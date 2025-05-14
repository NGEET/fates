module FatesRunningSummMod


  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_huge
  use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  
  implicit none
  private

  ! These are flags that specify how the averaging window works.
  ! Exponential moving average (EMA) windows have an arbitrary size and update frequency)
  ! and it is technically never reset, it just averages indefinitely.
  ! But hourly, six-hourly, daily, monthly and yearly fixed windows have pre-set
  ! window sizes associated with their namesake, and more importantly, they
  ! are zero'd at the beginning of the interval, and get equal average weighting
  ! over their construction period.

  
  integer, public, parameter :: moving_ema_window = 0   ! (exponential moving average)
  integer, public, parameter :: fixed_window = 1
  
  ! This type defines a type of mean. It does not
  ! define the variable, but it defines how
  ! often it is updated, how long its
  ! memory period is, and if it should be zero'd
  ! These are globally defined on the proc.
  
  type, public :: rsumm_def_type

     real(r8) :: mem_period   ! The total integration period (s)
     real(r8) :: up_period    ! The period between updates (s)
     integer  :: n_mem        ! How many updates per integration period?
     integer  :: method       ! Is this a fixed or moving window?

   contains
     
     procedure :: define
     
  end type rsumm_def_type

  
  ! This holds the time varying information for the mean
  ! which is instantiated on sites, patches, and cohorts
  
  type, public :: rsumm_type
     
     real(r8) :: c_mean  ! The current mean value, if this
                         ! is a moving window, it is the mean.
                         ! If this is a fixed window, it is only a partial mean
                         ! as the value uses equal update weights and is not
                         ! necessarily fully constructed.

     real(r8) :: l_mean  ! The latest reportable mean value
                         ! this value is actually the same
                         ! as c_mean for moving windows, and for fixed windows
                         ! it is the mean value when the time integration window
                         ! last completed.
     
     real(r8) :: c_minimum  ! The current minimum value. Depending on the time,
                            ! this will be a partial minimum as it needs to 
                            ! compare with all to fully determine the minimum

     real(r8) :: l_minimum  ! The latest reportable minimum value. This is the
                            ! minimum value for the last window period that was
                            ! fully assessed.
     
     real(r8) :: c_maximum  ! The current maximum value. Depending on the time,
                            ! this will be a partial minimum as it needs to 
                            ! compare with all to fully determine the maximum

     real(r8) :: l_maximum  ! The latest reportable maximum value. This is the
                            ! maximum value for the last window period that was
                            ! fully assessed.

     integer  :: c_index ! The number of values that have
                         ! been added to the mean so far
                         ! once this is >= n_mem then
                         ! the ema weight hits its cap

     ! This points to the global structure that
     ! defines the nature of this mean/avg
     type(rsumm_def_type), pointer :: def_type
     
   contains

     procedure :: GetMean
     procedure :: GetMin
     procedure :: GetMax
     procedure :: InitRSumm
     procedure :: UpdateRSumm
     procedure :: FuseRSumm
     procedure :: CopyFromDonor
     
  end type rsumm_type

  logical, parameter :: debug = .true.
  
  character(len=*), parameter, private :: sourcefile = &
         __FILE__


  ! Define the time methods that we want to have available to us
  
  class(rsumm_def_type), public, pointer :: ema_24hr   ! Exponential moving average - 24hr window
  class(rsumm_def_type), public, pointer :: fixed_24hr ! Fixed, 24-hour window
  class(rsumm_def_type), public, pointer :: ema_lpa    ! Exponential moving average - leaf photo acclimation
  class(rsumm_def_type), public, pointer :: ema_longterm  ! Exponential moving average - long-term leaf photo acclimation
  class(rsumm_def_type), public, pointer :: ema_60day  ! Exponential moving average, 60 day
                                                       ! Updated daily
  class(rsumm_def_type), public, pointer :: ema_storemem ! EMA used for smoothing N/C and P/C storage
  class(rsumm_def_type), public, pointer :: ema_sdlng_mort_par  ! EMA for seedling mort from light stress
  class(rsumm_def_type), public, pointer :: ema_sdlng2sap_par  ! EMA for seedling to sapling transition rates
                                                               ! based in par
  class(rsumm_def_type), public, pointer :: ema_sdlng_emerg_h2o  ! EMA for moisture-based seedling emergence
  class(rsumm_def_type), public, pointer :: ema_sdlng_mdd  ! EMA for seedling moisture deficit days 
  
  
  ! If we want to have different running mean specs based on
  ! pft or other types of constants
  type, public :: rsumm_arr_type
     class(rsumm_type), pointer :: p
  end type rsumm_arr_type
  

contains


  subroutine define(this,mem_period,up_period,method)
    !  This sub-routine sets parameters and method for this running summary structure

    class(rsumm_def_type) :: this

    real(r8),intent(in) :: mem_period ! The total integration period (s)
    real(r8),intent(in) :: up_period  ! The period between updates (s)
    integer,intent(in)  :: method     ! Is this a fixed or moving window?

    ! Check the memory and update periods
    if(debug) then
       if( abs(nint(mem_period/up_period)-mem_period/up_period) > nearzero ) then
          write(fates_log(), *) 'While defining a running mean definition'
          write(fates_log(), *) 'an update and memory period was specified'
          write(fates_log(), *) 'where the update period is not an exact fraction of the period'
          write(fates_log(), *) 'mem_period: ',mem_period
          write(fates_log(), *) 'up_period: ',up_period
          write(fates_log(), *) 'exiting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
       
    this%mem_period = mem_period
    this%up_period  = up_period
    this%method     = method
    this%n_mem      = nint(mem_period/up_period)

    return
  end subroutine define

  ! =====================================================================================
  
  function GetMean(this)
    ! Function to retrieve the running average for this variable.

    class(rsumm_type) :: this
    real(r8)          :: GetMean

    select case(this%def_type%method)
    case (moving_ema_window)
       if(this%c_index == 0 .and. debug) then
          write(fates_log(), *) 'attempting to get a running mean from a variable'
          write(fates_log(), *) 'that has not been given a value yet'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end select
    GetMean = this%l_mean

  end function GetMean

  ! =====================================================================================
  
  function GetMin(this)
    ! Function to retrieve the minimum value in this time window for this variable.

    class(rsumm_type) :: this
    real(r8)          :: GetMin

    select case(this%def_type%method)
    case (moving_ema_window)
       write(fates_log(), *) '   Attempting to get the minimum value from a variable'
       write(fates_log(), *) 'that is integrated using exponential moving averages.'
       call endrun(msg=errMsg(sourcefile, __LINE__))

    case (fixed_value)
       GetMin = this%l_minimum
    end select


  end function GetMin

  ! =====================================================================================
  
  function GetMax(this)
    ! Function to retrieve the maximum value in this time window for this variable.

    class(rsumm_type) :: this
    real(r8)          :: GetMax

    select case(this%def_type%method)
    case (moving_ema_window)
       write(fates_log(), *) '   Attempting to get the maximum value from a variable'
       write(fates_log(), *) 'that is integrated using exponential moving averages.'
       call endrun(msg=errMsg(sourcefile, __LINE__))

    case (fixed_value)
       GetMax = this%l_maximum
    end select


  end function GetMax

  ! =====================================================================================
  
  subroutine InitRSumm(this,rsumm_def,init_value,init_offset)
    !    Sub-routine to initialize the summary structure variables.
    ! If the initialization happens part-way through a fixed averaging window
    ! we need to account for this.  The current method moves the position
    ! index to match the offset, and then assumes that the init_value provided
    ! was a constant over the offset period.

    class(rsumm_type)             :: this
    type(rsumm_def_type), target  :: rsumm_def
    real(r8),intent(in),optional  :: init_value
    real(r8),intent(in),optional  :: init_offset
    
    ! Point to the averaging type
    this%def_type => rsumm_def

    select case(this%def_type%method)
    case (fixed_window)

       if(debug) then
          if(.not.(present(init_offset).and.present(init_value)) )then
             write(fates_log(), *) 'when initializing a temporal mean on a fixed window'
             write(fates_log(), *) 'there must be an initial value and a time offset'
             write(fates_log(), *) 'specified.'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          ! Check to see if the offset is an even increment of the update frequency
          if( abs(real(nint(init_offset/rsumm_def%up_period),r8)-(init_offset/rsumm_def%up_period)) > nearzero ) then
             write(fates_log(), *) 'when initializing a temporal mean on a fixed window'
             write(fates_log(), *) 'the time offset must be an increment of the update frequency'
             write(fates_log(), *) 'offset: ',init_offset
             write(fates_log(), *) 'up freq: ',rsumm_def%up_period 
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          if(init_offset<-nearzero) then
             write(fates_log(), *) 'offset must be positive: ',init_offset
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
       end if


       ! If the first value is offset, such that the we are a portion of the
       ! way through the window, we need to account for this. 
       this%c_index = modulo(nint(init_offset/rsumm_def%up_period),rsumm_def%n_mem)
       this%c_mean = real(this%c_index,r8)/real(rsumm_def%n_mem,r8)*init_value
       this%l_mean = init_value

       !    Minima and maxima are not weighted so they can both assigned the
       ! initial value directly.
       this%c_minimum = init_value
       this%l_minimum = init_value
       this%c_maximum = init_value
       this%l_maximum = init_value

    case (moving_ema_window)
       if (present(init_value)) then
          this%c_mean  = init_value
          this%l_mean  = init_value
          if(present(init_offset))then
             this%c_index = min(nint(init_offset/rsumm_def%up_period),rsumm_def%n_mem)
          else
             this%c_index = 1
          end if
       else
          this%c_mean  = nan
          this%l_mean  = nan
          this%c_index = 0
       end if

       ! We do not track minima and maxima when handling exponential moving averages.
       this%c_minimum = nan
       this%l_minimum = nan
       this%c_maximum = nan
       this%l_maximum = nan
    end select

    return
  end subroutine InitRSumm


  ! =====================================================================================

  
  subroutine CopyFromDonor(this, donor)
    ! Sub-routine to transfer running summary information from one structure to another

    class(rsumm_type) :: this             ! Receptor structure
    class(rsumm_type),intent(in) :: donor ! Donor structure

    if( .not.associated(this%def_type)) then
       write(fates_log(), *) 'Attempted to copy over running mean'
       write(fates_log(), *) 'info from a donor into a new structure'
       write(fates_log(), *) 'but the new structure did not have its'
       write(fates_log(), *) 'def_type pointer associated'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    this%c_mean    = donor%c_mean
    this%l_mean    = donor%l_mean
    this%c_minimum = donor%c_minimum
    this%l_minimum = donor%l_minimum
    this%c_maximum = donor%c_maximum
    this%l_maximum = donor%l_maximum
    this%c_index   = donor%c_index


    return
  end subroutine CopyFromDonor


  
  ! =====================================================================================
  
  subroutine UpdateRSumm(this, new_value)
    ! Sub-routine that updates the running summary metrics.

    class(rsumm_type) :: this
    real(r8)          :: new_value   ! The newest value added to the running mean
    real(r8)          :: wgt
    
    select case (this%def_type%method)
    case (moving_ema_window)

       this%c_index = min(this%def_type%n_mem,this%c_index + 1)

       if(this%c_index==1) then
          this%c_mean = new_value
       else
          wgt = 1._r8/real(this%c_index,r8)
          this%c_mean = this%c_mean*(1._r8-wgt) + wgt*new_value
       end if

       this%l_mean = this%c_mean

    case (fixed_window)

       ! If the last time we updated we had hit the
       ! end of the averaging memory period, and
       ! we are not using an indefinite running
       ! average, then zero things out
       
       this%c_index = this%c_index + 1
       wgt = this%def_type%up_period/this%def_type%mem_period
       this%c_mean = this%c_mean + new_value*wgt
       this%c_minimum = min(this%c_minimum,new_value)
       this%c_maximum = max(this%c_maximum,new_value)

       if(this%c_index == this%def_type%n_mem) then
          this%l_mean    = this%c_mean
          this%c_mean    = 0._r8
          this%c_minimum = fates_huge
          this%c_maximum = -fates_huge
          this%c_index = 0
       end if

    end select

    return
  end subroutine UpdateRSumm
  
  ! =====================================================================================
  
  subroutine FuseRSumm(this,donor,recip_wgt)
    ! Sub-routine that fuses data from two running summary structures.
    !
    ! Rules for fusion:
    ! If both entities have valid means already, then you simply use the
    ! weight provided to combine them.
    ! If this is a moving average, then update the index to be the larger of
    ! the two.
    ! if this is a fixed window, check that the index is the same between
    ! both.

    
    class(rsumm_type)          :: this
    class(rsumm_type), pointer :: donor
    real(r8),intent(in)        :: recip_wgt  ! Weighting factor for recipient (0-1)



    select case (this%def_type%method)
    case (fixed_window)
       if (this%c_index .ne. donor%c_index) then
          write(fates_log(), *) 'trying to fuse two fixed-window averages'
          write(fates_log(), *) 'that are at different points in the window?'
          write(fates_log(), *) 'c_mean', this%c_mean,  donor%c_mean
          write(fates_log(), *) 'l_mean', this%l_mean,  donor%l_mean
          write(fates_log(), *) 'c_index', this%c_index,  donor%c_index
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end select

    ! This last logic clause just simply prevents us from doing math
    ! on uninitialized values. If both are unintiailized, then
    ! leave the result as uninitialized
    if ( donor%c_index /= 0 ) then

       if ( this%c_index == 0 ) then
          this%c_mean    = donor%c_mean
          this%l_mean    = donor%l_mean
          this%c_minimum = donor%c_minimum
          this%l_minimum = donor%l_minimum
          this%c_maximum = donor%c_maximum
          this%l_maximum = donor%l_maximum
          this%c_index   = donor%c_index
       else
          ! Take the weighted mean between the two
          this%c_mean = this%c_mean*recip_wgt + donor%c_mean*(1._r8-recip_wgt)
          this%l_mean = this%l_mean*recip_wgt + donor%l_mean*(1._r8-recip_wgt)

          !    Minima and maxima are fused only for fixed window, otherwise they are
          ! not computed and must remain nan
          select case (this%def_type%method)
          case (fixed_window)
             this%c_minimum = this%c_minimum*recip_wgt + donor%c_minimum*(1._r8-recip_wgt)
             this%l_minimum = this%l_minimum*recip_wgt + donor%l_minimum*(1._r8-recip_wgt)
             this%c_maximum = this%c_maximum*recip_wgt + donor%c_maximum*(1._r8-recip_wgt)
             this%l_maximum = this%l_maximum*recip_wgt + donor%l_maximum*(1._r8-recip_wgt)
          end select

          ! Update the index to the larger of the two
          this%c_index = max(this%c_index,donor%c_index)
       end if

    end if
       
    return
  end subroutine FuseRSumm

  
end module FatesRunningSummMod
