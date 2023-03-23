module FatesRunningMeanMod


  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : r8 => fates_r8
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
  
  type, public :: rmean_def_type

     real(r8) :: mem_period   ! The total integration period (s)
     real(r8) :: up_period    ! The period between updates (s)
     integer  :: n_mem        ! How many updates per integration period?
     integer  :: method       ! Is this a fixed or moving window?

   contains
     
     procedure :: define
     
  end type rmean_def_type

  
  ! This holds the time varying information for the mean
  ! which is instantiated on sites, patches, and cohorts
  
  type, public :: rmean_type
     
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

     integer  :: c_index ! The number of values that have
                         ! been added to the mean so far
                         ! once this is >= n_mem then
                         ! the ema weight hits its cap

     ! This points to the global structure that
     ! defines the nature of this mean/avg
     type(rmean_def_type), pointer :: def_type
     
   contains

     procedure :: GetMean
     procedure :: InitRMean
     procedure :: UpdateRMean
     procedure :: FuseRMean
     procedure :: CopyFromDonor
     
  end type rmean_type


  logical, parameter :: debug = .true.
  
  character(len=*), parameter, private :: sourcefile = &
         __FILE__


  ! Define the time methods that we want to have available to us
  
  class(rmean_def_type), public, pointer :: ema_24hr   ! Exponential moving average - 24hr window
  class(rmean_def_type), public, pointer :: fixed_24hr ! Fixed, 24-hour window
  class(rmean_def_type), public, pointer :: ema_lpa    ! Exponential moving average - leaf photo acclimation
  class(rmean_def_type), public, pointer :: ema_longterm  ! Exponential moving average - long-term leaf photo acclimation
  class(rmean_def_type), public, pointer :: ema_60day  ! Exponential moving average, 60 day
                                                       ! Updated daily
  class(rmean_def_type), public, pointer :: ema_storemem ! EMA used for smoothing N/C and P/C storage

  ! If we want to have different running mean specs based on
  ! pft or other types of constants
  type, public :: rmean_arr_type
     class(rmean_def_type), pointer :: p
  end type rmean_arr_type
  
contains


  subroutine define(this,mem_period,up_period,method)

    class(rmean_def_type) :: this

    real(r8),intent(in) :: mem_period
    real(r8),intent(in) :: up_period
    integer,intent(in)  :: method

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

    class(rmean_type) :: this
    real(r8)          :: GetMean

    if(this%def_type%method .eq. moving_ema_window) then
       if(this%c_index == 0 .and. debug) then
          write(fates_log(), *) 'attempting to get a running mean from a variable'
          write(fates_log(), *) 'that has not been given a value yet'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    GetMean = this%l_mean

  end function GetMean
  
  ! =====================================================================================
  
  subroutine InitRMean(this,rmean_def,init_value,init_offset)

    class(rmean_type)             :: this
    type(rmean_def_type), target  :: rmean_def
    real(r8),intent(in),optional  :: init_value
    real(r8),intent(in),optional  :: init_offset

    ! If the initialization happens part-way through a fixed averaging window
    ! we need to account for this.  The current method moves the position
    ! index to match the offset, and then assumes that the init_value provided
    ! was a constant over the offset period.

    ! If the first value is offset, such that the we are a portion of the
    ! way through the window, we need to account for this. 
    
    ! Point to the averaging type
    this%def_type => rmean_def

    if(this%def_type%method .eq. fixed_window) then

       if(debug) then
          if(.not.(present(init_offset).and.present(init_value)) )then
             write(fates_log(), *) 'when initializing a temporal mean on a fixed window'
             write(fates_log(), *) 'there must be an initial value and a time offset'
             write(fates_log(), *) 'specified.'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          ! Check to see if the offset is an even increment of the update frequency
          if( abs(real(nint(init_offset/rmean_def%up_period),r8)-(init_offset/rmean_def%up_period)) > nearzero ) then
             write(fates_log(), *) 'when initializing a temporal mean on a fixed window'
             write(fates_log(), *) 'the time offset must be an inrement of the update frequency'
             write(fates_log(), *) 'offset: ',init_offset
             write(fates_log(), *) 'up freq: ',rmean_def%up_period 
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          if(init_offset<-nearzero) then
             write(fates_log(), *) 'offset must be positive: ',init_offset
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
       end if

       this%c_index = modulo(nint(init_offset/rmean_def%up_period),rmean_def%n_mem)
       this%c_mean = real(this%c_index,r8)/real(rmean_def%n_mem,r8)*init_value
       this%l_mean = init_value

    elseif(this%def_type%method .eq. moving_ema_window) then
       
       if(present(init_value))then
          this%c_mean  = init_value
          this%l_mean  = init_value
          if(present(init_offset))then
             this%c_index = min(nint(init_offset/rmean_def%up_period),rmean_def%n_mem)
          else
             this%c_index = 1
          end if
       else
          this%c_mean  = nan
          this%l_mean  = nan
          this%c_index = 0
       end if
       
    end if
    
    return
  end subroutine InitRMean


  ! =====================================================================================

  
  subroutine CopyFromDonor(this, donor)

    class(rmean_type) :: this
    class(rmean_type),intent(in) :: donor

    if( .not.associated(this%def_type)) then
       write(fates_log(), *) 'Attempted to copy over running mean'
       write(fates_log(), *) 'info from a donor into a new structure'
       write(fates_log(), *) 'but the new structure did not have its'
       write(fates_log(), *) 'def_type pointer associated'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    this%c_mean  = donor%c_mean
    this%l_mean  = donor%l_mean
    this%c_index = donor%c_index

    
    return
  end subroutine CopyFromDonor


  
  ! =====================================================================================
  
  subroutine UpdateRMean(this, new_value)

    class(rmean_type) :: this
    real(r8)          :: new_value   ! The newest value added to the running mean
    real(r8)          :: wgt
    
    if(this%def_type%method.eq.moving_ema_window) then

       this%c_index = min(this%def_type%n_mem,this%c_index + 1)
       
       if(this%c_index==1) then
          this%c_mean = new_value
       else
          wgt = 1._r8/real(this%c_index,r8)
          this%c_mean = this%c_mean*(1._r8-wgt) + wgt*new_value
       end if
       
       this%l_mean = this%c_mean
       
    else

       ! If the last time we updated we had hit the
       ! end of the averaging memory period, and
       ! we are not using an indefinite running
       ! average, then zero things out
       
       this%c_index = this%c_index + 1
       wgt = this%def_type%up_period/this%def_type%mem_period
       this%c_mean = this%c_mean + new_value*wgt

       if(this%c_index == this%def_type%n_mem) then
          this%l_mean  = this%c_mean
          this%c_mean  = 0._r8
          this%c_index = 0
       end if

       
    end if

    return
  end subroutine UpdateRmean
  
  ! =====================================================================================
  
  subroutine FuseRMean(this,donor,recip_wgt)

    ! Rules for fusion:
    ! If both entities have valid means already, then you simply use the
    ! weight provided to combine them.
    ! If this is a moving average, then update the index to be the larger of
    ! the two.
    ! if this is a fixed window, check that the index is the same between
    ! both.

    
    class(rmean_type)          :: this
    class(rmean_type), pointer :: donor
    real(r8),intent(in)        :: recip_wgt  ! Weighting factor for recipient (0-1)

    ! Unecessary
    !if(this%def_type%n_mem .ne. donor%def_type%n_mem) then
    !   write(fates_log(), *) 'memory size is somehow different during fusion'
    !   write(fates_log(), *) 'of two running mean variables: '!,this%name,donor%name
    !   call endrun(msg=errMsg(sourcefile, __LINE__))
    !end if

    if(this%def_type%method .eq. fixed_window ) then
       if (this%c_index .ne. donor%c_index) then
          write(fates_log(), *) 'trying to fuse two fixed-window averages'
          write(fates_log(), *) 'that are at different points in the window?'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    ! This last logic clause just simply prevents us from doing math
    ! on uninitialized values. If both are unintiailized, then
    ! leave the result as uninitialized
    if( .not. (donor%c_index==0) ) then

       if(this%c_index==0) then
          this%c_mean = donor%c_mean
          this%l_mean = donor%l_mean 
          this%c_index = donor%c_index
       else
          ! Take the weighted mean between the two
          this%c_mean = this%c_mean*recip_wgt + donor%c_mean*(1._r8-recip_wgt)
          this%l_mean = this%l_mean*recip_wgt + donor%l_mean*(1._r8-recip_wgt)
          ! Update the index to the larger of the two
          this%c_index = max(this%c_index,donor%c_index)
       end if

    end if
       
    return
  end subroutine FuseRMean

  
end module FatesRunningMeanMod
