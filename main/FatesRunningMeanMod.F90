module FatesRunningMeanMod


  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : r8 => fates_r8
  use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  
  implicit none
  private

  integer, parameter :: maxlen_varname = 8
  
  type, public :: rmean_type
     
     real(r8), allocatable    :: mem(:)        ! Array storing the memory of the mean
     real(r8)                 :: c_mean        ! The current mean value from the
                                               ! available memory, as of the last update
     integer                  :: c_index       ! The current memory index as per the
                                               ! last update
     integer                  :: n_mem         ! Total number of memory indices
     logical                  :: filled        ! Has enough time elapsed so all memory filled?
     !character(len=maxlen_varname) :: var_name ! A short name for this variable
     !                                          ! for diagnostic purposes
     
   contains

     procedure :: get_mean
     procedure :: InitRMean
     procedure :: UpdateRMean
     procedure :: FuseRMean
     
  end type rmean_type

  character(len=*), parameter, private :: sourcefile = &
         __FILE__
  
contains


  function get_mean(this)

    class(rmean_type)             :: this
    real(r8) :: get_mean

    if(this%c_index == 0) then
       write(fates_log(), *) 'attempting to get a running mean from a variable'
       write(fates_log(), *) 'that has not experienced any time to accumluate memory'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    get_mean = this%c_mean

  end function get_mean
  

  subroutine InitRMean(this,mem_period,up_period) !,init_value)

    class(rmean_type)             :: this
    !character(len=maxlen_varname) :: name       ! The name of the new variable
    real(r8)                      :: mem_period ! The period length in seconds that must be remembered
    real(r8)                      :: up_period  ! The period length in seconds that memory is updated
    ! (i.e. the resolution of the memory)
    !real(r8)                      :: init_value ! The first value to put in memory
    
    !this%name = name
    this%n_mem = nint(mem_period/up_period)
    
    if( abs(real(this%n_mem,r8)-mem_period/up_period) > nearzero ) then
       write(fates_log(), *) 'While defining a running mean variable: '!,this%var_name
       write(fates_log(), *) 'an update and total memory period was specified'
       write(fates_log(), *) 'where the update period is not an exact fraction of the period'
       write(fates_log(), *) 'mem_period: ',mem_period
       write(fates_log(), *) 'up_period: ',up_period
       write(fates_log(), *) 'exiting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Otherwise we allocate
    allocate(this%mem(this%n_mem))

    ! Initialize with nonsense numbers
    this%mem(:) = nan

    ! There are no memory spots filled with valid datapoints
    this%filled = .false.

    ! The current index of the memory is one
    this%c_index = 0
    
    !this%mem(1) = init_value

    
    return
  end subroutine InitRMean

  ! =====================================================================================
  
  subroutine UpdateRMean(this, new_value)

    class(rmean_type) :: this
    real(r8)          :: new_value   ! The newest value added to the running mean
    
    this%c_index = this%c_index + 1
    
    ! Update the index of the memory array
    ! and, determine if we have filled the memory yet
    ! If this index is greater than our memory slots,
    ! go back to the first
    if(this%c_index==this%n_mem) then
       this%filled  = .true.
    end if

    if(this%c_index>this%n_mem) this%c_index = 1
    
    this%mem(this%c_index) = new_value

    ! Update the running mean value. It will return a value
    ! if we have not filled in all the memory slots. To do this
    ! it will take a mean over what is available

    if(this%filled) then
       this%c_mean = sum(this%mem)/real(this%n_mem,r8)
    else
       this%c_mean = sum(this%mem(1:this%c_index))/real(this%c_index,r8)
    end if

    return
  end subroutine UpdateRmean
  
  ! =====================================================================================
  
  subroutine FuseRMean(this,donor,recip_wgt)

    ! When fusing the running mean of two entities, it is possible that they
    ! may have a different amount of memory spaces filled (at least in FATES). This
    ! is typical for newly created patches or cohorts, that litteraly just spawned
    ! So what generally happens is we walk backwards from the current memory index
    ! of both and take means where we can.  In places where one entity has more
    ! memory than the other, than we just use the value from the one that is there

    class(rmean_type)          :: this
    class(rmean_type), pointer :: donor
    real(r8),intent(in)        :: recip_wgt  ! Weighting factor for recipient (0-1)
    
    integer :: r_id   ! Loop index counter for the recipient (this)
    integer :: d_id   ! Loop index counter for the donor

    
    if(this%n_mem .ne. donor%n_mem) then
       write(fates_log(), *) 'memory size is somehow different during fusion'
       write(fates_log(), *) 'of two running mean variables: '!,this%name,donor%name
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if(this%filled .and. donor%filled) then

       ! Both are filled, take averages of each and be sure to use relative positions
       
       d_id = donor%c_index
       do r_id = 1,this%n_mem
          this%mem(r_id) = this%mem(r_id)*(recip_wgt) + &
               donor%mem(d_id)*(1._r8-recip_wgt)
          d_id = d_id+1
          if(d_id>donor%n_mem) d_id = 1
       end do

    elseif(this%filled .and. .not.donor%filled) then

       ! Here, we have only partial memory from the donor
       ! we we iterate through the donor's memory
       ! and average between the two in those spaces. Then
       ! we leave the rest untouched because we accept
       ! the values from the recipient. Also, keep the
       ! recipient's index.
       
       r_id = this%c_index
       do d_id = donor%c_index,1,-1
          this%mem(r_id) = this%mem(r_id)*(recip_wgt) + donor%mem(d_id)*(1._r8-recip_wgt)
          r_id = r_id - 1
          if(r_id==0) r_id = this%n_mem
       end do


    elseif(.not.this%filled .and. donor%filled) then

       ! Here we only have partial memory from the recipient
       ! so we iterate through the recipient's memory
       ! and average between the two.  Then we copy
       ! over the values from the donor for the indices that we
       ! didn't average because the donor has valid values
       
       d_id = donor%c_index
       do r_id = this%c_index,1,-1
          this%mem(r_id) = this%mem(r_id)*(recip_wgt) + donor%mem(d_id)*(1._r8-recip_wgt)
          d_id = d_id - 1
          if(d_id==0) d_id = donor%n_mem
       end do

       do r_id = this%n_mem,this%c_index+1,-1
          this%mem(r_id) =  donor%mem(d_id)
          d_id = d_id - 1
       end do
       ! Pass the current index of the donor since that was filled
       ! And also update the status to filled since it now
       ! has all memory filled from the donor
       this%c_index = donor%c_index
       this%filled  = .true.
       
    elseif(.not.this%filled .and. .not.donor%filled) then

       ! Here, neither is completely filled
       
       if( this%c_index>donor%c_index ) then

          ! In this case, leave all data as the recipient
          ! except for where there is donor. Keep the recipient's
          ! index and status, since it is larger and should remain unchanged
          
          r_id = this%c_index
          do d_id = donor%c_index,1,-1
             this%mem(r_id) = this%mem(r_id)*(recip_wgt) + donor%mem(d_id)*(1._r8-recip_wgt)
          end do
          
       else

          ! In this case, we do the same thing as the previous
          ! clause, but just switch the roles, then
          ! copy the donor to the recipient
          ! Also transfer the index from the donor, since that was
          ! higher and now reflects the filled memory
          
          d_id = donor%c_index
          do r_id = this%c_index,1,-1
             donor%mem(r_id) = this%mem(r_id)*(recip_wgt) + donor%mem(d_id)*(1._r8-recip_wgt)
          end do
          this%mem(:) = donor%mem
          this%c_index = donor%c_index
       end if
       
    end if

    ! Update the mean based on the fusion
    if(this%filled) then
       this%c_mean = sum(this%mem)/real(this%n_mem,r8)
    else
       if(this%c_index>0) this%c_mean = &
            sum(this%mem(1:this%c_index))/real(this%c_index,r8)
    end if

    
    return
  end subroutine FuseRMean

  
end module FatesRunningMeanMod
