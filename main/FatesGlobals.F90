module FatesGlobals
  ! NOTE(bja, 201608) This is a temporary hack module to store global
  ! data used inside fates. It's main use it to explicitly call out
  ! global data that needs to be dealt with, but doesn't have an
  ! immediately obvious home.

  use FatesConstantsMod         , only : r8 => fates_r8
  use TwoStreamMLPEMod          , only : TwoStreamLogInit
   
  implicit none
  private         ! By default everything is private

  integer :: fates_log_
  logical :: fates_global_verbose_

  ! Make public necessary subroutines and functions
  public :: FatesGlobalsInit
  public :: fates_log
  public :: fates_global_verbose
  public :: fates_endrun
  public :: FatesWarn
  public :: A2S
  public :: N2S
  public :: I2S
  public :: FatesReportTotalWarnings
  
  ! -------------------------------------------------------------------------------------
  ! Warning handling
  ! The objective here is to stop writing the same warning over and over again. After
  ! we've seen the same machine print out the same warning over and over again, we get
  ! the point and don't have to continue seeing the message.
  ! We also allow warnings to have their own unique or group identifier, which will
  ! make is to you only turn off warnings that are particularly chatty, and continue
  ! to allow warnings elsewhere that have not tripped as often.
  ! -------------------------------------------------------------------------------------
  
  integer, parameter :: max_ids = 200              ! Maximum number of unique warning ids
                                                   ! expand as necessary
  integer :: warn_counts(0:max_ids) = 0            ! Total number of times each id has warned
  integer, parameter :: max_warnings = 100         ! The maximum number of warnings before we
                                                   ! stop writing the warning
  logical :: warn_active(0:max_ids) = .true.       ! The current status of the warning. 
  logical, parameter :: warning_override = .false. ! If you really don't want any warnings
                                                   ! you can set this to true to avoid
                                                   ! printing any of these warnings to the log
                                                   ! It should also bypass the logicals bound inside
                                                   ! at the compiler level (?) and be faster

 
  
contains



  ! =====================================================================================

  subroutine FatesGlobalsInit(log_unit,global_verbose)

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    fates_log_ = log_unit
    fates_global_verbose_ = global_verbose

    call TwoStreamLogInit(log_unit)

  end subroutine FatesGlobalsInit

  ! =====================================================================================

  integer function fates_log()
    fates_log = fates_log_
  end function fates_log

  logical function fates_global_verbose()
    fates_global_verbose = fates_global_verbose_
  end function fates_global_verbose

  subroutine fates_endrun(msg, additional_msg) 

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    ! This subroutine was derived from CLM's
    ! endrun_vanilla() in abortutils.F90
    !
    use shr_sys_mod , only: shr_sys_abort
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: msg            ! string to be printed
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_sys_abort
    !-----------------------------------------------------------------------

    if (present(additional_msg)) then 
      write(fates_log(),*) 'ENDRUN: ', trim(additional_msg)
    else 
      write(fates_log(),*) 'ENDRUN: '
    end if
    
    call shr_sys_abort(msg)

  end subroutine fates_endrun

  ! =====================================================================================

  subroutine FatesWarn(msg,index)

    character(len=*), intent(in) :: msg      ! string to be printed
    integer,optional,intent(in)  :: index    ! warning index

    integer :: ind

    if(warning_override) return  ! Exit early if we are turning off warnings
    
    if(present(index))then
       ind = index
    else
       ind = 0
    end if

    ! Don't check if the index is within bounds, this routine could already
    ! be too expensive if this is in cohort loops
    warn_counts(ind) = warn_counts(ind) + 1

    if(warn_active(ind))then
       write(fates_log(),*) 'FATESWARN: '//trim(ADJUSTL(I2S(ind)))//' m: '//trim(msg)
       if(warn_counts(ind)> max_warnings) then
          warn_active(ind) = .false.
          write(fates_log(),*) 'FATESWARN: '//trim(ADJUSTL(I2S(ind)))//' has saturated messaging, no longer reporting'
       end if
    end if
    return
  end subroutine FatesWarn

  ! =====================================================================================


  subroutine FatesReportTotalWarnings()
    
    integer :: ind

    do ind = 0,max_ids
       
       if(warn_counts(ind)>0)then

          write(fates_log(),*) 'FATESWARN: '//trim(ADJUSTL(I2S(ind)))//' was triggered ',trim(ADJUSTL(I2S(warn_counts(ind))))//' times'
          
       end if
       
    end do

    
  end subroutine FatesReportTotalWarnings
  
  
  ! =====================================================================================

  function N2S(real_in) result(str)

    real(r8) :: real_in
    character(len=16) :: str

    !write(str,*) real_in
    write(str,'(E12.6)') real_in
    
  end function N2S

  ! =====================================================================================
  
  function I2S(int_in) result(str)

    integer :: int_in
    character(len=16) :: str

    !write(str,*) real_in
    write(str,'(I15)') int_in
    
  end function I2S
  
  ! =====================================================================================
  
  function A2S(reals_in) result(str)

    real(r8) :: reals_in(:)
    character(len=1024) :: str
    integer :: i, nreal
    
    str = ', '
    do i = 1,ubound(reals_in,1)
       str = trim(str)//', '//N2S(reals_in(i))
    end do
    
  end function A2S
  
end module FatesGlobals
