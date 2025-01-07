module TimingMod
  !
  !  DESCRIPTION
  !  this module implements a timer class for testing purposes
  !
  
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  private 
  
  type, public :: testing_timer 

    private
    real(r8) :: saved_time  ! saved time in [ms]

    contains

      procedure, public :: start_timer => start_timer_sub
      procedure, public :: elapsed_time => elapsed_time_fn

  end type testing_timer

  private :: start_timer_sub, elapsed_time_fn
  
  contains 
  
  !---------------------------------------------------------------------------------------
  
  subroutine start_timer_sub(this)
    !
    ! DESCRIPTION:
    ! Get and save the initial time
    !
    
    ! ARGUMENTS:
    class(testing_timer) :: this ! timer object

    ! LOCALS:
    integer, dimension(8) :: vals ! time value array

    ! get time
    call date_and_time(VALUES=vals)
    
    this%saved_time = 86400.D0*vals(3) + 3600.D0*vals(5) + 60.D0*vals(6) + vals(7) +     &
      0.001D0*vals(8)

  end subroutine start_timer_sub
  
  !---------------------------------------------------------------------------------------
  
  real(r8) function elapsed_time_fn(this)
    !
    ! DESCRIPTION:
    ! calculate elapsed time
    !
    
    ! ARGUMENTS:
    class(testing_timer) :: this ! timer object

    ! LOCALS:
    integer, dimension(8) :: vals         ! time value array
    real(r8)              :: current_time ! current time [ms]

    ! get time
    call date_and_time(VALUES=vals)
    current_time = 86400.D0*vals(3) + 3600.D0*vals(5) + 60.D0*vals(6) + vals(7) +        &
      0.001D0*vals(8)

    ! get elapsed time in seconds
    elapsed_time_fn = current_time - this%saved_time

  end function elapsed_time_fn
  
end module TimingMod