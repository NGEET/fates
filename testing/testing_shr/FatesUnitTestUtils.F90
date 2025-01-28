module FatesUnitTestUtils 
  
  ! DESCRIPTION:
  !  Miscellaneous methods to aid in unit testing
  !
  
  implicit none
  private
  
  public :: endrun_msg

contains

  !---------------------------------------------------------------------------------------
  
  function endrun_msg(msg)
    !
    ! DESCRIPTION:
    ! Gives the message thrown by shr_abort_abort, given a call to endrun(msg)
    !

    ! ARGUMENTS:
    character(len=:), allocatable :: endrun_msg  ! function result
    character(len=*), intent(in)  :: msg

    endrun_msg = 'ABORTED: '//trim(msg)

  end function endrun_msg
  
  !---------------------------------------------------------------------------------------
  
end module FatesUnitTestUtils