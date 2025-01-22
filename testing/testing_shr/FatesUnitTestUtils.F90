module FatesUnitTestUtils

  ! Miscellaneous utilities to aid unit testing

  implicit none
  private

  public :: endrun_msg ! Gives the message thrown by shr_abort_abort, given a call to endrun(msg)

contains

  function endrun_msg(msg)
    !
    ! DESCRIPTION:
    ! Gives the message thrown by shr_abort_abort, given a call to endrun(msg)
    !
    ! ARGUMENTS:
    character(len=:), allocatable :: endrun_msg  ! function result
    character(len=*), intent(in)  :: msg

    ! LOCALS:
    character(len=*), parameter :: subname = 'endrun_msg'

    endrun_msg = 'ENDRUN: '//trim(msg)

  end function endrun_msg

end module FatesUnitTestUtils