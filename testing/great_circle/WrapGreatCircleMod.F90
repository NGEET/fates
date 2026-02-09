module shr_log_mod
  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  
  public :: shr_log_errMsg
  
contains
  function shr_log_errMsg(source, line) result(ans)
    character(kind=c_char,len=*), intent(in) :: source
    integer(c_int), intent(in) :: line
    character(kind=c_char,len=4) :: cline  ! character version of int
    character(kind=c_char,len=128) :: ans
    
    write(cline,'(I4)') line
    ans = "source: " // trim(source) // " line: "// trim(cline)
    
  end function shr_log_errMsg
  
end module shr_log_mod

module FatesGlobals
  
  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use FatesConstantsMod, only : r8 => fates_r8

  integer :: stdo_unit = 6
  
contains
  
  integer function fates_log()
    fates_log = 6
  end function fates_log
  
  subroutine fates_endrun(msg)
    
    implicit none
    character(len=*), intent(in) :: msg    ! string to be printed
    
    write(stdo_unit,*) msg
    
    stop
    
  end subroutine fates_endrun
end module FatesGlobals
