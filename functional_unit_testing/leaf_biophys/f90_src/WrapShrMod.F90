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

module shr_sys_mod

  public :: shr_sys_abort

contains

  subroutine shr_sys_abort
    call exit(0)
  end subroutine shr_sys_abort
  
end module shr_sys_mod

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
  
  subroutine FatesWarn(msg,index)
    
    character(len=*), intent(in) :: msg      ! string to be printed
    integer,optional,intent(in)  :: index    ! warning index
    
    
    write(stdo_unit,*) index, msg
    
    stop
    
  end subroutine FatesWarn
  
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
    integer :: i
    
    str = ', '
    do i = 1,ubound(reals_in,1)
       str = trim(str)//', '//N2S(reals_in(i))
    end do
    
  end function A2S
end module FatesGlobals
