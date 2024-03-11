
! =======================================================================================
!
! This file is an alternative to key files in the fates
! filesystem. Noteably, we replace fates_r8 and fates_in
! with types that work with "ctypes".  This is
! a key step in working with python
!
! We also wrap FatesGlobals to reduce the dependancy
! cascade that it pulls in from shr_log_mod.
!
! =======================================================================================

module shr_log_mod

   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int

   contains

   function shr_log_errMsg(source, line) result(ans)
      character(kind=c_char,len=*), intent(in) :: source
      integer(c_int), intent(in) :: line
      character(kind=c_char,len=128) :: ans

      ans = "source: " // trim(source) // " line: "
   end function shr_log_errMsg

end module shr_log_mod


module FatesGlobals

   contains

   integer function fates_log()
      fates_log = -1
   end function fates_log

   subroutine fates_endrun(msg)

      implicit none
      character(len=*), intent(in) :: msg    ! string to be printed

      stop

   end subroutine fates_endrun

end module FatesGlobals
