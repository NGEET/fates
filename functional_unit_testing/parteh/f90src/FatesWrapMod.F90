! =======================================================================================
!
! This is the wrapper module that provides FATES data structures
!
! =======================================================================================

module EDTypesMod

  use iso_c_binding, only: fates_r8  => c_double
  use iso_c_binding, only: fates_int => c_int
  
  implicit none
  private

  integer(fates_int), parameter, public :: nlevleaf = 40
  real(fates_r8),     parameter, public :: dinc_ed = 1.0_fates_r8
  integer(fates_int), parameter, public :: nclmax   = 4

end module EDTypesMod


module shr_log_mod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int

  implicit none
  private

  ! Make necessary producers public
  public :: shr_log_errMsg

  contains

    function shr_log_errMsg(source, line) result(ans)
      character(kind=c_char,len=*), intent(in) :: source
      integer(c_int), intent(in) :: line
      character(kind=c_char,len=128) :: ans

      ans = "source: " // trim(source) // " line: "
    end function shr_log_errMsg
   
end module shr_log_mod

module FatesInterfaceMod

  use iso_c_binding, only: fates_r8  => c_double

  implicit none
  private

  real(fates_r8), parameter, public :: hlm_freq_day = 1.0_fates_r8/365.0_fates_r8
  
end module FatesInterfaceMod


module FatesGlobals

  implicit none
  private

  ! Make necessary producers public
  public :: fates_log
  public :: fates_endrun

  contains
    
    integer function fates_log()
      fates_log = 6      ! usually stdout
    end function fates_log
    
    subroutine fates_endrun(msg) 
      
      implicit none
      character(len=*), intent(in) :: msg    ! string to be printed
      stop
    end subroutine fates_endrun

end module FatesGlobals
