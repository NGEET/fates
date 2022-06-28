module FatesUtilsMod
  
  ! This module contains helper functions and subroutines which are general in nature.
  ! Think string parsing, timing, maybe numerics, etc.
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only      : fates_log

  implicit none
  private ! Modules are private by default

  ! Make public necessary subroutines and functions
  public :: check_hlm_list
  public :: check_var_real

contains
  
  
  function check_hlm_list(hlms,hlm_name) result(astatus)
    
    ! ---------------------------------------------------------------------------------
    ! This simple function compares a string of HLM tags to see if any of the names
    ! match the name of the currently active HLM. If any do, return true, if any
    ! don't, if any don't its a big secret.
    ! ---------------------------------------------------------------------------------
    
    character(len=*),intent(in) :: hlms
    character(len=*),intent(in) :: hlm_name

    integer :: index
    logical :: astatus

    astatus = .false.
    index = scan(trim(hlms),trim(hlm_name))
 
    if(index>0)then
       astatus=.true.
    end if
    return

  end function check_hlm_list
  
  ! =====================================================================================

  subroutine check_var_real(r8_var, var_name, return_code)

     real(r8),intent(in)         :: r8_var
     character(len=*),intent(in) :: var_name
     integer,intent(out)         :: return_code

     real(r8), parameter :: r8_type  = 1.0
     real(r8), parameter :: overflow  = huge(r8_type)
     real(r8), parameter :: underflow = tiny(r8_type)

     return_code = 0

     ! NaN check
     if (r8_var /= r8_var) then
        write(fates_log(),*) 'NaN detected, ',trim(var_name),': ',r8_var
        return_code = 1
     end if
     
     ! Overflow check (within 100th of max precision)
     if (abs(r8_var) > 0.01*overflow) then
        write(fates_log(),*) 'Nigh overflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 10
     end if

     ! Underflow check (within 100x of min precision)
     if (abs(r8_var) < 100.0_r8*underflow) then
        write(fates_log(),*) 'Nigh underflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 100
     end if


  end subroutine check_var_real
  
  !==========================================================================================!
  !   Function to compute the great circle distance between two points: the s suffix denotes !
  ! source point, and f denotes the destination - "forepoint"). The results are given in     !
  ! metres. The formula is intended to be accurate for both small and large distances and    !
  ! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
  ! close to the n*pi/2.                                                                     !
  !------------------------------------------------------------------------------------------!
  real function dist_gc(slons,slonf,slats,slatf)
  use consts_coms, only : erad    & ! intent(in)
                        , pio1808 ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real, intent(in) :: slons
   real, intent(in) :: slonf
   real, intent(in) :: slats
   real, intent(in) :: slatf
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)     :: lons
   real(kind=8)     :: lonf
   real(kind=8)     :: lats
   real(kind=8)     :: latf
   real(kind=8)     :: dlon
   real(kind=8)     :: dlat
   real(kind=8)     :: x
   real(kind=8)     :: y
   !---------------------------------------------------------------------------------------!

   !----- Convert the co-ordinates to double precision and to radians. --------------------!
   lons = dble(slons) * pio1808
   lonf = dble(slonf) * pio1808
   lats = dble(slats) * pio1808
   latf = dble(slatf) * pio1808
   dlon = lonf - lons
   dlat = latf - lats

   !----- Find the arcs. ------------------------------------------------------------------!
   x    = dsin(lats) * dsin(latf) + dcos(lats) * dcos(latf) * dcos(dlon)
   y    = dsqrt( (dcos(latf)*dsin(dlon)) * (dcos(latf)*dsin(dlon))                         &
               + (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))                  &
               * (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)) )

   !----- Convert the arcs to actual distance. --------------------------------------------!
   dist_gc = erad*sngl(datan2(y,x))

   return
   end function dist_gc


end module FatesUtilsMod
