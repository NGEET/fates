
! THIS IS A STRIPPED DOWN VERSION OF main/EDParamsMod.F90


module EDParamsMod
   !
   ! module that deals with reading the ED parameter file
   !
   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int
   use FatesConstantsMod, only : r8 => fates_r8

   implicit none
   private
   save

   integer(kind=c_int), parameter :: param_string_length = 32


   ! Hydraulics Control Parameters
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: hydr_kmax_rsurf1         !  maximum conducitivity for unit root surface 
                                                  !  soil to root direction (kg water/m2 root area/Mpa/s)
   character(kind=c_char,len=param_string_length),parameter,public :: hydr_name_kmax_rsurf1 = "fates_hydr_kmax_rsurf1"  
   
   real(r8),protected,public :: hydr_kmax_rsurf2         !  maximum conducitivity for unit root surface 
                                                  !  root to soil direciton (kg water/m2 root area/Mpa/s)
   character(kind=c_char,len=param_string_length),parameter,public :: hydr_name_kmax_rsurf2 = "fates_hydr_kmax_rsurf2" 

   real(r8),protected,public :: hydr_psi0          !  sapwood water potential at saturation (MPa)
   character(kind=c_char,len=param_string_length),parameter,public :: hydr_name_psi0 = "fates_hydr_psi0"

   real(r8),protected,public :: hydr_psicap        !  sapwood water potential at which capillary reserves exhausted (MPa)
   character(kind=c_char,len=param_string_length),parameter,public :: hydr_name_psicap = "fates_hydr_psicap"

   public :: EDParamsPySet
  
contains

  subroutine EDParamsPySet(rval,name)
    
    implicit none
    ! Arguments
    character(kind=c_char,len=*), intent(in) :: name
    real(r8),intent(in) :: rval
    
    if(trim(name) == trim(hydr_name_psi0))then
       hydr_psi0 = rval
    elseif(trim(name) == trim(hydr_name_psicap))then
       hydr_psicap = rval
    else
       print*,"ERROR in EDParamsPySet, uknown variable name: ",trim(name)
       stop
    end if
    
    return
  end subroutine EDParamsPySet
    


  
end module EDParamsMod
