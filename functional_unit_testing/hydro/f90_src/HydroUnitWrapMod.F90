
! This module holds hydro specific F90 code needed to run the unit tests
! This is stuff that we don't share with other unit tests, whereas UnitWrapMod.F90 
! can be built generically.

module HydroUnitWrapMod
   !
   ! module that deals with reading the ED parameter file
   !
   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesHydroWTFMod, only : wrf_type,wrf_type_vg,wrf_type_cch
   use FatesHydroWTFMod, only : wkf_type,wkf_type_vg,wkf_type_cch,wkf_type_tfs


   implicit none
   private
   save

   integer(kind=c_int), parameter :: param_string_length = 32

   integer, public, parameter :: van_genuchten      = 1
   integer, public, parameter :: campbell           = 2
   integer, public, parameter :: tfs                = 3


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


    class(wrf_type), pointer :: wrfs(:)   ! This holds all (soil and plant) water retention functions
    class(wkf_type), pointer :: wkfs(:)   ! 


  subroutine InitAllocWTFs(n_wrfs,n_wkfs)


      integer,intent(in) :: n_wrfs
      integer,intent(in) :: n_wkfs
      
      allocate(wrfs(n_wrfs))
      allocate(wkfs(n_wkfs))

      return
  end subroutine InitAllocWTFs


  subroutine SetWRF(index,itype,pvals)

      integer,intent(in)   :: index
      integer,intent(in)   :: itype
      real(r8), intent(in) :: pvals(:)

      class(wrf_type_vg), pointer :: wrf_vg
      class(wrf_type_cch), pointer :: wrf_cch


      if(itype == van_genuchten) then
          allocate(wrf_vg)
          wrfs(index) => wrf_vg
          wrf_vg%set_wrf_param_vg(pvals(1),pvals(2),pvals(3),pvals(4)) !alpha,psd,th_sat,th_res
      elseif(itype==campbell) then
          allocate(wrf_cch)
          wrfs(index) => wrf_cch
          wrf_cch%set_wrf_param_cch(pvals(1),pvals(2),pvals(3))  !th_sat,psi_sat,beta
      else
          print*,"UNKNOWN WRF"
          stop
      end if

      return
  end subroutine SetWRF

  subroutine SetWKF(index,itype,pvals)

      integer,intent(in)   :: index
      integer,intent(in)   :: itype
      real(r8), intent(in) :: pvals(:)

      class(wkf_type_vg), pointer :: wkf_vg
      class(wkf_type_cch), pointer :: wkf_cch


      if(itype == van_genuchten) then
          allocate(wkf_vg)
          wkfs(index) => wkf_vg
          wkf_vg%set_wkf_param_vg(pvals(1),pvals(2),pvals(3),pvals(4),pvals(5)) !alpha,psd,th_sat,th_res,tort
      elseif(itype==campbell) then
          allocate(wkf_cch)
          wkfs(index) => wkf_cch
          wkf_cch%set_wkf_param_cch(pvals(1),pvals(2),pvals(3)) !th_sat,psi_sat,beta
      elseif(itype==tfs) then
          allocate(wkf_tfs)
          wkfs(index) => wkf_tfs
          wkf_tfs%set_wkf_param_tfs(pvals(1),pvals(2),pvals(3)) !th_sat,p50,avuln
      else
          print*,"UNKNOWN WKF"
          stop
      end if

      return
  end subroutine SetWRF


  function WrapTHFromPSI(index,psi) result(th)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: th

      th = wrfs(index)%th_from_psi(psi)

      return
  end function WrapTHFromPSI


  function WrapPSIFromTH(index,th) result(psi)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: th
      real(r8) :: psi

      psi = wrfs(index)%psi_from_th(th)

  end function WrapPSIFromTH


  function WrapDPSIDTH(index,th) result(dpsidth)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: th
      real(r8) :: dpsidth

      dpsidth = wrfs(index)%dpsidth_from_th(th)

  end function WrapDPSIDTH

  
  function WrapDFTCDPSI(index,psi) result(dftcdpsi)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: dftcdpsi

      dftcdth = wrfs(index)%dftcdth_from_psi(psi)
      
  end function WrapDFTCDPSI


  function WrapFTCFromPSI(index,psi) result(ftc)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: ftc

      ftc = wrfs(index)%ftc_from_psi(psi)

      return
  end function WrapFTCFromPSI



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
    


  
end module HydroUnitWrapMod
