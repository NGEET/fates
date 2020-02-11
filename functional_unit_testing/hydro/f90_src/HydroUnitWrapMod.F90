
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
   use FatesHydroWTFMod, only : wrf_arr_type,wkf_arr_type,wrf_type_tfs

   implicit none
   public
   save

   integer(kind=c_int), parameter :: param_string_length = 32

   integer, public, parameter :: van_genuchten      = 1
   integer, public, parameter :: campbell           = 2
   integer, public, parameter :: tfs                = 3


   class(wrf_arr_type), public, pointer :: wrfs(:)   ! This holds all (soil and plant) water retention functions
   class(wkf_arr_type), public, pointer :: wkfs(:)   ! 

  
contains

  subroutine InitAllocWTFs(n_wrfs,n_wkfs)

      integer,intent(in) :: n_wrfs
      integer,intent(in) :: n_wkfs
      
      allocate(wrfs(n_wrfs))
      allocate(wkfs(n_wkfs))

      return
  end subroutine InitAllocWTFs

  ! =====================================================================================
  
  subroutine SetWRF(index,itype,npvals,pvals)

      ! The unit testing frameworks don't like assumed shape
      ! array arguments

      integer,intent(in)   :: index
      integer,intent(in)   :: itype
      integer,intent(in)   :: npvals
      real(r8), intent(in) :: pvals(npvals)

      class(wrf_type_vg), pointer :: wrf_vg
      class(wrf_type_cch), pointer :: wrf_cch
      class(wrf_type_tfs), pointer :: wrf_tfs

      print*,"ALLOCATING WRF",index,itype
      print*,pvals

      if(itype == van_genuchten) then
          allocate(wrf_vg)
          wrfs(index)%p => wrf_vg
          call wrf_vg%set_wrf_param(pvals) !alpha,psd,th_sat,th_res
      elseif(itype==campbell) then
          allocate(wrf_cch)
          wrfs(index)%p => wrf_cch
          call wrf_cch%set_wrf_param(pvals)  !th_sat,psi_sat,beta
      else
         allocate(wrf_tfs)
         wrfs(index)%p => wrf_tfs
         call wrf_tfs%set_wrf_param(pvals)
      end if

      return
  end subroutine SetWRF

  subroutine SetWKF(index,itype,npvals,pvals)

      integer,intent(in)   :: index
      integer,intent(in)   :: itype
      integer,intent(in)   :: npvals
      real(r8), intent(in) :: pvals(npvals)

      class(wkf_type_vg), pointer :: wkf_vg
      class(wkf_type_cch), pointer :: wkf_cch
      class(wkf_type_tfs), pointer :: wkf_tfs

      if(itype == van_genuchten) then
          allocate(wkf_vg)
          wkfs(index)%p => wkf_vg
          call wkf_vg%set_wkf_param(pvals) !alpha,psd,th_sat,th_res,tort
       elseif(itype==campbell) then
          allocate(wkf_cch)
          wkfs(index)%p => wkf_cch
          call wkf_cch%set_wkf_param(pvals) !th_sat,psi_sat,beta
      elseif(itype==tfs) then
          allocate(wkf_tfs)
          wkfs(index)%p => wkf_tfs
          call wkf_tfs%set_wkf_param(pvals) !th_sat,p50,avuln
      else
          print*,"UNKNOWN WKF"
          stop
      end if

      return
  end subroutine SetWKF


  function WrapTHFromPSI(index,psi) result(th)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: th

      th = wrfs(index)%p%th_from_psi(psi)

      return
  end function WrapTHFromPSI


  function WrapPSIFromTH(index,th) result(psi)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: th
      real(r8) :: psi

      psi = wrfs(index)%p%psi_from_th(th)

  end function WrapPSIFromTH


  function WrapDPSIDTH(index,th) result(dpsidth)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: th
      real(r8) :: dpsidth

      dpsidth = wrfs(index)%p%dpsidth_from_th(th)

  end function WrapDPSIDTH

  
  function WrapDFTCDPSI(index,psi) result(dftcdpsi)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: dftcdpsi

      dftcdpsi = wkfs(index)%p%dftcdpsi_from_psi(psi)
      
  end function WrapDFTCDPSI


  function WrapFTCFromPSI(index,psi) result(ftc)
      
      integer, intent(in) :: index
      real(r8),intent(in) :: psi
      real(r8) :: ftc

      ftc = wkfs(index)%p%ftc_from_psi(psi)

      return
  end function WrapFTCFromPSI

  
end module HydroUnitWrapMod
