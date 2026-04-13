module FatesRadiationDriveMod

  !-------------------------------------------------------------------------------------
  ! EDSurfaceRadiation
  !
  ! This module contains function and type definitions for all things related
  ! to radiative transfer in ED modules at the land surface.
  !
  !-------------------------------------------------------------------------------------

#include "shr_assert.h"

  use EDTypesMod        , only : ed_site_type
  use FatesPatchMod,      only : fates_patch_type
  ! [PORTED by Hui Tang: NVP cohort access for Beer's law radiation]
  use FatesCohortMod,     only : fates_cohort_type
  use EDParamsMod,        only : maxpft
  use EDParamsMod       , only : GetNVegLayers
  use EDParamsMod       , only : nvp_extinction_coeff  ! [PORTED by Hui Tang: NVP Beer's law k from parameter file]
  use FatesConstantsMod , only : r8 => fates_r8
  use FatesConstantsMod , only : fates_unset_r8
  use FatesConstantsMod , only : itrue
  use FatesConstantsMod , only : pi_const
  use FatesConstantsMod , only : nocomp_bareground
  use FatesConstantsMod , only : nearzero
  use FatesInterfaceTypesMod , only : bc_in_type
  use FatesInterfaceTypesMod , only : bc_out_type
  use FatesInterfaceTypesMod , only : numpft
  use FatesInterfaceTypesMod , only : hlm_radiation_model
  ! [PORTED by Hui Tang: NVP control flag and radiation model switch]
  use FatesInterfaceTypesMod , only : hlm_use_nvp
  use FatesInterfaceTypesMod , only : hlm_nvp_rad_model_ground
  use FatesRadiationMemMod, only : num_rad_stream_types
  use FatesRadiationMemMod, only : idirect, idiffuse
  use FatesRadiationMemMod, only : num_swb, ivis, inir, ipar
  use FatesRadiationMemMod, only : alb_ice, rho_snow, tau_snow
  use FatesRadiationMemMod, only : norman_solver
  use FatesRadiationMemMod, only : twostr_solver
  use TwoStreamMLPEMod, only : normalized_upper_boundary
  use FatesTwoStreamUtilsMod, only : FatesPatchFSun
  use FatesTwoStreamUtilsMod, only : CheckPatchRadiationBalance
  use FatesInterfaceTypesMod, only : hlm_hio_ignore_val
  use EDParamsMod        , only : dinc_vai,dlower_vai
  use EDParamsMod        , only : nclmax
  use EDParamsMod        , only : nlevleaf
  use EDCanopyStructureMod, only: calc_areaindex
  use FatesGlobals      , only : fates_log
  use FatesGlobals, only      : endrun => fates_endrun
  use EDPftvarcon,        only : EDPftvarcon_inst
  use FatesNormanRadMod,  only : PatchNormanRadiation
  ! [PORTED by Hui Tang: NVP PFT identification by zero stomatal intercept]
  use FatesLeafBiophysParamsMod, only : lb_params
  
  ! CIME globals
  use shr_log_mod       , only : errMsg => shr_log_errMsg

  implicit none

  private
  public :: FatesNormalizedCanopyRadiation  ! Surface albedo and two-stream fluxes
  public :: FatesSunShadeFracs
  public :: NVPBeerLawAbsorptance           ! Beer's law NVP absorptance (Approach A)

  logical :: debug = .false.  ! for debugging this module
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains

  subroutine FatesNormalizedCanopyRadiation(sites, bc_in, bc_out )

    ! Perform normalized (ie per unit downwelling radiative forcing) radiation
    ! scattering of the vegetation canopy.
    ! This call is normalized because the host wants an albedo for the next time
    ! step, but it does not have the absolute beam and diffuse forcing for the
    ! next step yet.
    ! However, with both Norman and Two stream, we save normalized scattering
    ! and absorption profiles amonst the vegetation, and that can
    ! be scaled by the forcing when we perform diagnostics, calculate heating
    ! rates (HLM side), and calculate absorbed leaf PAR for photosynthesis.

    !

    ! !ARGUMENTS:

    type(ed_site_type), intent(inout), target :: sites(:)      ! FATES site vector
    type(bc_in_type),   intent(in)            :: bc_in(:)
    type(bc_out_type),  intent(inout)         :: bc_out(:)

    ! !LOCAL VARIABLES:
    integer :: s                                   ! site loop counter
    integer :: nsites                              ! number of sites
    integer :: ifp                                 ! patch loop counter
    integer :: ib                                  ! radiation broad band counter
    integer :: cl, iv, icol, ft                    ! indices for canopy layer,leaf layer,
                                                   ! rad column and functional type
    integer :: nv                                  ! number of veg layers
    real(r8) :: area_frac                          ! area fraction for layer of interest
    real(r8) :: vai_top                            ! integrated (top-down) vegetation area
                                                   ! index at lop of layer
    real(r8) :: vai                                ! total VAI of the scattering element
    type(fates_patch_type), pointer :: currentPatch   ! patch pointer
    ! [PORTED by Hui Tang: NVP radiation (R3 ground-albedo modification, R4 Beer's law)]
    ! lai_nvp_pa already computed in EDCanopyStructureMod; only a brief cohort walk is
    ! needed here to find nvp_ft (the NVP PFT index) for the albedo lookup in R3.
    type(fates_cohort_type), pointer :: currentCohort ! cohort pointer (NVP PFT lookup only)
    integer  :: nvp_ft           ! NVP PFT index for EDPftvarcon albedo lookup
    real(r8) :: nvp_frac         ! NVP fractional coverage of patch [0-1]

    !-----------------------------------------------------------------------
    ! -------------------------------------------------------------------------------
    ! TODO (mv, 2014-10-29) the filter here is different than below
    ! this is needed to have the VOC's be bfb - this needs to be
    ! re-examined int he future
    ! RGK,2016-08-06: FATES is still incompatible with VOC emission module
    ! -------------------------------------------------------------------------------
    
    nsites = size(sites,dim=1)

    do s = 1, nsites

       ! Currently holding a copy of this at the site level for restarts
       sites(s)%coszen = bc_in(s)%coszen

       currentpatch => sites(s)%oldest_patch
       do while (associated(currentpatch))

          ifp = currentpatch%patchno
          
          if_bareground: if(currentpatch%nocomp_pft_label.ne.nocomp_bareground)then
             
             ! Initialize output boundary conditions with trivial assumption
             ! This matches CLM/ELM
             ! Albedo is perfect reflector, no flux into or through canopy
             bc_out(s)%albd_parb(ifp,:)            = 1._r8
             bc_out(s)%albi_parb(ifp,:)            = 1._r8
             bc_out(s)%fabi_parb(ifp,:)            = 0._r8
             bc_out(s)%fabd_parb(ifp,:)            = 0._r8
             bc_out(s)%ftdd_parb(ifp,:)            = 0._r8
             bc_out(s)%ftid_parb(ifp,:)            = 0._r8
             bc_out(s)%ftii_parb(ifp,:)            = 0._r8

             ! Zero diagnostics
             currentPatch%f_sun      (:,:,:) = 0._r8
             currentPatch%fabd_sun_z (:,:,:) = 0._r8
             currentPatch%fabd_sha_z (:,:,:) = 0._r8
             currentPatch%fabi_sun_z (:,:,:) = 0._r8
             currentPatch%fabi_sha_z (:,:,:) = 0._r8
             currentPatch%fabd       (:)     = 0._r8
             currentPatch%fabi       (:)     = 0._r8
             currentPatch%nrmlzd_parprof_pft_dir_z(:,:,:) = 0._r8
             currentPatch%nrmlzd_parprof_pft_dif_z(:,:,:) = 0._r8
             currentPatch%gnd_alb_dif(1:num_swb) = bc_in(s)%albgr_dif_rb(1:num_swb)
             currentPatch%gnd_alb_dir(1:num_swb) = bc_in(s)%albgr_dir_rb(1:num_swb)
             currentPatch%fcansno                = bc_in(s)%fcansno_pa(ifp)
             currentPatch%rad_error(:)           = hlm_hio_ignore_val

             ! [PORTED by Hui Tang: Step 8 R3 - modify ground albedo for NVP (no-snow case only)]
             ! Walk cohort list once to find the NVP PFT index for the albedo lookup.
             ! NVP is below snow, so ground albedo is only modified when frac_sno_eff_si == 0.
             nvp_frac = bc_out(s)%nvp_frac_pa(ifp)
             nvp_ft   = 0
             ! [PORTED by Hui Tang: R3 - find nvp_ft whenever NVP is present (needed for
             ! both no-snow ground albedo R3 and SNICAR omega output nvp_omega_pa)]
             if (hlm_use_nvp == itrue .and. nvp_frac > nearzero) then
                currentCohort => currentPatch%shortest
                do while (associated(currentCohort))
                   if (currentCohort%nvp_dz > nearzero) then
                      nvp_ft = currentCohort%pft
                      exit
                   end if
                   currentCohort => currentCohort%taller
                end do
                if (nvp_ft > 0) then
                   ! [PORTED by Hui Tang: set nvp_omega_pa for SNICAR layer-0 (both snow and no-snow)]
                   ! omega(ib) = rhol + taul: fraction of intercepted radiation that is scattered.
                   ! Used in SurfaceAlbedoMod to characterise NVP as a SNICAR pseudo-layer.
                   do ib = 1, num_swb
                      bc_out(s)%nvp_omega_pa(ifp,ib) = EDPftvarcon_inst%rhol(nvp_ft,ib) + &
                                                        EDPftvarcon_inst%taul(nvp_ft,ib)
                   end do
                   ! No-snow case only: modify ground albedo for Norman solver (R3)
                   ! Previous implementation (kept for reference):
                   ! if (bc_in(s)%frac_sno_eff_si <= 0._r8) then
                   !    do ib = 1, num_swb
                   !       currentPatch%gnd_alb_dir(ib) = ...
                   !       currentPatch%gnd_alb_dif(ib) = ...
                   !    end do
                   ! end if
                   ! [PORTED by Hui Tang: R3 applies only for Approach A (NVP as ground boundary)]
                   ! Approach B uses soil albedo as lower boundary; Norman solver handles NVP as
                   ! a canopy leaf layer, so ground albedo must remain the true soil albedo.
                   if (hlm_nvp_rad_model_ground == itrue .and. &
                        bc_in(s)%frac_sno_eff_si <= 0._r8) then
                      do ib = 1, num_swb
                         currentPatch%gnd_alb_dir(ib) = &
                              nvp_frac * EDPftvarcon_inst%rhol(nvp_ft,ib) + &
                              (1._r8 - nvp_frac) * bc_in(s)%albgr_dir_rb(ib)
                         currentPatch%gnd_alb_dif(ib) = &
                              nvp_frac * EDPftvarcon_inst%rhol(nvp_ft,ib) + &
                              (1._r8 - nvp_frac) * bc_in(s)%albgr_dif_rb(ib)
                      end do
                   end if
                end if
             end if
             
             if_zenith_flag: if( bc_in(s)%coszen>0._r8 )then
                
                select case(hlm_radiation_model)
                case(norman_solver)

                   call PatchNormanRadiation (currentPatch, &
                        bc_in(s)%coszen, &
                        bc_out(s)%albd_parb(ifp,:), &   ! Surface Albedo direct
                        bc_out(s)%albi_parb(ifp,:), &   ! Surface Albedo (indirect) diffuse
                        bc_out(s)%fabd_parb(ifp,:), &   ! Fraction direct absorbed by canopy per unit incident
                        bc_out(s)%fabi_parb(ifp,:), &   ! Fraction diffuse absorbed by canopy per unit incident
                        bc_out(s)%ftdd_parb(ifp,:), &   ! Down direct flux below canopy per unit direct at top
                        bc_out(s)%ftid_parb(ifp,:), &   ! Down diffuse flux below canopy per unit direct at top
                        bc_out(s)%ftii_parb(ifp,:))     ! Down diffuse flux below canopy per unit diffuse at top

                case(twostr_solver)

                   associate( twostr => currentPatch%twostr)

                     call twostr%CanopyPrep(currentPatch%fcansno) 
                     call twostr%ZenithPrep(sites(s)%coszen)

                     do ib = 1,num_swb

                        twostr%band(ib)%albedo_grnd_diff = currentPatch%gnd_alb_dif(ib)
                        twostr%band(ib)%albedo_grnd_beam = currentPatch%gnd_alb_dir(ib)

                        call twostr%Solve(ib,             &  ! in
                             normalized_upper_boundary,   &  ! in
                             1.0_r8,1.0_r8,               &  ! in
                             sites(s)%taulambda_2str,         &  ! inout (scratch)
                             sites(s)%omega_2str,             &  ! inout (scratch)
                             sites(s)%ipiv_2str,              &  ! inout (scratch)
                             bc_out(s)%albd_parb(ifp,ib), &  ! out
                             bc_out(s)%albi_parb(ifp,ib), &  ! out
                             currentPatch%rad_error(ib),  &  ! out
                             bc_out(s)%fabd_parb(ifp,ib), &  ! out
                             bc_out(s)%fabi_parb(ifp,ib), &  ! out
                             bc_out(s)%ftdd_parb(ifp,ib), &  ! out
                             bc_out(s)%ftid_parb(ifp,ib), &  ! out
                             bc_out(s)%ftii_parb(ifp,ib))

                        if(debug) then
                           currentPatch%twostr%band(ib)%Rbeam_atm = 1._r8
                           currentPatch%twostr%band(ib)%Rdiff_atm = 1._r8
                           call CheckPatchRadiationBalance(currentPatch, sites(s)%snow_depth, & 
                                ib, bc_out(s)%fabd_parb(ifp,ib),bc_out(s)%fabi_parb(ifp,ib))
                           currentPatch%twostr%band(ib)%Rbeam_atm = fates_unset_r8
                           currentPatch%twostr%band(ib)%Rdiff_atm = fates_unset_r8

                           if(bc_out(s)%fabi_parb(ifp,ib)>1.0 .or. bc_out(s)%fabd_parb(ifp,ib)>1.0)then
                              write(fates_log(),*) 'absorbed fraction > 1.0?'
                              write(fates_log(),*) ifp,ib,bc_out(s)%fabi_parb(ifp,ib),bc_out(s)%fabd_parb(ifp,ib)
                              call twostr%Dump(ib,lat=sites(s)%lat,lon=sites(s)%lon)
                              call endrun(msg=errMsg(sourcefile, __LINE__))
                           end if
                        end if
                     end do

                     ! Fill in the diagnostic arrays for normalized radiation profiles
                     do_cl: do cl = 1,twostr%n_lyr
                        do_icol: do icol = 1,twostr%n_col(cl)
                           ft = twostr%scelg(cl,icol)%pft
                           if_notair: if (ft>0) then
                              area_frac = twostr%scelg(cl,icol)%area
                              vai = twostr%scelg(cl,icol)%sai+twostr%scelg(cl,icol)%lai
                              nv = GetNVegLayers(vai)
                              do iv = 1, nv
                                 vai_top = dlower_vai(iv)
                                 currentPatch%nrmlzd_parprof_pft_dir_z(cl,ft,iv) = currentPatch%nrmlzd_parprof_pft_dir_z(cl,ft,iv) + &
                                      area_frac*twostr%GetRb(cl,icol,ivis,vai_top)
                                 currentPatch%nrmlzd_parprof_pft_dif_z(cl,ft,iv) = currentPatch%nrmlzd_parprof_pft_dif_z(cl,ft,iv) + &
                                      area_frac*twostr%GetRdDn(cl,icol,ivis,vai_top) + &
                                      area_frac*twostr%GetRdUp(cl,icol,ivis,vai_top)
                              end do
                           end if if_notair
                        end do do_icol
                     end do do_cl
                     
                   end associate
                end select

                ! [PORTED by Hui Tang: NVP absorptance for CLM energy balance when there when no-snow]
                ! Approach A: Stand-alone Beer's law on below-vascular-canopy fluxes (trd/tri). NVP not in
                !   Norman's canopy so ftdd/ftii are fluxes incident on NVP from above.
                ! Approach B: NVP IS in Norman's canopy; sum Norman's per-layer absorbed fractions
                !   (fabd_sun_z + fabd_sha_z) for NVP cohort layers. These are per m2 ground,
                !   area-weighted. k_nvp band-independent => apply PAR-band value to all swb bands.
                if (hlm_use_nvp == itrue) then
                   if (hlm_nvp_rad_model_ground == itrue) then
                      ! Approach A: Beer's law absorptance from below-vascular-canopy flux
                      call NVPBeerLawAbsorptance(nvp_frac, bc_out(s)%lai_nvp_pa(ifp), &
                           bc_out(s)%fabd_nvp_pa(ifp,:), bc_out(s)%fabi_nvp_pa(ifp,:))
                   else
                      ! Approach B: sum Norman per-layer output for NVP cohorts
                      bc_out(s)%fabd_nvp_pa(ifp,:) = 0._r8
                      bc_out(s)%fabi_nvp_pa(ifp,:) = 0._r8
                      if (nvp_frac > nearzero) then
                         do cl = 1, currentPatch%ncl_p
                            do ft = 1, numpft
                               if (lb_params%stomatal_intercept(ft) <= nearzero) then
                                  do iv = 1, currentPatch%nrad(cl,ft)
                                     bc_out(s)%fabd_nvp_pa(ifp,:) = bc_out(s)%fabd_nvp_pa(ifp,:) + &
                                          currentPatch%fabd_sun_z(cl,ft,iv) + currentPatch%fabd_sha_z(cl,ft,iv)
                                     bc_out(s)%fabi_nvp_pa(ifp,:) = bc_out(s)%fabi_nvp_pa(ifp,:) + &
                                          currentPatch%fabi_sun_z(cl,ft,iv) + currentPatch%fabi_sha_z(cl,ft,iv)
                                  end do
                               end if
                            end do
                         end do
                      end if
                   end if
                end if

             endif if_zenith_flag
          end if if_bareground
          currentPatch => currentPatch%younger
       end do
    end do
    
    return
  end subroutine FatesNormalizedCanopyRadiation

  ! ======================================================================================

  subroutine NVPBeerLawAbsorptance(nvp_frac, lai_nvp, fabd_nvp, fabi_nvp)
    ! [PORTED by Hui Tang: Beer's law NVP absorptance - extracted from FatesNormalizedCanopyRadiation]
    ! Computes Beer's law absorbed fraction per unit ground area for all radiation bands.
    ! k_nvp is band-independent so fabd_nvp == fabi_nvp.
    ! Used for CLM energy balance (sabg_lyr layer 0) and Approach A photosynthesis PAR.
    implicit none
    real(r8), intent(in)  :: nvp_frac              ! NVP fractional patch coverage [0-1]
    real(r8), intent(in)  :: lai_nvp               ! NVP thallus LAI [m2 thallus / m2 NVP crown]
    real(r8), intent(out) :: fabd_nvp(num_swb)     ! Beer's law direct absorptance per band [-]
    real(r8), intent(out) :: fabi_nvp(num_swb)     ! Beer's law diffuse absorptance per band [-]
    integer :: ib
    if (nvp_frac > nearzero) then
       do ib = 1, num_swb
          fabd_nvp(ib) = nvp_frac * (1._r8 - exp(-nvp_extinction_coeff * lai_nvp))
          fabi_nvp(ib) = fabd_nvp(ib)   ! band-independent: identical to direct
       end do
    else
       fabd_nvp(:) = 0._r8
       fabi_nvp(:) = 0._r8
    end if
  end subroutine NVPBeerLawAbsorptance

  ! ======================================================================================

  subroutine FatesSunShadeFracs(nsites, sites,bc_in,bc_out)

    implicit none

    ! Arguments
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    ! locals
    type (fates_patch_type),pointer :: cpatch   ! c"urrent" patch
    real(r8)          :: sunlai
    real(r8)          :: shalai
    real(r8)          :: elai
    integer           :: cl,ft
    integer           :: iv,ib
    integer           :: s
    integer           :: ifp
    integer           :: nv
    integer           :: icol
    ! Fraction of the canopy area associated with each pft and layer
    ! (used for weighting diagnostics)
    real(r8) :: area_vlpfcl(nlevleaf,maxpft,nclmax) 
    real(r8) :: vai_top,vai_bot
    real(r8) :: area_frac
    real(r8) :: Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac
    real(r8) :: vai
    logical  :: call_fail
    type(fates_patch_type), pointer :: fpatch ! patch pointer for failure reporting
    
    do s = 1,nsites

       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))

          ifp = cpatch%patchno
          
          if_bareground:if(cpatch%nocomp_pft_label.ne.nocomp_bareground)then !only for veg patches

             ! do not do albedo calculations for bare ground patch in SP mode
             
             ! Initialize diagnostics
             cpatch%ed_parsun_z(:,:,:) = 0._r8
             cpatch%ed_parsha_z(:,:,:) = 0._r8
             cpatch%ed_laisun_z(:,:,:) = 0._r8
             cpatch%ed_laisha_z(:,:,:) = 0._r8

             if_norm_twostr: if (hlm_radiation_model.eq.norman_solver) then

                sunlai = 0._r8
                shalai = 0._r8
                
                ! Loop over patches to calculate laisun_z and laisha_z for each layer.
                ! Derive canopy laisun, laisha, and fsun from layer sums.
                ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
                ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
                ! cpatch%f_sun is calculated in the surface_albedo routine...
                
                do cl = 1, cpatch%ncl_p
                   do ft = 1,numpft
                      do iv = 1,cpatch%nrad(cl,ft)
                         cpatch%ed_laisun_z(cl,ft,iv) = cpatch%elai_profile(cl,ft,iv) * &
                              cpatch%f_sun(cl,ft,iv)
                         cpatch%ed_laisha_z(cl,ft,iv) = cpatch%elai_profile(cl,ft,iv) * &
                              (1._r8 - cpatch%f_sun(cl,ft,iv))
                      end do
                      !needed for the VOC emissions, etc.
                      sunlai = sunlai + sum(cpatch%ed_laisun_z(cl,ft,1:cpatch%nrad(cl,ft)))
                      shalai = shalai + sum(cpatch%ed_laisha_z(cl,ft,1:cpatch%nrad(cl,ft)))
                   end do
                end do

                if(sunlai+shalai > 0._r8)then
                   bc_out(s)%fsun_pa(ifp) = sunlai / (sunlai+shalai)
                else
                   bc_out(s)%fsun_pa(ifp) = 0._r8
                endif
                
                if(debug)then
                   if(bc_out(s)%fsun_pa(ifp) > 1._r8)then
                      write(fates_log(),*) 'too much leaf area in profile',  bc_out(s)%fsun_pa(ifp), &
                           sunlai,shalai
                   endif
                end if
                
                elai = calc_areaindex(cpatch,'elai')
                
                bc_out(s)%laisun_pa(ifp) = elai*bc_out(s)%fsun_pa(ifp)
                bc_out(s)%laisha_pa(ifp) = elai*(1.0_r8-bc_out(s)%fsun_pa(ifp))
                
                ! Absorbed PAR profile through canopy
                ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
                ! are canopy integrated so that layer values equal big leaf values.
                
                do cl = 1, cpatch%ncl_p
                   do ft = 1,numpft
                      do iv = 1, cpatch%nrad(cl,ft)
                         
                         cpatch%ed_parsun_z(cl,ft,iv) = &
                              bc_in(s)%solad_parb(ifp,ipar)*cpatch%fabd_sun_z(cl,ft,iv) + &
                              bc_in(s)%solai_parb(ifp,ipar)*cpatch%fabi_sun_z(cl,ft,iv)
                         
                         cpatch%ed_parsha_z(cl,ft,iv) = &
                              bc_in(s)%solad_parb(ifp,ipar)*cpatch%fabd_sha_z(cl,ft,iv) + &
                              bc_in(s)%solai_parb(ifp,ipar)*cpatch%fabi_sha_z(cl,ft,iv)
                         
                      end do !iv
                   end do !ft
                end do !cl
                
                
                
             else  ! if_norm_twostr

                ! If there is no sun out, we have a trivial solution
                if_zenithflag: if( .not. sites(s)%coszen>0._r8 ) then

                   ! Initialize sun/shade fractions for times when zenith is not positive
                   bc_out(s)%laisun_pa(ifp) = 0._r8
                   bc_out(s)%laisha_pa(ifp) = calc_areaindex(cpatch,'elai')
                   bc_out(s)%fsun_pa(ifp)   = 0._r8
                   
                else

                   ! Two-stream 
                   ! -----------------------------------------------------------
                   do ib = 1,num_swb
                      cpatch%twostr%band(ib)%Rbeam_atm = bc_in(s)%solad_parb(ifp,ib)
                      cpatch%twostr%band(ib)%Rdiff_atm = bc_in(s)%solai_parb(ifp,ib)
                   end do
                   
                   area_vlpfcl(:,:,:) = 0._r8
                   cpatch%f_sun(:,:,:) = 0._r8
                   
                   call FatesPatchFSun(sites(s),cpatch,    &
                        bc_out(s)%fsun_pa(ifp),   &
                        bc_out(s)%laisun_pa(ifp), &
                        bc_out(s)%laisha_pa(ifp))
                   
                   associate(twostr => cpatch%twostr)
                     
                     do_cl: do cl = 1,twostr%n_lyr
                        do_icol: do icol = 1,twostr%n_col(cl)
                           
                           ft = twostr%scelg(cl,icol)%pft
                           if_notair: if (ft>0) then
                              area_frac = twostr%scelg(cl,icol)%area
                              vai = twostr%scelg(cl,icol)%sai+twostr%scelg(cl,icol)%lai

                              nv = GetNVegLayers(vai)

                              do iv = 1, nv
                                 
                                 vai_top = dlower_vai(iv)

                                 if(iv == nv) then
                                    vai_bot = twostr%scelg(cl,icol)%sai+twostr%scelg(cl,icol)%lai
                                 else
                                    vai_bot = dlower_vai(iv+1)
                                 end if
                                 
                                 call twostr%GetAbsRad(cl,icol,ipar,vai_top,vai_bot, &
                                      Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac,call_fail)

                                 if(call_fail) then
                                    write(fates_log(),*) 'patch failure:',cpatch%patchno,' of:'
                                    fpatch => sites(s)%oldest_patch
                                    do while (associated(fpatch))
                                       write(fates_log(),*) fpatch%patchno
                                       fpatch => fpatch%younger
                                    end do
                                    call twostr%Dump(ipar,lat=sites(s)%lat,lon=sites(s)%lon)
                                    call endrun(msg=errMsg(sourcefile, __LINE__))
                                 end if
                                 
                                 cpatch%f_sun(cl,ft,iv) = cpatch%f_sun(cl,ft,iv) + &
                                      area_frac*leaf_sun_frac
                                 cpatch%ed_parsun_z(cl,ft,iv) = cpatch%ed_parsun_z(cl,ft,iv) + &
                                      area_frac*(rd_abs_leaf*leaf_sun_frac + rb_abs_leaf)
                                 cpatch%ed_parsha_z(cl,ft,iv) = cpatch%ed_parsha_z(cl,ft,iv) + &
                                      area_frac*rd_abs_leaf*(1._r8-leaf_sun_frac)
                                 
                                 area_vlpfcl(iv,ft,cl) = area_vlpfcl(iv,ft,cl) + area_frac
                              end do
                           end if if_notair
                        end do do_icol
                        
                        do ft = 1,numpft
                           do_iv: do iv = 1,cpatch%nleaf(cl,ft)
                              if(area_vlpfcl(iv,ft,cl)<nearzero) exit do_iv
                              cpatch%f_sun(cl,ft,iv) = cpatch%f_sun(cl,ft,iv) / area_vlpfcl(iv,ft,cl)
                              cpatch%ed_parsun_z(cl,ft,iv) = cpatch%ed_parsun_z(cl,ft,iv) / area_vlpfcl(iv,ft,cl)
                              cpatch%ed_parsha_z(cl,ft,iv) = cpatch%ed_parsha_z(cl,ft,iv) / area_vlpfcl(iv,ft,cl)
                           end do do_iv
                        end do
                        
                     end do do_cl

                   end associate

                end if if_zenithflag
             endif if_norm_twostr
             
          end if if_bareground
          
          cpatch => cpatch%younger
       enddo


    enddo
    return

  end subroutine FatesSunShadeFracs

end module FatesRadiationDriveMod
