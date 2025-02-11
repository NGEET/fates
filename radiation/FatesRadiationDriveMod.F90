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
  use EDParamsMod,        only : maxpft
  use FatesConstantsMod , only : r8 => fates_r8
  use FatesConstantsMod , only : fates_unset_r8
  use FatesConstantsMod , only : itrue
  use FatesConstantsMod , only : pi_const
  use FatesConstantsMod , only : nocomp_bareground
  use FatesConstantsMod , only : nearzero
  use FatesInterfaceTypesMod , only : bc_in_type
  use FatesInterfaceTypesMod , only : bc_out_type
  use FatesInterfaceTypesMod , only : numpft
  use FatesRadiationMemMod, only : num_rad_stream_types
  use FatesRadiationMemMod, only : idirect, idiffuse
  use FatesRadiationMemMod, only : num_swb, ivis, inir, ipar
  use FatesRadiationMemMod, only : alb_ice, rho_snow, tau_snow
  use FatesRadiationMemMod, only : norman_solver
  use FatesRadiationMemMod, only : twostr_solver
  use EDParamsMod, only          : radiation_model
  use TwoStreamMLPEMod, only : normalized_upper_boundary
  use FatesTwoStreamUtilsMod, only : FatesPatchFSun
  use FatesTwoStreamUtilsMod, only : CheckPatchRadiationBalance
  use FatesInterfaceTypesMod        , only : hlm_hio_ignore_val
  use EDParamsMod        , only : dinc_vai,dlower_vai
  use EDParamsMod        , only : nclmax
  use EDParamsMod        , only : nlevleaf
  use EDCanopyStructureMod, only: calc_areaindex
  use FatesGlobals      , only : fates_log
  use FatesGlobals, only      : endrun => fates_endrun
  use EDPftvarcon,        only : EDPftvarcon_inst
  use FatesNormanRadMod,  only : PatchNormanRadiation
  
  ! CIME globals
  use shr_log_mod       , only : errMsg => shr_log_errMsg

  implicit none

  private
  public :: FatesNormalizedCanopyRadiation  ! Surface albedo and two-stream fluxes
  public :: FatesSunShadeFracs

  logical :: debug = .false.  ! for debugging this module
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  logical :: preserve_b4b = .true.
  
contains

  subroutine FatesNormalizedCanopyRadiation(nsites, sites, bc_in, bc_out )

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

    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)      ! FATES site vector
    type(bc_in_type),   intent(in)            :: bc_in(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)

    ! !LOCAL VARIABLES:
    integer :: s                                   ! site loop counter
    integer :: ifp                                 ! patch loop counter
    integer :: ib                                  ! radiation broad band counter
    type(fates_patch_type), pointer :: currentPatch   ! patch pointer

    !-----------------------------------------------------------------------
    ! -------------------------------------------------------------------------------
    ! TODO (mv, 2014-10-29) the filter here is different than below
    ! this is needed to have the VOC's be bfb - this needs to be
    ! re-examined int he future
    ! RGK,2016-08-06: FATES is still incompatible with VOC emission module
    ! -------------------------------------------------------------------------------


    do s = 1, nsites

       currentpatch => sites(s)%oldest_patch
       do while (associated(currentpatch))

          ! do not do albedo calculations for bare ground patch in SP mode

          ifp = currentpatch%patchno
          
          if_bareground: if(currentpatch%nocomp_pft_label.ne.nocomp_bareground)then
             
             ! Zero diagnostics
             currentPatch%f_sun      (:,:,:) = 0._r8
             currentPatch%fabd_sun_z (:,:,:) = 0._r8
             currentPatch%fabd_sha_z (:,:,:) = 0._r8
             currentPatch%fabi_sun_z (:,:,:) = 0._r8
             currentPatch%fabi_sha_z (:,:,:) = 0._r8
             currentPatch%fabd       (:)     = 0._r8
             currentPatch%fabi       (:)     = 0._r8
             currentPatch%nrmlzd_parprof_pft_dir_z(:,:,:,:) = 0._r8
             currentPatch%nrmlzd_parprof_pft_dif_z(:,:,:,:) = 0._r8
             
             currentPatch%rad_error(:)              = hlm_hio_ignore_val
             
             currentPatch%solar_zenith_flag         = bc_in(s)%filter_vegzen_pa(ifp)
             currentPatch%solar_zenith_angle        = bc_in(s)%coszen_pa(ifp)
             currentPatch%gnd_alb_dif(1:num_swb) = bc_in(s)%albgr_dif_rb(1:num_swb)
             currentPatch%gnd_alb_dir(1:num_swb) = bc_in(s)%albgr_dir_rb(1:num_swb)
             currentPatch%fcansno                   = bc_in(s)%fcansno_pa(ifp)
             
             if(radiation_model.eq.twostr_solver) then
                call currentPatch%twostr%CanopyPrep(bc_in(s)%fcansno_pa(ifp))
                call currentPatch%twostr%ZenithPrep(bc_in(s)%coszen_pa(ifp))
             end if

             if_zenith_flag: if(.not.currentPatch%solar_zenith_flag )then

                ! Sun below horizon, trivial solution
                ! Note (RGK-MLO): Investigate twilight mechanics for
                ! non-zero diffuse radiation when cosz<=0

                ! Temporarily turn off to preserve b4b
                ! preserve_b4b will be removed soon. This is kept here to prevent
                ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23) 
                if (.not.preserve_b4b) then
                   bc_out(s)%albd_parb(ifp,:)            = 1._r8
                   bc_out(s)%albi_parb(ifp,:)            = 1._r8
                   bc_out(s)%fabi_parb(ifp,:)            = 0._r8
                   bc_out(s)%fabd_parb(ifp,:)            = 0._r8
                   bc_out(s)%ftdd_parb(ifp,:)            = 0._r8
                   bc_out(s)%ftid_parb(ifp,:)            = 0._r8
                   bc_out(s)%ftii_parb(ifp,:)            = 0._r8
                end if
             else

                bc_out(s)%albd_parb(ifp,:)            = 0._r8  ! output HLM
                bc_out(s)%albi_parb(ifp,:)            = 0._r8  ! output HLM
                bc_out(s)%fabi_parb(ifp,:)            = 0._r8  ! output HLM
                bc_out(s)%fabd_parb(ifp,:)            = 0._r8  ! output HLM
                bc_out(s)%ftdd_parb(ifp,:)            = 1._r8 ! output HLM
                bc_out(s)%ftid_parb(ifp,:)            = 1._r8 ! output HLM
                bc_out(s)%ftii_parb(ifp,:)            = 1._r8 ! output HLM

                if_nrad: if (maxval(currentPatch%nrad(1,:))==0)then
                   ! there are no leaf layers in this patch. it is effectively bare ground.
                   bc_out(s)%fabd_parb(ifp,:) = 0.0_r8
                   bc_out(s)%fabi_parb(ifp,:) = 0.0_r8
                   currentPatch%rad_error(:)  = 0.0_r8

                   do ib = 1,num_swb
                      bc_out(s)%albd_parb(ifp,ib) = bc_in(s)%albgr_dir_rb(ib)
                      bc_out(s)%albi_parb(ifp,ib) = bc_in(s)%albgr_dif_rb(ib)
                      bc_out(s)%ftdd_parb(ifp,ib) = 1.0_r8
                      bc_out(s)%ftid_parb(ifp,ib) = 0.0_r8
                      bc_out(s)%ftii_parb(ifp,ib) = 1.0_r8
                   enddo

                else

                   select case(radiation_model)
                   case(norman_solver)

                      call PatchNormanRadiation (currentPatch, &
                           bc_out(s)%albd_parb(ifp,:), &   ! Surface Albedo direct
                           bc_out(s)%albi_parb(ifp,:), &   ! Surface Albedo (indirect) diffuse
                           bc_out(s)%fabd_parb(ifp,:), &   ! Fraction direct absorbed by canopy per unit incident
                           bc_out(s)%fabi_parb(ifp,:), &   ! Fraction diffuse absorbed by canopy per unit incident
                           bc_out(s)%ftdd_parb(ifp,:), &   ! Down direct flux below canopy per unit direct at top
                           bc_out(s)%ftid_parb(ifp,:), &   ! Down diffuse flux below canopy per unit direct at top
                           bc_out(s)%ftii_parb(ifp,:))     ! Down diffuse flux below canopy per unit diffuse at top

                   case(twostr_solver)

                        associate( twostr => currentPatch%twostr)

                        !call twostr%CanopyPrep(bc_in(s)%fcansno_pa(ifp))
                        !call twostr%ZenithPrep(bc_in(s)%coszen_pa(ifp))

                        do ib = 1,num_swb

                           twostr%band(ib)%albedo_grnd_diff = bc_in(s)%albgr_dif_rb(ib)
                           twostr%band(ib)%albedo_grnd_beam = bc_in(s)%albgr_dir_rb(ib)

                           call twostr%Solve(ib,             &  ! in
                                normalized_upper_boundary,   &  ! in
                                1.0_r8,1.0_r8,               &  ! in
                                sites(s)%taulambda_2str,        &  ! inout (scratch)
                                sites(s)%omega_2str,            &  ! inout (scratch)
                                sites(s)%ipiv_2str,             &  ! inout (scratch)
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
                      end associate

                   end select
                   
                end if if_nrad

             endif if_zenith_flag
          end if if_bareground

          currentPatch => currentPatch%younger
       end do       ! Loop linked-list patches
    enddo           ! Loop Sites

    return
  end subroutine FatesNormalizedCanopyRadiation

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
             cpatch%parprof_pft_dir_z(:,:,:) = 0._r8
             cpatch%parprof_pft_dif_z(:,:,:) = 0._r8

             bc_out(s)%fsun_pa(ifp) = 0._r8

             ! preserve_b4b will be removed soon. This is kept here to prevent
             ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23)
             if(.not.preserve_b4b)then
                bc_out(s)%laisun_pa(ifp) = 0._r8
                bc_out(s)%laisha_pa(ifp) = calc_areaindex(cpatch,'elai')
             end if

             sunlai  = 0._r8
             shalai  = 0._r8
             if_norm_twostr: if (radiation_model.eq.norman_solver) then
                
                ! Loop over patches to calculate laisun_z and laisha_z for each layer.
                ! Derive canopy laisun, laisha, and fsun from layer sums.
                ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
                ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
                
                ! cpatch%f_sun is calculated in the surface_albedo routine...
                
                do cl = 1, cpatch%ncl_p
                   do ft = 1,numpft
                      ! preserve_b4b will be removed soon. This is kept here to prevent
                      ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23)
                      if(.not.preserve_b4b) then
                         sunlai = sunlai + sum(cpatch%elai_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                              cpatch%f_sun(cl,ft,1:cpatch%nrad(cl,ft)))
                         shalai = shalai + sum(cpatch%elai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
                      else
                         do iv = 1,cpatch%nrad(cl,ft)
                            cpatch%ed_laisun_z(CL,ft,iv) = cpatch%elai_profile(CL,ft,iv) * &
                                 cpatch%f_sun(CL,ft,iv)

                            cpatch%ed_laisha_z(CL,ft,iv) = cpatch%elai_profile(CL,ft,iv) * &
                                 (1._r8 - cpatch%f_sun(CL,ft,iv))
                            
                         end do

                         !needed for the VOC emissions, etc.
                         sunlai = sunlai + sum(cpatch%ed_laisun_z(CL,ft,1:cpatch%nrad(CL,ft)))
                         shalai = shalai + sum(cpatch%ed_laisha_z(CL,ft,1:cpatch%nrad(CL,ft)))
                         
                      end if
                   end do
                end do
                ! preserve_b4b will be removed soon. This is kept here to prevent
                ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23)
                if(.not.preserve_b4b)then
                   shalai = shalai-sunlai
                end if
                
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
                
                ! Convert normalized radiation error units from fraction of radiation to W/m2
                do ib = 1,num_swb
                   cpatch%rad_error(ib) = cpatch%rad_error(ib) * &
                        (bc_in(s)%solad_parb(ifp,ib) + bc_in(s)%solai_parb(ifp,ib))
                end do
                
                ! output the actual PAR profiles through the canopy for diagnostic purposes
                do cl = 1, cpatch%ncl_p
                   do ft = 1,numpft
                      do iv = 1, cpatch%nrad(cl,ft)
                         cpatch%parprof_pft_dir_z(cl,ft,iv) = (bc_in(s)%solad_parb(ifp,ipar) * &
                              cpatch%nrmlzd_parprof_pft_dir_z(idirect,cl,ft,iv)) + &
                              (bc_in(s)%solai_parb(ifp,ipar) * &
                              cpatch%nrmlzd_parprof_pft_dir_z(idiffuse,cl,ft,iv))
                         
                         cpatch%parprof_pft_dif_z(cl,ft,iv) = (bc_in(s)%solad_parb(ifp,ipar) * &
                              cpatch%nrmlzd_parprof_pft_dif_z(idirect,cl,ft,iv)) + &
                              (bc_in(s)%solai_parb(ifp,ipar) * &
                              cpatch%nrmlzd_parprof_pft_dif_z(idiffuse,cl,ft,iv))
                         
                      end do ! iv
                   end do    ! ft
                end do       ! cl
                
             else

                ! If there is no sun out, we have a trivial solution
                if_zenithflag: if(cpatch%solar_zenith_flag ) then
                   
                   ! Two-stream 
                   ! -----------------------------------------------------------
                   do ib = 1,num_swb
                      cpatch%twostr%band(ib)%Rbeam_atm = bc_in(s)%solad_parb(ifp,ib)
                      cpatch%twostr%band(ib)%Rdiff_atm = bc_in(s)%solai_parb(ifp,ib)
                   end do
                   
                   area_vlpfcl(:,:,:) = 0._r8
                   cpatch%f_sun(:,:,:) = 0._r8
                   
                   call FatesPatchFSun(cpatch,    &
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
                              nv = minloc(dlower_vai, DIM=1, MASK=(dlower_vai>vai))
                              do iv = 1, nv
                                 
                                 vai_top = dlower_vai(iv)-dinc_vai(iv)
                                 vai_bot = min(dlower_vai(iv),twostr%scelg(cl,icol)%sai+twostr%scelg(cl,icol)%lai)
                                 
                                 cpatch%parprof_pft_dir_z(cl,ft,iv) = cpatch%parprof_pft_dir_z(cl,ft,iv) + &
                                      area_frac*twostr%GetRb(cl,icol,ivis,vai_top)
                                 cpatch%parprof_pft_dif_z(cl,ft,iv) = cpatch%parprof_pft_dif_z(cl,ft,iv) + &
                                      area_frac*twostr%GetRdDn(cl,icol,ivis,vai_top) + &
                                      area_frac*twostr%GetRdUp(cl,icol,ivis,vai_top)
                                 
                                 call twostr%GetAbsRad(cl,icol,ipar,vai_top,vai_bot, &
                                      Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac)
                                 
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
                              cpatch%parprof_pft_dir_z(cl,ft,iv) = &
                                   cpatch%parprof_pft_dir_z(cl,ft,iv) / area_vlpfcl(iv,ft,cl)
                              cpatch%parprof_pft_dif_z(cl,ft,iv) = &
                                   cpatch%parprof_pft_dif_z(cl,ft,iv) / area_vlpfcl(iv,ft,cl)
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
