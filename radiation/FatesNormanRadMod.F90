module FatesNormanRadMod

  !-------------------------------------------------------------------------------------
  ! EDSurfaceRadiation
  !
  ! This module contains function and type definitions for all things related
  ! to radiative transfer in ED modules at the land surface.
  !
  !-------------------------------------------------------------------------------------

#include "shr_assert.h"

  use EDTypesMod             , only : ed_site_type
  use FatesPatchMod          , only : fates_patch_type
  use EDParamsMod            , only : maxpft
  use FatesConstantsMod      , only : r8 => fates_r8
  use FatesConstantsMod      , only : itrue
  use FatesConstantsMod      , only : pi_const
  use FatesConstantsMod      , only : nocomp_bareground
  use FatesInterfaceTypesMod , only : bc_in_type
  use FatesInterfaceTypesMod , only : bc_out_type
  use FatesInterfaceTypesMod , only : numpft
  use EDParamsMod            , only : nclmax
  use EDParamsMod            , only : nlevleaf
  use FatesRadiationMemMod   , only : num_swb
  use FatesRadiationMemMod   , only : num_rad_stream_types
  use FatesRadiationmemMod   , only : idirect, idiffuse
  use FatesRadiationMemMod   , only : ivis, inir, ipar
  use EDCanopyStructureMod   , only : calc_areaindex
  use FatesGlobals           , only : fates_log
  use FatesGlobals           , only : endrun => fates_endrun
  use EDPftvarcon            , only : EDPftvarcon_inst

  ! CIME globals
  use shr_log_mod       , only : errMsg => shr_log_errMsg

  implicit none

  private
  public :: PatchNormanRadiation

  logical :: debug = .false.  ! for debugging this module
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !  real(r8), public  :: albice(num_swb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
  !       (/ 0.80_r8, 0.55_r8 /)

  !parameters of canopy snow reflectance model.
  ! the parameters in the 2-stream model are not directly analagous to those here
  ! and so they are stored here for now in common with the ice parameters above.
  ! in principle these could be moved to the parameter file.

  real(r8), public  :: albice(num_swb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.80_r8, 0.55_r8 /)
  real(r8), public  :: rho_snow(num_swb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.80_r8, 0.55_r8 /)
  real(r8), public  :: tau_snow(num_swb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.01_r8, 0.01_r8 /)
contains

  subroutine PatchNormanRadiation (currentPatch, &
       albd_parb_out, &   ! (ifp,ib)
       albi_parb_out, &   ! (ifp,ib)
       fabd_parb_out, &   ! (ifp,ib)
       fabi_parb_out, &   ! (ifp,ib)
       ftdd_parb_out, &   ! (ifp,ib)
       ftid_parb_out, &   ! (ifp,ib)
       ftii_parb_out)     ! (ifp,ib)

    ! -----------------------------------------------------------------------------------
    !
    ! This routine performs the Norman Radiation scattering for each patch.
    !
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! !ARGUMENTS:
    ! -----------------------------------------------------------------------------------

    type(fates_patch_type), intent(inout), target :: currentPatch
    real(r8), intent(inout) :: albd_parb_out(num_swb)
    real(r8), intent(inout) :: albi_parb_out(num_swb)
    real(r8), intent(inout) :: fabd_parb_out(num_swb)
    real(r8), intent(inout) :: fabi_parb_out(num_swb)
    real(r8), intent(inout) :: ftdd_parb_out(num_swb)
    real(r8), intent(inout) :: ftid_parb_out(num_swb)
    real(r8), intent(inout) :: ftii_parb_out(num_swb)

    ! Locals
    ! -----------------------------------------------------------------------------------

    integer  :: radtype, L, ft, j
    integer  :: iter                                          ! Iteration index
    integer  :: irep                                          ! Flag to exit iteration loop
    real(r8) :: sb
    real(r8) :: error                                         ! Error check
    real(r8) :: down_rad, up_rad                              ! Iterative solution do Dif_dn and Dif_up
    real(r8) :: ftweight(nclmax,maxpft,nlevleaf)
    real(r8) :: k_dir(maxpft)                              ! Direct beam extinction coefficient
    real(r8) :: tr_dir_z(nclmax,maxpft,nlevleaf)         ! Exponential transmittance of direct beam radiation through a single layer
    real(r8) :: tr_dif_z(nclmax,maxpft,nlevleaf)         ! Exponential transmittance of diffuse radiation through a single layer
    real(r8) :: weighted_dir_tr(nclmax)
    real(r8) :: weighted_fsun(nclmax)
    real(r8) :: weighted_dif_ratio(nclmax,num_swb)
    real(r8) :: weighted_dif_down(nclmax)
    real(r8) :: weighted_dif_up(nclmax)
    real(r8) :: refl_dif(nclmax,maxpft,nlevleaf,num_swb)  ! Term for diffuse radiation reflected by laye
    real(r8) :: tran_dif(nclmax,maxpft,nlevleaf,num_swb)  ! Term for diffuse radiation transmitted by layer
    real(r8) :: dif_ratio(nclmax,maxpft,nlevleaf,num_swb) ! Ratio of upward to forward diffuse fluxes
    real(r8) :: Dif_dn(nclmax,maxpft,nlevleaf)           ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)
    real(r8) :: Dif_up(nclmax,maxpft,nlevleaf)           ! Upward diffuse flux above canopy layer J (W/m**2 ground area)
    real(r8) :: lai_change(nclmax,maxpft,nlevleaf)       ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)

    real(r8) :: frac_lai                                ! Fraction of lai in each layer
    real(r8) :: frac_sai                                ! Fraction of sai in each layer
    real(r8) :: f_abs(nclmax,maxpft,nlevleaf,num_swb)    ! Fraction of light absorbed by surfaces.
    real(r8) :: rho_layer(nclmax,maxpft,nlevleaf,num_swb)! Weighted verage reflectance of layer
    real(r8) :: tau_layer(nclmax,maxpft,nlevleaf,num_swb)! Weighted average transmittance of layer
    real(r8) :: f_abs_leaf(nclmax,maxpft,nlevleaf,num_swb)
    real(r8) :: Abs_dir_z(maxpft,nlevleaf)
    real(r8) :: Abs_dif_z(maxpft,nlevleaf)
    real(r8) :: abs_rad(num_swb)                               !radiation absorbed by soil
    real(r8) :: tr_soili                                      ! Radiation transmitted to the soil surface.
    real(r8) :: tr_soild                                      ! Radiation transmitted to the soil surface.
    real(r8) :: phi1b(maxpft)                                 ! Radiation transmitted to the soil surface.
    real(r8) :: phi2b(maxpft)
    real(r8) :: laisum                                        ! cumulative lai+sai for canopy layer (at middle of layer)
    real(r8) :: angle

    real(r8),parameter :: tolerance = 0.000000001_r8


    integer, parameter :: max_diag_nlevleaf = 4
    integer, parameter :: diag_nlevleaf = min(nlevleaf,max_diag_nlevleaf)  ! for diagnostics, write a small number of leaf layers

    real(r8) :: denom
    real(r8) :: lai_reduction(nclmax)

    integer  :: fp,iv,s      ! array indices
    integer  :: ib               ! waveband number
    real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
    real(r8) :: gdir


    real(r8), parameter :: forc_dir(num_rad_stream_types) = (/ 1.0_r8, 0.0_r8 /)   ! These are binary switches used
    real(r8), parameter :: forc_dif(num_rad_stream_types) = (/ 0.0_r8, 1.0_r8 /)   ! to turn off and on radiation streams



    associate(&
         rhol         =>    EDPftvarcon_inst%rhol                     , & ! Input:  [real(r8) (:)   ] leaf reflectance: 1=vis, 2=nir
         rhos         =>    EDPftvarcon_inst%rhos                     , & ! Input:  [real(r8) (:)   ] stem reflectance: 1=vis, 2=nir
         taul         =>    EDPftvarcon_inst%taul                     , & ! Input:  [real(r8) (:)   ] leaf transmittance: 1=vis, 2=nir
         taus         =>    EDPftvarcon_inst%taus                     , & ! Input:  [real(r8) (:)   ] stem transmittance: 1=vis, 2=nir
         xl           =>    EDPftvarcon_inst%xl                       , & ! Input:  [real(r8) (:)   ] ecophys const - leaf/stem orientation index
         clumping_index  => EDPftvarcon_inst%clumping_index)



      ! Initialize local arrays

    weighted_dir_tr(:)   = 0._r8
    weighted_dif_down(:) = 0._r8
    weighted_dif_up(:)   = 0._r8

    tr_dir_z(:,:,:)      = 0._r8
    tr_dif_z(:,:,:)      = 0._r8
    lai_change(:,:,:)    = 0._r8
    Dif_up(:,:,:)        = 0._r8
    Dif_dn(:,:,:)        = 0._r8
    refl_dif(:,:,:,:)    = 0.0_r8
    tran_dif(:,:,:,:)    = 0.0_r8
    dif_ratio(:,:,:,:)   = 0.0_r8


    ! Initialize the ouput arrays
    ! ---------------------------------------------------------------------------------
    albd_parb_out(1:num_swb) = 0.0_r8
    albi_parb_out(1:num_swb) = 0.0_r8
    fabd_parb_out(1:num_swb) = 0.0_r8
    fabi_parb_out(1:num_swb) = 0.0_r8
    ftdd_parb_out(1:num_swb) = 1.0_r8
    ftid_parb_out(1:num_swb) = 1.0_r8
    ftii_parb_out(1:num_swb) = 1.0_r8

    ! Is this pft/canopy layer combination present in this patch?
    rho_layer(:,:,:,:)=0.0_r8
    tau_layer(:,:,:,:)=0.0_r8
    f_abs(:,:,:,:)=0.0_r8
    f_abs_leaf(:,:,:,:)=0._r8
    do L = 1,currentPatch%ncl_p
       do ft = 1,numpft
          currentPatch%canopy_mask(L,ft) = 0
          do  iv = 1, currentPatch%nrad(L,ft)
             if (currentPatch%canopy_area_profile(L,ft,iv) > 0._r8)then
                currentPatch%canopy_mask(L,ft) = 1

               if(currentPatch%elai_profile(L,ft,iv)+ currentPatch%esai_profile(L,ft,iv).gt.0.0_r8) then
                  frac_lai = currentPatch%elai_profile(L,ft,iv)/&
                       (currentPatch%elai_profile(L,ft,iv)+ currentPatch%esai_profile(L,ft,iv))
               else
                  frac_lai = 1.0_r8
               endif
               !frac_lai = 1.0_r8 ! make the same as previous codebase, in theory.
               frac_sai = 1.0_r8 - frac_lai

               ! layer level reflectance qualities
               do ib = 1,num_swb !vis, nir

                   rho_layer(L,ft,iv,ib)=frac_lai*rhol(ft,ib)+frac_sai*rhos(ft,ib)
                   tau_layer(L,ft,iv,ib)=frac_lai*taul(ft,ib)+frac_sai*taus(ft,ib)

                   ! adjust reflectance and transmittance for canopy snow
                   rho_layer(L,ft,iv,ib)=rho_layer(L,ft,iv,ib)*(1.0_r8- currentPatch%fcansno) &
                        + rho_snow(ib) * currentPatch%fcansno
                   tau_layer(L,ft,iv,ib)=tau_layer(L,ft,iv,ib)*(1.0_r8- currentPatch%fcansno) &
                        + tau_snow(ib) * currentPatch%fcansno

                   ! fraction of incoming light absorbed by leaves or stems.
                   f_abs(L,ft,iv,ib) = 1.0_r8 - tau_layer(L,ft,iv,ib) - rho_layer(L,ft,iv,ib)

                   ! the fraction of the vegetation absorbed light which is absorbed by leaves
                   f_abs_leaf(L,ft,iv,ib) = (1.0_r8- currentPatch%fcansno) * frac_lai* &
                   (1.0_r8 - rhol(ft,ib) - taul(ft,ib))/f_abs(L,ft,iv,ib)

                end do !ib
             endif
          end do !iv
       end do !ft
    end do !L


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Direct beam extinction coefficient, k_dir. PFT specific.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    cosz = max(0.001_r8, currentPatch%solar_zenith_angle ) !copied from previous radiation code...
    do ft = 1,numpft
       sb = (90._r8 - (acos(cosz)*180._r8/pi_const)) * (pi_const / 180._r8)
       phi1b(ft) = 0.5_r8 - 0.633_r8*xl(ft) - 0.330_r8*xl(ft)*xl(ft)
       phi2b(ft) = 0.877_r8 * (1._r8 - 2._r8*phi1b(ft)) !0 = horiz leaves, 1 - vert leaves.
       gdir = phi1b(ft) + phi2b(ft) * sin(sb)
       !how much direct light penetrates a singleunit of lai?
       k_dir(ft) = clumping_index(ft) * gdir / sin(sb)
    end do !FT




    !do this once for one unit of diffuse, and once for one unit of direct radiation
    do radtype = 1, num_rad_stream_types

       ! Extract information that needs to be provided by ED into local array.
       ! RGK: NOT SURE WHY WE NEED FTWEIGHT ...
       ! ------------------------------------------------------------------------------

       ftweight(:,:,:) = 0._r8
       do L = 1,currentPatch%NCL_p
          do ft = 1,numpft
             do  iv = 1, currentPatch%nrad(L,ft)
                !this is already corrected for area in CLAP
                ftweight(L,ft,iv) = currentPatch%canopy_area_profile(L,ft,iv)
             end do  !iv
          end do  !ft1
       end do  !L

       if(debug)then
          if (sum(ftweight(1,:,1))<0.999_r8)then
             write(fates_log(),*) 'canopy not full',ftweight(1,:,1)
          endif
          if (sum(ftweight(1,:,1))>1.0001_r8)then
             write(fates_log(),*) 'canopy too full',ftweight(1,:,1)
          endif
       end if
       
       do L = 1,currentPatch%NCL_p !start at the top canopy layer (1 is the top layer.)

          weighted_dir_tr(L)                 = 0.0_r8
          weighted_fsun(L)                   = 0._r8
          weighted_dif_ratio(L,1:num_swb) = 0._r8

          !Each canopy layer (canopy, understorey) has multiple 'parallel' pft's

          do ft =1,numpft

             if (currentPatch%canopy_mask(L,ft) == 1)then !only do calculation if there are the appropriate leaves.
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                ! Diffuse transmittance, tr_dif, do each layer with thickness elai_z.
                ! Estimated do nine sky angles in increments of 10 degrees
                ! PFT specific...
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                tr_dif_z(L,ft,:) = 0._r8
                do iv = 1,currentPatch%nrad(L,ft)
                   do j = 1,9
                      angle = (5._r8 + real(j - 1,r8) * 10._r8) * pi_const / 180._r8
                      gdir = phi1b(ft) + phi2b(ft) * sin(angle)
                      tr_dif_z(L,ft,iv) = tr_dif_z(L,ft,iv) + exp(-clumping_index(ft) * &
                           gdir / sin(angle) * &
                           (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv))) * &
                           sin(angle)*cos(angle)
                   end do

                   tr_dif_z(L,ft,iv) = tr_dif_z(L,ft,iv) * 2._r8 * (10._r8 * pi_const / 180._r8)

                end do


                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                ! Direct beam transmittance, tr_dir_z, uses cumulative LAI above layer J to give
                ! unscattered direct beam onto layer J. do each PFT section.
                ! This is just an  decay curve based on k_dir. (leaf & sun angle)
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                if (L==1)then
                   tr_dir_z(L,ft,1) = 1._r8
                else
                   tr_dir_z(L,ft,1)  = weighted_dir_tr(L-1)
                endif
                laisum = 0.00_r8
                !total direct beam getting to the bottom of the top canopy.
                do iv = 1,currentPatch%nrad(L,ft)
                   laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                   lai_change(L,ft,iv) = 0.0_r8
                   if (( ftweight(L,ft,iv+1)  >  0.0_r8 ) .and. ( ftweight(L,ft,iv+1)  <  ftweight(L,ft,iv) ))then
                      !where there is a partly empty leaf layer, some fluxes go straight through.
                      lai_change(L,ft,iv) = ftweight(L,ft,iv)-ftweight(L,ft,iv+1)
                   endif
                   if(debug)then
                      if (ftweight(L,ft,iv+1) - ftweight(L,ft,iv) > 1.e-10_r8)then
                         write(fates_log(),*) 'lower layer has more coverage. This is wrong' , &
                              ftweight(L,ft,iv),ftweight(L,ft,iv+1),ftweight(L,ft,iv+1)-ftweight(L,ft,iv)
                      endif
                   end if
                   
                   !n.b. in theory lai_change could be calculated daily in the ED code.
                   !This is light coming striaght through the canopy.
                   if (L==1)then
                      tr_dir_z(L,ft,iv+1) = exp(-k_dir(ft) * laisum)* &
                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                   else
                      tr_dir_z(L,ft,iv+1) = weighted_dir_tr(L-1)*exp(-k_dir(ft) * laisum)* &
                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                   endif

                   if (iv == 1)then
                      !this is the top layer.
                      tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                           ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                   else
                      !the lai_change(iv) affects the light incident on layer iv+2 not iv+1
                      ! light coming from the layer above (iv-1) goes through iv and onto iv+1.
                      if (lai_change(L,ft,iv-1) > 0.0_r8)then
                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv)* &
                              lai_change(L,ft,iv-1) / ftweight(L,ft,1)
                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv-1)* &
                              (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                      else
                         !account fot the light that comes striaght down from unfilled layers above.
                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                              ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                      endif
                   endif

                end do

                !add up all the weighted contributions from the different PFT columns.
                weighted_dir_tr(L) = weighted_dir_tr(L) + tr_dir_z(L,ft,currentPatch%nrad(L,ft)+1)*ftweight(L,ft,1)

                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                ! Sunlit and shaded fraction of leaf layer
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

                !laisum = 0._r8
                do iv = 1,currentPatch%nrad(L,ft)
                   ! Cumulative leaf area. Original code uses cumulative lai do layer.
                   ! Now use cumulative lai at center of layer.
                   ! Same as tr_dir_z calcualtions, but in the middle of the layer? FIX(RF,032414)-WHY?
                   if (iv  ==  1) then
                      laisum = 0.5_r8 * (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv))
                   else
                      laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                   end if


                   if (L == 1)then !top canopy layer
                      currentPatch%f_sun(L,ft,iv) = exp(-k_dir(ft) * laisum)* &
                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                   else
                      currentPatch%f_sun(L,ft,iv) = weighted_fsun(L-1)* exp(-k_dir(ft) * laisum)* &
                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                   endif

                   if ( iv > 1 ) then  ! becasue we are looking at this layer (not the next)
                      ! we only ever add fluxes if iv>1
                      if (lai_change(L,ft,iv-1) > 0.0_r8)then
                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                              currentPatch%f_sun(L,ft,iv) * &
                              lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                              currentPatch%f_sun(L,ft,iv-1) * &
                              (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                      else
                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                              currentPatch%f_sun(L,ft,iv-1) * &
                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                      endif
                   endif

                end do !iv

                weighted_fsun(L) = weighted_fsun(L) + currentPatch%f_sun(L,ft,currentPatch%nrad(L,ft))* &
                     ftweight(L,ft,1)

                ! instance where the first layer ftweight is used a proxy for the whole column. FTWA
                ! this is possibly a source of slight error. If we use the ftweight at the top of the PFT column,
                ! then we willl underestimate fsun, but if we use ftweight at the bottom of the column, we will
                ! underestimate it. Really, we should be tracking the release of direct light from the column as it tapers
                ! towards the ground. Is that necessary to get energy closure? It would be quite hard...
             endif !present.
          end do!pft loop
       end do !L


       do L = currentPatch%NCL_p,1, -1 !start at the bottom and work up.
          do ft = 1,numpft
             if (currentPatch%canopy_mask(L,ft) == 1)then

                !==============================================================================!
                ! Iterative solution do scattering
                !==============================================================================!

                do ib = 1,num_swb !vis, nir
                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                   ! Leaf scattering coefficient and terms do diffuse radiation reflected
                   ! and transmitted by a layer
                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

                   do iv = 1,currentPatch%nrad(L,ft)
                      !How much diffuse light is intercepted and then reflected?
                      refl_dif(L,ft,iv,ib) = (1._r8 - tr_dif_z(L,ft,iv)) * rho_layer(L,ft,iv,ib)
                      !How much diffuse light in this layer is transmitted?
                      tran_dif(L,ft,iv,ib) = (1._r8 - tr_dif_z(L,ft,iv)) * &
                           tau_layer(L,ft,iv,ib) + tr_dif_z(L,ft,iv)
                   end do

                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                   ! Ratio of upward to forward diffuse fluxes, dif_ratio
                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                   ! Soil diffuse reflectance (ratio of down to up radiation).
                   iv = currentPatch%nrad(L,ft) + 1
                   if (L  == currentPatch%NCL_p)then !nearest the soil
                      dif_ratio(L,ft,iv,ib) = currentPatch%gnd_alb_dif(ib)  !bc_in(s)%albgr_dif_rb(ib)
                   else
                      dif_ratio(L,ft,iv,ib) = weighted_dif_ratio(L+1,ib)
                   end if
                   ! Canopy layers, working upwardfrom soil with dif_ratio(iv+1) known
                   ! FIX(RF,032414) ray tracing eqution - need to find derivation of this...
                   ! for each unit going down, there are x units going up.
                   do iv = currentPatch%nrad(L,ft),1, -1
                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv+1,ib) * &
                           tran_dif(L,ft,iv,ib)*tran_dif(L,ft,iv,ib) / &
                           (1._r8 - dif_ratio(L,ft,iv+1,ib) * refl_dif(L,ft,iv,ib)) &
                           + refl_dif(L,ft,iv,ib)
                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) * &
                           ftweight(L,ft,iv)/ftweight(L,ft,1)
                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) + dif_ratio(L,ft,iv+1,ib) * &
                           (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                   end do
                   weighted_dif_ratio(L,ib) = weighted_dif_ratio(L,ib) + &
                        dif_ratio(L,ft,1,ib) * ftweight(L,ft,1)
                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                end do!num_swb
             endif ! currentPatch%canopy_mask
          end do!ft
       end do!L

       do ib = 1,num_swb

          currentPatch%rad_error(ib) = 0._r8
          
          Dif_dn(:,:,:) = 0.00_r8
          Dif_up(:,:,:) = 0.00_r8
          do L = 1, currentPatch%NCL_p !work down from the top of the canopy.
             weighted_dif_down(L) = 0._r8
             do ft = 1, numpft
                if (currentPatch%canopy_mask(L,ft) == 1)then
                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                   ! First estimates do downward and upward diffuse flux
                   !
                   ! Dif_dn =  forward diffuse flux onto layer J
                   ! Dif_up =  Upward diffuse flux above layer J
                   !
                   ! Solved here without direct beam radiation and using dif_ratio = Dif_up / Dif_dn
                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                   ! downward diffuse flux onto the top surface of the canopy

                   if (L == 1)then
                      Dif_dn(L,ft,1) = forc_dif(radtype)
                   else
                      Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                   end if
                   ! forward diffuse flux within the canopy and at soil, working forward through canopy
                   do iv = 1,currentPatch%nrad(L,ft)
                      denom = refl_dif(L,ft,iv,ib) *  dif_ratio(L,ft,iv,ib)
                      denom = 1._r8 - denom
                      Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) / &
                           denom *ftweight(L,ft,iv)/ftweight(L,ft,1)
                      if (iv > 1)then
                         if (lai_change(L,ft,iv-1) > 0.0_r8)then
                            !here we are thinking about whether the layer above had an laichange,
                            !but calculating the flux onto the layer below.
                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1)+ Dif_dn(L,ft,iv)* &
                                 lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1)+ Dif_dn(L,ft,iv-1)* &
                                 (ftweight(L,ft,1)-ftweight(L,ft,iv-1)/ftweight(L,ft,1))
                         else
                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1) + Dif_dn(L,ft,iv) * &
                                 (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                         endif
                      else
                         Dif_dn(L,ft,iv+1)    = Dif_dn(L,ft,iv+1) + Dif_dn(L,ft,iv) * &
                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                      endif
                   end do

                   weighted_dif_down(L) = weighted_dif_down(L) + Dif_dn(L,ft,currentPatch%nrad(L,ft)+1) * &
                        ftweight(L,ft,1)

                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                endif !present
             end do !ft
             if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is the the (incomplete) understorey?
                !Add on the radiation going through the canopy gaps.
                weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1)*(1.0-sum(ftweight(L,:,1)))
                !instance where the first layer ftweight is used a proxy for the whole column. FTWA
             endif
          end do !L

          do L = currentPatch%NCL_p,1 ,-1 !work up from the bottom.
             weighted_dif_up(L) = 0._r8
             do ft = 1, numpft
                if (currentPatch%canopy_mask(L,ft) == 1)then
                   !Bounce diffuse radiation off soil surface.
                   iv = currentPatch%nrad(L,ft) + 1
                   if (L==currentPatch%NCL_p)then !is this the bottom layer ?
                      Dif_up(L,ft,iv) = currentPatch%gnd_alb_dif(ib) * Dif_dn(L,ft,iv)
                   else
                      Dif_up(L,ft,iv) = weighted_dif_up(L+1)
                   end if
                   ! Upward diffuse flux within the canopy and above the canopy, working upward through canopy

                   do iv = currentPatch%nrad(L,ft), 1, -1
                      if (lai_change(L,ft,iv) > 0.0_r8)then
                         Dif_up(L,ft,iv) = dif_ratio(L,ft,iv,ib) * Dif_dn(L,ft,iv) * &
                              ftweight(L,ft,iv) / ftweight(L,ft,1)
                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * &
                              tran_dif(L,ft,iv,ib) * lai_change(L,ft,iv)/ftweight(L,ft,1)
                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * &
                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                         !nb is this the right constuction?
                         ! the radiation that hits the empty space is not reflected.
                      else
                         Dif_up(L,ft,iv) = dif_ratio(L,ft,iv,ib) * Dif_dn(L,ft,iv) * ftweight(L,ft,iv)
                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * (1.0_r8-ftweight(L,ft,iv))
                      endif
                   end do

                   weighted_dif_up(L) = weighted_dif_up(L) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                endif !present
             end do !ft
             if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                !Add on the radiation coming up through the canopy gaps.
                !diffuse to diffuse
                weighted_dif_up(L) = weighted_dif_up(L) +(1.0_r8-sum(ftweight(L,1:numpft,1))) * &
                     weighted_dif_down(L-1) * currentPatch%gnd_alb_dif(ib)
                !direct to diffuse
                weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(radtype) * &
                     weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1))) * currentPatch%gnd_alb_dir(ib)
             endif
          end do !L

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! 3. Iterative calculation of forward and upward diffuse fluxes, iNCL_puding
          ! scattered direct beam
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

          ! Flag to exit iteration loop: 0 = exit and 1 = iterate
          irep = 1
          ! Iteration loop
          iter = 0
          do while(irep ==1 .and. iter<50)

             iter = iter + 1
             irep = 0
             do L = 1,currentPatch%NCL_p !working from the top down
                weighted_dif_down(L) = 0._r8
                do ft =1,numpft
                   if (currentPatch%canopy_mask(L,ft) == 1)then
                      ! forward diffuse flux within the canopy and at soil, working forward through canopy
                      ! with Dif_up -from previous iteration-. Dif_dn(1) is the forward diffuse flux onto the canopy.
                      ! Note: down = forward flux onto next layer
                      if (L == 1)then !is this the top layer?
                         Dif_dn(L,ft,1) = forc_dif(radtype)
                      else
                         Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                      end if
                      down_rad = 0._r8

                      do iv = 1, currentPatch%nrad(L,ft)
                         ! down rad'n is the sum of the down and upwards reflected diffuse fluxes...
                         down_rad = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) + &
                              Dif_up(L,ft,iv+1) * refl_dif(L,ft,iv,ib)

                         !... plus the direct beam intercepted and intransmitted by this layer.
                         down_rad = down_rad + forc_dir(radtype) * tr_dir_z(L,ft,iv) * (1.00_r8 - &
                              exp(-k_dir(ft) *  (currentPatch%elai_profile(L,ft,iv)+ &
                              currentPatch%esai_profile(L,ft,iv))  )) * tau_layer(L,ft,iv,ib)


                         !... plus the direct beam intercepted and intransmitted by this layer.
                         ! modified to spread it out over the whole of incomplete layers.

                         down_rad = down_rad *(ftweight(L,ft,iv)/ftweight(L,ft,1))

                         if (iv > 1)then
                            if (lai_change(L,ft,iv-1) > 0.0_r8)then
                               down_rad = down_rad + Dif_dn(L,ft,iv)   * lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                               down_rad = down_rad + Dif_dn(L,ft,iv-1) * (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ &
                                    ftweight(L,ft,1)
                            else
                               down_rad = down_rad + Dif_dn(L,ft,iv)   * (ftweight(L,ft,1)-ftweight(L,ft,iv))/ &
                                    ftweight(L,ft,1)
                            endif
                         else
                            down_rad = down_rad + Dif_dn(L,ft,iv)   * (ftweight(L,ft,1)-ftweight(L,ft,iv))/ &
                                 ftweight(L,ft,1)
                         endif

                         !this is just Dif down, plus refl up, plus dir intercepted and turned into dif... ,
                         if (abs(down_rad - Dif_dn(L,ft,iv+1)) > tolerance)then
                            irep = 1
                         end if
                         Dif_dn(L,ft,iv+1) = down_rad

                      end do !iv

                      weighted_dif_down(L) = weighted_dif_down(L) + Dif_dn(L,ft,currentPatch%nrad(L,ft)+1) * &
                           ftweight(L,ft,1)

                   endif !present
                end do!ft
                if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                   weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1) * &
                        (1.0_r8-sum(ftweight(L,1:numpft,1)))
                end if
             end do ! do L loop

             do L = 1, currentPatch%NCL_p ! working from the top down.
                weighted_dif_up(L) = 0._r8
                do ft =1,numpft
                   if (currentPatch%canopy_mask(L,ft) == 1)then
                      ! Upward diffuse flux at soil or from lower canopy (forward diffuse and unscattered direct beam)
                      iv = currentPatch%nrad(L,ft) + 1
                      if (L==currentPatch%NCL_p)then  !In the bottom canopy layer, reflect off the soil
                         Dif_up(L,ft,iv) = Dif_dn(L,ft,iv) * currentPatch%gnd_alb_dif(ib) + &
                              forc_dir(radtype) * tr_dir_z(L,ft,iv) * currentPatch%gnd_alb_dir(ib)
                      else      !In the other canopy layers, reflect off the underlying vegetation.
                         Dif_up(L,ft,iv) =  weighted_dif_up(L+1)
                      end if

                      ! Upward diffuse flux within and above the canopy, working upward through canopy
                      ! with Dif_dn from previous interation.  Note: up = upward flux above current layer
                      do iv = currentPatch%nrad(L,ft),1,-1
                         !this is radiation up, by layer transmittance, by

                         !reflection of the lower layer,
                         up_rad = Dif_dn(L,ft,iv) * refl_dif(L,ft,iv,ib)
                         up_rad = up_rad + forc_dir(radtype) * tr_dir_z(L,ft,iv) * (1.00_r8 - exp(-k_dir(ft) * &
                              (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv))))* &
                              rho_layer(L,ft,iv,ib)
                         up_rad = up_rad + Dif_up(L,ft,iv+1) * tran_dif(L,ft,iv,ib)
                         up_rad = up_rad * ftweight(L,ft,iv)/ftweight(L,ft,1)
                         up_rad = up_rad + Dif_up(L,ft,iv+1) *(ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                         ! THE LOWER LAYER FLUX IS HOMOGENIZED, SO WE DON"T CONSIDER THE LAI_CHANGE HERE...

                         if (abs(up_rad - Dif_up(L,ft,iv)) > tolerance) then !are we close to the tolerance level?
                            irep = 1
                         end if
                         Dif_up(L,ft,iv) = up_rad

                      end do  !iv
                      weighted_dif_up(L) = weighted_dif_up(L) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                   end if !present
                end do!ft

                if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then  !is this the (incomplete) understorey?
                   !Add on the radiation coming up through the canopy gaps.
                   weighted_dif_up(L) = weighted_dif_up(L) +(1.0_r8-sum(ftweight(L,1:numpft,1))) * &
                        weighted_dif_down(L-1) * currentPatch%gnd_alb_dif(ib)
                   weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(radtype) * &
                        weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))*currentPatch%gnd_alb_dir(ib)
                end if
             end do!L
          end do ! do while over iter

          abs_rad(ib) = 0._r8
          tr_soili = 0._r8
          tr_soild = 0._r8

          do L = 1, currentPatch%NCL_p !working from the top down.
             abs_dir_z(:,:) = 0._r8
             abs_dif_z(:,:) = 0._r8
             do ft =1,numpft
                if (currentPatch%canopy_mask(L,ft) == 1)then
                   !==============================================================================!
                   ! Compute absorbed flux densities
                   !==============================================================================!

                   ! Absorbed direct beam and diffuse do leaf layers
                   do iv = 1, currentPatch%nrad(L,ft)
                      Abs_dir_z(ft,iv) = ftweight(L,ft,iv)* forc_dir(radtype) * tr_dir_z(L,ft,iv) * &
                           (1.00_r8 - exp(-k_dir(ft) *  (currentPatch%elai_profile(L,ft,iv)+ &
                           currentPatch%esai_profile(L,ft,iv)) )) * f_abs(L,ft,iv,ib)
                      Abs_dif_z(ft,iv) = ftweight(L,ft,iv)* ((Dif_dn(L,ft,iv) + &
                           Dif_up(L,ft,iv+1)) * (1.00_r8 - tr_dif_z(L,ft,iv)) * f_abs(L,ft,iv,ib))
                   end do

                   ! Absorbed direct beam and diffuse do soil
                   if (L == currentPatch%NCL_p)then
                      iv = currentPatch%nrad(L,ft) + 1
                      Abs_dif_z(ft,iv) = ftweight(L,ft,1)*Dif_dn(L,ft,iv) * (1.0_r8 - currentPatch%gnd_alb_dif(ib) )
                      Abs_dir_z(ft,iv) = ftweight(L,ft,1)*forc_dir(radtype) * &
                           tr_dir_z(L,ft,iv) * (1.0_r8 - currentPatch%gnd_alb_dir(ib)  )
                      tr_soild = tr_soild + ftweight(L,ft,1)*forc_dir(radtype) * tr_dir_z(L,ft,iv)
                      tr_soili = tr_soili + ftweight(L,ft,1)*Dif_dn(L,ft,iv)
                   end if

                   ! Absorbed radiation, shaded and sunlit portions of leaf layers
                   !here we get one unit of diffuse radiation... how much of
                   !it is absorbed?
                   if (ib == ivis) then ! only set the absorbed PAR for the visible light band.
                      do iv = 1, currentPatch%nrad(L,ft)
                         if (radtype==idirect) then
                            if ( debug ) then
                               write(fates_log(),*) 'EDsurfAlb 730 ',Abs_dif_z(ft,iv),currentPatch%f_sun(L,ft,iv)
                               write(fates_log(),*) 'EDsurfAlb 731 ', currentPatch%fabd_sha_z(L,ft,iv), &
                                    currentPatch%fabd_sun_z(L,ft,iv)
                            endif
                            currentPatch%fabd_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                 (1._r8 - currentPatch%f_sun(L,ft,iv))*f_abs_leaf(L,ft,iv,ib)
                            currentPatch%fabd_sun_z(L,ft,iv) =( Abs_dif_z(ft,iv) * &
                                 currentPatch%f_sun(L,ft,iv) + &
                                 Abs_dir_z(ft,iv))*f_abs_leaf(L,ft,iv,ib)
                         else
                            currentPatch%fabi_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                 (1._r8 - currentPatch%f_sun(L,ft,iv))*f_abs_leaf(L,ft,iv,ib)
                            currentPatch%fabi_sun_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                 currentPatch%f_sun(L,ft,iv)*f_abs_leaf(L,ft,iv,ib)
                         endif
                         if ( debug ) then
                            write(fates_log(),*) 'EDsurfAlb 740 ', currentPatch%fabd_sha_z(L,ft,iv), &
                                 currentPatch%fabd_sun_z(L,ft,iv)
                         endif
                      end do
                   endif ! ib


                   !==============================================================================!
                   ! Sum fluxes
                   !==============================================================================!
                   ! Solar radiation absorbed by ground
                   iv = currentPatch%nrad(L,ft) + 1
                   if (L==currentPatch%NCL_p)then
                      abs_rad(ib) = abs_rad(ib) +  (Abs_dir_z(ft,iv) + Abs_dif_z(ft,iv))
                   end if
                   ! Solar radiation absorbed by vegetation and sunlit/shaded leaves
                   do iv = 1,currentPatch%nrad(L,ft)
                      if (radtype == idirect)then
                         currentPatch%fabd(ib) = currentPatch%fabd(ib) + &
                              Abs_dir_z(ft,iv)+Abs_dif_z(ft,iv)
                         ! bc_out(s)%fabd_parb_out(ib) = currentPatch%fabd(ib)
                      else
                         currentPatch%fabi(ib) = currentPatch%fabi(ib) + Abs_dif_z(ft,iv)
                         ! bc_out(s)%fabi_parb_out(ib) = currentPatch%fabi(ib)
                      endif
                   end do

                   ! Albefor
                   if (L==1)then !top canopy layer.
                      if (radtype == idirect)then
                         albd_parb_out(ib) = albd_parb_out(ib) + &
                              Dif_up(L,ft,1) * ftweight(L,ft,1)
                      else
                         albi_parb_out(ib) = albi_parb_out(ib) + &
                              Dif_up(L,ft,1) * ftweight(L,ft,1)
                      end if
                   end if

                   ! pass normalized PAR profiles for use in diagnostic averaging for history fields
                   if (ib == ivis) then ! only diagnose PAR profiles for the visible band
                      do iv = 1, currentPatch%nrad(L,ft)
                         currentPatch%nrmlzd_parprof_pft_dir_z(radtype,L,ft,iv) = &
                              forc_dir(radtype) * tr_dir_z(L,ft,iv)
                         currentPatch%nrmlzd_parprof_pft_dif_z(radtype,L,ft,iv) = &
                              Dif_dn(L,ft,iv) + Dif_up(L,ft,iv)
                      end do
                   end if ! ib = visible
                end if ! present
             end do !ft
             if (radtype == idirect)then
                fabd_parb_out(ib) = currentPatch%fabd(ib)
             else
                fabi_parb_out(ib) = currentPatch%fabi(ib)
             endif


             !radiation absorbed from fluxes through unfilled part of lower canopy.
             if (currentPatch%NCL_p > 1.and.L == currentPatch%NCL_p)then
                abs_rad(ib) = abs_rad(ib) + weighted_dif_down(L-1) * &
                     (1.0_r8-sum(ftweight(L,1:numpft,1)))*(1.0_r8-currentPatch%gnd_alb_dif(ib) )
                abs_rad(ib) = abs_rad(ib) + forc_dir(radtype) * weighted_dir_tr(L-1) * &
                     (1.0_r8-sum(ftweight(L,1:numpft,1)))*(1.0_r8-currentPatch%gnd_alb_dir(ib) )
                tr_soili = tr_soili + weighted_dif_down(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))
                tr_soild = tr_soild + forc_dir(radtype) * weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))
             endif

             if (radtype == idirect)then
                currentPatch%tr_soil_dir(ib) = tr_soild
                currentPatch%tr_soil_dir_dif(ib) = tr_soili
                currentPatch%sabs_dir(ib)     = abs_rad(ib)
                ftdd_parb_out(ib)  = tr_soild
                ftid_parb_out(ib) =  tr_soili
             else
                currentPatch%tr_soil_dif(ib) = tr_soili
                currentPatch%sabs_dif(ib)     = abs_rad(ib)
                ftii_parb_out(ib) =  tr_soili
             end if

          end do!l


          !==============================================================================!
          ! Conservation check
          !==============================================================================!
          ! Total radiation balance: absorbed = incoming - outgoing

          if (radtype == idirect)then
             error = abs(currentPatch%sabs_dir(ib) - (currentPatch%tr_soil_dir(ib) * &
                  (1.0_r8-currentPatch%gnd_alb_dir(ib) ) + &
                  currentPatch%tr_soil_dir_dif(ib) * (1.0_r8-currentPatch%gnd_alb_dif(ib)     )))

             if(debug)then
                if ( abs(error) > 0.0001)then
                   write(fates_log(),*)'dir ground absorption error',error,currentPatch%sabs_dir(ib), &
                        currentPatch%tr_soil_dir(ib)* &
                        (1.0_r8-currentPatch%gnd_alb_dir(ib)  ),currentPatch%NCL_p,ib,sum(ftweight(1,1:numpft,1))
                   write(fates_log(),*) 'albedos',currentPatch%sabs_dir(ib) ,currentPatch%tr_soil_dir(ib), &
                        (1.0_r8-currentPatch%gnd_alb_dir(ib)  )
                   do ft =1,numpft
                      iv = currentPatch%nrad(1,ft) + 1
                      write(fates_log(),*) 'abs soil fluxes', Abs_dir_z(ft,iv),Abs_dif_z(ft,iv)
                   end do
                end if
             end if
             
          else
             if (debug) then
                if ( abs(currentPatch%sabs_dif(ib)-(currentPatch%tr_soil_dif(ib) * &
                     (1.0_r8-currentPatch%gnd_alb_dif(ib)  ))) > 0.0001_r8)then
                   write(fates_log(),*)'dif ground absorption error',currentPatch%sabs_dif(ib) , &
                        (currentPatch%tr_soil_dif(ib)* &
                        (1.0_r8-currentPatch%gnd_alb_dif(ib)  )),currentPatch%NCL_p,ib,sum(ftweight(1,1:numpft,1))
                endif
             end if
          endif

          if (radtype == idirect)then
             error = (forc_dir(radtype) + forc_dif(radtype)) - &
                  (fabd_parb_out(ib)  + albd_parb_out(ib) + currentPatch%sabs_dir(ib))
          else
             error = (forc_dir(radtype) + forc_dif(radtype)) - &
                  (fabi_parb_out(ib)  + albi_parb_out(ib) + currentPatch%sabs_dif(ib))
          endif

          ! ignore the current patch radiation error if the veg-covered fraction of the patch is really small
          if ( (currentPatch%total_canopy_area / currentPatch%area) .gt. tolerance ) then
             ! normalize rad error by the veg-covered fraction of the patch because that is
             ! the only part that this code applies to
             currentPatch%rad_error(ib) = currentPatch%rad_error(ib) + error &
                  * currentPatch%total_canopy_area / currentPatch%area

          endif

          lai_reduction(:) = 0.0_r8
          do L = 1, currentPatch%NCL_p
             do ft =1,numpft
                if (currentPatch%canopy_mask(L,ft) == 1)then
                   do iv = 1, currentPatch%nrad(L,ft)
                      if (lai_change(L,ft,iv) > 0.0_r8)then
                         lai_reduction(L) = max(lai_reduction(L),lai_change(L,ft,iv))
                      endif
                   enddo
                endif
             enddo
          enddo

          if (radtype == idirect)then
             !here we are adding a within-ED radiation scheme tolerance, and then adding the diffrence onto the albedo
             !it is important that the lower boundary for this is ~1000 times smaller than the tolerance in surface albedo.
             if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                albd_parb_out(ib) = albd_parb_out(ib) + error
                !this terms adds the error back on to the albedo. While this is partly inexcusable, it is
                ! in the medium term a solution that
                ! prevents the model from crashing with small and occasional energy balances issues.
                ! These are extremely difficult to debug, many have been solved already, leading
                ! to the complexity of this code, but where the system generates occasional errors, we
                ! will deal with them for now.
             end if
             
             if (abs(error)  >  0.15_r8)then
                if(debug)then
                   write(fates_log(),*) 'Large Dir Radn consvn error',error ,ib
                   write(fates_log(),*) 'diags', albd_parb_out(ib), ftdd_parb_out(ib), &
                        ftid_parb_out(ib), fabd_parb_out(ib)
                   write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   write(fates_log(),*) 'ftweight',ftweight(1,1:numpft,1:diag_nlevleaf)
                   write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                   write(fates_log(),*) 'ground albedo diffuse (ib)', currentPatch%gnd_alb_dir(ib)
                end if
                albd_parb_out(ib) = albd_parb_out(ib) + error
             end if
          else

             if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                albi_parb_out(ib) = albi_parb_out(ib) + error
             end if

             if (abs(error)  >  0.15_r8)then
                if(debug)then
                   write(fates_log(),*)  'lg Dif Radn consvn error',error ,ib
                   write(fates_log(),*) 'diags', albi_parb_out(ib), ftii_parb_out(ib), &
                        fabi_parb_out(ib)
                   !write(fates_log(),*) 'lai_change',lai_change(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   !write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   !write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   !write(fates_log(),*) 'ftweight',ftweight(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                   write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                   write(fates_log(),*) 'ground albedo diffuse (ib)', currentPatch%gnd_alb_dir(ib)
                   !write(fates_log(),*) 'rhol',rhol(1:numpft,:)
                   !write(fates_log(),*) 'ftw',sum(ftweight(1,1:numpft,1)),ftweight(1,1:numpft,1)
                   !write(fates_log(),*) 'present',currentPatch%canopy_mask(1,1:numpft)
                   !write(fates_log(),*) 'CAP',currentPatch%canopy_area_profile(1,1:numpft,1)
                end if
                albi_parb_out(ib) = albi_parb_out(ib) + error
             end if

             if (radtype == idirect)then
                error = (forc_dir(radtype) + forc_dif(radtype)) - &
                     (fabd_parb_out(ib)  + albd_parb_out(ib) + currentPatch%sabs_dir(ib))
             else
                error = (forc_dir(radtype) + forc_dif(radtype)) - &
                     (fabi_parb_out(ib)  + albi_parb_out(ib) + currentPatch%sabs_dif(ib))
             endif

             if(debug) then
                if (abs(error)  >  0.00000001_r8)then
                   write(fates_log(),*)  'there is still error after correction',error ,ib
                end if
             end if
             
          end if
       end do !num_swb

    enddo ! rad-type


  end associate
  return
end subroutine PatchNormanRadiation


end module FatesNormanRadMod
