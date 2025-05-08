module FATESPlantRespPhotosynthMod

  !-------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates the plant respiration and photosynthetic fluxes for the FATES model
  ! This code is similar to and was originally based off of the 'photosynthesis'
  ! subroutine in the CLM model.
  !
  ! Parameter for activation and deactivation energies were taken from:
  ! Activation energy, from:
  ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
  ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282
  ! High temperature deactivation, from:
  ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
  ! The factor "c" scales the deactivation to a value of 1.0 at 25C
  ! Photosynthesis and stomatal conductance parameters, from:
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  ! ------------------------------------------------------------------------------------

  ! !USES:

  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  use FatesGlobals,      only : FatesWarn,N2S,A2S,I2S
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : itrue
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
  use FatesConstantsMod, only : wm2_to_umolm2s
  use FatesConstantsMod, only : nocomp_bareground
  use FatesInterfaceTypesMod, only : hlm_use_planthydro
  use FatesInterfaceTypesMod, only : hlm_parteh_mode
  use FatesInterfaceTypesMod, only : numpft
  use FatesInterfaceTypesMod, only : nleafage
  use FatesUtilsMod,          only : QuadraticRoots => QuadraticRootsSridharachary
  use EDParamsMod,           only : maxpft
  use EDParamsMod,       only : nlevleaf
  use EDParamsMod,       only : nclmax
  use EDTypesMod,        only : do_fates_salinity
  use EDParamsMod,       only : q10_mr
  use FatesPatchMod,     only : fates_patch_type
  use FatesCohortMod,    only : fates_cohort_type
  use FatesConstantsMod, only : lmrmodel_ryan_1991
  use FatesConstantsMod, only : lmrmodel_atkin_etal_2017
  use PRTGenericMod,     only : prt_carbon_allom_hyp
  use PRTGenericMod,     only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,     only : carbon12_element
  use PRTGenericMod,     only : nitrogen_element
  use PRTGenericMod,     only : leaf_organ
  use PRTGenericMod,     only : fnrt_organ
  use PRTGenericMod,     only : sapw_organ
  use PRTGenericMod,     only : store_organ
  use PRTGenericMod,     only : repro_organ
  use PRTGenericMod,     only : struct_organ
  use EDParamsMod,       only : maintresp_nonleaf_baserate
  use PRTParametersMod,  only : prt_params
  use EDPftvarcon      , only : EDPftvarcon_inst
  use FatesRadiationMemMod, only : norman_solver,twostr_solver
  use FatesRadiationMemMod, only : ipar
  use FatesTwoStreamUtilsMod, only : FatesGetCohortAbsRad
  use FatesAllometryMod     , only : VegAreaLayer
  use FatesInterfaceTypesMod, only : hlm_radiation_model
  use FatesInterfaceTypesMod, only : hlm_maintresp_leaf_model
  use LeafBiophysicsMod, only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod, only : LeafHumidityStomaResis
  use LeafBiophysicsMod, only : GetCanopyGasParameters
  use LeafBiophysicsMod, only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod, only : LeafLayerMaintenanceRespiration_Atkin_etal_2017
  use LeafBiophysicsMod, only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod, only : LowstorageMainRespReduction
  use LeafBiophysicsMod, only : rsmax0
  use LeafBiophysicsMod, only : DecayCoeffVcmax
  use LeafBiophysicsMod, only : VeloToMolarCF
  use FatesRadiationMemMod, only : idirect
  
  ! CIME Globals
  use shr_log_mod , only      : errMsg => shr_log_errMsg

  implicit none
  private

  public :: FatesPlantRespPhotosynthDrive ! Called by the HLM-Fates interface

  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  character(len=1024) :: warn_msg   ! for defining a warning message

  logical   ::  debug = .false.
  !-------------------------------------------------------------------------------------

contains

  !--------------------------------------------------------------------------------------

  subroutine FatesPlantRespPhotosynthDrive (nsites, sites,bc_in,bc_out,dtime)

    ! -----------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    ! -----------------------------------------------------------------------------------

    ! !USES:
    use EDTypesMod        , only : ed_site_type
    use EDParamsMod       , only : dinc_vai
    use EDParamsMod       , only : dlower_vai
    use FatesInterfaceTypesMod , only : bc_in_type
    use FatesInterfaceTypesMod , only : bc_out_type
    use EDCanopyStructureMod, only : calc_areaindex
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : umol_per_mmol
   
    use FatesParameterDerivedMod, only : param_derived
    use FatesAllometryMod, only : bleaf, bstore_allom
    use FatesAllometryMod, only : storage_fraction_of_target
    use FatesAllometryMod, only : set_root_fraction
    use DamageMainMod, only : GetCrownReduction
    use FatesInterfaceTypesMod, only : hlm_use_tree_damage

    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    real(r8),intent(in)                     :: dtime


    ! LOCAL VARIABLES:
    ! -----------------------------------------------------------------------------------
    type (fates_patch_type) , pointer :: currentPatch
    type (fates_cohort_type), pointer :: currentCohort

    ! -----------------------------------------------------------------------------------
    ! These three arrays hold leaf-level biophysical rates that are calculated
    ! in one loop and then sent to the cohorts in another loop.  If hydraulics are
    ! on, we calculate a unique solution for each level-cohort-layer combination.
    ! If we are not using hydraulics, we calculate a unique solution for each
    ! level-pft-layer combination.  Thus the following three arrays are statically
    ! allocated for the maximum space of the two cases (numCohortsPerPatch)
    ! The "_z" suffix indicates these variables are discretized at the "leaf_layer"
    ! scale.
    ! Note: For these temporary arrays, we have the leaf layer dimension first
    ! and the canopy layer last. This order is chosen for efficiency. The arrays
    ! such as leaf area that are bound to the patch structure DO NOT follow this order
    ! as they are used in many other parts of the code with different looping, we
    ! are not modifying its order now.
    ! -----------------------------------------------------------------------------------

    ! leaf maintenance (dark) respiration [umol CO2/m**2/s]
    real(r8) :: lmr_z(nlevleaf,maxpft,nclmax)

    ! stomatal resistance [s/m]
    real(r8) :: rs_z(nlevleaf,maxpft,nclmax)

    ! net leaf photosynthesis averaged over sun and shade leaves. [umol CO2/m**2/s]
    real(r8) :: anet_av_z(nlevleaf,maxpft,nclmax)

    ! photsynthesis
    real(r8) :: psn_z(nlevleaf,maxpft,nclmax) 

    ! carbon 13 in newly assimilated carbon at leaf level
    real(r8) :: c13disc_z(nlevleaf,maxpft,nclmax) 
    
    ! Mask used to determine which leaf-layer biophysical rates have been
    ! used already
    logical :: rate_mask_z(nlevleaf,maxpft,nclmax)

    real(r8) :: vcmax_z                          ! leaf layer maximum rate of carboxylation
                                                 ! (umol co2/m**2/s)
    real(r8) :: jmax_z                           ! leaf layer maximum electron transport rate
                                                 ! (umol electrons/m**2/s)
    real(r8) :: kp_z                             ! leaf layer initial slope of CO2 response
                                                 ! curve (C4 plants)
    real(r8) :: mm_kco2                          ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2                           ! Michaelis-Menten constant for O2 (Pa)
    real(r8) :: co2_cpoint                       ! CO2 compensation point (Pa)
    real(r8) :: btran_eff                        ! effective transpiration wetness factor (0 to 1)
    real(r8) :: kn                               ! leaf nitrogen decay coefficient
    real(r8) :: gb_mol                           ! leaf boundary layer conductance (molar form: [umol /m**2/s])
    real(r8) :: nscaler                          ! leaf nitrogen scaling coefficient
    real(r8) :: rdark_scaler                     ! scaling coefficient for Atkin dark respiration
    real(r8) :: leaf_frac                        ! ratio of to leaf biomass to total alive biomass
    real(r8) :: tcsoi                            ! Temperature response function for root respiration.
    real(r8) :: tcwood                           ! Temperature response function for wood
    real(r8) :: patch_la                         ! exposed leaf area (patch scale)
    real(r8) :: live_stem_n                      ! Live stem (above-ground sapwood)
                                                 ! nitrogen content (kgN/plant)
    real(r8) :: live_croot_n                     ! Live coarse root (below-ground sapwood)
                                                 ! nitrogen content (kgN/plant)
    real(r8) :: sapw_c                           ! Sapwood carbon (kgC/plant)
    real(r8) :: store_c_target                   ! Target storage carbon (kgC/plant)
    real(r8) :: fnrt_c                           ! Fine root carbon (kgC/plant)
    real(r8) :: fnrt_n                           ! Fine root nitrogen content (kgN/plant)
    real(r8) :: leaf_c                           ! Leaf carbon (kgC/plant)
    real(r8) :: leaf_n                           ! leaf nitrogen content (kgN/plant)
    real(r8) :: g_sb_leaves                      ! Mean combined (stomata+boundary layer) leaf conductance [m/s]
                                                 ! over all of the patch's leaves.  The "sb" refers to the combined
                                                 ! "s"tomatal and "b"oundary layer.
                                                 ! This quantity is relevant on leaf surfaces. It does not
                                                 ! have units of /m2 leaf per say, but is implicitly on leaf surfaces
    real(r8) :: r_sb_leaves                      ! Mean leaf resistance over all the patch's leaves [s/m]
                                                 ! This is the direct reciprocal of g_sb_leaves
    real(r8) :: r_stomata                        ! Mean stomatal resistance across all leaves in the patch [s/m]
    real(r8) :: maintresp_reduction_factor       ! factor by which to reduce maintenance
                                                 ! respiration when storage pools are low
    real(r8) :: b_leaf                           ! leaf biomass kgC
    real(r8) :: frac                             ! storage pool as a fraction of target leaf biomass
                                                 ! over each cohort x layer.
    real(r8) :: cohort_eleaf_area                ! This is the effective leaf area [m2] reported by each cohort
    real(r8) :: lnc_top                          ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8) :: lmr25top                         ! canopy top leaf maint resp rate at 25C 
                                                 ! for this plant or pft (umol CO2/m**2/s)
    real(r8) :: lai_canopy_above                 ! the LAI in the canopy layers above the layer of interest
    real(r8) :: leaf_veg_frac                    ! fraction of vegetation area (leaf+stem) that is just leaf
    real(r8) :: cumulative_lai                   ! the cumulative LAI, top down, to the leaf layer of interest
    real(r8) :: leaf_psi                         ! leaf xylem matric potential [MPa] (only meaningful/used w/ hydro)
    real(r8) :: fnrt_mr_layer                    ! fine root maintenance respiation per layer [kgC/plant/s]
    integer  :: isunsha                          ! Index for differentiating sun and shade
    real(r8) :: fnrt_mr_nfix_layer               ! fineroot maintenance respiration
                                                 ! specifically for symbiotic fixation [kgC/plant/layer/s]
    real(r8) :: nfix_layer                       ! Nitrogen fixed in each layer this timestep [kgN/plant/layer/timestep]
    real(r8), allocatable :: rootfr_ft(:,:)      ! Root fractions per depth and PFT
    real(r8) :: agb_frac                         ! fraction of biomass aboveground
    real(r8) :: branch_frac                      ! fraction of aboveground woody biomass in branches
    real(r8) :: crown_reduction                  ! reduction in crown biomass from damage
    real(r8) :: sapw_c_bgw                       ! belowground sapwood
    real(r8) :: sapw_c_agw                       ! aboveground sapwood
    real(r8) :: sapw_c_undamaged                 ! the target sapwood of an undamaged tree
    real(r8) :: sapw_n                           ! sapwood nitrogen
    real(r8) :: sapw_n_bgw                       ! nitrogen in belowground portion of sapwood
    real(r8) :: sapw_n_agw                       ! nitrogen in aboveground portion of sapwood
    real(r8) :: sapw_n_undamaged                 ! nitrogen in sapwood of undamaged tree
    real(r8) :: rd_abs_leaf, rb_abs_leaf         ! Watts of diffuse and beam light absorbed by leaves over this
                                                 ! depth interval and ground footprint (m2)
    real(r8) :: r_abs_stem, r_abs_snow           ! Watts of light absorbed by stems and snow over this depth interval and
                                                 ! ground footprint 
    real(r8) :: rb_abs, rd_abs                   ! Total beam and diffuse radiation absorbed over this depth interval
                                                 ! and ground footprint
    real(r8) :: fsun                             ! sun-shade fraction
    real(r8) :: par_per_sunla, par_per_shala     ! PAR per sunlit and shaded leaf area [W/m2 leaf]
    real(r8) :: ac_utest, aj_utest               ! Gross rubisco and rubp limited assimilation (for unit tests)
    real(r8) :: ap_utest, co2_inter_c_utest      ! PEP limited assimilation, and intracellular co2 (for unit tests)
    real(r8),dimension(150) :: cohort_vaitop     ! The top-down integrated vegetation area index
                                                 ! (leaf+stem) at the top of the layer
    real(r8),dimension(150) :: cohort_vaibot     ! The top-down integrated vegetation area index
                                                 ! (leaf+stem) at the bottom of the layer
    real(r8),dimension(150) :: cohort_layer_elai ! exposed leaf area index of the layer
    real(r8),dimension(150) :: cohort_layer_esai ! exposed stem area index of the layer
    real(r8)               :: cohort_elai        ! exposed leaf area index of the cohort 
    real(r8)               :: cohort_esai        ! exposed stem area index of the cohort
    real(r8)               :: laisun,laisha      ! m2 of exposed or shaded leaves per m2 of crown
    real(r8)               :: leaf_area          ! m2 leaf per m2 footprint, for either sunlit or shaded leaf
    real(r8)               :: par_abs            ! absorbed PAR [umol photons/m2leaf/s]
    real(r8)               :: canopy_area        ! For Norman radiation, the fraction of crown area per
                                                 ! canopy area in each layer
                                                 ! they will sum to 1.0 in the fully closed canopy layers
    real(r8)               :: elai               ! effective leaf area index
    real(r8)               :: area_frac          ! this is either f_sun , or 1-f_sun, its for area weighting
    real(r8)               :: psn_ll             ! Leaf level photosyntheis
    real(r8)               :: gstoma_ll          ! leaf level stomatal conductance (separate for sun shade) [m/s]
    real(r8)               :: gstoma             ! stomatal conductance at leaf bin (sun/shade combined) [m/s]
    real(r8)               :: anet_ll            ! leaf level net assimilation   [umol CO2/m**2/s]
    real(r8)               :: c13disc_ll         ! leaf level c13 assimilation
    real(r8)               :: hydr_k_lwp         ! inner leaf humidity scaling coefficient [-]
    real(r8)               :: gs0                ! stomatal intercept, possibly scaled by btran depending on hypothesis
    real(r8)               :: gs1                ! stomatal slope, possibly scaled by btran depending on hypothesis
    real(r8)               :: gs2                ! optional btran scaling factor for Medlyn conductance only, instead
                                                 ! of applying to gs1, this would scale the whole non-intercept
                                                 ! portion of the conductance equation

    real(r8)               :: vmol_cf            ! velocity to molar conductance conversion (m/s) -> (umol/m2/s)
    
    ! -----------------------------------------------------------------------------------
    ! Keeping these two definitions in case they need to be added later
    !
    ! -----------------------------------------------------------------------------------
    !real(r8) :: psncanopy_pa  ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    !real(r8) :: lmrcanopy_pa  ! patch sunlit leaf maintenance respiration rate (umol CO2/m**2/s)

    integer  :: cl,s,iv,j,ps,ft,ifp ! indices
    integer  :: nv                  ! number of leaf layers
    integer  :: NCL_p               ! number of canopy layers in patch
    integer  :: iage                ! loop counter for leaf age classes
    integer  :: solve_iter          ! number of iterations required for photosynthesis solve
    
    ! Parameters
    ! Absolute convergence tolerance on solving intracellular CO2 concentration [Pa]

    real(r8), parameter :: ci_tol = 0.5_r8
    
    ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
    ! (gC/gN/s)
    ! ------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------
    ! Photosynthesis and stomatal conductance parameters, from:
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    ! -----------------------------------------------------------------------------------



    associate(  &
         slatop    => prt_params%slatop , &  ! specific leaf area at top of canopy,
                                             ! projected area basis [m^2/gC]
         woody     => prt_params%woody)  ! Is vegetation woody or not?

      do s = 1,nsites

         ! Multi-layer parameters scaled by leaf nitrogen profile.
         ! Loop through each canopy layer to calculate nitrogen profile using
         ! cumulative lai at the midpoint of the layer



         ! Pre-process some variables that are PFT dependent
         ! but not environmentally dependent
         ! ------------------------------------------------------------------------

         allocate(rootfr_ft(numpft, bc_in(s)%nlevsoil))

         do ft = 1,numpft
            call set_root_fraction(rootfr_ft(ft,:), ft, &
                 bc_in(s)%zi_sisl, &
                 bc_in(s)%max_rooting_depth_index_col)
         end do


         currentpatch => sites(s)%oldest_patch
         do while (associated(currentpatch))

            ifp = currentPatch%patchno
            
            if_notbare: if(currentpatch%nocomp_pft_label.ne.nocomp_bareground)then

               NCL_p = currentPatch%NCL_p

               ! Part I. Zero output boundary conditions
               ! ---------------------------------------------------------------------------
               bc_out(s)%rssun_pa(ifp)     = 0._r8
               bc_out(s)%rssha_pa(ifp)     = 0._r8

               g_sb_leaves = 0._r8
               patch_la    = 0._r8

               ! Part II. Filter out patches
               ! Patch level filter flag for photosynthesis calculations
               ! has a short memory, flags:
               ! 1 = patch has not been called
               ! 2 = patch is currently marked for photosynthesis
               ! 3 = patch has been called for photosynthesis already
               ! ---------------------------------------------------------------------------
               if_filter2: if(bc_in(s)%filter_photo_pa(ifp)==2)then


                  ! Part III. Calculate the number of sublayers for each pft and layer.
                  ! And then identify which layer/pft combinations have things in them.
                  ! Output:
                  ! currentPatch%ncan(:,:)
                  ! currentPatch%canopy_mask(:,:)
                  call UpdateCanopyNCanNRadPresent(currentPatch)


                  ! Part IV.  Identify some environmentally derived parameters:
                  !           These quantities are biologically irrelevant
                  !  Michaelis-Menten constant for CO2 (Pa)
                  !  Michaelis-Menten constant for O2 (Pa)
                  !  CO2 compensation point (Pa)
                  !  leaf boundary layer conductance of h20
                  !  constrained vapor pressure

                  call GetCanopyGasParameters(bc_in(s)%forc_pbot,       & ! in
                       bc_in(s)%oair_pa(ifp),    & ! in
                       bc_in(s)%t_veg_pa(ifp),   & ! in
                       mm_kco2,                  & ! out
                       mm_ko2,                   & ! out
                       co2_cpoint)

                  ! The host models use velocity based conductances and resistance
                  ! this is the factor that converts a conductance from
                  ! [m/s] to [umol/m2/s]
                  vmol_cf = VeloToMolarCF(bc_in(s)%forc_pbot,bc_in(s)%tgcm_pa(ifp))
                  
                  ! ------------------------------------------------------------------------
                  ! Part VI: Loop over all leaf layers.
                  ! The concept of leaf layers is a result of the radiative transfer scheme.
                  ! A leaf layer has uniform radiation environment.  Leaf layers are a group
                  ! of vegetation surfaces (stems and leaves) which inhabit the same
                  ! canopy-layer "CL", have the same functional type "ft" and within those
                  ! two partitions are further partitioned into vertical layers where
                  ! downwelling radiation attenuates in order.
                  ! In this phase we loop over the leaf layers and calculate the
                  ! photosynthesis and respiration of the layer (since all biophysical
                  ! properties are homogeneous).  After this step, we can loop through
                  ! our cohort list, associate each cohort with its list of leaf-layers
                  ! and transfer these quantities to the cohort.
                  ! With plant hydraulics, we must realize that photosynthesis and
                  ! respiration will be different for leaves of each cohort in the leaf
                  ! layers, as they will have there own hydraulic limitations.
                  ! NOTE: Only need to flush mask on the number of used pfts, not the whole
                  ! scratch space.
                  ! ------------------------------------------------------------------------
                  rate_mask_z(:,1:numpft,:) = .false.
                  psn_z(:,:,:)     = 0._r8
                  anet_av_z(:,:,:) = 0._r8
                  c13disc_z(:,:,:) = 0._r8
                  rs_z(:,:,:)      = 0._r8
                  lmr_z(:,:,:)     = 0._r8
                                    
                  if_any_cohorts: if(currentPatch%num_cohorts > 0)then

                     currentCohort => currentPatch%tallest
                     do_cohort_drive: do while (associated(currentCohort)) ! Cohort loop

                        ! Identify the canopy layer (cl), functional type (ft)
                        ! and the leaf layer (IV) for this cohort
                        ft = currentCohort%pft
                        cl = currentCohort%canopy_layer
                        
                        ! Calculate the cohort specific elai profile
                        ! And the top and bottom edges of the veg area index
                        ! of each layer bin are. Note, if the layers
                        ! sink below the ground snow line, then the effective
                        ! LAI and SAI start to shrink to zero, as well as
                        ! the difference between vaitop and vaibot.
                        if(currentCohort%treesai>0._r8)then
                           do iv = 1,currentCohort%nv
                              call VegAreaLayer(currentCohort%treelai, &
                                   currentCohort%treesai,              &
                                   currentCohort%height,               &
                                   iv,                                 &
                                   currentCohort%nv,                   &
                                   currentCohort%pft,                  &
                                   sites(s)%snow_depth,                &
                                   cohort_vaitop(iv),                  &
                                   cohort_vaibot(iv),                  & 
                                   cohort_layer_elai(iv),              &
                                   cohort_layer_esai(iv))
                           end do

                           cohort_elai = sum(cohort_layer_elai(1:currentCohort%nv))
                           cohort_esai = sum(cohort_layer_esai(1:currentCohort%nv))


                        else
                           cohort_layer_elai(:)   = 0._r8
                           cohort_layer_esai(:)   = 0._r8
                           cohort_vaitop(:) = 0._r8
                           cohort_vaibot(:) = 0._r8
                           cohort_elai = 0._r8
                           cohort_esai = 0._r8
                        end if

                        ! MLO. Assuming target to be related to leaf biomass when leaves are fully
                        ! flushed. But unsure whether this call is correct or not, shouldn't we get
                        ! the target value directly from the bstore_allom?
                        call bleaf(currentCohort%dbh,currentCohort%pft,&
                             currentCohort%crowndamage,currentCohort%canopy_trim,1.0_r8,store_c_target)
                        !                     call bstore_allom(currentCohort%dbh,currentCohort%pft, &
                        !                                       currentCohort%canopy_trim,store_c_target)

                        call storage_fraction_of_target(store_c_target, &
                             currentCohort%prt%GetState(store_organ, carbon12_element), &
                             frac)
                        call LowstorageMainRespReduction(frac,currentCohort%pft, &
                             maintresp_reduction_factor)

                        ! are there any leaves of this pft in this layer?
                        canopy_mask_if: if(currentPatch%canopy_mask(cl,ft) == 1)then

                           ! Loop over leaf-layers

                           leaf_veg_frac = currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
                           
                           leaf_layer_loop : do iv = 1,currentCohort%nv

                              ! ------------------------------------------------------------
                              ! If we are doing plant hydro-dynamics (or any run-type
                              ! where cohorts may generate different photosynthetic rates
                              ! of other cohorts in the same canopy-pft-layer combo),
                              ! we re-calculate the leaf biophysical rates for the
                              ! cohort-layer combo of interest.
                              ! but in the vanilla case, we only re-calculate if it has
                              ! not been done yet.
                              ! Other cases where we need to solve for every cohort
                              ! in every leaf layer:  nutrient dynamic mode, multiple leaf
                              ! age classes
                              ! ------------------------------------------------------------

                              rate_mask_if: if ( .not.rate_mask_z(iv,ft,cl) .or. &
                                   (hlm_use_planthydro.eq.itrue) .or. &
                                   (hlm_radiation_model .eq. twostr_solver ) .or. &
                                   (nleafage > 1) .or. &
                                   (hlm_parteh_mode .ne. prt_carbon_allom_hyp )   ) then

                                 
                                 ! These values are incremented, therefore since
                                 ! sometimes we re-do layers, we need to re-zero them as well
                                 ! since it is not an over-write
                                 psn_z(iv,ft,cl) = 0._r8
                                 anet_av_z(iv,ft,cl) = 0._r8
                                 c13disc_z(iv,ft,cl) = 0._r8
                                 
                                 if (hlm_use_planthydro.eq.itrue ) then

                                    btran_eff = currentCohort%co_hydr%btran 
                                    
                                    ! Find the cumulative LAI from the top of this cohort's crown to the
                                    ! center of the current veg layer. If this is the cohort's last layer
                                    ! then the mid-point is between the dlower and the total lai

                                    lai_canopy_above  = sum(currentPatch%canopy_layer_tlai(1:cl-1))
                                    
                                    if(iv == currentCohort%nv) then
                                       cumulative_lai = lai_canopy_above + leaf_veg_frac * (dlower_vai(iv)+0.5_r8*(currentCohort%treelai+currentCohort%treesai-dlower_vai(iv)))
                                    else
                                       cumulative_lai = lai_canopy_above + leaf_veg_frac * (dlower_vai(iv)+0.5_r8*dinc_vai(iv))
                                    end if

                                    leaf_psi = currentCohort%co_hydr%psi_ag(1)

                                 else

                                    btran_eff = currentPatch%btran_ft(ft)
                                    ! For consistency sake, we use total LAI here, and not exposed
                                    ! if the plant is under-snow, it will be effectively dormant for 
                                    ! the purposes of nscaler

                                    cumulative_lai = sum(currentPatch%canopy_layer_tlai(1:cl-1))  + &
                                         sum(currentPatch%tlai_profile(cl,ft,1:iv-1)) + &
                                         0.5*currentPatch%tlai_profile(cl,ft,iv)

                                    leaf_psi = fates_unset_r8

                                 end if

                                 if(do_fates_salinity)then
                                    btran_eff = btran_eff*currentPatch%bstress_sal_ft(ft)
                                 endif

                                 ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
                                 ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al 
                                 ! (2010) Biogeosciences, 7, 1833-1859

                                 kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
                                                        prt_params%leafn_vert_scaler_coeff1(ft), &
                                                        prt_params%leafn_vert_scaler_coeff2(ft))

                                 ! Scale for leaf nitrogen profile
                                 nscaler = exp(-kn * cumulative_lai)

                                 ! Leaf maintenance respiration to match the base rate used in CN
                                 ! but with the new temperature functions for C3 and C4 plants.

                                 ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
                                 ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s

                                 ! Then scale this value at the top of the canopy for canopy depth
                                 ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
                                 select case(hlm_parteh_mode)
                                 case (prt_carbon_allom_hyp)

                                    lnc_top  = prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ))/slatop(ft)

                                 case (prt_cnp_flex_allom_hyp)

                                    leaf_c  = currentCohort%prt%GetState(leaf_organ, carbon12_element)
                                    if( (leaf_c*slatop(ft)) > nearzero) then
                                       leaf_n  = currentCohort%prt%GetState(leaf_organ, nitrogen_element)
                                       lnc_top = leaf_n / (slatop(ft) * leaf_c )
                                    else
                                       lnc_top  = prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ))/slatop(ft)
                                    end if

                                    ! If one wants to break coupling with dynamic N conentrations,
                                    ! use the stoichiometry parameter
                                    ! lnc_top  = prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ))/slatop(ft)

                                 end select

                                 
                                 ! Part VII: Calculate dark respiration (leaf maintenance) for this layer

                                 select case (hlm_maintresp_leaf_model)

                                 case (lmrmodel_ryan_1991)

                                    call LeafLayerMaintenanceRespiration_Ryan_1991( lnc_top,     &  ! in
                                         nscaler,                  &  ! in
                                         ft,                       &  ! in
                                         bc_in(s)%t_veg_pa(ifp),   &  ! in
                                         lmr_z(iv,ft,cl))             ! out

                                 case (lmrmodel_atkin_etal_2017)

                                    ! This uses the relationship between leaf N and respiration from Atkin et al 
                                    ! for the top of the canopy, but then scales through the canopy based on a rdark_scaler.
                                    ! To assume proportionality with N through the canopy following Lloyd et al. 2010, use the
                                    ! default parameter value of 2.43, which results in the scaling of photosynthesis and respiration
                                    ! being proportional through the canopy. To have a steeper decrease in respiration than photosynthesis
                                    ! this number can be smaller. There is some observational evidence for this being the case
                                    ! in Lamour et al. 2023. 
                                    
                                    kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
                                         EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(ft), &
                                         EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(ft))

                                    rdark_scaler = exp(-kn * cumulative_lai)
                                    
                                    call LeafLayerMaintenanceRespiration_Atkin_etal_2017( lnc_top, &  ! in
                                         rdark_scaler,                       &  ! in
                                         ft,                                 &  ! in
                                         bc_in(s)%t_veg_pa(ifp),             &  ! in
                                         currentPatch%tveg_lpa%GetMean(),    &  ! in
                                         lmr_z(iv,ft,cl))                       ! out

                                 case default

                                    write (fates_log(),*)'error, incorrect leaf respiration model specified'
                                    call endrun(msg=errMsg(sourcefile, __LINE__))

                                 end select
                                 
                                 ! Pre-process PAR absorbed per unit leaf area for different schemes
                                 ! par_per_sunla = [W absorbed beam+diffuse radiation / m2 of sunlit leaves]
                                 ! par_per_shala = [W absorbed diffuse radiation / m2 of shaded leaves]
                                 ! fsun          = [m2 of sunlit leaves / m2 of total leaves]
                                 ! laisun:      m2 of exposed leaf, per m2 of crown. If this is the lowest layer
                                 !              for the pft/canopy group, than the m2 per crown is probably not
                                 !              as large as the layer above.
                                 ! ------------------------------------------------------------------

                                 if_radsolver: if(hlm_radiation_model.eq.norman_solver) then

                                    laisun = currentPatch%ed_laisun_z(cl,ft,iv)
                                    laisha = currentPatch%ed_laisha_z(cl,ft,iv)
                                    par_per_sunla = currentPatch%ed_parsun_z(cl,ft,iv)
                                    par_per_shala = currentPatch%ed_parsha_z(cl,ft,iv)
                                    canopy_area   = currentPatch%canopy_area_profile(cl,ft,iv)
                                    fsun = currentPatch%f_sun(cl,ft,iv)
                                    
                                 else    ! Two-stream

                                    if(cohort_layer_elai(iv) > nearzero .and. sites(s)%coszen>0._r8 ) then

                                       call FatesGetCohortAbsRad(currentPatch, currentCohort, ipar, &
                                            cohort_vaitop(iv), cohort_vaibot(iv), cohort_elai, cohort_esai, &
                                            rb_abs, rd_abs, rb_abs_leaf, rd_abs_leaf, fsun)

                                       ! rd_abs_leaf: Watts of diffuse light absorbed by leaves over this
                                       !              depth interval and ground footprint (m2)
                                       ! rd_abs_leaf*fsun  Watts of diffuse light absorbed by sunlit leaves
                                       !                   over this depth interval and ground footprint (m2)
                                       ! rb_abs_leaf       Watts of beam absorbed by sunlit leaves over this
                                       !                   depth interval and ground footprint (m2)
                                       ! cohort_layer_elai*fsun       Leaf area in sunlight within this interval and ground footprint
                                       ! cohort_layer_elai*(1-fsun)   Leaf area in shade within this interval and ground footprint

                                       laisun = (fsun*cohort_layer_elai(iv))
                                       laisha = ((1._r8 - fsun)*cohort_layer_elai(iv))
                                       if(fsun>nearzero) then
                                          par_per_sunla = (rd_abs_leaf*fsun + rb_abs_leaf)! / laisun
                                       else
                                          par_per_sunla = 0._r8
                                       end if
                                       par_per_shala = rd_abs_leaf*(1._r8-fsun) !/ laisha
                                       canopy_area = 1._r8 !currentPatch%canopy_area_profile(cl,ft,iv)
                                       
                                    else

                                       par_per_sunla = 0._r8
                                       par_per_shala = 0._r8
                                       laisun = 0.5_r8*cohort_layer_elai(iv)
                                       laisha = 0.5_r8*cohort_layer_elai(iv)
                                       canopy_area = 1._r8 !currentPatch%canopy_area_profile(cl,ft,iv)
                                       fsun = 0.5_r8 !avoid div0, should have no impact
                                       
                                    end if

                                 end if if_radsolver

                                 ! Perform photosynthesis calculations on sunlit and shaded leaves
                                 ! ---------------------------------------------------------------

                                 ! Calculate leaf boundary layer conductance in molar form [umol/m2/s]
                                 gb_mol = (1._r8/bc_in(s)%rb_pa(ifp)) * vmol_cf
                                 
                                 gstoma = 0._r8
                                 do_sunsha: do isunsha = 1,2

                                    ! Determine absorbed PAR per square meter of leaf
                                    ! If there is no leaf area perform a trivial solution

                                    if(isunsha == idirect) then
                                       leaf_area = laisun*canopy_area
                                       par_abs   = ConvertPar(leaf_area, par_per_sunla)
                                       area_frac = fsun
                                    else
                                       leaf_area = laisha*canopy_area
                                       par_abs   = ConvertPar(leaf_area, par_per_shala)
                                       area_frac = 1._r8 - fsun
                                    end if
                                    
                                    if( leaf_area < nearzero ) then

                                       ! Note: With no leaf area do not increment
                                       ! any fluxes. Assume a nominal conductance
                                       ! of maximum resistance
                                       gstoma = gstoma + area_frac/rsmax0
                                       
                                       cycle do_sunsha
                                    end if
                                    
                                    ! Part VII: Calculate (1) maximum rate of carboxylation (vcmax),
                                    ! (2) maximum electron transport rate, (3) triose phosphate
                                    ! utilization rate and (4) the initial slope of CO2 response curve
                                    ! (C4 plants). Earlier we calculated their base rates as dictated
                                    ! by their plant functional type and some simple scaling rules for
                                    ! nitrogen limitation baesd on canopy position (not prognostic).
                                    ! These rates are the specific rates used in the actual photosynthesis
                                    ! calculations that take localized environmental effects (temperature)
                                    ! into consideration.
                                    
                                    call LeafLayerBiophysicalRates(ft,        &  ! in
                                         currentCohort%vcmax25top,            &  ! in
                                         currentCohort%jmax25top,             &  ! in
                                         currentCohort%kp25top,               &  ! in
                                         nscaler,                             &  ! in
                                         bc_in(s)%t_veg_pa(ifp),              &  ! in
                                         bc_in(s)%dayl_factor_pa(ifp),        &  ! in
                                         currentPatch%tveg_lpa%GetMean(),     &  ! in
                                         currentPatch%tveg_longterm%GetMean(),&  ! in
                                         btran_eff,                           &  ! in
                                         vcmax_z,                             &  ! out
                                         jmax_z,                              &  ! out
                                         kp_z,                                &  ! out
                                         gs0,                                 &  ! out
                                         gs1,                                 &  ! out
                                         gs2 )                                   ! out


                                    if ( (hlm_use_planthydro.eq.itrue .and. EDPftvarcon_inst%hydr_k_lwp(ft)>nearzero) ) then
                                       hydr_k_lwp = EDPftvarcon_inst%hydr_k_lwp(ft)
                                    else
                                       hydr_k_lwp = 1._r8
                                    end if

                                    call LeafLayerPhotosynthesis(            & !
                                         par_abs,                            &  ! in
                                         ft,                                 &  ! in
                                         vcmax_z,                            &  ! in
                                         jmax_z,                             &  ! in
                                         kp_z,                               &  ! in
                                         gs0,                                &  ! in
                                         gs1,                                &  ! in
                                         gs2,                                &  ! in
                                         bc_in(s)%t_veg_pa(ifp),             &  ! in
                                         bc_in(s)%forc_pbot,                 &  ! in
                                         bc_in(s)%cair_pa(ifp),              &  ! in
                                         bc_in(s)%oair_pa(ifp),              &  ! in
                                         bc_in(s)%esat_tv_pa(ifp),           &  ! in
                                         gb_mol,                             &  ! in
                                         bc_in(s)%eair_pa(ifp),              &  ! in
                                         mm_kco2,                            &  ! in
                                         mm_ko2,                             &  ! in
                                         co2_cpoint,                         &  ! in
                                         lmr_z(iv,ft,cl),                    &  ! in
                                         ci_tol,                             &  ! in
                                         psn_ll,                             &  ! out
                                         gstoma_ll,                          &  ! out
                                         anet_ll,                            &  ! out
                                         c13disc_ll,                         &  ! out
                                         co2_inter_c_utest,                  &  ! out (unit tests)
                                         solve_iter)                            ! out performance tracking

                                    ! Average output quantities across sunlit and shaded leaves
                                    ! Convert from molar to velocity (umol /m**2/s) to (m/s)
                                    gstoma = gstoma + area_frac*(gstoma_ll / vmol_cf) 

                                    
                                    psn_z(iv,ft,cl) = psn_z(iv,ft,cl) + area_frac * psn_ll
                                    anet_av_z(iv,ft,cl) = anet_av_z(iv,ft,cl) + area_frac * anet_ll
                                    c13disc_z(iv,ft,cl) = c13disc_z(iv,ft,cl) + area_frac * c13disc_ll
                                    
                                 end do do_sunsha
                                 
                                 ! Stomatal resistance of the leaf-layer
                                 if ( (hlm_use_planthydro.eq.itrue .and. EDPftvarcon_inst%hydr_k_lwp(ft)>nearzero) ) then

                                    rs_z(iv,ft,cl) = LeafHumidityStomaResis(leaf_psi, EDPftvarcon_inst%hydr_k_lwp(ft), &
                                         bc_in(s)%t_veg_pa(ifp),bc_in(s)%cair_pa(ifp),bc_in(s)%forc_pbot, &
                                         bc_in(s)%rb_pa(ifp), gstoma, ft, bc_in(s)%esat_tv_pa(ifp) )
                                    
                                 else
                                    rs_z(iv,ft,cl)= 1._r8/gstoma
                                 end if
                                 
                                 rate_mask_z(iv,ft,cl) = .true.
                                 
                              end if rate_mask_if
                           end do leaf_layer_loop

                           
                           ! Zero cohort flux accumulators.
                           currentCohort%resp_m_tstep = 0.0_r8
                           currentCohort%gpp_tstep  = 0.0_r8
                           currentCohort%rdark      = 0.0_r8
                           currentCohort%ts_net_uptake = 0.0_r8
                           currentCohort%c13disc_clm = 0.0_r8

                           ! ---------------------------------------------------------------
                           ! Part VII: Transfer leaf flux rates (like maintenance respiration,
                           ! carbon assimilation and conductance) that are defined by the
                           ! leaf layer (which is area independent, ie /m2) onto each cohort
                           ! (where the rates become per cohort, ie /individual). Most likely
                           ! a sum over layers.
                           ! ---------------------------------------------------------------
                           nv = currentCohort%nv

                           if(hlm_radiation_model.eq.norman_solver) then

                              call ScaleLeafLayerFluxToCohort(nv,         & !in
                                   psn_z(1:nv,ft,cl),                     & !in
                                   lmr_z(1:nv,ft,cl),                     & !in
                                   rs_z(1:nv,ft,cl),                      & !in
                                   currentPatch%elai_profile(cl,ft,1:nv), & !in
                                   c13disc_z(1:nv,ft,cl),                 & !in
                                   currentCohort%c_area,                  & !in
                                   currentCohort%n,                       & !in
                                   bc_in(s)%rb_pa(ifp),                   & !in
                                   maintresp_reduction_factor,            & !in
                                   currentCohort%g_sb_laweight,           & !out
                                   currentCohort%gpp_tstep,               & !out
                                   currentCohort%rdark,                   & !out
                                   currentCohort%c13disc_clm,             & !out
                                   cohort_eleaf_area)                       !out

                           else

                              
                              
                              call ScaleLeafLayerFluxToCohort(nv,         & !in
                                   psn_z(1:nv,ft,cl),                     & !in
                                   lmr_z(1:nv,ft,cl),                     & !in
                                   rs_z(1:nv,ft,cl),                      & !in
                                   cohort_layer_elai(1:nv),               & !in
                                   c13disc_z(1:nv,ft,cl),                 & !in
                                   currentCohort%c_area,                  & !in
                                   currentCohort%n,                       & !in
                                   bc_in(s)%rb_pa(ifp),                   & !in
                                   maintresp_reduction_factor,            & !in
                                   currentCohort%g_sb_laweight,           & !out
                                   currentCohort%gpp_tstep,               & !out
                                   currentCohort%rdark,                   & !out
                                   currentCohort%c13disc_clm,             & !out
                                   cohort_eleaf_area)                       !out
                           end if

                              
                           ! Net Uptake does not need to be scaled, just transfer directly
                           currentCohort%ts_net_uptake(1:nv) = anet_av_z(1:nv,ft,cl) * umolC_to_kgC
                           
                        else
                           
                           ! In this case, the cohort had no leaves,
                           ! so no productivity,conductance, transpiration uptake
                           ! or dark respiration
                           cohort_eleaf_area       = 0.0_r8
                           currentCohort%gpp_tstep = 0.0_r8
                           currentCohort%rdark = 0.0_r8
                           currentCohort%g_sb_laweight = 0.0_r8
                           currentCohort%ts_net_uptake(:) = 0.0_r8

                        end if canopy_mask_if


                        ! ------------------------------------------------------------------
                        ! Part VIII: Calculate maintenance respiration in the sapwood and
                        ! fine root pools.
                        ! ------------------------------------------------------------------

                        ! Calculate the amount of nitrogen in the above and below ground
                        ! stem and root pools, used for maint resp
                        ! We are using the fine-root C:N ratio as an approximation for
                        ! the sapwood pools.
                        ! Units are in (kgN/plant)
                        ! ------------------------------------------------------------------

                        sapw_c   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
                        fnrt_c   = currentCohort%prt%GetState(fnrt_organ, carbon12_element)

                        if (hlm_use_tree_damage .eq. itrue) then

                           ! Crown damage currenly only reduces the aboveground portion of 
                           ! sapwood. Therefore we calculate the aboveground and the belowground portion 
                           ! sapwood for use in stem respiration. 
                           call GetCrownReduction(currentCohort%crowndamage, crown_reduction)

                        else
                           crown_reduction = 0.0_r8
                        end if

                        ! If crown reduction is zero, undamaged sapwood target will equal sapwood carbon
                        agb_frac = prt_params%allom_agb_frac(currentCohort%pft)
                        branch_frac = param_derived%branch_frac(currentCohort%pft)
                        sapw_c_undamaged = sapw_c / (1.0_r8 - (agb_frac * branch_frac * crown_reduction))

                        ! Undamaged below ground portion
                        sapw_c_bgw = sapw_c_undamaged * (1.0_r8 - agb_frac)

                        ! Damaged aboveground portion
                        sapw_c_agw = sapw_c - sapw_c_bgw                         


                        select case(hlm_parteh_mode)
                        case (prt_carbon_allom_hyp)

                           live_stem_n = sapw_c_agw * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))

                           live_croot_n = sapw_c_bgw * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))

                           fnrt_n = fnrt_c * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ))

                        case(prt_cnp_flex_allom_hyp)

                           live_stem_n = prt_params%allom_agb_frac(currentCohort%pft) * &
                                currentCohort%prt%GetState(sapw_organ, nitrogen_element)

                           live_croot_n = (1.0_r8-prt_params%allom_agb_frac(currentCohort%pft)) * &
                                currentCohort%prt%GetState(sapw_organ, nitrogen_element)


                           fnrt_n = currentCohort%prt%GetState(fnrt_organ, nitrogen_element)

                           if (hlm_use_tree_damage .eq. itrue) then

                              sapw_n = currentCohort%prt%GetState(sapw_organ, nitrogen_element)

                              sapw_n_undamaged = sapw_n / &
                                   (1.0_r8 - (agb_frac * branch_frac * crown_reduction))

                              sapw_n_bgw = sapw_n_undamaged * (1.0_r8 - agb_frac)
                              sapw_n_agw = sapw_n - sapw_n_bgw

                              live_croot_n = sapw_n_bgw

                              live_stem_n = sapw_n_agw

                           end if

                           ! If one wants to break coupling with dynamic N conentrations,
                           ! use the stoichiometry parameter
                           !
                           ! live_stem_n = prt_params%allom_agb_frac(currentCohort%pft) * &
                           !               sapw_c * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))
                           ! live_croot_n = (1.0_r8-prt_params%allom_agb_frac(currentCohort%pft)) * &
                           !               sapw_c * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))
                           ! fnrt_n = fnrt_c * prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ))


                        case default


                        end select

                        !------------------------------------------------------------------------------
                        ! Calculate Whole Plant Respiration
                        ! (this doesn't really need to be in this iteration at all, surely?)
                        ! Response: (RGK 12-2016): I think the positioning of these calls is
                        ! appropriate as of now.  Maintenance calculations in sapwood and roots
                        ! vary by cohort and with changing temperature at the minimum, and there are
                        ! no sub-pools chopping up those pools any finer that need to be dealt with.
                        !------------------------------------------------------------------------------

                        ! Live stem MR (kgC/plant/s) (above ground sapwood)
                        ! ------------------------------------------------------------------
                        if ( int(woody(ft)) == itrue) then
                           tcwood = q10_mr**((bc_in(s)%t_veg_pa(ifp)-tfrz - 20.0_r8)/10.0_r8)
                           ! kgC/s = kgN * kgC/kgN/s
                           currentCohort%livestem_mr  = live_stem_n * maintresp_nonleaf_baserate * tcwood * maintresp_reduction_factor
                        else
                           currentCohort%livestem_mr  = 0._r8
                        end if


                        ! Fine Root MR  (kgC/plant/s)
                        ! and calculate the N fixation rate as a function of the fixation-specific root respiration
                        ! for now use dev_arbitrary_pft as scaling term between 0 and 1 as additional increment of root respiration used for N fixation
                        ! ------------------------------------------------------------------
                        currentCohort%froot_mr = 0._r8
                        currentCohort%sym_nfix_tstep = 0._r8

                        ! n_fixation is integrated over the course of the day
                        ! this variable is zeroed at the end of the FATES dynamics sequence

                        do j = 1,bc_in(s)%nlevsoil
                           tcsoi  = q10_mr**((bc_in(s)%t_soisno_sl(j)-tfrz - 20.0_r8)/10.0_r8)

                           fnrt_mr_layer = fnrt_n * maintresp_nonleaf_baserate * tcsoi * rootfr_ft(ft,j) * maintresp_reduction_factor

                           ! calculate the cost of carbon for N fixation in each soil layer and calculate N fixation rate based on that [kgC / kgN]

                           call RootLayerNFixation(bc_in(s)%t_soisno_sl(j),ft,dtime,fnrt_mr_layer,fnrt_mr_nfix_layer,nfix_layer)

                           currentCohort%froot_mr = currentCohort%froot_mr + fnrt_mr_nfix_layer + fnrt_mr_layer 

                           currentCohort%sym_nfix_tstep = currentCohort%sym_nfix_tstep + nfix_layer


                        enddo

                        ! Coarse Root MR (kgC/plant/s) (below ground sapwood)
                        ! ------------------------------------------------------------------
                        if ( int(woody(ft)) == itrue) then
                           currentCohort%livecroot_mr = 0._r8
                           do j = 1,bc_in(s)%nlevsoil
                              ! Soil temperature used to adjust base rate of MR
                              tcsoi  = q10_mr**((bc_in(s)%t_soisno_sl(j)-tfrz - 20.0_r8)/10.0_r8)
                              currentCohort%livecroot_mr = currentCohort%livecroot_mr + &
                                   live_croot_n * maintresp_nonleaf_baserate * tcsoi * &
                                   rootfr_ft(ft,j) * maintresp_reduction_factor
                           enddo
                        else
                           currentCohort%livecroot_mr = 0._r8
                        end if


                        ! ------------------------------------------------------------------
                        ! Part IX: Perform some unit conversions (rate to integrated) and
                        ! calcualate some fluxes that are sums and nets of the base fluxes
                        ! ------------------------------------------------------------------

                        ! add on whole plant respiration values in kgC/indiv/s-1
                        currentCohort%resp_m_tstep = currentCohort%livestem_mr + &
                             currentCohort%livecroot_mr + &
                             currentCohort%froot_mr + &
                             currentCohort%rdark
                        
                        ! no drought response right now.. something like:
                        ! resp_m_tstep = resp_m_tstep * (1.0_r8 - currentPatch%btran_ft(currentCohort%pft) * &
                        !                    EDPftvarcon_inst%resp_drought_response(ft))

                        ! convert from kgC/indiv/s to kgC/indiv/timestep
                        currentCohort%resp_m_tstep  = currentCohort%resp_m_tstep  * dtime
                        currentCohort%gpp_tstep     = currentCohort%gpp_tstep * dtime
                        currentCohort%ts_net_uptake = currentCohort%ts_net_uptake * dtime
                        
                        ! save as a diagnostic the un-throttled maintenance respiration to be able to know how strong this is
                        currentCohort%resp_m_unreduced = currentCohort%resp_m_tstep / maintresp_reduction_factor
                        
                        ! Accumulate the combined conductance (stomatal+leaf boundary layer)
                        ! Note that currentCohort%g_sb_laweight is weighted by the leaf area
                        ! of each cohort and has units of [m/s] * [m2 leaf]

                        g_sb_leaves  = g_sb_leaves + currentCohort%g_sb_laweight

                        ! Accumulate the total effective leaf area from all cohorts
                        ! in this patch. Normalize by canopy area outside the loop
                        patch_la = patch_la + cohort_eleaf_area

                        currentCohort => currentCohort%shorter
                     enddo do_cohort_drive

                  end if if_any_cohorts

                  ! Normalize canopy total conductance by the effective LAI
                  ! The value here was integrated over each cohort x leaf layer
                  ! and was weighted by m2 of effective leaf area for each layer
                  
                  if_any_lai: if(patch_la>nearzero) then

                     ! Normalize the leaf-area weighted canopy conductance
                     ! The denominator is the total effective leaf area in the canopy,
                     ! units of [m/s]*[m2] / [m2] = [m/s]
                     
                     g_sb_leaves = g_sb_leaves / patch_la
                     
                     if_above_mincond: if( g_sb_leaves > (1._r8/rsmax0) ) then
                        
                        ! Combined mean leaf resistance is
                        ! the inverse of mean leaf conductance
                        r_sb_leaves  = 1.0_r8/g_sb_leaves
                        
                        if (r_sb_leaves<bc_in(s)%rb_pa(ifp)) then
                           write(fates_log(),*) 'Combined canopy resistance was somehow smaller than'
                           write(fates_log(),*) 'its boundary layer resistance component'
                           write(fates_log(),*) 'r_sb_leaves [s/m]: ',r_sb_leaves
                           write(fates_log(),*) 'bc_in(s)%rb_pa(ifp) [s/m]: ',bc_in(s)%rb_pa(ifp)
                           call endrun(msg=errMsg(sourcefile, __LINE__))
                        end if
                        
                        ! Mean leaf stomatal resistance for all patch leaves
                        r_stomata = (r_sb_leaves - bc_in(s)%rb_pa(ifp))
                        
                     else !if_above_mincond
                        
                        ! Here we prevent super high resistances
                        ! and use a nominal value when conductance is low
                        r_stomata = rsmax0
                        
                     end if if_above_mincond
                     
                     ! This will be multiplied by scaled by effective LAI in the host model
                     ! when it comes time to calculate a flux rate per unit ground
                     bc_out(s)%rssun_pa(ifp) = r_stomata
                     bc_out(s)%rssha_pa(ifp) = r_stomata
                     
                     ! This value is used for diagnostics, the molar form of conductance
                     ! is what is used in the field usually, so we track that form
                     ! vmol_cf :  s m**2/umol -> s/m (ideal gas conversion) [umol/m3]

                     currentPatch%c_stomata  = vmol_cf / r_stomata
                     
                  else !if_any_lai
                     
                     ! But this will prevent it from using an unintialized value
                     bc_out(s)%rssun_pa(ifp) = rsmax0
                     bc_out(s)%rssha_pa(ifp) = rsmax0
                     
                     ! This value is used for diagnostics, the molar form of conductance
                     ! is what is used in the field usually, so we track that form
                     currentPatch%c_stomata  = vmol_cf / rsmax0
                     
                  end if if_any_lai

                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_lblayer = vmol_cf / bc_in(s)%rb_pa(ifp)
                  
               end if if_filter2
               
            end if if_notbare

            currentPatch => currentPatch%younger
         end do

         deallocate(rootfr_ft)

      end do !site loop

    end associate
  end subroutine FatesPlantRespPhotosynthDrive

  ! ===========================================================================================


  subroutine RootLayerNFixation(t_soil,ft,dtime,fnrt_mr_layer,fnrt_mr_nfix_layer,nfix_layer)


    ! -------------------------------------------------------------------------------
    ! Symbiotic N Fixation is handled via Houlton et al 2008 and Fisher et al. 2010
    !
    ! A unifying framework for dinitrogen fixation in the terrestrial biosphere
    ! Benjamin Z. Houlton, Ying-Ping Wang, Peter M. Vitousek & Christopher B. Field 
    ! Nature volume 454, pages327–330 (2008)  https://doi.org/10.1038/nature07028
    !
    ! Carbon cost of plant nitrogen acquisition: A mechanistic, globally applicable model
    ! of plant nitrogen uptake, retranslocation, and fixation.  J. B. Fisher,S. Sitch,Y.
    ! Malhi,R. A. Fisher,C. Huntingford,S.-Y. Tan. Global Biogeochemical Cycles. March
    ! 2010 https://doi.org/10.1029/2009GB003621
    !
    ! ------------------------------------------------------------------------------


    real(r8),intent(in) :: t_soil              ! Temperature of the current soil layer [degC]
    integer,intent(in)  :: ft                  ! Functional type index
    real(r8),intent(in) :: dtime               ! Time step length [s]
    real(r8),intent(in) :: fnrt_mr_layer       ! Amount of maintenance respiration in the fine-roots
    ! for all non-fixation related processes [kgC/s]

    real(r8),intent(out) :: fnrt_mr_nfix_layer ! The added maintenance respiration due to nfixation
    ! to be added as a surcharge to non-fixation MR [kgC]
    real(r8),intent(out) :: nfix_layer         ! The amount of N fixed in this layer through
    ! symbiotic activity [kgN]

    real(r8) :: c_cost_nfix                    ! carbon cost of N fixation [kgC/kgN]
    real(r8) :: c_spent_nfix                   ! carbon spent on N fixation, per layer [kgC/plant/timestep]

    ! N fixation parameters from Houlton et al (2008) and Fisher et al (2010)
    real(r8), parameter :: s_fix = -6.25_r8 ! s parameter from FUN model (fisher et al 2010)
    real(r8), parameter :: a_fix = -3.62_r8 ! a parameter from Houlton et al. 2010 (a = -3.62 +/- 0.52)
    real(r8), parameter :: b_fix = 0.27_r8  ! b parameter from Houlton et al. 2010 (b = 0.27 +/-0.04)
    real(r8), parameter :: c_fix = 25.15_r8 ! c parameter from Houlton et al. 2010 (c = 25.15 +/- 0.66)

    ! Amount of C spent (as part of MR respiration) on symbiotic fixation [kgC/s]
    fnrt_mr_nfix_layer  = fnrt_mr_layer * prt_params%nfix_mresp_scfrac(ft)

    ! This is the unit carbon cost for nitrogen fixation. It is temperature dependant [kgC/kgN]
    c_cost_nfix = s_fix * (exp(a_fix + b_fix * (t_soil-tfrz) &
         * (1._r8 - 0.5_r8 * (t_soil-tfrz) / c_fix)) - 2._r8)

    ! Time integrated amount of carbon spent on fixation (in this layer) [kgC/plant/layer/tstep]
    c_spent_nfix = fnrt_mr_nfix_layer  * dtime

    ! Amount of nitrogen fixed in this layer [kgC/plant/layer/tstep]/[kgC/kgN] = [kgN/plant/layer/tstep]
    nfix_layer = c_spent_nfix / c_cost_nfix

    return
  end subroutine RootLayerNFixation


  ! =======================================================================================

  subroutine ScaleLeafLayerFluxToCohort(nv, &                         ! in
                                        psn_llz,     &                ! in 
                                        lmr_llz,     &                ! in 
                                        rs_llz,      &                ! in 
                                        elai_llz,    &                ! in 
                                        c13disc_llz, &                ! in 
                                        c_area,      &                ! in 
                                        nplant,      &                ! in 
                                        rb,          &                ! in 
                                        maintresp_reduction_factor, & ! in
                                        g_sb_laweight, &              ! out 
                                        gpp,         &                ! out 
                                        rdark,       &                ! out 
                                        c13disc_clm, &                ! out 
                                        cohort_eleaf_area )           ! out

    ! ------------------------------------------------------------------------------------
    ! This subroutine effectively integrates leaf carbon fluxes over the
    ! leaf layers to give cohort totals.
    ! Some arguments have the suffix "_llz".  This indicates that the vector
    ! is stratefied in the leaf-layer (ll)  dimension, and is a portion of the calling
    ! array which has the "_z" tag, thus "llz".
    ! ------------------------------------------------------------------------------------

    use FatesConstantsMod, only : umolC_to_kgC

    ! Arguments
    integer, intent(in)  :: nv               ! number of active leaf layers
    real(r8), intent(in) :: psn_llz(nv)      ! layer photosynthesis rate (GPP) [umolC/m2leaf/s]
    real(r8), intent(in) :: lmr_llz(nv)      ! layer dark respiration rate [umolC/m2leaf/s]
    real(r8), intent(in) :: rs_llz(nv)       ! leaf layer stomatal resistance [s/m]
    real(r8), intent(in) :: elai_llz(nv)     ! exposed LAI per layer [m2 leaf/ m2 pft footprint]
    real(r8), intent(in) :: c13disc_llz(nv)  ! leaf layer c13 discrimination, weighted mean
    real(r8), intent(in) :: c_area           ! crown area m2/m2
    real(r8), intent(in) :: nplant           ! indiv/m2
    real(r8), intent(in) :: rb               ! leaf boundary layer resistance (s/m)
    real(r8), intent(in) :: maintresp_reduction_factor  ! factor by which to reduce maintenance respiration
    real(r8), intent(out) :: g_sb_laweight     ! Combined conductance (stomatal + boundary layer) for the cohort
                                               ! weighted by leaf area [m/s]*[m2]
    real(r8), intent(out) :: gpp               ! GPP (kgC/indiv/s)
    real(r8), intent(out) :: rdark             ! Dark Leaf Respiration (kgC/indiv/s)
    real(r8), intent(out) :: cohort_eleaf_area ! Effective leaf area of the cohort [m2]
    real(r8), intent(out) :: c13disc_clm       ! unpacked Cohort level c13 discrimination
    real(r8)              :: sum_weight        ! sum of weight for unpacking d13c flux (c13disc_z) from
                                               ! (canopy_layer, pft, leaf_layer) matrix to cohort (c13disc_clm)

    ! GPP IN THIS SUBROUTINE IS A RATE. THE CALLING ARGUMENT IS GPP_TSTEP. AFTER THIS
    ! CALL THE RATE WILL BE MULTIPLIED BY THE INTERVAL TO GIVE THE INTEGRATED QUANT.

    ! Locals
    integer  :: il                       ! leaf layer index
    real(r8) :: cohort_layer_eleaf_area  ! the effective leaf area of the cohort's current layer [m2]

    cohort_eleaf_area = 0.0_r8
    g_sb_laweight     = 0.0_r8
    gpp               = 0.0_r8
    rdark             = 0.0_r8

    do il = 1, nv        ! Loop over the leaf layers this cohort participates in


       ! Cohort's total effective leaf area in this layer [m2]
       ! leaf area index of the layer [m2/m2 ground] * [m2 ground]
       ! elai_llz is the LAI for the whole PFT. Multiplying this by the ground
       ! area this cohort contributes, give the cohort's portion of the leaf
       ! area in this layer
       cohort_layer_eleaf_area = elai_llz(il) * c_area

       ! Increment the cohort's total effective leaf area [m2]
       cohort_eleaf_area       = cohort_eleaf_area + cohort_layer_eleaf_area

       ! Leaf conductance (stomatal and boundary layer)
       ! This should be the weighted average over the leaf surfaces.
       ! Since this is relevant to the stomata, its weighting should be based
       ! on total leaf area, and not really footprint area
       ! [m/s] * [m2 cohort's leaf layer]
       g_sb_laweight = g_sb_laweight + 1.0_r8/(rs_llz(il)+rb) * cohort_layer_eleaf_area

       ! GPP    [umolC/m2leaf/s] * [m2 leaf ] -> [umolC/s]
       gpp = gpp + psn_llz(il) * cohort_layer_eleaf_area

       ! Dark respiration
       ! [umolC/m2leaf/s] * [m2 leaf] 
       rdark = rdark + lmr_llz(il) * cohort_layer_eleaf_area

    end do



    if (nv > 1) then
       ! cohort%c13disc_clm as weighted mean of d13c flux at all related leave layers
       sum_weight = sum(psn_llz(1:nv-1) * elai_llz(1:nv-1))
       if (sum_weight .eq. 0.0_r8) then
          c13disc_clm = 0.0
       else
          c13disc_clm = sum(c13disc_llz(1:nv-1) * psn_llz(1:nv-1) * elai_llz(1:nv-1)) / sum_weight
       end if

    end if


    ! -----------------------------------------------------------------------------------
    ! We DO NOT normalize g_sb_laweight.
    ! The units that we are passing back are [m/s] * [m2 effective leaf]
    ! We will add these up over the whole patch, and then normalized
    ! by the patch's total leaf area in the calling routine
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Convert dark respiration and GPP from [umol/s] to [kgC/plant/s]
    ! Also, apply the maintenance respiration reduction factor
    ! -----------------------------------------------------------------------------------

    rdark     = rdark * umolC_to_kgC * maintresp_reduction_factor / nplant
    gpp       = gpp * umolC_to_kgC / nplant

    if ( debug ) then
       write(fates_log(),*) 'EDPhoto 816 ', gpp
       write(fates_log(),*) 'EDPhoto 817 ', psn_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 820 ', nv
       write(fates_log(),*) 'EDPhoto 821 ', elai_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 843 ', rdark
       write(fates_log(),*) 'EDPhoto 873 ', nv
       write(fates_log(),*) 'EDPhoto 874 ', cohort_eleaf_area
    endif

    return
  end subroutine ScaleLeafLayerFluxToCohort

  ! =====================================================================================

  subroutine UpdateCanopyNCanNRadPresent(currentPatch)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates two patch level quanities:
    ! currentPatch%ncan   and
    ! currentPatch%canopy_mask
    !
    ! currentPatch%ncan(:,:) is a two dimensional array that indicates
    ! the total number of leaf layers (including those that are not exposed to light)
    ! in each canopy layer and for each functional type.
    !
    ! currentPatch%nrad(:,:) is a two dimensional array that indicates
    ! the total number of EXPOSED leaf layers, but for all intents and purposes
    ! in the photosynthesis routine, this appears to be the same as %ncan...
    !
    ! currentPatch%canopy_mask(:,:) has the same dimensions, is binary, and
    ! indicates whether or not leaf layers are present (by evaluating the canopy area
    ! profile).
    ! ---------------------------------------------------------------------------------

    ! Arguments
    type(fates_patch_type), target :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    ! Locals
    integer :: cl  ! Canopy Layer Index
    integer :: ft  ! Function Type Index
    integer :: iv  ! index of the exposed leaf layer for each canopy layer and pft

    ! Loop through the cohorts in this patch, associate each cohort with a layer and PFT
    ! and use the cohort's memory of how many layer's it takes up to assign the maximum
    ! of the layer/pft index it is in
    ! ---------------------------------------------------------------------------------

    currentPatch%nleaf(:,:) = 0
    ! redo the canopy structure algorithm to get round a
    ! bug that is happening for site 125, FT13.
    currentCohort => currentPatch%tallest
    do while(associated(currentCohort))

       currentPatch%nleaf(currentCohort%canopy_layer,currentCohort%pft) = &
            max(currentPatch%nleaf(currentCohort%canopy_layer,currentCohort%pft), &
            currentCohort%NV)

       currentCohort => currentCohort%shorter

    enddo !cohort

    ! NRAD = NCAN ...
    currentPatch%nrad = currentPatch%nleaf

    ! Now loop through and identify which layer and pft combo has scattering elements
    do cl = 1,nclmax
       do ft = 1,numpft
          currentPatch%canopy_mask(cl,ft) = 0
          do iv = 1, currentPatch%nrad(cl,ft);
             if(currentPatch%canopy_area_profile(cl,ft,iv) > 0._r8)then
                currentPatch%canopy_mask(cl,ft) = 1
             end if
          end do !iv
       enddo !ft
    enddo !cl

    return
  end subroutine UpdateCanopyNCanNRadPresent
  
  ! =====================================================================================
  


  real(r8) function ConvertPar(leaf_area, par_wm2) result(par_umolm2s)
    !
    ! DESCRIPTION:
    ! Convert par from W/m2 to umol photons/m2leaf/s
    !

    ! ARGUMENTS:
    real(r8), intent(in) :: leaf_area ! leaf area index [m2 leaf / m2 ground]
    real(r8), intent(in) :: par_wm2   ! absorbed PAR [W/m2 ground]
    
    ! minimum Leaf area to solve, too little has shown instability
    real(r8), parameter :: min_la_to_solve = 0.0000000001_r8

    if (par_wm2 > nearzero .and. leaf_area > min_la_to_solve) then
       par_umolm2s = par_wm2/leaf_area*wm2_to_umolm2s
    else                 
       ! The radiative transfer schemes are imperfect
       ! they can sometimes generate negative values here if par or leaf area is 0.0
       par_umolm2s = 0.0_r8
    end if
     
  end function ConvertPar

end module FATESPlantRespPhotosynthMod
