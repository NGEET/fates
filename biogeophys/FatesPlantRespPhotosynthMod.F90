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
  use FatesConstantsMod, only : molar_mass_ratio_vapdry
  use FatesConstantsMod, only : molar_mass_water
  use FatesConstantsMod, only : rgas_J_K_mol
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
  use FatesConstantsMod, only : nocomp_bareground
  use FatesConstantsMod, only : photosynth_acclim_model_none
  use FatesConstantsMod, only : photosynth_acclim_model_kumarathunge_etal_2019
  use FatesInterfaceTypesMod, only : hlm_use_planthydro
  use FatesInterfaceTypesMod, only : hlm_parteh_mode
  use FatesInterfaceTypesMod, only : numpft
  use FatesInterfaceTypesMod, only : nleafage
  use FatesUtilsMod,          only : QuadraticRoots => QuadraticRootsSridharachary
  use EDParamsMod,           only : maxpft
  use EDParamsMod,       only : nlevleaf
  use EDParamsMod,       only : nclmax
  use PRTGenericMod,     only : max_nleafage
  use EDTypesMod,        only : do_fates_salinity
  use EDParamsMod,       only : q10_mr
  use FatesPatchMod,     only : fates_patch_type
  use FatesCohortMod,    only : fates_cohort_type
  use EDParamsMod,       only : maintresp_leaf_model
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
  use EDParamsMod,       only : stomatal_model
  use EDParamsMod,       only : stomatal_assim_model
  use EDParamsMod,       only : dayl_switch
  use EDParamsMod,       only : photo_tempsens_model
  use PRTParametersMod,  only : prt_params
  use EDPftvarcon      , only : EDPftvarcon_inst
  use TemperatureType,   only : temperature_type
  use FatesRadiationMemMod, only : norman_solver,twostr_solver
  use EDParamsMod,          only : radiation_model
  use FatesRadiationMemMod, only : ipar
  use FatesTwoStreamUtilsMod, only : FatesGetCohortAbsRad
  use FatesAllometryMod     , only : VegAreaLayer
  use FatesAllometryMod, only : decay_coeff_vcmax
  
  ! CIME Globals
  use shr_log_mod , only      : errMsg => shr_log_errMsg

  implicit none
  private

  public :: FatesPlantRespPhotosynthDrive ! Called by the HLM-Fates interface

  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  character(len=1024) :: warn_msg   ! for defining a warning message

  !-------------------------------------------------------------------------------------

  ! maximum stomatal resistance [s/m] (used across several procedures)
  real(r8),parameter :: rsmax0 =  2.e8_r8

  logical   ::  debug = .false.
  !-------------------------------------------------------------------------------------

  ! Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
  real(r8),parameter :: h2o_co2_stoma_diffuse_ratio = 1.6_r8

  ! Ratio of H2O/CO2 gass diffusion in the leaf boundary layer (approximate)
  real(r8),parameter :: h2o_co2_bl_diffuse_ratio = 1.4_r8

  ! Constants used to define C3 versus C4 photosynth pathways
  integer, parameter :: c3_path_index = 1
  integer, parameter :: c4_path_index = 0


  ! Constants used to define conductance models
  integer, parameter :: medlyn_model = 2
  integer, parameter :: ballberry_model = 1

  ! Alternatively, Gross Assimilation can be used to estimate
  ! leaf co2 partial pressure and therefore conductance. The default
  ! is to use anet
  integer, parameter :: net_assim_model = 1
  integer, parameter :: gross_assim_model = 2

  logical, parameter :: preserve_b4b = .true.

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
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol
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

    ! Photosynthesis [umol /m2 /s]
    real(r8) :: psn_z(nlevleaf,maxpft,nclmax)
    
    ! Mask used to determine which leaf-layer biophysical rates have been
    ! used already
    logical :: rate_mask_z(nlevleaf,maxpft,nclmax)

    real(r8) :: vcmax_z            ! leaf layer maximum rate of carboxylation
    ! (umol co2/m**2/s)
    real(r8) :: jmax_z             ! leaf layer maximum electron transport rate
    ! (umol electrons/m**2/s)
    real(r8) :: kp_z               ! leaf layer initial slope of CO2 response
    ! curve (C4 plants)
    real(r8) :: c13disc_z(nclmax,maxpft,nlevleaf) ! carbon 13 in newly assimilated carbon at leaf level

    real(r8) :: mm_kco2            ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2             ! Michaelis-Menten constant for O2 (Pa)
    real(r8) :: co2_cpoint         ! CO2 compensation point (Pa)
    real(r8) :: btran_eff          ! effective transpiration wetness factor (0 to 1)
    real(r8) :: stomatal_intercept_btran   ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8) :: kn                 ! leaf nitrogen decay coefficient
    real(r8) :: cf                 ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
    real(r8) :: gb_mol             ! leaf boundary layer conductance (molar form: [umol /m**2/s])
    real(r8) :: ceair              ! vapor pressure of air, constrained (Pa)
    real(r8) :: nscaler            ! leaf nitrogen scaling coefficient
    real(r8) :: leaf_frac          ! ratio of to leaf biomass to total alive biomass
    real(r8) :: tcsoi              ! Temperature response function for root respiration.
    real(r8) :: tcwood             ! Temperature response function for wood

    real(r8) :: patch_la           ! exposed leaf area (patch scale)
    real(r8) :: live_stem_n        ! Live stem (above-ground sapwood)
    ! nitrogen content (kgN/plant)
    real(r8) :: live_croot_n       ! Live coarse root (below-ground sapwood)
    ! nitrogen content (kgN/plant)
    real(r8) :: sapw_c             ! Sapwood carbon (kgC/plant)
    real(r8) :: store_c_target     ! Target storage carbon (kgC/plant)
    real(r8) :: fnrt_c             ! Fine root carbon (kgC/plant)
    real(r8) :: fnrt_n             ! Fine root nitrogen content (kgN/plant)
    real(r8) :: leaf_c             ! Leaf carbon (kgC/plant)
    real(r8) :: leaf_n             ! leaf nitrogen content (kgN/plant)
    real(r8) :: g_sb_leaves        ! Mean combined (stomata+boundary layer) leaf conductance [m/s]
    ! over all of the patch's leaves.  The "sb" refers to the combined
    ! "s"tomatal and "b"oundary layer.
    ! This quantity is relevant on leaf surfaces. It does not
    ! have units of /m2 leaf per say, but is implicitly on leaf surfaces
    real(r8) :: r_sb_leaves        ! Mean leaf resistance over all the patch's leaves [s/m]
    ! This is the direct reciprocal of g_sb_leaves
    real(r8) :: r_stomata          ! Mean stomatal resistance across all leaves in the patch [s/m]


    real(r8) :: maintresp_reduction_factor  ! factor by which to reduce maintenance
    ! respiration when storage pools are low
    real(r8) :: b_leaf             ! leaf biomass kgC
    real(r8) :: frac               ! storage pool as a fraction of target leaf biomass
    ! over each cohort x layer.
    real(r8) :: cohort_eleaf_area  ! This is the effective leaf area [m2] reported by each cohort
    real(r8) :: lnc_top            ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8) :: lmr25top           ! canopy top leaf maint resp rate at 25C 
    ! for this plant or pft (umol CO2/m**2/s)
    real(r8) :: leaf_inc           ! LAI-only portion of the vegetation increment of dinc_vai
    real(r8) :: lai_canopy_above   ! the LAI in the canopy layers above the layer of interest
    real(r8) :: lai_layers_above   ! the LAI in the leaf layers, within the current canopy,
    ! above the leaf layer of interest
    real(r8) :: lai_current        ! the LAI in the current leaf layer
    real(r8) :: cumulative_lai     ! the cumulative LAI, top down, to the leaf layer of interest
    real(r8) :: leaf_psi           ! leaf xylem matric potential [MPa] (only meaningful/used w/ hydro)
    real(r8) :: fnrt_mr_layer      ! fine root maintenance respiation per layer [kgC/plant/s]

    real(r8) :: fnrt_mr_nfix_layer ! fineroot maintenance respiration specifically for symbiotic fixation [kgC/plant/layer/s]
    real(r8) :: nfix_layer         ! Nitrogen fixed in each layer this timestep [kgN/plant/layer/timestep]
    real(r8), allocatable :: rootfr_ft(:,:)  ! Root fractions per depth and PFT

    real(r8) :: agb_frac              ! fraction of biomass aboveground
    real(r8) :: branch_frac           ! fraction of aboveground woody biomass in branches
    real(r8) :: crown_reduction       ! reduction in crown biomass from damage
    real(r8) :: sapw_c_bgw            ! belowground sapwood
    real(r8) :: sapw_c_agw            ! aboveground sapwood
    real(r8) :: sapw_c_undamaged      ! the target sapwood of an undamaged tree
    real(r8) :: sapw_n                ! sapwood nitrogen
    real(r8) :: sapw_n_bgw            ! nitrogen in belowground portion of sapwood
    real(r8) :: sapw_n_agw            ! nitrogen in aboveground portion of sapwood
    real(r8) :: sapw_n_undamaged      ! nitrogen in sapwood of undamaged tree
    real(r8) :: rd_abs_leaf, rb_abs_leaf, r_abs_stem, r_abs_snow, rb_abs, rd_abs
    real(r8) :: fsun
    real(r8) :: par_per_sunla, par_per_shala ! PAR per sunlit and shaded leaf area [W/m2 leaf]
    real(r8),dimension(75) :: cohort_vaitop
    real(r8),dimension(75) :: cohort_vaibot
    real(r8),dimension(75) :: cohort_layer_elai
    real(r8),dimension(75) :: cohort_layer_esai
    real(r8),dimension(75) :: cohort_layer_tlai
    real(r8),dimension(75) :: cohort_layer_tsai
    real(r8)               :: cohort_elai
    real(r8)               :: cohort_esai
    real(r8)               :: laisun,laisha
    real(r8)               :: canopy_area
    real(r8)               :: elai
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

    ! Parameters
    !
    ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
    ! (gC/gN/s)
    ! ------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------
    ! Photosynthesis and stomatal conductance parameters, from:
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    ! -----------------------------------------------------------------------------------



    associate(  &
         c3psn     => EDPftvarcon_inst%c3psn  , &
         slatop    => prt_params%slatop , &  ! specific leaf area at top of canopy,
                                ! projected area basis [m^2/gC]
         woody     => prt_params%woody,   &  ! Is vegetation woody or not?
         stomatal_intercept   => EDPftvarcon_inst%stomatal_intercept ) !Unstressed minimum stomatal conductance


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

               psn_z(:,:,:) = 0._r8
               
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
                  ! currentPatch%nleaf(:,:)
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
                       bc_in(s)%tgcm_pa(ifp),    & ! in
                       bc_in(s)%eair_pa(ifp),    & ! in
                       bc_in(s)%esat_tv_pa(ifp), & ! in
                       bc_in(s)%rb_pa(ifp),      & ! in
                       mm_kco2,                  & ! out
                       mm_ko2,                   & ! out
                       co2_cpoint,               & ! out
                       cf,                       & ! out
                       gb_mol,                   & ! out
                       ceair)                      ! out

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

                  if_any_cohorts: if(currentPatch%num_cohorts > 0.0)then
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
                        call lowstorage_maintresp_reduction(frac,currentCohort%pft, &
                             maintresp_reduction_factor)

                        ! are there any leaves of this pft in this layer?
                        canopy_mask_if: if(currentPatch%canopy_mask(cl,ft) == 1)then

                           ! Loop over leaf-layers
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
                                   (radiation_model .eq. twostr_solver ) .or. &
                                   (nleafage > 1) .or. &
                                   (hlm_parteh_mode .ne. prt_carbon_allom_hyp )   ) then

                                 if (hlm_use_planthydro.eq.itrue ) then

                                    stomatal_intercept_btran = max( cf/rsmax0,stomatal_intercept(ft)*currentCohort%co_hydr%btran )
                                    btran_eff = currentCohort%co_hydr%btran 

                                    ! dinc_vai(:) is the total vegetation area index of each "leaf" layer
                                    ! we convert to the leaf only portion of the increment
                                    ! ------------------------------------------------------
                                    leaf_inc    = dinc_vai(iv) * &
                                         currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)

                                    ! Now calculate the cumulative top-down lai of the current layer's midpoint
                                    lai_canopy_above  = sum(currentPatch%canopy_layer_tlai(1:cl-1)) 

                                    lai_layers_above  = (dlower_vai(iv) - dinc_vai(iv)) * &
                                         currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
                                    lai_current       = min(leaf_inc, currentCohort%treelai - lai_layers_above)
                                    cumulative_lai    = lai_canopy_above + lai_layers_above + 0.5*lai_current 

                                    leaf_psi = currentCohort%co_hydr%psi_ag(1)

                                 else

                                    stomatal_intercept_btran = max( cf/rsmax0,stomatal_intercept(ft)*currentPatch%btran_ft(ft) ) 

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

                                 kn = decay_coeff_vcmax(currentCohort%vcmax25top, &
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

                                 select case (maintresp_leaf_model)

                                 case (lmrmodel_ryan_1991)

                                    call LeafLayerMaintenanceRespiration_Ryan_1991( lnc_top,     &  ! in
                                         nscaler,                  &  ! in
                                         ft,                       &  ! in
                                         bc_in(s)%t_veg_pa(ifp),   &  ! in
                                         lmr_z(iv,ft,cl))             ! out

                                 case (lmrmodel_atkin_etal_2017)

                                    call LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top, &  ! in
                                         cumulative_lai,                     &  ! in
                                         currentCohort%vcmax25top,           &  ! in
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

                                 if_radsolver: if(radiation_model.eq.norman_solver) then

                                    laisun = currentPatch%ed_laisun_z(cl,ft,iv)
                                    laisha = currentPatch%ed_laisha_z(cl,ft,iv)
                                    par_per_sunla = currentPatch%ed_parsun_z(cl,ft,iv)
                                    par_per_shala = currentPatch%ed_parsha_z(cl,ft,iv)
                                    canopy_area   = currentPatch%canopy_area_profile(cl,ft,iv)
                                    fsun = currentPatch%f_sun(cl,ft,iv)
                                    
                                 else    ! Two-stream

                                    if(cohort_layer_elai(iv) > nearzero .and. currentPatch%solar_zenith_flag) then

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

                                 ! Part VII: Calculate (1) maximum rate of carboxylation (vcmax),
                                 ! (2) maximum electron transport rate, (3) triose phosphate
                                 ! utilization rate and (4) the initial slope of CO2 response curve
                                 ! (C4 plants). Earlier we calculated their base rates as dictated
                                 ! by their plant functional type and some simple scaling rules for
                                 ! nitrogen limitation baesd on canopy position (not prognostic).
                                 ! These rates are the specific rates used in the actual photosynthesis
                                 ! calculations that take localized environmental effects (temperature)
                                 ! into consideration.

                                 call LeafLayerBiophysicalRates(par_per_sunla, & ! in
                                      ft,                                 &  ! in
                                      currentCohort%vcmax25top,           &  ! in
                                      currentCohort%jmax25top,            &  ! in
                                      currentCohort%kp25top,              &  ! in
                                      nscaler,                            &  ! in
                                      bc_in(s)%t_veg_pa(ifp),             &  ! in
                                      bc_in(s)%dayl_factor_pa(ifp),       &  ! in
                                      currentPatch%tveg_lpa%GetMean(),    &  ! in
                                      currentPatch%tveg_longterm%GetMean(),&  ! in
                                      btran_eff,                          &  ! in
                                      vcmax_z,                            &  ! out
                                      jmax_z,                             &  ! out
                                      kp_z )                                 ! out

                                 ! Part IX: This call calculates the actual photosynthesis for the
                                 ! leaf layer, as well as the stomatal resistance and the net assimilated carbon.

                                 call LeafLayerPhotosynthesis(fsun,       &  ! in
                                      par_per_sunla,                      &  ! in
                                      par_per_shala,                      &  ! in
                                      laisun,                             &  ! in
                                      laisha,                             &  ! in
                                      canopy_area,                        &  ! in
                                      ft,                                 &  ! in
                                      vcmax_z,                            &  ! in
                                      jmax_z,                             &  ! in
                                      kp_z,                               &  ! in
                                      bc_in(s)%t_veg_pa(ifp),             &  ! in
                                      bc_in(s)%esat_tv_pa(ifp),           &  ! in
                                      bc_in(s)%forc_pbot,                 &  ! in
                                      bc_in(s)%cair_pa(ifp),              &  ! in
                                      bc_in(s)%oair_pa(ifp),              &  ! in
                                      btran_eff,                          &  ! in
                                      stomatal_intercept_btran,           &  ! in
                                      cf,                                 &  ! in
                                      gb_mol,                             &  ! in
                                      ceair,                              &  ! in
                                      mm_kco2,                            &  ! in
                                      mm_ko2,                             &  ! in
                                      co2_cpoint,                         &  ! in
                                      lmr_z(iv,ft,cl),                    &  ! in
                                      leaf_psi,                           &  ! in
                                      bc_in(s)%rb_pa(ifp),                &  ! in  
                                      psn_z(iv,ft,cl),                    &  ! out
                                      rs_z(iv,ft,cl),                     &  ! out
                                      anet_av_z(iv,ft,cl),                &  ! out
                                      c13disc_z(cl,ft,iv))                   ! out

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

                           ! Temporary bypass to preserve B4B behavior
                           if(radiation_model.eq.norman_solver) then

                              call ScaleLeafLayerFluxToCohort(nv,                                    & !in
                                   psn_z(1:nv,ft,cl),                     & !in
                                   lmr_z(1:nv,ft,cl),                     & !in
                                   rs_z(1:nv,ft,cl),                      & !in
                                   currentPatch%elai_profile(cl,ft,1:nv), & !in
                                   c13disc_z(cl, ft, 1:nv),               & !in
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

                              call ScaleLeafLayerFluxToCohort(nv,                                    & !in
                                   psn_z(1:nv,ft,cl),                     & !in
                                   lmr_z(1:nv,ft,cl),                     & !in
                                   rs_z(1:nv,ft,cl),                      & !in
                                   cohort_layer_elai(1:nv),               & !in
                                   c13disc_z(cl, ft, 1:nv),               & !in
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
                  ! preserve_b4b will be removed soon. This is kept here to prevent
                  ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23) 
                  if(preserve_b4b) then
                     patch_la = patch_la/ currentPatch%total_canopy_area
                  end if
                  
                  ! Normalize canopy total conductance by the effective LAI
                  ! The value here was integrated over each cohort x leaf layer
                  ! and was weighted by m2 of effective leaf area for each layer
                  
                  if_any_lai: if(patch_la>tiny(patch_la)) then

                     ! Normalize the leaf-area weighted canopy conductance
                     ! The denominator is the total effective leaf area in the canopy,
                     ! units of [m/s]*[m2] / [m2] = [m/s]
                     ! preserve_b4b will be removed soon. This is kept here to prevent
                     ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23) 
                     if_preserve_b4b3: if(preserve_b4b) then
                        elai     = calc_areaindex(currentPatch,'elai')
                        g_sb_leaves = g_sb_leaves / (elai*currentPatch%total_canopy_area)
                     else
                        g_sb_leaves = g_sb_leaves / max(0.1_r8*currentPatch%total_canopy_area,patch_la)
                     end if if_preserve_b4b3
                     
                     
                     if_above_mincond: if( g_sb_leaves > (1._r8/rsmax0) ) then
                        
                        ! Combined mean leaf resistance is the inverse of mean leaf conductance
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
                     currentPatch%c_stomata  = cf / r_stomata
                     
                  else !if_any_lai
                     
                     ! But this will prevent it from using an unintialized value
                     bc_out(s)%rssun_pa(ifp) = rsmax0
                     bc_out(s)%rssha_pa(ifp) = rsmax0
                     
                     ! This value is used for diagnostics, the molar form of conductance
                     ! is what is used in the field usually, so we track that form
                     currentPatch%c_stomata  = cf / rsmax0
                     
                  end if if_any_lai

                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_lblayer = cf / bc_in(s)%rb_pa(ifp)

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

subroutine LeafLayerPhotosynthesis(f_sun_lsl,         &  ! in
     parsun_lsl,        &  ! in
     parsha_lsl,        &  ! in
     laisun_lsl,        &  ! in
     laisha_lsl,        &  ! in
     canopy_area_lsl,   &  ! in
     ft,                &  ! in
     vcmax,             &  ! in
     jmax,              &  ! in
     co2_rcurve_islope, &  ! in
     veg_tempk,         &  ! in
     veg_esat,          &  ! in
     can_press,         &  ! in
     can_co2_ppress,    &  ! in
     can_o2_ppress,     &  ! in
     btran,             &  ! in
     stomatal_intercept_btran,  &  ! in
     cf,                &  ! in
     gb_mol,            &  ! in
     ceair,             &  ! in
     mm_kco2,           &  ! in
     mm_ko2,            &  ! in
     co2_cpoint,        &  ! in
     lmr,               &  ! in
     leaf_psi,          &  ! in
     rb,                &  ! in
     psn_out,           &  ! out
     rstoma_out,        &  ! out
     anet_av_out,       &  ! out
     c13disc_z)            ! out


  ! ------------------------------------------------------------------------------------
  ! This subroutine calculates photosynthesis and stomatal conductance within each leaf
  ! sublayer.
  ! A note on naming conventions: As this subroutine is called for every
  ! leaf-sublayer, many of the arguments are specific to that "leaf sub layer"
  ! (LSL), those variables are given a dimension tag "_lsl"
  ! Other arguments or variables may be indicative of scales broader than the LSL.
  ! ------------------------------------------------------------------------------------

  use EDParamsMod       , only : theta_cj_c3, theta_cj_c4


  ! Arguments
  ! ------------------------------------------------------------------------------------
  real(r8), intent(in) :: f_sun_lsl         !
  real(r8), intent(in) :: parsun_lsl        ! Absorbed PAR in sunlist leaves
  real(r8), intent(in) :: parsha_lsl        ! Absorved PAR in shaded leaves
  real(r8), intent(in) :: laisun_lsl        ! LAI in sunlit leaves
  real(r8), intent(in) :: laisha_lsl        ! LAI in shaded leaves
  real(r8), intent(in) :: canopy_area_lsl   !
  integer,  intent(in) :: ft                ! (plant) Functional Type Index
  real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
  real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
  real(r8), intent(in) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
  real(r8), intent(in) :: veg_tempk         ! vegetation temperature
  real(r8), intent(in) :: veg_esat          ! saturation vapor pressure at veg_tempk (Pa)

  ! Important Note on the following gas pressures.  This photosynthesis scheme will iteratively
  ! solve for the co2 partial pressure at the leaf surface (ie in the stomata). The reference
  ! point for these input values are NOT within that boundary layer that separates the stomata from
  ! the canopy air space.  The reference point for these is on the outside of that boundary
  ! layer.  This routine, which operates at the leaf scale, makes no assumptions about what the
  ! scale of the refernce is, it could be lower atmosphere, it could be within the canopy
  ! but most likely it is the closest value one can get to the edge of the leaf's boundary
  ! layer.  We use the convention "can_" because a reference point of within the canopy
  ! ia a best reasonable scenario of where we can get that information from.

  real(r8), intent(in) :: can_press       ! Air pressure NEAR the surface of the leaf (Pa)
  real(r8), intent(in) :: can_co2_ppress  ! Partial pressure of CO2 NEAR the leaf surface (Pa)
  real(r8), intent(in) :: can_o2_ppress   ! Partial pressure of O2 NEAR the leaf surface (Pa)
  real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1)
  real(r8), intent(in) :: stomatal_intercept_btran !water-stressed minimum stomatal conductance (umol H2O/m**2/s)
  real(r8), intent(in) :: cf              ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
  real(r8), intent(in) :: gb_mol          ! leaf boundary layer conductance (umol /m**2/s)
  real(r8), intent(in) :: ceair           ! vapor pressure of air, constrained (Pa)
  real(r8), intent(in) :: mm_kco2         ! Michaelis-Menten constant for CO2 (Pa)
  real(r8), intent(in) :: mm_ko2          ! Michaelis-Menten constant for O2 (Pa)
  real(r8), intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
  real(r8), intent(in) :: lmr             ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
  real(r8), intent(in) :: leaf_psi        ! Leaf water potential [MPa]
  real(r8), intent(in) :: rb              ! Boundary Layer resistance of leaf [s/m]
  
  real(r8), intent(out) :: psn_out        ! carbon assimilated in this leaf layer umolC/m2/s
  real(r8), intent(out) :: rstoma_out     ! stomatal resistance (1/gs_lsl) (s/m)
  real(r8), intent(out) :: anet_av_out    ! net leaf photosynthesis (umol CO2/m**2/s)
  ! averaged over sun and shade leaves.
  real(r8), intent(out) :: c13disc_z      ! carbon 13 in newly assimilated carbon

 

  
  ! Locals
  ! ------------------------------------------------------------------------
  integer :: c3c4_path_index    ! Index for which photosynthetic pathway
  ! is active.  C4 = 0,  C3 = 1
  integer :: sunsha             ! Index for differentiating sun and shade
  real(r8) :: gstoma            ! Stomatal Conductance of this leaf layer (m/s)
  real(r8) :: agross            ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
  real(r8) :: anet              ! net leaf photosynthesis (umol CO2/m**2/s)
  real(r8) :: a_gs              ! The assimilation (a) for calculating conductance (gs)
                                ! is either = to anet or agross
  real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
  real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
  real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
  real(r8) :: r1,r2             ! roots of quadratic equation
  real(r8) :: co2_inter_c       ! intercellular leaf CO2 (Pa)
  real(r8) :: co2_inter_c_old   ! intercellular leaf CO2 (Pa) (previous iteration)
  logical  :: loop_continue     ! Loop control variable
  integer  :: niter             ! iteration loop index
  real(r8) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
  real(r8) :: gs                ! leaf stomatal conductance (m/s)
  real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
  real(r8) :: gs_mol_err        ! gs_mol for error check
  real(r8) :: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: ap                ! product-limited (C3) or CO2-limited
  ! (C4) gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
  real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
  real(r8) :: init_co2_inter_c  ! First guess intercellular co2 specific to C path
  real(r8) :: term                 ! intermediate variable in Medlyn stomatal conductance model
  real(r8) :: vpd                  ! water vapor deficit in Medlyn stomatal model (KPa)


  ! Parameters
  ! ------------------------------------------------------------------------
  ! Fraction of light absorbed by non-photosynthetic pigments
  real(r8),parameter :: fnps = 0.15_r8

  ! term accounting that two photons are needed to fully transport a single 
  ! electron in photosystem 2
  real(r8), parameter :: photon_to_e = 0.5_r8

  ! Unit conversion of w/m2 to umol photons m-2 s-1
  real(r8), parameter :: wm2_to_umolm2s = 4.6_r8
  
  ! For plants with no leaves, a miniscule amount of conductance
  ! can happen through the stems, at a partial rate of cuticular conductance
  real(r8),parameter :: stem_cuticle_loss_frac = 0.1_r8

  ! empirical curvature parameter for electron transport rate
  real(r8),parameter :: theta_psii = 0.7_r8

  ! First guess on ratio between intercellular co2 and the atmosphere
  ! an iterator converges on actual
  real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
  real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

  ! quantum efficiency, used only for C4 (mol CO2 / mol photons) (index 0)
  real(r8),parameter,dimension(0:1) :: quant_eff = [0.05_r8,0.0_r8]

  ! empirical curvature parameter for ap photosynthesis co-limitation
  real(r8),parameter :: theta_ip = 0.999_r8

  ! minimum Leaf area to solve, too little has shown instability
  real(r8), parameter :: min_la_to_solve = 0.0000000001_r8
  
  associate( bb_slope  => EDPftvarcon_inst%bb_slope      ,& ! slope of BB relationship, unitless
       medlyn_slope=> EDPftvarcon_inst%medlyn_slope          , & ! Slope for Medlyn stomatal conductance model method, the unit is KPa^0.5
       stomatal_intercept=> EDPftvarcon_inst%stomatal_intercept )  !Unstressed minimum stomatal conductance, the unit is umol/m**2/s

  ! photosynthetic pathway: 0. = c4, 1. = c3
  c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ft))

  if (c3c4_path_index == c3_path_index) then
     init_co2_inter_c = init_a2l_co2_c3 * can_co2_ppress
  else
     init_co2_inter_c = init_a2l_co2_c4 * can_co2_ppress
  end if

  ! Part III: Photosynthesis and Conductance
  ! ----------------------------------------------------------------------------------

  if_daytime: if ( parsun_lsl <= 0._r8 ) then  ! night time
 
     anet_av_out = -lmr
     psn_out     = 0._r8

     ! The cuticular conductance already factored in maximum resistance as a bound
     ! no need to re-bound it

     rstoma_out = cf/stomatal_intercept_btran

     c13disc_z = 0.0_r8    !carbon 13 discrimination in night time carbon flux, note value of 1.0 is used in CLM

  else ! day time (a little bit more complicated ...)

     ! Is there leaf area? - (NV can be larger than 0 with only stem area if deciduous)
     
     if_leafarea: if ( laisun_lsl + laisha_lsl > 0._r8 ) then

        !Loop aroun shaded and unshaded leaves
        psn_out     = 0._r8    ! psn is accumulated across sun and shaded leaves.
        rstoma_out  = 0._r8    ! 1/rs is accumulated across sun and shaded leaves.
        anet_av_out = 0._r8
        gstoma  = 0._r8

        do  sunsha = 1,2
           ! Electron transport rate for C3 plants.
           ! Convert par from W/m2 to umol photons/m**2/s
           ! Convert from units of par absorbed per unit ground area to par
           ! absorbed per unit leaf area.

           if(sunsha == 1)then !sunlit
              if(( laisun_lsl * canopy_area_lsl) > min_la_to_solve)then

                 qabs = parsun_lsl / (laisun_lsl * canopy_area_lsl )
                 qabs = qabs * photon_to_e * (1._r8 - fnps) * wm2_to_umolm2s

              else
                 qabs = 0.0_r8
              end if
           else

              if( (parsha_lsl>nearzero) .and. (laisha_lsl * canopy_area_lsl) > min_la_to_solve  ) then

                 qabs = parsha_lsl / (laisha_lsl * canopy_area_lsl)
                 qabs = qabs * photon_to_e * (1._r8 - fnps) *  wm2_to_umolm2s
              else                 
                 ! The radiative transfer schemes are imperfect
                 ! they can sometimes generate negative values here
                 qabs = 0._r8
              end if
              
           end if

           !convert the absorbed par into absorbed par per m2 of leaf,
           ! so it is consistant with the vcmax and lmr numbers.
           aquad = theta_psii
           bquad = -(qabs + jmax)
           cquad = qabs * jmax
           call QuadraticRoots(aquad, bquad, cquad, r1, r2)
           je = min(r1,r2)

           ! Initialize intercellular co2
           co2_inter_c = init_co2_inter_c

           niter = 0
           loop_continue = .true.
           iter_loop: do while(loop_continue)
              ! Increment iteration counter. Stop if too many iterations
              niter = niter + 1

              ! Save old co2_inter_c
              co2_inter_c_old = co2_inter_c

              ! Photosynthesis limitation rate calculations
              if (c3c4_path_index == c3_path_index)then

                 ! C3: Rubisco-limited photosynthesis
                 ac = vcmax * max(co2_inter_c-co2_cpoint, 0._r8) / &
                      (co2_inter_c+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))

                 ! C3: RuBP-limited photosynthesis
                 aj = je * max(co2_inter_c-co2_cpoint, 0._r8) / &
                      (4._r8*co2_inter_c+8._r8*co2_cpoint)

                    ! Gross photosynthesis smoothing calculations. Co-limit ac and aj.
                    aquad = theta_cj_c3
                    bquad = -(ac + aj)
                    cquad = ac * aj
                    call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                    agross = min(r1,r2)

              else

                 ! C4: Rubisco-limited photosynthesis
                 ac = vcmax

                 ! C4: RuBP-limited photosynthesis
                 if(sunsha == 1)then !sunlit
                    !guard against /0's in the night.
                    if((laisun_lsl * canopy_area_lsl) > 0.0000000001_r8) then
                       aj = quant_eff(c3c4_path_index) * parsun_lsl * wm2_to_umolm2s
                       !convert from per cohort to per m2 of leaf)
                       aj = aj / (laisun_lsl * canopy_area_lsl)
                    else
                       aj = 0._r8
                    end if
                 else
                    aj = quant_eff(c3c4_path_index) * parsha_lsl * wm2_to_umolm2s
                    aj = aj / (laisha_lsl * canopy_area_lsl)
                 end if

                 ! C4: PEP carboxylase-limited (CO2-limited)
                 ap = co2_rcurve_islope * max(co2_inter_c, 0._r8) / can_press

                 ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap

                 aquad = theta_cj_c4
                 bquad = -(ac + aj)
                 cquad = ac * aj
                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 ai = min(r1,r2)

                 aquad = theta_ip
                 bquad = -(ai + ap)
                 cquad = ai * ap
                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 agross = min(r1,r2)

              end if
              
              ! Calculate anet, only exit iteration with negative anet when
              ! using anet in calculating gs this is version B  
              anet = agross  - lmr

              if ( stomatal_assim_model == gross_assim_model ) then
                 if ( stomatal_model == medlyn_model ) then
                    write (fates_log(),*) 'Gross Assimilation conductance is incompatible with the Medlyn model'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                 end if
                 a_gs = agross
              else
                 if (anet < 0._r8) then
                    loop_continue = .false.
                 end if
                 a_gs = anet
              end if

              ! With an <= 0, then gs_mol = stomatal_intercept_btran
              leaf_co2_ppress = can_co2_ppress- h2o_co2_bl_diffuse_ratio/gb_mol * a_gs * can_press 		   
              leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)

              ! A note about the use of the quadratic equations for calculating stomatal conductance
              ! ------------------------------------------------------------------------------------
              ! These two following models calculate the conductance between the intercellular leaf
              ! space and the leaf surface, not the canopy air space.  Transport between the leaf
              ! surface and the canopy air space is governed by the leaf boundary layer conductance.
              ! However, we need to estimate the properties at the surface of the leaf to solve for
              ! the stomatal conductance. We do this by using Fick's law (gradient resistance
              ! approximation of diffusion) to estimate the flux of water vapor across the
              ! leaf boundary layer, and balancing that with the flux across the stomata. It
              ! results in the following equation for leaf surface humidity:
              !
              ! e_s = (e_i g_s + e_c g_b)/(g_b + g_s)
              !
              ! The leaf surface humidity (e_s) becomes an expression of canopy humidity (e_c),
              ! intercellular humidity (e_i, which is the saturation humidity at leaf temperature),
              ! boundary layer conductance (g_b) (these are all known) and stomatal conductance
              ! (g_s) (this is still unknown).  This expression is substituted into the stomatal
              ! conductance equation. The resulting form of these equations becomes a quadratic.
              !
              ! For a detailed explanation, see the FATES technical note, section
              ! "1.11 Stomatal Conductance"
              !
              ! ------------------------------------------------------------------------------------
              
              
              if ( stomatal_model == medlyn_model ) then
                 !stomatal conductance calculated from Medlyn et al. (2011), the numerical &
                 !implementation was adapted from the equations in CLM5.0
                 vpd =  max((veg_esat - ceair), 50._r8) * 0.001_r8          !addapted from CLM5. Put some constraint on VPD
                 !when Medlyn stomatal conductance is being used, the unit is KPa. Ignoring the constraint will cause errors when model runs.
                 term = h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
                 aquad = 1.0_r8
                 bquad = -(2.0 * (stomatal_intercept_btran+ term) + (medlyn_slope(ft) * term)**2 / &
                      (gb_mol * vpd ))
                 cquad = stomatal_intercept_btran*stomatal_intercept_btran + &
                      (2.0*stomatal_intercept_btran + term * &
                      (1.0 - medlyn_slope(ft)* medlyn_slope(ft) / vpd)) * term

                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 gs_mol = max(r1,r2)

              else if ( stomatal_model == ballberry_model ) then         !stomatal conductance calculated from Ball et al. (1987)
                 aquad = leaf_co2_ppress
                 bquad = leaf_co2_ppress*(gb_mol - stomatal_intercept_btran) - bb_slope(ft) * a_gs * can_press
                 cquad = -gb_mol*(leaf_co2_ppress*stomatal_intercept_btran + &
                      bb_slope(ft)*anet*can_press * ceair/ veg_esat )

                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 gs_mol = max(r1,r2)
              end if
              
              ! Derive new estimate for co2_inter_c
              co2_inter_c = can_co2_ppress - anet * can_press * &
                   (h2o_co2_bl_diffuse_ratio*gs_mol+h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)

              ! Check for co2_inter_c convergence. Delta co2_inter_c/pair = mol/mol.
              ! Multiply by 10**6 to convert to umol/mol (ppm). Exit iteration if
              ! convergence criteria of +/- 1 x 10**-6 ppm is met OR if at least ten
              ! iterations (niter=10) are completed

              if ((abs(co2_inter_c-co2_inter_c_old)/can_press*1.e06_r8 <=  2.e-06_r8) &
                   .or. niter == 5) then
                 loop_continue = .false.
              end if
           end do iter_loop

           ! End of co2_inter_c iteration.  Check for an < 0, in which case gs_mol = bbb
           ! And Final estimates for leaf_co2_ppress and co2_inter_c 
           ! (needed for early exit of co2_inter_c iteration when an < 0)	 
           if (anet < 0._r8) then
              gs_mol = stomatal_intercept_btran
           end if

           ! Final estimates for leaf_co2_ppress and co2_inter_c
           leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb_mol * anet * can_press
           leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
           co2_inter_c = can_co2_ppress - anet * can_press * &
                (h2o_co2_bl_diffuse_ratio*gs_mol+h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)

           ! Convert gs_mol (umol /m**2/s) to gs (m/s) and then to rs (s/m)
           gs = gs_mol / cf

           ! estimate carbon 13 discrimination in leaf level carbon
           ! flux Liang WEI and Hang ZHOU 2018, based on
           ! Ubierna and Farquhar, 2014 doi:10.1111/pce.12346, using the simplified model:
           ! $\Delta ^{13} C = \alpha_s + (b - \alpha_s) \cdot \frac{C_i}{C_a}$
           ! just hard code b and \alpha_s for now, might move to parameter set in future
           ! b = 27.0 alpha_s = 4.4
           ! TODO, not considering C4 or CAM right now, may need to address this
           ! note co2_inter_c is intracelluar CO2, not intercelluar
           c13disc_z = 4.4_r8 + (27.0_r8 - 4.4_r8) * &
                min (can_co2_ppress, max (co2_inter_c, 0._r8)) / can_co2_ppress

           ! Accumulate total photosynthesis umol/m2 ground/s-1.
           ! weight per unit sun and sha leaves.
           if(sunsha == 1)then !sunlit
              psn_out     = psn_out + agross * f_sun_lsl
              anet_av_out = anet_av_out + anet * f_sun_lsl
              gstoma  = gstoma + 1._r8/(min(1._r8/gs, rsmax0)) * f_sun_lsl
           else
              psn_out = psn_out + agross * (1.0_r8-f_sun_lsl)
              anet_av_out = anet_av_out + anet * (1.0_r8-f_sun_lsl)
              gstoma  = gstoma + &
                   1._r8/(min(1._r8/gs, rsmax0)) * (1.0_r8-f_sun_lsl)
           end if

           ! Make sure iterative solution is correct
           if (gs_mol < 0._r8) then
              write (fates_log(),*)'Negative stomatal conductance:'
              write (fates_log(),*)'gs_mol= ',gs_mol
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if

           ! Compare with Medlyn model: gs_mol = 1.6*(1+m/sqrt(vpd)) * an/leaf_co2_ppress*p + b
           if ( stomatal_model == 2 ) then
              gs_mol_err = h2o_co2_stoma_diffuse_ratio*(1 + medlyn_slope(ft)/sqrt(vpd))*max(anet,0._r8)/leaf_co2_ppress*can_press + stomatal_intercept_btran
              ! Compare with Ball-Berry model: gs_mol = m * an * hs/leaf_co2_ppress*p + b
           else if ( stomatal_model == 1 ) then
              hs = (gb_mol*ceair + gs_mol* veg_esat ) / ((gb_mol+gs_mol)*veg_esat )
              gs_mol_err = bb_slope(ft)*max(anet, 0._r8)*hs/leaf_co2_ppress*can_press + stomatal_intercept_btran
           end if

           if (abs(gs_mol-gs_mol_err) > 1.e-01_r8) then
              warn_msg = 'Stomatal conductance error check - weak convergence: '//trim(N2S(gs_mol))//' '//trim(N2S(gs_mol_err))
              call FatesWarn(warn_msg,index=1)
           end if

        enddo !sunsha loop

        ! Stomatal resistance of the leaf-layer
        if ( (hlm_use_planthydro.eq.itrue .and. EDPftvarcon_inst%hydr_k_lwp(ft)>nearzero) ) then
           rstoma_out = LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press, veg_esat, &
                                               rb, gstoma, ft)
        else
           rstoma_out = 1._r8/gstoma
        end if
           
        
     else

        ! No leaf area. This layer is present only because of stems.
        ! Net assimilation is zero, not negative because there are
        ! no leaves to even respire
        ! (leaves are off, or have reduced to 0)

        psn_out     = 0._r8
        anet_av_out = 0._r8

        rstoma_out  = min(rsmax0,cf/(stem_cuticle_loss_frac*stomatal_intercept(ft)))
        c13disc_z = 0.0_r8

     end if if_leafarea !is there leaf area?


   end if if_daytime    ! night or day


   end associate
   return
  end subroutine LeafLayerPhotosynthesis

  ! =======================================================================================

  function LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press, veg_esat, &
       rb, gstoma, ft) result(rstoma_out)

    ! -------------------------------------------------------------------------------------
    ! This calculates inner leaf humidity as a function of mesophyll water potential 
    ! Adopted from  Vesala et al., 2017 https://www.frontiersin.org/articles/10.3389/fpls.2017.00054/full
    !
    ! Equation 1 in Vesala et al:
    ! lwp_star = wi/w0 = exp( k_lwp*leaf_psi*molar_mass_water/(rgas_J_k_mol * veg_tempk) )
    !
    ! Terms:
    ! leaf_psi: leaf water potential [MPa]
    ! k_lwp: inner leaf humidity scaling coefficient [-]
    ! rgas_J_K_mol: universal gas constant, [J/K/mol], 8.3144598
    ! molar_mass_water, molar mass of water, [g/mol]: 18.0
    !
    ! Unit conversions:
    ! 1 Pa = 1 N/m2 = 1 J/m3
    ! density of liquid water [kg/m3] = 1000
    ! 
    ! units of equation 1:  exp( [MPa]*[g/mol]/( [J/K/mol] * [K] ) )
    !                            [MJ/m3]*[g/mol]*[m3/kg]*[kg/g]*[J/MJ]  / ([J/mol])
    ! dimensionless:             [J/g]*[g/mol]/([J/mol])
    !
    ! Note: unit conversions drop out b/c [m3/kg]*[kg/g]*[J/MJ] = 1e-3*1.e-3*1e6 = 1.0
    !
    ! Junyan Ding 2021
    ! -------------------------------------------------------------------------------------

    ! Arguments
    real(r8) :: leaf_psi   ! Leaf water potential [MPa]
    real(r8) :: veg_tempk  ! Leaf temperature     [K]
    real(r8) :: ceair      ! vapor pressure of air, constrained [Pa]
    real(r8) :: can_press  ! Atmospheric pressure of canopy [Pa]
    real(r8) :: veg_esat   ! Saturated vapor pressure at veg surf [Pa]
    real(r8) :: rb         ! Leaf Boundary layer resistance [s/m]
    real(r8) :: gstoma     ! Stomatal Conductance of this leaf layer [m/s]
    integer  :: ft         ! Plant Functional Type
    real(r8) :: rstoma_out ! Total Stomatal resistance (stoma and BL) [s/m]

    ! Locals
    real(r8) :: k_lwp      ! Scaling coefficient for the ratio of leaf xylem
    ! water potential to mesophyll water potential
    real(r8) :: qs         ! Specific humidity [g/kg]
    real(r8) :: qsat       ! Saturation specific humidity  [g/kg]
    real(r8) :: qsat_adj   ! Adjusted saturation specific humidity  [g/kg]
    real(r8) :: lwp_star   ! leaf water potential scaling coefficient
    ! for inner leaf humidity, 0 means total dehydroted
    ! leaf, 1 means total saturated leaf

    ! Note: to disable this control, set k_lwp to zero, LWP_star will be 1
    k_lwp = EDPftvarcon_inst%hydr_k_lwp(ft)
    if (leaf_psi<0._r8) then
       lwp_star = exp(k_lwp*leaf_psi*molar_mass_water/(rgas_J_K_mol *veg_tempk))
    else 
       lwp_star = 1._r8
    end if

    ! compute specific humidity from vapor pressure
    ! q = molar_mass_ratio_vapdry*e/(can_press - (1-molar_mass_ratio_vapdry)*e) 
    ! source https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    ! now adjust inner leaf humidity by LWP_star

    qs = molar_mass_ratio_vapdry * ceair / (can_press - (1._r8-molar_mass_ratio_vapdry) * ceair)
    qsat = molar_mass_ratio_vapdry * veg_esat / (can_press - (1._r8-molar_mass_ratio_vapdry) * veg_esat)
    qsat_adj = qsat*lwp_star

    ! Adjusting gs (compute a virtual gs) that will be passed to host model

    if ( qsat_adj < qs ) then

       ! if inner leaf vapor pressure is less then or equal to that at leaf surface
       ! then set stomata resistance to be very large to stop the transpiration or back flow of vapor
       rstoma_out = rsmax0

    else

       rstoma_out = (qsat-qs)*( 1/gstoma + rb)/(qsat_adj - qs)-rb

    end if

    if (rstoma_out < nearzero ) then
       write (fates_log(),*) 'qsat:', qsat, 'qs:', qs
       write (fates_log(),*) 'LWP :', leaf_psi
       write (fates_log(),*) 'ceair:', ceair, 'veg_esat:', veg_esat            
       write (fates_log(),*) 'rstoma_out:', rstoma_out, 'rb:', rb  
       write (fates_log(),*) 'LWP_star', lwp_star 
       call endrun(msg=errMsg(sourcefile, __LINE__))                  
    end if

  end function LeafHumidityStomaResis


  ! =====================================================================================

  subroutine ScaleLeafLayerFluxToCohort(nv,          & ! in   currentCohort%nv
       psn_llz,     & ! in   %psn_z(1:currentCohort%nv,ft,cl)
       lmr_llz,     & ! in   lmr_z(1:currentCohort%nv,ft,cl)
       rs_llz,      & ! in   rs_z(1:currentCohort%nv,ft,cl)
       elai_llz,    & ! in   %elai_profile(cl,ft,1:currentCohort%nv)
       c13disc_llz, & ! in   c13disc_z(cl, ft, 1:currentCohort%nv)
       c_area,      & ! in   currentCohort%c_area
       nplant,      & ! in   currentCohort%n
       rb,          & ! in   bc_in(s)%rb_pa(ifp)
       maintresp_reduction_factor, & ! in
       g_sb_laweight, & ! out  currentCohort%g_sb_laweight [m/s] [m2-leaf]
       gpp,         &   ! out  currentCohort%gpp_tstep
       rdark,       &   ! out  currentCohort%rdark
       c13disc_clm, &   ! out currentCohort%c13disc_clm
       cohort_eleaf_area ) ! out [m2]

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

  function ft1_f(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft1_f

  ! =====================================================================================

  function fth_f(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temp function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temp function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temp inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth_f

  ! =====================================================================================

  function fth25_f(hd,se)result(ans)
    !
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES

    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temp function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25_f

  ! =====================================================================================

  subroutine UpdateCanopyNCanNRadPresent(currentPatch)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates two patch level quanities:
    ! currentPatch%nleaf   and
    ! currentPatch%canopy_mask
    !
    ! currentPatch%nleaf(:,:) is a two dimensional array that indicates
    ! the total number of leaf layers (including those that are not exposed to light)
    ! in each canopy layer and for each functional type.
    !
    ! currentPatch%nrad(:,:) is a two dimensional array that indicates
    ! the total number of EXPOSED leaf layers, but for all intents and purposes
    ! in the photosynthesis routine, this appears to be the same as %nleaf...
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
    do cl = 1,currentPatch%ncl_p
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

  ! ====================================================================================

  subroutine GetCanopyGasParameters(can_press, &
       can_o2_partialpress, &
       veg_tempk, &
       air_tempk, &
       air_vpress, &
       veg_esat,   &
       rb,        &
       mm_kco2,   &
       mm_ko2,    &
       co2_cpoint, &
       cf,         &
       gb_mol, &
       ceair)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the specific Michaelis Menten Parameters (pa) for CO2
    ! and O2, as well as the CO2 compentation point.
    ! ---------------------------------------------------------------------------------

    use FatesConstantsMod, only: umol_per_mol
    use FatesConstantsMod, only: mmol_per_mol
    use FatesConstantsMod, only: umol_per_kmol
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    ! Arguments
    real(r8), intent(in) :: can_press           ! Air pressure within the canopy (Pa)
    real(r8), intent(in) :: can_o2_partialpress ! Partial press of o2 in the canopy (Pa)
    real(r8), intent(in) :: veg_tempk           ! The temperature of the vegetation (K)
    real(r8), intent(in) :: air_tempk           ! Temperature of canopy air (K)
    real(r8), intent(in) :: air_vpress          ! Vapor pressure of canopy air (Pa)
    real(r8), intent(in) :: veg_esat            ! Saturated vapor pressure at veg surf (Pa)
    real(r8), intent(in) :: rb                  ! Leaf Boundary layer resistance (s/m)

    real(r8), intent(out) :: mm_kco2       ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(out) :: mm_ko2        !  Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(out) :: co2_cpoint    !  CO2 compensation point (Pa)
    real(r8), intent(out) :: cf            ! conversion factor between molar form and velocity form
    ! of conductance and resistance: [umol/m3]
    real(r8), intent(out) :: gb_mol        ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(out) :: ceair         ! vapor pressure of air, constrained (Pa)

    ! Locals
    real(r8) :: kc25                ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25                ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: sco                 ! relative specificity of rubisco
    real(r8) :: cp25                ! CO2 compensation point at 25C (Pa)

    ! ---------------------------------------------------------------------------------
    ! Intensive values (per mol of air)
    ! kc, ko, currentPatch, from: Bernacchi et al (2001)
    ! Plant, Cell and Environment 24:253-259
    ! ---------------------------------------------------------------------------------

    real(r8), parameter :: mm_kc25_umol_per_mol       = 404.9_r8
    real(r8), parameter :: mm_ko25_mmol_per_mol       = 278.4_r8
    real(r8), parameter :: co2_cpoint_umol_per_mol    = 42.75_r8

    ! Activation energy, from:
    ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
    ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430

    real(r8), parameter :: kcha    = 79430._r8  ! activation energy for kc (J/mol)
    real(r8), parameter :: koha    = 36380._r8  ! activation energy for ko (J/mol)
    real(r8), parameter :: cpha    = 37830._r8  ! activation energy for cp (J/mol)


    ! Derive sco from currentPatch and O2 using present-day O2 (0.209 mol/mol) and re-calculate
    ! currentPatch to account for variation in O2 using currentPatch = 0.5 O2 / sco

    ! FIXME (RGK 11-30-2016 THere are more constants here, but I don't have enough information
    ! about what they are or do, so I can't give them more descriptive names. Someone please
    ! fill this in when possible)

    kc25 = ( mm_kc25_umol_per_mol / umol_per_mol ) * can_press
    ko25 = ( mm_ko25_mmol_per_mol / mmol_per_mol ) * can_press
    sco  = 0.5_r8 * 0.209_r8 / (co2_cpoint_umol_per_mol / umol_per_mol )
    cp25 = 0.5_r8 * can_o2_partialpress / sco

    if( veg_tempk.gt.150_r8 .and. veg_tempk.lt.350_r8 )then
       mm_kco2       = kc25 * ft1_f(veg_tempk, kcha)
       mm_ko2         = ko25 * ft1_f(veg_tempk, koha)
       co2_cpoint     = cp25 * ft1_f(veg_tempk, cpha)
    else
       mm_kco2    = 1.0_r8
       mm_ko2     = 1.0_r8
       co2_cpoint = 1.0_r8
    end if

    ! ---------------------------------------------------------------------------------
    !
    ! cf is the conversion factor between molar form and velocity form
    ! of conductance and resistance: [umol/m3]
    !
    ! i.e.
    ! [m/s] * [umol/m3] -> [umol/m2/s]
    !
    ! Breakdown of the conversion factor: [ umol / m3 ]
    !
    ! Rgas [J /K /kmol]
    ! Air Potential Temperature [ K ]
    ! Canopy Pressure      [ Pa ]
    ! conversion: umol/kmol =  1e9
    !
    ! [ Pa * K * kmol umol/kmol  /  J K ] = [ Pa * umol / J ]
    ! since: 1 Pa = 1 N / m2
    ! [ Pa * umol / J ] = [ N * umol / J m2 ]
    ! since: 1 J = 1 N * m
    ! [ N * umol / J m2 ] = [ N * umol / N m3 ]
    ! [ umol / m3 ]
    !
    ! --------------------------------------------------------------------------------

    cf = can_press/(rgas * air_tempk )*umol_per_kmol
    gb_mol = (1._r8/ rb) * cf

    ! Constrain eair >= 0.05*esat_tv so that solution does not blow up. This ensures
    ! that hs does not go to zero. Also eair <= veg_esat so that hs <= 1
    ceair = min( max(air_vpress, 0.05_r8*veg_esat ),veg_esat )



    return
  end subroutine GetCanopyGasParameters

  ! ====================================================================================

  subroutine LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, &
       nscaler,       &
       ft,            &
       veg_tempk,     &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use EDPftvarcon      , only : EDPftvarcon_inst

    ! -----------------------------------------------------------------------
    ! Base maintenance respiration rate for plant tissues maintresp_leaf_ryan1991_baserate
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    ! Which is the default value of maintresp_nonleaf_baserate

    ! Arguments
    real(r8), intent(in)  :: lnc_top      ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8), intent(in)  :: nscaler      ! Scale for leaf nitrogen profile
    integer,  intent(in)  :: ft           ! (plant) Functional Type Index
    real(r8), intent(in)  :: veg_tempk    ! vegetation temperature
    real(r8), intent(out) :: lmr          ! Leaf Maintenance Respiration  (umol CO2/m**2/s)

    ! Locals
    real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top  ! canopy top leaf maint resp rate at 25C for this pft (umol CO2/m**2/s)
    integer :: c3c4_path_index    ! Index for which photosynthetic pathway

    ! Parameter
    real(r8), parameter :: lmrha = 46390._r8    ! activation energy for lmr (J/mol)
    real(r8), parameter :: lmrhd = 150650._r8   ! deactivation energy for lmr (J/mol)
    real(r8), parameter :: lmrse = 490._r8      ! entropy term for lmr (J/mol/K)
    real(r8), parameter :: lmrc = 1.15912391_r8 ! scaling factor for high
    ! temperature inhibition (25 C = 1.0)

    lmr25top = EDPftvarcon_inst%maintresp_leaf_ryan1991_baserate(ft) * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
    lmr25top = lmr25top * lnc_top / (umolC_to_kgC * g_per_kg)


    ! Part I: Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
    ! ----------------------------------------------------------------------------------
    lmr25 = lmr25top * nscaler

    ! photosynthetic pathway: 0. = c4, 1. = c3
    c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ft))

    if (c3c4_path_index == c3_path_index) then
       ! temperature sensitivity of C3 plants
       lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
            fth_f(veg_tempk, lmrhd, lmrse, lmrc)
    else
       ! temperature sensitivity of C4 plants
       lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
       lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
    endif

    ! Any hydrodynamic limitations could go here, currently none
    ! lmr = lmr * (nothing)

  end subroutine LeafLayerMaintenanceRespiration_Ryan_1991

  ! ====================================================================================   

  subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top, &
       cumulative_lai, &
       vcmax25top,     &
       ft,             &
       veg_tempk,      &
       tgrowth,        &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use FatesConstantsMod, only : lmr_b
    use FatesConstantsMod, only : lmr_c
    use FatesConstantsMod, only : lmr_TrefC
    use FatesConstantsMod, only : lmr_r_1
    use FatesConstantsMod, only : lmr_r_2
    use EDPftvarcon      , only : EDPftvarcon_inst

    ! Arguments
    real(r8), intent(in)  :: lnc_top          ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    integer,  intent(in)  :: ft               ! (plant) Functional Type Index
    real(r8), intent(in)  :: vcmax25top       ! top of canopy vcmax
    real(r8), intent(in)  :: cumulative_lai   ! cumulative lai above the current leaf layer
    real(r8), intent(in)  :: veg_tempk        ! vegetation temperature  (degrees K)
    real(r8), intent(in)  :: tgrowth          ! lagged vegetation temperature averaged over acclimation timescale (degrees K)
    real(r8), intent(out) :: lmr              ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
    
    ! Locals
    real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: r_0     ! base respiration rate, PFT-dependent (umol CO2/m**2/s)
    real(r8) :: r_t_ref ! acclimated ref respiration rate (umol CO2/m**2/s)
    real(r8) :: lmr25top  ! canopy top leaf maint resp rate at 25C for this pft (umol CO2/m**2/s)

    real(r8) :: rdark_scaler ! negative exponential scaling of rdark
    real(r8) :: kn           ! decay coefficient
   
    ! parameter values of r_0 as listed in Atkin et al 2017: (umol CO2/m**2/s) 
    ! Broad-leaved trees  1.7560
    ! Needle-leaf trees   1.4995
    ! Shrubs              2.0749
    ! C3 herbs/grasses    2.1956
    ! In the absence of better information, we use the same value for C4 grasses as C3 grasses.

    ! r_0 currently put into the EDPftvarcon_inst%dev_arbitrary_pft
    ! all figs in Atkin et al 2017 stop at zero Celsius so we will assume acclimation is fixed below that
    r_0 = EDPftvarcon_inst%maintresp_leaf_atkin2017_baserate(ft)

    ! This code uses the relationship between leaf N and respiration from Atkin et al 
    ! for the top of the canopy, but then scales through the canopy based on a rdark_scaler.
    ! To assume proportionality with N through the canopy following Lloyd et al. 2010, use the
    ! default parameter value of 2.43, which results in the scaling of photosynthesis and respiration
    ! being proportional through the canopy. To have a steeper decrease in respiration than photosynthesis
    ! this number can be smaller. There is some observational evidence for this being the case
    ! in Lamour et al. 2023. 

    kn = decay_coeff_vcmax(vcmax25top, &
                           EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(ft), &
                           EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(ft))

    rdark_scaler = exp(-kn * cumulative_lai)
    
    r_t_ref = max(0._r8, rdark_scaler * (r_0 + lmr_r_1 * lnc_top + lmr_r_2 * max(0._r8, (tgrowth - tfrz) )) )

    if (r_t_ref .eq. 0._r8) then
       warn_msg = 'Rdark is negative at this temperature and is capped at 0. tgrowth (C): '//trim(N2S(tgrowth-tfrz))//' pft: '//trim(I2S(ft))
       call FatesWarn(warn_msg,index=4)            
    end if

    lmr = r_t_ref * exp(lmr_b * (veg_tempk - tfrz - lmr_TrefC) + lmr_c * &
         ((veg_tempk-tfrz)**2 - lmr_TrefC**2))

  end subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017

  ! ====================================================================================

  subroutine LeafLayerBiophysicalRates( parsun_per_la, &
       ft,            &
       vcmax25top_ft, &
       jmax25top_ft, &
       co2_rcurve_islope25top_ft, &
       nscaler,    &
       veg_tempk,      &
       dayl_factor, &
       t_growth,   &
       t_home,     &
       btran, &
       vcmax, &
       jmax, &
       co2_rcurve_islope )

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the localized rates of several key photosynthesis
    ! rates.  By localized, we mean specific to the plant type and leaf layer,
    ! which factors in leaf physiology, as well as environmental effects.
    ! This procedure should be called prior to iterative solvers, and should
    ! have pre-calculated the reference rates for the pfts before this.
    !
    ! The output biophysical rates are:
    ! vcmax: maximum rate of carboxilation,
    ! jmax: maximum electron transport rate,
    ! co2_rcurve_islope: initial slope of CO2 response curve (C4 plants)
    ! ---------------------------------------------------------------------------------

    use EDPftvarcon         , only : EDPftvarcon_inst

    ! Arguments
    ! ------------------------------------------------------------------------------

    real(r8), intent(in) :: parsun_per_la   ! PAR absorbed per sunlit leaves for this layer
    integer,  intent(in) :: ft              ! (plant) Functional Type Index
    real(r8), intent(in) :: nscaler         ! Scale for leaf nitrogen profile
    real(r8), intent(in) :: vcmax25top_ft   ! canopy top maximum rate of carboxylation at 25C
    ! for this pft (umol CO2/m**2/s)
    real(r8), intent(in) :: jmax25top_ft    ! canopy top maximum electron transport rate at 25C
    ! for this pft (umol electrons/m**2/s)
    real(r8), intent(in) :: co2_rcurve_islope25top_ft ! initial slope of CO2 response curve
    ! (C4 plants) at 25C, canopy top, this pft
    real(r8), intent(in) :: veg_tempk           ! vegetation temperature
    real(r8), intent(in) :: dayl_factor         ! daylength scaling factor (0-1)
    real(r8), intent(in) :: t_growth            ! T_growth (short-term running mean temperature) (K)
    real(r8), intent(in) :: t_home              ! T_home (long-term running mean temperature) (K)
    real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1)

    real(r8), intent(out) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(out) :: jmax              ! maximum electron transport rate
    ! (umol electrons/m**2/s)
    real(r8), intent(out) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)

    ! Locals
    ! -------------------------------------------------------------------------------
    real(r8) :: vcmax25             ! leaf layer: maximum rate of carboxylation at 25C
    ! (umol CO2/m**2/s)
    real(r8) :: jmax25              ! leaf layer: maximum electron transport rate at 25C
    ! (umol electrons/m**2/s)
    real(r8) :: co2_rcurve_islope25 ! leaf layer: Initial slope of CO2 response curve
    ! (C4 plants) at 25C
    integer :: c3c4_path_index      ! Index for which photosynthetic pathway
    real(r8) :: dayl_factor_local   ! Local version of daylength factor

    ! Parameters
    ! ---------------------------------------------------------------------------------
    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: t_growth_celsius ! average growing temperature
    real(r8) :: t_home_celsius   ! average home temperature
    real(r8) :: jvr            ! ratio of Jmax25 / Vcmax25
    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)

    select case(photo_tempsens_model)
    case (photosynth_acclim_model_none) !No temperature acclimation
       vcmaxha = EDPftvarcon_inst%vcmaxha(FT)
       jmaxha  = EDPftvarcon_inst%jmaxha(FT)
       vcmaxhd = EDPftvarcon_inst%vcmaxhd(FT)
       jmaxhd  = EDPftvarcon_inst%jmaxhd(FT)
       vcmaxse = EDPftvarcon_inst%vcmaxse(FT)
       jmaxse  = EDPftvarcon_inst%jmaxse(FT)
    case (photosynth_acclim_model_kumarathunge_etal_2019) !Kumarathunge et al. temperature acclimation, Thome=30-year running mean
       t_growth_celsius = t_growth-tfrz
       t_home_celsius = t_home-tfrz
       vcmaxha = (42.6_r8 + (1.14_r8*t_growth_celsius))*1e3_r8 !J/mol
       jmaxha = 40.71_r8*1e3_r8 !J/mol
       vcmaxhd = 200._r8*1e3_r8 !J/mol
       jmaxhd = 200._r8*1e3_r8 !J/mol
       vcmaxse = (645.13_r8 - (0.38_r8*t_growth_celsius))
       jmaxse = 658.77_r8 - (0.84_r8*t_home_celsius) - 0.52_r8*(t_growth_celsius-t_home_celsius)
       jvr = 2.56_r8 - (0.0375_r8*t_home_celsius)-(0.0202_r8*(t_growth_celsius-t_home_celsius))
    case default
       write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    vcmaxc = fth25_f(vcmaxhd, vcmaxse)
    jmaxc  = fth25_f(jmaxhd, jmaxse)

    if(parsun_per_la <= 0._r8) then
       vcmax             = 0._r8
       jmax              = 0._r8
       co2_rcurve_islope = 0._r8
    else                                     ! day time

       ! update the daylength factor local variable if the switch is on
       if ( dayl_switch == itrue ) then
          dayl_factor_local = dayl_factor
       else
          dayl_factor_local = 1.0_r8
       endif

       ! Vcmax25top was already calculated to derive the nscaler function
       vcmax25 = vcmax25top_ft * nscaler * dayl_factor_local
       select case(photo_tempsens_model)
       case (photosynth_acclim_model_none)
          jmax25  = jmax25top_ft * nscaler * dayl_factor_local
       case (photosynth_acclim_model_kumarathunge_etal_2019) 
          jmax25 = vcmax25*jvr
       case default
          write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

       co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler

       ! Adjust for temperature
       ! photosynthetic pathway: 0. = c4, 1. = c3
       c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ft))

       if (c3c4_path_index == c3_path_index) then
          vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
       else
          vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
          vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
          vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
       end if

       jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)

       !q10 response of product limited psn.
       co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
    end if

    ! Adjust for water limitations
    vcmax = vcmax * btran

    return

  end subroutine LeafLayerBiophysicalRates

  subroutine lowstorage_maintresp_reduction(frac, pft, maintresp_reduction_factor)

    ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
    ! of this is that mortality of plants increases when storage is low because they are not able
    ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
    ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1,
    ! and the curvature of the function is determined by a parameter.

    ! Uses
    use EDPftvarcon         , only : EDPftvarcon_inst

    ! Arguments
    ! ------------------------------------------------------------------------------
    real(r8), intent(in) :: frac      ! ratio of storage to target leaf biomass
    integer,  intent(in) :: pft       ! what pft is this cohort?
    real(r8), intent(out) :: maintresp_reduction_factor  ! the factor by which to reduce maintenance respiration

    ! --------------------------------------------------------------------------------
    ! Parameters are at the PFT level:
    ! fates_maintresp_reduction_curvature controls the curvature of this.
    ! If this parameter is zero, then there is no reduction until the plant dies at storage = 0.
    ! If this parameter is one, then there is a linear reduction in respiration below the storage point.
    ! Intermediate values will give some (concave-downwards) curvature.
    !
    ! maintresp_reduction_intercept controls the maximum amount of throttling.
    ! zero means no throttling at any point, so it turns this mechanism off completely and so
    ! allows an entire cohort to die via negative carbon-induced termination mortality.
    ! one means complete throttling, so no maintenance respiration at all, when out of carbon.
    ! ---------------------------------------------------------------------------------

    if( frac .lt. 1._r8 )then
       if ( abs(EDPftvarcon_inst%maintresp_reduction_curvature(pft)-1._r8) > nearzero ) then
          maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
               EDPftvarcon_inst%maintresp_reduction_intercept(pft) * &
               (1._r8 - EDPftvarcon_inst%maintresp_reduction_curvature(pft)**frac) &
               / (1._r8-EDPftvarcon_inst%maintresp_reduction_curvature(pft))
       else  ! avoid nan answer for linear case
          maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
               EDPftvarcon_inst%maintresp_reduction_intercept(pft) * frac
       endif

    else
       maintresp_reduction_factor = 1._r8
    endif


  end subroutine lowstorage_maintresp_reduction

end module FATESPlantRespPhotosynthMod
