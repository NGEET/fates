module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure.
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : ifalse
  use FatesConstantsMod         , only : itrue
  use FatesConstantsMod         , only : fates_unset_int
  use FatesConstantsMod         , only : primaryland
  use FatesConstantsMod         , only : nearzero
  use FatesConstantsMod         , only : rsnbl_math_prec
  use EDTypesMod                , only : min_patch_area_forced
  use FatesConstantsMod         , only : n_landuse_cats
  use FatesConstantsMod         , only : is_crop
  use FatesConstantsMod         , only : fates_unset_r8
  use FatesConstantsMod         , only : nearzero, area_error_4, area_error_3
  use FatesGlobals              , only : endrun => fates_endrun
  use EDParamsMod               , only : nclmax
  use EDParamsMod               , only : regeneration_model
  use FatesGlobals              , only : fates_log
  use FatesInterfaceTypesMod    , only : hlm_is_restart
  use FatesInterfaceTypesMod    , only : hlm_current_tod
  use EDPftvarcon               , only : EDPftvarcon_inst
  use PRTParametersMod          , only : prt_params
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts
  use EDCohortDynamicsMod       , only : InitPRTObject
  use EDTypesMod                , only : set_patchno
  use EDPhysiologyMod           , only : calculate_sp_properties
  use ChecksBalancesMod         , only : SiteMassStock
  use FatesInterfaceTypesMod    , only : hlm_day_of_year
  use FatesRadiationMemMod      , only : num_swb
  use EDTypesMod                , only : ed_site_type
  use FatesPatchMod             , only : fates_patch_type
  use FatesCohortMod            , only : fates_cohort_type
  use EDTypesMod                , only : numWaterMem
  use EDTypesMod                , only : num_vegtemp_mem
  use EDTypesMod                , only : area, area_inv
  use EDTypesMod                , only : init_spread_near_bare_ground
  use EDTypesMod                , only : init_spread_inventory
  use FatesConstantsMod         , only : leaves_on
  use FatesConstantsMod         , only : leaves_off
  use FatesConstantsMod         , only : ihard_stress_decid
  use FatesConstantsMod         , only : isemi_stress_decid
  use PRTGenericMod             , only : num_elements
  use PRTGenericMod             , only : element_list
  use EDTypesMod                , only : phen_cstat_nevercold
  use EDTypesMod                , only : phen_cstat_iscold
  use EDTypesMod                , only : phen_dstat_timeoff
  use EDTypesMod                , only : phen_dstat_moistoff
  use EDTypesMod                , only : phen_cstat_notcold
  use EDTypesMod                , only : phen_dstat_moiston
  use FatesInterfaceTypesMod         , only : bc_in_type,bc_out_type
  use FatesInterfaceTypesMod         , only : hlm_use_planthydro
  use FatesInterfaceTypesMod         , only : hlm_use_inventory_init
  use FatesInterfaceTypesMod         , only : hlm_use_fixed_biogeog
  use FatesInterfaceTypesMod         , only : hlm_use_tree_damage
  use FatesInterfaceTypesMod         , only : hlm_use_sp
  use FatesInterfaceTypesMod         , only : hlm_use_luh
  use FatesInterfaceTypesMod         , only : numpft
  use FatesInterfaceTypesMod         , only : nleafage
  use FatesInterfaceTypesMod         , only : nlevsclass
  use FatesInterfaceTypesMod         , only : nlevcoage
  use FatesInterfaceTypesMod         , only : nlevdamage
  use FatesInterfaceTypesMod         , only : hlm_use_nocomp
  use FatesInterfaceTypesMod         , only : nlevage
  use FatesAllometryMod         , only : h2d_allom
  use FatesAllometryMod         , only : h_allom
  use FatesAllometryMod         , only : bagw_allom
  use FatesAllometryMod         , only : bbgw_allom
  use FatesAllometryMod         , only : bleaf
  use FatesAllometryMod         , only : bfineroot
  use FatesAllometryMod         , only : bsap_allom
  use FatesAllometryMod         , only : bdead_allom
  use FatesAllometryMod         , only : bstore_allom
  use FatesAllometryMod         , only : carea_allom
  use PRTGenericMod             , only : StorageNutrientTarget
  use FatesInterfaceTypesMod,      only : hlm_parteh_mode
  use PRTGenericMod,          only : prt_carbon_allom_hyp
  use PRTGenericMod,          only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,          only : prt_vartypes
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : nitrogen_element
  use PRTGenericMod,          only : phosphorus_element
  use PRTGenericMod,          only : SetState
  use FatesSizeAgeTypeIndicesMod,only : get_age_class_index
  use DamageMainMod,          only : undamaged_class
  use FatesConstantsMod,      only : n_term_mort_types
  use FatesInterfaceTypesMod, only : hlm_num_luh2_transitions
  use FatesConstantsMod,      only : nocomp_bareground_land, nocomp_bareground
  use FatesConstantsMod,      only : min_nocomp_pftfrac_perlanduse
  use EdTypesMod,             only : dump_site
  use SFNesterovMod,          only : nesterov_index


  ! CIME GLOBALS
  use shr_log_mod               , only : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : isnan => shr_infnan_isnan

  implicit none
  private

  logical   ::  debug = .false.
  integer :: istat           ! return status code
  character(len=255) :: smsg ! Message string for deallocation errors
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  public  :: zero_site
  public  :: init_site_vars
  public  :: init_patches
  public  :: set_site_properties
  private :: init_cohorts


  ! ============================================================================

contains

  ! ============================================================================

  subroutine init_site_vars( site_in, bc_in, bc_out )
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout) :: site_in
    type(bc_in_type),intent(in)       :: bc_in
    type(bc_out_type),intent(in)      :: bc_out
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
    integer :: el

    !
    allocate(site_in%term_nindivs_canopy(1:n_term_mort_types,1:nlevsclass,1:numpft))
    allocate(site_in%term_nindivs_ustory(1:n_term_mort_types,1:nlevsclass,1:numpft))
    allocate(site_in%demotion_rate(1:nlevsclass))
    allocate(site_in%promotion_rate(1:nlevsclass))
    allocate(site_in%imort_rate(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_canopy(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_ustory(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_cambial(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_crown(1:nlevsclass,1:numpft))
    allocate(site_in%growthflux_fusion(1:nlevsclass,1:numpft))
    allocate(site_in%mass_balance(1:num_elements))
    allocate(site_in%iflux_balance(1:num_elements))

    if (hlm_use_tree_damage .eq. itrue) then 
       allocate(site_in%term_nindivs_canopy_damage(1:nlevdamage, 1:nlevsclass, 1:numpft))
       allocate(site_in%term_nindivs_ustory_damage(1:nlevdamage, 1:nlevsclass, 1:numpft))
       allocate(site_in%imort_rate_damage(1:nlevdamage, 1:nlevsclass, 1:numpft))
       allocate(site_in%imort_cflux_damage(1:nlevdamage, 1:nlevsclass))
       allocate(site_in%term_cflux_canopy_damage(1:nlevdamage, 1:nlevsclass))
       allocate(site_in%term_cflux_ustory_damage(1:nlevdamage, 1:nlevsclass))
       allocate(site_in%fmort_rate_canopy_damage(1:nlevdamage, 1:nlevsclass, 1:numpft))
       allocate(site_in%fmort_rate_ustory_damage(1:nlevdamage, 1:nlevsclass, 1:numpft)) 
       allocate(site_in%fmort_cflux_canopy_damage(1:nlevdamage, 1:nlevsclass))
       allocate(site_in%fmort_cflux_ustory_damage(1:nlevdamage, 1:nlevsclass)) 
    else
       allocate(site_in%term_nindivs_canopy_damage(1,1,1))
       allocate(site_in%term_nindivs_ustory_damage(1,1,1))
       allocate(site_in%imort_rate_damage(1,1,1))
       allocate(site_in%imort_cflux_damage(1,1))
       allocate(site_in%term_cflux_canopy_damage(1,1))
       allocate(site_in%term_cflux_ustory_damage(1,1))
       allocate(site_in%fmort_rate_canopy_damage(1,1,1))
       allocate(site_in%fmort_rate_ustory_damage(1,1,1))
       allocate(site_in%fmort_cflux_canopy_damage(1,1))
       allocate(site_in%fmort_cflux_ustory_damage(1,1))
    end if

    allocate(site_in%term_carbonflux_canopy(1:n_term_mort_types,1:numpft))
    allocate(site_in%term_carbonflux_ustory(1:n_term_mort_types,1:numpft))
    allocate(site_in%imort_carbonflux(1:numpft))
    allocate(site_in%fmort_carbonflux_canopy(1:numpft))
    allocate(site_in%fmort_carbonflux_ustory(1:numpft))

    allocate(site_in%term_abg_flux(1:nlevsclass,1:numpft))
    allocate(site_in%imort_abg_flux(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_abg_flux(1:nlevsclass,1:numpft))

    site_in%nlevsoil   = bc_in%nlevsoil
    allocate(site_in%rootfrac_scr(site_in%nlevsoil))
    allocate(site_in%zi_soil(0:site_in%nlevsoil))
    allocate(site_in%dz_soil(site_in%nlevsoil))
    allocate(site_in%z_soil(site_in%nlevsoil))

    allocate(site_in%area_pft(1:numpft,1:n_landuse_cats))  
    allocate(site_in%landuse_vector_gt_min(1:n_landuse_cats))

    allocate(site_in%use_this_pft(1:numpft))
    allocate(site_in%area_by_age(1:nlevage))

    ! for CNP dynamics, track the mean l2fr of recruits
    ! for different pfts and canopy positions
    allocate(site_in%rec_l2fr(1:numpft,nclmax))


    ! SP mode
    allocate(site_in%sp_tlai(1:numpft))
    allocate(site_in%sp_tsai(1:numpft))
    allocate(site_in%sp_htop(1:numpft))

    ! Allocate site-level flux diagnostics
    ! -----------------------------------------------------------------------
    allocate(site_in%flux_diags%elem(1:num_elements))
    do el=1,num_elements
       allocate(site_in%flux_diags%elem(el)%surf_fine_litter_input(1:numpft))
       allocate(site_in%flux_diags%elem(el)%root_litter_input(1:numpft))
    end do
    allocate(site_in%flux_diags%nh4_uptake_scpf(numpft*nlevsclass))
    allocate(site_in%flux_diags%no3_uptake_scpf(numpft*nlevsclass))
    allocate(site_in%flux_diags%sym_nfix_scpf(numpft*nlevsclass))
    allocate(site_in%flux_diags%n_efflux_scpf(numpft*nlevsclass))
    allocate(site_in%flux_diags%p_uptake_scpf(numpft*nlevsclass))
    allocate(site_in%flux_diags%p_efflux_scpf(numpft*nlevsclass))
    
    

    ! Initialize the static soil
    ! arrays from the boundary (initial) condition
    site_in%zi_soil(:) = bc_in%zi_sisl(:)
    site_in%dz_soil(:) = bc_in%dz_sisl(:)
    site_in%z_soil(:)  = bc_in%z_sisl(:)

    ! Seed dispersal
    allocate(site_in%seed_in(1:numpft))
    allocate(site_in%seed_out(1:numpft))

    allocate(nesterov_index :: site_in%fireWeather)
    call site_in%fireWeather%Init()

  end subroutine init_site_vars

  ! ============================================================================
  subroutine zero_site( site_in )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout) ::  site_in
    !
    ! !LOCAL VARIABLES:
    integer :: el
    !----------------------------------------------------------------------

    site_in%oldest_patch     => null() ! pointer to oldest patch at the site
    site_in%youngest_patch   => null() ! pointer to yngest patch at the site


    ! PHENOLOGY

    site_in%cstatus          = fates_unset_int    ! are leaves in this pixel on or off?
    site_in%dstatus(:)       = fates_unset_int
    site_in%grow_deg_days    = nan  ! growing degree days
    site_in%snow_depth       = nan
    site_in%nchilldays       = fates_unset_int
    site_in%ncolddays        = fates_unset_int
    site_in%cleafondate      = fates_unset_int  ! doy of leaf on (cold)
    site_in%cleafoffdate     = fates_unset_int  ! doy of leaf off (cold)
    site_in%dleafondate(:)   = fates_unset_int  ! doy of leaf on (drought)
    site_in%dleafoffdate(:)  = fates_unset_int  ! doy of leaf off (drought)
    site_in%cndaysleafon     = fates_unset_int  ! days since leaf on (cold)
    site_in%cndaysleafoff    = fates_unset_int  ! days since leaf off (cold)
    site_in%dndaysleafon(:)  = fates_unset_int  ! days since leaf on (drought)
    site_in%dndaysleafoff(:) = fates_unset_int  ! days since leaf off (drought)
    site_in%elong_factor(:)  = nan              ! Elongation factor (0 - full abscission; 1 - fully flushed)

    site_in%liqvol_memory(:,:)  = nan
    site_in%smp_memory(:,:)  = nan
    site_in%vegtemp_memory(:) = nan              ! record of last 10 days temperature for senescence model.

    site_in%phen_model_date  = fates_unset_int

    ! Disturbance rates tracking
    site_in%primary_land_patchfusion_error = 0.0_r8
    site_in%disturbance_rates(:,:,:) = 0.0_r8
    site_in%landuse_transition_matrix(:,:) = 0.0_r8

    ! FIRE
    site_in%FDI              = 0.0_r8     ! daily fire danger index (0-1)
    site_in%NF               = 0.0_r8     ! daily lightning strikes per km2
    site_in%NF_successful    = 0.0_r8     ! daily successful iginitions per km2

    do el=1,num_elements
       ! Zero the state variables used for checking mass conservation
       call site_in%mass_balance(el)%ZeroMassBalState()
       call site_in%mass_balance(el)%ZeroMassBalFlux()
    end do

    call site_in%flux_diags%ZeroFluxDiags()
    
    ! This will be initialized in FatesSoilBGCFluxMod:PrepCH4BCs()
    ! It checks to see if the value is below -9000. If it is,
    ! it will assume the first value of the smoother is set
    site_in%ema_npp = -9999.9_r8

    ! termination and recruitment info
    site_in%term_nindivs_canopy(:,:,:) = 0._r8
    site_in%term_nindivs_ustory(:,:,:) = 0._r8
    site_in%term_crownarea_canopy = 0._r8
    site_in%term_crownarea_ustory = 0._r8
    site_in%imort_crownarea = 0._r8
    site_in%fmort_crownarea_canopy = 0._r8
    site_in%fmort_crownarea_ustory = 0._r8
    site_in%term_carbonflux_canopy(:,:) = 0._r8
    site_in%term_carbonflux_ustory(:,:) = 0._r8
    site_in%recruitment_rate(:) = 0._r8
    site_in%imort_rate(:,:) = 0._r8
    site_in%imort_carbonflux(:) = 0._r8
    site_in%fmort_rate_canopy(:,:) = 0._r8
    site_in%fmort_rate_ustory(:,:) = 0._r8
    site_in%fmort_carbonflux_canopy(:) = 0._r8
    site_in%fmort_carbonflux_ustory(:) = 0._r8
    site_in%fmort_rate_cambial(:,:) = 0._r8
    site_in%fmort_rate_crown(:,:) = 0._r8
    site_in%term_abg_flux(:,:) = 0._r8
    site_in%imort_abg_flux(:,:) = 0._r8
    site_in%fmort_abg_flux(:,:) = 0._r8

    ! fusoin-induced growth flux of individuals
    site_in%growthflux_fusion(:,:) = 0._r8

    ! demotion/promotion info
    site_in%demotion_rate(:) = 0._r8
    site_in%demotion_carbonflux = 0._r8
    site_in%promotion_rate(:) = 0._r8
    site_in%promotion_carbonflux = 0._r8

    ! damage transition info
    site_in%imort_rate_damage(:,:,:) = 0._r8
    site_in%term_nindivs_canopy_damage(:,:,:) = 0._r8
    site_in%term_nindivs_ustory_damage(:,:,:) = 0._r8
    site_in%imort_cflux_damage(:,:) = 0._r8
    site_in%term_cflux_canopy_damage(:,:) = 0._r8
    site_in%term_cflux_ustory_damage(:,:) = 0._r8
    site_in%crownarea_canopy_damage = 0._r8
    site_in%crownarea_ustory_damage = 0._r8
    site_in%fmort_rate_canopy_damage(:,:,:) = 0._r8
    site_in%fmort_rate_ustory_damage(:,:,:) = 0._r8
    site_in%fmort_cflux_canopy_damage(:,:) = 0._r8
    site_in%fmort_cflux_ustory_damage(:,:) = 0._r8

    ! Resources management (logging/harvesting, etc)
    site_in%resources_management%harvest_debt = 0.0_r8
    site_in%resources_management%harvest_debt_sec = 0.0_r8
    site_in%resources_management%trunk_product_site  = 0.0_r8

    ! canopy spread
    site_in%spread = 0._r8


    site_in%area_pft(:,:) = 0._r8
    site_in%area_bareground = 0._r8

    ! Seed dispersal
    site_in%seed_in(:) = 0.0_r8
    site_in%seed_out(:) = 0.0_r8

    site_in%use_this_pft(:) = fates_unset_int
    site_in%area_by_age(:) = 0._r8

    site_in%transition_landuse_from_off_to_on = .false.

  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties( nsites, sites,bc_in )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use EDParamsMod, only : crop_lu_pft_vector
    use EDParamsMod, only : max_nocomp_pfts_by_landuse
    !
    ! !ARGUMENTS

    integer, intent(in)                :: nsites
    type(ed_site_type) , intent(inout) :: sites(nsites)
    type(bc_in_type), intent(in)       :: bc_in(nsites)
    !
    ! !LOCAL VARIABLES:
    integer  :: s
    integer  :: cstat      ! cold status phenology flag
    real(r8) :: GDD
    integer  :: dstat      ! drought status phenology flag
    real(r8) :: liqvolmem
    real(r8) :: smpmem
    real(r8) :: elong_factor ! Elongation factor (0 - fully off; 1 - fully on)
    integer  :: cleafon    ! DOY for cold-decid leaf-on, initial guess
    integer  :: cleafoff   ! DOY for cold-decid leaf-off, initial guess
    integer  :: dleafoff   ! DOY for drought-decid leaf-off, initial guess
    integer  :: dleafon    ! DOY for drought-decid leaf-on, initial guess
    integer  :: cndleafon  ! days since leaf on  (cold), initial guess
    integer  :: cndleafoff ! days since leaf off  (cold), initial guess
    integer  :: dndleafon  ! days since leaf on  (drought), initial guess
    integer  :: dndleafoff ! days since leaf off (drought), initial guess
    integer  :: ft         ! PFT loop
    real(r8) :: sumarea    ! area of PFTs in nocomp mode.
    integer  :: hlm_pft    ! used in fixed biogeog mode
    integer  :: fates_pft  ! used in fixed biogeog mode
    integer  :: i_landusetype
    real(r8) :: temp_vec(numpft)  ! temporary vector
    integer  :: i_pftcount
    !----------------------------------------------------------------------


    ! If this is not a restart, we need to start with some reasonable
    ! starting points. If this is a restart, we leave the values
    ! as unset ints and reals, and let the restart values be read in
    ! after this routine

    if ( hlm_is_restart == ifalse ) then

       GDD      = 30.0_r8
       cleafon  = 100
       cleafoff = 300
       cndleafon  = 0
       cndleafoff = 0
       cstat    = phen_cstat_notcold     ! Leaves are on
       dstat    = phen_dstat_moiston     ! Leaves are on
       dleafoff = 300
       dleafon  = 100
       dndleafon  = 0
       dndleafoff = 0
       liqvolmem  = 0.5_r8 
       smpmem     = 0._r8
       elong_factor = 1._r8

       do s = 1,nsites
          sites(s)%nchilldays    = 0
          sites(s)%ncolddays     = 0        ! recalculated in phenology
                                            ! immediately, so yes this
                                            ! is memory-less, but needed
                                            ! for first value in history file
          sites(s)%phen_model_date         = 0
          sites(s)%cleafondate             = cleafon    - hlm_day_of_year
          sites(s)%cleafoffdate            = cleafoff   - hlm_day_of_year
          sites(s)%cndaysleafon            = cndleafon
          sites(s)%cndaysleafoff           = cndleafoff
          sites(s)%dleafoffdate (1:numpft) = dleafoff   - hlm_day_of_year
          sites(s)%dleafondate  (1:numpft) = dleafon    - hlm_day_of_year
          sites(s)%dndaysleafon (1:numpft) = dndleafon
          sites(s)%dndaysleafoff(1:numpft) = dndleafoff
          sites(s)%grow_deg_days   = GDD

          sites(s)%liqvol_memory(1:numWaterMem,1:numpft) = liqvolmem
          sites(s)%smp_memory(1:numWaterMem,1:numpft) = smpmem
          sites(s)%vegtemp_memory(1:num_vegtemp_mem) = 0._r8

          sites(s)%cstatus = cstat
          sites(s)%dstatus(1:numpft) = dstat
          sites(s)%elong_factor(1:numpft) = elong_factor

          sites(s)%NF         = 0.0_r8
          sites(s)%NF_successful  = 0.0_r8
          sites(s)%area_pft(:,:) = 0.0_r8

          do ft =  1,numpft
             sites(s)%rec_l2fr(ft,:) = prt_params%allom_l2fr(ft)
          end do

          ! Its difficult to come up with a resonable starting smoothing value, so
          ! we initialize on a cold-start to -1
          sites(s)%ema_npp = -9999._r8

          if(hlm_use_fixed_biogeog.eq.itrue)then

             use_fates_luh_if: if (hlm_use_luh .eq. itrue) then
                ! MAPPING OF FATES PFTs on to HLM_PFTs with land use
                ! add up the area associated with each FATES PFT
                ! where pft_areafrac_lu is the area of land in each HLM PFT and land use type (from surface dataset)
                ! hlm_pft_map is the area of that land in each FATES PFT (from param file)

                ! First check for NaNs in bc_in(s)%pft_areafrac_lu. If so, make everything bare ground.
                if ( .not. (any( isnan( bc_in(s)%pft_areafrac_lu (:,:) )) .or. isnan( bc_in(s)%baregroundfrac))) then
                   do i_landusetype = 1, n_landuse_cats
                      if (.not. is_crop(i_landusetype)) then
                         do hlm_pft = 1,size( EDPftvarcon_inst%hlm_pft_map,2)
                            do fates_pft = 1,numpft ! loop round all fates pfts for all hlm pfts
                               sites(s)%area_pft(fates_pft,i_landusetype) = sites(s)%area_pft(fates_pft,i_landusetype) + &
                                    EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft) * bc_in(s)%pft_areafrac_lu(hlm_pft,i_landusetype)
                            end do
                         end do !hlm_pft
                      else
                         ! for crops, we need to use different logic because the bc_in(s)%pft_areafrac_lu() information only exists for natural PFTs
                         sites(s)%area_pft(crop_lu_pft_vector(i_landusetype),i_landusetype) = 1._r8
                      endif
                   end do

                   sites(s)%area_bareground = bc_in(s)%baregroundfrac
                else
                   !if ( all( isnan( bc_in(s)%pft_areafrac_lu (:,:))) .and. isnan(bc_in(s)%baregroundfrac)) then
                      ! if given all NaNs, then make everything bare ground
                      sites(s)%area_bareground = 1._r8
                      sites(s)%area_pft(:,:) = 0._r8
                      write(fates_log(),*) 'Nan values for pftareafrac. dumping site info.'
                      call dump_site(sites(s))
                   !else
                   !   ! if only some things are NaN but not all, then something terrible has probably happened. crash.
                   !   write(fates_log(),*) 'some but, not all, of the data in the PFT by LU matrix at this site is NaN.'
                   !   write(fates_log(),*) 'recommend checking the dataset to see what has happened.'
                   !   call endrun(msg=errMsg(sourcefile, __LINE__))
                   !endif
                endif

             else
                ! MAPPING OF FATES PFTs on to HLM_PFTs
                ! add up the area associated with each FATES PFT
                ! where pft_areafrac is the area of land in each HLM PFT and (from surface dataset)
                ! hlm_pft_map is the area of that land in each FATES PFT (from param file)

                do hlm_pft = 1,size( EDPftvarcon_inst%hlm_pft_map,2)
                   do fates_pft = 1,numpft ! loop round all fates pfts for all hlm pfts
                      sites(s)%area_pft(fates_pft,primaryland) = sites(s)%area_pft(fates_pft,primaryland) + &
                           EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft) * bc_in(s)%pft_areafrac(hlm_pft)
                   end do
                   sites(s)%area_bareground = bc_in(s)%pft_areafrac(0)
                end do !hlm_pft

             endif use_fates_luh_if

             ! handle some edge cases
             do i_landusetype = 1, n_landuse_cats
                do ft =  1,numpft

                   ! remove tiny patches to prevent numerical errors in terminate patches
                   if (sites(s)%area_pft(ft, i_landusetype) .lt. min_nocomp_pftfrac_perlanduse &
                        .and. sites(s)%area_pft(ft, i_landusetype) .gt. nearzero) then
                      if(debug) write(fates_log(),*)  'removing small numbers in site%area_pft',s,ft,i_landusetype,sites(s)%area_pft(ft, i_landusetype)
                      sites(s)%area_pft(ft, i_landusetype)=0.0_r8
                   endif

                   ! if any areas are negative, then end run
                   if(sites(s)%area_pft(ft, i_landusetype).lt.0._r8)then
                      write(fates_log(),*) 'negative area',s,ft,i_landusetype,sites(s)%area_pft(ft, i_landusetype)
                      call endrun(msg=errMsg(sourcefile, __LINE__))
                   end if
                end do
             end do

             ! if in nocomp mode, and the number of nocomp PFTs of a given land use type is greater than the maximum number of patches
             ! allowed to be allocated for that land use type, then only keep the number of PFTs correspondign to the number of patches
             ! allowed on that land use type, starting with the PFTs with greatest area coverage and working down
             if (hlm_use_nocomp .eq. itrue) then
                do i_landusetype = 1, n_landuse_cats
                   ! count how many PFTs have areas greater than zero and compare to the number of patches allowed
                   if (COUNT(sites(s)%area_pft(:, i_landusetype) .gt. 0._r8) > max_nocomp_pfts_by_landuse(i_landusetype)) then
                      ! write current vector to log file
                      if(debug) write(fates_log(),*)  'too many PFTs for LU type ', i_landusetype, sites(s)%area_pft(:, i_landusetype)

                      ! start from largest area, put that PFT's area into a temp vector, and then work down to successively smaller-area PFTs,
                      ! at the end replace the original vector with the temp vector
                      temp_vec(:) = 0._r8
                      do i_pftcount = 1, max_nocomp_pfts_by_landuse(i_landusetype)
                         temp_vec(MAXLOC(sites(s)%area_pft(:, i_landusetype))) = &
                              sites(s)%area_pft(MAXLOC(sites(s)%area_pft(:, i_landusetype)), i_landusetype)
                         sites(s)%area_pft(MAXLOC(sites(s)%area_pft(:, i_landusetype)), i_landusetype) = 0._r8
                      end do
                      sites(s)%area_pft(:, i_landusetype) = temp_vec(:)

                      ! write adjusted vector to log file
                      if(debug) write(fates_log(),*)  'new PFT vector for LU type', i_landusetype, sites(s)%area_pft(:, i_landusetype)
                   endif
                end do
             end if

             ! re-normalize PFT area to ensure it sums to one for each (active) land use type
             ! for nocomp cases, track bare ground area as a separate quantity
             do i_landusetype = 1, n_landuse_cats
                sumarea = sum(sites(s)%area_pft(:,i_landusetype))
                if(sumarea.gt.nearzero)then
                      sites(s)%area_pft(:, i_landusetype) = sites(s)%area_pft(:, i_landusetype)/sumarea
                else
                   ! if no PFT area in primary lands, set bare ground fraction to one.
                   if ( i_landusetype .eq. primaryland) then
                      sites(s)%area_bareground = 1._r8
                      sites(s)%area_pft(:, i_landusetype) = 0._r8
                   endif
                end if
             end do
                
          end if !fixed biogeog

          do ft = 1,numpft
             ! Setting this to true ensures that all pfts
             ! are used for nocomp with no biogeog
             sites(s)%use_this_pft(ft) = itrue
             if(hlm_use_fixed_biogeog.eq.itrue)then
                if(any(sites(s)%area_pft(ft,:).gt.0.0_r8))then
                   sites(s)%use_this_pft(ft) = itrue
                else
                   sites(s)%use_this_pft(ft) = ifalse
                end if !area
             end if !SBG
          end do !ft

          ! need to set the minimum amount of allowable land-use fraction on a given site. this is a function of the minimum allowable patch size,
          ! and for nocomp simulations also the bare ground fraction and the minimum pft fraction for a given land-use type.
          if (hlm_use_nocomp .eq. itrue ) then
             if ( (1._r8 - sites(s)%area_bareground) .gt. nearzero) then
                sites(s)%min_allowed_landuse_fraction =   min_patch_area_forced / (AREA * min_nocomp_pftfrac_perlanduse * (1._r8 - sites(s)%area_bareground))
             else
                ! if all bare ground, shouldn't matter. but make it one anyway to really ignore land use (which should all be NaNs anyway)
                sites(s)%min_allowed_landuse_fraction = 1._r8
             endif
          else
             sites(s)%min_allowed_landuse_fraction =  min_patch_area_forced / AREA
          endif

       end do !site loop
    end if !restart

    return
  end subroutine set_site_properties


  ! ============================================================================
  subroutine init_patches( nsites, sites, bc_in)
    !
    ! !DESCRIPTION:
    ! initialize patches
    ! This may be call a near bare ground initialization, or it may
    ! load patches from an inventory.

    use FatesPlantHydraulicsMod, only : updateSizeDepRhizHydProps
    use FatesInventoryInitMod,   only : initialize_sites_by_inventory
    use FatesLandUseChangeMod,   only : GetLUHStatedata

    !
    ! !ARGUMENTS
    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    type(bc_in_type), intent(in)               :: bc_in(nsites)
    !
    ! !LOCAL VARIABLES:
    integer  :: s
    integer  :: el
    real(r8) :: age !notional age of this patch
    integer  :: ageclass
    real(r8) :: area_diff
    real(r8) :: area_error

    ! dummy locals
    real(r8) :: biomass_stock
    real(r8) :: litter_stock
    real(r8) :: seed_stock
    integer  :: n
    integer  :: start_patch
    integer  :: num_nocomp_pfts
    integer  :: nocomp_pft
    real(r8) :: newparea, newparea_withlanduse
    real(r8) :: total !check on area
    real(r8) :: litt_init  !invalid for satphen, 0 otherwise
    real(r8) :: old_carea
    logical  :: is_first_patch
    ! integer  :: n_luh_states
    ! integer  :: luh_state_counter
    real(r8) :: state_vector(n_landuse_cats)  ! [m2/m2]
    integer  :: i_lu, i_lu_state
    integer  :: n_active_landuse_cats
    integer  :: end_landuse_idx

    type(ed_site_type),  pointer :: sitep
    type(fates_patch_type), pointer :: newppft(:)
    type(fates_patch_type), pointer :: newp
    type(fates_cohort_type), pointer :: cohort
    type(fates_patch_type), pointer :: currentPatch

    ! List out some nominal patch values that are used for Near Bear Ground initializations
    ! as well as initializing inventory
    age                  = 0.0_r8
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    ! Two primary options, either a Near Bear Ground (NBG) or Inventory based cold-start
    ! ---------------------------------------------------------------------------------------------

    if ( hlm_use_inventory_init.eq.itrue ) then

       ! Initialize the site-level crown area spread factor (0-1)
       ! It is likely that closed canopy forest inventories
       ! have smaller spread factors than bare ground (they are crowded)
       do s = 1, nsites
          sites(s)%spread     = init_spread_inventory
       enddo

       call initialize_sites_by_inventory(nsites,sites,bc_in)


       ! For carbon balance checks, we need to initialize the
       ! total carbon stock
       do s = 1, nsites
          do el=1,num_elements
             call SiteMassStock(sites(s),el,sites(s)%mass_balance(el)%old_stock, &
                  biomass_stock,litter_stock,seed_stock)
             ! Initialize the integrated flux balance diagnostics
             ! No need to initialize the instantaneous states, those are re-calculated                                                  
             sites(s)%iflux_balance(el)%iflux_liveveg = &
                  (biomass_stock + seed_stock)*area_inv
             sites(s)%iflux_balance(el)%iflux_litter  = litter_stock * area_inv

          end do
          call set_patchno(sites(s),.false.,0)
       enddo
       
    else

       if(hlm_use_nocomp.eq.itrue)then
          num_nocomp_pfts = numpft
       else !default
          num_nocomp_pfts = 1
       end if !nocomp

       sites_loop: do s = 1, nsites
          sites(s)%sp_tlai(:) = 0._r8
          sites(s)%sp_tsai(:) = 0._r8
          sites(s)%sp_htop(:) = 0._r8

          ! Initialize the site-level crown area spread factor (0-1)
          ! It is likely that closed canopy forest inventories
          ! have smaller spread factors than bare ground (they are crowded)
          sites(s)%spread     = init_spread_near_bare_ground

          ! read in luh state data to determine initial land use types
          if (hlm_use_luh .eq. itrue) then

             ! Set the number of active land use categories to the maximum number
             ! This could be updated in the future to allow a variable number of
             ! categories based on which states are zero
             n_active_landuse_cats = n_landuse_cats
             call GetLUHStatedata(bc_in(s), state_vector)

             ! if the land use state vector is greater than the minimum value, set landuse_vector_gt_min flag to true
             ! otherwise set to false.
             do i_lu_state = 1, n_landuse_cats
                if (state_vector(i_lu_state) .gt. sites(s)%min_allowed_landuse_fraction) then
                   sites(s)%landuse_vector_gt_min(i_lu_state) = .true.
                else
                   sites(s)%landuse_vector_gt_min(i_lu_state) = .false.
                end if
             end do

          else
             ! If LUH2 data is not being used, we initialize with primarylands,
             ! i.e. array index equals '1'
             n_active_landuse_cats = primaryland
             state_vector(:) = 0._r8
             state_vector(primaryland) = 1._r8
          endif

          ! confirm that state vector sums to 1.
          if (abs(sum(state_vector(:))-1._r8) .gt. rsnbl_math_prec) then
             write(fates_log(),*) 'error that the state vector must sum to 1, but doesnt'
             write(fates_log(),*) 'sum(state_vector)', sum(state_vector)
             write(fates_log(),*) state_vector
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

          is_first_patch = .true.

          area_error = 0._r8
          ! first make a bare-ground patch if one is needed.
          make_bareground_patch_if: if (hlm_use_nocomp.eq.itrue .and. hlm_use_fixed_biogeog .eq.itrue) then

             newparea = area * sites(s)%area_bareground
             if (newparea  .gt. min_patch_area_forced) then
                
                allocate(newp)

                call newp%Create(age, newparea, nocomp_bareground_land, nocomp_bareground,     &
                     num_swb, numpft, sites(s)%nlevsoil, hlm_current_tod,      &
                     regeneration_model)

                ! set pointers for first patch (or only patch, if nocomp is false)
                newp%patchno = 1
                newp%younger => null()
                newp%older   => null()
                sites(s)%youngest_patch => newp
                sites(s)%oldest_patch   => newp
                is_first_patch = .false.

                ! Initialize the litter pools to zero, these
                ! pools will be populated by looping over the existing patches
                ! and transfering in mass
                if(hlm_use_sp.eq.itrue)then
                   litt_init = fates_unset_r8
                else
                   litt_init = 0._r8
                end if
                do el=1,num_elements
                   call newp%litter(el)%InitConditions(init_leaf_fines=litt_init, &
                        init_root_fines=litt_init, &
                        init_ag_cwd=litt_init, &
                        init_bg_cwd=litt_init, &
                        init_seed=litt_init,   &
                        init_seed_germ=litt_init)
                end do

             else
                area_error = area_error + newparea
             endif
          endif make_bareground_patch_if

          if (hlm_use_luh .eq. itrue) then
             end_landuse_idx = n_landuse_cats
          else
             end_landuse_idx = 1
          endif


          ! Next, create the non-bareground patches. We do this for either of two scenarios:
          ! If 1) we are not doing both nocomp & fixed-biogeo
          !    2) we are, but there is some non-zero bare-ground area
          not_all_bare_if: if( ((1._r8 - sites(s)%area_bareground) > nearzero) .or. &
                                     (.not.(hlm_use_nocomp.eq.itrue .and. hlm_use_fixed_biogeog.eq.itrue)) ) then
          
             ! now make one or more vegetated patches based on nocomp and land use logic
             luh_state_loop: do i_lu_state = 1, end_landuse_idx
                lu_state_present_if: if (state_vector(i_lu_state) .gt. nearzero) then
                   new_patch_nocomp_loop: do n = 1, num_nocomp_pfts
                      ! set the PFT index for patches if in nocomp mode.
                      if(hlm_use_nocomp.eq.itrue)then
                         nocomp_pft = n
                      else
                         nocomp_pft = fates_unset_int
                      end if
                      
                      if(hlm_use_nocomp.eq.itrue)then
                         ! In no competition mode, if we are using the fixed_biogeog filter
                         ! then each PFT has the area dictated  by the surface dataset.

                         ! If we are not using fixed biogeog model, each PFT gets the same area.
                         ! i.e. each grid cell is divided exactly into the number of FATES PFTs.

                         if(hlm_use_fixed_biogeog.eq.itrue)then
                            newparea = sites(s)%area_pft(nocomp_pft,i_lu_state) * area * state_vector(i_lu_state) &
                                 * (1._r8 - sites(s)%area_bareground)
                         else
                            newparea = area * state_vector(i_lu_state) / numpft
                         end if
                      else  ! The default case is initialized w/ one patch with the area of the whole site.
                         newparea = area * state_vector(i_lu_state)
                      end if  !nocomp mode

                      ! Stop patches being initilialized when PFT not present in nocomop mode
                      new_patch_area_gt_zero: if(newparea .gt. min_patch_area_forced) then 
                         allocate(newp)

                         call newp%Create(age, newparea, i_lu_state, nocomp_pft, &
                              num_swb, numpft, sites(s)%nlevsoil, hlm_current_tod, &
                              regeneration_model)

                         if (is_first_patch) then !is this the first patch?
                            ! set pointers for first patch (or only patch, if nocomp is false)
                            newp%patchno = 1
                            newp%younger => null()
                            newp%older   => null()
                            sites(s)%youngest_patch => newp
                            sites(s)%oldest_patch   => newp
                            is_first_patch = .false.
                         else
                            ! Set pointers for N>1 patches. Note this only happens when nocomp mode is on, or land use is on.
                            ! The new patch is the 'youngest' one, arbitrarily.
                            newp%patchno = nocomp_pft + (i_lu_state-1) * numpft
                            newp%older     => sites(s)%youngest_patch
                            newp%younger   => null()
                            sites(s)%youngest_patch%younger => newp
                            sites(s)%youngest_patch   => newp
                         end if

                         ! Initialize the litter pools to zero, these
                         ! pools will be populated by looping over the existing patches
                         ! and transfering in mass
                         if(hlm_use_sp.eq.itrue)then
                            litt_init = fates_unset_r8
                         else
                            litt_init = 0._r8
                         end if
                         do el=1,num_elements
                            call newp%litter(el)%InitConditions(init_leaf_fines=litt_init, &
                                 init_root_fines=litt_init, &
                                 init_ag_cwd=litt_init, &
                                 init_bg_cwd=litt_init, &
                                 init_seed=litt_init,   &
                                 init_seed_germ=litt_init)
                         end do

                         sitep => sites(s)
                         call init_cohorts(sitep, newp, bc_in(s))

                      else
                         area_error = area_error+ newparea
                      end if new_patch_area_gt_zero
                   end do new_patch_nocomp_loop
                end if lu_state_present_if
             end do luh_state_loop
          end if not_all_bare_if

          ! if we had to skip small patches above, resize things accordingly
          if ( area_error .gt. nearzero) then
             newp => sites(s)%oldest_patch
             do while (associated(newp))
                newp%area = newp%area * area/ (area - area_error)
                newp => newp%younger
             end do
          endif
          
          !check if the total area adds to the same as site area
          total = 0.0_r8
          newp => sites(s)%oldest_patch
          do while (associated(newp))
             total = total + newp%area
             newp => newp%younger
          end do
         
          area_diff = total - area
          if (abs(area_diff) > nearzero) then
             if (abs(area_diff) < area_error_4) then ! this is a precision error

                ! adjust areas of all patches so that they add up to total area
                newp => sites(s)%oldest_patch
                do while (associated(newp))
                   newp%area = newp%area * (area / total)
                   newp => newp%younger
                end do

             else !this is a big error not just a precision error.
                write(fates_log(),*) 'issue with patch area in EDinit', area_diff, total,sites(s)%lat,sites(s)%lon
                write(fates_log(),*) 'hlm_use_nocomp: ',hlm_use_nocomp
                write(fates_log(),*) 'hlm_use_fixed_biogeog: ',hlm_use_fixed_biogeog
                newp => sites(s)%oldest_patch
                do while (associated(newp))
                   write(fates_log(),*) newp%area, newp%nocomp_pft_label, newp%land_use_label
                   newp => newp%younger
                end do
                write(fates_log(),*) 'state_vector', state_vector
                write(fates_log(),*) 'area_error', area_error                
                write(fates_log(),*) 'area_bareground', sites(s)%area_bareground
                do i_lu_state = 1, end_landuse_idx
                   write(fates_log(),*) 'sites(s)%area_pft(:,i_lu_state)',i_lu_state, sites(s)%area_pft(:,i_lu_state)
                end do
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if  ! big error
          end if ! too much patch area
	  
          ! we might have messed up crown areas now - need to correct if SP mode
          if (hlm_use_sp .eq. itrue) then 
            newp => sites(s)%oldest_patch
            do while (associated(newp))
              cohort => newp%tallest
              do while (associated(cohort))
                if (abs(cohort%c_area - newp%area) < area_error_3) then ! correct if it's a very small error
                  old_carea = cohort%c_area 
                  cohort%c_area = cohort%c_area - (cohort%c_area - newp%area)
                  cohort%n = cohort%n*(cohort%c_area/old_carea)
                end if 
                cohort => cohort%shorter
              end do 
              newp => newp%younger
            end do
          end if 

          ! For carbon balance checks, we need to initialize the
          ! total carbon stock
          do el=1,num_elements
             call SiteMassStock(sites(s),el,sites(s)%mass_balance(el)%old_stock, &
                  biomass_stock,litter_stock,seed_stock)
             
             ! Initialize the integrated flux balance diagnostics
             ! No need to initialize the instantaneous states, those are re-calculated
             sites(s)%iflux_balance(el)%iflux_liveveg = &
                  (biomass_stock + seed_stock)*area_inv
             sites(s)%iflux_balance(el)%iflux_litter  = litter_stock * area_inv
             
          end do

          call set_patchno(sites(s),.false.,0)

       enddo sites_loop 
    end if

    ! zero all the patch fire variables for the first timestep
    do s = 1, nsites
       currentPatch => sites(s)%youngest_patch
       do while(associated(currentPatch))

          currentPatch%livegrass                  = 0._r8
          currentPatch%ros_front                  = 0._r8
          currentPatch%tau_l                      = 0._r8
          currentPatch%tfc_ros                    = 0._r8
          currentPatch%fi                         = 0._r8
          currentPatch%fire                       = 0
          currentPatch%fd                         = 0._r8
          currentPatch%ros_back                   = 0._r8
          currentPatch%scorch_ht(:)               = 0._r8
          currentPatch%frac_burnt                 = 0._r8
          currentPatch => currentPatch%older
       enddo
    enddo

    ! This sets the rhizosphere shells based on the plant initialization
    ! The initialization of the plant-relevant hydraulics variables
    ! were set from a call inside of the init_cohorts()->create_cohort() subroutine
    if (hlm_use_planthydro.eq.itrue) then
       do s = 1, nsites
          sitep => sites(s)
          call updateSizeDepRhizHydProps(sitep, bc_in(s))
       end do
    end if

    ! check to make sure there are no very tiny patches
    do s = 1, nsites
       currentPatch => sites(s)%youngest_patch
       do while(associated(currentPatch))
          if (currentPatch%area .lt. min_patch_area_forced) then
             write(fates_log(),*) 'edinit somehow making tiny patches',currentPatch%land_use_label, currentPatch%nocomp_pft_label, currentPatch%area 
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          currentPatch => currentPatch%older
       end do
    end do

    return
  end subroutine init_patches

  ! ============================================================================
 
   subroutine init_cohorts(site_in, patch_in, bc_in)
      !
      ! DESCRIPTION:
      ! initialize new cohorts on bare ground
      !

      ! ARGUMENTS
      type(ed_site_type),     intent(inout),  pointer  :: site_in
      type(fates_patch_type), intent(inout),  pointer  :: patch_in
      type(bc_in_type),       intent(in)               :: bc_in

      ! LOCAL VARIABLES:
      class(prt_vartypes),     pointer :: prt                   ! PARTEH object
      integer                          :: leaf_status           ! cohort phenology status [leaves on/off]
      integer                          :: pft                   ! index for PFT
      integer                          :: iage                  ! index for leaf age loop
      integer                          :: el                    ! index for element loop
      integer                          :: element_id            ! element index consistent with defs in PRTGeneric
      integer                          :: use_pft_local(numpft) ! determine whether this PFT is used for this patch and site
      integer                          :: crown_damage          ! crown damage class of the cohort [1 = undamaged, >1 = damaged] 
      real(r8)                         :: l2fr                  ! leaf to fineroot biomass ratio [kg kg-1]
      real(r8)                         :: canopy_trim           ! fraction of the maximum leaf biomass that we are targeting [0-1]
      real(r8)                         :: cohort_n              ! cohort density
      real(r8)                         :: dbh                   ! cohort dbh [cm]
      real(r8)                         :: height                ! cohort height [m]
      real(r8)                         :: c_area                ! cohort crown area [m2]
      real(r8)                         :: c_agw                 ! above ground (non-leaf) biomass [kgC]
      real(r8)                         :: c_bgw                 ! below ground (non-fineroot) biomss [kgC]
      real(r8)                         :: c_leaf                ! leaf biomass [kgC]
      real(r8)                         :: c_fnrt                ! fine root biomss [kgC]
      real(r8)                         :: c_sapw                ! sapwood biomass [kgC]
      real(r8)                         :: c_struct              ! structural (dead) biomass [kgC]
      real(r8)                         :: c_store               ! storage biomass [kgC]
      real(r8)                         :: a_sapw                ! sapwood area [m2]
      real(r8)                         :: m_struct              ! generic (any element) mass for structure [kg]
      real(r8)                         :: m_leaf                ! generic mass for leaf  [kg]
      real(r8)                         :: m_fnrt                ! generic mass for fine-root  [kg]
      real(r8)                         :: m_sapw                ! generic mass for sapwood [kg]
      real(r8)                         :: m_store               ! generic mass for storage [kg]
      real(r8)                         :: m_repro               ! generic mass for reproductive tissues [kg]
                                                                ! of all the organs in the recruits. Used for both [kg per plant] and [kg per cohort]
      real(r8)                         :: efleaf_coh         
      real(r8)                         :: effnrt_coh 
      real(r8)                         :: efstem_coh 
      real(r8)                         :: stem_drop_fraction    ! fraction of stem to absciss when leaves absciss
      real(r8)                         :: fnrt_drop_fraction    ! fraction of fine roots to absciss when leaves absciss
      integer, parameter               :: recruitstatus = 0     ! whether the newly created cohorts are recruited or initialized
      real(r8),parameter               :: zero_co_age = 0._r8   ! The age of a newly recruited cohort is zero
      !-------------------------------------------------------------------------------------

      patch_in%tallest  => null()
      patch_in%shortest => null()

      ! if any pfts are starting with large  size  then the whole  site needs a spread of 0
      do pft = 1, numpft
         if (EDPftvarcon_inst%initd(pft) < 0.0_r8) then   
            site_in%spread = init_spread_inventory
         end if
      end do

      
      ! Manage interactions of fixed biogeog (site level filter) and nocomp (patch level filter)
      ! Need to cover all potential biogeog x nocomp combinations
      ! 1. biogeog = false. nocomp = false: all PFTs on (DEFAULT)
      ! 2. biogeog = true.  nocomp = false: site level filter
      ! 3. biogeog = false. nocomp = true : patch level filter
      ! 4. biogeog = true.  nocomp = true : patch and site level filter
      ! in principle this could be a patch level variable.
      do pft = 1, numpft
         ! first turn every PFT ON, unless we are in a special case
         use_pft_local(pft) = itrue ! Case 1
         if (hlm_use_fixed_biogeog .eq. itrue) then !filter geographically
            use_pft_local(pft) = site_in%use_this_pft(pft) ! Case 2
            if (hlm_use_nocomp .eq. itrue .and. pft .ne. patch_in%nocomp_pft_label) then
               ! having set the biogeog filter as on or off, turn off all PFTs
               ! whose identity does not correspond to this patch label
               use_pft_local(pft) = ifalse ! Case 3
            endif
         else
            if (hlm_use_nocomp .eq. itrue .and. pft .ne. patch_in%nocomp_pft_label) then
               ! This case has all PFTs on their own patch everywhere
               use_pft_local(pft) = ifalse ! Case 4
            endif
         endif
      end do

      pft_loop: do pft =  1, numpft
         if_use_this_pft: if (use_pft_local(pft) .eq. itrue) then
            l2fr         = prt_params%allom_l2fr(pft)
            canopy_trim  = 1.0_r8
            crown_damage = 1  ! Assume no damage to begin with

            ! retrieve drop fraction of non-leaf tissues for phenology initialization
            fnrt_drop_fraction = prt_params%phen_fnrt_drop_fraction(pft)
            stem_drop_fraction = prt_params%phen_stem_drop_fraction(pft)

            ! initialize phenology variables
            if_spmode: if (hlm_use_sp == itrue) then 
               ! satellite phenology: do not override SP values with build-in phenology
               efleaf_coh = 1.0_r8 
               effnrt_coh = 1.0_r8
               efstem_coh = 1.0_r8 
               leaf_status = leaves_on
            else 
               ! use built-in phenology 
               if (prt_params%season_decid(pft) == itrue .and.                 &
                  any(site_in%cstatus == [phen_cstat_nevercold, phen_cstat_iscold])) then 
                  ! Cold deciduous, off season, assume complete abscission
                  efleaf_coh = 0.0_r8
                  effnrt_coh = 1.0_r8 - fnrt_drop_fraction 
                  efstem_coh = 1.0_r8 - stem_drop_fraction
                  leaf_status = leaves_off 
               else if (any(prt_params%stress_decid(pft) == [ihard_stress_decid, isemi_stress_decid])) then 
                  ! If the plant is drought deciduous, make sure leaf status is
                  ! always consistent with the leaf elongation factor.  For tissues
                  ! other than leaves, the actual drop fraction is a combination of the
                  ! elongation factor (e) and the drop fraction (x), which will ensure
                  ! that the remaining tissue biomass will be exactly e when x=1, and
                  ! exactly the original biomass when x = 0.
                  efleaf_coh = site_in%elong_factor(pft)
                  effnrt_coh = 1.0_r8 - (1.0_r8 - efleaf_coh)*fnrt_drop_fraction
                  efstem_coh = 1.0_r8 - (1.0_r8 - efleaf_coh)*stem_drop_fraction

                  if (efleaf_coh > 0.0_r8) then 
                     leaf_status = leaves_on 
                  else 
                     leaf_status = leaves_off 
                  end if 
               else 
                  ! Evergreens, or deciduous during growing season 
                  ! Assume leaves fully flushed 
                  efleaf_coh = 1.0_r8 
                  effnrt_coh = 1.0_r8 
                  efstem_coh = 1.0_r8 
                  leaf_status = leaves_on 
               end if 
            end if if_spmode 

            ! If positive EDPftvarcon_inst%initd is interpreted as initial recruit density.
            ! If negative EDPftvarcon_inst%initd is interpreted as initial dbh. 
            ! Dbh-initialization can only be used in nocomp mode.
            ! In the dbh-initialization case, we calculate crown area for a single tree and then calculate
            ! the density of plants  needed for a full canopy. 
            if_init_dens: if (EDPftvarcon_inst%initd(pft) > nearzero) then  ! interpret as initial density and calculate diameter

               cohort_n = EDPftvarcon_inst%initd(pft)*patch_in%area
               if (hlm_use_nocomp .eq. itrue) then !in nocomp mode we only have one PFT per patch
                  ! as opposed to numpft's. So we should up the initial density
                  ! to compensate (otherwise runs are very hard to compare)
                  ! this multiplies it by the number of PFTs there would have been in
                  ! the single shared patch in competition mode.
                  ! n.b. that this is the same as currentcohort%n = %initd(pft) &AREA
                  cohort_n = cohort_n*sum(site_in%use_this_pft)
               endif
               height = EDPftvarcon_inst%hgt_min(pft)

               ! h, dbh, leafc, n from SP values or from small initial size
               if (hlm_use_sp .eq. itrue) then
                  ! At this point, we do not know the bc_in values of tlai tsai and htop,
                  ! so this is initializing to an arbitrary value for the very first timestep.
                  ! Not sure if there's a way around this or not.
                  height = 0.5_r8
                  call calculate_SP_properties(height, 0.2_r8, 0.1_r8,         &
                     patch_in%area, pft, crown_damage, 1,                      &
                     EDPftvarcon_inst%vcmax25top(pft, 1), c_leaf, dbh,         &
                     cohort_n, c_area)
               else
                  ! calculate the plant diameter from height
                  call h2d_allom(height, pft, dbh)

                  ! Calculate the leaf biomass from allometry
                  ! (calculates a maximum first, then applies canopy trim)
                  call bleaf(dbh, pft, crown_damage, canopy_trim, efleaf_coh, c_leaf)
               endif  ! sp mode

            else ! interpret as initial diameter and calculate density 
               if (hlm_use_nocomp .eq. itrue) then
                  
                  dbh = abs(EDPftvarcon_inst%initd(pft))

                  ! calculate crown area of a single plant
                  call carea_allom(dbh, 1.0_r8, init_spread_inventory, pft, crown_damage,       &
                  c_area)

                  ! calculate initial density required to close canopy 
                  cohort_n = patch_in%area/c_area

                  ! Calculate the leaf biomass from allometry
                  ! (calculates a maximum first, then applies canopy trim)
                  call bleaf(dbh, pft, crown_damage, canopy_trim, efleaf_coh,  &
                       c_leaf)

                  ! calculate crown area of the cohort
                  call carea_allom(dbh, cohort_n, init_spread_inventory, pft, crown_damage,       &
                       c_area)

                  ! calculate height from diameter
                  call h_allom(dbh, pft, height)
 
               else
                  write(fates_log(),*) 'Negative fates_recruit_init_density can only be used in no comp mode'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               endif
            endif if_init_dens 

            ! calculate total above-ground biomass from allometry
            call bagw_allom(dbh, pft, crown_damage, efstem_coh, c_agw)

            ! calculate coarse root biomass from allometry
            call bbgw_allom(dbh, pft, efstem_coh, c_bgw)

            ! Calculate fine root biomass from allometry
            ! (calculates a maximum and then trimming value)
            call bfineroot(dbh, pft, canopy_trim, l2fr, effnrt_coh, c_fnrt)

            ! Calculate sapwood biomass
            call bsap_allom(dbh, pft, crown_damage, canopy_trim, efstem_coh,   &
               a_sapw, c_sapw)

            call bdead_allom(c_agw, c_bgw, c_sapw, pft, c_struct)
            call bstore_allom(dbh, pft, crown_damage, canopy_trim, c_store)

            if (debug) write(fates_log(),*) 'EDInitMod.F90 call create_cohort '

            ! --------------------------------------------------------------------------------
            ! Initialize the mass of every element in every organ of the organ
            ! --------------------------------------------------------------------------------

            prt => null()
            call InitPRTObject(prt)

            element_loop: do el = 1, num_elements

               element_id = element_list(el)
               ! If this is carbon12, then the initialization is straight forward
               ! otherwise, we use stoichiometric ratios
               select case(element_id)
               case(carbon12_element)
                  m_struct = c_struct
                  m_leaf   = c_leaf
                  m_fnrt   = c_fnrt
                  m_sapw   = c_sapw
                  m_store  = c_store
                  m_repro  = 0._r8
               case(nitrogen_element)
                  m_struct = c_struct*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(struct_organ))
                  m_leaf   = c_leaf*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ))
                  m_fnrt   = c_fnrt*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))
                  m_sapw   = c_sapw*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
                  m_repro  = 0._r8
                  m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
               case(phosphorus_element)

                  m_struct = c_struct*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(struct_organ))
                  m_leaf   = c_leaf*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(leaf_organ))
                  m_fnrt   = c_fnrt*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))
                  m_sapw   = c_sapw*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
                  m_repro  = 0._r8
                  m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
               end select

               select case(hlm_parteh_mode)
               case (prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp )
                  ! Put all of the leaf mass into the first bin
                  call SetState(prt, leaf_organ, element_id, m_leaf, 1)
                  do iage = 2,nleafage
                     call SetState(prt, leaf_organ, element_id, 0._r8, iage)
                  end do
                  call SetState(prt, fnrt_organ, element_id, m_fnrt)
                  call SetState(prt, sapw_organ, element_id, m_sapw)
                  call SetState(prt, store_organ, element_id, m_store)
                  call SetState(prt, struct_organ, element_id, m_struct)
                  call SetState(prt, repro_organ, element_id, m_repro)

               case default
                  write(fates_log(),*) 'Unspecified PARTEH module during create_cohort'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end select

            end do element_loop

            call prt%CheckInitialConditions()

            call create_cohort(site_in, patch_in, pft, cohort_n,               &
               height, zero_co_age, dbh, prt, efleaf_coh,                      &
               effnrt_coh, efstem_coh, leaf_status, recruitstatus,             &
               canopy_trim, c_area, 1, crown_damage, site_in%spread, bc_in)

         endif if_use_this_pft
      enddo pft_loop

      if (hlm_use_sp == ifalse) then
        call fuse_cohorts(site_in, patch_in,bc_in)
        call patch_in%SortCohorts()
      end if
      
      call patch_in%ValidateCohorts()

   end subroutine init_cohorts

  ! ======================================================================================

end module EDInitMod
