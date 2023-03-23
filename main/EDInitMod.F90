module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure.
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : ifalse
  use FatesConstantsMod         , only : itrue
  use FatesConstantsMod         , only : fates_unset_int
  use FatesConstantsMod         , only : primaryforest
  use FatesConstantsMod   , only : nearzero
  use FatesGlobals              , only : endrun => fates_endrun
  use EDTypesMod                , only : nclmax
  use FatesGlobals              , only : fates_log
  use FatesInterfaceTypesMod    , only : hlm_is_restart
  use EDPftvarcon               , only : EDPftvarcon_inst
  use PRTParametersMod          , only : prt_params
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDCohortDynamicsMod       , only : InitPRTObject
  use EDPatchDynamicsMod        , only : create_patch
  use EDPatchDynamicsMod        , only : set_patchno
  use EDPhysiologyMod           , only : assign_cohort_sp_properties
  use ChecksBalancesMod         , only : SiteMassStock
  use FatesInterfaceTypesMod    , only : hlm_day_of_year
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod                , only : numWaterMem
  use EDTypesMod                , only : num_vegtemp_mem
  use EDTypesMod                , only : maxpft
  use EDTypesMod                , only : AREA
  use EDTypesMod                , only : init_spread_near_bare_ground
  use EDTypesMod                , only : init_spread_inventory
  use EDTypesMod                , only : leaves_on
  use EDTypesMod                , only : leaves_off
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
  use FatesInterfaceTypesMod         , only : numpft
  use FatesInterfaceTypesMod         , only : nleafage
  use FatesInterfaceTypesMod         , only : nlevsclass
  use FatesInterfaceTypesMod         , only : nlevcoage
  use FatesInterfaceTypesMod         , only : nlevdamage
  use FatesInterfaceTypesMod         , only : hlm_use_nocomp
  use FatesInterfaceTypesMod         , only : nlevage

  use FatesAllometryMod         , only : h2d_allom
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

  ! CIME GLOBALS
  use shr_log_mod               , only : errMsg => shr_log_errMsg

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
    allocate(site_in%term_nindivs_canopy(1:nlevsclass,1:numpft))
    allocate(site_in%term_nindivs_ustory(1:nlevsclass,1:numpft))
    allocate(site_in%demotion_rate(1:nlevsclass))
    allocate(site_in%promotion_rate(1:nlevsclass))
    allocate(site_in%imort_rate(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_canopy(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_ustory(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_cambial(1:nlevsclass,1:numpft))
    allocate(site_in%fmort_rate_crown(1:nlevsclass,1:numpft))
    allocate(site_in%growthflux_fusion(1:nlevsclass,1:numpft))
    allocate(site_in%mass_balance(1:num_elements))
    allocate(site_in%flux_diags(1:num_elements))

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

    allocate(site_in%term_carbonflux_canopy(1:numpft))
    allocate(site_in%term_carbonflux_ustory(1:numpft))
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

    if (hlm_use_nocomp .eq. itrue .and. hlm_use_fixed_biogeog .eq. itrue) then
       allocate(site_in%area_pft(0:numpft))
    else  ! SP and nocomp require a bare-ground patch.
       allocate(site_in%area_pft(1:numpft))  
    endif

    allocate(site_in%use_this_pft(1:numpft))
    allocate(site_in%area_by_age(1:nlevage))

    ! for CNP dynamics, track the mean l2fr of recruits
    ! for different pfts and canopy positions
    allocate(site_in%rec_l2fr(1:numpft,nclmax))


    ! SP mode
    allocate(site_in%sp_tlai(1:numpft))
    allocate(site_in%sp_tsai(1:numpft))
    allocate(site_in%sp_htop(1:numpft))

    do el=1,num_elements
       allocate(site_in%flux_diags(el)%leaf_litter_input(1:numpft))
       allocate(site_in%flux_diags(el)%root_litter_input(1:numpft))
    end do

    ! Initialize the static soil
    ! arrays from the boundary (initial) condition

    site_in%zi_soil(:) = bc_in%zi_sisl(:)
    site_in%dz_soil(:) = bc_in%dz_sisl(:)
    site_in%z_soil(:)  = bc_in%z_sisl(:)

    !
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
    site_in%dstatus          = fates_unset_int
    site_in%grow_deg_days    = nan  ! growing degree days
    site_in%snow_depth       = nan
    site_in%nchilldays       = fates_unset_int
    site_in%ncolddays        = fates_unset_int
    site_in%cleafondate      = fates_unset_int
    site_in%cleafoffdate     = fates_unset_int
    site_in%dleafondate      = fates_unset_int
    site_in%dleafoffdate     = fates_unset_int
    site_in%water_memory(:)  = nan
    site_in%vegtemp_memory(:) = nan              ! record of last 10 days temperature for senescence model.

    site_in%phen_model_date  = fates_unset_int

    ! Disturbance rates tracking
    site_in%primary_land_patchfusion_error = 0.0_r8
    site_in%potential_disturbance_rates(:) = 0.0_r8
    site_in%disturbance_rates_secondary_to_secondary(:) = 0.0_r8
    site_in%disturbance_rates_primary_to_secondary(:) = 0.0_r8
    site_in%disturbance_rates_primary_to_primary(:) = 0.0_r8

    ! FIRE
    site_in%acc_ni           = 0.0_r8     ! daily nesterov index accumulating over time. time unlimited theoretically.
    site_in%FDI              = 0.0_r8     ! daily fire danger index (0-1)
    site_in%NF               = 0.0_r8     ! daily lightning strikes per km2
    site_in%NF_successful    = 0.0_r8     ! daily successful iginitions per km2

    do el=1,num_elements
       ! Zero the state variables used for checking mass conservation
       call site_in%mass_balance(el)%ZeroMassBalState()
       call site_in%mass_balance(el)%ZeroMassBalFlux()
       call site_in%flux_diags(el)%ZeroFluxDiags()
    end do


    ! termination and recruitment info
    site_in%term_nindivs_canopy(:,:) = 0._r8
    site_in%term_nindivs_ustory(:,:) = 0._r8
    site_in%term_crownarea_canopy = 0._r8
    site_in%term_crownarea_ustory = 0._r8
    site_in%imort_crownarea = 0._r8
    site_in%fmort_crownarea_canopy = 0._r8
    site_in%fmort_crownarea_ustory = 0._r8
    site_in%term_carbonflux_canopy(:) = 0._r8
    site_in%term_carbonflux_ustory(:) = 0._r8
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

    site_in%area_pft(:) = 0._r8
    site_in%use_this_pft(:) = fates_unset_int
    site_in%area_by_age(:) = 0._r8

  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties( nsites, sites,bc_in )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
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
    real(r8) :: acc_NI
    real(r8) :: watermem
    integer  :: cleafon    ! DOY for cold-decid leaf-on, initial guess
    integer  :: cleafoff   ! DOY for cold-decid leaf-off, initial guess
    integer  :: dleafoff   ! DOY for drought-decid leaf-off, initial guess
    integer  :: dleafon    ! DOY for drought-decid leaf-on, initial guess
    integer  :: ft         ! PFT loop
    real(r8) :: sumarea    ! area of PFTs in nocomp mode.
    integer  :: hlm_pft    ! used in fixed biogeog mode
    integer  :: fates_pft  ! used in fixed biogeog mode
    !----------------------------------------------------------------------


    ! If this is not a restart, we need to start with some reasonable
    ! starting points. If this is a restart, we leave the values
    ! as unset ints and reals, and let the restart values be read in
    ! after this routine

    if ( hlm_is_restart == ifalse ) then

       GDD      = 30.0_r8
       cleafon  = 100
       cleafoff = 300
       cstat    = phen_cstat_notcold     ! Leaves are on
       acc_NI   = 0.0_r8
       dstat    = phen_dstat_moiston     ! Leaves are on
       dleafoff = 300
       dleafon  = 100
       watermem = 0.5_r8

       do s = 1,nsites
          sites(s)%nchilldays    = 0
          sites(s)%ncolddays     = 0        ! recalculated in phenology
          ! immediately, so yes this
          ! is memory-less, but needed
          ! for first value in history file
          sites(s)%phen_model_date = 0
          sites(s)%cleafondate     = cleafon  - hlm_day_of_year
          sites(s)%cleafoffdate    = cleafoff - hlm_day_of_year
          sites(s)%dleafoffdate    = dleafoff - hlm_day_of_year
          sites(s)%dleafondate     = dleafon  - hlm_day_of_year
          sites(s)%grow_deg_days   = GDD

          sites(s)%water_memory(1:numWaterMem) = watermem
          sites(s)%vegtemp_memory(1:num_vegtemp_mem) = 0._r8

          sites(s)%cstatus = cstat
          sites(s)%dstatus = dstat

          sites(s)%acc_NI     = acc_NI
          sites(s)%NF         = 0.0_r8
          sites(s)%NF_successful  = 0.0_r8
          sites(s)%area_pft(:) = 0.0_r8

          do ft =  1,numpft
             sites(s)%rec_l2fr(ft,:) = prt_params%allom_l2fr(ft)
          end do

          ! Its difficult to come up with a resonable starting smoothing value, so
          ! we initialize on a cold-start to -1
          sites(s)%ema_npp = -9999._r8

          if(hlm_use_fixed_biogeog.eq.itrue)then
             ! MAPPING OF FATES PFTs on to HLM_PFTs
             ! add up the area associated with each FATES PFT
             ! where pft_areafrac is the area of land in each HLM PFT and (from surface dataset)
             ! hlm_pft_map is the area of that land in each FATES PFT (from param file)

             do hlm_pft = 1,size( EDPftvarcon_inst%hlm_pft_map,2)
                do fates_pft = 1,numpft ! loop round all fates pfts for all hlm pfts
                   sites(s)%area_pft(fates_pft) = sites(s)%area_pft(fates_pft) + &
                        EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft) * bc_in(s)%pft_areafrac(hlm_pft)
                end do
             end do !hlm_pft

             do ft =  1,numpft
                if(sites(s)%area_pft(ft).lt.0.01_r8.and.sites(s)%area_pft(ft).gt.0.0_r8)then
                   if(debug) write(fates_log(),*)  'removing small pft patches',s,ft,sites(s)%area_pft(ft)
                   sites(s)%area_pft(ft)=0.0_r8
                   ! remove tiny patches to prevent numerical errors in terminate patches
                endif
                if(sites(s)%area_pft(ft).lt.0._r8)then
                   write(fates_log(),*) 'negative area',s,ft,sites(s)%area_pft(ft)
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
                sites(s)%area_pft(ft)= sites(s)%area_pft(ft) * AREA ! rescale units to m2.
             end do

             ! re-normalize PFT area to ensure it sums to one.
             ! note that in areas of 'bare ground' (PFT 0 in CLM/ELM)
             ! the bare ground will no longer be proscribed and should emerge from FATES
             ! this may or may not be the right way to deal with this?

             if(hlm_use_nocomp.eq.ifalse)then ! when not in nocomp (i.e. or SP) mode, 
                ! subsume bare ground evenly into the existing patches.

                sumarea = sum(sites(s)%area_pft(1:numpft))
                do ft =  1,numpft
                   if(sumarea.gt.0._r8)then
                      sites(s)%area_pft(ft) = area * sites(s)%area_pft(ft)/sumarea
                   else
                      sites(s)%area_pft(ft) = area/numpft
                      ! in nocomp mode where there is only bare ground, we assign equal area to
                      ! all pfts and let the model figure out whether land should be bare or not.
                   end if
                end do !ft
             else ! for sp and nocomp mode, assert a bare ground patch if needed
                sumarea = sum(sites(s)%area_pft(1:numpft))

                ! In all the other FATES modes, bareground is the area in which plants
                ! do not grow of their own accord. In SP mode we assert that the canopy is full for
                ! each PFT patch. Thus, we also need to assert a bare ground area in
                ! order to not have all of the ground filled by leaves.

                ! Further to that, one could calculate bare ground as the remaining area when
                ! all fhe canopies are accounted for, but this means we don't pass balance checks
                ! on canopy are inside FATES, and so in SP mode, we define the bare groud
                ! patch as having a PFT identifier as zero.

                if(sumarea.lt.area)then !make some bare ground
                   sites(s)%area_pft(0) = area - sumarea
                end if
             end if !sp mode
          end if !fixed biogeog

          do ft = 1,numpft
             ! Setting this to true ensures that all pfts
             ! are used for nocomp with no biogeog
             sites(s)%use_this_pft(ft) = itrue
             if(hlm_use_fixed_biogeog.eq.itrue)then
                if(sites(s)%area_pft(ft).gt.0.0_r8)then
                   sites(s)%use_this_pft(ft) = itrue
                else
                   sites(s)%use_this_pft(ft) = ifalse
                end if !area
             end if !SBG
          end do !ft
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

    ! dummy locals
    real(r8) :: biomass_stock
    real(r8) :: litter_stock
    real(r8) :: seed_stock
    integer  :: n
    integer  :: start_patch
    integer  :: num_new_patches
    integer  :: nocomp_pft
    real(r8) :: newparea
    real(r8) :: tota !check on area
    integer  :: is_first_patch

    type(ed_site_type),  pointer :: sitep
    type(ed_patch_type), pointer :: newppft(:)
    type(ed_patch_type), pointer :: newp
    type(ed_patch_type), pointer :: currentPatch

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
          end do
       enddo

    else

       do s = 1, nsites
          sites(s)%sp_tlai(:) = 0._r8
          sites(s)%sp_tsai(:) = 0._r8
          sites(s)%sp_htop(:) = 0._r8

          ! Initialize the site-level crown area spread factor (0-1)
          ! It is likely that closed canopy forest inventories
          ! have smaller spread factors than bare ground (they are crowded)
          sites(s)%spread     = init_spread_near_bare_ground

          start_patch = 1   ! start at the first vegetated patch
          if(hlm_use_nocomp.eq.itrue)then
             num_new_patches = numpft
             if( hlm_use_fixed_biogeog .eq.itrue )then
                start_patch = 0 ! start at the bare ground patch
             endif
             !           allocate(newppft(numpft))
          else !default
             num_new_patches = 1
          end if !nocomp

          is_first_patch = itrue
          do n = start_patch, num_new_patches

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
                   newparea = sites(s)%area_pft(nocomp_pft)
                else
                   newparea = area / numpft
                end if
             else  ! The default case is initialized w/ one patch with the area of the whole site.
                newparea = area
             end if  !nocomp mode

             if(newparea.gt.0._r8)then ! Stop patches being initilialized when PFT not present in nocomop mode
                allocate(newp)

                call create_patch(sites(s), newp, age, newparea, primaryforest, nocomp_pft)

                if(is_first_patch.eq.itrue)then !is this the first patch?
                   ! set poointers for first patch (or only patch, if nocomp is false)
                   newp%patchno = 1
                   newp%younger => null()
                   newp%older   => null()
                   sites(s)%youngest_patch => newp
                   sites(s)%oldest_patch   => newp
                   is_first_patch = ifalse
                else
                   ! Set pointers for N>1 patches. Note this only happens when nocomp mode s on.
                   ! The new patch is the 'youngest' one, arbitrarily.
                   newp%patchno = nocomp_pft
                   newp%older     => sites(s)%youngest_patch
                   newp%younger   => null()
                   sites(s)%youngest_patch%younger => newp
                   sites(s)%youngest_patch   => newp
                end if

                ! Initialize the litter pools to zero, these
                ! pools will be populated by looping over the existing patches
                ! and transfering in mass
                do el=1,num_elements
                   call newp%litter(el)%InitConditions(init_leaf_fines=0._r8, &
                        init_root_fines=0._r8, &
                        init_ag_cwd=0._r8, &
                        init_bg_cwd=0._r8, &
                        init_seed=0._r8,   &
                        init_seed_germ=0._r8)
                end do

                sitep => sites(s)
                if(hlm_use_sp.eq.itrue)then
                   if(nocomp_pft.ne.0)then !don't initialize cohorts for SP bare ground patch
                      call init_cohorts(sitep, newp, bc_in(s))
                   end if
                else ! normal non SP case always call init cohorts
                   call init_cohorts(sitep, newp, bc_in(s))
                end if
             end if
          end do !no new patches

          !check if the total area adds to the same as site area
          tota = 0.0_r8
          newp => sites(s)%oldest_patch
          do while (associated(newp))
             tota=tota+newp%area
             newp=>newp%younger
          end do

          if(abs(tota-area).gt.nearzero*area)then
             if(abs(tota-area).lt.1.0e-10_r8)then ! this is a precision error
                if(sites(s)%oldest_patch%area.gt.(tota-area+nearzero))then
                   ! remove or add extra area
                   ! if the oldest patch has enough area, use that
                   sites(s)%oldest_patch%area = sites(s)%oldest_patch%area - (tota-area)
                   if(debug) write(fates_log(),*) 'fixing patch precision - oldest',s, tota-area
                else ! or otherwise take the area from the youngest patch.
                   sites(s)%youngest_patch%area = sites(s)%oldest_patch%area - (tota-area)
                   if(debug) write(fates_log(),*) 'fixing patch precision -youngest ',s, tota-area
                endif
             else !this is a big error not just a precision error.
                write(fates_log(),*) 'issue with patch area in EDinit',tota-area,tota
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif  ! big error
          end if ! too much patch area

          ! For carbon balance checks, we need to initialize the
          ! total carbon stock
          do el=1,num_elements
             call SiteMassStock(sites(s),el,sites(s)%mass_balance(el)%old_stock, &
                  biomass_stock,litter_stock,seed_stock)
          end do

          call set_patchno(sites(s))

       enddo !s
    end if

    ! zero all the patch fire variables for the first timestep
    do s = 1, nsites
       currentPatch => sites(s)%youngest_patch
       do while(associated(currentPatch))

          currentPatch%litter_moisture(:)         = 0._r8
          currentPatch%fuel_eff_moist             = 0._r8
          currentPatch%livegrass                  = 0._r8
          currentPatch%sum_fuel                   = 0._r8
          currentPatch%fuel_bulkd                 = 0._r8
          currentPatch%fuel_sav                   = 0._r8
          currentPatch%fuel_mef                   = 0._r8
          currentPatch%ros_front                  = 0._r8
          currentPatch%effect_wspeed              = 0._r8
          currentPatch%tau_l                      = 0._r8
          currentPatch%fuel_frac(:)               = 0._r8
          currentPatch%tfc_ros                    = 0._r8
          currentPatch%fi                         = 0._r8
          currentPatch%fire                       = 0
          currentPatch%fd                         = 0._r8
          currentPatch%ros_back                   = 0._r8
          currentPatch%scorch_ht(:)               = 0._r8
          currentPatch%frac_burnt                 = 0._r8
          currentPatch%burnt_frac_litter(:)       = 0._r8

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

    return
  end subroutine init_patches

  ! ============================================================================
  subroutine init_cohorts( site_in, patch_in, bc_in)
    !
    ! !DESCRIPTION:
    ! initialize new cohorts on bare ground
    !
    ! !USES:

    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout),  pointer  :: site_in
    type(ed_patch_type), intent(inout), pointer  :: patch_in
    type(bc_in_type), intent(in)                 :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type),pointer :: temp_cohort
    class(prt_vartypes),pointer  :: prt_obj
    integer  :: cstatus
    integer  :: pft
    integer  :: crowndamage ! which crown damage class
    integer  :: iage       ! index for leaf age loop
    integer  :: el         ! index for element loop
    integer  :: element_id ! element index consistent with defs in PRTGeneric
    integer  :: use_pft_local(numpft) ! determine whether this PFT is used for this patch and site.
    real(r8) :: c_agw      ! biomass above ground (non-leaf)     [kgC]
    real(r8) :: c_bgw      ! biomass below ground (non-fineroot) [kgC]
    real(r8) :: c_leaf     ! biomass in leaves [kgC]
    real(r8) :: c_fnrt     ! biomass in fine roots [kgC]
    real(r8) :: c_sapw     ! biomass in sapwood [kgC]
    real(r8) :: c_struct   ! biomass in structure (dead) [kgC]
    real(r8) :: c_store    ! biomass in storage [kgC]
    real(r8) :: a_sapw     ! area in sapwood (dummy) [m2]
    real(r8) :: m_struct   ! Generic (any element) mass for structure [kg]
    real(r8) :: m_leaf     ! Generic mass for leaf  [kg]
    real(r8) :: m_fnrt     ! Generic mass for fine-root  [kg]
    real(r8) :: m_sapw     ! Generic mass for sapwood [kg]
    real(r8) :: m_store    ! Generic mass for storage [kg]
    real(r8) :: m_repro    ! Generic mass for reproductive tissues [kg]
    real(r8) :: stem_drop_fraction

    integer, parameter :: rstatus = 0
    integer init

    real(r8) :: dummy_n    ! set cohort n to a dummy value of 1
 
    !----------------------------------------------------------------------

    patch_in%tallest  => null()
    patch_in%shortest => null()

    ! Manage interactions of fixed biogeog (site level filter) and
    ! nocomp (patch level filter)
    ! Need to cover all potential biogeog x nocomp combinations
    ! 1. biogeog = false. nocomp = false: all PFTs on (DEFAULT)
    ! 2. biogeog = true.  nocomp = false: site level filter
    ! 3. biogeog = false. nocomp = true : patch level filter
    ! 4. biogeog = true.  nocomp = true : patch and site level filter
    ! in principle this could be a patch level variable.
    do pft =  1,numpft
       ! Turn every PFT ON, unless we are in a special case.
       use_pft_local(pft) = itrue ! Case 1
       if(hlm_use_fixed_biogeog.eq.itrue)then !filter geographically
          use_pft_local(pft) = site_in%use_this_pft(pft) ! Case 2
          if(hlm_use_nocomp.eq.itrue.and.pft.ne.patch_in%nocomp_pft_label)then
             ! Having set the biogeog filter as on or off, turn off all PFTs
             ! whose identiy does not correspond to this patch label.
             use_pft_local(pft) = ifalse ! Case 3
          endif
       else
          if(hlm_use_nocomp.eq.itrue.and.pft.ne.patch_in%nocomp_pft_label)then
             ! This case has all PFTs on their own patch everywhere.
             use_pft_local(pft) = ifalse ! Case 4
          endif
       endif
    end do

    
    do pft =  1,numpft

       if(use_pft_local(pft).eq.itrue)then

          allocate(temp_cohort) ! temporary cohort
          temp_cohort%pft         = pft
          temp_cohort%l2fr = prt_params%allom_l2fr(pft)
          temp_cohort%canopy_trim = 1.0_r8
          temp_cohort%crowndamage = 1  ! Assume no damage to begin with
          
          ! If positive EDPftvarcon_inst%initd is interpreted as initial recruit density.
          ! If negative EDPftvarcon_inst%initd is interpreted as initial dbh. 
          ! Dbh-initialization can only be used in nocomp mode.
          ! In the dbh-initialization case, we calculate crown area for a single tree and then calculate
          ! the density of plants  needed for a full canopy. 

          if(EDPftvarcon_inst%initd(pft)>nearzero) then   ! interpret as initial density and calculate diameter

             temp_cohort%n           = EDPftvarcon_inst%initd(pft) * patch_in%area
             if(hlm_use_nocomp.eq.itrue)then !in nocomp mode we only have one PFT per patch
                ! as opposed to numpft's. So we should up the initial density
                ! to compensate (otherwise runs are very hard to compare)
                ! this multiplies it by the number of PFTs there would have been in
                ! the single shared patch in competition mode.
                ! n.b. that this is the same as currentcohort%n = %initd(pft) &AREA
                temp_cohort%n           =  temp_cohort%n * sum(site_in%use_this_pft)
             endif
             
             !  h,dbh,leafc,n from SP values or from small initial size.
             if(hlm_use_sp.eq.itrue)then
                init = itrue
                ! At this point, we do not know the bc_in values of tlai tsai and htop,
                ! so this is initializing to an arbitrary value for the very first timestep.
                ! Not sure if there's a way around this or not.
                call assign_cohort_SP_properties(temp_cohort, 0.5_r8,0.2_r8, 0.1_r8,patch_in%area,init,c_leaf)

             else
                temp_cohort%hite        = EDPftvarcon_inst%hgt_min(pft)
                ! Calculate the plant diameter from height
                call h2d_allom(temp_cohort%hite,pft,temp_cohort%dbh)

                ! Calculate the leaf biomass from allometry
                ! (calculates a maximum first, then applies canopy trim)
                call bleaf(temp_cohort%dbh,pft,temp_cohort%crowndamage, &
                     temp_cohort%canopy_trim,c_leaf)

             endif  ! sp mode

          else ! interpret as initial diameter and calculate density 
             if(hlm_use_nocomp .eq. itrue)then
                temp_cohort%dbh = abs(EDPftvarcon_inst%initd(pft))

                ! calculate crown area of a single plant
                dummy_n = 1.0_r8 ! make n=1 to get area of one tree
               
                call carea_allom(temp_cohort%dbh, dummy_n, init_spread_inventory, temp_cohort%pft, &
                     temp_cohort%crowndamage, temp_cohort%c_area)

                ! calculate initial density required to close canopy 
                temp_cohort%n  = patch_in%area / temp_cohort%c_area

                ! Calculate the leaf biomass from allometry
                ! (calculates a maximum first, then applies canopy trim)
                call bleaf(temp_cohort%dbh,pft,temp_cohort%crowndamage, &
                     temp_cohort%canopy_trim,c_leaf)

             else
                write(fates_log(),*) 'Negative fates_recruit_init_density can only be used in no comp mode'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif
          endif


          ! Calculate total above-ground biomass from allometry
          call bagw_allom(temp_cohort%dbh,pft,temp_cohort%crowndamage,c_agw)

          ! Calculate coarse root biomass from allometry
          call bbgw_allom(temp_cohort%dbh,pft,c_bgw)

          ! Calculate fine root biomass from allometry
          ! (calculates a maximum and then trimming value)
          call bfineroot(temp_cohort%dbh,pft,temp_cohort%canopy_trim,temp_cohort%l2fr,c_fnrt)

          ! Calculate sapwood biomass
          call bsap_allom(temp_cohort%dbh,pft,temp_cohort%crowndamage, &
               temp_cohort%canopy_trim,a_sapw,c_sapw)

          call bdead_allom( c_agw, c_bgw, c_sapw, pft, c_struct )

          call bstore_allom(temp_cohort%dbh, pft, temp_cohort%crowndamage, &
               temp_cohort%canopy_trim, c_store)

          cstatus = leaves_on

          stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(temp_cohort%pft)

          if(hlm_use_sp.eq.ifalse)then ! do not override SP vales with phenology

             if( prt_params%season_decid(pft) == itrue .and. &
                  any(site_in%cstatus == [phen_cstat_nevercold,phen_cstat_iscold])) then
                c_leaf = 0._r8
                c_sapw = (1.0_r8-stem_drop_fraction) * c_sapw
                c_struct  = (1.0_r8-stem_drop_fraction) * c_struct
                cstatus = leaves_off
             endif

             if ( prt_params%stress_decid(pft) == itrue .and. &
                  any(site_in%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff])) then
                c_leaf = 0._r8
                c_sapw = (1.0_r8-stem_drop_fraction) * c_sapw
                c_struct  = (1.0_r8-stem_drop_fraction) * c_struct
                cstatus = leaves_off
             endif

          end if ! SP mode

          if ( debug ) write(fates_log(),*) 'EDInitMod.F90 call create_cohort '

          temp_cohort%coage = 0.0_r8


          ! --------------------------------------------------------------------------------
          ! Initialize the mass of every element in every organ of the organ
          ! --------------------------------------------------------------------------------

          prt_obj => null()
          call InitPRTObject(prt_obj)

          do el = 1,num_elements

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

                m_struct = c_struct*prt_params%nitr_stoich_p1(pft,prt_params%organ_param_id(struct_organ))
                m_leaf   = c_leaf*prt_params%nitr_stoich_p1(pft,prt_params%organ_param_id(leaf_organ))
                m_fnrt   = c_fnrt*prt_params%nitr_stoich_p1(pft,prt_params%organ_param_id(fnrt_organ))
                m_sapw   = c_sapw*prt_params%nitr_stoich_p1(pft,prt_params%organ_param_id(sapw_organ))
                m_repro  = 0._r8
                m_store  = StorageNutrientTarget(pft,element_id,m_leaf,m_fnrt,m_sapw,m_struct)

             case(phosphorus_element)

                m_struct = c_struct*prt_params%phos_stoich_p1(pft,prt_params%organ_param_id(struct_organ))
                m_leaf   = c_leaf*prt_params%phos_stoich_p1(pft,prt_params%organ_param_id(leaf_organ))
                m_fnrt   = c_fnrt*prt_params%phos_stoich_p1(pft,prt_params%organ_param_id(fnrt_organ))
                m_sapw   = c_sapw*prt_params%phos_stoich_p1(pft,prt_params%organ_param_id(sapw_organ))
                m_repro  = 0._r8
                m_store  = StorageNutrientTarget(pft,element_id,m_leaf,m_fnrt,m_sapw,m_struct)

             end select

             select case(hlm_parteh_mode)
             case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp )

                ! Put all of the leaf mass into the first bin
                call SetState(prt_obj,leaf_organ, element_id,m_leaf,1)
                do iage = 2,nleafage
                   call SetState(prt_obj,leaf_organ, element_id,0._r8,iage)
                end do

                call SetState(prt_obj,fnrt_organ, element_id, m_fnrt)
                call SetState(prt_obj,sapw_organ, element_id, m_sapw)
                call SetState(prt_obj,store_organ, element_id, m_store)
                call SetState(prt_obj,struct_organ, element_id, m_struct)
                call SetState(prt_obj,repro_organ, element_id, m_repro)

             case default
                write(fates_log(),*) 'Unspecified PARTEH module during create_cohort'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end select

          end do

          call prt_obj%CheckInitialConditions()

          call create_cohort(site_in, patch_in, pft, temp_cohort%n, temp_cohort%hite, &
               temp_cohort%coage, temp_cohort%dbh, prt_obj, cstatus, rstatus,        &
               temp_cohort%canopy_trim, temp_cohort%c_area,1,temp_cohort%crowndamage, site_in%spread, bc_in)

          deallocate(temp_cohort, stat=istat, errmsg=smsg)
          if (istat/=0) then
             write(fates_log(),*) 'dealloc014: fail on deallocate(temp_cohort):'//trim(smsg)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

       endif !use_this_pft
    enddo !numpft

    ! (Keeping as an example)
    ! Pass patch level temperature to the new cohorts (this is a nominal 15C right now)
    !temp_cohort => patch_in%tallest
    !do while(associated(temp_cohort))
    !call temp_cohort%tveg_lpa%UpdateRmean(patch_in%tveg_lpa%GetMean())
    !temp_cohort => temp_cohort%shorter
    !end do

    call fuse_cohorts(site_in, patch_in,bc_in)
    call sort_cohorts(patch_in)


  end subroutine init_cohorts

  ! ===============================================================================================


end module EDInitMod
