module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : ifalse
  use FatesConstantsMod         , only : itrue
  use FatesConstantsMod         , only : fates_unset_int
  use FatesConstantsMod         , only : primaryforest
  use FatesGlobals              , only : endrun => fates_endrun
  use EDTypesMod                , only : nclmax
  use FatesGlobals              , only : fates_log
  use FatesInterfaceTypesMod         , only : hlm_is_restart
  use EDPftvarcon               , only : EDPftvarcon_inst
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDCohortDynamicsMod       , only : InitPRTObject
  use EDPatchDynamicsMod        , only : create_patch
  use ChecksBalancesMod         , only : SiteMassStock
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod                , only : numWaterMem
  use EDTypesMod                , only : num_vegtemp_mem
  use EDTypesMod                , only : maxpft
  use EDTypesMod                , only : AREA
  use EDTypesMod                , only : init_spread_near_bare_ground
  use EDTypesMod                , only : init_spread_inventory
  use EDTypesMod                , only : leaves_on
  use EDTypesMod                , only : leaves_off
  use EDTypesMod                , only : num_elements
  use EDTypesMod                , only : element_list
  use EDTypesMod                , only : phen_cstat_nevercold
  use EDTypesMod                , only : phen_cstat_iscold
  use EDTypesMod                , only : phen_dstat_timeoff
  use EDTypesMod                , only : phen_dstat_moistoff
  use EDTypesMod                , only : phen_cstat_notcold
  use EDTypesMod                , only : phen_dstat_moiston
  use EDTypesMod                , only : element_pos
  use FatesInterfaceTypesMod         , only : bc_in_type
  use FatesInterfaceTypesMod         , only : hlm_use_planthydro
  use FatesInterfaceTypesMod         , only : hlm_use_inventory_init
  use FatesInterfaceTypesMod         , only : numpft
  use FatesInterfaceTypesMod         , only : nleafage
  use FatesInterfaceTypesMod         , only : nlevsclass
  use FatesInterfaceTypesMod         , only : nlevcoage
  use FatesAllometryMod         , only : h2d_allom
  use FatesAllometryMod         , only : bagw_allom
  use FatesAllometryMod         , only : bbgw_allom
  use FatesAllometryMod         , only : bleaf
  use FatesAllometryMod         , only : bfineroot
  use FatesAllometryMod         , only : bsap_allom
  use FatesAllometryMod         , only : bdead_allom
  use FatesAllometryMod         , only : bstore_allom

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

  ! CIME GLOBALS
  use shr_log_mod               , only : errMsg => shr_log_errMsg

  implicit none
  private

  logical   ::  debug = .false.

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

  subroutine init_site_vars( site_in, bc_in )
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) :: site_in
    type(bc_in_type),intent(in)       :: bc_in
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
   
    site_in%nlevsoil   = bc_in%nlevsoil
    allocate(site_in%rootfrac_scr(site_in%nlevsoil))
    allocate(site_in%zi_soil(0:site_in%nlevsoil))
    allocate(site_in%dz_soil(site_in%nlevsoil))
    allocate(site_in%z_soil(site_in%nlevsoil))

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
    site_in%nchilldays       = fates_unset_int
    site_in%ncolddays        = fates_unset_int
    site_in%cleafondate      = fates_unset_int  ! doy of leaf on
    site_in%cleafoffdate     = fates_unset_int  ! doy of leaf off
    site_in%dleafondate      = fates_unset_int  ! doy of leaf on drought
    site_in%dleafoffdate     = fates_unset_int  ! doy of leaf on drought
    site_in%water_memory(:)  = nan
    site_in%vegtemp_memory(:) = nan              ! record of last 10 days temperature for senescence model.


    ! FIRE 
    site_in%acc_ni           = 0.0_r8     ! daily nesterov index accumulating over time. time unlimited theoretically.
    site_in%NF               = 0.0_r8     ! daily lightning strikes per km2 
    site_in%frac_burnt       = 0.0_r8     ! burn area read in from external file

    do el=1,num_elements
       ! Zero the state variables used for checking mass conservation
       call site_in%mass_balance(el)%ZeroMassBalState()
       call site_in%mass_balance(el)%ZeroMassBalFlux()
       call site_in%flux_diags(el)%ZeroFluxDiags()
    end do
       

    ! termination and recruitment info
    site_in%term_nindivs_canopy(:,:) = 0._r8
    site_in%term_nindivs_ustory(:,:) = 0._r8
    site_in%term_carbonflux_canopy = 0._r8
    site_in%term_carbonflux_ustory = 0._r8
    site_in%recruitment_rate(:) = 0._r8
    site_in%imort_rate(:,:) = 0._r8
    site_in%imort_carbonflux = 0._r8
    site_in%fmort_rate_canopy(:,:) = 0._r8
    site_in%fmort_rate_ustory(:,:) = 0._r8
    site_in%fmort_carbonflux_canopy = 0._r8
    site_in%fmort_carbonflux_ustory = 0._r8
    site_in%fmort_rate_cambial(:,:) = 0._r8
    site_in%fmort_rate_crown(:,:) = 0._r8

    ! fusoin-induced growth flux of individuals
    site_in%growthflux_fusion(:,:) = 0._r8

    ! demotion/promotion info
    site_in%demotion_rate(:) = 0._r8
    site_in%demotion_carbonflux = 0._r8
    site_in%promotion_rate(:) = 0._r8
    site_in%promotion_carbonflux = 0._r8
    
    ! Resources management (logging/harvesting, etc)
    site_in%resources_management%trunk_product_site  = 0.0_r8

    ! canopy spread
    site_in%spread = 0._r8

  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties( nsites, sites )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    

    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)

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

          sites(s)%cleafondate   = cleafon
          sites(s)%cleafoffdate  = cleafoff
          sites(s)%dleafoffdate  = dleafoff
          sites(s)%dleafondate   = dleafon
          sites(s)%grow_deg_days = GDD
          
          sites(s)%water_memory(1:numWaterMem) = watermem
          sites(s)%vegtemp_memory(1:num_vegtemp_mem) = 0._r8
          
          sites(s)%cstatus = cstat
          sites(s)%dstatus = dstat
          
          sites(s)%acc_NI     = acc_NI
          sites(s)%NF         = 0.0_r8         
          sites(s)%frac_burnt = 0.0_r8

          
       end do

    end if

    return
  end subroutine set_site_properties


  ! ============================================================================
  subroutine init_patches( nsites, sites, bc_in)
     !
     ! !DESCRIPTION:
     ! initialize patches
     ! This may be call a near bare ground initialization, or it may
     ! load patches from an inventory.

     !
     

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

     ! dummy locals
     real(r8) :: biomass_stock
     real(r8) :: litter_stock
     real(r8) :: seed_stock
     
     type(ed_site_type),  pointer :: sitep
     type(ed_patch_type), pointer :: newp

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

        !FIX(SPM,032414) clean this up...inits out of this loop
        do s = 1, nsites

           ! Initialize the site-level crown area spread factor (0-1)
           ! It is likely that closed canopy forest inventories
           ! have smaller spread factors than bare ground (they are crowded)
           sites(s)%spread     = init_spread_near_bare_ground

           allocate(newp)

           newp%patchno = 1
           newp%younger => null()
           newp%older   => null()

           sites(s)%youngest_patch => newp
           sites(s)%youngest_patch => newp
           sites(s)%oldest_patch   => newp

           ! make new patch...

           call create_patch(sites(s), newp, age, area, primaryforest)
           
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
           call init_cohorts(sitep, newp, bc_in(s))
           
           ! For carbon balance checks, we need to initialize the 
           ! total carbon stock
           do el=1,num_elements
              call SiteMassStock(sites(s),el,sites(s)%mass_balance(el)%old_stock, &
                   biomass_stock,litter_stock,seed_stock)
           end do
        enddo

     end if

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
    integer  :: iage       ! index for leaf age loop
    integer  :: el         ! index for element loop
    integer  :: element_id ! element index consistent with defs in PRTGeneric
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

    !----------------------------------------------------------------------

    patch_in%tallest  => null()
    patch_in%shortest => null()
    
    do pft =  1,numpft

       if(EDPftvarcon_inst%initd(pft)>1.0E-7) then

       allocate(temp_cohort) ! temporary cohort

       temp_cohort%pft         = pft
       temp_cohort%n           = EDPftvarcon_inst%initd(pft) * patch_in%area
       temp_cohort%hite        = EDPftvarcon_inst%hgt_min(pft)
       

       ! Calculate the plant diameter from height
       call h2d_allom(temp_cohort%hite,pft,temp_cohort%dbh)

       temp_cohort%canopy_trim = 1.0_r8

       ! Calculate total above-ground biomass from allometry
       call bagw_allom(temp_cohort%dbh,pft,c_agw)

       ! Calculate coarse root biomass from allometry
       call bbgw_allom(temp_cohort%dbh,pft,c_bgw)

       ! Calculate the leaf biomass from allometry
       ! (calculates a maximum first, then applies canopy trim)
       call bleaf(temp_cohort%dbh,pft,temp_cohort%canopy_trim,c_leaf)

       ! Calculate fine root biomass from allometry
       ! (calculates a maximum and then trimming value)
       call bfineroot(temp_cohort%dbh,pft,temp_cohort%canopy_trim,c_fnrt)

       ! Calculate sapwood biomass
       call bsap_allom(temp_cohort%dbh,pft,temp_cohort%canopy_trim,a_sapw,c_sapw)
       
       call bdead_allom( c_agw, c_bgw, c_sapw, pft, c_struct )

       call bstore_allom(temp_cohort%dbh, pft, temp_cohort%canopy_trim, c_store)

       temp_cohort%laimemory = 0._r8
       temp_cohort%sapwmemory = 0._r8
       temp_cohort%structmemory = 0._r8
       cstatus = leaves_on
       
       stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(temp_cohort%pft)
       
       if( EDPftvarcon_inst%season_decid(pft) == itrue .and. &
            any(site_in%cstatus == [phen_cstat_nevercold,phen_cstat_iscold])) then
         temp_cohort%laimemory = c_leaf
         temp_cohort%sapwmemory = c_sapw * stem_drop_fraction
         temp_cohort%structmemory = c_struct * stem_drop_fraction
         c_leaf = 0._r8
         c_sapw = (1.0_r8-stem_drop_fraction) * c_sapw
         c_struct  = (1.0_r8-stem_drop_fraction) * c_struct
         cstatus = leaves_off
       endif

       if ( EDPftvarcon_inst%stress_decid(pft) == itrue .and. &
            any(site_in%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff])) then
          temp_cohort%laimemory = c_leaf
          temp_cohort%sapwmemory = c_sapw * stem_drop_fraction
          temp_cohort%structmemory = c_struct * stem_drop_fraction
          c_leaf = 0._r8
          c_sapw = (1.0_r8-stem_drop_fraction) * c_sapw
          c_struct  = (1.0_r8-stem_drop_fraction) * c_struct
          cstatus = leaves_off
       endif

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
             
             m_struct = c_struct*EDPftvarcon_inst%prt_nitr_stoich_p2(pft,struct_organ)
             m_leaf   = c_leaf*EDPftvarcon_inst%prt_nitr_stoich_p2(pft,leaf_organ)
             m_fnrt   = c_fnrt*EDPftvarcon_inst%prt_nitr_stoich_p2(pft,fnrt_organ)
             m_sapw   = c_sapw*EDPftvarcon_inst%prt_nitr_stoich_p2(pft,sapw_organ)
             m_store  = c_store*EDPftvarcon_inst%prt_nitr_stoich_p2(pft,store_organ)
             m_repro  = 0._r8
             
          case(phosphorus_element)

             m_struct = c_struct*EDPftvarcon_inst%prt_phos_stoich_p2(pft,struct_organ)
             m_leaf   = c_leaf*EDPftvarcon_inst%prt_phos_stoich_p2(pft,leaf_organ)
             m_fnrt   = c_fnrt*EDPftvarcon_inst%prt_phos_stoich_p2(pft,fnrt_organ)
             m_sapw   = c_sapw*EDPftvarcon_inst%prt_phos_stoich_p2(pft,sapw_organ)
             m_store  = c_store*EDPftvarcon_inst%prt_phos_stoich_p2(pft,store_organ)
             m_repro  = 0._r8
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
            temp_cohort%coage, temp_cohort%dbh, prt_obj, temp_cohort%laimemory, &
            temp_cohort%sapwmemory, temp_cohort%structmemory, cstatus, rstatus,        &
             temp_cohort%canopy_trim, 1, site_in%spread, bc_in)

       deallocate(temp_cohort) ! get rid of temporary cohort

       endif

    enddo !numpft

    ! Zero the mass flux pools of the new cohorts
!    temp_cohort => patch_in%tallest
!    do while(associated(temp_cohort)) 
!       call temp_cohort%prt%ZeroRates()
!       temp_cohort => temp_cohort%shorter
!    end do

    call fuse_cohorts(site_in, patch_in,bc_in)
    call sort_cohorts(patch_in)

  end subroutine init_cohorts

  ! ===============================================================================================


end module EDInitMod
