module FatesSoilBGCFluxMod

  ! ============================================================================
  ! This module contains the routines that handle nutrient and carbon fluxes
  ! and states between the land-model's soil and FATES plants (uptake and plant
  ! characteristics for aquisition), and sending fragmented litter in FATES to
  ! the land-model's litter pool (which may include plant efflux).  
  ! ============================================================================

  use FatesConstantsMod, only    : r8 => fates_r8
  use FatesConstantsMod, only    : nearzero
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use PRTGenericMod     , only : num_elements
  use PRTGenericMod     , only : element_list
  use PRTGenericMod     , only : element_pos
  use PRTGenericMod     , only : prt_carbon_allom_hyp
  use PRTGenericMod     , only : prt_cnp_flex_allom_hyp
  use PRTGenericMod     , only : prt_vartypes
  use PRTGenericMod     , only : leaf_organ
  use PRTGenericMod     , only : sapw_organ, struct_organ
  use PRTGenericMod     , only : all_carbon_elements
  use PRTGenericMod     , only : carbon12_element
  use PRTGenericMod     , only : nitrogen_element
  use PRTGenericMod     , only : phosphorus_element
  use PRTGenericMod     , only : leaf_organ
  use PRTGenericMod     , only : fnrt_organ
  use PRTGenericMod     , only : sapw_organ
  use PRTGenericMod     , only : store_organ
  use PRTGenericMod     , only : repro_organ
  use PRTGenericMod     , only : struct_organ
  use PRTGenericMod     , only : SetState
  use PRTAllometricCNPMod,only : stoich_max
  use FatesAllometryMod, only : set_root_fraction
  use FatesAllometryMod , only : h_allom
  use FatesAllometryMod , only : h2d_allom
  use FatesAllometryMod , only : bagw_allom
  use FatesAllometryMod , only : bsap_allom
  use FatesAllometryMod , only : bleaf
  use FatesAllometryMod , only : bdead_allom
  use FatesAllometryMod , only : bstore_allom
  use FatesAllometryMod , only : bbgw_allom
  use FatesAllometryMod , only : carea_allom
  use EDTypesMod        , only : p_uptake_mode
  use EDTypesMod        , only : n_uptake_mode
  use EDTypesMod          , only : ed_site_type
  use EDTypesMod          , only : ed_patch_type
  use EDTypesMod          , only : ed_cohort_type
  use EDTypesMod          , only : AREA,AREA_INV
  use FatesInterfaceTypesMod, only    : bc_in_type
  use FatesInterfaceTypesMod, only    : bc_out_type
  use FatesInterfaceTypesMod, only    : numpft
  use FatesInterfaceTypesMod, only    : hlm_nu_com
  use FatesInterfaceTypesMod, only    : hlm_parteh_mode
  use FatesInterfaceTypesMod, only    : hlm_use_ch4
  use FatesConstantsMod , only : prescribed_p_uptake
  use FatesConstantsMod , only : prescribed_n_uptake
  use FatesConstantsMod , only : coupled_p_uptake
  use FatesConstantsMod , only : coupled_n_uptake
  use FatesConstantsMod, only    : days_per_sec
  use FatesConstantsMod, only    : g_per_kg
  use FatesConstantsMod, only    : kg_per_g
  use FatesConstantsMod, only    : fates_np_comp_scaling
  use FatesConstantsMod, only    : coupled_np_comp_scaling
  use FatesConstantsMod, only    : trivial_np_comp_scaling
  use FatesConstantsMod, only    : rsnbl_math_prec
  use FatesConstantsMod, only    : days_per_year
  use FatesConstantsMod, only    : sec_per_day
  use FatesConstantsMod, only    : years_per_day
  use FatesConstantsMod, only    : itrue
  use FatesLitterMod,        only : litter_type
  use FatesLitterMod    , only : ncwd
  use FatesLitterMod    , only : ndcmpy
  use FatesLitterMod    , only : ilabile
  use FatesLitterMod    , only : ilignin
  use FatesLitterMod    , only : icellulose
  use PRTParametersMod , only    : prt_params
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use FatesUtilsMod, only : check_var_real
  
  implicit none
  private

  public :: PrepCH4Bcs
  public :: PrepNutrientAquisitionBCs
  public :: UnPackNutrientAquisitionBCs
  public :: FluxIntoLitterPools

  
  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
       __FILE__


contains
  
  ! =====================================================================================

  function GetPlantDemand(ccohort,element_id) result(plant_demand)

    ! -----------------------------------------------------------------------------------
    ! This function calculates the plant's demand for a given nutrient
    ! based upon the need to fill its NPP demand and/or the need to
    ! get its tissues to their ideal stoichiometry ratios.
    ! This routine is used for informing BGC competition schemes, and
    ! for generating synthetic upake rates and also for calculating
    ! diagnostics.
    !
    ! THIS ROUTINE IS UNDERGOING MODIFICATIONS WILL CLEAN UP AFTER
    ! A DECENT FIRST HYPOTHESIS MANIFESTS
    ! -----------------------------------------------------------------------------------

    
    type(ed_cohort_type),intent(in) :: ccohort
    integer,intent(in)              :: element_id      ! Should match nitrogen_element or
                                                       ! phosphorus_element
    
    real(r8)              :: plant_demand ! Nutrient demand per plant [kg]
    real(r8)              :: plant_x      ! Total mass for element of interest [kg]
    real(r8)              :: plant_max_x  ! Maximum mass for element of interest [kg]
    integer               :: pft
    real(r8)              :: dbh
    real(r8)              :: fnrt_c

    real(r8), parameter :: smth_fac = 0.1_r8         ! Smoothing factor for updating
                                                     ! demand.
    real(r8), parameter :: init_demand_frac = 0.1_r8 ! Newly recruited plants
                                                     ! will specify a demand
                                                     ! based on their total nutrient
                                                     ! because they have not history
                                                     ! of need yet

    
    
    pft = ccohort%pft
    dbh = ccohort%dbh


    ! If the cohort has not experienced a day of integration
    ! (and thus any allocation yet), it has no deficit
    ! in its storage to drive any need, so it thus has no demand
    if(ccohort%isnew) then
       plant_demand = 0._r8
       return
    end if
       
    fnrt_c = ccohort%prt%GetState(fnrt_organ, carbon12_element)
    
    ! If the plant is not a newly recruited plant
    ! We use other methods of specifying nutrient demand
    ! -----------------------------------------------------------------------------------

    if(element_id.eq.nitrogen_element) then

       plant_demand = fnrt_c * EDPftvarcon_inst%eca_vmax_nh4(ccohort%pft) * sec_per_day
       
    elseif(element_id.eq.phosphorus_element) then

       plant_demand = fnrt_c * EDPftvarcon_inst%eca_vmax_p(ccohort%pft) * sec_per_day
       
    end if

  end function GetPlantDemand

  ! =====================================================================================
  
  
  subroutine UnPackNutrientAquisitionBCs(sites, bc_in)

    ! -----------------------------------------------------------------------------------
    ! The purpose of this routine is to recieve the nutrient uptake flux
    ! boundary conditions, and parse those fluxes out to the cohorts.
    !
    ! This routine should be called before FATES dynamics, particularly before
    ! any of the PARTEH code is called.  It is assumed that these uptake fluxes
    ! are being incremented each short-BGC timestep over the course of the day, and
    ! thus should be an integrated quantity, total nutrient uptake, over 1 day.
    ! At the end of this routine, after we have parsed the uptake to the cohorts,
    ! we can then zero out the input boundary conditions again, so they can be
    ! integrated.
    ! -----------------------------------------------------------------------------------


    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) :: sites(:)
    type(bc_in_type), intent(in)     :: bc_in(:)

    ! Locals
    integer                       :: nsites        ! number of sites
    integer                       :: s             ! site loop index
    integer                       :: j             ! soil layer
    integer                       :: icomp         ! competitor index
    integer                       :: id            ! decomp layer index
    integer                       :: pft           ! pft index
    type(ed_patch_type), pointer  :: cpatch        ! current patch pointer
    type(ed_cohort_type), pointer :: ccohort       ! current cohort pointer
    real(r8) :: fnrt_c                             ! fine-root carbon [kg]

    nsites = size(sites,dim=1)

    ! We can exit if this is a c-only simulation
    if(hlm_parteh_mode.eq.prt_carbon_allom_hyp) then
       ! These can now be zero'd
       do s = 1, nsites
          bc_in(s)%plant_nh4_uptake_flux(:,:) = 0._r8
          bc_in(s)%plant_no3_uptake_flux(:,:) = 0._r8
          bc_in(s)%plant_p_uptake_flux(:,:) = 0._r8
       end do
       return
    end if

    do s = 1, nsites

       ! If the plant is in "prescribed uptake mode"
       ! then we are not coupling with the soil bgc model.
       ! In this case, the bc_in structure is meaningless.
       ! Instead, we give the plants a parameterized fraction
       ! of their demand.  Routine GetPlantDemand() returns
       ! the plant demand.

       if (n_uptake_mode.eq.prescribed_n_uptake) then

          cpatch => sites(s)%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                pft = ccohort%pft
                ccohort%daily_n_demand = GetPlantDemand(ccohort,nitrogen_element)
                ccohort%daily_nh4_uptake = EDPftvarcon_inst%prescribed_nuptake(pft) * ccohort%daily_n_demand
                ccohort%daily_no3_uptake = 0._r8
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do

       elseif(n_uptake_mode.eq.coupled_n_uptake) then
          
          icomp = 0
          cpatch => sites(s)%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                icomp = icomp+1
                ccohort%daily_n_demand = GetPlantDemand(ccohort,nitrogen_element)
                ! N Uptake:  Convert g/m2/day -> kg/plant/day
                ccohort%daily_nh4_uptake = bc_in(s)%plant_nh4_uptake_flux(icomp,1)*kg_per_g*AREA/ccohort%n
                ccohort%daily_no3_uptake = bc_in(s)%plant_no3_uptake_flux(icomp,1)*kg_per_g*AREA/ccohort%n
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
             
       end if
       
       if (p_uptake_mode.eq.prescribed_p_uptake) then
          cpatch => sites(s)%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                pft = ccohort%pft
                ccohort%daily_p_demand = GetPlantDemand(ccohort,phosphorus_element)
                ccohort%daily_p_uptake = EDPftvarcon_inst%prescribed_puptake(pft) * ccohort%daily_p_demand
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
          
       elseif(p_uptake_mode.eq.coupled_p_uptake) then

          icomp = 0
          cpatch => sites(s)%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                icomp = icomp+1
                ccohort%daily_p_demand =  GetPlantDemand(ccohort,phosphorus_element)
                ! P Uptake:  Convert g/m2/day -> kg/plant/day
                ccohort%daily_p_uptake = bc_in(s)%plant_p_uptake_flux(icomp,1)*kg_per_g*AREA/ccohort%n
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
          
       end if
          
       ! These can now be zero'd
       bc_in(s)%plant_nh4_uptake_flux(:,:) = 0._r8
       bc_in(s)%plant_no3_uptake_flux(:,:) = 0._r8
       bc_in(s)%plant_p_uptake_flux(:,:) = 0._r8

    end do
    return
  end subroutine UnPackNutrientAquisitionBCs

  ! =====================================================================================

  subroutine PrepCH4BCs(csite,bc_in,bc_out)
    
    !
    ! This routine prepares the output boundary conditions for methane calculations
    ! in ELM/CLM.
    ! -----------------------------------------------------------------------------------
    
    
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) :: csite
    
    type(bc_out_type), intent(inout)  :: bc_out
    type(bc_in_type), intent(in)  :: bc_in
    type(ed_patch_type), pointer  :: cpatch        ! current patch pointer
    type(ed_cohort_type), pointer :: ccohort       ! current cohort pointer
    integer                       :: pft           ! plant functional type
    integer                       :: fp            ! patch index of the site
    real(r8) :: agnpp   ! Above ground daily npp
    real(r8) :: bgnpp   ! Below ground daily npp
    real(r8) :: plant_area ! crown area (m2) of all plants in patch
    real(r8) :: woody_area ! corwn area (m2) of woody plants in patch
    real(r8) :: fnrt_c  ! Fine root carbon [kg/plant]
    real(r8) :: sapw_net_alloc
    real(r8) :: store_net_alloc
    real(r8) :: fnrt_net_alloc
    real(r8) :: leaf_net_alloc
    real(r8) :: struct_net_alloc
    real(r8) :: repro_net_alloc

    ! Exit if we need not communicate with the hlm's ch4 module
    if(.not.(hlm_use_ch4==itrue)) return
    
    ! Initialize to zero
    bc_out%annavg_agnpp_pa(:) = 0._r8
    bc_out%annavg_bgnpp_pa(:) = 0._r8
    bc_out%annsum_npp_pa(:)   = 0._r8
    bc_out%rootfr_pa(:,:)  = 0._r8
    bc_out%frootc_pa(:)    = 0._r8
    bc_out%root_resp(:)  = 0._r8
    bc_out%woody_frac_aere_pa(:) = 0._r8
    
    fp = 0
    cpatch => csite%oldest_patch
    do while (associated(cpatch))

       ! Patch ordering when passing boundary conditions
       ! always goes from oldest to youngest, following
       ! the convention of EDPatchDynamics::set_patchno()
       
       fp    = fp + 1
       
       agnpp = 0._r8
       bgnpp = 0._r8
       woody_area = 0._r8
       plant_area = 0._r8
       
       ccohort => cpatch%tallest
       do while (associated(ccohort))

          ! For consistency, only apply calculations to non-new
          ! cohorts. New cohorts will not have respiration rates
          ! at this point in the call sequence.
          
          if(.not.ccohort%isnew) then
             
             pft   = ccohort%pft

             call set_root_fraction(csite%rootfrac_scr, pft, csite%zi_soil, &
                  bc_in%max_rooting_depth_index_col )

             fnrt_c   = ccohort%prt%GetState(fnrt_organ, carbon12_element)

             ! Fine root fraction over depth

             bc_out%rootfr_pa(fp,1:bc_in%nlevsoil) = &
                  bc_out%rootfr_pa(fp,1:bc_in%nlevsoil) + &
                  csite%rootfrac_scr(1:bc_in%nlevsoil)

             ! Fine root carbon, convert [kg/plant] -> [g/m2]
             bc_out%frootc_pa(fp) = &
                  bc_out%frootc_pa(fp) + &
                  fnrt_c*ccohort%n/cpatch%area * g_per_kg

             ! [kgC/day]
             sapw_net_alloc   = ccohort%prt%GetNetAlloc(sapw_organ, carbon12_element) * days_per_sec
             store_net_alloc  = ccohort%prt%GetNetAlloc(store_organ, carbon12_element) * days_per_sec
             leaf_net_alloc   = ccohort%prt%GetNetAlloc(leaf_organ, carbon12_element) * days_per_sec
             fnrt_net_alloc   = ccohort%prt%GetNetAlloc(fnrt_organ, carbon12_element) * days_per_sec
             struct_net_alloc = ccohort%prt%GetNetAlloc(struct_organ, carbon12_element) * days_per_sec
             repro_net_alloc  = ccohort%prt%GetNetAlloc(repro_organ, carbon12_element) * days_per_sec

             ! [kgC/plant/day] -> [gC/m2/s]
             agnpp = agnpp + ccohort%n/cpatch%area * (leaf_net_alloc + repro_net_alloc + &
                  prt_params%allom_agb_frac(pft)*(sapw_net_alloc+store_net_alloc+struct_net_alloc)) * g_per_kg

             ! [kgC/plant/day] -> [gC/m2/s]
             bgnpp = bgnpp + ccohort%n/cpatch%area * (fnrt_net_alloc  + &
                  (1._r8-prt_params%allom_agb_frac(pft))*(sapw_net_alloc+store_net_alloc+struct_net_alloc)) * g_per_kg

             ! (gC/m2/s) root respiration (fine root MR + total root GR)
             ! RGK: We do not save root respiration and average over the day. Until we do
             !      this is a best (bad) guess at fine root MR + total root GR
             !      (kgC/indiv/yr) -> gC/m2/s
             bc_out%root_resp(1:bc_in%nlevsoil) = bc_out%root_resp(1:bc_in%nlevsoil) + &
                  ccohort%resp_acc_hold*years_per_day*g_per_kg*days_per_sec* &
                  ccohort%n*area_inv*(1._r8-prt_params%allom_agb_frac(pft)) * csite%rootfrac_scr(1:bc_in%nlevsoil)

             if( prt_params%woody(pft)==itrue ) then
                woody_area = woody_area + ccohort%c_area
             end if
             plant_area = plant_area + ccohort%c_area
             
             
          end if
          
          ccohort => ccohort%shorter
       end do

       if( sum(bc_out%rootfr_pa(fp,1:bc_in%nlevsoil)) > nearzero) then
          bc_out%rootfr_pa(fp,1:bc_in%nlevsoil) = &
               bc_out%rootfr_pa(fp,1:bc_in%nlevsoil) / &
               sum(bc_out%rootfr_pa(fp,1:bc_in%nlevsoil)) 
       end if
          
       ! RGK: These averages should switch to the new patch averaging methods
       !      when available.  Right now we are not doing any time averaging
       !      because it would be mixing the memory of patches, which
       !      would be arguably worse than just using the instantaneous value
       
       ! gC/m2/s
       bc_out%annavg_agnpp_pa(fp) = agnpp
       bc_out%annavg_bgnpp_pa(fp) = bgnpp
       ! gc/m2/yr
       bc_out%annsum_npp_pa(fp) = (bgnpp+agnpp)*days_per_year*sec_per_day

       if(plant_area>nearzero) then
          bc_out%woody_frac_aere_pa(fp) = woody_area/plant_area
       end if
       
       cpatch => cpatch%younger
    end do

    return
  end subroutine PrepCH4BCs
  
  ! =====================================================================================

  subroutine PrepNutrientAquisitionBCs(csite, bc_in, bc_out)

    ! -----------------------------------------------------------------------------------
    ! This subroutine will generate the appropriate boundary condition output 
    ! structures, depending on:
    ! 1) Which soil-bgc competition method is active in the HLM
    ! 2) If nitrification/denitrification is turned on
    ! 3) Which competitor scaling type is used
    ! -----------------------------------------------------------------------------------

    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) :: csite
    type(bc_in_type), intent(in)      :: bc_in
    type(bc_out_type), intent(inout)  :: bc_out

    ! Locals
    integer                       :: icomp   ! competitor index
    integer                       :: j       ! soil layer index
    integer                       :: id      ! decomp index (might == j)
    integer                       :: pft     ! plant functional type
    type(ed_patch_type), pointer  :: cpatch  ! current patch pointer
    type(ed_cohort_type), pointer :: ccohort ! current cohort pointer
    real(r8) :: fnrt_c                       ! fine-root carbon [kg]
    real(r8) :: veg_rootc                    ! fine root carbon in each layer [g/m3]
    real(r8) :: decompmicc_layer             ! Microbial dedcomposer biomass for current layer

    real(r8), parameter :: decompmicc_lambda = 2.5_r8     ! Depth attenuation exponent for decomposer biomass
    real(r8), parameter :: decompmicc_zmax   = 7.0e-2_r8  ! Depth of maximum decomposer biomass


    ! Whether this is a trivial or coupled run,
    ! the following variables get initialized in the same way
    bc_out%veg_rootc(:,:) = 0._r8
    bc_out%ft_index(:)    = -1
    if(trim(hlm_nu_com).eq.'ECA')then
       bc_out%decompmicc(:)  = 0._r8
       bc_out%cn_scalar(:)   = 1._r8
       bc_out%cp_scalar(:)   = 1._r8
    end if

    if(fates_np_comp_scaling == trivial_np_comp_scaling) then
       if(trim(hlm_nu_com).eq.'RD')then
          bc_out%num_plant_comps = 1
          bc_out%ft_index(1)    = 1
          return
       end if
    end if

    ! For both the trivial case with ECA, and the coupled case
    ! we still need to calculate the root biomass and decompmicc
    ! arrays (the former for the latter when trivial). So we
    ! don't differentiate

    icomp = 0
    cpatch => csite%oldest_patch
    do while (associated(cpatch))
       ccohort => cpatch%tallest
       do while (associated(ccohort))

          if(fates_np_comp_scaling .eq. coupled_np_comp_scaling) then
             icomp = icomp+1
          else
             icomp = 1
          end if

          pft   = ccohort%pft
          bc_out%ft_index(icomp) = pft
           
          call set_root_fraction(csite%rootfrac_scr, pft, csite%zi_soil, &
               bc_in%max_rooting_depth_index_col )

          fnrt_c   = ccohort%prt%GetState(fnrt_organ, carbon12_element)

          ! Map the soil layers to the decomposition layers (which may be synonomous)
          ! veg_rootc in units:  [gC/m3] = [kgC/plant] * [plant/ha] * [ha/ 10k m2] * [1000 g / kg] * [1/m]

          do j = 1, bc_in%nlevdecomp
             id = bc_in%decomp_id(j)  ! Map from soil layer to decomp layer     
             veg_rootc = fnrt_c * ccohort%n * csite%rootfrac_scr(j) * AREA_INV * g_per_kg / csite%dz_soil(j)

             bc_out%veg_rootc(icomp,id) = bc_out%veg_rootc(icomp,id) + veg_rootc

             if(trim(hlm_nu_com).eq.'ECA')then

                ! We use a 3 parameter exponential attenuation function to estimate decomposer biomass
                ! The parameter EDPftvarcon_inst%decompmicc(pft) is the maximum amount found at depth
                ! decompmicc_zmax, and the profile attenuates with strength lambda

                decompmicc_layer = EDPftvarcon_inst%decompmicc(pft) * &
                     exp(-decompmicc_lambda*abs(csite%z_soil(j)-decompmicc_zmax))


                bc_out%decompmicc(id) = bc_out%decompmicc(id) + decompmicc_layer * veg_rootc
             end if

          end do
          ccohort => ccohort%shorter
       end do

       cpatch => cpatch%younger
    end do

    ! We calculate the decomposer microbial biomass by weighting with the
    ! root biomass. This is just the normalization step
    if(trim(hlm_nu_com).eq.'ECA')then
       do id = 1,bc_in%nlevdecomp
          bc_out%decompmicc(id) = bc_out%decompmicc(id) / &
               max(nearzero,sum(bc_out%veg_rootc(:,id),dim=1))
       end do
    end if

    if(fates_np_comp_scaling == coupled_np_comp_scaling) then
       bc_out%num_plant_comps = icomp
    else
       bc_out%num_plant_comps = 1
    end if
    
    return
  end subroutine PrepNutrientAquisitionBCs

  ! =====================================================================================

  subroutine FluxIntoLitterPools(csite, bc_in, bc_out)
    
    ! -----------------------------------------------------------------------------------
    ! Created by Charlie Koven and Rosie Fisher, 2014-2015
    ! take the flux out of the fragmenting litter pools and port into the decomposing 
    ! litter pools. 
    ! in this implementation, decomposing pools are assumed to be humus and non-flammable, 
    ! whereas fragmenting pools are assumed to be physically fragmenting but not 
    ! respiring. This is a simplification, but allows us to 
    !
    ! a) reconcile the need to track both chemical fractions (lignin, cellulose, labile) 
    !    and size fractions (trunk, branch, etc.)
    ! b) to impose a realistic delay on the surge of nutrients into the litter pools 
    !    when large CWD is added to the system via mortality
    !
    ! Because of the different subgrid structure, this subroutine includes the functionality
    ! that in the big-leaf BGC model, is calculated in SoilBiogeochemVerticalProfileMod
    !
    ! The ED code is resolved at a daily timestep, but all of the CN-BGC fluxes are passed 
    ! in as derivatives per second, and then accumulated in the CNStateUpdate routines.  
    ! One way of doing this is to pass back the CN fluxes per second, and keep them 
    ! constant for the whole day (making sure they are not overwritten.  This means that 
    ! the carbon gets passed back and forth between the photosynthesis code 
    ! (fast timestepping) to the ED code (slow timestepping), back to the BGC code 
    ! (fast timestepping).  This means that the state update for the litter pools and 
    ! for the CWD pools occurs at different timescales. 
    ! -----------------------------------------------------------------------------------

    
    use FatesInterfaceTypesMod, only : bc_in_type, bc_out_type
    use FatesInterfaceTypesMod, only : hlm_use_vertsoilc
    use FatesInterfaceTypesMod, only : hlm_numlevgrnd
    use FatesConstantsMod, only : itrue
    use FatesGlobals, only : endrun => fates_endrun
    use EDParamsMod , only : ED_val_cwd_flig, ED_val_cwd_fcel
   
    

    implicit none   

    ! !ARGUMENTS    
    type(ed_site_type) , intent(inout)         :: csite
    type(bc_in_type)   , intent(in)            :: bc_in
    type(bc_out_type)  , intent(inout),target  :: bc_out

    ! !LOCAL VARIABLES:
    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    real(r8), pointer              :: flux_cel_si(:)
    real(r8), pointer              :: flux_lab_si(:)
    real(r8), pointer              :: flux_lig_si(:)
    real(r8), pointer              :: efflux_ptr      ! Points to the current
                                                      ! element's root efflux                                         
    type(litter_type), pointer     :: litt
     
    real(r8) :: surface_prof(1:hlm_numlevgrnd) ! this array is used to distribute
                                               ! fragmented litter on the surface
                                               ! into the soil/decomposition
                                               ! layers. It exponentially decays
    real(r8) :: surface_prof_tot ! normalizes the surface_prof array
    integer  :: nlev_eff_soil    ! number of effective soil layers
    integer  :: nlev_eff_decomp  ! number of effective decomp layers
    real(r8) :: area_frac        ! fraction of site's area of current patch
    real(r8) :: z_decomp         ! Used for calculating depth midpoints of decomp layers
    integer  :: s                ! Site index
    integer  :: el               ! Element index (C,N,P,etc)
    integer  :: j                ! Soil layer index
    integer  :: id               ! Decomposition layer index
    integer  :: ic               ! CWD type index
    integer  :: ipft             ! PFT index

    ! NOTE(rgk, 201705) this parameter was brought over from SoilBiogeochemVerticalProfile
    ! how steep profile is for surface components (1/ e_folding depth) (1/m) 
    real(r8),  parameter :: surfprof_exp  = 10.

    ! This is the number of effective soil layers to transfer from
    nlev_eff_soil   = max(bc_in%max_rooting_depth_index_col, 1)
    
    ! The decomposition layers are most likely the exact same layers
    ! as the soil layers (same depths also), unless it is a simplified
    ! single layer case, where nlevdecomp = 1
    
    nlev_eff_decomp = min(bc_in%nlevdecomp,nlev_eff_soil)
    
    ! define a single shallow surface profile for surface additions 
    ! (leaves, stems, and N deposition). This sends the above ground
    ! mass into the soil pools using an exponential depth decay function.
    ! Since it is sending an absolute mass [kg] into variable layer
    ! widths, we multiply the profile by the layer width, so that
    ! wider layers get proportionally more.  After the masses
    ! are sent, each layer will normalize by depth.
    
    surface_prof(:) = 0._r8
    z_decomp = 0._r8
    do id = 1,nlev_eff_decomp
       z_decomp = z_decomp+0.5*bc_in%dz_decomp_sisl(id)
       surface_prof(id) = exp(-surfprof_exp * z_decomp) *  bc_in%dz_decomp_sisl(id)
       z_decomp = z_decomp+0.5*bc_in%dz_decomp_sisl(id)
    end do
    surface_prof_tot = sum(surface_prof)
    do id = 1,nlev_eff_decomp
       surface_prof(id) = surface_prof(id)/surface_prof_tot
    end do

    ! Loop over the different elements. 
    do el = 1, num_elements
       
       ! Zero out the boundary flux arrays
       ! Make a pointer to the cellulose, labile and lignan
       ! flux partitions.
       
       select case (element_list(el))
       case (carbon12_element)
          bc_out%litt_flux_cel_c_si(:) = 0.0_r8
          bc_out%litt_flux_lig_c_si(:) = 0.0_r8
          bc_out%litt_flux_lab_c_si(:) = 0.0_r8
          flux_cel_si => bc_out%litt_flux_cel_c_si(:)
          flux_lab_si => bc_out%litt_flux_lab_c_si(:)
          flux_lig_si => bc_out%litt_flux_lig_c_si(:)
          
       case (nitrogen_element) 
          bc_out%litt_flux_cel_n_si(:) = 0._r8
          bc_out%litt_flux_lig_n_si(:) = 0._r8
          bc_out%litt_flux_lab_n_si(:) = 0._r8
          flux_cel_si => bc_out%litt_flux_cel_n_si(:)
          flux_lab_si => bc_out%litt_flux_lab_n_si(:)
          flux_lig_si => bc_out%litt_flux_lig_n_si(:)
          
       case (phosphorus_element)
          bc_out%litt_flux_cel_p_si(:) = 0._r8
          bc_out%litt_flux_lig_p_si(:) = 0._r8
          bc_out%litt_flux_lab_p_si(:) = 0._r8
          flux_cel_si => bc_out%litt_flux_cel_p_si(:)
          flux_lab_si => bc_out%litt_flux_lab_p_si(:)
          flux_lig_si => bc_out%litt_flux_lig_p_si(:)
          
       end select

       currentPatch => csite%oldest_patch
       do while (associated(currentPatch))

          ! If there is any efflux (from stores overflowing)
          ! than pass that to the labile litter pool
          
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))
             if(.not.currentCohort%isnew)then
                if(element_list(el).eq.carbon12_element) then
                   efflux_ptr => currentCohort%daily_c_efflux
                elseif(element_list(el).eq.nitrogen_element) then
                   efflux_ptr => currentCohort%daily_n_efflux
                elseif(element_list(el).eq.phosphorus_element) then
                   efflux_ptr => currentCohort%daily_p_efflux
                end if

                call set_root_fraction(csite%rootfrac_scr, currentCohort%pft, csite%zi_soil, &
                     bc_in%max_rooting_depth_index_col )
                
                ! Unit conversion
                ! kg/plant/day * plant/ha * ha/m2 -> kg/m2/day
                ! Also, lets efflux out through the roots
                
                !do id = 1,nlev_eff_decomp
                !   flux_lab_si(id) = flux_lab_si(id) + &
                !        efflux_ptr * currentCohort%n* AREA_INV * surface_prof(id)
                !end do

                do id = 1,nlev_eff_soil
                   flux_lab_si(id) = flux_lab_si(id) + &
                        efflux_ptr * currentCohort%n* AREA_INV *csite%rootfrac_scr(id)
                end do
                
             end if
             currentCohort => currentCohort%shorter
          end do


          ! Set a pointer to the litter object
          ! for the current element on the current
          ! patch
          litt       => currentPatch%litter(el)
          area_frac  = currentPatch%area/area
          
          do ic = 1, ncwd

             do id = 1,nlev_eff_decomp
                flux_cel_si(id) = flux_cel_si(id) + &
                     litt%ag_cwd_frag(ic) * ED_val_cwd_fcel * area_frac * surface_prof(id)

                flux_lig_si(id) = flux_lig_si(id) + & 
                     litt%ag_cwd_frag(ic) * ED_val_cwd_flig * area_frac * surface_prof(id)
             end do

             do j = 1, nlev_eff_soil

                id = bc_in%decomp_id(j)  ! Map from soil layer to decomp layer

                flux_cel_si(id) = flux_cel_si(id) + &
                     litt%bg_cwd_frag(ic,j) * ED_val_cwd_fcel * area_frac

                flux_lig_si(id) = flux_lig_si(id) + &
                     litt%bg_cwd_frag(ic,j) * ED_val_cwd_flig * area_frac

             end do
          end do

          ! leaf and fine root fragmentation fluxes

          do id = 1,nlev_eff_decomp

             flux_lab_si(id) = flux_lab_si(id) + &
                  litt%leaf_fines_frag(ilabile) * area_frac* surface_prof(id)

             flux_cel_si(id) = flux_cel_si(id) + &
                  litt%leaf_fines_frag(icellulose) * area_frac* surface_prof(id)

             flux_lig_si(id) = flux_lig_si(id) + &
                  litt%leaf_fines_frag(ilignin) * area_frac* surface_prof(id)

          end do


          ! decaying seeds from the litter pool
          do ipft = 1,numpft
             do id = 1,nlev_eff_decomp

                flux_lab_si(id) = flux_lab_si(id) + &
                     (litt%seed_decay(ipft) + litt%seed_germ_decay(ipft)) * &
                     EDPftvarcon_inst%lf_flab(ipft) * area_frac* surface_prof(id)

                flux_cel_si(id) = flux_cel_si(id) + &
                     (litt%seed_decay(ipft) + litt%seed_germ_decay(ipft)) * &
                     EDPftvarcon_inst%lf_fcel(ipft) * area_frac* surface_prof(id)

                flux_lig_si(id) = flux_lig_si(id) + &
                     (litt%seed_decay(ipft) + litt%seed_germ_decay(ipft)) * &
                     EDPftvarcon_inst%lf_flig(ipft) * area_frac* surface_prof(id)
             end do
          end do


          do j = 1, nlev_eff_soil

             id = bc_in%decomp_id(j)
             flux_lab_si(id) = flux_lab_si(id) + &
                  litt%root_fines_frag(ilabile,j) * area_frac
             flux_cel_si(id) = flux_cel_si(id) + &
                  litt%root_fines_frag(icellulose,j) * area_frac
             flux_lig_si(id) = flux_lig_si(id) + &
                  litt%root_fines_frag(ilignin,j) * area_frac
          enddo

          currentPatch => currentPatch%younger
       end do

       ! Normalize all masses over the decomposition layer's depth
       ! Convert from kg/m2/day -> g/m3/s

       do id = 1,nlev_eff_decomp
          flux_cel_si(id) = days_per_sec * g_per_kg * &
               flux_cel_si(id) / bc_in%dz_decomp_sisl(id)
          flux_lig_si(id) = days_per_sec * g_per_kg * &
               flux_lig_si(id) / bc_in%dz_decomp_sisl(id)
          flux_lab_si(id) = days_per_sec * g_per_kg * &
               flux_lab_si(id) / bc_in%dz_decomp_sisl(id)
       end do

    end do  ! do elements

    return
  end subroutine FluxIntoLitterPools




end module FatesSoilBGCFluxMod
