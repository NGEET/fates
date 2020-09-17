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
  use FatesAllometryMod, only : set_root_fraction
  use FatesAllometryMod , only : h_allom
  use FatesAllometryMod , only : h2d_allom
  use FatesAllometryMod , only : bagw_allom
  use FatesAllometryMod , only : bsap_allom
  use FatesAllometryMod , only : bleaf
  use FatesAllometryMod , only : bfineroot
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
  use FatesConstantsMod , only : prescribed_p_uptake
  use FatesConstantsMod , only : prescribed_n_uptake
  use FatesConstantsMod , only : coupled_p_uptake
  use FatesConstantsMod , only : coupled_n_uptake
  use FatesConstantsMod, only    : days_per_sec
  use FatesConstantsMod, only    : g_per_kg
  use FatesConstantsMod, only    : kg_per_g
  use FatesConstantsMod, only    : fates_np_comp_scaling
  use FatesConstantsMod, only    : cohort_np_comp_scaling
  use FatesConstantsMod, only    : pft_np_comp_scaling
  use FatesConstantsMod, only    : rsnbl_math_prec
  use FatesLitterMod,        only : litter_type
  use FatesLitterMod    , only : ncwd
  use FatesLitterMod    , only : ndcmpy
  use FatesLitterMod    , only : ilabile
  use FatesLitterMod    , only : ilignin
  use FatesLitterMod    , only : icellulose
  use PRTParametersMod , only    : prt_params
  use EDPftvarcon      , only    : EDPftvarcon_inst
  
  implicit none
  private
  
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

    real(r8), parameter :: smth_fac = 0.8_r8         ! Smoothing factor for updating
                                                     ! demand.
    real(r8), parameter :: init_demand_frac = 0.1_r8 ! Newly recruited plants
                                                     ! will specify a demand
                                                     ! based on their total nutrient
                                                     ! because they have not history
                                                     ! of need yet

    
    
    pft = ccohort%pft
    dbh = ccohort%dbh


    ! If the cohort has not experienced a day of integration
    ! (and thus any allocation yet), we specify demand
    ! based purely on a fraction of its starting nutrient content
    if(ccohort%isnew) then

       if(element_id.eq.nitrogen_element) then
          plant_max_x = & 
               ccohort%prt%GetState(leaf_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,leaf_organ) + & 
               ccohort%prt%GetState(fnrt_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,fnrt_organ) + & 
               ccohort%prt%GetState(store_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,store_organ) + & 
               ccohort%prt%GetState(sapw_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,sapw_organ) + & 
               ccohort%prt%GetState(struct_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,struct_organ) + &
               ccohort%prt%GetState(repro_organ, carbon12_element)*prt_params%nitr_stoich_p2(pft,repro_organ)
       
       elseif(element_id.eq.phosphorus_element) then
          plant_max_x = &
               ccohort%prt%GetState(leaf_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,leaf_organ) + & 
               ccohort%prt%GetState(fnrt_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,fnrt_organ) + & 
               ccohort%prt%GetState(store_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,store_organ) + & 
               ccohort%prt%GetState(sapw_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,sapw_organ) + & 
               ccohort%prt%GetState(struct_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,struct_organ) + &
               ccohort%prt%GetState(repro_organ, carbon12_element)*prt_params%phos_stoich_p2(pft,repro_organ)
          
       end if

       plant_demand = init_demand_frac*plant_max_x
       return
    end if
       

    ! If the plant is not a newly recruited plant
    ! We use other methods of specifying nutrient demand
    ! -----------------------------------------------------------------------------------

    if(element_id.eq.nitrogen_element) then

       plant_demand = smth_fac*ccohort%daily_n_demand + (1._r8-smth_fac)*ccohort%daily_n_need2
       
    elseif(element_id.eq.phosphorus_element) then
 
       plant_demand = smth_fac*ccohort%daily_p_demand + (1._r8-smth_fac)*ccohort%daily_p_need2
       
    end if

  end function GetPlantDemand

  ! =====================================================================================
  
  
  subroutine UnPackNutrientAquisitionBCs(sites, bc_in)

    ! -----------------------------------------------------------------------------------
    ! The purpose of this routine is to recieve the nutrient uptake flux
    ! boundary conditions, and parse those fluxes out to the cohorts.
    !
    ! This routine should be called before FATES dynamics, particularly before
    ! any of the PARTEH code is called.  It is assumed that these utpake fluxes
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
    real(r8), allocatable :: fnrt_c_pft(:)            ! total mass of root for each PFT [kgC]


    nsites = size(sites,dim=1)

    
    ! Zero the uptake rates
    do s = 1, nsites
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             ccohort%daily_n_uptake = 0._r8
             ccohort%daily_p_uptake = 0._r8
             ccohort => ccohort%shorter
          end do
          cpatch => cpatch%younger
       end do
       
    end do
    
    ! We can exit if this is a c-only simulation
    if(hlm_parteh_mode.eq.prt_carbon_allom_hyp) then
       ! These can now be zero'd
       do s = 1, nsites
          bc_in(s)%plant_n_uptake_flux(:,:) = 0._r8
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
                ccohort%daily_n_uptake = &
                     EDPftvarcon_inst%prescribed_nuptake(pft) * &
                     GetPlantDemand(ccohort,nitrogen_element)
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
                ccohort%daily_p_uptake = &
                     EDPftvarcon_inst%prescribed_puptake(pft) * &
                     GetPlantDemand(ccohort,phosphorus_element)
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
       end if
       
       
       ! If nutrient competition is sent to the BGC model as PFTs
       ! and not as individual cohorts, we need to unravel the input
       ! boundary condition and send to cohort.  We do this downscaling
       ! by finding each cohort's fraction of total fine-root for the group
       
       if(n_uptake_mode.eq.coupled_n_uptake .or. p_uptake_mode.eq.coupled_p_uptake)then

          ! Note there are two scaling methods.  Either competition for
          ! N and/or P was performed by cohorts acting individually
          ! (cohort_np_comp_scaling) , or as PFTs (pft_np_comp_scaling)
          ! If we opt for the latter, then we assume that the nutrient
          ! uptake share of the cohort, matches the fraction of root
          ! mass it contributes to the group (PFT).
          
          if(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
             
             ! *Currently, all cohorts in a PFT have the same root
             ! fraction, so all we have to to is find its total mass fraction.
             
             fnrt_c_pft(:) = 0._r8
             cpatch => sites(s)%oldest_patch
             do while (associated(cpatch))
                ccohort => cpatch%tallest
                do while (associated(ccohort))
                   pft   = ccohort%pft
                   fnrt_c_pft(pft) = fnrt_c_pft(pft) + &
                        ccohort%prt%GetState(fnrt_organ, all_carbon_elements)*ccohort%n
                   ccohort => ccohort%shorter
                end do
                cpatch => cpatch%younger
             end do
             
          end if  ! end if(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
          
          ! --------------------------------------------------------------------------------
          ! Now that we have the arrays ready for downscaling (if needed)
          ! loop through all cohorts and acquire nutrient
          ! --------------------------------------------------------------------------------

          if(n_uptake_mode.eq.coupled_n_uptake) then
             icomp = 0
             cpatch => sites(s)%oldest_patch
             do while (associated(cpatch))
                ccohort => cpatch%tallest
                do while (associated(ccohort))
                   pft   = ccohort%pft
                   if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                      icomp = icomp+1

                      ! N Uptake:  Convert g/m2/day -> kg/plant/day
                      ccohort%daily_n_uptake = ccohort%daily_n_uptake + &
                           sum(bc_in(s)%plant_n_uptake_flux(icomp,:))*kg_per_g*AREA/ccohort%n

                   else
                      icomp = pft
                      ! Total fine-root carbon of the cohort [kgC/ha]
                      fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)*ccohort%n
                      
                      ! Loop through soil layers, add up the uptake this cohort gets from each layer
                      do id = 1,bc_in(s)%nlevdecomp
                         ccohort%daily_n_uptake = ccohort%daily_n_uptake + & 
                              bc_in(s)%plant_n_uptake_flux(icomp,id) * &
                              (fnrt_c/fnrt_c_pft(pft))*kg_per_g*AREA/ccohort%n
                      end do
                   end if
                   ccohort => ccohort%shorter
                end do
                cpatch => cpatch%younger
             end do
          end if

          if(p_uptake_mode.eq.coupled_p_uptake) then

             icomp = 0
             cpatch => sites(s)%oldest_patch
             do while (associated(cpatch))
                ccohort => cpatch%tallest
                do while (associated(ccohort))
                   pft   = ccohort%pft
                   if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                      icomp = icomp+1

                      ! P Uptake:  Convert g/m2/day -> kg/plant/day
                      ccohort%daily_p_uptake = ccohort%daily_p_uptake + &
                           sum(bc_in(s)%plant_p_uptake_flux(icomp,:))*kg_per_g*AREA/ccohort%n
                      
                   else
                      icomp = pft
                      ! Total fine-root carbon of the cohort [kgC/ha]
                      fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)*ccohort%n

                      ! Loop through soil layers, add up the uptake this cohort gets from each layer
                      do id = 1,bc_in(s)%nlevdecomp
                         ccohort%daily_p_uptake = ccohort%daily_p_uptake + & 
                              bc_in(s)%plant_p_uptake_flux(icomp,id) * &
                              (fnrt_c/fnrt_c_pft(pft))*kg_per_g*AREA/ccohort%n
                      end do
                   end if
                   ccohort => ccohort%shorter
                end do
                cpatch => cpatch%younger
             end do
          end if
          
          ! Free the temprorary arrays (if used)
          if(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
             deallocate(fnrt_c_pft)
          end if
          
       end if
       
       ! These can now be zero'd
       bc_in(s)%plant_n_uptake_flux(:,:) = 0._r8
       bc_in(s)%plant_p_uptake_flux(:,:) = 0._r8

    end do
    return
  end subroutine UnPackNutrientAquisitionBCs

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
    integer                       :: icomp         ! competitor index
    integer                       :: j             ! soil layer index
    integer                       :: id            ! decomp index (might == j)
    integer                       :: pft           ! plant functional type
    integer                       :: nlev_eff_soil ! number of active soil layers
    type(ed_patch_type), pointer  :: cpatch        ! current patch pointer
    type(ed_cohort_type), pointer :: ccohort       ! current cohort pointer
    real(r8) :: fnrt_c                             ! fine-root carbon [kg]

    ! Note that "leaf" in this context may contain the storage along with it

    real(r8) :: leaf_c                             ! leaf carbon [kg]
    real(r8) :: leaf_p                             ! leaf phosphorus [kg]
    real(r8) :: leaf_n                             ! leaf nitrogen [kg]
    real(r8) :: leaf_store_n
    real(r8) :: leaf_store_p
    real(r8) :: veg_rootc                          ! fine root carbon in each layer [g/m3]
    real(r8) :: dbh                                ! dbh (cm)
    real(r8) :: nc_leaf_ideal                      ! ideal leaf C:N ratio [kg/kg]
    real(r8) :: pc_leaf_ideal                      ! ideal leaf C:P ratio [kg/kg]
    real(r8) :: nc_leaf_actual                     ! actual leaf C:N ratio [kg/kg]
    real(r8) :: pc_leaf_actual                     ! actual leaf C:P ratio [kg/kg]
    real(r8) :: nc_leaf_min                        ! minimum leaf C:N ratio [kg/kg]
    real(r8) :: pc_leaf_min                        ! minimum leaf C:P ratio [kg/kg]
    real(r8) :: target_leaf_c                      ! maximum leaf C for this dbh [kg]
    real(r8) :: target_store_c                     ! maximum store C for this dbh [kg]
    real(r8) :: npp_n_demand                       ! Nitrogen needed to keep up with NPP  [kgN]
    real(r8) :: npp_p_demand                       ! Phosphorus needed to keep up with NPP [kgP]
    real(r8) :: deficit_n_demand                   ! Nitrogen needed to get stoich back to
                                                   ! optimal [kgN]
    real(r8) :: deficit_p_demand                   ! Phosphorus needed to get stoich back to
                                                   ! optimal [kgP]

    real(r8) :: nc_actual                          ! Actual N:C ratio of plant 
    real(r8) :: nc_min                             ! Minimum allowable N:C ratio to build tissue
    real(r8) :: nc_ideal                           ! Plant's ideal N:C ratio
    real(r8) :: cn_actual                          ! Actual C:N ratio of plant
    real(r8) :: cn_ideal                           ! Ideal C:N ratio of plant

    real(r8) :: pc_actual                          ! Actual P:C ratio of plant 
    real(r8) :: pc_min                             ! Minimum allowable P:C ratio to build tissue
    real(r8) :: pc_ideal                           ! Plant's ideal P:C ratio
    real(r8) :: cp_actual                          ! Actual C:P ratio of plant
    real(r8) :: cp_ideal                           ! Ideal C:P ratio of plant


    real(r8), parameter :: cn_stoich_var=0.2    ! variability of CN ratio (ECA)
    real(r8), parameter :: cp_stoich_var=0.4    ! variability of CP ratio (ECA)

    ! ECA: We are still testing different plant "neediness" functions, thus
    ! two methods
    !
    integer, parameter :: cnp_scalar_method1 = 1     ! RK's modification for FATES
    integer, parameter :: cnp_scalar_method2 = 2     ! Q. Zhu's original formulation
    integer, parameter :: cnp_scalar_method3 = 3     ! Force to be 1.0 all the time
    integer, parameter :: cnp_scalar_method  = cnp_scalar_method3


    real(r8) :: comp_per_pft(numpft) ! Competitors per PFT, used for averaging


    ! Run the trivial case where we do not have a nutrient model
    ! running in fates, send zero demands to the BGC model
    if((hlm_parteh_mode.ne.prt_cnp_flex_allom_hyp)) then! .or. &
!         (n_uptake_mode.eq.prescribed_n_uptake .and. p_uptake_mode.eq.prescribed_p_uptake)) then
       bc_out%n_plant_comps  = 1
       if(trim(hlm_nu_com).eq.'ECA')then
          bc_out%ft_index(1)    = 1
          bc_out%veg_rootc(1,:) = 0._r8
          bc_out%cn_scalar(1)   = 0._r8
          bc_out%cp_scalar(1)   = 0._r8
          bc_out%decompmicc(1)  = 0._r8
       elseif(trim(hlm_nu_com).eq.'RD') then
          bc_out%n_demand(1) = 0._r8
          bc_out%p_demand(1) = 0._r8
       end if
       return
    end if
    
    ! This is the number of effective soil layers to transfer from
    nlev_eff_soil   = max(bc_in%max_rooting_depth_index_col, 1)

    ! It is possible, through either symbiotic fixation, or through
    ! exudation, that FATES plants may provide a nutrient source term
    ! in the mineralized nitrogen (nh4) pool, same for phosphorous
    ! (CURRENTLY WE DO NOT HAVE THOSE ALROGITHMS)

    ! bc_out(s)%source_n(:) = 0._r8
    ! bc_out(s)%source_p(:) = 0._r8

    
    ! If we are using either un-coupled nutrient uptake, or RD competition
    ! in the soil BGC model, we must update the plant's demand
    ! -----------------------------------------------------------------------------------
    
    if( (trim(hlm_nu_com).eq.'RD') .or. &
        (n_uptake_mode .eq. prescribed_n_uptake) ) then

       cpatch => csite%oldest_patch
       do while (associated(cpatch))
          ccohort => cpatch%tallest
          do while (associated(ccohort))

             ccohort%daily_n_demand = GetPlantDemand(ccohort,nitrogen_element)
             ccohort%daily_p_demand = GetPlantDemand(ccohort,phosphorus_element)

             ccohort => ccohort%shorter
          end do
          
          cpatch => cpatch%younger
       end do
    end if

    
    ! ECA Specific Parameters
    ! --------------------------------------------------------------------------------
    if(trim(hlm_nu_com).eq.'ECA')then

       bc_out%veg_rootc(:,:) = 0._r8  ! Zero this, it will be incremented
       bc_out%cn_scalar(:)   = 0._r8
       bc_out%cp_scalar(:)   = 0._r8
       bc_out%decompmicc(:)  = 0._r8
       bc_out%ft_index(:)    = -1
       
       ! Loop over all patches and sum up the seed input for each PFT
       icomp = 0
       comp_per_pft(:) = 0     ! This counts how many competitors per
       ! pft, used for averaging
       cpatch => csite%oldest_patch
       do while (associated(cpatch))
          
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             
             pft   = ccohort%pft
             
             if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                icomp = icomp+1
             else
                icomp = pft
                comp_per_pft(pft) = comp_per_pft(pft) + 1
             end if
             
             call set_root_fraction(csite%rootfrac_scr, pft, csite%zi_soil)
             
             fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)
             
             ! Map the soil layers to the decomposition layers
             ! (which may be synonomous)
             ! veg_rootc in units:  [g/m3] = [kgC/plant] * [plant/ha] * [ha/ 10k m2] * [1000 g / kg] * [1/m]
             
             do j = 1, nlev_eff_soil
                id = bc_in%decomp_id(j)  ! Map from soil layer to decomp layer     
                veg_rootc = fnrt_c * ccohort%n * csite%rootfrac_scr(j) * AREA_INV * g_per_kg / csite%dz_soil(j)
                bc_out%veg_rootc(icomp,id) = bc_out%veg_rootc(icomp,id) + veg_rootc
                bc_out%decompmicc(id) = bc_out%decompmicc(id) + &
                     EDPftvarcon_inst%decompmicc(pft) * veg_rootc
             end do
             
             bc_out%ft_index(icomp) = pft

             ccohort => ccohort%shorter
          end do
          
          cpatch => cpatch%younger
       end do
       
       if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
          bc_out%n_plant_comps = icomp
       else
          bc_out%n_plant_comps = numpft
       end if

       ! We calculate the decomposer microbial biomass by weighting with the
       ! root biomass. This is just the normalization step
       do id = 1,bc_in%nlevdecomp
          bc_out%decompmicc(id) = bc_out%decompmicc(id) / &
               max(nearzero,sum(bc_out%veg_rootc(:,id),dim=1))
       end do

       ! ECA, in coupled mode affinity is diagnosed

       if(n_uptake_mode.eq.coupled_n_uptake) then
          icomp = 0
          cpatch => csite%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))

                pft   = ccohort%pft
                if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                   icomp = icomp+1
                else
                   icomp = pft
                end if
                
!                leaf_c = max(rsnbl_math_prec,ccohort%prt%GetState(leaf_organ, all_carbon_elements) + &
!                                             ccohort%prt%GetState(store_organ, all_carbon_elements))

                ! Target leaf biomass according to allometry and trimming
                call bleaf(ccohort%dbh,pft,ccohort%canopy_trim,target_leaf_c)
                call bstore_allom(ccohort%dbh,pft,ccohort%canopy_trim,target_store_c)
                
                leaf_store_n = max(rsnbl_math_prec,ccohort%prt%GetState(leaf_organ, nitrogen_element) + & 
                                                   ccohort%prt%GetState(store_organ, nitrogen_element))
                
                ! Calculate the ideal CN ratio for leaves and storage organs
                nc_ideal = ((target_leaf_c*prt_params%nitr_stoich_p2(pft,leaf_organ)) + &
                     (target_store_c*prt_params%nitr_stoich_p2(pft,store_organ))) / & 
                     (target_leaf_c+target_store_c)
                nc_min =  ((target_leaf_c*prt_params%nitr_stoich_p1(pft,leaf_organ)) + &
                     (target_store_c*prt_params%nitr_stoich_p1(pft,store_organ))) / & 
                     (target_leaf_c+target_store_c)
                
                nc_actual = max(leaf_store_n/(target_leaf_c+target_store_c),rsnbl_math_prec)

                select case(cnp_scalar_method)
                case(cnp_scalar_method1)
                   
                    !x0 = 0.5*(nc_ideal+nc_min)
                    
                    ! Fit the logistic shape parameter so that 95%tile of
                    ! nutrient concentration matches 95%tile of scalar
                    ! 0.95 = 1._r8/(1._r8 + exp(-logi_k*(  0.95*(nc_ideal-x0) )))
                    ! logi_k = -log(1._r8-0.95/0.95)/ (  0.95*(nc_ideal-x0) )
                    
                    !bc_out%cn_scalar(icomp) = 1._r8/(1._r8 + exp(-logi_k*(nc_actual-x0)))
                    
                    bc_out%cn_scalar(icomp) = bc_out%cn_scalar(icomp) + & 
                          min(1._r8,max(0._r8, & 
                          (nc_ideal - nc_actual + cn_stoich_var*nc_min) / & 
                          (nc_ideal - nc_min + cn_stoich_var*nc_min)))
                    
                case(cnp_scalar_method2)
                    cn_ideal = 1._r8/nc_ideal
                    cn_actual = 1._r8/nc_actual
                    bc_out%cn_scalar(icomp) = bc_out%cn_scalar(icomp) + & 
                          min(1._r8,max(0._r8, & 
                          (cn_actual - cn_ideal*(1._r8-cn_stoich_var))/(cn_ideal*cn_stoich_var)))
                    
                   
                case(cnp_scalar_method3)
                    ! Force cn_scalar to be maxed
                    bc_out%cn_scalar(icomp) = bc_out%cn_scalar(icomp) + 1.0_r8
                case default
                    write(fates_log(), *) 'undefined cnp_scalar mode'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                end select
                
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do

          ! Normalize the sum to a mean, if this is a PFT scale
          ! boundary flux
          if(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
             do icomp = 1, numpft
                bc_out%cn_scalar(icomp) = bc_out%cn_scalar(icomp)/real(comp_per_pft(icomp),r8)
             end do
          end if
          
       end if
       
       if(p_uptake_mode.eq.coupled_p_uptake) then

          
          ! ECA, in coupled mode affinity for P is diagnosed
          icomp = 0
          cpatch => csite%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))

                pft   = ccohort%pft
                if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                   icomp = icomp+1
                else
                   icomp = pft
                end if

                ! Target leaf biomass according to allometry and trimming
                call bleaf(ccohort%dbh,pft,ccohort%canopy_trim,target_leaf_c)
                call bstore_allom(ccohort%dbh,pft,ccohort%canopy_trim,target_store_c)
                
                leaf_store_p = max(rsnbl_math_prec,ccohort%prt%GetState(leaf_organ, phosphorus_element) + & 
                                                   ccohort%prt%GetState(store_organ, phosphorus_element))

                
                ! Calculate the ideal CN ratio for leaves and storage organs
                pc_ideal = ((target_leaf_c*prt_params%phos_stoich_p2(pft,leaf_organ)) + &
                     (target_store_c*prt_params%phos_stoich_p2(pft,store_organ))) / & 
                     (target_leaf_c+target_store_c)
                pc_min = ((target_leaf_c*prt_params%phos_stoich_p1(pft,leaf_organ)) + &
                     (target_store_c*prt_params%phos_stoich_p1(pft,store_organ))) / & 
                     (target_leaf_c+target_store_c)
                
                pc_actual = max(leaf_store_p/(target_leaf_c+target_store_c),rsnbl_math_prec)

                select case(cnp_scalar_method)
                case(cnp_scalar_method1)
                    bc_out%cp_scalar(icomp) = bc_out%cp_scalar(icomp) + & 
                          min(1._r8,max(0._r8, & 
                          (pc_ideal - pc_actual + cp_stoich_var*pc_min) / & 
                          (pc_ideal - pc_min + cp_stoich_var*pc_min)))
                    
                case(cnp_scalar_method2)
                    cp_actual = 1._r8/pc_actual
                    cp_ideal  = 1._r8/pc_ideal
                    bc_out%cp_scalar(icomp) = bc_out%cp_scalar(icomp) + & 
                          min(1._r8,max(0._r8, & 
                          (cp_actual - cp_ideal*(1._r8-cn_stoich_var))/(cp_ideal*cn_stoich_var)))
                    
                case(cnp_scalar_method3)
                    ! Force cn_scalar to be maxed
                    bc_out%cp_scalar(icomp) = bc_out%cp_scalar(icomp) + 1.0_r8
                    
                case default
                    write(fates_log(), *) 'undefined cnp_scalar mode'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                end select

                
                ccohort => ccohort%shorter
             end do
             
             cpatch => cpatch%younger
          end do
          
          if(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
             do icomp = 1, numpft
                bc_out%cp_scalar(icomp) = bc_out%cp_scalar(icomp)/real(comp_per_pft(icomp),r8)
             end do
          end if
          
       end if
       
    elseif(trim(hlm_nu_com).eq.'RD') then

       ! RD Specific Parameters
       ! --------------------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! Demand comes from two parts. (1) Some demand seeks to fullfill the N 
       ! or P requirements needed for the construction costs that match NPP. 
       ! We have to predict the next day based on the current day's NPP.  (2) 
       ! If a plant has flexible stoichiometry, it may currently have an N or 
       ! P concentration lower than optimal, and therefore demand seeks to 
       ! replace this deficit.
       ! --------------------------------------------------------------------
       ! [gX/m2/s]  convert [kgX/plant/day] * [plant/ha] * 
       !                                      [ha/10000 m2] * [1000 g/kg] * [1 day /86400 sec]
       
       bc_out%n_demand(:) = 0._r8
       bc_out%p_demand(:) = 0._r8
       
       if(n_uptake_mode.eq.coupled_n_uptake) then
          icomp = 0
          cpatch => csite%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                pft   = ccohort%pft
                dbh   = ccohort%dbh
                if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                   icomp = icomp+1
                else
                   icomp = pft
                end if
                bc_out%n_demand(icomp) = bc_out%n_demand(icomp) + &
                     ccohort%daily_n_demand*ccohort%n*AREA_INV*g_per_kg*days_per_sec
                     
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
       end if
       
       if(p_uptake_mode.eq.coupled_p_uptake) then
          icomp = 0
          cpatch => csite%oldest_patch
          do while (associated(cpatch))
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                pft   = ccohort%pft
                dbh   = ccohort%dbh
                if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
                   icomp = icomp+1
                else
                   icomp = pft
                end if
                bc_out%p_demand(icomp) = bc_out%p_demand(icomp) + & 
                     ccohort%daily_p_demand*ccohort%n*AREA_INV*g_per_kg*days_per_sec
                ccohort => ccohort%shorter
             end do
             cpatch => cpatch%younger
          end do
       end if

       if( (n_uptake_mode.eq.coupled_p_uptake) .or. &
           (p_uptake_mode.eq.coupled_p_uptake)) then
          if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
             bc_out%n_plant_comps = icomp
          else
             bc_out%n_plant_comps = numpft
          end if
          
       else
          bc_out%n_plant_comps = 1
          
       end if
       
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

    
    use FatesConstantsMod, only : sec_per_day
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
    integer  :: ft               ! PFT number
    integer  :: nlev_eff_soil    ! number of effective soil layers
    integer  :: nlev_eff_decomp  ! number of effective decomp layers
    real(r8) :: area_frac        ! fraction of site's area of current patch
    real(r8) :: z_decomp         ! Used for calculating depth midpoints of decomp layers
    integer  :: s                ! Site index
    integer  :: el               ! Element index (C,N,P,etc)
    integer  :: j                ! Soil layer index
    integer  :: id               ! Decomposition layer index
    integer  :: ic               ! CWD type index

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
          bc_out%litt_flux_cel_c_si(:) = 0._r8
          bc_out%litt_flux_lig_c_si(:) = 0._r8
          bc_out%litt_flux_lab_c_si(:) = 0._r8
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

          ! If we have prescribed boundary conditions
          ! we do not take N out of the BGC model's
          ! stores, so nor do we send any back
          if(n_uptake_mode.eq.prescribed_n_uptake) cycle
          
       case (phosphorus_element)
          bc_out%litt_flux_cel_p_si(:) = 0._r8
          bc_out%litt_flux_lig_p_si(:) = 0._r8
          bc_out%litt_flux_lab_p_si(:) = 0._r8
          flux_cel_si => bc_out%litt_flux_cel_p_si(:)
          flux_lab_si => bc_out%litt_flux_lab_p_si(:)
          flux_lig_si => bc_out%litt_flux_lig_p_si(:)

          ! If we have prescribed boundary conditions
          ! we do not take N out of the BGC model's
          ! stores, so nor do we send any back
          if(p_uptake_mode.eq.prescribed_p_uptake) cycle
          
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
                do id = 1,nlev_eff_decomp
                   flux_lab_si(id) = flux_lab_si(id) + &
                        efflux_ptr*currentCohort%n* AREA_INV * surface_prof(id)
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
