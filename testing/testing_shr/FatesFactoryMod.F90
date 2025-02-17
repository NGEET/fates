module FatesFactoryMod

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : leaves_on, leaves_off
  use FatesConstantsMod,           only : itrue
  use FatesConstantsMod,           only : ihard_stress_decid
  use FatesConstantsMod,           only : isemi_stress_decid
  use FatesConstantsMod,           only : primaryland
  use FatesConstantsMod,           only : sec_per_day, days_per_year
  use FatesGlobals,                only : fates_log
  use FatesGlobals,                only : endrun => fates_endrun
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use EDTypesMod,                  only : init_spread_inventory
  use FatesRadiationMemMod,        only : num_swb
  use EDParamsMod,                 only : ED_val_vai_top_bin_width
  use EDParamsMod,                 only : ED_val_vai_width_increase_factor
  use EDParamsMod,                 only : nlevleaf
  use EDParamsMod,                 only : dinc_vai
  use EDParamsMod,                 only : dlower_vai
  use EDParamsMod,                 only : nclmax
  use EDParamsMod,                 only : photo_temp_acclim_timescale
  use EDParamsMod,                 only : photo_temp_acclim_thome_time
  use FatesRunningMeanMod,         only : ema_24hr, fixed_24hr, ema_lpa, ema_longterm
  use FatesRunningMeanMod,         only : moving_ema_window, fixed_window
  use EDCohortDynamicsMod,         only : InitPRTObject
  use PRTParametersMod,            only : prt_params
  use PRTGenericMod,               only : element_pos
  use PRTGenericMod,               only : num_elements
  use PRTGenericMod,               only : element_list
  use PRTGenericMod,               only : SetState
  use PRTGenericMod,               only : prt_vartypes
  use PRTGenericMod,               only : leaf_organ
  use PRTGenericMod,               only : fnrt_organ
  use PRTGenericMod,               only : sapw_organ
  use PRTGenericMod,               only : store_organ
  use PRTGenericMod,               only : struct_organ
  use PRTGenericMod,               only : repro_organ
  use PRTGenericMod,               only : carbon12_element
  use PRTGenericMod,               only : nitrogen_element
  use PRTGenericMod,               only : phosphorus_element
  use PRTGenericMod,               only : prt_carbon_allom_hyp
  use PRTGenericMod,               only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,               only : StorageNutrientTarget
  use PRTAllometricCarbonMod,      only : InitPRTGlobalAllometricCarbon
  use FatesAllometryMod,           only : h_allom
  use FatesAllometryMod,           only : bagw_allom
  use FatesAllometryMod,           only : bbgw_allom
  use FatesAllometryMod,           only : bleaf
  use FatesAllometryMod,           only : bfineroot
  use FatesAllometryMod,           only : bsap_allom
  use FatesAllometryMod,           only : bdead_allom
  use FatesAllometryMod,           only : bstore_allom
  use FatesAllometryMod,           only : carea_allom
  use FatesInterfaceTypesMod,      only : hlm_parteh_mode
  use FatesInterfaceTypesMod,      only : nleafage
  use FatesSizeAgeTypeIndicesMod,  only : get_age_class_index
  use EDParamsMod,                 only : regeneration_model
  use SyntheticPatchTypes,         only : synthetic_patch_type
  use shr_log_mod,                 only : errMsg => shr_log_errMsg

  implicit none
  
  public :: GetSyntheticPatch
  public :: InitializeGlobals
  
  contains
  
  !---------------------------------------------------------------------------------------
  
  subroutine InitializeGlobals(step_size)
    !
    ! DESCRIPTION:
    ! Initialize globals needed for running factory
    
    ! ARGUMENTS:
    real(r8), intent(in) :: step_size ! step size to use
    
    ! LOCALS:
    integer :: i ! looping index
    
    ! initialize some values
    hlm_parteh_mode = prt_carbon_allom_hyp
    num_elements = 1
    allocate(element_list(num_elements))
    element_list(1) = carbon12_element
    element_pos(:) = 0
    element_pos(carbon12_element) = 1
    call InitPRTGlobalAllometricCarbon()
    
    allocate(ema_24hr)
    call ema_24hr%define(sec_per_day, step_size, moving_ema_window)
    allocate(fixed_24hr)
    call fixed_24hr%define(sec_per_day, step_size, fixed_window)
    allocate(ema_lpa)  
    call ema_lpa%define(photo_temp_acclim_timescale*sec_per_day, step_size,                &
      moving_ema_window)
    allocate(ema_longterm)  
    call ema_longterm%define(photo_temp_acclim_thome_time*days_per_year*sec_per_day,       &
      step_size, moving_ema_window)
      
    do i = 1, nlevleaf
      dinc_vai(i) = ED_val_vai_top_bin_width*ED_val_vai_width_increase_factor**(i-1)
    end do 
  
    do i = 1, nlevleaf
      dlower_vai(i) = sum(dinc_vai(1:i))
    end do
        
  end subroutine InitializeGlobals
  
  !---------------------------------------------------------------------------------------
  
  subroutine PRTFactory(prt, pft, c_struct, c_leaf, c_fnrt, c_sapw, c_store)
    !
    ! DESCRIPTION:
    ! Create a prt object
    
    ! ARGUMENTS:
    class(prt_vartypes), pointer, intent(inout) :: prt      ! PARTEH object
    integer,                      intent(in)    :: pft      ! plant functional type 
    real(r8),                     intent(in)    :: c_struct ! structural carbon [kgC]
    real(r8),                     intent(in)    :: c_leaf   ! leaf carbon [kgC]
    real(r8),                     intent(in)    :: c_fnrt   ! fine root carbon [kgC]
    real(r8),                     intent(in)    :: c_sapw   ! sapwood carbon [kgC]
    real(r8),                     intent(in)    :: c_store  ! storage carbon [kgC]
    
    ! LOCALS:
    integer  :: el         ! looping index
    integer  :: iage       ! looping index
    integer  :: element_id ! element id
    real(r8) :: m_struct   ! mass of structual biomass [kg]
    real(r8) :: m_leaf     ! mass of leaf biomass [kg]
    real(r8) :: m_fnrt     ! mass of fineroot biomass [kg]
    real(r8) :: m_sapw     ! mass of sapwood biomass [kg]
    real(r8) :: m_store    ! mass of storage biomass [kg]
    real(r8) :: m_repro    ! mass of reproductive tissue biomass [kg]
    
    prt => null()

    call InitPRTObject(prt)
        
    do el = 1, num_elements
      
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
          m_repro  = 0.0_r8
        case(nitrogen_element)
          m_struct = c_struct*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(struct_organ))
          m_leaf   = c_leaf*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ))
          m_fnrt   = c_fnrt*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))
          m_sapw   = c_sapw*prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
          m_repro  = 0.0_r8
          m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
        case(phosphorus_element)
          m_struct = c_struct*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(struct_organ))
          m_leaf   = c_leaf*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(leaf_organ))
          m_fnrt   = c_fnrt*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(fnrt_organ))
          m_sapw   = c_sapw*prt_params%phos_stoich_p1(pft, prt_params%organ_param_id(sapw_organ))
          m_repro  = 0.0_r8
          m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
      
      end select

      select case(hlm_parteh_mode)
        
        case (prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp)
          ! Put all of the leaf mass into the first bin
          call SetState(prt, leaf_organ, element_id, m_leaf, 1)
          do iage = 2, nleafage
            call SetState(prt, leaf_organ, element_id, 0.0_r8, iage)
          end do
          call SetState(prt, fnrt_organ, element_id, m_fnrt)
          call SetState(prt, sapw_organ, element_id, m_sapw)
          call SetState(prt, store_organ, element_id, m_store)
          call SetState(prt, struct_organ, element_id, m_struct)
          call SetState(prt, repro_organ, element_id, m_repro)

        case default
          write(fates_log(),*) 'Unspecified PARTEH module'
          call endrun(msg=errMsg(__FILE__, __LINE__))
      end select

    end do 

    call prt%CheckInitialConditions()
  
  end subroutine PRTFactory
   
  !---------------------------------------------------------------------------------------
  
  subroutine CohortFactory(cohort, pft, can_lai, dbh, number, crown_damage, status, age,         &
    canopy_trim, canopy_layer, elong_factor, patch_area)
    !
    ! DESCRIPTION:
    ! Create a FATES cohort
    !

    ! ARGUMENTS
    type(fates_cohort_type), target, intent(out)           :: cohort       ! cohort object
    integer,                          intent(in)           :: pft          ! plant functional type index
    real(r8),                         intent(in)           :: can_lai(:)   ! canopy lai of patch that cohort is on [m2/m2]
    real(r8),                         intent(in), optional :: dbh          ! diameter at breast height [cm]
    real(r8),                         intent(in), optional :: number       ! number density [/m2]
    integer,                          intent(in), optional :: crown_damage ! crown damage class
    integer,                          intent(in), optional :: status       ! growth status [leaves on/off]
    real(r8),                         intent(in), optional :: age          ! age [yr]
    real(r8),                         intent(in), optional :: canopy_trim  ! fraction of the maximum leaf biomass
    integer,                          intent(in), optional :: canopy_layer ! canopy layer
    real(r8),                         intent(in), optional :: elong_factor ! site-level elongation factor
    real(r8),                         intent(in), optional :: patch_area   ! patch area [m2]
    
    ! LOCALS:
    class(prt_vartypes), pointer :: prt                ! PARTEH object
    real(r8)                     :: dbh_local          ! local dbh [cm]
    real(r8)                     :: number_local       ! local number density [/m2]
    integer                      :: crown_damage_local ! local crown damage
    integer                      :: status_local       ! local phenology status
    integer                      :: canopy_layer_local ! local canopy layer
    real(r8)                     :: canopy_trim_local  ! local canopy trim
    real(r8)                     :: age_local          ! local age [yrs]
    real(r8)                     :: elong_fact_local   ! local elongation factor
    real(r8)                     :: patch_area_local   ! local patch area [m2]
    real(r8)                     :: height             ! height [m]
    real(r8)                     :: can_area           ! canopy area [m2]
    real(r8)                     :: c_struct           ! structural carbon [kgC]
    real(r8)                     :: c_leaf             ! leaf carbon [kgC]
    real(r8)                     :: c_fnrt             ! fine root carbon [kgC]
    real(r8)                     :: c_sapw             ! sapwood carbon [kgC]
    real(r8)                     :: c_store            ! storage carbon [kgC]
    real(r8)                     :: c_agw              ! aboveground biomass [kgC]
    real(r8)                     :: c_bgw              ! belowground biomass [kgC]
    real(r8)                     :: a_sapw             ! sapwood area [m2]
    real(r8)                     :: elongf_leaf        ! leaf elongation factor [fraction]
    real(r8)                     :: elongf_fnrt        ! fine-root "elongation factor" [fraction]
    real(r8)                     :: elongf_stem        ! stem "elongation factor" [fraction]
    
    ! CONSTANTS:
    real(r8), parameter :: dbh_default = 10.0_r8         ! default dbh [cm]
    real(r8), parameter :: canopy_trim_default = 1.0_r8  ! default canopy trim 
    integer,  parameter :: crown_damage_default = 1      ! default crown damage
    integer,  parameter :: status_default = leaves_on    ! default status
    real(r8), parameter :: age_default = 25.0_r8         ! default age [yrs]
    real(r8), parameter :: elong_factor_default = 1.0_r8 ! default elongation factor
    integer,  parameter :: canopy_layer_default = 1      ! default canopy layer
    real(r8), parameter :: patch_area_default = 100.0_r8 ! default patch area [m2]
    
    ! set local values
    if (present(dbh)) then
      dbh_local = dbh 
    else
      dbh_local = dbh_default
    end if
    
    if (present(crown_damage)) then 
      crown_damage_local = crown_damage
    else 
      crown_damage_local = crown_damage_default
    end if
    
    if (present(canopy_trim)) then 
      canopy_trim_local = canopy_trim
    else 
      canopy_trim_local = canopy_trim_default
    end if
    
    if (present(status)) then
      status_local = status 
    else
      status_local = status_default
    end if
    
    if (present(age)) then 
      age_local = age
    else 
      age_local = age_default 
    end if
    
    if (present(elong_factor)) then 
      elong_fact_local = elong_factor 
    else 
      elong_fact_local = elong_factor_default
    end if
    
    if (present(patch_area)) then 
      patch_area_local = patch_area 
    else 
      patch_area_local = patch_area_default 
    end if
    
    if (present(canopy_layer)) then 
      canopy_layer_local = canopy_layer 
    else 
      canopy_layer_local = canopy_layer_default 
    end if
    
    ! set leaf elongation factors
    if (prt_params%season_decid(pft) == itrue .and. status_local == leaves_off) then
      elongf_leaf = 0.0_r8
      elongf_fnrt = 1.0_r8 - prt_params%phen_fnrt_drop_fraction(pft) 
      elongf_stem = 1.0_r8 - prt_params%phen_stem_drop_fraction(pft)
    
    else if (any(prt_params%stress_decid(pft) == [ihard_stress_decid, isemi_stress_decid])) then
      elongf_leaf = elong_fact_local
      elongf_fnrt = 1.0_r8 - (1.0_r8 - elongf_leaf)*prt_params%phen_fnrt_drop_fraction(pft)
      elongf_stem = 1.0_r8 - (1.0_r8 - elongf_leaf)*prt_params%phen_stem_drop_fraction(pft)
      
      if (elongf_leaf > 0.0_r8) then 
        status_local = leaves_on
      else 
        status_local = leaves_off
      end if
    
    else
      elongf_leaf = 1.0_r8
      elongf_fnrt = 1.0_r8 
      elongf_stem = 1.0_r8 
    end if
    
    ! calculate allometric properties
    
    ! calculate or set number density
    if (present(number)) then
      
      number_local = number
    
    else 
      
      call carea_allom(dbh_local, 1.0_r8, init_spread_inventory, pft, crown_damage_local,  &
        can_area)
    
      ! calculate initial density required to close the canopy
      number_local = patch_area_local/can_area
      
    end if
    
    ! calculate leaf biomass
    call bleaf(dbh_local, pft, crown_damage_local, canopy_trim_local, elongf_leaf, c_leaf)
    
    ! recalculate crown area
    call carea_allom(dbh_local, number_local, init_spread_inventory, pft,                &
      crown_damage_local, can_area)
    
    ! calculate height  
    call h_allom(dbh_local, pft, height)
    
    ! calculate total above-ground biomass 
    call bagw_allom(dbh_local, pft, crown_damage_local, elongf_stem, c_agw)
    
    ! calculate coarse root biomass 
    call bbgw_allom(dbh_local, pft, elongf_stem, c_bgw)
    
    ! fine root biomass from allometry
    call bfineroot(dbh_local, pft, canopy_trim_local, prt_params%allom_l2fr(pft),        &
      elongf_fnrt, c_fnrt)
    
    ! sapwood biomass
    call bsap_allom(dbh_local, pft, crown_damage_local, canopy_trim_local, elongf_stem,  &
      a_sapw, c_sapw)
    
    ! structural biomass
    call bdead_allom(c_agw, c_bgw, c_sapw, pft, c_struct)
    
    ! storage biomass
    call bstore_allom(dbh_local, pft, crown_damage_local, canopy_trim_local, c_store)
        
    ! initialize the PRT object
    call PRTFactory(prt, pft, c_struct, c_leaf, c_fnrt, c_sapw, c_store)
    
    ! create the cohort
    call cohort%Create(prt, pft, number_local, height, age_local, dbh_local,             &
      status_local, canopy_trim_local, can_area, canopy_layer_local, crown_damage_local, &
      init_spread_inventory, can_lai, elongf_leaf, elongf_fnrt, elongf_stem)
  
  end subroutine CohortFactory
  
  !---------------------------------------------------------------------------------------
  
  subroutine PatchFactory(patch, age, area, num_swb, num_pft, num_levsoil,               &
    land_use_label, nocomp_pft, current_tod)
    !
    ! DESCRIPTION:
    ! Create a fates patch
    !
    
    ! ARGUMENTS:
    type(fates_patch_type), pointer, intent(out)           :: patch          ! patch object
    real(r8),                        intent(in)            :: age            ! patch age [yrs]
    real(r8),                        intent(in)            :: area           ! patch are [m2]
    integer,                         intent(in)            :: num_swb        ! number of shortwave bands
    integer,                         intent(in)            :: num_pft        ! number of pfts
    integer,                         intent(in)            :: num_levsoil    ! number of soil layers
    integer,                         intent(in), optional  :: land_use_label ! land use label
    integer,                         intent(in), optional  :: nocomp_pft     ! nocomp_pft label
    integer,                         intent(in), optional  :: current_tod    ! time of day [seconds past 0Z]
    
    ! LOCALS:
    integer :: land_use_label_local ! local land use label
    integer :: nocomp_pft_local     ! local nocomp pft label
    integer :: tod_local            ! local tod value
    
    ! CONSTANTS:
    integer :: land_use_label_default = primaryland ! default land use label
    integer :: nocomp_pft_default = 1               ! default nocomp pft label
    integer :: tod_default = 0                      ! default time of day
    
    ! set defaults if necessary
    if (present(land_use_label)) then 
      land_use_label_local = land_use_label
    else 
      land_use_label_local = land_use_label_default
    end if
    
    if (present(nocomp_pft)) then 
      nocomp_pft_local = nocomp_pft
    else 
      nocomp_pft_local = nocomp_pft_default
    end if
    
    if (present(current_tod)) then
      tod_local = current_tod
    else 
      tod_local = tod_default
    end if
    
    allocate(patch)
    call patch%Create(age, area, land_use_label_local, nocomp_pft_local, num_swb,        &
      num_pft, num_levsoil, tod_local, regeneration_model)
    
    patch%patchno = 1
    patch%younger => null()
    patch%older => null()
    patch%age_class = get_age_class_index(patch%age)

  end subroutine PatchFactory
  
  !---------------------------------------------------------------------------------------
  
  subroutine GetSyntheticPatch(patch_data, num_levsoil, patch)
    !
    ! DESCRIPTION:
    ! Create a synthetic patch based on input data
    !
    
    ! ARGUMETNS:
    type(synthetic_patch_type),          intent(in)  :: patch_data  ! synthetic patch data
    integer,                             intent(in)  :: num_levsoil ! number of soil layers
    type(fates_patch_type),     pointer, intent(out) :: patch       ! patch
    
    ! LOCALS:
    type(fates_cohort_type), pointer :: cohort          ! cohort object
    integer                          :: numpft          ! total number of pfts
    real(r8)                         :: can_lai(nclmax) ! canopy lai of plot
    real(r8)                         :: patch_age       ! patch age
    integer                          :: i               ! looping index
  
    numpft = size(prt_params%wood_density, dim=1)
    patch_age = maxval(patch_data%ages(:))
    can_lai(:) = 0.0_r8
    
    ! create the patch 
    call PatchFactory(patch, patch_age, patch_data%area, num_swb, numpft, num_levsoil)
    
    ! add cohorts
    do i = 1, patch_data%num_cohorts 
      allocate(cohort)
      call CohortFactory(cohort, patch_data%pft_ids(i), can_lai, dbh=patch_data%dbhs(i), &
        number=patch_data%densities(i)*patch_data%area, age=patch_data%ages(i),          &
        canopy_layer=patch_data%canopy_layers(i), patch_area=patch_data%area)
      
      call patch%InsertCohort(cohort)
      
    end do   
  
  end subroutine GetSyntheticPatch
  
  !---------------------------------------------------------------------------------------
  
  subroutine CreateTestPatchList(patch, heights, dbhs)
    !
    ! DESCRIPTION:
    ! Create a patch with a hard-coded cohort linked list
    ! Heights are supplied, optional dbhs
    ! Used for unit testing
    !
    
    ! ARGUMENTS:
    type(fates_patch_type), intent(out)          :: patch      ! patch object
    real(r8),               intent(in)           :: heights(:) ! hard-coded heights
    real(r8),               intent(in), optional :: dbhs(:)    ! optional hard-coded dbhs
    
    ! LOCALS:
    type(fates_cohort_type), pointer :: cohort, next_cohort ! cohort objects
    integer                          :: num_cohorts         ! number of cohorts to add to list 
    integer                          :: i                   ! looping index
    
    ! size of heights array must match sie of dbhs, if supplied
    if (present(dbhs)) then 
      if (size(heights) /= size(dbhs)) then 
        write(*, '(a)') "Size of heights array must match size of dbh array."
        stop
      end if 
    end if 
    
    num_cohorts = size(heights)
    
    ! initialize first cohort
    allocate(cohort)
    cohort%height = heights(1)
    if (present(dbhs)) cohort%dbh = dbhs(1)
    patch%shortest => cohort
    
    ! initialize the rest of the cohorts
    do i = 2, num_cohorts
      allocate(next_cohort)
      next_cohort%height = heights(i)
      if (present(dbhs)) next_cohort%dbh = dbhs(i)
      cohort%taller => next_cohort
      next_cohort%shorter => cohort
      cohort => next_cohort
    end do
    patch%tallest => cohort
  
  end subroutine CreateTestPatchList
  
end module FatesFactoryMod