module FatesFactoryMod

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : leaves_on, leaves_off
  use FatesCohortMod,      only : fates_cohort_type
  use EDTypesMod,          only : init_spread_inventory
  use EDCohortDynamicsMod, only : InitPRTObject
  use PRTParametersMod,    only : prt_params
  use PRTGenericMod,       only : num_elements
  use PRTGenericMod,       only : element_list
  use PRTGenericMod,       only : SetState
  use PRTGenericMod,       only : prt_vartypes
  use PRTGenericMod,       only : leaf_organ
  use PRTGenericMod,       only : fnrt_organ
  use PRTGenericMod,       only : sapw_organ
  use PRTGenericMod,       only : store_organ
  use PRTGenericMod,       only : struct_organ
  use PRTGenericMod,       only : repro_organ
  use PRTGenericMod,       only : carbon12_element
  use PRTGenericMod,       only : nitrogen_element
  use PRTGenericMod,       only : phosphorus_element
  use PRTGenericMod,       only : StorageNutrientTarget
  use FatesAllometryMod,   only : h_allom
  use FatesAllometryMod,   only : bagw_allom
  use FatesAllometryMod,   only : bbgw_allom
  use FatesAllometryMod,   only : bleaf
  use FatesAllometryMod,   only : bfineroot
  use FatesAllometryMod,   only : bsap_allom
  use FatesAllometryMod,   only : bdead_allom
  use FatesAllometryMod,   only : bstore_allom
  use FatesAllometryMod,   only : carea_allom
  
  implicit none
  
  public :: CohortFactory
  public :: PRTFactory
  
  ! CONSTANTS
  real(r8) :: patch_area_default = 100.0_r8 ! default patch area [m2]
  
  contains 
  
  subroutine PRTFactory(prt, pft, c_struct, c_leaf, c_fnrt, c_sapw, c_store)
    !
    ! DESCRIPTION:
    ! Create a mock-up of a prt object
    
    ! ARGUMENTS:
    class(prt_vartypes), pointer, intent(inout) :: prt ! PARTEH object
    
    ! LOCALS:
    integer :: el ! looping index
    
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
        case (prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp )
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
  
  subroutine CohortFactory(cohort, pft, dbh, crown_damage, status, age, canopy_trim,     &
    spread)
    !
    ! DESCRIPTION:
    ! Create a mock-up of a cohort
    !

    ! ARGUMENTS
    type(fates_cohort_type), pointer, intent(out)          :: cohort       ! cohort object
    integer,                          intent(in)           :: pft          ! plant functional type index
    real(r8),                         intent(in), optional :: dbh          ! diameter at breast height [cm]
    integer,                          intent(in), optional :: crown_damage ! crown damage class
    integer,                          intent(in), optional :: status       ! growth status [leaves on/off]
    real(r8),                         intent(in), optional :: age          ! age [yr]
    real(r8),                         intent(in), optional :: canopy_trim  ! fraction of the maximum leaf biomass
    real(r8),                         intent(in), optional :: spread       ! how spread crowns are in horizontal space
    real(r8),                         intent(in), optional :: elong_factor ! site-level elongation factor
    real(r8),                         intent(in), optional :: patch_area   ! patch area [m2]
    
    ! LOCALS:
    class(prt_vartypes), pointer :: prt                ! PARTEH object
    real(r8)                     :: dbh_local          ! local dbh [cm]
    integer                      :: crown_damage_local ! local crown damage
    integer                      :: status_local       ! local phenology status
    real(r8)                     :: canopy_trim_local  ! local canopy trim
    real(r8)                     :: age_local          ! local age [yrs]
    real(r8)                     :: elong_fact_local   ! local elongation factor
    real(r8)                     :: patch_area_local   ! local patch area [m2]
    real(r8)                     :: height             ! height [m]
    real(r8)                     :: number             ! number density [/m2]
    real(r8)                     :: can_area           ! canopy area [m2]
    real(r8)                     :: c_struct           ! structural carbon [kgC]
    real(r8)                     :: c_leaf             ! leaf carbon [kgC]
    real(r8)                     :: c_fnrt             ! fine root carbon [kgC]
    real(r8)                     :: c_sapw             ! sapwood carbon [kgC]
    real(r8)                     :: c_store            ! storage carbon [kgC]
    real(r8)                     :: c_agw              ! aboveground biomass [kgC]
    real(r8)                     :: c_bgw              ! belowground biomass [kgC]
    real(r8)                     :: a_sapw             ! sapwood area [m2]
    real(r8),                    :: elongf_leaf        ! leaf elongation factor [fraction]
    real(r8),                    :: elongf_fnrt        ! fine-root "elongation factor" [fraction]
    real(r8),                    :: elongf_stem        ! stem "elongation factor" [fraction]
    
    ! CONSTANTS:
    real(r8) :: dbh_default = 10.0_r8         ! default dbh [cm]
    real(r8) :: canopy_trim_default = 1.0_r8  ! default canopy trim 
    integer  :: crown_damage_default = 1      ! default crown damage
    integer  :: status_default = leaves_on    ! default status
    real(r8) :: age_default = 25.0_r8         ! default age [yrs]
    real(r8) :: elong_factor_default = 1.0_r8 ! default elongation factor
    
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
    
    ! set leaf elongation factors
    if (prt_params%season_decid(pft) == itrue .and. status_local == leaves_off) then
      elongf_leaf = 0.0_r8
      elongf_fnrt = 1.0_r8 - prt_params%phen_fnrt_drop_fraction(pft) 
      elongf_stem = 1.0_r8 - prt_params%phen_stem_drop_fraction(pft)
    
    else if (any(prt_params%stress_decid(pft) == [ihard_stress_decid, isemi_stress_decid]))
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
    
    ! calculate crown area of a single plant
    call carea_allom(dbh_local, 1.0_r8, init_spread_inventory, pft, crown_damage_local,  &
      can_area)
    
    ! calculate initial density required to close the canopy
    number = patch_area_local/can_area
    
    ! calculate leaf biomass
    call bleaf(dbh_local, pft, crown_damage_local, canopy_trim_local, efleaf_coh, c_leaf)
    
    ! recalculate crown area
    call carea_allom(dbh_local, number, init_spread_inventory, pft, crown_damage_local,  &
      can_area)
    
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
    
    ! allocate the cohort
    allocate(cohort)
    cohort%Create(prt, pft, number, height, age_local, dbh_local, status_local,          &
      canopy_trim_local, can_area, can_layer, crown_damage_local, spread, can_tlai,      &
      elongf_leaf, elongf_fnrt, elongf_stem)
  
  end subroutine CohortFactory
  
  !---------------------------------------------------------------------------------------
  
end module FatesFactoryMod