module FatesFactoryMod

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesCohortMod,      only : fates_cohort_type
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
  
  implicit none
  
  public :: CohortFactory
  public :: PRTFactory
  
  contains 
  
  subroutine PRTFactory(prt, pft, c_struct, c_leaf, c_fnrt, c_sapw, c_store)
    !
    ! DESCRIPTION:
    ! Create a mock-up of a prt object
    
    ! ARGUMENTS:
    class(prt_vartypes), pointer, intent(inout) :: prt ! PARTEH object
    
    ! LOCALS:
    integer :: el ! looping index
    
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
          m_repro  = 0._r8
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
          write(fates_log(),*) 'Unspecified PARTEH module during create_cohort'
          call endrun(msg=errMsg(__FILE__, __LINE__))
      end select

    end do 

    call prt%CheckInitialConditions()
  
  end subroutine PRTFactory
   
  
  subroutine CohortFactory(cohort, pft)
    !
    ! DESCRIPTION:
    ! Create a mock-up of a cohort
    !

    ! ARGUMENTS
    type(fates_cohort_type), pointer, intent(out)          :: cohort      ! cohort object
    integer,                          intent(in)           :: pft         ! plant funcitonal type index
    integer,                          intent(in), optional :: crowndamage ! crown damage class
    integer,                          intent(in), optional :: status      ! growth status [leaves on/off]
    real(r8),                         intent(in), optional :: dbh         ! diameter at breat height [cm]
    real(r8),                         intent(in), optional :: num         ! number of individuals [/m2]
    real(r8),                         intent(in), optional :: height      ! height [m]
    real(r8),                         intent(in), optional :: age         ! age [yr]
    real(r8),                         intent(in), optional :: ctrim       ! fraction of the maximum leaf biomass
    real(r8),                         intent(in), optional :: spread      ! how spread crowns are in horizontal space
    real(r8),                         intent(in), optional :: elongf_leaf ! leaf elongation factor [fraction]
    real(r8),                         intent(in), optional :: elongf_fnrt ! fine-root "elongation factor" [fraction]
    real(r8),                         intent(in), optional :: elongf_stem ! stem "elongation factor" [fraction]
    
    ! LOCALS:
    class(prt_vartypes), pointer :: prt      ! PARTEH object
    real(r8)                     :: can_area ! canopy area [m2]
    
      
    ! initialize the PRT object
    prt => null()
    call PRTFactory(prt, pft, c_struct, c_leaf, c_fnrt, c_sapw, c_store)
    
    ! allocate the cohort
    allocate(cohort)
    cohort%Create(prt, pft, num, height, age, dbh, status, ctrim, can_area, can_layer,   &
      crowndamage, spread, can_tlai, elongf_leaf, elongf_fnrt, elongf_stem)
  
  
  end subroutine CohortFactory
  
end module FatesFactoryMod