program FatesTestCohort

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use PRTGenericMod,               only : prt_carbon_allom_hyp
  use FatesInterfaceTypesMod,      only : hlm_parteh_mode
  use PRTParametersMod,            only : prt_params
  use FatesFactoryMod,             only : CohortFactory
  use PRTGenericMod,               only : element_list
  use PRTGenericMod,               only : element_pos
  use PRTGenericMod,               only : num_elements
  use PRTGenericMod,               only : carbon12_element
  use PRTAllometricCarbonMod,      only : InitPRTGlobalAllometricCarbon
  use EDParamsMod,                 only : ED_val_vai_top_bin_width
  use EDParamsMod,                 only : ED_val_vai_width_increase_factor
  use EDParamsMod,                 only : nlevleaf, nclmax
  use EDParamsMod,                 only : dinc_vai
  use EDParamsMod,                 only : dlower_vai
  
  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader    ! param reader instance
  character(len=:),                  allocatable :: param_file      ! input parameter file
  type(fates_cohort_type),           pointer     :: cohort          ! cohort
  real(r8)                                       :: can_lai(nclmax) ! canopy lai of plot
  integer                                        :: i               ! looping index

  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! initialize some values
  hlm_parteh_mode = prt_carbon_allom_hyp
  num_elements = 1
  allocate(element_list(num_elements))
  element_list(1) = carbon12_element
  element_pos(:) = 0
  element_pos(carbon12_element) = 1
  call InitPRTGlobalAllometricCarbon()

  do i = 1, nlevleaf
    dinc_vai(i) = ED_val_vai_top_bin_width*ED_val_vai_width_increase_factor**(i-1)
  end do 
  
  do i = 1, nlevleaf
    dlower_vai(i) = sum(dinc_vai(1:i))
  end do
  
  do i = 1, nclmax
    can_lai(i) = 1.0_r8
  end do
  
  call CohortFactory(cohort, 1, can_lai, canopy_layer=2)
  
end program FatesTestCohort

