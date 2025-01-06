program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use PRTGenericMod,               only : prt_carbon_allom_hyp
  use FatesInterfaceTypesMod,      only : hlm_parteh_mode
  use PRTParametersMod,            only : prt_params
  use FatesFactoryMod,             only : CohortFactory, PatchFactory
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
  use EDParamsMod,                 only : photo_temp_acclim_timescale
  use EDParamsMod,                 only : photo_temp_acclim_thome_time
  use FatesRadiationMemMod,        only : num_swb
  use FatesConstantsMod,           only : sec_per_day, days_per_year
  use FatesRunningMeanMod,         only : ema_24hr, fixed_24hr, ema_lpa, ema_longterm
  use FatesRunningMeanMod,         only : moving_ema_window, fixed_window
  use EDCohortDynamicsMod,         only : insert_cohort
  
  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader    ! param reader instance
  character(len=:),                  allocatable :: param_file      ! input parameter file
  type(fates_patch_type),            pointer     :: patch           ! patch
  type(fates_cohort_type),           pointer     :: cohort          ! cohort
  type(fates_cohort_type),           pointer     :: shorter_cohort  ! shorter cohort in linked list 
  type(fates_cohort_type),           pointer     :: taller_cohort   ! taller cohort in linked list
  real(r8)                                       :: can_lai(nclmax) ! canopy lai of plot
  integer                                        :: numpft          ! number of pfts (from parameter file)
  integer                                        :: i               ! looping index
  integer                                        :: tnull, snull

  ! CONSTANTS:
  integer  :: num_levsoil = 10      ! number of soil layers
  real(r8) :: step_size = 1800.0_r8 ! step-size [s]
  
  real(r8) :: age = 50.0_r8 
  real(r8) :: area = 500.0_r8 
  integer  :: num_cohorts = 2
  integer  :: pft = 2

  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  numpft = size(prt_params%wood_density, dim=1)
  
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
  
  ! --------------------------------------------------------------------------------------
  
  do i = 1, nclmax
    can_lai(i) = 0.0_r8
  end do

  call PatchFactory(patch, age, area, num_swb, numpft, num_levsoil)
  call CohortFactory(cohort, pft, can_lai, dbh=30.0_r8, number=7.5_r8, age=50.0_r8,      &
    patch_area=area)
    
  ! Put cohort at the right place in the linked list
  taller_cohort => patch%tallest
  shorter_cohort => patch%shortest

  if (associated(patch%tallest)) then
    tnull = 0
  else
    tnull = 1
    patch%tallest => cohort
  endif

  if (associated(patch%shortest)) then
    snull = 0
  else
    snull = 1
    patch%shortest => cohort
  endif
  call insert_cohort(patch, cohort, patch%tallest, patch%shortest, tnull, snull, &
     taller_cohort, shorter_cohort)

  if (associated(cohort)) deallocate(cohort)
  
  call CohortFactory(cohort, pft, can_lai, dbh=25.0_r8, number=7.5_r8, age=50.0_r8,      &
    patch_area=area)
  
    ! Put cohort at the right place in the linked list
  taller_cohort => patch%tallest
  shorter_cohort => patch%shortest

  if (associated(patch%tallest)) then
    tnull = 0
  else
    tnull = 1
    patch%tallest => cohort
  endif

  if (associated(patch%shortest)) then
    snull = 0
  else
    snull = 1
    patch%shortest => cohort
  endif
  call insert_cohort(patch, cohort, patch%tallest, patch%shortest, tnull, snull, &
     taller_cohort, shorter_cohort)

  if (associated(cohort)) deallocate(cohort)
    

end program FatesTestPatch

