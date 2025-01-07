program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use EDParamsMod,                 only : nclmax
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use PRTParametersMod,            only : prt_params  
  use FatesFactoryMod,             only : CohortFactory, PatchFactory, InitializeGlobals
  use FatesRadiationMemMod,        only : num_swb
  use EDCohortDynamicsMod,         only : insert_cohort, insert_cohort_2
  use TimingMod,                   only : testing_timer
  use RandMod,                     only : random0
  
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
  type(testing_timer)                            :: timer           ! timer object
  real(r8)                                       :: elapsed_time    ! elapsed time [ms]
  integer                                        :: n = 1000
  real(r8)                                       :: rand_dbh

  ! CONSTANTS:
  integer  :: num_levsoil = 10      ! number of soil layers
  real(r8) :: step_size = 1800.0_r8 ! step-size [s]
  
  real(r8) :: max_age 
  
  ! all subject to patch/biome type
  real(r8) :: age(6) = (/50.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, 50.0_r8/)
  real(r8) :: dbh(6) = (/30.0_r8, 25.0_r8, 10.0_r8, 4.0_r8, 50.0_r8, 30.0_r8/)
  real(r8) :: dens(6) = (/7.5_r8, 7.5_r8, 7.5_r8, 7.5_r8, 7.5_r8, 7.5_r8/)
  integer  :: can_layer(6) = (/1, 1, 1, 1, 1, 1/)
  integer  :: pft(6) = (/2, 2, 2, 2, 2, 2/)
  real(r8) :: area = 500.0_r8
  integer  :: num_cohorts = 6
  
  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! initialize some global data we need
  call InitializeGlobals(step_size)
  
  numpft = size(prt_params%wood_density, dim=1)

  call PatchFactory(patch, 50.0_r8, area, num_swb, numpft, num_levsoil)
  
  ! just set this to 0.0  
  can_lai(:) = 0.0_r8
  
  do i = 1, 100
    
    allocate(cohort)
    call CohortFactory(cohort, 2, can_lai, dbh=dbh(i), number=7.5_r8, age=50.0_r8,       &
      canopy_layer=1, patch_area=area)

    call insert_cohort_2(patch, cohort)
    
  end do
    
end program FatesTestPatch

