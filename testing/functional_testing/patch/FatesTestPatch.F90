program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use FatesFactoryMod,             only : InitializeGlobals, GetSyntheticPatch
  use SyntheticPatchTypes,         only : synthetic_patch_array_type
  
  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader ! param reader instance
  type(synthetic_patch_array_type)               :: patch_data   ! array of synthetic patches
  character(len=:),                  allocatable :: param_file   ! input parameter file
  type(fates_patch_type),            pointer     :: patch        ! patch
  type(fates_cohort_type),           pointer     :: cohort       ! cohort
  integer                                        :: i            ! patch array location
  
  real(r8) :: rand_heights(4) = (/20.0_r8, 10.0_r8, 1.0_r8, 10.0_r8/)

  ! CONSTANTS:
  integer,  parameter :: num_levsoil = 10      ! number of soil layers
  real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]
  
  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! initialize some global data we need
  call InitializeGlobals(step_size)
  
  ! get all the patch data
  call patch_data%GetSyntheticPatchData()
  
  i = patch_data%PatchDataPosition(patch_name='tropical')
  call GetSyntheticPatch(patch_data%patches(i), num_levsoil, patch)
  
  ! print out list in ascending order
  cohort => patch%shortest
  write(*,*) 'List cohorts in ascending order: '
  do while (associated(cohort))
    write (*,*) cohort%height, cohort%dbh
    cohort => cohort%taller
  end do
  write(*,*) ' '

  ! randomize the heights
  i = 1
  cohort => patch%shortest
  do while (associated(cohort))
    cohort%height = rand_heights(i)
    cohort => cohort%taller
    i = i + 1
  end do 
  cohort => patch%shortest
  write(*,*) 'Now they should be out of order: '
  do while (associated(cohort))
    write (*,*) cohort%height, cohort%dbh
    cohort => cohort%taller
  end do
  write(*,*) ' '
  
  call patch%SortCohorts2()
  
    ! print out list in ascending order
  cohort => patch%shortest
  write(*,*) 'Back in order?'
  do while (associated(cohort))
    write (*,*) cohort%height, cohort%dbh
    cohort => cohort%taller
  end do
  write(*,*) ' '
  
end program FatesTestPatch

