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
  type(fates_cohort_type), pointer               :: storebigcohort
  type(fates_cohort_type), pointer               :: storesmallcohort
  integer                                        :: tnull, snull

  ! CONSTANTS:
  integer,  parameter :: num_levsoil = 10      ! number of soil layers
  real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]
  integer,  parameter :: num_cohorts = 2
  real(r8), parameter :: height = 35.0
  
  ! insert 100 cohorts of the same height
  allocate(patch)
  
  do i = 1, num_cohorts
    
    allocate(cohort)
    cohort%height = height
    cohort%dbh = i
    
    storebigcohort => patch%tallest
    storesmallcohort => patch%shortest
  
    if(associated(patch%tallest)) then
        tnull = 0
      else
        tnull = 1
        patch%tallest => cohort
        cohort%taller => null()
      endif

      if(associated(patch%shortest))then
        snull = 0
      else
        snull = 1
        patch%shortest => cohort
        cohort%shorter => null()
      endif 
      
      call patch%InsertCohort(cohort)
      ! call patch%InsertCohort_old(cohort, patch%tallest,           &
      !  patch%shortest, tnull, snull, storebigcohort, storesmallcohort)
      
      ! patch%tallest  => storebigcohort
      ! patch%shortest => storesmallcohort
  end do
  
  cohort => patch%shortest
  do while (associated(cohort))
    print *, cohort%dbh
    cohort => cohort%taller
  end do
  
  print *, '------------------'
  
  call patch%SortCohorts()
  
  print *, '------------------'
  
  cohort => patch%shortest
  do while (associated(cohort))
    print *, cohort%dbh
    cohort => cohort%taller
  end do

  ! read in parameter file name from command line
  ! param_file = command_line_arg(1)
  
  ! ! read in parameter file
  ! call param_reader%Init(param_file)
  ! call param_reader%RetrieveParameters()
  
  ! ! initialize some global data we need
  ! call InitializeGlobals(step_size)
  
  ! ! get all the patch data
  ! call patch_data%GetSyntheticPatchData()
  
  ! i = patch_data%PatchDataPosition(patch_name='tropical')
  ! call GetSyntheticPatch(patch_data%patches(i), num_levsoil, patch)
  
  ! ! print out list in ascending order
  ! cohort => patch%shortest
  ! write(*,*) 'List cohorts in ascending order: '
  ! do while (associated(cohort))
  !   write (*,*) cohort%height, cohort%dbh
  !   cohort => cohort%taller
  ! end do
  ! write(*,*) ' '

end program FatesTestPatch

