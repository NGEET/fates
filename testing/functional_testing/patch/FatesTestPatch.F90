program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : itrue
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use FatesFactoryMod,             only : InitializeGlobals, GetSyntheticPatch
  use SyntheticPatchTypes,         only : synthetic_patch_array_type
  use EDCanopyStructureMod,        only : UpdatePatchLAI
  
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
  real(r8) :: site_spread
  real(r8) :: sitelevel_canopyarea, inc
  
  ! ! insert 100 cohorts of the same height
  ! allocate(patch)
  
  ! do i = 1, num_cohorts
    
  !   allocate(cohort)
  !   cohort%height = height
  !   cohort%dbh = i
    
  !   call patch%InsertCohort(cohort)

  ! end do
  
  ! cohort => patch%shortest
  ! do while (associated(cohort))
  !   print *, cohort%dbh
  !   cohort => cohort%taller
  ! end do
  
  ! print *, '------------------'
  
  ! call patch%SortCohorts()
  
  ! print *, '------------------'
  
  ! cohort => patch%shortest
  ! do while (associated(cohort))
  !   print *, cohort%dbh
  !   cohort => cohort%taller
  ! end do

  !read in parameter file name from command line
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
  write(*,*) patch%canopy_layer_tlai(:)
  
    ! print out list in ascending order
  cohort => patch%shortest
  write(*,*) 'Before UpdateLAI'
  do while (associated(cohort))
    write (*,*) cohort%dbh, cohort%treelai
    cohort => cohort%taller
  end do
  write(*,*) ' '
  
  call UpdatePatchLAI(patch)
  
  ! print out list in ascending order
  cohort => patch%shortest
  write(*,*) 'After UpdateLAI'
  do while (associated(cohort))
    write (*,*) cohort%dbh, cohort%treelai, cohort%nv
    cohort => cohort%taller
  end do
  write(*,*) ' '
  write(*,*) patch%canopy_layer_tlai(:)
  


end program FatesTestPatch

