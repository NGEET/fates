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
  type(fates_cohort_type),               pointer :: head, new_node
  real(r8) :: heights(8) = (/10.0_r8, 100.0_r8, 15.0_r8, 2.0_r8, 1.0_r8, 12.5001_r8, 20.0_r8, 0.5_r8/)
  
     ! create a hardcoded doubly linked list
      allocate(patch)
      
      allocate(head)
    
      head%height = heights(1)
      patch%shortest => head
    
      allocate(new_node)
      new_node%height = heights(2)
      head%taller => new_node
      
      allocate(new_node)
      new_node%height = heights(3)
      head%taller%taller => new_node
      new_node%shorter => head%taller
      
      allocate(new_node)
      new_node%height = heights(4)
      head%taller%taller%taller => new_node
      new_node%shorter => head%taller%taller
      
      allocate(new_node)
      new_node%height = heights(5)
      head%taller%taller%taller%taller => new_node
      new_node%shorter => head%taller%taller%taller
      
      allocate(new_node)
      new_node%height = heights(6)
      head%taller%taller%taller%taller%taller => new_node
      new_node%shorter => head%taller%taller%taller%taller
      
      allocate(new_node)
      new_node%height = heights(7)
      head%taller%taller%taller%taller%taller%taller => new_node
      new_node%shorter => head%taller%taller%taller%taller%taller
      
      allocate(new_node)
      new_node%height = heights(8)
      head%taller%taller%taller%taller%taller%taller%taller => new_node
      new_node%shorter => head%taller%taller%taller%taller%taller%taller
    
      patch%tallest => new_node
      
      ! sort cohorts
      call patch%SortCohorts()
      
      ! check backwards and forwards
      new_node => patch%shortest
      do while (associated(new_node))
        if (associated(new_node%taller)) then
          print *, new_node%height
        end if
        new_node => new_node%taller
      end do
      
      new_node => patch%tallest
      do while (associated(new_node))
        if (associated(new_node%shorter)) then 
        print *, new_node%height
        end if
        new_node => new_node%shorter
      end do

  ! ! CONSTANTS:
  ! integer,  parameter :: num_levsoil = 10      ! number of soil layers
  ! real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]
  
  ! ! read in parameter file name from command line
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

  ! ! randomize the heights
  ! i = 1
  ! cohort => patch%shortest
  ! do while (associated(cohort))
  !   cohort%height = rand_heights(i)
  !   cohort => cohort%taller
  !   i = i + 1
  ! end do 
  ! cohort => patch%shortest
  ! write(*,*) 'Now they should be out of order: '
  ! do while (associated(cohort))
  !   write (*,*) cohort%height, cohort%dbh
  !   cohort => cohort%taller
  ! end do
  ! write(*,*) ' '
  
  ! call patch%SortCohorts2()
  
  !   ! print out list in ascending order
  ! cohort => patch%shortest
  ! write(*,*) 'Back in order?'
  ! do while (associated(cohort))
  !   write (*,*) cohort%height, cohort%dbh
  !   cohort => cohort%taller
  ! end do
  ! write(*,*) ' '
  
end program FatesTestPatch

