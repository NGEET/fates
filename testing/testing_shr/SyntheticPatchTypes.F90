module SyntheticPatchTypes

  ! DESCRIPTION:
  !  Methods to create synthetic patch objects

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private
  
  integer, parameter :: chunk_size = 10
  
  ! patch data type to hold data about the synthetic patches
  type, public :: synthetic_patch_type
    character(len=:), allocatable :: patch_name        ! patch name, used for reference
    integer                       :: patch_id          ! patch id, used for reference
    integer                       :: num_cohorts       ! number of cohorts on patch
    real(r8)                      :: area              ! patch area [m2]
    real(r8), allocatable         :: ages(:)           ! cohort ages [yr]
    real(r8), allocatable         :: dbhs(:)           ! cohort dbhs [cm]
    real(r8), allocatable         :: densities(:)      ! cohort densities [/m2]
    integer,  allocatable         :: canopy_layers(:)  ! canopy layers
    integer,  allocatable         :: pft_ids(:)        ! pft ids
    
    contains 
      
      procedure :: InitSyntheticPatchData
  
  end type synthetic_patch_type
  
  ! --------------------------------------------------------------------------------------
  
  ! a class to just hold an array of these synthetic patches
  type, public :: synthetic_patch_array_type
  
    type(synthetic_patch_type), allocatable :: patches(:)  ! array of patches
    integer                                 :: num_patches ! total number of patches
    
    contains 
      
      procedure :: AddPatch
      procedure :: GetSyntheticPatchData
      procedure :: PatchDataPosition
  
  end type synthetic_patch_array_type
  
  ! --------------------------------------------------------------------------------------
  
  contains 
  
  subroutine InitSyntheticPatchData(this, patch_id, patch_name, area, ages, dbhs,        &
    densities, pft_ids, canopy_layers)
    !
    ! DESCRIPTION:
    ! Initializes a synthetic patch with input characteristics
    !
    
    ! ARGUMENTS:
    class(synthetic_patch_type), intent(inout) :: this             ! patch data to create
    integer,                     intent(in)    :: patch_id         ! patch id
    character(len=*),            intent(in)    :: patch_name       ! patch name 
    real(r8),                    intent(in)    :: area             ! patch area [m2]
    real(r8),                    intent(in)    :: ages(:)          ! cohort ages [yr]
    real(r8),                    intent(in)    :: dbhs(:)          ! cohort dbhs [cm]
    real(r8),                    intent(in)    :: densities(:)     ! cohort densities [/m2]
    integer,                     intent(in)    :: pft_ids(:)       ! pft ids
    integer,                     intent(in)    :: canopy_layers(:) ! canopy layers of cohorts
    
    ! LOCALS:
    integer :: num_cohorts ! number of cohorts on patch 
    integer :: i           ! looping index
    
    ! allocate arrays 
    num_cohorts = size(pft_ids)
    
    allocate(this%ages(num_cohorts))
    allocate(this%dbhs(num_cohorts))
    allocate(this%densities(num_cohorts))
    allocate(this%pft_ids(num_cohorts))
    allocate(this%canopy_layers(num_cohorts))
    
    ! set values
    this%patch_name = patch_name
    this%num_cohorts = num_cohorts
    this%patch_id = patch_id
    this%area = area
    
    do i = 1, num_cohorts 
      this%ages(i) = ages(i)
      this%dbhs(i) = dbhs(i)
      this%densities(i) = densities(i)
      this%pft_ids(i) = pft_ids(i)
      this%canopy_layers(i) = canopy_layers(i)
    end do 
        
  end subroutine InitSyntheticPatchData
  
  ! --------------------------------------------------------------------------------------
  
   subroutine AddPatch(this, patch_id, patch_name, area, ages, dbhs, densities, pft_ids, &
    canopy_layers)
    !
    ! DESCRIPTION:
    ! Adds a synthetic patch data to a dynamic array
    !
    
    ! ARGUMENTS:
    class(synthetic_patch_array_type), intent(inout) :: this              ! array of synthetic patches
    integer,                           intent(in)    :: patch_id          ! patch id
    character(len=*),                  intent(in)    :: patch_name        ! name of patch
    real(r8),                          intent(in)    :: area              ! patch area
    real(r8),                          intent(in)    :: ages(:)           ! cohort ages [yr]
    real(r8),                          intent(in)    :: dbhs(:)           ! cohort dbhs [cm]
    real(r8),                          intent(in)    :: densities(:)      ! cohort densities [/m2]
    integer,                           intent(in)    :: pft_ids(:)        ! pft ids
    integer,                           intent(in)    :: canopy_layers(:)  ! canopy layers
    
    ! LOCALS:
    type(synthetic_patch_type)              :: patch_data         ! synthetic patch data
    type(synthetic_patch_type), allocatable :: temporary_array(:) ! temporary array to hold data while re-allocating
    
    ! first make sure we have enough space in the array
    if (allocated(this%patches)) then
      ! already allocated to some size
      if (this%num_patches == size(this%patches)) then 
        ! need to add more space
        allocate(temporary_array(size(this%patches) + chunk_size))
        temporary_array(1:size(this%patches)) = this%patches
        call move_alloc(temporary_array, this%patches)
      end if 
      
      this%num_patches = this%num_patches + 1
  
    else 
      ! first element in array 
      allocate(this%patches(chunk_size))
      this%num_patches = 1
    end if 
    
    call patch_data%InitSyntheticPatchData(patch_id, patch_name, area, ages, dbhs,       &
      densities, pft_ids, canopy_layers)
    
    this%patches(this%num_patches) = patch_data
      
  end subroutine AddPatch
  
  ! --------------------------------------------------------------------------------------
  
  integer function PatchDataPosition(this, patch_id, patch_name)
    !
    ! DESCRIPTION:
    ! Returns the index of a desired synthetic patch data 
    !
    
    ! ARGUMENTS:
    class(synthetic_patch_array_type), intent(in)            :: this       ! array of patch data
    integer,                           intent(in), optional  :: patch_id   ! desired patch id
    character(len=*),                  intent(in), optional  :: patch_name ! desired patch name
        
    ! LOCALS:
    integer :: i ! looping index
    
    ! can't supply both
    if (present(patch_id) .and. present(patch_name)) then 
      write(*, '(a)') "Can only supply either a patch_id or a patch_name - not both"
      stop
    end if 
    
    do i = 1, this%num_patches
      if (present(patch_id)) then 
        if (this%patches(i)%patch_id == patch_id) then
          PatchDataPosition = i
          return
        end if
      else if (present(patch_name)) then 
        if (this%patches(i)%patch_name == patch_name) then
          PatchDataPosition = i
          return
        end if
      else 
        write(*, '(a)') "Must supply either a patch_id or a patch_name."
        stop
      end if
    end do
    write(*, '(a)') "Cannot find the synthetic patch type supplied"
    stop
  
  end function PatchDataPosition
  
  ! --------------------------------------------------------------------------------------
  
  subroutine GetSyntheticPatchData(this)
    !
    ! DESCRIPTION:
    ! Returns an array of hard-coded synthetic patch data
    ! 
    !
    
    ! ARGUMENTS:
    class(synthetic_patch_array_type), intent(inout) :: this ! array of synthetic patches
    
    call this%AddPatch(patch_id=1, patch_name='tropical', area=500.0_r8,                 &
      ages=(/100.0_r8, 80.0_r8, 40.0_r8, 20.0_r8/),                                      &
      dbhs=(/60.0_r8, 50.0_r8, 25.0_r8, 10.0_r8/),                                       &
      densities=(/0.005_r8, 0.008_r8, 0.02_r8, 0.017_r8/),                               &
      pft_ids=(/1, 1, 1, 1/),                                                            &
      canopy_layers=(/1, 1, 2, 2/))
    
    call this%AddPatch(patch_id=2, patch_name='evergreen', area=500.0_r8,                &
      ages=(/50.0_r8, 50.0_r8/),                                                         &
      dbhs=(/30.0_r8, 25.0_r8/),                                                         &
      densities=(/0.015_r8, 0.015_r8/),                                                  &
      pft_ids=(/2, 2/),                                                                  &
      canopy_layers=(/1, 1/))
      
    call this%AddPatch(patch_id=3, patch_name='savannah', area=500.0_r8,                 &
      ages=(/20.0_r8, 1.0_r8/),                                                          &
      dbhs=(/15.0_r8, 1.0_r8/),                                                          &
      densities=(/0.015_r8, 0.015_r8/),                                                  &
      pft_ids=(/5, 14/),                                                                 &
      canopy_layers=(/1, 2/))
      
    call this%AddPatch(patch_id=4, patch_name='grassland', area=500.0_r8,                &
      ages=(/1.0_r8, 2.0_r8/),                                                           &
      dbhs=(/1.0_r8, 1.0_r8/),                                                           &
      densities=(/0.015_r8, 0.015_r8/),                                                  &
      pft_ids=(/13, 13/),                                                                &
      canopy_layers=(/1, 1/))
      
    call this%AddPatch(patch_id=5, patch_name='temperate', area=500.0_r8,                &
      ages=(/80.0_r8, 50.0_r8, 20.0_r8, 5.0_r8/),                                        &
      dbhs=(/50.0_r8, 30.0_r8, 15.0_r8, 3.0_r8/),                                        &
      densities=(/0.005_r8, 0.01_r8, 0.015_r8, 0.005_r8/),                               &
      pft_ids=(/6, 2, 2, 9/),                                                            &
      canopy_layers=(/1, 1, 2, 2/))
    
  end subroutine GetSyntheticPatchData
  
  ! --------------------------------------------------------------------------------------
  
end module SyntheticPatchTypes
