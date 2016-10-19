module FatesHistoryDimensionMod

  use FatesConstantsMod, only : fates_short_string_length
  
  implicit none

    ! FIXME(bja, 2016-10) do these need to be strings, or can they be integer enumerations?
    character(*), parameter :: patch_r8 = 'PA_R8'
    character(*), parameter :: patch_ground_r8 = 'PA_GRND_R8'
    character(*), parameter :: patch_class_pft_r8 = 'PA_SCPF_R8'
    character(*), parameter :: site_r8 = 'SI_R8'
    character(*), parameter :: site_ground_r8 = 'SI_GRND_R8'
    character(*), parameter :: site_class_pft_r8 = 'SI_SCPF_R8'
    character(*), parameter :: patch_int = 'PA_INT'

    integer, parameter :: fates_num_dimension_types = 4
    character(*), parameter :: patch = 'patch'
    character(*), parameter :: column = 'column'
    character(*), parameter :: levgrnd = 'levgrnd'
    character(*), parameter :: levscpf = 'levscpf'

    ! patch = This is a structure that records where FATES patch boundaries
    ! on each thread point to in the host IO array, this structure
    ! is allocated by number of threads

    ! column = This is a structure that records where FATES column boundaries
    ! on each thread point to in the host IO array, this structure
    ! is allocated by number of threads

    ! ground = This is a structure that records the boundaries for the
    ! ground level (includes rock) dimension

    ! levscpf = This is a structure that records the boundaries for the
    ! number of size-class x pft dimension


  ! This structure is not allocated by thread, but the upper and lower boundaries
  ! of the dimension for each thread is saved in the clump_ entry
  type fates_history_dimension_type
     character(len=fates_short_string_length) :: name
     integer :: lower_bound
     integer :: upper_bound
     integer, allocatable :: clump_lower_bound(:) ! lower bound of thread's portion of HIO array
     integer, allocatable :: clump_upper_bound(:) ! upper bound of thread's portion of HIO array
   contains
     procedure, public :: Init => InitHistoryDimensions
     procedure, public :: SetThreadBounds => SetHistoryDimensionThreadBounds
  end type fates_history_dimension_type

contains

  ! =====================================================================================
  subroutine InitHistoryDimensions(this, name, num_threads, lower_bound, upper_bound)

    implicit none

    ! arguments
    class(fates_history_dimension_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: num_threads
    integer, intent(in) :: lower_bound
    integer, intent(in) :: upper_bound

    this%name = trim(name)
    this%lower_bound = lower_bound
    this%upper_bound = upper_bound

    allocate(this%clump_lower_bound(num_threads))
    this%clump_lower_bound(:) = -1

    allocate(this%clump_upper_bound(num_threads))
    this%clump_upper_bound(:) = -1

  end subroutine InitHistoryDimensions

  ! =====================================================================================

  subroutine SetHistoryDimensionThreadBounds(this, thread_index, lower_bound, upper_bound)

    implicit none

    class(fates_history_dimension_type), intent(inout) :: this
    integer, intent(in) :: thread_index
    integer, intent(in) :: lower_bound
    integer, intent(in) :: upper_bound

    this%clump_lower_bound(thread_index) = lower_bound
    this%clump_upper_bound(thread_index) = upper_bound

  end subroutine SetHistoryDimensionThreadBounds
  
end module FatesHistoryDimensionMod
