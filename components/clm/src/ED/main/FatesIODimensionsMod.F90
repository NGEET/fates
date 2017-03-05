module FatesIODimensionsMod

  use FatesConstantsMod, only : fates_short_string_length
  
  implicit none

    character(*), parameter :: cohort = 'cohort'
    character(*), parameter :: patch = 'patch'
    character(*), parameter :: column = 'column'
    character(*), parameter :: levgrnd = 'levgrnd'
    character(*), parameter :: levscpf = 'levscpf'
    character(*), parameter :: levscls = 'levscls'
    character(*), parameter :: levpft = 'levpft'
    character(*), parameter :: levage = 'levage'
    character(*), parameter :: levfuel = 'levfuel'
    character(*), parameter :: levcwdsc = 'levcwdsc'
    character(*), parameter :: levcan = 'levcan'
    character(*), parameter :: levcnlf = 'levcnlf'
    character(*), parameter :: levcnlfpft = 'lvcnlfpf'

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

    ! levscls = This is a structure that records the boundaries for the
    ! number of size-class dimension

    ! levpft = This is a structure that records the boundaries for the
    ! number of pft dimension

    ! levage = This is a structure that records the boundaries for the
    ! number of patch-age-class dimension

    ! levfuel = This is a structure that records the boundaries for the
    ! number of fuel-size-class dimension

    ! levcwdsc = This is a structure that records the boundaries for the
    ! number of coarse-woody-debris-size-class dimension

    ! levcan = This is a structure that records the boundaries for the
    ! number of canopy layer dimension

    ! levcnlf = This is a structure that records the boundaries for the
    ! number of cnanopy layer x leaf layer dimension

    ! levcnlfpft = This is a structure that records the boundaries for the
    ! number of canopy layer x leaf layer x pft dimension


    type, public :: fates_bounds_type
       integer :: patch_begin
       integer :: patch_end
       integer :: cohort_begin
       integer :: cohort_end
       integer :: column_begin          ! FATES does not have a "column" type
       integer :: column_end            ! we call this a "site" (rgk 11-2016)
       integer :: ground_begin
       integer :: ground_end
       integer :: sizepft_class_begin
       integer :: sizepft_class_end
       integer :: size_class_begin
       integer :: size_class_end
       integer :: pft_class_begin
       integer :: pft_class_end
       integer :: age_class_begin
       integer :: age_class_end
       integer :: fuel_begin
       integer :: fuel_end
       integer :: cwdsc_begin
       integer :: cwdsc_end
       integer :: can_begin
       integer :: can_end
       integer :: cnlf_begin
       integer :: cnlf_end
       integer :: cnlfpft_begin
       integer :: cnlfpft_end
    end type fates_bounds_type
    


  ! This structure is not allocated by thread, but the upper and lower boundaries
  ! of the dimension for each thread is saved in the clump_ entry
  type fates_io_dimension_type
     character(len=fates_short_string_length) :: name
     integer :: lower_bound
     integer :: upper_bound
     integer, allocatable :: clump_lower_bound(:) ! lower bound of thread's portion of HIO array
     integer, allocatable :: clump_upper_bound(:) ! upper bound of thread's portion of HIO array
   contains
     procedure, public :: Init
     procedure, public :: SetThreadBounds
  end type fates_io_dimension_type

contains

  ! =====================================================================================
  subroutine Init(this, name, num_threads, lower_bound, upper_bound)

    implicit none

    ! arguments
    class(fates_io_dimension_type), intent(inout) :: this
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

  end subroutine Init

  ! =====================================================================================

  subroutine SetThreadBounds(this, thread_index, lower_bound, upper_bound)

    implicit none

    class(fates_io_dimension_type), intent(inout) :: this
    integer, intent(in) :: thread_index
    integer, intent(in) :: lower_bound
    integer, intent(in) :: upper_bound

    this%clump_lower_bound(thread_index) = lower_bound
    this%clump_upper_bound(thread_index) = upper_bound

  end subroutine SetThreadBounds
  
end module FatesIODimensionsMod
