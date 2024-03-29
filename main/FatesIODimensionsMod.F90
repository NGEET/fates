module FatesIODimensionsMod

  use FatesConstantsMod, only : fates_short_string_length
  
  implicit none
  private

    ! The following dimension names must be replicated in
    ! CLM/ALMs histFileMod.F90 and 
    character(*), parameter, public :: levcapf = 'fates_levcapf'      ! matches histFileMod
    character(*), parameter, public :: levcacls = 'fates_levcacls'    ! matches histFileMod
   
    character(*), parameter, public  :: cohort = 'cohort'           ! matches clm_varcon
    character(*), parameter, public  :: column = 'column'           ! matches clm_varcon
    character(*), parameter, public  :: levsoil = 'levsoi'          ! matches clm_varcon
    character(*), parameter, public  :: levscag = 'fates_levscag'      ! matches histFileMod
    character(*), parameter, public  :: levscagpft = 'fates_levscagpf' ! matches histFileMod
    character(*), parameter, public  :: levagepft = 'fates_levagepft'  ! matches histFileMod
    character(*), parameter, public  :: levscpf = 'fates_levscpf'      ! matches histFileMod
    character(*), parameter, public  :: levscls = 'fates_levscls'      ! matches histFileMod
    character(*), parameter, public  :: levpft = 'fates_levpft'        ! matches histFileMod
    character(*), parameter, public  :: levage = 'fates_levage'        ! matches histFileMod
    character(*), parameter, public  :: levheight = 'fates_levheight'  ! matches histFileMod
    character(*), parameter, public  :: levfuel = 'fates_levfuel'      ! matches histFileMod
    character(*), parameter, public  :: levcwdsc = 'fates_levcwdsc'    ! matches histFileMod
    character(*), parameter, public  :: levcan = 'fates_levcan'        ! matches histFileMod
    character(*), parameter, public  :: levcnlf = 'fates_levcnlf'      ! matches histFileMod
    character(*), parameter, public  :: levcnlfpft = 'fates_levcnlfpf' ! matches histFileMod
    character(*), parameter, public  :: levclscpf = 'fates_levclscpf'   
    character(*), parameter, public  :: levcdsc = 'fates_levcdsc' ! matches histFileMod
    character(*), parameter, public  :: levcdpf = 'fates_levcdpf' ! matches histFileMod
    character(*), parameter, public  :: levcdam = 'fates_levcdam' ! matches histFileMod
    character(*), parameter, public  :: levagefuel = 'fates_levagefuel' ! matches histFileMod
    character(*), parameter, public  :: levelem =  'fates_levelem'
    character(*), parameter, public  :: levelpft = 'fates_levelpft'
    character(*), parameter, public  :: levelcwd = 'fates_levelcwd'
    character(*), parameter, public  :: levelage = 'fates_levelage'
    character(*), parameter, public  :: levlanduse = 'fates_levlanduse'
    character(*), parameter, public  :: levlulu = 'fates_levlulu'
    character(*), parameter, public  :: levlupft = 'fates_levlupft'
    
    ! column = This is a structure that records where FATES column boundaries
    ! on each thread point to in the host IO array, this structure
    ! is allocated by number of threads

    ! levsoil = This is a structure that records the boundaries for the
    ! soil level (includes rock) dimension

    ! levscpf = This is a structure that records the boundaries for the
    ! number of size-class x pft dimension

    ! levcapf = This is a structure that records the boundaries for the
    ! number of cohort-age-class x pft dimension

    ! levscls = This is a structure that records the boundaries for the
    ! number of size-class dimension

    ! levcacls = This is a structure that records the boundaries for the 
    ! number of cohort age class dimension

    ! levpft = This is a structure that records the boundaries for the
    ! number of pft dimension

    ! levage = This is a structure that records the boundaries for the
    ! number of patch-age-class dimension

    ! levheight = This is a structure that records the boundaries for the
    ! number of height dimension

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

    ! levcdsc = This is a structure that records the boundaries for the
    ! number of crown damage x size classes dimension

    ! levcdpf = This is a structure that records the boundaries for the
    ! number of crown damage x size classes x pft dimension

    ! levcdam = This is the structure that records the boundaries for the
    ! number of crown damage classes dimension
    
    ! levscag = This is a strcture that records the boundaries for the 
    ! number of size-classes x patch age

    ! levscagpft = This is a strcture that records the boundaries for the 
    ! number of size-classes x patch age x pft

    ! levagepft = This is a strcture that records the boundaries for the 
    ! number of patch age x pft

    ! levagefuel = This is a strcture that records the boundaries for the 
    ! number of patch age x fuel size class

    ! levclscpf = '' number of canopy layers x pft x size class
    
    ! levelem  = This records the boundaries for the number of elements
    ! levelpft = This records the boundaries for elements x pft
    ! levelcwd = This records the boundaries for element x cwd
    ! levelage = This records the boundaries for element x age

    ! levcdsc = This is a structure that records the boundaries for the
    ! number of crown damage x size classes dimension

    ! levcdpf = This is a structure that records the boundaries for the
    ! number of crown damage x size classes x pft dimension

    ! levcdam = This is the structure that records the boundaries for the
    ! number of crown damage classes dimension

    ! levlanduse = this is the structure that records the boundaries for the
    ! land use class dimension

    ! levlulu = this is the structure that records the boundaries for the
    ! (land use class) x (land use class) dimension

    ! levlupft = this is the structure that records the boundaries for the
    ! (land use class) x pft dimension

    type, public :: fates_bounds_type
       integer :: cohort_begin
       integer :: cohort_end
       integer :: column_begin          ! FATES does not have a "column" type
       integer :: column_end            ! we call this a "site" (rgk 11-2016)
       integer :: soil_begin
       integer :: soil_end
       integer :: sizeage_class_begin
       integer :: sizeage_class_end
       integer :: sizeagepft_class_begin
       integer :: sizeagepft_class_end
       integer :: agepft_class_begin
       integer :: agepft_class_end
       integer :: sizepft_class_begin
       integer :: sizepft_class_end
       integer :: coagepf_class_begin
       integer :: coagepf_class_end
       integer :: size_class_begin
       integer :: size_class_end
       integer :: coage_class_begin
       integer :: coage_class_end
       integer :: pft_class_begin
       integer :: pft_class_end
       integer :: age_class_begin
       integer :: age_class_end
       integer :: height_begin
       integer :: height_end
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
       integer :: cdsc_begin
       integer :: cdsc_end
       integer :: cdpf_begin
       integer :: cdpf_end
       integer :: cdam_begin
       integer :: cdam_end
       integer :: elem_begin
       integer :: elem_end
       integer :: elpft_begin
       integer :: elpft_end
       integer :: elcwd_begin
       integer :: elcwd_end
       integer :: elage_begin
       integer :: elage_end
       integer :: agefuel_begin
       integer :: agefuel_end
       integer :: clscpf_begin
       integer :: clscpf_end
       integer :: landuse_begin
       integer :: landuse_end
       integer :: lulu_begin
       integer :: lulu_end
       integer :: lupft_begin
       integer :: lupft_end
    end type fates_bounds_type
    


  ! This structure is not allocated by thread, but the upper and lower boundaries
  ! of the dimension for each thread is saved in the clump_ entry
  type, public :: fates_io_dimension_type
     character(len=fates_short_string_length) :: name
     integer :: lower_bound
     integer :: upper_bound
     integer, allocatable :: clump_lower_bound(:) ! lower bound of thread's portion of HIO array
     integer, allocatable :: clump_upper_bound(:) ! upper bound of thread's portion of HIO array
   contains
     procedure :: Init
     procedure :: SetThreadBounds
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
