module FatesHistoryVariableKindMod

  use FatesConstantsMod, only : fates_long_string_length
  use FatesGlobals, only : fates_log
  use FatesHistoryDimensionMod, only : fates_history_dimension_type

  implicit none

  ! NOTE(RGK, 2016) %active is not used yet. Was intended as a check on the HLM->FATES
  ! control parameter passing to ensure all active dimension types received all
  ! dimensioning specifications from the host, but we currently arent using those
  ! passing functions..

  ! This structure is not multi-threaded
  type fates_history_variable_kind_type
     character(len=fates_long_string_length) :: name ! String labelling this IO type
     integer              :: ndims       ! number of dimensions in this IO type
     integer, allocatable :: dimsize(:)  ! The size of each dimension
     logical, private :: active_
     integer :: dim1_index
     integer :: dim2_index

   contains

     procedure, public :: Init
     procedure, public :: set_active
     procedure, public :: is_active

  end type fates_history_variable_kind_type



contains

  ! ===================================================================================
  subroutine Init(this, name, num_dims)

    use FatesConstantsMod, only : fates_unset_int
    
    implicit none

    class(fates_history_variable_kind_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: num_dims
    
    this%name = trim(name)
    this%ndims = num_dims
    allocate(this%dimsize(this%ndims))
    this%dimsize(:) = fates_unset_int
    this%active_ = .false.
    this%dim1_index = fates_unset_int
    this%dim2_index = fates_unset_int
    
  end subroutine Init
  
 ! =======================================================================
 subroutine set_active(this)
   implicit none
   class(fates_history_variable_kind_type), intent(inout) :: this
   this%active_ = .true.
 end subroutine set_active

 logical function is_active(this)
   implicit none
   class(fates_history_variable_kind_type), intent(in) :: this
   is_active = this%active_
 end function is_active

  ! ====================================================================================

  function iotype_index(iotype_name, num_dim_kinds, dim_kinds) result(dk_index)

    ! argument
    character(len=*), intent(in) :: iotype_name
    integer, intent(in) :: num_dim_kinds
    type(fates_history_variable_kind_type), intent(in) :: dim_kinds(:)

    ! local
    integer :: dk_index

    do dk_index=1, num_dim_kinds
       if (trim(iotype_name) .eq. trim(dim_kinds(dk_index)%name)) then
          return
       end if
    end do
    write(fates_log(),*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    !end_run

  end function iotype_index
  
end module FatesHistoryVariableKindMod
