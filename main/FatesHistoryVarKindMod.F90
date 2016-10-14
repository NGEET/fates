module FatesHistoryVariableKindMod

  use FatesConstantsMod, only : fates_long_string_length
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
     logical              :: active
     type(fates_history_dimension_type), pointer :: dim1_ptr
     type(fates_history_dimension_type), pointer :: dim2_ptr

   contains

     procedure, public :: Init => InitVariableKind

  end type fates_history_variable_kind_type



contains

  ! ===================================================================================
  subroutine InitVariableKind(this, name, num_dims)

    use FatesConstantsMod, only : fates_unset_int
    
    implicit none

    class(fates_history_variable_kind_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: num_dims
    
    this%name = trim(name)
    this%ndims = num_dims
    allocate(this%dimsize(this%ndims))
    this%dimsize(:) = fates_unset_int
    this%active = .false.
    nullify(this%dim1_ptr)
    nullify(this%dim2_ptr)
    
  end subroutine InitVariableKind
  
  
  
end module FatesHistoryVariableKindMod
