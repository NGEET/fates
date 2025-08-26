module FatesInterfaceVariableTypeMod

  ! This module contains the type definition and associated type-bound procedures
  ! used to create an indexed list of associated HLM and FATES variables that are
  ! related across the application programming interface.
  ! This method is largely inspired by the FATES history infrastructure
  
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_long_string_length

  implicit none
  private
  

  type, public :: fates_interface_variable_type
    
    character(len=48)                       :: variable_name  ! variable common reference name
    character(len=fates_long_string_length) :: description    ! variable description
    
    logical :: active = .false.                               ! true if the variable is used by the host land model
    integer :: update_frequency                               ! this should facilitate a check that the update is being called in the correct subroutine 
   
    ! pointers to data (only one of these to be allocated per variable)
    ! TODO: make this polymorphic?
    integer, pointer  :: intscalar 
    integer, pointer  :: int1d(:)
    integer, pointer  :: int2d(:,:)
    integer, pointer  :: int3d(:,:,:)
    real(r8), pointer :: r8scalar
    real(r8), pointer :: r81d(:)
    real(r8), pointer :: r82d(:,:)
    real(r8), pointer :: r83d(:,:,:)
    
    contains
      procedure :: InitializeInterfaceVariables => Init
  end type fates_interface_variable_type
  
  contains
  
  subroutine InitializeInterfaceVariables(this, variable_name, description, active, &
                                          update_frequency)
                                          
    class(fates_interface_variable_type) :: this

    character(len=*), intent(in) :: variable_name  
    character(len=*), intent(in) :: description    
    logical, intent(in)          :: active
    integer, intent(in)          :: update_frequency
    
    this%variable_name = variable_name
    this%description = description
    this%active = active
    this%update_frequency = update_frequency
    
  end subroutine InitializeInterfaceVariables

end module FatesInterfaceVariableTypeMod