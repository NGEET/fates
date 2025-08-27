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
      procedure :: InitializeInterfaceVariable => Init
      procedure :: Register => RegisterInterfaceVariable_int_scalar
      
  end type fates_interface_variable_type
  
  contains
  
  ! ====================================================================================
  
    subroutine InitializeInterfaceVariable(this, variable_name)
                                            
      class(fates_interface_variable_type) :: this

      character(len=*), intent(in) :: variable_name  
      
      nullify(this%iscalar)
      nullify(this%int1d)
      nullify(this%int2d)
      nullify(this%int3d)
      nullify(this%rscalar)
      nullify(this%r81d)
      nullify(this%r82d)
      nullify(this%r83d)
      
      this%variable_name = variable_name
      this%active = .false.
      
    end subroutine InitializeInterfaceVariable

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable(this, data, active)
      
      class(fates_interface_variable_type) :: this

      class(*), pointer, intent(in) :: data
      logical, intent(in)           :: active
      
      this%data => data
      this%active = active
      
    end subroutine RegisterInterfaceVariable

end module FatesInterfaceVariableTypeMod