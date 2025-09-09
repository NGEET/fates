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

  ! Interface variable registry type
  type, public :: fates_interface_variable_type
    
    character(len=48) :: key            ! common registry key
    class(*), pointer :: data0d         ! scalar polymorphic data pointer
    class(*), pointer :: data1d(:)      ! 1D polymorphic data pointer
    class(*), pointer :: data2d(:,:)    ! 2D polymorphic data pointer
    class(*), pointer :: data3d(:,:,:)  ! 3D polymorphic data pointer
    logical           :: active         ! true if the variable is used by the host land model

    contains
      procedure :: Initialize => InitializeInterfaceVariable
      generic   :: Register => RegisterInterfaceVariable_1d, RegisterInterfaceVariable_2d

      procedure, private :: RegisterInterfaceVariable_1d
      procedure, private :: RegisterInterfaceVariable_2d
      
  end type fates_interface_variable_type
  
  contains
  
  ! ====================================================================================
  
    subroutine InitializeInterfaceVariable(this, key)
                                            
      class(fates_interface_variable_type) :: this

      character(len=*), intent(in) :: key
      
      this%data0d  => null()
      this%data1d  => null()
      this%data2d  => null()
      this%data3d  => null()
      this%key = key 
      this%active = .false.
      
    end subroutine InitializeInterfaceVariable

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_1d(this, data, active)
      
      class(fates_interface_variable_type) :: this

      class(*), target, intent(in) :: data(:)
      logical, intent(in)          :: active
      
      this%data1d => data(:)
      this%active = active
      
    end subroutine RegisterInterfaceVariable_1d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_2d(this, data, active)
      
      class(fates_interface_variable_type) :: this

      class(*), target, intent(in) :: data(:,:)
      logical, intent(in)          :: active
      
      this%data2d => data(:,:)
      this%active = active
      
    end subroutine RegisterInterfaceVariable_2d

  ! ====================================================================================

end module FatesInterfaceVariableTypeMod