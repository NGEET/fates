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

  integer, parameter, public :: subgrid_gridcell = 5
  integer, parameter, public :: subgrid_topounit = 4
  integer, parameter, public :: subgrid_landunit = 3
  integer, parameter, public :: subgrid_column   = 2
  integer, parameter, public :: subgrid_patch    = 1

  ! Interface registry variable type
  type, public :: fates_interface_variable_type
    
    character(len=48) :: key            ! common registry key
    class(*), pointer :: data0d         ! scalar polymorphic data pointer
    class(*), pointer :: data1d(:)      ! 1D polymorphic data pointer
    class(*), pointer :: data2d(:,:)    ! 2D polymorphic data pointer
    class(*), pointer :: data3d(:,:,:)  ! 3D polymorphic data pointer
    logical           :: active         ! true if the variable is used by the host land model
    integer           :: rank           ! rank of the variable (0, 1, 2, or 3)
    integer           :: subgrid        ! subgrid level (0 = gridcell, 1 = landunit, 2 = column, 3 = patch)

    contains
      procedure :: Initialize => InitializeInterfaceVariable
      procedure :: Update     => UpdateInterfaceVariable

      generic :: Register => RegisterInterfaceVariable_1d, RegisterInterfaceVariable_2d

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
    
    subroutine RegisterInterfaceVariable_1d(this, data, active, subgrid_index)
      
      class(fates_interface_variable_type) :: this

      class(*), target, intent(in) :: data(:)
      logical, intent(in)          :: active
      integer, intent(in)          :: subgrid_index

      this%data1d => data(:)
      this%active = active
      this%subgrid = subgrid_index

    end subroutine RegisterInterfaceVariable_1d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_2d(this, data, active, subgrid_index)
      
      class(fates_interface_variable_type) :: this

      class(*), target, intent(in)  :: data(:,:)
      logical, intent(in)           :: active
      integer, intent(in)           :: subgrid_index

      this%data2d => data(:,:)
      this%active = active
      this%subgrid = subgrid_index
      
    end subroutine RegisterInterfaceVariable_2d

  ! ====================================================================================
    
    subroutine UpdateInterfaceVariable(this, var, subgrid_indices)
      
      class(fates_interface_variable_type) :: this

      class(fates_interface_variable_type), intent(in) :: var
      integer, intent(in) :: subgrid_indices(:)

      ! This update method assumes that the first rank of the HLM data arrays
      ! corresponds to the subgrid level of the interface variable type.
      ! E.g. col_cf%w_scalar(c,1:nlevsoil) shows that the first rank is the column index.
      ! TODO: This should be held in an interface requirements document.
      


      ! TODO: add column index to the interface variable type to allow
      ! for appropriate slicing of input pointer array
      ! e.g.
      ! if (this%rank == 1)) then
      !   data_this => this%data1d
      ! else if (this%rank == 2) then
      !   data_this => this%data2d
      ! end if
      ! if (var%rank == 1)) then
      !   data_var => var%data1d
      ! else if (this%rank == 2) then
      !   data_var => var%data2d(this%col,:)
      ! end if
      ! data_this = data_var()
      ! This isn't exactly right, but you get the idea


    end subroutine UpdateInterfaceVariable

  ! ====================================================================================

end module FatesInterfaceVariableTypeMod