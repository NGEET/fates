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
      this%rank = rank(data)

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
      this%rank = rank(data) 
      
    end subroutine RegisterInterfaceVariable_2d

  ! ====================================================================================
    
    subroutine UpdateInterfaceVariable(this, var, subgrid_indices)
      
      class(fates_interface_variable_type) :: this             ! variable being updated
      class(fates_interface_variable_type), intent(in) :: var  ! variable update source
      integer, intent(in) :: subgrid_indices(:)                ! subgrid indices for the update source

      class(*), pointer :: data_var0d        => null()
      class(*), pointer :: data_var1d(:)     => null()
      class(*), pointer :: data_var2d(:,:)   => null()
      class(*), pointer :: data_var3d(:,:,:) => null()

      integer :: index  ! index for the subgrid level of the input interface variable 

      ! This update method assumes that the first rank of the HLM data arrays
      ! corresponds to the subgrid level of the interface variable type.
      ! E.g. col_cf%w_scalar(c,1:nlevsoil) shows that the first rank is the column index.
      ! TODO: This should be held in an interface requirements document.

      ! Get the subgrid index for the updating variable 
      index = subgrid_indices(var%subgrid_index)

      ! Update the data pointer based on the rank of the source variable while indexing
      ! into the appropriate subgrid level
      ! TODO: This assumes HLM->FATES direction; Validate this for FATES->HLM direction
      select case (var%rank)
        case(0)
          data_var0d => var%data0d
        case(1)
          data_var0d => var%data1d(index)
        case(2)
          data_var1d => var%data2d(index,:)
        case(3)
          data_var2d => var%data3d(index,:,:)
        case default
          call endrun(fates_log, 'FATES ERROR: Unsupported interface variable source rank in UpdateInterfaceVariable')
      end select

      ! Update the data pointer of the target variable based on its rank
      select case (this%rank)
        case(0)
          this%data0d = data_var0d
        case(1)
          this%data1d = data_var1d
        case(2)
          this%data2d = data_var2d
        case default
          call endrun(fates_log, 'FATES ERROR: Unsupported interface variable input rank in UpdateInterfaceVariable')
      end select

    end subroutine UpdateInterfaceVariable

  ! ====================================================================================

end module FatesInterfaceVariableTypeMod