module FatesInterfaceVariableTypeMod

  ! This module contains the type definition and associated type-bound procedures
  ! used to create an indexed list of associated HLM and FATES variables that are
  ! related across the HLM-FATES application programming interface (API).
  ! This method is largely inspired by the FATES history infrastructure

  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
 
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_long_string_length
  use FatesConstantsMod, only : fates_unset_int

  implicit none
  private

  ! Interface registry variable type
  ! This defined type holds the pointer to HLM or FATES variable that is
  ! associated with a common "key" name (defined in the interface parameter types).
  ! It also includes a set of data characterizing the type of data not defined by
  ! type or rank (e.g. HLM subgrid association, update frequency, etc).

  type, public :: fates_interface_variable_type
    
    character(len=48) :: key            ! common registry key
    class(*), pointer :: data0d         ! scalar polymorphic data pointer
    class(*), pointer :: data1d(:)      ! 1D polymorphic data pointer
    class(*), pointer :: data2d(:,:)    ! 2D polymorphic data pointer
    class(*), pointer :: data3d(:,:,:)  ! 3D polymorphic data pointer
    logical           :: active         ! true if the variable is used by the host land model
    logical           :: accumulate     ! If true, this variable should add the source to the target
    logical           :: zero_first     ! If true, zero the target variable before accumulation
    logical           :: is_last        ! True if the variable is associated with the last patch for the associated subgrid unit
    integer           :: bc_dir         ! 0 if bc_in, 1 if bc_out
    integer           :: data_rank      ! rank of the variable (0, 1, 2, or 3)
    integer           :: update_frequency ! frequency of updates 
    integer           :: subgrid_type   ! subgrid integer ID associated with this variable
    real(r8)          :: conversion_factor ! conversion factor to adjust units as necessary
    integer, allocatable :: data_size(:)   ! size of the first dimension of the variable

    contains
      procedure :: CheckBounds
      procedure :: Convert    => ConvertScaleInterfaceVariable
      procedure :: Initialize => InitializeInterfaceVariable
      procedure :: IsSubgridType
      procedure :: Normalize  => NormalizeLitterVariable
      procedure :: Update     => UpdateInterfaceVariable
      procedure :: Dump
      procedure :: SetLastState

      generic :: Register => RegisterInterfaceVariable_0d, &
                             RegisterInterfaceVariable_1d, &
                             RegisterInterfaceVariable_2d
      procedure, private :: RegisterInterfaceVariable_0d
      procedure, private :: RegisterInterfaceVariable_1d
      procedure, private :: RegisterInterfaceVariable_2d

      procedure, private :: CompareRegistryVariableSizes
      
  end type fates_interface_variable_type
  
  contains
  
  ! ====================================================================================
  
  subroutine CheckBounds(this, var)
  
    ! This procedure checks that the data bounds of the calling variable are consistent
    ! with the data bounds of the input argument variable

    class(fates_interface_variable_type), intent(in) :: this
    class(fates_interface_variable_type), intent(in) :: var

    ! Locals
    integer :: ub1, ub2
    integer :: lb1, lb2
    logical :: bounds_mismatch

    bounds_mismatch = .false.

    if (this%data_rank >= 1) then
      ub1 = ubound(this%data1d, dim=1)
      lb1 = lbound(this%data1d, dim=1)
      ub2 = ubound(var%data1d, dim=1)
      lb2 = lbound(var%data1d, dim=1)
      if (ub1 /= ub2 .or. lb1 /= lb2) then
        write(fates_log(),*) 'Dimension 1 bounds mismatch for variable: ', this%key
        write(fates_log(),*) 'Upper bounds: ', ub1, ', ', ub2
        write(fates_log(),*) 'Lower bounds: ', lb1, ', ', lb2
        bounds_mismatch = .true.
      end if
    else if (this%data_rank >= 2) then
      ub1 = ubound(this%data2d, dim=2)
      lb1 = lbound(this%data2d, dim=2)
      ub2 = ubound(var%data2d, dim=2)
      lb2 = lbound(var%data2d, dim=2)
      if (ub1 /= ub2 .or. lb1 /= lb2) then
        write(fates_log(),*) 'Dimension 2 bounds mismatch for variable: ', this%key
        write(fates_log(),*) 'Upper bounds: ', ub1, ', ', ub2
        write(fates_log(),*) 'Lower bounds: ', lb1, ', ', lb2
        bounds_mismatch = .true.
      end if
    end if

    if (bounds_mismatch) then
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine CheckBounds
  
  ! ====================================================================================
  
  subroutine ConvertScaleInterfaceVariable(this, scalar)
  
    ! This subroutine scales a given interface variable data.  This
    ! routine is provided primarily to allow for post accumulation conversion
    ! of units and scaling of data as necessary.  Note that this only provides
    ! scaling for real data types currently.

    class(fates_interface_variable_type), intent(inout) :: this
    real(r8), intent(in)                               :: scalar

    ! Locals
    integer :: i
    character(len=fates_long_string_length) :: err_msg = 'FATES Error: invalid interface type being scaled'

    select case (this%data_rank)
      case(0)
        select type(var => this%data0d)
          type is (real(r8))
            var = var * scalar
          class default
            write(fates_log(),*), err_msg
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case(1)
        select type(var => this%data1d)
          type is (real(r8))
            var = var * scalar
          class default
            write(fates_log(),*), err_msg
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case(2)
        select type(var => this%data2d)
          type is (real(r8))
            var = var * scalar
          class default
            write(fates_log(),*), err_msg
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case default
        write(fates_log(),*) 'FATES ERROR: Unsupported interface variable rank in ConvertScaleInterfaceVariable'
        call endrun(msg=errMsg(__FILE__, __LINE__))
      end select

  end subroutine ConvertScaleInterfaceVariable
  
  ! ====================================================================================
  
  subroutine Dump(this)
  
    ! This procedure will print output to the log file.  Developers should use this 
    ! procedure for diagnostic purposes when attempt to review interface variable
    ! data as it includes select type statements to resolve the polymorphic data types.

    class(fates_interface_variable_type), intent(in) :: this

    write(fates_log(),*) 'FATES Interface Variable Dump:'
    write(fates_log(),*) '  Key: ', this%key
    write(fates_log(),*) '  Active: ', this%active
    write(fates_log(),*) '  Accumulate: ', this%accumulate
    write(fates_log(),*) '  Data Rank: ', this%data_rank
    write(fates_log(),*) '  Data Size: ', this%data_size

    select case (this%data_rank)
      case(0)
        select type(var => this%data0d)
          type is (real(r8))
            write(fates_log(),*) '  Data (real): ', var
          type is (integer)
            write(fates_log(),*) '  Data (integer): ', var
          class default
            write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case(1)
        select type(var => this%data1d)
          type is (real(r8))
            write(fates_log(),*) '  Data (real): ', var
          type is (integer)
            write(fates_log(),*) '  Data (integer): ', var
          class default
            write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case(2)
        select type(var => this%data2d)
          type is (real(r8))
            write(fates_log(),*) '  Data (real): ', var
          type is (integer)
            write(fates_log(),*) '  Data (integer): ', var
          class default
            write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
            call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      case default
        write(fates_log(),*) 'FATES ERROR: Unsupported interface variable rank in Dump'
        call endrun(msg=errMsg(__FILE__, __LINE__))
      end select

  end subroutine Dump

  ! ====================================================================================
  
    subroutine InitializeInterfaceVariable(this, key, update_frequency, bc_dir)
    
      ! This procedure initializes an interface variable assigning its key
      ! update frequency and boundary condition direction.  It also unsets
      ! all other variables including the data pointers.
                                            
      class(fates_interface_variable_type), intent(inout) :: this
      character(len=*), intent(in) :: key
      integer, intent(in)          :: update_frequency
      integer, intent(in)          :: bc_dir

      
      allocate(this%data_size(3))

      ! Initialize components that are set later
      this%data_size = fates_unset_int
      this%data_rank = fates_unset_int
      this%data0d  => null()
      this%data1d  => null()
      this%data2d  => null()
      this%data3d  => null()
      this%active = .false.
      this%accumulate = .false.
      this%zero_first = .false.
      this%is_last = .false.
      this%conversion_factor = nan
      this%subgrid_type = fates_unset_int

      ! Initialize registry variable components that are updated at variable definition
      this%key = key 
      this%update_frequency = update_frequency
      this%bc_dir = bc_dir

    end subroutine InitializeInterfaceVariable
    
  ! ====================================================================================

  logical function IsSubgridType(this, subgrid_type)
  
    ! This procedure will check if the subgrid type of the interface variable
    ! matches the input argument subgrid type

    class(fates_interface_variable_type), intent(in) :: this
    integer, intent(in) :: subgrid_type

    IsSubgridType = (this%subgrid_type == subgrid_type)

  end function IsSubgridType

  ! ====================================================================================
  
    subroutine NormalizeLitterVariable(this, norm_var)
    
      ! This subroutine is specifically set up to normalize the litter flux interface
      ! variables by the decomposition thickness variable passed in as norm_factor.
      ! It could be extended later to handle other normalization needs as necessary.

      class(fates_interface_variable_type), intent(inout) :: this
      class(fates_interface_variable_type), intent(in)    :: norm_var  ! normalization interface variable

      ! Locals
      integer :: i
      character(len=fates_long_string_length) :: err_msg = 'FATES Error: invalid interface type being normalized'

      select case (this%data_rank)
        case(1)
          select type(var => this%data1d)
            type is (real(r8))
              select type(norm => norm_var%data1d)
                type is (real(r8))
                  do i = lbound(var, dim=1), ubound(var, dim=1)
                    var(i) = var(i) / norm(i)
                  end do
                class default
                  write(fates_log(),*), err_msg
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            class default
              write(fates_log(),*), err_msg
              call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
        case default
          write(fates_log(),*) 'FATES ERROR: Unsupported interface variable rank in NormalizeLitterVariable'
          call endrun(msg=errMsg(__FILE__, __LINE__))
        end select

    end subroutine NormalizeLitterVariable

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_0d(this, data, active, accumulate, is_first, subgrid_type, conversion_factor)
    
      ! This procedure registers the interface variable by associating the data pointer with 
      ! the scalar input data argument.  It also sets a number of required and optional arguments
      ! for metadata associated with the variable.

      class(fates_interface_variable_type), intent(inout) :: this
      class(*), target, intent(in) :: data
      logical, intent(in)          :: active
      logical, intent(in)          :: accumulate
      logical, intent(in)          :: is_first
      integer, intent(in)          :: subgrid_type
      real(r8), intent(in)         :: conversion_factor

      this%data0d => data
      this%active = active
      this%accumulate = accumulate
      this%zero_first = is_first
      this%data_rank = rank(data)
      this%conversion_factor = conversion_factor
      this%subgrid_type = subgrid_type

    end subroutine RegisterInterfaceVariable_0d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_1d(this, data, active, accumulate, is_first, subgrid_type, conversion_factor)

      ! This procedure registers the interface variable by associating the data pointer with 
      ! the 1D array input data argument.  It also sets a number of required and optional arguments
      ! for metadata associated with the variable.

      class(fates_interface_variable_type), intent(inout) :: this

      class(*), target, intent(in) :: data(:)
      logical, intent(in)          :: active
      logical, intent(in)          :: accumulate
      logical, intent(in)          :: is_first
      integer, intent(in)          :: subgrid_type  
      real(r8), intent(in)         :: conversion_factor

      this%data1d => data(:)
      this%active = active
      this%accumulate = accumulate
      this%zero_first = is_first
      this%data_rank = rank(data)
      this%data_size(1) = size(data, dim=1)
      this%conversion_factor = conversion_factor
      this%subgrid_type = subgrid_type

    end subroutine RegisterInterfaceVariable_1d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_2d(this, data, active, accumulate, is_first, subgrid_type, conversion_factor)

      ! This procedure registers the interface variable by associating the data pointer with 
      ! the 2D array input data argument.  It also sets a number of required and optional arguments
      ! for metadata associated with the variable.

      class(fates_interface_variable_type), intent(inout) :: this
      class(*), target, intent(in)  :: data(:,:)
      logical, intent(in)           :: active
      logical, intent(in)          :: accumulate
      logical, intent(in)          :: is_first
      integer, intent(in)          :: subgrid_type
      real(r8), intent(in)         :: conversion_factor

      this%data2d => data(:,:)
      this%active = active
      this%accumulate = accumulate
      this%zero_first = is_first
      this%data_rank = rank(data) 
      this%data_size(1) = size(data, dim=1)
      this%data_size(2) = size(data, dim=2)
      this%subgrid_type = subgrid_type
      this%conversion_factor = conversion_factor

    end subroutine RegisterInterfaceVariable_2d

  ! ====================================================================================

    subroutine SetLastState(this, last_state)
    
      ! This procedure sets the last state logical metadata for the calling interface variable

      class(fates_interface_variable_type), intent(inout) :: this
      logical, intent(in) :: last_state

      this%is_last = last_state

    end subroutine SetLastState

  ! ====================================================================================
    
    subroutine UpdateInterfaceVariable(this, var, is_last)
    
      ! This is the main procedure that drives the updates between the HLM and FATES.
      ! It updates the calling interface variable data with data from the argument
      ! interface variable.  Directionality of the update is from the argument to the
      ! caller.  The bulk of the code here is type selects to handle the polymorphic
      ! data pointers in the interface variables types.  
      ! This procedure utilizes metadata stored in the interface variables to determine
      ! how updates should be handled, for example whether or not to zero the calling
      ! interface variable data prior to adding the argument variable data to it.

      ! Arguments
      class(fates_interface_variable_type), intent(inout) :: this     ! variable being updated
      class(fates_interface_variable_type), intent(in)    :: var      ! variable update source
      logical, intent(out), optional                      :: is_last  ! true if last patch for subgrid unit   

      ! Locals
      character(len=fates_long_string_length) :: msg_mismatch = 'FATES ERROR: Mismatched interface variable types'

      ! Check that the dimensions of the source and target match
      call this%CompareRegistryVariableSizes(var)

      ! If the is_last argument is present, output the is_last flag state for this variable
      if (present(is_last)) is_last = this%is_last
        
      ! Update the data of the target variable using the source variable data pointer
      ! Make sure the types match for the polymorphic data to allow for copying from the 
      ! source to the target.
      ! Note that due to the use of polymorphic pointers, we must use select type constructs
      ! to determine the actual type of the data being pointed to allowing for type-safe assignment.
      ! This currently only supports real and integer types and no conversion between types 
      ! should be performed
      select case (this%data_rank)
        case(0)
          select type(dest => this%data0d)
            type is (real(r8))
              select type(source => var%data0d)
                type is (real(r8))
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0.0_r8
                    end if
                    dest = dest + source
                  else
                    dest = source
                  end if
                  ! Apply conversion factor if this is the last patch associated with the subgrid unit
                  if (this%is_last) then
                    dest = dest * var%conversion_factor
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data0d)
                type is (integer)
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0
                    end if
                    dest = dest + source 
                  else
                    dest = source 
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            class default
              write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
              call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

        case(1)

          select type(dest => this%data1d)
            type is (real(r8))
              select type(source => var%data1d)
                type is (real(r8))
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0.0_r8
                    end if
                    dest = dest + source
                  else
                    dest = source
                  end if
                  if (this%is_last) then
                    dest = dest * var%conversion_factor
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data1d)
                type is (integer)
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0
                    end if
                    dest = dest + source
                  else
                    dest = source
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            class default
              write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
              call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

        case(2)
          
          select type(dest => this%data2d)
            type is (real(r8))
              select type(source => var%data2d)
                type is (real(r8))
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0.0_r8
                    end if
                    dest = dest + source
                  else
                    dest = source
                  end if
                  if (this%is_last) then
                    dest = dest * var%conversion_factor
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data2d)
                type is (integer)
                  if (this%accumulate) then
                    if (this%zero_first) then
                      dest = 0
                    end if
                    dest = dest + source
                  else
                    dest = source
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end select
            class default
              write(fates_log(),*), 'FATES ERROR: Unsupported interface variable type'
              call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
        case default
          write(fates_log(),*) 'FATES ERROR: Unsupported interface variable target rank in UpdateInterfaceVariable'
          call endrun(msg=errMsg(__FILE__, __LINE__))
      end select

    end subroutine UpdateInterfaceVariable

  ! ====================================================================================
    
    subroutine CompareRegistryVariableSizes(this, var) 
    
      ! This procedure checks to make sure that the ranks and size of the
      ! data associated with the interface variables is the same.
      
      class(fates_interface_variable_type), intent(in) :: this ! variable being updated
      class(fates_interface_variable_type), intent(in) :: var  ! variable update source

      if (this%data_size(1) /= var%data_size(1) .or. & 
          this%data_size(2) /= var%data_size(2) .or. &
          this%data_size(3) /= var%data_size(3)) then

        write(fates_log(),*) 'FATES ERROR: Mismatched interface variable sizes in UpdateInterfaceVariable'

        if (this%data_rank == 1) then
          write(fates_log(),*) '  Target, size: ', this%key, this%data_size(1)
          write(fates_log(),*) '  Source, size: ', var%key, var%data_size(1)
        else if (this%data_rank == 2) then
          write(fates_log(),*) '  Target, size: ', this%key, this%data_size(1), this%data_size(2)
          write(fates_log(),*) '  Source, size: ', var%key, var%data_size(1), var%data_size(2)
        else if (this%data_rank == 3) then
          write(fates_log(),*) '  Target, size: ', this%key, this%data_size(1), this%data_size(2), this%data_size(3)
          write(fates_log(),*) '  Source, size: ', var%key, var%data_size(1), var%data_size(2), var%data_size(3)
        else 
          write(fates_log(),*) '  Unsupported interface variable rank in UpdateErrorMessage'
          write(fates_log(),*) '  Target key, rank: ', this%key, this%data_rank
          write(fates_log(),*) '  Source key, rank: ', var%key, var%data_rank
        end if

        call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end subroutine CompareRegistryVariableSizes
    
  ! ====================================================================================

end module FatesInterfaceVariableTypeMod