module FatesInterfaceVariableTypeMod

  ! This module contains the type definition and associated type-bound procedures
  ! used to create an indexed list of associated HLM and FATES variables that are
  ! related across the application programming interface.
  ! This method is largely inspired by the FATES history infrastructure

  use shr_log_mod           , only : errMsg => shr_log_errMsg
 
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_long_string_length
  use FatesConstantsMod, only : fates_unset_int

  implicit none
  private

  ! Interface registry variable type
  type, public :: fates_interface_variable_type
    
    character(len=48) :: key            ! common registry key
    class(*), pointer :: data0d         ! scalar polymorphic data pointer
    class(*), pointer :: data1d(:)      ! 1D polymorphic data pointer
    class(*), pointer :: data2d(:,:)    ! 2D polymorphic data pointer
    class(*), pointer :: data3d(:,:,:)  ! 3D polymorphic data pointer
    logical           :: active         ! true if the variable is used by the host land model
    logical           :: accumulate     ! If true, this variable should add the source to the target
    integer           :: bc_dir         ! 0 if bc_in, 1 if bc_out
    integer           :: data_rank      ! rank of the variable (0, 1, 2, or 3)
    integer           :: update_frequency ! frequency of updates 
    integer, allocatable :: data_size(:)   ! size of the first dimension of the variable

    contains
      procedure :: CheckBounds
      procedure :: Initialize => InitializeInterfaceVariable
      procedure :: Update     => UpdateInterfaceVariable
      procedure :: Dump

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
  
  subroutine Dump(this)

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

      ! Initialize registry variable components that are updated at variable definition
      this%key = key 
      this%update_frequency = update_frequency
      this%bc_dir = bc_dir

    end subroutine InitializeInterfaceVariable

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_0d(this, data, active, accumulate)

      class(fates_interface_variable_type), intent(inout) :: this
      class(*), target, intent(in) :: data
      logical, intent(in)          :: active
      logical, intent(in)          :: accumulate

      this%data0d => data
      this%active = active
      this%accumulate = accumulate
      this%data_rank = rank(data)

    end subroutine RegisterInterfaceVariable_0d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_1d(this, data, active, accumulate)

      class(fates_interface_variable_type), intent(inout) :: this

      class(*), target, intent(in) :: data(:)
      logical, intent(in)          :: active
      logical, intent(in)          :: accumulate

      this%data1d => data(:)
      this%active = active
      this%accumulate = accumulate
      this%data_rank = rank(data)
      this%data_size(1) = size(data, dim=1)

    end subroutine RegisterInterfaceVariable_1d

  ! ====================================================================================
    
    subroutine RegisterInterfaceVariable_2d(this, data, active, accumulate)

      class(fates_interface_variable_type), intent(inout) :: this
      class(*), target, intent(in)  :: data(:,:)
      logical, intent(in)           :: active
      logical, intent(in)          :: accumulate

      this%data2d => data(:,:)
      this%active = active
      this%accumulate = accumulate
      this%data_rank = rank(data) 
      this%data_size(1) = size(data, dim=1)
      this%data_size(2) = size(data, dim=2)

    end subroutine RegisterInterfaceVariable_2d

  ! ====================================================================================
    
    subroutine UpdateInterfaceVariable(this, var, scalar)

      ! Arguments
      class(fates_interface_variable_type), intent(inout) :: this ! variable being updated
      class(fates_interface_variable_type), intent(in)    :: var  ! variable update source
      real(r8), intent(in), optional                      :: scalar ! value to scale variable update 

      ! Locals
      real(r8) :: scalar_local
      character(len=fates_long_string_length) :: msg_mismatch = 'FATES ERROR: Mismatched interface variable types'

      ! Check if scalar is present and set default value to one
      ! Currently this assumes that the only real values are to be scaled
      if (present(scalar)) then
        scalar_local = scalar
      else
        scalar_local = 1.0_r8
      end if
      
      ! Check that the dimensions of the source and target match
      call this%CompareRegistryVariableSizes(var)
          
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
                    dest = dest + source * scalar_local
                  else
                    dest = source * scalar_local
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data0d)
                type is (integer)
                  if (this%accumulate) then
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
                    dest = dest + source * scalar_local
                  else
                    dest = source * scalar_local
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data1d)
                type is (integer)
                  if (this%accumulate) then
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
                    dest = dest + source * scalar_local
                  else
                    dest = source * scalar_local
                  end if
                class default
                  write(fates_log(),*), msg_mismatch 
                  call endrun(msg=errMsg(__FILE__, __LINE__))
              end select
            type is (integer)
              select type(source => var%data2d)
                type is (integer)
                  if (this%accumulate) then
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