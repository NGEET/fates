module FatesSynchronizedParamsMod

  ! NOTE: We currently do NOT use any "shared" or syncronized parameters
  !       between FATES and its hosts.  We previously shared q10 values.
  !       I will leave these values commented out instead of deleted
  !       to serve as a template for the possibility of future parameters.
  !       RGK 05-2019


  !-----------------------------------------------------------------------
  !
  ! !USES:
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  private ! Modules are private by default

  ! FatesSynchronizedParamsInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public  :: FatesSynchronizedParamsType
!      real(r8) :: Q10      ! temperature dependence
!      real(r8) :: froz_q10 ! separate q10 for frozen soil respiration rates
    contains

      ! Public member functions
      procedure :: RegisterParams
      procedure :: ReceiveParams

      ! Private member functions
      procedure, private :: Init
      procedure, private :: RegisterParamsScalar
      procedure, private :: ReceiveParamsScalar

  end type FatesSynchronizedParamsType

  type(FatesSynchronizedParamsType), public :: FatesSynchronizedParamsInst
  
  character(len=*), parameter :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------
  
contains

  subroutine Init(this)
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    class(FatesSynchronizedParamsType), intent(inout) :: this

!    this%Q10 = nan
!    this%froz_q10 = nan
    
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine RegisterParams(this, fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(FatesSynchronizedParamsType), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%Init()
    call this%RegisterParamsScalar(fates_params)
    
  end subroutine RegisterParams

  !-----------------------------------------------------------------------
  subroutine ReceiveParams(this, fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(FatesSynchronizedParamsType), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%ReceiveParamsScalar(fates_params)
    
  end subroutine ReceiveParams

  !-----------------------------------------------------------------------
  subroutine RegisterParamsScalar(this, fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_host_allpfts, dimension_shape_1d

    implicit none

    class(FatesSynchronizedParamsType), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_host_allpfts/)
    character(len=param_string_length) :: name

    call this%Init()


  end subroutine RegisterParamsScalar

  !-----------------------------------------------------------------------
  subroutine ReceiveParamsScalar(this, fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(FatesSynchronizedParamsType), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

!    name = 'q10_mr'
!    call fates_params%RetreiveParameter(name=name, &
!         data=this%Q10)

!    name = 'froz_q10'
!    call fates_params%RetreiveParameter(name=name, &
!         data=this%froz_q10)

  end subroutine ReceiveParamsScalar

end module FatesSynchronizedParamsMod
