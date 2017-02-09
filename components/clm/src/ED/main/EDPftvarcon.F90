module EDPftvarcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use clm_varpar  , only : mxpft, numrad, ivis, inir, nvariants
  use shr_kind_mod, only : r8 => shr_kind_r8

  use FatesGlobals, only : fates_log
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  integer, parameter, public :: lower_bound_pft = 0
  integer, parameter, public :: lower_bound_general = 1
  
  !ED specific variables. 
  type, public ::  EDPftvarcon_type
     real(r8), allocatable :: max_dbh            (:) ! maximum dbh at which height growth ceases...
     real(r8), allocatable :: freezetol          (:) ! minimum temperature tolerance...
     real(r8), allocatable :: wood_density       (:) ! wood density  g cm^-3  ...
     real(r8), allocatable :: alpha_stem         (:) ! live stem turnover rate. y-1
     real(r8), allocatable :: hgt_min            (:) ! sapling height m
     real(r8), allocatable :: cushion            (:) ! labile carbon storage target as multiple of leaf pool.
     real(r8), allocatable :: leaf_stor_priority (:) ! leaf turnover vs labile carbon use prioritisation. (1 = lose  leaves, 0 = use store).
     real(r8), allocatable :: leafwatermax       (:) ! degree to which respiration is limited by btran if btran = 0
     real(r8), allocatable :: rootresist         (:)
     real(r8), allocatable :: soilbeta           (:)
     real(r8), allocatable :: crown              (:)
     real(r8), allocatable :: bark_scaler        (:)
     real(r8), allocatable :: crown_kill         (:)
     real(r8), allocatable :: initd              (:)
     real(r8), allocatable :: sd_mort            (:)
     real(r8), allocatable :: seed_rain          (:)
     real(r8), allocatable :: BB_slope           (:)
     real(r8), allocatable :: root_long          (:) ! root longevity (yrs)
     real(r8), allocatable :: clone_alloc        (:) ! fraction of carbon balance allocated to clonal reproduction.
     real(r8), allocatable :: seed_alloc         (:) ! fraction of carbon balance allocated to seeds.
     real(r8), allocatable :: sapwood_ratio      (:) ! amount of sapwood per unit leaf carbon and m of height. gC/gC/m
     real(r8), allocatable :: dbh2h_m            (:) ! allocation parameter m from dbh to height
     real(r8), allocatable :: woody(:)
     real(r8), allocatable :: stress_decid(:)
     real(r8), allocatable :: season_decid(:)
     real(r8), allocatable :: evergreen(:)
     real(r8), allocatable :: froot_leaf(:)
     real(r8), allocatable :: slatop(:)
     real(r8), allocatable :: leaf_long(:)
     real(r8), allocatable :: roota_par(:)
     real(r8), allocatable :: rootb_par(:)
     real(r8), allocatable :: lf_flab(:)
     real(r8), allocatable :: lf_fcel(:)
     real(r8), allocatable :: lf_flig(:)
     real(r8), allocatable :: fr_flab(:)
     real(r8), allocatable :: fr_fcel(:)
     real(r8), allocatable :: fr_flig(:)
     real(r8), allocatable :: xl(:)
     real(r8), allocatable :: c3psn(:)
     real(r8), allocatable :: flnr(:)
     real(r8), allocatable :: fnitr(:)
     real(r8), allocatable :: leafcn(:)
     real(r8), allocatable :: frootcn(:)
     real(r8), allocatable :: smpso(:)
     real(r8), allocatable :: smpsc(:)
     real(r8), allocatable :: grperc(:) ! NOTE(bja, 2017-01) moved from EDParamsMod, was allocated as (maxPft=79), not (0:mxpft=78)!
     real(r8), allocatable :: rhol(:, :)
     real(r8), allocatable :: rhos(:, :)
     real(r8), allocatable :: taul(:, :)
     real(r8), allocatable :: taus(:, :)
     real(r8), allocatable :: rootprof_beta(:, :)
   contains
     procedure, public :: Init => EDpftconInit
     procedure, public :: Register
     procedure, public :: Receive
     procedure, private :: Register_PFT
     procedure, private :: Receive_PFT
     procedure, private :: Register_PFT_nvariants
     procedure, private :: Receive_PFT_nvariants
     procedure, private :: Register_PFT_numrad
     procedure, private :: Receive_PFT_numrad
  end type EDPftvarcon_type

  type(EDPftvarcon_type), public :: EDPftvarcon_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EDpftconrd ! Read and initialize vegetation (PFT) constants

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EDpftconInit(this)

    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this

  end subroutine EDpftconInit

  !-----------------------------------------------------------------------
  subroutine Register(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%Register_PFT(fates_params)
    call this%Register_PFT_numrad(fates_params)
    call this%Register_PFT_nvariants(fates_params)

  end subroutine Register

  !-----------------------------------------------------------------------
  subroutine Receive(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%Receive_PFT(fates_params)
    call this%Receive_PFT_numrad(fates_params)
    call this%Receive_PFT_nvariants(fates_params)

  end subroutine Receive

  !-----------------------------------------------------------------------
  subroutine Register_PFT(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)

    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
    !X!         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'max_dbh'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'freezetol'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'wood_density'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'alpha_stem'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'hgt_min'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'cushion'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'leaf_stor_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'leafwatermax'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'rootresist'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'soilbeta'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'crown'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'bark_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'crown_kill'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'initd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'sd_mort'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'seed_rain'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'BB_slope'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'root_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'clone_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'seed_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'sapwood_ratio'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'woody'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'stress_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'season_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'evergreen'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'froot_leaf'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'slatop'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'leaf_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'roota_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'rootb_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'lf_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'lf_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'lf_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fr_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fr_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fr_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'xl'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'c3psn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'flnr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fnitr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'leafcn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'frootcn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'smpso'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'smpsc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'grperc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)


  end subroutine Register_PFT

  !-----------------------------------------------------------------------
  subroutine Receive_PFT(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    name = 'max_dbh'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%max_dbh)

    name = 'freezetol'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%freezetol)

    name = 'wood_density'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%wood_density)

    name = 'alpha_stem'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%alpha_stem)

    name = 'hgt_min'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hgt_min)

    name = 'cushion'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%cushion)

    name = 'leaf_stor_priority'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leaf_stor_priority)

    name = 'leafwatermax'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leafwatermax)

    name = 'rootresist'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootresist)

    name = 'soilbeta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%soilbeta)

    name = 'crown'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown)

    name = 'bark_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%bark_scaler)

    name = 'crown_kill'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown_kill)

    name = 'initd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%initd)

    name = 'sd_mort'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%sd_mort)

    name = 'seed_rain'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_rain)

    name = 'BB_slope'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%BB_slope)

    name = 'root_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%root_long)

    name = 'clone_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%clone_alloc)

    name = 'seed_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_alloc)

    name = 'sapwood_ratio'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%sapwood_ratio)

    name = 'woody'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%woody)

    name = 'stress_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%stress_decid)

    name = 'season_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%season_decid)

    name = 'evergreen'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%evergreen)

    name = 'froot_leaf'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%froot_leaf)

    name = 'slatop'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%slatop)

    name = 'leaf_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leaf_long)

    name = 'roota_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%roota_par)

    name = 'rootb_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootb_par)

    name = 'lf_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flab)

    name = 'lf_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_fcel)

    name = 'lf_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flig)

    name = 'fr_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flab)

    name = 'fr_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_fcel)

    name = 'fr_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flig)

    name = 'xl'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%xl)

    name = 'c3psn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%c3psn)

    name = 'flnr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%flnr)

    name = 'fnitr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fnitr)

    name = 'leafcn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leafcn)

    name = 'frootcn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%frootcn)

    name = 'smpso'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpso)

    name = 'smpsc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpsc)

    name = 'grperc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%grperc)

  end subroutine Receive_PFT

  !-----------------------------------------------------------------------
  subroutine Register_PFT_numrad(this, fates_params)
    ! NOTE(bja, 2017-02) these are 2-d parameters, but they are
    ! currently stored in the parameter file as separate 1-d
    ! arrays. We have to register the parameters as 1-d arrays as they
    ! are on the parameter file. We store them as 2-d in the receive step.
    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)
    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)
    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
    !X!         dimension_names=dim_names)

    name = 'rholvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'rholnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'rhosvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'rhosnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'taulvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'taulnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'tausvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'tausnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)


  end subroutine Register_PFT_numrad

  !-----------------------------------------------------------------------
  subroutine Receive_PFT_numrad(this, fates_params)
    ! NOTE(bja, 2017-02) these are 2-d parameters, but they are
    ! currently stored in the parameter file as separate 1-d arrays.
    ! We can't allocate slices of arrays separately, so we have to
    ! manually allocate the memory here, retreive into a dummy array,
    ! and copy. All parameters in this subroutine are sized the same,
    ! so we can reused the dummy array. If someone wants to cleanup
    ! the input file, all this complexity can be removed.
    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length, max_dimensions

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    integer :: index
    integer :: dimension_shape
    integer :: dimension_sizes(max_dimensions)
    character(len=param_string_length) :: dimension_names(max_dimensions)
    logical :: is_host_param

    integer :: lower_bound_1, upper_bound_1, lower_bound_2, upper_bound_2
    real(r8), allocatable :: dummy_data(:)

    ! Fetch metadata from a representative variable. All variables
    ! called by this subroutine must be dimensioned the same way!
    name = 'rholvis'
    index = fates_params%FindIndex(name)
    call fates_params%GetMetaData(index, name, dimension_shape, dimension_sizes, dimension_names, is_host_param)
    lower_bound_1 = lower_bound_pft
    upper_bound_1 = lower_bound_pft + dimension_sizes(1) - 1
    lower_bound_2 = lower_bound_general
    upper_bound_2 = numrad

    allocate(dummy_data(lower_bound_1:upper_bound_1))

    !
    ! received rhol data
    !
    allocate(this%rhol(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'rholvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'rholnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received rhos data
    !
    allocate(this%rhos(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'rhosvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'rhosnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taul data
    !
    allocate(this%taul(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'taulvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'taulnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taus data
    !
    allocate(this%taus(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'tausvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taus(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'tausnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taus(lower_bound_1:upper_bound_1, inir) = dummy_data

  end subroutine Receive_PFT_numrad

  !-----------------------------------------------------------------------
  subroutine Register_PFT_nvariants(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_variants, dimension_name_pft, dimension_shape_2d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    ! NOTE(bja, 2017-01) initialization doesn't seem to work correctly
    ! if dim_names has a parameter qualifier.
    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_variants

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
    !X!         dimension_names=dim_names)

    name = 'rootprof_beta'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

  end subroutine Register_PFT_nvariants

  !-----------------------------------------------------------------------
  subroutine Receive_PFT_nvariants(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    name = 'rootprof_beta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootprof_beta)

  end subroutine Receive_PFT_nvariants

  !-----------------------------------------------------------------------
  subroutine EDpftconrd( ncid )
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t, ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    implicit none
    !
    type(file_desc_t), intent(inout) :: ncid   ! pio netCDF file id

    ! !LOCAL VARIABLES:

    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'EDpftconrd'              ! subroutine name

    !X!    call ncd_io('max_dbh',EDPftvarcon_inst%max_dbh, 'read', ncid, readvar=readv)
    !X!    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    !X!     call ncd_io('freezetol',EDPftvarcon_inst%freezetol, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('wood_density',EDPftvarcon_inst%wood_density, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('alpha_stem',EDPftvarcon_inst%alpha_stem, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('hgt_min',EDPftvarcon_inst%hgt_min, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('cushion',EDPftvarcon_inst%cushion, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('leaf_stor_priority',EDPftvarcon_inst%leaf_stor_priority, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('leafwatermax',EDPftvarcon_inst%leafwatermax, 'read', ncid, readvar=readv)
    !X!     if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('rootresist',EDPftvarcon_inst%rootresist,'read',     ncid, readvar=readv)
    !X!     if  ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )
    !X!
    !X!     call ncd_io('soilbeta',EDPftvarcon_inst%soilbeta,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('crown',EDPftvarcon_inst%crown,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('bark_scaler',EDPftvarcon_inst%bark_scaler,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('crown_kill',EDPftvarcon_inst%crown_kill,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('initd',EDPftvarcon_inst%initd,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('sd_mort',EDPftvarcon_inst%sd_mort,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('seed_rain',EDPftvarcon_inst%seed_rain,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('BB_slope',EDPftvarcon_inst%BB_slope,'read',         ncid, readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('root_long',EDPftvarcon_inst%root_long, 'read', ncid,  readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('seed_alloc',EDPftvarcon_inst%seed_alloc, 'read', ncid,  readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('clone_alloc',EDPftvarcon_inst%clone_alloc, 'read', ncid,  readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('sapwood_ratio',EDPftvarcon_inst%sapwood_ratio, 'read', ncid,  readvar=readv)
    !X!     if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')
    !X!
    !X!     call ncd_io('woody', EDPftvarcon_inst%woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('stress_decid', EDPftvarcon_inst%stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('season_decid', EDPftvarcon_inst%season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('evergreen', EDPftvarcon_inst%evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('froot_leaf', EDPftvarcon_inst%froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('slatop', EDPftvarcon_inst%slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('leaf_long', EDPftvarcon_inst%leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    !X!    call ncd_io('rootprof_beta', EDPftvarcon_inst%rootprof_beta, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    !X!     call ncd_io('roota_par', EDPftvarcon_inst%roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('rootb_par', EDPftvarcon_inst%rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('lf_flab', EDPftvarcon_inst%lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('lf_fcel', EDPftvarcon_inst%lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('lf_flig', EDPftvarcon_inst%lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('fr_flab', EDPftvarcon_inst%fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('fr_fcel', EDPftvarcon_inst%fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('fr_flig', EDPftvarcon_inst%fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    !X!    call ncd_io('rholvis', EDPftvarcon_inst%rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('rholnir', EDPftvarcon_inst%rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('rhosvis', EDPftvarcon_inst%rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('rhosnir', EDPftvarcon_inst% rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('taulvis', EDPftvarcon_inst%taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('taulnir', EDPftvarcon_inst%taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('tausvis', EDPftvarcon_inst%taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!    call ncd_io('tausnir', EDPftvarcon_inst%taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    !X!     call ncd_io('xl', EDPftvarcon_inst%xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('c3psn', EDPftvarcon_inst%c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('flnr', EDPftvarcon_inst%flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('fnitr', EDPftvarcon_inst%fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('leafcn', EDPftvarcon_inst%leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('frootcn', EDPftvarcon_inst%frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('smpso', EDPftvarcon_inst%smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('smpsc', EDPftvarcon_inst%smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    !X!
    !X!     call ncd_io('grperc', EDPftvarcon_inst%grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !X!     if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    ! HOLDING ON SEW ENSITIVITY-ANALYSIS PARAMETERS UNTIL MACHINE CONFIGS SET RGK/CX
    !    call ncd_io('dbh2h_m',EDPftvarcon_inst%dbh2h_m, 'read', ncid,  readvar=readv)
    !    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

  end subroutine EDpftconrd

end module EDPftvarcon

