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
     real(r8), allocatable :: dbh2h_m(:)
     real(r8), allocatable :: dbh2h_c(:)
     real(r8), allocatable :: dbh2bl_a(:)
     real(r8), allocatable :: dbh2bl_b(:)
     real(r8), allocatable :: dbh2bl_dbh2carea_expnt_diff(:)
     real(r8), allocatable :: dbh2bl_c(:)
     real(r8), allocatable :: dbh2bl_slascaler(:)
     real(r8), allocatable :: sai_scaler(:)
     real(r8), allocatable :: dbh2bd_a(:)
     real(r8), allocatable :: dbh2bd_b(:)
     real(r8), allocatable :: dbh2bd_c(:)
     real(r8), allocatable :: dbh2bd_d(:)
     real(r8), allocatable :: bmort(:)
     real(r8), allocatable :: hf_sm_threshold(:)
     real(r8), allocatable :: vcmaxha(:)
     real(r8), allocatable :: jmaxha(:)
     real(r8), allocatable :: tpuha(:)
     real(r8), allocatable :: vcmaxhd(:)
     real(r8), allocatable :: jmaxhd(:)
     real(r8), allocatable :: tpuhd(:)
     real(r8), allocatable :: vcmaxse(:)
     real(r8), allocatable :: jmaxse(:)
     real(r8), allocatable :: tpuse(:)
     real(r8), allocatable :: germination_timescale(:)
     real(r8), allocatable :: seed_decay_turnover(:)
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

    name = 'fates_max_dbh'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_freezetol'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_wood_density'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alpha_stem'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hgt_min'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cushion'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_stor_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leafwatermax'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_rootresist'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_soilbeta'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_crown'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_bark_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_crown_kill'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_initd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_sd_mort'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_rain'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_BB_slope'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_root_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_clone_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_sapwood_ratio'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_woody'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_stress_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_season_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_evergreen'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_froot_leaf'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_slatop'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_roota_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_rootb_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_xl'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_c3psn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_flnr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fnitr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leafcn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_frootcn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_smpso'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_smpsc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_grperc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2h_m'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2h_c'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bl_a'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bl_b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bl_dbh2carea_expnt_diff'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bl_c'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bl_slascaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_sai_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bd_a'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bd_b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bd_c'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_dbh2bd_d'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_bmort'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hf_sm_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_vcmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_jmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_tpuha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_vcmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_jmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_tpuhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_vcmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_jmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_tpuse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_germination_timescale'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_decay_turnover'
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

    name = 'fates_max_dbh'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%max_dbh)

    name = 'fates_freezetol'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%freezetol)

    name = 'fates_wood_density'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%wood_density)

    name = 'fates_alpha_stem'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%alpha_stem)

    name = 'fates_hgt_min'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hgt_min)

    name = 'fates_cushion'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%cushion)

    name = 'fates_leaf_stor_priority'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leaf_stor_priority)

    name = 'fates_leafwatermax'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leafwatermax)

    name = 'fates_rootresist'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootresist)

    name = 'fates_soilbeta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%soilbeta)

    name = 'fates_crown'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown)

    name = 'fates_bark_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%bark_scaler)

    name = 'fates_crown_kill'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown_kill)

    name = 'fates_initd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%initd)

    name = 'fates_sd_mort'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%sd_mort)

    name = 'fates_seed_rain'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_rain)

    name = 'fates_BB_slope'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%BB_slope)

    name = 'fates_root_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%root_long)

    name = 'fates_clone_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%clone_alloc)

    name = 'fates_seed_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_alloc)

    name = 'fates_sapwood_ratio'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%sapwood_ratio)

    name = 'fates_woody'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%woody)

    name = 'fates_stress_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%stress_decid)

    name = 'fates_season_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%season_decid)

    name = 'fates_evergreen'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%evergreen)

    name = 'fates_froot_leaf'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%froot_leaf)

    name = 'fates_slatop'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%slatop)

    name = 'fates_leaf_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leaf_long)

    name = 'fates_roota_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%roota_par)

    name = 'fates_rootb_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootb_par)

    name = 'fates_lf_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flab)

    name = 'fates_lf_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_fcel)

    name = 'fates_lf_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flig)

    name = 'fates_fr_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flab)

    name = 'fates_fr_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_fcel)

    name = 'fates_fr_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flig)

    name = 'fates_xl'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%xl)

    name = 'fates_c3psn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%c3psn)

    name = 'fates_flnr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%flnr)

    name = 'fates_fnitr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fnitr)

    name = 'fates_leafcn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leafcn)

    name = 'fates_frootcn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%frootcn)

    name = 'fates_smpso'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpso)

    name = 'fates_smpsc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpsc)

    name = 'fates_grperc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%grperc)

    name = 'fates_dbh2h_m'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2h_m)

    name = 'fates_dbh2h_c'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2h_c)

    name = 'fates_dbh2bl_a'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bl_a)

    name = 'fates_dbh2bl_b'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bl_b)

    name = 'fates_dbh2bl_dbh2carea_expnt_diff'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bl_dbh2carea_expnt_diff)

    name = 'fates_dbh2bl_c'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bl_c)

    name = 'fates_dbh2bl_slascaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bl_slascaler)

    name = 'fates_sai_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%sai_scaler)

    name = 'fates_dbh2bd_a'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bd_a)

    name = 'fates_dbh2bd_b'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bd_b)

    name = 'fates_dbh2bd_c'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bd_c)

    name = 'fates_dbh2bd_d'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh2bd_d)

    name = 'fates_bmort'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%bmort)

    name = 'fates_hf_sm_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hf_sm_threshold)

    name = 'fates_vcmaxha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxha)

    name = 'fates_jmaxha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxha)

    name = 'fates_tpuha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuha)

    name = 'fates_vcmaxhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxhd)

    name = 'fates_jmaxhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxhd)

    name = 'fates_tpuhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuhd)

    name = 'fates_vcmaxse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxse)

    name = 'fates_jmaxse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxse)

    name = 'fates_tpuse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuse)

    name = 'fates_germination_timescale'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%germination_timescale)

    name = 'fates_seed_decay_turnover'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_decay_turnover)

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

    name = 'fates_rholvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rholnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rhosvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rhosnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_taulvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_taulnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_tausvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_tausnir'
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
    name = 'fates_rholvis'
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
    
    name = 'fates_rholvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_rholnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received rhos data
    !
    allocate(this%rhos(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_rhosvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_rhosnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taul data
    !
    allocate(this%taul(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_taulvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_taulnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taus data
    !
    allocate(this%taus(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_tausvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taus(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_tausnir'
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

    name = 'fates_rootprof_beta'
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

    name = 'fates_rootprof_beta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootprof_beta)

  end subroutine Receive_PFT_nvariants

end module EDPftvarcon

