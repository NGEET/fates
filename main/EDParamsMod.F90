module EDParamsMod
   !
   ! module that deals with reading the ED parameter file
   !
   use shr_kind_mod      , only: r8 => shr_kind_r8
   use EDtypesMod        , only: maxPft
   use FatesParametersInterface, only : param_string_length
   
   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be 
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   
   real(r8),protected :: ED_val_grass_spread
   real(r8),protected :: ED_val_comp_excln
   real(r8),protected :: ED_val_stress_mort
   real(r8),protected :: ED_val_dispersal
   real(r8),protected :: ED_val_maxspread
   real(r8),protected :: ED_val_minspread
   real(r8),protected :: ED_val_init_litter
   real(r8),protected :: ED_val_nignitions
   real(r8),protected :: ED_val_understorey_death
   real(r8),protected :: ED_val_ag_biomass
   real(r8),protected :: ED_val_cwd_fcel
   real(r8),protected :: ED_val_cwd_flig
   real(r8),protected :: ED_val_bbopt_c3
   real(r8),protected :: ED_val_bbopt_c4
   real(r8),protected :: ED_val_base_mr_20
   real(r8),protected :: ED_val_phen_drought_threshold
   real(r8),protected :: ED_val_phen_doff_time
   real(r8),protected :: ED_val_phen_a
   real(r8),protected :: ED_val_phen_b
   real(r8),protected :: ED_val_phen_c
   real(r8),protected :: ED_val_phen_chiltemp
   real(r8),protected :: ED_val_phen_mindayson
   real(r8),protected :: ED_val_phen_ncolddayslim
   real(r8),protected :: ED_val_phen_coldtemp
   real(r8),protected :: ED_val_cohort_fusion_tol
   real(r8),protected :: ED_val_patch_fusion_tol
  
   character(len=param_string_length),parameter :: ED_name_grass_spread = "fates_grass_spread"
   character(len=param_string_length),parameter :: ED_name_comp_excln = "fates_comp_excln"
   character(len=param_string_length),parameter :: ED_name_stress_mort = "fates_stress_mort"
   character(len=param_string_length),parameter :: ED_name_dispersal = "fates_dispersal"
   character(len=param_string_length),parameter :: ED_name_maxspread = "fates_maxspread"
   character(len=param_string_length),parameter :: ED_name_minspread = "fates_minspread"
   character(len=param_string_length),parameter :: ED_name_init_litter = "fates_init_litter"
   character(len=param_string_length),parameter :: ED_name_nignitions = "fates_nfires"
   character(len=param_string_length),parameter :: ED_name_understorey_death = "fates_understorey_death"
   character(len=param_string_length),parameter :: ED_name_ag_biomass= "fates_ag_biomass"   
   character(len=param_string_length),parameter :: ED_name_cwd_fcel= "fates_cwd_fcel"   
   character(len=param_string_length),parameter :: ED_name_cwd_flig= "fates_cwd_flig"   
   character(len=param_string_length),parameter :: ED_name_bbopt_c3= "fates_bbopt_c3"   
   character(len=param_string_length),parameter :: ED_name_bbopt_c4= "fates_bbopt_c4"   
   character(len=param_string_length),parameter :: ED_name_base_mr_20= "fates_base_mr_20"   
   character(len=param_string_length),parameter :: ED_name_phen_drought_threshold= "fates_phen_drought_threshold"   
   character(len=param_string_length),parameter :: ED_name_phen_doff_time= "fates_phen_doff_time"   
   character(len=param_string_length),parameter :: ED_name_phen_a= "fates_phen_a"   
   character(len=param_string_length),parameter :: ED_name_phen_b= "fates_phen_b"   
   character(len=param_string_length),parameter :: ED_name_phen_c= "fates_phen_c"   
   character(len=param_string_length),parameter :: ED_name_phen_chiltemp= "fates_phen_chiltemp"   
   character(len=param_string_length),parameter :: ED_name_phen_mindayson= "fates_phen_mindayson"   
   character(len=param_string_length),parameter :: ED_name_phen_ncolddayslim= "fates_phen_ncolddayslim"   
   character(len=param_string_length),parameter :: ED_name_phen_coldtemp= "fates_phen_coldtemp"   
   character(len=param_string_length),parameter :: ED_name_cohort_fusion_tol= "fates_cohort_fusion_tol"   
   character(len=param_string_length),parameter :: ED_name_patch_fusion_tol= "fates_patch_fusion_tol"   
   
   public :: FatesParamsInit
   public :: FatesRegisterParams
   public :: FatesReceiveParams

  real(r8), protected :: fates_mortality_disturbance_fraction  = 0.5_r8 ! the fraction of canopy mortality that results in disturbance (i.e. transfer of area from new to old patch)                                                              
  
contains

  !-----------------------------------------------------------------------
  subroutine FatesParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    ED_val_grass_spread = nan
    ED_val_comp_excln = nan
    ED_val_stress_mort = nan
    ED_val_dispersal = nan
    ED_val_maxspread = nan
    ED_val_minspread = nan
    ED_val_init_litter = nan
    ED_val_nignitions = nan
    ED_val_understorey_death = nan
    ED_val_ag_biomass = nan
    ED_val_cwd_fcel = nan
    ED_val_cwd_flig = nan
    ED_val_bbopt_c3 = nan
    ED_val_bbopt_c4 = nan
    ED_val_base_mr_20 = nan
    ED_val_phen_drought_threshold = nan
    ED_val_phen_doff_time = nan
    ED_val_phen_a = nan
    ED_val_phen_b = nan
    ED_val_phen_c = nan
    ED_val_phen_chiltemp = nan
    ED_val_phen_mindayson = nan
    ED_val_phen_ncolddayslim = nan
    ED_val_phen_coldtemp = nan
    ED_val_cohort_fusion_tol = nan
    ED_val_patch_fusion_tol = nan

  end subroutine FatesParamsInit

  !-----------------------------------------------------------------------
  subroutine FatesRegisterParams(fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar1d, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_scalar1d/)

    call FatesParamsInit()

    call fates_params%RegisterParameter(name=ED_name_grass_spread, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_comp_excln, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_grass_spread, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_comp_excln, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_stress_mort, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_dispersal, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_maxspread, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_minspread, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_init_litter, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_nignitions, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_understorey_death, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_ag_biomass, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cwd_fcel, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cwd_flig, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_bbopt_c3, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_bbopt_c4, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_base_mr_20, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_drought_threshold, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_doff_time, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_a, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_b, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_c, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_chiltemp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_mindayson, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_ncolddayslim, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_coldtemp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cohort_fusion_tol, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_patch_fusion_tol, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

  end subroutine FatesRegisterParams

  
  !-----------------------------------------------------------------------
  subroutine FatesReceiveParams(fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetreiveParameter(name=ED_name_grass_spread, &
         data=ED_val_grass_spread)

    call fates_params%RetreiveParameter(name=ED_name_comp_excln, &
         data=ED_val_comp_excln)

    call fates_params%RetreiveParameter(name=ED_name_grass_spread, &
         data=ED_val_grass_spread)

    call fates_params%RetreiveParameter(name=ED_name_comp_excln, &
         data=ED_val_comp_excln)

    call fates_params%RetreiveParameter(name=ED_name_stress_mort, &
         data=ED_val_stress_mort)

    call fates_params%RetreiveParameter(name=ED_name_dispersal, &
         data=ED_val_dispersal)

    call fates_params%RetreiveParameter(name=ED_name_maxspread, &
         data=ED_val_maxspread)

    call fates_params%RetreiveParameter(name=ED_name_minspread, &
         data=ED_val_minspread)

    call fates_params%RetreiveParameter(name=ED_name_init_litter, &
         data=ED_val_init_litter)

    call fates_params%RetreiveParameter(name=ED_name_nignitions, &
         data=ED_val_nignitions)

    call fates_params%RetreiveParameter(name=ED_name_understorey_death, &
         data=ED_val_understorey_death)

    call fates_params%RetreiveParameter(name=ED_name_ag_biomass, &
         data=ED_val_ag_biomass)

    call fates_params%RetreiveParameter(name=ED_name_cwd_fcel, &
         data=ED_val_cwd_fcel)

    call fates_params%RetreiveParameter(name=ED_name_cwd_flig, &
         data=ED_val_cwd_flig)

    call fates_params%RetreiveParameter(name=ED_name_bbopt_c3, &
         data=ED_val_bbopt_c3)

    call fates_params%RetreiveParameter(name=ED_name_bbopt_c4, &
         data=ED_val_bbopt_c4)

    call fates_params%RetreiveParameter(name=ED_name_base_mr_20, &
         data=ED_val_base_mr_20)

    call fates_params%RetreiveParameter(name=ED_name_phen_drought_threshold, &
         data=ED_val_phen_drought_threshold)

    call fates_params%RetreiveParameter(name=ED_name_phen_doff_time, &
         data=ED_val_phen_doff_time)

    call fates_params%RetreiveParameter(name=ED_name_phen_a, &
         data=ED_val_phen_a)

    call fates_params%RetreiveParameter(name=ED_name_phen_b, &
         data=ED_val_phen_b)

    call fates_params%RetreiveParameter(name=ED_name_phen_c, &
         data=ED_val_phen_c)

    call fates_params%RetreiveParameter(name=ED_name_phen_chiltemp, &
         data=ED_val_phen_chiltemp)

    call fates_params%RetreiveParameter(name=ED_name_phen_mindayson, &
         data=ED_val_phen_mindayson)

    call fates_params%RetreiveParameter(name=ED_name_phen_ncolddayslim, &
         data=ED_val_phen_ncolddayslim)

    call fates_params%RetreiveParameter(name=ED_name_phen_coldtemp, &
         data=ED_val_phen_coldtemp)

    call fates_params%RetreiveParameter(name=ED_name_cohort_fusion_tol, &
         data=ED_val_cohort_fusion_tol)

    call fates_params%RetreiveParameter(name=ED_name_patch_fusion_tol, &
         data=ED_val_patch_fusion_tol)

  end subroutine FatesReceiveParams
  
end module EDParamsMod
