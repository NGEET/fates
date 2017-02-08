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
   real(r8),protected :: ED_val_nfires
   real(r8),protected :: ED_val_understorey_death
   real(r8),protected :: ED_val_profile_tol
   real(r8),protected :: ED_val_ag_biomass
  
   character(len=param_string_length),parameter :: ED_name_grass_spread = "grass_spread"
   character(len=param_string_length),parameter :: ED_name_comp_excln = "comp_excln"
   character(len=param_string_length),parameter :: ED_name_stress_mort = "stress_mort"
   character(len=param_string_length),parameter :: ED_name_dispersal = "dispersal"
   character(len=param_string_length),parameter :: ED_name_maxspread = "maxspread"
   character(len=param_string_length),parameter :: ED_name_minspread = "minspread"
   character(len=param_string_length),parameter :: ED_name_init_litter = "init_litter"
   character(len=param_string_length),parameter :: ED_name_nfires = "nfires"
   character(len=param_string_length),parameter :: ED_name_understorey_death = "understorey_death"
   character(len=param_string_length),parameter :: ED_name_profile_tol = "profile_tol"
   character(len=param_string_length),parameter :: ED_name_ag_biomass= "ag_biomass"   
   
   public :: FatesParamsInit
   public :: FatesRegisterParams
   public :: FatesReceiveParams
  
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
    ED_val_nfires = nan
    ED_val_understorey_death = nan
    ED_val_profile_tol = nan
    ED_val_ag_biomass = nan

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

    call fates_params%RegisterParameter(name=ED_name_nfires, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_understorey_death, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_profile_tol, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_ag_biomass, dimension_shape=dimension_shape_1d, &
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

    call fates_params%RetreiveParameter(name=ED_name_nfires, &
         data=ED_val_nfires)

    call fates_params%RetreiveParameter(name=ED_name_understorey_death, &
         data=ED_val_understorey_death)

    call fates_params%RetreiveParameter(name=ED_name_profile_tol, &
         data=ED_val_profile_tol)

    call fates_params%RetreiveParameter(name=ED_name_ag_biomass, &
         data=ED_val_ag_biomass)

  end subroutine FatesReceiveParams
  
end module EDParamsMod
