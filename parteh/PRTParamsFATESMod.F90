module PRTInitParamsFatesMod

  ! This is a FATES specific module for loading parameters through
  ! the CLM/ELM module system.

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : itrue,ifalse
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : years_per_day
  use FatesInterfaceTypesMod, only : hlm_parteh_mode
  use PRTParametersMod,  only : prt_params
  use PRTGenericMod,     only : num_organ_types
  use PRTGenericMod,     only : leaf_organ, fnrt_organ, store_organ
  use PRTGenericMod,     only : sapw_organ, struct_organ, repro_organ
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log 
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use PRTGenericMod,     only : prt_cnp_flex_allom_hyp,prt_carbon_allom_hyp
  
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  integer, parameter, public :: lower_bound_pft = 1
  integer, parameter, public :: lower_bound_general = 1

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PRTRegisterParams
  public :: PRTReceiveParams
  public :: PRTCheckParams
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine PRTRegisterParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params
    
    call PRTRegisterPFT(fates_params)
    call PRTRegisterPFTOrgans(fates_params)
    call PRTRegisterPFTLeafAge(fates_params)
    call Register_PFT_nvariants(fates_params)

  end subroutine PRTRegisterParams

  !-----------------------------------------------------------------------
  subroutine PRTReceiveParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call PRTReceivePFT(fates_params)
    call PRTReceivePFTOrgans(fates_params)
    call PRTReceivePFTLeafAge(fates_params)
    call Receive_PFT_nvariants(fates_params)
    
  end subroutine PRTReceiveParams

  !-----------------------------------------------------------------------
  subroutine PRTRegisterPFT(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)
    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)
    character(len=param_string_length) :: name


    name = 'fates_phen_stress_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_season_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_evergreen'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    
    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
    !X!         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fnrt_prof_a'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fnrt_prof_b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fnrt_prof_mode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_woody'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_wood_density'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_slamax'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_slatop'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_sai_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_seed_dbh_repro_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_stor_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_senleaf_long_fdrought'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_root_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_alloc_mature'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_c2b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_l2fr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_grperc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_dbh_maxheight'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_hmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_lmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_fmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_amode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_stmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_cmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_smode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_la_per_sa_int'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_la_per_sa_slp'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb_frac'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_d2bl1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2bl2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2bl3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_blca_expnt_diff'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2ca_coefficient_max'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2ca_coefficient_min'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb4'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_retrans_mode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_branch_turnover'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
  end subroutine PRTRegisterPFT
  
  !-----------------------------------------------------------------------
  
  subroutine PRTReceivePFT(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=prt_params%)

    name = 'fates_phen_stress_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%stress_decid)

    name = 'fates_phen_season_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%season_decid)

    name = 'fates_phen_evergreen'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%evergreen)
    
    name = 'fates_leaf_slamax'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%slamax)
    
    name = 'fates_leaf_slatop'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%slatop)

    name = 'fates_allom_sai_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_sai_scaler)

    name = 'fates_fnrt_prof_a'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_a)

    name = 'fates_fnrt_prof_b'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_b)

    name = 'fates_fnrt_prof_mode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_mode)
    
    name = 'fates_woody'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%woody)
    
    name = 'fates_wood_density'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%wood_density)
    
    name = 'fates_seed_dbh_repro_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%dbh_repro_threshold)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%cushion)

    name = 'fates_leaf_stor_priority'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%leaf_stor_priority)

    name = 'fates_senleaf_long_fdrought'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=prt_params%senleaf_long_fdrought)

    name = 'fates_root_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%root_long)

    name = 'fates_seed_alloc_mature'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%seed_alloc_mature)

    name = 'fates_seed_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%seed_alloc)

    name = 'fates_c2b'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%c2b)

    name = 'fates_grperc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%grperc)

    name = 'fates_allom_dbh_maxheight'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=prt_params%allom_dbh_maxheight)

    name = 'fates_allom_hmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_hmode)

    name = 'fates_allom_lmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_lmode)

    name = 'fates_allom_fmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_fmode)

    name = 'fates_allom_amode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_amode)

    name = 'fates_allom_stmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_stmode)

    name = 'fates_allom_cmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_cmode)

    name = 'fates_allom_smode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_smode)

    name = 'fates_allom_la_per_sa_int'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_la_per_sa_int)

    name = 'fates_allom_la_per_sa_slp'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_la_per_sa_slp)

    name = 'fates_allom_l2fr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_l2fr)

    name = 'fates_allom_agb_frac'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_agb_frac)

    name = 'fates_allom_d2h1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2h1)

    name = 'fates_allom_d2h2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2h2)

    name = 'fates_allom_d2h3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2h3)

    name = 'fates_allom_d2bl1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl1)

    name = 'fates_allom_d2bl2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl2)

    name = 'fates_allom_d2bl3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl3)

    name = 'fates_allom_blca_expnt_diff'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_blca_expnt_diff)

    name = 'fates_allom_d2ca_coefficient_max'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2ca_coefficient_max)

    name = 'fates_allom_d2ca_coefficient_min'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_d2ca_coefficient_min)

    name = 'fates_allom_agb1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_agb1)

    name = 'fates_allom_agb2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_agb2)

    name = 'fates_allom_agb3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_agb3)

    name = 'fates_allom_agb4'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%allom_agb4)
	 
    name = 'fates_branch_turnover'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%branch_long)
    
    name = 'fates_turnover_retrans_mode'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=prt_params%turnover_retrans_mode)


  end subroutine PRTReceivePFT

  !-----------------------------------------------------------------------

  subroutine PRTRegisterPFTLeafAge(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_leaf_age
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_2d
    
    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_leaf_age

    name = 'fates_leaf_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    return
  end subroutine PRTRegisterPFTLeafAge

  ! =====================================================================================
  
  subroutine Register_PFT_nvariants(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_variants, dimension_name_pft, dimension_shape_2d

    implicit none

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

  end subroutine Register_PFT_nvariants

  ! =====================================================================================
  
  subroutine Receive_PFT_nvariants(fates_params)

    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

  end subroutine Receive_PFT_nvariants
  
  ! =====================================================================================

  subroutine PRTReceivePFTLeafAge(fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name

     name = 'fates_leaf_long'
     call fates_params%RetreiveParameterAllocate(name=name, &
          data=prt_params%leaf_long)

     return
  end subroutine PRTReceivePFTLeafAge

  ! =====================================================================================

  
  subroutine PRTRegisterPFTOrgans(fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_prt_organs
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_2d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    ! NOTE(bja, 2017-01) initialization doesn't seem to work correctly
    ! if dim_names has a parameter qualifier.
    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_prt_organs

    name = 'fates_prt_nitr_stoich_p1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_nitr_stoich_p2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_phos_stoich_p1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_phos_stoich_p2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_alloc_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_carb_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)    

    name = 'fates_turnover_nitr_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_turnover_phos_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    
  end subroutine PRTRegisterPFTOrgans

  ! =====================================================================================

  subroutine PRTReceivePFTOrgans(fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name

     name = 'fates_prt_nitr_stoich_p1'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%nitr_stoich_p1)

     name = 'fates_prt_nitr_stoich_p2'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%nitr_stoich_p2)
     
     name = 'fates_prt_phos_stoich_p1'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%phos_stoich_p1)

     name = 'fates_prt_phos_stoich_p2'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%phos_stoich_p2)
    
     name = 'fates_prt_alloc_priority'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%alloc_priority)

     name = 'fates_turnover_carb_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%turnover_carb_retrans)

     name = 'fates_turnover_nitr_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%turnover_nitr_retrans)

     name = 'fates_turnover_phos_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=prt_params%turnover_phos_retrans)

  end subroutine PRTReceivePFTOrgans

  ! ===============================================================================================
  
  subroutine FatesReportPFTParams(is_master)
     
     ! Argument
     logical, intent(in) :: is_master  ! Only log if this is the master proc

     logical, parameter :: debug_report = .false.
     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft,ipft
     
     npft = size(prt_params%allom_hmode,1)
     
     if(debug_report .and. is_master) then
        
        if(npft>100)then
           write(fates_log(),*) 'you are trying to report pft parameters during initialization'
           write(fates_log(),*) 'but you have so many that it is over-running the format spec'
           write(fates_log(),*) 'simply bump up the muptiplier in parameter fmt0 shown above'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        write(fates_log(),*) '-----------  FATES PARTEH Parameters -----------------'
        write(fates_log(),fmt0) 'stress_decid = ',prt_params%stress_decid
        write(fates_log(),fmt0) 'season_decid = ',prt_params%season_decid
        write(fates_log(),fmt0) 'evergreen = ',prt_params%evergreen
        write(fates_log(),fmt0) 'wood_density = ',prt_params%wood_density
        write(fates_log(),fmt0) 'dbh max height = ',prt_params%allom_dbh_maxheight
        write(fates_log(),fmt0) 'dbh mature = ',prt_params%dbh_repro_threshold
        write(fates_log(),fmt0) 'cushion = ',prt_params%cushion
        write(fates_log(),fmt0) 'leaf_stor_priority = ',prt_params%leaf_stor_priority
        write(fates_log(),fmt0) 'root_long = ',prt_params%root_long
        write(fates_log(),fmt0) 'senleaf_long_fdrought = ',prt_params%senleaf_long_fdrought
        write(fates_log(),fmt0) 'seed_alloc_mature = ',prt_params%seed_alloc_mature
        write(fates_log(),fmt0) 'seed_alloc = ',prt_params%seed_alloc
        write(fates_log(),fmt0) 'slamax = ',prt_params%slamax
        write(fates_log(),fmt0) 'slatop = ',prt_params%slatop
        write(fates_log(),fmt0) 'allom_sai_scaler = ',prt_params%allom_sai_scaler
        write(fates_log(),fmt0) 'leaf_long = ',prt_params%leaf_long
        write(fates_log(),fmt0) 'grperc = ',prt_params%grperc
        write(fates_log(),fmt0) 'c2b = ',prt_params%c2b
        write(fates_log(),fmt0) 'branch_turnover = ',prt_params%branch_long
        write(fates_log(),fmt0) 'allom_hmode = ',prt_params%allom_hmode
        write(fates_log(),fmt0) 'allom_lmode = ',prt_params%allom_lmode
        write(fates_log(),fmt0) 'allom_fmode = ',prt_params%allom_fmode
        write(fates_log(),fmt0) 'allom_amode = ',prt_params%allom_amode
        write(fates_log(),fmt0) 'allom_cmode = ',prt_params%allom_cmode
        write(fates_log(),fmt0) 'allom_smode = ',prt_params%allom_smode
        write(fates_log(),fmt0) 'allom_la_per_sa_int = ',prt_params%allom_la_per_sa_int
        write(fates_log(),fmt0) 'allom_la_per_sa_slp = ',prt_params%allom_la_per_sa_slp
        write(fates_log(),fmt0) 'allom_l2fr = ',prt_params%allom_l2fr
        write(fates_log(),fmt0) 'allom_agb_frac = ',prt_params%allom_agb_frac
        write(fates_log(),fmt0) 'allom_d2h1 = ',prt_params%allom_d2h1
        write(fates_log(),fmt0) 'allom_d2h2 = ',prt_params%allom_d2h2
        write(fates_log(),fmt0) 'allom_d2h3 = ',prt_params%allom_d2h3
        write(fates_log(),fmt0) 'allom_d2bl1 = ',prt_params%allom_d2bl1
        write(fates_log(),fmt0) 'allom_d2bl2 = ',prt_params%allom_d2bl2
        write(fates_log(),fmt0) 'allom_d2bl3 = ',prt_params%allom_d2bl3
        write(fates_log(),fmt0) 'allom_blca_expnt_diff = ',prt_params%allom_blca_expnt_diff
        write(fates_log(),fmt0) 'allom_d2ca_coefficient_max = ',prt_params%allom_d2ca_coefficient_max
        write(fates_log(),fmt0) 'allom_d2ca_coefficient_min = ',prt_params%allom_d2ca_coefficient_min        
        write(fates_log(),fmt0) 'allom_agb1 = ',prt_params%allom_agb1
        write(fates_log(),fmt0) 'allom_agb2 = ',prt_params%allom_agb2
        write(fates_log(),fmt0) 'allom_agb3 = ',prt_params%allom_agb3
        write(fates_log(),fmt0) 'allom_agb4 = ',prt_params%allom_agb4
        write(fates_log(),fmt0) 'prt_nitr_stoich_p1 = ',prt_params%nitr_stoich_p1
        write(fates_log(),fmt0) 'prt_nitr_stoich_p2 = ',prt_params%nitr_stoich_p2
        write(fates_log(),fmt0) 'prt_phos_stoich_p1 = ',prt_params%phos_stoich_p1
        write(fates_log(),fmt0) 'prt_phos_stoich_p2 = ',prt_params%phos_stoich_p2
        write(fates_log(),fmt0) 'prt_alloc_priority = ',prt_params%alloc_priority
        write(fates_log(),fmt0) 'woody = ',prt_params%woody
        write(fates_log(),fmt0) 'roota_par = ',prt_params%fnrt_prof_a
        write(fates_log(),fmt0) 'rootb_par = ',prt_params%fnrt_prof_b
        write(fates_log(),fmt0) 'fnrt_prof_mode = ',prt_params%fnrt_prof_mode
        write(fates_log(),fmt0) 'turnover_carb_retrans = ',prt_params%turnover_carb_retrans
        write(fates_log(),fmt0) 'turnover_nitr_retrans = ',prt_params%turnover_nitr_retrans
        write(fates_log(),fmt0) 'turnover_phos_retrans = ',prt_params%turnover_phos_retrans
        write(fates_log(),*) '-------------------------------------------------'

     end if

  end subroutine FatesReportPFTParams

  ! =====================================================================================
  
  subroutine PRTCheckParams(is_master)

     ! ----------------------------------------------------------------------------------
     !
     ! This subroutine performs logical checks on user supplied parameters.  It cross
     ! compares various parameters and will fail if they don't make sense.  
     ! Examples:
     ! A tree can not be defined as both evergreen and deciduous.  A woody plant
     ! cannot have a structural biomass allometry intercept of 0, and a non-woody
     ! plant (grass) can't have a non-zero intercept...
     ! -----------------------------------------------------------------------------------


     ! Argument
     logical, intent(in) :: is_master    ! Only log if this is the master proc

     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft     ! number of PFTs
     integer :: ipft     ! pft index
     integer :: nleafage ! size of the leaf age class array
     integer :: iage     ! leaf age class index
     integer :: norgans  ! size of the plant organ dimension
     integer :: i, io    ! generic loop index and organ loop index

     
     integer, parameter,dimension(6) :: cnpflex_organs = &
          [leaf_organ, fnrt_organ, sapw_organ, store_organ, repro_organ, struct_organ]

     
     npft = size(prt_params%evergreen,1)

     ! Prior to performing checks copy grperc to the 
     ! organ dimensioned version

     norgans = size(prt_params%nitr_stoich_p1,2)

     if(.not.is_master) return




     if (norgans .ne. num_organ_types) then
        write(fates_log(),*) 'The size of the organ dimension for PRT parameters'
        write(fates_log(),*) 'as specified in the parameter file is incompatible.'
        write(fates_log(),*) 'All currently acceptable hypothesese are using'
        write(fates_log(),*) 'the full set of num_organ_types = ',num_organ_types
        write(fates_log(),*) 'The parameter file listed ',norgans
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if


     
     do ipft = 1,npft
        
        ! Check to see if evergreen, deciduous flags are mutually exclusive
        ! ----------------------------------------------------------------------------------

        if ( int(prt_params%evergreen(ipft) +    &
                 prt_params%season_decid(ipft) + &
                 prt_params%stress_decid(ipft)) .ne. 1 ) then
           
           write(fates_log(),*) 'PFT # ',ipft,' must be defined as having one of three'
           write(fates_log(),*) 'phenology habits, ie == 1'
           write(fates_log(),*) 'stress_decid: ',prt_params%stress_decid(ipft)
           write(fates_log(),*) 'season_decid: ',prt_params%season_decid(ipft)
           write(fates_log(),*) 'evergreen: ',prt_params%evergreen(ipft)
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if


        ! Check to see if mature and base seed allocation is greater than 1
        ! ----------------------------------------------------------------------------------
        if ( ( prt_params%seed_alloc(ipft) + &
               prt_params%seed_alloc_mature(ipft)) > 1.0_r8 ) then

           write(fates_log(),*) 'The sum of seed allocation from base and mature trees may'
           write(fates_log(),*) ' not exceed 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' seed_alloc: ',prt_params%seed_alloc(ipft)
           write(fates_log(),*) ' seed_alloc_mature: ',prt_params%seed_alloc_mature(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if

        ! Check if woody plants have a structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
        if ( ( prt_params%allom_agb1(ipft) <= tiny(prt_params%allom_agb1(ipft)) ) .and. &
             ( int(prt_params%woody(ipft)) .eq. 1 ) ) then

           write(fates_log(),*) 'Woody plants are expected to have a non-zero intercept'
           write(fates_log(),*) ' in the diameter to AGB allometry equations'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' allom_agb1: ',prt_params%allom_agb1(ipft)
           write(fates_log(),*) ' woody: ',int(prt_params%woody(ipft))
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if

        ! Check if non-woody plants have structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
!        if ( ( prt_params%allom_agb1(ipft) > tiny(prt_params%allom_agb1(ipft)) ) .and. &
!              ( int(prt_params%woody(ipft)) .ne. 1 ) ) then
!
!           write(fates_log(),*) 'Non-woody plants are expected to have a zero intercept'
!           write(fates_log(),*) ' in the diameter to AGB allometry equations'
!           write(fates_log(),*) ' This is because the definition of AGB (as far as allometry)'
!           write(fates_log(),*) ' is concerned, ignores leaf and fine-roots, and only contains'
!           write(fates_log(),*) ' woody tissues (sap and structural dead wood).'
!           write(fates_log(),*) ' PFT#: ',ipft
!           write(fates_log(),*) ' allom_agb1: ',prt_params%allom_agb1(ipft)
!           write(fates_log(),*) ' woody: ',int(prt_params%woody(ipft))
!           write(fates_log(),*) ' Aborting'
!           call endrun(msg=errMsg(sourcefile, __LINE__))
!
!        end if

        ! Check if leaf storage priority is between 0-1
        ! ----------------------------------------------------------------------------------
        
        if ( ( prt_params%leaf_stor_priority(ipft) < 0.0_r8 ) .or. &
             ( prt_params%leaf_stor_priority(ipft) > 1.0_r8 ) ) then

           write(fates_log(),*) 'Prioritization of carbon allocation to leaf'
           write(fates_log(),*) ' and root turnover replacement, must be between'
           write(fates_log(),*) ' 0 and 1'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) 'leaf_stor_priority: ',prt_params%leaf_stor_priority(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if


        ! Check re-translocations
        ! Seems reasonable to assume that sapwood, structure and reproduction
        ! should not be re-translocating mass upon turnover.
        ! Note to advanced users. Feel free to remove these checks...
        ! -------------------------------------------------------------------
        
        if ( (prt_params%turnover_carb_retrans(ipft,repro_organ) > nearzero) ) then
           write(fates_log(),*) ' Retranslocation of reproductive tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,repro_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,repro_organ) > nearzero) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,repro_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of reproductive tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,repro_organ)
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,repro_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,repro_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
           
        if ((prt_params%turnover_carb_retrans(ipft,sapw_organ) > nearzero)) then
           write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,sapw_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,sapw_organ) > nearzero) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,sapw_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        if ((prt_params%turnover_carb_retrans(ipft,struct_organ) > nearzero)) then
           write(fates_log(),*) ' Retranslocation of structural(dead) tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,struct_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,struct_organ) > nearzero) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,struct_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of structural(dead) tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,struct_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,struct_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Leaf retranslocation should be between 0 and 1
        if ( (prt_params%turnover_carb_retrans(ipft,leaf_organ) > 1.0_r8) .or. & 
             (prt_params%turnover_carb_retrans(ipft,leaf_organ) < 0.0_r8) ) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,leaf_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,leaf_organ) > 1.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,leaf_organ) > 1.0_r8) .or. &
               (prt_params%turnover_nitr_retrans(ipft,leaf_organ) < 0.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,leaf_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,leaf_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,leaf_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Fineroot retranslocation should be between 0-1
        if ((prt_params%turnover_carb_retrans(ipft,fnrt_organ) > 1.0_r8) .or. & 
            (prt_params%turnover_carb_retrans(ipft,fnrt_organ) < 0.0_r8)) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,fnrt_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,fnrt_organ) > 1.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,fnrt_organ) > 1.0_r8) .or. &
               (prt_params%turnover_nitr_retrans(ipft,fnrt_organ) < 0.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,fnrt_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,fnrt_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,fnrt_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Storage retranslocation should be between 0-1 (storage retrans seems weird, but who knows)
        if ((prt_params%turnover_carb_retrans(ipft,store_organ) > 1.0_r8) .or. & 
            (prt_params%turnover_carb_retrans(ipft,store_organ) < 0.0_r8)) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',prt_params%turnover_carb_retrans(ipft,store_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((prt_params%turnover_nitr_retrans(ipft,store_organ) > 1.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,store_organ) > 1.0_r8) .or. &
               (prt_params%turnover_nitr_retrans(ipft,store_organ) < 0.0_r8) .or.  & 
               (prt_params%turnover_phos_retrans(ipft,store_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,store_organ)
              write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,store_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Growth respiration
        !        if (parteh_mode .eq. prt_carbon_allom_hyp) then
        if ( ( prt_params%grperc(ipft) < 0.0_r8) .or. &
             ( prt_params%grperc(ipft) > 1.0_r8 ) ) then
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',prt_params%grperc(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
!        elseif(parteh_mode .eq. prt_cnp_flex_allom_hyp) then
!           if ( ( any(prt_params%grperc_organ(ipft,:) < 0.0_r8)) .or. &
!                ( any(prt_params%grperc_organ(ipft,:) >= 1.0_r8)) ) then
!              write(fates_log(),*) ' PFT#: ',ipft
!              write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',prt_params%grperc_organ(ipft,:)
!              write(fates_log(),*) ' Aborting'
!              call endrun(msg=errMsg(sourcefile, __LINE__))
!           end if
!        end if


        ! The first nitrogen stoichiometry is used in all cases
        if ( (any(prt_params%nitr_stoich_p1(ipft,:) < 0.0_r8)) .or. &
             (any(prt_params%nitr_stoich_p1(ipft,:) >= 1.0_r8))) then
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' N per C stoichiometry must bet between 0-1'
           write(fates_log(),*) prt_params%nitr_stoich_p1(ipft,:)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if


        if(hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then

           do i = 1,size(cnpflex_organs,dim=1)
              io = cnpflex_organs(i)
              if ( (prt_params%nitr_stoich_p1(ipft,io) < 0._r8) .or. &
                   (prt_params%nitr_stoich_p2(ipft,io) < 0._r8) .or. &
                   (prt_params%phos_stoich_p1(ipft,io) < 0._r8) .or. &
                   (prt_params%phos_stoich_p2(ipft,io) < 0._r8) .or. & 
                   (prt_params%nitr_stoich_p1(ipft,io) > 1._r8) .or. &
                   (prt_params%nitr_stoich_p2(ipft,io) > 1._r8) .or. &
                   (prt_params%phos_stoich_p1(ipft,io) > 1._r8) .or. &
                   (prt_params%phos_stoich_p2(ipft,io) > 1._r8) ) then
                 write(fates_log(),*) 'When the C,N,P allocation hypothesis with flexible'
                 write(fates_log(),*) 'stoichiometry is turned on (prt_cnp_flex_allom_hyp),'
                 write(fates_log(),*) 'all stoichiometries must be greater than or equal to zero,'
                 write(fates_log(),*) 'and less than 1 (probably way less than 1).'
                 write(fates_log(),*) 'Setting both p1 and p2 parameters to zero will turn'
                 write(fates_log(),*) 'off nutrient dynamics for the given species.'
                 write(fates_log(),*) 'You specified an organ/pft less than zero.'
                 write(fates_log(),*) 'PFT: ',ipft
                 write(fates_log(),*) 'organ index (see head of PRTGenericMod): ',io
                 write(fates_log(),*) 'nitr_stoich_p1: ',prt_params%nitr_stoich_p1(ipft,io)
                 write(fates_log(),*) 'nitr_stoich_p2: ',prt_params%phos_stoich_p1(ipft,io)
                 write(fates_log(),*) 'phos_stoich_p1: ',prt_params%nitr_stoich_p2(ipft,io)
                 write(fates_log(),*) 'phos_stoich_p2: ',prt_params%phos_stoich_p2(ipft,io)
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end do

           if ( any(prt_params%alloc_priority(ipft,:) < 0) .or. &
                any(prt_params%alloc_priority(ipft,:) > 6) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' Allocation priorities should be 0-6 for CNP flex hypothesis'
              write(fates_log(),*) prt_params%alloc_priority(ipft,:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if

        end if


        

        ! Check turnover time-scales
        
        nleafage = size(prt_params%leaf_long,dim=2)

        do iage = 1, nleafage

           if ( prt_params%leaf_long(ipft,iage)>nearzero ) then
              
              ! Check that leaf turnover doesn't exeed 1 day
              if ( (years_per_day / prt_params%leaf_long(ipft,iage)) > 1._r8 ) then
                 write(fates_log(),*) 'Leaf turnover time-scale is greater than 1 day!'
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',prt_params%leaf_long(ipft,iage),' [years]'
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
              
              ! Check to make sure that all other age-classes for this PFT also
              ! have non-zero entries, it wouldn't make sense otherwise
              if ( any(prt_params%leaf_long(ipft,:) <= nearzero) ) then
                 write(fates_log(),*) 'You specified a leaf_long that is zero or'
                 write(fates_log(),*) 'invalid for a particular age class.'
                 write(fates_log(),*) 'Yet, other age classes for this PFT are non-zero.'
                 write(fates_log(),*) 'this doesnt make sense.'
                 write(fates_log(),*) 'ipft = ',ipft
                 write(fates_log(),*) 'leaf_long(ipft,:) =  ',prt_params%leaf_long(ipft,:)
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if

           else
              if (prt_params%evergreen(ipft) .eq. itrue) then
                 write(fates_log(),*) 'You specified zero leaf turnover: '
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',prt_params%leaf_long(ipft,iage)
                 write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
                 write(fates_log(),*) 'that an evergreen would have leaf maintenance turnover'
                 write(fates_log(),*) 'disable this error if you are ok with this'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end if

        end do

        ! Check the turnover rates on the senescing leaf pool
        if ( prt_params%leaf_long(ipft,nleafage)>nearzero ) then
           
           ! Check that leaf turnover doesn't exeed 1 day
           if ( (years_per_day / &
                 (prt_params%leaf_long(ipft,nleafage) * &
                  prt_params%senleaf_long_fdrought(ipft))) > 1._r8 ) then
              write(fates_log(),*) 'Drought-senescent turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'leaf_long(ipft,nleafage)*senleaf_long_fdrought: ', &
                    prt_params%leaf_long(ipft,nleafage)*prt_params%senleaf_long_fdrought(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        if ( prt_params%senleaf_long_fdrought(ipft)<nearzero .or. &
             prt_params%senleaf_long_fdrought(ipft)>1._r8 ) then
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) must be greater than 0 '
           write(fates_log(),*) 'or less than or equal to 1.'
           write(fates_log(),*) 'Set this to 1 if you want no accelerated senescence turnover'
           write(fates_log(),*) 'ipft = ',ipft
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) = ',prt_params%senleaf_long_fdrought(ipft)
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
           

        if ( prt_params%root_long(ipft)>nearzero ) then
           
           ! Check that root turnover doesn't exeed 1 day
           if ( (years_per_day / prt_params%root_long(ipft)) > 1._r8 ) then
              write(fates_log(),*) 'Root turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',prt_params%root_long(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           
        else
           if (prt_params%evergreen(ipft) .eq. itrue) then
              write(fates_log(),*) 'You specified zero root turnover: '
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',prt_params%root_long(ipft)
              write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
              write(fates_log(),*) 'that an evergreen would have root maintenance turnover'
              write(fates_log(),*) 'disable this error if you are ok with this'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Check Branch turnover doesn't exceed one day
        if ( prt_params%branch_long(ipft)>nearzero ) then
           
           ! Check that branch turnover doesn't exeed 1 day
           if ( (years_per_day / prt_params%branch_long(ipft)) > 1._r8 ) then
              write(fates_log(),*) 'Branch turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'branch_long(ipft): ',prt_params%branch_long(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if


     end do


     return
   end subroutine PRTCheckParams

  
 end module PRTInitParamsFatesMod
