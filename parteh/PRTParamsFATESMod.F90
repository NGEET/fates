module PRTInitParamsFatesMod

  ! This is a FATES specific module for loading parametes through
  ! the CLM/ELM module system.

  
  use PRTParametersMod, only: prt_params
  use PRTGenericMod,  only : num_organ_types
  use PRTGenericMod,  only : leaf_organ, fnrt_organ, store_organ
  use PRTGenericMod,  only : sapw_organ, struct_organ, repro_organ
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log 
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  
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

    name = 'fates_roota_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_rootb_par'
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

    name = 'fates_roota_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_a)

    name = 'fates_rootb_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_b)
    
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

    name = 'fates_rootprof_beta'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

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

    name = 'fates_rootprof_beta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=prt_params%rootprof_beta)

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
        write(fates_log(),fmt0) 'rootprof_beta = ',prt_params%rootprof_beta
        write(fates_log(),fmt0) 'turnover_carb_retrans = ',prt_params%turnover_carb_retrans
        write(fates_log(),fmt0) 'turnover_nitr_retrans = ',prt_params%turnover_nitr_retrans
        write(fates_log(),fmt0) 'turnover_phos_retrans = ',prt_params%turnover_phos_retrans
        write(fates_log(),*) '-------------------------------------------------'

     end if

  end subroutine FatesReportPFTParams

end module PRTInitParamsFatesMod
