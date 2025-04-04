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
  use PRTGenericMod,     only : nitrogen_element, phosphorus_element
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log 
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use EDPftvarcon,       only : EDPftvarcon_inst
  use PRTGenericMod,     only : prt_cnp_flex_allom_hyp,prt_carbon_allom_hyp
  use FatesAllometryMod  , only : h_allom
  use FatesAllometryMod  , only : h2d_allom
  use FatesAllometryMod  , only : bagw_allom
  use FatesAllometryMod  , only : bsap_allom
  use FatesAllometryMod  , only : bleaf
  use FatesAllometryMod  , only : bfineroot
  use FatesAllometryMod  , only : bdead_allom
  use FatesAllometryMod  , only : bstore_allom
  use FatesAllometryMod  , only : bbgw_allom
  use FatesAllometryMod  , only : carea_allom
  use FatesAllometryMod  , only : CheckIntegratedAllometries
  use FatesAllometryMod, only : set_root_fraction
  use PRTGenericMod, only : StorageNutrientTarget
  use EDTypesMod,          only : init_recruit_trim
  use FatesConstantsMod,   only : ihard_stress_decid, isemi_stress_decid
  
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
  public :: PRTDerivedParams
  public :: NewRecruitTotalStoichiometry
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
    call PRTRegisterOrgan(fates_params)
    
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
    call PRTReceiveOrgan(fates_params)
    
  end subroutine PRTReceiveParams

  ! =====================================================================================
  
  subroutine PRTRegisterOrgan(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_prt_organs, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_prt_organs/)
    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_general /)
    character(len=param_string_length) :: name

    name = 'fates_alloc_organ_id'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
  end subroutine PRTRegisterOrgan
    
  ! =====================================================================================

  subroutine PRTReceiveOrgan(fates_params)

    ! Make sure to call this after PRTRegisterPFTOrgans
    
    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    real(r8), allocatable :: tmpreal(:)  ! Temporary variable to hold floats
    
    name = 'fates_alloc_organ_id'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%organ_id(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%organ_id)
    deallocate(tmpreal)
    
  end subroutine PRTReceiveOrgan
  
  ! =====================================================================================
  
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

    name = 'fates_phen_stem_drop_fraction'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_fnrt_drop_fraction'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_mindaysoff'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_drought_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_moist_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_fnrt_prof_a'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_fnrt_prof_b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_fnrt_prof_mode'
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
    
    name = 'fates_recruit_seed_dbh_repro_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alloc_store_priority_frac'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_senleaf_fdrought'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_fnrt'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leafn_vert_scaler_coeff1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leafn_vert_scaler_coeff2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_recruit_seed_alloc_mature'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_recruit_seed_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_trs_repro_alloc_a'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_trs_repro_alloc_b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_c2b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_l2fr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_pid_kd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_pid_ki'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_pid_kp'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_cnp_store_ovrflw_frac'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_nfix1'
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

    name = 'fates_allom_dmode'
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

    name = 'fates_allom_h2cd1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_h2cd2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_zroot_max_dbh'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_zroot_max_z'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_zroot_min_dbh'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_zroot_min_z'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_zroot_k'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
        
    name = 'fates_turnover_branch'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

   

    name = 'fates_cnp_nitr_store_ratio'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_phos_store_ratio'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
  end subroutine PRTRegisterPFT
  
  !-----------------------------------------------------------------------
  
  subroutine PRTReceivePFT(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    real(r8), allocatable :: tmpreal(:)  ! Temporary variable to hold floats
                                         ! that are converted to ints

    name = 'fates_phen_stress_decid'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%stress_decid(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%stress_decid)
    deallocate(tmpreal)
    
    name = 'fates_phen_season_decid'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%season_decid(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%season_decid)
    deallocate(tmpreal)
    
    name = 'fates_phen_evergreen'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%evergreen(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%evergreen)
    deallocate(tmpreal)

    name = 'fates_phen_stem_drop_fraction'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%phen_stem_drop_fraction)

    name = 'fates_phen_fnrt_drop_fraction'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%phen_fnrt_drop_fraction)

    name = 'fates_phen_mindaysoff'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%phen_doff_time)

    name = 'fates_phen_drought_threshold'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%phen_drought_threshold)

    name = 'fates_phen_moist_threshold'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%phen_moist_threshold)

    name = 'fates_leaf_slamax'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%slamax)
    
    name = 'fates_leaf_slatop'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%slatop)

    name = 'fates_allom_sai_scaler'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_sai_scaler)

    name = 'fates_allom_fnrt_prof_a'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_a)

    name = 'fates_allom_fnrt_prof_b'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_b)

    name = 'fates_allom_fnrt_prof_mode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%fnrt_prof_mode)

    name = 'fates_woody'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%woody(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%woody)
    deallocate(tmpreal)
    
    name = 'fates_wood_density'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%wood_density)
    
    name = 'fates_recruit_seed_dbh_repro_threshold'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%dbh_repro_threshold)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%cushion)

    name = 'fates_alloc_store_priority_frac'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%leaf_stor_priority)

    name = 'fates_turnover_senleaf_fdrought'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%senleaf_long_fdrought)

    name = 'fates_turnover_fnrt'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%root_long)

    name = 'fates_leafn_vert_scaler_coeff1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%leafn_vert_scaler_coeff1)

    name = 'fates_leafn_vert_scaler_coeff2'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%leafn_vert_scaler_coeff2)
    
    name = 'fates_recruit_seed_alloc_mature'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%seed_alloc_mature)
    
    name = 'fates_recruit_seed_alloc'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%seed_alloc)

    name = 'fates_trs_repro_alloc_a'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%repro_alloc_a)

    name = 'fates_trs_repro_alloc_b'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%repro_alloc_b)

    name = 'fates_c2b'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%c2b)

    name = 'fates_grperc'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%grperc)

    name = 'fates_allom_dbh_maxheight'
    call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%allom_dbh_maxheight)

    name = 'fates_allom_hmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_hmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_hmode)
    deallocate(tmpreal)

    name = 'fates_allom_lmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_lmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_lmode)
    deallocate(tmpreal)

    name = 'fates_allom_fmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_fmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_fmode)
    deallocate(tmpreal)

    name = 'fates_allom_amode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_amode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_amode)
    deallocate(tmpreal)

    name = 'fates_allom_stmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_stmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_stmode)
    deallocate(tmpreal)

    name = 'fates_allom_cmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_cmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_cmode)
    deallocate(tmpreal)

    name = 'fates_allom_smode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_smode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_smode)
    deallocate(tmpreal)

    name = 'fates_allom_dmode'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(prt_params%allom_dmode(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,prt_params%allom_dmode)
    deallocate(tmpreal)

    name = 'fates_allom_la_per_sa_int'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_la_per_sa_int)

    name = 'fates_allom_la_per_sa_slp'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_la_per_sa_slp)

    name = 'fates_allom_l2fr'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_l2fr)

    name = 'fates_cnp_pid_kp'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%pid_kp)

    name = 'fates_cnp_pid_ki'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%pid_ki)

    name = 'fates_cnp_pid_kd'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%pid_kd)
    
    name = 'fates_cnp_store_ovrflw_frac'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%store_ovrflw_frac)
    
    name = 'fates_cnp_nfix1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%nfix_mresp_scfrac)
    
    name = 'fates_allom_agb_frac'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_agb_frac)

    name = 'fates_allom_d2h1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2h1)

    name = 'fates_allom_d2h2'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2h2)

    name = 'fates_allom_d2h3'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2h3)

    name = 'fates_allom_d2bl1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl1)

    name = 'fates_allom_d2bl2'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl2)

    name = 'fates_allom_d2bl3'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2bl3)

    name = 'fates_allom_blca_expnt_diff'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_blca_expnt_diff)

    name = 'fates_allom_d2ca_coefficient_max'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2ca_coefficient_max)

    name = 'fates_allom_d2ca_coefficient_min'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_d2ca_coefficient_min)

    name = 'fates_allom_agb1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_agb1)

    name = 'fates_allom_agb2'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_agb2)

    name = 'fates_allom_agb3'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_agb3)

    name = 'fates_allom_agb4'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_agb4)

    name = 'fates_allom_h2cd1'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_h2cd1)

    name = 'fates_allom_h2cd2'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_h2cd2)

    name = 'fates_allom_zroot_max_dbh'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_zroot_max_dbh)
    
    name = 'fates_allom_zroot_max_z'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_zroot_max_z)
    
    name = 'fates_allom_zroot_min_dbh'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_zroot_min_dbh)
    
    name = 'fates_allom_zroot_min_z'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_zroot_min_z)
    
    name = 'fates_allom_zroot_k'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%allom_zroot_k)
    
    name = 'fates_turnover_branch'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%branch_long)
    
    name = 'fates_cnp_nitr_store_ratio'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%nitr_store_ratio)
    
    name = 'fates_cnp_phos_store_ratio'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=prt_params%phos_store_ratio)

    
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

    name = 'fates_turnover_leaf_canopy'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_leaf_ustory'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    return
  end subroutine PRTRegisterPFTLeafAge

  ! =====================================================================================

  subroutine ArrayNint(realarr,intarr)

    real(r8),intent(in)  :: realarr(:)
    integer,intent(out)  :: intarr(:)
    integer  :: i

    do i = 1,size(realarr,dim=1)
       intarr(i) = nint(realarr(i))
    end do
    
    return
  end subroutine ArrayNint
  
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
    !X!    call fates_params%RetrieveParameter(name=name, &
    !X!         data=this%)

  end subroutine Receive_PFT_nvariants
  
  ! =====================================================================================

  subroutine PRTReceivePFTLeafAge(fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name

     name = 'fates_turnover_leaf_canopy'
     call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%leaf_long)

     name = 'fates_turnover_leaf_ustory'
     call fates_params%RetrieveParameterAllocate(name=name, &
          data=prt_params%leaf_long_ustory)

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

    name = 'fates_stoich_nitr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_stoich_phos'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alloc_organ_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_cnp_turnover_nitr_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_cnp_turnover_phos_retrans'
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

     name = 'fates_stoich_nitr'
     call fates_params%RetrieveParameterAllocate(name=name, &
           data=prt_params%nitr_stoich_p1)

     name = 'fates_stoich_phos'
     call fates_params%RetrieveParameterAllocate(name=name, &
           data=prt_params%phos_stoich_p1)

     name = 'fates_alloc_organ_priority'
     call fates_params%RetrieveParameterAllocate(name=name, &
           data=prt_params%alloc_priority)

     name = 'fates_cnp_turnover_nitr_retrans'
     call fates_params%RetrieveParameterAllocate(name=name, &
           data=prt_params%turnover_nitr_retrans)

     name = 'fates_cnp_turnover_phos_retrans'
     call fates_params%RetrieveParameterAllocate(name=name, &
           data=prt_params%turnover_phos_retrans)

  end subroutine PRTReceivePFTOrgans

  ! ===============================================================================================
  
  subroutine FatesReportPFTParams(is_master)
     
     ! Argument
     logical, intent(in) :: is_master  ! Only log if this is the master proc

     logical, parameter :: debug_report = .false.
     character(len=15),parameter :: fmti = '(a,100(I12,1X))'
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
        write(fates_log(),fmti) 'stress_decid = ',prt_params%stress_decid
        write(fates_log(),fmti) 'season_decid = ',prt_params%season_decid
        write(fates_log(),fmti) 'evergreen = ',prt_params%evergreen
        write(fates_log(),fmt0) 'phen_fnrt_drop_fraction = ',prt_params%phen_fnrt_drop_fraction
        write(fates_log(),fmt0) 'phen_stem_drop_fraction = ',prt_params%phen_stem_drop_fraction
        write(fates_log(),fmt0) 'phen_doff_time = ',prt_params%phen_doff_time
        write(fates_log(),fmt0) 'phen_drought_threshold = ',prt_params%phen_drought_threshold
        write(fates_log(),fmt0) 'phen_moist_threshold = ',prt_params%phen_moist_threshold
        write(fates_log(),fmt0) 'wood_density = ',prt_params%wood_density
        write(fates_log(),fmt0) 'dbh max height = ',prt_params%allom_dbh_maxheight
        write(fates_log(),fmt0) 'dbh mature = ',prt_params%dbh_repro_threshold
        write(fates_log(),fmt0) 'cushion = ',prt_params%cushion
        write(fates_log(),fmt0) 'leaf_stor_priority = ',prt_params%leaf_stor_priority
        write(fates_log(),fmt0) 'root_long = ',prt_params%root_long
        write(fates_log(),fmt0) 'senleaf_long_fdrought = ',prt_params%senleaf_long_fdrought
        write(fates_log(),fmt0) 'seed_alloc_mature = ',prt_params%seed_alloc_mature
        write(fates_log(),fmt0) 'seed_alloc = ',prt_params%seed_alloc
        write(fates_log(),fmt0) 'repro_alloc_a = ',prt_params%repro_alloc_a
        write(fates_log(),fmt0) 'repro_alloc_b = ',prt_params%repro_alloc_b
        write(fates_log(),fmt0) 'slamax = ',prt_params%slamax
        write(fates_log(),fmt0) 'slatop = ',prt_params%slatop
        write(fates_log(),fmt0) 'allom_sai_scaler = ',prt_params%allom_sai_scaler
        write(fates_log(),fmt0) 'leaf_long = ',prt_params%leaf_long
        write(fates_log(),fmt0) 'leaf_long_ustory = ',prt_params%leaf_long_ustory
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
        write(fates_log(),fmt0) 'pid_kp = ',prt_params%pid_kp
        write(fates_log(),fmt0) 'pid_ki = ',prt_params%pid_ki
        write(fates_log(),fmt0) 'pid_kd = ',prt_params%pid_kd
        write(fates_log(),fmt0) 'store_ovrflw_frac = ',prt_params%store_ovrflw_frac
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
        write(fates_log(),fmt0) 'allom_h2cd1 = ',prt_params%allom_h2cd1
        write(fates_log(),fmt0) 'allom_h2cd2 = ',prt_params%allom_h2cd2

        write(fates_log(),fmt0) 'allom_zroot_max_dbh = ',prt_params%allom_zroot_max_dbh
        write(fates_log(),fmt0) 'allom_zroot_max_z = ',prt_params%allom_zroot_max_z
        write(fates_log(),fmt0) 'allom_zroot_min_dbh = ',prt_params%allom_zroot_min_dbh
        write(fates_log(),fmt0) 'allom_zroot_min_z = ',prt_params%allom_zroot_min_z
        write(fates_log(),fmt0) 'allom_zroot_k = ',prt_params%allom_zroot_k
        
        write(fates_log(),fmt0) 'stoich_nitr = ',prt_params%nitr_stoich_p1
        write(fates_log(),fmt0) 'stoich_phos = ',prt_params%phos_stoich_p1
        write(fates_log(),fmt0) 'alloc_organ_priority = ',prt_params%alloc_priority
        write(fates_log(),fmt0) 'woody = ',prt_params%woody
        write(fates_log(),fmt0) 'roota_par = ',prt_params%fnrt_prof_a
        write(fates_log(),fmt0) 'rootb_par = ',prt_params%fnrt_prof_b
        write(fates_log(),fmt0) 'fnrt_prof_mode = ',prt_params%fnrt_prof_mode
        write(fates_log(),fmt0) 'turnover_nitr_retrans = ',prt_params%turnover_nitr_retrans
        write(fates_log(),fmt0) 'turnover_phos_retrans = ',prt_params%turnover_phos_retrans
        write(fates_log(),fmti) 'organ_id = ',prt_params%organ_id
        write(fates_log(),fmt0) 'nitr_store_ratio = ',prt_params%nitr_store_ratio
        write(fates_log(),fmt0) 'phos_store_ratio = ',prt_params%phos_store_ratio
        write(fates_log(),fmt0) 'leafn_vert_scaler_coeff1 = ',prt_params%leafn_vert_scaler_coeff1
        write(fates_log(),fmt0) 'leafn_vert_scaler_coeff2 = ',prt_params%leafn_vert_scaler_coeff2
        write(fates_log(),*) '-------------------------------------------------'

     end if

  end subroutine FatesReportPFTParams

  ! =====================================================================================

  subroutine PRTDerivedParams()

    integer :: npft     ! number of PFTs
    integer :: ft       ! pft index
    integer :: norgans  ! number of organs in the parameter file
    integer :: i, io    ! generic loop index and organ loop index
    
    norgans = size(prt_params%organ_id,1)
    npft    = size(prt_params%evergreen,1)
    
    ! Set the reverse lookup map for organs to the parameter file index
    allocate(prt_params%organ_param_id(num_organ_types))
    
    ! Initialize them as invalid
    prt_params%organ_param_id(:) = -1
    
    do i = 1,norgans
       prt_params%organ_param_id(prt_params%organ_id(i)) = i
    end do

    return
  end subroutine PRTDerivedParams
    
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

     integer :: npft            ! number of PFTs
     integer :: ipft            ! pft index
     integer :: nleafage        ! size of the leaf age class array
     integer :: iage            ! leaf age class index
     integer :: norgans         ! size of the plant organ dimension
     integer :: i, io           ! generic loop index and organ loop index
     logical :: is_evergreen    ! Is the PFT evergreen
     logical :: is_season_decid ! Is the PFT cold-deciduous?
     logical :: is_stress_decid ! Is the PFT drought-deciduous?
     logical :: is_semi_decid   ! Is the PFT drought semi-deciduous?
     logical :: is_hmode_fine   ! Did the height allometry pass the check?
     integer :: nerror          ! Count number of errors. If this is not
                                !    zero by theend of the subroutine, stop 
                                !    the run.
     real(r8) :: height_crit    ! Critical height where crown depth equals height
     real(r8) :: height_max     ! Maximum height attainable by PFT.



     npft = size(prt_params%evergreen,1)

     ! Prior to performing checks copy grperc to the 
     ! organ dimensioned version

     norgans = size(prt_params%organ_id,1)

     if(.not.is_master) return

     ! Initialise nerror with zero. If anything is incorrectly set, nerror will be 
     ! positive, but we will hold on until all checks are performed before stopping
     ! the run.
     nerror = 0

     ! Initialise height allometry success flag to .true., and update it if there are
     ! inconsistencies.
     is_hmode_fine = .true.

     if( any(prt_params%organ_id(:)<1) .or. &
         any(prt_params%organ_id(:)>num_organ_types) ) then
        write(fates_log(),*) '---~---'
        write(fates_log(),*) 'prt_organ_ids should match the global ids'
        write(fates_log(),*) 'of organ types found in PRTGenericMod.F90'
        write(fates_log(),*) 'organ_ids: ',prt_params%organ_id(:)
        write(fates_log(),*) '---~---'
        write(fates_log(),*) ''
        write(fates_log(),*) ''
        nerror = nerror + 1
     end if

     ! Check to make sure the organ ids are valid if this is the
     ! cnp_flex_allom_hypothesis
     select case (hlm_parteh_mode)
     case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp)

         do io = 1,norgans
           if(prt_params%organ_id(io) == repro_organ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) 'with flexible cnp or c-only alloc hypotheses'
              write(fates_log(),*) 'reproductive tissues are a special case'
              write(fates_log(),*) 'and therefore should not be included in'
              write(fates_log(),*) 'the parameter file organ list'
              write(fates_log(),*) 'fates_prt_organ_id: ',prt_params%organ_id(:)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
           if(prt_params%organ_id(io) == store_organ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) 'with flexible cnp or c-only alloc hypotheses'
              write(fates_log(),*) 'storage is a special case'
              write(fates_log(),*) 'and therefore should not be included in'
              write(fates_log(),*) 'the parameter file organ list'
              write(fates_log(),*) 'fates_prt_organ_id: ',prt_params%organ_id(:)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if

        end do
     end select

     ! Make sure that the N fixation respiration surcharge fraction is
     ! between 0 and 1
     if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
        if(any(prt_params%nfix_mresp_scfrac(:)<0._r8) .or. any(prt_params%nfix_mresp_scfrac(:)>1.0_r8)) then
           write(fates_log(),*) '---~---'
           write(fates_log(),*) 'The N fixation surcharge nfix_mresp_sfrac (fates_nfix1) must be between 0-1.'
           write(fates_log(),*) 'here are the values: ',prt_params%nfix_mresp_scfrac(:)
           write(fates_log(),*) '---~---'
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end if
     end if

     
     pftloop: do ipft = 1,npft

        ! Check to see if evergreen, deciduous flags are mutually exclusive
        !      By the way, if these are mutually exclusive, shouldn't we define a
        !      single prt_params%leaf_phenology and a list of codes for the different
        !      types (i.e., ievergreen, iseason_decid, istress_hard, istress_semi, etc.)?
        ! ----------------------------------------------------------------------------------
        is_evergreen    = prt_params%evergreen(ipft)    == itrue
        is_season_decid = prt_params%season_decid(ipft) == itrue
        is_stress_decid = any(prt_params%stress_decid(ipft) == [ihard_stress_decid,isemi_stress_decid])
        is_semi_decid   = prt_params%stress_decid(ipft) == isemi_stress_decid

        if ( ( is_evergreen    .and. is_season_decid ) .or. &
             ( is_evergreen    .and. is_stress_decid ) .or. &
             ( is_season_decid .and. is_stress_decid ) ) then

           write(fates_log(),*) '---~---'
           write(fates_log(),*) 'PFT # ',ipft,' must be defined as having one of three'
           write(fates_log(),*) 'phenology habits, ie, only one of the flags below should'
           write(fates_log(),*) 'be different than ',ifalse
           write(fates_log(),*) 'stress_decid: ',prt_params%stress_decid(ipft)
           write(fates_log(),*) 'season_decid: ',prt_params%season_decid(ipft)
           write(fates_log(),*) 'evergreen: ',prt_params%evergreen(ipft)
           write(fates_log(),*) '---~---'
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end if


        ! When using the the drought semi-deciduous phenology, we must ensure that the lower
        ! and upper thresholds are consistent (i.e., that both are based on either soil
        ! water content or soil matric potential).
        if (is_semi_decid) then
           if ( prt_params%phen_drought_threshold(ipft)*prt_params%phen_moist_threshold(ipft) < 0._r8 ) then
              ! In case the product of the lower and upper thresholds is negative, the
              !    thresholds are inconsistent as both should be defined using the same 
              !    quantity.
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' When using drought semi-deciduous phenology,'
              write(fates_log(),*) '    the moist threshold must have the same sign as'
              write(fates_log(),*) '    the dry threshold.  Positive = soil water content [m3/m3],'
              write(fates_log(),*) '    Negative = soil matric potential [mm].'
              write(fates_log(),*) ' PFT                          = ',ipft
              write(fates_log(),*) ' Stress_decid                 = ',prt_params%stress_decid(ipft)
              write(fates_log(),*) ' fates_phen_drought_threshold = ',prt_params%phen_drought_threshold(ipft)
              write(fates_log(),*) ' fates_phen_moist_threshold   = ',prt_params%phen_moist_threshold  (ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1

           elseif ( prt_params%phen_drought_threshold(ipft) >= prt_params%phen_moist_threshold(ipft) ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' When using drought semi-deciduous phenology,'
              write(fates_log(),*) '   the moist threshold must be greater than the dry threshold.'
              write(fates_log(),*) '   By greater we mean more positive or less negative, and'
              write(fates_log(),*) '   they cannot be the identical.'
              write(fates_log(),*) ' PFT                          = ',ipft
              write(fates_log(),*) ' Stress_decid                 = ',prt_params%stress_decid(ipft)
              write(fates_log(),*) ' fates_phen_drought_threshold = ',prt_params%phen_drought_threshold(ipft)
              write(fates_log(),*) ' fates_phen_moist_threshold   = ',prt_params%phen_moist_threshold  (ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if

        ! For all deciduous PFTs, check that abscission fractions are all bounded.
        if (prt_params%evergreen(ipft) == ifalse) then
           ! Check if the fraction of fine roots to be actively abscised relative to leaf abscission
           ! is bounded between 0 and 1 (exactly 0 and 1 are acceptable).
           if ( ( prt_params%phen_fnrt_drop_fraction(ipft) < 0.0_r8 ) .or. &
                ( prt_params%phen_fnrt_drop_fraction(ipft) > 1.0_r8 ) ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' Abscission rate for fine roots must be between 0 and 1 for '
              write(fates_log(),*) ' deciduous PFTs.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' evergreen flag: (should be 0):',prt_params%evergreen(ipft)
              write(fates_log(),*) ' phen_fnrt_drop_fraction: ', prt_params%phen_fnrt_drop_fraction(ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if


           ! Check if the fraction of stems to be actively abscised relative to leaf abscission
           ! is bounded between 0 and 1 (exactly 0 and 1 are acceptable) when PFTs are non-woody, and
           ! that the fraction is zero for woody PFTs.  Stem abscission is a solution to avoid hydraulic
           ! problems when plant hydraulics is enabled.
           if ( ( prt_params%woody(ipft) == itrue )                     .and. &
                ( ( prt_params%phen_stem_drop_fraction(ipft) < 0.0_r8 ) .or.  &
                  ( prt_params%phen_stem_drop_fraction(ipft) > nearzero )   ) ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' Non-zero stem-drop fractions are not allowed for woody plants'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' part_params%woody:',prt_params%woody(ipft)
              write(fates_log(),*) ' phen_stem_drop_fraction: ', prt_params%phen_stem_drop_fraction(ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           elseif ( ( prt_params%phen_stem_drop_fraction(ipft) < 0.0_r8 ) .or. &
                    ( prt_params%phen_stem_drop_fraction(ipft) > 1.0_r8 ) ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' Deciduous non-wood plants must keep 0-100% of their stems'
              write(fates_log(),*) ' during the deciduous period.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' evergreen flag: (should be 0):',prt_params%evergreen(ipft)
              write(fates_log(),*) ' phen_stem_drop_fraction: ', prt_params%phen_stem_drop_fraction(ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if




        ! Check to see if mature and base seed allocation is greater than 1
        ! ----------------------------------------------------------------------------------
        if ( ( prt_params%seed_alloc(ipft) + &
               prt_params%seed_alloc_mature(ipft)) > 1.0_r8 ) then

           write(fates_log(),*) '---~---'
           write(fates_log(),*) 'The sum of seed allocation from base and mature trees may'
           write(fates_log(),*) ' not exceed 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' seed_alloc: ',prt_params%seed_alloc(ipft)
           write(fates_log(),*) ' seed_alloc_mature: ',prt_params%seed_alloc_mature(ipft)
           write(fates_log(),*) '---~---'
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end if

        ! Check if woody plants have a structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
        if ( ( prt_params%allom_agb1(ipft) <= tiny(prt_params%allom_agb1(ipft)) ) .and. &
             ( prt_params%woody(ipft) .eq. 1 ) ) then

           write(fates_log(),*) '---~---'
           write(fates_log(),*) 'Woody plants are expected to have a non-zero intercept'
           write(fates_log(),*) ' in the diameter to AGB allometry equations'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' allom_agb1: ',prt_params%allom_agb1(ipft)
           write(fates_log(),*) ' woody: ',prt_params%woody(ipft)
           write(fates_log(),*) '---~---'
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1

        end if

        ! Check if parameter 2 for dbh -> height is negative when allom_hmode is 2
        ! (Weibull function / Poorter et al. (2006)
        ! ML: FATES definition for parameter 2 is a bit unusual, which is why I added
        !     the check. Normally the minus sign is left outside the parameter for 
        !     Weibull functions.
        ! ----------------------------------------------------------------------------------
        if ( ( prt_params%allom_hmode(ipft) == 2     ) .and. &
             ( prt_params%allom_d2h2 (ipft) >  0._r8 ) ) then
           write(fates_log(),*) "---~---"
           write(fates_log(),*) " Incorrect settings for height allometry."
           write(fates_log(),*) ' PFT index:   ',ipft
           write(fates_log(),*) ' allom_hmode: ',prt_params%allom_hmode(ipft)
           write(fates_log(),*) ' allom_d2h2:  ',prt_params%allom_d2h2 (ipft)
           write(fates_log(),*) " Parameter ""allom_d2h2"" must be negative when using"
           write(fates_log(),*) "    allom_hmode = 2."
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
           
           ! Update flag so we do not run tests that depend on reasonable height allometry
           ! parameters. It is fine to bypass these additional tests because the code will
           ! stop due to the height allometry error.
           ! -------------------------------------------------------------------------------
           is_hmode_fine = .false.
        end if

        ! Make sure that the crown depth does not exceed plant height.
        ! ----------------------------------------------------------------------------------
        select_dmode_check: select case (prt_params%allom_dmode(ipft))
        case (1)
           ! Linear allometry
           ! -------------------------------------------------------------------------------
           if ( ( prt_params%allom_h2cd1 (ipft) <  nearzero ) .or. &
                ( prt_params%allom_h2cd1 (ipft) >  1._r8    ) ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) " Incorrect settings for crown depth allometry."
              write(fates_log(),*) ' PFT index:   ',ipft
              write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
              write(fates_log(),*) ' allom_h2cd1: ',prt_params%allom_h2cd1(ipft)
              write(fates_log(),*) " Parameter ""allom_h2cd1"" must be positive and <= 1"
              write(fates_log(),*) "    when allom_dmode = 1 (or allom_h2cd2 = 1)."
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        case (2)
           ! Log-linear allometry. Multiplier factor cannot be negative or zero.
           ! -------------------------------------------------------------------------------
           if ( prt_params%allom_h2cd1 (ipft) < nearzero ) then
              !   Calculations for the generic case require allom_h2cd1 to be positive. If
              ! not, issue an error.
              ! ----------------------------------------------------------------------------
              write(fates_log(),*) "---~---"
              write(fates_log(),*) " Incorrect settings for crown depth allometry."
              write(fates_log(),*) ' PFT index:   ',ipft
              write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
              write(fates_log(),*) ' allom_h2cd1: ',prt_params%allom_h2cd1(ipft)
              write(fates_log(),*) " Parameter ""allom_h2cd1"" must be positive when"
              write(fates_log(),*) "    allom_dmode = 2."
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           elseif ( abs(prt_params%allom_h2cd2 (ipft) - 1._r8 ) < nearzero ) then
              !    Special case when the log-linear equation reduces to linear. This must
              ! be checked separately to avoid singularities in the general case.
              ! ----------------------------------------------------------------------------
              if ( ( prt_params%allom_h2cd1 (ipft) < nearzero ) .or. &
                   ( prt_params%allom_h2cd1 (ipft) > 1._r8    ) ) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) " Incorrect settings for crown depth allometry."
                 write(fates_log(),*) ' PFT index:   ',ipft
                 write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
                 write(fates_log(),*) ' allom_h2cd1: ',prt_params%allom_h2cd1(ipft)
                 write(fates_log(),*) " Parameter ""allom_h2cd1"" must be positive and <= 1"
                 write(fates_log(),*) "    when allom_h2cd2 = 1 (or allom_dmode = 1)."
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if
           elseif (is_hmode_fine) then
              ! ----------------------------------------------------------------------------
              !    General log-linear case. Depending on the parameter values, crown depth 
              ! may exceed for very small or very large plants. The code has safeguards 
              ! to prevent this behaviour, but we must at least issue some warning.
              ! ----------------------------------------------------------------------------

              !   Find the critical height in which crown depth becomes height
              ! ----------------------------------------------------------------------------
              height_crit = prt_params%allom_h2cd1 (ipft) ** &
                           ( 1.0_r8 / (1.0_r8 - prt_params%allom_h2cd2 (ipft)) )

              !   Find the maximum height.
              ! ----------------------------------------------------------------------------
              call h_allom(prt_params%allom_dbh_maxheight(ipft),ipft,height_max)

              if ( ( prt_params%allom_h2cd2 (ipft)  <  1.0_r8     ) .and. &
                   ( EDPftvarcon_inst%hgt_min(ipft) < height_crit ) ) then
                 !   These parameters will cause the code to cap crown depth to height for
                 ! small plants. We print a warning message, but we do not stop the run.
                 ! -------------------------------------------------------------------------
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) "   WARNING!"
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) " Parameter settings will require capping crown"
                 write(fates_log(),*) "    depth to height for cohorts with height less"
                 write(fates_log(),*) "    than ""height_crit""."
                 write(fates_log(),*) " "
                 write(fates_log(),*) ' PFT index:   ',ipft
                 write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
                 write(fates_log(),*) ' allom_h2cd1: ',prt_params%allom_h2cd1(ipft)
                 write(fates_log(),*) ' allom_h2cd2: ',prt_params%allom_h2cd2(ipft)
                 write(fates_log(),*) ' height_crit: ',height_crit
                 write(fates_log(),*) ' height_min:  ',EDPftvarcon_inst%hgt_min(ipft)
                 write(fates_log(),*) " "
                 write(fates_log(),*) " To avoid this message, set ""allom_h2cd1"" and"
                 write(fates_log(),*) "    ""allom_h2cd2"" such that ""height_crit"" is"
                 write(fates_log(),*) "    less than ""height_min""."
                 write(fates_log(),*) " "
                 write(fates_log(),*) " height_crit = allom_h2cd1**(1/(1-allom_h2cd2))"
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ""
                 write(fates_log(),*) ""
              elseif ( ( prt_params%allom_h2cd2 (ipft) > 1.0_r8      ) .and. &
                       ( height_max                    > height_crit ) ) then
                 !   These parameters will cause the code to cap crown depth to height for
                 ! large plants. We print a warning message, but we do not stop the run.
                 ! -------------------------------------------------------------------------
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) "   WARNING!"
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) " Parameter settings will require capping crown"
                 write(fates_log(),*) "    depth to height for cohorts with height greater"
                 write(fates_log(),*) "    than ""height_crit""."
                 write(fates_log(),*) " "
                 write(fates_log(),*) ' PFT index:   ',ipft
                 write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
                 write(fates_log(),*) ' allom_h2cd1: ',prt_params%allom_h2cd1(ipft)
                 write(fates_log(),*) ' allom_h2cd2: ',prt_params%allom_h2cd2(ipft)
                 write(fates_log(),*) ' height_crit: ',height_crit
                 write(fates_log(),*) ' height_max:  ',height_max
                 write(fates_log(),*) " "
                 write(fates_log(),*) " To avoid this message, set ""allom_h2cd1"" and"
                 write(fates_log(),*) "    ""allom_h2cd2"" such that ""height_crit"" is"
                 write(fates_log(),*) "    greater than ""height_max""."
                 write(fates_log(),*) " "
                 write(fates_log(),*) " height_crit = allom_h2cd1**(1/(1-allom_h2cd2))"
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ""
                 write(fates_log(),*) ""
              end if
           end if
        case default
           write(fates_log(),*) "---~---"
           write(fates_log(),*) " Incorrect settings for crown depth allometry."
           write(fates_log(),*) ' PFT index:   ',ipft
           write(fates_log(),*) ' allom_dmode: ',prt_params%allom_dmode(ipft)
           write(fates_log(),*) " Parameter ""allom_dmode"" must be 1 or 2."
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end select select_dmode_check

        ! Check if non-woody plants have structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
!        if ( ( prt_params%allom_agb1(ipft) > tiny(prt_params%allom_agb1(ipft)) ) .and. &
!              ( iprt_params%woody(ipft) .ne. 1 ) ) then
!
!           write(fates_log(),*) "---~---"
!           write(fates_log(),*) 'Non-woody plants are expected to have a zero intercept'
!           write(fates_log(),*) ' in the diameter to AGB allometry equations'
!           write(fates_log(),*) ' This is because the definition of AGB (as far as allometry)'
!           write(fates_log(),*) ' is concerned, ignores leaf and fine-roots, and only contains'
!           write(fates_log(),*) ' woody tissues (sap and structural dead wood).'
!           write(fates_log(),*) ' PFT#: ',ipft
!           write(fates_log(),*) ' allom_agb1: ',prt_params%allom_agb1(ipft)
!           write(fates_log(),*) ' woody: ',prt_params%woody(ipft)
!           write(fates_log(),*) "---~---"
!           write(fates_log(),*) ''
!           write(fates_log(),*) ''
!           nerror = nerror + 1
!
!        end if

        ! Check if leaf storage priority is between 0-1
        ! ----------------------------------------------------------------------------------
        
        if ( ( prt_params%leaf_stor_priority(ipft) < 0.0_r8 ) .or. &
             ( prt_params%leaf_stor_priority(ipft) > 1.0_r8 ) ) then

           write(fates_log(),*) "---~---"
           write(fates_log(),*) 'Prioritization of carbon allocation to leaf'
           write(fates_log(),*) ' and root turnover replacement, must be between'
           write(fates_log(),*) ' 0 and 1'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) 'leaf_stor_priority: ',prt_params%leaf_stor_priority(ipft)
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1

        end if

        select case (hlm_parteh_mode)
        case (prt_cnp_flex_allom_hyp)

           ! Make sure nutrient storage fractions are positive
           if( prt_params%nitr_store_ratio(ipft) < 0._r8  ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'With parteh allometric CNP hypothesis'
              write(fates_log(),*) 'nitr_store_ratio must be > 0'
              write(fates_log(),*) 'PFT#: ',ipft
              write(fates_log(),*) 'nitr_store_ratio = ',prt_params%nitr_store_ratio(ipft)
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
           if( prt_params%phos_store_ratio(ipft) < 0._r8 ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'With parteh allometric CNP hypothesis'
              write(fates_log(),*) 'phos_store_ratio must be > 0'
              write(fates_log(),*) 'PFT#: ',ipft
              write(fates_log(),*) 'nitr_store_ratio = ',prt_params%phos_store_ratio(ipft)
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if

           do i = 1,norgans
              io = prt_params%organ_id(i)

              if(io == sapw_organ) then
                 if ((prt_params%turnover_nitr_retrans(ipft,i) > nearzero)) then
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
                    write(fates_log(),*) ' PFT#: ',ipft
                    write(fates_log(),*) ' nitrogen retrans: ',prt_params%turnover_nitr_retrans(ipft,i)
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ''
                    write(fates_log(),*) ''
                    nerror = nerror + 1
                 end if
                 if ((prt_params%turnover_phos_retrans(ipft,i) > nearzero)) then
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
                    write(fates_log(),*) ' PFT#: ',ipft
                    write(fates_log(),*) ' phosphorus retrans: ',prt_params%turnover_nitr_retrans(ipft,i)
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ''
                    write(fates_log(),*) ''
                    nerror = nerror + 1
                 end if
              elseif(io == struct_organ) then
                 if ((prt_params%turnover_nitr_retrans(ipft,i) > nearzero)) then
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ' Retranslocation of structural tissues should be zero.'
                    write(fates_log(),*) ' PFT#: ',ipft
                    write(fates_log(),*) ' carbon retrans: ',prt_params%turnover_nitr_retrans(ipft,i)
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ''
                    write(fates_log(),*) ''
                    nerror = nerror + 1
                 end if
                 if ((prt_params%turnover_phos_retrans(ipft,i) > nearzero)) then
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ' Retranslocation of structural tissues should be zero.'
                    write(fates_log(),*) ' PFT#: ',ipft
                    write(fates_log(),*) ' phosphorus retrans: ',prt_params%turnover_nitr_retrans(ipft,i)
                    write(fates_log(),*) "---~---"
                    write(fates_log(),*) ''
                    write(fates_log(),*) ''
                    nerror = nerror + 1
                 end if
              end if
              
              ! Otherwise, all other retranslocations should be between 0 and 1
              if ((prt_params%turnover_nitr_retrans(ipft,i) > 1.0_r8) .or.  & 
                   (prt_params%turnover_phos_retrans(ipft,i) > 1.0_r8) .or. &
                   (prt_params%turnover_nitr_retrans(ipft,i) < 0.0_r8) .or.  & 
                   (prt_params%turnover_phos_retrans(ipft,i) < 0.0_r8)) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ' Retranslocation should range from 0 to 1.'
                 write(fates_log(),*) ' PFT#: ',ipft
                 write(fates_log(),*) ' parameter file organ index: ',i,' global index: ',io
                 write(fates_log(),*) ' nitr: ',prt_params%turnover_nitr_retrans(ipft,i)
                 write(fates_log(),*) ' phos: ',prt_params%turnover_phos_retrans(ipft,i)
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if

           end do

        end select


        ! Growth respiration
        !        if (parteh_mode .eq. prt_carbon_allom_hyp) then
        if ( ( prt_params%grperc(ipft) < 0.0_r8) .or. &
             ( prt_params%grperc(ipft) > 1.0_r8 ) ) then
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',prt_params%grperc(ipft)
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end if
!        elseif(parteh_mode .eq. prt_cnp_flex_allom_hyp) then
!           if ( ( any(prt_params%grperc_organ(ipft,:) < 0.0_r8)) .or. &
!                ( any(prt_params%grperc_organ(ipft,:) >= 1.0_r8)) ) then
!              write(fates_log(),*) "---~---"
!              write(fates_log(),*) ' PFT#: ',ipft
!              write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',prt_params%grperc_organ(ipft,:)
!              write(fates_log(),*) "---~---"
!              write(fates_log(),*) ''
!              write(fates_log(),*) ''
!              nerror = nerror + 1
!           end if
!        end if

        select case (hlm_parteh_mode)
        case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp)
           ! The first nitrogen stoichiometry is used in all cases
           if ( (any(prt_params%nitr_stoich_p1(ipft,:) < 0.0_r8)) .or. &
                (any(prt_params%nitr_stoich_p1(ipft,:) >= 1.0_r8))) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' N per C stoichiometry must bet between 0-1'
              write(fates_log(),*) prt_params%nitr_stoich_p1(ipft,:)
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end select

        select case (hlm_parteh_mode)
        case (prt_cnp_flex_allom_hyp)

           do i = 1,norgans
              if ( (prt_params%nitr_stoich_p1(ipft,i) < 0._r8) .or. &
                   (prt_params%phos_stoich_p1(ipft,i) < 0._r8) .or. &
                   (prt_params%nitr_stoich_p1(ipft,i) > 1._r8) .or. &
                   (prt_params%phos_stoich_p1(ipft,i) > 1._r8) ) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) 'When the C,N,P allocation hypothesis with flexible'
                 write(fates_log(),*) 'stoichiometry is turned on (prt_cnp_flex_allom_hyp),'
                 write(fates_log(),*) 'all stoichiometries must be greater than or equal to zero,'
                 write(fates_log(),*) 'and less than 1 (probably way less than 1).'
                 write(fates_log(),*) 'You specified an organ/pft less than zero.'
                 write(fates_log(),*) 'PFT: ',ipft
                 write(fates_log(),*) 'organ index (see head of PRTGenericMod): ',io
                 write(fates_log(),*) 'nitr_stoich: ',prt_params%nitr_stoich_p1(ipft,i)
                 write(fates_log(),*) 'phos_stoich: ',prt_params%phos_stoich_p1(ipft,i)
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if
           end do

           if ( any(prt_params%alloc_priority(ipft,:) < 0) .or. &
                any(prt_params%alloc_priority(ipft,:) > 6) ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' Allocation priorities should be 0-6 for CNP flex hypothesis'
              write(fates_log(),*) prt_params%alloc_priority(ipft,:)
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if

        end select


        ! Check turnover time-scales
        
        nleafage = size(prt_params%leaf_long,dim=2)

        do iage = 1, nleafage

           if ( prt_params%leaf_long(ipft,iage)>nearzero ) then
              
              ! Check that leaf turnover doesn't exeed 1 day
              if ( (years_per_day / prt_params%leaf_long(ipft,iage)) > 1._r8 ) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) 'Leaf turnover time-scale is greater than 1 day!'
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',prt_params%leaf_long(ipft,iage),' [years]'
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if
              
              ! Check to make sure that all other age-classes for this PFT also
              ! have non-zero entries, it wouldn't make sense otherwise
              if ( any(prt_params%leaf_long(ipft,:) <= nearzero) ) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) 'You specified a leaf_long that is zero or'
                 write(fates_log(),*) 'invalid for a particular age class.'
                 write(fates_log(),*) 'Yet, other age classes for this PFT are non-zero.'
                 write(fates_log(),*) 'this doesnt make sense.'
                 write(fates_log(),*) 'ipft = ',ipft
                 write(fates_log(),*) 'leaf_long(ipft,:) =  ',prt_params%leaf_long(ipft,:)
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if

           else
              if (prt_params%evergreen(ipft) .eq. itrue) then
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) 'You specified zero leaf turnover: '
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',prt_params%leaf_long(ipft,iage)
                 write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
                 write(fates_log(),*) 'that an evergreen would have leaf maintenance turnover'
                 write(fates_log(),*) 'disable this error if you are ok with this'
                 write(fates_log(),*) "---~---"
                 write(fates_log(),*) ''
                 write(fates_log(),*) ''
                 nerror = nerror + 1
              end if
           end if

        end do

        ! Check the turnover rates on the senescing leaf pool
        if ( prt_params%leaf_long(ipft,nleafage)>nearzero ) then
           
           ! Check that leaf turnover doesn't exeed 1 day
           if ( (years_per_day / &
                 (prt_params%leaf_long(ipft,nleafage) * &
                  prt_params%senleaf_long_fdrought(ipft))) > 1._r8 ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'Drought-senescent turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'leaf_long(ipft,nleafage)*senleaf_long_fdrought: ', &
                    prt_params%leaf_long(ipft,nleafage)*prt_params%senleaf_long_fdrought(ipft),' [years]'
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if
        
        if ( prt_params%senleaf_long_fdrought(ipft)<nearzero .or. &
             prt_params%senleaf_long_fdrought(ipft)>1._r8 ) then
           write(fates_log(),*) "---~---"
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) must be greater than 0 '
           write(fates_log(),*) 'or less than or equal to 1.'
           write(fates_log(),*) 'Set this to 1 if you want no accelerated senescence turnover'
           write(fates_log(),*) 'ipft = ',ipft
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) = ',prt_params%senleaf_long_fdrought(ipft)
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
        end if
           

        if ( prt_params%root_long(ipft)>nearzero ) then
           
           ! Check that root turnover doesn't exeed 1 day
           if ( (years_per_day / prt_params%root_long(ipft)) > 1._r8 ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'Root turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',prt_params%root_long(ipft),' [years]'
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
           
        else
           if (prt_params%evergreen(ipft) .eq. itrue) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'You specified zero root turnover: '
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',prt_params%root_long(ipft)
              write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
              write(fates_log(),*) 'that an evergreen would have root maintenance turnover'
              write(fates_log(),*) 'disable this error if you are ok with this'
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if
        
        ! Check Branch turnover doesn't exceed one day
        if ( prt_params%branch_long(ipft)>nearzero ) then
           
           ! Check that branch turnover doesn't exeed 1 day
           if ( (years_per_day / prt_params%branch_long(ipft)) > 1._r8 ) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'Branch turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'branch_long(ipft): ',prt_params%branch_long(ipft),' [years]'
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if


     end do pftloop


     ! If any error was found, abort. We add a single point to abort the run after all
     ! checks so users can get all the errors and address them in one go (as opposed to
     ! multiple submissions).
     if (nerror > 0) then
        write(fates_log(),*) 'One or more parameter errors found. Aborting.'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     return
   end subroutine PRTCheckParams

   ! ====================================================================================
   
   function NewRecruitTotalStoichiometry(ft,l2fr,element_id) result(recruit_stoich)

     ! ----------------------------------------------------------------------------------
     ! This function calculates the total N:C or P:C ratio for a newly recruited plant
     ! It does this by first identifying the dbh of a new plant, then uses
     ! allometry to calculate the starting amount of carbon, and then uses
     ! the stoichiometry parameters to determine the proportional mass of N or P
     !
     ! This process only has to be called once, and is then stored in parameter
     ! constants for each PFT.  These values are used for determining nutrient
     ! fluxes into seed pools (on plant), and also from germinated seed polls (on ground)
     ! into new recruits.
     ! ----------------------------------------------------------------------------------

     integer,intent(in)  :: ft
     integer,intent(in)  :: element_id
     real(r8),intent(in) :: l2fr
     real(r8)            :: recruit_stoich  ! nutrient to carbon ratio of recruit

     real(r8) :: dbh         ! dbh of the new recruit [cm]
     real(r8) :: c_leaf      ! target leaf biomass [kgC]
     real(r8) :: c_fnrt      ! target fine root biomass [kgC]
     real(r8) :: c_sapw      ! target sapwood biomass [kgC]
     real(r8) :: a_sapw      ! target sapwood cross section are [m2] (dummy)
     real(r8) :: c_agw       ! target Above ground biomass [kgC]
     real(r8) :: c_bgw       ! target Below ground biomass [kgC]
     real(r8) :: c_struct    ! target Structural biomass [kgc]
     real(r8) :: c_store     ! target Storage biomass [kgC]
     real(r8) :: c_total     ! total target carbon
     real(r8) :: nutr_total  ! total target nutrient

     integer, parameter :: not_damaged = 1 ! this is also in MainDamageMod, here for dependency purposes

     ! For all tissues, assume tissues of new recruits will be fully flushed.
     call h2d_allom(EDPftvarcon_inst%hgt_min(ft),ft,dbh)
     call bleaf(dbh,ft,not_damaged,init_recruit_trim, 1.0_r8, c_leaf)
     call bfineroot(dbh,ft,init_recruit_trim,l2fr, 1.0_r8,c_fnrt)
     call bsap_allom(dbh,ft,not_damaged,init_recruit_trim, 1.0_r8,a_sapw, c_sapw)
     call bagw_allom(dbh,ft,not_damaged, 1.0_r8,c_agw)
     call bbgw_allom(dbh,ft, 1.0_r8,c_bgw)
     call bdead_allom(c_agw,c_bgw,c_sapw,ft,c_struct)
     call bstore_allom(dbh,ft,not_damaged,init_recruit_trim,c_store)

     ! Total carbon in a newly recruited plant
     c_total = c_leaf + c_fnrt + c_sapw + c_struct + c_store

     ! Total nutrient in a newly recruited plant
     select case(element_id)
     case(nitrogen_element)

        nutr_total = &
             c_struct*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(struct_organ)) + &
             c_leaf*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)) + &
             c_fnrt*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)) + & 
             c_sapw*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)) + &
             StorageNutrientTarget(ft, element_id, &
             c_leaf*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)), &
             c_fnrt*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)), &
             c_sapw*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)), &
             c_struct*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(struct_organ)))

     case(phosphorus_element)

        nutr_total = &
             c_struct*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(struct_organ)) + &
             c_leaf*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)) + &
             c_fnrt*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)) + & 
             c_sapw*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)) + &
             StorageNutrientTarget(ft, element_id, &
             c_leaf*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)), &
             c_fnrt*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)), &
             c_sapw*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)), &
             c_struct*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(struct_organ)))


     end select

     recruit_stoich = nutr_total/c_total


     return
   end function NewRecruitTotalStoichiometry
  
 end module PRTInitParamsFatesMod
