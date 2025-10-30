module PRTInitParamsFatesMod

  ! This is a FATES specific module for loading parameters through
  ! the CLM/ELM module system.

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : itrue
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
  use FatesConstantsMod,   only : ievergreen
  use FatesConstantsMod,   only : isemi_stress_decid
  use JSONParameterUtilsMod,only : params_type,param_type
  
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PRTCheckParams
  public :: PRTDerivedParams
  public :: TransferParamsPRT
  public :: NewRecruitTotalStoichiometry
  !-----------------------------------------------------------------------

contains

  ! =====================================================================================

  subroutine TransferParamsPRT(pstruct)
    
    type(params_type) :: pstruct         ! Data structure containing all parameters and dimensions
    type(param_type),pointer :: param_p  ! Pointer to one specific parameter

    integer                  :: num_pft
    integer                  :: num_organ
    integer                  :: num_ageclass

    num_ageclass = pstruct%GetDimSizeFromName('fates_leafage_class')
    num_organ = pstruct%GetDimSizeFromName('fates_plant_organs')
    num_pft   = pstruct%GetDimSizeFromName('fates_pft')
    
    param_p => pstruct%GetParamFromName('fates_alloc_organ_id')
    allocate(prt_params%organ_id(num_organ))
    prt_params%organ_id(:) = param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_leaf_habit')
    allocate(prt_params%phen_leaf_habit(num_pft))
    prt_params%phen_leaf_habit(:) = param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_stem_drop_fraction')
    allocate(prt_params%phen_stem_drop_fraction(num_pft))
    prt_params%phen_stem_drop_fraction(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_fnrt_drop_fraction')
    allocate(prt_params%phen_fnrt_drop_fraction(num_pft))
    prt_params%phen_fnrt_drop_fraction(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_mindaysoff')
    allocate(prt_params%phen_doff_time(num_pft))
    prt_params%phen_doff_time(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_drought_threshold')
    allocate(prt_params%phen_drought_threshold(num_pft))
    prt_params%phen_drought_threshold(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_phen_moist_threshold')
    allocate(prt_params%phen_moist_threshold(num_pft))
    prt_params%phen_moist_threshold(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_leaf_slamax')
    allocate(prt_params%slamax(num_pft))
    prt_params%slamax(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_slatop')
    allocate(prt_params%slatop(num_pft))
    prt_params%slatop(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_sai_scaler')
    allocate(prt_params%allom_sai_scaler(num_pft))
    prt_params%allom_sai_scaler(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_fnrt_prof_a')
    allocate(prt_params%fnrt_prof_a(num_pft))
    prt_params%fnrt_prof_a(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_fnrt_prof_b')
    allocate(prt_params%fnrt_prof_b(num_pft))
    prt_params%fnrt_prof_b(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_fnrt_prof_mode')
    allocate(prt_params%fnrt_prof_mode(num_pft))
    prt_params%fnrt_prof_mode(:) = param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_woody')
    allocate(prt_params%woody(num_pft))
    prt_params%woody(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_wood_density')
    allocate(prt_params%wood_density(num_pft))
    prt_params%wood_density(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_seed_dbh_repro_threshold')
    allocate(prt_params%dbh_repro_threshold(num_pft))
    prt_params%dbh_repro_threshold(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_alloc_storage_cushion')
    allocate(prt_params%cushion(num_pft))
    prt_params%cushion(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_alloc_store_priority_frac')
    allocate(prt_params%leaf_stor_priority(num_pft))
    prt_params%leaf_stor_priority(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_turnover_senleaf_fdrought')
    allocate(prt_params%(num_pft))
    prt_params%senleaf_long_fdrought(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_turnover_fnrt')
    allocate(prt_params%(num_pft))
    prt_params%root_long(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_leafn_vert_scaler_coeff1')
    allocate(prt_params%(num_pft))
    prt_params%leafn_vert_scaler_coeff1(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_leafn_vert_scaler_coeff2')
    allocate(prt_params%(num_pft))
    prt_params%leafn_vert_scaler_coeff2(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_seed_alloc_mature')
    allocate(prt_params%(num_pft))
    prt_params%seed_alloc_mature(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_seed_alloc')
    allocate(prt_params%(num_pft))
    prt_params%seed_alloc(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_trs_repro_alloc_a')
    allocate(prt_params%(num_pft))
    prt_params%repro_alloc_a(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_trs_repro_alloc_b')
    allocate(prt_params%(num_pft))
    prt_params%repro_alloc_b(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_c2b')
    allocate(prt_params%(num_pft))
    prt_params%c2b(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_grperc')
    allocate(prt_params%(num_pft))
    prt_params%grperc(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_dbh_maxheight')
    allocate(prt_params%(num_pft))
    prt_params%allom_dbh_maxheight(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_hmode')
    allocate(prt_params%allom_hmode(num_pft))
    prt_params%allom_hmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_lmode')
    allocate(prt_params%allom_lmode(num_pft))
    prt_params%allom_lmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_fmode')
    allocate(prt_params%allom_fmode(num_pft))
    prt_params%allom_fmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_amode')
    allocate(prt_params%allom_amode(num_pft))
    prt_params%allom_amode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_stmode')
    allocate(prt_params%allom_stmode(num_pft))
    prt_params%allom_stmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_cmode')
    allocate(prt_params%allom_cmode(num_pft))
    prt_params%allom_cmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_smode')
    allocate(prt_params%allom_smode(num_pft))
    prt_params%allom_smode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_dmode')
    allocate(prt_params%allom_dmode(num_pft))
    prt_params%allom_dmode(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_la_per_sa_int')
    allocate(prt_params%allom_la_per_sa_int(num_pft))
    prt_params%allom_la_per_sa_int(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_la_per_sa_slp')
    allocate(prt_params%allom_la_per_sa_slp(num_pft))
    prt_params%allom_la_per_sa_slp(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_l2fr')
    allocate(prt_params%allom_l2fr(num_pft))
    prt_params%allom_l2fr(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_cnp_pid_kp')
    allocate(prt_params%pid_kp(num_pft))
    prt_params%pid_kp(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_cnp_pid_ki')
    allocate(prt_params%pid_ki(num_pft))
    prt_params%pid_ki(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_cnp_pid_kd')
    allocate(prt_params%pid_kd(num_pft))
    prt_params%pid_kd(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_store_ovrflw_frac')
    allocate(prt_params%store_ovrflw_frac(num_pft))
    prt_params%store_ovrflw_frac(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_nfix1')
    allocate(prt_params%nfix_mresp_scfrac(num_pft))
    prt_params%nfix_mresp_scfrac(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_agb_frac')
    allocate(prt_params%allom_agb_frac(num_pft))
    prt_params%allom_agb_frac(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2h1')
    allocate(prt_params%allom_d2h1(num_pft))
    prt_params%allom_d2h1(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2h2')
    allocate(prt_params%allom_d2h2(num_pft))
    prt_params%allom_d2h2(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2h3')
    allocate(prt_params%allom_d2h3(num_pft))
    prt_params%allom_d2h3(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2bl1')
    allocate(prt_params%allom_d2bl1(num_pft))
    prt_params%allom_d2bl1(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2bl2')
    allocate(prt_params%allom_d2bl2(num_pft))
    prt_params%allom_d2bl2(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2bl3')
    allocate(prt_params%allom_d2bl3(num_pft))
    prt_params%allom_d2bl3(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_blca_expnt_diff')
    allocate(prt_params%allom_blca_expnt_diff(num_pft))
    prt_params%allom_blca_expnt_diff(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2ca_coefficient_max')
    allocate(prt_params%allom_d2ca_coefficient_max(num_pft))
    prt_params%allom_d2ca_coefficient_max(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_d2ca_coefficient_min')
    allocate(prt_params%allom_d2ca_coefficient_min(num_pft))
    prt_params%allom_d2ca_coefficient_min(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_agb1')
    allocate(prt_params%allom_agb1(num_pft))
    prt_params%allom_agb1(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_agb2')
    allocate(prt_params%allom_agb2(num_pft))
    prt_params%allom_agb2(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_agb3')
    allocate(prt_params%allom_agb3(num_pft))
    prt_params%allom_agb3(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_agb4')
    allocate(prt_params%allom_agb4(num_pft))
    prt_params%allom_agb4(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_h2cd1')
    allocate(prt_params%allom_h2cd1(num_pft))
    prt_params%allom_h2cd1(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_h2cd2')
    allocate(prt_params%allom_h2cd2(num_pft))
    prt_params%allom_h2cd2(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_allom_zroot_max_dbh')
    allocate(prt_params%allom_zroot_max_dbh(num_pft))
    prt_params%allom_zroot_max_dbh(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_zroot_max_z')
    allocate(prt_params%allom_zroot_max_z(num_pft))
    prt_params%allom_zroot_max_z(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_zroot_min_dbh')
    allocate(prt_params%allom_zroot_min_dbh(num_pft))
    prt_params%allom_zroot_min_dbh(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_zroot_min_z')
    allocate(prt_params%allom_zroot_min_z(num_pft))
    prt_params%allom_zroot_min_z(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_zroot_k')
    allocate(prt_params%allom_zroot_k(num_pft))
    prt_params%allom_zroot_k(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_turnover_branch')
    allocate(prt_params%branch_long(num_pft))
    prt_params%branch_long(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_nitr_store_ratio')
    allocate(prt_params%nitr_store_ratio(num_pft))
    prt_params%nitr_store_ratio(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_phos_store_ratio')
    allocate(prt_params%phos_store_ratio(num_pft))
    prt_params%phos_store_ratio(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_turnover_leaf_canopy')
    allocate(prt_params%leaf_long(num_pft,num_ageclass))
    prt_params%leaf_long(:,:) = param_p%r_data_2d(:)

    param_p => pstruct%GetParamFromName('fates_turnover_leaf_ustory')
    allocate(prt_params%leaf_long_ustory(num_pft,num_ageclass))
    prt_params%leaf_long_ustory(:,:) = param_p%r_data_2d(:)

    param_p => pstruct%GetParamFromName('fates_stoich_nitr')
    allocate(prt_params%nitr_stoich_p1(num_pft,num_organ))
    prt_params%nitr_stoich_p1(:) = param_p%r_data_2d(:)
    
    param_p => pstruct%GetParamFromName('fates_stoich_phos')
    allocate(prt_params%phos_stoich_p1(num_pft,num_organ))
    prt_params%phos_stoich_p1(:) = param_p%r_data_2d(:)
    
    param_p => pstruct%GetParamFromName('fates_alloc_organ_priority')
    allocate(prt_params%alloc_priority(num_pft,num_organ))
    prt_params%alloc_priority(:) = param_p%i_data_2d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_turnover_nitr_retrans')
    allocate(prt_params%turnover_nitr_retrans(num_pft,num_organ))
    prt_params%turnover_nitr_retrans(:) = param_p%r_data_2d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_turnover_phos_retrans')
    allocate(prt_params%turnover_phos_retrans(num_pft,num_organ))
    prt_params%turnover_phos_retrans(:) = param_p%r_data_2d(:)

    return
  end subroutine TransferParamsPRT

  ! =====================================================================================
  
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
        write(fates_log(),fmti) 'phen_leaf_habit = ',prt_params%phen_leaf_habit
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
    npft    = size(prt_params%phen_leaf_habit,1)
    
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
     ! A woody plant cannot have a structural biomass allometry intercept of 0, and a 
     ! non-woody plant (grass) can't have a non-zero intercept...
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
     logical :: is_hmode_fine   ! Did the height allometry pass the check?
     integer :: nerror          ! Count number of errors. If this is not
                                !    zero by theend of the subroutine, stop 
                                !    the run.
     real(r8) :: height_crit    ! Critical height where crown depth equals height
     real(r8) :: height_max     ! Maximum height attainable by PFT.



     npft = size(prt_params%phen_leaf_habit,1)

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

        ! When using the the drought semi-deciduous phenology, we must ensure that the lower
        ! and upper thresholds are consistent (i.e., that both are based on either soil
        ! water content or soil matric potential).
        if (prt_params%phen_leaf_habit(ipft) == isemi_stress_decid) then
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
              write(fates_log(),*) ' phen_leaf_habit              = ',prt_params%phen_leaf_habit(ipft)
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
              write(fates_log(),*) ' phen_leaf_habit              = ',prt_params%phen_leaf_habit(ipft)
              write(fates_log(),*) ' fates_phen_drought_threshold = ',prt_params%phen_drought_threshold(ipft)
              write(fates_log(),*) ' fates_phen_moist_threshold   = ',prt_params%phen_moist_threshold  (ipft)
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
           end if
        end if

        ! For all deciduous PFTs, check that abscission fractions are all bounded.
        if (prt_params%phen_leaf_habit(ipft) /= ievergreen) then
           ! Check if the fraction of fine roots to be actively abscised relative to leaf abscission
           ! is bounded between 0 and 1 (exactly 0 and 1 are acceptable).
           if ( ( prt_params%phen_fnrt_drop_fraction(ipft) < 0.0_r8 ) .or. &
                ( prt_params%phen_fnrt_drop_fraction(ipft) > 1.0_r8 ) ) then
              write(fates_log(),*) '---~---'
              write(fates_log(),*) ' Abscission rate for fine roots must be between 0 and 1 for '
              write(fates_log(),*) ' deciduous PFTs.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' phen_leaf_habit: ',prt_params%phen_leaf_habit(ipft)
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
              write(fates_log(),*) ' phen_leaf_habit: ',prt_params%phen_leaf_habit(ipft)
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

           elseif (prt_params%phen_leaf_habit(ipft) == ievergreen) then
              write(fates_log(),*) "---~---"
              write(fates_log(),*) 'You specified zero leaf turnover: '
              write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
              write(fates_log(),*) 'phen_leaf_habit: ',prt_params%phen_leaf_habit(ipft)
              write(fates_log(),*) 'leaf_long(ipft,iage): ',prt_params%leaf_long(ipft,iage)
              write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
              write(fates_log(),*) 'that an evergreen would have leaf maintenance turnover'
              write(fates_log(),*) 'disable this error if you are ok with this'
              write(fates_log(),*) "---~---"
              write(fates_log(),*) ''
              write(fates_log(),*) ''
              nerror = nerror + 1
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
        elseif (prt_params%phen_leaf_habit(ipft) == ievergreen) then
           write(fates_log(),*) "---~---"
           write(fates_log(),*) 'You specified zero root turnover: '
           write(fates_log(),*) 'ipft: ',ipft
           write(fates_log(),*) 'phen_leaf_habit: ',prt_params%phen_leaf_habit(ipft)
           write(fates_log(),*) 'root_long(ipft): ',prt_params%root_long(ipft)
           write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
           write(fates_log(),*) 'that an evergreen would have root maintenance turnover'
           write(fates_log(),*) 'disable this error if you are ok with this'
           write(fates_log(),*) "---~---"
           write(fates_log(),*) ''
           write(fates_log(),*) ''
           nerror = nerror + 1
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
