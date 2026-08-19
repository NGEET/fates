module FatesUnitTestPRTParametersMod

  use FatesConstantsMod, only: r8 => fates_r8
  use FatesInterfaceTypesMod, only: numpft, hlm_use_planthydro, nleafage
  use EDParamsMod, only: hydr_htftype_node
  use PRTParametersMod, only: prt_params
  use FatesParameterDerivedMod, only : param_derived
  use EDPftvarconMockMod, only: init_mock_edpftvarcon

  implicit none

contains

  subroutine init_mock_prt_params(npft_in)
    integer, intent(in), optional :: npft_in
    integer :: npft, max_pft
    
    if (present(npft_in)) then
       npft = npft_in
    else
       npft = 1
    end if
    max_pft = npft
    
    numpft = npft
    hlm_use_planthydro = 1
    nleafage = 1
    
    if (allocated(prt_params%phen_leaf_habit)) then
       if (size(prt_params%phen_leaf_habit) /= npft) then
          call cleanup_mock_prt_params()
       end if
    end if
    
    if (.not. allocated(hydr_htftype_node)) then
      allocate(hydr_htftype_node(4))
      hydr_htftype_node(:) = 2
    end if
    
    if (.not. allocated(prt_params%phen_leaf_habit)) then
      allocate(prt_params%phen_leaf_habit(npft))
      allocate(prt_params%phen_fnrt_drop_fraction(npft))
      allocate(prt_params%phen_stem_drop_fraction(npft))
      allocate(prt_params%phen_drought_threshold(npft))
      allocate(prt_params%phen_moist_threshold(npft))
      allocate(prt_params%phen_doff_time(npft))
      allocate(prt_params%senleaf_long_fdrought(max_pft))
      allocate(prt_params%leaf_long(max_pft, max_pft))
      allocate(prt_params%leaf_long_ustory(max_pft, max_pft))
      allocate(prt_params%root_long(max_pft))
      allocate(prt_params%branch_long(max_pft))
      allocate(prt_params%turnover_nitr_retrans(max_pft, max_pft))
      allocate(prt_params%turnover_phos_retrans(max_pft, max_pft))
      allocate(prt_params%leafn_vert_scaler_coeff1(max_pft))
      allocate(prt_params%leafn_vert_scaler_coeff2(max_pft))
      allocate(prt_params%grperc(max_pft))
      allocate(prt_params%nitr_stoich_p1(max_pft, max_pft))
      allocate(prt_params%phos_stoich_p1(max_pft, max_pft))
      allocate(prt_params%nitr_store_ratio(max_pft))
      allocate(prt_params%phos_store_ratio(max_pft))
      allocate(prt_params%organ_id(max_pft))
      allocate(prt_params%alloc_priority(max_pft, max_pft))
      allocate(prt_params%cushion(max_pft))
      allocate(prt_params%leaf_stor_priority(max_pft))
      allocate(prt_params%dbh_repro_threshold(max_pft))
      allocate(prt_params%seed_alloc_mature(max_pft))
      allocate(prt_params%seed_alloc(max_pft))
      allocate(prt_params%repro_alloc_a(max_pft))
      allocate(prt_params%repro_alloc_b(max_pft))
      allocate(prt_params%organ_param_id(max_pft))
      allocate(prt_params%fnrt_prof_mode(max_pft))
      allocate(prt_params%fnrt_prof_a(max_pft))
      allocate(prt_params%fnrt_prof_b(max_pft))
      allocate(prt_params%c2b(max_pft))
      allocate(prt_params%wood_density(max_pft))
      allocate(prt_params%woody(max_pft))
      allocate(prt_params%slamax(max_pft))
      allocate(prt_params%slatop(max_pft))
      allocate(prt_params%allom_sai_scaler(max_pft))
      allocate(prt_params%allom_dbh_maxheight(max_pft))
      allocate(prt_params%allom_hmode(max_pft))
      allocate(prt_params%allom_lmode(max_pft))
      allocate(prt_params%allom_fmode(max_pft))
      allocate(prt_params%allom_amode(max_pft))
      allocate(prt_params%allom_cmode(max_pft))
      allocate(prt_params%allom_smode(max_pft))
      allocate(prt_params%allom_stmode(max_pft))
      allocate(prt_params%allom_dmode(max_pft))
      allocate(prt_params%allom_la_per_sa_int(max_pft))
      allocate(prt_params%allom_la_per_sa_slp(max_pft))
      allocate(prt_params%allom_l2fr(max_pft))
      allocate(prt_params%allom_agb_frac(max_pft))
      allocate(prt_params%allom_d2h1(max_pft))
      allocate(prt_params%allom_d2h2(max_pft))
      allocate(prt_params%allom_d2h3(max_pft))
      allocate(prt_params%allom_d2bl1(max_pft))
      allocate(prt_params%allom_d2bl2(max_pft))
      allocate(prt_params%allom_d2bl3(max_pft))
      allocate(prt_params%allom_blca_expnt_diff(max_pft))
      allocate(prt_params%allom_d2ca_coefficient_max(max_pft))
      allocate(prt_params%allom_d2ca_coefficient_min(max_pft))
      allocate(prt_params%allom_agb1(max_pft))
      allocate(prt_params%allom_agb2(max_pft))
      allocate(prt_params%allom_agb3(max_pft))
      allocate(prt_params%allom_agb4(max_pft))
      allocate(prt_params%allom_h2cd1(max_pft))
      allocate(prt_params%allom_h2cd2(max_pft))
      allocate(prt_params%allom_zroot_max_dbh(max_pft))
      allocate(prt_params%allom_zroot_max_z(max_pft))
      allocate(prt_params%allom_zroot_min_dbh(max_pft))
      allocate(prt_params%allom_zroot_min_z(max_pft))
      allocate(prt_params%allom_zroot_k(max_pft))
      allocate(prt_params%pid_kp(max_pft))
      allocate(prt_params%pid_ki(max_pft))
      allocate(prt_params%pid_kd(max_pft))
      allocate(prt_params%store_ovrflw_frac(max_pft))
      allocate(prt_params%nfix_mresp_scfrac(max_pft))
    end if

      ! --- PHENOLOGY DEFAULTS ---
      ! Evergreen habit provides continuous canopy leaf area without requiring seasonal driver triggers
      prt_params%phen_leaf_habit = 1
      ! Complete loss upon senescence isolates turnover rates from partial retention mechanics
      prt_params%phen_fnrt_drop_fraction = 1.0_r8
      prt_params%phen_stem_drop_fraction = 1.0_r8
      ! Unit thresholds normalize drought/moisture phenology scaling
      prt_params%phen_drought_threshold = 1.0_r8
      prt_params%phen_moist_threshold = 1.0_r8
      prt_params%phen_doff_time = 1.0_r8
      prt_params%senleaf_long_fdrought = 1.0_r8

      ! --- ORGAN LONGEVITY & TURNOVER ---
      ! 1-year baseline lifespan establishes standard annual turnover for temperate broadleaf vegetation
      prt_params%leaf_long = 1.0_r8
      prt_params%leaf_long_ustory = 1.0_r8
      prt_params%root_long = 1.0_r8
      prt_params%branch_long = 1.0_r8
      ! Complete nutrient retranslocation isolates carbon pool turnover from nutrient limitation feedback
      prt_params%turnover_nitr_retrans = 1.0_r8
      prt_params%turnover_phos_retrans = 1.0_r8

      ! --- CANOPY & STOICHIOMETRY SCALERS ---
      prt_params%leafn_vert_scaler_coeff1 = 1.0_r8
      prt_params%leafn_vert_scaler_coeff2 = 1.0_r8
      prt_params%grperc = 1.0_r8
      prt_params%nitr_stoich_p1 = 1.0_r8
      prt_params%phos_stoich_p1 = 1.0_r8
      prt_params%nitr_store_ratio = 1.0_r8
      prt_params%phos_store_ratio = 1.0_r8

      ! --- ALLOCATION & REPRODUCTION DEFAULTS ---
      prt_params%organ_id = 1
      prt_params%alloc_priority = 1
      prt_params%cushion = 1.0_r8
      prt_params%leaf_stor_priority = 1.0_r8
      prt_params%dbh_repro_threshold = 1.0_r8
      prt_params%seed_alloc_mature = 1.0_r8
      prt_params%seed_alloc = 1.0_r8
      prt_params%repro_alloc_a = 1.0_r8
      prt_params%repro_alloc_b = 1.0_r8
      prt_params%organ_param_id = 1

      ! --- ROOT PROFILE PARAMETERS ---
      ! Zeng 2001 exponential root profile mode (mode 1)
      prt_params%fnrt_prof_mode = 1.0_r8
      ! Exponential depth decay coefficients (Zeng 2001) establish realistic root density profile with depth
      prt_params%fnrt_prof_a = 7.0_r8
      prt_params%fnrt_prof_b = 2.0_r8

      ! --- BIOMASS & ALLOMETRY CONVERSION FACTORS ---
      ! 50% carbon content per unit dry plant biomass (2.0 g biomass / g C)
      prt_params%c2b = 2.0_r8
      ! Temperate broadleaf sapwood density default (0.5 g/cm3)
      prt_params%wood_density = 0.5_r8
      ! Woody plant functional type flag (1 = woody tree, 0 = non-woody)
      prt_params%woody = 1
      ! Specific leaf area bounds canopy light interception (SLA top/max)
      prt_params%slamax = 0.02_r8
      prt_params%slatop = 0.015_r8
      prt_params%allom_sai_scaler = 1.0_r8
      prt_params%allom_dbh_maxheight = 1.0_r8

      ! Allometry component modes (1 = standard FATES allometric equations)
      prt_params%allom_hmode = 1
      prt_params%allom_lmode = 1
      prt_params%allom_fmode = 1
      prt_params%allom_amode = 1
      prt_params%allom_cmode = 1
      prt_params%allom_smode = 1
      prt_params%allom_stmode = 1
      prt_params%allom_dmode = 1

      prt_params%allom_la_per_sa_int = 1.0_r8
      prt_params%allom_la_per_sa_slp = 1.0_r8
      prt_params%allom_l2fr = 1.0_r8
      ! 60% aboveground / 40% belowground biomass allocation ratio (Saldarriaga et al. 1988)
      prt_params%allom_agb_frac = 0.6_r8
      prt_params%allom_d2h1 = 1.0_r8
      prt_params%allom_d2h2 = 1.0_r8
      prt_params%allom_d2h3 = 1.0_r8
      prt_params%allom_d2bl1 = 0.1_r8
      prt_params%allom_d2bl2 = 1.0_r8
      prt_params%allom_d2bl3 = 1.0_r8
      prt_params%allom_blca_expnt_diff = 1.0_r8
      prt_params%allom_d2ca_coefficient_max = 1.0_r8
      prt_params%allom_d2ca_coefficient_min = 1.0_r8
      prt_params%allom_agb1 = 1.0_r8
      prt_params%allom_agb2 = 1.0_r8
      prt_params%allom_agb3 = 1.0_r8
      prt_params%allom_agb4 = 1.0_r8
      prt_params%allom_h2cd1 = 1.0_r8
      prt_params%allom_h2cd2 = 1.0_r8

      ! --- ROOTING DEPTH SCALING PARAMETERS ---
      ! Rooting depth scaling bounds established by Jackson et al. 1996
      prt_params%allom_zroot_max_dbh = 100.0_r8
      prt_params%allom_zroot_min_dbh = 1.0_r8
      prt_params%allom_zroot_max_z = 2.0_r8
      prt_params%allom_zroot_min_z = 0.5_r8
      prt_params%allom_zroot_k = 0.05_r8

      ! --- STORAGE & CONTROLLER DEFAULTS ---
      prt_params%pid_kp = 1.0_r8
      prt_params%pid_ki = 1.0_r8
      prt_params%pid_kd = 1.0_r8
      prt_params%store_ovrflw_frac = 1.0_r8
      prt_params%nfix_mresp_scfrac = 1.0_r8

    call init_mock_edpftvarcon()
    call param_derived%Init(npft)
  end subroutine init_mock_prt_params

  subroutine cleanup_mock_prt_params()
    ! Deallocate mock prt_params arrays to support clean test tearDown()
    if (allocated(hydr_htftype_node)) deallocate(hydr_htftype_node)

    if (allocated(prt_params%phen_leaf_habit)) then
      deallocate(prt_params%phen_leaf_habit)
      deallocate(prt_params%phen_fnrt_drop_fraction)
      deallocate(prt_params%phen_stem_drop_fraction)
      deallocate(prt_params%phen_drought_threshold)
      deallocate(prt_params%phen_moist_threshold)
      deallocate(prt_params%phen_doff_time)
      deallocate(prt_params%senleaf_long_fdrought)
      deallocate(prt_params%leaf_long)
      deallocate(prt_params%leaf_long_ustory)
      deallocate(prt_params%root_long)
      deallocate(prt_params%branch_long)
      deallocate(prt_params%turnover_nitr_retrans)
      deallocate(prt_params%turnover_phos_retrans)
      deallocate(prt_params%leafn_vert_scaler_coeff1)
      deallocate(prt_params%leafn_vert_scaler_coeff2)
      deallocate(prt_params%grperc)
      deallocate(prt_params%nitr_stoich_p1)
      deallocate(prt_params%phos_stoich_p1)
      deallocate(prt_params%nitr_store_ratio)
      deallocate(prt_params%phos_store_ratio)
      deallocate(prt_params%organ_id)
      deallocate(prt_params%alloc_priority)
      deallocate(prt_params%cushion)
      deallocate(prt_params%leaf_stor_priority)
      deallocate(prt_params%dbh_repro_threshold)
      deallocate(prt_params%seed_alloc_mature)
      deallocate(prt_params%seed_alloc)
      deallocate(prt_params%repro_alloc_a)
      deallocate(prt_params%repro_alloc_b)
      deallocate(prt_params%organ_param_id)
      deallocate(prt_params%fnrt_prof_mode)
      deallocate(prt_params%fnrt_prof_a)
      deallocate(prt_params%fnrt_prof_b)
      deallocate(prt_params%c2b)
      deallocate(prt_params%wood_density)
      deallocate(prt_params%woody)
      deallocate(prt_params%slamax)
      deallocate(prt_params%slatop)
      deallocate(prt_params%allom_sai_scaler)
      deallocate(prt_params%allom_dbh_maxheight)
      deallocate(prt_params%allom_hmode)
      deallocate(prt_params%allom_lmode)
      deallocate(prt_params%allom_fmode)
      deallocate(prt_params%allom_amode)
      deallocate(prt_params%allom_cmode)
      deallocate(prt_params%allom_smode)
      deallocate(prt_params%allom_stmode)
      deallocate(prt_params%allom_dmode)
      deallocate(prt_params%allom_la_per_sa_int)
      deallocate(prt_params%allom_la_per_sa_slp)
      deallocate(prt_params%allom_l2fr)
      deallocate(prt_params%allom_agb_frac)
      deallocate(prt_params%allom_d2h1)
      deallocate(prt_params%allom_d2h2)
      deallocate(prt_params%allom_d2h3)
      deallocate(prt_params%allom_d2bl1)
      deallocate(prt_params%allom_d2bl2)
      deallocate(prt_params%allom_d2bl3)
      deallocate(prt_params%allom_blca_expnt_diff)
      deallocate(prt_params%allom_d2ca_coefficient_max)
      deallocate(prt_params%allom_d2ca_coefficient_min)
      deallocate(prt_params%allom_agb1)
      deallocate(prt_params%allom_agb2)
      deallocate(prt_params%allom_agb3)
      deallocate(prt_params%allom_agb4)
      deallocate(prt_params%allom_h2cd1)
      deallocate(prt_params%allom_h2cd2)
      deallocate(prt_params%allom_zroot_max_dbh)
      deallocate(prt_params%allom_zroot_max_z)
      deallocate(prt_params%allom_zroot_min_dbh)
      deallocate(prt_params%allom_zroot_min_z)
      deallocate(prt_params%allom_zroot_k)
      deallocate(prt_params%pid_kp)
      deallocate(prt_params%pid_ki)
      deallocate(prt_params%pid_kd)
      deallocate(prt_params%store_ovrflw_frac)
      deallocate(prt_params%nfix_mresp_scfrac)
    end if
  end subroutine cleanup_mock_prt_params

end module FatesUnitTestPRTParametersMod
