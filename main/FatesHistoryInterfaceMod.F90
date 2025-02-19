module FatesHistoryInterfaceMod

  use FatesConstantsMod        , only : r8 => fates_r8
  use FatesConstantsMod        , only : fates_avg_flag_length
  use FatesConstantsMod        , only : fates_short_string_length
  use FatesConstantsMod        , only : fates_long_string_length
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : calloc_abs_error
  use FatesConstantsMod        , only : mg_per_kg
  use FatesConstantsMod        , only : pi_const
  use FatesConstantsMod        , only : nearzero
  use FatesConstantsMod        , only : t_water_freeze_k_1atm
  use FatesConstantsMod        , only : n_term_mort_types
  use FatesConstantsMod        , only : i_term_mort_type_cstarv
  use FatesConstantsMod        , only : i_term_mort_type_canlev
  use FatesConstantsMod        , only : i_term_mort_type_numdens
  use FatesConstantsMod        , only : nocomp_bareground_land
  use FatesConstantsMod        , only : nocomp_bareground
  use FatesGlobals             , only : fates_log
  use FatesGlobals             , only : endrun => fates_endrun
  use EDParamsMod              , only : nclmax, maxpft
  use FatesConstantsMod        , only : ican_upper
  use PRTGenericMod            , only : element_pos
  use PRTGenericMod            , only : num_elements
  use PRTGenericMod            , only : prt_cnp_flex_allom_hyp
  use EDTypesMod               , only : site_fluxdiags_type
  use EDTypesMod               , only : elem_diag_type
  use EDtypesMod               , only : ed_site_type
  use FatesCohortMod           , only : fates_cohort_type
  use FatesPatchMod            , only : fates_patch_type
  use EDtypesMod               , only : AREA
  use EDtypesMod               , only : AREA_INV
  use EDTypesMod               , only : numWaterMem
  use EDTypesMod               , only : num_vegtemp_mem
  use PRTGenericMod            , only : element_list
  use FatesIOVariableKindMod   , only : group_dyna_simple, group_dyna_complx
  use FatesIOVariableKindMod   , only : group_hifr_simple, group_hifr_complx
  use FatesIOVariableKindMod   , only : group_hydr_simple, group_hydr_complx
  use FatesIOVariableKindMod   , only : group_nflx_simple, group_nflx_complx
  use FatesConstantsMod        , only : N_DIST_TYPES
  use FatesConstantsMod        , only : dtype_ifall
  use FatesConstantsMod        , only : dtype_ifire
  use FatesConstantsMod        , only : dtype_ilog
  use FatesIODimensionsMod     , only : fates_io_dimension_type
  use FatesIOVariableKindMod   , only : fates_io_variable_kind_type
  use FatesIOVariableKindMod   , only : site_int
  use FatesHistoryVariableType , only : fates_history_variable_type
  use FatesInterfaceTypesMod        , only : hlm_hio_ignore_val
  use FatesInterfaceTypesMod        , only : hlm_use_planthydro
  use FatesInterfaceTypesMod        , only : hlm_use_ed_st3
  use FatesInterfaceTypesMod        , only : hlm_use_cohort_age_tracking
  use FatesInterfaceTypesMod        , only : hlm_use_tree_damage
  use FatesInterfaceTypesMod        , only : nlevdamage
  use FatesInterfaceTypesMod        , only : numpft
  use FatesInterfaceTypesMod        , only : hlm_freq_day
  use FatesInterfaceTypesMod        , only : hlm_parteh_mode
  use FatesInterfaceTypesMod        , only : hlm_use_sp
  use EDParamsMod              , only : ED_val_comp_excln
  use EDParamsMod              , only : ED_val_phen_coldtemp
  use EDParamsMod                   , only : nlevleaf
  use EDParamsMod               , only : ED_val_history_height_bin_edges
  use EDParamsMod               , only : ED_val_history_ageclass_bin_edges
  use FatesInterfaceTypesMod        , only : nlevsclass, nlevage
  use FatesInterfaceTypesMod        , only : nlevheight
  use FatesInterfaceTypesMod        , only : bc_in_type
  use FatesInterfaceTypesMod        , only : bc_out_type
  use FatesInterfaceTypesMod        , only : hlm_model_day
  use FatesInterfaceTypesMod        , only : nlevcoage
  use FatesInterfaceTypesMod        , only : hlm_use_nocomp
  use FatesInterfaceTypesMod        , only : hlm_use_fixed_biogeog
  use FatesRadiationMemMod          , only : ivis,inir
  use FatesInterfaceTypesMod        , only : hlm_hist_level_hifrq,hlm_hist_level_dynam
  use FatesIOVariableKindMod, only : site_r8, site_soil_r8, site_size_pft_r8
  use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
  use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
  use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
  use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
  use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
  use FatesIOVariableKindMod, only : site_height_r8, site_agefuel_r8
  use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
  use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8, site_clscpf_r8
  use FatesIOVariableKindMod, only : site_cdpf_r8, site_cdsc_r8, site_cdam_r8
  use FatesIOVariableKindMod, only : site_landuse_r8, site_lulu_r8, site_lupft_r8
  use FatesConstantsMod   , only : n_landuse_cats
  use FatesAllometryMod             , only : CrownDepth
  use FatesAllometryMod             , only : bstore_allom, bsap_allom
  use FatesAllometryMod             , only : set_root_fraction

  use EDPftvarcon              , only : EDPftvarcon_inst
  use PRTParametersMod         , only : prt_params

  ! CIME Globals
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use shr_infnan_mod           , only : isnan => shr_infnan_isnan

  use FatesConstantsMod        , only : g_per_kg
  use FatesConstantsMod        , only : kg_per_g
  use FatesConstantsMod        , only : ha_per_m2
  use FatesConstantsMod        , only : days_per_sec
  use FatesConstantsMod        , only : sec_per_day
  use FatesConstantsMod        , only : days_per_sec
  use FatesConstantsMod        , only : days_per_year
  use FatesConstantsMod        , only : years_per_day
  use FatesConstantsMod        , only : m2_per_km2
  use FatesConstantsMod        , only : J_per_kJ
  use FatesConstantsMod        , only : m2_per_ha
  use FatesConstantsMod        , only : ha_per_m2
  use FatesConstantsMod        , only : m_per_cm
  use FatesConstantsMod        , only : m_per_mm
  use FatesConstantsMod        , only : sec_per_min
  use FatesConstantsMod        , only : umol_per_mol,mol_per_umol
  use FatesConstantsMod        , only : pa_per_mpa
  use FatesConstantsMod        , only : dens_fresh_liquid_water
  use FatesConstantsMod        , only : grav_earth
  use FatesLitterMod           , only : litter_type
  use FatesConstantsMod        , only : secondaryland
  use FatesConstantsMod        , only : primaryland

  use PRTGenericMod            , only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod            , only : struct_organ, store_organ, repro_organ
  use PRTGenericMod            , only : carbon12_element
  use PRTGenericMod            , only : nitrogen_element, phosphorus_element
  use PRTGenericMod            , only : prt_carbon_allom_hyp
  use PRTAllometricCNPMod      , only : stoich_max,stoich_growth_min
  use FatesSizeAgeTypeIndicesMod, only : get_layersizetype_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_age_class_index

  use FatesFuelClassesMod , only : num_fuel_classes
  use FatesLitterMod      , only : ncwd
  use FatesConstantsMod   , only : ican_upper
  use FatesConstantsMod   , only : ican_ustory
  use FatesSizeAgeTypeIndicesMod, only : get_sizeage_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_sizeagepft_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_agepft_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_agefuel_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_height_index
  use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_cdamagesize_class_index
  use FatesSizeAgeTypeIndicesMod, only : get_cdamagesizepft_class_index
  use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index

  implicit none
  private          ! By default everything is private

  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.  Definitions are not provided, for an explanation of
  ! the variable go to its registry.  (IH_ signifies "index history")
  !
  ! Because of the complex sub-gridscale structure of FATES, in which multiple patches and cohorts
  ! exist within a gridcell, along with vertical gradients within and between canopy layers, as well
  ! as distinct classes such as PFTs or fuel size bins, there are multiple different dimensions in
  ! which it is possible to output history variables to better understand what's going on.
  !
  ! a key point is that, while the number of patches or cohorts can in principle be large, and
  ! the age and size indices of a given patch or cohort can be finely resolved, we collapse these
  ! continuously varying indices into bins of time-invariant width for the purposes of history
  ! outputting.  This is because a given patch or cohort may not persist across a given interval
  ! of history averaging, so it is better to output all patches of cohorts whose index is within
  ! a given interval along the size or age bin.
  !
  ! Another particularity of the issue of FATES shifting its subgrid structure frequently
  ! and possibly having multiple (or zero) patches or cohorts within a given bin is that, if you
  ! want to output an average quantities across some dimension, such as a mean carbon flux across
  ! patch area of a given age, in general it is better to output both the numerator and denominator
  ! of the averaging calculation separately, rather than the average itself, and then calculate
  ! the average in post-processing. So, e.g. this means outputting both the patch area and the
  ! product of the flux within each patch and the patch area as separate variables.  Doing this
  ! allows conservation even when the weights are changing rapidly and simplifies the logic when
  ! the number of patches or cohorts may be anywhere from zero to a large number.
  !
  ! So what this means is that anything that is disaggregated at the patch area requires
  ! outputting the patch age distribution (in units of patch area / site area) as the denominator
  ! of the average and then calculating the numerator of the average as XXX times the patch
  ! area so (so in units of XXX * patch area / site area). For cohort-level quantities,
  ! this requires outputting the number density (in units of individuals per site area), etc.
  !
  ! For reference, some standardized abbreviations of the FATES dimensions are listed here:
  ! scls = size-class dimension
  ! cacls = cohort age-class dimension
  ! pft  = the pft dimension
  ! age  = the age bin dimension
  ! height = the height bin dimension
  ! cwdsc  = the coarse woody debris size class dimension
  !
  ! Since the netcdf interface can only handle variables with a certain number of dimensions,
  ! we have create some "multiplexed" dimensions that combine two or more dimensions into a
  ! single dimension.  Examples of these are the following:
  ! scpf = size class x PFT
  ! cacpf = cohort age class x PFT
  ! cnlf = canopy layer x leaf layer
  ! cnlfpft = canopy layer x leaf layer x PFT
  ! scag = size class bin x age bin
  ! scagpft = size class bin x age bin x PFT
  ! agepft  = age bin x PFT
  ! agefuel = age bin x fuel size class


  ! A recipe for adding a new history variable to this module:
  ! (1) decide what time frequency it makes sense to update the variable at, and what dimension(s)
  !     you want to output the variable on
  ! (2) add the ih_ integer variable in the immediately following section of the module.
  !     use the suffix as outlined above for the dimension you are using.
  ! (3) define a corresponding hio_ variable by associating it to the ih_ variable
  !     in the associate section of the subroutine that corresponds to the time-updating
  !     frequency that you've chosen
  !     (i.e. if half-hourly, then work in subroutine update_history_prod; if daily,
  !     then work in subroutine update_history_dyn)
  ! (4) within that subroutine, add the logic that passes the information from the
  !     fates-native variable (possibly on a patch or cohort structure) to the history
  !     hio_ variable that you've associated to.
  ! (5) add the variable name, metadata, units, dimension, updating frequency, the ih_ variable
  !     index, etc via a call to the set_history_var method in the subroutine define_history_vars.
  !

  ! --STEPS-- TO CONVERT TO HISTORY LEVELS:
  !           SPLIT INTO HIFRQ AND DYNAMICS
  !           GO UP IN ORDER


  ! Indices to 1D Patch variables

  integer :: ih_storec_si
  integer :: ih_storectfrac_si
  integer :: ih_storectfrac_canopy_scpf
  integer :: ih_storectfrac_ustory_scpf
  integer :: ih_leafc_si
  integer :: ih_sapwc_si
  integer :: ih_fnrtc_si
  integer :: ih_fnrtc_sl
  integer :: ih_reproc_si
  integer :: ih_totvegc_si

  ! Nutrient relevant diagnostics (CNP)
  ! ---------------------------------------------------------------
  ! These are active if if(any(element_list(:)==nitrogen_element))
  integer :: ih_storen_si
  integer :: ih_leafn_si
  integer :: ih_sapwn_si
  integer :: ih_fnrtn_si
  integer :: ih_repron_si
  integer :: ih_totvegn_si
  integer :: ih_storentfrac_si
  integer :: ih_totvegn_scpf
  integer :: ih_leafn_scpf
  integer :: ih_fnrtn_scpf
  integer :: ih_storen_scpf
  integer :: ih_sapwn_scpf
  integer :: ih_repron_scpf
  integer :: ih_storentfrac_canopy_scpf
  integer :: ih_storentfrac_understory_scpf

  ! These are active if if(any(element_list(:)==phosphorus_element))
  integer :: ih_storep_si
  integer :: ih_leafp_si
  integer :: ih_sapwp_si
  integer :: ih_fnrtp_si
  integer :: ih_reprop_si
  integer :: ih_totvegp_si
  integer :: ih_storeptfrac_si
  integer :: ih_totvegp_scpf
  integer :: ih_leafp_scpf
  integer :: ih_fnrtp_scpf
  integer :: ih_reprop_scpf
  integer :: ih_storep_scpf
  integer :: ih_sapwp_scpf
  integer :: ih_storeptfrac_canopy_scpf
  integer :: ih_storeptfrac_understory_scpf

  integer :: ih_l2fr_si
  integer :: ih_l2fr_clscpf
  integer :: ih_recl2fr_canopy_pf
  integer :: ih_recl2fr_ustory_pf
  
  integer :: ih_nh4uptake_scpf
  integer :: ih_no3uptake_scpf
  integer :: ih_puptake_scpf
  integer :: ih_nh4uptake_si
  integer :: ih_no3uptake_si
  integer :: ih_puptake_si
  integer :: ih_nefflux_si
  integer :: ih_pefflux_si
  integer :: ih_nefflux_scpf
  integer :: ih_pefflux_scpf
  integer :: ih_nfix_si
  integer :: ih_nfix_scpf
  integer :: ih_ndemand_si
  integer :: ih_ndemand_scpf
  integer :: ih_pdemand_si
  integer :: ih_pdemand_scpf
  
  integer :: ih_trimming_si
  integer :: ih_area_plant_si
  integer :: ih_area_trees_si
  integer :: ih_litter_in_elem
  integer :: ih_litter_out_elem
  integer :: ih_seed_bank_elem
  integer :: ih_fates_fraction_si
  integer :: ih_litter_in_si            ! carbon only
  integer :: ih_litter_out_si           ! carbon only
  integer :: ih_seed_bank_si            ! carbon only
  integer :: ih_seeds_in_si             ! carbon only
  integer :: ih_seeds_in_local_si       ! carbon only
  integer :: ih_ungerm_seed_bank_si        ! carbon only
  integer :: ih_seedling_pool_si    ! carbon only
  integer :: ih_ba_weighted_height_si
  integer :: ih_ca_weighted_height_si
  integer :: ih_seeds_in_local_elem
  integer :: ih_seeds_in_extern_elem
  integer :: ih_seed_decay_elem
  integer :: ih_seed_germ_elem

  integer :: ih_fines_ag_elem
  integer :: ih_fines_bg_elem
  integer :: ih_cwd_ag_elem
  integer :: ih_cwd_bg_elem
  integer :: ih_cwd_elcwd
  integer :: ih_burn_flux_elem

  ! Size-class x PFT mass states

  integer :: ih_bstor_canopy_si_scpf
  integer :: ih_bstor_understory_si_scpf
  integer :: ih_bleaf_canopy_si_scpf
  integer :: ih_bleaf_understory_si_scpf
  !  Size-class x PFT LAI states
  integer :: ih_lai_canopy_si_scpf
  integer :: ih_lai_understory_si_scpf
  !  Size-class x PFT LAI states
  integer :: ih_crownarea_canopy_si_scpf
  integer :: ih_crownarea_understory_si_scpf
  
  integer :: ih_totvegc_scpf
  integer :: ih_leafc_scpf
  integer :: ih_fnrtc_scpf
  integer :: ih_storec_scpf
  integer :: ih_sapwc_scpf
  integer :: ih_reproc_scpf
  integer :: ih_bdead_si
  integer :: ih_balive_si
  integer :: ih_agb_si
  integer :: ih_npp_si
  integer :: ih_gpp_si
  integer :: ih_aresp_si
  integer :: ih_maint_resp_si
  integer :: ih_growth_resp_si
  integer :: ih_excess_resp_si
  integer :: ih_ar_canopy_si
  integer :: ih_gpp_canopy_si
  integer :: ih_ar_understory_si
  integer :: ih_gpp_understory_si
  integer :: ih_canopy_biomass_si
  integer :: ih_understory_biomass_si
  integer :: ih_maint_resp_unreduced_si

  integer :: ih_primaryland_fusion_error_si

  ! land-use-resolved variables
  integer :: ih_area_si_landuse
  integer :: ih_biomass_si_landuse
  integer :: ih_burnedarea_si_landuse
  integer :: ih_gpp_si_landuse
  integer :: ih_npp_si_landuse

  ! land use by land use variables
  integer :: ih_disturbance_rate_si_lulu
  integer :: ih_transition_matrix_si_lulu
  
  integer :: ih_fire_disturbance_rate_si
  integer :: ih_logging_disturbance_rate_si
  integer :: ih_fall_disturbance_rate_si
  integer :: ih_harvest_debt_si
  integer :: ih_harvest_debt_sec_si
  integer :: ih_harvest_woodprod_carbonflux_si
  integer :: ih_luchange_woodprod_carbonflux_si
  
  ! Indices to site by size-class by age variables
  integer :: ih_nplant_si_scag
  integer :: ih_nplant_canopy_si_scag
  integer :: ih_nplant_understory_si_scag
  integer :: ih_ddbh_canopy_si_scag
  integer :: ih_ddbh_understory_si_scag
  integer :: ih_mortality_canopy_si_scag
  integer :: ih_mortality_understory_si_scag

  ! Indices to site by size-class by age by pft variables
  integer :: ih_nplant_si_scagpft

  ! Indices to site by patch age by pft variables
  integer :: ih_biomass_si_agepft
  integer :: ih_npp_si_agepft
  integer :: ih_scorch_height_si_agepft

  ! Indices to (site) variables
  integer :: ih_tveg24_si
  integer :: ih_tlongterm_si
  integer :: ih_tgrowth_si
  integer :: ih_tveg_si
  integer :: ih_nep_si
  integer :: ih_hr_si

  integer :: ih_c_stomata_si
  integer :: ih_c_lblayer_si
  integer :: ih_vis_rad_err_si
  integer :: ih_nir_rad_err_si
  integer :: ih_fire_c_to_atm_si
  integer :: ih_interr_liveveg_elem
  integer :: ih_interr_litter_elem
  integer :: ih_cbal_err_fates_si
  integer :: ih_err_fates_elem

  integer :: ih_npatches_si
  integer :: ih_ncohorts_si
  integer :: ih_demotion_carbonflux_si
  integer :: ih_promotion_carbonflux_si
  integer :: ih_canopy_mortality_carbonflux_si
  integer :: ih_understory_mortality_carbonflux_si
  integer :: ih_canopy_mortality_crownarea_si
  integer :: ih_understory_mortality_crownarea_si
  integer :: ih_canopy_spread_si
  integer :: ih_npp_leaf_si
  integer :: ih_npp_seed_si
  integer :: ih_npp_stem_si
  integer :: ih_npp_froot_si
  integer :: ih_npp_croot_si
  integer :: ih_npp_stor_si
  integer :: ih_leaf_mr_si
  integer :: ih_froot_mr_si
  integer :: ih_livestem_mr_si
  integer :: ih_livecroot_mr_si
  integer :: ih_woodproduct_si
  integer :: ih_h2oveg_si
  integer :: ih_h2oveg_dead_si
  integer :: ih_h2oveg_recruit_si
  integer :: ih_h2oveg_growturn_err_si
  integer :: ih_h2oveg_hydro_err_si
  integer :: ih_lai_si
  integer :: ih_elai_si
  
  integer :: ih_site_cstatus_si
  integer :: ih_gdd_si
  integer :: ih_site_nchilldays_si
  integer :: ih_site_ncolddays_si
  integer :: ih_cleafoff_si
  integer :: ih_cleafon_si

  integer :: ih_nesterov_fire_danger_si
  integer :: ih_fire_nignitions_si
  integer :: ih_fire_fdi_si
  integer :: ih_fire_intensity_area_product_si
  integer :: ih_spitfire_ros_si
  integer :: ih_effect_wspeed_si
  integer :: ih_tfc_ros_si
  integer :: ih_fire_intensity_si
  integer :: ih_fire_area_si
  integer :: ih_fire_fuel_bulkd_si
  integer :: ih_fire_fuel_eff_moist_si
  integer :: ih_fire_fuel_sav_si
  integer :: ih_fire_fuel_mef_si
  integer :: ih_sum_fuel_si
  integer :: ih_fragmentation_scaler_sl

  integer :: ih_nplant_si_scpf
  integer :: ih_gpp_si_scpf
  integer :: ih_npp_totl_si_scpf
  integer :: ih_npp_leaf_si_scpf
  integer :: ih_npp_seed_si_scpf
  integer :: ih_npp_fnrt_si_scpf
  integer :: ih_npp_bgsw_si_scpf
  integer :: ih_npp_bgdw_si_scpf
  integer :: ih_npp_agsw_si_scpf
  integer :: ih_npp_agdw_si_scpf
  integer :: ih_npp_stor_si_scpf
  integer :: ih_grazing_si

  integer :: ih_mortality_canopy_si_scpf
  integer :: ih_mortality_understory_si_scpf
  integer :: ih_m3_mortality_canopy_si_scpf
  integer :: ih_m3_mortality_understory_si_scpf
  integer :: ih_nplant_canopy_si_scpf
  integer :: ih_nplant_understory_si_scpf
  integer :: ih_ddbh_canopy_si_scpf
  integer :: ih_ddbh_understory_si_scpf
  integer :: ih_gpp_canopy_si_scpf
  integer :: ih_gpp_understory_si_scpf
  integer :: ih_ar_canopy_si_scpf
  integer :: ih_ar_understory_si_scpf

  integer :: ih_ddbh_si_scpf
  integer :: ih_growthflux_si_scpf
  integer :: ih_growthflux_fusion_si_scpf
  integer :: ih_ba_si_scpf
  integer :: ih_agb_si_scpf
  integer :: ih_m1_si_scpf
  integer :: ih_m2_si_scpf
  integer :: ih_m3_si_scpf
  integer :: ih_m4_si_scpf
  integer :: ih_m5_si_scpf
  integer :: ih_m6_si_scpf
  integer :: ih_m7_si_scpf
  integer :: ih_m8_si_scpf
  integer :: ih_m9_si_scpf
  integer :: ih_m10_si_scpf
  integer :: ih_m11_si_scpf

  integer :: ih_crownfiremort_si_scpf
  integer :: ih_cambialfiremort_si_scpf

  integer :: ih_abg_mortality_cflux_si_scpf
  integer :: ih_abg_productivity_cflux_si_scpf

  integer :: ih_m10_si_capf
  integer :: ih_nplant_si_capf

  integer :: ih_ar_si_scpf
  integer :: ih_ar_grow_si_scpf
  integer :: ih_ar_maint_si_scpf
  integer :: ih_ar_darkm_si_scpf
  integer :: ih_ar_agsapm_si_scpf
  integer :: ih_ar_crootm_si_scpf
  integer :: ih_ar_frootm_si_scpf

  integer :: ih_c13disc_si_scpf

  ! indices to (site x scls [size class bins]) variables
  integer :: ih_ba_si_scls
  integer :: ih_nplant_si_scls
  integer :: ih_nplant_canopy_si_scls
  integer :: ih_nplant_understory_si_scls
  integer :: ih_lai_canopy_si_scls
  integer :: ih_lai_understory_si_scls
  integer :: ih_sai_canopy_si_scls
  integer :: ih_sai_understory_si_scls
  integer :: ih_mortality_canopy_si_scls
  integer :: ih_mortality_understory_si_scls
  integer :: ih_m3_mortality_canopy_si_scls
  integer :: ih_m3_mortality_understory_si_scls

  integer :: ih_demotion_rate_si_scls
  integer :: ih_promotion_rate_si_scls
  integer :: ih_trimming_canopy_si_scls
  integer :: ih_trimming_understory_si_scls
  integer :: ih_crown_area_canopy_si_scls
  integer :: ih_crown_area_understory_si_scls
  integer :: ih_ddbh_canopy_si_scls
  integer :: ih_ddbh_understory_si_scls
  integer :: ih_agb_si_scls
  integer :: ih_biomass_si_scls
  integer :: ih_mortality_canopy_secondary_si_scls

  ! mortality vars
  integer :: ih_m1_si_scls
  integer :: ih_m2_si_scls
  integer :: ih_m3_si_scls
  integer :: ih_m4_si_scls
  integer :: ih_m5_si_scls
  integer :: ih_m6_si_scls
  integer :: ih_m7_si_scls
  integer :: ih_m8_si_scls
  integer :: ih_m9_si_scls
  integer :: ih_m10_si_scls

  integer :: ih_m10_si_cacls
  integer :: ih_nplant_si_cacls

  ! lots of non-default diagnostics for understanding canopy versus understory carbon balances
  integer :: ih_rdark_canopy_si_scls
  integer :: ih_livestem_mr_canopy_si_scls
  integer :: ih_livecroot_mr_canopy_si_scls
  integer :: ih_froot_mr_canopy_si_scls
  integer :: ih_resp_g_canopy_si_scls
  integer :: ih_resp_m_canopy_si_scls
  integer :: ih_leaf_md_canopy_si_scls
  integer :: ih_root_md_canopy_si_scls
  integer :: ih_carbon_balance_canopy_si_scls
  integer :: ih_bstore_md_canopy_si_scls
  integer :: ih_bdead_md_canopy_si_scls
  integer :: ih_bsw_md_canopy_si_scls
  integer :: ih_seed_prod_canopy_si_scls
  integer :: ih_npp_leaf_canopy_si_scls
  integer :: ih_npp_fnrt_canopy_si_scls
  integer :: ih_npp_sapw_canopy_si_scls
  integer :: ih_npp_dead_canopy_si_scls
  integer :: ih_npp_seed_canopy_si_scls
  integer :: ih_npp_stor_canopy_si_scls

  integer :: ih_rdark_understory_si_scls
  integer :: ih_livestem_mr_understory_si_scls
  integer :: ih_livecroot_mr_understory_si_scls
  integer :: ih_froot_mr_understory_si_scls
  integer :: ih_resp_g_understory_si_scls
  integer :: ih_resp_m_understory_si_scls
  integer :: ih_leaf_md_understory_si_scls
  integer :: ih_root_md_understory_si_scls
  integer :: ih_carbon_balance_understory_si_scls
  integer :: ih_bsw_md_understory_si_scls
  integer :: ih_bdead_md_understory_si_scls
  integer :: ih_bstore_md_understory_si_scls
  integer :: ih_seed_prod_understory_si_scls
  integer :: ih_npp_leaf_understory_si_scls
  integer :: ih_npp_fnrt_understory_si_scls
  integer :: ih_npp_sapw_understory_si_scls
  integer :: ih_npp_dead_understory_si_scls
  integer :: ih_npp_seed_understory_si_scls
  integer :: ih_npp_stor_understory_si_scls

  integer :: ih_yesterdaycanopylevel_canopy_si_scls
  integer :: ih_yesterdaycanopylevel_understory_si_scls

  ! indices to (site x pft) variables
  integer :: ih_biomass_si_pft
  integer :: ih_leafbiomass_si_pft
  integer :: ih_storebiomass_si_pft
  integer :: ih_nindivs_si_pft
  integer :: ih_recruitment_si_pft
  integer :: ih_recruitment_cflux_si_pft
  integer :: ih_mortality_si_pft
  integer :: ih_mortality_carbonflux_si_pft
  integer :: ih_hydraulicmortality_carbonflux_si_pft
  integer :: ih_cstarvmortality_carbonflux_si_pft
  integer :: ih_firemortality_carbonflux_si_pft
  integer :: ih_cstarvmortality_continuous_carbonflux_si_pft
  integer :: ih_crownarea_si_pft
  integer :: ih_canopycrownarea_si_pft
  integer :: ih_crownarea_si_cnlf
  integer :: ih_gpp_si_pft
  integer :: ih_npp_si_pft
  integer :: ih_site_dstatus_si_pft
  integer :: ih_dleafoff_si_pft
  integer :: ih_dleafon_si_pft
  integer :: ih_meanliqvol_si_pft
  integer :: ih_meansmp_si_pft
  integer :: ih_elong_factor_si_pft
  integer :: ih_nocomp_pftpatchfraction_si_pft
  integer :: ih_nocomp_pftnpatches_si_pft
  integer :: ih_nocomp_pftburnedarea_si_pft
  integer :: ih_seeds_out_gc_si_pft
  integer :: ih_seeds_in_gc_si_pft

  ! indices to (site x patch-age) variables
  integer :: ih_area_si_age
  integer :: ih_lai_si_age
  integer :: ih_canopy_area_si_age
  integer :: ih_gpp_si_age
  integer :: ih_npp_si_age
  integer :: ih_ncl_si_age
  integer :: ih_npatches_si_age
  integer :: ih_zstar_si_age
  integer :: ih_biomass_si_age
  integer :: ih_c_stomata_si_age
  integer :: ih_c_lblayer_si_age
  integer :: ih_agesince_anthrodist_si_age
  integer :: ih_secondarylands_area_si_age
  integer :: ih_primarylands_area_si_age
  integer :: ih_area_burnt_si_age
  ! integer :: ih_fire_rate_of_spread_front_si_age
  integer :: ih_fire_intensity_si_age
  integer :: ih_fire_sum_fuel_si_age

  ! indices to (site x height) variables
  integer :: ih_canopy_height_dist_si_height
  integer :: ih_leaf_height_dist_si_height

  ! Indices to hydraulics variables

  integer :: ih_errh2o_scpf
  integer :: ih_tran_scpf

  !  integer :: ih_h2osoi_si_scagpft  ! hijacking the scagpft dimension instead of creating a new shsl dimension
  integer :: ih_sapflow_scpf
  integer :: ih_sapwood_area_scpf
  integer :: ih_sapflow_si
  integer :: ih_iterh1_scpf
  integer :: ih_iterh2_scpf
  integer :: ih_supsub_scpf
  integer :: ih_ath_scpf
  integer :: ih_tth_scpf
  integer :: ih_sth_scpf
  integer :: ih_lth_scpf
  integer :: ih_awp_scpf
  integer :: ih_twp_scpf
  integer :: ih_swp_scpf
  integer :: ih_lwp_scpf
  integer :: ih_aflc_scpf
  integer :: ih_tflc_scpf
  integer :: ih_sflc_scpf
  integer :: ih_lflc_scpf
  integer :: ih_btran_scpf

  ! Hydro: Soil water states
  integer :: ih_rootwgt_soilvwc_si
  integer :: ih_rootwgt_soilvwcsat_si
  integer :: ih_rootwgt_soilmatpot_si

  ! Hydro: Soil water state by layer
  integer :: ih_soilmatpot_sl
  integer :: ih_soilvwc_sl
  integer :: ih_soilvwcsat_sl

  ! Hydro: Root water Uptake rates
  integer :: ih_rootuptake_si
  integer :: ih_rootuptake_sl
  integer :: ih_rootuptake0_scpf
  integer :: ih_rootuptake10_scpf
  integer :: ih_rootuptake50_scpf
  integer :: ih_rootuptake100_scpf


  ! indices to (site x fuel class) variables
  integer :: ih_litter_moisture_si_fuel
  integer :: ih_burnt_frac_litter_si_fuel
  integer :: ih_fuel_amount_si_fuel

  ! indices to (site x cwd size class) variables
  integer :: ih_cwd_ag_si_cwdsc
  integer :: ih_cwd_bg_si_cwdsc
  integer :: ih_cwd_ag_in_si_cwdsc
  integer :: ih_cwd_bg_in_si_cwdsc
  integer :: ih_cwd_ag_out_si_cwdsc
  integer :: ih_cwd_bg_out_si_cwdsc

  ! indices to (site x [canopy layer x leaf layer]) variables
  integer :: ih_parsun_z_si_cnlf
  integer :: ih_parsha_z_si_cnlf
  integer :: ih_laisun_z_si_cnlf
  integer :: ih_laisha_z_si_cnlf
  integer :: ih_ts_net_uptake_si_cnlf
  integer :: ih_crownarea_clll
  integer :: ih_parprof_dir_si_cnlf
  integer :: ih_parprof_dif_si_cnlf

  ! indices to (site x [canopy layer x leaf layer x pft]) variables
  integer :: ih_parsun_z_si_cnlfpft
  integer :: ih_parsha_z_si_cnlfpft
  integer :: ih_laisun_clllpf
  integer :: ih_laisha_clllpf
  integer :: ih_parprof_dir_si_cnlfpft
  integer :: ih_parprof_dif_si_cnlfpft
  integer :: ih_crownfrac_clllpf

  ! indices to site x crown damage variables
  ! site x crown damage x pft x sizeclass
  ! site x crown damage x size class
  integer :: ih_nplant_si_cdpf
  integer :: ih_nplant_canopy_si_cdpf
  integer :: ih_nplant_understory_si_cdpf
  integer :: ih_mortality_si_cdpf
  integer :: ih_mortality_canopy_si_cdpf
  integer :: ih_mortality_understory_si_cdpf
  integer :: ih_m3_si_cdpf
  integer :: ih_m11_si_cdpf
  integer :: ih_m3_mortality_canopy_si_cdpf
  integer :: ih_m3_mortality_understory_si_cdpf
  integer :: ih_m11_mortality_canopy_si_cdpf
  integer :: ih_m11_mortality_understory_si_cdpf
  integer :: ih_ddbh_si_cdpf
  integer :: ih_ddbh_canopy_si_cdpf
  integer :: ih_ddbh_understory_si_cdpf

  ! crownarea damaged
  integer :: ih_crownarea_canopy_damage_si
  integer :: ih_crownarea_ustory_damage_si

  ! indices to (site x canopy layer) variables
  integer :: ih_parsun_si_can
  integer :: ih_parsha_si_can
  integer :: ih_laisun_si_can
  integer :: ih_laisha_si_can
  integer :: ih_crownarea_cl

  ! indices to (patch age x fuel size class) variables
  integer :: ih_fuel_amount_age_fuel

  ! The number of variable dim/kind types we have defined (static)

  integer, parameter, public :: fates_history_num_dimensions = 50
  integer, parameter, public :: fates_history_num_dim_kinds = 50

  type, public :: fates_history_interface_type

     ! Instance of the list of history output varialbes
     type(fates_history_variable_type), allocatable :: hvars(:)
     integer, private :: num_history_vars_

     ! Instanteat one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(fates_io_variable_kind_type) :: dim_kinds(fates_history_num_dim_kinds)

     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure is
     ! allocated by number of threads. This could be dynamically
     ! allocated, but is unlikely to change...?
     type(fates_io_dimension_type) :: dim_bounds(fates_history_num_dimensions)

     !! THESE WERE EXPLICITLY PRIVATE WHEN TYPE WAS PUBLIC
     integer, private :: column_index_, levsoil_index_, levscpf_index_
     integer, private :: levscls_index_, levpft_index_, levage_index_
     integer, private :: levfuel_index_, levcwdsc_index_, levscag_index_
     integer, private :: levcan_index_, levcnlf_index_, levcnlfpft_index_
     integer, private :: levcdpf_index_, levcdsc_index_, levcdam_index_ 
     integer, private :: levscagpft_index_, levagepft_index_
     integer, private :: levheight_index_, levagefuel_index_
     integer, private :: levelem_index_, levelpft_index_
     integer, private :: levelcwd_index_, levelage_index_
     integer, private :: levcacls_index_, levcapf_index_
     integer, private :: levclscpf_index_
     integer, private :: levlanduse_index_, levlulu_index_, levlupft_index_

   contains

     procedure :: Init
     procedure :: SetThreadBoundsEach
     procedure :: initialize_history_vars
     procedure :: assemble_history_output_types

     procedure :: update_history_dyn
     procedure :: update_history_dyn1
     procedure :: update_history_dyn2
     procedure :: update_history_hifrq
     procedure :: update_history_hifrq1
     procedure :: update_history_hifrq2
     procedure :: update_history_hydraulics
     procedure :: update_history_nutrflux

     ! 'get' methods used by external callers to access private read only data

     procedure :: num_history_vars
     procedure :: column_index
     procedure :: levsoil_index
     procedure :: levscpf_index
     procedure :: levscls_index
     procedure :: levcapf_index
     procedure :: levcacls_index
     procedure :: levpft_index
     procedure :: levage_index
     procedure :: levfuel_index
     procedure :: levcwdsc_index
     procedure :: levcan_index
     procedure :: levcnlf_index
     procedure :: levcnlfpft_index
     procedure :: levcdpf_index
     procedure :: levcdsc_index
     procedure :: levcdam_index
     procedure :: levscag_index
     procedure :: levscagpft_index
     procedure :: levagepft_index
     procedure :: levheight_index
     procedure :: levelem_index
     procedure :: levelpft_index
     procedure :: levelcwd_index
     procedure :: levelage_index
     procedure :: levagefuel_index
     procedure :: levclscpf_index
     procedure :: levlanduse_index
     procedure :: levlulu_index
     procedure :: levlupft_index

     ! private work functions
     procedure, private :: define_history_vars
     procedure, private :: set_history_var
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: set_column_index
     procedure, private :: set_levsoil_index
     procedure, private :: set_levscpf_index
     procedure, private :: set_levcacls_index
     procedure, private :: set_levcapf_index
     procedure, private :: set_levscls_index
     procedure, private :: set_levpft_index
     procedure, private :: set_levage_index
     procedure, private :: set_levfuel_index
     procedure, private :: set_levcwdsc_index
     procedure, private :: set_levcan_index
     procedure, private :: set_levcnlf_index
     procedure, private :: set_levcnlfpft_index
     procedure, private :: set_levcdpf_index
     procedure, private :: set_levcdsc_index
     procedure, private :: set_levcdam_index
     procedure, private :: set_levscag_index
     procedure, private :: set_levscagpft_index
     procedure, private :: set_levagepft_index
     procedure, private :: set_levheight_index
     procedure, private :: set_levagefuel_index
     procedure, private :: set_levclscpf_index
     procedure, private :: set_levlanduse_index
     procedure, private :: set_levlulu_index
     procedure, private :: set_levlupft_index
     procedure, private :: set_levelem_index
     procedure, private :: set_levelpft_index
     procedure, private :: set_levelcwd_index
     procedure, private :: set_levelage_index

     procedure, public :: flush_hvars
     procedure, public :: zero_site_hvars
     procedure, public :: flush_all_hvars

  end type fates_history_interface_type

  character(len=*), parameter :: sourcefile = &
       __FILE__


  ! The instance of the type

  type(fates_history_interface_type), public :: fates_hist


contains

  ! ======================================================================

  subroutine Init(this, num_threads, fates_bounds)

    use FatesIODimensionsMod, only : column, levsoil, levscpf
    use FatesIODimensionsMod, only : levscls, levpft, levage
    use FatesIODimensionsMod, only : levcacls, levcapf
    use FatesIODimensionsMod, only : levfuel, levcwdsc, levscag
    use FatesIODimensionsMod, only : levscagpft, levagepft
    use FatesIODimensionsMod, only : levcan, levcnlf, levcnlfpft
    use FatesIODimensionsMod, only : fates_bounds_type
    use FatesIODimensionsMod, only : levheight, levagefuel
    use FatesIODimensionsMod, only : levelem, levelpft
    use FatesIODimensionsMod, only : levelcwd, levelage, levclscpf
    use FatesIODimensionsMod, only : levcdpf, levcdsc, levcdam
    use FatesIODimensionsMod, only : levlanduse, levlulu, levlupft

    implicit none

    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    dim_count = dim_count + 1
    call this%set_levsoil_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levsoil, num_threads, &
         fates_bounds%soil_begin, fates_bounds%soil_end)

    dim_count = dim_count + 1
    call this%set_levscpf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscpf, num_threads, &
         fates_bounds%sizepft_class_begin, fates_bounds%sizepft_class_end)

    dim_count = dim_count + 1
    call this%set_levscls_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscls, num_threads, &
         fates_bounds%size_class_begin, fates_bounds%size_class_end)

    dim_count = dim_count + 1
    call this%set_levcacls_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcacls, num_threads, &
         fates_bounds%coage_class_begin, fates_bounds%coage_class_end)

    dim_count = dim_count + 1
    call this%set_levcapf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcapf, num_threads, &
         fates_bounds%coagepf_class_begin, fates_bounds%coagepf_class_end)

    dim_count = dim_count + 1
    call this%set_levpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levpft, num_threads, &
         fates_bounds%pft_class_begin, fates_bounds%pft_class_end)

    dim_count = dim_count + 1
    call this%set_levage_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levage, num_threads, &
         fates_bounds%age_class_begin, fates_bounds%age_class_end)

    dim_count = dim_count + 1
    call this%set_levfuel_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levfuel, num_threads, &
         fates_bounds%fuel_begin, fates_bounds%fuel_end)

    dim_count = dim_count + 1
    call this%set_levcwdsc_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcwdsc, num_threads, &
         fates_bounds%cwdsc_begin, fates_bounds%cwdsc_end)

    dim_count = dim_count + 1
    call this%set_levcan_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcan, num_threads, &
         fates_bounds%can_begin, fates_bounds%can_end)

    dim_count = dim_count + 1
    call this%set_levcnlf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlf, num_threads, &
         fates_bounds%cnlf_begin, fates_bounds%cnlf_end)

    dim_count = dim_count + 1
    call this%set_levcnlfpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlfpft, num_threads, &
         fates_bounds%cnlfpft_begin, fates_bounds%cnlfpft_end)

    dim_count = dim_count + 1
    call this%set_levcdpf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcdpf, num_threads, &
         fates_bounds%cdpf_begin, fates_bounds%cdpf_end)

    dim_count = dim_count + 1
    call this%set_levcdsc_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcdsc, num_threads, &
         fates_bounds%cdsc_begin, fates_bounds%cdsc_end)

    dim_count = dim_count + 1
    call this%set_levcdam_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcdam, num_threads, &
         fates_bounds%cdam_begin, fates_bounds%cdam_end)

    dim_count = dim_count + 1
    call this%set_levscag_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscag, num_threads, &
         fates_bounds%sizeage_class_begin, fates_bounds%sizeage_class_end)

    dim_count = dim_count + 1
    call this%set_levscagpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscagpft, num_threads, &
         fates_bounds%sizeagepft_class_begin, fates_bounds%sizeagepft_class_end)

    dim_count = dim_count + 1
    call this%set_levagepft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levagepft, num_threads, &
         fates_bounds%agepft_class_begin, fates_bounds%agepft_class_end)

    dim_count = dim_count + 1
    call this%set_levheight_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levheight, num_threads, &
         fates_bounds%height_begin, fates_bounds%height_end)

    dim_count = dim_count + 1
    call this%set_levelem_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelem, num_threads, &
         fates_bounds%elem_begin, fates_bounds%elem_end)

    dim_count = dim_count + 1
    call this%set_levelpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelpft, num_threads, &
         fates_bounds%elpft_begin, fates_bounds%elpft_end)

    dim_count = dim_count + 1
    call this%set_levelcwd_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelcwd, num_threads, &
         fates_bounds%elcwd_begin, fates_bounds%elcwd_end)

    dim_count = dim_count + 1
    call this%set_levelage_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelage, num_threads, &
         fates_bounds%elage_begin, fates_bounds%elage_end)

    dim_count = dim_count + 1
    call this%set_levagefuel_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levagefuel, num_threads, &
         fates_bounds%agefuel_begin, fates_bounds%agefuel_end)

    dim_count = dim_count + 1
    call this%set_levclscpf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levclscpf, num_threads, &
         fates_bounds%clscpf_begin, fates_bounds%clscpf_end)

    dim_count = dim_count + 1
    call this%set_levlanduse_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levlanduse, num_threads, &
         fates_bounds%landuse_begin, fates_bounds%landuse_end)

    dim_count = dim_count + 1
    call this%set_levlulu_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levlulu, num_threads, &
         fates_bounds%lulu_begin, fates_bounds%lulu_end)

    dim_count = dim_count + 1
    call this%set_levlupft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levlupft, num_threads, &
         fates_bounds%lupft_begin, fates_bounds%lupft_end)

  end subroutine Init

  ! ======================================================================
  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)

    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index

    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)

    index = this%levsoil_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%soil_begin, thread_bounds%soil_end)

    index = this%levscpf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%sizepft_class_begin, thread_bounds%sizepft_class_end)

    index = this%levscls_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%size_class_begin, thread_bounds%size_class_end)

    index = this%levcacls_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%coage_class_begin, thread_bounds%coage_class_end)

    index = this%levcapf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%coagepf_class_begin, thread_bounds%coagepf_class_end)

    index = this%levpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%pft_class_begin, thread_bounds%pft_class_end)

    index = this%levage_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%age_class_begin, thread_bounds%age_class_end)

    index = this%levfuel_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%fuel_begin, thread_bounds%fuel_end)

    index = this%levcwdsc_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cwdsc_begin, thread_bounds%cwdsc_end)

    index = this%levcan_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%can_begin, thread_bounds%can_end)

    index = this%levcnlf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cnlf_begin, thread_bounds%cnlf_end)

    index = this%levcnlfpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cnlfpft_begin, thread_bounds%cnlfpft_end)

    index = this%levcdpf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cdpf_begin, thread_bounds%cdpf_end)

    index = this%levcdsc_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cdsc_begin, thread_bounds%cdsc_end)

    index = this%levcdam_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cdam_begin, thread_bounds%cdam_end)

    index = this%levscag_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%sizeage_class_begin, thread_bounds%sizeage_class_end)

    index = this%levscagpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%sizeagepft_class_begin, thread_bounds%sizeagepft_class_end)

    index = this%levagepft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%agepft_class_begin, thread_bounds%agepft_class_end)

    index = this%levheight_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%height_begin, thread_bounds%height_end)

    index = this%levelem_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elem_begin, thread_bounds%elem_end)

    index = this%levelpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elpft_begin, thread_bounds%elpft_end)

    index = this%levelcwd_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elcwd_begin, thread_bounds%elcwd_end)

    index = this%levelage_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elage_begin, thread_bounds%elage_end)

    index = this%levagefuel_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%agefuel_begin, thread_bounds%agefuel_end)

    index = this%levclscpf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%clscpf_begin, thread_bounds%clscpf_end)

    index = this%levlanduse_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%landuse_begin, thread_bounds%landuse_end)

    index = this%levlulu_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%lulu_begin, thread_bounds%lulu_end)

    index = this%levlupft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%lupft_begin, thread_bounds%lupft_end)


  end subroutine SetThreadBoundsEach

  ! ===================================================================================
  subroutine assemble_history_output_types(this)



    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(site_r8, 1, this%column_index())

    call this%set_dim_indices(site_soil_r8, 1, this%column_index())
    call this%set_dim_indices(site_soil_r8, 2, this%levsoil_index())

    call this%set_dim_indices(site_size_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_pft_r8, 2, this%levscpf_index())

    call this%set_dim_indices(site_size_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_r8, 2, this%levscls_index())

    call this%set_dim_indices(site_coage_r8, 1, this%column_index())
    call this%set_dim_indices(site_coage_r8, 2, this%levcacls_index())

    call this%set_dim_indices(site_coage_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_coage_pft_r8, 2, this%levcapf_index())

    call this%set_dim_indices(site_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_pft_r8, 2, this%levpft_index())

    call this%set_dim_indices(site_age_r8, 1, this%column_index())
    call this%set_dim_indices(site_age_r8, 2, this%levage_index())

    call this%set_dim_indices(site_fuel_r8, 1, this%column_index())
    call this%set_dim_indices(site_fuel_r8, 2, this%levfuel_index())

    call this%set_dim_indices(site_cwdsc_r8, 1, this%column_index())
    call this%set_dim_indices(site_cwdsc_r8, 2, this%levcwdsc_index())

    call this%set_dim_indices(site_can_r8, 1, this%column_index())
    call this%set_dim_indices(site_can_r8, 2, this%levcan_index())

    call this%set_dim_indices(site_cnlf_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlf_r8, 2, this%levcnlf_index())

    call this%set_dim_indices(site_cnlfpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlfpft_r8, 2, this%levcnlfpft_index())

    call this%set_dim_indices(site_cdpf_r8, 1, this%column_index())
    call this%set_dim_indices(site_cdpf_r8, 2, this%levcdpf_index())

    call this%set_dim_indices(site_cdsc_r8, 1, this%column_index())
    call this%set_dim_indices(site_cdsc_r8, 2, this%levcdsc_index())

    call this%set_dim_indices(site_cdam_r8, 1, this%column_index())
    call this%set_dim_indices(site_cdam_r8, 2, this%levcdam_index())

    call this%set_dim_indices(site_scag_r8, 1, this%column_index())
    call this%set_dim_indices(site_scag_r8, 2, this%levscag_index())

    call this%set_dim_indices(site_scagpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_scagpft_r8, 2, this%levscagpft_index())

    call this%set_dim_indices(site_agepft_r8, 1, this%column_index())
    call this%set_dim_indices(site_agepft_r8, 2, this%levagepft_index())

    call this%set_dim_indices(site_height_r8, 1, this%column_index())
    call this%set_dim_indices(site_height_r8, 2, this%levheight_index())

    call this%set_dim_indices(site_elem_r8, 1, this%column_index())
    call this%set_dim_indices(site_elem_r8, 2, this%levelem_index())

    call this%set_dim_indices(site_elpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_elpft_r8, 2, this%levelpft_index())

    call this%set_dim_indices(site_elcwd_r8, 1, this%column_index())
    call this%set_dim_indices(site_elcwd_r8, 2, this%levelcwd_index())

    call this%set_dim_indices(site_elage_r8, 1, this%column_index())
    call this%set_dim_indices(site_elage_r8, 2, this%levelage_index())

    call this%set_dim_indices(site_agefuel_r8, 1, this%column_index())
    call this%set_dim_indices(site_agefuel_r8, 2, this%levagefuel_index())

    call this%set_dim_indices(site_clscpf_r8, 1, this%column_index())
    call this%set_dim_indices(site_clscpf_r8, 2, this%levclscpf_index())

    call this%set_dim_indices(site_landuse_r8, 1, this%column_index())
    call this%set_dim_indices(site_landuse_r8, 2, this%levlanduse_index())

    call this%set_dim_indices(site_lulu_r8, 1, this%column_index())
    call this%set_dim_indices(site_lulu_r8, 2, this%levlulu_index())

    call this%set_dim_indices(site_lupft_r8, 1, this%column_index())
    call this%set_dim_indices(site_lupft_r8, 2, this%levlupft_index())

  end subroutine assemble_history_output_types

  ! ===================================================================================

  subroutine set_dim_indices(this, dk_name, idim, dim_index)

    use FatesIOVariableKindMod , only : iotype_index

    implicit none

    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)     :: dk_name
    integer, intent(in)              :: idim  ! dimension index
    integer, intent(in) :: dim_index


    ! local
    integer :: ityp

    ityp = iotype_index(trim(dk_name), fates_history_num_dim_kinds, this%dim_kinds)

    ! First check to see if the dimension is allocated
    if (this%dim_kinds(ityp)%ndims < idim) then
       write(fates_log(), *) 'Trying to define dimension size to a dim-type structure'
       write(fates_log(), *) 'but the dimension index does not exist'
       write(fates_log(), *) 'type: ',dk_name,' ndims: ',this%dim_kinds(ityp)%ndims,' input dim:',idim
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (idim == 1) then
       this%dim_kinds(ityp)%dim1_index = dim_index
    else if (idim == 2) then
       this%dim_kinds(ityp)%dim2_index = dim_index
    end if

    ! With the map, we can set the dimension size
    this%dim_kinds(ityp)%dimsize(idim) = this%dim_bounds(dim_index)%upper_bound - &
         this%dim_bounds(dim_index)%lower_bound + 1

  end subroutine set_dim_indices

  ! =======================================================================
  subroutine set_column_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%column_index_ = index
  end subroutine set_column_index

  integer function column_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    column_index = this%column_index_
  end function column_index

  ! =======================================================================
  subroutine set_levsoil_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levsoil_index_ = index
  end subroutine set_levsoil_index

  integer function levsoil_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levsoil_index = this%levsoil_index_
  end function levsoil_index

  ! =======================================================================
  subroutine set_levscpf_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levscpf_index_ = index
  end subroutine set_levscpf_index

  integer function levscpf_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscpf_index = this%levscpf_index_
  end function levscpf_index

  ! =======================================================================
  subroutine set_levscls_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levscls_index_ = index
  end subroutine set_levscls_index

  integer function levscls_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscls_index = this%levscls_index_
  end function levscls_index

  !=========================================================================
  subroutine set_levcacls_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcacls_index_ = index
  end subroutine set_levcacls_index

  integer function levcacls_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcacls_index = this%levcacls_index_
  end function levcacls_index

  !=========================================================================
  subroutine set_levcapf_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcapf_index_ = index
  end subroutine set_levcapf_index

  integer function levcapf_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcapf_index = this%levcapf_index_
  end function levcapf_index

  ! =======================================================================
  subroutine set_levpft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levpft_index_ = index
  end subroutine set_levpft_index

  integer function levpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levpft_index = this%levpft_index_
  end function levpft_index

  ! =======================================================================
  subroutine set_levage_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levage_index_ = index
  end subroutine set_levage_index

  integer function levage_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levage_index = this%levage_index_
  end function levage_index

  ! =======================================================================
  subroutine set_levfuel_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levfuel_index_ = index
  end subroutine set_levfuel_index

  integer function levfuel_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levfuel_index = this%levfuel_index_
  end function levfuel_index

  ! =======================================================================
  subroutine set_levcwdsc_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcwdsc_index_ = index
  end subroutine set_levcwdsc_index

  integer function levcwdsc_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcwdsc_index = this%levcwdsc_index_
  end function levcwdsc_index

  ! =======================================================================
  subroutine set_levcan_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcan_index_ = index
  end subroutine set_levcan_index

  integer function levcan_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcan_index = this%levcan_index_
  end function levcan_index

  ! =======================================================================
  subroutine set_levcnlf_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcnlf_index_ = index
  end subroutine set_levcnlf_index

  integer function levcnlf_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcnlf_index = this%levcnlf_index_
  end function levcnlf_index

  ! =======================================================================
  subroutine set_levcnlfpft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcnlfpft_index_ = index
  end subroutine set_levcnlfpft_index

  integer function levcnlfpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcnlfpft_index = this%levcnlfpft_index_
  end function levcnlfpft_index

  ! =======================================================================
  subroutine set_levcdpf_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcdpf_index_ = index
  end subroutine set_levcdpf_index

  integer function levcdpf_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcdpf_index = this%levcdpf_index_
  end function levcdpf_index

  ! =======================================================================
  subroutine set_levcdsc_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcdsc_index_ = index
  end subroutine set_levcdsc_index

  integer function levcdsc_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcdsc_index = this%levcdsc_index_
  end function levcdsc_index

  ! =======================================================================
  subroutine set_levcdam_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levcdam_index_ = index
  end subroutine set_levcdam_index

  integer function levcdam_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levcdam_index = this%levcdam_index_
  end function levcdam_index

  ! ======================================================================================
  subroutine set_levscag_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levscag_index_ = index
  end subroutine set_levscag_index

  integer function levscag_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscag_index = this%levscag_index_
  end function levscag_index

  ! ======================================================================================
  subroutine set_levscagpft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levscagpft_index_ = index
  end subroutine set_levscagpft_index

  integer function levscagpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscagpft_index = this%levscagpft_index_
  end function levscagpft_index

  ! ======================================================================================
  subroutine set_levagepft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levagepft_index_ = index
  end subroutine set_levagepft_index

  integer function levagepft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levagepft_index = this%levagepft_index_
  end function levagepft_index

  ! ======================================================================================
  subroutine set_levheight_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levheight_index_ = index
  end subroutine set_levheight_index

  integer function levheight_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levheight_index = this%levheight_index_
  end function levheight_index

  ! ======================================================================================

  subroutine set_levelem_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levelem_index_ = index
  end subroutine set_levelem_index

  integer function levelem_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelem_index = this%levelem_index_
  end function levelem_index

  ! ======================================================================================

  subroutine set_levelpft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levelpft_index_ = index
  end subroutine set_levelpft_index

  integer function levelpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelpft_index = this%levelpft_index_
  end function levelpft_index

  ! ======================================================================================

  subroutine set_levelcwd_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levelcwd_index_ = index
  end subroutine set_levelcwd_index

  integer function levelcwd_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelcwd_index = this%levelcwd_index_
  end function levelcwd_index

  ! ======================================================================================

  subroutine set_levelage_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levelage_index_ = index
  end subroutine set_levelage_index

  integer function levelage_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelage_index = this%levelage_index_
  end function levelage_index

  ! ======================================================================================

  subroutine set_levagefuel_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levagefuel_index_ = index
  end subroutine set_levagefuel_index

  integer function levagefuel_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levagefuel_index = this%levagefuel_index_
  end function levagefuel_index
  ! ======================================================================================

  subroutine set_levclscpf_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levclscpf_index_ = index
  end subroutine set_levclscpf_index

  integer function levclscpf_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levclscpf_index = this%levclscpf_index_
  end function levclscpf_index

  ! ======================================================================================

  subroutine set_levlanduse_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levlanduse_index_ = index
  end subroutine set_levlanduse_index

  integer function levlanduse_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levlanduse_index = this%levlanduse_index_
  end function levlanduse_index

  ! ======================================================================================

  subroutine set_levlulu_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levlulu_index_ = index
  end subroutine set_levlulu_index

  integer function levlulu_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levlulu_index = this%levlulu_index_
  end function levlulu_index

  ! ======================================================================================

  subroutine set_levlupft_index(this, index)
    implicit none
    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%levlupft_index_ = index
  end subroutine set_levlupft_index

  integer function levlupft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levlupft_index = this%levlupft_index_
  end function levlupft_index

  ! =====================================================================================

  subroutine zero_site_hvars(this, currentSite, upfreq_in)

    ! This routine zero's a history diagnostic variable
    ! but only zero's on fates sites
    ! This should be called prior to filling the variable
    ! and after they have been flushed to the ignore value

    class(fates_history_interface_type)    :: this        ! hvars_interface instance
    integer, intent(in)                    :: upfreq_in   !
    type(ed_site_type), intent(in), target :: currentSite ! site instance

    integer :: ivar     ! history variable index
    integer :: ndims    ! number of dimensions

    do ivar=1,ubound(this%hvars,1)
       if (this%hvars(ivar)%upfreq == upfreq_in) then 

          ndims = this%dim_kinds(this%hvars(ivar)%dim_kinds_index)%ndims

          if(trim(this%dim_kinds(this%hvars(ivar)%dim_kinds_index)%name) == site_int)then
             write(fates_log(),*)'add in zeroing provision for SI_INT'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          if(ndims==1) then
             this%hvars(ivar)%r81d(currentSite%h_gid) = 0._r8
          elseif(ndims==2) then
             this%hvars(ivar)%r82d(currentSite%h_gid,:) = 0._r8
          elseif(ndims==3) then
             this%hvars(ivar)%r83d(currentSite%h_gid,:,:) = 0._r8
          end if
       end if
    end do

    return
  end subroutine zero_site_hvars


  ! ======================================================================================

  subroutine flush_all_hvars(this,nc)

    ! A wrapper to flush all active history
    ! groups to their flush value
    
    class(fates_history_interface_type)    :: this
    integer,intent(in)                     :: nc
    
    if(hlm_hist_level_hifrq>0) then
       call this%flush_hvars(nc,upfreq_in=group_hifr_simple)
       if (hlm_use_planthydro.eq.itrue) call this%flush_hvars(nc,upfreq_in=group_hydr_simple)
       if(hlm_hist_level_hifrq>1) then
          call this%flush_hvars(nc,upfreq_in=group_hifr_complx)
          if (hlm_use_planthydro.eq.itrue) call this%flush_hvars(nc,upfreq_in=group_hydr_complx)
       end if
    end if

    if(hlm_hist_level_dynam>0) then
       call this%flush_hvars(nc,upfreq_in=group_dyna_simple)
       call this%flush_hvars(nc,upfreq_in=group_nflx_simple)
       if(hlm_hist_level_dynam>1) then
          call this%flush_hvars(nc,upfreq_in=group_dyna_complx)
          call this%flush_hvars(nc,upfreq_in=group_nflx_complx)
       end if
    end if
    
    return
  end subroutine flush_all_hvars
  
  ! ======================================================================================

  subroutine flush_hvars(this,nc,upfreq_in)

    class(fates_history_interface_type)        :: this
    integer,intent(in)                     :: nc
    integer,intent(in)                     :: upfreq_in
    integer                      :: ivar
    integer                      :: lb1,ub1,lb2,ub2

    do ivar=1,ubound(this%hvars,1)
       if (this%hvars(ivar)%upfreq == upfreq_in) then ! Only flush variables with update on dynamics step
          call this%hvars(ivar)%HFlush(nc, this%dim_bounds, this%dim_kinds)
       end if
    end do

  end subroutine flush_hvars


  ! =====================================================================================

  subroutine set_history_var(this, vname, units, long, use_default, avgflag, vtype, &
       hlms, upfreq, ivar, initialize, index, flush_to_zero)

    use FatesUtilsMod, only     : check_hlm_list
    use FatesInterfaceTypesMod, only : hlm_name

    implicit none

    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)  :: vname
    character(len=*), intent(in)  :: units
    character(len=*), intent(in)  :: long
    character(len=*), intent(in)  :: use_default
    character(len=*), intent(in)  :: avgflag
    character(len=*), intent(in)  :: vtype
    character(len=*), intent(in)  :: hlms
    integer, intent(in)           :: upfreq
    logical, intent(in)           :: initialize
    integer, intent(inout)        :: ivar
    integer, intent(inout)        :: index  ! This is the index for the variable of
    ! interest that is associated with an
    ! explict name (for fast reference during update)
    ! A zero is passed back when the variable is
    ! not used
    logical, intent(in), optional :: flush_to_zero

    ! locals
    integer   :: ub1, lb1, ub2, lb2    ! Bounds for allocating the var
    integer   :: ityp
    real(r8)  :: flushval
    logical   :: write_var


    ! Flushing to the ignore val coerces all FATES diagnostics to be
    ! relevant only on FATES sites. This way we do not average zero's
    ! at locations not on FATES columns
    ! We make one exception to this rule, for the fates_fraction variable.  That way
    ! we can always know what fraction of the gridcell FATES is occupying.

    flushval = hlm_hio_ignore_val
    if (present(flush_to_zero)) then
       if (flush_to_zero) then
          flushval = 0.0_r8
       endif
    endif

    write_var = check_hlm_list(trim(hlms), trim(hlm_name))
    if( write_var ) then
       ivar  = ivar+1
       index = ivar

       if (initialize) then
          call this%hvars(ivar)%Init(vname, units, long, use_default,          &
               vtype, avgflag, flushval, upfreq, fates_history_num_dim_kinds,  &
               this%dim_kinds, this%dim_bounds)
       end if
    else
       index = 0
    end if

    return
  end subroutine set_history_var

  ! ====================================================================================

  subroutine init_dim_kinds_maps(this)

    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : site_r8, site_soil_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_height_r8, site_agefuel_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8, site_clscpf_r8
    use FatesIOVariableKindMod, only : site_cdpf_r8, site_cdsc_r8, site_cdam_r8

    implicit none

    ! Arguments
    class(fates_history_interface_type), intent(inout) :: this


    integer :: index

    index = 1
    ! 1d Site
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! site x soil
    index = index + 1
    call this%dim_kinds(index)%Init(site_soil_r8, 2)

    ! site x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_pft_r8, 2)

    ! site x size-class
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_r8, 2)

    ! site x cohort age-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_coage_pft_r8, 2)

    ! site x cohort age-class
    index = index + 1
    call this%dim_kinds(index)%Init(site_coage_r8, 2)

    ! site x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_pft_r8, 2)

    ! site x patch-age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_age_r8, 2)

    ! site x fuel size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_fuel_r8, 2)

    ! site x cwd size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cwdsc_r8, 2)

    ! site x can class
    index = index + 1
    call this%dim_kinds(index)%Init(site_can_r8, 2)

    ! site x cnlf class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlf_r8, 2)

    ! site x cnlfpft class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlfpft_r8, 2)

    ! site x crown damage x pft x size class 
    index = index + 1
    call this%dim_kinds(index)%Init(site_cdpf_r8, 2)

    ! site x crown damage x size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cdsc_r8, 2)

    ! site x crown damage
    index = index + 1
    call this%dim_kinds(index)%Init(site_cdam_r8, 2)

    ! site x size-class x age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_scag_r8, 2)

    ! site x size-class x age class x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_scagpft_r8, 2)

    ! site x age class x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_agepft_r8, 2)

    ! site x height
    index = index + 1
    call this%dim_kinds(index)%Init(site_height_r8, 2)

    ! site x elemenet
    index = index + 1
    call this%dim_kinds(index)%Init(site_elem_r8, 2)

    ! site x element x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_elpft_r8, 2)

    ! site x element x cwd
    index = index + 1
    call this%dim_kinds(index)%Init(site_elcwd_r8, 2)

    ! site x element x age
    index = index + 1
    call this%dim_kinds(index)%Init(site_elage_r8, 2)

    ! site x age x fuel size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_agefuel_r8, 2)

    ! site x age x fuel size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_clscpf_r8, 2)

    ! site x land use class
    index = index + 1
    call this%dim_kinds(index)%Init(site_landuse_r8, 2)

    ! site x land use x land use class
    index = index + 1
    call this%dim_kinds(index)%Init(site_lulu_r8, 2)

    ! site x land use x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_lupft_r8, 2)

    ! FIXME(bja, 2016-10) assert(index == fates_history_num_dim_kinds)
  end subroutine init_dim_kinds_maps

  ! =======================================================================

  subroutine update_history_nutrflux(this,csite)

    ! TODO IN FUTURE PR:
    !       CHANGE THIS NAME FROM NUTRFLUX TO
    !       DYNAM_FLUX AND THE EXISISTING DYNAMICS
    !       TO DYNAM_STATE. MOVE FLUX DIAGNOSTIC
    !       TO THIS ROUTINE, AND USE THIS ROUTINE
    !       TO TRACK ALL DYNAMICS FLUX DIAGNOSTICS
    !       AND CALL THIS AFTER DISTURBANCE RATES
    !       HAVE BEEN CALLED, BUT BEFORE PATCHES
    !       HAVE BEEN SPLIT AND SPAWNED

    
    ! Update history diagnostics for nutrient dynamics variables.
    ! This is a separate routine because we like to handle these
    ! things before patches are reshuffled during disturbance, and
    ! thus this is called immediately after PARTEH allocation
    ! These diagnostics must be zero'd at the beginning
    ! of the dynamics call (not here, because this is a
    ! being called at the cohort level)

    ! Arguments
    class(fates_history_interface_type) :: this
    type(ed_site_type), intent(in)      :: csite

    type(fates_patch_type), pointer     :: cpatch
    type(fates_cohort_type), pointer    :: ccohort
    integer :: iclscpf   ! layer x size x pft class index
    integer :: iscpf     ! Size x pft class index
    integer :: io_si     ! site's global index in the history vector
    integer :: el        ! element loop index
    integer :: ft        ! pft loop index
    real(r8):: uconv     ! combined unit conversion factor
    real(r8) :: fnrt_c   ! cohort fine-root c

    ! Process variables with time-space dimensions only
    ! ---------------------------------------------------------------------------------------------

    if_dynam1: if(hlm_hist_level_dynam>0) then
       
         ! history site index
         io_si  = csite%h_gid

         ! zero nutrient fluxes
         call this%zero_site_hvars(csite,upfreq_in=group_nflx_simple)
         
         cpatch => csite%youngest_patch
         do while(associated(cpatch))

            ccohort => cpatch%shortest
            do while(associated(ccohort))

               ! If this is a new cohort, do not make diagnostics
               if(ccohort%isnew) then
                  ccohort => ccohort%taller
                  cycle
               end if

               ! unit conversion factor to get x/plant/day -> x/m2/sec
               uconv = ccohort%n * ha_per_m2 * days_per_sec

               fnrt_c   = ccohort%prt%GetState(fnrt_organ, carbon12_element)

               ! Loop over the different elements. 
               do el = 1, num_elements

                  select case (element_list(el))
                  case (carbon12_element)

                     ! Excess carbon respired
                     this%hvars(ih_excess_resp_si)%r81d(io_si) = &
                          this%hvars(ih_excess_resp_si)%r81d(io_si) + &
                          ccohort%resp_excess_hold*uconv/days_per_year

                  case (nitrogen_element) 

                     ! Mineralized uptake of NH4, NO3
                     fates_hist%hvars(ih_nh4uptake_si)%r81d(io_si) =       &
                          fates_hist%hvars(ih_nh4uptake_si)%r81d(io_si)  + &
                          ccohort%daily_nh4_uptake*uconv

                     fates_hist%hvars(ih_no3uptake_si)%r81d(io_si) =       &
                          fates_hist%hvars(ih_no3uptake_si)%r81d(io_si)  + &
                          ccohort%daily_no3_uptake*uconv

                     ! Symbiotic Fixation
                     fates_hist%hvars(ih_nfix_si)%r81d(io_si) = &
                          fates_hist%hvars(ih_nfix_si)%r81d(io_si) + &
                          ccohort%sym_nfix_daily*uconv

                     ! Efflux/exudation
                     this%hvars(ih_nefflux_si)%r81d(io_si) = &
                          this%hvars(ih_nefflux_si)%r81d(io_si) + & 
                          ccohort%daily_n_efflux*uconv

                     ! Demand
                     this%hvars(ih_ndemand_si)%r81d(io_si) = &
                          this%hvars(ih_ndemand_si)%r81d(io_si) + &
                          ccohort%daily_n_demand*uconv
                     
                  case (phosphorus_element)

                     ! Mineralized uptake of PO4
                     fates_hist%hvars(ih_puptake_si)%r81d(io_si) =       &
                          fates_hist%hvars(ih_puptake_si)%r81d(io_si)  + &
                          ccohort%daily_p_gain*uconv

                     ! Efflux
                     this%hvars(ih_pefflux_si)%r81d(io_si) = &
                          this%hvars(ih_pefflux_si)%r81d(io_si) + &
                          ccohort%daily_p_efflux*uconv

                     ! Demand
                     this%hvars(ih_pdemand_si)%r81d(io_si) = &
                          this%hvars(ih_pdemand_si)%r81d(io_si) + & 
                          ccohort%daily_p_demand*uconv
                     
                  end select
               end do

               ccohort => ccohort%taller
            end do

            cpatch => cpatch%older
         end do

    end if if_dynam1

    ! Process multiplexed variables
    ! ---------------------------------------------------------------------------------------------
    
    if_dynam2: if(hlm_hist_level_dynam>1) then

         ! history site index
         io_si  = csite%h_gid

         call this%zero_site_hvars(csite,upfreq_in=group_nflx_complx)
         
         cpatch => csite%youngest_patch
         do while(associated(cpatch))

            ccohort => cpatch%shortest
            do while(associated(ccohort))

               ! If this is a new cohort, do not make diagnostics
               if(ccohort%isnew) then
                  ccohort => ccohort%taller
                  cycle
               end if

               ! size class index
               iscpf = ccohort%size_by_pft_class

               ! layer by size by pft index
               iclscpf = get_layersizetype_class_index(ccohort%canopy_layer,ccohort%dbh,ccohort%pft)

               ! unit conversion factor to get x/plant/day -> x/m2/sec
               uconv = ccohort%n * ha_per_m2 * days_per_sec

               ! Loop over the different elements. 
               do el = 1, num_elements

                  select case (element_list(el))

                  case (nitrogen_element) 

                     ! Mineralized uptake of NH4, NO3
                     fates_hist%hvars(ih_nh4uptake_scpf)%r82d(io_si,iscpf) =           &
                          fates_hist%hvars(ih_nh4uptake_scpf)%r82d(io_si,iscpf) +      &
                          ccohort%daily_nh4_uptake*uconv

                     fates_hist%hvars(ih_no3uptake_scpf)%r82d(io_si,iscpf) =           &
                          fates_hist%hvars(ih_no3uptake_scpf)%r82d(io_si,iscpf) +      &
                          ccohort%daily_no3_uptake*uconv

                     ! Fixation
                     fates_hist%hvars(ih_nfix_scpf)%r82d(io_si,iscpf) =           &
                          fates_hist%hvars(ih_nfix_scpf)%r82d(io_si,iscpf) +      &
                          ccohort%sym_nfix_daily*uconv

                     ! Efflux/exudation
                     this%hvars(ih_nefflux_scpf)%r82d(io_si,iscpf) = &
                          this%hvars(ih_nefflux_scpf)%r82d(io_si,iscpf) + &
                          ccohort%daily_n_efflux*uconv

                     ! Demand
                     this%hvars(ih_ndemand_scpf)%r82d(io_si,iscpf) = &
                          this%hvars(ih_ndemand_scpf)%r82d(io_si,iscpf) + &
                          ccohort%daily_n_demand*uconv
                     
                  case (phosphorus_element)

                     ! Mineralized uptake of PO4
                     fates_hist%hvars(ih_puptake_scpf)%r82d(io_si,iscpf) =             &
                          fates_hist%hvars(ih_puptake_scpf)%r82d(io_si,iscpf) +        &
                          ccohort%daily_p_gain*uconv

                     ! Efflux
                     this%hvars(ih_pefflux_scpf)%r82d(io_si,iscpf) = &
                          this%hvars(ih_pefflux_scpf)%r82d(io_si,iscpf) + & 
                          ccohort%daily_p_efflux*uconv
                     
                     ! Demand
                     this%hvars(ih_pdemand_scpf)%r82d(io_si,iscpf) = &
                          this%hvars(ih_pdemand_scpf)%r82d(io_si,iscpf) + &
                          ccohort%daily_p_demand*uconv
                     
                  end select
               end do

               ccohort => ccohort%taller
            end do

            cpatch => cpatch%older
         end do

    end if if_dynam2
    
    return
  end subroutine update_history_nutrflux

  ! ====================================================================================

  subroutine update_history_dyn(this,nc,nsites,sites,bc_in)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after Ecosystem Dynamics have been processed.
    ! This is the general routine that will call the single or multi-dimensional
    ! routines if they are called for by the user
    ! ---------------------------------------------------------------------------------

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    ! Locals

    ! If we don't have dynamics turned on, we just abort these diagnostics

    ! There is future work slated to split dynamics diagnostics into those
    ! related to states, and those related to fluxes. States should be fine
    ! to report in ST3 mode.
    
    if (hlm_use_ed_st3.eq.itrue) return

    if(hlm_hist_level_dynam>0) then
       call update_history_dyn1(this,nc,nsites,sites,bc_in)
       if(hlm_hist_level_dynam>1) then
          call update_history_dyn2(this,nc,nsites,sites,bc_in)
       end if
    end if

    
    
    return
  end subroutine update_history_dyn

  ! =========================================================================

  subroutine update_history_dyn1(this,nc,nsites,sites,bc_in)



    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)

    type(fates_cohort_type), pointer :: ccohort
    type(fates_patch_type),  pointer :: cpatch
    type(elem_diag_type), pointer :: elflux_diags_c ! Pointer to site level carbon fluxes
    type(litter_type), pointer :: litt     ! Generic pointer to any litter pool

    integer  :: s                  ! site counter
    integer  :: io_si              ! site's index in the history output array space
    integer  :: el                 ! element index
    integer  :: ft                 ! pft index
    real(r8) :: site_ba            ! Site basal area used for weighting
    real(r8) :: cohort_ba          ! Cohort basal area
    real(r8) :: site_ca            ! Site crown area used for weighting
    real(r8) :: store_max          ! Maximum storage capacity for carbon and nutrients
    real(r8) :: sapw_m             ! Sapwood mass (elemental, c,n or p) [kg/plant]
    real(r8) :: struct_m           ! Structural mass ""
    real(r8) :: leaf_m             ! Leaf mass ""
    real(r8) :: fnrt_m             ! Fineroot mass ""
    real(r8) :: store_m            ! Storage mass ""
    real(r8) :: alive_m            ! Alive biomass (sap+leaf+fineroot+repro+storage) ""
    real(r8) :: total_m            ! Total vegetation mass
    real(r8) :: repro_m            ! Total reproductive mass (on plant) ""
    real(r8) :: sapw_m_turnover    ! sapwood turnover rate [kg/yr]
    real(r8) :: store_m_turnover   ! storage turnover rate [kg/yr]
    real(r8) :: leaf_m_turnover    ! leaf turnover rate [kg/yr]
    real(r8) :: fnrt_m_turnover    ! fine-root turnover rate [kg/yr]
    real(r8) :: struct_m_turnover  ! structural turnover rate [kg/yr]
    real(r8) :: sapw_m_net_alloc   ! mass allocated to sapwood [kg/yr]
    real(r8) :: store_m_net_alloc  ! mass allocated to storage [kg/yr]
    real(r8) :: leaf_m_net_alloc   ! mass allocated to leaf [kg/yr]
    real(r8) :: fnrt_m_net_alloc   ! mass allocated to fine-root [kg/yr]
    real(r8) :: struct_m_net_alloc ! mass allocated to structure [kg/yr]
    real(r8) :: repro_m_net_alloc  ! mass allocated to reproduction [kg/yr]
    real(r8) :: leaf_herbivory     ! mass of leaves eaten by herbivores [kg/yr]
    real(r8) :: n_perm2            ! abundance per m2
    real(r8) :: area_frac  ! Fraction of area for this patch
    
    associate( hio_npatches_si         => this%hvars(ih_npatches_si)%r81d, &
         hio_ncohorts_si         => this%hvars(ih_ncohorts_si)%r81d, &
         hio_trimming_si         => this%hvars(ih_trimming_si)%r81d, &
         hio_area_plant_si       => this%hvars(ih_area_plant_si)%r81d, &
         hio_area_trees_si  => this%hvars(ih_area_trees_si)%r81d, &
         hio_fates_fraction_si   => this%hvars(ih_fates_fraction_si)%r81d, &
         hio_ba_weighted_height_si  => this%hvars(ih_ba_weighted_height_si)%r81d, &
         hio_ca_weighted_height_si  => this%hvars(ih_ca_weighted_height_si)%r81d, &
         hio_canopy_spread_si    => this%hvars(ih_canopy_spread_si)%r81d, &
         hio_nesterov_fire_danger_si => this%hvars(ih_nesterov_fire_danger_si)%r81d, &
         hio_fire_nignitions_si => this%hvars(ih_fire_nignitions_si)%r81d, &
         hio_fire_fdi_si => this%hvars(ih_fire_fdi_si)%r81d, &
         hio_spitfire_ros_si     => this%hvars(ih_spitfire_ros_si)%r81d, &
         hio_tfc_ros_si          => this%hvars(ih_tfc_ros_si)%r81d, &
         hio_effect_wspeed_si    => this%hvars(ih_effect_wspeed_si)%r81d, &
         hio_fire_intensity_si   => this%hvars(ih_fire_intensity_si)%r81d, &
         hio_fire_intensity_area_product_si => this%hvars(ih_fire_intensity_area_product_si)%r81d, &
         hio_fire_area_si        => this%hvars(ih_fire_area_si)%r81d, &
         hio_fire_fuel_bulkd_si  => this%hvars(ih_fire_fuel_bulkd_si)%r81d, &
         hio_fire_fuel_eff_moist_si => this%hvars(ih_fire_fuel_eff_moist_si)%r81d, &
         hio_fire_fuel_sav_si    => this%hvars(ih_fire_fuel_sav_si)%r81d, &
         hio_fire_fuel_mef_si    => this%hvars(ih_fire_fuel_mef_si)%r81d, &
         hio_sum_fuel_si         => this%hvars(ih_sum_fuel_si)%r81d,  &
         hio_litter_in_si        => this%hvars(ih_litter_in_si)%r81d, &
         hio_litter_out_si       => this%hvars(ih_litter_out_si)%r81d, &
         hio_npp_si              => this%hvars(ih_npp_si)%r81d, &
         hio_aresp_si            => this%hvars(ih_aresp_si)%r81d, &
         hio_growth_resp_si      => this%hvars(ih_growth_resp_si)%r81d, &
         hio_seed_bank_si        => this%hvars(ih_seed_bank_si)%r81d, &
         hio_ungerm_seed_bank_si => this%hvars(ih_ungerm_seed_bank_si)%r81d, &
         hio_seedling_pool_si    => this%hvars(ih_seedling_pool_si)%r81d, &
         hio_seeds_in_si         => this%hvars(ih_seeds_in_si)%r81d, &
         hio_seeds_in_local_si   => this%hvars(ih_seeds_in_local_si)%r81d, &
         hio_bdead_si            => this%hvars(ih_bdead_si)%r81d, &
         hio_balive_si           => this%hvars(ih_balive_si)%r81d, &
         hio_agb_si              => this%hvars(ih_agb_si)%r81d, &
         hio_canopy_biomass_si   => this%hvars(ih_canopy_biomass_si)%r81d, &
         hio_ustory_biomass_si   => this%hvars(ih_understory_biomass_si)%r81d, &
         hio_primaryland_fusion_error_si    => this%hvars(ih_primaryland_fusion_error_si)%r81d, &
         hio_fire_disturbance_rate_si      => this%hvars(ih_fire_disturbance_rate_si)%r81d, &
         hio_logging_disturbance_rate_si   => this%hvars(ih_logging_disturbance_rate_si)%r81d, &
         hio_fall_disturbance_rate_si      => this%hvars(ih_fall_disturbance_rate_si)%r81d, &
         hio_harvest_debt_si     => this%hvars(ih_harvest_debt_si)%r81d, &
         hio_harvest_debt_sec_si => this%hvars(ih_harvest_debt_sec_si)%r81d, &
         hio_npp_leaf_si         => this%hvars(ih_npp_leaf_si)%r81d, &
         hio_npp_seed_si         => this%hvars(ih_npp_seed_si)%r81d, &
         hio_npp_stem_si         => this%hvars(ih_npp_stem_si)%r81d, &
         hio_npp_froot_si        => this%hvars(ih_npp_froot_si)%r81d, &
         hio_npp_croot_si        => this%hvars(ih_npp_croot_si)%r81d, &
         hio_npp_stor_si         => this%hvars(ih_npp_stor_si)%r81d, &
         hio_grazing_si          => this%hvars(ih_grazing_si)%r81d, &
         hio_canopy_mortality_crownarea_si     => this%hvars(ih_canopy_mortality_crownarea_si)%r81d, &
         hio_ustory_mortality_crownarea_si => this%hvars(ih_understory_mortality_crownarea_si)%r81d, &
         hio_fire_c_to_atm_si  => this%hvars(ih_fire_c_to_atm_si)%r81d, &
         hio_demotion_carbonflux_si        => this%hvars(ih_demotion_carbonflux_si)%r81d, &
         hio_promotion_carbonflux_si       => this%hvars(ih_promotion_carbonflux_si)%r81d, &
         hio_canopy_mortality_carbonflux_si     => this%hvars(ih_canopy_mortality_carbonflux_si)%r81d, &
         hio_ustory_mortality_carbonflux_si => this%hvars(ih_understory_mortality_carbonflux_si)%r81d, &
         hio_woodproduct_si                 => this%hvars(ih_woodproduct_si)%r81d, &
         hio_gdd_si                           => this%hvars(ih_gdd_si)%r81d, &
         hio_site_ncolddays_si                => this%hvars(ih_site_ncolddays_si)%r81d, &
         hio_site_nchilldays_si               => this%hvars(ih_site_nchilldays_si)%r81d, &
         hio_site_cstatus_si                  => this%hvars(ih_site_cstatus_si)%r81d, &
         hio_cleafoff_si                      => this%hvars(ih_cleafoff_si)%r81d, &
         hio_cleafon_si                       => this%hvars(ih_cleafon_si)%r81d, &
         hio_cbal_err_fates_si                => this%hvars(ih_cbal_err_fates_si)%r81d, &
         hio_tveg24                           => this%hvars(ih_tveg24_si)%r81d, &
         hio_tlongterm                        => this%hvars(ih_tlongterm_si)%r81d, &
         hio_tgrowth                          => this%hvars(ih_tgrowth_si)%r81d, &
         hio_lai_si                           => this%hvars(ih_lai_si)%r81d, &
         hio_elai_si                          => this%hvars(ih_elai_si)%r81d, &
         hio_harvest_woodprod_carbonflux_si   => this%hvars(ih_harvest_woodprod_carbonflux_si)%r81d, &
         hio_luchange_woodprod_carbonflux_si => this%hvars(ih_luchange_woodprod_carbonflux_si)%r81d)

      ! ---------------------------------------------------------------------------------
      ! Loop through the FATES scale hierarchy and fill the history IO arrays
      ! ---------------------------------------------------------------------------------

      siteloop: do s = 1,nsites

         io_si  = sites(s)%h_gid

         site_ba = 0._r8
         site_ca = 0._r8

         call this%zero_site_hvars(sites(s),upfreq_in=group_dyna_simple)
         
         ! set the fates fraction to one, since it is zero on non-fates columns, &
         ! the average is the total gridcell fates fraction
         hio_fates_fraction_si(io_si) = 1._r8

         ! Total carbon model error [kgC/day -> kgC/s]
         hio_cbal_err_fates_si(io_si) = &
              sites(s)%mass_balance(element_pos(carbon12_element))%err_fates / sec_per_day

         
         
         ! Total carbon lost to atmosphere from burning (kgC/site/day -> kgC/m2/s)
         hio_fire_c_to_atm_si(io_si) = &
              sites(s)%mass_balance(element_pos(carbon12_element))%burn_flux_to_atm * &
              ha_per_m2 * days_per_sec

         ! damage variables - site level - this needs to be OUT of the patch loop 
         if(hlm_use_tree_damage .eq. itrue) then

            this%hvars(ih_crownarea_canopy_damage_si)%r81d(io_si) = &
                 this%hvars(ih_crownarea_canopy_damage_si)%r81d(io_si) + &
                 sites(s)%crownarea_canopy_damage  * days_per_year * 1 / m2_per_ha

            this%hvars(ih_crownarea_ustory_damage_si)%r81d(io_si) = &
                 this%hvars(ih_crownarea_ustory_damage_si)%r81d(io_si) + &
                 sites(s)%crownarea_ustory_damage  * days_per_year * 1 / m2_per_ha

         end if

         ! Canopy spread index (0-1)
         hio_canopy_spread_si(io_si) = sites(s)%spread

         ! Update the site status for cold deciduous (drought deciduous is now PFT dependent)
         hio_site_cstatus_si(io_si)   = real(sites(s)%cstatus,r8)

         ! Number of chill days and cold days
         hio_site_nchilldays_si(io_si) = real(sites(s)%nchilldays,r8)
         hio_site_ncolddays_si(io_si)  = real(sites(s)%ncolddays,r8)

         ! Growing degree-days
         hio_gdd_si(io_si) = sites(s)%grow_deg_days

         ! Model days elapsed since leaf on/off for cold-deciduous
         hio_cleafoff_si(io_si) = real(sites(s)%phen_model_date - sites(s)%cleafoffdate,r8)
         hio_cleafon_si(io_si)  = real(sites(s)%phen_model_date - sites(s)%cleafondate,r8)

         ! track total wood product accumulation at the site level
         hio_woodproduct_si(io_si) = sites(s)%resources_management%trunk_product_site &
              * AREA_INV

         ! site-level fire variables:

         ! Nesterov index (unitless)
         hio_nesterov_fire_danger_si(io_si) = sites(s)%fireWeather%fire_weather_index
         
         hio_effect_wspeed_si(io_si) = sites(s)%fireWeather%effective_windspeed/sec_per_min

         ! number of ignitions [#/km2/day -> #/m2/s]
         hio_fire_nignitions_si(io_si) = sites(s)%NF_successful / m2_per_km2 /  &
              sec_per_day

         ! Fire danger index (FDI) (0-1)
         hio_fire_fdi_si(io_si) = sites(s)%FDI

         ! If hydraulics are turned on, track the error terms associated with
         ! dynamics [kg/m2]
         if(hlm_use_planthydro.eq.itrue)then
            this%hvars(ih_h2oveg_dead_si)%r81d(io_si)         = sites(s)%si_hydr%h2oveg_dead
            this%hvars(ih_h2oveg_recruit_si)%r81d(io_si)      = sites(s)%si_hydr%h2oveg_recruit
            this%hvars(ih_h2oveg_growturn_err_si)%r81d(io_si) = sites(s)%si_hydr%h2oveg_growturn_err
         end if

         hio_harvest_debt_si(io_si) = sites(s)%resources_management%harvest_debt
         hio_harvest_debt_sec_si(io_si) = sites(s)%resources_management%harvest_debt_sec

         ! error in primary lands from patch fusion [m2 m-2 day-1] -> [m2 m-2 yr-1]
         hio_primaryland_fusion_error_si(io_si) = sites(s)%primary_land_patchfusion_error * days_per_year

         ! output site-level disturbance rates [m2 m-2 day-1] -> [m2 m-2 yr-1] - TO DO rework this

         hio_fire_disturbance_rate_si(io_si) = &
              sum(sites(s)%disturbance_rates(dtype_ifire,1:n_landuse_cats,1:n_landuse_cats)) * &
              days_per_year

         hio_logging_disturbance_rate_si(io_si) = &
              sum(sites(s)%disturbance_rates(dtype_ilog,1:n_landuse_cats,1:n_landuse_cats)) * &
              days_per_year

         hio_fall_disturbance_rate_si(io_si) = &
              sum(sites(s)%disturbance_rates(dtype_ifall,1:n_landuse_cats,1:n_landuse_cats)) * &
              days_per_year

         hio_harvest_woodprod_carbonflux_si(io_si) = AREA_INV * &
              sum(sites(s)%mass_balance(element_pos(carbon12_element))%wood_product_harvest(1:numpft))

         hio_luchange_woodprod_carbonflux_si(io_si) = AREA_INV * &
              sum(sites(s)%mass_balance(element_pos(carbon12_element))%wood_product_landusechange(1:numpft))
          

         ! carbon flux associated with mortality of trees dying by fire
         hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
              sum(sites(s)%fmort_carbonflux_canopy(:)) / g_per_kg

         
         hio_ustory_mortality_carbonflux_si(io_si) = hio_ustory_mortality_carbonflux_si(io_si) + &
              sum(sites(s)%fmort_carbonflux_ustory(:)) / g_per_kg

         ! treat carbon flux from imort the same way
         hio_ustory_mortality_carbonflux_si(io_si) = hio_ustory_mortality_carbonflux_si(io_si) + &
              sum(sites(s)%imort_carbonflux(:))

         ! convert kg C / ha / day to kgc / m2 / sec
         hio_demotion_carbonflux_si(io_si) = sites(s)%demotion_carbonflux * ha_per_m2 * days_per_sec
         hio_promotion_carbonflux_si(io_si) = sites(s)%promotion_carbonflux * ha_per_m2 * days_per_sec
         !
         ! mortality-associated carbon fluxes

         hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
              sum(sites(s)%term_carbonflux_canopy(:,:)) * days_per_sec * ha_per_m2
         
         hio_ustory_mortality_carbonflux_si(io_si) = hio_ustory_mortality_carbonflux_si(io_si) + &
              sum(sites(s)%term_carbonflux_ustory(:,:)) * days_per_sec * ha_per_m2

         ! add site level mortality counting to crownarea diagnostic
         hio_canopy_mortality_crownarea_si(io_si) = hio_canopy_mortality_crownarea_si(io_si) + &
              sites(s)%fmort_crownarea_canopy + &
              sites(s)%term_crownarea_canopy * days_per_year

         hio_ustory_mortality_crownarea_si(io_si) = hio_ustory_mortality_crownarea_si(io_si) + &
              sites(s)%fmort_crownarea_ustory + &
              sites(s)%term_crownarea_ustory * days_per_year + &
              sites(s)%imort_crownarea


         elflux_diags_c => sites(s)%flux_diags%elem(element_pos(carbon12_element))

         hio_litter_in_si(io_si) = (sum(elflux_diags_c%cwd_ag_input(:)) + &
              sum(elflux_diags_c%cwd_bg_input(:)) + &
              sum(elflux_diags_c%surf_fine_litter_input(:)) + &
              sum(elflux_diags_c%root_litter_input(:))) * &
              AREA_INV * days_per_sec

         ! Loop through patches to sum up diagonistics
         cpatch => sites(s)%oldest_patch
         patchloop: do while(associated(cpatch))

            ! Increment the number of patches per site
            hio_npatches_si(io_si) = hio_npatches_si(io_si) + 1._r8

            hio_lai_si(io_si) = hio_lai_si(io_si) + sum( cpatch%canopy_area_profile(:,:,:) * cpatch%tlai_profile(:,:,:) ) * &
                 cpatch%total_canopy_area * AREA_INV
            
            hio_elai_si(io_si) = hio_elai_si(io_si) + sum( cpatch%canopy_area_profile(:,:,:) * cpatch%elai_profile(:,:,:) ) * &
                 cpatch%total_canopy_area * AREA_INV
            
            ! 24hr veg temperature
            hio_tveg24(io_si) = hio_tveg24(io_si) + &
                 (cpatch%tveg24%GetMean()- t_water_freeze_k_1atm)*cpatch%area*AREA_INV

            ! long-term veg temperature
            hio_tlongterm(io_si) = hio_tlongterm(io_si) + &
                 (cpatch%tveg_longterm%GetMean()- t_water_freeze_k_1atm)*cpatch%area*AREA_INV

            ! long-term running mean veg temperature (tgrowth)
            hio_tgrowth(io_si) = hio_tgrowth(io_si) + &
                 (cpatch%tveg_lpa%GetMean()- t_water_freeze_k_1atm)*cpatch%area*AREA_INV

            ! Canopy trimming - degree to which canopy expansion is limited by leaf economics (0-1)
            if(associated(cpatch%tallest))then
               hio_trimming_si(io_si) = hio_trimming_si(io_si) + cpatch%tallest%canopy_trim * cpatch%area * AREA_INV
            endif

            ! area occupied by plants and trees [m2/m2]
            hio_area_plant_si(io_si) = hio_area_plant_si(io_si) + min(cpatch%total_canopy_area,cpatch%area) * AREA_INV
            hio_area_trees_si(io_si) = hio_area_trees_si(io_si) + min(cpatch%total_tree_area,cpatch%area) * AREA_INV

            ! Patch specific variables that are already calculated
            ! These things are all duplicated. Should they all be converted to LL or array structures RF?
            ! define scalar to counteract the patch albedo scaling logic for conserved quantities

            ! Update Fire Variables
            hio_spitfire_ros_si(io_si)         = hio_spitfire_ros_si(io_si) + cpatch%ROS_front * cpatch%area * AREA_INV / sec_per_min
            hio_tfc_ros_si(io_si)              = hio_tfc_ros_si(io_si) + cpatch%TFC_ROS * cpatch%area * AREA_INV
            hio_fire_intensity_si(io_si)       = hio_fire_intensity_si(io_si) + cpatch%FI * cpatch%area * AREA_INV * J_per_kJ
            hio_fire_area_si(io_si)            = hio_fire_area_si(io_si) + cpatch%frac_burnt * cpatch%area * AREA_INV / sec_per_day
            hio_fire_fuel_bulkd_si(io_si)      = hio_fire_fuel_bulkd_si(io_si) + cpatch%fuel%bulk_density_notrunks * cpatch%area * AREA_INV
            hio_fire_fuel_eff_moist_si(io_si)  = hio_fire_fuel_eff_moist_si(io_si) + cpatch%fuel%average_moisture_notrunks * cpatch%area * AREA_INV
            hio_fire_fuel_sav_si(io_si)        = hio_fire_fuel_sav_si(io_si) + cpatch%fuel%SAV_notrunks * cpatch%area * AREA_INV / m_per_cm
            hio_fire_fuel_mef_si(io_si)        = hio_fire_fuel_mef_si(io_si) + cpatch%fuel%MEF_notrunks * cpatch%area * AREA_INV
            hio_sum_fuel_si(io_si)             = hio_sum_fuel_si(io_si) + cpatch%fuel%non_trunk_loading * cpatch%area * AREA_INV

            hio_fire_intensity_area_product_si(io_si) = hio_fire_intensity_area_product_si(io_si) + &
                 cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV * J_per_kJ

            litt => cpatch%litter(element_pos(carbon12_element))

            area_frac = cpatch%area * AREA_INV

            ! Sum up all output fluxes (fragmentation) kgC/m2/day -> kgC/m2/s
            hio_litter_out_si(io_si) = hio_litter_out_si(io_si) + &
                 (sum(litt%leaf_fines_frag(:)) + &
                 sum(litt%root_fines_frag(:,:)) + &
                 sum(litt%ag_cwd_frag(:)) + &
                 sum(litt%bg_cwd_frag(:,:)) + &
                 sum(litt%seed_decay(:)) + &
                 sum(litt%seed_germ_decay(:))) * &
                 area_frac * days_per_sec

            ! Sum up total seed bank (germinated and ungerminated)
            hio_seed_bank_si(io_si) = hio_seed_bank_si(io_si) + &
                 (sum(litt%seed(:))+sum(litt%seed_germ(:))) * &
                 area_frac

            ! Sum up total seed bank (just ungerminated)
            hio_ungerm_seed_bank_si(io_si) = hio_ungerm_seed_bank_si(io_si) + &
                 sum(litt%seed(:)) * area_frac

            ! Sum up total seedling pool  
            hio_seedling_pool_si(io_si) = hio_seedling_pool_si(io_si) + &
                 sum(litt%seed_germ(:)) * area_frac

            ! Sum up the input flux into the seed bank (local and external)
            hio_seeds_in_si(io_si) = hio_seeds_in_si(io_si) + &
                 (sum(litt%seed_in_local(:)) + sum(litt%seed_in_extern(:))) * &
                 area_frac * days_per_sec

            hio_seeds_in_local_si(io_si) = hio_seeds_in_local_si(io_si) + &
                 sum(litt%seed_in_local(:)) * &
                 area_frac * days_per_sec

            ! loop through cohorts on patch
            ccohort => cpatch%shortest
            cohortloop: do while(associated(ccohort))

               ft = ccohort%pft
               n_perm2 = ccohort%n * AREA_INV

               ! Increment the number of cohorts per site
               hio_ncohorts_si(io_si) = hio_ncohorts_si(io_si) + 1._r8

               ! Update biomass components
               ! Mass pools [kg]
               elloop: do el = 1, num_elements

                  sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
                  struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
                  leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
                  fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
                  store_m  = ccohort%prt%GetState(store_organ, element_list(el))
                  repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))

                  alive_m  = leaf_m + fnrt_m + sapw_m
                  total_m  = alive_m + store_m + struct_m

                  ! Plant multi-element states and fluxes
                  ! Zero states, and set the fluxes
                  if( element_list(el).eq.carbon12_element )then

                     ! mass in different tissues [kg/ha] -> [kg/m2]
                     this%hvars(ih_storec_si)%r81d(io_si) =                       &
                          this%hvars(ih_storec_si)%r81d(io_si) + ccohort%n *        &
                          store_m / m2_per_ha
                     this%hvars(ih_leafc_si)%r81d(io_si) =                        &
                          this%hvars(ih_leafc_si)%r81d(io_si) + ccohort%n *         &
                          leaf_m / m2_per_ha
                     this%hvars(ih_fnrtc_si)%r81d(io_si) =                        &
                          this%hvars(ih_fnrtc_si)%r81d(io_si) + ccohort%n *         &
                          fnrt_m / m2_per_ha
                     this%hvars(ih_reproc_si)%r81d(io_si) =                       &
                          this%hvars(ih_reproc_si)%r81d(io_si)+ ccohort%n *         &
                          repro_m / m2_per_ha
                     this%hvars(ih_sapwc_si)%r81d(io_si) =                        &
                          this%hvars(ih_sapwc_si)%r81d(io_si) + ccohort%n *         &
                          sapw_m / m2_per_ha
                     this%hvars(ih_totvegc_si)%r81d(io_si) =                      &
                          this%hvars(ih_totvegc_si)%r81d(io_si)+ ccohort%n *        &
                          total_m / m2_per_ha

                     call bstore_allom(ccohort%dbh,ccohort%pft,ccohort%crowndamage,ccohort%canopy_trim, &
                          store_max)

                     this%hvars(ih_storectfrac_si)%r81d(io_si)  = &
                          this%hvars(ih_storectfrac_si)%r81d(io_si) + ccohort%n * store_max/m2_per_ha

                     hio_bdead_si(io_si) = hio_bdead_si(io_si)  + n_perm2 * struct_m
                     hio_balive_si(io_si) = hio_balive_si(io_si) + n_perm2 * alive_m

                     hio_agb_si(io_si) = hio_agb_si(io_si) + n_perm2 *            &
                          ( leaf_m + (sapw_m + struct_m + store_m) * prt_params%allom_agb_frac(ccohort%pft) )

                     if( hlm_parteh_mode == prt_cnp_flex_allom_hyp) then
                        this%hvars(ih_l2fr_si)%r81d(io_si) = &
                             this%hvars(ih_l2fr_si)%r81d(io_si) + &
                             ccohort%l2fr *ccohort%n * fnrt_m / m2_per_ha
                     else
                        this%hvars(ih_l2fr_si)%r81d(io_si) = &
                             this%hvars(ih_l2fr_si)%r81d(io_si) + &
                             prt_params%allom_l2fr(ft) *ccohort%n * fnrt_m / m2_per_ha
                     end if

                  elseif(element_list(el).eq.nitrogen_element)then

                     store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                     this%hvars(ih_storen_si)%r81d(io_si)  =                      &
                          this%hvars(ih_storen_si)%r81d(io_si) + ccohort%n *        &
                          store_m / m2_per_ha
                     this%hvars(ih_storentfrac_si)%r81d(io_si)  =                 &
                          this%hvars(ih_storentfrac_si)%r81d(io_si) + ccohort%n *   &
                          store_max / m2_per_ha
                     this%hvars(ih_leafn_si)%r81d(io_si)   =                      &
                          this%hvars(ih_leafn_si)%r81d(io_si) + ccohort%n *         &
                          leaf_m / m2_per_ha
                     this%hvars(ih_fnrtn_si)%r81d(io_si)   =                      &
                          this%hvars(ih_fnrtn_si)%r81d(io_si) + ccohort%n *         &
                          fnrt_m / m2_per_ha
                     this%hvars(ih_repron_si)%r81d(io_si)  =                      &
                          this%hvars(ih_repron_si)%r81d(io_si) + ccohort%n *        &
                          repro_m / m2_per_ha
                     this%hvars(ih_sapwn_si)%r81d(io_si)   =                      &
                          this%hvars(ih_sapwn_si)%r81d(io_si) + ccohort%n *         &
                          sapw_m / m2_per_ha
                     this%hvars(ih_totvegn_si)%r81d(io_si) =                      &
                          this%hvars(ih_totvegn_si)%r81d(io_si) + ccohort%n *       &
                          total_m / m2_per_ha

                  elseif(element_list(el).eq.phosphorus_element) then

                     store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                     this%hvars(ih_storep_si)%r81d(io_si)  =                      &
                          this%hvars(ih_storep_si)%r81d(io_si) + ccohort%n *        &
                          store_m / m2_per_ha
                     this%hvars(ih_storeptfrac_si)%r81d(io_si)  =                 &
                          this%hvars(ih_storeptfrac_si)%r81d(io_si) + ccohort%n *   &
                          store_max / m2_per_ha
                     this%hvars(ih_leafp_si)%r81d(io_si)   =                      &
                          this%hvars(ih_leafp_si)%r81d(io_si) + ccohort%n *         &
                          leaf_m / m2_per_ha
                     this%hvars(ih_fnrtp_si)%r81d(io_si)   =                      &
                          this%hvars(ih_fnrtp_si)%r81d(io_si) + ccohort%n *         &
                          fnrt_m / m2_per_ha
                     this%hvars(ih_reprop_si)%r81d(io_si)  =                      &
                          this%hvars(ih_reprop_si)%r81d(io_si) + ccohort%n *        &
                          repro_m / m2_per_ha
                     this%hvars(ih_sapwp_si)%r81d(io_si)   =                      &
                          this%hvars(ih_sapwp_si)%r81d(io_si) + ccohort%n *         &
                          sapw_m / m2_per_ha
                     this%hvars(ih_totvegp_si)%r81d(io_si) =                      &
                          this%hvars(ih_totvegp_si)%r81d(io_si)+ ccohort%n *        &
                          total_m / m2_per_ha
                  end if
               end do elloop

               ! FLUXES ---
               ! Flux Variables (cohorts must had experienced a day before any of these values
               ! have any meaning, otherwise they are just inialization values
               notnew: if( .not.(ccohort%isnew) ) then

                  hio_npp_si(io_si) = hio_npp_si(io_si) + &
                       ccohort%npp_acc_hold * n_perm2 / days_per_year / sec_per_day

                  hio_growth_resp_si(io_si) =  hio_growth_resp_si(io_si) + &
                       ccohort%resp_g_acc_hold * n_perm2 / days_per_year / sec_per_day

                  hio_aresp_si(io_si) = hio_aresp_si(io_si) + &
                       (ccohort%resp_g_acc_hold + ccohort%resp_m_acc_hold + &
                       ccohort%resp_excess_hold) * n_perm2 / days_per_year / sec_per_day

                  ! Turnover pools [kgC/day] * [day/yr] = [kgC/yr]
                  sapw_m_turnover   = ccohort%prt%GetTurnover(sapw_organ, carbon12_element) * days_per_year
                  store_m_turnover  = ccohort%prt%GetTurnover(store_organ, carbon12_element) * days_per_year
                  leaf_m_turnover   = ccohort%prt%GetTurnover(leaf_organ, carbon12_element) * days_per_year
                  fnrt_m_turnover   = ccohort%prt%GetTurnover(fnrt_organ, carbon12_element) * days_per_year
                  struct_m_turnover = ccohort%prt%GetTurnover(struct_organ, carbon12_element) * days_per_year

                  ! Net change from allocation and transport [kgC/day] * [day/yr] = [kgC/yr]
                  sapw_m_net_alloc   = ccohort%prt%GetNetAlloc(sapw_organ, carbon12_element) * days_per_year
                  store_m_net_alloc  = ccohort%prt%GetNetAlloc(store_organ, carbon12_element) * days_per_year
                  leaf_m_net_alloc   = ccohort%prt%GetNetAlloc(leaf_organ, carbon12_element) * days_per_year
                  fnrt_m_net_alloc   = ccohort%prt%GetNetAlloc(fnrt_organ, carbon12_element) * days_per_year
                  struct_m_net_alloc = ccohort%prt%GetNetAlloc(struct_organ, carbon12_element) * days_per_year
                  repro_m_net_alloc  = ccohort%prt%GetNetAlloc(repro_organ, carbon12_element) * days_per_year

                  ! ecosystem-level, organ-partitioned NPP/allocation fluxes
                  ! [kgC/yr] -> [kgC/sec]
                  hio_npp_leaf_si(io_si) = hio_npp_leaf_si(io_si) +               &
                       leaf_m_net_alloc * n_perm2 / days_per_year / sec_per_day
                  hio_npp_seed_si(io_si) = hio_npp_seed_si(io_si) +               &
                       repro_m_net_alloc * n_perm2 / days_per_year / sec_per_day
                  hio_npp_stem_si(io_si) = hio_npp_stem_si(io_si) +               &
                       (sapw_m_net_alloc + struct_m_net_alloc) * n_perm2 *          &
                       (prt_params%allom_agb_frac(ccohort%pft)) /                   &
                       days_per_year / sec_per_day
                  hio_npp_froot_si(io_si) = hio_npp_froot_si(io_si) +             &
                       fnrt_m_net_alloc * n_perm2 / days_per_year / sec_per_day
                  hio_npp_croot_si(io_si) = hio_npp_croot_si(io_si) +             &
                       (sapw_m_net_alloc + struct_m_net_alloc) * n_perm2 *          &
                       (1._r8-prt_params%allom_agb_frac(ccohort%pft)) /             &
                       days_per_year / sec_per_day
                  hio_npp_stor_si(io_si) = hio_npp_stor_si(io_si) +               &
                       store_m_net_alloc * n_perm2 / days_per_year / sec_per_day

                  leaf_herbivory   = ccohort%prt%GetHerbivory(leaf_organ, carbon12_element) * days_per_year  !cdkcdk
                  hio_grazing_si(io_si) = hio_grazing_si(io_si) + leaf_herbivory * n_perm2 / days_per_year / sec_per_day

                  ! Woody State Variables (basal area growth increment)
                  if ( prt_params%woody(ft) == itrue) then

                     cohort_ba = 0.25_r8*pi_const*((ccohort%dbh/100.0_r8)**2.0_r8)*ccohort%n

                     hio_ba_weighted_height_si(io_si) = hio_ba_weighted_height_si(io_si) + &
                          ccohort%height * cohort_ba

                     site_ba = site_ba + cohort_ba

                  end if

                  ! THIS NEEDS TO BE NORMALIZED (RGK)
                  hio_ca_weighted_height_si(io_si) = hio_ca_weighted_height_si(io_si) + &
                       ccohort%height * ccohort%c_area / m2_per_ha

                  site_ca = site_ca + ccohort%c_area / m2_per_ha

                  ! RGK - CANOPY/USTORY BIOMASS IS NOT A FLUX, NEED NOT BE CONDITIONED BY isnew
                  ! ----------------------------------------------------------------------------------
                  if (ccohort%canopy_layer .eq. 1) then
                     hio_canopy_biomass_si(io_si) = hio_canopy_biomass_si(io_si) + n_perm2 * total_m

                     hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
                          (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                          ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * &
                          total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2 + &
                          (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_m * &
                          ccohort%n * ha_per_m2

                     hio_canopy_mortality_crownarea_si(io_si) = hio_canopy_mortality_crownarea_si(io_si) + &
                          (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                          ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * &
                          ccohort%c_area  + &
                          (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                          ccohort%c_area * sec_per_day * days_per_year

                  else
                     hio_ustory_biomass_si(io_si) = hio_ustory_biomass_si(io_si) + n_perm2 * total_m

                     hio_ustory_mortality_carbonflux_si(io_si) = hio_ustory_mortality_carbonflux_si(io_si) + &
                          (ccohort%bmort + ccohort%hmort + ccohort%cmort +   &
                          ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * &
                          total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2 + &
                          (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_m * &
                          ccohort%n * ha_per_m2

                     hio_ustory_mortality_crownarea_si(io_si) = hio_ustory_mortality_crownarea_si(io_si) + &
                          (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                          ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * &
                          ccohort%c_area  + &
                          (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                          ccohort%c_area * sec_per_day * days_per_year

                  end if

               end if notnew

               ccohort => ccohort%taller
            enddo cohortloop ! cohort loop

            cpatch => cpatch%younger
         end do patchloop !patch loop

         ! Perform any necessary normalizations
         ! ----------------------------------------------------------------------------------------

         ! Normalize crown-area weighted height
         if(site_ca>nearzero)then
            hio_ca_weighted_height_si(io_si) = hio_ca_weighted_height_si(io_si)/site_ca
         end if

         ! divide basal-area-weighted height by basal area to get mean
         if ( site_ba .gt. nearzero ) then
            hio_ba_weighted_height_si(io_si) = hio_ba_weighted_height_si(io_si)/site_ba
         endif
         
         elloop2: do el = 1, num_elements
            if( element_list(el).eq.carbon12_element )then
               if( this%hvars(ih_storectfrac_si)%r81d(io_si)>nearzero ) then
                  this%hvars(ih_storectfrac_si)%r81d(io_si) = this%hvars(ih_storec_si)%r81d(io_si) / &
                       this%hvars(ih_storectfrac_si)%r81d(io_si)
               end if
            elseif( element_list(el).eq.nitrogen_element )then
               if( this%hvars(ih_storentfrac_si)%r81d(io_si)>nearzero ) then
                  this%hvars(ih_storentfrac_si)%r81d(io_si)  = this%hvars(ih_storen_si)%r81d(io_si) / &
                       this%hvars(ih_storentfrac_si)%r81d(io_si)
               end if
            elseif( element_list(el).eq.phosphorus_element )then
               if( this%hvars(ih_storeptfrac_si)%r81d(io_si)>nearzero ) then
                  this%hvars(ih_storeptfrac_si)%r81d(io_si) = this%hvars(ih_storep_si)%r81d(io_si) / &
                       this%hvars(ih_storeptfrac_si)%r81d(io_si)
               end if
            end if
         end do elloop2


         if(this%hvars(ih_fnrtc_si)%r81d(io_si)>nearzero)then
            this%hvars(ih_l2fr_si)%r81d(io_si) = this%hvars(ih_l2fr_si)%r81d(io_si) / &
                 this%hvars(ih_fnrtc_si)%r81d(io_si)
         else
            this%hvars(ih_l2fr_si)%r81d(io_si) = hlm_hio_ignore_val
         end if
         
         ! zero the site-level termination carbon flux variable
         sites(s)%term_carbonflux_canopy(:,:) = 0._r8
         sites(s)%term_carbonflux_ustory(:,:) = 0._r8
         sites(s)%crownarea_canopy_damage = 0._r8
         sites(s)%crownarea_ustory_damage = 0._r8

      end do siteloop

    end associate
    return
  end subroutine update_history_dyn1

  ! =========================================================================================

  subroutine update_history_dyn2(this,nc,nsites,sites,bc_in)

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)

    type(fates_cohort_type), pointer :: ccohort
    type(fates_patch_type),  pointer :: cpatch
    type(litter_type), pointer :: litt_c   ! Pointer to the carbon12 litter pool
    type(litter_type), pointer :: litt     ! Generic pointer to any litter pool
    integer  :: s                  ! site counter
    integer  :: ipa2           ! patch index matching host model array space
    integer  :: io_si              ! site's index in the history output array space
    integer  :: el                 ! element index
    integer  :: ft                 ! pft index
    real(r8) :: site_ba            ! Site basal area used for weighting
    integer  :: model_day_int      ! Integer model day since simulation start
    real(r8) :: store_max          ! Maximum storage capacity for carbon and nutrients
    real(r8) :: sapw_m             ! Sapwood mass (elemental, c,n or p) [kg/plant]
    real(r8) :: struct_m           ! Structural mass ""
    real(r8) :: leaf_m             ! Leaf mass ""
    real(r8) :: fnrt_m             ! Fineroot mass ""
    real(r8) :: store_m            ! Storage mass ""
    real(r8) :: alive_m            ! Alive biomass (sap+leaf+fineroot+repro+storage) ""
    real(r8) :: total_m            ! Total vegetation mass
    real(r8) :: repro_m            ! Total reproductive mass (on plant) ""
    real(r8) :: sapw_m_turnover    ! sapwood turnover rate [kg/yr]
    real(r8) :: store_m_turnover   ! storage turnover rate [kg/yr]
    real(r8) :: leaf_m_turnover    ! leaf turnover rate [kg/yr]
    real(r8) :: fnrt_m_turnover    ! fine-root turnover rate [kg/yr]
    real(r8) :: struct_m_turnover  ! structural turnover rate [kg/yr]
    real(r8) :: sapw_m_net_alloc   ! mass allocated to sapwood [kg/yr]
    real(r8) :: store_m_net_alloc  ! mass allocated to storage [kg/yr]
    real(r8) :: leaf_m_net_alloc   ! mass allocated to leaf [kg/yr]
    real(r8) :: fnrt_m_net_alloc   ! mass allocated to fine-root [kg/yr]
    real(r8) :: struct_m_net_alloc ! mass allocated to structure [kg/yr]
    real(r8) :: repro_m_net_alloc  ! mass allocated to reproduction [kg/yr]
    real(r8) :: n_perm2            ! abundance per m2
    integer  :: ageclass_since_anthrodist  ! what is the equivalent age class for
                                           ! time-since-anthropogenic-disturbance of secondary forest
    real(r8) :: area_frac  ! Fraction of area for this patch
    real(r8) :: frac_canopy_in_bin  ! fraction of a leaf's canopy that is within a given height bin
    real(r8) :: binbottom,bintop    ! edges of height bins
    integer  :: height_bin_max, height_bin_min   ! which height bin a given cohort's canopy is in
    integer  :: ican, ileaf, cnlf_indx  ! iterators for leaf and canopy level
    integer  :: elcwd, i_cwd            ! combined index of element and pft or cwd
    integer  :: i_scpf,i_pft,i_scls     ! iterators for scpf, pft, and scls dims
    integer  :: i_cacls, i_capf      ! iterators for cohort age and cohort age x pft
    integer  :: i_fuel            ! iterators for fuel dims
    integer  :: i_heightbin  ! iterator for height bins
    integer  :: iagepft     ! age x pft index
    integer  :: ilyr      ! Soil index for nlevsoil
    integer  :: iscag        ! size-class x age index
    integer  :: iscagpft     ! size-class x age x pft index
    integer  :: icdpf, icdsc, icdam ! iterators for the crown damage level
    integer  :: i_agefuel     ! age x fuel size class index
    real(r8) :: gpp_cached    ! gpp from previous timestep, for c13 discrimination
    real(r8) :: crown_depth   ! Depth of the crown [m]
    real(r8) :: gpp_cached_scpf(numpft*nlevsclass)  ! variable used to cache gpp value in previous time step; for C13 discrimination
    real(r8) :: storen_canopy_scpf(numpft*nlevsclass)
    real(r8) :: storen_understory_scpf(numpft*nlevsclass)
    real(r8) :: storep_canopy_scpf(numpft*nlevsclass)
    real(r8) :: storep_understory_scpf(numpft*nlevsclass)
    real(r8) :: storec_canopy_scpf(numpft*nlevsclass)
    real(r8) :: storec_understory_scpf(numpft*nlevsclass)
    real(r8) :: a_sapw ! sapwood area [m^2]
    real(r8) :: c_sapw ! sapwood biomass [kgC]

    integer  :: i_dist, j_dist

    type(elem_diag_type), pointer :: elflux_diags
    type(elem_diag_type), pointer :: elflux_diags_c


    real(r8), parameter :: reallytalltrees = 1000.   ! some large number (m)

    associate( hio_biomass_si_pft      => this%hvars(ih_biomass_si_pft)%r82d, &
         hio_leafbiomass_si_pft  => this%hvars(ih_leafbiomass_si_pft)%r82d, &
         hio_storebiomass_si_pft => this%hvars(ih_storebiomass_si_pft)%r82d, &
         hio_nindivs_si_pft      => this%hvars(ih_nindivs_si_pft)%r82d, &
         hio_recruitment_si_pft  => this%hvars(ih_recruitment_si_pft)%r82d, &
         hio_recruitment_cflux_si_pft  => this%hvars(ih_recruitment_cflux_si_pft)%r82d, &
         hio_seeds_out_gc_si_pft => this%hvars(ih_seeds_out_gc_si_pft)%r82d, &
         hio_seeds_in_gc_si_pft  => this%hvars(ih_seeds_in_gc_si_pft)%r82d, &
         hio_mortality_si_pft    => this%hvars(ih_mortality_si_pft)%r82d, &
         hio_mortality_carbonflux_si_pft  => this%hvars(ih_mortality_carbonflux_si_pft)%r82d, &
         hio_cstarvmortality_carbonflux_si_pft  => this%hvars(ih_cstarvmortality_carbonflux_si_pft)%r82d, &
         hio_hydraulicmortality_carbonflux_si_pft  => this%hvars(ih_hydraulicmortality_carbonflux_si_pft)%r82d, &
         hio_firemortality_carbonflux_si_pft  => this%hvars(ih_firemortality_carbonflux_si_pft)%r82d, &
         hio_crownarea_si_pft    => this%hvars(ih_crownarea_si_pft)%r82d, &
         hio_canopycrownarea_si_pft  => this%hvars(ih_canopycrownarea_si_pft)%r82d, &
         hio_gpp_si_pft  => this%hvars(ih_gpp_si_pft)%r82d, &
         hio_npp_si_pft  => this%hvars(ih_npp_si_pft)%r82d, &
         hio_fragmentation_scaler_sl  => this%hvars(ih_fragmentation_scaler_sl)%r82d,  &
         hio_litter_in_elem      => this%hvars(ih_litter_in_elem)%r82d, &
         hio_litter_out_elem     => this%hvars(ih_litter_out_elem)%r82d, &
         hio_seed_bank_elem      => this%hvars(ih_seed_bank_elem)%r82d, &
         hio_seeds_in_local_elem => this%hvars(ih_seeds_in_local_elem)%r82d, &
         hio_seed_in_extern_elem => this%hvars(ih_seeds_in_extern_elem)%r82d, &
         hio_seed_decay_elem     => this%hvars(ih_seed_decay_elem)%r82d, &
         hio_seed_germ_elem      => this%hvars(ih_seed_germ_elem)%r82d, &
         hio_gpp_si_scpf         => this%hvars(ih_gpp_si_scpf)%r82d, &
         hio_npp_totl_si_scpf    => this%hvars(ih_npp_totl_si_scpf)%r82d, &
         hio_npp_leaf_si_scpf    => this%hvars(ih_npp_leaf_si_scpf)%r82d, &
         hio_npp_seed_si_scpf    => this%hvars(ih_npp_seed_si_scpf)%r82d, &
         hio_npp_fnrt_si_scpf    => this%hvars(ih_npp_fnrt_si_scpf)%r82d, &
         hio_npp_bgsw_si_scpf    => this%hvars(ih_npp_bgsw_si_scpf)%r82d, &
         hio_npp_bgdw_si_scpf    => this%hvars(ih_npp_bgdw_si_scpf)%r82d, &
         hio_npp_agsw_si_scpf    => this%hvars(ih_npp_agsw_si_scpf)%r82d, &
         hio_npp_agdw_si_scpf    => this%hvars(ih_npp_agdw_si_scpf)%r82d, &
         hio_npp_stor_si_scpf    => this%hvars(ih_npp_stor_si_scpf)%r82d, &
         hio_bstor_canopy_si_scpf      => this%hvars(ih_bstor_canopy_si_scpf)%r82d, &
         hio_bstor_understory_si_scpf  => this%hvars(ih_bstor_understory_si_scpf)%r82d, &
         hio_bleaf_canopy_si_scpf      => this%hvars(ih_bleaf_canopy_si_scpf)%r82d, &
         hio_bleaf_understory_si_scpf  => this%hvars(ih_bleaf_understory_si_scpf)%r82d, &
         hio_lai_canopy_si_scpf        => this%hvars(ih_lai_canopy_si_scpf)%r82d, &
         hio_lai_understory_si_scpf    => this%hvars(ih_lai_understory_si_scpf)%r82d, &
         hio_crownarea_canopy_si_scpf     => this%hvars(ih_crownarea_canopy_si_scpf)%r82d, &
         hio_crownarea_understory_si_scpf => this%hvars(ih_crownarea_understory_si_scpf)%r82d, &
         hio_mortality_canopy_si_scpf         => this%hvars(ih_mortality_canopy_si_scpf)%r82d, &
         hio_mortality_canopy_secondary_si_scls      => this%hvars(ih_mortality_canopy_secondary_si_scls)%r82d, &
         hio_mortality_understory_si_scpf     => this%hvars(ih_mortality_understory_si_scpf)%r82d, &
         hio_m3_mortality_canopy_si_scpf      => this%hvars(ih_m3_mortality_canopy_si_scpf)%r82d, &
         hio_m3_mortality_understory_si_scpf  => this%hvars(ih_m3_mortality_understory_si_scpf)%r82d, &
         hio_m3_mortality_canopy_si_scls    => this%hvars(ih_m3_mortality_canopy_si_scls)%r82d, &
         hio_m3_mortality_understory_si_scls => this%hvars(ih_m3_mortality_understory_si_scls)%r82d, &
         hio_nplant_canopy_si_scpf     => this%hvars(ih_nplant_canopy_si_scpf)%r82d, &
         hio_nplant_understory_si_scpf => this%hvars(ih_nplant_understory_si_scpf)%r82d, &
         hio_ddbh_canopy_si_scpf       => this%hvars(ih_ddbh_canopy_si_scpf)%r82d, &
         hio_ddbh_understory_si_scpf   => this%hvars(ih_ddbh_understory_si_scpf)%r82d, &
         hio_ddbh_canopy_si_scls       => this%hvars(ih_ddbh_canopy_si_scls)%r82d, &
         hio_ddbh_understory_si_scls   => this%hvars(ih_ddbh_understory_si_scls)%r82d, &
         hio_gpp_canopy_si_scpf        => this%hvars(ih_gpp_canopy_si_scpf)%r82d, &
         hio_gpp_understory_si_scpf    => this%hvars(ih_gpp_understory_si_scpf)%r82d, &
         hio_ar_canopy_si_scpf         => this%hvars(ih_ar_canopy_si_scpf)%r82d, &
         hio_ar_understory_si_scpf     => this%hvars(ih_ar_understory_si_scpf)%r82d, &
         hio_ddbh_si_scpf        => this%hvars(ih_ddbh_si_scpf)%r82d, &
         hio_growthflux_si_scpf        => this%hvars(ih_growthflux_si_scpf)%r82d, &
         hio_growthflux_fusion_si_scpf        => this%hvars(ih_growthflux_fusion_si_scpf)%r82d, &
         hio_ba_si_scpf          => this%hvars(ih_ba_si_scpf)%r82d, &
         hio_agb_si_scpf         => this%hvars(ih_agb_si_scpf)%r82d, &
         hio_nplant_si_scpf      => this%hvars(ih_nplant_si_scpf)%r82d, &
         hio_nplant_si_capf      => this%hvars(ih_nplant_si_capf)%r82d, &
         hio_m1_si_scpf          => this%hvars(ih_m1_si_scpf)%r82d, &
         hio_m2_si_scpf          => this%hvars(ih_m2_si_scpf)%r82d, &
         hio_m3_si_scpf          => this%hvars(ih_m3_si_scpf)%r82d, &
         hio_m4_si_scpf          => this%hvars(ih_m4_si_scpf)%r82d, &
         hio_m5_si_scpf          => this%hvars(ih_m5_si_scpf)%r82d, &
         hio_m6_si_scpf          => this%hvars(ih_m6_si_scpf)%r82d, &
         hio_m7_si_scpf          => this%hvars(ih_m7_si_scpf)%r82d, &
         hio_m8_si_scpf          => this%hvars(ih_m8_si_scpf)%r82d, &
         hio_m9_si_scpf          => this%hvars(ih_m9_si_scpf)%r82d, &
         hio_m10_si_scpf         => this%hvars(ih_m10_si_scpf)%r82d, &
         hio_m10_si_capf         => this%hvars(ih_m10_si_capf)%r82d, &
         hio_crownfiremort_si_scpf     => this%hvars(ih_crownfiremort_si_scpf)%r82d, &
         hio_cambialfiremort_si_scpf   => this%hvars(ih_cambialfiremort_si_scpf)%r82d, &
         hio_abg_mortality_cflux_si_scpf    => this%hvars(ih_abg_mortality_cflux_si_scpf)%r82d, &
         hio_abg_productivity_cflux_si_scpf => this%hvars(ih_abg_productivity_cflux_si_scpf)%r82d, &
         hio_burn_flux_elem    => this%hvars(ih_burn_flux_elem)%r82d, &
         hio_m1_si_scls          => this%hvars(ih_m1_si_scls)%r82d, &
         hio_m2_si_scls          => this%hvars(ih_m2_si_scls)%r82d, &
         hio_m3_si_scls          => this%hvars(ih_m3_si_scls)%r82d, &
         hio_m4_si_scls          => this%hvars(ih_m4_si_scls)%r82d, &
         hio_m5_si_scls          => this%hvars(ih_m5_si_scls)%r82d, &
         hio_m6_si_scls          => this%hvars(ih_m6_si_scls)%r82d, &
         hio_m7_si_scls          => this%hvars(ih_m7_si_scls)%r82d, &
         hio_m8_si_scls          => this%hvars(ih_m8_si_scls)%r82d, &
         hio_m9_si_scls          => this%hvars(ih_m9_si_scls)%r82d, &
         hio_m10_si_scls         => this%hvars(ih_m10_si_scls)%r82d, &
         hio_m10_si_cacls        => this%hvars(ih_m10_si_cacls)%r82d)

      ! Break up associates for NAG compilers
      associate(hio_c13disc_si_scpf     => this%hvars(ih_c13disc_si_scpf)%r82d, &
           hio_cwd_elcwd           => this%hvars(ih_cwd_elcwd)%r82d, &
           hio_cwd_ag_elem         => this%hvars(ih_cwd_ag_elem)%r82d, &
           hio_cwd_bg_elem         => this%hvars(ih_cwd_bg_elem)%r82d, &
           hio_fines_ag_elem       => this%hvars(ih_fines_ag_elem)%r82d, &
           hio_fines_bg_elem       => this%hvars(ih_fines_bg_elem)%r82d, &
           hio_ba_si_scls          => this%hvars(ih_ba_si_scls)%r82d, &
           hio_agb_si_scls          => this%hvars(ih_agb_si_scls)%r82d, &
           hio_biomass_si_scls          => this%hvars(ih_biomass_si_scls)%r82d, &
           hio_nplant_si_scls         => this%hvars(ih_nplant_si_scls)%r82d, &
           hio_nplant_si_cacls        => this%hvars(ih_nplant_si_cacls)%r82d, &
           hio_nplant_canopy_si_scls         => this%hvars(ih_nplant_canopy_si_scls)%r82d, &
           hio_nplant_understory_si_scls     => this%hvars(ih_nplant_understory_si_scls)%r82d, &
           hio_lai_canopy_si_scls         => this%hvars(ih_lai_canopy_si_scls)%r82d, &
           hio_lai_understory_si_scls     => this%hvars(ih_lai_understory_si_scls)%r82d, &
           hio_sai_canopy_si_scls         => this%hvars(ih_sai_canopy_si_scls)%r82d, &
           hio_sai_understory_si_scls     => this%hvars(ih_sai_understory_si_scls)%r82d, &
           hio_mortality_canopy_si_scls      => this%hvars(ih_mortality_canopy_si_scls)%r82d, &
           hio_mortality_understory_si_scls  => this%hvars(ih_mortality_understory_si_scls)%r82d, &
           hio_demotion_rate_si_scls         => this%hvars(ih_demotion_rate_si_scls)%r82d, &
           hio_promotion_rate_si_scls        => this%hvars(ih_promotion_rate_si_scls)%r82d, &
           hio_trimming_canopy_si_scls         => this%hvars(ih_trimming_canopy_si_scls)%r82d, &
           hio_trimming_understory_si_scls     => this%hvars(ih_trimming_understory_si_scls)%r82d, &
           hio_crown_area_canopy_si_scls         => this%hvars(ih_crown_area_canopy_si_scls)%r82d, &
           hio_crown_area_understory_si_scls     => this%hvars(ih_crown_area_understory_si_scls)%r82d, &
           hio_leaf_md_canopy_si_scls           => this%hvars(ih_leaf_md_canopy_si_scls)%r82d, &
           hio_root_md_canopy_si_scls           => this%hvars(ih_root_md_canopy_si_scls)%r82d, &
           hio_carbon_balance_canopy_si_scls    => this%hvars(ih_carbon_balance_canopy_si_scls)%r82d, &
           hio_bsw_md_canopy_si_scls            => this%hvars(ih_bsw_md_canopy_si_scls)%r82d, &
           hio_bdead_md_canopy_si_scls          => this%hvars(ih_bdead_md_canopy_si_scls)%r82d, &
           hio_bstore_md_canopy_si_scls         => this%hvars(ih_bstore_md_canopy_si_scls)%r82d, &
           hio_seed_prod_canopy_si_scls         => this%hvars(ih_seed_prod_canopy_si_scls)%r82d, &
           hio_npp_leaf_canopy_si_scls          => this%hvars(ih_npp_leaf_canopy_si_scls)%r82d, &
           hio_npp_fnrt_canopy_si_scls         => this%hvars(ih_npp_fnrt_canopy_si_scls)%r82d, &
           hio_npp_sapw_canopy_si_scls           => this%hvars(ih_npp_sapw_canopy_si_scls)%r82d, &
           hio_npp_dead_canopy_si_scls         => this%hvars(ih_npp_dead_canopy_si_scls)%r82d, &
           hio_npp_seed_canopy_si_scls         => this%hvars(ih_npp_seed_canopy_si_scls)%r82d, &
           hio_npp_stor_canopy_si_scls         => this%hvars(ih_npp_stor_canopy_si_scls)%r82d, &
           hio_leaf_md_understory_si_scls       => this%hvars(ih_leaf_md_understory_si_scls)%r82d, &
           hio_root_md_understory_si_scls       => this%hvars(ih_root_md_understory_si_scls)%r82d, &
           hio_carbon_balance_understory_si_scls=> this%hvars(ih_carbon_balance_understory_si_scls)%r82d, &
           hio_bstore_md_understory_si_scls     => this%hvars(ih_bstore_md_understory_si_scls)%r82d, &
           hio_bsw_md_understory_si_scls        => this%hvars(ih_bsw_md_understory_si_scls)%r82d, &
           hio_bdead_md_understory_si_scls      => this%hvars(ih_bdead_md_understory_si_scls)%r82d, &
           hio_seed_prod_understory_si_scls     => this%hvars(ih_seed_prod_understory_si_scls)%r82d, &
           hio_npp_leaf_understory_si_scls      => this%hvars(ih_npp_leaf_understory_si_scls)%r82d, &
           hio_npp_fnrt_understory_si_scls     => this%hvars(ih_npp_fnrt_understory_si_scls)%r82d, &
           hio_npp_sapw_understory_si_scls       => this%hvars(ih_npp_sapw_understory_si_scls)%r82d, &
           hio_npp_dead_understory_si_scls     => this%hvars(ih_npp_dead_understory_si_scls)%r82d, &
           hio_npp_seed_understory_si_scls     => this%hvars(ih_npp_seed_understory_si_scls)%r82d, &
           hio_npp_stor_understory_si_scls     => this%hvars(ih_npp_stor_understory_si_scls)%r82d, &
           hio_nplant_si_scagpft                => this%hvars(ih_nplant_si_scagpft)%r82d, &
           hio_npp_si_agepft                    => this%hvars(ih_npp_si_agepft)%r82d, &
           hio_biomass_si_agepft                => this%hvars(ih_biomass_si_agepft)%r82d, &
           hio_scorch_height_si_agepft          => this%hvars(ih_scorch_height_si_agepft)%r82d, &
           hio_yesterdaycanopylevel_canopy_si_scls     => this%hvars(ih_yesterdaycanopylevel_canopy_si_scls)%r82d, &
           hio_yesterdaycanopylevel_understory_si_scls => this%hvars(ih_yesterdaycanopylevel_understory_si_scls)%r82d, &
           hio_area_si_age         => this%hvars(ih_area_si_age)%r82d, &
           hio_lai_si_age          => this%hvars(ih_lai_si_age)%r82d, &
           hio_canopy_area_si_age  => this%hvars(ih_canopy_area_si_age)%r82d, &
           hio_ncl_si_age          => this%hvars(ih_ncl_si_age)%r82d, &
           hio_npatches_si_age     => this%hvars(ih_npatches_si_age)%r82d, &
           hio_zstar_si_age        => this%hvars(ih_zstar_si_age)%r82d, &
           hio_biomass_si_age        => this%hvars(ih_biomass_si_age)%r82d, &
           hio_npp_si_age                     => this%hvars(ih_npp_si_age)%r82d, &
           hio_npp_si_landuse                 => this%hvars(ih_npp_si_landuse)%r82d, &
           hio_agesince_anthrodist_si_age     => this%hvars(ih_agesince_anthrodist_si_age)%r82d, &
           hio_secondarylands_area_si_age    => this%hvars(ih_secondarylands_area_si_age)%r82d, &
           hio_primarylands_area_si_age      => this%hvars(ih_primarylands_area_si_age)%r82d, &
           hio_area_si_landuse               => this%hvars(ih_area_si_landuse)%r82d, &
           hio_biomass_si_landuse            => this%hvars(ih_biomass_si_landuse)%r82d, &
           hio_burnedarea_si_landuse         => this%hvars(ih_burnedarea_si_landuse)%r82d, &
           hio_area_burnt_si_age              => this%hvars(ih_area_burnt_si_age)%r82d, &
                                ! hio_fire_rate_of_spread_front_si_age  => this%hvars(ih_fire_rate_of_spread_front_si_age)%r82d, &
           hio_fire_intensity_si_age          => this%hvars(ih_fire_intensity_si_age)%r82d, &
           hio_fire_sum_fuel_si_age           => this%hvars(ih_fire_sum_fuel_si_age)%r82d, &
           hio_burnt_frac_litter_si_fuel      => this%hvars(ih_burnt_frac_litter_si_fuel)%r82d, &
           hio_fuel_amount_si_fuel            => this%hvars(ih_fuel_amount_si_fuel)%r82d, &
           hio_fuel_amount_age_fuel            => this%hvars(ih_fuel_amount_age_fuel)%r82d, &
           hio_canopy_height_dist_si_height   => this%hvars(ih_canopy_height_dist_si_height)%r82d, &
           hio_leaf_height_dist_si_height     => this%hvars(ih_leaf_height_dist_si_height)%r82d, &
           hio_litter_moisture_si_fuel        => this%hvars(ih_litter_moisture_si_fuel)%r82d, &
           hio_cwd_ag_si_cwdsc                  => this%hvars(ih_cwd_ag_si_cwdsc)%r82d, &
           hio_cwd_bg_si_cwdsc                  => this%hvars(ih_cwd_bg_si_cwdsc)%r82d, &
           hio_cwd_ag_in_si_cwdsc               => this%hvars(ih_cwd_ag_in_si_cwdsc)%r82d, &
           hio_cwd_bg_in_si_cwdsc               => this%hvars(ih_cwd_bg_in_si_cwdsc)%r82d, &
           hio_cwd_ag_out_si_cwdsc              => this%hvars(ih_cwd_ag_out_si_cwdsc)%r82d, &
           hio_cwd_bg_out_si_cwdsc              => this%hvars(ih_cwd_bg_out_si_cwdsc)%r82d, &
           hio_crownarea_si_cnlf                => this%hvars(ih_crownarea_si_cnlf)%r82d, &
           hio_crownarea_cl                 => this%hvars(ih_crownarea_cl)%r82d, &
           hio_nplant_si_scag                   => this%hvars(ih_nplant_si_scag)%r82d, &
           hio_nplant_canopy_si_scag            => this%hvars(ih_nplant_canopy_si_scag)%r82d, &
           hio_nplant_understory_si_scag        => this%hvars(ih_nplant_understory_si_scag)%r82d, &
           hio_ddbh_canopy_si_scag              => this%hvars(ih_ddbh_canopy_si_scag)%r82d, &
           hio_ddbh_understory_si_scag          => this%hvars(ih_ddbh_understory_si_scag)%r82d, &
           hio_mortality_canopy_si_scag         => this%hvars(ih_mortality_canopy_si_scag)%r82d, &
           hio_mortality_understory_si_scag     => this%hvars(ih_mortality_understory_si_scag)%r82d )

        ! Break up associates for NAG compilers
        associate( hio_site_dstatus_si_pft              => this%hvars(ih_site_dstatus_si_pft)%r82d, &
             hio_dleafoff_si_pft                  => this%hvars(ih_dleafoff_si_pft)%r82d, &
             hio_dleafon_si_pft                   => this%hvars(ih_dleafon_si_pft)%r82d, &
             hio_meanliqvol_si_pft                => this%hvars(ih_meanliqvol_si_pft)%r82d, &
             hio_meansmp_si_pft                   => this%hvars(ih_meansmp_si_pft)%r82d, &
             hio_elong_factor_si_pft              => this%hvars(ih_elong_factor_si_pft)%r82d, &
             hio_nplant_si_scag                   => this%hvars(ih_nplant_si_scag)%r82d, &
             hio_nplant_canopy_si_scag            => this%hvars(ih_nplant_canopy_si_scag)%r82d, &
             hio_nplant_understory_si_scag        => this%hvars(ih_nplant_understory_si_scag)%r82d, &
             hio_disturbance_rate_si_lulu         => this%hvars(ih_disturbance_rate_si_lulu)%r82d, &
             hio_cstarvmortality_continuous_carbonflux_si_pft  => this%hvars(ih_cstarvmortality_continuous_carbonflux_si_pft)%r82d, &
             hio_transition_matrix_si_lulu      => this%hvars(ih_transition_matrix_si_lulu)%r82d, &
             hio_sapwood_area_scpf              => this%hvars(ih_sapwood_area_scpf)%r82d)

          model_day_int = nint(hlm_model_day)

          ! ---------------------------------------------------------------------------------
          ! Loop through the FATES scale hierarchy and fill the history IO arrays
          ! ---------------------------------------------------------------------------------

          
          siteloop: do s = 1,nsites

             io_si  = sites(s)%h_gid

             ! C13 will not get b4b restarts on the first day because
             ! there is no mechanism to remember the previous day's values
             ! through a restart. This should be added with the next refactor
             gpp_cached_scpf(:) = hio_gpp_si_scpf(io_si,:)
             
             call this%zero_site_hvars(sites(s),upfreq_in=group_dyna_complx)
             
             ! These are weighting factors
             storen_canopy_scpf(:) = 0._r8
             storen_understory_scpf(:) = 0._r8
             storep_canopy_scpf(:) = 0._r8
             storep_understory_scpf(:) = 0._r8
             storec_canopy_scpf(:) = 0._r8
             storec_understory_scpf(:) = 0._r8
             
             do i_dist = 1, n_landuse_cats
                do j_dist = 1, n_landuse_cats
                   
                   ! roll up disturbance rates in land-use x land-use array into a single dimension
                   hio_disturbance_rate_si_lulu(io_si, i_dist+n_landuse_cats*(j_dist-1)) = &
                        sum(sites(s)%disturbance_rates(1:n_dist_types,i_dist, j_dist)) * &
                        days_per_year

                   ! get the land use transition matrix and output that to history.
                   ! (mainly a sanity check, can maybe remove before integration)
                   hio_transition_matrix_si_lulu(io_si, i_dist+n_landuse_cats*(j_dist-1)) = &
                        sites(s)%landuse_transition_matrix(i_dist, j_dist)
                end do
             end do
             
             do el = 1, num_elements

                if((hlm_use_ed_st3 .eq. ifalse) .and. (hlm_use_sp .eq. ifalse)) then

                   ! Total model error [kg/day -> kg/s]  (all elements)
                   this%hvars(ih_err_fates_elem)%r82d(io_si,el) = sites(s)%mass_balance(el)%err_fates / sec_per_day
                   
                   this%hvars(ih_interr_liveveg_elem)%r82d(io_si,el) =  sites(s)%flux_diags%elem(el)%err_liveveg

                   this%hvars(ih_interr_litter_elem)%r82d(io_si,el) = sites(s)%flux_diags%elem(el)%err_litter
                   
                end if
                   
                ! Total element lost to atmosphere from burning (kg/site/day -> kg/m2/s)
                hio_burn_flux_elem(io_si,el) = &
                     sites(s)%mass_balance(el)%burn_flux_to_atm * ha_per_m2 *           &
                     days_per_sec

             end do


             ! Update drought deciduous information (now separated by PFT).
             do ft = 1,numpft
                ! Update the site-PFT status for drought deciduous
                hio_site_dstatus_si_pft(io_si,ft) = real(sites(s)%dstatus(ft),r8)

                ! Model days elapsed since leaf off/on for drought deciduous
                hio_dleafoff_si_pft(io_si,ft)     = real(sites(s)%dndaysleafon (ft),r8)
                hio_dleafon_si_pft(io_si,ft)      = real(sites(s)%dndaysleafoff(ft),r8)

                ! Leaf elongation factor (0 means fully abscissed, 1 means fully flushed).
                hio_elong_factor_si_pft(io_si,ft) = sites(s)%elong_factor(ft)

                if(model_day_int>numWaterMem)then
                   ! Mean liquid water content (m3/m3) used for drought phenology
                   hio_meanliqvol_si_pft(io_si,ft) = &
                        sum(sites(s)%liqvol_memory(1:numWaterMem,ft))/real(numWaterMem,r8)

                   ! Mean soil matric potential (Pa) used for drought phenology
                   hio_meansmp_si_pft(io_si,ft) = &
                        sum(sites(s)%smp_memory(1:numWaterMem,ft))/real(numWaterMem,r8) &
                        * dens_fresh_liquid_water * grav_earth * m_per_mm
                end if
             end do

             ! Loop through patches to sum up diagonistics
             cpatch => sites(s)%oldest_patch
             patchloop: do while(associated(cpatch))


                cpatch%age_class  = get_age_class_index(cpatch%age)

                ! Increment the fractional area in each age class bin
                hio_area_si_age(io_si,cpatch%age_class) = hio_area_si_age(io_si,cpatch%age_class) &
                     + cpatch%area * AREA_INV

                ! ignore land use info on nocomp bareground (where landuse label = 0)
                if (cpatch%land_use_label .gt. nocomp_bareground_land) then 
                   hio_area_si_landuse(io_si, cpatch%land_use_label) = &
                        hio_area_si_landuse(io_si, cpatch%land_use_label) &
                        + cpatch%area * AREA_INV

                   hio_burnedarea_si_landuse(io_si, cpatch%land_use_label) = &
                        hio_burnedarea_si_landuse(io_si, cpatch%land_use_label) + &
                        cpatch%frac_burnt * cpatch%area * AREA_INV / sec_per_day
                end if
                   
                ! Increment some patch-age-resolved diagnostics
                hio_lai_si_age(io_si,cpatch%age_class) = hio_lai_si_age(io_si,cpatch%age_class) &
                     + sum(cpatch%tlai_profile(:,:,:) * cpatch%canopy_area_profile(:,:,:) ) * cpatch%total_canopy_area
                
                hio_ncl_si_age(io_si,cpatch%age_class) = hio_ncl_si_age(io_si,cpatch%age_class) &
                     + cpatch%ncl_p * cpatch%area
                
                hio_npatches_si_age(io_si,cpatch%age_class) = hio_npatches_si_age(io_si,cpatch%age_class) + 1._r8



                if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
                   hio_zstar_si_age(io_si,cpatch%age_class) = hio_zstar_si_age(io_si,cpatch%age_class) &
                        + cpatch%zstar * cpatch%area * AREA_INV
                endif

                ! some diagnostics on secondary forest area and its age distribution
                if ( cpatch%land_use_label .eq. secondaryland ) then

                   ageclass_since_anthrodist = get_age_class_index(cpatch%age_since_anthro_disturbance)

                   hio_agesince_anthrodist_si_age(io_si,ageclass_since_anthrodist) = &
                        hio_agesince_anthrodist_si_age(io_si,ageclass_since_anthrodist)  &
                        + cpatch%area * AREA_INV

                   hio_secondarylands_area_si_age(io_si,cpatch%age_class) = &
                        hio_secondarylands_area_si_age(io_si,cpatch%age_class) & 
                        + cpatch%area * AREA_INV

                else if ( cpatch%land_use_label .eq. primaryland) then
                   hio_primarylands_area_si_age(io_si,cpatch%age_class) = &
                        hio_primarylands_area_si_age(io_si,cpatch%age_class) & 
                        + cpatch%area * AREA_INV

                endif



                ! patch-age-resolved fire variables
                do ft = 1,numpft
                   ! for scorch height, weight the value by patch area within any
                   ! given age class - in the event that there is more than one
                   ! patch per age class.
                   iagepft = cpatch%age_class + (ft-1) * nlevage
                   hio_scorch_height_si_agepft(io_si,iagepft) = hio_scorch_height_si_agepft(io_si,iagepft) + &
                        cpatch%Scorch_ht(ft) * cpatch%area

                   ! and also pft-labeled patch areas in the event that we are in nocomp mode
                   if ( hlm_use_nocomp .eq. itrue .and. cpatch%nocomp_pft_label .eq. ft) then 
                      this%hvars(ih_nocomp_pftpatchfraction_si_pft)%r82d(io_si,ft) = &
                           this%hvars(ih_nocomp_pftpatchfraction_si_pft)%r82d(io_si,ft) + cpatch%area * AREA_INV

                      this%hvars(ih_nocomp_pftnpatches_si_pft)%r82d(io_si,ft) = &
                           this%hvars(ih_nocomp_pftnpatches_si_pft)%r82d(io_si,ft) + 1._r8

                      this%hvars(ih_nocomp_pftburnedarea_si_pft)%r82d(io_si,ft) = &
                           this%hvars(ih_nocomp_pftburnedarea_si_pft)%r82d(io_si,ft) + &
                           cpatch%frac_burnt * cpatch%area * AREA_INV / sec_per_day
                   endif

                end do

                ! fractional area burnt [frac/day] -> [frac/sec]
                hio_area_burnt_si_age(io_si,cpatch%age_class) = hio_area_burnt_si_age(io_si,cpatch%age_class) + &
                     cpatch%frac_burnt * cpatch%area * AREA_INV / sec_per_day

                ! hio_fire_rate_of_spread_front_si_age(io_si, cpatch%age_class) = hio_fire_rate_of_spread_si_age(io_si, cpatch%age_class) + &
                !    cpatch%ros_front * cpatch*frac_burnt * cpatch%area * AREA_INV

                ! Fire intensity weighted by burned fraction [kJ/m/s] -> [J/m/s]
                hio_fire_intensity_si_age(io_si, cpatch%age_class) = hio_fire_intensity_si_age(io_si, cpatch%age_class) + &
                     cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV * J_per_kJ

                ! Fuel sum [kg/m2]
                hio_fire_sum_fuel_si_age(io_si, cpatch%age_class) = hio_fire_sum_fuel_si_age(io_si, cpatch%age_class) +  &
                     cpatch%fuel%non_trunk_loading * cpatch%area * AREA_INV



                ! loop through cohorts on patch
                ccohort => cpatch%shortest
                cohortloop: do while(associated(ccohort))

                   ft = ccohort%pft

                   ! get indices for size class x pft and cohort age x pft
                   ! size class is the fastest changing dimension
                   call sizetype_class_index(ccohort%dbh, ccohort%pft,                &
                        ccohort%size_class, ccohort%size_by_pft_class)
                   ! cohort age is the fastest changing dimension
                   call coagetype_class_index(ccohort%coage, ccohort%pft,             &
                        ccohort%coage_class, ccohort%coage_by_pft_class)

                   n_perm2 = ccohort%n * AREA_INV

                   hio_canopy_area_si_age(io_si,cpatch%age_class) = hio_canopy_area_si_age(io_si,cpatch%age_class) &
                        + ccohort%c_area * AREA_INV

                   ! calculate leaf height distribution, assuming leaf area is evenly distributed thru crown depth
                   call CrownDepth(ccohort%height,ft,crown_depth)
                   height_bin_max = get_height_index(ccohort%height)
                   height_bin_min = get_height_index(ccohort%height - crown_depth)
                   do i_heightbin = height_bin_min, height_bin_max
                      binbottom = ED_val_history_height_bin_edges(i_heightbin)
                      if (i_heightbin .eq. nlevheight) then
                         bintop = reallytalltrees
                      else
                         bintop = ED_val_history_height_bin_edges(i_heightbin+1)
                      endif
                      ! what fraction of a cohort's crown is in this height bin?
                      frac_canopy_in_bin = (min(bintop,ccohort%height) - &
                           max(binbottom,ccohort%height-crown_depth)) / &
                           (crown_depth)

                      hio_leaf_height_dist_si_height(io_si,i_heightbin) = &
                           hio_leaf_height_dist_si_height(io_si,i_heightbin) + &
                           ccohort%c_area * AREA_INV * ccohort%treelai * frac_canopy_in_bin

                      ! if ( ( ccohort%c_area * AREA_INV * ccohort%treelai * frac_canopy_in_bin) .lt. 0._r8) then
                      !    write(fates_log(),*) ' negative hio_leaf_height_dist_si_height:'
                      !    write(fates_log(),*) '   c_area, treelai, frac_canopy_in_bin:', ccohort%c_area, ccohort%treelai, frac_canopy_in_bin
                      ! endif
                   end do

                   if (ccohort%canopy_layer .eq. 1) then
                      ! calculate the area of canopy that is within each height bin
                      hio_canopy_height_dist_si_height(io_si,height_bin_max) = &
                           hio_canopy_height_dist_si_height(io_si,height_bin_max) + ccohort%c_area * AREA_INV
                   endif

                   call set_root_fraction(sites(s)%rootfrac_scr, ccohort%pft, sites(s)%zi_soil, &
                        bc_in(s)%max_rooting_depth_index_col )

                   ! Update biomass components
                   ! Mass pools [kg]
                   elloop: do el = 1, num_elements

                      sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
                      struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
                      leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
                      fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
                      store_m  = ccohort%prt%GetState(store_organ, element_list(el))
                      repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
                      alive_m  = leaf_m + fnrt_m + sapw_m
                      total_m  = alive_m + store_m + struct_m

                      i_scpf = ccohort%size_by_pft_class


                      ! Plant multi-element states and fluxes
                      ! Zero states, and set the fluxes
                      if( element_list(el).eq.carbon12_element )then

                         call bstore_allom(ccohort%dbh,ccohort%pft,ccohort%crowndamage,ccohort%canopy_trim, store_max)
                         
                         ! Determine the root carbon biomass in kg/m3
                         ! [kg/m3] = [kg/plant] * [plant/ha] / [m3/ha] * [fraction] / [m]

                         do ilyr = 1,sites(s)%nlevsoil
                            this%hvars(ih_fnrtc_sl)%r82d(io_si,ilyr) = this%hvars(ih_fnrtc_sl)%r82d(io_si,ilyr) + &
                                 fnrt_m * ccohort%n / area * sites(s)%rootfrac_scr(ilyr) / sites(s)%dz_soil(ilyr)
                         end do

                         ! Update PFT partitioned biomass components
                         hio_leafbiomass_si_pft(io_si,ft) = hio_leafbiomass_si_pft(io_si,ft) + &
                              (ccohort%n * AREA_INV) * leaf_m

                         hio_storebiomass_si_pft(io_si,ft) = hio_storebiomass_si_pft(io_si,ft) + &
                              (ccohort%n * AREA_INV) * store_m

                         hio_nindivs_si_pft(io_si,ft) = hio_nindivs_si_pft(io_si,ft) + &
                              ccohort%n * AREA_INV

                         if(ccohort%isnew) then
                            hio_recruitment_cflux_si_pft(io_si, ft) = hio_recruitment_cflux_si_pft(io_si, ft) + &
                                 (ccohort%n * AREA_INV) * total_m * days_per_year
                         end if

                         hio_biomass_si_pft(io_si, ft) = hio_biomass_si_pft(io_si, ft) + &
                              (ccohort%n * AREA_INV) * total_m

                         ! update total biomass per age bin
                         hio_biomass_si_age(io_si,cpatch%age_class) = hio_biomass_si_age(io_si,cpatch%age_class) &
                              + total_m * ccohort%n * AREA_INV

                         ! biomass by land use type
                         hio_biomass_si_landuse(io_si, cpatch%land_use_label) = &
                              hio_biomass_si_landuse(io_si, cpatch%land_use_label) &
                              + total_m * ccohort%n * AREA_INV

                         if (ccohort%canopy_layer .eq. 1) then
                            storec_canopy_scpf(i_scpf) = &
                                 storec_canopy_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storectfrac_canopy_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storectfrac_canopy_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max
                         else
                            storec_understory_scpf(i_scpf) = &
                                 storec_understory_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storectfrac_ustory_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storectfrac_ustory_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max  
                         end if


                      elseif(element_list(el).eq.nitrogen_element)then

                         store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                         if (ccohort%canopy_layer .eq. 1) then
                            storen_canopy_scpf(i_scpf) = &
                                 storen_canopy_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storentfrac_canopy_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storentfrac_canopy_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max
                         else
                            storen_understory_scpf(i_scpf) = &
                                 storen_understory_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storentfrac_understory_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storentfrac_understory_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max  
                         end if

                      elseif(element_list(el).eq.phosphorus_element) then

                         store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                         if (ccohort%canopy_layer .eq. 1) then
                            storep_canopy_scpf(i_scpf) = &
                                 storep_canopy_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storeptfrac_canopy_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storeptfrac_canopy_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max
                         else
                            storep_understory_scpf(i_scpf) = &
                                 storep_understory_scpf(i_scpf) + ccohort%n * store_m
                            this%hvars(ih_storeptfrac_understory_scpf)%r82d(io_si,i_scpf) = & 
                                 this%hvars(ih_storeptfrac_understory_scpf)%r82d(io_si,i_scpf) + &
                                 ccohort%n * store_max
                         end if


                      end if
                   end do elloop

                   ! Update PFT crown area
                   hio_crownarea_si_pft(io_si, ft) = hio_crownarea_si_pft(io_si, ft) + &
                        ccohort%c_area * AREA_INV

                   if (ccohort%canopy_layer .eq. 1) then
                      ! Update PFT canopy crown area
                      hio_canopycrownarea_si_pft(io_si, ft) = hio_canopycrownarea_si_pft(io_si, ft) + &
                           ccohort%c_area * AREA_INV
                   end if

                   ! Site by Size-Class x PFT (SCPF)
                   ! ------------------------------------------------------------------------

                   ! Flux Variables (cohorts must had experienced a day before any of these values
                   ! have any meaning, otherwise they are just inialization values
                   notnew: if( .not.(ccohort%isnew) ) then

                      ! update pft-resolved NPP and GPP fluxes
                      hio_gpp_si_pft(io_si, ft) = hio_gpp_si_pft(io_si, ft) + &
                           ccohort%gpp_acc_hold * n_perm2 / (days_per_year* sec_per_day)

                      hio_npp_si_pft(io_si, ft) = hio_npp_si_pft(io_si, ft) + &
                           ccohort%npp_acc_hold * n_perm2 / (days_per_year*sec_per_day)

                      if (cpatch%land_use_label .gt. nocomp_bareground_land) then
                         hio_npp_si_landuse(io_si,cpatch%land_use_label) = hio_npp_si_landuse(io_si,cpatch%land_use_label) &
                              + ccohort%npp_acc_hold * n_perm2 / (days_per_year*sec_per_day)
                      end if

                      ! Turnover pools [kgC/day] * [day/yr] = [kgC/yr]
                      sapw_m_turnover   = ccohort%prt%GetTurnover(sapw_organ, carbon12_element) * days_per_year
                      store_m_turnover  = ccohort%prt%GetTurnover(store_organ, carbon12_element) * days_per_year
                      leaf_m_turnover   = ccohort%prt%GetTurnover(leaf_organ, carbon12_element) * days_per_year
                      fnrt_m_turnover   = ccohort%prt%GetTurnover(fnrt_organ, carbon12_element) * days_per_year
                      struct_m_turnover = ccohort%prt%GetTurnover(struct_organ, carbon12_element) * days_per_year

                      ! Net change from allocation and transport [kgC/day] * [day/yr] = [kgC/yr]
                      sapw_m_net_alloc   = ccohort%prt%GetNetAlloc(sapw_organ, carbon12_element) * days_per_year
                      store_m_net_alloc  = ccohort%prt%GetNetAlloc(store_organ, carbon12_element) * days_per_year
                      leaf_m_net_alloc   = ccohort%prt%GetNetAlloc(leaf_organ, carbon12_element) * days_per_year
                      fnrt_m_net_alloc   = ccohort%prt%GetNetAlloc(fnrt_organ, carbon12_element) * days_per_year
                      struct_m_net_alloc = ccohort%prt%GetNetAlloc(struct_organ, carbon12_element) * days_per_year
                      repro_m_net_alloc  = ccohort%prt%GetNetAlloc(repro_organ, carbon12_element) * days_per_year



                      associate( scpf => ccohort%size_by_pft_class,                   &
                           scls => ccohort%size_class,                          &
                           cacls => ccohort%coage_class,                        &
                           capf => ccohort%coage_by_pft_class,                  &
                           cdam => ccohort%crowndamage)

                        ! convert [kgC/plant/year] -> [kgC/m2/s]
                        hio_gpp_si_scpf(io_si,scpf) = hio_gpp_si_scpf(io_si,scpf) +     &
                             n_perm2*ccohort%gpp_acc_hold / (days_per_year*sec_per_day)
                        
                        hio_npp_totl_si_scpf(io_si,scpf) = hio_npp_totl_si_scpf(io_si,scpf) + &
                             ccohort%npp_acc_hold * n_perm2 / (days_per_year*sec_per_day)

                        hio_npp_leaf_si_scpf(io_si,scpf) = hio_npp_leaf_si_scpf(io_si,scpf) + &
                             leaf_m_net_alloc*n_perm2 / (days_per_year*sec_per_day)
                        
                        hio_npp_fnrt_si_scpf(io_si,scpf) = hio_npp_fnrt_si_scpf(io_si,scpf) + &
                             fnrt_m_net_alloc*n_perm2 / (days_per_year*sec_per_day)
                        
                        hio_npp_bgsw_si_scpf(io_si,scpf) = hio_npp_bgsw_si_scpf(io_si,scpf) + &
                             sapw_m_net_alloc*n_perm2*(1._r8-prt_params%allom_agb_frac(ccohort%pft)) / &
                             (days_per_year*sec_per_day)
                        
                        hio_npp_agsw_si_scpf(io_si,scpf) = hio_npp_agsw_si_scpf(io_si,scpf) + &
                             sapw_m_net_alloc*n_perm2*prt_params%allom_agb_frac(ccohort%pft) / &
                             (days_per_year*sec_per_day)
                        
                        hio_npp_bgdw_si_scpf(io_si,scpf) = hio_npp_bgdw_si_scpf(io_si,scpf) + &
                             struct_m_net_alloc*n_perm2*(1._r8-prt_params%allom_agb_frac(ccohort%pft)) / &
                             (days_per_year*sec_per_day)
                        
                        hio_npp_agdw_si_scpf(io_si,scpf) = hio_npp_agdw_si_scpf(io_si,scpf) + &
                             struct_m_net_alloc*n_perm2*prt_params%allom_agb_frac(ccohort%pft) / &
                             (days_per_year*sec_per_day)
                        
                        hio_npp_seed_si_scpf(io_si,scpf) = hio_npp_seed_si_scpf(io_si,scpf) + &
                             repro_m_net_alloc*n_perm2 / (days_per_year*sec_per_day)
                        
                        hio_npp_stor_si_scpf(io_si,scpf) = hio_npp_stor_si_scpf(io_si,scpf) + &
                             store_m_net_alloc*n_perm2 / (days_per_year*sec_per_day)

                        ! Woody State Variables (basal area growth increment)
                        if ( prt_params%woody(ft) == itrue) then

                           ! basal area  [m2/m2]
                           hio_ba_si_scpf(io_si,scpf) = hio_ba_si_scpf(io_si,scpf) + &
                                0.25_r8*pi_const*((ccohort%dbh/100.0_r8)**2.0_r8)*ccohort%n / m2_per_ha

                           ! also by size class only
                           hio_ba_si_scls(io_si,scls) = hio_ba_si_scls(io_si,scls) + &
                                0.25_r8*pi_const*((ccohort%dbh/100.0_r8)**2.0_r8)*       &
                                ccohort%n / m2_per_ha

                           ! growth increment
                           hio_ddbh_si_scpf(io_si,scpf) = hio_ddbh_si_scpf(io_si,scpf) + &
                                ccohort%ddbhdt*ccohort%n / m2_per_ha * m_per_cm

                           call bsap_allom(ccohort%dbh,ccohort%pft,ccohort%crowndamage,&
                                ccohort%canopy_trim, ccohort%efstem_coh, a_sapw, c_sapw)

                           ! sapwood area [m2/m2]
                           hio_sapwood_area_scpf(io_si,scpf)  = hio_sapwood_area_scpf(io_si,scpf)+ &
                                a_sapw*ccohort%n / m2_per_ha

                        end if

                        ! mortality sums [#/m2]
                        hio_m1_si_scpf(io_si,scpf) = hio_m1_si_scpf(io_si,scpf) +       &
                             ccohort%bmort*ccohort%n / m2_per_ha
                        hio_m2_si_scpf(io_si,scpf) = hio_m2_si_scpf(io_si,scpf) +       &
                             ccohort%hmort*ccohort%n / m2_per_ha
                        hio_m3_si_scpf(io_si,scpf) = hio_m3_si_scpf(io_si,scpf) +       &
                             ccohort%cmort*ccohort%n / m2_per_ha

                        hio_m7_si_scpf(io_si,scpf) = hio_m7_si_scpf(io_si,scpf) +       &
                             (ccohort%lmort_direct + ccohort%lmort_collateral +           &
                             ccohort%lmort_infra) * ccohort%n / m2_per_ha

                        hio_m8_si_scpf(io_si,scpf) = hio_m8_si_scpf(io_si,scpf) +       &
                             ccohort%frmort*ccohort%n / m2_per_ha
                        hio_m9_si_scpf(io_si,scpf) = hio_m9_si_scpf(io_si,scpf) +       &
                             ccohort%smort*ccohort%n / m2_per_ha

                        if (hlm_use_cohort_age_tracking .eq.itrue) then
                           hio_m10_si_scpf(io_si,scpf) = hio_m10_si_scpf(io_si,scpf) +  &
                                ccohort%asmort*ccohort%n / m2_per_ha
                           hio_m10_si_capf(io_si,capf) = hio_m10_si_capf(io_si,capf) +  &
                                ccohort%asmort*ccohort%n / m2_per_ha
                           hio_m10_si_scls(io_si,scls) = hio_m10_si_scls(io_si,scls) +  &
                                ccohort%asmort*ccohort%n / m2_per_ha
                           hio_m10_si_cacls(io_si,cacls) = hio_m10_si_cacls(io_si,cacls)+ &
                                ccohort%asmort*ccohort%n / m2_per_ha
                        end if




                        hio_m1_si_scls(io_si,scls) = hio_m1_si_scls(io_si,scls) + ccohort%bmort*ccohort%n / m2_per_ha
                        hio_m2_si_scls(io_si,scls) = hio_m2_si_scls(io_si,scls) + ccohort%hmort*ccohort%n / m2_per_ha
                        hio_m3_si_scls(io_si,scls) = hio_m3_si_scls(io_si,scls) + ccohort%cmort*ccohort%n / m2_per_ha
                        hio_m7_si_scls(io_si,scls) = hio_m7_si_scls(io_si,scls) + &
                             (ccohort%lmort_direct+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n / m2_per_ha
                        hio_m8_si_scls(io_si,scls) = hio_m8_si_scls(io_si,scls) + &
                             ccohort%frmort*ccohort%n / m2_per_ha
                        hio_m9_si_scls(io_si,scls) = hio_m9_si_scls(io_si,scls) + ccohort%smort*ccohort%n / m2_per_ha

                        !C13 discrimination
                        if(abs(gpp_cached_scpf(scpf)-hlm_hio_ignore_val)>nearzero .and. &
                              (gpp_cached_scpf(scpf) + ccohort%gpp_acc_hold) > 0.0_r8) then
                           
                           gpp_cached = gpp_cached_scpf(scpf)*days_per_year*sec_per_day
                           
                           hio_c13disc_si_scpf(io_si,scpf) = ((hio_c13disc_si_scpf(io_si,scpf) * gpp_cached) + &
                                (ccohort%c13disc_acc * ccohort%gpp_acc_hold)) / (gpp_cached + ccohort%gpp_acc_hold)
                        else
                           hio_c13disc_si_scpf(io_si,scpf) = 0.0_r8
                        end if

                        ! number density [/m2]
                        hio_nplant_si_scpf(io_si,scpf) = hio_nplant_si_scpf(io_si,scpf) + ccohort%n / m2_per_ha

                        ! number density along the cohort age dimension
                        if (hlm_use_cohort_age_tracking .eq.itrue) then
                           hio_nplant_si_capf(io_si,capf) = hio_nplant_si_capf(io_si,capf) + ccohort%n / m2_per_ha
                           hio_nplant_si_cacls(io_si,cacls) = hio_nplant_si_cacls(io_si,cacls)+ccohort%n / m2_per_ha
                        end if

                        ! damage variables - cohort level 
                        if(hlm_use_tree_damage .eq. itrue) then

                           icdpf = get_cdamagesizepft_class_index(ccohort%dbh, ccohort%crowndamage, ccohort%pft)

                           this%hvars(ih_mortality_si_cdpf)%r82d(io_si,icdpf) = &
                                this%hvars(ih_mortality_si_cdpf)%r82d(io_si,icdpf) + &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort + &
                                ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                ccohort%n * sec_per_day * days_per_year / m2_per_ha

                           ! crown damage by size by pft
                           this%hvars(ih_nplant_si_cdpf)%r82d(io_si, icdpf) = &
                                this%hvars(ih_nplant_si_cdpf)%r82d(io_si, icdpf) + ccohort%n / m2_per_ha
                           this%hvars(ih_m3_si_cdpf)%r82d(io_si, icdpf) = &
                                this%hvars(ih_m3_si_cdpf)%r82d(io_si, icdpf) + &
                                ccohort%cmort * ccohort%n / m2_per_ha

                           ! mortality
                           this%hvars(ih_m11_si_scpf)%r82d(io_si,scpf) = &
                                this%hvars(ih_m11_si_scpf)%r82d(io_si,scpf) + &
                                ccohort%dgmort*ccohort%n / m2_per_ha
                           this%hvars(ih_m11_si_cdpf)%r82d(io_si,icdpf) = &
                                this%hvars(ih_m11_si_cdpf)%r82d(io_si,icdpf) + &
                                ccohort%dgmort*ccohort%n / m2_per_ha

                           this%hvars(ih_ddbh_si_cdpf)%r82d(io_si,icdpf) = &
                                this%hvars(ih_ddbh_si_cdpf)%r82d(io_si,icdpf) + &
                                ccohort%ddbhdt*ccohort%n / m2_per_ha * m_per_cm

                        end if

                        ! Carbon only metrics
                        sapw_m   = ccohort%prt%GetState(sapw_organ, carbon12_element)
                        struct_m = ccohort%prt%GetState(struct_organ, carbon12_element)
                        leaf_m   = ccohort%prt%GetState(leaf_organ, carbon12_element)
                        fnrt_m   = ccohort%prt%GetState(fnrt_organ, carbon12_element)
                        store_m  = ccohort%prt%GetState(store_organ, carbon12_element)
                        repro_m  = ccohort%prt%GetState(repro_organ, carbon12_element)
                        alive_m  = leaf_m + fnrt_m + sapw_m
                        total_m  = alive_m + store_m + struct_m

                        hio_mortality_carbonflux_si_pft(io_si,ccohort%pft) = hio_mortality_carbonflux_si_pft(io_si,ccohort%pft) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                             ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * &
                             total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2 + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_m * &
                             ccohort%n * ha_per_m2


                        hio_hydraulicmortality_carbonflux_si_pft(io_si,ccohort%pft) = &
                             hio_hydraulicmortality_carbonflux_si_pft(io_si,ccohort%pft) + &
                             ccohort%hmort * total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2

                        hio_cstarvmortality_carbonflux_si_pft(io_si,ccohort%pft) = &
                             hio_cstarvmortality_carbonflux_si_pft(io_si,ccohort%pft) + &
                             ccohort%cmort * total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2

                        hio_cstarvmortality_continuous_carbonflux_si_pft(io_si,ccohort%pft) = &
                             hio_cstarvmortality_continuous_carbonflux_si_pft(io_si,ccohort%pft) + &
                             ccohort%cmort * total_m * ccohort%n * days_per_sec * years_per_day * ha_per_m2
                        
                        ! Aboveground mortality
                        hio_abg_mortality_cflux_si_scpf(io_si,scpf) = hio_abg_mortality_cflux_si_scpf(io_si,scpf) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * &
                             ( (sapw_m + struct_m + store_m ) * prt_params%allom_agb_frac(ccohort%pft) + &
                             leaf_m ) * ccohort%n * days_per_sec * years_per_day * ha_per_m2 + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                             ( (sapw_m + struct_m + store_m ) * prt_params%allom_agb_frac(ccohort%pft) + &
                             leaf_m ) * ccohort%n * ha_per_m2

                        ! Aboveground woody productivity
                        hio_abg_productivity_cflux_si_scpf(io_si,scpf) = hio_abg_productivity_cflux_si_scpf(io_si,scpf) + &
                             ( (sapw_m_net_alloc + struct_m_net_alloc + store_m_net_alloc) * prt_params%allom_agb_frac(ccohort%pft) + &
                             leaf_m_net_alloc ) * n_perm2 / &
                             days_per_year / sec_per_day


                        ! number density by size and biomass
                        hio_agb_si_scls(io_si,scls) = hio_agb_si_scls(io_si,scls) + &
                             total_m * ccohort%n * prt_params%allom_agb_frac(ccohort%pft) * AREA_INV

                        hio_agb_si_scpf(io_si,scpf) = hio_agb_si_scpf(io_si,scpf) + &
                             total_m * ccohort%n * prt_params%allom_agb_frac(ccohort%pft) * AREA_INV

                        hio_biomass_si_scls(io_si,scls) = hio_biomass_si_scls(io_si,scls) + &
                             total_m * ccohort%n * AREA_INV

                        ! age-resolved cohort-based areas

                        hio_npp_si_age(io_si,cpatch%age_class) = hio_npp_si_age(io_si,cpatch%age_class) + &
                             ccohort%n * ccohort%npp_acc_hold * AREA_INV / days_per_year / sec_per_day

                        ! update size-class x patch-age related quantities

                        iscag = get_sizeage_class_index(ccohort%dbh,cpatch%age)

                        hio_nplant_si_scag(io_si,iscag) = hio_nplant_si_scag(io_si,iscag) + ccohort%n / m2_per_ha

                        hio_nplant_si_scls(io_si,scls) = hio_nplant_si_scls(io_si,scls) + ccohort%n / m2_per_ha


                        ! update size, age, and PFT - indexed quantities
                        iscagpft = get_sizeagepft_class_index(ccohort%dbh,cpatch%age,ccohort%pft)

                        hio_nplant_si_scagpft(io_si,iscagpft) = hio_nplant_si_scagpft(io_si,iscagpft) + ccohort%n / m2_per_ha

                        ! update age and PFT - indexed quantities
                        iagepft = get_agepft_class_index(cpatch%age,ccohort%pft)

                        hio_npp_si_agepft(io_si,iagepft) = hio_npp_si_agepft(io_si,iagepft) + &
                             ccohort%n * ccohort%npp_acc_hold * AREA_INV / days_per_year / sec_per_day

                        hio_biomass_si_agepft(io_si,iagepft) = hio_biomass_si_agepft(io_si,iagepft) + &
                             total_m * ccohort%n * AREA_INV

                        ! update SCPF/SCLS- and canopy/subcanopy- partitioned quantities
                        canlayer: if (ccohort%canopy_layer .eq. 1) then
                           hio_nplant_canopy_si_scag(io_si,iscag) = hio_nplant_canopy_si_scag(io_si,iscag) + ccohort%n / m2_per_ha
                           hio_mortality_canopy_si_scag(io_si,iscag) = hio_mortality_canopy_si_scag(io_si,iscag) + &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                                ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha
                           hio_ddbh_canopy_si_scag(io_si,iscag) = hio_ddbh_canopy_si_scag(io_si,iscag) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha
                           hio_bstor_canopy_si_scpf(io_si,scpf) = hio_bstor_canopy_si_scpf(io_si,scpf) + &
                                store_m * ccohort%n / m2_per_ha
                           hio_bleaf_canopy_si_scpf(io_si,scpf) = hio_bleaf_canopy_si_scpf(io_si,scpf) + &
                                leaf_m * ccohort%n / m2_per_ha
                           hio_lai_canopy_si_scpf(io_si,scpf) = hio_lai_canopy_si_scpf(io_si,scpf) + &
                                ccohort%treelai*ccohort%c_area * AREA_INV
                           
                           hio_crownarea_canopy_si_scpf(io_si,scpf) = hio_crownarea_canopy_si_scpf(io_si,scpf) + &
                                ccohort%c_area * AREA_INV
                           !hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &
                           !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                           !     ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                           hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort + &
                                ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                ccohort%n * sec_per_day * days_per_year / m2_per_ha

                           hio_m3_mortality_canopy_si_scpf(io_si,scpf) = hio_m3_mortality_canopy_si_scpf(io_si,scpf) + &
                                ccohort%cmort * ccohort%n / m2_per_ha

                           hio_nplant_canopy_si_scpf(io_si,scpf) = hio_nplant_canopy_si_scpf(io_si,scpf) + ccohort%n / m2_per_ha
                           hio_nplant_canopy_si_scls(io_si,scls) = hio_nplant_canopy_si_scls(io_si,scls) + ccohort%n / m2_per_ha
                           hio_lai_canopy_si_scls(io_si,scls) = hio_lai_canopy_si_scls(io_si,scls) + &
                                ccohort%treelai*ccohort%c_area * AREA_INV
                           hio_sai_canopy_si_scls(io_si,scls) = hio_sai_canopy_si_scls(io_si,scls) + &
                                ccohort%treesai*ccohort%c_area * AREA_INV
                           hio_trimming_canopy_si_scls(io_si,scls) = hio_trimming_canopy_si_scls(io_si,scls) + &
                                ccohort%n * ccohort%canopy_trim / m2_per_ha
                           hio_crown_area_canopy_si_scls(io_si,scls) = hio_crown_area_canopy_si_scls(io_si,scls) + &
                                ccohort%c_area * AREA_INV
                           hio_gpp_canopy_si_scpf(io_si,scpf) = hio_gpp_canopy_si_scpf(io_si,scpf) +  &
                                n_perm2*ccohort%gpp_acc_hold / days_per_year / sec_per_day
                           hio_ar_canopy_si_scpf(io_si,scpf) = hio_ar_canopy_si_scpf(io_si,scpf) + &
                                n_perm2*(ccohort%resp_m_acc_hold + ccohort%resp_g_acc_hold + &
                                ccohort%resp_excess_hold) / days_per_year / sec_per_day
                           ! growth increment
                           hio_ddbh_canopy_si_scpf(io_si,scpf) = hio_ddbh_canopy_si_scpf(io_si,scpf) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha
                           hio_ddbh_canopy_si_scls(io_si,scls) = hio_ddbh_canopy_si_scls(io_si,scls) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha

                           ! sum of all mortality
                           hio_mortality_canopy_si_scls(io_si,scls) = hio_mortality_canopy_si_scls(io_si,scls) + &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort +   &
                                ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                ccohort%n * sec_per_day * days_per_year / m2_per_ha

                           hio_m3_mortality_canopy_si_scls(io_si,scls) = hio_m3_mortality_canopy_si_scls(io_si,scls) + &
                                ccohort%cmort * ccohort%n / m2_per_ha

                           hio_carbon_balance_canopy_si_scls(io_si,scls) = hio_carbon_balance_canopy_si_scls(io_si,scls) + &
                                ccohort%n * ccohort%npp_acc_hold / m2_per_ha / days_per_year / sec_per_day

                           ! damage variables - canopy
                           if(hlm_use_tree_damage .eq. itrue) then

                              ! carbon starvation mortality in the canopy by size x damage x pft 
                              this%hvars(ih_m3_mortality_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_m3_mortality_canopy_si_cdpf)%r82d(io_si,icdpf)+&
                                   ccohort%cmort * ccohort%n / m2_per_ha

                              ! damage mortality in the canopy by size x damage x pft 
                              this%hvars(ih_m11_mortality_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_m11_mortality_canopy_si_cdpf)%r82d(io_si,icdpf)+&
                                   ccohort%dgmort * ccohort%n / m2_per_ha

                              this%hvars(ih_mortality_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_mortality_canopy_si_cdpf)%r82d(io_si,icdpf)+ &
                                   (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort +  ccohort%smort + &
                                   ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                   (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                   ccohort%n * sec_per_day * days_per_year / m2_per_ha

                              ! nplants by damage 
                              this%hvars(ih_nplant_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_nplant_canopy_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%n / m2_per_ha

                              ! growth rate by damage x size x pft in the canopy
                              this%hvars(ih_ddbh_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_ddbh_canopy_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%ddbhdt*ccohort%n / m2_per_ha * m_per_cm

                           end if ! end if damage


                           hio_leaf_md_canopy_si_scls(io_si,scls) = hio_leaf_md_canopy_si_scls(io_si,scls) + &
                                leaf_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_root_md_canopy_si_scls(io_si,scls) = hio_root_md_canopy_si_scls(io_si,scls) + &
                                fnrt_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bsw_md_canopy_si_scls(io_si,scls) = hio_bsw_md_canopy_si_scls(io_si,scls) + &
                                sapw_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bstore_md_canopy_si_scls(io_si,scls) = hio_bstore_md_canopy_si_scls(io_si,scls) + &
                                store_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bdead_md_canopy_si_scls(io_si,scls) = hio_bdead_md_canopy_si_scls(io_si,scls) + &
                                struct_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_seed_prod_canopy_si_scls(io_si,scls) = hio_seed_prod_canopy_si_scls(io_si,scls) + &
                                ccohort%seed_prod * ccohort%n / m2_per_ha / sec_per_day

                           hio_npp_leaf_canopy_si_scls(io_si,scls) = hio_npp_leaf_canopy_si_scls(io_si,scls) + &
                                leaf_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_fnrt_canopy_si_scls(io_si,scls) = hio_npp_fnrt_canopy_si_scls(io_si,scls) + &
                                fnrt_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_sapw_canopy_si_scls(io_si,scls) = hio_npp_sapw_canopy_si_scls(io_si,scls) + &
                                sapw_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_dead_canopy_si_scls(io_si,scls) = hio_npp_dead_canopy_si_scls(io_si,scls) + &
                                struct_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_seed_canopy_si_scls(io_si,scls) = hio_npp_seed_canopy_si_scls(io_si,scls) + &
                                repro_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_stor_canopy_si_scls(io_si,scls) = hio_npp_stor_canopy_si_scls(io_si,scls) + &
                                store_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day

                           hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) = &
                                hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) + &
                                ccohort%canopy_layer_yesterday * ccohort%n / m2_per_ha


                        else canlayer
                           hio_nplant_understory_si_scag(io_si,iscag) = hio_nplant_understory_si_scag(io_si,iscag) + ccohort%n / m2_per_ha
                           hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                                ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha
                           hio_ddbh_understory_si_scag(io_si,iscag) = hio_ddbh_understory_si_scag(io_si,iscag) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha
                           hio_bstor_understory_si_scpf(io_si,scpf) = hio_bstor_understory_si_scpf(io_si,scpf) + &
                                store_m * ccohort%n / m2_per_ha
                           hio_bleaf_understory_si_scpf(io_si,scpf) = hio_bleaf_understory_si_scpf(io_si,scpf) + &
                                leaf_m  * ccohort%n / m2_per_ha

                           hio_lai_understory_si_scpf(io_si,scpf) = hio_lai_understory_si_scpf(io_si,scpf) + &
                                ccohort%treelai*ccohort%c_area  * AREA_INV
                           hio_crownarea_understory_si_scpf(io_si,scpf) = hio_crownarea_understory_si_scpf(io_si,scpf) + &
                                ccohort%c_area * AREA_INV
                           
                           !hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                           !    (ccohort%bmort + ccohort%hmort + ccohort%cmort +
                           !      ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                           hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                                ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                ccohort%n * sec_per_day * days_per_year / m2_per_ha

                           hio_m3_mortality_understory_si_scpf(io_si,scpf) = hio_m3_mortality_understory_si_scpf(io_si,scpf) + &
                                ccohort%cmort * ccohort%n / m2_per_ha

                           if(cpatch%land_use_label .eq. secondaryland) then
                              hio_mortality_canopy_secondary_si_scls(io_si,scls) = hio_mortality_canopy_secondary_si_scls(io_si,scls) + &
                                   (ccohort%bmort + ccohort%hmort + ccohort%cmort +   &
                                   ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n / m2_per_ha + &
                                   (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                   ccohort%n * sec_per_day * days_per_year / m2_per_ha
                           end if

                           hio_nplant_understory_si_scpf(io_si,scpf) = hio_nplant_understory_si_scpf(io_si,scpf) + ccohort%n / m2_per_ha
                           hio_nplant_understory_si_scls(io_si,scls) = hio_nplant_understory_si_scls(io_si,scls) + ccohort%n / m2_per_ha
                           hio_lai_understory_si_scls(io_si,scls) = hio_lai_understory_si_scls(io_si,scls) + &
                                ccohort%treelai*ccohort%c_area  * AREA_INV
                           hio_sai_understory_si_scls(io_si,scls) = hio_sai_understory_si_scls(io_si,scls) + &
                                ccohort%treelai*ccohort%c_area  * AREA_INV
                           hio_trimming_understory_si_scls(io_si,scls) = hio_trimming_understory_si_scls(io_si,scls) + &
                                ccohort%n * ccohort%canopy_trim / m2_per_ha
                           hio_crown_area_understory_si_scls(io_si,scls) = hio_crown_area_understory_si_scls(io_si,scls) + &
                                ccohort%c_area * AREA_INV
                           hio_gpp_understory_si_scpf(io_si,scpf)      = hio_gpp_understory_si_scpf(io_si,scpf)      + &
                                n_perm2*ccohort%gpp_acc_hold / days_per_year / sec_per_day
                           hio_ar_understory_si_scpf(io_si,scpf)      = hio_ar_understory_si_scpf(io_si,scpf)      + &
                                n_perm2*(ccohort%resp_m_acc_hold + ccohort%resp_g_acc_hold + &
                                ccohort%resp_excess_hold) / days_per_year / sec_per_day

                           ! growth increment
                           hio_ddbh_understory_si_scpf(io_si,scpf) = hio_ddbh_understory_si_scpf(io_si,scpf) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha
                           hio_ddbh_understory_si_scls(io_si,scls) = hio_ddbh_understory_si_scls(io_si,scls) + &
                                ccohort%ddbhdt*ccohort%n * m_per_cm / m2_per_ha

                           ! sum of all mortality
                           hio_mortality_understory_si_scls(io_si,scls) = hio_mortality_understory_si_scls(io_si,scls) + &
                                (ccohort%bmort + ccohort%hmort + ccohort%cmort +   &
                                ccohort%frmort + ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                ccohort%n * sec_per_day * days_per_year / m2_per_ha

                           hio_m3_mortality_understory_si_scls(io_si,scls) = hio_m3_mortality_understory_si_scls(io_si,scls) + &
                                ccohort%cmort * ccohort%n / m2_per_ha



                           hio_carbon_balance_understory_si_scls(io_si,scls) = hio_carbon_balance_understory_si_scls(io_si,scls) + &
                                ccohort%npp_acc_hold * ccohort%n / m2_per_ha / days_per_year / sec_per_day

                           ! damage variables - understory
                           if(hlm_use_tree_damage .eq. itrue) then

                              ! carbon mortality in the understory by damage x size x pft
                              this%hvars(ih_m3_mortality_understory_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_m3_mortality_understory_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%cmort * ccohort%n / m2_per_ha

                              ! damage in the understory by damage x size x pft
                              this%hvars(ih_m11_mortality_understory_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_m11_mortality_understory_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%dgmort * ccohort%n / m2_per_ha

                              ! total mortality of understory cohorts by damage x size x pft
                              this%hvars(ih_mortality_understory_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_mortality_understory_si_cdpf)%r82d(io_si,icdpf) + &
                                   (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort + &
                                   ccohort%smort + ccohort%asmort + ccohort%dgmort) * ccohort%n / m2_per_ha + &
                                   (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                                   ccohort%n * sec_per_day * days_per_year / m2_per_ha

                              this%hvars(ih_nplant_understory_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_nplant_understory_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%n / m2_per_ha

                              ! growth rate by size x damage x pft  - understory
                              this%hvars(ih_ddbh_understory_si_cdpf)%r82d(io_si,icdpf) = &
                                   this%hvars(ih_ddbh_understory_si_cdpf)%r82d(io_si,icdpf) + &
                                   ccohort%ddbhdt*ccohort%n / m2_per_ha * m_per_cm

                           end if ! end if damage

                           hio_leaf_md_understory_si_scls(io_si,scls) = hio_leaf_md_understory_si_scls(io_si,scls) + &
                                leaf_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_root_md_understory_si_scls(io_si,scls) = hio_root_md_understory_si_scls(io_si,scls) + &
                                fnrt_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bsw_md_understory_si_scls(io_si,scls) = hio_bsw_md_understory_si_scls(io_si,scls) + &
                                sapw_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bstore_md_understory_si_scls(io_si,scls) = hio_bstore_md_understory_si_scls(io_si,scls) + &
                                store_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_bdead_md_understory_si_scls(io_si,scls) = hio_bdead_md_understory_si_scls(io_si,scls) + &
                                struct_m_turnover * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_seed_prod_understory_si_scls(io_si,scls) = hio_seed_prod_understory_si_scls(io_si,scls) + &
                                ccohort%seed_prod * ccohort%n / m2_per_ha  / sec_per_day

                           hio_npp_leaf_understory_si_scls(io_si,scls) = hio_npp_leaf_understory_si_scls(io_si,scls) + &
                                leaf_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_fnrt_understory_si_scls(io_si,scls) = hio_npp_fnrt_understory_si_scls(io_si,scls) + &
                                fnrt_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_sapw_understory_si_scls(io_si,scls) = hio_npp_sapw_understory_si_scls(io_si,scls) + &
                                sapw_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_dead_understory_si_scls(io_si,scls) = hio_npp_dead_understory_si_scls(io_si,scls) + &
                                struct_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_seed_understory_si_scls(io_si,scls) = hio_npp_seed_understory_si_scls(io_si,scls) + &
                                repro_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day
                           hio_npp_stor_understory_si_scls(io_si,scls) = hio_npp_stor_understory_si_scls(io_si,scls) + &
                                store_m_net_alloc * ccohort%n / m2_per_ha / days_per_year / sec_per_day

                           hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) = &
                                hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) + &
                                ccohort%canopy_layer_yesterday * ccohort%n / m2_per_ha
                        endif canlayer
                        !
                        !
                        ccohort%canopy_layer_yesterday = real(ccohort%canopy_layer, r8)
                        !
                        ! growth flux of individuals into a given bin
                        ! track the actual growth here, the virtual growth from fusion lower down
                        if ( (scls - ccohort%size_class_lasttimestep ) .gt. 0) then
                           do i_scls = ccohort%size_class_lasttimestep + 1, scls
                              i_scpf = (ccohort%pft-1)*nlevsclass+i_scls
                              hio_growthflux_si_scpf(io_si,i_scpf) = hio_growthflux_si_scpf(io_si,i_scpf) + &
                                   ccohort%n * days_per_year / m2_per_ha
                           end do
                        end if
                        ccohort%size_class_lasttimestep = scls

                      end associate
                   else notnew ! i.e. cohort%isnew
                      !
                      ! if cohort is new, track its growth flux into the first size bin
                      i_scpf = (ccohort%pft-1)*nlevsclass+1
                      hio_growthflux_si_scpf(io_si,i_scpf) =                          &
                           hio_growthflux_si_scpf(io_si,i_scpf) + ccohort%n *           &
                           days_per_year / m2_per_ha
                      ccohort%size_class_lasttimestep = 1

                   end if notnew

                   ! resolve some canopy area profiles, both total and of occupied leaves
                   ican = ccohort%canopy_layer
                   !
                   hio_crownarea_cl(io_si, ican) = hio_crownarea_cl(io_si, ican) + ccohort%c_area / AREA
                   !
                   do ileaf=1,ccohort%nv
                      cnlf_indx = ileaf + (ican-1) * nlevleaf
                      hio_crownarea_si_cnlf(io_si, cnlf_indx) = hio_crownarea_si_cnlf(io_si, cnlf_indx) + &
                           ccohort%c_area / AREA
                   end do

                   ccohort => ccohort%taller
                enddo cohortloop ! cohort loop



                do ilyr = 1,sites(s)%nlevsoil
                   hio_fragmentation_scaler_sl(io_si,ilyr) = hio_fragmentation_scaler_sl(io_si,ilyr) + cpatch%fragmentation_scaler(ilyr) * cpatch%area * AREA_INV
                end do

                do i_fuel = 1, num_fuel_classes

                   i_agefuel = get_agefuel_class_index(cpatch%age,i_fuel)
                   hio_fuel_amount_age_fuel(io_si,i_agefuel) = hio_fuel_amount_age_fuel(io_si,i_agefuel) + &
                        cpatch%fuel%frac_loading(i_fuel) * cpatch%fuel%non_trunk_loading * cpatch%area * AREA_INV

                   hio_litter_moisture_si_fuel(io_si, i_fuel) = hio_litter_moisture_si_fuel(io_si, i_fuel) + &
                        cpatch%fuel%effective_moisture(i_fuel) * cpatch%area * AREA_INV

                   hio_fuel_amount_si_fuel(io_si, i_fuel) = hio_fuel_amount_si_fuel(io_si, i_fuel) + &
                        cpatch%fuel%frac_loading(i_fuel) * cpatch%fuel%non_trunk_loading * cpatch%area * AREA_INV

                   hio_burnt_frac_litter_si_fuel(io_si, i_fuel) = hio_burnt_frac_litter_si_fuel(io_si, i_fuel) + &
                        cpatch%fuel%frac_burnt(i_fuel) * cpatch%frac_burnt * cpatch%area * AREA_INV
                end do




                ! Update Litter Flux Variables

                litt_c       => cpatch%litter(element_pos(carbon12_element))


                do i_cwd = 1, ncwd

                   hio_cwd_ag_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_si_cwdsc(io_si, i_cwd) + &
                        litt_c%ag_cwd(i_cwd)*cpatch%area * AREA_INV
                   hio_cwd_bg_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_si_cwdsc(io_si, i_cwd) + &
                        sum(litt_c%bg_cwd(i_cwd,:)) * cpatch%area * AREA_INV

                   hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) + &
                        litt_c%ag_cwd_frag(i_cwd)*cpatch%area * AREA_INV /              &
                        days_per_year / sec_per_day

                   hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) + &
                        sum(litt_c%bg_cwd_frag(i_cwd,:)) * cpatch%area * AREA_INV / &
                        days_per_year / sec_per_day

                end do

                cpatch => cpatch%younger
             end do patchloop !patch loop




             ! divide so-far-just-summed but to-be-averaged patch-age-class
             ! variables by patch-age-class area to get mean values
             do ipa2 = 1, nlevage
                if (hio_area_si_age(io_si, ipa2) .gt. nearzero) then
                   hio_lai_si_age(io_si, ipa2) = hio_lai_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
                   hio_ncl_si_age(io_si, ipa2) = hio_ncl_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
                   do ft = 1, numpft
                      iagepft = ipa2 + (ft-1) * nlevage
                      hio_scorch_height_si_agepft(io_si, iagepft) = &
                           hio_scorch_height_si_agepft(io_si, iagepft) / (hio_area_si_age(io_si, ipa2)*AREA)
                   enddo
                else
                   hio_lai_si_age(io_si, ipa2) = 0._r8
                   hio_ncl_si_age(io_si, ipa2) = 0._r8
                endif
             end do



             ! pass the cohort termination mortality as a flux to the history, and then reset the termination mortality buffer
             ! note there are various ways of reporting the total mortality, so pass to these as well
             do i_pft = 1, numpft
                hio_cstarvmortality_carbonflux_si_pft(io_si,i_pft) = &
                     hio_cstarvmortality_carbonflux_si_pft(io_si,i_pft) + &
                     (sites(s)%term_carbonflux_ustory(i_term_mort_type_cstarv,i_pft) + &
                     sites(s)%term_carbonflux_canopy(i_term_mort_type_cstarv,i_pft)) * days_per_sec * ha_per_m2
             end do
             do ft = 1, numpft
                do i_scls = 1,nlevsclass
                   i_scpf = (ft-1)*nlevsclass + i_scls
                   !
                   ! termination mortality. sum of canopy and understory indices
                   
                   ! move carbon starvation-related termination mortality to the carbon starvation mortality type and only consider
                   ! the other two types of termination mortality here.
                   
                   hio_m6_si_scpf(io_si,i_scpf) = &
                        (sum(sites(s)%term_nindivs_canopy(i_term_mort_type_canlev:n_term_mort_types,i_scls,ft)) + &
                        sum(sites(s)%term_nindivs_ustory(i_term_mort_type_canlev:n_term_mort_types,i_scls,ft))) *              &
                        days_per_year / m2_per_ha
                   
                   hio_m6_si_scls(io_si,i_scls) = hio_m6_si_scls(io_si,i_scls) +  &
                        (sum(sites(s)%term_nindivs_canopy(i_term_mort_type_canlev:n_term_mort_types,i_scls,ft)) +              &
                        sum(sites(s)%term_nindivs_ustory(i_term_mort_type_canlev:n_term_mort_types,i_scls,ft))) *              &
                        days_per_year / m2_per_ha
                   !
                   ! add the carbon starvation-related termination mortality to the carbon starvation diagnostics
                   hio_m3_si_scpf(io_si,i_scpf) = hio_m3_si_scpf(io_si,i_scpf) +                               &
                        (sites(s)%term_nindivs_canopy(i_term_mort_type_cstarv,i_scls,ft) +              &
                        sites(s)%term_nindivs_ustory(i_term_mort_type_cstarv,i_scls,ft)) *              &
                        days_per_year / m2_per_ha     

                   hio_m3_si_scls(io_si,i_scls) = hio_m3_si_scls(io_si,i_scls) +                               &
                        (sites(s)%term_nindivs_canopy(i_term_mort_type_cstarv,i_scls,ft) +              &
                        sites(s)%term_nindivs_ustory(i_term_mort_type_cstarv,i_scls,ft)) *              &
                        days_per_year / m2_per_ha
                   
                   ! add c-starve termination mortality to canopy and understory M3 mortality (N/m^2/yr)
                   hio_m3_mortality_canopy_si_scpf(io_si,i_scpf) = &
                        hio_m3_mortality_canopy_si_scpf(io_si,i_scpf) + &
                        sites(s)%term_nindivs_canopy(i_term_mort_type_cstarv,i_scls,ft) * &
                        days_per_year / m2_per_ha

                   hio_m3_mortality_understory_si_scpf(io_si,i_scpf) = &
                        hio_m3_mortality_understory_si_scpf(io_si,i_scpf) + &
                        sites(s)%term_nindivs_ustory(i_term_mort_type_cstarv,i_scls,ft) * &
                        days_per_year / m2_per_ha
                   
                   hio_m3_mortality_canopy_si_scls(io_si,i_scls) = &
                        hio_m3_mortality_canopy_si_scls(io_si,i_scls) + &
                        sites(s)%term_nindivs_canopy(i_term_mort_type_cstarv,i_scls,ft) * &
                        days_per_year / m2_per_ha
                   
                   hio_m3_mortality_understory_si_scls(io_si,i_scls) = &
                        hio_m3_mortality_understory_si_scls(io_si,i_scls) + &
                        sites(s)%term_nindivs_ustory(i_term_mort_type_cstarv,i_scls,ft) * &
                        days_per_year / m2_per_ha
                   
                   !
                   ! add termination mortality to canopy and understory mortality
                   hio_mortality_canopy_si_scls(io_si,i_scls) = hio_mortality_canopy_si_scls(io_si,i_scls) + &
                        sum(sites(s)%term_nindivs_canopy(:,i_scls,ft)) * days_per_year / m2_per_ha
                   
                   hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                        sum(sites(s)%term_nindivs_ustory(:,i_scls,ft)) * days_per_year / m2_per_ha

                   hio_mortality_canopy_si_scpf(io_si,i_scpf) = hio_mortality_canopy_si_scpf(io_si,i_scpf) + &
                        sum(sites(s)%term_nindivs_canopy(:,i_scls,ft)) * days_per_year / m2_per_ha

                   hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                        sum(sites(s)%term_nindivs_ustory(:,i_scls,ft)) * days_per_year / m2_per_ha


                   !
                   ! imort on its own
                   hio_m4_si_scpf(io_si,i_scpf) = sites(s)%imort_rate(i_scls, ft) / m2_per_ha
                   hio_m4_si_scls(io_si,i_scls) = hio_m4_si_scls(io_si,i_scls) + sites(s)%imort_rate(i_scls, ft) / m2_per_ha
                   !
                   ! add imort to other mortality terms. consider imort as understory mortality even if it happens in
                   ! cohorts that may have been promoted as part of the patch creation, and use the pre-calculated site-level
                   ! values to avoid biasing the results by the dramatically-reduced number densities in cohorts that are subject to imort
                   hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                        sites(s)%imort_rate(i_scls, ft) / m2_per_ha
                   hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                        sites(s)%imort_rate(i_scls, ft) / m2_per_ha
                   !
                   iscag = i_scls ! since imort is by definition something that only happens in newly disturbed patches, treat as such
                   hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                        sites(s)%imort_rate(i_scls, ft) / m2_per_ha

                   ! fire mortality from the site-level diagnostic rates
                   hio_m5_si_scpf(io_si,i_scpf) = (sites(s)%fmort_rate_canopy(i_scls, ft) + &
                        sites(s)%fmort_rate_ustory(i_scls, ft)) / m2_per_ha
                   hio_m5_si_scls(io_si,i_scls) = hio_m5_si_scls(io_si,i_scls) + &
                        (sites(s)%fmort_rate_canopy(i_scls, ft) +              &
                        sites(s)%fmort_rate_ustory(i_scls, ft)) / m2_per_ha
                   !
                   hio_crownfiremort_si_scpf(io_si,i_scpf) = sites(s)%fmort_rate_crown(i_scls, ft) / m2_per_ha
                   hio_cambialfiremort_si_scpf(io_si,i_scpf) = sites(s)%fmort_rate_cambial(i_scls, ft) / m2_per_ha
                   !
                   ! fire components of overall canopy and understory mortality
                   hio_mortality_canopy_si_scpf(io_si,i_scpf) = hio_mortality_canopy_si_scpf(io_si,i_scpf) + &
                        sites(s)%fmort_rate_canopy(i_scls, ft) / m2_per_ha
                   hio_mortality_canopy_si_scls(io_si,i_scls) = hio_mortality_canopy_si_scls(io_si,i_scls) + &
                        sites(s)%fmort_rate_canopy(i_scls, ft) / m2_per_ha

                   ! the fire mortality rates for each layer are total dead, since the usable
                   ! output will then normalize by the counts, we are allowed to sum over layers
                   hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                        sites(s)%fmort_rate_ustory(i_scls, ft) / m2_per_ha

                   hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                        sites(s)%fmort_rate_ustory(i_scls, ft) / m2_per_ha

                   !
                   ! for scag variables, also treat as happening in the newly-disurbed patch

                   hio_mortality_canopy_si_scag(io_si,iscag) = hio_mortality_canopy_si_scag(io_si,iscag) + &
                        sites(s)%fmort_rate_canopy(i_scls, ft) / m2_per_ha
                   hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                        sites(s)%fmort_rate_ustory(i_scls, ft) / m2_per_ha

                   ! while in this loop, pass the fusion-induced growth rate flux to history
                   hio_growthflux_fusion_si_scpf(io_si,i_scpf) = hio_growthflux_fusion_si_scpf(io_si,i_scpf) + &
                        sites(s)%growthflux_fusion(i_scls, ft) * days_per_year / m2_per_ha

                end do
             end do



             do ft = 1, numpft
                hio_mortality_carbonflux_si_pft(io_si,ft) = hio_mortality_carbonflux_si_pft(io_si,ft) + &
                     (sites(s)%fmort_carbonflux_canopy(ft) + &
                     sites(s)%fmort_carbonflux_ustory(ft) ) / g_per_kg + &
                     sites(s)%imort_carbonflux(ft) + &
                     sum(sites(s)%term_carbonflux_ustory(:,ft)) * days_per_sec * ha_per_m2 + &
                     sum(sites(s)%term_carbonflux_canopy(:,ft)) * days_per_sec * ha_per_m2

                hio_firemortality_carbonflux_si_pft(io_si,ft) = sites(s)%fmort_carbonflux_canopy(ft) / g_per_kg
             end do

             ! add imort and fmort to aboveground woody mortality 
             do ft = 1, numpft
                do i_scls = 1,nlevsclass
                   i_scpf = (ft-1)*nlevsclass + i_scls
                   hio_abg_mortality_cflux_si_scpf(io_si,i_scpf) = hio_abg_mortality_cflux_si_scpf(io_si,i_scpf) + &
                        (sites(s)%fmort_abg_flux(i_scls,ft) / g_per_kg ) + &
                        sites(s)%imort_abg_flux(i_scls,ft)  +  &
                        (sites(s)%term_abg_flux(i_scls,ft)  * days_per_sec * ha_per_m2 ) 
                end do
             end do


             if(hlm_use_tree_damage .eq. itrue) then

                do ft = 1, numpft
                   do icdam = 1, nlevdamage
                      do i_scls = 1,nlevsclass

                         icdsc = (icdam-1)*nlevsclass + i_scls
                         icdpf = (icdam-1)*nlevsclass + i_scls + &
                              (ft-1) * nlevsclass * nlevdamage

                         this%hvars(ih_mortality_si_cdpf)%r82d(io_si, icdpf) = &
                              this%hvars(ih_mortality_si_cdpf)%r82d(io_si, icdpf) + &
                              ( (sites(s)%term_nindivs_canopy_damage(icdam, i_scls, ft) * days_per_year) + &
                              (sites(s)%term_nindivs_ustory_damage(icdam, i_scls, ft) * days_per_year) + &
                              sites(s)%imort_rate_damage(icdam, i_scls, ft) + & 
                              sites(s)%fmort_rate_canopy_damage(icdam, i_scls, ft) + &
                              sites(s)%fmort_rate_ustory_damage(icdam, i_scls, ft) ) / m2_per_ha

                         this%hvars(ih_mortality_canopy_si_cdpf)%r82d(io_si,icdpf) = &
                              this%hvars(ih_mortality_canopy_si_cdpf)%r82d(io_si,icdpf) + &
                              ( sites(s)%term_nindivs_canopy_damage(icdam,i_scls,ft) * days_per_year + &
                              sites(s)%fmort_rate_canopy_damage(icdam, i_scls, ft) )/ m2_per_ha

                         this%hvars(ih_mortality_understory_si_cdpf)%r82d(io_si,icdpf) = &
                              this%hvars(ih_mortality_understory_si_cdpf)%r82d(io_si,icdpf) + &
                              ( sites(s)%term_nindivs_ustory_damage(icdam, i_scls,ft) * days_per_year + &
                              sites(s)%imort_rate_damage(icdam, i_scls, ft) + &
                              sites(s)%fmort_rate_ustory_damage(icdam, i_scls, ft) )/ m2_per_ha

                      end do
                   end do
                end do
             end if
             sites(s)%term_nindivs_canopy(:,:,:) = 0._r8
             sites(s)%term_nindivs_ustory(:,:,:) = 0._r8
             sites(s)%imort_carbonflux(:) = 0._r8
             sites(s)%imort_rate(:,:) = 0._r8
             sites(s)%fmort_rate_canopy(:,:) = 0._r8
             sites(s)%fmort_rate_ustory(:,:) = 0._r8
             sites(s)%fmort_carbonflux_canopy(:) = 0._r8
             sites(s)%fmort_carbonflux_ustory(:) = 0._r8
             sites(s)%fmort_rate_cambial(:,:) = 0._r8
             sites(s)%fmort_rate_crown(:,:) = 0._r8
             sites(s)%growthflux_fusion(:,:) = 0._r8
             sites(s)%fmort_abg_flux(:,:) = 0._r8
             sites(s)%imort_abg_flux(:,:) = 0._r8
             sites(s)%term_abg_flux(:,:) = 0._r8

             sites(s)%imort_rate_damage(:,:,:) = 0.0_r8
             sites(s)%term_nindivs_canopy_damage(:,:,:) = 0.0_r8
             sites(s)%term_nindivs_ustory_damage(:,:,:) = 0.0_r8
             sites(s)%imort_cflux_damage(:,:) = 0._r8
             sites(s)%term_cflux_canopy_damage(:,:) = 0._r8
             sites(s)%term_cflux_ustory_damage(:,:) = 0._r8
             sites(s)%fmort_rate_canopy_damage(:,:,:) = 0._r8
             sites(s)%fmort_rate_ustory_damage(:,:,:) = 0._r8
             sites(s)%fmort_cflux_canopy_damage(:,:) = 0._r8
             sites(s)%fmort_cflux_ustory_damage(:,:) = 0._r8

             ! pass the recruitment rate as a flux to the history, and then reset the recruitment buffer
             do ft = 1, numpft
                ! pass the recruitment rate as a flux to the history, and then reset the recruitment buffer
                hio_recruitment_si_pft(io_si,ft) = sites(s)%recruitment_rate(ft) * days_per_year / m2_per_ha

                ! Gridcell output and inputs
                hio_seeds_out_gc_si_pft(io_si,ft) = sites(s)%seed_out(ft)
                hio_seeds_in_gc_si_pft(io_si,ft) = sites(s)%seed_in(ft)
             end do
             sites(s)%recruitment_rate(:) = 0._r8

             ! summarize all of the mortality fluxes by PFT
             do ft = 1, numpft
                do i_scls = 1,nlevsclass
                   i_scpf = (ft-1)*nlevsclass + i_scls

                   hio_mortality_si_pft(io_si,ft) = hio_mortality_si_pft(io_si,ft) + &
                        hio_m1_si_scpf(io_si,i_scpf) + &
                        hio_m2_si_scpf(io_si,i_scpf) + &
                        hio_m3_si_scpf(io_si,i_scpf) + &
                        hio_m4_si_scpf(io_si,i_scpf) + &
                        hio_m5_si_scpf(io_si,i_scpf) + &
                        hio_m6_si_scpf(io_si,i_scpf) + &
                        hio_m7_si_scpf(io_si,i_scpf) + &
                        hio_m8_si_scpf(io_si,i_scpf) + &
                        hio_m9_si_scpf(io_si,i_scpf) + &
                        hio_m10_si_scpf(io_si,i_scpf)

                   if(hlm_use_tree_damage .eq. itrue) then
                      hio_mortality_si_pft(io_si, ft) = hio_mortality_si_pft(io_si,ft) + &
                           this%hvars(ih_m11_si_scpf)%r82d(io_si,i_scpf)
                   end if

                end do
             end do

             ! ------------------------------------------------------------------------------
             ! Some carbon only litter diagnostics (legacy)
             ! ------------------------------------------------------------------------------

             elflux_diags_c => sites(s)%flux_diags%elem(element_pos(carbon12_element))

             ! ------------------------------------------------------------------------------
             ! Diagnostics discretized by element type
             ! ------------------------------------------------------------------------------

             do el = 1, num_elements

                elflux_diags => sites(s)%flux_diags%elem(el)

                ! Sum up all input litter fluxes (above below, fines, cwd) [kg/ha/day]
                hio_litter_in_elem(io_si, el) = (sum(elflux_diags%cwd_ag_input(:)) +    &
                     sum(elflux_diags%cwd_bg_input(:)) + sum(elflux_diags%surf_fine_litter_input(:)) + &
                     sum(elflux_diags%root_litter_input(:))) / m2_per_ha / sec_per_day


                ! Plant multi-element states and fluxes
                ! Zero states, and set the fluxes
                if(element_list(el).eq.carbon12_element)then
                   this%hvars(ih_totvegc_scpf)%r82d(io_si,:) = 0._r8
                   this%hvars(ih_leafc_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_fnrtc_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_sapwc_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_storec_scpf)%r82d(io_si,:)  = 0._r8
                   this%hvars(ih_reproc_scpf)%r82d(io_si,:)  = 0._r8

                elseif(element_list(el).eq.nitrogen_element)then

                   this%hvars(ih_totvegn_scpf)%r82d(io_si,:) = 0._r8
                   this%hvars(ih_leafn_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_fnrtn_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_sapwn_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_storen_scpf)%r82d(io_si,:)  = 0._r8
                   this%hvars(ih_repron_scpf)%r82d(io_si,:)  = 0._r8

                elseif(element_list(el).eq.phosphorus_element)then
                   this%hvars(ih_totvegp_scpf)%r82d(io_si,:) = 0._r8
                   this%hvars(ih_leafp_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_fnrtp_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_sapwp_scpf)%r82d(io_si,:)   = 0._r8
                   this%hvars(ih_storep_scpf)%r82d(io_si,:)  = 0._r8
                   this%hvars(ih_reprop_scpf)%r82d(io_si,:)  = 0._r8



                end if

                cpatch => sites(s)%oldest_patch
                do while(associated(cpatch))

                   litt => cpatch%litter(el)

                   area_frac = cpatch%area * AREA_INV

                   ! Sum up all output fluxes (fragmentation)
                   hio_litter_out_elem(io_si,el) = hio_litter_out_elem(io_si,el) + &
                        (sum(litt%leaf_fines_frag(:)) + &
                        sum(litt%root_fines_frag(:,:)) + &
                        sum(litt%ag_cwd_frag(:)) + &
                        sum(litt%bg_cwd_frag(:,:)) + &
                        sum(litt%seed_decay(:)) + &
                        sum(litt%seed_germ_decay(:))) * cpatch%area / m2_per_ha / sec_per_day

                   hio_seed_bank_elem(io_si,el) = hio_seed_bank_elem(io_si,el) + &
                        sum(litt%seed(:)) * cpatch%area / m2_per_ha

                   hio_seed_germ_elem(io_si,el) = hio_seed_germ_elem(io_si,el) + &
                        sum(litt%seed_germ(:)) *  cpatch%area / m2_per_ha

                   hio_seed_decay_elem(io_si,el) = hio_seed_decay_elem(io_si,el) + &
                        sum(litt%seed_decay(:) + litt%seed_germ_decay(:) ) *       &
                        cpatch%area / m2_per_ha / sec_per_day

                   hio_seeds_in_local_elem(io_si,el) = hio_seeds_in_local_elem(io_si,el) + &
                        sum(litt%seed_in_local(:)) *  cpatch%area / m2_per_ha / sec_per_day

                   hio_seed_in_extern_elem(io_si,el) = hio_seed_in_extern_elem(io_si,el) + &
                        sum(litt%seed_in_extern(:)) * cpatch%area / m2_per_ha / sec_per_day

                   ! Litter State Variables
                   hio_cwd_ag_elem(io_si,el) = hio_cwd_ag_elem(io_si,el) + &
                        sum(litt%ag_cwd(:)) * cpatch%area / m2_per_ha

                   hio_cwd_bg_elem(io_si,el) = hio_cwd_bg_elem(io_si,el) + &
                        sum(litt%bg_cwd(:,:)) * cpatch%area / m2_per_ha

                   hio_fines_ag_elem(io_si,el) = hio_fines_ag_elem(io_si,el) + &
                        sum(litt%leaf_fines(:)) * cpatch%area / m2_per_ha

                   hio_fines_bg_elem(io_si,el) = hio_fines_bg_elem(io_si,el) + &
                        sum(litt%root_fines(:,:)) * cpatch%area / m2_per_ha

                   do i_cwd=1,ncwd
                      elcwd = (el-1)*ncwd+i_cwd
                      hio_cwd_elcwd(io_si,elcwd) = hio_cwd_elcwd(io_si,elcwd) +   &
                           (litt%ag_cwd(i_cwd) + sum(litt%bg_cwd(i_cwd,:))) *        &
                           cpatch%area / m2_per_ha

                   end do

                   ! Load Mass States
                   ccohort => cpatch%tallest
                   do while(associated(ccohort))

                      sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
                      struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
                      leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
                      fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
                      store_m  = ccohort%prt%GetState(store_organ, element_list(el))
                      repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
                      total_m  = sapw_m+struct_m+leaf_m+fnrt_m+store_m+repro_m


                      i_scpf = ccohort%size_by_pft_class

                      if(element_list(el).eq.carbon12_element)then
                         this%hvars(ih_totvegc_scpf)%r82d(io_si,i_scpf) =             &
                              this%hvars(ih_totvegc_scpf)%r82d(io_si,i_scpf) +          &
                              total_m * ccohort%n / m2_per_ha
                         this%hvars(ih_leafc_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_leafc_scpf)%r82d(io_si,i_scpf) +            &
                              leaf_m * ccohort%n / m2_per_ha
                         this%hvars(ih_fnrtc_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_fnrtc_scpf)%r82d(io_si,i_scpf) +            &
                              fnrt_m * ccohort%n / m2_per_ha
                         this%hvars(ih_sapwc_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_sapwc_scpf)%r82d(io_si,i_scpf) +            &
                              sapw_m * ccohort%n / m2_per_ha
                         this%hvars(ih_storec_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_storec_scpf)%r82d(io_si,i_scpf) +           &
                              store_m * ccohort%n / m2_per_ha
                         this%hvars(ih_reproc_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_reproc_scpf)%r82d(io_si,i_scpf) +           &
                              repro_m * ccohort%n / m2_per_ha
                      elseif(element_list(el).eq.nitrogen_element)then

                         store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                         this%hvars(ih_totvegn_scpf)%r82d(io_si,i_scpf) =             &
                              this%hvars(ih_totvegn_scpf)%r82d(io_si,i_scpf) +          &
                              total_m * ccohort%n / m2_per_ha
                         this%hvars(ih_leafn_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_leafn_scpf)%r82d(io_si,i_scpf) +            &
                              leaf_m * ccohort%n / m2_per_ha
                         this%hvars(ih_fnrtn_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_fnrtn_scpf)%r82d(io_si,i_scpf) +            &
                              fnrt_m * ccohort%n / m2_per_ha
                         this%hvars(ih_sapwn_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_sapwn_scpf)%r82d(io_si,i_scpf) +            &
                              sapw_m * ccohort%n / m2_per_ha
                         this%hvars(ih_storen_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_storen_scpf)%r82d(io_si,i_scpf) +           &
                              store_m * ccohort%n / m2_per_ha
                         this%hvars(ih_repron_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_repron_scpf)%r82d(io_si,i_scpf) +           &
                              repro_m * ccohort%n / m2_per_ha

                      elseif(element_list(el).eq.phosphorus_element)then

                         store_max = ccohort%prt%GetNutrientTarget(element_list(el),store_organ,stoich_growth_min)

                         this%hvars(ih_totvegp_scpf)%r82d(io_si,i_scpf) =             &
                              this%hvars(ih_totvegp_scpf)%r82d(io_si,i_scpf) +          &
                              total_m * ccohort%n / m2_per_ha
                         this%hvars(ih_leafp_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_leafp_scpf)%r82d(io_si,i_scpf) +            &
                              leaf_m * ccohort%n / m2_per_ha
                         this%hvars(ih_fnrtp_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_fnrtp_scpf)%r82d(io_si,i_scpf) +            &
                              fnrt_m * ccohort%n / m2_per_ha
                         this%hvars(ih_sapwp_scpf)%r82d(io_si,i_scpf) =               &
                              this%hvars(ih_sapwp_scpf)%r82d(io_si,i_scpf) +            &
                              sapw_m * ccohort%n / m2_per_ha
                         this%hvars(ih_storep_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_storep_scpf)%r82d(io_si,i_scpf) +           &
                              store_m * ccohort%n / m2_per_ha
                         this%hvars(ih_reprop_scpf)%r82d(io_si,i_scpf) =              &
                              this%hvars(ih_reprop_scpf)%r82d(io_si,i_scpf) +           &
                              repro_m * ccohort%n / m2_per_ha

                      end if

                      ccohort => ccohort%shorter
                   end do ! end cohort loop

                   cpatch => cpatch%younger
                end do ! end patch loop

             end do ! end element loop

             do ft = 1, numpft
                do i_scls = 1,nlevsclass
                   i_scpf = (ft-1)*nlevsclass + i_scls

                   if( this%hvars(ih_storectfrac_canopy_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                      this%hvars(ih_storectfrac_canopy_scpf)%r82d(io_si,i_scpf) = &
                           storec_canopy_scpf(i_scpf) / &
                           this%hvars(ih_storectfrac_canopy_scpf)%r82d(io_si,i_scpf)
                   end if
                   if( this%hvars(ih_storectfrac_ustory_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                      this%hvars(ih_storectfrac_ustory_scpf)%r82d(io_si,i_scpf) = &
                           storec_understory_scpf(i_scpf) / &
                           this%hvars(ih_storectfrac_ustory_scpf)%r82d(io_si,i_scpf)
                   end if

                end do
             end do

             do el = 1, num_elements

                if(element_list(el).eq.nitrogen_element)then

                   do ft = 1, numpft
                      do i_scls = 1,nlevsclass
                         i_scpf = (ft-1)*nlevsclass + i_scls

                         if( this%hvars(ih_storentfrac_canopy_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                            this%hvars(ih_storentfrac_canopy_scpf)%r82d(io_si,i_scpf) = &
                                 storen_canopy_scpf(i_scpf) / &
                                 this%hvars(ih_storentfrac_canopy_scpf)%r82d(io_si,i_scpf)
                         end if
                         if( this%hvars(ih_storentfrac_understory_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                            this%hvars(ih_storentfrac_understory_scpf)%r82d(io_si,i_scpf) = &
                                 storen_understory_scpf(i_scpf) / &
                                 this%hvars(ih_storentfrac_understory_scpf)%r82d(io_si,i_scpf)
                         end if

                      end do
                   end do
                elseif(element_list(el).eq.phosphorus_element)then

                   do ft = 1, numpft
                      do i_scls = 1,nlevsclass
                         i_scpf = (ft-1)*nlevsclass + i_scls

                         if( this%hvars(ih_storeptfrac_canopy_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                            this%hvars(ih_storeptfrac_canopy_scpf)%r82d(io_si,i_scpf) = &
                                 storep_canopy_scpf(i_scpf) / &
                                 this%hvars(ih_storeptfrac_canopy_scpf)%r82d(io_si,i_scpf)
                         end if
                         if( this%hvars(ih_storeptfrac_understory_scpf)%r82d(io_si,i_scpf)>nearzero ) then
                            this%hvars(ih_storeptfrac_understory_scpf)%r82d(io_si,i_scpf) = &
                                 storep_understory_scpf(i_scpf) / &
                                 this%hvars(ih_storeptfrac_understory_scpf)%r82d(io_si,i_scpf)
                         end if

                      end do
                   end do
                end if
             end do

             ! pass demotion rates and associated carbon fluxes to history
             do i_scls = 1,nlevsclass
                hio_demotion_rate_si_scls(io_si,i_scls) = sites(s)%demotion_rate(i_scls) * days_per_year / m2_per_ha
                hio_promotion_rate_si_scls(io_si,i_scls) = sites(s)%promotion_rate(i_scls) * days_per_year / m2_per_ha
             end do

             ! add the site-level disturbance-associated cwd and litter input fluxes to thir respective flux fields

             do i_cwd = 1, ncwd
                hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) + &
                     elflux_diags_c%cwd_ag_input(i_cwd) / days_per_year / sec_per_day

                hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) + &
                     elflux_diags_c%cwd_bg_input(i_cwd) / days_per_year / sec_per_day

             end do
             
             do ft = 1,numpft
                this%hvars(ih_recl2fr_canopy_pf)%r82d(io_si,ft) = sites(s)%rec_l2fr(ft,1)
                this%hvars(ih_recl2fr_ustory_pf)%r82d(io_si,ft) = sites(s)%rec_l2fr(ft,2)
             end do
             
          enddo siteloop ! site loop

        end associate
      end associate
    end associate

    return
  end subroutine update_history_dyn2

  ! ===============================================================================================

  subroutine update_history_hifrq(this,nc,nsites,sites,bc_in,bc_out,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! at the model time-step frequency.
    ! This is the general routine that will call the single or multi-dimensional
    ! routines if they are called for by the user
    ! ---------------------------------------------------------------------------------
    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    type(bc_out_type)       , intent(in)            :: bc_out(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    if(hlm_hist_level_hifrq>0) then
       call update_history_hifrq1(this,nc,nsites,sites,bc_in,bc_out,dt_tstep)
       if(hlm_hist_level_hifrq>1) then
          call update_history_hifrq2(this,nc,nsites,sites,bc_in,bc_out,dt_tstep)
       end if
    end if


    return
  end subroutine update_history_hifrq

  subroutine update_history_hifrq1(this,nc,nsites,sites,bc_in,bc_out,dt_tstep)

    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    type(bc_out_type)       , intent(in)            :: bc_out(nsites)
    real(r8)                , intent(in)            :: dt_tstep

    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: age_class  ! class age index
    real(r8) :: site_area_veg_inv  ! inverse canopy area of the site (1/m2)
    real(r8) :: site_area_rad_inv   ! inverse canopy area of site for only
    ! patches that called the solver
    real(r8) :: dt_tstep_inv        ! inverse timestep (1/sec)
    real(r8) :: n_perm2             ! number of plants per square meter
    real(r8) :: sum_area_rad        ! sum of patch canopy areas
    real(r8),allocatable :: age_area_rad(:)

    type(fates_patch_type),pointer  :: cpatch
    type(fates_cohort_type),pointer :: ccohort


    associate( hio_gpp_si                   => this%hvars(ih_gpp_si)%r81d, &
         hio_npp_si                   => this%hvars(ih_npp_si)%r81d, &
         hio_aresp_si                 => this%hvars(ih_aresp_si)%r81d, &
         hio_maint_resp_si            => this%hvars(ih_maint_resp_si)%r81d, &
         hio_growth_resp_si           => this%hvars(ih_growth_resp_si)%r81d, &
         hio_c_stomata_si             => this%hvars(ih_c_stomata_si)%r81d, &
         hio_c_lblayer_si             => this%hvars(ih_c_lblayer_si)%r81d, &
         hio_vis_rad_err_si           => this%hvars(ih_vis_rad_err_si)%r81d, &
         hio_nir_rad_err_si           => this%hvars(ih_nir_rad_err_si)%r81d, &
         hio_nep_si                   => this%hvars(ih_nep_si)%r81d, &
         hio_hr_si                    => this%hvars(ih_hr_si)%r81d, &
         hio_gpp_canopy_si            => this%hvars(ih_gpp_canopy_si)%r81d, &
         hio_ar_canopy_si             => this%hvars(ih_ar_canopy_si)%r81d, &
         hio_gpp_understory_si        => this%hvars(ih_gpp_understory_si)%r81d, &
         hio_ar_understory_si         => this%hvars(ih_ar_understory_si)%r81d, &
         hio_leaf_mr_si               => this%hvars(ih_leaf_mr_si)%r81d, &
         hio_froot_mr_si              => this%hvars(ih_froot_mr_si)%r81d, &
         hio_livecroot_mr_si          => this%hvars(ih_livecroot_mr_si)%r81d, &
         hio_livestem_mr_si           => this%hvars(ih_livestem_mr_si)%r81d, &
         hio_maint_resp_unreduced_si  => this%hvars(ih_maint_resp_unreduced_si)%r81d, &
         hio_tveg                     => this%hvars(ih_tveg_si)%r81d)


      ! THIS CAN BE REMOVED WHEN BOTH CTSM AND E3SM CALL FLUSH_ALL_HVARS
      ! THIS IS NOT A LIABILITY, IT IS JUST REDUNDANT
      call this%flush_hvars(nc,upfreq_in=group_hifr_simple)
      
      dt_tstep_inv = 1.0_r8/dt_tstep

      allocate(age_area_rad(size(ED_val_history_ageclass_bin_edges,1)+1))

      do_sites: do s = 1,nsites

         call this%zero_site_hvars(sites(s), upfreq_in=group_hifr_simple)

         io_si  = sites(s)%h_gid

         hio_nep_si(io_si) = -bc_in(s)%tot_het_resp * kg_per_g
         hio_hr_si(io_si)  =  bc_in(s)%tot_het_resp * kg_per_g

         ! Diagnostics that are only incremented if we called the radiation solver
         ! We do not call the radiation solver if
         ! a) there is no vegetation
         ! b) there is no light! (ie cos(zenith) ~= 0)
         age_area_rad(:) = 0._r8
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            ! We initialize the solver error to the ignore value
            ! in the radiation driver. It is only modified if the
            ! solver was called. The solver will be called for NIR
            ! if VIS is called, and likewise the same for conservation
            ! error. So the check on VIS solve error will catch all.
            if( abs(cpatch%rad_error(ivis))>nearzero ) then
               age_class = get_age_class_index(cpatch%age)
               age_area_rad(age_class) = age_area_rad(age_class) + cpatch%total_canopy_area
            end if
            cpatch => cpatch%younger
         end do

         sum_area_rad = sum(age_area_rad(:))

         if_anyrad: if(sum_area_rad<nearzero)then
            hio_vis_rad_err_si(io_si)          = hlm_hio_ignore_val
            hio_nir_rad_err_si(io_si)          = hlm_hio_ignore_val
         else

            hio_vis_rad_err_si(io_si)          = 0._r8
            hio_nir_rad_err_si(io_si)          = 0._r8

            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))
               if( abs(cpatch%rad_error(ivis))>nearzero ) then
                  age_class = get_age_class_index(cpatch%age)

                  hio_vis_rad_err_si(io_si) = hio_vis_rad_err_si(io_si) + &
                       cpatch%rad_error(ivis)*cpatch%total_canopy_area/sum_area_rad
                  hio_nir_rad_err_si(io_si) = hio_nir_rad_err_si(io_si) + &
                       cpatch%rad_error(inir)*cpatch%total_canopy_area/sum_area_rad

               end if
               cpatch => cpatch%younger
            end do
         end if if_anyrad

         ! Diagnostics that are only relevant if there is vegetation present on this site
         ! ie, non-zero canopy area

         if (hlm_use_nocomp .eq. itrue .and. hlm_use_fixed_biogeog .eq. itrue) then
            site_area_veg_inv = area - sites(s)%area_bareground * area
         else
            site_area_veg_inv = 0._r8
            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))
               site_area_veg_inv = site_area_veg_inv + cpatch%total_canopy_area
               cpatch => cpatch%younger
            end do !patch loop
         end if
         
         if_veg_area: if(site_area_veg_inv < nearzero) then

            hio_c_stomata_si(io_si) = hlm_hio_ignore_val
            hio_c_lblayer_si(io_si) = hlm_hio_ignore_val
            hio_tveg(io_si)         = hlm_hio_ignore_val

            exit if_veg_area

         else

            site_area_veg_inv = 1._r8/site_area_veg_inv

            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))

               hio_c_stomata_si(io_si) = hio_c_stomata_si(io_si) + &
                    cpatch%c_stomata * cpatch%total_canopy_area * mol_per_umol * site_area_veg_inv

               hio_c_lblayer_si(io_si) = hio_c_lblayer_si(io_si) + &
                    cpatch%c_lblayer * cpatch%total_canopy_area * mol_per_umol * site_area_veg_inv

               ! Only accumulate the instantaneous vegetation temperature for vegetated patches
               if (cpatch%nocomp_pft_label.ne.nocomp_bareground)then
                  hio_tveg(io_si) = hio_tveg(io_si) + &
                       (bc_in(s)%t_veg_pa(cpatch%patchno) - t_water_freeze_k_1atm) * &
                       cpatch%total_canopy_area * site_area_veg_inv
               end if

               ccohort => cpatch%shortest
               do while(associated(ccohort))

                  n_perm2   = ccohort%n * AREA_INV

                  if_notnew: if ( .not. ccohort%isnew ) then

                     ! scale up cohort fluxes to the site level
                     ! these fluxes have conversions of [kg/plant/timestep] -> [kg/m2/s]
                     
                     ! Net Ecosystem Production [kgC/m2/s]. Use yesterday's growth respiration
                     hio_nep_si(io_si) = hio_nep_si(io_si) + &
                          (ccohort%gpp_tstep-ccohort%resp_m_tstep) * n_perm2 * dt_tstep_inv - &
                          (ccohort%resp_g_acc_hold+ccohort%resp_excess_hold) * n_perm2 / days_per_year / sec_per_day

                     hio_gpp_si(io_si) = hio_gpp_si(io_si) + &
                          ccohort%gpp_tstep * n_perm2 * dt_tstep_inv

                     hio_maint_resp_si(io_si) = hio_maint_resp_si(io_si) + &
                          ccohort%resp_m_tstep * n_perm2 * dt_tstep_inv

                     hio_maint_resp_unreduced_si(io_si) = hio_maint_resp_unreduced_si(io_si) + &
                          ccohort%resp_m_unreduced * n_perm2 * dt_tstep_inv

                     ! Maintenance respiration of different organs
                     hio_leaf_mr_si(io_si) = hio_leaf_mr_si(io_si) + ccohort%rdark &
                          * n_perm2
                     hio_froot_mr_si(io_si) = hio_froot_mr_si(io_si) + ccohort%froot_mr &
                          * n_perm2
                     hio_livecroot_mr_si(io_si) = hio_livecroot_mr_si(io_si) + ccohort%livecroot_mr &
                          * n_perm2
                     hio_livestem_mr_si(io_si) = hio_livestem_mr_si(io_si) + ccohort%livestem_mr &
                          * n_perm2

                     ! accumulate fluxes on canopy- and understory- separated fluxes
                     ! these fluxes have conversions of [kg/plant/timestep] -> [kg/m2/s]
                     if (ccohort%canopy_layer .eq. 1) then

                        hio_gpp_canopy_si(io_si) = hio_gpp_canopy_si(io_si) + &
                             ccohort%gpp_tstep * n_perm2 * dt_tstep_inv

                     else

                        hio_gpp_understory_si(io_si) = hio_gpp_understory_si(io_si) + &
                             ccohort%gpp_tstep * n_perm2 * dt_tstep_inv

                     end if

                  end if if_notnew
                  ccohort => ccohort%taller
               end do
               cpatch => cpatch%younger
            end do
         end if if_veg_area
      end do do_sites

      deallocate(age_area_rad)

    end associate
    return
  end subroutine update_history_hifrq1

  ! ===============================================================================================

  subroutine update_history_hifrq2(this,nc,nsites,sites,bc_in,bc_out,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays for multi-dimension arrays
    ! that change rapidly.  This is an expensive call, the model will probably run
    ! much faster if the user is not using any of these diagnostics.
    ! ---------------------------------------------------------------------------------

    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    type(bc_out_type)       , intent(in)            :: bc_out(nsites)
    real(r8)                , intent(in)            :: dt_tstep

    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_area_by_age(nlevage)  ! patch area in each bin for normalizing purposes
    real(r8) :: canopy_area_by_age(nlevage) ! canopy area in each bin for normalizing purposes
    real(r8) :: site_area_veg_inv           ! 1/area of the site that is not bare-ground 
    integer  :: ipa2     ! patch incrementer
    integer  :: clllpf_indx, cnlf_indx, ipft, ican, ileaf ! more iterators and indices
    real(r8) :: clllpf_area  ! area footprint (m2) for the current cl x ll x pft bin
    real(r8) :: clll_area    ! area footprint (m2) for the cl x ll bin (ie adds up pfts in parallel)
    real(r8) :: cl_area      ! total weight of all ll x pft bins in the canopy layer

    type(fates_patch_type),pointer  :: cpatch
    type(fates_cohort_type),pointer :: ccohort
    real(r8) :: dt_tstep_inv          ! Time step in frequency units (/s)

    associate( hio_ar_si_scpf                      => this%hvars(ih_ar_si_scpf)%r82d, &
         hio_ar_grow_si_scpf                 => this%hvars(ih_ar_grow_si_scpf)%r82d, &
         hio_ar_maint_si_scpf                => this%hvars(ih_ar_maint_si_scpf)%r82d, &
         hio_ar_agsapm_si_scpf               => this%hvars(ih_ar_agsapm_si_scpf)%r82d, &
         hio_ar_darkm_si_scpf                => this%hvars(ih_ar_darkm_si_scpf)%r82d, &
         hio_ar_crootm_si_scpf               => this%hvars(ih_ar_crootm_si_scpf)%r82d, &
         hio_ar_frootm_si_scpf               => this%hvars(ih_ar_frootm_si_scpf)%r82d, &
         hio_rdark_canopy_si_scls            => this%hvars(ih_rdark_canopy_si_scls)%r82d, &
         hio_livestem_mr_canopy_si_scls      => this%hvars(ih_livestem_mr_canopy_si_scls)%r82d, &
         hio_livecroot_mr_canopy_si_scls     => this%hvars(ih_livecroot_mr_canopy_si_scls)%r82d, &
         hio_froot_mr_canopy_si_scls         => this%hvars(ih_froot_mr_canopy_si_scls)%r82d, &
         hio_resp_g_canopy_si_scls           => this%hvars(ih_resp_g_canopy_si_scls)%r82d, &
         hio_resp_m_canopy_si_scls           => this%hvars(ih_resp_m_canopy_si_scls)%r82d, &
         hio_rdark_understory_si_scls        => this%hvars(ih_rdark_understory_si_scls)%r82d, &
         hio_livestem_mr_understory_si_scls  => this%hvars(ih_livestem_mr_understory_si_scls)%r82d, &
         hio_livecroot_mr_understory_si_scls => this%hvars(ih_livecroot_mr_understory_si_scls)%r82d, &
         hio_froot_mr_understory_si_scls     => this%hvars(ih_froot_mr_understory_si_scls)%r82d, &
         hio_resp_g_understory_si_scls       => this%hvars(ih_resp_g_understory_si_scls)%r82d, &
         hio_resp_m_understory_si_scls       => this%hvars(ih_resp_m_understory_si_scls)%r82d, &
         hio_gpp_si_age                      => this%hvars(ih_gpp_si_age)%r82d, &
         hio_gpp_si_landuse                  => this%hvars(ih_gpp_si_landuse)%r82d, &
         hio_c_stomata_si_age                => this%hvars(ih_c_stomata_si_age)%r82d, &
         hio_c_lblayer_si_age                => this%hvars(ih_c_lblayer_si_age)%r82d, &
         hio_parsun_z_si_cnlf                => this%hvars(ih_parsun_z_si_cnlf)%r82d, &
         hio_parsha_z_si_cnlf                => this%hvars(ih_parsha_z_si_cnlf)%r82d, &
         hio_ts_net_uptake_si_cnlf           => this%hvars(ih_ts_net_uptake_si_cnlf)%r82d, &
         hio_parsun_z_si_cnlfpft             => this%hvars(ih_parsun_z_si_cnlfpft)%r82d, &
         hio_parsha_z_si_cnlfpft             => this%hvars(ih_parsha_z_si_cnlfpft)%r82d, &
         hio_laisun_z_si_cnlf                => this%hvars(ih_laisun_z_si_cnlf)%r82d, &
         hio_laisha_z_si_cnlf                => this%hvars(ih_laisha_z_si_cnlf)%r82d, &
         hio_laisun_clllpf                   => this%hvars(ih_laisun_clllpf)%r82d, &
         hio_laisha_clllpf                   => this%hvars(ih_laisha_clllpf)%r82d, &
         hio_crownfrac_clllpf                => this%hvars(ih_crownfrac_clllpf)%r82d, &
         hio_parprof_dir_si_cnlf             => this%hvars(ih_parprof_dir_si_cnlf)%r82d, &
         hio_parprof_dif_si_cnlf             => this%hvars(ih_parprof_dif_si_cnlf)%r82d, &
         hio_parprof_dir_si_cnlfpft          => this%hvars(ih_parprof_dir_si_cnlfpft)%r82d, &
         hio_parprof_dif_si_cnlfpft          => this%hvars(ih_parprof_dif_si_cnlfpft)%r82d, &
         hio_parsun_si_can                   => this%hvars(ih_parsun_si_can)%r82d, &
         hio_parsha_si_can                   => this%hvars(ih_parsha_si_can)%r82d, &
         hio_laisun_si_can                    => this%hvars(ih_laisun_si_can)%r82d, &
         hio_laisha_si_can                    => this%hvars(ih_laisha_si_can)%r82d )


      ! THIS CAN BE REMOVED WHEN BOTH CTSM AND E3SM CALL FLUSH_ALL_HVARS
      ! THIS IS NOT A LIABILITY, IT IS JUST REDUNDANT 
      call this%flush_hvars(nc,upfreq_in=group_hifr_complx)
      
      dt_tstep_inv = 1.0_r8/dt_tstep

      do_sites: do s = 1,nsites

         call this%zero_site_hvars(sites(s), upfreq_in=group_hifr_complx)
         
         site_area_veg_inv = 0._r8
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            site_area_veg_inv = site_area_veg_inv + cpatch%total_canopy_area
            cpatch => cpatch%younger
         end do !patch loop

         ! If there is no vegetation, go to the next site
         if(site_area_veg_inv < nearzero) cycle do_sites

         site_area_veg_inv = 1._r8/site_area_veg_inv

         io_si  = sites(s)%h_gid

         patch_area_by_age(1:nlevage) = 0._r8
         canopy_area_by_age(1:nlevage) = 0._r8

         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))

            patch_area_by_age(cpatch%age_class)  = &
                 patch_area_by_age(cpatch%age_class) + cpatch%area

            canopy_area_by_age(cpatch%age_class) = &
                 canopy_area_by_age(cpatch%age_class) + cpatch%total_canopy_area

            ! Canopy resitance terms
            hio_c_stomata_si_age(io_si,cpatch%age_class) = &
                 hio_c_stomata_si_age(io_si,cpatch%age_class) + &
                 cpatch%c_stomata * cpatch%total_canopy_area * mol_per_umol

            hio_c_lblayer_si_age(io_si,cpatch%age_class) = &
                 hio_c_lblayer_si_age(io_si,cpatch%age_class) + &
                 cpatch%c_lblayer * cpatch%total_canopy_area * mol_per_umol

            ccohort => cpatch%shortest
            do while(associated(ccohort))

               n_perm2   = ccohort%n * AREA_INV

               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  associate( scpf => ccohort%size_by_pft_class, &
                       scls => ccohort%size_class )

                    ! Total AR (kgC/m2/s) = (kgC/plant/step) / (s/step) * (plant/m2)
                    hio_ar_si_scpf(io_si,scpf)    =   hio_ar_si_scpf(io_si,scpf) + &
                         (ccohort%resp_m_tstep*dt_tstep_inv) * n_perm2 + &
                         (ccohort%resp_g_acc_hold + ccohort%resp_excess_hold)* n_perm2 / days_per_year / sec_per_day

                    ! Growth AR (kgC/m2/s)   ! CDK: this should be daily
                    hio_ar_grow_si_scpf(io_si,scpf) = hio_ar_grow_si_scpf(io_si,scpf) + &
                         ccohort%resp_g_acc_hold * n_perm2 / days_per_year / sec_per_day

                    ! Maint AR (kgC/m2/s)
                    hio_ar_maint_si_scpf(io_si,scpf) = hio_ar_maint_si_scpf(io_si,scpf) + &
                         (ccohort%resp_m_tstep*dt_tstep_inv) * n_perm2

                    ! Maintenance AR partition variables are stored as rates (kgC/plant/s)
                    ! (kgC/m2/s) = (kgC/plant/s) * (plant/m2)
                    hio_ar_agsapm_si_scpf(io_si,scpf) = hio_ar_agsapm_si_scpf(io_si,scpf) + &
                         ccohort%livestem_mr * n_perm2

                    ! (kgC/m2/s) = (kgC/plant/s) * (plant/m2)
                    hio_ar_darkm_si_scpf(io_si,scpf) = hio_ar_darkm_si_scpf(io_si,scpf) + &
                         ccohort%rdark * n_perm2

                    ! (kgC/m2/s) = (kgC/plant/s) * (plant/m2)
                    hio_ar_crootm_si_scpf(io_si,scpf) = hio_ar_crootm_si_scpf(io_si,scpf) + &
                         ccohort%livecroot_mr * n_perm2

                    ! (kgC/m2/s) = (kgC/plant/s) * (plant/m2)
                    hio_ar_frootm_si_scpf(io_si,scpf) = hio_ar_frootm_si_scpf(io_si,scpf) + &
                         ccohort%froot_mr * n_perm2

                    ! accumulate fluxes per patch age bin
                    hio_gpp_si_age(io_si,cpatch%age_class) = hio_gpp_si_age(io_si,cpatch%age_class) &
                         + ccohort%gpp_tstep * ccohort%n * dt_tstep_inv

                    if (cpatch%land_use_label .gt. nocomp_bareground_land) then
                       hio_gpp_si_landuse(io_si,cpatch%land_use_label) = hio_gpp_si_landuse(io_si,cpatch%land_use_label) &
                            + ccohort%gpp_tstep * ccohort%n * dt_tstep_inv
                    end if

                    ! accumulate fluxes on canopy- and understory- separated fluxes
                    if (ccohort%canopy_layer .eq. 1) then

                       ! size-resolved respiration fluxes are in kg C / m2 / s
                       hio_rdark_canopy_si_scls(io_si,scls) = hio_rdark_canopy_si_scls(io_si,scls) + &
                            ccohort%rdark  * ccohort%n * ha_per_m2
                       hio_livestem_mr_canopy_si_scls(io_si,scls) = hio_livestem_mr_canopy_si_scls(io_si,scls) + &
                            ccohort%livestem_mr  * ccohort%n * ha_per_m2
                       hio_livecroot_mr_canopy_si_scls(io_si,scls) = hio_livecroot_mr_canopy_si_scls(io_si,scls) + &
                            ccohort%livecroot_mr  * ccohort%n * ha_per_m2
                       hio_froot_mr_canopy_si_scls(io_si,scls) = hio_froot_mr_canopy_si_scls(io_si,scls) + &
                            ccohort%froot_mr  * ccohort%n * ha_per_m2
                       hio_resp_g_canopy_si_scls(io_si,scls) = hio_resp_g_canopy_si_scls(io_si,scls) + &
                            ccohort%resp_g_acc_hold * n_perm2 / days_per_year / sec_per_day
                       hio_resp_m_canopy_si_scls(io_si,scls) = hio_resp_m_canopy_si_scls(io_si,scls) + &
                            ccohort%resp_m_tstep  * ccohort%n * dt_tstep_inv * ha_per_m2
                    else

                       ! size-resolved respiration fluxes are in kg C / m2 / s
                       hio_rdark_understory_si_scls(io_si,scls) = hio_rdark_understory_si_scls(io_si,scls) + &
                            ccohort%rdark  * ccohort%n * ha_per_m2
                       hio_livestem_mr_understory_si_scls(io_si,scls) = hio_livestem_mr_understory_si_scls(io_si,scls) + &
                            ccohort%livestem_mr  * ccohort%n  * ha_per_m2
                       hio_livecroot_mr_understory_si_scls(io_si,scls) = hio_livecroot_mr_understory_si_scls(io_si,scls) + &
                            ccohort%livecroot_mr  * ccohort%n  * ha_per_m2
                       hio_froot_mr_understory_si_scls(io_si,scls) = hio_froot_mr_understory_si_scls(io_si,scls) + &
                            ccohort%froot_mr  * ccohort%n  * ha_per_m2
                       hio_resp_g_understory_si_scls(io_si,scls) = hio_resp_g_understory_si_scls(io_si,scls) + &
                            ccohort%resp_g_acc_hold * n_perm2 / days_per_year / sec_per_day
                       hio_resp_m_understory_si_scls(io_si,scls) = hio_resp_m_understory_si_scls(io_si,scls) + &
                            ccohort%resp_m_tstep  * ccohort%n * dt_tstep_inv  * ha_per_m2
                    endif
                  end associate
               endif

!!! canopy leaf carbon balance
               ican = ccohort%canopy_layer
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) = hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%ts_net_uptake(ileaf) * dt_tstep_inv * ccohort%c_area * area_inv
               end do

               ccohort => ccohort%taller
            enddo ! cohort loop


            ! summarize radiation profiles through the canopy
            ! --------------------------------------------------------------------

            do_pft1: do ipft=1,numpft
               do_canlev1: do ican=1,cpatch%ncl_p
                  do_leaflev1: do ileaf=1,cpatch%nleaf(ican,ipft)

                     ! calculate where we are on multiplexed dimensions
                     clllpf_indx = ileaf + (ican-1) * nlevleaf + (ipft-1) * nlevleaf * nclmax
                     cnlf_indx = ileaf + (ican-1) * nlevleaf

                     ! canopy_area_profile is the fraction of the total canopy area that
                     ! is occupied by this bin.  If you add up the top leaf layer bins in the
                     ! top canopy layers, for all pfts, that should equal to 1

                     clllpf_area = cpatch%canopy_area_profile(ican,ipft,ileaf)*cpatch%total_canopy_area

                     ! Canopy by leaf by pft level diagnostics
                     ! -------------------------------------------------------------------
                     hio_parsun_z_si_cnlfpft(io_si,clllpf_indx) = hio_parsun_z_si_cnlfpft(io_si,clllpf_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * clllpf_area

                     hio_parsha_z_si_cnlfpft(io_si,clllpf_indx) = hio_parsha_z_si_cnlfpft(io_si,clllpf_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * clllpf_area

                     ! elai_profile is the m2 of leaf inside the m2 of bin.

                     hio_laisun_clllpf(io_si, clllpf_indx) = hio_laisun_clllpf(io_si, clllpf_indx) + &
                          cpatch%elai_profile(ican,ipft,ileaf)*cpatch%f_sun(ican,ipft,ileaf)*clllpf_area

                     hio_laisha_clllpf(io_si,clllpf_indx) = hio_laisha_clllpf(io_si,clllpf_indx) + &
                          cpatch%elai_profile(ican,ipft,ileaf)*(1._r8-cpatch%f_sun(ican,ipft,ileaf))*clllpf_area

                     hio_parprof_dir_si_cnlfpft(io_si,clllpf_indx) = hio_parprof_dir_si_cnlfpft(io_si,clllpf_indx) + &
                          cpatch%parprof_pft_dir_z(ican,ipft,ileaf) * clllpf_area

                     hio_parprof_dif_si_cnlfpft(io_si,clllpf_indx) = hio_parprof_dif_si_cnlfpft(io_si,clllpf_indx) + &
                          cpatch%parprof_pft_dif_z(ican,ipft,ileaf) * clllpf_area

                     ! The fractional area of Canopy layer and PFTs can be used
                     ! do upscale the CLLLPF properties
                     hio_crownfrac_clllpf(io_si,clllpf_indx) = hio_crownfrac_clllpf(io_si,clllpf_indx) + &
                          clllpf_area


                     ! Canopy by leaf layer (mean across pfts) level diagnostics
                     ! ----------------------------------------------------------------------------
                     hio_parprof_dir_si_cnlf(io_si,cnlf_indx) = hio_parprof_dir_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%parprof_pft_dir_z(ican,ipft,ileaf) * clllpf_area

                     hio_parprof_dif_si_cnlf(io_si,cnlf_indx) = hio_parprof_dif_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%parprof_pft_dif_z(ican,ipft,ileaf) * clllpf_area

                     hio_parsun_z_si_cnlf(io_si,cnlf_indx) = hio_parsun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * clllpf_area

                     hio_parsha_z_si_cnlf(io_si,cnlf_indx) = hio_parsha_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * clllpf_area

                     hio_laisun_z_si_cnlf(io_si,cnlf_indx) = hio_laisun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%f_sun(ican,ipft,ileaf)*clllpf_area

                     hio_laisha_z_si_cnlf(io_si,cnlf_indx) = hio_laisha_z_si_cnlf(io_si,cnlf_indx) + &
                          (1._r8-cpatch%f_sun(ican,ipft,ileaf))*clllpf_area

                     ! Canopy mean diagnostics
                     ! --------------------------------------------------------------

                     hio_parsun_si_can(io_si,ican) = hio_parsun_si_can(io_si,ican) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * clllpf_area
                     hio_parsha_si_can(io_si,ican) = hio_parsha_si_can(io_si,ican) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * clllpf_area

                     hio_laisun_si_can(io_si,ican) = hio_laisun_si_can(io_si,ican) + &
                          cpatch%f_sun(ican,ipft,ileaf)*cpatch%elai_profile(ican,ipft,ileaf) * clllpf_area
                     hio_laisha_si_can(io_si,ican) = hio_laisha_si_can(io_si,ican) + &
                          (1._r8-cpatch%f_sun(ican,ipft,ileaf))*cpatch%elai_profile(ican,ipft,ileaf) * clllpf_area


                  end do do_leaflev1
               end do do_canlev1
            end do do_pft1

            cpatch => cpatch%younger
         end do !patch loop

         ! Normalize the radiation multiplexed diagnostics
         ! Set values that dont have canopy elements to ignore
         ! ----------------------------------------------------------------------------

         do_ican2: do ican = 1,nclmax

            cl_area = 0._r8
            do_ileaf2: do ileaf = 1,nlevleaf

               clll_area = 0._r8
               do_ipft2: do ipft = 1,numpft

                  clllpf_indx = ileaf + (ican-1) * nlevleaf + (ipft-1) * nlevleaf * nclmax
                  if( hio_crownfrac_clllpf(io_si,clllpf_indx)<nearzero)then
                     hio_parsun_z_si_cnlfpft(io_si,clllpf_indx) = hlm_hio_ignore_val
                     hio_parsha_z_si_cnlfpft(io_si,clllpf_indx) = hlm_hio_ignore_val
                     hio_laisun_clllpf(io_si,clllpf_indx) = hlm_hio_ignore_val
                     hio_laisha_clllpf(io_si,clllpf_indx) = hlm_hio_ignore_val
                     hio_parprof_dir_si_cnlfpft(io_si,clllpf_indx) = hlm_hio_ignore_val
                     hio_parprof_dif_si_cnlfpft(io_si,clllpf_indx) = hlm_hio_ignore_val
                  else

                     hio_parsun_z_si_cnlfpft(io_si,clllpf_indx) = &
                          hio_parsun_z_si_cnlfpft(io_si,clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     hio_parsha_z_si_cnlfpft(io_si,clllpf_indx) = &
                          hio_parsha_z_si_cnlfpft(io_si,clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     hio_laisun_clllpf(io_si, clllpf_indx) = &
                          hio_laisun_clllpf(io_si, clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     hio_laisha_clllpf(io_si,clllpf_indx) = &
                          hio_laisha_clllpf(io_si,clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     hio_parprof_dir_si_cnlfpft(io_si,clllpf_indx) = &
                          hio_parprof_dir_si_cnlfpft(io_si,clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     hio_parprof_dif_si_cnlfpft(io_si,clllpf_indx) = &
                          hio_parprof_dif_si_cnlfpft(io_si,clllpf_indx)/hio_crownfrac_clllpf(io_si,clllpf_indx)

                     clll_area = clll_area + hio_crownfrac_clllpf(io_si,clllpf_indx)
                     cl_area   = cl_area + hio_crownfrac_clllpf(io_si,clllpf_indx)

                     ! Convert from total m2 to fraction of the site
                     hio_crownfrac_clllpf(io_si,clllpf_indx) = hio_crownfrac_clllpf(io_si,clllpf_indx)*site_area_veg_inv                     
                  end if
               end do do_ipft2

               cnlf_indx = ileaf + (ican-1) * nlevleaf

               if(clll_area<nearzero)then
                  hio_parprof_dir_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
                  hio_parprof_dif_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
                  hio_parsun_z_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
                  hio_parsha_z_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
                  hio_laisun_z_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
                  hio_laisha_z_si_cnlf(io_si,cnlf_indx) = hlm_hio_ignore_val
               else

                  hio_parprof_dir_si_cnlf(io_si,cnlf_indx) = &
                       hio_parprof_dir_si_cnlf(io_si,cnlf_indx)/clll_area
                  hio_parprof_dif_si_cnlf(io_si,cnlf_indx) = &
                       hio_parprof_dif_si_cnlf(io_si,cnlf_indx)/clll_area
                  hio_parsun_z_si_cnlf(io_si,cnlf_indx) = &
                       hio_parsun_z_si_cnlf(io_si,cnlf_indx)/clll_area
                  hio_parsha_z_si_cnlf(io_si,cnlf_indx) = &
                       hio_parsha_z_si_cnlf(io_si,cnlf_indx)/clll_area
                  hio_laisun_z_si_cnlf(io_si,cnlf_indx) = &
                       hio_laisun_z_si_cnlf(io_si,cnlf_indx)/clll_area
                  hio_laisha_z_si_cnlf(io_si,cnlf_indx) = &
                       hio_laisha_z_si_cnlf(io_si,cnlf_indx)/clll_area
               end if
            end do do_ileaf2

            if(cl_area<nearzero)then

               hio_parsun_si_can(io_si,ican) = hlm_hio_ignore_val
               hio_parsha_si_can(io_si,ican) = hlm_hio_ignore_val
               hio_laisun_si_can(io_si,ican) = hlm_hio_ignore_val
               hio_laisha_si_can(io_si,ican) = hlm_hio_ignore_val

            else
               ! Since these are integrated metrics, ie absorbed over depth
               ! and total leaf over depth, we just want to normalize by the
               ! the area of the footprint.  The weightings they had
               ! recieved were always in m2 (ie the footprint of the bin)
               hio_parsun_si_can(io_si,ican) = hio_parsun_si_can(io_si,ican) * site_area_veg_inv
               hio_parsha_si_can(io_si,ican) = hio_parsha_si_can(io_si,ican) * site_area_veg_inv
               hio_laisun_si_can(io_si,ican) = hio_laisun_si_can(io_si,ican) * site_area_veg_inv
               hio_laisha_si_can(io_si,ican) = hio_laisha_si_can(io_si,ican) * site_area_veg_inv
            end if

         end do do_ican2


         ! Normalize age stratified diagnostics
         ! ----------------------------------------------------------------
         do ipa2 = 1, nlevage
            if (patch_area_by_age(ipa2) .gt. nearzero) then
               hio_gpp_si_age(io_si, ipa2) = &
                    hio_gpp_si_age(io_si, ipa2) / (patch_area_by_age(ipa2))
            else
               hio_gpp_si_age(io_si, ipa2) = 0._r8
            endif

            ! Normalize resistance diagnostics
            if (canopy_area_by_age(ipa2) .gt. nearzero) then
               hio_c_stomata_si_age(io_si,ipa2) = &
                    hio_c_stomata_si_age(io_si,ipa2) / canopy_area_by_age(ipa2)

               hio_c_lblayer_si_age(io_si,ipa2) = &
                    hio_c_lblayer_si_age(io_si,ipa2) / canopy_area_by_age(ipa2)
            else
               hio_c_stomata_si_age(io_si,ipa2) = 0._r8
               hio_c_lblayer_si_age(io_si,ipa2) = 0._r8
            end if

         end do

      enddo do_sites ! site loop

    end associate

  end subroutine update_history_hifrq2

  ! =====================================================================================

  subroutine update_history_hydraulics(this,nc,nsites,sites,bc_in,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------

    use FatesHydraulicsMemMod, only : ed_cohort_hydr_type, nshell
    use FatesHydraulicsMemMod, only : ed_site_hydr_type

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    real(r8)                , intent(in)            :: dt_tstep

    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ft               ! functional type index
    !    integer  :: io_shsl  ! The combined "SH"ell "S"oil "L"ayer index in the IO array
    real(r8) :: ncohort_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
    ! should be "hio_nplant_si_scpf"
    real(r8) :: nplant_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
    ! should be "hio_nplant_si_scpf"
    real(r8) :: number_fraction
    real(r8) :: number_fraction_rate
    real(r8) :: mean_aroot
    integer  :: ipa2     ! patch incrementer
    integer  :: ix       ! histogram x (count) bin index
    integer  :: iscpf    ! index of the scpf group
    integer  :: ipft     ! index of the pft loop
    integer  :: iscls    ! index of the size-class loop
    integer  :: k        ! rhizosphere shell index
    integer  :: j        ! rhizosphere (ie root) layer index
    integer  :: j_bc     ! Soil layer index (ie boundary condition grid index)
    integer  :: j_t,j_b  ! top and bottom soil layer matching current rhiz layer
    integer  :: nlevrhiz ! number of rhizosphere layers
    integer  :: nlevsoil ! number of soil layers
    real(r8) :: mean_soil_vwc    ! mean soil volumetric water content [m3/m3]
    real(r8) :: mean_soil_vwcsat ! mean soil saturated volumetric water content [m3/m3]
    real(r8) :: mean_soil_matpot ! mean soil water potential [MPa]
    real(r8) :: layer_areaweight ! root area weighting factor for each soil layer
    real(r8) :: areaweight       ! root area weighting factor for column
    real(r8) :: vwc              ! volumetric water content of layer [m3/m3] = theta
    real(r8) :: vwc_sat          ! saturated water content of layer [m3/m3]
    real(r8) :: psi              ! matric potential of soil layer
    real(r8) :: depth_frac       ! fraction of rhizosphere layer depth occupied by current soil layer
    character(2) :: fmt_char
    type(fates_patch_type),pointer  :: cpatch
    type(fates_cohort_type),pointer :: ccohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: site_hydr
    real(r8) :: per_dt_tstep          ! Time step in frequency units (/s)
    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! Should this be 365.25?

    if(hlm_use_planthydro.eq.ifalse) return

    per_dt_tstep = 1._r8 / dt_tstep

    if_hifrq0: if(hlm_hist_level_hifrq>0) then

       ! THIS CAN BE REMOVED WHEN BOTH CTSM AND E3SM CALL FLUSH_ALL_HVARS
       ! THIS IS NOT A LIABILITY, IT IS JUST REDUNDANT 
       call this%flush_hvars(nc,upfreq_in=group_hydr_simple)
       
       associate(   hio_h2oveg_hydro_err_si   => this%hvars(ih_h2oveg_hydro_err_si)%r81d, &
            hio_rootwgt_soilvwc_si    => this%hvars(ih_rootwgt_soilvwc_si)%r81d, &
            hio_rootwgt_soilvwcsat_si => this%hvars(ih_rootwgt_soilvwcsat_si)%r81d, &
            hio_rootwgt_soilmatpot_si => this%hvars(ih_rootwgt_soilmatpot_si)%r81d, &
            hio_sapflow_si        => this%hvars(ih_sapflow_si)%r81d, &
            hio_rootuptake_si         => this%hvars(ih_rootuptake_si)%r81d,  &
            hio_h2oveg_si         => this%hvars(ih_h2oveg_si)%r81d )

         do s = 1,nsites

            call this%zero_site_hvars(sites(s),upfreq_in=group_hydr_simple)

            site_hydr => sites(s)%si_hydr
            nlevrhiz = site_hydr%nlevrhiz
            nlevsoil = bc_in(s)%nlevsoil
            io_si  = sites(s)%h_gid

            hio_h2oveg_si(io_si)              = site_hydr%h2oveg
            hio_h2oveg_hydro_err_si(io_si)    = site_hydr%h2oveg_hydro_err
            hio_rootuptake_si(io_si) = sum(site_hydr%rootuptake_sl,dim=1)

            ! Get column means of some soil diagnostics, these are weighted
            ! by the amount of fine-root surface area in each layer
            ! --------------------------------------------------------------------

            mean_soil_vwc    = 0._r8
            mean_soil_matpot = 0._r8
            mean_soil_vwcsat = 0._r8
            areaweight       = 0._r8

            do j=1,nlevrhiz

               j_t = site_hydr%map_r2s(j,1) ! top soil layer matching rhiz layer
               j_b = site_hydr%map_r2s(j,2) ! bottom soil layer matching rhiz layer

               do j_bc = j_t,j_b

                  vwc     = bc_in(s)%h2o_liqvol_sl(j_bc)
                  psi     = site_hydr%wrf_soil(j)%p%psi_from_th(vwc) ! MLO: Any reason for not using smp_sl?
                  ! cap capillary pressure
                  ! psi = max(-1e5_r8,psi) Removing cap as that is inconstistent
                  !                        with model internals and physics. Should
                  !                        implement caps inside the functions
                  !                        if desired. (RGK 12-2021)
                  vwc_sat = bc_in(s)%watsat_sl(j_bc)
                  depth_frac = bc_in(s)%dz_sisl(j_bc)/site_hydr%dz_rhiz(j)

                  ! If there are any roots, we use root weighting
                  if(sum(site_hydr%l_aroot_layer(:),dim=1) > nearzero) then
                     layer_areaweight = site_hydr%l_aroot_layer(j)*depth_frac*pi_const*site_hydr%rs1(j)**2.0

                     ! If there are no roots, we use depth weighting
                  else
                     layer_areaweight = bc_in(s)%dz_sisl(j_bc)
                  endif

                  areaweight       = areaweight + layer_areaweight
                  mean_soil_vwc    = mean_soil_vwc + vwc*layer_areaweight
                  mean_soil_vwcsat = mean_soil_vwcsat + vwc_sat*layer_areaweight
                  mean_soil_matpot = mean_soil_matpot + psi*layer_areaweight

               end do
            end do

            hio_rootwgt_soilvwc_si(io_si)    = mean_soil_vwc/areaweight
            hio_rootwgt_soilvwcsat_si(io_si) = mean_soil_vwcsat/areaweight
            hio_rootwgt_soilmatpot_si(io_si) = mean_soil_matpot/areaweight  * pa_per_mpa

            ! calculate site sapflow
            do ipft = 1, numpft
               do iscls = 1,nlevsclass
                  hio_sapflow_si(io_si) = hio_sapflow_si(io_si) + &
                       site_hydr%sapflow_scpf(iscls, ipft) * ha_per_m2
               end do
            end do
            
            
         end do
       end associate
    end if if_hifrq0

    if_hifrq1: if(hlm_hist_level_hifrq>1) then

       associate( hio_errh2o_scpf  => this%hvars(ih_errh2o_scpf)%r82d, &
            hio_tran_scpf         => this%hvars(ih_tran_scpf)%r82d, &
            hio_sapflow_scpf      => this%hvars(ih_sapflow_scpf)%r82d, &
            hio_iterh1_scpf       => this%hvars(ih_iterh1_scpf)%r82d, &
            hio_iterh2_scpf       => this%hvars(ih_iterh2_scpf)%r82d, &
            hio_ath_scpf          => this%hvars(ih_ath_scpf)%r82d, &
            hio_tth_scpf          => this%hvars(ih_tth_scpf)%r82d, &
            hio_sth_scpf          => this%hvars(ih_sth_scpf)%r82d, &
            hio_lth_scpf          => this%hvars(ih_lth_scpf)%r82d, &
            hio_awp_scpf          => this%hvars(ih_awp_scpf)%r82d, &
            hio_twp_scpf          => this%hvars(ih_twp_scpf)%r82d, &
            hio_swp_scpf          => this%hvars(ih_swp_scpf)%r82d, &
            hio_lwp_scpf          => this%hvars(ih_lwp_scpf)%r82d, &
            hio_aflc_scpf          => this%hvars(ih_aflc_scpf)%r82d, &
            hio_tflc_scpf          => this%hvars(ih_tflc_scpf)%r82d, &
            hio_sflc_scpf          => this%hvars(ih_sflc_scpf)%r82d, &
            hio_lflc_scpf          => this%hvars(ih_lflc_scpf)%r82d, &
            hio_btran_scpf        => this%hvars(ih_btran_scpf)%r82d, &
            
            hio_nplant_si_scpf    => this%hvars(ih_nplant_si_scpf)%r82d, &
            hio_nplant_si_capf    => this%hvars(ih_nplant_si_capf)%r82d, &
            hio_soilmatpot_sl         => this%hvars(ih_soilmatpot_sl)%r82d, &
            hio_soilvwc_sl            => this%hvars(ih_soilvwc_sl)%r82d, &
            hio_soilvwcsat_sl         => this%hvars(ih_soilvwcsat_sl)%r82d, &
            hio_rootuptake_sl         => this%hvars(ih_rootuptake_sl)%r82d, &
            hio_rootuptake0_scpf      => this%hvars(ih_rootuptake0_scpf)%r82d, &
            hio_rootuptake10_scpf     => this%hvars(ih_rootuptake10_scpf)%r82d, &
            hio_rootuptake50_scpf     => this%hvars(ih_rootuptake50_scpf)%r82d, &
            hio_rootuptake100_scpf    => this%hvars(ih_rootuptake100_scpf)%r82d )

         ! THIS CAN BE REMOVED WHEN BOTH CTSM AND E3SM CALL FLUSH_ALL_HVARS
         ! THIS IS NOT A LIABILITY, IT IS JUST REDUNDANT 
         call this%flush_hvars(nc,upfreq_in=group_hydr_complx)
         
         do s = 1,nsites
            
            site_hydr => sites(s)%si_hydr
            nlevrhiz = site_hydr%nlevrhiz
            nlevsoil = bc_in(s)%nlevsoil
            io_si  = sites(s)%h_gid

            call this%zero_site_hvars(sites(s),upfreq_in=group_hydr_complx)
            
            hio_rootuptake_sl(io_si,1:nlevsoil) = site_hydr%rootuptake_sl(1:nlevsoil)

            ! Get column means of some soil diagnostics, these are weighted
            ! by the amount of fine-root surface area in each layer
            ! --------------------------------------------------------------------

            mean_soil_vwc    = 0._r8
            mean_soil_matpot = 0._r8
            mean_soil_vwcsat = 0._r8
            areaweight       = 0._r8

            do j=1,nlevrhiz

               j_t = site_hydr%map_r2s(j,1) ! top soil layer matching rhiz layer
               j_b = site_hydr%map_r2s(j,2) ! bottom soil layer matching rhiz layer

               do j_bc = j_t,j_b

                  vwc     = bc_in(s)%h2o_liqvol_sl(j_bc)
                  psi     = site_hydr%wrf_soil(j)%p%psi_from_th(vwc) ! MLO: Any reason for not using smp_sl?
                  ! cap capillary pressure
                  ! psi = max(-1e5_r8,psi) Removing cap as that is inconstistent
                  !                        with model internals and physics. Should
                  !                        implement caps inside the functions
                  !                        if desired. (RGK 12-2021)
                  vwc_sat = bc_in(s)%watsat_sl(j_bc)
                  depth_frac = bc_in(s)%dz_sisl(j_bc)/site_hydr%dz_rhiz(j)

                  ! If there are any roots, we use root weighting
                  if(sum(site_hydr%l_aroot_layer(:),dim=1) > nearzero) then
                     layer_areaweight = site_hydr%l_aroot_layer(j)*depth_frac*pi_const*site_hydr%rs1(j)**2.0

                     ! If there are no roots, we use depth weighting
                  else
                     layer_areaweight = bc_in(s)%dz_sisl(j_bc)
                  endif

                  areaweight       = areaweight + layer_areaweight
                  mean_soil_vwc    = mean_soil_vwc + vwc*layer_areaweight
                  mean_soil_vwcsat = mean_soil_vwcsat + vwc_sat*layer_areaweight
                  mean_soil_matpot = mean_soil_matpot + psi*layer_areaweight

                  hio_soilmatpot_sl(io_si,j_bc) = psi * pa_per_mpa
                  hio_soilvwc_sl(io_si,j_bc)    = vwc
                  hio_soilvwcsat_sl(io_si,j_bc) = vwc_sat

               end do
            end do

            ! Normalization counters
            nplant_scpf(:) = 0._r8
            ncohort_scpf(:) = 0._r8
            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))
               ccohort => cpatch%shortest
               do while(associated(ccohort))
                  if ( .not. ccohort%isnew ) then
                     ! Calculate index for the scpf class
                     iscpf = ccohort%size_by_pft_class
                     nplant_scpf(iscpf) = nplant_scpf(iscpf) + ccohort%n
                     ncohort_scpf(iscpf) = ncohort_scpf(iscpf) + 1._r8
                  end if
                  ccohort => ccohort%taller
               enddo ! cohort loop
               cpatch => cpatch%younger
            end do !patch loop

            do ipft = 1, numpft
               do iscls = 1,nlevsclass
                  iscpf = (ipft-1)*nlevsclass + iscls
                  hio_sapflow_scpf(io_si,iscpf)       = site_hydr%sapflow_scpf(iscls, ipft) * ha_per_m2
                  hio_rootuptake0_scpf(io_si,iscpf)   = site_hydr%rootuptake0_scpf(iscls,ipft) * ha_per_m2
                  hio_rootuptake10_scpf(io_si,iscpf)  = site_hydr%rootuptake10_scpf(iscls,ipft) * ha_per_m2
                  hio_rootuptake50_scpf(io_si,iscpf)  = site_hydr%rootuptake50_scpf(iscls,ipft) * ha_per_m2
                  hio_rootuptake100_scpf(io_si,iscpf) = site_hydr%rootuptake100_scpf(iscls,ipft) * ha_per_m2
               end do
            end do

            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))

               ccohort => cpatch%shortest
               do while(associated(ccohort))

                  ccohort_hydr => ccohort%co_hydr

                  if ( .not. ccohort%isnew ) then

                     ! Calculate index for the scpf class
                     iscpf = ccohort%size_by_pft_class

                     ! scale up cohort fluxes to their sites
                     number_fraction_rate = (ccohort%n / nplant_scpf(iscpf)) * per_dt_tstep

                     ! scale cohorts to mean quantity
                     number_fraction = (ccohort%n / nplant_scpf(iscpf))

                     hio_errh2o_scpf(io_si,iscpf) = hio_errh2o_scpf(io_si,iscpf) + &
                          ccohort_hydr%errh2o * number_fraction_rate ! [kg/indiv/s]

                     hio_tran_scpf(io_si,iscpf) = hio_tran_scpf(io_si,iscpf) + &
                          (ccohort_hydr%qtop) * number_fraction_rate ! [kg/indiv/s]

                     hio_iterh1_scpf(io_si,iscpf)          = hio_iterh1_scpf(io_si,iscpf) + &
                          ccohort_hydr%iterh1/ncohort_scpf(iscpf)

                     hio_iterh2_scpf(io_si,iscpf)          = hio_iterh2_scpf(io_si,iscpf) + &
                          ccohort_hydr%iterh2/ncohort_scpf(iscpf)

                     mean_aroot = sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                          sum(ccohort_hydr%v_aroot_layer(:))

                     hio_ath_scpf(io_si,iscpf)             = hio_ath_scpf(io_si,iscpf) + &
                          mean_aroot * number_fraction      ! [m3 m-3]

                     hio_tth_scpf(io_si,iscpf)             = hio_tth_scpf(io_si,iscpf) + &
                          ccohort_hydr%th_troot  * number_fraction         ! [m3 m-3]

                     hio_sth_scpf(io_si,iscpf)             = hio_sth_scpf(io_si,iscpf) + &
                          ccohort_hydr%th_ag(2)  * number_fraction        ! [m3 m-3]

                     hio_lth_scpf(io_si,iscpf)             =  hio_lth_scpf(io_si,iscpf) + &
                          ccohort_hydr%th_ag(1)  * number_fraction        ! [m3 m-3]

                     mean_aroot = sum(ccohort_hydr%psi_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                          sum(ccohort_hydr%v_aroot_layer(:))

                     hio_awp_scpf(io_si,iscpf)             = hio_awp_scpf(io_si,iscpf) + &
                          mean_aroot * number_fraction * pa_per_mpa ! [Pa]

                     hio_twp_scpf(io_si,iscpf)             = hio_twp_scpf(io_si,iscpf) + &
                          ccohort_hydr%psi_troot  * number_fraction * pa_per_mpa     ! [Pa]

                     hio_swp_scpf(io_si,iscpf)             = hio_swp_scpf(io_si,iscpf) + &
                          ccohort_hydr%psi_ag(2)  * number_fraction * pa_per_mpa     ! [Pa]

                     hio_lwp_scpf(io_si,iscpf)             = hio_lwp_scpf(io_si,iscpf) + &
                          ccohort_hydr%psi_ag(1)  * number_fraction * pa_per_mpa      ! [Pa]

                     mean_aroot = sum(ccohort_hydr%ftc_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                          sum(ccohort_hydr%v_aroot_layer(:))
                     hio_aflc_scpf(io_si,iscpf)             = hio_aflc_scpf(io_si,iscpf) + &
                          mean_aroot   * number_fraction

                     hio_tflc_scpf(io_si,iscpf)             = hio_tflc_scpf(io_si,iscpf) + &
                          ccohort_hydr%ftc_troot  * number_fraction

                     hio_sflc_scpf(io_si,iscpf)             = hio_sflc_scpf(io_si,iscpf) + &
                          ccohort_hydr%ftc_ag(2)  * number_fraction

                     hio_lflc_scpf(io_si,iscpf)             = hio_lflc_scpf(io_si,iscpf) + &
                          ccohort_hydr%ftc_ag(1)  * number_fraction

                     hio_btran_scpf(io_si,iscpf)           = hio_btran_scpf(io_si,iscpf) + &
                          ccohort_hydr%btran  * number_fraction        ! [-]

                  endif

                  ccohort => ccohort%taller
               enddo ! cohort loop

               cpatch => cpatch%younger
            end do !patch loop

            if((hlm_use_ed_st3.eq.ifalse) ) then
               do iscpf=1,nlevsclass*numpft
                  if ((abs(hio_nplant_si_scpf(io_si, iscpf)-(nplant_scpf(iscpf)*ha_per_m2)) > 1.0E-8_r8) .and. &
                       (hio_nplant_si_scpf(io_si, iscpf) .ne. hlm_hio_ignore_val)) then
                     write(fates_log(),*) 'numpft:',numpft
                     write(fates_log(),*) 'nlevsclass:',nlevsclass
                     write(fates_log(),*) 'scpf:',iscpf
                     write(fates_log(),*) 'io_si:',io_si
                     write(fates_log(),*) 'hio_nplant_si_scpf:',hio_nplant_si_scpf(io_si, iscpf)
                     write(fates_log(),*) 'nplant_scpf:',nplant_scpf(iscpf)
                     write(fates_log(),*) 'nplant check on hio_nplant_si_scpf fails during hydraulics history updates'
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end if
               end do
            end if

         end do
       end associate

    end if if_hifrq1

  end subroutine update_history_hydraulics


  ! ====================================================================================
  integer function num_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(in) :: this

    num_history_vars = this%num_history_vars_

  end function num_history_vars

  ! ====================================================================================

  subroutine initialize_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    ! Determine how many of the history IO variables registered in FATES
    ! are going to be allocated
    call this%define_history_vars(initialize_variables=.false.)

    ! Allocate the list of history output variable objects
    allocate(this%hvars(this%num_history_vars()))

    ! construct the object that defines all of the IO variables
    call this%define_history_vars(initialize_variables=.true.)

  end subroutine initialize_history_vars

  ! ====================================================================================

  subroutine define_history_vars(this, initialize_variables)

    ! ---------------------------------------------------------------------------------
    !
    !                    REGISTRY OF HISTORY OUTPUT VARIABLES
    !
    ! This subroutine is called in two contexts, either in count mode or initialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
    !
    ! Note 1 there are different ways you can flush or initialize the output fields.
    ! If you flush to a native type, (such as zero), the entire slab which covers
    ! indices which may not be relevant to FATES, are flushed to this value.  So
    ! in that case, lakes and crops that are not controlled by FATES will zero'd
    ! and when values are scaled up to the land-grid, the zero's for non FATES will
    ! be included.  This is good and correct if nothing is there.
    !
    ! But, what if crops exist in the host model and occupy a fraction of the land-surface
    ! shared with natural vegetation? In that case, you want to flush your arrays
    ! with a value that the HLM treats as "do not average"
    !
    ! If your HLM makes use of, and you want, INTEGER OUTPUT, pass the flushval as
    ! a real.  The applied flush value will use the NINT() intrinsic function
    ! ---------------------------------------------------------------------------------

    use FatesIOVariableKindMod, only : site_r8, site_soil_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_pft_r8, site_coage_r8
    use FatesIOVariableKindMod, only : site_height_r8, site_agefuel_r8
    use FatesInterfaceTypesMod, only : hlm_use_planthydro

    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_cdsc_r8, site_cdpf_r8, site_cdam_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8, site_clscpf_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8


    implicit none

    class(fates_history_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?

    integer :: ivar
    character(len=10) :: tempstring

    ivar=0

    ! Variable names should start with the 'FATES_' prefix and end with a suffix
    ! depending on how it is indexed (i.e. the dimension):
    ! site                     (site_r8)        : no suffix
    ! cohort age               (site_coage_r8)  : AC
    ! patch age                (site_age_r8)    : AP
    ! canopy layer             (site_can_r8)    : CL
    ! coarse woody debris size (site_cwdsc_r8)  : DC
    ! element                  (site_elem_r8)   : EL
    ! leaf layer                                : LL
    ! fuel class               (site_fuel_r8)   : FC
    ! height                   (site_height_r8) : HT
    ! plant functional type    (site_pft_r8)    : PF
    ! soil layer               (site_soil_r8)   : SL
    ! cohort size              (site_size_r8)   : SZ
    ! cohort crown damage      (site_cd_r8)     : CD

    ! Multiple dimensions should have multiple two-code suffixes:
    ! cohort age x pft                (site_cooage_r8)   : ACPF
    ! patch age x fuel class          (site_agefuel_r8)  : APFC
    ! patch age x pft                 (site_agepft_r8)   : APPF
    ! canopy layer x leaf layer       (site_cnlf_r8)     : CLLL
    ! canopy layer x leaf layer x pft (site_cnlfpft_r8)  : CLLLPF
    ! element x cwd size              (site_elcwd_r8)    : ELDC
    ! cohort size x patch age         (site_scag_r8)     : SZAP
    ! cohort size x patch age x pft   (site_scagpft_r8)  : SZAPPF
    ! cohort size x pft               (site_size_pft_r8) : SZPF
    ! canopy layer x size x pft       (site_clscpf_r8)   : CLSZPF (NOT ACTIVE)
    ! cohort size x crown damage       (site_cdsc_r8)     : SZCD
    ! cohort size x crown damage x pft (site_cdpf_r8)     : CDPF


    if_dyn0: if(hlm_hist_level_dynam>0) then

       ! Site level counting variables
       call this%set_history_var(vname='FATES_NPATCHES', units='',                &
            long='total number of patches per site', use_default='active',        &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_npatches_si)

       call this%set_history_var(vname='FATES_NCOHORTS', units='',                &
            long='total number of cohorts per site', use_default='active',        &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_ncohorts_si)

       ! Patch variables
       call this%set_history_var(vname='FATES_TRIMMING', units='1',               &
            long='degree to which canopy expansion is limited by leaf economics (0-1)', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_trimming_si)

       call this%set_history_var(vname='FATES_AREA_PLANTS', units='m2 m-2',       &
            long='area occupied by all plants per m2 land area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,      &
            initialize=initialize_variables, index=ih_area_plant_si)

       call this%set_history_var(vname='FATES_AREA_TREES', units='m2 m-2',        &
            long='area occupied by woody plants per m2 land area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_area_trees_si)

       call this%set_history_var(vname='FATES_FRACTION', units='m2 m-2',          &
            long='total gridcell fraction which FATES is running over', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fates_fraction_si, flush_to_zero=.true.)

       call this%set_history_var(vname='FATES_BA_WEIGHTED_HEIGHT', units='m',        &
            long='basal area-weighted mean height of woody plants', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_ba_weighted_height_si)

       call this%set_history_var(vname='FATES_CA_WEIGHTED_HEIGHT', units='m',        &
            long='crown area-weighted mean height of canopy plants', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_ca_weighted_height_si)

       call this%set_history_var(vname='FATES_COLD_STATUS', units='',             &
            long='site-level cold status, 0=not cold-dec, 1=too cold for leaves, 2=not too cold',  &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index=ih_site_cstatus_si)

       call this%set_history_var(vname='FATES_GDD', units='degree_Celsius',       &
            long='site-level growing degree days', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables, index=ih_gdd_si)

       call this%set_history_var(vname='FATES_NCHILLDAYS', units = 'days',        &
            long='site-level number of chill days', use_default='active',         &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_site_nchilldays_si)

       call this%set_history_var(vname='FATES_NCOLDDAYS', units = 'days',         &
            long='site-level number of cold days', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_site_ncolddays_si)

       call this%set_history_var(vname='FATES_DAYSINCE_COLDLEAFOFF',              &
            units='days', long='site-level days elapsed since cold leaf drop',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_cleafoff_si)

       call this%set_history_var(vname='FATES_DAYSINCE_COLDLEAFON',               &
            units='days', long='site-level days elapsed since cold leaf flush',   &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_cleafon_si)

       call this%set_history_var(vname='FATES_CANOPY_SPREAD', units='',           &
            long='scaling factor (0-1) between tree basal area and canopy area',  &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_canopy_spread_si)

       call this%set_history_var(vname='FATES_LAI', units='m2 m-2',               &
            long='total leaf area index per m2 land area',                        &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_lai_si)
       
       call this%set_history_var(vname='FATES_ELAI', units='m2 m-2',               &
            long='exposed (non snow-occluded) leaf area index per m2 land area',                        &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_elai_si)
       
       call this%set_history_var(vname='FATES_WOOD_PRODUCT', units='kg m-2',      &
            long='total wood product from logging in kg carbon per m2 land area', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',   &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_woodproduct_si)

       ! Fire Variables

       call this%set_history_var(vname='FATES_NESTEROV_INDEX', units='',          &
            long='nesterov fire danger index', use_default='active',              &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_nesterov_fire_danger_si)

       call this%set_history_var(vname='FATES_IGNITIONS',                         &
            units='m-2 s-1',                                                      &
            long='number of successful fire ignitions per m2 land area per second',  &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_nignitions_si)

       call this%set_history_var(vname='FATES_FDI', units='1',                    &
            long='Fire Danger Index (probability that an ignition will lead to a fire)', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_fdi_si)

       call this%set_history_var(vname='FATES_ROS', units='m s-1',                &
            long='fire rate of spread in meters per second',                      &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_spitfire_ros_si)

       call this%set_history_var(vname='FATES_EFFECT_WSPEED', units='m s-1',      &
            long ='effective wind speed for fire spread in meters per second',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_effect_wspeed_si)

       call this%set_history_var(vname='FATES_FUELCONSUMED', units='kg m-2',      &
            long ='total fuel consumed in kg carbon per m2 land area',            &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_tfc_ros_si)

       call this%set_history_var(vname='FATES_FIRE_INTENSITY',                    &
            units='J m-1 s-1',                                                    &
            long='spitfire surface fireline intensity in J per m per second',     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_intensity_si)

       call this%set_history_var(vname='FATES_FIRE_INTENSITY_BURNFRAC',           &
            units='J m-1 s-1',                                                    &
            long='product of surface fire intensity and burned area fraction -- divide by FATES_BURNFRAC to get area-weighted mean intensity', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_intensity_area_product_si)

       call this%set_history_var(vname='FATES_BURNFRAC', units='s-1',             &
            long='burned area fraction per second', use_default='active',         &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',                           &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_area_si)

       call this%set_history_var(vname='FATES_FUEL_MEF', units='m3 m-3',          &
            long='fuel moisture of extinction (volumetric)',                      &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index=ih_fire_fuel_mef_si)

       call this%set_history_var(vname='FATES_FUEL_BULKD',                        &
            units='kg m-3', long='fuel bulk density in kg per m3',                &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_fire_fuel_bulkd_si )

       call this%set_history_var(vname='FATES_FUEL_EFF_MOIST', units='m3 m-3',    &
            long='spitfire fuel moisture (volumetric)', use_default='active',     &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,      &
            initialize=initialize_variables, index = ih_fire_fuel_eff_moist_si)

       call this%set_history_var(vname='FATES_FUEL_SAV', units='m-1',             &
            long='spitfire fuel surface area to volume ratio',                    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_fire_fuel_sav_si)

       call this%set_history_var(vname='FATES_FUEL_AMOUNT', units='kg m-2',       &
            long='total ground fuel related to FATES_ROS (omits 1000hr fuels) in kg C per m2 land area',   &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_sum_fuel_si)
       ! Litter Variables

       call this%set_history_var(vname='FATES_LITTER_IN', units='kg m-2 s-1',     &
            long='litter flux in kg carbon per m2 per second',                    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_litter_in_si)

       call this%set_history_var(vname='FATES_LITTER_OUT', units='kg m-2 s-1',    &
            long='litter flux out in kg carbon (exudation, fragmentation, seed decay)',   &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_litter_out_si)

       call this%set_history_var(vname='FATES_SEED_BANK', units='kg m-2',         &
            long='total seed mass of all PFTs in kg carbon per m2 land area',     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_seed_bank_si)

       call this%set_history_var(vname='FATES_UNGERM_SEED_BANK', units='kg m-2',         &
            long='ungerminated seed mass of all PFTs in kg carbon per m2 land area',     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_ungerm_seed_bank_si)

       call this%set_history_var(vname='FATES_SEEDLING_POOL', units='kg m-2',         &
            long='total seedling (ie germinated seeds) mass of all PFTs in kg carbon per m2 land area',     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_seedling_pool_si)

       call this%set_history_var(vname='FATES_SEEDS_IN', units='kg m-2 s-1',      &
            long='seed production rate in kg carbon per m2 second',               &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_seeds_in_si)

       call this%set_history_var(vname='FATES_SEEDS_IN_LOCAL', units='kg m-2 s-1',      &
            long='local seed production rate in kg carbon per m2 second',               &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_seeds_in_local_si)

       call this%set_history_var(vname='FATES_STOREC', units='kg m-2',            &
            long='total biomass in live plant storage in kg carbon per m2 land area', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_storec_si)

       call this%set_history_var(vname='FATES_STOREC_TF', units='kg kg-1',         &
            long='Storage C fraction of target', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,   &
            ivar=ivar, initialize=initialize_variables, index = ih_storectfrac_si )

       call this%set_history_var(vname='FATES_VEGC', units='kg m-2',              &
            long='total biomass in live plants in kg carbon per m2 land area',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_totvegc_si)

       call this%set_history_var(vname='FATES_SAPWOODC', units='kg m-2',          &
            long='total biomass in live plant sapwood in kg carbon per m2',       &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_sapwc_si)

       call this%set_history_var(vname='FATES_LEAFC', units='kg m-2',             &
            long='total biomass in live plant leaves in kg carbon per m2',        &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_leafc_si)

       call this%set_history_var(vname='FATES_FROOTC', units='kg m-2',            &
            long='total biomass in live plant fine roots in kg carbon per m2',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_fnrtc_si)

       call this%set_history_var(vname='FATES_REPROC', units='kg m-2',            &
            long='total biomass in live plant reproductive tissues in kg carbon per m2', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_reproc_si)

       call this%set_history_var(vname='FATES_NPP', units='kg m-2 s-1',           &
            long='net primary production in kg carbon per m2 per second',         &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables, index = ih_npp_si)

       call this%set_history_var(vname='FATES_AUTORESP', units='kg m-2 s-1',     &
            long='autotrophic respiration in kg carbon per m2 per second',        &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables, index = ih_aresp_si)

       call this%set_history_var(vname='FATES_GROWTH_RESP', units='kg m-2 s-1',   &
            long='growth respiration in kg carbon per m2 per second',             &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_growth_resp_si)


       ! Output specific to the chemical species dynamics used (parteh)
       call this%set_history_var(vname='FATES_L2FR', units='kg kg-1',                   &
            long='The leaf to fineroot biomass multiplier for target allometry', & 
            use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,    &
            ivar=ivar, initialize=initialize_variables, index = ih_l2fr_si)

       nitrogen_active_if0: if(any(element_list(:)==nitrogen_element)) then

          call this%set_history_var(vname='FATES_NH4UPTAKE', units='kg m-2 s-1',  &
               long='ammonium uptake rate by plants in kg NH4 per m2 per second', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_nh4uptake_si)

          call this%set_history_var(vname='FATES_NO3UPTAKE', units='kg m-2 s-1',  &
               long='nitrate uptake rate by plants in kg NO3 per m2 per second',  &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_no3uptake_si)

          call this%set_history_var(vname='FATES_NEFFLUX', units='kg m-2 s-1',    &
               long='nitrogen effluxed from plant in kg N per m2 per second (unused)', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_nefflux_si)

          call this%set_history_var(vname='FATES_NDEMAND', units='kg m-2 s-1',      &
               long='plant nitrogen need (algorithm dependent) in kg N per m2 per second', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_ndemand_si)
          
          call this%set_history_var(vname='FATES_NFIX_SYM', units='kg m-2 s-1',      &
               long='symbiotic dinitrogen fixation in kg N per m2 per second', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_nfix_si)
          
          call this%set_history_var(vname='FATES_STOREN', units='kg m-2',         &
               long='total nitrogen in live plant storage', use_default='active', &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,              &
               ivar=ivar, initialize=initialize_variables, index = ih_storen_si)

          call this%set_history_var(vname='FATES_STOREN_TF', units='1',           &
               long='storage N fraction of target', use_default='active',         &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,   &
               initialize=initialize_variables, index = ih_storentfrac_si)

          call this%set_history_var(vname='FATES_VEGN', units='kg m-2',           &
               long='total nitrogen in live plants', use_default='active',        &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,   &
               initialize=initialize_variables, index = ih_totvegn_si)

          call this%set_history_var(vname='FATES_SAPWOODN', units='kg m-2',       &
               long='total nitrogen in live plant sapwood', use_default='active', &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,              &
               ivar=ivar, initialize=initialize_variables, index = ih_sapwn_si)

          call this%set_history_var(vname='FATES_LEAFN', units='kg m-2',          &
               long='total nitrogen in live plant leaves', use_default='active',  &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,   &
               initialize=initialize_variables, index = ih_leafn_si)

          call this%set_history_var(vname='FATES_FROOTN', units='kg m-2',         &
               long='total nitrogen in live plant fine-roots',                    &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_fnrtn_si)

          call this%set_history_var(vname='FATES_REPRON', units='kg m-2',         &
               long='total nitrogen in live plant reproductive tissues',          &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_repron_si)
       end if nitrogen_active_if0

       phosphorus_active_if0: if(any(element_list(:)==phosphorus_element)) then
          call this%set_history_var(vname='FATES_STOREP', units='kg m-2',         &
               long='total phosphorus in live plant storage',                     &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_storep_si)

          call this%set_history_var(vname='FATES_STOREP_TF', units='1',           &
               long='storage P fraction of target', use_default='active',         &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,              &
               ivar=ivar, initialize=initialize_variables,                        &
               index = ih_storeptfrac_si)

          call this%set_history_var(vname='FATES_VEGP', units='kg m-2',           &
               long='total phosphorus in live plants', use_default='active',      &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,              &
               ivar=ivar, initialize=initialize_variables, index = ih_totvegp_si)

          call this%set_history_var(vname='FATES_SAPWOODP', units='kg m-2',       &
               long='Total phosphorus in live plant sapwood', use_default='active', &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,   &
               initialize=initialize_variables, index = ih_sapwp_si)

          call this%set_history_var(vname='FATES_LEAFP', units='kg m-2',          &
               long='total phosphorus in live plant leaves',                      &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_leafp_si)

          call this%set_history_var(vname='FATES_FROOTP', units='kg m-2',         &
               long='total phosphorus in live plant fine roots',                  &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_fnrtp_si)

          call this%set_history_var(vname='FATES_REPROP', units='kg m-2',         &
               long='total phosphorus in live plant reproductive tissues',        &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_reprop_si)

          call this%set_history_var(vname='FATES_PUPTAKE', units='kg m-2 s-1',    &
               long='mineralized phosphorus uptake rate of plants in kg P per m2 per second', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_puptake_si)

          call this%set_history_var(vname='FATES_PEFFLUX', units='kg m-2 s-1',    &
               long='phosphorus effluxed from plant in kg P per m2 per second (unused)', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_pefflux_si)

          call this%set_history_var(vname='FATES_PDEMAND', units='kg m-2 s-1',      &
               long='plant phosphorus need (algorithm dependent) in kg P per m2 per second', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_pdemand_si)
          
       end if phosphorus_active_if0

       call this%set_history_var(vname='FATES_STRUCTC', units='kg m-2',           &
            long='structural biomass in kg carbon per m2 land area',              &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_bdead_si)

       call this%set_history_var(vname='FATES_NONSTRUCTC', units='kg m-2',        &
            long='non-structural biomass (sapwood + leaf + fineroot) in kg carbon per m2', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_balive_si)

       call this%set_history_var(vname='FATES_VEGC_ABOVEGROUND', units='kg m-2',  &
            long='aboveground biomass in kg carbon per m2 land area',             &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_agb_si)

       call this%set_history_var(vname='FATES_CANOPY_VEGC', units='kg m-2',       &
            long='biomass of canopy plants in kg carbon per m2 land area',        &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_canopy_biomass_si)

       call this%set_history_var(vname='FATES_USTORY_VEGC', units='kg m-2',   &
            long='biomass of understory plants in kg carbon per m2 land area',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_understory_biomass_si)

       ! disturbance rates

       call this%set_history_var(vname='FATES_PRIMARY_PATCHFUSION_ERR',           &
            units='m2 m-2 yr-1',                                                  &
            long='error in total primary lands associated with patch fusion',     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_primaryland_fusion_error_si)

       call this%set_history_var(vname='FATES_DISTURBANCE_RATE_FIRE',             &
            units='m2 m-2 yr-1', long='disturbance rate from fire',               &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_fire_disturbance_rate_si)

       call this%set_history_var(vname='FATES_DISTURBANCE_RATE_LOGGING',          &
            units='m2 m-2 yr-1', long='disturbance rate from logging',            &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_logging_disturbance_rate_si)

       call this%set_history_var(vname='FATES_DISTURBANCE_RATE_TREEFALL',         &
            units='m2 m-2 yr-1', long='disturbance rate from treefall',           &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_fall_disturbance_rate_si)
       
       call this%set_history_var(vname='FATES_HARVEST_WOODPROD_C_FLUX',           &
            units='kg m-2 yr-1',                                                  &
            long='harvest-associated wood product carbon flux in kg C per m2 per year', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_harvest_woodprod_carbonflux_si)
       
       call this%set_history_var(vname='FATES_LUCHANGE_WOODPROD_C_FLUX',     &
            units='kg m-2 yr-1',                                                  &
            long='land-use-change-associated wood product carbon flux in kg C per m2 per year', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_luchange_woodprod_carbonflux_si)
       
       call this%set_history_var(vname='FATES_TVEG24', units='degree_Celsius', &
            long='fates 24-hr running mean vegetation temperature by site', &
            use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, &
            ivar=ivar, initialize=initialize_variables, index = ih_tveg24_si )

       call this%set_history_var(vname='FATES_TLONGTERM', units='degree_Celsius', &
            long='fates 30-year running mean vegetation temperature by site', &
            use_default='inactive', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, &
            ivar=ivar, initialize=initialize_variables, index = ih_tlongterm_si )

       call this%set_history_var(vname='FATES_TGROWTH', units='degree_Celsius', &
            long='fates long-term running mean vegetation temperature by site', &
            use_default='inactive', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple, &
            ivar=ivar, initialize=initialize_variables, index = ih_tgrowth_si )

       call this%set_history_var(vname='FATES_HARVEST_DEBT', units='kg C',                   &
            long='Accumulated carbon failed to be harvested',  use_default='active',     &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,   &
            ivar=ivar, initialize=initialize_variables, index = ih_harvest_debt_si )

       call this%set_history_var(vname='FATES_HARVEST_DEBT_SEC', units='kg C',                   &
            long='Accumulated carbon failed to be harvested from secondary patches',  use_default='active',     &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_dyna_simple,   &
            ivar=ivar, initialize=initialize_variables, index = ih_harvest_debt_sec_si )

       ! Nutrient flux variables (dynamics call frequency)
       ! ----------------------------------------------------
       call this%set_history_var(vname='FATES_EXCESS_RESP', units='kg m-2 s-1',    &
            long='respiration of un-allocatable carbon gain', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_nflx_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_excess_resp_si)

       
       ! slow carbon fluxes associated with mortality from or transfer betweeen canopy and understory

       call this%set_history_var(vname='FATES_DEMOTION_CARBONFLUX',               &
            units = 'kg m-2 s-1',                                                &
            long='demotion-associated biomass carbon flux from canopy to understory in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_demotion_carbonflux_si)

       call this%set_history_var(vname='FATES_PROMOTION_CARBONFLUX',              &
            units = 'kg m-2 s-1',                                                &
            long='promotion-associated biomass carbon flux from understory to canopy in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_promotion_carbonflux_si)

       call this%set_history_var(vname='FATES_MORTALITY_CFLUX_CANOPY',            &
            units = 'kg m-2 s-1',                                                &
            long='flux of biomass carbon from live to dead pools from mortality of canopy plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_canopy_mortality_carbonflux_si)

       call this%set_history_var(vname='FATES_MORTALITY_CFLUX_USTORY',            &
            units = 'kg m-2 s-1',                                                &
            long='flux of biomass carbon from live to dead pools from mortality of understory plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_understory_mortality_carbonflux_si)

       call this%set_history_var(vname='MORTALITY_CROWNAREA_CANOPY',              &
            units = 'm2/ha/year',                                                &
            long='Crown area of canopy trees that died',                         &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_canopy_mortality_crownarea_si )

       call this%set_history_var(vname='MORTALITY_CROWNAREA_UNDERSTORY',          &
            units = 'm2/ha/year',                                                &
            long='Crown aera of understory trees that died',                     &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_understory_mortality_crownarea_si )

       call this%set_history_var(vname='FATES_FIRE_CLOSS', units='kg m-2 s-1',    &
            long='carbon loss to atmosphere from fire in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_fire_c_to_atm_si)

       call this%set_history_var(vname='FATES_CBALANCE_ERROR',                    &
            units='kg s-1',                                                       &
            long='total carbon error in kg carbon per second',                    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_cbal_err_fates_si)


       
       call this%set_history_var(vname='FATES_LEAF_ALLOC', units='kg m-2 s-1',    &
            long='allocation to leaves in kg carbon per m2 per second',          &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_leaf_si)

       call this%set_history_var(vname='FATES_SEED_ALLOC', units='kg m-2 s-1',    &
            long='allocation to seeds in kg carbon per m2 per second',           &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_seed_si)

       call this%set_history_var(vname='FATES_STEM_ALLOC', units='kg m-2 s-1',    &
            long='allocation to stem in kg carbon per m2 per second',            &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_stem_si)

       call this%set_history_var(vname='FATES_FROOT_ALLOC', units='kg m-2 s-1',   &
            long='allocation to fine roots in kg carbon per m2 per second',      &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_froot_si)

       call this%set_history_var(vname='FATES_CROOT_ALLOC', units='kg m-2 s-1',   &
            long='allocation to coarse roots in kg carbon per m2 per second',    &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_croot_si)

       call this%set_history_var(vname='FATES_STORE_ALLOC', units='kg m-2 s-1',   &
            long='allocation to storage tissues in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_npp_stor_si)

       call this%set_history_var(vname='FATES_GRAZING', units='kg m-2 s-1',    &
            long='grazing by herbivores of leaves in kg carbon per m2 per second',          &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_grazing_si)

       hydro_active_if: if(hlm_use_planthydro.eq.itrue) then
          call this%set_history_var(vname='FATES_VEGH2O_DEAD', units = 'kg m-2',  &
               long='cumulative water stored in dead biomass due to mortality',  &
               use_default='inactive', avgflag='A', vtype=site_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,                              &
               initialize=initialize_variables, index = ih_h2oveg_dead_si)

          call this%set_history_var(vname='FATES_VEGH2O_RECRUIT',                 &
               units = 'kg m-2', long='amount of water in new recruits',         &
               use_default='inactive', avgflag='A', vtype=site_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,                              &
               initialize=initialize_variables, index = ih_h2oveg_recruit_si)

          call this%set_history_var(vname='FATES_VEGH2O_GROWTURN_ERR',            &
               units = 'kg m-2',                                                 &
               long='cumulative net borrowed (+) or lost (-) from water storage due to combined growth & turnover', &
               use_default='inactive', avgflag='A', vtype=site_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_simple, ivar=ivar,                              &
               initialize=initialize_variables, index = ih_h2oveg_growturn_err_si)
       end if hydro_active_if

       if_crowndamage1: if(hlm_use_tree_damage .eq. itrue) then

          call this%set_history_var(vname='FATES_CROWNAREA_CANOPY_CD', units = 'm2 m-2 yr-1',         &
               long='crownarea lost to damage each year', use_default='inactive',   &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables, index = ih_crownarea_canopy_damage_si )
          
          call this%set_history_var(vname='FATES_CROWNAREA_USTORY_CD', units = 'm2 m-2 yr-1',         &
               long='crownarea lost to damage each year', use_default='inactive',   &
               avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
               upfreq=group_dyna_simple, ivar=ivar, initialize=initialize_variables, index = ih_crownarea_ustory_damage_si )

       end if if_crowndamage1
       

       if_dyn1: if(hlm_hist_level_dynam>1) then

          call this%set_history_var(vname='FATES_NPP_LU', units='kg m-2 s-1',        &
               long='net primary productivity by land use type in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_landuse_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_si_landuse)

          call this%set_history_var(vname='FATES_PATCHAREA_LU', units='m2 m-2',      &
               long='patch area by land use type', use_default='active',  &
               avgflag='A', vtype=site_landuse_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, &
               ivar=ivar, initialize=initialize_variables, index=ih_area_si_landuse)

          call this%set_history_var(vname='FATES_VEGC_LU', units='kg m-2',      &
               long='Vegetation Carbon by land use type', use_default='active',  &
               avgflag='A', vtype=site_landuse_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, &
               ivar=ivar, initialize=initialize_variables, index=ih_biomass_si_landuse)

          call this%set_history_var(vname='FATES_BURNEDAREA_LU', units='s-1',      &
               long='burned area by land use type', use_default='active',  &
               avgflag='A', vtype=site_landuse_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, &
               ivar=ivar, initialize=initialize_variables, index=ih_burnedarea_si_landuse)

          call this%set_history_var(vname='FATES_TRANSITION_MATRIX_LULU', units='m2 m-2 yr-1',      &
               long='land use transition matrix', use_default='active',  &
               avgflag='A', vtype=site_lulu_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,  &
               initialize=initialize_variables, index=ih_transition_matrix_si_lulu)
       
          call this%set_history_var(vname='FATES_DISTURBANCE_RATE_MATRIX_LULU', units='m2 m-2 yr-1',   &
               long='disturbance rates by land use type x land use type matrix', use_default='active', &
               avgflag='A', vtype=site_lulu_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,  &
               initialize=initialize_variables, index=ih_disturbance_rate_si_lulu)
          
          call this%set_history_var(vname='FATES_VEGC_PF', units='kg m-2',           &
               long='total PFT-level biomass in kg of carbon per land area',         &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_biomass_si_pft)

          call this%set_history_var(vname='FATES_RECRUITMENT_CFLUX_PF', units='kg m-2 yr-1',  &
               long='total PFT-level biomass of new recruits in kg of carbon per land area',         &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=1, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_recruitment_cflux_si_pft)

          call this%set_history_var(vname='FATES_LEAFC_PF', units='kg m-2',          &
               long='total PFT-level leaf biomass in kg carbon per m2 land area',    &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_leafbiomass_si_pft)

          call this%set_history_var(vname='FATES_STOREC_PF', units='kg m-2',         &
               long='total PFT-level stored biomass in kg carbon per m2 land area',  &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_storebiomass_si_pft)

          call this%set_history_var(vname='FATES_CROWNAREA_PF',  units='m2 m-2',     &
               long='total PFT-level crown area per m2 land area',                   &
               use_default='active', avgflag='A', vtype=site_pft_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_crownarea_si_pft)

          call this%set_history_var(vname='FATES_CANOPYCROWNAREA_PF',                &
               units='m2 m-2', long='total PFT-level canopy-layer crown area per m2 land area', &
               use_default='active', avgflag='A', vtype=site_pft_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_canopycrownarea_si_pft)

          call this%set_history_var(vname='FATES_GPP_PF', units='kg m-2 s-1',        &
               long='total PFT-level GPP in kg carbon per m2 land area per second',  &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_gpp_si_pft)

          call this%set_history_var(vname='FATES_NPP_PF', units='kg m-2 s-1',       &
               long='total PFT-level NPP in kg carbon per m2 land area per second',  &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_npp_si_pft)

          call this%set_history_var(vname='FATES_NPLANT_PF', units='m-2',           &
               long='total PFT-level number of individuals per m2 land area',        &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_nindivs_si_pft)

          call this%set_history_var(vname='FATES_RECRUITMENT_PF',                    &
               units='m-2 yr-1',                                                     &
               long='PFT-level recruitment rate in number of individuals per m2 land area per year',  &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_recruitment_si_pft)

          call this%set_history_var(vname='FATES_SEEDS_IN_GRIDCELL_PF',                    &
               units='kg',                                                      &
               long='Site-level seed mass input from neighboring gridcells per pft',  &
               use_default='inactive', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_seeds_in_gc_si_pft)

          call this%set_history_var(vname='FATES_SEEDS_OUT_GRIDCELL_PF',                    &
               units='kg',                                                      &
               long='Site-level seed mass output to neighboring gridcells per pft',  &
               use_default='inactive', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_seeds_out_gc_si_pft)

          call this%set_history_var(vname='FATES_MORTALITY_PF', units='m-2 yr-1',    &
               long='PFT-level mortality rate in number of individuals per m2 land area per year', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_mortality_si_pft)

          !MLO - Drought-deciduous phenology variables are now defined for each PFT.
          call this%set_history_var(vname='FATES_DROUGHT_STATUS_PF',                     &
               units='',                                                                &
               long='PFT-level drought status, <2 too dry for leaves, >=2 not too dry', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_site_dstatus_si_pft)

          call this%set_history_var(vname='FATES_DAYSINCE_DROUGHTLEAFOFF_PF',           &
               units='days', long='PFT-level days elapsed since drought leaf drop',     &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_dleafoff_si_pft)

          call this%set_history_var(vname='FATES_DAYSINCE_DROUGHTLEAFON_PF',            &
               units='days',                                                            &
               long='PFT-level days elapsed since drought leaf flush',                  &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_dleafon_si_pft)

          call this%set_history_var(vname='FATES_MEANLIQVOL_DROUGHTPHEN_PF',            &
               units='m3 m-3',                                                          &
               long='PFT-level mean liquid water volume for drought phenolgy',          &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_meanliqvol_si_pft)

          call this%set_history_var(vname='FATES_MEANSMP_DROUGHTPHEN_PF',               &
               units='Pa',                                                              &
               long='PFT-level mean soil matric potential for drought phenology',       &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_meansmp_si_pft)

          call this%set_history_var(vname='FATES_ELONG_FACTOR_PF',                      &
               units='1',                                                               &
               long='PFT-level mean elongation factor (partial flushing/abscission)',   &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM',    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                    &
               index=ih_elong_factor_si_pft)

          nocomp_if: if (hlm_use_nocomp .eq. itrue) then
             call this%set_history_var(vname='FATES_NOCOMP_NPATCHES_PF', units='',      &
                  long='number of patches per PFT (nocomp-mode-only)',                  &
                  use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
                  index=ih_nocomp_pftnpatches_si_pft)

             call this%set_history_var(vname='FATES_NOCOMP_PATCHAREA_PF', units='m2 m-2',&
                  long='total patch area allowed per PFT (nocomp-mode-only)',           &
                  use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
                  index=ih_nocomp_pftpatchfraction_si_pft)

             call this%set_history_var(vname='FATES_NOCOMP_BURNEDAREA_PF', units='s-1', &
                  long='total burned area of PFT-labeled patch area (nocomp-mode-only)',&
                  use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
                  index=ih_nocomp_pftburnedarea_si_pft)
          endif nocomp_if

          ! patch age class variables
          call this%set_history_var(vname='FATES_PATCHAREA_AP', units='m2 m-2',      &
               long='patch area by age bin per m2 land area', use_default='active',  &
               avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,  &
               initialize=initialize_variables, index=ih_area_si_age)

          call this%set_history_var(vname='FATES_LAI_AP', units='m2 m-2',            &
               long='total leaf area index by age bin per m2 land area',                   &
               use_default='active', avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_lai_si_age)



          call this%set_history_var(vname='FATES_CANOPYAREA_AP', units='m2 m-2',     &
               long='canopy area by age bin per m2 land area', use_default='active', &
               avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,  &
               initialize=initialize_variables, index=ih_canopy_area_si_age)

          call this%set_history_var(vname='FATES_NCL_AP', units='',                  &
               long='number of canopy levels by age bin',                            &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_ncl_si_age)

          call this%set_history_var(vname='FATES_NPATCH_AP', units='',               &
               long='number of patches by age bin', use_default='inactive',          &
               avgflag='A', vtype=site_age_r8, hlms='CLM:ALM',                       &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_npatches_si_age)

          if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
             tempstring = 'active'
          else
             tempstring = 'inactive'
          endif

          call this%set_history_var(vname='FATES_ZSTAR_AP', units='m',               &
               long='product of zstar and patch area by age bin (divide by FATES_PATCHAREA_AP to get mean zstar)', &
               use_default=trim(tempstring), avgflag='A', vtype=site_age_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_zstar_si_age)

          call this%set_history_var(vname='FATES_CANOPYAREA_HT', units='m2 m-2',     &
               long='canopy area height distribution',                               &
               use_default='active', avgflag='A', vtype=site_height_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_canopy_height_dist_si_height)

          call this%set_history_var(vname='FATES_LEAFAREA_HT', units='m2 m-2',       &
               long='leaf area height distribution', use_default='active',           &
               avgflag='A', vtype=site_height_r8, hlms='CLM:ALM',                    &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_leaf_height_dist_si_height)

          call this%set_history_var(vname='FATES_VEGC_AP', units='kg m-2',           &
               long='total biomass within a given patch age bin in kg carbon per m2 land area', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_biomass_si_age)

          call this%set_history_var(vname='FATES_SECONDARY_ANTHRODISTAGE_AP',          &
               units='m2 m-2',                                                       &
               long='secondary forest patch area age distribution since anthropogenic disturbance', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_agesince_anthrodist_si_age)

          call this%set_history_var(vname='FATES_SECONDARY_AREA_AP',                &
               units='m2 m-2',                                                       &
               long='secondary forest patch area age distribution since any kind of disturbance', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_secondarylands_area_si_age)

          call this%set_history_var(vname='FATES_PRIMARY_AREA_AP',                &
               units='m2 m-2',                                                       &
               long='primary forest patch area age distribution since any kind of disturbance', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index=ih_primarylands_area_si_age)

          call this%set_history_var(vname='FATES_FRAGMENTATION_SCALER_SL', units='', &
               long='factor (0-1) by which litter/cwd fragmentation proceeds relative to max rate by soil layer',  &
               use_default='active', avgflag='A', vtype=site_soil_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_fragmentation_scaler_sl)

          call this%set_history_var(vname='FATES_FUEL_MOISTURE_FC', units='m3 m-3',  &
               long='spitfire fuel class-level fuel moisture (volumetric)',          &
               use_default='active', avgflag='A', vtype=site_fuel_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_litter_moisture_si_fuel)

          call this%set_history_var(vname='FATES_FUEL_AMOUNT_FC', units='kg m-2',    &
               long='spitfire fuel-class level fuel amount in kg carbon per m2 land area', &
               use_default='active', avgflag='A', vtype=site_fuel_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_fuel_amount_si_fuel)

          call this%set_history_var(vname='FATES_FUEL_AMOUNT_APFC', units='kg m-2',  &
               long='spitfire fuel quantity in each age x fuel class in kg carbon per m2 land area', &
               use_default='inactive', avgflag='A', vtype=site_agefuel_r8,           &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_fuel_amount_age_fuel)

          call this%set_history_var(vname='FATES_BURNFRAC_AP', units='s-1',          &
               long='spitfire fraction area burnt (per second) by patch age',        &
               use_default='active', avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index = ih_area_burnt_si_age)

          call this%set_history_var(vname='FATES_FIRE_INTENSITY_BURNFRAC_AP',        &
               units='J m-1 s-1', &
               long='product of fire intensity and burned fraction, resolved by patch age (so divide by FATES_BURNFRAC_AP to get burned-area-weighted-average intensity)', &
               use_default='active', avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index = ih_fire_intensity_si_age)

          call this%set_history_var(vname='FATES_FUEL_AMOUNT_AP', units='kg m-2',    &
               long='spitfire ground fuel (kg carbon per m2) related to FATES_ROS (omits 1000hr fuels) within each patch age bin (divide by FATES_PATCHAREA_AP to get fuel per unit area of that-age patch)', &
               use_default='active', avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index = ih_fire_sum_fuel_si_age)

          call this%set_history_var(vname='FATES_FUEL_BURNT_BURNFRAC_FC', units='1', &
               long='product of fraction (0-1) of fuel burnt and burnt fraction (divide by FATES_BURNFRAC to get burned-area-weighted mean fraction fuel burnt)', &
               use_default='active', avgflag='A', vtype=site_fuel_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_burnt_frac_litter_si_fuel)

          call this%set_history_var(vname='FATES_LITTER_IN_EL', units='kg m-2 s-1',  &
               long='litter flux in in kg element per m2 per second',                &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_litter_in_elem)

          call this%set_history_var(vname='FATES_LITTER_OUT_EL', units='kg m-2 s-1', &
               long='litter flux out (exudation, fragmentation and seed decay) in kg element', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_litter_out_elem)

          call this%set_history_var(vname='FATES_SEED_BANK_EL', units='kg m-2',      &
               long='element-level total seed mass of all PFTs in kg element per m2', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_seed_bank_elem)

          call this%set_history_var(vname='FATES_SEEDS_IN_LOCAL_EL',                 &
               units='kg m-2 s-1',                                                   &
               long='within-site, element-level seed production rate in kg element per m2 per second', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_seeds_in_local_elem)

          call this%set_history_var(vname='FATES_SEEDS_IN_EXTERN_EL',                &
               units='kg m-2 s-1', long='external seed influx rate in kg element per m2 per second', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_seeds_in_extern_elem)

          call this%set_history_var(vname='FATES_SEED_GERM_EL', units='kg m-2',  &
               long='element-level total germinated seed mass of all PFTs in kg element per m2', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_seed_germ_elem)

          call this%set_history_var(vname='FATES_SEED_DECAY_EL', units='kg m-2 s-1', &
               long='seed mass decay (germinated and un-germinated) in kg element per m2 per second', &
               use_default='active', avgflag='A', vtype=site_elem_r8,                &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_seed_decay_elem)

          ! SITE LEVEL CARBON STATE VARIABLES

          call this%set_history_var(vname='FATES_STOREC_TF_USTORY_SZPF', units='kg kg-1',         &
               long='Storage C fraction of target by size x pft, in the understory', use_default='inactive',          &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,   &
               ivar=ivar, initialize=initialize_variables, index = ih_storectfrac_ustory_scpf )


          call this%set_history_var(vname='FATES_STOREC_TF_CANOPY_SZPF', units='kg kg-1',         &
               long='Storage C fraction of target by size x pft, in the canopy', use_default='inactive',          &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,   &
               ivar=ivar, initialize=initialize_variables, index = ih_storectfrac_canopy_scpf )


          call this%set_history_var(vname='FATES_FROOTC_SL', units='kg m-3',                        &
               long='Total carbon in live plant fine-roots over depth', use_default='active', &
               avgflag='A', vtype=site_soil_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,         &
               ivar=ivar, initialize=initialize_variables, index = ih_fnrtc_sl )


          ! Output specific to the chemical species dynamics used (parteh)
          call this%set_history_var(vname='FATES_L2FR_CANOPY_REC_PF', units='kg kg-1', &
               long='The leaf to fineroot biomass multiplier for recruits (canopy)',   & 
               use_default='active', &
               avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,    &
               ivar=ivar, initialize=initialize_variables, index = ih_recl2fr_canopy_pf)
          
          call this%set_history_var(vname='FATES_L2FR_USTORY_REC_PF', units='kg kg-1',                   &
               long='The leaf to fineroot biomass multiplier for recruits (understory)', & 
               use_default='active', &
               avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,    &
               ivar=ivar, initialize=initialize_variables, index = ih_recl2fr_ustory_pf)
          
             !call this%set_history_var(vname='FATES_L2FR_CLSZPF', units='kg kg-1',                   &
             !     long='The leaf to fineroot biomass multiplier for target allometry', & 
             !     use_default='inactive', &
             !     avgflag='A', vtype=site_clscpf_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,    &
             !     ivar=ivar, initialize=initialize_variables, index = ih_l2fr_clscpf)

          nitrogen_active_if1: if(any(element_list(:)==nitrogen_element)) then

             call this%set_history_var(vname='FATES_NH4UPTAKE_SZPF',                 &
                  units='kg m-2 s-1',                                                &
                  long='ammonium uptake rate by plants by size-class x pft in kg NH4 per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_nh4uptake_scpf)

             call this%set_history_var(vname='FATES_NO3UPTAKE_SZPF',                 &
                  units='kg m-2 s-1',                                                &
                  long='nitrate uptake rate by plants by size-class x pft in kg NO3 per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_no3uptake_scpf)

             call this%set_history_var(vname='FATES_NEFFLUX_SZPF', units='kg m-2 s-1', &
                  long='nitrogen efflux, root to soil, by size-class x pft in kg N per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_nefflux_scpf)

             call this%set_history_var(vname='FATES_NDEMAND_SZPF', units='kg m-2 s-1', &
                  long='plant N need (algorithm dependent), by size-class x pft in kg N per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_ndemand_scpf)
             
             call this%set_history_var(vname='FATES_NFIX_SYM_SZPF', units='kg m-2 s-1', &
                  long='symbiotic dinitrogen fixation, by size-class x pft in kg N per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_nfix_scpf)
             
             call this%set_history_var(vname='FATES_VEGN_SZPF', units='kg m-2',      &
                  long='total (live) vegetation nitrogen mass by size-class x pft in kg N per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_totvegn_scpf)

             call this%set_history_var(vname='FATES_LEAFN_SZPF', units='kg m-2',     &
                  long='leaf nitrogen mass by size-class x pft in kg N per m2',      &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_leafn_scpf)

             call this%set_history_var(vname='FATES_FROOTN_SZPF', units='kg m-2',    &
                  long='fine-root nitrogen mass by size-class x pft in kg N per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_fnrtn_scpf)

             call this%set_history_var(vname='FATES_SAPWOODN_SZPF', units='kg m-2',  &
                  long='sapwood nitrogen mass by size-class x pft in kg N per m2',   &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_sapwn_scpf)

             call this%set_history_var(vname='FATES_STOREN_SZPF', units='kg m-2',    &
                  long='storage nitrogen mass by size-class x pft in kg N per m2',   &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_storen_scpf)

             call this%set_history_var(vname='FATES_STOREN_TF_CANOPY_SZPF',          &
                  units='1',                                                         &
                  long='storage nitrogen fraction (0-1) of target, in canopy, by size-class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_storentfrac_canopy_scpf)

             call this%set_history_var(vname='FATES_STOREN_TF_USTORY_SZPF',      &
                  units='1',                                                         &
                  long='storage nitrogen fraction (0-1) of target, in understory, by size-class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables,                                   &
                  index = ih_storentfrac_understory_scpf)

             call this%set_history_var(vname='FATES_REPRON_SZPF', units='kg m-2',    &
                  long='reproductive nitrogen mass (on plant) by size-class x pft in kg N per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_repron_scpf)

          end if nitrogen_active_if1

          phosphorus_active_if1: if(any(element_list(:)==phosphorus_element)) then

             call this%set_history_var(vname='FATES_VEGP_SZPF', units='kg m-2',      &
                  long='total (live) vegetation phosphorus mass by size-class x pft in kg P per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_totvegp_scpf)

             call this%set_history_var(vname='FATES_LEAFP_SZPF', units='kg m-2', &
                  long='leaf phosphorus mass by size-class x pft', use_default='inactive', &
                  avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_leafp_scpf )

             call this%set_history_var(vname='FATES_FROOTP_SZPF', units='kg m-2',    &
                  long='fine-root phosphorus mass by size-class x pft in kg P per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_fnrtp_scpf)

             call this%set_history_var(vname='FATES_SAPWOODP_SZPF', units='kg m-2',  &
                  long='sapwood phosphorus mass by size-class x pft in kg P per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_sapwp_scpf)

             call this%set_history_var(vname='FATES_STOREP_SZPF', units='kg m-2',    &
                  long='storage phosphorus mass by size-class x pft in kg P per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_storep_scpf)

             call this%set_history_var(vname='FATES_STOREP_TF_CANOPY_SZPF',          &
                  units='1',                                                         &
                  long='storage phosphorus fraction (0-1) of target, in canopy, by size-class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_storeptfrac_canopy_scpf)

             call this%set_history_var(vname='FATES_STOREP_TF_USTORY_SZPF',          &
                  units='1',                                                         &
                  long='storage phosphorus fraction (0-1) of target, in understory, by size-class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables,                                   &
                  index = ih_storeptfrac_understory_scpf)

             call this%set_history_var(vname='FATES_REPROP_SZPF', units='kg m-2',    &
                  long='reproductive phosphorus mass (on plant) by size-class x pft in kg P per m2', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_reprop_scpf)

             call this%set_history_var(vname='FATES_PUPTAKE_SZPF',                   &
                  units='kg m-2 s-1',                                                &
                  long='phosphorus uptake rate by plants, by size-class x pft in kg P per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_puptake_scpf)

             call this%set_history_var(vname='FATES_PEFFLUX_SZPF',                   &
                  units='kg m-2 s-1',                                                &
                  long='phosphorus efflux, root to soil, by size-class x pft in kg P per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_pefflux_scpf)
             
             call this%set_history_var(vname='FATES_PDEMAND_SZPF', units='kg m-2 s-1', &
                  long='plant P need (algorithm dependent), by size-class x pft in kg P per m2 per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_nflx_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_pdemand_scpf)
             
          end if phosphorus_active_if1

          call this%set_history_var(vname='FATES_CROWNAREA_CLLL', units='m2 m-2',    &
               long='area fraction of the total ground occupied by each canopy-leaf layer', &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_crownarea_si_cnlf)

          call this%set_history_var(vname='FATES_CROWNAREA_CL', units='m2 m-2',      &
               long='area fraction of the canopy footprint occupied by each canopy-leaf layer', use_default='active',   &
               avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,             &
               ivar=ivar, initialize=initialize_variables, index = ih_crownarea_cl)


          call this%set_history_var(vname='FATES_MORTALITY_CFLUX_PF', units='kg m-2 s-1',    &
               long='PFT-level flux of biomass carbon from live to dead pool from mortality', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_mortality_carbonflux_si_pft)

          call this%set_history_var(vname='FATES_MORTALITY_FIRE_CFLUX_PF', units='kg m-2 s-1',    &
               long='PFT-level flux of biomass carbon from live to dead pool from fire mortality', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_firemortality_carbonflux_si_pft)

          call this%set_history_var(vname='FATES_MORTALITY_HYDRO_CFLUX_PF', units='kg m-2 s-1',    &
               long='PFT-level flux of biomass carbon from live to dead pool from hydraulic failure mortality', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_hydraulicmortality_carbonflux_si_pft)

          call this%set_history_var(vname='FATES_MORTALITY_CSTARV_CFLUX_PF', units='kg m-2 s-1',    &
               long='PFT-level flux of biomass carbon from live to dead pool from carbon starvation mortality (both continuous and termination)', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_cstarvmortality_carbonflux_si_pft)

          call this%set_history_var(vname='FATES_MORT_CSTARV_CONT_CFLUX_PF', units='kg m-2 s-1',    &
               long='PFT-level flux of biomass carbon from live to dead pool from carbon starvation mortality (Continuous-only, without termination)', &
               use_default='active', avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_cstarvmortality_continuous_carbonflux_si_pft)
          
          call this%set_history_var(vname='FATES_ABOVEGROUND_MORT_SZPF', units='kg m-2 s-1',    &
               long='Aboveground flux of carbon from AGB to necromass due to mortality', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_abg_mortality_cflux_si_scpf)

          call this%set_history_var(vname='FATES_ABOVEGROUND_PROD_SZPF', units='kg m-2 s-1',    &
               long='Aboveground carbon productivity', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index=ih_abg_productivity_cflux_si_scpf)
          ! size class by age dimensioned variables

          call this%set_history_var(vname='FATES_NPLANT_SZAP', units = 'm-2',        &
               long='number of plants per m2 in each size x age class',             &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_si_scag)

          call this%set_history_var(vname='FATES_NPLANT_CANOPY_SZAP', units = 'm-2', &
               long='number of plants per m2 in canopy in each size x age class',   &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_canopy_si_scag)

          call this%set_history_var(vname='FATES_NPLANT_USTORY_SZAP',                &
               units = 'm-2',                                                       &
               long='number of plants per m2 in understory in each size x age class', &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_understory_si_scag)

          call this%set_history_var(vname='FATES_DDBH_CANOPY_SZAP',                  &
               units = 'm m-2 yr-1',                                                &
               long='growth rate of canopy plants in meters DBH per m2 per year in canopy in each size x age class', &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_canopy_si_scag)

          call this%set_history_var(vname='FATES_DDBH_USTORY_SZAP',                  &
               units = 'm m-2 yr-1',                                                &
               long='growth rate of understory plants in meters DBH per m2 per year in each size x age class', &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_understory_si_scag)

          call this%set_history_var(vname='FATES_MORTALITY_CANOPY_SZAP',             &
               units = 'm-2 yr-1',                                                  &
               long='mortality rate of canopy plants in number of plants per m2 per year in each size x age class', &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_mortality_canopy_si_scag)

          call this%set_history_var(vname='FATES_MORTALITY_USTORY_SZAP',             &
               units = 'm-2 yr-1',                                                  &
               long='mortality rate of understory plants in number of plants per m2 per year in each size x age class', &
               use_default='inactive', avgflag='A', vtype=site_scag_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_mortality_understory_si_scag)

          ! size x age x pft dimensioned

          call this%set_history_var(vname='FATES_NPLANT_SZAPPF',units = 'm-2',       &
               long='number of plants per m2 in each size x age x pft class',       &
               use_default='inactive', avgflag='A', vtype=site_scagpft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_si_scagpft)

          ! age x pft dimensioned
          call this%set_history_var(vname='FATES_NPP_APPF',units = 'kg m-2 s-1',     &
               long='NPP per PFT in each age bin in kg carbon per m2 per second',   &
               use_default='inactive', avgflag='A', vtype=site_agepft_r8,           &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_npp_si_agepft)

          call this%set_history_var(vname='FATES_NPP_AP', units='kg m-2 s-1',        &
               long='net primary productivity by age bin in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_si_age)

          call this%set_history_var(vname='FATES_VEGC_APPF',units = 'kg m-2',        &
               long='biomass per PFT in each age bin in kg carbon per m2',          &
               use_default='inactive', avgflag='A', vtype=site_agepft_r8,           &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_biomass_si_agepft)

          call this%set_history_var(vname='FATES_SCORCH_HEIGHT_APPF',units = 'm',    &
               long='SPITFIRE flame Scorch Height (calculated per PFT in each patch age bin)', &
               use_default='inactive', avgflag='A', vtype=site_agepft_r8,           &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_scorch_height_si_agepft)


          ! Carbon Flux (grid dimension x scpf) (THESE ARE DEFAULT INACTIVE!!!
          !                                     (BECAUSE THEY TAKE UP SPACE!!!
          ! ===================================================================================

          call this%set_history_var(vname='FATES_GPP_SZPF', units='kg m-2 s-1',      &
               long='gross primary production by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_gpp_si_scpf)

          call this%set_history_var(vname='FATES_GPP_CANOPY_SZPF',                   &
               units='kg m-2 s-1',                                                  &
               long='gross primary production of canopy plants by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_gpp_canopy_si_scpf)

          call this%set_history_var(vname='FATES_AUTORESP_CANOPY_SZPF',             &
               units='kg m-2 s-1',                                                  &
               long='autotrophic respiration of canopy plants by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_canopy_si_scpf)

          call this%set_history_var(vname='FATES_GPP_USTORY_SZPF',               &
               units='kg m-2 s-1',                                                  &
               long='gross primary production of understory plants by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_gpp_understory_si_scpf)

          call this%set_history_var(vname='FATES_AUTORESP_USTORY_SZPF',         &
               units='kg m-2 s-1',                                                  &
               long='autotrophic respiration of understory plants by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_understory_si_scpf)

          call this%set_history_var(vname='FATES_NPP_SZPF', units='kg m-2 s-1',      &
               long='total net primary production by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_npp_totl_si_scpf)

          call this%set_history_var(vname='FATES_LEAF_ALLOC_SZPF', units='kg m-2 s-1', &
               long='allocation to leaves by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_npp_leaf_si_scpf)

          call this%set_history_var(vname='FATES_SEED_ALLOC_SZPF', units='kg m-2 s-1',  &
               long='allocation to seeds by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                  &
               initialize=initialize_variables, index = ih_npp_seed_si_scpf)

          call this%set_history_var(vname='FATES_FROOT_ALLOC_SZPF',                   &
               units='kg m-2 s-1',                                                   &
               long='allocation to fine roots by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                  &
               initialize=initialize_variables, index = ih_npp_fnrt_si_scpf)

          call this%set_history_var(vname='FATES_BGSAPWOOD_ALLOC_SZPF',               &
               units='kg m-2 s-1',                                                   &
               long='allocation to below-ground sapwood by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_bgsw_si_scpf)

          call this%set_history_var(vname='FATES_BGSTRUCT_ALLOC_SZPF', units='kg m-2 s-1', &
               long='allocation to below-ground structural (deadwood) by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_bgdw_si_scpf)

          call this%set_history_var(vname='FATES_AGSAPWOOD_ALLOC_SZPF',               &
               units='kg m-2 s-1',                                                   &
               long='allocation to above-ground sapwood by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_agsw_si_scpf)

          call this%set_history_var(vname = 'FATES_AGSTRUCT_ALLOC_SZPF',              &
               units='kg m-2 s-1',                                                   &
               long='allocation to above-ground structural (deadwood) by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_agdw_si_scpf)

          call this%set_history_var(vname = 'FATES_STORE_ALLOC_SZPF',                 &
               units='kg m-2 s-1',                                                   &
               long='allocation to storage C by pft/size in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_stor_si_scpf)

          call this%set_history_var(vname='FATES_DDBH_SZPF', units = 'm m-2 yr-1',   &
               long='diameter growth increment by pft/size', use_default='inactive', &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',                 &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                &
               index = ih_ddbh_si_scpf)

          call this%set_history_var(vname='FATES_GROWTHFLUX_SZPF',                   &
               units = 'm-2 yr-1',                                                  &
               long='flux of individuals into a given size class bin via growth and recruitment', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_growthflux_si_scpf)

          call this%set_history_var(vname='FATES_GROWTHFLUX_FUSION_SZPF',            &
               units = 'm-2 yr-1',                                                  &
               long='flux of individuals into a given size class bin via fusion',   &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_growthflux_fusion_si_scpf)

          call this%set_history_var(vname='FATES_DDBH_CANOPY_SZPF',                  &
               units = 'm m-2 yr-1',                                                &
               long='diameter growth increment by pft/size',                        &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_canopy_si_scpf)

          call this%set_history_var(vname='FATES_DDBH_USTORY_SZPF',              &
               units = 'm m-2 yr-1',                                                &
               long='diameter growth increment by pft/size',                        &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_understory_si_scpf)

          call this%set_history_var(vname='FATES_BASALAREA_SZPF', units = 'm2 m-2',  &
               long='basal area by pft/size', use_default='inactive',               &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,       &
               ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scpf)

          call this%set_history_var(vname='FATES_SAPWOOD_AREA_SZPF', units = 'm2 m-2',  &
               long='sapwood area by pft/size', use_default='inactive',               &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,       &
               ivar=ivar, initialize=initialize_variables, index = ih_sapwood_area_scpf)

          call this%set_history_var(vname='FATES_VEGC_ABOVEGROUND_SZPF',             &
               units = 'kg m-2',                                                     &
               long='aboveground biomass by pft/size in kg carbon per m2',           &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_agb_si_scpf)

          call this%set_history_var(vname='FATES_NPLANT_SZPF', units = 'm-2',        &
               long='stem number density by pft/size', use_default='inactive',      &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',                 &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                &
               index = ih_nplant_si_scpf)

          call this%set_history_var(vname='FATES_NPLANT_ACPF', units = 'm-2',        &
               long='stem number density by pft and age class',                      &
               use_default='inactive', avgflag='A', vtype=site_coage_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_nplant_si_capf)

          call this%set_history_var(vname='FATES_MORTALITY_BACKGROUND_SZPF',         &
               units = 'm-2 yr-1',                                                  &
               long='background mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m1_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_HYDRAULIC_SZPF',          &
               units = 'm-2 yr-1',                                                  &
               long='hydraulic mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m2_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_CSTARV_SZPF',             &
               units = 'm-2 yr-1',                                                  &
               long='carbon starvation mortality by pft/size in number of plants per m2 per year (both continous and termination)', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_IMPACT_SZPF',             &
               units = 'm-2 yr-1',                                                  &
               long='impact mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m4_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_FIRE_SZPF',               &
               units = 'm-2 yr-1',                                                  &
               long='fire mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m5_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_CROWNSCORCH_SZPF',        &
               units = 'm-2 yr-1',                                                  &
               long='fire mortality from crown scorch by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_crownfiremort_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_CAMBIALBURN_SZPF',        &
               units = 'm-2 yr-1',                                                  &
               long='fire mortality from cambial burn by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cambialfiremort_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_TERMINATION_SZPF',        &
               units = 'm-2 yr-1',                                                  &
               long='termination mortality (excluding C-starvation) by pft/size in number pf plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m6_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_LOGGING_SZPF',            &
               units = 'm-2 yr-1',                                                  &
               long='logging mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive',           &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,       &
               ivar=ivar, initialize=initialize_variables, index = ih_m7_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_FREEZING_SZPF',           &
               units = 'm-2 yr-1',                                                  &
               long='freezing mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m8_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_SENESCENCE_SZPF',         &
               units = 'm-2 yr-1',                                                 &
               long='senescence mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m9_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_AGESCEN_SZPF',            &
               units = 'm-2 yr-1',                                                   &
               long='age senescence mortality by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype =site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_m10_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_AGESCEN_ACPF',            &
               units='m-2 yr-1',                                                     &
               long='age senescence mortality by pft/cohort age in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype =site_coage_pft_r8,        &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index =ih_m10_si_capf)

          call this%set_history_var(vname='FATES_MORTALITY_CANOPY_SZPF',             &
               units = 'm-2 yr-1',                                                  &
               long='total mortality of canopy plants by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_mortality_canopy_si_scpf)

          call this%set_history_var(vname='FATES_M3_MORTALITY_CANOPY_SZPF',          &
               units = 'm-2 yr-1',                                &
               long='C starvation mortality of canopy plants by pft/size',          &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_mortality_canopy_si_scpf )

          call this%set_history_var(vname='FATES_M3_MORTALITY_USTORY_SZPF',          &
               units = 'm-2 yr-1',                                                   &
               long='C starvation mortality of understory plants by pft/size',      &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_mortality_understory_si_scpf )


          call this%set_history_var(vname='FATES_C13DISC_SZPF', units = 'per mil',   &
               long='C13 discrimination by pft/size',use_default='inactive',         &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',                  &
               upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,                 &
               index = ih_c13disc_si_scpf)
          
          call this%set_history_var(vname='FATES_STOREC_CANOPY_SZPF', units = 'kg m-2', &
               long='biomass in storage pools of canopy plants by pft/size in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bstor_canopy_si_scpf)

          call this%set_history_var(vname='FATES_LEAFC_CANOPY_SZPF',                 &
               units = 'kg m-2',                                                    &
               long='biomass in leaves of canopy plants by pft/size in kg carbon per m2', &
               use_default='inactive', &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,       &
               ivar=ivar, initialize=initialize_variables,                          &
               index = ih_bleaf_canopy_si_scpf)

          call this%set_history_var(vname='FATES_LAI_CANOPY_SZPF',                   &
               units = 'm2 m-2',                                                    &
               long='Leaf area index (LAI) of canopy plants by pft/size',           &
               use_default='inactive', &
               avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,       &
               ivar=ivar, initialize=initialize_variables,                          &
               index = ih_lai_canopy_si_scpf )          

          call this%set_history_var(vname='FATES_CROWNAREA_CANOPY_SZPF',             &
               units = 'm2 m-2',                                                    &
               long='Total crown area of canopy plants by pft/size',                &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_crownarea_canopy_si_scpf )

          call this%set_history_var(vname='FATES_CROWNAREA_USTORY_SZPF',             &
               units = 'm2 m-2',                                                    &
               long='Total crown area of understory plants by pft/size',            &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_crownarea_understory_si_scpf )
          
          call this%set_history_var(vname='FATES_NPLANT_CANOPY_SZPF', units = 'm-2', &
               long='number of canopy plants by size/pft per m2',                   &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_canopy_si_scpf)

          call this%set_history_var(vname='FATES_MORTALITY_USTORY_SZPF',         &
               units = 'm-2 yr-1',                                                  &
               long='total mortality of understory plants by pft/size in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_mortality_understory_si_scpf)

          call this%set_history_var(vname='FATES_STOREC_USTORY_SZPF',            &
               units = 'kg m-2',                                                    &
               long='biomass in storage pools of understory plants by pft/size in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bstor_understory_si_scpf)

          call this%set_history_var(vname='FATES_LEAFC_USTORY_SZPF',             &
               units = 'kg m-2',                                                    &
               long='biomass in leaves of understory plants by pft/size in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bleaf_understory_si_scpf)

          call this%set_history_var(vname='FATES_LAI_USTORY_SZPF',                   &
               units = 'm2 m-2',                                                    &
               long='Leaf area index (LAI) of understory plants by pft/size',       &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_lai_understory_si_scpf )

          call this%set_history_var(vname='FATES_NPLANT_USTORY_SZPF',            &
               units = 'm-2',                                                       &
               long='density of understory plants by pft/size in number of plants per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_understory_si_scpf)

          call this%set_history_var(vname='FATES_CWD_ABOVEGROUND_DC', units='kg m-2', &
               long='debris class-level aboveground coarse woody debris stocks in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_ag_si_cwdsc)

          call this%set_history_var(vname='FATES_CWD_BELOWGROUND_DC', units='kg m-2', &
               long='debris class-level belowground coarse woody debris stocks in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_bg_si_cwdsc)

          call this%set_history_var(vname='FATES_CWD_ABOVEGROUND_IN_DC',             &
               units='kg m-2 s-1',                                                  &
               long='debris class-level aboveground coarse woody debris input in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_ag_in_si_cwdsc)

          call this%set_history_var(vname='FATES_CWD_BELOWGROUND_IN_DC',             &
               units='kg m-2 s-1',                                                  &
               long='debris class-level belowground coarse woody debris input in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_bg_in_si_cwdsc)

          call this%set_history_var(vname='FATES_CWD_ABOVEGROUND_OUT_DC',            &
               units='kg m-2 s-1',                                                  &
               long='debris class-level aboveground coarse woody debris output in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_ag_out_si_cwdsc)

          call this%set_history_var(vname='FATES_CWD_BELOWGROUND_OUT_DC',            &
               units='kg m-2 s-1',                                                  &
               long='debris class-level belowground coarse woody debris output in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_cwdsc_r8,            &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_bg_out_si_cwdsc)


          ! size-class only variables

          call this%set_history_var(vname='FATES_DDBH_CANOPY_SZ',                    &
               units = 'm m-2 yr-1', long='diameter growth increment by size of canopy plants', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_canopy_si_scls)

          call this%set_history_var(vname='FATES_DDBH_USTORY_SZ',                &
               units = 'm m-2 yr-1', long='diameter growth increment by size of understory plants', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ddbh_understory_si_scls)

          call this%set_history_var(vname='FATES_YESTCANLEV_CANOPY_SZ',              &
               units = 'm-2',                                                       &
               long='yesterdays canopy level for canopy plants by size class in number of plants per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_yesterdaycanopylevel_canopy_si_scls)

          call this%set_history_var(vname='FATES_YESTCANLEV_USTORY_SZ',          &
               units = 'm-2',                                                       &
               long='yesterdays canopy level for understory plants by size class in number of plants per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_yesterdaycanopylevel_understory_si_scls)

          call this%set_history_var(vname='FATES_BASALAREA_SZ', units = 'm2 m-2',    &
               long='basal area by size class', use_default='active',               &
               avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,           &
               ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scls)

          call this%set_history_var(vname='FATES_VEGC_ABOVEGROUND_SZ',               &
               units = 'kg m-2',                                                   &
               long='aboveground biomass by size class in kg carbon per m2',        &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_agb_si_scls)

          call this%set_history_var(vname='FATES_VEGC_SZ', units = 'kg m-2',         &
               long='total biomass by size class in kg carbon per m2',              &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_biomass_si_scls)

          call this%set_history_var(vname='FATES_DEMOTION_RATE_SZ',                  &
               units = 'm-2 yr-1',                                                  &
               long='demotion rate from canopy to understory by size class in number of plants per m2 per year', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_demotion_rate_si_scls)

          call this%set_history_var(vname='FATES_PROMOTION_RATE_SZ',                 &
               units = 'm-2 yr-1',                                                  &
               long='promotion rate from understory to canopy by size class',       &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_promotion_rate_si_scls)

          call this%set_history_var(vname='FATES_NPLANT_CANOPY_SZ',                  &
               units = 'm-2',               &
               long='number of canopy plants per m2 by size class',                 &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_canopy_si_scls)

          call this%set_history_var(vname='FATES_LAI_CANOPY_SZ', units = 'm2 m-2',   &
               long='leaf area index (LAI) of canopy plants by size class',         &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_lai_canopy_si_scls)

          call this%set_history_var(vname='FATES_SAI_CANOPY_SZ', units = 'm2 m-2',   &
               long='stem area index (SAI) of canopy plants by size class',         &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_sai_canopy_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_CANOPY_SZ',               &
               units = 'm-2 yr-1',                                                  &
               long='total mortality of canopy trees by size class in number of plants per m2', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_mortality_canopy_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_CANOPY_SE_SZ',               &
               units = 'm-2 yr-1',                                                  &
               long='total mortality of canopy trees by size class in number of plants per m2, secondary patches', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_mortality_canopy_secondary_si_scls)

          call this%set_history_var(vname='FATES_NPLANT_USTORY_SZ',              &
               units = 'm-2',                                                       &
               long='number of understory plants per m2 by size class',             &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_nplant_understory_si_scls)

          call this%set_history_var(vname='FATES_M3_MORTALITY_CANOPY_SZ',            &
               units = 'm-2 yr-1',                                                   &
               long='C starvation mortality of canopy plants by size',              &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_mortality_canopy_si_scls )

          call this%set_history_var(vname='FATES_M3_MORTALITY_USTORY_SZ',            &
               units = 'm-2 yr-1',                                                   &
               long='C starvation mortality of understory plants by size',          &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_mortality_understory_si_scls )

          call this%set_history_var(vname='FATES_LAI_USTORY_SZ',                 &
               units = 'm2 m-2',                                                    &
               long='leaf area index (LAI) of understory plants by size class',     &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_lai_understory_si_scls)

          call this%set_history_var(vname='FATES_SAI_USTORY_SZ',                 &
               units = 'm2 m-2',                                                    &
               long='stem area index (SAI) of understory plants by size class',     &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_sai_understory_si_scls)

          call this%set_history_var(vname='FATES_NPLANT_SZ', units = 'm-2',          &
               long='number of plants per m2 by size class', use_default='active',  &
               avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', upfreq=group_dyna_complx,           &
               ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scls)

          call this%set_history_var(vname='FATES_NPLANT_AC', units = 'm-2',          &
               long='number of plants per m2 by cohort age class',                   &
               use_default='active', avgflag='A', vtype=site_coage_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                  &
               initialize=initialize_variables, index = ih_nplant_si_cacls)

          call this%set_history_var(vname='FATES_MORTALITY_BACKGROUND_SZ',           &
               units = 'm-2 yr-1',                                                  &
               long='background mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m1_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_HYDRAULIC_SZ',            &
               units = 'm-2 yr-1',                                                  &
               long='hydraulic mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m2_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_CSTARV_SZ',               &
               units = 'm-2 yr-1',                                                  &
               long='carbon starvation mortality by size in number of plants per m2 per year (both continous and termination)', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m3_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_IMPACT_SZ',               &
               units = 'm-2 yr-1',                                                  &
               long='impact mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m4_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_FIRE_SZ',                 &
               units = 'm-2 yr-1',                                                  &
               long='fire mortality by size in number of plants per m2 per year',   &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m5_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_TERMINATION_SZ',          &
               units = 'm-2 yr-1',                                                  &
               long='termination mortality (excluding C-starvation) by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m6_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_LOGGING_SZ',              &
               units = 'm-2 yr-1',                                                  &
               long='logging mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m7_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_FREEZING_SZ',             &
               units = 'm-2 yr-1',                                                  &
               long='freezing mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m8_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_SENESCENCE_SZ',           &
               units = 'm-2 yr-1',                                                  &
               long='senescence mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m9_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_AGESCEN_SZ',              &
               units = 'm-2 yr-1',                                                  &
               long='age senescence mortality by size in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m10_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_AGESCEN_AC',              &
               units = 'm-2 yr-1',                                                  &
               long='age senescence mortality by cohort age in number of plants per m2 per year', &
               use_default='active', avgflag='A', vtype=site_coage_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_m10_si_cacls)

          call this%set_history_var(vname='FATES_NPP_CANOPY_SZ', units = 'kg m-2 s-1', &
               long='NPP of canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_carbon_balance_canopy_si_scls)

          call this%set_history_var(vname='FATES_NPP_USTORY_SZ', units = 'kg m-2 s-1', &
               long='NPP of understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_carbon_balance_understory_si_scls)

          call this%set_history_var(vname='FATES_MORTALITY_USTORY_SZ',           &
               units = 'm-2 yr-1',                                                  &
               long='total mortality of understory trees by size class in individuals per m2 per year', &
               use_default='active', avgflag='A', vtype=site_size_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_mortality_understory_si_scls)

          call this%set_history_var(vname='FATES_TRIMMING_CANOPY_SZ', units = 'm-2', &
               long='trimming term of canopy plants weighted by plant density, by size class', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_trimming_canopy_si_scls)

          call this%set_history_var(vname='FATES_TRIMMING_USTORY_SZ',            &
               units = 'm-2',                                                       &
               long='trimming term of understory plants weighted by plant density, by size class', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_trimming_understory_si_scls)

          call this%set_history_var(vname='FATES_CROWNAREA_CANOPY_SZ', units = 'm2 m-2', &
               long='total crown area of canopy plants by size class',              &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_crown_area_canopy_si_scls)

          call this%set_history_var(vname='FATES_CROWNAREA_USTORY_SZ', units = 'm2 m-2', &
               long='total crown area of understory plants by size class',          &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_crown_area_understory_si_scls)

          call this%set_history_var(vname='FATES_LEAFCTURN_CANOPY_SZ',               &
               units = 'kg m-2 s-1',                                                &
               long='leaf turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_leaf_md_canopy_si_scls)

          call this%set_history_var(vname='FATES_FROOTCTURN_CANOPY_SZ',              &
               units = 'kg m-2 s-1',                                                &
               long='fine root turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_root_md_canopy_si_scls)

          call this%set_history_var(vname='FATES_STORECTURN_CANOPY_SZ',              &
               units = 'kg m-2 s-1',                                                &
               long='storage turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bstore_md_canopy_si_scls)

          call this%set_history_var(vname='FATES_STRUCTCTURN_CANOPY_SZ',             &
               units = 'kg m-2 s-1',                                                &
               long='structural C turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bdead_md_canopy_si_scls)

          call this%set_history_var(vname='FATES_SAPWOODCTURN_CANOPY_SZ',            &
               units = 'kg m-2 s-1',                                                &
               long='sapwood turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bsw_md_canopy_si_scls)

          call this%set_history_var(vname='FATES_SEED_PROD_CANOPY_SZ',               &
               units = 'kg m-2 s-1',                                                &
               long='seed production of canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_seed_prod_canopy_si_scls)

          call this%set_history_var(vname='FATES_LEAF_ALLOC_CANOPY_SZ',               &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to leaves for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_leaf_canopy_si_scls)

          call this%set_history_var(vname='FATES_FROOT_ALLOC_CANOPY_SZ',              &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to fine root C for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_fnrt_canopy_si_scls)

          call this%set_history_var(vname='FATES_SAPWOOD_ALLOC_CANOPY_SZ',            &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to sapwood C for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_sapw_canopy_si_scls)

          call this%set_history_var(vname='FATES_STRUCT_ALLOC_CANOPY_SZ',             &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to structural C for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_dead_canopy_si_scls)

          call this%set_history_var(vname='FATES_SEED_ALLOC_CANOPY_SZ',               &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to reproductive C for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_seed_canopy_si_scls)

          call this%set_history_var(vname='FATES_STORE_ALLOC_CANOPY_SZ',              &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to storage C for canopy plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_stor_canopy_si_scls)

          call this%set_history_var(vname='FATES_LEAFCTURN_USTORY_SZ',           &
               units = 'kg m-2 s-1',                                                 &
               long='leaf turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_leaf_md_understory_si_scls)

          call this%set_history_var(vname='FATES_FROOTCTURN_USTORY_SZ',          &
               units = 'kg m-2 s-1',                                                &
               long='fine root turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_root_md_understory_si_scls)

          call this%set_history_var(vname='FATES_STORECTURN_USTORY_SZ',              &
               units = 'kg m-2 s-1',                                                &
               long='storage C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_bstore_md_understory_si_scls)

          call this%set_history_var(vname='FATES_STRUCTCTURN_USTORY_SZ',         &
               units = 'kg m-2 s-1',                                                &
               long='structural C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_bdead_md_understory_si_scls)

          call this%set_history_var(vname='FATES_SAPWOODCTURN_USTORY_SZ',        &
               units = 'kg m-2 s-1',                                                &
               long='sapwood C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_bsw_md_understory_si_scls)

          call this%set_history_var(vname='FATES_SEED_PROD_USTORY_SZ',           &
               units = 'kg m-2 s-1',                                                &
               long='seed production of understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_seed_prod_understory_si_scls)

          call this%set_history_var(vname='FATES_LEAF_ALLOC_USTORY_SZ',           &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to leaves for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_leaf_understory_si_scls)

          call this%set_history_var(vname='FATES_FROOT_ALLOC_USTORY_SZ',          &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to fine roots for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_fnrt_understory_si_scls)

          call this%set_history_var(vname='FATES_SAPWOOD_ALLOC_USTORY_SZ',        &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to sapwood C for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_sapw_understory_si_scls)

          call this%set_history_var(vname='FATES_STRUCT_ALLOC_USTORY_SZ',         &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to structural C for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_dead_understory_si_scls)

          call this%set_history_var(vname='FATES_SEED_ALLOC_USTORY_SZ',           &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to reproductive C for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_seed_understory_si_scls)

          call this%set_history_var(vname='FATES_STORE_ALLOC_USTORY_SZ',          &
               units = 'kg m-2 s-1',                                                 &
               long='allocation to storage C for understory plants by size class in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_npp_stor_understory_si_scls)

          ! CROWN DAMAGE VARIABLES
          if_crowndamage: if(hlm_use_tree_damage .eq. itrue) then 

             call this%set_history_var(vname='FATES_NPLANT_CDPF', units = 'm-2',     &
                  long='N. plants per damage x size x pft class', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_cdpf )

             call this%set_history_var(vname='FATES_NPLANT_CANOPY_CDPF', units = 'm-2',     &
                  long='N. plants per damage x size x pft class', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_cdpf )

             call this%set_history_var(vname='FATES_NPLANT_USTORY_CDPF', units = 'm-2',     &
                  long='N. plants in the understory per damage x size x pft class', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_cdpf )

             call this%set_history_var(vname='FATES_M3_CDPF', units = 'm-2 yr-1',          &
                  long='carbon starvation mortality by damaage/pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m3_si_cdpf )

             call this%set_history_var(vname='FATES_M11_SZPF', units = 'm-2 yr-1',         &
                  long='damage mortality by pft/size',use_default='inactive', &
                  avgflag='A', vtype =site_size_pft_r8, hlms='CLM:ALM',      &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m11_si_scpf )

             call this%set_history_var(vname='FATES_M11_CDPF', units = 'm-2 yr-1',          &
                  long='damage mortality by damaage/pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m11_si_cdpf )

             call this%set_history_var(vname='FATES_MORTALITY_CDPF', units = 'm-2 yr-1',          &
                  long='mortality by damage class by size by pft', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',   &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_mortality_si_cdpf )

             call this%set_history_var(vname='FATES_M3_MORTALITY_CANOPY_CDPF', units = 'm-2 yr-1',          &
                  long='C starvation mortality of canopy plants by damage/pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m3_mortality_canopy_si_cdpf )

             call this%set_history_var(vname='FATES_M3_MORTALITY_USTORY_CDPF', units = 'm-2 yr-1',          &
                  long='C starvation mortality of understory plants by pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m3_mortality_understory_si_cdpf )

             call this%set_history_var(vname='FATES_M11_MORTALITY_CANOPY_CDPF', units = 'm-2 yr-1',          &
                  long='damage mortality of canopy plants by damage/pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m11_mortality_canopy_si_cdpf )

             call this%set_history_var(vname='FATES_M11_MORTALITY_USTORY_CDPF', units = 'm-2 yr-1',          &
                  long='damage mortality of understory plants by pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_m11_mortality_understory_si_cdpf )

             call this%set_history_var(vname='FATES_MORTALITY_CANOPY_CDPF', units = 'm-2 yr-1',          &
                  long='mortality of canopy plants by damage/pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',   &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_cdpf )

             call this%set_history_var(vname='FATES_MORTALITY_USTORY_CDPF', units = 'm-2 yr-1',          &
                  long='mortality of understory plants by pft/size', use_default='inactive', &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_cdpf )

             call this%set_history_var(vname='FATES_DDBH_CDPF', units = 'm m-2 yr-1',               &
                  long='ddbh annual increment growth by damage x size pft', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',    &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_si_cdpf )

             call this%set_history_var(vname='FATES_DDBH_CANOPY_CDPF', units = 'm m-2 yr-1',               &
                  long='ddbh annual canopy increment growth by damage x size pft', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',   &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_cdpf )

             call this%set_history_var(vname='FATES_DDBH_USTORY_CDPF', units = 'm m-2 yr-1',               &
                  long='ddbh annual understory increment growth by damage x size pft', use_default='inactive',   &
                  avgflag='A', vtype=site_cdpf_r8, hlms='CLM:ALM',     &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_cdpf )

          end if if_crowndamage


          call this%set_history_var(vname='FATES_FIRE_FLUX_EL', units='kg m-2 s-1',  &
               long='loss to atmosphere from fire by element in kg element per m2 per s', &
               use_default='active', avgflag='A', vtype=site_elem_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_burn_flux_elem)

          if((hlm_use_ed_st3 .eq. ifalse) .and. (hlm_use_sp .eq. ifalse))then
             call this%set_history_var(vname='FATES_ERROR_EL', units='kg s-1',          &
                  long='total mass-balance error in kg per second by element',          &
                  use_default='active', avgflag='A', vtype=site_elem_r8,                &
                  hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
                  index = ih_err_fates_elem)
             
             call this%set_history_var(vname='FATES_INTERR_LIVEVEG_EL',units='kg m-2',             &
                  long='Bias error between integrated flux and (minus) state in live vegetation ', &
                  use_default='active', avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM',           &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,            &
                  index = ih_interr_liveveg_elem)
             
             call this%set_history_var(vname='FATES_INTERR_LITTER_EL',units='kg m-2',             &
                  long='Bias error between integrated flux and (minus) state in litter ', &
                  use_default='active', avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM',           &
                  upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables,            &
                  index = ih_interr_litter_elem)
          end if
             
          call this%set_history_var(vname='FATES_LITTER_AG_FINE_EL', units='kg m-2', &
               long='mass of aboveground litter in fines (leaves, nonviable seed) by element', &
               use_default='active', avgflag='A', vtype=site_elem_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_fines_ag_elem)

          call this%set_history_var(vname='FATES_LITTER_BG_FINE_EL', units='kg m-2', &
               long='mass of belowground litter in fines (fineroots) by element',   &
               use_default='active', avgflag='A', vtype=site_elem_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_fines_bg_elem)

          call this%set_history_var(vname='FATES_LITTER_BG_CWD_EL', units='kg m-2',  &
               long='mass of belowground litter in coarse woody debris (coarse roots) by element', &
               use_default='active', avgflag='A', vtype=site_elem_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_bg_elem)

          call this%set_history_var(vname='FATES_LITTER_AG_CWD_EL', units='kg m-2',  &
               long='mass of aboveground litter in coarse woody debris (trunks/branches/twigs) by element', &
               use_default='active', avgflag='A', vtype=site_elem_r8,               &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_ag_elem)

          call this%set_history_var(vname='FATES_LITTER_CWD_ELDC', units='kg m-2',   &
               long='total mass of litter in coarse woody debris by element and coarse woody debris size', &
               use_default='active', avgflag='A', vtype=site_elcwd_r8,              &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_cwd_elcwd)

          ! Mass states C/N/P SCPF dimensions
          ! CARBON
          call this%set_history_var(vname='FATES_VEGC_SZPF', units='kg m-2',         &
               long='total vegetation biomass in live plants by size-class x pft in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_totvegc_scpf)

          call this%set_history_var(vname='FATES_LEAFC_SZPF', units='kg m-2',        &
               long='leaf carbon mass by size-class x pft in kg carbon per m2',      &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_leafc_scpf)

          call this%set_history_var(vname='FATES_FROOTC_SZPF', units='kg m-2',       &
               long='fine-root carbon mass by size-class x pft in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_fnrtc_scpf)

          call this%set_history_var(vname='FATES_SAPWOODC_SZPF', units='kg m-2',     &
               long='sapwood carbon mass by size-class x pft in kg carbon per m2',   &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_sapwc_scpf)

          call this%set_history_var(vname='FATES_STOREC_SZPF', units='kg m-2',       &
               long='storage carbon mass by size-class x pft in kg carbon per m2',   &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_storec_scpf)

          call this%set_history_var(vname='FATES_REPROC_SZPF', units='kg m-2',       &
               long='reproductive carbon mass (on plant) by size-class x pft in kg carbon per m2', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,          &
               hlms='CLM:ALM', upfreq=group_dyna_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_reproc_scpf)



       end if if_dyn1
    end if if_dyn0

    !HERE


    if_hifrq0: if(hlm_hist_level_hifrq>0) then

       ! Canopy Resistance

       call this%set_history_var(vname='FATES_STOMATAL_COND',                     &
            units='mol m-2 s-1', long='mean stomatal conductance',                &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_c_stomata_si)

       call this%set_history_var(vname='FATES_LBLAYER_COND', units='mol m-2 s-1', &
            long='mean leaf boundary layer conductance', use_default='active',    &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM',  upfreq=group_hifr_simple,                &
            ivar=ivar, initialize=initialize_variables, index = ih_c_lblayer_si)

       ! Temperature

       call this%set_history_var(vname='FATES_TVEG', units='degree_Celsius', &
            long='fates instantaneous mean vegetation temperature by site', &
            use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_hifr_simple, &
            ivar=ivar, initialize=initialize_variables, index = ih_tveg_si )

       call this%set_history_var(vname='FATES_VIS_RAD_ERROR', units='-',          &
            long='mean two-stream solver error for VIS', use_default='active',            &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_hifr_simple,                 &
            ivar=ivar, initialize=initialize_variables, index = ih_vis_rad_err_si)

       call this%set_history_var(vname='FATES_NIR_RAD_ERROR', units='-',          &
            long='mean two-stream solver error for NIR', use_default='active',            &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', upfreq=group_hifr_simple,                 &
            ivar=ivar, initialize=initialize_variables, index = ih_nir_rad_err_si)

       ! Ecosystem Carbon Fluxes (updated rapidly, upfreq=group_hifr_simple)

       call this%set_history_var(vname='FATES_GPP', units='kg m-2 s-1',           &
            long='gross primary production in kg carbon per m2 per second',       &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables, index = ih_gpp_si)

       call this%set_history_var(vname='FATES_MAINT_RESP', units='kg m-2 s-1',    &
            long='maintenance respiration in kg carbon per m2 land area per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_maint_resp_si)

       call this%set_history_var(vname='FATES_MAINT_RESP_UNREDUCED', units='kg m-2 s-1',    &
            long='diagnostic maintenance respiration if the low-carbon-storage reduction is ignored', &
            use_default='unactive', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_maint_resp_unreduced_si)

       ! fast fluxes separated canopy/understory
       call this%set_history_var(vname='FATES_GPP_CANOPY', units='kg m-2 s-1',    &
            long='gross primary production of canopy plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_gpp_canopy_si)

       call this%set_history_var(vname='FATES_AUTORESP_CANOPY',                  &
            units='kg m-2 s-1',                                                   &
            long='autotrophic respiration of canopy plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_ar_canopy_si)

       call this%set_history_var(vname='FATES_GPP_USTORY',                        &
            units='kg m-2 s-1',                                                   &
            long='gross primary production of understory plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_gpp_understory_si)

       call this%set_history_var(vname='FATES_AUTORESP_USTORY',                   &
            units='kg m-2 s-1',                                                   &
            long='autotrophic respiration of understory plants in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                 &
            index = ih_ar_understory_si)



       call this%set_history_var(vname='FATES_LEAFMAINTAR',                       &
            units = 'kg m-2 s-1',                                                &
            long='leaf maintenance autotrophic respiration in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_leaf_mr_si)

       call this%set_history_var(vname='FATES_FROOTMAINTAR',                      &
            units = 'kg m-2 s-1',                                                &
            long='fine root maintenance autotrophic respiration in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_froot_mr_si)

       call this%set_history_var(vname='FATES_CROOTMAINTAR',                      &
            units = 'kg m-2 s-1',                                                &
            long='live coarse root maintenance autotrophic respiration in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_livecroot_mr_si)

       call this%set_history_var(vname='FATES_LSTEMMAINTAR',                      &
            units = 'kg m-2 s-1',                                                &
            long='live stem maintenance autotrophic respiration in kg carbon per m2 per second', &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_livestem_mr_si)

       call this%set_history_var(vname='FATES_NEP', units='kg m-2 s-1',           &
            long='net ecosystem production in kg carbon per m2 per second',      &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',    &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables,                &
            index = ih_nep_si)

       call this%set_history_var(vname='FATES_HET_RESP', units='kg m-2 s-1',      &
            long='heterotrophic respiration in kg carbon per m2 per second',      &
            use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',     &
            upfreq=group_hifr_simple, ivar=ivar, initialize=initialize_variables, index = ih_hr_si)

       hydro_active_if0: if(hlm_use_planthydro.eq.itrue) then
          call this%set_history_var(vname='FATES_SAPFLOW', units='kg m-2 s-1',    &
               long='areal sap flow rate in kg per m2 ground area per second',               &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM', &
               upfreq=group_hydr_simple, ivar=ivar, initialize=initialize_variables,             &
               index = ih_sapflow_si)
          
          call this%set_history_var(vname='FATES_ROOTWGT_SOILVWC', units='m3 m-3', &
               long='soil volumetric water content, weighted by root area',       &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_hydr_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_rootwgt_soilvwc_si)

          call this%set_history_var(vname='FATES_ROOTWGT_SOILVWCSAT',             &
               units='m3 m-3',                                                    &
               long='soil saturated volumetric water content, weighted by root area', &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_hydr_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_rootwgt_soilvwcsat_si)

          call this%set_history_var(vname='FATES_ROOTWGT_SOILMATPOT', units='Pa', &
               long='soil matric potential, weighted by root area',               &
               use_default='active', avgflag='A', vtype=site_r8, hlms='CLM:ALM',  &
               upfreq=group_hydr_simple, ivar=ivar, initialize=initialize_variables,              &
               index = ih_rootwgt_soilmatpot_si)

          call this%set_history_var(vname='FATES_ROOTUPTAKE', units='kg m-2 s-1', &
               long='root water uptake rate', use_default='active', avgflag='A',  &
               vtype=site_r8, hlms='CLM:ALM', upfreq=group_hydr_simple, ivar=ivar,                &
               initialize=initialize_variables, index = ih_rootuptake_si)


          call this%set_history_var(vname='FATES_VEGH2O', units = 'kg m-2',       &
               long='water stored inside vegetation tissues (leaf, stem, roots)', &
               use_default='inactive', avgflag='A', vtype=site_r8,               &
               hlms='CLM:ALM', upfreq=group_hydr_simple, ivar=ivar,                              &
               initialize=initialize_variables, index = ih_h2oveg_si)

          call this%set_history_var(vname='FATES_VEGH2O_HYDRO_ERR',               &
               units = 'kg m-2',                                                 &
               long='cumulative net borrowed (+) from plant_stored_h2o due to plant hydrodynamics', &
               use_default='inactive', avgflag='A', vtype=site_r8,               &
               hlms='CLM:ALM', upfreq=group_hydr_simple, ivar=ivar,                              &
               initialize=initialize_variables, index = ih_h2oveg_hydro_err_si)
       end if hydro_active_if0

       !HERE

       if_hifrq1: if(hlm_hist_level_hifrq>1) then

          ! This next group are multidimensional variables that are updated
          ! over the short timestep. We turn off these variables when we want
          ! to save time (and some space)

          call this%set_history_var(vname='FATES_GPP_AP', units='kg m-2 s-1',        &
               long='gross primary productivity by age bin in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_gpp_si_age)

          call this%set_history_var(vname='FATES_GPP_LU', units='kg m-2 s-1',        &
               long='gross primary productivity by land use type in kg carbon per m2 per second', &
               use_default='inactive', avgflag='A', vtype=site_landuse_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_gpp_si_landuse)

          call this%set_history_var(vname='FATES_RDARK_USTORY_SZ',               &
               units = 'kg m-2 s-1',                                                &
               long='dark respiration for understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_rdark_understory_si_scls)

          call this%set_history_var(vname='FATES_LSTEMMAINTAR_USTORY_SZ',        &
               units = 'kg m-2 s-1',                                                &
               long='live stem maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_livestem_mr_understory_si_scls)

          call this%set_history_var(vname='FATES_CROOTMAINTAR_USTORY_SZ',        &
               units = 'kg m-2 s-1',                                                &
               long='live coarse root maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_livecroot_mr_understory_si_scls)

          call this%set_history_var(vname='FATES_FROOTMAINTAR_USTORY_SZ',        &
               units = 'kg m-2 s-1',                                                &
               long='fine root maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_froot_mr_understory_si_scls)

          call this%set_history_var(vname='FATES_GROWAR_USTORY_SZ',              &
               units = 'kg m-2 s-1',                                                &
               long='growth autotrophic respiration of understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_resp_g_understory_si_scls)

          call this%set_history_var(vname='FATES_MAINTAR_USTORY_SZ',             &
               units = 'kg m-2 s-1',                                                &
               long='maintenance autotrophic respiration of understory plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM',                                                      &
               upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables,                &
               index = ih_resp_m_understory_si_scls)

          call this%set_history_var(vname='FATES_RDARK_CANOPY_SZ',                   &
               units = 'kg m-2 s-1',                                                &
               long='dark respiration for canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_rdark_canopy_si_scls)

          call this%set_history_var(vname='FATES_CROOTMAINTAR_CANOPY_SZ',            &
               units = 'kg m-2 s-1',                                                &
               long='live coarse root maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_livecroot_mr_canopy_si_scls)

          call this%set_history_var(vname='FATES_FROOTMAINTAR_CANOPY_SZ',            &
               units = 'kg m-2 s-1',                                                &
               long='live coarse root maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_froot_mr_canopy_si_scls)

          call this%set_history_var(vname='FATES_GROWAR_CANOPY_SZ',                  &
               units = 'kg m-2 s-1',                                                 &
               long='growth autotrophic respiration of canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_resp_g_canopy_si_scls)

          call this%set_history_var(vname='FATES_MAINTAR_CANOPY_SZ',                 &
               units = 'kg m-2 s-1',                                                &
               long='maintenance autotrophic respiration of canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_resp_m_canopy_si_scls)

          call this%set_history_var(vname='FATES_LSTEMMAINTAR_CANOPY_SZ',            &
               units = 'kg m-2 s-1',                                                &
               long='live stem maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size', &
               use_default='inactive', avgflag='A', vtype=site_size_r8,             &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables,                                     &
               index = ih_livestem_mr_canopy_si_scls)

          call this%set_history_var(vname='FATES_AUTORESP_SZPF',                     &
               units = 'kg m-2 s-1',                                                &
               long='total autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_si_scpf)

          call this%set_history_var(vname='FATES_GROWAR_SZPF',                       &
               units = 'kg m-2 s-1',                                                &
               long='growth autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_grow_si_scpf)

          call this%set_history_var(vname='FATES_MAINTAR_SZPF',                      &
               units = 'kg m-2 s-1',          &
               long='maintenance autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_maint_si_scpf)

          call this%set_history_var(vname='FATES_RDARK_SZPF',                        &
               units = 'kg m-2 s-1',                                                &
               long='dark portion of maintenance autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_darkm_si_scpf)

          call this%set_history_var(vname='FATES_AGSAPMAINTAR_SZPF',                 &
               units = 'kg m-2 s-1',                                                &
               long='above-ground sapwood maintenance autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_agsapm_si_scpf)

          call this%set_history_var(vname='FATES_BGSAPMAINTAR_SZPF',                 &
               units = 'kg m-2 s-1',                                                &
               long='below-ground sapwood maintenance autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_crootm_si_scpf)

          call this%set_history_var(vname='FATES_FROOTMAINTAR_SZPF',                 &
               units = 'kg m-2 s-1',                                                &
               long='fine root maintenance autotrophic respiration in kg carbon per m2 per second by pft/size', &
               use_default='inactive', avgflag='A', vtype=site_size_pft_r8,         &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_ar_frootm_si_scpf)

          call this%set_history_var(vname='FATES_PARSUN_CLLL', units='W m-2',      &
               long='PAR absorbed in the sun by each canopy and leaf layer',         &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                  &
               initialize=initialize_variables, index = ih_parsun_z_si_cnlf)

          call this%set_history_var(vname='FATES_PARSHA_CLLL', units='W m-2',      &
               long='PAR absorbed in the shade by each canopy and leaf layer',       &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar,                                  &
               initialize=initialize_variables, index = ih_parsha_z_si_cnlf)

          call this%set_history_var(vname='FATES_PARSUN_CLLLPF', units='W m-2',    &
               long='PAR absorbed in the sun by each canopy, leaf, and PFT',         &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parsun_z_si_cnlfpft)

          call this%set_history_var(vname='FATES_PARSHA_CLLLPF', units='W m-2',    &
               long='PAR absorbed in the shade by each canopy, leaf, and PFT',       &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parsha_z_si_cnlfpft)

          call this%set_history_var(vname='FATES_PARSUN_CL', units='W m-2',        &
               long='PAR absorbed by sunlit leaves in each canopy layer', &
               use_default='inactive', avgflag='A', vtype=site_can_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parsun_si_can )

          call this%set_history_var(vname='FATES_PARSHA_CL', units='W m-2',        &
               long='PAR absorbed by shaded leaves in each canopy layer', &
               use_default='inactive', avgflag='A', vtype=site_can_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parsha_si_can)

          call this%set_history_var(vname='FATES_LAISUN_CLLL', units='m2 m-2',     &
               long='LAI in the sun by each canopy and leaf layer',                  &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisun_z_si_cnlf)

          call this%set_history_var(vname='FATES_LAISHA_CLLL', units='m2 m-2',     &
               long='LAI in the shade by each canopy and leaf layer',                &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisha_z_si_cnlf)

          call this%set_history_var(vname='FATES_LAISUN_CLLLPF', units='m2 m-2',   &
               long='Sunlit leaf area by each canopy, leaf, and PFT',                  &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisun_clllpf)

          call this%set_history_var(vname='FATES_LAISHA_CLLLPF', units='m2 m-2',            &
               long='Shaded leaf area by each canopy, leaf, and PFT',            &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisha_clllpf)

          call this%set_history_var(vname='FATES_PARPROF_DIR_CLLLPF', units='W m-2', &
               long='radiative profile of direct PAR through each canopy, leaf, and PFT', &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parprof_dir_si_cnlfpft)

          call this%set_history_var(vname='FATES_PARPROF_DIF_CLLLPF', units='W m-2', &
               long='radiative profile of diffuse PAR through each canopy, leaf, and PFT', &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parprof_dif_si_cnlfpft)

          call this%set_history_var(vname='FATES_LAISUN_CL', units='m2 m-2',     &
               long='LAI of sunlit leaves by canopy layer',     &
               use_default='inactive', avgflag='A', vtype=site_can_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisun_si_can)

          call this%set_history_var(vname='FATES_LAISHA_CL', units='m2 m-2',     &
               long='LAI of shaded leaves by canopy layer',   &
               use_default='inactive', avgflag='A', vtype=site_can_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_laisha_si_can)

          call this%set_history_var(vname='FATES_PARPROF_DIR_CLLL', units='W m-2',   &
               long='radiative profile of direct PAR through each canopy and leaf layer (averaged across PFTs)', &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parprof_dir_si_cnlf)

          call this%set_history_var(vname='FATES_PARPROF_DIF_CLLL', units='W m-2',   &
               long='radiative profile of diffuse PAR through each canopy and leaf layer (averaged across PFTs)', &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_parprof_dif_si_cnlf)

          ! canopy-resolved fluxes and structure

          call this%set_history_var(vname='FATES_NET_C_UPTAKE_CLLL',                 &
               units='kg m-2 s-1',                                                   &
               long='net carbon uptake in kg carbon per m2 per second by each canopy and leaf layer per unit ground area (i.e. divide by CROWNAREA_CLLL to make per leaf area)', &
               use_default='inactive', avgflag='A', vtype=site_cnlf_r8,              &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_ts_net_uptake_si_cnlf)

          call this%set_history_var(vname='FATES_CROWNFRAC_CLLLPF', units='m2 m-2', &
               long='area fraction of the canopy footprint occupied by each canopy-leaf-pft layer', &
               use_default='inactive', avgflag='A', vtype=site_cnlfpft_r8,           &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_crownfrac_clllpf)

          call this%set_history_var(vname='FATES_LBLAYER_COND_AP',                   &
               units='mol m-2 s-1',                                                  &
               long='mean leaf boundary layer conductance - by patch age',           &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM',  upfreq=group_hifr_complx, ivar=ivar,                                 &
               initialize=initialize_variables, index = ih_c_lblayer_si_age)

          ! Canopy resistance
          call this%set_history_var(vname='FATES_STOMATAL_COND_AP',                  &
               units='mol m-2 s-1', long='mean stomatal conductance - by patch age', &
               use_default='inactive', avgflag='A', vtype=site_age_r8,               &
               hlms='CLM:ALM', upfreq=group_hifr_complx, ivar=ivar, initialize=initialize_variables, &
               index = ih_c_stomata_si_age)

          ! PLANT HYDRAULICS

          hydro_active_if1: if(hlm_use_planthydro.eq.itrue) then

             call this%set_history_var(vname='FATES_ERRH2O_SZPF', units='kg s-1',    &
                  long='mean individual water balance error in kg per individual per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_errh2o_scpf)

             call this%set_history_var(vname='FATES_TRAN_SZPF', units='kg s-1',      &
                  long='mean individual transpiration rate in kg per individual per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_tran_scpf)

             call this%set_history_var(vname='FATES_SAPFLOW_SZPF', units='kg m-2 s-1', &
                  long='areal sap flow rate dimensioned by size x pft in kg per m2 ground area per second', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_sapflow_scpf)


             call this%set_history_var(vname='FATES_ITERH1_SZPF', units='count indiv-1 step-1', &
                  long='water balance error iteration diagnostic 1', &
                  use_default='inactive', &
                  avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',     &
                  upfreq=group_hydr_complx, ivar=ivar, initialize=initialize_variables, index = ih_iterh1_scpf )

             call this%set_history_var(vname='FATES_ITERH2_SZPF', units='count indiv-1 step-1', &
                  long='water balance error iteration diagnostic 2', &
                  use_default='inactive', &
                  avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM',     &
                  upfreq=group_hydr_complx, ivar=ivar, initialize=initialize_variables, index = ih_iterh2_scpf )

             call this%set_history_var(vname='FATES_ABSROOT_H2O_SZPF',               &
                  units='m3 m-3',                                                   &
                  long='absorbing volumetric root water content by size class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_ath_scpf)

             call this%set_history_var(vname='FATES_TRANSROOT_H2O_SZPF',             &
                  units='m3 m-3',                                                   &
                  long='transporting volumetric root water content by size class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index =  ih_tth_scpf)

             call this%set_history_var(vname='FATES_STEM_H2O_SZPF', units='m3 m-3',   &
                  long='stem volumetric water content by size class x pft',         &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_sth_scpf)

             call this%set_history_var(vname='FATES_LEAF_H2O_SZPF', units='m3 m-3',   &
                  long='leaf volumetric water content by size class x pft',         &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_lth_scpf)

             call this%set_history_var(vname='FATES_ABSROOT_H2OPOT_SZPF', units='Pa',   &
                  long='absorbing root water potential by size class x pft',        &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_awp_scpf)

             call this%set_history_var(vname='FATES_TRANSROOT_H2OPOT_SZPF',          &
                  units='Pa', long='transporting root water potential by size class x pft', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_twp_scpf)

             call this%set_history_var(vname='FATES_STEM_H2OPOT_SZPF', units='Pa',   &
                  long='stem water potential by size class x pft',                  &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_swp_scpf)

             call this%set_history_var(vname='FATES_LEAF_H2OPOT_SZPF', units='Pa',   &
                  long='leaf water potential by size class x pft',                  &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_lwp_scpf)

             call this%set_history_var(vname='FATES_ABSROOT_CONDFRAC_SZPF', units='1',   &
                  long='absorbing root fraction (0-1) of condutivity by size class x pft', &
                  use_default='active', avgflag='A', vtype=site_size_pft_r8,        &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_aflc_scpf)

             call this%set_history_var(vname='FATES_TRANSROOT_CONDFRAC_SZPF', units='1', &
                  long='transporting root fraction (0-1) of condutivity by size class x pft', &
                  use_default='active', avgflag='A', vtype=site_size_pft_r8,        &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_tflc_scpf)

             call this%set_history_var(vname='FATES_STEM_CONDFRAC_SZPF', units='1',   &
                  long='stem water fraction (0-1) of condutivity by size class x pft', &
                  use_default='active', avgflag='A', vtype=site_size_pft_r8,        &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_sflc_scpf)

             call this%set_history_var(vname='FATES_LEAF_CONDFRAC_SZPF', units='1',   &
                  long='leaf water fraction (0-1) of condutivity by size class x pft', &
                  use_default='active', avgflag='A', vtype=site_size_pft_r8,        &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_lflc_scpf)

             call this%set_history_var(vname='FATES_BTRAN_SZPF', units='1',          &
                  long='mean individual level BTRAN by size class x pft',           &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,      &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                              &
                  initialize=initialize_variables, index = ih_btran_scpf)

             call this%set_history_var(vname='FATES_SOILMATPOT_SL', units='Pa',      &
                  long='soil water matric potenial by soil layer',                   &
                  use_default='inactive', avgflag='A', vtype=site_soil_r8,           &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_soilmatpot_sl)

             call this%set_history_var(vname='FATES_SOILVWC_SL', units='m3 m-3',     &
                  long='soil volumetric water content by soil layer',                &
                  use_default='inactive', avgflag='A', vtype=site_soil_r8,           &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_soilvwc_sl)

             call this%set_history_var(vname='FATES_SOILVWCSAT_SL', units='m3 m-3',  &
                  long='soil saturated volumetric water content by soil layer',      &
                  use_default='inactive', avgflag='A', vtype=site_soil_r8,           &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_soilvwcsat_sl)

             call this%set_history_var(vname='FATES_ROOTUPTAKE_SL',                  &
                  units='kg m-2 s-1',                                               &
                  long='root water uptake rate by soil layer',                       &
                  use_default='inactive', avgflag='A', vtype=site_soil_r8,           &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_rootuptake_sl)

             call this%set_history_var(vname='FATES_ROOTUPTAKE0_SZPF',               &
                  units='kg m-2 m-1 s-1',                                            &
                  long='root water uptake from 0 to to 10 cm depth, by plant size x pft ', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_rootuptake0_scpf)

             call this%set_history_var(vname='FATES_ROOTUPTAKE10_SZPF',              &
                  units='kg m-2 m-1 s-1',                                            &
                  long='root water uptake from 10 to to 50 cm depth, by plant size x pft ', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_rootuptake10_scpf)

             call this%set_history_var(vname='FATES_ROOTUPTAKE50_SZPF',              &
                  units='kg m-2 m-1 s-1',                                            &
                  long='root water uptake from 50 to to 100 cm depth, by plant size x pft ', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_rootuptake50_scpf)

             call this%set_history_var(vname='FATES_ROOTUPTAKE100_SZPF',             &
                  units='kg m-2 m-1 s-1',                                            &
                  long='root water uptake below 100 cm depth, by plant size x pft ', &
                  use_default='inactive', avgflag='A', vtype=site_size_pft_r8,       &
                  hlms='CLM:ALM', upfreq=group_hydr_complx, ivar=ivar,                               &
                  initialize=initialize_variables, index = ih_rootuptake100_scpf)

          end if hydro_active_if1

       end if if_hifrq1


    end if if_hifrq0

    ! Must be last thing before return
    this%num_history_vars_ = ivar

  end subroutine define_history_vars


  ! ====================================================================================
  ! DEPRECATED, TRANSITIONAL OR FUTURE CODE SECTION
  ! ====================================================================================

  !subroutine set_fates_hio_str(tag,iotype_name, iostr_val)

  !       ! Arguments
  !       character(len=*), intent(in)           :: tag
  !       character(len=*), optional, intent(in) :: iotype_name
  !       integer, optional, intent(in)         :: iostr_val

  !       ! local variables
  !       logical              :: all_set
  !       integer,  parameter  :: unset_int = -999
  !       real(r8), parameter  :: unset_double = -999.9
  !       integer              :: ityp, idim

  !       select case (trim(tag))
  !       case('flush_to_unset')
  !          write(*, *) ''
  !          write(*, *) 'Flushing FATES IO types prior to transfer from host'
  !          do ityp=1,ubound(iovar_str, 1)
  !             iovar_str(ityp)%dimsize = unset_int
  !             iovar_str(ityp)%active  = .false.
  !          end do

  !       case('check_allset')
  !          do ityp=1,ubound(iovar_str, 1)
  !             write(*, *) 'Checking to see if ',iovar_str(ityp)%name, ' IO communicators were sent to FATES'
  !             if(iovar_str(ityp)%active)then
  !                if(iovar_str(ityp)%offset .eq. unset_int) then
  !                   write(*, *) 'FATES offset information of IO type:', iovar_str(ityp)%name
  !                   write(*, *) 'was never set'
  !                   ! end_run('MESSAGE')
  !                end if
  !                do idim=1, iovar_str(ityp)%ndims
  !                   if(iovar_str(ityp)%dimsize(idim) .eq. unset_int) then
  !                      write(*, *) 'FATES dimension information of IO type:', iovar_str(ityp)%name
  !                      write(*, *) 'was never set'
  !                      ! end_run('MESSAGE')
  !                   end if
  !                end do
  !             end if
  !          end do
  !          write(*, *) 'Checked. All history IO specifications properly sent to FATES.'
  !       case default

  !          ! Must have two arguments if this is not a check or flush
  !          if(present(iostr_val) .and. present(iotype_name))then
  !
  !             ! Tag in this case is dimsize or offset
  !             select case (trim(tag))
  !
  !             case('offset')
  !                ityp=iotype_index(trim(iotype_name))
  !                iovar_str(ityp)%offset = iostr_val
  !                write(*, *) 'Transfering offset for IOTYPE',iotype_name, ' to FATES'

  !             case('dimsize1')
  !                ityp=iotype_index(trim(iotype_name))
  !                iovar_str(ityp)%dimsize(1) = iostr_val
  !                write(*, *) 'Transfering 1st dimension size for IOTYPE',iotype_name, ' to FATES'

  !             case('dimsize2')
  !                ityp=iotype_index(trim(iotype_name))
  !                if(ubound(iovar_str(ityp)%dimsize, 1)==1)then
  !                   write(fates_log(), *) 'Transfering second dimensional bound to unallocated space'
  !                   write(fates_log(), *) 'type:', iotype_name
  !                   ! end_run
  !                end if
  !                iovar_str(ityp)%dimsize(2) = iostr_val
  !                write(*, *) 'Transfering 2nd dimension size for IOTYPE',iotype_name, ' to FATES'

  !             case('dimsize3')
  !                ityp=iotype_index(trim(iotype_name))
  !                if(ubound(iovar_str(ityp)%dimsize, 1)<3)then
  !                   write(fates_log(), *) 'Transfering third dimensional bound to unallocated space'
  !                   write(fates_log(), *) 'type:', iotype_name
  !                   ! end_run
  !                end if
  !                iovar_str(ityp)%dimsize(3) = iostr_val
  !                write(*, *) 'Transfering 3rd dimension size for IOTYPE',iotype_name, ' to FATES'

  !             case default
  !                write(*, *) 'IO parameter not recognized:', trim(tag)
  !                ! end_run
  !             end select
  !          else
  !             write(*, *) 'no value was provided for the tag'
  !          end if
  !
  !       end select
  !       return
  !     end subroutine set_fates_hio_str



end module FatesHistoryInterfaceMod
