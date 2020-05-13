module FatesHistoryInterfaceMod

  use FatesConstantsMod        , only : r8 => fates_r8
  use FatesConstantsMod        , only : fates_avg_flag_length
  use FatesConstantsMod        , only : fates_short_string_length
  use FatesConstantsMod        , only : fates_long_string_length
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : calloc_abs_error
  use FatesConstantsMod        , only : mg_per_kg
  use FatesConstantsMod        , only : pi_const
  use FatesGlobals             , only : fates_log
  use FatesGlobals             , only : endrun => fates_endrun
  use EDTypesMod               , only : nclmax
  use EDTypesMod               , only : ican_upper
  use EDTypesMod               , only : element_pos
  use EDTypesMod               , only : num_elements
  use EDTypesMod               , only : site_fluxdiags_type
  use EDtypesMod               , only : ed_site_type
  use EDtypesMod               , only : ed_cohort_type
  use EDtypesMod               , only : ed_patch_type  
  use EDtypesMod               , only : AREA
  use EDtypesMod               , only : AREA_INV
  use EDTypesMod               , only : numWaterMem
  use EDTypesMod               , only : num_vegtemp_mem
  use EDTypesMod               , only : site_massbal_type
  use EDTypesMod               , only : element_list
  use FatesIODimensionsMod     , only : fates_io_dimension_type
  use FatesIOVariableKindMod   , only : fates_io_variable_kind_type
  use FatesHistoryVariableType , only : fates_history_variable_type
  use FatesInterfaceTypesMod        , only : hlm_hio_ignore_val
  use FatesInterfaceTypesMod        , only : hlm_use_planthydro
  use FatesInterfaceTypesMod        , only : hlm_use_ed_st3
  use FatesInterfaceTypesMod        , only : hlm_use_cohort_age_tracking
  use FatesInterfaceTypesMod        , only : numpft
  use FatesInterfaceTypesMod        , only : hlm_freq_day
  use EDParamsMod              , only : ED_val_comp_excln
  use EDParamsMod              , only : ED_val_phen_coldtemp
  use FatesInterfaceTypesMod        , only : nlevsclass, nlevage
  use FatesInterfaceTypesMod        , only : nlevheight
  use FatesInterfaceTypesMod        , only : bc_in_type
  use FatesInterfaceTypesMod        , only : hlm_model_day
  use FatesInterfaceTypesMod        , only : nlevcoage

  ! FIXME(bja, 2016-10) need to remove CLM dependancy 
  use EDPftvarcon              , only : EDPftvarcon_inst

  ! CIME Globals
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use shr_infnan_mod           , only : isnan => shr_infnan_isnan
  use FatesConstantsMod        , only : g_per_kg
  use FatesConstantsMod        , only : ha_per_m2
  use FatesConstantsMod        , only : days_per_sec
  use FatesConstantsMod        , only : sec_per_day
  use FatesConstantsMod        , only : days_per_year
  use FatesConstantsMod        , only : years_per_day
  use FatesLitterMod           , only : litter_type
  use FatesConstantsMod        , only : secondaryforest

  use PRTGenericMod            , only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod            , only : struct_organ, store_organ, repro_organ
  use PRTGenericMod            , only : all_carbon_elements
  use PRTGenericMod            , only : carbon12_element


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
  
  ! Indices to 1D Patch variables

  integer :: ih_trimming_si
  integer :: ih_area_plant_si
  integer :: ih_area_trees_si

  integer :: ih_cwd_elcwd

  integer :: ih_litter_in_si    ! carbon only
  integer :: ih_litter_out_si   ! carbon only
  integer :: ih_seed_bank_si    ! carbon only
  integer :: ih_seeds_in_si     ! carbon only

  integer :: ih_litter_in_elem
  integer :: ih_litter_out_elem
  integer :: ih_seed_bank_elem
  integer :: ih_seeds_in_local_elem
  integer :: ih_seeds_in_extern_elem
  integer :: ih_seed_decay_elem
  integer :: ih_seed_germ_elem

  integer :: ih_fines_ag_elem
  integer :: ih_fines_bg_elem
  integer :: ih_cwd_ag_elem
  integer :: ih_cwd_bg_elem
  integer :: ih_burn_flux_elem

  integer :: ih_daily_temp
  integer :: ih_daily_rh
  integer :: ih_daily_prec
 
  integer :: ih_bstore_si
  integer :: ih_bdead_si
  integer :: ih_balive_si
  integer :: ih_bleaf_si
  integer :: ih_bsapwood_si
  integer :: ih_bfineroot_si
  integer :: ih_btotal_si
  integer :: ih_agb_si
  integer :: ih_npp_si
  integer :: ih_gpp_si
  integer :: ih_aresp_si
  integer :: ih_maint_resp_si
  integer :: ih_growth_resp_si
  integer :: ih_ar_canopy_si
  integer :: ih_gpp_canopy_si
  integer :: ih_ar_understory_si
  integer :: ih_gpp_understory_si
  integer :: ih_canopy_biomass_si
  integer :: ih_understory_biomass_si

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

  integer :: ih_nep_si

  integer :: ih_c_stomata_si
  integer :: ih_c_lblayer_si

  integer :: ih_fire_c_to_atm_si


  integer :: ih_cbal_err_fates_si
  integer :: ih_err_fates_si

  integer :: ih_npatches_si
  integer :: ih_ncohorts_si
  integer :: ih_demotion_carbonflux_si
  integer :: ih_promotion_carbonflux_si
  integer :: ih_canopy_mortality_carbonflux_si
  integer :: ih_understory_mortality_carbonflux_si
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
  integer :: ih_fraction_secondary_forest_si
  integer :: ih_biomass_secondary_forest_si
  integer :: ih_woodproduct_si
  integer :: ih_h2oveg_si
  integer :: ih_h2oveg_dead_si
  integer :: ih_h2oveg_recruit_si
  integer :: ih_h2oveg_growturn_err_si
  integer :: ih_h2oveg_pheno_err_si
  integer :: ih_h2oveg_hydro_err_si


  
  integer :: ih_site_cstatus_si
  integer :: ih_site_dstatus_si
  integer :: ih_gdd_si
  integer :: ih_site_nchilldays_si
  integer :: ih_site_ncolddays_si
  integer :: ih_cleafoff_si
  integer :: ih_cleafon_si
  integer :: ih_dleafoff_si
  integer :: ih_dleafon_si
  integer :: ih_meanliqvol_si

  integer :: ih_nesterov_fire_danger_si
  integer :: ih_fire_intensity_area_product_si
  integer :: ih_spitfire_ros_si
  integer :: ih_fire_ros_area_product_si
  integer :: ih_effect_wspeed_si
  integer :: ih_tfc_ros_si
  integer :: ih_tfc_ros_area_product_si
  integer :: ih_fire_intensity_si
  integer :: ih_fire_area_si
  integer :: ih_fire_fuel_bulkd_si
  integer :: ih_fire_fuel_eff_moist_si
  integer :: ih_fire_fuel_sav_si
  integer :: ih_fire_fuel_mef_si
  integer :: ih_sum_fuel_si

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
  
  integer :: ih_bstor_canopy_si_scpf
  integer :: ih_bstor_understory_si_scpf
  integer :: ih_bleaf_canopy_si_scpf
  integer :: ih_bleaf_understory_si_scpf
  integer :: ih_mortality_canopy_si_scpf
  integer :: ih_mortality_understory_si_scpf
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
  integer :: ih_crownfiremort_si_scpf
  integer :: ih_cambialfiremort_si_scpf

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
  integer :: ih_mortality_si_pft
  integer :: ih_crownarea_si_pft


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
  integer :: ih_secondaryforest_area_si_age
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
  integer :: ih_fabd_sun_si_cnlf
  integer :: ih_fabd_sha_si_cnlf
  integer :: ih_fabi_sun_si_cnlf
  integer :: ih_fabi_sha_si_cnlf
  integer :: ih_ts_net_uptake_si_cnlf
  integer :: ih_crownarea_si_cnlf
  integer :: ih_parprof_dir_si_cnlf
  integer :: ih_parprof_dif_si_cnlf

  ! indices to (site x [canopy layer x leaf layer x pft]) variables
  integer :: ih_parsun_z_si_cnlfpft
  integer :: ih_parsha_z_si_cnlfpft
  integer :: ih_laisun_z_si_cnlfpft
  integer :: ih_laisha_z_si_cnlfpft
  integer :: ih_fabd_sun_si_cnlfpft
  integer :: ih_fabd_sha_si_cnlfpft
  integer :: ih_fabi_sun_si_cnlfpft
  integer :: ih_fabi_sha_si_cnlfpft
  integer :: ih_parprof_dir_si_cnlfpft
  integer :: ih_parprof_dif_si_cnlfpft

  ! indices to (site x canopy layer) variables
  integer :: ih_parsun_top_si_can
  integer :: ih_parsha_top_si_can
  integer :: ih_laisun_top_si_can
  integer :: ih_laisha_top_si_can
  integer :: ih_fabd_sun_top_si_can
  integer :: ih_fabd_sha_top_si_can
  integer :: ih_fabi_sun_top_si_can
  integer :: ih_fabi_sha_top_si_can
  integer :: ih_crownarea_si_can

  ! The number of variable dim/kind types we have defined (static)

  integer, parameter, public :: fates_history_num_dimensions = 50
  integer, parameter, public :: fates_history_num_dim_kinds = 50

  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type, public :: iovar_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: patch1_index(:) ! maps site index to the HIO patch 1st position
  end type iovar_map_type


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
     
     type(iovar_map_type), pointer :: iovar_map(:)

   
     !! THESE WERE EXPLICITLY PRIVATE WHEN TYPE WAS PUBLIC
     integer, private :: patch_index_, column_index_, levgrnd_index_, levscpf_index_
     integer, private :: levscls_index_, levpft_index_, levage_index_
     integer, private :: levfuel_index_, levcwdsc_index_, levscag_index_
     integer, private :: levcan_index_, levcnlf_index_, levcnlfpft_index_
     integer, private :: levscagpft_index_, levagepft_index_
     integer, private :: levheight_index_
     integer, private :: levelem_index_, levelpft_index_
     integer, private :: levelcwd_index_, levelage_index_
     integer, private :: levcacls_index_, levcapf_index_

     
   contains
     
     procedure :: Init
     procedure :: SetThreadBoundsEach
     procedure :: initialize_history_vars
     procedure :: assemble_history_output_types
     
     procedure :: update_history_dyn
     procedure :: update_history_prod
     procedure :: update_history_cbal
     procedure :: update_history_hydraulics

     ! 'get' methods used by external callers to access private read only data

     procedure :: num_history_vars
     procedure :: patch_index
     procedure :: column_index
     procedure :: levgrnd_index
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
     procedure :: levscag_index
     procedure :: levscagpft_index
     procedure :: levagepft_index
     procedure :: levheight_index
     procedure :: levelem_index
     procedure :: levelpft_index
     procedure :: levelcwd_index
     procedure :: levelage_index

     ! private work functions
     procedure, private :: define_history_vars
     procedure, private :: set_history_var
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: flush_hvars

     procedure, private :: set_patch_index
     procedure, private :: set_column_index
     procedure, private :: set_levgrnd_index
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
     procedure, private :: set_levscag_index
     procedure, private :: set_levscagpft_index
     procedure, private :: set_levagepft_index
     procedure, private :: set_levheight_index
     
     procedure, private :: set_levelem_index
     procedure, private :: set_levelpft_index
     procedure, private :: set_levelcwd_index
     procedure, private :: set_levelage_index


  end type fates_history_interface_type
   
  character(len=*), parameter :: sourcefile = &
         __FILE__

contains

  ! ======================================================================
  
  subroutine Init(this, num_threads, fates_bounds)

    use FatesIODimensionsMod, only : patch, column, levgrnd, levscpf
    use FatesIODimensionsMod, only : levscls, levpft, levage
    use FatesIODimensionsMod, only : levcacls, levcapf
    use FatesIODimensionsMod, only : levfuel, levcwdsc, levscag
    use FatesIODimensionsMod, only : levscagpft, levagepft
    use FatesIODimensionsMod, only : levcan, levcnlf, levcnlfpft
    use FatesIODimensionsMod, only : fates_bounds_type
    use FatesIODimensionsMod, only : levheight
    use FatesIODimensionsMod, only : levelem, levelpft
    use FatesIODimensionsMod, only : levelcwd, levelage

    implicit none

    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_patch_index(dim_count)
    call this%dim_bounds(dim_count)%Init(patch, num_threads, &
         fates_bounds%patch_begin, fates_bounds%patch_end)

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    dim_count = dim_count + 1
    call this%set_levgrnd_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levgrnd, num_threads, &
         fates_bounds%ground_begin, fates_bounds%ground_end)

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
    

    ! FIXME(bja, 2016-10) assert(dim_count == FatesHistorydimensionmod::num_dimension_types)

    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%iovar_map(num_threads))
    
  end subroutine Init

  ! ======================================================================
  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)

    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index
    
    index = this%patch_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%patch_begin, thread_bounds%patch_end)

    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)

    index = this%levgrnd_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%ground_begin, thread_bounds%ground_end)

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
    

    


    
  end subroutine SetThreadBoundsEach
  
  ! ===================================================================================
  subroutine assemble_history_output_types(this)

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8

   implicit none

    class(fates_history_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(patch_r8, 1, this%patch_index())

    call this%set_dim_indices(site_r8, 1, this%column_index())

    call this%set_dim_indices(patch_ground_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(site_ground_r8, 1, this%column_index())
    call this%set_dim_indices(site_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(patch_size_pft_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_size_pft_r8, 2, this%levscpf_index())

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
       stop
       !end_run
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
 subroutine set_patch_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%patch_index_ = index
 end subroutine set_patch_index

 integer function patch_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   patch_index = this%patch_index_
 end function patch_index

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
 subroutine set_levgrnd_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levgrnd_index_ = index
 end subroutine set_levgrnd_index

 integer function levgrnd_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levgrnd_index = this%levgrnd_index_
 end function levgrnd_index

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

 subroutine flush_hvars(this,nc,upfreq_in)
 
   class(fates_history_interface_type)        :: this
   integer,intent(in)                     :: nc
   integer,intent(in)                     :: upfreq_in
   integer                      :: ivar
   integer                      :: lb1,ub1,lb2,ub2

   do ivar=1,ubound(this%hvars,1)
      if (this%hvars(ivar)%upfreq == upfreq_in) then ! Only flush variables with update on dynamics step
         call this%hvars(ivar)%flush(nc, this%dim_bounds, this%dim_kinds)
         
      end if
   end do
   
end subroutine flush_hvars

  
  ! =====================================================================================
   
  subroutine set_history_var(this, vname, units, long, use_default, avgflag, vtype, &
       hlms, flushval, upfreq, ivar, initialize, index)

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
    real(r8), intent(in)          :: flushval ! IF THE TYPE IS AN INT WE WILL round with NINT
    integer, intent(in)           :: upfreq
    logical, intent(in) :: initialize
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used
    

    ! locals
    integer :: ub1, lb1, ub2, lb2    ! Bounds for allocating the var
    integer :: ityp

    logical :: write_var

    write_var = check_hlm_list(trim(hlms), trim(hlm_name))
    if( write_var ) then
       ivar  = ivar+1
       index = ivar    
       
       if (initialize) then
          call this%hvars(ivar)%Init(vname, units, long, use_default, &
               vtype, avgflag, flushval, upfreq, &
               fates_history_num_dim_kinds, this%dim_kinds, this%dim_bounds)
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
    ! PA_R8   : 1D patch scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8

    implicit none
    
    ! Arguments
    class(fates_history_interface_type), intent(inout) :: this
       

    integer :: index

    ! 1d Patch
    index = 1
    call this%dim_kinds(index)%Init(patch_r8, 1)

    ! 1d Site
    index = index + 1
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! patch x ground
    index = index + 1
    call this%dim_kinds(index)%Init(patch_ground_r8, 2)

    ! patch x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(patch_size_pft_r8, 2)

    ! site x ground
    index = index + 1
    call this%dim_kinds(index)%Init(site_ground_r8, 2)

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


    ! FIXME(bja, 2016-10) assert(index == fates_history_num_dim_kinds)
  end subroutine init_dim_kinds_maps

 ! =======================================================================

  subroutine update_history_cbal(this,nc,nsites,sites,bc_in,dtime)

     use EDtypesMod          , only : ed_site_type
      

     ! Arguments
     class(fates_history_interface_type)             :: this
     integer                 , intent(in)            :: nc   ! clump index
     integer                 , intent(in)            :: nsites
     type(ed_site_type)      , intent(inout), target :: sites(nsites)
     type(bc_in_type)        , intent(in)            :: bc_in(nsites)
     real(r8)                , intent(in)            :: dtime   ! Time-step (s)
     
     ! Locals
     integer  :: s        ! The local site index
     integer  :: io_si     ! The site index of the IO array
     real(r8) :: inv_dtime  ! inverse of dtime (faster math)
     type(ed_cohort_type), pointer  :: ccohort ! current cohort
     type(ed_patch_type) , pointer  :: cpatch ! current patch
     
     associate( hio_nep_si => this%hvars(ih_nep_si)%r81d )
       
       ! ---------------------------------------------------------------------------------
       ! Flush arrays to values defined by %flushval (see registry entry in
       ! subroutine define_history_vars()
       ! ---------------------------------------------------------------------------------

       call this%flush_hvars(nc,upfreq_in=3)        

       inv_dtime = 1._r8/dtime
       
       do s = 1,nsites
           
           io_si  = this%iovar_map(nc)%site_index(s)

           hio_nep_si(io_si) = -bc_in(s)%tot_het_resp ! (gC/m2/s)
           
           cpatch => sites(s)%oldest_patch
           do while(associated(cpatch))
               ccohort => cpatch%shortest
               do while(associated(ccohort))
               
                   ! Add up the total Net Ecosystem Production
                   ! for this timestep.  [gC/m2/s]
                   hio_nep_si(io_si) = hio_nep_si(io_si) + &
                        (ccohort%gpp_tstep - ccohort%resp_tstep) * &
                        g_per_kg * ccohort%n * area_inv * inv_dtime
                   ccohort => ccohort%taller
               end do
               cpatch => cpatch%younger
           end do
       end do
      end associate

   end subroutine update_history_cbal
   

  ! ====================================================================================
  
  subroutine update_history_dyn(this,nc,nsites,sites)
    
    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after Ecosystem Dynamics have been processed.
    ! ---------------------------------------------------------------------------------
    

    use EDtypesMod          , only : nfsc
    use FatesLitterMod      , only : ncwd
    use EDtypesMod          , only : ican_upper
    use EDtypesMod          , only : ican_ustory
    use FatesSizeAgeTypeIndicesMod, only : get_sizeage_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_sizeagepft_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_agepft_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_height_index
    use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
    use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
    use EDTypesMod        , only : nlevleaf
    use EDParamsMod,           only : ED_val_history_height_bin_edges

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    
    ! Locals
    type(litter_type), pointer         :: litt_c   ! Pointer to the carbon12 litter pool
    type(litter_type), pointer         :: litt     ! Generic pointer to any litter pool
    type(site_fluxdiags_type), pointer :: flux_diags
    type(site_fluxdiags_type), pointer :: flux_diags_c
    type(site_massbal_type), pointer :: site_mass

    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa, ipa2 ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    integer  :: cwd
    integer  :: elcwd, elpft            ! combined index of element and pft or cwd
    integer  :: i_scpf,i_pft,i_scls     ! iterators for scpf, pft, and scls dims
    integer  :: i_cacls, i_capf      ! iterators for cohort age and cohort age x pft
    integer  :: i_cwd,i_fuel            ! iterators for cwd and fuel dims
    integer  :: iscag        ! size-class x age index
    integer  :: iscagpft     ! size-class x age x pft index
    integer  :: iagepft      ! age x pft index
    integer  :: ican, ileaf, cnlf_indx  ! iterators for leaf and canopy level
    integer  :: height_bin_max, height_bin_min   ! which height bin a given cohort's canopy is in
    integer  :: i_heightbin  ! iterator for height bins
    integer  :: el           ! Loop index for elements
    integer  :: model_day_int ! integer model day from reference 
    integer  :: ageclass_since_anthrodist  ! what is the equivalent age class for
                                           ! time-since-anthropogenic-disturbance of secondary forest

    
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: dbh         ! diameter ("at breast height")
    real(r8) :: coage       ! cohort age 
    real(r8) :: npp_partition_error ! a check that the NPP partitions sum to carbon allocation
    real(r8) :: frac_canopy_in_bin  ! fraction of a leaf's canopy that is within a given height bin
    real(r8) :: binbottom,bintop    ! edges of height bins
    
    real(r8) :: gpp_cached ! variable used to cache gpp value in previous time step; for C13 discrimination

    ! The following are all carbon states, turnover and net allocation flux variables
    ! the organs of relevance should be self explanatory
    real(r8) :: sapw_c
    real(r8) :: struct_c
    real(r8) :: leaf_c
    real(r8) :: fnrt_c
    real(r8) :: store_c
    real(r8) :: alive_c
    real(r8) :: total_c
    real(r8) :: sapw_c_turnover
    real(r8) :: store_c_turnover
    real(r8) :: leaf_c_turnover
    real(r8) :: fnrt_c_turnover
    real(r8) :: struct_c_turnover
    real(r8) :: sapw_c_net_alloc
    real(r8) :: store_c_net_alloc
    real(r8) :: leaf_c_net_alloc
    real(r8) :: fnrt_c_net_alloc
    real(r8) :: struct_c_net_alloc
    real(r8) :: repro_c_net_alloc
    real(r8) :: area_frac

    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    real(r8), parameter :: reallytalltrees = 1000.   ! some large number (m)
    
    integer :: tmp

    associate( hio_npatches_si         => this%hvars(ih_npatches_si)%r81d, &
               hio_ncohorts_si         => this%hvars(ih_ncohorts_si)%r81d, &
               hio_trimming_si         => this%hvars(ih_trimming_si)%r81d, &
               hio_area_plant_si       => this%hvars(ih_area_plant_si)%r81d, &
               hio_area_trees_si  => this%hvars(ih_area_trees_si)%r81d, & 
               hio_canopy_spread_si    => this%hvars(ih_canopy_spread_si)%r81d, &
               hio_biomass_si_pft      => this%hvars(ih_biomass_si_pft)%r82d, &
               hio_leafbiomass_si_pft  => this%hvars(ih_leafbiomass_si_pft)%r82d, &
               hio_storebiomass_si_pft => this%hvars(ih_storebiomass_si_pft)%r82d, &
               hio_nindivs_si_pft      => this%hvars(ih_nindivs_si_pft)%r82d, &
               hio_recruitment_si_pft  => this%hvars(ih_recruitment_si_pft)%r82d, &
               hio_mortality_si_pft    => this%hvars(ih_mortality_si_pft)%r82d, &
               hio_crownarea_si_pft    => this%hvars(ih_crownarea_si_pft)%r82d, &
               hio_nesterov_fire_danger_si => this%hvars(ih_nesterov_fire_danger_si)%r81d, &
               hio_spitfire_ros_si     => this%hvars(ih_spitfire_ros_si)%r81d, &
               hio_fire_ros_area_product_si=> this%hvars(ih_fire_ros_area_product_si)%r81d, &
               hio_tfc_ros_si          => this%hvars(ih_tfc_ros_si)%r81d, &
               hio_tfc_ros_area_product_si => this%hvars(ih_tfc_ros_area_product_si)%r81d, &
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
               hio_seed_bank_si        => this%hvars(ih_seed_bank_si)%r81d, &
               hio_seeds_in_si         => this%hvars(ih_seeds_in_si)%r81d, &
               hio_litter_in_elem      => this%hvars(ih_litter_in_elem)%r82d, &
               hio_litter_out_elem     => this%hvars(ih_litter_out_elem)%r82d, &
               hio_seed_bank_elem      => this%hvars(ih_seed_bank_elem)%r82d, &
               hio_seeds_in_local_elem => this%hvars(ih_seeds_in_local_elem)%r82d, &
               hio_seed_in_extern_elem => this%hvars(ih_seeds_in_extern_elem)%r82d, & 
               hio_seed_decay_elem     => this%hvars(ih_seed_decay_elem)%r82d, &
               hio_seed_germ_elem      => this%hvars(ih_seed_germ_elem)%r82d, &

               hio_bstore_si           => this%hvars(ih_bstore_si)%r81d, &
               hio_bdead_si            => this%hvars(ih_bdead_si)%r81d, &
               hio_balive_si           => this%hvars(ih_balive_si)%r81d, &
               hio_bleaf_si            => this%hvars(ih_bleaf_si)%r81d, &
               hio_bsapwood_si         => this%hvars(ih_bsapwood_si)%r81d, &
               hio_bfineroot_si        => this%hvars(ih_bfineroot_si)%r81d, &
               hio_btotal_si           => this%hvars(ih_btotal_si)%r81d, &
               hio_agb_si              => this%hvars(ih_agb_si)%r81d, &
               hio_canopy_biomass_si   => this%hvars(ih_canopy_biomass_si)%r81d, &
               hio_understory_biomass_si   => this%hvars(ih_understory_biomass_si)%r81d, &
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
               hio_npp_leaf_si         => this%hvars(ih_npp_leaf_si)%r81d, &
               hio_npp_seed_si         => this%hvars(ih_npp_seed_si)%r81d, &
               hio_npp_stem_si         => this%hvars(ih_npp_stem_si)%r81d, &
               hio_npp_froot_si        => this%hvars(ih_npp_froot_si)%r81d, &
               hio_npp_croot_si        => this%hvars(ih_npp_croot_si)%r81d, &
               hio_npp_stor_si         => this%hvars(ih_npp_stor_si)%r81d, &
               hio_bstor_canopy_si_scpf      => this%hvars(ih_bstor_canopy_si_scpf)%r82d, &
               hio_bstor_understory_si_scpf  => this%hvars(ih_bstor_understory_si_scpf)%r82d, &
               hio_bleaf_canopy_si_scpf      => this%hvars(ih_bleaf_canopy_si_scpf)%r82d, &
               hio_bleaf_understory_si_scpf  => this%hvars(ih_bleaf_understory_si_scpf)%r82d, &
               hio_mortality_canopy_si_scpf         => this%hvars(ih_mortality_canopy_si_scpf)%r82d, &
               hio_mortality_understory_si_scpf     => this%hvars(ih_mortality_understory_si_scpf)%r82d, &
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

               hio_fire_c_to_atm_si  => this%hvars(ih_fire_c_to_atm_si)%r81d, &
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
               hio_m10_si_cacls        => this%hvars(ih_m10_si_cacls)%r82d, &
              
	       hio_c13disc_si_scpf     => this%hvars(ih_c13disc_si_scpf)%r82d, &                    

               hio_cwd_elcwd           => this%hvars(ih_cwd_elcwd)%r82d, &
               hio_cwd_ag_elem         => this%hvars(ih_cwd_ag_elem)%r82d, &
               hio_cwd_bg_elem         => this%hvars(ih_cwd_bg_elem)%r82d, &
               hio_fines_ag_elem       => this%hvars(ih_fines_bg_elem)%r82d, &
               hio_fines_bg_elem       => this%hvars(ih_fines_ag_elem)%r82d, &
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
               hio_demotion_carbonflux_si        => this%hvars(ih_demotion_carbonflux_si)%r81d, &
               hio_promotion_rate_si_scls        => this%hvars(ih_promotion_rate_si_scls)%r82d, &
               hio_trimming_canopy_si_scls         => this%hvars(ih_trimming_canopy_si_scls)%r82d, &
               hio_trimming_understory_si_scls     => this%hvars(ih_trimming_understory_si_scls)%r82d, &
               hio_crown_area_canopy_si_scls         => this%hvars(ih_crown_area_canopy_si_scls)%r82d, &
               hio_crown_area_understory_si_scls     => this%hvars(ih_crown_area_understory_si_scls)%r82d, &
               hio_promotion_carbonflux_si       => this%hvars(ih_promotion_carbonflux_si)%r81d, &
               hio_canopy_mortality_carbonflux_si     => this%hvars(ih_canopy_mortality_carbonflux_si)%r81d, &
               hio_understory_mortality_carbonflux_si => this%hvars(ih_understory_mortality_carbonflux_si)%r81d, &
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
               hio_fraction_secondary_forest_si   => this%hvars(ih_fraction_secondary_forest_si)%r81d, &
               hio_biomass_secondary_forest_si    => this%hvars(ih_biomass_secondary_forest_si)%r81d, &
               hio_woodproduct_si                 => this%hvars(ih_woodproduct_si)%r81d, &
               hio_agesince_anthrodist_si_age     => this%hvars(ih_agesince_anthrodist_si_age)%r82d, &
               hio_secondaryforest_area_si_age    => this%hvars(ih_secondaryforest_area_si_age)%r82d, &
               hio_area_burnt_si_age              => this%hvars(ih_area_burnt_si_age)%r82d, &
               ! hio_fire_rate_of_spread_front_si_age  => this%hvars(ih_fire_rate_of_spread_front_si_age)%r82d, &
               hio_fire_intensity_si_age          => this%hvars(ih_fire_intensity_si_age)%r82d, &
               hio_fire_sum_fuel_si_age           => this%hvars(ih_fire_sum_fuel_si_age)%r82d, &
               hio_burnt_frac_litter_si_fuel      => this%hvars(ih_burnt_frac_litter_si_fuel)%r82d, &
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
               hio_crownarea_si_can                 => this%hvars(ih_crownarea_si_can)%r82d, &
               hio_nplant_si_scag                   => this%hvars(ih_nplant_si_scag)%r82d, &
               hio_nplant_canopy_si_scag            => this%hvars(ih_nplant_canopy_si_scag)%r82d, &
               hio_nplant_understory_si_scag        => this%hvars(ih_nplant_understory_si_scag)%r82d, &
               hio_ddbh_canopy_si_scag              => this%hvars(ih_ddbh_canopy_si_scag)%r82d, &
               hio_ddbh_understory_si_scag          => this%hvars(ih_ddbh_understory_si_scag)%r82d, &
               hio_mortality_canopy_si_scag         => this%hvars(ih_mortality_canopy_si_scag)%r82d, &
               hio_mortality_understory_si_scag     => this%hvars(ih_mortality_understory_si_scag)%r82d, &
               hio_site_cstatus_si                  => this%hvars(ih_site_cstatus_si)%r81d, &
               hio_site_dstatus_si                  => this%hvars(ih_site_dstatus_si)%r81d, &
               hio_gdd_si                           => this%hvars(ih_gdd_si)%r81d, &
               hio_site_ncolddays_si                => this%hvars(ih_site_ncolddays_si)%r81d, &
               hio_site_nchilldays_si               => this%hvars(ih_site_nchilldays_si)%r81d, &
               hio_cleafoff_si                      => this%hvars(ih_cleafoff_si)%r81d, &
               hio_cleafon_si                       => this%hvars(ih_cleafon_si)%r81d, &
               hio_dleafoff_si                      => this%hvars(ih_dleafoff_si)%r81d, &
               hio_dleafon_si                       => this%hvars(ih_dleafoff_si)%r81d, &
               hio_meanliqvol_si                    => this%hvars(ih_meanliqvol_si)%r81d, &
               hio_cbal_err_fates_si                => this%hvars(ih_cbal_err_fates_si)%r81d, &
               hio_err_fates_si                     => this%hvars(ih_err_fates_si)%r82d )

               
      ! ---------------------------------------------------------------------------------
      ! Flush arrays to values defined by %flushval (see registry entry in
      ! subroutine define_history_vars()
      ! ---------------------------------------------------------------------------------
      call this%flush_hvars(nc,upfreq_in=1)


      ! If we don't have dynamics turned on, we just abort these diagnostics
      if (hlm_use_ed_st3.eq.itrue) return

      model_day_int = nint(hlm_model_day)

      ! ---------------------------------------------------------------------------------
      ! Loop through the FATES scale hierarchy and fill the history IO arrays
      ! ---------------------------------------------------------------------------------
      
      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1

         ! Total carbon model error [kgC/day -> mgC/day]
         hio_cbal_err_fates_si(io_si) = &
               sites(s)%mass_balance(element_pos(carbon12_element))%err_fates * mg_per_kg

         ! Total carbon lost to atmosphere from burning (kgC/site/day -> gC/m2/s)
         hio_fire_c_to_atm_si(io_si) = &
              sites(s)%mass_balance(element_pos(carbon12_element))%burn_flux_to_atm * &
              g_per_kg * ha_per_m2 * days_per_sec

         ! Total model error [kg/day -> mg/day]  (all elements)
         do el = 1, num_elements
             site_mass => sites(s)%mass_balance(el)
             hio_err_fates_si(io_si,el) = site_mass%err_fates * mg_per_kg

             ! Total element lost to atmosphere from burning (kg/site/day -> g/m2/s)
             hio_burn_flux_elem(io_si,el) = &
                  sites(s)%mass_balance(el)%burn_flux_to_atm * &
                  g_per_kg * ha_per_m2 * days_per_sec

         end do

         hio_canopy_spread_si(io_si)        = sites(s)%spread

         ! Update the site statuses (stati?)
         hio_site_cstatus_si(io_si)   = real(sites(s)%cstatus,r8)
         hio_site_dstatus_si(io_si)   = real(sites(s)%dstatus,r8)

         !count number of days for leaves off
         hio_site_nchilldays_si(io_si) = real(sites(s)%nchilldays,r8)
         hio_site_ncolddays_si(io_si)  = real(sites(s)%ncolddays,r8)

            
         hio_gdd_si(io_si)      = sites(s)%grow_deg_days
         hio_cleafoff_si(io_si) = real(model_day_int - sites(s)%cleafoffdate,r8)
         hio_cleafon_si(io_si)  = real(model_day_int - sites(s)%cleafondate,r8)
         hio_dleafoff_si(io_si) = real(model_day_int - sites(s)%dleafoffdate,r8)
         hio_dleafon_si(io_si)  = real(model_day_int - sites(s)%dleafondate,r8)

         if(model_day_int>numWaterMem)then
            hio_meanliqvol_si(io_si) = &
                 sum(sites(s)%water_memory(1:numWaterMem))/real(numWaterMem,r8)
         end if

         ! track total wood product accumulation at the site level
         hio_woodproduct_si(io_si)          = sites(s)%resources_management%trunk_product_site &
              * AREA_INV * g_per_kg
         
         ! site-level fire variables
         hio_nesterov_fire_danger_si(io_si) = sites(s)%acc_NI

         ! If hydraulics are turned on, track the error terms
         ! associated with dynamics

         if(hlm_use_planthydro.eq.itrue)then
            this%hvars(ih_h2oveg_dead_si)%r81d(io_si)         = sites(s)%si_hydr%h2oveg_dead
            this%hvars(ih_h2oveg_recruit_si)%r81d(io_si)      = sites(s)%si_hydr%h2oveg_recruit
            this%hvars(ih_h2oveg_growturn_err_si)%r81d(io_si) = sites(s)%si_hydr%h2oveg_growturn_err
            this%hvars(ih_h2oveg_pheno_err_si)%r81d(io_si)    = sites(s)%si_hydr%h2oveg_pheno_err
         end if

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            ! Increment the number of patches per site
            hio_npatches_si(io_si) = hio_npatches_si(io_si) + 1._r8

            cpatch%age_class  = get_age_class_index(cpatch%age)

            ! Increment the fractional area in each age class bin
            hio_area_si_age(io_si,cpatch%age_class) = hio_area_si_age(io_si,cpatch%age_class) &
                 + cpatch%area * AREA_INV

            ! Increment some patch-age-resolved diagnostics
            hio_lai_si_age(io_si,cpatch%age_class) = hio_lai_si_age(io_si,cpatch%age_class) &
                  + sum(cpatch%tlai_profile(:,:,:)) * cpatch%area

            hio_ncl_si_age(io_si,cpatch%age_class) = hio_ncl_si_age(io_si,cpatch%age_class) &
                  + cpatch%ncl_p * cpatch%area
            hio_npatches_si_age(io_si,cpatch%age_class) = hio_npatches_si_age(io_si,cpatch%age_class) + 1._r8
            if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
               hio_zstar_si_age(io_si,cpatch%age_class) = hio_zstar_si_age(io_si,cpatch%age_class) &
                    + cpatch%zstar * cpatch%area * AREA_INV
            endif

            ! some diagnostics on secondary forest area and its age distribution
            if ( cpatch%anthro_disturbance_label .eq. secondaryforest ) then
               hio_fraction_secondary_forest_si(io_si) = hio_fraction_secondary_forest_si(io_si) + &
                    cpatch%area * AREA_INV
               
               ageclass_since_anthrodist = get_age_class_index(cpatch%age_since_anthro_disturbance)
               
               hio_agesince_anthrodist_si_age(io_si,ageclass_since_anthrodist) = &
                    hio_agesince_anthrodist_si_age(io_si,ageclass_since_anthrodist)  &
                    + cpatch%area * AREA_INV

               hio_secondaryforest_area_si_age(io_si,cpatch%age_class) = &
                    hio_secondaryforest_area_si_age(io_si,cpatch%age_class)  &
                    + cpatch%area * AREA_INV
            endif
            
            !!! patch-age-resolved fire variables
            do i_pft = 1,numpft
               ! for scorch height, weight the value by patch area within any given age calss (in the event that there is
               ! more than one patch per age class.
               iagepft = cpatch%age_class + (i_pft-1) * nlevage
               hio_scorch_height_si_agepft(io_si,iagepft) = hio_scorch_height_si_agepft(io_si,iagepft) + &
                    cpatch%Scorch_ht(i_pft) * cpatch%area
            end do

            hio_area_burnt_si_age(io_si,cpatch%age_class) = hio_area_burnt_si_age(io_si,cpatch%age_class) + &
                 cpatch%frac_burnt * cpatch%area * AREA_INV

            ! hio_fire_rate_of_spread_front_si_age(io_si, cpatch%age_class) = hio_fire_rate_of_spread_si_age(io_si, cpatch%age_class) + &
            !      cpatch%ros_front * cpatch*frac_burnt * cpatch%area * AREA_INV

            hio_fire_intensity_si_age(io_si, cpatch%age_class) = hio_fire_intensity_si_age(io_si, cpatch%age_class) + &
                 cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV

            hio_fire_sum_fuel_si_age(io_si, cpatch%age_class) = hio_fire_sum_fuel_si_age(io_si, cpatch%age_class) + &
                 cpatch%sum_fuel * cpatch%area * AREA_INV
             
            if(associated(cpatch%tallest))then
               hio_trimming_si(io_si) = hio_trimming_si(io_si) + cpatch%tallest%canopy_trim * cpatch%area * AREA_INV
            endif
            
            hio_area_plant_si(io_si) = hio_area_plant_si(io_si) + min(cpatch%total_canopy_area,cpatch%area) * AREA_INV

            hio_area_trees_si(io_si) = hio_area_trees_si(io_si) + min(cpatch%total_tree_area,cpatch%area) * AREA_INV
            
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ft = ccohort%pft

               call sizetype_class_index(ccohort%dbh, ccohort%pft, ccohort%size_class, ccohort%size_by_pft_class)
               call coagetype_class_index(ccohort%coage, ccohort%pft, &
                                          ccohort%coage_class, ccohort%coage_by_pft_class)
              
               ! Increment the number of cohorts per site
               hio_ncohorts_si(io_si) = hio_ncohorts_si(io_si) + 1._r8
               
               n_perm2   = ccohort%n * AREA_INV
                              
               hio_canopy_area_si_age(io_si,cpatch%age_class) = hio_canopy_area_si_age(io_si,cpatch%age_class) &
                    + ccohort%c_area * AREA_INV

               ! calculate leaf height distribution, assuming leaf area is evenly distributed thru crown depth
               height_bin_max = get_height_index(ccohort%hite)
               height_bin_min = get_height_index(ccohort%hite * (1._r8 - EDPftvarcon_inst%crown(ft)))
               do i_heightbin = height_bin_min, height_bin_max
                  binbottom = ED_val_history_height_bin_edges(i_heightbin)
                  if (i_heightbin .eq. nlevheight) then
                     bintop = reallytalltrees
                  else
                     bintop = ED_val_history_height_bin_edges(i_heightbin+1)
                  endif
                  ! what fraction of a cohort's crown is in this height bin?
                  frac_canopy_in_bin = (min(bintop,ccohort%hite) - &
                       max(binbottom,ccohort%hite * (1._r8 - EDPftvarcon_inst%crown(ft)))) / &
                       (ccohort%hite * EDPftvarcon_inst%crown(ft))
                  !
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

               ! Update biomass components


               ! Mass pools [kgC]
               sapw_c   = ccohort%prt%GetState(sapw_organ, all_carbon_elements)
               struct_c = ccohort%prt%GetState(struct_organ, all_carbon_elements)
               leaf_c   = ccohort%prt%GetState(leaf_organ, all_carbon_elements)
               fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)
               store_c  = ccohort%prt%GetState(store_organ, all_carbon_elements)
 
               alive_c  = leaf_c + fnrt_c + sapw_c
               total_c  = alive_c + store_c + struct_c

               hio_bleaf_si(io_si)     = hio_bleaf_si(io_si)  + n_perm2 * leaf_c  * g_per_kg
               hio_bstore_si(io_si)    = hio_bstore_si(io_si) + n_perm2 * store_c * g_per_kg
               hio_bdead_si(io_si)     = hio_bdead_si(io_si)  + n_perm2 * struct_c * g_per_kg
               hio_balive_si(io_si)    = hio_balive_si(io_si) + n_perm2 * alive_c * g_per_kg

               hio_bsapwood_si(io_si)  = hio_bsapwood_si(io_si)   + n_perm2 * sapw_c  * g_per_kg
               hio_bfineroot_si(io_si) = hio_bfineroot_si(io_si) + n_perm2 * fnrt_c * g_per_kg
               hio_btotal_si(io_si)    = hio_btotal_si(io_si) + n_perm2 * total_c * g_per_kg

               hio_agb_si(io_si)       = hio_agb_si(io_si) + n_perm2 * g_per_kg * &
                    ( leaf_c + (sapw_c + struct_c + store_c) * EDPftvarcon_inst%allom_agb_frac(ccohort%pft) )

               ! Update PFT partitioned biomass components
               hio_leafbiomass_si_pft(io_si,ft) = hio_leafbiomass_si_pft(io_si,ft) + &
                    (ccohort%n * AREA_INV) * leaf_c     * g_per_kg
             
               hio_storebiomass_si_pft(io_si,ft) = hio_storebiomass_si_pft(io_si,ft) + &
                    (ccohort%n * AREA_INV) * store_c   * g_per_kg
               
               hio_nindivs_si_pft(io_si,ft) = hio_nindivs_si_pft(io_si,ft) + &
                    ccohort%n * AREA_INV

               hio_biomass_si_pft(io_si, ft) = hio_biomass_si_pft(io_si, ft) + &
                    (ccohort%n * AREA_INV) * total_c * g_per_kg

               ! Update PFT crown area
               hio_crownarea_si_pft(io_si, ft) = hio_crownarea_si_pft(io_si, ft) + &
                    ccohort%c_area 

               ! update total biomass per age bin
               hio_biomass_si_age(io_si,cpatch%age_class) = hio_biomass_si_age(io_si,cpatch%age_class) &
                    + total_c * ccohort%n * AREA_INV

               ! track the total biomass on all secondary lands
               if ( cpatch%anthro_disturbance_label .eq. secondaryforest ) then
                  hio_biomass_secondary_forest_si(io_si) = hio_biomass_secondary_forest_si(io_si) + &
                       total_c * ccohort%n * AREA_INV
               endif

               ! Site by Size-Class x PFT (SCPF) 
               ! ------------------------------------------------------------------------

               dbh = ccohort%dbh !-0.5*(1./365.25)*ccohort%ddbhdt

               ! Flux Variables (cohorts must had experienced a day before any of these values
               ! have any meaning, otherwise they are just inialization values
               if( .not.(ccohort%isnew) ) then

                  ! Turnover pools [kgC/day] * [day/yr] = [kgC/yr]
                  sapw_c_turnover   = ccohort%prt%GetTurnover(sapw_organ, all_carbon_elements) * days_per_year
                  store_c_turnover  = ccohort%prt%GetTurnover(store_organ, all_carbon_elements) * days_per_year
                  leaf_c_turnover   = ccohort%prt%GetTurnover(leaf_organ, all_carbon_elements) * days_per_year
                  fnrt_c_turnover   = ccohort%prt%GetTurnover(fnrt_organ, all_carbon_elements) * days_per_year
                  struct_c_turnover = ccohort%prt%GetTurnover(struct_organ, all_carbon_elements) * days_per_year
                  
                  ! Net change from allocation and transport [kgC/day] * [day/yr] = [kgC/yr]
                 sapw_c_net_alloc   = ccohort%prt%GetNetAlloc(sapw_organ, all_carbon_elements) * days_per_year
                  store_c_net_alloc  = ccohort%prt%GetNetAlloc(store_organ, all_carbon_elements) * days_per_year
                  leaf_c_net_alloc   = ccohort%prt%GetNetAlloc(leaf_organ, all_carbon_elements) * days_per_year
                  fnrt_c_net_alloc   = ccohort%prt%GetNetAlloc(fnrt_organ, all_carbon_elements) * days_per_year
                  struct_c_net_alloc = ccohort%prt%GetNetAlloc(struct_organ, all_carbon_elements) * days_per_year
                  repro_c_net_alloc  = ccohort%prt%GetNetAlloc(repro_organ, all_carbon_elements) * days_per_year

                  ! ecosystem-level, organ-partitioned NPP/allocation fluxes                                                                                                         
                  hio_npp_leaf_si(io_si) = hio_npp_leaf_si(io_si) + leaf_c_net_alloc * n_perm2
                  hio_npp_seed_si(io_si) = hio_npp_seed_si(io_si) + repro_c_net_alloc * n_perm2
                  hio_npp_stem_si(io_si) = hio_npp_stem_si(io_si) + (sapw_c_net_alloc + struct_c_net_alloc) * n_perm2 * &
                       (EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                  hio_npp_froot_si(io_si) = hio_npp_froot_si(io_si) + fnrt_c_net_alloc * n_perm2
                  hio_npp_croot_si(io_si) = hio_npp_croot_si(io_si) + (sapw_c_net_alloc + struct_c_net_alloc) * n_perm2 * &
                       (1._r8-EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                  hio_npp_stor_si(io_si) = hio_npp_stor_si(io_si) + store_c_net_alloc * n_perm2
                  
                  associate( scpf => ccohort%size_by_pft_class, &

                       scls => ccohort%size_class, &
                       cacls => ccohort%coage_class, &
                       capf => ccohort%coage_by_pft_class)
     
			     
		    gpp_cached = hio_gpp_si_scpf(io_si,scpf)
      
                    hio_gpp_si_scpf(io_si,scpf)      = hio_gpp_si_scpf(io_si,scpf)      + &
                                                       n_perm2*ccohort%gpp_acc_hold  ! [kgC/m2/yr]
                    hio_npp_totl_si_scpf(io_si,scpf) = hio_npp_totl_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_acc_hold *n_perm2
                    
                    
                    hio_npp_leaf_si_scpf(io_si,scpf) = hio_npp_leaf_si_scpf(io_si,scpf) + &
                                                       leaf_c_net_alloc*n_perm2
                    hio_npp_fnrt_si_scpf(io_si,scpf) = hio_npp_fnrt_si_scpf(io_si,scpf) + &
                                                       fnrt_c_net_alloc*n_perm2
                    hio_npp_bgsw_si_scpf(io_si,scpf) = hio_npp_bgsw_si_scpf(io_si,scpf) + &
                                                       sapw_c_net_alloc*n_perm2*           &
                                                       (1._r8-EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                    hio_npp_agsw_si_scpf(io_si,scpf) = hio_npp_agsw_si_scpf(io_si,scpf) + &
                                                       sapw_c_net_alloc*n_perm2*           &
                                                       EDPftvarcon_inst%allom_agb_frac(ccohort%pft)
                    hio_npp_bgdw_si_scpf(io_si,scpf) = hio_npp_bgdw_si_scpf(io_si,scpf) + &
                                                       struct_c_net_alloc*n_perm2*         &
                                                       (1._r8-EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                    hio_npp_agdw_si_scpf(io_si,scpf) = hio_npp_agdw_si_scpf(io_si,scpf) + &
                                                       struct_c_net_alloc*n_perm2*         &
                                                       EDPftvarcon_inst%allom_agb_frac(ccohort%pft)
                    hio_npp_seed_si_scpf(io_si,scpf) = hio_npp_seed_si_scpf(io_si,scpf) + &
                                                       repro_c_net_alloc*n_perm2
                    hio_npp_stor_si_scpf(io_si,scpf) = hio_npp_stor_si_scpf(io_si,scpf) + &
                                                       store_c_net_alloc*n_perm2

                    ! Woody State Variables (basal area growth increment)
                    if (EDPftvarcon_inst%woody(ft) == 1) then

                       ! basal area  [m2/ha]
                       hio_ba_si_scpf(io_si,scpf) = hio_ba_si_scpf(io_si,scpf) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n

                       ! also by size class only
                       hio_ba_si_scls(io_si,scls) = hio_ba_si_scls(io_si,scls) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n

                       ! growth increment
                       hio_ddbh_si_scpf(io_si,scpf) = hio_ddbh_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n

                    end if

                    hio_m1_si_scpf(io_si,scpf) = hio_m1_si_scpf(io_si,scpf) + ccohort%bmort*ccohort%n
                    hio_m2_si_scpf(io_si,scpf) = hio_m2_si_scpf(io_si,scpf) + ccohort%hmort*ccohort%n
                    hio_m3_si_scpf(io_si,scpf) = hio_m3_si_scpf(io_si,scpf) + ccohort%cmort*ccohort%n
                    hio_m7_si_scpf(io_si,scpf) = hio_m7_si_scpf(io_si,scpf) + &
                         (ccohort%lmort_direct+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n
                    hio_m8_si_scpf(io_si,scpf) = hio_m8_si_scpf(io_si,scpf) + ccohort%frmort*ccohort%n
                    hio_m9_si_scpf(io_si,scpf) = hio_m9_si_scpf(io_si,scpf) + ccohort%smort*ccohort%n
                    
                    if (hlm_use_cohort_age_tracking .eq.itrue) then
                       hio_m10_si_scpf(io_si,scpf) = hio_m10_si_scpf(io_si,scpf) + ccohort%asmort*ccohort%n
                       hio_m10_si_capf(io_si,capf) = hio_m10_si_capf(io_si,capf) + ccohort%asmort*ccohort%n
                       hio_m10_si_scls(io_si,scls) = hio_m10_si_scls(io_si,scls) + ccohort%asmort*ccohort%n
                       hio_m10_si_cacls(io_si,cacls) = hio_m10_si_cacls(io_si,cacls)+ &
                            ccohort%asmort*ccohort%n
                    end if
                    
                    hio_m1_si_scls(io_si,scls) = hio_m1_si_scls(io_si,scls) + ccohort%bmort*ccohort%n
                    hio_m2_si_scls(io_si,scls) = hio_m2_si_scls(io_si,scls) + ccohort%hmort*ccohort%n
                    hio_m3_si_scls(io_si,scls) = hio_m3_si_scls(io_si,scls) + ccohort%cmort*ccohort%n
                    hio_m7_si_scls(io_si,scls) = hio_m7_si_scls(io_si,scls) + &
                         (ccohort%lmort_direct+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n
                    hio_m8_si_scls(io_si,scls) = hio_m8_si_scls(io_si,scls) + &
                         ccohort%frmort*ccohort%n
                    hio_m9_si_scls(io_si,scls) = hio_m9_si_scls(io_si,scls) + ccohort%smort*ccohort%n
                  
                  
                   
                    !C13 discrimination
                    if(gpp_cached + ccohort%gpp_acc_hold > 0.0_r8)then
                       hio_c13disc_si_scpf(io_si,scpf) = ((hio_c13disc_si_scpf(io_si,scpf) * gpp_cached) + &
                            (ccohort%c13disc_acc * ccohort%gpp_acc_hold)) / (gpp_cached + ccohort%gpp_acc_hold)
                    else
                       hio_c13disc_si_scpf(io_si,scpf) = 0.0_r8
                    endif

                    ! number density [/ha]
                    hio_nplant_si_scpf(io_si,scpf) = hio_nplant_si_scpf(io_si,scpf) + ccohort%n

                    ! number density along the cohort age dimension
                    if (hlm_use_cohort_age_tracking .eq.itrue) then
                       hio_nplant_si_capf(io_si,capf) = hio_nplant_si_capf(io_si,capf) + ccohort%n
                       hio_nplant_si_cacls(io_si,cacls) = hio_nplant_si_cacls(io_si,cacls)+ccohort%n
                    end if
                    
                    ! number density by size and biomass
                    hio_agb_si_scls(io_si,scls) = hio_agb_si_scls(io_si,scls) + &
                         total_c * ccohort%n * EDPftvarcon_inst%allom_agb_frac(ccohort%pft) * AREA_INV

                    hio_agb_si_scpf(io_si,scpf) = hio_agb_si_scpf(io_si,scpf) + &
                         total_c * ccohort%n * EDPftvarcon_inst%allom_agb_frac(ccohort%pft) * AREA_INV

                    hio_biomass_si_scls(io_si,scls) = hio_biomass_si_scls(io_si,scls) + &
                          total_c * ccohort%n * AREA_INV

                    ! update size-class x patch-age related quantities

                    iscag = get_sizeage_class_index(ccohort%dbh,cpatch%age)
                    
                    hio_nplant_si_scag(io_si,iscag) = hio_nplant_si_scag(io_si,iscag) + ccohort%n

                    hio_nplant_si_scls(io_si,scls) = hio_nplant_si_scls(io_si,scls) + ccohort%n
                    
                  
                    ! update size, age, and PFT - indexed quantities

                    iscagpft = get_sizeagepft_class_index(ccohort%dbh,cpatch%age,ccohort%pft)
                    
                    hio_nplant_si_scagpft(io_si,iscagpft) = hio_nplant_si_scagpft(io_si,iscagpft) + ccohort%n

                    ! update age and PFT - indexed quantities

                    iagepft = get_agepft_class_index(cpatch%age,ccohort%pft)
                    
                    hio_npp_si_agepft(io_si,iagepft) = hio_npp_si_agepft(io_si,iagepft) + &
                         ccohort%n * ccohort%npp_acc_hold * AREA_INV

                    hio_biomass_si_agepft(io_si,iagepft) = hio_biomass_si_agepft(io_si,iagepft) + &
                          total_c * ccohort%n * AREA_INV

                    ! update SCPF/SCLS- and canopy/subcanopy- partitioned quantities
                    if (ccohort%canopy_layer .eq. 1) then
                       hio_nplant_canopy_si_scag(io_si,iscag) = hio_nplant_canopy_si_scag(io_si,iscag) + ccohort%n
                       hio_mortality_canopy_si_scag(io_si,iscag) = hio_mortality_canopy_si_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n
                       hio_ddbh_canopy_si_scag(io_si,iscag) = hio_ddbh_canopy_si_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_canopy_si_scpf(io_si,scpf) = hio_bstor_canopy_si_scpf(io_si,scpf) + &
                             store_c * ccohort%n
                       hio_bleaf_canopy_si_scpf(io_si,scpf) = hio_bleaf_canopy_si_scpf(io_si,scpf) + &
                             leaf_c * ccohort%n

                       hio_canopy_biomass_si(io_si) = hio_canopy_biomass_si(io_si) + n_perm2 * total_c * g_per_kg

                       !hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &
                       !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                       ! ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                       hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &

                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort + & 
                            ccohort%smort + ccohort%asmort) * ccohort%n + &
			    (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                            ccohort%n * sec_per_day * days_per_year

                       hio_nplant_canopy_si_scpf(io_si,scpf) = hio_nplant_canopy_si_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_canopy_si_scls(io_si,scls) = hio_nplant_canopy_si_scls(io_si,scls) + ccohort%n
                       hio_lai_canopy_si_scls(io_si,scls) = hio_lai_canopy_si_scls(io_si,scls) + &
                                                            ccohort%treelai*ccohort%c_area * AREA_INV
                       hio_sai_canopy_si_scls(io_si,scls) = hio_sai_canopy_si_scls(io_si,scls) + &
                                                            ccohort%treesai*ccohort%c_area * AREA_INV
                       hio_trimming_canopy_si_scls(io_si,scls) = hio_trimming_canopy_si_scls(io_si,scls) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_canopy_si_scls(io_si,scls) = hio_crown_area_canopy_si_scls(io_si,scls) + &
                            ccohort%c_area
                       hio_gpp_canopy_si_scpf(io_si,scpf)      = hio_gpp_canopy_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_canopy_si_scpf(io_si,scpf)      = hio_ar_canopy_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold
                       ! growth increment
                       hio_ddbh_canopy_si_scpf(io_si,scpf) = hio_ddbh_canopy_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_canopy_si_scls(io_si,scls) = hio_ddbh_canopy_si_scls(io_si,scls) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_canopy_si_scls(io_si,scls) = hio_mortality_canopy_si_scls(io_si,scls) + &

                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                             ccohort%n * sec_per_day * days_per_year

                       hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * &
                            total_c * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                            (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_c * &
                            ccohort%n * g_per_kg * ha_per_m2
                       

                       hio_carbon_balance_canopy_si_scls(io_si,scls) = hio_carbon_balance_canopy_si_scls(io_si,scls) + &
                             ccohort%n * ccohort%npp_acc_hold
                       
                      
                       hio_leaf_md_canopy_si_scls(io_si,scls) = hio_leaf_md_canopy_si_scls(io_si,scls) + &
                             leaf_c_turnover * ccohort%n
                       hio_root_md_canopy_si_scls(io_si,scls) = hio_root_md_canopy_si_scls(io_si,scls) + &
                             fnrt_c_turnover * ccohort%n
                       hio_bsw_md_canopy_si_scls(io_si,scls) = hio_bsw_md_canopy_si_scls(io_si,scls) + &
                             sapw_c_turnover * ccohort%n
                       hio_bstore_md_canopy_si_scls(io_si,scls) = hio_bstore_md_canopy_si_scls(io_si,scls) + &
                             store_c_turnover * ccohort%n
                       hio_bdead_md_canopy_si_scls(io_si,scls) = hio_bdead_md_canopy_si_scls(io_si,scls) + &
                             struct_c_turnover * ccohort%n
                       hio_seed_prod_canopy_si_scls(io_si,scls) = hio_seed_prod_canopy_si_scls(io_si,scls) + &
                             ccohort%seed_prod * ccohort%n

                       hio_npp_leaf_canopy_si_scls(io_si,scls) = hio_npp_leaf_canopy_si_scls(io_si,scls) + &
                             leaf_c_net_alloc * ccohort%n
                       hio_npp_fnrt_canopy_si_scls(io_si,scls) = hio_npp_fnrt_canopy_si_scls(io_si,scls) + &
                             fnrt_c_net_alloc * ccohort%n
                       hio_npp_sapw_canopy_si_scls(io_si,scls) = hio_npp_sapw_canopy_si_scls(io_si,scls) + &
                             sapw_c_net_alloc * ccohort%n
                       hio_npp_dead_canopy_si_scls(io_si,scls) = hio_npp_dead_canopy_si_scls(io_si,scls) + &
                             struct_c_net_alloc * ccohort%n
                       hio_npp_seed_canopy_si_scls(io_si,scls) = hio_npp_seed_canopy_si_scls(io_si,scls) + &
                             repro_c_net_alloc * ccohort%n
                       hio_npp_stor_canopy_si_scls(io_si,scls) = hio_npp_stor_canopy_si_scls(io_si,scls) + &
                             store_c_net_alloc * ccohort%n
                       
                       hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) = &
                            hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    else
                       hio_nplant_understory_si_scag(io_si,iscag) = hio_nplant_understory_si_scag(io_si,iscag) + ccohort%n
                       hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n
                       hio_ddbh_understory_si_scag(io_si,iscag) = hio_ddbh_understory_si_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_understory_si_scpf(io_si,scpf) = hio_bstor_understory_si_scpf(io_si,scpf) + &
                             store_c * ccohort%n
                       hio_bleaf_understory_si_scpf(io_si,scpf) = hio_bleaf_understory_si_scpf(io_si,scpf) + &
                             leaf_c  * ccohort%n
                       hio_understory_biomass_si(io_si) = hio_understory_biomass_si(io_si) + &
                             n_perm2 * total_c * g_per_kg
                       !hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                        !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + 
                       !      ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                       hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
			    (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                            ccohort%n * sec_per_day * days_per_year

                       hio_nplant_understory_si_scpf(io_si,scpf) = hio_nplant_understory_si_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_understory_si_scls(io_si,scls) = hio_nplant_understory_si_scls(io_si,scls) + ccohort%n
                       hio_lai_understory_si_scls(io_si,scls) = hio_lai_understory_si_scls(io_si,scls) + &
                                                                ccohort%treelai*ccohort%c_area  * AREA_INV
                       hio_sai_understory_si_scls(io_si,scls) = hio_sai_understory_si_scls(io_si,scls) + &
                                                                ccohort%treelai*ccohort%c_area  * AREA_INV
                       hio_trimming_understory_si_scls(io_si,scls) = hio_trimming_understory_si_scls(io_si,scls) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_understory_si_scls(io_si,scls) = hio_crown_area_understory_si_scls(io_si,scls) + &
                            ccohort%c_area
                       hio_gpp_understory_si_scpf(io_si,scpf)      = hio_gpp_understory_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_understory_si_scpf(io_si,scpf)      = hio_ar_understory_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold

                       ! growth increment
                       hio_ddbh_understory_si_scpf(io_si,scpf) = hio_ddbh_understory_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_understory_si_scls(io_si,scls) = hio_ddbh_understory_si_scls(io_si,scls) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_understory_si_scls(io_si,scls) = hio_mortality_understory_si_scls(io_si,scls) + &

                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                             ccohort%n * sec_per_day * days_per_year
                       
                       hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * &
                             total_c * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_c * &
                             ccohort%n * g_per_kg * ha_per_m2

                       hio_carbon_balance_understory_si_scls(io_si,scls) = hio_carbon_balance_understory_si_scls(io_si,scls) + &
                             ccohort%npp_acc_hold * ccohort%n

                       hio_leaf_md_understory_si_scls(io_si,scls) = hio_leaf_md_understory_si_scls(io_si,scls) + &
                            leaf_c_turnover * ccohort%n
                       hio_root_md_understory_si_scls(io_si,scls) = hio_root_md_understory_si_scls(io_si,scls) + &
                            fnrt_c_turnover * ccohort%n
                       hio_bsw_md_understory_si_scls(io_si,scls) = hio_bsw_md_understory_si_scls(io_si,scls) + &
                             sapw_c_turnover * ccohort%n
                       hio_bstore_md_understory_si_scls(io_si,scls) = hio_bstore_md_understory_si_scls(io_si,scls) + &
                             store_c_turnover * ccohort%n
                       hio_bdead_md_understory_si_scls(io_si,scls) = hio_bdead_md_understory_si_scls(io_si,scls) + &
                             struct_c_turnover * ccohort%n
                       hio_seed_prod_understory_si_scls(io_si,scls) = hio_seed_prod_understory_si_scls(io_si,scls) + &
                            ccohort%seed_prod * ccohort%n

                       hio_npp_leaf_understory_si_scls(io_si,scls) = hio_npp_leaf_understory_si_scls(io_si,scls) + &
                             leaf_c_net_alloc * ccohort%n
                       hio_npp_fnrt_understory_si_scls(io_si,scls) = hio_npp_fnrt_understory_si_scls(io_si,scls) + &
                             fnrt_c_net_alloc * ccohort%n
                       hio_npp_sapw_understory_si_scls(io_si,scls) = hio_npp_sapw_understory_si_scls(io_si,scls) + &
                             sapw_c_net_alloc * ccohort%n
                       hio_npp_dead_understory_si_scls(io_si,scls) = hio_npp_dead_understory_si_scls(io_si,scls) + &
                             struct_c_net_alloc * ccohort%n
                       hio_npp_seed_understory_si_scls(io_si,scls) = hio_npp_seed_understory_si_scls(io_si,scls) + &
                             repro_c_net_alloc * ccohort%n
                       hio_npp_stor_understory_si_scls(io_si,scls) = hio_npp_stor_understory_si_scls(io_si,scls) + &
                             store_c_net_alloc * ccohort%n
                       
                       hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) = &
                            hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    endif
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
                               ccohort%n * days_per_year
                       end do
                    end if
                    ccohort%size_class_lasttimestep = scls


                    !
                  end associate
               else  ! i.e. cohort%isnew
                  !
                  ! if cohort is new, track its growth flux into the first size bin
                  i_scpf = (ccohort%pft-1)*nlevsclass+1
                  hio_growthflux_si_scpf(io_si,i_scpf) = hio_growthflux_si_scpf(io_si,i_scpf) + ccohort%n * days_per_year
                  ccohort%size_class_lasttimestep = 1
                                   
               end if

               ! resolve some canopy area profiles, both total and of occupied leaves
               ican = ccohort%canopy_layer
               !
               hio_crownarea_si_can(io_si, ican) = hio_crownarea_si_can(io_si, ican) + ccohort%c_area / AREA
               !
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_crownarea_si_cnlf(io_si, cnlf_indx) = hio_crownarea_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%c_area / AREA
               end do
               
               ccohort => ccohort%taller
            enddo ! cohort loop
            
            ! Patch specific variables that are already calculated
            ! These things are all duplicated. Should they all be converted to LL or array structures RF? 
            ! define scalar to counteract the patch albedo scaling logic for conserved quantities
                        
            ! Update Fire Variables
            hio_spitfire_ros_si(io_si)         = hio_spitfire_ros_si(io_si) + cpatch%ROS_front * cpatch%area * AREA_INV
            hio_fire_ros_area_product_si(io_si)= hio_fire_ros_area_product_si(io_si) + &
                 cpatch%frac_burnt * cpatch%ROS_front * cpatch%area * AREA_INV
            hio_effect_wspeed_si(io_si)        = hio_effect_wspeed_si(io_si) + cpatch%effect_wspeed * cpatch%area * AREA_INV
            hio_tfc_ros_si(io_si)              = hio_tfc_ros_si(io_si) + cpatch%TFC_ROS * cpatch%area * AREA_INV
            hio_tfc_ros_area_product_si(io_si) = hio_tfc_ros_area_product_si(io_si) + &
                 cpatch%frac_burnt * cpatch%TFC_ROS * cpatch%area * AREA_INV
            hio_fire_intensity_si(io_si)       = hio_fire_intensity_si(io_si) + cpatch%FI * cpatch%area * AREA_INV
            hio_fire_area_si(io_si)            = hio_fire_area_si(io_si) + cpatch%frac_burnt * cpatch%area * AREA_INV
            hio_fire_fuel_bulkd_si(io_si)      = hio_fire_fuel_bulkd_si(io_si) + cpatch%fuel_bulkd * cpatch%area * AREA_INV
            hio_fire_fuel_eff_moist_si(io_si)  = hio_fire_fuel_eff_moist_si(io_si) + cpatch%fuel_eff_moist * cpatch%area * AREA_INV
            hio_fire_fuel_sav_si(io_si)        = hio_fire_fuel_sav_si(io_si) + cpatch%fuel_sav * cpatch%area * AREA_INV
            hio_fire_fuel_mef_si(io_si)        = hio_fire_fuel_mef_si(io_si) + cpatch%fuel_mef * cpatch%area * AREA_INV
            hio_sum_fuel_si(io_si)             = hio_sum_fuel_si(io_si) + cpatch%sum_fuel * g_per_kg * cpatch%area * AREA_INV
            
            do i_fuel = 1,nfsc
               hio_litter_moisture_si_fuel(io_si, i_fuel) = hio_litter_moisture_si_fuel(io_si, i_fuel) + &
                    cpatch%litter_moisture(i_fuel) * cpatch%area * AREA_INV

               hio_burnt_frac_litter_si_fuel(io_si, i_fuel) = hio_burnt_frac_litter_si_fuel(io_si, i_fuel) + &
                    cpatch%burnt_frac_litter(i_fuel) * cpatch%frac_burnt * cpatch%area * AREA_INV
            end do


            hio_fire_intensity_area_product_si(io_si) = hio_fire_intensity_area_product_si(io_si) + &
                 cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV

            ! Update Litter Flux Variables

            litt_c       => cpatch%litter(element_pos(carbon12_element))
            flux_diags_c => sites(s)%flux_diags(element_pos(carbon12_element))
                         
            do i_cwd = 1, ncwd

                hio_cwd_ag_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_si_cwdsc(io_si, i_cwd) + &
                      litt_c%ag_cwd(i_cwd)*cpatch%area * AREA_INV * g_per_kg
                hio_cwd_bg_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_si_cwdsc(io_si, i_cwd) + &
                      sum(litt_c%bg_cwd(i_cwd,:)) * cpatch%area * AREA_INV * g_per_kg
                
                hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) + &
                      litt_c%ag_cwd_frag(i_cwd)*cpatch%area * AREA_INV * g_per_kg
                
                hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) + &
                      sum(litt_c%bg_cwd_frag(i_cwd,:)) * cpatch%area * AREA_INV * g_per_kg

            end do

            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         ! divide so-far-just-summed but to-be-averaged patch-age-class variables by patch-age-class area to get mean values
         do ipa2 = 1, nlevage
            if (hio_area_si_age(io_si, ipa2) .gt. tiny) then
               hio_lai_si_age(io_si, ipa2) = hio_lai_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
               hio_ncl_si_age(io_si, ipa2) = hio_ncl_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
               do i_pft = 1, numpft
                  iagepft = ipa2 + (i_pft-1) * nlevage
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
            do i_scls = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_scls
               !
               ! termination mortality. sum of canopy and understory indices
               hio_m6_si_scpf(io_si,i_scpf) = (sites(s)%term_nindivs_canopy(i_scls,i_pft) + &
                                               sites(s)%term_nindivs_ustory(i_scls,i_pft)) * days_per_year

               hio_m6_si_scls(io_si,i_scls) = hio_m6_si_scls(io_si,i_scls) + &
                     (sites(s)%term_nindivs_canopy(i_scls,i_pft) + &
                      sites(s)%term_nindivs_ustory(i_scls,i_pft)) * days_per_year
                     

               !
               ! add termination mortality to canopy and understory mortality
               hio_mortality_canopy_si_scls(io_si,i_scls) = hio_mortality_canopy_si_scls(io_si,i_scls) + &
                    sites(s)%term_nindivs_canopy(i_scls,i_pft) * days_per_year

               hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                    sites(s)%term_nindivs_ustory(i_scls,i_pft) * days_per_year

               hio_mortality_canopy_si_scpf(io_si,i_scpf) = hio_mortality_canopy_si_scpf(io_si,i_scpf) + &
                     sites(s)%term_nindivs_canopy(i_scls,i_pft) * days_per_year

               hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                     sites(s)%term_nindivs_ustory(i_scls,i_pft) * days_per_year

               !
               ! imort on its own
               hio_m4_si_scpf(io_si,i_scpf) = sites(s)%imort_rate(i_scls, i_pft)
               hio_m4_si_scls(io_si,i_scls) = hio_m4_si_scls(io_si,i_scls) + sites(s)%imort_rate(i_scls, i_pft)
               !
               ! add imort to other mortality terms. consider imort as understory mortality even if it happens in 
               ! cohorts that may have been promoted as part of the patch creation, and use the pre-calculated site-level 
               ! values to avoid biasing the results by the dramatically-reduced number densities in cohorts that are subject to imort
               hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                    sites(s)%imort_rate(i_scls, i_pft)
               hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                    sites(s)%imort_rate(i_scls, i_pft)
               !
               iscag = i_scls ! since imort is by definition something that only happens in newly disturbed patches, treat as such
               hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                    sites(s)%imort_rate(i_scls, i_pft)

               ! fire mortality from the site-level diagnostic rates
               hio_m5_si_scpf(io_si,i_scpf) = sites(s)%fmort_rate_canopy(i_scls, i_pft) + &
                     sites(s)%fmort_rate_ustory(i_scls, i_pft)
               hio_m5_si_scls(io_si,i_scls) = hio_m5_si_scls(io_si,i_scls) + &
                     sites(s)%fmort_rate_canopy(i_scls, i_pft) +  sites(s)%fmort_rate_ustory(i_scls, i_pft)
               !
               hio_crownfiremort_si_scpf(io_si,i_scpf) = sites(s)%fmort_rate_crown(i_scls, i_pft)
               hio_cambialfiremort_si_scpf(io_si,i_scpf) = sites(s)%fmort_rate_cambial(i_scls, i_pft)
               !
               ! fire components of overall canopy and understory mortality
               hio_mortality_canopy_si_scpf(io_si,i_scpf) = hio_mortality_canopy_si_scpf(io_si,i_scpf) + &
                    sites(s)%fmort_rate_canopy(i_scls, i_pft)
               hio_mortality_canopy_si_scls(io_si,i_scls) = hio_mortality_canopy_si_scls(io_si,i_scls) + &
                    sites(s)%fmort_rate_canopy(i_scls, i_pft)

               ! the fire mortality rates for each layer are total dead, since the usable
               ! output will then normalize by the counts, we are allowed to sum over layers
               hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                     sites(s)%fmort_rate_ustory(i_scls, i_pft)

               hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                     sites(s)%fmort_rate_ustory(i_scls, i_pft)

               !
               ! carbon flux associated with mortality of trees dying by fire
               hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
                     sites(s)%fmort_carbonflux_canopy
               
               hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
                     sites(s)%fmort_carbonflux_ustory
               
               !
               ! for scag variables, also treat as happening in the newly-disurbed patch

               hio_mortality_canopy_si_scag(io_si,iscag) = hio_mortality_canopy_si_scag(io_si,iscag) + &
                    sites(s)%fmort_rate_canopy(i_scls, i_pft)
               hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                    sites(s)%fmort_rate_ustory(i_scls, i_pft)

               ! while in this loop, pass the fusion-induced growth rate flux to history
               hio_growthflux_fusion_si_scpf(io_si,i_scpf) = hio_growthflux_fusion_si_scpf(io_si,i_scpf) + &
                    sites(s)%growthflux_fusion(i_scls, i_pft) * days_per_year
            end do
         end do
         !
         
         ! treat carbon flux from imort the same way
         hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
              sites(s)%imort_carbonflux
         !
         sites(s)%term_nindivs_canopy(:,:) = 0._r8
         sites(s)%term_nindivs_ustory(:,:) = 0._r8
         sites(s)%imort_carbonflux = 0._r8
         sites(s)%imort_rate(:,:) = 0._r8
         sites(s)%fmort_rate_canopy(:,:) = 0._r8
         sites(s)%fmort_rate_ustory(:,:) = 0._r8
         sites(s)%fmort_carbonflux_canopy = 0._r8
         sites(s)%fmort_carbonflux_ustory = 0._r8
         sites(s)%fmort_rate_cambial(:,:) = 0._r8
         sites(s)%fmort_rate_crown(:,:) = 0._r8
         sites(s)%growthflux_fusion(:,:) = 0._r8

         ! pass the recruitment rate as a flux to the history, and then reset the recruitment buffer
         do i_pft = 1, numpft
            hio_recruitment_si_pft(io_si,i_pft) = sites(s)%recruitment_rate(i_pft) * days_per_year
         end do
         sites(s)%recruitment_rate(:) = 0._r8

         ! summarize all of the mortality fluxes by PFT
         do i_pft = 1, numpft
            do i_scls = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_scls

               hio_mortality_si_pft(io_si,i_pft) = hio_mortality_si_pft(io_si,i_pft) + &
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

            end do
         end do
         
         ! ------------------------------------------------------------------------------
         ! Some carbon only litter diagnostics (legacy)
         ! ------------------------------------------------------------------------------

         flux_diags => sites(s)%flux_diags(element_pos(carbon12_element))

         hio_litter_in_si(io_si) = (sum(flux_diags%cwd_ag_input(:)) + &
              sum(flux_diags%cwd_bg_input(:)) + &
              sum(flux_diags%leaf_litter_input(:)) + &
              sum(flux_diags%root_litter_input(:))) * &
              g_per_kg * AREA_INV * days_per_sec

         hio_litter_out_si(io_si) = 0._r8
         hio_seed_bank_si(io_si)  = 0._r8
         hio_seeds_in_si(io_si)   = 0._r8

         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            litt => cpatch%litter(element_pos(carbon12_element))
            
            area_frac = cpatch%area * AREA_INV
            
            ! Sum up all output fluxes (fragmentation) kgC/m2/day -> gC/m2/s
            hio_litter_out_si(io_si) = hio_litter_out_si(io_si) + &
                 (sum(litt%leaf_fines_frag(:)) + &
                 sum(litt%root_fines_frag(:,:)) + &
                 sum(litt%ag_cwd_frag(:)) + &
                 sum(litt%bg_cwd_frag(:,:))) * &
                 area_frac * g_per_kg * days_per_sec

            ! Sum up total seed bank (germinated and ungerminated)
            hio_seed_bank_si(io_si) = hio_seed_bank_si(io_si) + &
                 (sum(litt%seed(:))+sum(litt%seed_germ(:))) * &
                 area_frac * g_per_kg * days_per_sec

            ! Sum up the input flux into the seed bank (local and external)
            hio_seeds_in_si(io_si) = hio_seeds_in_si(io_si) + &
                 (sum(litt%seed_in_local(:)) + sum(litt%seed_in_extern(:))) * &
                 area_frac * g_per_kg * days_per_sec
            
            cpatch => cpatch%younger
         end do
         

         ! ------------------------------------------------------------------------------
         ! Diagnostics discretized by element type
         ! ------------------------------------------------------------------------------

         hio_cwd_elcwd(io_si,:)   = 0._r8

         
         do el = 1, num_elements
            
            flux_diags => sites(s)%flux_diags(el)
            
            ! Sum up all input litter fluxes (above below, fines, cwd)
            hio_litter_in_elem(io_si, el) =  & 
                 sum(flux_diags%cwd_ag_input(:)) + & 
                 sum(flux_diags%cwd_bg_input(:)) + &
                 sum(flux_diags%leaf_litter_input(:)) + &
                 sum(flux_diags%root_litter_input(:))

            hio_cwd_ag_elem(io_si,el)         = 0._r8
            hio_cwd_bg_elem(io_si,el)         = 0._r8
            hio_fines_ag_elem(io_si,el)       = 0._r8
            hio_fines_bg_elem(io_si,el)       = 0._r8

            hio_seed_bank_elem(io_si,el)      = 0._r8
            hio_seed_germ_elem(io_si,el)      = 0._r8
            hio_seed_decay_elem(io_si,el)     = 0._r8
            hio_seeds_in_local_elem(io_si,el) = 0._r8
            hio_seed_in_extern_elem(io_si,el) = 0._r8
            hio_litter_out_elem(io_si,el)     = 0._r8
            
            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))

               litt => cpatch%litter(el)

               area_frac = cpatch%area * AREA_INV

               ! Sum up all output fluxes (fragmentation)
               hio_litter_out_elem(io_si,el) = hio_litter_out_elem(io_si,el) + &
                    (sum(litt%leaf_fines_frag(:)) + &
                     sum(litt%root_fines_frag(:,:)) + &
                     sum(litt%ag_cwd_frag(:)) + & 
                     sum(litt%bg_cwd_frag(:,:))) * area_frac

               hio_seed_bank_elem(io_si,el) = hio_seed_bank_elem(io_si,el) + & 
                    sum(litt%seed(:)) * area_frac

               hio_seed_germ_elem(io_si,el) = hio_seed_germ_elem(io_si,el) + &
                    sum(litt%seed_germ(:)) * area_frac
                    
               hio_seed_decay_elem(io_si,el) = hio_seed_decay_elem(io_si,el) + & 
                    sum(litt%seed_decay(:)) * area_frac

               hio_seeds_in_local_elem(io_si,el) = hio_seeds_in_local_elem(io_si,el) + & 
                    sum(litt%seed_in_local(:)) * area_frac

               hio_seed_in_extern_elem(io_si,el) = hio_seed_in_extern_elem(io_si,el) + & 
                    sum(litt%seed_in_extern(:)) * area_frac

               ! Litter State Variables
               hio_cwd_ag_elem(io_si,el) = hio_cwd_ag_elem(io_si,el) + &
                     sum(litt%ag_cwd(:)) * area_frac
               
               hio_cwd_bg_elem(io_si,el) = hio_cwd_bg_elem(io_si,el) + &
                     sum(litt%bg_cwd(:,:)) * area_frac
               
               hio_fines_ag_elem(io_si,el) = hio_fines_ag_elem(io_si,el) + & 
                     sum(litt%leaf_fines(:)) * area_frac
               
               hio_fines_bg_elem(io_si,el) = hio_fines_bg_elem(io_si,el) + &
                     sum(litt%root_fines(:,:)) * area_frac


               do cwd=1,ncwd
                   elcwd = (el-1)*ncwd+cwd
                   hio_cwd_elcwd(io_si,elcwd) = hio_cwd_elcwd(io_si,elcwd) + & 
                         (litt%ag_cwd(cwd) + sum(litt%bg_cwd(cwd,:))) * area_frac

               end do

                    
               cpatch => cpatch%younger
            end do

         end do
         



         
         ! pass demotion rates and associated carbon fluxes to history
         do i_scls = 1,nlevsclass
            hio_demotion_rate_si_scls(io_si,i_scls) = sites(s)%demotion_rate(i_scls) * days_per_year
            hio_promotion_rate_si_scls(io_si,i_scls) = sites(s)%promotion_rate(i_scls) * days_per_year
         end do
         !
         ! convert kg C / ha / day to gc / m2 / sec
         hio_demotion_carbonflux_si(io_si) = sites(s)%demotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         hio_promotion_carbonflux_si(io_si) = sites(s)%promotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         !
         ! mortality-associated carbon fluxes
         
         hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
               sites(s)%term_carbonflux_canopy * g_per_kg * days_per_sec * ha_per_m2
         
         hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
               sites(s)%term_carbonflux_ustory * g_per_kg * days_per_sec * ha_per_m2

         ! and zero the site-level termination carbon flux variable
         sites(s)%term_carbonflux_canopy = 0._r8
         sites(s)%term_carbonflux_ustory = 0._r8
         !

         ! add the site-level disturbance-associated cwd and litter input fluxes to thir respective flux fields

         do i_cwd = 1, ncwd
             hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) + &
                   flux_diags_c%cwd_ag_input(i_cwd) * g_per_kg
             
             hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) + &
                   flux_diags_c%cwd_bg_input(i_cwd) * g_per_kg

         end do

         ! and reset the disturbance-related field buffers

         do el = 1, num_elements
             call sites(s)%flux_diags(el)%ZeroFluxDiags()
         end do

      enddo ! site loop
      
    end associate

    return
  end subroutine update_history_dyn
 
 ! ======================================================================================

 subroutine update_history_prod(this,nc,nsites,sites,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use EDTypesMod          , only : nclmax, nlevleaf
    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_area_by_age(nlevage) ! patch area in each bin for normalizing purposes
    real(r8) :: canopy_area_by_age(nlevage) ! canopy area in each bin for normalizing purposes
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    integer  :: ipa2     ! patch incrementer
    integer :: cnlfpft_indx, cnlf_indx, ipft, ican, ileaf ! more iterators and indices
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    real(r8) :: per_dt_tstep          ! Time step in frequency units (/s)

    associate( hio_gpp_si         => this%hvars(ih_gpp_si)%r81d, &
               hio_npp_si         => this%hvars(ih_npp_si)%r81d, &
               hio_aresp_si       => this%hvars(ih_aresp_si)%r81d, &
               hio_maint_resp_si  => this%hvars(ih_maint_resp_si)%r81d, &
               hio_growth_resp_si => this%hvars(ih_growth_resp_si)%r81d, &
               hio_c_stomata_si   => this%hvars(ih_c_stomata_si)%r81d, &
               hio_c_lblayer_si   => this%hvars(ih_c_lblayer_si)%r81d, &
               hio_ar_si_scpf     => this%hvars(ih_ar_si_scpf)%r82d, &
               hio_ar_grow_si_scpf   => this%hvars(ih_ar_grow_si_scpf)%r82d, &
               hio_ar_maint_si_scpf  => this%hvars(ih_ar_maint_si_scpf)%r82d, &
               hio_ar_agsapm_si_scpf => this%hvars(ih_ar_agsapm_si_scpf)%r82d, &
               hio_ar_darkm_si_scpf  => this%hvars(ih_ar_darkm_si_scpf)%r82d, &
               hio_ar_crootm_si_scpf => this%hvars(ih_ar_crootm_si_scpf)%r82d, &
               hio_ar_frootm_si_scpf => this%hvars(ih_ar_frootm_si_scpf)%r82d, &
               hio_gpp_canopy_si     => this%hvars(ih_gpp_canopy_si)%r81d, &
               hio_ar_canopy_si      => this%hvars(ih_ar_canopy_si)%r81d, &
               hio_gpp_understory_si => this%hvars(ih_gpp_understory_si)%r81d, &
               hio_ar_understory_si  => this%hvars(ih_ar_understory_si)%r81d, &
               hio_rdark_canopy_si_scls             => this%hvars(ih_rdark_canopy_si_scls)%r82d, &
               hio_livestem_mr_canopy_si_scls       => this%hvars(ih_livestem_mr_canopy_si_scls)%r82d, &
               hio_livecroot_mr_canopy_si_scls      => this%hvars(ih_livecroot_mr_canopy_si_scls)%r82d, &
               hio_froot_mr_canopy_si_scls          => this%hvars(ih_froot_mr_canopy_si_scls)%r82d, &
               hio_resp_g_canopy_si_scls            => this%hvars(ih_resp_g_canopy_si_scls)%r82d, &
               hio_resp_m_canopy_si_scls            => this%hvars(ih_resp_m_canopy_si_scls)%r82d, &
               hio_rdark_understory_si_scls         => this%hvars(ih_rdark_understory_si_scls)%r82d, &
               hio_livestem_mr_understory_si_scls   => this%hvars(ih_livestem_mr_understory_si_scls)%r82d, &
               hio_livecroot_mr_understory_si_scls  => this%hvars(ih_livecroot_mr_understory_si_scls)%r82d, &
               hio_froot_mr_understory_si_scls      => this%hvars(ih_froot_mr_understory_si_scls)%r82d, &
               hio_resp_g_understory_si_scls        => this%hvars(ih_resp_g_understory_si_scls)%r82d, &
               hio_resp_m_understory_si_scls        => this%hvars(ih_resp_m_understory_si_scls)%r82d, &
               hio_leaf_mr_si         => this%hvars(ih_leaf_mr_si)%r81d, &
               hio_froot_mr_si        => this%hvars(ih_froot_mr_si)%r81d, &
               hio_livecroot_mr_si    => this%hvars(ih_livecroot_mr_si)%r81d, &
               hio_livestem_mr_si     => this%hvars(ih_livestem_mr_si)%r81d, &
               hio_gpp_si_age         => this%hvars(ih_gpp_si_age)%r82d, &
               hio_npp_si_age         => this%hvars(ih_npp_si_age)%r82d, &
               hio_c_stomata_si_age   => this%hvars(ih_c_stomata_si_age)%r82d, &
               hio_c_lblayer_si_age   => this%hvars(ih_c_lblayer_si_age)%r82d, &
               hio_parsun_z_si_cnlf     => this%hvars(ih_parsun_z_si_cnlf)%r82d, &
               hio_parsha_z_si_cnlf     => this%hvars(ih_parsha_z_si_cnlf)%r82d, &
               hio_ts_net_uptake_si_cnlf => this%hvars(ih_ts_net_uptake_si_cnlf)%r82d, &
               hio_parsun_z_si_cnlfpft  => this%hvars(ih_parsun_z_si_cnlfpft)%r82d, &
               hio_parsha_z_si_cnlfpft  => this%hvars(ih_parsha_z_si_cnlfpft)%r82d, &
               hio_laisun_z_si_cnlf     => this%hvars(ih_laisun_z_si_cnlf)%r82d, &
               hio_laisha_z_si_cnlf     => this%hvars(ih_laisha_z_si_cnlf)%r82d, &
               hio_laisun_z_si_cnlfpft  => this%hvars(ih_laisun_z_si_cnlfpft)%r82d, &
               hio_laisha_z_si_cnlfpft  => this%hvars(ih_laisha_z_si_cnlfpft)%r82d, &
               hio_laisun_top_si_can     => this%hvars(ih_laisun_top_si_can)%r82d, &
               hio_laisha_top_si_can     => this%hvars(ih_laisha_top_si_can)%r82d, &
               hio_fabd_sun_si_cnlfpft  => this%hvars(ih_fabd_sun_si_cnlfpft)%r82d, &
               hio_fabd_sha_si_cnlfpft  => this%hvars(ih_fabd_sha_si_cnlfpft)%r82d, &
               hio_fabi_sun_si_cnlfpft  => this%hvars(ih_fabi_sun_si_cnlfpft)%r82d, &
               hio_fabi_sha_si_cnlfpft  => this%hvars(ih_fabi_sha_si_cnlfpft)%r82d, &
               hio_fabd_sun_si_cnlf  => this%hvars(ih_fabd_sun_si_cnlf)%r82d, &
               hio_fabd_sha_si_cnlf  => this%hvars(ih_fabd_sha_si_cnlf)%r82d, &
               hio_fabi_sun_si_cnlf  => this%hvars(ih_fabi_sun_si_cnlf)%r82d, &
               hio_fabi_sha_si_cnlf  => this%hvars(ih_fabi_sha_si_cnlf)%r82d, &
               hio_parprof_dir_si_cnlf  => this%hvars(ih_parprof_dir_si_cnlf)%r82d, &
               hio_parprof_dif_si_cnlf  => this%hvars(ih_parprof_dif_si_cnlf)%r82d, &
               hio_parprof_dir_si_cnlfpft  => this%hvars(ih_parprof_dir_si_cnlfpft)%r82d, &
               hio_parprof_dif_si_cnlfpft  => this%hvars(ih_parprof_dif_si_cnlfpft)%r82d, &
               hio_fabd_sun_top_si_can  => this%hvars(ih_fabd_sun_top_si_can)%r82d, &
               hio_fabd_sha_top_si_can  => this%hvars(ih_fabd_sha_top_si_can)%r82d, &
               hio_fabi_sun_top_si_can  => this%hvars(ih_fabi_sun_top_si_can)%r82d, &
               hio_fabi_sha_top_si_can  => this%hvars(ih_fabi_sha_top_si_can)%r82d, &
               hio_parsun_top_si_can     => this%hvars(ih_parsun_top_si_can)%r82d, &
               hio_parsha_top_si_can     => this%hvars(ih_parsha_top_si_can)%r82d &
               )


      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=2)

      per_dt_tstep = 1.0_r8/dt_tstep

      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1
         
         ipa = 0
         cpatch => sites(s)%oldest_patch

         patch_area_by_age(1:nlevage) = 0._r8
         canopy_area_by_age(1:nlevage) = 0._r8

         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            patch_area_by_age(cpatch%age_class)  = &
                 patch_area_by_age(cpatch%age_class) + cpatch%area

            canopy_area_by_age(cpatch%age_class) = &
                 canopy_area_by_age(cpatch%age_class) + cpatch%total_canopy_area

            ! Canopy resitance terms
            hio_c_stomata_si_age(io_si,cpatch%age_class) = &
                 hio_c_stomata_si_age(io_si,cpatch%age_class) + &
                 cpatch%c_stomata * cpatch%total_canopy_area
            
            hio_c_lblayer_si_age(io_si,cpatch%age_class) = &
                 hio_c_lblayer_si_age(io_si,cpatch%age_class) + &
                 cpatch%c_lblayer * cpatch%total_canopy_area
            
            hio_c_stomata_si(io_si) = hio_c_stomata_si(io_si) + &
                 cpatch%c_stomata * cpatch%total_canopy_area
            
            hio_c_lblayer_si(io_si) = hio_c_lblayer_si(io_si) + &
                 cpatch%c_lblayer * cpatch%total_canopy_area

            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               n_perm2   = ccohort%n * AREA_INV
               
               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  associate( scpf => ccohort%size_by_pft_class, &
                             scls => ccohort%size_class )
                    
                  ! scale up cohort fluxes to the site level
                  hio_npp_si(io_si) = hio_npp_si(io_si) + &
                        ccohort%npp_tstep * g_per_kg * n_perm2 * per_dt_tstep
                  hio_gpp_si(io_si) = hio_gpp_si(io_si) + &
                        ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep
                  hio_aresp_si(io_si) = hio_aresp_si(io_si) + &
                        ccohort%resp_tstep * g_per_kg * n_perm2 * per_dt_tstep
                  hio_growth_resp_si(io_si) = hio_growth_resp_si(io_si) + &
                        ccohort%resp_g * g_per_kg * n_perm2 * per_dt_tstep
                  hio_maint_resp_si(io_si) = hio_maint_resp_si(io_si) + &
                        ccohort%resp_m * g_per_kg * n_perm2 * per_dt_tstep
                  
                  ! aggregate MR fluxes to the site level
                  hio_leaf_mr_si(io_si) = hio_leaf_mr_si(io_si) + ccohort%rdark &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_froot_mr_si(io_si) = hio_froot_mr_si(io_si) + ccohort%froot_mr &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_livecroot_mr_si(io_si) = hio_livecroot_mr_si(io_si) + ccohort%livecroot_mr &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_livestem_mr_si(io_si) = hio_livestem_mr_si(io_si) + ccohort%livestem_mr &
                       * n_perm2 *  sec_per_day * days_per_year

                  ! Total AR (kgC/m2/yr) = (kgC/plant/step) / (s/step) * (plant/m2) * (s/yr)
                  hio_ar_si_scpf(io_si,scpf)    =   hio_ar_si_scpf(io_si,scpf) + &
                        (ccohort%resp_tstep/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Growth AR (kgC/m2/yr)
                  hio_ar_grow_si_scpf(io_si,scpf) = hio_ar_grow_si_scpf(io_si,scpf) + &
                        (ccohort%resp_g/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Maint AR (kgC/m2/yr)
                  hio_ar_maint_si_scpf(io_si,scpf) = hio_ar_maint_si_scpf(io_si,scpf) + &
                        (ccohort%resp_m/dt_tstep) * n_perm2 * sec_per_day * days_per_year
                  
                  ! Maintenance AR partition variables are stored as rates (kgC/plant/s)
                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_agsapm_si_scpf(io_si,scpf) = hio_ar_agsapm_si_scpf(io_si,scpf) + &
                        ccohort%livestem_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_darkm_si_scpf(io_si,scpf) = hio_ar_darkm_si_scpf(io_si,scpf) + &
                        ccohort%rdark * n_perm2 *  sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_crootm_si_scpf(io_si,scpf) = hio_ar_crootm_si_scpf(io_si,scpf) + &
                        ccohort%livecroot_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_frootm_si_scpf(io_si,scpf) = hio_ar_frootm_si_scpf(io_si,scpf) + &
                        ccohort%froot_mr * n_perm2  * sec_per_day * days_per_year


                  ! accumulate fluxes per patch age bin
                  hio_gpp_si_age(io_si,cpatch%age_class) = hio_gpp_si_age(io_si,cpatch%age_class) &
                       + ccohort%gpp_tstep * ccohort%n * g_per_kg * per_dt_tstep
                  hio_npp_si_age(io_si,cpatch%age_class) = hio_npp_si_age(io_si,cpatch%age_class) &
                       + ccohort%npp_tstep * ccohort%n * g_per_kg * per_dt_tstep

                  ! accumulate fluxes on canopy- and understory- separated fluxes
                  if (ccohort%canopy_layer .eq. 1) then
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_canopy_si(io_si) = hio_gpp_canopy_si(io_si) + &
                          ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     hio_ar_canopy_si(io_si) = hio_ar_canopy_si(io_si) + &
                          ccohort%resp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_canopy_si_scls(io_si,scls) = hio_rdark_canopy_si_scls(io_si,scls) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_canopy_si_scls(io_si,scls) = hio_livestem_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_canopy_si_scls(io_si,scls) = hio_livecroot_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_canopy_si_scls(io_si,scls) = hio_froot_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_canopy_si_scls(io_si,scls) = hio_resp_g_canopy_si_scls(io_si,scls) + &
                          ccohort%resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_canopy_si_scls(io_si,scls) = hio_resp_m_canopy_si_scls(io_si,scls) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  else
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_understory_si(io_si) = hio_gpp_understory_si(io_si) + &
                          ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     hio_ar_understory_si(io_si) = hio_ar_understory_si(io_si) + &
                          ccohort%resp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_understory_si_scls(io_si,scls) = hio_rdark_understory_si_scls(io_si,scls) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_understory_si_scls(io_si,scls) = hio_livestem_mr_understory_si_scls(io_si,scls) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_understory_si_scls(io_si,scls) = hio_livecroot_mr_understory_si_scls(io_si,scls) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_understory_si_scls(io_si,scls) = hio_froot_mr_understory_si_scls(io_si,scls) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_understory_si_scls(io_si,scls) = hio_resp_g_understory_si_scls(io_si,scls) + &
                          ccohort%resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_understory_si_scls(io_si,scls) = hio_resp_m_understory_si_scls(io_si,scls) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  endif
                end associate
               endif

               !!! canopy leaf carbon balance
               ican = ccohort%canopy_layer
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) = hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%ts_net_uptake(ileaf) * g_per_kg * per_dt_tstep * ccohort%c_area / AREA
               end do

               ccohort => ccohort%taller
            enddo ! cohort loop

            ! summarize radiation profiles through the canopy
            do ipft=1,numpft
               do ican=1,nclmax         !  cpatch%ncl_p  ?
                  do ileaf=1,nlevleaf   !  cpatch%ncan(ican,ipft) ?
                     ! calculate where we are on multiplexed dimensions
                     cnlfpft_indx = ileaf + (ican-1) * nlevleaf + (ipft-1) * nlevleaf * nclmax 
                     cnlf_indx = ileaf + (ican-1) * nlevleaf
                     !
                     ! first do all the canopy x leaf x pft calculations
                     hio_parsun_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_parsun_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_parsha_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_laisun_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_laisha_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sun_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sha_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sun_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sha_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_parprof_dir_si_cnlfpft(io_si,cnlfpft_indx) = hio_parprof_dir_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%parprof_pft_dir_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parprof_dif_si_cnlfpft(io_si,cnlfpft_indx) = hio_parprof_dif_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%parprof_pft_dif_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     ! summarize across all PFTs
                     hio_parsun_z_si_cnlf(io_si,cnlf_indx) = hio_parsun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_si_cnlf(io_si,cnlf_indx) = hio_parsha_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_si_cnlf(io_si,cnlf_indx) = hio_laisun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_si_cnlf(io_si,cnlf_indx) = hio_laisha_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_si_cnlf(io_si,cnlf_indx) = hio_fabd_sun_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_si_cnlf(io_si,cnlf_indx) = hio_fabd_sha_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_si_cnlf(io_si,cnlf_indx) = hio_fabi_sun_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_si_cnlf(io_si,cnlf_indx) = hio_fabi_sha_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV

                  end do
                  !
                  ! summarize just the top leaf level across all PFTs, for each canopy level
                  hio_parsun_top_si_can(io_si,ican) = hio_parsun_top_si_can(io_si,ican) + &
                       cpatch%ed_parsun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_parsha_top_si_can(io_si,ican) = hio_parsha_top_si_can(io_si,ican) + &
                       cpatch%ed_parsha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_laisun_top_si_can(io_si,ican) = hio_laisun_top_si_can(io_si,ican) + &
                       cpatch%ed_laisun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_laisha_top_si_can(io_si,ican) = hio_laisha_top_si_can(io_si,ican) + &
                       cpatch%ed_laisha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_fabd_sun_top_si_can(io_si,ican) = hio_fabd_sun_top_si_can(io_si,ican) + &
                       cpatch%fabd_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabd_sha_top_si_can(io_si,ican) = hio_fabd_sha_top_si_can(io_si,ican) + &
                       cpatch%fabd_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sun_top_si_can(io_si,ican) = hio_fabi_sun_top_si_can(io_si,ican) + &
                       cpatch%fabi_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sha_top_si_can(io_si,ican) = hio_fabi_sha_top_si_can(io_si,ican) + &
                       cpatch%fabi_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
               end do
            end do

            ! PFT-mean radiation profiles
            do ican=1,nclmax
               do ileaf=1,nlevleaf
                  ! calculate where we are on multiplexed dimensions
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  !
                  hio_parprof_dir_si_cnlf(io_si,cnlf_indx) = hio_parprof_dir_si_cnlf(io_si,cnlf_indx) + &
                       cpatch%parprof_dir_z(ican,ileaf) * cpatch%area * AREA_INV
                  hio_parprof_dif_si_cnlf(io_si,cnlf_indx) = hio_parprof_dif_si_cnlf(io_si,cnlf_indx) + &
                       cpatch%parprof_dif_z(ican,ileaf) * cpatch%area * AREA_INV
               end do
            end do

            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         do ipa2 = 1, nlevage
            if (patch_area_by_age(ipa2) .gt. tiny) then
               hio_gpp_si_age(io_si, ipa2) = hio_gpp_si_age(io_si, ipa2) / (patch_area_by_age(ipa2))
               hio_npp_si_age(io_si, ipa2) = hio_npp_si_age(io_si, ipa2) / (patch_area_by_age(ipa2))
            else
               hio_gpp_si_age(io_si, ipa2) = 0._r8
               hio_npp_si_age(io_si, ipa2) = 0._r8
            endif

            ! Normalize resistance diagnostics
            if (canopy_area_by_age(ipa2) .gt. tiny) then
               hio_c_stomata_si_age(io_si,ipa2) = &
                    hio_c_stomata_si_age(io_si,ipa2) / canopy_area_by_age(ipa2)

               hio_c_lblayer_si_age(io_si,ipa2) = &
                    hio_c_lblayer_si_age(io_si,ipa2) / canopy_area_by_age(ipa2)
            else
               hio_c_stomata_si_age(io_si,ipa2) = 0._r8
               hio_c_lblayer_si_age(io_si,ipa2) = 0._r8
            end if
            
         end do
         
         ! Normalize resistance diagnostics
         if ( sum(canopy_area_by_age(1:nlevage)) .gt. tiny) then
            hio_c_stomata_si(io_si) = hio_c_stomata_si(io_si) / sum(canopy_area_by_age(1:nlevage))
            hio_c_lblayer_si(io_si) = hio_c_lblayer_si(io_si) / sum(canopy_area_by_age(1:nlevage))
         else
            hio_c_stomata_si(io_si) = 0._r8
            hio_c_lblayer_si(io_si) = 0._r8
         end if
 
      enddo ! site loop

    end associate
 
  end subroutine update_history_prod

  ! =====================================================================================

  subroutine update_history_hydraulics(this,nc,nsites,sites,bc_in,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use FatesHydraulicsMemMod, only : ed_cohort_hydr_type, nshell
    use FatesHydraulicsMemMod, only : ed_site_hydr_type
    use EDTypesMod           , only : maxpft

    
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
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: ft               ! functional type index
    integer  :: scpf
!    integer  :: io_shsl  ! The combined "SH"ell "S"oil "L"ayer index in the IO array
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    real(r8) :: ncohort_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
                                                   ! should be "hio_nplant_si_scpf"
    real(r8) :: number_fraction
    real(r8) :: number_fraction_rate
    real(r8) :: mean_aroot
    integer  :: ipa2     ! patch incrementer
    integer  :: iscpf    ! index of the scpf group
    integer  :: ipft     ! index of the pft loop
    integer  :: iscls    ! index of the size-class loop
    integer  :: k        ! rhizosphere shell index
    integer  :: jsoil    ! soil layer index
    integer  :: jrhiz    ! rhizosphere layer index
    integer  :: jr1, jr2 ! Rhizosphere top and bottom layers
    integer  :: nlevrhiz ! number of rhizosphere layers
    real(r8) :: mean_soil_vwc    ! mean soil volumetric water content [m3/m3]
    real(r8) :: mean_soil_vwcsat ! mean soil saturated volumetric water content [m3/m3]
    real(r8) :: mean_soil_matpot ! mean soil water potential [MPa]
    real(r8) :: layer_areaweight ! root area weighting factor for each soil layer
    real(r8) :: areaweight       ! root area weighting factor for column
    real(r8) :: vwc              ! volumetric water content of layer [m3/m3] = theta
    real(r8) :: vwc_sat          ! saturated water content of layer [m3/m3]
    real(r8) :: psi              ! matric potential of soil layer
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: site_hydr

    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! Should this be 365.25?

    
    if(hlm_use_planthydro.eq.ifalse) return

    associate( hio_errh2o_scpf  => this%hvars(ih_errh2o_scpf)%r82d, &
          hio_tran_scpf         => this%hvars(ih_tran_scpf)%r82d, &
          hio_sapflow_scpf      => this%hvars(ih_sapflow_scpf)%r82d, &
          hio_sapflow_si        => this%hvars(ih_sapflow_si)%r81d, & 
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
          hio_h2oveg_si         => this%hvars(ih_h2oveg_si)%r81d, &
          hio_nplant_si_scpf    => this%hvars(ih_nplant_si_scpf)%r82d, &
          hio_nplant_si_capf    => this%hvars(ih_nplant_si_capf)%r82d, &
          hio_h2oveg_hydro_err_si   => this%hvars(ih_h2oveg_hydro_err_si)%r81d, &
          hio_rootwgt_soilvwc_si    => this%hvars(ih_rootwgt_soilvwc_si)%r81d, &
          hio_rootwgt_soilvwcsat_si => this%hvars(ih_rootwgt_soilvwcsat_si)%r81d, & 
          hio_rootwgt_soilmatpot_si => this%hvars(ih_rootwgt_soilmatpot_si)%r81d, &
          hio_soilmatpot_sl         => this%hvars(ih_soilmatpot_sl)%r82d, &
          hio_soilvwc_sl            => this%hvars(ih_soilvwc_sl)%r82d, &
          hio_soilvwcsat_sl         => this%hvars(ih_soilvwcsat_sl)%r82d, &
          hio_rootuptake_si         => this%hvars(ih_rootuptake_si)%r81d, &
          hio_rootuptake_sl         => this%hvars(ih_rootuptake_sl)%r82d, &
          hio_rootuptake0_scpf      => this%hvars(ih_rootuptake0_scpf)%r82d, &
          hio_rootuptake10_scpf     => this%hvars(ih_rootuptake10_scpf)%r82d, &
          hio_rootuptake50_scpf     => this%hvars(ih_rootuptake50_scpf)%r82d, &
          hio_rootuptake100_scpf    => this%hvars(ih_rootuptake100_scpf)%r82d )

      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=4)
      
      do s = 1,nsites

         site_hydr => sites(s)%si_hydr
         nlevrhiz = site_hydr%nlevrhiz
         jr1 = site_hydr%i_rhiz_t
         jr2 = site_hydr%i_rhiz_b

         io_si  = this%iovar_map(nc)%site_index(s)
         
         hio_h2oveg_si(io_si)              = site_hydr%h2oveg
         hio_h2oveg_hydro_err_si(io_si)    = site_hydr%h2oveg_hydro_err

         ncohort_scpf(:) = 0.0_r8  ! Counter for normalizing weighting 
                                   ! factors for cohort mean propoerties
                                   ! This is actually used as a check
                                   ! on hio_nplant_si_scpf
         
         ! Get column means of some soil diagnostics, these are weighted
         ! by the amount of fine-root surface area in each layer
         ! --------------------------------------------------------------------
         
         mean_soil_vwc    = 0._r8
         mean_soil_matpot = 0._r8
         mean_soil_vwcsat = 0._r8
         areaweight       = 0._r8
         
         do jrhiz=1,nlevrhiz
            
            jsoil = jrhiz + jr1-1
            vwc     = bc_in(s)%h2o_liqvol_sl(jsoil)
            psi     = site_hydr%wrf_soil(jrhiz)%p%psi_from_th(vwc)
            vwc_sat = bc_in(s)%watsat_sl(jsoil)
            layer_areaweight = site_hydr%l_aroot_layer(jrhiz)*pi_const*site_hydr%rs1(jrhiz)**2.0
            mean_soil_vwc    = mean_soil_vwc + vwc*layer_areaweight
            mean_soil_vwcsat = mean_soil_vwcsat + vwc_sat*layer_areaweight
            mean_soil_matpot = mean_soil_matpot + psi*layer_areaweight
            areaweight       = areaweight + layer_areaweight

            hio_soilmatpot_sl(io_si,jsoil) = psi
            hio_soilvwc_sl(io_si,jsoil)    = vwc
            hio_soilvwcsat_sl(io_si,jsoil) = vwc_sat
            
         end do
         
         hio_rootwgt_soilvwc_si(io_si)    = mean_soil_vwc/areaweight
         hio_rootwgt_soilvwcsat_si(io_si) = mean_soil_vwcsat/areaweight
         hio_rootwgt_soilmatpot_si(io_si) = mean_soil_matpot/areaweight
         
         hio_rootuptake_si(io_si) = sum(site_hydr%rootuptake_sl,dim=1)
         hio_rootuptake_sl(io_si,:) = 0._r8
         hio_rootuptake_sl(io_si,jr1:jr2) = site_hydr%rootuptake_sl(1:nlevrhiz)
         hio_rootuptake_si(io_si) = sum(site_hydr%sapflow_scpf)
         
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               if ( .not. ccohort%isnew ) then
                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  ncohort_scpf(iscpf) = ncohort_scpf(iscpf) + ccohort%n
               end if
               ccohort => ccohort%taller
            enddo ! cohort loop
            cpatch => cpatch%younger
         end do !patch loop

         do ipft = 1, numpft
            do iscls = 1,nlevsclass
               iscpf = (ipft-1)*nlevsclass + iscls
               hio_sapflow_scpf(io_si,iscpf)       = site_hydr%sapflow_scpf(iscls, ipft)
               hio_rootuptake0_scpf(io_si,iscpf)   = site_hydr%rootuptake0_scpf(iscls,ipft)
               hio_rootuptake10_scpf(io_si,iscpf)  = site_hydr%rootuptake10_scpf(iscls,ipft)
               hio_rootuptake50_scpf(io_si,iscpf)  = site_hydr%rootuptake50_scpf(iscls,ipft)
               hio_rootuptake100_scpf(io_si,iscpf) = site_hydr%rootuptake100_scpf(iscls,ipft)
            end do
         end do

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            ccohort => cpatch%shortest
            do while(associated(ccohort))

               ccohort_hydr => ccohort%co_hydr
               
               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  
                  ! scale up cohort fluxes to their sites
                  number_fraction_rate = (ccohort%n / ncohort_scpf(iscpf))/dt_tstep
                  
                  ! scale cohorts to mean quantity
                  number_fraction = (ccohort%n / ncohort_scpf(iscpf))
                  
                  hio_errh2o_scpf(io_si,iscpf) = hio_errh2o_scpf(io_si,iscpf) + &
                        ccohort_hydr%errh2o * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_tran_scpf(io_si,iscpf) = hio_tran_scpf(io_si,iscpf) + &
                        (ccohort_hydr%qtop) * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_iterh1_scpf(io_si,iscpf)          = hio_iterh1_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh1  * number_fraction             ! [-]
                  
                  hio_iterh2_scpf(io_si,iscpf)          = hio_iterh2_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh2 * number_fraction             ! [-]

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
                       mean_aroot * number_fraction     ! [MPa]
                  
                  hio_twp_scpf(io_si,iscpf)             = hio_twp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_troot  * number_fraction       ! [MPa]
                  
                  hio_swp_scpf(io_si,iscpf)             = hio_swp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_ag(2)  * number_fraction       ! [MPa]
                  
                  hio_lwp_scpf(io_si,iscpf)             = hio_lwp_scpf(io_si,iscpf) + &
                       ccohort_hydr%psi_ag(1)  * number_fraction       ! [MPa]

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
            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         if(hlm_use_ed_st3.eq.ifalse) then
            do scpf=1,nlevsclass*numpft
               if( abs(hio_nplant_si_scpf(io_si, scpf)-ncohort_scpf(scpf)) > 1.0E-8_r8 ) then
                  write(fates_log(),*) 'numpft:',numpft
                  write(fates_log(),*) 'nlevsclass:',nlevsclass
                  write(fates_log(),*) 'scpf:',scpf
                  write(fates_log(),*) 'io_si:',io_si
                  write(fates_log(),*) 'hio_nplant_si_scpf:',hio_nplant_si_scpf(io_si, scpf)
                  write(fates_log(),*) 'ncohort_scpf:',ncohort_scpf(scpf)
                  write(fates_log(),*) 'nplant check on hio_nplant_si_scpf fails during hydraulics history updates'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
            end do
         end if

      enddo ! site loop

    end associate
 
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
    ! This subroutine is called in two contexts, either in count mode or inialize mode
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

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8    
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_pft_r8, site_coage_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesInterfaceTypesMod     , only : hlm_use_planthydro
    
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8


    implicit none
    
    class(fates_history_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?

    integer :: ivar
    character(len=10) :: tempstring 
    
    ivar=0
    
    ! Site level counting variables
    call this%set_history_var(vname='ED_NPATCHES', units='none',                &
         long='Total number of ED patches per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches_si)

    call this%set_history_var(vname='ED_NCOHORTS', units='none',                &
         long='Total number of ED cohorts per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_ncohorts_si)
    
    ! Patch variables
    call this%set_history_var(vname='TRIMMING', units='none',                   &
         long='Degree to which canopy expansion is limited by leaf economics',  & 
         use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_trimming_si)
    
    call this%set_history_var(vname='AREA_PLANT', units='m2',                   &
         long='area occupied by all plants', use_default='active',              &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_plant_si)
    
    call this%set_history_var(vname='AREA_TREES', units='m2',                   &
         long='area occupied by woody plants', use_default='active',            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_trees_si)

    call this%set_history_var(vname='SITE_COLD_STATUS', units='0,1,2', &
          long='Site level cold status, 0=not cold-dec, 1=too cold for leaves, 2=not-too cold',  &
          use_default='active',                                                  &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
          ivar=ivar, initialize=initialize_variables, index = ih_site_cstatus_si )

    call this%set_history_var(vname='SITE_DROUGHT_STATUS', units='0,1,2,3', &
          long='Site level drought status, <2 too dry for leaves, >=2 not-too dry', &
          use_default='active',                                                  &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
          ivar=ivar, initialize=initialize_variables, index = ih_site_dstatus_si)

    call this%set_history_var(vname='SITE_GDD', units='degC',  &
         long='site level growing degree days',                &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_gdd_si)
    
    call this%set_history_var(vname='SITE_NCHILLDAYS', units = 'days', &
         long='site level number of chill days', &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_site_nchilldays_si)

    call this%set_history_var(vname='SITE_NCOLDDAYS', units = 'days', &
         long='site level number of cold days', &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_site_ncolddays_si)

    call this%set_history_var(vname='SITE_DAYSINCE_COLDLEAFOFF', units='days', &
         long='site level days elapsed since cold leaf drop', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_cleafoff_si)

    call this%set_history_var(vname='SITE_DAYSINCE_COLDLEAFON', units='days', &
         long='site level days elapsed since cold leaf flush', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_cleafon_si) 

    call this%set_history_var(vname='SITE_DAYSINCE_DROUGHTLEAFOFF', units='days', &
         long='site level days elapsed since drought leaf drop', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_dleafoff_si)
    
    call this%set_history_var(vname='SITE_DAYSINCE_DROUGHTLEAFON', units='days', &
         long='site level days elapsed since drought leaf flush', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_dleafon_si)

    call this%set_history_var(vname='SITE_MEANLIQVOL_DROUGHTPHEN', units='m3/m3', &
         long='site level mean liquid water volume for drought phen', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_meanliqvol_si)

    call this%set_history_var(vname='CANOPY_SPREAD', units='0-1',               &
         long='Scaling factor between tree basal area and canopy area',         &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_spread_si)

    call this%set_history_var(vname='PFTbiomass', units='gC/m2',                   &
         long='total PFT level biomass', use_default='active',                     &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_pft )

    call this%set_history_var(vname='PFTleafbiomass', units='gC/m2',              &
         long='total PFT level leaf biomass', use_default='active',                &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_leafbiomass_si_pft )

    call this%set_history_var(vname='PFTstorebiomass',  units='gC/m2',            &
         long='total PFT level stored biomass', use_default='active',              &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_storebiomass_si_pft )

    call this%set_history_var(vname='PFTcrownarea',  units='m2/ha',            &
         long='total PFT level crown area', use_default='inactive',              &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_si_pft )
    
    call this%set_history_var(vname='PFTnindivs',  units='indiv / m2',            &
         long='total PFT level number of individuals', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_nindivs_si_pft )

    call this%set_history_var(vname='RECRUITMENT',  units='indiv/ha/yr',            &
         long='Rate of recruitment by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_recruitment_si_pft )

    call this%set_history_var(vname='MORTALITY',  units='indiv/ha/yr',            &
         long='Rate of total mortality by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_mortality_si_pft )

    ! patch age class variables
    call this%set_history_var(vname='PATCH_AREA_BY_AGE', units='m2/m2',             &
         long='patch area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_area_si_age )

    call this%set_history_var(vname='LAI_BY_AGE', units='m2/m2',                   &
         long='leaf area index by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_lai_si_age )

    call this%set_history_var(vname='CANOPY_AREA_BY_AGE', units='m2/m2',             &
         long='canopy area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_area_si_age )
    
    call this%set_history_var(vname='NCL_BY_AGE', units='--',                   &
         long='number of canopy levels by age bin', use_default='inactive',             &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_ncl_si_age )

    call this%set_history_var(vname='NPATCH_BY_AGE', units='--',                   &
         long='number of patches by age bin', use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches_si_age )

    if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
       tempstring = 'active'
    else
       tempstring = 'inactive'
    endif
    call this%set_history_var(vname='ZSTAR_BY_AGE', units='m',                   &
         long='product of zstar and patch area by age bin (divide by PATCH_AREA_BY_AGE to get mean zstar)', &
         use_default=trim(tempstring),                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_zstar_si_age )

    call this%set_history_var(vname='CANOPY_HEIGHT_DIST', units='m2/m2',                   &
         long='canopy height distribution', use_default='active',                     &
         avgflag='A', vtype=site_height_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_height_dist_si_height )

    call this%set_history_var(vname='LEAF_HEIGHT_DIST', units='m2/m2',                   &
         long='leaf height distribution', use_default='active',                     &
         avgflag='A', vtype=site_height_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_leaf_height_dist_si_height )

    call this%set_history_var(vname='BIOMASS_BY_AGE', units='kgC/m2',                   &
         long='Total Biomass within a given patch age bin', &
         use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_age )

    ! Secondary forest area and age diagnostics

    call this%set_history_var(vname='SECONDARY_FOREST_FRACTION', units='m2/m2', &
         long='Secondary forest fraction', use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fraction_secondary_forest_si )

    call this%set_history_var(vname='WOOD_PRODUCT', units='gC/m2', &
         long='Total wood product from logging', use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_woodproduct_si )

    call this%set_history_var(vname='SECONDARY_FOREST_BIOMASS', units='kgC/m2', &
         long='Biomass on secondary lands (per total site area, mult by SECONDARY_FOREST_FRACTION to get per secondary forest area)',&
         use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_secondary_forest_si )

    call this%set_history_var(vname='SECONDARY_AREA_AGE_ANTHRO_DIST', units='m2/m2', &
         long='Secondary forest patch area age distribution since anthropgenic disturbance', &
         use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_agesince_anthrodist_si_age )

    call this%set_history_var(vname='SECONDARY_AREA_PATCH_AGE_DIST', units='m2/m2', &
         long='Secondary forest patch area age distribution since any kind of disturbance', &
         use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_secondaryforest_area_si_age )


    ! Fire Variables

    call this%set_history_var(vname='FIRE_NESTEROV_INDEX', units='none',       &
         long='nesterov_fire_danger index', use_default='active',               &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_nesterov_fire_danger_si)

    call this%set_history_var(vname='FIRE_ROS', units='m/min',                 &
         long='fire rate of spread m/min', use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_spitfire_ros_si)

    call this%set_history_var(vname='FIRE_ROS_AREA_PRODUCT', units='m/min',                 &
         long='product of fire rate of spread (m/min) and burned area (fraction)--divide by FIRE_AREA to get burned-area-weighted-mean ROS', use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_ros_area_product_si)

    call this%set_history_var(vname='EFFECT_WSPEED', units='none',             &
         long ='effective windspeed for fire spread', use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_effect_wspeed_si )

    call this%set_history_var(vname='FIRE_TFC_ROS', units='kgC/m2',              &
         long ='total fuel consumed', use_default='active',                     &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_tfc_ros_si )

    call this%set_history_var(vname='FIRE_TFC_ROS_AREA_PRODUCT', units='kgC/m2',              &
         long ='product of total fuel consumed and burned area--divide by FIRE_AREA to get burned-area-weighted-mean TFC', use_default='active',                     &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_tfc_ros_area_product_si )

    call this%set_history_var(vname='FIRE_INTENSITY', units='kJ/m/s',          &
         long='spitfire fire intensity: kJ/m/s', use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_si )

    call this%set_history_var(vname='FIRE_INTENSITY_AREA_PRODUCT', units='kJ/m/s',          &
         long='spitfire product of fire intensity and burned area (divide by FIRE_AREA to get area-weighted mean intensity)', &
         use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_area_product_si )

    call this%set_history_var(vname='FIRE_AREA', units='fraction',             &
         long='spitfire fire area burn fraction', use_default='active',                    &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_area_si )

    call this%set_history_var(vname='FIRE_FUEL_MEF', units='m',                &
         long='spitfire fuel moisture',  use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_mef_si )

    call this%set_history_var(vname='FIRE_FUEL_BULKD', units='kg biomass/m3',              &
         long='spitfire fuel bulk density',  use_default='active',              &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_bulkd_si )

    call this%set_history_var(vname='FIRE_FUEL_EFF_MOIST', units='m',          &
         long='spitfire fuel moisture', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_eff_moist_si )

    call this%set_history_var(vname='FIRE_FUEL_SAV', units='per m',                &
         long='spitfire fuel surface/volume ',  use_default='active',           &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_sav_si )

    call this%set_history_var(vname='SUM_FUEL', units='gC m-2',                &
         long='total ground fuel related to ros (omits 1000hr fuels)',          & 
         use_default='active',                                                  & 
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_sum_fuel_si )

    call this%set_history_var(vname='FUEL_MOISTURE_NFSC', units='-',                &
         long='spitfire size-resolved fuel moisture', use_default='active',       &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_moisture_si_fuel )

    call this%set_history_var(vname='AREA_BURNT_BY_PATCH_AGE', units='m2/m2', &
         long='spitfire area burnt by patch age (divide by patch_area_by_age to get burnt fraction by age)', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_area_burnt_si_age )

    call this%set_history_var(vname='FIRE_INTENSITY_BY_PATCH_AGE', units='kJ/m/2', &
         long='product of fire intensity and burned area, resolved by patch age (so divide by AREA_BURNT_BY_PATCH_AGE to get burned-area-weighted-average intensity', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_si_age )

    call this%set_history_var(vname='SUM_FUEL_BY_PATCH_AGE', units='gC / m2 of site area', &
         long='spitfire ground fuel related to ros (omits 1000hr fuels) within each patch age bin (divide by patch_area_by_age to get fuel per unit area of that-age patch)', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_sum_fuel_si_age )

    call this%set_history_var(vname='BURNT_LITTER_FRAC_AREA_PRODUCT', units='fraction', &
         long='product of fraction of fuel burnt and burned area (divide by FIRE_AREA to get burned-area-weighted mean fraction fuel burnt)', &
         use_default='active', &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_burnt_frac_litter_si_fuel )


    ! Litter Variables

    call this%set_history_var(vname='LITTER_IN', units='gC m-2 s-1',           &
         long='FATES litter flux in',  use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_in_si )

    call this%set_history_var(vname='LITTER_OUT', units='gC m-2 s-1',          &
         long='FATES litter flux out',  use_default='active',                  & 
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_out_si )

    call this%set_history_var(vname='SEED_BANK', units='gC m-2',               &
         long='Total Seed Mass of all PFTs',  use_default='active',             &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_bank_si )

    call this%set_history_var(vname='SEEDS_IN', units='gC m-2 s-1',            &
         long='Seed Production Rate',  use_default='active',                    &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_si )

    call this%set_history_var(vname='LITTER_IN_ELEM', units='kg m-2 d-1',         &
         long='FATES litter flux in',  use_default='active',                      &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_in_elem )

    call this%set_history_var(vname='LITTER_OUT_ELEM', units='kg m-2 d-1',         &
         long='FATES litter flux out (fragmentation only)',  use_default='active',                      & 
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_out_elem )

    call this%set_history_var(vname='SEED_BANK_ELEM', units='kg m-2',             &
         long='Total Seed Mass of all PFTs',  use_default='active',               &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_bank_elem )

    call this%set_history_var(vname='SEEDS_IN_LOCAL_ELEM', units='kg m-2 d-1',     &
         long='Within Site Seed Production Rate',  use_default='active',           &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_local_elem )

    call this%set_history_var(vname='SEEDS_IN_EXTERN_ELEM', units='kg m-2 d-1',     &
         long='External Seed Influx Rate',  use_default='active',                   &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_extern_elem )

    call this%set_history_var(vname='SEED_GERM_ELEM', units='kg m-2 d-1',          &
         long='Seed mass converted into new cohorts', use_default='active',        &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_germ_elem )

    call this%set_history_var(vname='SEED_DECAY', units='kg m-2 d-1',           &
         long='Seed mass decay', use_default='active',                          &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_decay_elem )
    
    call this%set_history_var(vname='ED_bstore', units='gC m-2',                  &
         long='Storage biomass', use_default='active',                          &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bstore_si )

    call this%set_history_var(vname='ED_bdead', units='gC m-2',                   &
         long='Dead (structural) biomass (live trees, not CWD)',                &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bdead_si )

    call this%set_history_var(vname='ED_balive', units='gC m-2',                  &
         long='Live biomass', use_default='active',                             &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_balive_si )

    call this%set_history_var(vname='ED_bleaf', units='gC m-2',                   &
         long='Leaf biomass',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bleaf_si )

    call this%set_history_var(vname='ED_bsapwood', units='gC m-2',                 &
         long='Sapwood biomass',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bsapwood_si )    

    call this%set_history_var(vname='ED_bfineroot', units='gC m-2',                 &
         long='Fine root biomass',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bfineroot_si )    

    call this%set_history_var(vname='ED_biomass', units='gC m-2',                  &
         long='Total biomass',  use_default='active',                           &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_btotal_si )

    call this%set_history_var(vname='AGB', units='gC m-2',                  &
         long='Aboveground biomass',  use_default='active',                           &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_agb_si )

    call this%set_history_var(vname='BIOMASS_CANOPY', units='gC m-2',                   &
         long='Biomass of canopy plants',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_biomass_si )

    call this%set_history_var(vname='BIOMASS_UNDERSTORY', units='gC m-2',                   &
         long='Biomass of understory plants',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_understory_biomass_si )

    ! Canopy Resistance 

    call this%set_history_var(vname='C_STOMATA', units='umol m-2 s-1',                   &
         long='mean stomatal conductance', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_stomata_si )

    call this%set_history_var(vname='C_LBLAYER', units='umol m-2 s-1',                   &
         long='mean leaf boundary layer conductance', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_lblayer_si )


    ! Ecosystem Carbon Fluxes (updated rapidly, upfreq=2)

    call this%set_history_var(vname='NPP', units='gC/m^2/s',                &
         long='net primary production',  use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_si )

    call this%set_history_var(vname='GPP', units='gC/m^2/s',                   &
         long='gross primary production',  use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_si )

    call this%set_history_var(vname='AR', units='gC/m^2/s',                 &
         long='autotrophic respiration', use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_aresp_si )

    call this%set_history_var(vname='GROWTH_RESP', units='gC/m^2/s',           &
         long='growth respiration', use_default='active',                       &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_growth_resp_si )

    call this%set_history_var(vname='MAINT_RESP', units='gC/m^2/s',            &
         long='maintenance respiration', use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_maint_resp_si )

    ! Canopy resistance 

    call this%set_history_var(vname='C_STOMATA_BY_AGE', units='umol m-2 s-1',                   &
         long='mean stomatal conductance - by patch age', use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_stomata_si_age )

    call this%set_history_var(vname='C_LBLAYER_BY_AGE', units='umol m-2 s-1',                   &
         long='mean leaf boundary layer conductance - by patch age', use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_lblayer_si_age )

    ! fast fluxes by age bin
    call this%set_history_var(vname='NPP_BY_AGE', units='gC/m^2/s',                   &
         long='net primary productivity by age bin', use_default='inactive',           &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_si_age )

    call this%set_history_var(vname='GPP_BY_AGE', units='gC/m^2/s',                   &
         long='gross primary productivity by age bin', use_default='inactive',         &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_si_age )

    ! fast fluxes separated canopy/understory
    call this%set_history_var(vname='GPP_CANOPY', units='gC/m^2/s',                   &
         long='gross primary production of canopy plants',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy_si )

    call this%set_history_var(vname='AR_CANOPY', units='gC/m^2/s',                 &
         long='autotrophic respiration of canopy plants', use_default='active',       &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy_si )

    call this%set_history_var(vname='GPP_UNDERSTORY', units='gC/m^2/s',                   &
         long='gross primary production of understory plants',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory_si )

    call this%set_history_var(vname='AR_UNDERSTORY', units='gC/m^2/s',                 &
         long='autotrophic respiration of understory plants', use_default='active',       &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_understory_si )


    ! fast radiative fluxes resolved through the canopy
    call this%set_history_var(vname='PARSUN_Z_CNLF', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_si_cnlf )

    call this%set_history_var(vname='PARSHA_Z_CNLF', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_si_cnlf )

    call this%set_history_var(vname='PARSUN_Z_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_si_cnlfpft )

    call this%set_history_var(vname='PARSHA_Z_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_si_cnlfpft )

    call this%set_history_var(vname='PARSUN_Z_CAN', units='W/m2',                 &
         long='PAR absorbed in the sun by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_top_si_can )

    call this%set_history_var(vname='PARSHA_Z_CAN', units='W/m2',                 &
         long='PAR absorbed in the shade by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_top_si_can )

    call this%set_history_var(vname='LAISUN_Z_CNLF', units='m2/m2',                 &
         long='LAI in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_si_cnlf )

    call this%set_history_var(vname='LAISHA_Z_CNLF', units='m2/m2',                 &
         long='LAI in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_si_cnlf )

    call this%set_history_var(vname='LAISUN_Z_CNLFPFT', units='m2/m2',                 &
         long='LAI in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_si_cnlfpft )

    call this%set_history_var(vname='LAISHA_Z_CNLFPFT', units='m2/m2',                 &
         long='LAI in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_si_cnlfpft )

    call this%set_history_var(vname='LAISUN_TOP_CAN', units='m2/m2',                 &
         long='LAI in the sun by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_top_si_can )

    call this%set_history_var(vname='LAISHA_TOP_CAN', units='m2/m2',                 &
         long='LAI in the shade by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_top_si_can )

    call this%set_history_var(vname='FABD_SUN_CNLFPFT', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_si_cnlfpft )

    call this%set_history_var(vname='FABD_SHA_CNLFPFT', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_si_cnlfpft )

    call this%set_history_var(vname='FABI_SUN_CNLFPFT', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_si_cnlfpft )

    call this%set_history_var(vname='FABI_SHA_CNLFPFT', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_si_cnlfpft )

    call this%set_history_var(vname='FABD_SUN_CNLF', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_si_cnlf )

    call this%set_history_var(vname='FABD_SHA_CNLF', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_si_cnlf )

    call this%set_history_var(vname='FABI_SUN_CNLF', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_si_cnlf )

    call this%set_history_var(vname='FABI_SHA_CNLF', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_si_cnlf )

    call this%set_history_var(vname='PARPROF_DIR_CNLFPFT', units='W/m2',                 &
         long='Radiative profile of direct PAR through each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dir_si_cnlfpft )

    call this%set_history_var(vname='PARPROF_DIF_CNLFPFT', units='W/m2',                 &
         long='Radiative profile of diffuse PAR through each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dif_si_cnlfpft )

    call this%set_history_var(vname='PARPROF_DIR_CNLF', units='W/m2',                 &
         long='Radiative profile of direct PAR through each canopy and leaf layer (averaged across PFTs)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dir_si_cnlf )

    call this%set_history_var(vname='PARPROF_DIF_CNLF', units='W/m2',                 &
         long='Radiative profile of diffuse PAR through each canopy and leaf layer (averaged across PFTs)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dif_si_cnlf )

    call this%set_history_var(vname='FABD_SUN_TOPLF_BYCANLAYER', units='fraction',                 &
         long='sun fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_top_si_can )

    call this%set_history_var(vname='FABD_SHA_TOPLF_BYCANLAYER', units='fraction',                 &
         long='shade fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_top_si_can )

    call this%set_history_var(vname='FABI_SUN_TOPLF_BYCANLAYER', units='fraction',                 &
         long='sun fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_top_si_can )

    call this%set_history_var(vname='FABI_SHA_TOPLF_BYCANLAYER', units='fraction',                 &
         long='shade fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_top_si_can )

    !!! canopy-resolved fluxes and structure
    call this%set_history_var(vname='NET_C_UPTAKE_CNLF', units='gC/m2/s',                 &
         long='net carbon uptake by each canopy and leaf layer per unit ground area (i.e. divide by CROWNAREA_CNLF to make per leaf area)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ts_net_uptake_si_cnlf )

    call this%set_history_var(vname='CROWNAREA_CNLF', units='m2/m2',                 &
         long='total crown area that is occupied by leaves in each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_si_cnlf )

    call this%set_history_var(vname='CROWNAREA_CAN', units='m2/m2',                 &
         long='total crown area in each canopy layer', &
         use_default='active',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_si_can )

    ! slow carbon fluxes associated with mortality from or transfer betweeen canopy and understory
    call this%set_history_var(vname='DEMOTION_CARBONFLUX', units = 'gC/m2/s',               &
          long='demotion-associated biomass carbon flux from canopy to understory', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_carbonflux_si )

    call this%set_history_var(vname='PROMOTION_CARBONFLUX', units = 'gC/m2/s',               &
          long='promotion-associated biomass carbon flux from understory to canopy', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_carbonflux_si )

    call this%set_history_var(vname='MORTALITY_CARBONFLUX_CANOPY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of canopy plants', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_canopy_mortality_carbonflux_si )

    call this%set_history_var(vname='MORTALITY_CARBONFLUX_UNDERSTORY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of understory plants',use_default='active',&
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_understory_mortality_carbonflux_si )

    ! size class by age dimensioned variables
    call this%set_history_var(vname='NPLANT_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in each size x age class', use_default='active',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scag )

    call this%set_history_var(vname='NPLANT_CANOPY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in canopy in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scag )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in understory in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scag )

    call this%set_history_var(vname='DDBH_CANOPY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of canopy plantsnumber of plants per hectare in canopy in each size x age class', &
          use_default='inactive', avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scag )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of understory plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scag )

    call this%set_history_var(vname='MORTALITY_CANOPY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of canopy plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scag )

    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of understory plantsin each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scag )

    ! size x age x pft dimensioned
    call this%set_history_var(vname='NPLANT_SCAGPFT',units = 'plants/ha',               &
          long='number of plants per hectare in each size x age x pft class', use_default='inactive',   &
          avgflag='A', vtype=site_scagpft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scagpft )

    ! age x pft dimensioned
    call this%set_history_var(vname='NPP_AGEPFT',units = 'kgC/m2/yr',               &
          long='NPP per PFT in each age bin', use_default='inactive',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_si_agepft )

    call this%set_history_var(vname='BIOMASS_AGEPFT',units = 'kg C / m2',               &
          long='biomass per PFT in each age bin', use_default='inactive',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_agepft )

    call this%set_history_var(vname='SCORCH_HEIGHT',units = 'm',               &
          long='SPITFIRE Flame Scorch Height (calculated per PFT in each patch age bin)', &
          use_default='active',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_scorch_height_si_agepft )


    ! Carbon Flux (grid dimension x scpf) (THESE ARE DEFAULT INACTIVE!!!
    !                                     (BECAUSE THEY TAKE UP SPACE!!!
    ! ===================================================================================

    call this%set_history_var(vname='GPP_SCPF', units='kgC/m2/yr',            &
          long='gross primary production by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_si_scpf )

    call this%set_history_var(vname='GPP_CANOPY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of canopy plants by pft/size ', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy_si_scpf )

    call this%set_history_var(vname='AR_CANOPY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of canopy plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy_si_scpf )

    call this%set_history_var(vname='GPP_UNDERSTORY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory_si_scpf )

    call this%set_history_var(vname='AR_UNDERSTORY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_understory_si_scpf )

    call this%set_history_var(vname='NPP_SCPF', units='kgC/m2/yr',            &
          long='total net primary production by pft/size', use_default='inactive',       &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_totl_si_scpf )

    call this%set_history_var(vname='NPP_LEAF_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into leaves by pft/size', use_default='inactive',               &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_si_scpf )

   call this%set_history_var(vname='NPP_SEED_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into seeds by pft/size', use_default='inactive',                &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_si_scpf )

   call this%set_history_var(vname='NPP_FNRT_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into fine roots by pft/size', use_default='inactive',           &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_si_scpf )

   call this%set_history_var(vname='NPP_BGSW_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into below-ground sapwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgsw_si_scpf )

   call this%set_history_var(vname='NPP_BGDW_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into below-ground deadwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgdw_si_scpf )

   call this%set_history_var(vname='NPP_AGSW_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into above-ground sapwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agsw_si_scpf )

   call this%set_history_var(vname = 'NPP_AGDW_SCPF', units='kgC/m2/yr',    &
         long='NPP flux into above-ground deadwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agdw_si_scpf )

   call this%set_history_var(vname = 'NPP_STOR_SCPF', units='kgC/m2/yr',    &
         long='NPP flux into storage by pft/size', use_default='inactive',              &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_si_scpf )

    call this%set_history_var(vname='DDBH_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_si_scpf )

    call this%set_history_var(vname='GROWTHFLUX_SCPF', units = 'n/yr/ha',         &
          long='flux of individuals into a given size class bin via growth and recruitment',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_growthflux_si_scpf )

    call this%set_history_var(vname='GROWTHFLUX_FUSION_SCPF', units = 'n/yr/ha',         &
          long='flux of individuals into a given size class bin via fusion',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_growthflux_fusion_si_scpf )

    call this%set_history_var(vname='DDBH_CANOPY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scpf )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scpf )

    call this%set_history_var(vname='BA_SCPF', units = 'm2/ha',               &
          long='basal area by pft/size', use_default='inactive',   &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scpf )

    call this%set_history_var(vname='AGB_SCPF', units = 'kgC/m2', &
         long='Aboveground biomass by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_agb_si_scpf ) 

    call this%set_history_var(vname='NPLANT_SCPF', units = 'N/ha',         &
          long='stem number density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scpf )

    call this%set_history_var(vname='NPLANT_CAPF', units = 'N/ha',       &
         long='stem number density by pft/coage', use_default='inactive', &
         avgflag='A', vtype=site_coage_pft_r8, hlms='CLM:ALM',flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_capf )

    call this%set_history_var(vname='M1_SCPF', units = 'N/ha/yr',          &
          long='background mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m1_si_scpf )
    
    call this%set_history_var(vname='M2_SCPF', units = 'N/ha/yr',          &
          long='hydraulic mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m2_si_scpf )

    call this%set_history_var(vname='M3_SCPF', units = 'N/ha/yr',          &
          long='carbon starvation mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m3_si_scpf )

    call this%set_history_var(vname='M4_SCPF', units = 'N/ha/yr',          &
          long='impact mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m4_si_scpf )

    call this%set_history_var(vname='M5_SCPF', units = 'N/ha/yr',          &
          long='fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m5_si_scpf )

    call this%set_history_var(vname='CROWNFIREMORT_SCPF', units = 'N/ha/yr',          &
          long='crown fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crownfiremort_si_scpf )

    call this%set_history_var(vname='CAMBIALFIREMORT_SCPF', units = 'N/ha/yr',          &
          long='cambial fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cambialfiremort_si_scpf )

    call this%set_history_var(vname='M6_SCPF', units = 'N/ha/yr',          &
          long='termination mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m6_si_scpf )

    call this%set_history_var(vname='M7_SCPF', units = 'N/ha/event',               &
          long='logging mortality by pft/size',use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m7_si_scpf )

    call this%set_history_var(vname='M8_SCPF', units = 'N/ha/yr',          &
          long='freezing mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m8_si_scpf )

    call this%set_history_var(vname='M9_SCPF', units = 'N/ha/yr',          &
          long='senescence mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m9_si_scpf )

    call this%set_history_var(vname='M10_SCPF', units = 'N/ha/yr',         &
         long='age senescence mortality by pft/size',use_default='inactive', &
         avgflag='A', vtype =site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_si_scpf )
    
    call this%set_history_var(vname='M10_CAPF',units='N/ha/yr',         &
         long='age senescence mortality by pft/cohort age',use_default='inactive', &
         avgflag='A', vtype =site_coage_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,         &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index =ih_m10_si_capf )

    call this%set_history_var(vname='MORTALITY_CANOPY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scpf )

    call this%set_history_var(vname='C13disc_SCPF', units = 'per mil',               &
         long='C13 discrimination by pft/size',use_default='inactive',           &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_c13disc_si_scpf ) 

    call this%set_history_var(vname='BSTOR_CANOPY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_canopy_si_scpf )

    call this%set_history_var(vname='BLEAF_CANOPY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_canopy_si_scpf )

    call this%set_history_var(vname='NPLANT_CANOPY_SCPF', units = 'N/ha',         &
          long='stem number of canopy plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scpf )

    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scpf )

    call this%set_history_var(vname='BSTOR_UNDERSTORY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_understory_si_scpf )

    call this%set_history_var(vname='BLEAF_UNDERSTORY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_understory_si_scpf )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCPF', units = 'N/ha',         &
          long='stem number of understory plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scpf )

    call this%set_history_var(vname='CWD_AG_CWDSC', units='gC/m^2', &
          long='size-resolved AG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_CWDSC', units='gC/m^2', &
          long='size-resolved BG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_si_cwdsc )

    call this%set_history_var(vname='CWD_AG_IN_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_in_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_IN_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_in_si_cwdsc )

    call this%set_history_var(vname='CWD_AG_OUT_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_out_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_OUT_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_out_si_cwdsc )

    ! Size structured diagnostics that require rapid updates (upfreq=2)

    call this%set_history_var(vname='AR_SCPF',units = 'kgC/m2/yr',          &
          long='total autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_si_scpf )
    
    call this%set_history_var(vname='AR_GROW_SCPF',units = 'kgC/m2/yr',          &
          long='growth autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_grow_si_scpf )

    call this%set_history_var(vname='AR_MAINT_SCPF',units = 'kgC/m2/yr',          &
          long='maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_maint_si_scpf )

    call this%set_history_var(vname='AR_DARKM_SCPF',units = 'kgC/m2/yr',          &
          long='dark portion of maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_darkm_si_scpf )

    call this%set_history_var(vname='AR_AGSAPM_SCPF',units = 'kgC/m2/yr',          &
          long='above-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_agsapm_si_scpf )
    
    call this%set_history_var(vname='AR_CROOTM_SCPF',units = 'kgC/m2/yr',          &
          long='below-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_crootm_si_scpf )

    call this%set_history_var(vname='AR_FROOTM_SCPF',units = 'kgC/m2/yr',          &
          long='fine root maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_frootm_si_scpf )

    ! size-class only variables

    call this%set_history_var(vname='DDBH_CANOPY_SCLS', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scls )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCLS', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scls )

    call this%set_history_var(vname='YESTERDAYCANLEV_CANOPY_SCLS', units = 'indiv/ha',               &
          long='Yesterdays canopy level for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_canopy_si_scls )

    call this%set_history_var(vname='YESTERDAYCANLEV_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='Yesterdays canopy level for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_understory_si_scls )

    call this%set_history_var(vname='BA_SCLS', units = 'm2/ha',               &
          long='basal area by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scls )

    call this%set_history_var(vname='AGB_SCLS', units = 'kgC/m2',               &
          long='Aboveground biomass by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_agb_si_scls )

    call this%set_history_var(vname='BIOMASS_SCLS', units = 'kgC/m2',               &
          long='Total biomass by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_scls )

    call this%set_history_var(vname='DEMOTION_RATE_SCLS', units = 'indiv/ha/yr',               &
          long='demotion rate from canopy to understory by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_rate_si_scls )

    call this%set_history_var(vname='PROMOTION_RATE_SCLS', units = 'indiv/ha/yr',               &
          long='promotion rate from understory to canopy by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_rate_si_scls )

    call this%set_history_var(vname='NPLANT_CANOPY_SCLS', units = 'indiv/ha',               &
          long='number of canopy plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scls )
  
    call this%set_history_var(vname='LAI_CANOPY_SCLS', units = 'm2/m2',               &
          long='Leaf are index (LAI) by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_lai_canopy_si_scls )

    call this%set_history_var(vname='SAI_CANOPY_SCLS', units = 'm2/m2',               &
          long='stem area index(SAI) by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sai_canopy_si_scls )

    call this%set_history_var(vname='MORTALITY_CANOPY_SCLS', units = 'indiv/ha/yr',               &
          long='total mortality of canopy trees by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scls )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scls )

    call this%set_history_var(vname='LAI_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_lai_understory_si_scls )

    call this%set_history_var(vname='SAI_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sai_understory_si_scls )

    call this%set_history_var(vname='NPLANT_SCLS', units = 'indiv/ha',               &
          long='number of plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scls )

    call this%set_history_var(vname='NPLANT_CACLS', units = 'indiv/ha',          &
         long='number of plants by coage class', use_default='active',   &
         avgflag='A', vtype=site_coage_r8, hlms='CLM:ALM', flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_cacls )

    call this%set_history_var(vname='M1_SCLS', units = 'N/ha/yr',          &
          long='background mortality by size', use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m1_si_scls )
    
    call this%set_history_var(vname='M2_SCLS', units = 'N/ha/yr',          &
          long='hydraulic mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m2_si_scls )

    call this%set_history_var(vname='M3_SCLS', units = 'N/ha/yr',          &
          long='carbon starvation mortality by size', use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m3_si_scls )

    call this%set_history_var(vname='M4_SCLS', units = 'N/ha/yr',          &
          long='impact mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m4_si_scls )

    call this%set_history_var(vname='M5_SCLS', units = 'N/ha/yr',          &
          long='fire mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m5_si_scls )

    call this%set_history_var(vname='M6_SCLS', units = 'N/ha/yr',          &
          long='termination mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m6_si_scls )

    call this%set_history_var(vname='M7_SCLS', units = 'N/ha/event',               &
          long='logging mortality by size',use_default='active',           &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m7_si_scls )

    call this%set_history_var(vname='M8_SCLS', units = 'N/ha/event',               &
          long='freezing mortality by size',use_default='active',           &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m8_si_scls )

    call this%set_history_var(vname='M9_SCLS', units = 'N/ha/yr',              &
          long='senescence mortality by size',use_default='active',         &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m9_si_scls )

    call this%set_history_var(vname='M10_SCLS', units = 'N/ha/yr',              &
          long='age senescence mortality by size',use_default='active',         &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,     &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_si_scls ) 

    call this%set_history_var(vname='M10_CACLS', units = 'N/ha/yr',             &
          long='age senescence mortality by cohort age',use_default='active',      &
          avgflag='A', vtype=site_coage_r8, hlms='CLM:ALM', flushval=0.0_r8,     &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_si_cacls )

    call this%set_history_var(vname='CARBON_BALANCE_CANOPY_SCLS', units = 'kg C / ha / yr', &
          long='CARBON_BALANCE for canopy plants by size class', use_default='inactive',    &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_canopy_si_scls )

    call this%set_history_var(vname='CARBON_BALANCE_UNDERSTORY_SCLS', units = 'kg C / ha / yr', &
          long='CARBON_BALANCE for understory plants by size class', use_default='inactive',    &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_understory_si_scls )
    
    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCLS', units = 'indiv/ha/yr',               &
          long='total mortality of understory trees by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scls )

    call this%set_history_var(vname='TRIMMING_CANOPY_SCLS', units = 'indiv/ha',               &
          long='trimming term of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_canopy_si_scls )

    call this%set_history_var(vname='TRIMMING_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='trimming term of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_understory_si_scls )

    call this%set_history_var(vname='CROWN_AREA_CANOPY_SCLS', units = 'm2/ha',               &
          long='total crown area of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_canopy_si_scls )

    call this%set_history_var(vname='CROWN_AREA_UNDERSTORY_SCLS', units = 'm2/ha',               &
          long='total crown area of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_understory_si_scls )

    call this%set_history_var(vname='LEAF_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LEAF_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_canopy_si_scls )
    
    call this%set_history_var(vname='ROOT_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='ROOT_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_canopy_si_scls )

    call this%set_history_var(vname='BSTORE_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='BSTORE_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstore_md_canopy_si_scls )

    call this%set_history_var(vname='BDEAD_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='BDEAD_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bdead_md_canopy_si_scls )

    call this%set_history_var(vname='BSW_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='BSW_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bsw_md_canopy_si_scls )

    call this%set_history_var(vname='SEED_PROD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='SEED_PROD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_LEAF_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_LEAF for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=-999.9_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_FROOT_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_FROOT for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_BSW_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BSW for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_sapw_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_BDEAD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BDEAD for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_dead_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_BSEED_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BSEED for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_canopy_si_scls )
    
   call this%set_history_var(vname='NPP_STORE_CANOPY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_STORE for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_canopy_si_scls )
    
    call this%set_history_var(vname='LEAF_MR', units = 'kg C / m2 / yr',               &
          long='RDARK (leaf maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_leaf_mr_si )
    
    call this%set_history_var(vname='FROOT_MR', units = 'kg C / m2 / yr',               &
          long='fine root maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_si )
    
    call this%set_history_var(vname='LIVECROOT_MR', units = 'kg C / m2 / yr',               &
          long='live coarse root maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_si )
    
    call this%set_history_var(vname='LIVESTEM_MR', units = 'kg C / m2 / yr',               &
          long='live stem maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_si )
    
    call this%set_history_var(vname='RDARK_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RDARK for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_canopy_si_scls )
    
    call this%set_history_var(vname='LIVESTEM_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_canopy_si_scls )
    
    call this%set_history_var(vname='LIVECROOT_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_canopy_si_scls )
    
    call this%set_history_var(vname='FROOT_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='FROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_canopy_si_scls )
    
    call this%set_history_var(vname='RESP_G_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_G for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_canopy_si_scls )
    
    call this%set_history_var(vname='RESP_M_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_M for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_canopy_si_scls )

    call this%set_history_var(vname='LEAF_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LEAF_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_understory_si_scls )
    
    call this%set_history_var(vname='ROOT_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='ROOT_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_understory_si_scls )

    call this%set_history_var(vname='BSTORE_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='BSTORE_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstore_md_understory_si_scls )
    
    call this%set_history_var(vname='BDEAD_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='BDEAD_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bdead_md_understory_si_scls )

    call this%set_history_var(vname='BSW_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='BSW_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bsw_md_understory_si_scls )
    
    call this%set_history_var(vname='SEED_PROD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='SEED_PROD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_understory_si_scls )
    
   call this%set_history_var(vname='NPP_LEAF_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_LEAF for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_understory_si_scls )
    
   call this%set_history_var(vname='NPP_FROOT_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_FROOT for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_understory_si_scls )
    
   call this%set_history_var(vname='NPP_BSW_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BSW for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_sapw_understory_si_scls )
    
   call this%set_history_var(vname='NPP_BDEAD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BDEAD for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_dead_understory_si_scls )
    
   call this%set_history_var(vname='NPP_BSEED_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_BSEED for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_understory_si_scls )
    
   call this%set_history_var(vname='NPP_STORE_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
         long='NPP_STORE for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_understory_si_scls )
    
    call this%set_history_var(vname='RDARK_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RDARK for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_understory_si_scls )
    
    call this%set_history_var(vname='LIVESTEM_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_understory_si_scls )
    
    call this%set_history_var(vname='LIVECROOT_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_understory_si_scls )
    
    call this%set_history_var(vname='FROOT_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='FROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_understory_si_scls )
    
    call this%set_history_var(vname='RESP_G_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_G for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_understory_si_scls )
    
    call this%set_history_var(vname='RESP_M_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_M for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_understory_si_scls )


    ! CARBON BALANCE VARIABLES THAT DEPEND ON HLM BGC INPUTS

    call this%set_history_var(vname='NEP', units='gC/m^2/s', &
          long='net ecosystem production', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_nep_si )

    call this%set_history_var(vname='Fire_Closs', units='gC/m^2/s', &
          long='ED/SPitfire Carbon loss to atmosphere', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fire_c_to_atm_si )
   
    call this%set_history_var(vname='FIRE_FLUX', units='g/m^2/s', &
          long='ED-spitfire loss to atmosphere of elements', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_burn_flux_elem )
   
    call this%set_history_var(vname='CBALANCE_ERROR_FATES', units='mgC/day',  &
         long='total carbon error, FATES', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cbal_err_fates_si )

    call this%set_history_var(vname='ERROR_FATES', units='mg/day',  &
         long='total error, FATES mass-balance', use_default='active', &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_err_fates_si )

    call this%set_history_var(vname='LITTER_FINES_AG_ELEM', units='kg/m^2', &
          long='mass of above ground  litter in fines (leaves,nonviable seed)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fines_ag_elem )

    call this%set_history_var(vname='LITTER_FINES_BG_ELEM', units='kg/m^2', &
          long='mass of below ground litter in fines (fineroots)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fines_bg_elem )

    call this%set_history_var(vname='LITTER_CWD_BG_ELEM', units='kg/m^2', &
          long='mass of below ground litter in CWD (coarse roots)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_elem )

    call this%set_history_var(vname='LITTER_CWD_AG_ELEM', units='kg/m^2', &
          long='mass of above ground litter in CWD (trunks/branches/twigs)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_elem )

    call this%set_history_var(vname='LITTER_CWD', units='kg/m^2', &
          long='total mass of litter in CWD', use_default='active', &
          avgflag='A', vtype=site_elcwd_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_elcwd )

    ! organ-partitioned NPP / allocation fluxes
    call this%set_history_var(vname='NPP_LEAF', units='kgC/m2/yr',       &
          long='NPP flux into leaves', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_si )

    call this%set_history_var(vname='NPP_SEED', units='kgC/m2/yr',       &
          long='NPP flux into seeds', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_si )

    call this%set_history_var(vname='NPP_STEM', units='kgC/m2/yr',       &
          long='NPP flux into stem', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stem_si )

    call this%set_history_var(vname='NPP_FROOT', units='kgC/m2/yr',       &
          long='NPP flux into fine roots', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_froot_si )

    call this%set_history_var(vname='NPP_CROOT', units='kgC/m2/yr',       &
          long='NPP flux into coarse roots', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_croot_si )

    call this%set_history_var(vname='NPP_STOR', units='kgC/m2/yr',       &
          long='NPP flux into storage tissues', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_si )


    ! PLANT HYDRAULICS

    if(hlm_use_planthydro.eq.itrue) then
       
       call this%set_history_var(vname='FATES_ERRH2O_SCPF', units='kg/indiv/s', &
             long='mean individual water balance error', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_errh2o_scpf )

       call this%set_history_var(vname='FATES_TRAN_SCPF', units='kg/indiv/s', &
             long='mean individual transpiration rate', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_tran_scpf )

       call this%set_history_var(vname='FATES_SAPFLOW_SCPF', units='kg/ha/s', &
             long='areal sap flow rate dimensioned by size x pft', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sapflow_scpf )

       call this%set_history_var(vname='FATES_SAPFLOW_SI', units='kg/ha/s', &
             long='areal sap flow rate', use_default='active', &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sapflow_si )

       
       call this%set_history_var(vname='FATES_ITERH1_SCPF', units='count/indiv/step', &
             long='number of outer iterations required to achieve tolerable water balance error', &
             use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh1_scpf )
       
       call this%set_history_var(vname='FATES_ITERH2_SCPF', units='count/indiv/step', &
             long='number of inner iterations required to achieve tolerable water balance error', &
             use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh2_scpf )
       
       call this%set_history_var(vname='FATES_ATH_SCPF', units='m3 m-3', &
             long='absorbing root water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_ath_scpf )
       
       call this%set_history_var(vname='FATES_TTH_SCPF', units='m3 m-3', &
             long='transporting root water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index =  ih_tth_scpf )
       
       call this%set_history_var(vname='FATES_STH_SCPF', units='m3 m-3', &
             long='stem water contenet', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sth_scpf )
       
       call this%set_history_var(vname='FATES_LTH_SCPF', units='m3 m-3', &
             long='leaf water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lth_scpf )

       call this%set_history_var(vname='FATES_AWP_SCPF', units='MPa', &
             long='absorbing root water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_awp_scpf )
       
       call this%set_history_var(vname='FATES_TWP_SCPF', units='MPa', &
             long='transporting root water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_twp_scpf )
       
       call this%set_history_var(vname='FATES_SWP_SCPF', units='MPa', &
             long='stem water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_swp_scpf )
       
       call this%set_history_var(vname='FATES_LWP_SCPF', units='MPa', &
             long='leaf water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lwp_scpf )
 
       call this%set_history_var(vname='FATES_AFLC_SCPF', units='fraction', &
             long='absorbing root fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_aflc_scpf )
       
       call this%set_history_var(vname='FATES_TFLC_SCPF', units='fraction', &
             long='transporting root fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_tflc_scpf )
       
       call this%set_history_var(vname='FATES_SFLC_SCPF', units='fraction', &
             long='stem water fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sflc_scpf )
       
       call this%set_history_var(vname='FATES_LFLC_SCPF', units='fraction', &
             long='leaf fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lflc_scpf )
       
       call this%set_history_var(vname='FATES_BTRAN_SCPF', units='unitless', &
             long='mean individual level btran', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_btran_scpf )
       
       call this%set_history_var(vname='FATES_ROOTWGT_SOILVWC_SI', units='m3 m-3', &
            long='soil volumetric water content, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilvwc_si )

       call this%set_history_var(vname='FATES_ROOTWGT_SOILVWCSAT_SI', units='m3 m-3', &
            long='soil saturated volumetric water content, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilvwcsat_si )
       
       call this%set_history_var(vname='FATES_ROOTWGT_SOILMATPOT_SI', units='MPa', &
            long='soil matric potential, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilmatpot_si )
       
       call this%set_history_var(vname='FATES_SOILMATPOT_SL', units='MPa', &
            long='soil water matric potenial by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilmatpot_sl )
       
       call this%set_history_var(vname='FATES_SOILVWC_SL', units='m3 m-3', &
            long='soil volumetric water content by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilvwc_sl )
       
       call this%set_history_var(vname='FATES_SOILVWCSAT_SL', units='m3 m-3', &
            long='soil saturated volumetric water content by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilvwcsat_sl )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE_SI', units='kg ha-1 s-1', &
            long='root water uptake rate', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake_si )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE_SL', units='kg ha-1 s-1', &
            long='root water uptake rate by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake_sl )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE0_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 0 to to 10 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake0_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE10_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 10 to to 50 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake10_scpf )

       call this%set_history_var(vname='FATES_ROOTUPTAKE50_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 50 to to 100 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake50_scpf )

       call this%set_history_var(vname='FATES_ROOTUPTAKE100_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake below 100 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake100_scpf )

       call this%set_history_var(vname='H2OVEG', units = 'kg/m2',               &
             long='water stored inside vegetation tissues (leaf, stem, roots)', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_si )

       call this%set_history_var(vname='H2OVEG_DEAD', units = 'kg/m2',               &
             long='cumulative plant_stored_h2o in dead biomass due to mortality', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_dead_si )

       call this%set_history_var(vname='H2OVEG_RECRUIT', units = 'kg/m2',               &
             long='amount of water in new recruits', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_recruit_si )
    
       call this%set_history_var(vname='H2OVEG_GROWTURN_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) or lost (-) from plant_stored_h2o due to combined growth & turnover', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_growturn_err_si )
    
       call this%set_history_var(vname='H2OVEG_PHENO_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) from plant_stored_h2o due to leaf emergence', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_pheno_err_si )
     
       call this%set_history_var(vname='H2OVEG_HYDRO_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) from plant_stored_h2o due to plant hydrodynamics', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_hydro_err_si )
    end if

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
