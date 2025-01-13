module EDParamsMod

   !
   ! module that deals with reading the ED parameter file
   !

   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : nearzero
   use FatesParametersInterface, only : param_string_length
   use FatesGlobals        , only : fates_log
   use FatesGlobals        , only : endrun => fates_endrun
   use FatesConstantsMod,    only : fates_unset_r8
   use FatesConstantsMod,    only : cstarvation_model_lin
   use FatesConstantsMod,    only : n_landuse_cats

   ! CIME Globals
   use shr_log_mod         , only : errMsg => shr_log_errMsg

   implicit none
   private
   save

   !
   ! this is what the user can use for the actual values
   !

   real(r8),protected, public :: vai_top_bin_width           ! width in VAI units of uppermost leaf+stem
                                                             ! layer scattering element in each canopy layer [m2/m2]
   real(r8),protected, public :: vai_width_increase_factor   ! factor by which each leaf+stem scattering element
                                                             ! increases in VAI width (1 = uniform spacing)
   real(r8),protected, public :: photo_temp_acclim_timescale ! Length of the window for the exponential moving average (ema)
                                                             ! of vegetation temperature used in photosynthesis and respiration
                                                             ! temperature acclimation [days]
   real(r8),protected, public :: photo_temp_acclim_thome_time ! Length of the window for the long-term exponential moving average (ema)
                                                              ! of vegetation temperature used in photosynthesis 
                                                              ! T_home term in Kumarathunge parameterization [years]
   integer,protected, public :: maintresp_leaf_model  ! switch for choosing between leaf maintenance
                                                      ! respiration model. 1=Ryan (1991), 2=Atkin et al (2017)
   real(r8),protected, public :: sdlng_emerg_h2o_timescale !Length of the window for the exponential moving
                                                                 !average of smp used to calculate seedling emergence
   real(r8),protected, public :: sdlng_mort_par_timescale !Length of the window for the exponential moving average 
                                                                !of par at the seedling layer used to calculate 
                                                                !seedling mortality
   real(r8),protected, public :: sdlng_mdd_timescale !Length of the window for the exponential moving average
                                                           ! of moisture deficit days used to calculate seedling mortality
   real(r8),protected, public :: sdlng2sap_par_timescale !Length of the window for the exponential 
                                                               !moving average of par at the seedling layer used to 
                                                               !calculate seedling to sapling transition rates
   integer,protected, public :: photo_tempsens_model  ! switch for choosing the model that defines the temperature
                                                      ! sensitivity of photosynthetic parameters (vcmax, jmax).
                                                      ! 1=non-acclimating, 2=Kumarathunge et al., 2019

   integer,protected, public :: radiation_model       ! Switch betrween Norman (1) and Two-stream (2) radiation models

   integer,protected, public :: mort_cstarvation_model ! Switch for carbon starvation mortality:
                                                       ! 1 -- Linear model
                                                       ! 2 -- Exponential model

   real(r8),protected, public :: fates_mortality_disturbance_fraction ! the fraction of canopy mortality that results in disturbance
   real(r8),protected, public :: ED_val_comp_excln                    ! weighting factor for canopy layer exclusion and promotion
   real(r8),protected, public :: ED_val_vai_top_bin_width             ! width in VAI units of uppermost leaf+stem layer scattering element
   real(r8),protected, public :: ED_val_vai_width_increase_factor     ! factor by which each leaf+stem scattering element increases in VAI width
   real(r8),protected, public :: ED_val_nignitions                    ! number of annual ignitions per square km
   real(r8),protected, public :: ED_val_understorey_death             ! fraction of plants in understorey cohort impacted by overstorey tree-fall
   real(r8),protected, public :: ED_val_cwd_fcel                      ! Cellulose fraction for CWD
   real(r8),protected, public :: ED_val_cwd_flig                      ! Lignin fraction of coarse woody debris
   real(r8),protected, public :: maintresp_nonleaf_baserate           ! Base maintenance respiration rate for plant tissues
   real(r8),protected, public :: ED_val_phen_a                        ! GDD accumulation function, intercept parameter: gdd_thesh = a + b exp(c*ncd)
   real(r8),protected, public :: ED_val_phen_b                        ! GDD accumulation function, multiplier parameter: gdd_thesh = a + b exp(c*ncd)
   real(r8),protected, public :: ED_val_phen_c                        ! GDD accumulation function, exponent parameter: gdd_thesh = a + b exp(c*ncd)
   real(r8),protected, public :: ED_val_phen_chiltemp                 ! chilling day counting threshold for vegetation
   real(r8),protected, public :: ED_val_phen_mindayson                ! day threshold compared against days since leaves became on-allometry
   real(r8),protected, public :: ED_val_phen_ncolddayslim             ! day threshold exceedance for temperature leaf-drop
   real(r8),protected, public :: ED_val_phen_coldtemp                 ! vegetation temperature exceedance that flags a cold-day for leaf-drop
   real(r8),protected, public :: ED_val_cohort_size_fusion_tol        ! minimum fraction in difference in dbh between cohorts
   real(r8),protected, public :: ED_val_cohort_age_fusion_tol         ! minimum fraction in differece in cohort age between cohorts
   real(r8),protected, public :: ED_val_patch_fusion_tol              ! minimum fraction in difference in profiles between patches
   real(r8),protected, public :: ED_val_canopy_closure_thresh         ! site-level canopy closure point where trees take on forest (narrow) versus savannah (wide) crown allometry
   integer,protected, public  :: stomatal_model                       ! switch for choosing between stomatal conductance models, 1 for Ball-Berry, 2 for Medlyn
   integer,protected, public  :: dayl_switch                          ! switch for turning on or off day length factor scaling for photosynthetic parameters
   integer,protected, public  :: regeneration_model                   ! Switch for choosing between regeneration models:
                                                                      ! (1) for Fates default
                                                                      ! (2) for the Tree Recruitment Scheme (Hanbury-Brown et al., 2022)
                                                                      ! (3) for the Tree Recruitment Scheme without seedling dynamics
   
   
   logical,protected, public :: active_crown_fire        ! flag, 1=active crown fire 0=no active crown fire
   character(len=param_string_length),parameter :: fates_name_active_crown_fire = "fates_fire_active_crown_fire"

   real(r8), protected, public :: cg_strikes             ! fraction of cloud to ground lightning strikes (0-1)
   character(len=param_string_length),parameter :: fates_name_cg_strikes="fates_fire_cg_strikes"

   ! empirical curvature parameters for ac, aj photosynthesis co-limitation, c3 and c4 plants respectively
   real(r8),protected,public  :: theta_cj_c3    ! Empirical curvature parameter for ac, aj photosynthesis co-limitation in c3 plants
   real(r8),protected,public  :: theta_cj_c4    ! Empirical curvature parameter for ac, aj photosynthesis co-limitation in c4 plants

     ! Global identifier of how nutrients interact with the host land model
  ! either they are fully coupled, or they generate uptake rates synthetically
  ! in prescribed mode. In the latter, there is both NO mass removed from the HLM's soil
  ! BGC N and P pools, and there is also none removed.

   integer, public :: n_uptake_mode
   integer, public :: p_uptake_mode

   real(r8), parameter, public :: soil_tfrz_thresh = -2.0_r8 ! Soil temperature threshold below which hydraulic failure mortality is off (non-hydro only) in degrees C
   
   integer, parameter, public :: nclmax = 2   ! Maximum number of canopy layers (used only for scratch arrays)
                                              ! We would make this even higher, but making this
                                              ! a little lower keeps the size down on some output arrays
                                              ! For large arrays at patch level we use dynamic allocation

   ! parameters that govern the VAI (LAI+SAI) bins used in radiative transfer code
   integer, parameter, public :: nlevleaf = 30   ! number of leaf+stem layers in each canopy layer

   real(r8), public :: dinc_vai(nlevleaf)   = fates_unset_r8 ! VAI bin widths array
   real(r8), public :: dlower_vai(nlevleaf) = fates_unset_r8 ! lower edges of VAI bins
 
   integer, parameter, public :: maxpft = 16      ! maximum number of PFTs allowed
   
   real(r8),protected,public  :: q10_mr     ! Q10 for respiration rate (for soil fragmenation and plant respiration)    (unitless)
   real(r8),protected,public  :: q10_froz   ! Q10 for frozen-soil respiration rates (for soil fragmentation)            (unitless)

   ! grazing parameters
   real(r8),protected,public :: landuse_grazing_carbon_use_eff
   real(r8),protected,public :: landuse_grazing_maxheight
   real(r8),protected,public :: landuse_grazing_nitrogen_use_eff
   real(r8),protected,public :: landuse_grazing_phosphorus_use_eff
   real(r8),protected,public :: landuse_grazing_rate(n_landuse_cats)

   ! Unassociated pft dimensioned free parameter that developers can use for testing arbitrary new hypotheses
   ! (THIS PARAMETER IS UNUSED, FEEL FREE TO USE IT FOR WHATEVER PURPOSE YOU LIKE. WE CAN
   !  HELP MIGRATE YOUR USAGE OF THE PARMETER TO A PERMANENT HOME LATER)
   real(r8),protected,public  :: dev_arbitrary  ! Unassociated free parameter that developers can use for testing arbitrary new hypotheses
   character(len=param_string_length),parameter,public :: name_dev_arbitrary = "fates_dev_arbitrary"
   
   ! parameters whose size is defined in the parameter file
   real(r8),protected,allocatable,public :: ED_val_history_sizeclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_ageclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_height_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_coageclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_damage_bin_edges(:)
   
   ! Switch that defines the current pressure-volume and pressure-conductivity model
   ! to be used at each node (compartment/organ)
   ! 1  = Christofferson et al. 2016 (TFS),   2 = Van Genuchten 1980
   
   character(len=param_string_length),parameter,public :: ED_name_sdlng_emerg_h2o_timescale = "fates_trs_seedling_emerg_h2o_timescale"
   character(len=param_string_length),parameter,public :: ED_name_sdlng_mort_par_timescale = "fates_trs_seedling_mort_par_timescale"
   character(len=param_string_length),parameter,public :: ED_name_sdlng_mdd_timescale = "fates_trs_seedling_mdd_timescale"
   character(len=param_string_length),parameter,public :: ED_name_sdlng2sap_par_timescale = "fates_trs_seedling2sap_par_timescale"
   integer, protected,allocatable,public :: hydr_htftype_node(:)
   character(len=param_string_length),parameter,public :: ED_name_photo_temp_acclim_timescale = "fates_leaf_photo_temp_acclim_timescale"
   character(len=param_string_length),parameter,public :: ED_name_photo_temp_acclim_thome_time = "fates_leaf_photo_temp_acclim_thome_time"
   character(len=param_string_length),parameter,public :: name_photo_tempsens_model = "fates_leaf_photo_tempsens_model"
   character(len=param_string_length),parameter,public :: name_maintresp_model = "fates_maintresp_leaf_model"
   character(len=param_string_length),parameter,public :: name_radiation_model = "fates_rad_model"
   character(len=param_string_length),parameter,public :: ED_name_hydr_htftype_node = "fates_hydro_htftype_node"
   character(len=param_string_length),parameter,public :: ED_name_mort_disturb_frac = "fates_mort_disturb_frac"
   character(len=param_string_length),parameter,public :: ED_name_mort_cstarvation_model = "fates_mort_cstarvation_model"
   character(len=param_string_length),parameter,public :: ED_name_comp_excln = "fates_comp_excln"
   character(len=param_string_length),parameter,public :: ED_name_vai_top_bin_width = "fates_vai_top_bin_width"
   character(len=param_string_length),parameter,public :: ED_name_vai_width_increase_factor = "fates_vai_width_increase_factor"
   character(len=param_string_length),parameter,public :: ED_name_nignitions = "fates_fire_nignitions"
   character(len=param_string_length),parameter,public :: ED_name_understorey_death = "fates_mort_understorey_death"
   character(len=param_string_length),parameter,public :: ED_name_cwd_fcel= "fates_frag_cwd_fcel"   
   character(len=param_string_length),parameter,public :: ED_name_cwd_flig= "fates_frag_cwd_flig"   
   character(len=param_string_length),parameter,public :: fates_name_maintresp_nonleaf_baserate= "fates_maintresp_nonleaf_baserate"
   character(len=param_string_length),parameter,public :: ED_name_phen_a= "fates_phen_gddthresh_a"   
   character(len=param_string_length),parameter,public :: ED_name_phen_b= "fates_phen_gddthresh_b"   
   character(len=param_string_length),parameter,public :: ED_name_phen_c= "fates_phen_gddthresh_c"   
   character(len=param_string_length),parameter,public :: ED_name_phen_chiltemp= "fates_phen_chilltemp"   
   character(len=param_string_length),parameter,public :: ED_name_phen_mindayson= "fates_phen_mindayson"
   character(len=param_string_length),parameter,public :: ED_name_phen_ncolddayslim= "fates_phen_ncolddayslim"   
   character(len=param_string_length),parameter,public :: ED_name_phen_coldtemp= "fates_phen_coldtemp"   
   character(len=param_string_length),parameter,public :: ED_name_cohort_size_fusion_tol= "fates_cohort_size_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_cohort_age_fusion_tol = "fates_cohort_age_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_patch_fusion_tol= "fates_patch_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_canopy_closure_thresh= "fates_canopy_closure_thresh"      
   character(len=param_string_length),parameter,public :: ED_name_stomatal_model= "fates_leaf_stomatal_model"
   character(len=param_string_length),parameter,public :: ED_name_dayl_switch= "fates_daylength_factor_switch"
   character(len=param_string_length),parameter,public :: ED_name_regeneration_model= "fates_regeneration_model"

   character(len=param_string_length),parameter,public :: name_theta_cj_c3 = "fates_leaf_theta_cj_c3"
   character(len=param_string_length),parameter,public :: name_theta_cj_c4 = "fates_leaf_theta_cj_c4"
   
   character(len=param_string_length),parameter :: fates_name_q10_mr="fates_q10_mr"
   character(len=param_string_length),parameter :: fates_name_q10_froz="fates_q10_froz"

   ! non-scalar parameter names
   character(len=param_string_length),parameter,public :: ED_name_history_sizeclass_bin_edges= "fates_history_sizeclass_bin_edges"      
   character(len=param_string_length),parameter,public :: ED_name_history_ageclass_bin_edges= "fates_history_ageclass_bin_edges"      
   character(len=param_string_length),parameter,public :: ED_name_history_height_bin_edges= "fates_history_height_bin_edges"
   character(len=param_string_length),parameter,public :: ED_name_history_coageclass_bin_edges = "fates_history_coageclass_bin_edges"
   character(len=param_string_length),parameter,public :: ED_name_history_damage_bin_edges = "fates_history_damage_bin_edges"
   character(len=param_string_length),parameter,public :: ED_name_crop_lu_pft_vector = "fates_landuse_crop_lu_pft_vector"
   character(len=param_string_length),parameter,public :: ED_name_maxpatches_by_landuse = "fates_maxpatches_by_landuse"
   character(len=param_string_length),parameter,public :: ED_name_max_nocomp_pfts_by_landuse = "fates_max_nocomp_pfts_by_landuse"


   ! Hydraulics Control Parameters (ONLY RELEVANT WHEN USE_FATES_HYDR = TRUE)
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: hydr_kmax_rsurf1         !  maximum conducitivity for unit root surface 
                                                  !  soil to root direction (kg water/m2 root area/Mpa/s)
   character(len=param_string_length),parameter,public :: hydr_name_kmax_rsurf1 = "fates_hydro_kmax_rsurf1"  
   
   real(r8),protected,public :: hydr_kmax_rsurf2         !  maximum conducitivity for unit root surface 
                                                  !  root to soil direciton (kg water/m2 root area/Mpa/s)
   character(len=param_string_length),parameter,public :: hydr_name_kmax_rsurf2 = "fates_hydro_kmax_rsurf2" 

   real(r8),protected,public :: hydr_psi0          !  sapwood water potential at saturation (MPa)
   character(len=param_string_length),parameter,public :: hydr_name_psi0 = "fates_hydro_psi0"

   real(r8),protected,public :: hydr_psicap        !  sapwood water potential at which capillary reserves exhausted (MPa)
   character(len=param_string_length),parameter,public :: hydr_name_psicap = "fates_hydro_psicap"

   
   ! Switch that defines which hydraulic solver to use
   ! 1 = Taylor solution that solves plant fluxes with 1 layer
   !     sequentially placing solution on top of previous layer solves
   ! 2 = Picard solution that solves all fluxes in a plant and
   !     the soil simultaneously, 2D: soil x (root + shell)
   ! 3 = Newton-Raphson (Deprecated) solution that solves all fluxes in a plant and
   !     the soil simultaneously, 2D: soil x (root + shell)
   
   integer,protected,public :: hydr_solver        !  switch designating hydraulics numerical solver
   character(len=param_string_length),parameter,public :: hydr_name_solver = "fates_hydro_solver"
   
   !Soil BGC parameters, mostly used for testing FATES when not coupled to the dynamics bgc hlm
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: bgc_soil_salinity ! site-level soil salinity for FATES when not coupled to dynamic soil BGC of salinity
   character(len=param_string_length),parameter,public :: bgc_name_soil_salinity= "fates_soil_salinity"      

   ! Switch designating whether to use net or gross assimilation in the stomata model
   integer, protected, public :: stomatal_assim_model
   character(len=param_string_length), parameter, public :: stomatal_assim_name = "fates_leaf_stomatal_assim_model"

   ! Integer code that options how damage events are structured
   integer, protected, public :: damage_event_code
   character(len=param_string_length), parameter, public :: damage_name_event_code = "fates_damage_event_code"

   integer,protected,public :: damage_canopy_layer_code  ! Code that changes whether damage affects canopy trees (1), understory trees (2)
   character(len=param_string_length),parameter,public :: damage_name_canopy_layer_code = "fates_damage_canopy_layer_code"
   
   ! Maximum allowable primary and secondary patches
   ! These values are USED FOR ALLOCATIONS IN BOTH FATES AND CLM/ELM!!!!
   ! The number of patches specified in the parameter file may be over-written.
   ! For instance, in SP mode, we want the same number of primary patches as the number of PFTs
   ! in the fates parameter file, and zero secondary.
   ! thus they are not protected here.
   
   integer, public :: maxpatches_by_landuse(n_landuse_cats)
   integer, public :: max_nocomp_pfts_by_landuse(n_landuse_cats)
   integer, public :: maxpatch_total

   ! which crops can be grown on a given crop land use type
   integer,protected,public :: crop_lu_pft_vector(n_landuse_cats)

   ! Maximum allowable cohorts per patch
   integer, protected, public :: max_cohort_per_patch
   character(len=param_string_length), parameter, public :: maxcohort_name = "fates_maxcohort"

   
   
   ! Logging Control Parameters (ONLY RELEVANT WHEN USE_FATES_LOGGING = TRUE)
   ! ----------------------------------------------------------------------------------------------

   real(r8),protected,public :: logging_dbhmin              ! Minimum dbh at which logging is applied (cm)
                                                            ! Typically associated with harvesting
   character(len=param_string_length),parameter,public :: logging_name_dbhmin = "fates_landuse_logging_dbhmin"

   real(r8),protected,public :: logging_dbhmax              ! Maximum dbh at which logging is applied (cm)
                                                            ! Typically associated with fire suppression
   character(len=param_string_length),parameter,public :: logging_name_dbhmax = "fates_landuse_logging_dbhmax"


   real(r8),protected,public :: logging_collateral_frac     ! Ratio of collateral mortality to direct logging mortality
   character(len=param_string_length),parameter,public :: logging_name_collateral_frac = "fates_landuse_logging_collateral_frac"

   real(r8),protected,public :: logging_coll_under_frac ! Fraction of understory plants that die when logging disturbance
                                                 ! is generated
   character(len=param_string_length),parameter,public :: logging_name_coll_under_frac = "fates_landuse_logging_coll_under_frac"
   
   real(r8),protected,public :: logging_direct_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter,public :: logging_name_direct_frac = "fates_landuse_logging_direct_frac"

   real(r8),protected,public :: logging_mechanical_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter,public :: logging_name_mechanical_frac = "fates_landuse_logging_mechanical_frac"

   real(r8),protected,public :: logging_event_code          ! Code that options how logging events are structured 
   character(len=param_string_length),parameter,public :: logging_name_event_code = "fates_landuse_logging_event_code"
   
   real(r8),protected,public :: logging_dbhmax_infra        ! "Tree diameter, above which infrastructure from logging does not impact damage or mortality.
   character(len=param_string_length),parameter,public :: logging_name_dbhmax_infra = "fates_landuse_logging_dbhmax_infra"
   
   real(r8),protected,public :: logging_export_frac        ! "fraction of trunk product being shipped offsite, the 
                                                    ! leftovers will be left onsite as large CWD
   character(len=param_string_length),parameter,public :: logging_name_export_frac ="fates_landuse_logging_export_frac"   

   ! grazing-related parameters
   character(len=param_string_length),parameter,public :: name_landuse_grazing_rate               = "fates_landuse_grazing_rate"
   character(len=param_string_length),parameter,public :: name_landuse_grazing_carbon_use_eff     = "fates_landuse_grazing_carbon_use_eff"
   character(len=param_string_length),parameter,public :: name_landuse_grazing_maxheight          = "fates_landuse_grazing_maxheight"
   character(len=param_string_length),parameter,public :: name_landuse_grazing_nitrogen_use_eff   = "fates_landuse_grazing_nitrogen_use_eff"
   character(len=param_string_length),parameter,public :: name_landuse_grazing_phosphorus_use_eff = "fates_landuse_grazing_phosphorus_use_eff"

   real(r8),protected,public :: eca_plant_escalar  ! scaling factor for plant fine root biomass to 
                                               ! calculate nutrient carrier enzyme abundance (ECA)

   

   
   character(len=param_string_length),parameter,public :: eca_name_plant_escalar = "fates_cnp_eca_plant_escalar"

   public :: FatesParamsInit
   public :: FatesRegisterParams
   public :: FatesReceiveParams
   public :: FatesReportParams
  
contains

  !-----------------------------------------------------------------------
  subroutine FatesParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    vai_top_bin_width                     = nan
    vai_width_increase_factor             = nan
    photo_temp_acclim_timescale           = nan
    sdlng_emerg_h2o_timescale             = nan
    sdlng_mort_par_timescale              = nan
    sdlng_mdd_timescale                   = nan
    sdlng2sap_par_timescale               = nan
    photo_temp_acclim_thome_time          = nan
    photo_tempsens_model                  = -9
    maintresp_leaf_model                  = -9
    radiation_model                       = -9
    fates_mortality_disturbance_fraction  = nan
    mort_cstarvation_model                = -9
    ED_val_comp_excln                     = nan
    ED_val_vai_top_bin_width              = nan
    ED_val_vai_width_increase_factor      = nan
    ED_val_nignitions                     = nan
    ED_val_understorey_death              = nan
    ED_val_cwd_fcel                       = nan
    ED_val_cwd_flig                       = nan
    maintresp_nonleaf_baserate            = nan
    ED_val_phen_a                         = nan
    ED_val_phen_b                         = nan
    ED_val_phen_c                         = nan
    ED_val_phen_chiltemp                  = nan
    ED_val_phen_mindayson                 = nan
    ED_val_phen_ncolddayslim              = nan
    ED_val_phen_coldtemp                  = nan
    ED_val_cohort_size_fusion_tol         = nan
    ED_val_cohort_age_fusion_tol          = nan
    ED_val_patch_fusion_tol               = nan
    ED_val_canopy_closure_thresh          = nan
    stomatal_model                        = -9
    dayl_switch                           = -9
    regeneration_model                    = -9
    stomatal_assim_model                  = -9
    max_cohort_per_patch                  = -9
    hydr_kmax_rsurf1                      = nan
    hydr_kmax_rsurf2                      = nan
    hydr_psi0                             = nan
    hydr_psicap                           = nan
    hydr_solver                           = -9
    bgc_soil_salinity                     = nan
    logging_dbhmin                        = nan
    logging_dbhmax                        = nan
    logging_collateral_frac               = nan
    logging_direct_frac                   = nan
    logging_mechanical_frac               = nan
    logging_event_code                    = nan
    logging_dbhmax_infra                  = nan
    logging_export_frac                   = nan
    eca_plant_escalar                     = nan
    q10_mr                                = nan
    q10_froz                              = nan
    theta_cj_c3                           = nan
    theta_cj_c4                           = nan
    dev_arbitrary                         = nan
    damage_event_code                     = -9
    damage_canopy_layer_code              = -9
    landuse_grazing_carbon_use_eff        = nan
    landuse_grazing_nitrogen_use_eff      = nan
    landuse_grazing_phosphorus_use_eff    = nan
    landuse_grazing_maxheight             = nan
    landuse_grazing_rate(:)               = nan

  end subroutine FatesParamsInit

  !-----------------------------------------------------------------------
  subroutine FatesRegisterParams(fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_1d
    use FatesParametersInterface, only : dimension_name_history_size_bins, dimension_name_history_age_bins
    use FatesParametersInterface, only : dimension_name_history_height_bins, dimension_name_hydr_organs
    use FatesParametersInterface, only : dimension_name_history_coage_bins, dimension_name_history_damage_bins
    use FatesParametersInterface, only : dimension_shape_scalar, dimension_name_landuse


    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)
    character(len=param_string_length), parameter :: dim_names_sizeclass(1) = (/dimension_name_history_size_bins/)
    character(len=param_string_length), parameter :: dim_names_ageclass(1) = (/dimension_name_history_age_bins/)
    character(len=param_string_length), parameter :: dim_names_height(1) = (/dimension_name_history_height_bins/)
    character(len=param_string_length), parameter :: dim_names_coageclass(1) = (/dimension_name_history_coage_bins/)
    character(len=param_string_length), parameter :: dim_names_hydro_organs(1) = (/dimension_name_hydr_organs/)
    character(len=param_string_length), parameter :: dim_names_damageclass(1)= (/dimension_name_history_damage_bins/)
    character(len=param_string_length), parameter :: dim_names_landuse(1)= (/dimension_name_landuse/)
    
    call FatesParamsInit()

    call fates_params%RegisterParameter(name=ED_name_photo_temp_acclim_timescale, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_sdlng_emerg_h2o_timescale, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_sdlng_mort_par_timescale, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_sdlng_mdd_timescale, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_sdlng2sap_par_timescale, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_photo_temp_acclim_thome_time, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_photo_tempsens_model,dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_radiation_model,dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=name_maintresp_model,dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=name_theta_cj_c3, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=name_theta_cj_c4, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=ED_name_mort_disturb_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_mort_cstarvation_model, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_comp_excln, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_vai_top_bin_width, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_vai_width_increase_factor, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_nignitions, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_understorey_death, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cwd_fcel, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cwd_flig, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=fates_name_maintresp_nonleaf_baserate, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_a, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_b, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_c, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_chiltemp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_mindayson, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_ncolddayslim, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_coldtemp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cohort_size_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cohort_age_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_patch_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_canopy_closure_thresh, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_stomatal_model, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_dayl_switch, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_regeneration_model, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
	 
    call fates_params%RegisterParameter(name=stomatal_assim_name, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=maxcohort_name, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=hydr_name_solver, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=hydr_name_kmax_rsurf1, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=hydr_name_kmax_rsurf2, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=hydr_name_psi0, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=hydr_name_psicap, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=bgc_name_soil_salinity, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar) 

    call fates_params%RegisterParameter(name=logging_name_dbhmin, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_dbhmax, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_collateral_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_coll_under_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_direct_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_mechanical_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_event_code, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_dbhmax_infra, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_export_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=eca_name_plant_escalar, dimension_shape=dimension_shape_scalar, & 
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=fates_name_q10_mr, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=fates_name_q10_froz, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_dev_arbitrary, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=damage_name_event_code, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=damage_name_canopy_layer_code, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_landuse_grazing_carbon_use_eff, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_landuse_grazing_maxheight, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_landuse_grazing_nitrogen_use_eff, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=name_landuse_grazing_phosphorus_use_eff, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    ! non-scalar parameters

    call fates_params%RegisterParameter(name=ED_name_hydr_htftype_node, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_hydro_organs)
    
    call fates_params%RegisterParameter(name=ED_name_history_sizeclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_sizeclass)

    call fates_params%RegisterParameter(name=ED_name_history_ageclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_ageclass)

    call fates_params%RegisterParameter(name=ED_name_history_height_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_height)

    call fates_params%RegisterParameter(name=fates_name_active_crown_fire, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=fates_name_cg_strikes, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=ED_name_history_coageclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_coageclass)

    call fates_params%RegisterParameter(name=ED_name_history_damage_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_damageclass)

    call fates_params%RegisterParameter(name=ED_name_crop_lu_pft_vector, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_landuse)

    call fates_params%RegisterParameter(name=ED_name_maxpatches_by_landuse, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_landuse)

    call fates_params%RegisterParameter(name=ED_name_max_nocomp_pfts_by_landuse, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_landuse)

    call fates_params%RegisterParameter(name=name_landuse_grazing_rate, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_landuse)
  end subroutine FatesRegisterParams

  
  !-----------------------------------------------------------------------
  subroutine FatesReceiveParams(fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar
    use FatesConstantsMod, only: primaryland, secondaryland, rangeland, pastureland, cropland

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    real(r8) :: tmpreal ! local real variable for changing type on read
    real(r8), allocatable :: hydr_htftype_real(:)
    real(r8), allocatable :: tmp_vector_by_landuse1(:)  ! local real vector for changing type on read
    real(r8), allocatable :: tmp_vector_by_landuse2(:)  ! local real vector for changing type on read
    real(r8), allocatable :: tmp_vector_by_landuse3(:)  ! local real vector for changing type on read
    real(r8), allocatable :: tmp_vector_by_landuse4(:)  ! local real vector for changing type on read

    call fates_params%RetrieveParameter(name=ED_name_photo_temp_acclim_timescale, &
         data=photo_temp_acclim_timescale)

    call fates_params%RetrieveParameter(name=ED_name_sdlng_emerg_h2o_timescale, &
         data=sdlng_emerg_h2o_timescale)

    call fates_params%RetrieveParameter(name=ED_name_sdlng_mort_par_timescale, &
         data=sdlng_mort_par_timescale)

    call fates_params%RetrieveParameter(name=ED_name_sdlng_mdd_timescale, &
         data=sdlng_mdd_timescale)

    call fates_params%RetrieveParameter(name=ED_name_sdlng2sap_par_timescale, &
         data=sdlng2sap_par_timescale)
    
    call fates_params%RetrieveParameter(name=ED_name_photo_temp_acclim_thome_time, &
         data=photo_temp_acclim_thome_time)

    call fates_params%RetrieveParameter(name=name_photo_tempsens_model, &
         data=tmpreal)
    photo_tempsens_model = nint(tmpreal)

    call fates_params%RetrieveParameter(name=name_radiation_model, &
         data=tmpreal)
    radiation_model = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=name_maintresp_model, &
         data=tmpreal)
    maintresp_leaf_model = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=ED_name_mort_disturb_frac, &
          data=fates_mortality_disturbance_fraction)

    call fates_params%RetrieveParameter(name=ED_name_mort_cstarvation_model, &
         data=tmpreal)
    mort_cstarvation_model = nint(tmpreal)
        
    call fates_params%RetrieveParameter(name=ED_name_comp_excln, &
         data=ED_val_comp_excln)

    call fates_params%RetrieveParameter(name=ED_name_vai_top_bin_width, &
         data=ED_val_vai_top_bin_width)

    call fates_params%RetrieveParameter(name=ED_name_vai_width_increase_factor, &
         data=ED_val_vai_width_increase_factor)

    call fates_params%RetrieveParameter(name=ED_name_nignitions, &
         data=ED_val_nignitions)

    call fates_params%RetrieveParameter(name=ED_name_understorey_death, &
         data=ED_val_understorey_death)

    call fates_params%RetrieveParameter(name=ED_name_cwd_fcel, &
         data=ED_val_cwd_fcel)

    call fates_params%RetrieveParameter(name=ED_name_cwd_flig, &
         data=ED_val_cwd_flig)

    call fates_params%RetrieveParameter(name=fates_name_maintresp_nonleaf_baserate, &
         data=maintresp_nonleaf_baserate)

    call fates_params%RetrieveParameter(name=ED_name_phen_a, &
         data=ED_val_phen_a)

    call fates_params%RetrieveParameter(name=ED_name_phen_b, &
         data=ED_val_phen_b)

    call fates_params%RetrieveParameter(name=ED_name_phen_c, &
         data=ED_val_phen_c)

    call fates_params%RetrieveParameter(name=ED_name_phen_chiltemp, &
         data=ED_val_phen_chiltemp)

    call fates_params%RetrieveParameter(name=ED_name_phen_mindayson, &
         data=ED_val_phen_mindayson)

    call fates_params%RetrieveParameter(name=ED_name_phen_ncolddayslim, &
         data=ED_val_phen_ncolddayslim)

    call fates_params%RetrieveParameter(name=ED_name_phen_coldtemp, &
         data=ED_val_phen_coldtemp)

    call fates_params%RetrieveParameter(name=ED_name_cohort_size_fusion_tol, &
         data=ED_val_cohort_size_fusion_tol)

    call fates_params%RetrieveParameter(name=ED_name_cohort_age_fusion_tol, &
         data=ED_val_cohort_age_fusion_tol)

    call fates_params%RetrieveParameter(name=ED_name_patch_fusion_tol, &
         data=ED_val_patch_fusion_tol)
    
    call fates_params%RetrieveParameter(name=ED_name_canopy_closure_thresh, &
         data=ED_val_canopy_closure_thresh)

    call fates_params%RetrieveParameter(name=ED_name_stomatal_model, &
         data=tmpreal)
    stomatal_model = nint(tmpreal)

    call fates_params%RetrieveParameter(name=ED_name_dayl_switch, &
         data=tmpreal)
    dayl_switch = nint(tmpreal)

    call fates_params%RetrieveParameter(name=ED_name_regeneration_model, &
         data=tmpreal)
    regeneration_model = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=stomatal_assim_name, &
         data=tmpreal)
    stomatal_assim_model = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=maxcohort_name, &
         data=tmpreal)
    max_cohort_per_patch = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=hydr_name_kmax_rsurf1, &
          data=hydr_kmax_rsurf1)

    call fates_params%RetrieveParameter(name=hydr_name_kmax_rsurf2, &
          data=hydr_kmax_rsurf2)	 
    
    call fates_params%RetrieveParameter(name=hydr_name_psi0, &
          data=hydr_psi0)

    call fates_params%RetrieveParameter(name=hydr_name_psicap, &
          data=hydr_psicap)

    call fates_params%RetrieveParameter(name=hydr_name_solver, &
         data=tmpreal)
    hydr_solver = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=bgc_name_soil_salinity, &
          data=bgc_soil_salinity)	  

    call fates_params%RetrieveParameter(name=logging_name_dbhmin, &
          data=logging_dbhmin)

    call fates_params%RetrieveParameter(name=logging_name_dbhmax, &
          data=logging_dbhmax)

    call fates_params%RetrieveParameter(name=logging_name_collateral_frac, &
          data=logging_collateral_frac)
    
    call fates_params%RetrieveParameter(name=logging_name_coll_under_frac, &
          data=logging_coll_under_frac)

    call fates_params%RetrieveParameter(name=logging_name_direct_frac, &
          data=logging_direct_frac)

    call fates_params%RetrieveParameter(name=logging_name_mechanical_frac, &
          data=logging_mechanical_frac)
    
    call fates_params%RetrieveParameter(name=logging_name_event_code, &
          data=logging_event_code)

    call fates_params%RetrieveParameter(name=logging_name_dbhmax_infra, &
          data=logging_dbhmax_infra)

    call fates_params%RetrieveParameter(name=logging_name_export_frac, &
          data=logging_export_frac)

    call fates_params%RetrieveParameter(name=eca_name_plant_escalar, &
          data=eca_plant_escalar)

    call fates_params%RetrieveParameter(name=name_theta_cj_c3, &
          data=theta_cj_c3)

     call fates_params%RetrieveParameter(name=name_theta_cj_c4, &
          data=theta_cj_c4)
     
    call fates_params%RetrieveParameter(name=fates_name_q10_mr, &
          data=q10_mr)
    
    call fates_params%RetrieveParameter(name=fates_name_q10_froz, &
         data=q10_froz)
    
    call fates_params%RetrieveParameter(name=name_dev_arbitrary, &
         data=dev_arbitrary)

    call fates_params%RetrieveParameter(name=fates_name_active_crown_fire, & 
          data=tmpreal)
    active_crown_fire = (abs(tmpreal-1.0_r8)<nearzero)

    call fates_params%RetrieveParameter(name=fates_name_cg_strikes, &
          data=cg_strikes)

    call fates_params%RetrieveParameter(name=damage_name_event_code, &
         data=tmpreal)
    damage_event_code = nint(tmpreal)
    
    call fates_params%RetrieveParameter(name=damage_name_canopy_layer_code, &
         data=tmpreal)
    damage_canopy_layer_code = nint(tmpreal)
    
    ! parameters that are arrays of size defined within the params file and thus need allocating as well
    call fates_params%RetrieveParameterAllocate(name=ED_name_history_sizeclass_bin_edges, &
          data=ED_val_history_sizeclass_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_history_ageclass_bin_edges, &
          data=ED_val_history_ageclass_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_history_height_bin_edges, &
          data=ED_val_history_height_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_history_coageclass_bin_edges, &
         data=ED_val_history_coageclass_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_history_damage_bin_edges, &
         data=ED_val_history_damage_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_crop_lu_pft_vector, &
         data=tmp_vector_by_landuse1)

    crop_lu_pft_vector(:) = nint(tmp_vector_by_landuse1(:))
    deallocate(tmp_vector_by_landuse1)

    call fates_params%RetrieveParameterAllocate(name=ED_name_maxpatches_by_landuse, &
         data=tmp_vector_by_landuse2)

    maxpatches_by_landuse(:) = nint(tmp_vector_by_landuse2(:))
    maxpatch_total = sum(maxpatches_by_landuse(:))
    deallocate(tmp_vector_by_landuse2)

    call fates_params%RetrieveParameterAllocate(name=ED_name_max_nocomp_pfts_by_landuse, &
         data=tmp_vector_by_landuse3)

    max_nocomp_pfts_by_landuse(:) = nint(tmp_vector_by_landuse3(:))
    deallocate(tmp_vector_by_landuse3)

    call fates_params%RetrieveParameterAllocate(name=ED_name_hydr_htftype_node, &
         data=hydr_htftype_real)
    allocate(hydr_htftype_node(size(hydr_htftype_real)))
    hydr_htftype_node(:) = nint(hydr_htftype_real(:))
    deallocate(hydr_htftype_real)

    call fates_params%RetrieveParameter(name=name_landuse_grazing_carbon_use_eff, &
         data=landuse_grazing_carbon_use_eff)

    call fates_params%RetrieveParameter(name=name_landuse_grazing_nitrogen_use_eff, &
         data=landuse_grazing_nitrogen_use_eff)

    call fates_params%RetrieveParameter(name=name_landuse_grazing_phosphorus_use_eff, &
         data=landuse_grazing_phosphorus_use_eff)

    call fates_params%RetrieveParameter(name=name_landuse_grazing_maxheight, &
         data=landuse_grazing_maxheight)

    call fates_params%RetrieveParameterAllocate(name=name_landuse_grazing_rate, &
         data=tmp_vector_by_landuse4)

    landuse_grazing_rate(:) = tmp_vector_by_landuse4(:)

    deallocate(tmp_vector_by_landuse4)

  end subroutine FatesReceiveParams
  
  ! =====================================================================================

  subroutine FatesReportParams(is_master)

     logical,intent(in) :: is_master

     character(len=32),parameter :: fmt0 = '(a,(F12.4))'
     character(len=32),parameter :: fmti = '(a,(I4))'
     logical, parameter :: debug_report = .false.
     
     if(debug_report .and. is_master) then
        
        write(fates_log(),*) '-----------  FATES Scalar Parameters -----------------'
        write(fates_log(),fmt0) 'vai_top_bin_width = ',vai_top_bin_width
        write(fates_log(),fmt0) 'vai_width_increase_factor = ',vai_width_increase_factor
        write(fates_log(),fmt0) 'photo_temp_acclim_timescale = ',photo_temp_acclim_timescale
        write(fates_log(),fmt0) 'sdlng_emerg_h2o_timescale = ', sdlng_emerg_h2o_timescale
        write(fates_log(),fmt0) 'sdlng_mort_par_timescale = ', sdlng_mort_par_timescale
        write(fates_log(),fmt0) 'sdlng_mdd_timescale = ', sdlng_mdd_timescale
        write(fates_log(),fmt0) 'sdlng2sap_par_timescale = ', sdlng2sap_par_timescale
        write(fates_log(),fmt0) 'photo_temp_acclim_timescale (days) = ',photo_temp_acclim_timescale
        write(fates_log(),fmt0) 'photo_temp_acclim_thome_time (years) = ',photo_temp_acclim_thome_time
        write(fates_log(),fmti) 'hydr_htftype_node = ',hydr_htftype_node
        write(fates_log(),fmt0) 'fates_mortality_disturbance_fraction = ',fates_mortality_disturbance_fraction
        write(fates_log(),fmt0) 'ED_val_comp_excln = ',ED_val_comp_excln
        write(fates_log(),fmt0) 'ED_val_vai_top_bin_width = ',ED_val_vai_top_bin_width
        write(fates_log(),fmt0) 'ED_val_vai_width_increase_factor = ',ED_val_vai_width_increase_factor
        write(fates_log(),fmt0) 'ED_val_nignitions = ',ED_val_nignitions
        write(fates_log(),fmt0) 'ED_val_understorey_death = ',ED_val_understorey_death
        write(fates_log(),fmt0) 'ED_val_cwd_fcel = ',ED_val_cwd_fcel
        write(fates_log(),fmt0) 'ED_val_cwd_flig = ',ED_val_cwd_flig
        write(fates_log(),fmt0) 'fates_maintresp_nonleaf_baserate = ', maintresp_nonleaf_baserate
        write(fates_log(),fmt0) 'ED_val_phen_a = ',ED_val_phen_a
        write(fates_log(),fmt0) 'ED_val_phen_b = ',ED_val_phen_b
        write(fates_log(),fmt0) 'ED_val_phen_c = ',ED_val_phen_c
        write(fates_log(),fmt0) 'ED_val_phen_chiltemp = ',ED_val_phen_chiltemp
        write(fates_log(),fmt0) 'ED_val_phen_mindayson = ',ED_val_phen_mindayson
        write(fates_log(),fmt0) 'ED_val_phen_ncolddayslim = ',ED_val_phen_ncolddayslim
        write(fates_log(),fmt0) 'ED_val_phen_coldtemp = ',ED_val_phen_coldtemp
        write(fates_log(),fmt0) 'ED_val_cohort_size_fusion_tol = ',ED_val_cohort_size_fusion_tol
        write(fates_log(),fmt0) 'ED_val_cohort_age_fusion_tol = ',ED_val_cohort_age_fusion_tol
        write(fates_log(),fmt0) 'ED_val_patch_fusion_tol = ',ED_val_patch_fusion_tol
        write(fates_log(),fmt0) 'ED_val_canopy_closure_thresh = ',ED_val_canopy_closure_thresh
        write(fates_log(),fmt0) 'regeneration_model = ',regeneration_model
        write(fates_log(),fmt0) 'dayl_switch = ',dayl_switch      
        write(fates_log(),fmt0) 'stomatal_model = ',stomatal_model
        write(fates_log(),fmt0) 'stomatal_assim_model = ',stomatal_assim_model            
        write(fates_log(),fmt0) 'hydro_kmax_rsurf1 = ',hydr_kmax_rsurf1
        write(fates_log(),fmt0) 'hydro_kmax_rsurf2 = ',hydr_kmax_rsurf2  
        write(fates_log(),fmt0) 'hydro_psi0 = ',hydr_psi0
        write(fates_log(),fmt0) 'hydro_psicap = ',hydr_psicap
        write(fates_log(),fmt0) 'hydro_solver = ',hydr_solver
        write(fates_log(),fmt0) 'bgc_soil_salinity = ', bgc_soil_salinity
        write(fates_log(),fmt0) 'logging_dbhmin = ',logging_dbhmin
        write(fates_log(),fmt0) 'logging_dbhmax = ',logging_dbhmax
        write(fates_log(),fmt0) 'logging_collateral_frac = ',logging_collateral_frac
        write(fates_log(),fmt0) 'logging_coll_under_frac = ',logging_coll_under_frac
        write(fates_log(),fmt0) 'logging_direct_frac = ',logging_direct_frac
        write(fates_log(),fmt0) 'logging_mechanical_frac = ',logging_mechanical_frac
        write(fates_log(),fmt0) 'logging_event_code = ',logging_event_code
        write(fates_log(),fmt0) 'logging_dbhmax_infra = ',logging_dbhmax_infra
        write(fates_log(),fmt0) 'eca_plant_escalar = ',eca_plant_escalar
        write(fates_log(),fmt0) 'q10_mr = ',q10_mr
        write(fates_log(),fmt0) 'q10_froz = ',q10_froz
        write(fates_log(),fmt0) 'cg_strikes = ',cg_strikes
        write(fates_log(),'(a,L2)') 'active_crown_fire = ',active_crown_fire
        write(fates_log(),fmt0) 'damage_event_code = ',damage_event_code
        write(fates_log(),fmt0) 'damage_canopy_layer_code = ', damage_canopy_layer_code
	write(fates_log(),fmt0) 'landuse_grazing_carbon_use_eff = ', landuse_grazing_carbon_use_eff
        write(fates_log(),fmt0) 'name_landuse_grazing_nitrogen_use_eff = ', name_landuse_grazing_nitrogen_use_eff
        write(fates_log(),fmt0) 'name_landuse_grazing_phosphorus_use_eff = ', name_landuse_grazing_phosphorus_use_eff
        write(fates_log(),fmt0) 'name_landuse_grazing_maxheight = ', name_landuse_grazing_maxheight
        write(fates_log(),fmt0) 'name_landuse_grazing_rate(:) = ', name_landuse_grazing_rate(:)
        write(fates_log(),*) '------------------------------------------------------'

     end if

  end subroutine FatesReportParams

  
end module EDParamsMod
