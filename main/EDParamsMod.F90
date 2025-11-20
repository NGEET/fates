module EDParamsMod

  !
  ! Things related to FATES parameter Constants
  !
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : itrue
  use FatesGlobals        , only : fates_log
  use FatesGlobals        , only : endrun => fates_endrun
  use FatesConstantsMod,    only : fates_unset_r8
  use FatesConstantsMod,    only : n_landuse_cats
  use JSONParameterUtilsMod,only : params_type,param_type
  
  ! CIME Globals
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  
  implicit none
  private
  save

   real(r8),protected, public :: vai_top_bin_width                    ! width in VAI units of uppermost leaf+stem
                                                                      ! layer scattering element in each canopy layer [m2/m2]
   real(r8),protected, public :: vai_width_increase_factor            ! factor by which each leaf+stem scattering element
                                                                      ! increases in VAI width (1 = uniform spacing)
   real(r8),protected, public :: photo_temp_acclim_timescale          ! Length of the window for the exponential moving average (ema)
                                                                      ! of vegetation temperature used in photosynthesis and respiration
                                                                      ! temperature acclimation [days]
   real(r8),protected, public :: photo_temp_acclim_thome_time         ! Length of the window for the long-term exponential moving average (ema)
                                                                      ! of vegetation temperature used in photosynthesis 
                                                                      ! T_home term in Kumarathunge parameterization [years]
   real(r8),protected, public :: sdlng_emerg_h2o_timescale            ! Length of the window for the exponential moving
                                                                      ! average of smp used to calculate seedling emergence
   real(r8),protected, public :: sdlng_mort_par_timescale             ! Length of the window for the exponential moving average 
                                                                      ! of par at the seedling layer used to calculate 
                                                                      ! seedling mortality
   real(r8),protected, public :: sdlng_mdd_timescale                  ! Length of the window for the exponential moving average
                                                                      ! of moisture deficit days used to calculate seedling mortality
   real(r8),protected, public :: sdlng2sap_par_timescale              ! Length of the window for the exponential 
                                                                      ! moving average of par at the seedling layer used to 
                                                                      ! calculate seedling to sapling transition rates
   real(r8),protected, public :: mortality_disturbance_fraction       ! the fraction of canopy mortality that results in disturbance
   real(r8),protected, public :: comp_excln_exp                       ! weighting factor (exponent) for canopy layer exclusion and promotion
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

   logical,protected, public :: active_crown_fire                     ! flag, 1=active crown fire 0=no active crown fire

   real(r8), protected, public :: cg_strikes                          ! fraction of cloud to ground lightning strikes (0-1)

                                                                      ! Global identifier of how nutrients interact with the host land model
                                                                      ! either they are fully coupled, or they generate uptake rates synthetically
                                                                      ! in prescribed mode. In the latter, there is both NO mass removed from the HLM's soil
                                                                      ! BGC N and P pools, and there is also none removed.

   integer, public :: n_uptake_mode
   integer, public :: p_uptake_mode

   real(r8), parameter, public :: soil_tfrz_thresh = -2.0_r8          ! Soil temperature threshold below which hydraulic failure mortality is off (non-hydro only) in degrees C
   
   integer, parameter, public :: nclmax = 3   ! Maximum number of canopy layers allowed
                                              ! We would make this even higher, but making this
                                              ! a little lower keeps the size down on some output arrays
                                              ! For large arrays at patch level we use dynamic allocation

                                                                      ! parameters that govern the VAI (LAI+SAI) bins used in radiative transfer code
   integer, parameter, public :: nlevleaf = 30                        ! number of leaf+stem layers in each canopy layer

   real(r8), public :: dinc_vai(nlevleaf)   = fates_unset_r8          ! VAI bin widths array
   real(r8), public :: dlower_vai(nlevleaf) = fates_unset_r8          ! numericaly (not vertically) lower edges of VAI bins
                                                                      ! starting with zero in the first index, the last bin
                                                                      ! is assumed to be bounded, but a user can override this
                                                                      ! if change a local parameter vai_capping in tree_lai()
                                                                      ! in the allometry module
   
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
   
   ! parameters whose size is defined in the parameter file
   real(r8),protected,allocatable,public :: ED_val_history_sizeclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_ageclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_height_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_coageclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_damage_bin_edges(:)
   
   ! Hydraulics Control Parameters (ONLY RELEVANT WHEN USE_FATES_HYDR = TRUE)
   ! ----------------------------------------------------------------------------------------------

   ! Switch that defines the hydraulic transfer functions for each organ.
   ! campbell_type           = 3
   ! smooth1_campbell_type   = 31
   ! smooth2_campbell_type   = 32
   ! tfs_type                = 1
   ! van Genuchten 1980 model = 2 
   integer, protected,allocatable,public :: hydr_htftype_node(:) 
   
   real(r8),protected,public :: hydr_kmax_rsurf1         !  maximum conducitivity for unit root surface 
                                                         !  soil to root direction (kg water/m2 root area/Mpa/s)
   
   real(r8),protected,public :: hydr_kmax_rsurf2         !  maximum conducitivity for unit root surface 
                                                         !  root to soil direciton (kg water/m2 root area/Mpa/s)

   real(r8),protected,public :: hydr_psi0          !  sapwood water potential at saturation (MPa)

   real(r8),protected,public :: hydr_psicap        !  sapwood water potential at which capillary reserves exhausted (MPa)

   !Soil BGC parameters, mostly used for testing FATES when not coupled to the dynamics bgc hlm
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: bgc_soil_salinity ! site-level soil salinity for FATES when not coupled to dynamic soil BGC of salinity

   ! Integer code that options how damage events are structured
   integer, protected, public :: damage_event_code

   integer,protected,public :: damage_canopy_layer_code  ! Code that changes whether damage affects canopy trees (1), understory trees (2)
   
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
   
   
   ! Logging Control Parameters (ONLY RELEVANT WHEN USE_FATES_LOGGING = TRUE)
   ! ----------------------------------------------------------------------------------------------

   real(r8),protected,public :: logging_dbhmin              ! Minimum dbh at which logging is applied (cm)
                                                            ! Typically associated with harvesting

   real(r8),protected,public :: logging_dbhmax              ! Maximum dbh at which logging is applied (cm)
                                                            ! Typically associated with fire suppression

   real(r8),protected,public :: logging_collateral_frac     ! Ratio of collateral mortality to direct logging mortality

   real(r8),protected,public :: logging_coll_under_frac ! Fraction of understory plants that die when logging disturbance
                                                 ! is generated
   
   real(r8),protected,public :: logging_direct_frac         ! Fraction of stems logged per event

   real(r8),protected,public :: logging_mechanical_frac         ! Fraction of stems logged per event

   real(r8),protected,public :: logging_event_code          ! Code that options how logging events are structured 
   
   real(r8),protected,public :: logging_dbhmax_infra        ! "Tree diameter, above which infrastructure from logging does not impact damage or mortality.
   
   real(r8),protected,public :: logging_export_frac        ! "fraction of trunk product being shipped offsite, the 
                                                    ! leftovers will be left onsite as large CWD

   real(r8),protected,public :: eca_plant_escalar  ! scaling factor for plant fine root biomass to 
                                               ! calculate nutrient carrier enzyme abundance (ECA)

   public :: TransferParamsGeneric
   public :: FatesReportParams
   public :: GetNVegLayers

   
 contains


   function GetNVegLayers(treevai) result(nv)

     real(r8) :: treevai  ! The LAI+SAI of the cohort (m2/m2)
     integer  :: nv

     nv = count(treevai .gt. dlower_vai(:))

   end function GetNVegLayers
     
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
    mortality_disturbance_fraction        = nan
    comp_excln_exp                        = nan
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
    max_cohort_per_patch                  = -9
    hydr_kmax_rsurf1                      = nan
    hydr_kmax_rsurf2                      = nan
    hydr_psi0                             = nan
    hydr_psicap                           = nan
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
    dev_arbitrary                         = nan
    damage_event_code                     = -9
    damage_canopy_layer_code              = -9
    landuse_grazing_carbon_use_eff        = nan
    landuse_grazing_nitrogen_use_eff      = nan
    landuse_grazing_phosphorus_use_eff    = nan
    landuse_grazing_maxheight             = nan
    landuse_grazing_rate(:)               = nan

  end subroutine FatesParamsInit

  ! =====================================================================================

  subroutine TransferParamsGeneric(pstruct)

    ! -----------------------------------------------------------------------------------
    ! Transfer non-specific/generic parameter values from the data structure "pstruct"
    ! to named primitive data structures
    ! -----------------------------------------------------------------------------------

    implicit none

    type(params_type) :: pstruct         ! Data structure containing all parameters and dimensions
    type(param_type),pointer :: param_p  ! Pointer to one specific parameter

    integer :: num_hydr_organ
    
    call FatesParamsInit()

    num_hydr_organ = pstruct%GetDimSizeFromName('fates_hydr_organs')
    
    param_p => pstruct%GetParamFromName("fates_fire_active_crown_fire")
    active_crown_fire = (param_p%i_data_scalar == itrue)
    
    param_p => pstruct%GetParamFromName("fates_fire_cg_strikes")
    cg_strikes = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_dev_arbitrary")
    dev_arbitrary = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_leaf_photo_temp_acclim_timescale")
    photo_temp_acclim_timescale = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_leaf_photo_temp_acclim_thome_time")
    photo_temp_acclim_thome_time = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_trs_seedling_emerg_h2o_timescale")
    sdlng_emerg_h2o_timescale = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_trs_seedling_mort_par_timescale")
    sdlng_mort_par_timescale = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_trs_seedling_mdd_timescale")
    sdlng_mdd_timescale = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_trs_seedling2sap_par_timescale")
    sdlng2sap_par_timescale = param_p%r_data_scalar
   
    param_p => pstruct%GetParamFromName("fates_mort_disturb_frac")
    mortality_disturbance_fraction = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_comp_excln")
    comp_excln_exp = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_vai_top_bin_width")
    vai_top_bin_width = param_p%r_data_scalar
   
    param_p => pstruct%GetParamFromName("fates_vai_width_increase_factor")
    vai_width_increase_factor = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_fire_nignitions")
    ED_val_nignitions = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_mort_understorey_death")
    ED_val_understorey_death = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_frag_cwd_fcel")
    ED_val_cwd_fcel = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_frag_cwd_flig")
    ED_val_cwd_flig = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_maintresp_nonleaf_baserate")
    maintresp_nonleaf_baserate = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_gddthresh_a")
    ED_val_phen_a = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_gddthresh_b")
    ED_val_phen_a = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_gddthresh_c")
    ED_val_phen_a = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_chilltemp")
    ED_val_phen_chiltemp = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_mindayson")
    ED_val_phen_mindayson = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_ncolddayslim")
    ED_val_phen_ncolddayslim = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_phen_coldtemp")
    ED_val_phen_coldtemp = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_cohort_size_fusion_tol")
    ED_val_cohort_size_fusion_tol = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_cohort_age_fusion_tol")
    ED_val_cohort_age_fusion_tol = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_patch_fusion_tol")
    ED_val_patch_fusion_tol = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_canopy_closure_thresh")
    ED_val_canopy_closure_thresh = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_q10_mr")
    q10_mr = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_q10_froz")
    q10_froz = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_history_sizeclass_bin_edges")
    allocate(ED_val_history_sizeclass_bin_edges(size(param_p%r_data_1d,dim=1)))
    ED_val_history_sizeclass_bin_edges(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_history_ageclass_bin_edges")
    allocate(ED_val_history_ageclass_bin_edges(size(param_p%r_data_1d,dim=1)))
    ED_val_history_ageclass_bin_edges(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_history_height_bin_edges")
    allocate(ED_val_history_height_bin_edges(size(param_p%r_data_1d,dim=1)))
    ED_val_history_height_bin_edges(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_history_coageclass_bin_edges")
    allocate(ED_val_history_coageclass_bin_edges(size(param_p%r_data_1d,dim=1)))
    ED_val_history_coageclass_bin_edges(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_history_damage_bin_edges")
    allocate(ED_val_history_damage_bin_edges(size(param_p%r_data_1d,dim=1)))
    ED_val_history_damage_bin_edges(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_landuse_crop_lu_pft_vector")
    crop_lu_pft_vector(:) =  param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName("fates_maxpatches_by_landuse")
    maxpatches_by_landuse(:) = param_p%i_data_1d(:)
    maxpatch_total = sum(maxpatches_by_landuse(:))
    
    param_p => pstruct%GetParamFromName("fates_max_nocomp_pfts_by_landuse")
    max_nocomp_pfts_by_landuse(:) = param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName("fates_hydro_htftype_node")
    allocate(hydr_htftype_node(num_hydr_organ))
    hydr_htftype_node(:) = param_p%i_data_1d(:)
    
    ! if use_hydro
    param_p => pstruct%GetParamFromName("fates_hydro_kmax_rsurf1")
    hydr_kmax_rsurf1 = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_hydro_kmax_rsurf2")
    hydr_kmax_rsurf2 = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_hydro_psi0")
    hydr_psi0 = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_hydro_psicap")
    hydr_psicap = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_soil_salinity")
    bgc_soil_salinity = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_damage_event_code")
    damage_event_code = param_p%i_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_damage_canopy_layer_code")
    damage_canopy_layer_code = param_p%i_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_maxcohort")
    max_cohort_per_patch = param_p%i_data_scalar

    param_p => pstruct%GetParamFromName("fates_landuse_logging_dbhmin")
    logging_dbhmin = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_dbhmax")
    logging_dbhmax = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_collateral_frac")
    logging_collateral_frac = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_coll_under_frac")
    logging_coll_under_frac = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_direct_frac")
    logging_direct_frac = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_mechanical_frac")
    logging_mechanical_frac = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_event_code")
    logging_event_code = param_p%i_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_dbhmax_infra")
    logging_dbhmax_infra = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_logging_export_frac")
    logging_export_frac  = param_p%r_data_scalar

    param_p => pstruct%GetParamFromName("fates_landuse_grazing_rate")
    landuse_grazing_rate(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName("fates_landuse_grazing_carbon_use_eff")
    landuse_grazing_carbon_use_eff = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_grazing_maxheight")
    landuse_grazing_maxheight = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_grazing_nitrogen_use_eff")
    landuse_grazing_nitrogen_use_eff = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_landuse_grazing_phosphorus_use_eff")
    landuse_grazing_phosphorus_use_eff = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_cnp_eca_plant_escalar")
    eca_plant_escalar = param_p%r_data_scalar

    return
  end subroutine TransferParamsGeneric


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
        write(fates_log(),fmt0) 'mortality_disturbance_fraction = ',mortality_disturbance_fraction
        write(fates_log(),fmt0) 'comp_excln_exp = ',comp_excln_exp
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
        write(fates_log(),fmt0) 'hydro_kmax_rsurf1 = ',hydr_kmax_rsurf1
        write(fates_log(),fmt0) 'hydro_kmax_rsurf2 = ',hydr_kmax_rsurf2  
        write(fates_log(),fmt0) 'hydro_psi0 = ',hydr_psi0
        write(fates_log(),fmt0) 'hydro_psicap = ',hydr_psicap
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
        write(fates_log(),fmt0) 'landuse_grazing_carbon_use_eff = ',landuse_grazing_carbon_use_eff
        write(fates_log(),fmt0) 'landuse_grazing_nitrogen_use_eff = ', landuse_grazing_nitrogen_use_eff
        write(fates_log(),fmt0) 'landuse_grazing_phosphorus_use_eff = ', landuse_grazing_phosphorus_use_eff
        write(fates_log(),fmt0) 'landuse_grazing_maxheight = ', landuse_grazing_maxheight
        write(fates_log(),fmt0) 'landuse_grazing_rate(:) = ', landuse_grazing_rate(:)
        write(fates_log(),*) '------------------------------------------------------'

     end if

   end subroutine FatesReportParams

  
end module EDParamsMod
