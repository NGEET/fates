module FatesConstantsMod

  ! This module is used to define global _immutable_ data. Everything in this module 
  ! must have the parameter attribute.

  implicit none
  private 

  ! KINDS
  integer, parameter, public :: fates_r8  = selected_real_kind(12) ! 8 byte real
  integer, parameter, public :: fates_int = selected_int_kind(8)   ! 4 byte int

  ! STRING LENGTH
  integer, parameter, public :: fates_avg_flag_length     = 3
  integer, parameter, public :: fates_short_string_length = 32
  integer, parameter, public :: fates_long_string_length  = 199

  ! INITIALIZING AND TESTING UNSET VALUES
  integer,        parameter, public :: fates_unset_int       = -9999           ! unset integers
  real(fates_r8), parameter, public :: fates_unset_r8        = -1.e36_fates_r8 ! unset r8s
  real(fates_r8), parameter, public :: fates_check_param_set = 9.9e32_fates_r8 ! unset in parameter file 

  ! BOOLS
  integer, parameter, public :: itrue  = 1 ! integer equivalent of true  (in case some compilers dont auto convert)
  integer, parameter, public :: ifalse = 0 ! integer equivalent of false (in case come compilers dont auto convert)

  ! INDEXING AND LABELS
  integer, parameter, public :: n_dbh_bins              = 6  ! no. of dbh bins used when comparing patches
  integer, parameter, public :: maxpft                  = 16 ! maximum number of PFTs allowed
                                                             !  the parameter file may determine that fewer
                                                             !  are used, but this helps allocate scratch
                                                             !  space and output arrays.
  integer, parameter, public  :: numWaterMem            = 10 ! watermemory saved as site level var
  integer, parameter, public  :: numlevsoil_max         = 30 ! scratch space used for static arrays
                                                             !  The actual number of soil layers should not exceed this
  integer, parameter, public :: num_vegtemp_mem        = 10 ! window of time over which we track temp for cold sensecence (days)
  integer, parameter, public :: n_rad_stream_types      = 2  ! number of radiation streams used (direct/diffuse)
  integer, parameter, public :: idirect                 = 1  ! array index for direct radiation
  integer, parameter, public :: idiffuse                = 2  ! array index for diffuse radiation
  integer, parameter, public :: n_dist_types            = 3  ! disturbance modes 1) tree-fall, 2) fire, 3) logging
  integer, parameter, public :: leaves_on               = 2  ! flag specifying that a deciduous plant has leaves
                                                             !   and should be allocating to them as well
  integer, parameter, public :: leaves_off              = 1  ! flag specifying that a deciduous plant has dropped
                                                             !   its leaves and should not be trying to allocate
                                                             !   towards any growth
  integer, parameter, public :: phen_cstat_nevercold    = 0  ! site has not experienced a cold period over a large number
                                                             !  of days, leaves are dropped and flagged as non-cold region
  integer, parameter, public :: phen_cstat_iscold       = 1  ! site is in a cold-state where leaves should have fallen
  integer, parameter, public :: phen_cstat_notcold      = 2  ! site is in a warm-state where leaves are allowed to flush
  integer, parameter, public :: phen_dstat_timeoff      = 0  ! site has off due to time exceedance (drought phenology)
  integer, parameter, public :: phen_dstat_moistoff     = 1  ! site has off due to moisture avail (drought phenology)
  integer, parameter, public :: phen_dstat_moiston      = 2  ! site has on due to moisture avail (drought phenology)
  integer, parameter, public :: phen_dstat_timeon       = 3  ! site has on due to time exceedance (drought phenology)
  integer, parameter, public :: ican_upper              = 1  ! nominal index for the upper canopy
  integer, parameter, public :: ican_ustory             = 2  ! nominal index for diagnostics that refer to understory layers 
                                                             !  (all layers that are not the top canopy layer)
  integer, parameter, public :: nocomp_bareground       = 0  ! bareground label for no competition mode
  integer, parameter, public :: dtype_ifall             = 1  ! index for naturally occuring tree-fall generated event
  integer, parameter, public :: dtype_ifire             = 2  ! index for fire generated disturbance event
  integer, parameter, public :: dtype_ilog              = 3  ! index for logging generated disturbance event
  integer, parameter, public :: prescribed_p_uptake     = 1  ! flag specifying prescribed phosophorous uptake/turnover 
  integer, parameter, public :: coupled_p_uptake        = 2  ! flag specifying coupled phosophorous uptake/turnover from host model
  integer, parameter, public :: prescribed_n_uptake     = 1  ! flag specifying prescribed nitrogen uptake/turnover 
  integer, parameter, public :: coupled_n_uptake        = 2  ! flag specifying coupled nitrogen uptake/turnover from host model
  integer, parameter, public :: coupled_np_comp_scaling = 1  ! flag signaling at least 1 chemical element (ie N or P)
  integer, parameter, public :: trivial_np_comp_scaling = 2  ! flag definition indicates that either
                                                             !  nutrients are turned off in FATES, or, that the
                                                             !  plants are not coupled with below ground chemistry. In
                                                             !  this situation, we send token boundary condition information
  
  ! integer labels for specifying photosynthetic acclimation model
  integer, parameter, public :: photosynth_acclim_model_none                   = 1
  integer, parameter, public :: photosynth_acclim_model_kumarathunge_etal_2019 = 2

  ! integer labels for specifying leaf maintenance respiration models
  integer, parameter, public :: lmrmodel_ryan_1991         = 1
  integer, parameter, public :: lmrmodel_atkin_etal_2017   = 2

  ! labels for patch disturbance history
  integer, parameter, public :: n_anthro_disturbance_categories = 2
  integer, parameter, public :: primaryforest                   = 1
  integer, parameter, public :: secondaryforest                 = 2

  ! integer labels for specifying harvest units
  integer, parameter, public :: hlm_harvest_area_fraction = 1 ! harvest by area
  integer, parameter, public :: hlm_harvest_carbon        = 2 ! harvest based on carbon extracted

  ! ERROR TOLERANCES

  ! Allowable error in carbon allocations, should be applied to estimates
  ! of carbon conservation in units of kgC/plant.  This gives an effective
  ! error tolerance of 1 microgram.
  real(fates_r8), parameter, public :: calloc_abs_error = 1.0e-9_fates_r8

  ! Rounding errors seem to hover around 1e-15 for the gnu compiler
  ! when not applying compiler directives for safe math.  An example
  ! of this is taking a vector of numbers, dividing through by their sum,
  ! multiplying each by their original sum, and then seeing if their addition
  ! matches the original sum.  Other simple examples of rounding errors
  ! are simply changing the orders:  a*b*c .ne. a*c*b
  ! This value here is used as an error expectation comparison
  ! for multiplication/division procedures, also allowing for 3 orders
  ! of magnitude of buffer error (ie instead of 1e-15)
  real(fates_r8), parameter, public :: rsnbl_math_prec = 1.0e-12_fates_r8

  ! This is the precision of 8byte reals (~1e-308)
  real(fates_r8), parameter, public :: tinyr8 = tiny(1.0_fates_r8)

  ! We mostly use this in place of logical comparisons
  ! between reals with zero, as the chances are their
  ! precisions are preventing perfect zero in comparison
  real(fates_r8), parameter, public :: nearzero = 1.0e-30_fates_r8

  ! UNIT CONVERSIONS
  real(fates_r8), parameter, public :: umolC_to_kgC    = 12.0E-9_fates_r8 ! umols of Carbon to kg of Carbon
                                                                          !  (1 mol = 12 g)
                                                                          !  we do not use umolC_per_kg because it is
                                                                          !  a non-terminating decimal
  real(fates_r8), parameter, public :: mg_per_kg       = 1.0e6_fates_r8   ! milligrams per kilgrams
  real(fates_r8), parameter, public :: g_per_kg        = 1000.0_fates_r8  ! grams per kilogram
  real(fates_r8), parameter, public :: kg_per_g        = 0.001_fates_r8   ! kilograms per gram
  real(fates_r8), parameter, public :: mg_per_g        = 1000.0_fates_r8  ! milligrams per gram
  real(fates_r8), parameter, public :: kg_per_Megag    = 1000.0_fates_r8  ! kilgrams per Megagram
  real(fates_r8), parameter, public :: umol_per_mmol   = 1000.0_fates_r8  ! micromoles per millimole
  real(fates_r8), parameter, public :: mmol_per_mol    = 1000.0_fates_r8  ! millimoles per mole
  real(fates_r8), parameter, public :: umol_per_mol    = 1.0E6_fates_r8   ! micromoles per mole
  real(fates_r8), parameter, public :: mol_per_umol    = 1.0E-6_fates_r8  ! moles per micromole
  real(fates_r8), parameter, public :: umol_per_kmol   = 1.0E9_fates_r8   ! micromoles per kilomole
  real(fates_r8), parameter, public :: m_per_mm        = 1.0E-3_fates_r8  ! meters per millimeter
  real(fates_r8), parameter, public :: mm_per_m        = 1.0E3_fates_r8   ! millimeters per meter
  real(fates_r8), parameter, public :: m_per_cm        = 1.0E-2_fates_r8  ! meters per centimeter
  real(fates_r8), parameter, public :: m2_per_ha       = 1.0e4_fates_r8   ! meters squared per hectare
  real(fates_r8), parameter, public :: ha_per_m2       = 1.0e-4_fates_r8  ! hectares per meters squared
  real(fates_r8), parameter, public :: m2_per_km2      = 1.0e6_fates_r8   ! meters squared per kilometers squared
  real(fates_r8), parameter, public :: cm2_per_m2      = 10000.0_fates_r8 ! centimeters squared per meters squared
  real(fates_r8), parameter, public :: m3_per_mm3      = 1.0E-9_fates_r8  ! meters cubed per millimeters cubed
  real(fates_r8), parameter, public :: m3_per_cm3      = 1.0E-6_fates_r8  ! meters cubed per centimeters cubed
  real(fates_r8), parameter, public :: cm3_per_m3      = 1.0E6_fates_r8   ! centimeters cubed per meters cubed
  real(fates_r8), parameter, public :: sec_per_min     = 60.0_fates_r8    ! seconds per minute
  real(fates_r8), parameter, public :: sec_per_day     = 86400.0_fates_r8 ! seconds per day
  real(fates_r8), parameter, public :: days_per_year   = 365.00_fates_r8  ! days per year (assume HLM uses 365 day calendar)
                                                                          !   if we need to link to 365.25 calendared HLM,
                                                                          !   rewire to pass through interface
  real(fates_r8), parameter, public :: months_per_year = 12.0_fates_r8    ! months per year
  real(fates_r8), parameter, public :: J_per_kJ        = 1000.0_fates_r8  ! Jules per kiloJoules
  real(fates_r8), parameter, public :: pa_per_mpa      = 1.e6_fates_r8    ! megapascals to pascals
  real(fates_r8), parameter, public :: mpa_per_pa      = 1.e-6_fates_r8   ! pascals to megapascals
  real(fates_r8), parameter, public :: days_per_sec    = 1.0_fates_r8/sec_per_day   ! days per second
  real(fates_r8), parameter, public :: years_per_day   = 1.0_fates_r8/days_per_year ! years per day

  ! PHYSICAL CONSTANTS
  real(fates_r8), parameter, public :: rgas_J_K_kmol           = 8314.4598_fates_r8 ! universal gas constant [J/K/kmol]
  real(fates_r8), parameter, public :: rgas_J_K_mol            = 8.3144598_fates_r8 ! universal gas constant [J/k/mol]
  real(fates_r8), parameter, public :: t_water_freeze_k_1atm   = 273.15_fates_r8    ! freezing point of water at 1 atm [K]
  real(fates_r8), parameter, public :: t_water_freeze_k_triple = 273.16_fates_r8    ! freezing point of water at triple point [K]
  real(fates_r8), parameter, public :: dens_fresh_liquid_water = 1.0E3_fates_r8     ! density of fresh liquid water [kg/m3]
  real(fates_r8), parameter, public :: molar_mass_water        = 18.0_fates_r8      ! molar mass of water [g/mol]
  real(fates_r8), parameter, public :: molar_mass_ratio_vapdry = 0.622_fates_r8     ! approximate molar mass of water vapor to dry air (-)
  real(fates_r8), parameter, public :: grav_earth              = 9.8_fates_r8       ! gravity constant on earth [m/s]

  ! NUMERICAL INQUIRY
  real(fates_r8), parameter, public :: fates_huge = huge(g_per_kg)
  real(fates_r8), parameter, public :: fates_tiny = tiny(g_per_kg)

  ! GEOMETRIC CONSTANTS
  real(fates_r8), parameter, public :: pi_const = 3.14159265359_fates_r8 ! pi

  ! MODES
  logical, parameter, public :: do_fates_salinity     = .false. ! Flag to turn on/off salinity effects on the effective "btran"
                                                                !   btran stress function
  logical, parameter, public :: homogenize_seed_pfts  = .false. ! special mode to cause PFTs to create seed mass of all currently-existing PFTs

  ! MODEL PARAMETERS
  real(fates_r8), parameter, public :: secondary_age_threshold = 94._fates_r8     ! less than this value is young secondary land
                                                                                  !   based on average age of global
                                                                                  !   secondary 1900s land in hurtt-2011 [years]
  real(fates_r8), parameter, public :: init_recruit_trim       = 0.8_fates_r8     ! the initial trimming value that
                                                                                  !   new recruits start with
  real(fates_r8), parameter, public :: area                    = 10000.0_fates_r8 ! notional area of simulated forest [m2]
  real(fates_r8), parameter, public :: area_inv                = 1/area           ! inverse of the notion area (faster math) [/m2]
  real(fates_r8), parameter, public :: init_site_GDD           = 30.0_fates_r8    ! initial growing degree-days for sites (non-restart)
  integer,        parameter, public :: init_cleafon            = 100
  integer,        parameter, public :: init_cleafoff           = 300
  integer,        parameter, public :: init_dleafon            = 100
  integer,        parameter, public :: init_dleafoff           = 300
  real(fates_r8), parameter, public :: init_watermem           = 0.5_fates_r8     ! initial moisture in the site water_memory array

  ! PATCH FUSION 
  real(fates_r8), parameter, public :: force_patchfuse_min_biomass                 = 0.005_fates_r8 ! min biomass below which to force-fuse patches [kg/m2 patch area]
  real(fates_r8), parameter, public :: patchfusion_dbhbin_loweredges(n_dbh_bins)   =       &        ! array of bin lower edges for comparing patches [cm]
                                                                    (/0._fates_r8, 5._fates_r8, 20._fates_r8, 50._fates_r8, 100._fates_r8, 150._fates_r8/)                 
  real(fates_r8), parameter, public :: patch_fusion_tolerance_relaxation_increment = 1.1_fates_r8   ! amount by which to increment patch fusion threshold
  real(fates_r8), parameter, public :: max_age_of_second_oldest_patch              = 200._fates_r8  ! age above which to combine all patches [years]

  ! COHORT FUSION AND TERMINATION
  real(fates_r8), parameter, public :: hitemax               = 30.0_fates_r8   ! max dbh value used in hgt profile comparison [cm]
  integer,        parameter, public :: n_hite_bins           = 60              ! no. of hite bins used to distribute LAI
  real(fates_r8), parameter, public :: min_npm2              = 1.0E-7_fates_r8 ! minimum cohort number density before termination [/m2]
  real(fates_r8), parameter, public :: min_patch_area        = 0.01_fates_r8   ! smallest allowable patch area before termination [m2]
  real(fates_r8), parameter, public :: min_patch_area_forced = 0.0001_fates_r8 ! patch termination will not fuse the youngest patch
                                                                               !   if the area is less than min_patch_area.
                                                                               !   however, it is allowed to fuse the youngest patch
                                                                               !   if the fusion area is less than min_patch_area_forced [m2]
  real(fates_r8), parameter, public :: min_nppatch    = min_npm2*min_patch_area ! minimum number of cohorts per patch [/patch]
  real(fates_r8), parameter, public :: min_n_safemath = 1.0E-12_fates_r8       ! in some cases, we want to immediately remove super small
                                                                                ! number densities of cohorts to prevent FPEs

  ! This is the community level amount of spread expected in nearly-bare-ground
  ! and inventory starting modes.
  ! These are used to initialize only. These values will scale between
  ! the PFT defined maximum and minimum crown area scaing parameters.
  !
  ! A value of 1 indicates that
  ! plants should have crown areas at maximum spread for their size and PFT.
  ! A value of 0 means that they have the least amount of spread for their
  ! size and PFT.
  real(fates_r8), parameter, public :: init_spread_near_bare_ground = 1.0_fates_r8
  real(fates_r8), parameter, public :: init_spread_inventory        = 0.0_fates_r8

end module FatesConstantsMod
