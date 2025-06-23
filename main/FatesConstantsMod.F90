module FatesConstantsMod

  ! This module is used to define global _immutable_ data. Everything in
  ! this module must have the parameter attribute.

  implicit none
  private ! Modules are private by default

  ! kinds
  integer, parameter, public :: fates_r8 = selected_real_kind(12) ! 8 byte real
  integer, parameter, public :: fates_int = selected_int_kind(8)  ! 4 byte int

  ! string lengths
  integer, parameter, public :: fates_avg_flag_length = 3
  integer, parameter, public :: fates_short_string_length = 32
  integer, parameter, public :: fates_long_string_length = 199

  ! Used to initialize and test unset integers
  integer, parameter, public :: fates_unset_int = -9999

  ! Used to initialize and test unset r8s
  real(fates_r8), parameter, public :: fates_unset_r8 = -1.e36_fates_r8

  ! Used to check if a parameter was specified in the parameter file (or left as _)
  real(fates_r8), parameter, public :: fates_check_param_set = 9.9e32_fates_r8

  ! Integer equivalent of true  (in case some compilers dont auto convert)
  integer, parameter, public :: itrue = 1

  ! Integer equivalent of false (in case come compilers dont auto convert)
  integer, parameter, public :: ifalse = 0

  ! the parameter file may determine that fewer
  ! are used, but this helps allocate scratch
  ! space and output arrays.

  integer , parameter, public       :: N_DBH_BINS = 6  ! no. of dbh bins used when comparing patches
  real(fates_r8), parameter, public :: patchfusion_dbhbin_loweredges(N_DBH_BINS) = &
  (/0._fates_r8, 5._fates_r8, 20._fates_r8, 50._fates_r8, 100._fates_r8, 150._fates_r8/) ! array of bin lower edges for comparing patches

  
  real(fates_r8), parameter, public :: min_vai_bin_sum = 5.0_fates_r8   ! The sum of vai increments used to discretize the canopy vertically                                                                                                                           ! must be larger than this number. 
  
  integer , parameter, public :: N_DIST_TYPES = 4          ! Disturbance Modes 1) tree-fall, 2) fire, 3) logging, 4) land-use change
  integer , parameter, public :: dtype_ifall  = 1          ! index for naturally occuring tree-fall generated event
  integer , parameter, public :: dtype_ifire  = 2          ! index for fire generated disturbance event
  integer , parameter, public :: dtype_ilog   = 3          ! index for logging generated disturbance event
  integer , parameter, public :: dtype_ilandusechange = 4  ! index for land use change disturbance (not including logging)

  ! Labels for patch land use type information
  integer, parameter, public :: n_landuse_cats = 5
  integer, parameter, public :: primaryland = 1
  integer, parameter, public :: secondaryland = 2
  integer, parameter, public :: rangeland = 3
  integer, parameter, public :: pastureland = 4
  integer, parameter, public :: cropland = 5
  logical, parameter, dimension(n_landuse_cats), public :: is_crop = [.false., .false.,.false.,.false.,.true.]
  integer, parameter, public :: n_crop_lu_types = 1

  ! Bareground nocomp land use label
  integer, parameter, public :: nocomp_bareground_land = 0  ! not a real land use type, only for labeling any bare-ground nocomp patches

  ! Bareground nocomp PFT label for no competition mode
  integer, parameter, public :: nocomp_bareground = 0

  integer, parameter, public :: leaves_on  = 2  ! Flag specifying that a deciduous plant has leaves
                                                ! and should be allocating to them as well
  integer, parameter, public :: leaves_off = 1  ! Flag specifying that a deciduous plant has dropped
                                                ! its leaves and should not be trying to allocate
                                                ! towards any growth.
  integer, parameter, public :: leaves_shedding = 3  ! Flag specifying that a deciduous plant has leaves
                                                     ! but is shedding them (partial shedding). This plant
                                                     ! should not allocate carbon towards growth or 
                                                     ! reproduction.
integer, parameter, public :: ihard_stress_decid = 1 ! If the PFT is stress (drought) deciduous,
                                                     !  this flag is used to tell that the PFT
                                                     !  is a "hard" deciduous (i.e., the plant
                                                     !  has only two statuses, the plant either
                                                     !  sheds all leaves when it's time, or seeks
                                                     !  to flush the leaves back to allometry 
                                                     !  when conditions improve.
integer, parameter, public :: isemi_stress_decid = 2 ! If the PFT is stress (drought) deciduous,
                                                     !  this flag is used to tell that the PFT
                                                     !  is a semi-deciduous (i.e., the plant
                                                     !  can downregulate the amount of leaves
                                                     !  relative to the allometry based on 
                                                     !  soil moisture conditions. It can still
                                                     !  shed all leaves if conditions are very
                                                     !  dry.

  integer, parameter, public :: ican_upper = 1  ! nominal index for the upper canopy
  integer, parameter, public :: ican_ustory = 2 ! nominal index for diagnostics that refer to understory layers 
                                                !  (all layers that are not the top canopy layer)

  ! Flags specifying how phosphorous uptake and turnover interacts
  ! with the host model.
  integer, public, parameter :: prescribed_p_uptake = 1
  integer, public, parameter :: coupled_p_uptake    = 2

  ! Flags specifying how nitrogen uptake and turnover interacts
  ! with the host model.
  integer, public, parameter :: prescribed_n_uptake = 1
  integer, public, parameter :: coupled_n_uptake    = 2
  integer, public, parameter :: coupled_np_comp_scaling = 1 ! This flag signals that at least 1 chemical element (ie N or P)
  
  !Flags specifying how tree regeneration works
  
  integer, public, parameter :: TRS_no_seedling_dyn = 3                          ! Constant defining the Tree Recruitment
                                                                                 ! Scheme switch. This value turns on 
                                                                                 ! size-based reproductive allocation 
                                                                                 ! and allocation to non-seed 
                                                                                 ! reproductive biomass, but does not turn 
                                                                                 ! on seedling dynamics.
  integer, public, parameter :: TRS_regeneration = 2                             ! Constant defining the Tree Recruitment
                                                                                 ! Scheme switch. Turns on full TRS.
  integer, public, parameter :: default_regeneration = 1                         ! Constant defining FATES's default 
                                                                                 ! regeneration scheme switch.
  real(fates_r8), public, parameter :: min_max_dbh_for_trees = 15._fates_r8      ! If pfts have a max dbh less 
                                                                                 ! than this value FATES 
                                                                                 ! will use the default regeneration scheme.
                                                                                 ! Avoids TRS for shrubs / grasses.

  integer, public, parameter :: trivial_np_comp_scaling = 2 ! This flag definition indicates that either
                                                            ! nutrients are turned off in FATES, or, that the
                                                            ! plants are not coupled with below ground chemistry. In
                                                            ! this situation, we send token boundary condition information.


  ! This flag specifies the scaling of how we present
  ! nutrient competitors to the HLM's soil BGC model
  
  integer, public :: fates_np_comp_scaling = fates_unset_int

  real(fates_r8), parameter, public :: secondary_age_threshold = 94._fates_r8 ! less than this value is young secondary land
                                                            ! based on average age of global
                                                            ! secondary 1900s land in hurtt-2011

  ! integer labels for specifying harvest units
  integer, parameter, public :: photosynth_acclim_model_none = 1
  integer, parameter, public :: photosynth_acclim_model_kumarathunge_etal_2019 = 2

  ! integer labels for specifying harvest units
  integer, parameter, public :: hlm_harvest_area_fraction = 1 ! Code for harvesting by area
  integer, parameter, public :: hlm_harvest_carbon = 2 ! Code for harvesting based on carbon extracted.

  ! integer labels for specifying harvest debt status
  integer, parameter, public :: fates_no_harvest_debt = 0
  integer, parameter, public :: fates_with_harvest_debt = 1
  integer, parameter, public :: fates_bypass_harvest_debt = 2  ! Do not calculate harvest debt for area based harvest

  ! integer labels for specifying leaf maintenance respiration models
  integer, parameter, public :: lmrmodel_ryan_1991         = 1
  integer, parameter, public :: lmrmodel_atkin_etal_2017   = 2

  ! integer labels for specifying carbon starvation model
  integer, parameter, public :: cstarvation_model_lin = 1 ! Linear scaling
  integer, parameter, public :: cstarvation_model_exp = 2 ! Exponential scaling

  ! Error Tolerances

  ! Allowable error in carbon allocations, should be applied to estimates
  ! of carbon conservation in units of kgC/plant.  This gives an effective
  ! error tolerance of 1 microgram.
  real(fates_r8), parameter, public :: calloc_abs_error = 1.0e-9_fates_r8

  ! area tolerance checks
  real(fates_r8), parameter, public :: area_error_1     = 1.0e-16_fates_r8 ! error tolerance for area checks (canopy, patch)
  real(fates_r8), parameter, public :: area_error_2     = 1.0e-12_fates_r8 ! error tolerance for tree lai checks
  real(fates_r8), parameter, public :: area_error_3     = 10.e-9_fates_r8  ! error tolerance for area checks (canopy, patch)
  real(fates_r8), parameter, public :: area_error_4     = 1.0e-10_fates_r8 ! error tolerance for area checks

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

  ! in nocomp simulations, what is the minimum PFT fraction for any given land use type?
  real(fates_r8), parameter, public :: min_nocomp_pftfrac_perlanduse = 0.01_fates_r8

  ! This is the precision of 8byte reals (~1e-308)
  real(fates_r8), parameter, public :: tinyr8 = tiny(1.0_fates_r8)

  ! We mostly use this in place of logical comparisons
  ! between reals with zero, as the chances are their
  ! precisions are preventing perfect zero in comparison
  real(fates_r8), parameter, public :: nearzero = 1.0e-30_fates_r8

  ! Unit conversion constants:

  ! Conversion factor umols of Carbon -> kg of Carbon (1 mol = 12g)
  ! We do not use umolC_per_kg because it is a non-terminating decimal

  real(fates_r8), parameter, public :: umolC_to_kgC = 12.0E-9_fates_r8

  ! Conversion factor: miligrams per kilogram
  real(fates_r8), parameter, public :: mg_per_kg    = 1.0e6_fates_r8

  ! Conversion factor: grams per kilograms
  real(fates_r8), parameter, public :: g_per_kg = 1000.0_fates_r8

  ! Conversion factor: kilograms per gram
  real(fates_r8), parameter, public :: kg_per_g = 0.001_fates_r8

  ! Conversion factor: miligrams per grams
  real(fates_r8), parameter, public :: mg_per_g = 1000.0_fates_r8

  ! Conversion factor: kilograms per Megagram
  real(fates_r8), parameter, public :: kg_per_Megag = 1000.0_fates_r8

  ! Conversion factor: micromoles per milimole
  real(fates_r8), parameter, public :: umol_per_mmol = 1000.0_fates_r8

  ! Conversion factor: milimoles per mole
  real(fates_r8), parameter, public :: mmol_per_mol = 1000.0_fates_r8

  ! Conversion factor: micromoles per mole
  real(fates_r8), parameter, public :: umol_per_mol = 1.0E6_fates_r8

  ! Conversion factor: moles per micromole
  real(fates_r8), parameter, public :: mol_per_umol = 1.0E-6_fates_r8
  
  ! Conversion factor: umols per kilomole
  real(fates_r8), parameter, public :: umol_per_kmol = 1.0E9_fates_r8

  ! Conversion factor: meters per milimeter
  real(fates_r8), parameter, public :: m_per_mm = 1.0E-3_fates_r8

  ! Conversion factor: milimeters per meter
  real(fates_r8), parameter, public :: mm_per_m = 1.0E3_fates_r8

  ! Conversion factor: millimeters per centimeter (ahb added this 7/7/2021)
  real(fates_r8), parameter, public :: mm_per_cm = 10.0_fates_r8
  
  ! Conversion factor: meters per centimeter
  real(fates_r8), parameter, public :: m_per_cm = 1.0E-2_fates_r8

  ! Conversion factor: m2 per ha
  real(fates_r8), parameter, public :: m2_per_ha = 1.0e4_fates_r8

  ! Conversion factor: m2 per km2
  real(fates_r8), parameter, public :: m2_per_km2 = 1.0e6_fates_r8

  ! Conversion factor: cm2 per m2
  real(fates_r8), parameter, public :: cm2_per_m2 = 10000.0_fates_r8

  ! Conversion factor: m3 per mm3
  real(fates_r8), parameter, public :: m3_per_mm3 = 1.0E-9_fates_r8

  ! Conversion factor: cubic meters per cubic cm
  real(fates_r8), parameter, public :: m3_per_cm3 = 1.0E-6_fates_r8

  real(fates_r8), parameter, public :: cm3_per_m3 = 1.0E6_fates_r8

  ! Conversion factor :: ha per m2
  real(fates_r8), parameter, public :: ha_per_m2 = 1.0e-4_fates_r8

  ! Conversion: seconds per minute
  real(fates_r8), parameter, public :: sec_per_min = 60.0_fates_r8

  ! Conversion: seconds per day
  real(fates_r8), parameter, public :: sec_per_day = 86400.0_fates_r8

  ! Conversion: megajoules per joule
  real(fates_r8), parameter, public :: megajoules_per_joule = 1.0E-6_fates_r8
 
  
  ! Conversion: days per second
  real(fates_r8), parameter, public :: days_per_sec = 1.0_fates_r8/86400.0_fates_r8

  ! Conversion: days per year. assume HLM uses 365 day calendar.
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter, public :: days_per_year = 365.00_fates_r8

  ! Integer version of days per year.
  integer, parameter, public :: ndays_per_year = nint(days_per_year)

  ! Conversion: years per day. assume HLM uses 365 day calendar.
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter, public :: years_per_day = 1.0_fates_r8/365.00_fates_r8

  ! Conversion: months per year
  real(fates_r8), parameter, public :: months_per_year = 12.0_fates_r8

  ! Conversion: Joules per kiloJoules
  real(fates_r8), parameter, public :: J_per_kJ = 1000.0_fates_r8

  ! Physical constants
  
  ! dewpoint calculation
  real(fates_r8), parameter, public :: dewpoint_a = 17.62_fates_r8
  real(fates_r8), parameter, public :: dewpoint_b = 243.12_fates_r8 ![degrees C]

  ! universal gas constant [J/K/kmol]
  real(fates_r8), parameter, public :: rgas_J_K_kmol          = 8314.4598_fates_r8

  ! universal gas constant [J/k/mol]
  real(fates_r8), parameter, public :: rgas_J_K_mol           = 8.3144598_fates_r8

  ! freezing point of water at 1 atm (K)
  real(fates_r8), parameter, public :: t_water_freeze_k_1atm   = 273.15_fates_r8

  ! freezing point of water at triple point (K)
  real(fates_r8), parameter, public :: t_water_freeze_k_triple = 273.16_fates_r8

  ! Density of fresh liquid water (kg/m3)
  real(fates_r8), parameter, public :: dens_fresh_liquid_water = 1.0E3_fates_r8

  ! Molar mass of water (g/mol)
  real(fates_r8), parameter, public :: molar_mass_water = 18.0_fates_r8

  ! Approximate molar mass of water vapor to dry air (-)
  real(fates_r8), parameter, public :: molar_mass_ratio_vapdry= 0.622_fates_r8
  
  ! Gravity constant on earth [m/s]
  real(fates_r8), parameter, public :: grav_earth = 9.8_fates_r8

  ! Megapascals to pascals
  real(fates_r8), parameter, public :: pa_per_mpa = 1.e6_fates_r8

  ! Pascals to megapascals
  real(fates_r8), parameter, public :: mpa_per_pa = 1.e-6_fates_r8

  ! Conversion: megapascals per mm H2O suction
  real(fates_r8), parameter, public :: mpa_per_mm_suction = dens_fresh_liquid_water * &
                                       grav_earth * 1.0E-9_fates_r8

  ! For numerical inquiry
  real(fates_r8), parameter, public :: fates_huge = huge(g_per_kg)

  real(fates_r8), parameter, public :: fates_tiny = tiny(g_per_kg)
  
  ! Geodesy constants (WGS 84)
  real(fates_r8), parameter, public :: earth_radius_eq = 6378137.0_fates_r8                       ! equitorial radius, earth [m]
  real(fates_r8), parameter, public :: earth_flattening = 1.0_fates_r8 / 298.257223563_fates_r8  ! flattening [non-dimensional]
  

  ! Geometric Constants

  ! PI
  real(fates_r8), parameter, public :: pi_const = 3.14159265359_fates_r8
  real(fates_r8), parameter, public :: rad_per_deg = pi_const/180.0_fates_r8

  ! Rdark constants from Atkin et al., 2017 https://doi.org/10.1007/978-3-319-68703-2_6
  ! and Heskel et al., 2016 https://doi.org/10.1073/pnas.1520282113
  real(fates_r8), parameter, public :: lmr_b = 0.1012_fates_r8       ! (degrees C**-1)

  real(fates_r8), parameter, public :: lmr_c = -0.0005_fates_r8      ! (degrees C**-2)

  real(fates_r8), parameter, public :: lmr_TrefC = 25._fates_r8      ! (degrees C)

  real(fates_r8), parameter, public :: lmr_r_1 = 0.2061_fates_r8     ! (umol CO2/m**2/s / (gN/(m2 leaf))) 

  real(fates_r8), parameter, public :: lmr_r_2 = -0.0402_fates_r8    ! (umol CO2/m**2/s/degree C)
  
  ! some integers related to termination mortality
  integer, parameter, public :: n_term_mort_types = 3
  integer, parameter, public :: i_term_mort_type_cstarv = 1
  integer, parameter, public :: i_term_mort_type_canlev = 2
  integer, parameter, public :: i_term_mort_type_numdens = 3

end module FatesConstantsMod
