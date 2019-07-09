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

  ! Integer equivalent of true  (in case some compilers dont auto convert)
  integer, parameter, public :: itrue = 1

  ! Integer equivalent of false (in case come compilers dont auto convert)
  integer, parameter, public :: ifalse = 0

  ! Labels for patch disturbance history
  integer, parameter, public :: n_anthro_disturbance_categories = 2
  integer, parameter, public :: primaryforest = 1
  integer, parameter, public :: secondaryforest = 2
  
  ! Error Tolerances

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

  ! Unit conversion constants:

  ! Conversion factor umols of Carbon -> kg of Carbon (1 mol = 12g)
  ! We do not use umolC_per_kg because it is a non-terminating decimal
  real(fates_r8), parameter, public :: umolC_to_kgC = 12.0E-9_fates_r8
  
  ! Conversion factor: grams per kilograms
  real(fates_r8), parameter, public :: g_per_kg = 1000.0_fates_r8
  
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

  ! Conversion factor: umols per kilomole
  real(fates_r8), parameter, public :: umol_per_kmol = 1.0E9_fates_r8

  ! Conversion factor: m2 per ha
  real(fates_r8), parameter, public :: m2_per_ha = 1.0e4_fates_r8

  ! Conversion factor: cm2 per m2
  real(fates_r8), parameter, public :: cm2_per_m2 = 10000.0_fates_r8

  ! Conversion factor :: ha per m2
  real(fates_r8), parameter, public :: ha_per_m2 = 1.0e-4_fates_r8

  ! Conversion: seconds per minute
  real(fates_r8), parameter, public :: sec_per_min = 60.0_fates_r8

  ! Conversion: seconds per day
  real(fates_r8), parameter, public :: sec_per_day = 86400.0_fates_r8

  ! Conversion: days per second
  real(fates_r8), parameter, public :: days_per_sec = 1.0_fates_r8/86400.0_fates_r8

  ! Conversion: days per year. assume HLM uses 365 day calendar.  
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter, public :: days_per_year = 365.00_fates_r8
  
  ! Conversion: years per day. assume HLM uses 365 day calendar.  
  ! If we need to link to 365.25-day-calendared HLM, rewire to pass through interface
  real(fates_r8), parameter, public :: years_per_day = 1.0_fates_r8/365.00_fates_r8
  
  ! Physical constants

  ! universal gas constant [J/K/kmol]
  real(fates_r8), parameter, public :: rgas_J_K_kmol          = 8314.4598_fates_r8

  ! freezing point of water at 1 atm (K)
  real(fates_r8), parameter, public :: t_water_freeze_k_1atm   = 273.15_fates_r8     

  ! freezing point of water at triple point (K)
  real(fates_r8), parameter, public :: t_water_freeze_k_triple = 273.16_fates_r8      

  ! Density of fresh liquid water (kg/m3)
  real(fates_r8), parameter, public :: dens_fresh_liquid_water = 1.0E3_fates_r8

  ! Gravity constant on earth [m/s]
  real(fates_r8), parameter, public :: grav_earth = 9.8_fates_r8


  ! For numerical inquiry
  real(fates_r8), parameter, public :: fates_huge = huge(g_per_kg)

  real(fates_r8), parameter, public :: fates_tiny = tiny(g_per_kg)

  ! Geometric Constants
  
  ! PI
  real(fates_r8), parameter, public :: pi_const = 3.14159265359_fates_r8

end module FatesConstantsMod
