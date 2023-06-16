module EDTypesMod

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesGlobals,          only : endrun => fates_endrun
  use FatesConstantsMod,     only : ifalse
  use FatesConstantsMod,     only : itrue
  use FatesGlobals,          only : fates_log
  use FatesHydraulicsMemMod, only : ed_cohort_hydr_type
  use FatesHydraulicsMemMod, only : ed_site_hydr_type
  use PRTGenericMod,         only : prt_vartypes
  use PRTGenericMod,         only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,         only : repro_organ, store_organ, struct_organ
  use PRTGenericMod,         only : prt_carbon_allom_hyp
  use PRTGenericMod,         only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,         only : num_organ_types
  use PRTGenericMod,         only : num_elements
  use PRTGenericMod,         only : element_list
  use PRTGenericMod,         only : num_element_types
  use PRTGenericMod,         only : carbon12_element
  use FatesLitterMod,        only : litter_type
  use FatesLitterMod,        only : ncwd, NFSC
  use FatesConstantsMod,     only : n_anthro_disturbance_categories
  use FatesConstantsMod,     only : days_per_year
  use FatesConstantsMod,     only : fates_unset_r8
  use FatesRunningMeanMod,   only : rmean_type
  use FatesInterfaceTypesMod,only : bc_in_type
  use FatesInterfaceTypesMod,only : bc_out_type
  use FatesInterfaceTypesMod,only : hlm_parteh_mode
  use FatesCohortMod,        only : fates_cohort_type
  use FatesPatchMod,         only : fates_patch_type
  use EDParamsMod,           only : maxSWb, nclmax, nlevleaf
  use FatesConstantsMod,     only : maxpft
  use FatesConstantsMod,     only : n_dbh_bins, n_dist_types
  use shr_log_mod,           only : errMsg => shr_log_errMsg

  implicit none
  private               ! By default everything is private
  save
              
  real(r8), parameter, public :: init_recruit_trim = 0.8_r8    ! This is the initial trimming value that
                                                               ! new recruits start with

  ! -------------------------------------------------------------------------------------
  ! Radiation parameters
  ! These should be part of the radiation module, but since we only have one option
  ! this is ok for now. (RGK 04-2018)
  ! -------------------------------------------------------------------------------------

  integer, parameter, public :: n_rad_stream_types = 2    ! The number of radiation streams used (direct/diffuse)
 
  integer, parameter, public :: idirect   = 1             ! This is the array index for direct radiation
  integer, parameter, public :: idiffuse  = 2             ! This is the array index for diffuse radiation


 ! Flag to turn on/off salinity effects on the effective "btran"
  ! btran stress function.

  logical, parameter, public :: do_fates_salinity = .false.


  ! This is the community level amount of spread expected in nearly-bare-ground
  ! and inventory starting modes.
  ! These are used to initialize only. These values will scale between
  ! the PFT defined maximum and minimum crown area scaing parameters.
  !
  ! A value of 1 indicates that
  ! plants should have crown areas at maximum spread for their size and PFT.
  ! A value of 0 means that they have the least amount of spread for their
  ! size and PFT.
  
  real(r8), parameter, public :: init_spread_near_bare_ground = 1.0_r8
  real(r8), parameter, public :: init_spread_inventory        = 0.0_r8


  ! MODEL PARAMETERS

  real(r8), parameter, public :: area                 = 10000.0_r8 ! Notional area of simulated forest m2
  real(r8), parameter, public :: area_inv             = 1.0e-4_r8  ! Inverse of the notion area (faster math)

  integer, parameter, public  :: numWaterMem          = 10         ! watermemory saved as site level var

  integer, parameter, public  :: numlevsoil_max       = 30         ! This is scratch space used for static arrays
                                                                   ! The actual number of soil layers should not exceed this


  ! BIOLOGY/BIOGEOCHEMISTRY        
  integer , parameter, public :: num_vegtemp_mem      = 10         ! Window of time over which we track temp for cold sensecence (days)
  integer , parameter, public :: dtype_ifall          = 1          ! index for naturally occuring tree-fall generated event
  integer , parameter, public :: dtype_ifire          = 2          ! index for fire generated disturbance event
  integer , parameter, public :: dtype_ilog           = 3          ! index for logging generated disturbance event

  ! Phenology status flag definitions (cold type is cstat, dry type is dstat)

  integer, parameter, public :: phen_cstat_nevercold = 0        ! This (location/plant) has not experienced a cold period over a large number
                                                        ! of days, leaves are dropped and flagged as non-cold region
  integer, parameter, public :: phen_cstat_iscold    = 1        ! This (location/plant) is in a cold-state where leaves should have fallen
  integer, parameter, public :: phen_cstat_notcold   = 2        ! This site is in a warm-state where leaves are allowed to flush

  integer, parameter, public :: phen_dstat_timeoff   = 0       ! Leaves off due to time exceedance (drought phenology)
  integer, parameter, public :: phen_dstat_moistoff  = 1       ! Leaves off due to moisture avail  (drought phenology)
  integer, parameter, public :: phen_dstat_moiston   = 2       ! Leaves on due to moisture avail   (drought phenology)
  integer, parameter, public :: phen_dstat_timeon    = 3       ! Leaves on due to time exceedance  (drought phenology)

  ! PATCH FUSION 
  real(r8), parameter, public :: force_patchfuse_min_biomass = 0.005_r8   ! min biomass (kg / m2 patch area) below which to force-fuse patches
  real(r8), parameter, public :: patchfusion_dbhbin_loweredges(N_DBH_BINS) = &
       (/0._r8, 5._r8, 20._r8, 50._r8, 100._r8, 150._r8/)                 ! array of bin lower edges for comparing patches
  real(r8), parameter, public :: patch_fusion_tolerance_relaxation_increment = 1.1_r8 ! amount by which to increment patch fusion threshold
  real(r8), parameter, public :: max_age_of_second_oldest_patch = 200._r8 ! age in years above which to combine all patches

  ! COHORT FUSION
  real(r8), parameter, public :: HITEMAX              = 30.0_r8    ! max dbh value used in hgt profile comparison 
  integer , parameter, public :: N_HITE_BINS          = 60         ! no. of hite bins used to distribute LAI

  ! COHORT TERMINATION

  real(r8), parameter, public :: min_npm2       = 1.0E-7_r8               ! minimum cohort number density per m2 before termination
  real(r8), parameter, public :: min_patch_area = 0.01_r8                 ! smallest allowable patch area before termination
  real(r8), parameter, public :: min_patch_area_forced = 0.0001_r8        ! patch termination will not fuse the youngest patch
                                                                          ! if the area is less than min_patch_area.
                                                                          ! however, it is allowed to fuse the youngest patch
                                                                          ! if the fusion area is less than min_patch_area_forced

  real(r8), parameter, public :: min_nppatch    = min_npm2*min_patch_area ! minimum number of cohorts per patch (min_npm2*min_patch_area)
  real(r8), parameter, public :: min_n_safemath = 1.0E-12_r8              ! in some cases, we want to immediately remove super small
                                                                          ! number densities of cohorts to prevent FPEs

  ! special mode to cause PFTs to create seed mass of all currently-existing PFTs
  logical, parameter, public :: homogenize_seed_pfts  = .false.
  character(len=*), parameter, private :: sourcefile = __FILE__

  !************************************
  !** Resources management type      **
  ! YX
  !************************************
  type, public :: ed_resources_management_type
    
     real(r8) ::  trunk_product_site                       ! Actual  trunk product at site level KgC/site
     real(r8) ::  harvest_debt                             ! the amount of kgC per site that did not successfully harvested 
     real(r8) ::  harvest_debt_sec                         ! the amount of kgC per site from secondary patches that did
                                                           ! not successfully harvested

     !debug variables
     real(r8) ::  delta_litter_stock                       ! kgC/site = kgC/ha
     real(r8) ::  delta_biomass_stock                      ! kgC/site
     real(r8) ::  delta_individual                         ! 
  
  end type ed_resources_management_type

  ! =====================================================================================

end module EDTypesMod
