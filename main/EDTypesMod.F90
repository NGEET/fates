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

  type, public :: site_fluxdiags_type

     ! ----------------------------------------------------------------------------------
     ! Diagnostics for fluxes into the litter pool from plants
     ! these fluxes are the total from 
     ! (1) turnover from living plants
     ! (2) mass transfer from non-disturbance inducing mortality events
     ! (3) mass transfer from disturbance inducing mortality events
     ! [kg / ha / day]
     ! ---------------------------------------------------------------------------------

     real(r8) :: cwd_ag_input(1:ncwd)               
     real(r8) :: cwd_bg_input(1:ncwd)               
     real(r8),allocatable :: leaf_litter_input(:)
     real(r8),allocatable :: root_litter_input(:)
     
   contains

     procedure :: ZeroFluxDiags
     
  end type site_fluxdiags_type

  ! ====================================================================================

  type, public ::  site_massbal_type

     ! ----------------------------------------------------------------------------------
     ! This type is used for accounting purposes to ensure that we are not
     ! loosing or creating mass. This type is supposed to be allocated for each element 
     ! we simulate (e.g. carbon12_element, etc)
     ! Note that the unit of "site", is nominally equivalent to 1 hectare
     !
     ! This set of mass checks are for INCREMENTAL checks during the dynamics step.
     ! ----------------------------------------------------------------------------------
     
     real(r8) :: old_stock    ! remember biomass stock from last time  [Kg/site]
     real(r8) :: err_fates    ! Total mass balance error for FATES processes     [kg/site]


     ! ----------------------------------------------------------------------------------
     ! Group 3: Components of the total site level mass fluxes
     ! ----------------------------------------------------------------------------------

     real(r8) :: gpp_acc          ! Accumulated gross primary productivity [kg/site/day]
     real(r8) :: aresp_acc        ! Accumulated autotrophic respiration [kg/site/day]

     real(r8) :: net_root_uptake  ! Net uptake of carbon or nutrients through the roots [kg/site/day]
                                  ! could include exudation, and for N this also includes symbiotic
                                  ! fixation

     real(r8) :: seed_in          ! Total mass of external seed rain into fates site [kg/site/day]
                                  ! This is from external grid-cells or from user parameterization
                                  ! (user param seed rain, or dispersal model)
     real(r8) :: seed_out         ! Total mass of seeds exported outside of fates site [kg/site/day]
                                  ! (this is not used currently, placeholder, rgk feb-2019)

     real(r8) :: frag_out         ! Litter and coarse woody debris fragmentation flux [kg/site/day]

     real(r8) :: wood_product          ! Total mass exported as wood product [kg/site/day]
     real(r8) :: burn_flux_to_atm      ! Total mass burned and exported to the atmosphere [kg/site/day]

     real(r8) :: flux_generic_in       ! Used for prescribed or artificial input fluxes
                                       ! and initialization [kg/site/day]
     real(r8) :: flux_generic_out      ! Used for prescribed or artificial output fluxes
                                       ! for instance when prescribed physiology is on
     real(r8) :: patch_resize_err      ! This is the amount of mass gained (or loss when negative)
                                       ! due to re-sizing patches when area math starts to lose
                                       ! precision

   contains

     procedure :: ZeroMassBalState
     procedure :: ZeroMassBalFlux
     
  end type site_massbal_type
  

  contains
      
    ! =====================================================================================

    subroutine ZeroFluxDiags(this)
      
      class(site_fluxdiags_type) :: this
      
      this%cwd_ag_input(:)      = 0._r8
      this%cwd_bg_input(:)      = 0._r8
      this%leaf_litter_input(:) = 0._r8
      this%root_litter_input(:) = 0._r8
      
      return
    end subroutine ZeroFluxDiags

    ! =====================================================================================
    
    subroutine ZeroMassBalState(this)
      
      class(site_massbal_type) :: this
      
      this%old_stock = 0._r8
      this%err_fates = 0._r8
      
      return
    end subroutine ZeroMassBalState
    
    subroutine ZeroMassBalFlux(this)
      
      class(site_massbal_type) :: this

      this%gpp_acc           = 0._r8
      this%aresp_acc         = 0._r8
      this%net_root_uptake   = 0._r8
      this%seed_in           = 0._r8
      this%seed_out          = 0._r8
      this%frag_out          = 0._r8
      this%wood_product      = 0._r8
      this%burn_flux_to_atm  = 0._r8
      this%flux_generic_in   = 0._r8
      this%flux_generic_out  = 0._r8
      this%patch_resize_err  = 0._r8

      return
  end subroutine ZeroMassBalFlux
   
  ! =====================================================================================
  
  
end module EDTypesMod
