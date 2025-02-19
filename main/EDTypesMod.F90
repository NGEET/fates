module EDTypesMod

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesGlobals,          only : endrun => fates_endrun
  use FatesConstantsMod,     only : ifalse
  use FatesConstantsMod,     only : itrue
  use FatesConstantsMod,     only : nocomp_bareground_land
  use FatesConstantsMod,     only : nocomp_bareground
  use FatesConstantsMod,     only : secondaryland
  use FatesConstantsMod,     only : secondary_age_threshold
  use FatesConstantsMod,     only : nearzero
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
  use FatesLitterMod,        only : ncwd
  use FatesConstantsMod,     only : days_per_year
  use FatesRunningMeanMod,   only : rmean_type,rmean_arr_type
  use FatesConstantsMod,     only : fates_unset_r8
  use FatesInterfaceTypesMod,only : bc_in_type
  use FatesInterfaceTypesMod,only : bc_out_type
  use FatesConstantsMod    , only : n_landuse_cats
  use FatesInterfaceTypesMod,only : hlm_parteh_mode
  use FatesCohortMod,        only : fates_cohort_type
  use FatesPatchMod,         only : fates_patch_type
  use EDParamsMod,           only : nclmax, nlevleaf, maxpft
  use FatesConstantsMod,     only : n_dbh_bins, n_dist_types
  use shr_log_mod,           only : errMsg => shr_log_errMsg
  use SFFireWeatherMod,      only : fire_weather

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

  ! Phenology status flag definitions (cold type is cstat, dry type is dstat)

  integer, parameter, public :: phen_cstat_nevercold = 0        ! This (location/plant) has not experienced a cold period over a large number
                                                        ! of days, leaves are dropped and flagged as non-cold region
  integer, parameter, public :: phen_cstat_iscold    = 1        ! This (location/plant) is in a cold-state where leaves should have fallen
  integer, parameter, public :: phen_cstat_notcold   = 2        ! This site is in a warm-state where leaves are allowed to flush

  integer, parameter, public :: phen_dstat_timeoff   = 0       ! Leaves off due to time exceedance (drought phenology)
  integer, parameter, public :: phen_dstat_moistoff  = 1       ! Leaves off due to moisture avail  (drought phenology)
  integer, parameter, public :: phen_dstat_moiston   = 2       ! Leaves on due to moisture avail   (drought phenology)
  integer, parameter, public :: phen_dstat_timeon    = 3       ! Leaves on due to time exceedance  (drought phenology)
  integer, parameter, public :: phen_dstat_pshed     = 4 ! Leaves partially abscissing       (drought phenology)

  ! PATCH FUSION 
  real(r8), parameter, public :: force_patchfuse_min_biomass = 0.005_r8   ! min biomass (kg / m2 patch area) below which to force-fuse patches
  real(r8), parameter, public :: patch_fusion_tolerance_relaxation_increment = 1.1_r8 ! amount by which to increment patch fusion threshold
  real(r8), parameter, public :: max_age_of_second_oldest_patch = 200._r8 ! age in years above which to combine all patches

  ! COHORT FUSION
  real(r8), parameter, public :: HEIGHTMAX           = 30.0_r8    ! max dbh value used in hgt profile comparison 
  integer , parameter, public :: N_HEIGHT_BINS       = 60         ! no. of height bins used to distribute LAI

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

  type, public :: elem_diag_type

     ! ----------------------------------------------------------------------------------
     ! Diagnostics of fluxes
     ! These act as an intermediary to write fluxes to the history
     ! file after number densities of plants have changed. They also
     ! allow the history flux diagnostics to be rebuilt during restart
     !
     ! Litter fluxes are the total from 
     ! (1) turnover from living plants
     ! (2) mass transfer from non-disturbance inducing mortality events
     ! (3) mass transfer from disturbance inducing mortality events
     ! [kg / ha / day]
     ! ---------------------------------------------------------------------------------
     
     real(r8) :: cwd_ag_input(1:ncwd)               
     real(r8) :: cwd_bg_input(1:ncwd)               
     real(r8),allocatable :: surf_fine_litter_input(:)
     real(r8),allocatable :: root_litter_input(:)

     real(r8) :: tot_seed_turnover ! decay of living seed bank to
                                   ! fragmented litter [kg/m2/day]
     real(r8) :: exported_harvest  ! mass of harvested vegetation exported and not sent to litter [kg/m2/day]
     real(r8) :: burned_liveveg    ! Amount of mass burned from living plants [kg/m2/day]

     ! Integrated Error Terms ( Int. Flux - State ) 
     
     real(r8) :: err_liveveg       ! Error from comparing [state-integrated flux]
                                   ! in live vegetation [kg/m2]
     real(r8) :: err_litter        ! Net change in litter [kg/m2]
     
     
  end type elem_diag_type

  ! -------------------------------------------------------------------------

  type, public :: site_ifluxbal_type

     ! The combination of living vegetation and litter accounts
     ! for all of the mass that is tracked by FATES. We use these
     ! data structures to ensure that an instantaneous assessment
     ! of the mass of live vegetation and litter, is the same
     ! as the initial condition plus the integrated fluxes in and
     ! out of those pools over the duration of the simulation.
     
     ! Mass in living vegetation, this includes:
     ! All organs on living plants, including non respiring tissues
     ! and heartwood.
     ! The live seed pool
     
     real(r8) :: state_liveveg   ! Assessed instanteously [kg/m2]
     real(r8) :: iflux_liveveg   ! Integrated daily       [kg/m2]

     ! Mass in all litter tracked by FATES
     ! This includes only the unfragmented litter, litter that
     ! has fragmented is tracked by the host models
     
     real(r8) :: state_litter    ! Assessed instantaneously [kg/m2]
     real(r8) :: iflux_litter    ! Integrated daily         [kg/m2]

  end type site_ifluxbal_type

  ! -------------------------------------------------------------------------
  
  type, public :: site_fluxdiags_type


     ! This is for all diagnostics that are uniform over all elements (C,N,P)
     
     type(elem_diag_type), pointer :: elem(:)

     ! This variable is slated as to-do, but the fluxdiags type needs
     ! to be refactored first. Currently this type is allocated
     ! by chemical species (ie C, N or P). GPP is C, but not N or P (RGK 0524)
     ! Previous day GPP [kgC/m2/year], partitioned by size x pft
     !real(r8),allocatable :: gpp_prev_scpf(:)

     real(r8) :: npp          ! kg m-2 day-1
     
     ! Nutrient Flux Diagnostics
     
     real(r8) :: resp_excess  ! plant carbon respired due to carbon overflow
                              ! this happens when nutrients are limiting construction
                              ! of new tissues kg m-2 s-1
     real(r8) :: nh4_uptake   ! plant nh4 uptake, kg m-2 s-1
     real(r8) :: no3_uptake   ! plant no3 uptake, kg m-2 s-1
     real(r8) :: sym_nfix     ! plant N uptake via symbiotic fixation kg m-2 s-1
     real(r8) :: n_efflux     ! efflux of unusable N from plant to soil labile pool kg m-2 s-1
     real(r8) :: p_uptake     ! po4 uptake, kg m-2 s-1
     real(r8) :: p_efflux     ! efflux of unusable P from plant to soil labile pool kg m-2 s-1

     ! Size by PFT delineated nutrient flux diagnostics (same units as above)
     ! These are only allocated if both complex history diagnostics, and the
     ! species of interest (N or P) is requested
     
     real(r8),allocatable :: nh4_uptake_scpf(:)
     real(r8),allocatable :: no3_uptake_scpf(:)
     real(r8),allocatable :: sym_nfix_scpf(:)
     real(r8),allocatable :: n_efflux_scpf(:)
     real(r8),allocatable :: p_uptake_scpf(:)
     real(r8),allocatable :: p_efflux_scpf(:)


     
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

     real(r8) :: wood_product_harvest(maxpft)    ! Total mass exported as wood product from wood harvest [kg/site/day]

     real(r8) :: wood_product_landusechange(maxpft)    ! Total mass exported as wood product from land use change [kg/site/day]

     real(r8) :: burn_flux_to_atm      ! Total mass burned and exported to the atmosphere [kg/site/day]

     real(r8) :: flux_generic_in       ! Used for prescribed or artificial input fluxes
                                       ! and initialization [kg/site/day]
     real(r8) :: flux_generic_out      ! Used for prescribed or artificial output fluxes
                                       ! for instance when prescribed physiology is on
     real(r8) :: patch_resize_err      ! This is the amount of mass gained (or loss when negative)
                                       ! due to re-sizing patches when area math starts to lose
                                       ! precision

     real(r8) :: herbivory_flux_out    ! loss of element due to grazing (and/or browsing) by herbivores
     
   contains

     procedure :: ZeroMassBalState
     procedure :: ZeroMassBalFlux
     
  end type site_massbal_type
  

  !************************************
  !** Site type structure           **
  !************************************

  type, public :: ed_site_type
     
     ! POINTERS  
     type (fates_patch_type), pointer :: oldest_patch => null()   ! pointer to oldest patch at the site  
     type (fates_patch_type), pointer :: youngest_patch => null() ! pointer to yngest patch at the site
     
     ! Resource management
     type (ed_resources_management_type) :: resources_management ! resources_management at the site 

     ! If this simulation uses shared memory then the sites need to know what machine
     ! index they are on. This index is (currently) only used to identify the sites
     ! position in history output fields
     !integer :: clump_id 

     ! Global index of this site in the history output file
     integer :: h_gid
     
     ! INDICES 
     real(r8) ::  lat                                          ! latitude:  degrees 
     real(r8) ::  lon                                          ! longitude: degrees 

     ! Fixed Biogeography mode inputs
     real(r8), allocatable :: area_PFT(:,:)                    ! Area allocated to individual PFTs, indexed by land use class  [ha/ha of non-bareground area]
     real(r8) :: area_bareground                               ! Area allocated to bare ground in nocomp configurations (corresponds to HLM PFT 0) [ha/ha]

     integer, allocatable  :: use_this_pft(:)                  ! Is area_PFT > 0 ? (1=yes, 0=no)

     ! Total area of patches in each age bin [m2]
     real(r8), allocatable :: area_by_age(:)

     
     ! Nutrient relevant 
     real(r8), allocatable :: rec_l2fr(:,:) ! A running mean of the l2fr's for the newly
                                            ! recruited, pft x canopy_layer
     real(r8) :: ema_npp                    ! An exponential moving average of NPP [gC/m2/year]
                                            ! The lengthscale is hard-coded "ema_npp_tcale"
                                            ! in FatesSoilBGCFluxMod. Used solely to inform bc_out%ema_npp
                                            ! which is used for fixation

     
     ! Two-stream scratch arrays
     real(r8), allocatable :: omega_2str(:,:)   ! This is the matrix that is inverted to solve
                                                ! the linear system of equations in the two-stream
                                                ! radiation module. This array will grow
                                                ! and shrink depending on how many scattering
                                                ! elements there are. This matrix is square,
                                                ! and needs to be larger than 2 x number-of-elements
                                                ! for each patch on the site

     real(r8), allocatable :: taulambda_2str(:) ! These are the coefficients of the two-stream
                                                ! linear system of equations (ie the unknowns, "lambda")
                                                ! As well as the left-side (constants, "tau"). Since
                                                ! the LAPACK solver dgesv uses the latter
                                                ! as the argument and over-writes, we only
                                                ! need one array

     integer, allocatable :: ipiv_2str(:)       ! pivot indices for the lapack 2str solver

     ! SP mode target PFT level variables
     real(r8), allocatable :: sp_tlai(:)                      ! target TLAI per FATES pft
     real(r8), allocatable :: sp_tsai(:)                      ! target TSAI per FATES pft
     real(r8), allocatable :: sp_htop(:)                      ! target HTOP per FATES pft
     
     ! Instantaneous Mass Balance (allocation for each element)

     type(site_massbal_type), pointer :: mass_balance(:)

     ! Integrated Mass Balance checks, i.e. do the integrated
     ! fluxes match the state?  One is for live vegetation,
     ! one if for litter

     type(site_ifluxbal_type), pointer :: iflux_balance(:)
     
     
     ! Flux diagnostics
     ! These are used to write history output based on fluxes
     ! that are calculated before mortality and cohort/patch restructuring
     ! but are written after patch structure. This structure also allows for
     ! accurate restarting of the history      

     type(site_fluxdiags_type) :: flux_diags
     

     ! PHENOLOGY 
     real(r8) ::  grow_deg_days                                ! Phenology growing degree days
     real(r8) ::  snow_depth                                   ! site-level snow depth (used for ELAI/TLAI calcs)

     integer  ::  cstatus                                      ! are leaves in this pixel on or off for cold decid
                                                               ! 0 = this site has not experienced a cold period over at least
                                                               !     400 days, leaves are dropped and flagged as non-cold region
                                                               ! 1 = this site is in a cold-state where leaves should have fallen
                                                               ! 2 = this site is in a warm-state where leaves are allowed to flush
     integer  ::  dstatus(maxpft)                              ! are leaves in this pixel on or off for drought decid
                                                               ! 0 = leaves off due to time exceedance
                                                               ! 1 = leaves off due to moisture avail
                                                               ! 2 = leaves on due to moisture avail
                                                               ! 3 = leaves on due to time exceedance
                                                               ! 4 = leaves partially on (ED2-like phenology)
     integer  ::  nchilldays                                   ! num chilling days: (for botta gdd trheshold calculation)
     integer  ::  ncolddays                                    ! num cold days: (must exceed threshold to drop leaves)
     real(r8) ::  vegtemp_memory(num_vegtemp_mem)              ! record of last 10 days temperature for senescence model. deg C
     integer  ::  cleafondate                                  ! model date (day integer) of leaf on (cold):-
     integer  ::  cleafoffdate                                 ! model date (day integer) of leaf off (cold):-
     integer  ::  cndaysleafon                                 ! number of days since leaf on period started (cold)
     integer  ::  cndaysleafoff                                ! number of days since leaf off period started (cold)
     integer  ::  dleafondate(maxpft)                          ! model date (day integer) of leaf on drought:-
     integer  ::  dleafoffdate(maxpft)                         ! model date (day integer) of leaf off drought:-
     integer  ::  dndaysleafon(maxpft)                         ! number of days since leaf on period started (drought)
     integer  ::  dndaysleafoff(maxpft)                        ! number of days since leaf off period started (drought)
     real(r8) ::  elong_factor(maxpft)                         ! Elongation factor (ED2-like phenology). This is zero when leaves are
                                                               ! completely off, and one when they are completely flushed.
     integer  ::  phen_model_date                              ! current model date (day integer)
                                                               ! this date stays continuous when
                                                               ! in runs that are restarted, regardless of
                                                               ! the conditions of restart

     real(r8) ::  liqvol_memory(numWaterMem,maxpft)            ! last 10 days of soil liquid water volume (drought phenology)
     real(r8) ::  smp_memory(numWaterMem,maxpft)               ! last 10 days of soil matric potential (drought phenology)


     ! FIRE
     real(r8) ::  wind                                         ! daily wind in m/min for Spitfire units 
     real(r8) ::  fdi                                          ! daily probability an ignition event will start a fire
     real(r8) ::  NF                                           ! daily ignitions in km2
     real(r8) ::  NF_successful                                ! daily ignitions in km2 that actually lead to fire
     class(fire_weather), pointer :: fireWeather               ! fire weather object

     ! PLANT HYDRAULICS
     type(ed_site_hydr_type), pointer :: si_hydr

     ! Soil Layering

     integer :: nlevsoil                      ! Number of soil layers in this site
     real(r8), allocatable :: zi_soil(:)      ! interface level below a "z" level (m)
                                              ! this contains a zero index for surface.
     real(r8), allocatable :: dz_soil(:)      ! layer thickness (m)
     real(r8), allocatable :: z_soil(:)       ! layer depth (m)
     real(r8), allocatable :: rootfrac_scr(:) ! This is just allocated scratch space to hold
                                              ! root fractions. Since root fractions may be dependent
                                              ! on cohort properties, and we do not want to store this infromation
                                              ! on each cohort, we do not keep root fractions in
                                              ! memory, and instead calculate them on demand.
                                              ! This array is allocated over the number of soil
                                              ! layers for each site, and save allocating deallocating.
                                              ! NOTE: THIS SCRATCH SPACE WOULD NOT BE THREAD-SAFE
                                              ! IF WE FORK ON PATCHES


     ! Mineralized nutrient flux from veg to the soil, via multiple mechanisms
     ! inluding symbiotic fixation, or other 

     !real(r8) :: allocatable :: minn_flux_out  ! kg/ha/day
     !real(r8) :: allocatable :: minp_flux_out  ! kg/ha/day

     
     ! DIAGNOSTICS

     ! TERMINATION, RECRUITMENT, DEMOTION, and DISTURBANCE
     
     real(r8) :: term_crownarea_canopy                 ! crownarea from termination mortality, per canopy level
     real(r8) :: term_crownarea_ustory                 ! crownarea from termination mortality, per canopy level    
 
     real(r8) :: imort_crownarea                       ! crownarea of individuals killed due to impact mortality per year. [m2 day]

     real(r8) :: fmort_crownarea_canopy                ! crownarea of canopy indivs killed due to fire per year. [m2/sec]
     real(r8) :: fmort_crownarea_ustory                ! crownarea of understory indivs killed due to fire per year [m2/sec] 

     real(r8), allocatable :: term_nindivs_canopy(:,:,:)   ! number of canopy individuals that were in cohorts which 
                                                           ! were terminated this timestep, by termination type, size x pft
     real(r8), allocatable :: term_nindivs_ustory(:,:,:)   ! number of understory individuals that were in cohorts which 
                                                           ! were terminated this timestep, by termination type, size x pft

     real(r8), allocatable :: term_carbonflux_canopy(:,:)  ! carbon flux from live to dead pools associated 
                                                           ! with termination mortality, by termination type and pft. [kgC/ha/day]
     real(r8), allocatable :: term_carbonflux_ustory(:,:)  ! carbon flux from live to dead pools associated 
                                                         ! with termination mortality, by termination type and pft.  [kgC/ha/day]    
     real(r8), allocatable :: imort_carbonflux(:)        ! biomass of individuals killed due to impact mortality per year, by pft. [kgC/m2/sec]
     real(r8), allocatable :: fmort_carbonflux_canopy(:) ! biomass of canopy indivs killed due to fire per year, by pft. [gC/m2/sec]
     real(r8), allocatable :: fmort_carbonflux_ustory(:) ! biomass of understory indivs killed due to fire per year, by pft [gC/m2/sec] 

     real(r8), allocatable :: term_abg_flux(:,:)          ! aboveground biomass lost due to termination mortality x size x pft
     real(r8), allocatable :: imort_abg_flux(:,:)         ! aboveground biomass lost due to impact mortality x size x pft [kgC/m2/sec]
     real(r8), allocatable :: fmort_abg_flux(:,:)         ! aboveground biomass lost due to fire mortality x size x pft


     real(r8) :: demotion_carbonflux                     ! biomass of demoted individuals from canopy to understory [kgC/ha/day]
     real(r8) :: promotion_carbonflux                    ! biomass of promoted individuals from understory to canopy [kgC/ha/day]
     real(r8) :: recruitment_rate(1:maxpft)              ! number of individuals that were recruited into new cohorts
     real(r8), allocatable :: demotion_rate(:)           ! rate of individuals demoted from canopy to understory per FATES timestep
     real(r8), allocatable :: promotion_rate(:)          ! rate of individuals promoted from understory to canopy per FATES timestep
     real(r8), allocatable :: imort_rate(:,:)            ! rate of individuals killed due to impact mortality per year.  on size x pft array
     

     real(r8), allocatable :: fmort_rate_canopy(:,:)     ! rate of canopy individuals killed due to fire mortality per year.  
                                                         ! on size x pft array  (1:nlevsclass,1:numpft)
     real(r8), allocatable :: fmort_rate_ustory(:,:)     ! rate of understory individuals killed due to fire mortality per year.  
                                                         ! on size x pft array  (1:nlevsclass,1:numpft)
    
     real(r8), allocatable :: fmort_rate_cambial(:,:)    ! rate of individuals killed due to fire mortality 
                                                         ! from cambial damage per year.  on size x pft array
     real(r8), allocatable :: fmort_rate_crown(:,:)      ! rate of individuals killed due to fire mortality 
                                                         ! from crown damage per year.  on size x pft array

     real(r8), allocatable :: imort_rate_damage(:,:,:)     ! number of individuals per damage class that die from impact mortality
     real(r8), allocatable :: term_nindivs_canopy_damage(:,:,:) ! number of individuals per damage class that die from termination mortality - canopy
     real(r8), allocatable :: term_nindivs_ustory_damage(:,:,:) ! number of individuals per damage class that die from termination mortality - canopy
     real(r8), allocatable :: fmort_rate_canopy_damage(:,:,:) ! number of individuals per damage class that die from fire - canopy
     real(r8), allocatable :: fmort_rate_ustory_damage(:,:,:) ! number of individuals per damage class that die from fire - ustory
     real(r8), allocatable :: fmort_cflux_canopy_damage(:,:) ! cflux per damage class that die from fire - canopy
     real(r8), allocatable :: fmort_cflux_ustory_damage(:,:) ! cflux per damage class that die from fire - ustory
     real(r8), allocatable :: imort_cflux_damage(:,:)         ! carbon flux from impact mortality by damage class [kgC/m2/sec]
     real(r8), allocatable :: term_cflux_canopy_damage(:,:)          ! carbon flux from termination mortality by damage class
     real(r8), allocatable :: term_cflux_ustory_damage(:,:)          ! carbon flux from termination mortality by damage class

     real(r8), allocatable :: growthflux_fusion(:,:)     ! rate of individuals moving into a given size class bin
                                                         ! due to fusion in a given day. on size x pft array 


     real(r8)              :: crownarea_canopy_damage                 ! crown area of canopy that is damaged annually  
     real(r8)              :: crownarea_ustory_damage                 ! crown area of understory that is damaged annually
     
     ! Canopy Spread
     real(r8) ::  spread                                          ! dynamic canopy allometric term [unitless]

     ! Seed dispersal
      real(r8), allocatable :: seed_out(:)                               ! amount of seed leaving the site [kg/site/day]
      real(r8), allocatable :: seed_in(:)                                ! amount of seed dispersed into the site from neighbouring cells  [kg/site/day]

     ! site-level variables to keep track of the disturbance rates, both actual and "potential"
     real(r8) :: disturbance_rates(N_DIST_TYPES,n_landuse_cats, n_landuse_cats)  ! actual disturbance rates for each disturbance type  [m2/m2/day]
     real(r8) :: primary_land_patchfusion_error             ! error term in total area of primary patches associated with patch fusion [m2/m2/day]
     real(r8) :: landuse_transition_matrix(n_landuse_cats, n_landuse_cats) ! land use transition matrix as read in from HLM and aggregated to FATES land use types [m2/m2/year]

     real(r8) :: min_allowed_landuse_fraction             ! minimum amount of land-use type below which the resulting patches would be too small [m2/m2]
     logical, allocatable :: landuse_vector_gt_min(:)     ! is the land use state vector for each land use type greater than the minimum below which we ignore?
     logical :: transition_landuse_from_off_to_on         ! special flag to use only when reading restarts, which triggers procedure to initialize land use

     contains

       procedure, public :: get_current_landuse_statevector
       procedure, public :: get_secondary_young_fraction

  end type ed_site_type
  
  ! Make public necessary subroutines and functions
  public :: dump_site
  public :: CalculateTreeGrassAreaSite
  public :: set_patchno
  
contains

  ! ============================================================================

  subroutine set_patchno( currentSite, check , call_id)

    !
    ! !DESCRIPTION:
    ! Give patches an order number from the oldest to youngest. 
    ! Oldest patches start with an index of 1.
    ! Special case: For no-comp runs, we treat the bare-ground
    ! patch as index 0.
    
    type(ed_site_type),intent(in) :: currentSite
    logical,intent(in) :: check     ! If true, we are checking order, not setting
    integer,intent(in) :: call_id   ! An index used for testing
    type(fates_patch_type), pointer :: currentPatch
    integer patchno

    !---------------------------------------------------------------------
    
    patchno = 1
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       if(currentPatch%nocomp_pft_label.eq.nocomp_bareground)then
          ! for bareground patch, we make the patch number 0
          if(check .and. currentPatch%patchno.ne.0)then
             write(fates_log(),*)'nocomp patch numbering is not correct:',currentPatch%patchno,'call_id:',call_id
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          currentPatch%patchno = 0
       else
          if(check .and. currentPatch%patchno.ne.patchno) then
             write(fates_log(),*)'patch numbering is not correct:',currentPatch%patchno,patchno,'call_id:',call_id
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          currentPatch%patchno = patchno
          patchno = patchno + 1
       endif
       currentPatch => currentPatch%younger
    enddo
    
    return
  end subroutine set_patchno
    
  ! =====================================================================================

    subroutine ZeroFluxDiags(this)
      
      class(site_fluxdiags_type) :: this
      integer :: el
      
      do el = 1,num_elements
         this%elem(el)%cwd_ag_input(:)      = 0._r8
         this%elem(el)%cwd_bg_input(:)      = 0._r8
         this%elem(el)%surf_fine_litter_input(:) = 0._r8
         this%elem(el)%root_litter_input(:) = 0._r8
         this%elem(el)%burned_liveveg       = 0._r8
         this%elem(el)%tot_seed_turnover   = 0._r8
         this%elem(el)%exported_harvest   = 0._r8
         this%elem(el)%err_liveveg   = 0._r8
         this%elem(el)%err_litter   = 0._r8

      end do

     this%npp = 0._r8
     this%resp_excess = 0._r8
     this%nh4_uptake  = 0._r8
     this%no3_uptake  = 0._r8
     this%sym_nfix = 0._r8
     this%n_efflux = 0._r8
     this%p_uptake = 0._r8
     this%p_efflux = 0._r8
     
     this%nh4_uptake_scpf(:) = 0._r8
     this%no3_uptake_scpf(:) = 0._r8
     this%sym_nfix_scpf(:) = 0._r8
     this%n_efflux_scpf(:) = 0._r8
     this%p_uptake_scpf(:) = 0._r8
     this%p_efflux_scpf(:) = 0._r8
      
     ! We don't zero gpp_prev_scpf because this is not
     ! incremented like others, it is assigned at the end
     ! of the daily history write process
     
     
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
      this%wood_product_harvest(:)        = 0._r8
      this%wood_product_landusechange(:)  = 0._r8
      this%burn_flux_to_atm  = 0._r8
      this%flux_generic_in   = 0._r8
      this%flux_generic_out  = 0._r8
      this%patch_resize_err  = 0._r8
      this%herbivory_flux_out= 0._r8

      return
  end subroutine ZeroMassBalFlux
   
  ! =====================================================================================

  subroutine CalculateTreeGrassAreaSite(csite, tree_fraction, grass_fraction, bare_fraction)
    !
    !  DESCRIPTION:
    !  Calculates total grass, tree, and bare fractions for a site

    ! ARGUMENTS:
    type(ed_site_type), intent(inout) :: csite          ! site object
    real(r8),           intent(out)   :: tree_fraction  ! total site tree fraction
    real(r8),           intent(out)   :: grass_fraction ! total site grass fraction
    real(r8),           intent(out)   :: bare_fraction  ! total site bare fraction

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch ! patch object
    
    tree_fraction = 0.0_r8
    grass_fraction = 0.0_r8
    
    currentPatch => csite%oldest_patch
    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        call currentPatch%UpdateTreeGrassArea()
        tree_fraction = tree_fraction + currentPatch%total_tree_area/AREA
        grass_fraction = grass_fraction + currentPatch%total_grass_area/AREA
      end if 
      currentPatch => currentPatch%younger
    end do

    ! if cover > 1.0, grasses are under the trees
    grass_fraction = min(grass_fraction, 1.0_r8 - tree_fraction)
    bare_fraction = 1.0_r8 - tree_fraction - grass_fraction

  end subroutine CalculateTreeGrassAreaSite

  !---------------------------------------------------------------------------------------

  subroutine dump_site(csite) 

   type(ed_site_type),intent(in),target :: csite


   ! EDTypes is 

   write(fates_log(),*) '----------------------------------------'
   write(fates_log(),*) ' Site Coordinates                       '
   write(fates_log(),*) '----------------------------------------'
   write(fates_log(),*) 'latitude                    = ', csite%lat
   write(fates_log(),*) 'longitude                   = ', csite%lon
   write(fates_log(),*) '----------------------------------------'
   return

  end subroutine dump_site

  ! =====================================================================================

  function get_current_landuse_statevector(this) result(current_state_vector)

     !
     ! !DESCRIPTION:
     !  Calculate how much of a site is each land use category.
     !  this does not include bare ground when nocomp + fixed biogeography is on,
     !  so will not sum to one in that case. otherwise it will sum to one.
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     class(ed_site_type) :: this
     real(r8)            :: current_state_vector(n_landuse_cats)

     ! !LOCAL VARIABLES:
     type(fates_patch_type), pointer :: currentPatch

     current_state_vector(:) = 0._r8

     currentPatch => this%oldest_patch
     do while (associated(currentPatch))
        if (currentPatch%land_use_label .gt. nocomp_bareground_land) then
           current_state_vector(currentPatch%land_use_label) = &
                current_state_vector(currentPatch%land_use_label) + &
                currentPatch%area/AREA
        end if
        currentPatch => currentPatch%younger
     end do

   end function get_current_landuse_statevector

   ! =====================================================================================

   function get_secondary_young_fraction(this) result(secondary_young_fraction)

     !
     ! !DESCRIPTION:
     !  Calculate how much of the secondary area is "young", i.e. below the age threshold.
     !  If no seconday patch area at all, return -1.
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     class(ed_site_type) :: this
     real(r8)            :: secondary_young_fraction
     real(r8)            :: secondary_young_area
     real(r8)            :: secondary_old_area

     ! !LOCAL VARIABLES:
     type(fates_patch_type), pointer :: currentPatch

     secondary_young_area = 0._r8
     secondary_old_area = 0._r8

     currentPatch => this%oldest_patch
     do while (associated(currentPatch))
        if (currentPatch%land_use_label .eq. secondaryland) then
           if ( currentPatch%age .ge. secondary_age_threshold ) then
              secondary_old_area = secondary_old_area + currentPatch%area
           else
              secondary_young_area = secondary_young_area + currentPatch%area
           end if
        end if
        currentPatch => currentPatch%younger
     end do

     if ( (secondary_young_area + secondary_old_area) .gt. nearzero ) then
        secondary_young_fraction = secondary_young_area / (secondary_young_area + secondary_old_area)
     else
        secondary_young_fraction = -1._r8
     endif

   end function get_secondary_young_fraction

end module EDTypesMod
