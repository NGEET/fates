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
  use EDParamsMod,           only : maxSWb, nclmax, nlevleaf
  use shr_log_mod,           only : errMsg => shr_log_errMsg

  
  implicit none
  private               ! By default everything is private
  save

  integer, parameter, public :: maxpft = 16               ! maximum number of PFTs allowed
                                                          ! the parameter file may determine that fewer
                                                          ! are used, but this helps allocate scratch
                                                          ! space and output arrays.
                                                  
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
  integer , parameter, public :: N_DIST_TYPES         = 3          ! Disturbance Modes 1) tree-fall, 2) fire, 3) logging
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
  integer , parameter, public :: N_DBH_BINS           = 6                 ! no. of dbh bins used when comparing patches
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
  !** Patch type structure           **
  !************************************

  type, public :: ed_patch_type

     ! POINTERS
     type (fates_cohort_type), pointer :: tallest => null()           ! pointer to patch's tallest cohort    
     type (fates_cohort_type), pointer :: shortest => null()          ! pointer to patch's shortest cohort
     type (ed_patch_type),  pointer :: older => null()             ! pointer to next older patch   
     type (ed_patch_type),  pointer :: younger => null()           ! pointer to next younger patch      

     !INDICES
     integer  :: patchno                                           ! unique number given to each new patch created for tracking

     ! PATCH INFO
     real(r8) ::  age                                              ! average patch age: years                   
     integer  ::  age_class                                        ! age class of the patch for history binning purposes
     real(r8) ::  area                                             ! patch area: m2  
     integer  ::  countcohorts                                     ! Number of cohorts in patch
     integer  ::  ncl_p                                            ! Number of occupied canopy layers
     integer  ::  anthro_disturbance_label                         ! patch label for anthropogenic disturbance classification
     real(r8) ::  age_since_anthro_disturbance                     ! average age for secondary forest since last anthropogenic disturbance


     ! Running means
     !class(rmean_type), pointer :: t2m                          ! Place-holder for 2m air temperature (variable window-size)
     class(rmean_type), pointer :: tveg24                        ! 24-hour mean vegetation temperature (K)
     class(rmean_type), pointer :: tveg_lpa                      ! Running mean of vegetation temperature at the
                                                                 ! leaf photosynthesis acclimation timescale [K]
     class(rmean_type), pointer :: tveg_longterm                ! Long-Term Running mean of vegetation temperature at the
                                                                 ! leaf photosynthesis acclimation timescale [K] (i.e T_home)

     integer  ::  nocomp_pft_label                               ! Where nocomp is active, use this label for patch ID.
                                                                 ! Each patch ID corresponds to a pft number since each
                                                                 ! patch has only one pft.  Bareground patches are given
                                                                 ! a zero integer as a label.
                                                                 ! If nocomp is not active this is set to unset.
                                                                 ! This is set in create_patch as an argument
                                                                 ! to that procedure.


     ! LEAF ORGANIZATION
     real(r8) ::  pft_agb_profile(maxpft,n_dbh_bins)            ! binned above ground biomass, for patch fusion: KgC/m2
     real(r8) ::  canopy_layer_tlai(nclmax)                     ! total leaf area index of each canopy layer
                                                                ! used to determine attenuation of parameters during
                                                                ! photosynthesis m2 veg / m2 of canopy area (patch without bare ground)
     real(r8) ::  total_canopy_area                             ! area that is covered by vegetation : m2
     real(r8) ::  total_tree_area                               ! area that is covered by woody vegetation : m2
     real(r8) ::  zstar                                         ! height of smallest canopy tree -- only meaningful in "strict PPA" mode

     real(r8) :: c_stomata                                      ! Mean stomatal conductance of all leaves in the patch   [umol/m2/s]
     real(r8) :: c_lblayer                                      ! Mean boundary layer conductance of all leaves in the patch [umol/m2/s]
      
                                                                ! UNITS for the ai profiles
                                                                ! [ m2 leaf / m2 contributing crown footprints]
     real(r8) ::  tlai_profile(nclmax,maxpft,nlevleaf)          ! total   leaf area in each canopy layer, pft, and leaf layer. 
     real(r8) ::  elai_profile(nclmax,maxpft,nlevleaf)          ! exposed leaf area in each canopy layer, pft, and leaf layer
     real(r8) ::  tsai_profile(nclmax,maxpft,nlevleaf)          ! total   stem area in each canopy layer, pft, and leaf layer
     real(r8) ::  esai_profile(nclmax,maxpft,nlevleaf)          ! exposed stem area in each canopy layer, pft, and leaf layer
     real(r8) ::  radiation_error                               ! radiation error (w/m2)
     real(r8) ::  layer_height_profile(nclmax,maxpft,nlevleaf)
     real(r8) ::  canopy_area_profile(nclmax,maxpft,nlevleaf)   ! fraction of crown area per canopy area in each layer
                                                                ! they will sum to 1.0 in the fully closed canopy layers
                                                                ! but only in leaf-layers that contain contributions
                                                                ! from all cohorts that donate to canopy_area


     ! layer, pft, and leaf layer:-
     integer  ::  canopy_mask(nclmax,maxpft)                    ! is there any of this pft in this canopy layer?      
     integer  ::  nrad(nclmax,maxpft)                           ! number of exposed leaf layers for each canopy layer and pft
     integer  ::  ncan(nclmax,maxpft)                           ! number of total   leaf layers for each canopy layer and pft

     !RADIATION FLUXES      
     real(r8) :: fcansno                                        ! Fraction of canopy covered in snow

     logical  ::  solar_zenith_flag                             ! integer flag specifying daylight (based on zenith angle)
     real(r8) ::  solar_zenith_angle                            ! solar zenith angle (radians)

     real(r8) ::  gnd_alb_dif(maxSWb)                           ! ground albedo for diffuse rad, both bands (fraction)
     real(r8) ::  gnd_alb_dir(maxSWb)                           ! ground albedo for direct rad, both bands (fraction)
     
     real(r8) ::  fabd_sun_z(nclmax,maxpft,nlevleaf)            ! sun fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabd_sha_z(nclmax,maxpft,nlevleaf)            ! shade fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sun_z(nclmax,maxpft,nlevleaf)            ! sun fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sha_z(nclmax,maxpft,nlevleaf)            ! shade fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-

     real(r8) ::  ed_laisun_z(nclmax,maxpft,nlevleaf)           ! amount of LAI in the sun   in each canopy layer, 
     ! pft, and leaf layer. m2/m2
     real(r8) ::  ed_laisha_z(nclmax,maxpft,nlevleaf)           ! amount of LAI in the shade in each canopy layer,
     real(r8) ::  ed_parsun_z(nclmax,maxpft,nlevleaf)           ! PAR absorbed  in the sun   in each canopy layer,
     real(r8) ::  ed_parsha_z(nclmax,maxpft,nlevleaf)           ! PAR absorbed  in the shade in each canopy layer,
     real(r8) ::  f_sun(nclmax,maxpft,nlevleaf)                 ! fraction of leaves in the sun in each canopy layer, pft, 

     ! radiation profiles for comparison against observations

     ! normalized direct photosynthetically active radiation profiles by 
     ! incident type (direct/diffuse at top of canopy),leaf,pft,leaf (unitless)
     real(r8) ::  nrmlzd_parprof_pft_dir_z(n_rad_stream_types,nclmax,maxpft,nlevleaf)  

     ! normalized diffuse photosynthetically active radiation profiles by 
     ! incident type (direct/diffuse at top of canopy),leaf,pft,leaf (unitless)
     real(r8) ::  nrmlzd_parprof_pft_dif_z(n_rad_stream_types,nclmax,maxpft,nlevleaf)  

     ! normalized direct photosynthetically active radiation profiles by 
     ! incident type (direct/diffuse at top of canopy),leaf,leaf (unitless) 
     real(r8) ::  nrmlzd_parprof_dir_z(n_rad_stream_types,nclmax,nlevleaf)         

     ! normalized diffuse photosynthetically active radiation profiles by 
     ! incident type (direct/diffuse at top of canopy),leaf,leaf (unitless) 
     real(r8) ::  nrmlzd_parprof_dif_z(n_rad_stream_types,nclmax,nlevleaf)
         
     real(r8) ::  parprof_pft_dir_z(nclmax,maxpft,nlevleaf)   ! direct-beam PAR profile through canopy, by canopy,PFT,leaf level (w/m2)
     real(r8) ::  parprof_pft_dif_z(nclmax,maxpft,nlevleaf)   ! diffuse     PAR profile through canopy, by canopy,PFT,leaf level (w/m2)
     real(r8) ::  parprof_dir_z(nclmax,nlevleaf)              ! direct-beam PAR profile through canopy, by canopy,leaf level (w/m2)
     real(r8) ::  parprof_dif_z(nclmax,nlevleaf)              ! diffuse     PAR profile through canopy, by canopy,leaf level (w/m2)

     ! and leaf layer. m2/m2
     real(r8),allocatable ::  tr_soil_dir(:)                              ! fraction of incoming direct  radiation that (cm_numSWb)
     ! is transmitted to the soil as direct
     real(r8),allocatable ::  tr_soil_dif(:)                              ! fraction of incoming diffuse radiation that 
     ! is transmitted to the soil as diffuse
     real(r8),allocatable ::  tr_soil_dir_dif(:)                          ! fraction of incoming direct  radiation that 
     ! is transmitted to the soil as diffuse
     real(r8),allocatable ::  fab(:)                                      ! fraction of incoming total   radiation that is absorbed by the canopy
     real(r8),allocatable ::  fabd(:)                                     ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8),allocatable ::  fabi(:)                                     ! fraction of incoming diffuse radiation that is absorbed by the canopy
     real(r8),allocatable ::  sabs_dir(:)                                 ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8),allocatable ::  sabs_dif(:)                                 ! fraction of incoming diffuse radiation that is absorbed by the canopy


     ! PHOTOSYNTHESIS       

     real(r8) ::  psn_z(nclmax,maxpft,nlevleaf)               ! carbon assimilation in each canopy layer, pft, and leaf layer. umolC/m2/s

     ! ROOTS
     real(r8) ::  btran_ft(maxpft)                ! btran calculated seperately for each PFT:-
     real(r8) ::  bstress_sal_ft(maxpft)          ! bstress from salinity calculated seperately for each PFT:-   


     ! These two variables are only used for external seed rain currently.
     real(r8) :: nitr_repro_stoich(maxpft)        ! The NC ratio of a new recruit in this patch
     real(r8) :: phos_repro_stoich(maxpft)        ! The PC ratio of a new recruit in this patch

     
     ! DISTURBANCE 
     real(r8) ::  disturbance_rates(n_dist_types)                  ! disturbance rate from 1) mortality 
                                                                   !                       2) fire: fraction/day 
                                                                   !                       3) logging mortatliy
     real(r8) ::  fract_ldist_not_harvested                        ! fraction of logged area that is canopy trees that weren't harvested


     ! Litter and Coarse Woody Debris

     type(litter_type), pointer :: litter(:)  ! Litter (leaf,fnrt,CWD and seeds) for different elements

     real(r8),allocatable :: fragmentation_scaler(:)            ! Scale rate of litter fragmentation based on soil layer. 0 to 1.

     !FUEL CHARECTERISTICS
     real(r8) ::  sum_fuel                                         ! total ground fuel related to ros (omits 1000hr fuels): KgC/m2
     real(r8) ::  fuel_frac(nfsc)                                  ! fraction of each litter class in the ros_fuel:-.  
     real(r8) ::  livegrass                                        ! total aboveground grass biomass in patch.  KgC/m2
     real(r8) ::  fuel_bulkd                                       ! average fuel bulk density of the ground fuel. kgBiomass/m3
                                                                   ! (incl. live grasses. omits 1000hr fuels). KgC/m3
     real(r8) ::  fuel_sav                                         ! average surface area to volume ratio of the ground fuel. cm-1
                                                                   ! (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_mef                                         ! average moisture of extinction factor 
                                                                   ! of the ground fuel (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_eff_moist                                   ! effective avearage fuel moisture content of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels)
     real(r8) ::  litter_moisture(nfsc)

     ! FIRE SPREAD
     real(r8) ::  ros_front                                        ! rate of forward  spread of fire: m/min
     real(r8) ::  ros_back                                         ! rate of backward spread of fire: m/min
     real(r8) ::  effect_wspeed                                    ! windspeed modified by fraction of relative grass and tree cover: m/min
     real(r8) ::  tau_l                                            ! Duration of lethal heating: mins
     real(r8) ::  fi                                               ! average fire intensity of flaming front:  kj/m/s or kw/m
     integer  ::  fire                                             ! Is there a fire? 1=yes 0=no
     real(r8) ::  fd                                               ! fire duration: mins

     ! FIRE EFFECTS     
     real(r8) ::  scorch_ht(maxpft)                                ! scorch height: m 
     real(r8) ::  frac_burnt                                       ! fraction burnt: frac patch/day  
     real(r8) ::  tfc_ros                                          ! total intensity-relevant fuel consumed - no trunks.  KgC/m2 of burned ground/day
     real(r8) ::  burnt_frac_litter(nfsc)                          ! fraction of each litter pool burned, conditional on it being burned


     ! PLANT HYDRAULICS   (not currently used in hydraulics RGK 03-2018)  
     ! type(ed_patch_hydr_type) , pointer :: pa_hydr              ! All patch hydraulics data, see FatesHydraulicsMemMod.F90


  end type ed_patch_type

  
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
  

  !************************************
  !** Site type structure           **
  !************************************

  type, public :: ed_site_type
     
     ! POINTERS  
     type (ed_patch_type), pointer :: oldest_patch => null()   ! pointer to oldest patch at the site  
     type (ed_patch_type), pointer :: youngest_patch => null() ! pointer to yngest patch at the site
     
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
     real(r8), allocatable :: area_PFT(:)                      ! Area allocated to individual PFTs    
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

     
     
     ! SP mode target PFT level variables
     real(r8), allocatable :: sp_tlai(:)                      ! target TLAI per FATES pft
     real(r8), allocatable :: sp_tsai(:)                      ! target TSAI per FATES pft
     real(r8), allocatable :: sp_htop(:)                      ! target HTOP per FATES pft
     
     ! Mass Balance (allocation for each element)

     type(site_massbal_type), pointer :: mass_balance(:)

     ! Flux diagnostics (allocation for each element)

     type(site_fluxdiags_type), pointer :: flux_diags(:)

     ! PHENOLOGY 
     real(r8) ::  grow_deg_days                                ! Phenology growing degree days
     real(r8) ::  snow_depth                                   ! site-level snow depth (used for ELAI/TLAI calcs)

     integer  ::  cstatus                                      ! are leaves in this pixel on or off for cold decid
                                                               ! 0 = this site has not experienced a cold period over at least
                                                               !     400 days, leaves are dropped and flagged as non-cold region
                                                               ! 1 = this site is in a cold-state where leaves should have fallen
                                                               ! 2 = this site is in a warm-state where leaves are allowed to flush
     integer  ::  dstatus                                      ! are leaves in this pixel on or off for drought decid
                                                               ! 0 = leaves off due to time exceedance
                                                               ! 1 = leaves off due to moisture avail
                                                               ! 2 = leaves on due to moisture avail
                                                               ! 3 = leaves on due to time exceedance
     integer  ::  nchilldays                                   ! num chilling days: (for botta gdd trheshold calculation)
     integer  ::  ncolddays                                    ! num cold days: (must exceed threshold to drop leaves)
     real(r8) ::  vegtemp_memory(num_vegtemp_mem)              ! record of last 10 days temperature for senescence model. deg C
     integer  ::  cleafondate                                  ! model date (day integer) of leaf on (cold):-
     integer  ::  cleafoffdate                                 ! model date (day integer) of leaf off (cold):-
     integer  ::  dleafondate                                  ! model date (day integer) of leaf on drought:-
     integer  ::  dleafoffdate                                 ! model date (day integer) of leaf off drought:-
     integer  ::  phen_model_date                              ! current model date (day integer)
                                                               ! this date stays continuous when
                                                               ! in runs that are restarted, regardless of
                                                               ! the conditions of restart

     real(r8) ::  water_memory(numWaterMem)                             ! last 10 days of soil moisture memory...


     ! FIRE
     real(r8) ::  wind                                         ! daily wind in m/min for Spitfire units 
     real(r8) ::  acc_ni                                       ! daily nesterov index accumulating over time.
     real(r8) ::  fdi                                          ! daily probability an ignition event will start a fire
     real(r8) ::  NF                                           ! daily ignitions in km2
     real(r8) ::  NF_successful                                ! daily ignitions in km2 that actually lead to fire

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

     real(r8), allocatable :: term_nindivs_canopy(:,:) ! number of canopy individuals that were in cohorts which 
                                                       ! were terminated this timestep, on size x pft
     real(r8), allocatable :: term_nindivs_ustory(:,:) ! number of understory individuals that were in cohorts which 
                                                       ! were terminated this timestep, on size x pft

     real(r8), allocatable :: term_carbonflux_canopy(:)  ! carbon flux from live to dead pools associated 
                                                         ! with termination mortality, per canopy level
     real(r8), allocatable :: term_carbonflux_ustory(:)  ! carbon flux from live to dead pools associated 
                                                         ! with termination mortality, per canopy level    
     real(r8), allocatable :: imort_carbonflux(:)        ! biomass of individuals killed due to impact mortality per year. [kgC/ha/day]
     real(r8), allocatable :: fmort_carbonflux_canopy(:) ! biomass of canopy indivs killed due to fire per year. [gC/m2/sec]
     real(r8), allocatable :: fmort_carbonflux_ustory(:) ! biomass of understory indivs killed due to fire per year [gC/m2/sec] 

     real(r8), allocatable :: term_abg_flux(:,:)          ! aboveground biomass lost due to termination mortality x size x pft
     real(r8), allocatable :: imort_abg_flux(:,:)         ! aboveground biomass lost due to impact mortality x size x pft
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
     real(r8), allocatable :: imort_cflux_damage(:,:)         ! carbon flux from impact mortality by damage class
     real(r8), allocatable :: term_cflux_canopy_damage(:,:)          ! carbon flux from termination mortality by damage class
     real(r8), allocatable :: term_cflux_ustory_damage(:,:)          ! carbon flux from termination mortality by damage class

     real(r8), allocatable :: growthflux_fusion(:,:)     ! rate of individuals moving into a given size class bin
                                                         ! due to fusion in a given day. on size x pft array 


     real(r8)              :: crownarea_canopy_damage                 ! crown area of canopy that is damaged annually  
     real(r8)              :: crownarea_ustory_damage                 ! crown area of understory that is damaged annually
     
     ! Canopy Spread
     real(r8) ::  spread                                          ! dynamic canopy allometric term [unitless]

     ! site-level variables to keep track of the disturbance rates, both actual and "potential"
     real(r8) :: disturbance_rates_primary_to_primary(N_DIST_TYPES)      ! actual disturbance rates from primary patches to primary patches [m2/m2/day]
     real(r8) :: disturbance_rates_primary_to_secondary(N_DIST_TYPES)    ! actual disturbance rates from primary patches to secondary patches [m2/m2/day]
     real(r8) :: disturbance_rates_secondary_to_secondary(N_DIST_TYPES)  ! actual disturbance rates from secondary patches to secondary patches [m2/m2/day]
     real(r8) :: potential_disturbance_rates(N_DIST_TYPES)               ! "potential" disturb rates (i.e. prior to the "which is most" logic) [m2/m2/day]
     real(r8) :: primary_land_patchfusion_error                          ! error term in total area of primary patches associated with patch fusion [m2/m2/day]
     
  end type ed_site_type

  ! Make public necessary subroutines and functions
  public :: val_check_ed_vars
  public :: dump_site
  public :: dump_patch
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
  
  subroutine val_check_ed_vars(currentPatch,var_aliases,return_code)

     ! ----------------------------------------------------------------------------------
     ! Perform numerical checks on variables of interest.
     ! The input string is of the form:  'VAR1_NAME:VAR2_NAME:VAR3_NAME'
     ! ----------------------------------------------------------------------------------


     use FatesUtilsMod,only : check_hlm_list
     use FatesUtilsMod,only : check_var_real

     ! Arguments
     type(ed_patch_type),intent(in), target :: currentPatch
     character(len=*),intent(in)            :: var_aliases
     integer,intent(out)                    :: return_code ! return 0 for all fine
                                                           ! return 1 if a nan detected
                                                           ! return 10+ if an overflow
                                                           ! return 100% if an underflow
     ! Locals
     type(fates_cohort_type), pointer          :: currentCohort

     
     ! Check through a registry of variables to check
     
     if ( check_hlm_list(trim(var_aliases),'co_n') ) then

        currentCohort => currentPatch%shortest
        do while(associated(currentCohort))
           call check_var_real(currentCohort%n,'cohort%n',return_code)
           if(.not.(return_code.eq.0)) then
              call dump_patch(currentPatch)
              call currentCohort%dump()
              return
           end if
           currentCohort => currentCohort%taller
        end do
     end if
     
     if ( check_hlm_list(trim(var_aliases),'co_dbh') ) then

        currentCohort => currentPatch%shortest
        do while(associated(currentCohort))        
           call check_var_real(currentCohort%dbh,'cohort%dbh',return_code)
           if(.not.(return_code.eq.0)) then
              call dump_patch(currentPatch)
              call currentCohort%dump()
              return
           end if
           currentCohort => currentCohort%taller
        end do
     end if

     if ( check_hlm_list(trim(var_aliases),'pa_area') ) then

        call check_var_real(currentPatch%area,'patch%area',return_code)
        if(.not.(return_code.eq.0)) then
           call dump_patch(currentPatch)
           return
        end if
     end if
     


     return
  end subroutine val_check_ed_vars

  ! =====================================================================================

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


  subroutine dump_patch(cpatch)

     type(ed_patch_type),intent(in),target :: cpatch

     ! locals
     integer :: el  ! element loop counting index

     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) ' Dumping Patch Information              '
     write(fates_log(),*) ' (omitting arrays)                      '
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) 'pa%patchno            = ',cpatch%patchno
     write(fates_log(),*) 'pa%age                = ',cpatch%age
     write(fates_log(),*) 'pa%age_class          = ',cpatch%age_class
     write(fates_log(),*) 'pa%area               = ',cpatch%area
     write(fates_log(),*) 'pa%countcohorts       = ',cpatch%countcohorts
     write(fates_log(),*) 'pa%ncl_p              = ',cpatch%ncl_p
     write(fates_log(),*) 'pa%total_canopy_area  = ',cpatch%total_canopy_area
     write(fates_log(),*) 'pa%total_tree_area    = ',cpatch%total_tree_area
     write(fates_log(),*) 'pa%zstar              = ',cpatch%zstar
     write(fates_log(),*) 'pa%solar_zenith_flag  = ',cpatch%solar_zenith_flag
     write(fates_log(),*) 'pa%solar_zenith_angle = ',cpatch%solar_zenith_angle
     write(fates_log(),*) 'pa%gnd_alb_dif        = ',cpatch%gnd_alb_dif(:)
     write(fates_log(),*) 'pa%gnd_alb_dir        = ',cpatch%gnd_alb_dir(:)
     write(fates_log(),*) 'pa%c_stomata          = ',cpatch%c_stomata
     write(fates_log(),*) 'pa%c_lblayer          = ',cpatch%c_lblayer
     write(fates_log(),*) 'pa%disturbance_rates  = ',cpatch%disturbance_rates(:)
     write(fates_log(),*) 'pa%anthro_disturbance_label = ',cpatch%anthro_disturbance_label
     write(fates_log(),*) '----------------------------------------'
     do el = 1,num_elements
        write(fates_log(),*) 'element id: ',element_list(el)
        write(fates_log(),*) 'seed mass: ',sum(cpatch%litter(el)%seed)
        write(fates_log(),*) 'seed germ mass: ',sum(cpatch%litter(el)%seed_germ)
        write(fates_log(),*) 'leaf fines(pft): ',sum(cpatch%litter(el)%leaf_fines)
        write(fates_log(),*) 'root fines(pft,sl): ',sum(cpatch%litter(el)%root_fines)
        write(fates_log(),*) 'ag_cwd(c): ',sum(cpatch%litter(el)%ag_cwd)
        write(fates_log(),*) 'bg_cwd(c,sl): ',sum(cpatch%litter(el)%bg_cwd)
     end do

     return

  end subroutine dump_patch

  ! =====================================================================================
  
  

end module EDTypesMod
