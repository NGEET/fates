module EDTypesMod

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesConstantsMod,     only : ifalse
  use FatesConstantsMod,     only : itrue
  use FatesGlobals,          only : fates_log
  use FatesHydraulicsMemMod, only : ed_cohort_hydr_type
  use FatesHydraulicsMemMod, only : ed_site_hydr_type
  use PRTGenericMod,         only : prt_vartypes
  use PRTGenericMod,         only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,         only : repro_organ, store_organ, struct_organ
  use PRTGenericMod,         only : all_carbon_elements
  use PRTGenericMod,         only : num_element_types
  use FatesLitterMod,        only : litter_type
  use FatesLitterMod,        only : ncwd
  use FatesConstantsMod,     only : n_anthro_disturbance_categories
  use FatesConstantsMod,     only : days_per_year
  
  implicit none
  private               ! By default everything is private
  save

  integer, parameter, public :: maxPatchesPerSite  = 14   ! maximum number of patches to live on a site
  integer, parameter, public :: maxPatchesPerSite_by_disttype(n_anthro_disturbance_categories)  = &
                                                     (/ 10, 4 /)  !!! MUST SUM TO maxPatchesPerSite !!!
  integer,  public :: maxCohortsPerPatch = 100            ! maximum number of cohorts per patch
  
  integer, parameter, public :: nclmax = 2                ! Maximum number of canopy layers
  integer, parameter, public :: ican_upper = 1            ! Nominal index for the upper canopy
  integer, parameter, public :: ican_ustory = 2           ! Nominal index for diagnostics that refer
                                                          ! to understory layers (all layers that
                                                          ! are not the top canopy layer)

  integer, parameter, public :: nlevleaf = 30             ! number of leaf layers in canopy layer
  integer, parameter, public :: maxpft = 15               ! maximum number of PFTs allowed
                                                          ! the parameter file may determine that fewer
                                                          ! are used, but this helps allocate scratch
                                                          ! space and output arrays.
                                                  
  integer, parameter, public :: max_nleafage = 4          ! This is the maximum number of leaf age pools, 
                                                          ! used for allocating scratch space


  ! -------------------------------------------------------------------------------------
  ! Radiation parameters
  ! These should be part of the radiation module, but since we only have one option
  ! this is ok for now. (RGK 04-2018)
  ! -------------------------------------------------------------------------------------


  integer, parameter, public :: n_rad_stream_types = 2    ! The number of radiation streams used (direct/diffuse)
 
  integer, parameter, public :: idirect   = 1             ! This is the array index for direct radiation
  integer, parameter, public :: idiffuse  = 2             ! This is the array index for diffuse radiation


  ! TODO: we use this cp_maxSWb only because we have a static array q(size=2) of
  ! land-ice abledo for vis and nir.  This should be a parameter, which would
  ! get us on track to start using multi-spectral or hyper-spectral (RGK 02-2017)

  integer, parameter, public :: maxSWb = 2      ! maximum number of broad-bands in the
                                                ! shortwave spectrum cp_numSWb <= cp_maxSWb
                                                ! this is just for scratch-array purposes
                                                ! if cp_numSWb is larger than this value
                                                ! simply bump this number up as needed

  integer, parameter, public :: ivis = 1        ! This is the array index for short-wave
                                                ! radiation in the visible spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod
  integer, parameter, public :: inir = 2        ! This is the array index for short-wave
                                                ! radiation in the near-infrared spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod

  integer, parameter, public :: ipar = ivis     ! The photosynthetically active band
                                                ! can be approximated to be equal to the visible band


  integer, parameter, public :: leaves_on  = 2  ! Flag specifying that a deciduous plant has leaves
                                                ! and should be allocating to them as well
  integer, parameter, public :: leaves_off = 1  ! Flag specifying that a deciduous plant has dropped
                                                ! its leaves and should not be trying to allocate
                                                ! towards any growth.

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
  real(r8), parameter, public :: dinc_ed              = 1.0_r8     ! size of VAI bins (LAI+SAI)  [CHANGE THIS NAME WITH NEXT INTERFACE
                                                           ! UPDATE]
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


  ! SPITFIRE     

  integer,  parameter, public :: NFSC                 = NCWD+2     ! number fuel size classes  (4 cwd size classes, leaf litter, and grass)
  integer,  parameter, public :: tw_sf                = 1          ! array index of twig pool for spitfire
  integer,  parameter, public :: lb_sf                = 3          ! array index of large branch pool for spitfire
  integer,  parameter, public :: tr_sf                = 4          ! array index of dead trunk pool for spitfire
  integer,  parameter, public :: dl_sf                = 5          ! array index of dead leaf pool for spitfire (dead grass and dead leaves)
  integer,  parameter, public :: lg_sf                = 6          ! array index of live grass pool for spitfire

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

  character*4 yearchar                    

  ! special mode to cause PFTs to create seed mass of all currently-existing PFTs
  logical, parameter, public :: homogenize_seed_pfts  = .false.

  

  ! Global identifiers for which elements we are using (apply mostly to litter)

  integer, public              :: num_elements          ! This is the number of elements in this simulation
                                                        ! e.g. (C,N,P,K, etc)
  integer, allocatable, public :: element_list(:)       ! This vector holds the element ids that are found
                                                        ! in PRTGenericMod.F90. examples are carbon12_element
                                                        ! nitrogen_element, etc.

  integer, public :: element_pos(num_element_types)       ! This is the reverse lookup
                                                        ! for element types. Pick an element
                                                        ! global index, and it gives you
                                                        ! the position in the element_list

  !************************************
  !** COHORT type structure          **
  !************************************
  type, public :: ed_cohort_type

     ! POINTERS
     type (ed_cohort_type) , pointer :: taller   => null()       ! pointer to next tallest cohort     
     type (ed_cohort_type) , pointer :: shorter  => null()       ! pointer to next shorter cohort     
     type (ed_patch_type)  , pointer :: patchptr => null()       ! pointer to patch that cohort is in


     
     ! Multi-species, multi-organ Plant Reactive Transport (PRT)
     ! Contains carbon and nutrient state variables for various plant organs

     class(prt_vartypes), pointer :: prt

     ! VEGETATION STRUCTURE
     integer  ::  pft                                    ! pft number
     real(r8) ::  n                                      ! number of individuals in cohort per 'area' (10000m2 default)
     real(r8) ::  dbh                                    ! dbh: cm
     real(r8) ::  coage                                  ! cohort age in years
     real(r8) ::  hite                                   ! height: meters
     integer  ::  indexnumber                            ! unique number for each cohort. (within clump?)
     real(r8) ::  laimemory                              ! target leaf biomass- set from previous year: kGC per indiv
     real(r8) ::  sapwmemory                             ! target sapwood biomass- set from previous year: kGC per indiv
     real(r8) ::  structmemory                           ! target structural biomass- set from previous year: kGC per indiv
     integer  ::  canopy_layer                           ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)
     real(r8) ::  canopy_layer_yesterday                 ! recent canopy status of cohort
                                                         ! (1 = canopy, 2 = understorey, etc.)  
                                                         ! real to be conservative during fusion

     real(r8) ::  lai                                    ! leaf area index of cohort: m2 leaf area of entire cohort per m2 of canopy area of a patch
     real(r8) ::  sai                                    ! stem area index of cohort: m2 leaf area of entire cohort per m2 of canopy area of a patch
     real(r8) ::  g_sb_laweight                          ! Total conductance (stomata+boundary layer) of the cohort, weighted by its leaf area [m/s]*[m2]
     real(r8) ::  canopy_trim                            ! What is the fraction of the maximum leaf biomass that we are targeting? :-
     real(r8) ::  leaf_cost                              ! How much does it cost to maintain leaves: kgC/m2/year-1
     real(r8) ::  excl_weight                            ! How much of this cohort is demoted each year, as a proportion of all cohorts:-
     real(r8) ::  prom_weight                            ! How much of this cohort is promoted each year, as a proportion of all cohorts:-
     integer  ::  nv                                     ! Number of leaf layers: -
     integer  ::  status_coh                             ! growth status of plant  (2 = leaves on , 1 = leaves off)
     real(r8) ::  c_area                                 ! areal extent of canopy (m2)
     real(r8) ::  treelai                                ! lai of an individual within cohort leaf area (m2) / crown area (m2)
     real(r8) ::  treesai                                ! stem area index of an indiv. within cohort: stem area (m2) / crown area (m2)
     logical  ::  isnew                                  ! flag to signify a new cohort, new cohorts have not experienced
                                                         ! npp or mortality and should therefore not be fused or averaged
     integer  ::  size_class                             ! An index that indicates which diameter size bin the cohort currently resides in
                                                         ! this is used for history output. We maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls
     integer  ::  coage_class                            ! An index that indicates which age bin the cohort currently resides in 
                                                         ! used for history output.
     integer  ::  size_by_pft_class                      ! An index that indicates the cohorts position of the joint size-class x functional
                                                         ! type classification. We also maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls
     integer  ::  coage_by_pft_class                     ! An index that indicates the cohorts position of the join cohort age class x PFT 
     integer ::  size_class_lasttimestep                 ! size class of the cohort at the last time step

     ! CARBON FLUXES 
     
     ! ----------------------------------------------------------------------------------
     ! NPP, GPP and RESP: Instantaneous, accumulated and accumulated-hold types.*
     ! 
     ! _tstep:    The instantaneous estimate that is calculated at each rapid plant biophysics
     !            time-step (ie photosynthesis, sub-hourly). (kgC/indiv/timestep)
     ! _acc:      The accumulation of the _tstep variable from the beginning to ending of
     !            the dynamics time-scale.  This variable is zero'd during initialization and
     !            after the dynamics call-sequence is completed.  (kgC/indiv/day)
     ! _acc_hold: While _acc is zero'd after the dynamics call sequence and then integrated, 
     !            _acc_hold "holds" the integrated value until the next time dynamics is 
     !            called. This is necessary for restarts. This variable also has units
     !            converted to a useful rate (kgC/indiv/yr)
     ! ----------------------------------------------------------------------------------

     real(r8) ::  gpp_tstep          ! Gross Primary Production (see above *)
     real(r8) ::  gpp_acc
     real(r8) ::  gpp_acc_hold

     real(r8) ::  npp_tstep          ! Net Primary Production (see above *)
     real(r8) ::  npp_acc
     real(r8) ::  npp_acc_hold

     real(r8) ::  resp_tstep         ! Autotrophic respiration (see above *)
     real(r8) ::  resp_acc
     real(r8) ::  resp_acc_hold
     
     ! carbon 13c discrimination
     real(r8) ::  c13disc_clm         ! carbon 13 discrimination in new synthesized carbon: part-per-mil, at each indiv/timestep
     real(r8) ::  c13disc_acc         ! carbon 13 discrimination in new synthesized carbon: part-per-mil, at each indiv/day, at the end of a day


     ! The following four biophysical rates are assumed to be
     ! at the canopy top, at reference temp 25C, and based on the 
     ! leaf age weighted average of the PFT parameterized values. The last
     ! condition is why it is dynamic and tied to the cohort

     real(r8) :: vcmax25top  ! Maximum carboxylation at the cohort's top 
                             ! at reference temperature (25C).
     real(r8) :: jmax25top   ! canopy top: maximum electron transport 
                             ! rate at 25C (umol electrons/m**2/s)
     real(r8) :: tpu25top    ! canopy top: triose phosphate utilization
                             ! rate at 25C (umol CO2/m**2/s)
     real(r8) :: kp25top     ! canopy top: initial slope of CO2 response
                             ! curve (C4 plants) at 25C



     real(r8) ::  ts_net_uptake(nlevleaf)              ! Net uptake of leaf layers: kgC/m2/timestep
     real(r8) ::  year_net_uptake(nlevleaf)            ! Net uptake of leaf layers: kgC/m2/year

     ! RESPIRATION COMPONENTS
     real(r8) ::  rdark                                  ! Dark respiration: kgC/indiv/s
     real(r8) ::  resp_g                                 ! Growth respiration:  kgC/indiv/timestep
     real(r8) ::  resp_m                                 ! Maintenance respiration:  kgC/indiv/timestep 
     real(r8) ::  livestem_mr                            ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (Above ground)
     real(r8) ::  livecroot_mr                           ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (below ground)
     real(r8) ::  froot_mr                               ! Live fine root   maintenance respiration: kgC/indiv/s

     !MORTALITY
     real(r8) ::  dmort                                  ! proportional mortality rate. (year-1)

     ! Mortality Rate Partitions
     real(r8) ::  bmort                                  ! background mortality rate        n/year
     real(r8) ::  cmort                                  ! carbon starvation mortality rate n/year
     real(r8) ::  hmort                                  ! hydraulic failure mortality rate n/year
     real(r8) ::  frmort                                 ! freezing mortality               n/year
     real(r8) ::  smort                                  ! senesence mortality              n/year
     real(r8) ::  asmort                                 ! age senescence mortality         n/year
     
      ! Logging Mortality Rate 
      ! Yi Xu & M. Huang
     real(r8) ::  lmort_direct                           ! directly logging rate            fraction /per logging activity
     real(r8) ::  lmort_collateral                       ! collaterally damaged rate        fraction /per logging activity
     real(r8) ::  lmort_infra                            ! mechanically damaged rate        fraction /per logging activity
     real(r8) ::  l_degrad                               ! rate of trees that are not killed but suffer from forest degradation
                                                         ! (i.e. they are moved to newly-anthro-disturbed secondary 
                                                         !  forest patch).  fraction /per logging activity

     real(r8) :: seed_prod                               ! diagnostic seed production rate [kgC/plant/day]

     ! NITROGEN POOLS      
     ! ----------------------------------------------------------------------------------
     ! Nitrogen pools are not prognostic in the current implementation.
     ! They are diagnosed during photosynthesis using a simple C2N parameter. Local values
     ! used in that routine.
     ! ----------------------------------------------------------------------------------

     ! GROWTH DERIVIATIVES
     real(r8) ::  dndt                                   ! time derivative of cohort size  : n/year
     real(r8) ::  dhdt                                   ! time derivative of height       : m/year
     real(r8) ::  ddbhdt                                 ! time derivative of dbh          : cm/year
     real(r8) ::  dbdeaddt                               ! time derivative of dead biomass         : KgC/year

     ! FIRE
     real(r8) ::  fraction_crown_burned                  ! proportion of crown affected by fire:-
     real(r8) ::  cambial_mort                           ! probability that trees dies due to cambial char 
                                                         ! (conditional on the tree being subjected to the fire)
     real(r8) ::  crownfire_mort                         ! probability of tree post-fire mortality 
                                                         ! due to crown scorch (conditional on the tree being subjected to the fire)
     real(r8) ::  fire_mort                              ! post-fire mortality from cambial and crown damage assuming two are independent:-

     ! Hydraulics
     type(ed_cohort_hydr_type), pointer :: co_hydr       ! All cohort hydraulics data, see FatesHydraulicsMemMod.F90

  end type ed_cohort_type

  !************************************
  !** Patch type structure           **
  !************************************

  type, public :: ed_patch_type

     ! POINTERS
     type (ed_cohort_type), pointer :: tallest => null()           ! pointer to patch's tallest cohort    
     type (ed_cohort_type), pointer :: shortest => null()          ! pointer to patch's shortest cohort
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
     real(r8) ::  btran_ft(maxpft)                              ! btran calculated seperately for each PFT:-
     real(r8) ::  bstress_sal_ft(maxpft)                        ! bstress from salinity calculated seperately for each PFT:-   
     

     ! DISTURBANCE 
     real(r8) ::  disturbance_rates(n_dist_types)                  ! disturbance rate from 1) mortality 
                                                                   !                       2) fire: fraction/day 
                                                                   !                       3) logging mortatliy
     real(r8) ::  disturbance_rate                                 ! larger effective disturbance rate: fraction/day
     integer  ::  disturbance_mode                                 ! index identifying which disturbance was applied
                                                                   ! can be one of: dtype_ifall, dtype_ilog or dtype_ifire
     real(r8) ::  fract_ldist_not_harvested                        ! fraction of logged area that is canopy trees that weren't harvested


     ! Litter and Coarse Woody Debris

     type(litter_type), pointer :: litter(:)  ! Litter (leaf,fnrt,CWD and seeds) for different elements

     real(r8) :: fragmentation_scaler          ! Scale rate of litter fragmentation. 0 to 1.

     real(r8) ::  repro(maxpft)                                 ! allocation to reproduction per PFT : KgC/m2

     !FUEL CHARECTERISTICS
     real(r8) ::  sum_fuel                                         ! total ground fuel related to ros (omits 1000hr fuels): KgC/m2
     real(r8) ::  fuel_frac(nfsc)                                  ! fraction of each litter class in the ros_fuel:-.  
     real(r8) ::  livegrass                                        ! total aboveground grass biomass in patch.  KgC/m2
     real(r8) ::  fuel_bulkd                                       ! average fuel bulk density of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels). KgC/m3
     real(r8) ::  fuel_sav                                         ! average surface area to volume ratio of the ground fuel 
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
     real(r8) ::  frac_burnt                                       ! fraction burnt: frac gridcell/day  
     real(r8) ::  tfc_ros                                          ! total fuel consumed - no trunks.  KgC/m2/day
     real(r8) ::  burnt_frac_litter(nfsc)                          ! fraction of each litter pool burned:-


     ! PLANT HYDRAULICS   (not currently used in hydraulics RGK 03-2018)  
     ! type(ed_patch_hydr_type) , pointer :: pa_hydr              ! All patch hydraulics data, see FatesHydraulicsMemMod.F90


  end type ed_patch_type

  
  !************************************
  !** Resources management type      **
  ! YX
  !************************************
  type, public :: ed_resources_management_type
    
     real(r8) ::  trunk_product_site                       ! Actual  trunk product at site level KgC/site

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
                                  ! (if carbon most likely exudation, if even active)

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



     ! INDICES 
     real(r8) ::  lat                                          ! latitude:  degrees 
     real(r8) ::  lon                                          ! longitude: degrees 
     
     ! Mass Balance (allocation for each element)

     type(site_massbal_type), pointer :: mass_balance(:)

     ! Flux diagnostics (allocation for each element)

     type(site_fluxdiags_type), pointer :: flux_diags(:)

     ! PHENOLOGY 
     real(r8) ::  grow_deg_days                                ! Phenology growing degree days

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

     real(r8) ::  water_memory(numWaterMem)                             ! last 10 days of soil moisture memory...


     ! FIRE
     real(r8) ::  wind                                         ! daily wind in m/min for Spitfire units 
     real(r8) ::  acc_ni                                       ! daily nesterov index accumulating over time.
     real(r8) ::  fdi                                          ! daily probability an ignition event will start a fire
     real(r8) ::  NF                                           ! daily ignitions in km2
     real(r8) ::  frac_burnt                                   ! fraction of area burnt in this day.

     ! PLANT HYDRAULICS
     type(ed_site_hydr_type), pointer :: si_hydr

     ! Soil Layering

     integer :: nlevsoil                      ! Number of soil layers in this site
     real(r8), allocatable :: zi_soil(:)      ! interface level below a "z" level (m)
                                              ! this contains a zero index for surface.
     real(r8), allocatable :: dz_soil(:)      ! layer thickness (m)
     real(r8), allocatable :: z_soil(:)       ! layer depth (m)
     real(r8), allocatable :: rootfrac_scr(:) ! This is just allocated scratch space to hold
                                              ! root fractions. Since root fractions may be dependant
                                              ! on cohort properties, and we do not want to store this infromation
                                              ! on each cohort, we do not keep root fractions in
                                              ! memory, and instead calculate them on demand.
                                              ! This array is allocated over the number of soil
                                              ! layers for each site, and save allocating deallocating.
                                              ! NOTE: THIS SCRATCH SPACE WOULD NOT BE THREAD-SAFE
                                              ! IF WE FORK ON PATCHES

     
     ! DIAGNOSTICS

     ! TERMINATION, RECRUITMENT, DEMOTION, and DISTURBANCE
     
     real(r8), allocatable :: term_nindivs_canopy(:,:) ! number of canopy individuals that were in cohorts which 
                                                       ! were terminated this timestep, on size x pft
     real(r8), allocatable :: term_nindivs_ustory(:,:) ! number of understory individuals that were in cohorts which 
                                                       ! were terminated this timestep, on size x pft
 
     real(r8) :: term_carbonflux_canopy                ! carbon flux from live to dead pools associated 
                                                       ! with termination mortality, per canopy level
     real(r8) :: term_carbonflux_ustory                ! carbon flux from live to dead pools associated 
                                                       ! with termination mortality, per canopy level    
     real(r8) :: demotion_carbonflux                             ! biomass of demoted individuals from canopy to understory [kgC/ha/day]
     real(r8) :: promotion_carbonflux                            ! biomass of promoted individuals from understory to canopy [kgC/ha/day]
     real(r8) :: imort_carbonflux                                ! biomass of individuals killed due to impact mortality per year. [kgC/ha/day]
     real(r8) :: fmort_carbonflux_canopy                         ! biomass of canopy indivs killed due to fire per year. [gC/m2/sec]
     real(r8) :: fmort_carbonflux_ustory                         ! biomass of understory indivs killed due to fire per year [gC/m2/sec] 

     real(r8) :: recruitment_rate(1:maxpft)            ! number of individuals that were recruited into new cohorts
     real(r8), allocatable :: demotion_rate(:)         ! rate of individuals demoted from canopy to understory per FATES timestep
    
     real(r8), allocatable :: promotion_rate(:)                  ! rate of individuals promoted from understory to canopy per FATES timestep
     
     real(r8), allocatable :: imort_rate(:,:)                    ! rate of individuals killed due to impact mortality per year.  on size x pft array
     

     real(r8), allocatable :: fmort_rate_canopy(:,:)             ! rate of canopy individuals killed due to fire mortality per year.  
                                                                 ! on size x pft array  (1:nlevsclass,1:numpft)
     real(r8), allocatable :: fmort_rate_ustory(:,:)             ! rate of understory individuals killed due to fire mortality per year.  
                                                                 ! on size x pft array  (1:nlevsclass,1:numpft)
    
     real(r8), allocatable :: fmort_rate_cambial(:,:)            ! rate of individuals killed due to fire mortality 
                                                                 ! from cambial damage per year.  on size x pft array
     real(r8), allocatable :: fmort_rate_crown(:,:)              ! rate of individuals killed due to fire mortality 
                                                                 ! from crown damage per year.  on size x pft array

     real(r8), allocatable :: growthflux_fusion(:,:)             ! rate of individuals moving into a given size class bin
     ! due to fusion in a given day. on size x pft array 



     ! Canopy Spread
     real(r8) ::  spread                                          ! dynamic canopy allometric term [unitless]
     
  end type ed_site_type

  ! Make public necessary subroutines and functions
  public :: val_check_ed_vars
  public :: dump_site
  public :: dump_patch
  public :: dump_cohort
  public :: dump_cohort_hydr

  contains


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
     type(ed_cohort_type), pointer          :: currentCohort

     
     ! Check through a registry of variables to check
     
     if ( check_hlm_list(trim(var_aliases),'co_n') ) then

        currentCohort => currentPatch%shortest
        do while(associated(currentCohort))
           call check_var_real(currentCohort%n,'cohort%n',return_code)
           if(.not.(return_code.eq.0)) then
              call dump_patch(currentPatch)
              call dump_cohort(currentCohort)
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
              call dump_cohort(currentCohort)
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
     write(fates_log(),*) 'pa%disturbance_rate   = ',cpatch%disturbance_rate
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
  
  subroutine dump_cohort(ccohort)


     type(ed_cohort_type),intent(in),target :: ccohort
     
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) ' Dumping Cohort Information             '
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) 'co%pft                    = ', ccohort%pft
     write(fates_log(),*) 'co%n                      = ', ccohort%n                         
     write(fates_log(),*) 'co%dbh                    = ', ccohort%dbh                                        
     write(fates_log(),*) 'co%hite                   = ', ccohort%hite
     write(fates_log(),*) 'co%coage                  = ', ccohort%coage
     write(fates_log(),*) 'co%laimemory              = ', ccohort%laimemory
     write(fates_log(),*) 'co%sapwmemory             = ', ccohort%sapwmemory
     write(fates_log(),*) 'co%structmemory           = ', ccohort%structmemory
     
     write(fates_log(),*) 'leaf carbon               = ', ccohort%prt%GetState(leaf_organ,all_carbon_elements) 
     write(fates_log(),*) 'fineroot carbon           = ', ccohort%prt%GetState(fnrt_organ,all_carbon_elements) 
     write(fates_log(),*) 'sapwood carbon            = ', ccohort%prt%GetState(sapw_organ,all_carbon_elements) 
     write(fates_log(),*) 'structural (dead) carbon  = ', ccohort%prt%GetState(struct_organ,all_carbon_elements) 
     write(fates_log(),*) 'storage carbon            = ', ccohort%prt%GetState(store_organ,all_carbon_elements) 
     write(fates_log(),*) 'reproductive carbon       = ', ccohort%prt%GetState(repro_organ,all_carbon_elements) 

     write(fates_log(),*) 'co%lai                    = ', ccohort%lai                         
     write(fates_log(),*) 'co%sai                    = ', ccohort%sai  
     write(fates_log(),*) 'co%g_sb_laweight          = ', ccohort%g_sb_laweight
     write(fates_log(),*) 'co%leaf_cost              = ', ccohort%leaf_cost
     write(fates_log(),*) 'co%canopy_layer           = ', ccohort%canopy_layer
     write(fates_log(),*) 'co%canopy_layer_yesterday = ', ccohort%canopy_layer_yesterday
     write(fates_log(),*) 'co%nv                     = ', ccohort%nv
     write(fates_log(),*) 'co%status_coh             = ', ccohort%status_coh
     write(fates_log(),*) 'co%canopy_trim            = ', ccohort%canopy_trim
     write(fates_log(),*) 'co%excl_weight            = ', ccohort%excl_weight               
     write(fates_log(),*) 'co%prom_weight            = ', ccohort%prom_weight               
     write(fates_log(),*) 'co%size_class             = ', ccohort%size_class
     write(fates_log(),*) 'co%size_by_pft_class      = ', ccohort%size_by_pft_class
     write(fates_log(),*) 'co%coage_class            = ', ccohort%coage_class
     write(fates_log(),*) 'co%coage_by_pft_class     = ', ccohort%coage_by_pft_class
     write(fates_log(),*) 'co%gpp_acc_hold           = ', ccohort%gpp_acc_hold
     write(fates_log(),*) 'co%gpp_acc                = ', ccohort%gpp_acc
     write(fates_log(),*) 'co%gpp_tstep              = ', ccohort%gpp_tstep
     write(fates_log(),*) 'co%npp_acc_hold           = ', ccohort%npp_acc_hold
     write(fates_log(),*) 'co%npp_tstep              = ', ccohort%npp_tstep
     write(fates_log(),*) 'co%npp_acc                = ', ccohort%npp_acc
     write(fates_log(),*) 'co%resp_tstep             = ', ccohort%resp_tstep
     write(fates_log(),*) 'co%resp_acc               = ', ccohort%resp_acc
     write(fates_log(),*) 'co%resp_acc_hold          = ', ccohort%resp_acc_hold
     write(fates_log(),*) 'co%rdark                  = ', ccohort%rdark
     write(fates_log(),*) 'co%resp_m                 = ', ccohort%resp_m
     write(fates_log(),*) 'co%resp_g                 = ', ccohort%resp_g
     write(fates_log(),*) 'co%livestem_mr            = ', ccohort%livestem_mr
     write(fates_log(),*) 'co%livecroot_mr           = ', ccohort%livecroot_mr
     write(fates_log(),*) 'co%froot_mr               = ', ccohort%froot_mr
     write(fates_log(),*) 'co%dmort                  = ', ccohort%dmort
     write(fates_log(),*) 'co%treelai                = ', ccohort%treelai
     write(fates_log(),*) 'co%treesai                = ', ccohort%treesai
     write(fates_log(),*) 'co%c_area                 = ', ccohort%c_area
     write(fates_log(),*) 'co%cmort                  = ', ccohort%cmort
     write(fates_log(),*) 'co%bmort                  = ', ccohort%bmort
     write(fates_log(),*) 'co%smort                  = ', ccohort%smort
     write(fates_log(),*) 'co%asmort                 = ', ccohort%asmort
     write(fates_log(),*) 'co%hmort                  = ', ccohort%hmort
     write(fates_log(),*) 'co%frmort                 = ', ccohort%frmort
     write(fates_log(),*) 'co%asmort                 = ', ccohort%asmort
     write(fates_log(),*) 'co%isnew                  = ', ccohort%isnew
     write(fates_log(),*) 'co%dndt                   = ', ccohort%dndt
     write(fates_log(),*) 'co%dhdt                   = ', ccohort%dhdt
     write(fates_log(),*) 'co%ddbhdt                 = ', ccohort%ddbhdt
     write(fates_log(),*) 'co%dbdeaddt               = ', ccohort%dbdeaddt
     write(fates_log(),*) 'co%fraction_crown_burned  = ', ccohort%fraction_crown_burned
     write(fates_log(),*) 'co%fire_mort              = ', ccohort%fire_mort
     write(fates_log(),*) 'co%crownfire_mort         = ', ccohort%crownfire_mort
     write(fates_log(),*) 'co%cambial_mort           = ', ccohort%cambial_mort
     write(fates_log(),*) 'co%size_class             = ', ccohort%size_class
     write(fates_log(),*) 'co%size_by_pft_class      = ', ccohort%size_by_pft_class
     if (associated(ccohort%co_hydr) ) then
        call dump_cohort_hydr(ccohort)
     endif 
     write(fates_log(),*) '----------------------------------------'
     return
  end subroutine dump_cohort

  ! =====================================================================================
  
  subroutine dump_cohort_hydr(ccohort)


     type(ed_cohort_type),intent(in),target :: ccohort
     type(ed_cohort_hydr_type), pointer :: ccohort_hydr
     ccohort_hydr => ccohort%co_hydr
     
     write(fates_log(),*) '--------------------------------------------'
     write(fates_log(),*) ' Dumping Cohort Plant Hydraulic Information '
     write(fates_log(),*) 'ccohort_hydr%th_aroot(:) = ', ccohort_hydr%th_aroot(:)
     write(fates_log(),*) 'ccohort_hydr%v_aroot_layer_init(:) = ', ccohort_hydr%v_aroot_layer_init(:)
     write(fates_log(),*) 'ccohort_hydr%v_aroot_layer(:) = ', ccohort_hydr%v_aroot_layer(:)
     write(fates_log(),*) '--------------------------------------------'
     return
  end subroutine dump_cohort_hydr

end module EDTypesMod
