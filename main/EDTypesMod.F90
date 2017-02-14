module EDTypesMod

  use shr_kind_mod , only : r8 => shr_kind_r8;
  use decompMod    , only : bounds_type 
  use clm_varpar   , only : nlevgrnd, mxpft
  use domainMod    , only : domain_type
  use shr_sys_mod  , only : shr_sys_flush

  implicit none
  save

  !SWITCHES THAT ARE READ IN
  integer         RESTART                                  ! restart flag, 1= read initial system state 0 = bare ground

  ! MODEL PARAMETERS
  real(r8)            :: timestep_secs                     ! subdaily timestep in seconds (e.g. 1800 or 3600) 
  
  real(r8), parameter :: AREA                 = 10000.0_r8 ! Notional area of simulated forest m2
  integer  doy

  integer, parameter  :: invalidValue         = -9999      ! invalid value for gcells,
  ! cohorts, and patches

  ! for setting number of patches per gridcell and number of cohorts per patch
  ! for I/O and converting to a vector

  integer, parameter :: maxPatchesPerCol      = 10          !
  integer, parameter :: maxCohortsPerPatch    = 160         !
  integer, parameter :: cohorts_per_col       = 1600        ! This is the max number of individual items one can store per 

                                                           ! each grid cell and effects the striding in the ED restart 
                                                           ! data as some fields are arrays where each array is
                                                           ! associated with one cohort

  integer, parameter :: numWaterMem           = 10         ! watermemory saved as site level var

  ! BIOLOGY/BIOGEOCHEMISTRY        
  integer , parameter :: INTERNAL_RECRUITMENT = 1          ! internal recruitment fla  1=yes  
  integer , parameter :: EXTERNAL_RECRUITMENT = 0          ! external recruitment flag 1=yes  
  integer , parameter :: SENES                = 10         ! Window of time over which we track temp for cold sensecence (days)
  real(r8), parameter :: DINC_ED              = 1.0_r8     ! size of LAI bins. 
  integer , parameter :: N_DIST_TYPES         = 2          ! number of disturbance types (mortality, fire)
  integer , parameter :: numpft_ed            = 2          ! number of PFTs used in ED. 
  integer , parameter :: maxPft               = 79         ! max number of PFTs potentially used by CLM 


  ! SPITFIRE     
  integer , parameter :: NLSC                 = 6          ! number carbon compartments in above ground litter array 
  integer , parameter :: NFSC                 = 6          ! number fuel size classes  
  integer , parameter :: N_EF                 = 7          ! number of emission factors. One per trace gas or aerosol species.
  integer,  parameter :: NCWD                 = 4          ! number of coarse woody debris pools
  integer,  parameter :: lg_sf                = 6          ! array index of live grass pool for spitfire
  integer,  parameter :: dg_sf                = 1          ! array index of dead grass pool for spitfire
  integer,  parameter :: tr_sf                = 5          ! array index of dead trunk pool for spitfire
  integer,  parameter :: lb_sf                = 4          ! array index of lrge branch pool for spitfire 
  real(r8), parameter :: fire_threshold       = 35.0_r8    ! threshold for fires that spread or go out. KWm-2

  ! COHORT FUSION          
  real(r8), parameter :: FUSETOL              = 0.05_r8     ! min fractional difference in dbh between cohorts

  ! PATCH FUSION 
  real(r8), parameter :: patchfusion_profile_tolerance = 0.05_r8  ! minimum fraction in difference in profiles between patches
  real(r8), parameter :: NTOL                 = 0.05_r8    ! min plant density for hgt bin to be used in height profile comparisons 
  real(r8), parameter :: HITEMAX              = 30.0_r8    ! max dbh value used in hgt profile comparison 
  real(r8), parameter :: DBHMAX               = 150.0_r8   ! max dbh value used in hgt profile comparison 
  integer , parameter :: N_HITE_BINS          = 60         ! no. of hite bins used to distribute LAI
  integer , parameter :: N_DBH_BINS           = 5          ! no. of dbh bins used when comparing patches


  real(r8), parameter :: min_npm2       = 1.0d-5   ! minimum cohort number density per m2 before termination
  real(r8), parameter :: min_patch_area = 0.001_r8 ! smallest allowable patch area before termination
  real(r8), parameter :: min_nppatch    = 1.0d-8   ! minimum number of cohorts per patch (min_npm2*min_patch_area)
  real(r8), parameter :: min_n_safemath = 1.0d-15  ! in some cases, we want to immediately remove super small
                                                   ! number densities of cohorts to prevent FPEs, this is usually
                                                   ! just relevant in the first day after recruitment

  character*4 yearchar                    

  !the lower limit of the size classes of ED cohorts
  !0-10,10-20...
  integer, parameter :: nlevsclass_ed = 13    ! Number of dbh size classes for size structure analysis
                                              ! |0-1,1-2,2-3,3-4,4-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80,80-90,90-100,100+|
!  real(r8), parameter, dimension(16) ::  sclass_ed  = (/0.0_r8,1.0_r8,2.0_r8,3.0_r8,4.0_r8,5.0_r8,10.0_r8,20.0_r8,30.0_r8,40.0_r8, &
!                                                       50.0_r8,60.0_r8,70.0_r8,80.0_r8,90.0_r8,100.0_r8/)

  real(r8), parameter, dimension(nlevsclass_ed) ::  sclass_ed  = (/0.0_r8,5.0_r8,10.0_r8,15.0_r8,20.0_r8,30.0_r8,40.0_r8, &
                                                       50.0_r8,60.0_r8,70.0_r8,80.0_r8,90.0_r8,100.0_r8/)

  integer, parameter :: nlevage_ed = 7  ! Number of patch-age classes for age structured analyses
  real(r8), parameter, dimension(nlevage_ed) ::  ageclass_ed  = (/0.0_r8,1.0_r8,2._r8,5.0_r8,10.0_r8,20.0_r8,50.0_r8/)
  

 !  integer, parameter :: nlevsclass_ed = 17
 !  real(r8), parameter, dimension(17) ::  sclass_ed  = (/0.1_r8, 5.0_r8,10.0_r8,15.0_r8,20.0_r8,25.0_r8, & 
 !                                                       30.0_r8,35.0_r8,40.0_r8,45.0_r8,50.0_r8,55.0_r8, &
 !                                                       60.0_r8,65.0_r8,70.0_r8,75.0_r8,80.0_r8/)

  integer, parameter :: nlevmclass_ed = 5      ! nlev "mortality" classes in ED
                                               ! Number of ways to die
                                               ! (background,hydraulic,carbon,impact,fire)

  character(len = 10), parameter,dimension(5) :: char_list = (/"background","hydraulic ","carbon    ","impact    ","fire      "/)


  ! These vectors are used for history output mapping
  real(r8) ,allocatable :: levsclass_ed(:) ! The lower bound on size classes for ED trees. This 
                                           ! is used really for IO into the
                                           ! history tapes. It gets copied from
                                           ! the parameter array sclass_ed.
  integer , allocatable :: pft_levscpf_ed(:)
  integer , allocatable :: scls_levscpf_ed(:) 
  real(r8), allocatable :: levage_ed(:) 
  integer , allocatable :: levpft_ed(:) 

  
  ! Control Parameters (cp_)            
  ! -------------------------------------------------------------------------------------

  ! These parameters are dictated by FATES internals
  
  integer, parameter :: cp_nclmax = 2       ! Maximum number of canopy layers

  integer, parameter :: cp_nlevcan = 40     ! number of leaf layers in canopy layer
  
  integer, parameter :: cp_maxSWb = 2       ! maximum number of broad-bands in the
                                            ! shortwave spectrum cp_numSWb <= cp_maxSWb
                                            ! this is just for scratch-array purposes
                                            ! if cp_numSWb is larger than this value
                                            ! simply bump this number up as needed
  ! These parameters are dictated by the host model or driver
  
  integer :: cp_numSWb       ! Number of broad-bands in the short-wave radiation
                             ! specturm to track 
                             ! (typically 2 as a default, VIS/NIR, in ED variants <2016)

  integer :: cp_numlevgrnd   ! Number of ground layers
  integer :: cp_numlevsoil   ! Number of soil layers

  ! Number of GROUND layers for the purposes of biogeochemistry; can be either 1 
  ! or the total number of soil layers (includes bedrock)
  integer :: cp_numlevdecomp_full  

  ! Number of SOIL layers for the purposes of biogeochemistry; can be either 1 
  ! or the total number of soil layers
  integer :: cp_numlevdecomp

  ! This character string passed by the HLM is used during the processing of IO
  ! data, so that FATES knows which IO variables it should prepare.  For instance
  ! ATS, ALM and CLM will only want variables specficially packaged for them.
  ! This string will dictate which filter is enacted.
  character(len=16) :: cp_hlm_name

  ! This value can be flushed to history diagnostics, such that the
  ! HLM will interpret that the value should not be included in the average.
  real(r8) :: cp_hio_ignore_val


  ! Is this the master processor, typically useful for knowing if 
  ! the current machine should be printing out messages to the logs or terminals
  ! 1 = TRUE (is master) 0 = FALSE (is not master)
  integer :: cp_masterproc


  ! Module switches (this will be read in one day)
  ! This variable only exists now to serve as a place holder
  !!!!!!!!!! THIS SHOULD NOT BE SET TO TRUE !!!!!!!!!!!!!!!!!
  logical,parameter :: use_fates_plant_hydro = .false.

  !************************************
  !** COHORT type structure          **
  !************************************
  type ed_cohort_type

     ! POINTERS
     type (ed_cohort_type) , pointer :: taller   => null()       ! pointer to next tallest cohort     
     type (ed_cohort_type) , pointer :: shorter  => null()       ! pointer to next shorter cohort     
     type (ed_patch_type)  , pointer :: patchptr => null()       ! pointer to patch that cohort is in
     type (ed_site_type)   , pointer :: siteptr  => null()       ! pointer to site that cohort is in

     ! VEGETATION STRUCTURE
     integer  ::  pft                                    ! pft number
     real(r8) ::  n                                      ! number of individuals in cohort per 'area' (10000m2 default)
     real(r8) ::  dbh                                    ! dbh: cm
     real(r8) ::  hite                                   ! height: meters
     integer  ::  indexnumber                            ! unique number for each cohort. (within clump?)
     real(r8) ::  balive                                 ! total living biomass: kGC per indiv
     real(r8) ::  bdead                                  ! dead biomass:  kGC per indiv
     real(r8) ::  bstore                                 ! stored carbon: kGC per indiv
     real(r8) ::  laimemory                              ! target leaf biomass- set from previous year: kGC per indiv
     integer  ::  canopy_layer                           ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)
     real(r8) ::  b                                      ! total biomass: kGC per indiv
     real(r8) ::  bsw                                    ! sapwood in stem and roots: kGC per indiv
     real(r8) ::  bl                                     ! leaf biomass: kGC per indiv
     real(r8) ::  br                                     ! fine root biomass: kGC per indiv
     real(r8) ::  lai                                    ! leaf area index of cohort   m2/m2
     real(r8) ::  sai                                    ! stem area index of cohort   m2/m2
     real(r8) ::  gscan                                  ! Stomatal resistance of cohort. 
     real(r8) ::  canopy_trim                            ! What is the fraction of the maximum leaf biomass that we are targeting? :-
     real(r8) ::  leaf_cost                              ! How much does it cost to maintain leaves: kgC/m2/year-1
     real(r8) ::  excl_weight                            ! How much of this cohort is demoted each year, as a proportion of all cohorts:-
     real(r8) ::  prom_weight                            ! How much of this cohort is promoted each year, as a proportion of all cohorts:-
     integer  ::  nv                                     ! Number of leaf layers: -
     integer  ::  status_coh                             ! growth status of plant  (2 = leaves on , 1 = leaves off)
     real(r8) ::  c_area                                 ! areal extent of canopy (m2)
     real(r8) ::  treelai                                ! lai of tree (total leaf area (m2) / canopy area (m2)
     real(r8) ::  treesai                                ! stem area index of tree (total stem area (m2) / canopy area (m2)
     logical  ::  isnew                                  ! flag to signify a new cohort, new cohorts have not experienced
                                                         ! npp or mortality and should therefore not be fused or averaged
     integer  ::  size_class                             ! An index that indicates which diameter size bin the cohort currently resides in
                                                         ! this is used for history output. We maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls
     integer  ::  size_by_pft_class                      ! An index that indicates the cohorts position of the joint size-class x functional
                                                         ! type classification. We also maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls


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

     ! Net Primary Production Partitions

     real(r8) ::  npp_leaf                               ! NPP into leaves (includes replacement of turnover):  KgC/indiv/day
     real(r8) ::  npp_froot                              ! NPP into fine roots (includes replacement of turnover):  KgC/indiv/day
     real(r8) ::  npp_bsw                                ! NPP into sapwood: KgC/indiv/day
     real(r8) ::  npp_bdead                              ! NPP into deadwood (structure):  KgC/indiv/day
     real(r8) ::  npp_bseed                              ! NPP into seeds: KgC/indiv/day
     real(r8) ::  npp_store                              ! NPP into storage: KgC/indiv/day

     real(r8) ::  ts_net_uptake(cp_nlevcan)              ! Net uptake of leaf layers: kgC/m2/s
     real(r8) ::  year_net_uptake(cp_nlevcan)            ! Net uptake of leaf layers: kgC/m2/year

     ! RESPIRATION COMPONENTS
     real(r8) ::  rdark                                  ! Dark respiration: kgC/indiv/s
     real(r8) ::  resp_g                                 ! Growth respiration:  kgC/indiv/timestep
     real(r8) ::  resp_m                                 ! Maintenance respiration:  kgC/indiv/timestep 
     real(r8) ::  livestem_mr                            ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (Above ground)
     real(r8) ::  livecroot_mr                           ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (below ground)
     real(r8) ::  froot_mr                               ! Live fine root   maintenance respiration: kgC/indiv/s

     ! ALLOCATION
     real(r8) ::  md                                     ! plant maintenance demand: kgC/indiv/year
     real(r8) ::  leaf_md                                ! leaf  maintenance demand: kgC/indiv/year
     real(r8) ::  root_md                                ! root  maintenance demand: kgC/indiv/year
     real(r8) ::  carbon_balance                         ! carbon remaining for growth and storage: kg/indiv/year
     real(r8) ::  seed_prod                              ! reproduction seed and clonal: KgC/indiv/year
     real(r8) ::  leaf_litter                            ! leaf litter from phenology: KgC/m2
     real(r8) ::  woody_turnover                         ! amount of wood lost each day: kgC/indiv/year. Currently set to zero.

     !MORTALITY
     real(r8) ::  dmort                                  ! proportional mortality rate. (year-1)

     ! Mortality Rate Partitions
     real(r8) ::  bmort                                  ! background mortality rate        n/year
     real(r8) ::  cmort                                  ! carbon starvation mortality rate n/year
     real(r8) ::  hmort                                  ! hydraulic failure mortality rate n/year
     real(r8) ::  imort                                  ! mortality from impacts by others n/year
     real(r8) ::  fmort                                  ! fire mortality                   n/year

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
     real(r8) ::  dbalivedt                              ! time derivative of total living biomass : KgC/year
     real(r8) ::  dbdeaddt                               ! time derivative of dead biomass         : KgC/year
     real(r8) ::  dbstoredt                              ! time derivative of stored biomass       : KgC/year
     real(r8) ::  storage_flux                           ! flux from npp into bstore               : KgC/year

     ! FIRE
     real(r8) ::  cfa                                    ! proportion of crown affected by fire:-
     real(r8) ::  cambial_mort                           ! probability that trees dies due to cambial char:-
     real(r8) ::  crownfire_mort                         ! probability of tree post-fire mortality due to crown scorch:-
     real(r8) ::  fire_mort                              ! post-fire mortality from cambial and crown damage assuming two are independent:-

  end type ed_cohort_type

  !************************************
  !** Patch type structure           **
  !************************************

  type ed_patch_type

     ! POINTERS
     type (ed_cohort_type), pointer :: tallest => null()           ! pointer to patch's tallest cohort    
     type (ed_cohort_type), pointer :: shortest => null()          ! pointer to patch's shortest cohort
     type (ed_patch_type),  pointer :: older => null()             ! pointer to next older patch   
     type (ed_patch_type),  pointer :: younger => null()           ! pointer to next younger patch      
     type (ed_site_type),   pointer :: siteptr => null()           ! pointer to the site that the patch is in

     !INDICES
     integer  :: patchno                                           ! unique number given to each new patch created for tracking

     ! INTERF-TODO: THIS VARIABLE SHOULD BE REMOVED
     integer  :: clm_pno                                           ! clm patch number (index of p vector)

     ! PATCH INFO
     real(r8) ::  age                                              ! average patch age: years                   
     integer  ::  age_class                                        ! age class of the patch for history binning purposes
     real(r8) ::  area                                             ! patch area: m2  
     integer  ::  countcohorts                                     ! Number of cohorts in patch
     integer  ::  ncl_p                                            ! Number of occupied canopy layers

     ! LEAF ORGANIZATION
     real(r8) ::  spread(cp_nclmax)                                   ! dynamic ratio of dbh to canopy area: cm/m2
     real(r8) ::  pft_agb_profile(numpft_ed,n_dbh_bins)            ! binned above ground biomass, for patch fusion: KgC/m2
     real(r8) ::  canopy_layer_lai(cp_nclmax)                         ! lai that is shading this canopy layer: m2/m2 
     real(r8) ::  total_canopy_area                                ! area that is covered by vegetation : m2
     real(r8) ::  total_tree_area                                  ! area that is covered by woody vegetation : m2
     real(r8) ::  canopy_area                                      ! area that is covered by vegetation : m2 (is this different to total_canopy_area?
     real(r8) ::  bare_frac_area                                   ! bare soil in this patch expressed as a fraction of the total soil surface.
     real(r8) ::  lai                                              ! leaf area index of patch

     real(r8) ::  tlai_profile(cp_nclmax,numpft_ed,cp_nlevcan)        ! total   leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  elai_profile(cp_nclmax,numpft_ed,cp_nlevcan)        ! exposed leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  tsai_profile(cp_nclmax,numpft_ed,cp_nlevcan)        ! total   stem area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  esai_profile(cp_nclmax,numpft_ed,cp_nlevcan)        ! exposed stem area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  layer_height_profile(cp_nclmax,numpft_ed,cp_nlevcan)
     real(r8) ::  canopy_area_profile(cp_nclmax,numpft_ed,cp_nlevcan) ! fraction of canopy in each canopy 
     ! layer, pft, and leaf layer:-
     integer  ::  present(cp_nclmax,numpft_ed)                        ! is there any of this pft in this canopy layer?      
     integer  ::  nrad(cp_nclmax,numpft_ed)                           ! number of exposed leaf layers for each canopy layer and pft
     integer  ::  ncan(cp_nclmax,numpft_ed)                           ! number of total   leaf layers for each canopy layer and pft

     !RADIATION FLUXES      
     real(r8) ::  fabd_sun_z(cp_nclmax,numpft_ed,cp_nlevcan)          ! sun fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabd_sha_z(cp_nclmax,numpft_ed,cp_nlevcan)          ! shade fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sun_z(cp_nclmax,numpft_ed,cp_nlevcan)          ! sun fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sha_z(cp_nclmax,numpft_ed,cp_nlevcan)          ! shade fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-

     real(r8) ::  ed_laisun_z(cp_nclmax,numpft_ed,cp_nlevcan)         ! amount of LAI in the sun   in each canopy layer, 
     ! pft, and leaf layer. m2/m2
     real(r8) ::  ed_laisha_z(cp_nclmax,numpft_ed,cp_nlevcan)         ! amount of LAI in the shade in each canopy layer,
     real(r8) ::  ed_parsun_z(cp_nclmax,numpft_ed,cp_nlevcan)         ! PAR absorbed  in the sun   in each canopy layer,
     real(r8) ::  ed_parsha_z(cp_nclmax,numpft_ed,cp_nlevcan)         ! PAR absorbed  in the shade in each canopy layer,
     real(r8) ::  f_sun(cp_nclmax,numpft_ed,cp_nlevcan)               ! fraction of leaves in the sun in each canopy layer, pft, 

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


     !SEED BANK
     real(r8) :: seeds_in(numpft_ed)                               ! seed production KgC/m2/year
     real(r8) :: seed_decay(numpft_ed)                             ! seed decay in KgC/m2/year
     real(r8) :: seed_germination(numpft_ed)                       ! germination rate of seed pool in KgC/m2/year

     ! PHOTOSYNTHESIS       
     real(r8) ::  psn_z(cp_nclmax,numpft_ed,cp_nlevcan)               ! carbon assimilation in each canopy layer, pft, and leaf layer. umolC/m2/s

     ! ROOTS
     real(r8), allocatable ::  rootfr_ft(:,:)                      ! root fraction of each PFT in each soil layer:-
     real(r8), allocatable ::  rootr_ft(:,:)                       ! fraction of water taken from each PFT and soil layer:-
     real(r8) ::  btran_ft(numpft_ed)                              ! btran calculated seperately for each PFT:-   

     ! DISTURBANCE 
     real(r8) ::  disturbance_rates(n_dist_types)                  ! disturbance rate from 1) mortality and 2) fire: fraction/day
     real(r8) ::  disturbance_rate                                 ! larger effective disturbance rate: fraction/day

     ! LITTER AND COARSE WOODY DEBRIS 
     ! Pools of litter (non respiring) 
     real(r8) ::  cwd_ag(ncwd)                                     ! above ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  cwd_bg(ncwd)                                     ! below ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  leaf_litter(numpft_ed)                           ! above ground leaf litter that does not respire. KgC/m2
     real(r8) ::  root_litter(numpft_ed)                           ! below ground fine root litter that does not respire. KgC/m2

     ! Fluxes of litter (non respiring) 
     real(r8) :: fragmentation_scaler                              ! Scale rate of litter fragmentation. 0 to 1.
     real(r8) :: cwd_ag_in(ncwd)                                   ! Flux into CWD_AG from turnover and mortality KgC/m2/y
     real(r8) :: cwd_bg_in(ncwd)                                   ! Flux into cwd_bg from root turnover and mortality KgC/m2/y
     real(r8) :: cwd_ag_out(ncwd)                                  ! Flux out of AG CWD into AG litter KgC/m2/y
     real(r8) :: cwd_bg_out(ncwd)                                  ! Flux out of BG CWD into BG litter KgC/m2/


     real(r8) :: leaf_litter_in(numpft_ed)                         ! Flux in  to AG leaf litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: leaf_litter_out(numpft_ed)                        ! Flux out of AG leaf litter from fragmentation KgC/m2/y
     real(r8) :: root_litter_in(numpft_ed)                         ! Flux in  to BG root litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: root_litter_out(numpft_ed)                        ! Flux out of BG root from fragmentation KgC/m2/y

     ! Derivatives of litter (non respiring) 
     real(r8) ::  dcwd_AG_dt(ncwd)                                 ! rate of change of above ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dcwd_BG_dt(ncwd)                                 ! rate of change of below ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dleaf_litter_dt(numpft_ed)                       ! rate of change of leaf litter in each size class: KgC/m2/year. 
     real(r8) ::  droot_litter_dt(numpft_ed)                       ! rate of change of root litter in each size class: KgC/m2/year. 

     real(r8) ::  canopy_mortality_woody_litter                    ! flux of wood litter in to litter pool: KgC/m2/year
     real(r8) ::  canopy_mortality_leaf_litter(numpft_ed)          ! flux in to  leaf litter from tree death: KgC/m2/year
     real(r8) ::  canopy_mortality_root_litter(numpft_ed)          ! flux in to froot litter  from tree death: KgC/m2/year

     real(r8) ::  repro(numpft_ed)                                 ! allocation to reproduction per PFT : KgC/m2

     !FUEL CHARECTERISTICS
     real(r8) ::  sum_fuel                                         ! total ground fuel related to ros (omits 1000hr fuels): KgC/m2
     real(r8) ::  fuel_frac(ncwd+2)                                ! fraction of each litter class in the ros_fuel:-.  
     real(r8) ::  livegrass                                        ! total aboveground grass biomass in patch.  KgC/m2
     real(r8) ::  fuel_bulkd                                       ! average fuel bulk density of the ground fuel 
     ! (incl. live grasses. omits 1000hr fuels). KgC/m3
     real(r8) ::  fuel_sav                                         ! average surface area to volume ratio of the ground fuel 
     ! (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_mef                                         ! average moisture of extinction factor 
     ! of the ground fuel (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_eff_moist                                   ! effective avearage fuel moisture content of the ground fuel 
     ! (incl. live grasses. omits 1000hr fuels)
     real(r8) ::  litter_moisture(ncwd+2)

     ! FIRE SPREAD
     real(r8) ::  ros_front                                        ! rate of forward  spread of fire: m/min
     real(r8) ::  ros_back                                         ! rate of backward spread of fire: m/min
     real(r8) ::  effect_wspeed                                    ! windspeed modified by fraction of relative grass and tree cover: m/min
     real(r8) ::  tau_l                                            ! Duration of lethal heating: mins
     real(r8) ::  fi                                               ! average fire intensity of flaming front:  kj/m/s or kw/m
     integer  ::  fire                                             ! Is there a fire? 1=yes 0=no
     real(r8) ::  fd                                               ! fire duration: mins
     real(r8) ::  nf                                               ! number of fires initiated daily: n/gridcell/day
     real(r8) ::  sh                                               ! average scorch height: m 

     ! FIRE EFFECTS     
     real(r8) ::  ab                                               ! area burnt:  m2/day
     real(r8) ::  frac_burnt                                       ! fraction burnt: frac gridcell/day  
     real(r8) ::  tfc_ros                                          ! total fuel consumed - no trunks.  KgC/m2/day
     real(r8) ::  burnt_frac_litter(nfsc)                          ! fraction of each litter pool burned:-

   contains

     procedure, public :: set_root_fraction

  end type ed_patch_type

  !************************************
  !** Site type structure           **
  !************************************

  type ed_site_type

     ! POINTERS  
     type (ed_patch_type), pointer :: oldest_patch => null()   ! pointer to oldest patch at the site  
     type (ed_patch_type), pointer :: youngest_patch => null() ! pointer to yngest patch at the site

     ! INDICES 
     real(r8) ::  lat                                          ! latitude:  degrees 
     real(r8) ::  lon                                          ! longitude: degrees 

     ! CARBON BALANCE       
     real(r8) :: flux_in                                      ! for carbon balance purpose. C coming into biomass pool:  KgC/site
     real(r8) :: flux_out                                     ! for carbon balance purpose. C leaving ED pools  KgC/site
     real(r8) :: old_stock                                    ! for accounting purposes, remember biomass stock from last time:  KgC/site
     real(r8) :: npp                                          ! used for calculating NEP and NBP during BGC summarization phase
     real(r8) :: nep                                          ! Net ecosystem production, i.e. fast-timescale carbon balance that 
                                                              ! does not include disturbance [gC/m2/s]
     real(r8) :: nbp                                          ! Net biosphere production, i.e. slow-timescale carbon balance that 
                                                              ! integrates to total carbon change [gC/m2/s]
     real(r8) :: tot_seed_rain_flux                           ! [gC/m2/s] total flux of carbon from seed rain
     real(r8) :: fire_c_to_atm                                ! total fire carbon loss to atmosphere [gC/m2/s]
     real(r8) :: ed_litter_stock                              ! litter in [gC/m2]
     real(r8) :: cwd_stock                                    ! coarse woody debris [gC/m2]
     real(r8) :: biomass_stock                                ! total biomass at the column level in [gC / m2]
     real(r8) :: totfatesc                                    ! Total FATES carbon at the site, including vegetation, CWD, seeds, 
                                                              ! and FATES portion of litter [gC/m2] 
     real(r8) :: totbgcc                                      ! Total BGC carbon at the site, including litter, and soil pools [gC/m2] 
     real(r8) :: totecosysc                                   ! Total ecosystem C at the site, including vegetation, 
                                                              ! CWD, litter (from HLM and FATES), and soil pools [gC/m2]

     real(r8) :: totfatesc_old                                ! Total FATES C at the site from last call to balance check [gC/m2]
     real(r8) :: totbgcc_old                                  ! Total BGC C at the site from last call to balance check [gC/m2] 
     real(r8) :: totecosysc_old                               ! Total ecosystem C at the site from last call to balance check [gC/m2]
     
     real(r8) :: fates_to_bgc_this_ts                         ! total flux of carbon from FATES to BGC models on current timestep [gC/m2/s] 
     real(r8) :: fates_to_bgc_last_ts                         ! total flux of carbon from FATES to BGC models on previous timestep [gC/m2/s] 

     real(r8) :: cbal_err_fates                               ! [gC/m2/s]  total carbon balance error for FATES processes
     real(r8) :: cbal_err_bgc                                 ! [gC/m2/s]  total carbon balance error for BGC (HLM) processes
     real(r8) :: cbal_err_tot                                 ! [gC/m2/s]  total carbon balance error for all land processes

     real(r8) :: nep_timeintegrated                           ! Net ecosystem production accumulated over model time-steps [gC/m2]
     real(r8) :: hr_timeintegrated                            ! Heterotrophic respiration accumulated over model time-steps [gC/m2]
     real(r8) :: npp_timeintegrated                           ! Net primary production accumulated over model time-steps [gC/m2]
     real(r8) :: nbp_integrated                               ! Net biosphere production accumulated over model time-steps [gC/m2]


     ! DISTURBANCE
     real(r8) ::  disturbance_mortality                        ! site level disturbance rates from mortality.
     real(r8) ::  disturbance_fire                             ! site level disturbance rates from fire.  
     integer  ::  dist_type                                    ! disturbance dist_type id.
     real(r8) ::  disturbance_rate                             ! site total dist rate

     ! PHENOLOGY 
     real(r8) ::  ED_GDD_site                                  ! ED Phenology growing degree days.
     integer  ::  status                                       ! are leaves in this pixel on or off for cold decid
     integer  ::  dstatus                                      ! are leaves in this pixel on or off for drought decid
     real(r8) ::  ncd                                          ! no chilling days:-
     real(r8) ::  last_n_days(senes)                           ! record of last 10 days temperature for senescence model. deg C
     integer  ::  leafondate                                   ! doy of leaf on:-
     integer  ::  leafoffdate                                  ! doy of leaf off:-
     integer  ::  dleafondate                                  ! doy of leaf on drought:-
     integer  ::  dleafoffdate                                 ! doy of leaf on drought:-
     real(r8) ::  water_memory(10)                             ! last 10 days of soil moisture memory...

     !SEED BANK
     real(r8) :: seed_bank(numpft_ed)                              ! seed pool in KgC/m2/year
     real(r8) :: dseed_dt(numpft_ed)
     real(r8) :: seed_rain_flux(numpft_ed)                         ! flux of seeds from exterior KgC/m2/year (needed for C balance purposes)

     ! FIRE 
     real(r8) ::  acc_ni                                       ! daily nesterov index accumulating over time.
     real(r8) ::  ab                                           ! daily burnt area: m2
     real(r8) ::  frac_burnt                                   ! fraction of soil burnt in this day.
     real(r8) ::  total_burn_flux_to_atm                       ! total carbon burnt to the atmosphere in this day. KgC/site
     real(r8) ::  cwd_ag_burned(ncwd)
     real(r8) ::  leaf_litter_burned(numpft_ed)

  end type ed_site_type

  !************************************
  !** Userdata type structure       **
  !************************************

!  type userdata
!     integer  ::   cohort_number            ! Counts up the number of cohorts which have been made.
!     integer  ::   n_sub                    ! num of substeps in year 
!     real(r8) ::   deltat                   ! fraction of year used for each timestep (1/N_SUB)
!     integer  ::   time_period              ! Within year timestep (1:N_SUB) day of year
!     integer  ::   restart_year             ! Which year of simulation are we starting in? 
!  end type userdata
!  type(userdata), public, target :: udata   ! THIS WAS NOT THREADSAFE
  !-------------------------------------------------------------------------------------!

  public :: ed_hist_scpfmaps

contains

  !-------------------------------------------------------------------------------------!
  subroutine ed_hist_scpfmaps
    ! This subroutine allocates and populates the variables
    ! that define the mapping of variables in history files in the "scpf" format
    ! back to
    ! its respective size-class "sc" and pft "pf"

    integer :: i
    integer :: isc
    integer :: ipft

    allocate( levsclass_ed(1:nlevsclass_ed   ))
    allocate( pft_levscpf_ed(1:nlevsclass_ed*mxpft))
    allocate(scls_levscpf_ed(1:nlevsclass_ed*mxpft))
    allocate( levpft_ed(1:mxpft   ))
    allocate( levage_ed(1:nlevage_ed   ))

    ! Fill the IO array of plant size classes
    ! For some reason the history files did not like
    ! a hard allocation of sclass_ed
    levsclass_ed(:) = sclass_ed(:)
    
    levage_ed(:) = ageclass_ed(:)

    ! make pft array
    do ipft=1,mxpft
       levpft_ed(ipft) = ipft
    end do

    ! Fill the IO arrays that match pft and size class to their combined array
    i=0
    do ipft=1,mxpft
       do isc=1,nlevsclass_ed
          i=i+1
          pft_levscpf_ed(i) = ipft
          scls_levscpf_ed(i) = isc
       end do
    end do

  end subroutine ed_hist_scpfmaps

  !-------------------------------------------------------------------------------------!
  function map_clmpatch_to_edpatch(site, clmpatch_number) result(edpatch_pointer)
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(in), target :: site
    integer, intent(in) :: clmpatch_number 
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: edpatch_pointer
    !----------------------------------------------------------------------
    
    ! There is a one-to-one mapping between edpatches and clmpatches. To obtain
    ! this mapping - the following is computed elsewhere in the code base
    ! (1) what is the weight respective to the column of clmpatch? 
    !     dynEDMod determines this via the following logic
    !        if (clm_patch%is_veg(p) .or. clm_patch%is_bareground(p)) then
    !           clm_patch%wtcol(p) = clm_patch%wt_ed(p)
    !        else
    !           clm_patch%wtcol(p)  = 0.0_r8 
    !        end if
    ! (2) is the clmpatch active? 
    !     subgridWeightsMod uses the following logic (in routine is_active_p) to determine if
    !     clmpatch_number is active ( this is a shortened version of the logic to capture
    !     only the essential parts relevent here)
    !         if (clmpatch%wtcol(p) > 0._r8) is_active_p = .true.

    edpatch_pointer => site%oldest_patch    
    do while ( clmpatch_number /= edpatch_pointer%clm_pno )
       edpatch_pointer => edpatch_pointer%younger
    end do

  end function map_clmpatch_to_edpatch

  !-------------------------------------------------------------------------------------!
  subroutine set_root_fraction( this , depth_gl)
    !
    ! !DESCRIPTION:
    !  Calculates the fractions of the root biomass in each layer for each pft. 
    !
    ! !USES:
    use PatchType   , only : clmpatch => patch
    use EDPftvarcon   , only : EDPftvarcon_inst
    !
    ! !ARGUMENTS    
    class(ed_patch_type) :: this
    real(r8),intent(in)  :: depth_gl(0:cp_numlevgrnd)
    !
    ! !LOCAL VARIABLES:
    integer :: lev,p,c,ft
    !----------------------------------------------------------------------

    do ft = 1,numpft_ed 
       do lev = 1, cp_numlevgrnd
          this%rootfr_ft(ft,lev) = 0._r8
       enddo

       do lev = 1, cp_numlevsoil-1
          this%rootfr_ft(ft,lev) = .5_r8*( &
                 exp(-EDPftvarcon_inst%roota_par(ft) * depth_gl(lev-1))  &
               + exp(-EDPftvarcon_inst%rootb_par(ft) * depth_gl(lev-1))  &
               - exp(-EDPftvarcon_inst%roota_par(ft) * depth_gl(lev))    &
               - exp(-EDPftvarcon_inst%rootb_par(ft) * depth_gl(lev)))
       end do
    end do

  end subroutine set_root_fraction


  ! =====================================================================================
  
 


end module EDTypesMod
