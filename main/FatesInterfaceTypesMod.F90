module FatesInterfaceTypesMod
  
  use FatesConstantsMod   , only : r8 => fates_r8
  use FatesConstantsMod   , only : itrue,ifalse
  use FatesConstantsMod   , only : fates_unset_int
  use FatesGlobals        , only : fates_global_verbose
  use FatesGlobals        , only : fates_log
  use FatesGlobals        , only : endrun => fates_endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)

  use FatesInterfaceVariableTypeMod, only : fates_interface_variable_type
  
  implicit none

   private        ! By default everything is private

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   ! -------------------------------------------------------------------------------------
   ! Parameters that are dictated by the Host Land Model
   ! THESE ARE NOT DYNAMIC. SHOULD BE SET ONCE DURING INTIALIZATION.
   ! -------------------------------------------------------------------------------------

  
   integer, public :: hlm_numSWb  ! Number of broad-bands in the short-wave radiation
                                             ! specturm to track 
                                             ! (typically 2 as a default, VIS/NIR, in ED variants <2016)

   integer, public :: hlm_ivis    ! The HLMs assumption of the array index associated with the 
                                             ! visible portion of the spectrum in short-wave radiation arrays

   integer, public :: hlm_inir    ! The HLMs assumption of the array index associated with the 
                                             ! NIR portion of the spectrum in short-wave radiation arrays

   integer, public :: hlm_maxlevsoil   ! Max number of soil layers


   integer, public :: hlm_is_restart   ! Is the HLM signalling that this is a restart
                                                  ! type simulation?
                                                  ! 1=TRUE, 0=FALSE
   
   character(len=16), public :: hlm_name ! This character string passed by the HLM
                                                    ! is used during the processing of IO data, 
                                                    ! so that FATES knows which IO variables it 
                                                    ! should prepare.  For instance
                                                    ! ATS, ALM and CLM will only want variables 
                                                    ! specficially packaged for them.
                                                    ! This string sets which filter is enacted.

   character(len=16), public :: hlm_decomp ! This string defines which soil decomposition
                                           ! scheme is active
                                           ! expected values are one of CENTURY,MIMICS,CTC
   
   
   character(len=16), public :: hlm_nu_com ! This string defines which soil
                                                      ! nutrient competition scheme is in use.
                                                      ! current options with
                                                      ! E3SM: RD, ECA
                                                      ! CESM: NONE
                                                      ! ATS: ?
                                                      ! NORESM: ?
   

   integer, public :: hlm_nitrogen_spec   ! This flag signals which nitrogen
                                                     ! species are active if any:
                                                     ! 0: none
                                                     ! 1: nh4 only
                                                     ! 2: nh4 and no3

   integer, public :: hlm_phosphorus_spec ! Signals if phosphorous is turned on in the HLM
                                                     ! 0: none
                                                     ! 1: p is on

   real(r8), public :: hlm_stepsize        ! The step-size of the host land model (s)
                                           ! moreover, this is the shortest main-model timestep
                                           ! at which fates will be called on the main model integration loop
  
   real(r8), public :: hlm_hio_ignore_val  ! This value can be flushed to history 
                                                      ! diagnostics, such that the
                                                      ! HLM will interpret that the value should not 
                                                      ! be included in the average.
   
   integer, public :: hlm_masterproc  ! Is this the master processor, typically useful
                                                 ! for knowing if the current machine should be 
                                                 ! printing out messages to the logs or terminals
                                                 ! 1 = TRUE (is master) 0 = FALSE (is not master)

   integer, public :: hlm_ipedof      ! The HLM pedotransfer index
                                                 ! this is only used by the plant hydraulics
                                                 ! submodule to check and/or enable consistency
                                                 ! between the pedotransfer functions of the HLM
                                                 ! and how it moves and stores water in its
                                                 ! rhizosphere shells

   integer, public :: hlm_parteh_mode   ! This flag signals which Plant Allocation and Reactive
                                                   ! Transport (exensible) Hypothesis (PARTEH) to use

   integer, public :: hlm_seeddisp_cadence ! This flag signals at what cadence to disperse seeds across gridcells
                                                   ! 0 => no seed dispersal
                                                   ! 1, 2, 3 => daily, monthly, yearly dispersal

   integer, public :: hlm_use_ch4       ! This flag signals whether the methane model in ELM/CLM is
                                        ! active, and therefore whether or not boundary conditions
                                        ! need to be prepped
   
   integer, public :: hlm_use_vertsoilc ! This flag signals whether or not the 
                                                   ! host model is using vertically discretized
                                                   ! soil carbon
                                                   ! 1 = TRUE,  0 = FALSE
   
   integer, public :: hlm_spitfire_mode  ! Flag to signal SPITFIRE mode
                                         ! See namelist_definition_clm4_5.xml
                                         ! ignitions: 1=constant, >1=external data sources (lightning and/or anthropogenic)

   integer, public :: hlm_use_managed_fire    ! Flag to enable managed fire mode.  Requires spitfire to be on.

   integer, public :: hlm_use_lu_harvest      ! This flag signals whether or not to use
                                                         ! harvest data from the hlm
                                                         ! 0 = do not use lu harvest from hlm
                                                         ! 1 = use lu harvest from hlm  
                                                         ! If 1, it automatically sets
                                                         ! hlm_use_logging to 1

   integer, public :: hlm_num_lu_harvest_cats    ! number of hlm harvest categories (e.g. primary forest harvest, secondary young forest harvest, etc.)
                                                         ! this is the first dimension of:
                                                         ! harvest_rates in dynHarvestMod
                                                         ! bc_in%hlm_harvest_rates and bc_in%hlm_harvest_catnames

   integer, public :: hlm_use_luh                   ! flag to signal whether or not to use luh2 drivers
   integer, public :: hlm_use_potentialveg          ! flag to signal whether or not to use potential vegetation only
                                                    ! (i.e., no land use and instead force all lands to be primary)
   integer, public :: hlm_num_luh2_states           ! number of land use state types provided in LUH2 forcing dataset

   integer, public :: hlm_num_luh2_transitions      ! number of land use transition types provided in LUH2 forcing dataset

   integer, public :: hlm_sf_nofire_def               ! Definition of a no-fire case for hlm_spitfire_mode
   integer, public :: hlm_sf_scalar_lightning_def     ! Definition of a scalar-lightning case for hlm_spitfire_mode
   integer, public :: hlm_sf_successful_ignitions_def ! Definition of a successful-ignition dataset case for hlm_spitfire_mode
   integer, public :: hlm_sf_anthro_ignitions_def      ! Definition of an anthropogenic-ignition dataset case for hlm_spitfire_mode
   

   integer, public :: hlm_use_logging       ! This flag signals whether or not to use
                                                       ! the logging module
                                                         ! If hlm_use_lu_harvest is zero,
                                                         ! then logging is determined by
                                                         ! the fates parameter file
                                                         ! If hlm_use_lu_harvest is non-zero,
                                                         ! then this flag is automatically
                                                         ! set to 1 and logging is determined
                                                         ! by the lu harvest input from the hlm

   integer, public :: hlm_use_planthydro    ! This flag signals whether or not to use
                                                       ! plant hydraulics (bchristo/xu methods)
                                                       ! 1 = TRUE, 0 = FALSE
                                                       ! THIS IS CURRENTLY NOT SUPPORTED 

   integer, public :: hlm_use_cohort_age_tracking ! This flag signals whether or not to use
                                                  ! cohort age tracking. 1 = TRUE, 0 = FALSE


   integer, public :: hlm_use_tree_damage         ! This flag signals whether or not to turn on the
                                                  ! tree damage module

   integer, public :: hlm_hydr_solver             ! Switch that defines which hydraulic solver to use
                                                  ! 1 = Taylor solution that solves plant fluxes with 1 layer
                                                  !     sequentially placing solution on top of previous layer solves
                                                  ! 2 = Picard solution that solves all fluxes in a plant and
                                                  !     the soil simultaneously, 2D: soil x (root + shell)
                                                  ! 3 = Newton-Raphson (Deprecated) solution that solves all fluxes in a plant and
                                                  !     the soil simultaneously, 2D: soil x (root + shell)

   integer, public :: hlm_maintresp_leaf_model    ! switch for choosing between leaf maintenance
                                                  ! respiration model. 1=Ryan (1991), 2=Atkin et al (2017)

   integer, public :: hlm_mort_cstarvation_model  ! Switch for carbon starvation mortality:
                                                  ! 1 -- Linear model
                                                  ! 2 -- Exponential model

   integer, public :: hlm_radiation_model         ! Switch for radiation model
                                                  ! Norman (1) and Two-stream (2)

   integer, public :: hlm_regeneration_model      ! Switch for choosing between regeneration models:
                                                  ! (1) for Fates default
                                                  ! (2) for the Tree Recruitment Scheme (Hanbury-Brown et al., 2022)
                                                  ! (3) for the Tree Recruitment Scheme without seedling dynamics

   integer, public :: hlm_use_ed_st3              ! This flag signals whether or not to use
                                                  ! (ST)atic (ST)and (ST)ructure mode (ST3)
                                                  ! Essentially, this gives us the ability
                                                  ! to turn off "dynamics", ie growth, disturbance
                                                  ! recruitment and mortality.
                                                  ! (EXPERIMENTAL!!!!! - RGK 07-2017)
                                                  ! 1 = TRUE, 0 = FALSE
                                                  ! default should be FALSE (dynamics on)
                                                  ! cannot be true with prescribed_phys

   integer, public :: hlm_use_ed_prescribed_phys ! This flag signals whether or not to use
                                                            ! prescribed physiology, somewhat the opposite
                                                            ! to ST3, in this case can turn off
                                                            ! fast processes like photosynthesis and respiration
                                                            ! and prescribe NPP
                                                            ! (NOT CURRENTLY IMPLEMENTED - PLACEHOLDER)
                                                            ! 1 = TRUE, 0 = FALSE
                                                            ! default should be FALSE (biophysics on)
                                                            ! cannot be true with st3 mode

   integer, public :: hlm_use_inventory_init     ! Initialize this simulation from
                                                            ! an inventory file. If this is toggled on
                                                            ! an inventory control file must be specified
                                                            ! as well.
                                                            ! 1 = TRUE, 0 = FALSE
   
   character(len=256), public :: hlm_inventory_ctrl_file ! This is the full path to the
                                                                    ! inventory control file that
                                                                    ! specifieds the availabel inventory datasets
                                                                    ! there locations and their formats
                                                                    ! This need only be defined when
                                                                    ! hlm_use_inventory_init = 1

  integer, public ::  hlm_use_fixed_biogeog                         !  Flag to use FATES fixed biogeography mode
                                                                    !  1 = TRUE, 0 = FALSE 

  integer, public ::  hlm_use_nocomp                                !  Flag to use FATES no competition mode
                                                                    !  1 = TRUE, 0 = FALSE

  integer, public ::  hlm_use_sp                                    !  Flag to use FATES satellite phenology (LAI) mode
                                                                    !  1 = TRUE, 0 = FALSE

  
  ! Flag specifying what types of history fields to allocate and prepare
  ! The "_dynam" refers to history fields that can be updated on the dynamics (daily) step
  ! THe "_hifrq" refers to history fields that can be updated on the model (high-frequency) step
  ! 0 = no output
  ! 1 = site-level averages only
  ! 2 = allow the second dimension
  
  integer, public :: hlm_hist_level_dynam                           
                                                                    
  integer, public :: hlm_hist_level_hifrq
  
   ! -------------------------------------------------------------------------------------
   ! Parameters that are dictated by FATES and known to be required knowledge
   !  needed by the HLMs
   ! -------------------------------------------------------------------------------------

   ! Variables mostly used for dimensioning host land model (HLM) array spaces
   
   integer, public :: fates_maxElementsPerPatch ! maxElementsPerPatch is the value that is ultimately
                                                           ! used to set the size of the largest arrays necessary
                                                           ! in things like restart files (probably hosted by the 
                                                           ! HLM). The size of these arrays are not a parameter
                                                           ! because it is simply the maximum of several different
                                                           ! dimensions. It is possible that this would be the
                                                           ! maximum number of cohorts per patch, but
                                                           ! but it could be other things.

   integer, public :: fates_maxElementsPerSite  ! This is the max number of individual items one can store per 
                                                           ! each grid cell and effects the striding in the ED restart 
                                                           ! data as some fields are arrays where each array is
                                                           ! associated with one cohort

   ! The number of patches that FATES wants the HLM to allocate
   ! for non sp/no-comp, this is the same number of patches that
   ! fates tracks, plus 1 for the bare-ground. For sp/nocomp, it is the
   ! maximum between the number fates tracks, and the numper of PFTs+CFTs
   ! in the surface file.
   ! for an SP simulation, if there are more PFTs and CFTs in the surface
   ! dataset than the number of PFTs in FATES, we have to allocate with
   ! the prior so that we can hold the LAI data
   integer, public :: fates_maxPatchesPerSite

   integer, public :: max_comp_per_site         ! This is the maximum number of nutrient aquisition
                                                           ! competitors that will be generated on each site
   
   
   integer, public :: fates_dispersal_kernel_mode   ! Flag to signal the type of kernel used for grid cell seed dispersal

   integer, parameter, public :: fates_dispersal_kernel_exponential = 1  ! exponential dispersal kernel
   integer, parameter, public :: fates_dispersal_kernel_exppower = 2     ! exponential power (ExP) dispersal kernel
   integer, parameter, public :: fates_dispersal_kernel_logsech = 3      ! logistic-sech (LogS) dispersal kernel

   integer, parameter, public :: fates_dispersal_cadence_none = 0     ! no dispersal (use seed rain only)
   integer, parameter, public :: fates_dispersal_cadence_daily = 1    ! Disperse seeds daily
   integer, parameter, public :: fates_dispersal_cadence_monthly = 2  ! Disperse seeds monthly
   integer, parameter, public :: fates_dispersal_cadence_yearly = 3   ! Disperse seeds yearly
   
   integer, parameter :: hlm_subgrid_levels = 5  ! The number of subgrid hierarchy levels that the HLM
                                                 ! Including the gridcell level, ELM = 5, CLM = 4

   ! Subgrid levels for HLM-FATES interface variable
   integer, parameter, public :: subgrid_gridcell_index = 5
   integer, parameter, public :: subgrid_topounit_index = 4
   integer, parameter, public :: subgrid_landunit_index = 3
   integer, parameter, public :: subgrid_column_index = 2
   integer, parameter, public :: subgrid_patch_index = 1

    ! Registry keys parameters
    character(len=*), parameter, public :: hlm_fates_soil_level = 'soil_level_number'
    character(len=*), parameter, public :: hlm_fates_decomp_frac_moisture = 'decomp_frac_moisture'
    character(len=*), parameter, public :: hlm_fates_decomp_frac_temperature = 'decomp_frac_temperature'
    
    ! Registry update frequency parameters
    integer, parameter :: registry_update_init = 1       ! variable only needs to be updated during initialization
    integer, parameter :: registry_update_daily = 2      ! variable needs to be updated daily
    integer, parameter :: registry_update_timestep = 3   ! variable needs to be updated at each timestep
    integer, parameter :: registry_update_types_num = 3  ! number of update frequency types

   ! -------------------------------------------------------------------------------------
   ! These vectors are used for history output mapping
   ! CLM/ALM have limited support for multi-dimensional history output arrays.
   ! FATES structure and composition is multi-dimensional, so we end up "multi-plexing"
   ! multiple dimensions into one dimension.  These new dimensions need definitions,
   ! mapping to component dimensions, and definitions for those component dimensions as
   ! well.
   ! -------------------------------------------------------------------------------------
   
   real(r8), public, allocatable :: fates_hdim_levcoage(:)         ! cohort age class lower bound dimension
   integer , public, allocatable :: fates_hdim_pfmap_levcapf(:)    ! map of pfts into cohort age class x pft dimension
   integer , public, allocatable :: fates_hdim_camap_levcapf(:)    ! map of cohort age class into cohort age x pft dimension
   real(r8), public, allocatable :: fates_hdim_levsclass(:)        ! plant size class lower bound dimension
   integer , public, allocatable :: fates_hdim_pfmap_levscpf(:)    ! map of pfts into size-class x pft dimension
   integer , public, allocatable :: fates_hdim_scmap_levscpf(:)    ! map of size-class into size-class x pft dimension
   real(r8), public, allocatable :: fates_hdim_levage(:)           ! patch age lower bound dimension
   real(r8), public, allocatable :: fates_hdim_levheight(:)        ! height lower bound dimension
   integer , public, allocatable :: fates_hdim_levpft(:)           ! plant pft dimension
   integer , public, allocatable :: fates_hdim_levlanduse(:)       ! land use label dimension
   integer , public, allocatable :: fates_hdim_levfuel(:)          ! fire fuel size class (fsc) dimension
   integer , public, allocatable :: fates_hdim_levcwdsc(:)         ! cwd class dimension
   integer , public, allocatable :: fates_hdim_levcan(:)           ! canopy-layer dimension 
   real(r8), public, allocatable :: fates_hdim_levleaf(:)          ! leaf-layer dimension, integrated VAI [m2/m2]
   integer , public, allocatable :: fates_hdim_levelem(:)              ! element dimension
   integer , public, allocatable :: fates_hdim_canmap_levcnlf(:)   ! canopy-layer map into the canopy-layer x leaf-layer dim
   integer , public, allocatable :: fates_hdim_lfmap_levcnlf(:)    ! leaf-layer map into the can-layer x leaf-layer dimension
   integer , public, allocatable :: fates_hdim_canmap_levcnlfpf(:) ! can-layer map into the can-layer x pft x leaf-layer dim
   integer , public, allocatable :: fates_hdim_lfmap_levcnlfpf(:)  ! leaf-layer map into the can-layer x pft x leaf-layer dim
   integer , public, allocatable :: fates_hdim_pftmap_levcnlfpf(:) ! pft map into the canopy-layer x pft x leaf-layer dim
   integer , public, allocatable :: fates_hdim_scmap_levscag(:)    ! map of size-class into size-class x patch age dimension
   integer , public, allocatable :: fates_hdim_agmap_levscag(:)    ! map of patch-age into size-class x patch age dimension
   integer , public, allocatable :: fates_hdim_scmap_levscagpft(:)     ! map of size-class into size-class x patch age x pft dimension
   integer , public, allocatable :: fates_hdim_agmap_levscagpft(:)     ! map of patch-age into size-class x patch age x pft dimension
   integer , public, allocatable :: fates_hdim_pftmap_levscagpft(:)    ! map of pft into size-class x patch age x pft dimension
   integer , public, allocatable :: fates_hdim_agmap_levagepft(:)      ! map of patch-age into patch age x pft dimension
   integer , public, allocatable :: fates_hdim_pftmap_levagepft(:)     ! map of pft into patch age x pft dimension
   integer , public, allocatable :: fates_hdim_agmap_levagefuel(:)     ! map of patch-age into patch age x fsc dimension
   integer , public, allocatable :: fates_hdim_fscmap_levagefuel(:)    ! map of fuel size-class into patch age x fsc dimension
   
   integer , public, allocatable :: fates_hdim_elmap_levelpft(:)       ! map of elements in the element x pft dimension
   integer , public, allocatable :: fates_hdim_elmap_levelcwd(:)       ! map of elements in the element x cwd dimension
   integer , public, allocatable :: fates_hdim_elmap_levelage(:)       ! map of elements in the element x age dimension
   integer , public, allocatable :: fates_hdim_pftmap_levelpft(:)       ! map of pfts in the element x pft dimension
   integer , public, allocatable :: fates_hdim_cwdmap_levelcwd(:)       ! map of cwds in the element x cwd dimension
   integer , public, allocatable :: fates_hdim_agemap_levelage(:)       ! map of ages in the element x age dimension
   
   integer , public, allocatable :: fates_hdim_pftmap_levcdpf(:)   ! map of pfts into size x crowndamage x pft dimension
   integer , public, allocatable :: fates_hdim_cdmap_levcdpf(:)    ! map of crowndamage into size x crowndamage x pft
   integer , public, allocatable :: fates_hdim_scmap_levcdpf(:)    ! map of size into size x crowndamage x pft
   integer , public, allocatable :: fates_hdim_cdmap_levcdsc(:)    ! map of crowndamage into size x crowndamage
   integer , public, allocatable :: fates_hdim_scmap_levcdsc(:)    ! map of size into size x crowndamage
   integer , public, allocatable :: fates_hdim_levdamage(:)        ! plant damage class lower bound dimension   
   
   ! ------------------------------------------------------------------------------------
   !                              DYNAMIC BOUNDARY CONDITIONS
   ! ------------------------------------------------------------------------------------


   ! -------------------------------------------------------------------------------------
   ! Scalar Timing Variables
   ! It is assumed that all of the sites on a given machine will be synchronous.
   ! It is also assumed that the HLM will control time.
   ! -------------------------------------------------------------------------------------
   integer, public  :: hlm_current_year    ! Current year
   integer, public  :: hlm_current_month   ! month of year
   integer, public  :: hlm_current_day     ! day of month
   integer, public  :: hlm_current_tod     ! time of day (seconds past 0Z)
   integer, public  :: hlm_current_date    ! time of day (seconds past 0Z)
   integer, public  :: hlm_reference_date  ! YYYYMMDD
   real(r8), public :: hlm_model_day       ! elapsed days between current date and ref
   integer, public  :: hlm_day_of_year     ! The integer day of the year
   integer, public  :: hlm_days_per_year   ! The HLM controls time, some HLMs may 
                                                      ! include a leap
   real(r8), public :: hlm_freq_day        ! fraction of year for daily time-step 
                                                      ! (1/days_per_year_, this is a frequency
   

   ! -------------------------------------------------------------------------------------
   !
   ! Constant parameters that are dictated by the fates parameter file
   !
   ! -------------------------------------------------------------------------------------

   integer, public :: numpft           ! The total number of PFTs defined in the simulation
   integer, public :: nlevsclass       ! The total number of cohort size class bins output to history
   integer, public :: nlevage          ! The total number of patch age bins output to history
   integer, public :: nlevheight       ! The total number of height bins output to history
   integer, public :: nlevcoage        ! The total number of cohort age bins output to history 
   integer, public :: nleafage         ! The total number of leaf age classes
   integer, public :: nlevdamage       ! The total number of damage classes
   
   ! -------------------------------------------------------------------------------------
   ! Structured Boundary Conditions (SITE/PATCH SCALE)
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! Naming conventions:   _si  means site dimensions (scalar in that case)
   !                       _pa  means patch dimensions
   !                       _rb  means radiation band
   !                       _sl  means soil layer
   !                       _sisl means site x soil layer
   ! ------------------------------------------------------------------------------------

   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches


      ! Soil layer structure

      integer              :: nlevsoil           ! the number of soil layers in this column
      integer              :: nlevdecomp         ! the number of soil layers in the column
                                                 ! that are biogeochemically active
      real(r8),allocatable :: zi_sisl(:)         ! interface level below a "z" level (m)
                                                 ! this contains a zero index for surface.
      real(r8),allocatable :: dz_sisl(:)         ! layer thickness (m)
      real(r8),allocatable :: z_sisl(:)          ! layer depth (m)

      ! Decomposition Layer Structure
      real(r8), allocatable :: dz_decomp_sisl(:) ! This should match dz_sisl(), unless
                                                 ! only one layer is chosen, in that
                                                 ! case, it has its own depth, which
                                                 ! has traditionally been 1 meter

      integer,allocatable  :: decomp_id(:)       ! The decomposition layer index that each
                                                 ! soil layer maps to. This will either
                                                 ! be equivalent (ie integer ascending)
                                                 ! Or, all will be 1.

      ! Decomposition fractions
      real(r8),allocatable :: w_scalar_sisl(:)   ! fraction by which decomposition is limited by moisture availability
      real(r8),allocatable :: t_scalar_sisl(:)   ! fraction by which decomposition is limited by temperature
      
      ! Fire Model

      ! 24-hour lightning or ignitions [#/km2/day]
      real(r8),allocatable :: lightning24(:)

      ! Population density [#/km2]
      real(r8),allocatable :: pop_density(:)
      
      ! Average precipitation over the last 24 hours [mm/s]
      real(r8), allocatable :: precip24_pa(:)
      
      ! Average relative humidity over past 24 hours [-]
      real(r8), allocatable :: relhumid24_pa(:)

      ! Patch 24-hour running mean of wind (m/s ?)
      real(r8), allocatable :: wind24_pa(:)

      ! Radiation variables for calculating sun/shade fractions
      ! ---------------------------------------------------------------------------------

      ! Downwelling direct beam radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solad_parb(:,:)  

      ! Downwelling diffuse (I-ndirect) radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solai_parb(:,:)

      
      ! Nutrient input fluxes (these are integrated fluxes over the day, most
      !                        likely calculated over shorter dynamics steps,
      !                        and then incremented until the end of the day)
      !
      ! Note 1: If these are indexed by COHORT, they don't also need to be indexed
      !         by decomposition layer. So it is allocated with 2nd dim=1.
      ! Note 2: Has it's own zero'ing call
      real(r8), pointer :: plant_nh4_uptake_flux(:,:) ! Ammonium uptake flux for
                                                      ! each competitor [gN/m2/day]

      real(r8), pointer :: plant_no3_uptake_flux(:,:) ! Nitrate uptake flux for
                                                      ! each competitor [gN/m2/day]
      
      real(r8), pointer :: plant_p_uptake_flux(:,:)   ! Phosphorus input flux for
                                                      ! each competitor [gP/m2/day]


      ! Photosynthesis variables
      ! ---------------------------------------------------------------------------------

      ! Patch level filter flag for photosynthesis calculations
      ! has a short memory, flags:
      ! 1 = patch has not been called
      ! 2 = patch is currently marked for photosynthesis
      ! 3 = patch has been called for photosynthesis at least once
      integer, allocatable  :: filter_photo_pa(:)

      ! atmospheric pressure (Pa)
      real(r8)              :: forc_pbot             

      ! daylength scaling factor (0-1)
      real(r8), allocatable :: dayl_factor_pa(:)
      
      ! saturation vapor pressure at t_veg (Pa)
      real(r8), allocatable :: esat_tv_pa(:)

      ! vapor pressure of canopy air (Pa)
      real(r8), allocatable :: eair_pa(:)

      ! Atmospheric O2 partial pressure (Pa)
      real(r8), allocatable :: oair_pa(:)

      ! Atmospheric CO2 partial pressure (Pa)
      real(r8), allocatable :: cair_pa(:)

      ! boundary layer resistance (s/m)
      real(r8), allocatable :: rb_pa(:)

      ! vegetation temperature (Kelvin)
      real(r8), allocatable :: t_veg_pa(:)
             
      ! air temperature at agcm reference height (kelvin)
      real(r8), allocatable :: tgcm_pa(:)

      ! soil temperature (Kelvin)
      real(r8), allocatable :: t_soisno_sl(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------

      ! Cosine of the zenith angle (0-1) - site level
      real(r8) :: coszen
      
      ! fraction of canopy that is covered in snow
      real(r8), allocatable :: fcansno_pa(:)
       
      ! Abledo of the ground for direct radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dir_rb(:)

      ! Albedo of the ground for diffuse radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dif_rb(:)
      
      ! LitterFlux Boundaries
      ! the index of the deepest model soil level where roots may be
      ! due to permafrost or bedrock constraints
      integer  :: max_rooting_depth_index_col

      ! BGC Accounting

      real(r8) :: tot_het_resp  ! total heterotrophic respiration  (gC/m2/s)
      real(r8) :: tot_somc      ! total soil organic matter carbon (gc/m2)
      real(r8) :: tot_litc      ! total litter carbon tracked in the HLM (gc/m2)

      ! Canopy Structure

      real(r8) :: snow_depth_si    ! Depth of snow in snowy areas of site (m)
      real(r8) :: frac_sno_eff_si  ! Fraction of ground covered by snow (0-1)

      ! Hydrology variables for BTRAN
      ! ---------------------------------------------------------------------------------

      ! Soil suction potential of layers in each site, negative, [mm]
      real(r8), allocatable :: smp_sl(:)
      
      !soil salinity of layers in each site [ppt]
      real(r8), allocatable :: salinity_sl(:)

      ! Effective porosity = porosity - vol_ic, of layers in each site [-]
      real(r8), allocatable :: eff_porosity_sl(:)

      ! volumetric soil water at saturation (porosity)
      real(r8), allocatable :: watsat_sl(:)

      ! Temperature of soil layers [K]
      real(r8), allocatable :: tempk_sl(:)

      ! Liquid volume in soil layer (m3/m3)
      real(r8), allocatable :: h2o_liqvol_sl(:)

      ! Site level filter for uptake response functions
      logical               :: filter_btran


      ! ALL HYDRO DATA STRUCTURES SHOULD NOW BE ALLOCATED ON RHIZOSPHERE LEVELS
      
      ! Plant-Hydro
      ! ---------------------------------------------------------------------------------
      
      real(r8),allocatable :: qflx_transp_pa(:)    ! Transpiration flux as dictated by the HLM's
                                                   ! canopy solver. [mm H2O/s] [+ into root]
      real(r8),allocatable :: swrad_net_pa(:)      ! Net absorbed shortwave radiation (W/m2)
      real(r8),allocatable :: lwrad_net_pa(:)      ! Net absorbed longwave radiation (W/m2)
      real(r8),allocatable :: watsat_sisl(:)       ! volumetric soil water at saturation (porosity)
      real(r8),allocatable :: watres_sisl(:)       ! volumetric residual soil water
      real(r8),allocatable :: sucsat_sisl(:)       ! minimum soil suction (mm) 
      real(r8),allocatable :: bsw_sisl(:)          ! Clapp and Hornberger "b"
      real(r8),allocatable :: hksat_sisl(:)        ! hydraulic conductivity at saturation (mm H2O /s)
      real(r8),allocatable :: h2o_liq_sisl(:)      ! Liquid water mass in each layer (kg/m2)
      real(r8) :: smpmin_si                        ! restriction for min of soil potential (mm)

      ! Land use
      ! ---------------------------------------------------------------------------------
      real(r8),allocatable :: hlm_harvest_rates(:)    ! annual harvest rate per cat from hlm for a site
      character(len=64), allocatable :: hlm_harvest_catnames(:)  ! names of hlm_harvest d1
      real(r8),allocatable :: hlm_luh_states(:)
      character(len=64),allocatable :: hlm_luh_state_names(:)
      real(r8),allocatable :: hlm_luh_transitions(:)
      character(len=64),allocatable :: hlm_luh_transition_names(:)


      integer :: hlm_harvest_units  ! what units are the harvest rates specified in? [area vs carbon]
    
      real(r8) :: site_area    ! Actual area of current site [m2], only used in carbon-based harvest

      ! Fixed biogeography mode 
      real(r8), allocatable :: pft_areafrac(:)          ! Fractional area of the FATES column occupied by each PFT

      ! Fixed biogeography mode with land use active
      real(r8), allocatable :: pft_areafrac_lu(:,:)     ! Fractional area occupied by each PFT on each land use type
      real(r8) :: baregroundfrac                        ! fractional area held as bare-ground

    
     ! Satellite Phenology (SP) input variables.  (where each patch only has one PFT)
     ! ---------------------------------------------------------------------------------
     real(r8),allocatable :: hlm_sp_tlai(:)  ! Interpolated daily total LAI (leaf area index) input from HLM per patch/pft 
     real(r8),allocatable :: hlm_sp_tsai(:)  ! Interpolated sailt total SAI (stem area index) input from HLM per patch/pft
     real(r8),allocatable :: hlm_sp_htop(:)  ! Interpolated daily canopy vegetation height    input from HLM per patch/pft
     
     contains
     
      procedure :: Initialize => InitializeBCIn

   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

      ! Sunlit canopy LAI
      real(r8),allocatable :: laisun_pa(:)
      
      ! Shaded canopy LAI
      real(r8),allocatable :: laisha_pa(:)
      
      ! Logical stating whether a soil layer can have water uptake by plants
      ! The only condition right now is that liquid water exists
      ! The name (suction) is used to indicate that soil suction should be calculated
      logical, allocatable :: active_suction_sl(:)

      ! Effective fraction of roots in each soil layer 
      real(r8), allocatable :: rootr_pasl(:,:)

      ! Integrated (vertically) transpiration wetness factor (0 to 1) 
      ! (diagnostic, should not be used by HLM)
      real(r8), allocatable :: btran_pa(:)

      ! Sunlit canopy resistance [s/m]
      real(r8), allocatable :: rssun_pa(:)

      ! Shaded canopy resistance [s/m]
      real(r8), allocatable :: rssha_pa(:)

      ! leaf photosynthesis (umol CO2 /m**2/ s)
      ! (NOT CURRENTLY USED, PLACE-HOLDER)
      !real(r8), allocatable :: psncanopy_pa(:)

      ! leaf maintenance respiration rate (umol CO2/m**2/s) 
      ! (NOT CURRENTLY USED, PLACE-HOLDER)
      !real(r8), allocatable :: lmrcanopy_pa(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------

      ! Note: We initialize and default the radiatioon balance to assume that the
      ! canopy is invisible, and the soil absorbs all radiation.
      
      ! Surface albedo (direct) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albd_parb(:,:)
      
      ! Surface albedo (diffuse) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albi_parb(:,:)                 
      
      ! Flux absorbed by canopy per unit direct flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabd_parb(:,:) 
      
      ! Flux absorbed by canopy per unit diffuse flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabi_parb(:,:)

      ! Down direct flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftdd_parb(:,:)

      ! Down diffuse flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftid_parb(:,:)
      
      ! Down diffuse flux below canopy per unit diffuse flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftii_parb(:,:)


      ! Mass fluxes to BGC from fragmentation of litter into decomposing pools
      
      real(r8), allocatable :: litt_flux_cel_c_si(:) ! cellulose carbon litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lig_c_si(:) ! lignin carbon litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_c_si(:) ! labile carbon litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_cel_n_si(:) ! cellulose nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lig_n_si(:) ! lignin nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_n_si(:) ! labile nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_cel_p_si(:) ! cellulose phosphorus litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lig_p_si(:) ! lignin phosphorus litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_p_si(:) ! labile phosphorus litter, fates->BGC g/m3/s

      
      ! MIMICS Boundary Conditions
      ! -----------------------------------------------------------------------------------
      real(r8) :: litt_flux_ligc_per_n  ! lignin carbon per total nitrogen
                                        ! in the fragmentation flux, per square meter [g/g]
      

      ! Nutrient competition boundary conditions
      ! (These are all pointer allocations, this is because the host models
      !  will point to these arrays)
      ! ---------------------------------------------------------------------------------

      integer               :: num_plant_comps ! Number of unique competitors

      real(r8), allocatable :: source_nh4(:) ! FATES generated source of ammonium to the mineralized N pool
                                             ! in the BGC model [gN/m3]
      real(r8), allocatable :: source_p(:)   ! FATES generated source of phosphorus to mineralized P
                                             ! pool in the BGC model [gP/m3]
      
      real(r8), pointer :: veg_rootc(:,:)    ! Total fine-root carbon of each competitor
                                             ! [gC/m3 of site area]  
                                             ! (maxcohort_per_site x nlevdecomp)
      real(r8), pointer :: decompmicc(:)     ! Microbial decomposer biomass [gc/m3] 
                                             ! (numpft x nledecomp_full)
      integer, pointer :: ft_index(:)        ! functional type index of each competitor
                                             ! (maxcohort_per_site)
      real(r8), pointer :: cn_scalar(:)      ! C:N scaling factor for root n uptake 
                                             ! kinetics (exact meaning differs between
                                             ! soil BGC hypotheses)
      real(r8), pointer :: cp_scalar(:)      ! C:P scaling factor for root p uptake
                                             ! kinetics (exact meaning differs between
                                             ! soil BGC hypotheses)




      ! RD Nutrient Boundary Conditions
      ! ---------------------------------------------------------------------------------

      !real(r8), pointer :: n_demand(:)       ! Nitrogen demand from each competitor
      !                                       ! for use in ELMs CTC/RD [g/m2/s] 
      !real(r8), pointer :: p_demand(:)       ! Phosophorus demand from each competitor
      !                                       ! for use in ELMs CTC/RD [g/m2/s] 




      
      ! CH4 Boundary Conditions
      ! -----------------------------------------------------------------------------------
      real(r8), pointer :: annavg_agnpp_pa(:)    ! annual average patch npp above ground (gC/m2/s)
      real(r8), pointer :: annavg_bgnpp_pa(:)    ! annual average patch npp below ground (gC/m2/s)
      real(r8), pointer :: annsum_npp_pa(:)      ! annual sum patch npp (gC/m2/yr)
      real(r8), pointer :: frootc_pa(:)          ! Carbon in fine roots (gC/m2)
      real(r8), pointer :: root_resp(:)          ! (gC/m2/s) root respiration (fine root MR + total root GR)
      real(r8), pointer :: rootfr_pa(:,:)        ! Rooting fraction with depth
      real(r8), pointer :: woody_frac_aere_pa(:) ! Woody plant fraction (by crown area) of all plants
                                                 ! used for calculating patch-level aerenchyma porosity
      
      real(r8)          :: ema_npp               ! site-level NPP smoothed over time, see PrepCH4BCs()
                                                 ! used for N fixation in ELM/CLM right now
      ! Canopy Structure

      real(r8), allocatable :: elai_pa(:)  ! exposed leaf area index
      real(r8), allocatable :: esai_pa(:)  ! exposed stem area index
      real(r8), allocatable :: tlai_pa(:)  ! total leaf area index
      real(r8), allocatable :: tsai_pa(:)  ! total stem area index
      real(r8), allocatable :: htop_pa(:)  ! top of the canopy [m]
      real(r8), allocatable :: hbot_pa(:)  ! bottom of canopy? [m]

      real(r8), allocatable :: z0m_pa(:)   ! roughness length [m]
      real(r8), allocatable :: displa_pa(:) ! displacement height [m]
      real(r8), allocatable :: dleaf_pa(:)  ! leaf characteristic dimension/width/diameter [m]

      real(r8), allocatable :: canopy_fraction_pa(:) ! Area fraction of each patch in the site
                                                     ! Use most likely for weighting
                                                     ! This is currently the projected canopy
                                                     ! area of each patch [0-1]

      real(r8), allocatable :: frac_veg_nosno_alb_pa(:) ! This is not really a fraction
                                                        ! this is actually binary based on if any
                                                        ! vegetation in the patch is exposed.
                                                        ! [0,1]

     integer, allocatable :: nocomp_pft_label_pa(:) ! in nocomp and SP mode, each patch has a PFT identity. 

      ! FATES Hydraulics


      
      real(r8) :: plant_stored_h2o_si         ! stored water in LIVE+DEAD vegetation (kg/m2 H2O)
                                              ! Assuming density of 1Mg/m3 ~= mm/m2 H2O
                                              ! This must be set and transfered prior to clm_drv()
                                              ! following the calls to ed_update_site()
                                              ! ed_update_site() is called during both the restart
                                              ! and coldstart process

      real(r8),allocatable :: qflx_soil2root_sisl(:)   ! Water flux from soil into root by site and soil layer
                                                       ! [mm H2O/s] [+ into root]
      
      real(r8),allocatable :: qflx_ro_sisl(:)          ! Water flux runoff generated by
                                                       ! root to soil flux super-saturating the soils
                                                       ! This does seem unlikely, but we need accomodate
                                                       ! small fluxes for various reasons
                                                       ! [mm H2O/s]

      ! FATES LULCC
      real(r8) :: hrv_deadstemc_to_prod10c   ! Harvested C flux to 10-yr wood product pool [Site-Level, gC m-2 s-1]
      real(r8) :: hrv_deadstemc_to_prod100c  ! Harvested C flux to 100-yr wood product pool [Site-Level, gC m-2 s-1]
      real(r8) :: gpp_site  ! Site level GPP, for NBP diagnosis in HLM [Site-Level, gC m-2 s-1]
      real(r8) :: ar_site   ! Site level Autotrophic Resp, for NBP diagnosis in HLM [Site-Level, gC m-2 s-1]

      ! direct carbon loss to atm pathways
      real(r8) :: grazing_closs_to_atm_si    ! Loss of carbon to atmosphere via grazing [Site-Level, gC m-2 s-1]
      real(r8) :: fire_closs_to_atm_si       ! Loss of carbon to atmosphere via burning (includes burning from land use change) [Site-Level, gC m-2 s-1]

      ! summary carbon stock variables
      real(r8) :: veg_c_si                   ! Total vegetation carbon [Site-Level, gC m-2]
      real(r8) :: litter_cwd_c_si            ! Total litter plus CWD carbon [Site-Level, gC m-2]
      real(r8) :: seed_c_si                  ! Total seed carbon [Site-Level, gC m-2]

   end type bc_out_type


   ! This type holds parameter constants
   ! These parameter constants only need to specified once, and never modified again.
   ! After re-factoring this module to split the procedures from the data-types
   ! we can then set the datatypes as protected.

   type, public :: bc_pconst_type

       ! Nutrient competition boundary conditions for ECA hypothesis
       ! Note, these "could" be stored globaly for each machine, saving them on
       ! each column is inefficient. Each of these are dimensioned by PFT.
       
       integer           :: max_plant_comps

       real(r8), pointer :: vmax_nh4(:)
       real(r8), pointer :: vmax_no3(:)
       real(r8), pointer :: vmax_p(:)
       
       real(r8), pointer :: eca_km_nh4(:)
       real(r8), pointer :: eca_km_no3(:)
       real(r8), pointer :: eca_km_p(:)
       
       real(r8), pointer :: eca_km_ptase(:)     
       real(r8), pointer :: eca_vmax_ptase(:)
       real(r8), pointer :: eca_alpha_ptase(:)
       real(r8), pointer :: eca_lambda_ptase(:)
       real(r8)          :: eca_plant_escalar

       integer, pointer  :: j_uptake(:)         ! Mapping between decomposition
                                                ! layers and the uptake layers
                                                ! in FATES (is either incrementally
                                                ! increasing, or all 1s)

   end type bc_pconst_type

   ! Base type to be extended for the API registry
   type, public :: fates_interface_registry_base_type

    ! Container array of interface variables indexed by key
    type(fates_interface_variable_type), allocatable :: hlm_vars(:) 
    type(fates_interface_variable_type), allocatable :: fates_vars(:) 
    
    ! Array of keys associated with the interface variables
    character(len=48), allocatable :: key(:)

    ! Variable regsitry metadata
    integer :: num_api_vars                        ! number of variables in the registry
    integer :: num_api_vars_update_init           ! number of variables that update only at initialization
    integer :: num_api_vars_update_daily          ! number of variables that update daily

    ! Array of update frequency values for each variable index
    integer, allocatable :: update_frequency(:)

    ! Arrays that hold the indices of variables based on update frequency
    integer, allocatable :: filter_init(:)       ! index of variables that update only at initialization
    integer, allocatable :: filter_daily(:)      ! index of variables that update daily
    ! integer, allocatable :: filter_timestep(:)   ! index of variables that update at each timestep

    ! Subgrid index data
    integer, private :: gidx
    integer, private :: tidx
    integer, private :: lidx
    integer, private :: cidx
    integer, private :: sidx
    integer, private :: hpidx
    integer, private :: fpidx

    contains

      procedure :: InitializeInterfaceRegistry
      procedure :: InitializeInterfaceVariables
      procedure :: SetSubgridIndices
      procedure :: Update => UpdateInterfaceVariables

      generic :: Register => RegisterInterfaceVariables_0d, & 
                             RegisterInterfaceVariables_1d, &
                             RegisterInterfaceVariables_2d
      procedure, private :: RegisterInterfaceVariables_0d
      procedure, private :: RegisterInterfaceVariables_1d
      procedure, private :: RegisterInterfaceVariables_2d

      procedure, private :: DefineInterfaceRegistry
      procedure, private :: DefineInterfaceVariable
      procedure, private :: SetFilterMapArrays
      procedure, private :: GetRegistryIndex
      procedure, private :: GetRegistryKey

   end type fates_interface_registry_base_type

   public :: ZeroBCOutCarbonFluxes
   
 contains
  
       
  ! ======================================================================================
 
  subroutine InitializeBCIn(this)
    
    class(bc_in_type), intent(inout) :: this
    
    ! Allocate the boundary condition variables
    allocate(this%w_scalar_sisl(this%nlevsoil))
    allocate(this%t_scalar_sisl(this%nlevsoil))
    
    ! Unset the values
    this%w_scalar_sisl = nan
    this%t_scalar_sisl = nan
    
  end subroutine InitializeBCIn
 
  ! ======================================================================================

  subroutine ZeroBCOutCarbonFluxes(bc_out)

    ! !ARGUMENTS
    type(bc_out_type), intent(inout)   :: bc_out

    bc_out%grazing_closs_to_atm_si = nan    ! set via site_mass%burn_flux
    bc_out%fire_closs_to_atm_si    = nan    ! set via site_mass%herbivory_flux_out
    bc_out%gpp_site                = 0._r8
    bc_out%ar_site                 = 0._r8

  end subroutine ZeroBCOutCarbonFluxes

  ! ======================================================================================

  subroutine InitializeInterfaceRegistry(this)

    ! This initializes the interface registry

    class(fates_interface_registry_base_type), intent(inout) :: this

    logical :: initialize

    ! Initialize registry counters
    this%num_api_vars = 0
    this%num_api_vars_update_init = 0
    this%num_api_vars_update_daily = 0
    
    ! First count up the keys defined in the registry and the registry counters
    call this%DefineInterfaceRegistry(initialize=.false.)

    ! Allocate the registry variables arrays
    allocate(this%fates_vars(this%num_api_vars))
    allocate(this%hlm_vars(this%num_api_vars))
    allocate(this%key(this%num_api_vars))
    allocate(this%update_frequency(this%num_api_vars))
    
    ! Unset the allocatables not including the interface variables
    this%update_frequency(:) = fates_unset_int
    this%filter_init(:) = fates_unset_int
    this%filter_daily(:) = fates_unset_int

    ! Now initialize the registry keys
    call this%DefineInterfaceRegistry(initialize=.true.)

    ! Allocate the index filter maps
    allocate(this%filter_init(this%num_api_vars_update_init))
    allocate(this%filter_daily(this%num_api_vars_update_daily))
    
    ! Set filter map arrays
    call this%SetFilterMapArrays()

  end subroutine InitializeInterfaceRegistry

  ! ======================================================================================
  
  subroutine DefineInterfaceRegistry(this, initialize)

    ! This procedure defines the list of common names to be associated with FATES and HLM
    ! variables.

    class(fates_interface_registry_base_type), intent(inout) :: this

    logical, intent(in) :: initialize  ! false = count up the keys in the registry

    integer :: index   ! Index to be incremented for each call to DefineInterfaceVariable()
    
    ! Initialize the index
    index = 0

    ! Define the interface registry names and indices
    ! Variables that only need to be updated during initialization, such as dimensions
    call this%DefineInterfaceVariable(key=hlm_fates_soil_level, initialize=initialize, index=index, &
                                      update_frequency=registry_update_init)
                                      
    
    ! Variables that need to be updated daily
    call this%DefineInterfaceVariable(key=hlm_fates_decomp_frac_moisture, initialize=initialize, index=index)
    call this%DefineInterfaceVariable(key=hlm_fates_decomp_frac_temperature, initialize=initialize, index=index)

  end subroutine DefineInterfaceRegistry

  ! ======================================================================================

  subroutine DefineInterfaceVariable(this, key, initialize, index, update_frequency)

    class(fates_interface_registry_base_type), intent(inout) :: this

    character(len=*), intent(in)  :: key
    logical, intent(in)           :: initialize
    integer, intent(inout)        :: index
    integer, intent(in), optional :: update_frequency

    ! Local variables
    integer :: index_type
    integer :: update_frequency_local
    
    ! Increment the index
    index = index + 1  

    ! If not initializing, increment the registry count variables, otherwise initialize the variable at the correct index
    if (initialize) then 

      ! Initialize the variable
      if (present(update_frequency)) then
        select case (update_frequency)
        case (registry_update_init)
          update_frequency_local = registry_update_init
        case (registry_update_daily)
          update_frequency_local = registry_update_daily
        case default
          write(fates_log(),*) 'ERROR: Unrecognized update frequency in DefineInterfaceVariable(): ', update_frequency
          call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      else
        ! Default to daily update frequency
        update_frequency_local = registry_update_daily
      end if

      ! Set the key for each index
      this%key(index) = key
      
      ! Set the update frequency array values
      this%update_frequency(index) = update_frequency_local
      
      ! Initialize the interface variables and pass the key and update frequency to each for metadata
      call this%hlm_vars(index)%Initialize(key, update_frequency_local)
      call this%fates_vars(index)%Initialize(key, update_frequency_local)

    ! Not initializing, just counting the variables
    else 

      ! Increment the total API count
      this%num_api_vars = this%num_api_vars + 1

      ! Increment the count for the update frequency counts, defaulting to daily if not specified
      if (present(update_frequency)) then
        select case (update_frequency)
        case (registry_update_init)
          this%num_api_vars_update_init = this%num_api_vars_update_init + 1
        case (registry_update_daily)
          this%num_api_vars_update_daily = this%num_api_vars_update_daily + 1
        case default
          write(fates_log(),*) 'ERROR: Unrecognized update frequency in DefineInterfaceVariable(): ', update_frequency
          call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
      else
        ! Default to daily update frequency
        this%num_api_vars_update_daily = this%num_api_vars_update_daily + 1
      end if
      
    end if

  end subroutine DefineInterfaceVariable

  ! ======================================================================================
  
  subroutine SetSubgridIndices(this, gridcell, topounit, landunit, column, hlmpatch)
    
    class(fates_interface_registry_base_type), intent(inout) :: this
    integer, intent(in) :: gridcell
    integer, intent(in) :: topounit
    integer, intent(in) :: landunit
    integer, intent(in) :: column
    integer, intent(in) :: hlmpatch
    
    this%gidx = gridcell
    this%tidx = topounit
    this%lidx = landunit
    this%cidx = column
    this%hpidx = hlmpatch
    
  end subroutine SetSubgridIndices
  
  ! ======================================================================================
  
  integer function GetColumnIndex(this) result(cidx)
  
    class(fates_interface_registry_base_type), intent(inout) :: this
    
    cidx = this%cidx
    
  end function GetColumnIndex
  
  ! ======================================================================================
  
  subroutine SetFilterMapArrays(this)

    class(fates_interface_registry_base_type), intent(inout) :: this

    integer :: index
    integer :: count_init
    integer :: count_daily

    ! Initialize counters
    count_init = 0
    count_daily = 0
    
    ! Iterate over all registered variables and populate the filter maps accordingly
    do index = 1, this%num_api_vars
      if (this%update_frequency(index) == registry_update_init) then
        count_init = count_init + 1
        this%filter_init(count_init) = index
      else if (this%update_frequency(index) == registry_update_daily) then
        count_daily = count_daily + 1
        this%filter_daily(count_daily) = index
      else
        write(fates_log(),*) 'ERROR: Unrecognized update frequency in SetFilterMapArrays(): ', this%update_frequency(index)
        call endrun(msg=errMsg(__FILE__, __LINE__))
      end if
    end do
    
    ! Check that the counts match the expected sizes
    if (count_init /= this%num_api_vars_update_init .or. &
        count_daily /= this%num_api_vars_update_daily) then
          
      write(fates_log(),*) 'ERROR: Mismatch in initialization counts in SetFilterMapArrays(): '
      write(fates_log(),*) '  count_init = ', count_init, ' expected = ', this%num_api_vars_update_init
      write(fates_log(),*) '  count_daily = ', count_daily, ' expected = ', this%num_api_vars_update_daily  
      call endrun(msg=errMsg(__FILE__, __LINE__))

    end if

  end subroutine SetFilterMapArrays
  
  ! ======================================================================================

  subroutine RegisterInterfaceVariables_0d(this, key, data, hlm_flag)

    ! This procedure is called by the to associate a data variable
    ! with a particular registry key

    ! Arguments
    class(fates_interface_registry_base_type), intent(inout) :: this
    class(*), target, intent(in)  :: data          ! data to be associated with key
    character(len=*), intent(in)  :: key           ! variable registry key 
    logical, intent(in)           :: hlm_flag      ! Is the variable being register from the HLM?
    
    ! Get index from registry key and associate the given data pointer
    if (hlm_flag) then
      call this%hlm_vars(this%GetRegistryIndex(key))%Register(data, active=.true.)
    else
      call this%fates_vars(this%GetRegistryIndex(key))%Register(data, active=.true.)
    end if


  end subroutine RegisterInterfaceVariables_0d

  ! ======================================================================================

  subroutine RegisterInterfaceVariables_1d(this, key, data, hlm_flag)

    ! This procedure is called by the to associate a data variable
    ! with a particular registry key

    class(fates_interface_registry_base_type), intent(inout) :: this
    class(*), target, intent(in)  :: data(:)    ! data to be associated with key
    character(len=*), intent(in)  :: key        ! variable registry key 
    logical, intent(in)           :: hlm_flag   ! Is the variable being register from the HLM?

    ! Get index from registry key and associate the given data pointer
    if (hlm_flag) then
      call this%hlm_vars(this%GetRegistryIndex(key))%Register(data(:), active=.true.)
    else
      call this%fates_vars(this%GetRegistryIndex(key))%Register(data(:), active=.true.)
    end if

  end subroutine RegisterInterfaceVariables_1d

  ! ======================================================================================

  subroutine RegisterInterfaceVariables_2d(this, key, data, hlm_flag)

    ! This procedure is called by the to associate a data variable
    ! with a particular registry key

    class(fates_interface_registry_base_type), intent(inout) :: this
    class(*), target, intent(in)  :: data(:,:)  ! data to be associated with key
    character(len=*), intent(in)  :: key        ! variable registry key 
    logical, intent(in)           :: hlm_flag   ! Is the variable being register from the HLM?

    ! Get index from registry key and associate the given data pointer
    if (hlm_flag) then
      call this%hlm_vars(this%GetRegistryIndex(key))%Register(data(:,:), active=.true.)
    else
      call this%fates_vars(this%GetRegistryIndex(key))%Register(data(:,:), active=.true.)
    end if

  end subroutine RegisterInterfaceVariables_2d

  ! ======================================================================================
  
  subroutine InitializeInterfaceVariables(this)

    ! Arguments
    class(fates_interface_registry_base_type), intent(inout) :: this  ! registry being initialized

    ! Locals
    integer :: i  ! initialization iterator
    integer :: j  ! variable index
    
    ! Update the boundary conditions necessary during initialization only
    do i = 1, this%num_api_vars_update_init
      
      ! Get the variable index from the init filter
      j = this%filter_init(i)
      
      ! Update the variables
      call this%fates_vars(j)%Update(hlm_vars(j))

    end do
    
  end subroutine
  
  ! ======================================================================================

  subroutine UpdateInterfaceVariables(this)

    class(fates_interface_registry_base_type), intent(inout) :: this

    integer :: i

    ! Iterate over all registered variables
    do i = 1, this%num_api_vars
      this%fates_vars(i)%Update(hlm_vars(i))
    end do

  end subroutine UpdateInterfaceVariables

  ! ======================================================================================

  integer function GetRegistryIndex(this, key) result(index)

    ! This procedure returns the index associated with the key provided

    class(fates_interface_registry_base_type), intent(in) :: this

    character(len=*), intent(in) :: key    ! variable registry key to search

    integer :: ivar  ! Iterator

    ! Iterate over the registry until the associated key is found
    do ivar = 1, this%num_api_vars
      if (this%key(ivar) == key) then
        index = ivar
        return
      end if
    end do

  end function GetRegistryIndex

  ! ======================================================================================

  function GetRegistryKey(this, index) result(key)

    ! This procedure returns the index associated with the key provided

    class(fates_interface_registry_base_type), intent(in) :: this

    integer, intent(in) :: index    ! variable registry index
    character(len=:), allocatable :: key

    key = this%key(index)

  end function GetRegistryKey

  ! ======================================================================================
  
end module FatesInterfaceTypesMod
