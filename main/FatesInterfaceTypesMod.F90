module FatesInterfaceTypesMod
  
  use FatesConstantsMod   , only : r8 => fates_r8
  use FatesConstantsMod   , only : itrue,ifalse
  use FatesGlobals        , only : fates_global_verbose
  use FatesGlobals        , only : fates_log
  use FatesGlobals        , only : endrun => fates_endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
  use EDTypesMod          , only : ed_site_type
  
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


   integer, public :: hlm_numlevgrnd   ! Number of ground layers
                                                  ! NOTE! SOIL LAYERS ARE NOT A GLOBAL, THEY 
                                                  ! ARE VARIABLE BY SITE

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
   
   integer, public :: hlm_max_patch_per_site ! The HLM needs to exchange some patch
                                                        ! level quantities with FATES
                                                        ! FATES does not dictate those allocations
                                                        ! since it happens pretty early in
                                                        ! the model initialization sequence.
                                                        ! So we want to at least query it,
                                                        ! compare it to our maxpatchpersite,
                                                        ! and gracefully halt if we are over-allocating

   integer, public :: hlm_parteh_mode   ! This flag signals which Plant Allocation and Reactive
                                                   ! Transport (exensible) Hypothesis (PARTEH) to use


   integer, public :: hlm_use_vertsoilc ! This flag signals whether or not the 
                                                   ! host model is using vertically discretized
                                                   ! soil carbon
                                                   ! 1 = TRUE,  0 = FALSE
   
   integer, public :: hlm_use_spitfire  ! This flag signals whether or not to use SPITFIRE
                                                   ! 1 = TRUE, 0 = FALSE


   integer, public :: hlm_use_logging       ! This flag signals whether or not to use
                                                       ! the logging module

   integer, public :: hlm_use_planthydro    ! This flag signals whether or not to use
                                                       ! plant hydraulics (bchristo/xu methods)
                                                       ! 1 = TRUE, 0 = FALSE
                                                       ! THIS IS CURRENTLY NOT SUPPORTED 

   integer, public :: hlm_use_cohort_age_tracking ! This flag signals whether or not to use
                                                             ! cohort age tracking. 1 = TRUE, 0 = FALSE

   integer, public :: hlm_use_ed_st3        ! This flag signals whether or not to use
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
   integer , public, allocatable :: fates_hdim_levfuel(:)          ! fire fuel class dimension
   integer , public, allocatable :: fates_hdim_levcwdsc(:)         ! cwd class dimension
   integer , public, allocatable :: fates_hdim_levcan(:)           ! canopy-layer dimension 
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
   
   integer , public, allocatable :: fates_hdim_elmap_levelpft(:)       ! map of elements in the element x pft dimension
   integer , public, allocatable :: fates_hdim_elmap_levelcwd(:)       ! map of elements in the element x cwd dimension
   integer , public, allocatable :: fates_hdim_elmap_levelage(:)       ! map of elements in the element x age dimension
   integer , public, allocatable :: fates_hdim_pftmap_levelpft(:)       ! map of pfts in the element x pft dimension
   integer , public, allocatable :: fates_hdim_cwdmap_levelcwd(:)       ! map of cwds in the element x cwd dimension
   integer , public, allocatable :: fates_hdim_agemap_levelage(:)       ! map of ages in the element x age dimension

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
      

      ! Vegetation Dynamics
      ! ---------------------------------------------------------------------------------

      ! Patch 24 hour vegetation temperature [K]
      real(r8),allocatable :: t_veg24_pa(:)  
      
      ! Fire Model

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
      
      ! Filter for vegetation patches with a positive zenith angle (daylight)
      logical, allocatable :: filter_vegzen_pa(:)

      ! Cosine of the zenith angle (0-1), by patch
      ! Note RGK: It does not seem like the code would currently generate
      !           different zenith angles for different patches (nor should it)
      !           I am leaving it at this scale for simplicity.  Patches should
      !           have no spacially variable information
      real(r8), allocatable :: coszen_pa(:)
      
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

      ! Temperature of ground layers [K]
      real(r8), allocatable :: tempk_sl(:)

      ! Liquid volume in ground layer (m3/m3)
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
      
   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

      ! Sunlit canopy LAI
      real(r8),allocatable :: laisun_pa(:)
      
      ! Shaded canopy LAI
      real(r8),allocatable :: laisha_pa(:)
      
      ! Logical stating whether a ground layer can have water uptake by plants
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
      real(r8), allocatable :: litt_flux_lig_c_si(:) ! lignan carbon litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_c_si(:) ! labile carbon litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_cel_n_si(:) ! cellulose nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lig_n_si(:) ! lignan nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_n_si(:) ! labile nitrogen litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_cel_p_si(:) ! cellulose phosphorus litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lig_p_si(:) ! lignan phosphorus litter, fates->BGC g/m3/s
      real(r8), allocatable :: litt_flux_lab_p_si(:) ! labile phosphorus litter, fates->BGC g/m3/s
      
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


   end type bc_out_type


   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), pointer :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)


   end type fates_interface_type


 contains
   
   ! ====================================================================================
   
   
    
  end module FatesInterfaceTypesMod
