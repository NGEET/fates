module FatesSiteMod 

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : fates_unset_int
  use FatesConstantsMod,           only : ifalse
  use FatesConstantsMod,           only : itrue
  use FatesConstantsMod,           only : numWaterMem, num_vegtemp_mem
  use EDParamsMod,                 only : nclmax
  use FatesGlobals,                only : fates_log
  use FatesPatchMod,               only : fates_patch_type
  use FatesResourcesManagementMod, only : fates_resources_management_type
  use FatesMassBalTypeMod,         only : site_massbal_type
  
  use FatesConstantsMod,           only : N_DIST_TYPES
  use FatesHydraulicsMemMod,       only : ed_site_hydr_type
  use PRTGenericMod,               only : num_elements
  use FatesConstantsMod,           only : maxpft
  use FatesLitterMod,              only : ncwd
  use FatesInterfaceTypesMod,      only : nlevage, nlevdamage, nlevsclass, numpft
  use FatesInterfaceTypesMod,      only : hlm_use_nocomp, hlm_use_fixed_biogeog
  use FatesInterfaceTypesMod,      only : hlm_use_tree_damage

  use shr_infnan_mod,              only : nan => shr_infnan_nan, assignment(=)

  implicit none
  private

  type, public :: fates_site_type
     
    ! POINTERS  
    type(fates_patch_type),    pointer :: oldest_patch => null()   ! pointer to oldest patch on site  
    type(fates_patch_type),    pointer :: youngest_patch => null() ! pointer to youngest patch on site
    type(site_massbal_type),   pointer :: mass_balance(:)          ! mass balance (allocation for each element)
    type(site_fluxdiags_type), pointer :: flux_diags(:)            ! flux diagnostics (allocation for each element)
    type(ed_site_hydr_type),   pointer :: si_hydr                  ! plant hydraulics

    !-------------------------------------------------------------------------------------
     
    ! RESOURCE MANAGEMENT
    type(fates_resources_management_type) :: resources_management ! resources_management at the site

    !-------------------------------------------------------------------------------------

    ! INDICES
    integer  :: h_gid ! global index of this site in the history output file
    real(r8) :: lat   ! latitude [degrees] 
    real(r8) :: lon   ! longitude [degrees] 

    !-------------------------------------------------------------------------------------

    ! FIXED BIOGEOGRAPHY MODE INPUTS
    real(r8), allocatable :: area_PFT(:)     ! area allocated to individual PFTs    
    integer,  allocatable :: use_this_pft(:) ! is area_PFT > 0?

    !-------------------------------------------------------------------------------------

    ! SP MODE TARGET PFT-LEVEL VARIABLES
    real(r8), allocatable :: sp_tlai(:) ! target TLAI per FATES pft
    real(r8), allocatable :: sp_tsai(:) ! target TSAI per FATES pft
    real(r8), allocatable :: sp_htop(:) ! target HTOP per FATES pft

    !-------------------------------------------------------------------------------------

    ! SITE CHARACTERISTICS
    real(r8), allocatable :: area_by_age(:) ! total area of patches in each age bin [m2]
    real(r8), allocatable :: rec_l2fr(:,:)  ! running mean of the l2fr's for the newly
                                            !   recruited, pft x canopy_layer
    real(r8) :: ema_npp                     ! exponential moving average of NPP [gC/m2/year]
                                            !   The length scale is hard-coded "ema_npp_tcale"
                                            !   in FatesSoilBGCFluxMod. Used solely to inform bc_out%ema_npp
                                            !   which is used for fixation
    real(r8) :: snow_depth                  ! site-level snow depth [m] (used for ELAI/TLAI calcs)
    real(r8) :: spread                      ! canopy spread, dynamic canopy allometric term [unitless]

    !-------------------------------------------------------------------------------------

    ! PHENOLOGY 
    real(r8) :: grow_deg_days                   ! growing degree days
    integer  :: nchilldays                      ! number of chilling days (for botta gdd threshhold calculation)
    integer  :: ncolddays                       ! number of cold days (must exceed threshold to drop leaves)
    integer  :: cstatus                         ! are leaves in this pixel on or off for cold deciduous
                                                !   0 = this site has not experienced a cold period over at least
                                                !     400 days, leaves are dropped and flagged as non-cold region
                                                !   1 = this site is in a cold-state where leaves should have fallen
                                                !   2 = this site is in a warm-state where leaves are allowed to flush
    integer  :: dstatus                         ! are leaves in this pixel on or off for drought deciduous
                                                !   0 = leaves off due to time exceedance
                                                !   1 = leaves off due to moisture availibility
                                                !   2 = leaves on due to moisture avail
                                                !   3 = leaves on due to time exceedance
    integer  :: cleafondate                     ! model date of leaf on (cold) [day integer]
    integer  :: cleafoffdate                    ! model date of leaf off (cold) [day integer]
    integer  :: dleafondate                     ! model date of leaf on (drought) [day integer]
    integer  :: dleafoffdate                    ! model date of leaf off (drought) [day integer]
    integer  :: phen_model_date                 ! current model date [day integer]
                                                !   this date stays continuous when
                                                !   in runs that are restarted, regardless of
                                                !   the conditions of restart
    real(r8) :: water_memory(numWaterMem)       ! last numwatermem days of soil moisture memory [m3/m3]
    real(r8) :: vegtemp_memory(num_vegtemp_mem) ! last num_vegtemp_mem days temperature for senescence model [deg C]

    !-------------------------------------------------------------------------------------

    ! FIRE
    real(r8) :: wind          ! daily wind speed [m/min]
    real(r8) :: acc_ni        ! daily nesterov index accumulating over time
    real(r8) :: fdi           ! daily probability an ignition event will start a fire [0-1]
    real(r8) :: NF            ! daily ignitions [/km2]
    real(r8) :: NF_successful ! daily ignitions that actually lead to fire [/km2]

    !-------------------------------------------------------------------------------------

    ! SOIL LAYERS
    integer               :: nlevsoil        ! number of soil layers in this site
    real(r8), allocatable :: zi_soil(:)      ! interface level below a "z" level [m]
                                             !   this contains a zero index for surface.
    real(r8), allocatable :: dz_soil(:)      ! layer thickness [m]
    real(r8), allocatable :: z_soil(:)       ! layer depth [m]
    real(r8), allocatable :: rootfrac_scr(:) ! This is just allocated scratch space to hold
                                             !   root fractions. Since root fractions may be dependent
                                             !   on cohort properties, and we do not want to store this infromation
                                             !   on each cohort, we do not keep root fractions in
                                             !   memory, and instead calculate them on demand.
                                             !   This array is allocated over the number of soil
                                             !   layers for each site, and save allocating deallocating.
                                             !   NOTE: THIS SCRATCH SPACE WOULD NOT BE THREAD-SAFE
                                             !   IF WE FORK ON PATCHES

    !-------------------------------------------------------------------------------------

    ! TERMINATION, RECRUITMENT, DEMOTION, and DISTURBANCE DIAGNOSTICS
    real(r8)              :: term_crownarea_canopy             ! crown area of individuals killed due to termination mortality for canopy [m2]
    real(r8)              :: term_crownarea_ustory             ! crown area of individuals killed due to termination mortality for understory [m2] 
    real(r8)              :: imort_crownarea                   ! crown area of individuals killed due to impact mortality per year [m2/day]
    real(r8)              :: fmort_crownarea_canopy            ! crown area of canopy indivs killed due to fire per year [m2/sec]
    real(r8)              :: fmort_crownarea_ustory            ! crown area of understory indivs killed due to fire per year [m2/sec] 
    real(r8), allocatable :: term_nindivs_canopy(:,:)          ! number of canopy individuals that were in cohorts which 
                                                               !   were terminated this timestep (size x pft) [/day]
    real(r8), allocatable :: term_nindivs_ustory(:,:)          ! number of understory individuals that were in cohorts which 
                                                               !   were terminated this timestep (size x pft) [/day]
    real(r8), allocatable :: term_carbonflux_canopy(:)         ! carbon flux from live to dead pools associated 
                                                               !   with termination mortality, per canopy level. [kgC/ha/day]
    real(r8), allocatable :: term_carbonflux_ustory(:)         ! carbon flux from live to dead pools associated 
                                                               !  with termination mortality, per canopy level.  [kgC/ha/day]    
    real(r8), allocatable :: imort_carbonflux(:)               ! biomass of individuals killed due to impact mortality per year. [kgC/m2/sec]
    real(r8), allocatable :: fmort_carbonflux_canopy(:)        ! biomass of canopy indivs killed due to fire per year. [gC/m2/sec]
    real(r8), allocatable :: fmort_carbonflux_ustory(:)        ! biomass of understory indivs killed due to fire per year [gC/m2/sec] 
    real(r8), allocatable :: term_abg_flux(:,:)                ! aboveground biomass lost due to termination mortality x size x pft
    real(r8), allocatable :: imort_abg_flux(:,:)               ! aboveground biomass lost due to impact mortality x size x pft [kgC/m2/sec]
    real(r8), allocatable :: fmort_abg_flux(:,:)               ! aboveground biomass lost due to fire mortality x size x pft
    real(r8)              :: demotion_carbonflux               ! biomass of demoted individuals from canopy to understory [kgC/ha/day]
    real(r8)              :: promotion_carbonflux              ! biomass of promoted individuals from understory to canopy [kgC/ha/day]
    real(r8)              :: recruitment_rate(1:maxpft)        ! number of individuals that were recruited into new cohorts
    real(r8), allocatable :: demotion_rate(:)                  ! rate of individuals demoted from canopy to understory per FATES timestep
    real(r8), allocatable :: promotion_rate(:)                 ! rate of individuals promoted from understory to canopy per FATES timestep
    real(r8), allocatable :: imort_rate(:,:)                   ! rate of individuals killed due to impact mortality per year.  on size x pft array
    real(r8), allocatable :: fmort_rate_canopy(:,:)            ! rate of canopy individuals killed due to fire mortality per year.  
                                                               !  on size x pft array  (1:nlevsclass,1:numpft)
    real(r8), allocatable :: fmort_rate_ustory(:,:)            ! rate of understory individuals killed due to fire mortality per year.  
                                                               !  on size x pft array  (1:nlevsclass,1:numpft)
    real(r8), allocatable :: fmort_rate_cambial(:,:)           ! rate of individuals killed due to fire mortality 
                                                               !  from cambial damage per year.  on size x pft array
    real(r8), allocatable :: fmort_rate_crown(:,:)             ! rate of individuals killed due to fire mortality 
                                                               !  from crown damage per year.  on size x pft array
    real(r8), allocatable :: imort_rate_damage(:,:,:)          ! number of individuals per damage class that die from impact mortality
    real(r8), allocatable :: term_nindivs_canopy_damage(:,:,:) ! number of individuals per damage class that die from termination mortality - canopy
    real(r8), allocatable :: term_nindivs_ustory_damage(:,:,:) ! number of individuals per damage class that die from termination mortality - canopy
    real(r8), allocatable :: fmort_rate_canopy_damage(:,:,:)   ! number of individuals per damage class that die from fire - canopy
    real(r8), allocatable :: fmort_rate_ustory_damage(:,:,:)   ! number of individuals per damage class that die from fire - ustory
    real(r8), allocatable :: fmort_cflux_canopy_damage(:,:)    ! cflux per damage class that die from fire - canopy
    real(r8), allocatable :: fmort_cflux_ustory_damage(:,:)    ! cflux per damage class that die from fire - ustory
    real(r8), allocatable :: imort_cflux_damage(:,:)           ! carbon flux from impact mortality by damage class [kgC/m2/sec]
    real(r8), allocatable :: term_cflux_canopy_damage(:,:)     ! carbon flux from termination mortality by damage class
    real(r8), allocatable :: term_cflux_ustory_damage(:,:)     ! carbon flux from termination mortality by damage class [kgC]
    real(r8), allocatable :: growthflux_fusion(:,:)            ! rate of individuals moving into a given size class bin
                                                               !  due to fusion in a given day (size x pft)
    real(r8)              :: crownarea_canopy_damage           ! crown area of canopy that is damaged annually  
    real(r8)              :: crownarea_ustory_damage           ! crown area of understory that is damaged annually

    !-------------------------------------------------------------------------------------

    ! ACTUAL AND POTENTIAL DISTURBANCE RATES
    real(r8) :: disturbance_rates_primary_to_primary(N_DIST_TYPES)     ! actual disturbance rates from primary patches to primary patches [m2/m2/day]
    real(r8) :: disturbance_rates_primary_to_secondary(N_DIST_TYPES)   ! actual disturbance rates from primary patches to secondary patches [m2/m2/day]
    real(r8) :: disturbance_rates_secondary_to_secondary(N_DIST_TYPES) ! actual disturbance rates from secondary patches to secondary patches [m2/m2/day]
    real(r8) :: potential_disturbance_rates(N_DIST_TYPES)              ! "potential" disturb rates (i.e. prior to the "which is most" logic) [m2/m2/day]
    real(r8) :: primary_land_patchfusion_error                         ! error term in total area of primary patches associated with patch fusion [m2/m2/day]
    
    !-------------------------------------------------------------------------------------

    ! Mineralized nutrient flux from veg to the soil, via multiple mechanisms
    ! inluding symbiotic fixation, or other 

    !real(r8) :: allocatable :: minn_flux_out  ! kg/ha/day
    !real(r8) :: allocatable :: minp_flux_out  ! kg/ha/day

    !-------------------------------------------------------------------------------------

    contains 

    procedure :: Init
    procedure :: NanValues
    procedure :: ZeroValues
    procedure :: Dump

  end type fates_site_type

  !=======================================================================================

  type, public :: site_fluxdiags_type

    ! ------------------------------------------------------------------------------------
    ! Diagnostics for fluxes into the litter pool from plants
    ! these fluxes are the total from 
    ! (1) turnover from living plants
    ! (2) mass transfer from non-disturbance inducing mortality events
    ! (3) mass transfer from disturbance inducing mortality events
    ! [kg / ha / day]
    ! ------------------------------------------------------------------------------------

    real(r8)              :: cwd_ag_input(1:ncwd)               
    real(r8)              :: cwd_bg_input(1:ncwd)               
    real(r8), allocatable :: leaf_litter_input(:)
    real(r8), allocatable :: root_litter_input(:)
  
  contains

    procedure :: ZeroFluxDiags
  
  end type site_fluxdiags_type

  !=======================================================================================

  contains

    subroutine Init(this, num_levsoil, zi_sisl, dz_sisl, z_sisl)
      !
      !  DESCRIPTION:
      !  Initialize variables for the site type
      !

      ! ARGUMENTS:
      class(fates_site_type), intent(inout) :: this        ! site object
      integer,               intent(in)     :: num_levsoil ! the number of soil layers in this column
      real(r8),               intent(in)    :: zi_sisl(:)  ! soil interface level below a "z" level [m]
      real(r8),               intent(in)    :: dz_sisl(:)  ! soil layer thickness [m]
      real(r8),               intent(in)    :: z_sisl(:)   ! soil layer depth [m]

      ! LOCALS:
      integer :: el ! looping index

      ! allocate arrays
      allocate(this%term_nindivs_canopy(1:nlevsclass,1:numpft))
      allocate(this%term_nindivs_ustory(1:nlevsclass,1:numpft))
      allocate(this%demotion_rate(1:nlevsclass))
      allocate(this%promotion_rate(1:nlevsclass))
      allocate(this%imort_rate(1:nlevsclass,1:numpft))
      allocate(this%fmort_rate_canopy(1:nlevsclass,1:numpft))
      allocate(this%fmort_rate_ustory(1:nlevsclass,1:numpft))
      allocate(this%fmort_rate_cambial(1:nlevsclass,1:numpft))
      allocate(this%fmort_rate_crown(1:nlevsclass,1:numpft))
      allocate(this%growthflux_fusion(1:nlevsclass,1:numpft))
      allocate(this%mass_balance(1:num_elements))
      allocate(this%flux_diags(1:num_elements))

      if (hlm_use_tree_damage .eq. itrue) then 
        allocate(this%term_nindivs_canopy_damage(1:nlevdamage,1:nlevsclass,1:numpft))
        allocate(this%term_nindivs_ustory_damage(1:nlevdamage,1:nlevsclass,1:numpft))
        allocate(this%imort_rate_damage(1:nlevdamage,1:nlevsclass,1:numpft))
        allocate(this%imort_cflux_damage(1:nlevdamage,1:nlevsclass))
        allocate(this%term_cflux_canopy_damage(1:nlevdamage,1:nlevsclass))
        allocate(this%term_cflux_ustory_damage(1:nlevdamage,1:nlevsclass))
        allocate(this%fmort_rate_canopy_damage(1:nlevdamage,1:nlevsclass,1:numpft))
        allocate(this%fmort_rate_ustory_damage(1:nlevdamage,1:nlevsclass,1:numpft)) 
        allocate(this%fmort_cflux_canopy_damage(1:nlevdamage,1:nlevsclass))
        allocate(this%fmort_cflux_ustory_damage(1:nlevdamage,1:nlevsclass)) 
      else
        allocate(this%term_nindivs_canopy_damage(1,1,1))
        allocate(this%term_nindivs_ustory_damage(1,1,1))
        allocate(this%imort_rate_damage(1,1,1))
        allocate(this%imort_cflux_damage(1,1))
        allocate(this%term_cflux_canopy_damage(1,1))
        allocate(this%term_cflux_ustory_damage(1,1))
        allocate(this%fmort_rate_canopy_damage(1,1,1))
        allocate(this%fmort_rate_ustory_damage(1,1,1))
        allocate(this%fmort_cflux_canopy_damage(1,1))
        allocate(this%fmort_cflux_ustory_damage(1,1))
      end if

      allocate(this%term_carbonflux_canopy(1:numpft))
      allocate(this%term_carbonflux_ustory(1:numpft))
      allocate(this%imort_carbonflux(1:numpft))
      allocate(this%fmort_carbonflux_canopy(1:numpft))
      allocate(this%fmort_carbonflux_ustory(1:numpft))

      allocate(this%term_abg_flux(1:nlevsclass,1:numpft))
      allocate(this%imort_abg_flux(1:nlevsclass,1:numpft))
      allocate(this%fmort_abg_flux(1:nlevsclass,1:numpft))

      allocate(this%rootfrac_scr(num_levsoil))
      allocate(this%zi_soil(0:num_levsoil))
      allocate(this%dz_soil(num_levsoil))
      allocate(this%z_soil(num_levsoil))

      if (hlm_use_nocomp .eq. itrue .and. hlm_use_fixed_biogeog .eq. itrue) then
        ! SP and nocomp require a bare-ground patch.
        allocate(this%area_pft(0:numpft))
      else  
        allocate(this%area_pft(1:numpft))  
      endif

      allocate(this%use_this_pft(1:numpft))
      allocate(this%area_by_age(1:nlevage))

      ! for CNP dynamics, track the mean l2fr of recruits
      ! for different pfts and canopy positions
      allocate(this%rec_l2fr(1:numpft,nclmax))

      ! SP mode
      allocate(this%sp_tlai(1:numpft))
      allocate(this%sp_tsai(1:numpft))
      allocate(this%sp_htop(1:numpft))

      do el = 1, num_elements
        allocate(this%flux_diags(el)%leaf_litter_input(1:numpft))
        allocate(this%flux_diags(el)%root_litter_input(1:numpft))
      end do

      ! initialize all values to nan
      call this%NanValues()

      ! set some values to zero
      call this%ZeroValues()

      ! initialize the static soil arrays from the boundary (initial) condition
      this%nlevsoil = num_levsoil
      this%zi_soil(:) = zi_sisl(:)
      this%dz_soil(:) = dz_sisl(:)
      this%z_soil(:)  = z_sisl(:)

    end subroutine Init

    !=====================================================================================

    subroutine NanValues(this)
      !
      !  DESCRIPTION:
      !  Sets all values in site to nan
      !

      ! ARGUMENTS:
      class(fates_site_type), intent(inout) :: this

      ! set pointers to null
      this%oldest_patch    => null()   
      this%youngest_patch  => null()
      nullify(this%oldest_patch)
      nullify(this%youngest_patch)

      ! INDICES
      !this%h_gid                                      = fates_unset_int
      !this%lat                                        = nan 
      !this%lon                                        = nan 

      ! FIXED BIOGEOGRAPHY MODE INPUTS
      this%area_PFT(:)                                = nan
      this%use_this_pft(:)                            = fates_unset_int

      ! SP MODE TARGET PFT-LEVEL VARIABLES
      this%sp_tlai(:)                                 = nan
      this%sp_tsai(:)                                 = nan 
      this%sp_htop(:)                                 = nan 

      ! SITE CHARACTERISTICS
      this%area_by_age(:)                              = nan 
      this%rec_l2fr(:,:)                               = nan
      !this%ema_npp                                     = nan 
      this%snow_depth                                  = nan 
      this%spread                                      = nan 

      ! PHENOLOGY 
      this%grow_deg_days                               = nan                
      this%nchilldays                                  = fates_unset_int                 
      this%ncolddays                                   = fates_unset_int                  
      this%cstatus                                     = fates_unset_int                       
      this%dstatus                                     = fates_unset_int                      
      this%cleafondate                                 = fates_unset_int          
      this%cleafoffdate                                = fates_unset_int              
      this%dleafondate                                 = fates_unset_int 
      this%dleafoffdate                                = fates_unset_int             
      this%phen_model_date                             = fates_unset_int 
      this%water_memory(:)                             = nan      
      this%vegtemp_memory(:)                           = nan 

      ! FIRE
      this%wind                                        = nan          
      this%acc_ni                                      = nan     
      this%fdi                                         = nan          
      this%NF                                          = nan        
      this%NF_successful                               = nan 

      ! SOIL LAYERS
      this%nlevsoil                                    = fates_unset_int 
      this%zi_soil(:)                                  = nan       
      this%dz_soil(:)                                  = nan     
      this%z_soil(:)                                   = nan     
      this%rootfrac_scr(:)                             = nan 

      ! TERMINATION, RECRUITMENT, DEMOTION, and DISTURBANCE DIAGNOSTICS
      this%term_crownarea_canopy                       = nan       
      this%term_crownarea_ustory                       = nan            
      this%imort_crownarea                             = nan                  
      this%fmort_crownarea_canopy                      = nan         
      this%fmort_crownarea_ustory                      = nan           
      this%term_nindivs_canopy(:,:)                    = nan                                 
      this%term_nindivs_ustory(:,:)                    = nan                                                                      
      this%term_carbonflux_canopy(:)                   = nan                                                               
      this%term_carbonflux_ustory(:)                   = nan                                                        
      this%imort_carbonflux(:)                         = nan               
      this%fmort_carbonflux_canopy(:)                  = nan         
      this%fmort_carbonflux_ustory(:)                  = nan         
      this%term_abg_flux(:,:)                          = nan              
      this%imort_abg_flux(:,:)                         = nan               
      this%fmort_abg_flux(:,:)                         = nan               
      this%demotion_carbonflux                         = nan              
      this%promotion_carbonflux                        = nan             
      this%recruitment_rate(:)                         = nan         
      this%demotion_rate(:)                            = nan                  
      this%promotion_rate(:)                           = nan                  
      this%imort_rate(:,:)                             = nan                    
      this%fmort_rate_canopy(:,:)                      = nan                                                                         
      this%fmort_rate_ustory(:,:)                      = nan                                                                         
      this%fmort_rate_cambial(:,:)                     = nan                                                             
      this%fmort_rate_crown(:,:)                       = nan                   
      this%imort_rate_damage(:,:,:)                    = nan           
      this%term_nindivs_canopy_damage(:,:,:)           = nan  
      this%term_nindivs_ustory_damage(:,:,:)           = nan  
      this%fmort_rate_canopy_damage(:,:,:)             = nan     
      this%fmort_rate_ustory_damage(:,:,:)             = nan    
      this%fmort_cflux_canopy_damage(:,:)              = nan     
      this%fmort_cflux_ustory_damage(:,:)              = nan     
      this%imort_cflux_damage(:,:)                     = nan            
      this%term_cflux_canopy_damage(:,:)               = nan   
      this%term_cflux_ustory_damage(:,:)               = nan      
      this%growthflux_fusion(:,:)                      = nan        
      this%crownarea_canopy_damage                     = nan             
      this%crownarea_ustory_damage                     = nan       

      ! ACTUAL AND POTENTIAL DISTURBANCE RATES
      this%disturbance_rates_primary_to_primary(:)     = nan     
      this%disturbance_rates_primary_to_secondary(:)   = nan     
      this%disturbance_rates_secondary_to_secondary(:) = nan  
      this%potential_disturbance_rates(:)              = nan               
      this%primary_land_patchfusion_error              = nan     

    end subroutine NanValues

    !=====================================================================================

    subroutine ZeroValues(this)
      !
      ! DESCRIPTION:
      !     sets specific variables in the site to zero
      !

      ! ARGUMENTS:
      class(fates_site_type), intent(inout) :: this

      ! LOCALS:
      integer :: el ! looping index

      ! zero the state variables used for checking mass conservation
      do el = 1, num_elements
        call this%mass_balance(el)%ZeroMassBalState()
        call this%mass_balance(el)%ZeroMassBalFlux()
        call this%flux_diags(el)%ZeroFluxDiags()
      end do

      ! zero resources management variables
      call this%resources_management%ZeroVals()

      this%primary_land_patchfusion_error              = 0.0_r8
      this%potential_disturbance_rates(:)              = 0.0_r8
      this%disturbance_rates_secondary_to_secondary(:) = 0.0_r8
      this%disturbance_rates_primary_to_secondary(:)   = 0.0_r8
      this%disturbance_rates_primary_to_primary(:)     = 0.0_r8
      this%acc_ni                                      = 0.0_r8    
      this%FDI                                         = 0.0_r8     
      this%NF                                          = 0.0_r8     
      this%NF_successful                               = 0.0_r8     
      this%term_nindivs_canopy(:,:)                    = 0.0_r8
      this%term_nindivs_ustory(:,:)                    = 0.0_r8
      this%term_crownarea_canopy                       = 0.0_r8
      this%term_crownarea_ustory                       = 0.0_r8
      this%imort_crownarea                             = 0.0_r8
      this%fmort_crownarea_canopy                      = 0.0_r8
      this%fmort_crownarea_ustory                      = 0.0_r8
      this%term_carbonflux_canopy(:)                   = 0.0_r8
      this%term_carbonflux_ustory(:)                   = 0.0_r8
      this%recruitment_rate(:)                         = 0.0_r8
      this%imort_rate(:,:)                             = 0.0_r8
      this%imort_carbonflux(:)                         = 0.0_r8
      this%fmort_rate_canopy(:,:)                      = 0.0_r8
      this%fmort_rate_ustory(:,:)                      = 0.0_r8
      this%fmort_carbonflux_canopy(:)                  = 0.0_r8
      this%fmort_carbonflux_ustory(:)                  = 0.0_r8
      this%fmort_rate_cambial(:,:)                     = 0.0_r8
      this%fmort_rate_crown(:,:)                       = 0.0_r8
      this%term_abg_flux(:,:)                          = 0.0_r8
      this%imort_abg_flux(:,:)                         = 0.0_r8
      this%fmort_abg_flux(:,:)                         = 0.0_r8
      this%growthflux_fusion(:,:)                      = 0.0_r8
      this%demotion_rate(:)                            = 0.0_r8
      this%demotion_carbonflux                         = 0.0_r8
      this%promotion_rate(:)                           = 0.0_r8
      this%promotion_carbonflux                        = 0.0_r8
      this%imort_rate_damage(:,:,:)                    = 0.0_r8
      this%term_nindivs_canopy_damage(:,:,:)           = 0.0_r8
      this%term_nindivs_ustory_damage(:,:,:)           = 0.0_r8
      this%imort_cflux_damage(:,:)                     = 0.0_r8
      this%term_cflux_canopy_damage(:,:)               = 0.0_r8
      this%term_cflux_ustory_damage(:,:)               = 0.0_r8
      this%crownarea_canopy_damage                     = 0.0_r8
      this%crownarea_ustory_damage                     = 0.0_r8
      this%fmort_rate_canopy_damage(:,:,:)             = 0.0_r8
      this%fmort_rate_ustory_damage(:,:,:)             = 0.0_r8
      this%fmort_cflux_canopy_damage(:,:)              = 0.0_r8
      this%fmort_cflux_ustory_damage(:,:)              = 0.0_r8
      this%spread                                      = 0.0_r8
      this%area_pft(:)                                 = 0.0_r8
      this%area_by_age(:)                              = 0.0_r8

    end subroutine ZeroValues
  
    ! ====================================================================================

    subroutine Dump(this) 
      !
      !  DESCRIPTION:
      !  Print out relevant information for a site
      !
  
      ! ARGUMENTS:
      class(fates_site_type), intent(in) :: this ! site object
  
      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) ' Site Coordinates                       '
      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) 'latitude                    = ', this%lat
      write(fates_log(),*) 'longitude                   = ', this%lon
      write(fates_log(),*) '----------------------------------------'
  
    end subroutine Dump

    !=====================================================================================

    subroutine ZeroFluxDiags(this)
      !
      !  DESCRIPTION:
      !  Zero all variables for the flux diags type
      !
      
      ! ARGUMENTS:
      class(site_fluxdiags_type) :: this ! flux diags object
      
      this%cwd_ag_input(:)      = 0.0_r8
      this%cwd_bg_input(:)      = 0.0_r8
      this%leaf_litter_input(:) = 0.0_r8
      this%root_litter_input(:) = 0.0_r8
      
    end subroutine ZeroFluxDiags

    !=====================================================================================

end module FatesSiteMod
