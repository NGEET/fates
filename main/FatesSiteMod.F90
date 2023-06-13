module FatesSiteMod 

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesPatchMod,     only : fates_patch_type
  use EDTypesMod,        only : ed_resources_management_type
  use EDTypesMod,        only : site_fluxdiags_type
  use EDTypesMod,        only : site_massbal_type

  implicit none
  private

  type, public :: fates_site_type
     
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
                                                         ! with termination mortality, per canopy level. [kgC/ha/day]
     real(r8), allocatable :: term_carbonflux_ustory(:)  ! carbon flux from live to dead pools associated 
                                                         ! with termination mortality, per canopy level.  [kgC/ha/day]    
     real(r8), allocatable :: imort_carbonflux(:)        ! biomass of individuals killed due to impact mortality per year. [kgC/m2/sec]
     real(r8), allocatable :: fmort_carbonflux_canopy(:) ! biomass of canopy indivs killed due to fire per year. [gC/m2/sec]
     real(r8), allocatable :: fmort_carbonflux_ustory(:) ! biomass of understory indivs killed due to fire per year [gC/m2/sec] 

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

     ! site-level variables to keep track of the disturbance rates, both actual and "potential"
     real(r8) :: disturbance_rates_primary_to_primary(N_DIST_TYPES)      ! actual disturbance rates from primary patches to primary patches [m2/m2/day]
     real(r8) :: disturbance_rates_primary_to_secondary(N_DIST_TYPES)    ! actual disturbance rates from primary patches to secondary patches [m2/m2/day]
     real(r8) :: disturbance_rates_secondary_to_secondary(N_DIST_TYPES)  ! actual disturbance rates from secondary patches to secondary patches [m2/m2/day]
     real(r8) :: potential_disturbance_rates(N_DIST_TYPES)               ! "potential" disturb rates (i.e. prior to the "which is most" logic) [m2/m2/day]
     real(r8) :: primary_land_patchfusion_error                          ! error term in total area of primary patches associated with patch fusion [m2/m2/day]
     
  end type fates_site_type

end module FatesSiteMod