module FatesPatchMod

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesCohortMod,      only : fates_cohort_type
  use FatesRunningMeanMod, only : rmean_type
  use FatesLitterMod,      only : nfsc
  use FatesLitterMod,      only : litter_type
  use EDParamsMod,         only : maxSWb, nlevleaf, nclmax
  use FatesConstantsMod,   only : n_dbh_bins, maxpft, n_dist_types
  use FatesConstantsMod,   only : n_rad_stream_types

  implicit none
  private

  type, public :: fates_patch_type

    ! POINTERS
    type (fates_cohort_type), pointer :: tallest => null()  ! pointer to patch's tallest cohort    
    type (fates_cohort_type), pointer :: shortest => null() ! pointer to patch's shortest cohort
    type (fates_patch_type),  pointer :: older => null()    ! pointer to next older patch   
    type (fates_patch_type),  pointer :: younger => null()  ! pointer to next younger patch
  
    !:.........................................................................:

    ! INDICES
    integer  :: patchno          ! unique number given to each new patch created for tracking
    integer  :: nocomp_pft_label ! when nocomp is active, use this label for patch ID
                                 !   each patch ID corresponds to a pft number since each
                                 !   patch has only one pft.  Bareground patches are given
                                 !   a zero integer as a label. If nocomp is not active this 
                                 !   is set to unset. This is set in patch%create as an argument
                                 !   to that procedure.

    !:.........................................................................:

    ! PATCH INFO
    real(r8) :: age                          ! average patch age [years]                  
    integer  :: age_class                    ! age class of the patch for history binning purposes
    real(r8) :: area                         ! patch area [m2]
    integer  :: countcohorts                 ! number of cohorts in patch
    integer  :: ncl_p                        ! number of occupied canopy layers
    integer  :: anthro_disturbance_label     ! patch label for anthropogenic disturbance classification
    real(r8) :: age_since_anthro_disturbance ! average age for secondary forest since last anthropogenic disturbance [years]

    !:.........................................................................:

    ! RUNNING MEANS
    !class(rmean_type), pointer :: t2m          ! place-holder for 2m air temperature (variable window-size)
    class(rmean_type), pointer :: tveg24        ! 24-hour mean vegetation temperature [K]
    class(rmean_type), pointer :: tveg_lpa      ! running mean of vegetation temperature at the
                                                !   leaf photosynthesis acclimation timescale [K]
    class(rmean_type), pointer :: tveg_longterm ! long-term running mean of vegetation temperature at the
                                                !   leaf photosynthesis acclimation timescale [K] (i.e T_home)

    !:.........................................................................:

    ! LEAF ORGANIZATION
    real(r8) :: pft_agb_profile(maxpft,n_dbh_bins)          ! binned aboveground biomass, for patch fusion [kgC/m2]
    real(r8) :: canopy_layer_tlai(nclmax)                   ! total leaf area index of each canopy layer [m2 veg/m2 canopy area]
                                                            !   (patch without bare ground)
                                                            !   used to determine attenuation of parameters during photosynthesis
    real(r8) :: total_canopy_area                           ! area that is covered by vegetation [m2]
    real(r8) :: total_tree_area                             ! area that is covered by woody vegetation [m2]
    real(r8) :: zstar                                       ! height of smallest canopy tree, only meaningful in "strict PPA" mode [m]
    real(r8) :: elai_profile(nclmax,maxpft,nlevleaf)        ! exposed leaf area in each canopy layer, pft, and leaf layer [m2 leaf/m2 contributing crown area]
    real(r8) :: esai_profile(nclmax,maxpft,nlevleaf)        ! exposed stem area in each canopy layer, pft, and leaf layer [m2 leaf/m2 contributing crown area]
    real(r8) :: tlai_profile(nclmax,maxpft,nlevleaf)
    real(r8) :: tsai_profile(nclmax,maxpft,nlevleaf)
    real(r8) :: canopy_area_profile(nclmax,maxpft,nlevleaf) ! fraction of crown area per canopy area in each layer
                                                            !   they will sum to 1.0 in the fully closed canopy layers
                                                            !   but only in leaf-layers that contain contributions
                                                            !   from all cohorts that donate to canopy_area
    integer  :: canopy_mask(nclmax,maxpft)                  ! is there any of this pft in this canopy layer?      
    integer  :: nrad(nclmax,maxpft)                         ! number of exposed leaf layers for each canopy layer and pft
    integer  :: ncan(nclmax,maxpft)                         ! number of total   leaf layers for each canopy layer and pft
    real(r8) :: c_stomata                                   ! mean stomatal conductance of all leaves in the patch   [umol/m2/s]
    real(r8) :: c_lblayer                                   ! mean boundary layer conductance of all leaves in the patch [umol/m2/s]
    real(r8) ::  layer_height_profile(nclmax,maxpft,nlevleaf)

    !:.........................................................................:

    real(r8) :: psn_z(nclmax,maxpft,nlevleaf) 
    real(r8) :: nrmlzd_parprof_pft_dir_z(n_rad_stream_types,nclmax,maxpft,nlevleaf)
    real(r8) :: nrmlzd_parprof_pft_dif_z(n_rad_stream_types,nclmax,maxpft,nlevleaf)
    real(r8) :: nrmlzd_parprof_dir_z(n_rad_stream_types,nclmax,nlevleaf)
    real(r8) :: nrmlzd_parprof_dif_z(n_rad_stream_types,nclmax,nlevleaf)

    !:.........................................................................:

    ! RADIATION
    real(r8) :: radiation_error                           ! radiation error [W/m2] 
    real(r8) :: fcansno                                   ! fraction of canopy covered in snow [0-1]
    logical  :: solar_zenith_flag                         ! integer flag specifying daylight (based on zenith angle)
    real(r8) :: solar_zenith_angle                        ! solar zenith angle [radians]
    real(r8) :: gnd_alb_dif(maxSWb)                       ! ground albedo for diffuse rad, both bands [0-1]
    real(r8) :: gnd_alb_dir(maxSWb)                       ! ground albedo for direct rad, both bands [0-1]

    ! organized by canopy layer, pft, and leaf layer
    real(r8) :: fabd_sun_z(nclmax,maxpft,nlevleaf)        ! sun fraction of direct light absorbed [0-1]
    real(r8) :: fabd_sha_z(nclmax,maxpft,nlevleaf)        ! shade fraction of direct light absorbed [0-1]
    real(r8) :: fabi_sun_z(nclmax,maxpft,nlevleaf)        ! sun fraction of indirect light absorbed [0-1]
    real(r8) :: fabi_sha_z(nclmax,maxpft,nlevleaf)        ! shade fraction of indirect light absorbed [0-1]
    real(r8) :: ed_parsun_z(nclmax,maxpft,nlevleaf)       ! PAR absorbed in the sun [W/m2]   
    real(r8) :: ed_parsha_z(nclmax,maxpft,nlevleaf)       ! PAR absorbed in the shade [W/m2]
    real(r8) :: ed_laisun_z(nclmax,maxpft,nlevleaf)
    real(r8) :: ed_laisha_z(nclmax,maxpft,nlevleaf)
    real(r8) :: f_sun(nclmax,maxpft,nlevleaf)             ! fraction of leaves in the sun [0-1]

    ! radiation profiles for comparison against observations
    real(r8) :: parprof_pft_dir_z(nclmax,maxpft,nlevleaf) ! direct-beam PAR profile through canopy, by canopy, PFT, leaf level [W/m2]
    real(r8) :: parprof_pft_dif_z(nclmax,maxpft,nlevleaf) ! diffuse     PAR profile through canopy, by canopy, PFT, leaf level [W/m2]
    real(r8) :: parprof_dir_z(nclmax,nlevleaf)            ! direct-beam PAR profile through canopy, by canopy, leaf level [W/m2]
    real(r8) :: parprof_dif_z(nclmax,nlevleaf)            ! diffuse     PAR profile through canopy, by canopy, leaf level [W/m2]
    
    real(r8), allocatable :: tr_soil_dir(:)               ! fraction of incoming direct radiation transmitted to the soil as direct, by numSWB [0-1]
    real(r8), allocatable :: tr_soil_dif(:)               ! fraction of incoming diffuse radiation that is transmitted to the soil as diffuse [0-1]
    real(r8), allocatable :: tr_soil_dir_dif(:)           ! fraction of incoming direct radiation that is transmitted to the soil as diffuse [0-1]
    real(r8), allocatable :: fab(:)                       ! fraction of incoming total   radiation that is absorbed by the canopy
    real(r8), allocatable :: fabd(:)                      ! fraction of incoming direct  radiation that is absorbed by the canopy
    real(r8), allocatable :: fabi(:)                      ! fraction of incoming diffuse radiation that is absorbed by the canopy
    real(r8), allocatable :: sabs_dir(:)                  ! fraction of incoming direct  radiation that is absorbed by the canopy
    real(r8), allocatable :: sabs_dif(:)                  ! fraction of incoming diffuse radiation that is absorbed by the canopy

    !:.........................................................................:

    ! ROOTS
    real(r8) :: btran_ft(maxpft)       ! btran calculated seperately for each PFT
    real(r8) :: bstress_sal_ft(maxpft) ! bstress from salinity calculated seperately for each PFT
    
    !:.........................................................................:

    ! EXTERNAL SEED RAIN
    real(r8) :: nitr_repro_stoich(maxpft) ! The NC ratio of a new recruit in this patch
    real(r8) :: phos_repro_stoich(maxpft) ! The PC ratio of a new recruit in this patch

    !:.........................................................................:
    
   ! DISTURBANCE 
    real(r8) :: disturbance_rates(n_dist_types) ! disturbance rate [0-1/day] from 1) mortality 
                                                !                                 2) fire 
                                                !                                 3) logging mortatliy
    real(r8) :: fract_ldist_not_harvested       ! fraction of logged area that is canopy trees that weren't harvested [0-1]

    !:.........................................................................:

    ! LITTER AND COARSE WOODY DEBRIS
    type(litter_type), pointer :: litter(:)               ! litter (leaf,fnrt,CWD and seeds) for different elements
    real(r8), allocatable      :: fragmentation_scaler(:) ! scale rate of litter fragmentation based on soil layer [0-1]

    !:.........................................................................:

    ! FUELS AND FIRE
    ! fuel characteristics
    real(r8) :: sum_fuel                ! total ground fuel related to ROS (omits 1000 hr fuels) [kgC/m2]
    real(r8) :: fuel_frac(nfsc)         ! fraction of each litter class in the ros_fuel [0-1]
    real(r8) :: livegrass               ! total aboveground grass biomass in patch [kgC/m2]
    real(r8) :: fuel_bulkd              ! average fuel bulk density of the ground fuel. [kg/m3]
                                        ! (incl. live grasses, omits 1000hr fuels)
    real(r8) :: fuel_sav                ! average surface area to volume ratio of the ground fuel [cm-1]
                                        ! (incl. live grasses, omits 1000hr fuels)
    real(r8) :: fuel_mef                ! average moisture of extinction factor 
                                        ! of the ground fuel (incl. live grasses, omits 1000hr fuels)
    real(r8) :: fuel_eff_moist          ! effective avearage fuel moisture content of the ground fuel 
                                        ! (incl. live grasses. omits 1000hr fuels)
    real(r8) :: litter_moisture(nfsc)   ! moisture of litter [m3/m3]

    ! fire spread
    real(r8) :: ros_front               ! rate of forward  spread of fire [m/min]
    real(r8) :: ros_back                ! rate of backward spread of fire [m/min]
    real(r8) :: effect_wspeed           ! windspeed modified by fraction of relative grass and tree cover [m/min]
    real(r8) :: tau_l                   ! duration of lethal heating [min]
    real(r8) :: fi                      ! average fire intensity of flaming front [kJ/m/s] or [kW/m]
    integer  :: fire                    ! is there a fire? [1=yes; 0=no]
    real(r8) :: fd                      ! fire duration [min]

    ! fire effects      
    real(r8) :: scorch_ht(maxpft)       ! scorch height [m] 
    real(r8) :: frac_burnt              ! fraction burnt [0-1/day]  
    real(r8) :: tfc_ros                 ! total intensity-relevant fuel consumed - no trunks [kgC/m2 of burned ground/day]
    real(r8) :: burnt_frac_litter(nfsc) ! fraction of each litter pool burned, conditional on it being burned [0-1]

  end type fates_patch_type

end module FatesPatchMod