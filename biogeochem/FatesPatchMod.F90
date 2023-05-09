module FatesPatchMod

  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : fates_unset_r8
  use FatesConstantsMod,   only : fates_unset_int
  use FatesConstantsMod,   only : primaryforest, secondaryforest
  use FatesGlobals,        only : fates_log
  use FatesGlobals,        only : endrun => fates_endrun
  use FatesUtilsMod,       only : check_hlm_list
  use FatesUtilsMod,       only : check_var_real
  use FatesCohortMod,      only : fates_cohort_type
  use FatesRunningMeanMod, only : rmean_type
  use FatesLitterMod,      only : nfsc
  use FatesLitterMod,      only : litter_type
  use PRTGenericMod,       only : num_elements
  use PRTGenericMod,       only : element_list
  use EDParamsMod,         only : maxSWb, nlevleaf, nclmax
  use FatesConstantsMod,   only : n_dbh_bins, maxpft, n_dist_types
  use FatesConstantsMod,   only : n_rad_stream_types
  use FatesConstantsMod,   only : t_water_freeze_k_1atm
  use FatesRunningMeanMod, only : ema_24hr, fixed_24hr, ema_lpa, ema_longterm

  use shr_infnan_mod,      only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod,         only : errMsg => shr_log_errMsg

  implicit none
  private

  ! for error message writing
  character(len=*), parameter :: sourcefile = __FILE__

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
    !:.........................................................................:

    real(r8) :: layer_height_profile(nclmax,maxpft,nlevleaf)
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
    real(r8)              :: sum_fuel                ! total ground fuel related to ROS (omits 1000 hr fuels) [kgC/m2]
    real(r8)              :: fuel_frac(nfsc)         ! fraction of each litter class in the ros_fuel [0-1]
    real(r8)              :: livegrass               ! total aboveground grass biomass in patch [kgC/m2]
    real(r8)              :: fuel_bulkd              ! average fuel bulk density of the ground fuel. [kg/m3]
                                                       ! (incl. live grasses, omits 1000hr fuels)
    real(r8)              :: fuel_sav                ! average surface area to volume ratio of the ground fuel [cm-1]
                                                       ! (incl. live grasses, omits 1000hr fuels)
    real(r8)              :: fuel_mef                ! average moisture of extinction factor 
                                                       ! of the ground fuel (incl. live grasses, omits 1000hr fuels)
    real(r8)              :: fuel_eff_moist          ! effective avearage fuel moisture content of the ground fuel 
                                                       ! (incl. live grasses. omits 1000hr fuels)
    real(r8)              :: litter_moisture(nfsc)   ! moisture of litter [m3/m3]

    ! fire spread
    real(r8)              :: ros_front               ! rate of forward  spread of fire [m/min]
    real(r8)              :: ros_back                ! rate of backward spread of fire [m/min]
    real(r8)              :: effect_wspeed           ! windspeed modified by fraction of relative grass and tree cover [m/min]
    real(r8)              :: tau_l                   ! duration of lethal heating [min]
    real(r8)              :: fi                      ! average fire intensity of flaming front [kJ/m/s] or [kW/m]
    integer               :: fire                    ! is there a fire? [1=yes; 0=no]
    real(r8)              :: fd                      ! fire duration [min]

    ! fire effects      
    real(r8)              :: scorch_ht(maxpft)       ! scorch height [m] 
    real(r8)              :: frac_burnt              ! fraction burnt [0-1/day]  
    real(r8)              :: tfc_ros                 ! total intensity-relevant fuel consumed - no trunks [kgC/m2 of burned ground/day]
    real(r8)              :: burnt_frac_litter(nfsc) ! fraction of each litter pool burned, conditional on it being burned [0-1]

    !:.........................................................................:
    
    ! PLANT HYDRAULICS (not currently used in hydraulics RGK 03-2018)  
    ! type(ed_patch_hydr_type), pointer :: pa_hydr ! All patch hydraulics data, see FatesHydraulicsMemMod.F90
    
    contains

      procedure :: init
      procedure :: nan_values 
      procedure :: zero_values
      procedure :: init_running_means
      procedure :: init_litter
      procedure :: create
      procedure :: free_memory
      procedure :: dump
      procedure :: check_vars

  end type fates_patch_type

  contains 

    subroutine init(this, numSWb, nlevsoil)
      !
      !  DESCRIPTION:
      !  Initialize a new patch - allocate arrays and set values to nan and/or 
      !     0.0
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this        ! patch object
      integer,                 intent(in)    :: numSWb      ! number of shortwave broad-bands to track
      integer,                 intent(in)    :: nlevsoil    ! number of soil layers

      ! allocate arrays 
      allocate(this%tr_soil_dir(numSWb))
      allocate(this%tr_soil_dif(numSWb))
      allocate(this%tr_soil_dir_dif(numSWb))
      allocate(this%fab(numSWb))
      allocate(this%fabd(numSWb))
      allocate(this%fabi(numSWb))
      allocate(this%sabs_dir(numSWb))
      allocate(this%sabs_dif(numSWb))
      allocate(this%fragmentation_scaler(nlevsoil))

      ! initialize all values to nan
      call this%nan_values()

      ! zero values that should be zeroed
      call this%zero_values()

    end subroutine init 

    !:.........................................................................:

    subroutine nanValues(this)
      !
      !  DESCRIPTION:
      !  Sets all values in patch to nan
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this ! patch object

      ! set pointers to null
      this%tallest  => null()   
      this%shortest => null()  
      this%older    => null()  
      this%younger  => null()
      nullify(this%tallest)
      nullify(this%shortest)
      nullify(this%older)
      nullify(this%younger)

      ! INDICES
      this%patchno                      = fates_unset_int
      this%nocomp_pft_label             = fates_unset_int                  
  
      ! PATCH INFO
      this%age                          = nan                          
      this%age_class                    = fates_unset_int
      this%area                         = nan    
      this%countcohorts                 = fates_unset_int 
      this%ncl_p                        = fates_unset_int
      this%anthro_disturbance_label     = fates_unset_int
      this%age_since_anthro_disturbance = nan
      
      ! LEAF ORGANIZATION
      this%pft_agb_profile(:,:)         = nan
      this%canopy_layer_tlai(:)         = nan               
      this%total_canopy_area            = nan
      this%total_tree_area              = nan 
      this%zstar                        = nan 
      this%elai_profile(:,:,:)          = nan 
      this%esai_profile(:,:,:)          = nan   
      this%tlai_profile(:,:,:)          = nan 
      this%tsai_profile(:,:,:)          = nan 
      this%canopy_area_profile(:,:,:)   = nan  
      this%canopy_mask(:,:)             = fates_unset_int
      this%nrad(:,:)                    = fates_unset_int
      this%ncan(:,:)                    = fates_unset_int
      this%c_stomata                    = nan 
      this%c_lblayer                    = nan
      this%layer_height_profile(:,:,:)  = nan
      
      this%psn_z(:,:,:)                 = nan 
      this%nrmlzd_parprof_pft_dir_z(:,:,:,:) = nan
      this%nrmlzd_parprof_pft_dif_z(:,:,:,:) = nan
      this%nrmlzd_parprof_dir_z(:,:,:)  = nan
      this%nrmlzd_parprof_dir_z(:,:,:) = nan

      ! RADIATION
      this%radiation_error              = nan 
      this%fcansno                      = nan 
      this%solar_zenith_flag            = .false. 
      this%solar_zenith_angle           = nan 
      this%gnd_alb_dif(:)               = nan 
      this%gnd_alb_dir(:)               = nan
      this%fabd_sun_z(:,:,:)            = nan 
      this%fabd_sha_z(:,:,:)            = nan 
      this%fabi_sun_z(:,:,:)            = nan 
      this%fabi_sha_z(:,:,:)            = nan  
      this%ed_laisun_z(:,:,:)           = nan 
      this%ed_laisha_z(:,:,:)           = nan 
      this%ed_parsun_z(:,:,:)           = nan 
      this%ed_parsha_z(:,:,:)           = nan 
      this%f_sun(:,:,:)                 = nan
      this%parprof_pft_dir_z(:,:,:)     = nan 
      this%parprof_pft_dif_z(:,:,:)     = nan
      this%parprof_dir_z(:,:)           = nan
      this%parprof_dif_z(:,:)           = nan
      this%tr_soil_dir(:)               = nan    
      this%tr_soil_dif(:)               = nan    
      this%tr_soil_dir_dif(:)           = nan
      this%fab(:)                       = nan    
      this%fabd(:)                      = nan    
      this%fabi(:)                      = nan
      this%sabs_dir(:)                  = nan 
      this%sabs_dif(:)                  = nan 
      
      ! ROOTS
      this%btran_ft(:)                  = nan 
      this%bstress_sal_ft(:)            = nan 

      ! EXTERNAL SEED RAIN 
      this%nitr_repro_stoich(:)         = nan 
      this%phos_repro_stoich(:)         = nan 
  
      ! DISTURBANCE 
      this%disturbance_rates(:)         = nan
      this%fract_ldist_not_harvested    = nan 

      ! LITTER AND COARSE WOODY DEBRIS
      this%fragmentation_scaler(:)      = nan 
  
      ! FUELS AND FIRE
      this%sum_fuel                     = nan 
      this%fuel_frac(:)                 = nan 
      this%livegrass                    = nan 
      this%fuel_bulkd                   = nan 
      this%fuel_sav                     = nan
      this%fuel_mef                     = nan 
      this%fuel_eff_moist               = nan 
      this%litter_moisture(:)           = nan
      this%ros_front                    = nan
      this%ros_back                     = nan   
      this%effect_wspeed                = nan    
      this%tau_l                        = nan
      this%fi                           = nan 
      this%fire                         = fates_unset_int
      this%fd                           = nan 
      this%scorch_ht(:)                 = nan 
      this%frac_burnt                   = nan
      this%tfc_ros                      = nan    
      this%burnt_frac_litter(:)         = nan    
  
    end subroutine nan_values

    !:.........................................................................:

    subroutine zero_values(this)
      !
      ! DESCRIPTION:
      !  sets specific variables in patch to zero
      !  these should only be values that are incremented, so that we can
      !  catch all other uninitialized variables with nans

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this
          
      ! LEAF ORGANIZATION
      this%canopy_layer_tlai(:)              = 0.0_r8
      this%total_tree_area                   = 0.0_r8  
      this%zstar                             = 0.0_r8
      this%elai_profile(:,:,:)               = 0.0_r8
      this%c_stomata                         = 0.0_r8 
      this%c_lblayer                         = 0.0_r8
      this%psn_z(:,:,:)                      = 0.0_r8
      this%nrmlzd_parprof_pft_dir_z(:,:,:,:) = 0.0_r8
      this%nrmlzd_parprof_pft_dif_z(:,:,:,:) = 0.0_r8
      this%nrmlzd_parprof_dir_z(:,:,:)       = 0.0_r8
      this%nrmlzd_parprof_dif_z(:,:,:)       = 0.0_r8

      ! RADIATION
      this%radiation_error                   = 0.0_r8
      this%fabd_sun_z(:,:,:)                 = 0.0_r8 
      this%fabd_sha_z(:,:,:)                 = 0.0_r8 
      this%fabi_sun_z(:,:,:)                 = 0.0_r8 
      this%fabi_sha_z(:,:,:)                 = 0.0_r8  
      this%ed_parsun_z(:,:,:)                = 0.0_r8 
      this%ed_parsha_z(:,:,:)                = 0.0_r8 
      this%ed_laisun_z(:,:,:)                = 0.0_r8
      this%ed_laisha_z(:,:,:)                = 0.0_r8 
      this%f_sun                             = 0.0_r8
      this%tr_soil_dir_dif(:)                = 0.0_r8
      this%fab(:)                            = 0.0_r8
      this%fabi(:)                           = 0.0_r8
      this%fabd(:)                           = 0.0_r8
      this%sabs_dir(:)                       = 0.0_r8
      this%sabs_dif(:)                       = 0.0_r8

      ! ROOTS
      this%btran_ft(:)                       = 0.0_r8

      ! DISTURBANCE 
      this%disturbance_rates(:)              = 0.0_r8 
      this%fract_ldist_not_harvested         = 0.0_r8

      ! LITTER AND COARSE WOODY DEBRIS
      this%fragmentation_scaler(:)           = 0.0_r8

      ! FIRE
      this%sum_fuel                          = 0.0_r8
      this%fuel_frac(:)                      = 0.0_r8
      this%livegrass                         = 0.0_r8
      this%fuel_bulkd                        = 0.0_r8
      this%fuel_sav                          = 0.0_r8
      this%fuel_mef                          = 0.0_r8
      this%fuel_eff_moist                    = 0.0_r8
      this%litter_moisture(:)                = 0.0_r8
      this%ros_front                         = 0.0_r8
      this%ros_back                          = 0.0_r8
      this%effect_wspeed                     = 0.0_r8
      this%tau_l                             = 0.0_r8
      this%fi                                = 0.0_r8
      this%fd                                = 0.0_r8
      this%scorch_ht(:)                      = 0.0_r8  
      this%frac_burnt                        = 0.0_r8  
      this%tfc_ros                           = 0.0_r8
      this%burnt_frac_litter(:)              = 0.0_r8

    end subroutine zero_values

    !:.........................................................................:

    subroutine init_running_means(this, current_tod)
      !
      ! DESCRIPTION:
      ! set initial values for patch running means
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this        ! patch object
      integer,                 intent(in)    :: current_tod ! time of day [seconds past 0Z]

      ! PARAMETERS:
      ! Until bc's are pointed to by sites give veg a default temp [K]
      real(r8), parameter :: temp_init_veg = 15._r8 + t_water_freeze_k_1atm

      allocate(this%tveg24)
      allocate(this%tveg_lpa)
      allocate(this%tveg_longterm)

      ! set initial values for running means
      call this%tveg24%InitRMean(fixed_24hr, init_value=temp_init_veg,         &
        init_offset=real(current_tod, r8))
      call this%tveg_lpa%InitRmean(ema_lpa, init_value=temp_init_veg)
      call this%tveg_longterm%InitRmean(ema_longterm, init_value=temp_init_veg)

    end subroutine init_running_means

    !:.........................................................................:

    subroutine init_litter(this, numpft, nlevsoil)
      !
      ! DESCRIPTION:
      ! set initial values for litter
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this     ! patch object
      integer,                 intent(in)    :: numpft   ! number of pfts to simulate
      integer,                 intent(in)    :: nlevsoil ! number of soil layers

      ! LOCALS:
      integer :: el ! looping index

      allocate(this%litter(num_elements))

      do el = 1, num_elements
        call this%litter(el)%InitAllocate(numpft, nlevsoil, element_list(el))
        call this%litter(el)%ZeroFlux()
        call this%litter(el)%InitConditions(init_leaf_fines=fates_unset_r8,  &
          init_root_fines=fates_unset_r8, init_ag_cwd=fates_unset_r8,        &
          init_bg_cwd=fates_unset_r8, init_seed=fates_unset_r8,              &
          init_seed_germ=fates_unset_r8)
     end do

    end subroutine init_litter

    !:.........................................................................:

    subroutine create(this, age, areap, label, nocomp_pft, numSWb, numpft,     &
      nlevsoil, current_tod) 
      !
      ! DESCRIPTION:
      ! create a new patch with input and default values
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this        ! patch object
      real(r8),                intent(in)    :: age         ! notional age of this patch in years
      real(r8),                intent(in)    :: areap       ! initial area of this patch in m2. 
      integer,                 intent(in)    :: label       ! anthropogenic disturbance label
      integer,                 intent(in)    :: nocomp_pft  ! no-competition mode pft label
      integer,                 intent(in)    :: numSWb      ! number of shortwave broad-bands to track
      integer,                 intent(in)    :: numpft      ! number of pfts to simulate
      integer,                 intent(in)    :: nlevsoil    ! number of soil layers
      integer,                 intent(in)    :: current_tod ! time of day [seconds past 0Z]
    
      ! initialize patch
      ! sets all values to nan, then some values to zero
      call this%init(numpft, numSWb, nlevsoil)

      ! initialize running means for patch
      call this%init_running_means(current_tod)
      
      ! initialize litter
      call this%init_litter(numpft, nlevsoil)
    
      ! assign known patch attributes 
      this%age       = age   
      this%age_class = 1
      this%area      = areap 

      ! assign anthropgenic disturbance category and label
      this%anthro_disturbance_label = label
      if (label .eq. secondaryforest) then
        this%age_since_anthro_disturbance = age
      else
        this%age_since_anthro_disturbance = fates_unset_r8
      endif
      this%nocomp_pft_label = nocomp_pft

      this%tr_soil_dir(:) = 1.0_r8
      this%tr_soil_dif(:) = 1.0_r8
      this%NCL_p          = 1

    end subroutine create

    !:.........................................................................:

    subroutine free_memory(this)
      !
      ! DESCRIPTION:
      ! deallocate the allocatable memory associated with this patch
      ! this DOES NOT deallocate the patch structure itself
      !
  
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this
  
      ! LOCALS:
      type(fates_cohort_type), pointer :: ccohort ! current cohort
      type(fates_cohort_type), pointer :: ncohort ! next cohort
      integer                          :: el      ! loop counter for elements
      integer                          :: istat   ! return status code
      character(len=255)               :: smsg    ! message string for deallocation errors
      
      ! first deallocate the cohorts
      ccohort => this%shortest
      do while(associated(ccohort))
        ncohort => ccohort%taller
        call ccohort%free_memory()
        deallocate(ccohort, stat=istat, errmsg=smsg)
        if (istat /= 0) then
          write(fates_log(),*) 'dealloc007: fail on deallocate(cchort):'//trim(smsg)
          call endrun(msg=errMsg(sourcefile, __LINE__))
        endif
        ccohort => ncohort
      end do
  
      ! deallocate all litter objects
      do el=1,num_elements
        call this%litter(el)%DeallocateLitt()
      end do
      deallocate(this%litter, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc008: fail on deallocate(this%litter):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      
      ! deallocate the allocatable arrays
      deallocate(this%tr_soil_dir,              & 
                 this%tr_soil_dif,              & 
                 this%tr_soil_dir_dif,          & 
                 this%fab,                      &
                 this%fabd,                     &
                 this%fabi,                     &
                 this%sabs_dir,                 &
                 this%sabs_dif,                 &
                 this%fragmentation_scaler,     &
                 stat=istat, errmsg=smsg)

      if (istat/=0) then
        write(fates_log(),*) 'dealloc009: fail on deallocate patch vectors:'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      
      ! deallocate running means
      deallocate(this%tveg24, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc010: fail on deallocate(this%tveg24):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      deallocate(this%tveg_lpa, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc011: fail on deallocate(this%tveg_lpa):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      deallocate(this%tveg_longterm, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc012: fail on deallocate(this%tveg_longterm):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      
    end subroutine free_memory

    !:.........................................................................:

    subroutine dump(this)
      !
      ! DESCRIPTION:
      ! print attributes of a patch
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(in) :: this ! patch object

      ! LOCALS:
      integer :: el  ! element loop counting index

      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) ' Dumping Patch Information              '
      write(fates_log(),*) ' (omitting arrays)                      '
      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) 'pa%patchno            = ',this%patchno
      write(fates_log(),*) 'pa%age                = ',this%age
      write(fates_log(),*) 'pa%age_class          = ',this%age_class
      write(fates_log(),*) 'pa%area               = ',this%area
      write(fates_log(),*) 'pa%countcohorts       = ',this%countcohorts
      write(fates_log(),*) 'pa%ncl_p              = ',this%ncl_p
      write(fates_log(),*) 'pa%total_canopy_area  = ',this%total_canopy_area
      write(fates_log(),*) 'pa%total_tree_area    = ',this%total_tree_area
      write(fates_log(),*) 'pa%zstar              = ',this%zstar
      write(fates_log(),*) 'pa%solar_zenith_flag  = ',this%solar_zenith_flag
      write(fates_log(),*) 'pa%solar_zenith_angle = ',this%solar_zenith_angle
      write(fates_log(),*) 'pa%gnd_alb_dif        = ',this%gnd_alb_dif(:)
      write(fates_log(),*) 'pa%gnd_alb_dir        = ',this%gnd_alb_dir(:)
      write(fates_log(),*) 'pa%c_stomata          = ',this%c_stomata
      write(fates_log(),*) 'pa%c_lblayer          = ',this%c_lblayer
      write(fates_log(),*) 'pa%disturbance_rates  = ',this%disturbance_rates(:)
      write(fates_log(),*) 'pa%anthro_disturbance_label = ',this%anthro_disturbance_label
      write(fates_log(),*) '----------------------------------------'

      do el = 1, num_elements
        write(fates_log(),*) 'element id: ',element_list(el)
        write(fates_log(),*) 'seed mass: ',sum(this%litter(el)%seed)
        write(fates_log(),*) 'seed germ mass: ',sum(this%litter(el)%seed_germ)
        write(fates_log(),*) 'leaf fines(pft): ',sum(this%litter(el)%leaf_fines)
        write(fates_log(),*) 'root fines(pft,sl): ',sum(this%litter(el)%root_fines)
        write(fates_log(),*) 'ag_cwd(c): ',sum(this%litter(el)%ag_cwd)
        write(fates_log(),*) 'bg_cwd(c,sl): ',sum(this%litter(el)%bg_cwd)
      end do

    end subroutine dump

    !:.........................................................................:

    subroutine check_vars(this, var_aliases, return_code)
      !
      ! DESCRIPTION:
      ! perform numerical checks on patch variables of interest
      ! The input string is of the form:  'VAR1_NAME:VAR2_NAME:VAR3_NAME'
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(in)  :: this        ! patch object
      character(len=*),        intent(in)  :: var_aliases
      integer,                 intent(out) :: return_code ! return 0 for all fine
                                                          ! return 1 if a nan detected
                                                          ! return 10+ if an overflow
                                                          ! return 100% if an underflow
      ! LOCALS:
      type(fates_cohort_type), pointer :: currentCohort

      ! Check through a registry of variables to check

      if (check_hlm_list(trim(var_aliases), 'co_n')) then
        currentCohort => this%shortest
        do while(associated(currentCohort))
            call check_var_real(currentCohort%n, 'cohort%n', return_code)
            if (.not.(return_code .eq. 0)) then
                call this%dump()
                call currentCohort%dump()
                return
            end if
            currentCohort => currentCohort%taller
        end do
      end if
  
      if (check_hlm_list(trim(var_aliases), 'co_dbh')) then
        currentCohort => this%shortest
        do while(associated(currentCohort))        
            call check_var_real(currentCohort%dbh, 'cohort%dbh', return_code)
            if (.not. (return_code .eq. 0)) then
              call this%dump()
              call currentCohort%dump()
              return
            end if
            currentCohort => currentCohort%taller
        end do
      end if

      if (check_hlm_list(trim(var_aliases), 'pa_area')) then
        call check_var_real(this%area, 'patch%area', return_code)
        if (.not. (return_code .eq. 0)) then
            call this%dump()
            return
        end if
      end if

end subroutine check_vars

!:.........................................................................:

end module FatesPatchMod