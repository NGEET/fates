module FatesPatchMod

  use FatesConstantsMod,      only : r8 => fates_r8
  use FatesConstantsMod,      only : fates_unset_r8
  use FatesConstantsMod,      only : fates_unset_int
  use FatesConstantsMod,      only : primaryland, secondaryland
  use FatesConstantsMod,      only : n_landuse_cats
  use FatesConstantsMod,      only : TRS_regeneration
  use FatesConstantsMod,      only : itrue, ifalse
  use FatesGlobals,           only : fates_log
  use FatesGlobals,           only : endrun => fates_endrun
  use FatesUtilsMod,          only : check_hlm_list
  use FatesUtilsMod,          only : check_var_real
  use FatesCohortMod,         only : fates_cohort_type
  use FatesRunningMeanMod,    only : rmean_type, rmean_arr_type
  use FatesLitterMod,         only : litter_type
  use FatesFuelMod,           only : fuel_type
  use PRTGenericMod,          only : num_elements
  use PRTGenericMod,          only : element_list
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : struct_organ, leaf_organ, sapw_organ
  use PRTParametersMod,       only : prt_params
  use FatesConstantsMod,      only : nocomp_bareground
  use EDParamsMod,            only : nlevleaf, nclmax, maxpft
  use FatesConstantsMod,      only : n_dbh_bins, n_dist_types
  use FatesConstantsMod,      only : t_water_freeze_k_1atm
  use FatesRunningMeanMod,    only : ema_24hr, fixed_24hr, ema_lpa, ema_longterm
  use FatesRunningMeanMod,    only : ema_sdlng_emerg_h2o, ema_sdlng_mort_par
  use FatesRunningMeanMod,    only : ema_sdlng2sap_par, ema_sdlng_mdd
  use TwoStreamMLPEMod,       only : twostream_type
  use FatesRadiationMemMod,   only : num_swb
  use FatesRadiationMemMod,   only : num_rad_stream_types
  use FatesInterfaceTypesMod, only : hlm_hio_ignore_val
  use FatesInterfaceTypesMod, only : numpft
  use shr_infnan_mod,         only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod,            only : errMsg => shr_log_errMsg

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
  
    !---------------------------------------------------------------------------

    ! INDICES
    integer  :: patchno          ! unique number given to each new patch created for tracking
    integer  :: nocomp_pft_label ! when nocomp is active, use this label for patch ID
                                 !   each patch ID corresponds to a pft number since each
                                 !   patch has only one pft.  Bareground patches are given
                                 !   a zero integer as a label. If nocomp is not active this 
                                 !   is set to unset. This is set in patch%Create as an argument
                                 !   to that procedure.

    !---------------------------------------------------------------------------

    ! PATCH INFO
    real(r8) :: age                          ! average patch age [years]                  
    integer  :: age_class                    ! age class of the patch for history binning purposes
    real(r8) :: area                         ! patch area [m2]
    integer  :: num_cohorts                  ! number of cohorts in patch
    integer  :: ncl_p                        ! number of occupied canopy layers
    integer  :: land_use_label               ! patch label for land use classification (primaryland, secondaryland, etc)
    real(r8) :: age_since_anthro_disturbance ! average age for secondary forest since last anthropogenic disturbance [years]
    logical  :: changed_landuse_this_ts      ! logical flag to track patches that have just undergone land use change [only used with nocomp and land use change]

    !---------------------------------------------------------------------------

    ! RUNNING MEANS
    !class(rmean_type),    pointer :: t2m                  ! place-holder for 2m air temperature (variable window-size)
    class(rmean_type),     pointer :: tveg24               ! 24-hour mean vegetation temperature [K]
    class(rmean_type),     pointer :: tveg_lpa             ! running mean of vegetation temperature at the
                                                           !   leaf photosynthesis acclimation timescale [K]
    class(rmean_type),     pointer :: tveg_longterm        ! long-term running mean of vegetation temperature at the
                                                           !   leaf photosynthesis acclimation timescale [K] (i.e T_home)
    class(rmean_type),     pointer :: seedling_layer_par24 ! 24-hour mean of photosynthetically active radiation at seedling layer [W/m2]
    class(rmean_arr_type), pointer :: sdlng_emerg_smp(:)   ! running mean of soil matric potential at the seedling
                                                           !   rooting depth at the H2O seedling emergence timescale (see sdlng_emerg_h2o_timescale parameter)
    class(rmean_type),     pointer :: sdlng_mort_par       ! running mean of photosythetically active radiation
                                                           ! at the seedling layer and at the par-based seedling  
                                                           ! mortality timescale (sdlng_mort_par_timescale)
    class(rmean_arr_type), pointer :: sdlng_mdd(:)         ! running mean of moisture deficit days
                                                           ! at the seedling layer and at the mdd-based seedling  
                                                           ! mortality timescale (sdlng_mdd_timescale) 
                                                           ! (sdlng2sap_par_timescale)
    class(rmean_type), pointer :: sdlng2sap_par            ! running mean of photosythetically active radiation
                                                           ! at the seedling layer and at the par-based seedling  
                                                           ! to sapling transition timescale 
                                                           ! (sdlng2sap_par_timescale)
    
    !---------------------------------------------------------------------------

    ! LEAF ORGANIZATION
    real(r8) :: pft_agb_profile(maxpft,n_dbh_bins)          ! binned aboveground biomass, for patch fusion [kgC/m2]
    real(r8) :: canopy_layer_tlai(nclmax)                   ! total leaf area index of each canopy layer [m2 veg/m2 canopy area]
                                                              !   (patch without bare ground)
                                                              !   used to determine attenuation of parameters during photosynthesis
    real(r8) :: total_canopy_area                           ! area that is covered by vegetation [m2]
    real(r8) :: total_tree_area                             ! area that is covered by woody vegetation [m2]
    real(r8) :: total_grass_area                            ! area that is covered by non-woody vegetation [m2]
    real(r8) :: zstar                                       ! height of smallest canopy tree, only meaningful in "strict PPA" mode [m]

    ! exposed leaf area in each canopy layer, pft, and leaf layer [m2 leaf/m2 contributing crown area]
    real(r8), allocatable :: elai_profile(:,:,:)  ! nclmax,maxpft,nlevleaf)
    ! exposed stem area in each canopy layer, pft, and leaf layer [m2 leaf/m2 contributing crown area]      
    real(r8), allocatable :: esai_profile(:,:,:)  ! nclmax,maxpft,nlevleaf)
    ! total leaf area (includes that which is under snow-pack) 
    real(r8), allocatable :: tlai_profile(:,:,:)  ! nclmax,maxpft,nlevleaf)
    ! total stem area (includes that which is under snow-pack)
    real(r8), allocatable :: tsai_profile(:,:,:)  ! nclmax,maxpft,nlevleaf)
    
    real(r8), allocatable :: canopy_area_profile(:,:,:) ! nclmax,maxpft,nlevleaf) ! fraction of crown area per canopy area in each layer
                                                        !   they will sum to 1.0 in the fully closed canopy layers
                                                        !   but only in leaf-layers that contain contributions
                                                        !   from all cohorts that donate to canopy_area

    integer  :: canopy_mask(nclmax,maxpft)                  ! is there any of this pft in this canopy layer?      
    integer  :: nrad(nclmax,maxpft)                         ! number of exposed vegetation layers for each canopy layer and pft
    integer  :: nleaf(nclmax,maxpft)                        ! number of total leaf layers for each canopy layer and pft
    real(r8) :: c_stomata                                   ! mean stomatal conductance of all leaves in the patch   [umol/m2/s]
    real(r8) :: c_lblayer                                   ! mean boundary layer conductance of all leaves in the patch [umol/m2/s]
    
    real(r8),allocatable :: nrmlzd_parprof_pft_dir_z(:,:,:,:) !num_rad_stream_types,nclmax,maxpft,nlevleaf)
    real(r8),allocatable :: nrmlzd_parprof_pft_dif_z(:,:,:,:) !num_rad_stream_types,nclmax,maxpft,nlevleaf)

    !---------------------------------------------------------------------------

    ! RADIATION
    real(r8) :: rad_error(num_swb)                        ! radiation consv error by band [W/m2]
    real(r8) :: fcansno                                   ! fraction of canopy covered in snow [0-1]
    logical  :: solar_zenith_flag                         ! integer flag specifying daylight (based on zenith angle)
    real(r8) :: solar_zenith_angle                        ! solar zenith angle [radians]
    real(r8) :: gnd_alb_dif(num_swb)                      ! ground albedo for diffuse rad, both bands [0-1]
    real(r8) :: gnd_alb_dir(num_swb)                      ! ground albedo for direct rad, both bands [0-1]
    
    
    ! organized by canopy layer, pft, and leaf layer
    real(r8),allocatable :: fabd_sun_z(:,:,:)    !nclmax,maxpft,nlevleaf)        ! sun fraction of direct light absorbed [0-1]
    real(r8),allocatable :: fabd_sha_z(:,:,:)    !nclmax,maxpft,nlevleaf)        ! shade fraction of direct light absorbed [0-1]
    real(r8),allocatable :: fabi_sun_z(:,:,:)    !nclmax,maxpft,nlevleaf)        ! sun fraction of indirect light absorbed [0-1]
    real(r8),allocatable :: fabi_sha_z(:,:,:)    !nclmax,maxpft,nlevleaf)        ! shade fraction of indirect light absorbed [0-1]
    real(r8),allocatable :: ed_parsun_z(:,:,:)   !nclmax,maxpft,nlevleaf)       ! PAR absorbed in the sun [W/m2]   
    real(r8),allocatable :: ed_parsha_z(:,:,:)   !nclmax,maxpft,nlevleaf)       ! PAR absorbed in the shade [W/m2]
    real(r8),allocatable :: f_sun(:,:,:)         !nclmax,maxpft,nlevleaf)             ! fraction of leaves in the sun [0-1]
    real(r8),allocatable :: ed_laisun_z(:,:,:)   !nclmax,maxpft,nlevleaf)
    real(r8),allocatable :: ed_laisha_z(:,:,:)   !nclmax,maxpft,nlevleaf)

    
    ! radiation profiles for comparison against observations
    real(r8),allocatable :: parprof_pft_dir_z(:,:,:)   !nclmax,maxpft,nlevleaf) ! direct-beam PAR profile through canopy, by canopy, PFT, leaf level [W/m2]
    real(r8),allocatable :: parprof_pft_dif_z(:,:,:)   !nclmax,maxpft,nlevleaf) ! diffuse     PAR profile through canopy, by canopy, PFT, leaf level [W/m2]
    
    real(r8), allocatable :: tr_soil_dir(:)               ! fraction of incoming direct radiation transmitted to the soil as direct, by numSWB [0-1]
    real(r8), allocatable :: tr_soil_dif(:)               ! fraction of incoming diffuse radiation that is transmitted to the soil as diffuse [0-1]
    real(r8), allocatable :: tr_soil_dir_dif(:)           ! fraction of incoming direct radiation that is transmitted to the soil as diffuse [0-1]
    real(r8), allocatable :: fab(:)                       ! fraction of incoming total   radiation that is absorbed by the canopy
    real(r8), allocatable :: fabd(:)                      ! fraction of incoming direct  radiation that is absorbed by the canopy
    real(r8), allocatable :: fabi(:)                      ! fraction of incoming diffuse radiation that is absorbed by the canopy
    real(r8), allocatable :: sabs_dir(:)                  ! fraction of incoming direct  radiation that is absorbed by the canopy
    real(r8), allocatable :: sabs_dif(:)                  ! fraction of incoming diffuse radiation that is absorbed by the canopy

    ! Twostream data structures
    type(twostream_type) :: twostr                        ! This holds all two-stream data and procedures
   
    
    !---------------------------------------------------------------------------

    ! ROOTS
    real(r8) :: btran_ft(maxpft)       ! btran calculated seperately for each PFT
    real(r8) :: bstress_sal_ft(maxpft) ! bstress from salinity calculated seperately for each PFT
    
    !---------------------------------------------------------------------------

    ! EXTERNAL SEED RAIN
    real(r8) :: nitr_repro_stoich(maxpft) ! The NC ratio of a new recruit in this patch
    real(r8) :: phos_repro_stoich(maxpft) ! The PC ratio of a new recruit in this patch

    !---------------------------------------------------------------------------
    
   ! DISTURBANCE 
    real(r8) :: disturbance_rates(n_dist_types)           ! disturbance rate [0-1/day] from 1) mortality
                                                          !                                 2) fire
                                                          !                                 3) logging mortatliy
                                                          !                                 4) land use change
    real(r8) :: landuse_transition_rates(n_landuse_cats)  ! land use tranision rate
    real(r8) :: fract_ldist_not_harvested                 ! fraction of logged area that is canopy trees that weren't harvested [0-1]

    !---------------------------------------------------------------------------

    ! LITTER AND COARSE WOODY DEBRIS
    type(litter_type), pointer :: litter(:)               ! litter (leaf,fnrt,CWD and seeds) for different elements
    type(fuel_type),   pointer :: fuel                    ! fuel class 
    real(r8), allocatable      :: fragmentation_scaler(:) ! scale rate of litter fragmentation based on soil layer [0-1]

    !---------------------------------------------------------------------------

    ! FUELS AND FIRE
    ! fuel characteristics
    real(r8)              :: livegrass               ! total aboveground grass biomass in patch [kgC/m2]

    ! fire spread
    real(r8)              :: ros_front               ! rate of forward  spread of fire [m/min]
    real(r8)              :: ros_back                ! rate of backward spread of fire [m/min]
    real(r8)              :: tau_l                   ! duration of lethal heating [min]
    real(r8)              :: fi                      ! average fire intensity of flaming front [kJ/m/s] or [kW/m]
    integer               :: fire                    ! is there a fire? [1=yes; 0=no]
    real(r8)              :: fd                      ! fire duration [min]
    real(r8)              :: frac_burnt              ! fraction of patch burnt by fire

    ! fire effects      
    real(r8)              :: scorch_ht(maxpft)       ! scorch height [m] 
    real(r8)              :: tfc_ros                 ! total intensity-relevant fuel consumed - no trunks [kgC/m2 of burned ground/day]
    !---------------------------------------------------------------------------
    
    ! PLANT HYDRAULICS (not currently used in hydraulics RGK 03-2018)  
    ! type(ed_patch_hydr_type), pointer :: pa_hydr ! All patch hydraulics data, see FatesHydraulicsMemMod.F90
    
    contains

      procedure :: Init
      procedure :: NanValues
      procedure :: NanDynamics
      procedure :: ZeroValues
      procedure :: ZeroDynamics
      procedure :: ReAllocateDynamics
      procedure :: InitRunningMeans
      procedure :: InitLitter
      procedure :: Create
      procedure :: CountCohorts
      procedure :: ValidateCohorts
      procedure :: InsertCohort
      procedure :: SortCohorts
      procedure :: UpdateTreeGrassArea
      procedure :: UpdateLiveGrass
      procedure :: FreeMemory
      procedure :: Dump
      procedure :: CheckVars

  end type fates_patch_type

  contains 

    !===========================================================================

    subroutine Init(this, num_swb, num_levsoil)
      !
      !  DESCRIPTION:
      !  Initialize a new patch - allocate arrays and set values to nan and/or 0.0
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this        ! patch object
      integer,                 intent(in)    :: num_swb     ! number of shortwave broad-bands to track
      integer,                 intent(in)    :: num_levsoil ! number of soil layers

      ! allocate arrays 
      allocate(this%tr_soil_dir(num_swb))
      allocate(this%tr_soil_dif(num_swb))
      allocate(this%tr_soil_dir_dif(num_swb))
      allocate(this%fab(num_swb))
      allocate(this%fabd(num_swb))
      allocate(this%fabi(num_swb))
      allocate(this%sabs_dir(num_swb))
      allocate(this%sabs_dif(num_swb))
      allocate(this%fragmentation_scaler(num_levsoil))

      ! initialize all values to nan
      call this%NanValues()

      ! zero values that should be zeroed
      call this%ZeroValues()

    end subroutine Init 

    !===========================================================================

    subroutine ReAllocateDynamics(this)

      ! ------------------------------------------------------------------------
      ! Perform allocations and re-allocations of potentially large patch arrays
      !
      ! Note that this routine is called at the very end of dynamics, after trimming
      ! and after canopy structure. This is important because we want to allocate
      ! the array spaces for leaf area and radiation scattering based on the most
      ! updated values.  Note, that this routine is also called before we re-initialize
      ! radiation scattering during restarts.
      ! ------------------------------------------------------------------------

      ! arguments
      class(fates_patch_type), intent(inout) :: this        ! patch object

      ! locals
      logical  :: re_allocate              ! Should we re-allocate the patch arrays?
      integer  :: prev_nveg                ! Previous number of vegetation layers
      integer  :: nveg                     ! Number of vegetation layers
      integer  :: ncan                     ! Number of canopy layers
      integer  :: prev_ncan                ! Number of canopy layers previously
                                           ! as defined in the allocation space
    
      ncan = this%ncl_p
      nveg = maxval(this%nleaf(:,:))

      ! Assume we will need to allocate, unless the
      ! arrays already are allocated and require the same size
      re_allocate = .true.

      ! If the large patch arrays are not new, deallocate them
      if(allocated(this%elai_profile)) then

         prev_ncan = ubound(this%tlai_profile,1)
         prev_nveg = ubound(this%tlai_profile,3)

         ! We re-allocate if the number of canopy layers has changed, or
         ! if the number of vegetation layers is larger than previously.
         ! However, we also re-allocate if the number of vegetation layers
         ! is not just smaller than previously allocated, but A GOOD BIT smaller
         ! than previously allocated.  Why?
         ! We do this so that we are not always re-allocating.

         if( prev_ncan .ne. ncan .or. (nveg>prev_nveg) .or. (nveg<prev_nveg-2) ) then

            deallocate(this%tlai_profile)
            deallocate(this%tsai_profile)
            deallocate(this%elai_profile)
            deallocate(this%esai_profile)
            deallocate(this%f_sun)
            deallocate(this%fabd_sun_z)
            deallocate(this%fabd_sha_z)
            deallocate(this%fabi_sun_z)
            deallocate(this%fabi_sha_z)
            deallocate(this%nrmlzd_parprof_pft_dir_z)
            deallocate(this%nrmlzd_parprof_pft_dif_z)
            deallocate(this%ed_parsun_z)
            deallocate(this%ed_parsha_z)
            deallocate(this%ed_laisun_z)
            deallocate(this%ed_laisha_z)
            deallocate(this%parprof_pft_dir_z)
            deallocate(this%parprof_pft_dif_z)
            deallocate(this%canopy_area_profile)
         else
            ! The number of canopy layers has not changed
            ! no need to deallocate or reallocate
            re_allocate = .false.
         end if

      end if

      ! Allocate dynamic patch arrays

      if(re_allocate) then

         ! Add a little bit of buffer to the nveg
         ! so it doesn't need to be reallocated as much
         nveg = nveg + 1
         
         allocate(this%tlai_profile(ncan,numpft,nveg))
         allocate(this%tsai_profile(ncan,numpft,nveg))
         allocate(this%elai_profile(ncan,numpft,nveg))
         allocate(this%esai_profile(ncan,numpft,nveg))
         allocate(this%canopy_area_profile(ncan,numpft,nveg))
         allocate(this%f_sun(ncan,numpft,nveg))
         allocate(this%fabd_sun_z(ncan,numpft,nveg))
         allocate(this%fabd_sha_z(ncan,numpft,nveg))
         allocate(this%fabi_sun_z(ncan,numpft,nveg))
         allocate(this%fabi_sha_z(ncan,numpft,nveg))
         allocate(this%nrmlzd_parprof_pft_dir_z(num_rad_stream_types,ncan,numpft,nveg))
         allocate(this%nrmlzd_parprof_pft_dif_z(num_rad_stream_types,ncan,numpft,nveg))
         allocate(this%ed_parsun_z(ncan,numpft,nveg))
         allocate(this%ed_parsha_z(ncan,numpft,nveg))
         allocate(this%ed_laisun_z(ncan,numpft,nveg))
         allocate(this%ed_laisha_z(ncan,numpft,nveg))
         allocate(this%parprof_pft_dir_z(ncan,numpft,nveg))
         allocate(this%parprof_pft_dif_z(ncan,numpft,nveg))
      end if

      return
    end subroutine ReAllocateDynamics
    
    !===========================================================================
    
    subroutine NanDynamics(this)
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this ! patch object

      this%elai_profile(:,:,:)          = nan 
      this%esai_profile(:,:,:)          = nan   
      this%tlai_profile(:,:,:)          = nan 
      this%tsai_profile(:,:,:)          = nan
      this%canopy_area_profile(:,:,:)   = nan  
      this%nrmlzd_parprof_pft_dir_z(:,:,:,:) = nan
      this%nrmlzd_parprof_pft_dif_z(:,:,:,:) = nan

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
      
    end subroutine NanDynamics

    !===========================================================================
    
    subroutine NanValues(this)
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
      this%num_cohorts                  = fates_unset_int 
      this%ncl_p                        = fates_unset_int
      this%land_use_label               = fates_unset_int
      this%age_since_anthro_disturbance = nan
      
      ! LEAF ORGANIZATION
      this%pft_agb_profile(:,:)         = nan
      this%canopy_layer_tlai(:)         = nan               
      this%total_canopy_area            = nan
      this%total_tree_area              = nan
      this%total_grass_area             = nan
      this%zstar                        = nan 

     
      this%canopy_mask(:,:)             = fates_unset_int
      this%nrad(:,:)                    = fates_unset_int
      this%nleaf(:,:)                   = fates_unset_int
      this%c_stomata                    = nan 
      this%c_lblayer                    = nan
      

      
      ! RADIATION
      this%rad_error(:)                 = nan
      this%fcansno                      = nan 
      this%solar_zenith_flag            = .false. 
      this%solar_zenith_angle           = nan 
      this%gnd_alb_dif(:)               = nan 
      this%gnd_alb_dir(:)               = nan

      

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

      ! LAND USE
      this%landuse_transition_rates(:)  = nan

      ! LITTER AND COARSE WOODY DEBRIS
      this%fragmentation_scaler(:)      = nan 
  
      ! FUELS AND FIRE
      this%livegrass                    = nan 
      this%ros_front                    = nan
      this%ros_back                     = nan   
      this%tau_l                        = nan
      this%fi                           = nan 
      this%fire                         = fates_unset_int
      this%fd                           = nan 
      this%scorch_ht(:)                 = nan 
      this%tfc_ros                      = nan
      this%frac_burnt                   = nan
      
    end subroutine NanValues

    !===========================================================================

    subroutine ZeroDynamics(this)
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this

      this%f_sun(:,:,:) = 0._r8
      this%fabd_sun_z(:,:,:) = 0._r8
      this%fabi_sun_z(:,:,:) = 0._r8
      this%fabd_sha_z(:,:,:) = 0._r8
      this%fabi_sha_z(:,:,:) = 0._r8
      this%nrmlzd_parprof_pft_dir_z(:,:,:,:) = 0._r8
      this%nrmlzd_parprof_pft_dif_z(:,:,:,:) = 0._r8

      ! Added
      this%elai_profile(:,:,:)          = 0._r8
      this%esai_profile(:,:,:)          = 0._r8
      this%tlai_profile(:,:,:)          = 0._r8
      this%tsai_profile(:,:,:)          = 0._r8
      this%canopy_area_profile(:,:,:)   = 0._r8
      
      this%ed_laisun_z(:,:,:)           = 0._r8
      this%ed_laisha_z(:,:,:)           = 0._r8
      this%ed_parsun_z(:,:,:)           = 0._r8
      this%ed_parsha_z(:,:,:)           = 0._r8
      this%parprof_pft_dir_z(:,:,:)     = 0._r8
      this%parprof_pft_dif_z(:,:,:)     = 0._r8
      
    end subroutine ZeroDynamics
    
    !===========================================================================
    
    subroutine ZeroValues(this)
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
      this%total_grass_area                  = 0.0_r8
      this%zstar                             = 0.0_r8
      
      this%c_stomata                         = 0.0_r8 
      this%c_lblayer                         = 0.0_r8

      
      ! RADIATION
      this%rad_error(:)                      = 0.0_r8
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

      ! LAND USE
      this%landuse_transition_rates(:)       = 0.0_r8

      ! LITTER AND COARSE WOODY DEBRIS
      this%fragmentation_scaler(:)           = 0.0_r8

      ! FIRE
      this%livegrass                         = 0.0_r8
      this%ros_front                         = 0.0_r8
      this%ros_back                          = 0.0_r8
      this%tau_l                             = 0.0_r8
      this%fi                                = 0.0_r8
      this%fd                                = 0.0_r8
      this%scorch_ht(:)                      = 0.0_r8  
      this%tfc_ros                           = 0.0_r8
      this%frac_burnt                        = 0.0_r8

    end subroutine ZeroValues

    !===========================================================================

    subroutine InitRunningMeans(this, current_tod, regeneration_model, numpft)
      !
      ! DESCRIPTION:
      ! set initial values for patch running means
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this               ! patch object
      integer,                 intent(in)    :: current_tod        ! time of day [seconds past 0Z]
      integer,                 intent(in)    :: regeneration_model ! regeneration model type
      integer,                 intent(in)    :: numpft             ! number of pfts on patch

      ! PARAMETERS:
      ! Until bc's are pointed to by sites give veg a default temp [K]
      real(r8), parameter :: temp_init_veg     = 15._r8 + t_water_freeze_k_1atm
      real(r8), parameter :: init_seedling_par = 5.0_r8      ! arbitrary initialization for seedling layer [MJ m-2 d-1]
      real(r8), parameter :: init_seedling_smp = -26652.0_r8 ! abitrary initialization of smp [mm]
      integer             :: pft                             ! pft looping index

      allocate(this%tveg24)
      allocate(this%tveg_lpa)
      allocate(this%tveg_longterm)

      ! set initial values for running means
      call this%tveg24%InitRMean(fixed_24hr, init_value=temp_init_veg,         &
        init_offset=real(current_tod, r8))
      call this%tveg_lpa%InitRmean(ema_lpa, init_value=temp_init_veg)
      call this%tveg_longterm%InitRmean(ema_longterm, init_value=temp_init_veg)

      if (regeneration_model == TRS_regeneration) then
        allocate(this%seedling_layer_par24)
        allocate(this%sdlng_mdd(numpft))
        allocate(this%sdlng_emerg_smp(numpft))
        allocate(this%sdlng_mort_par)
        allocate(this%sdlng2sap_par)

        call this%seedling_layer_par24%InitRMean(fixed_24hr,                   &
          init_value=init_seedling_par, init_offset=real(current_tod, r8))
        call this%sdlng_mort_par%InitRMean(ema_sdlng_mort_par,                 &
          init_value=temp_init_veg)
        call this%sdlng2sap_par%InitRMean(ema_sdlng2sap_par,                   &
          init_value=init_seedling_par)

        do pft = 1,numpft
          allocate(this%sdlng_mdd(pft)%p)
          allocate(this%sdlng_emerg_smp(pft)%p)

          call this%sdlng_mdd(pft)%p%InitRMean(ema_sdlng_mdd,             &
            init_value=0.0_r8)
          call this%sdlng_emerg_smp(pft)%p%InitRMean(ema_sdlng_emerg_h2o, &
            init_value=init_seedling_smp)
        end do
     end if

    end subroutine InitRunningMeans

    !===========================================================================

    subroutine InitLitter(this, num_pft, num_levsoil)
      !
      ! DESCRIPTION:
      ! set initial values for litter
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this        ! patch object
      integer,                 intent(in)    :: num_pft     ! number of pfts to simulate
      integer,                 intent(in)    :: num_levsoil ! number of soil layers

      ! LOCALS:
      integer :: el ! looping index

      allocate(this%litter(num_elements))

      do el = 1, num_elements
        call this%litter(el)%InitAllocate(num_pft, num_levsoil, element_list(el))
        call this%litter(el)%ZeroFlux()
        call this%litter(el)%InitConditions(init_leaf_fines=fates_unset_r8,  &
          init_root_fines=fates_unset_r8, init_ag_cwd=fates_unset_r8,        &
          init_bg_cwd=fates_unset_r8, init_seed=fates_unset_r8,              &
          init_seed_germ=fates_unset_r8)
     end do

    end subroutine InitLitter

    !===========================================================================

    subroutine Create(this, age, area, land_use_label, nocomp_pft, num_swb, num_pft,    &
      num_levsoil, current_tod, regeneration_model) 
      !
      ! DESCRIPTION:
      ! create a new patch with input and default values
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this               ! patch object
      real(r8),                intent(in)    :: age                ! notional age of this patch in years
      real(r8),                intent(in)    :: area               ! initial area of this patch in m2. 
      integer,                 intent(in)    :: land_use_label     ! land use label
      integer,                 intent(in)    :: nocomp_pft         ! no-competition mode pft label
      integer,                 intent(in)    :: num_swb            ! number of shortwave broad-bands to track
      integer,                 intent(in)    :: num_pft            ! number of pfts to simulate
      integer,                 intent(in)    :: num_levsoil        ! number of soil layers
      integer,                 intent(in)    :: current_tod        ! time of day [seconds past 0Z]
      integer,                 intent(in)    :: regeneration_model ! regeneration model version
    
      ! initialize patch
      ! sets all values to nan, then some values to zero
      call this%Init(num_swb, num_levsoil)

      ! initialize running means for patch
      call this%InitRunningMeans(current_tod, regeneration_model, num_pft)
      
      ! initialize litter
      call this%InitLitter(num_pft, num_levsoil)

      ! initialize fuel
      allocate(this%fuel)
      call this%fuel%Init()

      this%twostr%scelg => null()  ! The radiation module will check if this
                                   ! is associated, since it is not, it will then
                                   ! initialize and allocate
      
      ! assign known patch attributes 
      this%age       = age   
      this%age_class = 1
      this%area      = area 

      ! assign anthropgenic disturbance category and label
      this%land_use_label = land_use_label
      if (land_use_label .eq. secondaryland) then
        this%age_since_anthro_disturbance = age
      else
        this%age_since_anthro_disturbance = fates_unset_r8
      endif
      this%nocomp_pft_label = nocomp_pft

      this%tr_soil_dir(:) = 1.0_r8
      this%tr_soil_dif(:) = 1.0_r8
      this%NCL_p          = 1

      this%changed_landuse_this_ts = .false.

    end subroutine Create

    !===========================================================================

    subroutine UpdateTreeGrassArea(this)
      !
      ! DESCRIPTION:
      ! calculate and update the total tree area and grass area (by canopy) on patch
      !

      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this ! patch object 

      ! LOCALS:
      type(fates_cohort_type), pointer :: currentCohort ! cohort object
      real(r8)                         :: tree_area     ! treed area of patch [m2]
      real(r8)                         :: grass_area    ! grass area of patch [m2]

      if (this%nocomp_pft_label /= nocomp_bareground) then 
        tree_area = 0.0_r8
        grass_area = 0.0_r8
        
        currentCohort => this%tallest
        do while(associated(currentCohort))
          if (prt_params%woody(currentCohort%pft) == itrue) then
            tree_area = tree_area + currentCohort%c_area
          else
            grass_area = grass_area + currentCohort%c_area
          end if
          currentCohort => currentCohort%shorter
        end do
        
        this%total_tree_area = min(tree_area, this%area)
        this%total_grass_area = min(grass_area, this%area)
      end if 

    end subroutine UpdateTreeGrassArea

    !===========================================================================

    subroutine UpdateLiveGrass(this)
      !
      ! DESCRIPTION:
      ! Calculates the sum of live grass biomass [kgC/m2] on a patch
    
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this ! patch
      
      ! LOCALS:
      real(r8)                         :: live_grass    ! live grass [kgC/m2]
      type(fates_cohort_type), pointer :: currentCohort ! cohort type

      live_grass = 0.0_r8
      currentCohort => this%tallest
      do while(associated(currentCohort))
          ! for grasses sum all aboveground tissues
          if (prt_params%woody(currentCohort%pft) == ifalse) then 
            live_grass = live_grass +                                      &
              (currentCohort%prt%GetState(leaf_organ, carbon12_element) +  &
              currentCohort%prt%GetState(sapw_organ, carbon12_element) +   &
              currentCohort%prt%GetState(struct_organ, carbon12_element))* &
              currentCohort%n/this%area
        endif
        currentCohort => currentCohort%shorter
      enddo
      
      this%livegrass = live_grass

    end subroutine UpdateLiveGrass

    !===========================================================================

    subroutine FreeMemory(this, regeneration_model, numpft)
      !
      ! DESCRIPTION:
      ! deallocate the allocatable memory associated with this patch
      ! this DOES NOT deallocate the patch structure itself
      !
  
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout) :: this
      integer,                 intent(in)    :: regeneration_model
      integer,                 intent(in)    :: numpft
  
      ! LOCALS:
      type(fates_cohort_type), pointer :: ccohort ! current cohort
      type(fates_cohort_type), pointer :: ncohort ! next cohort
      integer                          :: el      ! loop counter for elements
      integer                          :: pft     ! loop counter for pfts
      integer                          :: istat   ! return status code
      character(len=255)               :: smsg    ! message string for deallocation errors
      
      ! first deallocate the cohorts
      ccohort => this%shortest
      do while(associated(ccohort))
        ncohort => ccohort%taller
        call ccohort%FreeMemory()
        deallocate(ccohort, stat=istat, errmsg=smsg)
        if (istat /= 0) then
          write(fates_log(),*) 'dealloc007: fail on deallocate(cchort):'//trim(smsg)
          call endrun(msg=errMsg(sourcefile, __LINE__))
        endif
        ccohort => ncohort
      end do

      ! Deallocate Radiation scattering elements
      if(associated(this%twostr%scelg)) then
         call this%twostr%DeallocTwoStream()
      end if
      
      ! deallocate all litter objects
      do el=1,num_elements
        call this%litter(el)%DeallocateLitt()
      end do
      deallocate(this%litter, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc008: fail on deallocate(this%litter):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif

      deallocate(this%fuel, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc009: fail on deallocate patch fuel:'//trim(smsg)
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

      ! These arrays are allocated via a call from EDCanopyStructureMod
      ! while determining how many canopy and leaf layers the patch has.
      ! Its possible that patches may be spawned and destroyed before
      ! ever reaching that routine, thus we must check to see
      ! if the they are already allocated.
      if(allocated(this%elai_profile)) then
         deallocate(this%tlai_profile)
         deallocate(this%tsai_profile)
         deallocate(this%elai_profile)
         deallocate(this%esai_profile)
         deallocate(this%f_sun)
         deallocate(this%fabd_sun_z)
         deallocate(this%fabd_sha_z)
         deallocate(this%fabi_sun_z)
         deallocate(this%fabi_sha_z)
         deallocate(this%nrmlzd_parprof_pft_dir_z)
         deallocate(this%nrmlzd_parprof_pft_dif_z)
         deallocate(this%ed_parsun_z)
         deallocate(this%ed_parsha_z)
         deallocate(this%ed_laisun_z)
         deallocate(this%ed_laisha_z)
         deallocate(this%parprof_pft_dir_z)
         deallocate(this%parprof_pft_dif_z)
         deallocate(this%canopy_area_profile)
      end if
      
      if (istat/=0) then
        write(fates_log(),*) 'dealloc010: fail on deallocate patch vectors:'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      
      ! deallocate running means
      deallocate(this%tveg24, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc011: fail on deallocate(this%tveg24):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      deallocate(this%tveg_lpa, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc012: fail on deallocate(this%tveg_lpa):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      deallocate(this%tveg_longterm, stat=istat, errmsg=smsg)
      if (istat/=0) then
        write(fates_log(),*) 'dealloc013: fail on deallocate(this%tveg_longterm):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif

      if (regeneration_model == TRS_regeneration) then 
        deallocate(this%seedling_layer_par24)
        deallocate(this%sdlng_mort_par)
        deallocate(this%sdlng2sap_par)
        do pft = 1, numpft 
          deallocate(this%sdlng_mdd(pft)%p)
        end do 
        deallocate(this%sdlng_mdd)
        do pft = 1, numpft 
          deallocate(this%sdlng_emerg_smp(pft)%p)
        end do 
        deallocate(this%sdlng_emerg_smp)
      end if 
      
    end subroutine FreeMemory

    !===========================================================================
    
    subroutine InsertCohort(this, cohort)
      !
      ! DESCRIPTION:
      ! Inserts a cohort into a patch's linked list structure
      !
      
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout), target  :: this   ! patch 
      type(fates_cohort_type), intent(inout), pointer :: cohort ! cohort to insert
      
      ! LOCALS:
      type(fates_cohort_type), pointer :: temp_cohort1, temp_cohort2 ! temporary cohorts to store pointers
      
      ! validate the cohort before insertion
      if (.not. associated(cohort)) then
        call endrun(msg="cohort is not allocated",                                       &
          additional_msg=errMsg(sourcefile, __LINE__))
          return
      end if
      
      ! check for inconsistent list state
      if ((.not. associated(this%shortest) .and. associated(this%tallest)) .or.          &
        (associated(this%shortest) .and. .not. associated(this%tallest))) then
        call endrun(msg="inconsistent list state",                                       &
          additional_msg=errMsg(sourcefile, __LINE__))
          return
      end if
    
      ! nothing in the list - add to head
      if (.not. associated(this%shortest)) then
        this%shortest => cohort
        this%tallest  => this%shortest
        cohort%taller => null()
        cohort%shorter => null()
        return
      end if 
        
      ! shortest - add to front of list
      if (cohort%height < this%shortest%height) then
        temp_cohort1 => this%shortest         ! save current shortest in temporary pointer
        cohort%taller => this%shortest        ! attach cohort to list
        this%shortest => cohort               ! cohort is now the shortest 
        this%shortest%shorter => null()       ! nullify new head's "shorter" pointer
        temp_cohort1%shorter => this%shortest ! new head is previous head's 'shorter'
        return
      end if 

      ! tallest - add to end
      if (cohort%height >= this%tallest%height) then
        this%tallest%taller => cohort        ! attach cohort to end of list
        temp_cohort1 => this%tallest         ! store current tallest in temporary pointer
        this%tallest => cohort               ! cohort is now the tallest
        this%tallest%shorter => temp_cohort1 ! new tail is previous tails's 'taller'
        this%tallest%taller => null()        ! nullify new tails's "taller" pointer
        return
      end if 

      ! traverse list to find where to put cohort
      temp_cohort1 => this%shortest
      temp_cohort2 => temp_cohort1%taller
      do while (associated(temp_cohort2))
        
        ! validate list structure before insertion
        
        if (associated(temp_cohort1%taller) .and.                                        &
          .not. associated(temp_cohort1%taller%shorter, temp_cohort1)) then 
          call endrun(msg="corrupted list structure",                                    &
            additional_msg=errMsg(sourcefile, __LINE__))
            return
        end if 
        
        if ((cohort%height >= temp_cohort1%height) .and. (cohort%height < temp_cohort2%height)) then 
          ! add cohort here
          cohort%taller => temp_cohort2
          temp_cohort1%taller => cohort
          cohort%shorter => temp_cohort1
          temp_cohort2%shorter => cohort
          exit
        end if
        temp_cohort1 => temp_cohort2
        temp_cohort2 => temp_cohort2%taller
      end do
    
    end subroutine InsertCohort
    
    !===========================================================================
    
    subroutine ValidateCohorts(this)
      !
      ! DESCRIPTION:
      ! Validates a patch's cohort linked list
      !
      
      ! ARGUMENTS:
      class(fates_patch_type), intent(in), target :: this ! patch
           
      ! LOCALS:
      type(fates_cohort_type), pointer :: currentCohort                 ! cohort object
      integer                          :: forward_count, backward_count ! forwards and backwards counts of cohorts

      ! check initial conditions
      if (.not. associated(this%shortest) .and. .not. associated(this%tallest)) then
          ! validation passed - empty list
          return
      else if (.not. associated(this%shortest) .or. .not. associated(this%tallest)) then
          call endrun(msg="one of shortest or tallest is null",                          &
            additional_msg=errMsg(sourcefile, __LINE__))
            return
      end if
      
      ! initialize counts
      forward_count = 0
      backward_count = 0
      
      ! traverse taller chain
      currentCohort => this%shortest
      do while (associated(currentCohort))
        forward_count = forward_count + 1

        ! validate cohort
        if (associated(currentCohort%taller)) then
          if (.not. associated(currentCohort%taller%shorter, currentCohort)) then
              call endrun(msg="mismatch in patch's taller chain",                        &
                additional_msg=errMsg(sourcefile, __LINE__))
                return
          end if
        else
          if (.not. associated(currentCohort, this%tallest)) then
            call endrun(msg="cohort list does not end at tallest",                       &
              additional_msg=errMsg(sourcefile, __LINE__))
              return
          end if
        end if
          currentCohort => currentCohort%taller
      end do
      
      ! traverse shorter chain
      currentCohort => this%tallest
      do while (associated(currentCohort))
        backward_count = backward_count + 1

        if (associated(currentCohort%shorter)) then
          if (.not. associated(currentCohort%shorter%taller, currentCohort)) then
            call endrun(msg="mismatch in patch's shorter chain",                         &
              additional_msg=errMsg(sourcefile, __LINE__))
              return
          end if
        else
          if (.not. associated(currentCohort, this%shortest)) then
            call endrun(msg="cohort list does not start at shortest",                    &
              additional_msg=errMsg(sourcefile, __LINE__))
              return
          end if
        end if
        currentCohort => currentCohort%shorter
      end do
      
      ! check consistency between forward and backward counts
      if (forward_count /= backward_count) then
        call endrun(msg="forward and backward traversal counts do not match",            &
          additional_msg=errMsg(sourcefile, __LINE__))
          return
      end if
        
    end subroutine ValidateCohorts
    
    !===========================================================================
    
    subroutine CountCohorts(this)
      !
      ! DESCRIPTION:
      ! Counts the number of a cohorts in a patch's linked list and updates 
      ! the this%num_cohorts attribute
      !
      
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout), target :: this ! patch 
      
      ! LOCALS:
      type(fates_cohort_type), pointer :: currentCohort ! cohort object
      integer                          :: cohort_count  ! count of cohorts 
      
      cohort_count = 0
      currentCohort => this%shortest
      do while (associated(currentCohort))
        cohort_count = cohort_count + 1
        currentCohort => currentCohort%taller
      end do
      
      this%num_cohorts = cohort_count
          
    end subroutine CountCohorts
    
    !===========================================================================
    
    subroutine SortCohorts(this)
      !
      ! DESCRIPTION: sort cohorts in patch's linked list
      ! uses insertion sort to build a new list
      !
    
      ! ARGUMENTS:
      class(fates_patch_type), intent(inout), target :: this ! patch
      
      ! LOCALS:
      type(fates_cohort_type), pointer :: currentCohort
      type(fates_cohort_type), pointer :: nextCohort
      
      ! check for inconsistent list state
      if (.not. associated(this%shortest) .and. .not. associated(this%tallest)) then
          ! empty list
          return
      else if (.not. associated(this%shortest) .or. .not. associated(this%tallest)) then
          call endrun(msg="inconsistent list state",                                     &
            additional_msg=errMsg(sourcefile, __LINE__))
          return
      end if
      
      ! hold on to current linked list so we don't lose it
      currentCohort => this%shortest
      
      ! reset the current list: we'll build it incrementally
      this%shortest => null()
      this%tallest => null()
      
      ! insert each cohort
      do while (associated(currentCohort))
        ! store the next cohort to sort
        nextCohort => currentCohort%taller
        call this%InsertCohort(currentCohort)
        currentCohort => nextCohort
      end do
    
    end subroutine SortCohorts
  
    !===========================================================================

    subroutine Dump(this)
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
      write(fates_log(),*) 'pa%num_cohorts        = ',this%num_cohorts
      write(fates_log(),*) 'pa%ncl_p              = ',this%ncl_p
      write(fates_log(),*) 'pa%total_canopy_area  = ',this%total_canopy_area
      write(fates_log(),*) 'pa%total_tree_area    = ',this%total_tree_area
      write(fates_log(),*) 'pa%total_grass_area   = ',this%total_grass_area
      write(fates_log(),*) 'pa%zstar              = ',this%zstar
      write(fates_log(),*) 'pa%solar_zenith_flag  = ',this%solar_zenith_flag
      write(fates_log(),*) 'pa%solar_zenith_angle = ',this%solar_zenith_angle
      write(fates_log(),*) 'pa%gnd_alb_dif        = ',this%gnd_alb_dif(:)
      write(fates_log(),*) 'pa%gnd_alb_dir        = ',this%gnd_alb_dir(:)
      write(fates_log(),*) 'pa%c_stomata          = ',this%c_stomata
      write(fates_log(),*) 'pa%c_lblayer          = ',this%c_lblayer
      write(fates_log(),*) 'pa%disturbance_rates  = ',this%disturbance_rates(:)
      write(fates_log(),*) 'pa%land_use_label     = ',this%land_use_label
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

    end subroutine Dump

    !===========================================================================

    subroutine CheckVars(this, var_aliases, return_code)
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
                call this%Dump()
                call currentCohort%Dump()
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
              call this%Dump()
              call currentCohort%Dump()
              return
            end if
            currentCohort => currentCohort%taller
        end do
      end if

      if (check_hlm_list(trim(var_aliases), 'pa_area')) then
        call check_var_real(this%area, 'patch%area', return_code)
        if (.not. (return_code .eq. 0)) then
            call this%Dump()
            return
        end if
      end if

    end subroutine CheckVars

    !===========================================================================  

end module FatesPatchMod
