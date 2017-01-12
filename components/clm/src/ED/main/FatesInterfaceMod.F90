module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------
   ! Used CLM Modules
   ! INTERF-TODO:  NO CLM MODULES SHOULD BE ACCESSIBLE BY THE FATES
   ! PUBLIC API!!!!
   ! ------------------------------------------------------------------------------------

   use EDtypesMod            , only : ed_site_type
   use EDtypesMod            , only : maxPatchesPerCol
   use EDtypesMod            , only : cp_nclmax
   use EDtypesMod            , only : cp_numSWb
   use EDtypesMod            , only : cp_numlevgrnd
   use EDtypesMod            , only : cp_maxSWb
   use EDtypesMod            , only : cp_numlevdecomp
   use EDtypesMod            , only : cp_numlevdecomp_full
   use EDtypesMod            , only : cp_hlm_name
   use EDtypesMod            , only : cp_hio_ignore_val
   use EDtypesMod            , only : cp_numlevsoil
   use EDtypesMod            , only : cp_masterproc
   use FatesConstantsMod     , only : r8 => fates_r8

   implicit none

   ! ------------------------------------------------------------------------------------
   ! Notes on types
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! Naming conventions:   _gl  means ground layer dimensions
   !                       _si  means site dimensions (scalar in that case)
   !                       _pa  means patch dimensions
   !                       _rb  means radiation band
   ! ------------------------------------------------------------------------------------
   
   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches

      ! Timing Variables
      integer  :: current_year    ! Current year
      integer  :: current_month   ! month of year
      integer  :: current_day     ! day of month
      integer  :: current_tod     ! time of day (seconds past 0Z)
      integer  :: current_date    ! time of day (seconds past 0Z)
      integer  :: reference_date  ! YYYYMMDD
      real(r8) :: model_day       ! elapsed days between current date and reference
                                  ! uses ESMF functions, so prefered to pass it in as
                                  ! argument rather than calculate directly

      ! Vegetation Dynamics
      ! ---------------------------------------------------------------------------------

      ! The site level 24 hour vegetation temperature is used for various purposes during vegetation 
      ! dynamics.  However, we are currently using the bare ground patch's value [K]
      ! TO-DO: Get some consensus on the correct vegetation temperature used for phenology.
      ! It is possible that the bare-ground value is where the average is being stored.
      ! (RGK-01-2017)
      real(r8)             :: t_veg24_si

      ! Patch 24 hour vegetation temperature [K]
      real(r8),allocatable :: t_veg24_pa(:)  
      
      ! NOTE: h2osoi_vol_si is used to update surface water memory
      ! CLM/ALM may be using "waterstate%h2osoi_vol_col" on the first index (coli,1)
      ! to inform this. I think this should be re-evaluated (RGK 01/2017)
      ! Site volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
      real(r8) :: h2osoi_vol_si 

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

      ! Hydrology variables for BTRAN
      ! ---------------------------------------------------------------------------------

      ! Soil suction potential of layers in each site, negative, [mm]
      real(r8), allocatable :: smp_gl(:)

      ! Effective porosity = porosity - vol_ic, of layers in each site [-]
      real(r8), allocatable :: eff_porosity_gl(:)

      ! volumetric soil water at saturation (porosity)
      real(r8), allocatable :: watsat_gl(:)

      ! Temperature of ground layers [K]
      real(r8), allocatable :: tempk_gl(:)

      ! Liquid volume in ground layer
      real(r8), allocatable :: h2o_liqvol_gl(:)

      ! Site level filter for uptake response functions
      logical               :: filter_btran

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
      real(r8), allocatable :: t_soisno_gl(:)

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

      ! Ground Layer Structure
      real(r8),allocatable :: depth_gl(:)      ! Depth in vertical direction of ground layers
                                   ! Interface level below a "z" level (m) (1:cp_nlevgrnd) 

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
      logical, allocatable :: active_suction_gl(:)

      ! Effective fraction of roots in each soil layer 
      real(r8), allocatable :: rootr_pagl(:,:)

      ! Integrated (vertically) transpiration wetness factor (0 to 1) 
      ! (diagnostic, should not be used by HLM)
      real(r8), allocatable :: btran_pa(:)

      ! Sunlit canopy resistance [s/m]
      real(r8), allocatable :: rssun_pa(:)

      ! Shaded canopy resistance [s/m]
      real(r8), allocatable :: rssha_pa(:)

      ! Canopy conductance [mmol m-2 s-1]
      real(r8), allocatable :: gccanopy_pa(:)

      ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s)
      real(r8), allocatable :: psncanopy_pa(:)

      ! patch sunlit leaf maintenance respiration rate (umol CO2/m**2/s) 
      real(r8), allocatable :: lmrcanopy_pa(:)

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


      ! litterfall fluxes of C from FATES patches to BGC columns

      ! total labile    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lab_c_col(:)      

      !total cellulose litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_cel_c_col(:)      
      
      !total lignin    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lig_c_col(:)      

      ! Canopy Structure

      real(r8), allocatable :: elai_pa(:)  ! exposed leaf area index
      real(r8), allocatable :: esai_pa(:)  ! exposed stem area index
      real(r8), allocatable :: tlai_pa(:)  ! total leaf area index
      real(r8), allocatable :: tsai_pa(:)  ! total stem area index
      real(r8), allocatable :: htop_pa(:)  ! top of the canopy [m]
      real(r8), allocatable :: hbot_pa(:)  ! bottom of canopy? [m]

      real(r8), allocatable :: canopy_fraction_pa(:) ! Area fraction of each patch in the site
                                                     ! Use most likely for weighting
                                                     ! This is currently the projected canopy
                                                     ! area of each patch [0-1]

      real(r8), allocatable :: frac_veg_nosno_alb_pa(:) ! This is not really a fraction
                                                        ! this is actually binary based on if any
                                                        ! vegetation in the patch is exposed.
                                                        ! [0,1]

      
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

   contains
      
      procedure, public :: zero_bcs

   end type fates_interface_type

   public :: FatesInterfaceInit
   public :: set_fates_ctrlparms


contains

   ! ====================================================================================
  subroutine FatesInterfaceInit(log_unit, global_verbose)

    use FatesGlobals, only : FatesGlobalsInit

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    call FatesGlobalsInit(log_unit, global_verbose)

  end subroutine FatesInterfaceInit

   ! ====================================================================================

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
      ! Incrementally walk through linked list and deallocate
      
      
      
      ! Deallocate the site list
!      deallocate (this%sites)
      
      return
   end subroutine fates_clean


   ! ====================================================================================
   

   subroutine allocate_bcin(bc_in)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      
      ! Allocate input boundaries
      
      ! Vegetation Dynamics
      allocate(bc_in%t_veg24_pa(numPatchesPerCol))

      allocate(bc_in%wind24_pa(numPatchesPerCol))
      allocate(bc_in%relhumid24_pa(numPatchesPerCol))
      allocate(bc_in%precip24_pa(numPatchesPerCol))
      
      ! Radiation
      allocate(bc_in%solad_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_in%solai_parb(maxPatchesPerCol,cp_numSWb))
      
      ! Hydrology
      allocate(bc_in%smp_gl(cp_numlevgrnd))
      allocate(bc_in%eff_porosity_gl(cp_numlevgrnd))
      allocate(bc_in%watsat_gl(cp_numlevgrnd))
      allocate(bc_in%tempk_gl(cp_numlevgrnd))
      allocate(bc_in%h2o_liqvol_gl(cp_numlevgrnd))

      ! Photosynthesis
      allocate(bc_in%filter_photo_pa(maxPatchesPerCol))
      allocate(bc_in%dayl_factor_pa(maxPatchesPerCol))
      allocate(bc_in%esat_tv_pa(maxPatchesPerCol))
      allocate(bc_in%eair_pa(maxPatchesPerCol))
      allocate(bc_in%oair_pa(maxPatchesPerCol))
      allocate(bc_in%cair_pa(maxPatchesPerCol))
      allocate(bc_in%rb_pa(maxPatchesPerCol))
      allocate(bc_in%t_veg_pa(maxPatchesPerCol))
      allocate(bc_in%tgcm_pa(maxPatchesPerCol))
      allocate(bc_in%t_soisno_gl(cp_numlevgrnd))

      ! Canopy Radiation
      allocate(bc_in%filter_vegzen_pa(maxPatchesPerCol))
      allocate(bc_in%coszen_pa(maxPatchesPerCol))
      allocate(bc_in%albgr_dir_rb(cp_numSWb))
      allocate(bc_in%albgr_dif_rb(cp_numSWb))

      ! Carbon Balance Checking
      ! (snow-depth and snow fraction are site level and not vectors)
      
      ! Ground layer structure
      allocate(bc_in%depth_gl(0:cp_numlevgrnd))

      return
   end subroutine allocate_bcin
   
   subroutine allocate_bcout(bc_out)

      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_out_type), intent(inout) :: bc_out
      
      
      ! Radiation
      allocate(bc_out%fsun_pa(maxPatchesPerCol))
      allocate(bc_out%laisun_pa(maxPatchesPerCol))
      allocate(bc_out%laisha_pa(maxPatchesPerCol))
      
      ! Hydrology
      allocate(bc_out%active_suction_gl(cp_numlevgrnd))
      allocate(bc_out%rootr_pagl(maxPatchesPerCol,cp_numlevgrnd))
      allocate(bc_out%btran_pa(maxPatchesPerCol))
      
      ! Photosynthesis
      allocate(bc_out%rssun_pa(maxPatchesPerCol))
      allocate(bc_out%rssha_pa(maxPatchesPerCol))
      allocate(bc_out%gccanopy_pa(maxPatchesPerCol))
      allocate(bc_out%lmrcanopy_pa(maxPatchesPerCol))
      allocate(bc_out%psncanopy_pa(maxPatchesPerCol))
      
      ! Canopy Radiation
      allocate(bc_out%albd_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%albi_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%fabd_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%fabi_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftdd_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftid_parb(maxPatchesPerCol,cp_numSWb))
      allocate(bc_out%ftii_parb(maxPatchesPerCol,cp_numSWb))

      ! biogeochemistry
      allocate(bc_out%FATES_c_to_litr_lab_c_col(cp_numlevdecomp_full))        
      allocate(bc_out%FATES_c_to_litr_cel_c_col(cp_numlevdecomp_full))
      allocate(bc_out%FATES_c_to_litr_lig_c_col(cp_numlevdecomp_full))

      ! Canopy Structure
      allocate(bc_out%elai_pa(maxPatchesPerCol))
      allocate(bc_out%esai_pa(maxPatchesPerCol))
      allocate(bc_out%tlai_pa(maxPatchesPerCol))
      allocate(bc_out%tsai_pa(maxPatchesPerCol))
      allocate(bc_out%htop_pa(maxPatchesPerCol))
      allocate(bc_out%hbot_pa(maxPatchesPerCol))
      allocate(bc_out%canopy_fraction_pa(maxPatchesPerCol))
      allocate(bc_out%frac_veg_nosno_alb_pa(maxPatchesPerCol))


      return
   end subroutine allocate_bcout

   ! ====================================================================================

   subroutine zero_bcs(this,s)

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s

      ! Input boundaries
      this%bc_in(s)%current_year   = 0
      this%bc_in(s)%current_month  = 0 
      this%bc_in(s)%current_day    = 0 
      this%bc_in(s)%current_tod    = 0 
      this%bc_in(s)%current_date   = 0
      this%bc_in(s)%reference_date = 0 
      this%bc_in(s)%model_day      = 0.0_r8
      this%bc_in(s)%t_veg24_si     = 0.0_r8
      this%bc_in(s)%t_veg24_pa(:)  = 0.0_r8
      this%bc_in(s)%h2osoi_vol_si  = 0.0_r8
      this%bc_in(s)%precip24_pa(:) = 0.0_r8
      this%bc_in(s)%relhumid24_pa(:) = 0.0_r8
      this%bc_in(s)%wind24_pa(:)     = 0.0_r8

      this%bc_in(s)%solad_parb(:,:)     = 0.0_r8
      this%bc_in(s)%solai_parb(:,:)     = 0.0_r8
      this%bc_in(s)%smp_gl(:)           = 0.0_r8
      this%bc_in(s)%eff_porosity_gl(:)  = 0.0_r8
      this%bc_in(s)%watsat_gl(:)        = 0.0_r8
      this%bc_in(s)%tempk_gl(:)         = 0.0_r8
      this%bc_in(s)%h2o_liqvol_gl(:)    = 0.0_r8
      this%bc_in(s)%filter_vegzen_pa(:) = .false.
      this%bc_in(s)%coszen_pa(:)        = 0.0_r8
      this%bc_in(s)%albgr_dir_rb(:)     = 0.0_r8
      this%bc_in(s)%albgr_dif_rb(:)     = 0.0_r8
      this%bc_in(s)%max_rooting_depth_index_col = 0
      this%bc_in(s)%tot_het_resp        = 0.0_r8
      this%bc_in(s)%tot_somc            = 0.0_r8 
      this%bc_in(s)%tot_litc            = 0.0_r8
      this%bc_in(s)%snow_depth_si       = 0.0_r8
      this%bc_in(s)%frac_sno_eff_si     = 0.0_r8
      this%bc_in(s)%depth_gl(:)         = 0.0_r8
      
      ! Output boundaries
      this%bc_out(s)%active_suction_gl(:) = .false.
      this%bc_out(s)%fsun_pa(:)      = 0.0_r8
      this%bc_out(s)%laisun_pa(:)    = 0.0_r8
      this%bc_out(s)%laisha_pa(:)    = 0.0_r8
      this%bc_out(s)%rootr_pagl(:,:) = 0.0_r8
      this%bc_out(s)%btran_pa(:)     = 0.0_r8

      this%bc_out(s)%FATES_c_to_litr_lab_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_cel_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_lig_c_col(:) = 0.0_r8

      this%bc_out(s)%rssun_pa(:)     = 0.0_r8
      this%bc_out(s)%rssha_pa(:)     = 0.0_r8
      this%bc_out(s)%gccanopy_pa(:)  = 0.0_r8
      this%bc_out(s)%psncanopy_pa(:) = 0.0_r8
      this%bc_out(s)%lmrcanopy_pa(:) = 0.0_r8

      this%bc_out(s)%albd_parb(:,:) = 0.0_r8
      this%bc_out(s)%albi_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabd_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabi_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftdd_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftid_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftii_parb(:,:) = 0.0_r8

      this%bc_out(s)%elai_pa(:) = 0.0_r8
      this%bc_out(s)%esai_pa(:) = 0.0_r8
      this%bc_out(s)%tlai_pa(:) = 0.0_r8
      this%bc_out(s)%tsai_pa(:) = 0.0_r8
      this%bc_out(s)%htop_pa(:) = 0.0_r8
      this%bc_out(s)%hbot_pa(:) = 0.0_r8
      this%bc_out(s)%canopy_fraction_pa(:) = 0.0_r8
      this%bc_out(s)%frac_veg_nosno_alb_pa(:) = 0.0_r8
      
      return
    end subroutine zero_bcs
   
   ! ==================================================================================== 

   subroutine set_fates_ctrlparms(tag,ival,rval,cval)
      
      ! ---------------------------------------------------------------------------------
      ! INTERF-TODO:  NEED ALLOWANCES FOR REAL AND CHARACTER ARGS..
      ! Certain model control parameters and dimensions used by FATES are dictated by 
      ! the the driver or the host mode. To see which parameters should be filled here
      ! please also look at the ctrl_parms_type in FATESTYpeMod, in the section listing
      ! components dictated by the host model.
      !
      ! Some important points:
      ! 1. Calls to this function are likely from the clm_fates module in the HLM.
      ! 2. The calls should be preceeded by a flush function.
      ! 3. All values in ctrl_parm (FATESTypesMod.F90) that are classified as 
      !    'dictated by the HLM' must be listed in this subroutine
      ! 4. Should look like this:
      ! 
      ! call set_fates_ctrlparms('flush_to_unset')
      ! call set_fates_ctrlparms('num_sw_bbands',numrad)  ! or other variable
      ! ...
      ! call set_fates_ctrlparms('num_lev_ground',nlevgrnd)   ! or other variable
      ! call set_fates_ctrlparms('check_allset') 
      !
      ! RGK-2016
      ! ---------------------------------------------------------------------------------

      use FatesGlobals, only : fates_log, fates_global_verbose

      ! Arguments
      integer, optional, intent(in)         :: ival
      real(r8), optional, intent(in)        :: rval
      character(len=*),optional, intent(in) :: cval
      character(len=*),intent(in)           :: tag
      
      ! local variables
      logical              :: all_set
      integer,  parameter  :: unset_int = -999
      real(r8), parameter  :: unset_double = -999.9
      
      
      select case (trim(tag))
      case('flush_to_unset')
         if (fates_global_verbose()) then
            write(fates_log(), *) 'Flushing FATES control parameters prior to transfer from host'
         end if
         cp_numSwb     = unset_int
         cp_numlevgrnd = unset_int
         cp_numlevsoil = unset_int
         cp_numlevdecomp_full = unset_int
         cp_numlevdecomp      = unset_int
         cp_hlm_name          = 'unset'
         cp_hio_ignore_val    = unset_double
         cp_masterproc        = unset_int

      case('check_allset')
         
         if(cp_numSWb .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_masterproc .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES parameter unset: cp_masterproc'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numSWb > cp_maxSWb) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES sets a maximum number of shortwave bands'
               write(fates_log(), *) 'for some scratch-space, cp_maxSWb'
               write(fates_log(), *) 'it defaults to 2, but can be increased as needed'
               write(fates_log(), *) 'your driver or host model is intending to drive'
               write(fates_log(), *) 'FATES with:',cp_numSWb,' bands.'
               write(fates_log(), *) 'please increase cp_maxSWb in EDTypes to match'
               write(fates_log(), *) 'or exceed this value'
            end if
            ! end_run('MESSAGE')
         end if

         if(cp_numlevgrnd .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numlevsoil .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numlevdecomp_full .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp_full'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(cp_numlevdecomp .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(trim(cp_hlm_name) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_name'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if( abs(cp_hio_ignore_val-unset_double)<1e-10 ) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hio_ignore'
            end if
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if (fates_global_verbose()) then
            write(fates_log(), *) 'Checked. All control parameters sent to FATES.'
         end if

         
      case default

         if(present(ival))then
            select case (trim(tag))

            case('masterproc')
               cp_masterproc = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering masterproc = ',ival,' to FATES'
               end if

            case('num_sw_bbands')
               cp_numSwb = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_sw_bbands = ',ival,' to FATES'
               end if
               
            case('num_lev_ground')
               cp_numlevgrnd = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_ground = ',ival,' to FATES'
               end if

            case('num_lev_soil')
               cp_numlevsoil = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_ground = ',ival,' to FATES'
               end if

            case('num_levdecomp_full')
               cp_numlevdecomp_full = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_levdecomp_full = ',ival,' to FATES'
               end if
            
            case('num_levdecomp')
               cp_numlevdecomp = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_levdecomp = ',ival,' to FATES'
               end if

            case default
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
            
         end if
         
         if(present(rval))then
            select case (trim(tag))
            case ('hio_ignore_val')
               cp_hio_ignore_val = rval
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hio_ignore_val = ',rval,' to FATES'
               end if
            case default
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
         end if

         if(present(cval))then
            select case (trim(tag))
               
            case('hlm_name')
               cp_hlm_name = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the HLM name = ',trim(cval)
               end if

            case default
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
         end if

      end select
            
      return
   end subroutine set_fates_ctrlparms



end module FatesInterfaceMod
