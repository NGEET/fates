module FatesLitterMod

  ! -------------------------------------------------------------------------------------
  ! This module contains methods and type definitions for all things litter.
  ! "litter" means all organic material that is no longer associated with a live plant.
  ! Also, in FATES we only track "un-fragmented" and "un-decomposed" litter. This
  ! is a decision of pragmatism, as FATES is not a soil decomposition model, yet FATES
  ! does need to retain litter for fire calculations.  Therefore, we retain
  ! undecomposed litter for a period of time in FATES, until it fragments and is passed
  ! to another model to handle deocomposition.
  ! 
  ! This encompasses:  1) "Coarse Woody Debris"
  !                    2) fine materials leaves, roots etc 
  !                       (sometimes exclusively refered to as litter)
  !                    3) Reproductive materials (seeds, nuts, fruits)
  !
  ! Important point:   THESE POOLS DO NOT CONTAIN DECOMPOSING MATTER !!!!!
  !
  ! Another Important Point: We track the fine litter by its "decomposability" pool.
  !                          However, we don't actually apply any differential
  !                          turnover rates based on these pools, we are just
  !                          differentiating, tracking and preserving them to be 
  !                          passed in the correct partitions to the BGC model. 
  !                          Their partitions are a PFT parameter.
  ! 
  ! -------------------------------------------------------------------------------------


   ! To-do:
   ! 8) In CWD_IN, add the flux diagnostics, then remove the
   !    patch level rate in the history code

   
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : i4 => fates_int
   use FatesConstantsMod, only : nearzero
   use FatesConstantsMod, only : calloc_abs_error
   use FatesConstantsMod, only : fates_unset_r8
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log 
   use shr_log_mod      , only : errMsg => shr_log_errMsg

   implicit none
   private

   public :: adjust_SF_CWD_frac

   integer, public, parameter :: ncwd  = 4    ! number of coarse woody debris pools 
                                              ! (twig,s branch,l branch, trunk)

   integer, public, parameter :: ndcmpy = 3   ! number of "decomposability" pools in
                                              ! fines (lignin, cellulose, labile)

   integer, public, parameter :: ilabile     = 1   ! Array index for labile portion
   integer, public, parameter :: icellulose  = 2   ! Array index for cellulose portion
   integer, public, parameter :: ilignin     = 3   ! Array index for the lignin portion

    ! SPITFIRE     

  integer,  parameter, public :: NFSC                 = NCWD+2     ! number fuel size classes  (4 cwd size classes, leaf litter, and grass)
  integer,  parameter, public :: tw_sf                = 1          ! array index of twig pool for spitfire
  integer,  parameter, public :: lb_sf                = 3          ! array index of large branch pool for spitfire
  integer,  parameter, public :: tr_sf                = 4          ! array index of dead trunk pool for spitfire
  integer,  parameter, public :: dl_sf                = 5          ! array index of dead leaf pool for spitfire (dead grass and dead leaves)
  integer,  parameter, public :: lg_sf                = 6          ! array index of live grass pool for spitfire



   type, public ::  litter_type
      
      ! This object is allocated for each element (C, N, P, etc) that we wish to track.

      integer              :: element_id            ! This element ID should
                                                    ! be associated with the element
                                                    ! types listed in parteh/PRTGenericMod.F90

      ! ---------------------------------------------------------------------------------
      ! Prognostic variables (litter and coarse woody debris)
      ! Note that we do not track the fines (leaf/fine-root debris) by PFT. We track them 
      ! by their decomposing pools (i.e. chemical fraction).  This is the same dimensioning
      ! that gets passed back to the external BGC model, and saves a lot of space.
      ! ---------------------------------------------------------------------------------


      real(r8)             :: ag_cwd(ncwd)          ! above ground coarse wood debris (cwd)         [kg/m2]
      real(r8),allocatable :: bg_cwd(:,:)           ! below ground coarse wood debris (cwd x soil)  [kg/m2]
      real(r8),allocatable :: leaf_fines(:)         ! above ground leaf litter (dcmpy)              [kg/m2]
      real(r8),allocatable :: root_fines(:,:)       ! below ground fine root litter (dcmpy x soil)  [kg/m2]
      
      real(r8),allocatable :: seed(:)               ! the seed pool (viable)    (pft) [kg/m2]
      real(r8),allocatable :: seed_germ(:)          ! the germinated seed pool  (pft) [kg/m2]


      ! ---------------------------------------------------------------------------------
      ! Fluxes in - dying trees / seed rain  (does not include disturbance fluxes)
      ! ---------------------------------------------------------------------------------

      real(r8)             ::  ag_cwd_in(ncwd)      ! (cwd)        [kg/m2/day]
      real(r8),allocatable ::  bg_cwd_in(:,:)       ! (cwd x soil) [kg/m2/day]
      real(r8),allocatable ::  leaf_fines_in(:)     ! (dcmpy)       [kg/m2/day]
      real(r8),allocatable ::  root_fines_in(:,:)   ! (dcmpy x soil [kg/m2/day]

      real(r8),allocatable ::  seed_in_local(:)     ! (pft)        [kg/m2/day] (from local sources)
      real(r8),allocatable ::  seed_in_extern(:)    ! (pft)        [kg/m2/day] (from outside cell)

                                                    
      ! ---------------------------------------------------------------------------------
      ! Fluxes out - fragmentation, seed decay (does not include disturbance)
      ! ---------------------------------------------------------------------------------

      real(r8)             ::  ag_cwd_frag(ncwd)    ! above ground cwd fragmentation flux   [kg/m2/day]
      real(r8),allocatable ::  bg_cwd_frag(:,:)     ! below ground cwd fragmentation flux   [kg/m2/day]
      real(r8),allocatable ::  leaf_fines_frag(:)   ! above ground fines fragmentation flux [kg/m2/day]
      real(r8),allocatable ::  root_fines_frag(:,:) ! kg/m2/day

      real(r8), allocatable :: seed_decay(:)      ! decay of viable seeds to litter     [kg/m2/day]
      real(r8), allocatable :: seed_germ_decay(:) ! decay of germinated seeds to litter [kg/m2/day]
      real(r8), allocatable :: seed_germ_in(:)    ! flux from viable to germinated seed [kg/m2/day]

    contains
      
      procedure,non_overridable :: InitAllocate
      procedure,non_overridable :: DeallocateLitt
      procedure,non_overridable :: InitConditions
      procedure,non_overridable :: FuseLitter
      procedure,non_overridable :: CopyLitter
      procedure,non_overridable :: ZeroFlux
      procedure,non_overridable :: GetTotalLitterMass
      
   end type litter_type

   ! Part 3: Public extended types

   character(len=*), parameter, private :: sourcefile = __FILE__

contains

  subroutine FuseLitter(this,self_area,donor_area,donor_litt)

    ! -----------------------------------------------------------------------------------
    ! The litter pools are all area normalized.  This routine
    ! will use area weighting to determine the resulting
    ! litter density per area of all the pools. Essentially
    ! summing up the total mass by multiplying each component
    ! area, and then normalizing by the new total.
    ! -----------------------------------------------------------------------------------


    class(litter_type) :: this
    real(r8),intent(in)           :: self_area
    real(r8),intent(in)           :: donor_area
    type(litter_type),intent(in)  :: donor_litt

    ! locals
    integer  :: nlevsoil        ! number of soil layers
    integer  :: c               ! cwd index
    integer  :: pft             ! pft index
    integer  :: ilyr            ! soil layer index
    integer  :: dcmpy           ! dcmpyical pool index
    integer  :: npft            ! number of PFTs
    real(r8) :: self_weight     ! weighting of the recieving litter pool
    real(r8) :: donor_weight    ! weighting of the donating litter pool
    

    nlevsoil = size(this%bg_cwd,dim=2)
    npft     = size(this%seed,dim=1)

    self_weight  = self_area /(donor_area+self_area)
    donor_weight = 1._r8 - self_weight

    
    do c=1,ncwd
       this%ag_cwd(c)      = this%ag_cwd(c) *self_weight +  &
                             donor_litt%ag_cwd(c) * donor_weight
       this%ag_cwd_in(c)   = this%ag_cwd_in(c) *self_weight + &
                             donor_litt%ag_cwd_in(c) * donor_weight
       this%ag_cwd_frag(c) = this%ag_cwd_frag(c) *self_weight + &
                             donor_litt%ag_cwd_frag(c) * donor_weight
       do ilyr = 1,nlevsoil
          this%bg_cwd(c,ilyr)      = this%bg_cwd(c,ilyr) * self_weight + &
                                     donor_litt%bg_cwd(c,ilyr) * donor_weight
          this%bg_cwd_in(c,ilyr)   = this%bg_cwd_in(c,ilyr) * self_weight + &
                                     donor_litt%bg_cwd_in(c,ilyr) * donor_weight
          this%bg_cwd_frag(c,ilyr) = this%bg_cwd_frag(c,ilyr) * self_weight + &
                                     donor_litt%bg_cwd_frag(c,ilyr) * donor_weight
       end do

    end do

    
    do pft=1,npft
       
       this%seed(pft)            = this%seed(pft) * self_weight + &
                                   donor_litt%seed(pft) * donor_weight
       this%seed_germ(pft)       = this%seed_germ(pft) * self_weight + &
                                   donor_litt%seed_germ(pft) * donor_weight
       
       this%seed_in_local(pft)   = this%seed_in_local(pft) * self_weight + &
                                   donor_litt%seed_in_local(pft) * donor_weight
       this%seed_in_extern(pft)  = this%seed_in_extern(pft) * self_weight + &
                                   donor_litt%seed_in_extern(pft) * donor_weight
       
       this%seed_decay(pft)      = this%seed_decay(pft) * self_weight + &
                                   donor_litt%seed_decay(pft) * donor_weight
       this%seed_germ_decay(pft) = this%seed_germ_decay(pft) * self_weight + &
                                   donor_litt%seed_germ_decay(pft) * donor_weight
       this%seed_germ_in(pft)    = this%seed_germ_in(pft) * self_weight + &
                                   donor_litt%seed_germ_in(pft) * donor_weight
   end do


   do dcmpy=1,ndcmpy

       this%leaf_fines(dcmpy)      = this%leaf_fines(dcmpy) * self_weight + &
                                   donor_litt%leaf_fines(dcmpy) * donor_weight
       this%leaf_fines_in(dcmpy)   = this%leaf_fines_in(dcmpy) * self_weight + &
                                   donor_litt%leaf_fines_in(dcmpy) * donor_weight
       this%leaf_fines_frag(dcmpy) = this%leaf_fines_frag(dcmpy) * self_weight + &
                                   donor_litt%leaf_fines_frag(dcmpy) * donor_weight

       do ilyr=1,nlevsoil
           this%root_fines(dcmpy,ilyr)     = this%root_fines(dcmpy,ilyr) * self_weight + &
                                            donor_litt%root_fines(dcmpy,ilyr) * donor_weight
          this%root_fines_in(dcmpy,ilyr)   = this%root_fines_in(dcmpy,ilyr) * self_weight + &
                                            donor_litt%root_fines_in(dcmpy,ilyr) * donor_weight
          this%root_fines_frag(dcmpy,ilyr) = this%root_fines_frag(dcmpy,ilyr) * self_weight + &
                                            donor_litt%root_fines_frag(dcmpy,ilyr) * donor_weight
       end do
    end do

    return
  end subroutine FuseLitter

  ! =====================================================================================

  subroutine CopyLitter(this,donor_litt)

    ! This might not need to ever be called.  When a new patch is created
    ! from disturbance, litter initialization is handled elsewhere (EDPatchDynamics)


    class(litter_type) :: this
    type(litter_type),intent(in) :: donor_litt


    this%ag_cwd(:)      = donor_litt%ag_cwd(:)
    this%ag_cwd_in(:)   = donor_litt%ag_cwd_in(:)
    this%ag_cwd_frag(:) = donor_litt%ag_cwd_frag(:)
    
    this%bg_cwd(:,:)      = donor_litt%bg_cwd(:,:)
    this%bg_cwd_in(:,:)   = donor_litt%bg_cwd_in(:,:)
    this%bg_cwd_frag(:,:) = donor_litt%bg_cwd_frag(:,:)

    this%leaf_fines(:)    = donor_litt%leaf_fines(:)
    this%seed(:)          = donor_litt%seed(:)
    this%seed_germ(:)     = donor_litt%seed_germ(:)
    this%leaf_fines_in(:) = donor_litt%leaf_fines_in(:)
    this%seed_in_local(:) = donor_litt%seed_in_local(:)
    
    this%seed_in_extern(:)    = donor_litt%seed_in_extern(:)
    this%leaf_fines_frag(:)   = donor_litt%leaf_fines_frag(:)
    
    this%seed_decay(:)        = donor_litt%seed_decay(:)
    this%seed_germ_decay(:)   = donor_litt%seed_germ_decay(:)
    this%seed_germ_in(:)      = donor_litt%seed_germ_in(:)
    this%root_fines(:,:)      = donor_litt%root_fines(:,:)
    this%root_fines_in(:,:)   = donor_litt%root_fines_in(:,:)
    this%root_fines_frag(:,:) = donor_litt%root_fines_frag(:,:)

    return
  end subroutine CopyLitter

  ! =====================================================================================

  subroutine InitAllocate(this,numpft,numlevsoil,element_id)

    class(litter_type) :: this
    integer,intent(in)  :: numpft     ! number of plant functional types
    integer,intent(in)  :: numlevsoil ! number of soil layers
    integer,intent(in)  :: element_id ! PARTEH compliant element index

    this%element_id = element_id

    allocate(this%bg_cwd_in(ncwd,numlevsoil))
    allocate(this%bg_cwd(ncwd,numlevsoil))
    allocate(this%bg_cwd_frag(ncwd,numlevsoil))

    allocate(this%leaf_fines(ndcmpy))
    allocate(this%root_fines(ndcmpy,numlevsoil))
    allocate(this%leaf_fines_in(ndcmpy))
    allocate(this%root_fines_in(ndcmpy,numlevsoil))
    allocate(this%leaf_fines_frag(ndcmpy))
    allocate(this%root_fines_frag(ndcmpy,numlevsoil))

    allocate(this%seed_in_local(numpft))
    allocate(this%seed_in_extern(numpft))
    allocate(this%seed(numpft))
    allocate(this%seed_germ(numpft))
    allocate(this%seed_germ_in(numpft))
    allocate(this%seed_germ_decay(numpft))
    allocate(this%seed_decay(numpft))

    ! Initialize everything to a nonsense flag
    this%ag_cwd(:)            = fates_unset_r8
    this%bg_cwd(:,:)          = fates_unset_r8
    this%leaf_fines(:)        = fates_unset_r8
    this%root_fines(:,:)      = fates_unset_r8
    this%seed(:)              = fates_unset_r8
    this%seed_germ(:)         = fates_unset_r8

    this%ag_cwd_in(:)         = fates_unset_r8
    this%bg_cwd_in(:,:)       = fates_unset_r8
    this%leaf_fines_in(:)     = fates_unset_r8
    this%root_fines_in(:,:)   = fates_unset_r8
    this%seed_in_local(:)     = fates_unset_r8
    this%seed_in_extern(:)    = fates_unset_r8

    this%ag_cwd_frag(:)       = fates_unset_r8
    this%bg_cwd_frag(:,:)     = fates_unset_r8
    this%leaf_fines_frag(:)   = fates_unset_r8
    this%root_fines_frag(:,:) = fates_unset_r8

    this%seed_decay(:)        = fates_unset_r8
    this%seed_germ_decay(:)   = fates_unset_r8
    this%seed_germ_in(:)      = fates_unset_r8

    return
  end subroutine InitAllocate

  ! =====================================================================================

  subroutine InitConditions(this, &
                            init_leaf_fines, &
                            init_root_fines, &
                            init_ag_cwd,     &
                            init_bg_cwd,     &
                            init_seed,       &
                            init_seed_germ)
    
    ! This procedure initialized litter pools.  This does not allow initialization
    ! of each soil layer depth, or decomposability pool. This is meant for
    ! uniform initializations. This is used for cold-starts, but is not
    ! used in restarts.  For patch fusion, this routine is used to zero the pools
    ! before accumulating debris from multiple patches.


    class(litter_type) :: this
    real(r8),intent(in) :: init_leaf_fines
    real(r8),intent(in) :: init_root_fines
    real(r8),intent(in) :: init_ag_cwd
    real(r8),intent(in) :: init_bg_cwd
    real(r8),intent(in) :: init_seed
    real(r8),intent(in) :: init_seed_germ
    
    this%ag_cwd(:)              = init_ag_cwd
    this%bg_cwd(:,:)            = init_bg_cwd
    this%leaf_fines(:)          = init_leaf_fines
    this%root_fines(:,:)        = init_root_fines
    this%seed(:)                = init_seed
    this%seed_germ(:)           = init_seed_germ

    return
  end subroutine InitConditions
  
  ! =====================================================================================
  
  subroutine DeallocateLitt(this)
    
    class(litter_type) :: this

    deallocate(this%bg_cwd)
    deallocate(this%leaf_fines)
    deallocate(this%root_fines)
    deallocate(this%seed)
    deallocate(this%seed_germ)

    deallocate(this%bg_cwd_in)
    deallocate(this%leaf_fines_in)
    deallocate(this%root_fines_in)
    deallocate(this%seed_in_local)
    deallocate(this%seed_in_extern)
    
    deallocate(this%bg_cwd_frag)
    deallocate(this%leaf_fines_frag)
    deallocate(this%root_fines_frag)
   
    deallocate(this%seed_decay)
    deallocate(this%seed_germ_decay)
    deallocate(this%seed_germ_in)

    return
  end subroutine DeallocateLitt
  
  ! =====================================================================================

  subroutine ZeroFlux(this)
    
    class(litter_type) :: this
    
    this%ag_cwd_in(:)         = 0._r8
    this%bg_cwd_in(:,:)       = 0._r8
    this%leaf_fines_in(:)     = 0._r8
    this%root_fines_in(:,:)   = 0._r8
    this%seed_in_local(:)     = 0._r8
    this%seed_in_extern(:)    = 0._r8

    this%ag_cwd_frag(:)       = 0._r8
    this%bg_cwd_frag(:,:)     = 0._r8
    this%leaf_fines_frag(:)   = 0._r8
    this%root_fines_frag(:,:) = 0._r8
    
    this%seed_germ_in(:)      = 0._r8
    this%seed_decay(:)        = 0._r8
    this%seed_germ_decay(:)   = 0._r8


    return
  end subroutine ZeroFlux

  ! ===================================================

  function GetTotalLitterMass(this) result(total_mass)
    
    class(litter_type) :: this
    real(r8) :: total_mass
    
    total_mass = sum(this%ag_cwd) + &
                 sum(this%bg_cwd) + &
                 sum(this%root_fines) + &
                 sum(this%leaf_fines) + & 
                 sum(this%seed) + & 
                 sum(this%seed_germ)
    
    return
  end function GetTotalLitterMass
  
  ! =====================================================

  subroutine adjust_SF_CWD_frac(dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)
     
     !DESCRIPTION
     !Adjust  the partitioning of struct + sawp into cwd pools based on 
     !cohort dbh. This avoids struct and sapw from small cohorts going to
     !fuel classes that are larger than their dbh. Here, struct + sapw are sent to the cwd
     !class consistent with fuel class diameter thresholds (Fosberg et al., 1971;
     !Rothermel, 1983)

     !ARGUMENTS
     real(r8), intent(in)               :: dbh !dbh of cohort [cm]
     integer, intent(in)                :: ncwd !number of cwd pools
     real(r8), intent(in)               :: SF_val_CWD_frac(:) !fates parameter specifying the
                                                              !fraction of struct + sapw going
                                                              !to each CWD class
     real(r8), intent(out)              :: SF_val_CWD_frac_adj(:) !Updated cwd paritions
     !
     !LOCAL VARIABLES
     !These diameter ranges are based on work by Fosberg et al., 1971 
     
     real(r8), parameter :: lb_max_diam        = 7.6 !max diameter [cm] for large branch (100 hr fuel)
     real(r8), parameter :: sb_max_diam        = 2.5 !max diameter [cm] for small branch (10 hr fuel)
     real(r8), parameter :: twig_max_diam      = 0.6 !max diameter [cm] for twig (1 hr fuel)
     !------------------------------------------------------------------------------------

     
     SF_val_CWD_frac_adj = SF_val_CWD_frac
     
     !If dbh is > 7.6 cm diameter then the main stem is 1,000 hr fuel and we don't change
     !how biomass is partitioned among cwd classes.
     if (dbh > lb_max_diam) then
        return

     !When dbh is greater than the max size of a small branch (10 hr fuel) but less than or 
     !equal to the max size of a large branch (e.g. saplings where dbh > 2.5 cm and < 7.5 cm)
     !we redistribute the biomass that would have gone to 1,000 hr fuel among the smaller cwd classes.
     !This keeps the biomass proportions among the combustible classes the same as the scenario 
     !where a fraction of struct + sawp goes to 1,000 hr fuel.
     else if (dbh > sb_max_diam .and. dbh .le. lb_max_diam) then
        SF_val_CWD_frac_adj(ncwd) = 0.0
        SF_val_CWD_frac_adj(ncwd-1) = SF_val_CWD_frac(ncwd-1) / (1.0_r8 - SF_val_CWD_frac(ncwd) )
	SF_val_CWD_frac_adj(ncwd-2) = SF_val_CWD_frac(ncwd-2) / (1.0_r8 - SF_val_CWD_frac(ncwd) )
	SF_val_CWD_frac_adj(ncwd-3) = SF_val_CWD_frac(ncwd-3) / (1.0_r8 - SF_val_CWD_frac(ncwd) )

     !When dbh is greater than the max size of a twig (1 hr fuel) but less than or 
     !equal to the max size of a small branch (10 hr fuel) we redistribute the biomass among the smaller
     !pools.
     else if (dbh > twig_max_diam .and. dbh .le. sb_max_diam) then
        SF_val_CWD_frac_adj(ncwd) = 0.0
        SF_val_CWD_frac_adj(ncwd-1) = 0.0
        SF_val_CWD_frac_adj(ncwd-2) = SF_val_CWD_frac(ncwd-2) / (1.0_r8 - (SF_val_CWD_frac(ncwd) + &
                                      SF_val_CWD_frac(ncwd-1)))
	SF_val_CWD_frac_adj(ncwd-3) = SF_val_CWD_frac(ncwd-3) / (1.0_r8 - (SF_val_CWD_frac(ncwd) + &
                                      SF_val_CWD_frac(ncwd-1)))
                                              
     !If dbh is less than or equal to the max size of a twig we send all 
     !biomass to twigs.
     else if (dbh .le. twig_max_diam) then
        SF_val_CWD_frac_adj(ncwd) = 0.0
        SF_val_CWD_frac_adj(ncwd-1) = 0.0
	SF_val_CWD_frac_adj(ncwd-2) = 0.0
	SF_val_CWD_frac_adj(ncwd-3) = sum(SF_val_CWD_frac)

     endif 
  end subroutine adjust_SF_CWD_frac
  
end module FatesLitterMod
