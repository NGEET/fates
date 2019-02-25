module FatesLittMod

  ! -------------------------------------------------------------------------------------
  ! This module contains methods and type definitions for all things litter.
  ! "litter" is used here in the broadest sense.  It means all organic material
  ! that is no longer associated with a live plant.
  !
  ! This encompasses:  1) "Coarse Woody Debris"
  !                    2) fine materials leaves, roots etc 
  !                       (sometimes exclusively refered to as litter)
  !                    3) Reproductive materials (seeds, nuts, fruits)
  !
  ! Important point:   THESE POOLS DO NOT CONTAIN DECOMPOSING MATTER !!!!!
  ! 
  ! Continued: These litter pools will fragment, and then be passed to a 
  ! soil-biogeochemical model (NOT FATES) to handle decomposition.
  ! -------------------------------------------------------------------------------------


   ! To-do:
   ! 0) Update the initialization sequence to use
   !    the array style, and also set the element_id
   ! 1) Add the state integration step
   ! 2) Evaluate how we are handling disturbance
   ! 3) Add root fraction variable to each cohort?
   ! 4) Tracking site level mass balances, where should the new litter structure 
   !    help in tracking that sort of thing? 
   ! 5) Add a routine that performs all of the site-level mass balance
   !    integrations from daily sources?

   
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : i4 => fates_int
   use FatesConstantsMod, only : nearzero
   use FatesConstantsMod, only : calloc_abs_error
   use FatesConstantsMod, only : fates_unset_real

   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log 
   use shr_log_mod      , only : errMsg => shr_log_errMsg

   implicit none
   private
   
   type litt_vartype

      
      ! This object is allocated for each element (C, N, P, etc) that we wish to track.

      integer              :: element_id               ! This element ID should
                                                       ! be associated with the element
                                                       ! types listed in parteh/PRTGenericMod.F90

      ! Prognostic variables (litter and coarse woody debris)

      real(r8)             ::  ag_cwd(ncwd)            ! above ground coarse wood debris kg/m2
      real(r8),allocatable ::  bg_cwd(:,:)             ! below ground coarse wood debris kg/m2
      real(r8),allocatable ::  leaf_fines(:)           ! above ground leaf litter        kg/m2
      real(r8),allocatable ::  root_fines(:,:)         ! below ground fine root litter   kg/m2

      real(r8),allocatable ::  seed(:)                 ! the seed pool (viable)          kg/m2


      ! Fluxes in - dying trees / seed rain
      
      real(r8)             ::  ag_cwd_in(ncwd)         ! kg/m2/day
      real(r8),allocatable ::  bg_cwd_in(:,:)          ! kg/m2/day
      real(r8),allocatable ::  leaf_fines_in(:)        ! kg/m2/day
      real(r8),allocatable ::  root_fines_in(:,:)      ! kg/m2/day

      real(r8),allocatable ::  seed_in_local(:)        ! kg/m2/day (from local sources)
      real(r8),allocatable ::  seed_in_extern(:)       ! kg/m2/day (from outside cell)

      ! Fluxes out - fragmentation
      
      real(r8)             ::  ag_cwd_frag(ncwd)     ! kg/m2/day
      real(r8),allocatable ::  bg_cwd_frag(:,:)      ! kg/m2/day
      real(r8),allocatable ::  leaf_fines_frag(:)    ! kg/m2/day
      real(r8),allocatable ::  root_fines_frag(:,:)  ! kg/m2/day

      ! Fluxes out from burning are not tracked here
      ! because this process changes the size of the patch
      ! as well, and the unit fluxes that would be tracked
      ! here are convoluted
      
      !real(r8)             ::  ag_cwd_burn_atm(ncwd)     ! kg/m2/day
      !real(r8),allocatable ::  bg_cwd_burn_atm(:,:)      ! kg/m2/day
      !real(r8),allocatable ::  leaf_fines_burn_atm(:)    ! kg/m2/day
      !real(r8),allocatable ::  root_fines_burn_atm(:,:)  ! kg/m2/day


      ! Fluxes out - germination
      real(r8),allocatable :: seed_germ(:)           ! kg/m2/day


      ! Flux out from exporting (harvesting, seed eflux?)
      real(r8) :: exported                             ! kg/m2/day


      ! Flux (in/out ... transfer) from seed to litter
      real(r8),allocatable :: seed_decay(:)            ! kg/m2/day


      


    contains
      
      procedure,non_overridable :: InitAllocate
      procedure,non_overridable :: DeallocateLitt
      procedure,non_overridable :: InitConditions
      procedure,non_overridable :: ZeroFlux
      
   end type litt_vartype

   ! Part 3: Public extended types

   character(len=*), parameter, private :: sourcefile = __FILE__

   
contains

  subroutine InitAllocate(this,numpft,numlevsoil)

    class(litt_vartype) :: this
    integer,intent(in)  :: numpft   ! number of plant functional types
    integer,intent(in)  :: numlevsoil ! number of soil layers

    allocate(this%bg_cwd(ncwd,numlevsoil))
    allocate(this%leaf_fines(numpft))
    allocate(this%root_fines(numpft,numlevsoil))
    allocate(this%seed(numpft))

    allocate(this%bg_cwd_in(ncwd,numlevsoil))
    allocate(this%leaf_fines_in(numpft))
    allocate(this%root_fines_in(numpft,numlevsoil))
    allocate(this%seed_in_local(numpft))
    allocate(this%seed_in_extern(numpft))

    allocate(this%bg_cwd_frag(ncwd,numlevsoil))
    allocate(this%leaf_fines_frag(numpft))
    allocate(this%root_fines_frag(numpft,numlevsoil))

    allocate(this%bg_cwd_burn(ncwd,numlevsoil))
    allocate(this%leaf_fines_burn(numpft))
    allocate(this%root_fines_burn(numpft,numlevsoil))

    allocate(this%seed_germ(numpft))
    allocate(this%seed_decay(numpft))


    ! Initialize everything to a nonsense flag
    this%ag_cwd(:)               = fates_unset_real
    this%bg_cwd(:,:)             = fates_unset_real
    this%leaf_fines(:)          = fates_unset_real
    this%root_fines(:,:)        = fates_unset_real
    this%seed(:)                 = fates_unset_real

    this%ag_cwd_in(:)            = fates_unset_real
    this%bg_cwd_in(:,:)          = fates_unset_real
    this%leaf_fines_in(:)       = fates_unset_real
    this%root_fines_in(:,:)     = fates_unset_real
    this%seed_in_local(:)            = fates_unset_real
    this%seed_in_extern(:)      = fates_unset_real

    this%ag_cwd_frag(:)        = fates_unset_real
    this%bg_cwd_frag(:,:)      = fates_unset_real
    this%leaf_fines_frag(:)   = fates_unset_real
    this%root_fines_frag(:,:) = fates_unset_real

    this%ag_cwd_burn(:)        = fates_unset_real
    this%bg_cwd_burn(:,:)      = fates_unset_real
    this%leaf_fines_burn(:)   = fates_unset_real
    this%root_fines_burn(:,:) = fates_unset_real

    this%exported               = fates_unset_real

    this%seed_germ(:)         = fates_unset_real
    this%seed_decay(:)        = fates_unset_real

    return
  end subroutine InitAllocate

  ! =====================================================================================

  subroutine InitConditions(this, &
                            init_leaf_fines, &
                            init_root_fines, &
                            init_ag_cwd,     &
                            init_bg_cwd,     &
                            init_seed)
    
    class(litt_vartype) :: this
    real(r8),intent(in) :: init_leaf_fines
    real(r8),intent(in) :: init_root_fines
    real(r8),intent(in) :: init_ag_cwd
    real(r8),intent(in) :: init_bg_cwd
    real(r8),intent(in) :: init_seed
    
    this%ag_cwd(:)              = init_leaf_fines
    this%bg_cwd(:,:)            = init_root_fines
    this%leaf_fines(:)          = init_ag_cwd
    this%root_fines(:,:)        = init_bg_cwd
    this%seed(:)                = init_seed
    
    return
  end subroutine InitConditions
  
  ! =====================================================================================
  
  subroutine DeallocateLitt(this)
    
    class(litt_vartype) :: this

    deallocate(this%bg_cwd)
    deallocate(this%leaf_fines)
    deallocate(this%root_fines)
    deallocate(this%seed)

    deallocate(this%bg_cwd_in)
    deallocate(this%leaf_fines_in)
    deallocate(this%root_fines_in)
    deallocate(this%seed_in_local)
    deallocate(this%seed_in_extern)
    
    deallocate(this%bg_cwd_frag)
    deallocate(this%leaf_fines_frag)
    deallocate(this%root_fines_frag)
   
    deallocate(this%bg_cwd_burn)
    deallocate(this%leaf_fines_burn)
    deallocate(this%root_fines_burn)
    
    deallocate(this%seed_decay)
    deallocate(this%seed_germ)

    return
  end subroutine DeallocateLitt
  
  ! =====================================================================================

  subroutine ZeroFlux(this)
    
    class(litt_vartype) :: this
    
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
    
    this%ag_cwd_burn(:)       = 0._r8
    this%bg_cwd_burn(:,:)     = 0._r8
    this%leaf_fines_burn(:)   = 0._r8
    this%root_fines_burn(:,:) = 0._r8

    this%seed_germ(:)         = 0._r8
    this%seed_decay(:)        = 0._r8

    this%exported             = 0._r8
    
    return
  end subroutine ZeroFlux



  
end module FatesLittMod
