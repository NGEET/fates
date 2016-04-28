module clmfates_interfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library with the CLM/ALM/ATS/etc model driver.  
   ! All connections between the two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see FATES memory 
   ! structures.
   ! The routines here, that call FATES library routines, will not pass any types defined
   ! by the driving land model (DLM).
   ! 
   ! either native type arrays (int,real,log, etc) or packed into ED boundary condition
   ! structures.
   !
   ! Note that CLM/ALM does use Shared Memory Parallelism (SMP), where processes such as 
   ! the update of state variables are forked.  However, IO is not assumed to be 
   ! threadsafe and therefore memory spaces reserved for IO must be continuous vectors,
   ! and moreover they must be pushed/pulled from history IO for each individual 
   ! bounds_proc memory space as a unit.
   !
   ! Therefore, the state variables in the clm_fates communicator is vectorized by
   ! threadcount, and the IO communication arrays are not.
   !
   !
   ! Conventions:
   ! keep line widths within 90 spaces
   ! DLM acronym = Driving Land Model
   !
   ! -------------------------------------------------------------------------------------

   !  use ed_driver_interface, only: 
   
   ! Used CLM Modules
   use PatchType         , only : patch
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use decompMod         , only : bounds_type
   use WaterStateType    , only : waterstate_type
   use CanopyStateType   , only : canopystate_type
   use TemperatureType   , only : temperature_type
   use SoilStateType     , only : soilstate_type
   use clm_varctl        , only : iulog
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use clm_time_manager  , only : is_restart
   use ncdio_pio         , only : file_desc_t


   ! Used ED Modules
   use EDtypesMod            , only : ed_patch_type, ed_site_type, numpft_ed
   use EDtypesMod            , only : map_clmpatch_to_edpatch
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs
   use EDMainMod             , only : ed_driver
   use EDRestVectorMod       , only : EDRest
   use EDInitMod             , only : ed_init_sites
   use EDMainMod             , only : ed_update_site

   implicit none


   type, public :: fates_thread_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated
      ! with a grid-cell, this is intended to be migrated to columns
      ! prev:  type(ed_site_type)::ed_allsites_inst
      type(ed_site_type), allocatable :: site_inst(:)
      
      ! INTERF-TODO ADD THE DLM->FATES BOUNDARY CONDITION CLASS
      ! These are boundary condition variables populated by the DLM
      ! type(fates_bc_type) :: dlm2fates
      
   contains

      ! Procedures for initializing FATES threaded memory and communicators
      procedure, public :: thread_init
      procedure, public :: thread_clean
      procedure, public :: site_init
      procedure, public :: thread_restart
      procedure, public :: canopy_sunshade_fracs

   end type fates_thread_type
   
   ! ------------------------------------------------------------------------------------
   
   type, public :: dlm_fates_interface_type
      
      !      private
      

      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_thread_type), allocatable :: thread (:)
      
      ! fates2dlm_inst (previously called "clm_ed_inst") contains types and variables
      ! that are passed back to the driving land model, ie fates-to-dlm.  
      ! usefull to a calling model.  In this case DLM means "Driving Land Model"
      ! prev:  type(ed_clm_type)::ed_clm_inst
      type(ed_clm_type) :: fates2dlm_inst   

      
      ! These are phenology relevant variables (strange how phenology gets
      ! its own subclass)
      ! prev: type(ed_phenology_type)::ed_phenology_inst
      ! SOON TO BE DEPRECATED - PHENOLOGY TO BE MOVED INTO ed_sites(:)
      type(ed_phenology_type)   :: phen_inst  


      ! INTERF-TODO: we will need a new bounding type (maybe?)
      ! depending on how we structure the memory
      ! We will likely have a fates_bcs (boundary conditions) type
      ! And this type could take many forms, but probably be pointers
      ! to various CLM/ALM/ATS/etc types
      ! This is a structure of pointers that maps between
      ! the calling model and the fates model
      !      type(fates_bounds_type) :: bound_fate


   contains

     procedure, public :: init
     procedure, public :: fates2dlm_link
     procedure, public :: dynamics_driv
     
  end type dlm_fates_interface_type

contains

   subroutine thread_init(this,bounds_clump)
      

      implicit none
      
      ! Input Arguments
      class(fates_thread_type), intent(inout) :: this
      type(bounds_type),intent(in)            :: bounds_clump 
      
      
      ! Initialize the mapping elements between FATES and the DLM
      
      ! These bounds are for a single clump (thread)
      allocate (this%site_inst(bounds_clump%begg:bounds_clump%endg))
      
      return
   end subroutine thread_init

   ! ------------------------------------------------------------------------------------

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine thread_clean(this,bounds_clump)
      
      implicit none
      
      ! Input Arguments
      class(fates_thread_type), intent(inout) :: this
      type(bounds_type),intent(in)            :: bounds_clump 
      
      ! Incrementally walk through linked list and deallocate
      
      ! Deallocate the site list
      deallocate (this%site_inst)
      
      return
   end subroutine thread_clean

   ! ------------------------------------------------------------------------------------
   
   subroutine site_init(this,bounds_clump)
         
      ! CLM:  called from initialize2()
      ! ALM:  ??
      
      ! Input Arguments
      class(fates_thread_type), intent(inout) :: this
      type(bounds_type),intent(in)            :: bounds_clump
      
      ! locals
      integer :: g
      
      ! Initialize  (INTERF-TODO THIS ROUTINE CALLS CLM STUFF-MIGRATE CODE TO HERE)
      call ed_init_sites( bounds_clump,                                               &
            this%site_inst(bounds_clump%begg:bounds_clump%endg) )
      
      do g = bounds_clump%begg,bounds_clump%endg
         if (this%site_inst(g)%istheresoil) then
            call ed_update_site(this%site_inst(g))
         end if
      end do
      
      return
   end subroutine site_init
   
   ! ------------------------------------------------------------------------------------

   subroutine thread_restart(this, bounds_clump, ncid, flag )
      
      implicit none
      class(fates_thread_type), intent(inout)     :: this
      type(bounds_type),intent(in)                :: bounds_clump
      type(file_desc_t)       , intent(inout)     :: ncid    ! netcdf id
      character(len=*)        , intent(in)        :: flag    !'read' or 'write'
      
      call EDRest( bounds_clump, this%site_inst(bounds_clump%begg:bounds_clump%endg), &
            ncid, flag )
      return
   end subroutine thread_restart
   
   ! ------------------------------------------------------------------------------------
   
   subroutine canopy_sunshade_fracs(this,filter_nourbanp, num_nourbanp, &
         atm2lnd_inst,canopystate_inst)
         
         
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      
      implicit none
      
      ! Input Arguments
      class(fates_thread_type), intent(inout) :: this
      
      ! patch filter for non-urban points
      integer, intent(in),dimension(:)     :: filter_nourbanp
      
      ! number of patches in non-urban points in patch  filter
      integer, intent(in)                  :: num_nourbanp       
      
      ! direct and diffuse downwelling radiation (W/m2)
      type(atm2lnd_type),intent(in)        :: atm2lnd_inst
      
      ! Input/Output Arguments to CLM
      type(canopystate_type),intent(inout) :: canopystate_inst
      
      ! Local Variables
      integer  :: fp                          ! non-urban filter patch index
      integer  :: p                           ! patch index
      integer  :: g                           ! grid cell index
      integer, parameter :: ipar = 1          ! The band index for PAR
      type(ed_patch_type), pointer :: cpatch  ! c"urrent" patch
      
      associate( forc_solad => atm2lnd_inst%forc_solad_grc, &
            forc_solai => atm2lnd_inst%forc_solai_grc, &
            fsun       => canopystate_inst%fsun_patch)
        
        do fp = 1,num_nourbanp
           
           p = filter_nourbanp(fp)
           g = patch%gridcell(p)
           
           if ( patch%is_veg(p) ) then 
              cpatch => map_clmpatch_to_edpatch(this%site_inst(g), p) 
              
              call ED_SunShadeFracs(cpatch,forc_solad(g,ipar),forc_solai(g,ipar),fsun(p))
              
           endif
           
        end do
      end associate
      return
   end subroutine canopy_sunshade_fracs
    
   subroutine init(this,bounds_proc)
    
    ! ---------------------------------------------------------------------------------
    ! This initializes the dlm_fates_interface_type 
    !
    ! site_inst is the root of the ED state hierarchy (instantaneous info on 
    ! the state of the ecosystem).  As such, it governs the connection points between
    ! the host (which also dictates its allocation) and its patch structures.
    !
    ! site_inst may associate with different scales in different models. In
    ! CLM, it is being designed to relate to column scale.
    !
    ! This global may become relegated to this module. 
    !
    ! Note: CLM/ALM currently wants site_inst to be allocated even if ed
    ! is not turned on
    ! ---------------------------------------------------------------------------------
    
    implicit none
    
    ! Input Arguments
    class(dlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                   :: bounds_proc
    
    
    ! Initialize the FATES communicators with the DLM
    ! This involves to stages
    ! 1) allocate the vectors dlm_inst
    ! 2) add the history variables defined in clm_inst to the history machinery
    call this%fates2dlm_inst%Init(bounds_proc)
    
    ! Initialize ED phenology variables
    ! This also involves two stages
    ! 1) allocate the vectors in phen_inst
    ! 2) add the phenology history variables to the history machinery
    call this%phen_inst%Init(bounds_proc)
    
    ! -----------------------------------------------------------------------------------
    ! Initialization of the state-threads is handled by the calling subroutine
    ! clm_instInit
    ! -----------------------------------------------------------------------------------
    
  end subroutine init
  
  ! -----------------------------------------------------------------------------------
  

  
  ! -----------------------------------------------------------------------------------
  
  subroutine fates2dlm_link(this,bounds_clump,  nc, waterstate_inst, canopystate_inst)
    
    ! CLM:  called from initialize2()
    ! ALM:  ??
    
    ! Input Arguments
    class(dlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    type(waterstate_type)   , intent(inout)     :: waterstate_inst
    type(canopystate_type)  , intent(inout)     :: canopystate_inst
    integer, intent(in)                         :: nc
    
    call this%fates2dlm_inst%ed_clm_link( bounds_clump,                      &
          this%thread(nc)%site_inst(bounds_clump%begg:bounds_clump%endg),    &
          this%phen_inst,                                                    &
          waterstate_inst,                                                   &
          canopystate_inst)
    
    
    return
 end subroutine fates2dlm_link
  
 ! -------------------------------------------------------------------------------------
  
 subroutine dynamics_driv(this, nc, bounds_clump,                                      &
       atm2lnd_inst, soilstate_inst, temperature_inst,                                 &
       waterstate_inst, canopystate_inst)
    
    ! This wrapper is called daily from clm_driver
    ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
    ! ed_driver is not a dlm_fates_inst_type procedure because we need an extra step 
    ! to process array bounding information 
    
    implicit none
    class(dlm_fates_interface_type), intent(inout)     :: this
    type(bounds_type),intent(in)                :: bounds_clump
    type(atm2lnd_type)      , intent(in)        :: atm2lnd_inst
    type(soilstate_type)    , intent(in)        :: soilstate_inst
    type(temperature_type)  , intent(in)        :: temperature_inst
    integer                 , intent(in)        :: nc
    type(waterstate_type)   , intent(inout)     :: waterstate_inst
    type(canopystate_type)  , intent(inout)     :: canopystate_inst
    
    ! INTERF-TODO: REMOVE ED_DRIVER ARGUMENTS OF CLM STUCTURED TYPES AND
    ! REPLACE THEM WITH FATES_BC TYPES WITH A BOUNDS MAPPING SCHEME
    
    ! RENAME TO fates_driver()
    call ed_driver( bounds_clump,                                                        &
          this%thread(nc)%site_inst(bounds_clump%begg:bounds_clump%endg),                &
          this%fates2dlm_inst,                                                           &
          this%phen_inst,                                                                &
          atm2lnd_inst, soilstate_inst, temperature_inst,                                &
          waterstate_inst, canopystate_inst)
    
    return
 end subroutine dynamics_driv
  
 ! -------------------------------------------------------------------------------------
 
!  THESE WRAPPERS MAY COME IN HANDY, KEEPING FOR NOW

!  subroutine set_fates2dlm_inst(this,bounds_clump, setval_scalar) 
!
!    ! This is a simple wrapper to flush some FATES -> CLM communicators
!
!    implicit none
!    class(dlm_fates_interface_type), intent(inout) :: this
!    type(bounds_type),intent(in)                :: bounds_clump
!    real(r8),intent(in) :: setval_scalar      ! This value will flush all 
!    
!    call this%fates2dlm_inst%SetValues( bounds_clump, setval_scalar )
!  end subroutine set_fates2dlm_inst
! -------------------------------------------------------------------------------------
!  subroutine phen_accvars_init(this,bounds_clump)
!
!    implicit none
!    class(dlm_fates_interface_type), intent(inout) :: this
!    type(bounds_type),intent(in)                :: bounds_clump
!
!    call this%phen_inst%initAccVars(bounds_clump)
!
!    return
!  end subroutine phen_accvars_init
  
  ! -------------------------------------------------------------------------------------







  ! ---------------------------------------------------------------------------------------
  

    
    
  end module clmfates_interfaceMod
