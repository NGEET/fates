module CLMFatesInterfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library/API with the CLM/ALM/ATS/etc model driver.  
   ! All connections between the two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see FATES memory 
   ! structures.
   ! The routines here, that call FATES library routines, will not pass any types defined
   ! by the driving land model (HLM).
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
   use clm_time_manager  , only : get_days_per_year, get_curr_date
   use clm_time_manager  , only : get_ref_date, timemgr_datediff 
   use spmdMod           , only : masterproc

   ! Used FATES Modules
   use FatesInterfaceMod     , only : fates_interface_type
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDTypesMod            , only : udata
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site

   implicit none

   type, public :: hlm_fates_interface_type
      
      !      private
      

      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_interface_type), allocatable :: fates (:)
      
      ! fates2hlm_inst (previously called "clm_ed_inst") contains types and variables
      ! that are passed back to the driving land model, ie fates-to-hlm.  
      ! usefull to a calling model.  In this case HLM means "Hosting Land Model"
      ! prev:  type(ed_clm_type)::ed_clm_inst
      type(ed_clm_type) :: fates2hlm_inst   

      
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
      procedure, public :: fates2hlm_link
      procedure, public :: dynamics_driv
      
   end type hlm_fates_interface_type

contains
   
   subroutine init(this,bounds_proc)
      
      ! ---------------------------------------------------------------------------------
      ! This initializes the dlm_fates_interface_type 
      !
      ! sites is the root of the ED state hierarchy (instantaneous info on 
      ! the state of the ecosystem).  As such, it governs the connection points between
      ! the host (which also dictates its allocation) and its patch structures.
      !
      ! sites may associate with different scales in different models. In
      ! CLM, it is being designed to relate to column scale.
      !
      ! This global may become relegated to this module. 
      !
      ! Note: CLM/ALM currently wants sites to be allocated even if ed
      ! is not turned on
      ! ---------------------------------------------------------------------------------
      
      implicit none
      
      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_proc
      
      
      ! Initialize the FATES communicators with the HLM
      ! This involves to stages
      ! 1) allocate the vectors hlm_inst
      ! 2) add the history variables defined in clm_inst to the history machinery
      call this%fates2hlm_inst%Init(bounds_proc)
      
      ! Initialize ED phenology variables
      ! This also involves two stages
      ! 1) allocate the vectors in phen_inst
      ! 2) add the phenology history variables to the history machinery
      call this%phen_inst%Init(bounds_proc)
      
      ! ---------------------------------------------------------------------------------
      ! Initialization of the state-threads is handled by the calling subroutine
      ! clm_instInit
      ! ---------------------------------------------------------------------------------
      
   end subroutine init
   
   ! ------------------------------------------------------------------------------------
   
   subroutine fates2hlm_link(this, bounds_clump, nc, waterstate_inst, canopystate_inst)
      
      ! CLM:  called from initialize2()
      ! ALM:  ??
      
      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_clump
      type(waterstate_type)   , intent(inout)        :: waterstate_inst
      type(canopystate_type)  , intent(inout)        :: canopystate_inst
      integer, intent(in)                            :: nc
      
      call this%fates2hlm_inst%ed_clm_link( bounds_clump,                  &
            this%fates(nc)%sites(bounds_clump%begg:bounds_clump%endg), &
            this%phen_inst,                                                &
            waterstate_inst,                                               &
            canopystate_inst)
      
      
      return
   end subroutine fates2hlm_link
  
   ! ------------------------------------------------------------------------------------
   
   subroutine dynamics_driv(this, nc, bounds_clump,      &
         atm2lnd_inst, soilstate_inst, temperature_inst, &
         waterstate_inst, canopystate_inst)
    
      ! This wrapper is called daily from clm_driver
      ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
      ! ed_driver is not a hlm_fates_inst_type procedure because we need an extra step 
      ! to process array bounding information 
      
      implicit none
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_clump
      type(atm2lnd_type)      , intent(in)           :: atm2lnd_inst
      type(soilstate_type)    , intent(in)           :: soilstate_inst
      type(temperature_type)  , intent(in)           :: temperature_inst
      integer                 , intent(in)           :: nc
      type(waterstate_type)   , intent(inout)        :: waterstate_inst
      type(canopystate_type)  , intent(inout)        :: canopystate_inst

      ! !LOCAL VARIABLES:
      real(r8) :: dayDiff                  ! day of run
      integer  :: dayDiffInt               ! integer of day of run
      integer  :: g                        ! gridcell  
      integer  :: yr                       ! year (0, ...)
      integer  :: mon                      ! month (1, ..., 12)
      integer  :: day                      ! day of month (1, ..., 31)
      integer  :: sec                      ! seconds of the day
      integer  :: ncdate                   ! current date
      integer  :: nbdate                   ! base date (reference date)
      !-----------------------------------------------------------------------

      
      ! ---------------------------------------------------------------------------------
      ! INTERF-TODO: REMOVE ED_DRIVER ARGUMENTS OF CLM STUCTURED TYPES AND
      ! REPLACE THEM WITH FATES_BC TYPES WITH ITS OWN MAPPING SCHEME
      ! ALSO, NOTE THAT THE ED_DYNAMICS IS A MODULE OF FATES NOW
      ! ie:
      ! fates(nc)%fatesbc%leaf_temp <=> canopystate_inst%
      !
      ! call this%fates(nc)%ed_driver(this%fates(nc)%site,    &
      !                               this%fates(nc)%fatesbc)
      ! ---------------------------------------------------------------------------------
      

      call this%fates2hlm_inst%SetValues( bounds_clump, 0._r8 )

      ! timing statements. 
      udata%n_sub = get_days_per_year()
            udata%deltat = 1.0_r8/dble(udata%n_sub) !for working out age of patches in years        
      if(udata%time_period == 0)then             
         udata%time_period = udata%n_sub
      endif
      
      call get_curr_date(yr, mon, day, sec)
      ncdate = yr*10000 + mon*100 + day
      call get_ref_date(yr, mon, day, sec)
      nbdate = yr*10000 + mon*100 + day
      
      call timemgr_datediff(nbdate, 0, ncdate, sec, dayDiff)
      
      dayDiffInt = floor(dayDiff)
      udata%time_period = mod( dayDiffInt , udata%n_sub )
      

      ! TODO-INTEF: PROCEDURE FOR CONVERTING CLM/ALM FIELDS TO MODEL BOUNDARY
      ! CONDITIONS. IE. 


      ! where most things happen
      do g = bounds_clump%begg,bounds_clump%endg
         if (this%fates(nc)%sites(g)%istheresoil) then
            call ed_ecosystem_dynamics(this%fates(nc)%sites(g), &
                  this%fates2hlm_inst,  &
                  this%phen_inst, atm2lnd_inst, &
                  soilstate_inst, temperature_inst, waterstate_inst)
            
            call ed_update_site(this%fates(nc)%sites(g))
         endif
      enddo

      ! link to CLM/ALM structures
      call this%fates2hlm_inst%ed_clm_link( bounds_clump,                  &
            this%fates(nc)%sites(bounds_clump%begg:bounds_clump%endg),     &
            this%phen_inst,                                                &
            waterstate_inst,                                               &
            canopystate_inst)


      if (masterproc) then
         write(iulog, *) 'clm: leaving ED model', bounds_clump%begg, &
                                                  bounds_clump%endg, dayDiffInt
      end if

      
      return
   end subroutine dynamics_driv
   
   ! ------------------------------------------------------------------------------------
   !  THESE WRAPPERS MAY COME IN HANDY, KEEPING FOR NOW
   !  subroutine set_fates2hlm_inst(this,bounds_clump, setval_scalar) 
   !
   !    ! This is a simple wrapper to flush some FATES -> CLM communicators
   !
   !    implicit none
   !    class(hlm_fates_interface_type), intent(inout) :: this
   !    type(bounds_type),intent(in)                :: bounds_clump
   !    real(r8),intent(in) :: setval_scalar      ! This value will flush all 
   !    
   !    call this%fates2hlm_inst%SetValues( bounds_clump, setval_scalar )
   !  end subroutine set_fates2hlm_inst
   ! ------------------------------------------------------------------------------------
   !  subroutine phen_accvars_init(this,bounds_clump)
   !
   !    implicit none
   !    class(hlm_fates_interface_type), intent(inout) :: this
   !    type(bounds_type),intent(in)                :: bounds_clump
   !
   !    call this%phen_inst%initAccVars(bounds_clump)
   !
   !    return
   !  end subroutine phen_accvars_init
   ! ------------------------------------------------------------------------------------
    
end module CLMFatesInterfaceMod
