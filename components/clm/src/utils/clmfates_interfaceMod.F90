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
   use clm_varpar        , only : numpft
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use clm_time_manager  , only : is_restart
   use ncdio_pio         , only : file_desc_t
   use clm_time_manager  , only : get_days_per_year, get_curr_date
   use clm_time_manager  , only : get_ref_date, timemgr_datediff 
   use spmdMod           , only : masterproc
   use decompMod         , only : get_proc_bounds, get_proc_clumps, get_clump_bounds
   use ColumnType        , only : col
   use LandunitType      , only : lun
   use landunit_varcon   , only : istsoil

   ! Used FATES Modules
   use FatesInterfaceMod     , only : fates_interface_type
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDTypesMod            , only : udata
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site
   use EDPftVarcon           , only : EDpftvarcon_inst
   use EDEcophysConType      , only : EDecophysconInit
   implicit none

   type, private :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:) 

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

      !



   end type f2hmap_type
   

   type, public :: hlm_fates_interface_type
      
      !      private
      

      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_interface_type), allocatable :: fates (:)
      

      ! This memory structure is used to map fates sites
      ! into the host model.  Currently, the FATES site
      ! and its column number matching are its only members

      type(f2hmap_type), allocatable  :: f2hmap(:)

      ! fates2hlm (previously called "clm_ed_inst") contains types and variables
      ! that are passed back to the driving land model, ie fates-to-hlm.  
      ! usefull to a calling model.  In this case HLM means "Hosting Land Model"
      ! prev:  type(ed_clm_type)::ed_clm_inst

      type(ed_clm_type) :: fates2hlm  

      
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
   
   subroutine init(this,bounds_proc, use_ed)
      
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
      logical,intent(in)                             :: use_ed   ! NEEDS TO BE PASSED (FOR NOW)
                                                                 ! BC THE FATES SITE VECTORS
                                                                 ! NEED TO BE GENERATED
                                                                 ! FOR NON-ED AS WELL.  SO
                                                                 ! ONLY PART OF THIS MAY BE OPERATIVE
      ! local variables
      integer                                        :: nclumps   ! Number of threads
      integer                                        :: nc        ! thread index
      integer                                        :: s         ! FATES site index
      integer                                        :: c         ! HLM column index
      integer, allocatable                           :: collist (:)
      type(bounds_type)                              :: bounds_clump

      if (use_ed) then
         
         ! Initialize the FATES communicators with the HLM
         ! This involves to stages
         ! 1) allocate the vectors
         ! 2) add the history variables defined in clm_inst to the history machinery
         call this%fates2hlm%Init(bounds_proc)
         
         ! Initialize ED phenology variables
         ! This also involves two stages
         ! 1) allocate the vectors in phen_inst
         ! 2) add the phenology history variables to the history machinery
         call this%phen_inst%Init(bounds_proc)
         
         
         call EDecophysconInit( EDpftvarcon_inst, numpft )

      end if
         

      nclumps = get_proc_clumps()
      allocate(this%fates(nclumps))
      allocate(this%f2hmap(nclumps))

      do nc = 1,nclumps
         
         call get_clump_bounds(nc, bounds_clump)
         nmaxcol = bounds_clump%endc - bounds_clump%begc + 1
         allocate(collist(1:nmaxcol))
         
         ! Allocate the mapping that points columns to FATES sites, 0 is NA
         allocate(self%f2hmap(nc)%hsites(bounds_clump%begc:bounds_clump%endc))

         ! Initialize all columns with a zero index, which indicates no FATES site
         self%f2hmap(nc)%hsites(:) = 0

         s = 0
         do c = bounds_clump%begc,bounds_clump%endc
            l = col%landunit(c)
               
            ! These are the key constraints that determine if this column
            ! will have a FATES site associated with it
            if (col%active(c) .and. lun%itype(l) == istsoil ) then 
               s = s + 1
               collist(s) = c
               self%f2hmap(nc)%hsites(c) = s
            endif
            
         enddo
         
         ! Allocate vectors that match FATES sites with HLM columns
         allocate(self%f2hmap(nc)%fcolumn(s))

         ! Assign the h2hmap indexing
         self%f2hmap(nc)%fcolumn(1:s)         =  collist(1:s)
         
         ! Deallocate the temporary arrays
         deallocate(collist)
         
         ! Set the number of FATES sites
         this%fates(nc)%nsites = s

         ! Allocate the FATES sites
         allocate (this%fates(nc)%sites(s))

!         call this%fates(nc)%init()
      end do

      
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
      
      call this%fates2hlm%ed_clm_link( bounds_clump,                  &
            this%fates(nc)%sites, &
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
      

      call this%fates2hlm%SetValues( bounds_clump, 0._r8 )

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
      do s = 1,this%fates(nc)%nsites

!         if (this%fates(nc)%sites(g)%istheresoil) then
            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates2hlm,                                  &
                  this%phen_inst, atm2lnd_inst,                    &
                  soilstate_inst, temperature_inst, waterstate_inst)
            
            call ed_update_site(this%fates(nc)%sites(s))
!         endif
      enddo

      ! link to CLM/ALM structures
      call this%fates2hlm%ed_clm_link( bounds_clump,                  &
            this%fates(nc)%sites,     &
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

   subroutine init_restart(this, ncid, flag, waterstate_inst, canopystate_inst )

      implicit none

      ! Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(file_desc_t)              , intent(inout) :: ncid    ! netcdf id
      character(len=*)               , intent(in)    :: flag    !'read' or 'write'
      type(waterstate_type)          , intent(in)    :: waterstate_inst
      type(canopystate_type)         , intent(in)    :: canopystate_inst

      ! Locals
      type(bounds_type) :: bounds_clump
      integer           :: nc
      integer           :: nclumps

      nclumps = get_proc_clumps()
      do nc = 1, nclumps
         call get_clump_bounds(nc, bounds_clump)

         call EDRest( bounds_clump, this%fates(nc)%sites, this%fates(nc)%nsites, &
               this%f2hmap(nc)%fcolumn, ncid, flag )

         if ( trim(flag) == 'read' ) then
            
            call this%fates2hlm%ed_clm_link( bounds_clump,                       &
                  this%fates(nc)%sites,                                          &
                  this%phen_inst,                                                &
                  waterstate_inst,                                               &
                  canopystate_inst)


         end if
      end do

      call clm_fates%fates2hlm%restart(bounds, ncid, flag)
      
      return
   end subroutine init_restart



   ! ------------------------------------------------------------------------------------
   !  THESE WRAPPERS MAY COME IN HANDY, KEEPING FOR NOW
   !  subroutine set_fates2hlm(this,bounds_clump, setval_scalar) 
   !
   !    ! This is a simple wrapper to flush some FATES -> CLM communicators
   !
   !    implicit none
   !    class(hlm_fates_interface_type), intent(inout) :: this
   !    type(bounds_type),intent(in)                :: bounds_clump
   !    real(r8),intent(in) :: setval_scalar      ! This value will flush all 
   !    
   !    call this%fates2hlm%SetValues( bounds_clump, setval_scalar )
   !  end subroutine set_fates2hlm
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
