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
   ! INTERF-TODO: NEED AN INVALID R8 SETTING FOR FATES
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
   use EnergyFluxType    , only : energyflux_type
   use SoilStateType     , only : soilstate_type
   use PhotosynthesisMod , only : photosyns_type
   use clm_varctl        , only : iulog, use_ed
   use clm_varcon        , only : tfrz
   use clm_varcon        , only : spval
   use clm_varpar        , only : numpft,            &
                                  numrad,            &
                                  nlevgrnd,          &
                                  nlevdecomp,        &
                                  nlevdecomp_full
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
   use clm_time_manager  , only : is_restart
   use ncdio_pio         , only : file_desc_t
   use clm_time_manager  , only : get_days_per_year, &
                                  get_curr_date
   use clm_time_manager  , only : get_ref_date,      &
                                  timemgr_datediff 
   use spmdMod           , only : masterproc
   use decompMod         , only : get_proc_bounds,   &
                                  get_proc_clumps,   &
                                  get_clump_bounds
   use GridCellType      , only : grc
   use ColumnType        , only : col
   use LandunitType      , only : lun
   use landunit_varcon   , only : istsoil
   use abortutils        , only : endrun
   use shr_log_mod       , only : errMsg => shr_log_errMsg    

   ! Used FATES Modules
   use FatesInterfaceMod     , only : fates_interface_type, &
                                      set_fates_ctrlparms,  &
                                      allocate_bcin,        &
                                      allocate_bcout
   
   use HistoryIOMod          , only : fates_hio_interface_type

   use EDCLMLinkMod          , only : ed_clm_type
   use EDTypesMod            , only : udata
   use EDTypesMod            , only : ed_patch_type
   use EDtypesMod            , only : numPatchesPerCol
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site
   use EDInitMod             , only : zero_site
   use EDInitMod             , only : init_patches
   use EDInitMod             , only : set_site_properties
   use EDPftVarcon           , only : EDpftvarcon_inst
   use EDEcophysConType      , only : EDecophysconInit
   use EDRestVectorMod       , only : EDRest
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs, ED_Norman_Radiation
   use EDBtranMod            , only : btran_ed, &
                                      get_active_suction_layers

   use EDPhotosynthesisMod   , only : Photosynthesis_ED
   use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
   use EDPhysiologyMod       , only: flux_into_litter_pools

   implicit none

   type, public :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:) 

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

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

      ! fates_hio is the interface class for the history output
      type(fates_hio_interface_type) :: fates_hio

   contains
      
      procedure, public :: init
      procedure, public :: init_allocate
      procedure, public :: check_hlm_active
      procedure, public :: init_restart
      procedure, public :: init_coldstart
      procedure, public :: dynamics_driv
      procedure, public :: wrap_sunfrac
      procedure, public :: wrap_btran
      procedure, public :: wrap_photosynthesis
      procedure, public :: wrap_accumulatefluxes
      procedure, public :: prep_canopyfluxes
      procedure, public :: wrap_canopy_radiation
      procedure, private :: wrap_litter_fluxout
      procedure, private :: init_history_io

   end type hlm_fates_interface_type


   logical :: DEBUG  = .false.

contains
   
   ! ====================================================================================

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

      if (use_ed) then
         
         ! Initialize the FATES communicators with the HLM
         ! This involves to stages
         ! 1) allocate the vectors
         ! 2) add the history variables defined in clm_inst to the history machinery
         call this%fates2hlm%Init(bounds_proc)
                  
         call EDecophysconInit( EDpftvarcon_inst, numpft )

      end if
         
      if(DEBUG)then
         write(iulog,*) 'Entering clm_fates%init'
      end if


      nclumps = get_proc_clumps()
      allocate(this%fates(nclumps))
      allocate(this%f2hmap(nclumps))


      ! ---------------------------------------------------------------------------------
      ! Send dimensions and other model controling parameters to FATES.  These
      ! are obviously only those parameters that are dictated by the host
      ! ---------------------------------------------------------------------------------
      
      ! Force FATES parameters that are recieve type, to the unset value
      call set_fates_ctrlparms('flush_to_unset')
      
      ! Send parameters individually
      call set_fates_ctrlparms('num_sw_bbands',ival=numrad)
      call set_fates_ctrlparms('num_lev_ground',ival=nlevgrnd)
      call set_fates_ctrlparms('num_levdecomp',ival=nlevdecomp)
      call set_fates_ctrlparms('num_levdecomp_full',ival=nlevdecomp_full)
      call set_fates_ctrlparms('hlm_name',cval='CLM')
      call set_fates_ctrlparms('hio_ignore_val',rval=spval)

      ! Check through FATES parameters to see if all have been set
      call set_fates_ctrlparms('check_allset')


      if(DEBUG)then
         write(iulog,*) 'clm_fates%init():  allocating for ',nclumps,' threads'
      end if

      return
   end subroutine init

   ! ====================================================================================
   
   subroutine init_allocate(this)
      
      implicit none
      
      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      ! local variables
      integer                                        :: nclumps   ! Number of threads
      integer                                        :: nc        ! thread index
      integer                                        :: s         ! FATES site index
      integer                                        :: c         ! HLM column index
      integer                                        :: l         ! HLM LU index
      integer, allocatable                           :: collist (:)
      type(bounds_type)                              :: bounds_clump
      type(bounds_type)                              :: bounds_proc
      integer                                        :: nmaxcol

      if(DEBUG)then
         write(iulog,*) 'Entering clm_fates%init_allocate'
      end if

      nclumps = get_proc_clumps()

      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,nmaxcol,s,c,l,collist)
      do nc = 1,nclumps
         
         call get_clump_bounds(nc, bounds_clump)
         nmaxcol = bounds_clump%endc - bounds_clump%begc + 1

         allocate(collist(1:nmaxcol))
         
         ! Allocate the mapping that points columns to FATES sites, 0 is NA
         allocate(this%f2hmap(nc)%hsites(bounds_clump%begc:bounds_clump%endc))

         ! Initialize all columns with a zero index, which indicates no FATES site
         this%f2hmap(nc)%hsites(:) = 0

         s = 0
         do c = bounds_clump%begc,bounds_clump%endc
            l = col%landunit(c)
               
            ! These are the key constraints that determine if this column
            ! will have a FATES site associated with it
            
            

            ! INTERF-TODO: WE HAVE NOT FILTERED OUT FATES SITES ON INACTIVE COLUMNS.. YET
            ! NEED A RUN-TIME ROUTINE THAT CLEARS AND REWRITES THE SITE LIST

            if (lun%itype(l) == istsoil ) then
               s = s + 1
               collist(s) = c
               this%f2hmap(nc)%hsites(c) = s
               if(DEBUG)then
                  write(iulog,*) 'clm_fates%init(): thread',nc,': found column',c,'with lu',l
                  write(iulog,*) 'LU type:', lun%itype(l)
               end if
            endif
            
         enddo

         if(DEBUG)then
            write(iulog,*) 'clm_fates%init(): thread',nc,': allocated ',s,' sites'
         end if

         ! Allocate vectors that match FATES sites with HLM columns
         ! RGK: Sites and fcolumns are forced as args during clm_driv() as of 6/4/2016
         ! We may have to give these a dummy allocation of 1, which should
         ! not be a problem since we always iterate on nsites.

         allocate(this%f2hmap(nc)%fcolumn(s))

         ! Assign the h2hmap indexing
         this%f2hmap(nc)%fcolumn(1:s)         =  collist(1:s)
         
         ! Deallocate the temporary arrays
         deallocate(collist)
         
         ! Set the number of FATES sites
         this%fates(nc)%nsites = s

         ! Allocate the FATES sites
         allocate (this%fates(nc)%sites(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (in)
         allocate(this%fates(nc)%bc_in(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (out)
         allocate(this%fates(nc)%bc_out(this%fates(nc)%nsites))

         ! Allocate and Initialize the Boundary Condition Arrays
         ! These are staticaly allocated at maximums, so
         ! No information about the patch or cohort
         ! structure is needed at this step
         
         do s = 1, this%fates(nc)%nsites
            call allocate_bcin(this%fates(nc)%bc_in(s))
            call allocate_bcout(this%fates(nc)%bc_out(s))
            call this%fates(nc)%zero_bcs(s)
         end do


         if( this%fates(nc)%nsites == 0 ) then
            write(iulog,*) 'Clump ',nc,' had no valid FATES sites'
            write(iulog,*) 'This will likely cause problems until code is improved'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

      end do
      !$OMP END PARALLEL DO
      
      call get_proc_bounds(bounds_proc)
      call this%init_history_io(bounds_proc)

   end subroutine init_allocate
   
  
   ! ------------------------------------------------------------------------------------
   

   subroutine check_hlm_active(this, nc, bounds_clump)

      
      implicit none
      class(hlm_fates_interface_type), intent(inout) :: this
      integer                                        :: nc
      type(bounds_type),intent(in)                   :: bounds_clump
      
      ! local variables
      integer :: c

      ! FATES-TODO: THIS SHOULD BE CHANGED TO DO RE-ALLOCATION
      ! INSTEAD OF FAILURE
      
      do c = bounds_clump%begc,bounds_clump%endc

         ! FATES ACTIVE BUT HLM IS NOT
         if(this%f2hmap(nc)%hsites(c)>0 .and. .not.col%active(c)) then

            
            write(iulog,*) 'INACTIVE COLUMN WITH ACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(msg=errMsg(__FILE__, __LINE__))

         elseif (this%f2hmap(nc)%hsites(c)==0 .and. col%active(c)) then
            
            write(iulog,*) 'ACTIVE COLUMN WITH INACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end do



   end subroutine check_hlm_active

   ! ------------------------------------------------------------------------------------

   subroutine dynamics_driv(this, nc, bounds_clump,      &
         atm2lnd_inst, soilstate_inst, temperature_inst, &
         waterstate_inst, canopystate_inst, soilbiogeochem_carbonflux_inst)
    
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
      type(soilbiogeochem_carbonflux_type), intent(out) :: soilbiogeochem_carbonflux_inst

      ! !LOCAL VARIABLES:
      real(r8) :: dayDiff                  ! day of run
      integer  :: dayDiffInt               ! integer of day of run
      integer  :: s                        ! site
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

            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates2hlm,                                  &
                  atm2lnd_inst,                                    &
                  soilstate_inst, temperature_inst, waterstate_inst)
            
            call ed_update_site(this%fates(nc)%sites(s))

      enddo

      call wrap_litter_fluxout(this, nc, bounds_clump, canopystate_inst, soilbiogeochem_carbonflux_inst)
      
      
      ! link to CLM/ALM structures
      call this%fates2hlm%ed_clm_link( bounds_clump,               &
            this%fates(nc)%sites,                                  &
            this%fates(nc)%nsites,                                 &
            this%f2hmap(nc)%fcolumn,                               &
            waterstate_inst,                                       &
            canopystate_inst)


      call this%fates_hio%update_history_dyn( nc, &
           this%fates(nc)%sites,       &
           this%fates(nc)%nsites,      &
           this%f2hmap(nc)%fcolumn)


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
      type(waterstate_type)          , intent(inout) :: waterstate_inst
      type(canopystate_type)         , intent(inout) :: canopystate_inst

      ! Locals
      type(bounds_type) :: bounds_clump
      integer           :: nc
      integer           :: nclumps

      nclumps = get_proc_clumps()
      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
      do nc = 1, nclumps
         if (this%fates(nc)%nsites>0) then
            call get_clump_bounds(nc, bounds_clump)
            
            call EDRest( bounds_clump,                                             &
                 this%fates(nc)%sites,                                             &
                 this%fates(nc)%nsites,                                            &
                 this%f2hmap(nc)%fcolumn, ncid, flag )
            
            if ( trim(flag) == 'read' ) then
               
               call this%fates2hlm%ed_clm_link( bounds_clump,                      &
                    this%fates(nc)%sites,                                          &
                    this%fates(nc)%nsites,                                         &
                    this%f2hmap(nc)%fcolumn,                                       &
                    waterstate_inst,                                               &
                    canopystate_inst)
               
               call this%fates_hio%update_history_dyn( nc, &
                    this%fates(nc)%sites,       &
                    this%fates(nc)%nsites,      &
                    this%f2hmap(nc)%fcolumn)

            end if
         end if
         call this%fates2hlm%restart(bounds_clump, ncid, flag)
      end do
      !$OMP END PARALLEL DO
      
      return
   end subroutine init_restart

   ! ====================================================================================

   subroutine init_coldstart(this, waterstate_inst, canopystate_inst)

     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     type(waterstate_type)          , intent(inout) :: waterstate_inst
     type(canopystate_type)         , intent(inout) :: canopystate_inst

     ! locals
     integer                                        :: nclumps
     integer                                        :: nc
     type(bounds_type)                              :: bounds_clump
     ! locals
     integer :: s
     integer :: c
     integer :: g


     nclumps = get_proc_clumps()

     !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s,c,g)
     do nc = 1, nclumps
        
        if ( this%fates(nc)%nsites>0 ) then

           call get_clump_bounds(nc, bounds_clump)

           do s = 1,this%fates(nc)%nsites
              call zero_site(this%fates(nc)%sites(s))
              c = this%f2hmap(nc)%fcolumn(s)
              g = col%gridcell(c)
              this%fates(nc)%sites(s)%lat = grc%latdeg(g)  
              this%fates(nc)%sites(s)%lon = grc%londeg(g)
           end do
           
           call set_site_properties(this%fates(nc)%sites, this%fates(nc)%nsites)
           
           call init_patches(this%fates(nc)%sites, this%fates(nc)%nsites)

           do s = 1,this%fates(nc)%nsites
              call ed_update_site(this%fates(nc)%sites(s))
           end do
           
           call this%fates2hlm%ed_clm_link( bounds_clump,           &
                this%fates(nc)%sites,                               &
                this%fates(nc)%nsites,                              &
                this%f2hmap(nc)%fcolumn,                            &
                waterstate_inst,                                    &
                canopystate_inst)

           call this%fates_hio%update_history_dyn( nc, &
                this%fates(nc)%sites,       &
                this%fates(nc)%nsites,      &
                this%f2hmap(nc)%fcolumn)


        end if
     end do
     !$OMP END PARALLEL DO
     return
   end subroutine init_coldstart

   ! ======================================================================================
   
   subroutine wrap_sunfrac(this,nc,atm2lnd_inst,canopystate_inst)
         
      
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      
      implicit none
      
      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      
      integer, intent(in)                  :: nc
      
      ! direct and diffuse downwelling radiation (W/m2)
      type(atm2lnd_type),intent(in)        :: atm2lnd_inst
      
      ! Input/Output Arguments to CLM
      type(canopystate_type),intent(inout) :: canopystate_inst
      
      ! Local Variables
      integer  :: p                           ! global index of the host patch
      integer  :: g                           ! global index of the host gridcell
      integer  :: c                           ! global index of the host column

      integer  :: s                           ! FATES site index
      integer  :: ifp                         ! FATEs patch index
                                              ! this is the order increment of patch
                                              ! on the site
      
      type(ed_patch_type), pointer :: cpatch  ! c"urrent" patch  INTERF-TODO: SHOULD
                                              ! BE HIDDEN AS A FATES PRIVATE
      
      associate( forc_solad => atm2lnd_inst%forc_solad_grc, &
                 forc_solai => atm2lnd_inst%forc_solai_grc, &
                 fsun       => canopystate_inst%fsun_patch, &
                 laisun     => canopystate_inst%laisun_patch, &               
                 laisha     => canopystate_inst%laisha_patch )

        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! The sun-shade calculations are performed only on FATES patches
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           g = col%gridcell(c)

           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
           !do ifp = 1, this%fates(nc)%bc_in(s)%npatches

              p = ifp+col%patchi(c)

              this%fates(nc)%bc_in(s)%solad_parb(ifp,:) = forc_solad(g,:)
              this%fates(nc)%bc_in(s)%solai_parb(ifp,:) = forc_solai(g,:)

           end do
        end do

        ! -------------------------------------------------------------------------------
        ! Call FATES public function to calculate internal sun/shade structures
        ! as well as total patch sun/shade fraction output boundary condition
        ! -------------------------------------------------------------------------------

        call ED_SunShadeFracs(this%fates(nc)%sites,  &
                              this%fates(nc)%nsites, &
                              this%fates(nc)%bc_in,  &
                              this%fates(nc)%bc_out)

        ! -------------------------------------------------------------------------------
        ! Transfer the FATES output boundary condition for canopy sun/shade fraction
        ! to the HLM
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
              p = ifp+col%patchi(c)
              fsun(p)   = this%fates(nc)%bc_out(s)%fsun_pa(ifp)
              laisun(p) = this%fates(nc)%bc_out(s)%laisun_pa(ifp)
              laisha(p) = this%fates(nc)%bc_out(s)%laisha_pa(ifp)
           end do
        end do

      end associate
      return
   end subroutine wrap_sunfrac
   
   ! ===================================================================================

   subroutine prep_canopyfluxes(this, nc, fn, filterp, photosyns_inst)

     ! ----------------------------------------------------------------------
     ! the main function for calculating photosynthesis is called within a
     ! loop based on convergence.  Some intitializations, including 
     ! canopy resistance must be intitialized before the loop
     ! ----------------------------------------------------------------------
    
     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     integer, intent(in)                            :: nc
     integer, intent(in)                            :: fn
     integer, intent(in)                            :: filterp(fn)
     type(photosyns_type),intent(inout)             :: photosyns_inst
     ! locals
     integer                                        :: f,p,c,s
     ! parameters
     integer,parameter                              :: rsmax0 = 2.e4_r8

     if (.not.use_ed) return
     
     do s = 1, this%fates(nc)%nsites
        ! filter flag == 1 means that this patch has not been called for photosynthesis
        this%fates(nc)%bc_in(s)%filter_photo_pa(:) = 1
     end do
  end subroutine prep_canopyfluxes

   ! ====================================================================================
   
   subroutine wrap_btran(this,nc,fn,filterc,soilstate_inst, waterstate_inst, &
                         temperature_inst, energyflux_inst,  &
                         soil_water_retention_curve)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates btran for FATES, this will be an input boundary
      ! condition for FATES photosynthesis/transpiration.
      !
      ! This subroutine also calculates rootr
      ! 
      ! ---------------------------------------------------------------------------------

      use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

      implicit none
      
      ! Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      integer                , intent(in)            :: nc
      integer                , intent(in)            :: fn
      integer                , intent(in)            :: filterc(fn) ! This is a list of
                                                                        ! columns with exposed veg
      type(soilstate_type)   , intent(inout)         :: soilstate_inst
      type(waterstate_type)  , intent(in)            :: waterstate_inst
      type(temperature_type) , intent(in)            :: temperature_inst
      type(energyflux_type)  , intent(inout)         :: energyflux_inst
      class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

      ! local variables
      real(r8) :: smp_node ! Soil suction potential, negative, [mm]
      real(r8) :: s_node
      integer  :: s
      integer  :: c
      integer  :: j
      integer  :: ifp
      integer  :: p
      
      associate(& 
         sucsat      => soilstate_inst%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm) 
         watsat      => soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         bsw         => soilstate_inst%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b" 
         eff_porosity => soilstate_inst%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice       
         t_soisno    => temperature_inst%t_soisno_col       , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         h2osoi_liqvol => waterstate_inst%h2osoi_liqvol_col , & ! Input: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR
         btran       => energyflux_inst%btran_patch         , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1) 
         btran2       => energyflux_inst%btran2_patch       , & ! Output: [real(r8) (:)   ]  
         rresis      => energyflux_inst%rresis_patch        , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd) 
         rootr       => soilstate_inst%rootr_patch          & ! Output: [real(r8) (:,:) ]  Fraction of water uptake in each layer
         )

        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! Critical step: a filter is being passed in that dictates which columns have
        ! exposed vegetation (above snow).  This is necessary, because various hydrologic
        ! variables like h2osoi_liqvol are not calculated and will have uninitialized
        ! values outside this list.
        !
        ! bc_in(s)%filter_btran      (this is in, but is also used in this subroutine)
        !
        ! We also filter a second time within this list by determining which soil layers
        ! have conditions for active uptake based on soil moisture and temperature. This
        ! must be determined by FATES (science stuff).  But the list of layers and patches
        ! needs to be passed back to the interface, because it then needs to request
        ! suction on these layers via CLM/ALM functions.  We cannot wide-swath calculate
        ! this on all layers, because values with no moisture or low temps will generate
        ! unstable values and cause sigtraps.
        ! -------------------------------------------------------------------------------
        
        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)

           ! Check to see if this column is in the exposed veg filter
           if( any(filterc==c) )then
              
              this%fates(nc)%bc_in(s)%filter_btran = .true.
              do j = 1,nlevgrnd
                 this%fates(nc)%bc_in(s)%tempk_gl(j)         = t_soisno(c,j)
                 this%fates(nc)%bc_in(s)%h2o_liqvol_gl(j)    = h2osoi_liqvol(c,j)
                 this%fates(nc)%bc_in(s)%eff_porosity_gl(j)  = eff_porosity(c,j)
                 this%fates(nc)%bc_in(s)%watsat_gl(j)        = watsat(c,j)
              end do

           else
              this%fates(nc)%bc_in(s)%filter_btran = .false.
              this%fates(nc)%bc_in(s)%tempk_gl(:)         = -999._r8
              this%fates(nc)%bc_in(s)%h2o_liqvol_gl(:)    = -999._r8
              this%fates(nc)%bc_in(s)%eff_porosity_gl(:)  = -999._r8
              this%fates(nc)%bc_in(s)%watsat_gl(:)        = -999._r8
           end if

        end do

        ! -------------------------------------------------------------------------------
        ! This function evaluates the ground layer to determine if
        ! root water uptake can happen, and soil suction should even
        ! be calculated.  We ask FATES for a boundary condition output
        ! logical because we don't want science calculations in the interface
        ! yet... hydrology (suction calculation) is provided by the host
        ! so we need fates to tell us where to calculate suction
        ! but not calculate it itself. Yeah, complicated, but thats life.
        ! -------------------------------------------------------------------------------
        call get_active_suction_layers(this%fates(nc)%sites,  &
                                       this%fates(nc)%nsites, &
                                       this%fates(nc)%bc_in,  &
                                       this%fates(nc)%bc_out)

        ! Now that the active layers of water uptake have been decided by fates
        ! Calculate the suction that is passed back to fates
        ! Note that the filter_btran is unioned with active_suction_gl

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           do j = 1,nlevgrnd
              if(this%fates(nc)%bc_out(s)%active_suction_gl(j)) then
                 s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)
                 call soil_water_retention_curve%soil_suction(c,j,s_node, soilstate_inst, smp_node)
                 this%fates(nc)%bc_in(s)%smp_gl(j)           = smp_node
              end if
           end do
        end do
        
        ! -------------------------------------------------------------------------------
        ! Suction and active uptake layers calculated, lets calculate uptake (btran)
        ! This will calculate internals, as well as output boundary conditions: 
        ! btran, rootr
        ! -------------------------------------------------------------------------------

        call btran_ed(this%fates(nc)%sites,  &
                      this%fates(nc)%nsites, &
                      this%fates(nc)%bc_in,  &
                      this%fates(nc)%bc_out)


        ! -------------------------------------------------------------------------------
        ! Convert output BC's
        ! For CLM/ALM this wrapper provides return variables that should
        ! be similar to that of calc_root_moist_stress().  However,
        ! CLM/ALM-FATES simulations will no make use of rresis, btran or btran2
        ! outside of FATES. We do not have code in place to calculate btran2 or
        ! rresis right now, so we force to bad.  We have btran calculated so we
        ! pass it in case people want diagnostics.  rootr is actually the only
        ! variable that will be used, as it is needed to help distribute the
        ! the transpiration sink to the appropriate layers. (RGK)
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
              
              p = ifp+col%patchi(c)
              
              do j = 1,nlevgrnd
                 
                 rresis(p,j) = -999.9  ! We do not calculate this correctly
                 ! it should not thought of as valid output until we decide to.
                 rootr(p,j)  = this%fates(nc)%bc_out(s)%rootr_pagl(ifp,j)
                 btran(p)    = this%fates(nc)%bc_out(s)%btran_pa(ifp)
                 btran2(p)   = -999.9  ! Not available, force to nonsense
                 
              end do
           end do
        end do
      end associate
      return
   end subroutine wrap_btran

   ! ====================================================================================
   
   subroutine wrap_photosynthesis(this, nc, bounds, fn, filterp, &
         esat_tv, eair, oair, cair, rb, dayl_factor,             &
         atm2lnd_inst, temperature_inst, canopystate_inst, photosyns_inst)
   
    use shr_log_mod       , only : errMsg => shr_log_errMsg
    use abortutils        , only : endrun
    use decompMod         , only : bounds_type
    use clm_time_manager  , only : get_step_size
    use clm_varcon        , only : rgas, tfrz, namep  
    use clm_varpar        , only : nlevsoi, mxpft
    use clm_varctl        , only : iulog
    use pftconMod         , only : pftcon
    use perf_mod          , only : t_startf, t_stopf
    use PatchType         , only : patch
    use quadraticMod      , only : quadratic
    use EDParamsMod       , only : ED_val_grperc
    use EDSharedParamsMod , only : EDParamsShareInst
    use EDTypesMod        , only : numpft_ed, dinc_ed
    use EDtypesMod        , only : ed_patch_type, ed_cohort_type, ed_site_type, numpft_ed, numPatchesPerCol
    use EDEcophysContype  , only : EDecophyscon
   
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    integer                , intent(in)            :: nc                          ! clump index
    type(bounds_type)      , intent(in)            :: bounds
    integer                , intent(in)            :: fn                          ! size of pft filter
    integer                , intent(in)            :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)            :: esat_tv(bounds%begp: )      ! saturation vapor pressure at t_veg (Pa)
    real(r8)               , intent(in)            :: eair( bounds%begp: )        ! vapor pressure of canopy air (Pa)
    real(r8)               , intent(in)            :: oair( bounds%begp: )        ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)            :: cair( bounds%begp: )        ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)            :: rb( bounds%begp: )          ! boundary layer resistance (s/m)
    real(r8)               , intent(in)            :: dayl_factor( bounds%begp: ) ! scalar (0-1) for daylength
    type(atm2lnd_type)     , intent(in)            :: atm2lnd_inst
    type(temperature_type) , intent(in)            :: temperature_inst
    type(canopystate_type) , intent(inout)         :: canopystate_inst
    type(photosyns_type)   , intent(inout)         :: photosyns_inst

    integer                                        :: s,c,p,ifp,j,icp
    real(r8)                                       :: dtime

    call t_startf('edpsn')
    associate(&
          t_soisno  => temperature_inst%t_soisno_col , &
          t_veg     => temperature_inst%t_veg_patch  , &
          tgcm      => temperature_inst%thm_patch    , &
          forc_pbot => atm2lnd_inst%forc_pbot_downscaled_col, &
          rssun     => photosyns_inst%rssun_patch  , &
          rssha     => photosyns_inst%rssha_patch,   &
          psnsun    => photosyns_inst%psnsun_patch,  &
          psnsha    => photosyns_inst%psnsha_patch)
      
      do s = 1, this%fates(nc)%nsites
         
         c = this%f2hmap(nc)%fcolumn(s)
         
         do j = 1,nlevsoi
            this%fates(nc)%bc_in(s)%t_soisno_gl(j)   = t_soisno(c,j)  ! soil temperature (Kelvin)
         end do
         this%fates(nc)%bc_in(s)%forc_pbot           = forc_pbot(c)   ! atmospheric pressure (Pa)

         do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
            
            p = ifp+col%patchi(c)

            ! Check to see if this patch is in the filter
            ! Note that this filter is most likely changing size, and getting smaller
            ! and smaller as more patch have converged on solution
            if( any(filterp==p) )then

               ! This filter is flushed to 1 before the canopyflux stability iterator
               ! It is set to status 2 if it is an active patch within the iterative loop
               ! After photosynthesis is called, it is upgraded to 3 if it was called.
               ! After all iterations we can evaluate which patches have a final flag
               ! of 3 to check if we missed any.

               this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 2
               this%fates(nc)%bc_in(s)%dayl_factor_pa(ifp) = dayl_factor(p) ! scalar (0-1) for daylength
               this%fates(nc)%bc_in(s)%esat_tv_pa(ifp)     = esat_tv(p)     ! saturation vapor pressure at t_veg (Pa)
               this%fates(nc)%bc_in(s)%eair_pa(ifp)        = eair(p)        ! vapor pressure of canopy air (Pa)
               this%fates(nc)%bc_in(s)%oair_pa(ifp)        = oair(p)        ! Atmospheric O2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%cair_pa(ifp)        = cair(p)        ! Atmospheric CO2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%rb_pa(ifp)          = rb(p)          ! boundary layer resistance (s/m)
               this%fates(nc)%bc_in(s)%t_veg_pa(ifp)       = t_veg(p)       ! vegetation temperature (Kelvin)     
               this%fates(nc)%bc_in(s)%tgcm_pa(ifp)        = tgcm(p)        ! air temperature at agcm reference height (kelvin)
            end if
         end do
      end do

      dtime = get_step_size()
      
      ! Call photosynthesis
      
      call  Photosynthesis_ED (this%fates(nc)%sites,  &
                               this%fates(nc)%nsites, &
                               this%fates(nc)%bc_in,  &
                               this%fates(nc)%bc_out, &
                               dtime)

      ! Perform a double check to see if all patches on naturally vegetated columns
      ! were activated for photosynthesis
      ! ---------------------------------------------------------------------------------
      do icp = 1,fn
         p = filterp(icp)
         c = patch%column(p)
         s = this%f2hmap(nc)%hsites(c)
         ! do if structure here and only pass natveg columns
         ifp = p-col%patchi(c)
         if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 2)then
            write(iulog,*) 'Not all patches on the natveg column in the photosynthesis'
            write(iulog,*) 'filter ran photosynthesis'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         else
            this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 3
            rssun(p) = this%fates(nc)%bc_out(s)%rssun_pa(ifp)
            rssha(p) = this%fates(nc)%bc_out(s)%rssha_pa(ifp)
            
            ! These fields are marked with a bad-value flag
            photosyns_inst%psnsun_patch(p)   = spval
            photosyns_inst%psnsha_patch(p)   = spval
         end if
      end do
      
    end associate
    call t_stopf('edpsn')
    return
 end subroutine wrap_photosynthesis

 ! ======================================================================================

 subroutine wrap_accumulatefluxes(this, nc, fn, filterp)

   use clm_time_manager  , only : get_step_size
   
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   integer                , intent(in)            :: nc                   ! clump index
   integer                , intent(in)            :: fn                   ! size of pft filter
   integer                , intent(in)            :: filterp(fn)          ! pft filter
   
   integer                                        :: s,c,p,ifp,icp
   real(r8)                                       :: dtime

    ! Run a check on the filter
    do icp = 1,fn
       p = filterp(icp)
       c = patch%column(p)
       s = this%f2hmap(nc)%hsites(c)
       ifp = p-col%patchi(c)
       if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 3)then
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do

    dtime = get_step_size()
    call  AccumulateFluxes_ED(this%fates(nc)%sites,  &
                               this%fates(nc)%nsites, &
                               this%fates(nc)%bc_in,  &
                               this%fates(nc)%bc_out, &
                               dtime)

    
    call this%fates_hio%update_history_prod(nc, &
                               this%fates(nc)%sites,  &
                               this%fates(nc)%nsites, &
                               this%f2hmap(nc)%fcolumn, &
                               dtime)

    return

 end subroutine wrap_accumulatefluxes

 ! ======================================================================================

 subroutine wrap_canopy_radiation(this, bounds_clump, nc, &
         filter_vegsol, num_vegsol, coszen, surfalb_inst)


    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)             :: bounds_clump
    ! filter for vegetated pfts with coszen>0
    integer            , intent(in)            :: nc ! clump index
    integer            , intent(in)            :: filter_vegsol(num_vegsol)    
    integer            , intent(in)            :: num_vegsol                 
    ! cosine solar zenith angle for next time step
    real(r8)           , intent(in)            :: coszen( bounds_clump%begp: )        
    type(surfalb_type) , intent(inout)         :: surfalb_inst 
    
    ! locals
    integer                                    :: s,c,p,ifp,icp

    associate(&
         albgrd_col   =>    surfalb_inst%albgrd_col         , & !in
         albgri_col   =>    surfalb_inst%albgri_col         , & !in
         albd         =>    surfalb_inst%albd_patch         , & !out
         albi         =>    surfalb_inst%albi_patch         , & !out
         fabd         =>    surfalb_inst%fabd_patch         , & !out
         fabi         =>    surfalb_inst%fabi_patch         , & !out
         ftdd         =>    surfalb_inst%ftdd_patch         , & !out
         ftid         =>    surfalb_inst%ftid_patch         , & !out
         ftii         =>    surfalb_inst%ftii_patch)            !out

    do s = 1, this%fates(nc)%nsites

       c = this%f2hmap(nc)%fcolumn(s)
       do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
          
          p = ifp+col%patchi(c)
          
          if( any(filter_vegsol==p) )then
    
             this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) = .true.
             this%fates(nc)%bc_in(s)%coszen_pa(ifp)  = coszen(p)
             this%fates(nc)%bc_in(s)%albgr_dir_rb(:) = albgrd_col(c,:)
             this%fates(nc)%bc_in(s)%albgr_dif_rb(:) = albgri_col(c,:)

          else
             
             this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) = .false.

          end if

       end do
    end do

    call ED_Norman_Radiation(this%fates(nc)%sites,  &
                             this%fates(nc)%nsites, &
                             this%fates(nc)%bc_in,  &
                             this%fates(nc)%bc_out)
    
    ! Pass FATES BC's back to HLM
    ! -----------------------------------------------------------------------------------
    do icp = 1,num_vegsol
       p = filter_vegsol(icp)
       c = patch%column(p)
       s = this%f2hmap(nc)%hsites(c)
       ! do if structure here and only pass natveg columns
       ifp = p-col%patchi(c)

       if(.not.this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) )then
          write(iulog,*) 'Not all patches on the natveg column were passed to canrad'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       else
          albd(p,:) = this%fates(nc)%bc_out(s)%albd_parb(ifp,:)
          albi(p,:) = this%fates(nc)%bc_out(s)%albi_parb(ifp,:)
          fabd(p,:) = this%fates(nc)%bc_out(s)%fabd_parb(ifp,:)
          fabi(p,:) = this%fates(nc)%bc_out(s)%fabi_parb(ifp,:)
          ftdd(p,:) = this%fates(nc)%bc_out(s)%ftdd_parb(ifp,:)
          ftid(p,:) = this%fates(nc)%bc_out(s)%ftid_parb(ifp,:)
          ftii(p,:) = this%fates(nc)%bc_out(s)%ftii_parb(ifp,:)
       end if
    end do
    
  end associate

    return
 end subroutine wrap_canopy_radiation

 ! ======================================================================================
 
 subroutine wrap_litter_fluxout(this, nc, bounds_clump, canopystate_inst, soilbiogeochem_carbonflux_inst)
     
    implicit none
    
    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    integer                , intent(in)            :: nc
    type(bounds_type),intent(in)                   :: bounds_clump
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    type(soilbiogeochem_carbonflux_type), intent(out) :: soilbiogeochem_carbonflux_inst
    
    ! local variables
    integer :: s, c
    
    
    ! process needed input boundary conditions to define rooting profiles
    ! call subroutine to aggregate ED litter output fluxes and package them for handing across interface
    ! process output into the dimensions that the BGC model wants (column, depth, and litter fractions)
    
    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = canopystate_inst%altmax_lastyear_indx_col(c)
    end do
    
    call flux_into_litter_pools(this%fates(nc)%sites,  &
                                this%fates(nc)%nsites, &
                                this%fates(nc)%bc_in,  &
                                this%fates(nc)%bc_out)
    
    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       soilbiogeochem_carbonflux_inst%FATES_c_to_litr_lab_c_col(c,:) = this%fates(nc)%bc_out(s)%FATES_c_to_litr_lab_c_col(:)
       soilbiogeochem_carbonflux_inst%FATES_c_to_litr_cel_c_col(c,:) = this%fates(nc)%bc_out(s)%FATES_c_to_litr_cel_c_col(:)
       soilbiogeochem_carbonflux_inst%FATES_c_to_litr_lig_c_col(c,:) = this%fates(nc)%bc_out(s)%FATES_c_to_litr_lig_c_col(:)
    end do


 end subroutine wrap_litter_fluxout

 ! ======================================================================================

 subroutine init_history_io(this,bounds_proc)

   use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
   use EDtypesMod , only : nlevsclass_ed
   use clm_varpar , only : mxpft, nlevgrnd

   ! Arguments
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type),intent(in)                   :: bounds_proc  ! Currently "proc"
   
   
   ! Locals
   type(bounds_type)                              :: bounds_clump
   integer :: nvar  ! number of IO variables found
   integer :: ivar  ! variable index 1:nvar
   integer :: nc    ! thread counter 1:nclumps
   integer :: nclumps ! number of threads on this proc
   integer :: s     ! FATES site index
   integer :: c     ! ALM/CLM column index
   character(len=32) :: dim2name
   
   ! This routine initializes the types of output variables
   ! not the variables themselves, just the types
   ! ---------------------------------------------------------------------------------

   if(.not.use_ed) return
   
   !associate(hio => this%fates_hio)
   
   nclumps = get_proc_clumps()

   ! ------------------------------------------------------------------------------------
   ! PART I: Set FATES DIMENSIONING INFORMATION
   !       
   ! -------------------------------------------------------------------------------
   ! Those who wish add variables that require new dimensions, please
   ! see FATES: HistoryIOMod.F90.  Dimension types are defined at the top of the
   ! module, and a new explicitly named instance of that type should be created.
   ! With this new dimension, a new output type/kind can contain that dimension.
   ! A new type/kind can be added to the iovar_dk structure, which defines its members
   ! in created in init_iovar_dk_maps().  Make sure to increase the size of n_iovar_dk.
   ! A type/kind of output is defined by the data type (ie r8,int,..)
   ! and the dimensions.  Keep in mind that 3D variables (or 4D if you include time)
   ! are not really supported in CLM/ALM right now.  There are ways around this
   ! limitations by creating combined dimensions, for instance the size+pft dimension
   ! "scpf"
   ! ------------------------------------------------------------------------------------
   
   call this%fates_hio%dim_init(this%fates_hio%iopa_dim,'patch',nclumps,bounds_proc%begp,bounds_proc%endp)
   call this%fates_hio%dim_init(this%fates_hio%iosi_dim,'column',nclumps,bounds_proc%begc,bounds_proc%endc)
   call this%fates_hio%dim_init(this%fates_hio%iogrnd_dim,'levgrnd',nclumps,1,nlevgrnd)
   call this%fates_hio%dim_init(this%fates_hio%ioscpf_dim,'levscpf',nclumps,1,nlevsclass_ed*mxpft)

   ! Allocate the mapping between FATES indices and the IO indices
   allocate(this%fates_hio%iovar_map(nclumps))

   
   ! Define the bounds on the first dimension for each thread
   !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s,c)
   do nc = 1,nclumps
      
      call get_clump_bounds(nc, bounds_clump)
      
      ! thread bounds for patch
      call this%fates_hio%set_dim_thread_bounds(this%fates_hio%iopa_dim,nc,bounds_clump%begp,bounds_clump%endp)
      call this%fates_hio%set_dim_thread_bounds(this%fates_hio%iosi_dim,nc,bounds_clump%begc,bounds_clump%endc)
      call this%fates_hio%set_dim_thread_bounds(this%fates_hio%iogrnd_dim,nc,1,nlevgrnd)
      call this%fates_hio%set_dim_thread_bounds(this%fates_hio%ioscpf_dim,nc,1,nlevsclass_ed*mxpft)

      ! ------------------------------------------------------------------------------------
      ! PART I.5: SET SOME INDEX MAPPINGS SPECIFICALLY FOR SITE<->COLUMN AND PATCH 
      ! ------------------------------------------------------------------------------------

      allocate(this%fates_hio%iovar_map(nc)%site_index(this%fates(nc)%nsites))
      allocate(this%fates_hio%iovar_map(nc)%patch1_index(this%fates(nc)%nsites))
      
      do s=1,this%fates(nc)%nsites
         c = this%f2hmap(nc)%fcolumn(s)
         this%fates_hio%iovar_map(nc)%site_index(s)   = c
         this%fates_hio%iovar_map(nc)%patch1_index(s) = col%patchi(c)+1
      end do
      
   end do
   !$OMP END PARALLEL DO
   
   ! ------------------------------------------------------------------------------------
   ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
   ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
   ! ------------------------------------------------------------------------------------

   call this%fates_hio%init_iovar_dk_maps()
   
   call this%fates_hio%set_dim_ptrs(dk_name='PA_R8',idim=1,dim_target=this%fates_hio%iopa_dim)

   call this%fates_hio%set_dim_ptrs(dk_name='SI_R8',idim=1,dim_target=this%fates_hio%iosi_dim)

   call this%fates_hio%set_dim_ptrs(dk_name='PA_GRND_R8',idim=1,dim_target=this%fates_hio%iopa_dim)
   call this%fates_hio%set_dim_ptrs(dk_name='PA_GRND_R8',idim=2,dim_target=this%fates_hio%iogrnd_dim)

   call this%fates_hio%set_dim_ptrs(dk_name='SI_GRND_R8',idim=1,dim_target=this%fates_hio%iosi_dim)
   call this%fates_hio%set_dim_ptrs(dk_name='SI_GRND_R8',idim=2,dim_target=this%fates_hio%iogrnd_dim)

   call this%fates_hio%set_dim_ptrs(dk_name='PA_SCPF_R8',idim=1,dim_target=this%fates_hio%iopa_dim)
   call this%fates_hio%set_dim_ptrs(dk_name='PA_SCPF_R8',idim=2,dim_target=this%fates_hio%ioscpf_dim)

   call this%fates_hio%set_dim_ptrs(dk_name='SI_SCPF_R8',idim=1,dim_target=this%fates_hio%iosi_dim)
   call this%fates_hio%set_dim_ptrs(dk_name='SI_SCPF_R8',idim=2,dim_target=this%fates_hio%ioscpf_dim)

   
   ! ------------------------------------------------------------------------------------
   ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
   ! HLM ACCORDING TO THEIR TYPES
   ! ------------------------------------------------------------------------------------
   
   ! Determine how many of the history IO variables registered in FATES
   ! are going to be allocated
   
   call this%fates_hio%define_history_vars('count',nvar)
   this%fates_hio%n_hvars = nvar
   
   ! Allocate the list of history output variable objects
   allocate(this%fates_hio%hvars(nvar))
   
   ! construct the object that defines all of the IO variables
   call this%fates_hio%define_history_vars('initialize')
   
   do ivar = 1,nvar
      
      associate( vname    => this%fates_hio%hvars(ivar)%vname, &
                 vunits   => this%fates_hio%hvars(ivar)%units,   &
                 vlong    => this%fates_hio%hvars(ivar)%long, &
                 vdefault => this%fates_hio%hvars(ivar)%use_default, &
                 vavgflag => this%fates_hio%hvars(ivar)%avgflag,  &
                 ioname   => this%fates_hio%hvars(ivar)%iovar_dk_ptr%name )
        
 
        
        select case(trim(ioname))
        case('PA_R8')
           call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_patch=this%fates_hio%hvars(ivar)%r81d,    &
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)
           
        case('SI_R8')
           call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=this%fates_hio%hvars(ivar)%r81d,      & 
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)

        case('PA_GRND_R8')
           dim2name = this%fates_hio%hvars(ivar)%iovar_dk_ptr%dim2_ptr%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         & ! <--- addfld2d
                              type2d=trim(dim2name),                        & ! <--- type2d
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_patch=this%fates_hio%hvars(ivar)%r82d,    & 
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)
           
        case('PA_SCPF_R8')
           dim2name = this%fates_hio%hvars(ivar)%iovar_dk_ptr%dim2_ptr%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
                              type2d=trim(dim2name),                        &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_patch=this%fates_hio%hvars(ivar)%r82d,    & 
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)
        case('SI_GRND_R8')
           dim2name = this%fates_hio%hvars(ivar)%iovar_dk_ptr%dim2_ptr%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
                              type2d=trim(dim2name),                        &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=this%fates_hio%hvars(ivar)%r82d,      & 
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)
        case('SI_SCPF_R8')
           dim2name = this%fates_hio%hvars(ivar)%iovar_dk_ptr%dim2_ptr%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
                              type2d=trim(dim2name),                        &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=this%fates_hio%hvars(ivar)%r82d,      & 
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)

        case default
           write(iulog,*) 'A FATES iotype was created that was not registerred'
           write(iulog,*) 'in CLM.:',trim(ioname)
           call endrun(msg=errMsg(__FILE__, __LINE__))
        end select
          
      end associate
   end do
   return
 end subroutine init_history_io


end module CLMFatesInterfaceMod
