module clmed_interfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the ED library with the CLM model driver.  All connections between the 
   ! two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see ED memory structures.
   ! The routines here, that call ED library routines, will only pass CLM arguments as
   ! either native type arrays (int,real,log, etc) or packed into ED boundary condition
   ! structures.
   !
   ! Conventions:
   ! all communicator functions start with CLMEDInterf_
   ! keep line widths within 90 spaces
   ! host model data structures are passed as simple vectors when possible with the
   ! intention that this script will be more portable to a similar host model.  Although
   ! this current method does break with the current convention of calling subroutines
   ! from clm_driv(), where data structure groups are passed, "ie waterstate_type".
   !
   ! -------------------------------------------------------------------------------------

   !  use ed_driver_interface, only: 
   
   ! Used CLM Modules
   use PatchType         , only : patch
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use decompMod         , only : bounds_type
   use WaterStateType    , only : waterstate_type
   use CanopyStateType   , only : canopystate_type
   use clm_varctl        , only : iulog
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use clm_time_manager  , only : is_restart
   
   ! Used ED Modules
   use EDtypesMod            , only : ed_patch_type, ed_site_type, numpft_ed
   use EDtypesMod            , only : map_clmpatch_to_edpatch
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDInitMod             , only : ed_init
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs
   
   implicit none
   
   

   type, public :: lsm_ed_interface_type

      private

      ! Previously many of these sub-classes had "ed_" naming conventions
      ! These prefixes are now implied as part of the lsm_ed_interface_type
      ! and are no longer necessary

      ! This is the root of the ED hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated
      ! with a grid-cell, this is intended to be migrated to columns

      ! prev:  type(ed_site_type)::ed_allsites_inst
      type(ed_site_type), allocatable :: site_inst(:)

      
      ! These are the communicator variables that are populated by ED, and are
      ! usefull to a calling model.  In this case DLM means "Driving Land Model"
      
      ! prev:  type(ed_clm_type)::ed_clm_inst
      type(dlm_type)    :: dlm_inst   

      
      ! These are phenology relevant variables (strange how phenology gets
      ! its own subclass)

      ! prev: type(ed_phenology_type)::ed_phenology_inst
      type(phen_type)   :: phen_inst  

   contains

      procedure, public :: Init
      procedure, public :: AllocateEDSites
      procedure, public :: CanopySunShadeFracs

   end type lsm_ed_interface_type


   ! FOR CLM: type lsm_ed_interface_type describes "clm_ed", which is defined in clm_instMod
   ! FOR ALM: NA

contains

   
   subroutine Init(this,bounds)

      ! ---------------------------------------------------------------------------------
      ! This initializes the lsm_ed_interface_type 
      !
      ! ed_allsites_inst is the root of the ED state hierarchy (instantaneous info on 
      ! the state of the ecosystem).  As such, it governs the connection points between
      ! the host (which also dictates its allocation) and its patch structures.
      !
      ! ed_allsites_inst may associate with different scales in different models. In
      ! CLM, it is being designed to relate to column scale.
      !
      ! This global may become relegated to this module. 
      !
      ! Note: CLM/ALM currently wants ed_allsites_inst to be allocated even if ed
      ! is not turned on
      ! ---------------------------------------------------------------------------------

      implicit none

      ! Input Arguments
      class(lsm_ed_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                :: bounds 

      
      ! Initialize the mapping elements between ED and the DLM

      ! ed_bounds = func(bounds)

      ! These bounds are for a single clump (thread)
      allocate (this%ed_site_inst(bounds%begg:bounds%endg))

      ! Initialize the ed_communicators with the DLM
      ! This involves to stages
      ! 1) allocate the vectors dlm_inst
      ! 2) add the history variables defined in clm_inst to the history machinery
      this%dlm_inst%Init(bounds)
      
      ! Initialize ED phenology variables
      ! This also involves two stages
      ! 1) allocate the vectors in phen_inst
      ! 2) add the phenology history variables to the history machinery
      this%phen_inst%Init(bounds)


   end subroutine Init


   subroutine StateInit(self,bounds_clump)

      
      ! CLM:  called from initialize2()
      ! ALM:  ??

      ! Input Arguments
      class(lsm_ed_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                :: bounds_clump

      ! locals
      integer :: g

      ! It is assumed that bounds is bounds_clump

      if ( .not. is_restart() ) then

         ! Initialize  (THIS ROUTINE CALLS CLM STUFF-MIGRATE CODE TO HERE)
         call ed_init_sites( bounds_clump, self%sites_inst(bounds_clump%begg:bounds_clump%endg) )

         do g = bounds%begg,bounds%endg
            if (this%sites_inst(g)%istheresoil) then
               call ed_update_site(this%sites_inst(g))
            end if
         end do
      endif
      


   end subroutine StateInit


   subroutine DLMInit(self,bounds,waterstate_inst, canopystate_inst)

      ! CLM:  called from initialize2()
      ! ALM:  ??

      ! Input Arguments
      class(lsm_ed_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                :: bounds_clump
      type(waterstate_type)   , intent(inout)     :: waterstate_inst
      type(canopystate_type)  , intent(inout)     :: canopystate_inst

      if ( .not. is_restart() ) then

       self%dlm_link%ed_clm_link( bounds,                                     &
                                  self%sites_inst(bounds%begg:bounds%endg),   &
                                  self%phen_inst,                             &
                                  waterstate_inst,                            &
                                  canopystate_inst)
    endif

   subroutine AllocateEDSites(this,bounds)
   


      
      return
   end subroutine AllocateEDSites

   ! -------------------------------------------------------------------------------------
    
   !CLMEDInterf_ed_init
   subroutine InitializeED(this,bounds, ed_clm_inst, ed_phenology_inst,  &
                                  waterstate_inst, canopystate_inst   )

      ! ----------------------------------------------------------------------------------
      ! This subroutine is simply a wrapper right now that calls ed_init
      ! As we start to pull CLM stuff out of ed_init, that CLM stuff will be moved into 
      ! this routine and ED things may stay in the ed_init routine or likewise be moved 
      ! here.  Currently, the only thing this subroutine does is wrap in 
      ! ed_allsites_inst as an ED global
      ! ----------------------------------------------------------------------------------
      
      implicit none
      ! !ARGUMENTS    
      class(lsm_ed_interface_type), intent(inout)     :: this
      type(bounds_type)       , intent(in)            :: bounds  ! clump bounds
      type(ed_clm_type)       , intent(inout)         :: ed_clm_inst
      type(ed_phenology_type) , intent(inout)         :: ed_phenology_inst
      type(waterstate_type)   , intent(inout)         :: waterstate_inst
      type(canopystate_type)  , intent(inout)         :: canopystate_inst
      
      ! Note bounds and the two "state" structures should not pass this routine
      ! but that will be staged for a later PR (RGK)

      call ed_init( bounds, ed_allsites_inst(bounds%begg:bounds%endg), ed_clm_inst, &
            ed_phenology_inst, waterstate_inst, canopystate_inst)


      


      return
   end subroutine InitializeED

   ! -------------------------------------------------------------------------------------

   subroutine InitEDSubTypes(this,bounds)

      class(lsm_ed_interface_type), intent(inout)     :: this


      this%phen_inst%initAccVars(bounds)



   end subroutine InitEDTypes


  
   subroutine CanopySunShadeFracs(this,filter_nourbanp, num_nourbanp, &
                                              atm2lnd_inst,canopystate_inst)


      ! ----------------------------------------------------------------------------------
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      ! ----------------------------------------------------------------------------------

      implicit none

      ! Input Arguments
      class(lsm_ed_interface_type), intent(inout) :: this
      
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
              cpatch => map_clmpatch_to_edpatch(this%ed_allsites_inst(g), p) 
              
              call ED_SunShadeFracs(cpatch,forc_solad(g,ipar),forc_solai(g,ipar),fsun(p))

           endif
           
        end do
      end associate
      return
   end subroutine CanopySunShadeFracs





end module clmed_interfaceMod
