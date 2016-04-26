module clmfates_interfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library with the (CLM)/ALM/ATS/etc model driver.  
   ! All connections between the two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see FATES memory 
   ! structures.
   ! The routines here, that call FATES library routines, will only pass CLM arguments as
   ! either native type arrays (int,real,log, etc) or packed into ED boundary condition
   ! structures.
   !
   ! Conventions:
   ! keep line widths within 90 spaces
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
   
   ! Used ED Modules
   use EDtypesMod            , only : ed_patch_type, ed_site_type, numpft_ed
   use EDtypesMod            , only : map_clmpatch_to_edpatch
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs
   
   implicit none

   type, public :: dlm_fates_interface_type
      
      private
      
      ! Change in naming conventions:
      ! Previously, many of these sub-classes had "ed_" naming conventions
      ! However, 
      ! 1) the model has gone through a name change and is now called
      ! FATES (however, FATES will continue to use core ED (Moorcroft et al. 2001)
      ! mechanics.
      ! 2) the submodules and memory structurs do not need to have the name
      ! of the FATES model, as it is now implied as part of this interface.


      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated
      ! with a grid-cell, this is intended to be migrated to columns
      ! prev:  type(ed_site_type)::ed_allsites_inst
      type(ed_site_type), allocatable :: site_inst(:)

      
      ! These are the communicator variables that are populated by ED/FATES, and are
      ! usefull to a calling model.  In this case DLM means "Driving Land Model"
      ! prev:  type(ed_clm_type)::ed_clm_inst
      type(ed_clm_type) :: fates2dlm_inst   

      
      ! These are phenology relevant variables (strange how phenology gets
      ! its own subclass)
      ! prev: type(ed_phenology_type)::ed_phenology_inst
      type(phen_type)   :: phen_inst  


      ! INTERF-TODO: we will need a new bounding type (maybe?)
      ! depending on how we structure the memory
      ! We will likely have a fates_bcs (boundary conditions) type
      ! And this type could take many forms, but probably be pointers
      ! to various CLM/ALM/ATS/etc types
      ! This is a structure of pointers that maps between
      ! the calling model and the fates model
      !      type(fates_bounds_type) :: bound_fate


   contains
     
     ! Procedures for initializing FATES memory and communicators
     procedure, public :: init
     procedure, public :: site_init            
     procedure, public :: fates2dlm_init
     procedure, public :: phen_accvars_init
     
     ! Run-time procedures
     procedure, public :: dynamics_driv         ! the daily FATES timestep driver
     procedure, public :: set_fates2dlm_inst    ! wrapper for fates2clm_inst%SetValues
     procedure, public :: canopy_sunshade_fracs ! wrapper to calculate sun/shade fracs
     
  end type dlm_fates_interface_type


  ! FOR CLM: type dlm_fates_interface_type describes "clm_fates"
  !          which is defined in clm_instMod
  ! FOR ALM: type dlm_fates_interface_type describes ?

contains
    
  subroutine init(this,bounds_clump)
    
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
    type(bounds_type),intent(in)                :: bounds_clump 
    
    
    ! Initialize the mapping elements between FATES and the DLM
    
    ! These bounds are for a single clump (thread)
    ! CHANGE NAMING CONVENTION FROM ED->FATES (RGK 04-25-2016)
    allocate (this%site_inst(bounds_clump%begg:bounds_clump%endg))
    
    ! Initialize the FATES communicators with the DLM
    ! This involves to stages
    ! 1) allocate the vectors dlm_inst
    ! 2) add the history variables defined in clm_inst to the history machinery
    this%fates2dlm_inst%Init(bounds_clump)
    
    ! Initialize ED phenology variables
    ! This also involves two stages
    ! 1) allocate the vectors in phen_inst
    ! 2) add the phenology history variables to the history machinery
    this%phen_inst%Init(bounds_clump)
    
    
  end subroutine init
  
  ! -----------------------------------------------------------------------------------
  
  subroutine site_init(self,bounds_clump)
    
    
    ! CLM:  called from initialize2()
    ! ALM:  ??
    
    ! Input Arguments
    class(lsm_ed_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    
    ! locals
    integer :: g
    
    ! It is assumed that bounds is bounds_clump
    
    ! Initialize  (INTERF-TODO THIS ROUTINE CALLS CLM STUFF-MIGRATE CODE TO HERE)
    call ed_init_sites( bounds_clump,                                               &
         self%site_inst(bounds_clump%begg:bounds_clump%endg) )
    
    do g = bounds_clump%begg,bounds_clump%endg
       if (this%site_inst(g)%istheresoil) then
          call ed_update_site(this%site_inst(g))
       end if
    end do
    
    return
  end subroutine site_init
  
  ! -----------------------------------------------------------------------------------
  
  subroutine fates2dlm_init(self,bounds_clump,waterstate_inst, canopystate_inst)
    
    ! CLM:  called from initialize2()
    ! ALM:  ??
    
    ! Input Arguments
    class(lsm_ed_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    type(waterstate_type)   , intent(inout)     :: waterstate_inst
    type(canopystate_type)  , intent(inout)     :: canopystate_inst
    
    call self%fates2dlm_inst%ed_clm_link( bounds_clump,    &
         self%site_inst(bounds_clump%begg:bounds_clump%endg),    &
         self%phen_inst,                             &
         waterstate_inst,                            &
         canopystate_inst)
    
      
    return
  end subroutine fates2dlm_init
  
  ! -------------------------------------------------------------------------------------
  
  subroutine dynamics_driv(this, bounds_clump,                                          &
       atm2lnd_inst, soilstate_inst, temperature_inst,                                  &
       waterstate_inst, canopystate_inst)

    ! This wrapper is called daily from clm_driver
    ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
    ! ed_driver is not a dlm_fates_inst_type procedure because we need an extra step 
    ! to process array bounding information 
    
    implicit none
    class(lsm_ed_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    type(atm2lnd_type)      , intent(in)            :: atm2lnd_inst
    type(soilstate_type)    , intent(in)            :: soilstate_inst
    type(temperature_type)  , intent(in)            :: temperature_inst
    type(waterstate_type)   , intent(inout)         :: waterstate_inst
    type(canopystate_type)  , intent(inout)         :: canopystate_inst

    ! INTERF-TODO: REMOVE ED_DRIVER ARGUMENTS OF CLM STUCTURED TYPES AND
    ! REPLACE THEM WITH FATES_BC TYPES WITH A BOUNDS MAPPING SCHEME
    
    ! RENAME TO fates_driver()
    call ed_driver( bounds_clump,                                                       &
         this%ed_allsites_inst(bounds_clump%begg:bounds_clump%endg),                    &
         this%dlm_inst,                                                                 &
         this%phen_inst,                                                                &
         atm2lnd_inst, soilstate_inst, temperature_inst,                                &
         waterstate_inst, canopystate_inst)
    
    return
  end subroutine dynamics_driv
  
  ! -------------------------------------------------------------------------------------
  
  subroutine set_fates2dlm_inst(this,bounds_clump, setval_scalar) 

    ! This is a simple wrapper to flush some FATES -> CLM communicators

    implicit none
    class(lsm_ed_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    real(r8),intent(in) :: setval_scalar      ! This value will flush all 
    
    call this%fates2clm_inst%SetValues( bounds_clump, setval_scalar )


  end subroutine set_fates2dlm_inst

  ! -------------------------------------------------------------------------------------

  subroutine phen_accvars_init(this,bounds_clump)

    implicit none
    class(lsm_ed_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                :: bounds_clump

    call this%phen_inst%initAccVars(bounds_clump)

    return
  end subroutine phen_accvars_init
  
  ! ---------------------------------------------------------------------------------------
  
  subroutine canopy_sunshade_fracs(this,filter_nourbanp, num_nourbanp, &
                                              atm2lnd_inst,canopystate_inst)


      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.


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
    end subroutine canopy_sunshade_fracs
    
    
  end module clmfates_interfaceMod
