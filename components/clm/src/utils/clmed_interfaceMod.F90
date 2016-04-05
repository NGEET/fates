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
   
   ! Used ED Modules
   use EDtypesMod            , only : ed_patch_type, ed_site_type, numpft_ed
   use EDtypesMod            , only : ed_allsites_inst
   use EDtypesMod            , only : map_clmpatch_to_edpatch
   use EDCLMLinkMod          , only : ed_clm_type
   use EDPhenologyType       , only : ed_phenology_type
   use EDInitMod             , only : ed_init
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs
   
   implicit none
   
   
contains
   
   
   subroutine CLMEDInterf_AllocateAllSites(bounds,use_ed)
   
      ! ---------------------------------------------------------------------------------
      ! This subroutine allocates the ed_allsites_inst global.
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
      type(bounds_type),intent(in)    :: bounds 
      logical, intent(in)             :: use_ed   ! does the model want ed turned on?

      allocate (ed_allsites_inst(bounds%begg:bounds%endg))

!      if (use_ed) then
!         call ed_clm_inst%Init(bounds)
!         call ed_phenology_inst%Init(bounds)
!         call EDecophysconInit( EDpftvarcon_inst, numpft)
!      end if

      
      return
   end subroutine CLMEDInterf_AllocateAllSites

   ! -------------------------------------------------------------------------------------
    

   subroutine CLMEDInterf_ed_init(bounds, ed_clm_inst, ed_phenology_inst,  &
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
   end subroutine CLMEDInterf_ed_init

   ! -------------------------------------------------------------------------------------

  
   subroutine CLMEDInterf_CanopySunShadeFracs(bounds,forc_solad, forc_solai,  &
                                              filter_nourbanp, num_nourbanp,  fsun)      
      

      ! ----------------------------------------------------------------------------------
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      ! ----------------------------------------------------------------------------------

      implicit none

      ! Input Arguments

      type(bounds_type)       , intent(in)            :: bounds  ! clump bounds

      ! direct beam radiation (W/m**2) ! => atm2lnd_inst%forc_solad_grc
      real(r8),intent(in) :: forc_solad(bounds%begg:,:)
      
      ! diffuse radiation (W/m**2)     ! => atm2lnd_inst%forc_solai_grc
      real(r8),intent(in) :: forc_solai(bounds%begg:,:)
      
      ! patch filter for non-urban points
      integer, intent(in),dimension(:)    :: filter_nourbanp

      ! number of patches in non-urban points in patch  filter
      integer, intent(in)                 :: num_nourbanp       

      ! Output Arguments to CLM
      ! sunlit fraction of canopy  ! => canopystate_inst%fsun_patch  
      real(r8),intent(inout) :: fsun(bounds%begp:)

      ! Local Variables
      integer  :: fp                          ! non-urban filter patch index
      integer  :: p                           ! patch index
      integer  :: g                           ! grid cell index
      integer, parameter :: ipar = 1          ! The band index for PAR
      type(ed_patch_type), pointer :: cpatch  ! c"urrent" patch

      ! global ed_allsites_inst(g) should be allocated begg:endg
      ! if this run is using clumps (SMP?) the patch%gricell pointers
      ! should reflect that.
      
      do fp = 1,num_nourbanp
         
         p = filter_nourbanp(fp)
         g = patch%gridcell(p)

         if ( patch%is_veg(p) ) then 
            cpatch => map_clmpatch_to_edpatch(ed_allsites_inst(g), p) 
            
            call ED_SunShadeFracs(cpatch,forc_solad(g,ipar),forc_solai(g,ipar),fsun(p))
            

         endif
         
      end do
      return
   end subroutine CLMEDInterf_CanopySunShadeFracs

   
end module clmed_interfaceMod
