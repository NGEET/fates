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
   use decompMod         , only : bounds_type
   use CanopyStateType   , only : canopystate_type
   use atm2lndType       , only : atm2lnd_type
   use ncdio_pio         , only : file_desc_t
   use PatchType         , only : patch
   use ColumnType        , only : col
   ! ------------------------------------------------------------------------------------

   use EDtypesMod            , only : ed_patch_type, ed_site_type, numpft_ed
   use EDtypesMod            , only : map_clmpatch_to_edpatch
   use EDSurfaceRadiationMod , only : ED_SunShadeFracs
   use EDInitMod             , only : ed_init_sites
   use EDMainMod             , only : ed_update_site
   use EDRestVectorMod       , only : EDRest
   
   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated
      ! with a grid-cell, this is intended to be migrated to columns
      ! prev:  type(ed_site_type)::ed_allsites_inst

      integer                         :: nsites

      type(ed_site_type), allocatable :: sites(:)
      
      ! INTERF-TODO ADD THE DLM->FATES BOUNDARY CONDITION CLASS
      ! These are boundary condition variables populated by the DLM
      ! type(fates_bc_type) :: fatesbc
      
   contains
      
      ! Procedures for initializing FATES threaded memory and communicators
      procedure, public :: init
      procedure, public :: fates_clean
      procedure, public :: site_init
      procedure, public :: fates_restart
      procedure, public :: canopy_sunshade_fracs

   end type fates_interface_type

contains

!   subroutine init(this,bounds_clump)
!
!      implicit none
!      
!      ! Input Arguments
!      class(fates_interface_type), intent(inout) :: this
!
!      ! INTERF-TODO:  AS THE FATES PUBLIC API- BOUNDS CLUMP WILL NOT BE ALLOWED
!      ! IN HERE FOR MUCH LONGER.
!      type(bounds_type),intent(in)            :: bounds_clump 
!      
!      ! Initialize the mapping elements between FATES and the DLM
!      
!      ! These bounds are for a single clump (thread)
!      allocate (this%sites(this%nsites))
!      
!      return
!   end subroutine init
   
   ! ------------------------------------------------------------------------------------

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this,bounds_clump)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)               :: bounds_clump 
      
      ! Incrementally walk through linked list and deallocate
      
      ! Deallocate the site list
      deallocate (this%sites)
      
      return
   end subroutine fates_clean

   ! ------------------------------------------------------------------------------------

   subroutine init_coldstart(this,fcolumn)
         
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in)                        :: fcolumn(this%nsites)
      
      ! locals
      integer :: s
      integer :: c
      integer :: g
      
      do s = 1,this%nsites

         call zero_site(this%sites(s))
         
         c = fcolumn(s)
         g = col%gridcell(c)  ! TODO-INTERF: col% and grc% should not be accessible here

         this%sites(s)%lat = grc%latdeg(g)  
         this%sites(s)%lon = grc%londeg(g)

      end do

      call set_site_properties(this%sites,this%nsites)

      call init_patches(this%sites, this%nsites)

      do s = 1,this%nsites
         call ed_update_site(this%sites(s))
      end do
      
      return
   end subroutine init_coldstart
   
   ! ------------------------------------------------------------------------------------
   
   subroutine init_restart(this, bounds_clump, ncid, flag )
      
      implicit none
      class(fates_interface_type), intent(inout)  :: this
      type(bounds_type),intent(in)                :: bounds_clump
      type(file_desc_t)       , intent(inout)     :: ncid    ! netcdf id
      character(len=*)        , intent(in)        :: flag    !'read' or 'write'
      
      call EDRest( bounds_clump, this%sites, this%nsites,
            ncid, flag )
      return
   end subroutine init_restart

   ! ------------------------------------------------------------------------------------
   
   subroutine canopy_sunshade_fracs(this ,filter_nourbanp, num_nourbanp, &
         atm2lnd_inst,canopystate_inst)
         
      
      ! TODO-INTERF: THIS ROUTINE NEEDS TO BE WRAPPED BY A CLM_FATES CALL
      !              IN THAT CALL THE BOUNDARY CONDITIONS SHOULD BE PREPPED
      !              SO THAT THIS CALL DOES NOT HAVE CLM TYPES HERE

      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
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
              cpatch => map_clmpatch_to_edpatch(this%sites(g), p) 
              
              call ED_SunShadeFracs(cpatch,forc_solad(g,ipar),forc_solai(g,ipar),fsun(p))
              
           endif
           
        end do
      end associate
      return
   end subroutine canopy_sunshade_fracs

   


end module FatesInterfaceMod
