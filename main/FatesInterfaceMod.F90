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
   ! ------------------------------------------------------------------------------------

   use EDtypesMod            , only : ed_site_type,      &
                                      numPatchesPerCol
   use shr_kind_mod          , only : r8 => shr_kind_r8  ! INTERF-TODO: REMOVE THIS
   
   
   ! ------------------------------------------------------------------------------------
   ! Certain dimension information used by FATES is dictated by the the driver
   ! or the host model, or perhaps may be some compromise between what FATES will want
   ! in a best-case scenario and what space the driver/host will allow based on its
   ! memory constraints (most-likely due to IO)
   ! ------------------------------------------------------------------------------------
   
   type, private :: fates_dims_type
      
      integer :: numSWBands   ! Maximum number of broad-bands in the short-wave radiation
                              ! specturm to track 
                              ! (typically 2 as a default, VIS/NIR, in ED variants <2016)
      
      
   end type fates_dims_type



   ! ------------------------------------------------------------------------------------
   ! Notes on types
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! ------------------------------------------------------------------------------------
   
   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches

      ! Downwelling direct beam radiation (patch,broad-band) [W/m2?]
      real(r8),allocatable :: solad_pa(:,:)  

      ! Downwelling diffuse (I-ndirect) radiation (patch,broad-band) [W/m2]
      real(r8),allocatable :: solai_pa(:,:)


   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

   end type bc_out_type


   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), allocatable :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)

   contains
      
      procedure, public :: allocate_bcs
      procedure, public :: zero_bcs

   end type fates_interface_type

   ! ------------------------------------------------------------------------------------
   ! Dimension information is independent of which clump it is on
   ! and since these are typically read only variables and never updated on the clump
   ! we need not attach these to the threaded instances of the fates_interface_type
   ! ------------------------------------------------------------------------------------
   
   type(fates_dims_type) :: fates_dims



contains

   ! ====================================================================================

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
      ! Incrementally walk through linked list and deallocate
      
      
      
      ! Deallocate the site list
      deallocate (this%sites)
      
      return
   end subroutine fates_clean


   ! ====================================================================================
   

   subroutine allocate_bcs(this,s)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s
      
      ! Allocate input boundaries
      
      allocate(this%bc_in(s)%solad_pa(numPatchesPerCol,fates_dims%numSWBands))
      allocate(this%bc_in(s)%solai_pa(numPatchesPerCol,fates_dims%numSWBands))
      
      ! Allocate output boundaries
      
      allocate(this%bc_out(s)%fsun_pa(numPatchesPerCol))

      return
   end subroutine allocate_bcs

   ! ====================================================================================

   subroutine zero_bcs(this,s)

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s

      ! Input boundaries

      this%bc_in(s)%solad_pa(:,:) = 0.0_r8
      this%bc_in(s)%solai_pa(:,:) = 0.0_r8
      
      ! Output boundaries
      
      this%bc_out(s)%fsun_pa(:) = 0.0_r8


      return
   end subroutine zero_bcs


end module FatesInterfaceMod
