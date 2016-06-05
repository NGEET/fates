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

   use EDtypesMod            , only : ed_site_type
   
   
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
!      procedure, public :: fates_clean

   end type fates_interface_type

contains

   ! ------------------------------------------------------------------------------------

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

end module FatesInterfaceMod
