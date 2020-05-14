module ChecksBalancesMod

   use shr_kind_mod,      only : r8 => shr_kind_r8
   use shr_const_mod,     only : SHR_CONST_CDAY
   use EDtypesMod,        only : ed_site_type
   use EDTypesMod,        only : ed_patch_type
   use EDTypesMod,        only : ed_cohort_type
   use EDTypesMod,        only : AREA
   use EDTypesMod,        only : site_massbal_type
   use EDTypesMod,        only : num_elements
   use EDTypesMod,        only : element_list
   use FatesInterfaceTypesMod, only : numpft
   use FatesConstantsMod, only : g_per_kg
   use FatesInterfaceTypesMod, only : bc_in_type
   use FatesLitterMod,    only : litter_type
   use FatesLitterMod,    only : ncwd
   use FatesLitterMod,    only : ndcmpy
   use PRTGenericMod,     only : all_carbon_elements
   use PRTGenericMod,     only : leaf_organ
   use PRTGenericMod,     only : fnrt_organ
   use PRTGenericMod,     only : sapw_organ
   use PRTGenericMod,     only : store_organ
   use PRTGenericMod,     only : repro_organ
   use PRTGenericMod,     only : struct_organ
   use FatesGlobals,      only : fates_log 
   use shr_log_mod,       only : errMsg => shr_log_errMsg
   use FatesGlobals,      only : endrun => fates_endrun

   implicit none
   
   private
   public :: SiteMassStock
   public :: PatchMassStock

   character(len=*), parameter, private :: sourcefile = &
        __FILE__

contains
   
  ! ==============================================================================================
  
  subroutine SiteMassStock(currentSite,el,total_stock,biomass_stock,litter_stock,seed_stock)
  
     type(ed_site_type),intent(inout),target :: currentSite
     integer,intent(in)                      :: el           ! This is the element index
                                                             ! in FATES (not the parteh global id)
     real(r8),intent(out)                    :: total_stock    ! kg
     real(r8),intent(out)                    :: litter_stock   ! kg
     real(r8),intent(out)                    :: biomass_stock  ! kg
     real(r8),intent(out)                    :: seed_stock     ! kg
     type(ed_patch_type), pointer            :: currentPatch
     type(ed_cohort_type), pointer           :: currentCohort
     real(r8)                                :: patch_biomass  ! kg
     real(r8)                                :: patch_seed     ! kg
     real(r8)                                :: patch_litter   ! kg
       
     litter_stock  = 0.0_r8
     biomass_stock = 0.0_r8
     seed_stock    = 0.0_r8

     currentPatch => currentSite%oldest_patch 
     do while(associated(currentPatch))

        call PatchMassStock(currentPatch,el,patch_biomass,patch_seed,patch_litter)
        litter_stock  = litter_stock + patch_litter
        biomass_stock = biomass_stock + patch_biomass
        seed_stock    = seed_stock + patch_seed

        currentPatch => currentPatch%younger
     enddo !end patch loop
     
     total_stock = biomass_stock + seed_stock + litter_stock

     return
  end subroutine SiteMassStock

  ! =====================================================================================

  subroutine PatchMassStock(currentPatch,el,live_stock,seed_stock,litter_stock)

      ! ---------------------------------------------------------------------------------
      ! Sum up the mass of the different stocks on a patch for each element
      ! ---------------------------------------------------------------------------------
      type(ed_patch_type),intent(inout),target :: currentPatch
      integer,intent(in)                       :: el
      real(r8),intent(out)                     :: live_stock
      real(r8),intent(out)                     :: seed_stock
      real(r8),intent(out)                     :: litter_stock

      type(litter_type), pointer            :: litt           ! litter object
      type(ed_cohort_type), pointer         :: currentCohort
      integer                               :: element_id

      litt => currentPatch%litter(el)
      element_id = element_list(el)

      ! Total non-seed litter in [kg]
      litter_stock = currentPatch%area * &
            (sum(litt%ag_cwd)                  + &
            sum(litt%bg_cwd) + &
            sum(litt%leaf_fines)              + &
            sum(litt%root_fines))
      
        ! Total mass of viable seeds in [kg]
      seed_stock = currentPatch%area * &
            (sum(litt%seed) + sum(litt%seed_germ))

      ! Total mass on living plants
      live_stock = 0._r8
      currentCohort => currentPatch%tallest
      do while(associated(currentCohort))
          live_stock = live_stock + &
                (currentCohort%prt%GetState(struct_organ,element_id) + &
                currentCohort%prt%GetState(sapw_organ,element_id) + &
                currentCohort%prt%GetState(leaf_organ,element_id) + &
                currentCohort%prt%GetState(fnrt_organ,element_id) + &
                currentCohort%prt%GetState(store_organ,element_id) + &
                currentCohort%prt%GetState(repro_organ,element_id) ) &
                * currentCohort%n
          currentCohort => currentCohort%shorter
      enddo !end cohort loop 

      return
  end subroutine PatchMassStock


  
  ! =====================================================================================
  
  subroutine CheckLitterPools(currentSite,bc_in)

    ! -----------------------------------------------------------------------------------
    !
    ! This subroutine checks that the litter pools do not have weird values.
    ! Currently, we are only checking for negatives.  
    !
    ! This is not a carbon balance check.
    ! 
    ! We will usually keep this routine turned off, unless we are actively 
    ! debugging.
    !
    ! -----------------------------------------------------------------------------------

    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(bc_in_type), intent(in)               :: bc_in
    
    ! Local variables
    type(ed_patch_type), pointer :: currentPatch
    type(litter_type), pointer   :: litt    ! Litter object
    integer :: el                           ! Litter element loop index
    integer :: element_id                   ! parteh consistent litter index
    integer :: c                            ! CWD loop index
    integer :: ilyr                         ! soil layer index
    integer :: pft                          ! pft index
    integer :: dcmpy                        ! decomposability index
    integer :: numlevsoil                   ! number of soil layers

    ! We only really run this scheme if we think things are really broken.  
    ! The balance checks should be our first line of defense that are
    ! always on.

    logical, parameter :: do_litter_checks = .true.

    
    if(.not.do_litter_checks) return  

    numlevsoil = bc_in%nlevsoil
    
    do el = 1, num_elements
       
       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))
          
          litt => currentPatch%litter(el)
          element_id = litt%element_id

          do c = 1,ncwd
             if(litt%ag_cwd(c)<0._r8) then
                write(fates_log(),*) 'In pool: ',c
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative AG CWD: ',litt%ag_cwd(c)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
                
             do ilyr = 1,numlevsoil
                if(litt%bg_cwd(c,ilyr)<0._r8) then
                   write(fates_log(),*) 'In pool: ',c
                   write(fates_log(),*) 'Soil layer: ',ilyr
                   write(fates_log(),*) 'Element id: ',element_id
                   write(fates_log(),*) 'Negative BG CWD: ',litt%bg_cwd(c,ilyr)
                   write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
                
             end do
          end do
          
          do pft = 1,numpft
             
             if(litt%seed(pft)<0._r8) then
                write(fates_log(),*) 'For PFT: ',pft
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative seed pool: ',litt%seed(pft)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             

          end do

          do dcmpy = 1,ndcmpy
             
             if(litt%leaf_fines(dcmpy)<0._r8)then
                write(fates_log(),*) 'For PFT: ',pft
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative leaf fine litter: ',litt%leaf_fines(dcmpy)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             do ilyr = 1,numlevsoil
                if(litt%root_fines(dcmpy,ilyr)<0._r8)then
                   write(fates_log(),*) 'For PFT: ',dcmpy
                   write(fates_log(),*) 'Soil layer: ',ilyr
                   write(fates_log(),*) 'Element id: ',element_id
                   write(fates_log(),*) 'Negative leaf fine litter: ',litt%root_fines(dcmpy,ilyr)
                   write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
             end do
             
          end do
          currentPatch => currentPatch%older
       end do
    end do
    
    return
  end subroutine CheckLitterPools




   
end module ChecksBalancesMod
