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
   use FatesInterfaceMod, only : numpft
   use FatesConstantsMod, only : g_per_kg
   use FatesInterfaceMod, only : bc_in_type
   use FatesLitterMod,    only : litter_type
   use FatesLitterMod,    only : ncwd
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
   public :: InitMassBalanceColdStart

   character(len=*), parameter, private :: sourcefile = &
        __FILE__

contains
   
  subroutine InitMassBalanceColdStart(sites)
    
    ! arguments
    type(ed_site_type), intent(inout), target :: sites(:)

    ! locals
    type(site_massbal_type),pointer :: site_mass
    integer :: el
    integer :: s
    real(r8) :: total_stock
    real(r8) :: biomass_stock
    real(r8) :: litter_stock
    real(r8) :: seed_stock
    
    do s = 1,size(sites(:),dim=1)
    
       do el = 1,num_elements
      
          site_mass => sites(s)%mass_balance(el)
          
          call SiteMassStock(sites(s),el,total_stock,biomass_stock,litter_stock,seed_stock)
          
          site_mass%old_stock     = total_stock
          
       end do
       
    end do
    
    return
  end subroutine InitMassBalanceColdStart  

  ! ==============================================================================================

  subroutine SiteMassStock(currentSite,i_element,total_stock,biomass_stock,litter_stock,seed_stock)
     
     type(ed_site_type),intent(inout),target :: currentSite
     integer,intent(in)                      :: i_element    ! This is the element index
                                                             ! in FATES (not the parteh global id)
     real(r8),intent(out)                    :: total_stock    ! kg
     real(r8),intent(out)                    :: litter_stock   ! kg
     real(r8),intent(out)                    :: biomass_stock  ! kg
     real(r8),intent(out)                    :: seed_stock     ! kg
     integer                                 :: element_id     ! parteh element id
     type(litter_type), pointer              :: litt           ! litter object
     type(ed_patch_type), pointer            :: currentPatch
     type(ed_cohort_type), pointer           :: currentCohort
     
     litter_stock  = 0.0_r8
     biomass_stock = 0.0_r8
     seed_stock    = 0.0_r8

     element_id = element_list(i_element)

     currentPatch => currentSite%oldest_patch 
     do while(associated(currentPatch))

        litt => currentPatch%litter(i_element)

        ! Total non-seed litter in [kg]
        litter_stock = litter_stock + currentPatch%area * &
             (sum(litt%ag_cwd)                  + &
              sum(litt%bg_cwd) + &
              sum(litt%leaf_fines)              + &
              sum(litt%root_fines))

        ! Total mass of viable seeds in [kg]
        seed_stock = seed_stock + currentPatch%area * &
             (sum(litt%seed) + sum(litt%seed_germ))

        ! Total mass on living plants
        currentCohort => currentPatch%tallest
        do while(associated(currentCohort))
           biomass_stock =  biomass_stock + &
                 (currentCohort%prt%GetState(struct_organ,element_id) + &
                 currentCohort%prt%GetState(sapw_organ,element_id) + &
                 currentCohort%prt%GetState(leaf_organ,element_id) + &
                 currentCohort%prt%GetState(fnrt_organ,element_id) + &
                 currentCohort%prt%GetState(store_organ,element_id) + &
                 currentCohort%prt%GetState(repro_organ,element_id) ) &
                 * currentCohort%n
           currentCohort => currentCohort%shorter
        enddo !end cohort loop 
        currentPatch => currentPatch%younger
     enddo !end patch loop
     
     total_stock = biomass_stock + seed_stock + litter_stock

     return
  end subroutine SiteMassStock
  
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
             
             if(litt%leaf_fines(pft)<0._r8)then
                write(fates_log(),*) 'For PFT: ',pft
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative leaf fine litter: ',litt%leaf_fines(pft)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             do ilyr = 1,numlevsoil
                if(litt%root_fines(pft,ilyr)<0._r8)then
                   write(fates_log(),*) 'For PFT: ',pft
                   write(fates_log(),*) 'Soil layer: ',ilyr
                   write(fates_log(),*) 'Element id: ',element_id
                   write(fates_log(),*) 'Negative leaf fine litter: ',litt%root_fines(pft,ilyr)
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
