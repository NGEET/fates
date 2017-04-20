
module EDLoggingMortalityMod

! ====================================================================================
!  Purpose: 1. create logging mortalities: 
!           (a)logging mortality (size class level)
!           (b)collateral mortality (size class level)
!           (c)infrastructure mortality (cohort level)

!           2. move the logged trunk fluxes from live into product pool 
!              
!           3. move logging-associated mortality fluxes from live to CWD

!           4. keep carbon balance (in ed_total_balance_check)
!
!  Yi Xu
!  Date: 2017
! ====================================================================================

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use EDTypesMod       , only : ed_cohort_type, ed_patch_type,ncwd, ed_site_type, ed_resources_management_type
  use clm_varctl       , only : logging_time
  use pftconMod        , only : pftcon

  implicit none
  private
  
  public :: LoggingMortality_rates
  public :: No_LoggingMortality_rates
  public :: logging_litter_fluxes
  public :: Logging_threshold
  
contains
 
 subroutine Logging_threshold (minimum_diameter_logging, threshold_sizeclass)
    ! Calculate the threshold of size class for logging
	
	
	real(r8), intent(in)  :: minimum_diameter_logging
	integer, intent(out) :: threshold_sizeclass 
	

	if (minimum_diameter_logging >= 0.0_r8 .and. minimum_diameter_logging< 5.0_r8) then 
		threshold_sizeclass = 1
	elseif (minimum_diameter_logging>=5.0_r8 .and. minimum_diameter_logging<10.0_r8) then
		threshold_sizeclass = 2
	elseif (minimum_diameter_logging>=10.0_r8 .and. minimum_diameter_logging<15.0_r8) then
		threshold_sizeclass = 3
	elseif (minimum_diameter_logging>=15.0_r8 .and. minimum_diameter_logging<20.0_r8) then
		threshold_sizeclass = 4
	elseif (minimum_diameter_logging>=20.0_r8 .and. minimum_diameter_logging<30.0_r8) then
		threshold_sizeclass = 5
	elseif (minimum_diameter_logging>=30.0_r8 .and. minimum_diameter_logging<40.0_r8) then
		threshold_sizeclass = 6
	elseif (minimum_diameter_logging>=40.0_r8 .and. minimum_diameter_logging<50.0_r8) then
		threshold_sizeclass = 7
	elseif (minimum_diameter_logging>=50.0_r8 .and. minimum_diameter_logging<60.0_r8) then
		threshold_sizeclass = 8
	elseif (minimum_diameter_logging>=60.0_r8 .and. minimum_diameter_logging<70.0_r8) then
		threshold_sizeclass = 9
	elseif (minimum_diameter_logging>=70.0_r8 .and. minimum_diameter_logging<80.0_r8) then
		threshold_sizeclass = 10
	elseif (minimum_diameter_logging>=80.0_r8 .and. minimum_diameter_logging<90.0_r8) then
		threshold_sizeclass = 11
	elseif (minimum_diameter_logging>=90.0_r8 .and. minimum_diameter_logging<100.0_r8) then
		threshold_sizeclass = 12
	elseif (minimum_diameter_logging>100.0_r8) then
		threshold_sizeclass = 13
	end if
	
 
 end subroutine Logging_threshold
 

 
 subroutine LoggingMortality_rates( site_in, pft_i, size_class, threshold_sizeclass, lmort_logging,lmort_collateral,lmort_infra )

    
   integer, intent(in) :: size_class  !An index that indicates which diameter size bin the cohort currently resides in
   integer, intent(in) :: threshold_sizeclass ! threshold of size class for logging 
   integer, intent(in) :: pft_i
   
   real(r8),intent(out) :: lmort_logging ! logging mortality_rates, share same size class
   real(r8),intent(out) :: lmort_collateral  ! logging impact mortality_rates, share same size class
   real(r8),intent(out) :: lmort_infra  ! infrastructure mortality_rates, share same size class
  
   real(r8)  :: adjustment ! adjustment for mortality rates
   type(ed_site_type) , intent(inout), target :: site_in
   
   
   adjustment=1.0_r8
   ! Read in rules of logging mortalities
   
   
  

   site_in%resouces_management%logging_collatoral_mortality_rate = site_in%resouces_management%logging_ratio*site_in%resouces_management%fraction_trees_logged          ! collaterally damaged rate        %/per logging activity
   site_in%resouces_management%logging_infrastructure_mortality_rate = site_in%resouces_management%logging_ratio*site_in%resouces_management%fraction_trees_logged		! mechanically damaged rate        %/per logging activity
	

 
   
   if (logging_time) then 
   
      if(pftcon%woody(pft_i) == 1)then ! only set logging rates for trees
		  ! Log trees whose DBH > 50cm 
		  ! Pass logging rates to cohort level 
		  if (size_class>=threshold_sizeclass) then
			lmort_logging=site_in%resouces_management%fraction_trees_logged*adjustment
		  end if
		  lmort_collateral=site_in%resouces_management%logging_collatoral_mortality_rate*adjustment
		  lmort_infra=site_in%resouces_management%logging_infrastructure_mortality_rate*adjustment
	   else
		  lmort_logging=0.0_r8
		  lmort_collateral=0.0_r8
		  lmort_infra=0.0_r8
	   end if
	   
	   
   else 

	  lmort_logging=0.0_r8
	  lmort_collateral=0.0_r8
	  lmort_infra=0.0_r8
   end if 
   
   
   
 end subroutine LoggingMortality_rates
 
 subroutine No_LoggingMortality_rates(lmort_logging,lmort_collateral,lmort_infra )
	 real(r8),intent(out) :: lmort_logging ! logging mortality_rates, share same size class
     real(r8),intent(out) :: lmort_collateral  ! logging impact mortality_rates, share same size class
	 real(r8),intent(out) :: lmort_infra  ! infrastructure mortality_rates, share same size class
   
 	  lmort_logging=0.0_r8
	  lmort_collateral=0.0_r8
	  lmort_infra=0.0_r8
	  
 end subroutine
 
! ============================================================================
  subroutine logging_litter_fluxes(cp_Site, cp_target, new_patch_target, patch_site_areadis)
	
	
	!  DESCRIPTION:
    !  Carbon going from ongoing mortality into CWD pools. 
	!  This module includes all the mortalities except fire mortality and logging mortalities
    !  Purpose: 
    !	1) move logging-associated mortality fluxes from live to CWD
	!   2) move the logging fluxes from live into product pool 
	!   3) carbon balance check 
	!  E.g,:
	!  Remove trunk of logged trees from litter/CWD
	!  Add other parts of logged trees and all parts of collaterally and mechanically damaged trees into litter  
	

	!USES:
	use EDParamsMod,  only : ED_val_ag_biomass, ED_val_understorey_death
    use SFParamsMod,  only : SF_val_cwd_frac
	use EDtypesMod, only : area
	use EDtypesMod, only : ed_site_type, ed_patch_type,ed_cohort_type
	use pftconMod, only : pftcon
	
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(inout), target :: cp_Site 
    type(ed_patch_type) , intent(inout), target :: cp_target 
    type(ed_patch_type) , intent(inout), target :: new_patch_target

    real(r8)            , intent(in)            :: patch_site_areadis
	
	!LOCAL VARIABLES:
    real(r8) :: cwd_litter_density
    real(r8) :: litter_area ! area over which to distribute this litter. 
	
	
	type(ed_site_type)  , pointer :: currentSite
	type(ed_patch_type) , pointer :: currentPatch 
	type(ed_patch_type) , pointer :: new_patch
	type(ed_cohort_type), pointer :: currentCohort

   
    real(r8) :: bcroot               ! amount of below ground coarse root per cohort  kgC. (goes into CWD_BG)
    real(r8) :: bstem                ! amount of above ground stem biomass per cohort  kgC.(goes into CWG_AG)
	real(r8) :: dead_tree_density    ! # trees killed by all logging per m2
	real(r8) :: damaged_dead_tree_density    ! # trees killed by collateral and mechanical damages per m2
	
	real(r8) :: logging_density      ! # logged trees per m2 (logged trunk)
	
	real(r8) :: np_mult           !Fraction of the new patch which came from the current patch (and so needs the same litter) 
    integer :: p,c
	 
    currentSite  => cp_Site
    currentPatch => cp_target
    new_patch => new_patch_target
	

	 
	if (currentPatch%logging==1) then ! 
		
		currentCohort => currentPatch%shortest
		do while(associated(currentCohort))       
		   p = currentCohort%pft
		   if(currentPatch%disturbance_rates(3) > currentPatch%disturbance_rates(2) .and. currentPatch%disturbance_rates(3) > currentPatch%disturbance_rates(1) )then !mortality is dominant disturbance 
			  if(pftcon%woody(p) == 1)then   !woody=1 trees only, woody=0 grass, canopy_layer=0 woody=1 small trees
			     !******************************************/
				! PART 1) Put twig, small branch and large branch of dead trees biomass from logging and damage into above ground litter pool
				!         Put below ground dead trees biomass from logging and damage into below ground litter pool
				! This happens before the plant number of selective logging have been updated (need to check other codes), 
				! but the effects of other mortalities due on plant number have been updated in ED_integrate_state_variable
			  
				 !************************************/ 
				 ! Number of trees that died because of the selective logging, per m2 of ground. 

				 ! total amount of above-ground tree biomass in patch. kgC/m2 in SFMainMod.F90 
				 !(currentCohort%bl+ED_val_ag_biomass* (currentCohort%bsw + currentCohort%bdead))*currentCohort%n
				 ! So remove leaf biomass, stem biomass per tree is (currentCohort%bsw + currentCohort%bdead) * ED_val_ag_biomass  
					 
				 ! stem biomass per tree
				  bstem  = (currentCohort%bsw + currentCohort%bdead) * ED_val_ag_biomass     
				 
				 ! coarse root biomass per tree
				  bcroot = (currentCohort%bsw + currentCohort%bdead) * (1.0_r8 - ED_val_ag_biomass)
				  
				 ! density of dead trees per m2. Beware /AREA= 10000.0_r8 ! Notional area of simulated forest m2
				 dead_tree_density  = (currentCohort%lmort_logging + currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n / AREA  
				 
				 ! direct logged dead tree per m2 Beware /AREA= 10000.0_r8 ! Notional area of simulated forest m2
				 damaged_dead_tree_density = (currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n / AREA  
				 
				 
				  !new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + dead_tree_density * (currentCohort%bl) * (1.0_r8-currentCohort%cfa)
				  new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + dead_tree_density * (currentCohort%bl) 
				  new_patch%root_litter(p) = new_patch%root_litter(p) + dead_tree_density * (currentCohort%br+currentCohort%bstore)
			 
				  !currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + dead_tree_density * (currentCohort%bl) * (1.0_r8-currentCohort%cfa)
				  currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + dead_tree_density * (currentCohort%bl)
				  currentPatch%root_litter(p) = currentPatch%root_litter(p) + dead_tree_density * (currentCohort%br+currentCohort%bstore)
      
				  ! track as diagnostic fluxes
				  currentSite%leaf_litter_diagnostic_input_carbonflux(p) = currentSite%leaf_litter_diagnostic_input_carbonflux(p) + &
                  (currentCohort%bl) *  (currentCohort%lmort_logging + currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n/ AREA
				  
                  currentSite%root_litter_diagnostic_input_carbonflux(p) = currentSite%root_litter_diagnostic_input_carbonflux(p) + &
                  (currentCohort%br+currentCohort%bstore) * (currentCohort%lmort_logging + currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n/ AREA

				  
	              !above ground coarse woody debris of twigs, small branches and large branches are from logged and damaged trees 
				   do c = 1,3
						new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem
						currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem 
						
						! track as diagnostic fluxes
						currentSite%CWD_AG_diagnostic_input_carbonflux(c) = currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
						 SF_val_CWD_frac(c) * bstem * (currentCohort%lmort_logging + currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n / AREA
						
						
				   enddo
				   

				   ! above ground coarse woody debris of trunk are only from damaged trees.The logged trees of this part are transported off-site. 
			       do c = 4,4
						new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + damaged_dead_tree_density * SF_val_CWD_frac(c) * bstem
						currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + damaged_dead_tree_density * SF_val_CWD_frac(c) * bstem
						
						!*******************************************************************/
						! PART 2 Save logged trunk into product tool
						! PART 3 Pass product tool into flux out for carbon balance check 
						!*******************************************************************/
						
						
						 logging_density = (currentCohort%lmort_logging * currentCohort%n) / AREA  
						
						!currentPatch%trunk_product unit is kGC/m2
						 currentPatch%trunk_product = currentPatch%trunk_product + logging_density * SF_val_CWD_frac(c) * bstem

						 currentSite%flux_out = currentSite%flux_out + logging_density * SF_val_CWD_frac(c) * bstem * AREA
						 
						 
						! track as diagnostic fluxes
						 currentSite%CWD_AG_diagnostic_input_carbonflux(c) = currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
						 SF_val_CWD_frac(c) * bstem * (currentCohort%lmort_logging * currentCohort%n) / AREA
						
						
						
				   enddo
			      
				  !below ground coarse woody debris of all sizes are from logged and damaged trees
				   do c = 1,ncwd
						new_patch%cwd_bg(c) = new_patch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot
						currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot

						! track as diagnostic fluxes
						currentSite%CWD_BG_diagnostic_input_carbonflux(c) = currentSite%CWD_BG_diagnostic_input_carbonflux(c) + &
						 SF_val_CWD_frac(c) * bcroot * (currentCohort%lmort_logging + currentCohort%lmort_collateral + currentCohort%lmort_infra )* currentCohort%n/ AREA  
						
						
				   enddo
				   

			   end if
			   
			end if
		
			   currentCohort => currentCohort%taller      

			   
		enddo !currentCohort
		
	
	end if
	
	
  end subroutine logging_litter_fluxes
 
end module EDLoggingMortalityMod