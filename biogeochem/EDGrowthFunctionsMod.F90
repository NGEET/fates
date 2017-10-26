module EDGrowthFunctionsMod

  ! ============================================================================
  ! Functions that control the trajectory of plant growth. 
  ! Ideally these would all use parameters that are fed in from the parameter file. 
  ! At present, there is only a single allocation trajectory. 
  ! ============================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals     , only : fates_log
  use EDPftvarcon        , only : EDPftvarcon_inst
  use EDTypesMod       , only : ed_cohort_type, nlevleaf, dinc_ed
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesAllometryMod, only : bleaf

  implicit none
  private

  public ::  tree_lai
  public ::  tree_sai
  public ::  mortality_rates

  logical :: DEBUG_growth = .false.

  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  real(r8) function tree_lai( cohort_in )

    ! ============================================================================
    !  LAI of individual trees is a function of the total leaf area and the total canopy area.   
    ! ============================================================================

    type(ed_cohort_type), intent(inout) :: cohort_in       

    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC

    if( cohort_in%bl  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(fates_log(),*) 'problem in treelai',cohort_in%bl,cohort_in%pft
    endif

    if( cohort_in%status_coh  ==  2 ) then ! are the leaves on? 
       slat = 1000.0_r8 * EDPftvarcon_inst%slatop(cohort_in%pft) ! m2/g to m2/kg
       leafc_per_unitarea = cohort_in%bl/(cohort_in%c_area/cohort_in%n) !KgC/m2
       if(leafc_per_unitarea > 0.0_r8)then
          tree_lai = leafc_per_unitarea * slat  !kg/m2 * m2/kg = unitless LAI 
       else
          tree_lai = 0.0_r8
       endif
    else
       tree_lai = 0.0_r8
    endif !status
    cohort_in%treelai = tree_lai

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treelai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much lai' , cohort_in%treelai , cohort_in%pft , nlevleaf * dinc_ed
    endif

    return

  end function tree_lai
  
  ! ============================================================================

  real(r8) function tree_sai( cohort_in )

    ! ============================================================================
    !  SAI of individual trees is a function of the total dead biomass per unit canopy area.   
    ! ============================================================================

    type(ed_cohort_type), intent(inout) :: cohort_in       

    real(r8) :: bdead_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: sai_scaler     

    sai_scaler = EDPftvarcon_inst%allom_sai_scaler(cohort_in%pft) 

    if( cohort_in%bdead  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(fates_log(),*) 'problem in treesai',cohort_in%bdead,cohort_in%pft
    endif

    bdead_per_unitarea = cohort_in%bdead/(cohort_in%c_area/cohort_in%n) !KgC/m2
    tree_sai = bdead_per_unitarea * sai_scaler !kg/m2 * m2/kg = unitless LAI 
   
    cohort_in%treesai = tree_sai

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treesai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much sai' , cohort_in%treesai , cohort_in%pft , nlevleaf * dinc_ed
    endif

    return

  end function tree_sai
  
  ! ============================================================================

  subroutine mortality_rates( cohort_in,cmort,hmort,bmort )

    ! ============================================================================
    !  Calculate mortality rates as a function of carbon storage       
    ! ============================================================================

    use EDParamsMod,  only : ED_val_stress_mort
    use FatesInterfaceMod,  only : hlm_use_ed_prescribed_phys

    type (ed_cohort_type), intent(in) :: cohort_in
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: b_leaf ! leaf biomass kgC
    real(r8) :: hf_sm_threshold    ! hydraulic failure soil moisture threshold 


    if (hlm_use_ed_prescribed_phys .eq. ifalse) then

    ! 'Background' mortality (can vary as a function of density as in ED1.0 and ED2.0, but doesn't here for tractability) 
    bmort = EDPftvarcon_inst%bmort(cohort_in%pft) 

    ! Proxy for hydraulic failure induced mortality. 
    hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)

    if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= hf_sm_threshold)then 
       hmort = ED_val_stress_mort
     else
       hmort = 0.0_r8
     endif 
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then
       call bleaf(cohort_in%dbh,cohort_in%hite,cohort_in%pft,cohort_in%canopy_trim,b_leaf)
       if( b_leaf > 0._r8 .and. cohort_in%bstore <= b_leaf )then
          frac = cohort_in%bstore/ b_leaf
          cmort = max(0.0_r8,ED_val_stress_mort*(1.0_r8 - frac))
        else
          cmort = 0.0_r8
       endif

    else
       write(fates_log(),*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer
    endif

    !mortality_rates = bmort + hmort + cmort

    else ! i.e. hlm_use_ed_prescribed_phys is true
       if ( cohort_in%canopy_layer .eq. 1) then
          bmort = EDPftvarcon_inst%prescribed_mortality_canopy(cohort_in%pft)
       else
          bmort = EDPftvarcon_inst%prescribed_mortality_understory(cohort_in%pft)
       endif
       cmort = 0._r8
       hmort = 0._r8
    endif

 end subroutine mortality_rates

! ============================================================================

end module EDGrowthFunctionsMod
