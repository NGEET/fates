module EDMortalityFunctionsMod

  ! ============================================================================
  ! Functions that control mortality.
  ! ============================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals     , only : fates_log
  use EDPftvarcon      , only : EDPftvarcon_inst
  use EDTypesMod       , only : ed_cohort_type
  use FatesConstantsMod, only : itrue,ifalse
  use FatesAllometryMod, only : bleaf
  use EDParamsMod      , only : ED_val_stress_mort
  use FatesInterfaceMod, only : hlm_use_ed_prescribed_phys

  implicit none
  private


  public ::  mortality_rates

  logical :: DEBUG_growth = .false.

  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains



  subroutine mortality_rates( cohort_in,cmort,hmort,bmort )

    ! ============================================================================
    !  Calculate mortality rates as a function of carbon storage       
    ! ============================================================================

   

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

end module EDMortalityFunctionsMod
