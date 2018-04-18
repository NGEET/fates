module EDMortalityFunctionsMod

  ! ============================================================================
  ! Functions that control mortality.
  ! ============================================================================

   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesGlobals          , only : fates_log
   use EDPftvarcon           , only : EDPftvarcon_inst
   use EDTypesMod            , only : ed_cohort_type
   use EDTypesMod            , only : ed_site_type
   use EDTypesMod            , only : ed_patch_type
   use FatesConstantsMod     , only : itrue,ifalse
   use FatesAllometryMod     , only : bleaf
   use FatesAllometryMod     , only : storage_fraction_of_target
   use FatesInterfaceMod     , only : bc_in_type
   use FatesInterfaceMod     , only : hlm_use_ed_prescribed_phys
   use FatesInterfaceMod     , only : hlm_freq_day
   use EDLoggingMortalityMod , only : LoggingMortality_frac
   use EDParamsMod           , only : fates_mortality_disturbance_fraction
   use FatesInterfaceMod     , only : bc_in_type


   implicit none
   private
   
   
   public :: mortality_rates
   public :: Mortality_Derivative
   
   logical :: DEBUG_growth = .false.
   
   ! ============================================================================
   ! 10/30/09: Created by Rosie Fisher
   ! 02/20/18: Refactored Ryan Knox
   ! ============================================================================


contains



  subroutine mortality_rates( cohort_in,bc_in,cmort,hmort,bmort,frmort )

    ! ============================================================================
    !  Calculate mortality rates from carbon storage, hydraulic cavitation, 
    !  background and freezing
    ! ============================================================================
    
    use FatesConstantsMod,  only : tfrz => t_water_freeze_k_1atm
   

    type (ed_cohort_type), intent(in) :: cohort_in 
    type (bc_in_type), intent(in) :: bc_in
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality
    real(r8),intent(out) :: frmort ! freezing stress mortality

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: b_leaf ! target leaf biomass kgC
    real(r8) :: hf_sm_threshold    ! hydraulic failure soil moisture threshold 
    real(r8) :: temp_dep           ! Temp. function (freezing mortality)
    real(r8) :: temp_in_C          ! Daily averaged temperature in Celcius
    real(r8),parameter :: frost_mort_buffer = 5.0_r8  ! 5deg buffer for freezing mortality

    logical, parameter :: test_zero_mortality = .false. ! Developer test which
                                                        ! may help to debug carbon imbalances
                                                        ! and the like

    if (hlm_use_ed_prescribed_phys .eq. ifalse) then

    ! 'Background' mortality (can vary as a function of 
    !  density as in ED1.0 and ED2.0, but doesn't here for tractability) 
    bmort = EDPftvarcon_inst%bmort(cohort_in%pft) 

    ! Proxy for hydraulic failure induced mortality. 
    hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)

    if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= hf_sm_threshold)then 
       hmort = EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
     else
       hmort = 0.0_r8
     endif 
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then
       call bleaf(cohort_in%dbh,cohort_in%pft,cohort_in%canopy_trim,b_leaf)
       call storage_fraction_of_target(b_leaf, cohort_in%bstore, frac)
       if( frac .lt. 1._r8) then
          cmort = max(0.0_r8,EDPftvarcon_inst%mort_scalar_cstarvation(cohort_in%pft) * &
               (1.0_r8 - frac))
       else
          cmort = 0.0_r8
       endif

    else
       write(fates_log(),*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer
    endif


    !    Mortality due to cold and freezing stress (frmort), based on ED2 and:           
    !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions 
    !           of land-use change, CO2 fertilization, and climate variability to the    
    !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              
    !           doi: 10.1111/j.1365-2486.2006.01254.x                                    

    temp_in_C = bc_in%t_veg24_si - tfrz
    temp_dep  = max(0.0,min(1.0,1.0 - (temp_in_C - &
                EDPftvarcon_inst%freezetol(cohort_in%pft))/frost_mort_buffer) )
    frmort    = EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft) * temp_dep


    !mortality_rates = bmort + hmort + cmort

    else ! i.e. hlm_use_ed_prescribed_phys is true
       if ( cohort_in%canopy_layer .eq. 1) then
          bmort = EDPftvarcon_inst%prescribed_mortality_canopy(cohort_in%pft)
       else
          bmort = EDPftvarcon_inst%prescribed_mortality_understory(cohort_in%pft)
       endif
       cmort  = 0._r8
       hmort  = 0._r8
       frmort = 0._r8
    endif

    if (test_zero_mortality) then
       cmort = 0.0_r8
       hmort = 0.0_r8
       frmort = 0.0_r8
       bmort = 0.0_r8
    end if
       
    return
 end subroutine mortality_rates

 ! ============================================================================

 subroutine Mortality_Derivative( currentSite, currentCohort, bc_in)

    !
    ! !DESCRIPTION:
    ! Calculate the change in number density per unit time from the contributing
    ! rates.  These rates are not disturbance-inducing rates (that is handled
    ! elsewhere).
    !
    ! !USES:

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_cohort_type),intent(inout), target :: currentCohort
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)
    real(r8) :: frmort   ! freezing mortality rate (fraction per year)
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    integer  :: ipft              ! local copy of the pft index
    !----------------------------------------------------------------------

    ipft = currentCohort%pft
    
    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    call mortality_rates(currentCohort,bc_in,cmort,hmort,bmort,frmort)
    call LoggingMortality_frac(ipft, currentCohort%dbh, &
                               currentCohort%lmort_direct,                       &
                               currentCohort%lmort_collateral,                    &
                               currentCohort%lmort_infra )

    if (currentCohort%canopy_layer > 1)then 
       
       ! Include understory logging mortality rates not associated with disturbance
       dndt_logging = (currentCohort%lmort_direct     + &
                       currentCohort%lmort_collateral + &
                       currentCohort%lmort_infra)/hlm_freq_day

       currentCohort%dndt = -1.0_r8 * (cmort+hmort+bmort+frmort+dndt_logging) * currentCohort%n
    else
       currentCohort%dndt = -(1.0_r8 - fates_mortality_disturbance_fraction) &
            * (cmort+hmort+bmort+frmort) * currentCohort%n
    endif

    return

 end subroutine Mortality_Derivative

end module EDMortalityFunctionsMod
