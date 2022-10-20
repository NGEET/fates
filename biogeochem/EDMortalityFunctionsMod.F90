module EDMortalityFunctionsMod

  ! ============================================================================
  ! Functions that control mortality.
  ! ============================================================================

   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesGlobals          , only : fates_log
   use FatesGlobals          , only : endrun => fates_endrun
   use FatesGlobals          , only : fates_log
   use EDPftvarcon           , only : EDPftvarcon_inst
   use EDTypesMod            , only : ed_cohort_type
   use EDTypesMod            , only : ed_site_type
   use EDTypesMod            , only : ed_patch_type
   use EDTypesMod            , only : leaves_on
   use FatesConstantsMod     , only : itrue,ifalse
   use FatesAllometryMod     , only : bleaf
   use FatesAllometryMod     , only : storage_fraction_of_target
   use FatesInterfaceTypesMod     , only : bc_in_type
   use FatesInterfaceTypesMod     , only : hlm_use_ed_prescribed_phys
   use FatesInterfaceTypesMod     , only : hlm_freq_day
   use FatesInterfaceTypesMod     , only : hlm_use_planthydro
   use FatesInterfaceTypesMod     , only : hlm_model_day
   use FatesInterfaceTypesMod     , only : hlm_use_frosthard
   use FatesInterfaceTypesMod     , only : hlm_current_year
   use FatesInterfaceTypesMod     , only : hlm_current_month
   use FatesInterfaceTypesMod     , only : hlm_current_day
   use PRTParametersMod           , only : prt_params
   use EDLoggingMortalityMod , only : LoggingMortality_frac
   use EDParamsMod           , only : fates_mortality_disturbance_fraction

   use PRTGenericMod,          only : all_carbon_elements
   use PRTGenericMod,          only : store_organ
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   
   implicit none
   private
   

   logical, parameter :: debug = .false.
   character(len=*), parameter, private :: sourcefile = &
        __FILE__
   
   public :: mortality_rates
   public :: Mortality_Derivative
   public :: ExemptTreefallDist
   public :: Hardening_scheme
   
   
   ! ============================================================================
   ! 10/30/09: Created by Rosie Fisher
   ! 02/20/18: Refactored Ryan Knox
   ! ============================================================================


contains



  subroutine mortality_rates( currentSite,cohort_in,bc_in,cmort,hmort,bmort,frmort,smort,asmort )

    ! ============================================================================
    !  Calculate mortality rates from carbon storage, hydraulic cavitation, 
    !  background and freezing and size and age dependent senescence
    ! ============================================================================
    
    use FatesConstantsMod,  only : tfrz => t_water_freeze_k_1atm
    use FatesInterfaceTypesMod        , only : hlm_hio_ignore_val   
    use FatesConstantsMod,  only : fates_check_param_set
    
    type (ed_cohort_type), intent(in) :: cohort_in 
    type (ed_site_type), intent(inout), target  :: currentSite
    type (bc_in_type), intent(in) :: bc_in
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality
    real(r8),intent(out) :: frmort ! freezing stress mortality
    real(r8),intent(out) :: smort  ! size dependent senescence term
    real(r8),intent(out) :: asmort ! age dependent senescence term 

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: leaf_c_target      ! target leaf biomass kgC
    real(r8) :: store_c
    real(r8) :: hf_sm_threshold    ! hydraulic failure soil moisture threshold 
    real(r8) :: hf_flc_threshold   ! hydraulic failure fractional loss of conductivity threshold
    real(r8) :: mort_ip_size_senescence ! inflection point for increase in mortality with dbh 
    real(r8) :: mort_r_size_senescence  ! rate of mortality increase with dbh in senesence term
    real(r8) :: mort_ip_age_senescence ! inflection point for increase in mortality with age
    real(r8) :: mort_r_age_senescence ! rate of mortality increase with age in senescence term 
    real(r8) :: temp_dep_fraction  ! Temp. function (freezing mortality)
    real(r8) :: temp_in_C          ! Daily averaged temperature in Celcius
    real(r8) :: min_fmc_ag         ! minimum fraction of maximum conductivity for aboveground
    real(r8) :: min_fmc_tr         ! minimum fraction of maximum conductivity for transporting root
    real(r8) :: min_fmc_ar         ! minimum fraction of maximum conductivity for absorbing root
    real(r8) :: min_fmc            ! minimum fraction of maximum conductivity for whole plant
    real(r8) :: flc                ! fractional loss of conductivity 
    real(r8) :: Tmin               ! minimum daily average temperature
    real(r8) :: hard_hydr          ! reduced hydro mortality if hardened 
    real(r8) :: max_h                      !maximum hardiness level
    real(r8), parameter :: min_h = -2.0_r8 !minimum hardiness level
    real(r8), parameter :: frost_mort_buffer = 5.0_r8  ! 5deg buffer for freezing mortality
    logical, parameter :: test_zero_mortality = .false. ! Developer test which
                                                        ! may help to debug carbon imbalances
                                                        ! and the like
     
   ! Size Dependent Senescence
    ! rate (r) and inflection point (ip) define the increase in mortality rate with dbh
    mort_r_size_senescence = EDPftvarcon_inst%mort_r_size_senescence(cohort_in%pft)
    mort_ip_size_senescence = EDPftvarcon_inst%mort_ip_size_senescence(cohort_in%pft)
    
    ! if param values have been set then calculate smort
    if ( mort_ip_size_senescence < fates_check_param_set ) then 
       smort = 1.0_r8 / ( 1.0_r8 + exp( -1.0_r8 * mort_r_size_senescence * &
            (cohort_in%dbh - mort_ip_size_senescence) ) ) 
    else
       smort = 0.0_r8
    end if

    ! if param values have been set then calculate asmort

    

    mort_r_age_senescence = EDPftvarcon_inst%mort_r_age_senescence(cohort_in%pft)
    mort_ip_age_senescence = EDPftvarcon_inst%mort_ip_age_senescence(cohort_in%pft)
    
    if ( mort_ip_age_senescence < fates_check_param_set ) then
       ! Age Dependent Senescence
       ! rate and inflection point define the change in mortality with age
       
       asmort = 1.0_r8 / (1.0_r8 + exp(-1.0_r8 * mort_r_age_senescence * &
            (cohort_in%coage - mort_ip_age_senescence ) ) )
    else
       asmort = 0.0_r8
    end if


    
if (hlm_use_ed_prescribed_phys .eq. ifalse) then

    ! 'Background' mortality (can vary as a function of 
    !  density as in ED1.0 and ED2.0, but doesn't here for tractability) 
       
       bmort = EDPftvarcon_inst%bmort(cohort_in%pft)

    ! Proxy for hydraulic failure induced mortality. 
    hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)
    hf_flc_threshold = EDPftvarcon_inst%hf_flc_threshold(cohort_in%pft)
    if(hlm_use_planthydro.eq.itrue)then
     !note the flc is set as the fraction of max conductivity in hydro
     min_fmc_ag = minval(cohort_in%co_hydr%ftc_ag(:))
     min_fmc_tr = cohort_in%co_hydr%ftc_troot
     min_fmc_ar = minval(cohort_in%co_hydr%ftc_aroot(:))
     min_fmc = min(min_fmc_ag, min_fmc_tr)
     min_fmc = min(min_fmc, min_fmc_ar)
     flc = 1.0_r8-min_fmc
     if(flc >= hf_flc_threshold .and. hf_flc_threshold < 1.0_r8 )then
       if (hlm_use_hydrohard .eq. itrue .and. currentSite%hard_level2(cohort_in%pft) < -3._r8) then  
         max_h=min(max(EDPftvarcon_inst%freezetol(cohort_in%pft),max(currentSite%hardtemp,-60._r8)-10._r8),min_h)
         hard_hydr=EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)*( ((currentSite%hard_level2(cohort_in%pft)-max_h)/(min_h-max_h) )  *0.5_r8+0.5_r8)
       else 
         hard_hydr=EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
       end if
       hmort = (flc-hf_flc_threshold)/(1.0_r8-hf_flc_threshold) * hard_hydr
     else
       hmort = 0.0_r8
     endif      
    else
     if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= hf_sm_threshold)then 
       hmort = EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
     else
       hmort = 0.0_r8
     endif
    endif 
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then

       call bleaf(cohort_in%dbh,cohort_in%pft,cohort_in%canopy_trim,leaf_c_target)
       store_c = cohort_in%prt%GetState(store_organ,all_carbon_elements)

       call storage_fraction_of_target(leaf_c_target, store_c, frac)
       if( frac .lt. 1._r8) then
          cmort = max(0.0_r8,EDPftvarcon_inst%mort_scalar_cstarvation(cohort_in%pft) * &
               (1.0_r8 - frac))
       else
          cmort = 0.0_r8
       endif

    else
       write(fates_log(),*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    !-------------------------------------------------------------------------------- 
    !    Mortality due to cold and freezing stress (frmort), based on ED2 and:           
    !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions 
    !           of land-use change, CO2 fertilization, and climate variability to the    
    !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              
    !           doi: 10.1111/j.1365-2486.2006.01254.x                                    

    if (hlm_use_frosthard .eq. itrue) then !Hardiness level dependent frost mortality
       Tmin=bc_in%tmin24_si-273.15_r8
       temp_dep_fraction  = max(0.0_r8, min(1.0_r8,(-Tmin + &
       max(currentSite%hard_level2(cohort_in%pft),EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft)))/frost_mort_buffer))
       if (nint(hlm_model_day)<185) then
          temp_dep_fraction=0._r8
       endif
       frmort    = EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft)*temp_dep_fraction
    else
       temp_in_C = bc_in%t_veg24_pa(ifp) - tfrz
       temp_dep_fraction  = max(0.0_r8, min(1.0_r8, 1.0_r8 - (temp_in_C - &
                            EDPftvarcon_inst%freezetol(cohort_in%pft))/frost_mort_buffer) )
       frmort    = EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft) * temp_dep_fraction
    endif


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
       smort = 0.0_r8
       asmort = 0.0_r8
    end if
       
    return
 end subroutine mortality_rates

 ! ============================================================================

 subroutine Mortality_Derivative( currentSite, currentCohort, bc_in, frac_site_primary)

    !
    ! !DESCRIPTION:
    ! Calculate the change in number density per unit time from the contributing
    ! rates.  These rates are not disturbance-inducing rates (that is handled
    ! elsewhere).
    !
    ! !USES:

    use FatesInterfaceTypesMod, only : hlm_freq_day
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_cohort_type),intent(inout), target :: currentCohort
    type(bc_in_type), intent(in)               :: bc_in
    real(r8), intent(in)                       :: frac_site_primary
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)
    real(r8) :: frmort   ! freezing mortality rate (fraction per year)
    real(r8) :: smort    ! size dependent senescence mortality rate (fraction per year)
    real(r8) :: asmort   ! age dependent senescence mortality rate (fraction per year)
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    integer  :: ipft              ! local copy of the pft index
    !----------------------------------------------------------------------

    ipft = currentCohort%pft
    
    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    call mortality_rates(currentSite,currentCohort,bc_in,cmort,hmort,bmort,frmort,smort, asmort)
    call LoggingMortality_frac(ipft, currentCohort%dbh, currentCohort%canopy_layer, &
                               currentCohort%lmort_direct,                       &
                               currentCohort%lmort_collateral,                    &
                               currentCohort%lmort_infra,                        &
                               currentCohort%l_degrad, &
                               bc_in%hlm_harvest_rates, &
                               bc_in%hlm_harvest_catnames, &
                               bc_in%hlm_harvest_units, &
                               currentCohort%patchptr%anthro_disturbance_label, &
                               currentCohort%patchptr%age_since_anthro_disturbance, &
                               frac_site_primary)

    
    

    if (currentCohort%canopy_layer > 1)then 
       ! Include understory logging mortality rates not associated with disturbance
       dndt_logging = (currentCohort%lmort_direct     + &
            currentCohort%lmort_collateral + &
            currentCohort%lmort_infra)/hlm_freq_day

       
       currentCohort%dndt = -1.0_r8 * &
            (cmort+hmort+bmort+frmort+smort+asmort + dndt_logging) &
            * currentCohort%n
    else

       ! Mortality from logging in the canopy is ONLY disturbance generating, don't
       ! update number densities via non-disturbance inducing death

       ! for plants whose death is not considered disturbance (i.e. grasses),
       ! need to include all of their mortality here rather than part of it here
       ! and part in disturbance routine.

       currentCohort%dndt= -(cmort+hmort+bmort+frmort+smort+asmort) * currentCohort%n
       if ( .not. ExemptTreefallDist(currentCohort)) then
          currentCohort%dndt = (1.0_r8-fates_mortality_disturbance_fraction) * currentCohort%dndt
       endif

    endif

    return

 end subroutine Mortality_Derivative
! ============================================================================

  subroutine Hardening_scheme(currentSite,currentPatch,cohort_in,bc_in )
    !
    ! !DESCRIPTION:
    ! Hardening module from Rammig et al. 2010. 
    ! Controls the unrealistic root water release by reducing root confuctivity, controls freezing mortality

    ! !ARGUMENTS
    type (ed_site_type), intent(inout), target  :: currentSite
    type (ed_patch_type), intent(inout), target  :: currentPatch
    type (ed_cohort_type), intent(inout) :: cohort_in
    type (bc_in_type), intent(in) :: bc_in


    ! !LOCAL VARIABLES:
    real(r8) :: Tmean                        ! Daily average temperature °C
    real(r8) :: Tmin                         ! Daily minimum temperature °C
    real(r8) :: max_h                        ! Maximum hardiness level
    real(r8) :: max_h_dehard                 ! Maximum hardiness level for dehardening function
    real(r8), parameter :: min_h = -2.0_r8   ! Minimum hardiness level from Bigras for Picea abies (°C)
    real(r8) :: target_h                     ! Target hardiness 
    real(r8) :: rate_h                       ! Hardening rate     
    real(r8) :: rate_dh                      ! Dehardening rate
    real(r8) :: hard_level_prev              ! Temporary variable for the previous time-step hardiness level
    real(r8) :: dayl_thresh                  ! Threshold for the onset of the hardening only period
    integer  :: ipft                         ! pft index

    Tmean=bc_in%temp24_si-273.15_r8
    if (Tmean<-200.0_r8) then
       Tmean=0.0_r8
    endif
    Tmin=bc_in%tmin24_si-273.15_r8

    max_h=min(max(EDPftvarcon_inst%freezetol(cohort_in%pft),max(currentSite%hardtemp,-60._r8)-10._r8),min_h)

    !Calculation of the target hardiness
    if (Tmean <= max_h/1.5_r8) then !
       target_h=max_h
    else if (Tmean>= (6.0_r8-max_h/6._r8)) then
       target_h=min_h
    else
       target_h = -sin((3.14159_r8*(0.5_r8+(Tmean-max_h/1.5_r8)/(-max_h/1.5_r8+(6.0_r8-max_h/6._r8)))))*(min_h-max_h)/2._r8-(min_h-max_h)/2._r8+min_h
    end if
    !Calculation of the hardening rate
    if (Tmean <= max_h/2) then
       rate_h=(max_h-min_h)/-31.11_r8+0.1_r8
    else if (Tmean >= 20.0_r8) then
       rate_h=0.1_r8
    else 
       rate_h = sin((3.14159_r8*(0.5+(Tmean-max_h/2._r8)/(-max_h/2._r8+20._r8))))*((max_h-min_h)/-62.22_r8)+(((max_h-min_h)/-62.22_r8)+0.1_r8)
    end if
    !Calculation of the dehardening rate
    if (max_h>-5.111_r8) then !this if steatement so that there is always possibility for 0.5 dehardening rate at 12 degreesC.
       max_h_dehard=-5.111_r8
    else 
       max_h_dehard=max_h
    end if
    ipft = cohort_in%pft
    if (prt_params%season_decid(ipft) == itrue) then
     if (Tmean <= 5.0_r8) then
       rate_dh=0.0_r8
     else if (Tmean >= 15.0_r8) then
       rate_dh=5.0_r8*(max_h_dehard-min_h)/-31.11_r8
     else
       rate_dh=(Tmean-5.0_r8)*((max_h_dehard-min_h)/-62.22_r8)
     end if      
    else
     if (Tmean <= 2.5_r8) then
       rate_dh=0.0_r8
     else if (Tmean >= 12.5_r8) then
       rate_dh=5.0_r8*(max_h_dehard-min_h)/-31.11_r8
     else
       rate_dh=(Tmean-2.5_r8)*((max_h_dehard-min_h)/-62.22_r8)
     end if   
    end if

    dayl_thresh= 42000.0_r8 + ( (-30.0_r8 - max(-60.0_r8,min(0.0_r8,currentSite%hardtemp)) )/15.0_r8) * 4500.0_r8  
    !Hardening calculation
    cohort_in%hard_level_prev = cohort_in%hard_level

    if (cohort_in%hard_level_prev + rate_dh > min_h) then
       cohort_in%hard_level = min_h
    else if (cohort_in%hard_level_prev >= target_h) then
       cohort_in%hard_level = cohort_in%hard_level_prev - rate_h
    else if (cohort_in%hard_level_prev < target_h) then
       cohort_in%hard_level = cohort_in%hard_level_prev + rate_dh
    end if
    if (bc_in%dayl_si <= dayl_thresh .and. bc_in%dayl_si < bc_in%prev_dayl_si) then
       if (cohort_in%hard_level_prev >= target_h) then
       cohort_in%hard_level = cohort_in%hard_level_prev - rate_h
       else 
          cohort_in%hard_level = cohort_in%hard_level_prev ! now dehardening is not possible in autumn but the hardiness level is also allowed to remain steady
       end if
    else if (bc_in%dayl_si > dayl_thresh .and. bc_in%dayl_si < bc_in%prev_dayl_si) then
          cohort_in%hard_level = min_h
    end if  
    ipft = cohort_in%pft
    if (prt_params%season_decid(ipft) == itrue .and. cohort_in%status_coh == leaves_on .and. bc_in%dayl_si > bc_in%prev_dayl_si .and. &
        (nint(hlm_model_day) >= currentSite%cleafondate .or. nint(hlm_model_day) >= currentSite%dleafondate)) then         
       cohort_in%hard_level = min_h
    end if
    if (cohort_in%hard_level > min_h) then
       cohort_in%hard_level = min_h
    end if
    if (cohort_in%hard_level < max_h) then
       cohort_in%hard_level = max_h
    end if 

    return

  end subroutine Hardening_scheme
  
 ! ============================================================================

 function ExemptTreefallDist(ccohort) result(is_exempt)

   use PRTParametersMod      , only : prt_params

   ! ============================================================================
   !  Determine whether or not to consider some fraction of a cohort's crown
   !  area as disturbed patch area when individuals of that cohort die.
   ! ============================================================================

   ! Arguments
   type(ed_cohort_type),intent(in), target :: ccohort

   logical :: is_exempt ! if true, then treat all mortality from this cohort as non-disturbance-generating

   !! current logic is only to exempt non-woody plants

   if ( prt_params%woody(ccohort%pft) == ifalse ) then
      is_exempt = .true.
   else
      is_exempt = .false.
   endif

   return

 end function ExemptTreefallDist

end module EDMortalityFunctionsMod
