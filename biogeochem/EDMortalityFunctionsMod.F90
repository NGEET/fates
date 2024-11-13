module EDMortalityFunctionsMod

  ! ============================================================================
  ! Functions that control mortality.
  ! ============================================================================

   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesGlobals          , only : fates_log
   use FatesGlobals          , only : endrun => fates_endrun
   use FatesGlobals          , only : fates_log
   use EDPftvarcon           , only : EDPftvarcon_inst
   use FatesCohortMod        , only : fates_cohort_type
   use EDTypesMod            , only : ed_site_type
   use EDParamsMod,            only : maxpft
   use EDParamsMod           , only : mort_cstarvation_model
   use FatesConstantsMod     , only : itrue,ifalse
   use FatesConstantsMod     , only : cstarvation_model_lin
   use FatesConstantsMod     , only : cstarvation_model_exp
   use FatesConstantsMod     , only : nearzero
   use FatesConstantsMod     , only : ihard_stress_decid
   use FatesConstantsMod     , only : isemi_stress_decid
   use FatesConstantsMod     , only : leaves_off
   use FatesAllometryMod     , only : bleaf
   use FatesAllometryMod     , only : storage_fraction_of_target
   use FatesInterfaceTypesMod     , only : bc_in_type
   use FatesInterfaceTypesMod     , only : hlm_use_ed_prescribed_phys
   use FatesInterfaceTypesMod     , only : hlm_freq_day
   use FatesInterfaceTypesMod     , only : hlm_use_planthydro
   use FatesInterfaceTypesMod     , only : hlm_use_tree_damage
   use EDLoggingMortalityMod , only : LoggingMortality_frac
   use EDParamsMod           , only : fates_mortality_disturbance_fraction
   use FatesConstantsMod     , only : n_landuse_cats
   use PRTGenericMod,          only : carbon12_element
   use PRTGenericMod,          only : store_organ
   use PRTParametersMod      , only : prt_params
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   
   implicit none
   private
   

   logical, parameter :: debug = .false.
   character(len=*), parameter, private :: sourcefile = &
        __FILE__
   
   public :: mortality_rates
   public :: Mortality_Derivative
   public :: ExemptTreefallDist
   
   
   ! ============================================================================
   ! 10/30/09: Created by Rosie Fisher
   ! 02/20/18: Refactored Ryan Knox
   ! ============================================================================


contains

  subroutine mortality_rates( cohort_in,bc_in, btran_ft, mean_temp,             &
      cmort,hmort,bmort, frmort,smort,asmort,dgmort )

    ! ============================================================================
    !  Calculate mortality rates from carbon storage, hydraulic cavitation, 
    !  background and freezing and size and age dependent senescence
    ! ============================================================================
    
    use FatesConstantsMod,      only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod,      only : fates_check_param_set
    use DamageMainMod,          only : GetDamageMortality
    use EDParamsmod,            only : soil_tfrz_thresh
    
    type (fates_cohort_type), intent(in) :: cohort_in 
    type (bc_in_type), intent(in) :: bc_in
    real(r8), intent(in)          :: btran_ft(maxpft) 
    real(r8), intent(in)          :: mean_temp
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality
    real(r8),intent(out) :: frmort ! freezing stress mortality
    real(r8),intent(out) :: smort  ! size dependent senescence term
    real(r8),intent(out) :: asmort ! age dependent senescence term
    real(r8),intent(out) :: dgmort ! damage dependent mortality

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: target_leaf_c      ! target leaf biomass for the current trim status and
                                   ! damage class [kgC]
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
    logical  :: is_decid_dormant   ! Flag to signal that the cohort is deciduous and dormant


    real(r8), parameter :: frost_mort_buffer = 5.0_r8  ! 5deg buffer for freezing mortality
    logical, parameter :: test_zero_mortality = .false. ! Developer test which
                                                        ! may help to debug carbon imbalances
                                                        ! and the like

    ! Check if the PFT is deciduous and leaves are completely abscised.  If this is the case,
    ! we prevent hydraulic failure mortality to occur as plants are leafless. For now, 
    ! semi-deciduous plants with senescing leaves may still die of hydraulic failure, but in
    ! the future we could accelerate senescence to avoid mortality. Note that both drought 
    ! deciduous and cold deciduous are considered here to be consistent with the idea that
    ! plants without leaves cannot die of hydraulic failure.
    is_decid_dormant =                                                            & !
       ( prt_params%stress_decid(cohort_in%pft) == ihard_stress_decid .or.        & ! Drought deciduous
         prt_params%stress_decid(cohort_in%pft) == isemi_stress_decid .or.        & ! Semi-deciduous
         prt_params%season_decid(cohort_in%pft) == itrue                  ) .and. & ! Cold deciduous
       ( cohort_in%status_coh == leaves_off )                                     ! ! Fully abscised
    
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

    ! Damage dependent mortality
    if (hlm_use_tree_damage .eq. itrue) then
       call GetDamageMortality(cohort_in%crowndamage, cohort_in%pft, dgmort)
    else
       dgmort = 0.0_r8
    end if

    
    if (hlm_use_ed_prescribed_phys .eq. ifalse) then

       ! 'Background' mortality (can vary as a function of 
       !  density as in ED1.0 and ED2.0, but doesn't here for tractability) 

       bmort = EDPftvarcon_inst%bmort(cohort_in%pft)

       ! Proxy for hydraulic failure induced mortality.
       hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)
       hf_flc_threshold = EDPftvarcon_inst%hf_flc_threshold(cohort_in%pft)

       if (hlm_use_planthydro == itrue) then
          !note the flc is set as the fraction of max conductivity in hydro
          min_fmc_ag = minval(cohort_in%co_hydr%ftc_ag(:))
          min_fmc_tr = cohort_in%co_hydr%ftc_troot
          min_fmc_ar = minval(cohort_in%co_hydr%ftc_aroot(:))
          min_fmc = min(min_fmc_ag, min_fmc_tr)
          min_fmc = min(min_fmc, min_fmc_ar)
          flc = 1.0_r8-min_fmc
          if(flc >= hf_flc_threshold .and. hf_flc_threshold < 1.0_r8 )then 
             hmort = (flc-hf_flc_threshold)/(1.0_r8-hf_flc_threshold) * &
                  EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
          else
             hmort = 0.0_r8
          endif

       else
          ! When FATES-Hydro is off, hydraulic failure mortality occurs only when btran
          ! falls below a threshold and plants have leaves.
          if ( (.not. is_decid_dormant) .and. &
               ( btran_ft(cohort_in%pft) <= hf_sm_threshold ) .and. &
               ( ( minval(bc_in%t_soisno_sl) - tfrz ) > soil_tfrz_thresh ) ) then
             hmort = EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
          else
             hmort = 0.0_r8
          end if
       end if

       ! Carbon Starvation induced mortality.
       if ( cohort_in%dbh  >  0._r8 ) then

          ! Find the current ratio between storage biomass and leaf biomass, which will be
          ! used to define carbon starvation mortality.  The reference leaf biomass is 
          ! always for when plants are fully flushed (but accounting for damage and 
          ! trimming).
          call bleaf(cohort_in%dbh,cohort_in%pft,cohort_in%crowndamage,cohort_in%canopy_trim, &
               1.0_r8, target_leaf_c)
          store_c = cohort_in%prt%GetState(store_organ,carbon12_element)
          call storage_fraction_of_target(target_leaf_c, store_c, frac)

          ! Select the carbon starvation mortality model (linear or exponential)s.
          select case (mort_cstarvation_model)
          case (cstarvation_model_lin)
             ! Linear model. Carbon starvation mortality will be zero when fraction of
             ! storage is greater than or equal to mort_upthresh_cstarvation, and will
             ! increase to the maximum mortality (mort_scalar_cstarvation) when frac = 0.
             cmort = EDPftvarcon_inst%mort_scalar_cstarvation(cohort_in%pft) * &
                max(0.0_r8, (EDPftvarcon_inst%mort_upthresh_cstarvation(cohort_in%pft)-frac) / &
                            EDPftvarcon_inst%mort_upthresh_cstarvation(cohort_in%pft) )

          case (cstarvation_model_exp)
             ! Exponential model.  Maximum carbon starvation mortality 
             ! (mort_scalar_cstarvation) occurs when frac=0. Parameter 
             ! mort_upthresh_cstarvation controls the the e-folding factor for frac. The
             ! smaller the mort_upthresh_cstarvation, the faster the mortality will decay.
             cmort = EDPftvarcon_inst%mort_scalar_cstarvation(cohort_in%pft) * &
                     exp(- frac / EDPftvarcon_inst%mort_upthresh_cstarvation(cohort_in%pft))

          case default
              write(fates_log(),*) &
                 'Invalid carbon starvation model (',mort_cstarvation_model,').'
              call endrun(msg=errMsg(sourcefile, __LINE__))
          end select

          ! Make sure the mortality is set to zero when tiny.
          if (cmort <= nearzero) then
             cmort = 0.0_r8
          end if

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

    temp_in_C = mean_temp - tfrz
    
       temp_dep_fraction  = max(0.0_r8, min(1.0_r8, 1.0_r8 - (temp_in_C - &
            EDPftvarcon_inst%freezetol(cohort_in%pft))/frost_mort_buffer) )
       frmort    = EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft) * temp_dep_fraction

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
       dgmort = 0.0_r8
    end if

    return
 end subroutine mortality_rates

 ! ============================================================================

 subroutine Mortality_Derivative( currentSite, currentCohort, bc_in, btran_ft, &
      mean_temp, land_use_label, age_since_anthro_disturbance,       &
      current_fates_landuse_state_vector, harvestable_forest_c, harvest_tag)

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
    type(fates_cohort_type),intent(inout), target :: currentCohort
    type(bc_in_type), intent(in)               :: bc_in
    real(r8),         intent(in)               :: btran_ft(maxpft)
    real(r8),         intent(in)               :: mean_temp
    integer,          intent(in)               :: land_use_label
    real(r8),         intent(in)               :: age_since_anthro_disturbance
    real(r8),         intent(in)               :: current_fates_landuse_state_vector(n_landuse_cats)
    
    real(r8), intent(in) :: harvestable_forest_c(:)   ! total carbon available for logging, kgC site-1
    integer, intent(out) :: harvest_tag(:)    ! tag to record the harvest status
                                              ! for the calculation of harvest debt in C-based
                                              ! harvest mode
                                              ! 0 - successful;
                                              ! 1 - unsuccessful since not enough carbon
                                              ! 2 - not applicable
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)
    real(r8) :: frmort   ! freezing mortality rate (fraction per year)
    real(r8) :: smort    ! size dependent senescence mortality rate (fraction per year)
    real(r8) :: asmort   ! age dependent senescence mortality rate (fraction per year)
    real(r8) :: dgmort   ! damage mortality (fraction per year)
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    integer  :: ipft              ! local copy of the pft index
   
   !----------------------------------------------------------------------

    ipft = currentCohort%pft
    
    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    call mortality_rates(currentCohort,bc_in,btran_ft, mean_temp,              &
      cmort,hmort,bmort,frmort, smort, asmort, dgmort)
    call LoggingMortality_frac(currentSite, bc_in, ipft, currentCohort%dbh, currentCohort%canopy_layer, &
                               currentCohort%lmort_direct,                       &
                               currentCohort%lmort_collateral,                    &
                               currentCohort%lmort_infra,                        &
                               currentCohort%l_degrad, &
                               bc_in%hlm_harvest_rates, &
                               bc_in%hlm_harvest_catnames, &
                               bc_in%hlm_harvest_units, &
                               land_use_label, &
                               age_since_anthro_disturbance, &
                               current_fates_landuse_state_vector, harvestable_forest_c, harvest_tag)

    if (currentCohort%canopy_layer > 1)then 
       ! Include understory logging mortality rates not associated with disturbance
       dndt_logging = (currentCohort%lmort_direct     + &
            currentCohort%lmort_collateral + &
            currentCohort%lmort_infra)/hlm_freq_day

       
       currentCohort%dndt = -1.0_r8 * &
            (cmort+hmort+bmort+frmort+smort+asmort+dgmort + dndt_logging) &
            * currentCohort%n
    else

       ! Mortality from logging in the canopy is ONLY disturbance generating, don't
       ! update number densities via non-disturbance inducing death
       ! for plants whose death is not considered disturbance (i.e. grasses),
       ! need to include all of their mortality here rather than part of it here
       ! and part in disturbance routine.

       currentCohort%dndt= -(cmort+hmort+bmort+frmort+smort+asmort+dgmort) * currentCohort%n
       if ( .not. ExemptTreefallDist(currentCohort)) then
          currentCohort%dndt = (1.0_r8-fates_mortality_disturbance_fraction) * currentCohort%dndt
       endif

    endif

    return

 end subroutine Mortality_Derivative

 ! ============================================================================

 function ExemptTreefallDist(ccohort) result(is_exempt)

   use PRTParametersMod      , only : prt_params

   ! ============================================================================
   !  Determine whether or not to consider some fraction of a cohort's crown
   !  area as disturbed patch area when individuals of that cohort die.
   ! ============================================================================

   ! Arguments
   type(fates_cohort_type),intent(in), target :: ccohort

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
