module DamageMainMod

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : i4 => fates_int
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : years_per_day
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use EDPftvarcon           , only : EDPftvarcon_inst
  use EDParamsMod           , only : damage_event_code
  use EDtypesMod            , only : ed_site_type
  use EDtypesMod            , only : ed_patch_type
  use EDtypesMod            , only : ed_cohort_type
  use EDtypesMod            , only : AREA
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : SetState
  use FatesInterfaceTypesMod, only : hlm_current_day
  use FatesInterfaceTypesMod, only : hlm_current_month
  use FatesInterfaceTypesMod, only : hlm_current_year
  use FatesInterfaceTypesMod, only : hlm_model_day
  use FatesInterfaceTypesMod , only : hlm_day_of_year
  
  
  
  implicit none
  private

  logical, protected :: damage_time  ! if true then damage occurs during current time step

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
  public :: GetCrownReduction
  public :: GetDamageFrac
  public :: IsItDamageTime
  public :: damage_time
  public :: GetDamageMortality
  public :: DamageRecovery
  
  logical :: debug = .false.  ! for debugging


  ! The following is the special classification for undamaged plants
  ! and is used in contexts where cohort%damageclass is used. This is
  ! to flag to the user that an undamaged plant is assumed in those contexts
  
  integer, parameter, public :: undamaged_class = 1

  
  ! ============================================================================
  ! ============================================================================

contains

  subroutine DamageRecovery(currentCohort,recoveryCohort)

    !---------------------------------------------------------------------------
    ! JN March 2021
    ! At this point it is possible that damaged cohorts have reached their
    ! target allometries. There is a choice now - if they have excess carbon,
    ! they can use it to grow along their reduced allometric targets  - i.e.
    ! dbh and all carbon pools grow out together. OR they can use excess carbon to
    ! jump to a lower damage class by changing their target allometry and growing 
    ! to meet new C pools for same dbh.
    !
    ! d = damage class
    ! --------------------------------------------------------------------------
    
    type(ed_cohort_type) :: currentCohort
    type(ed_cohort_type), pointer :: recoveryCohort



    if (crowndamage > 1 .and. carbon_balance > calloc_abs_error) then

       if(damage_recovery_scalar > 0.0_r8) then 
       ! 1. What is excess carbon?
       ! carbon_balance
          
       !  2. What is biomass required to go from current damage level to next damage level?

       ! mass of this damage class
       mass_d = (sum(leaf_c(1:nleafage)) + fnrt_c + store_c + sapw_c + struct_c ) 
  
       ! Target sapwood biomass according to allometry and trimming [kgC]
       call bsap_allom(dbh,ipft, crowndamage-1, canopy_trim,sapw_area,targetn_sapw_c)
       ! Target total above ground biomass in woody/fibrous tissues  [kgC]
       call bagw_allom(dbh,ipft, crowndamage-1, targetn_agw_c)
       ! Target total below ground biomass in woody/fibrous tissues [kgC] 
       call bbgw_allom(dbh,ipft,targetn_bgw_c)
       ! Target total dead (structrual) biomass [kgC]
       call bdead_allom( targetn_agw_c, targetn_bgw_c, targetn_sapw_c, ipft, targetn_struct_c)
       ! Target fine-root biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
       call bfineroot(dbh,ipft,canopy_trim,targetn_fnrt_c)
       ! Target storage carbon [kgC,kgC/cm]
       call bstore_allom(dbh,ipft,crowndamage-1, canopy_trim,targetn_store_c)
       ! Target leaf biomass according to allometry and trimming
       if(leaf_status==2) then
          call bleaf(dbh,ipft,crowndamage-1, canopy_trim,targetn_leaf_c)
       else
          targetn_leaf_c = 0._r8
       end if


       mass_dminus1 = (max(sum(leaf_c), targetn_leaf_c) + max(fnrt_c, targetn_fnrt_c) + &
            max(store_c, targetn_store_c) + max(sapw_c, targetn_sapw_c) + &
            max(struct_c, targetn_struct_c)) 

       ! Carbon needed to get from current mass to allometric target mass of next damage class up
       recovery_demand = mass_dminus1 - mass_d
       
       ! 3. How many trees can get there with excess carbon?
       max_recover_n =  carbon_balance * n / recovery_demand 

       ! 4. Use the scalar to decide how many to recover
       n_recover = max_recover_n * damage_recovery_scalar

       ! carbon balance needs to be updated

       ! there is a special case where damage_recovery_scalar = 1, but
       ! max_recover_n > n (i.e. there is more carbon than needed for all
       ! individuals to recover to the next damage class.
       ! in this case we can cheat, by making n_recover 0 and simply
       ! allowing the donor cohort to recover and then go through
       ! prt - will this work though? if they are not anywhere near allometry?
       

       if(damage_recovery_scalar .eq. 1.0_r8 .and. max_recover_n > n) then
          n_recover = 0.0_r8
          crowndamage = crowndamage - 1
          ! call prt from within itself here? 
       else
          carbon_balance = (n * carbon_balance - (recovery_demand * n_recover)) /(n-n_recover)
       end if

       ! we reduce number density here and continue on with daily prt for the
       ! part of the cohort that is not recovering - staying fixed on its
       ! current reduced allometries
       n = n - n_recover

       ! Outside of parteh we will copy the cohort and allow the
       ! recovery portion to change allometric targets.
       
       end if ! end if some recovery is permited
    end if ! end if crowndamage 
    !------------------------------------------------------------------------------------




    
    if(currentCohort%crowndamage > 1) then

       ! N is inout boundary condition so has now been updated. The difference must
       ! go to a new cohort
       n_recover = n_old - currentCohort%n

       if(n_recover > nearzero) then

          allocate(nc)
          if(hlm_use_planthydro .eq. itrue) call InitHydrCohort(CurrentSite,nc)
          ! Initialize the PARTEH object and point to the
          ! correct boundary condition fields
          nc%prt => null()
          call InitPRTObject(nc%prt)
          call InitPRTBoundaryConditions(nc)          
          !    call zero_cohort(nc)  
          call copy_cohort(currentCohort, nc)

          nc%n = n_recover
          nc%crowndamage = currentCohort%crowndamage - 1

          ! Need to adjust the crown area which is NOT on a per individual basis
          nc%c_area = nc%n/n_old * currentCohort%c_area
          currentCohort%c_area = currentCohort%c_area - nc%c_area

          ! This new cohort spends carbon balance on growing out pools
          ! (but not dbh) to reach new allometric targets
          ! This was already calculated within parteh - this cohort should just
          ! be able to hit allometric targets of one damage class down
          call nc%prt%DamageRecovery()

          ! at this point we need to update fluxes or this cohort will
          ! fail its mass conservation checks

          sapw_c   = nc%prt%GetState(sapw_organ, all_carbon_elements)
          struct_c = nc%prt%GetState(struct_organ, all_carbon_elements)
          leaf_c   = nc%prt%GetState(leaf_organ, all_carbon_elements)
          fnrt_c   = nc%prt%GetState(fnrt_organ, all_carbon_elements)
          store_c  = nc%prt%GetState(store_organ, all_carbon_elements)
          repro_c  = nc%prt%GetState(repro_organ, all_carbon_elements)
          nc_carbon = sapw_c + struct_c + leaf_c + fnrt_c + store_c + repro_c

          cc_sapw_c   = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
          cc_struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
          cc_leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
          cc_fnrt_c   = currentCohort%prt%GetState(fnrt_organ, all_carbon_elements)
          cc_store_c  = currentCohort%prt%GetState(store_organ, all_carbon_elements)
          cc_repro_c  = currentCohort%prt%GetState(repro_organ, all_carbon_elements)
          cc_carbon = cc_sapw_c + cc_struct_c + cc_leaf_c + cc_fnrt_c + cc_store_c + cc_repro_c


          call PRTDamageRecoveryFluxes(nc%prt, leaf_organ, leaf_c0, leaf_c, cc_leaf_c)
          call PRTDamageRecoveryFluxes(nc%prt, repro_organ, repro_c0, repro_c, cc_repro_c)
          call PRTDamageRecoveryFluxes(nc%prt, sapw_organ, sapw_c0, sapw_c, cc_sapw_c)
          call PRTDamageRecoveryFluxes(nc%prt, struct_organ, struct_c0, struct_c, cc_struct_c)
          call PRTDamageRecoveryFluxes(nc%prt, store_organ, store_c0, store_c, cc_store_c)
          call PRTDamageRecoveryFluxes(nc%prt, fnrt_organ, fnrt_c0, fnrt_c, cc_fnrt_c)

          ! update crown area
          call carea_allom(nc%dbh, nc%n, currentSite%spread, nc%pft, nc%crowndamage, nc%c_area)
          call carea_allom(currentCohort%dbh, currentCohort%n, currentSite%spread, &
               currentCohort%pft, currentCohort%crowndamage, currentCohort%c_area)



          !----------- Insert copy into linked list ----------------------! 
          nc%shorter => currentCohort
          if(associated(currentCohort%taller))then
             nc%taller => currentCohort%taller
             currentCohort%taller%shorter => nc
          else
             currentPatch%tallest => nc    
             nc%taller => null()
          endif
          currentCohort%taller => nc

       end if ! end if greater than nearzero

    end if ! end if crowndamage > 1




    return
  end subroutine DamageRecovery


  subroutine IsItDamageTime(is_master, currentSite)

    !----------------------------------------------------------------------------
    ! This subroutine determines whether damage should occur (it is called daily)
    ! This is almost an exact replica of the IsItLoggingTime subroutine
    !-----------------------------------------------------------------------------


    integer, intent(in) :: is_master
    type(ed_site_type), intent(inout), target :: currentSite

    integer :: icode     ! Integer equivalent of the event code (parameter file only allows reals)
    integer :: damage_date  ! Day of month for damage extracted from event code
    integer :: damage_month ! Month of year for damage extracted from event code
    integer :: damage_year  ! Year for damage extracted from event code
    integer :: model_day_int  ! Model day
    
    character(len=64) :: fmt = '(a,i2.2,a,i2.2,a,i4.4)'

    damage_time = .false.
    icode = int(damage_event_code)

    model_day_int = nint(hlm_model_day)
    
    if(icode .eq. 1) then
       ! Damage is turned off 
       damage_time = .false.

    else if(icode .eq. 2) then
       ! Damage event on first time step 
       if(model_day_int .eq.1) then
          damage_time = .true.
       end if

    else if(icode .eq. 3) then
       ! Damage event every day - not sure this is recommended as it will result in a very large
       ! number of cohorts 
       damage_time = .true.

    else if(icode .eq. 4) then
       ! Damage event once a month
       if(hlm_current_day.eq.1 ) then
          damage_time = .true.
       end if

    else if(icode < 0 .and. icode > -366) then
       ! Damage event every year on a specific day of the year
       if(hlm_day_of_year .eq. abs(icode) ) then
          damage_time = .true.
       end if

    else if(icode > 10000 ) then
       ! Specific Event: YYYYMMDD
       damage_date  = icode - int(100* floor(real(icode)/100))
       damage_year  = floor(real(icode)/10000)
       damage_month = floor(real(icode)/100) - damage_year*100

       if(hlm_current_day .eq. damage_date .and. &
            hlm_current_month .eq. damage_month .and. &
            hlm_current_year .eq. damage_year ) then
          damage_time = .true.
       end if

    else
       ! Bad damage event flag
       write(fates_log(),*) 'An invalid damage code was specified in fates_params'
       write(fates_log(),*) 'Check DamageMainMod.F90:IsItDamageTime()'
       write(fates_log(),*) 'for a breakdown of the valide codes and change'
       write(fates_log(),*) 'fates_damage_event_code in the file accordingly.'
       write(fates_log(),*) 'exiting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if(damage_time .and. (is_master.eq.itrue) ) then
       write(fates_log(),fmt) 'Damage Event Enacted on date: ', &
            hlm_current_month,'-', hlm_current_day,'-',hlm_current_year
    end if
  
    return
    
  end subroutine IsItDamageTime
  
  !----------------------------------------------------------------------------

  subroutine GetDamageFrac(cc_cd, nc_cd, pft, dist_frac)


    ! given current cohort damage class find the fraction of individuals
    ! going to the new damage class.
    ! Consults a look up table of transitions from param derived. 
    
    ! USES
    use FatesInterfaceTypesMod, only : nlevdamage
    use FatesConstantsMod, only : years_per_day
    use FatesParameterDerivedMod, only : param_derived

       
    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    integer, intent(in) :: pft
    real(r8), intent(out) :: dist_frac             ! probability of current cohort moving to new damage level

    dist_frac = param_derived%damage_transitions(cc_cd, nc_cd, pft) !* years_per_day
    ! (if damage is occuring annually don't do this)
    

  end subroutine GetDamageFrac

  !-------------------------------------------------------    
  
  subroutine GetCrownReduction(crowndamage, crown_reduction)

    !------------------------------------------------------------------                                                                     
    ! This function takes the crown damage class of a cohort (integer)
    ! and returns the fraction of the crown that is lost
    ! Since crowndamage class = 1 means no damage, we subtract one     
    ! before multiplying by 0.2                                    
    ! Therefore, first damage class is 20% loss of crown, second 40% etc.                                                                   
    !-------------------------------------------------------------------                                                                    
    use FatesInterfaceTypesMod     , only : nlevdamage

    integer(i4), intent(in)   :: crowndamage
    real(r8),    intent(out)  :: crown_reduction

    ! local variables
    real(r8) :: class_width

    class_width = 1.0_r8/nlevdamage
    crown_reduction = min(1.0_r8, (real(crowndamage) - 1.0_r8) * class_width)

    return
  end subroutine GetCrownReduction


  !----------------------------------------------------------------------------------------


  subroutine GetDamageMortality(crowndamage,pft, dgmort)

    use FatesInterfaceTypesMod     , only : nlevdamage
    use EDPftvarcon                , only : EDPftvarcon_inst
    
    integer(i4), intent(in) :: crowndamage
    integer(i4), intent(in) :: pft
    real(r8),    intent(out) :: dgmort

    ! local variables
    real(r8) :: damage_mort_p1
    real(r8) :: damage_mort_p2
    real(r8) :: class_width
    real(r8) :: crown_loss
    
    class_width = 1.0_r8/real(nlevdamage)
    
    ! parameter to determine slope of exponential
    damage_mort_p1 = EDPftvarcon_inst%damage_mort_p1(pft)
    damage_mort_p2 = EDPftvarcon_inst%damage_mort_p2(pft)

    ! make damage mortality a function of crownloss and not crowndamage
    ! class so that it doesn't need to be re-parameterised if the number
    ! of damage classes change.
    crown_loss = min(1.0_r8, (real(crowndamage) - 1.0_r8) * class_width)
    
    if (crowndamage .eq. 1 ) then
       dgmort = 0.0_r8
    else
       dgmort = 1.0_r8 / (1.0_r8 + exp(-1.0_r8 * damage_mort_p2 * &
            (crown_loss - damage_mort_p1) ) )

    end if

    return
  end subroutine GetDamageMortality
  !----------------------------------------------------------------------------------------

 
end module DamageMainMod

