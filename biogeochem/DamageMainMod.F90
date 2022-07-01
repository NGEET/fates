module DamageMainMod

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : i4 => fates_int
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : years_per_day
  use FatesConstantsMod     , only : nearzero
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use EDPftvarcon           , only : EDPftvarcon_inst
  use EDParamsMod           , only : damage_event_code
  use EDTypesMod            , only : ed_site_type
  use EDTypesMod            , only : ed_patch_type
  use EDTypesMod            , only : ed_cohort_type
  use EDTypesMod            , only : AREA
  use EDTypesMod            , only : leaves_on
  use PRTGenericMod,          only : num_elements
  use PRTGenericMod,          only : element_list
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : nitrogen_element
  use PRTGenericMod,          only : phosphorus_element
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
  use PRTGenericMod,          only : StorageNutrientTarget
  use PRTParametersMod,       only : prt_params
  use FatesInterfaceTypesMod, only : hlm_current_day
  use FatesInterfaceTypesMod, only : hlm_current_month
  use FatesInterfaceTypesMod, only : hlm_current_year
  use FatesInterfaceTypesMod, only : hlm_model_day
  use FatesInterfaceTypesMod, only : hlm_day_of_year
  use FatesInterfaceTypesMod, only : hlm_use_planthydro
  use FatesAllometryMod,      only : bsap_allom
  use FatesAllometryMod,      only : bagw_allom
  use FatesAllometryMod,      only : bbgw_allom
  use FatesAllometryMod,      only : bdead_allom
  use FatesAllometryMod,      only : bfineroot
  use FatesAllometryMod,      only : bstore_allom
  use FatesAllometryMod,      only : bleaf
  use EDCohortDynamicsMod,   only : copy_cohort
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use EDCohortDynamicsMod  , only : InitPRTObject
  use EDCohortDynamicsMod  , only : InitPRTBoundaryConditions
  
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

  subroutine DamageRecovery(csite,cpatch,ccohort,newly_recovered)

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

    type(ed_site_type)   :: csite            ! Site of the current cohort
    type(ed_patch_type)  :: cpatch           ! patch of the current cohort
    type(ed_cohort_type),pointer :: ccohort  ! Current (damaged) cohort
    logical              :: newly_recovered  ! true if we create a new cohort

    ! locals
    type(ed_cohort_type), pointer :: rcohort ! New cohort that recovers by
                                             ! having an lower damage class
    real(r8) :: sapw_area                    ! sapwood area
    real(r8) :: target_sapw_c,target_sapw_m  ! sapwood mass, C and N/P
    real(r8) :: target_agw_c                 ! target above ground wood
    real(r8) :: target_bgw_c                    ! target below ground wood
    real(r8) :: target_struct_c,target_struct_m ! target structural C and N/P
    real(r8) :: target_fnrt_c,target_fnrt_m     ! target fine-root C and N/P
    real(r8) :: target_leaf_c,target_leaf_m     ! target leaf C and N/P
    real(r8) :: target_store_c,target_store_m   ! target storage C and N/P
    real(r8) :: target_repro_m                  ! target reproductive C/N/P
    real(r8) :: leaf_m,fnrt_m,sapw_m            ! actual masses in organs C/N/P
    real(r8) :: struct_m,store_m,repro_m        ! actual masses in organs C/N/P
    real(r8) :: mass_d                          ! intermediate term for nplant_recover
    real(r8) :: mass_dminus1                    ! intermediate term for nplant_recover
    real(r8) :: available_m                     ! available mass that can be used to 
                                                ! improve damage class
    real(r8) :: recovery_demand                 ! amount of mass needed to get to 
                                                ! get to the target of the next damage class
    real(r8) :: max_recover_nplant              ! max number of plants that could get to
                                                ! target of next class
    real(r8) :: nplant_recover                  ! number of plants in cohort that will
                                                ! recover to the next class
    integer  :: el                                ! element loop counter
    
    associate(dbh => ccohort%dbh, &
         ipft => ccohort%pft, &
         canopy_trim => ccohort%canopy_trim)

      ! If we are currently undamaged, no recovery
      ! necessary, do nothing and return a null pointer
      ! If the damage_recovery_scalar is zero, which
      ! would be an unusual testing case, but possible,
      ! then no recovery is possible, do nothing and
      ! return a null pointer
      if ((ccohort%crowndamage == undamaged_class) .or. &
           (EDPftvarcon_inst%damage_recovery_scalar(ipft) < nearzero) ) then
         newly_recovered = .false.
         return
      end if

      
      ! If we have not returned, then this cohort both has
      ! a damaged status, and the ability to recover from that damage
      ! -----------------------------------------------------------------

      ! To determine recovery, the first priority is to determine how much
      ! resources (C,N,P) are required to recover the plant to the target
      ! pool sizes of the next (less) damage class
      
      ! Target sapwood biomass according to allometry and trimming [kgC]
      call bsap_allom(dbh,ipft, ccohort%crowndamage-1, canopy_trim,sapw_area,target_sapw_c)
      ! Target total above ground biomass in woody/fibrous tissues  [kgC]
      call bagw_allom(dbh,ipft, ccohort%crowndamage-1, target_agw_c)
      ! Target total below ground biomass in woody/fibrous tissues [kgC] 
      call bbgw_allom(dbh,ipft,target_bgw_c)
      ! Target total dead (structrual) biomass [kgC]
      call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)
      ! Target fine-root biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
      call bfineroot(dbh,ipft,canopy_trim,target_fnrt_c)
      ! Target storage carbon [kgC,kgC/cm]
      call bstore_allom(dbh,ipft,ccohort%crowndamage-1, canopy_trim,target_store_c)
      ! Target leaf biomass according to allometry and trimming
      if(ccohort%status_coh==leaves_on) then
         call bleaf(dbh,ipft,ccohort%crowndamage-1, canopy_trim,target_leaf_c)
      else
         target_leaf_c = 0._r8
      end if

      ! We will be taking the number of recovering plants
      ! based on minimum of available resources for C/N/P (initialize high)
      nplant_recover = 1.e10_r8
      
      do el=1,num_elements
         
         ! Actual mass of chemical species in the organs
         leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
         store_m  = ccohort%prt%GetState(store_organ, element_list(el))
         sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
         fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
         struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
         repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
         
         ! Target mass of chemical species in organs, based on stature,
         ! allometry and stoichiometry parameters
         select case (element_list(el))
         case (carbon12_element)
            target_store_m  = target_store_c
            target_leaf_m   = target_leaf_c
            target_fnrt_m   = target_fnrt_c
            target_struct_m = target_struct_c
            target_sapw_m   = target_sapw_c
            target_repro_m  = 0._r8
            available_m     = ccohort%npp_acc
         case (nitrogen_element) 
            target_struct_m = target_struct_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
            target_leaf_m = target_leaf_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
            target_fnrt_m = target_fnrt_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
            target_sapw_m = target_sapw_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
            target_repro_m  = 0._r8
            target_store_m = StorageNutrientTarget(ipft, element_list(el), &
                 target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            ! For nutrients, all uptake is immediately put into storage, so just swap
            ! them and assume storage is what is available, but needs to be filled up
            available_m     = store_m
            store_m         = 0._r8
         case (phosphorus_element)
            target_struct_m = target_struct_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
            target_leaf_m = target_leaf_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
            target_fnrt_m = target_fnrt_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
            target_sapw_m = target_sapw_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
            target_repro_m  = 0._r8
            target_store_m = StorageNutrientTarget(ipft, element_list(el), &
                 target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            ! For nutrients, all uptake is immediately put into storage, so just swap
            ! them and assume storage is what is available, but needs to be filled up
            available_m     = store_m
            store_m         = 0._r8
         end select
         
         ! 1. What is excess carbon?
         ! carbon_balance
         
         !  2. What is biomass required to go from current
         !     damage level to next damage level?
         
         ! mass of this damage class
         mass_d = leaf_m + store_m + sapw_m + fnrt_m + struct_m + repro_m
         
         mass_dminus1 = max(leaf_m, target_leaf_m) + max(fnrt_m, target_fnrt_m) + &
              max(store_m, target_store_m) + max(sapw_m, target_sapw_m) + &
              max(struct_m, target_struct_m)
         
         ! Mass needed to get from current mass to allometric
         ! target mass of next damage class up
         recovery_demand = mass_dminus1 - mass_d
         
         ! 3. How many trees can get there with excess carbon?
         max_recover_nplant =  available_m * ccohort%n / recovery_demand 
         
         ! 4. Use the scalar to decide how many to recover
         nplant_recover = min(nplant_recover,max(0._r8,max_recover_nplant * &
                              EDPftvarcon_inst%damage_recovery_scalar(ipft) ))
         
      end do
          
      ! there is a special case where damage_recovery_scalar = 1, but
      ! max_recover_nplant > n (i.e. there is more carbon than needed for all
      ! individuals to recover to the next damage class.
      ! in this case we can cheat, by making n_recover 0 and simply
      ! allowing the donor cohort to recover and then go through
      ! prt - will this work though? if they are not anywhere near allometry?
      
      if( abs(EDPftvarcon_inst%damage_recovery_scalar(ipft)-1._r8) < nearzero .and. &
           nplant_recover > ccohort%n) then
         nplant_recover = 0.0_r8
         ccohort%crowndamage = ccohort%crowndamage - 1
      end if

      if(nplant_recover < nearzero) then

         newly_recovered = .false.
         return
         
      else
         newly_recovered = .true.
         allocate(rcohort)
         if(hlm_use_planthydro .eq. itrue) call InitHydrCohort(csite,rcohort)
         ! Initialize the PARTEH object and point to the
         ! correct boundary condition fields
         rcohort%prt => null()
         call InitPRTObject(rcohort%prt)
         call InitPRTBoundaryConditions(rcohort)
         call copy_cohort(ccohort, rcohort)

         rcohort%n = nplant_recover
          
         rcohort%crowndamage = ccohort%crowndamage - 1
         
         ! Need to adjust the crown area which is NOT on a per individual basis
         call carea_allom(dbh,rcohort%n,csite%spread,ipft,rcohort%crowndamage,rcohort%c_area)
         !rcohort%n/n_old * ccohort%c_area
         !ccohort%c_area = ccohort%c_area - rcohort%c_area

         ! Update properties of the un-recovered (donor) cohort
         ccohort%n = ccohort%n - rcohort%n
         ccohort%c_area = ccohort%c_area * ccohort%n / (ccohort%n+rcohort%n)

         !----------- Insert copy into linked list ----------------------!
         ! This subroutine is called within a loop in EDMain that
         ! proceeds short to tall. We want the newly created cohort
         ! to have an opportunity to experience the list, so we add
         ! it in the list in a position taller than the current cohort
         ! --------------------------------------------------------------!
         
         rcohort%shorter => ccohort
         if(associated(ccohort%taller))then
            rcohort%taller => ccohort%taller
            ccohort%taller%shorter => rcohort
         else
            cpatch%tallest => rcohort    
            rcohort%taller => null()
         endif
         ccohort%taller => rcohort
         
      end if ! end if greater than nearzero

    end associate
    
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
       damage_date  = icode - int(100* floor(real(icode,r8)/100._r8))
       damage_year  = floor(real(icode,r8)/10000._r8)
       damage_month = floor(real(icode,r8)/100._r8) - damage_year*100

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
    crown_reduction = min(1.0_r8, (real(crowndamage,r8) - 1.0_r8) * class_width)

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
    
    class_width = 1.0_r8/real(nlevdamage,r8)
    
    ! parameter to determine slope of exponential
    damage_mort_p1 = EDPftvarcon_inst%damage_mort_p1(pft)
    damage_mort_p2 = EDPftvarcon_inst%damage_mort_p2(pft)

    ! make damage mortality a function of crownloss and not crowndamage
    ! class so that it doesn't need to be re-parameterised if the number
    ! of damage classes change.
    crown_loss = min(1.0_r8, (real(crowndamage,r8) - 1.0_r8) * class_width)
    
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

