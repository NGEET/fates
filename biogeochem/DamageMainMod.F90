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

  logical, protected :: DamageTime  ! if true then damage occurs during current time step

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
  public :: GetCrownReduction
  public :: GetDamageFrac
  public :: IsItDamageTime
  public :: DamageTime
  public :: GetDamageMortality
  
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains


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

    DamageTime = .false.
    icode = int(damage_event_code)

    model_day_int = nint(hlm_model_day)
    
    if(icode .eq. 1) then
       ! Damage is turned off 
       DamageTime = .false.

    else if(icode .eq. 2) then
       ! Damage event on first time step 
       if(model_day_int .eq.1) then
          DamageTime = .true.
       end if

    else if(icode .eq. 3) then
       ! Damage event every day - not sure this is recommended as it will result in a very large
       ! number of cohorts 
       DamageTime = .true.

    else if(icode .eq. 4) then
       ! Damage event once a month
       if(hlm_current_day.eq.1 ) then
          DamageTime = .true.
       end if

    else if(icode < 0 .and. icode > -366) then
       ! Damage event every year on a specific day of the year
       if(hlm_day_of_year .eq. abs(icode) ) then
          DamageTime = .true.
       end if

    else if(icode > 10000 ) then
       ! Specific Event: YYYYMMDD
       damage_date  = icode - int(100* floor(real(icode)/100))
       damage_year  = floor(real(icode)/10000)
       damage_month = floor(real(icode)/100) - damage_year*100

       if(hlm_current_day .eq. damage_date .and. &
            hlm_current_month .eq. damage_month .and. &
            hlm_current_year .eq. damage_year ) then
          DamageTime = .true.
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

    if(DamageTime .and. (is_master.eq.itrue) ) then
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

