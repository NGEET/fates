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
  use EDParamsMod           , only : ED_val_history_damage_bin_edges
  use EDTypesMod            , only : ed_site_type
  use EDTypesMod            , only : ed_patch_type
  use EDTypesMod            , only : ed_cohort_type
  use EDTypesMod            , only : AREA
  use FatesInterfaceTypesMod, only : hlm_current_day
  use FatesInterfaceTypesMod, only : hlm_current_month
  use FatesInterfaceTypesMod, only : hlm_current_year
  use FatesInterfaceTypesMod, only : hlm_model_day
  use FatesInterfaceTypesMod, only : hlm_day_of_year
  
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
  
  
  logical :: debug = .false.  ! for debugging


  ! The following is the special classification for undamaged plants
  ! and is used in contexts where cohort%damageclass is used. This is
  ! to flag to the user that an undamaged plant is assumed in those contexts
  
  integer, parameter, public :: undamaged_class = 1

  
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

    damage_time = .false.
    icode = int(damage_event_code)

    model_day_int = int(hlm_model_day)
    
    if(icode .eq. 1) then
       ! Damage is turned off 
       damage_time = .false.

    else if(icode .eq. 2) then
       ! Damage event on first time step 
       if(model_day_int .eq.1) then
          damage_time = .true.
       end if

    else if(icode .eq. 3) then
       ! Damage event every day - this is not recommended as it will result in a very large
       ! number of cohorts which will likely be terminated 
       damage_time = .true.

    else if(icode .eq. 4) then
       ! Damage event once a month
       if(hlm_current_day.eq.1 ) then
          damage_time = .true.
       end if

    else if(icode < 0 .and. icode > -366) then
       ! Damage event every year on a specific day of the year
       ! specified as negative day of year
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
       write(fates_log(),*) 'for a breakdown of the valid codes and change'
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


    ! Given the current cohort damage class find the fraction of individuals
    ! going to the new damage class.
    ! This subroutine consults a look up table of transitions from param derived. 
    
    ! USES
    use FatesParameterDerivedMod, only : param_derived

       
    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    integer, intent(in) :: pft                     ! plant functional type
    real(r8), intent(out) :: dist_frac             ! fraction of current cohort moving to
                                                   ! new damage level

    dist_frac = param_derived%damage_transitions(cc_cd, nc_cd, pft) 
    
  end subroutine GetDamageFrac

  !-------------------------------------------------------    
  
  subroutine GetCrownReduction(crowndamage, crown_reduction)

    !------------------------------------------------------------------
    ! This subroutine takes the crown damage class of a cohort (integer)
    ! and returns the fraction of the crown that is lost.                                                                  
    !-------------------------------------------------------------------       

    integer(i4), intent(in)   :: crowndamage        ! crown damage class of the cohort
    real(r8),    intent(out)  :: crown_reduction    ! fraction of crown lost from damage

    crown_reduction = ED_val_history_damage_bin_edges(crowndamage)/100.0_r8
    
    return
  end subroutine GetCrownReduction


  !----------------------------------------------------------------------------------------


  subroutine GetDamageMortality(crowndamage,pft, dgmort)

    !------------------------------------------------------------------
    ! This subroutine calculates damage-dependent mortality. 
    ! Not all damage related mortality will be captured by mechanisms in FATES
    ! (e.g. carbon starvation mortality). Damage could also lead to damage 
    ! due to unrepresented mechanisms such as pathogens or increased 
    ! vulnerability to wind throws. This function captures mortality due to 
    ! those unrepresented mechanisms. 
    !------------------------------------------------------------------

    use EDPftvarcon                , only : EDPftvarcon_inst
    
    integer(i4), intent(in) :: crowndamage       ! crown damage class of the cohort
    integer(i4), intent(in) :: pft               ! plant functional type
    real(r8),    intent(out) :: dgmort           ! mortality directly associated with damage

    ! local variables
    real(r8) :: damage_mort_p1       ! inflection  point of the damage mortalty relationship
    real(r8) :: damage_mort_p2       ! rate parameter for the damage mortality relationship
    real(r8) :: crown_loss           ! fraction of  crown lost 
        
    damage_mort_p1 = EDPftvarcon_inst%damage_mort_p1(pft)
    damage_mort_p2 = EDPftvarcon_inst%damage_mort_p2(pft)

    ! make damage mortality a function of crownloss and not crowndamage
    ! class so that it doesn't need to be re-parameterised if the number
    ! of damage classes change.
    crown_loss = ED_val_history_damage_bin_edges(crowndamage)/100.0_r8

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

