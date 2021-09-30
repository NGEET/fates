module DamageMainMod

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : i4 => fates_int
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : years_per_day
  use FatesGlobals          , only : fates_log

  use EDPftvarcon           , only : EDPftvarcon_inst

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

  implicit none
  private

  logical, protected :: damage_time  ! if true then damage occurs during current time step

  public :: get_crown_reduction
  public :: get_damage_frac
  public :: is_it_damage_time
  public :: damage_time
  public :: get_damage_mortality
  
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains


  subroutine is_it_damage_time(is_master, currentSite)

    !----------------------------------------------------------------------------
    ! This subroutine determines whether damage should occur (it is called daily)
    !-----------------------------------------------------------------------------

    use FatesInterfaceTypesMod , only : hlm_day_of_year

    integer, intent(in) :: is_master
    type(ed_site_type), intent(inout), target :: currentSite
    
    
    damage_time = .false.

    if (hlm_day_of_year .eq. 1) then
       damage_time = .true.
    end if
   
  end subroutine is_it_damage_time
  
  !----------------------------------------------------------------------------

  subroutine get_damage_frac(cc_cd, nc_cd, pft, dist_frac)


    ! given current cohort damage class find the fraction of individuals
    ! going to the new damage class.
    ! Consults a look up table of transitions from param derived. 
    
    ! USES
    use FatesInterfaceTypesMod, only : ncrowndamage
    use FatesConstantsMod, only : years_per_day
    use FatesParameterDerivedMod, only : param_derived

       
    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    integer, intent(in) :: pft
    real(r8), intent(out) :: dist_frac             ! probability of current cohort moving to new damage level

    dist_frac = param_derived%damage_transitions(cc_cd, nc_cd, pft) !* years_per_day (if damage is occuring annually don't do this)
    
    
  end subroutine get_damage_frac

  !-------------------------------------------------------    
  
  subroutine get_crown_reduction(crowndamage, crown_reduction)

    !------------------------------------------------------------------                                                                     
    ! This function takes the crown damage class of a cohort (integer)
    ! and returns the fraction of the crown that is lost
    ! Since crowndamage class = 1 means no damage, we subtract one     
    ! before multiplying by 0.2                                    
    ! Therefore, first damage class is 20% loss of crown, second 40% etc.                                                                   
    !-------------------------------------------------------------------                                                                    
    use FatesInterfaceTypesMod     , only : ncrowndamage

    integer(i4), intent(in)   :: crowndamage
    real(r8),    intent(out)  :: crown_reduction

    ! local variables
    real(r8) :: class_width

    class_width = 1.0_r8/ncrowndamage
    crown_reduction = min(1.0_r8, (real(crowndamage) - 1.0_r8) * class_width)

    return
  end subroutine get_crown_reduction


  !----------------------------------------------------------------------------------------


  subroutine get_damage_mortality(crowndamage,pft, dgmort)

    use FatesInterfaceTypesMod     , only : ncrowndamage
    use EDPftvarcon                , only : EDPftvarcon_inst
    
    integer(i4), intent(in) :: crowndamage
    integer(i4), intent(in) :: pft
    real(r8),    intent(out) :: dgmort

    ! local variables
    real(r8) :: damage_mort_p1
    real(r8) :: damage_mort_p2
    real(r8) :: class_width
    real(r8) :: crown_loss
    
    class_width = 1.0_r8/real(ncrowndamage)
    
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
  end subroutine get_damage_mortality
  !----------------------------------------------------------------------------------------

 
end module DamageMainMod

