module FatesTestLeafPhotoMod
  !
  ! DESCRIPTION:
  ! Helper methods for running leaf-level photosynthesis
  !
  
  use FatesConstantsMod,      only : r8 => fates_r8
  use FatesInterfaceTypesMod, only : hlm_maintresp_leaf_model
  use FatesConstantsMod,      only : lmrmodel_ryan_1991
  use FatesConstantsMod,      only : lmrmodel_atkin_etal_2017
  use PRTParametersMod,       only : prt_params
  use PRTGenericMod,          only : leaf_organ
  use LeafBiophysicsMod,      only : GetCanopyGasParameters
  use LeafBiophysicsMod,      only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod,      only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,      only : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysicsMod,      only : LeafLayerMaintenanceRespiration_Atkin_etal_2017
  
  implicit none
  private
  
  ! convergence tolerance for intracellular CO2 [Pa]
  real(r8), public, parameter :: ci_tol = 0.5_r8

  type, public :: leaf_capacity_type
    ! One leaf layer's photosynthetic capacity and dark respiration, evaluated at a
    ! specific set of leaf conditions
    real(r8) :: vcmax ! maximum rate of carboxylation [umolCO2/m2 leaf/s]
    real(r8) :: jmax  ! maximum electron transport rate [umol electrons/m2 leaf/s]
    real(r8) :: kp    ! initial slope of the CO2 response curve, C4 plants [umol/m2 leaf/s]
    real(r8) :: gs0   ! effective stomatal conductance intercept [umol H2O/m2 leaf/s]
    real(r8) :: gs1   ! effective stomatal conductance slope
    real(r8) :: gs2   ! alternative btran term applied to the whole non-intercept side of the Medlyn conductance equation
    real(r8) :: lmr   ! leaf maintenance (dark) respiration [umolCO2/m2 leaf/s]
  end type leaf_capacity_type

  public :: LeafNitrogenContent
  public :: EvaluateLeafPhotosynthesis
  
  contains 
  
  ! ==========================================================================
  
  function LeafNitrogenContent(pft) result(lnc_top)
    !
    ! DESCRIPTION:
    ! Leaf nitrogen content at the canopy top [gN/m2 leaf], needed by both leaf
    ! maintenance respiration models. Requires PRTDerivedParams() to have
    ! already populated prt_params%organ_param_id

    ! ARGUMENTS:
    integer, intent(in) :: pft ! plant functional type index

    ! LOCALS:
    real(r8) :: lnc_top ! leaf N content at the canopy top [gN/m2 leaf]

    lnc_top = prt_params%nitr_stoich_p1(pft, prt_params%organ_param_id(leaf_organ)) / &
      prt_params%slatop(pft)

  end function LeafNitrogenContent
  
  ! ==========================================================================
    
  subroutine EvaluateLeafPhotosynthesis(pft, par_abs, veg_tempk, t_growth, t_home, &
    can_press, can_co2_ppress, can_o2_ppress, veg_esat, can_vpress, gb, nscaler,   &
    rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top,     &
    agross, anet, gs, ci)
    !
    ! DESCRIPTION:
    ! Evaluates leaf-level photosynthesis at arbitrary prescribed driver
    ! conditions. This reproduces the full current production call sequence:
    !
    !   GetCanopyGasParameters -> LeafLayerBiophysicalRates ->
    !   LeafLayerMaintenanceRespiration_* -> LeafLayerPhotosynthesis
    !
    ! lb_params' model switches (electron_transport_model, stomatal_model,
    ! stomatal_assim_model, photo_tempsens_model) and hlm_maintresp_leaf_model
    ! are HLM-namelist-controlled in production and so are not set by any call
    ! in this module. The calling driver must set them explicitly
  
    ! ARGUMENTS:
    integer,  intent(in)  :: pft            ! plant functional type index
    real(r8), intent(in)  :: par_abs        ! absorbed PAR per unit leaf area [umol photons/m2 leaf/s]
    real(r8), intent(in)  :: veg_tempk      ! instantaneous leaf temperature [K]
    real(r8), intent(in)  :: t_growth       ! 10-day running-mean growth temperature [K]
    real(r8), intent(in)  :: t_home         ! long-term running-mean home temperature [K]
    real(r8), intent(in)  :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8), intent(in)  :: veg_esat       ! saturation vapor pressure at veg_tempk [Pa]
    real(r8), intent(in)  :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8), intent(in)  :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8), intent(in)  :: nscaler        ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8), intent(in)  :: rdark_scaler   ! leaf respiration vertical-scaling factor [0-1], Atkin only
    real(r8), intent(in)  :: dayl_factor    ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8), intent(in)  :: btran          ! soil moisture stress factor [0-1]
    real(r8), intent(in)  :: vcmax25top     ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8), intent(in)  :: jmax25top      ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8), intent(in)  :: kp25top        ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8), intent(in)  :: lnc_top        ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    real(r8), intent(out) :: agross         ! gross photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: anet           ! net photosynthesis [umolC/m2/s]
    real(r8), intent(out) :: gs             ! leaf stomatal conductance [umol H2O/m2/s]
    real(r8), intent(out) :: ci             ! intracellular CO2 [Pa]

    ! LOCALS:
    type(leaf_capacity_type) :: cap         ! leaf biophysical capacity/dark respiration at these conditions
    real(r8)                 :: mm_kco2     ! Michaelis-Menten constant for CO2 [Pa]
    real(r8)                 :: mm_ko2      ! Michaelis-Menten constant for O2 [Pa]
    real(r8)                 :: co2_cpoint  ! Michaelis-Menten constants for CO2 compensation point [Pa]
    real(r8)                 :: c13disc     ! carbon-13 discrimination (unused diagnostic here)
    integer                  :: solve_iter  ! Ci-solver iteration count (unused diagnostic here)

    call GetCanopyGasParameters(can_press, can_o2_ppress, veg_tempk, mm_kco2, mm_ko2, co2_cpoint)

    call LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler,             &
      rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top,  &
      cap)

    call LeafLayerPhotosynthesis(par_abs, pft, cap%vcmax, cap%jmax, cap%kp,       &
      cap%gs0, cap%gs1, cap%gs2, veg_tempk, can_press, can_co2_ppress,            &
      can_o2_ppress, veg_esat, gb, can_vpress, mm_kco2, mm_ko2, co2_cpoint,       &
      cap%lmr, ci_tol, agross, gs, anet, c13disc, ci, solve_iter)
      
  end subroutine EvaluateLeafPhotosynthesis
  
  ! ==========================================================================
  
  subroutine LeafLayerCapacity(pft, veg_tempk, t_growth, t_home, nscaler,         &
    rdark_scaler, dayl_factor, btran, vcmax25top, jmax25top, kp25top, lnc_top, cap)
    !
    ! DESCRIPTION:
    ! One leaf layer's photosynthetic capacity and dark respiration at the given
    ! leaf conditions
    !
    ! The leaf maintenance respiration model is selected by
    ! hlm_maintresp_leaf_model, which is HLM-namelist-controlled in production
    ! and so must be set explicitly by the calling driver
    !
    ! t_growth doubles as the Atkin et al. (2017) acclimation temperature, the
    ! same quantity production passes there (currentPatch%tveg_lpa%GetMean())
    !

    ! ARGUMENTS:
    integer,                  intent(in)  :: pft          ! plant functional type index
    real(r8),                 intent(in)  :: veg_tempk    ! instantaneous leaf temperature [K]
    real(r8),                 intent(in)  :: t_growth     ! 10-day running-mean growth temperature [K]
    real(r8),                 intent(in)  :: t_home       ! long-term running-mean home temperature [K]
    real(r8),                 intent(in)  :: nscaler      ! leaf nitrogen vertical-scaling factor [0-1]
    real(r8),                 intent(in)  :: rdark_scaler ! leaf respiration vertical-scaling factor [0-1], Atkin only
    real(r8),                 intent(in)  :: dayl_factor  ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8),                 intent(in)  :: btran        ! soil moisture stress factor [0-1]
    real(r8),                 intent(in)  :: vcmax25top   ! reference (25C, canopy-top) maximum carboxylation rate [umol/m2/s]
    real(r8),                 intent(in)  :: jmax25top    ! reference (25C, canopy-top) maximum electron transport rate [umol/m2/s]
    real(r8),                 intent(in)  :: kp25top      ! reference (25C, canopy-top) initial slope of C4 CO2 response [umol/m2/s]
    real(r8),                 intent(in)  :: lnc_top      ! leaf N content at the canopy top [gN/m2 leaf] (see LeafNitrogenContent)
    type(leaf_capacity_type), intent(out) :: cap          ! this layer's capacity/dark respiration at these conditions

    call LeafLayerBiophysicalRates(pft, vcmax25top, jmax25top, kp25top, nscaler,   &
      veg_tempk, dayl_factor, t_growth, t_home, btran, cap%vcmax, cap%jmax,       &
      cap%kp, cap%gs0, cap%gs1, cap%gs2)

    select case (hlm_maintresp_leaf_model)

      case (lmrmodel_ryan_1991)

        call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, nscaler, pft,      &
          veg_tempk, cap%lmr)

      case (lmrmodel_atkin_etal_2017)

        call LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top,              &
          rdark_scaler, pft, veg_tempk, t_growth, cap%lmr)

      case default

        write(*,*) 'LeafLayerCapacity: unrecognized leaf respiration model: ',     &
          hlm_maintresp_leaf_model
        error stop

    end select

  end subroutine LeafLayerCapacity
  
end module FatesTestLeafPhotoMod