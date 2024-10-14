module LeafBiophysicsMod

  !-------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! 
  ! This module contains routines for leaf-level biophyics.  These
  ! routines are all called with primitive arguments to facilitate
  ! use accross models, with the exception of an internally defined
  ! set of constants associated with plant functional type.
  !
  ! Assumptions:
  !
  ! Units are always in micro-moles (umol), square meters (m2) and seconds
  !
  ! Calculations are per-square-meter of leaf tissue.
  !
  ! 
  !
  ! ------------------------------------------------------------------------------------

  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use shr_sys_mod,       only : shr_sys_abort
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  use FatesGlobals,      only : FatesWarn,N2S,A2S,I2S
  use FatesConstantsMod, only : itrue
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : molar_mass_ratio_vapdry
  use FatesConstantsMod, only : molar_mass_water
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
  use FatesConstantsMod, only : wm2_to_umolm2s
  use FatesConstantsMod, only : nocomp_bareground
  use FatesConstantsMod, only : lmrmodel_ryan_1991
  use FatesConstantsMod, only : lmrmodel_atkin_etal_2017
  use FatesConstantsMod, only : kpa_per_pa
  use FatesConstantsMod, only : umol_per_kmol
  use FatesUtilsMod,     only : QuadraticRoots => QuadraticRootsSridharachary
  use FatesConstantsMod, only : rgas_J_K_kmol
  use FatesConstantsMod, only : rgas_J_K_mol
  use FatesConstantsMod, only : g_per_kg
  
  implicit none
  private

  public :: LeafLayerPhotosynthesis
  public :: LeafHumidityStomaResis
  public :: GetCanopyGasParameters
  public :: GetStomatalInterceptBtran
  public :: LeafLayerMaintenanceRespiration_Ryan_1991
  public :: LeafLayerMaintenanceRespiration_Atkin_etal_2017
  public :: LeafLayerBiophysicalRates
  public :: LowstorageMainRespReduction
  public :: GetConstrainedVPress
  public :: DecayCoeffVcmax
  public :: StomatalVaporPressureFromLWP
  public :: QSat
  public :: AgrossRubiscoC3
  public :: AgrossRuBPC3
  public :: AgrossRuBPC4
  public :: AgrossPEPC4
  public :: StomatalCondMedlyn
  public :: StomatalCondBallBerry
  public :: VeloToMolarCF
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  character(len=1024) :: warn_msg   ! for defining a warning message

  !-------------------------------------------------------------------------------------

  ! maximum stomatal resistance [s/m] (used across several procedures)
  real(r8),public, parameter :: rsmax0 =  2.e8_r8

  ! minimum allowable stomatal conductance for numerics purposes [m/s]
  real(r8),parameter :: gsmin0 = 1._r8/rsmax0

  ! minimum allowable conductance for numerics purposes at 20C (293.15K)
  ! and 1 standard atmosphere (101325 Pa) [umol/m2/s]
  ! this follows the same math as  FatesPlantRespPhotosynthMod:FetMolarVeloCF()
  real(r8),parameter :: gsmin0_20C1A_mol = gsmin0 * 101325.0_r8/(rgas_J_K_kmol * 293.15 )*umol_per_kmol

  ! Set this to true to perform debugging
  logical,parameter   ::  debug = .true.

  ! Set this to true if you want to assume zero resistance between
  ! leaf surface and boundary layer, which then assumes leaf surface
  ! humidity is the same as humidity on the other side of the boundary
  ! layer, essentially bypassing the use
  ! of Fick's law while calculating stomatal conductance.
  logical, parameter :: zero_bl_resist = .false.

  ! Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
  real(r8),parameter :: h2o_co2_stoma_diffuse_ratio = 1.6_r8

  ! Ratio of H2O/CO2 gass diffusion in the leaf boundary layer (approximate)
  real(r8),parameter :: h2o_co2_bl_diffuse_ratio = 1.4_r8

  ! Constants used to define C3 versus C4 photosynth pathways
  integer, parameter :: c3_path_index = 1
  integer, parameter :: c4_path_index = 0


  ! Constants used to define conductance models
  integer, parameter :: medlyn_model = 2
  integer, parameter :: ballberry_model = 1

  ! Alternatively, Gross Assimilation can be used to estimate
  ! leaf co2 partial pressure and therefore conductance. The default
  ! is to use anet
  integer, parameter :: net_assim_model = 1
  integer, parameter :: gross_assim_model = 2

  ! Constants defining the photosynthesis temperature acclimation model
  integer, parameter :: photosynth_acclim_model_none = 1
  integer, parameter :: photosynth_acclim_model_kumarathunge_etal_2019 = 2


  ! Rdark constants from Atkin et al., 2017 https://doi.org/10.1007/978-3-319-68703-2_6
  ! and Heskel et al., 2016 https://doi.org/10.1073/pnas.1520282113
  real(r8), parameter :: lmr_b = 0.1012_r8       ! (degrees C**-1)
  real(r8), parameter :: lmr_c = -0.0005_r8      ! (degrees C**-2)
  real(r8), parameter :: lmr_TrefC = 25._r8      ! (degrees C)
  ! These two are public for error checking during parameter read-in
  real(r8), parameter, public :: lmr_r_1 = 0.2061_r8     ! (umol CO2/m**2/s / (gN/(m2 leaf))) 
  real(r8), parameter, public :: lmr_r_2 = -0.0402_r8    ! (umol CO2/m**2/s/degree C)

  ! Fraction of light absorbed by non-photosynthetic pigments
  real(r8),parameter :: fnps = 0.15_r8
  
  ! term accounting that two photons are needed to fully transport a single 
  ! electron in photosystem 2
  real(r8), parameter :: photon_to_e = 0.5_r8
  
  ! empirical curvature parameter for electron transport rate
  real(r8),parameter :: theta_psii = 0.7_r8
  
  ! These are parameter constants read in externally
  ! some are differentiated by PFT, others are not
  ! -------------------------------------------------------------------------------------

  type, public :: leafbiophys_params_type

     integer, allocatable :: c3psn(:)                             ! Photosynthetic pathway index (C3=1,C4=0,etc)
     real(r8),allocatable :: medlyn_slope(:)                      ! Stomatal Slope, Medlyn, e.g. g1 [-]
     real(r8),allocatable :: bb_slope(:)                          ! Stomatal Slope, Ball-Berry, e.g. g1 [-]
     real(r8),allocatable :: stomatal_intercept(:)                ! Stomatal int, BB or Medlyn, e.g. g0, [-]
     real(r8),allocatable :: maintresp_leaf_ryan1991_baserate(:)  ! Base maintenance resp rate M.Ryan 1991 [gC gN-1 s-1]
     real(r8),allocatable :: maintresp_leaf_atkin2017_baserate(:) ! Base maintenance resp rate Atkin 2017 [umol CO2 m-2 s-1]
     real(r8),allocatable :: maintresp_reduction_curvature(:)     ! curvature of MR reduction as f(carbon storage),
                                                                  ! 1=linear, 0=very curved
     real(r8),allocatable :: maintresp_reduction_intercept(:)     ! intercept of MR reduction as f(carbon storage),
                                                                  ! 0=no throttling, 1=max throttling
     real(r8),allocatable :: maintresp_reduction_upthresh (:)     ! Upper threshold for storage biomass (relative 
                                                                  ! to leaf biomass) above which MR is not reduced
     real(r8),allocatable :: vcmaxha(:)                           ! activation energy for vcmax (J/mol)
     real(r8),allocatable :: jmaxha(:)                            ! activation energy for jmax (J/mol)
     real(r8),allocatable :: vcmaxhd(:)                           ! deactivation energy for vcmax (J/mol)
     real(r8),allocatable :: jmaxhd(:)                            ! deactivation energy for jmax (J/mol)
     real(r8),allocatable :: vcmaxse(:)                           ! entropy term for vcmax (J/mol/K)
     real(r8),allocatable :: jmaxse(:)                            ! entropy term for jmax (J/mol/K)
     real(r8)             :: theta_cj_c3                          ! Empirical curvature parameter for ac and 
                                                                  ! aj photosynthesis co-limitation in c3 plants
     real(r8)             :: theta_cj_c4                          ! Empirical curvature parameter for ac and aj
                                                                  ! photosynthesis co-limitation in c4 plants
     integer              :: dayl_switch                          ! switch for turning on or off day length factor
                                                                  ! scaling for photosynthetic parameters
     integer              :: photo_tempsens_model                 ! switch for choosing the model that defines the temperature
                                                                  ! sensitivity of photosynthetic parameters (vcmax, jmax).
                                                                  ! 1=non-acclimating, 2=Kumarathunge et al., 2019
     integer              :: stomatal_assim_model                 ! Switch designating whether to use net or
                                                                  ! gross assimilation in the stomata model
     integer              :: stomatal_model                       ! switch for choosing between stomatal conductance models,
                                                                  ! for Ball-Berry, 2 for Medlyn
     
     ! -------------------------------------------------------------------------------------
     ! Note the omission of several parameter constants:
     !
     ! Vcmax at 25C, canopy top
     ! Jmax at 25C, at canopy
     ! Kp (initial co2 response slope) at 25C, canopy top
     !
     ! These parameters are omitted because some models (like FATES) have these as functions
     ! of leaf age. By setting these as arguments in our routines, we add a degree of
     ! flexibility. So these routines can be run on mixed-age leaves if desired.
     !--------------------------------------------------------------------------------------

  end type leafbiophys_params_type


  type(leafbiophys_params_type),public :: lb_params



  ! A possible sequence of calls for leaf biophysics is as follows:
  ! 1) determine any gas quantities or parameters that are derived and
  !    are applicable to a super-set of leaf-layers (like MM, and compensation points)
  ! 2) loop over discrete portions of leaves, perhaps differentiated by spatial position and pft
  ! 3) determine if this particular leaf is present and has been solved
  ! 4) solve for leaf maintenance respiration (e.g. dark respiration) 
  ! 5) update leaf-level rates, such as vcmax, jmax
  ! 6) solve for photosynthesis
  ! 7) solve for maintenance respiration of other tissues (other library?)

contains
  
  subroutine StomatalCondMedlyn(anet,ft,veg_esat,can_vpress,stomatal_intercept_btran,leaf_co2_ppress,can_press,gb,gs)

    ! Input
    real(r8), intent(in) :: anet                     ! net leaf photosynthesis (umol CO2/m**2/s)
    integer, intent(in)  :: ft                       ! plant functional type index
    real(r8), intent(in) :: veg_esat                 ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press                ! Air pressure NEAR the surface of the leaf (Pa)
    real(r8), intent(in) :: gb                       ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress               ! vapor pressure of canopy air (Pa)
    real(r8), intent(in) :: stomatal_intercept_btran ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: leaf_co2_ppress          ! CO2 partial pressure at the leaf surface (Pa)
                                                     ! at the boundary layer interface with the leaf

    ! Output
    real(r8),intent(out) :: gs                       ! leaf stomatal conductance (umol H2O/m**2/s)
    

    ! locals
    real(r8) :: vpd                                  ! water vapor deficit in Medlyn stomatal model [KPa]
    real(r8) :: term                                 ! intermediate term used to simplify equations
    real(r8) :: aquad,bquad,cquad                    ! quadradic solve terms
    real(r8) :: r1,r2                                ! quadradic solve roots

    ! Evaluate trival solution, if there is no positive net assimiolation
    ! the stomatal conductance is the intercept conductance
    if (anet <= nearzero) then
       gs = stomatal_intercept_btran
       return
    end if

    
   
    
    ! stomatal conductance calculated from Medlyn et al. (2011), the numerical &
    ! implementation was adapted from the equations in CLM5.0 [kPa]

    vpd =  max((veg_esat - can_vpress), 50._r8) * kpa_per_pa       !addapted from CLM5. Put some constraint on VPD

    if(zero_bl_resist) then

       ! We assume zero resistance in the leaf boundary layer, and that humidity at
       ! the leaf surface is equal to humidity outside the boundary layer
       gs = h2o_co2_stoma_diffuse_ratio*(1._r8 + lb_params%medlyn_slope(ft)/sqrt(vpd)) * &
            anet/leaf_co2_ppress*can_press + stomatal_intercept_btran

    else

       term = h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
       aquad = 1.0_r8
       bquad = -(2.0 * (stomatal_intercept_btran+ term) + (lb_params%medlyn_slope(ft) * term)**2 / &
            (gb * vpd ))
       cquad = stomatal_intercept_btran*stomatal_intercept_btran + &
            (2.0*stomatal_intercept_btran + term * &
            (1.0 - lb_params%medlyn_slope(ft)* lb_params%medlyn_slope(ft) / vpd)) * term
       
       call QuadraticRoots(aquad, bquad, cquad, r1, r2)
       gs = max(r1,r2)
       
    end if
    
    return
  end subroutine StomatalCondMedlyn

  ! =======================================================================================

  subroutine StomatalCondBallBerry(a_gs,ft,veg_esat,can_vpress,stomatal_intercept_btran,leaf_co2_ppress,can_press, gb, gs)


                                                     ! Input
    integer, intent(in)  :: ft                       ! plant functional type index
    real(r8), intent(in) :: veg_esat                 ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press                ! Air pressure NEAR the surface of the leaf (Pa)
    real(r8), intent(in) :: gb                       ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress               ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in) :: stomatal_intercept_btran ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: leaf_co2_ppress          ! CO2 partial pressure at leaf surface (Pa)
    real(r8), intent(in) :: a_gs                     ! The assimilation (a) for calculating conductance (gs)
                                                     ! is either = to anet or agross

                                                     ! Output
    real(r8), intent(out) :: gs                      ! leaf stomatal conductance (umol H2O/m**2/s)

                                                     ! locals
    real(r8) :: hs                                   ! vapor pressure over saturation vapor pressure ratio
    real(r8) :: ceair                                ! constrained vapor pressure (Pa)
    real(r8) :: aquad,bquad,cquad                    ! quadradic solve terms
    real(r8) :: r1,r2                                ! quadradic solve roots

    if (a_gs <= nearzero) then
       gs = stomatal_intercept_btran
       return
    end if
    
    ! Apply a constraint to the vapor pressure
    ceair = GetConstrainedVPress(can_vpress,veg_esat)
    
    if(zero_bl_resist) then

       hs = (ceair/ veg_esat)  
       gs = lb_params%bb_slope(ft)*a_gs*hs/leaf_co2_ppress*can_press + stomatal_intercept_btran

    else
       aquad = leaf_co2_ppress
       bquad = leaf_co2_ppress*(gb - stomatal_intercept_btran) - lb_params%bb_slope(ft) * a_gs * can_press
       cquad = -gb*(leaf_co2_ppress * stomatal_intercept_btran + &
            lb_params%bb_slope(ft) * a_gs * can_press * ceair/ veg_esat )

       call QuadraticRoots(aquad, bquad, cquad, r1, r2)
       gs = max(r1,r2)
    end if


    !print*,gs

    return
  end subroutine StomatalCondBallBerry

  ! =====================================================================================
  
  function AgrossRubiscoC3(vcmax,co2_inter_c,can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2) result(ac)

    ! Input
    real(r8) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8) :: co2_inter_c       ! intercellular leaf CO2 (Pa)
    real(r8) :: co2_cpoint        ! CO2 compensation point (Pa)
    real(r8) :: can_o2_ppress     ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8) :: mm_kco2           ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2            ! Michaelis-Menten constant for O2 (Pa)
    
    ! Output
    real(r8) :: ac               ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    
    ac = vcmax * max(co2_inter_c-co2_cpoint, 0._r8) / &
         (co2_inter_c+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))
    
  end function AgrossRubiscoC3
  
  ! =====================================================================================

  function AgrossRuBPC3(par_abs_wm2,jmax,co2_inter_c,co2_cpoint ) result(aj)

    ! Input
    real(r8) :: par_abs_wm2           ! Absorbed PAR per leaf area [W/m2]
    real(r8) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: co2_inter_c      ! intercellular leaf CO2 (Pa)
    real(r8) :: co2_cpoint        ! CO2 compensation point (Pa)

    ! Output
    real(r8) :: aj               ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)

    ! locals
    real(r8) :: par_abs_umol      ! Absorbed PAR that gets to the photocenters,
                                  ! converted to umol photons/m2/s
    real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    
    ! Electron transport rate for C3 plants.
    ! Convert absorbed photon density [umol/m2 leaf /s] to
    ! that absorbed only by the photocenters (fnps) and also
    ! convert from photon energy into electron transport rate (photon_to_e)
    
    par_abs_umol = par_abs_wm2*photon_to_e*(1.0_r8 - fnps)
    
    ! convert the absorbed par into absorbed par per m2 of leaf,
    ! so it is consistant with the vcmax and lmr numbers.
    aquad = theta_psii
    bquad = -(par_abs_umol + jmax)
    cquad = par_abs_umol * jmax
    call QuadraticRoots(aquad, bquad, cquad, r1, r2)
    je = min(r1,r2)
    
    aj = je * max(co2_inter_c-co2_cpoint, 0._r8) / &
                 (4._r8*co2_inter_c+8._r8*co2_cpoint)

    
  end function AgrossRuBPC3
  
  ! =======================================================================================

  function AgrossRuBPC4(par_abs_wm2) result(aj)

    real(r8) :: par_abs_wm2      ! Absorbed PAR per leaf area [W/m2]
    real(r8) :: aj               ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)

    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(r8),parameter :: c4_quant_eff = 0.05_r8
    
    aj = c4_quant_eff*par_abs_wm2*photon_to_e*(1.0_r8 - fnps)
    
  end function AgrossRuBPC4

  ! =======================================================================================

  function AgrossPEPC4(co2_inter_c,co2_rcurve_islope,can_press) result(ap)

    real(r8) :: co2_inter_c       ! intercellular leaf CO2 (Pa)
    real(r8) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
    real(r8) :: can_press         ! Air pressure near the surface of the leaf (Pa)
    real(r8) :: ap                ! product-limited (C3) or CO2-limited
                                  ! (C4) gross photosynthesis (umol CO2/m**2/s)

    ap = co2_rcurve_islope * max(co2_inter_c, 0._r8) / can_press
    
  end function AgrossPEPC4

  ! =======================================================================================
  
  subroutine LeafLayerPhotosynthesis(par_abs,           &  ! in
       leaf_area,         &  ! in
       ft,                &  ! in
       vcmax,             &  ! in
       jmax,              &  ! in
       co2_rcurve_islope, &  ! in
       veg_tempk,         &  ! in
       can_press,         &  ! in
       can_co2_ppress,    &  ! in
       can_o2_ppress,     &  ! in
       btran,             &  ! in
       gb,                &  ! in
       can_vpress,        &  ! in
       mm_kco2,           &  ! in
       mm_ko2,            &  ! in
       co2_cpoint,        &  ! in
       lmr,               &  ! in
       leaf_psi,          &  ! in   (currently dummy)
       agross_out,        &  ! out
       gs_out,            &  ! out
       anet_out,          &  ! out
       c13disc_out,       &  ! out
       ac,                &  ! out
       aj,                &  ! out
       ap,                &  ! out
       co2_intercel)                   ! out


    ! ------------------------------------------------------------------------------------
    ! This subroutine calculates photosynthesis and stomatal conductance within each leaf
    ! sublayer.
    ! A note on naming conventions: As this subroutine is called for every
    ! leaf-sublayer, many of the arguments are specific to that "leaf sub layer"
    ! (LSL), those variables are given a dimension tag "_lsl"
    ! Other arguments or variables may be indicative of scales broader than the LSL.
    ! ------------------------------------------------------------------------------------

    ! Arguments
    ! ------------------------------------------------------------------------------------
    real(r8), intent(in) :: par_abs           ! Absorbed PAR per leaf area [W/m2]
    real(r8), intent(in) :: leaf_area         ! leaf area per ground area [m2/m2] 
    integer,  intent(in) :: ft                ! (plant) Functional Type Index
    real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in) :: veg_tempk         ! vegetation temperature
    real(r8), intent(in) :: can_press         ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in) :: can_co2_ppress    ! Partial pressure of CO2 near the leaf surface (Pa)
    real(r8), intent(in) :: can_o2_ppress     ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8), intent(in) :: btran             ! transpiration wetness factor (0 to 1)
    real(r8), intent(in) :: gb                ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress        ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in) :: mm_kco2           ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(in) :: mm_ko2            ! Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(in) :: co2_cpoint        ! CO2 compensation point (Pa)
    real(r8), intent(in) :: lmr               ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
    real(r8), intent(in) :: leaf_psi          ! Leaf water potential [MPa]

    real(r8), intent(out) :: agross_out       ! gross photosynthesis (umolC/m2/s)
    real(r8), intent(out) :: gs_out           ! leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out) :: anet_out         ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: c13disc_out      ! carbon 13 in newly assimilated carbon
    real(r8), intent(out) :: ac               ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: aj               ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: ap               ! product-limited (C3) or CO2-limited
                                              ! (C4) gross photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: co2_intercel      ! intercellular leaf CO2 (Pa)

    ! Important Note on the gas pressures as input arguments.  This photosynthesis scheme will iteratively
    ! solve for the co2 partial pressure at the leaf surface (ie in the stomata). The reference
    ! point for these input values are NOT within that boundary layer that separates the stomata from
    ! the canopy air space.  The reference point for these is on the outside of that boundary
    ! layer.  This routine, which operates at the leaf scale, makes no assumptions about what the
    ! scale of the refernce is, it could be lower atmosphere, it could be within the canopy
    ! but most likely it is the closest value one can get to the edge of the leaf's boundary
    ! layer.  We use the convention "can_" because a reference point of within the canopy
    ! ia a best reasonable scenario of where we can get that information from.


    
    ! Locals
    ! ------------------------------------------------------------------------
    
    real(r8) :: a_gs              ! The assimilation (a) for calculating conductance (gs)
                                  ! is either = to anet or agross
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    
    real(r8) :: co2_intercel_old   ! intercellular leaf CO2 (Pa) (previous iteration)
    logical  :: loop_continue     ! Loop control variable
    integer  :: niter             ! iteration loop index
 
    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: stomatal_intercept_btran ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8) :: veg_esat                             ! Saturation vapor pressure at leaf-surface [Pa]
    real(r8) :: veg_qs                               ! DUMMY, specific humidity at leaf-surface [kg/kg]
       
    ! Parameters
    ! ------------------------------------------------------------------------
    

    ! For plants with no leaves, a miniscule amount of conductance
    ! can happen through the stems, at a partial rate of cuticular conductance
    ! THIS IS NOT USED
    !real(r8),parameter :: stem_cuticle_loss_frac = 0.1_r8

   

    ! First guess on ratio between intercellular co2 and the atmosphere
    ! an iterator converges on actual
    real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
    real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

    ! When iteratively solving for intercellular co2 concentration, this
    ! is the maximum tolerable change to accept convergence [Pa]
    real(r8),parameter :: co2_intercel_tol = 1._r8

    ! Maximum number of iterations on intercelluar co2 solver until is quits
    integer, parameter :: max_iters = 10

    ! empirical curvature parameter for ap photosynthesis co-limitation
    real(r8),parameter :: theta_ip = 0.999_r8

    ! Set diagnostics as un-initialized
    ac         = -999._r8
    aj         = -999._r8
    ap         = -999._r8
    
    ! Trivial solution - No biomass - no photosynthesis, no conductance, no respiration, no nothin
    ! --------------------------------------------------------------------------------------------

    if(  leaf_area < nearzero ) then
       agross_out   = 0._r8
       gs_out    = 0._r8 
       anet_out  = 0._r8
       c13disc_out = 0._r8
       return
    end if

    ! Find the stomatal conductance intercept

    stomatal_intercept_btran = max(gsmin0_20C1A_mol,lb_params%stomatal_intercept(ft)*btran)
    
    ! Less, but still trivial solution - biomass, but no light, no photosynthesis
    ! Stomatal conductance is the intercept of the conductance functions
    ! ---------------------------------------------------------------------------------------------
    if (par_abs < nearzero ) then
       anet_out   = -lmr
       agross_out    = 0._r8
       gs_out     = stomatal_intercept_btran
       c13disc_out  = 0.0_r8
       return
    end if


    ! Not trivial solution, some biomass and some light
    ! Initialize first guess of intercellular co2 conc [Pa]
    if (lb_params%c3psn(ft) == c3_path_index) then
       co2_intercel = init_a2l_co2_c3 * can_co2_ppress
    else
       co2_intercel = init_a2l_co2_c4 * can_co2_ppress
    end if

    ! Perform iterative solution to converge on net assimilation,
    ! stomatal conductance and intercellular co2 concentration
    ! this is not a gradient method, it merely re-runs calculations
    ! based off of updated intercellular co2 concentration until
    ! there is minimal deviation from last attempt

    niter = 0
    loop_continue = .true.
    iter_loop: do while(loop_continue)

       ! Increment iteration counter. Stop if too many iterations
       niter = niter + 1

       ! Save old co2_intercel
       co2_intercel_old = co2_intercel

       ! Photosynthesis limitation rate calculations
       if (lb_params%c3psn(ft) == c3_path_index)then

          ! C3: Rubisco-limited photosynthesis
          ac = AgrossRubiscoC3(vcmax,co2_intercel,can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2)

          ! C3: RuBP-limited photosynthesis
          aj = AgrossRuBPC3(par_abs,jmax,co2_intercel,co2_cpoint )

          ! Gross photosynthesis smoothing calculations. Co-limit ac and aj.
          ! RGK: We can remove this smoothing, right? theta is always nearly 1...?
          aquad = lb_params%theta_cj_c3
          bquad = -(ac + aj)
          cquad = ac * aj
          call QuadraticRoots(aquad, bquad, cquad, r1, r2)
          agross_out = min(r1,r2)

          ! RGK: agross_out = min(ac,aj)
          
       else

          ! C4: Rubisco-limited photosynthesis
          ac = vcmax

          ! C4: RuBP-limited photosynthesis
          aj = AgrossRuBPC4(par_abs)
            
          ! C4: PEP carboxylase-limited (CO2-limited)
          ap = AgrossPEPC4(co2_intercel,co2_rcurve_islope,can_press)
                      
          ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap

          aquad = lb_params%theta_cj_c4
          bquad = -(ac + aj)
          cquad = ac * aj
          call QuadraticRoots(aquad, bquad, cquad, r1, r2)
          ai = min(r1,r2)

          aquad = theta_ip
          bquad = -(ai + ap)
          cquad = ai * ap
          call QuadraticRoots(aquad, bquad, cquad, r1, r2)
          agross_out = min(r1,r2)
          
          ! RGK:  agross_out = minval([ac,aj,ap])
          
       end if

       ! Calculate anet
       anet_out = agross_out  - lmr

       if (  lb_params%stomatal_assim_model == gross_assim_model ) then
          if ( lb_params%stomatal_model == medlyn_model ) then
             write (fates_log(),*) 'Gross Assimilation conductance is incompatible with the Medlyn model'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          a_gs = agross_out
       else

          ! RGK: I'm not sure why we break out here.  Anet can be negative when the plant
          ! is performing photosynthesis, but respiration is larger than production, I see
          ! no reason to exit the solver early...
          if (anet_out < 0._r8) then
             loop_continue = .false.
          end if
          
          a_gs = anet_out
       end if

       ! A note about the use of the quadratic equations for calculating stomatal conductance
       ! ------------------------------------------------------------------------------------
       ! These two following models calculate the conductance between the intercellular leaf
       ! space and the leaf surface, not the canopy air space.  Transport between the leaf
       ! surface and the canopy air space is governed by the leaf boundary layer conductance.
       ! However, we need to estimate the properties at the surface of the leaf to solve for
       ! the stomatal conductance. We do this by using Fick's law (gradient resistance
       ! approximation of diffusion) to estimate the flux of water vapor across the
       ! leaf boundary layer, and balancing that with the flux across the stomata. It
       ! results in the following equation for leaf surface humidity:
       !
       ! e_s = (e_i g_s + e_c g_b)/(g_b + g_s)
       !
       ! The leaf surface humidity (e_s) becomes an expression of canopy humidity (e_c),
       ! intercellular humidity (e_i, which is the saturation humidity at leaf temperature),
       ! boundary layer conductance (g_b) (these are all known) and stomatal conductance
       ! (g_s) (this is still unknown).  This expression is substituted into the stomatal
       ! conductance equation. The resulting form of these equations becomes a quadratic.
       !
       ! For a detailed explanation, see the FATES technical note, section
       ! "1.11 Stomatal Conductance"
       !
       ! ------------------------------------------------------------------------------------

       leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb * anet_out * can_press 		   

       ! RGK: THIS 1e-6 CAP IS UNREASONBLY LOW. Units are Pascals. If the leaf internal CO2 is below 1
       ! pascal even, that is incredibly low. It should be some moderate fraction of
       ! atmospheric CO2 (which is roughly 40ish circa 2024), something that is much closer
       ! to that value than even 25% to 200%
       ! leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
       
       ! Determine saturation vapor pressure at the leaf surface, from temp and atm-pressure
       call QSat(veg_tempk, can_press, veg_qs, veg_esat)

       ! RGK: Should this be moved here?
       ! veg_esat = StomatalVaporPressureFromLWP(leaf_psi, k_lwp, veg_tempk, can_press, can_vpress)
       
       if ( lb_params%stomatal_model == medlyn_model ) then
          call StomatalCondMedlyn(anet_out,ft,veg_esat,can_vpress,stomatal_intercept_btran, &
                                  leaf_co2_ppress,can_press,gb,gs_out)
       else
          call StomatalCondBallBerry(a_gs,ft,veg_esat,can_vpress,stomatal_intercept_btran, &
                                     leaf_co2_ppress,can_press,gb,gs_out)
       end if

       ! Ensure that the conductance is above some nominal incredibly small value
       ! gs_out = max(gs_out,gsmin0_20C1A_mol)

       if(debug) then
          if(gs_out < stomatal_intercept_btran) then
             write(fates_log(),*) 'stomatal conductance is below the intercept, shouldnt be:',gs_out
             write(fates_log(),*) 'input:', lb_params%stomatal_model,anet_out,a_gs,ft,veg_esat,can_vpress
             write(fates_log(),*) 'input:',stomatal_intercept_btran,leaf_co2_ppress,can_press,gb
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if

       ! Derive new estimate for co2_intercel
       co2_intercel = can_co2_ppress - anet_out * can_press * &
            (h2o_co2_bl_diffuse_ratio*gs_out+h2o_co2_stoma_diffuse_ratio*gb) / (gb*gs_out)

       ! Check for co2_intercel convergence.
       ! The tolerance is if the new solution is less than 1 Pascal within
       ! the previous solution.  Typical values for the range are 20-50,
       ! with values less than atmospheric during production, and values
       ! greater than atmospheric during dark respiration

       if( abs(co2_intercel-co2_intercel_old) < co2_intercel_tol .or. niter >= max_iters) then
          loop_continue = .false.
       end if
    end do iter_loop

    ! estimate carbon 13 discrimination in leaf level carbon
    ! flux Liang WEI and Hang ZHOU 2018, based on
    ! Ubierna and Farquhar, 2014 doi:10.1111/pce.12346, using the simplified model:
    ! $\Delta ^{13} C = \alpha_s + (b - \alpha_s) \cdot \frac{C_i}{C_a}$
    ! just hard code b and \alpha_s for now, might move to parameter set in future
    ! b = 27.0 alpha_s = 4.4
    ! TODO, not considering C4 or CAM right now, may need to address this
    ! note co2_intercel is intracelluar CO2, not intercelluar
    c13disc_out = 4.4_r8 + (27.0_r8 - 4.4_r8) * &
         min (can_co2_ppress, max (co2_intercel, 0._r8)) / can_co2_ppress

    return
  end subroutine LeafLayerPhotosynthesis

  ! =======================================================================================

  function LeafHumidityStomaResis(leaf_psi, k_lwp, veg_tempk, can_vpress, can_press, &
       rb, gstoma, ft) result(rstoma_out)

    ! -------------------------------------------------------------------------------------
    ! This calculates inner leaf humidity as a function of mesophyll water potential 
    ! Adopted from  Vesala et al., 2017 https://www.frontiersin.org/articles/10.3389/fpls.2017.00054/full
    !
    ! Equation 1 in Vesala et al:
    ! lwp_star = wi/w0 = exp( k_lwp*leaf_psi*molar_mass_water/(rgas * veg_tempk) )
    !
    ! Terms:
    ! leaf_psi: leaf water potential [MPa]
    ! k_lwp: inner leaf humidity scaling coefficient [-]
    ! rgas: universal gas constant, [J/K/mol]
    ! molar_mass_water, molar mass of water, [g/mol]: 18.0
    !
    ! Unit conversions:
    ! 1 Pa = 1 N/m2 = 1 J/m3
    ! density of liquid water [kg/m3] = 1000
    ! 
    ! units of equation 1:  exp( [MPa]*[g/mol]/( [J/K/mol] * [K] ) )
    !                            [MJ/m3]*[g/mol]*[m3/kg]*[kg/g]*[J/MJ]  / ([J/mol])
    ! dimensionless:             [J/g]*[g/mol]/([J/mol])
    !
    ! Note: unit conversions drop out b/c [m3/kg]*[kg/g]*[J/MJ] = 1e-3*1.e-3*1e6 = 1.0
    !
    !
    ! RGK: Not clear to me, why this is not coupled into the stomatal conductance schemes.
    ! this should affect Anet...No?
    !   
    ! Junyan Ding 2021
    ! -------------------------------------------------------------------------------------

    ! Arguments
    real(r8) :: leaf_psi   ! Leaf water potential [MPa]
    real(r8) :: k_lwp      ! Scaling coefficient for the ratio of leaf xylem (user parameter)
    real(r8) :: veg_tempk  ! Leaf temperature     [K]
    real(r8) :: can_vpress ! vapor pressure of air (unconstrained) [Pa]
    real(r8) :: can_press  ! Atmospheric pressure of canopy [Pa]
    
    real(r8) :: rb         ! Leaf Boundary layer resistance [s/m]
    real(r8) :: gstoma     ! Stomatal Conductance of this leaf layer [m/s]
    integer  :: ft         ! Plant Functional Type
    real(r8) :: rstoma_out ! Total Stomatal resistance (stoma and BL) [s/m]
    
    ! Locals
    real(r8) :: ceair      ! vapor pressure of air, constrained [Pa]
                           ! water potential to mesophyll water potential
    real(r8) :: qs         ! Specific humidity [g/kg]
    real(r8) :: veg_esat   ! Saturated vapor pressure at veg surf [Pa]
    real(r8) :: qsat_alt   ! Saturation specific humidity  [g/kg]
    real(r8) :: qsat_loc   ! Saturation specific humidity  [g/kg]
    real(r8) :: qsat_adj   ! Adjusted saturation specific humidity  [g/kg]
    real(r8) :: lwp_star   ! leaf water potential scaling coefficient
    ! for inner leaf humidity, 0 means total dehydroted
    ! leaf, 1 means total saturated leaf

    ! Note: if k_lwp is an argument, but also a user defined parameter.
    ! If it is zero, LWP_star will be 1. 
    
    if (leaf_psi<0._r8) then
       lwp_star = exp(k_lwp*leaf_psi*molar_mass_water/(rgas_J_K_mol*veg_tempk))
    else 
       lwp_star = 1._r8
    end if
    
    call QSat(veg_tempk, can_press, qsat_alt, veg_esat)

    qsat_alt = qsat_alt * g_per_kg
    
    ceair = GetConstrainedVPress(can_vpress,veg_esat)
    
    ! compute specific humidity from vapor pressure
    ! q = molar_mass_ratio_vapdry*e/(can_press - (1-molar_mass_ratio_vapdry)*e) 
    ! source https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    ! now adjust inner leaf humidity by LWP_star

    qs = molar_mass_ratio_vapdry * ceair / (can_press - (1._r8-molar_mass_ratio_vapdry) * ceair)
    qsat_loc = molar_mass_ratio_vapdry * veg_esat / (can_press - (1._r8-molar_mass_ratio_vapdry) * veg_esat)
    qsat_adj = qsat_loc*lwp_star

    if(debug .and. (abs(qsat_loc-qsat_alt) > 1.e-2)) then
       write (fates_log(),*) 'qsat from QSat():', qsat_alt
       write (fates_log(),*) 'qsat localy :', qsat_loc
       write (fates_log(),*) 'values of qsat are too different'
       call endrun(msg=errMsg(sourcefile, __LINE__))  
    end if
    
    ! Adjusting gs (compute a virtual gs) that will be passed to host model

    if ( qsat_adj < qs ) then

       ! if inner leaf vapor pressure is less then or equal to that at leaf surface
       ! then set stomata resistance to be very large to stop the transpiration or back flow of vapor
       rstoma_out = rsmax0

    else

       rstoma_out = (qsat_loc-qs)*( 1/gstoma + rb)/(qsat_adj - qs)-rb

    end if

    if (rstoma_out < nearzero ) then
       write (fates_log(),*) 'qsat:', qsat_loc, 'qs:', qs
       write (fates_log(),*) 'LWP :', leaf_psi
       write (fates_log(),*) 'pft :', ft
       write (fates_log(),*) 'ceair:', ceair, 'veg_esat:', veg_esat            
       write (fates_log(),*) 'rstoma_out:', rstoma_out, 'rb:', rb  
       write (fates_log(),*) 'LWP_star', lwp_star 
       call endrun(msg=errMsg(sourcefile, __LINE__))                  
    end if

  end function LeafHumidityStomaResis

  ! =====================================================================================
  
  function StomatalVaporPressureFromLWP(leaf_psi, k_lwp, veg_tempk, &
                                        can_press, can_vpress) result(stoma_esat)

    ! -------------------------------------------------------------------------------------
    ! This calculates inner leaf humidity as a function of mesophyll water potential 
    ! Adopted from  Vesala et al., 2017
    ! https://www.frontiersin.org/articles/10.3389/fpls.2017.00054/full
    !
    ! Equation 1 in Vesala et al:
    ! lwp_star = wi/w0 = exp( k_lwp*leaf_psi*molar_mass_water/(rgas * veg_tempk) )
    !
    ! Terms:
    ! leaf_psi: leaf water potential [MPa]
    ! k_lwp: inner leaf humidity scaling coefficient [-]
    ! rgas: universal gas constant, [J/K/mol], 8.3144598
    ! molar_mass_water, molar mass of water, [g/mol]: 18.0
    !
    ! Junyan Ding 2021
    ! Adapted by Ryan Knox to use alternate functions like QSat() or inside
    ! Medlyn/BB solvers
    ! -------------------------------------------------------------------------------------

    ! Input
    real(r8) :: leaf_psi   ! Leaf water potential [MPa]
    real(r8) :: k_lwp      ! Scaling coefficient for the ratio of leaf xylem (user parameter)
    real(r8) :: veg_tempk  ! Leaf temperature     [K]
    real(r8) :: can_press  ! Atmospheric Pressure in the canopy [Pa]
    real(r8) :: can_vpress ! Actual Vapor Pressure in the canopy [Pa]
    
    ! Output
    real(r8) :: stoma_esat ! The vapor pressure at the surface of the stomata [Pa]

    
    real(r8) :: lwp_star   ! leaf water potential scaling coefficient
    real(r8) :: qs,es      ! saturation specific humidity [kg/kg] and vapor pressure [Pa]
    real(r8) :: stoma_qsat ! The specific humidity at the surface of the stomata [kg/kg]

    
    ! for inner leaf humidity:
    !    0 means total dehydrated leaf
    !    1 means total saturated leaf

    ! Note: if k_lwp is zero, LWP_star will be 1. 
    
    if (leaf_psi<0._r8) then
       lwp_star = exp(k_lwp*leaf_psi*molar_mass_water/(rgas_J_K_mol*veg_tempk))
    else 
       lwp_star = 1._r8
    end if

    call QSat(veg_tempk, can_press, qs, es)

    ! Calculate the reduced specific humidity based off of the scaled
    ! down saturated specific humidity at the leaf surface [kg/kg]
    stoma_qsat = lwp_star * qs

    ! Convert humidity to a vapor pressure [kg/kg] -> [Pa]
    ! qs = es * 0.622 / (p - 0.378*es)
    ! qs = 0.622 / (p/es - 0.378)
    ! qs (p/es - 0.378) = 0.622 
    ! qs*p/es  = 0.622 + qs*0.378
    ! 1/es  = (0.622 + qs*0.378)/(qs*p)
    ! es = (qs*p)/(0.622 + qs*0.378)

    stoma_esat = (stoma_qsat*can_press)/(0.622_r8 + stoma_qsat*0.378_r8)

    ! We don't allow reverse transpiration, prevent the stomatal surface
    ! vapor pressure from being lower than the actual!

    stoma_esat = max(can_vpress,stoma_esat)
    
    
  end function StomatalVaporPressureFromLWP


  
  ! =====================================================================================

  function ft1_f(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! This function, along with the fth1_f() function, scales the Vcmax25 and Jmax25
    ! constants so they are adjusted for leaf temperature.
    ! This equation does not have a name, but is equation 9.10 in the CLM5.0 technical
    ! manual.
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES
    

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas_J_K_mol*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft1_f

  ! =====================================================================================

  function fth_f(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! This function, along with the ft1_f() function, scales the Vcmax25 and Jmax25
    ! constants so they are adjusted for leaf temperature.
    ! This equation does not have a name, but is equation 9.10 in the CLM5.0 technical
    ! manual.
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temp function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temp function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temp inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas_J_K_mol*tl) ) )

    
    
    return
  end function fth_f

  ! =====================================================================================

  function fth25_f(hd,se)result(ans)
    !
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES


    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temp function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas_J_K_mol*(tfrz+25._r8)) )

    return
  end function fth25_f

  ! =====================================================================================

  subroutine GetCanopyGasParameters(can_press, &
       can_o2_partialpress, &
       veg_tempk, &
       mm_kco2,   &
       mm_ko2,    &
       co2_cpoint)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the specific Michaelis Menten Parameters (pa) for CO2
    ! and O2, as well as the CO2 compentation point.
    ! ---------------------------------------------------------------------------------

    use FatesConstantsMod, only: umol_per_mol
    use FatesConstantsMod, only: mmol_per_mol
    use FatesConstantsMod, only: umol_per_kmol

    ! Arguments
    real(r8), intent(in) :: can_press           ! Air pressure within the canopy (Pa)
    real(r8), intent(in) :: can_o2_partialpress ! Partial press of o2 in the canopy (Pa)
    real(r8), intent(in) :: veg_tempk           ! The temperature of the vegetation (K)

    real(r8), intent(out) :: mm_kco2       ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(out) :: mm_ko2        !  Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(out) :: co2_cpoint    !  CO2 compensation point (Pa)
                                           ! of conductance and resistance: [umol/m3]

    ! Locals
    real(r8) :: kc25                ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25                ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: sco                 ! relative specificity of rubisco
    real(r8) :: cp25                ! CO2 compensation point at 25C (Pa)

    ! ---------------------------------------------------------------------------------
    ! Intensive values (per mol of air)
    ! kc, ko, currentPatch, from: Bernacchi et al (2001)
    ! Plant, Cell and Environment 24:253-259
    ! ---------------------------------------------------------------------------------

    real(r8), parameter :: mm_kc25_umol_per_mol       = 404.9_r8
    real(r8), parameter :: mm_ko25_mmol_per_mol       = 278.4_r8
    real(r8), parameter :: co2_cpoint_umol_per_mol    = 42.75_r8

    ! Activation energy, from:
    ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
    ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430

    real(r8), parameter :: kcha    = 79430._r8  ! activation energy for kc (J/mol)
    real(r8), parameter :: koha    = 36380._r8  ! activation energy for ko (J/mol)
    real(r8), parameter :: cpha    = 37830._r8  ! activation energy for cp (J/mol)


    ! Derive sco from currentPatch and O2 using present-day O2 (0.209 mol/mol) and re-calculate
    ! currentPatch to account for variation in O2 using currentPatch = 0.5 O2 / sco

    ! FIXME (RGK 11-30-2016 THere are more constants here, but I don't have enough information
    ! about what they are or do, so I can't give them more descriptive names. Someone please
    ! fill this in when possible)

    kc25 = ( mm_kc25_umol_per_mol / umol_per_mol ) * can_press
    ko25 = ( mm_ko25_mmol_per_mol / mmol_per_mol ) * can_press
    sco  = 0.5_r8 * 0.209_r8 / (co2_cpoint_umol_per_mol / umol_per_mol )
    cp25 = 0.5_r8 * can_o2_partialpress / sco

    if( veg_tempk.gt.150_r8 .and. veg_tempk.lt.350_r8 )then
       mm_kco2       = kc25 * ft1_f(veg_tempk, kcha)
       mm_ko2         = ko25 * ft1_f(veg_tempk, koha)
       co2_cpoint     = cp25 * ft1_f(veg_tempk, cpha)
    else
       mm_kco2    = 1.0_r8
       mm_ko2     = 1.0_r8
       co2_cpoint = 1.0_r8
    end if

    return
  end subroutine GetCanopyGasParameters


  ! ====================================================================================
  
  subroutine LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, &
       nscaler,       &
       ft,            &
       veg_tempk,     &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : umolC_to_kgC
    

    ! -----------------------------------------------------------------------
    ! Base maintenance respiration rate for plant tissues maintresp_leaf_ryan1991_baserate
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    ! Which is the default value of maintresp_nonleaf_baserate

    ! Arguments
    real(r8), intent(in)  :: lnc_top      ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8), intent(in)  :: nscaler      ! Scale for leaf nitrogen profile
    integer,  intent(in)  :: ft           ! (plant) Functional Type Index
    real(r8), intent(in)  :: veg_tempk    ! vegetation temperature
    real(r8), intent(out) :: lmr          ! Leaf Maintenance Respiration  (umol CO2/m**2/s)

    ! Locals
    real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top  ! canopy top leaf maint resp rate at 25C for this pft (umol CO2/m**2/s)

    ! Parameter
    real(r8), parameter :: lmrha = 46390._r8    ! activation energy for lmr (J/mol)
    real(r8), parameter :: lmrhd = 150650._r8   ! deactivation energy for lmr (J/mol)
    real(r8), parameter :: lmrse = 490._r8      ! entropy term for lmr (J/mol/K)
    real(r8), parameter :: lmrc = 1.15912391_r8 ! scaling factor for high
    ! temperature inhibition (25 C = 1.0)

    lmr25top = lb_params%maintresp_leaf_ryan1991_baserate(ft) * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
    lmr25top = lmr25top * lnc_top / (umolC_to_kgC * g_per_kg)


    ! Part I: Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
    ! ----------------------------------------------------------------------------------
    lmr25 = lmr25top * nscaler


    if (lb_params%c3psn(ft) == c3_path_index) then
       ! temperature sensitivity of C3 plants
       lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
            fth_f(veg_tempk, lmrhd, lmrse, lmrc)
    else
       ! temperature sensitivity of C4 plants
       lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
       lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
    endif

    ! Any hydrodynamic limitations could go here, currently none
    ! lmr = lmr * (nothing)

  end subroutine LeafLayerMaintenanceRespiration_Ryan_1991

  ! ====================================================================================   

  subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top, &
       rdark_scaler,     &
       ft,             &
       veg_tempk,      &
       tgrowth,        &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm

    ! Arguments
    real(r8), intent(in)  :: lnc_top          ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8), intent(in)  :: rdark_scaler     ! Decay coefficient for vertical scaling of respiration. See note 1 below
    integer,  intent(in)  :: ft               ! (plant) Functional Type Index
    real(r8), intent(in)  :: veg_tempk        ! vegetation temperature  (degrees K)
    real(r8), intent(in)  :: tgrowth          ! lagged vegetation temperature averaged over acclimation timescale (degrees K)
    real(r8), intent(out) :: lmr              ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
    
    
    ! Locals
    real(r8) :: r_0     ! base respiration rate, PFT-dependent (umol CO2/m**2/s)
    real(r8) :: r_t_ref ! acclimated ref respiration rate (umol CO2/m**2/s)
  
    ! parameter values of r_0 as listed in Atkin et al 2017: (umol CO2/m**2/s) 
    ! Broad-leaved trees  1.7560
    ! Needle-leaf trees   1.4995
    ! Shrubs              2.0749
    ! C3 herbs/grasses    2.1956
    ! In the absence of better information, we use the same value for C4 grasses as C3 grasses.

    ! r_0 currently put into the EDPftvarcon_inst%dev_arbitrary_pft
    ! all figs in Atkin et al 2017 stop at zero Celsius so we will assume acclimation is fixed below that
    r_0 = lb_params%maintresp_leaf_atkin2017_baserate(ft)

    ! Note 1:
    ! This code uses the relationship between leaf N and respiration from Atkin et al 
    ! for the top of the canopy, but then scales through the canopy based on a rdark_scaler.
    ! To assume proportionality with N through the canopy following Lloyd et al. 2010, use the
    ! default parameter value of 2.43, which results in the scaling of photosynthesis and respiration
    ! being proportional through the canopy. To have a steeper decrease in respiration than photosynthesis
    ! this number can be smaller. There is some observational evidence for this being the case
    ! in Lamour et al. 2023.

    r_t_ref = max(0._r8, rdark_scaler * (r_0 + lmr_r_1 * lnc_top + lmr_r_2 * max(0._r8, (tgrowth - tfrz) )) )

    if (r_t_ref .eq. 0._r8) then
       warn_msg = 'Rdark is negative at this temperature and is capped at 0. tgrowth (C): ' &
            //trim(N2S(tgrowth-tfrz))//' pft: '//trim(I2S(ft))
       call FatesWarn(warn_msg,index=4)            
    end if

    lmr = r_t_ref * exp(lmr_b * (veg_tempk - tfrz - lmr_TrefC) + lmr_c * &
         ((veg_tempk-tfrz)**2 - lmr_TrefC**2))

  end subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017

  ! ====================================================================================

  subroutine LeafLayerBiophysicalRates( ft,            &
       vcmax25top_ft, &
       jmax25top_ft, &
       co2_rcurve_islope25top_ft, &
       nscaler,    &
       veg_tempk,      &
       dayl_factor, &
       t_growth,   &
       t_home,     &
       btran, &
       vcmax, &
       jmax, &
       co2_rcurve_islope )

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the localized rates of several key photosynthesis
    ! rates.  By localized, we mean specific to the plant type and leaf layer,
    ! which factors in leaf physiology, as well as environmental effects.
    ! This procedure should be called prior to iterative solvers, and should
    ! have pre-calculated the reference rates for the pfts before this.
    !
    ! The output biophysical rates are:
    ! vcmax: maximum rate of carboxilation,
    ! jmax: maximum electron transport rate,
    ! co2_rcurve_islope: initial slope of CO2 response curve (C4 plants)
    ! ---------------------------------------------------------------------------------

    ! Arguments
    ! ------------------------------------------------------------------------------

    integer,  intent(in) :: ft                        ! (plant) Functional Type Index
    real(r8), intent(in) :: nscaler                   ! Scale for leaf nitrogen profile
    real(r8), intent(in) :: vcmax25top_ft             ! canopy top maximum rate of carboxylation at 25C
                                                      ! for this pft (umol CO2/m**2/s)
    real(r8), intent(in) :: jmax25top_ft              ! canopy top maximum electron transport rate at 25C
                                                      ! for this pft (umol electrons/m**2/s)
    real(r8), intent(in) :: co2_rcurve_islope25top_ft ! initial slope of CO2 response curve
                                                      ! (C4 plants) at 25C, canopy top, this pft
    real(r8), intent(in) :: veg_tempk                 ! vegetation temperature
    real(r8), intent(in) :: dayl_factor               ! daylength scaling factor (0-1)
    real(r8), intent(in) :: t_growth                  ! T_growth (short-term running mean temperature) (K)
    real(r8), intent(in) :: t_home                    ! T_home (long-term running mean temperature) (K)
    real(r8), intent(in) :: btran                     ! transpiration wetness factor (0 to 1)
    
    real(r8), intent(out) :: vcmax                    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(out) :: jmax                     ! maximum electron transport rate
                                                      ! (umol electrons/m**2/s)
    real(r8), intent(out) :: co2_rcurve_islope        ! initial slope of CO2 response curve (C4 plants)

    ! Locals
    ! -------------------------------------------------------------------------------
    real(r8) :: vcmax25             ! leaf layer: maximum rate of carboxylation at 25C
    ! (umol CO2/m**2/s)
    real(r8) :: jmax25              ! leaf layer: maximum electron transport rate at 25C
    ! (umol electrons/m**2/s)
    real(r8) :: co2_rcurve_islope25 ! leaf layer: Initial slope of CO2 response curve
    ! (C4 plants) at 25C
    real(r8) :: dayl_factor_local   ! Local version of daylength factor

    ! Parameters
    ! ---------------------------------------------------------------------------------
    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: t_growth_celsius ! average growing temperature
    real(r8) :: t_home_celsius   ! average home temperature
    real(r8) :: jvr            ! ratio of Jmax25 / Vcmax25
    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)


    ! update the daylength factor local variable if the switch is on
    if ( lb_params%dayl_switch == itrue ) then
       dayl_factor_local = dayl_factor
    else
       dayl_factor_local = 1.0_r8
    endif
    
    select case(lb_params%photo_tempsens_model)
    case (photosynth_acclim_model_none) !No temperature acclimation
       vcmaxha = lb_params%vcmaxha(FT)
       jmaxha  = lb_params%jmaxha(FT)
       vcmaxhd = lb_params%vcmaxhd(FT)
       jmaxhd  = lb_params%jmaxhd(FT)
       vcmaxse = lb_params%vcmaxse(FT)
       jmaxse  = lb_params%jmaxse(FT)
       
    case (photosynth_acclim_model_kumarathunge_etal_2019) !Kumarathunge et al. temperature acclimation, Thome=30-year running mean
       t_growth_celsius = t_growth-tfrz
       t_home_celsius = t_home-tfrz
       vcmaxha = (42.6_r8 + (1.14_r8*t_growth_celsius))*1e3_r8 !J/mol
       jmaxha = 40.71_r8*1e3_r8 !J/mol
       vcmaxhd = 200._r8*1e3_r8 !J/mol
       jmaxhd = 200._r8*1e3_r8 !J/mol
       vcmaxse = (645.13_r8 - (0.38_r8*t_growth_celsius))
       jmaxse = 658.77_r8 - (0.84_r8*t_home_celsius) - 0.52_r8*(t_growth_celsius-t_home_celsius)
       jvr = 2.56_r8 - (0.0375_r8*t_home_celsius)-(0.0202_r8*(t_growth_celsius-t_home_celsius))

       
    case default
       write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
       write (fates_log(),*)'lb_params%photo_tempsens_model: ',lb_params%photo_tempsens_model
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    vcmaxc = fth25_f(vcmaxhd, vcmaxse)
    jmaxc  = fth25_f(jmaxhd, jmaxse)
    
    ! Vcmax25top was already calculated to derive the nscaler function
    vcmax25 = vcmax25top_ft * nscaler * dayl_factor_local
    
    select case( lb_params%photo_tempsens_model)
    case (photosynth_acclim_model_none)
       jmax25  = jmax25top_ft * nscaler * dayl_factor_local
    case (photosynth_acclim_model_kumarathunge_etal_2019) 
       jmax25 = vcmax25*jvr
    case default
       write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
       write (fates_log(),*)'lb_params%photo_tempsens_model:',lb_params%photo_tempsens_model
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    
    co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler
    
    ! Adjust for temperature
    ! photosynthetic pathway: 0. = c4, 1. = c3
    
    if (lb_params%c3psn(ft) == c3_path_index) then
       vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
    else
       vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
       vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
       vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
    end if
    
    jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)
    
    !q10 response of product limited psn.
    co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
 
 
    ! Adjust for water limitations
    vcmax = vcmax * btran

    return

  end subroutine LeafLayerBiophysicalRates

  ! =====================================================================================
  
  function DecayCoeffVcmax(vcmax25top,slope_param,intercept_param) result(decay_coeff_vcmax)
    
    ! ---------------------------------------------------------------------------------
    ! This function estimates the decay coefficient used to estimate vertical
    ! attenuation of properties in the canopy.
    !
    ! Decay coefficient (kn) is a function of vcmax25top for each pft.
    !
    ! Currently, this decay is applied to vcmax attenuation, SLA (optionally)
    ! and leaf respiration (optionally w/ Atkin)
    !
    ! ---------------------------------------------------------------------------------
    
    !ARGUMENTS

    real(r8),intent(in) :: vcmax25top
    real(r8),intent(in) :: slope_param      ! multiplies vcmax25top
    real(r8),intent(in) :: intercept_param  ! adds to vcmax25top

    real(r8) :: decay_coeff_vcmax
    
    !LOCAL VARIABLES
    ! -----------------------------------------------------------------------------------
    
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
    ! kn = 0.11. Here, we derive kn from vcmax25 as in Lloyd et al 
    ! (2010) Biogeosciences, 7, 1833-1859
    ! This function is also used to vertically scale leaf maintenance
    ! respiration.
    
    decay_coeff_vcmax = exp(slope_param * vcmax25top - intercept_param)
    
    return
  end function DecayCoeffVcmax
  
  ! =====================================================================================

  subroutine LowstorageMainRespReduction(frac, pft, maintresp_reduction_factor)

    ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
    ! of this is that mortality of plants increases when storage is low because they are not able
    ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
    ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1,
    ! and the curvature of the function is determined by a parameter.

    ! Arguments
    ! ------------------------------------------------------------------------------
    real(r8), intent(in) :: frac      ! ratio of storage to target leaf biomass
    integer,  intent(in) :: pft       ! what pft is this cohort?
    real(r8), intent(out) :: maintresp_reduction_factor  ! the factor by which to reduce maintenance respiration

    ! --------------------------------------------------------------------------------
    ! Parameters are at the PFT level:
    ! fates_maintresp_reduction_curvature controls the curvature of this.
    ! If this parameter is zero, then there is no reduction until the plant dies at storage = 0.
    ! If this parameter is one, then there is a linear reduction in respiration below the storage point.
    ! Intermediate values will give some (concave-downwards) curvature.
    !
    ! maintresp_reduction_intercept controls the maximum amount of throttling.
    ! zero means no throttling at any point, so it turns this mechanism off completely and so
    ! allows an entire cohort to die via negative carbon-induced termination mortality.
    ! one means complete throttling, so no maintenance respiration at all, when out of carbon.
    ! ---------------------------------------------------------------------------------

    if( frac .lt. 1._r8 )then
       if ( abs(lb_params%maintresp_reduction_curvature(pft)-1._r8) > nearzero ) then
          maintresp_reduction_factor = (1._r8 - lb_params%maintresp_reduction_intercept(pft)) + &
               lb_params%maintresp_reduction_intercept(pft) * &
               (1._r8 - lb_params%maintresp_reduction_curvature(pft)**frac) &
               / (1._r8-lb_params%maintresp_reduction_curvature(pft))
       else  ! avoid nan answer for linear case
          maintresp_reduction_factor = (1._r8 - lb_params%maintresp_reduction_intercept(pft)) + &
               lb_params%maintresp_reduction_intercept(pft) * frac
       endif

    else
       maintresp_reduction_factor = 1._r8
    endif


  end subroutine LowstorageMainRespReduction



  ! =====================================================================================

  real(r8) function GetConstrainedVPress(air_vpress,veg_esat) result(ceair)

    ! Return a constrained vapor pressure [Pa]
    
    real(r8) :: air_vpress   ! vapor pressure of the air (unconstrained) [Pa]
    real(r8) :: veg_esat     ! saturated vapor pressure [Pa]

    ! Constrain eair >= 0.05*esat_tv so that solution does not blow up. This ensures
    ! that hs (ie air_vpress / veg_esat)  in the Ball-Berry equation does not go to zero.
    ! Also eair <= veg_esat so that hs <= 1
    
    ceair = min( max(air_vpress, 0.05_r8*veg_esat ),veg_esat )
    
  end function GetConstrainedVPress

   ! ====================================================================================
  
  real(r8) function GetStomatalInterceptBtran(ft,btran) result(stomatal_intercept_btran)

    ! This helper function returns the stomatal intercept for conductance equations.
    ! The core mechanic here is to multiply the user parameter by the moisture limitation
    ! function "btran". It also places a lower limit on this based off of the minimum
    ! allowable conductance.
    
    ! Arguments
    integer  :: ft
    real(r8) :: btran

    stomatal_intercept_btran = max(gsmin0_20C1A_mol,lb_params%stomatal_intercept(ft)*btran )
    
  end function GetStomatalInterceptBtran

  ! =====================================================================================

  subroutine QSat (T, p, qs, es, qsdT, esdT)
    !
    ! !DESCRIPTION:
    !
    ! THIS IS AN EXACT CLONE OF QSat routine in QSatMod.F90 from CTSM, tag ctsm5.2.028
    ! https://github.com/ESCOMP/CTSM/blob/ctsm5.2.028/src/biogeophys/QSatMod.F90
    ! 
    ! Computes saturation mixing ratio and (optionally) the change in saturation mixing
    ! ratio with respect to temperature. Mixing ratio and specific humidity are
    ! approximately equal and can be treated as the same.
    ! Reference:  Polynomial approximations from:
    !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
    !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
    
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out), optional :: es       ! vapor pressure (pa)
    real(r8), intent(out), optional :: qsdT     ! d(qs)/d(T)
    real(r8), intent(out), optional :: esdT     ! d(es)/d(T)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: es_local    ! local version of es (in case es is not present)
    real(r8) :: esdT_local  ! local version of esdT (in case esdT is not present)
    real(r8) :: td,vp,vp1,vp2

    ! For water vapor (temperature range 0C-100C)
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
    ! For derivative:water vapor
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
    ! For ice (temperature range -75C-0C)
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
    ! For derivative:ice
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
    
    
    !-----------------------------------------------------------------------

    td = min(100.0_r8, max(-75.0_r8, T - tfrz ))

    if (td >= 0.0_r8) then
       es_local = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
       es_local = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    endif

    es_local = es_local * 100._r8            ! pa
    vp    = 1.0_r8   / (p - 0.378_r8*es_local)
    vp1   = 0.622_r8 * vp
    qs    = es_local * vp1             ! kg/kg
    if (present(es)) then
       es = es_local
    end if

    if (present(qsdT) .or. present(esdT)) then
       if (td >= 0.0_r8) then
          esdT_local = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
               + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
       else
          esdT_local = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
               + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
       end if

       esdT_local = esdT_local * 100._r8            ! pa/K
       vp2 = vp1 * vp
       if (present(qsdT)) then
          qsdT = esdT_local * vp2 * p         ! 1 / K
       end if
       if (present(esdT)) then
          esdT = esdT_local
       end if
    end if

  end subroutine QSat

  ! =====================================================================================
  
  real(r8) function VeloToMolarCF(press,tempk) result(cf)

    ! ---------------------------------------------------------------------------------
    !
    ! "cf" is the conversion factor between molar form and velocity form
    ! of conductance and resistance. The units on this factor are: [umol/m3]
    ! This uses the ideal gas law. This routine is necessary because the
    ! photosynthesis module uses units of moles, while the land-energy
    ! balance solvers in the host models use units of velocity.
    !
    ! i.e.
    ! [m/s] * [umol/m3] -> [umol/m2/s]
    !
    ! Breakdown of the conversion factor: [ umol / m3 ]
    !
    ! Rgas [J /K /kmol]
    ! Air Potential Temperature [ K ]
    ! Air Pressure      [ Pa ]
    ! conversion: umol/kmol =  1e9
    !
    ! [ Pa * K * kmol umol/kmol  /  J K ] = [ Pa * umol / J ]
    ! since: 1 Pa = 1 N / m2
    ! [ Pa * umol / J ] = [ N * umol / J m2 ]
    ! since: 1 J = 1 N * m
    ! [ N * umol / J m2 ] = [ N * umol / N m3 ]
    ! [ umol / m3 ]
    !
    ! --------------------------------------------------------------------------------

    ! Arguments
    real(r8) :: press    ! air pressure at point of interest [Pa]
    real(r8) :: tempk    ! temperature at point of interest  [K]
    
    cf = press/(rgas_J_K_kmol * tempk )*umol_per_kmol
    
  end function VeloToMolarCF

  ! =====================================================================================
  
end module LeafBiophysicsMod
