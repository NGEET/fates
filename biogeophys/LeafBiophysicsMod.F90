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
  use FatesConstantsMod, only : nocomp_bareground
  use FatesConstantsMod, only : lmrmodel_ryan_1991
  use FatesConstantsMod, only : lmrmodel_atkin_etal_2017
  use FatesConstantsMod, only : kpa_per_pa
  use FatesConstantsMod, only : umol_per_kmol
  use FatesUtilsMod,     only : QuadraticRoots => QuadraticRootsSridharachary
  use FatesConstantsMod, only : rgas_J_K_kmol
  use FatesConstantsMod, only : rgas_J_K_mol
  use FatesConstantsMod, only : g_per_kg
  use FatesConstantsMod, only : umolC_to_kgC
  
  implicit none
  private

  public :: LeafLayerPhotosynthesis
  public :: LeafHumidityStomaResis
  public :: GetCanopyGasParameters
  public :: LeafLayerMaintenanceRespiration_Ryan_1991
  public :: LeafLayerMaintenanceRespiration_Atkin_etal_2017
  public :: LeafLayerBiophysicalRates
  public :: LowstorageMainRespReduction
  public :: GetConstrainedVPress
  public :: DecayCoeffVcmax
  public :: QSat
  public :: AgrossRubiscoC3
  public :: AgrossRuBPC3
  public :: AgrossRuBPC4
  public :: AgrossPEPC4
  public :: StomatalCondMedlyn
  public :: StomatalCondBallBerry
  public :: VeloToMolarCF
  public :: CiMinMax
  public :: CiFunc
  public :: CiBisection
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  character(len=1024) :: warn_msg   ! for defining a warning message

  !-------------------------------------------------------------------------------------

  ! maximum stomatal resistance [s/m]
  real(r8),public, parameter :: rsmax0 =  2.e8_r8
  
  ! minimum allowable conductance for numerics purposes at 20C (293.15K)
  ! and 1 standard atmosphere (101325 Pa) [umol/m2/s]
  ! rgas_J_K_kmol 8314.4598
  ! this follows the same math as  FatesPlantRespPhotosynthMod:FetMolarVeloCF()
  ! 101325.0_r8/(8314.4598 * 293.15 )*1.e9/2.e8  ~ 0.2 [umol/m2/s]
  real(r8),parameter :: gsmin0 = 101325.0_r8/(rgas_J_K_kmol * 293.15 )*umol_per_kmol/rsmax0

  ! minimum allowable conductance for getting reasonable
  ! bounds on extreme equations [umol/m2/s]
  ! real(r8),parameter :: gsmin0 = 50._r8
  
  ! Set this to true to perform debugging
  logical,parameter   ::  debug = .false.

  ! Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
  real(r8),parameter :: h2o_co2_stoma_diffuse_ratio = 1.6_r8

  ! Ratio of H2O/CO2 gass diffusion in the leaf boundary layer (approximate)
  real(r8),parameter :: h2o_co2_bl_diffuse_ratio = 1.4_r8

  ! Constants used to define C3 versus C4 photosynth pathways
  integer, public, parameter :: c3_path_index = 1
  integer, public, parameter :: c4_path_index = 0


  ! Constants used to define conductance models
  integer, parameter :: medlyn_model = 2
  integer, parameter :: ballberry_model = 1

  ! Alternatively, Gross Assimilation can be used to estimate
  ! leaf co2 partial pressure and therefore conductance. The default
  ! is to use anet
  integer, parameter :: net_assim_model = 1
  integer, parameter :: gross_assim_model = 2

  ! Constants defining the photosynthesis temperature acclimation model
  integer, parameter :: photosynth_acclim_model_none = 0
  integer, parameter :: photosynth_acclim_model_kumarathunge_etal_2019 = 1

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


  ! curvature parameter for quadratic smoothing of C4 gross assimilation
  real(r8),parameter :: theta_ip_c4 = 0.95_r8  !0.95 is from Collatz 91, 0.999 was api 36
  real(r8),parameter :: theta_cj_c4 = 0.98_r8  !0.98 from Collatz 91,  0.099 was api 36

  ! This bypasses a potential fix to assimilation in BB
  logical, parameter :: bb_agsfix_bypass = .true.

  ! Some observations have indicated that vcmax and jmax
  ! never really drop down to 0, or close to it. Setting
  ! this to true will cap vcmax and jmax to a fraction
  ! of its value at 25C. 
  logical, parameter :: do_mincap_vcjmax = .false.
  real(r8),parameter :: min_vcmax_frac = 0.10_r8
  real(r8),parameter :: min_jmax_frac  = 0.10_r8

  
  ! For plants with no leaves, a miniscule amount of conductance
  ! can happen through the stems, at a partial rate of cuticular conductance
  real(r8),parameter :: stem_cuticle_loss_frac = 0.1_r8


  ! The stomatal slope can either be scaled by btran or not. FATES had
  ! a precedent of using this into 2024, but discussions here: https://github.com/NGEET/fates/issues/719
  ! suggest we should try other hypotheses

  integer, parameter,public :: btran_on_gs_none = 0       ! do not apply btran to conductance terms
  integer, parameter,public :: btran_on_gs_gs0 = 1        ! apply btran to stomatal intercept (API 36.1)
  integer, parameter,public :: btran_on_gs_gs1 = 2        ! apply btran to stomatal slope
  integer, parameter,public :: btran_on_gs_gs01 = 3       ! apply btran to both stomatal slope and intercept
  integer, parameter,public :: btran_on_gs_gs2  = 4       ! apply btran to the whole non-intercept portion
                                                          ! of conductance equation. (NOTE! for Medlyn,
                                                          ! this is different than btran_on_gs_gs1,
                                                          ! for Ball-Berry, this is the SAME as
                                                          ! btran_on_gs_gs1
  integer, parameter,public :: btran_on_gs_gs02 = 5       ! same as btran_on_gs_gs2, but also apply to the intercept
  
  integer, parameter,public :: btran_on_ag_none  = 0      ! do not apply btran to vcmax or jmax
  integer, parameter,public :: btran_on_ag_vcmax = 1      ! apply btran to vcmax (API 36.1)
  integer, parameter,public :: btran_on_ag_vcmax_jmax = 2 ! apply btran to vcmax and jmax
  
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
     integer              :: dayl_switch                          ! switch for turning on or off day length factor
                                                                  ! scaling for photosynthetic parameters
     integer              :: photo_tempsens_model                 ! switch for choosing the model that defines the temperature
                                                                  ! sensitivity of photosynthetic parameters (vcmax, jmax).
                                                                  ! 0=non-acclimating, 1=Kumarathunge et al., 2019
     integer              :: stomatal_assim_model                 ! Switch designating whether to use net or
                                                                  ! gross assimilation in the stomata model
     integer              :: stomatal_model                       ! switch for choosing between stomatal conductance models,
                                                                  ! for Ball-Berry, 2 for Medlyn
     integer,allocatable :: stomatal_btran_model(:)               ! index for how btran effects conductance
                                                                  ! 0: btran does not scale the stomatal slope or intercept
                                                                  ! 1: btran scales the stomatal intercept only
                                                                  ! 2: btran scales the stomatal slope only
                                                                  ! 3: btran scales both stomatal slope and intercept
     integer,allocatable :: agross_btran_model(:)                 ! index for how btran scales gross assimilation processes
                                                                  ! 0: btran does not scale vcmax or jmax
                                                                  ! 1: btran scales only vcmax
                                                                  ! 2: btran scales both vcmax and jmax
     
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
  
  subroutine StomatalCondMedlyn(anet,veg_esat,can_vpress,gs0,gs1,gs2,leaf_co2_ppress,can_press,gb,gs)

    ! -----------------------------------------------------------------------------------
    ! Calculate stomatal conductance (gs [umol/m2/s]) via the Medlyn approach, described in:
    ! Medlyn et al. Reconciling the optimal and empirical approaches to modelling stomatal
    ! conductance.  Global Change Biology (2011) 17, 2134–2144, doi: 10.1111/j.1365-2486.2010.02375.x
    !
    ! Original implementation in FATES is described in:
    ! Li, Q. and Serbin, S. P. and Lamour, J. and Davidson, K. J. and Ely, K. S. and Rogers, A.
    !  Implementation and evaluation of the unified stomatal optimization approach in the Functionally
    !  Assembled Terrestrial Ecosystem Simulator (FATES). Geoscientific Model Development, 15(11), 2022.
    !  doi: 10.5194/gmd-15-4313-2022.
    !
    ! The implementation, and adaptation to include a leaf boundary layer in serial
    ! resitance was taken from CLM5.0
    ! -----------------------------------------------------------------------------------

    
    ! Input
    real(r8), intent(in) :: anet            ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(in) :: veg_esat        ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press       ! Air pressure NEAR the surface of the leaf (Pa)
    real(r8), intent(in) :: gb              ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress      ! vapor pressure of canopy air (Pa)
    real(r8), intent(in) :: gs0,gs1,gs2     ! stomatal intercept, and two different
                                            ! mutually exclusive slopes
                                            ! gs1 is a slope that is possibly multiplied by btran
                                            ! gs2 is an alternate location for btran scaling
                                            ! that applies to the non-intercept portion of
                                            ! the equation
    real(r8), intent(in) :: leaf_co2_ppress ! CO2 partial pressure at the leaf surface (Pa)
                                            ! at the boundary layer interface with the leaf

    ! Output
    real(r8),intent(out) :: gs                       ! leaf stomatal conductance (umol H2O/m**2/s)

    ! locals
    real(r8) :: vpd                                  ! water vapor deficit in Medlyn stomatal model [KPa]
    real(r8) :: term                                 ! intermediate term used to simplify equations
    real(r8) :: aquad,bquad,cquad                    ! quadradic solve terms
    real(r8) :: r1,r2                                ! quadradic solve roots
    real(r8),parameter :: min_vpd_pa = 50._r8        ! Minimum allowable vpd [Pa]
    real(r8),parameter :: anet_min = 0.001_r8
    logical :: err
    
    ! Evaluate trival solution, if there is no positive net assimiolation
    ! the stomatal conductance is the intercept conductance
    if (anet <= anet_min) then
       gs = gs0
       return
    end if

    ! Trivial case (gs2 near 0)
    if(gs2<0.01_r8) then
       gs = gs0
       return
    end if
    
    ! Trivial case (gs1 near 0)
    if(gs1<0.01_r8)then
       gs = gs0 + h2o_co2_stoma_diffuse_ratio * anet/(leaf_co2_ppress/can_press)
       return
    end if
    
    vpd =  max((veg_esat - can_vpress), min_vpd_pa) * kpa_per_pa
    term = gs2 * h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
    aquad = 1.0_r8
    bquad = -(2.0 * (gs0 + term) + (gs1 * term)**2 / (gb * vpd ))
    cquad = gs0*gs0 +  (2.0*gs0 + term * &
         (1.0 - gs1*gs1 / vpd)) * term
    
    call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)
    gs = max(r1,r2)

    if(debug)then
       if(err)then
          write(fates_log(),*) "medquadfail:",anet,veg_esat,can_vpress,gs0,gs1,gs2,leaf_co2_ppress,can_press
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    ! RGK: re-derived solution to check units.
    ! For 50*10*10*5*50 operations in unit tests the original method is
    ! about 2% faster. Alternative derivation below.
    
    ! vpd =  max((veg_esat - can_vpress), min_vpd_pa)
    ! term = h2o_co2_stoma_diffuse_ratio * anet * can_press / leaf_co2_ppress
    ! term2 = pa_per_kpa*(gs1*term)**2._r8
    ! aquad = gb*vpd
    ! bquad = -2*gb*vpd*(gs0 + term) - term2
    ! cquad = gb*vpd*(gs0 + term)**2._r8 - gb*term2

    ! call QuadraticRoots(aquad, bquad, cquad, r3, r4)

    ! if( abs(max(r1,r2)-max(r3,r4))>1.e-5 ) then
    !    print*,"Math check failed",r1,r2,r3,r4,anet,veg_esat,can_vpress,gs0,gs1,leaf_co2_ppress,can_press,gb
    !    stop
    ! end if
       
    return
  end subroutine StomatalCondMedlyn

  ! =======================================================================================

  subroutine StomatalCondBallBerry(a_gs,a_net,veg_esat,can_vpress,gs0,gs1,leaf_co2_ppress,can_press, gb, gs)


    ! Input
    real(r8), intent(in) :: veg_esat                 ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press                ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in) :: gb                       ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress               ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in) :: gs0,gs1                  ! Stomatal intercept and slope (umol H2O/m**2/s)
    real(r8), intent(in) :: leaf_co2_ppress          ! CO2 partial pressure at leaf surface (Pa)
    real(r8), intent(in) :: a_gs                     ! The assimilation (a) for calculating conductance (gs)
                                                     ! is either = to anet or agross
    real(r8), intent(in) :: a_net                    ! Net assimilation rate of co2 (umol/m2/s)
                                                     ! Output
    real(r8), intent(out) :: gs                      ! leaf stomatal conductance (umol H2O/m**2/s)

                                                     ! locals
    real(r8) :: ceair                                ! constrained vapor pressure (Pa)
    real(r8) :: aquad,bquad,cquad                    ! quadradic solve terms
    real(r8) :: r1,r2                                ! quadradic solve roots
    logical :: err


    if (a_gs <= nearzero) then
       gs = gs0
       return
    end if
       
    ! Trivial case (gs1 near 0)
    if(gs1<0.01_r8)then
       gs = gs0
       return
    end if
        
    ! Apply a constraint to the vapor pressure
    ceair = GetConstrainedVPress(can_vpress,veg_esat)
    
    aquad = leaf_co2_ppress
    bquad = leaf_co2_ppress*(gb - gs0) - gs1 * a_gs * can_press
    
    if(bb_agsfix_bypass) then
       cquad = -gb*(leaf_co2_ppress * gs0 + gs1 * a_net * can_press * ceair/ veg_esat )
    else
       cquad = -gb*(leaf_co2_ppress * gs0 + gs1 * a_gs* can_press * ceair/ veg_esat )
    end if
    
    call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

    if(debug)then
       if(err)then
          write(fates_log(),*) "bbquadfail:",a_net,a_gs,veg_esat,can_vpress,gs0,gs1,leaf_co2_ppress,can_press
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    gs = max(r1,r2)
    
    return
  end subroutine StomatalCondBallBerry

  ! =====================================================================================
  
  function AgrossRubiscoC3(vcmax,ci,can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2) result(ac)

    ! Input
    real(r8) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8) :: ci                ! intracellular leaf CO2 (Pa)
    real(r8) :: co2_cpoint        ! CO2 compensation point (Pa)
    real(r8) :: can_o2_ppress     ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8) :: mm_kco2           ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2            ! Michaelis-Menten constant for O2 (Pa)
    
    ! Output
    real(r8) :: ac               ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    
    ac = vcmax * max(ci-co2_cpoint, 0._r8) / &
         (ci+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))
    
  end function AgrossRubiscoC3
  
  ! =====================================================================================

  function GetJe(par_abs,jmax) result(je)

    ! Input
    real(r8) :: par_abs           ! Absorbed PAR per leaf area [umol photons/m**2/s]
    real(r8) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    real(r8) :: jpar              ! absorbed photons in photocenters as an electron
                                  ! rate (umol electrons/m**2/s)
    logical :: err

    ! Electron transport rate for C3 plants.
    ! Convert absorbed photon density [umol/m2 leaf /s] to
    ! that absorbed only by the photocenters (fnps) and also
    ! convert from photon energy into electron transport rate (photon_to_e)
    jpar = par_abs*photon_to_e*(1.0_r8 - fnps)
    
    ! convert the absorbed par into absorbed par per m2 of leaf,
    ! so it is consistant with the vcmax and lmr numbers.
    aquad = theta_psii
    bquad = -(jpar + jmax)
    cquad = jpar * jmax
    call QuadraticRoots(aquad, bquad, cquad, r1, r2, err)

    if(debug)then
       if(err)then
          write(fates_log(),*) "jequadfail:",par_abs,jpar,jmax
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    je = min(r1,r2)
    
  end function GetJe

  ! =====================================================================================
  
  function AgrossRuBPC3(par_abs,jmax,ci,co2_cpoint) result(aj)

    ! Input
    real(r8) :: par_abs    ! Absorbed PAR per leaf area [umol photons/m2leaf/s ]
    real(r8) :: jmax       ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: ci         ! intracellular leaf CO2 (Pa)
    real(r8) :: co2_cpoint ! CO2 compensation point (Pa)

    ! Output
    real(r8) :: aj         ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)

    ! locals
    real(r8) :: je         ! actual electron transport rate (umol electrons/m**2/s)
    
    ! Get the smoothed (quadratic between J and Jmax) electron transport rate
    je = GetJe(par_abs,jmax)

    
    aj = je * max(ci-co2_cpoint, 0._r8) / &
                 (4._r8*ci+8._r8*co2_cpoint)
    
    
  end function AgrossRuBPC3
  
  ! =======================================================================================

  function AgrossRuBPC4(par_abs) result(aj)

    real(r8) :: par_abs ! Absorbed PAR per leaf area [umol photons/m2leaf/s ]
    real(r8) :: aj      ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)

    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(r8),parameter :: c4_quant_eff = 0.05_r8
    
    aj = c4_quant_eff*par_abs
    
  end function AgrossRuBPC4

  ! =======================================================================================

  function AgrossPEPC4(ci,kp,can_press) result(ap)

    real(r8) :: ci                ! intracellular leaf CO2 (Pa)
    real(r8) :: kp                ! initial co2 response slope 
    real(r8) :: can_press         ! Air pressure near the surface of the leaf (Pa)
    real(r8) :: ap                ! PEP limited gross assimilation rate (umol co2/m2/s)     
    
    ap = kp * max(ci, 0._r8) / can_press
    
  end function AgrossPEPC4

  ! =======================================================================================


  subroutine CiMinMax(ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
       can_co2_ppress,can_o2_ppress,can_press,lmr,par_abs, &
       gb,gs0,ci_min,ci_max)

    ! This routine is used to find the first values of Ci that are the bounds
    ! for the bisection algorithm. It finds the values associated with minimum
    ! and maximum conductance when equating the source and sink limitations
    ! on net photosynthesis.

    ! input
    integer, intent(in)  :: ft             ! plant functional type index
    real(r8), intent(in) :: vcmax          ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in) :: jmax           ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: kp             ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in) :: co2_cpoint     ! CO2 compensation point (Pa)
    real(r8), intent(in) :: mm_kco2        ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(in) :: mm_ko2         ! Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(in) :: can_co2_ppress ! Partial pressure of CO2 near the leaf surface (Pa)
    real(r8), intent(in) :: can_o2_ppress  ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8), intent(in) :: can_press      ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in) :: lmr            ! leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in) :: par_abs        ! par absorbed per unit leaf area [umol photons/m2leaf/s ]
    real(r8), intent(in) :: gb             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: gs0            ! stomatal intercept
    
    ! output
    real(r8), intent(out) :: ci_max  ! Intracellular Co2 at maximum conductance [Pa]
    real(r8), intent(out) :: ci_min  ! Intracellular CO2 at minimum conductance [Pa]
    
    ! Intermediate terms
    real(r8), dimension(3) :: ci  ! Ci for Rubisco,RuBP and PEP respectively
    real(r8), dimension(3) :: ag  ! Agross for Rubisco, Rubp and PEP

    ! These are compound terms used to solve the equation that balances
    ! net assimilation with the gradient flux equation
    real(r8) :: a,b,c,d,e,f,g
    real(r8) :: ap,ac,aj,ai
    real(r8) :: Je                ! Electron tranport rate (umol e/m2/s)
    real(r8) :: rmin ,rmax        ! Maximum and minimum resistance [s/umol h2o/m2]
    real(r8) :: aquad,bquad,cquad,r1,r2
    logical :: err
    
    ! Minimum possible resistance (with a little buffer)
    rmin = 0.75_r8*h2o_co2_bl_diffuse_ratio/gb
    ! Maximum possible resistance (with a littel buffer)
    rmax = 1.25_r8*(h2o_co2_bl_diffuse_ratio/gb + h2o_co2_stoma_diffuse_ratio/gs0)
    
    if (lb_params%c3psn(ft) == c4_path_index)then

       ! Maximum conductance (minimum resistance)
       
       ag(1) = vcmax
       ag(2) = AgrossRuBPC4(par_abs)

       ! C4: Rubisco-limited photosynthesis
       ac = vcmax
       
       ! C4: RuBP-limited photosynthesis
       aj = AgrossRuBPC4(par_abs)
       
       aquad = theta_cj_c4
       bquad = -(ac + aj)
       cquad = ac * aj
       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c41quadfail_minmax1:",par_abs,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ai = min(r1,r2)

       a = rmin*can_press
       
       aquad = theta_ip_c4/a**2.0_r8 + kp/(can_press*a)
       bquad = - theta_ip_c4*2.0_r8*can_co2_ppress/(a*a) &
               - theta_ip_c4*2._r8*lmr/a &
               - kp*lmr/can_press &
               + ai/a &
               - kp*can_co2_ppress/(can_press*a) &
               + ai*kp/can_press
       cquad = theta_ip_c4*can_co2_ppress*can_co2_ppress/(a*a) &
             + 2._r8*lmr*can_co2_ppress*theta_ip_c4/a  &
             + theta_ip_c4*lmr*lmr &
             - ai*can_co2_ppress/a &
             - ai*lmr

       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c41quadfail_minmax2:",par_abs,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ci(3) = max(r1,r2)
       
       !! C4: PEP carboxylase-limited (CO2-limited)
       ap = AgrossPEPC4(ci(3),kp,can_press)

       aquad = theta_ip_c4
       bquad = -(ai + ap)
       cquad = ai * ap
       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c42quadfail:",par_abs,ci,kp,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ag(:) = min(r1,r2)

       if(debug) then
          if( abs(ci(3)-(can_co2_ppress - (ag(3)-lmr)*can_press*rmin))>0.001_r8) then
             write(fates_log(),*) "c4 ci conv check fail"
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       
       ci(1) = can_co2_ppress - (ag(1)-lmr)*can_press*rmin
       ci(2) = can_co2_ppress - (ag(1)-lmr)*can_press*rmin
       
       ci_max = ci(minloc(ag,DIM=1))

       ! Minimum conductance (maximum resistance)

       a = rmax*can_press
       
       aquad = theta_ip_c4/a**2.0_r8 + kp/(can_press*a)
       bquad = - theta_ip_c4*2.0_r8*can_co2_ppress/(a*a) &
               - theta_ip_c4*2._r8*lmr/a &
               - kp*lmr/can_press &
               + ai/a &
               - kp*can_co2_ppress/(can_press*a) &
               + ai*kp/can_press
       cquad = theta_ip_c4*can_co2_ppress*can_co2_ppress/(a*a) &
             + 2._r8*lmr*can_co2_ppress*theta_ip_c4/a  &
             + theta_ip_c4*lmr*lmr &
             - ai*can_co2_ppress/a &
             - ai*lmr

       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c41quadfail_minmax2:",par_abs,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ci(3) = max(r1,r2)

       ap = AgrossPEPC4(ci(3),kp,can_press)

       aquad = theta_ip_c4
       bquad = -(ai + ap)
       cquad = ai * ap
       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c42quadfail:",par_abs,ci,kp,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ag(:) = min(r1,r2)
       
       ci(1) = can_co2_ppress - (ag(1)-lmr)*can_press*rmax
       ci(2) = can_co2_ppress - (ag(1)-lmr)*can_press*rmax

       if(debug) then
          if( abs(ci(3)-(can_co2_ppress - (ag(3)-lmr)*can_press*rmax))>0.001_r8) then
             write(fates_log(),*) "c4 ci conv check fail"
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       ci_min = ci(minloc(ag,DIM=1))

       return
    end if
       
    ! Get the maximum e tranport rate for when we solve for RuBP (twice)
    Je = GetJe(par_abs,jmax)
    
    ! Find ci at maximum conductance (1/inf = 0)
    
    a = can_co2_ppress
    b = can_press*rmin
    c = vcmax
    d = vcmax*co2_cpoint
    e = 1._r8
    f = mm_kco2*(1._r8+can_o2_ppress / mm_ko2 )
    g = lmr
    ci(1) = CiFromAnetDiffGrad(a,b,c,d,e,f,g)
    ag(1) = AgrossRubiscoC3(vcmax,ci(1),can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2)
    
    
    ! check the math
    if(debug)then
       if ( abs((can_co2_ppress-ci(1))/b - (ag(1)-lmr)  ) > 1.e-3_r8 ) then
          write(fates_log(),*) 'incorrect ci, max cond:',ci(1),(can_co2_ppress-ci(1))/b,ag(1)-lmr
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    c = Je
    d = Je*co2_cpoint
    e = 4._r8
    f = 8._r8*co2_cpoint
    g = lmr
    ci(2) = CiFromAnetDiffGrad(a,b,c,d,e,f,g)
    ag(2) = AgrossRuBPC3(par_abs,jmax,ci(2),co2_cpoint)
    
    if(debug)then
       if ( abs((can_co2_ppress-ci(2))/b - (ag(2)-lmr))  > 1.e-3_r8 ) then
          write(fates_log(),*)'incorrect ci, max cond:',ci(2),(can_co2_ppress-ci(2))/b,ag(2)-lmr
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    ! Take the Ci that matches the minimizing gross assimilation
    ci_max = ci(minloc(ag(1:2),DIM=1))
    
    ! Find ci at minimum conductance (1/g0) (max conductance)

    a = can_co2_ppress
    b = can_press*rmax
    c = vcmax
    d = vcmax*co2_cpoint
    e = 1._r8
    f = mm_kco2*(1._r8+can_o2_ppress / mm_ko2 )
    g = lmr
    ci(1) = CiFromAnetDiffGrad(a,b,c,d,e,f,g) 
    ag(1) = AgrossRubiscoC3(vcmax,ci(1),can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2)
    
    ! check the math
    if(debug)then
       if ( abs((can_co2_ppress-ci(1))/b - (ag(1)-lmr)  ) > 1.e-3_r8 ) then
          write(fates_log(),*)'incorrect ci, min cond:',ci(1),(can_co2_ppress-ci(1))/b,ag(1)-lmr
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    c = Je
    d = Je*co2_cpoint
    e = 4._r8
    f = 8._r8*co2_cpoint
    g = lmr
    ci(2) = CiFromAnetDiffGrad(a,b,c,d,e,f,g)
    ag(2) = AgrossRuBPC3(par_abs,jmax,ci(2),co2_cpoint)
    
    if(debug)then
       if ( abs((can_co2_ppress-ci(2))/b -(ag(2)-lmr))  > 1.e-3_r8 ) then
          write(fates_log(),*)'incorrect ci, min cond:',ci(2),(can_co2_ppress-ci(2))/b,ag(2)-lmr
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    ci_min = ci(minloc(ag(1:2),DIM=1))
    
    return
  end subroutine CiMinMax

  ! =====================================================================================

  function CiFromAnetDiffGrad(a,b,c,d,e,f,g) result(ci)

    ! When the equation for net photosynthesis is coupled with the diffusion
    ! flux equation, for both Rubisco and RuBP limited assimilation
    ! the form is like so:
    !
    ! (a-ci)/b = (c*ci - d)/(e*ci + f) - g 
    !
    ! The expression below simply solves for ci
    ! This function is called to find endpoints, where conductance is maximized
    ! and minimized, to perform a binary search
    
    real(r8) :: a,b,c,d,e,f,g     ! compound terms to solve the coupled Anet
                                  ! and diffusive flux gradient equations
    real(r8) :: ci                ! intracellular co2 [Pa]
    real(r8) :: r1,r2             ! roots for quadratic
    real(r8) :: aquad,bquad,cquad
    logical  :: err

    aquad = -(e/b)
    bquad = (e*a-f)/b + e*g - c
    cquad = (f*a)/b + f*g + d
    call QuadraticRoots(aquad, bquad, cquad, r1, r2, err)

    ci = max(r1,r2)

  end function CiFromAnetDiffGrad
  
  ! =======================================================================================

  subroutine CiFunc(ci, &
       ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
       can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
       gs0, gs1, gs2, &
       anet,agross,gs,fval)

    ! -----------------------------------------------------------------------------------
    ! DESCRIPTION:
    ! The photosynthesis solver must solve 3 equations simultaneously,
    ! Anet, Stomatal Conductance and the update to Ci.  To solve this,
    ! we search values of Ci that satisfy the equations, and create an evaluation
    ! function "fval" that is the difference between the input value and the predicted
    ! value of ci.
    !
    ! fval = ci_input - ci_updated
    !
    ! The updated ci:
    ! ci =  ca - (1.37rb+1.65rs))*patm*anet
    !
    ! Thus:
    ! fval = ci_input - [ca-(1.37rb+1.65rs))*patm*anet]
    !
    ! This method borrows from Jinyun Tang's 2011 code for CLM
    ! -----------------------------------------------------------------------------------

    
    real(r8), intent(in)  :: ci             ! Input (trial) intracellular leaf CO2 (Pa)
    integer, intent(in)   :: ft             ! plant functional type index
    real(r8), intent(in)  :: vcmax          ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in)  :: jmax           ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)  :: kp             ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in)  :: mm_kco2        ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(in)  :: mm_ko2         ! Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(in)  :: co2_cpoint     ! CO2 compensation point (Pa)
    real(r8), intent(in)  :: can_press      ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in)  :: can_co2_ppress ! Partial pressure of CO2 near the leaf surface (Pa)
    real(r8), intent(in)  :: can_o2_ppress  ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8), intent(in)  :: can_vpress     ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in)  :: lmr            ! leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)  :: par_abs        ! par absorbed per unit lai [umol photons/m2leaf/s ]
    real(r8), intent(in)  :: gb             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in)  :: veg_tempk      ! vegetation temperature
    real(r8), intent(in)  :: veg_esat       !
    real(r8), intent(in)  :: gs0            ! conductance intercept (umol H20/m2/s)
    real(r8), intent(in)  :: gs1            ! conductance slope (could be multiplied by btran)
    real(r8), intent(in)  :: gs2            ! for Medlyn only: either 1 or btran
    real(r8), intent(out) :: anet           ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: agross         ! gross leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: gs             ! stomatal conductance (umol h2o/m2/s)
    real(r8), intent(out) :: fval           ! ci_input - ci_updated  (Pa)
    
    !real(r8) :: veg_esat          ! Saturation vapor pressure at leaf-surface [Pa]
    real(r8) :: veg_qs            ! DUMMY, specific humidity at leaf-surface [kg/kg]
    real(r8) :: a_gs              ! The assimilation (a) for calculating conductance (gs)
                                  ! is either = to anet or agross
    real(r8) :: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: ap                ! product-limited (C3) or CO2-limited
                                  ! (C4) gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: ai                ! For C4, smoothed co-limited assimilation of Rubisco/RuBP
    real(r8) :: aquad,bquad,cquad ! a,b,c terms in the quadratic equation
    real(r8) :: r1,r2             ! The roots from the quadratic
    logical  :: err
    !------------------------------------------------------------------------------

    ! Photosynthesis limitation rate calculations
    if (lb_params%c3psn(ft) == c3_path_index)then
       
       ! C3: Rubisco-limited photosynthesis
       ac = AgrossRubiscoC3(vcmax,ci,can_o2_ppress,co2_cpoint,mm_kco2,mm_ko2)
       
       ! C3: RuBP-limited photosynthesis
       aj = AgrossRuBPC3(par_abs,jmax,ci,co2_cpoint )

       ! Take the minimum, no smoothing
       agross = min(ac,aj)
       
    else
       
       ! C4: Rubisco-limited photosynthesis
       ac = vcmax

       ! C4: RuBP-limited photosynthesis
       aj = AgrossRuBPC4(par_abs)
                 
       ! C4: PEP carboxylase-limited (CO2-limited)
       ap = AgrossPEPC4(ci,kp,can_press)

       aquad = theta_cj_c4
       bquad = -(ac + aj)
       cquad = ac * aj
       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c41quadfail:",par_abs,ci,kp,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if

       ai = min(r1,r2)

       aquad = theta_ip_c4
       bquad = -(ai + ap)
       cquad = ai * ap
       call QuadraticRoots(aquad, bquad, cquad, r1, r2,err)

       if(debug)then
          if(err)then
             write(fates_log(),*) "c42quadfail:",par_abs,ci,kp,can_press,vcmax
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if
       
       agross = min(r1,r2)
       
    end if

    ! Calculate anet
    anet = agross  - lmr
    
    if (  lb_params%stomatal_assim_model == gross_assim_model ) then
       if ( lb_params%stomatal_model == medlyn_model ) then
          write (fates_log(),*) 'Gross Assimilation conductance is incompatible with the Medlyn model'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       a_gs = agross
    else
       a_gs = anet
    end if

    ! A note about the use of the quadratic equations for calculating stomatal conductance
    ! ------------------------------------------------------------------------------------
    ! These two following models calculate the conductance between the intracellular leaf
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
    ! intracellular humidity (e_i, which is the saturation humidity at leaf temperature),
    ! boundary layer conductance (g_b) (these are all known) and stomatal conductance
    ! (g_s) (this is still unknown).  This expression is substituted into the stomatal
    ! conductance equation. The resulting form of these equations becomes a quadratic.
    !
    ! For a detailed explanation, see the FATES technical note, section
    ! "1.11 Stomatal Conductance"
    !
    ! ------------------------------------------------------------------------------------

    leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb * a_gs * can_press
    leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
    
    ! Determine saturation vapor pressure at the leaf surface, from temp and atm-pressure
    !call QSat(veg_tempk, can_press, veg_qs, veg_esat)
    
    if ( lb_params%stomatal_model == medlyn_model ) then
       call StomatalCondMedlyn(anet,veg_esat,can_vpress,gs0,gs1,gs2, &
            leaf_co2_ppress,can_press,gb,gs)
    else
       call StomatalCondBallBerry(a_gs,anet,veg_esat,can_vpress,gs0,gs1, &
            leaf_co2_ppress,can_press,gb,gs)
    end if
    
    ! Derive new estimate for ci
    ! ci = can_co2_ppress - anet * can_press * &
    ! (h2o_co2_bl_diffuse_ratio/gb + h2o_co2_stoma_diffuse_ratio/gs_out)
    
    ! fval = ci_input - ci_predicted
    fval = ci - (can_co2_ppress - anet * can_press * &
         (h2o_co2_bl_diffuse_ratio*gs + h2o_co2_stoma_diffuse_ratio*gb)/(gb*gs))

  end subroutine CiFunc

  ! ====================================================================================
  
  subroutine CiBisection(ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
       can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
       gs0,gs1,gs2,ci_tol, &
       anet,agross,gs,ci,solve_iter)

    ! -----------------------------------------------------------------------------------
    !
    ! Co-solve for photosynthesis and stomatal conductance via a bisectional search for
    ! intracellular CO2. This is an inefficient yet robust search method that we resort
    ! to if the recursive solver cannot find an answer within its convergence criteria.
    ! -----------------------------------------------------------------------------------
    
    integer, intent(in)    :: ft             ! plant functional type index
    real(r8), intent(in)   :: vcmax          ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in)   :: jmax           ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)   :: kp             ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in)   :: mm_kco2        ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(in)   :: mm_ko2         ! Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(in)   :: co2_cpoint     ! CO2 compensation point (Pa)
    real(r8), intent(in)   :: can_press      ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in)   :: can_co2_ppress ! Partial pressure of CO2 near the leaf surface (Pa)
    real(r8), intent(in)   :: can_o2_ppress  ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8), intent(in)   :: can_vpress     ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in)   :: lmr            ! leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)   :: par_abs        ! par absorbed per unit lai [umol photons/m2leaf/s ]
    real(r8), intent(in)   :: gb             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in)   :: veg_tempk      ! vegetation temperature [Kelvin]
    real(r8), intent(in)   :: veg_esat
    real(r8), intent(in)   :: gs0,gs1,gs2    ! stomatal intercept and slope and alternative btran
    real(r8), intent(in)   :: ci_tol         ! Convergence tolerance for solutions for intracellular CO2 (Pa)
    real(r8), intent(out)  :: anet           ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out)  :: agross         ! gross leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out)  :: gs             ! stomatal conductance (umol h2o/m2/s)
    real(r8), intent(out)  :: ci             ! Input (trial) intracellular leaf CO2 (Pa)
    integer, intent(inout) :: solve_iter     ! number of bisections required

    ! With bisection, we need to keep track of three different ci values at any given time
    ! The high, the low and the bisection.

    real(r8) :: ci_h, fval_h
    real(r8) :: ci_l, fval_l
    real(r8) :: ci_b, fval_b

    logical :: loop_continue   ! Continue bisecting until tolerance reached

    ! Maximum number of iterations on intracelluar co2 solver until is quits
    integer, parameter :: max_iters = 200

    ! Find the starting points (end-points) for bisection
    ! We dont need stomatal slope because we just want the two extremes
    ! which is the intercept and infinite conductance
    call CiMinMax(ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
         can_co2_ppress,can_o2_ppress,can_press,lmr,par_abs, &
         gb,gs0,ci_l,ci_h)
    
    call CiFunc(ci_h, &
         ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
         can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
         gs0,gs1,gs2, &
         anet,agross,gs,fval_h)
    
    call CiFunc(ci_l, &
         ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
         can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
         gs0,gs1,gs2, &
         anet,agross,gs,fval_l)

    ! It is necessary that our starting points are on opposite sides of the root
    if( nint(fval_h/abs(fval_h)) .eq. nint(fval_l/abs(fval_l)) ) then
       
       ! Try an exteremly large bisection range, if this doesn't work, then
       ! fail the run
       ci_h = 0.000001_r8
       call CiFunc(ci_h, &
            ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
            can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
            gs0,gs1,gs2, &
            anet,agross,gs,fval_h)

       ci_l = 2000._r8
       call CiFunc(ci_l, &
         ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
         can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk, veg_esat, &
         gs0,gs1,gs2, &
         anet,agross,gs,fval_l)
       
       ! It is necessary that our starting points are on opposite sides of the root
       if( nint(fval_h/abs(fval_h)) .eq. nint(fval_l/abs(fval_l)) ) then
          write(fates_log(),*)'While attempting bisection for Ci calculations,'
          write(fates_log(),*)'the two starting values for Ci were on the same'
          write(fates_log(),*)'side of the root. Try increasing and decreasing'
          write(fates_log(),*)'init_ci_high and init_ci_low respectively'
          write(fates_log(),*) "ci_h=",ci_h,"fval_h=",fval_h,"ci_l=",ci_l,"fval_l=",fval_l
          write(fates_log(),*) "ft= ",ft,"is c3psn:",lb_params%c3psn(ft) == c3_path_index
          write(fates_log(),*) "vcmax=",vcmax,"jmax=",jmax,"kp=",kp
          write(fates_log(),*) "co2_cpoint=",co2_cpoint,"mm_kco2=",mm_kco2,"mm_ko2=",mm_ko2
          write(fates_log(),*) "can_co2_ppress=",can_co2_ppress,"can_o2_ppress=",can_o2_ppress,"can_press=",can_press
          write(fates_log(),*) "can_vpress=",can_vpress,"lmr=",lmr,"par_abs=",par_abs,"gb=",gb
          write(fates_log(),*) "veg_tempk=",veg_tempk,"gs0=",gs0,"gs1=",gs1,"gs2=",gs2,"ci_tol=",ci_tol
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    loop_continue = .true.
    ci_b = 0.5*(ci_l+ci_h)
    bi_iter_loop: do while(loop_continue)

       solve_iter = solve_iter + 1

       call CiFunc(ci_b, &
            ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
            can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
            gs0,gs1,gs2, &
            anet,agross,gs,fval_b)

       if(abs(fval_b)<ci_tol) then
          !loop_continue = .false.
          exit bi_iter_loop
       end if
       
       if( solve_iter == max_iters) then
          write (fates_log(),*) 'Ci bisection during photosynthesis failed'
          write (fates_log(),*) 'try increasing tolerance or widening the starting points'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       
       ! Determining which side of the interval we want to retain
       if( nint(fval_b/abs(fval_b)) .eq. nint(fval_h/abs(fval_h)) ) then
          fval_h = fval_b
          ci_h   = ci_b
       else
          fval_l = fval_b
          ci_l   = ci_b
       end if

       ci_b = 0.5*(ci_l+ci_h)

    end do bi_iter_loop
    
    ! The Ci value that is constent with the last call to CiFunc()
    ! is the bisection point minus the difference
    
    ci = ci_b - fval_b

    return
  end subroutine CiBisection
  
  ! =====================================================================================
  
  subroutine LeafLayerPhotosynthesis( &
       par_abs,           &  ! in
       ft,                &  ! in
       vcmax,             &  ! in
       jmax,              &  ! in
       kp,                &  ! in
       gs0,               &  ! in
       gs1,               &  ! in
       gs2,               &  ! in
       veg_tempk,         &  ! in
       can_press,         &  ! in
       can_co2_ppress,    &  ! in
       can_o2_ppress,     &  ! in
       veg_esat,          &  ! in
       gb,                &  ! in
       can_vpress,        &  ! in
       mm_kco2,           &  ! in
       mm_ko2,            &  ! in
       co2_cpoint,        &  ! in
       lmr,               &  ! in
       ci_tol,            &  ! in
       agross,            &  ! out
       gs,                &  ! out
       anet,              &  ! out
       c13disc,           &  ! out
       ci,                &  ! out
       solve_iter)           ! out


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
    real(r8), intent(in) :: par_abs           ! Absorbed PAR per leaf area [umol photons/m2 leaf/s]
    integer,  intent(in) :: ft                ! (plant) Functional Type Index
    real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: kp                ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in) :: gs0               ! effective stomatal conductance intercept
    real(r8), intent(in) :: gs1               ! effective stomatal conductance slope, applies to 
    real(r8), intent(in) :: gs2               ! for Medlyn conductance only, alternative btran
                                              ! term that applies to non-intercept side of equation
    real(r8), intent(in) :: veg_tempk         ! vegetation temperature
    real(r8), intent(in) :: can_press         ! Air pressure near the surface of the leaf (Pa)
    real(r8), intent(in) :: can_co2_ppress    ! Partial pressure of CO2 near the leaf surface (Pa)
    real(r8), intent(in) :: can_o2_ppress     ! Partial pressure of O2 near the leaf surface (Pa)
    real(r8), intent(in) :: veg_esat          ! saturation vapor pressure at vegetation
    real(r8), intent(in) :: gb                ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: can_vpress        ! vapor pressure of the canopy air (Pa)
    real(r8), intent(in) :: mm_kco2           ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(in) :: mm_ko2            ! Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(in) :: co2_cpoint        ! CO2 compensation point (Pa)
    real(r8), intent(in) :: lmr               ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
    real(r8), intent(in) :: ci_tol            ! Convergence tolerance for solutions for intracellular CO2 (Pa)
    real(r8), intent(out) :: agross           ! gross photosynthesis (umolC/m2/s)
    real(r8), intent(out) :: gs               ! leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out) :: anet             ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8), intent(out) :: c13disc          ! carbon 13 in newly assimilated carbon
    real(r8), intent(out) :: ci               ! intracellular leaf CO2 (Pa)
    integer,  intent(out) :: solve_iter       ! Number of iterations required for the solve
    
    ! Important Note on the gas pressures as input arguments.  This photosynthesis scheme will iteratively
    ! solve for the co2 partial pressure at the leaf surface (ie in the stomata). The reference
    ! point for these input values are NOT within that boundary layer that separates the stomata from
    ! the canopy air space.  The reference point for these is on the outside of that boundary
    ! layer.  This routine, which operates at the leaf scale, makes no assumptions about what the
    ! scale of the reference is, it could be lower atmosphere, it could be within the canopy
    ! but most likely it is the closest value one can get to the edge of the leaf's boundary
    ! layer.  We use the convention "can_" because a reference point of within the canopy
    ! ia a best reasonable scenario of where we can get that information from.
    
    ! Locals
    ! ------------------------------------------------------------------------
    real(r8) :: ci0               ! Local perturbation values of intracellular co2 [Pa]
    real(r8) :: fval              ! Change in intracellular co2 through
                                  ! one iteration of a solve (ci_init - ci) [Pa]
    logical  :: loop_continue     ! Loop control variable
      
    ! Parameters
    ! ------------------------------------------------------------------------
    

    ! First guess on ratio between intracellular co2 and the atmosphere
    ! an iterator converges on actual
    real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
    real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

    ! For testing, it is useful to force the bisection method
    logical, parameter :: force_bisection = .false.

    ! Maximum number of iterations on intracelluar co2 solver until is quits
    integer, parameter :: max_iters = 10

    ! Assume a trival solve until we encounter both leaf and light
    solve_iter = 0

    ! Less, but still trivial solution - biomass, but no light, no photosynthesis
    ! Stomatal conductance is the intercept of the conductance functions
    ! ---------------------------------------------------------------------------------------------
    if (par_abs < nearzero ) then
       anet    = -lmr
       agross  = 0._r8
       gs      = gs0
       c13disc = 0.0_r8
       return
    end if
    
    ! Not trivial solution, some biomass and some light
    ! Initialize first guess of intracellular co2 conc [Pa]
    if (lb_params%c3psn(ft) == c3_path_index) then
       ci0 = init_a2l_co2_c3 * can_co2_ppress
    else
       ci0 = init_a2l_co2_c4 * can_co2_ppress
    end if

    loop_continue = .true.
    iter_loop: do while(loop_continue)
       
       ! Increment iteration counter. Stop if too many iterations
       solve_iter = solve_iter + 1

       ! This function calculates the difference between
       ! the guessed ci (ie ci0) and an updated value and returns
       ! anet, agross, gs and the difference in ci (fval)
       ! We encapsulate this process in the CiFunc() because
       ! this is the objective funtion that we want the root of
       ! in various search methods such as NR, Secant and bisection
       if(.not.force_bisection)then
          call CiFunc(ci0, &
               ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
               can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
               gs0,gs1,gs2, &
               anet,agross,gs,fval)
          
          ! fval = ci_input - ci_predicted
          ! ci_predicted = ci_input - fval
          ci = ci0 - fval
          
          ! In main, ci_tol = 2*can_press
          ! Special convergence requirement to satisfy B4B with main
          
          if (abs(fval) <= ci_tol ) then
             loop_continue = .false.
             exit iter_loop
          end if

       end if
       
       ci0 = ci
       
       if( solve_iter == max_iters .or. force_bisection) then
          call CiBisection( &
               ft,vcmax,jmax,kp,co2_cpoint,mm_kco2,mm_ko2, &
               can_co2_ppress,can_o2_ppress,can_press,can_vpress,lmr,par_abs,gb,veg_tempk,veg_esat, &
               gs0,gs1,gs2,ci_tol, &
               anet,agross,gs,ci,solve_iter)
          loop_continue = .false.
          exit iter_loop
       end if
       
    end do iter_loop

    
    ! $\Delta ^{13} C = \alpha_s + (b - \alpha_s) \cdot \frac{C_i}{C_a}$
    ! just hard code b and \alpha_s for now, might move to parameter set in future
    ! b = 27.0 alpha_s = 4.4
    ! TODO, not considering C4 or CAM right now, may need to address this

    c13disc = 4.4_r8 + (27.0_r8 - 4.4_r8) * &
         min (can_co2_ppress, max (ci, 0._r8)) / can_co2_ppress

    return
  end subroutine LeafLayerPhotosynthesis

  ! =======================================================================================

  function LeafHumidityStomaResis(leaf_psi, k_lwp, veg_tempk, can_vpress, can_press, &
       rb, gstoma, ft, veg_esat) result(rstoma_out)

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
    real(r8) :: veg_esat   ! Saturation vapor pressure at veg surface [Pa]
    
    ! Locals
    real(r8) :: ceair      ! vapor pressure of air, constrained [Pa]
                           ! water potential to mesophyll water potential
    real(r8) :: qs         ! Specific humidity [g/kg]
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
    
    ! call QSat(veg_tempk, can_press, qsat_alt, veg_esat)

    qsat_alt = qsat_alt * g_per_kg
    
    ceair = GetConstrainedVPress(can_vpress,veg_esat)
    
    ! compute specific humidity from vapor pressure
    ! q = molar_mass_ratio_vapdry*e/(can_press - (1-molar_mass_ratio_vapdry)*e) 
    ! source https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    ! now adjust inner leaf humidity by LWP_star

    qs = molar_mass_ratio_vapdry * ceair / (can_press - (1._r8-molar_mass_ratio_vapdry) * ceair)
    qsat_loc = molar_mass_ratio_vapdry * veg_esat / (can_press - (1._r8-molar_mass_ratio_vapdry) * veg_esat)
    qsat_adj = qsat_loc*lwp_star

    if(debug)then
       if (abs(qsat_loc-qsat_alt) > 1.e-2) then
          write (fates_log(),*) 'qsat from QSat():', qsat_alt
          write (fates_log(),*) 'qsat localy :', qsat_loc
          write (fates_log(),*) 'values of qsat are too different'
          call endrun(msg=errMsg(sourcefile, __LINE__))  
       end if
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
       kp25_ft, &
       nscaler,    &
       veg_tempk,      &
       dayl_factor, &
       t_growth,   &
       t_home,     &
       btran, &
       vcmax, &
       jmax, &
       kp, &
       gs0, &
       gs1, &
       gs2)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the localized values of several key photosynthesis
    ! rates.  By localized, we mean specific to the plant type and leaf layer,
    ! which factors in leaf physiology, as well as environmental effects.
    ! This procedure should be called prior to iterative solvers, and should
    ! have pre-calculated the reference rates for the pfts before this.
    !
    ! The output biophysical rates are:
    ! vcmax: maximum rate of carboxilation,
    ! jmax: maximum electron transport rate,
    ! kp: initial slope of CO2 response curve (C4 plants)
    ! gs0,gs1,gs2: Stomatal conductance slopes and intercepts (multiplied by btran)
    ! ---------------------------------------------------------------------------------

    ! Arguments
    ! ------------------------------------------------------------------------------

    integer,  intent(in) :: ft                        ! (plant) Functional Type Index
    real(r8), intent(in) :: nscaler                   ! Scale for leaf nitrogen profile
    real(r8), intent(in) :: vcmax25top_ft             ! canopy top maximum rate of carboxylation at 25C
                                                      ! for this pft (umol CO2/m**2/s)
    real(r8), intent(in) :: jmax25top_ft              ! canopy top maximum electron transport rate at 25C
                                                      ! for this pft (umol electrons/m**2/s)
    real(r8), intent(in) :: kp25_ft                   ! initial slope of CO2 response curve
                                                      ! (C4 plants) at 25C, canopy top, this pft
    real(r8), intent(in) :: veg_tempk                 ! vegetation temperature
    real(r8), intent(in) :: dayl_factor               ! daylength scaling factor (0-1)
    real(r8), intent(in) :: t_growth                  ! T_growth (short-term running mean temperature) (K)
    real(r8), intent(in) :: t_home                    ! T_home (long-term running mean temperature) (K)
    real(r8), intent(in) :: btran                     ! transpiration wetness factor (0 to 1)
    real(r8), intent(out) :: vcmax                    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(out) :: jmax                     ! maximum electron transport rate
                                                      ! (umol electrons/m**2/s)
    real(r8), intent(out) :: kp                       ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(out) :: gs0                      ! effective stomatal intercept
    real(r8), intent(out) :: gs1                      ! effective stomatal slope
    real(r8), intent(out) :: gs2                      ! alternative btran term applied to Medlyn
                                                      ! conductance on whole non-intercept side of equation
    
    ! Locals
    ! -------------------------------------------------------------------------------
    real(r8) :: vcmax25           ! leaf layer: maximum rate of carboxylation at 25C
                                  ! (umol CO2/m**2/s)
    real(r8) :: jmax25            ! leaf layer: maximum electron transport rate at 25C
                                  ! (umol electrons/m**2/s)
    real(r8) :: dayl_factor_local ! Local version of daylength factor

    ! Parameters
    ! ---------------------------------------------------------------------------------
    real(r8) :: vcmaxha          ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha           ! activation energy for jmax (J/mol)
    real(r8) :: vcmaxhd          ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd           ! deactivation energy for jmax (J/mol)
    real(r8) :: vcmaxse          ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse           ! entropy term for jmax (J/mol/K)
    real(r8) :: t_growth_celsius ! average growing temperature
    real(r8) :: t_home_celsius   ! average home temperature
    real(r8) :: jvr              ! ratio of Jmax25 / Vcmax25
    real(r8) :: vcmaxc           ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc            ! scaling factor for high temperature inhibition (25 C = 1.0)


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
    
    ! Adjust for temperature
    ! photosynthetic pathway: 0. = c4, 1. = c3
    
    if (lb_params%c3psn(ft) == c3_path_index) then
       vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
       kp = -9999._r8
    else
       vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
       vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
       vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
       kp = kp25_ft * nscaler * 2._r8**((min(veg_tempk,310._r8)-(tfrz+25._r8))/10._r8)
    end if

    jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)
 
    ! Adjust various rates for water limitations
    ! -----------------------------------------------------------------------------------

    ! Apply water limitations to Vcmax
    if(lb_params%agross_btran_model(ft) .ne. btran_on_ag_none) then
       vcmax = vcmax * btran
    end if

    ! Apply water limitations to Jmax
    if(lb_params%agross_btran_model(ft) .eq. btran_on_ag_vcmax_jmax) then
       jmax = jmax * btran
    end if

    ! Make sure that vcmax and jmax do not drop below a lower
    ! threshold, see where the constants are defined for an explanation.
    ! -----------------------------------------------------------------------------------
    if( do_mincap_vcjmax ) then
       
       vcmax = max(min_vcmax_frac*vcmax25top_ft,vcmax)

       jmax = max(min_jmax_frac*jmax25top_ft,jmax)
       
    end if

    ! Apply water limitations to stomatal intercept (hypothesis dependent)

    if(lb_params%stomatal_btran_model(ft)==btran_on_gs_gs0  .or. &
       lb_params%stomatal_btran_model(ft)==btran_on_gs_gs01 .or. &
       lb_params%stomatal_btran_model(ft)==btran_on_gs_gs02 )then

       gs0 = max(gsmin0,lb_params%stomatal_intercept(ft)*btran)
       
    else
       gs0 = max(gsmin0,lb_params%stomatal_intercept(ft))
    end if

    ! Apply water limitations to stomatal slope (hypothesis dependent)
    gs2 = 1._r8 
    if(lb_params%stomatal_btran_model(ft)==btran_on_gs_gs1 .or. &
       lb_params%stomatal_btran_model(ft)==btran_on_gs_gs01)then
       if(lb_params%stomatal_model.eq.medlyn_model ) then
          gs1 = lb_params%medlyn_slope(ft)*btran
       else
          gs1 = lb_params%bb_slope(ft)*btran
       end if
    elseif(lb_params%stomatal_btran_model(ft)==btran_on_gs_gs2 .or. &
           lb_params%stomatal_btran_model(ft)==btran_on_gs_gs02)then
       if(lb_params%stomatal_model.eq.medlyn_model ) then
          gs2 = btran
          gs1 = lb_params%medlyn_slope(ft)
       else
          gs1 = lb_params%bb_slope(ft)
       end if
    else
       if(lb_params%stomatal_model.eq.medlyn_model ) then
          gs1 = lb_params%medlyn_slope(ft)
       else
          gs1 = lb_params%bb_slope(ft)
       end if
    end if

    

    
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

  subroutine LowstorageMainRespReduction(frac, ft, maintresp_reduction_factor)

    ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
    ! of this is that mortality of plants increases when storage is low because they are not able
    ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
    ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1,
    ! and the curvature of the function is determined by a parameter.

    ! Arguments
    ! ------------------------------------------------------------------------------
    real(r8), intent(in) :: frac                        ! ratio of storage to target leaf biomass
    integer,  intent(in) :: ft                          ! what pft is this cohort?
    real(r8), intent(out) :: maintresp_reduction_factor ! the factor by which to reduce maintenance respiration

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
       if ( abs(lb_params%maintresp_reduction_curvature(ft)-1._r8) > nearzero ) then
          maintresp_reduction_factor = (1._r8 - lb_params%maintresp_reduction_intercept(ft)) + &
               lb_params%maintresp_reduction_intercept(ft) * &
               (1._r8 - lb_params%maintresp_reduction_curvature(ft)**frac) &
               / (1._r8-lb_params%maintresp_reduction_curvature(ft))
       else  ! avoid nan answer for linear case
          maintresp_reduction_factor = (1._r8 - lb_params%maintresp_reduction_intercept(ft)) + &
               lb_params%maintresp_reduction_intercept(ft) * frac
       endif

    else
       maintresp_reduction_factor = 1._r8
    endif


  end subroutine LowstorageMainRespReduction

  ! =====================================================================================

  real(r8) function GetConstrainedVPress(air_vpress,veg_esat) result(ceair)

    ! -----------------------------------------------------------------------------------
    ! Return a constrained vapor pressure [Pa]
    !
    ! Vapor pressure may not be greater than saturation vapor pressure
    ! and may not be less than 1% of saturation vapor pressure
    !
    ! This function is used to make sure that Ball-Berry conductance is well
    ! behaved.
    ! -----------------------------------------------------------------------------------
    
    real(r8) :: air_vpress              ! vapor pressure of the air (unconstrained) [Pa]
    real(r8) :: veg_esat                ! saturated vapor pressure [Pa]
    real(r8) :: min_frac_esat = 0.05_r8 ! We don't allow vapor pressures
                                        ! below this fraction amount of saturation vapor pressure

    ceair = min( max(air_vpress, min_frac_esat*veg_esat ),veg_esat )
    
  end function GetConstrainedVPress

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
    ! photosynthesis module uses units of micromoles, while the land-energy
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
