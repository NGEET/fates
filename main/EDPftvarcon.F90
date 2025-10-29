module EDPftvarcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:

  use FatesRadiationMemMod, only: num_swb,ivis,inir
  use FatesRadiationMemMod, only: norman_solver,twostr_solver
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : itrue, ifalse
  use PRTParametersMod, only : prt_params
  use LeafBiophysicsMod, only : lb_params
  use LeafBiophysicsMod, only : btran_on_gs_gs02,btran_on_ag_vcmax_jmax
  use LeafBiophysicsMod, only : lmr_r_1, lmr_r_2
  use LeafBiophysicsMod, only : c3_path_index, c4_path_index
  use FatesGlobals,   only : fates_log
  use FatesGlobals,   only : endrun => fates_endrun
  use FatesLitterMod, only : ilabile,icellulose,ilignin
  use PRTGenericMod,  only : leaf_organ, fnrt_organ, store_organ
  use PRTGenericMod,  only : sapw_organ, struct_organ, repro_organ
  use PRTGenericMod,  only : prt_cnp_flex_allom_hyp,prt_carbon_allom_hyp
  use FatesInterfaceTypesMod, only : hlm_nitrogen_spec, hlm_phosphorus_spec
  use FatesInterfaceTypesMod, only : hlm_parteh_mode
  use FatesInterfaceTypesMod, only : hlm_nu_com
  use FatesConstantsMod   , only : ievergreen
  use FatesConstantsMod   , only : prescribed_p_uptake
  use FatesConstantsMod   , only : prescribed_n_uptake
  use FatesConstantsMod   , only : coupled_p_uptake
  use FatesConstantsMod   , only : coupled_n_uptake
  use FatesConstantsMod   , only : default_regeneration
  use FatesConstantsMod   , only : TRS_regeneration
  use FatesConstantsMod   , only : TRS_no_seedling_dyn

   ! CIME Globals
  use shr_log_mod ,   only : errMsg => shr_log_errMsg

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  integer, parameter, public :: lower_bound_pft = 1
  integer, parameter, public :: lower_bound_general = 1

  !ED specific variables.
  type, public ::  EDPftvarcon_type

     real(r8), allocatable :: freezetol(:)           ! minimum temperature tolerance
     real(r8), allocatable :: hgt_min(:)             ! sapling height m
     real(r8), allocatable :: dleaf(:)               ! leaf characteristic dimension length (m)
     real(r8), allocatable :: z0mr(:)                ! ratio of roughness length of vegetation to height (-)
     real(r8), allocatable :: displar(:)             ! ratio of displacement height to canopy top height
     real(r8), allocatable :: bark_scaler(:)         ! scaler from dbh to bark thickness. For fire model.
     real(r8), allocatable :: crown_kill(:)          ! scaler on fire death. For fire model.
     real(r8), allocatable :: initd(:)               ! initial seedling density

     real(r8), allocatable :: seed_suppl(:)          ! seeds that come from outside the gridbox.

     real(r8), allocatable :: lf_flab(:)             ! Leaf litter labile fraction [-]
     real(r8), allocatable :: lf_fcel(:)             ! Leaf litter cellulose fraction [-]
     real(r8), allocatable :: lf_flig(:)             ! Leaf litter lignin fraction [-]
     real(r8), allocatable :: fr_flab(:)             ! Fine-root litter labile fraction [-]
     real(r8), allocatable :: fr_fcel(:)             ! Fine-root litter cellulose fraction [-]
     real(r8), allocatable :: fr_flig(:)             ! Fine-root litter lignatn fraction [-]
     real(r8), allocatable :: xl(:)                  ! Leaf-stem orientation index
     real(r8), allocatable :: clumping_index(:)      ! factor describing how much self-occlusion
                                                     ! of leaf scattering elements
                                                     ! decreases light interception
     real(r8), allocatable :: smpso(:)               ! Soil water potential at full stomatal opening
                                                     ! (non-HYDRO mode only) [mm]
     real(r8), allocatable :: smpsc(:)               ! Soil water potential at full stomatal closure
                                                     ! (non-HYDRO mode only) [mm]

     real(r8), allocatable :: maintresp_leaf_vert_scaler_coeff1(:) ! leaf maintenance respiration decrease through the canopy param 1
                                                                   ! only with Atkin et al. 2017 respiration model
     real(r8), allocatable :: maintresp_leaf_vert_scaler_coeff2(:) ! leaf maintenance respiration decrease through the canopy param 2
                                                                   ! only with Atkin et al. 2017 respiraiton model 
     real(r8), allocatable :: bmort(:)
     real(r8), allocatable :: mort_ip_size_senescence(:) ! inflection point of dbh dependent senescence
     real(r8), allocatable :: mort_r_size_senescence(:)  ! rate of change in mortality with dbh
     real(r8), allocatable :: mort_ip_age_senescence(:)  ! inflection point of age dependent senescence
     real(r8), allocatable :: mort_r_age_senescence(:)   ! rate of change in mortality with age
     real(r8), allocatable :: mort_scalar_coldstress(:)  ! maximum mortality rate from cold stress
     real(r8), allocatable :: mort_scalar_cstarvation(:) ! maximum mortality rate from carbon starvation
     real(r8), allocatable :: mort_scalar_hydrfailure(:) ! maximum mortality rate from hydraulic failure
     real(r8), allocatable :: mort_upthresh_cstarvation(:) ! threshold for storage biomass (relative to target leaf biomass) above which carbon starvation is zero
     real(r8), allocatable :: hf_sm_threshold(:)         ! soil moisture (btran units) at which drought mortality begins for non-hydraulic model
     real(r8), allocatable :: hf_flc_threshold(:)        ! plant fractional loss of conductivity at which drought mortality begins for hydraulic model

     real(r8), allocatable :: germination_rate(:)        ! Fraction of seed mass germinating per year (yr-1)
     real(r8), allocatable :: seed_decay_rate(:)         ! Fraction of seed mass (both germinated and
                                                         ! ungerminated), decaying per year    (yr-1)
     real(r8), allocatable :: seed_dispersal_pdf_scale(:)  ! Seed dispersal scale parameter, Bullock et al. (2017)
     real(r8), allocatable :: seed_dispersal_pdf_shape(:)  ! Seed dispersal shape parameter, Bullock et al. (2017)
     real(r8), allocatable :: seed_dispersal_max_dist(:) ! Maximum seed dispersal distance parameter (m)
     real(r8), allocatable :: seed_dispersal_fraction(:) ! Fraction of seed rain to disperse, per pft

     real(r8), allocatable :: repro_frac_seed(:)         ! fraciton of reproductive carbon that is seed
     real(r8), allocatable :: a_emerg(:)                 ! mean fraction of seed bank emerging [day-1]
     real(r8), allocatable :: b_emerg(:)                 ! seedling emergence sensitivity to soil moisture
     real(r8), allocatable :: par_crit_germ(:)           ! critical light level for germination [MJ m2-1 day-1]
     real(r8), allocatable :: seedling_psi_emerg(:)      ! critical soil moisture for seedling emergence [mm h2o suction]
     real(r8), allocatable :: seedling_psi_crit(:)       ! critical soil moisture initiating seedling stress
     real(r8), allocatable :: seedling_light_rec_a(:)    ! coefficient in light-based seedling to sapling transition rate
     real(r8), allocatable :: seedling_light_rec_b(:)    ! coefficient in light-based seedling to sapling transition rate
     real(r8), allocatable :: seedling_mdd_crit(:)       ! critical moisture deficit day accumulation for seedling moisture-based
                                                         ! seedling mortality to begin
     real(r8), allocatable :: seedling_h2o_mort_a(:)     ! coefficient in moisture-based seedling mortality
     real(r8), allocatable :: seedling_h2o_mort_b(:)     ! coefficient in moisture-based seedling mortality
     real(r8), allocatable :: seedling_h2o_mort_c(:)     ! coefficient in moisture-based seedling mortality
     real(r8), allocatable :: seedling_root_depth(:)     ! rooting depth of seedlings [m]
     real(r8), allocatable :: seedling_light_mort_a(:)   ! light-based seedling mortality coefficient
     real(r8), allocatable :: seedling_light_mort_b(:)   ! light-based seedling mortality coefficient
     real(r8), allocatable :: background_seedling_mort(:)! background seedling mortality [yr-1]

     real(r8), allocatable :: trim_limit(:)              ! Limit to reductions in leaf area w stress (m2/m2)
     real(r8), allocatable :: trim_inc(:)                ! Incremental change in trimming function   (m2/m2)
     real(r8), allocatable :: rhol(:, :)                 ! Leaf reflectance; second dim: 1 = vis, 2 = nir
     real(r8), allocatable :: rhos(:, :)                 ! Stem reflectance; second dim: 1 = vis, 2 = nir
     real(r8), allocatable :: taul(:, :)                 ! Leaf transmittance; second dim: 1 = vis, 2 = nir
     real(r8), allocatable :: taus(:, :)                 ! Stem transmittance; second dim: 1 = vis, 2 = nir

     ! Fire Parameters (No PFT vector capabilities in their own routines)
     ! See fire/SFParamsMod.F90 for bulk of fire parameters
     ! -------------------------------------------------------------------------------------------
     real(r8), allocatable :: fire_alpha_SH(:)      ! spitfire parameter, alpha scorch height
                                                    ! Equation 16 Thonicke et al 2010

     ! Non-PARTEH Allometry Parameters
     ! --------------------------------------------------------------------------------------------


     real(r8), allocatable :: allom_frbstor_repro(:)  ! fraction of bstrore for reproduction after mortality

     ! Prescribed Physiology Mode Parameters
     real(r8), allocatable :: prescribed_npp_canopy(:)           ! this is only for the
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_npp_understory(:)       ! this is only for the
                                                                 ! prescribed physiology mode
     real(r8), allocatable :: prescribed_mortality_canopy(:)     ! this is only for the
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_mortality_understory(:) ! this is only for the
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_recruitment(:)          ! this is only for the
                                                                 ! prescribed_physiology_mode


     ! Damage Parameters

     real(r8), allocatable :: damage_frac(:)             ! Fraction of each cohort damaged per year
     real(r8), allocatable :: damage_mort_p1(:)          ! Inflection point for damage mortality function
     real(r8), allocatable :: damage_mort_p2(:)          ! Rate parameter for damage mortality function
     real(r8), allocatable :: damage_recovery_scalar(:)  ! what fraction of cohort gets to recover

     ! Nutrient Aquisition (ECA & RD)


     real(r8), allocatable :: decompmicc(:)             ! microbial decomposer biomass gC/m3
                                                        ! on root surface

     real(r8), allocatable :: vmax_nh4(:) ! maximum production rate for plant NH4 uptake   [gN/gC/s]
     real(r8), allocatable :: vmax_no3(:) ! maximum production rate for plant NO3 uptake   [gN/gC/s]
                                          ! For ECA: these rates will be applied separately to
                                          ! draw from mineralized nh4 and no3 pools independantly.
                                          ! For RD: these rates will be added, to construct a total
                                          ! N demand, which will be applied to NH4 and then NO3
                                          ! sequentially
     real(r8), allocatable :: vmax_p(:)   ! maximum production rate for plant p uptake     [gP/gC/s]

     

     
     ! ECA Parameters: See Zhu et al. Multiple soil nutrient competition between plants,
     !                     microbes, and mineral surfaces: model development, parameterization,
     !                     and example applications in several tropical forests.  Biogeosciences,
     !                     13, pp.341-363, 2016.
     ! KM: Michaeles-Menten half-saturation constants for ECA (plant–enzyme affinity)
     ! VMAX: Product of the reaction-rate and enzyme abundance for each PFT in ECA
     ! Note*: units of [gC] is grams carbon of fine-root
     
     real(r8), allocatable :: eca_km_nh4(:)   ! half-saturation constant for plant nh4 uptake  [gN/m3]
     
     real(r8), allocatable :: eca_km_no3(:)   ! half-saturation constant for plant no3 uptake  [gN/m3]
    
     real(r8), allocatable :: eca_km_p(:)     ! half-saturation constant for plant p uptake    [gP/m3]
     
     real(r8), allocatable :: eca_km_ptase(:)     ! half-saturation constant for biochemical P production [gP/m3]
     real(r8), allocatable :: eca_vmax_ptase(:)   ! maximum production rate for biochemical P prod        [gP/gC/s]
     real(r8), allocatable :: eca_alpha_ptase(:)  ! Fraction of min P generated from ptase activity
                                                  ! that is immediately sent to the plant [/]
     real(r8), allocatable :: eca_lambda_ptase(:) ! critical value for Ptase that incurs
                                                  ! biochemical production, fraction based how much
                                                  ! more in need a plant is for P versus N [/]

     ! Phenology related things

     real(r8), allocatable :: phenflush_fraction(:)       ! Maximum fraction of storage carbon used to flush leaves
                                                          ! on bud-burst [kgC/kgC]
     real(r8), allocatable :: phen_cold_size_threshold(:) ! stem/leaf drop occurs on DBH size of decidious non-woody
                                                          ! (coastal grass) plants larger than the threshold value

     ! Nutrient Aquisition parameters
     real(r8), allocatable :: prescribed_nuptake(:)   ! If there is no soil BGC model active,
                                                      ! prescribe an uptake rate for nitrogen, this is the fraction of plant demand

     real(r8), allocatable :: prescribed_puptake(:)   ! If there is no soil BGC model active,
                                                      ! prescribe an uptake rate for phosphorus
                                                      ! This is the fraction of plant demand

     ! Unassociated pft dimensioned free parameter that
     ! developers can use for testing arbitrary new hypothese
     real(r8), allocatable :: dev_arbitrary_pft(:)

     ! Parameters dimensioned by PFT and leaf age
     real(r8), allocatable :: vcmax25top(:,:)             ! maximum carboxylation rate of Rub. at 25C,
                                                          ! canopy top [umol CO2/m^2/s].  Dimensioned by
                                                          ! leaf age-class
     ! Plant Hydraulic Parameters
     ! ---------------------------------------------------------------------------------------------

     ! PFT Dimension
     real(r8), allocatable :: hydr_p_taper(:)       ! xylem taper exponent
     real(r8), allocatable :: hydr_rs2(:)           ! absorbing root radius (m)
     real(r8), allocatable :: hydr_srl(:)           ! specific root length (m g-1)
     real(r8), allocatable :: hydr_rfrac_stem(:)    ! fraction of total tree resistance from troot to canopy
     real(r8), allocatable :: hydr_avuln_gs(:)      ! shape parameter for stomatal control of water vapor exiting leaf
     real(r8), allocatable :: hydr_p50_gs(:)        ! water potential at 50% loss of stomatal conductance
     real(r8), allocatable :: hydr_k_lwp(:)         ! inner leaf humidity scaling coefficient 

     ! PFT x Organ Dimension  (organs are: 1=leaf, 2=stem, 3=transporting root, 4=absorbing root)
     ! ----------------------------------------------------------------------------------

     ! Van Genuchten PV PK curves
     real(r8), allocatable :: hydr_vg_alpha_node(:,:)   ! capilary length parameter in van Genuchten model
     real(r8), allocatable :: hydr_vg_m_node(:,:)       ! pore size distribution, m in van Genuchten 1980 model, range (0,1)
     real(r8), allocatable :: hydr_vg_n_node(:,:)       ! pore size distribution, n in van Genuchten 1980 model, range >2

     ! TFS PV-PK curves
     real(r8), allocatable :: hydr_avuln_node(:,:)  ! xylem vulernability curve shape parameter
     real(r8), allocatable :: hydr_p50_node(:,:)    ! xylem water potential at 50% conductivity loss (MPa)
     real(r8), allocatable :: hydr_epsil_node(:,:)  ! bulk elastic modulus (MPa)
     real(r8), allocatable :: hydr_pitlp_node(:,:)  ! turgor loss point (MPa)
     real(r8), allocatable :: hydr_fcap_node(:,:)   ! fraction of (1-resid_node) that is capillary in source
     real(r8), allocatable :: hydr_pinot_node(:,:)  ! osmotic potential at full turgor
     real(r8), allocatable :: hydr_kmax_node(:,:)   ! maximum xylem conductivity per unit conducting xylem area

     ! Parameters for both VG and TFS PV-PK curves
     real(r8), allocatable :: hydr_resid_node(:,:)  ! residual fraction (fraction)
     real(r8), allocatable :: hydr_thetas_node(:,:) ! saturated water content (cm3/cm3)

     ! Table that maps HLM pfts to FATES pfts for fixed biogeography mode
     ! The values are area fractions
     real(r8), allocatable :: hlm_pft_map(:,:)

     ! Land-use and land-use change related PFT parameters
     real(r8), allocatable :: harvest_pprod10(:)              ! fraction of harvest wood product that goes to 10-year product pool (remainder goes to 100-year pool)
     real(r8), allocatable :: landusechange_frac_burned(:)    ! fraction of land use change-generated and not-exported material that is burned (the remainder goes to litter)
     real(r8), allocatable :: landusechange_frac_exported(:)  ! fraction of land use change-generated wood material that is exported to wood product (the remainder is either burned or goes to litter)
     real(r8), allocatable :: landusechange_pprod10(:)        ! fraction of land use change wood product that goes to 10-year product pool (remainder goes to 100-year pool)

     ! Grazing
     real(r8), allocatable :: landuse_grazing_palatability(:) ! Relative intensity of leaf grazing/browsing per PFT (unitless 0-1)

  end type EDPftvarcon_type

  type(EDPftvarcon_type), public :: EDPftvarcon_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: TransferParamsPFT
  public :: FatesReportPFTParams
  public :: FatesCheckParams
  public :: GetDecompyFrac
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine TransferParamsPFT(pstruct)

    implicit none

    type(params_type) :: pstruct         ! Data structure containing all parameters and dimensions
    type(param_type),pointer :: param_p  ! Pointer to one specific parameter
    integer                  :: numpft
    integer                  :: num_hlm_pft
    integer                  :: num_ageclass
    integer                  :: num_hydrorgan
    
    numpft       = pstruct%GetDimSizeFromName('fates_pft')
    num_hlm_pft  = pstruct%GetDimSizeFromName('fates_hlm_pftno')
    num_ageclass = pstruct%GetDimSizeFromName('fates_leafage_class')
    num_hydrorgan = pstruct%GetDimSizeFromName('fates_hydr_organs')
    
    ! Section 1: 1D PFT dimension
    ! --------------------------------------------------------------------------
    
    param_p => pstruct%GetParamFromName('fates_mort_freezetol')
    allocate(EDPftvarcon_inst%freezetol(numpft))
    EDPftvarcon_inst%freezetol(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_height_min')
    allocate(EDPftvarcon_inst%hgt_min(numpft))
    EDPftvarcon_inst%hgt_min(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_fire_bark_scaler')
    allocate(EDPftvarcon_inst%bark_scaler(numpft))
    EDPftvarcon_inst%bark_scaler(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_fire_crown_kill')
    allocate(EDPftvarcon_inst%crown_kill(numpft))
    EDPftvarcon_inst%crown_kill(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_init_density')
    allocate(EDPftvarcon_inst%initd(numpft))
    EDPftvarcon_inst%initd(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_seed_supplement')
    allocate(EDPftvarcon_inst%seed_suppl(numpft))
    EDPftvarcon_inst%seed_suppl(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_leaf_flab')
    allocate(EDPftvarcon_inst%lf_flab(numpft))
    EDPftvarcon_inst%lf_flab(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_leaf_fcel')
    allocate(EDPftvarcon_inst%lf_fcel(numpft))
    EDPftvarcon_inst%lf_fcel(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_leaf_flig')
    allocate(EDPftvarcon_inst%lf_flig(numpft))
    EDPftvarcon_inst%lf_flig(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_fnrt_flab')
    allocate(EDPftvarcon_inst%fr_flab(numpft))
    EDPftvarcon_inst%fr_flab(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_fnrt_fcel')
    allocate(EDPftvarcon_inst%fr_fcel(numpft))
    EDPftvarcon_inst%fr_fcel(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_fnrt_flig')
    allocate(EDPftvarcon_inst%fr_flig(numpft))
    EDPftvarcon_inst%fr_flig(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_rad_leaf_xl')
    allocate(EDPftvarcon_inst%xl(numpft))
    EDPftvarcon_inst%xl(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_rad_leaf_clumping_index')
    allocate(EDPftvarcon_inst%clumping_index(numpft))
    EDPftvarcon_inst%clumping_index(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_nonhydro_smpso')
    allocate(EDPftvarcon_inst%smpso(numpft))
    EDPftvarcon_inst%smpso(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_nonhydro_smpsc')
    allocate(EDPftvarcon_inst%smpsc(numpft))
    EDPftvarcon_inst%smpsc(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_leaf_vert_scaler_coeff1')
    allocate(EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(numpft))
    EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_leaf_vert_scaler_coeff2')
    allocate(EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(numpft))
    EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_prescribed_npp_canopy')
    allocate(EDPftvarcon_inst%prescribed_npp_canopy(numpft))
    EDPftvarcon_inst%prescribed_npp_canopy(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_prescribed_npp_understory')
    allocate(EDPftvarcon_inst%prescribed_npp_understory(numpft))
    EDPftvarcon_inst%prescribed_npp_understory(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_prescribed_canopy')
    allocate(EDPftvarcon_inst%prescribed_mortality_canopy(numpft))
    EDPftvarcon_inst%prescribed_mortality_canopy(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_prescribed_understory')
    allocate(EDPftvarcon_inst%prescribed_mortality_understory(numpft))
    EDPftvarcon_inst%prescribed_mortality_understory(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_prescribed_rate')
    allocate(EDPftvarcon_inst%prescribed_recruitment(numpft))
    EDPftvarcon_inst%prescribed_recruitment(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_damage_frac')
    allocate(EDPftvarcon_inst%damage_frac(numpft))
    EDPftvarcon_inst%damage_frac(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_damage_mort_p1')
    allocate(EDPftvarcon_inst%damage_mort_p1(numpft))
    EDPftvarcon_inst%damage_mort_p1(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_damage_mort_p2')
    allocate(EDPftvarcon_inst%damage_mort_p2(numpft))
    EDPftvarcon_inst%damage_mort_p2(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_damage_recovery_scalar')
    allocate(EDPftvarcon_inst%damage_recovery_scalar(numpft))
    EDPftvarcon_inst%damage_recovery_scalar(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_fire_alpha_SH')
    allocate(EDPftvarcon_inst%fire_alpha_SH(numpft))
    EDPftvarcon_inst%fire_alpha_SH(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_allom_frbstor_repro')
    allocate(EDPftvarcon_inst%allom_frbstor_repro(numpft))
    EDPftvarcon_inst%allom_frbstor_repro(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_p_taper')
    allocate(EDPftvarcon_inst%hydr_p_taper(numpft))
    EDPftvarcon_inst%hydr_p_taper(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_rs2')
    allocate(EDPftvarcon_inst%hydr_rs2(numpft))
    EDPftvarcon_inst%hydr_rs2(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_srl')
    allocate(EDPftvarcon_inst%hydr_srl(numpft))
    EDPftvarcon_inst%hydr_srl(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_rfrac_stem')
    allocate(EDPftvarcon_inst%hydr_rfrac_stem(numpft))
    EDPftvarcon_inst%hydr_rfrac_stem(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_avuln_gs')
    allocate(EDPftvarcon_inst%hydr_avuln_gs(numpft))
    EDPftvarcon_inst%hydr_avuln_gs(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_p50_gs')
    allocate(EDPftvarcon_inst%hydr_p50_gs(numpft))
    EDPftvarcon_inst%hydr_p50_gs(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_k_lwp')
    allocate(EDPftvarcon_inst%hydr_k_lwp(numpft))
    EDPftvarcon_inst%hydr_k_lwp(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_bmort')
    allocate(EDPftvarcon_inst%bmort(numpft))
    EDPftvarcon_inst%bmort(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_scalar_coldstress')
    allocate(EDPftvarcon_inst%mort_scalar_coldstress(numpft))
    EDPftvarcon_inst%mort_scalar_coldstress(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_scalar_cstarvation')
    allocate(EDPftvarcon_inst%mort_scalar_cstarvation(numpft))
    EDPftvarcon_inst%mort_scalar_cstarvation(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_scalar_hydrfailure')
    allocate(EDPftvarcon_inst%mort_scalar_hydrfailure(numpft))
    EDPftvarcon_inst%mort_scalar_hydrfailure(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_upthresh_cstarvation')
    allocate(EDPftvarcon_inst%mort_upthresh_cstarvation(numpft))
    EDPftvarcon_inst%mort_upthresh_cstarvation(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_ip_size_senescence')
    allocate(EDPftvarcon_inst%mort_ip_size_senescence(numpft))
    EDPftvarcon_inst%mort_ip_size_senescence(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_r_size_senescence')
    allocate(EDPftvarcon_inst%mort_r_size_senescence(numpft))
    EDPftvarcon_inst%mort_r_size_senescence(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_ip_age_senescence')
    allocate(EDPftvarcon_inst%mort_ip_age_senescence(numpft))
    EDPftvarcon_inst%mort_ip_age_senescence(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_r_age_senescence')
    allocate(EDPftvarcon_inst%mort_r_age_senescence(numpft))
    EDPftvarcon_inst%mort_r_age_senescence(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_scalar_coldstress')
    allocate(EDPftvarcon_inst%mort_scalar_coldstress(numpft))
    EDPftvarcon_inst%mort_scalar_coldstress(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_scalar_cstarvation')
    allocate(EDPftvarcon_inst%mort_scalar_cstarvation(numpft))
    EDPftvarcon_inst%mort_scalar_cstarvation(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_upthresh_cstarvation')
    allocate(EDPftvarcon_inst%mort_upthresh_cstarvation(numpft))
    EDPftvarcon_inst%mort_upthresh_cstarvation(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_hf_sm_threshold')
    allocate(EDPftvarcon_inst%hf_sm_threshold(numpft))
    EDPftvarcon_inst%hf_sm_threshold(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_mort_hf_flc_threshold')
    allocate(EDPftvarcon_inst%hf_flc_threshold(numpft))
    EDPftvarcon_inst%hf_flc_threshold(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_recruit_seed_germination_rate')
    allocate(EDPftvarcon_inst%germination_rate(numpft))
    EDPftvarcon_inst%germination_rate(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_repro_frac_seed')
    allocate(EDPftvarcon_inst%repro_frac_seed(numpft))
    EDPftvarcon_inst%repro_frac_seed(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_a_emerg')
    allocate(EDPftvarcon_inst%a_emerg(numpft))
    EDPftvarcon_inst%a_emerg(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_b_emerg')
    allocate(EDPftvarcon_inst%b_emerg(numpft))
    EDPftvarcon_inst%b_emerg(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_par_crit_germ')
    allocate(EDPftvarcon_inst%par_crit_germ(numpft))
    EDPftvarcon_inst%par_crit_germ(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_psi_emerg')
    allocate(EDPftvarcon_inst%seedling_psi_emerg(numpft))
    EDPftvarcon_inst%seedling_psi_emerg(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_psi_crit')
    allocate(EDPftvarcon_inst%seedling_psi_crit(numpft))
    EDPftvarcon_inst%seedling_psi_crit(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_light_rec_a')
    allocate(EDPftvarcon_inst%seedling_light_rec_a(numpft))
    EDPftvarcon_inst%seedling_light_rec_a(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_light_rec_b')
    allocate(EDPftvarcon_inst%seedling_light_rec_b(numpft))
    EDPftvarcon_inst%seedling_light_rec_b(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_mdd_crit')
    allocate(EDPftvarcon_inst%seedling_mdd_crit(numpft))
    EDPftvarcon_inst%seedling_mdd_crit(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_h2o_mort_a')
    allocate(EDPftvarcon_inst%seedling_h2o_mort_a(numpft))
    EDPftvarcon_inst%seedling_h2o_mort_a(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_h2o_mort_b')
    allocate(EDPftvarcon_inst%seedling_h2o_mort_b(numpft))
    EDPftvarcon_inst%seedling_h2o_mort_b(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_h2o_mort_c')
    allocate(EDPftvarcon_inst%seedling_h2o_mort_c(numpft))
    EDPftvarcon_inst%seedling_h2o_mort_c(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_root_depth')
    allocate(EDPftvarcon_inst%seedling_root_depth(numpft))
    EDPftvarcon_inst%seedling_root_depth(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_light_mort_a')
    allocate(EDPftvarcon_inst%seedling_light_mort_a(numpft))
    EDPftvarcon_inst%seedling_light_mort_a(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_light_mort_b')
    allocate(EDPftvarcon_inst%seedling_light_mort_b(numpft))
    EDPftvarcon_inst%seedling_light_mort_b(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trs_seedling_background_mort')
    allocate(EDPftvarcon_inst%background_seedling_mort(numpft))
    EDPftvarcon_inst%background_seedling_mort(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_frag_seed_decay_rate')
    allocate(EDPftvarcon_inst%seed_decay_rate(numpft))
    EDPftvarcon_inst%seed_decay_rate(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_seed_dispersal_pdf_scale')
    allocate(EDPftvarcon_inst%seed_dispersal_pdf_scale(numpft))
    EDPftvarcon_inst%seed_dispersal_pdf_scale(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_seed_dispersal_pdf_shape')
    allocate(EDPftvarcon_inst%seed_dispersal_pdf_shape(numpft))
    EDPftvarcon_inst%seed_dispersal_pdf_shape(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_seed_dispersal_max_dist')
    allocate(EDPftvarcon_inst%seed_dispersal_max_dist(numpft))
    EDPftvarcon_inst%seed_dispersal_max_dist(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_seed_dispersal_fraction')
    allocate(EDPftvarcon_inst%seed_dispersal_fraction(numpft))         
    EDPftvarcon_inst%seed_dispersal_fraction(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trim_limit')
    allocate(EDPftvarcon_inst%trim_limit(numpft))
    EDPftvarcon_inst%trim_limit(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_trim_inc')
    allocate(EDPftvarcon_inst%trim_inc(numpft))
    EDPftvarcon_inst%trim_inc(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_turb_leaf_diameter')
    allocate(EDPftvarcon_inst%dleaf(numpft))
    EDPftvarcon_inst%dleaf(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_turb_z0mr')
    allocate(EDPftvarcon_inst%z0mr(numpft))
    EDPftvarcon_inst%z0mr(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_turb_displar')
    allocate(EDPftvarcon_inst%displar(numpft))
    EDPftvarcon_inst%displar(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_phen_flush_fraction')
    allocate(EDPftvarcon_inst%phenflush_fraction(numpft))
    EDPftvarcon_inst%phenflush_fraction(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_phen_cold_size_threshold')
    allocate(EDPftvarcon_inst%phen_cold_size_threshold(numpft))
    EDPftvarcon_inst%phen_cold_size_threshold(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_prescribed_nuptake')
    allocate(EDPftvarcon_inst%prescribed_nuptake(numpft))
    EDPftvarcon_inst%prescribed_nuptake(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_prescribed_puptake')
    allocate(EDPftvarcon_inst%prescribed_puptake(numpft))
    EDPftvarcon_inst%prescribed_puptake(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_dev_arbitrary_pft')
    allocate(EDPftvarcon_inst%dev_arbitrary_pft(numpft))
    EDPftvarcon_inst%dev_arbitrary_pft(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_decompmicc')
    allocate(EDPftvarcon_inst%decompmicc(numpft))
    EDPftvarcon_inst%decompmicc(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_km_nh4')
    allocate(EDPftvarcon_inst%eca_km_nh4(numpft))
    EDPftvarcon_inst%eca_km_nh4(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_vmax_nh4')
    allocate(EDPftvarcon_inst%vmax_nh4(numpft))
    EDPftvarcon_inst%vmax_nh4(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_km_no3')
    allocate(EDPftvarcon_inst%eca_km_no3(numpft))
    EDPftvarcon_inst%eca_km_no3(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_vmax_no3')
    allocate(EDPftvarcon_inst%vmax_no3(numpft))
    EDPftvarcon_inst%vmax_no3(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_km_p')
    allocate(EDPftvarcon_inst%eca_km_p(numpft))
    EDPftvarcon_inst%eca_km_p(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_vmax_p')
    allocate(EDPftvarcon_inst%vmax_p(numpft))
    EDPftvarcon_inst%vmax_p(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_km_ptase')
    allocate(EDPftvarcon_inst%eca_km_ptase(numpft))
    EDPftvarcon_inst%eca_km_ptase(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_vmax_ptase')
    allocate(EDPftvarcon_inst%eca_vmax_ptase(numpft))
    EDPftvarcon_inst%eca_vmax_ptase(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_alpha_ptase')
    allocate(EDPftvarcon_inst%eca_alpha_ptase(numpft))
    EDPftvarcon_inst%eca_alpha_ptase(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_cnp_eca_lambda_ptase')
    allocate(EDPftvarcon_inst%eca_lambda_ptase(numpft))
    EDPftvarcon_inst%eca_lambda_ptase(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_landuse_harvest_pprod10')
    allocate(EDPftvarcon_inst%harvest_pprod10(numpft))
    EDPftvarcon_inst%harvest_pprod10(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_landuse_luc_frac_burned')
    allocate(EDPftvarcon_inst%landusechange_frac_burned(numpft))
    EDPftvarcon_inst%landusechange_frac_burned(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_landuse_luc_frac_exported')
    allocate(EDPftvarcon_inst%landusechange_frac_exported(numpft))
    EDPftvarcon_inst%landusechange_frac_exported(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_landuse_luc_pprod10')
    allocate(EDPftvarcon_inst%landusechange_pprod10(numpft))
    EDPftvarcon_inst%landusechange_pprod10(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_landuse_grazing_palatability')
    allocate(EDPftvarcon_inst%landuse_grazing_palatability(numpft))
    EDPftvarcon_inst%landuse_grazing_palatability(:) = param_p%r_data_1d(:)

    ! Section 2: 2D PFT x HLM-PFT dimension
    ! --------------------------------------------------------------------------
    param_p => pstruct%GetParamFromName('fates_hlm_pft_map')
    allocate(EDPftvarcon_inst%hlm_pft_map(numpft,num_hlm_pft))
    EDPftvarcon_inst%hlm_pft_map(:,:) = param_p%r_data_2d(:,:)

    ! Section 3: 2D PFT x vis/nir dimension
    ! --------------------------------------------------------------------------
    
    allocate(EDPftvarcon_inst%rhol(numpft,2))

    param_p => pstruct%GetParamFromName('fates_rad_leaf_rhovis')
    EDPftvarcon_inst%rhol(:,ivis) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_rad_leaf_rhonir')
    EDPftvarcon_inst%rhol(:,inir) = param_p%r_data_1d(:)

    allocate(EDPftvarcon_inst%rhos(numpft,2))

    param_p => pstruct%GetParamFromName('fates_rad_stem_rhovis')
    EDPftvarcon_inst%rhos(:,ivis) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_rad_stem_rhonir')
    EDPftvarcon_inst%rhos(:,inir) = param_p%r_data_1d(:)
    
    allocate(EDPftvarcon_inst%taul(numpft,2))
    
    param_p => pstruct%GetParamFromName('fates_rad_leaf_tauvis')
    EDPftvarcon_inst%taul(:,ivis) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_rad_leaf_taunir')
    EDPftvarcon_inst%taul(:,nir) = param_p%r_data_1d(:)

    allocate(EDPftvarcon_inst%taus(numpft,2))
    
    param_p => pstruct%GetParamFromName('fates_rad_stem_tauvis')
    EDPftvarcon_inst%taus(:,ivis) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_rad_stem_taunir')
    EDPftvarcon_inst%taus(:,nir) = param_p%r_data_1d(:)

    ! Section 4: 2D PFT x leaf age dimension
    ! --------------------------------------------------------------------------
    
    param_p => pstruct%GetParamFromName('fates_leaf_vcmax25top')
    allocate(EDPftvarcon_inst%vcmax25top(numpft,num_ageclass))
    EDPftvarcon_inst%vcmax25top(:,:) = param_p%r_data_2d(:,:)

    ! Section 5: 2D PFT x hydro organ dimension
    ! --------------------------------------------------------------------------

    param_p => pstruct%GetParamFromName('fates_hydro_vg_alpha_node')
    allocate(EDPftvarcon_inst%hydr_vg_alpha_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_alpha_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_vg_m_node')
    allocate(EDPftvarcon_inst%hydr_vg_m_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_m_node(:,:) = param_p%r_data_2d(:,:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_vg_n_node')
    allocate(EDPftvarcon_inst%hydr_vg_n_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_n_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_avuln_node')
    allocate(EDPftvarcon_inst%hydr_avuln_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_avuln_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_p50_node')
    allocate(EDPftvarcon_inst%%hydr_p50_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%%hydr_p50_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_thetas_node')
    allocate(EDPftvarcon_inst%hydr_thetas_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_thetas_node(:,:) = param_p%r_data_2d(:,:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_epsil_node')
    allocate(EDPftvarcon_inst%hydr_epsil_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_epsil_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_pitlp_node')
    allocate(EDPftvarcon_inst%hydr_pitlp_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_pitlp_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_resid_node')
    allocate(EDPftvarcon_inst%hydr_resid_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_resid_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_fcap_node')
    allocate(EDPftvarcon_inst%hydr_fcap_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_fcap_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_pinot_node')
    allocate(EDPftvarcon_inst%hydr_pinot_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_pinot_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_kmax_node')
    allocate(EDPftvarcon_inst%hydr_kmax_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_kmax_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_vg_alpha_node')
    allocate(EDPftvarcon_inst%hydr_vg_alpha_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_alpha_node(:,:) = param_p%r_data_2d(:,:)
    
    param_p => pstruct%GetParamFromName('fates_hydro_vg_m_node')
    allocate(EDPftvarcon_inst%hydr_vg_m_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_m_node(:,:) = param_p%r_data_2d(:,:)

    param_p => pstruct%GetParamFromName('fates_hydro_vg_n_node')
    allocate(EDPftvarcon_inst%hydr_vg_n_node(numpft,num_hydrorgan))
    EDPftvarcon_inst%hydr_vg_n_node(:,:) = param_p%r_data_2d(:,:)


    return
  end subroutine TransferParamsPFT

  ! ===============================================================================================

  subroutine FatesReportPFTParams(is_master)

     ! Argument
     logical, intent(in) :: is_master  ! Only log if this is the master proc

     logical, parameter :: debug_report = .false.
     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft,ipft

     npft = size(EDPftvarcon_inst%initd,1)

     if(debug_report .and. is_master) then

        if(npft>100)then
           write(fates_log(),*) 'you are trying to report pft parameters during initialization'
           write(fates_log(),*) 'but you have so many that it is over-running the format spec'
           write(fates_log(),*) 'simply bump up the muptiplier in parameter fmt0 shown above'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        write(fates_log(),*) '-----------  FATES PFT Parameters -----------------'
        write(fates_log(),fmt0) 'freezetol = ',EDPftvarcon_inst%freezetol
        write(fates_log(),fmt0) 'hgt_min = ',EDPftvarcon_inst%hgt_min
        write(fates_log(),fmt0) 'dleaf = ',EDPftvarcon_inst%dleaf
        write(fates_log(),fmt0) 'z0mr = ',EDPftvarcon_inst%z0mr
        write(fates_log(),fmt0) 'displar = ',EDPftvarcon_inst%displar
        write(fates_log(),fmt0) 'bark_scaler = ',EDPftvarcon_inst%bark_scaler
        write(fates_log(),fmt0) 'crown_kill = ',EDPftvarcon_inst%crown_kill
        write(fates_log(),fmt0) 'initd = ',EDPftvarcon_inst%initd
        write(fates_log(),fmt0) 'seed_suppl = ',EDPftvarcon_inst%seed_suppl
        write(fates_log(),fmt0) 'lf_flab = ',EDPftvarcon_inst%lf_flab
        write(fates_log(),fmt0) 'lf_fcel = ',EDPftvarcon_inst%lf_fcel
        write(fates_log(),fmt0) 'lf_flig = ',EDPftvarcon_inst%lf_flig
        write(fates_log(),fmt0) 'fr_flab = ',EDPftvarcon_inst%fr_flab
        write(fates_log(),fmt0) 'fr_fcel = ',EDPftvarcon_inst%fr_fcel
        write(fates_log(),fmt0) 'fr_flig = ',EDPftvarcon_inst%fr_flig
        write(fates_log(),fmt0) 'xl = ',EDPftvarcon_inst%xl
        write(fates_log(),fmt0) 'clumping_index = ',EDPftvarcon_inst%clumping_index
        
        write(fates_log(),fmt0) 'vcmax25top = ',EDPftvarcon_inst%vcmax25top
        write(fates_log(),fmt0) 'smpso = ',EDPftvarcon_inst%smpso
        write(fates_log(),fmt0) 'smpsc = ',EDPftvarcon_inst%smpsc
        write(fates_log(),fmt0) 'bmort = ',EDPftvarcon_inst%bmort
        write(fates_log(),fmt0) 'mort_ip_size_senescence = ', EDPftvarcon_inst%mort_ip_size_senescence
        write(fates_log(),fmt0) 'mort_r_size_senescence = ', EDPftvarcon_inst%mort_r_size_senescence
        write(fates_log(),fmt0) 'mort_ip_age_senescence = ', EDPftvarcon_inst%mort_ip_age_senescence
        write(fates_log(),fmt0) 'mort_r_age_senescence = ', EDPftvarcon_inst%mort_r_age_senescence
        write(fates_log(),fmt0) 'mort_scalar_coldstress = ',EDPftvarcon_inst%mort_scalar_coldstress
        write(fates_log(),fmt0) 'mort_scalar_cstarvation = ',EDPftvarcon_inst%mort_scalar_cstarvation
        write(fates_log(),fmt0) 'mort_scalar_hydrfailure = ',EDPftvarcon_inst%mort_scalar_hydrfailure
        write(fates_log(),fmt0) 'mort_upthresh_cstarvation = ',EDPftvarcon_inst%mort_upthresh_cstarvation
        write(fates_log(),fmt0) 'hf_sm_threshold = ',EDPftvarcon_inst%hf_sm_threshold
        write(fates_log(),fmt0) 'hf_flc_threshold = ',EDPftvarcon_inst%hf_flc_threshold
 
        write(fates_log(),fmt0) 'germination_timescale = ',EDPftvarcon_inst%germination_rate
        write(fates_log(),fmt0) 'seed_decay_turnover = ',EDPftvarcon_inst%seed_decay_rate
        write(fates_log(),fmt0) 'seed_dispersal_pdf_scale = ',EDPftvarcon_inst%seed_dispersal_pdf_scale
        write(fates_log(),fmt0) 'seed_dispersal_pdf_shape = ',EDPftvarcon_inst%seed_dispersal_pdf_shape
        write(fates_log(),fmt0) 'seed_dispersal_max_dist = ',EDPftvarcon_inst%seed_dispersal_max_dist
        write(fates_log(),fmt0) 'seed_dispersal_fraction = ',EDPftvarcon_inst%seed_dispersal_fraction        
        write(fates_log(),fmt0) 'repro_frac_seed = ',EDPftvarcon_inst%repro_frac_seed
        write(fates_log(),fmt0) 'a_emerg = ',EDPftvarcon_inst%a_emerg
        write(fates_log(),fmt0) 'b_emerg = ',EDPftvarcon_inst%b_emerg
        write(fates_log(),fmt0) 'par_crit_germ = ',EDPftvarcon_inst%par_crit_germ
        write(fates_log(),fmt0) 'seedling_psi_emerg = ',EDPftvarcon_inst%seedling_psi_emerg
        write(fates_log(),fmt0) 'seedling_psi_crit = ',EDPftvarcon_inst%seedling_psi_crit
        write(fates_log(),fmt0) 'seedling_mdd_crit = ',EDPftvarcon_inst%seedling_mdd_crit
        write(fates_log(),fmt0) 'seedling_light_rec_a = ',EDPftvarcon_inst%seedling_light_rec_a
        write(fates_log(),fmt0) 'seedling_light_rec_b = ',EDPftvarcon_inst%seedling_light_rec_b
        write(fates_log(),fmt0) 'background_seedling_mort = ',EDPftvarcon_inst%background_seedling_mort
        write(fates_log(),fmt0) 'seedling_root_depth = ',EDPftvarcon_inst%seedling_root_depth        
        write(fates_log(),fmt0) 'seedling_h2o_mort_a = ',EDPftvarcon_inst%seedling_h2o_mort_a        
        write(fates_log(),fmt0) 'seedling_h2o_mort_b = ',EDPftvarcon_inst%seedling_h2o_mort_b        
        write(fates_log(),fmt0) 'seedling_h2o_mort_c = ',EDPftvarcon_inst%seedling_h2o_mort_c        
        write(fates_log(),fmt0) 'trim_limit = ',EDPftvarcon_inst%trim_limit
        write(fates_log(),fmt0) 'trim_inc = ',EDPftvarcon_inst%trim_inc
        write(fates_log(),fmt0) 'rhol = ',EDPftvarcon_inst%rhol
        write(fates_log(),fmt0) 'rhos = ',EDPftvarcon_inst%rhos
        write(fates_log(),fmt0) 'taul = ',EDPftvarcon_inst%taul
        write(fates_log(),fmt0) 'taus = ',EDPftvarcon_inst%taus
        write(fates_log(),fmt0) 'phen_flush_fraction',EDpftvarcon_inst%phenflush_fraction
        write(fates_log(),fmt0) 'phen_cold_size_threshold = ',EDPftvarcon_inst%phen_cold_size_threshold
        write(fates_log(),fmt0) 'fire_alpha_SH = ',EDPftvarcon_inst%fire_alpha_SH
        write(fates_log(),fmt0) 'allom_frbstor_repro = ',EDPftvarcon_inst%allom_frbstor_repro
        write(fates_log(),fmt0) 'hydro_p_taper = ',EDPftvarcon_inst%hydr_p_taper
        write(fates_log(),fmt0) 'hydro_rs2 = ',EDPftvarcon_inst%hydr_rs2
        write(fates_log(),fmt0) 'hydro_srl = ',EDPftvarcon_inst%hydr_srl
        write(fates_log(),fmt0) 'hydro_rfrac_stem = ',EDPftvarcon_inst%hydr_rfrac_stem
        write(fates_log(),fmt0) 'hydro_avuln_gs = ',EDPftvarcon_inst%hydr_avuln_gs
        write(fates_log(),fmt0) 'hydro_k_lwp = ',EDPftvarcon_inst%hydr_k_lwp
        write(fates_log(),fmt0) 'hydro_p50_gs = ',EDPftvarcon_inst%hydr_p50_gs
        write(fates_log(),fmt0) 'hydro_avuln_node = ',EDPftvarcon_inst%hydr_avuln_node
        write(fates_log(),fmt0) 'hydro_p50_node = ',EDPftvarcon_inst%hydr_p50_node
        write(fates_log(),fmt0) 'hydro_thetas_node = ',EDPftvarcon_inst%hydr_thetas_node
        write(fates_log(),fmt0) 'hydro_epsil_node = ',EDPftvarcon_inst%hydr_epsil_node
        write(fates_log(),fmt0) 'hydro_pitlp_node = ',EDPftvarcon_inst%hydr_pitlp_node
        write(fates_log(),fmt0) 'hydro_resid_node = ',EDPftvarcon_inst%hydr_resid_node
        write(fates_log(),fmt0) 'hydro_fcap_node = ',EDPftvarcon_inst%hydr_fcap_node
        write(fates_log(),fmt0) 'hydro_pinot_node = ',EDPftvarcon_inst%hydr_pinot_node
        write(fates_log(),fmt0) 'hydro_kmax_node = ',EDPftvarcon_inst%hydr_kmax_node
        write(fates_log(),fmt0) 'hlm_pft_map = ', EDPftvarcon_inst%hlm_pft_map
        write(fates_log(),fmt0) 'hydro_vg_alpha_node  = ',EDPftvarcon_inst%hydr_vg_alpha_node
        write(fates_log(),fmt0) 'hydro_vg_m_node  = ',EDPftvarcon_inst%hydr_vg_m_node
        write(fates_log(),fmt0) 'hydro_vg_n_node  = ',EDPftvarcon_inst%hydr_vg_n_node
        write(fates_log(),fmt0) 'maintresp_leaf_vert_scaler_coeff1 = ',EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1
        write(fates_log(),fmt0) 'maintresp_leaf_vert_scaler_coeff2 = ',EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2
        write(fates_log(),*) '-------------------------------------------------'

     end if

  end subroutine FatesReportPFTParams


  ! =====================================================================================

  subroutine FatesCheckParams(is_master)

     ! ----------------------------------------------------------------------------------
     !
     ! This subroutine performs logical checks on user supplied parameters.  It cross
     ! compares various parameters and will fail if they don't make sense.
     ! Examples:
     ! A woody plant cannot have a structural biomass allometry intercept of 0, and a 
     ! non-woody plant (grass) can't have a non-zero intercept...
     ! -----------------------------------------------------------------------------------
    use FatesConstantsMod  , only : fates_check_param_set
    use FatesConstantsMod  , only : itrue, ifalse
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    
    use EDParamsMod        , only : logging_mechanical_frac, logging_collateral_frac
    use EDParamsMod        , only : logging_direct_frac,logging_export_frac
    use FatesInterfaceTypesMod, only : hlm_use_fixed_biogeog,hlm_use_sp, hlm_name
    use FatesInterfaceTypesMod, only : hlm_use_inventory_init
    use FatesInterfaceTypesMod, only : hlm_use_nocomp
    use EDParamsMod        , only : max_nocomp_pfts_by_landuse, maxpatches_by_landuse
    use FatesConstantsMod  , only : n_landuse_cats

     ! Argument
     logical, intent(in) :: is_master    ! Only log if this is the master proc

     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft     ! number of PFTs
     integer :: ipft     ! pft index
     integer :: nleafage ! size of the leaf age class array
     integer :: iage     ! leaf age class index
     integer :: norgans  ! size of the plant organ dimension
     integer  :: hlm_pft    ! used in fixed biogeog mode
     integer  :: fates_pft  ! used in fixed biogeog mode
     integer  :: i_lu    ! land use index

     real(r8) :: sumarea    ! area of PFTs in nocomp mode.
     real(r8) :: neg_lmr_temp ! temperature at which lmr would got negative 
     real(r8) :: r_0 ! base respiartion rate, PFT-dependent
     real(r8) :: lnc_top ! leaf nitrogen content at top of canopy
     
     
     npft = size(EDPftvarcon_inst%freezetol,1)

     if(.not.is_master) return

     select case (hlm_parteh_mode)
     case (prt_cnp_flex_allom_hyp)

        ! Check to see if either RD/ECA/MIC is turned on

        if (.not.( (trim(hlm_nu_com).eq.'RD') .or. (trim(hlm_nu_com).eq.'ECA'))) then
           write(fates_log(),*) 'FATES PARTEH with allometric flexible CNP must have'
           write(fates_log(),*) 'a valid BGC model enabled: RD,ECA,MIC or SYN'
           write(fates_log(),*) 'nu_comp: ',trim(hlm_nu_com)
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        
        ! If nitrogen is turned on, check to make sure there are valid ammonium
        ! parameters
        if(hlm_nitrogen_spec>0)then
           if (trim(hlm_nu_com).eq.'ECA') then

              if(any(EDpftvarcon_inst%eca_km_nh4(:)<0._r8) ) then
                 write(fates_log(),*) 'ECA with nitrogen is turned on'
                 write(fates_log(),*) 'bad ECA km value(s) for nh4: ',EDpftvarcon_inst%eca_km_nh4(:)
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if

              if(hlm_nitrogen_spec==2)then
                 if(any(EDpftvarcon_inst%eca_km_no3(:)<0._r8)) then
                    write(fates_log(),*) 'ECA with nit/denitr is turned on'
                    write(fates_log(),*) 'bad ECA km value(s) for no3: ',EDpftvarcon_inst%eca_km_no3(:)
                    write(fates_log(),*) 'Aborting'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                 end if
              end if
           end if
        end if
        
        ! If any PFTs are specified as either prescribed N or P uptake
        ! then they all must be !
        
        if (any(EDPftvarcon_inst%prescribed_nuptake(:) < -nearzero ) .or. &
             any(EDPftvarcon_inst%prescribed_nuptake(:) > 10._r8 ) ) then
           write(fates_log(),*) 'Negative values for EDPftvarcon_inst%prescribed_nuptake(:)'
           write(fates_log(),*) 'are not allowed. Reasonable ranges for this parameter are zero'
           write(fates_log(),*) 'to something slightly larger than 1, so we set a cap at 10.'
           write(fates_log(),*) 'Set to zero to turn off and use coupled nutrients.'
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        elseif (any(abs(EDPftvarcon_inst%prescribed_nuptake(:)) > nearzero )) then
           if(.not.all(abs(EDPftvarcon_inst%prescribed_nuptake(:)) > nearzero )) then
              write(fates_log(),*) 'If any PFTs are specified as having prescribed N'
              write(fates_log(),*) 'uptake, then they must all. Note, prescribed'
              write(fates_log(),*) 'rates are associated with any value abs(x)>nearzero'
              write(fates_log(),*) 'EDPftvarcon_inst%prescribed_nuptake(:):', &
                   EDPftvarcon_inst%prescribed_nuptake(:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! We are using a simple phosphatase model right now. There is
        ! no critical value (lambda) , and there is no preferential uptake (alpha).
        ! Make sure these parameters are both set to 0.

        if ((hlm_phosphorus_spec>0) .and. (trim(hlm_nu_com).eq.'ECA')) then
           if (any(abs(EDPftvarcon_inst%eca_lambda_ptase(:)) > nearzero ) ) then
              write(fates_log(),*) 'Critical Values for phosphatase in ECA are not'
              write(fates_log(),*) 'enabled right now. Please set fates_eca_lambda_ptase = 0'
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           if (any(abs(EDPftvarcon_inst%eca_alpha_ptase(:)) > nearzero ) ) then
              write(fates_log(),*) 'There is no preferential plant uptake of P from phosphatase'
              write(fates_log(),*) 'enabled right now. Please set fates_eca_alpha_ptase = 0'
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
     case (prt_carbon_allom_hyp)
        ! No additional checks needed for now.
        continue

     case default

        write(fates_log(),*) 'FATES Plant Allocation and Reactive Transport has'
        write(fates_log(),*) 'only 2 modules supported, allometric carbon and CNP.'
        write(fates_log(),*) 'fates_parteh_mode must be set to 1 or 2 in the namelist'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end select

     ! logging parameters, make sure they make sense
     if ( (logging_mechanical_frac + logging_collateral_frac + logging_direct_frac) .gt. 1._r8) then
        write(fates_log(),*) 'the sum of logging_mechanical_frac + logging_collateral_frac + logging_direct_frac'
        write(fates_log(),*) 'must be less than or equal to 1.'
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     endif

     ! Same for phosphorus
     if (any(EDPftvarcon_inst%prescribed_puptake(:) < -nearzero ) .or. &
          any(EDPftvarcon_inst%prescribed_puptake(:) > 10._r8 )) then
        write(fates_log(),*) 'Negative values for EDPftvarcon_inst%prescribed_puptake(:)'
        write(fates_log(),*) 'are not allowed. Reasonable ranges for this parameter are zero'
        write(fates_log(),*) 'to something slightly larger than 1, so we set a cap at 10.'
        write(fates_log(),*) 'Set to zero or unset to turn off and use coupled nutrients.'
        write(fates_log(),*) ' Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     elseif (any(abs(EDPftvarcon_inst%prescribed_puptake(:)) > nearzero )) then
        if(.not.all(abs(EDPftvarcon_inst%prescribed_puptake(:)) > nearzero )) then
           write(fates_log(),*) 'If any PFTs are specified as having prescribed P'
           write(fates_log(),*) 'uptake, then they must all. Note, prescribed'
           write(fates_log(),*) 'rates are associated with any value abs(x)>nearzero'
           write(fates_log(),*) 'EDPftvarcon_inst%prescribed_puptake(:):', &
                EDPftvarcon_inst%prescribed_puptake(:)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
     end if


     

     do ipft = 1,npft

        ! xl must be between -0.4 and 0.6 according to Bonan (2019) doi:10.1017/9781107339217 pg. 238
        !-----------------------------------------------------------------------------------
        if (EDPftvarcon_inst%xl(ipft) < -0.4 .or. EDPftvarcon_inst%xl(ipft) > 0.6) then
          write(fates_log(),*) 'fates_rad_leaf_xl for pft ', ipft, ' is outside the allowed range of -0.4 to 0.6'
          write(fates_log(),*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
        end if 

        ! Check that the seed dispersal parameters are all set if one of them is set
        !-----------------------------------------------------------------------------------
        if (( EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft) < fates_check_param_set ) .and. &
            (( EDPftvarcon_inst%seed_dispersal_max_dist(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_fraction(ipft) > fates_check_param_set ))) then

        write(fates_log(),*) 'Seed dispersal is on per fates_seed_dispersal_pdf_scale being set'
        write(fates_log(),*) 'Please provide values for all other seed_dispersal parameters'
        write(fates_log(),*) 'See Bullock et al. (2017) for reasonable values'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        if (( EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft) < fates_check_param_set ) .and. &
            (( EDPftvarcon_inst%seed_dispersal_max_dist(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_fraction(ipft) > fates_check_param_set ))) then

        write(fates_log(),*) 'Seed dispersal is on per fates_seed_dispersal_pdf_shape being set'
        write(fates_log(),*) 'Please provide values for all other seed_dispersal parameters'
        write(fates_log(),*) 'See Bullock et al. (2017) for reasonable values'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        
        if (( EDPftvarcon_inst%seed_dispersal_max_dist(ipft) < fates_check_param_set ) .and. &
            ((  EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft) > fates_check_param_set  ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_fraction(ipft) > fates_check_param_set ))) then

        write(fates_log(),*) 'Seed dispersal is on per seed_dispersal_max_dist being set'
        write(fates_log(),*) 'Please provide values for all other seed_dispersal parameters'
        write(fates_log(),*) 'See Bullock et al. (2017) for reasonable values'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        
        if (( EDPftvarcon_inst%seed_dispersal_fraction(ipft) < fates_check_param_set ) .and. &
        ((  EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft) > fates_check_param_set  ) .or. &
         (  EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft) > fates_check_param_set  ) .or. &
         (  EDPftvarcon_inst%seed_dispersal_max_dist(ipft) > fates_check_param_set ))) then

         write(fates_log(),*) 'Seed dispersal is on per seed_dispersal_fraction being set'
         write(fates_log(),*) 'Please provide values for all other seed_dispersal parameters'
         write(fates_log(),*) 'See Bullock et al. (2017) for reasonable values'
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
        
        ! Check that parameter ranges for the seed dispersal fraction make sense
        !-----------------------------------------------------------------------------------
        if (( EDPftvarcon_inst%seed_dispersal_fraction(ipft) < fates_check_param_set ) .and. &
            ((  EDPftvarcon_inst%seed_dispersal_fraction(ipft) > 1.0_r8 ) .or. &
             (  EDPftvarcon_inst%seed_dispersal_fraction(ipft) < 0.0_r8 ))) then

        write(fates_log(),*) 'Seed dispersal is on per seed_dispersal_fraction being set'
        write(fates_log(),*) 'Please make sure the fraction value is between 0 and 1'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        
        ! Check that parameter ranges for age-dependent mortality make sense
        !-----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%mort_ip_age_senescence(ipft) < fates_check_param_set ) .and. &
             (  EDPftvarcon_inst%mort_r_age_senescence(ipft) > fates_check_param_set ) ) then

           write(fates_log(),*) 'Age-dependent mortality is on'
           write(fates_log(),*) 'Please also set mort_r_age_senescence'
           write(fates_log(),*) 'Sensible values are between 0.03-0.06'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        ! Check that parameter ranges for size-dependent mortality make sense
        !-----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%mort_ip_size_senescence(ipft) < fates_check_param_set ) .and. &
             (  EDPftvarcon_inst%mort_r_size_senescence(ipft) > fates_check_param_set ) ) then

           write(fates_log(),*) 'Size-dependent mortality is on'
           write(fates_log(),*) 'Please also set mort_r_size_senescence'
           write(fates_log(),*) 'Sensible values are between 0.03-0.06'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        ! Check that parameter ranges for size-dependent mortality make sense
        !-----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%mort_ip_size_senescence(ipft) < 0.0_r8 ) .or. &
           ( EDPftvarcon_inst%mort_r_size_senescence(ipft) < 0.0_r8 ) ) then

           write(fates_log(),*) 'Either mort_ip_size_senescence or mort_r_size_senescence'
           write(fates_log(),*) 'is negative which makes no biological sense.'
           write(fates_log(),*) 'Sensible values for ip are between 1 and 300?'
           write(fates_log(),*) 'Sensible values for r are between 0.03-0.06'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if


        ! Check that parameter ranges for size-dependent mortality make sense
        !-----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%mort_ip_size_senescence(ipft) < 0.0_r8 ) .or. &
           ( EDPftvarcon_inst%mort_r_size_senescence(ipft) < 0.0_r8 ) ) then

           write(fates_log(),*) 'Either mort_ip_size_senescence or mort_r_size_senescence'
           write(fates_log(),*) 'is negative which makes no biological sense.'
           write(fates_log(),*) 'Sensible values for ip are between 1 and 300?'
           write(fates_log(),*) 'Sensible values for r are between 0.03-0.06'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        ! Check if the fraction of storage used for flushing deciduous trees
        ! is greater than zero, and less than or equal to 1.
        if (prt_params%phen_leaf_habit(ipft) /= ievergreen) then
           if ( ( EDPftvarcon_inst%phenflush_fraction(ipft) < nearzero ) .or. &
                ( EDPFtvarcon_inst%phenflush_fraction(ipft) > 1 ) ) then

              write(fates_log(),*) ' Deciduous plants must flush some storage carbon'
              write(fates_log(),*) ' on bud-burst. If phenflush_fraction is not greater than 0'
              write(fates_log(),*) ' it will not be able to put out any leaves. Plants need leaves.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' phen_leaf_habit: (evergreen should be ',ievergreen,'):', &
                                        int(prt_params%phen_leaf_habit(ipft))
              write(fates_log(),*) ' phenflush_fraction: ', EDPFtvarcon_inst%phenflush_fraction(ipft)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Check if freezing tolerance is within reasonable bounds
        ! ----------------------------------------------------------------------------------

        if ( ( EDPftvarcon_inst%freezetol(ipft) > 60.0_r8 ) .or. &
             ( EDPFtvarcon_inst%freezetol(ipft) < -273.1_r8 ) ) then

           write(fates_log(),*) 'Freezing tolerance was set to a strange value'
           write(fates_log(),*) ' Units should be degrees celcius. It cannot'
           write(fates_log(),*) ' be less than absolute zero, and we check to see'
           write(fates_log(),*) ' if it is greater than 60C, which would be ludicrous as well'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' freezetol: ', EDPFtvarcon_inst%freezetol(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if

        ! Check that in initial density is not equal to zero in a cold-start run
        !-----------------------------------------------------------------------------------
        
        if ( hlm_use_inventory_init == ifalse .and. & 
             abs( EDPftvarcon_inst%initd(ipft) ) < nearzero ) then
          
           write(fates_log(),*) ' In a cold start run initial density cannot be zero.'
           write(fates_log(),*) ' For a bare ground run set to initial recruit density.'
           write(fates_log(),*) ' If no-comp is on it is possible to initialize with larger  '
           write(fates_log(),*) ' plants by setting fates_recruit_init_density to a negative number'
           write(fates_log(),*) ' which will be interpreted as (absolute) initial dbh. '
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if
           

        ! Check to make sure that if a grass sapwood allometry is used, it is not
        ! a woody plant.
        if ( ( prt_params%allom_smode(ipft)==2 ) .and. (prt_params%woody(ipft)==itrue) ) then
           write(fates_log(),*) 'Allometry mode 2 is a mode that is only appropriate'
           write(fates_log(),*) 'for a grass functional type. Sapwood allometry is set with'
           write(fates_log(),*) 'fates_allom_smode in the parameter file. Woody versus non woody'
           write(fates_log(),*) 'plants are set via fates_woody in the parameter file.'
           write(fates_log(),*) 'Current settings for pft number: ',ipft
           write(fates_log(),*) 'fates_woody: true'
           write(fates_log(),*) 'fates_allom_smode: ',prt_params%allom_smode(ipft)
           write(fates_log(),*) 'Please correct this discrepancy before re-running. Aborting.'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        ! Check if fraction of storage to reproduction is between 0-1
        ! ----------------------------------------------------------------------------------

        if ( ( EDPftvarcon_inst%allom_frbstor_repro(ipft) < 0.0_r8 ) .or. &
             ( EDPftvarcon_inst%allom_frbstor_repro(ipft) > 1.0_r8 ) ) then

           write(fates_log(),*) 'fraction of storage to reproduction'
           write(fates_log(),*) ' after plants die, must be between'
           write(fates_log(),*) ' 0 and 1'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' allom_frbstor_repro: ',EDPftvarcon_inst%allom_frbstor_repro(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if


        ! Check if photosynthetic pathway is neither C3/C4
        ! ----------------------------------------------------------------------------------
        if(.not.any(lb_params%c3psn(ipft) == [c3_path_index,c4_path_index])) then

           write(fates_log(),*) ' Two photosynthetic pathways are currently supported'
           write(fates_log(),*) ' C4 plants have c3psn = ',c3_path_index
           write(fates_log(),*) ' C3 plants have c3psn = ',c4_path_index
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' c3psn(pft): ',lb_params%c3psn(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if

        if( hlm_use_fixed_biogeog .eq. itrue ) then
           ! check that the host-fates PFT map adds to one along HLM dimension so that all the HLM area
           ! goes to a FATES PFT.  Each FATES PFT can get < or > 1 of an HLM PFT.
           do hlm_pft = 1,size( EDPftvarcon_inst%hlm_pft_map,2)
              sumarea = sum(EDPftvarcon_inst%hlm_pft_map(1:npft,hlm_pft))
              if(abs(sumarea-1.0_r8).gt.nearzero)then
                 write(fates_log(),*) 'The distribution of this host land model PFT :',hlm_pft
                 write(fates_log(),*) 'into FATES PFTs, does not add up to 1.0.'
                 write(fates_log(),*) 'Error is:',sumarea-1.0_r8
                 write(fates_log(),*) 'and the hlm_pft_map is:', EDPftvarcon_inst%hlm_pft_map(1:npft,hlm_pft)
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end do !hlm_pft
        end if
        
     end do !ipft


     ! if nocomp is enabled, check to make sure the max number of nocomp PFTs per land use is
     ! less than or equal to the max number of patches per land use. (unless this is an
     ! SP run, then all PFTS are tracked on the primary LU and the others are allocated
     ! zero patch space

     if ( hlm_use_nocomp .eq. itrue .and. hlm_use_sp.eq.ifalse) then
        do i_lu = 1, n_landuse_cats
           if (max_nocomp_pfts_by_landuse(i_lu) .gt. maxpatches_by_landuse(i_lu)) then
              write(fates_log(),*) 'The max number of nocomp PFTs must all be less than or equal to the number of patches, for a given land use type'
              write(fates_log(),*) 'land use index:',i_lu
              write(fates_log(),*) 'max_nocomp_pfts_by_landuse(i_lu):', max_nocomp_pfts_by_landuse(i_lu)
              write(fates_log(),*) 'maxpatches_by_landuse(i_lu):', maxpatches_by_landuse(i_lu)
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end do
     endif
     

     ! Check the temperature at which Rdark would become negative for each PFT -
     ! given their parameters
     !------------------------------------------------------------------------------------
     do ipft = 1,npft
        
        r_0 = lb_params%maintresp_leaf_atkin2017_baserate(ipft)

        lnc_top = prt_params%nitr_stoich_p1(ipft, prt_params%organ_param_id(leaf_organ))
        
        ! From LeafLayerMaintenanceRespiration_Atkin_etal_2017
        ! r_t_ref = nscaler * (r_0 + r_1 * lnc_top + r_2 * max(0._r8, (tgrowth - tfrz) ))

        ! find temperature at which whole term is negative
        neg_lmr_temp = ( -1._r8 * (  r_0  + lmr_r_1 * lnc_top ) ) / lmr_r_2

        write(fates_log(),*)  'PFT  ',  ipft
        write(fates_log(),*)  'will have  negative Rdark at ', neg_lmr_temp, 'degrees C' 
        write(fates_log(),*)  'with these values of slatop, nitrogen stoichiometry and' 
        write(fates_log(),*)  'maintresp_leaf_atkin2017_baserate.'
        write(fates_log(),*)  'See LeafLayerMaintenanceRespiration_Atkin_etal_2017 in '
        write(fates_log(),*)  'FatesPlantRespPhotosynthMod'  
     
     end do ! ipft

     ! Check to make sure the btran limitation models are within expected ranges

     do ipft = 1,npft
        
        if( lb_params%stomatal_btran_model(ipft) < 0 .or. &
            lb_params%stomatal_btran_model(ipft) > btran_on_gs_gs02 ) then

           write(fates_log(),*)  'PFT  ',  ipft
           write(fates_log(),*)  'Undefined fates_leaf_stomatal_btran_model = ',lb_params%stomatal_btran_model
           write(fates_log(),*)  'See biogeophys/LeafbiophysicsMod.F90 btran_on_gs_* for model types'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        if( lb_params%agross_btran_model(ipft) < 0 .or. &
            lb_params%agross_btran_model(ipft) > btran_on_ag_vcmax_jmax ) then
           write(fates_log(),*)  'PFT  ',  ipft
           write(fates_log(),*)  'Undefined fates_leaf_agross_btran_model = ',lb_params%agross_btran_model
           write(fates_log(),*)  'See biogeophys/LeafbiophysicsMod.F90 btran_on_ag_* for model types'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

     end do
     
     
     

!!    ! Checks for HYDRO
!!    if( hlm_use_planthydro == itrue ) then
!!
!!        do ipft=1,numpft
!!
!!            ! Calculate fine-root density and see if the result
!!            ! is reasonable.
!!            ! kg/m3
!!
!!            dens_aroot = 1._r8/(g_per_kg*pi_const*EDPftvarcon_inst%hydr_rs2(ipft)**2.0_r8*EDPftvarcon_inst%hydr_srl(ipft))
!!
!!            if(rho_aroot>max_dens_aroot .or. dens_aroot<min_dens_aroot)then
!!                write(fates_log(),*) 'The two parameters controling fine root radius'
!!                write(fates_log(),*) 'and specific root length, have generated'
!!                write(fates_log(),*) 'a strange root density.'
!!                call endrun(msg=errMsg(sourcefile, __LINE__))
!!            end if
!!
!!        end if
!!    end do

     return
  end subroutine FatesCheckParams


  ! =====================================================================================

  function GetDecompyFrac(pft,organ_id,dcmpy) result(decompy_frac)


      ! This simple routine matches the correct decomposibility pool's
      ! material fraction with the pft parameter data.

      integer, intent(in) :: pft
      integer, intent(in) :: organ_id
      integer, intent(in) :: dcmpy
      real(r8) :: decompy_frac

      ! Decomposability for leaves
      if(organ_id == leaf_organ)then
         select case(dcmpy)
         case(ilabile)
            decompy_frac=EDPftvarcon_inst%lf_flab(pft)
         case(icellulose)
            decompy_frac=EDPftvarcon_inst%lf_fcel(pft)
         case(ilignin)
            decompy_frac=EDPftvarcon_inst%lf_flig(pft)
         case default
            write(fates_log(),*) 'Unknown decompositiblity pool index: ',dcmpy
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select
      ! Decomposability for fine-roots
      elseif(organ_id == fnrt_organ) then
         select case(dcmpy)
         case(ilabile)
            decompy_frac=EDPftvarcon_inst%fr_flab(pft)
         case(icellulose)
            decompy_frac=EDPftvarcon_inst%fr_fcel(pft)
         case(ilignin)
            decompy_frac=EDPftvarcon_inst%fr_flig(pft)
         case default
            write(fates_log(),*) 'Unknown decompositiblity pool index: ',dcmpy
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select
      else
         write(fates_log(),*) 'Unknown parteh organ index: ',organ_id
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      return
  end function GetDecompyFrac


end module EDPftvarcon

