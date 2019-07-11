module EDPftvarcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use EDTypesMod  ,   only : maxSWb, ivis, inir
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : itrue, ifalse
  use FatesConstantsMod, only : years_per_day
  use FatesGlobals,   only : fates_log
  use FatesGlobals,   only : endrun => fates_endrun

  use PRTGenericMod,  only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,  only : prt_carbon_allom_hyp
  use PRTGenericMod,  only : num_organ_types
  use PRTGenericMod,  only : leaf_organ, fnrt_organ, store_organ
  use PRTGenericMod,  only : sapw_organ, struct_organ, repro_organ

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
     real(r8), allocatable :: wood_density(:)        ! wood density  g cm^-3  ...
     real(r8), allocatable :: hgt_min(:)             ! sapling height m
     real(r8), allocatable :: dbh_repro_threshold(:) ! diameter at which mature plants shift allocation
     real(r8), allocatable :: dleaf(:)               ! leaf characteristic dimension length (m)
     real(r8), allocatable :: z0mr(:)                ! ratio of roughness length of vegetation to height (-) 
     real(r8), allocatable :: displar(:)             ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: cushion(:)             ! labile carbon storage target as multiple of leaf pool.
     real(r8), allocatable :: leaf_stor_priority(:)  ! leaf turnover vs labile carbon use prioritisation
                                                     ! (1 = lose  leaves, 0 = use store).
     real(r8), allocatable :: crown(:)               ! fraction of the height of the plant 
                                                     ! that is occupied by crown. For fire model. 
     real(r8), allocatable :: bark_scaler(:)         ! scaler from dbh to bark thickness. For fire model.
     real(r8), allocatable :: crown_kill(:)          ! scaler on fire death. For fire model. 
     real(r8), allocatable :: initd(:)               ! initial seedling density 
     real(r8), allocatable :: seed_suppl(:)           ! seeds that come from outside the gridbox.
     real(r8), allocatable :: BB_slope(:)            ! ball berry slope parameter
     
     real(r8), allocatable :: seed_alloc_mature(:)   ! fraction of carbon balance allocated to 
                                                     ! clonal reproduction.
     real(r8), allocatable :: seed_alloc(:)          ! fraction of carbon balance allocated to seeds.
     real(r8), allocatable :: c2b(:)                 ! Carbon to biomass multiplier [kg/kgC]
     real(r8), allocatable :: woody(:)               ! Does the plant have wood?      (1=yes, 0=no)

                                                     ! The following three PFT classes 
                                                     ! are mutually exclusive
     real(r8), allocatable :: stress_decid(:)        ! Is the plant stress deciduous? (1=yes, 0=no)
     real(r8), allocatable :: season_decid(:)        ! Is the plant seasonally deciduous (1=yes, 0=no)
     real(r8), allocatable :: evergreen(:)           ! Is the plant an evergreen (1=yes, 0=no)

     real(r8), allocatable :: slamax(:)              ! Maximum specific leaf area of plant (at bottom) [m2/gC]
     real(r8), allocatable :: slatop(:)              ! Specific leaf area at canopy top [m2/gC]
     
     real(r8), allocatable :: roota_par(:)           ! Normalized Root profile scaling parameter A
     real(r8), allocatable :: rootb_par(:)           ! Normalized root profile scaling parameter B
     real(r8), allocatable :: lf_flab(:)             ! Leaf litter labile fraction [-]
     real(r8), allocatable :: lf_fcel(:)             ! Leaf litter cellulose fraction [-]
     real(r8), allocatable :: lf_flig(:)             ! Leaf litter lignan fraction [-]
     real(r8), allocatable :: fr_flab(:)             ! Fine-root litter labile fraction [-]
     real(r8), allocatable :: fr_fcel(:)             ! Fine-root litter cellulose fraction [-]
     real(r8), allocatable :: fr_flig(:)             ! Fine-root litter lignatn fraction [-]
     real(r8), allocatable :: xl(:)                  ! Leaf-stem orientation index
     real(r8), allocatable :: clumping_index(:)      ! factor describing how much self-occlusion 
                                                     ! of leaf scattering elements 
                                                     ! decreases light interception
     real(r8), allocatable :: c3psn(:)               ! index defining the photosynthetic 
                                                     ! pathway C4 = 0,  C3 = 1
    
     real(r8), allocatable :: smpso(:)               ! Soil water potential at full stomatal opening 
                                                     ! (non-HYDRO mode only) [mm]
     real(r8), allocatable :: smpsc(:)               ! Soil water potential at full stomatal closure 
                                                     ! (non-HYDRO mode only) [mm]


     real(r8), allocatable :: maintresp_reduction_curvature(:) ! curvature of MR reduction as f(carbon storage), 
                                                               ! 1=linear, 0=very curved
     real(r8), allocatable :: maintresp_reduction_intercept(:) ! intercept of MR reduction as f(carbon storage), 
                                                               ! 0=no throttling, 1=max throttling
     real(r8), allocatable :: bmort(:)
     real(r8), allocatable :: mort_scalar_coldstress(:)
     real(r8), allocatable :: mort_scalar_cstarvation(:)
     real(r8), allocatable :: mort_scalar_hydrfailure(:)
     real(r8), allocatable :: hf_sm_threshold(:)
     real(r8), allocatable :: hf_flc_threshold(:)
     real(r8), allocatable :: vcmaxha(:)
     real(r8), allocatable :: jmaxha(:)
     real(r8), allocatable :: tpuha(:)
     real(r8), allocatable :: vcmaxhd(:)
     real(r8), allocatable :: jmaxhd(:)
     real(r8), allocatable :: tpuhd(:)
     real(r8), allocatable :: vcmaxse(:)
     real(r8), allocatable :: jmaxse(:)
     real(r8), allocatable :: tpuse(:)
     real(r8), allocatable :: germination_timescale(:)
     real(r8), allocatable :: seed_decay_turnover(:)
     
     real(r8), allocatable :: trim_limit(:)              ! Limit to reductions in leaf area w stress (m2/m2)
     real(r8), allocatable :: trim_inc(:)                ! Incremental change in trimming function   (m2/m2)
     real(r8), allocatable :: rhol(:, :)
     real(r8), allocatable :: rhos(:, :)
     real(r8), allocatable :: taul(:, :)
     real(r8), allocatable :: taus(:, :)
     real(r8), allocatable :: rootprof_beta(:, :)

     ! Fire Parameters (No PFT vector capabilities in their own routines)
     ! See fire/SFParamsMod.F90 for bulk of fire parameters
     ! -------------------------------------------------------------------------------------------
     real(r8), allocatable :: fire_alpha_SH(:)      ! spitfire parameter, alpha scorch height
                                                    ! Equation 16 Thonicke et al 2010

     ! Allometry Parameters
     ! -------------------------------------------------------------------------------------------- 
     real(r8), allocatable :: allom_dbh_maxheight(:) ! dbh at which height growth ceases
     
     real(r8), allocatable :: allom_hmode(:)        ! height allometry function type
     real(r8), allocatable :: allom_lmode(:)        ! maximum leaf allometry function type
     real(r8), allocatable :: allom_fmode(:)        ! maximum root allometry function type
     real(r8), allocatable :: allom_amode(:)        ! AGB allometry function type
     real(r8), allocatable :: allom_cmode(:)        ! Coarse root allometry function type
     real(r8), allocatable :: allom_smode(:)        ! sapwood allometry function type
     real(r8), allocatable :: allom_stmode(:)       ! storage allometry functional type 
                                                    ! (HARD-CODED FOR TIME BEING, RGK 11-2017)
     real(r8), allocatable :: allom_la_per_sa_int(:)   ! Leaf area to sap area conversion, intercept 
                                                    ! (sapwood area / leaf area) [cm2/m2]
     real(r8), allocatable :: allom_la_per_sa_slp(:)   ! Leaf area to sap area conversion, slope 
                                                    ! (sapwood area / leaf area / diameter) [cm2/m2/cm]
     real(r8), allocatable :: allom_l2fr(:)         ! Fine root biomass per leaf biomass ratio [kgC/kgC]
     real(r8), allocatable :: allom_agb_frac(:)     ! Fraction of stem above ground [-]
     real(r8), allocatable :: allom_d2h1(:)         ! Parameter 1 for d2h allometry (intercept, or "c")
     real(r8), allocatable :: allom_d2h2(:)         ! Parameter 2 for d2h allometry (slope, or "m")
     real(r8), allocatable :: allom_d2h3(:)         ! Parameter 3 for d2h allometry (optional)
     real(r8), allocatable :: allom_d2bl1(:)        ! Parameter 1 for d2bl allometry (intercept)
     real(r8), allocatable :: allom_d2bl2(:)        ! Parameter 2 for d2bl allometry (slope)
     real(r8), allocatable :: allom_d2bl3(:)           ! Parameter 3 for d2bl allometry (optional)
     real(r8), allocatable :: allom_sai_scaler(:)      ! 
     real(r8), allocatable :: allom_blca_expnt_diff(:) ! Any difference in the exponent between the leaf
                                                       ! biomass and crown area scaling
     real(r8), allocatable :: allom_d2ca_coefficient_max(:)  ! upper (savanna) value for crown 
                                                             ! area to dbh coefficient
     real(r8), allocatable :: allom_d2ca_coefficient_min(:)  ! lower (closed-canopy forest) value for crown 
                                                             ! area to dbh coefficient
     real(r8), allocatable :: allom_agb1(:)         ! Parameter 1 for agb allometry
     real(r8), allocatable :: allom_agb2(:)         ! Parameter 2 for agb allometry
     real(r8), allocatable :: allom_agb3(:)         ! Parameter 3 for agb allometry
     real(r8), allocatable :: allom_agb4(:)         ! Parameter 3 for agb allometry
     
     real(r8), allocatable :: allom_frbstor_repro(:) ! fraction of bstrore for reproduction after mortality

     ! Prescribed Physiology Mode Parameters
     real(r8), allocatable :: prescribed_npp_canopy(:)           ! this is only for the special 
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_npp_understory(:)       ! this is only for the special 
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_mortality_canopy(:)     ! this is only for the special
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_mortality_understory(:) ! this is only for the special 
                                                                 ! prescribed_physiology_mode
     real(r8), allocatable :: prescribed_recruitment(:)          ! this is only for the special 
                                                                 ! prescribed_physiology_mode

     
     ! Plant Reactive Transport (allocation)

     real(r8), allocatable :: grperc(:)                 ! Growth respiration per unit Carbon gained
                                                        ! One value for whole plant
                                                        ! ONLY parteh_mode == 1  [kg/kg]

     real(r8), allocatable :: prt_grperc_organ(:,:)     ! Unit growth respiration (pft x organ) [kg/kg]
                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                        ! THIS IS NOT READ IN BY THE PARAMETER FILE
                                                        ! THIS IS JUST FILLED BY GRPERC.  WE KEEP THIS
                                                        ! PARAMETER FOR HYPOTHESIS TESTING (ADVANCED USE)
                                                        ! IT HAS THE PRT_ TAG BECAUSE THIS PARAMETER
                                                        ! IS USED INSIDE PARTEH, WHILE GRPERC IS APPLIED
                                                        ! IN THE LEAF BIOPHYSICS SCHEME
                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     real(r8), allocatable :: prt_nitr_stoich_p1(:,:)   ! Parameter 1 for nitrogen stoichiometry (pft x organ) 
     real(r8), allocatable :: prt_nitr_stoich_p2(:,:)   ! Parameter 2 for nitrogen stoichiometry (pft x organ) 
     real(r8), allocatable :: prt_phos_stoich_p1(:,:)   ! Parameter 1 for phosphorous stoichiometry (pft x organ) 
     real(r8), allocatable :: prt_phos_stoich_p2(:,:)   ! Parameter 2 for phosphorous stoichiometry (pft x organ) 
     real(r8), allocatable :: prt_alloc_priority(:,:)   ! Allocation priority for each organ (pft x organ) [integer 0-6]

     ! Turnover related things

     real(r8), allocatable :: phenflush_fraction(:)       ! Maximum fraction of storage carbon used to flush leaves
                                                          ! on bud-burst [kgC/kgC]
     real(r8), allocatable :: phen_cold_size_threshold(:) ! stem/leaf drop occurs on DBH size of decidious non-woody 
                                                          !(coastal grass) plants larger than the threshold value 							  
     real(r8), allocatable :: phen_stem_drop_fraction(:)  ! Fraction of stem dropped/senescened for decidious 
                                                          ! non-woody (grass) plants						  

     real(r8), allocatable :: senleaf_long_fdrought(:)    ! Multiplication factor for leaf longevity of senescent 
                                                          ! leaves during drought( 1.0 indicates no change)

     real(r8), allocatable :: root_long(:)                ! root turnover time (longevity) (pft)             [yr]
     real(r8), allocatable :: branch_turnover(:)          ! Turnover time for branchfall on live trees (pft) [yr]
     real(r8), allocatable :: turnover_retrans_mode(:)    ! Retranslocation method (pft)

     real(r8), allocatable :: turnover_carb_retrans(:,:)  ! carbon re-translocation fraction (pft x organ)
     real(r8), allocatable :: turnover_nitr_retrans(:,:)  ! nitrogen re-translocation fraction (pft x organ)
     real(r8), allocatable :: turnover_phos_retrans(:,:)  ! phosphorous re-translocation fraction (pft x organ)

     ! Parameters dimensioned by PFT and leaf age
     real(r8), allocatable :: leaf_long(:,:)              ! Leaf turnover time (longevity) (pft x age-class)
                                                          ! If there is >1 class, it is the longevity from
                                                          ! one class to the next [yr]
     
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

     ! PFT x Organ Dimension  (organs are: 1=leaf, 2=stem, 3=transporting root, 4=absorbing root)
     real(r8), allocatable :: hydr_avuln_node(:,:)  ! xylem vulernability curve shape parameter 
     real(r8), allocatable :: hydr_p50_node(:,:)    ! xylem water potential at 50% conductivity loss (MPa)
     real(r8), allocatable :: hydr_thetas_node(:,:) ! saturated water content (cm3/cm3)
     real(r8), allocatable :: hydr_epsil_node(:,:)  ! bulk elastic modulus (MPa)
     real(r8), allocatable :: hydr_pitlp_node(:,:)  ! turgor loss point (MPa)
     real(r8), allocatable :: hydr_resid_node(:,:)  ! residual fraction (fraction)
     real(r8), allocatable :: hydr_fcap_node(:,:)   ! fraction of (1-resid_node) that is capillary in source
     real(r8), allocatable :: hydr_pinot_node(:,:)  ! osmotic potential at full turgor
     real(r8), allocatable :: hydr_kmax_node(:,:)   ! maximum xylem conductivity per unit conducting xylem area
     
   contains
     procedure, public :: Init => EDpftconInit
     procedure, public :: Register
     procedure, public :: Receive
     procedure, private :: Register_PFT
     procedure, private :: Receive_PFT
     procedure, private :: Register_PFT_nvariants
     procedure, private :: Receive_PFT_nvariants
     procedure, private :: Register_PFT_hydr_organs
     procedure, private :: Receive_PFT_hydr_organs 
     procedure, private :: Register_PFT_prt_organs
     procedure, private :: Receive_PFT_prt_organs
     procedure, private :: Register_PFT_leafage
     procedure, private :: Receive_PFT_leafage
     procedure, private :: Register_PFT_numrad
     procedure, private :: Receive_PFT_numrad
  end type EDPftvarcon_type

  type(EDPftvarcon_type), public :: EDPftvarcon_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FatesReportPFTParams
  public :: FatesCheckParams
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EDpftconInit(this)

    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this

  end subroutine EDpftconInit

  !-----------------------------------------------------------------------
  subroutine Register(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%Register_PFT(fates_params)
    call this%Register_PFT_numrad(fates_params)
    call this%Register_PFT_nvariants(fates_params)
    call this%Register_PFT_hydr_organs(fates_params)
    call this%Register_PFT_prt_organs(fates_params)
    call this%Register_PFT_leafage(fates_params)
    
  end subroutine Register

  !-----------------------------------------------------------------------
  subroutine Receive(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    call this%Receive_PFT(fates_params)
    call this%Receive_PFT_numrad(fates_params)
    call this%Receive_PFT_nvariants(fates_params)
    call this%Receive_PFT_hydr_organs(fates_params)
    call this%Receive_PFT_prt_organs(fates_params)
    call this%Receive_PFT_leafage(fates_params)

  end subroutine Receive

  !-----------------------------------------------------------------------
  subroutine Register_PFT(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)

    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
    !X!         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_dbh_repro_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_mort_freezetol'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_wood_density'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_recruit_hgt_min'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_stor_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fire_crown_depth_frac'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fire_bark_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fire_crown_kill'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_recruit_initd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_suppl'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_BB_slope'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_senleaf_long_fdrought'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_root_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_alloc_mature'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_alloc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_c2b'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_woody'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_stress_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_season_decid'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_evergreen'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_l2fr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_slamax'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_slatop'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_roota_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_rootb_par'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_lf_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_flab'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_fcel'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_fr_flig'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_xl'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)


    name = 'fates_leaf_clumping_index'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_c3psn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_smpso'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_smpsc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_grperc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_reduction_curvature'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_maintresp_reduction_intercept'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prescribed_npp_canopy'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prescribed_npp_understory'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prescribed_mortality_canopy'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prescribed_mortality_understory'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prescribed_recruitment'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
       
    name = 'fates_fire_alpha_SH'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_dbh_maxheight'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_hmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_lmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_fmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_amode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_stmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_cmode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_smode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_la_per_sa_int'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_la_per_sa_slp'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb_frac'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2h3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_allom_d2bl1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2bl2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2bl3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_blca_expnt_diff'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2ca_coefficient_max'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_d2ca_coefficient_min'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_sai_scaler'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb3'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_allom_agb4'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
	 
    name = 'fates_allom_frbstor_repro'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_p_taper'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_rs2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_srl'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_rfrac_stem'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_hydr_avuln_gs'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_hydr_p50_gs'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_mort_bmort'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
	  
    name = 'fates_mort_scalar_coldstress'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_mort_scalar_cstarvation'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_mort_scalar_hydrfailure'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_mort_hf_sm_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
	 
    name = 'fates_mort_hf_flc_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
	 
    name = 'fates_leaf_vcmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_tpuha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_vcmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_tpuhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_vcmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_tpuse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_germination_timescale'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_seed_decay_turnover'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_retrans_mode'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_branch_turnover'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_trim_limit'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_trim_inc'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_leaf_diameter'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_z0mr'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_displar'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phenflush_fraction'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
	  
    name = 'fates_phen_cold_size_threshold'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_phen_stem_drop_fraction'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
  end subroutine Register_PFT

  !-----------------------------------------------------------------------
  subroutine Receive_PFT(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    name = 'fates_seed_dbh_repro_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dbh_repro_threshold)

    name = 'fates_mort_freezetol'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%freezetol)

    name = 'fates_wood_density'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%wood_density)

    name = 'fates_recruit_hgt_min'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hgt_min)

    name = 'fates_alloc_storage_cushion'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%cushion)

    name = 'fates_leaf_stor_priority'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%leaf_stor_priority)

    name = 'fates_fire_crown_depth_frac'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown)

    name = 'fates_fire_bark_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%bark_scaler)

    name = 'fates_fire_crown_kill'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%crown_kill)

    name = 'fates_recruit_initd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%initd)

    name = 'fates_seed_suppl'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_suppl)

    name = 'fates_leaf_BB_slope'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%BB_slope)

    name = 'fates_senleaf_long_fdrought'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%senleaf_long_fdrought)

    name = 'fates_root_long'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%root_long)

    name = 'fates_seed_alloc_mature'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_alloc_mature)

    name = 'fates_seed_alloc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_alloc)

    name = 'fates_c2b'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%c2b)

    name = 'fates_woody'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%woody)

    name = 'fates_phen_stress_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%stress_decid)

    name = 'fates_phen_season_decid'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%season_decid)

    name = 'fates_phen_evergreen'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%evergreen)

    name = 'fates_leaf_slamax'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%slamax)

    name = 'fates_leaf_slatop'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%slatop)

    name = 'fates_roota_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%roota_par)

    name = 'fates_rootb_par'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootb_par)

    name = 'fates_lf_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flab)

    name = 'fates_lf_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_fcel)

    name = 'fates_lf_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%lf_flig)

    name = 'fates_fr_flab'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flab)

    name = 'fates_fr_fcel'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_fcel)

    name = 'fates_fr_flig'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fr_flig)

    name = 'fates_leaf_xl'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%xl)

    name = 'fates_leaf_clumping_index'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%clumping_index)

    name = 'fates_leaf_c3psn'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%c3psn)

    name = 'fates_smpso'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpso)

    name = 'fates_smpsc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%smpsc)

    name = 'fates_grperc'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%grperc)

    name = 'fates_maintresp_reduction_curvature'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%maintresp_reduction_curvature)

    name = 'fates_maintresp_reduction_intercept'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%maintresp_reduction_intercept)

    name = 'fates_prescribed_npp_canopy'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%prescribed_npp_canopy)

    name = 'fates_prescribed_npp_understory'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%prescribed_npp_understory)

    name = 'fates_prescribed_mortality_canopy'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%prescribed_mortality_canopy)

    name = 'fates_prescribed_mortality_understory'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%prescribed_mortality_understory)

    name = 'fates_prescribed_recruitment'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%prescribed_recruitment)

    name = 'fates_fire_alpha_SH'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%fire_alpha_SH)

    name = 'fates_allom_dbh_maxheight'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%allom_dbh_maxheight)

    name = 'fates_allom_hmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_hmode)

    name = 'fates_allom_lmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_lmode)

    name = 'fates_allom_fmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_fmode)

    name = 'fates_allom_amode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_amode)

    name = 'fates_allom_stmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_stmode)

    name = 'fates_allom_cmode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_cmode)

    name = 'fates_allom_smode'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_smode)

    name = 'fates_allom_la_per_sa_int'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_la_per_sa_int)

    name = 'fates_allom_la_per_sa_slp'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_la_per_sa_slp)

    name = 'fates_allom_l2fr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_l2fr)

    name = 'fates_allom_agb_frac'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_agb_frac)

    name = 'fates_allom_d2h1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2h1)

    name = 'fates_allom_d2h2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2h2)

    name = 'fates_allom_d2h3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2h3)

    name = 'fates_allom_d2bl1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2bl1)

    name = 'fates_allom_d2bl2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2bl2)

    name = 'fates_allom_d2bl3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2bl3)

    name = 'fates_allom_blca_expnt_diff'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_blca_expnt_diff)

    name = 'fates_allom_d2ca_coefficient_max'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2ca_coefficient_max)

    name = 'fates_allom_d2ca_coefficient_min'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_d2ca_coefficient_min)

    name = 'fates_allom_sai_scaler'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_sai_scaler)

    name = 'fates_allom_agb1'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_agb1)

    name = 'fates_allom_agb2'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_agb2)

    name = 'fates_allom_agb3'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_agb3)

    name = 'fates_allom_agb4'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_agb4)
	 
    name = 'fates_allom_frbstor_repro'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%allom_frbstor_repro) 

    name = 'fates_hydr_p_taper'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_p_taper)

    name = 'fates_hydr_rs2'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_rs2)
    
    name = 'fates_hydr_srl'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_srl)
    
    name = 'fates_hydr_rfrac_stem'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_rfrac_stem)

    name = 'fates_hydr_avuln_gs'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_avuln_gs)
    
    name = 'fates_hydr_p50_gs'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%hydr_p50_gs)

    name = 'fates_mort_bmort'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%bmort)
	 	 
    name = 'fates_mort_scalar_coldstress'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%mort_scalar_coldstress)

    name = 'fates_mort_scalar_cstarvation'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%mort_scalar_cstarvation)

    name = 'fates_mort_scalar_hydrfailure'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%mort_scalar_hydrfailure)

    name = 'fates_mort_hf_sm_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hf_sm_threshold)
	 
    name = 'fates_mort_hf_flc_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%hf_flc_threshold)
	 
    name = 'fates_leaf_vcmaxha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxha)

    name = 'fates_leaf_jmaxha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxha)

    name = 'fates_leaf_tpuha'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuha)

    name = 'fates_leaf_vcmaxhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxhd)

    name = 'fates_leaf_jmaxhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxhd)

    name = 'fates_leaf_tpuhd'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuhd)

    name = 'fates_leaf_vcmaxse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%vcmaxse)

    name = 'fates_leaf_jmaxse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%jmaxse)

    name = 'fates_leaf_tpuse'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%tpuse)

    name = 'fates_seed_germination_timescale'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%germination_timescale)

    name = 'fates_seed_decay_turnover'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%seed_decay_turnover)

    name = 'fates_branch_turnover'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%branch_turnover)
    
    name = 'fates_turnover_retrans_mode'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%turnover_retrans_mode)

    name = 'fates_trim_limit'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%trim_limit)

    name = 'fates_trim_inc'
    call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%trim_inc)

    name = 'fates_leaf_diameter'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%dleaf)

    name = 'fates_z0mr'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%z0mr)

    name = 'fates_displar'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%displar)

    name = 'fates_phenflush_fraction'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%phenflush_fraction)
	  	 
    name = 'fates_phen_cold_size_threshold'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%phen_cold_size_threshold)
	 
    name = 'fates_phen_stem_drop_fraction'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%phen_stem_drop_fraction)	 
	 

  end subroutine Receive_PFT

  !-----------------------------------------------------------------------
  subroutine Register_PFT_numrad(this, fates_params)
    ! NOTE(bja, 2017-02) these are 2-d parameters, but they are
    ! currently stored in the parameter file as separate 1-d
    ! arrays. We have to register the parameters as 1-d arrays as they
    ! are on the parameter file. We store them as 2-d in the receive step.
    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)
    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)
    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
    !X!         dimension_names=dim_names)

    name = 'fates_rholvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rholnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rhosvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_rhosnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_taulvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_taulnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_tausvis'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    name = 'fates_tausnir'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)


  end subroutine Register_PFT_numrad

  !-----------------------------------------------------------------------
  subroutine Receive_PFT_numrad(this, fates_params)
    ! NOTE(bja, 2017-02) these are 2-d parameters, but they are
    ! currently stored in the parameter file as separate 1-d arrays.
    ! We can't allocate slices of arrays separately, so we have to
    ! manually allocate the memory here, retreive into a dummy array,
    ! and copy. All parameters in this subroutine are sized the same,
    ! so we can reused the dummy array. If someone wants to cleanup
    ! the input file, all this complexity can be removed.
    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length, max_dimensions

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    integer :: index
    integer :: dimension_shape
    integer :: dimension_sizes(max_dimensions)
    character(len=param_string_length) :: dimension_names(max_dimensions)
    logical :: is_host_param

    integer :: lower_bound_1, upper_bound_1, lower_bound_2, upper_bound_2
    real(r8), allocatable :: dummy_data(:)

    ! Fetch metadata from a representative variable. All variables
    ! called by this subroutine must be dimensioned the same way!
    name = 'fates_rholvis'
    index = fates_params%FindIndex(name)
    call fates_params%GetMetaData(index, name, dimension_shape, dimension_sizes, dimension_names, is_host_param)
    lower_bound_1 = lower_bound_pft
    upper_bound_1 = lower_bound_pft + dimension_sizes(1) - 1
    lower_bound_2 = lower_bound_general
    upper_bound_2 = maxSWb      ! When we have radiation parameters read in as a vector
                                ! We will compare the vector dimension size that we
                                ! read-in to the parameterized size that fates expects

    allocate(dummy_data(lower_bound_1:upper_bound_1))

    !
    ! received rhol data
    !
    allocate(this%rhol(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_rholvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_rholnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhol(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received rhos data
    !
    allocate(this%rhos(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_rhosvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_rhosnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%rhos(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taul data
    !
    allocate(this%taul(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_taulvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_taulnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taul(lower_bound_1:upper_bound_1, inir) = dummy_data

    !
    ! received taus data
    !
    allocate(this%taus(lower_bound_1:upper_bound_1, lower_bound_2:upper_bound_2))
    
    name = 'fates_tausvis'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taus(lower_bound_1:upper_bound_1, ivis) = dummy_data

    name = 'fates_tausnir'
    call fates_params%RetreiveParameter(name=name, &
         data=dummy_data)
    this%taus(lower_bound_1:upper_bound_1, inir) = dummy_data

  end subroutine Receive_PFT_numrad

  !-----------------------------------------------------------------------
  subroutine Register_PFT_nvariants(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_variants, dimension_name_pft, dimension_shape_2d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    ! NOTE(bja, 2017-01) initialization doesn't seem to work correctly
    ! if dim_names has a parameter qualifier.
    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_variants

    !X!    name = ''
    !X!    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
    !X!         dimension_names=dim_names)

    name = 'fates_rootprof_beta'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

  end subroutine Register_PFT_nvariants

  !-----------------------------------------------------------------------
  subroutine Receive_PFT_nvariants(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length) :: name

    !X!    name = ''
    !X!    call fates_params%RetreiveParameter(name=name, &
    !X!         data=this%)

    name = 'fates_rootprof_beta'
    call fates_params%RetreiveParameterAllocate(name=name, &
         data=this%rootprof_beta)

  end subroutine Receive_PFT_nvariants

  ! -----------------------------------------------------------------------

  subroutine Register_PFT_leafage(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_leaf_age
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_2d
    
    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_leaf_age

    name = 'fates_leaf_long'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_leaf_vcmax25top'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    

    return
  end subroutine Register_PFT_leafage

  ! =====================================================================================

  subroutine Receive_PFT_leafage(this, fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(EDPftvarcon_type), intent(inout) :: this
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name

     name = 'fates_leaf_long'
     call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%leaf_long)
     
     name = 'fates_leaf_vcmax25top'
     call fates_params%RetreiveParameterAllocate(name=name, &
          data=this%vcmax25top)

     return
   end subroutine Receive_PFT_leafage

  ! =====================================================================================

  
  subroutine Register_PFT_prt_organs(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_prt_organs
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_2d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    ! NOTE(bja, 2017-01) initialization doesn't seem to work correctly
    ! if dim_names has a parameter qualifier.
    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_prt_organs

    name = 'fates_prt_nitr_stoich_p1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_nitr_stoich_p2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_phos_stoich_p1'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_phos_stoich_p2'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_prt_alloc_priority'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_turnover_carb_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)    

    name = 'fates_turnover_nitr_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_turnover_phos_retrans'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    
 end subroutine Register_PFT_prt_organs

  ! =====================================================================================

  subroutine Receive_PFT_prt_organs(this, fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(EDPftvarcon_type), intent(inout) :: this
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name

     name = 'fates_prt_nitr_stoich_p1'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%prt_nitr_stoich_p1)

     name = 'fates_prt_nitr_stoich_p2'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%prt_nitr_stoich_p2)
     
     name = 'fates_prt_phos_stoich_p1'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%prt_phos_stoich_p1)

     name = 'fates_prt_phos_stoich_p2'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%prt_phos_stoich_p2)
    
     name = 'fates_prt_alloc_priority'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%prt_alloc_priority)

     name = 'fates_turnover_carb_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%turnover_carb_retrans)

     name = 'fates_turnover_nitr_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%turnover_nitr_retrans)

     name = 'fates_turnover_phos_retrans'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%turnover_phos_retrans)

  end subroutine Receive_PFT_prt_organs

  ! -----------------------------------------------------------------------
  
  subroutine Register_PFT_hydr_organs(this, fates_params)

    use FatesParametersInterface, only : fates_parameters_type, param_string_length
    use FatesParametersInterface, only : max_dimensions, dimension_name_hydr_organs
    use FatesParametersInterface, only : dimension_name_pft, dimension_shape_2d

    implicit none

    class(EDPftvarcon_type), intent(inout) :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    integer, parameter :: dim_lower_bound(2) = (/ lower_bound_pft, lower_bound_general /)
    character(len=param_string_length) :: dim_names(2)
    character(len=param_string_length) :: name

    ! NOTE(bja, 2017-01) initialization doesn't seem to work correctly
    ! if dim_names has a parameter qualifier.
    dim_names(1) = dimension_name_pft
    dim_names(2) = dimension_name_hydr_organs

    name = 'fates_hydr_avuln_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_hydr_p50_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_thetas_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_epsil_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_pitlp_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_resid_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_fcap_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_hydr_pinot_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_hydr_kmax_node'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_2d, &
          dimension_names=dim_names, lower_bounds=dim_lower_bound)
    

  end subroutine Register_PFT_hydr_organs

  !-----------------------------------------------------------------------

  subroutine Receive_PFT_hydr_organs(this, fates_params)
     
     use FatesParametersInterface, only : fates_parameters_type
     use FatesParametersInterface, only : param_string_length
     
     implicit none
     
     class(EDPftvarcon_type), intent(inout) :: this
     class(fates_parameters_type), intent(inout) :: fates_params
     
     character(len=param_string_length) :: name
     
     name = 'fates_hydr_avuln_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_avuln_node)

     name = 'fates_hydr_p50_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_p50_node)

     name = 'fates_hydr_thetas_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_thetas_node)
     
     name = 'fates_hydr_epsil_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_epsil_node)
     
     name = 'fates_hydr_pitlp_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_pitlp_node)
     
     name = 'fates_hydr_resid_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_resid_node)
     
     name = 'fates_hydr_fcap_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_fcap_node)

     name = 'fates_hydr_pinot_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_pinot_node)

     name = 'fates_hydr_kmax_node'
     call fates_params%RetreiveParameterAllocate(name=name, &
           data=this%hydr_kmax_node)

  end subroutine Receive_PFT_hydr_organs

  ! ===============================================================================================
  
  subroutine FatesReportPFTParams(is_master)
     
     ! Argument
     logical, intent(in) :: is_master  ! Only log if this is the master proc

     logical, parameter :: debug_report = .false.
     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft,ipft
     
     npft = size(EDPftvarcon_inst%evergreen,1)
     
     if(debug_report .and. is_master) then
        
        if(npft>100)then
           write(fates_log(),*) 'you are trying to report pft parameters during initialization'
           write(fates_log(),*) 'but you have so many that it is over-running the format spec'
           write(fates_log(),*) 'simply bump up the muptiplier in parameter fmt0 shown above'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        write(fates_log(),*) '-----------  FATES PFT Parameters -----------------'
        write(fates_log(),fmt0) 'dbh max height = ',EDPftvarcon_inst%allom_dbh_maxheight
        write(fates_log(),fmt0) 'dbh mature = ',EDPftvarcon_inst%dbh_repro_threshold
        write(fates_log(),fmt0) 'freezetol = ',EDPftvarcon_inst%freezetol
        write(fates_log(),fmt0) 'wood_density = ',EDPftvarcon_inst%wood_density
        write(fates_log(),fmt0) 'hgt_min = ',EDPftvarcon_inst%hgt_min
        write(fates_log(),fmt0) 'dleaf = ',EDPftvarcon_inst%dleaf
        write(fates_log(),fmt0) 'z0mr = ',EDPftvarcon_inst%z0mr
        write(fates_log(),fmt0) 'displar = ',EDPftvarcon_inst%displar
        write(fates_log(),fmt0) 'cushion = ',EDPftvarcon_inst%cushion
        write(fates_log(),fmt0) 'leaf_stor_priority = ',EDPftvarcon_inst%leaf_stor_priority
        write(fates_log(),fmt0) 'crown = ',EDPftvarcon_inst%crown
        write(fates_log(),fmt0) 'bark_scaler = ',EDPftvarcon_inst%bark_scaler
        write(fates_log(),fmt0) 'crown_kill = ',EDPftvarcon_inst%crown_kill
        write(fates_log(),fmt0) 'initd = ',EDPftvarcon_inst%initd
        write(fates_log(),fmt0) 'seed_suppl = ',EDPftvarcon_inst%seed_suppl
        write(fates_log(),fmt0) 'BB_slope = ',EDPftvarcon_inst%BB_slope
        write(fates_log(),fmt0) 'root_long = ',EDPftvarcon_inst%root_long
        write(fates_log(),fmt0) 'senleaf_long_fdrought = ',EDPftvarcon_inst%senleaf_long_fdrought
        write(fates_log(),fmt0) 'seed_alloc_mature = ',EDPftvarcon_inst%seed_alloc_mature
        write(fates_log(),fmt0) 'seed_alloc = ',EDPftvarcon_inst%seed_alloc
        write(fates_log(),fmt0) 'woody = ',EDPftvarcon_inst%woody
        write(fates_log(),fmt0) 'stress_decid = ',EDPftvarcon_inst%stress_decid
        write(fates_log(),fmt0) 'season_decid = ',EDPftvarcon_inst%season_decid
        write(fates_log(),fmt0) 'evergreen = ',EDPftvarcon_inst%evergreen
        write(fates_log(),fmt0) 'slamax = ',EDPftvarcon_inst%slamax
        write(fates_log(),fmt0) 'slatop = ',EDPftvarcon_inst%slatop        
        write(fates_log(),fmt0) 'leaf_long = ',EDPftvarcon_inst%leaf_long
        write(fates_log(),fmt0) 'roota_par = ',EDPftvarcon_inst%roota_par
        write(fates_log(),fmt0) 'rootb_par = ',EDPftvarcon_inst%rootb_par
        write(fates_log(),fmt0) 'lf_flab = ',EDPftvarcon_inst%lf_flab
        write(fates_log(),fmt0) 'lf_fcel = ',EDPftvarcon_inst%lf_fcel
        write(fates_log(),fmt0) 'lf_flig = ',EDPftvarcon_inst%lf_flig
        write(fates_log(),fmt0) 'fr_flab = ',EDPftvarcon_inst%fr_flab
        write(fates_log(),fmt0) 'fr_fcel = ',EDPftvarcon_inst%fr_fcel
        write(fates_log(),fmt0) 'fr_flig = ',EDPftvarcon_inst%fr_flig
        write(fates_log(),fmt0) 'xl = ',EDPftvarcon_inst%xl
        write(fates_log(),fmt0) 'clumping_index = ',EDPftvarcon_inst%clumping_index
        write(fates_log(),fmt0) 'c3psn = ',EDPftvarcon_inst%c3psn
        write(fates_log(),fmt0) 'vcmax25top = ',EDPftvarcon_inst%vcmax25top
        write(fates_log(),fmt0) 'smpso = ',EDPftvarcon_inst%smpso
        write(fates_log(),fmt0) 'smpsc = ',EDPftvarcon_inst%smpsc
        write(fates_log(),fmt0) 'grperc = ',EDPftvarcon_inst%grperc
        write(fates_log(),fmt0) 'c2b = ',EDPftvarcon_inst%c2b
        write(fates_log(),fmt0) 'bmort = ',EDPftvarcon_inst%bmort
	write(fates_log(),fmt0) 'mort_scalar_coldstress = ',EDPftvarcon_inst%mort_scalar_coldstress
        write(fates_log(),fmt0) 'mort_scalar_cstarvation = ',EDPftvarcon_inst%mort_scalar_cstarvation
        write(fates_log(),fmt0) 'mort_scalar_hydrfailure = ',EDPftvarcon_inst%mort_scalar_hydrfailure
        write(fates_log(),fmt0) 'hf_sm_threshold = ',EDPftvarcon_inst%hf_sm_threshold
        write(fates_log(),fmt0) 'hf_flc_threshold = ',EDPftvarcon_inst%hf_flc_threshold
        write(fates_log(),fmt0) 'vcmaxha = ',EDPftvarcon_inst%vcmaxha
        write(fates_log(),fmt0) 'jmaxha = ',EDPftvarcon_inst%jmaxha
        write(fates_log(),fmt0) 'tpuha = ',EDPftvarcon_inst%tpuha
        write(fates_log(),fmt0) 'vcmaxhd = ',EDPftvarcon_inst%vcmaxhd
        write(fates_log(),fmt0) 'jmaxhd = ',EDPftvarcon_inst%jmaxhd
        write(fates_log(),fmt0) 'tpuhd = ',EDPftvarcon_inst%tpuhd
        write(fates_log(),fmt0) 'vcmaxse = ',EDPftvarcon_inst%vcmaxse
        write(fates_log(),fmt0) 'jmaxse = ',EDPftvarcon_inst%jmaxse
        write(fates_log(),fmt0) 'tpuse = ',EDPftvarcon_inst%tpuse
        write(fates_log(),fmt0) 'germination_timescale = ',EDPftvarcon_inst%germination_timescale
        write(fates_log(),fmt0) 'seed_decay_turnover = ',EDPftvarcon_inst%seed_decay_turnover
        write(fates_log(),fmt0) 'branch_turnover = ',EDPftvarcon_inst%branch_turnover
        write(fates_log(),fmt0) 'trim_limit = ',EDPftvarcon_inst%trim_limit
        write(fates_log(),fmt0) 'trim_inc = ',EDPftvarcon_inst%trim_inc
        write(fates_log(),fmt0) 'rhol = ',EDPftvarcon_inst%rhol
        write(fates_log(),fmt0) 'rhos = ',EDPftvarcon_inst%rhos
        write(fates_log(),fmt0) 'taul = ',EDPftvarcon_inst%taul 
        write(fates_log(),fmt0) 'taus = ',EDPftvarcon_inst%taus
        write(fates_log(),fmt0) 'phenflush_fraction',EDpftvarcon_inst%phenflush_fraction
	write(fates_log(),fmt0) 'phen_cold_size_threshold = ',EDPftvarcon_inst%phen_cold_size_threshold
        write(fates_log(),fmt0) 'phen_stem_drop_fraction',EDpftvarcon_inst%phen_stem_drop_fraction
        write(fates_log(),fmt0) 'rootprof_beta = ',EDPftvarcon_inst%rootprof_beta
        write(fates_log(),fmt0) 'fire_alpha_SH = ',EDPftvarcon_inst%fire_alpha_SH
        write(fates_log(),fmt0) 'allom_hmode = ',EDPftvarcon_inst%allom_hmode
        write(fates_log(),fmt0) 'allom_lmode = ',EDPftvarcon_inst%allom_lmode
        write(fates_log(),fmt0) 'allom_fmode = ',EDPftvarcon_inst%allom_fmode
        write(fates_log(),fmt0) 'allom_amode = ',EDPftvarcon_inst%allom_amode
        write(fates_log(),fmt0) 'allom_cmode = ',EDPftvarcon_inst%allom_cmode
        write(fates_log(),fmt0) 'allom_smode = ',EDPftvarcon_inst%allom_smode
        write(fates_log(),fmt0) 'allom_la_per_sa_int = ',EDPftvarcon_inst%allom_la_per_sa_int
        write(fates_log(),fmt0) 'allom_la_per_sa_slp = ',EDPftvarcon_inst%allom_la_per_sa_slp
        write(fates_log(),fmt0) 'allom_l2fr = ',EDPftvarcon_inst%allom_l2fr
        write(fates_log(),fmt0) 'allom_agb_frac = ',EDPftvarcon_inst%allom_agb_frac
        write(fates_log(),fmt0) 'allom_d2h1 = ',EDPftvarcon_inst%allom_d2h1
        write(fates_log(),fmt0) 'allom_d2h2 = ',EDPftvarcon_inst%allom_d2h2
        write(fates_log(),fmt0) 'allom_d2h3 = ',EDPftvarcon_inst%allom_d2h3
        write(fates_log(),fmt0) 'allom_d2bl1 = ',EDPftvarcon_inst%allom_d2bl1
        write(fates_log(),fmt0) 'allom_d2bl2 = ',EDPftvarcon_inst%allom_d2bl2
        write(fates_log(),fmt0) 'allom_d2bl3 = ',EDPftvarcon_inst%allom_d2bl3
        write(fates_log(),fmt0) 'allom_sai_scaler = ',EDPftvarcon_inst%allom_sai_scaler
        write(fates_log(),fmt0) 'allom_blca_expnt_diff = ',EDPftvarcon_inst%allom_blca_expnt_diff
        write(fates_log(),fmt0) 'allom_d2ca_coefficient_max = ',EDPftvarcon_inst%allom_d2ca_coefficient_max
        write(fates_log(),fmt0) 'allom_d2ca_coefficient_min = ',EDPftvarcon_inst%allom_d2ca_coefficient_min        
        write(fates_log(),fmt0) 'allom_agb1 = ',EDPftvarcon_inst%allom_agb1
        write(fates_log(),fmt0) 'allom_agb2 = ',EDPftvarcon_inst%allom_agb2
        write(fates_log(),fmt0) 'allom_agb3 = ',EDPftvarcon_inst%allom_agb3
        write(fates_log(),fmt0) 'allom_agb4 = ',EDPftvarcon_inst%allom_agb4
	write(fates_log(),fmt0) 'allom_frbstor_repro = ',EDPftvarcon_inst%allom_frbstor_repro
        write(fates_log(),fmt0) 'hydr_p_taper = ',EDPftvarcon_inst%hydr_p_taper
        write(fates_log(),fmt0) 'hydr_rs2 = ',EDPftvarcon_inst%hydr_rs2
        write(fates_log(),fmt0) 'hydr_srl = ',EDPftvarcon_inst%hydr_srl
        write(fates_log(),fmt0) 'hydr_rfrac_stem = ',EDPftvarcon_inst%hydr_rfrac_stem
        write(fates_log(),fmt0) 'hydr_avuln_gs = ',EDPftvarcon_inst%hydr_avuln_gs
        write(fates_log(),fmt0) 'hydr_p50_gs = ',EDPftvarcon_inst%hydr_p50_gs
        write(fates_log(),fmt0) 'hydr_avuln_node = ',EDPftvarcon_inst%hydr_avuln_node
        write(fates_log(),fmt0) 'hydr_p50_node = ',EDPftvarcon_inst%hydr_p50_node
        write(fates_log(),fmt0) 'hydr_thetas_node = ',EDPftvarcon_inst%hydr_thetas_node 
        write(fates_log(),fmt0) 'hydr_epsil_node = ',EDPftvarcon_inst%hydr_epsil_node
        write(fates_log(),fmt0) 'hydr_pitlp_node = ',EDPftvarcon_inst%hydr_pitlp_node
        write(fates_log(),fmt0) 'hydr_resid_node = ',EDPftvarcon_inst%hydr_resid_node
        write(fates_log(),fmt0) 'hydr_fcap_node = ',EDPftvarcon_inst%hydr_fcap_node
        write(fates_log(),fmt0) 'hydr_pinot_node = ',EDPftvarcon_inst%hydr_pinot_node
        write(fates_log(),fmt0) 'hydr_kmax_node = ',EDPftvarcon_inst%hydr_kmax_node
        
        
        write(fates_log(),fmt0) 'prt_nitr_stoich_p1 = ',EDPftvarcon_inst%prt_nitr_stoich_p1
        write(fates_log(),fmt0) 'prt_nitr_stoich_p2 = ',EDPftvarcon_inst%prt_nitr_stoich_p2
        write(fates_log(),fmt0) 'prt_phos_stoich_p1 = ',EDPftvarcon_inst%prt_phos_stoich_p1
        write(fates_log(),fmt0) 'prt_phos_stoich_p2 = ',EDPftvarcon_inst%prt_phos_stoich_p2
        write(fates_log(),fmt0) 'prt_grperc_organ   = ',EDPftvarcon_inst%prt_grperc_organ
        write(fates_log(),fmt0) 'prt_alloc_priority = ',EDPftvarcon_inst%prt_alloc_priority

        write(fates_log(),fmt0) 'turnover_carb_retrans = ',EDPftvarcon_inst%turnover_carb_retrans
        write(fates_log(),fmt0) 'turnover_nitr_retrans = ',EDPftvarcon_inst%turnover_nitr_retrans
        write(fates_log(),fmt0) 'turnover_phos_retrans = ',EDPftvarcon_inst%turnover_phos_retrans

        write(fates_log(),*) '-------------------------------------------------'

     end if

  end subroutine FatesReportPFTParams


  ! =====================================================================================

  subroutine FatesCheckParams(is_master, parteh_mode)

     ! ----------------------------------------------------------------------------------
     !
     ! This subroutine performs logical checks on user supplied parameters.  It cross
     ! compares various parameters and will fail if they don't make sense.  
     ! Examples:
     ! A tree can not be defined as both evergreen and deciduous.  A woody plant
     ! cannot have a structural biomass allometry intercept of 0, and a non-woody
     ! plant (grass) can't have a non-zero intercept...
     ! -----------------------------------------------------------------------------------


     ! Argument
     logical, intent(in) :: is_master    ! Only log if this is the master proc
     integer, intent(in) :: parteh_mode  ! argument for nl flag hlm_parteh_mode

     character(len=32),parameter :: fmt0 = '(a,100(F12.4,1X))'

     integer :: npft     ! number of PFTs
     integer :: ipft     ! pft index
     integer :: nleafage ! size of the leaf age class array
     integer :: iage     ! leaf age class index
     integer :: norgans  ! size of the plant organ dimension

     npft = size(EDPftvarcon_inst%evergreen,1)

     ! Prior to performing checks copy grperc to the 
     ! organ dimensioned version

     norgans = size(EDPftvarcon_inst%prt_nitr_stoich_p1,2)
     allocate(EDPftvarcon_inst%prt_grperc_organ(npft,norgans))
     do ipft = 1,npft
        EDPftvarcon_inst%prt_grperc_organ(ipft,1:norgans) = EDPftvarcon_inst%grperc(ipft)
     end do

     
     if(.not.is_master) return

     if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
        write(fates_log(),*) 'FATES Plant Allocation and Reactive Transport'
        write(fates_log(),*) 'with flexible target stoichiometry for NP and'
        write(fates_log(),*) 'allometrically constrianed C is still under development'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        
     elseif (parteh_mode .ne. prt_carbon_allom_hyp) then
        
        write(fates_log(),*) 'FATES Plant Allocation and Reactive Transport has'
        write(fates_log(),*) 'only 1 module supported, allometric carbon only.'
        write(fates_log(),*) 'fates_parteh_mode must be set to 1 in the namelist'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if


     if (norgans .ne. num_organ_types) then
        write(fates_log(),*) 'The size of the organ dimension for PRT parameters'
        write(fates_log(),*) 'as specified in the parameter file is incompatible.'
        write(fates_log(),*) 'All currently acceptable hypothesese are using'
        write(fates_log(),*) 'the full set of num_organ_types = ',num_organ_types
        write(fates_log(),*) 'The parameter file listed ',norgans
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     
     do ipft = 1,npft
        
        ! Check to see if evergreen, deciduous flags are mutually exclusive
        ! ----------------------------------------------------------------------------------

        if ( int(EDPftvarcon_inst%evergreen(ipft) +    &
                 EDPftvarcon_inst%season_decid(ipft) + &
                 EDPftvarcon_inst%stress_decid(ipft)) .ne. 1 ) then
           
           write(fates_log(),*) 'PFT # ',ipft,' must be defined as having one of three'
           write(fates_log(),*) 'phenology habits, ie == 1'
           write(fates_log(),*) 'stress_decid: ',EDPftvarcon_inst%stress_decid(ipft)
           write(fates_log(),*) 'season_decid: ',EDPftvarcon_inst%season_decid(ipft)
           write(fates_log(),*) 'evergreen: ',EDPftvarcon_inst%evergreen(ipft)
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if


        ! Check to see if mature and base seed allocation is greater than 1
        ! ----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%seed_alloc(ipft) + &
               EDPftvarcon_inst%seed_alloc_mature(ipft)) > 1.0_r8 ) then

           write(fates_log(),*) 'The sum of seed allocation from base and mature trees may'
           write(fates_log(),*) ' not exceed 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' seed_alloc: ',EDPftvarcon_inst%seed_alloc(ipft)
           write(fates_log(),*) ' seed_alloc_mature: ',EDPftvarcon_inst%seed_alloc_mature(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if

        ! Check if woody plants have a structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
        if ( ( EDPftvarcon_inst%allom_agb1(ipft) <= tiny(EDPftvarcon_inst%allom_agb1(ipft)) ) .and. &
             ( int(EDPftvarcon_inst%woody(ipft)) .eq. 1 ) ) then

           write(fates_log(),*) 'Woody plants are expected to have a non-zero intercept'
           write(fates_log(),*) ' in the diameter to AGB allometry equations'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' allom_agb1: ',EDPftvarcon_inst%allom_agb1(ipft)
           write(fates_log(),*) ' woody: ',int(EDPftvarcon_inst%woody(ipft))
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if

        ! Check if non-woody plants have structural biomass (agb) intercept
        ! ----------------------------------------------------------------------------------
!        if ( ( EDPftvarcon_inst%allom_agb1(ipft) > tiny(EDPftvarcon_inst%allom_agb1(ipft)) ) .and. &
!              ( int(EDPftvarcon_inst%woody(ipft)) .ne. 1 ) ) then
!
!           write(fates_log(),*) 'Non-woody plants are expected to have a zero intercept'
!           write(fates_log(),*) ' in the diameter to AGB allometry equations'
!           write(fates_log(),*) ' This is because the definition of AGB (as far as allometry)'
!           write(fates_log(),*) ' is concerned, ignores leaf and fine-roots, and only contains'
!           write(fates_log(),*) ' woody tissues (sap and structural dead wood).'
!           write(fates_log(),*) ' PFT#: ',ipft
!           write(fates_log(),*) ' allom_agb1: ',EDPftvarcon_inst%allom_agb1(ipft)
!           write(fates_log(),*) ' woody: ',int(EDPftvarcon_inst%woody(ipft))
!           write(fates_log(),*) ' Aborting'
!           call endrun(msg=errMsg(sourcefile, __LINE__))
!
!        end if


        ! Check if the fraction of storage used for flushing deciduous trees
        ! is greater than zero, and less than or equal to 1.

        if ( int(EDPftvarcon_inst%evergreen(ipft)) .ne. 1 ) then 
           if ( ( EDPftvarcon_inst%phenflush_fraction(ipft) < nearzero ) .or. &
                ( EDPFtvarcon_inst%phenflush_fraction(ipft) > 1 ) ) then
              
              write(fates_log(),*) ' Deciduous plants must flush some storage carbon'
              write(fates_log(),*) ' on bud-burst. If phenflush_fraction is not greater than 0'
              write(fates_log(),*) ' it will not be able to put out any leaves. Plants need leaves.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' evergreen flag: (shold be 0):',int(EDPftvarcon_inst%evergreen(ipft))
              write(fates_log(),*) ' phenflush_fraction: ', EDPFtvarcon_inst%phenflush_fraction(ipft)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           if ( ( EDPftvarcon_inst%phen_stem_drop_fraction(ipft) < 0.0_r8 ) .or. &
                ( EDPFtvarcon_inst%phen_stem_drop_fraction(ipft) > 1 ) ) then
              
              write(fates_log(),*) ' Deciduous non-wood plants must keep 0-100% of their stems'
              write(fates_log(),*) ' during the dedicous period.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' evergreen flag: (shold be 0):',int(EDPftvarcon_inst%evergreen(ipft))
              write(fates_log(),*) ' phen_stem_drop_fraction: ', EDPFtvarcon_inst%phen_stem_drop_fraction(ipft)
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

        ! Check if leaf storage priority is between 0-1
        ! ----------------------------------------------------------------------------------
        
        if ( ( EDPftvarcon_inst%leaf_stor_priority(ipft) < 0.0_r8 ) .or. &
             ( EDPftvarcon_inst%leaf_stor_priority(ipft) > 1.0_r8 ) ) then

           write(fates_log(),*) 'Prioritization of carbon allocation to leaf'
           write(fates_log(),*) ' and root turnover replacement, must be between'
           write(fates_log(),*) ' 0 and 1'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) 'leaf_stor_priority: ',EDPftvarcon_inst%leaf_stor_priority(ipft)
           write(fates_log(),*) ' Aborting'
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
        
        if ( ( EDPftvarcon_inst%c3psn(ipft) < 0.0_r8 ) .or. &
             ( EDPftvarcon_inst%c3psn(ipft) > 1.0_r8 ) ) then

           write(fates_log(),*) ' Two photosynthetic pathways are currently supported'
           write(fates_log(),*) ' C4 plants have c3psn = 0'
           write(fates_log(),*) ' C3 plants have c3psn = 1'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' c3psn(pft): ',EDPftvarcon_inst%c3psn(ipft)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))

        end if


        ! Check re-translocations
        ! Seems reasonable to assume that sapwood, structure and reproduction
        ! should not be re-translocating mass upon turnover.
        ! Note to advanced users. Feel free to remove these checks...
        ! -------------------------------------------------------------------
        
        if ( (EDPftvarcon_inst%turnover_carb_retrans(ipft,repro_organ) > nearzero) ) then
           write(fates_log(),*) ' Retranslocation of reproductive tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,repro_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,repro_organ) > nearzero) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,repro_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of reproductive tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,repro_organ)
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,repro_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,repro_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
           
        if ((EDPftvarcon_inst%turnover_carb_retrans(ipft,sapw_organ) > nearzero)) then
           write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,sapw_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,sapw_organ) > nearzero) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,sapw_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of sapwood tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,sapw_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        if ((EDPftvarcon_inst%turnover_carb_retrans(ipft,struct_organ) > nearzero)) then
           write(fates_log(),*) ' Retranslocation of structural(dead) tissues should be zero.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,struct_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,struct_organ) > nearzero) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,struct_organ) > nearzero) ) then
              write(fates_log(),*) ' Retranslocation of structural(dead) tissues should be zero.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,struct_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,struct_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Leaf retranslocation should be between 0 and 1
        if ( (EDPftvarcon_inst%turnover_carb_retrans(ipft,leaf_organ) > 1.0_r8) .or. & 
             (EDPftvarcon_inst%turnover_carb_retrans(ipft,leaf_organ) < 0.0_r8) ) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,leaf_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,leaf_organ) > 1.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,leaf_organ) > 1.0_r8) .or. &
               (EDPftvarcon_inst%turnover_nitr_retrans(ipft,leaf_organ) < 0.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,leaf_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,leaf_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,leaf_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Fineroot retranslocation should be between 0-1
        if ((EDPftvarcon_inst%turnover_carb_retrans(ipft,fnrt_organ) > 1.0_r8) .or. & 
            (EDPftvarcon_inst%turnover_carb_retrans(ipft,fnrt_organ) < 0.0_r8)) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,fnrt_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,fnrt_organ) > 1.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,fnrt_organ) > 1.0_r8) .or. &
               (EDPftvarcon_inst%turnover_nitr_retrans(ipft,fnrt_organ) < 0.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,fnrt_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,fnrt_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,fnrt_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Storage retranslocation should be between 0-1 (storage retrans seems weird, but who knows)
        if ((EDPftvarcon_inst%turnover_carb_retrans(ipft,store_organ) > 1.0_r8) .or. & 
            (EDPftvarcon_inst%turnover_carb_retrans(ipft,store_organ) < 0.0_r8)) then
           write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' carbon: ',EDPftvarcon_inst%turnover_carb_retrans(ipft,store_organ)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ((EDPftvarcon_inst%turnover_nitr_retrans(ipft,store_organ) > 1.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,store_organ) > 1.0_r8) .or. &
               (EDPftvarcon_inst%turnover_nitr_retrans(ipft,store_organ) < 0.0_r8) .or.  & 
               (EDPftvarcon_inst%turnover_phos_retrans(ipft,store_organ) < 0.0_r8)) then
              write(fates_log(),*) ' Retranslocation of leaf tissues should be between 0 and 1.'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' nitr: ',EDPftvarcon_inst%turnover_nitr_retrans(ipft,store_organ)
              write(fates_log(),*) ' phos: ',EDPftvarcon_inst%turnover_phos_retrans(ipft,store_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Growth respiration
        if (parteh_mode .eq. prt_carbon_allom_hyp) then
           if ( ( EDPftvarcon_inst%grperc(ipft) < 0.0_r8) .or. &
                ( EDPftvarcon_inst%grperc(ipft) > 1.0_r8 ) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',EDPftvarcon_inst%grperc(ipft)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        elseif(parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ( ( any(EDPftvarcon_inst%prt_grperc_organ(ipft,:) < 0.0_r8)) .or. &
                ( any(EDPftvarcon_inst%prt_grperc_organ(ipft,:) >= 1.0_r8)) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' Growth respiration must be between 0 and 1: ',EDPftvarcon_inst%prt_grperc_organ(ipft,:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Stoichiometric Ratios

        ! Firstly, the seed production and germination models cannot handle nutrients. So 
        ! we assume (for now) that seeds do not have nutrients (parteh_mode = 1 is c only)
        if(parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ( (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ) > nearzero) .or. &
                 (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ) < -nearzero) .or. & 
                 (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ) > nearzero) .or. &
                 (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ) < -nearzero) .or. &
                 (EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,repro_organ) > nearzero) .or. &
                 (EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,repro_organ) < -nearzero) .or. & 
                 (EDPftvarcon_inst%prt_phos_stoich_p2(ipft,repro_organ) > nearzero) .or. &
                 (EDPftvarcon_inst%prt_phos_stoich_p2(ipft,repro_organ) < -nearzero) ) then
              write(fates_log(),*) 'N & P should be zero in reproductive tissues'
              write(fates_log(),*) 'until nutrients are coupled into recruitment'
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ)
              write(fates_log(),*) EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ)
              write(fates_log(),*) EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,repro_organ)
              write(fates_log(),*) EDPftvarcon_inst%prt_phos_stoich_p2(ipft,repro_organ)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! The first nitrogen stoichiometry is used in all cases
        if ( (any(EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,:) < 0.0_r8)) .or. &
             (any(EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,:) >= 1.0_r8))) then
           write(fates_log(),*) ' PFT#: ',ipft
           write(fates_log(),*) ' N per C stoichiometry must bet between 0-1'
           write(fates_log(),*) EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,:)
           write(fates_log(),*) ' Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        if(parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if( (any(EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,:) < 0.0_r8)) .or. &
               (any(EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,:) >= 1.0_r8)) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' N per C stoichiometry must bet between 0-1'
              write(fates_log(),*) EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if

        ! Stoichiometric Ratios
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ( (any(EDPftvarcon_inst%prt_phos_stoich_p1(ipft,:) < 0.0_r8)) .or. &
                (any(EDPftvarcon_inst%prt_phos_stoich_p1(ipft,:) >= 1.0_r8)) .or. &
                (any(EDPftvarcon_inst%prt_phos_stoich_p2(ipft,:) < 0.0_r8)) .or. &
                (any(EDPftvarcon_inst%prt_phos_stoich_p2(ipft,:) >= 1.0_r8)) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' P per C stoichiometry must bet between 0-1'
              write(fates_log(),*) EDPftvarcon_inst%prt_phos_stoich_p1(ipft,:)
              write(fates_log(),*) EDPftvarcon_inst%prt_phos_stoich_p2(ipft,:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        if (parteh_mode .eq. prt_cnp_flex_allom_hyp) then
           if ( any(EDPftvarcon_inst%prt_alloc_priority(ipft,:) < 0) .or. &
                any(EDPftvarcon_inst%prt_alloc_priority(ipft,:) > 6) ) then
              write(fates_log(),*) ' PFT#: ',ipft
              write(fates_log(),*) ' Allocation priorities should be 0-6 for H1'
              write(fates_log(),*) EDPftvarcon_inst%prt_alloc_priority(ipft,:)
              write(fates_log(),*) ' Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if


        ! Check turnover time-scales
        
        nleafage = size(EDPftvarcon_inst%leaf_long,dim=2)

        do iage = 1, nleafage

           if ( EDPftvarcon_inst%leaf_long(ipft,iage)>nearzero ) then
              
              ! Check that leaf turnover doesn't exeed 1 day
              if ( (years_per_day / EDPftvarcon_inst%leaf_long(ipft,iage)) > 1._r8 ) then
                 write(fates_log(),*) 'Leaf turnover time-scale is greater than 1 day!'
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',EDPftvarcon_inst%leaf_long(ipft,iage),' [years]'
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
              
              ! Check to make sure that all other age-classes for this PFT also
              ! have non-zero entries, it wouldn't make sense otherwise
              if ( any(EDPftvarcon_inst%leaf_long(ipft,:) <= nearzero) ) then
                 write(fates_log(),*) 'You specified a leaf_long that is zero or'
                 write(fates_log(),*) 'invalid for a particular age class.'
                 write(fates_log(),*) 'Yet, other age classes for this PFT are non-zero.'
                 write(fates_log(),*) 'this doesnt make sense.'
                 write(fates_log(),*) 'ipft = ',ipft
                 write(fates_log(),*) 'leaf_long(ipft,:) =  ',EDPftvarcon_inst%leaf_long(ipft,:)
                 write(fates_log(),*) 'Aborting'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if

           else
              if (EDPftvarcon_inst%evergreen(ipft) .eq. itrue) then
                 write(fates_log(),*) 'You specified zero leaf turnover: '
                 write(fates_log(),*) 'ipft: ',ipft,' iage: ',iage
                 write(fates_log(),*) 'leaf_long(ipft,iage): ',EDPftvarcon_inst%leaf_long(ipft,iage)
                 write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
                 write(fates_log(),*) 'that an evergreen would have leaf maintenance turnover'
                 write(fates_log(),*) 'disable this error if you are ok with this'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end if

        end do

        ! Check the turnover rates on the senescing leaf pool
        if ( EDPftvarcon_inst%leaf_long(ipft,nleafage)>nearzero ) then
           
           ! Check that leaf turnover doesn't exeed 1 day
           if ( (years_per_day / &
                 (EDPftvarcon_inst%leaf_long(ipft,nleafage) * &
                  EDPftvarcon_inst%senleaf_long_fdrought(ipft))) > 1._r8 ) then
              write(fates_log(),*) 'Drought-senescent turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'leaf_long(ipft,nleafage)*senleaf_long_fdrought: ', &
                    EDPftvarcon_inst%leaf_long(ipft,nleafage)*EDPftvarcon_inst%senleaf_long_fdrought(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        if ( EDPftvarcon_inst%senleaf_long_fdrought(ipft)<nearzero .or. &
             EDPftvarcon_inst%senleaf_long_fdrought(ipft)>1._r8 ) then
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) must be greater than 0 '
           write(fates_log(),*) 'or less than or equal to 1.'
           write(fates_log(),*) 'Set this to 1 if you want no accelerated senescence turnover'
           write(fates_log(),*) 'ipft = ',ipft
           write(fates_log(),*) 'senleaf_long_fdrought(ipft) = ',EDPftvarcon_inst%senleaf_long_fdrought(ipft)
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
           

        if ( EDPftvarcon_inst%root_long(ipft)>nearzero ) then
           
           ! Check that root turnover doesn't exeed 1 day
           if ( (years_per_day / EDPftvarcon_inst%root_long(ipft)) > 1._r8 ) then
              write(fates_log(),*) 'Root turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',EDPftvarcon_inst%root_long(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           
        else
           if (EDPftvarcon_inst%evergreen(ipft) .eq. itrue) then
              write(fates_log(),*) 'You specified zero root turnover: '
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'root_long(ipft): ',EDPftvarcon_inst%root_long(ipft)
              write(fates_log(),*) 'yet this is an evergreen PFT, and it only makes sense'
              write(fates_log(),*) 'that an evergreen would have root maintenance turnover'
              write(fates_log(),*) 'disable this error if you are ok with this'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if
        
        ! Check Branch turnover doesn't exceed one day
        if ( EDPftvarcon_inst%branch_turnover(ipft)>nearzero ) then
           
           ! Check that branch turnover doesn't exeed 1 day
           if ( (years_per_day / EDPftvarcon_inst%branch_turnover(ipft)) > 1._r8 ) then
              write(fates_log(),*) 'Branch turnover time-scale is greater than 1 day!'
              write(fates_log(),*) 'ipft: ',ipft
              write(fates_log(),*) 'branch_turnover(ipft): ',EDPftvarcon_inst%branch_turnover(ipft),' [years]'
              write(fates_log(),*) 'Aborting'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end if


     end do


     return
  end subroutine FatesCheckParams




end module EDPftvarcon

