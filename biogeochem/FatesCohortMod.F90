module FatesCohortMod

  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesConstantsMod,          only : fates_unset_int
  use FatesConstantsMod,          only : ifalse, itrue
  use FatesConstantsMod,          only : nearzero
  use FatesConstantsMod,          only : ican_upper, ican_ustory
  use EDParamsMod,                only : nlevleaf
  use FatesGlobals,               only : endrun => fates_endrun
  use FatesGlobals,               only : fates_log
  use PRTGenericMod,              only : max_nleafage
  use PRTGenericMod,              only : prt_vartypes
  use PRTGenericMod,              only : prt_carbon_allom_hyp
  use PRTGenericMod,              only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,              only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,              only : repro_organ, store_organ, struct_organ
  use PRTGenericMod,              only : carbon12_element
  use PRTParametersMod,           only : prt_params
  use FatesParameterDerivedMod,   only : param_derived
  use FatesHydraulicsMemMod,      only : ed_cohort_hydr_type
  use FatesInterfaceTypesMod,     only : hlm_parteh_mode
  use FatesInterfaceTypesMod,     only : hlm_use_sp
  use FatesInterfaceTypesMod,     only : hlm_use_planthydro
  use FatesInterfaceTypesMod,     only : nleafage
  use EDPftvarcon,                only : EDPftvarcon_inst
  use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
  use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
  use FatesAllometryMod,          only : carea_allom, tree_lai_sai
  use PRTAllometricCarbonMod,     only : ac_bc_inout_id_dbh, ac_bc_inout_id_netdc
  use PRTAllometricCarbonMod,     only : ac_bc_in_id_cdamage, ac_bc_in_id_pft
  use PRTAllometricCarbonMod,     only : ac_bc_in_id_ctrim, ac_bc_in_id_lstat
  use PRTAllometricCarbonMod,     only : ac_bc_in_id_efleaf
  use PRTAllometricCarbonMod,     only : ac_bc_in_id_effnrt
  use PRTAllometricCarbonMod,     only : ac_bc_in_id_efstem
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_pft, acnp_bc_in_id_ctrim
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_lstat, acnp_bc_in_id_netdc
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_netdc, acnp_bc_in_id_nc_repro
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_pc_repro, acnp_bc_in_id_cdamage
  use PRTAllometricCNPMod,        only : acnp_bc_inout_id_dbh, acnp_bc_inout_id_resp_excess
  use PRTAllometricCNPMod,        only : acnp_bc_inout_id_l2fr, acnp_bc_inout_id_cx_int
  use PRTAllometricCNPMod,        only : acnp_bc_inout_id_emadcxdt, acnp_bc_inout_id_cx0
  use PRTAllometricCNPMod,        only : acnp_bc_inout_id_netdn, acnp_bc_inout_id_netdp
  use PRTAllometricCNPMod,        only : acnp_bc_out_id_cefflux, acnp_bc_out_id_nefflux
  use PRTAllometricCNPMod,        only : acnp_bc_out_id_pefflux, acnp_bc_out_id_limiter
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_efleaf
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_effnrt
  use PRTAllometricCNPMod,        only : acnp_bc_in_id_efstem

  use shr_infnan_mod,             only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod,                only : errMsg => shr_log_errMsg

  implicit none
  private
  
  ! PARAMETERS
  character(len=*), parameter, private :: sourcefile = __FILE__

  ! FATES COHORT TYPE
  type, public :: fates_cohort_type

    ! POINTERS
    type (fates_cohort_type), pointer :: taller   => null() ! pointer to next tallest cohort     
    type (fates_cohort_type), pointer :: shorter  => null() ! pointer to next shorter cohort
    
    !---------------------------------------------------------------------------

    ! Multi-species, multi-organ Plant Reactive Transport (PRT)
    ! Contains carbon and nutrient state variables for various plant organs
    class(prt_vartypes), pointer :: prt
    real(r8)                     :: l2fr ! leaf to fineroot biomass ratio [kg root / kg leaf]
                                            ! (this is constant in carbon only simulationss, and 
                                            ! is set by the allom_l2fr parameter).  
                                            ! For nutrient enabled simulations, this is dynamic.  
                                            ! In cold-start simulations, the allom_l2fr 
                                            ! parameter sets the starter value. 

    !---------------------------------------------------------------------------

    ! VEGETATION STRUCTURE
    integer  :: pft                     ! pft index
    real(r8) :: n                       ! number of individuals in cohort per 'area' (10000m2 default) [/m2]
    real(r8) :: dbh                     ! diameter at breast height [cm]
    real(r8) :: coage                   ! age [years]
    real(r8) :: height                  ! height [m]
    integer  :: indexnumber             ! unique number for each cohort (within clump?)
    integer  :: canopy_layer            ! canopy status of cohort [1 = canopy, 2 = understorey, etc.]
    real(r8) :: canopy_layer_yesterday  ! recent canopy status of cohort [1 = canopy, 2 = understorey, etc.]
                                        !   real to be conservative during fusion
    integer  :: crowndamage             ! crown damage class of the cohort [1 = undamaged, >1 = damaged]                     
    real(r8) :: g_sb_laweight           ! total conductance (stomata + boundary layer) of the cohort
                                        !   weighted by its leaf area [m/s]*[m2]
    real(r8) :: canopy_trim             ! fraction of the maximum leaf biomass that we are targeting [0-1]
    real(r8) :: leaf_cost               ! how much does it cost to maintain leaves [kgC/m2/year]
    real(r8) :: excl_weight             ! how much of this cohort is demoted each year, as a proportion of all cohorts
    real(r8) :: prom_weight             ! how much of this cohort is promoted each year, as a proportion of all cohorts
    integer  :: nv                      ! number of leaf layers
    integer  :: status_coh              ! growth status of plant  [2 = leaves on , 1 = leaves off]
    real(r8) :: efleaf_coh              ! elongation factor for leaves [fraction]
    real(r8) ::  effnrt_coh             ! elongation factor for fine roots [fraction]
    real(r8) ::  efstem_coh             ! elongation factor for stem [fraction]
                                        !   for all the elongation factors, 0 means fully abscissed, and 
                                        !   1 means fully flushed.
    real(r8) :: c_area                  ! areal extent of canopy [m2]
    real(r8) :: treelai                 ! lai of an individual within cohort leaf area [m2 leaf area/m2 crown area]
    real(r8) :: treesai                 ! stem area index of an individual within cohort [m2 stem area/m2 crown area]
    logical  :: isnew                   ! flag to signify a new cohort - new cohorts have not experienced
                                        !   npp or mortality and should therefore not be fused or averaged
    integer  :: size_class              ! index that indicates which diameter size bin the cohort currently resides in
                                        !   this is used for history output. We maintain this in the main cohort memory
                                        !   because we don't want to continually re-calculate the cohort's position when
                                        !   performing size diagnostics at high-frequency calls
    integer  :: coage_class             ! index that indicates which age bin the cohort currently resides in 
                                        !   (used for history output)
    integer  :: size_by_pft_class       ! index that indicates the cohorts position of the joint size-class x functional
                                        !   type classification. We also maintain this in the main cohort memory
                                        !   because we don't want to continually re-calculate the cohort's position when
                                        !   performing size diagnostics at high-frequency calls
    integer  :: coage_by_pft_class      ! index that indicates the cohorts position of the join cohort age class x PFT 
    integer  :: size_class_lasttimestep ! size class of the cohort at the last time step

    !---------------------------------------------------------------------------

    ! CARBON AND NUTRIENT FLUXES 

    ! --------------------------------------------------------------------------
    ! NPP, GPP and RESP: Instantaneous, accumulated and accumulated-hold types*
    ! 
    ! _tstep:    The instantaneous estimate that is calculated at each rapid plant biophysics
    !            time-step (ie photosynthesis, sub-hourly) [kgC/indiv/timestep]
    ! _acc:      The accumulation of the _tstep variable from the beginning to ending of
    !            the dynamics time-scale.  This variable is zero'd during initialization and
    !            after the dynamics call-sequence is completed.  [kgC/indiv/day]
    ! _acc_hold: While _acc is zero'd after the dynamics call sequence and then integrated, 
    !            _acc_hold "holds" the integrated value until the next time dynamics is 
    !            called. This is useful because growth and excess respiration
    !            are calculated once daily, but we want to remove the average
    !            flux from the daily NEP signal, so we remove it from the next day.
    !            The hold variables are also useful for rebuilding history on restart.
    !            Units converted to a useful rate [kgC/indiv/yr]
    ! --------------------------------------------------------------------------

    real(r8) :: gpp_tstep                 ! Gross Primary Production (see above *)
    real(r8) :: gpp_acc
    real(r8) :: gpp_acc_hold

    real(r8) :: npp_acc
    real(r8) :: npp_acc_hold

    real(r8) :: resp_m_tstep              ! Maintenance respiration (see above *)
    real(r8) :: resp_m_acc
    real(r8) :: resp_m_acc_hold
    real(r8) :: resp_g_acc_hold

    real(r8) :: c13disc_clm               ! carbon 13 discrimination in new synthesized carbon at each indiv/timestep [ppm]
    real(r8) :: c13disc_acc               ! carbon 13 discrimination in new synthesized carbon at each indiv/day
                                          !   at the end of a day [ppm]

    ! The following four biophysical rates are assumed to be at the canopy top, at reference temp 25degC, 
    ! and based on the leaf age weighted average of the PFT parameterized values. 
    ! The last condition is why it is dynamic and tied to the cohort

    real(r8) :: vcmax25top                ! maximum carboxylation at canopy top and 25degC [umol CO2/m2/s]
    real(r8) :: jmax25top                 ! maximum electron transport rate at canopy top and 25degC [umol electrons/m2/s]
    real(r8) :: tpu25top                  ! triose phosphate utilization rate at canopy top and 25degC [umol CO2/m2/s]
    real(r8) :: kp25top                   ! initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: ts_net_uptake(nlevleaf)   ! net uptake of leaf layers [kgC/m2/timestep]
    real(r8) :: year_net_uptake(nlevleaf) ! net uptake of leaf layers [kgC/m2/year]

    ! used for CNP
    integer  :: cnp_limiter               ! which element is limiting growth [0 = none, 1 = C, 2 = N, 3 = P]
    real(r8) :: cx_int                    ! time integration of the log of the relative carbon storage over relative nutrient
    real(r8) :: ema_dcxdt                 ! derivative of the log of the relative carbon storage over relative nutrient
    real(r8) :: cx0                       !   value on the previous time-step of log of the relative carbon storage over 
                                          !   relative nutrient
    real(r8) :: nc_repro                  ! N:C ratio of a new recruit, used also for defining reproductive stoich
    real(r8) :: pc_repro                  ! P:C ratio of a new recruit

    ! Nutrient Fluxes (if N, P, etc. are turned on)
    real(r8) :: daily_nh4_uptake          ! integrated daily uptake of mineralized ammonium through competitive acquisition 
                                          !   in soil [kgN/plant/day]
    real(r8) :: daily_no3_uptake          ! integrated daily uptake of mineralized nitrate through competitive acquisition 
                                          !   in soil [kgN/plant/day]

    real(r8) :: sym_nfix_daily            ! accumulated symbiotic N fixation from the roots [kgN/indiv/day]
    real(r8) :: sym_nfix_tstep            ! symbiotic N fixation from the roots for the time-step [kgN/indiv/timestep]

    real(r8) :: daily_n_gain              ! sum of fixation and uptake of mineralized NH4/NO3 in solution as well as 
                                          !   symbiotic fixation
    real(r8) :: daily_p_gain              ! integrated daily uptake of mineralized P through competitive acquisition 
                                          !   in soil [kgP/plant/day]

    real(r8) :: daily_c_efflux            ! daily mean efflux of excess carbon from roots into labile pool [kgC/plant/day]
    real(r8) :: daily_n_efflux            ! daily mean efflux of excess nitrogen from roots into labile pool [kgN/plant/day]
    real(r8) :: daily_p_efflux            ! daily mean efflux of excess phophorus from roots into labile pool [kgP/plant/day]

    real(r8) :: daily_n_demand            ! daily amount of N demanded by the plant [kgN/plant/day]
    real(r8) :: daily_p_demand            ! daily amount of P demanded by the plant [kgN/plant/day]

    real(r8) :: seed_prod                 ! diagnostic seed production rate [kgC/plant/day]

    !---------------------------------------------------------------------------

    integer :: twostr_col  ! The column index in the two-stream solution that this cohort is part of
    
    ! RESPIRATION COMPONENTS
    real(r8) :: resp_excess_hold ! respiration of excess carbon [kgC/indiv/yr]
                                 ! note: this is flagged "hold" because it is calculated
                                 ! at the end of the day (dynamics) but is used
                                 ! on the following day (like growth respiration)
                                 ! to aid in reporting a more accurate sub-daily
                                 ! NEP

    real(r8) :: rdark            ! dark respiration [kgC/indiv/s]
    real(r8) :: resp_m_unreduced ! diagnostic-only unreduced maintenance respiration [kgC/indiv/timestep]
    real(r8) :: livestem_mr      ! aboveground live stem maintenance respiration [kgC/indiv/s]
    real(r8) :: livecroot_mr     ! belowground live stem maintenance respiration [kgC/indiv/s]
    real(r8) :: froot_mr         ! live fine root maintenance respiration [kgC/indiv/s]

    !---------------------------------------------------------------------------

    ! DAMAGE
    real(r8) :: branch_frac ! fraction of aboveground woody biomass in branches [0-1]

    !---------------------------------------------------------------------------

    ! MORTALITY
    real(r8) :: dmort            ! proportional mortality rate [/year]

    ! Mortality Rate Partitions
    real(r8) :: bmort            ! background mortality rate [indiv/year]
    real(r8) :: cmort            ! carbon starvation mortality rate [indiv/year]
    real(r8) :: hmort            ! hydraulic failure mortality rate [indiv/year]
    real(r8) :: frmort           ! freezing mortality rate [indiv/year]
    real(r8) :: smort            ! senesence mortality [indiv/year]
    real(r8) :: asmort           ! age senescence mortality [indiv/year]
    real(r8) :: dgmort           ! damage mortality [indiv/year]

    ! Logging Mortality Rate 
    ! Yi Xu & M. Huang
    real(r8) :: lmort_direct     ! directly logging rate [fraction/logging activity]
    real(r8) :: lmort_collateral ! collaterally damaged rate [fraction/logging activity]
    real(r8) :: lmort_infra      ! mechanically damaged rate [fraction/logging activity]
    real(r8) :: l_degrad         ! rate of trees that are not killed but suffer from forest degradation
                                 !  (i.e. they are moved to newly-anthro-disturbed secondary 
                                 !  forest patch)  [fraction/logging activity]

    !---------------------------------------------------------------------------

    ! NITROGEN POOLS      
    ! --------------------------------------------------------------------------
    ! Nitrogen pools are not prognostic in the current implementation.
    ! They are diagnosed during photosynthesis using a simple C2N parameter. 
    ! Local values are used in that routine.
    ! --------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    ! GROWTH DERIVIATIVES
    real(r8) :: dndt      ! time derivative of cohort size [n/year]
    real(r8) :: dhdt      ! time derivative of height [m/year]
    real(r8) :: ddbhdt    ! time derivative of dbh [cm/year]
    real(r8) :: dbdeaddt  ! time derivative of dead biomass [kgC/year]

    !---------------------------------------------------------------------------

    ! FIRE
    real(r8) ::  fraction_crown_burned ! proportion of crown affected by fire [0-1]
    real(r8) ::  cambial_mort          ! probability that trees dies due to cambial charring [0-1]
                                       !  (conditional on the tree being subjected to the fire)
    real(r8) ::  crownfire_mort        ! probability of tree post-fire mortality from crown scorch [0-1]
                                       !  (conditional on the tree being subjected to the fire)
    real(r8) ::  fire_mort             ! post-fire mortality from cambial and crown damage assuming two are independent [0-1]

    !---------------------------------------------------------------------------

    ! HYDRAULICS
    type(ed_cohort_hydr_type), pointer :: co_hydr ! all cohort hydraulics data, see FatesHydraulicsMemMod.F90

    contains
    
    procedure :: Init
    procedure :: NanValues
    procedure :: ZeroValues
    procedure :: Create
    procedure :: Copy
    procedure :: FreeMemory
    procedure :: CanUpperUnder
    procedure :: InitPRTBoundaryConditions
    procedure :: UpdateCohortBioPhysRates
    procedure :: Dump

  end type fates_cohort_type

  contains

    !===========================================================================

    subroutine Init(this, prt)
      !
      !  DESCRIPTION:
      !  Create new cohort and set default values for all variables
      !
  
      ! ARGUMENTS:
      class(fates_cohort_type), intent(inout)          :: this
      class(prt_vartypes),      intent(inout), pointer :: prt ! allocated PARTEH object
  
      call this%NanValues()  ! make everything in the cohort not-a-number
      call this%ZeroValues() ! zero things that need to be zeroed
      
      ! point to the PARTEH object
      this%prt => prt
  
      ! The PARTEH cohort object should be allocated and already
      ! initialized in this routine.
      call this%prt%CheckInitialConditions()
    
      ! new cohorts do not have mortality rates, nor have they moved any
      ! carbon when they are created.  They will bias our statistics
      ! until they have experienced a full day.  We need a newly recruited flag.
      ! This flag will be set to false after it has experienced
      ! growth, disturbance and mortality.
      this%isnew = .true.
  
    end subroutine Init
  
    !===========================================================================

    subroutine NanValues(this)
      !
      ! DESCRIPTION:
      !  make all the cohort variables NaN or unset so they aren't used before defined
      !
   
      ! ARGUMENTS:
      class(fates_cohort_type), intent(inout) :: this
    
      ! set pointers to null
      this%taller => null()  
      this%shorter => null()
      this%prt => null()
      this%co_hydr => null()
      nullify(this%taller)
      nullify(this%shorter)
      nullify(this%prt)
      nullify(this%co_hydr)
   
      ! VEGETATION STRUCTURE
      this%l2fr                    = nan 
      this%pft                     = fates_unset_int  
      this%n                       = nan
      this%dbh                     = nan
      this%coage                   = nan 
      this%height                  = nan 
      this%indexnumber             = fates_unset_int
      this%canopy_layer            = fates_unset_int
      this%canopy_layer_yesterday  = nan  
      this%crowndamage             = fates_unset_int
      this%g_sb_laweight           = nan
      this%canopy_trim             = nan
      this%leaf_cost               = nan
      this%excl_weight             = nan
      this%prom_weight             = nan 
      this%nv                      = fates_unset_int  
      this%status_coh              = fates_unset_int
      this%efleaf_coh              = nan
      this%effnrt_coh              = nan 
      this%efstem_coh              = nan
      this%c_area                  = nan 
      this%treelai                 = nan
      this%treesai                 = nan
      this%isnew                   = .false.
      this%size_class              = fates_unset_int
      this%coage_class             = fates_unset_int
      this%size_by_pft_class       = fates_unset_int
      this%coage_by_pft_class      = fates_unset_int
      this%size_class_lasttimestep = fates_unset_int
   
      ! CARBON AND NUTRIENT FLUXES 
      this%gpp_tstep               = nan
      this%gpp_acc                 = nan
      this%gpp_acc_hold            = nan
      this%npp_acc                 = nan 
      this%npp_acc_hold            = nan
      this%resp_m_tstep            = nan 
      this%resp_m_acc              = nan 
      this%resp_m_acc_hold         = nan
      this%resp_g_acc_hold         = nan
      this%c13disc_clm             = nan
      this%c13disc_acc             = nan
      this%vcmax25top              = nan
      this%jmax25top               = nan
      this%tpu25top                = nan
      this%kp25top                 = nan
      this%year_net_uptake(:)      = nan 
      this%ts_net_uptake(:)        = nan
      this%cnp_limiter             = fates_unset_int
      this%cx_int                  = nan
      this%ema_dcxdt               = nan
      this%cx0                     = nan
      this%nc_repro                = nan
      this%pc_repro                = nan
      this%daily_nh4_uptake        = nan
      this%daily_no3_uptake        = nan
      this%sym_nfix_daily          = nan
      this%sym_nfix_tstep          = nan
      this%daily_n_gain            = nan
      this%daily_p_gain            = nan
      this%daily_c_efflux          = nan
      this%daily_n_efflux          = nan
      this%daily_p_efflux          = nan
      this%daily_n_demand          = nan
      this%daily_p_demand          = nan
      this%seed_prod               = nan
   
      ! RESPIRATION COMPONENTS
      this%rdark                   = nan
      this%resp_m_unreduced        = nan 
      this%resp_excess_hold        = nan 
      this%livestem_mr             = nan 
      this%livecroot_mr            = nan 
      this%froot_mr                = nan 
   
      ! DAMAGE
      this%branch_frac             = nan 
   
      ! MORTALITY
      this%dmort                   = nan
      this%bmort                   = nan 
      this%cmort                   = nan
      this%frmort                  = nan 
      this%smort                   = nan 
      this%asmort                  = nan 
      this%dgmort                  = nan 
      this%lmort_direct            = nan
      this%lmort_collateral        = nan
      this%lmort_infra             = nan
      this%l_degrad                = nan
   
      ! GROWTH DERIVATIVES
      this%dndt                    = nan 
      this%dhdt                    = nan 
      this%ddbhdt                  = nan
      this%dbdeaddt                = nan 
   
      ! FIRE
      this%fraction_crown_burned   = nan 
      this%cambial_mort            = nan 
      this%crownfire_mort          = nan 
      this%fire_mort               = nan 
      
    end subroutine NanValues
   
    !===========================================================================
   
    subroutine ZeroValues(this)
      !
      ! DESCRIPTION:
      ! Zero variables that need to be accounted for if this cohort is altered 
      ! before they are defined.
      !
      ! ARGUMENTS
      class(fates_cohort_type), intent(inout) :: this
      
      this%g_sb_laweight           = 0._r8
   
      this%leaf_cost               = 0._r8
      this%excl_weight             = 0._r8
      this%prom_weight             = 0._r8
      this%nv                      = 0
      this%status_coh              = 0
      this%efleaf_coh              = 0.0_r8
      this%effnrt_coh              = 0.0_r8 
      this%efstem_coh              = 0.0_r8
    
      this%treesai                 = 0._r8
      this%size_class              = 1
      this%coage_class             = 1
   
      this%size_class_lasttimestep = 0
      this%gpp_tstep               = 0._r8
      this%gpp_acc                 = 0._r8
      this%npp_acc                 = 0._r8
      this%resp_m_tstep            = 0._r8
      this%resp_m_acc              = 0._r8

      ! do not zero these, they are not built
      ! so more appropriate to leave unzerod
      ! to prevent uninitialized use
      ! this%gpp_acc_hold            = nan
      ! this%npp_acc_hold            = nan
      ! this%resp_m_acc_hold         = nan
      ! this%resp_g_acc_hold         = nan
      ! this%resp_excess_hold        = nan
      
      this%c13disc_clm             = 0._r8
      this%c13disc_acc             = 0._r8
   
      this%ts_net_uptake(:)        = 0._r8
      this%year_net_uptake(:)      = 999._r8 ! this needs to be 999, or trimming of new cohorts will break.
   
      this%daily_nh4_uptake        = 0._r8
      this%daily_no3_uptake        = 0._r8
      
      ! fixation is also integrated over the course of the day and must be 
      !     zeroed upon creation and after plant resource allocation
      this%sym_nfix_daily          = 0._r8
      this%daily_n_gain            = 0._r8
      this%daily_p_gain            = 0._r8
   
      ! daily nutrient fluxes are INTEGRATED over the course of the day.  
      !     These variables MUST be zerod upon creation AND after allocation. 
      !     These variables exist in carbon-only mode but are not used.
      this%daily_c_efflux          = 0._r8
      this%daily_n_efflux          = 0._r8
      this%daily_p_efflux          = 0._r8
   
      ! initialize these as negative
      this%daily_n_demand          = -9._r8
      this%daily_p_demand          = -9._r8
      this%seed_prod               = 0._r8
      this%rdark                   = 0._r8
      this%resp_m_unreduced        = 0._r8
      this%livestem_mr             = 0._r8
      this%livecroot_mr            = 0._r8
      this%froot_mr                = 0._r8
   
      this%dmort                   = 0._r8
      this%lmort_direct            = 0._r8
      this%lmort_collateral        = 0._r8
      this%lmort_infra             = 0._r8
      this%l_degrad                = 0._r8
      this%fraction_crown_burned   = 0._r8
      this%cambial_mort            = 0._r8
      this%crownfire_mort          = 0._r8
      this%fire_mort               = 0._r8
    
    end subroutine ZeroValues
   
    !===========================================================================

    subroutine Create(this, prt, pft, nn, height, coage, dbh, status,            &
      ctrim, carea, clayer, crowndamage, spread, can_tlai, elongf_leaf,        &
      elongf_fnrt, elongf_stem)
      !
      ! DESCRIPTION:
      ! set up values for a newly created cohort
      
      ! ARGUMENTS
      class(fates_cohort_type), intent(inout), target  :: this             ! cohort object
      class(prt_vartypes),      intent(inout), pointer :: prt              ! The allocated PARTEH object
      integer,                  intent(in)             :: pft              ! cohort Plant Functional Type
      integer,                  intent(in)             :: crowndamage      ! cohort damage class 
      integer,                  intent(in)             :: clayer           ! canopy status of cohort [canopy/understory]
      integer,                  intent(in)             :: status           ! growth status of cohort [leaves on/off]
      real(r8),                 intent(in)             :: nn               ! number of individuals in cohort [/m2]
      real(r8),                 intent(in)             :: height           ! cohort height [m]
      real(r8),                 intent(in)             :: coage            ! cohort age [yr]
      real(r8),                 intent(in)             :: dbh              ! cohort diameter at breat height [cm]
      real(r8),                 intent(in)             :: ctrim            ! fraction of the maximum leaf biomass 
      real(r8),                 intent(in)             :: spread           ! how spread crowns are in horizontal space
      real(r8),                 intent(in)             :: carea            ! area of cohort, for SP mode [m2]
      real(r8),                 intent(in)             :: can_tlai(:)      ! patch-level total LAI of each canopy layer
      real(r8),                 intent(in)             :: elongf_leaf      ! leaf elongation factor [fraction]
      real(r8),                 intent(in)             :: elongf_fnrt      ! fine-root "elongation factor" [fraction]
      real(r8),                 intent(in)             :: elongf_stem      ! stem "elongation factor" [fraction]

      ! LOCAL VARIABLES:
      integer  :: iage        ! loop counter for leaf age classes
      real(r8) :: leaf_c      ! total leaf carbon [kgC]
      real(r8) :: treesai     ! stem area index within crown [m2/m2]
      
      ! initialize cohort
      call this%Init(prt)
      
      ! set values
      this%pft          = pft
      this%crowndamage  = crowndamage
      this%canopy_layer = clayer
      this%canopy_layer_yesterday = real(clayer, r8)
      this%status_coh   = status
      this%n            = nn
      this%height       = height
      this%dbh          = dbh
      this%coage        = coage
      this%canopy_trim  = ctrim
      this%efleaf_coh   = elongf_leaf
      this%effnrt_coh   = elongf_fnrt
      this%efstem_coh   = elongf_stem

      ! This routine may be called during restarts, and at this point in the call sequence
      ! the actual cohort data is unknown, as this is really only used for allocation
      ! In these cases, testing if things like biomass are reasonable is premature
      ! However, in this part of the code, we will pass in nominal values for size, number and type
      if (this%dbh <= 0._r8 .or. this%n == 0._r8 .or. this%pft == 0) then
        write(fates_log(),*) 'FATES: something is zero in cohort%Create',      &
          this%dbh, this%n, this%pft
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif

      ! Initialize the leaf to fineroot biomass ratio.
      ! For C-only, this will stay constant, for nutrient-enabled this will be
      ! dynamic.  In both cases, new cohorts are initialized with the minimum. 
      ! This works in the nutrient enabled case because cohorts are also 
      ! initialized with full stores, which match with minimum fineroot biomass
      this%l2fr = prt_params%allom_l2fr(pft)

      if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
        this%cx_int      = 0._r8  ! Assume balanced N,P/C stores ie log(1) = 0
        this%cx0         = 0._r8  ! Assume balanced N,P/C stores ie log(1) = 0
        this%ema_dcxdt   = 0._r8  ! Assume unchanged dCX/dt
        this%cnp_limiter = 0      ! Assume limitations are unknown
      end if

      ! This sets things like vcmax25top, that depend on the leaf age fractions 
      ! (which are defined by PARTEH)
      call this%UpdateCohortBioPhysRates()

      ! calculate size classes
      call sizetype_class_index(this%dbh, this%pft, this%size_class,           &
        this%size_by_pft_class)

      ! If cohort age tracking is off we call this here once, just so everything
      ! is in the first bin. This makes it easier to copy and terminate cohorts 
      ! later.
      ! We don't need to update this ever if cohort age tracking is off
      call coagetype_class_index(this%coage, this%pft, this%coage_class,       &
        this%coage_by_pft_class)

      ! asssign or calculate canopy extent and depth
      if (hlm_use_sp .eq. ifalse) then
        call carea_allom(this%dbh, this%n, spread, this%pft, this%crowndamage, &
          this%c_area)
      else
        ! set this from previously precision-controlled value in SP mode
        this%c_area = carea 
      endif

      ! Query PARTEH for the leaf carbon [kg]
      leaf_c = this%prt%GetState(leaf_organ, carbon12_element)

      call tree_lai_sai(leaf_c, this%pft, this%c_area, this%n,           &
           this%canopy_layer, can_tlai, this%vcmax25top, this%dbh, this%crowndamage,          &
           this%canopy_trim, this%efstem_coh, 2, this%treelai, treesai)

      if (hlm_use_sp .eq. ifalse) then
         this%treesai = treesai
      end if
     

      call this%InitPRTBoundaryConditions()

    end subroutine Create

    !===========================================================================

    subroutine Copy(this, copyCohort) 
      !
      ! DESCRIPTION:
      ! copies all the variables in one cohort into a new cohort
      !

      ! ARGUMENTS
      class(fates_cohort_type), intent(in)    :: this       ! old cohort 
      class(fates_cohort_type), intent(inout) :: copyCohort ! new cohort

      copyCohort%indexnumber = fates_unset_int
      
      ! POINTERS
      copyCohort%taller  => NULL() 
      copyCohort%shorter => NULL() 

      ! PRT
      call copyCohort%prt%CopyPRTVartypes(this%prt)
      copyCohort%l2fr                    = this%l2fr
      
      ! VEGETATION STRUCTURE
      copyCohort%pft                     = this%pft
      copyCohort%n                       = this%n
      copyCohort%dbh                     = this%dbh
      copyCohort%coage                   = this%coage
      copyCohort%height                  = this%height
      copyCohort%canopy_layer            = this%canopy_layer
      copyCohort%canopy_layer_yesterday  = this%canopy_layer_yesterday
      copyCohort%crowndamage             = this%crowndamage
      copyCohort%g_sb_laweight           = this%g_sb_laweight
      copyCohort%canopy_trim             = this%canopy_trim
      copyCohort%leaf_cost               = this%leaf_cost
      copyCohort%excl_weight             = this%excl_weight
      copyCohort%prom_weight             = this%prom_weight
      copyCohort%nv                      = this%nv
      copyCohort%status_coh              = this%status_coh
      copyCohort%efleaf_coh              = this%efleaf_coh 
      copyCohort%effnrt_coh              = this%effnrt_coh 
      copyCohort%efstem_coh              = this%efstem_coh
      copyCohort%c_area                  = this%c_area
      copyCohort%treelai                 = this%treelai
      copyCohort%treesai                 = this%treesai
      copyCohort%isnew                   = this%isnew
      copyCohort%size_class              = this%size_class
      copyCohort%coage_class             = this%coage_class
      copyCohort%size_by_pft_class       = this%size_by_pft_class
      copyCohort%coage_by_pft_class      = this%coage_by_pft_class
      copyCohort%size_class_lasttimestep = this%size_class_lasttimestep

      ! CARBON AND NUTRIENT FLUXES
      copyCohort%gpp_tstep               = this%gpp_tstep
      copyCohort%gpp_acc                 = this%gpp_acc
      copyCohort%gpp_acc_hold            = this%gpp_acc_hold
      copyCohort%npp_acc                 = this%npp_acc
      copyCohort%npp_acc_hold            = this%npp_acc_hold
      copyCohort%resp_m_tstep            = this%resp_m_tstep
      copyCohort%resp_m_acc              = this%resp_m_acc
      copyCohort%resp_m_acc_hold         = this%resp_m_acc_hold
      copyCohort%resp_g_acc_hold         = this%resp_g_acc_hold
      copyCohort%c13disc_clm             = this%c13disc_clm
      copyCohort%c13disc_acc             = this%c13disc_acc
      copyCohort%vcmax25top              = this%vcmax25top
      copyCohort%jmax25top               = this%jmax25top
      copyCohort%tpu25top                = this%tpu25top
      copyCohort%kp25top                 = this%kp25top
      copyCohort%ts_net_uptake           = this%ts_net_uptake
      copyCohort%year_net_uptake         = this%year_net_uptake
      copyCohort%cnp_limiter             = this%cnp_limiter

      if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then 
        copyCohort%cx_int                  = this%cx_int
        copyCohort%ema_dcxdt               = this%ema_dcxdt
        copyCohort%cx0                     = this%cx0
      end if 

      copyCohort%nc_repro                = this%nc_repro
      copyCohort%daily_nh4_uptake        = this%daily_nh4_uptake
      copyCohort%daily_no3_uptake        = this%daily_no3_uptake
      copyCohort%sym_nfix_daily          = this%sym_nfix_daily
      copyCohort%sym_nfix_tstep          = this%sym_nfix_tstep
      copyCohort%daily_n_gain            = this%daily_n_gain
      copyCohort%daily_p_gain            = this%daily_p_gain
      copyCohort%daily_c_efflux          = this%daily_c_efflux
      copyCohort%daily_n_efflux          = this%daily_n_efflux
      copyCohort%daily_p_efflux          = this%daily_p_efflux
      copyCohort%daily_n_demand          = this%daily_n_demand
      copyCohort%daily_p_demand          = this%daily_p_demand
      copyCohort%seed_prod               = this%seed_prod

      ! RESPIRATION COMPONENTS
      copyCohort%rdark                   = this%rdark
      copyCohort%resp_m_unreduced        = this%resp_m_unreduced
      copyCohort%resp_excess_hold        = this%resp_excess_hold
      copyCohort%livestem_mr             = this%livestem_mr
      copyCohort%livecroot_mr            = this%livecroot_mr
      copyCohort%froot_mr                = this%froot_mr

      ! DAMAGE
      copyCohort%branch_frac             = this%branch_frac

      ! MORTALITY
      copyCohort%dmort                   = this%dmort
      copyCohort%bmort                   = this%bmort
      copyCohort%cmort                   = this%cmort
      copyCohort%hmort                   = this%hmort
      copyCohort%frmort                  = this%frmort
      copyCohort%smort                   = this%smort
      copyCohort%asmort                  = this%asmort
      copyCohort%dgmort                  = this%dgmort
      copyCohort%lmort_direct            = this%lmort_direct
      copyCohort%lmort_collateral        = this%lmort_collateral
      copyCohort%lmort_infra             = this%lmort_infra
      copyCohort%l_degrad                = this%l_degrad

      ! GROWTH DERIVATIVES
      copyCohort%dndt                    = this%dndt
      copyCohort%dhdt                    = this%dhdt
      copyCohort%ddbhdt                  = this%ddbhdt
      copyCohort%dbdeaddt                = this%dbdeaddt

      ! FIRE
      copyCohort%fraction_crown_burned   = this%fraction_crown_burned
      copyCohort%cambial_mort            = this%cambial_mort
      copyCohort%crownfire_mort          = this%crownfire_mort
      copyCohort%fire_mort               = this%fire_mort

      ! HYDRAULICS
      if (hlm_use_planthydro .eq. itrue) then
        call copyCohort%co_hydr%CopyCohortHydraulics(this%co_hydr)
      endif

    end subroutine Copy

    !===========================================================================

    subroutine FreeMemory(this)
      !
      ! DESCRIPTION:
      ! deallocates all dynamic memory and objects within the cohort structure
      ! DOES NOT deallocate the cohort structure itself
      !

      ! ARGUMENTS
      class(fates_cohort_type), intent(inout) :: this ! cohort object

      ! LOCALS:
      integer            :: istat ! return status code
      character(len=255) :: smsg  ! error message
 
      ! at this point, nothing should be pointing to current cohort
      if (hlm_use_planthydro .eq. itrue) then
        call this%co_hydr%DeAllocateHydrCohortArrays()
        deallocate(this%co_hydr)
      end if
 
      ! deallocate the cohort's PRT structures
      call this%prt%DeallocatePRTVartypes()
 
      ! Deallocate the PRT object
      deallocate(this%prt, stat=istat, errmsg=smsg)
      if (istat /= 0) then
        write(fates_log(),*) 'dealloc002: fail in deallocate(currentCohort%prt):'//trim(smsg)
        call endrun(msg=errMsg(sourcefile, __LINE__))
      endif

    end subroutine FreeMemory

    !===========================================================================
  
    subroutine InitPRTBoundaryConditions(this)      
      !
      ! DESCRIPTION:
      ! Set the boundary conditions that flow in an out of the PARTEH
      ! allocation hypotheses.  Each of these calls to "RegsterBC" are simply
      ! setting pointers.
      ! For instance, if the hypothesis wants to know what
      ! the DBH of the plant is, then we pass in the dbh as an argument (copyCohort%dbh),
      ! and also tell it which boundary condition we are talking about (which is
      ! defined by an integer index (ac_bc_inout_id_dbh)
      !
      ! Again, elaborated Example:
      ! "ac_bc_inout_id_dbh" is the unique integer that defines the object index
      ! for the allometric carbon "ac" boundary condition "bc" for DBH "dbh"
      ! that is classified as input and output "inout".
      ! See PRTAllometricCarbonMod.F90 to track its usage.
      ! bc_rval is used as the optional argument identifyer to specify a real
      ! value boundary condition.
      ! bc_ival is used as the optional argument identifyer to specify an integer
      ! value boundary condition.
      
      ! ARGUMENTS:
      class(fates_cohort_type), intent(inout), target :: this
      
      select case(hlm_parteh_mode)
      case (prt_carbon_allom_hyp)
   
        ! Register boundary conditions for the Carbon Only Allometric Hypothesis
  
        call this%prt%RegisterBCInOut(ac_bc_inout_id_dbh, bc_rval=this%dbh)
        call this%prt%RegisterBCInOut(ac_bc_inout_id_netdc, bc_rval=this%npp_acc)
        call this%prt%RegisterBCIn(ac_bc_in_id_cdamage, bc_ival=this%crowndamage)
        call this%prt%RegisterBCIn(ac_bc_in_id_pft, bc_ival=this%pft)
        call this%prt%RegisterBCIn(ac_bc_in_id_ctrim, bc_rval=this%canopy_trim)
        call this%prt%RegisterBCIn(ac_bc_in_id_lstat, bc_ival=this%status_coh)
        call this%prt%RegisterBCIn(ac_bc_in_id_efleaf, bc_rval = this%efleaf_coh)
        call this%prt%RegisterBCIn(ac_bc_in_id_effnrt, bc_rval = this%effnrt_coh)
        call this%prt%RegisterBCIn(ac_bc_in_id_efstem, bc_rval = this%efstem_coh)
        
      case (prt_cnp_flex_allom_hyp)
   
        ! Register boundary conditions for the CNP Allometric Hypothesis
   
        call this%prt%RegisterBCIn(acnp_bc_in_id_pft, bc_ival=this%pft)
        call this%prt%RegisterBCIn(acnp_bc_in_id_ctrim, bc_rval=this%canopy_trim)
        call this%prt%RegisterBCIn(acnp_bc_in_id_lstat, bc_ival=this%status_coh)
        call this%prt%RegisterBCIn(acnp_bc_in_id_efleaf, bc_rval = this%efleaf_coh)
        call this%prt%RegisterBCIn(acnp_bc_in_id_effnrt, bc_rval = this%effnrt_coh)
        call this%prt%RegisterBCIn(acnp_bc_in_id_efstem, bc_rval = this%efstem_coh)
        call this%prt%RegisterBCIn(acnp_bc_in_id_netdc, bc_rval=this%npp_acc)
  
        call this%prt%RegisterBCIn(acnp_bc_in_id_nc_repro, bc_rval=this%nc_repro)
        call this%prt%RegisterBCIn(acnp_bc_in_id_pc_repro, bc_rval=this%pc_repro)
        call this%prt%RegisterBCIn(acnp_bc_in_id_cdamage, bc_ival=this%crowndamage)
        
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_dbh, bc_rval=this%dbh)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_resp_excess, bc_rval=this%resp_excess_hold)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_l2fr, bc_rval=this%l2fr)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_cx_int, bc_rval=this%cx_int)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_emadcxdt, bc_rval=this%ema_dcxdt)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_cx0, bc_rval=this%cx0)
        
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_netdn, bc_rval=this%daily_n_gain)
        call this%prt%RegisterBCInOut(acnp_bc_inout_id_netdp, bc_rval=this%daily_p_gain)
        
        call this%prt%RegisterBCOut(acnp_bc_out_id_cefflux, bc_rval=this%daily_c_efflux)
        call this%prt%RegisterBCOut(acnp_bc_out_id_nefflux, bc_rval=this%daily_n_efflux)
        call this%prt%RegisterBCOut(acnp_bc_out_id_pefflux, bc_rval=this%daily_p_efflux)
        call this%prt%RegisterBCOut(acnp_bc_out_id_limiter, bc_ival=this%cnp_limiter)
        
      case DEFAULT
   
        write(fates_log(),*) 'You specified an unknown PRT module'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
   
      end select
   
    end subroutine InitPRTBoundaryConditions
   
    !===========================================================================

    subroutine UpdateCohortBioPhysRates(this)
      !
      ! DESCRIPTION:
      ! Update the four key biophysical rates of leaves based on the changes 
      ! in a cohort's leaf age proportions.
      !
      ! This should be called after growth.  Growth occurs
      ! after turnover and damage states are applied to the tree.
      ! Therefore, following growth, the leaf mass fractions
      ! of different age classes are unchanged until the next day.

      ! ARGUMENTS
      class(fates_cohort_type), intent(inout) :: this ! cohort object

      ! LOCAL VARIABLES
      real(r8) :: frac_leaf_aclass(max_nleafage)  ! fraction of leaves in each age-class
      integer  :: iage                            ! loop index for leaf ages
      integer  :: ipft                            ! plant functional type index

      ! First, calculate the fraction of leaves in each age class
      ! It is assumed that each class has the same proportion across leaf layers
      do iage = 1, nleafage
         frac_leaf_aclass(iage) = this%prt%GetState(leaf_organ,                &
          carbon12_element, iage)
      end do

      ! If there are leaves, then perform proportional weighting on the four rates
      ! We assume that leaf age does not effect the specific leaf area, so the mass
      ! fractions are applicable to these rates

      ipft = this%pft

      if (sum(frac_leaf_aclass(1:nleafage)) > nearzero .and.                   &
        hlm_use_sp .eq. ifalse) then

        frac_leaf_aclass(1:nleafage) = frac_leaf_aclass(1:nleafage)/           &
          sum(frac_leaf_aclass(1:nleafage))

        this%vcmax25top = sum(EDPftvarcon_inst%vcmax25top(ipft, 1:nleafage)*   &
          frac_leaf_aclass(1:nleafage))

        this%jmax25top = sum(param_derived%jmax25top(ipft, 1:nleafage)*        &
          frac_leaf_aclass(1:nleafage))

        this%tpu25top = sum(param_derived%tpu25top(ipft, 1:nleafage)*          &
          frac_leaf_aclass(1:nleafage))

        this%kp25top = sum(param_derived%kp25top(ipft, 1:nleafage)*            &
          frac_leaf_aclass(1:nleafage))

      else if (hlm_use_sp .eq. itrue) then
          
        this%vcmax25top = EDPftvarcon_inst%vcmax25top(ipft, 1)
        this%jmax25top = param_derived%jmax25top(ipft, 1)
        this%tpu25top = param_derived%tpu25top(ipft, 1)
        this%kp25top = param_derived%kp25top(ipft, 1)

      else
      
        this%vcmax25top = 0._r8
        this%jmax25top  = 0._r8
        this%tpu25top   = 0._r8
        this%kp25top    = 0._r8

      end if

    end subroutine UpdateCohortBioPhysRates

    !===========================================================================

    function CanUpperUnder(this) result(can_position)
      !
      ! DESCRIPTION:
      ! This simple function is used to determine if a cohort's crown position 
      !   is in the upper portion (ie the canopy) or the understory.  This 
      !   differentiation is only used for diagnostic purposes.  Functionally, 
      !   the model uses the canopy layer position, which may have more than 
      !   two layers at any given time. Utlimately, every plant that is not in 
      !   the top layer (canopy), is considered understory.
      !
   
      ! ARGUMENTS:
      class(fates_cohort_type) :: this         ! current cohort of interest
      integer                  :: can_position ! canopy position 
      
      if (this%canopy_layer == 1)then
        can_position = ican_upper
      else
        can_position = ican_ustory
      end if
      
    end function CanUpperUnder
   
    !===========================================================================

    subroutine Dump(this)
      !
      !  DESCRIPTION:
      !  Print out attributes of a cohort
      !
   
      ! ARGUMENTS:
      class(fates_cohort_type), intent(in), target :: this
   
      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) ' Dumping Cohort Information             '
      write(fates_log(),*) '----------------------------------------'
      write(fates_log(),*) 'cohort%pft                    = ', this%pft
      write(fates_log(),*) 'cohort%n                      = ', this%n                         
      write(fates_log(),*) 'cohort%dbh                    = ', this%dbh                                        
      write(fates_log(),*) 'cohort%height                 = ', this%height
      write(fates_log(),*) 'cohort%crowndamage            = ', this%crowndamage
      write(fates_log(),*) 'cohort%coage                  = ', this%coage
      write(fates_log(),*) 'cohort%l2fr                   = ', this%l2fr
      write(fates_log(),*) 'leaf carbon                   = ', this%prt%GetState(leaf_organ,carbon12_element) 
      write(fates_log(),*) 'fineroot carbon               = ', this%prt%GetState(fnrt_organ,carbon12_element) 
      write(fates_log(),*) 'sapwood carbon                = ', this%prt%GetState(sapw_organ,carbon12_element) 
      write(fates_log(),*) 'structural (dead) carbon      = ', this%prt%GetState(struct_organ,carbon12_element) 
      write(fates_log(),*) 'storage carbon                = ', this%prt%GetState(store_organ,carbon12_element) 
      write(fates_log(),*) 'reproductive carbon           = ', this%prt%GetState(repro_organ,carbon12_element) 
      write(fates_log(),*) 'cohort%g_sb_laweight          = ', this%g_sb_laweight
      write(fates_log(),*) 'cohort%leaf_cost              = ', this%leaf_cost
      write(fates_log(),*) 'cohort%canopy_layer           = ', this%canopy_layer
      write(fates_log(),*) 'cohort%canopy_layer_yesterday = ', this%canopy_layer_yesterday
      write(fates_log(),*) 'cohort%nv                     = ', this%nv
      write(fates_log(),*) 'cohort%status_coh             = ', this%status_coh
      write(fates_log(),*) 'co%status_coh                 = ', this%status_coh
      write(fates_log(),*) 'co%efleaf_coh                 = ', this%efleaf_coh
      write(fates_log(),*) 'co%effnrt_coh                 = ', this%effnrt_coh
      write(fates_log(),*) 'co%efstem_coh                 = ', this%efstem_coh
      write(fates_log(),*) 'cohort%canopy_trim            = ', this%canopy_trim
      write(fates_log(),*) 'cohort%excl_weight            = ', this%excl_weight               
      write(fates_log(),*) 'cohort%prom_weight            = ', this%prom_weight               
      write(fates_log(),*) 'cohort%size_class             = ', this%size_class
      write(fates_log(),*) 'cohort%size_by_pft_class      = ', this%size_by_pft_class
      write(fates_log(),*) 'cohort%coage_class            = ', this%coage_class
      write(fates_log(),*) 'cohort%coage_by_pft_class     = ', this%coage_by_pft_class
      write(fates_log(),*) 'cohort%gpp_acc_hold           = ', this%gpp_acc_hold
      write(fates_log(),*) 'cohort%gpp_acc                = ', this%gpp_acc
      write(fates_log(),*) 'cohort%gpp_tstep              = ', this%gpp_tstep
      write(fates_log(),*) 'cohort%npp_acc_hold           = ', this%npp_acc_hold
      write(fates_log(),*) 'cohort%npp_acc                = ', this%npp_acc
      write(fates_log(),*) 'cohort%resp_m_tstep           = ', this%resp_m_tstep
      write(fates_log(),*) 'cohort%resp_m_acc             = ', this%resp_m_acc
      write(fates_log(),*) 'cohort%resp_m_acc_hold        = ', this%resp_m_acc_hold
      write(fates_log(),*) 'cohort%resp_g_acc_hold        = ', this%resp_g_acc_hold
      write(fates_log(),*) 'cohort%rdark                  = ', this%rdark
      write(fates_log(),*) 'cohort%livestem_mr            = ', this%livestem_mr
      write(fates_log(),*) 'cohort%livecroot_mr           = ', this%livecroot_mr
      write(fates_log(),*) 'cohort%froot_mr               = ', this%froot_mr
      write(fates_log(),*) 'cohort%dgmort                 = ', this%dgmort
      write(fates_log(),*) 'cohort%treelai                = ', this%treelai
      write(fates_log(),*) 'cohort%treesai                = ', this%treesai
      write(fates_log(),*) 'cohort%c_area                 = ', this%c_area
      write(fates_log(),*) 'cohort%cmort                  = ', this%cmort
      write(fates_log(),*) 'cohort%bmort                  = ', this%bmort
      write(fates_log(),*) 'cohort%smort                  = ', this%smort
      write(fates_log(),*) 'cohort%asmort                 = ', this%asmort
      write(fates_log(),*) 'cohort%dgmort                 = ', this%dgmort
      write(fates_log(),*) 'cohort%hmort                  = ', this%hmort
      write(fates_log(),*) 'cohort%frmort                 = ', this%frmort
      write(fates_log(),*) 'cohort%asmort                 = ', this%asmort
      write(fates_log(),*) 'cohort%lmort_direct           = ', this%lmort_direct
      write(fates_log(),*) 'cohort%lmort_collateral       = ', this%lmort_collateral
      write(fates_log(),*) 'cohort%lmort_infra            = ', this%lmort_infra
      write(fates_log(),*) 'cohort%isnew                  = ', this%isnew
      write(fates_log(),*) 'cohort%dndt                   = ', this%dndt
      write(fates_log(),*) 'cohort%dhdt                   = ', this%dhdt
      write(fates_log(),*) 'cohort%ddbhdt                 = ', this%ddbhdt
      write(fates_log(),*) 'cohort%dbdeaddt               = ', this%dbdeaddt
      write(fates_log(),*) 'cohort%fraction_crown_burned  = ', this%fraction_crown_burned
      write(fates_log(),*) 'cohort%fire_mort              = ', this%fire_mort
      write(fates_log(),*) 'cohort%crownfire_mort         = ', this%crownfire_mort
      write(fates_log(),*) 'cohort%cambial_mort           = ', this%cambial_mort
      write(fates_log(),*) 'cohort%size_class             = ', this%size_class
      write(fates_log(),*) 'cohort%size_by_pft_class      = ', this%size_by_pft_class
   
      if (associated(this%co_hydr)) call this%co_hydr%Dump()
   
      write(fates_log(),*) '----------------------------------------'
   
    end subroutine Dump
   
    !===========================================================================
   
end module FatesCohortMod
