module PRTAllometricCNPMod

   ! ------------------------------------------------------------------------------------
   !
   ! This module contains all of the specific functions and types for
   ! Plant Allocation and Reactive Transport Extensible Hypotheses (PARTEH)
   !
   ! Carbon-Nitrogen-Phosphorus (CNP) Prioritized Allometric Allocations
   !
   ! Ryan Knox Aug 2018
   !
   ! ------------------------------------------------------------------------------------

  use PRTGenericMod , only  : prt_global_type
  use PRTGenericMod , only  : prt_global
  use PRTGenericMod , only  : prt_vartypes
  use PRTGenericMod , only  : carbon12_element
  use PRTGenericMod , only  : nitrogen_element
  use PRTGenericMod , only  : phosphorus_element
  use PRTGenericMod , only  : max_nleafage
  use PRTGenericMod , only  : l2fr_min
  use PRTGenericMod , only  : leaf_organ
  use PRTGenericMod , only  : fnrt_organ
  use PRTGenericMod , only  : sapw_organ
  use PRTGenericMod , only  : store_organ
  use PRTGenericMod , only  : repro_organ
  use PRTGenericMod , only  : struct_organ
  use PRTGenericMod , only  : num_organ_types
  use PRTGenericMod , only  : prt_cnp_flex_allom_hyp
  use PRTGenericMod , only  : StorageNutrientTarget
  
  use FatesAllometryMod   , only : bleaf
  use FatesAllometryMod   , only : bsap_allom
  use FatesAllometryMod   , only : bfineroot
  use FatesAllometryMod   , only : bstore_allom
  use FatesAllometryMod   , only : bdead_allom
  use FatesAllometryMod   , only : bbgw_allom
  use FatesAllometryMod   , only : bagw_allom
  use FatesAllometryMod   , only : h_allom
  use FatesAllometryMod   , only : CheckIntegratedAllometries

  use FatesGlobals        , only : endrun => fates_endrun
  use FatesGlobals        , only : fates_log
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use FatesConstantsMod   , only : r8 => fates_r8
  use FatesConstantsMod   , only : i4 => fates_int
  use FatesConstantsMod   , only : calloc_abs_error
  use FatesConstantsMod   , only : rsnbl_math_prec
  use FatesConstantsMod   , only : years_per_day
  use FatesConstantsMod   , only : mm_per_cm
  use FatesIntegratorsMod , only : RKF45
  use FatesIntegratorsMod , only : Euler
  use FatesConstantsMod   , only : calloc_abs_error
  use FatesConstantsMod   , only : nearzero
  use FatesConstantsMod   , only : itrue
  use FatesConstantsMod   , only : fates_unset_r8
  use FatesConstantsMod   , only : fates_unset_int
  use FatesConstantsMod   , only : sec_per_day
  use FatesConstantsMod   , only : TRS_regeneration
  use FatesConstantsMod   , only : default_regeneration
  use FatesConstantsMod   , only : TRS_no_seedling_dyn
  use FatesConstantsMod   , only : min_max_dbh_for_trees
  use PRTParametersMod    , only : prt_params
  use FatesConstantsMod   , only : leaves_on,leaves_off
  use FatesConstantsMod   , only : leaves_shedding
  use EDParamsMod         , only : p_uptake_mode
  use EDParamsMod         , only : n_uptake_mode
  use FatesConstantsMod   , only : prescribed_p_uptake
  use FatesConstantsMod   , only : prescribed_n_uptake
  use EDPftvarcon, only : EDPftvarcon_inst
  use EDParamsMod       , only : regeneration_model


  
  implicit none
  private

  ! -------------------------------------------------------------------------------------
  !
  ! Define the state variables for this specific hypothesis. Give them units and define
  ! the indices that correspond with the generic classifications of PRT variables
  !
  ! -------------------------------------------------------------------------------------

  integer, parameter :: leaf_c_id   = 1       ! leaf carbon index
  integer, parameter :: fnrt_c_id   = 2       ! fine-root carbon index
  integer, parameter :: sapw_c_id   = 3       ! sapwood carbon index
  integer, parameter :: store_c_id  = 4       ! storage carbon index
  integer, parameter :: repro_c_id  = 5       ! reproductive carbon index
  integer, parameter :: struct_c_id = 6       ! structural carbon index

  integer, parameter :: leaf_n_id   = 7
  integer, parameter :: fnrt_n_id   = 8
  integer, parameter :: sapw_n_id   = 9
  integer, parameter :: store_n_id  = 10
  integer, parameter :: repro_n_id  = 11
  integer, parameter :: struct_n_id = 12

  integer, parameter :: leaf_p_id   = 13
  integer, parameter :: fnrt_p_id   = 14
  integer, parameter :: sapw_p_id   = 15
  integer, parameter :: store_p_id  = 16
  integer, parameter :: repro_p_id  = 17
  integer, parameter :: struct_p_id = 18

  ! Total number of state variables
  integer, parameter :: num_vars    = 18
  
  
  ! Global identifiers for the two stoichiometry values
  integer,public, parameter :: stoich_growth_min = 1     ! Flag for stoichiometry associated with
                                                         ! minimum needed for growth
  
  ! This is deprecated until a reasonable hypothesis is in place (RGK)
  integer,public, parameter :: stoich_max = 2            ! Flag for stoichiometry associated with 
                                                         ! maximum for that organ
  
  ! This is the ordered list of organs used in this module
  ! -------------------------------------------------------------------------------------

  integer, parameter                        :: num_organs = 6
  
  ! Converting from local to global organ id
  integer, parameter,dimension(num_organs) :: l2g_organ_list = & 
     [leaf_organ, fnrt_organ, sapw_organ, store_organ, repro_organ, struct_organ]


  
  ! These are local indices associated with organs and quantities
  ! that can be integrated (namely, growth respiration during stature growth
  ! and dbh)
  
  integer, parameter :: leaf_id = 1
  integer, parameter :: fnrt_id = 2
  integer, parameter :: sapw_id = 3
  integer, parameter :: store_id = 4
  integer, parameter :: repro_id = 5
  integer, parameter :: struct_id = 6
  integer, parameter :: dbh_id = 7

  integer, parameter :: num_intgr_vars = 7


  
  ! -------------------------------------------------------------------------------------
  ! Input/Output Boundary Indices (These are public, and therefore
  !       each boundary condition across all modules must
  !       have a unique name !!!!)
  ! They are used in the routine, and also changed in the routine before
  ! being passed back
  ! -------------------------------------------------------------------------------------


  integer, public, parameter :: acnp_bc_inout_id_dbh        = 1  ! Plant DBH
  integer, public, parameter :: acnp_bc_inout_id_resp_excess = 2  ! Respiration of excess storage
  integer, public, parameter :: acnp_bc_inout_id_l2fr       = 3  ! leaf 2 fineroot scalar, this
                                                                 ! is dynamic with CNP
  integer, public, parameter :: acnp_bc_inout_id_netdn    = 4 ! Index for the net daily NH4 input BC
  integer, public, parameter :: acnp_bc_inout_id_netdp    = 5 ! Index for the net daily P input BC
  integer, public, parameter :: acnp_bc_inout_id_cx_int   = 6 ! Index for the EMA log storage ratio max(N,P)/C
  integer, public, parameter :: acnp_bc_inout_id_cx0      = 7 ! Index for the previous step's log storage ratio max(N,P)/C
  integer, public, parameter :: acnp_bc_inout_id_emadcxdt = 8 ! Index for the EMA log storage ratio derivative d max(NP)/C dt
  integer, public, parameter :: num_bc_inout              = 8

  ! -------------------------------------------------------------------------------------
  ! Input only Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: acnp_bc_in_id_pft      =  1 ! Index for the PFT input BC
  integer, public, parameter :: acnp_bc_in_id_ctrim    =  2 ! Index for the canopy trim function
  integer, public, parameter :: acnp_bc_in_id_lstat    =  3 ! phenology status logical
  integer, public, parameter :: acnp_bc_in_id_netdc    =  4 ! Index for the net daily C input BC
  integer, public, parameter :: acnp_bc_in_id_nc_repro =  5
  integer, public, parameter :: acnp_bc_in_id_pc_repro =  6
  integer, public, parameter :: acnp_bc_in_id_cdamage  =  7 ! Index for the crowndamage input BC
  integer, public, parameter :: acnp_bc_in_id_efleaf   =  8 ! Leaf elongation factor
  integer, public, parameter :: acnp_bc_in_id_effnrt   =  9 ! Fine-root "elongation factor"
  integer, public, parameter :: acnp_bc_in_id_efstem   = 10 ! Stem "elongation factor"
  integer, parameter         :: num_bc_in              = 10

  ! -------------------------------------------------------------------------------------
  ! Output Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: acnp_bc_out_id_cefflux = 1  ! Daily exudation of C  [kg]
  integer, public, parameter :: acnp_bc_out_id_nefflux = 2  ! Daily exudation of N  [kg]
  integer, public, parameter :: acnp_bc_out_id_pefflux = 3  ! Daily exudation of P  [kg]
  integer, public, parameter :: acnp_bc_out_id_limiter = 4  ! The minimum of the Nutrient ratio over c ratio
  
  integer, parameter         :: num_bc_out             = 4  ! Total number of


  ! Indices for parameters passed to the integrator
  integer,private, parameter :: intgr_parm_ctrim   = 1
  integer,private, parameter :: intgr_parm_pft     = 2
  integer,private, parameter :: intgr_parm_l2fr    = 3
  integer,private, parameter :: intgr_parm_cdamage = 4
  integer,private, parameter :: intgr_parm_efleaf  = 5
  integer,private, parameter :: intgr_parm_effnrt  = 6
  integer,private, parameter :: intgr_parm_efstem  = 7
  integer,private, parameter :: num_intgr_parm     = 7
  
  ! -------------------------------------------------------------------------------------
  ! Define the size of the coorindate vector.  For this hypothesis, there is only
  ! one pool per each species x organ combination, except for leaves (WHICH HAVE AGE)
  ! icd refers to the first index in the leaf array, which is the youngest.  Growth
  ! and allocation only happens in the youngest bin by definition
  ! -------------------------------------------------------------------------------------
  integer, parameter :: icd = 1

  ! These constants define different methods of dealing with excess carbon at
  ! the end of the allocation process, assuming that N or P is limiting growth.
  ! You can either exude (send to soil), retain (grow storage without limit),
  ! or burn (as respiration to the atm).
  integer, parameter :: exude_c_store_overflow = 1
  integer, parameter :: retain_c_store_overflow = 2
  integer, parameter :: burn_c_store_overflow = 3
  integer, parameter :: store_c_overflow = burn_c_store_overflow

  ! These constants define if/how growth is limited by
  ! one of the 3 chemical species, 0 indicates there is some
  ! degree of co limitation
  integer, parameter  :: cnp_limited = 0
  integer, parameter  :: c_limited = 1
  integer, parameter  :: n_limited = 2
  integer, parameter  :: p_limited = 3

  ! Flags to select using the equivalent carbon method of co-limitation,
  ! or to just grow with available carbon and let it fix itself on the
  ! next step
  
  integer, parameter  :: grow_lim_conly = 1  ! Just use C to decide stature on this step
  integer, parameter  :: grow_lim_estNP = 2  ! Estimate equivalent C from N and P
  integer, parameter  :: grow_lim_type = grow_lim_estNP
    
  ! Following growth, if desired, you can prioritize that 
  ! reproductive tissues get balanced CNP
  logical, parameter :: prioritize_repro_nutr_growth = .true.
  
  ! If this parameter is true, then the fine-root l2fr optimization
  ! scheme will remove biomass from roots without restriction if the
  ! l2fr is getting smaller.
  logical, parameter :: use_unrestricted_contraction = .true.

  ! -------------------------------------------------------------------------------------
  ! This is the core type that holds this specific
  ! plant reactive transport (PRT) module
  ! -------------------------------------------------------------------------------------
  type, public, extends(prt_vartypes) :: cnp_allom_prt_vartypes

  contains

     procedure :: DailyPRT     => DailyPRTAllometricCNP
     procedure :: FastPRT      => FastPRTAllometricCNP
     procedure :: GetNutrientTarget => GetNutrientTargetCNP
     
     ! Extended functions specific to Allometric CNP
     procedure :: CNPPrioritizedReplacement
     procedure :: CNPStatureGrowth
     procedure :: EstimateGrowthNC
     procedure :: CNPAdjustFRootTargets
     procedure :: CNPAllocateRemainder
     procedure :: GetDeficit
     procedure :: TrimFineRoot
  end type cnp_allom_prt_vartypes


  ! ------------------------------------------------------------------------------------
  !
  ! This next class is an extention of the base instance that maps state variables
  !      to the outside model.
  !
  ! This is the instance of the mapping table and variable definitions
  ! this is only allocated once per node
  ! ------------------------------------------------------------------------------------

   class(prt_global_type), public, target, allocatable :: prt_global_acnp

   character(len=*), parameter, private :: sourcefile = __FILE__
   logical, parameter :: debug = .false.
   
   public :: InitPRTGlobalAllometricCNP
   

contains


  subroutine InitPRTGlobalAllometricCNP()

     ! ----------------------------------------------------------------------------------
     ! Initialize and populate the general mapping table that
     ! organizes the specific variables in this module to
     ! pre-ordained groups, so they can be used to inform
     ! the rest of the model
     !
     ! This routine is not part of the sp_pool_vartypes class
     ! because it is the same for all plants and we need not
     ! waste memory on it.
     ! -----------------------------------------------------------------------------------

     integer :: nleafage
     integer :: istat
     character(len=255) :: smsg

     allocate(prt_global_acnp, stat=istat, errmsg=smsg)
     if (istat/=0) call endrun(msg='allocate stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
     allocate(prt_global_acnp%state_descriptor(num_vars), stat=istat, errmsg=smsg)
     if (istat/=0) call endrun(msg='allocate stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
     
     prt_global_acnp%hyp_name = 'Allometric Flexible C+N+P'

     prt_global_acnp%hyp_id = prt_cnp_flex_allom_hyp

     call prt_global_acnp%ZeroGlobal()

     ! The number of leaf age classes can be determined from the parameter file,
     ! notably the size of the leaf-longevity parameter's second dimension.
     ! This is the same value in FatesInterfaceMod.F90

     nleafage = size(prt_params%leaf_long,dim=2)

     if(nleafage>max_nleafage) then
        write(fates_log(),*) 'The allometric carbon PARTEH hypothesis'
        write(fates_log(),*) 'sets a maximum number of leaf age classes'
        write(fates_log(),*) 'used for scratch space. The model wants'
        write(fates_log(),*) 'exceed that. Simply increase max_nleafage'
        write(fates_log(),*) 'found in parteh/PRTAllometricCarbonMod.F90'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     call prt_global_acnp%RegisterVarInGlobal(leaf_c_id,'Leaf Carbon','leaf_c',leaf_organ,carbon12_element,nleafage)
     call prt_global_acnp%RegisterVarInGlobal(fnrt_c_id,'Fine Root Carbon','fnrt_c',fnrt_organ,carbon12_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(sapw_c_id,'Sapwood Carbon','sapw_c',sapw_organ,carbon12_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(store_c_id,'Storage Carbon','store_c',store_organ,carbon12_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(struct_c_id,'Structural Carbon','struct_c',struct_organ,carbon12_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(repro_c_id,'Reproductive Carbon','repro_c',repro_organ,carbon12_element,icd)

     call prt_global_acnp%RegisterVarInGlobal(leaf_n_id,'Leaf Nitrogen','leaf_n',leaf_organ,nitrogen_element,nleafage)
     call prt_global_acnp%RegisterVarInGlobal(fnrt_n_id,'Fine Root Nitrogen','fnrt_n',fnrt_organ,nitrogen_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(sapw_n_id,'Sapwood Nitrogen','sapw_n',sapw_organ,nitrogen_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(store_n_id,'Storage Nitrogen','store_n',store_organ,nitrogen_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(struct_n_id,'Structural Nitrogen','struct_n',struct_organ,nitrogen_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(repro_n_id,'Reproductive Nitrogen','repro_n',repro_organ,nitrogen_element,icd)

     call prt_global_acnp%RegisterVarInGlobal(leaf_p_id,'Leaf Phosphorus','leaf_p',leaf_organ,phosphorus_element,nleafage)
     call prt_global_acnp%RegisterVarInGlobal(fnrt_p_id,'Fine Root Phosphorus','fnrt_p',fnrt_organ,phosphorus_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(sapw_p_id,'Sapwood Phosphorus','sapw_p',sapw_organ,phosphorus_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(store_p_id,'Storage Phosphorus','store_p',store_organ,phosphorus_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(struct_p_id,'Structural Phosphorus','struct_p',struct_organ,phosphorus_element,icd)
     call prt_global_acnp%RegisterVarInGlobal(repro_p_id,'Reproductive Phosphorus','repro_p',repro_organ,phosphorus_element,icd)


     ! Set some of the array sizes for input and output boundary conditions
     prt_global_acnp%num_bc_in    = num_bc_in
     prt_global_acnp%num_bc_out   = num_bc_out
     prt_global_acnp%num_bc_inout = num_bc_inout
     prt_global_acnp%num_vars     = num_vars

     ! Have the global generic pointer, point to this hypothesis' object
     prt_global => prt_global_acnp

     return
  end subroutine InitPRTGlobalAllometricCNP


  ! =====================================================================================


  subroutine DailyPRTAllometricCNP(this,phase)

    class(cnp_allom_prt_vartypes) :: this
    integer,intent(in)           :: phase         ! the phase splits the routine into parts
    ! note that phasing is used primarily to
    ! accomodate the damage module. Damage
    ! and nutrient cycling are not yet compatable though
    ! hence, we simply return from any phase but phase 1

    
    ! Pointers to in-out bcs
    real(r8),pointer :: dbh          ! Diameter at breast height [cm]
    real(r8),pointer :: resp_excess   ! Respiration of any un-allocatable C
    real(r8),pointer :: l2fr         ! Leaf to fineroot ratio of target biomass
    
    ! Input only bcs
    integer          :: ipft         ! Plant Functional Type index
    real(r8)         :: c_gain       ! Daily carbon balance for this cohort [kgC]
    real(r8),pointer :: n_gain       ! Daily nitrogen uptake through fine-roots [kgN]
    real(r8),pointer :: p_gain       ! Daily phosphorus uptake through fine-roots [kgN]
    real(r8)         :: canopy_trim  ! The canopy trimming function [0-1]
    integer          :: crown_damage ! which crown damage clas
    real(r8)         :: elongf_leaf  ! Leaf elongation factor [0-1]
    real(r8)         :: elongf_fnrt  ! Fine-root "elongation factor" [0-1]
    real(r8)         :: elongf_stem  ! Stem "elongation factor" [0-1]

    ! Pointers to output bcs
    real(r8),pointer :: c_efflux   ! Total plant efflux of carbon (kgC)
    real(r8),pointer :: n_efflux   ! Total plant efflux of nitrogen (kgN)
    real(r8),pointer :: p_efflux   ! Total plant efflux of phosphorus (kgP)

    ! Allometry targets (kg/plant) and (kg/cm/plant)
    real(r8), dimension(num_organ_types) :: target_c, target_dcdd
    ! Initial states (for accounting) (kg/plant)
    real(r8), dimension(num_organ_types) :: state_c0, state_n0, state_p0
    ! Allometry partial targets
    real(r8) :: agw_c_target,agw_dcdd_target
    real(r8) :: bgw_c_target,bgw_dcdd_target
    real(r8) :: sapw_area
    real(r8) :: store_flux
    integer :: i       ! generic organ loop index
    integer :: i_org   ! organ index
    integer :: i_var   ! variable index

    ! These are daily mass gains, frozen in time, not drawn from, and thus
    ! these are only used for evaluating mass balancing at the end
    real(r8) :: dbh0
    real(r8) :: c_gain0
    real(r8) :: n_gain0
    real(r8) :: p_gain0
    real(r8) :: resp_excess0

    ! Used for mass checking, total mass allocated based
    ! on change in the states, should match gain0's
    real(r8) :: allocated_c
    real(r8) :: allocated_n
    real(r8) :: allocated_p
    real(r8) :: sum_c ! error checking sum


    ! Phasing is only used to accomodate the
    ! damage module. Since this is incompatible with CNP
    ! Ignore all subsequent calls after the first
    if (phase.ne.1) return
    
    
    ! In/out boundary conditions
    resp_excess => this%bc_inout(acnp_bc_inout_id_resp_excess)%rval
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    dbh0        =  dbh
    l2fr        => this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    n_gain      => this%bc_inout(acnp_bc_inout_id_netdn)%rval
    p_gain      => this%bc_inout(acnp_bc_inout_id_netdp)%rval


    ! Assume that there is no other source of excess respiration
    ! so it is safe to zero it. In the third stage we will
    ! decide if this should be updated
    resp_excess  = 0._r8
    resp_excess0 = resp_excess

    
    ! integrator variables

    ! Copy the input only boundary conditions into readable local variables
    ! We don't use pointers, because inputs should be intent in only
    ! Also, we save the initial values of many of these BC's
    ! for checking and resetting if needed
    ! -----------------------------------------------------------------------------------
    c_gain      = this%bc_in(acnp_bc_in_id_netdc)%rval
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    crown_damage = this%bc_in(acnp_bc_in_id_cdamage)%ival
    elongf_leaf  = this%bc_in(acnp_bc_in_id_efleaf)%rval
    elongf_fnrt  = this%bc_in(acnp_bc_in_id_effnrt)%rval
    elongf_stem  = this%bc_in(acnp_bc_in_id_efstem)%rval

    ! If either n or p uptake is in prescribed mode
    ! set the gains to something massive. 1 kilo of pure
    ! nutrient should be wayyy more than enough
    if(n_uptake_mode.eq.prescribed_n_uptake) then
       n_gain = 1.e3_r8
    end if
    if(p_uptake_mode.eq.prescribed_p_uptake) then
       p_gain = 1.e3_r8
    end if
    
    n_gain0      = n_gain
    p_gain0      = p_gain
    c_gain0      = c_gain


    ! Calculate Carbon allocation targets
    ! -----------------------------------------------------------------------------------

    ! Set carbon targets based on the plant's current stature
    target_c(:) = fates_unset_r8
    target_dcdd(:) = fates_unset_r8
    call bsap_allom(dbh,ipft,crown_damage,canopy_trim,elongf_stem,sapw_area,target_c(sapw_organ),target_dcdd(sapw_organ))
    call bagw_allom(dbh,ipft,crown_damage,elongf_stem,agw_c_target,agw_dcdd_target)
    call bbgw_allom(dbh,ipft,elongf_stem,bgw_c_target,bgw_dcdd_target)
    call bdead_allom(agw_c_target,bgw_c_target,target_c(sapw_organ),ipft,target_c(struct_organ), &
                     agw_dcdd_target,bgw_dcdd_target,target_dcdd(sapw_organ),target_dcdd(struct_organ))
    call bleaf(dbh,ipft,crown_damage,canopy_trim, elongf_leaf, target_c(leaf_organ), target_dcdd(leaf_organ))
    call bfineroot(dbh,ipft,canopy_trim, l2fr, elongf_fnrt, target_c(fnrt_organ), target_dcdd(fnrt_organ))
    call bstore_allom(dbh,ipft,crown_damage, canopy_trim, target_c(store_organ), target_dcdd(store_organ))
    target_c(repro_organ) = 0._r8
    target_dcdd(repro_organ) = 0._r8

    ! ===================================================================================
    ! Step 1: Evaluate nutrient storage in the plant. Depending on how low
    ! these stores are, we will move proportionally more or less of the daily carbon
    ! gain to increase the target fine-root biomass, fill up to target
    ! and then attempt to get them up to stoichiometry targets.
    ! ===================================================================================

    
    
    ! Remember the original C,N,P states to help with final
    ! evaluation of how much was allocated
    ! -----------------------------------------------------------------------------------

    do i = 1,num_organs
       i_org = l2g_organ_list(i) ! global index from PRTGeneric
       i_var = prt_global%sp_organ_map(i_org,carbon12_element)
       state_c0(i_org)  = this%variables(i_var)%val(1)
       i_var = prt_global%sp_organ_map(i_org,nitrogen_element)
       state_n0(i_org)  = this%variables(i_var)%val(1)
       i_var = prt_global%sp_organ_map(i_org,phosphorus_element)
       state_p0(i_org)  =  this%variables(i_var)%val(1)
    end do
    
    
    ! Output only boundary conditions
    c_efflux    => this%bc_out(acnp_bc_out_id_cefflux)%rval;  c_efflux = 0._r8
    n_efflux    => this%bc_out(acnp_bc_out_id_nefflux)%rval;  n_efflux = 0._r8
    p_efflux    => this%bc_out(acnp_bc_out_id_pefflux)%rval;  p_efflux = 0._r8

    ! ===================================================================================
    ! Step 0.  Transfer all stored nutrient into the daily uptake pool. Also
    !          transfer C storage that is above the target (ie transfer overflow)
    ! ===================================================================================

    ! Put overflow storage into the net daily pool
    store_flux = max(0._r8, this%variables(store_c_id)%val(1) - target_c(store_organ))
    c_gain = c_gain + store_flux
    this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) - store_flux

    n_gain = n_gain + sum(this%variables(store_n_id)%val(:))
    this%variables(store_n_id)%val(:) = 0._r8

    p_gain = p_gain + sum(this%variables(store_p_id)%val(:))
    this%variables(store_p_id)%val(:) = 0._r8
    
    
    ! ===================================================================================
    ! Step 2.  Prioritized allocation to replace tissues from turnover, and/or pay
    ! any un-paid maintenance respiration from storage.
    ! ===================================================================================
    
    call this%CNPPrioritizedReplacement(c_gain, n_gain, p_gain, target_c)
       
    sum_c = 0._r8
    do i = 1,num_organs
       i_org = l2g_organ_list(i)
       i_var = prt_global%sp_organ_map(i_org,carbon12_element)
       sum_c = sum_c+this%variables(i_var)%val(1)
    end do
    if( abs((c_gain0-c_gain) - &
            (sum_c-sum(state_c0(:),dim=1))) >calloc_abs_error ) then
       write(fates_log(),*) 'Carbon not balancing I'
       do i  = 1,num_organs
          i_org = l2g_organ_list(i)
          write(fates_log(),*) 'c: ',this%variables(prt_global%sp_organ_map(i_org,carbon12_element))%val(1)
       end do
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    ! ===================================================================================
    ! Step 3. Grow out the stature of the plant by allocating to tissues beyond
    ! current targets. 
    ! Attempts have been made to get all pools and species closest to allometric
    ! targets based on prioritized relative demand and allometry functions.
    ! ===================================================================================

    call this%CNPStatureGrowth(c_gain, n_gain, p_gain, target_c, target_dcdd)

    sum_c = 0._r8
    do i = 1,num_organs
       i_org = l2g_organ_list(i)
       i_var = prt_global%sp_organ_map(i_org,carbon12_element)
       sum_c = sum_c+this%variables(i_var)%val(1)
    end do
    if( abs((c_gain0-c_gain) - &
            (sum_c-sum(state_c0(:),dim=1))) >calloc_abs_error ) then
       write(fates_log(),*) 'Carbon not balancing II'
       do i  = 1,num_organs
          i_org = l2g_organ_list(i)
          write(fates_log(),*) 'c: ',this%variables(prt_global%sp_organ_map(i_org,carbon12_element))%val(1)
       end do
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! ===================================================================================
    ! Step 3. 
    ! At this point, at least 1 of the 3 resources have been used up.
    ! Allocate the remaining resources, or as a last resort, efflux them.
    ! ===================================================================================

    call this%CNPAllocateRemainder(c_gain, n_gain, p_gain, &
         c_efflux, n_efflux, p_efflux,target_c,target_dcdd)


    if(n_uptake_mode.ne.prescribed_n_uptake) then
       if( abs(n_gain) > 0.1_r8*calloc_abs_error) then
          write(fates_log(),*) 'Allocation scheme should had used up all mass gain pools'
          write(fates_log(),*) 'Any mass that cannot be allocated should be effluxed'
          write(fates_log(),*) 'n_gain: ',n_gain
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    if(p_uptake_mode.ne.prescribed_p_uptake) then
       if( abs(p_gain) > 0.01_r8*calloc_abs_error) then
          write(fates_log(),*) 'Allocation scheme should had used up all mass gain pools'
          write(fates_log(),*) 'Any mass that cannot be allocated should be effluxed'
          write(fates_log(),*) 'p_gain: ',p_gain
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    if( abs(c_gain) > calloc_abs_error) then
       write(fates_log(),*) 'Allocation scheme should had used up all mass gain pools'
       write(fates_log(),*) 'Any mass that cannot be allocated should be effluxed'
       write(fates_log(),*) 'c_gain: ',c_gain
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Perform a final tally on what was used (allocated)
    ! Since this is also a check against what was available
    ! we include what is lost through respiration of excess storage
    
    allocated_c = (resp_excess-resp_excess0) + c_efflux
    allocated_n = n_efflux
    allocated_p = p_efflux
    
    ! Update the allocation flux diagnostic arrays for each 3 elements
    do i = 1,num_organs

       i_org = l2g_organ_list(i)

       i_var = prt_global%sp_organ_map(i_org,carbon12_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (this%variables(i_var)%val(1) - state_c0(i_org))

       allocated_c = allocated_c + (this%variables(i_var)%val(1) - state_c0(i_org))
       
       i_var = prt_global%sp_organ_map(i_org,nitrogen_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (this%variables(i_var)%val(1) - state_n0(i_org))

       allocated_n = allocated_n + (this%variables(i_var)%val(1) - state_n0(i_org))
       
       i_var = prt_global%sp_organ_map(i_org,phosphorus_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (this%variables(i_var)%val(1) - state_p0(i_org))

       allocated_p = allocated_p + (this%variables(i_var)%val(1) - state_p0(i_org))
       
    end do
    
    if(debug) then

       ! Error Check: Do a final balance between how much mass
       ! we had to work with, and how much was allocated
       
       if ( abs(allocated_c - (c_gain0-c_gain)) > calloc_abs_error .or. & 
            abs(allocated_n - (n_gain0-n_gain)) > calloc_abs_error .or. &
            abs(allocated_p - (p_gain0-p_gain)) > calloc_abs_error ) then
          write(fates_log(),*) 'CNP allocation scheme did not balance mass.'
          write(fates_log(),*) 'c_gain0: ',c_gain0,' allocated_c: ',allocated_c,resp_excess,resp_excess0,c_efflux
          write(fates_log(),*) 'n_gain0: ',n_gain0,' allocated_n: ',allocated_n
          write(fates_log(),*) 'p_gain0: ',p_gain0,' allocated_p: ',allocated_p

          do i = 1,num_organs
             i_org = l2g_organ_list(i)
             i_var = prt_global%sp_organ_map(i_org,carbon12_element)
             write(fates_log(),*) i_org, this%variables(i_var)%val(1)-state_c0(i_org)
          end do
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    ! IF this was prescribed, then we dictate the uptake
    ! and pass that back as an output, otherwise
    ! we set the gains to what we started with so that
    ! it can be used again for mass balance checking and diagnostics
    
    if(n_uptake_mode.eq.prescribed_n_uptake) then
       n_gain = n_gain0-n_gain
    else
       n_gain = n_gain0
    end if
    if(p_uptake_mode.eq.prescribed_p_uptake) then
       p_gain = p_gain0-p_gain
    else
       p_gain = p_gain0
    end if



    ! If fine-roots are allocated above their
    ! target (perhaps with some buffer, but perhaps not)
    ! then 
    call this%TrimFineRoot()

    return
  end subroutine DailyPRTAllometricCNP

  
  function SafeLog(val) result(logval)

    ! The log functions used to transform storage ratios
    ! need not be large. Even a ratio of 10 is sending a strong signal to the
    ! root adaptation algorithm to change course pretty strongly. We set
    ! bounds of e3 here to prevent numerical overflows and underflows
    
    real(r8) :: val
    real(r8) :: logval
    real(r8), parameter :: safelog_min = 0.001_r8  !Don't pass anything smaller to a log
    real(r8), parameter :: safelog_max = 1000._r8

    logval = log(max(safelog_min,min(safelog_max,val)))
    
  end function SafeLog

  
  ! =====================================================================================
  
  subroutine CNPAdjustFRootTargets(this, target_c, target_dcdd)

    class(cnp_allom_prt_vartypes) :: this
    real(r8)                      :: target_c(:)
    real(r8)                      :: target_dcdd(:)

    real(r8), pointer :: l2fr           ! leaf to fineroot target biomass scaler
    integer           :: ipft           ! PFT index
    real(r8), pointer :: dbh
    real(r8)          :: canopy_trim
    integer  :: leaf_status
    integer, pointer :: limiter
    real(r8) :: elongf_fnrt
    real(r8) :: store_c_max, store_c_act
    real(r8) :: store_nut_max, store_nut_act
    real(r8) :: l2fr_delta
    real(r8) :: cn_ratio, cp_ratio ! ratio of relative C storage over relative N or P storage
    real(r8) :: dcxdt_ratio        ! log change (derivative) of the maximum of the N/C and P/C storage ratio
    real(r8) :: cx_logratio        ! log Maximum of the C/N and C/P storage ratio
    real(r8), pointer :: cx_int    ! Integration of the cx_logratio 
    real(r8), pointer :: cx0       ! The log of the cx ratio from previous time-step
    real(r8), pointer :: ema_dcxdt ! the EMA of the change in log storage ratio

    real(r8), parameter :: pid_drv_wgt = 1._r8/20._r8   ! n-day smoothing of the derivative
                                                        ! of the process function in the PID controller

    leaf_status = this%bc_in(acnp_bc_in_id_lstat)%ival
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    elongf_fnrt = this%bc_in(acnp_bc_in_id_effnrt)%rval
    l2fr        => this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim =  this%bc_in(acnp_bc_in_id_ctrim)%rval
    cx_int      => this%bc_inout(acnp_bc_inout_id_cx_int)%rval
    cx0         => this%bc_inout(acnp_bc_inout_id_cx0)%rval
    ema_dcxdt   => this%bc_inout(acnp_bc_inout_id_emadcxdt)%rval
    limiter     => this%bc_out(acnp_bc_out_id_limiter)%ival
    ! Abort if leaves are off
    if(leaf_status.eq.leaves_off) return


    ! Step 1: Determine the process function for the controller.  Generally, this is
    ! some indicator about the relative health of the plant in terms of carbon versus
    ! nutrient.  There are a few ways to cast this function, but right now we are using
    ! the relative amount of Carbon storage (actual/maximum) divided by the relative amount
    ! of nutrient (actual/maximum). We take the natural log of this ratio. And then we take
    ! maximum of the two quotients that use nitrogen and phosphorus.
    ! -----------------------------------------------------------------------------------

    store_c_max = target_c(store_organ)

    store_c_act = max(0.001_r8*store_c_max,this%GetState(store_organ, carbon12_element) + &
         this%bc_in(acnp_bc_in_id_netdc)%rval)

    if(n_uptake_mode.eq.prescribed_n_uptake)then
       cn_ratio = -1._r8
    else

       ! Calculate the relative nitrogen storage fraction,
       ! over the relative carbon storage fraction.

       store_nut_max = this%GetNutrientTarget(nitrogen_element,store_organ,stoich_growth_min)

       store_nut_act = max(0.001_r8*store_nut_max, &
            this%GetState(store_organ, nitrogen_element) + &
            this%bc_inout(acnp_bc_inout_id_netdn)%rval)

       cn_ratio = (store_c_act/store_c_max)/(store_nut_act/store_nut_max)

    end if

    if(p_uptake_mode.eq.prescribed_p_uptake)then
       cp_ratio = -1._r8
    else

       ! Calculate the relative phosphorus storage fraction,
       ! over the relative carbon storage fraction.

       store_nut_max = this%GetNutrientTarget(phosphorus_element,store_organ,stoich_growth_min) 

       store_nut_act = max(0.001_r8*store_nut_max, &
            this%GetState(store_organ, phosphorus_element) + &
            this%bc_inout(acnp_bc_inout_id_netdp)%rval)

       cp_ratio = (store_c_act/store_c_max)/(store_nut_act/store_nut_max)

    end if

    ! Use the limiting nutrient species
    if( (n_uptake_mode.eq.prescribed_n_uptake) .and. &
         (p_uptake_mode.eq.prescribed_p_uptake) )then
       cx_int = 0._r8
       ema_dcxdt = 0._r8
       cx0 = 0.0_r8
       return
    else

       if (n_uptake_mode.eq.prescribed_n_uptake) then
          cx_logratio = SafeLog(cp_ratio)
       elseif (p_uptake_mode.eq.prescribed_p_uptake) then
          cx_logratio = SafeLog(cn_ratio)
       else
          cx_logratio = SafeLog(max(cp_ratio,cn_ratio))
       end if

       ! If cx_logratio has just crossed zero, then
       ! reset the integrator. This will be true if
       ! the sign of the current ratio is different than
       ! the sign of the previous

       cx_int = cx_int + cx_logratio

       ! Reset the integrator if its sign changes
       if( abs(cx_logratio)>nearzero .and. abs(cx0)>nearzero) then
          if( abs(cx_logratio/abs(cx_logratio) - cx0/abs(cx0)) > nearzero ) then
             cx_int = cx_logratio
          end if
       end if

       dcxdt_ratio = cx_logratio-cx0

       ema_dcxdt = pid_drv_wgt*dcxdt_ratio + (1._r8-pid_drv_wgt)*ema_dcxdt

       cx0 = cx_logratio


    end if

    l2fr_delta = prt_params%pid_kp(ipft)*cx_logratio + &
         prt_params%pid_ki(ipft)*cx_int + &
         prt_params%pid_kd(ipft)*ema_dcxdt

    ! Apply the delta, also, avoid generating incredibly small l2fr's,
    ! super small l2frs will occur in plants that perpetually get almost
    ! now carbon gain, such as newly recruited plants in a dark understory

    l2fr = max(l2fr_min, l2fr + l2fr_delta)

    ! Find the updated target fineroot biomass
    call bfineroot(dbh,ipft,canopy_trim, l2fr, elongf_fnrt, target_c(fnrt_organ),target_dcdd(fnrt_organ))

    return
  end subroutine CNPAdjustFRootTargets

  ! =====================================================================================

  subroutine TrimFineRoot(this)
    
    ! The following section allows forceful turnover of fine-roots if a new L2FR is generated
    ! that is lower than the previous l2fr. The maintenance turnover (background) rate
    ! will automatically accomodate a lower l2fr, but if the change is large it will
    ! not keep pace.  Note 1: however, that the algorithm for calculating l2fr will prevent
    ! large drops in l2fr (unless that safegaurd is removed). Note 2: this section may also
    ! generate mass check errors in the main CNPAllocation routine, this is because the "val" is
    ! changing but the net_allocated is not reciprocating, which is expected. 
    
    ! Keep a buffer above the L2FR in the hopes that natural turnover will catch
    ! up.
    class(cnp_allom_prt_vartypes) :: this
       
    real(r8) :: fnrt_flux_c
    real(r8) :: turn_flux_c
    real(r8) :: store_flux_c
    real(r8) :: nc_fnrt
    real(r8) :: pc_fnrt
    real(r8) :: target_fnrt_c
    real(r8),parameter :: nday_buffer = 0._r8
    real(r8),parameter :: fnrt_opt_eff = 0._r8  ! If we want to transfer resources to storage
    
    if(.not.use_unrestricted_contraction)return
    
    associate( ipft         => this%bc_in(acnp_bc_in_id_pft)%ival,  &
         l2fr         => this%bc_inout(acnp_bc_inout_id_l2fr)%rval, &
         dbh          => this%bc_inout(acnp_bc_inout_id_dbh)%rval,  &
         elongf_fnrt  => this%bc_in(acnp_bc_in_id_effnrt)%rval,     &
         canopy_trim  => this%bc_in(acnp_bc_in_id_ctrim)%rval)

      ! Find the updated target fineroot biomass
      call bfineroot(dbh,ipft,canopy_trim, l2fr, elongf_fnrt, target_fnrt_c)

      fnrt_flux_c = max(0._r8,this%variables(fnrt_c_id)%val(1)*(1._r8-nday_buffer*(years_per_day / prt_params%root_long(ipft))) - target_fnrt_c )

      if(fnrt_flux_c>nearzero) then

         turn_flux_c = (1._r8 - fnrt_opt_eff)*fnrt_flux_c
         store_flux_c = fnrt_opt_eff*fnrt_flux_c

         nc_fnrt = this%variables(fnrt_n_id)%val(1)/this%variables(fnrt_c_id)%val(1)
         pc_fnrt = this%variables(fnrt_p_id)%val(1)/this%variables(fnrt_c_id)%val(1)

         this%variables(fnrt_c_id)%val(1)        = this%variables(fnrt_c_id)%val(1) - fnrt_flux_c
         this%variables(fnrt_c_id)%turnover(1)   = this%variables(fnrt_c_id)%turnover(1) + turn_flux_c
         this%variables(fnrt_c_id)%net_alloc(1)  = this%variables(fnrt_c_id)%net_alloc(1) - store_flux_c
         this%variables(store_c_id)%val(1)       = this%variables(store_c_id)%val(1) + store_flux_c
         this%variables(store_c_id)%net_alloc(1) = this%variables(store_c_id)%net_alloc(1) + store_flux_c

         this%variables(fnrt_n_id)%val(1)        = this%variables(fnrt_n_id)%val(1) - fnrt_flux_c * nc_fnrt
         this%variables(fnrt_n_id)%turnover(1)   = this%variables(fnrt_n_id)%turnover(1) + turn_flux_c * nc_fnrt
         this%variables(fnrt_n_id)%net_alloc(1)  = this%variables(fnrt_n_id)%net_alloc(1) - store_flux_c * nc_fnrt
         this%variables(store_n_id)%val(1)       = this%variables(store_n_id)%val(1) + store_flux_c * nc_fnrt
         this%variables(store_n_id)%net_alloc(1) = this%variables(store_n_id)%net_alloc(1) + store_flux_c * nc_fnrt

         this%variables(fnrt_p_id)%val(1)        = this%variables(fnrt_p_id)%val(1) - fnrt_flux_c * pc_fnrt
         this%variables(fnrt_p_id)%turnover(1)   = this%variables(fnrt_p_id)%turnover(1) + turn_flux_c * pc_fnrt
         this%variables(fnrt_p_id)%net_alloc(1)  = this%variables(fnrt_p_id)%net_alloc(1) - store_flux_c * pc_fnrt
         this%variables(store_p_id)%val(1)       = this%variables(store_p_id)%val(1) + store_flux_c * pc_fnrt
         this%variables(store_p_id)%net_alloc(1) = this%variables(store_p_id)%net_alloc(1) + store_flux_c * pc_fnrt

      end if
    end associate
    return
  end subroutine TrimFineRoot
  
  ! =====================================================================================
    
  subroutine CNPPrioritizedReplacement(this,c_gain, n_gain, p_gain, target_c)
    
      
    ! -----------------------------------------------------------------------------------
    ! Alternative allocation hypothesis for the prioritized replacement phase.
    ! This is more similar to the current (04/2020) carbon only hypothesis.
    ! -----------------------------------------------------------------------------------
    
    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain
    real(r8), intent(inout) :: n_gain
    real(r8), intent(inout) :: p_gain
    real(r8), intent(in)    :: target_c(:)      ! Indexed by global organ (from PRTGenericMod)

    integer  :: n_curpri_org

    integer, dimension(num_organs) :: curpri_org   ! organ ID's of the current priority level
    real(r8), dimension(num_organs) :: deficit_c ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_n ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_p ! Deficit to get to target from current    [kg]
    
    integer  :: i, ii, i_org             ! Loop indices (mostly for organs)
    integer  :: i_var                    ! variable index
    integer  :: i_pri                    ! loop index for priority
    integer  :: ipft                     ! Plant functional type index of this plant
    integer  :: leaf_status              ! Is this plant in a leaf on or off status?
    real(r8) :: canopy_trim              ! trim factor for maximum leaf biomass
    real(r8) :: target_n                 ! Target mass of N for a given organ [kg]
    real(r8) :: target_p                 ! Target mass of P for a given organ [kg]
    real(r8) :: elongf_leaf              ! Leaf elongation factor
    real(r8) :: elongf_fnrt              ! Fine-root "elongation factor"
    real(r8) :: elongf_stem              ! Stem "elongation factor"
    integer  :: priority_code            ! Index for priority level of each organ
    real(r8) :: sum_c_demand             ! Carbon demanded to bring tissues up to allometry (kg)
    real(r8) :: store_below_target       ! The amount of storage that is less than the target (kg)
    real(r8) :: store_target_fraction    ! The fraction of actual storage carbon over the target (kg)
    real(r8) :: store_demand             ! Based on the target fraction, an exponential function defining
                                         ! how much carbon we should try to put back into storage
    real(r8) :: store_c_flux             ! The amount of C we draw from gains to give back to storage (kg)
    real(r8) :: sum_c_flux               ! The flux to bring tissues up to allometry (kg)
    real(r8) :: c_flux                   ! carbon flux into an arbitrary pool (kg)
    integer :: n_max_priority    ! Maximum possible number of priority levels is
                                 ! the total number organs plus 1, which allows
                                 ! each organ to have its own level, and ignore
                                 ! the specialized priority 1
    
    
    leaf_status     = this%bc_in(acnp_bc_in_id_lstat)%ival
    elongf_leaf     = this%bc_in(acnp_bc_in_id_efleaf)%rval
    elongf_fnrt     = this%bc_in(acnp_bc_in_id_effnrt)%rval
    elongf_stem     = this%bc_in(acnp_bc_in_id_efstem)%rval
    ipft            = this%bc_in(acnp_bc_in_id_pft)%ival
    canopy_trim     = this%bc_in(acnp_bc_in_id_ctrim)%rval

    
    n_max_priority = maxval(prt_params%organ_param_id(:))
    if(n_max_priority>10 .or. n_max_priority<0)then
       write(fates_log(),*) 'was unable to interpret prt_params%organ_param_id'
       write(fates_log(),*) 'for cnp allocation, there should be non-zero values <10'
       write(fates_log(),*) 'your values: ',prt_params%organ_param_id(:)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! -----------------------------------------------------------------------------------
    ! Notes on indexes:
    !
    ! i_org: this is the index that matches the global organ indices found in PRTGenericMod
    !
    ! prt_params%alloc_priority is a parameter array that only holds a subset of the
    ! organs.  For instance it does not include reproductive and storage, because those
    ! organs are special.  We can find out the global organ index (ie i_org) from the
    ! parameter array prt_params%organ_id
    !
    ! i and ii are just local indices used to iterate through sub-groups of organs.
    ! for instance, i will generally iterate through the curpri_org, which is an
    ! array of the global organ indices (i_org) for the current subset of organs that
    ! should be allocated at the current priority level
    !
    ! -----------------------------------------------------------------------------------

    
    ! -----------------------------------------------------------------------------------
    ! Preferential transfer of available carbon and nutrients into the highest
    ! priority pools, and maintenance respiration. We will loop through the available
    ! pools, and identify if that pool is part of the highest transfer priority.
    ! If it is, then we track the variable ids associated with that pool for each CNP
    ! species.  It "should" work fine if there are NO priority=1 pools...
    ! -----------------------------------------------------------------------------------

    curpri_org(:) = fates_unset_int    ! reset "current-priority" organ ids
    i = 0
    do ii = 1, size(prt_params%organ_id,1) 

       ! universal organ index from PRTGenericMod
       i_org = prt_params%organ_id(ii)

       ! Don't allow allocation to leaves if they are in an "off" status.
       ! Also, dont allocate to replace turnover if this is not evergreen
       ! (this prevents accidental re-flushing on the day they drop)
       if( ( any(leaf_status == [leaves_off,leaves_shedding]) .or. &
             (prt_params%evergreen(ipft) /= itrue) ) &
            .and. (i_org == leaf_organ)) cycle

       ! The priority code associated with this organ
       priority_code = int(prt_params%alloc_priority(ipft, ii))

       ! 1 is the highest priority code possible
       if( priority_code == 1 ) then
          i = i + 1
          curpri_org(i) = i_org
          deficit_c(i)  = max(0._r8,this%GetDeficit(carbon12_element,i_org,target_c(i_org)))
       end if

    end do

      
    ! Number of pools in the current priority level
    n_curpri_org = i

    ! -----------------------------------------------------------------------------------
    ! The high-priority pools, and their associated variable
    ! indices have been identified.
    !
    ! Let us now calculate the fluxes into these priority pools
    ! The first step is to replace just their maintenance turnover
    ! -----------------------------------------------------------------------------------

    sum_c_demand = 0._r8
    do i = 1,n_curpri_org
       i_org = curpri_org(i)
       i_var       = prt_global%sp_organ_map(i_org,carbon12_element)
       sum_c_demand = sum_c_demand + prt_params%leaf_stor_priority(ipft) * &
            sum(this%variables(i_var)%turnover(:))
    end do

    sum_c_flux = max(0._r8,min(sum_c_demand,this%variables(store_c_id)%val(1)+c_gain))
    
    if (sum_c_flux> nearzero ) then
       
       ! We pay this even if we don't have the carbon
       ! Just don't pay so much carbon that storage+carbon_balance can't pay for it
       
       do i = 1,n_curpri_org

          i_org = curpri_org(i)
          i_var = prt_global%sp_organ_map(i_org,carbon12_element)

          c_flux = sum_c_flux*(prt_params%leaf_stor_priority(ipft) * &
               sum(this%variables(i_var)%turnover(:))/sum_c_demand)
          
          ! Add carbon to the pool
          this%variables(i_var)%val(1) = this%variables(i_var)%val(1) + c_flux
          
          ! Remove from daily  carbon gain
          c_gain = c_gain - c_flux

       end do
    end if

    ! Determine nutrient demand and make tansfers
    do i = 1, n_curpri_org
       
       i_org = curpri_org(i)
       
       ! Update the nitrogen deficits
       ! Note that the nitrogen target is tied to the stoichiometry of the growing pool only (pos = 1)

       target_n = this%GetNutrientTarget(nitrogen_element,i_org,stoich_growth_min)
       deficit_n(i) = max(0.0_r8, target_n - this%GetState(i_org, nitrogen_element,1))
       
       ! Update the phosphorus deficits (which are based off of carbon actual..)
       ! Note that the phsophorus target is tied to the stoichiometry of thegrowing pool only (also)
       target_p = this%GetNutrientTarget(phosphorus_element,i_org,stoich_growth_min)
       deficit_p(i) = max(0.0_r8, target_p - this%GetState(i_org, phosphorus_element,1))

    end do
    
    ! Allocate nutrients at this priority level
    ! Nitrogen
    call ProportionalNutrAllocation(this,deficit_n(1:n_curpri_org), &
         n_gain, nitrogen_element, curpri_org(1:n_curpri_org))
    
    ! Phosphorus
    call ProportionalNutrAllocation(this,deficit_p(1:n_curpri_org), &
         p_gain, phosphorus_element, curpri_org(1:n_curpri_org))
    
    ! -----------------------------------------------------------------------------------
    ! IV. if carbon balance is negative, re-coup the losses from storage
    !       if it is positive, give some love to storage carbon
    ! -----------------------------------------------------------------------------------
    
    if( c_gain < 0.0_r8 ) then

       ! Storage will have to pay for any negative gains
       store_c_flux           = -c_gain
       c_gain                 = c_gain                + store_c_flux

       this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) - store_c_flux
       
    else

       ! This is just a cap, don't fill up more than is needed (shouldn't even apply)
       store_below_target     = max(target_c(store_organ) - this%variables(store_c_id)%val(1),0._r8)
       
       ! This is the desired need for carbon
       store_target_fraction  = max(this%variables(store_c_id)%val(1)/target_c(store_organ),0._r8)
       store_demand           = max(c_gain*(exp(-1.*store_target_fraction**4._r8) - exp( -1.0_r8 )),0._r8)

       ! The flux is the (positive) minimum of all three
       store_c_flux           = min(store_below_target,store_demand)
       
       c_gain                 = c_gain  - store_c_flux

       this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) + store_c_flux
   
   end if
   
    
    ! -----------------------------------------------------------------------------------
    !  If carbon is still available, allocate to remaining high
    !        carbon balance is guaranteed to be >=0 beyond this point
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Bring all pools, in priority order, up to allometric targets if possible
    ! Repeat priority order 1 as well.
    ! -----------------------------------------------------------------------------------
    
    priority_loop: do i_pri = 1, n_max_priority
       
       curpri_org(:) = fates_unset_int    ! "current-priority" organ indices

       i = 0

       ! Storage has a special hard-coded priority level of 2
       if( i_pri == 2 ) then
          curpri_org(1) = store_organ
          i=1
       end if
          
       ! Loop over all organs in the CNP routine, which
       do ii = 1, size(prt_params%organ_id,1) 

          ! universal organ index from PRTGenericMod
          i_org = prt_params%organ_id(ii)
          
          ! The priority code associated with this organ
          
          priority_code = int(prt_params%alloc_priority(ipft,ii))
          
          ! Don't allow allocation to leaves if they are in an "off" status.
          ! (this prevents accidental re-flushing on the day they drop)
          if( any(leaf_status == [leaves_off,leaves_shedding]) .and. &
              (i_org == leaf_organ) ) cycle

          if( priority_code == i_pri ) then
             i = i + 1
             curpri_org(i) = i_org
          end if
       end do

       n_curpri_org  = i

       do i = 1,n_curpri_org
          i_org = curpri_org(i)
          deficit_c(i) = max(0._r8,this%GetDeficit(carbon12_element,i_org,target_c(i_org)))
       end do
          
       ! Bring carbon up to target first, this order is required
       ! because we need to know the resulting carbon concentrations
       ! before  we set the allometric targets for the nutrients

       sum_c_demand = 0._r8
       do i=1,n_curpri_org
          i_org = curpri_org(i)
          sum_c_demand = sum_c_demand + deficit_c(i)
       end do
       
       sum_c_flux = min(c_gain, sum_c_demand)
       
       ! Transfer carbon into pools if there is any
       if (sum_c_flux>nearzero) then
          do i = 1, n_curpri_org
             
             i_org = curpri_org(i)

             c_flux  =  sum_c_flux*deficit_c(i)/sum_c_demand
             
             ! Update the carbon pool
             i_var = prt_global%sp_organ_map(i_org,carbon12_element)
             this%variables(i_var)%val(1) = this%variables(i_var)%val(1) + c_flux
             
             ! Update carbon pools deficit
             deficit_c(i) = max(0._r8,deficit_c(i) - c_flux)
             
             ! Reduce the carbon gain
             c_gain      = c_gain - c_flux
             
          end do
       end if

       ! Determine nutrient demand and make tansfers
       do i = 1, n_curpri_org
          
          i_org = curpri_org(i)

          ! Update the nitrogen deficits
          ! Note that the nitrogen target is tied to the stoichiometry of thegrowing pool only
          target_n = this%GetNutrientTarget(nitrogen_element,i_org,stoich_growth_min)
          deficit_n(i) = max(0.0_r8, target_n - this%GetState(i_org, nitrogen_element,1) )
             
          ! Update the phosphorus deficits (which are based off of carbon actual..)
          ! Note that the phsophorus target is tied to the stoichiometry of thegrowing pool only (also)
          target_p = this%GetNutrientTarget(phosphorus_element,i_org,stoich_growth_min)
          deficit_p(i) = max(0.0_r8, target_p - this%GetState(i_org, phosphorus_element,1) )

       end do
          
       ! Allocate nutrients at this priority level Nitrogen
       call ProportionalNutrAllocation(this,deficit_n(1:n_curpri_org), &
            n_gain, nitrogen_element, curpri_org(1:n_curpri_org))
       
       ! Phosphorus
       call ProportionalNutrAllocation(this,deficit_p(1:n_curpri_org), &
            p_gain, phosphorus_element, curpri_org(1:n_curpri_org))
       

    end do priority_loop


    
    
    return
  end subroutine CNPPrioritizedReplacement

  
  ! =====================================================================================
    
  subroutine CNPStatureGrowth(this,c_gain, n_gain, p_gain, &
                              target_c, target_dcdd)
    
    
    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain      ! Total daily C gain that remains to be used
    real(r8), intent(inout) :: n_gain      ! Total N available for allocation
                                           ! (new uptake + storage)
    real(r8), intent(inout) :: p_gain      ! Total P available for allocation
                                           ! (new uptake + storage)
    real(r8), intent(in) :: target_c(:)    ! target carbon mass for each organ (before growth)
    real(r8), intent(in) :: target_dcdd(:) ! target carbon mass derivative (wrt dbh) before growth)


    real(r8), pointer :: dbh
    integer           :: ipft
    integer, pointer  :: limiter                 ! Integer flagging which (C,N,P) is limiting
    real(r8)          :: canopy_trim             ! fraction of crown trimmed
    integer           :: crown_damage            ! Damage status level
    real(r8)          :: elongf_leaf             ! Elongation factor (leaves)
    real(r8)          :: elongf_fnrt             ! Elongation factor (fine roots)
    real(r8)          :: elongf_stem             ! Elongation factor (woods)
    real(r8)          :: leaf_status             ! leaves on or off?
    real(r8)          :: l2fr                    ! leaf to fineroot allometry multiplier
    integer  :: i, ii                            ! organ index loops (masked and unmasked)
    integer  :: i_org                            ! global organ index
    real(r8) :: total_dcostdd                    ! Total carbon transferred to all pools for unit growth
    logical  :: step_pass                        ! flag stating if the integration sub-steps passed checks
    real(r8) :: totalC                           ! total carbon sent to integrator (kg)
    real(r8) :: deltaC                           ! trial value for substep change in carbon (kg)
    real(r8) :: cdeficit                         ! carbon deficit from target
    integer  :: ierr                             ! error flag for allometric growth step
    integer  :: nsteps                           ! number of sub-steps
    real(r8) :: avg_nc,avg_pc                    ! Estimated average N/C and P/C ratios of
                                                 ! allocated carbon during stature growth
    real(r8) :: repro_c_frac                     ! Fraction of C allocated to reproduction
                                                 ! at current stature (dbh) [/]
    real(r8) :: sum_c_flux                       ! Sum of the carbon allocated, as reported
                                                 ! by the ODE solver. [kg]
    real(r8) :: c_flux_adj                       ! Adjustment to total carbon flux during stature growth
                                                 ! intended to correct integration error (kg/kg)
    real(r8) :: c_flux                           ! Carbon flux from the gain pool to an organ (kgC)
    real(r8) :: n_flux,p_flux
    real(r8) :: c_gstature                       ! Carbon reserved for stature growth (kg)
    real(r8) :: target_n                         ! Target mass of N for a given organ [kg]
    real(r8) :: target_p                         ! Target mass of P for a given organ [kg]
    real(r8) :: sum_n_demand                     ! Total N deficit to overcome after C stature growth [kg]
    real(r8) :: sum_p_demand                     ! Total P deficit to overcome after C stature growth [kg]
    real(r8), dimension(num_organs) :: deficit_n ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_p ! Deficit to get to target from current    [kg]
    integer,dimension(num_organs) :: mask_organs ! This works with "state_mask", the list
                                                 ! of organs (local ids) in the mask
    integer,dimension(num_organs) :: mask_gorgans ! List of organ global indices in the mask
    integer                       :: n_mask_organs
    
    ! Integrator error checking
    integer  :: i_var
    integer  :: nbins
    real(r8) :: dbh_tp1
    real(r8) :: leafc_tp1
    real(r8) :: fnrtc_tp1
    real(r8) :: sapwc_tp1
    real(r8) :: storec_tp1
    real(r8) :: structc_tp1
    real(r8) :: leaf_c_target_tp1
    real(r8) :: fnrt_c_target_tp1
    real(r8) :: sapw_c_target_tp1
    real(r8) :: agw_c_target_tp1
    real(r8) :: bgw_c_target_tp1
    real(r8) :: struct_c_target_tp1
    real(r8) :: store_c_target_tp1
    real(r8) :: sapw_area
    
    ! Integegrator variables
    ! These are not global because we want a unique instance for each time the routine is called
    ! ----------------------------------------------------------------------------------------
    
    real(r8),dimension(num_intgr_vars) :: state_array      ! Vector of carbon pools passed to integrator
    real(r8),dimension(num_intgr_vars) :: state_array_out  ! Vector of carbon pools passed back from integrator
    logical,dimension(num_intgr_vars)  :: state_mask       ! Mask of active pools during integration
    integer , parameter :: max_substeps = 300            ! Maximum allowable iterations
    real(r8), parameter :: max_trunc_error = 1.0_r8        ! Maximum allowable truncation error
    integer,  parameter :: ODESolve = 2                    ! 1=RKF45,  2=Euler

    real(r8)            :: intgr_params(num_intgr_parm)

    
    real(r8) :: neq_cgain, peq_cgain  ! N and P equivalent c_gain spent on growth
    real(r8) :: cnp_gain              ! used as a check to see efficiency of limited growth
    
    

    leaf_status = this%bc_in(acnp_bc_in_id_lstat)%ival
    elongf_leaf = this%bc_in(acnp_bc_in_id_efleaf)%rval
    elongf_fnrt = this%bc_in(acnp_bc_in_id_effnrt)%rval
    elongf_stem = this%bc_in(acnp_bc_in_id_efstem)%rval
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    crown_damage = this%bc_in(acnp_bc_in_id_cdamage)%ival
    limiter     => this%bc_out(acnp_bc_out_id_limiter)%ival
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    l2fr        = this%bc_inout(acnp_bc_inout_id_l2fr)%rval ! This variable is not updated in this
                                                            ! routine, and is therefore not a pointer
    
    if( c_gain <= calloc_abs_error ) then
       limiter = c_limited
       if((n_gain <= 0.1_r8*calloc_abs_error) .or. &
          (p_gain <= 0.02_r8*calloc_abs_error)) limiter = cnp_limited
    else
       if(n_gain <= 0.1_r8*calloc_abs_error) limiter = n_limited
       if(p_gain <= 0.02_r8*calloc_abs_error) limiter = p_limited
    end if

    limiter = 0
    
    ! If any of these resources is essentially tapped out,
    ! then there is no point in performing growth
    ! It also seems impossible that we would be in a leaf-off status
    ! and have enough carbon to grow stature, but its possible that
    ! a plant had a productive last day before the phenology scheme
    ! signaled a drop. If this is the case, we can't grow stature
    ! cause that would force the leaves back on, so just leave.

    if( c_gain <= calloc_abs_error .or. &
        n_gain <= 0.1_r8*calloc_abs_error .or. &
        p_gain <= 0.02_r8*calloc_abs_error .or. &
        any(leaf_status == [leaves_off,leaves_shedding]) ) then
       return
    end if
       
    intgr_params(:)                  = fates_unset_r8
    intgr_params(intgr_parm_ctrim)   = this%bc_in(acnp_bc_in_id_ctrim)%rval
    intgr_params(intgr_parm_pft)     = real(this%bc_in(acnp_bc_in_id_pft)%ival,r8)
    intgr_params(intgr_parm_l2fr)    = this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    intgr_params(intgr_parm_cdamage) = real(this%bc_in(acnp_bc_in_id_cdamage)%ival,r8)
    intgr_params(intgr_parm_efleaf)  = this%bc_in(acnp_bc_in_id_efleaf)%rval
    intgr_params(intgr_parm_effnrt)  = this%bc_in(acnp_bc_in_id_effnrt)%rval
    intgr_params(intgr_parm_efstem)  = this%bc_in(acnp_bc_in_id_efstem)%rval
    state_mask(:) = .false.
    mask_organs(:) = fates_unset_int
    mask_gorgans(:) = fates_unset_int
    
    ! Go through and flag the integrating variables as either pools that
    ! are growing in this iteration, or not.   At this point, if carbon for growth
    ! remains, it means that all pools are up to, or above the target.  If
    ! it is above, that is because of numerical integration errors, or fusion.
    ! In that case, we flag that pool to not be included in stature growth. It will
    ! catch up with the other pools in the next stature growth steps.

    ii = 0
    do i = 1, num_organs

       i_org = l2g_organ_list(i)
       
       cdeficit = this%GetDeficit(carbon12_element,i_org,target_c(i_org))
       
       if ( cdeficit > calloc_abs_error ) then
          ! In this case, we somehow still have carbon to play with,
          ! yet one of the pools is below its current target
          ! gracefully fail
          write(fates_log(),*) 'A carbon pool has reached the stature growth step'
          write(fates_log(),*) 'yet its deficit is too large to integrate '
          write(fates_log(),*) 'organ: ',i_org
          write(fates_log(),*) 'carbon gain: ',c_gain
          write(fates_log(),*) 'leaves status:', leaf_status
          write(fates_log(),*) cdeficit, target_c(i_org)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       elseif( (-cdeficit) > calloc_abs_error ) then
          ! In this case, we are above our target (ie negative deficit (fusion?))
          ! and if so, this pool does not have to grow and will
          ! catch up in the next step
          state_mask(i)            = .false.    ! flag index of state variable
       else
          ! In this case, the pool is close enough to the target
          ! to be flagged for growth
          state_mask(i)            = .true.    ! flag index of state variable

          ! Reproduction is a special case, don't add it to the
          ! list of organs... yet
          if (i_org.ne.repro_organ) then
             ii=ii+1
             mask_organs(ii) = i
             mask_gorgans(ii) = i_org
          end if
          
       end if
    end do

    n_mask_organs = ii

    if(debug) then
       if ( .not.any(state_mask(1:num_organs)) ) then
          write(fates_log(),*) 'No organs seem to have carbon masses that are'
          write(fates_log(),*) 'on allometry. Apparently, all are above allometry'
          write(fates_log(),*) 'which should be impossible. We allow for all but'
          write(fates_log(),*) 'one pool to be above allometry. Structure in woody'
          write(fates_log(),*) 'plants, and roots in grasses are not allowed above target.'
          write(fates_log(),*) 'pft: ',ipft
          write(fates_log(),*) 'dbh: ',dbh
          write(fates_log(),*) 'c targets: ',target_c(:)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    ! fraction of carbon going towards reproduction. reproductive carbon is
    ! just different from the other pools. It is not based on proportionality,
    ! so its mask is set differently.  We (inefficiently) just included
    ! reproduction in the previous loop, but oh well, we over-write now.

    ! If the TRS is switched off, or if the plant is a shrub or grass
    ! then we use FATES's default reproductive allocation.
    ! We designate a plant a shrub or grass if its dbh at maximum height
    ! is less than 15 cm

    if ( regeneration_model == default_regeneration .or. &
         prt_params%allom_dbh_maxheight(ipft) < min_max_dbh_for_trees ) then

       if (dbh <= prt_params%dbh_repro_threshold(ipft)) then
          repro_c_frac = prt_params%seed_alloc(ipft)
       else
          repro_c_frac = prt_params%seed_alloc(ipft) + prt_params%seed_alloc_mature(ipft)
       end if
    
    ! If the TRS is switched on (with or w/o seedling dynamics) then reproductive allocation is
    ! a pft-specific function of dbh. This allows for the representation of different
    ! reproductive schedules (Wenk and Falster, 2015)
    else if ( any(regeneration_model == [TRS_regeneration, TRS_no_seedling_dyn]) .and. &
                  prt_params%allom_dbh_maxheight(ipft) > min_max_dbh_for_trees ) then

       repro_c_frac = prt_params%seed_alloc(ipft) * &
       (exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm) / &
       (1 + exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm)))

    else
       
       write(fates_log(),*) 'unknown seed allocation and regeneration model, exiting'
       write(fates_log(),*) 'regeneration_model: ',regeneration_model
       call endrun(msg=errMsg(sourcefile, __LINE__))
       
    end if ! regeneration switch 


    if(repro_c_frac>nearzero)then
       state_mask(repro_id)            = .true.
       n_mask_organs = n_mask_organs + 1
       mask_organs(n_mask_organs)  = repro_id
       mask_gorgans(n_mask_organs) = repro_organ
    else
       state_mask(repro_id)            = .false.
    end if

    ! Calculate the total CARBON allocation rate per diameter increment
    ! Include the growth respiration costs "total_dcostdd" in that estimate
    ! ALso include the non-respiring cost, which is needed to project reproductive
    ! costs.
    ! --------------------------------------------------------------------------------

    ! First objective is to find the extrapolated proportions of carbon going to
    ! each pool. This has nothing to do with carbon conservation, it is just used
    ! to make a rough prediction of how much nutrient is needed to match carbon
    
    total_dcostdd = 0._r8

    do i = 1, n_mask_organs
       i_org = mask_gorgans(ii)
       total_dcostdd = total_dcostdd + target_dcdd(i_org)
    end do
    
    ! We can either proceed with stature growth by using all of the carbon
    ! available, or we can try to estimate the limitations of N and P
    ! and thereby reduce the amount of C we are willing to use to try
    ! and match what is available in n_gain and p_gain.  Note that the
    ! c-only option does allow limitations to eventually occur, because
    ! it is assumed that in the dynamics call that follows this one, there may
    ! or not be enough N or P to reach the stature growth step at all.

    if(grow_lim_type == grow_lim_conly) then
       c_gstature = c_gain
       limiter = 0
    elseif (grow_lim_type == grow_lim_estNP) then

       call EstimateGrowthNC(this,target_c,target_dcdd,state_mask,avg_nc,avg_pc)

       neq_cgain = n_gain/avg_nc
       peq_cgain = p_gain/avg_pc
       
       if(c_gain<neq_cgain) then

          if(c_gain < peq_cgain) then
             limiter = c_limited
             c_gstature = c_gain
             cnp_gain = c_gain
          else
             limiter = p_limited
             c_gstature = peq_cgain
             cnp_gain = p_gain
          end if
          
       else
          if(neq_cgain < peq_cgain) then
             limiter = n_limited
             c_gstature = neq_cgain
             cnp_gain   = n_gain
          else
             limiter = p_limited
             c_gstature = peq_cgain
             cnp_gain = p_gain
          end if
          
       end if
       
    end if

    
   if_stature_growth: if(c_gstature > nearzero) then

       ! Initialize the adaptive integrator arrays and flags
       ! -----------------------------------------------------------------------------------
       
      if(ODESolve == 2) then
         this%ode_opt_step = c_gstature
      end if
      
      ! If this flag is set to 0, then
      ! we have a successful integration
      ierr             = 1
      nsteps           = 0
      totalC           = c_gstature
      
      ! Fill the state array with element masses for each organ
      do i = 1, num_organs
         i_org = l2g_organ_list(i)
         i_var = prt_global%sp_organ_map(i_org,carbon12_element)
         state_array(i) = this%variables(i_var)%val(1)
      end do
      
      state_mask(dbh_id)       = .true.
      state_array(dbh_id)      = dbh

      do_solve_check: do while( ierr .ne. 0 )
         
         deltaC = min(totalC,this%ode_opt_step)
         if(ODESolve == 1) then
            
            call RKF45(AllomCNPGrowthDeriv,state_array,state_mask,deltaC,totalC, &
                 max_trunc_error,intgr_params,state_array_out,this%ode_opt_step,step_pass)
            
         elseif(ODESolve == 2) then
            
            call Euler(AllomCNPGrowthDeriv,state_array,state_mask, &
                 deltaC,totalC,intgr_params,state_array_out)

            ! Here we check to see if the solution is reasonably
            ! close to allometry, we also have to add up all leaf bins
            ! for this check.
            
            leafc_tp1 = state_array_out(leaf_id)
            i_var     = prt_global%sp_organ_map(leaf_organ,carbon12_element)
            nbins     = prt_global%state_descriptor(i_var)%num_pos
            do i = 2,nbins
               leafc_tp1 = leafc_tp1 + this%variables(i_var)%val(i)
            end do

            call CheckIntegratedAllometries(state_array_out(dbh_id),ipft,crown_damage,canopy_trim,  &
                 elongf_leaf, elongf_fnrt, elongf_stem, l2fr, &
                 leafc_tp1, state_array_out(fnrt_id), state_array_out(sapw_id), &
                 state_array_out(store_id), state_array_out(struct_id), &
                 state_mask(leaf_id), state_mask(fnrt_id), state_mask(sapw_id), &
                 state_mask(store_id),state_mask(struct_id),  max_trunc_error, step_pass)

            if(step_pass)  then
               this%ode_opt_step = deltaC
            else
               this%ode_opt_step = 0.5_r8*deltaC
            end if
            
         else
            write(fates_log(),*) 'An integrator was chosen that DNE'
            write(fates_log(),*) 'ODESolve = ',ODESolve
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         nsteps = nsteps + 1
         
         if(step_pass)  then
            totalC            = totalC - deltaC
            state_array(:)    = state_array_out(:)
         end if
          
         ! TotalC should eventually be whittled down to near-zero
         ! --------------------------------------------------------------------------------
         if_completed_solve:  if( (totalC < calloc_abs_error) )then
            
            ierr           = 0
            
            ! Sum up the total flux predicted by the integrator,
            ! which SHOULD be c_gstature, except
            ! for integration errors.  To make carbon
            ! perfectly preserved, we calculate this bias
            ! and make a linear (proportional) correction to all pools.
            
            sum_c_flux = 0.0_r8
            do ii = 1, n_mask_organs
               i     = mask_organs(ii)
               i_org = mask_gorgans(ii)
               i_var = prt_global%sp_organ_map(i_org,carbon12_element)
               sum_c_flux = sum_c_flux + (state_array(i) - this%variables(i_var)%val(1))
            end do
            
            ! This is a correction factor that forces
            ! mass conservation
            c_flux_adj = c_gstature/sum_c_flux
            
            do ii = 1, n_mask_organs
               
               i     = mask_organs(ii)
               i_org = mask_gorgans(ii)
               i_var = prt_global%sp_organ_map(i_org,carbon12_element)
               
               ! Calculate adjusted flux
               c_flux   = (state_array(i) - this%variables(i_var)%val(1))*c_flux_adj
               
               ! update the carbon pool (in all pools flux goes into the first pool)
               this%variables(i_var)%val(1) = this%variables(i_var)%val(1) + c_flux
               
               ! Remove carbon from the daily gain
               c_gain = c_gain - c_flux
               
            end do
            
            ! Update dbh
            dbh = state_array(dbh_id)       
            
         else

            if_step_exceedance: if (nsteps > max_substeps ) then
               
               write(fates_log(),*) 'CNP Plant Growth Integrator could not find'
               write(fates_log(),*) 'a solution in less than ',max_substeps,' tries'
               write(fates_log(),*) 'Aborting'
               write(fates_log(),*) 'mask: ',state_mask
               write(fates_log(),*) 'smallest deltaC',this%ode_opt_step
               write(fates_log(),*) 'totalC',totalC,c_gain,neq_cgain,peq_cgain
               write(fates_log(),*) 'pft: ',ipft
               write(fates_log(),*) 'trim: ',canopy_trim
               write(fates_log(),*) 'l2fr: ',l2fr
               write(fates_log(),*) 'dbh: ',dbh
               write(fates_log(),*) 'elongf_leaf: ',elongf_leaf
               write(fates_log(),*) 'elongf_fnrt: ',elongf_fnrt
               write(fates_log(),*) 'elongf_stem: ',elongf_stem
               write(fates_log(),*) 'dCleaf_dd: ',target_dcdd(leaf_organ)
               write(fates_log(),*) 'dCfnrt_dd: ',target_dcdd(fnrt_organ)
               write(fates_log(),*) 'dCstore_dd: ',target_dcdd(store_organ)
               write(fates_log(),*) 'dCsapw_dd: ',target_dcdd(sapw_organ)
               write(fates_log(),*) 'dCstruct_dd: ',target_dcdd(struct_organ)
               write(fates_log(),*) 'repro c frac: ',repro_c_frac
               dbh_tp1     = state_array_out(dbh_id)
               leafc_tp1   = state_array_out(leaf_id)
               fnrtc_tp1   = state_array_out(fnrt_id)
               sapwc_tp1   = state_array_out(sapw_id)
               storec_tp1  = state_array_out(store_id)
               structc_tp1 = state_array_out(struct_id)
               
               call bleaf(dbh_tp1,ipft,crown_damage,canopy_trim, elongf_leaf, leaf_c_target_tp1)
               call bfineroot(dbh_tp1,ipft,canopy_trim,l2fr, elongf_fnrt, fnrt_c_target_tp1)
               call bsap_allom(dbh_tp1,ipft,crown_damage,canopy_trim, elongf_stem, sapw_area,sapw_c_target_tp1)
               call bagw_allom(dbh_tp1,ipft,crown_damage, elongf_stem, agw_c_target_tp1)
               call bbgw_allom(dbh_tp1,ipft, elongf_stem, bgw_c_target_tp1)
               call bdead_allom(agw_c_target_tp1,bgw_c_target_tp1, sapw_c_target_tp1, ipft, struct_c_target_tp1)
               call bstore_allom(dbh_tp1,ipft,crown_damage,canopy_trim,store_c_target_tp1)

               write(fates_log(),*) 'leaf_c: ',leafc_tp1, leaf_c_target_tp1,leafc_tp1-leaf_c_target_tp1
               write(fates_log(),*) 'fnrt_c: ',fnrtc_tp1, fnrt_c_target_tp1,fnrtc_tp1- fnrt_c_target_tp1
               write(fates_log(),*) 'sapw_c: ',sapwc_tp1, sapw_c_target_tp1 ,sapwc_tp1- sapw_c_target_tp1
               write(fates_log(),*) 'store_c: ',storec_tp1, store_c_target_tp1,storec_tp1- store_c_target_tp1
               write(fates_log(),*) 'struct_c: ',structc_tp1, struct_c_target_tp1,structc_tp1- struct_c_target_tp1
               call endrun(msg=errMsg(sourcefile, __LINE__))

            end if if_step_exceedance
               
         end if if_completed_solve
         
      end do do_solve_check

      ! Prioritize nutrient transfer to the reproductive pool
      ! Note, that if we do not keep reproductive tissues on stoichiometry, the seed
      ! pool for that pft will be off stoichiometry, and one of C,N or P will limit
      ! recruitment. Per the current model formulation, new recruits are forced to
      ! have their maximum stoichiometry in each organ. The total stoichiometry
      ! of the recruits should match the stoichiometry of the seeds

      if(prioritize_repro_nutr_growth)then

         target_n = this%GetNutrientTarget(nitrogen_element,repro_organ,stoich_growth_min)
         deficit_n(1) = this%GetDeficit(nitrogen_element,repro_organ,target_n)
         n_flux = max(0._r8,min(n_gain,deficit_n(1)))
         
         target_p = this%GetNutrientTarget(phosphorus_element,repro_organ,stoich_growth_min)
         deficit_p(1) = this%GetDeficit(phosphorus_element,repro_organ,target_p)
         p_flux = max(0._r8,min(p_gain,deficit_p(1)))
         
         this%variables(repro_n_id)%val(1) = this%variables(repro_n_id)%val(1) + n_flux
         this%variables(repro_p_id)%val(1) = this%variables(repro_p_id)%val(1) + p_flux

         n_gain = n_gain - n_flux
         p_gain = p_gain - p_flux
         
      end if
      
      ! -----------------------------------------------------------------------------------
      ! Nutrient Fluxes proportionally to each pool (these should be fully actualized)
      ! (this also removes from the gain pools)
      ! -----------------------------------------------------------------------------------

      sum_n_demand = 0._r8   ! For error checking
      sum_p_demand = 0._r8   ! For error checking
      do ii = 1, n_mask_organs

         i     = mask_organs(ii)
         i_org = mask_gorgans(ii)
         
         target_n = this%GetNutrientTarget(nitrogen_element,i_org,stoich_growth_min)
         target_p = this%GetNutrientTarget(phosphorus_element,i_org,stoich_growth_min)

         deficit_n(ii) = this%GetDeficit(nitrogen_element,i_org,target_n)
         sum_n_demand = sum_n_demand+max(0._r8,deficit_n(ii))
         
         deficit_p(ii) = this%GetDeficit(phosphorus_element,i_org,target_p)
         sum_p_demand = sum_p_demand+max(0._r8,deficit_p(ii))
         
      end do

      ! TODO: mask_organs should be a vector of global organs
      
      ! Nitrogen
      call ProportionalNutrAllocation(this,deficit_n(1:n_mask_organs), & 
           n_gain, nitrogen_element,mask_gorgans(1:n_mask_organs))
      
      ! Phosphorus
      call ProportionalNutrAllocation(this,deficit_p(1:n_mask_organs), &
           p_gain, phosphorus_element,mask_gorgans(1:n_mask_organs))

    end if if_stature_growth

    return
  end subroutine CNPStatureGrowth
  
  ! =====================================================================================

  subroutine CNPAllocateRemainder(this, c_gain,n_gain,p_gain, &
                                  c_efflux, n_efflux, p_efflux, &
                                  target_c,target_dcdd)

    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain
    real(r8), intent(inout) :: n_gain
    real(r8), intent(inout) :: p_gain 
    real(r8), intent(inout) :: c_efflux
    real(r8), intent(inout) :: n_efflux
    real(r8), intent(inout) :: p_efflux
    real(r8)  :: target_c(:)
    real(r8) :: target_dcdd(:)
    
    integer  :: i
    real(r8), dimension(num_organs) :: deficit_n
    real(r8), dimension(num_organs) :: deficit_p
    real(r8) :: target_n
    real(r8) :: target_p
    real(r8) :: store_c_target   ! Target amount of C in storage including "overflow" [kgC]
    real(r8) :: total_c_flux     ! Total C flux from gains into storage and growth R [kgC]
    real(r8), pointer :: dbh
    real(r8), pointer :: resp_excess
    integer           :: ipft
    integer, pointer  :: limiter
    real(r8)          :: canopy_trim
    integer           :: crown_damage
    
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    resp_excess => this%bc_inout(acnp_bc_inout_id_resp_excess)%rval
    limiter     => this%bc_out(acnp_bc_out_id_limiter)%ival
    crown_damage = this%bc_in(acnp_bc_in_id_cdamage)%ival
    
    ! -----------------------------------------------------------------------------------
    ! If nutrients are still available, then we can bump up the values in the pools
    !  towards the OPTIMAL target values.
    ! -----------------------------------------------------------------------------------

    do i = 1, num_organs
       
       ! Update the nitrogen and phosphorus deficits
       target_n = this%GetNutrientTarget(nitrogen_element,l2g_organ_list(i),stoich_growth_min)
       target_p = this%GetNutrientTarget(phosphorus_element,l2g_organ_list(i),stoich_growth_min)

       if(l2g_organ_list(i)==store_organ)then
          target_n = target_n * (1._r8 + prt_params%store_ovrflw_frac(ipft))
          target_p = target_p * (1._r8 + prt_params%store_ovrflw_frac(ipft))
       end if

       deficit_n(i) = max(0._r8,this%GetDeficit(nitrogen_element,l2g_organ_list(i),target_n))
       deficit_p(i) = max(0._r8,this%GetDeficit(phosphorus_element,l2g_organ_list(i),target_p))
          
    end do
    
    ! -----------------------------------------------------------------------------------
    ! Nutrient Fluxes proportionally to each pool (these should be fully actualized)
    ! (this also removes from the gain pools)
    ! -----------------------------------------------------------------------------------
    
    ! Nitrogen
    call ProportionalNutrAllocation(this,deficit_n(1:num_organs), & 
         n_gain, nitrogen_element, l2g_organ_list(1:num_organs))
    
    ! Phosphorus
    call ProportionalNutrAllocation(this,deficit_p(1:num_organs), &
         p_gain, phosphorus_element, l2g_organ_list(1:num_organs))


    ! This routine updates the l2fr (leaf 2 fine-root multiplier) variable
    ! It will also update the target
    call this%CNPAdjustFRootTargets(target_c,target_dcdd)
        
    ! -----------------------------------------------------------------------------------
    ! If carbon is still available, lets cram some into storage overflow
    ! We will do this last, because we wanted the non-overflow storage
    ! value to draw minimum and optimal nutrient fluxes
    ! -----------------------------------------------------------------------------------

    if(c_gain>calloc_abs_error) then

       if(store_c_overflow == retain_c_store_overflow)then
          
          total_c_flux = c_gain
          ! Transfer excess carbon into storage overflow
          this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) + total_c_flux
          c_gain              = c_gain - total_c_flux
          
       elseif(store_c_overflow == burn_c_store_overflow) then

          ! Update carbon based allometric targets
          call bstore_allom(dbh,ipft,crown_damage,canopy_trim, store_c_target)

          ! Allow some overflow
          store_c_target = store_c_target * (1._r8 + prt_params%store_ovrflw_frac(ipft))
          
          total_c_flux = max(0._r8,min(c_gain, store_c_target - this%variables(store_c_id)%val(1) ))
          
          ! Transfer excess carbon INTO storage overflow
          this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) + total_c_flux
          c_gain              = c_gain - total_c_flux

          resp_excess = resp_excess + c_gain
          c_gain      = 0._r8
          
       elseif(store_c_overflow == exude_c_store_overflow)then
                 
          ! Update carbon based allometric targets
          call bstore_allom(dbh,ipft,crown_damage,canopy_trim, store_c_target)
          
          ! Estimate the overflow
          store_c_target = store_c_target * (1._r8 + prt_params%store_ovrflw_frac(ipft))
          
          total_c_flux = max(0.0, min(c_gain, store_c_target - this%variables(store_c_id)%val(1)))
          ! Transfer excess carbon into storage overflow
          this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) + total_c_flux
          c_gain = c_gain - total_c_flux
          
       end if

    end if

    ! If we had some poor numerical precision resulting
    ! in negative gains, use storage to get them back to zero
    ! they should be very very small
    if(c_gain<-nearzero) then
       this%variables(store_c_id)%val(1) = this%variables(store_c_id)%val(1) + c_gain
       c_gain                            = 0
    end if
    if(n_gain<-nearzero) then
       this%variables(store_n_id)%val(1) = this%variables(store_n_id)%val(1) + n_gain
       n_gain                            = 0
    end if
    if(p_gain<-nearzero) then
       this%variables(store_p_id)%val(1) = this%variables(store_p_id)%val(1) + p_gain
       p_gain                            = 0
    end if
    
    

    ! Figure out what to do with excess carbon and nutrients
    ! 1) excude through roots cap at 0 to flush out imprecisions
    ! -----------------------------------------------------------------------------------

    ! If either n or p uptake is in prescribed mode
    ! don't efflux anything, we will use the remainder
    ! n_gain and p_gain to specify the demand as what was used
    ! and what was uptaken
    
    if(n_uptake_mode.eq.prescribed_n_uptake) then
       n_efflux  = 0._r8
    else
       n_efflux = n_gain
       n_gain   = 0._r8
    end if
    
    if(p_uptake_mode.eq.prescribed_p_uptake) then
       p_efflux = 0._r8
    else
       p_efflux = p_gain
       p_gain   = 0._r8
    end if

    c_efflux = c_gain
    c_gain   = 0.0_r8
    
    

    nullify(dbh)

    return
  end subroutine CNPAllocateRemainder

  ! =====================================================================================

  subroutine FastPRTAllometricCNP(this)

    implicit none
    class(cnp_allom_prt_vartypes) :: this     ! this class

    ! This routine does nothing, because
    ! we currently don't have any fast-timestep processes
    ! Think of this as a stub.

    return
  end subroutine FastPRTAllometricCNP

  ! =====================================================================================



  ! =====================================================================================
  
  function GetDeficit(this,element_id,organ_id,target_m) result(deficit_m)

    class(cnp_allom_prt_vartypes) :: this
    integer,intent(in)            :: element_id
    integer,intent(in)            :: organ_id
    real(r8),intent(in)           :: target_m

    integer  :: i_var
    real(r8) :: deficit_m
    
    i_var = prt_global%sp_organ_map(organ_id,element_id)

    if(element_id.eq.carbon12_element) then
       deficit_m = target_m - sum(this%variables(i_var)%val(:),dim=1)
    else
       deficit_m = target_m - this%variables(i_var)%val(1)
    end if
       
    return
  end function GetDeficit
  
  ! =====================================================================================

  function GetNutrientTargetCNP(this,element_id,organ_id,stoich_mode) result(target_m)
    
    class(cnp_allom_prt_vartypes) :: this
    integer, intent(in)           :: element_id
    integer, intent(in)           :: organ_id
    integer, intent(in),optional  :: stoich_mode
    real(r8)                      :: target_m    ! Target amount of nutrient for this organ [kg]
    
    real(r8)         :: target_c
    real(r8),pointer :: dbh
    real(r8)         :: canopy_trim
    real(r8)         :: l2fr
    integer          :: ipft
    integer          :: i_cvar
    integer          :: crown_damage
    real(r8)         :: sapw_area
    real(r8)         :: leaf_c_target,fnrt_c_target
    real(r8)         :: sapw_c_target,agw_c_target
    real(r8)         :: bgw_c_target,struct_c_target
    real(r8)         :: elongf_leaf
    real(r8)         :: elongf_fnrt
    real(r8)         :: elongf_stem
    real(r8)         :: nc_repro,pc_repro


    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    elongf_leaf = this%bc_in(acnp_bc_in_id_efleaf)%rval
    elongf_fnrt = this%bc_in(acnp_bc_in_id_effnrt)%rval
    elongf_stem = this%bc_in(acnp_bc_in_id_efstem)%rval
    i_cvar      = prt_global%sp_organ_map(organ_id,carbon12_element)
    l2fr        = this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    nc_repro    = this%bc_in(acnp_bc_in_id_nc_repro)%rval
    pc_repro    = this%bc_in(acnp_bc_in_id_pc_repro)%rval
    crown_damage = this%bc_in(acnp_bc_in_id_cdamage)%ival
    
    ! Storage of nutrients are assumed to have different compartments than
    ! for carbon, and thus their targets are not associated with a tissue
    ! but is more represented as a fraction of the maximum amount of nutrient
    ! that can be bound in non-reproductive tissues
    
    if(organ_id == store_organ) then

       call bleaf(dbh,ipft,crown_damage,canopy_trim, elongf_leaf, leaf_c_target)
       call bfineroot(dbh,ipft,canopy_trim,l2fr, elongf_fnrt, fnrt_c_target)
       call bsap_allom(dbh,ipft,crown_damage,canopy_trim, elongf_stem, sapw_area,sapw_c_target)
       call bagw_allom(dbh,ipft,crown_damage, elongf_stem, agw_c_target)
       call bbgw_allom(dbh,ipft, elongf_stem, bgw_c_target)
       call bdead_allom(agw_c_target,bgw_c_target, sapw_c_target, ipft, struct_c_target)

       ! Target for storage is a fraction of the sum target of all
       ! non-reproductive organs

       if( element_id == nitrogen_element) then
          
          target_m = StorageNutrientTarget(ipft, element_id, &
               leaf_c_target*prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ)), &
               fnrt_c_target*prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ)), &
               sapw_c_target*prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ)), & 
               struct_c_target*prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(struct_organ)))
       else
          
          target_m = StorageNutrientTarget(ipft, element_id, &
               leaf_c_target*prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ)), &
               fnrt_c_target*prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ)), &
               sapw_c_target*prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ)), & 
               struct_c_target*prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(struct_organ)))

       end if

       ! This is only called during phase 3, remainder and allows
       ! us to have some overflow to avoid exudation/efflux if possible
       if( stoich_mode == stoich_max ) then
          target_m = target_m*(1._r8 + prt_params%store_ovrflw_frac(ipft))
       end if
       
    elseif(organ_id == repro_organ) then

       target_c = this%variables(i_cvar)%val(1)
       if( element_id == nitrogen_element) then
          target_m = target_c * nc_repro
       else
          target_m = target_c * pc_repro
       end if
       
    else


       if(.not.present(stoich_mode))then
          write(fates_log(),*) 'Must specify if nutrient target is growthmin or max'
          write(fates_log(),*) 'for non-reproductive and non-storage organs'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       
       ! In all cases, we want the first index because for non-leaves
       ! that is the only index, and for leaves, that is the newly
       ! growing index.
       target_c = this%variables(i_cvar)%val(1)
    
       if( stoich_mode == stoich_growth_min ) then
          if( element_id == nitrogen_element) then
             target_m = target_c * prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(organ_id))
          else
             target_m = target_c * prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(organ_id))
          end if
       elseif( stoich_mode == stoich_max ) then
          !if( element_id == nitrogen_element) then
          !   target_m = target_c * prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(organ_id))
          !else
          !   target_m = target_c * prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(organ_id))
          !end if
          write(fates_log(),*) 'invalid stoichiometry mode specified while getting'
          write(fates_log(),*) 'nutrient targets'
          write(fates_log(),*) 'stoich_mode: ',stoich_mode
          call endrun(msg=errMsg(sourcefile, __LINE__))
       else
          write(fates_log(),*) 'invalid stoichiometry mode specified while getting'
          write(fates_log(),*) 'nutrient targets'
          write(fates_log(),*) 'stoich_mode: ',stoich_mode
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    nullify(dbh)

    return
  end function GetNutrientTargetCNP


  
  ! =====================================================================================
  
  subroutine ProportionalNutrAllocation(this,deficit_m, gain_m, element_id, list)

    ! -----------------------------------------------------------------------------------
    ! This routine allocates nutrients to a set of organs based on proportional
    ! need.  It is assumed that the input arrays are not sparse, and the fluxes
    ! are based purely off of there deficit from some ideal state.
    ! Note: this may or may not be called inside some preferential organ filter.
    ! -----------------------------------------------------------------------------------

    class(cnp_allom_prt_vartypes) :: this
    real(r8),intent(inout) :: deficit_m(:) ! Nutrient mass deficit of species
                                           ! over set of organs
    integer, intent(in)    :: list(:)      ! List of organ indices from PRTGenericMod
    real(r8),intent(inout) :: gain_m       ! Total nutrient mass gain to
                                           ! work with
    integer,intent(in) :: element_id       ! Element global index

    ! locals
    integer :: num_organs
    integer :: i,i_org
    integer :: i_var
    real(r8) :: flux
    real(r8) :: sum_deficit
    real(r8) :: sum_flux

    num_organs = size(list,dim=1)
       
    sum_deficit = 0._r8
    do i = 1, num_organs
       i_org = list(i)
       sum_deficit = sum_deficit + max(0._r8,deficit_m(i))
    end do
    
    if (sum_deficit>nearzero) then
       
       sum_flux = min(gain_m, sum_deficit)
       
       do i = 1, num_organs
          i_org = list(i)
          
          flux  = sum_flux * max(0._r8,deficit_m(i))/sum_deficit

          i_var = prt_global%sp_organ_map(i_org,element_id)
          this%variables(i_var)%val(1) = this%variables(i_var)%val(1) + flux
          
          deficit_m(i) = deficit_m(i) - flux
          gain_m      = gain_m - flux
          
      end do
      
    end if
    
    if(debug) then
       if(gain_m < -calloc_abs_error) then
          write(fates_log(),*) 'Somehow we have negative nutrient gain'
          write(fates_log(),*) 'during proportional allocation'
          write(fates_log(),*) 'gain: ',gain_m,'element: ',element_id
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    return
  end subroutine ProportionalNutrAllocation

  ! =====================================================================================

  function AllomCNPGrowthDeriv(l_state_array,l_state_mask,cbalance,intgr_params) result(dCdx)

      ! ---------------------------------------------------------------------------------
      ! This function calculates the derivatives for the carbon pools
      ! relative to the amount of carbon balance.
      ! This function is based off of allometry.
      !
      ! Important Note:  While this routine is carbon-only, the total carbon balance
      ! that governs how much each pool is integrated, was likely limited by nutrient
      ! availability.
      !
      ! ---------------------------------------------------------------------------------

      ! Arguments
      real(r8),intent(in), dimension(:) :: l_state_array      ! Vector of carbon pools
                                                        ! dbh,leaf,root,sap,store,dead
      logical,intent(in), dimension(:)  :: l_state_mask ! logical mask of active pools
                                                        ! some may be turned off
      real(r8),intent(in)               :: cbalance     ! The carbon balance of the
                                                        ! partial step (independant var)

      real(r8), intent(in),dimension(:) :: intgr_params  ! Generic Array used to pass
                                                        ! parameters into this function


      ! Return Value
      real(r8),dimension(lbound(l_state_array,dim=1):ubound(l_state_array,dim=1)) :: dCdx

      ! locals
      integer  :: ipft             ! PFT index
      real(r8) :: canopy_trim      ! Canopy trimming function (boundary condition [0-1]
      integer  :: crown_damage     ! Damage class
      real(r8) :: l2fr             ! leaf to fineroot biomass multiplier 
      real(r8) :: leaf_c_target    ! target leaf biomass, dummy var (kgC)
      real(r8) :: fnrt_c_target    ! target fine-root biomass, dummy var (kgC)
      real(r8) :: sapw_c_target    ! target sapwood biomass, dummy var (kgC)
      real(r8) :: agw_c_target     ! target aboveground wood, dummy var (kgC)
      real(r8) :: bgw_c_target     ! target belowground wood, dummy var (kgC)
      real(r8) :: store_c_target   ! target storage, dummy var (kgC)
      real(r8) :: struct_c_target  ! target structural biomas, dummy var (kgC)
      real(r8) :: sapw_area
      real(r8) :: leaf_dcdd_target      ! target leaf biomass derivative wrt d, (kgC/cm)
      real(r8) :: fnrt_dcdd_target      ! target fine-root biomass derivative wrt d, (kgC/cm)
      real(r8) :: sapw_dcdd_target      ! target sapwood biomass derivative wrt d, (kgC/cm)
      real(r8) :: agw_dcdd_target       ! target AG wood biomass derivative wrt d, (kgC/cm)
      real(r8) :: bgw_dcdd_target       ! target BG wood biomass derivative wrt d, (kgC/cm)
      real(r8) :: store_dcdd_target     ! target storage biomass derivative wrt d, (kgC/cm)
      real(r8) :: struct_dcdd_target    ! target structural biomass derivative wrt d, (kgC/cm)
      real(r8) :: repro_fraction        ! fraction of carbon balance directed towards reproduction (kgC/kgC)
      real(r8) :: total_dcostdd         ! carbon cost for non-reproductive pools per unit increment of dbh
      real(r8) :: elongf_leaf           ! Leaf elongation factor (0-1)
      real(r8) :: elongf_fnrt           ! Fine-root "elongation factor" (0-1)
      real(r8) :: elongf_stem           ! Stem "elongation factor" (0-1)


      associate( dbh         => l_state_array(dbh_id),      &
                 leaf_c      => l_state_array(leaf_id),   &
                 fnrt_c      => l_state_array(fnrt_id),   &
                 sapw_c      => l_state_array(sapw_id),   &
                 store_c     => l_state_array(store_id),  &
                 struct_c    => l_state_array(struct_id), &
                 repro_c     => l_state_array(repro_id),  &
                 mask_dbh    => l_state_mask(dbh_id),       &
                 mask_leaf   => l_state_mask(leaf_id),    &
                 mask_fnrt   => l_state_mask(fnrt_id),    &
                 mask_sapw   => l_state_mask(sapw_id),    &
                 mask_store  => l_state_mask(store_id),   &
                 mask_struct => l_state_mask(struct_id),  &
                 mask_repro  => l_state_mask(repro_id) )

        canopy_trim  = intgr_params(intgr_parm_ctrim)
        ipft         = int(intgr_params(intgr_parm_pft))
        l2fr         = intgr_params(intgr_parm_l2fr)
        crown_damage = int(intgr_params(intgr_parm_cdamage))
        elongf_leaf  = intgr_params(intgr_parm_efleaf)
        elongf_fnrt  = intgr_params(intgr_parm_effnrt)
        elongf_stem  = intgr_params(intgr_parm_efstem)

        call bleaf(dbh,ipft,crown_damage,canopy_trim, elongf_leaf, leaf_c_target,leaf_dcdd_target)
        call bfineroot(dbh,ipft,canopy_trim,l2fr, elongf_fnrt, fnrt_c_target,fnrt_dcdd_target)
        call bsap_allom(dbh,ipft,crown_damage,canopy_trim, elongf_stem, sapw_area,sapw_c_target,sapw_dcdd_target)
        call bagw_allom(dbh,ipft,crown_damage, elongf_stem,agw_c_target,agw_dcdd_target)
        call bbgw_allom(dbh,ipft, elongf_stem,bgw_c_target,bgw_dcdd_target)
        call bdead_allom(agw_c_target,bgw_c_target, sapw_c_target, ipft, struct_c_target, &
                         agw_dcdd_target, bgw_dcdd_target, sapw_dcdd_target, struct_dcdd_target)
        call bstore_allom(dbh,ipft,crown_damage,canopy_trim,store_c_target,store_dcdd_target)

        if (mask_repro) then

           ! If the TRS is switched off then we use FATES's default reproductive allocation.
           if ( regeneration_model == default_regeneration .or. &
                prt_params%allom_dbh_maxheight(ipft) < min_max_dbh_for_trees ) then ! The Tree Recruitment Scheme 
                                                                             ! is only for trees
              if (dbh <= prt_params%dbh_repro_threshold(ipft)) then
                 repro_fraction = prt_params%seed_alloc(ipft)
              else
                 repro_fraction = prt_params%seed_alloc(ipft) + prt_params%seed_alloc_mature(ipft)
              end if
   
           ! If the TRS is switched on (with or w/o seedling dynamics) then reproductive allocation is
           ! a pft-specific function of dbh. This allows for the representation of different
           ! reproductive schedules (Wenk and Falster, 2015)
           else if ( any(regeneration_model == [TRS_regeneration, TRS_no_seedling_dyn]) .and. &
                     prt_params%allom_dbh_maxheight(ipft) > min_max_dbh_for_trees ) then

              repro_fraction = prt_params%seed_alloc(ipft) * &
              (exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm) / &
              (1 + exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm)))
           else
              write(fates_log(),*) 'unknown seed allocation and regeneration model, exiting'
              write(fates_log(),*) 'regeneration_model: ',regeneration_model
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if ! regeneration switch 
          
        else ! mask repro
           repro_fraction = 0._r8
        end if !mask repro

        total_dcostdd = 0._r8
        if (mask_struct) then
           total_dcostdd = total_dcostdd + struct_dcdd_target
        end if
        if (mask_leaf) then
           total_dcostdd = total_dcostdd + leaf_dcdd_target
        end if
        if (mask_fnrt) then
           total_dcostdd = total_dcostdd + fnrt_dcdd_target
        end if
        if (mask_sapw) then
           total_dcostdd = total_dcostdd + sapw_dcdd_target
        end if
        if (mask_store) then
           total_dcostdd = total_dcostdd + store_dcdd_target
        end if

        dCdx(:) = 0.0_r8

        ! It is possible that with some asymptotic, or hard
        ! capped allometries, that all growth rates reach zero.
        ! In this case, if there is carbon, give it to reproduction

        if(total_dcostdd > nearzero) then

           if (mask_struct) then
              dCdx(struct_id) =                  struct_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           end if
           if (mask_leaf) then
              dCdx(leaf_id)  =                   leaf_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           end if
           if (mask_fnrt) then
              dCdx(fnrt_id)  =                   fnrt_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           end if
           if (mask_sapw) then
              dCdx(sapw_id)  =                   sapw_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           end if
           if (mask_store) then
              dCdx(store_id) =                   store_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           end if
           if (mask_repro) then
              dCdx(repro_id)  =                  repro_fraction
           end if

           if( abs(sum(dCdx,dim=1)-1.0_r8)>rsnbl_math_prec ) then
              write(fates_log(),*) 'dCdx should sum to 1'
              write(fates_log(),*) 'dCdx: ',dCdx
              write(fates_log(),*) 'repro fraction: ',repro_fraction
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           
           dCdx(dbh_id)      = (1.0_r8/total_dcostdd)*(1.0_r8 - repro_fraction)
           
        else

           if(repro_fraction<nearzero) then
              write(fates_log(),*) 'A plant has reached a stature where'
              write(fates_log(),*) 'allometry is dictating that no more growth'
              write(fates_log(),*) 'should be occuring in any of its non-reporductive pools'
              write(fates_log(),*) 'However, at this point, the fraction sent to reproduction'
              write(fates_log(),*) 'is also zero, which is a conflict. The plant needs'
              write(fates_log(),*) 'to spend its carbon somewhere no?'
              write(fates_log(),*) 'dbh: ',dbh
              write(fates_log(),*) 'pft: ',ipft
              write(fates_log(),*) 'repro fraction: ',repro_fraction
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
              
           dCdx(repro_id) = 1._r8
           
        end if

      end associate

      return
   end function AllomCNPGrowthDeriv

   ! =====================================================================================
   
   
   subroutine EstimateGrowthNC(this,target_c,target_dcdd,state_mask,avg_nc,avg_pc)

     ! This routine predicts the effective nutrient/carbon allocation ratio
     ! for the forthcoming growth step. This helps the growth step predict
     ! which element will be limiting, and reduce the amount of carbon
     ! used to make the step.

     class(cnp_allom_prt_vartypes) :: this
     real(r8)            :: target_c(:)
     real(r8)            :: target_dcdd(:)
     logical             :: state_mask(:)
     real(r8)            :: avg_nc  ! Average N:C ratio
     real(r8)            :: avg_pc  ! Average P:C ratio

     real(r8) :: repro_c_frac
     real(r8) :: total_w      ! Weight (dC/dd) for the ratios
     real(r8) :: store_nc
     real(r8) :: store_pc
     real(r8) :: repro_w,leaf_w,fnrt_w,sapw_w,struct_w,store_w
     
     associate(dbh    => this%bc_inout(acnp_bc_inout_id_dbh)%rval, & 
          ipft        => this%bc_in(acnp_bc_in_id_pft)%ival, &
          nc_repro    => this%bc_in(acnp_bc_in_id_nc_repro)%rval, &
          pc_repro    => this%bc_in(acnp_bc_in_id_pc_repro)%rval)
     
     if(state_mask(repro_id)) then
        
        ! If the TRS is switched off then we use FATES's default reproductive allocation.
        if ( regeneration_model == default_regeneration .or. &
             prt_params%allom_dbh_maxheight(ipft) < min_max_dbh_for_trees ) then ! The Tree Recruitment Scheme 
                                                                                 ! is only for trees
           if (dbh <= prt_params%dbh_repro_threshold(ipft)) then
              repro_c_frac = prt_params%seed_alloc(ipft)
           else
              repro_c_frac = prt_params%seed_alloc(ipft) + prt_params%seed_alloc_mature(ipft)
           end if
   
        ! If the TRS is switched on (with or w/o seedling dynamics) then reproductive allocation is
        ! a pft-specific function of dbh. This allows for the representation of different
        ! reproductive schedules (Wenk and Falster, 2015)
        else if ( any(regeneration_model == [TRS_regeneration, TRS_no_seedling_dyn]) .and. &
                  prt_params%allom_dbh_maxheight(ipft) > min_max_dbh_for_trees ) then

           repro_c_frac = prt_params%seed_alloc(ipft) * &
           (exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm) / &
           (1 + exp(prt_params%repro_alloc_b(ipft) + prt_params%repro_alloc_a(ipft)*dbh*mm_per_cm)))
        else
           write(fates_log(),*) 'unknown seed allocation and regeneration model, exiting'
           write(fates_log(),*) 'regeneration_model: ',regeneration_model
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if ! regeneration switch 

     else ! state mask
        repro_c_frac = 0._r8
     end if ! state mask
        
     ! Estimate the total weight
     total_w = 0._r8
     avg_nc = 0._r8
     avg_pc = 0._r8
     
     if(state_mask(leaf_id)) then
        leaf_w = target_dcdd(leaf_organ) * (1._r8 - repro_c_frac)
        total_w = total_w + leaf_w
        avg_nc = avg_nc + leaf_w * prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
        avg_pc = avg_pc + leaf_w * prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
     end if
     if(state_mask(fnrt_id)) then
        fnrt_w = target_dcdd(fnrt_organ) * (1._r8 - repro_c_frac)
        total_w = total_w + fnrt_w
        avg_nc = avg_nc + fnrt_w * prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
        avg_pc = avg_pc + fnrt_w * prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
     end if
     if(state_mask(sapw_id)) then
        sapw_w = target_dcdd(sapw_organ) * (1._r8 - repro_c_frac)
        total_w = total_w + sapw_w
        avg_nc = avg_nc + sapw_w * prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
        avg_pc = avg_pc + sapw_w * prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
     end if
     if(state_mask(struct_id)) then
        struct_w = target_dcdd(struct_organ) * (1._r8 - repro_c_frac)
        total_w = total_w + struct_w
        avg_nc = avg_nc + struct_w * prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
        avg_pc = avg_pc + struct_w * prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
     end if
     if(state_mask(store_id)) then
        store_w = target_dcdd(store_organ) * (1._r8 - repro_c_frac)
        total_w = total_w + store_w
        store_nc = this%GetNutrientTarget(nitrogen_element,store_organ,stoich_growth_min) / target_c(store_organ)
        store_pc = this%GetNutrientTarget(phosphorus_element,store_organ,stoich_growth_min) / target_c(store_organ)
        avg_nc = avg_nc + store_w * store_nc
        avg_pc = avg_pc + store_w * store_pc
     end if

     if(state_mask(repro_id)) then

        ! repro_w = total*repro_c_frac
        ! repro_w = (total_w + repro_w)*repro_c_frac = total_w*repro_c_frac + repro_w*repro_c_frac
        ! repro_w * (1 - repro_c_frac) = total_w*repro_c_frac
        ! repro_w = total_w * repro_c_frac/(1-repro_c_frac)

        if(1._r8 - repro_c_frac < nearzero) then
           repro_w = repro_c_frac 
        else
           repro_w = total_w * repro_c_frac/(1._r8 - repro_c_frac)
        end if
        
        total_w = total_w  + repro_w
        avg_nc = avg_nc + repro_w * nc_repro
        avg_pc = avg_pc + repro_w * pc_repro
     end if

     avg_nc = avg_nc / total_w
     avg_pc = avg_pc / total_w

   end associate
     
   return
 end subroutine EstimateGrowthNC

end module PRTAllometricCNPMod
