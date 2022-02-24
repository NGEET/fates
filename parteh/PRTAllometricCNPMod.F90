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
  
  use PRTGenericMod , only  : leaf_organ
  use PRTGenericMod , only  : fnrt_organ
  use PRTGenericMod , only  : sapw_organ
  use PRTGenericMod , only  : store_organ
  use PRTGenericMod , only  : repro_organ
  use PRTGenericMod , only  : struct_organ
  use PRTGenericMod , only  : all_organs
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
  use FatesIntegratorsMod , only : RKF45
  use FatesIntegratorsMod , only : Euler
  use FatesConstantsMod   , only : calloc_abs_error
  use FatesConstantsMod   , only : nearzero
  use FatesConstantsMod   , only : itrue
  use FatesConstantsMod   , only : fates_unset_r8
  use FatesConstantsMod   , only : fates_unset_int
  use FatesConstantsMod   , only : sec_per_day
  use PRTParametersMod    , only : prt_params
  use EDTypesMod          , only : leaves_on,leaves_off
  
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
  integer,public, parameter :: stoich_max = 2            ! Flag for stoichiometry associated with 
                                                  ! maximum for that organ

  
  ! This is the ordered list of organs used in this module
  ! -------------------------------------------------------------------------------------

  integer, parameter                        :: num_organs = 6
  
  ! Converting from local to global organ id
  integer, parameter,dimension(num_organs) :: organ_list = & 
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
  integer, public, parameter :: acnp_bc_inout_id_rmaint_def = 2  ! Index for any accumulated
                                                                 ! maintenance respiration deficit
  integer, public, parameter :: acnp_bc_inout_id_l2fr       = 3  ! leaf 2 fineroot scalar, this
                                                                 ! is dynamic with CNP
  integer, public, parameter :: num_bc_inout                = 3

  ! -------------------------------------------------------------------------------------
  ! Input only Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: acnp_bc_in_id_pft      = 1 ! Index for the PFT input BC
  integer, public, parameter :: acnp_bc_in_id_ctrim    = 2 ! Index for the canopy trim function
  integer, public, parameter :: acnp_bc_in_id_lstat    = 3 ! phenology status logical
  integer, public, parameter :: acnp_bc_in_id_netdc    = 4 ! Index for the net daily C input BC
  integer, public, parameter :: acnp_bc_in_id_netdnh4  = 5 ! Index for the net daily NH4 input BC
  integer, public, parameter :: acnp_bc_in_id_netdno3  = 6 ! Index for the net daily NO3 input BC
  integer, public, parameter :: acnp_bc_in_id_netdp    = 7 ! Index for the net daily P input BC
  !integer, public, parameter :: acnp_bc_in_id_l2fr_ema = 7 ! Index for the moving average ema
  !integer, public, parameter :: acnp_bc_in_id_ncs_ema = 9  ! Index for N/C storage ratio (EMA)
  !integer, public, parameter :: acnp_bc_in_id_pcs_ema = 10 ! Index for P/C storage ratio (EMA)
  
  ! 0=leaf off, 1=leaf on
  integer, parameter         :: num_bc_in             = 7

  ! -------------------------------------------------------------------------------------
  ! Output Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: acnp_bc_out_id_cefflux = 1  ! Daily exudation of C  [kg]
  integer, public, parameter :: acnp_bc_out_id_nefflux = 2  ! Daily exudation of N  [kg]
  integer, public, parameter :: acnp_bc_out_id_pefflux = 3  ! Daily exudation of P  [kg]
  
  integer, parameter         :: num_bc_out                = 3  ! Total number of


  ! Indices for parameters passed to the integrator
  integer,private, parameter :: intgr_parm_ctrim = 1
  integer,private, parameter :: intgr_parm_pft   = 2
  integer,private, parameter :: intgr_parm_l2fr  = 3
  integer,private, parameter :: num_intgr_parm   = 3
  
  ! -------------------------------------------------------------------------------------
  ! Define the size of the coorindate vector.  For this hypothesis, there is only
  ! one pool per each species x organ combination, except for leaves (WHICH HAVE AGE)
  ! -------------------------------------------------------------------------------------
  integer, parameter :: icd = 1


  real(r8), parameter :: store_overflow_frac = 0.15      ! The fraction above target allowed in storage

  logical, parameter :: force_store_c_overflow = .true.

  ! User may want to attempt matching results with the
  ! C-only allocation module. If so, then set reproduce_conly
  ! and make sure both fnrt and leaf are set to the highest
  ! priority order, sapwood and storage are set to the
  ! second highest, and then structure is last. When this is
  ! flagged as true, it changes the logic in the first allocation
  ! phase, to give first dibs to leaves, even though they are
  ! in the same priority group as fineroots.

  logical, parameter :: reproduce_conly = .false.

  ! Definitions for the regulation functions. These typically translate
  ! a storage fraction into a scalar used to regulate resources somehow
  
  integer, parameter :: regulate_linear = 1
  integer, parameter :: regulate_logi   = 2
  integer, parameter :: regulate_CN_logi = 3
  integer, parameter :: regulate_CN_dfdd = 4
  integer, parameter :: regulate_CN_ema  = 5


  real(r8), public, parameter :: fnrt_adapt_tscl = 100._r8   ! Fine-root adaptation timescale (days)
                                                             ! or, how many days it takes
                                                             ! for a doubling or halving of the l2fr



  
  ! Array of pointers are difficult in F90
  ! This structure is a necessary intermediate 
  type :: parray_type
     real(r8), pointer :: ptr
  end type parray_type

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
     procedure :: CNPAdjustFRootTargets
     procedure :: CNPAllocateRemainder
     procedure :: GetDeficit
     procedure :: StorageRegulator
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

     allocate(prt_global_acnp)
     allocate(prt_global_acnp%state_descriptor(num_vars))

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


  subroutine DailyPRTAllometricCNP(this)

    class(cnp_allom_prt_vartypes) :: this
    
    ! Pointers to in-out bcs
    real(r8),pointer :: dbh          ! Diameter at breast height [cm]
    real(r8),pointer :: maint_r_def  ! Current maintenance respiration deficit [kgC]
    real(r8),pointer :: l2fr         ! Leaf to fineroot ratio of target biomass
    
    ! Input only bcs
    integer  :: ipft        ! Plant Functional Type index
    real(r8) :: c_gain      ! Daily carbon balance for this cohort [kgC]
    real(r8) :: n_gain      ! Daily nitrogen uptake through fine-roots [kgN]
    real(r8) :: p_gain      ! Daily phosphorus uptake through fine-roots [kgN]
    real(r8) :: canopy_trim ! The canopy trimming function [0-1]
    !real(r8) :: l2fr_ema    ! Mean (EMA) l2fr
    
    ! Pointers to output bcs
    real(r8),pointer :: c_efflux   ! Total plant efflux of carbon (kgC)
    real(r8),pointer :: n_efflux   ! Total plant efflux of nitrogen (kgN)
    real(r8),pointer :: p_efflux   ! Total plant efflux of phosphorus (kgP)
    real(r8),pointer :: growth_r   ! Total plant growth respiration this step (kgC)

    ! These are pointers to the state variables, rearranged in organ dimensioned
    ! arrays.  This is useful because we loop through organs so often
    type(parray_type),pointer :: state_c(:)   ! State array for carbon, by organ [kg]
    type(parray_type),pointer :: state_n(:)   ! State array for N, by organ [kg]
    type(parray_type),pointer :: state_p(:)   ! State array for P, by organ [kg]

    integer :: i_org   ! organ index
    integer :: i_var   ! variable index

    ! Agruments for allometry functions, that are not in the target_c array
    real(r8) :: agw_c_target,agw_dcdd_target
    real(r8) :: bgw_c_target,bgw_dcdd_target
    real(r8) :: sapw_area
    integer  :: cnp_limiter
    real(r8) :: max_store_n
    ! These arrays hold various support variables dimensioned by organ
    ! Zero suffix indicates the initial state values at the beginning of the routine
    ! _unl suffix indicates values used for tracking nutrient need (ie unlimited)
    ! target is the target masses associated with the plant stature, and
    ! also the derivative of c wrt diameter at current diameter
    real(r8), dimension(num_organs) :: target_c, target_dcdd
    real(r8), dimension(num_organs) :: state_c0, state_n0, state_p0

    ! These are daily mass gains, frozen in time, not drawn from, and thus
    ! these are only used for evaluating mass balancing at the end
    real(r8) :: dbh0
    real(r8) :: c_gain0
    real(r8) :: n_gain0
    real(r8) :: p_gain0
    real(r8) :: maint_r_def0

    ! Used for mass checking, total mass allocated based
    ! on change in the states, should match gain0's
    real(r8) :: allocated_c
    real(r8) :: allocated_n
    real(r8) :: allocated_p
    real(r8) :: target_n,target_p
    real(r8) :: sum_c ! error checking sum



    ! ===================================================================================
    ! Step 1: Evaluate nutrient storage in the plant. Depending on how low
    ! these stores are, we will move proportionally more or less of the daily carbon
    ! gain to increase the target fine-root biomass, fill up to target
    ! and then attempt to get them up to stoichiometry targets.
    ! ===================================================================================

    ! This routine actually just updates the l2fr variable
    call this%CNPAdjustFRootTargets()

    
    ! integrator variables

    ! Copy the input only boundary conditions into readable local variables
    ! We don't use pointers, because inputs should be intent in only
    ! Also, we save the initial values of many of these BC's
    ! for checking and resetting if needed
    ! -----------------------------------------------------------------------------------
    c_gain      = this%bc_in(acnp_bc_in_id_netdc)%rval; c_gain0      = c_gain
    n_gain      = this%bc_in(acnp_bc_in_id_netdnh4)%rval + &
                  this%bc_in(acnp_bc_in_id_netdno3)%rval
    n_gain0      = n_gain
    p_gain      = this%bc_in(acnp_bc_in_id_netdp)%rval; p_gain0      = p_gain
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    
    
    ! Output only boundary conditions
    c_efflux    => this%bc_out(acnp_bc_out_id_cefflux)%rval;  c_efflux = 0._r8
    n_efflux    => this%bc_out(acnp_bc_out_id_nefflux)%rval;  n_efflux = 0._r8
    p_efflux    => this%bc_out(acnp_bc_out_id_pefflux)%rval;  p_efflux = 0._r8

    ! In/out boundary conditions
    maint_r_def => this%bc_inout(acnp_bc_inout_id_rmaint_def)%rval; maint_r_def0 = maint_r_def
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval;        dbh0         = dbh
    l2fr        => this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    
    
    ! If more than 1 leaf age bin is present, this
    ! call advances leaves in their age, but does
    ! not actually remove any biomass from the plant
    
    call this%AgeLeaves(ipft,sec_per_day)

    
    ! Set all of the per-organ pointer arrays
    ! Note: Since growth only happens in the 1st leaf bin, we only
    ! point to that bin.  However, we need to account for all bins
    ! when we calculate the deficit
    
    allocate(state_c(num_organs))
    allocate(state_n(num_organs))
    allocate(state_p(num_organs))

    ! Set carbon targets based on the plant's current stature
    target_c(:) = fates_unset_r8
    target_dcdd(:) = fates_unset_r8
    call bsap_allom(dbh,ipft,canopy_trim,sapw_area,target_c(sapw_id),target_dcdd(sapw_id)  )
    call bagw_allom(dbh,ipft,agw_c_target,agw_dcdd_target)
    call bbgw_allom(dbh,ipft,bgw_c_target,bgw_dcdd_target)
    call bdead_allom(agw_c_target,bgw_c_target, target_c(sapw_id), ipft, target_c(struct_id), &
                     agw_dcdd_target, bgw_dcdd_target, target_dcdd(sapw_id), target_dcdd(struct_id))
    call bleaf(dbh,ipft,canopy_trim, target_c(leaf_id), target_dcdd(leaf_id))
    call bfineroot(dbh,ipft,canopy_trim, l2fr, target_c(fnrt_id), target_dcdd(fnrt_id))
    call bstore_allom(dbh,ipft,canopy_trim, target_c(store_id), target_dcdd(store_id))
    target_c(repro_id) = 0._r8
    target_dcdd(repro_id) = 0._r8

    ! Initialize the the state, and keep a record of this state
    ! as we may actuall run the allocation process twice, and
    ! will need this state to both reset, and measure total
    ! mass fluxes
    do i_org = 1,num_organs

       i_var = prt_global%sp_organ_map(organ_list(i_org),carbon12_element)
       state_c(i_org)%ptr => this%variables(i_var)%val(1)
       state_c0(i_org)  = this%variables(i_var)%val(1)

       i_var = prt_global%sp_organ_map(organ_list(i_org),nitrogen_element)
       state_n(i_org)%ptr => this%variables(i_var)%val(1)
       state_n0(i_org)  = this%variables(i_var)%val(1)

       i_var = prt_global%sp_organ_map(organ_list(i_org),phosphorus_element)
       state_p(i_org)%ptr => this%variables(i_var)%val(1)
       state_p0(i_org)  =  this%variables(i_var)%val(1)
       
    end do

    
    call bfineroot(dbh,ipft,canopy_trim, l2fr, target_c(fnrt_id), target_dcdd(fnrt_id))
    
    ! ===================================================================================
    ! Step 0.  Transfer all stored nutrient into the daily uptake pool.
    ! ===================================================================================

    !i_var = prt_global%sp_organ_map(store_organ,carbon12_element)
    !c_gain = c_gain + max(0._r8,sum(this%variables(i_var)%val(:))-target_c(store_id))
    !this%variables(i_var)%val(1) = this%variables(i_var)%val(1)-max(0._r8,sum(this%variables(i_var)%val(:))-target_c(store_id))
    
    i_var = prt_global%sp_organ_map(store_organ,nitrogen_element)
    n_gain = n_gain + sum(this%variables(i_var)%val(:))
    this%variables(i_var)%val(:) = 0._r8

    i_var = prt_global%sp_organ_map(store_organ,phosphorus_element)
    p_gain = p_gain + sum(this%variables(i_var)%val(:))
    this%variables(i_var)%val(:) = 0._r8
    
    ! ===================================================================================
    ! Step 2.  Prioritized allocation to replace tissues from turnover, and/or pay
    ! any un-paid maintenance respiration from storage.
    ! ===================================================================================
    
    call this%CNPPrioritizedReplacement(maint_r_def, c_gain, n_gain, p_gain, &
             state_c, state_n, state_p, target_c)

    
    sum_c = 0._r8
    do i_org = 1,num_organs
       sum_c = sum_c+state_c(i_org)%ptr
    end do
    if( abs((c_gain0-c_gain) - &
            (sum_c-sum(state_c0(:),dim=1)+(maint_r_def0-maint_r_def))) >calloc_abs_error ) then
       write(fates_log(),*) 'Carbon not balancing I'
       do i_org = 1,num_organs
          write(fates_log(),*) 'state_c: ',state_c(i_org)%ptr,state_c0(i_org)
       end do
       write(fates_log(),*) maint_r_def0-maint_r_def
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    ! ===================================================================================
    ! Step 3. Grow out the stature of the plant by allocating to tissues beyond
    ! current targets. 
    ! Attempts have been made to get all pools and species closest to allometric
    ! targets based on prioritized relative demand and allometry functions.
    ! ===================================================================================

    call this%CNPStatureGrowth(c_gain, n_gain, p_gain, &
         state_c, state_n, state_p, target_c, target_dcdd, cnp_limiter)
    
    sum_c = 0._r8
    do i_org = 1,num_organs
       sum_c = sum_c+state_c(i_org)%ptr
    end do
    if( abs((c_gain0-c_gain) - &
            (sum_c-sum(state_c0(:),dim=1)+(maint_r_def0-maint_r_def))) >calloc_abs_error ) then
       write(fates_log(),*) 'Carbon not balanceing II'
       do i_org = 1,num_organs
          write(fates_log(),*) 'state_c: ',state_c(i_org)%ptr,state_c0(i_org)
       end do
       write(fates_log(),*) maint_r_def0-maint_r_def
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! ===================================================================================
    ! Step 3. 
    ! At this point, at least 1 of the 3 resources have been used up.
    ! Allocate the remaining resources, or as a last resort, efflux them.
    ! ===================================================================================
    
    call this%CNPAllocateRemainder(c_gain, n_gain, p_gain,  &
         state_c, state_n, state_p, c_efflux, n_efflux, p_efflux)

    ! Error Check: Make sure that the mass gains are completely used up
    if( abs(c_gain) > calloc_abs_error .or. &
        abs(n_gain) > 0.1_r8*calloc_abs_error .or. &
        abs(p_gain) > 0.02_r8*calloc_abs_error ) then
       write(fates_log(),*) 'Allocation scheme should had used up all mass gain pools'
       write(fates_log(),*) 'Any mass that cannot be allocated should be effluxed'
       write(fates_log(),*) 'c_gain: ',c_gain
       write(fates_log(),*) 'n_gain: ',n_gain
       write(fates_log(),*) 'p_gain: ',p_gain
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if


    ! Perform a final tally on what was used (allocated)
    ! Since this is also a check against what was available
    ! we include maintenance pay-back and efflux to the "allocated"
    ! pool to make sure everything balances.
    
    allocated_c = (maint_r_def0-maint_r_def) + c_efflux
    allocated_n = n_efflux
    allocated_p = p_efflux
    
    ! Update the allocation flux diagnostic arrays for each 3 elements
    do i_org = 1,num_organs
       
       i_var = prt_global%sp_organ_map(organ_list(i_org),carbon12_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (state_c(i_org)%ptr - state_c0(i_org))

       allocated_c = allocated_c + (state_c(i_org)%ptr - state_c0(i_org))
       
       i_var = prt_global%sp_organ_map(organ_list(i_org),nitrogen_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (state_n(i_org)%ptr - state_n0(i_org))

       allocated_n = allocated_n + (state_n(i_org)%ptr - state_n0(i_org))
       
       i_var = prt_global%sp_organ_map(organ_list(i_org),phosphorus_element)
       this%variables(i_var)%net_alloc(1) = &
            this%variables(i_var)%net_alloc(1) + (state_p(i_org)%ptr - state_p0(i_org))

       allocated_p = allocated_p + (state_p(i_org)%ptr - state_p0(i_org))
       
    end do
    
    if(debug) then

       ! Error Check: Do a final balance between how much mass
       ! we had to work with, and how much was allocated
       
       if ( abs(allocated_c - c_gain0) > calloc_abs_error .or. & 
            abs(allocated_n - n_gain0) > calloc_abs_error .or. &
            abs(allocated_p - p_gain0) > calloc_abs_error ) then
          write(fates_log(),*) 'CNP allocation scheme did not balance mass.'
          write(fates_log(),*) 'c_gain0: ',c_gain0,' allocated_c: ',allocated_c
          write(fates_log(),*) 'n_gain0: ',n_gain0,' allocated_n: ',allocated_n
          write(fates_log(),*) 'p_gain0: ',p_gain0,' allocated_p: ',allocated_p

          do i_org = 1,num_organs
             write(fates_log(),*) i_org, state_c(i_org)%ptr-state_c0(i_org)
          end do
          write(fates_log(),*) (maint_r_def0-maint_r_def), c_efflux
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    !target_n = this%GetNutrientTarget(nitrogen_element,store_organ,stoich_max)
    !target_p = this%GetNutrientTarget(phosphorus_element,store_organ,stoich_max)
    
    deallocate(state_c)
    deallocate(state_n)
    deallocate(state_p)
    
    return
  end subroutine DailyPRTAllometricCNP

  ! =====================================================================================
  subroutine CNPAdjustFRootTargets(this)

    class(cnp_allom_prt_vartypes) :: this

    real(r8), pointer :: l2fr           ! leaf to fineroot target biomass scaler
    integer           :: ipft           ! PFT index
    real(r8), pointer :: dbh
    real(r8)          :: canopy_trim
    !!real(r8)          :: l2fr_ema       ! Moving average L2FR (EMA)

    real(r8) :: l2fr_actual
    real(r8) :: fnrt_c_target ! Target fineroot C before we change l2fr
    real(r8) :: leaf_c_target
    real(r8) :: n_regulator   ! Nitrogen storage regulation function scaler
    real(r8) :: p_regulator   ! Phosphorus storage regulation function scaler
    real(r8) :: np_regulator  ! Combined NP storage regulation function scaler
    real(r8) :: turnfrac      ! The factional amount of root biomass that may
                              ! remain after a day of root turnover without
                              ! replacement
    real(r8) :: fnrt_frac     ! fine-root's current fraction of the target
    real(r8) :: loss_flux_c
    real(r8) :: loss_flux_n
    real(r8) :: loss_flux_p
    real(r8) :: fnrt_c_above_target
    
    integer, parameter :: regulate_type = regulate_CN_ema

    ipft         =  this%bc_in(acnp_bc_in_id_pft)%ival
    l2fr         => this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    dbh          => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim  = this%bc_in(acnp_bc_in_id_ctrim)%rval
    !!l2fr_ema     = this%bc_in(acnp_bc_in_id_l2fr_ema)%rval
    
    associate( l2fr_min => prt_params%allom_l2fr_min(ipft), &
         l2fr_max => prt_params%allom_l2fr_max(ipft))

      call this%StorageRegulator(nitrogen_element, regulate_type,n_regulator)
      call this%StorageRegulator(phosphorus_element, regulate_type,p_regulator)

      ! We take the maximum here, because the maximum is reflective of the
      ! element with the lowest storage, which is the limiting element
      
      np_regulator = max(n_regulator,p_regulator)

      ! Update the leaf-to-fineroot ratio used
      ! to set fine-root biomass allometry

      if(regulate_type == regulate_CN_logi)then

        l2fr = l2fr_min + max(0._r8,min(1.0_r8,np_regulator))*(l2fr_max-l2fr_min)

         
      elseif(regulate_type == regulate_CN_ema) then

         ! To prevent the target l2fr from diverging too far from the
         ! actual l2fr, create some constraints.
         l2fr_actual = this%GetState(fnrt_organ, carbon12_element)/this%GetState(leaf_organ, carbon12_element)

         ! Only update L2FR if some leaves are out
         call bleaf(dbh,ipft,canopy_trim,leaf_c_target)
         if(this%GetState(leaf_organ, carbon12_element)/leaf_c_target>0.5_r8) then
            l2fr = l2fr * np_regulator
         end if
         
         
      elseif(regulate_type == regulate_CN_dfdd)then

         ! Also, we aren't allowed to increase root biomass target if
         ! we are very low on root biomass relative to the target
         ! And we aren't allowed to reduce the target if we are above the target
         
         call bfineroot(dbh,ipft,canopy_trim,l2fr,fnrt_c_target)
         fnrt_frac = this%GetState(fnrt_organ, carbon12_element)/fnrt_c_target

         turnfrac = (years_per_day / prt_params%root_long(ipft))
         
         ! If there is low root compared to the max, don't allow cap growth
         if(fnrt_frac < 1._r8-1.5_r8*turnfrac)then
            np_regulator = min(np_regulator,1.0)
         end if
         
         ! If there is high root compared to the max, don't allow cap decrease
         if(fnrt_frac > 1.0_r8+1.5_r8*turnfrac)then
            np_regulator = max(np_regulator,1.0)
         end if

         ! Don't allow us to drop l2fr more than what the maximum loss to turnover
         ! would be for one day.
         ! this will prevent the algorithm from snowballing.  This is doubly important
         ! because if C is low compared to N or P, then the plant is probably
         ! not very productive, and will not be growing. A growing plant can reach
         ! equilibrium root mass more quickly. (might be unnecessary given
         ! the growth caps prior to this...?)
         
         l2fr = l2fr * min(max(np_regulator,1._r8-2._r8*turnfrac),1._r8+2._r8*turnfrac)
         
      end if
      
      ! Find the updated target fineroot biomass
      call bfineroot(dbh,ipft,canopy_trim, l2fr, fnrt_c_target)

      fnrt_c_above_target = max(0._r8,this%GetState(fnrt_organ, carbon12_element) - fnrt_c_target)

      ! Allow no reabsorption? (any reabsorption of nutrients will further push the N/C or P/C imbalance
      ! and if we are dropping roots, its because we had excess nutrient compared to carbon anyway
      ! Stop the positive feedback

      ! Since this is really a stop gap, we want to allow allocation to handle most of the
      ! fine-root adaption, and only have this kick in when the target starts to drift significantly
      ! from the actual

      loss_flux_c = 0._r8
      
      !loss_flux_c = fnrt_c_above_target*max(fnrt_c_above_target/fnrt_c_target-0.1_r8,0._r8)
      
      
      loss_flux_n = loss_flux_c*this%variables(fnrt_n_id)%val(1)/this%variables(fnrt_c_id)%val(1)
      this%variables(fnrt_n_id)%val(1) = this%variables(fnrt_n_id)%val(1) - loss_flux_n
      this%variables(fnrt_n_id)%turnover(1) = this%variables(fnrt_n_id)%turnover(1) + loss_flux_p

      loss_flux_p = loss_flux_c*this%variables(fnrt_p_id)%val(1)/this%variables(fnrt_c_id)%val(1)
      this%variables(fnrt_p_id)%val(1) = this%variables(fnrt_p_id)%val(1) - loss_flux_p
      this%variables(fnrt_p_id)%turnover(1) = this%variables(fnrt_p_id)%turnover(1) + loss_flux_p
      
      this%variables(fnrt_c_id)%val(1) = this%variables(fnrt_c_id)%val(1) - loss_flux_c
      this%variables(fnrt_c_id)%turnover(1) = this%variables(fnrt_c_id)%turnover(1) + loss_flux_c

      

    end associate

    return
  end subroutine CNPAdjustFRootTargets

  ! =====================================================================================
    
  subroutine CNPPrioritizedReplacement(this, & 
       maint_r_deficit, c_gain, n_gain, p_gain, &
       state_c, state_n, state_p, target_c)

      
    ! -----------------------------------------------------------------------------------
    ! Alternative allocation hypothesis for the prioritized replacement phase.
    ! This is more similar to the current (04/2020) carbon only hypothesis.
    ! -----------------------------------------------------------------------------------
    
    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain
    real(r8), intent(inout) :: n_gain
    real(r8), intent(inout) :: p_gain
    real(r8), intent(inout) :: maint_r_deficit  ! Not currently used
    type(parray_type) :: state_c(:)     ! State array for carbon, by organ [kg]
    type(parray_type) :: state_n(:)     ! State array for N, by organ [kg]
    type(parray_type) :: state_p(:)     ! State array for P, by organ [kg]
    real(r8), intent(in) :: target_c(:)
    
    integer  :: n_curpri_org
    integer, dimension(num_organs) :: curpri_org         ! C variable ID's of the current priority level
    real(r8), dimension(num_organs) :: deficit_c ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_n ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_p ! Deficit to get to target from current    [kg]
    integer  :: i, ii, i_org             ! Loop indices (mostly for organs)
    integer  :: i_cvar                   ! variable index
    integer  :: i_pri                    ! loop index for priority
    integer  :: ipft                     ! Plant functional type index of this plant
    integer  :: leaf_status              ! Is this plant in a leaf on or off status?
    real(r8) :: dbh                      ! DBH [cm]
    real(r8) :: canopy_trim              ! trim factor for maximum leaf biomass
    real(r8) :: target_n                 ! Target mass of N for a given organ [kg]
    real(r8) :: target_p                 ! Target mass of P for a given organ [kg]
    integer  :: priority_code            ! Index for priority level of each organ
    real(r8) :: sum_c_demand            ! Carbon demanded to bring tissues up to allometry (kg)
    real(r8) :: sum_n_deficit            ! The nitrogen deficit of all pools for given priority level (kg)
    real(r8) :: sum_p_deficit            ! The phosphorus deficit of all pools for given priority level (kg)
    real(r8) :: store_below_target
    real(r8) :: store_target_fraction
    real(r8) :: store_demand
    real(r8) :: store_c_flux
    real(r8) :: sum_c_flux               ! The flux to bring tissues up to allometry (kg)
    real(r8) :: sum_n_flux               ! The flux of nitrogen ""  (kg)
    real(r8) :: sum_p_flux               ! The flux of phosphorus "" (Kg)
    real(r8) :: c_flux                   ! carbon flux into an arbitrary pool (kg)
    real(r8) :: gr_flux                  ! carbon flux to fulfill growth respiration of an arbitrary pool (kg)
    real(r8) :: n_flux                   ! nitrogen flux into  an arbitrary pool (kg)
    real(r8) :: p_flux                   ! phosphorus flux into an arbitrary pool (kg)
    real(r8) :: maint_r_def_flux         ! Flux into maintenance respiration during priority 1 allocation
    real(r8) :: c_gain_flux              ! Flux used to pay back negative carbon gain (from storage) (kgC)
    real(r8) :: sapw_area
    integer, parameter  :: n_max_priority = num_organs + 1 ! Maximum possible number of priority levels is
                                                           ! the total number organs plus 1, which allows
                                                           ! each organ to have its own level, and ignore
                                                           ! the specialized priority 1
    
    leaf_status     = this%bc_in(acnp_bc_in_id_lstat)%ival
    ipft            = this%bc_in(acnp_bc_in_id_pft)%ival
    canopy_trim     = this%bc_in(acnp_bc_in_id_ctrim)%rval
    
    ! -----------------------------------------------------------------------------------
    ! Preferential transfer of available carbon and nutrients into the highest
    ! priority pools, and maintenance respiration. We will loop through the available
    ! pools, and identify if that pool is part of the highest transfer priority.
    ! If it is, then we track the variable ids associated with that pool for each CNP
    ! species.  It "should" work fine if there are NO priority=1 pools...
    ! -----------------------------------------------------------------------------------
    
    curpri_org(:) = fates_unset_int    ! reset "current-priority" organ ids
    i = 0
    do ii = 1, num_organs
       
       deficit_c(ii) = max(0._r8,this%GetDeficit(carbon12_element,organ_list(ii),target_c(ii)))

       ! The following logic bars any organs that were not given allocation priority
       if( prt_params%organ_param_id(organ_list(ii)) < 1 ) cycle

       ! The priority code associated with this organ
       priority_code = int(prt_params%alloc_priority(ipft, prt_params%organ_param_id(organ_list(ii))))
          
       ! Don't allow allocation to leaves if they are in an "off" status.
       ! Also, dont allocate to replace turnover if this is not evergreen
       ! (this prevents accidental re-flushing on the day they drop)
       if( ((leaf_status.eq.leaves_off) .or. (prt_params%evergreen(ipft) .ne. itrue)) &
            .and. (organ_list(ii).eq.leaf_organ)) cycle
       
       ! 1 is the highest priority code possible
       if( priority_code == 1 ) then
          i = i + 1
          curpri_org(i) = ii
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
       ii = curpri_org(i)

       i_cvar       = prt_global%sp_organ_map(organ_list(ii),carbon12_element)
       sum_c_demand = sum_c_demand + prt_params%leaf_stor_priority(ipft) * &
            sum(this%variables(i_cvar)%turnover(:))
          
    end do


    sum_c_flux = max(0._r8,min(sum_c_demand,state_c(store_id)%ptr+c_gain))
    
    if (sum_c_flux> nearzero ) then
       
       ! We pay this even if we don't have the carbon
       ! Just don't pay so much carbon that storage+carbon_balance can't pay for it
       
       do ii = 1,n_curpri_org
          i = curpri_org(ii)
          
          i_cvar = prt_global%sp_organ_map(organ_list(i),carbon12_element)

          c_flux = sum_c_flux*(prt_params%leaf_stor_priority(ipft) * &
               sum(this%variables(i_cvar)%turnover(:))/sum_c_demand)

          ! Add carbon to the pool
          state_c(i)%ptr = state_c(i)%ptr + c_flux
          
          ! Remove from daily  carbon gain
          c_gain = c_gain - c_flux

       end do
    end if

    ! Determine nutrient demand and make tansfers
    do i = 1, n_curpri_org
       
       i_org = curpri_org(i)
       
       ! Update the nitrogen deficits
       ! Note that the nitrogen target is tied to the stoichiometry of thegrowing pool only

       target_n = this%GetNutrientTarget(nitrogen_element,organ_list(i_org),stoich_growth_min)
       deficit_n(i_org) = max(0.0_r8, target_n - state_n(i_org)%ptr )
       
       ! Update the phosphorus deficits (which are based off of carbon actual..)
       ! Note that the phsophorus target is tied to the stoichiometry of thegrowing pool only (also)
       target_p = this%GetNutrientTarget(phosphorus_element,organ_list(i_org),stoich_growth_min)
       deficit_p(i_org) = max(0.0_r8, target_p - state_p(i_org)%ptr )

    end do
    
    ! Allocate nutrients at this priority level
    ! Nitrogen
    call ProportionalNutrAllocation(state_n, deficit_n, &
         n_gain, nitrogen_element, curpri_org(1:n_curpri_org))
    
    ! Phosphorus
    call ProportionalNutrAllocation(state_p, deficit_p, &
         p_gain, phosphorus_element, curpri_org(1:n_curpri_org))
    

    
    ! -----------------------------------------------------------------------------------
    ! IV. if carbon balance is negative, re-coup the losses from storage
    !       if it is positive, give some love to storage carbon
    ! -----------------------------------------------------------------------------------
    
    if( c_gain < 0.0_r8 ) then

       ! Storage will have to pay for any negative gains
       store_c_flux           = -c_gain
       c_gain                 = c_gain                + store_c_flux
       state_c(store_id)%ptr  = state_c(store_id)%ptr - store_c_flux
       
    else

       ! This is just a cap, don't fill up more than is needed (shouldn't even apply)
       store_below_target     = max(target_c(store_id) - state_c(store_id)%ptr,0._r8)
       
       ! This is the desired need for carbon
       store_target_fraction  = max(state_c(store_id)%ptr/target_c(store_id),0._r8)
       
       store_demand           = max(c_gain*(exp(-1.*store_target_fraction**4._r8) - exp( -1.0_r8 )),0._r8)

       ! The flux is the (positive) minimum of all three
       store_c_flux           = min(store_below_target,store_demand)
       
       c_gain                 = c_gain  - store_c_flux
       state_c(store_id)%ptr  = state_c(store_id)%ptr + store_c_flux
       
   
   end if
   
    
    ! -----------------------------------------------------------------------------------
    !  If carbon is still available, allocate to remaining high
    !        carbon balance is guaranteed to be >=0 beyond this point
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Bring all pools, in priority order, up to allometric targets if possible
    ! -----------------------------------------------------------------------------------
    
    priority_loop: do i_pri = 1, n_max_priority
       
       curpri_org(:) = fates_unset_int    ! "current-priority" organ indices

       i = 0
       do ii = 1, num_organs
          
          ! The priority code associated with this organ
          ! Storage has a special hard-coded priority level of 2
          ! Note that it is also implicitly part of step 1

          if( organ_list(ii).eq.store_organ ) then
             priority_code = 2
          else
             if( prt_params%organ_param_id(organ_list(ii)) <1 ) then
                priority_code = -1
             else
                priority_code = int(prt_params%alloc_priority(ipft,prt_params%organ_param_id(organ_list(ii))))
             end if
          end if
       
          
          ! Don't allow allocation to leaves if they are in an "off" status.
          ! (this prevents accidental re-flushing on the day they drop)
          if((leaf_status.eq.leaves_off) .and. (organ_list(ii).eq.leaf_organ)) cycle
          
          
          if( priority_code == i_pri ) then
             deficit_c(ii) = max(0._r8,this%GetDeficit(carbon12_element,organ_list(ii),target_c(ii)))
             i = i + 1
             curpri_org(i) = ii
          end if
       end do

       ! Bring carbon up to target first, this order is required
       ! because we need to know the resulting carbon concentrations
       ! before  we set the allometric targets for the nutrients

       n_curpri_org  = i


       sum_c_demand = 0._r8
       do i=1,n_curpri_org
          i_org = curpri_org(i)
          sum_c_demand = sum_c_demand + deficit_c(i_org)
       end do
       
       sum_c_flux = min(c_gain, sum_c_demand)
       
       ! Transfer carbon into pools if there is any
       if (sum_c_flux>nearzero) then
          do i = 1, n_curpri_org
             
             i_org = curpri_org(i)

             c_flux  =  sum_c_flux*deficit_c(i_org)/sum_c_demand
             
             ! Update the carbon pool
             state_c(i_org)%ptr = state_c(i_org)%ptr + c_flux
             
             ! Update carbon pools deficit
             deficit_c(i_org) = max(0._r8,deficit_c(i_org) - c_flux)
             
             ! Reduce the carbon gain
             c_gain      = c_gain - c_flux
             
          end do
       end if

       ! Determine nutrient demand and make tansfers
       do i = 1, n_curpri_org
          
          i_org = curpri_org(i)

          ! Update the nitrogen deficits
          ! Note that the nitrogen target is tied to the stoichiometry of thegrowing pool only
          target_n = this%GetNutrientTarget(nitrogen_element,organ_list(i_org),stoich_growth_min)
          deficit_n(i_org) = max(0.0_r8, target_n - state_n(i_org)%ptr )
             
          ! Update the phosphorus deficits (which are based off of carbon actual..)
          ! Note that the phsophorus target is tied to the stoichiometry of thegrowing pool only (also)
          target_p = this%GetNutrientTarget(phosphorus_element,organ_list(i_org),stoich_growth_min)
          deficit_p(i_org) = max(0.0_r8, target_p - state_p(i_org)%ptr )

       end do
          
       ! Allocate nutrients at this priority level Nitrogen
       call ProportionalNutrAllocation(state_n, deficit_n, &
            n_gain, nitrogen_element, curpri_org(1:n_curpri_org))
       
       ! Phosphorus
       call ProportionalNutrAllocation(state_p, deficit_p, &
            p_gain, phosphorus_element, curpri_org(1:n_curpri_org))
       

    end do priority_loop


    
    
    return
  end subroutine CNPPrioritizedReplacement

  
  ! =====================================================================================
    
  subroutine CNPStatureGrowth(this,c_gain, n_gain, p_gain, &
                              state_c, state_n, state_p,           &
                              target_c, target_dcdd, cnp_limiter)
    
    
    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain      ! Total daily C gain that remains to be used
    real(r8), intent(inout) :: n_gain      ! Total N available for allocation
                                           ! (new uptake + storage)
    real(r8), intent(inout) :: p_gain      ! Total P available for allocation
                                           ! (new uptake + storage)
    type(parray_type) :: state_c(:)       ! State array for carbon, by organ [kg]
    type(parray_type) :: state_n(:)       ! State array for N, by organ [kg]
    type(parray_type) :: state_p(:)       ! State array for P, by organ [kg]
    real(r8), intent(in) :: target_c(:)
    real(r8), intent(in) :: target_dcdd(:)
    integer, intent(out) :: cnp_limiter
    
    real(r8), pointer :: dbh
    integer           :: ipft
    real(r8)          :: canopy_trim
    real(r8)          :: leaf_status
    real(r8)          :: l2fr
    
    integer  :: i, ii                            ! organ index loops (masked and unmasked)
    integer  :: istep                            ! outer step iteration loop
    real(r8) :: grow_c_from_c                    ! carbon transferred into tissues
    real(r8) :: grow_c_from_n                    ! carbon needed to match N transfers to tissues
    real(r8) :: grow_c_from_p                    ! carbon needed to match P transfers to tissues
    real(r8) :: total_dcostdd                    ! Total carbon transferred to all pools for unit growth
    logical  :: step_pass                        ! flag stating if the integration sub-steps passed checks
    real(r8) :: totalC                           ! total carbon sent to integrator (kg)
    real(r8) :: deltaC                           ! trial value for substep change in carbon (kg)
    real(r8) :: cdeficit                         ! carbon deficit from target
    integer  :: ierr                             ! error flag for allometric growth step
    integer  :: nsteps                           ! number of sub-steps
    real(r8) :: repro_c_frac                     ! Fraction of C allocated to reproduction
                                                 ! at current stature (dbh) [/]
    real(r8) :: sum_c_flux                       ! Sum of the carbon allocated, as reported
                                                 ! by the ODE solver. [kg]
    real(r8) :: np_limit
    real(r8) :: n_match
    real(r8) :: p_match
    real(r8) :: c_flux_adj                       ! Adjustment to total carbon flux during stature growth
                                                 ! intended to correct integration error (kg/kg)
    real(r8) :: c_flux                           ! Carbon flux from the gain pool to an organ (kgC)
    real(r8) :: gr_flux                          ! Growth respiration flux for the current transaction (kgC)
    real(r8) :: c_gstature                       ! Carbon reserved for stature growth (kg)
    real(r8) :: target_n                         ! Target mass of N for a given organ [kg]
    real(r8) :: target_p                         ! Target mass of P for a given organ [kg]
    real(r8) :: sum_n_demand                     ! Total N deficit to overcome after C stature growth [kg]
    real(r8) :: sum_p_demand                     ! Total P deficit to overcome after C stature growth [kg]
    real(r8), dimension(num_organs) :: frac_c    ! Fraction of C going towards each pool
                                                 ! (only used when calculating which species limits)
    real(r8), dimension(num_organs) :: deficit_n ! Deficit to get to target from current    [kg]
    real(r8), dimension(num_organs) :: deficit_p ! Deficit to get to target from current    [kg]
    integer,dimension(num_organs) :: mask_organs ! This works with "state_mask", the list
                                                 ! of organs in the mask
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

    integer, parameter  :: grow_lim_type = 3 ! Dev flag for growth limitation algorithm
                                             ! 1 = tries to calculate equivalent carbon
                                             ! 2 = modification of 1
                                             ! 3 = don't limit, and assume nutrient limitations will prevent calling
                                             !     of this step on the next cycle if they exist
    
    integer, parameter  :: c_limited = 1
    integer, parameter  :: n_limited = 2
    integer, parameter  :: p_limited = 3
    

    leaf_status = this%bc_in(acnp_bc_in_id_lstat)%ival
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    l2fr        = this%bc_inout(acnp_bc_inout_id_l2fr)%rval ! This variable is not updated in this
                                                            ! routine, and is therefore not a pointer
    
    cnp_limiter = 0
    
    ! If any of these resources is essentially tapped out,
    ! then there is no point in performing growth
    ! It also seems impossible that we would be in a leaf-off status
    ! and have enough carbon to grow stature, but its possible that
    ! a plant had a productive last day before the phenology scheme
    ! signaled a drop. If this is the case, we can't grow stature
    ! cause that would force the leaves back on, so just leave.

    
    if( c_gain <= calloc_abs_error  .or. &
         leaf_status.eq.leaves_off  .or. &
         n_gain <= 0.1_r8*calloc_abs_error .or. &
         p_gain <= 0.02_r8*calloc_abs_error ) then
       return
    end if
    
       
    intgr_params(:)              = fates_unset_r8
    intgr_params(intgr_parm_ctrim) = this%bc_in(acnp_bc_in_id_ctrim)%rval
    intgr_params(intgr_parm_pft)   = real(this%bc_in(acnp_bc_in_id_pft)%ival)
    intgr_params(intgr_parm_l2fr)  = this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    
    state_mask(:) = .false.
    mask_organs(:) = fates_unset_int
    
    ! Go through and flag the integrating variables as either pools that
    ! are growing in this iteration, or not.   At this point, if carbon for growth
    ! remains, it means that all pools are up to, or above the target.  If
    ! it is above, that is because of numerical integration errors, or fusion.
    ! In that case, we flag that pool to not be included in stature growth. It will
    ! catch up with the other pools in the next stature growth steps.

    ii = 0
    do i = 1, num_organs

       cdeficit = this%GetDeficit(carbon12_element,organ_list(i),target_c(i))
       
       if ( cdeficit > calloc_abs_error ) then
          ! In this case, we somehow still have carbon to play with,
          ! yet one of the pools is below its current target
          ! gracefully fail
          write(fates_log(),*) 'A carbon pool has reached the stature growth step'
          write(fates_log(),*) 'yet its deficit is too large to integrate '
          write(fates_log(),*) 'organ: ',i
          write(fates_log(),*) 'carbon gain: ',c_gain
          write(fates_log(),*) 'leaves status:', leaf_status
          write(fates_log(),*) cdeficit, target_c(i), state_c(i)%ptr
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
          if (organ_list(i).ne.repro_organ) then
             ii=ii+1
             mask_organs(ii) = i
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
          write(fates_log(),*) 'c state1 : ',state_c(1)%ptr
          write(fates_log(),*) 'c targets: ',target_c(1:num_organs)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    ! fraction of carbon going towards reproduction. reproductive carbon is
    ! just different from the other pools. It is not based on proportionality,
    ! so its mask is set differently.  We (inefficiently) just included
    ! reproduction in the previous loop, but oh well, we over-write now.
    
    if (dbh <= prt_params%dbh_repro_threshold(ipft)) then
       repro_c_frac = prt_params%seed_alloc(ipft)
    else
       repro_c_frac = prt_params%seed_alloc(ipft) + prt_params%seed_alloc_mature(ipft)
    end if
    
    if(repro_c_frac>nearzero)then
       state_mask(repro_id)            = .true.
       ii = ii + 1
       n_mask_organs = ii
       mask_organs(ii) = repro_id
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

    do ii = 1, n_mask_organs
       i = mask_organs(ii)
       total_dcostdd = total_dcostdd + target_dcdd(i)
    end do

    frac_c(:) = 0._r8
    do ii = 1, n_mask_organs
       i = mask_organs(ii)
       frac_c(i) = target_dcdd(i)/total_dcostdd * (1.0_r8 - repro_c_frac)
    end do
    frac_c(repro_id) = repro_c_frac

    if(debug) then
       if ( abs(sum(frac_c,dim=1)-1._r8)>rsnbl_math_prec ) then
          write(fates_log(),*) 'predicted carbon allocation fractions dont sum to 1?'
          write(fates_log(),*) 'frac_c(:):',frac_c
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

       
    ! No mathematical co-limitation of growth
    ! This assumes that limitations will prevent
    ! organs from allowing the growth step to even occur
    ! and thus from an algorithmic level limit growth
    
    c_gstature = c_gain

    
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
         state_array(i) = state_c(i)%ptr
      end do
      
      state_mask(dbh_id)       = .true.
      state_array(dbh_id)      = dbh
      
      totalC                   = c_gstature

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
            
            call CheckIntegratedAllometries(state_array_out(dbh_id),ipft,canopy_trim, l2fr,  &
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
               i = mask_organs(ii)
               sum_c_flux = sum_c_flux + (state_array(i) - state_c(i)%ptr)
            end do
            
            ! This is a correction factor that forces
            ! mass conservation
            c_flux_adj = c_gstature/sum_c_flux
            
            do ii = 1, n_mask_organs
               
               i = mask_organs(ii)
               
               ! Calculate adjusted flux
               c_flux   = (state_array(i) - state_c(i)%ptr)*c_flux_adj
               
               ! update the carbon pool (in all pools flux goes into the first pool)
               state_c(i)%ptr = state_c(i)%ptr + c_flux
               
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
               write(fates_log(),*) 'totalC',totalC
               write(fates_log(),*) 'pft: ',ipft
               write(fates_log(),*) 'dbh: ',dbh
               write(fates_log(),*) 'dCleaf_dd: ',target_dcdd(leaf_id)
               write(fates_log(),*) 'dCfnrt_dd: ',target_dcdd(fnrt_id)
               write(fates_log(),*) 'dCstore_dd: ',target_dcdd(store_id)
               write(fates_log(),*) 'dCsapw_dd: ',target_dcdd(sapw_id)
               write(fates_log(),*) 'dCstruct_dd: ',target_dcdd(struct_id)
               write(fates_log(),*) 'repro c frac: ',repro_c_frac
               dbh_tp1     = state_array_out(dbh_id)
               leafc_tp1   = state_array_out(leaf_id)
               fnrtc_tp1   = state_array_out(fnrt_id)
               sapwc_tp1   = state_array_out(sapw_id)
               storec_tp1  = state_array_out(store_id)
               structc_tp1 = state_array_out(struct_id)
               
               call bleaf(dbh_tp1,ipft,canopy_trim,leaf_c_target_tp1)
               call bfineroot(dbh_tp1,ipft,canopy_trim,l2fr,fnrt_c_target_tp1)
               call bsap_allom(dbh_tp1,ipft,canopy_trim,sapw_area,sapw_c_target_tp1)
               call bagw_allom(dbh_tp1,ipft,agw_c_target_tp1)
               call bbgw_allom(dbh_tp1,ipft,bgw_c_target_tp1)
               call bdead_allom(agw_c_target_tp1,bgw_c_target_tp1, sapw_c_target_tp1, ipft, struct_c_target_tp1)
               call bstore_allom(dbh_tp1,ipft,canopy_trim,store_c_target_tp1)
               
               write(fates_log(),*) 'leaf_c: ',leafc_tp1, leaf_c_target_tp1,leafc_tp1-leaf_c_target_tp1
               write(fates_log(),*) 'fnrt_c: ',fnrtc_tp1, fnrt_c_target_tp1,fnrtc_tp1- fnrt_c_target_tp1
               write(fates_log(),*) 'sapw_c: ',sapwc_tp1, sapw_c_target_tp1 ,sapwc_tp1- sapw_c_target_tp1
               write(fates_log(),*) 'store_c: ',storec_tp1, store_c_target_tp1,storec_tp1- store_c_target_tp1
               write(fates_log(),*) 'struct_c: ',structc_tp1, struct_c_target_tp1,structc_tp1- struct_c_target_tp1
               write(fates_log(),*) 'sapw_c_t0: ',state_c(sapw_id)%ptr, target_c(sapw_id)
               
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

      !!target_n = this%GetNutrientTarget(nitrogen_element,repro_organ,stoich_growth_min)
      !!deficit_n(repro_id) = this%GetDeficit(nitrogen_element,repro_organ,target_n)

      !!target_p = this%GetNutrientTarget(phosphorus_element,repro_organ,stoich_growth_min)
      !!deficit_p(repro_id) = this%GetDeficit(phosphorus_element,repro_organ,target_p)
      
      ! Nitrogen for
      !!call ProportionalNutrAllocation(state_n, deficit_n, n_gain, nitrogen_element,[repro_id])
      
      ! Phosphorus
      !!call ProportionalNutrAllocation(state_p, deficit_p, p_gain, phosphorus_element,[repro_id])

      
      ! -----------------------------------------------------------------------------------
      ! Nutrient Fluxes proportionally to each pool (these should be fully actualized)
      ! (this also removes from the gain pools)
      ! -----------------------------------------------------------------------------------

      sum_n_demand = 0._r8   ! For error checking
      sum_p_demand = 0._r8   ! For error checking
      do ii = 1, n_mask_organs

         i = mask_organs(ii)

         target_n = this%GetNutrientTarget(nitrogen_element,organ_list(i),stoich_growth_min)
         target_p = this%GetNutrientTarget(phosphorus_element,organ_list(i),stoich_growth_min)

         deficit_n(i) = this%GetDeficit(nitrogen_element,organ_list(i),target_n)
         sum_n_demand = sum_n_demand+max(0._r8,deficit_n(i))
         
         deficit_p(i) = this%GetDeficit(phosphorus_element,organ_list(i),target_p)
         sum_p_demand = sum_p_demand+max(0._r8,deficit_p(i))
         
      end do
       
       ! Nitrogen
       call ProportionalNutrAllocation(state_n,deficit_n, & 
            n_gain, nitrogen_element,mask_organs(1:n_mask_organs))
       
       ! Phosphorus
       call ProportionalNutrAllocation(state_p, deficit_p, &
            p_gain, phosphorus_element,mask_organs(1:n_mask_organs))
       
       
    end if if_stature_growth
    
    return
  end subroutine CNPStatureGrowth
  
  ! =====================================================================================

  subroutine CNPAllocateRemainder(this,c_gain, n_gain, p_gain, &
                                  state_c, state_n, state_p, c_efflux, n_efflux, p_efflux)

    class(cnp_allom_prt_vartypes) :: this
    real(r8), intent(inout) :: c_gain
    real(r8), intent(inout) :: n_gain
    real(r8), intent(inout) :: p_gain
    type(parray_type) :: state_c(:)       ! State array for carbon, by organ [kg]
    type(parray_type) :: state_n(:)       ! State array for N, by organ [kg]
    type(parray_type) :: state_p(:)       ! State array for P, by organ [kg]
    real(r8), intent(inout) :: c_efflux
    real(r8), intent(inout) :: n_efflux
    real(r8), intent(inout) :: p_efflux

    integer  :: i
    real(r8), dimension(num_organs) :: deficit_n
    real(r8), dimension(num_organs) :: deficit_p
    real(r8) :: target_n
    real(r8) :: target_p
    real(r8) :: store_c_target   ! Target amount of C in storage including "overflow" [kgC]
    real(r8) :: total_c_flux     ! Total C flux from gains into storage and growth R [kgC]
    real(r8) :: growth_r_flux    ! Growth respiration for filling storage [kgC]
    real(r8) :: store_m_flux     ! Flux into storage [kg]
    integer, dimension(num_organs),parameter :: all_organs = [1,2,3,4,5,6]
    real(r8), pointer :: dbh
    integer           :: ipft
    real(r8)          :: canopy_trim
    

    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival


    ! -----------------------------------------------------------------------------------
    ! If nutrients are still available, then we can bump up the values in the pools
    !  towards the OPTIMAL target values.
    ! -----------------------------------------------------------------------------------

    do i = 1, num_organs
       
       ! Update the nitrogen and phosphorus deficits
       target_n = this%GetNutrientTarget(nitrogen_element,organ_list(i),stoich_max)
       deficit_n(i) = max(0._r8,this%GetDeficit(nitrogen_element,organ_list(i),target_n))
       
       target_p = this%GetNutrientTarget(phosphorus_element,organ_list(i),stoich_max)
       deficit_p(i) = max(0._r8,this%GetDeficit(phosphorus_element,organ_list(i),target_p))
          
    end do

    ! -----------------------------------------------------------------------------------
    ! Nutrient Fluxes proportionally to each pool (these should be fully actualized)
    ! (this also removes from the gain pools)
    ! -----------------------------------------------------------------------------------
    
    ! Nitrogen
    call ProportionalNutrAllocation(state_n(1:num_organs), & 
         deficit_n(1:num_organs), & 
         n_gain, nitrogen_element, all_organs)
    
    ! Phosphorus
    call ProportionalNutrAllocation(state_p(1:num_organs), &
         deficit_p(1:num_organs), &
         p_gain, phosphorus_element, all_organs)


    ! Optional hypothesis (
    ! If any N or P is still hanging around, put it in storage

    !state_n(store_id)%ptr = state_n(store_id)%ptr + n_gain
    !state_p(store_id)%ptr = state_p(store_id)%ptr + p_gain

    !n_gain = 0._r8
    !p_gain = 0._r8
    
    
    ! -----------------------------------------------------------------------------------
    ! If carbon is still available, lets cram some into storage overflow
    ! We will do this last, because we wanted the non-overflow storage
    ! value to draw minimum and optimal nutrient fluxes
    ! -----------------------------------------------------------------------------------

    if(c_gain>calloc_abs_error) then


!       select(c_overlow_method)
!    case(store_c_overflow)
       
!    case(efflux_c_overflow)
       
!    case(burn_c_overflow)
       
!    end if


    
       if(force_store_c_overflow)then

          total_c_flux = c_gain
       else
          
          ! Update carbon based allometric targets
          call bstore_allom(dbh,ipft,canopy_trim, store_c_target)
          
          ! Estimate the overflow
          store_c_target = store_c_target * (1.0_r8 + store_overflow_frac)
          
          total_c_flux = min(c_gain,max(0.0, (store_c_target - state_c(store_id)%ptr)))
          
       end if
       ! Transfer excess carbon into storage overflow
       state_c(store_id)%ptr = state_c(store_id)%ptr + total_c_flux
       c_gain              = c_gain - total_c_flux
       
          

    end if

    

    ! Figure out what to do with excess carbon and nutrients
    ! 1) excude through roots cap at 0 to flush out imprecisions
    ! -----------------------------------------------------------------------------------

    c_efflux = max(0.0_r8,c_gain)
    n_efflux = max(0.0_r8,n_gain)
    p_efflux = max(0.0_r8,p_gain)


    c_gain = 0.0_r8
    n_gain = 0.0_r8
    p_gain = 0.0_r8

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
    real(r8)         :: sapw_area
    real(r8)         :: leaf_c_target,fnrt_c_target
    real(r8)         :: sapw_c_target,agw_c_target
    real(r8)         :: bgw_c_target,struct_c_target
    
    dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    canopy_trim = this%bc_in(acnp_bc_in_id_ctrim)%rval
    ipft        = this%bc_in(acnp_bc_in_id_pft)%ival
    i_cvar      = prt_global%sp_organ_map(organ_id,carbon12_element)
    l2fr        = this%bc_inout(acnp_bc_inout_id_l2fr)%rval
    
    ! Storage of nutrients are assumed to have different compartments than
    ! for carbon, and thus their targets are not associated with a tissue
    ! but is more represented as a fraction of the maximum amount of nutrient
    ! that can be bound in non-reproductive tissues
    
    if(organ_id == store_organ) then

       call bleaf(dbh,ipft,canopy_trim,leaf_c_target)
       call bfineroot(dbh,ipft,canopy_trim,l2fr,fnrt_c_target)
       call bsap_allom(dbh,ipft,canopy_trim,sapw_area,sapw_c_target)
       call bagw_allom(dbh,ipft,agw_c_target)
       call bbgw_allom(dbh,ipft,bgw_c_target)
       call bdead_allom(agw_c_target,bgw_c_target, sapw_c_target, ipft, struct_c_target)

       ! Target for storage is a fraction of the sum target of all
       ! non-reproductive organs

       if( element_id == nitrogen_element) then
          
          target_m = StorageNutrientTarget(ipft, element_id, &
               leaf_c_target*prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(leaf_organ)), &
               fnrt_c_target*prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(fnrt_organ)), &
               sapw_c_target*prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(sapw_organ)), & 
               struct_c_target*prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(struct_organ)))
       else
          
          target_m = StorageNutrientTarget(ipft, element_id, &
               leaf_c_target*prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(leaf_organ)), &
               fnrt_c_target*prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(fnrt_organ)), &
               sapw_c_target*prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(sapw_organ)), & 
               struct_c_target*prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(struct_organ)))
          
       end if

       ! Hard-code the growth minimum storage stoichiometry to 75% of maximum
       if( stoich_mode == stoich_growth_min ) then
          target_m = target_m*0.25_r8
       end if

    elseif(organ_id == repro_organ) then

       target_c = this%variables(i_cvar)%val(1)
       if( element_id == nitrogen_element) then
          target_m = target_c * prt_params%nitr_recr_stoich(ipft)
       else
          target_m = target_c * prt_params%phos_recr_stoich(ipft)
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
          if( element_id == nitrogen_element) then
             target_m = target_c * prt_params%nitr_stoich_p2(ipft,prt_params%organ_param_id(organ_id))
          else
             target_m = target_c * prt_params%phos_stoich_p2(ipft,prt_params%organ_param_id(organ_id))
          end if
       else
          write(fates_log(),*) 'invalid stoichiometry mode specified while getting'
          write(fates_log(),*) 'nutrient targets'
          write(fates_log(),*) 'stoich_mode: ',stoich_mode
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if
    
    return
  end function GetNutrientTargetCNP


  
  ! =====================================================================================
  
  subroutine ProportionalNutrAllocation(state_m, deficit_m, gain_m, element_id, list)

    ! -----------------------------------------------------------------------------------
    ! This routine allocates nutrients to a set of organs based on proportional
    ! need.  It is assumed that the input arrays are not sparse, and the fluxes
    ! are based purely off of there deficit from some ideal state.
    ! Note: this may or may not be called inside some preferential organ filter.
    ! -----------------------------------------------------------------------------------

    type(parray_type) :: state_m(:)        ! Current mass of nutrient
                                           ! of arbitrary species
                                           ! over some arbitrary set of organs
    real(r8),intent(inout) :: deficit_m(:) ! Nutrient mass deficit of species
                                           ! over set of organs
    integer, intent(in)    :: list(:)      ! List of indices if sparse
    real(r8),intent(inout) :: gain_m       ! Total nutrient mass gain to
                                           ! work with
    integer,intent(in) :: element_id       ! Element global index (for debugging)

    ! locals
    integer :: num_organs
    integer :: i,ii
    real(r8) :: flux
    real(r8) :: sum_deficit
    real(r8) :: sum_flux

    num_organs = size(list,dim=1)
       
    sum_deficit = 0._r8
    do ii = 1, num_organs
       i = list(ii)
       sum_deficit = sum_deficit + max(0._r8,deficit_m(i))
    end do
    
    if (sum_deficit>nearzero) then
       
       sum_flux = min(gain_m, sum_deficit)
       
       do ii = 1, num_organs
          i = list(ii)
          
          flux        = sum_flux * max(0._r8,deficit_m(i))/sum_deficit
          state_m(i)%ptr  = state_m(i)%ptr + flux
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
      real(r8) :: total_dcdd_target     ! target total (not reproductive) biomass derivative wrt d, (kgC/cm)
      real(r8) :: repro_fraction        ! fraction of carbon balance directed towards reproduction (kgC/kgC)
      real(r8) :: total_dcostdd         ! carbon cost for non-reproductive pools per unit increment of dbh
      

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

        canopy_trim = intgr_params(intgr_parm_ctrim)
        ipft        = int(intgr_params(intgr_parm_pft))
        l2fr        = intgr_params(intgr_parm_l2fr)

        call bleaf(dbh,ipft,canopy_trim,leaf_c_target,leaf_dcdd_target)
        call bfineroot(dbh,ipft,canopy_trim,l2fr,fnrt_c_target,fnrt_dcdd_target)
        call bsap_allom(dbh,ipft,canopy_trim,sapw_area,sapw_c_target,sapw_dcdd_target)
        call bagw_allom(dbh,ipft,agw_c_target,agw_dcdd_target)
        call bbgw_allom(dbh,ipft,bgw_c_target,bgw_dcdd_target)
        call bdead_allom(agw_c_target,bgw_c_target, sapw_c_target, ipft, struct_c_target, &
                         agw_dcdd_target, bgw_dcdd_target, sapw_dcdd_target, struct_dcdd_target)
        call bstore_allom(dbh,ipft,canopy_trim,store_c_target,store_dcdd_target)

        if (mask_repro) then
           ! fraction of carbon going towards reproduction
           if (dbh <= prt_params%dbh_repro_threshold(ipft)) then
              repro_fraction = prt_params%seed_alloc(ipft)
           else
              repro_fraction = prt_params%seed_alloc(ipft) + prt_params%seed_alloc_mature(ipft)
           end if
        else
           repro_fraction = 0._r8
        end if

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

   ! ====================================================================================

   subroutine TargetAllometryCheck(bleaf,bfroot,bsap,bstore,bdead, &
                                   bt_leaf,bt_froot,bt_sap,bt_store,bt_dead, &
                                   grow_leaf,grow_froot,grow_sapw,grow_store)

     ! Arguments
     real(r8),intent(in) :: bleaf   !actual
     real(r8),intent(in) :: bfroot
     real(r8),intent(in) :: bsap
     real(r8),intent(in) :: bstore
     real(r8),intent(in) :: bdead
     real(r8),intent(in) :: bt_leaf   !target
     real(r8),intent(in) :: bt_froot
     real(r8),intent(in) :: bt_sap
     real(r8),intent(in) :: bt_store
     real(r8),intent(in) :: bt_dead
     logical,intent(out) :: grow_leaf  !growth flag
     logical,intent(out) :: grow_froot
     logical,intent(out) :: grow_sapw
     logical,intent(out) :: grow_store

     if( (bt_leaf - bleaf)>calloc_abs_error) then
        write(fates_log(),*) 'leaves are not on-allometry at the growth step'
        write(fates_log(),*) 'exiting',bleaf,bt_leaf
        call endrun(msg=errMsg(sourcefile, __LINE__))
     elseif( (bleaf - bt_leaf)>calloc_abs_error) then
        ! leaf is above allometry, ignore
        grow_leaf = .false.
     else
        grow_leaf = .true.
     end if

     if( (bt_froot - bfroot)>calloc_abs_error) then
        write(fates_log(),*) 'fineroots are not on-allometry at the growth step'
        write(fates_log(),*) 'exiting',bfroot, bt_froot
        call endrun(msg=errMsg(sourcefile, __LINE__))
     elseif( ( bfroot-bt_froot)>calloc_abs_error ) then
        grow_froot = .false.
     else
        grow_froot = .true.
     end if

     if( (bt_sap - bsap)>calloc_abs_error) then
        write(fates_log(),*) 'sapwood is not on-allometry at the growth step'
        write(fates_log(),*) 'exiting',bsap, bt_sap
        call endrun(msg=errMsg(sourcefile, __LINE__))
     elseif( ( bsap-bt_sap)>calloc_abs_error ) then
        grow_sapw = .false.
     else
        grow_sapw = .true.
     end if

     if( (bt_store - bstore)>calloc_abs_error) then
        write(fates_log(),*) 'storage is not on-allometry at the growth step'
        write(fates_log(),*) 'exiting',bstore,bt_store
        call endrun(msg=errMsg(sourcefile, __LINE__))
     elseif( ( bstore-bt_store)>calloc_abs_error ) then
        grow_store = .false.
     else
        grow_store = .true.
     end if

     if( (bt_dead - bdead)>calloc_abs_error) then
        write(fates_log(),*) 'structure not on-allometry at the growth step'
        write(fates_log(),*) 'exiting',bdead,bt_dead
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if
   end subroutine TargetAllometryCheck

   
   ! =====================================================================================
   
   
   subroutine StorageRegulator(this,element_id,regulate_type,c_scalar)

     ! -----------------------------------------------------------------------------------
     ! This function evaluates the storage of either N or P, and returns
     ! a scalar that is used to regulate fine-root biomass in some way. Depending
     ! on the type of method, this may either be an absolute scalar on biomass
     ! or it may be a rate of change on that scalar. 
     ! -----------------------------------------------------------------------------------


     ! Arguments (in)
     class(cnp_allom_prt_vartypes) :: this
     integer,intent(in)  :: element_id              ! element id consistent with parteh/PRTGenericMod.F90
     integer,intent(in)  :: regulate_type


     ! Arguments (out)
     real(r8) :: c_scalar

     ! Locals
     real(r8) :: store_frac  ! Current nutrient storage relative to max
     real(r8) :: store_max   ! Maximum nutrient storable by plant
     real(r8) :: store_c     ! Current storage carbon
     real(r8) :: store_c_max ! Current maximum storage carbon
     integer  :: icode       ! real variable checking code
     real(r8) :: store_x
     integer  :: i_var
     real(r8), parameter :: c_eq_offset = 0.95 ! This shifts the center-point
                                               ! of the N:C or P:C storage equlibrium
                                               ! by multiplying the C term. If its less than 1 it
                                               ! shifts left and great than one it shifts right.
                                               ! It should shift left to help mitigate wasted N and P
                                               ! storage overflow
     
     ! For N/C logistic
     real(r8) :: logi_k          ! logistic function k
     real(r8) :: store_x0        ! storage fraction inflection point
     real(r8) :: logi_min        ! minimum cn_scalar for logistic
     real(r8) :: l2fr_delta_max
     real(r8) :: l2fr_delta_min
     real(r8) :: l2fr_actual
     real(r8) :: leaf_c, fnrt_c
     real(r8) :: leaf_c_target
     real(r8) :: log_nc_frac
     real(r8) :: store_c_frac
     real(r8) :: c_gain
     real(r8) :: fnrt_c_target
     
     ! This is the storage fraction where downregulation starts if using
     ! a linear function
     real(r8), parameter :: store_frac0 = 0.85_r8

     real(r8), parameter :: c_max = 1.0_r8   ! Maximum allowable result of the function
     real(r8), parameter :: c_min = 0.0_r8   ! Minimum allowable result of the function

     ! This fraction governs
     ! how much carbon from daily gains + storage overflow, is allowed to
     ! be spent on growing out roots. This inludes getting roots
     ! back on allometry before growing out
     
     real(r8), parameter :: max_l2fr_cgain_frac = 0.5_r8
     
     associate(dbh         => this%bc_inout(acnp_bc_inout_id_dbh)%rval, & 
               canopy_trim => this%bc_in(acnp_bc_in_id_ctrim)%rval, &
               ipft        => this%bc_in(acnp_bc_in_id_pft)%ival, & 
               l2fr        => this%bc_inout(acnp_bc_inout_id_l2fr)%rval)
              
       store_max = this%GetNutrientTarget(element_id,store_organ,stoich_max)

       ! Storage fractions could more than the target, depending on the
       ! hypothesis and functions involved, but should typically be 0-1
       ! The cap of 2 is for numerics and preventing weird math
       store_frac = max(0.01_r8,min(1.0_r8,this%GetState(store_organ, element_id)/store_max))
       
       if(regulate_type == regulate_linear) then

          c_scalar = min(c_max,max(c_min,1.0 - (store_frac - store_frac0)/(1.0_r8-store_frac0)))

       elseif(regulate_type == regulate_logi) then

          logi_k   = 30.0_r8
          store_x0 = 0.7_r8
          logi_min = 0.0_r8
          
          ! In this method, we define the c_scalar term
          ! with a logistic function that goes to 1 (full need)
          ! as the plant's nutrien storage hits a low threshold
          ! and goes to 0, no demand, as the plant's nutrient
          ! storage approaches it's maximum holding capacity

          c_scalar = max(c_min,min(c_max,logi_min + (1.0_r8-logi_min)/(1.0_r8 + exp(logi_k*(store_frac-store_x0)))))

          !call check_var_real(c_scalar,'c_scalar',icode)
          !if (icode .ne. 0) then
          !   write(fates_log(),*) 'c_scalar is invalid, element: ',element_id
          !   write(fates_log(),*) 'ending'
          !   call endrun(msg=errMsg(sourcefile, __LINE__))
          !endif

       elseif(regulate_type == regulate_CN_logi .or. regulate_type == regulate_CN_ema) then

          logi_k   = 1.0_r8
          store_x0 = 0.0_r8
          logi_min = 0.0_r8

          ! Only update L2FR if some leaves are out
          call bleaf(dbh,ipft,canopy_trim,leaf_c_target)
          if(this%GetState(leaf_organ, carbon12_element)/leaf_c_target>0.01_r8) then
             store_max = this%GetNutrientTarget(element_id,store_organ,stoich_max)
             
             ! Storage fractions could be more than the target, depending on the
             ! hypothesis and functions involved, but should typically be 0-1
             ! The cap of 2 is for numerics and preventing weird math
             
             store_frac = max(0.01_r8,min(2.0_r8,this%GetState(store_organ, element_id)/store_max))
             
             
             call bstore_allom(dbh,ipft,canopy_trim,store_c_max)
             
             ! Since we don't dump storage carbon
             ! these stores can actually get pretty large, so the cap of 10x is numerically
             ! feasable, and should also minimize stress on the logistic function
             store_c_frac = max(0.01_r8,min(10.0_r8,c_eq_offset*(this%GetState(store_organ, carbon12_element)/store_c_max)))
             
             
             ! -----------------------------------------------------------------------------
             ! To decide the upper limit on expanding root growth, we perform a carbon
             ! balance. Note that if we are growing roots out more, than we have proportionaly
             ! more C compared to other resources. Specifically, we want to limit root growth
             ! such that allocation to roots can't exceed a certain fraction of the daily
             ! available carbon.  This fraction is "max_l2fr_cgain_frac".
             ! Additional notes.  When calculating the "allocation to roots", we consider
             ! both the carbon necessary to get the roots "on allometry" plux the carbon
             ! necessary to expand them.
             !
             !  (l2fr_delta_max*target_fnrt_c + target_fnrt_c-actual_fnrt_c )/c_gain
             !               < max_l2fr_cgain_frac
             !
             ! or
             !
             ! l2fr_delta_max*target_fnrt_c < max_l2fr_cgain_frac * (c_gain -
             !                                                       (target_fnrt_c-actual_fnrt_c) -
             !                                                       (target_leaf_c-actual_leaf_c))
             !
             ! ------------------------------------------------------------------------------
             
             call bfineroot(dbh,ipft,canopy_trim, l2fr, fnrt_c_target)
             call bleaf(dbh,ipft,canopy_trim,leaf_c_target)
             
             ! If there is overflow storage, add this to the gain
             c_gain = this%bc_in(acnp_bc_in_id_netdc)%rval + &
                  max(0._r8,this%GetState(store_organ, carbon12_element)-store_c_max)

             fnrt_c = this%GetState(fnrt_organ, carbon12_element)
             leaf_c = this%GetState(leaf_organ, carbon12_element)
             
             l2fr_delta_max = max_l2fr_cgain_frac / fnrt_c_target * &
                  (c_gain - (fnrt_c_target-fnrt_c) - (leaf_c_target-leaf_c)) 

             ! This value could be negative if there is no gain, or less gain
             ! than what can replace leaf/root, just ensure the multiplier is GT 1
             
             l2fr_delta_max = max(1._r8,l2fr_delta_max)
             

             ! Second constraint, folding timescale
             ! 2.0 = l2fr_delta_max^frnt_adapt_tscl
             l2fr_delta_max = min(l2fr_delta_max, 2._r8**(1._r8/fnrt_adapt_tscl))
             
             
             l2fr_actual =  this%GetState(fnrt_organ, carbon12_element) / &
                  this%GetState(leaf_organ, carbon12_element)

             ! Constrain change in l2fr minimum to be no more than what is lost 
             ! in turnover for a day
             l2fr_delta_min = 1._r8-(years_per_day / prt_params%root_long(ipft))

             ! Second constraint, folding timescale
             l2fr_delta_min = max(l2fr_delta_min, 0.5_r8**(1._r8/fnrt_adapt_tscl))

             
             log_nc_frac = log( store_frac / store_c_frac )
             
             ! This is a logistic between -1 and 1
             c_scalar = 2._r8*max(0._r8, &
                  min(1._r8,logi_min + (1._r8-logi_min)/(1._r8 + exp(logi_k*(log_nc_frac-store_x0)))))-1.0_r8
             
             if(c_scalar>0.0_r8)then
                c_scalar = 1._r8 + c_scalar*(l2fr_delta_max-1._r8)
             else
                c_scalar = 1._r8 + c_scalar*(1._r8-l2fr_delta_min)
             end if
             
          else
             c_scalar = 1._r8
          end if
          
          
       elseif(regulate_type == regulate_CN_dfdd) then

          store_c = this%GetState(store_organ, carbon12_element)
          call bstore_allom(dbh,ipft,canopy_trim,store_c_max)

          ! Fraction of N per fraction of C
          ! If this is greater than 1, then we have more N in storage than
          ! we have C, so we downregulate. If this is less than 1, then
          ! we have less N in storage than we have C, so up-regulate

         
          
          store_frac = max(0.01_r8,store_frac) / max(0.01_r8,c_eq_offset*(store_c/store_c_max))

          c_scalar = 1._r8 - 0.02_r8*log(store_frac)

          !print*,element_id,store_frac,c_scalar
          
       end if

     end associate


   end subroutine StorageRegulator



end module PRTAllometricCNPMod
