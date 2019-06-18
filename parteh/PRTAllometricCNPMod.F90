module PRTAllometricCNPMod

   ! ------------------------------------------------------------------------------------
   !
   ! This module contains all of the specific functions and types for
   ! Plant Allocation and Reactive Transport Extensible Hypotheses (PARTEH)
   ! Carbon-Nitrogen-Phosphorus Prioritized Allometric Allocations
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

  use PRTGenericMod , only  : leaf_organ
  use PRTGenericMod , only  : fnrt_organ
  use PRTGenericMod , only  : sapw_organ
  use PRTGenericMod , only  : store_organ
  use PRTGenericMod , only  : repro_organ
  use PRTGenericMod , only  : struct_organ
  use PRTGenericMod , only  : all_organs
  use PRTGenericMod , only  : prt_cnp_flex_allom_hyp

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
  use FatesIntegratorsMod , only : RKF45
  use FatesIntegratorsMod , only : Euler
  use EDPftvarcon         , only : EDPftvarcon_inst
  use FatesConstantsMod   , only : calloc_abs_error
  use FatesConstantsMod   , only : nearzero
  use FatesConstantsMod   , only : itrue
  use FatesConstantsMod   , only : years_per_day

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


  ! This is the ordered list of organs used in this module
  ! -------------------------------------------------------------------------------------

  integer, parameter                        :: num_organs = 6
  integer, parameter, dimension(num_organs) :: organ_list = &
       [leaf_organ, fnrt_organ, sapw_organ, store_organ, repro_organ, struct_organ]
  

  ! -------------------------------------------------------------------------------------
  ! Set up the array used in the integrator, these are PRIVATE
  ! Total number of variables passed to the integrator
  ! this is all carbon pools, plus dbh and growth respiration
  ! *Note the absence of N+P. N+P availability will effect
  ! how much stature growth occurs, but it effects the integration
  ! bounds, instead of being integrated as prognostic variables
  ! themselves.
  ! -------------------------------------------------------------------------------------
  
  integer, parameter :: intgr_leaf_c_id   = 1    ! matches the main state order
  integer, parameter :: intgr_fnrt_c_id   = 2    ! ""
  integer, parameter :: intgr_sapw_c_id   = 3
  integer, parameter :: intgr_store_c_id  = 4
  integer, parameter :: intgr_repro_c_id  = 5
  integer, parameter :: intgr_struct_c_id = 6
  integer, parameter :: intgr_leaf_gr_id   = 7   ! matches the main state order
  integer, parameter :: intgr_fnrt_gr_id   = 8
  integer, parameter :: intgr_sapw_gr_id   = 9
  integer, parameter :: intgr_store_gr_id  = 10
  integer, parameter :: intgr_repro_gr_id  = 11
  integer, parameter :: intgr_struct_gr_id = 12
  integer, parameter :: intgr_dbh_id       = 13
  integer, parameter :: num_intgr_vars     = 13

  
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
  integer, public, parameter :: acnp_bc_inout_id_netdc      = 3  ! Index for the net daily C input BC
  integer, public, parameter :: acnp_bc_inout_id_netdn      = 4  ! Index for the net daily N input BC
  integer, public, parameter :: acnp_bc_inout_id_netdp      = 5  ! Index for the net daily P input BC
  integer, public, parameter :: num_bc_inout                = 5

  ! -------------------------------------------------------------------------------------
  ! Input only Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: acnp_bc_in_id_pft     = 1     ! Index for the PFT input BC
  integer, public, parameter :: acnp_bc_in_id_ctrim   = 2     ! Index for the canopy trim function
  integer, public, parameter :: acnp_bc_in_id_leafon  = 3     ! phenology status logical
                                                              ! 0=leaf off, 1=leaf on
 
  integer, parameter         :: num_bc_in             = 3

  ! -------------------------------------------------------------------------------------
  ! Output Boundary Indices (These are public)
  ! -------------------------------------------------------------------------------------
  
  integer, public, parameter :: acnp_bc_out_id_rootcexude = 1  ! Daily exudation of C  [kg]
  integer, public, parameter :: acnp_bc_out_id_rootnexude = 2  ! Daily exudation of N  [kg]
  integer, public, parameter :: acnp_bc_out_id_rootpexude = 3  ! Daily exudation of P  [kg]
  integer, public, parameter :: acnp_bc_out_id_growresp   = 4
  integer, parameter         :: num_bc_out                = 4  ! Total number of 


  ! -------------------------------------------------------------------------------------
  ! Define the size of the coorindate vector.  For this hypothesis, there is only
  ! one pool per each species x organ combination, except for leaves (WHICH HAVE AGE)
  ! -------------------------------------------------------------------------------------
  integer, parameter :: icd = 1


  ! This is the maximum number of leaf age pools  (used for allocating scratch space)
  integer, parameter         :: max_nleafage  = 10


  ! -------------------------------------------------------------------------------------
  ! This is the core type that holds this specific
  ! plant reactive transport (PRT) module
  ! -------------------------------------------------------------------------------------


  type, public, extends(prt_vartypes) :: cnp_allom_prt_vartypes
     
  contains
     
     procedure :: DailyPRT     => DailyPRTAllometricCNP
     procedure :: FastPRT      => FastPRTAllometricCNP
     
  end type cnp_allom_prt_vartypes

   
   ! ------------------------------------------------------------------------------------
   !
   ! This next class is an extention of the base instance that maps state variables
   !      to the outside model.
   !
   ! ------------------------------------------------------------------------------------
   
   character(len=*), parameter, private :: sourcefile = __FILE__

   ! This is the instance of the mapping table and variable definitions
   ! this is only allocated once per node

   class(prt_global_type), public, target, allocatable :: prt_global_acnp

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

     nleafage = size(EDPftvarcon_inst%leaf_long,dim=2)

     if(nleafage>max_nleafage) then
        write(fates_log(),*) 'The allometric carbon PARTEH hypothesis'
        write(fates_log(),*) 'sets a maximum number of leaf age classes'
        write(fates_log(),*) 'used for scratch space. The model wants'
        write(fates_log(),*) 'exceed that. Simply increase max_nleafage'
        write(fates_log(),*) 'found in parteh/PRTAllometricCarbonMod.F90'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     ! These mappings help define and classify the different variables in the global sense

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

    real(r8),pointer :: dbh              ! Diameter at breast height [cm]
    real(r8),pointer :: carbon_gain      ! Daily carbon balance for this cohort [kgC]
    real(r8),pointer :: nitrogen_gain    ! Daily nitrogen uptake through fine-roots [kgN]
    real(r8),pointer :: phosphorus_gain ! Daily phosphorus uptake through fine-roots [kgN]
    real(r8),pointer :: maint_r_deficit  ! Current maintenance respiration deficit [kgC]


    real(r8) :: canopy_trim              ! The canopy trimming function [0-1]

    real(r8) :: grow_c_from_c            ! carbon transferred into tissues
    real(r8) :: grow_c_from_n            ! carbon needed to match N transfers to tissues
    real(r8) :: grow_c_from_p            ! carbon needed to match P transfers to tissues

    real(r8) :: agw_c_target             ! Target for above-ground carbon biomass (not a state)
    real(r8) :: bgw_c_target             ! Target for below-ground carbon biomass (not a state)
    real(r8) :: agw_dcdd_target          ! Target for above-ground carbon biomass derivative (not a state)
    real(r8) :: bgw_dcdd_target          ! Target for below-ground carbon biomass derivative (not a state)

    real(r8) :: total_dcostdd            ! Total carbon transferred to all pools for unit growth
    real(r8) :: struct_c_frac            ! Approx frac of C sent to structure during stature growth 
    real(r8) :: leaf_c_frac              ! Approx frac of C sent to leaves during stature growth 
    real(r8) :: fnrt_c_frac              ! Approx frac of C sent to fineroots during stature growth 
    real(r8) :: sapw_c_frac              ! Approx frac of C sent to sapwood during stature growth 
    real(r8) :: store_c_frac             ! Approx frac of C sent to storage during stature growth 
    real(r8) :: repro_c_frac             ! Approx frac of C sent to reproduction during stature growth 

    real(r8) :: fnrt_c_flux              ! Unadjusted C flux into fineroots during stature growth (kg)
    real(r8) :: leaf_c_flux              ! Unadjusted C Flux into leaves during stature growth (kg)
    real(r8) :: sapw_c_flux              ! Unadjusted C flux into sapwood during stature growth (kg)
    real(r8) :: struct_c_flux            ! Unadjusted C flux into structure during stature growth (kg)
    real(r8) :: store_c_flux             ! Unadjusted C flux into storage during stature growth (kg)
    real(r8) :: repro_c_flux             ! Unadjusted C flux into reproduction during stature growth (kg)
    real(r8) :: growth_resp_flux         ! Unadjusted C flux into growth respiration (kg)

    real(r8) :: maint_r_def_flux         ! Flux into maintenance respiration during priority 1 allocation
    
    real(r8) :: fnrt_dcdd_target         ! target derivative of fineroot carbon (kg C) / cm
    real(r8) :: leaf_dcdd_target         ! target derivative of leaf carbon (kg C) / cm
    real(r8) :: sapw_dcdd_target         ! target derivative of sapwood carbon (kg C) / cm
    real(r8) :: struct_dcdd_target       ! target derivative of structure carbon (kg C) / cm
    real(r8) :: store_dcdd_target        ! target derivative of storage carbon (kg C) / cm

    

    real(r8) :: carbon_gain0             ! Carbon gain before any allocations (input) kgC
    real(r8) :: maint_r_deficit0         ! Accumulated maintenance respiration deficit 
                                         ! prior to allocations (kgC)
    real(r8) :: nitrogen_error           ! Total nitrogen imbalance kg
    real(r8) :: rel_nitrogen_error       ! Relative nitrogen imbalance kg
    real(r8) :: nitrogen_gain0           ! Nitrogen gain before any allocations
    real(r8) :: phosphorus_error        ! Total phosphorus imbalance (kg)
    real(r8) :: rel_phosphorus_error    ! Relative phosphorus imbalance (kg)
    real(r8) :: phosphorus_gain0        ! Total phosphorus gains before allocation (kg)


    real(r8) :: c_flux_adj               ! Adjustment to total carbon flux during stature growth
                                         ! intended to correct integration error (kg/kg)

    real(r8) :: carbon_gain_flux         ! Flux of carbon out of the carbon gains pool (kg)
    real(r8) :: carbon_gstature          ! Carbon reserved for stature growth (kg)

    real(r8) :: p_m                      ! Maintenance prioritization scalar 
                                         ! 1.0 means maintenance deficit repayments
                                         ! are completely prioritized over priority 1 tissues
                                         ! 0.5 is equal payment, and 0.0 is deprioritized

    integer  :: priority_code            ! Index for priority level of each organ
    real(r8) :: sum_c_demand             ! Carbon demanded to bring tissues up to allometry (kg) 
    real(r8) :: sum_n_demand             ! The nitrogen demand of all pools for given priority level (kg)
    real(r8) :: sum_p_demand             ! The phosphorus demand of all pools for given priority level (kg) 
    real(r8) :: mo_sum_c_flux            ! The total flux going to both organs and maint resp at priority 1 (kg)
    real(r8) :: mo_sum_c_frac            ! The fraction of demanded carbon during priority 1, to both organs and
                                         ! maintenance respiration that can be paid (kg/kg)
    real(r8) :: sum_c_flux               ! The flux to bring tissues up to allometry (kg)
    real(r8) :: sum_n_flux               ! The flux of nitrogen ""  (kg)
    real(r8) :: sum_p_flux               ! The flux of phosphorus "" (Kg)  
    real(r8) :: c_flux                   ! carbon flux into an arbitrary pool (kg)
    real(r8) :: gr_flux                  ! carbon flux to fulfill growth respiration of an arbitrary pool (kg)
    real(r8) :: n_flux                   ! nitrogen flux into  an arbitrary pool (kg)
    real(r8) :: p_flux                   ! phosphorus flux into an arbitrary pool (kg)
    real(r8) :: redist_c_flux            ! The redistributed flux of carbon during priority 1 transfer (kg)
    real(r8) :: store_c_transferable     ! The amount of storage carbon that can be transferred 
                                         ! during priority 1 allocation (kg)
    real(r8) :: sapw_area                ! dummy var for sapwood area (m2)

    integer  :: leaves_on                ! phenological status of plant (2=leaves on, 1=inactive)
    integer  :: ipft                     ! Plant Functional Type index
                                         ! integrator variables
    logical  :: step_pass                ! flag stating if the integration sub-steps passed checks
    real(r8) :: totalC                   ! total carbon sent to integrator (kg)
    real(r8) :: deltaC                   ! trial value for substep change in carbon (kg)
    real(r8) :: cdeficit                 ! carbon deficit from target
    integer  :: ierr                     ! error flag for allometric growth step
    integer  :: nsteps                   ! number of sub-steps
    integer  :: istep                    ! current substep index
    real(r8) :: hite_out                 ! dummy height variable

    integer  :: nleafage                 ! Number of ACTUAL leaf age bins
    integer  :: i                        ! Generic loop counter (mostly used for organ counting)
    integer  :: ii                       ! Generic loop counter (mostly used for organ counting)
    integer  :: i_age                    ! Generic loop counter for leaf age bins
    integer  :: i_var                    ! Index for generic variable
    integer  :: i_cvar                   ! Index for a carbon variable
    integer  :: i_nvar                   ! Index for a nitrogen variable
    integer  :: i_pvar                   ! Index for a phosphorus variable
    integer  :: i_gorgan                 ! Index for organ using the global (generic)
                                         ! indexing system
    integer  :: i_priority               ! Index for the current priority level


    real(r8),dimension(max_nleafage) :: leaf_c0
    real(r8),dimension(max_nleafage) :: leaf_n0
    real(r8),dimension(max_nleafage) :: leaf_p0

    real(r8) :: leaf_age_flux_frac    ! Fraction of leaf mass in each age bin that is transferred
                                      ! to the next
                                      ! Initial value of carbon used to determine net flux
    real(r8) :: fnrt_c0               ! during this routine
    real(r8) :: sapw_c0               ! ""   
    real(r8) :: store_c0              ! ""
    real(r8) :: repro_c0              ! ""
    real(r8) :: struct_c0             ! ""
    real(r8) :: fnrt_n0               ! ""
    real(r8) :: sapw_n0               ! ""   
    real(r8) :: store_n0              ! ""
    real(r8) :: repro_n0              ! ""
    real(r8) :: struct_n0             ! ""
    real(r8) :: fnrt_p0               ! ""
    real(r8) :: sapw_p0               ! ""   
    real(r8) :: store_p0              ! ""
    real(r8) :: repro_p0              ! ""
    real(r8) :: struct_p0             ! ""


    real(r8), dimension(num_vars)  :: v_target             ! Target quantity for each state variable [kg]
    real(r8), dimension(num_vars)  :: v_demand             ! Demand to get to target from current    [kg]
    

    integer  :: num_organs_curpri                          ! Number of pools at current priority level


    real(r8), dimension(num_organs) :: r_g                 ! Unit growth respriation for each pool [kg/kg]
    real(r8), dimension(num_organs) :: growth_resp         ! Total growth respiration for each pool [kg]
    real(r8) :: total_growth_respiration                   ! Total plant respiration (kgC)

    integer, dimension(num_organs) :: curpri_c_ids         ! C variable ID's of the current priority level
    integer, dimension(num_organs) :: curpri_n_ids         ! N variable ID's of the current priority level
    integer, dimension(num_organs) :: curpri_p_ids         ! P variable ID's of the current priority level

    ! Integegrator variables
    ! These are not global because we want a unique instance for each time the routine is called
    ! ----------------------------------------------------------------------------------------

    real(r8),dimension(num_intgr_vars) :: state_array      ! Vector of carbon pools passed to integrator
    real(r8),dimension(num_intgr_vars) :: state_array_out  ! Vector of carbon pools passed back from integrator
    logical,dimension(num_intgr_vars)  :: state_mask       ! Mask of active pools during integration

    integer, parameter  :: n_max_priority = num_organs + 1 ! Maximum possible number of priority levels is
                                                           ! the total number organs plus 1, which allows
                                                           ! each organ to have its own level, and ignore
                                                           ! the specialized priority 1

    integer , parameter :: max_substeps = 300              ! Maximum allowable iterations
    real(r8), parameter :: max_trunc_error = 1.0_r8        ! Maximum allowable truncation error
    integer,  parameter :: ODESolve = 2                    ! 1=RKF45,  2=Euler

    real(r8), parameter :: store_overflow_frac = 0.15      ! The fraction above target allowed in storage
                                                           ! to accept overflow, or extra carbon at end of allocation

    ! This is a local array containing the boundary conditions
    ! we need this (for now at least) because the integration layer needs things
    ! packed into simple types

    real(r8) ::  intgr_params(num_bc_in)


    ! Create a large number of aliases used for convenience

    associate( & 

          ! Carbon State
          leaf_c   => this%variables(leaf_c_id)%val, &
          fnrt_c   => this%variables(fnrt_c_id)%val(icd), &
          sapw_c   => this%variables(sapw_c_id)%val(icd), &
          store_c  => this%variables(store_c_id)%val(icd), &
          repro_c  => this%variables(repro_c_id)%val(icd), &
          struct_c => this%variables(struct_c_id)%val(icd), &

          ! Carbon Targets
          leaf_c_target   => v_target(leaf_c_id), &
          fnrt_c_target   => v_target(fnrt_c_id), &
          sapw_c_target   => v_target(sapw_c_id), &
          store_c_target  => v_target(store_c_id), &
          repro_c_target  => v_target(repro_c_id), &
          struct_c_target => v_target(struct_c_id), &
          
          ! Carbon Demand
          leaf_c_demand   => v_demand(leaf_c_id), &
          fnrt_c_demand   => v_demand(fnrt_c_id), &
          sapw_c_demand   => v_demand(sapw_c_id), &
          store_c_demand  => v_demand(store_c_id), &
          repro_c_demand  => v_demand(repro_c_id), &
          struct_c_demand => v_demand(struct_c_id), &

          ! Nitrogen State
          leaf_n   => this%variables(leaf_n_id)%val, &
          fnrt_n   => this%variables(fnrt_n_id)%val(icd), &
          sapw_n   => this%variables(sapw_n_id)%val(icd), &
          store_n  => this%variables(store_n_id)%val(icd), &
          repro_n  => this%variables(repro_n_id)%val(icd), &
          struct_n => this%variables(struct_n_id)%val(icd), &

          ! Phosphorus State
          leaf_p   => this%variables(leaf_p_id)%val, &
          fnrt_p   => this%variables(fnrt_p_id)%val(icd), &
          sapw_p   => this%variables(sapw_p_id)%val(icd), &
          store_p  => this%variables(store_p_id)%val(icd), &
          repro_p  => this%variables(repro_p_id)%val(icd), &
          struct_p => this%variables(struct_p_id)%val(icd), &

          ! This is the list of organs
          organ_map => prt_global%organ_map, &
          
          ! This is the list of organ x element
          sp_organ_map => prt_global%sp_organ_map)


    ! Copy the in/out boundary conditions as pointers
    ! -----------------------------------------------------------------------------------

    dbh                               => this%bc_inout(acnp_bc_inout_id_dbh)%rval
    carbon_gain                       => this%bc_inout(acnp_bc_inout_id_netdc)%rval
    maint_r_deficit                   => this%bc_inout(acnp_bc_inout_id_rmaint_def)%rval
    nitrogen_gain                     => this%bc_inout(acnp_bc_inout_id_netdn)%rval
    phosphorus_gain                  => this%bc_inout(acnp_bc_inout_id_netdp)%rval

    ! Copy the in boundary conditions into readable local variables  
    ! We don't use pointers, because inputs should be intent in only
    ! -----------------------------------------------------------------------------------

    canopy_trim                       = this%bc_in(acnp_bc_in_id_ctrim)%rval
    leaves_on                         = this%bc_in(acnp_bc_in_id_leafon)%ival
    ipft                              = this%bc_in(acnp_bc_in_id_pft)%ival

    r_g(leaf_organ)                    = EDPftvarcon_inst%prt_grperc_organ(ipft,leaf_organ)
    r_g(fnrt_organ)                    = EDPftvarcon_inst%prt_grperc_organ(ipft,fnrt_organ)
    r_g(sapw_organ)                    = EDPftvarcon_inst%prt_grperc_organ(ipft,sapw_organ)
    r_g(store_organ)                   = EDPftvarcon_inst%prt_grperc_organ(ipft,store_organ)
    r_g(struct_organ)                  = EDPftvarcon_inst%prt_grperc_organ(ipft,struct_organ)
    r_g(repro_organ)                   = EDPftvarcon_inst%prt_grperc_organ(ipft,repro_organ)

    intgr_params(:)                    = -9.9e32_r8
    intgr_params(acnp_bc_in_id_ctrim) = this%bc_in(acnp_bc_in_id_ctrim)%rval
    intgr_params(acnp_bc_in_id_pft)   = real(this%bc_in(acnp_bc_in_id_pft)%ival)

    ! All growth respiration is handled internally, and is tallied in this
    ! daily routine.  Initialize it here.
    total_growth_respiration = 0.0_r8
    
    
    ! Initialize some fields at the beginning of the routine for balance checking later on
    carbon_gain0      = carbon_gain
    nitrogen_gain0    = nitrogen_gain
    phosphorus_gain0 = phosphorus_gain
    maint_r_deficit0   = maint_r_deficit


    ! Initialize Targets and Demands to uninitialized flag
    v_target(:) = -999.9_r8
    v_demand(:) = -999.9_r8

    ! Initialize growth respiration to zero
    growth_resp(:) = 0.0_r8

    ! -----------------------------------------------------------------------------------
    ! I. Remember the values for the state variables at the beginning of this
    ! routines. We will then use that to determine their net allocation and reactive
    ! transport flux "%net_alloc" at the end.
    ! -----------------------------------------------------------------------------------

    nleafage = prt_global%state_descriptor(leaf_c_id)%num_pos ! Number of leaf age class

    leaf_c0(1:nleafage) = leaf_c(1:nleafage)  ! Set initial leaf carbon 
    leaf_n0(1:nleafage) = leaf_n(1:nleafage)
    leaf_p0(1:nleafage) = leaf_p(1:nleafage)
    fnrt_c0 = fnrt_c                          ! Set initial fine-root carbon
    sapw_c0 = sapw_c                          ! Set initial sapwood carbon
    store_c0 = store_c                        ! Set initial storage carbon 
    repro_c0 = repro_c                        ! Set initial reproductive carbon
    struct_c0 = struct_c                      ! Set initial structural carbon
    fnrt_n0 = fnrt_n                          ! Set initial fine-root carbon
    sapw_n0 = sapw_n                          ! Set initial sapwood carbon
    store_n0 = store_n                        ! Set initial storage carbon 
    repro_n0 = repro_n                        ! Set initial reproductive carbon
    struct_n0 = struct_n                      ! Set initial structural carbon
    fnrt_p0 = fnrt_p                          ! Set initial fine-root carbon
    sapw_p0 = sapw_p                          ! Set initial sapwood carbon
    store_p0 = store_p                        ! Set initial storage carbon 
    repro_p0 = repro_p                        ! Set initial reproductive carbon
    struct_p0 = struct_p                      ! Set initial structural carbon

    ! -----------------------------------------------------------------------------------
    ! If we have more than one leaf age classification, allow
    ! some leaf biomass to transition to the older classes.  NOTE! This is not handling
    ! losses due to turnover (ie. flux from the oldest senescing class). This is only
    ! internal.
    ! -----------------------------------------------------------------------------------

    if(nleafage>1) then
       do i_age = 1,nleafage-1
          if (EDPftvarcon_inst%leaf_long(ipft,i_age)>nearzero) then

             leaf_age_flux_frac = years_per_day / EDPftvarcon_inst%leaf_long(ipft,i_age)

             leaf_c(i_age)    = leaf_c(i_age)   - leaf_c0(i_age) * leaf_age_flux_frac
             leaf_c(i_age+1)  = leaf_c(i_age+1) + leaf_c0(i_age) * leaf_age_flux_frac

             leaf_n(i_age)    = leaf_n(i_age)   - leaf_n0(i_age) * leaf_age_flux_frac
             leaf_n(i_age+1)  = leaf_n(i_age+1) + leaf_n0(i_age) * leaf_age_flux_frac

             leaf_p(i_age)    = leaf_p(i_age)   - leaf_p0(i_age) * leaf_age_flux_frac
             leaf_p(i_age+1)  = leaf_p(i_age+1) + leaf_p0(i_age) * leaf_age_flux_frac

          end if
       end do
    end if

   
    ! -----------------------------------------------------------------------------------
    ! I. Calculate target size of the carbon components of each organ
    ! -----------------------------------------------------------------------------------
    
    ! Target sapwood biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
    call bsap_allom(dbh,ipft,canopy_trim,sapw_area,sapw_c_target, sapw_dcdd_target)

    ! Target total above ground deriv. biomass in woody/fibrous tissues  [kgC, kgC/cm]
    call bagw_allom(dbh,ipft,agw_c_target,agw_dcdd_target)

    ! Target total below ground deriv. biomass in woody/fibrous tissues [kgC, kgC/cm] 
    call bbgw_allom(dbh,ipft,bgw_c_target, bgw_dcdd_target)

    ! Target total dead (structrual) biomass and deriv. [kgC, kgC/cm]
    call bdead_allom( agw_c_target, bgw_c_target, sapw_c_target, ipft, struct_c_target, &
                      agw_dcdd_target, bgw_dcdd_target, sapw_dcdd_target, struct_dcdd_target)
    
    ! Target leaf biomass according to allometry and trimming
    call bleaf(dbh,ipft,canopy_trim, leaf_c_target, leaf_dcdd_target)

    ! Target fine-root biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
    call bfineroot(dbh,ipft,canopy_trim, fnrt_c_target, fnrt_dcdd_target)

    ! Target storage carbon [kgC,kgC/cm]
    call bstore_allom(dbh,ipft,canopy_trim, store_c_target, store_dcdd_target)

    ! -----------------------------------------------------------------------------------
    ! Calculate the allometric deficits for all pools
    ! Note that we are factoring in the growth respiration taxes on the carbon pools
    ! 
    ! Capping to positive values is necessary, some pools may be above targets
    ! because of either numerical integration errors, or in the case of storage, 
    ! we may had been limited in the last growth step by nutrients, and had 
    ! overflow.
    ! -----------------------------------------------------------------------------------
    
    leaf_c_demand   = max(0.0_r8,(leaf_c_target - sum(leaf_c)) * (1.0_r8 + r_g(leaf_organ)))
    fnrt_c_demand   = max(0.0_r8,(fnrt_c_target - fnrt_c) * (1.0_r8 + r_g(fnrt_organ)))
    sapw_c_demand   = max(0.0_r8,(sapw_c_target - sapw_c) * (1.0_r8 + r_g(sapw_organ)))
    struct_c_demand = max(0.0_r8,(struct_c_target - struct_c) * (1.0_r8 + r_g(struct_organ)))
    store_c_demand  = max(0.0_r8,(store_c_target - store_c) * (1.0_r8 + r_g(store_organ)))

    ! Make sure that reproductive demand is always zero until stature growth step. We only
    ! allocate a portion during that part of the algorithm
    repro_c_demand  = 0.0_r8


    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    ! Preferential transfer of available carbon and nutrients into the highest 
    ! priority pools, and maintenance respiration. We will loop through the available
    ! pools, and identify if that pool is part of the highest transfer priority.
    ! If it is, then we track the variable ids associated with that pool for each CNP
    ! species.  It "should" work fine if there are NO priority=1 pools...
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

    curpri_c_ids(:) = 0    ! "current-priority" carbon variable id's
    curpri_n_ids(:) = 0    ! "current-priority" nitrogen variable id's
    curpri_p_ids(:) = 0    ! "current-priority" phosphorus variable id's
    
    
    i = 0
    do ii = 1, num_organs
       
       ! The priority code associated with this organ
       priority_code = int(EDPftvarcon_inst%prt_alloc_priority(ipft, organ_list(ii)))
       
       ! 1 is the highest priority code possible
       if( priority_code == 1 ) then
          i = i + 1
          curpri_c_ids(i) = sp_organ_map(organ_list(ii),carbon12_element)
          curpri_n_ids(i) = sp_organ_map(organ_list(ii),nitrogen_element)
          curpri_p_ids(i) = sp_organ_map(organ_list(ii),phosphorus_element)
       end if
    end do
    
    ! Number of pools in the current priority level
    num_organs_curpri = i
    
    ! -----------------------------------------------------------------------------------
    ! The high-priority pools, and their associated variable
    ! indices have been identified.
    ! 
    ! Let us now calculate the fluxes into these priority pools and replace
    ! maintenance respiration demand.  
    !
    ! Note: We will only allow transfers during this stage, 
    ! into tissues if the plant has a leaf-on phenology status.  We only 
    ! apply that filter here, because this is the only place we will
    ! use storage carbon to replace tissues, otherwise we use carbon
    ! gains which would be zero anyway if leaves are off.
    ! -----------------------------------------------------------------------------------
    
    sum_c_demand = 0.0_r8
    if(leaves_on == itrue) then
       do i = 1, num_organs_curpri
          i_cvar       = curpri_c_ids(i)
          sum_c_demand = sum_c_demand + v_demand(i_cvar)
       end do
    end if
    
    if( (sum_c_demand+maint_r_deficit)>nearzero .or. carbon_gain < 0.0_r8   ) then
    
       ! the actual amount of carbon that will be used to fulfill the priority flux
       ! to both organ pools and maintenance respiration.
       ! All of these pools should be positive in this next calculation.
       
       if (carbon_gain >= 0.0_r8) then

          ! The transferable storage to leaf/fine-root/maintenance is a linear function
          ! of availability
          store_c_transferable   = store_c*min(1.0_r8,store_c / store_c_target)

          mo_sum_c_flux   = min(sum_c_demand + maint_r_deficit, &
               store_c_transferable + carbon_gain )
          
          ! The fraction of the carbon that that will actually be transferred per
          ! the total needed
          mo_sum_c_frac   = mo_sum_c_flux / (sum_c_demand + maint_r_deficit)
          
          p_m = EDPftvarcon_inst%leaf_stor_priority(ipft)
          if( p_m >= 0.5) then
             redist_c_flux  = min( (p_m - 0.5_r8)/0.5_r8 * mo_sum_c_frac * maint_r_deficit , &
                  (1.0_r8 - mo_sum_c_frac) * sum_c_demand )
             
             sum_c_flux  = mo_sum_c_frac * sum_c_demand + redist_c_flux
             maint_r_def_flux = mo_sum_c_frac * maint_r_deficit   - redist_c_flux
             
          else
             redist_c_flux  = min( (0.5_r8 - p_m)/0.5_r8 * mo_sum_c_frac * sum_c_demand, &
                  (1.0_r8 - mo_sum_c_frac) * maint_r_deficit ) 
             
             sum_c_flux  = mo_sum_c_frac * sum_c_demand - redist_c_flux
             maint_r_def_flux = mo_sum_c_frac * maint_r_deficit   + redist_c_flux
             
          end if
          
       else

          ! In some model configurations, we may pass in a negative
          ! carbon balance, a mode where we are not tracking maintenance
          ! respiration deficit, and instead, we are paying it by
          ! creating something of a negative NPP for the the day.
          ! In this csae, we assume that the maintenance respiration,
          ! and therefore the negative carbon gain, has already
          ! factored in storage carbon availability. So we do not impoase
          ! any limitations here, other than a check to make sure the
          ! pool isnt depleted into a negative value. ie.. all of it
          ! is transferrable.
          
          ! note that in cases where we instead track a maintenance
          ! respiration deficit, and pay for it in this routine,
          ! we will never pass in a negative carbon gain, and thus
          ! the transferrable storage carbon is calculated here.

          carbon_gain_flux = -carbon_gain
          
          if ( carbon_gain_flux > store_c ) then
             write(fates_log(),*) 'There is not enough transferable storage carbon'
             write(fates_log(),*) 'to pay back negative daily carbon gain.'
             write(fates_log(),*) 'daily carbon_gain:',carbon_gain
             write(fates_log(),*) 'transferrable carbon:',store_c_transferable
             write(fates_log(),*) 'Difference gain+transferable:',carbon_gain+store_c_transferable
             write(fates_log(),*) 'Exiting.'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          ! When carbon gain is negative, the carbon transferred from storage
          ! to pay back negative gains is prioritized.
          ! Then we see what is left, and are willing to transfer out of storage
          ! That value is sent to priority tissues.
          
          store_c_transferable = max(0.0_r8,store_c-carbon_gain_flux) * &
                                 min(1.0_r8,(store_c-carbon_gain_flux) / store_c_target)
          sum_c_flux           = min(store_c_transferable,sum_c_demand)


       end if
       
       ! -----------------------------------------------------------------------------------
       ! The flux fractions of carbon have been calculated, now transfer
       ! their respective pools and out of carbon gain and storage
       ! -----------------------------------------------------------------------------------
       
       do i = 1, num_organs_curpri
          
          ! The local index for this variable
          i_cvar   = curpri_c_ids(i)
          i_nvar   = curpri_n_ids(i)
          i_pvar   = curpri_p_ids(i)
          
          ! The global index for this pool (for pulling parameters)
          i_gorgan = prt_global%state_descriptor(i_cvar)%organ_id
          
          ! The total carbon drawn for this pool
          if(sum_c_demand > nearzero) then
             c_flux   = sum_c_flux * v_demand(i_cvar) / sum_c_demand
          else
             c_flux   = 0.0_r8
          end if
          
          ! Transfer the carbon into the pool
          ! Note* If this is a variable with age bins (ie leaves)
          ! we are also using bin 1 for all fluxes at this stage (growth bin)
          ! which is also "icd"

          this%variables(i_cvar)%val(icd) = this%variables(i_cvar)%val(icd) + &
                                            c_flux/(1.0_r8 + r_g(i_gorgan))
          
          ! Transfer the carbon into the pool growth respiration cost
          growth_resp(i_cvar) = growth_resp(i_cvar) + & 
                c_flux * r_g(i_gorgan)/(1.0_r8 + r_g(i_gorgan))

          ! Update the carbon demand (based on all bins
          v_demand(i_cvar) = max(0.0_r8,(v_target(i_cvar) - &
               sum(this%variables(i_cvar)%val(:))) * (1.0_r8 + r_g(i_gorgan)))
          
          
          ! Update the nitrogen demands (which are based off of carbon actual..)
          ! Note that the nitrogen target is tied to the stoichiometry of thegrowing pool only
          v_target(i_nvar) = this%variables(i_cvar)%val(icd)*EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,i_gorgan)
          v_demand(i_nvar) = max(0.0_r8, v_target(i_nvar) - this%variables(i_nvar)%val(icd))
          
          
          ! Update the phosphorus demands (which are based off of carbon actual..)
          ! Note that the phsophorus target is tied to the stoichiometry of thegrowing pool only (also)
          v_target(i_pvar) = this%variables(i_cvar)%val(icd)*EDPftvarcon_inst%prt_phos_stoich_p1(ipft,i_gorgan)
          v_demand(i_pvar) = max(0.0_r8, v_target(i_pvar) - this%variables(i_pvar)%val(icd))
          
       end do
    
       ! Update the maintenance demand (this needs a tarrif.  the maintenance deficit 
       ! is not playing fair, and is making us the laughing stock of the plant allocation world
       
       if(carbon_gain >= 0.0_r8) then
          
          maint_r_deficit  = maint_r_deficit - maint_r_def_flux
          
          ! Carbon flux coming from daily gains
          carbon_gain_flux = min(carbon_gain, mo_sum_c_flux)
          
          carbon_gain = carbon_gain - carbon_gain_flux

          ! Draw from storage (if necessary)
          store_c_flux =  max(0.0_r8,  mo_sum_c_flux - carbon_gain_flux)

       else

          ! Case where we have negative carbon gains (we aren't tracking
          ! maintenance R deficit in this case)

          carbon_gain = carbon_gain + carbon_gain_flux

          store_c_flux = carbon_gain_flux + sum_c_flux

       end if
       
       store_c     = store_c     - store_c_flux
       
       ! Update storeage demand
       i_gorgan = prt_global%state_descriptor(store_c_id)%organ_id
       v_demand(store_c_id) = &
            max(0.0_r8,(v_target(store_c_id) - &
            this%variables(store_c_id)%val(icd)) * (1.0_r8 + r_g(i_gorgan)))
       


    end if  ! if( (sum_c_demand+maint_r_deficit)>nearzero ) then

    ! -----------------------------------------------------------------------------------
    ! Calculate highest priority nutrient (NP) fluxes
    ! -----------------------------------------------------------------------------------

    sum_n_demand = 0.0_r8
    sum_p_demand = 0.0_r8
    do i = 1, num_organs_curpri
       ! The local index for this variable
       i_nvar   = curpri_n_ids(i)
       i_pvar   = curpri_p_ids(i)
       sum_n_demand = sum_n_demand + v_demand(i_nvar)
       sum_p_demand = sum_p_demand + v_demand(i_pvar)
    end do

    if (sum_n_demand>nearzero) then
       sum_n_flux = min(nitrogen_gain, sum_n_demand) 
       do i = 1, num_organs_curpri
          i_nvar   = curpri_n_ids(i)
          n_flux   = sum_n_flux * v_demand(i_nvar)/sum_n_demand

          ! Update the state (all mass is placed in 1st age bin if exists)
          this%variables(i_nvar)%val(icd) = this%variables(i_nvar)%val(icd) + n_flux

          ! Update the demand
          v_demand(i_nvar) = max(0.0_r8,v_target(i_nvar) - this%variables(i_nvar)%val(icd))

       end do
    end if
    if (sum_p_demand>nearzero) then
       sum_p_flux = min(phosphorus_gain, sum_p_demand) 
       do i = 1, num_organs_curpri
          i_pvar = curpri_p_ids(i)
          p_flux = sum_p_flux * v_demand(i_pvar)/sum_p_demand

          ! Update the state
          this%variables(i_pvar)%val(icd) = this%variables(i_pvar)%val(icd) + p_flux
          
          ! Update the demand
          v_demand(i_pvar) = max(0.0_r8,v_target(i_pvar) - this%variables(i_pvar)%val(icd))

       end do
    end if

    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    !
    ! Bring all pools, in priority order, up to allometric targets if possible
    ! We have already done priority 1, so start with 2
    !
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

    do i_priority = 2, n_max_priority

       curpri_c_ids(:) = 0    ! "current-priority" carbon variable id's
       curpri_n_ids(:) = 0    ! "current-priority" nitrogen variable id's
       curpri_p_ids(:) = 0    ! "current-priority" phosphorus variable id's

       i = 0
       do ii = 1, num_organs
          
          ! The priority code associated with this organ
          priority_code = int(EDPftvarcon_inst%prt_alloc_priority(ipft, organ_list(ii)))
          
          ! 1 is the highest priority code possible
          if( priority_code == i_priority ) then
             i = i + 1
             curpri_c_ids(i) = sp_organ_map(organ_list(ii),carbon12_element)
             curpri_n_ids(i) = sp_organ_map(organ_list(ii),nitrogen_element)
             curpri_p_ids(i) = sp_organ_map(organ_list(ii),phosphorus_element)
          end if
       end do
       
       ! Bring carbon up to target first, this order is required
       ! because we need to know the resulting carbon concentrations
       ! before  we set the allometric targets for the nutrients

       num_organs_curpri  = i
       sum_c_demand = 0.0_r8
       
       do i = 1, num_organs_curpri
          ! The local index for this variable
          i_cvar   = curpri_c_ids(i)
          sum_c_demand = sum_c_demand + v_demand(i_cvar)
       end do
       
       sum_c_flux = min(carbon_gain, sum_c_demand)

       do i = 1, num_organs_curpri
          i_cvar   = curpri_c_ids(i)
          i_nvar   = curpri_n_ids(i)
          i_pvar   = curpri_p_ids(i)

          if (sum_c_demand>nearzero) then
             c_flux = sum_c_flux * v_demand(i_cvar)/sum_c_demand
          else
             c_flux = 0.0_r8
          end if

          i_gorgan = prt_global%state_descriptor(i_cvar)%organ_id

          ! Update the carbon pool
          this%variables(i_cvar)%val(icd) = this%variables(i_cvar)%val(icd) + c_flux/(1.0_r8 + r_g(i_gorgan))

          ! Estimate the growth respiration associated with this allocation
          growth_resp(i_cvar) =  growth_resp(i_cvar) + &
                c_flux * r_g(i_gorgan)/(1.0_r8 + r_g(i_gorgan))

          
          ! Update all targets and demands now that the carbon pool has been updated

          ! carbon demand
          v_demand(i_cvar) = max(0.0_r8,(v_target(i_cvar) - &
               sum(this%variables(i_cvar)%val(:)))*(1.0_r8 + r_g(i_gorgan)))
          
          ! nitrogen demand (Note that nitrogen demand is only relevant in the growing pool)

          v_target(i_nvar) = this%variables(i_cvar)%val(icd) * &
               EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,i_gorgan)

          v_demand(i_nvar) = max(0.0_r8,v_target(i_nvar) - this%variables(i_nvar)%val(icd))
          
          ! phosphorus demand (note that phosphorus demand is only relevant in the growing pool)
          v_target(i_pvar) = this%variables(i_cvar)%val(icd) * &
               EDPftvarcon_inst%prt_phos_stoich_p1(ipft,i_gorgan)
          
          v_demand(i_pvar) = max(0.0_r8,v_target(i_pvar) - this%variables(i_pvar)%val(icd))
          
          carbon_gain      = carbon_gain - c_flux

       end do

       ! Bring nutrients up to target for this priority level

       sum_n_demand = 0.0_r8
       sum_p_demand = 0.0_r8
       
       do i = 1, num_organs_curpri
          ! The local index for this variable
          i_nvar   = curpri_n_ids(i)
          i_pvar   = curpri_p_ids(i)
          sum_n_demand = sum_n_demand + v_demand(i_nvar)
          sum_p_demand = sum_p_demand + v_demand(i_pvar)
       end do
       
       sum_n_flux = min(nitrogen_gain, sum_n_demand) 
       sum_p_flux = min(phosphorus_gain, sum_p_demand)

       do i = 1, num_organs_curpri
          i_cvar   = curpri_c_ids(i)
          i_nvar   = curpri_n_ids(i)
          i_pvar   = curpri_p_ids(i)

          if (sum_n_demand>nearzero) then
             n_flux = sum_n_flux * v_demand(i_nvar)/sum_n_demand
          else
             n_flux = 0.0_r8
          end if
          if (sum_p_demand>nearzero) then
             p_flux = sum_p_flux * v_demand(i_pvar)/sum_p_demand
          else
             p_flux = 0.0_r8
          end if

          i_gorgan = prt_global%state_descriptor(i_cvar)%organ_id

          this%variables(i_nvar)%val(icd) = this%variables(i_nvar)%val(icd) + n_flux
          this%variables(i_pvar)%val(icd) = this%variables(i_pvar)%val(icd) + p_flux

          ! update nutrient demands, but not targets (carbon already changed, targets already set)
          v_demand(i_nvar) = max(0.0_r8,v_target(i_nvar) - this%variables(i_nvar)%val(icd))
          v_demand(i_pvar) = max(0.0_r8,v_target(i_pvar) - this%variables(i_pvar)%val(icd))
          
          nitrogen_gain    = nitrogen_gain - n_flux
          phosphorus_gain = phosphorus_gain - p_flux

       end do

    end do
    
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    ! Attempts have been made to get all pools and species closest to allometric
    ! targets based on prioritized relative demand and allometry functions.
    ! Now, grow tissues beyond current targets.  Growth is limited by how far each pool
    ! is allowed to extend beyond their nutrient-carbon ratios.  Some pools may allow
    ! some slip room, some may allow none.  This is specified as the fraction of mass
    ! deficit allowed, relative to the current target value.
    ! parteh_n_growth_regulation
    ! 0 = the pool must be at or above the current stoichiometric target to allow growth
    ! 1 = the pool is allowed to grow up to and until its the nutrient concentration
    !     deficit has dropped to 100% of the current target.
    ! 2.5 = similar to 1, but 250% of its current target
    !
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    !
    ! For carbon:    there is a throttle down on how much can be transferred
    !                based on how far off stoichiometry the nutrients are.
    ! For nutrients: Follows after the carbon step.  Is again, just based off of
    !                proportionality, or, relative demand.
    ! -----------------------------------------------------------------------------------

    ! Lets forcast each carbon pools relative needs by extrapolating the 
    ! slope at the current point.  We will use these slopes to anticipate
    ! which nutrient pool may limit growth
    ! The remaining carbon nugget may or may not be down-regulated based
    ! on nutrient limitations

    ! the carbon balance sub-step (deltaC) will be halved and tried again
    ! -----------------------------------------------------------------------------------
    
    if( carbon_gain > nearzero ) then
       
       state_mask(:) = .false.

       ! Go through and flag the integrating variables as either pools that 
       ! are growing in this iteration, or not.   At this point, if carbon for growth
       ! remains, it means that all pools are up to, or above the target.  If
       ! it is above, that is because of numerical integration errors, or fusion.
       ! In that case, we flag that pool to not be included in stature growth. It will
       ! catch up with the other pools in the next stature growth steps.

       do i = 1, num_organs
          
          i_cvar = sp_organ_map(organ_list(i),carbon12_element)
          
          cdeficit =  v_target(i_cvar) - sum(this%variables(i_cvar)%val(:))

          if ( cdeficit > calloc_abs_error ) then
             ! In this case, we somehow still have carbon to play with,
             ! yet one of the pools is below its current target
             ! gracefully fail
             write(fates_log(),*) 'A carbon pool has reached the stature growth step'
             write(fates_log(),*) 'yet its deficit is too large to integrate '
             write(fates_log(),*) 'organ: ',i
             write(fates_log(),*) cdeficit, v_target(i_cvar),sum(this%variables(i_cvar)%val(:))
             call endrun(msg=errMsg(sourcefile, __LINE__))
             
          elseif( (-cdeficit) > calloc_abs_error ) then
             ! In this case, we are above our target (ie negative deficit (fusion?))
             ! and if so, this pool does not have to grow and will
             ! catch up in the next step
             state_mask(i)            = .false.    ! flag index of state variable
             state_mask(i+num_organs) = .false.    ! flag index for its growth respiration too
             
          else
             ! In this case, the pool is close enough to the target
             ! to be flagged for growth

             state_mask(i)            = .true.    ! flag index of state variable
             state_mask(i+num_organs) = .true.    ! flag index for its growth respiration too
          end if
       end do


       ! fraction of carbon going towards reproduction reproductive carbon is 
       ! just different from the other pools. It is not based on proportionality, 
       ! so its mask is set differently.  We (inefficiently) just included
       ! reproduction in the previous loop, but oh well, we over-write now.
       
       if (dbh <= EDPftvarcon_inst%dbh_repro_threshold(ipft)) then
          repro_c_frac = EDPftvarcon_inst%seed_alloc(ipft)
       else
          repro_c_frac = EDPftvarcon_inst%seed_alloc(ipft) + EDPftvarcon_inst%seed_alloc_mature(ipft)
       end if
       
       if(repro_c_frac>nearzero)then
          state_mask(repro_organ)            = .true.
          state_mask(repro_organ+num_organs) = .true.
       else
          state_mask(repro_organ)            = .false.
          state_mask(repro_organ+num_organs) = .false.
       end if

       ! Calculate the total CARBON allocation rate per diameter increment
       ! Include the growth respiration costs "total_dcostdd" in that estimate
       ! ALso include the non-respiring cost, which is needed to project reproductive
       ! costs.
       ! --------------------------------------------------------------------------------
       
       ! First objective is to find the extrapolated proportions of carbon going to
       ! each pool. This has nothing to do with carbon conservation, it is just used
       ! to make a rough prediction of how much nutrient is needed to match carbon

       total_dcostdd = struct_dcdd_target
       if (state_mask(leaf_c_id)) then
          total_dcostdd = total_dcostdd + leaf_dcdd_target
       end if
       if (state_mask(fnrt_c_id)) then
          total_dcostdd = total_dcostdd + fnrt_dcdd_target
       end if
       if (state_mask(sapw_c_id)) then
          total_dcostdd = total_dcostdd + sapw_dcdd_target
       end if
       if (state_mask(store_c_id)) then
          total_dcostdd = total_dcostdd + store_dcdd_target
       end if

       struct_c_frac = struct_dcdd_target/total_dcostdd * (1.0_r8 - repro_c_frac)
       if (state_mask(leaf_c_id)) then
          leaf_c_frac = leaf_dcdd_target/total_dcostdd * (1.0_r8 - repro_c_frac)
       else
          leaf_c_frac = 0.0_r8
       end if
       if (state_mask(intgr_fnrt_c_id)) then
          fnrt_c_frac = fnrt_dcdd_target/total_dcostdd * (1.0_r8 - repro_c_frac)
       else
          fnrt_c_frac = 0.0_r8
       end if
       if (state_mask(intgr_sapw_c_id)) then
          sapw_c_frac = sapw_dcdd_target/total_dcostdd * (1.0_r8 - repro_c_frac)
       else
          sapw_c_frac = 0.0_r8
       end if
       if (state_mask(intgr_store_c_id)) then
          store_c_frac = store_dcdd_target/total_dcostdd * (1.0_r8 - repro_c_frac)
       else
          store_c_frac = 0.0_r8
       end if
       
       ! Calculate an approximation of the total amount of carbon that would be needed
       ! to match the amount of each nutrient used.  We also add in the amount of nutrient
       ! that may or may-not exist above each pool's minimum stoichiometry...
       ! --------------------------------------------------------------------------------
       
       grow_c_from_c = carbon_gain * struct_c_frac / (1.0_r8 + r_g(struct_organ))
       grow_c_from_n = nitrogen_gain * struct_c_frac / EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,struct_organ)
       grow_c_from_p = phosphorus_gain * struct_c_frac / EDPftvarcon_inst%prt_phos_stoich_p1(ipft,struct_organ)

       if (state_mask(intgr_leaf_c_id))  then
          grow_c_from_c = grow_c_from_c + carbon_gain * leaf_c_frac / (1.0_r8 + r_g(leaf_organ))
          grow_c_from_n = grow_c_from_n + nitrogen_gain * leaf_c_frac / EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,leaf_organ)
          grow_c_from_p = grow_c_from_p + phosphorus_gain * leaf_c_frac /  EDPftvarcon_inst%prt_phos_stoich_p1(ipft,leaf_organ)

          ! Nutrient pools may already be above the current minimum target, add that head start onto the equiavelnt
          ! value so the "head-start" quantity is factored in...
          grow_c_from_n = grow_c_from_n + &
                max(0.0_r8,this%variables(leaf_n_id)%val(icd) - v_target(leaf_n_id)) / &
                EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,leaf_organ)
          grow_c_from_p = grow_c_from_p + &
                max(0.0_r8,this%variables(leaf_p_id)%val(icd) - v_target(leaf_p_id)) / &
                EDPftvarcon_inst%prt_phos_stoich_p1(ipft,leaf_organ)

       end if
       if (state_mask(intgr_fnrt_c_id)) then
          grow_c_from_c = grow_c_from_c + carbon_gain * fnrt_c_frac / (1.0_r8 + r_g(fnrt_organ))
          grow_c_from_n = grow_c_from_n + nitrogen_gain * fnrt_c_frac / EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,fnrt_organ)
          grow_c_from_p = grow_c_from_p + phosphorus_gain * fnrt_c_frac / EDPftvarcon_inst%prt_phos_stoich_p1(ipft,fnrt_organ)

          grow_c_from_n = grow_c_from_n + &
                max(0.0_r8,this%variables(fnrt_n_id)%val(icd) - v_target(fnrt_n_id)) / &
                EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,fnrt_organ)
          grow_c_from_p = grow_c_from_p + &
                max(0.0_r8,this%variables(fnrt_p_id)%val(icd) - v_target(fnrt_p_id)) / &
                EDPftvarcon_inst%prt_phos_stoich_p1(ipft,fnrt_organ)

       end if
       if (state_mask(intgr_sapw_c_id)) then
          grow_c_from_c = grow_c_from_c + carbon_gain * sapw_c_frac / (1.0_r8 + r_g(sapw_organ))
          grow_c_from_n = grow_c_from_n + nitrogen_gain * sapw_c_frac / EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,sapw_organ)
          grow_c_from_p = grow_c_from_p + phosphorus_gain * sapw_c_frac / EDPftvarcon_inst%prt_phos_stoich_p1(ipft,sapw_organ)

          grow_c_from_n = grow_c_from_n + &
                max(0.0_r8,this%variables(sapw_n_id)%val(icd) - v_target(sapw_n_id)) / &
                EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,sapw_organ)
          grow_c_from_p = grow_c_from_p + &
                max(0.0_r8,this%variables(sapw_p_id)%val(icd) - v_target(sapw_p_id)) / &
                EDPftvarcon_inst%prt_phos_stoich_p1(ipft,sapw_organ)


       end if
       if (state_mask(intgr_store_c_id)) then
          grow_c_from_c = grow_c_from_c + carbon_gain * store_c_frac / &
                          (1.0_r8 + r_g(store_organ))
          grow_c_from_n = grow_c_from_n + nitrogen_gain * store_c_frac / &
                          EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,store_organ)
          grow_c_from_p = grow_c_from_p + phosphorus_gain * store_c_frac / &
                          EDPftvarcon_inst%prt_phos_stoich_p1(ipft,store_organ)

          grow_c_from_n = grow_c_from_n + &
                max(0.0_r8,this%variables(store_n_id)%val(icd) - v_target(store_n_id)) / &
                EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,store_organ)
          grow_c_from_p = grow_c_from_p + &
                max(0.0_r8,this%variables(store_p_id)%val(icd) - v_target(store_p_id)) / &
                EDPftvarcon_inst%prt_phos_stoich_p1(ipft,store_organ)


       end if
       if (state_mask(intgr_repro_c_id)) then
          grow_c_from_c = grow_c_from_c + carbon_gain * repro_c_frac / (1.0_r8 + r_g(repro_organ))
          grow_c_from_n = grow_c_from_n + nitrogen_gain * repro_c_frac / EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ)
          grow_c_from_p = grow_c_from_p + phosphorus_gain * repro_c_frac / EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ)

          ! The amount of existing reproductive material has no bearing on how much 
          ! is allocated.

       end if
       
       ! --------------------------------------------------------------------------------
       ! We limit growth to align with the species would motivate the least flux of 
       ! carbon into growing tissues to match.  This is only an approximation of how much 
       ! growth we get out of each, and they don't have to be perfect.  As long as we 
       ! don't use more carbon than we have (we wont) and if we use the actual numerical 
       ! integrator in the trasfer step, the nutrients will be transferred linearly in 
       ! the next step. if they dip slightly above or below their target allometries, 
       ! its no big deal.
       ! --------------------------------------------------------------------------------
       
       ! Also, if all stoichiometries are all 0.. then set grow_c_from_x to sufficiently
       ! high that it does not limit growth

       if (  (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,leaf_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,fnrt_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,sapw_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,store_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,struct_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ) < nearzero ) ) then
          grow_c_from_n = 1.0e10_r8
       end if
       if (  (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,leaf_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,fnrt_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,sapw_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,store_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,struct_organ) < nearzero ) .and. &
             (EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ) < nearzero ) ) then
          grow_c_from_p = 1.0e10_r8
       end if

       if(grow_c_from_c > nearzero) then
          carbon_gstature = carbon_gain * min(grow_c_from_c, grow_c_from_n, grow_c_from_p)/grow_c_from_c
       else
          write(fates_log(),*) 'Somehow grow_c_from_c is near zero',grow_c_from_c
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if


       if(carbon_gstature > nearzero) then
          
          ! Initialize the adaptive integrator arrays and flags
          ! -----------------------------------------------------------------------------------
          ierr             = 1
          totalC           = carbon_gstature
          nsteps           = 0
          
          ! Transfer carbon variables into the integrator
          do i = 1, num_organs
             i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
             ! Carbon in the actual tissue
             state_array(i) = sum(this%variables(i_cvar)%val(:))
             ! Respired carbon (should already be zero'd)
             state_array(i+num_organs) = 0.0_r8
          end do
          
          ! Transfer DBH into the integrator
          state_mask(intgr_dbh_id)       = .true.
          state_array(intgr_dbh_id)      = dbh
          
          if(ODESolve == 2) then
             this%ode_opt_step = totalC
          end if
          
          do while( ierr .ne. 0 )
             
             deltaC = min(totalC,this%ode_opt_step)
             if(ODESolve == 1) then
                call RKF45(AllomCNPGrowthDeriv,state_array,state_mask,deltaC,totalC, &
                     max_trunc_error,intgr_params,state_array_out,this%ode_opt_step,step_pass)
                
             elseif(ODESolve == 2) then
                call Euler(AllomCNPGrowthDeriv,state_array,state_mask,deltaC,totalC,intgr_params,state_array_out)
                !  step_pass = .true.
                call CheckIntegratedAllometries(state_array_out(intgr_dbh_id),ipft,canopy_trim,  &
                     state_array_out(intgr_leaf_c_id), state_array_out(intgr_fnrt_c_id), state_array_out(intgr_sapw_c_id), &
                     state_array_out(intgr_store_c_id), state_array_out(intgr_struct_c_id), &
                     state_mask(intgr_leaf_c_id), state_mask(intgr_fnrt_c_id), state_mask(intgr_sapw_c_id), &
                     state_mask(intgr_store_c_id),state_mask(intgr_struct_c_id),  max_trunc_error, step_pass)
                if(step_pass)  then
                   this%ode_opt_step = deltaC
                else
                   this%ode_opt_step = 0.5*deltaC
                end if
             else
                write(fates_log(),*) 'An integrator was chosen that DNE'
                write(fates_log(),*) 'ODESolve = ',ODESolve
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             nsteps = nsteps + 1
             
             if (step_pass) then ! If true, then step is accepted
                totalC         = totalC - deltaC
                state_array(:) = state_array_out(:)
             end if
             
             if(nsteps > max_substeps ) then
                write(fates_log(),*) 'CNP Plant Growth Integrator could not find'
                write(fates_log(),*) 'a solution in less than ',max_substeps,' tries'
                write(fates_log(),*) 'Aborting'
                write(fates_log(),*) 'carbon_balance',carbon_gain
                write(fates_log(),*) 'deltaC',deltaC
                write(fates_log(),*) 'totalC',totalC
                
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             ! TotalC should eventually be whittled down to near zero
             ! At that point, update the actual states
             ! --------------------------------------------------------------------------------
             if( (totalC < calloc_abs_error) .and. (step_pass) )then
                
                ierr           = 0

                ! Sum up the total flux predicted by the integrator,
                ! which SHOULD be carbon_gstature, except
                ! for integration errors.  To make carbon
                ! perfectly preserved, we calculate this bias
                ! and make a linear (proportional) correction to all pools.
                
                sum_c_flux = 0.0_r8
                do i = 1, num_organs
                   i_cvar = sp_organ_map(organ_list(i),carbon12_element)
                   if(state_mask(i_cvar))then
                      ! Add in the carbon flux
                      sum_c_flux = sum_c_flux + (state_array(i_cvar) - sum(this%variables(i_cvar)%val(:)))
                      ! Add in the growth respiration flux
                      sum_c_flux = sum_c_flux + state_array(i_cvar+num_organs)
                   end if
                end do
                
                c_flux_adj = carbon_gstature/sum_c_flux
                
                sum_n_demand = 0.0_r8
                sum_p_demand = 0.0_r8
                
                do i = 1, num_organs

                   i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
                   i_nvar   = sp_organ_map(organ_list(i),nitrogen_element)
                   i_pvar   = sp_organ_map(organ_list(i),phosphorus_element)
                   
                   i_gorgan = organ_list(i)
                   
                   ! Calculate adjusted flux
                   c_flux   = (state_array(i_cvar) - sum(this%variables(i_cvar)%val(:)))*c_flux_adj
                   gr_flux  = state_array(i_cvar+num_organs)*c_flux_adj
                   
                   ! update the carbon pool (in all pools flux goes into the first pool)
                   this%variables(i_cvar)%val(icd) = this%variables(i_cvar)%val(icd) + c_flux
                   
                   ! Track the growth respiration
                   growth_resp(i_cvar) =  growth_resp(i_cvar) + gr_flux
                   
                   ! update the nitrogen target and demand for the MINIMUM stoichiometry
                   ! Calculate this on the growing bin only

                   v_target(i_nvar) = &
                        this%variables(i_cvar)%val(icd)*EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,i_gorgan)
                   
                   v_demand(i_nvar) = &
                        max(0.0_r8,v_target(i_nvar) - this%variables(i_nvar)%val(icd))
                   
                   ! Update the phosphorus target and demand for the MINIMUM stoichiometry
                   ! Calculate this on the growing bin only

                   v_target(i_pvar) = &
                        this%variables(i_cvar)%val(icd)*EDPftvarcon_inst%prt_phos_stoich_p1(ipft,i_gorgan)
                   
                   v_demand(i_pvar) = &
                        max(0.0_r8,v_target(i_pvar) - this%variables(i_pvar)%val(icd))
                   
                   ! Sum up the total demands for each N&P for all pools
                   sum_n_demand = sum_n_demand + v_demand(i_nvar)
                   sum_p_demand = sum_p_demand + v_demand(i_pvar)
                   
                end do
                
                ! Determine the total flux from N&P pool
                sum_n_flux = min(nitrogen_gain, sum_n_demand) 
                sum_p_flux = min(phosphorus_gain, sum_p_demand)
                
                do i = 1, num_organs

                   i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
                   i_nvar   = sp_organ_map(organ_list(i),nitrogen_element)
                   i_pvar   = sp_organ_map(organ_list(i),phosphorus_element)
                   i_gorgan = organ_list(i)
                   
                   ! Transport nitrogren flux into the pool
                   if(sum_n_demand>nearzero) this%variables(i_nvar)%val(icd) = this%variables(i_nvar)%val(icd) + &
                        sum_n_flux * v_demand(i_nvar)  / sum_n_demand
                   
                   ! Transport phosphorus flux into the pool
                   if(sum_p_demand>nearzero) this%variables(i_pvar)%val(icd) = this%variables(i_pvar)%val(icd) + &
                        sum_p_flux * v_demand(i_pvar)  / sum_p_demand
                   
                end do
                
                ! --------------------------------------------------------------------------
                ! Remove fluxes from the C&N&P gains pool
                ! --------------------------------------------------------------------------
                
                carbon_gain      = carbon_gain - carbon_gstature
                nitrogen_gain    = nitrogen_gain - sum_n_flux
                phosphorus_gain = phosphorus_gain - sum_p_flux
                

                ! --------------------------------------------------------------------------
                ! Update the non-pool integrated variables
                ! --------------------------------------------------------------------------
                
                dbh                = state_array(intgr_dbh_id)
   
                
             end if
          end do
          
       end if
       
    end if
    
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    !
    ! If nutrients are still available, then we can bump up the values in the pools
    !  towards the OPTIMAL target values.
    !
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

    sum_n_demand = 0.0_r8
    sum_p_demand = 0.0_r8
    do i = 1, num_organs
       
       i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
       i_nvar   = sp_organ_map(organ_list(i),nitrogen_element)
       i_pvar   = sp_organ_map(organ_list(i),phosphorus_element)
       i_gorgan = organ_list(i)
       
       ! This is the optimal nitrogen target and demand
       v_target(i_nvar) = this%variables(i_cvar)%val(icd) * &
            EDPftvarcon_inst%prt_nitr_stoich_p2(ipft,i_gorgan) 
       
       v_demand(i_nvar) = max(0.0_r8, v_target(i_nvar) - this%variables(i_nvar)%val(icd))

       sum_n_demand = sum_n_demand + v_demand(i_nvar)


       ! This is the optimal phosphorus target and demand
       v_target(i_pvar) = this%variables(i_cvar)%val(icd) * &
            EDPftvarcon_inst%prt_phos_stoich_p2(ipft,i_gorgan) 

       v_demand(i_pvar) = max(0.0_r8, v_target(i_pvar) - this%variables(i_pvar)%val(icd))
       
       sum_p_demand = sum_p_demand + v_demand(i_pvar)
       
    end do

    sum_n_flux =  min(nitrogen_gain,sum_n_demand)
    sum_p_flux =  min(phosphorus_gain,sum_p_demand)

    do i = 1, num_organs
       
       i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
       i_nvar   = sp_organ_map(organ_list(i),nitrogen_element)
       i_pvar   = sp_organ_map(organ_list(i),phosphorus_element)

       if(sum_n_demand>nearzero) then
          n_flux = sum_n_flux * v_demand(i_nvar) / sum_n_demand
       else
          n_flux = 0.0_r8
       end if
       if(sum_p_demand>nearzero) then
          p_flux = sum_p_flux * v_demand(i_pvar) / sum_p_demand
       else
          p_flux = 0.0_r8
       end if
       
       this%variables(i_nvar)%val(icd) = this%variables(i_nvar)%val(icd) + n_flux
       this%variables(i_pvar)%val(icd) = this%variables(i_pvar)%val(icd) + p_flux

    end do
    
    nitrogen_gain    = nitrogen_gain - sum_n_flux
    phosphorus_gain = phosphorus_gain - sum_p_flux

    ! -----------------------------------------------------------------------------------
    ! If carbon is still available, lets cram some into storage overflow
    ! We will do this last, because we wanted the non-overflow storage
    ! value to draw minimum and optimal nutrient fluxes
    ! -----------------------------------------------------------------------------------
    
    if(carbon_gain>nearzero) then

       ! Update carbon based allometric targets
       call bstore_allom(dbh,ipft,canopy_trim, store_c_target, store_dcdd_target)
       
       ! Estimate the overflow
       store_c_target = store_c_target * (1.0_r8 + store_overflow_frac)
       
       store_c_demand = max(0.0,(store_c_target - &
            this%variables(store_c_id)%val(icd))*(1.0 + r_g(store_organ)))
       
       store_c_flux   = min(carbon_gain,store_c_demand)

       ! Transfer excess carbon into storage overflow
       this%variables(store_c_id)%val(icd) = &
            this%variables(store_c_id)%val(icd) + store_c_flux / (1.0_r8 + r_g(store_organ))

       growth_resp(store_c_id) =  growth_resp(store_c_id) + &
             store_c_flux * r_g(store_organ) / (1.0_r8 + r_g(store_organ))
       
       carbon_gain = carbon_gain - store_c_flux
       

    end if

    ! Sum up growth respiration
    total_growth_respiration = 0.0_r8
    do i = 1, num_organs
       i_cvar   = sp_organ_map(organ_list(i),carbon12_element)
       total_growth_respiration = total_growth_respiration + growth_resp(i_cvar)
    end do
    
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    !
    ! Figure out what to do with excess carbon and nutrients
    !
    ! 1) excude through roots cap at 0 to flush out imprecisions
    ! 
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    

    this%bc_out(acnp_bc_out_id_rootcexude)%rval = max(0.0_r8,carbon_gain)
    this%bc_out(acnp_bc_out_id_rootnexude)%rval = max(0.0_r8,nitrogen_gain)
    this%bc_out(acnp_bc_out_id_rootpexude)%rval = max(0.0_r8,phosphorus_gain)
    this%bc_out(acnp_bc_out_id_growresp)%rval   = total_growth_respiration

    carbon_gain      = 0.0_r8
    nitrogen_gain    = 0.0_r8
    phosphorus_gain = 0.0_r8


    ! Update the diagnostic on daily rate of change

    do i_age = 1,nleafage
       this%variables(leaf_c_id)%net_alloc(i_age) = &
            this%variables(leaf_c_id)%net_alloc(i_age) + &
            (leaf_c(i_age) - leaf_c0(i_age))

       this%variables(leaf_n_id)%net_alloc(i_age) = &
            this%variables(leaf_n_id)%net_alloc(i_age) + &
            (leaf_n(i_age) - leaf_n0(i_age))
       
       this%variables(leaf_p_id)%net_alloc(i_age) = &
            this%variables(leaf_p_id)%net_alloc(i_age) + &
            (leaf_p(i_age) - leaf_p0(i_age))

    end do

    this%variables(fnrt_c_id)%net_alloc(icd) = &
         this%variables(fnrt_c_id)%net_alloc(icd) + (fnrt_c - fnrt_c0)
    
    this%variables(sapw_c_id)%net_alloc(icd) = &
         this%variables(sapw_c_id)%net_alloc(icd) + (sapw_c - sapw_c0)
    
    this%variables(store_c_id)%net_alloc(icd) = &
         this%variables(store_c_id)%net_alloc(icd) + (store_c - store_c0)
    
    this%variables(repro_c_id)%net_alloc(icd) = &
         this%variables(repro_c_id)%net_alloc(icd) + (repro_c - repro_c0)
    
    this%variables(struct_c_id)%net_alloc(icd) = &
         this%variables(struct_c_id)%net_alloc(icd) + (struct_c - struct_c0)

    this%variables(fnrt_n_id)%net_alloc(icd) = &
         this%variables(fnrt_n_id)%net_alloc(icd) + (fnrt_n - fnrt_n0)
    
    this%variables(sapw_n_id)%net_alloc(icd) = &
         this%variables(sapw_n_id)%net_alloc(icd) + (sapw_n - sapw_n0)
    
    this%variables(store_n_id)%net_alloc(icd) = &
         this%variables(store_n_id)%net_alloc(icd) + (store_n - store_n0)
    
    this%variables(repro_n_id)%net_alloc(icd) = &
         this%variables(repro_n_id)%net_alloc(icd) + (repro_n - repro_n0)
    
    this%variables(struct_n_id)%net_alloc(icd) = &
         this%variables(struct_n_id)%net_alloc(icd) + (struct_n - struct_n0)

    this%variables(fnrt_p_id)%net_alloc(icd) = &
         this%variables(fnrt_p_id)%net_alloc(icd) + (fnrt_p - fnrt_p0)
    
    this%variables(sapw_p_id)%net_alloc(icd) = &
         this%variables(sapw_p_id)%net_alloc(icd) + (sapw_p - sapw_p0)
    
    this%variables(store_p_id)%net_alloc(icd) = &
         this%variables(store_p_id)%net_alloc(icd) + (store_p - store_p0)
    
    this%variables(repro_p_id)%net_alloc(icd) = &
         this%variables(repro_p_id)%net_alloc(icd) + (repro_p - repro_p0)
    
    this%variables(struct_p_id)%net_alloc(icd) = &
         this%variables(struct_p_id)%net_alloc(icd) + (struct_p - struct_p0)
          
  end associate

  return
  end subroutine DailyPRTAllometricCNP
  
  ! =====================================================================================

  subroutine FastPRTAllometricCNP(this)
      
    implicit none
    class(cnp_allom_prt_vartypes) :: this     ! this class
    
    ! This routine does nothing, becaues
    ! we currently don't have any fast-timestep processes
    ! Think of this as a stub.
    
    return
  end subroutine FastPRTAllometricCNP

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


      associate( dbh         => l_state_array(intgr_dbh_id),      &
                 leaf_c      => l_state_array(intgr_leaf_c_id),   &
                 fnrt_c      => l_state_array(intgr_fnrt_c_id),   &
                 sapw_c      => l_state_array(intgr_sapw_c_id),   &
                 store_c     => l_state_array(intgr_store_c_id),  &
                 struct_c    => l_state_array(intgr_struct_c_id), &
                 repro_c     => l_state_array(intgr_repro_c_id),  &
                 mask_dbh    => l_state_mask(intgr_dbh_id),       &
                 mask_leaf   => l_state_mask(intgr_leaf_c_id),    &  
                 mask_fnrt   => l_state_mask(intgr_fnrt_c_id),    &
                 mask_sapw   => l_state_mask(intgr_sapw_c_id),    &
                 mask_store  => l_state_mask(intgr_store_c_id),   &
                 mask_struct => l_state_mask(intgr_struct_c_id),  &
                 mask_repro  => l_state_mask(intgr_repro_c_id) )


        canopy_trim = intgr_params(acnp_bc_in_id_ctrim)
        ipft        = int(intgr_params(acnp_bc_in_id_pft))

        if(dbh>huge(dbh)) then
           print*,"BIG D IN DERIV:",dbh
           stop
        end if

        call bleaf(dbh,ipft,canopy_trim,leaf_c_target,leaf_dcdd_target)
        call bfineroot(dbh,ipft,canopy_trim,fnrt_c_target,fnrt_dcdd_target)
        call bsap_allom(dbh,ipft,canopy_trim,sapw_area,sapw_c_target,sapw_dcdd_target)
        call bagw_allom(dbh,ipft,agw_c_target,agw_dcdd_target)
        call bbgw_allom(dbh,ipft,bgw_c_target,bgw_dcdd_target)
        call bdead_allom(agw_c_target,bgw_c_target, sapw_c_target, ipft, struct_c_target, &
                         agw_dcdd_target, bgw_dcdd_target, sapw_dcdd_target, struct_dcdd_target) 
        call bstore_allom(dbh,ipft,canopy_trim,store_c_target,store_dcdd_target)
        
        ! fraction of carbon going towards reproduction
        if (dbh <= EDPftvarcon_inst%dbh_repro_threshold(ipft)) then
           repro_fraction = EDPftvarcon_inst%seed_alloc(ipft)
        else
           repro_fraction = EDPftvarcon_inst%seed_alloc(ipft) + EDPftvarcon_inst%seed_alloc_mature(ipft)
        end if

        total_dcostdd = struct_dcdd_target * (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,struct_organ))

        if (mask_leaf) then
           total_dcostdd = total_dcostdd + leaf_dcdd_target * (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,leaf_organ))
        end if
        if (mask_fnrt) then
           total_dcostdd = total_dcostdd + fnrt_dcdd_target * (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,fnrt_organ))
        end if
        if (mask_sapw) then
           total_dcostdd = total_dcostdd + sapw_dcdd_target * (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,sapw_organ))
        end if
        if (mask_store) then
           total_dcostdd = total_dcostdd + store_dcdd_target * (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,store_organ))
        end if
        
        dCdx(:) = 0.0_r8

        ! It is possible that with some asymptotic, or hard
        ! capped allometries, that all growth rates reach zero.
        ! In this case, if there is carbon, give it to reproduction

        if(total_dcostdd > nearzero) then
           
           dCdx(intgr_struct_c_id) = struct_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
           dCdx(intgr_struct_gr_id) = struct_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction) * &
                EDPftvarcon_inst%prt_grperc_organ(ipft,struct_organ)
           
           if (mask_leaf) then
              dCdx(intgr_leaf_c_id)  = leaf_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
              dCdx(intgr_leaf_gr_id) = (leaf_dcdd_target/total_dcostdd) *  (1.0_r8 - repro_fraction) * &
                   EDPftvarcon_inst%prt_grperc_organ(ipft,leaf_organ)
           end if
           
           if (mask_fnrt) then
              dCdx(intgr_fnrt_c_id)  = fnrt_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
              dCdx(intgr_fnrt_gr_id) = (fnrt_dcdd_target/total_dcostdd) * (1.0_r8 - repro_fraction) * &
                   EDPftvarcon_inst%prt_grperc_organ(ipft,fnrt_organ)
           end if
           
           if (mask_sapw) then
              dCdx(intgr_sapw_c_id)  = sapw_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
              dCdx(intgr_sapw_gr_id) = (sapw_dcdd_target/total_dcostdd) * (1.0_r8 - repro_fraction) * &
                   EDPftvarcon_inst%prt_grperc_organ(ipft,sapw_organ)
           end if
           
           if (mask_store) then
              dCdx(intgr_store_c_id)  = store_dcdd_target/total_dcostdd * (1.0_r8 - repro_fraction)
              dCdx(intgr_store_gr_id) = (store_dcdd_target/total_dcostdd) * (1.0_r8 - repro_fraction) * &
                   EDPftvarcon_inst%prt_grperc_organ(ipft,store_organ)
           end if
           
           if (mask_repro) then
              dCdx(intgr_repro_c_id)  = repro_fraction / (1.0_r8 + EDPftvarcon_inst%prt_grperc_organ(ipft,repro_organ) )
              dCdx(intgr_repro_gr_id) = repro_fraction * EDPftvarcon_inst%prt_grperc_organ(ipft,repro_organ) / &
                   (1.0_r8 +  EDPftvarcon_inst%prt_grperc_organ(ipft,repro_organ) )
           end if
           
           dCdx(intgr_dbh_id)      = (1.0_r8/total_dcostdd)*(1.0_r8 - repro_fraction)

           
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

   subroutine FastPRTAC(this)
      
      implicit none
      class(cnp_allom_prt_vartypes) :: this     ! this class
      
      ! This routine does nothing, because in the carbon only allometric RT model
      ! we currently don't have any fast-timestep processes
      ! Think of this as a stub.
      
      



      return
    end subroutine FastPRTAC


end module PRTAllometricCNPMod
  
