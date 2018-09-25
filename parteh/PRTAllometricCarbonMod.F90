module PRTAllometricCarbonMod

   ! ------------------------------------------------------------------------------------
   !
   ! This module contains all of the specific functions and types for
   ! Plant Allocation and Reactive Transport Extensible Hypotheses (PARTEH)
   ! CARBON only, allometric growth hypothesis
   ! 
   ! Ryan Knox Apr 2018
   !
   ! TO-DO: THE MAPPING TABLES SHOULD BE PROTECTED STATUS. TEST ADDING THIS AFTER 1ST
   ! SUCCESFULL RUN
   !
   ! ------------------------------------------------------------------------------------

  use PRTGenericMod , only  : prt_instance_type
  use PRTGenericMod , only  : prt_vartype
  use PRTGenericMod , only  : prt_vartypes
  use PRTGenericMod , only  : carbon12_species
  use PRTGenericMod , only  : leaf_organ
  use PRTGenericMod , only  : fnrt_organ
  use PRTGenericMod , only  : sapw_organ
  use PRTGenericMod , only  : store_organ
  use PRTGenericMod , only  : repro_organ
  use PRTGenericMod , only  : struct_organ

  use PRTLossFluxesMod    , only : PRTMaintTurnover

  use FatesInterfaceMod   , only : hlm_freq_day
  use FatesAllometryMod   , only : bleaf
  use FatesAllometryMod   , only : bsap_allom
  use FatesAllometryMod   , only : bfineroot
  use FatesAllometryMod   , only : bstore_allom
  use FatesAllometryMod   , only : bdead_allom
  use FatesAllometryMod   , only : bbgw_allom
  use FatesAllometryMod   , only : bagw_allom
  use FatesAllometryMod   , only : h_allom
  use FatesAllometryMod   , only : CheckIntegratedAllometries
  use FatesAllometryMod   , only : StructureResetOfDH

  use FatesGlobals        , only : endrun => fates_endrun
  use FatesGlobals        , only : fates_log
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use FatesConstantsMod   , only : r8 => fates_r8
  use FatesConstantsMod   , only : i4 => fates_int
  use FatesIntegratorsMod , only : RKF45
  use FatesIntegratorsMod , only : Euler
  use EDPftvarcon         , only : EDPftvarcon_inst
  use FatesConstantsMod   , only : calloc_abs_error
  use FatesConstantsMod   , only : nearzero
  use FatesConstantsMod   , only : itrue
  !  use PARTEHUtilitiesMod  , only : MaintenanceTurnover

  implicit none
  private

  ! -------------------------------------------------------------------------------------
  !
  ! Define the state variables for this specific hypothesis. Give them units and define 
  ! the indices that correspond with the generic classifications of PRT variables
  !
  ! -------------------------------------------------------------------------------------

  integer, parameter :: leaf_c_id   = 1
  integer, parameter :: fnrt_c_id   = 2
  integer, parameter :: sapw_c_id   = 3
  integer, parameter :: store_c_id  = 4
  integer, parameter :: repro_c_id  = 5
  integer, parameter :: struct_c_id = 6
  integer, parameter :: ac_num_vars = 6  ! Number of PRT variables
  
  integer, parameter :: dbh_id          = 7  ! This is just used for the integrator
  integer, parameter :: n_integration_vars = 7


  ! -------------------------------------------------------------------------------------
  ! Boundary Conditions
  ! -------------------------------------------------------------------------------------
  ! Input Boundary Indices (These are public, and therefore
  !       each boundary condition across all modules must
  !       have a unique name !!!!)
  ! -------------------------------------------------------------------------------------

  integer, public, parameter :: ac_bc_inout_id_dbh   = 1   ! Plant DBH
  integer, public, parameter :: ac_bc_inout_id_netdc = 2   ! Index for the net daily C input BC
  integer, parameter         :: num_bc_inout         = 2


  integer, public, parameter :: ac_bc_in_id_pft   = 1   ! Index for the PFT input BC
  integer, public, parameter :: ac_bc_in_id_ctrim = 2   ! Index for the canopy trim function
  integer, parameter         :: num_bc_in         = 2


  ! -------------------------------------------------------------------------------------
  ! Define the size of the coorindate vector.  For this hypothesis, there is only
  ! one pool per each species x organ combination.
  ! -------------------------------------------------------------------------------------
  integer, parameter         :: icd               = 1   ! Only 1 coordinate per variable


  ! -------------------------------------------------------------------------------------
  ! This is the core type that holds this specific
  ! plant reactive transport (PRT) module
  ! -------------------------------------------------------------------------------------

  type callom_prt_vartype

     real(r8) :: allom_deficit   ! Deficit of plant WRT allometric target

  end type callom_prt_vartype



  type, public, extends(prt_vartypes) :: callom_prt_vartypes

     type(callom_prt_vartype),allocatable :: aux_variables(:)

   contains

     procedure :: DailyPRT     => DailyPRTAC
     procedure :: FastPRT      => FastPRTAC
     procedure :: InitAllocate => InitAllocateAC

   end type callom_prt_vartypes
   
   ! ------------------------------------------------------------------------------------
   !
   ! This next class is an extention of the base instance that maps state variables
   !      to the outside model.
   !
   ! ------------------------------------------------------------------------------------
   
   character(len=*), parameter, private :: sourcefile = __FILE__


   ! This is the instance of the mapping table and variable definitions
   ! this is only allocated once per node
   class(prt_instance_type), target, allocatable :: prt_instance_ac


   public :: InitPRTInstanceAC


contains
  
 
  subroutine InitPRTInstanceAC()

     ! ----------------------------------------------------------------------------------
     ! Initialize and populate the general mapping table that 
     ! organizes the specific variables in this module to
     ! pre-ordained groups, so they can be used to inform
     ! the rest of the model
     ! -----------------------------------------------------------------------------------

     allocate(prt_instance_ac)
     allocate(prt_instance_ac%state_descriptor(ac_num_vars))

     prt_instance_ac%hyp_name = 'Allometric Carbon Only'
     
     call prt_instance_ac%ZeroInstance()

     ! Populate the array
     ! This is a carbon only scheme, no isotopes, so should be simple
     ! The "indices array" max not exceed max_types_per_sp_organ
     ! If that array limit is not large enough for new hypothesis
     ! simply increase it.  It will not use much memory or increase loop sizes
     

     call prt_instance_ac%InitInstance(leaf_c_id,"Leaf Carbon","leaf_c",leaf_organ,carbon12_species)
     call prt_instance_ac%InitInstance(fnrt_c_id,"Fine Root Carbon","fnrt_c",fnrt_organ,carbon12_species)
     call prt_instance_ac%InitInstance(sapw_c_id,"Sapwood Carbon","sapw_c",sapw_organ,carbon12_species)
     call prt_instance_ac%InitInstance(store_c_id,"Storage Carbon","store_c",store_organ,carbon12_species)
     call prt_instance_ac%InitInstance(struct_c_id,"Structural Carbon","struct_c",struct_organ,carbon12_species)
     call prt_instance_ac%InitInstance(repro_c_id,"Reproductive Carbon","repro_c",repro_organ,carbon12_species)
     
     return
   end subroutine InitPRTInstanceAC


   ! =====================================================================================
   

   subroutine InitAllocateAC(this)
    
     ! ----------------------------------------------------------------------------------
     ! This initialization is called everytime a plant/cohort
     ! is newly recruited.  This simply sets-up, allocates
     ! and sets some initialization values
     ! ----------------------------------------------------------------------------------

     class(callom_prt_vartypes)                :: this     ! this class
        
     integer :: ivar


     ! Set the instance pointer to the correct instance
     ! ----------------------------------------------------------------------------------

     this%prt_instance => prt_instance_ac

     
     ! Allocate the boundar condition arrays and flush them to no-data flags
     ! ----------------------------------------------------------------------------------

     allocate(this%bc_in(num_bc_in))
     allocate(this%bc_inout(num_bc_inout))

     
     ! Allocate the state variables
     allocate(this%variables(ac_num_vars))
     
     do ivar = 1, ac_num_vars

        this%variables(ivar)%num_pos = icd
        allocate(this%variables(ivar)%val(icd))
        allocate(this%variables(ivar)%val0(icd))
        allocate(this%variables(ivar)%turnover(icd))
        allocate(this%variables(ivar)%dvaldt(icd))

     end do

     ! Initialize the optimum step size as very large.

     this%ode_opt_step = 1e6_r8

     
     return
   end subroutine InitAllocateAC
  

  ! =====================================================================================
  

  subroutine DailyPRTAC(this)

    
    ! The class is the only argument, input and output bc's are globals
    class(callom_prt_vartypes)   :: this          ! this class

    ! -----------------------------------------------------------------------------------
    ! These are local copies of the in/out boundary condition structure
    ! -----------------------------------------------------------------------------------

    real(r8),pointer :: dbh               ! Diameter at breast height [cm]
                                          ! this local will point to both in and out bc's
    real(r8),pointer :: carbon_balance    ! Daily carbon balance for this cohort [kgC]

    ! These are local copies of the input only boundary conditions
    real(r8) :: canopy_trim           ! The canopy trimming function [0-1]
    integer  :: ipft                  ! Plant Functional Type index

    ! -----------------------------------------------------------------------------------
    ! Local copies of output boundary conditions
    ! -----------------------------------------------------------------------------------
    
    real(r8) :: target_leaf_c     ! target leaf carbon [kgC]
    real(r8) :: target_fnrt_c     ! target fine-root carbon [kgC]
    real(r8) :: target_sapw_c     ! target sapwood carbon [kgC]
    real(r8) :: target_store_c    ! target storage carbon [kgC]
    real(r8) :: target_agw_c      ! target above ground carbon in woody tissues [kgC]
    real(r8) :: target_bgw_c      ! target below ground carbon in woody tissues [kgC]
    real(r8) :: target_struct_c   ! target structural carbon [kgC]

    real(r8) :: leaf_below_target     ! fineroot biomass below target amount [kgC]
    real(r8) :: fnrt_below_target     ! fineroot biomass below target amount [kgC]
    real(r8) :: sapw_below_target      ! sapwood biomass below target amount [kgC]
    real(r8) :: store_below_target    ! storage biomass below target amount [kgC]
    real(r8) :: struct_below_target   ! dead (structural) biomass below target amount [kgC]
    real(r8) :: total_below_target    ! total biomass below the allometric target [kgC]

    real(r8) :: flux_adj              ! adjustment made to growth flux term to minimize error [kgC]
    real(r8) :: store_target_fraction ! ratio between storage and leaf biomass when on allometry [kgC]

    real(r8) :: leaf_c_demand  ! leaf carbon that is demanded to replace maintenance turnover [kgC]
    real(r8) :: fnrt_c_demand  ! fineroot carbon that is demanded to replace 
                                      ! maintenance turnover [kgC]
    real(r8) :: total_c_demand ! total carbon that is demanded to replace maintenance turnover [kgC]
    real(r8) :: sapw_area             ! dummy sapwood area
    logical  :: step_pass             ! Did the integration step pass?

    real(r8) :: leaf_c_flux
    real(r8) :: fnrt_c_flux             
    real(r8) :: sapw_c_flux
    real(r8) :: store_c_flux
    real(r8) :: repro_c_flux
    real(r8) :: struct_c_flux

    logical  :: grow_leaf             ! Are leaves at allometric target and should be grown?
    logical  :: grow_fnrt             ! Are fine-roots at allometric target and should be grown?
    logical  :: grow_sapw             ! Is sapwood at allometric target and should be grown?
    logical  :: grow_store            ! Is storage at allometric target and should be grown?

                                      ! integrator variables
    real(r8) :: deltaC                ! trial value for substep
    integer  :: ierr                  ! error flag for allometric growth step
    integer  :: nsteps                ! number of sub-steps
    integer  :: istep                 ! current substep index
    real(r8) :: totalC                ! total carbon allocated over alometric growth step
    real(r8) :: hite_out              ! dummy height variable

    integer  :: i_var                 ! local index for iterating state variables


    ! Integegrator variables

    real(r8),dimension(n_integration_vars) :: c_pool       ! Vector of carbon pools passed to integrator
    real(r8),dimension(n_integration_vars) :: c_pool_out   ! Vector of carbon pools passed back from integrator
    logical,dimension(n_integration_vars)  :: c_mask       ! Mask of active pools during integration

    real(r8), parameter :: cbal_prec = 1.0e-15_r8   ! Desired precision in carbon balance
    integer , parameter :: max_substeps = 300       ! Maximum allowable iterations
    real(r8), parameter :: max_trunc_error = 1.0_r8 ! Maximum allowable truncation error
    integer,  parameter :: ODESolve = 2             ! 1=RKF45,  2=Euler

    real(r8) ::  intgr_params(num_bc_in)


    ! This is a local array containing the boundary conditions
    ! we need this (for now at least) because the integration layer needs things
    ! packed into simple types


    ! This array is used to hold parameters that must be passed through
    ! a generic integrator to the derivative functions

    associate( & 
          leaf_c   => this%variables(leaf_c_id)%val(icd), &
          fnrt_c   => this%variables(fnrt_c_id)%val(icd), &
          sapw_c   => this%variables(sapw_c_id)%val(icd), &
          store_c  => this%variables(store_c_id)%val(icd), &
          repro_c  => this%variables(repro_c_id)%val(icd), &
          struct_c => this%variables(struct_c_id)%val(icd))


    ! ===================================================================================
    !
    ! !!!! CALCULATIONS THAT SHOULD NOW BE OUTSIDE OF THIS ROUTINE !!!!
    ! WE USED TO SET THE "ISNEW" FLAG HERE
    ! MAKE SURE THAT IT IS SET AFTER THIS ROUTINE IS CALLED
    ! IT SHOULD BE CALLED FOR ANY INSTANCE, NOT JUST THIS
    ! currentSite%flux_in = currentSite%flux_in + currentCohort%npp_acc * currentCohort%n
    ! 
    ! ===================================================================================

    ! Copy the boundary conditions into readable local variables  
    ! We don't use pointers, because inputs should be intent in only

    dbh                             => this%bc_inout(ac_bc_inout_id_dbh)%rval
    carbon_balance                  => this%bc_inout(ac_bc_inout_id_netdc)%rval

    canopy_trim                     = this%bc_in(ac_bc_in_id_ctrim)%rval
    ipft                            = this%bc_in(ac_bc_in_id_pft)%ival

    intgr_params(:)                 = -9.9e32_r8
    intgr_params(ac_bc_in_id_ctrim) = this%bc_in(ac_bc_in_id_ctrim)%rval
    intgr_params(ac_bc_in_id_pft)   = real(this%bc_in(ac_bc_in_id_pft)%ival)


    ! -----------------------------------------------------------------------------------
    ! I. Calculate target size of the biomass compartment for a given dbh.   
    ! -----------------------------------------------------------------------------------

    
    ! Target sapwood biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
    !call bsap_allom(dbh,ipft,canopy_trim,sapw_area,target_sapw_c)
    call bsap_allom(dbh,ipft,canopy_trim,target_sapw_c)
    

    ! Target total above ground deriv. biomass in woody/fibrous tissues  [kgC, kgC/cm]
    call bagw_allom(dbh,ipft,target_agw_c)

    ! Target total below ground deriv. biomass in woody/fibrous tissues [kgC, kgC/cm] 
    call bbgw_allom(dbh,ipft,target_bgw_c)

    ! Target total dead (structrual) biomass and deriv. [kgC, kgC/cm]
    call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)

    


    ! ------------------------------------------------------------------------------------
    ! If structure is larger than target, then we need to correct some integration errors
    ! by slightly increasing dbh to match it.
    ! For grasses, if leaf biomass is larger than target, then we reset dbh to match
    ! -----------------------------------------------------------------------------------
    if( (( struct_c - target_struct_c ) > calloc_abs_error) .and. &
          (EDPftvarcon_inst%woody(ipft) == itrue) ) then

       call StructureResetOfDH( struct_c, ipft, &
             canopy_trim, dbh, hite_out )

       ! Target sapwood biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
       !call bsap_allom(dbh,ipft,canopy_trim,sapw_area,target_sapw_c)
       call bsap_allom(dbh,ipft,canopy_trim,target_sapw_c)

       ! Target total above ground deriv. biomass in woody/fibrous tissues  [kgC, kgC/cm]
       call bagw_allom(dbh,ipft,target_agw_c)
       
       ! Target total below ground deriv. biomass in woody/fibrous tissues [kgC, kgC/cm] 
       call bbgw_allom(dbh,ipft,target_bgw_c)
       
       ! Target total dead (structrual) biomass and deriv. [kgC, kgC/cm]
       call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)

    end if
    
    ! Target leaf biomass according to allometry and trimming
    call bleaf(dbh,ipft,canopy_trim,target_leaf_c)

    ! Target fine-root biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
    call bfineroot(dbh,ipft,canopy_trim,target_fnrt_c)

    ! Target storage carbon [kgC,kgC/cm]
    call bstore_allom(dbh,ipft,canopy_trim,target_store_c)

    
    ! -----------------------------------------------------------------------------------
    !  Set memory of the old state variables for comparison
    ! -----------------------------------------------------------------------------------

    do i_var = 1,ac_num_vars
       this%variables(i_var)%val0(icd) = this%variables(i_var)%val(icd)
    end do

    ! -----------------------------------------------------------------------------------
    ! II. Call maintenance turnover
    ! This will increment %turnover and decrease %val
    ! ----------------------------------------------------------------------------------

    call PRTMaintTurnover(this,ipft)

    ! -----------------------------------------------------------------------------------
    ! III.  Prioritize some amount of carbon to replace leaf/root turnover
    !         Make sure it isnt a negative payment, and either pay what is available
    !         or forcefully pay from storage. 
    ! -----------------------------------------------------------------------------------
    
    leaf_c_demand   = max(0.0_r8,(target_leaf_c - leaf_c))
    fnrt_c_demand   = max(0.0_r8,(target_fnrt_c - fnrt_c))

    total_c_demand = leaf_c_demand + fnrt_c_demand
    
    if (total_c_demand> nearzero) then

       ! If we are testing b4b, then we pay this even if we don't have the carbon
       ! Just don't pay so much carbon that storage+carbon_balance can't pay for it
       leaf_c_flux = min(leaf_c_demand, &
                         max(0.0_r8,(store_c+carbon_balance)* &
                         (leaf_c_demand/total_c_demand)))
       
       carbon_balance = carbon_balance - leaf_c_flux
       leaf_c         = leaf_c + leaf_c_flux

       ! If we are testing b4b, then we pay this even if we don't have the carbon
       fnrt_c_flux = min(fnrt_c_demand, &
                         max(0.0_r8, (store_c+carbon_balance)* &
                         (fnrt_c_demand/total_c_demand)))

       carbon_balance = carbon_balance - fnrt_c_flux
       fnrt_c         = fnrt_c + fnrt_c_flux

    end if

    ! -----------------------------------------------------------------------------------
    ! IV. if carbon balance is negative, re-coup the losses from storage
    !       if it is positive, give some love to storage carbon
    ! -----------------------------------------------------------------------------------

    if( carbon_balance < 0.0_r8 ) then
       
       store_c_flux           = carbon_balance
       carbon_balance         = carbon_balance - store_c_flux
       store_c                = store_c + store_c_flux

    else

       store_below_target     = max(target_store_c - store_c,0.0_r8)
       store_target_fraction  = max(0.0_r8, store_c/target_store_c )

       store_c_flux           = min(store_below_target,carbon_balance * &
                                max(exp(-1.*store_target_fraction**4._r8) - exp( -1.0_r8 ),0.0_r8))

       carbon_balance         = carbon_balance - store_c_flux
       store_c                = store_c + store_c_flux

    end if

    ! -----------------------------------------------------------------------------------
    ! V.  If carbon is still available, prioritize some allocation to replace
    !        the rest of the leaf/fineroot deficit
    !        carbon balance is guaranteed to be >=0 beyond this point
    ! -----------------------------------------------------------------------------------
    
    leaf_c_demand   = max(0.0_r8,(target_leaf_c - leaf_c))
    fnrt_c_demand   = max(0.0_r8,(target_fnrt_c - fnrt_c))

    total_c_demand = leaf_c_demand + fnrt_c_demand
    
    if( (carbon_balance > nearzero ) .and. (total_c_demand>nearzero)) then

       leaf_c_flux    = min(leaf_c_demand, &
                        carbon_balance*(leaf_c_demand/total_c_demand))
       carbon_balance = carbon_balance - leaf_c_flux
       leaf_c         = leaf_c + leaf_c_flux
       
       fnrt_c_flux    = min(fnrt_c_demand, &
                            carbon_balance*(fnrt_c_demand/total_c_demand))
       carbon_balance = carbon_balance - fnrt_c_flux
       fnrt_c         = fnrt_c + fnrt_c_flux

    end if

    ! -----------------------------------------------------------------------------------
    ! VI.  If carbon is still available, we try to push all live 
    !        pools back towards allometry. But only upwards, if fusion happened
    !        to generate some pools above allometric target, don't reduce the pool,
    !        just ignore it until the rest of the plant grows to meet it.
    ! -----------------------------------------------------------------------------------
    if( carbon_balance > nearzero ) then

       leaf_below_target  = max(target_leaf_c - leaf_c,0.0_r8)
       fnrt_below_target  = max(target_fnrt_c - fnrt_c,0.0_r8)
       sapw_below_target  = max(target_sapw_c - sapw_c,0.0_r8)
       store_below_target = max(target_store_c - store_c,0.0_r8)
       
       total_below_target = leaf_below_target + fnrt_below_target + &
                            sapw_below_target + store_below_target
    
       if ( total_below_target > nearzero ) then
          
          if( total_below_target > carbon_balance) then
             leaf_c_flux  = carbon_balance * leaf_below_target/total_below_target
             fnrt_c_flux = carbon_balance * fnrt_below_target/total_below_target
             sapw_c_flux  = carbon_balance * sapw_below_target/total_below_target
             store_c_flux = carbon_balance * store_below_target/total_below_target
          else
             leaf_c_flux  = leaf_below_target
             fnrt_c_flux = fnrt_below_target
             sapw_c_flux = sapw_below_target
             store_c_flux = store_below_target
          end if

          carbon_balance               = carbon_balance - leaf_c_flux
          leaf_c                       = leaf_c + leaf_c_flux
          
          carbon_balance               = carbon_balance - fnrt_c_flux
          fnrt_c                      = fnrt_c + fnrt_c_flux
          
          carbon_balance               = carbon_balance - sapw_c_flux
          sapw_c                       = sapw_c + sapw_c_flux
          
          carbon_balance               = carbon_balance - store_c_flux
          store_c                      = store_c  +  store_c_flux
          
       end if
    end if
    
    ! -----------------------------------------------------------------------------------
    ! VIII.  If carbon is still available, replenish the structural pool to get
    !           back on allometry
    ! -----------------------------------------------------------------------------------

    if( carbon_balance > nearzero ) then
    
       struct_below_target  = max(target_struct_c - struct_c ,0.0_r8)
    
       if ( struct_below_target > 0.0_r8) then
       
          struct_c_flux          = min(carbon_balance,struct_below_target)
          carbon_balance         = carbon_balance - struct_c_flux
          struct_c               = struct_c + struct_c_flux
          
       end if

    end if
    
    ! -----------------------------------------------------------------------------------
    ! IX.  If carbon is yet still available ...
    !        Our pools are now either on allometry or above (from fusion).
    !        We we can increment those pools at or below,
    !        including structure and reproduction according to their rates
    !        Use an adaptive euler integration. If the error is not nominal,
    !        the carbon balance sub-step (deltaC) will be halved and tried again
    ! -----------------------------------------------------------------------------------
    
    if( carbon_balance > nearzero ) then
       
       ! This routine checks that actual carbon is not below that targets. It does
       ! allow actual pools to be above the target, and in these cases, it sends
       ! a false on the "grow_<>" flag, allowing the plant to grow into these pools.
       ! It also checks to make sure that structural biomass is not above the target.
       if ( EDPftvarcon_inst%woody(ipft) == itrue ) then
          call TargetAllometryCheck(leaf_c, fnrt_c, sapw_c, &
                                    store_c, struct_c,       &
                                    target_leaf_c, target_fnrt_c, &
                                    target_sapw_c, target_store_c, target_struct_c, &
                                    grow_leaf, grow_fnrt, grow_sapw, grow_store)
       else
          grow_leaf  = .true.
          grow_fnrt = .true.
          grow_sapw  = .true.
          grow_store = .true.
       end if

       ! Initialize the adaptive integrator arrays and flags
       ! -----------------------------------------------------------------------------------
       ierr             = 1
       totalC           = carbon_balance
       nsteps           = 0
       
       c_pool(:) = 0.0_r8
       c_mask(:) = .false.

       c_pool(leaf_c_id)   = leaf_c
       c_pool(fnrt_c_id)   = fnrt_c
       c_pool(sapw_c_id)   = sapw_c
       c_pool(store_c_id)  = store_c
       c_pool(struct_c_id) = struct_c
       c_pool(repro_c_id)  = repro_c
       c_pool(dbh_id)      = dbh
       
       c_mask(leaf_c_id)   = grow_leaf
       c_mask(fnrt_c_id)   = grow_fnrt
       c_mask(sapw_c_id)   = grow_sapw
       c_mask(store_c_id)  = grow_store
       c_mask(struct_c_id) = .true.                ! Always increment dead on growth step
       c_mask(repro_c_id)  = .true.                ! Always calculate reproduction on growth
       c_mask(dbh_id)      = .true.                ! Always increment dbh on growth step
       
       if(ODESolve == 2) then
          this%ode_opt_step = totalC
       end if
       
       do while( ierr .ne. 0 )
          
          deltaC = min(totalC,this%ode_opt_step)
          if(ODESolve == 1) then
             call RKF45(AllomCGrowthDeriv,c_pool,c_mask,deltaC,totalC, &
                   max_trunc_error,intgr_params,c_pool_out,this%ode_opt_step,step_pass)
             
          elseif(ODESolve == 2) then
             call Euler(AllomCGrowthDeriv,c_pool,c_mask,deltaC,totalC,intgr_params,c_pool_out)
             !  step_pass = .true.
             call CheckIntegratedAllometries(c_pool_out(dbh_id),ipft,canopy_trim,  &
                   c_pool_out(leaf_c_id), c_pool_out(fnrt_c_id), c_pool_out(sapw_c_id), &
                   c_pool_out(store_c_id), c_pool_out(struct_c_id), &
                   c_mask(leaf_c_id), c_mask(fnrt_c_id), c_mask(sapw_c_id), &
                   c_mask(store_c_id),c_mask(struct_c_id),  max_trunc_error, step_pass)
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
             totalC    = totalC - deltaC
             c_pool(:) = c_pool_out(:)
          end if
          
          if(nsteps > max_substeps ) then
             write(fates_log(),*) 'Plant Growth Integrator could not find'
             write(fates_log(),*) 'a solution in less than ',max_substeps,' tries'
             write(fates_log(),*) 'Aborting'
             write(fates_log(),*) 'carbon_balance',carbon_balance
             write(fates_log(),*) 'deltaC',deltaC
             write(fates_log(),*) 'totalC',totalC
             write(fates_log(),*) 'leaf:',grow_leaf,target_leaf_c,target_leaf_c - leaf_c
             write(fates_log(),*) 'fnrt:',grow_fnrt,target_fnrt_c,target_fnrt_c - fnrt_c
             write(fates_log(),*) 'sap:',grow_sapw,target_sapw_c, target_sapw_c - sapw_c
             write(fates_log(),*) 'store:',grow_store,target_store_c,target_store_c - store_c
             write(fates_log(),*) 'dead:',target_struct_c,target_struct_c - struct_c
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          ! TotalC should eventually be whittled down to near zero
          ! At that point, update the actual states
          ! --------------------------------------------------------------------------------
          if( (totalC < calloc_abs_error) .and. (step_pass) )then

             ierr           = 0
             leaf_c_flux    = c_pool(leaf_c_id)   - leaf_c
             fnrt_c_flux    = c_pool(fnrt_c_id)   - fnrt_c
             sapw_c_flux    = c_pool(sapw_c_id)   - sapw_c
             store_c_flux   = c_pool(store_c_id)  - store_c
             struct_c_flux  = c_pool(struct_c_id) - struct_c
             repro_c_flux   = c_pool(repro_c_id)  - repro_c
             
             ! Make an adjustment to flux partitions to make it match remaining c balance
             flux_adj       = carbon_balance/(leaf_c_flux+fnrt_c_flux+sapw_c_flux + &
                                              store_c_flux+struct_c_flux+repro_c_flux)

             
             leaf_c_flux    = leaf_c_flux*flux_adj
             fnrt_c_flux    = fnrt_c_flux*flux_adj
             sapw_c_flux    = sapw_c_flux*flux_adj
             store_c_flux   = store_c_flux*flux_adj
             struct_c_flux  = struct_c_flux*flux_adj
             repro_c_flux   = repro_c_flux*flux_adj
             
             carbon_balance = carbon_balance - leaf_c_flux
             leaf_c         = leaf_c + leaf_c_flux
             
             carbon_balance = carbon_balance - fnrt_c_flux
             fnrt_c         = fnrt_c + fnrt_c_flux
             
             carbon_balance = carbon_balance - sapw_c_flux
             sapw_c         = sapw_c + sapw_c_flux
             
             carbon_balance = carbon_balance - store_c_flux
             store_c        = store_c + store_c_flux
             
             carbon_balance = carbon_balance - struct_c_flux
             struct_c       = struct_c + struct_c_flux
             
             carbon_balance = carbon_balance - repro_c_flux
             repro_c        = repro_c  + repro_c_flux
             
             dbh            = c_pool(dbh_id)

             ! THESE HAVE TO BE SET OUTSIDE OF THIS ROUTINE
             !!          cohort%seed_prod = cohort%seed_prod + brepro_flux / hlm_freq_day
             !!          cohort%dhdt      = (h_sub-cohort%hite)/hlm_freq_day
             !!          cohort%ddbhdt    = (dbh_sub-dbh_in)/hlm_freq_day
             
             if( abs(carbon_balance)>calloc_abs_error ) then
                write(fates_log(),*) 'carbon conservation error while integrating pools'
                write(fates_log(),*) 'along alometric curve'
                write(fates_log(),*) 'carbon_balance = ',carbon_balance,totalC
                write(fates_log(),*) 'exiting'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
          end if
       end do
    end if

  end associate
  
  return
  end subroutine DailyPRTAC
  
  ! =====================================================================================
  
  function AllomCGrowthDeriv(c_pools,c_mask,cbalance,intgr_params) result(dCdx)

      ! ---------------------------------------------------------------------------------
      ! This function calculates the derivatives for the carbon pools
      ! relative to the amount of carbon balance.  This function is based completely
      ! off of allometry, and assumes that there are no other species (ie nutrients) that
      ! govern allocation.
      ! ---------------------------------------------------------------------------------
      
      ! Arguments
      real(r8),intent(in), dimension(:) :: c_pools      ! Vector of carbon pools
                                                        ! dbh,leaf,root,sap,store,dead
      logical,intent(in), dimension(:)  :: c_mask       ! logical mask of active pools
                                                        ! some may be turned off
      real(r8),intent(in)               :: cbalance     ! The carbon balance of the
                                                        ! partial step (independant var)
                                             
      real(r8), intent(in),dimension(:) :: intgr_params  ! Generic Array used to pass 
                                                         ! parameters into this function


      ! Return Value
      real(r8),dimension(lbound(c_pools,dim=1):ubound(c_pools,dim=1)) :: dCdx 

      ! locals
      integer  :: ipft       ! PFT index
      real(r8) :: canopy_trim    ! Canopy trimming function (boundary condition [0-1]
      real(r8) :: ct_leaf    ! target leaf biomass, dummy var (kgC)
      real(r8) :: ct_fnrt   ! target fine-root biomass, dummy var (kgC)
      real(r8) :: ct_sap     ! target sapwood biomass, dummy var (kgC)
      real(r8) :: ct_agw     ! target aboveground wood, dummy var (kgC)
      real(r8) :: ct_bgw     ! target belowground wood, dummy var (kgC)
      real(r8) :: ct_store   ! target storage, dummy var (kgC)
      real(r8) :: ct_dead    ! target structural biomas, dummy var (kgC)
      real(r8) :: sapw_area             ! dummy sapwood area
      real(r8) :: ct_dleafdd     ! target leaf biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dfnrtdd    ! target fine-root biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dsapdd      ! target sapwood biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dagwdd      ! target AG wood biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dbgwdd      ! target BG wood biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dstoredd    ! target storage biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_ddeaddd     ! target structural biomass derivative wrt d, (kgC/cm)
      real(r8) :: ct_dtotaldd    ! target total (not reproductive) biomass derivative wrt d, (kgC/cm)
      real(r8) :: repro_fraction ! fraction of carbon balance directed towards reproduction (kgC/kgC)


      associate( dbh    => c_pools(dbh_id), &
                 cleaf  => c_pools(leaf_c_id), &
                 cfnrt  => c_pools(fnrt_c_id), &
                 csap   => c_pools(sapw_c_id), &
                 cstore => c_pools(store_c_id), &
                 cdead  => c_pools(struct_c_id), &
                 crepro => c_pools(repro_c_id), &    ! Unused (memoryless)
                 mask_dbh  => c_mask(dbh_id), &    ! Unused (dbh always grows)
                 mask_leaf => c_mask(leaf_c_id), &  
                 mask_fnrt => c_mask(fnrt_c_id), &
                 mask_sap  => c_mask(sapw_c_id), &
                 mask_store => c_mask(store_c_id), &
                 mask_dead  => c_mask(struct_c_id),  & ! Unused (dead always grows)
                 mask_repro => c_mask(repro_c_id) )

        canopy_trim = intgr_params(ac_bc_in_id_ctrim)
        ipft        = nint(intgr_params(ac_bc_in_id_pft))

        if(dbh>huge(dbh)) then
           print*,"BIG D IN DERIV:",dbh
           stop
        end if

        call bleaf(dbh,ipft,canopy_trim,ct_leaf,ct_dleafdd)
        call bfineroot(dbh,ipft,canopy_trim,ct_fnrt,ct_dfnrtdd)
        !call bsap_allom(dbh,ipft,canopy_trim,sapw_area,ct_sap,ct_dsapdd)
        call bsap_allom(dbh,ipft,canopy_trim,ct_sap,ct_dsapdd)
        call bagw_allom(dbh,ipft,ct_agw,ct_dagwdd)
        call bbgw_allom(dbh,ipft,ct_bgw,ct_dbgwdd)
        call bdead_allom(ct_agw,ct_bgw, ct_sap, ipft, ct_dead, &
                         ct_dagwdd, ct_dbgwdd, ct_dsapdd, ct_ddeaddd)
        call bstore_allom(dbh,ipft,canopy_trim,ct_store,ct_dstoredd)
        
        ! fraction of carbon going towards reproduction
        if (dbh <= EDPftvarcon_inst%dbh_repro_threshold(ipft)) then ! cap on leaf biomass
           repro_fraction = EDPftvarcon_inst%seed_alloc(ipft)
        else
           repro_fraction = EDPftvarcon_inst%seed_alloc(ipft) + EDPftvarcon_inst%seed_alloc_mature(ipft)
        end if

        dCdx = 0.0_r8

        ct_dtotaldd = ct_ddeaddd
        if (mask_leaf)  ct_dtotaldd = ct_dtotaldd + ct_dleafdd 
        if (mask_fnrt) ct_dtotaldd = ct_dtotaldd + ct_dfnrtdd 
        if (mask_sap)   ct_dtotaldd = ct_dtotaldd + ct_dsapdd
        if (mask_store) ct_dtotaldd = ct_dtotaldd + ct_dstoredd

        ! It is possible that with some asymptotic, or hard
        ! capped allometries, that all growth rates reach zero.
        ! In this case, if there is carbon, give it to reproduction

        if(ct_dtotaldd<=tiny(ct_dtotaldd))then

           dCdx(struct_c_id)  = 0.0_r8
           dCdx(dbh_id)    = 0.0_r8
           dCdx(leaf_c_id)  = 0.0_r8
           dCdx(fnrt_c_id) = 0.0_r8
           dCdx(sapw_c_id)   = 0.0_r8
           dCdx(store_c_id) = 0.0_r8
           dCdx(repro_c_id) = 1.0_r8

        else

           dCdx(struct_c_id) = (ct_ddeaddd/ct_dtotaldd)*(1.0_r8-repro_fraction)   
           dCdx(dbh_id)   = (1.0_r8/ct_dtotaldd)*(1.0_r8-repro_fraction)
        
           if (mask_leaf) then
              dCdx(leaf_c_id) = (ct_dleafdd/ct_dtotaldd)*(1.0_r8-repro_fraction)
           else
              dCdx(leaf_c_id) = 0.0_r8
           end if
           
           if (mask_fnrt) then
              dCdx(fnrt_c_id) = (ct_dfnrtdd/ct_dtotaldd)*(1.0_r8-repro_fraction)
           else
              dCdx(fnrt_c_id) = 0.0_r8
           end if
           
           if (mask_sap) then
              dCdx(sapw_c_id) = (ct_dsapdd/ct_dtotaldd)*(1.0_r8-repro_fraction)
           else
              dCdx(sapw_c_id) = 0.0_r8
           end if
           
           if (mask_store) then
              dCdx(store_c_id) = (ct_dstoredd/ct_dtotaldd)*(1.0_r8-repro_fraction)
           else
              dCdx(store_c_id) = 0.0_r8
           end if
           
           dCdx(repro_c_id) = repro_fraction

        end if
        
      end associate

      return
   end function AllomCGrowthDeriv

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
      class(callom_prt_vartypes) :: this     ! this class
      
      ! This routine does nothing, because in the carbon only allometric RT model
      ! we currently don't have any fast-timestep processes
      ! Think of this as a stub.
      
      



      return
    end subroutine FastPRTAC


end module PRTAllometricCarbonMod
  
