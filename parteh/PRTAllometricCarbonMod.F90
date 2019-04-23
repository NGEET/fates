module PRTAllometricCarbonMod

   ! ------------------------------------------------------------------------------------
   !
   ! This module contains all of the specific functions and types for
   ! Plant Allocation and Reactive Transport Extensible Hypotheses (PARTEH)
   ! CARBON only, allometric growth hypothesis
   ! 
   ! Adapted from code originally in ED, by Rosie Fisher and Paul Moorcroft
   ! This refactor written by : Ryan Knox Apr 2018
   !
   ! ------------------------------------------------------------------------------------

  use PRTGenericMod , only  : prt_global_type
  use PRTGenericMod , only  : prt_global
  use PRTGenericMod , only  : prt_vartype
  use PRTGenericMod , only  : prt_vartypes
  use PRTGenericMod , only  : carbon12_element
  use PRTGenericMod , only  : leaf_organ
  use PRTGenericMod , only  : fnrt_organ
  use PRTGenericMod , only  : sapw_organ
  use PRTGenericMod , only  : store_organ
  use PRTGenericMod , only  : repro_organ
  use PRTGenericMod , only  : struct_organ
  use PRTGenericMod , only  : un_initialized
  use PRTGenericMod , only  : prt_carbon_allom_hyp

  use FatesAllometryMod   , only : bleaf
  use FatesAllometryMod   , only : bsap_allom
  use FatesAllometryMod   , only : bfineroot
  use FatesAllometryMod   , only : bstore_allom
  use FatesAllometryMod   , only : bdead_allom
  use FatesAllometryMod   , only : bbgw_allom
  use FatesAllometryMod   , only : bagw_allom
  use FatesAllometryMod   , only : h_allom
  use FatesAllometryMod   , only : CheckIntegratedAllometries
  use FatesAllometryMod   , only : ForceDBH

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
  use FatesConstantsMod   , only : years_per_day


  implicit none
  private

  ! -------------------------------------------------------------------------------------
  !
  ! Define the state variables for this specific hypothesis.
  !
  ! -------------------------------------------------------------------------------------

  integer, parameter :: leaf_c_id   = 1   ! Unique object index for leaf carbon
  integer, parameter :: fnrt_c_id   = 2   ! Unique object index for fine-root carbon
  integer, parameter :: sapw_c_id   = 3   ! Unique object index for sapwood carbon 
  integer, parameter :: store_c_id  = 4   ! Unique object index for storage carbon
  integer, parameter :: repro_c_id  = 5   ! Unique object index for reproductive carbon
  integer, parameter :: struct_c_id = 6   ! Unique object index for structural carbon
  integer, parameter :: num_vars = 6      ! THIS MUST MATCH THE LARGEST INDEX ABOVE
  
  
  ! For this hypothesis, we integrate dbh along with the other 6. Since this
  ! is a boundary condition, we do not add it to the state array, but we do want
  ! to include it with the integrator array.

  integer, parameter :: dbh_id             = 7  ! This is just used for the integrator
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
  integer, parameter         :: num_bc_inout         = 2   ! Number of in & output boundary conditions


  integer, public, parameter :: ac_bc_in_id_pft   = 1   ! Index for the PFT input BC
  integer, public, parameter :: ac_bc_in_id_ctrim = 2   ! Index for the canopy trim function
  integer, parameter         :: num_bc_in         = 2   ! Number of input boundary condition

  ! THere are no purely output boundary conditions
  integer, parameter         :: num_bc_out        = 0   ! Number of purely output boundary condtions

  ! -------------------------------------------------------------------------------------
  ! Define the size of the coorindate vector.  For this hypothesis, there is only
  ! one pool per each species x organ combination.
  ! -------------------------------------------------------------------------------------
  integer, parameter         :: icd               = 1   ! Only 1 coordinate per variable

  
  ! This is the maximum number of leaf age pools  (used for allocating scratch space)
  integer, parameter         :: max_nleafage  = 10

  ! -------------------------------------------------------------------------------------
  ! This is the core type that holds this specific
  ! plant reactive transport (PRT) module
  ! -------------------------------------------------------------------------------------


  type, public, extends(prt_vartypes) :: callom_prt_vartypes

   contains

     procedure :: DailyPRT     => DailyPRTAllometricCarbon
     procedure :: FastPRT      => FastPRTAllometricCarbon

   end type callom_prt_vartypes
   
   ! ------------------------------------------------------------------------------------
   !
   ! This next class is an extention of the base instance that maps state variables
   !      to the outside model.
   !
   ! ------------------------------------------------------------------------------------
   
   character(len=*), parameter, private :: sourcefile = __FILE__


   ! This is the instance of the mapping table and variable definitions
   ! this is only allocated once per node.  This should be read-only
   ! everywhere in the code, except for where it is populated in this init routine
   ! below.

   class(prt_global_type), public, target, allocatable :: prt_global_ac


   public :: InitPRTGlobalAllometricCarbon


contains
  
 
  subroutine InitPRTGlobalAllometricCarbon()

     ! ----------------------------------------------------------------------------------
     ! Initialize and populate the object that holds the descriptions of the variables,
     ! and contains the mappings of each variable to the pre-ordained organ
     ! and species list, and the number of boundary conditions of each 3 types.
     !
     ! This is called very early on in the call sequence of the model, and should occur
     ! before any plants start being initialized.  These mapping tables must 
     ! exist before that happens.  This initialization only happens once on each
     ! machine, and the mapping will be read-only, and a global thing. This step
     ! is not initializing the data structures bound to the plants.
     !
     ! There are two mapping tables.  One mapping table is a 2d array organized
     ! by organ and species, that contains the variable index:
     ! 
     ! prt_global%sp_organ_map
     !
     ! The other mapping table is similar, but it is a 1D array, a list of the organs.
     ! And each of these the in turn points to a list of the indices associated
     ! with that organ.  This is useful when you want to do lots of stuff to a specified
     ! organ. 
     ! 
     ! prt_global%organ_map
     !
     ! IMPORTANT NOTE:  Once this object is populated, we can use this to properly
     ! allocate the "prt_vartypes_type" objects that attached to each plant. That process
     ! is handled by generic functions, and does not need to be written in each hypothesis.
     ! 
     ! -----------------------------------------------------------------------------------

     integer :: nleafage

     allocate(prt_global_ac)
     
     ! The "state descriptor" object holds things like the names, the symbols, the units
     ! of each variable. By putting it in an object, we can loop through them when
     ! doing things like reading/writing history and restarts

     allocate(prt_global_ac%state_descriptor(num_vars))

     prt_global_ac%hyp_name = 'Allometric Carbon Only'
     
     prt_global_ac%hyp_id = prt_carbon_allom_hyp

     ! Set mapping tables to zero
     call prt_global_ac%ZeroGlobal()

     
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

     ! Register the variables. Each variable must be associated with a global identifier
     ! for an organ and species.  Leaves are a little special in that they are discretized
     ! further by age class.  Although the code works fine if this collapses to 1.

     call prt_global_ac%RegisterVarInGlobal(leaf_c_id,"Leaf Carbon","leaf_c",leaf_organ,carbon12_element,nleafage)
     call prt_global_ac%RegisterVarInGlobal(fnrt_c_id,"Fine Root Carbon","fnrt_c",fnrt_organ,carbon12_element,icd)
     call prt_global_ac%RegisterVarInGlobal(sapw_c_id,"Sapwood Carbon","sapw_c",sapw_organ,carbon12_element,icd)
     call prt_global_ac%RegisterVarInGlobal(store_c_id,"Storage Carbon","store_c",store_organ,carbon12_element,icd)
     call prt_global_ac%RegisterVarInGlobal(struct_c_id,"Structural Carbon","struct_c",struct_organ,carbon12_element,icd)
     call prt_global_ac%RegisterVarInGlobal(repro_c_id,"Reproductive Carbon","repro_c",repro_organ,carbon12_element,icd)
     
     ! Set some of the array sizes for input and output boundary conditions
     prt_global_ac%num_bc_in    = num_bc_in
     prt_global_ac%num_bc_out   = num_bc_out
     prt_global_ac%num_bc_inout = num_bc_inout
     prt_global_ac%num_vars     = num_vars

     ! Have the global generic pointer, point to this hypothesis' object
     prt_global => prt_global_ac


     return
  end subroutine InitPRTGlobalAllometricCarbon

  ! =====================================================================================
  

  subroutine DailyPRTAllometricCarbon(this)

    ! -----------------------------------------------------------------------------------
    !
    ! This is the main routine that handles allocation associated with the 1st
    ! hypothesis;  carbon only, and growth governed by allometry
    ! 
    ! This routine is explained in the technical documentation in detail.
    !
    ! Some points:
    ! 1) dbh, while not a PARTEH "state variable", is passed in from FATES (or other
    !    model), is integrated along with the mass based state variables, and then
    !    passed back to the ecosystem model. It is a "inout" style boundary condition.
    !
    ! 2) It is assumed that both growth respiration, and maintenance respiration
    !    costs have already been paid, and therefore the "carbon_balance" boundary
    !    condition is the net carbon gained by the plant over the coarse of the day.
    !    Think of "daily integrated NPP".
    ! 
    ! 3) This routine will completely spend carbon_balance if it enters as a positive 
    !    value, or replace carbon balance (using storage) if it enters as a negative value.
    !    
    ! 4) It is assumed that the ecosystem model calling this routine has ensured that
    !    the net amount of negative carbon is no greater than that which can be replaced
    !    by storage.  This routine will crash gracefully if that is not true.
    !
    ! 5) Leaves and fine-roots are given top priority, but just to replace maintenance 
    !    turnover. This can also draw from strorage.
    ! 
    ! 6) Storage is given next available carbon gain, either to push up to zero, 
    !    or to use it to top off stores.
    !
    ! 7) Third priority is then given to leaves and fine-roots again, but can only use
    !    carbon gain. Also, this transfer will attempt to get pools up to allometry.
    ! 
    ! 8) Fourth priority is to bring other live pools up to allometry, and then structure.
    ! 
    ! 9) Finally, if carbon is yet still available, it will grow all pools out concurrently
    !    including some to reproduction.
    !
    ! ----------------------------------------------------------------------------------

    
    ! The class is the only argument
    class(callom_prt_vartypes)   :: this          ! this class

    ! -----------------------------------------------------------------------------------
    ! These are local copies of the in/out boundary condition structure
    ! -----------------------------------------------------------------------------------

    real(r8),pointer :: dbh            ! Diameter at breast height [cm]
                                       ! this local will point to both in and out bc's
    real(r8),pointer :: carbon_balance ! Daily carbon balance for this cohort [kgC]

    real(r8) :: canopy_trim            ! The canopy trimming function [0-1]
    integer  :: ipft                   ! Plant Functional Type index


    real(r8) :: target_leaf_c         ! target leaf carbon [kgC]
    real(r8) :: target_fnrt_c         ! target fine-root carbon [kgC]
    real(r8) :: target_sapw_c         ! target sapwood carbon [kgC]
    real(r8) :: target_store_c        ! target storage carbon [kgC]
    real(r8) :: target_agw_c          ! target above ground carbon in woody tissues [kgC]
    real(r8) :: target_bgw_c          ! target below ground carbon in woody tissues [kgC]
    real(r8) :: target_struct_c       ! target structural carbon [kgC]

    real(r8) :: sapw_area             ! dummy var, x-section area of sapwood [m2]

    real(r8) :: leaf_below_target     ! fineroot biomass below target amount [kgC]
    real(r8) :: fnrt_below_target     ! fineroot biomass below target amount [kgC]
    real(r8) :: sapw_below_target     ! sapwood biomass below target amount [kgC]
    real(r8) :: store_below_target    ! storage biomass below target amount [kgC]
    real(r8) :: struct_below_target   ! dead (structural) biomass below target amount [kgC]
    real(r8) :: total_below_target    ! total biomass below the allometric target [kgC]

    real(r8) :: flux_adj              ! adjustment made to growth flux term to minimize error [kgC]
    real(r8) :: store_target_fraction ! ratio between storage and leaf biomass when on allometry [kgC]

    real(r8) :: leaf_c_demand         ! leaf carbon that is demanded to replace maintenance turnover [kgC]
    real(r8) :: fnrt_c_demand         ! fineroot carbon that is demanded to replace 
                                      ! maintenance turnover [kgC]
    real(r8) :: total_c_demand        ! total carbon that is demanded to replace maintenance turnover [kgC]
    logical  :: step_pass             ! Did the integration step pass?

    real(r8) :: leaf_c_flux           ! Transfer into leaves at various stages [kgC]
    real(r8) :: fnrt_c_flux           ! Transfer into fine-roots at various stages [kgC]
    real(r8) :: sapw_c_flux           ! Transfer into sapwood at various stages [kgC]
    real(r8) :: store_c_flux          ! Transfer into storage at various stages [kgC]
    real(r8) :: repro_c_flux          ! Transfer into reproduction at the final stage [kgC]
    real(r8) :: struct_c_flux         ! Transfer into structure at various stages [kgC]

    real(r8),dimension(max_nleafage) :: leaf_c0 

                                      ! Initial value of carbon used to determine net flux
    real(r8) :: fnrt_c0               ! during this routine
    real(r8) :: sapw_c0               ! ""   
    real(r8) :: store_c0              ! ""
    real(r8) :: repro_c0              ! ""
    real(r8) :: struct_c0             ! ""

    logical  :: grow_struct
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

    integer  :: i_var                 ! index for iterating state variables
    integer  :: i_age                 ! index for iterating leaf ages
    integer  :: nleafage              ! number of leaf age classifications
    real(r8) :: leaf_age_flux         ! carbon mass flux between leaf age classification pools


    ! Integegrator variables c_pool is "mostly" carbon variables, it also includes
    ! dbh...
    ! -----------------------------------------------------------------------------------

    real(r8),dimension(n_integration_vars) :: c_pool     ! Vector of carbon pools passed to integrator
    real(r8),dimension(n_integration_vars) :: c_pool_out ! Vector of carbon pools passed back from integrator
    logical,dimension(n_integration_vars)  :: c_mask     ! Mask of active pools during integration

    integer , parameter :: max_substeps = 300            ! Maximum allowable iterations

    real(r8), parameter :: max_trunc_error = 1.0_r8      ! Maximum allowable truncation error

    integer,  parameter :: ODESolve = 2                  ! 1=RKF45,  2=Euler

    integer,  parameter :: iexp_leaf = 1                 ! index 1 is the expanding (i.e. youngest)
                                                         ! leaf age class, and therefore
                                                         ! all new allocation goes into that pool

    real(r8) ::  intgr_params(num_bc_in)                 ! The boundary conditions to this routine,
                                                         ! are pressed into an array that is also
                                                         ! passed to the integrators

    associate( & 
          leaf_c   => this%variables(leaf_c_id)%val, &
          fnrt_c   => this%variables(fnrt_c_id)%val(icd), &
          sapw_c   => this%variables(sapw_c_id)%val(icd), &
          store_c  => this%variables(store_c_id)%val(icd), &
          repro_c  => this%variables(repro_c_id)%val(icd), &
          struct_c => this%variables(struct_c_id)%val(icd))


    ! -----------------------------------------------------------------------------------
    ! 0.
    ! Copy the boundary conditions into readable local variables.
    ! We don't use pointers for bc's that ar "in" only, only "in-out" and "out"
    ! -----------------------------------------------------------------------------------

    dbh                             => this%bc_inout(ac_bc_inout_id_dbh)%rval
    carbon_balance                  => this%bc_inout(ac_bc_inout_id_netdc)%rval

    canopy_trim                     = this%bc_in(ac_bc_in_id_ctrim)%rval
    ipft                            = this%bc_in(ac_bc_in_id_pft)%ival

    intgr_params(:)                 = un_initialized
    intgr_params(ac_bc_in_id_ctrim) = this%bc_in(ac_bc_in_id_ctrim)%rval
    intgr_params(ac_bc_in_id_pft)   = real(this%bc_in(ac_bc_in_id_pft)%ival)
    
    ! -----------------------------------------------------------------------------------
    ! I. Remember the values for the state variables at the beginning of this
    ! routines. We will then use that to determine their net allocation and reactive
    ! transport flux "%net_alloc" at the end.
    ! -----------------------------------------------------------------------------------

    nleafage = prt_global%state_descriptor(leaf_c_id)%num_pos ! Number of leaf age class

    leaf_c0(1:nleafage) = leaf_c(1:nleafage)  ! Set initial leaf carbon 
    fnrt_c0 = fnrt_c                          ! Set initial fine-root carbon
    sapw_c0 = sapw_c                          ! Set initial sapwood carbon
    store_c0 = store_c                        ! Set initial storage carbon 
    repro_c0 = repro_c                        ! Set initial reproductive carbon
    struct_c0 = struct_c                      ! Set initial structural carbon

    
    ! -----------------------------------------------------------------------------------
    ! If we have more than one leaf age classification, allow
    ! some leaf biomass to transition to the older classes.  NOTE! This is not handling
    ! losses due to turnover (ie. flux from the oldest senescing class). This is only
    ! internal.
    ! (rgk 12-15-2018: Have Chonggang confirm that aging should not be restricted
    ! to evergreens)
    ! -----------------------------------------------------------------------------------

    if(nleafage>1) then
       do i_age = 1,nleafage-1
          if (EDPftvarcon_inst%leaf_long(ipft,i_age)>nearzero) then
             leaf_age_flux   = leaf_c0(i_age) * years_per_day / EDPftvarcon_inst%leaf_long(ipft,i_age)
             leaf_c(i_age)   = leaf_c(i_age) - leaf_age_flux
             leaf_c(i_age+1) = leaf_c(i_age+1) + leaf_age_flux
          end if
       end do
    end if
    

    ! -----------------------------------------------------------------------------------
    ! II. Calculate target size of the biomass compartment for a given dbh.   
    ! -----------------------------------------------------------------------------------
    
    ! Target sapwood biomass according to allometry and trimming [kgC]
    call bsap_allom(dbh,ipft,canopy_trim,sapw_area,target_sapw_c)
    
    ! Target total above ground biomass in woody/fibrous tissues  [kgC]
    call bagw_allom(dbh,ipft,target_agw_c)
    
    ! Target total below ground biomass in woody/fibrous tissues [kgC] 
    call bbgw_allom(dbh,ipft,target_bgw_c)
    
    ! Target total dead (structrual) biomass [kgC]
    call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)
    
    ! Target leaf biomass according to allometry and trimming
    call bleaf(dbh,ipft,canopy_trim,target_leaf_c)
    
    ! Target fine-root biomass and deriv. according to allometry and trimming [kgC, kgC/cm]
    call bfineroot(dbh,ipft,canopy_trim,target_fnrt_c)
    
    ! Target storage carbon [kgC,kgC/cm]
    call bstore_allom(dbh,ipft,canopy_trim,target_store_c)


    ! -----------------------------------------------------------------------------------
    ! III.  Prioritize some amount of carbon to replace leaf/root turnover
    !         Make sure it isnt a negative payment, and either pay what is available
    !         or forcefully pay from storage. 
    ! -----------------------------------------------------------------------------------
    
    if( EDPftvarcon_inst%evergreen(ipft) ==1 ) then
       leaf_c_demand   = max(0.0_r8, &
             EDPftvarcon_inst%leaf_stor_priority(ipft)*sum(this%variables(leaf_c_id)%turnover(:)))
    else
       leaf_c_demand   = 0.0_r8
    end if
    
    fnrt_c_demand   = max(0.0_r8, &
          EDPftvarcon_inst%leaf_stor_priority(ipft)*this%variables(fnrt_c_id)%turnover(icd))

    total_c_demand = leaf_c_demand + fnrt_c_demand
    
    if (total_c_demand> nearzero ) then

       ! We pay this even if we don't have the carbon
       ! Just don't pay so much carbon that storage+carbon_balance can't pay for it

       leaf_c_flux = min(leaf_c_demand, &
                         max(0.0_r8,(store_c+carbon_balance)* &
                         (leaf_c_demand/total_c_demand)))
       
       ! Add carbon to the youngest age pool (i.e iexp_leaf = index 1)
       carbon_balance    = carbon_balance   - leaf_c_flux
       leaf_c(iexp_leaf) = leaf_c(iexp_leaf) + leaf_c_flux

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
    
    leaf_c_demand   = max(0.0_r8,(target_leaf_c - sum(leaf_c(1:nleafage))))
    fnrt_c_demand   = max(0.0_r8,(target_fnrt_c - fnrt_c))

    total_c_demand = leaf_c_demand + fnrt_c_demand
    
    if( (carbon_balance > nearzero ) .and. (total_c_demand>nearzero)) then

       leaf_c_flux    = min(leaf_c_demand, &
                        carbon_balance*(leaf_c_demand/total_c_demand))
       carbon_balance = carbon_balance - leaf_c_flux
       leaf_c(iexp_leaf) = leaf_c(iexp_leaf) + leaf_c_flux
       
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

       leaf_below_target  = max(target_leaf_c - sum(leaf_c(1:nleafage)),0.0_r8)
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
          leaf_c(iexp_leaf)            = leaf_c(iexp_leaf) + leaf_c_flux
          
          carbon_balance               = carbon_balance - fnrt_c_flux
          fnrt_c                       = fnrt_c + fnrt_c_flux
          
          carbon_balance               = carbon_balance - sapw_c_flux
          sapw_c                       = sapw_c + sapw_c_flux
          
          carbon_balance               = carbon_balance - store_c_flux
          store_c                      = store_c  +  store_c_flux
          
       end if
    end if
    
    ! -----------------------------------------------------------------------------------
    ! VII.  If carbon is still available, replenish the structural pool to get
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
    ! VIII.  If carbon is yet still available ...
    !        Our pools are now either on allometry or above (from fusion).
    !        We we can increment those pools at or below,
    !        including structure and reproduction according to their rates
    !        Use an adaptive euler integration. If the error is not nominal,
    !        the carbon balance sub-step (deltaC) will be halved and tried again
    !
    ! Note that we compare against calloc_abs_error here because it is possible
    ! that all the carbon was effectively used up, but a miniscule amount
    ! remains due to numerical precision (ie -20 or so), so even though
    ! the plant has not been brought to be "on allometry", it thinks it has carbon
    ! left to allocate, and thus it must be on allometry when its not.
    ! -----------------------------------------------------------------------------------
    
    if( carbon_balance > calloc_abs_error ) then
       
       ! This routine checks that actual carbon is not below that targets. It does
       ! allow actual pools to be above the target, and in these cases, it sends
       ! a false on the "grow_<>" flag, allowing the plant to grow into these pools.
       ! It also checks to make sure that structural biomass is not above the target.

       if( (target_store_c - store_c)>calloc_abs_error) then
          write(fates_log(),*) 'storage is not on-allometry at the growth step'
          write(fates_log(),*) 'exiting'
          write(fates_log(),*) 'cbal: ',carbon_balance
          write(fates_log(),*) 'near-zero',nearzero
          write(fates_log(),*) 'store_c: ',store_c
          write(fates_log(),*) 'target c: ',target_store_c
          write(fates_log(),*) 'store_c0:', store_c0
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       

       call TargetAllometryCheck(sum(leaf_c(1:nleafage)), fnrt_c, sapw_c, &
                                 store_c, struct_c,       &
                                 target_leaf_c, target_fnrt_c, &
                                 target_sapw_c, target_store_c, target_struct_c, &
                                 grow_struct, grow_leaf, grow_fnrt, grow_sapw, grow_store)

       ! --------------------------------------------------------------------------------
       ! The numerical integration of growth requires that the instantaneous state
       ! variables are passed in as an array.  We call it "c_pool".
       !
       ! Initialize the adaptive integrator arrays and flags
       ! --------------------------------------------------------------------------------

       ierr             = 1
       totalC           = carbon_balance
       nsteps           = 0
       
       c_pool(:) = 0.0_r8                        ! Zero state variable array
       c_mask(:) = .false.                       ! This mask tells the integrator
                                                 ! which indices are active. Its possible
                                                 ! that due to fusion, or previous numerical
                                                 ! truncation errors, that one of these pools
                                                 ! may be larger than its target! We check
                                                 ! this, and if true, then we flag that
                                                 ! pool to be ignored. c_mask(i) = .false.
                                                 ! For grasses, since they don't grow very 
                                                 ! large and thus won't accumulate such large
                                                 ! errors, we always mask as true.

       c_pool(leaf_c_id)   = sum(leaf_c(1:nleafage))
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
       c_mask(struct_c_id) = grow_struct
       c_mask(repro_c_id)  = .true.                ! Always calculate reproduction on growth
       c_mask(dbh_id)      = .true.                ! Always increment dbh on growth step
       

       ! When using the Euler method, we keep things simple.  We always try
       ! to make the first integration step to span the entirety of the integration
       ! window for the independent variable (available carbon)

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
             
             ! When integrating along the allometric curve, we have the luxury of perfect
             ! hindsite.  Ie, after we have made our step, we can see if the amount
             ! of each carbon we have matches the target associated with the new dbh.
             ! The following call evaluates how close we are to the allometically defined
             ! targets. If we are too far (governed by max_trunc_error), then we
             ! pass back the pass/fail flag (step_pass) as false.  If false, then
             ! we halve the step-size, and then retry.  If that step was fine, then
             ! we remember the current step size as a good next guess.
             
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
             write(fates_log(),*) 'An integrator was chosen that does not exist'
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
             write(fates_log(),*) 'leaf:',grow_leaf,target_leaf_c,target_leaf_c - sum(leaf_c(:))
             write(fates_log(),*) 'fnrt:',grow_fnrt,target_fnrt_c,target_fnrt_c - fnrt_c
             write(fates_log(),*) 'sap:',grow_sapw,target_sapw_c, target_sapw_c - sapw_c
             write(fates_log(),*) 'store:',grow_store,target_store_c,target_store_c - store_c
             write(fates_log(),*) 'dead:',target_struct_c,target_struct_c - struct_c
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          !
          ! TotalC should eventually be whittled down to near zero.
          ! The solvers are not perfect, so we can't expect it to be perfectly zero.
          ! Note that calloc_abs_error is 1e-9, which is really small (1 microgram of carbon)
          ! yet also six orders of magnitude greater than typical rounding errors (~1e-15).
  
          ! At that point, update the actual states
          ! --------------------------------------------------------------------------------
          if( (totalC < calloc_abs_error) .and. (step_pass) )then

             ierr           = 0
             leaf_c_flux    = c_pool(leaf_c_id)   - sum(leaf_c(1:nleafage))
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
             
             carbon_balance    = carbon_balance - leaf_c_flux
             leaf_c(iexp_leaf) = leaf_c(iexp_leaf) + leaf_c_flux
             
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

    ! Track the net allocations and transport from this routine

    do i_age = 1,nleafage
       this%variables(leaf_c_id)%net_alloc(i_age) = &
             this%variables(leaf_c_id)%net_alloc(i_age) + &
             (leaf_c(i_age) - leaf_c0(i_age))
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



   end associate
  
   return
  end subroutine DailyPRTAllometricCarbon
  
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
      ! Change in carbon (each pool) per change in total allocatable carbon (kgC/kgC)
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
      real(r8) :: sapw_area      ! dummy sapwood area
      real(r8) :: ct_dleafdd     ! target leaf biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dfnrtdd     ! target fine-root biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dsapdd      ! target sapwood biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dagwdd      ! target AG wood biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dbgwdd      ! target BG wood biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dstoredd    ! target storage biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_ddeaddd     ! target structural biomass derivative wrt diameter, (kgC/cm)
      real(r8) :: ct_dtotaldd    ! target total (not reproductive) biomass derivative wrt diameter, (kgC/cm)
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
        ipft        = int(intgr_params(ac_bc_in_id_pft))

        call bleaf(dbh,ipft,canopy_trim,ct_leaf,ct_dleafdd)
        call bfineroot(dbh,ipft,canopy_trim,ct_fnrt,ct_dfnrtdd)
        call bsap_allom(dbh,ipft,canopy_trim,sapw_area,ct_sap,ct_dsapdd)

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
                                   grow_dead,grow_leaf,grow_froot,grow_sapw,grow_store)

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
     logical,intent(out) :: grow_dead
       
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
     elseif( (bdead-bt_dead)> calloc_abs_error) then
        grow_dead = .false.
     else
        grow_dead = .true.
     end if
     

     return
   end subroutine TargetAllometryCheck

   ! =====================================================================================

   subroutine FastPRTAllometricCarbon(this)
      
      implicit none
      class(callom_prt_vartypes) :: this     ! this class
      
      ! This routine does nothing, because in the carbon only allometric RT model
      ! we currently don't have any fast-timestep processes
      ! Think of this as a stub.


      return
   end subroutine FastPRTAllometricCarbon


end module PRTAllometricCarbonMod
  
