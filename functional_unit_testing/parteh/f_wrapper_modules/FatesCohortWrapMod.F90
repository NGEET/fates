! =======================================================================================
!
! This is the wrapper module that provides FATES data structures
!
! =======================================================================================

module FatesCohortWrapMod
   
  use iso_c_binding, only : r8 => c_double
  use iso_c_binding, only : i4 => c_int
  use iso_c_binding, only : c_char
  use FatesAllometryMod,     only : bleaf
  use FatesAllometryMod,     only : bfineroot
  use FatesAllometryMod,     only : bsap_allom
  use FatesAllometryMod,     only : bagw_allom
  use FatesAllometryMod,     only : bbgw_allom
  use FatesAllometryMod,     only : bdead_allom
  use FatesAllometryMod,     only : bstore_allom
  use FatesAllometryMod,     only : h2d_allom
  use FatesAllometryMod,     only : tree_lai
  use FatesAllometryMod,     only : carea_allom

  use EDPftvarcon,            only : EDPftvarcon_inst

  use PRTGenericMod,          only : prt_vartypes
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : nitrogen_element
  use PRTGenericMod,          only : phosphorous_element
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : SetState
  use PRTGenericMod,          only : prt_global

  use PRTAllometricCarbonMod, only : callom_prt_vartypes
  use PRTAllometricCarbonMod, only : ac_bc_inout_id_netdc
  use PRTAllometricCarbonMod, only : ac_bc_in_id_pft
  use PRTAllometricCarbonMod, only : ac_bc_in_id_ctrim
  use PRTAllometricCarbonMod, only : ac_bc_inout_id_dbh
  use PRTAllometricCarbonMod, only : prt_global_ac

  use PRTLossFluxesMod, only : PRTMaintTurnover
  use PRTLossFluxesMod, only : PRTDeciduousTurnover
  use PRTLossFluxesMod, only : PRTPhenologyFlush

  use PRTAllometricCNPMod,    only : prt_global_acnp
  use PRTAllometricCNPMod,    only : cnp_allom_prt_vartypes
  use PRTAllometricCNPMod,    only : acnp_bc_inout_id_dbh
  use PRTAllometricCNPMod,    only : acnp_bc_inout_id_netdc
  use PRTAllometricCNPMod,    only : acnp_bc_inout_id_rmaint_def
  use PRTAllometricCNPMod,    only : acnp_bc_inout_id_netdn
  use PRTAllometricCNPMod,    only : acnp_bc_inout_id_netdp
  use PRTAllometricCNPMod,    only : acnp_bc_in_id_ctrim
  use PRTAllometricCNPMod,    only : acnp_bc_in_id_pft
  use PRTAllometricCNPMod,    only : acnp_bc_in_id_status
  use PRTAllometricCNPMod,    only : acnp_bc_out_id_rootcexude
  use PRTAllometricCNPMod,    only : acnp_bc_out_id_rootnexude
  use PRTAllometricCNPMod,    only : acnp_bc_out_id_rootpexude
  use PRTAllometricCNPMod,    only : acnp_bc_out_id_growresp


  use FatesConstantsMod   , only : nearzero

  use EDTypesMod            , only : nclmax

  use FatesGlobals          , only : endrun => fates_endrun
  use FatesGlobals          , only : fates_log
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  
  implicit none
  private ! Modules are private by default
  
  type, public :: ed_cohort_type
     
     integer  :: pft             ! pft number
     integer  :: parteh_model        ! The PARTEH allocation hypothesis used
     real(r8) :: dbh             ! dbh: cm
     integer  :: status_coh      ! leaf status 1=off, 2=on
     real(r8) :: canopy_trim     ! Trimming function for the canopy
     
     real(r8) :: dhdt            ! time derivative of height       : m/year
     real(r8) :: ddbhdt          ! time derivative of dbh          : cm/year

     real(r8) :: daily_carbon_gain        ! 
     real(r8) :: daily_nitrogen_gain      !
     real(r8) :: daily_phosphorous_gain   !
     real(r8) :: daily_r_grow             !
     real(r8) :: daily_r_maint            !
     real(r8) :: daily_r_maint_demand          !
     real(r8) :: accum_r_maint_deficit  !
     real(r8) :: carbon_root_exudate           !
     real(r8) :: nitrogen_root_exudate         !
     real(r8) :: phosphorous_root_exudate      !

     
     ! Multi-species, multi-pool Reactive Transport 
     class(prt_vartypes), pointer :: prt
     
  end type ed_cohort_type
  
  ! Global Instances
  
  type(ed_cohort_type), pointer, public       :: cohort_array(:)
  integer, public :: numcohort

  character(len=*), parameter, private :: sourcefile = __FILE__

  ! Make necessary procedures public
  public :: CohortInitAlloc
  public :: CohortPySet
  public :: WrapDailyPRT
  public :: WrapQueryVars
  public :: WrapQueryDiagnostics
  
contains
  
  subroutine CohortInitAlloc(numcohorts)
    
    ! Arguments
    integer(i4), intent(in) :: numcohorts
    
    ! Locals
    integer(i4)                   :: ico
    type(ed_cohort_type), pointer :: ccohort
    
    
    allocate(cohort_array(numcohorts))
    
    do ico = 1,numcohorts
       ccohort               => cohort_array(ico)
       ccohort%parteh_model              = -1
       ccohort%pft                      = -9
       ccohort%dbh                      = -999.9_r8
       ccohort%status_coh               = -1
       ccohort%canopy_trim              = -999.9_r8
       ccohort%dhdt                     = -999.9_r8
       ccohort%ddbhdt                   = -999.9_r8
       ccohort%daily_carbon_gain        = -999.9_r8
       ccohort%daily_nitrogen_gain      = -999.9_r8
       ccohort%daily_phosphorous_gain   = -999.9_r8
       ccohort%daily_r_grow             = -999.9_r8
       ccohort%daily_r_maint            = -999.9_r8
       ccohort%daily_r_maint_demand     = -999.9_r8
       ccohort%accum_r_maint_deficit    = -999.9_r8
       ccohort%carbon_root_exudate      = -999.9_r8
       ccohort%nitrogen_root_exudate    = -999.9_r8
       ccohort%phosphorous_root_exudate = -999.9_r8
    end do

    return
  end subroutine CohortInitAlloc
  
  ! =====================================================================================  
  
  subroutine CohortPySet(ipft,hgt_min,canopy_trim)
    
    implicit none
    ! Arguments
    integer(i4),intent(in) :: ipft
    real(r8),intent(in)    :: hgt_min
    real(r8),intent(in)    :: canopy_trim
 
    ! Locals

    type(ed_cohort_type), pointer :: ccohort   ! Current cohort
    real(r8) :: leaf_c
    real(r8) :: fnrt_c
    real(r8) :: sapw_c
    real(r8) :: agw_c
    real(r8) :: bgw_c
    real(r8) :: struct_c
    real(r8) :: repro_c
    real(r8) :: store_c

    real(r8) :: sapw_area ! dummy area cross-sec

    real(r8) :: leaf_n
    real(r8) :: fnrt_n
    real(r8) :: sapw_n
    real(r8) :: struct_n
    real(r8) :: repro_n
    real(r8) :: store_n
    real(r8) :: leaf_p
    real(r8) :: fnrt_p
    real(r8) :: sapw_p
    real(r8) :: struct_p
    real(r8) :: repro_p
    real(r8) :: store_p


    class(callom_prt_vartypes), pointer :: callom_prt
    class(cnp_allom_prt_vartypes), pointer :: cnpallom_prt

    
    ccohort => cohort_array(ipft)
    

    ccohort%pft          = int(ipft)
    ccohort%parteh_model = int(EDPftvarcon_inst%parteh_model(ipft))

    call h2d_allom(hgt_min,ipft,ccohort%dbh)
    ccohort%canopy_trim = canopy_trim

    ! Use allometry to compute initial values
    
    ! Leaf biomass (carbon)
    call bleaf(ccohort%dbh, ipft, canopy_trim, leaf_c)
    
    ! Fine-root biomass (carbon)
    call bfineroot(ccohort%dbh, ipft, canopy_trim, fnrt_c)
    
    ! Sapwood biomass (carbon)
    call bsap_allom(ccohort%dbh, ipft, canopy_trim, sapw_area, sapw_c)
    
    ! Above ground woody biomass (carbon)
    call bagw_allom(ccohort%dbh, ipft, agw_c)
    
    ! Below ground woody biomass (carbon)
    call bbgw_allom(ccohort%dbh, ipft, bgw_c)
    
    ! Total structural biomass (carbon)
    call bdead_allom(agw_c, bgw_c, sapw_c, ipft, struct_c) 
    
    ! Target storage carbon [kgC,kgC/cm]
    call bstore_allom(ccohort%dbh, ipft, canopy_trim, store_c)
    
    repro_c = 0.0_r8

    
    select case(ccohort%parteh_model)
    case (1)
       prt_global => prt_global_ac
       allocate(callom_prt)
       ccohort%prt => callom_prt

    case(2)
       prt_global => prt_global_acnp
       allocate(cnpallom_prt)
       ccohort%prt => cnpallom_prt
       
    case DEFAULT
       write(fates_log(),*) 'You specified an unknown PRT module'
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    call ccohort%prt%InitPRTVartype()

    select case(ccohort%parteh_model)
    case (1)

       call SetState(ccohort%prt,leaf_organ, carbon12_element, leaf_c)
       call SetState(ccohort%prt,fnrt_organ, carbon12_element, fnrt_c)
       call SetState(ccohort%prt,sapw_organ, carbon12_element, sapw_c)
       call SetState(ccohort%prt,store_organ, carbon12_element, store_c)
       call SetState(ccohort%prt,struct_organ , carbon12_element, struct_c)
       call SetState(ccohort%prt,repro_organ , carbon12_element, repro_c)

       call ccohort%prt%RegisterBCInOut(ac_bc_inout_id_dbh,bc_rval = ccohort%dbh)
       call ccohort%prt%RegisterBCInOut(ac_bc_inout_id_netdc,bc_rval = ccohort%daily_carbon_gain)

       call ccohort%prt%RegisterBCIn(ac_bc_in_id_pft,bc_ival = ccohort%pft)
       call ccohort%prt%RegisterBCIn(ac_bc_in_id_ctrim,bc_rval = ccohort%canopy_trim)

    case (2)

       ! Initializing with the target stoichiometric ratios
       ! (OR you can initialize with the minimum ratios too.... p2)
       leaf_n = leaf_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,leaf_organ)
       fnrt_n = fnrt_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,fnrt_organ)
       sapw_n = sapw_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,sapw_organ)
       store_n = store_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,store_organ)
       struct_n = struct_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,struct_organ)
       repro_n = repro_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,repro_organ)
       
       leaf_p = leaf_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,leaf_organ)
       fnrt_p = fnrt_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,fnrt_organ)
       sapw_p = sapw_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,sapw_organ)
       store_p = store_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,store_organ)
       struct_p = struct_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,struct_organ)
       repro_p = repro_c * EDPftvarcon_inst%prt_phos_stoich_p1(ipft,repro_organ)

       ccohort%accum_r_maint_deficit = 0.0_r8
       
       call SetState(ccohort%prt,leaf_organ, carbon12_element, leaf_c)
       call SetState(ccohort%prt,fnrt_organ, carbon12_element, fnrt_c)
       call SetState(ccohort%prt,sapw_organ, carbon12_element, sapw_c)
       call SetState(ccohort%prt,store_organ, carbon12_element, store_c)
       call SetState(ccohort%prt,struct_organ , carbon12_element, struct_c)
       call SetState(ccohort%prt,repro_organ , carbon12_element, repro_c)

       call SetState(ccohort%prt,leaf_organ, nitrogen_element, leaf_n)
       call SetState(ccohort%prt,fnrt_organ, nitrogen_element, fnrt_n)
       call SetState(ccohort%prt,sapw_organ, nitrogen_element, sapw_n)
       call SetState(ccohort%prt,store_organ, nitrogen_element, store_n)
       call SetState(ccohort%prt,struct_organ , nitrogen_element, struct_n)
       call SetState(ccohort%prt,repro_organ , nitrogen_element, repro_n)

       call SetState(ccohort%prt,leaf_organ, phosphorous_element, leaf_p)
       call SetState(ccohort%prt,fnrt_organ, phosphorous_element, fnrt_p)
       call SetState(ccohort%prt,sapw_organ, phosphorous_element, sapw_p)
       call SetState(ccohort%prt,store_organ, phosphorous_element, store_p)
       call SetState(ccohort%prt,struct_organ , phosphorous_element, struct_p)
       call SetState(ccohort%prt,repro_organ , phosphorous_element, repro_p)
       
       ! Register In/Out Boundary Conditions
       call ccohort%prt%RegisterBCInOut(acnp_bc_inout_id_dbh,bc_rval   = ccohort%dbh)
       call ccohort%prt%RegisterBCInOut(acnp_bc_inout_id_netdc,bc_rval = ccohort%daily_carbon_gain)
       call ccohort%prt%RegisterBCInOut(acnp_bc_inout_id_netdn,bc_rval = ccohort%daily_nitrogen_gain)
       call ccohort%prt%RegisterBCInOut(acnp_bc_inout_id_netdp,bc_rval = ccohort%daily_phosphorous_gain)
       call ccohort%prt%RegisterBCInOut(acnp_bc_inout_id_rmaint_def, bc_rval = ccohort%accum_r_maint_deficit) 

       ! Register Input only BC's
       call ccohort%prt%RegisterBCIn(acnp_bc_in_id_pft,bc_ival   = ccohort%pft)
       call ccohort%prt%RegisterBCIn(acnp_bc_in_id_ctrim,bc_rval = ccohort%canopy_trim)
       call ccohort%prt%RegisterBCIn(acnp_bc_in_id_status,bc_ival = ccohort%status_coh)

       ! Register Output Boundary Conditions
       call ccohort%prt%RegisterBCOut(acnp_bc_out_id_rootcexude,bc_rval = ccohort%carbon_root_exudate)
       call ccohort%prt%RegisterBCOut(acnp_bc_out_id_rootnexude,bc_rval = ccohort%nitrogen_root_exudate)
       call ccohort%prt%RegisterBCOut(acnp_bc_out_id_rootpexude,bc_rval = ccohort%phosphorous_root_exudate)
       call ccohort%prt%RegisterBCOut(acnp_bc_out_id_growresp,bc_rval = ccohort%daily_r_grow )

    
    end select
    
    call ccohort%prt%CheckInitialConditions()

    
  end subroutine CohortPySet

  ! =====================================================================================

  subroutine WrapDailyPRT(ipft,daily_carbon_gain,canopy_trim,flush_c,drop_frac_c,leaf_status, &
       daily_nitrogen_gain, daily_phosphorous_gain,daily_r_maint_demand )
    
    implicit none
    ! Arguments
    integer(i4),intent(in) :: ipft
    real(r8),intent(in)    :: daily_carbon_gain
    real(r8),intent(in)    :: canopy_trim
    real(r8),intent(in)    :: flush_c
    real(r8),intent(in)    :: drop_frac_c
    integer,intent(in)     :: leaf_status
    real(r8), intent(in), optional :: daily_nitrogen_gain
    real(r8), intent(in), optional :: daily_phosphorous_gain
    real(r8), intent(in), optional :: daily_r_maint_demand

    type(ed_cohort_type), pointer :: ccohort
    

    ccohort               => cohort_array(ipft)

    ccohort%status_coh = leaf_status

    ! Zero the rate of change and the turnover arrays

    call ccohort%prt%ZeroRates()

    call PRTDeciduousTurnover(ccohort%prt, ipft, leaf_organ , drop_frac_c)
    
    call PRTPhenologyFlush(ccohort%prt, ipft, leaf_organ, flush_c)

    call PRTMaintTurnover(ccohort%prt, ipft)

    select case(int(ccohort%parteh_model))
    case (1)
       prt_global => prt_global_ac
       ccohort%daily_carbon_gain = daily_carbon_gain

       call ccohort%prt%DailyPRT()

       ccohort%daily_r_grow = 0.0_r8
       ccohort%carbon_root_exudate = 0.0_r8

    case (2)
       prt_global => prt_global_acnp
       ccohort%daily_carbon_gain      = daily_carbon_gain
       ccohort%daily_nitrogen_gain    = daily_nitrogen_gain
       ccohort%daily_phosphorous_gain = daily_phosphorous_gain
       ccohort%accum_r_maint_deficit  = ccohort%accum_r_maint_deficit + &
                                        daily_r_maint_demand 

       call ccohort%prt%DailyPRT()


    case DEFAULT
       write(fates_log(),*) 'You specified an unknown PRT module'
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
 

    call ccohort%prt%CheckMassConservation(ipft,1)

    
    
    return
  end subroutine WrapDailyPRT
  
  ! =====================================================================================
  
  subroutine WrapQueryVars(ipft,leaf_area,crown_area,agb,store_c,target_leaf_c)
    
    implicit none
    ! Arguments
    integer(i4),intent(in) :: ipft
    real(r8),intent(out)    :: leaf_area
    real(r8),intent(out)    :: crown_area
    real(r8),intent(out)    :: agb
    real(r8),intent(out)    :: store_c
    real(r8),intent(out)    :: target_leaf_c

    real(r8) :: leaf_c
    type(ed_cohort_type), pointer :: ccohort

    real(r8),parameter :: nplant = 1.0_r8
    real(r8),parameter :: site_spread = 1.0_r8

    real(r8), parameter, dimension(nclmax) :: canopy_lai = [0.0_r8,0.0_r8,0.0_r8,0.0_r8]
    integer, parameter  :: cl1 = 1
    
    ccohort     => cohort_array(ipft)


    select case(int(ccohort%parteh_model))
    case (1)
       prt_global => prt_global_ac
    case (2)
       prt_global => prt_global_acnp
    end select
    

    leaf_c  = ccohort%prt%GetState(leaf_organ, all_carbon_elements )
    store_c = ccohort%prt%GetState(store_organ, all_carbon_elements )
    
    call carea_allom(ccohort%dbh,nplant,site_spread,ipft,crown_area)

    leaf_area = crown_area*tree_lai(leaf_c, ipft, crown_area, nplant, cl1, canopy_lai) 

    call bagw_allom(ccohort%dbh,ipft,agb)

    call bleaf(ccohort%dbh,ipft, ccohort%canopy_trim, target_leaf_c)


    return
 end subroutine WrapQueryVars
  
 
 subroutine WrapQueryDiagnostics(ipft, dbh, &
                                 leaf_c, fnrt_c, sapw_c, store_c, struct_c, repro_c, &
                                 leaf_cturn, fnrt_cturn, sapw_cturn, store_cturn, struct_cturn, &
                                 leaf_n, fnrt_n, sapw_n, store_n, struct_n, repro_n, &
                                 leaf_nturn, fnrt_nturn, sapw_nturn, store_nturn, struct_nturn, &
                                 leaf_p, fnrt_p, sapw_p, store_p, struct_p, repro_p, &
                                 leaf_pturn, fnrt_pturn, sapw_pturn, store_pturn, struct_pturn, &
                                 crown_area, &
                                 carbon_root_exudate, nitrogen_root_exudate, phosphorous_root_exudate, &
                                 growth_resp )
    
    implicit none
    ! Arguments
    integer(i4),intent(in) :: ipft
    real(r8),intent(out)   :: dbh

    real(r8),intent(out)   :: leaf_c
    real(r8),intent(out)   :: fnrt_c
    real(r8),intent(out)   :: sapw_c
    real(r8),intent(out)   :: store_c
    real(r8),intent(out)   :: struct_c
    real(r8),intent(out)   :: repro_c
    real(r8),intent(out)   :: leaf_cturn
    real(r8),intent(out)   :: fnrt_cturn
    real(r8),intent(out)   :: sapw_cturn
    real(r8),intent(out)   :: store_cturn
    real(r8),intent(out)   :: struct_cturn

    real(r8),intent(out)   :: leaf_n
    real(r8),intent(out)   :: fnrt_n
    real(r8),intent(out)   :: sapw_n
    real(r8),intent(out)   :: store_n
    real(r8),intent(out)   :: struct_n
    real(r8),intent(out)   :: repro_n
    real(r8),intent(out)   :: leaf_nturn
    real(r8),intent(out)   :: fnrt_nturn
    real(r8),intent(out)   :: sapw_nturn
    real(r8),intent(out)   :: store_nturn
    real(r8),intent(out)   :: struct_nturn

    real(r8),intent(out)   :: leaf_p
    real(r8),intent(out)   :: fnrt_p
    real(r8),intent(out)   :: sapw_p
    real(r8),intent(out)   :: store_p
    real(r8),intent(out)   :: struct_p
    real(r8),intent(out)   :: repro_p
    real(r8),intent(out)   :: leaf_pturn
    real(r8),intent(out)   :: fnrt_pturn
    real(r8),intent(out)   :: sapw_pturn
    real(r8),intent(out)   :: store_pturn
    real(r8),intent(out)   :: struct_pturn


    real(r8),intent(out)   :: carbon_root_exudate
    real(r8),intent(out)   :: nitrogen_root_exudate
    real(r8),intent(out)   :: phosphorous_root_exudate
    real(r8),intent(out)   :: growth_resp
    
    real(r8),intent(out)   :: crown_area
    type(ed_cohort_type), pointer :: ccohort
    real(r8),parameter :: nplant = 1.0_r8
    real(r8),parameter :: site_spread = 1.0_r8

    ccohort => cohort_array(ipft)   
    
    select case(int(ccohort%parteh_model))
    case (1)
       prt_global => prt_global_ac
    case (2)
       prt_global => prt_global_acnp
    end select
    
    dbh    = ccohort%dbh

    leaf_c = ccohort%prt%GetState(organ_id=leaf_organ, species_id=all_carbon_elements)
    fnrt_c = ccohort%prt%GetState(organ_id=fnrt_organ, species_id=all_carbon_elements)
    sapw_c = ccohort%prt%GetState(organ_id=sapw_organ, species_id=all_carbon_elements)
    store_c = ccohort%prt%GetState(organ_id=store_organ, species_id=all_carbon_elements)
    struct_c = ccohort%prt%GetState(organ_id=struct_organ, species_id=all_carbon_elements)
    repro_c = ccohort%prt%GetState(organ_id=repro_organ, species_id=all_carbon_elements)

    leaf_cturn = ccohort%prt%GetTurnover(organ_id=leaf_organ, species_id=all_carbon_elements)
    fnrt_cturn = ccohort%prt%GetTurnover(organ_id=fnrt_organ, species_id=all_carbon_elements)
    sapw_cturn = ccohort%prt%GetTurnover(organ_id=sapw_organ, species_id=all_carbon_elements)
    store_cturn = ccohort%prt%GetTurnover(organ_id=store_organ, species_id=all_carbon_elements)
    struct_cturn = ccohort%prt%GetTurnover(organ_id=struct_organ, species_id=all_carbon_elements)

    leaf_n = ccohort%prt%GetState(organ_id=leaf_organ, species_id=nitrogen_element)
    fnrt_n = ccohort%prt%GetState(organ_id=fnrt_organ, species_id=nitrogen_element)
    sapw_n = ccohort%prt%GetState(organ_id=sapw_organ, species_id=nitrogen_element)
    store_n = ccohort%prt%GetState(organ_id=store_organ, species_id=nitrogen_element)
    struct_n = ccohort%prt%GetState(organ_id=struct_organ, species_id=nitrogen_element)
    repro_n = ccohort%prt%GetState(organ_id=repro_organ, species_id=nitrogen_element)
    
    leaf_nturn = ccohort%prt%GetTurnover(organ_id=leaf_organ, species_id=nitrogen_element)
    fnrt_nturn = ccohort%prt%GetTurnover(organ_id=fnrt_organ, species_id=nitrogen_element)
    sapw_nturn = ccohort%prt%GetTurnover(organ_id=sapw_organ, species_id=nitrogen_element)
    store_nturn = ccohort%prt%GetTurnover(organ_id=store_organ, species_id=nitrogen_element)
    struct_nturn = ccohort%prt%GetTurnover(organ_id=struct_organ, species_id=nitrogen_element)

    leaf_p = ccohort%prt%GetState(organ_id=leaf_organ, species_id=phosphorous_element)
    fnrt_p = ccohort%prt%GetState(organ_id=fnrt_organ, species_id=phosphorous_element)
    sapw_p = ccohort%prt%GetState(organ_id=sapw_organ, species_id=phosphorous_element)
    store_p = ccohort%prt%GetState(organ_id=store_organ, species_id=phosphorous_element)
    struct_p = ccohort%prt%GetState(organ_id=struct_organ, species_id=phosphorous_element)
    repro_p = ccohort%prt%GetState(organ_id=repro_organ, species_id=phosphorous_element)
    
    leaf_pturn = ccohort%prt%GetTurnover(organ_id=leaf_organ, species_id=phosphorous_element)
    fnrt_pturn = ccohort%prt%GetTurnover(organ_id=fnrt_organ, species_id=phosphorous_element)
    sapw_pturn = ccohort%prt%GetTurnover(organ_id=sapw_organ, species_id=phosphorous_element)
    store_pturn = ccohort%prt%GetTurnover(organ_id=store_organ, species_id=phosphorous_element)
    struct_pturn = ccohort%prt%GetTurnover(organ_id=struct_organ, species_id=phosphorous_element)

    growth_resp = ccohort%daily_r_grow

    call carea_allom(ccohort%dbh,nplant,site_spread,ipft,crown_area)

    carbon_root_exudate = ccohort%carbon_root_exudate
    nitrogen_root_exudate = ccohort%nitrogen_root_exudate
    phosphorous_root_exudate = ccohort%phosphorous_root_exudate

    return
 end subroutine WrapQueryDiagnostics



   
end module FatesCohortWrapMod
