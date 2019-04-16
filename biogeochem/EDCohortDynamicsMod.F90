module EDCohortDynamicsMod
  !
  ! !DESCRIPTION:
  ! Cohort stuctures in ED. 
  !
  ! !USES: 
  use FatesGlobals          , only : endrun => fates_endrun
  use FatesGlobals          , only : fates_log
  use FatesInterfaceMod     , only : hlm_freq_day
  use FatesInterfaceMod     , only : bc_in_type
  use FatesInterfaceMod     , only : hlm_use_planthydro
  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : fates_unset_int
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : fates_unset_r8
  use FatesConstantsMod     , only : nearzero
  use FatesConstantsMod     , only : calloc_abs_error
  use FatesInterfaceMod     , only : hlm_days_per_year
  use FatesInterfaceMod     , only : nleafage
  use EDPftvarcon           , only : EDPftvarcon_inst
  use FatesParameterDerivedMod, only : param_derived
  use EDTypesMod            , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod            , only : nclmax
  use EDTypesMod            , only : ncwd
  use EDTypesMod            , only : maxCohortsPerPatch
  use EDTypesMod            , only : AREA
  use EDTypesMod            , only : min_npm2, min_nppatch
  use EDTypesMod            , only : min_n_safemath
  use EDTypesMod            , only : nlevleaf
  use EDTypesMod            , only : equal_leaf_aclass
  use EDTypesMod            , only : first_leaf_aclass
  use EDTypesMod            , only : nan_leaf_aclass
  use EDTypesMod            , only : max_nleafage
  use EDTypesMod            , only : ican_upper
  use FatesInterfaceMod      , only : hlm_use_planthydro
  use FatesInterfaceMod      , only : hlm_parteh_mode
  use FatesPlantHydraulicsMod, only : FuseCohortHydraulics
  use FatesPlantHydraulicsMod, only : CopyCohortHydraulics
  use FatesPlantHydraulicsMod, only : updateSizeDepTreeHydProps
  use FatesPlantHydraulicsMod, only : initTreeHydStates
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesPlantHydraulicsMod, only : DeallocateHydrCohort
  use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
  use FatesPlantHydraulicsMod, only : UpdateTreeHydrNodes
  use FatesPlantHydraulicsMod, only : UpdateTreeHydrLenVolCond
  use FatesPlantHydraulicsMod, only : SavePreviousCompartmentVolumes
  use FatesPlantHydraulicsMod, only : ConstrainRecruitNumber
  use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
  use FatesAllometryMod  , only : bleaf
  use FatesAllometryMod  , only : bfineroot
  use FatesAllometryMod  , only : bsap_allom
  use FatesAllometryMod  , only : bagw_allom
  use FatesAllometryMod  , only : bbgw_allom
  use FatesAllometryMod  , only : bdead_allom
  use FatesAllometryMod  , only : h_allom
  use FatesAllometryMod  , only : carea_allom
  use FatesAllometryMod  , only : ForceDBH
  use FatesAllometryMod  , only : tree_lai, tree_sai

  use PRTGenericMod,          only : prt_carbon_allom_hyp   
  use PRTGenericMod,          only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,          only : InitPRTVartype
  use PRTGenericMod,          only : prt_vartypes
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
  use PRTGenericMod,          only : SetState

  use PRTAllometricCarbonMod, only : callom_prt_vartypes
  use PRTAllometricCarbonMod, only : ac_bc_inout_id_netdc
  use PRTAllometricCarbonMod, only : ac_bc_in_id_pft
  use PRTAllometricCarbonMod, only : ac_bc_in_id_ctrim
  use PRTAllometricCarbonMod, only : ac_bc_inout_id_dbh

  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)  

  ! CIME globals
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  !
  implicit none
  private
  !
  public :: create_cohort
  public :: zero_cohort
  public :: nan_cohort
  public :: terminate_cohorts
  public :: fuse_cohorts
  public :: insert_cohort
  public :: sort_cohorts
  public :: copy_cohort
  public :: count_cohorts
  public :: InitPRTCohort
  public :: UpdateCohortBioPhysRates
  public :: EvaluateAndCorrectDBH

  logical, parameter :: debug  = .false. ! local debug flag

  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  integer, parameter, private :: conserve_crownarea_and_number_not_dbh = 1
  integer, parameter, private :: conserve_dbh_and_number_not_crownarea = 2

  integer, parameter, private :: cohort_fusion_conservation_method = conserve_crownarea_and_number_not_dbh
  
  ! 10/30/09: Created by Rosie Fisher
  !-------------------------------------------------------------------------------------!

contains

  !-------------------------------------------------------------------------------------!

  subroutine create_cohort(currentSite, patchptr, pft, nn, hite, dbh, bleaf, bfineroot, &
                           bsap, bdead, bstore, laimemory, status, recruitstatus,ctrim, &
                           clayer, spread, leaf_aclass_init, bc_in)

    !
    ! !DESCRIPTION:
    ! create new cohort
    ! There are 4 places this is called
    ! 1) Initializing new cohorts at the beginning of a cold-start simulation
    ! 2) Initializing new recruits during dynamics
    ! 3) Initializing new cohorts at the beginning of a inventory read
    ! 4) Initializing new cohorts during restart
    !
    ! It is assumed that in the first 3, this is called with a reasonable amount of starter information.
    !
    ! !USES:
    !
    ! !ARGUMENTS    

    type(ed_site_type), intent(inout),   target :: currentSite
    type(ed_patch_type), intent(inout), pointer :: patchptr
    integer,  intent(in)   :: pft                        ! Cohort Plant Functional Type
    integer,  intent(in)   :: clayer                     ! canopy status of cohort 
                                                         ! (1 = canopy, 2 = understorey, etc.)
    integer,  intent(in)   :: status                     ! growth status of plant  
                                                         ! (2 = leaves on , 1 = leaves off)
    integer,  intent(in)   :: recruitstatus              ! recruit status of plant  
                                                         ! (1 = recruitment , 0 = other)
    real(r8), intent(in)   :: nn                         ! number of individuals in cohort 
                                                         ! per 'area' (10000m2 default)
    real(r8), intent(in)   :: hite                       ! height: meters
    real(r8), intent(in)   :: dbh                        ! dbh: cm
    real(r8), intent(in)   :: bleaf                      ! biomass in leaves: kgC
    real(r8), intent(in)   :: bfineroot                  ! biomass in fineroots: kgC
    real(r8), intent(in)   :: bsap                       ! biomass in sapwood: kgC
    real(r8), intent(in)   :: bdead                      ! total dead biomass: kGC per indiv
    real(r8), intent(in)   :: bstore                     ! stored carbon: kGC per indiv
    real(r8), intent(in)   :: laimemory                  ! target leaf biomass- set from 
                                                         ! previous year: kGC per indiv
    real(r8), intent(in)   :: ctrim                      ! What is the fraction of the maximum 
                                                         ! leaf biomass that we are targeting?
    real(r8), intent(in)   :: spread                     ! The community assembly effects how 
                                                         ! spread crowns are in horizontal space
    integer,  intent(in)   :: leaf_aclass_init           ! how to initialized the leaf age class
                                                         ! distribution
    integer :: iage                                      ! loop counter for leaf age classes
    type(bc_in_type), intent(in) :: bc_in                ! External boundary conditions
     
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: new_cohort         ! Pointer to New Cohort structure.
    type(ed_cohort_type), pointer :: storesmallcohort 
    type(ed_cohort_type), pointer :: storebigcohort   
    real(r8) :: frac_leaf_aclass(max_nleafage)   ! Fraction of leaves in each age-class
    integer  :: tnull,snull                      ! are the tallest and shortest cohorts allocate
    integer :: nlevsoi_hyd                       ! number of hydraulically active soil layers 

    !----------------------------------------------------------------------

    allocate(new_cohort)

    call nan_cohort(new_cohort)  ! Make everything in the cohort not-a-number
    call zero_cohort(new_cohort) ! Zero things that need to be zeroed. 

    !**********************/
    ! Define cohort state variable
    !**********************/
 
    new_cohort%indexnumber  = fates_unset_int ! Cohort indexing was not thread-safe, setting
                                              ! bogus value for the time being (RGK-012017)

    new_cohort%patchptr     => patchptr

    new_cohort%pft          = pft     
    new_cohort%status_coh   = status
    new_cohort%n            = nn
    new_cohort%hite         = hite
    new_cohort%dbh          = dbh
    new_cohort%canopy_trim  = ctrim
    new_cohort%canopy_layer = clayer
    new_cohort%canopy_layer_yesterday = real(clayer, r8)
    new_cohort%laimemory    = laimemory

    
    ! All newly initialized cohorts start off with an assumption
    ! about leaf age (depending on what is calling the initialization
    ! of this cohort

    if(leaf_aclass_init .eq. equal_leaf_aclass) then
       frac_leaf_aclass(1:nleafage) = 1._r8 / real(nleafage,r8)
    elseif(leaf_aclass_init .eq. first_leaf_aclass) then
       frac_leaf_aclass(1:nleafage) = 0._r8
       frac_leaf_aclass(1)          = 1._r8
    elseif(leaf_aclass_init .eq. nan_leaf_aclass) then
       frac_leaf_aclass(1:nleafage) = nan
    else
       write(fates_log(),*) 'An unknown leaf age distribution was'
       write(fates_log(),*) 'requested during create cohort'
       write(fates_log(),*) 'leaf_aclass_init: ',leaf_aclass_init
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Initialize the Plant allocative Reactive Transport (PaRT) module
    ! Choose from one of the extensible hypotheses (EH)
    ! -----------------------------------------------------------------------------------

    call InitPRTCohort(new_cohort)

    ! The initialization allocates memory, but the boundary and initial
    ! contitions must be set. All new cohorts go through create_cohort()
    ! so this should be the only place this is called.  Alternatively
    ! cohorts can be copied and fused, but special routines handle that.
    ! -----------------------------------------------------------------------------------

    select case(hlm_parteh_mode)
    case (prt_carbon_allom_hyp)

       do iage = 1,nleafage
          call SetState(new_cohort%prt,leaf_organ, carbon12_element, &
                bleaf*frac_leaf_aclass(iage),iage)
       end do
       call SetState(new_cohort%prt,fnrt_organ, carbon12_element, bfineroot)
       call SetState(new_cohort%prt,sapw_organ, carbon12_element, bsap)
       call SetState(new_cohort%prt,store_organ, carbon12_element, bstore)
       call SetState(new_cohort%prt,struct_organ , carbon12_element, bdead)
       call SetState(new_cohort%prt,repro_organ , carbon12_element, 0.0_r8)

    end select
    

    ! This call cycles through the initial conditions, and makes sure that they
    ! are all initialized.
    ! -----------------------------------------------------------------------------------

    call new_cohort%prt%CheckInitialConditions()

    ! This sets things like vcmax25top, that depend on the
    ! leaf age fractions
    call UpdateCohortBioPhysRates(new_cohort)

    call sizetype_class_index(new_cohort%dbh,new_cohort%pft, &
                              new_cohort%size_class,new_cohort%size_by_pft_class)


    ! This routine may be called during restarts, and at this point in the call sequence
    ! the actual cohort data is unknown, as this is really only used for allocation
    ! In these cases, testing if things like biomass are reasonable is pre-mature
    ! However, in this part of the code, we will pass in nominal values for size, number and type
    
    if (new_cohort%dbh <= 0.0_r8 .or. new_cohort%n == 0._r8 .or. new_cohort%pft == 0 ) then
       write(fates_log(),*) 'ED: something is zero in create_cohort', &
                             new_cohort%dbh,new_cohort%n, &
                             new_cohort%pft
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Assign canopy extent and depth
    call carea_allom(new_cohort%dbh,new_cohort%n,spread,new_cohort%pft,new_cohort%c_area)

    new_cohort%treelai = tree_lai(bleaf, new_cohort%pft, new_cohort%c_area,    &
                                  new_cohort%n, new_cohort%canopy_layer,               &
                                  patchptr%canopy_layer_tlai,new_cohort%vcmax25top )    

    new_cohort%treesai = tree_sai(new_cohort%pft, new_cohort%dbh, new_cohort%canopy_trim,   &
                                  new_cohort%c_area, new_cohort%n, new_cohort%canopy_layer, &
                                  patchptr%canopy_layer_tlai, new_cohort%treelai,new_cohort%vcmax25top,2 )  

    new_cohort%lai     = new_cohort%treelai * new_cohort%c_area/patchptr%area


    ! Put cohort at the right place in the linked list
    storebigcohort   => patchptr%tallest
    storesmallcohort => patchptr%shortest 

    if (associated(patchptr%tallest)) then
       tnull = 0
    else
       tnull = 1
       patchptr%tallest => new_cohort
    endif

    if (associated(patchptr%shortest)) then
       snull = 0
    else
       snull = 1
       patchptr%shortest => new_cohort 
    endif

    ! Recuits do not have mortality rates, nor have they moved any
    ! carbon when they are created.  They will bias our statistics
    ! until they have experienced a full day.  We need a newly recruited flag.
    ! This flag will be set to false after it has experienced 
    ! growth, disturbance and mortality.
    new_cohort%isnew = .true.

    if( hlm_use_planthydro.eq.itrue ) then

       nlevsoi_hyd = currentSite%si_hydr%nlevsoi_hyd

       ! This allocates array spaces
       call InitHydrCohort(currentSite,new_cohort)

       ! This calculates node heights
       call UpdateTreeHydrNodes(new_cohort%co_hydr,new_cohort%pft, &
                                new_cohort%hite,nlevsoi_hyd,bc_in)

       ! This calculates volumes, lengths and max conductances
       call UpdateTreeHydrLenVolCond(new_cohort,nlevsoi_hyd,bc_in)
       
       ! Since this is a newly initialized plant, we set the previous compartment-size
       ! equal to the ones we just calculated.
       call SavePreviousCompartmentVolumes(new_cohort%co_hydr)
       
       ! This comes up with starter suctions and then water contents
       ! based on the soil values
       call initTreeHydStates(currentSite,new_cohort, bc_in)

       if(recruitstatus==1)then
          new_cohort%co_hydr%is_newly_recruited = .true.

          ! If plant hydraulics is active, we must constrain the
          ! number density of the new recruits based on the moisture
          ! available to be subsumed in the new plant tissues.
          ! So we go through the process of pre-initializing the hydraulic
          ! states in the temporary cohort, to calculate this new number density

          call ConstrainRecruitNumber(currentSite,new_cohort, bc_in)
       endif

    endif
    
    call insert_cohort(new_cohort, patchptr%tallest, patchptr%shortest, tnull, snull, &
         storebigcohort, storesmallcohort)

    patchptr%tallest  => storebigcohort 
    patchptr%shortest => storesmallcohort

  end subroutine create_cohort

  ! -------------------------------------------------------------------------------------

  subroutine InitPRTCohort(new_cohort)

     ! ----------------------------------------------------------------------------------
     ! This subroutine simply allocates and attaches the correct PRT object.
     ! The call to InitPRTVartype() performs the allocation of the variables
     ! and boundary conditions inside the object.  It also initializes
     ! all values as unitialized (large bogus values).
     !
     ! Each PARTEH allocation hypothesis has different expectations of boundary conditions.
     ! These are specified by pointers to values in the host model. Because these
     ! are pointers, they just need to be set once when the prt object is first initalized.
     ! The calls below to "RegisterBCINOut", "RegisterBCIn" and "RegisterBCOut" are
     ! setting those pointers.
     ! -----------------------------------------------------------------------------------

     !
     ! !ARGUMENTS    
     type(ed_cohort_type), intent(inout), target  :: new_cohort
     type(callom_prt_vartypes), pointer :: callom_prt


     ! Allocate the PRT class object
     ! Each hypothesis has a different object which is an extension
     ! of the base class.

     select case(hlm_parteh_mode)
     case (prt_carbon_allom_hyp)
        
        allocate(callom_prt)
        new_cohort%prt => callom_prt
     
     case DEFAULT

        write(fates_log(),*) 'You specified an unknown PRT module'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))

     end select
     
     ! This is the call to allocate the data structures in the PRT object
     ! This call will be extended to each specific class.

     call new_cohort%prt%InitPRTVartype()


     ! Set the boundary conditions that flow in an out of the PARTEH
     ! allocation hypotheses.  These are pointers in the PRT objects that
     ! point to values outside in the FATES model.

     ! Example:
     ! "ac_bc_inout_id_dbh" is the unique integer that defines the object index
     ! for the allometric carbon "ac" boundary condition "bc" for DBH "dbh"
     ! that is classified as input and output "inout".
     ! See PRTAllometricCarbonMod.F90 to track its usage.
     ! bc_rval is used as the optional argument identifyer to specify a real
     ! value boundary condition.
     ! bc_ival is used as the optional argument identifyer to specify an integer
     ! value boundary condition.
     

     select case(hlm_parteh_mode)
     case (prt_carbon_allom_hyp)

        ! Register boundary conditions for the Carbon Only Allometric Hypothesis

        call new_cohort%prt%RegisterBCInOut(ac_bc_inout_id_dbh,bc_rval = new_cohort%dbh)
        call new_cohort%prt%RegisterBCInOut(ac_bc_inout_id_netdc,bc_rval = new_cohort%npp_acc)
        call new_cohort%prt%RegisterBCIn(ac_bc_in_id_pft,bc_ival = new_cohort%pft)
        call new_cohort%prt%RegisterBCIn(ac_bc_in_id_ctrim,bc_rval = new_cohort%canopy_trim)

     end select


     return
  end subroutine InitPRTCohort

  !-------------------------------------------------------------------------------------!

  subroutine nan_cohort(cc_p)
    !
    ! !DESCRIPTION:
    !  Make all the cohort variables NaN so they aren't used before defined.   
    !
    ! !USES:

    use FatesConstantsMod, only : fates_unset_int

    !
    ! !ARGUMENTS    
    type (ed_cohort_type), intent(inout), target  :: cc_p
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type)   , pointer             :: currentCohort
    !----------------------------------------------------------------------

    currentCohort => cc_p

    currentCohort%taller      => null()       ! pointer to next tallest cohort     
    currentCohort%shorter     => null()       ! pointer to next shorter cohort     
    currentCohort%patchptr    => null()       ! pointer to patch that cohort is in

    nullify(currentCohort%taller) 
    nullify(currentCohort%shorter) 
    nullify(currentCohort%patchptr) 

    ! VEGETATION STRUCTURE
    currentCohort%pft                = fates_unset_int  ! pft number                           
    currentCohort%indexnumber        = fates_unset_int  ! unique number for each cohort. (within clump?)
    currentCohort%canopy_layer       = fates_unset_int  ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)   
    currentCohort%canopy_layer_yesterday       = nan  ! recent canopy status of cohort (1 = canopy, 2 = understorey, etc.)   
    currentCohort%NV                 = fates_unset_int  ! Number of leaf layers: -
    currentCohort%status_coh         = fates_unset_int  ! growth status of plant  (2 = leaves on , 1 = leaves off)
    currentCohort%size_class         = fates_unset_int  ! size class index
    currentCohort%size_class_lasttimestep = fates_unset_int  ! size class index
    currentCohort%size_by_pft_class  = fates_unset_int  ! size by pft classification index

    currentCohort%n                  = nan ! number of individuals in cohort per 'area' (10000m2 default)     
    currentCohort%dbh                = nan ! 'diameter at breast height' in cm                            
    currentCohort%hite               = nan ! height: meters                   
    currentCohort%laimemory          = nan ! target leaf biomass- set from previous year: kGC per indiv
    currentCohort%lai                = nan ! leaf area index of cohort   m2/m2      
    currentCohort%sai                = nan ! stem area index of cohort   m2/m2
    currentCohort%g_sb_laweight      = nan ! Total leaf conductance of cohort (stomata+blayer) weighted by leaf-area [m/s]*[m2]
    currentCohort%canopy_trim        = nan ! What is the fraction of the maximum leaf biomass that we are targeting? :-
    currentCohort%leaf_cost          = nan ! How much does it cost to maintain leaves: kgC/m2/year-1
    currentCohort%excl_weight        = nan ! How much of this cohort is demoted each year, as a proportion of all cohorts:-
    currentCohort%prom_weight        = nan ! How much of this cohort is promoted each year, as a proportion of all cohorts:-
    currentCohort%c_area             = nan ! areal extent of canopy (m2)
    currentCohort%treelai            = nan ! lai of tree (total leaf area (m2) / canopy area (m2)
    currentCohort%treesai            = nan ! stem area index of tree (total stem area (m2) / canopy area (m2)

    currentCohort%vcmax25top = nan 
    currentCohort%jmax25top  = nan 
    currentCohort%tpu25top   = nan 
    currentCohort%kp25top    = nan 

    ! CARBON FLUXES 
    currentCohort%gpp_acc_hold       = nan ! GPP:  kgC/indiv/year
    currentCohort%gpp_tstep          = nan ! GPP:  kgC/indiv/timestep
    currentCohort%gpp_acc            = nan ! GPP:  kgC/indiv/day         
    currentCohort%npp_acc_hold       = nan ! NPP:  kgC/indiv/year
    currentCohort%npp_tstep          = nan ! NPP:  kGC/indiv/timestep
    currentCohort%npp_acc            = nan ! NPP:  kgC/indiv/day  
    currentCohort%year_net_uptake(:) = nan ! Net uptake of individual leaf layers kgC/m2/year
    currentCohort%ts_net_uptake(:)   = nan ! Net uptake of individual leaf layers kgC/m2/s
    currentCohort%resp_acc_hold      = nan ! RESP: kgC/indiv/year
    currentCohort%resp_tstep         = nan ! RESP: kgC/indiv/timestep
    currentCohort%resp_acc           = nan ! RESP: kGC/cohort/day
    
    currentCohort%c13disc_clm        = nan ! C13 discrimination, per mil at indiv/timestep
    currentCohort%c13disc_acc        = nan ! C13 discrimination, per mil at indiv/timestep at indiv/daily at the end of a day

    !RESPIRATION
    currentCohort%rdark              = nan
    currentCohort%resp_m             = nan ! Maintenance respiration.  kGC/cohort/year
    currentCohort%resp_g             = nan ! Growth respiration.       kGC/cohort/year
    currentCohort%livestem_mr        = nan ! Live stem maintenance respiration. kgC/indiv/s-1 
    currentCohort%livecroot_mr       = nan ! Coarse root maintenance respiration. kgC/indiv/s-1 
    currentCohort%froot_mr           = nan ! Fine root maintenance respiration. kgC/indiv/s-1 

    ! ALLOCATION
    currentCohort%dmort              = nan ! proportional mortality rate. (year-1)

    ! logging
    currentCohort%lmort_direct       = nan
    currentCohort%lmort_infra        = nan
    currentCohort%lmort_collateral   = nan
    currentCohort%l_degrad           = nan


    currentCohort%seed_prod          = nan ! reproduction seed and clonal: KgC/indiv/year
    currentCohort%c_area             = nan ! areal extent of canopy (m2)
    currentCohort%treelai            = nan ! lai of tree (total leaf area (m2) / canopy area (m2)
    currentCohort%treesai            = nan ! stem area index of tree (total stem area (m2) / canopy area (m2)


    ! VARIABLES NEEDED FOR INTEGRATION 
    currentCohort%dndt               = nan ! time derivative of cohort size 
    currentCohort%dhdt               = nan ! time derivative of height 
    currentCohort%ddbhdt             = nan ! time derivative of dbh 

    ! FIRE
    currentCohort%fraction_crown_burned = nan ! proportion of crown affected by fire
    currentCohort%cambial_mort          = nan ! probability that trees dies due to cambial char P&R (1986)
    currentCohort%crownfire_mort        = nan ! probability of tree post-fire mortality due to crown scorch
    currentCohort%fire_mort             = nan ! post-fire mortality from cambial and crown damage assuming two are independent

  end subroutine nan_cohort

  !-------------------------------------------------------------------------------------!

  subroutine zero_cohort(cc_p)
    !
    ! !DESCRIPTION:
    ! Zero variables that need to be accounted for if 
    ! this cohort is altered before they are defined.       
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type (ed_cohort_type), intent(inout), target  :: cc_p
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type)   , pointer             :: currentCohort
    !----------------------------------------------------------------------

    currentCohort => cc_p

    currentCohort%NV                 = 0    
    currentCohort%status_coh         = 0    
    currentCohort%rdark              = 0._r8
    currentCohort%resp_m             = 0._r8 
    currentCohort%resp_g             = 0._r8
    currentCohort%livestem_mr        = 0._r8
    currentCohort%livecroot_mr       = 0._r8
    currentCohort%froot_mr           = 0._r8
    currentCohort%fire_mort          = 0._r8 
    currentcohort%npp_acc            = 0._r8
    currentcohort%gpp_acc            = 0._r8
    currentcohort%resp_acc           = 0._r8
    currentcohort%npp_tstep          = 0._r8
    currentcohort%gpp_tstep          = 0._r8
    currentcohort%resp_tstep         = 0._r8
    currentcohort%resp_acc_hold      = 0._r8

    currentcohort%year_net_uptake(:) = 999._r8 ! this needs to be 999, or trimming of new cohorts will break. 
    currentcohort%ts_net_uptake(:)   = 0._r8
    currentcohort%seed_prod          = 0._r8
    currentcohort%fraction_crown_burned = 0._r8 
    currentCohort%size_class            = 1
    currentCohort%size_class_lasttimestep = 0
    currentcohort%npp_acc_hold       = 0._r8 
    currentcohort%gpp_acc_hold       = 0._r8  
    currentcohort%dmort              = 0._r8 
    currentcohort%g_sb_laweight      = 0._r8 
    currentcohort%treesai            = 0._r8  
    currentCohort%lmort_direct       = 0._r8
    currentCohort%lmort_infra        = 0._r8
    currentCohort%lmort_collateral   = 0._r8
    currentCohort%l_degrad           = 0._r8    
    currentCohort%leaf_cost          = 0._r8
    currentcohort%excl_weight        = 0._r8
    currentcohort%prom_weight        = 0._r8
    currentcohort%crownfire_mort     = 0._r8
    currentcohort%cambial_mort       = 0._r8
    currentCohort%c13disc_clm        = 0._r8 
    currentCohort%c13disc_acc        = 0._r8
    
  end subroutine zero_cohort

  !-------------------------------------------------------------------------------------!
  subroutine terminate_cohorts( currentSite, currentPatch, level )
    !
    ! !DESCRIPTION:
    ! terminates cohorts when they get too small      
    !
    ! !USES:
    use SFParamsMod, only : SF_val_CWD_frac
    !
    ! !ARGUMENTS    
    type (ed_site_type) , intent(inout), target :: currentSite
    type (ed_patch_type), intent(inout), target :: currentPatch
    integer             , intent(in)            :: level

    ! Important point regarding termination levels.  Termination is typically
    ! called after fusion.  We do this so that we can re-capture the biomass that would
    ! otherwise be lost from termination.  The biomass of a fused plant remains in the
    ! live pool.  However, some plant number densities can be so low that they 
    ! can cause numerical instabilities.  Thus, we call terminate_cohorts at level=1
    ! before fusion to get rid of these cohorts that are so incredibly sparse, and then
    ! terminate the remainder at level 2 for various other reasons.

    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_cohort_type) , pointer :: shorterCohort
    type (ed_cohort_type) , pointer :: tallerCohort

    real(r8) :: leaf_c    ! leaf carbon [kg]
    real(r8) :: store_c   ! storage carbon [kg]
    real(r8) :: sapw_c    ! sapwood carbon [kg]
    real(r8) :: fnrt_c    ! fineroot carbon [kg]
    real(r8) :: repro_c   ! reproductive carbon [kg]
    real(r8) :: struct_c  ! structural carbon [kg]

    integer :: terminate   ! do we terminate (1) or not (0) 
    integer :: c           ! counter for litter size class. 
    integer :: levcan      ! canopy level
    !----------------------------------------------------------------------


    currentCohort => currentPatch%shortest
    do while (associated(currentCohort))

       terminate = 0 
       tallerCohort => currentCohort%taller

       leaf_c  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
       store_c = currentCohort%prt%GetState(store_organ, all_carbon_elements)
       sapw_c  = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
       fnrt_c  = currentCohort%prt%GetState(fnrt_organ, all_carbon_elements)
       struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
       repro_c  = currentCohort%prt%GetState(repro_organ, all_carbon_elements)

       ! Check if number density is so low is breaks math (level 1)
       if (currentcohort%n <  min_n_safemath .and. level == 1) then
         terminate = 1
	 if ( debug ) then
             write(fates_log(),*) 'terminating cohorts 0',currentCohort%n/currentPatch%area,currentCohort%dbh
         endif
       endif

       ! The rest of these are only allowed if we are not dealing with a recruit (level 2)
       if (.not.currentCohort%isnew .and. level == 2) then

         ! Not enough n or dbh
         if  (currentCohort%n/currentPatch%area <= min_npm2 .or.	&  !
              currentCohort%n <= min_nppatch .or. &
              (currentCohort%dbh < 0.00001_r8 .and. store_c < 0._r8) ) then 
            terminate = 1

            if ( debug ) then
               write(fates_log(),*) 'terminating cohorts 1',currentCohort%n/currentPatch%area,currentCohort%dbh
            endif
         endif

         ! Outside the maximum canopy layer
         if (currentCohort%canopy_layer > nclmax ) then 
           terminate = 1
           if ( debug ) then
             write(fates_log(),*) 'terminating cohorts 2', currentCohort%canopy_layer
           endif
         endif

         ! live biomass pools are terminally depleted
         if ( ( sapw_c+leaf_c+fnrt_c ) < 1e-10_r8  .or.  &
               store_c  < 1e-10_r8) then 
            terminate = 1  
            if ( debug ) then
              write(fates_log(),*) 'terminating cohorts 3', &
                    sapw_c,leaf_c,fnrt_c,store_c
            endif
         endif

         ! Total cohort biomass is negative
         if ( ( struct_c+sapw_c+leaf_c+fnrt_c+store_c ) < 0._r8) then
            terminate = 1
            if ( debug ) then
            write(fates_log(),*) 'terminating cohorts 4', & 
                  struct_c,sapw_c,leaf_c,fnrt_c,store_c

         endif

         endif 
      endif    !  if (.not.currentCohort%isnew .and. level == 2) then

      if (terminate == 1) then 
         
          ! preserve a record of the to-be-terminated cohort for mortality accounting
          levcan = currentCohort%canopy_layer

          if( hlm_use_planthydro == itrue ) &
             call AccumulateMortalityWaterStorage(currentSite,currentCohort,currentCohort%n)

          if(levcan==ican_upper) then
             currentSite%term_nindivs_canopy(currentCohort%size_class,currentCohort%pft) = &
                   currentSite%term_nindivs_canopy(currentCohort%size_class,currentCohort%pft) + currentCohort%n
 
             currentSite%term_carbonflux_canopy = currentSite%term_carbonflux_canopy + &
                   currentCohort%n * (struct_c+sapw_c+leaf_c+fnrt_c+store_c+repro_c)
          else
             currentSite%term_nindivs_ustory(currentCohort%size_class,currentCohort%pft) = &
                   currentSite%term_nindivs_ustory(currentCohort%size_class,currentCohort%pft) + currentCohort%n
 
             currentSite%term_carbonflux_ustory = currentSite%term_carbonflux_ustory + &
                   currentCohort%n * (struct_c+sapw_c+leaf_c+fnrt_c+store_c+repro_c)
          end if

          !put the litter from the terminated cohorts straight into the fragmenting pools
          if (currentCohort%n.gt.0.0_r8) then
             do c=1,ncwd

                currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + currentCohort%n*(struct_c+sapw_c) / &
                     currentPatch%area &
                     * SF_val_CWD_frac(c) * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) 
                currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + currentCohort%n*(struct_c+sapw_c) / &
                     currentPatch%area &
                     * SF_val_CWD_frac(c) * (1.0_r8 -  EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) 
             enddo
             
             currentPatch%leaf_litter(currentCohort%pft) = currentPatch%leaf_litter(currentCohort%pft) + currentCohort%n* &
                  (leaf_c)/currentPatch%area

             currentPatch%root_litter(currentCohort%pft) = currentPatch%root_litter(currentCohort%pft) + currentCohort%n* &
                  (fnrt_c+store_c)/currentPatch%area 


             ! keep track of the above fluxes at the site level as a CWD/litter input flux (in kg / site-m2 / yr)
             do c=1,ncwd
                currentSite%CWD_AG_diagnostic_input_carbonflux(c)  = currentSite%CWD_AG_diagnostic_input_carbonflux(c) &
                     + currentCohort%n*(struct_c + sapw_c) * &
                     SF_val_CWD_frac(c) * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * hlm_days_per_year / AREA
                currentSite%CWD_BG_diagnostic_input_carbonflux(c)  = currentSite%CWD_BG_diagnostic_input_carbonflux(c) &
                     + currentCohort%n*(struct_c + sapw_c) * &
                     SF_val_CWD_frac(c) * (1.0_r8 -  EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))  * hlm_days_per_year / AREA
             enddo
             
             currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                  currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) +  &
                  currentCohort%n * (leaf_c) * hlm_days_per_year  / AREA
             currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                  currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) + &
                  currentCohort%n * (fnrt_c + store_c) * hlm_days_per_year  / AREA

          end if
          
          ! Zero out the state pools
          call SetState(currentCohort%prt,leaf_organ,carbon12_element,0.0_r8)
          call SetState(currentCohort%prt,fnrt_organ,carbon12_element,0.0_r8)
          call SetState(currentCohort%prt,sapw_organ,carbon12_element,0.0_r8)
          call SetState(currentCohort%prt,struct_organ,carbon12_element,0.0_r8)
          call SetState(currentCohort%prt,repro_organ,carbon12_element,0.0_r8)
          call SetState(currentCohort%prt,store_organ,carbon12_element,0.0_r8)

          ! Set pointers and remove the current cohort from the list
          shorterCohort => currentCohort%shorter
          
          if (.not. associated(tallerCohort)) then
             currentPatch%tallest => shorterCohort
             if(associated(shorterCohort)) shorterCohort%taller => null()
          else 
             tallerCohort%shorter => shorterCohort

          endif
          
          if (.not. associated(shorterCohort)) then
             currentPatch%shortest => tallerCohort
             if(associated(tallerCohort)) tallerCohort%shorter => null()
          else 
             shorterCohort%taller => tallerCohort
          endif
          
          ! At this point, nothing should be pointing to current Cohort
          if (hlm_use_planthydro.eq.itrue) call DeallocateHydrCohort(currentCohort)

          ! Deallocate the cohort's PRT structure
          call currentCohort%prt%DeallocatePRTVartypes()
          deallocate(currentCohort%prt)

          deallocate(currentCohort)
          nullify(currentCohort)
          
       endif
       currentCohort => tallerCohort
    enddo

  end subroutine terminate_cohorts

  !-------------------------------------------------------------------------------------!

  subroutine fuse_cohorts(currentSite, currentPatch, bc_in)  

     !
     ! !DESCRIPTION:
     ! Join similar cohorts to reduce total number            
     !
     ! !USES:
     use EDParamsMod , only :  ED_val_cohort_fusion_tol
     !
     ! !ARGUMENTS   
     type (ed_site_type), intent(inout),  target :: currentSite 
     type (ed_patch_type), intent(inout), target :: currentPatch
     type (bc_in_type), intent(in)               :: bc_in
     !

     ! !LOCAL VARIABLES:
     type (ed_cohort_type) , pointer :: currentCohort
     type (ed_cohort_type) , pointer :: nextc
     type (ed_cohort_type) , pointer :: nextnextc

     type (ed_cohort_type) , pointer :: shorterCohort
     type (ed_cohort_type) , pointer :: tallerCohort

     integer  :: i  
     integer  :: fusion_took_place
     integer  :: iterate    ! do we need to keep fusing to get below maxcohorts?
     integer  :: nocohorts
     real(r8) :: newn
     real(r8) :: diff
     real(r8) :: leaf_c_next   ! Leaf carbon * plant density of current (for weighting)
     real(r8) :: leaf_c_curr   ! Leaf carbon * plant density of next (for weighting)
     real(r8) :: leaf_c_target 
     real(r8) :: dynamic_fusion_tolerance
     real(r8) :: dbh
     real(r8) :: leaf_c             ! leaf carbon [kg]

     integer  :: largersc, smallersc, sc_i        ! indices for tracking the growth flux caused by fusion
     real(r8) :: larger_n, smaller_n

     logical, parameter :: fuse_debug = .false.   ! This debug is over-verbose
                                                 ! and gets its own flag

     !----------------------------------------------------------------------

     !set initial fusion tolerance
     dynamic_fusion_tolerance = ED_val_cohort_fusion_tol

     !This needs to be a function of the canopy layer, because otherwise, at canopy closure
     !the number of cohorts doubles and very dissimilar cohorts are fused together
     !because c_area and biomass are non-linear with dbh, this causes several mass inconsistancies
     !in theory, all of this routine therefore causes minor losses of C and area, but these are below 
     !detection limit normally. 

     iterate = 1
     fusion_took_place = 0   

     !---------------------------------------------------------------------!
     !  Keep doing this until nocohorts <= maxcohorts                         !
     !---------------------------------------------------------------------!
     
     if (associated(currentPatch%shortest)) then  
        do while(iterate == 1)
           
           currentCohort => currentPatch%tallest
           
           ! The following logic continues the loop while the current cohort is not the shortest cohort
           ! if they point to the same target (ie equivalence), then the loop ends.
           ! This loop is different than the simple "continue while associated" loop in that
           ! it omits the last cohort (because it has already been compared by that point)
           
           do while ( .not.associated(currentCohort,currentPatch%shortest) )

              nextc => currentPatch%tallest

              do while (associated(nextc))
                 nextnextc => nextc%shorter
                 diff = abs((currentCohort%dbh - nextc%dbh)/(0.5*(currentCohort%dbh + nextc%dbh)))  

                 !Criteria used to divide up the height continuum into different cohorts.

                 if (diff < dynamic_fusion_tolerance) then

                    ! Don't fuse a cohort with itself!
                    if (.not.associated(currentCohort,nextc) ) then

                       if (currentCohort%pft == nextc%pft) then              

                          ! check cohorts in same c. layer. before fusing

                          if (currentCohort%canopy_layer == nextc%canopy_layer) then 

                             ! Note: because newly recruited cohorts that have not experienced
                             ! a day yet will have un-known flux quantities or change rates
                             ! we don't want them fusing with non-new cohorts.  We allow them
                             ! to fuse with other new cohorts to keep the total number of cohorts
                             ! down.

                             if( currentCohort%isnew.eqv.nextc%isnew ) then

                                newn = currentCohort%n + nextc%n

                                fusion_took_place = 1         

                                if ( fuse_debug .and. currentCohort%isnew ) then
                                   write(fates_log(),*) 'Fusing Two Cohorts'
                                   write(fates_log(),*) 'newn: ',newn
                                   write(fates_log(),*) 'Cohort I, Cohort II' 
                                   write(fates_log(),*) 'n:',currentCohort%n,nextc%n
                                   write(fates_log(),*) 'isnew:',currentCohort%isnew,nextc%isnew
                                   write(fates_log(),*) 'laimemory:',currentCohort%laimemory,nextc%laimemory
                                   write(fates_log(),*) 'hite:',currentCohort%hite,nextc%hite
                                   write(fates_log(),*) 'dbh:',currentCohort%dbh,nextc%dbh
                                   write(fates_log(),*) 'pft:',currentCohort%pft,nextc%pft
                                   write(fates_log(),*) 'canopy_trim:',currentCohort%canopy_trim,nextc%canopy_trim
                                   write(fates_log(),*) 'canopy_layer_yesterday:', &
                                         currentCohort%canopy_layer_yesterday,nextc%canopy_layer_yesterday
                                   do i=1, nlevleaf
                                      write(fates_log(),*) 'leaf level: ',i,'year_net_uptake', &
                                            currentCohort%year_net_uptake(i),nextc%year_net_uptake(i)
                                   end do
                                end if


                                ! Fuse all mass pools
                                call currentCohort%prt%WeightedFusePRTVartypes(nextc%prt, &
                                                                               currentCohort%n/newn )

                                ! Leaf biophysical rates (use leaf mass weighting)
                                ! -----------------------------------------------------------------
                                call UpdateCohortBioPhysRates(currentCohort)

                                currentCohort%laimemory   = (currentCohort%n*currentCohort%laimemory   &
                                      + nextc%n*nextc%laimemory)/newn

                                currentCohort%canopy_trim = (currentCohort%n*currentCohort%canopy_trim &
                                      + nextc%n*nextc%canopy_trim)/newn
				
                                ! c13disc_acc calculation; weighted mean by GPP
                                if ((currentCohort%n * currentCohort%gpp_acc + nextc%n * nextc%gpp_acc) .eq. 0.0_r8) then
                                    currentCohort%c13disc_acc = 0.0_r8
                                else  
                                    currentCohort%c13disc_acc = (currentCohort%n * currentCohort%gpp_acc * currentCohort%c13disc_acc +   &
                                          nextc%n * nextc%gpp_acc * nextc%c13disc_acc)/    &
                                          (currentCohort%n * currentCohort%gpp_acc + nextc%n * nextc%gpp_acc)
                                endif
                                
                                select case(cohort_fusion_conservation_method)
                                   !
                                   ! -----------------------------------------------------------------
                                   ! Because cohort fusion is an unavoidable but non-physical process,
                                   ! and because of the various nonlinear allometric relationships,
                                   ! it isn't possible to simultaneously conserve all of the allometric
                                   ! relationships during cohort fusion.  We will always conserve carbon,
                                   ! but there are choices to made about what else to conserve or not.
                                   ! In particular, there is a choice to be made of conservation amongst
                                   ! the number density, stem diameter, and crown area. Below,
                                   ! some different conservation relationships can be chosen during fusion.
                                   ! -----------------------------------------------------------------
                                   !
                                case(conserve_crownarea_and_number_not_dbh)
                                   !
                                   ! -----------------------------------------------------------------
                                   ! conserve total crown area during the fusion step, and then calculate
                                   ! dbh of the fused cohort as that which conserves both crown area and
                                   ! the dbh to crown area allometry.  dbh will be updated in the next
                                   ! growth step in the (likely) event that dbh to structural iomass
                                   ! allometry is exceeded. if using a capped crown area allometry and
                                   ! above the cap, then calculate as the weighted average of fusing
                                   ! cohorts' dbh
                                   ! -----------------------------------------------------------------
                                   !
                                   call carea_allom(currentCohort%dbh,currentCohort%n, &
                                         currentSite%spread,currentCohort%pft,&
                                         currentCohort%c_area,inverse=.false.)
                                   
                                   call carea_allom(nextc%dbh,nextc%n, &
                                         currentSite%spread,nextc%pft,&
                                         nextc%c_area,inverse=.false.)
                                   
                                   currentCohort%c_area = currentCohort%c_area + nextc%c_area

                                   !
                                   call carea_allom(dbh,newn,currentSite%spread,currentCohort%pft,&
                                        currentCohort%c_area,inverse=.true.)
                                   !
                                   if (abs(dbh-fates_unset_r8)<nearzero) then
                                      currentCohort%dbh = (currentCohort%n*currentCohort%dbh         &
                                           + nextc%n*nextc%dbh)/newn

                                      if( EDPftvarcon_inst%woody(currentCohort%pft) == itrue ) then

                                          call ForceDBH( currentCohort%pft, currentCohort%canopy_trim, &
                                               currentCohort%dbh, currentCohort%hite, &
                                               bdead = currentCohort%prt%GetState(struct_organ,all_carbon_elements))

                                      end if
                                      !
                                      call carea_allom(currentCohort%dbh,newn,currentSite%spread,currentCohort%pft,&
                                            currentCohort%c_area,inverse=.false.)
                                      
                                   else
                                      currentCohort%dbh = dbh
                                   endif

                                   !
                                   call h_allom(currentCohort%dbh,currentCohort%pft,currentCohort%hite)
                                   !
                                case(conserve_dbh_and_number_not_crownarea)
                                   !
                                   ! -----------------------------------------------------------------
                                   ! Here we conserve the mean stem diameter of the trees in the cohorts
                                   ! rather than the crown area of the cohort
                                   ! -----------------------------------------------------------------
                                   !
                                   currentCohort%dbh         = (currentCohort%n*currentCohort%dbh         &
                                        + nextc%n*nextc%dbh)/newn
                                   !
                                   call h_allom(currentCohort%dbh,currentCohort%pft,currentCohort%hite)
                                   !
                                   ! -----------------------------------------------------------------
                                   ! If fusion pushed structural biomass to be larger than
                                   ! the allometric target value derived by diameter, we
                                   ! then increase diameter and height until the allometric
                                   ! target matches actual bdead. (if it is the other way around
                                   ! we then just let the carbon pools grow to fill out allometry)
                                   ! -----------------------------------------------------------------
                                   !
                                   if( EDPftvarcon_inst%woody(currentCohort%pft) == itrue ) then
                                      call ForceDBH( currentCohort%pft, currentCohort%canopy_trim, &
                                           currentCohort%dbh, currentCohort%hite, &
                                           bdead = currentCohort%prt%GetState(struct_organ,all_carbon_elements))

                                   end if
                                   !
                                   call carea_allom(currentCohort%dbh,newn,currentSite%spread,currentCohort%pft,&
                                        currentCohort%c_area,inverse=.false.)
                                   !
                                case default
                                    write(fates_log(),*) 'FATES: Invalid choice for cohort_fusion_conservation_method'
                                   call endrun(msg=errMsg(sourcefile, __LINE__))
                                end select


                                ! If fusion forces the actual leaf biomass to be unreasonably
                                ! greater than the target (ie 25%), reset the DBH

!                                call bleaf(currentCohort%dbh,currentCohort%pft, &
!                                     currentCohort%canopy_trim,leaf_c_target)
                                
                                leaf_c = currentCohort%prt%GetState(leaf_organ,all_carbon_elements)

!                                if (leaf_c > leaf_c_target*1.25_r8) then
!                                   call ForceDBH( currentCohort%pft, currentCohort%canopy_trim, &
!                                        currentCohort%dbh, currentCohort%hite, &
!                                        bl = leaf_c)
!                                   call carea_allom(currentCohort%dbh,newn,currentSite%spread,currentCohort%pft, &
!                                        currentCohort%c_area,inverse=.false.)
!                                end if


                                currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, newn, &
                                               currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, &
                                               currentCohort%vcmax25top)
                                currentCohort%treesai = tree_sai(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_trim, &
                                               currentCohort%c_area, newn, currentCohort%canopy_layer, &
                                               currentPatch%canopy_layer_tlai, currentCohort%treelai,currentCohort%vcmax25top,1 ) 



                                call sizetype_class_index(currentCohort%dbh,currentCohort%pft, &
                                      currentCohort%size_class,currentCohort%size_by_pft_class)
				      

                                if(hlm_use_planthydro.eq.itrue) then			  					  				  
                                    call FuseCohortHydraulics(currentSite,currentCohort,nextc,bc_in,newn)				    
                                endif

                                ! recent canopy history
                                currentCohort%canopy_layer_yesterday  = (currentCohort%n*currentCohort%canopy_layer_yesterday  + &
                                      nextc%n*nextc%canopy_layer_yesterday)/newn


                                ! keep track of the size class bins so that we can monitor growth fluxes
                                ! compare the values.  if they are the same, then nothing needs to be done. if not, track the diagnostic flux
                                if (currentCohort%size_class_lasttimestep .ne. nextc%size_class_lasttimestep ) then
                                   !
                                   ! keep track of which was which, irresespective of which cohort they were in
                                   if (currentCohort%size_class_lasttimestep .gt. nextc%size_class_lasttimestep) then
                                      largersc = currentCohort%size_class_lasttimestep
                                      smallersc = nextc%size_class_lasttimestep
                                      larger_n = currentCohort%n
                                      smaller_n = nextc%n
                                   else
                                      largersc = nextc%size_class_lasttimestep
                                      smallersc = currentCohort%size_class_lasttimestep
                                      larger_n = nextc%n
                                      smaller_n = currentCohort%n
                                   endif
                                   !
                                   ! it is possible that fusion has caused cohorts separated by at least two size bin deltas to join.  
                                   ! so slightly complicated to keep track of because the resulting cohort could be in one of the old bins or in between
                                   ! structure as a loop to handle the general case
                                   !
                                   ! first the positive growth case
                                   do sc_i = smallersc + 1, currentCohort%size_class
                                      currentSite%growthflux_fusion(sc_i, currentCohort%pft) = &
                                           currentSite%growthflux_fusion(sc_i, currentCohort%pft) + smaller_n
                                   end do
                                   !
                                   ! next the negative growth case
                                   do sc_i = currentCohort%size_class + 1, largersc
                                      currentSite%growthflux_fusion(sc_i, currentCohort%pft) = &
                                           currentSite%growthflux_fusion(sc_i, currentCohort%pft) - larger_n
                                   end do
                                   ! now that we've tracked the change flux.  reset the memory of the prior timestep
                                   currentCohort%size_class_lasttimestep = currentCohort%size_class
                                endif
                                   
                                ! Flux and biophysics variables have not been calculated for recruits we just default to 
                                ! their initization values, which should be the same for eahc
                                
                                if ( .not.currentCohort%isnew) then

                                   currentCohort%seed_prod      = (currentCohort%n*currentCohort%seed_prod + &
                                         nextc%n*nextc%seed_prod)/newn
                                   currentCohort%gpp_acc        = (currentCohort%n*currentCohort%gpp_acc     + &
                                         nextc%n*nextc%gpp_acc)/newn
                                   currentCohort%npp_acc        = (currentCohort%n*currentCohort%npp_acc     + &
                                         nextc%n*nextc%npp_acc)/newn
                                   currentCohort%resp_acc       = (currentCohort%n*currentCohort%resp_acc    + &
                                         nextc%n*nextc%resp_acc)/newn
                                   currentCohort%resp_acc_hold  = &
                                         (currentCohort%n*currentCohort%resp_acc_hold + &
                                         nextc%n*nextc%resp_acc_hold)/newn
                                   currentCohort%npp_acc_hold   = &
                                         (currentCohort%n*currentCohort%npp_acc_hold + &
                                         nextc%n*nextc%npp_acc_hold)/newn
                                   currentCohort%gpp_acc_hold   = &
                                         (currentCohort%n*currentCohort%gpp_acc_hold + &
                                         nextc%n*nextc%gpp_acc_hold)/newn

                                   currentCohort%dmort          = (currentCohort%n*currentCohort%dmort       + &
                                         nextc%n*nextc%dmort)/newn

                                   currentCohort%fire_mort      = (currentCohort%n*currentCohort%fire_mort   + &
                                         nextc%n*nextc%fire_mort)/newn

                                   ! mortality diagnostics
                                   currentCohort%cmort = (currentCohort%n*currentCohort%cmort + nextc%n*nextc%cmort)/newn
                                   currentCohort%hmort = (currentCohort%n*currentCohort%hmort + nextc%n*nextc%hmort)/newn
                                   currentCohort%bmort = (currentCohort%n*currentCohort%bmort + nextc%n*nextc%bmort)/newn
                                   currentCohort%frmort = (currentCohort%n*currentCohort%frmort + nextc%n*nextc%frmort)/newn

                                   ! logging mortality, Yi Xu
                                   currentCohort%lmort_direct = (currentCohort%n*currentCohort%lmort_direct + &
                                         nextc%n*nextc%lmort_direct)/newn
                                   currentCohort%lmort_collateral = (currentCohort%n*currentCohort%lmort_collateral + &
                                         nextc%n*nextc%lmort_collateral)/newn
                                   currentCohort%lmort_infra = (currentCohort%n*currentCohort%lmort_infra + &
                                         nextc%n*nextc%lmort_infra)/newn
                                   currentCohort%l_degrad = (currentCohort%n*currentCohort%l_degrad + &
                                         nextc%n*nextc%l_degrad)/newn
                                   
                                   ! biomass and dbh tendencies
                                   currentCohort%ddbhdt     = (currentCohort%n*currentCohort%ddbhdt  + &
                                         nextc%n*nextc%ddbhdt)/newn

                                   do i=1, nlevleaf     
                                      if (currentCohort%year_net_uptake(i) == 999._r8 .or. nextc%year_net_uptake(i) == 999._r8) then
                                         currentCohort%year_net_uptake(i) = &
                                               min(nextc%year_net_uptake(i),currentCohort%year_net_uptake(i))
                                      else
                                         currentCohort%year_net_uptake(i) = (currentCohort%n*currentCohort%year_net_uptake(i) + &
                                               nextc%n*nextc%year_net_uptake(i))/newn                
                                      endif
                                   enddo
                                   
                                end if !(currentCohort%isnew)

                                currentCohort%n = newn     

                                ! Set pointers and remove the current cohort from the list
                                
                                shorterCohort => nextc%shorter
                                tallerCohort  => nextc%taller
                                
                                if (.not. associated(tallerCohort)) then
                                   currentPatch%tallest => shorterCohort
                                   if(associated(shorterCohort)) shorterCohort%taller => null()
                                else 
                                   tallerCohort%shorter => shorterCohort
                                endif
                                
                                if (.not. associated(shorterCohort)) then
                                   currentPatch%shortest => tallerCohort
                                   if(associated(tallerCohort)) tallerCohort%shorter => null()
                                else 
                                   shorterCohort%taller => tallerCohort
                                endif
                                
                                ! At this point, nothing should be pointing to current Cohort
				! update hydraulics quantities that are functions of hite & biomasses
				! deallocate the hydro structure of nextc
                                if (hlm_use_planthydro.eq.itrue) then				    
				    call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
				          currentCohort%pft,currentCohort%c_area)
                                    leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                                    currentCohort%treelai = tree_lai(leaf_c,             &
                                       currentCohort%pft, currentCohort%c_area, currentCohort%n, &
                                       currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, &
                                       currentCohort%vcmax25top  )			    
				   call updateSizeDepTreeHydProps(currentSite,currentCohort, bc_in)  				   
				   call DeallocateHydrCohort(nextc)
				endif

                                ! Deallocate the cohort's PRT structure
                                call nextc%prt%DeallocatePRTVartypes()
                                deallocate(nextc%prt)
                                deallocate(nextc)
                                nullify(nextc)

                             endif ! if( currentCohort%isnew.eqv.nextc%isnew ) then
                          endif !canopy layer
                       endif !pft
                    endif  !index no. 
                 endif !diff   
                 
                 nextc => nextnextc
                 
              enddo !end checking nextc cohort loop

              ! Ususally we always point to the next cohort. But remember ...
              ! this loop exits when current becomes the shortest, not when
              ! it finishes and becomes the null pointer.  If there is no
              ! shorter cohort, then it is shortest, and will exit
              ! Note also that it is possible that it entered here as the shortest
              ! which is possible if nextc was the shortest and was removed.

              if (associated (currentCohort%shorter)) then
                 currentCohort => currentCohort%shorter
              endif
              
           enddo !end currentCohort cohort loop

           !---------------------------------------------------------------------!
           ! Is the number of cohorts larger than the maximum?                   !
           !---------------------------------------------------------------------!   
           nocohorts = 0
           currentCohort => currentPatch%tallest
           do while(associated(currentCohort))
              nocohorts = nocohorts + 1
              currentCohort => currentCohort%shorter
           enddo

           if (nocohorts > maxCohortsPerPatch) then
              iterate = 1
              !---------------------------------------------------------------------!
              ! Making profile tolerance larger means that more fusion will happen  !
              !---------------------------------------------------------------------!        
              dynamic_fusion_tolerance = dynamic_fusion_tolerance * 1.1_r8

              !write(fates_log(),*) 'maxcohorts exceeded',dynamic_fusion_tolerance

           else

              iterate = 0
        endif

        if ( dynamic_fusion_tolerance .gt. 100._r8) then
              ! something has gone terribly wrong and we need to report what
              write(fates_log(),*) 'exceeded reasonable expectation of cohort fusion.'
              currentCohort => currentPatch%tallest
              nocohorts = 0
              do while(associated(currentCohort))
                 write(fates_log(),*) 'cohort ', nocohorts, currentCohort%dbh, currentCohort%canopy_layer, currentCohort%n
                 nocohorts = nocohorts + 1
                 currentCohort => currentCohort%shorter
              enddo
              call endrun(msg=errMsg(sourcefile, __LINE__))
           endif

        enddo !do while nocohorts>maxcohorts

     endif ! patch. 

     if (fusion_took_place == 1) then  ! if fusion(s) occured sort cohorts 
        call sort_cohorts(currentPatch)
     endif

  end subroutine fuse_cohorts

!-------------------------------------------------------------------------------------!

  subroutine sort_cohorts(patchptr)  
    ! ============================================================================
    !                 sort cohorts into the correct order   DO NOT CHANGE THIS IT WILL BREAK
    ! ============================================================================

    type(ed_patch_type) , intent(inout), target :: patchptr

    type(ed_patch_type) , pointer :: current_patch
    type(ed_cohort_type), pointer :: current_c, next_c
    type(ed_cohort_type), pointer :: shortestc, tallestc 
    type(ed_cohort_type), pointer :: storesmallcohort 
    type(ed_cohort_type), pointer :: storebigcohort   
    integer :: snull,tnull

    current_patch => patchptr
    tallestc  => NULL()
    shortestc => NULL()
    storebigcohort   => null()
    storesmallcohort => null()
    current_c => current_patch%tallest 

    do while (associated(current_c))  
       next_c => current_c%shorter
       tallestc  => storebigcohort 
       shortestc => storesmallcohort   
       if (associated(tallestc)) then
          tnull = 0
       else
          tnull = 1
          tallestc => current_c
       endif

       if (associated(shortestc)) then
          snull = 0
       else
          snull = 1
          shortestc => current_c
       endif

       call insert_cohort(current_c, tallestc, shortestc, tnull, snull, storebigcohort, storesmallcohort)

       current_patch%tallest  => storebigcohort 
       current_patch%shortest => storesmallcohort
       current_c => next_c

    enddo

  end subroutine sort_cohorts

  !-------------------------------------------------------------------------------------!
  subroutine insert_cohort(pcc, ptall, pshort, tnull, snull, storebigcohort, storesmallcohort)
    !
    ! !DESCRIPTION:
    ! Insert cohort into linked list                  
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_cohort_type) , intent(inout), target          :: pcc
    type(ed_cohort_type) , intent(inout), target          :: ptall
    type(ed_cohort_type) , intent(inout), target          :: pshort
    integer              , intent(in)                     :: tnull
    integer              , intent(in)                     :: snull
    type(ed_cohort_type) , intent(inout),pointer,optional :: storesmallcohort ! storage of the smallest cohort for insertion routine
    type(ed_cohort_type) , intent(inout),pointer,optional :: storebigcohort   ! storage of the largest cohort for insertion routine 
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: current
    type(ed_cohort_type), pointer :: tallptr, shortptr, icohort
    type(ed_cohort_type), pointer :: ptallest, pshortest 
    real(r8) :: tsp
    integer :: tallptrnull,exitloop
    !----------------------------------------------------------------------

    currentPatch => pcc%patchptr
    ptallest => ptall
    pshortest => pshort

    if (tnull == 1) then
       ptallest => null()
    endif
    if (snull == 1) then
       pshortest => null()
    endif

    icohort => pcc ! assign address to icohort local name  
    !place in the correct place in the linked list of heights 
    !begin by finding cohort that is just taller than the new cohort 
    tsp = icohort%hite

    current => pshortest
    exitloop = 0
    !starting with shortest tree on the grid, find tree just  
    !taller than tree being considered and return its pointer 
    if (associated(current)) then
       do while (associated(current).and.exitloop == 0)
          if (current%hite < tsp) then
             current => current%taller   
          else
             exitloop = 1 
          endif
       enddo
    endif

    if (associated(current)) then
       tallptr => current
       tallptrnull = 0
    else
       tallptr => null()
       tallptrnull = 1
    endif

    !new cohort is tallest 
    if (.not.associated(tallptr)) then  
       !new shorter cohort to the new cohort is the old tallest cohort 
       shortptr => ptallest

       !new cohort is tallest cohort and next taller remains null 
       ptallest => icohort
       if (present(storebigcohort)) then
          storebigcohort => icohort
       end if
       currentPatch%tallest => icohort 
       icohort%patchptr%tallest => icohort  
       !new cohort is not tallest 
    else
       !next shorter cohort to new cohort is the next shorter cohort 
       !to the cohort just taller than the new cohort 
       shortptr => tallptr%shorter

       !new cohort becomes the next shorter cohort to the cohort 
       !just taller than the new cohort 
       tallptr%shorter => icohort
    endif

    !new cohort is shortest 
    if (.not.associated(shortptr)) then
       !next shorter reamins null 
       !cohort is placed at the bottom of the list 
       pshortest => icohort
       if (present(storesmallcohort)) then
          storesmallcohort => icohort 
       end if
       currentPatch%shortest => icohort  
       icohort%patchptr%shortest => icohort 
    else
       !new cohort is not shortest and becomes next taller cohort 
       !to the cohort just below it as defined in the previous block 
       shortptr%taller => icohort
    endif

    ! assign taller and shorter links for the new cohort 
    icohort%taller => tallptr
    if (tallptrnull == 1) then 
       icohort%taller=> null()
    endif
    icohort%shorter => shortptr

  end subroutine insert_cohort

  !-------------------------------------------------------------------------------------!
  subroutine copy_cohort( currentCohort,copyc )
    !
    ! !DESCRIPTION:
    ! Copies all the variables in one cohort into another empty cohort                                    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_cohort_type), intent(inout) , target ::  copyc         ! New cohort argument.
    type(ed_cohort_type), intent(in)    , target ::  currentCohort ! Old cohort argument.
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::  n,o           ! New and old cohort pointers
    !----------------------------------------------------------------------

    o => currentCohort
    n => copyc

    n%indexnumber     = fates_unset_int
    
    ! VEGETATION STRUCTURE
    n%pft             = o%pft
    n%n               = o%n                         
    n%dbh             = o%dbh                                        
    n%hite            = o%hite
    n%laimemory       = o%laimemory
    n%lai             = o%lai                         
    n%sai             = o%sai  
    n%g_sb_laweight   = o%g_sb_laweight
    n%leaf_cost       = o%leaf_cost
    n%canopy_layer    = o%canopy_layer
    n%canopy_layer_yesterday    = o%canopy_layer_yesterday
    n%nv              = o%nv
    n%status_coh      = o%status_coh
    n%canopy_trim     = o%canopy_trim
    n%excl_weight     = o%excl_weight               
    n%prom_weight     = o%prom_weight               
    n%size_class      = o%size_class
    n%size_class_lasttimestep      = o%size_class_lasttimestep
    n%size_by_pft_class = o%size_by_pft_class

    ! This transfers the PRT objects over.
    call n%prt%CopyPRTVartypes(o%prt)

    ! Leaf biophysical rates
    n%vcmax25top = o%vcmax25top
    n%jmax25top  = o%jmax25top
    n%tpu25top   = o%tpu25top
    n%kp25top    = o%kp25top 

    ! CARBON FLUXES
    n%gpp_acc_hold    = o%gpp_acc_hold
    n%gpp_acc         = o%gpp_acc
    n%gpp_tstep       = o%gpp_tstep

    n%npp_acc_hold    = o%npp_acc_hold
    n%npp_tstep       = o%npp_tstep
    n%npp_acc         = o%npp_acc

    if ( debug ) write(fates_log(),*) 'EDcohortDyn Ia ',o%npp_acc
    if ( debug ) write(fates_log(),*) 'EDcohortDyn Ib ',o%resp_acc

    n%resp_tstep      = o%resp_tstep
    n%resp_acc        = o%resp_acc
    n%resp_acc_hold   = o%resp_acc_hold
    n%year_net_uptake = o%year_net_uptake
    n%ts_net_uptake   = o%ts_net_uptake

    ! C13 discrimination
    n%c13disc_clm   = o%c13disc_clm
    n%c13disc_acc   = o%c13disc_acc

    !RESPIRATION
    n%rdark           = o%rdark
    n%resp_m          = o%resp_m
    n%resp_g          = o%resp_g
    n%livestem_mr     = o%livestem_mr
    n%livecroot_mr    = o%livecroot_mr
    n%froot_mr        = o%froot_mr
 
    ! ALLOCATION
    n%dmort           = o%dmort
    n%seed_prod       = o%seed_prod
    n%treelai         = o%treelai
    n%treesai         = o%treesai
    n%c_area          = o%c_area

    ! Mortality diagnostics
    n%cmort = o%cmort
    n%bmort = o%bmort
    n%hmort = o%hmort
    n%frmort = o%frmort

    ! logging mortalities, Yi Xu
    n%lmort_direct     =o%lmort_direct
    n%lmort_collateral =o%lmort_collateral
    n%lmort_infra      =o%lmort_infra
    n%l_degrad         =o%l_degrad    

    ! Flags
    n%isnew = o%isnew

    ! VARIABLES NEEDED FOR INTEGRATION 
    n%dndt            = o%dndt
    n%dhdt            = o%dhdt
    n%ddbhdt          = o%ddbhdt

    if ( debug ) write(fates_log(),*) 'EDCohortDyn dpstoredt ',o%dbstoredt

    ! FIRE 
    n%fraction_crown_burned = o%fraction_crown_burned
    n%fire_mort             = o%fire_mort
    n%crownfire_mort        = o%crownfire_mort
    n%cambial_mort          = o%cambial_mort

    ! Plant Hydraulics
    
    if( hlm_use_planthydro.eq.itrue ) then
      call CopyCohortHydraulics(n,o)
    endif

    ! indices for binning
    n%size_class      = o%size_class
    n%size_class_lasttimestep      = o%size_class_lasttimestep
    n%size_by_pft_class   = o%size_by_pft_class

    !Pointers
    n%taller          => NULL()     ! pointer to next tallest cohort     
    n%shorter         => NULL()     ! pointer to next shorter cohort     
    n%patchptr        => o%patchptr ! pointer to patch that cohort is in 

  end subroutine copy_cohort

  !-------------------------------------------------------------------------------------!
  function count_cohorts( currentPatch ) result ( backcount )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), target :: currentPatch      !new site
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::currentCohort   !new patch
    integer backcount
    !----------------------------------------------------------------------

    currentCohort => currentPatch%shortest

    currentPatch%countcohorts = 0
    do while (associated(currentCohort)) 
       currentPatch%countcohorts = currentPatch%countcohorts + 1 
       currentCohort => currentCohort%taller  
    enddo

    backcount = 0
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort)) 
       backcount = backcount + 1
       currentCohort => currentCohort%shorter    
    enddo

    if (backcount /= currentPatch%countcohorts) then
       write(fates_log(),*) 'problem with linked list, not symmetrical' 
    endif

  end function count_cohorts

  ! ===================================================================================

  subroutine UpdateCohortBioPhysRates(currentCohort)

       ! --------------------------------------------------------------------------------
       ! This routine updates the four key biophysical rates of leaves
       ! based on the changes in a cohort's leaf age proportions
       !
       ! This should be called after growth.  Growth occurs
       ! after turnover and damage states are applied to the tree.
       ! Therefore, following growth, the leaf mass fractions
       ! of different age classes are unchanged until the next day.
       ! --------------------------------------------------------------------------------

       type(ed_cohort_type),intent(inout) :: currentCohort
       
       
       real(r8) :: frac_leaf_aclass(max_nleafage)  ! Fraction of leaves in each age-class
       integer  :: iage                            ! loop index for leaf ages
       integer  :: ipft                            ! plant functional type index

       ! First, calculate the fraction of leaves in each age class
       ! It is assumed that each class has the same proportion
       ! across leaf layers

       do iage = 1, nleafage
          frac_leaf_aclass(iage) = &
                currentCohort%prt%GetState(leaf_organ, all_carbon_elements,iage)
       end do

       ! If there are leaves, then perform proportional weighting on the four rates
       ! We assume that leaf age does not effect the specific leaf area, so the mass
       ! fractions are applicable to these rates
       
       if(sum(frac_leaf_aclass(1:nleafage))>nearzero) then

          ipft = currentCohort%pft

          frac_leaf_aclass(1:nleafage) =  frac_leaf_aclass(1:nleafage) / &
                sum(frac_leaf_aclass(1:nleafage))
          
          currentCohort%vcmax25top = sum(EDPftvarcon_inst%vcmax25top(ipft,1:nleafage) * &
                frac_leaf_aclass(1:nleafage))
          
          currentCohort%jmax25top  = sum(param_derived%jmax25top(ipft,1:nleafage) * &
                frac_leaf_aclass(1:nleafage))
          
          currentCohort%tpu25top   = sum(param_derived%tpu25top(ipft,1:nleafage) * &
                frac_leaf_aclass(1:nleafage))
          
          currentCohort%kp25top    = sum(param_derived%kp25top(ipft,1:nleafage) * & 
                frac_leaf_aclass(1:nleafage))

       else
          
          currentCohort%vcmax25top = 0._r8          
          currentCohort%jmax25top  = 0._r8
          currentCohort%tpu25top   = 0._r8
          currentCohort%kp25top    = 0._r8

       end if


       return
    end subroutine UpdateCohortBioPhysRates

  
  ! ============================================================================


  subroutine EvaluateAndCorrectDBH(currentCohort,delta_dbh,delta_hite)

    ! -----------------------------------------------------------------------------------
    ! If the current diameter of a plant is somehow less than what is allometrically 
    ! consistent with stuctural biomass (or, in the case of grasses, leaf biomass) 
    ! then correct (increase) the dbh to match that.
    ! -----------------------------------------------------------------------------------

    ! argument
    type(ed_cohort_type),intent(inout) :: currentCohort
    real(r8),intent(out)               :: delta_dbh
    real(r8),intent(out)               :: delta_hite
    
    ! locals
    real(r8) :: dbh
    real(r8) :: canopy_trim
    integer  :: ipft
    real(r8) :: sapw_area
    real(r8) :: target_sapw_c
    real(r8) :: target_agw_c
    real(r8) :: target_bgw_c
    real(r8) :: target_struct_c
    real(r8) :: target_leaf_c
    real(r8) :: struct_c
    real(r8) :: hite_out
    real(r8) :: leaf_c
    
    dbh  = currentCohort%dbh
    ipft = currentCohort%pft
    canopy_trim = currentCohort%canopy_trim

    delta_dbh   = 0._r8
    delta_hite  = 0._r8
    
    if( EDPftvarcon_inst%woody(ipft) == itrue) then

       struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
    
       ! Target sapwood biomass according to allometry and trimming [kgC]
       call bsap_allom(dbh,ipft,canopy_trim,sapw_area,target_sapw_c)
       
       ! Target total above ground biomass in woody/fibrous tissues  [kgC]
       call bagw_allom(dbh,ipft,target_agw_c)
       
       ! Target total below ground biomass in woody/fibrous tissues [kgC] 
       call bbgw_allom(dbh,ipft,target_bgw_c)
       
       ! Target total dead (structrual) biomass [kgC]
       call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)
       
       ! ------------------------------------------------------------------------------------
       ! If structure is larger than target, then we need to correct some integration errors
       ! by slightly increasing dbh to match it.
       ! For grasses, if leaf biomass is larger than target, then we reset dbh to match
       ! -----------------------------------------------------------------------------------
       
       if( (struct_c - target_struct_c ) > calloc_abs_error ) then
          call ForceDBH( ipft, canopy_trim, dbh, hite_out, bdead=struct_c )
          delta_dbh = dbh - currentCohort%dbh 
          delta_hite = hite_out - currentCohort%hite
          currentCohort%dbh  = dbh
          currentCohort%hite = hite_out
       end if
       
    else

       ! This returns the sum of leaf carbon over all (age) bins
       leaf_c  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)

       ! Target leaf biomass according to allometry and trimming
       call bleaf(dbh,ipft,canopy_trim,target_leaf_c)

       if( ( leaf_c - target_leaf_c ) > calloc_abs_error ) then
          call ForceDBH( ipft, canopy_trim, dbh, hite_out, bl=leaf_c )
          delta_dbh = dbh - currentCohort%dbh 
          delta_hite = hite_out - currentCohort%hite
          currentCohort%dbh = dbh
          currentCohort%hite = hite_out
       end if
       
    end if
    return
  end subroutine EvaluateAndCorrectDBH
  

end module EDCohortDynamicsMod
