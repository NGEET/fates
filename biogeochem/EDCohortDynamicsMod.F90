Module EDCohortDynamicsMod
  !
  ! DESCRIPTION:
  ! Cohort stuctures in FATES
  !
  
  ! USES:
  use FatesGlobals          , only : endrun => fates_endrun
  use FatesGlobals          , only : fates_log
  use FatesInterfaceTypesMod     , only : bc_in_type
  use FatesInterfaceTypesMod     , only : hlm_use_planthydro
  use FatesInterfaceTypesMod     , only : hlm_use_cohort_age_tracking
  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : fates_unset_r8
  use FatesConstantsMod     , only : nearzero
  use FatesConstantsMod     , only : calloc_abs_error
  use FatesInterfaceTypesMod     , only : nleafage
  use SFParamsMod           , only : SF_val_CWD_frac
  use EDPftvarcon           , only : EDPftvarcon_inst
  use EDPftvarcon           , only : GetDecompyFrac
  use PRTParametersMod      , only : prt_params
  use EDTypesMod            , only : ed_site_type
  use FatesPatchMod,          only : fates_patch_type
  use FatesCohortMod       , only : fates_cohort_type
  use EDParamsMod            , only : nclmax
  use PRTGenericMod         , only : element_list
  use PRTGenericMod         , only : StorageNutrientTarget
  use FatesLitterMod        , only : ncwd
  use FatesLitterMod        , only : ndcmpy
  use FatesLitterMod        , only : litter_type
  use FatesLitterMod        , only : adjust_SF_CWD_frac
  use EDParamsMod           , only : max_cohort_per_patch
  use EDTypesMod            , only : min_npm2, min_nppatch
  use EDTypesMod            , only : min_n_safemath
  use EDParamsMod            , only : nlevleaf
  use FatesConstantsMod     , only : ican_upper
  use EDTypesMod            , only : elem_diag_type
  use PRTGenericMod         , only : num_elements
  use FatesConstantsMod     , only : leaves_off
  use FatesConstantsMod     , only : leaves_shedding
  use FatesConstantsMod     , only : ihard_stress_decid
  use FatesConstantsMod     , only : isemi_stress_decid
  use EDParamsMod           , only : ED_val_cohort_age_fusion_tol
  use FatesInterfaceTypesMod      , only : hlm_use_planthydro
  use FatesInterfaceTypesMod      , only : hlm_parteh_mode
  use FatesPlantHydraulicsMod, only : FuseCohortHydraulics
  use FatesPlantHydraulicsMod, only : UpdateSizeDepPlantHydProps
  use FatesPlantHydraulicsMod, only : InitPlantHydStates
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
  use FatesPlantHydraulicsMod, only : UpdatePlantHydrNodes
  use FatesPlantHydraulicsMod, only : UpdatePlantHydrLenVol
  use FatesPlantHydraulicsMod, only : UpdatePlantKmax
  use FatesPlantHydraulicsMod, only : SavePreviousCompartmentVolumes
  use FatesPlantHydraulicsMod, only : ConstrainRecruitNumber
  use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
  use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
  use FatesAllometryMod  , only : bleaf
  use FatesAllometryMod  , only : bfineroot
  use FatesAllometryMod  , only : bsap_allom
  use FatesAllometryMod  , only : bagw_allom
  use FatesAllometryMod  , only : bbgw_allom
  use FatesAllometryMod  , only : bdead_allom
  use FatesAllometryMod  , only : h_allom
  use FatesAllometryMod  , only : carea_allom
  use FatesAllometryMod  , only : bstore_allom
  use FatesAllometryMod  , only : ForceDBH
  use FatesAllometryMod    , only : set_root_fraction
  use PRTGenericMod,          only : prt_carbon_allom_hyp
  use PRTGenericMod,          only : prt_cnp_flex_allom_hyp
  use PRTGenericMod,          only : prt_vartypes
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : nitrogen_element
  use PRTGenericMod,          only : phosphorus_element
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTAllometricCarbonMod, only : callom_prt_vartypes
  use PRTAllometricCNPMod,    only : cnp_allom_prt_vartypes
  use DamageMainMod,          only : undamaged_class
  use FatesConstantsMod,      only : i_term_mort_type_cstarv
  use FatesConstantsMod,      only : i_term_mort_type_canlev
  use FatesConstantsMod,      only : i_term_mort_type_numdens

  use shr_infnan_mod,         only : nan => shr_infnan_nan, assignment(=)  
  use shr_log_mod,            only : errMsg => shr_log_errMsg

  !
  implicit none
  private
  !
  public :: create_cohort
  public :: terminate_cohorts
  public :: terminate_cohort
  public :: fuse_cohorts
  public :: InitPRTObject
  public :: SendCohortToLitter
  public :: EvaluateAndCorrectDBH
  public :: DamageRecovery
  
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
  subroutine create_cohort(currentSite, patchptr, pft, nn, height, coage, dbh,     &
    prt, elongf_leaf, elongf_fnrt, elongf_stem, status, recruitstatus, ctrim,   &
    carea, clayer, crowndamage, spread, bc_in)

    !
    ! DESCRIPTION:
    ! create new cohort
    ! There are 4 places this is called
    ! 1) Initializing new cohorts at the beginning of a cold-start simulation
    ! 2) Initializing new recruits during dynamics
    ! 3) Initializing new cohorts at the beginning of a inventory read
    ! 4) Initializing new cohorts during restart
    !
    ! It is assumed that in the first 3, this is called with a reasonable amount of starter information.
    !

    ! ARGUMENTS:
    type(ed_site_type),     intent(inout), target  :: currentSite   ! site object
    type(fates_patch_type), intent(inout), pointer :: patchptr      ! pointer to patch object
    integer,                intent(in)             :: pft           ! cohort Plant Functional Type
    integer,                intent(in)             :: crowndamage   ! cohort damage class
    integer,                intent(in)             :: clayer        ! canopy status of cohort [1=canopy; 2=understorey]
    integer,                intent(in)             :: status        ! growth status of plant [1=leaves off; 2=leaves on]
    integer,                intent(in)             :: recruitstatus ! recruit status of plant [1 = recruitment , 0 = other]                
    real(r8),               intent(in)             :: nn            ! number of individuals in cohort [/m2]
    real(r8),               intent(in)             :: height        ! cohort height [m]
    real(r8),               intent(in)             :: coage         ! cohort age [m]
    real(r8),               intent(in)             :: dbh           ! cohort diameter at breast height [cm]
    real(r8),               intent(in)             :: elongf_leaf   ! leaf elongation factor [fraction] - 0: fully abscissed; 1: fully flushed
    real(r8),               intent(in)             :: elongf_fnrt   ! fine-root "elongation factor" [fraction]
    real(r8),               intent(in)             :: elongf_stem   ! stem "elongation factor" [fraction]
    class(prt_vartypes),    intent(inout), pointer :: prt           ! allocated PARTEH object
    real(r8),               intent(in)             :: ctrim         ! fraction of the maximum leaf biomass we are targeting
    real(r8),               intent(in)             :: spread        ! how spread crowns are in horizontal space
    real(r8),               intent(in)             :: carea         ! area of cohort - ONLY USED IN SP MODE [m2]
    type(bc_in_type),       intent(in)             :: bc_in         ! external boundary conditions

    ! LOCAL VARIABLES:
    type(fates_cohort_type), pointer :: newCohort        ! pointer to New Cohort structure.
    real(r8)                         :: rmean_temp       ! running mean temperature
    integer                          :: nlevrhiz         ! number of rhizosphere layers

    !----------------------------------------------------------------------

    ! create new cohort
    allocate(newCohort)
    call newCohort%Create(prt, pft, nn, height, coage, dbh, status, ctrim, carea,            &
      clayer, crowndamage, spread, patchptr%canopy_layer_tlai, elongf_leaf, elongf_fnrt,    &
      elongf_stem)
      
    ! Allocate running mean functions

    !  (Keeping as an example)
    !! allocate(newCohort%tveg_lpa)
    !! call newCohort%tveg_lpa%InitRMean(ema_lpa,init_value=patchptr%tveg_lpa%GetMean())

    if (hlm_use_planthydro .eq. itrue) then

      nlevrhiz = currentSite%si_hydr%nlevrhiz

      ! This allocates array spaces
      call InitHydrCohort(currentSite, newCohort)

      ! zero out the water balance error
      newCohort%co_hydr%errh2o = 0._r8

      ! This calculates node heights
      call UpdatePlantHydrNodes(newCohort, newCohort%pft, &
        newCohort%height,currentSite%si_hydr)

      ! This calculates volumes and lengths
      call UpdatePlantHydrLenVol(newCohort,currentSite%si_hydr)

      ! This updates the Kmax's of the plant's compartments
      call UpdatePlantKmax(newCohort%co_hydr,newCohort,currentSite%si_hydr)

      ! Since this is a newly initialized plant, we set the previous compartment-size
      ! equal to the ones we just calculated.
      call SavePreviousCompartmentVolumes(newCohort%co_hydr)

      ! This comes up with starter suctions and then water contents
      ! based on the soil values
      call InitPlantHydStates(currentSite,newCohort)

      if(recruitstatus==1)then

        newCohort%co_hydr%is_newly_recruited = .true.

          ! If plant hydraulics is active, we must constrain the
          ! number density of the new recruits based on the moisture
          ! available to be subsumed in the new plant tissues.
          ! So we go through the process of pre-initializing the hydraulic
          ! states in the temporary cohort, to calculate this new number density
          rmean_temp = patchptr%tveg24%GetMean()
          call ConstrainRecruitNumber(currentSite, newCohort, patchptr,           &
            bc_in, rmean_temp)

      endif

    endif

    call patchptr%InsertCohort(newCohort)

  end subroutine create_cohort

   ! ------------------------------------------------------------------------------------!

  subroutine InitPRTObject(prt)

    ! -----------------------------------------------------------------------------------
    !
    ! This routine allocates the PARTEH object that is associated with each cohort.
    ! The argument that is passed in is a pointer that is then associated with this
    ! newly allocated object.
    ! The object that is allocated is the specific extended class for the hypothesis
    ! of choice.
    ! Following this, the object and its internal mappings are initialized.
    ! This routine does NOT set any of the initial conditions, or boundary conditions
    ! such as the organ/element masses.  Those are handled after this call.
    !
    ! -----------------------------------------------------------------------------------

    ! Argument
    class(prt_vartypes), pointer :: prt

    ! Potential Extended types
    class(callom_prt_vartypes), pointer :: c_allom_prt
    class(cnp_allom_prt_vartypes), pointer :: cnp_allom_prt


    select case(hlm_parteh_mode)
    case (prt_carbon_allom_hyp)

        allocate(c_allom_prt)
        prt => c_allom_prt

    case (prt_cnp_flex_allom_hyp)

       allocate(cnp_allom_prt)
       prt => cnp_allom_prt

    case DEFAULT

        write(fates_log(),*) 'You specified an unknown PRT module'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))

    end select

    ! This is the call to allocate the data structures in the PRT object
    ! This call will be extended to each specific class.

    call prt%InitPRTVartype()


    return
  end subroutine InitPRTObject

  !-------------------------------------------------------------------------------------!

  subroutine terminate_cohorts( currentSite, currentPatch, level , call_index, bc_in)
    !
    ! !DESCRIPTION:
    ! terminates all cohorts when they get too small
    !
    ! !USES:
    
    !
    ! !ARGUMENTS
    type (ed_site_type) , intent(inout) :: currentSite
    type (fates_patch_type), intent(inout) :: currentPatch
    integer             , intent(in)    :: level
    integer                             :: call_index
    type(bc_in_type), intent(in)        :: bc_in

    ! Important point regarding termination levels.  Termination is typically
    ! called after fusion.  We do this so that we can re-capture the biomass that would
    ! otherwise be lost from termination.  The biomass of a fused plant remains in the
    ! live pool.  However, some plant number densities can be so low that they
    ! can cause numerical instabilities.  Thus, we call terminate_cohorts at level=1
    ! before fusion to get rid of these cohorts that are so incredibly sparse, and then
    ! terminate the remainder at level 2 for various other reasons.

    !
    ! !LOCAL VARIABLES:
    type (fates_cohort_type) , pointer :: currentCohort
    type (fates_cohort_type) , pointer :: shorterCohort
    type (fates_cohort_type) , pointer :: tallerCohort

    real(r8) :: leaf_c    ! leaf carbon [kg]
    real(r8) :: store_c   ! storage carbon [kg]
    real(r8) :: sapw_c    ! sapwood carbon [kg]
    real(r8) :: fnrt_c    ! fineroot carbon [kg]
    real(r8) :: repro_c   ! reproductive carbon [kg]
    real(r8) :: struct_c  ! structural carbon [kg]
    integer :: terminate  ! do we terminate (itrue) or not (ifalse)
    integer :: istat      ! return status code
    character(len=255) :: smsg
    integer :: termination_type
    !----------------------------------------------------------------------

    currentCohort => currentPatch%shortest
    do while (associated(currentCohort))

       terminate = ifalse
       termination_type = 0
       tallerCohort => currentCohort%taller

       leaf_c  = currentCohort%prt%GetState(leaf_organ, carbon12_element)
       store_c = currentCohort%prt%GetState(store_organ, carbon12_element)
       sapw_c  = currentCohort%prt%GetState(sapw_organ, carbon12_element)
       fnrt_c  = currentCohort%prt%GetState(fnrt_organ, carbon12_element)
       struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
       repro_c  = currentCohort%prt%GetState(repro_organ, carbon12_element)

       ! Check if number density is so low is breaks math (level 1)
       if (currentcohort%n <  min_n_safemath .and. level == 1) then
          terminate = itrue
          termination_type = i_term_mort_type_numdens
          if ( debug ) then
             write(fates_log(),*) 'terminating cohorts 0',currentCohort%n/currentPatch%area,currentCohort%dbh,currentCohort%pft,call_index
          endif
       endif

       ! The rest of these are only allowed if we are not dealing with a recruit (level 2)
       if (.not.currentCohort%isnew .and. level == 2) then

         ! Not enough n or dbh
         if  (currentCohort%n/currentPatch%area <= min_npm2 .or.	&  !
              currentCohort%n <= min_nppatch .or. &
              (currentCohort%dbh < 0.00001_r8 .and. store_c < 0._r8) ) then
            terminate = itrue
            termination_type = i_term_mort_type_numdens
            if ( debug ) then
               write(fates_log(),*) 'terminating cohorts 1',currentCohort%n/currentPatch%area,currentCohort%dbh,currentCohort%pft,call_index
            endif
         endif

         ! Outside the maximum canopy layer
         if (currentCohort%canopy_layer > nclmax ) then
           terminate = itrue
           termination_type = i_term_mort_type_canlev
           if ( debug ) then
             write(fates_log(),*) 'terminating cohorts 2', currentCohort%canopy_layer,currentCohort%pft,call_index
           endif
         endif

         ! live biomass pools are terminally depleted
         if ( ( sapw_c+leaf_c+fnrt_c ) < 1e-10_r8  .or.  &
               store_c  < 1e-10_r8) then
            terminate = itrue
            termination_type = i_term_mort_type_cstarv
            if ( debug ) then
              write(fates_log(),*) 'terminating cohorts 3', &
                    sapw_c,leaf_c,fnrt_c,store_c,currentCohort%pft,call_index
            endif
         endif

         ! Total cohort biomass is negative
         if ( ( struct_c+sapw_c+leaf_c+fnrt_c+store_c ) < 0._r8) then
            terminate = itrue
            termination_type = i_term_mort_type_cstarv
            if ( debug ) then
               write(fates_log(),*) 'terminating cohorts 4', &
                    struct_c,sapw_c,leaf_c,fnrt_c,store_c,currentCohort%pft,call_index
            endif

        endif
      endif    !  if (.not.currentCohort%isnew .and. level == 2) then

      if (terminate == itrue) then
         call terminate_cohort(currentSite, currentPatch, currentCohort, bc_in, termination_type)
         deallocate(currentCohort, stat=istat, errmsg=smsg)
         if (istat/=0) then
            write(fates_log(),*) 'dealloc001: fail on terminate_cohorts:deallocate(currentCohort):'//trim(smsg)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
      endif
      currentCohort => tallerCohort
   enddo

  end subroutine terminate_cohorts

  !-------------------------------------------------------------------------------------!
  subroutine terminate_cohort(currentSite, currentPatch, currentCohort, bc_in, termination_type)
   !
   ! !DESCRIPTION:
   ! Terminates an individual cohort and updates the site-level
   ! updates the carbon flux and nuber of individuals appropriately
   !
   ! !USES:
   !
   ! !ARGUMENTS
   type (ed_site_type)  , intent(inout), target :: currentSite
   type (fates_patch_type) , intent(inout), target :: currentPatch
   type (fates_cohort_type), intent(inout), target :: currentCohort
   type(bc_in_type), intent(in)                :: bc_in
   integer, intent(in)                         :: termination_type

   ! !LOCAL VARIABLES:
   type (fates_cohort_type) , pointer :: shorterCohort
   type (fates_cohort_type) , pointer :: tallerCohort

   real(r8) :: leaf_c    ! leaf carbon [kg]
   real(r8) :: store_c   ! storage carbon [kg]
   real(r8) :: sapw_c    ! sapwood carbon [kg]
   real(r8) :: fnrt_c    ! fineroot carbon [kg]
   real(r8) :: repro_c   ! reproductive carbon [kg]
   real(r8) :: struct_c  ! structural carbon [kg]
   integer :: terminate  ! do we terminate (itrue) or not (ifalse)
   integer :: c           ! counter for litter size class.
   integer :: levcan      ! canopy level
   
   !----------------------------------------------------------------------

   ! check termination_type; it should not be 0
   if (termination_type == 0) then
      write(fates_log(),*) 'termination_type=0'
      call endrun(msg=errMsg(sourcefile, __LINE__))
   endif
   
   leaf_c  = currentCohort%prt%GetState(leaf_organ, carbon12_element)
   store_c = currentCohort%prt%GetState(store_organ, carbon12_element)
   sapw_c  = currentCohort%prt%GetState(sapw_organ, carbon12_element)
   fnrt_c  = currentCohort%prt%GetState(fnrt_organ, carbon12_element)
   struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
   repro_c  = currentCohort%prt%GetState(repro_organ, carbon12_element)
   
   ! preserve a record of the to-be-terminated cohort for mortality accounting
   levcan = currentCohort%canopy_layer

   if( hlm_use_planthydro == itrue ) &
      call AccumulateMortalityWaterStorage(currentSite,currentCohort,currentCohort%n)

   ! Update the site-level carbon flux and individuals count for the appropriate canopy layer
   if(levcan==ican_upper) then
      currentSite%term_nindivs_canopy(termination_type,currentCohort%size_class,currentCohort%pft) = &
            currentSite%term_nindivs_canopy(termination_type,currentCohort%size_class,currentCohort%pft) + currentCohort%n

      currentSite%term_carbonflux_canopy(termination_type,currentCohort%pft) = currentSite%term_carbonflux_canopy(termination_type,currentCohort%pft) + &
            currentCohort%n * (struct_c+sapw_c+leaf_c+fnrt_c+store_c+repro_c)
   else
      currentSite%term_nindivs_ustory(termination_type,currentCohort%size_class,currentCohort%pft) = &
            currentSite%term_nindivs_ustory(termination_type,currentCohort%size_class,currentCohort%pft) + currentCohort%n

      currentSite%term_carbonflux_ustory(termination_type,currentCohort%pft) = currentSite%term_carbonflux_ustory(termination_type,currentCohort%pft) + &
            currentCohort%n * (struct_c+sapw_c+leaf_c+fnrt_c+store_c+repro_c)
   end if

   currentSite%term_abg_flux(currentCohort%size_class, currentCohort%pft) = &
        currentSite%term_abg_flux(currentCohort%size_class, currentCohort%pft) + &
        currentCohort%n * ( (struct_c+sapw_c+store_c) * prt_params%allom_agb_frac(currentCohort%pft) + &
        leaf_c )
  

   ! put the litter from the terminated cohorts
   ! straight into the fragmenting pools

   if (currentCohort%n.gt.0.0_r8) then
      call SendCohortToLitter(currentSite,currentPatch, &
           currentCohort,currentCohort%n,bc_in)
   end if

   ! Set pointers and deallocate the current cohort from the list
   shorterCohort => currentCohort%shorter
   tallerCohort => currentCohort%taller

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

   call currentCohort%FreeMemory()

 end subroutine terminate_cohort 
  
  ! =====================================================================================

  subroutine SendCohortToLitter(csite,cpatch,ccohort,nplant,bc_in)

    ! -----------------------------------------------------------------------------------
    ! This routine transfers the existing mass in all pools and all elements
    ! on a vegetation cohort, into the litter pool.
    !
    ! Important: (1) This IS NOT turnover, this is not a partial transfer.
    !            (2) This is from a select number of plants in the cohort. ie this is
    !                not a "whole-sale" sending of all plants to litter.
    !            (3) This does not affect the PER PLANT mass pools, so
    !                do not update any PARTEH structures.
    !            (4) The change in plant number density (due to death or termination)
    !                IS NOT handled here.
    !            (5) This routine is NOT used for disturbance, mostly
    !                because this routine assumes a cohort lands in its patch
    !                Whereas the disturbance scheme does NOT assume that.
    ! -----------------------------------------------------------------------------------

    ! Arguments
    type (ed_site_type)   , target  :: csite
    type (fates_patch_type)  , target  :: cpatch
    type (fates_cohort_type) , target  :: ccohort
    real(r8)                        :: nplant     ! Number (absolute)
                                                  ! of plants to transfer
    type(bc_in_type), intent(in)    :: bc_in

    type(litter_type), pointer        :: litt       ! Litter object for each element
    type(elem_diag_type),pointer :: elflux_diags

    real(r8) :: leaf_m    ! leaf mass [kg]
    real(r8) :: store_m   ! storage mass [kg]
    real(r8) :: sapw_m    ! sapwood mass [kg]
    real(r8) :: fnrt_m    ! fineroot mass [kg]
    real(r8) :: repro_m   ! reproductive mass [kg]
    real(r8) :: struct_m  ! structural mass [kg]
    real(r8) :: plant_dens! plant density [/m2]
    real(r8) :: dcmpy_frac! fraction of mass going to each decomposability partition
    integer  :: el        ! loop index for elements
    integer  :: c         ! loop index for CWD
    integer  :: pft       ! pft index of the cohort
    integer  :: crowndamage ! the crown damage class of the cohort
    integer  :: sl        ! loop index for soil layers
    integer  :: dcmpy     ! loop index for decomposability
    real(r8) :: SF_val_CWD_frac_adj(4) !Updated wood partitioning to CWD based on dbh
    !----------------------------------------------------------------------

    pft = ccohort%pft

    plant_dens = nplant/cpatch%area

    call set_root_fraction(csite%rootfrac_scr, pft, csite%zi_soil, &
         bc_in%max_rooting_depth_index_col)

    do el=1,num_elements

       store_m  = ccohort%prt%GetState(store_organ, element_list(el))
       fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
       repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
       if (prt_params%woody(ccohort%pft) == itrue) then
          leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
          sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
          struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
       else
          leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el)) + &
               ccohort%prt%GetState(sapw_organ, element_list(el)) + &
               ccohort%prt%GetState(struct_organ, element_list(el))
          sapw_m   = 0._r8
          struct_m = 0._r8
       endif

       litt => cpatch%litter(el)
       elflux_diags => csite%flux_diags%elem(el)

       !adjust how wood is partitioned between the cwd classes based on cohort dbh
       call adjust_SF_CWD_frac(ccohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)

       do c=1,ncwd

          ! above ground CWD
          litt%ag_cwd(c) = litt%ag_cwd(c) + plant_dens * &
               (struct_m+sapw_m)  * SF_val_CWD_frac_adj(c) * &
               prt_params%allom_agb_frac(pft)

          ! below ground CWD
          do sl=1,csite%nlevsoil
             litt%bg_cwd(c,sl) = litt%bg_cwd(c,sl) + plant_dens * &
                  (struct_m+sapw_m) * SF_val_CWD_frac_adj(c) * &
                  (1.0_r8 - prt_params%allom_agb_frac(pft)) * &
                  csite%rootfrac_scr(sl)
          enddo

          ! above ground
          elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
                (struct_m+sapw_m) * SF_val_CWD_frac_adj(c) * &
                prt_params%allom_agb_frac(pft) * nplant

          ! below ground
          elflux_diags%cwd_bg_input(c)  = elflux_diags%cwd_bg_input(c) + &
                (struct_m + sapw_m) * SF_val_CWD_frac_adj(c) * &
                (1.0_r8 - prt_params%allom_agb_frac(pft)) * nplant

       enddo

       do dcmpy=1,ndcmpy
           dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)

           litt%leaf_fines(dcmpy) = litt%leaf_fines(dcmpy) + &
                 plant_dens * (leaf_m+repro_m) * dcmpy_frac

           dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
           do sl=1,csite%nlevsoil
               litt%root_fines(dcmpy,sl) = litt%root_fines(dcmpy,sl) + &
                     plant_dens * (fnrt_m+store_m) * csite%rootfrac_scr(sl) * dcmpy_frac
           end do

       end do

       elflux_diags%surf_fine_litter_input(pft) = &
             elflux_diags%surf_fine_litter_input(pft) +  &
             (leaf_m+repro_m) * nplant
       elflux_diags%root_litter_input(pft) = &
             elflux_diags%root_litter_input(pft) +  &
             (fnrt_m+store_m) * nplant


    end do

    return
  end subroutine SendCohortToLitter

  !--------------------------------------------------------------------------------------

  subroutine fuse_cohorts(currentSite, currentPatch, bc_in)

     !
     ! !DESCRIPTION:
     ! Join similar cohorts to reduce total number
     !
     ! !USES:
     use EDParamsMod , only :  ED_val_cohort_size_fusion_tol
     use EDParamsMod , only :  ED_val_cohort_age_fusion_tol
     use FatesInterfaceTypesMod , only :  hlm_use_cohort_age_tracking
     use FatesConstantsMod , only : itrue
     use FatesConstantsMod, only : days_per_year
     

     !
     ! !ARGUMENTS
     type (ed_site_type), intent(inout)           :: currentSite
     type (fates_patch_type), intent(inout), pointer :: currentPatch
     type (bc_in_type), intent(in)                :: bc_in
     !

     ! !LOCAL VARIABLES:
     type (fates_cohort_type) , pointer :: currentCohort
     type (fates_cohort_type) , pointer :: nextc
     type (fates_cohort_type) , pointer :: nextnextc

     type (fates_cohort_type) , pointer :: shorterCohort
     type (fates_cohort_type) , pointer :: tallerCohort

     integer  :: i
     integer  :: fusion_took_place
     integer  :: iterate    ! do we need to keep fusing to get below maxcohorts?
     integer  :: nocohorts
     real(r8) :: newn
     real(r8) :: diff
     real(r8) :: coage_diff
     real(r8) :: leaf_c_next   ! Leaf carbon * plant density of current (for weighting)
     real(r8) :: leaf_c_curr   ! Leaf carbon * plant density of next (for weighting)
     real(r8) :: leaf_c_target
     real(r8) :: dynamic_size_fusion_tolerance
     real(r8) :: dynamic_age_fusion_tolerance
     real(r8) :: dbh
     real(r8) :: leaf_c             ! leaf carbon [kg]
     real(r8) :: target_storec      ! Target storage C
     
     integer  :: largersc, smallersc, sc_i        ! indices for tracking the growth flux caused by fusion
     real(r8) :: larger_n, smaller_n
     integer  :: oldercacls, youngercacls, cacls_i ! indices for tracking the age flux caused by fusion
     real(r8) :: older_n, younger_n
     real(r8) :: crown_reduction

     logical, parameter :: fuse_debug = .false.   ! This debug is over-verbose
                                                 ! and gets its own flag
     integer  :: istat         ! return status code
     character(len=255) :: smsg

     !----------------------------------------------------------------------

     !set initial fusion tolerance (in cm)
     dynamic_size_fusion_tolerance = ED_val_cohort_size_fusion_tol
     ! set the cohort age fusion tolerance (in fraction of years)
     dynamic_age_fusion_tolerance = ED_val_cohort_age_fusion_tol


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
                 diff = abs((currentCohort%dbh - nextc%dbh)/(0.5_r8*(currentCohort%dbh + nextc%dbh)))

                 !Criteria used to divide up the height continuum into different cohorts.

                 if (diff < dynamic_size_fusion_tolerance) then

                    ! Only fuse if the cohorts are within x years of each other
                    ! if they are the same age we make diff 0- to avoid errors divding by zero
                    !NB if cohort age tracking is off then the age of both should be 0
                    ! and hence the age fusion criterion is met
                    if (abs(currentCohort%coage - nextc%coage)<nearzero ) then
                       coage_diff = 0.0_r8
                    else
                       coage_diff = abs((currentCohort%coage - nextc%coage)/ &
                            (0.5_r8*(currentCohort%coage + nextc%coage)))
                    end if

                    if (coage_diff <= dynamic_age_fusion_tolerance ) then

                       ! Don't fuse a cohort with itself!
                       if (.not.associated(currentCohort,nextc) ) then

                          if (currentCohort%pft == nextc%pft) then

                             ! check cohorts have same damage class before fusing
                             if (currentCohort%crowndamage == nextc%crowndamage) then

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
                                      write(fates_log(),*) 'height:',currentCohort%height,nextc%height
                                      write(fates_log(),*) 'coage:',currentCohort%coage,nextc%coage
                                      write(fates_log(),*) 'dbh:',currentCohort%dbh,nextc%dbh
                                      write(fates_log(),*) 'pft:',currentCohort%pft,nextc%pft
                                      write(fates_log(),*) 'crowndamage:',currentCohort%crowndamage,nextc%crowndamage
                                      write(fates_log(),*) 'canopy_trim:',currentCohort%canopy_trim,nextc%canopy_trim
                                      write(fates_log(),*) 'canopy_layer_yesterday:', &
                                           currentCohort%canopy_layer_yesterday,nextc%canopy_layer_yesterday
                                      do i=1, nlevleaf
                                         write(fates_log(),*) 'leaf level: ',i,'year_net_uptake', &
                                              currentCohort%year_net_uptake(i),nextc%year_net_uptake(i)
                                      end do
                                   end if

                                   

                                      
                                   ! new cohort age is weighted mean of two cohorts
                                   currentCohort%coage = &
                                        (currentCohort%coage * (currentCohort%n/(currentCohort%n + nextc%n))) + &
                                        (nextc%coage * (nextc%n/(currentCohort%n + nextc%n)))

                                   ! update the cohort age again
                                   if (hlm_use_cohort_age_tracking .eq.itrue) then
                                      call coagetype_class_index(currentCohort%coage, currentCohort%pft, &
                                           currentCohort%coage_class, currentCohort%coage_by_pft_class)
                                   end if

                                   ! Fuse all mass pools
                                   call currentCohort%prt%WeightedFusePRTVartypes(nextc%prt, &
                                        currentCohort%n/newn )

                                   ! Leaf biophysical rates (use leaf mass weighting)
                                   ! -----------------------------------------------------------------
                                   call currentCohort%UpdateCohortBioPhysRates()
                                   
                                   currentCohort%l2fr = (currentCohort%n*currentCohort%l2fr &
                                        + nextc%n*nextc%l2fr)/newn

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
                                           currentCohort%crowndamage, &
                                           currentCohort%c_area,inverse=.false.)

                                      call carea_allom(nextc%dbh,nextc%n, &
                                           currentSite%spread,nextc%pft,&
                                           nextc%crowndamage, &
                                           nextc%c_area,inverse=.false.)

                                      currentCohort%c_area = currentCohort%c_area + nextc%c_area

                                      !
                                      dbh = currentCohort%dbh
                                      call carea_allom(dbh,newn,currentSite%spread,currentCohort%pft,&
                                           currentCohort%crowndamage,currentCohort%c_area,inverse=.true.)
                                      !
                                      if (abs(dbh-fates_unset_r8)<nearzero) then
                                         currentCohort%dbh = (currentCohort%n*currentCohort%dbh         &
                                              + nextc%n*nextc%dbh)/newn

                                         if( prt_params%woody(currentCohort%pft) == itrue ) then

                                            call ForceDBH( currentCohort%pft, currentCohort%crowndamage, & 
                                                 currentCohort%canopy_trim, &
                                                 currentCohort%efleaf_coh, currentCohort%efstem_coh, &
                                                 currentCohort%dbh, currentCohort%height, &
                                                 bdead = currentCohort%prt%GetState(struct_organ,carbon12_element))

                                         end if
                                         !
                                         call carea_allom(currentCohort%dbh,newn,currentSite%spread,currentCohort%pft,&
                                              currentCohort%crowndamage, currentCohort%c_area,inverse=.false.)

                                      else
                                         currentCohort%dbh = dbh
                                      endif

                                      !
                                      call h_allom(currentCohort%dbh,currentCohort%pft,currentCohort%height)
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
                                      call h_allom(currentCohort%dbh,currentCohort%pft,currentCohort%height)
                                      !
                                      ! -----------------------------------------------------------------
                                      ! If fusion pushed structural biomass to be larger than
                                      ! the allometric target value derived by diameter, we
                                      ! then increase diameter and height until the allometric
                                      ! target matches actual bdead. (if it is the other way around
                                      ! we then just let the carbon pools grow to fill out allometry)
                                      ! -----------------------------------------------------------------
                                      !
                                      if( prt_params%woody(currentCohort%pft) == itrue ) then
                                         call ForceDBH( currentCohort%pft, currentCohort%crowndamage, & 
                                              currentCohort%canopy_trim, &
                                              currentCohort%efleaf_coh, currentCohort%efstem_coh, &
                                              currentCohort%dbh, currentCohort%height, &
                                              bdead = currentCohort%prt%GetState(struct_organ,carbon12_element))

                                      end if
                                      !
                                      call carea_allom(currentCohort%dbh,newn,currentSite%spread,currentCohort%pft,&
                                           currentCohort%crowndamage, currentCohort%c_area,inverse=.false.)
                                      !
                                   case default
                                      write(fates_log(),*) 'FATES: Invalid choice for cohort_fusion_conservation_method'
                                      call endrun(msg=errMsg(sourcefile, __LINE__))
                                   end select

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
                                   ! their initization values, which should be the same for each

                                   if ( .not.currentCohort%isnew) then
                                      currentCohort%seed_prod      = (currentCohort%n*currentCohort%seed_prod + &
                                           nextc%n*nextc%seed_prod)/newn
                                      currentCohort%gpp_acc        = (currentCohort%n*currentCohort%gpp_acc     + &
                                           nextc%n*nextc%gpp_acc)/newn
                                      currentCohort%npp_acc        = (currentCohort%n*currentCohort%npp_acc     + &
                                           nextc%n*nextc%npp_acc)/newn
                                      currentCohort%resp_m_acc       = (currentCohort%n*currentCohort%resp_m_acc    + &
                                           nextc%n*nextc%resp_m_acc)/newn
                                      currentCohort%resp_m_acc_hold  = &
                                           (currentCohort%n*currentCohort%resp_m_acc_hold + &
                                           nextc%n*nextc%resp_m_acc_hold)/newn
                                      currentCohort%resp_g_acc_hold  = &
                                           (currentCohort%n*currentCohort%resp_g_acc_hold + &
                                           nextc%n*nextc%resp_g_acc_hold)/newn
                                      currentCohort%npp_acc_hold   = &
                                           (currentCohort%n*currentCohort%npp_acc_hold + &
                                           nextc%n*nextc%npp_acc_hold)/newn
                                      currentCohort%gpp_acc_hold   = &
                                           (currentCohort%n*currentCohort%gpp_acc_hold + &
                                           nextc%n*nextc%gpp_acc_hold)/newn

                                      currentCohort%resp_excess_hold = &
                                           (currentCohort%n*currentCohort%resp_excess_hold + &
                                           nextc%n*nextc%resp_excess_hold)/newn

                                      currentCohort%dmort          = (currentCohort%n*currentCohort%dmort       + &
                                           nextc%n*nextc%dmort)/newn

                                      currentCohort%fire_mort      = (currentCohort%n*currentCohort%fire_mort   + &
                                           nextc%n*nextc%fire_mort)/newn

                                      ! mortality diagnostics
                                      currentCohort%cmort = (currentCohort%n*currentCohort%cmort + nextc%n*nextc%cmort)/newn
                                      currentCohort%hmort = (currentCohort%n*currentCohort%hmort + nextc%n*nextc%hmort)/newn
                                      currentCohort%bmort = (currentCohort%n*currentCohort%bmort + nextc%n*nextc%bmort)/newn
                                      currentCohort%smort = (currentCohort%n*currentCohort%smort + nextc%n*nextc%smort)/newn
                                      currentCohort%asmort = (currentCohort%n*currentCohort%asmort + nextc%n*nextc%asmort)/newn
                                      currentCohort%frmort = (currentCohort%n*currentCohort%frmort + nextc%n*nextc%frmort)/newn

                                      ! Nutrients
                                      if(hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then

                                         if(nextc%n > currentCohort%n) currentCohort%cnp_limiter = nextc%cnp_limiter

                                         currentCohort%cx_int = (currentCohort%n*currentCohort%cx_int + &
                                              nextc%n*nextc%cx_int)/newn
                                         currentCohort%ema_dcxdt = (currentCohort%n*currentCohort%ema_dcxdt + &
                                              nextc%n*nextc%ema_dcxdt)/newn
                                         currentCohort%cx0 = (currentCohort%n*currentCohort%cx0 + &
                                              nextc%n*nextc%cx0)/newn

                                         ! These variables do not need to be rescaled because they
                                         ! are written to history immediately after calculation
                                         
                                         currentCohort%daily_nh4_uptake = (currentCohort%n*currentCohort%daily_nh4_uptake + &
                                              nextc%n*nextc%daily_nh4_uptake)/newn
                                         currentCohort%daily_no3_uptake = (currentCohort%n*currentCohort%daily_no3_uptake + &
                                              nextc%n*nextc%daily_no3_uptake)/newn
                                         currentCohort%sym_nfix_daily = (currentCohort%n*currentCohort%sym_nfix_daily + &
                                              nextc%n*nextc%sym_nfix_daily)/newn
                                         currentCohort%daily_n_gain = (currentCohort%n*currentCohort%daily_n_gain + &
                                              nextc%n*nextc%daily_n_gain)/newn
                                         currentCohort%daily_p_gain = (currentCohort%n*currentCohort%daily_p_gain + &
                                              nextc%n*nextc%daily_p_gain)/newn
                                         currentCohort%daily_p_demand = (currentCohort%n*currentCohort%daily_p_demand + &
                                              nextc%n*nextc%daily_p_demand)/newn
                                         currentCohort%daily_n_demand = (currentCohort%n*currentCohort%daily_n_demand + &
                                              nextc%n*nextc%daily_n_demand)/newn
                                         currentCohort%daily_c_efflux = (currentCohort%n*currentCohort%daily_c_efflux + &
                                              nextc%n*nextc%daily_c_efflux)/newn
                                         currentCohort%daily_n_efflux = (currentCohort%n*currentCohort%daily_n_efflux + &
                                              nextc%n*nextc%daily_n_efflux)/newn
                                         currentCohort%daily_p_efflux = (currentCohort%n*currentCohort%daily_p_efflux + &
                                              nextc%n*nextc%daily_p_efflux)/newn
                                      end if
                                         

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
                                   ! update hydraulics quantities that are functions of height & biomasses
                                   ! deallocate the hydro structure of nextc
                                   if (hlm_use_planthydro.eq.itrue) then
                                      call UpdateSizeDepPlantHydProps(currentSite,currentCohort, bc_in)
                                   endif

                                   call nextc%FreeMemory()
                                   deallocate(nextc, stat=istat, errmsg=smsg)
                                   if (istat/=0) then
                                      write(fates_log(),*) 'dealloc003: fail on deallocate(nextc):'//trim(smsg)
                                      call endrun(msg=errMsg(sourcefile, __LINE__))
                                   endif

                                endif ! if( currentCohort%isnew.eqv.nextc%isnew ) then
                             endif !canopy layer
                             endif ! crowndamage 
                          endif !pft
                       endif  !index no.
                    endif  ! cohort age diff
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


           if ( hlm_use_cohort_age_tracking .eq.itrue) then
              if ( nocohorts > max_cohort_per_patch ) then
                 iterate = 1
                 !---------------------------------------------------------------------!
                 ! Making profile tolerance larger means that more fusion will happen  !
                 !---------------------------------------------------------------------!
                 dynamic_size_fusion_tolerance = dynamic_size_fusion_tolerance * 1.1_r8
                 dynamic_age_fusion_tolerance = dynamic_age_fusion_tolerance * 1.1_r8

              else

                 iterate = 0
              endif

           else

              if (nocohorts > max_cohort_per_patch) then
                 iterate = 1
                 !---------------------------------------------------------------------!
                 ! Making profile tolerance larger means that more fusion will happen  !
                 !---------------------------------------------------------------------!
                 dynamic_size_fusion_tolerance = dynamic_size_fusion_tolerance * 1.1_r8

              else

                 iterate = 0
              endif
           end if


        if ( dynamic_size_fusion_tolerance .gt. 100._r8) then
              ! something has gone terribly wrong and we need to report what
              write(fates_log(),*) 'exceeded reasonable expectation of cohort fusion.'
              currentCohort => currentPatch%tallest
              nocohorts = 0
              do while(associated(currentCohort))
                 write(fates_log(),*) 'cohort ', nocohorts, currentCohort%dbh,&
                      currentCohort%coage, currentCohort%canopy_layer, currentCohort%n
                 nocohorts = nocohorts + 1
                 currentCohort => currentCohort%shorter
              enddo
              call endrun(msg=errMsg(sourcefile, __LINE__))
           endif

        enddo !do while nocohorts>maxcohorts

     endif ! patch.
     
     if (fusion_took_place == 1) then  ! if fusion(s) occured sort cohorts
        call currentPatch%SortCohorts()
        call currentPatch%ValidateCohorts()
     endif
   
  end subroutine fuse_cohorts

!-------------------------------------------------------------------------------------!

  subroutine EvaluateAndCorrectDBH(currentCohort,delta_dbh,delta_height)

    ! -----------------------------------------------------------------------------------
    ! If the current diameter of a plant is somehow less than what is allometrically
    ! consistent with stuctural biomass (or, in the case of grasses, leaf biomass)
    ! then correct (increase) the dbh to match that.
    ! -----------------------------------------------------------------------------------
    
    ! argument
    type(fates_cohort_type),intent(inout) :: currentCohort
    real(r8),intent(out)               :: delta_dbh
    real(r8),intent(out)               :: delta_height

    ! locals
    real(r8) :: dbh
    real(r8) :: canopy_trim
    integer  :: ipft
    integer  :: icrowndamage
    real(r8) :: sapw_area
    real(r8) :: target_sapw_c
    real(r8) :: target_agw_c
    real(r8) :: target_bgw_c
    real(r8) :: target_struct_c
    real(r8) :: target_leaf_c
    real(r8) :: struct_c
    real(r8) :: height_out
    real(r8) :: leaf_c
    real(r8) :: crown_reduction
    real(r8) :: elongf_leaf
    real(r8) :: elongf_stem
    
    dbh  = currentCohort%dbh
    ipft = currentCohort%pft
    icrowndamage = currentCohort%crowndamage
    canopy_trim = currentCohort%canopy_trim
    elongf_leaf = currentCohort%efleaf_coh
    elongf_stem = currentCohort%efstem_coh

    delta_dbh   = 0._r8
    delta_height  = 0._r8

    if( prt_params%woody(currentCohort%pft) == itrue) then

       struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)

       ! Target sapwood biomass according to allometry, trimming and phenology [kgC]
       call bsap_allom(dbh,ipft,icrowndamage,canopy_trim, elongf_stem, sapw_area,target_sapw_c)
       
       ! Target total above ground biomass in woody/fibrous tissues
       ! according to allometry, trimming and phenology [kgC]
       call bagw_allom(dbh,ipft, icrowndamage, elongf_stem, target_agw_c)
       
       ! Target total below ground biomass in woody/fibrous tissues
       ! according to allometry, trimming and phenology [kgC]
       call bbgw_allom(dbh,ipft, elongf_stem, target_bgw_c)

       ! Target total dead (structrual) biomass [kgC]
       call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)

       ! ------------------------------------------------------------------------------------
       ! If structure is larger than target, then we need to correct some integration errors
       ! by slightly increasing dbh to match it.
       ! For grasses, if leaf biomass is larger than target, then we reset dbh to match
       ! -----------------------------------------------------------------------------------

       if( (struct_c - target_struct_c ) > calloc_abs_error ) then

          call ForceDBH( ipft,icrowndamage,canopy_trim, elongf_leaf, elongf_stem, &
               dbh, height_out, bdead=struct_c)

          delta_dbh = dbh - currentCohort%dbh 
          delta_height = height_out - currentCohort%height
          currentCohort%dbh  = dbh
          currentCohort%height = height_out
       end if

    else

       ! This returns the sum of leaf carbon over all (age) bins
       leaf_c  = currentCohort%prt%GetState(leaf_organ, carbon12_element)

       ! Target leaf biomass according to allometry, trimming and phenology
       call bleaf(dbh,ipft,icrowndamage, canopy_trim, elongf_leaf, target_leaf_c)

       if( ( leaf_c - target_leaf_c ) > calloc_abs_error ) then
          call ForceDBH( ipft, icrowndamage, canopy_trim, elongf_leaf, elongf_stem, &
               dbh, height_out, bl=leaf_c )
          delta_dbh = dbh - currentCohort%dbh
          delta_height = height_out - currentCohort%height
          currentCohort%dbh = dbh
          currentCohort%height = height_out
       end if

    end if
    return
  end subroutine EvaluateAndCorrectDBH

  !------------------------------------------------------------------------------------

  subroutine DamageRecovery(csite,cpatch,ccohort,newly_recovered)

    !---------------------------------------------------------------------------
    ! JN March 2021
    ! At this point it is possible that damaged cohorts have reached their
    ! target allometries. There is a choice now - if they have excess carbon,
    ! they can use it to grow along their reduced allometric targets  - i.e.
    ! dbh and all carbon pools grow out together. OR they can use excess carbon to
    ! jump to a lower damage class by changing their target allometry and growing 
    ! to meet new C pools for same dbh.
    !
    ! d = damage class
    ! --------------------------------------------------------------------------

    type(ed_site_type)   :: csite            ! Site of the current cohort
    type(fates_patch_type)  :: cpatch           ! patch of the current cohort
    type(fates_cohort_type),pointer :: ccohort  ! Current (damaged) cohort
    logical              :: newly_recovered  ! true if we create a new cohort

    ! locals
    type(fates_cohort_type), pointer :: rcohort ! New cohort that recovers by
                                             ! having a lower damage class
    real(r8) :: sapw_area                    ! sapwood area
    real(r8) :: target_sapw_c,target_sapw_m  ! sapwood mass, C and N/P
    real(r8) :: target_agw_c                 ! target above ground wood
    real(r8) :: target_bgw_c                    ! target below ground wood
    real(r8) :: target_struct_c,target_struct_m ! target structural C and N/P
    real(r8) :: target_fnrt_c,target_fnrt_m     ! target fine-root C and N/P
    real(r8) :: target_leaf_c,target_leaf_m     ! target leaf C and N/P
    real(r8) :: target_store_c,target_store_m   ! target storage C and N/P
    real(r8) :: target_repro_m                  ! target reproductive C/N/P
    real(r8) :: leaf_m,fnrt_m,sapw_m            ! actual masses in organs C/N/P
    real(r8) :: struct_m,store_m,repro_m        ! actual masses in organs C/N/P
    real(r8) :: mass_d                          ! intermediate term for nplant_recover
    real(r8) :: mass_dminus1                    ! intermediate term for nplant_recover
    real(r8) :: available_m                     ! available mass that can be used to 
                                                ! improve damage class
    real(r8) :: recovery_demand                 ! amount of mass needed to get to 
                                                ! the target of the next damage class
    real(r8) :: max_recover_nplant              ! max number of plants that could get to
                                                ! target of next class
    real(r8) :: nplant_recover                  ! number of plants in cohort that will
                                                ! recover to the next class
    integer  :: el                                ! element loop counter

    logical  :: is_hydecid_dormant    ! Flag to signal that the cohort is drought deciduous and dormant
    logical  :: is_sedecid_dormant    ! Flag to signal this is a deciduous PFT
    
    associate(dbh => ccohort%dbh, &
         ipft => ccohort%pft, &
         canopy_trim => ccohort%canopy_trim, &
         elongf_leaf => ccohort%efleaf_coh,  &
         elongf_fnrt => ccohort%effnrt_coh,  &
         elongf_stem => ccohort%efstem_coh)

      ! If we are currently undamaged, no recovery
      ! necessary, do nothing and return a null pointer
      ! If the damage_recovery_scalar is zero, which
      ! would be an unusual testing case, but possible,
      ! then no recovery is possible, do nothing and
      ! return a null pointer
      if ((ccohort%crowndamage == undamaged_class) .or. &
           (EDPftvarcon_inst%damage_recovery_scalar(ipft) < nearzero) ) then
         newly_recovered = .false.
         return
      end if


      !--- Set some logical flags to simplify "if" blocks
      is_hydecid_dormant = &
         any(prt_params%stress_decid(ipft) == [ihard_stress_decid,isemi_stress_decid] ) &
         .and. any(ccohort%status_coh == [leaves_off,leaves_shedding] )
      is_sedecid_dormant = &
         ( prt_params%season_decid(ipft) == itrue ) &
         .and. any(ccohort%status_coh == [leaves_off,leaves_shedding] )

      ! If plants are drought deciduous and are losing or lost all leaves, they cannot
      ! allocate carbon to any growth or recovery. Return a null pointer and wait until
      ! the growing season.
      if (is_hydecid_dormant) then
         newly_recovered = .false.
         return
      end if


      ! If we have not returned, then this cohort both has
      ! a damaged status, and the ability to recover from that damage
      ! -----------------------------------------------------------------

      ! To determine recovery, the first priority is to determine how much
      ! resources (C,N,P) are required to recover the plant to the target
      ! pool sizes of the next (less) damage class

      ! Target sapwood biomass according to allometry, trimming and phenology [kgC]
      call bsap_allom(dbh,ipft, ccohort%crowndamage-1, canopy_trim, elongf_stem, &
           sapw_area,target_sapw_c)
      ! Target total above ground biomass in woody/fibrous tissues
      ! according to allometry, trimming and phenology [kgC]
      call bagw_allom(dbh,ipft, ccohort%crowndamage-1, elongf_stem, target_agw_c)
      ! Target total below ground biomass in woody/fibrous tissues
      ! according to allometry, trimming and phenology [kgC]
      call bbgw_allom(dbh,ipft, elongf_stem, target_bgw_c)
      ! Target total dead (structrual) biomass [kgC]
      call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, ipft, target_struct_c)
      ! Target fine-root biomass according to allometry, trimming and phenology [kgC]
      call bfineroot(dbh,ipft,canopy_trim,ccohort%l2fr, elongf_fnrt, target_fnrt_c)
      ! Target storage carbon [kgC]
      call bstore_allom(dbh,ipft,ccohort%crowndamage-1, canopy_trim,target_store_c)
      ! Target leaf biomass according to allometry, trimming and phenology [kgC]
      call bleaf(dbh,ipft,ccohort%crowndamage-1, canopy_trim, elongf_leaf, target_leaf_c)


      ! If plants are cold deciduous, we do not let them recover leaves, but we allow
      ! them to recover other tissues. This is due to back-compatibility, but we may
      ! want to revisit this later.
      if (is_sedecid_dormant) then
         target_leaf_c   = 0._r8
      end if

      ! We will be taking the number of recovering plants
      ! based on minimum of available resources for C/N/P (initialize high)
      nplant_recover = 1.e10_r8
      
      do el=1,num_elements
         
         ! Actual mass of chemical species in the organs
         leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
         store_m  = ccohort%prt%GetState(store_organ, element_list(el))
         sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
         fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
         struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
         repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
         
         ! Target mass of chemical species in organs, based on stature,
         ! allometry and stoichiometry parameters
         select case (element_list(el))
         case (carbon12_element)
            target_store_m  = target_store_c
            target_leaf_m   = target_leaf_c
            target_fnrt_m   = target_fnrt_c
            target_struct_m = target_struct_c
            target_sapw_m   = target_sapw_c
            target_repro_m  = 0._r8
            available_m     = ccohort%npp_acc
         case (nitrogen_element) 
            target_struct_m = target_struct_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
            target_leaf_m = target_leaf_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
            target_fnrt_m = target_fnrt_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
            target_sapw_m = target_sapw_c * &
                 prt_params%nitr_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
            target_repro_m  = 0._r8
            target_store_m = StorageNutrientTarget(ipft, element_list(el), &
                 target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            ! For nutrients, all uptake is immediately put into storage, so just swap
            ! them and assume storage is what is available, but needs to be filled up
            available_m     = store_m
            store_m         = 0._r8
         case (phosphorus_element)
            target_struct_m = target_struct_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(struct_organ))
            target_leaf_m = target_leaf_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(leaf_organ))
            target_fnrt_m = target_fnrt_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(fnrt_organ))
            target_sapw_m = target_sapw_c * &
                 prt_params%phos_stoich_p1(ipft,prt_params%organ_param_id(sapw_organ))
            target_repro_m  = 0._r8
            target_store_m = StorageNutrientTarget(ipft, element_list(el), &
                 target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            ! For nutrients, all uptake is immediately put into storage, so just swap
            ! them and assume storage is what is available, but needs to be filled up
            available_m     = store_m
            store_m         = 0._r8
         end select
         
         ! 1. What is excess carbon?
         ! carbon_balance
         
         !  2. What is biomass required to go from current
         !     damage level to next damage level?
         
         ! mass of this damage class
         mass_d = leaf_m + store_m + sapw_m + fnrt_m + struct_m + repro_m
         
         mass_dminus1 = max(leaf_m, target_leaf_m) + max(fnrt_m, target_fnrt_m) + &
              max(store_m, target_store_m) + max(sapw_m, target_sapw_m) + &
              max(struct_m, target_struct_m)
         
         ! Mass needed to get from current mass to allometric
         ! target mass of next damage class up
         recovery_demand = mass_dminus1 - mass_d
         
         ! 3. How many trees can get there with excess carbon?
         max_recover_nplant =  available_m * ccohort%n / recovery_demand 
         
         ! 4. Use the scalar to decide how many to recover
         nplant_recover = min(nplant_recover,min(ccohort%n,max(0._r8,max_recover_nplant * &
                              EDPftvarcon_inst%damage_recovery_scalar(ipft) )))
         
      end do
          
      if(nplant_recover < nearzero) then

         newly_recovered = .false.
         return
         
      else
         newly_recovered = .true.
         allocate(rcohort)
         if(hlm_use_planthydro .eq. itrue) call InitHydrCohort(csite,rcohort)
         ! Initialize the PARTEH object and point to the
         ! correct boundary condition fields
         rcohort%prt => null()
         call InitPRTObject(rcohort%prt)
         call rcohort%InitPRTBoundaryConditions()
         call ccohort%Copy(rcohort)

         rcohort%n = nplant_recover
          
         rcohort%crowndamage = ccohort%crowndamage - 1
         
         ! Need to adjust the crown area which is NOT on a per individual basis
         call carea_allom(dbh,rcohort%n,csite%spread,ipft,rcohort%crowndamage,rcohort%c_area)
        
         ! Update properties of the un-recovered (donor) cohort
         ccohort%n = ccohort%n - rcohort%n
         ccohort%c_area = ccohort%c_area * ccohort%n / (ccohort%n+rcohort%n)

         !----------- Insert copy into linked list ----------------------!
         ! This subroutine is called within a loop in EDMain that
         ! proceeds short to tall. We want the newly created cohort
         ! to have an opportunity to experience the list, so we add
         ! it in the list in a position taller than the current cohort
         ! --------------------------------------------------------------!
         
         rcohort%shorter => ccohort
         if(associated(ccohort%taller))then
            rcohort%taller => ccohort%taller
            ccohort%taller%shorter => rcohort
         else
            cpatch%tallest => rcohort    
            rcohort%taller => null()
         endif
         ccohort%taller => rcohort
         
      end if ! end if greater than nearzero

    end associate
    
    return
  end subroutine DamageRecovery




!:.........................................................................:

  

end module EDCohortDynamicsMod
