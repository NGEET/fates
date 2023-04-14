module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   use EDTypesMod                , only : ed_site_type
   use EDTypesMod                , only : dinc_vai
   use EDTypesMod                , only : dlower_vai
   use EDParamsMod               , only : ED_val_vai_top_bin_width
   use EDParamsMod               , only : ED_val_vai_width_increase_factor
   use EDParamsMod               , only : ED_val_history_damage_bin_edges
   use EDParamsMod               , only : maxpatch_total
   use EDParamsMod               , only : maxpatch_primary
   use EDParamsMod               , only : maxpatch_secondary
   use EDParamsMod               , only : max_cohort_per_patch
   use EDTypesMod                , only : maxSWb
   use EDTypesMod                , only : ivis
   use EDTypesMod                , only : inir
   use EDTypesMod                , only : nclmax
   use EDTypesMod                , only : nlevleaf
   use EDTypesMod                , only : maxpft
   use EDTypesMod                , only : do_fates_salinity
   use EDTypesMod                , only : numWaterMem
   use EDTypesMod                , only : numlevsoil_max
   use EDTypesMod                , only : ed_site_type
   use EDTypesMod                , only : ed_patch_type
   use EDTypesMod                , only : ed_cohort_type
   use EDTypesMod                , only : area_inv
   use EDTypesMod                , only : num_vegtemp_mem
   use FatesConstantsMod         , only : r8 => fates_r8
   use FatesConstantsMod         , only : itrue,ifalse
   use FatesConstantsMod         , only : nearzero
   use FatesConstantsMod         , only : sec_per_day
   use FatesConstantsMod         , only : days_per_year
   use FatesGlobals              , only : fates_global_verbose
   use FatesGlobals              , only : fates_log
   use FatesGlobals              , only : endrun => fates_endrun
   use FatesConstantsMod             , only : fates_unset_r8
   use FatesLitterMod            , only : ncwd
   use FatesLitterMod            , only : ndcmpy
   use EDPftvarcon               , only : FatesReportPFTParams
   use EDPftvarcon               , only : FatesCheckParams
   use EDPftvarcon               , only : EDPftvarcon_inst
   use SFParamsMod               , only : SpitFireCheckParams
   use EDParamsMod               , only : FatesReportParams
   use EDParamsMod               , only : bgc_soil_salinity
   use FatesPlantHydraulicsMod   , only : InitHydroGlobals
   use EDParamsMod               , only : photo_temp_acclim_timescale
   use EDParamsMod               , only : photo_temp_acclim_thome_time
   use EDParamsMod               , only : ED_val_history_sizeclass_bin_edges
   use EDParamsMod               , only : ED_val_history_ageclass_bin_edges
   use EDParamsMod               , only : ED_val_history_height_bin_edges
   use EDParamsMod               , only : ED_val_history_coageclass_bin_edges
   use CLMFatesParamInterfaceMod , only : FatesReadParameters
   use EDTypesMod                , only : p_uptake_mode
   use EDTypesMod                , only : n_uptake_mode
   use EDTypesMod                , only : ed_site_type
   use FatesConstantsMod         , only : prescribed_p_uptake
   use FatesConstantsMod         , only : prescribed_n_uptake
   use FatesConstantsMod         , only : coupled_p_uptake
   use FatesConstantsMod         , only : coupled_n_uptake
   use FatesConstantsMod         , only : fates_np_comp_scaling
   use FatesConstantsMod         , only : coupled_np_comp_scaling
   use FatesConstantsMod         , only : trivial_np_comp_scaling
   use PRTGenericMod             , only : num_elements
   use PRTGenericMod             , only : element_list
   use PRTGenericMod             , only : element_pos
   use EDParamsMod               , only : eca_plant_escalar
   use PRTGenericMod             , only : prt_carbon_allom_hyp
   use PRTGenericMod             , only : prt_cnp_flex_allom_hyp
   use PRTGenericMod             , only : carbon12_element
   use PRTGenericMod             , only : nitrogen_element
   use PRTGenericMod             , only : phosphorus_element
   use PRTGenericMod             , only : num_organ_types
   use PRTGenericMod             , only : leaf_organ, fnrt_organ, store_organ
   use PRTGenericMod             , only : sapw_organ, struct_organ, repro_organ
   use PRTParametersMod          , only : prt_params
   use PRTInitParamsFatesMod     , only : PRTCheckParams, PRTDerivedParams
   use PRTAllometricCarbonMod    , only : InitPRTGlobalAllometricCarbon
   use PRTAllometricCNPMod       , only : InitPRTGlobalAllometricCNP
   use FatesRunningMeanMod       , only : ema_24hr
   use FatesRunningMeanMod       , only : fixed_24hr
   use FatesRunningMeanMod       , only : ema_lpa
   use FatesRunningMeanMod       , only : ema_longterm
   use FatesRunningMeanMod       , only : ema_60day
   use FatesRunningMeanMod       , only : moving_ema_window
   use FatesRunningMeanMod       , only : fixed_window
   use FatesHistoryInterfaceMod  , only : fates_hist
   use FatesHydraulicsMemMod     , only : nshell
   use FatesHydraulicsMemMod     , only : nlevsoi_hyd_max
   
   ! CIME Globals
   use shr_log_mod               , only : errMsg => shr_log_errMsg
   use shr_infnan_mod            , only : nan => shr_infnan_nan, assignment(=)

   ! Just use everything from FatesInterfaceTypesMod, this is
   ! its sister code
   use FatesInterfaceTypesMod

   implicit none

   private

   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), pointer :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)


      ! These are parameter constants that FATES may need to provide a host model
      ! We have other methods of reading in input parameters. Since these
      ! are parameter constants, we don't need them allocated over every site,one
      ! instance is fine.
      
      type(bc_pconst_type) :: bc_pconst
      

   end type fates_interface_type
   
   

   character(len=*), parameter :: sourcefile = &
        __FILE__

   ! Make public necessary subroutines and functions
   public :: FatesInterfaceInit
   public :: set_fates_ctrlparms
   public :: SetFatesTime
   public :: SetFatesGlobalElements1
   public :: SetFatesGlobalElements2
   public :: FatesReportParameters
   public :: allocate_bcin
   public :: allocate_bcout
   public :: allocate_bcpconst
   public :: set_bcpconst
   public :: zero_bcs
   public :: set_bcs
   public :: UpdateFatesRMeansTStep
   public :: InitTimeAveragingGlobals
   
contains

  ! ====================================================================================
  subroutine FatesInterfaceInit(log_unit,global_verbose)
    
    use FatesGlobals, only : FatesGlobalsInit
    
    implicit none
    
    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    call FatesGlobalsInit(log_unit,global_verbose)
    
  end subroutine FatesInterfaceInit

  ! ====================================================================================
  
  ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
  subroutine fates_clean(this)
      
    implicit none
    
    ! Input Arguments
    class(fates_interface_type), intent(inout) :: this
    
    ! Incrementally walk through linked list and deallocate
    
    
      
    ! Deallocate the site list
    !      deallocate (this%sites)
      
    return
  end subroutine fates_clean
  

  ! ====================================================================================

   
  subroutine allocate_bcpconst(bc_pconst,nlevdecomp)
    
    type(bc_pconst_type), intent(inout) :: bc_pconst
    integer             , intent(in)    :: nlevdecomp 

    allocate(bc_pconst%vmax_nh4(numpft))
    allocate(bc_pconst%vmax_no3(numpft))    
    allocate(bc_pconst%vmax_p(numpft))
    allocate(bc_pconst%eca_km_nh4(numpft))
    allocate(bc_pconst%eca_km_no3(numpft))
    allocate(bc_pconst%eca_km_p(numpft))
    allocate(bc_pconst%eca_km_ptase(numpft))
    allocate(bc_pconst%eca_vmax_ptase(numpft))
    allocate(bc_pconst%eca_alpha_ptase(numpft))
    allocate(bc_pconst%eca_lambda_ptase(numpft))
    allocate(bc_pconst%j_uptake(nlevdecomp))
    
    return
  end subroutine allocate_bcpconst
  
  ! ====================================================================================
  
  subroutine set_bcpconst(bc_pconst,nlevdecomp)

    type(bc_pconst_type), intent(inout) :: bc_pconst
    integer             , intent(in)    :: nlevdecomp 
    integer                             :: j
    
    bc_pconst%vmax_nh4(1:numpft)         = EDPftvarcon_inst%vmax_nh4(1:numpft)
    bc_pconst%vmax_no3(1:numpft)         = EDPftvarcon_inst%vmax_no3(1:numpft)
    bc_pconst%vmax_p(1:numpft)           = EDPftvarcon_inst%vmax_p(1:numpft)
    
    bc_pconst%eca_km_nh4(1:numpft)       = EDPftvarcon_inst%eca_km_nh4(1:numpft)
    bc_pconst%eca_km_no3(1:numpft)       = EDPftvarcon_inst%eca_km_no3(1:numpft)
    bc_pconst%eca_km_p(1:numpft)         = EDPftvarcon_inst%eca_km_p(1:numpft)
    bc_pconst%eca_km_ptase(1:numpft)     = EDPftvarcon_inst%eca_km_ptase(1:numpft)
    bc_pconst%eca_vmax_ptase(1:numpft)   = EDPftvarcon_inst%eca_vmax_ptase(1:numpft)
    bc_pconst%eca_alpha_ptase(1:numpft)  = EDPftvarcon_inst%eca_alpha_ptase(1:numpft) 
    bc_pconst%eca_lambda_ptase(1:numpft) = EDPftvarcon_inst%eca_lambda_ptase(1:numpft)
    bc_pconst%eca_plant_escalar          = eca_plant_escalar
    
    return
  end subroutine set_bcpconst

  ! ====================================================================================
   
  subroutine zero_bcs(fates,s)

    type(fates_interface_type), intent(inout) :: fates
    integer, intent(in) :: s
    
    ! Input boundaries
    
    fates%bc_in(s)%lightning24(:)      = 0.0_r8
    fates%bc_in(s)%pop_density(:)      = 0.0_r8
    fates%bc_in(s)%precip24_pa(:)      = 0.0_r8
    fates%bc_in(s)%relhumid24_pa(:)    = 0.0_r8
    fates%bc_in(s)%wind24_pa(:)        = 0.0_r8
     
    fates%bc_in(s)%solad_parb(:,:)     = 0.0_r8
    fates%bc_in(s)%solai_parb(:,:)     = 0.0_r8
    fates%bc_in(s)%smp_sl(:)           = 0.0_r8
    fates%bc_in(s)%eff_porosity_sl(:)  = 0.0_r8
    fates%bc_in(s)%watsat_sl(:)        = 0.0_r8
    fates%bc_in(s)%tempk_sl(:)         = 0.0_r8
    fates%bc_in(s)%h2o_liqvol_sl(:)    = 0.0_r8
    fates%bc_in(s)%filter_vegzen_pa(:) = .false.
    fates%bc_in(s)%coszen_pa(:)        = 0.0_r8
    fates%bc_in(s)%fcansno_pa(:)       = 0.0_r8
    fates%bc_in(s)%albgr_dir_rb(:)     = 0.0_r8
    fates%bc_in(s)%albgr_dif_rb(:)     = 0.0_r8
    fates%bc_in(s)%max_rooting_depth_index_col = 0
    fates%bc_in(s)%tot_het_resp        = 0.0_r8
    fates%bc_in(s)%tot_somc            = 0.0_r8 
    fates%bc_in(s)%tot_litc            = 0.0_r8
    fates%bc_in(s)%snow_depth_si       = 0.0_r8
    fates%bc_in(s)%frac_sno_eff_si     = 0.0_r8
    fates%bc_in(s)%w_scalar_sisl(:)    = 0.0_r8
    fates%bc_in(s)%t_scalar_sisl(:)    = 0.0_r8
    
    if(do_fates_salinity)then
       fates%bc_in(s)%salinity_sl(:)   = 0.0_r8
    endif
    
    if (hlm_use_planthydro.eq.itrue) then
       
       fates%bc_in(s)%qflx_transp_pa(:) = 0.0_r8
       fates%bc_in(s)%swrad_net_pa(:) = 0.0_r8
       fates%bc_in(s)%lwrad_net_pa(:) = 0.0_r8
       fates%bc_in(s)%watsat_sisl(:) = 0.0_r8
       fates%bc_in(s)%watres_sisl(:) = 0.0_r8
       fates%bc_in(s)%sucsat_sisl(:) = 0.0_r8
       fates%bc_in(s)%bsw_sisl(:) = 0.0_r8
       fates%bc_in(s)%hksat_sisl(:) = 0.0_r8
    end if

    
    ! Output boundaries
    fates%bc_out(s)%active_suction_sl(:) = .false.
    fates%bc_out(s)%fsun_pa(:)      = 0.0_r8
    fates%bc_out(s)%laisun_pa(:)    = 0.0_r8
    fates%bc_out(s)%laisha_pa(:)    = 0.0_r8
    fates%bc_out(s)%rootr_pasl(:,:) = 0.0_r8
    fates%bc_out(s)%btran_pa(:)     = 0.0_r8

    ! MIMIC litter quality, always initialize to unset
    fates%bc_out(s)%litt_flux_ligc_per_n = fates_unset_r8

    
    ! Fates -> BGC fragmentation mass fluxes
    select case(hlm_parteh_mode) 
    case(prt_carbon_allom_hyp)
       fates%bc_out(s)%litt_flux_cel_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lig_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lab_c_si(:) = 0._r8
    case(prt_cnp_flex_allom_hyp) 
       
       fates%bc_in(s)%plant_nh4_uptake_flux(:,:) = 0._r8
       fates%bc_in(s)%plant_no3_uptake_flux(:,:) = 0._r8
       fates%bc_in(s)%plant_p_uptake_flux(:,:) = 0._r8
       fates%bc_out(s)%source_p(:)           = 0._r8
       fates%bc_out(s)%source_nh4(:)         = 0._r8
       fates%bc_out(s)%litt_flux_cel_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lig_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lab_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_cel_n_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lig_n_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lab_n_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_cel_p_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lig_p_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lab_p_si(:) = 0._r8
       
    case default
       write(fates_log(), *) 'An unknown parteh hypothesis was passed'
       write(fates_log(), *) 'while zeroing output boundary conditions'
       write(fates_log(), *) 'hlm_parteh_mode: ',hlm_parteh_mode
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    
    fates%bc_out(s)%rssun_pa(:)     = 0.0_r8
    fates%bc_out(s)%rssha_pa(:)     = 0.0_r8
    
    fates%bc_out(s)%albd_parb(:,:) = 0.0_r8
    fates%bc_out(s)%albi_parb(:,:) = 0.0_r8
    fates%bc_out(s)%fabd_parb(:,:) = 0.0_r8
    fates%bc_out(s)%fabi_parb(:,:) = 0.0_r8
    fates%bc_out(s)%ftdd_parb(:,:) = 0.0_r8
    fates%bc_out(s)%ftid_parb(:,:) = 0.0_r8
    fates%bc_out(s)%ftii_parb(:,:) = 0.0_r8
    
    fates%bc_out(s)%elai_pa(:)   = 0.0_r8
    fates%bc_out(s)%esai_pa(:)   = 0.0_r8
    fates%bc_out(s)%tlai_pa(:)   = 0.0_r8
    fates%bc_out(s)%tsai_pa(:)   = 0.0_r8
    fates%bc_out(s)%htop_pa(:)   = 0.0_r8
    fates%bc_out(s)%hbot_pa(:)   = 0.0_r8
    fates%bc_out(s)%displa_pa(:) = 0.0_r8
    fates%bc_out(s)%z0m_pa(:)    = 0.0_r8
    fates%bc_out(s)%dleaf_pa(:)   = 0.0_r8
    fates%bc_out(s)%nocomp_pft_label_pa(:) = 0
    
    fates%bc_out(s)%canopy_fraction_pa(:) = 0.0_r8
    fates%bc_out(s)%frac_veg_nosno_alb_pa(:) = 0.0_r8
    
    if (hlm_use_planthydro.eq.itrue) then
       fates%bc_out(s)%qflx_soil2root_sisl(:) = 0.0_r8
       fates%bc_out(s)%qflx_ro_sisl(:)        = 0.0_r8
    end if
    fates%bc_out(s)%plant_stored_h2o_si = 0.0_r8

    ! Land Use realated
    fates%bc_out(s)%gpp_site = 0.0_r8
    fates%bc_out(s)%ar_site = 0.0_r8
    fates%bc_out(s)%hrv_deadstemc_to_prod10c = 0.0_r8
    fates%bc_out(s)%hrv_deadstemc_to_prod100c = 0.0_r8
    
    return
  end subroutine zero_bcs

  ! ===========================================================================

   subroutine allocate_bcin(bc_in, nlevsoil_in, nlevdecomp_in, num_lu_harvest_cats,natpft_lb,natpft_ub)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      integer,intent(in)              :: nlevsoil_in
      integer,intent(in)              :: nlevdecomp_in
      integer,intent(in)              :: num_lu_harvest_cats
      integer,intent(in)              :: natpft_lb,natpft_ub ! dimension bounds of the array holding surface file pft data
      
      ! Allocate input boundaries

      bc_in%nlevsoil   = nlevsoil_in

      if(nlevsoil_in > numlevsoil_max) then
         write(fates_log(), *) 'The number of soil layers imposed by the host model'
         write(fates_log(), *) 'is larger than what we have allocated in our static'
         write(fates_log(), *) 'arrays. Please increase the size of numlevsoil_max'
         write(fates_log(), *) 'found in EDTypesMod.F90'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      bc_in%nlevdecomp = nlevdecomp_in

      if (hlm_use_vertsoilc == itrue) then
         if(bc_in%nlevdecomp .ne. bc_in%nlevsoil) then
            write(fates_log(), *) 'The host has signaled a vertically resolved'
            write(fates_log(), *) 'soil decomposition model. Therefore, the '
            write(fates_log(), *) 'total number of soil layers should equal the'
            write(fates_log(), *) 'total number of decomposition layers.'
            write(fates_log(), *) 'nlevdecomp: ',bc_in%nlevdecomp
            write(fates_log(), *) 'nlevsoil: ',bc_in%nlevsoil
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      else
         if(bc_in%nlevdecomp .ne. 1)then
            write(fates_log(), *) 'The host has signaled a non-vertically resolved'
            write(fates_log(), *) 'soil decomposition model. Therefore, the '
            write(fates_log(), *) 'total number of decomposition layers should be 1.'
            write(fates_log(), *) 'nlevdecomp: ',bc_in%nlevdecomp
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      end if

      ! Plant Nutrient Aquisition variables
      ! If we are up-scaling to PFT, then we need to pass bach PFTxlayer
      ! if we don't, then there is ambiguity in the uptake. If we
      ! do not upscale to PFT, then we can simply send back the
      ! uptake for each cohort, and don't need to allocate by layer
      ! Allocating differently could save a lot of memory and time

      if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp) then
         allocate(bc_in%plant_nh4_uptake_flux(max_comp_per_site,1))
         allocate(bc_in%plant_no3_uptake_flux(max_comp_per_site,1))
         allocate(bc_in%plant_p_uptake_flux(max_comp_per_site,1))
      else
         allocate(bc_in%plant_nh4_uptake_flux(1,1))
         allocate(bc_in%plant_no3_uptake_flux(1,1))
         allocate(bc_in%plant_p_uptake_flux(1,1))
      end if


      allocate(bc_in%zi_sisl(0:nlevsoil_in))
      allocate(bc_in%dz_sisl(nlevsoil_in))
      allocate(bc_in%z_sisl(nlevsoil_in))
      allocate(bc_in%decomp_id(nlevsoil_in))
      allocate(bc_in%dz_decomp_sisl(nlevdecomp_in))
      allocate(bc_in%w_scalar_sisl(nlevsoil_in))
      allocate(bc_in%t_scalar_sisl(nlevsoil_in))

      ! Lightning (or successful ignitions) and population density
      ! Fire related variables
      allocate(bc_in%lightning24(maxpatch_total))
      allocate(bc_in%pop_density(maxpatch_total))
      allocate(bc_in%wind24_pa(maxpatch_total))
      allocate(bc_in%relhumid24_pa(maxpatch_total))
      allocate(bc_in%precip24_pa(maxpatch_total))
      
      ! Radiation
      allocate(bc_in%solad_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_in%solai_parb(maxpatch_total,hlm_numSWb))
      
      ! Hydrology
      allocate(bc_in%smp_sl(nlevsoil_in))
      allocate(bc_in%eff_porosity_sl(nlevsoil_in))
      allocate(bc_in%watsat_sl(nlevsoil_in))
      allocate(bc_in%tempk_sl(nlevsoil_in))
      allocate(bc_in%h2o_liqvol_sl(nlevsoil_in))
      
      !BGC
      if(do_fates_salinity) then
         allocate(bc_in%salinity_sl(nlevsoil_in))
      endif

      
      
      ! Photosynthesis
      allocate(bc_in%filter_photo_pa(maxpatch_total))
      allocate(bc_in%dayl_factor_pa(maxpatch_total))
      allocate(bc_in%esat_tv_pa(maxpatch_total))
      allocate(bc_in%eair_pa(maxpatch_total))
      allocate(bc_in%oair_pa(maxpatch_total))
      allocate(bc_in%cair_pa(maxpatch_total))
      allocate(bc_in%rb_pa(maxpatch_total))
      allocate(bc_in%t_veg_pa(maxpatch_total))
      allocate(bc_in%tgcm_pa(maxpatch_total))
      allocate(bc_in%t_soisno_sl(nlevsoil_in))

      ! Canopy Radiation
      allocate(bc_in%filter_vegzen_pa(maxpatch_total))
      allocate(bc_in%coszen_pa(maxpatch_total))
      allocate(bc_in%fcansno_pa(maxpatch_total))
      allocate(bc_in%albgr_dir_rb(hlm_numSWb))
      allocate(bc_in%albgr_dif_rb(hlm_numSWb))

      ! Plant-Hydro BC's
      if (hlm_use_planthydro.eq.itrue) then

         allocate(bc_in%qflx_transp_pa(maxpatch_total))
         allocate(bc_in%swrad_net_pa(maxpatch_total))
         allocate(bc_in%lwrad_net_pa(maxpatch_total))
         
         allocate(bc_in%watsat_sisl(nlevsoil_in))
         allocate(bc_in%watres_sisl(nlevsoil_in))
         allocate(bc_in%sucsat_sisl(nlevsoil_in))
         allocate(bc_in%bsw_sisl(nlevsoil_in))
         allocate(bc_in%hksat_sisl(nlevsoil_in))
         allocate(bc_in%h2o_liq_sisl(nlevsoil_in)); bc_in%h2o_liq_sisl = nan
      end if

      ! Land use

      ! harvest flag denote data from hlm,
      ! while the logging flag signifies only that logging is occurring (which could just be FATES logging)
      if (hlm_use_lu_harvest .gt. 0) then
         allocate(bc_in%hlm_harvest_rates(num_lu_harvest_cats))
         allocate(bc_in%hlm_harvest_catnames(num_lu_harvest_cats))
      else ! LoggingMortality_frac needs these passed to it regardless of harvest
         allocate(bc_in%hlm_harvest_rates(0))
         allocate(bc_in%hlm_harvest_catnames(0))
      end if

      allocate(bc_in%pft_areafrac(natpft_lb:natpft_ub))

      ! Variables for SP mode. 
      if(hlm_use_sp.eq.itrue) then
        allocate(bc_in%hlm_sp_tlai(natpft_lb:natpft_ub))
        allocate(bc_in%hlm_sp_tsai(natpft_lb:natpft_ub))     
        allocate(bc_in%hlm_sp_htop(natpft_lb:natpft_ub))
      end if 
      return
   end subroutine allocate_bcin

   ! ====================================================================================
   
   subroutine allocate_bcout(bc_out, nlevsoil_in, nlevdecomp_in)

      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_out_type), intent(inout) :: bc_out
      integer,intent(in)               :: nlevsoil_in
      integer,intent(in)               :: nlevdecomp_in
      
      ! Radiation
      allocate(bc_out%fsun_pa(maxpatch_total))
      allocate(bc_out%laisun_pa(maxpatch_total))
      allocate(bc_out%laisha_pa(maxpatch_total))
      
      ! Hydrology
      allocate(bc_out%active_suction_sl(nlevsoil_in))
      allocate(bc_out%rootr_pasl(maxpatch_total,nlevsoil_in))
      allocate(bc_out%btran_pa(maxpatch_total))
      
      ! Photosynthesis

      allocate(bc_out%rssun_pa(maxpatch_total))
      allocate(bc_out%rssha_pa(maxpatch_total))
      
      ! Canopy Radiation
      allocate(bc_out%albd_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%albi_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%fabd_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%fabi_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%ftdd_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%ftid_parb(maxpatch_total,hlm_numSWb))
      allocate(bc_out%ftii_parb(maxpatch_total,hlm_numSWb))


      ! We allocate the boundary conditions to the BGC
      ! model, regardless of what scheme we use. The BGC
      ! model in ELM allocates all species C,N,P even if they
      ! are not turned on. Also, it is feasible that the
      ! one would want to allow soil BGC nutrient dynamics
      ! to proceed even if we are not passing source fluxes
      ! or uptake from FATES.
      ! When FATES does not have nutrients enabled, these
      ! arrays are indexed by 1.
      
      !if(trim(hlm_nu_com).eq.'RD') then
      !   allocate(bc_out%n_demand(max_comp_per_site))
      !   allocate(bc_out%p_demand(max_comp_per_site))
      !end if

      ! Used in both
      allocate(bc_out%veg_rootc(max_comp_per_site,nlevdecomp_in))
      allocate(bc_out%ft_index(max_comp_per_site))
         
      if(trim(hlm_nu_com).eq.'ECA') then
         allocate(bc_out%decompmicc(nlevdecomp_in))
         allocate(bc_out%cn_scalar(max_comp_per_site))
         allocate(bc_out%cp_scalar(max_comp_per_site))
      end if

      ! Include the bare-ground patch for these patch-level boundary conditions
      ! (it will always be zero for all of these)
      if(hlm_use_ch4.eq.itrue) then
         allocate(bc_out%annavg_agnpp_pa(0:maxpatch_total));bc_out%annavg_agnpp_pa(:)=nan
         allocate(bc_out%annavg_bgnpp_pa(0:maxpatch_total));bc_out%annavg_bgnpp_pa(:)=nan
         allocate(bc_out%annsum_npp_pa(0:maxpatch_total));bc_out%annsum_npp_pa(:)=nan
         allocate(bc_out%frootc_pa(0:maxpatch_total));bc_out%frootc_pa(:)=nan
         allocate(bc_out%root_resp(nlevsoil_in));bc_out%root_resp(:)=nan
         allocate(bc_out%woody_frac_aere_pa(0:maxpatch_total));bc_out%woody_frac_aere_pa(:)=nan
         allocate(bc_out%rootfr_pa(0:maxpatch_total,nlevsoil_in))
         bc_out%rootfr_pa(:,:)=nan

         ! Give the bare-ground root fractions a nominal fraction of unity over depth
         bc_out%rootfr_pa(0,1:nlevsoil_in)=1._r8/real(nlevsoil_in,r8)
      end if

      bc_out%ema_npp = nan
      
      
      ! Fates -> BGC fragmentation mass fluxes
      select case(hlm_parteh_mode) 
      case(prt_carbon_allom_hyp)
         allocate(bc_out%litt_flux_cel_c_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lig_c_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lab_c_si(nlevdecomp_in))
      case(prt_cnp_flex_allom_hyp) 

         
         allocate(bc_out%litt_flux_cel_c_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lig_c_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lab_c_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_cel_n_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lig_n_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lab_n_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_cel_p_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lig_p_si(nlevdecomp_in))
         allocate(bc_out%litt_flux_lab_p_si(nlevdecomp_in))

         allocate(bc_out%source_nh4(nlevdecomp_in))
         allocate(bc_out%source_p(nlevdecomp_in))

      case default
         write(fates_log(), *) 'An unknown parteh hypothesis was passed'
         write(fates_log(), *) 'to the site level output boundary conditions'
         write(fates_log(), *) 'hlm_parteh_mode: ',hlm_parteh_mode
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select


      ! Canopy Structure
      allocate(bc_out%elai_pa(maxpatch_total))
      allocate(bc_out%esai_pa(maxpatch_total))
      allocate(bc_out%tlai_pa(maxpatch_total))
      allocate(bc_out%tsai_pa(maxpatch_total))
      allocate(bc_out%htop_pa(maxpatch_total))
      allocate(bc_out%hbot_pa(maxpatch_total))
      allocate(bc_out%dleaf_pa(maxpatch_total))

      allocate(bc_out%displa_pa(maxpatch_total))
      allocate(bc_out%z0m_pa(maxpatch_total))

      allocate(bc_out%canopy_fraction_pa(maxpatch_total))
      allocate(bc_out%frac_veg_nosno_alb_pa(maxpatch_total))

      allocate(bc_out%nocomp_pft_label_pa(maxpatch_total))

      ! Plant-Hydro BC's
      if (hlm_use_planthydro.eq.itrue) then
         allocate(bc_out%qflx_soil2root_sisl(nlevsoil_in))
         allocate(bc_out%qflx_ro_sisl(nlevsoil_in))
      end if

      return
   end subroutine allocate_bcout

   ! ====================================================================================

   subroutine set_bcs(bc_in)

       ! --------------------------------------------------------------------------------
       !
       ! This subroutine is called directly from the HLM to set boundary condition not yet 
       !     functional from hlm. This allows flexibility for model testing.
       !
       ! This subroutine MUST BE CALLED AFTER the FATES PFT parameter file has been read in,
       ! and the EDPftvarcon_inst structure has been made.
       ! This subroutine must ALSO BE CALLED BEFORE the history file dimensions
       ! are set.
       ! 
       ! --------------------------------------------------------------------------------
      implicit none
      type(bc_in_type), intent(inout) :: bc_in

      ! Input boundaries
      ! Warning: these "z" type variables
      ! are written only once at the beginning
      ! so THIS ROUTINE SHOULD NOT BE CALLED AFTER
      ! INITIALIZATION
      if(do_fates_salinity)then
         bc_in%salinity_sl(:)     = bgc_soil_salinity
      endif
      
    end subroutine set_bcs

    ! ===================================================================================
    
    subroutine SetFatesGlobalElements1(use_fates,surf_numpft,surf_numcft)

       ! --------------------------------------------------------------------------------
       !
       ! This is the first FATES routine that is called.
       !
       ! spmode,biogeog and nocomp mode flags have been passed prior to this call
       ! --------------------------------------------------------------------------------


      implicit none
      
      logical,intent(in) :: use_fates    ! Is fates turned on?
      integer,intent(in) :: surf_numpft  ! Number of PFTs in surface dataset
      integer,intent(in) :: surf_numcft  ! Number of CFTs in surface dataset

      integer :: fates_numpft  ! Number of PFTs tracked in FATES
      
      if (use_fates) then
         
         ! Self explanatory, read the fates parameter file
         call FatesReadParameters()

         fates_numpft = size(prt_params%wood_density,dim=1)
         
         if(hlm_use_sp==itrue)then

            ! For an SP run we also just use the primary patches
            ! to hold all PFTs.  So create the same number of
            ! patches as the number of PFTs

            maxpatch_primary   = fates_numpft
            maxpatch_secondary = 0
            maxpatch_total     = fates_numpft
            
            ! If this is an SP run, we actually need enough patches on the
            ! CLM/ELM side of the code to hold the LAI data.  This
            ! number may be larger than what fates requires.  Of course
            ! we may have multiple PFTs in the surface datafile mapping
            ! to FATES.  The surf_numpft includes the bare ground.
            ! maxpatch_total does not include the bare ground (so add 1)
            
            fates_maxPatchesPerSite = max(surf_numpft+surf_numcft,maxpatch_total+1)

         else

            ! If we are using fixed biogeography or no-comp then we
            ! can also apply those constraints to maxpatch_primary and secondary
            ! and that value will match fates_maxPatchesPerSite
            
            if(hlm_use_nocomp==itrue) then

               maxpatch_primary = max(maxpatch_primary,fates_numpft)
               maxpatch_total = maxpatch_primary + maxpatch_secondary
               !if(maxpatch_primary<fates_numpft)then
               !   write(fates_log(),*) 'warning: lower number of patches than pfts'
               !   write(fates_log(),*) 'this may become a problem in nocomp mode'
               !end if
            end if

            ! maxpatch_total does not include the bare ground (so add 1)
            fates_maxPatchesPerSite = maxpatch_total+1
            
         end if
             
      end if

    end subroutine SetFatesGlobalElements1

    ! ====================================================================================
    
    subroutine SetFatesGlobalElements2(use_fates)

      ! --------------------------------------------------------------------------------
      !
      ! This is the second FATES routine that is called.
      !
      ! --------------------------------------------------------------------------------

      logical,intent(in) :: use_fates    ! Is fates turned on?
      integer :: i

      if (use_fates) then

         if(lbound(prt_params%wood_density(:),dim=1) .eq. 0 ) then
            numpft = size(prt_params%wood_density,dim=1)-1
         elseif(lbound(prt_params%wood_density(:),dim=1) .eq. 1 ) then
            numpft = size(prt_params%wood_density,dim=1)
         else
            write(fates_log(), *) 'While assessing the number of FATES PFTs,'
            write(fates_log(), *) 'it was found that the lower bound was neither 0 or 1?'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(numpft>maxpft) then
            write(fates_log(), *) 'The number of PFTs dictated by the FATES parameter file'
            write(fates_log(), *) 'is larger than the maximum allowed. Increase the FATES parameter constant'
            write(fates_log(), *) 'FatesInterfaceMod.F90:maxpft accordingly'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         ! Identify the number of leaf age-classes
         
         if( (lbound(prt_params%leaf_long(:,:),dim=2) .eq. 0) .or. &
             (ubound(prt_params%leaf_long(:,:),dim=2) .eq. 0) ) then
            write(fates_log(), *) 'While assessing the number of FATES leaf age classes,'
            write(fates_log(), *) 'The second dimension of leaf_long was 0?'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         else
            nleafage = size(prt_params%leaf_long,dim=2)
         end if

         
         ! These values are used to define the restart file allocations and general structure
         ! of memory for the cohort arrays
         if(hlm_use_sp.eq.itrue) then
            fates_maxElementsPerPatch = maxSWb
         else
            fates_maxElementsPerPatch = max(maxSWb,max_cohort_per_patch, ndcmpy*hlm_maxlevsoil ,ncwd*hlm_maxlevsoil)
         end if
         
         fates_maxElementsPerSite = max(fates_maxPatchesPerSite * fates_maxElementsPerPatch, &
              numWatermem, num_vegtemp_mem, num_elements, nlevsclass*numpft)

         if(hlm_use_planthydro==itrue)then
            fates_maxElementsPerSite = max(fates_maxElementsPerSite, nshell*nlevsoi_hyd_max )
         end if
         
         
         ! Set the maximum number of nutrient aquisition competitors per site
         ! This is used to set array sizes for the boundary conditions.
         ! Note: since BGC code may be active even when no nutrients
         ! present, we still need to allocate things when no nutrients


         if (any(abs(EDPftvarcon_inst%prescribed_nuptake(:)) > nearzero )) then
            n_uptake_mode = prescribed_n_uptake
         else
            n_uptake_mode = coupled_n_uptake
         end if

         if (any(abs(EDPftvarcon_inst%prescribed_puptake(:)) > nearzero )) then
            p_uptake_mode = prescribed_p_uptake
         else
            p_uptake_mode = coupled_p_uptake
         end if
         
         if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp ) then

            if((p_uptake_mode==coupled_p_uptake) .or. (n_uptake_mode==coupled_n_uptake))then
               max_comp_per_site = fates_maxElementsPerSite
               fates_np_comp_scaling = coupled_np_comp_scaling
            else
               max_comp_per_site = 1
               fates_np_comp_scaling = trivial_np_comp_scaling
            end if

         else
            max_comp_per_site = 1
            fates_np_comp_scaling = trivial_np_comp_scaling
         end if
            
         ! calculate the bin edges for radiative transfer calculations
         ! VAI bin widths array 
         do i = 1,nlevleaf
            dinc_vai(i) = ED_val_vai_top_bin_width * ED_val_vai_width_increase_factor ** (i-1)
         end do

         ! lower edges of VAI bins       
         do i = 1,nlevleaf
            dlower_vai(i) = sum(dinc_vai(1:i))
         end do

         ! Identify number of size and age class bins for history output
         ! assume these arrays are 1-indexed
         nlevsclass = size(ED_val_history_sizeclass_bin_edges,dim=1)
         nlevage = size(ED_val_history_ageclass_bin_edges,dim=1)
         nlevheight = size(ED_val_history_height_bin_edges,dim=1)
         nlevcoage = size(ED_val_history_coageclass_bin_edges,dim=1)
         nlevdamage = size(ED_val_history_damage_bin_edges, dim=1)
         
         ! do some checks on the size, age, and height bin arrays to make sure they make sense:
         ! make sure that all start at zero, and that both are monotonically increasing
         if ( ED_val_history_sizeclass_bin_edges(1) .ne. 0._r8 ) then
            write(fates_log(), *) 'size class bins specified in parameter file must start at zero'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
         if ( ED_val_history_ageclass_bin_edges(1) .ne. 0._r8 ) then
            write(fates_log(), *) 'age class bins specified in parameter file must start at zero'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
         if ( ED_val_history_height_bin_edges(1) .ne. 0._r8 ) then
            write(fates_log(), *) 'height class bins specified in parameter file must start at zero'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
         do i = 2,nlevsclass
            if ( (ED_val_history_sizeclass_bin_edges(i) - ED_val_history_sizeclass_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'size class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do
         do i = 2,nlevage
            if ( (ED_val_history_ageclass_bin_edges(i) - ED_val_history_ageclass_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'age class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do
         do i = 2,nlevheight
            if ( (ED_val_history_height_bin_edges(i) - ED_val_history_height_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'height class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do
         do i = 2,nlevcoage
            if ( (ED_val_history_coageclass_bin_edges(i) - ED_val_history_coageclass_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'cohort age class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do

         ! Initialize Hydro globals 
         ! (like water retention functions)
         ! this needs to know the number of PFTs, which is
         ! determined in that call
         call InitHydroGlobals()
   
         ! Initialize the Plant Allocation and Reactive Transport
         ! global functions and mapping tables
         ! Also associate the elements defined in PARTEH with a list in FATES
         ! "element_list" is useful because it allows the fates side of the code
         ! to loop through elements, and call the correct PARTEH interfaces
         ! automatically.
         call InitPARTEHGlobals()
         
         
         ! Set Various Mapping Arrays used in history output as well
         ! These will not be used if use_ed or use_fates is false
         call fates_history_maps()

       

      else
         ! If we are not using FATES, the cohort dimension is still
         ! going to be initialized, lets set it to the smallest value
         ! possible so that the dimensioning info takes up little space

         fates_maxElementsPerPatch = 1
      
         fates_maxElementsPerSite = 1
         

      end if
    end subroutine SetFatesGlobalElements2

    ! ======================================================================

    subroutine InitTimeAveragingGlobals()
      
      ! Instantiate the time-averaging method globals
      ! NOTE: It may be possible in the future that the HLM model timesteps
      ! are dynamic in time or space, in that case, these would no longer
      ! be global constants.

      allocate(ema_24hr)
      call ema_24hr%define(sec_per_day, hlm_stepsize, moving_ema_window)
      allocate(fixed_24hr)
      call fixed_24hr%define(sec_per_day, hlm_stepsize, fixed_window)
      allocate(ema_lpa)  ! note that this parameter has units of days
      call ema_lpa%define(photo_temp_acclim_timescale*sec_per_day, &
           hlm_stepsize,moving_ema_window)
      allocate(ema_longterm)  ! note that this parameter has units of years
      call ema_longterm%define(photo_temp_acclim_thome_time*days_per_year*sec_per_day, & 
           hlm_stepsize,moving_ema_window)
      
      !allocate(ema_60day)
      !call ema_60day%define(prt_params%fnrt_adapt_tscl*sec_per_day,sec_per_day,moving_ema_window)
      !class(rmean_arr_type), pointer :: ema_fnrt_tscale(:)
      !rmean_arr_type
      
      
      return
    end subroutine InitTimeAveragingGlobals

      
    ! ======================================================================
    
    subroutine InitPARTEHGlobals()

     ! Initialize the Plant Allocation and Reactive Transport
     ! global functions and mapping tables
     ! Also associate the elements defined in PARTEH with a list in FATES
     ! "element_list" is useful because it allows the fates side of the code
     ! to loop through elements, and call the correct PARTEH interfaces
     ! automatically.
     
     select case(hlm_parteh_mode)
     case(prt_carbon_allom_hyp)

        num_elements = 1
        allocate(element_list(num_elements))
        element_list(1) = carbon12_element
        element_pos(:) = 0
        element_pos(carbon12_element) = 1

        call InitPRTGlobalAllometricCarbon()

     case(prt_cnp_flex_allom_hyp)
        
        num_elements = 3
        allocate(element_list(num_elements))
        element_list(1) = carbon12_element
        element_list(2) = nitrogen_element
        element_list(3) = phosphorus_element
        element_pos(:)  = 0
        element_pos(carbon12_element)   = 1
        element_pos(nitrogen_element)   = 2
        element_pos(phosphorus_element) = 3

        call InitPRTGlobalAllometricCNP()
        
     case DEFAULT
        write(fates_log(),*) 'You specified an unknown PRT module'
        write(fates_log(),*) 'Check your setting for fates_parteh_mode'
        write(fates_log(),*) 'in the CLM namelist. The only valid value now is 1'
        write(fates_log(),*) 'Aborting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
       
    end select

   end subroutine InitPARTEHGlobals

   !==============================================================================================
    
    subroutine fates_history_maps
       
       use EDTypesMod, only : NFSC
       use EDTypesMod, only : nclmax
       use EDTypesMod, only : nlevleaf
       use EDParamsMod, only : ED_val_history_sizeclass_bin_edges
       use EDParamsMod, only : ED_val_history_ageclass_bin_edges
       use EDParamsMod, only : ED_val_history_height_bin_edges
       use EDParamsMod, only : ED_val_history_coageclass_bin_edges

       ! ------------------------------------------------------------------------------------------
       ! This subroutine allocates and populates the variables
       ! that define the mapping of variables in history files in multiplexed dimensions like
       ! the "scpf" format
       ! back to
       ! their respective single component dimensions, like size-class "sc" and pft "pf"
       ! ------------------------------------------------------------------------------------------

       integer :: i
       integer :: isc
       integer :: ipft
       integer :: icwd
       integer :: ifuel
       integer :: ican
       integer :: icdam
       integer :: ileaf
       integer :: iage
       integer :: iheight
       integer :: icoage
       integer :: iel

       allocate( fates_hdim_levsclass(1:nlevsclass   ))
       allocate( fates_hdim_pfmap_levscpf(1:nlevsclass*numpft))
       allocate( fates_hdim_scmap_levscpf(1:nlevsclass*numpft))
       allocate( fates_hdim_levpft(1:numpft   ))
       allocate( fates_hdim_levfuel(1:NFSC   ))
       allocate( fates_hdim_levcwdsc(1:NCWD   ))
       allocate( fates_hdim_levage(1:nlevage   ))
       allocate( fates_hdim_levheight(1:nlevheight   ))
       allocate( fates_hdim_levcoage(1:nlevcoage ))
       allocate( fates_hdim_pfmap_levcapf(1:nlevcoage*numpft))
       allocate( fates_hdim_camap_levcapf(1:nlevcoage*numpft))

       allocate( fates_hdim_levcan(nclmax))
       allocate( fates_hdim_levelem(num_elements))
       allocate( fates_hdim_levleaf(nlevleaf))
       allocate( fates_hdim_canmap_levcnlf(nlevleaf*nclmax))
       allocate( fates_hdim_lfmap_levcnlf(nlevleaf*nclmax))
       allocate( fates_hdim_canmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_lfmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_pftmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_scmap_levscag(nlevsclass * nlevage ))
       allocate( fates_hdim_agmap_levscag(nlevsclass * nlevage ))
       allocate( fates_hdim_scmap_levscagpft(nlevsclass * nlevage * numpft))
       allocate( fates_hdim_agmap_levscagpft(nlevsclass * nlevage * numpft))
       allocate( fates_hdim_pftmap_levscagpft(nlevsclass * nlevage * numpft))
       allocate( fates_hdim_agmap_levagepft(nlevage * numpft))
       allocate( fates_hdim_pftmap_levagepft(nlevage * numpft))
       allocate( fates_hdim_agmap_levagefuel(nlevage * nfsc))
       allocate( fates_hdim_fscmap_levagefuel(nlevage * nfsc))

       allocate( fates_hdim_elmap_levelpft(num_elements*numpft))
       allocate( fates_hdim_elmap_levelcwd(num_elements*ncwd))
       allocate( fates_hdim_elmap_levelage(num_elements*nlevage))
       allocate( fates_hdim_pftmap_levelpft(num_elements*numpft))
       allocate( fates_hdim_cwdmap_levelcwd(num_elements*ncwd))
       allocate( fates_hdim_agemap_levelage(num_elements*nlevage))

       allocate( fates_hdim_levdamage(1:nlevdamage ))
       allocate( fates_hdim_scmap_levcdsc(nlevsclass*nlevdamage))
       allocate( fates_hdim_cdmap_levcdsc(nlevsclass*nlevdamage))
       allocate( fates_hdim_scmap_levcdpf(nlevsclass*nlevdamage * numpft))
       allocate( fates_hdim_cdmap_levcdpf(nlevsclass*nlevdamage * numpft))
       allocate( fates_hdim_pftmap_levcdpf(nlevsclass*nlevdamage * numpft))
       
       ! Fill the IO array of plant size classes
       fates_hdim_levsclass(:) = ED_val_history_sizeclass_bin_edges(:)
       fates_hdim_levage(:) = ED_val_history_ageclass_bin_edges(:)
       fates_hdim_levheight(:) = ED_val_history_height_bin_edges(:)
       fates_hdim_levcoage(:) = ED_val_history_coageclass_bin_edges(:)
       fates_hdim_levleaf(:) = dlower_vai(:)
       fates_hdim_levdamage(:) = ED_val_history_damage_bin_edges(:)
       
       ! make pft array
       do ipft=1,numpft
          fates_hdim_levpft(ipft) = ipft
       end do

       ! make fuel array
       do ifuel=1,NFSC
          fates_hdim_levfuel(ifuel) = ifuel
       end do

       ! make cwd array
       do icwd=1,NCWD
          fates_hdim_levcwdsc(icwd) = icwd
       end do

       ! make canopy array
       do ican = 1,nclmax
          fates_hdim_levcan(ican) = ican
       end do

       ! Make an element array, each index is the PARTEH global identifier index

       do iel = 1, num_elements
           fates_hdim_levelem(iel) = element_list(iel)
       end do
       
       i = 0
       do iel = 1, num_elements
           do ipft=1,numpft
               i = i+1
               fates_hdim_elmap_levelpft(i)  = iel
               fates_hdim_pftmap_levelpft(i) = ipft
           end do
       end do
       
       i = 0
       do iel = 1, num_elements
           do icwd = 1, ncwd
               i = i+1
               fates_hdim_elmap_levelcwd(i)  = iel
               fates_hdim_cwdmap_levelcwd(i) = icwd
           end do
       end do
       
       i = 0
       do iel = 1, num_elements
           do iage=1,nlevage
               i = i+1
               fates_hdim_elmap_levelage(i) = iel
               fates_hdim_agemap_levelage(i) = iage
           end do
       end do

       ! Fill the IO arrays that match pft and size class to their combined array
       i=0
       do ipft=1,numpft
          do isc=1,nlevsclass
             i=i+1
             fates_hdim_pfmap_levscpf(i) = ipft
             fates_hdim_scmap_levscpf(i) = isc
          end do
       end do

       i=0
       do ipft=1,numpft
          do icoage=1,nlevcoage
             i=i+1
             fates_hdim_pfmap_levcapf(i) = ipft
             fates_hdim_camap_levcapf(i) = icoage
          end do
       end do

       i=0
       do ican=1,nclmax
          do ileaf=1,nlevleaf
             i=i+1
             fates_hdim_canmap_levcnlf(i) = ican
             fates_hdim_lfmap_levcnlf(i) = ileaf
          end do
       end do

       i=0
       do iage=1,nlevage
          do isc=1,nlevsclass
             i=i+1
             fates_hdim_scmap_levscag(i) = isc
             fates_hdim_agmap_levscag(i) = iage
          end do
       end do

       i=0
       do icdam=1,nlevdamage
          do isc=1,nlevsclass
             i=i+1
             fates_hdim_scmap_levcdsc(i) = isc
             fates_hdim_cdmap_levcdsc(i) = icdam
          end do
       end do

       i=0
       do ipft=1,numpft
          do icdam=1,nlevdamage
             do isc=1,nlevsclass
                i=i+1
                fates_hdim_scmap_levcdpf(i) = isc
                fates_hdim_cdmap_levcdpf(i) = icdam
                fates_hdim_pftmap_levcdpf(i) = ipft
             end do
          end do
       end do
       
       i=0
       do ipft=1,numpft
          do ican=1,nclmax
             do ileaf=1,nlevleaf
                i=i+1
                fates_hdim_canmap_levcnlfpf(i) = ican
                fates_hdim_lfmap_levcnlfpf(i) = ileaf
                fates_hdim_pftmap_levcnlfpf(i) = ipft
             end do
          end do
       end do

       i=0
       do ipft=1,numpft
          do iage=1,nlevage
             do isc=1,nlevsclass
                i=i+1
                fates_hdim_scmap_levscagpft(i) = isc
                fates_hdim_agmap_levscagpft(i) = iage
                fates_hdim_pftmap_levscagpft(i) = ipft
             end do
          end do
       end do

       i=0
       do ipft=1,numpft
          do iage=1,nlevage
             i=i+1
             fates_hdim_agmap_levagepft(i) = iage
             fates_hdim_pftmap_levagepft(i) = ipft
          end do
       end do

       i=0
       do iage=1,nlevage
          do ifuel=1,NFSC
             i=i+1
             fates_hdim_agmap_levagefuel(i) = iage
             fates_hdim_fscmap_levagefuel(i) = ifuel
          end do
       end do

    end subroutine fates_history_maps

    ! ===================================================================================

    subroutine SetFatesTime(current_year_in, current_month_in, &
                          current_day_in, current_tod_in, &
                          current_date_in, reference_date_in, &
                          model_day_in, day_of_year_in, &
                          days_per_year_in, freq_day_in)

     ! This subroutine should be called directly from the HLM
     
     integer,  intent(in) :: current_year_in
     integer,  intent(in) :: current_month_in
     integer,  intent(in) :: current_day_in
     integer,  intent(in) :: current_tod_in
     integer,  intent(in) :: current_date_in
     integer,  intent(in) :: reference_date_in
     real(r8), intent(in) :: model_day_in
     integer,  intent(in) :: day_of_year_in
     integer,  intent(in) :: days_per_year_in
     real(r8), intent(in) :: freq_day_in

     hlm_current_year   = current_year_in
     hlm_current_month  = current_month_in
     hlm_current_day    = current_day_in
     hlm_current_tod    = current_tod_in
     hlm_current_date   = current_date_in
     hlm_reference_date = reference_date_in
     hlm_model_day      = model_day_in
     hlm_day_of_year    = day_of_year_in
     hlm_days_per_year  = days_per_year_in
     hlm_freq_day       = freq_day_in

  end subroutine SetFatesTime

  ! ==================================================================================== 

  subroutine set_fates_ctrlparms(tag,ival,rval,cval)
      
      ! ---------------------------------------------------------------------------------
      ! Certain model control parameters and dimensions used by FATES are dictated by 
      ! the the driver or the host mode. To see which parameters should be filled here
      ! please also look at the ctrl_parms_type in FATESTYpeMod, in the section listing
      ! components dictated by the host model.
      !
      ! Some important points:
      ! 1. Calls to this function are likely from the clm_fates module in the HLM.
      ! 2. The calls should be preceeded by a flush function.
      ! 3. All values in ctrl_parm (FATESTypesMod.F90) that are classified as 
      !    'dictated by the HLM' must be listed in this subroutine
      ! 4. Should look like this:
      ! 
      ! call set_fates_ctrlparms('flush_to_unset')
      ! call set_fates_ctrlparms('num_sw_bbands',numrad)  ! or other variable
      ! ...
      ! call set_fates_ctrlparms('num_lev_soil',nlevsoi)   ! or other variable
      ! call set_fates_ctrlparms('check_allset') 
      !
      ! RGK-2016
      ! ---------------------------------------------------------------------------------
      use FatesConstantsMod, only : fates_check_param_set
    
    
      ! Arguments
      integer, optional, intent(in)         :: ival
      real(r8), optional, intent(in)        :: rval
      character(len=*),optional, intent(in) :: cval
      character(len=*),intent(in)           :: tag
      
      ! local variables
      logical              :: all_set
      integer,  parameter  :: unset_int = -999
      real(r8), parameter  :: unset_double = -999.9
      
      
      select case (trim(tag))
      case('flush_to_unset')
         if (fates_global_verbose()) then
            write(fates_log(), *) 'Flushing FATES control parameters prior to transfer from host'
         end if

         hlm_numSWb     = unset_int
         hlm_inir       = unset_int
         hlm_ivis       = unset_int
         hlm_is_restart = unset_int
         hlm_maxlevsoil   = unset_int
         hlm_name         = 'unset'
         hlm_hio_ignore_val   = unset_double
         hlm_masterproc   = unset_int
         hlm_ipedof       = unset_int
         hlm_nu_com      = 'unset'
         hlm_decomp      = 'unset'
         hlm_nitrogen_spec = unset_int
         hlm_use_tree_damage = unset_int
         hlm_phosphorus_spec = unset_int
         hlm_use_ch4       = unset_int
         hlm_use_vertsoilc = unset_int
         hlm_parteh_mode   = unset_int
         hlm_spitfire_mode = unset_int
         hlm_sf_nofire_def = unset_int
         hlm_sf_scalar_lightning_def = unset_int
         hlm_sf_successful_ignitions_def = unset_int
         hlm_sf_anthro_ignitions_def = unset_int
         hlm_use_planthydro = unset_int
         hlm_use_lu_harvest   = unset_int
         hlm_num_lu_harvest_cats   = unset_int
         hlm_use_cohort_age_tracking = unset_int
         hlm_use_logging   = unset_int
         hlm_use_ed_st3    = unset_int
         hlm_use_ed_prescribed_phys = unset_int
         hlm_use_fixed_biogeog = unset_int
         hlm_use_nocomp = unset_int   
         hlm_use_sp = unset_int
         hlm_use_inventory_init = unset_int
         hlm_inventory_ctrl_file = 'unset'

      case('check_allset')
         
         if(hlm_numSWb .eq. unset_int) then
            write(fates_log(), *) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_masterproc .eq. unset_int) then
            write(fates_log(), *) 'FATES parameter unset: hlm_masterproc'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numSWb > maxSWb) then
            write(fates_log(), *) 'FATES sets a maximum number of shortwave bands'
            write(fates_log(), *) 'for some scratch-space, maxSWb'
            write(fates_log(), *) 'it defaults to 2, but can be increased as needed'
            write(fates_log(), *) 'your driver or host model is intending to drive'
            write(fates_log(), *) 'FATES with:',hlm_numSWb,' bands.'
            write(fates_log(), *) 'please increase maxSWb in EDTypes to match'
            write(fates_log(), *) 'or exceed this value'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (  .not.((hlm_use_planthydro.eq.1).or.(hlm_use_planthydro.eq.0))    ) then
            write(fates_log(), *) 'The FATES namelist planthydro flag must be 0 or 1, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         elseif (hlm_use_planthydro.eq.1 ) then
               write(fates_log(), *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(fates_log(), *) ''
               write(fates_log(), *) ' use_fates_planthydro is an      EXPERIMENTAL FEATURE        '
               write(fates_log(), *) ' please see header of fates/biogeophys/FatesHydraulicsMod.F90'
               write(fates_log(), *) ' for more information.'
               write(fates_log(), *) ''
               write(fates_log(), *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if

         if ( (hlm_use_lu_harvest .lt. 0).or.(hlm_use_lu_harvest .gt. 1) ) then
            write(fates_log(), *) 'The FATES lu_harvest flag must be 0 or 1,  exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( (hlm_num_lu_harvest_cats .lt. 0) ) then
            write(fates_log(), *) 'The FATES number of hlm harvest cats must be >= 0, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( .not.((hlm_use_logging .eq.1).or.(hlm_use_logging.eq.0))    ) then
            write(fates_log(), *) 'The FATES namelist use_logging flag must be 0 or 1, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         if ( ( ANY(EDPftvarcon_inst%mort_ip_age_senescence < fates_check_param_set )) .and. &
           (hlm_use_cohort_age_tracking .eq.0 ) ) then
           write(fates_log(),*) 'Age dependent mortality cannot be on if'
           write(fates_log(),*) 'cohort age tracking is off.'
           write(fates_log(),*) 'Set use_fates_cohort_age_tracking = .true.'
           write(fates_log(),*) 'in FATES namelist options'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

         if (  .not.((hlm_use_ed_st3.eq.1).or.(hlm_use_ed_st3.eq.0))    ) then
            write(fates_log(), *) 'The FATES namelist stand structure flag must be 0 or 1, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if (  .not.((hlm_use_ed_prescribed_phys.eq.1).or.(hlm_use_ed_prescribed_phys.eq.0))    ) then
            write(fates_log(), *) 'The FATES namelist prescribed physiology flag must be 0 or 1, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( hlm_use_ed_prescribed_phys.eq.1 .and. hlm_use_ed_st3.eq.1 ) then
            write(fates_log(), *) 'FATES ST3 and prescribed physiology cannot both be turned on.'
            write(fates_log(), *) 'Review the namelist entries, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( hlm_use_inventory_init.eq.1  .and. hlm_use_cohort_age_tracking .eq.1) then
            write(fates_log(), *) 'Fates inventory init cannot be used with age dependent mortality'
            write(fates_log(), *) 'Set use_fates_cohort_age_tracking to 0 or turn off inventory init'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (  .not.((hlm_use_inventory_init.eq.1).or.(hlm_use_inventory_init.eq.0))    ) then
            write(fates_log(), *) 'The FATES NL inventory flag must be 0 or 1, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(trim(hlm_inventory_ctrl_file) .eq. 'unset') then
            write(fates_log(),*) 'namelist entry for fates inventory control file is unset, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ivis .ne. ivis) then
            write(fates_log(), *) 'FATES assumption about the index of visible shortwave'
            write(fates_log(), *) 'radiation is different from the HLM, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(hlm_inir .ne. inir) then
            write(fates_log(), *) 'FATES assumption about the index of NIR shortwave'
            write(fates_log(), *) 'radiation is different from the HLM, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_is_restart .eq. unset_int) then
            write(fates_log(), *) 'FATES parameter unset: hlm_is_restart, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_maxlevsoil .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: hlm_maxlevsoil, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_name) .eq. 'unset') then
            write(fates_log(),*) 'FATES dimension/parameter unset: hlm_name, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_decomp) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_decomp, exiting'
               write(fates_log(),*) 'valid: MIMICS, CENTURY, CTC'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if( .not. ((trim(hlm_decomp) .eq. 'MIMICS') .or. &
              (trim(hlm_decomp) .eq. 'CENTURY') .or. &
              (trim(hlm_decomp) .eq. 'CTC') .or. &
              (trim(hlm_decomp) .eq. 'NONE'))   ) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_decomp, exiting'
               write(fates_log(),*) 'valid: NONE, MIMICS, CENTURY, CTC, yours: ',trim(hlm_decomp)
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_nu_com) .eq. 'unset') then
            write(fates_log(),*) 'FATES dimension/parameter unset: hlm_nu_com, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_tree_damage .eq. unset_int) then
            write(fates_log(),*) 'FATES dimension/parameter unset: hlm_use_tree_damage, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         else
            if((hlm_use_tree_damage .eq. itrue) .and. &
                 (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp))then
               write(fates_log(),*) 'FATES tree damage (use_fates_tree_damage = .true.) is not'
               write(fates_log(),*) '(yet) compatible with CNP allocation (fates_parteh_mode = 2)'
               call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

            
         end if

         if(hlm_nitrogen_spec .eq. unset_int) then
            write(fates_log(),*) 'FATES parameters unset: hlm_nitrogen_spec, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_phosphorus_spec .eq. unset_int) then
            write(fates_log(),*) 'FATES parameters unset: hlm_phosphorus_spec, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if( abs(hlm_hio_ignore_val-unset_double)<1e-10 ) then
            write(fates_log(),*) 'FATES dimension/parameter unset: hio_ignore'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ipedof .eq. unset_int) then
            write(fates_log(), *) 'index for the HLMs pedotransfer function unset: hlm_ipedof, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_parteh_mode .eq. unset_int) then
            write(fates_log(), *) 'switch deciding which plant reactive transport model to use is unset, hlm_parteh_mode, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_ch4 .eq. unset_int) then
            write(fates_log(), *) 'switch for the HLMs CH4 module unset: hlm_use_ch4, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_vertsoilc .eq. unset_int) then
            write(fates_log(), *) 'switch for the HLMs soil carbon discretization unset: hlm_use_vertsoilc, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_spitfire_mode .eq. unset_int) then
            write(fates_log(), *) 'switch for SPITFIRE unset: hlm_spitfire_mode, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_nofire_def .eq. unset_int) then
            write(fates_log(), *) 'definition of no-fire mode unset: hlm_sf_nofire_def, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_scalar_lightning_def .eq. unset_int) then
            write(fates_log(), *) 'definition of scalar lightning mode unset: hlm_sf_scalltng_def, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_successful_ignitions_def .eq. unset_int) then
            write(fates_log(), *) 'definition of successful ignition mode unset: hlm_sf_successful, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_anthro_ignitions_def .eq. unset_int) then
            write(fates_log(), *) 'definition of anthro-ignition mode unset: hlm_sf_anthig_def, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_name).eq.'CLM' .and. hlm_parteh_mode .eq. 2) then
            if( sum(abs(EDPftvarcon_inst%prescribed_puptake(:)))<nearzero .and. &
                sum(abs(EDPftvarcon_inst%prescribed_nuptake(:)))<nearzero) then
               write(fates_log(), *) 'PARTEH hypothesis 2 is only viable with forced'
               write(fates_log(), *) 'boundary conditions for CLM (currently).'
               write(fates_log(), *) 'prescribed_puptake or prescribed_nuptake must > 0'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end if
         
        if(hlm_use_fixed_biogeog.eq.unset_int) then
           if(fates_global_verbose()) then
             write(fates_log(), *) 'switch for fixed biogeog unset: hlm_use_fixed_biogeog, exiting'
           end if
           call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_nocomp.eq.unset_int) then
            write(fates_log(), *) 'switch for no competition mode. '
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_sp.eq.unset_int) then
            write(fates_log(), *) 'switch for SP mode. '
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_cohort_age_tracking .eq. unset_int) then
            write(fates_log(), *) 'switch for cohort_age_tracking  unset: hlm_use_cohort_age_tracking, exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_sp.eq.itrue.and.hlm_use_nocomp.eq.ifalse)then
            write(fates_log(), *) 'SP cannot be on if nocomp mode is off. Exiting. '
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_sp.eq.itrue.and.hlm_use_fixed_biogeog.eq.ifalse)then
            write(fates_log(), *) 'SP cannot be on if fixed biogeog mode is off. Exiting. '
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (fates_global_verbose()) then
            write(fates_log(), *) 'Checked. All control parameters sent to FATES.'
         end if
         
      case default

         if(present(ival))then
            select case (trim(tag))

            case('masterproc')
               hlm_masterproc = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering masterproc = ',ival,' to FATES'
               end if

            case('num_sw_bbands')
               hlm_numSwb = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_sw_bbands = ',ival,' to FATES'
               end if
               
            case('vis_sw_index')
               hlm_ivis = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering index associated with visible SW rad = ',ival,' to FATES'
               end if
            
            case('nir_sw_index')
               hlm_inir = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering index associated with NIR SW rad = ',ival,' to FATES'
               end if

            case('is_restart')
               hlm_is_restart = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering flag signaling restart / not-restart = ',ival,' to FATES'
               end if

            case('num_lev_soil')
               hlm_maxlevsoil = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_soil = ',ival,' to FATES'
               end if

            case('soilwater_ipedof')
               hlm_ipedof = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_ipedof = ',ival,' to FATES'
               end if

            case('use_tree_damage')
               hlm_use_tree_damage = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_tree_damage = ',ival,' to FATES'
               end if
               
            case('nitrogen_spec')
               hlm_nitrogen_spec = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_nitrogen_spec = ',ival,' to FATES'
               end if

            case('phosphorus_spec')
               hlm_phosphorus_spec = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_phosphorus_spec = ',ival,' to FATES'
               end if

            case('use_ch4')
               hlm_use_ch4 = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_ch4 = ',ival,' to FATES'
               end if
               
            case('use_vertsoilc')
               hlm_use_vertsoilc = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_vertsoilc= ',ival,' to FATES'
               end if
               
            case('parteh_mode')
               hlm_parteh_mode = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_parteh_mode= ',ival,' to FATES'
               end if

            case('spitfire_mode')
               hlm_spitfire_mode = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_spitfire_mode =',ival,' to FATES'
              end if
              
           case('sf_nofire_def')
               hlm_sf_nofire_def = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_sf_nofire_def =',ival,' to FATES'
               end if

           case('sf_scalar_lightning_def')
               hlm_sf_scalar_lightning_def = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_sf_scalar_lightning_def =',ival,' to FATES'
               end if

           case('sf_successful_ignitions_def')
               hlm_sf_successful_ignitions_def = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_sf_successful_ignition_def =',ival,' to FATES'
               end if

           case('sf_anthro_ignitions_def')
               hlm_sf_anthro_ignitions_def = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_sf_anthro_ignition_def =',ival,' to FATES'
               end if

               
            case('use_fixed_biogeog')
                hlm_use_fixed_biogeog = ival
               if (fates_global_verbose()) then
                   write(fates_log(),*) 'Transfering hlm_use_fixed_biogeog= ',ival,' to FATES'
               end if
            
            case('use_nocomp')
                hlm_use_nocomp = ival
               if (fates_global_verbose()) then
                   write(fates_log(),*) 'Transfering hlm_use_nocomp= ',ival,' to FATES'
               end if

            case('use_sp')
               hlm_use_sp = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_sp= ',ival,' to FATES'
               end if

            case('use_planthydro')
               hlm_use_planthydro = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_planthydro= ',ival,' to FATES'
               end if

            case('use_lu_harvest')
               hlm_use_lu_harvest = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_lu_harvest= ',ival,' to FATES'
               end if

            case('num_lu_harvest_cats')
               hlm_num_lu_harvest_cats = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_num_lu_harvest_cats= ',ival,' to FATES'
               end if

            case('use_cohort_age_tracking')
               hlm_use_cohort_age_tracking = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_cohort_age_tracking= ',ival,' to FATES'
               end if

            case('use_logging')
               hlm_use_logging = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_logging= ',ival,' to FATES'
               end if

            case('use_ed_st3')
               hlm_use_ed_st3 = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_ed_st3= ',ival,' to FATES'
               end if

            case('use_ed_prescribed_phys')
               hlm_use_ed_prescribed_phys = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_ed_prescribed_phys= ',ival,' to FATES'
               end if

            case('use_inventory_init')
               hlm_use_inventory_init = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_inventory_init= ',ival,' to FATES'
               end if

            case default
               write(fates_log(), *) 'fates NL tag not recognized:',trim(tag)
               !! call endrun(msg=errMsg(sourcefile, __LINE__))
            end select
            
         end if
         
         if(present(rval))then
            select case (trim(tag))
            case ('hio_ignore_val')
               hlm_hio_ignore_val = rval
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hio_ignore_val = ',rval,' to FATES'
               end if
            case default
               write(fates_log(),*) 'fates NL tag not recognized:',trim(tag)
               !! call endrun(msg=errMsg(sourcefile, __LINE__))
            end select
         end if

         if(present(cval))then
            select case (trim(tag))
               
            case('hlm_name')
               hlm_name = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the HLM name = ',trim(cval)
               end if

            case('nu_com')
               hlm_nu_com = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the nutrient competition name = ',trim(cval)
               end if

            case('decomp_method')
               hlm_decomp = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the decomp method name = ',trim(cval)
               end if

            case('inventory_ctrl_file')
               hlm_inventory_ctrl_file = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the name of the inventory control file = ',trim(cval)
               end if
               
            case default
               write(fates_log(),*) 'fates NL tag not recognized:',trim(tag)
               !! call endrun(msg=errMsg(sourcefile, __LINE__))
            end select
         end if

      end select
            
      return
   end subroutine set_fates_ctrlparms

   ! ====================================================================================

   subroutine FatesReportParameters(masterproc)
      
      ! -----------------------------------------------------
      ! Simple parameter reporting functions
      ! A debug like print flag is contained in each routine
      ! -----------------------------------------------------

      logical,intent(in) :: masterproc

      call FatesReportPFTParams(masterproc)
      call FatesReportParams(masterproc)
      call FatesCheckParams(masterproc)    ! Check general fates parameters
      call PRTDerivedParams()              ! Update PARTEH derived constants
      call PRTCheckParams(masterproc)      ! Check PARTEH parameters
      call SpitFireCheckParams(masterproc)
      

      
      return
   end subroutine FatesReportParameters

   ! =====================================================================================

   subroutine UpdateFatesRMeansTStep(sites,bc_in)

     ! In this routine, we update any FATES buffers where
     ! we calculate running means. It is assumed that this buffer is updated
     ! on the model time-step.

     type(ed_site_type), intent(inout) :: sites(:)
     type(bc_in_type), intent(in)      :: bc_in(:)
     
     type(ed_patch_type),  pointer :: cpatch
     type(ed_cohort_type), pointer :: ccohort
     integer :: s, ifp, io_si

     do s = 1,size(sites,dim=1)

        ifp=0
        cpatch => sites(s)%oldest_patch
        do while(associated(cpatch))
           if (cpatch%patchno .ne. 0) then
           ifp=ifp+1
           call cpatch%tveg24%UpdateRMean(bc_in(s)%t_veg_pa(ifp))
           call cpatch%tveg_lpa%UpdateRMean(bc_in(s)%t_veg_pa(ifp))
           call cpatch%tveg_longterm%UpdateRMean(bc_in(s)%t_veg_pa(ifp))

           
           !ccohort => cpatch%tallest
           !do while (associated(ccohort))
           !   call ccohort%tveg_lpa%UpdateRMean(bc_in(s)%t_veg_pa(ifp))
           !   ccohort => ccohort%shorter
           !end do
           end if
           
           cpatch => cpatch%younger
        enddo
     end do

     return
   end subroutine UpdateFatesRMeansTStep
      
 end module FatesInterfaceMod
