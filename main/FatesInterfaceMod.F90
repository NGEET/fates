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
   use EDTypesMod                , only : maxPatchesPerSite
   use EDTypesMod                , only : maxCohortsPerPatch
   use EDTypesMod                , only : maxSWb
   use EDTypesMod                , only : ivis
   use EDTypesMod                , only : inir
   use EDTypesMod                , only : nclmax
   use EDTypesMod                , only : nlevleaf
   use EDTypesMod                , only : maxpft
   use EDTypesMod                , only : do_fates_salinity
   use EDTypesMod                , only : numWaterMem
   use EDTypesMod                , only : numlevsoil_max
   use FatesConstantsMod         , only : r8 => fates_r8
   use FatesConstantsMod         , only : itrue,ifalse
   use FatesConstantsMod         , only : nearzero
   use FatesGlobals              , only : fates_global_verbose
   use FatesGlobals              , only : fates_log
   use FatesGlobals              , only : endrun => fates_endrun
   use FatesLitterMod            , only : ncwd
   use FatesLitterMod            , only : ndcmpy
   use EDPftvarcon               , only : FatesReportPFTParams
   use EDPftvarcon               , only : FatesCheckParams
   use EDPftvarcon               , only : EDPftvarcon_inst
   use SFParamsMod               , only : SpitFireCheckParams
   use EDParamsMod               , only : FatesReportParams
   use EDParamsMod               , only : bgc_soil_salinity
   use FatesPlantHydraulicsMod   , only : InitHydroGlobals
   use EDParamsMod               , only : ED_val_history_sizeclass_bin_edges
   use EDParamsMod               , only : ED_val_history_ageclass_bin_edges
   use EDParamsMod               , only : ED_val_history_height_bin_edges
   use EDParamsMod               , only : ED_val_history_coageclass_bin_edges
   use CLMFatesParamInterfaceMod , only : FatesReadParameters
   use EDTypesMod                , only : p_uptake_mode
   use EDTypesMod                , only : n_uptake_mode
   use FatesConstantsMod         , only : prescribed_p_uptake
   use FatesConstantsMod         , only : prescribed_n_uptake
   use FatesConstantsMod         , only : coupled_p_uptake
   use FatesConstantsMod         , only : coupled_n_uptake
   use FatesConstantsMod         , only : fates_np_comp_scaling
   use FatesConstantsMod         , only : cohort_np_comp_scaling
   use FatesConstantsMod         , only : pft_np_comp_scaling
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
   use PRTInitParamsFatesMod     , only : PRTCheckParams
   use PRTAllometricCarbonMod    , only : InitPRTGlobalAllometricCarbon
   use PRTAllometricCNPMod       , only : InitPRTGlobalAllometricCNP


   ! CIME Globals
   use shr_log_mod               , only : errMsg => shr_log_errMsg
   use shr_infnan_mod            , only : nan => shr_infnan_nan, assignment(=)

   ! Just use everything from FatesInterfaceTypesMod, this is
   ! its sister code
   use FatesInterfaceTypesMod


   implicit none

   private

   character(len=*), parameter :: sourcefile = &
        __FILE__

   ! Make public necessary subroutines and functions
   public :: FatesInterfaceInit
   public :: set_fates_ctrlparms
   public :: SetFatesTime
   public :: SetFatesGlobalElements
   public :: FatesReportParameters
   public :: allocate_bcin
   public :: allocate_bcout
   public :: allocate_bcpconst
   public :: set_bcpconst
   public :: zero_bcs
   public :: set_bcs

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
    
    allocate(bc_pconst%eca_km_nh4(numpft))
    allocate(bc_pconst%eca_vmax_nh4(numpft))
    allocate(bc_pconst%eca_km_no3(numpft))
    allocate(bc_pconst%eca_vmax_no3(numpft))      
    allocate(bc_pconst%eca_km_p(numpft))
    allocate(bc_pconst%eca_vmax_p(numpft))
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
    
    bc_pconst%eca_km_nh4(1:numpft)       = EDPftvarcon_inst%eca_km_nh4(1:numpft)
    bc_pconst%eca_vmax_nh4(1:numpft)     = EDPftvarcon_inst%eca_vmax_nh4(1:numpft)
    bc_pconst%eca_km_no3(1:numpft)       = EDPftvarcon_inst%eca_km_no3(1:numpft)
    bc_pconst%eca_vmax_no3(1:numpft)     = EDPftvarcon_inst%eca_vmax_no3(1:numpft)
    bc_pconst%eca_km_p(1:numpft)         = EDPftvarcon_inst%eca_km_p(1:numpft)
    bc_pconst%eca_vmax_p(1:numpft)       = EDPftvarcon_inst%eca_vmax_p(1:numpft)
    bc_pconst%eca_km_ptase(1:numpft)     = EDPftvarcon_inst%eca_km_ptase(1:numpft)
    bc_pconst%eca_vmax_ptase(1:numpft)   = EDPftvarcon_inst%eca_vmax_ptase(1:numpft)
    bc_pconst%eca_alpha_ptase(1:numpft)  = EDPftvarcon_inst%eca_alpha_ptase(1:numpft) 
    bc_pconst%eca_lambda_ptase(1:numpft) = EDPftvarcon_inst%eca_lambda_ptase(1:numpft)
    bc_pconst%eca_plant_escalar          = eca_plant_escalar
    if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
       bc_pconst%j_uptake(1:nlevdecomp)    = 1
    else
       do j=1,nlevdecomp
          bc_pconst%j_uptake(j) = j
       end do
    end if
    
    return
  end subroutine set_bcpconst

  ! ====================================================================================
   
  subroutine zero_bcs(fates,s)

    type(fates_interface_type), intent(inout) :: fates
    integer, intent(in) :: s
    
    ! Input boundaries
    
    fates%bc_in(s)%t_veg24_pa(:)  = 0.0_r8
    fates%bc_in(s)%precip24_pa(:) = 0.0_r8
    fates%bc_in(s)%relhumid24_pa(:) = 0.0_r8
    fates%bc_in(s)%wind24_pa(:)     = 0.0_r8

    fates%bc_in(s)%lightning24(:)      = 0.0_r8
    fates%bc_in(s)%pop_density(:)      = 0.0_r8
    fates%bc_in(s)%solad_parb(:,:)     = 0.0_r8
    fates%bc_in(s)%solai_parb(:,:)     = 0.0_r8
    fates%bc_in(s)%smp_sl(:)           = 0.0_r8
    fates%bc_in(s)%eff_porosity_sl(:)  = 0.0_r8
    fates%bc_in(s)%watsat_sl(:)        = 0.0_r8
    fates%bc_in(s)%tempk_sl(:)         = 0.0_r8
    fates%bc_in(s)%h2o_liqvol_sl(:)    = 0.0_r8
    fates%bc_in(s)%filter_vegzen_pa(:) = .false.
    fates%bc_in(s)%coszen_pa(:)        = 0.0_r8
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

    ! Fates -> BGC fragmentation mass fluxes
    select case(hlm_parteh_mode) 
    case(prt_carbon_allom_hyp)
       fates%bc_out(s)%litt_flux_cel_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lig_c_si(:) = 0._r8
       fates%bc_out(s)%litt_flux_lab_c_si(:) = 0._r8
    case(prt_cnp_flex_allom_hyp) 
       
       fates%bc_in(s)%plant_n_uptake_flux(:,:) = 0._r8
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
    
    fates%bc_out(s)%canopy_fraction_pa(:) = 0.0_r8
    fates%bc_out(s)%frac_veg_nosno_alb_pa(:) = 0.0_r8
    
    if (hlm_use_planthydro.eq.itrue) then
       fates%bc_out(s)%qflx_soil2root_sisl(:) = 0.0_r8
       fates%bc_out(s)%qflx_ro_sisl(:)        = 0.0_r8
    end if
    fates%bc_out(s)%plant_stored_h2o_si = 0.0_r8
    
    return
  end subroutine zero_bcs

  ! ===========================================================================

   subroutine allocate_bcin(bc_in, nlevsoil_in, nlevdecomp_in, num_lu_harvest_cats)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      integer,intent(in)              :: nlevsoil_in
      integer,intent(in)              :: nlevdecomp_in
      integer,intent(in)              :: num_lu_harvest_cats
      
      ! Allocate input boundaries

      bc_in%nlevsoil   = nlevsoil_in

      if(nlevsoil_in > numlevsoil_max) then
         write(fates_log(), *) 'The number of soil layers imposed by the host model'
         write(fates_log(), *) 'is larger than what we have allocated in our static'
         write(fates_log(), *) 'arrays. Please increase the size of numlevsoil_max'
         write(fates_log(), *) 'found in EDTypesMod.F90'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if( (nlevsoil_in*ndcmpy) > fates_maxElementsPerPatch .or. &
          (nlevsoil_in*ncwd) > fates_maxElementsPerPatch) then
          write(fates_log(), *) 'The restart files require that space is allocated'
          write(fates_log(), *) 'to accomodate the multi-dimensional patch arrays'
          write(fates_log(), *) 'that are nlevsoil*numpft and nlevsoil*ncwd'
          write(fates_log(), *) 'fates_maxElementsPerPatch = ',fates_maxElementsPerPatch
          write(fates_log(), *) 'nlevsoil = ',nlevsoil_in
          write(fates_log(), *) 'dcmpy = ',ndcmpy
          write(fates_log(), *) 'ncwd  = ',ncwd
          write(fates_log(), *) 'numpft*nlevsoil = ',nlevsoil_in*numpft
          write(fates_log(), *) 'ncwd*nlevsoil = ',ncwd * nlevsoil_in
          write(fates_log(), *) 'To increase max_elements, change numlevsoil_max'
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
         if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
            allocate(bc_in%plant_n_uptake_flux(max_comp_per_site,1))
            allocate(bc_in%plant_p_uptake_flux(max_comp_per_site,1))
         else
            allocate(bc_in%plant_n_uptake_flux(max_comp_per_site,bc_in%nlevdecomp))
            allocate(bc_in%plant_p_uptake_flux(max_comp_per_site,bc_in%nlevdecomp))
         end if
      else
         allocate(bc_in%plant_n_uptake_flux(1,1))
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
      allocate(bc_in%lightning24(maxPatchesPerSite))
      allocate(bc_in%pop_density(maxPatchesPerSite))

      ! Vegetation Dynamics
      allocate(bc_in%t_veg24_pa(maxPatchesPerSite))

      allocate(bc_in%wind24_pa(maxPatchesPerSite))
      allocate(bc_in%relhumid24_pa(maxPatchesPerSite))
      allocate(bc_in%precip24_pa(maxPatchesPerSite))
      
      ! Radiation
      allocate(bc_in%solad_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_in%solai_parb(maxPatchesPerSite,hlm_numSWb))
      
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
      allocate(bc_in%filter_photo_pa(maxPatchesPerSite))
      allocate(bc_in%dayl_factor_pa(maxPatchesPerSite))
      allocate(bc_in%esat_tv_pa(maxPatchesPerSite))
      allocate(bc_in%eair_pa(maxPatchesPerSite))
      allocate(bc_in%oair_pa(maxPatchesPerSite))
      allocate(bc_in%cair_pa(maxPatchesPerSite))
      allocate(bc_in%rb_pa(maxPatchesPerSite))
      allocate(bc_in%t_veg_pa(maxPatchesPerSite))
      allocate(bc_in%tgcm_pa(maxPatchesPerSite))
      allocate(bc_in%t_soisno_sl(nlevsoil_in))

      ! Canopy Radiation
      allocate(bc_in%filter_vegzen_pa(maxPatchesPerSite))
      allocate(bc_in%coszen_pa(maxPatchesPerSite))
      allocate(bc_in%albgr_dir_rb(hlm_numSWb))
      allocate(bc_in%albgr_dif_rb(hlm_numSWb))

      ! Plant-Hydro BC's
      if (hlm_use_planthydro.eq.itrue) then

         allocate(bc_in%qflx_transp_pa(maxPatchesPerSite))
         allocate(bc_in%swrad_net_pa(maxPatchesPerSite))
         allocate(bc_in%lwrad_net_pa(maxPatchesPerSite))
         
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

      allocate(bc_in%pft_areafrac(maxpft))

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
      allocate(bc_out%fsun_pa(maxPatchesPerSite))
      allocate(bc_out%laisun_pa(maxPatchesPerSite))
      allocate(bc_out%laisha_pa(maxPatchesPerSite))
      
      ! Hydrology
      allocate(bc_out%active_suction_sl(nlevsoil_in))
      allocate(bc_out%rootr_pasl(maxPatchesPerSite,nlevsoil_in))
      allocate(bc_out%btran_pa(maxPatchesPerSite))
      
      ! Photosynthesis

      allocate(bc_out%rssun_pa(maxPatchesPerSite))
      allocate(bc_out%rssha_pa(maxPatchesPerSite))
      
      ! Canopy Radiation
      allocate(bc_out%albd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%albi_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%fabd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%fabi_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftdd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftid_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftii_parb(maxPatchesPerSite,hlm_numSWb))


      ! We allocate the boundary conditions to the BGC
      ! model, regardless of what scheme we use. The BGC
      ! model in ELM allocates all species C,N,P even if they
      ! are not turned on. Also, it is feasible that the
      ! one would want to allow soil BGC nutrient dynamics
      ! to proceed even if we are not passing source fluxes
      ! or uptake from FATES.
      ! When FATES does not have nutrients enabled, these
      ! arrays are indexed by 1.
      
      if(trim(hlm_nu_com).eq.'RD') then
         allocate(bc_out%n_demand(max_comp_per_site))
         allocate(bc_out%p_demand(max_comp_per_site))
      end if

      if(trim(hlm_nu_com).eq.'ECA') then
         allocate(bc_out%veg_rootc(max_comp_per_site,nlevdecomp_in))
         allocate(bc_out%decompmicc(nlevdecomp_in))
         allocate(bc_out%ft_index(max_comp_per_site))
         allocate(bc_out%cn_scalar(max_comp_per_site))
         allocate(bc_out%cp_scalar(max_comp_per_site))
      end if

      
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
      allocate(bc_out%elai_pa(maxPatchesPerSite))
      allocate(bc_out%esai_pa(maxPatchesPerSite))
      allocate(bc_out%tlai_pa(maxPatchesPerSite))
      allocate(bc_out%tsai_pa(maxPatchesPerSite))
      allocate(bc_out%htop_pa(maxPatchesPerSite))
      allocate(bc_out%hbot_pa(maxPatchesPerSite))
      allocate(bc_out%dleaf_pa(maxPatchesPerSite))

      allocate(bc_out%displa_pa(maxPatchesPerSite))
      allocate(bc_out%z0m_pa(maxPatchesPerSite))

      allocate(bc_out%canopy_fraction_pa(maxPatchesPerSite))
      allocate(bc_out%frac_veg_nosno_alb_pa(maxPatchesPerSite))

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
    
    subroutine SetFatesGlobalElements(use_fates)

       ! --------------------------------------------------------------------------------
       !
       ! This is the first FATES routine that is called.
       !
       ! This subroutine MUST BE CALLED AFTER the FATES PFT parameter file has been read in,
       ! and the EDPftvarcon_inst structure has been made.
       ! This subroutine MUST BE CALLED AFTER NL VARIABLES ARE READ (ie hlm_parteh_mode,etc)
       ! This subroutine must ALSO BE CALLED BEFORE the history file dimensions
       ! are set.
       ! 
       ! This routine requires no information from the HLM. This routine is responsible
       ! for generating the globals that are required by the HLM that are entirely
       ! FATES derived.
       !
       ! --------------------------------------------------------------------------------


      implicit none
      
      logical,intent(in) :: use_fates    ! Is fates turned on?
      integer :: i
      
      if (use_fates) then
         
         ! first read the non-PFT parameters
         call FatesReadParameters()

         ! Identify the number of PFTs by evaluating a pft array
         ! Using wood density as that is not expected to be deprecated any time soon

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

         if ( hlm_use_cohort_age_tracking .eq. itrue) then
            maxCohortsPerPatch = 300
         else
            maxCohortsPerPatch = 100
         end if
         
         ! These values are used to define the restart file allocations and general structure
         ! of memory for the cohort arrays
         
         fates_maxElementsPerPatch = max(maxCohortsPerPatch, ndcmpy*numlevsoil_max ,ncwd*numlevsoil_max)

         if (maxPatchesPerSite * fates_maxElementsPerPatch <  numWaterMem) then
            write(fates_log(), *) 'By using such a tiny number of maximum patches and maximum cohorts'
            write(fates_log(), *) ' this could create problems for indexing in restart files'
            write(fates_log(), *) ' The multiple of the two has to be greater than numWaterMem'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         fates_maxElementsPerSite = maxPatchesPerSite * fates_maxElementsPerPatch


         ! Set the maximum number of nutrient aquisition competitors per site
         ! This is used to set array sizes for the boundary conditions.
         ! Note: since BGC code may be active even when no nutrients
         ! present, we still need to allocate things when no nutrients

         if (hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp ) then
            if(fates_np_comp_scaling.eq.cohort_np_comp_scaling) then
               max_comp_per_site = fates_maxElementsPerSite
            elseif(fates_np_comp_scaling.eq.pft_np_comp_scaling) then
               max_comp_per_site = numpft
            else
               write(fates_log(), *) 'An unknown nutrient competitor scaling method was chosen?'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         else
            max_comp_per_site = 1
         end if
            


         ! Identify number of size and age class bins for history output
         ! assume these arrays are 1-indexed
         nlevsclass = size(ED_val_history_sizeclass_bin_edges,dim=1)
         nlevage = size(ED_val_history_ageclass_bin_edges,dim=1)
         nlevheight = size(ED_val_history_height_bin_edges,dim=1)
         nlevcoage = size(ED_val_history_coageclass_bin_edges,dim=1)

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


    end subroutine SetFatesGlobalElements

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


       ! Fill the IO array of plant size classes
       fates_hdim_levsclass(:) = ED_val_history_sizeclass_bin_edges(:)
       fates_hdim_levage(:) = ED_val_history_ageclass_bin_edges(:)
       fates_hdim_levheight(:) = ED_val_history_height_bin_edges(:)
       fates_hdim_levcoage(:) = ED_val_history_coageclass_bin_edges(:)

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
      ! call set_fates_ctrlparms('num_lev_ground',nlevgrnd)   ! or other variable
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
         hlm_numlevgrnd   = unset_int
         hlm_name         = 'unset'
         hlm_hio_ignore_val   = unset_double
         hlm_masterproc   = unset_int
         hlm_ipedof       = unset_int
         hlm_nu_com      = 'unset'
         hlm_nitrogen_spec = unset_int
         hlm_phosphorus_spec = unset_int
         hlm_max_patch_per_site = unset_int
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
         !hlm_use_nocomp = unset_int    ! future reduced complexity mode
         hlm_use_inventory_init = unset_int
         hlm_inventory_ctrl_file = 'unset'

      case('check_allset')
         
         if(hlm_numSWb .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_masterproc .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES parameter unset: hlm_masterproc'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numSWb > maxSWb) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES sets a maximum number of shortwave bands'
               write(fates_log(), *) 'for some scratch-space, maxSWb'
               write(fates_log(), *) 'it defaults to 2, but can be increased as needed'
               write(fates_log(), *) 'your driver or host model is intending to drive'
               write(fates_log(), *) 'FATES with:',hlm_numSWb,' bands.'
               write(fates_log(), *) 'please increase maxSWb in EDTypes to match'
               write(fates_log(), *) 'or exceed this value'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (  .not.((hlm_use_planthydro.eq.1).or.(hlm_use_planthydro.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist planthydro flag must be 0 or 1, exiting'
            end if
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
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES lu_harvest flag must be 0 or 1,  exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( (hlm_num_lu_harvest_cats .lt. 0) ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES number of hlm harvest cats must be >= 0, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( .not.((hlm_use_logging .eq.1).or.(hlm_use_logging.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist use_logging flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         if ( ( ANY(EDPftvarcon_inst%mort_ip_age_senescence < fates_check_param_set )) .and. &
           (hlm_use_cohort_age_tracking .eq.0 ) ) then

           write(fates_log(),*) 'Age dependent mortality cannot be on if'
           write(fates_log(),*) 'cohort age tracking is off.'
           write(fates_log(),*) 'Set hlm_use_cohort_age_tracking = .true.'
           write(fates_log(),*) 'in FATES namelist options'
           write(fates_log(),*) 'Aborting'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
         

         if (  .not.((hlm_use_ed_st3.eq.1).or.(hlm_use_ed_st3.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist stand structure flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if (  .not.((hlm_use_ed_prescribed_phys.eq.1).or.(hlm_use_ed_prescribed_phys.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist prescribed physiology flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( hlm_use_ed_prescribed_phys.eq.1 .and. hlm_use_ed_st3.eq.1 ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES ST3 and prescribed physiology cannot both be turned on.'
               write(fates_log(), *) 'Review the namelist entries, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( hlm_use_inventory_init.eq.1  .and. hlm_use_cohort_age_tracking .eq.1) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'Fates inventory init cannot be used with age dependent mortality'
               write(fates_log(), *) 'Set hlm_use_cohort_age_tracking to 0 or turn off inventory init'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         

         
         if (  .not.((hlm_use_inventory_init.eq.1).or.(hlm_use_inventory_init.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES NL inventory flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(trim(hlm_inventory_ctrl_file) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'namelist entry for fates inventory control file is unset, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ivis .ne. ivis) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES assumption about the index of visible shortwave'
               write(fates_log(), *) 'radiation is different from the HLM, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(hlm_inir .ne. inir) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES assumption about the index of NIR shortwave'
               write(fates_log(), *) 'radiation is different from the HLM, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_is_restart .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES parameter unset: hlm_is_restart, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numlevgrnd .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_name) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_name, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_nu_com) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_nu_com, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(hlm_nitrogen_spec .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES parameters unset: hlm_nitrogen_spec, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_phosphorus_spec .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES parameters unset: hlm_phosphorus_spec, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if( abs(hlm_hio_ignore_val-unset_double)<1e-10 ) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hio_ignore'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ipedof .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'index for the HLMs pedotransfer function unset: hlm_ipedof, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         if(hlm_max_patch_per_site .eq. unset_int ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'the number of patch-space per site unset: hlm_max_patch_per_site, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         elseif(hlm_max_patch_per_site < maxPatchesPerSite ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES is trying to allocate space for more patches per site, than the HLM has space for.'
               write(fates_log(), *) 'hlm_max_patch_per_site (HLM side): ', hlm_max_patch_per_site
               write(fates_log(), *) 'maxPatchesPerSite (FATES side): ', maxPatchesPerSite
               write(fates_log(), *)
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_parteh_mode .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch deciding which plant reactive transport model to use is unset, hlm_parteh_mode, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_vertsoilc .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch for the HLMs soil carbon discretization unset: hlm_use_vertsoilc, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_spitfire_mode .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch for SPITFIRE unset: hlm_spitfire_mode, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_nofire_def .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'definition of no-fire mode unset: hlm_sf_nofire_def, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_scalar_lightning_def .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'definition of scalar lightning mode unset: hlm_sf_scalltng_def, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_successful_ignitions_def .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'definition of successful ignition mode unset: hlm_sf_successful, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if(hlm_sf_anthro_ignitions_def .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'definition of anthro-ignition mode unset: hlm_sf_anthig_def, exiting'
            end if
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
             write(fates_log(), *) 'switch for fixed biogeog unset: him_use_fixed_biogeog, exiting'
           end if
           call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

        ! Future reduced complexity mode   
        !if(hlm_use_nocomp.eq.unset_int) then
        !      if(fates_global_verbose()) then
        !     write(fates_log(), *) 'switch for no competition mode. '
        !    end if
        !   call endrun(msg=errMsg(sourcefile, __LINE__))
        ! end if

         if(hlm_use_cohort_age_tracking .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch for cohort_age_tracking  unset: hlm_use_cohort_age_tracking, exiting'
            end if
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

            case('num_lev_ground')
               hlm_numlevgrnd = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_ground = ',ival,' to FATES'
               end if

            case('soilwater_ipedof')
               hlm_ipedof = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_ipedof = ',ival,' to FATES'
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

               
            case('max_patch_per_site')
               hlm_max_patch_per_site = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_max_patch_per_site = ',ival,' to FATES'
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
            
            ! Future reduced complexity mode   
            !case('use_nocomp')
            !    hlm_use_nocomp = ival
            !   if (fates_global_verbose()) then
            !       write(fates_log(),*) 'Transfering hlm_use_nocomp= ',ival,' to FATES'
            !   end if


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
               if (fates_global_verbose()) then
                   write(fates_log(), *) 'tag not recognized:',trim(tag)
                   call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
               ! end_run
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
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
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

            case('inventory_ctrl_file')
               hlm_inventory_ctrl_file = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the name of the inventory control file = ',trim(cval)
               end if
               
            case default
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
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
      call PRTCheckParams(masterproc)      ! Check PARTEH parameters
      call SpitFireCheckParams(masterproc)

      return
   end subroutine FatesReportParameters

  ! =====================================================================================




end module FatesInterfaceMod
