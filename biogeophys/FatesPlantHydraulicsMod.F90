module FatesPlantHydraulicsMod

  ! ==============================================================================================
  ! This module contains the relevant code for plant hydraulics. Currently, one hydraulics module
  ! is available.  Other methods of estimating plant hydraulics may become available in future
  ! releases.  For now, please cite the following reference if this module is used to generate
  ! published research:
  !     
  ! Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N. M., Galbraith, D. R., Baker, 
  !   T. R., Kruijt, B., Rowland, L., Fisher, R. A., Binks, O. J., Sevanto, S., Xu, C., Jansen, 
  !   S., Choat, B., Mencuccini, M., McDowell, N. G., Meir, P. Linking hydraulic traits to 
  !   tropical forest function in a size-structured and trait-driven model (TFS~v.1-Hydro). 
  !   Geoscientific Model Development, 9(11), 2016, pp: 4227-4255, 
  !   https://www.geosci-model-dev.net/9/4227/2016/, DOI = 10.5194/gmd-9-4227-2016.
  !
  ! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
  !
  !  PLANT HYDRAULICS IS AN EXPERIMENTAL OPTION THAT IS STILL UNDERGOING TESTING.  
  ! 
  ! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
  !
  ! ==============================================================================================

  use FatesGlobals, only      : endrun => fates_endrun
  use FatesGlobals, only      : fates_log

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_huge
  use FatesConstantsMod, only : denh2o => dens_fresh_liquid_water
  use FatesConstantsMod, only : grav_earth
  use FatesConstantsMod, only : ifalse, itrue
  use FatesConstantsMod, only : pi_const
  use FatesConstantsMod, only : cm2_per_m2
  use FatesConstantsMod, only : g_per_kg
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : mpa_per_pa
  use FatesConstantsMod, only : m_per_mm
  use FatesConstantsMod, only : mg_per_kg
  use FatesConstantsMod, only : pa_per_mpa
  use FatesConstantsMod, only : rsnbl_math_prec
  use FatesConstantsMod, only : m3_per_mm3
  use FatesConstantsMod, only : cm3_per_m3
  use FatesConstantsMod, only : kg_per_g
  use FatesConstantsMod, only : fates_unset_r8

  use EDParamsMod       , only : hydr_kmax_rsurf1
  use EDParamsMod       , only : hydr_kmax_rsurf2
  use EDParamsMod       , only : hydr_psi0
  use EDParamsMod       , only : hydr_psicap
  
  use EDTypesMod        , only : ed_site_type
  use EDTypesMod        , only : ed_patch_type
  use EDTypesMod        , only : ed_cohort_type
  use EDTypesMod        , only : AREA_INV
  use EDTypesMod        , only : AREA
  use EDTypesMod        , only : leaves_on

  use FatesInterfaceTypesMod  , only : bc_in_type
  use FatesInterfaceTypesMod  , only : bc_out_type
  use FatesInterfaceTypesMod  , only : hlm_use_planthydro
  use FatesInterfaceTypesMod  , only : hlm_ipedof
  use FatesInterfaceTypesMod  , only : numpft
  use FatesInterfaceTypesMod  , only : nlevsclass

  use FatesAllometryMod, only    : bleaf
  use FatesAllometryMod, only    : bsap_allom
  use FatesAllometryMod, only    : CrownDepth
  use FatesAllometryMod , only   : set_root_fraction
  use FatesHydraulicsMemMod, only: use_2d_hydrosolve
  use FatesHydraulicsMemMod, only: ed_site_hydr_type
  use FatesHydraulicsMemMod, only: ed_cohort_hydr_type
  use FatesHydraulicsMemMod, only: n_hypool_plant
  use FatesHydraulicsMemMod, only: n_hypool_leaf
  use FatesHydraulicsMemMod, only: n_hypool_tot
  use FatesHydraulicsMemMod, only: n_hypool_stem
  use FatesHydraulicsMemMod, only: n_hypool_troot
  use FatesHydraulicsMemMod, only: n_hypool_aroot
  use FatesHydraulicsMemMod, only: n_plant_media
  use FatesHydraulicsMemMod, only: nshell
  use FatesHydraulicsMemMod, only: n_hypool_ag
  use FatesHydraulicsMemMod, only: stomata_p_media
  use FatesHydraulicsMemMod, only: leaf_p_media
  use FatesHydraulicsMemMod, only: stem_p_media
  use FatesHydraulicsMemMod, only: troot_p_media
  use FatesHydraulicsMemMod, only: aroot_p_media
  use FatesHydraulicsMemMod, only: rhiz_p_media
  use FatesHydraulicsMemMod, only: nlevsoi_hyd_max
  use FatesHydraulicsMemMod, only: cohort_recruit_water_layer
  use FatesHydraulicsMemMod, only: recruit_water_avail_layer  
  use FatesHydraulicsMemMod, only: rwccap, rwcft
  use FatesHydraulicsMemMod, only: ignore_layer1
  
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,          only : store_organ, repro_organ, struct_organ

  use clm_time_manager  , only : get_step_size, get_nstep

  use EDPftvarcon, only : EDPftvarcon_inst

  use FatesHydroWTFMod, only : wrf_arr_type
  use FatesHydroWTFMod, only : wkf_arr_type
  use FatesHydroWTFMod, only : wrf_type, wrf_type_vg, wrf_type_cch, wrf_type_tfs
  use FatesHydroWTFMod, only : wkf_type, wkf_type_vg, wkf_type_cch, wkf_type_tfs


  ! CIME Globals
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : isnan => shr_infnan_isnan


  implicit none


  ! 1=leaf, 2=stem, 3=troot, 4=aroot
  ! Several of these may be better transferred to the parameter file in due time (RGK)

  integer, public :: use_ed_planthydraulics    =  1      ! 0 => use vanilla btran
  ! 1 => use BC hydraulics; 
  ! 2 => use CX hydraulics

  ! The following options are temporarily unavailable (RGK 09-06-19)
  ! ----------------------------------------------------------------------------------

  ! logical, public :: do_dqtopdth_leaf          = .false.  ! should a nonzero dqtopdth_leaf
  ! term be applied to the plant 
  ! hydraulics numerical solution?
  ! logical, public :: do_dyn_xylemrefill        = .false.  ! should the dynamics of xylem refilling 
  ! (i.e., non-instantaneous) be considered 
  ! within plant hydraulics?
  ! logical, public :: do_kbound_upstream        = .true.  ! should the hydraulic conductance at the 
  ! boundary between nodes be taken to be a
  ! function of the upstream loss of 
  ! conductivity (flc)?

  ! DO NOT TURN THIS ON. LEAVING THIS ONLY IF THE HLMS START HAVING
  ! TROUBLE RESPONDING TO SUPERSATURATION
  logical :: purge_supersaturation = .false. ! If for some reason the roots force water
                                             ! into a saturated soil layer, or push it slightly
                                             ! past saturation, should we attempt to help
                                             ! fix the situation by assigning some
                                             ! of the water to a runoff term?
  

  logical, public :: do_growthrecruiteffects   = .true. ! should size- or root length-dependent 
                                                        ! hydraulic properties and states be 
                                                        ! updated every day when trees grow or 
                                                        ! when recruitment happens?

  ! If this is set to true, then the conductance over a path between nodes, is defined
  ! by the side of the path with higher potential only. 
  logical, parameter     :: do_upstream_k = .true.


  
  logical :: do_parallel_stem = .true. ! If this mode is active, we treat the conduit through
                                       ! the plant (in 1D solves) as closed from root layer
                                       ! to the stomata. The effect of this, is that
                                       ! the conductances through stem and leaf surfaces
                                       ! are reduced by the fraction of active root
                                       ! conductance, and for each soil-layer, integration
                                       ! proceeds over the entire time-step.



  real(r8), parameter :: thsat_buff = 0.001_r8 ! Ensure that this amount of buffer
                                               ! is left between soil moisture and saturation [m3/m3]
                                               ! (if we are going to help purge super-saturation)
                                             
  logical,parameter :: debug = .false.          ! flag to report warning in hydro


  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  integer, public, parameter :: van_genuchten_type      = 1
  integer, public, parameter :: campbell_type           = 2
  integer, public, parameter :: tfs_type                = 3
  
  integer, parameter :: plant_wrf_type = tfs_type
  integer, parameter :: plant_wkf_type = tfs_type
  integer, parameter :: soil_wrf_type  = campbell_type
  integer, parameter :: soil_wkf_type  = campbell_type
  
  
  ! Define the global object that holds the water retention functions
  ! for plants of each different porous media type, and plant functional type
  
  class(wrf_arr_type),pointer :: wrf_plant(:,:)
  
  ! Define the global object that holds the water conductance functions
  ! for plants of each different porous media type, and plant functional type
  
  class(wkf_arr_type), pointer :: wkf_plant(:,:)

  ! Testing parameters for Van Genuchten soil WRTs
  ! unused unless van_genuchten_type is selected, also
  ! it would be much better to use the native parameters passed in
  ! from the HLM's soil model
  real(r8), parameter :: alpha_vg  = 0.001_r8
  real(r8), parameter :: th_sat_vg = 0.65_r8
  real(r8), parameter :: th_res_vg = 0.15_r8
  real(r8), parameter :: psd_vg    = 2.7_r8
  real(r8), parameter :: tort_vg   = 0.5_r8

  ! The maximum allowable water balance error over a plant-soil continuum
  ! for a given step [kgs] (0.1 mg)
  real(r8), parameter :: max_wb_step_err = 1.e-7_r8 

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: AccumulateMortalityWaterStorage
  public :: RecruitWaterStorage
  public :: hydraulics_drive
  public :: InitHydrSites
  public :: HydrSiteColdStart
  public :: BTranForHLMDiagnosticsFromCohortHydr
  public :: InitHydrCohort
  public :: DeallocateHydrCohort
  public :: UpdateH2OVeg
  public :: CopyCohortHydraulics
  public :: FuseCohortHydraulics
  public :: UpdateSizeDepPlantHydProps
  public :: UpdateSizeDepPlantHydStates
  public :: UpdatePlantPsiFTCFromTheta
  public :: InitPlantHydStates
  public :: UpdateSizeDepRhizHydProps
  public :: UpdateSizeDepRhizHydStates
  public :: RestartHydrStates
  public :: SavePreviousCompartmentVolumes
  public :: SavePreviousRhizVolumes
  public :: UpdatePlantHydrNodes
  public :: UpdatePlantHydrLenVol
  public :: UpdatePlantKmax
  public :: ConstrainRecruitNumber
  public :: InitHydroGlobals

  !------------------------------------------------------------------------------
  ! 01/18/16: Created by Brad Christoffersen
  ! 02/xx/17: Refactoring by Ryan Knox and Brad Christoffersen
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine hydraulics_drive( nsites, sites, bc_in,bc_out,dtime )

    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    real(r8),intent(in)                     :: dtime


    select case (use_ed_planthydraulics)

    case (1)

       call FillDrainRhizShells(nsites, sites, bc_in, bc_out )
       call hydraulics_BC(nsites, sites,bc_in,bc_out,dtime )

    case (2)

       !call Hydraulics_CX()

    case DEFAULT

    end select

  end subroutine Hydraulics_Drive

  ! =====================================================================================

  subroutine RestartHydrStates(sites,nsites,bc_in,bc_out)

    ! It is assumed that the following state variables have been read in by
    ! the restart machinery. 
    !
    ! co_hydr%th_ag
    ! co_hydr%th_troot
    ! co_hydr%th_aroot
    ! si_hydr%r_node_shell
    ! si_hydr%v_shell
    ! si_hydr%h2osoi_liqvol_shell
    ! si_hydr%l_aroot_layer
    !
    ! The goal of this subroutine is to call
    ! the correct sequence of hydraulics initializations to repopulate
    ! information that relies on these key states, as well as other vegetation
    ! states such as carbon pools and plant geometry.

    integer                     , intent(in)            :: nsites
    type(ed_site_type)          , intent(inout), target :: sites(nsites)
    type(bc_in_type)            , intent(in)            :: bc_in(nsites)
    type(bc_out_type)           , intent(inout)         :: bc_out(nsites)

    ! locals
    ! ----------------------------------------------------------------------------------
    ! LL pointers
    type(ed_patch_type),pointer                         :: cpatch      ! current patch
    type(ed_cohort_type),pointer                        :: ccohort     ! current cohort
    type(ed_cohort_hydr_type),pointer                   :: ccohort_hydr
    type(ed_site_hydr_type),pointer                     :: csite_hydr   
    integer                                             :: s           ! site loop counter
    integer                                             :: j           ! soil layer index
    integer                                             :: j_bc        ! soil layer index of boundary condition
    class(wrf_type_vg), pointer                         :: wrf_vg
    class(wkf_type_vg), pointer                         :: wkf_vg
    class(wrf_type_cch), pointer                        :: wrf_cch
    class(wkf_type_cch), pointer                        :: wkf_cch

    do s = 1,nsites
       csite_hydr=>sites(s)%si_hydr
       
       cpatch => sites(s)%oldest_patch
       do while(associated(cpatch))

          ccohort => cpatch%shortest
          do while(associated(ccohort))  

             ccohort_hydr => ccohort%co_hydr

             ! This calculates node heights
             call UpdatePlantHydrNodes(ccohort_hydr,ccohort%pft,ccohort%hite, &
                  sites(s)%si_hydr)

             ! This calculates volumes and lengths
             call UpdatePlantHydrLenVol(ccohort,csite_hydr)

             ! This updates the Kmax's of the plant's compartments
             call UpdatePlantKmax(ccohort_hydr,ccohort,sites(s)%si_hydr)

             ! Since this is a newly initialized plant, we set the previous compartment-size
             ! equal to the ones we just calculated.
             call SavePreviousCompartmentVolumes(ccohort_hydr) 

             ccohort => ccohort%taller
          enddo

          cpatch => cpatch%younger
       end do

       sites(s)%si_hydr%l_aroot_layer_init(:)  = fates_unset_r8
       sites(s)%si_hydr%r_node_shell_init(:,:) = fates_unset_r8
       sites(s)%si_hydr%v_shell_init(:,:)      = fates_unset_r8

       ! --------------------------------------------------------------------------------
       ! Initialize water transfer functions
       ! which include both water retention functions (WRFs)
       ! as well as the water conductance (K) functions (WKFs)
       ! But, this is only for soil!
       ! --------------------------------------------------------------------------------
       ! Initialize the Water Retention Functions
       ! -----------------------------------------------------------------------------------
       
       select case(soil_wrf_type)
       case(van_genuchten_type)
          do j=1,sites(s)%si_hydr%nlevrhiz
             j_bc = j+csite_hydr%i_rhiz_t-1
             allocate(wrf_vg)
             sites(s)%si_hydr%wrf_soil(j)%p => wrf_vg
             call wrf_vg%set_wrf_param([alpha_vg, psd_vg, bc_in(s)%watsat_sisl(j_bc), th_res_vg])
          end do
       case(campbell_type)
          do j=1,sites(s)%si_hydr%nlevrhiz
             j_bc = j+csite_hydr%i_rhiz_t-1
             allocate(wrf_cch)
             sites(s)%si_hydr%wrf_soil(j)%p => wrf_cch
             call wrf_cch%set_wrf_param([bc_in(s)%watsat_sisl(j_bc), &  
                  (-1.0_r8)*bc_in(s)%sucsat_sisl(j_bc)*denh2o*grav_earth*mpa_per_pa*m_per_mm , & 
                                      bc_in(s)%bsw_sisl(j_bc)])
          end do
       case(tfs_type)
          write(fates_log(),*) 'TFS water retention curves not available for soil'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       
       ! -----------------------------------------------------------------------------------
       ! Initialize the Water Conductance (K) Functions
       ! -----------------------------------------------------------------------------------
       
       select case(soil_wkf_type)
       case(van_genuchten_type)
          do j=1,sites(s)%si_hydr%nlevrhiz
             j_bc = j+csite_hydr%i_rhiz_t-1
             allocate(wkf_vg)
             sites(s)%si_hydr%wkf_soil(j)%p => wkf_vg
             call wkf_vg%set_wkf_param([alpha_vg, psd_vg, th_sat_vg, th_res_vg, tort_vg])
          end do
       case(campbell_type)
          do j=1,sites(s)%si_hydr%nlevrhiz
             j_bc = j+csite_hydr%i_rhiz_t-1
             allocate(wkf_cch)
             sites(s)%si_hydr%wkf_soil(j)%p => wkf_cch
             call wkf_cch%set_wkf_param([bc_in(s)%watsat_sisl(j_bc), &
                  (-1.0_r8)*bc_in(s)%sucsat_sisl(j_bc)*denh2o*grav_earth*mpa_per_pa*m_per_mm , & 
                  bc_in(s)%bsw_sisl(j_bc)])
          end do
       case(tfs_type)
          write(fates_log(),*) 'TFS conductance not used in soil'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       



       ! Update static quantities related to the rhizosphere
       call UpdateSizeDepRhizVolLenCon(sites(s), bc_in(s))

       ! We update the "initial" values of the rhizosphere after
       ! the previous call to make sure that the conductances are updated
       ! Now we set the prevous to the current so that the water states
       ! are not perturbed
       call SavePreviousRhizVolumes(sites(s))

    end do

    call UpdateH2OVeg(nsites,sites,bc_out)

    return
  end subroutine RestartHydrStates

  ! ====================================================================================

  subroutine InitPlantHydStates(site, cohort)

    ! REQUIRED INPUTS:
    !
    !  ccohort_hydr%z_node_troot(:)
    !  ccohort_hydr%z_node_aroot
    !  ccohort_hydr%z_node_ag
    ! 
    ! !DESCRIPTION: 
    !
    ! !USES:

    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target   :: site   ! current site pointer
    type(ed_cohort_type), intent(inout), target :: cohort ! current cohort pointer
    !
    ! !LOCAL VARIABLES:
    type(ed_site_hydr_type), pointer   :: site_hydr
    type(ed_cohort_hydr_type), pointer :: cohort_hydr
    integer  :: j,k                                 ! layer and node indices
    integer  :: ft                                  ! functional type index
    real(r8) :: psi_rhiz1                           ! pressure in first rhizosphere shell [MPa]
    real(r8) :: dz                                  ! depth of the current layer [m]
    real(r8) :: h_aroot_mean                        ! minimum total potential of absorbing roots
    real(r8), parameter :: psi_aroot_init = -0.2_r8 ! Initialize aroots with -0.2 MPa
    real(r8), parameter :: dh_dz = 0.02_r8          ! amount to decrease downstream
                                                    ! compartment total potentials [MPa/meter]
    
    ! In init mode = 1, set absorbing roots to -0.2 MPa
    !              = 2, use soil as starting point, match total potentials
    !                   and then reduce plant compartment total potential by 1KPa
    !                   for transporting root node, match the lowest total potential
    !                   in absorbing roots
    integer, parameter :: init_mode = 2
    
    
    site_hydr   => site%si_hydr
    cohort_hydr => cohort%co_hydr
    ft          =  cohort%pft

    ! Set abosrbing root

    if(init_mode == 2) then
       
!       h_aroot_mean = 0._r8

       do j=1, site_hydr%nlevrhiz
          
          ! Match the potential of the absorbing root to the inner rhizosphere shell
          cohort_hydr%psi_aroot(j) = site_hydr%wrf_soil(j)%p%psi_from_th(site_hydr%h2osoi_liqvol_shell(j,1))

          ! Calculate the mean total potential (include height) of absorbing roots
!          h_aroot_mean = h_aroot_mean + cohort_hydr%psi_aroot(j) + mpa_per_pa*denh2o*grav_earth*(-site_hydr%zi_rhiz(j))
          
          cohort_hydr%th_aroot(j) = wrf_plant(aroot_p_media,ft)%p%th_from_psi(cohort_hydr%psi_aroot(j))
          cohort_hydr%ftc_aroot(j) = wkf_plant(aroot_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_aroot(j))
       end do
       
    else
       
       do j=1, site_hydr%nlevrhiz
          cohort_hydr%psi_aroot(j) = psi_aroot_init
          ! Calculate the mean total potential (include height) of absorbing roots
!          h_aroot_mean = h_aroot_mean + cohort_hydr%psi_aroot(j) + mpa_per_pa*denh2o*grav_earth*(-site_hydr%zi_rhiz(j))
          cohort_hydr%th_aroot(j) = wrf_plant(aroot_p_media,ft)%p%th_from_psi(cohort_hydr%psi_aroot(j))
          cohort_hydr%ftc_aroot(j) = wkf_plant(aroot_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_aroot(j))
       end do
    end if
    
    !h_aroot_mean = h_aroot_mean/real(site_hydr%nlevrhiz,r8)
    h_aroot_mean = minval(cohort_hydr%psi_aroot(:) + mpa_per_pa*denh2o*grav_earth*(-site_hydr%zi_rhiz(:)))

    ! initialize plant water potentials with slight potential gradient (or zero) (dh/dz = C)
    ! the assumption is made here that initial conditions for soil water will 
    ! be in (or at least close to) hydrostatic equilibrium as well, so that
    ! it doesn't matter which absorbing root layer the transporting root water


    ! Set the transporting root to be in equilibrium with mean potential
    ! of the absorbing roots, minus any gradient we add

    cohort_hydr%psi_troot = h_aroot_mean - &
         mpa_per_pa*denh2o*grav_earth*cohort_hydr%z_node_troot - dh_dz

    cohort_hydr%th_troot = wrf_plant(troot_p_media,ft)%p%th_from_psi(cohort_hydr%psi_troot)
    cohort_hydr%ftc_troot = wkf_plant(troot_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_troot)


    ! working our way up a tree, assigning water potentials that are in
    ! hydrostatic equilibrium (minus dh_dz offset) with the water potential immediately below
    dz = cohort_hydr%z_node_ag(n_hypool_ag) - cohort_hydr%z_node_troot

    cohort_hydr%psi_ag(n_hypool_ag) = cohort_hydr%psi_troot - &
         mpa_per_pa*denh2o*grav_earth*dz - dh_dz


    cohort_hydr%th_ag(n_hypool_ag) = wrf_plant(stem_p_media,ft)%p%th_from_psi(cohort_hydr%psi_ag(n_hypool_ag))
    cohort_hydr%ftc_ag(n_hypool_ag) = wkf_plant(stem_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_ag(n_hypool_ag))

    
    do k=n_hypool_ag-1, 1, -1
       dz = cohort_hydr%z_node_ag(k) - cohort_hydr%z_node_ag(k+1)
       cohort_hydr%psi_ag(k) = cohort_hydr%psi_ag(k+1) - &
            mpa_per_pa*denh2o*grav_earth*dz - &
            dh_dz

       cohort_hydr%th_ag(k) = wrf_plant(site_hydr%pm_node(k),ft)%p%th_from_psi(cohort_hydr%psi_ag(k))
       cohort_hydr%ftc_ag(k) = wkf_plant(site_hydr%pm_node(k),ft)%p%ftc_from_psi(cohort_hydr%psi_ag(k))
    end do

    cohort_hydr%errh2o_growturn_ag(:)    = 0.0_r8
    cohort_hydr%errh2o_growturn_troot    = 0.0_r8
    cohort_hydr%errh2o_growturn_aroot    = 0.0_r8
    cohort_hydr%errh2o_pheno_ag(:)       = 0.0_r8
    cohort_hydr%errh2o_pheno_troot       = 0.0_r8
    cohort_hydr%errh2o_pheno_aroot       = 0.0_r8

    !initialize cohort-level btran

    cohort_hydr%btran = wkf_plant(stomata_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_ag(1))


    !flc_gs_from_psi(cohort_hydr%psi_ag(1),cohort%pft)
    
    ! We do allow for positive pressures.
    ! But starting off with positive pressures is something we try to avoid
    if ( (cohort_hydr%psi_troot>0.0_r8) .or. &
          any(cohort_hydr%psi_ag(:)>0._r8) .or. &
         any(cohort_hydr%psi_aroot(:)>0._r8) ) then
       write(fates_log(),*) 'Initialized plant compartments with positive pressure?'
       write(fates_log(),*) 'psi troot: ',cohort_hydr%psi_troot
       write(fates_log(),*) 'psi ag(:): ',cohort_hydr%psi_ag(:)
       write(fates_log(),*) 'psi_aroot(:): ',cohort_hydr%psi_aroot(:)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    

  end subroutine InitPlantHydStates
  
  ! =====================================================================================

  subroutine UpdatePlantPsiFTCFromTheta(ccohort,csite_hydr)
    
    ! This subroutine updates the potential and the fractional
    ! of total conductivity based on the relative water
    ! content
    ! Arguments
    type(ed_cohort_type),intent(inout), target :: ccohort
    type(ed_site_hydr_type),intent(in), target :: csite_hydr

    ! Locals
    integer :: ft  ! Plant functional type
    integer :: k   ! loop index for compartments
    integer :: j   ! Loop index for soil layers
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr


    
    ccohort_hydr => ccohort%co_hydr
    ft = ccohort%pft
    
    ! Update Psi and FTC in above-ground compartments
    ! -----------------------------------------------------------------------------------
    do k = 1,n_hypool_leaf
        ccohort_hydr%psi_ag(k) = wrf_plant(leaf_p_media,ft)%p%psi_from_th(ccohort_hydr%th_ag(k)) 
        ccohort_hydr%ftc_ag(k) = wkf_plant(leaf_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(k))
    end do

    do k = n_hypool_leaf+1, n_hypool_ag 
       ccohort_hydr%psi_ag(k) = wrf_plant(stem_p_media,ft)%p%psi_from_th(ccohort_hydr%th_ag(k))
       ccohort_hydr%ftc_ag(k) = wkf_plant(stem_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(k))
    end do

    ! Update the Psi and FTC for the transporting root compartment
    ccohort_hydr%psi_troot = wrf_plant(troot_p_media,ft)%p%psi_from_th(ccohort_hydr%th_troot)
    ccohort_hydr%ftc_troot = wkf_plant(troot_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_troot)

    ! Update the Psi and FTC for the absorbing roots
    do j = 1, csite_hydr%nlevrhiz
       ccohort_hydr%psi_aroot(j) = wrf_plant(aroot_p_media,ft)%p%psi_from_th(ccohort_hydr%th_aroot(j)) 
       ccohort_hydr%ftc_aroot(j) = wkf_plant(aroot_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_aroot(j))
    end do

    return
  end subroutine UpdatePlantPsiFTCFromTheta


  ! =====================================================================================


  subroutine UpdatePlantHydrNodes(ccohort_hydr,ft,plant_height,csite_hydr)

    ! --------------------------------------------------------------------------------
    ! This subroutine calculates the nodal heights critical to hydraulics in the plant
    !
    ! Inputs:  Plant height
    !          Plant functional type
    !          Number of soil hydraulic layers
    !
    ! Outputs: cohort_hydr%z_node_ag(:)      
    !                     %z_lower_ag(:)
    !                     %z_upper_ag(:)
    !                     %z_node_troot
    !                     %z_node_aroot(:)
    ! --------------------------------------------------------------------------------

    ! Arguments
    type(ed_cohort_hydr_type), intent(inout) :: ccohort_hydr 
    integer,intent(in)                       :: ft            ! plant functional type index
    real(r8), intent(in)                     :: plant_height  ! [m]
    type(ed_site_hydr_type), intent(in)      :: csite_hydr


    ! Locals
    integer  :: nlevrhiz      ! number of rhizosphere layers
    real(r8) :: roota         ! root profile parameter a zeng2001_crootfr
    real(r8) :: rootb         ! root profile parameter b zeng2001_crootfr
    real(r8) :: crown_depth   ! crown depth for the plant [m]
    real(r8) :: dz_canopy     ! discrete crown depth intervals [m]
    real(r8) :: z_stem        ! the height of the plants stem below crown [m]
    real(r8) :: dz_stem       ! vertical stem discretization                           [m]
    real(r8) :: dcumul_rf     ! cumulative root distribution discretization            [-]
    real(r8) :: cumul_rf      ! cumulative root distribution where depth is determined [-]
    real(r8) :: z_cumul_rf    ! depth at which cumul_rf occurs                         [m]
    integer  :: k             ! Loop counter for compartments

    ! Crown Nodes
    ! in special case where n_hypool_leaf = 1, the node height of the canopy
    ! water pool is 1/2 the distance from the bottom of the canopy to the top of the tree
    roota                      = EDPftvarcon_inst%fnrt_prof_a(ft)
    rootb                      = EDPftvarcon_inst%fnrt_prof_b(ft)
    nlevrhiz                   = csite_hydr%nlevrhiz
    call CrownDepth(plant_height,crown_depth)

    dz_canopy                  = crown_depth / real(n_hypool_leaf,r8)
    do k=1,n_hypool_leaf
       ccohort_hydr%z_lower_ag(k)   = plant_height - dz_canopy*real(k,r8)
       ccohort_hydr%z_node_ag(k)    = ccohort_hydr%z_lower_ag(k) + 0.5_r8*dz_canopy
       ccohort_hydr%z_upper_ag(k)   = ccohort_hydr%z_lower_ag(k) + dz_canopy
    enddo


    ! Stem Nodes
    ! in special case where n_hypool_stem = 1, the node height of the stem water pool is
    ! 1/2 the height from the ground to the bottom of the canopy
    z_stem                     = plant_height - crown_depth
    dz_stem                    = z_stem / real(n_hypool_stem,r8)
    do k=n_hypool_leaf+1,n_hypool_ag
       ccohort_hydr%z_upper_ag(k)   = real(n_hypool_stem - (k - 1 - n_hypool_leaf),r8)*dz_stem
       ccohort_hydr%z_node_ag(k)    = ccohort_hydr%z_upper_ag(k) - 0.5_r8*dz_stem
       ccohort_hydr%z_lower_ag(k)   = ccohort_hydr%z_upper_ag(k) - dz_stem
    enddo

    ! Transporting Root Node depth [m] (negative from surface)

    call bisect_rootfr(roota, rootb, 0._r8, 1.E10_r8, &
         0.001_r8, 0.001_r8, 0.5_r8, z_cumul_rf)
    z_cumul_rf =  min(z_cumul_rf, abs(csite_hydr%zi_rhiz(nlevrhiz)))
    ccohort_hydr%z_node_troot = -z_cumul_rf
    
    return
  end subroutine UpdatePlantHydrNodes

  ! =====================================================================================

  subroutine SavePreviousCompartmentVolumes(ccohort_hydr)

    type(ed_cohort_hydr_type),intent(inout) :: ccohort_hydr


    ! Saving the current compartment volumes into an "initial" save-space
    ! allows us to see how the compartments change size when plants
    ! change size and effect water contents

    ccohort_hydr%v_ag_init(:)          =  ccohort_hydr%v_ag(:)
    ccohort_hydr%v_troot_init          =  ccohort_hydr%v_troot
    ccohort_hydr%v_aroot_layer_init(:) =  ccohort_hydr%v_aroot_layer(:)

    return
  end subroutine SavePreviousCompartmentVolumes

  ! =====================================================================================

  subroutine UpdateSizeDepPlantHydProps(currentSite,ccohort,bc_in)


    ! DESCRIPTION: Updates absorbing root length (total and its vertical distribution)
    ! as well as the consequential change in the size of the 'representative' rhizosphere
    ! shell radii, volumes, and compartment volumes of plant tissues

    ! !USES:
    use shr_sys_mod        , only : shr_sys_abort

    ! ARGUMENTS:
    type(ed_site_type)     , intent(in)             :: currentSite ! Site stuff
    type(ed_cohort_type)   , intent(inout)          :: ccohort     ! current cohort pointer
    type(bc_in_type)       , intent(in)             :: bc_in       ! Boundary Conditions

    ! Locals
    integer                            :: nlevrhiz             ! Number of total soil layers
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    integer                            :: ft

    nlevrhiz                =  currentSite%si_hydr%nlevrhiz
    ccohort_hydr               => ccohort%co_hydr
    ft                         =  ccohort%pft

    ! Save the current vegetation compartment volumes into
    ! a save space so that it can be compared with the updated quantity.

    call SavePreviousCompartmentVolumes(ccohort_hydr)

    ! This updates all of the z_node positions
    call UpdatePlantHydrNodes(ccohort_hydr,ft,ccohort%hite,currentSite%si_hydr)

    ! This updates plant compartment volumes, lengths and 
    ! maximum conductances. Make sure for already
    ! initialized vegetation, that SavePreviousCompartment
    ! volumes, and UpdatePlantHydrNodes is called prior to this.
    call UpdatePlantHydrLenVol(ccohort,currentSite%si_hydr)

    ! This updates the Kmax's of the plant's compartments
    call UpdatePlantKmax(ccohort_hydr,ccohort,currentsite%si_hydr)


  end subroutine UpdateSizeDepPlantHydProps

  ! =====================================================================================

  subroutine UpdatePlantHydrLenVol(ccohort,site_hydr)

    ! -----------------------------------------------------------------------------------
    ! This subroutine calculates two attributes of a plant:
    ! 1) the volumes of storage compartments in the plants
    ! 2) the lenghts of the organs 
    ! These are not dependent on the hydraulic state of the
    ! plant, it is more about the structural characteristics and how much biomass
    ! is present in the different tissues.
    !
    ! Inputs, plant geometries, plant carbon pools, z_node values
    !
    ! -----------------------------------------------------------------------------------

    ! Arguments
    type(ed_cohort_type),intent(inout)  :: ccohort
    type(ed_site_hydr_type),intent(in)  :: site_hydr 
    
    type(ed_cohort_hydr_type),pointer :: ccohort_hydr     ! Plant hydraulics structure
    integer  :: j,k
    integer  :: ft                           ! Plant functional type index 
    real(r8) :: roota                        ! root profile parameter a zeng2001_crootfr
    real(r8) :: rootb                        ! root profile parameter b zeng2001_crootfr
    real(r8) :: leaf_c                       ! Current amount of leaf carbon in the plant                            [kg]
    real(r8) :: leaf_c_target                ! Target leaf carbon (with some conditions) [kgC]
    real(r8) :: fnrt_c                       ! Current amount of fine-root carbon in the plant                       [kg]
    real(r8) :: sapw_c                       ! Current amount of sapwood carbon in the plant                         [kg]
    real(r8) :: struct_c                     ! Current amount of structural carbon in the plant                      [kg]
    real(r8) :: woody_bg_c                   ! belowground woody biomass in carbon units                             [kgC/indiv]
    real(r8) :: z_stem                       ! the height of the plants stem below crown [m]
    real(r8) :: sla                          ! specific leaf area                                                    [cm2/g]
    real(r8) :: v_aroot_tot                  ! total compartment volume of all absorbing roots for cohort            [m3]
    real(r8) :: l_aroot_tot                  ! total length of absorbing roots for cohrot                            [m]
    real(r8) :: denleaf                      ! leaf dry mass per unit fresh leaf volume                              [kg/m3]     
    real(r8) :: a_sapwood                    ! sapwood area                                                          [m2]
    real(r8) :: a_sapwood_target             ! sapwood cross-section area at reference height, at target biomass     [m2]
    real(r8) :: sapw_c_target                ! sapwood carbon, at target                                             [kgC]
    real(r8) :: v_sapwood                    ! sapwood volume                                                        [m3]
    real(r8) :: v_troot                      ! transporting root volume                                              [m3/indiv]
    real(r8) :: rootfr                       ! mass fraction of roots in each layer                                  [kg/kg]
    real(r8) :: crown_depth                  ! Depth of the plant's crown [m]
    real(r8) :: norm                         ! total root fraction used <1
    integer  :: nlevrhiz                     ! number of rhizosphere levels


    ! We allow the transporting root to donate a fraction of its volume to the absorbing
    ! roots to help mitigate numerical issues due to very small volumes. This is the
    ! fraction the transporting roots donate to those layers
    real(r8), parameter :: t2aroot_vol_donate_frac = 0.65_r8

    real(r8), parameter :: min_leaf_frac = 0.1_r8   ! Fraction of maximum leaf carbon that
                                                    ! we set as our lower cap on leaf volume
    real(r8), parameter :: min_trim      = 0.1_r8   ! The lower cap on trimming function used
                                                    ! to estimate maximum leaf carbon
    
    ccohort_hydr => ccohort%co_hydr
    ft           = ccohort%pft
    nlevrhiz = site_hydr%nlevrhiz
    leaf_c   = ccohort%prt%GetState(leaf_organ, all_carbon_elements)
    sapw_c   = ccohort%prt%GetState(sapw_organ, all_carbon_elements)
    fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)
    struct_c = ccohort%prt%GetState(struct_organ, all_carbon_elements)

    roota    = EDPftvarcon_inst%fnrt_prof_a(ft)
    rootb    = EDPftvarcon_inst%fnrt_prof_b(ft)

    ! Leaf Volumes
    ! -----------------------------------------------------------------------------------

    ! NOTE: SLATOP currently does not use any vertical scaling functions
    ! but that may not be so forever. ie sla = slatop (RGK-082017)
    ! m2/gC * cm2/m2 -> cm2/gC
    
    sla                        = EDPftvarcon_inst%slatop(ft) * cm2_per_m2 
    
    ! empirical regression data from leaves at Caxiuana (~ 8 spp)
    denleaf                    = -2.3231_r8*sla/EDPftvarcon_inst%c2b(ft) + 781.899_r8    
 
    ! Leaf volumes
    ! Note: Leaf volumes of zero is problematic for two reasons.  Zero volumes create
    ! numerical difficulties, and they could also create problems when a leaf is trying
    ! to re-flush.
    ! Therefore, if the leaf is in an "off" status, then we do not update the leaf
    ! volume.  This way the volume is where it was when it dropped, and this is consistent
    ! with the theory that leaf water potentials drive growth and re-flushing, not the
    ! other way around.  However, it is possible that we may have recruits with an
    ! "off" status (due to external seed rain active during a dry or cold season). If a
    ! cohort is newly created, we must give it a starting volume.
    ! We also place a lower bound on how low the leaf volume is allowed to go, which is 10%
    ! of the plant's carrying capacity.

    
    ! [kgC] * [kg/kgC] / [kg/m3] -> [m3]


    ! Get the target, or rather, maximum leaf carrying capacity of plant
    ! Lets also avoid super-low targets that have very low trimming functions

    call bleaf(ccohort%dbh,ccohort%pft,max(ccohort%canopy_trim,min_trim),leaf_c_target)

    if( (ccohort%status_coh == leaves_on) .or. ccohort_hydr%is_newly_recruited ) then
       ccohort_hydr%v_ag(1:n_hypool_leaf) = max(leaf_c,min_leaf_frac*leaf_c_target) * &
            EDPftvarcon_inst%c2b(ft) / denleaf/ real(n_hypool_leaf,r8)
    end if
    
    ! Step sapwood volume
    ! -----------------------------------------------------------------------------------

    ! BOC...may be needed for testing/comparison w/ v_sapwood 
    ! kg  / ( g cm-3 * cm3/m3 * kg/g ) -> m3    
    ! v_stem       = b_stem_biom / (EDPftvarcon_inst%wood_density(ft) * kg_per_g * cm3_per_m3 ) 

    ! calculate the sapwood cross-sectional area
    call bsap_allom(ccohort%dbh,ccohort%pft,ccohort%canopy_trim,a_sapwood_target,sapw_c_target)

    ! uncomment this if you want to use
    ! the actual sapwood, which may be lower than target due to branchfall.
    a_sapwood = a_sapwood_target  ! * sapw_c / sapw_c_target

    ! alternative cross section calculation
    ! a_sapwood    = a_leaf_tot / ( 0.001_r8 + 0.025_r8 * ccohort%hite ) * 1.e-4_r8

    call CrownDepth(ccohort%hite,crown_depth)
    z_stem       = ccohort%hite - crown_depth
    v_sapwood    = a_sapwood * z_stem    ! + 0.333_r8*a_sapwood*crown_depth
    ccohort_hydr%v_ag(n_hypool_leaf+1:n_hypool_ag) = v_sapwood / n_hypool_stem


    ! Determine belowground biomass as a function of total (sapwood, heartwood, 
    ! leaf, fine root) biomass then subtract out the fine root biomass to get 
    ! coarse (transporting) root biomass

    woody_bg_c = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(ft)) * (sapw_c + struct_c)
    
    v_troot                    = woody_bg_c * EDPftvarcon_inst%c2b(ft) / & 
          (EDPftvarcon_inst%wood_density(ft)*kg_per_g*cm3_per_m3)
    
    
    ! Estimate absorbing root total length (all layers)
    ! SRL is in m/g
    ! [m] = [kgC]*1000[g/kg]*[kg/kgC]*[m/g]
    ! ------------------------------------------------------------------------------
    l_aroot_tot        = fnrt_c*g_per_kg*EDPftvarcon_inst%c2b(ft)*EDPftvarcon_inst%hydr_srl(ft)
    
    
    ! Estimate absorbing root volume (all layers)
    ! ------------------------------------------------------------------------------
    v_aroot_tot        = pi_const * (EDPftvarcon_inst%hydr_rs2(ft)**2._r8) * l_aroot_tot

    ! The transporting root donates some of its volume
    ! to the layer-by-layer absorbing root (which is now a hybrid compartment)
    ! ------------------------------------------------------------------------------
    ccohort_hydr%v_troot = (1._r8-t2aroot_vol_donate_frac) * v_troot
    
    ! Partition the total absorbing root lengths and volumes into the active soil layers
    ! We have a condition, where we may ignore the first layer
    ! ------------------------------------------------------------------------------
    
    norm = 1._r8 - &
          zeng2001_crootfr(roota, rootb,site_hydr%zi_rhiz(1)-site_hydr%dz_rhiz(1), site_hydr%zi_rhiz(nlevrhiz))
    
    do j=1,nlevrhiz
        
        rootfr = norm*(zeng2001_crootfr(roota, rootb, site_hydr%zi_rhiz(j),site_hydr%zi_rhiz(nlevrhiz)) - &
              zeng2001_crootfr(roota, rootb, site_hydr%zi_rhiz(j)-site_hydr%dz_rhiz(j),site_hydr%zi_rhiz(nlevrhiz)))
        
        ccohort_hydr%l_aroot_layer(j) = rootfr*l_aroot_tot
        
        ! This is a hybrid absorbing root and transporting root volume
        ccohort_hydr%v_aroot_layer(j) = rootfr*(v_aroot_tot + t2aroot_vol_donate_frac*v_troot)

    end do
       
    return
  end subroutine UpdatePlantHydrLenVol

  ! =====================================================================================

  subroutine UpdateSizeDepPlantHydStates(currentSite,ccohort)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    use FatesUtilsMod  , only : check_var_real

    ! !ARGUMENTS:
    type(ed_site_type)    , intent(in)    :: currentSite ! Site stuff
    type(ed_cohort_type)   , intent(inout) :: ccohort
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type),pointer :: csite_hydr
    integer  :: j,k,FT                       ! indices
    integer  :: err_code = 0
    real(r8) :: th_ag_uncorr(      n_hypool_ag) ! uncorrected aboveground water content[m3 m-3]
    real(r8) :: th_troot_uncorr                 ! uncorrected transporting root water content[m3 m-3]
    real(r8) :: th_aroot_uncorr(currentSite%si_hydr%nlevrhiz)    ! uncorrected absorbing root water content[m3 m-3] 
    real(r8), parameter :: small_theta_num = 1.e-7_r8  ! avoids theta values equalling thr or ths         [m3 m-3]

    integer :: nstep !number of time steps
    !-----------------------------------------------------------------------

    ccohort_hydr => ccohort%co_hydr
    FT      =  cCohort%pft
    
    associate(pm_node => currentSite%si_hydr%pm_node)
    
    ! MAYBE ADD A NAN CATCH?  If UpdateSizeDepPlantHydProps() was not called twice prior to the first
    ! time this routine is called for a new cohort, then v_ag_init(k) will be a nan.
    ! It should be ok, but may be vulnerable if code is changed (RGK 02-2017)

    ! UPDATE WATER CONTENTS (assume water for growth comes from within tissue itself
    ! -- apply water mass conservation)

    do k=1,n_hypool_leaf
       if( ccohort_hydr%v_ag(k) > nearzero ) then
          th_ag_uncorr(k)    = ccohort_hydr%th_ag(k)   * &
               ccohort_hydr%v_ag_init(k) /ccohort_hydr%v_ag(k)
          ccohort_hydr%th_ag(k) = constrain_water_contents(th_ag_uncorr(k), small_theta_num, ft, pm_node(k))
       else
          th_ag_uncorr(k)    = ccohort_hydr%th_ag(k) 
       end if
    end do

    do k=n_hypool_leaf+1,n_hypool_ag
       th_ag_uncorr(k)    = ccohort_hydr%th_ag(k)   * &
            ccohort_hydr%v_ag_init(k) /ccohort_hydr%v_ag(k)
       ccohort_hydr%th_ag(k) = constrain_water_contents(th_ag_uncorr(k), small_theta_num, ft, pm_node(k))
    enddo

    th_troot_uncorr = ccohort_hydr%th_troot * ccohort_hydr%v_troot_init /ccohort_hydr%v_troot
    ccohort_hydr%th_troot = constrain_water_contents(th_troot_uncorr, small_theta_num, ft, pm_node(3))


    ccohort_hydr%errh2o_growturn_aroot = 0._r8
    do j=1,currentSite%si_hydr%nlevrhiz
       th_aroot_uncorr(j) = ccohort_hydr%th_aroot(j) * &
            ccohort_hydr%v_aroot_layer_init(j)/ccohort_hydr%v_aroot_layer(j)
       ccohort_hydr%th_aroot(j) = constrain_water_contents(th_aroot_uncorr(j), small_theta_num, ft, pm_node(4))
       ccohort_hydr%errh2o_growturn_aroot = ccohort_hydr%errh2o_growturn_aroot + & 
            denh2o*cCohort%n*AREA_INV*(ccohort_hydr%th_aroot(j)-th_aroot_uncorr(j))*ccohort_hydr%v_aroot_layer(j)
    enddo

    ! Storing mass balance error
    ! + means water created; - means water destroyed
    ccohort_hydr%errh2o_growturn_ag(:) = denh2o*cCohort%n*AREA_INV*ccohort_hydr%v_ag(:) * & 
         (ccohort_hydr%th_ag(:)-th_ag_uncorr(:))
    ccohort_hydr%errh2o_growturn_troot = denh2o*cCohort%n*AREA_INV*ccohort_hydr%v_troot * &
         (ccohort_hydr%th_troot-th_troot_uncorr)

    csite_hydr =>currentSite%si_hydr
    csite_hydr%h2oveg_growturn_err = csite_hydr%h2oveg_growturn_err + &
         sum(ccohort_hydr%errh2o_growturn_ag(:)) + & 
         ccohort_hydr%errh2o_growturn_troot      + &
         ccohort_hydr%errh2o_growturn_aroot


    ! UPDATES OF WATER POTENTIALS ARE DONE PRIOR TO RICHARDS' SOLUTION WITHIN FATESPLANTHYDRAULICSMOD.F90
    end associate

  end subroutine UpdateSizeDepPlantHydStates

  ! =====================================================================================

  function constrain_water_contents(th_uncorr, delta, ft, pm_type) result(th_corr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: th_uncorr ! uncorrected water content (m3 m-3)
    real(r8) , intent(in) :: delta
    integer , intent(in)  :: ft
    integer , intent(in)  :: pm_type
    !
    ! !Local:
    real(r8) :: thr                    ! residual water content (m3 m-3)
    real(r8) :: ths                    ! saturated water content (m3 m-3)
    !
    ! !RESULT
    real(r8) :: th_corr                ! corrected water content
    !
    !------------------------------------------------------------------------
    ths     = EDPftvarcon_inst%hydr_thetas_node(ft,pm_type)
    thr     = EDPftvarcon_inst%hydr_resid_node(ft,pm_type)
    th_corr = max((thr+delta),min((ths-delta),th_uncorr))

    return

  end function constrain_water_contents

  ! =====================================================================================

  subroutine CopyCohortHydraulics(newCohort, oldCohort)

    ! Arguments
    type(ed_cohort_type), intent(inout), target :: newCohort
    type(ed_cohort_type), intent(inout), target :: oldCohort

    ! Locals
    type(ed_cohort_hydr_type), pointer :: ncohort_hydr
    type(ed_cohort_hydr_type), pointer :: ocohort_hydr


    ncohort_hydr => newCohort%co_hydr
    ocohort_hydr => oldCohort%co_hydr

    ! Node heights
    ncohort_hydr%z_node_ag             = ocohort_hydr%z_node_ag
    ncohort_hydr%z_upper_ag            = ocohort_hydr%z_upper_ag
    ncohort_hydr%z_lower_ag            = ocohort_hydr%z_lower_ag
    ncohort_hydr%z_node_troot          = ocohort_hydr%z_node_troot

    ! Compartment kmax's
    ncohort_hydr%kmax_petiole_to_leaf  = ocohort_hydr%kmax_petiole_to_leaf
    ncohort_hydr%kmax_stem_lower       = ocohort_hydr%kmax_stem_lower
    ncohort_hydr%kmax_stem_upper       = ocohort_hydr%kmax_stem_upper
    ncohort_hydr%kmax_troot_upper      = ocohort_hydr%kmax_troot_upper
    ncohort_hydr%kmax_troot_lower      = ocohort_hydr%kmax_troot_lower
    ncohort_hydr%kmax_aroot_upper      = ocohort_hydr%kmax_aroot_upper
    ncohort_hydr%kmax_aroot_lower      = ocohort_hydr%kmax_aroot_lower
    ncohort_hydr%kmax_aroot_radial_in  = ocohort_hydr%kmax_aroot_radial_in
    ncohort_hydr%kmax_aroot_radial_out = ocohort_hydr%kmax_aroot_radial_out

    ! Compartment volumes
    ncohort_hydr%v_ag_init             = ocohort_hydr%v_ag_init
    ncohort_hydr%v_ag                  = ocohort_hydr%v_ag
    ncohort_hydr%v_troot_init          = ocohort_hydr%v_troot_init
    ncohort_hydr%v_troot               = ocohort_hydr%v_troot
    ncohort_hydr%v_aroot_layer_init    = ocohort_hydr%v_aroot_layer_init
    ncohort_hydr%v_aroot_layer         = ocohort_hydr%v_aroot_layer
    ncohort_hydr%l_aroot_layer         = ocohort_hydr%l_aroot_layer

    ! State Variables
    ncohort_hydr%th_ag                 = ocohort_hydr%th_ag
    ncohort_hydr%th_troot              = ocohort_hydr%th_troot
    ncohort_hydr%th_aroot              = ocohort_hydr%th_aroot
    ncohort_hydr%psi_ag                = ocohort_hydr%psi_ag
    ncohort_hydr%psi_troot             = ocohort_hydr%psi_troot
    ncohort_hydr%psi_aroot             = ocohort_hydr%psi_aroot
    ncohort_hydr%ftc_ag                = ocohort_hydr%ftc_ag
    ncohort_hydr%ftc_troot             = ocohort_hydr%ftc_troot
    ncohort_hydr%ftc_aroot             = ocohort_hydr%ftc_aroot

    ! Other
    ncohort_hydr%btran                 = ocohort_hydr%btran
    ncohort_hydr%supsub_flag           = ocohort_hydr%supsub_flag
    ncohort_hydr%iterh1                = ocohort_hydr%iterh1
    ncohort_hydr%iterh2                = ocohort_hydr%iterh2
    ncohort_hydr%iterlayer             = ocohort_hydr%iterlayer
    ncohort_hydr%errh2o                = ocohort_hydr%errh2o
    ncohort_hydr%errh2o_growturn_ag    = ocohort_hydr%errh2o_growturn_ag
    ncohort_hydr%errh2o_pheno_ag       = ocohort_hydr%errh2o_pheno_ag
    ncohort_hydr%errh2o_growturn_troot = ocohort_hydr%errh2o_growturn_troot
    ncohort_hydr%errh2o_pheno_troot    = ocohort_hydr%errh2o_pheno_troot
    ncohort_hydr%errh2o_growturn_aroot = ocohort_hydr%errh2o_growturn_aroot
    ncohort_hydr%errh2o_pheno_aroot    = ocohort_hydr%errh2o_pheno_aroot

    ! BC PLANT HYDRAULICS - flux terms
    ncohort_hydr%qtop                  = ocohort_hydr%qtop

    ncohort_hydr%is_newly_recruited    = ocohort_hydr%is_newly_recruited

  end subroutine CopyCohortHydraulics

  ! =====================================================================================
  subroutine FuseCohortHydraulics(currentSite,currentCohort, nextCohort, bc_in, newn)


    type(ed_cohort_type), intent(inout), target :: currentCohort ! current cohort
    type(ed_cohort_type), intent(inout), target :: nextCohort    ! next (donor) cohort
    type(ed_site_type), intent(inout), target :: currentSite    ! current site

    type(bc_in_type), intent(in)                :: bc_in
    real(r8), intent(in)                        :: newn

    ! !LOCAL VARIABLES:
    type(ed_site_hydr_type), pointer :: site_hydr
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr  ! current cohort hydraulics derived type
    type(ed_cohort_hydr_type), pointer :: ncohort_hydr  ! donor (next) cohort hydraulics d type
    integer  :: j,k                                     ! indices
    integer  :: ft

    site_hydr => currentSite%si_hydr

    ccohort_hydr => currentCohort%co_hydr
    ncohort_hydr => nextCohort%co_hydr

    ccohort_hydr%th_ag(:)    = (currentCohort%n*ccohort_hydr%th_ag(:)    + &
         nextCohort%n*ncohort_hydr%th_ag(:))/newn
    ccohort_hydr%th_troot    = (currentCohort%n*ccohort_hydr%th_troot    + &
         nextCohort%n*ncohort_hydr%th_troot)/newn
    ccohort_hydr%th_aroot(:) = (currentCohort%n*ccohort_hydr%th_aroot(:) + &
         nextCohort%n*ncohort_hydr%th_aroot(:))/newn
    ccohort_hydr%supsub_flag = 0

    ! Only save the iteration counters for the worse of the two cohorts
    if(ncohort_hydr%iterh1 > ccohort_hydr%iterh1)then
       ccohort_hydr%iterh1      = ncohort_hydr%iterh1
       ccohort_hydr%iterh2      = ncohort_hydr%iterh2
       ccohort_hydr%iterlayer   = ncohort_hydr%iterlayer
    end if

    ft = currentCohort%pft
    do k=1,n_hypool_leaf
       ccohort_hydr%psi_ag(k) = wrf_plant(leaf_p_media,ft)%p%psi_from_th(ccohort_hydr%th_ag(k)) 
       ccohort_hydr%ftc_ag(k) = wkf_plant(leaf_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(k)) 
    end do

    do k = n_hypool_leaf+1,n_hypool_ag
       ccohort_hydr%psi_ag(k) = wrf_plant(stem_p_media,ft)%p%psi_from_th(ccohort_hydr%th_ag(k)) 
       ccohort_hydr%ftc_ag(k) = wkf_plant(stem_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(k))
    end do

    ccohort_hydr%psi_troot = wrf_plant(troot_p_media,ft)%p%psi_from_th(ccohort_hydr%th_troot) 
    ccohort_hydr%ftc_troot = wkf_plant(troot_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_troot)

    do j=1,site_hydr%nlevrhiz
       ccohort_hydr%psi_aroot(j) = wrf_plant(aroot_p_media,ft)%p%psi_from_th(ccohort_hydr%th_aroot(j)) 
       ccohort_hydr%ftc_aroot(j) = wkf_plant(aroot_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_aroot(j))
    end do


    ccohort_hydr%btran = wkf_plant(stomata_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(1))

    ccohort_hydr%qtop     = (currentCohort%n*ccohort_hydr%qtop     + &
         nextCohort%n*ncohort_hydr%qtop)/newn

    ccohort_hydr%errh2o                   = (currentCohort%n*ccohort_hydr%errh2o                   + &
         nextCohort%n*ncohort_hydr%errh2o)/newn
    ccohort_hydr%errh2o_growturn_ag(:)    = (currentCohort%n*ccohort_hydr%errh2o_growturn_ag(:)    + &
         nextCohort%n*ncohort_hydr%errh2o_growturn_ag(:))/newn
    ccohort_hydr%errh2o_pheno_ag(:)       = (currentCohort%n*ccohort_hydr%errh2o_pheno_ag(:)       + &
         nextCohort%n*ncohort_hydr%errh2o_pheno_ag(:))/newn
    ccohort_hydr%errh2o_growturn_troot    = (currentCohort%n*ccohort_hydr%errh2o_growturn_troot    + &
         nextCohort%n*ncohort_hydr%errh2o_growturn_troot)/newn
    ccohort_hydr%errh2o_pheno_troot       = (currentCohort%n*ccohort_hydr%errh2o_pheno_troot       + &
         nextCohort%n*ncohort_hydr%errh2o_pheno_troot)/newn
    ccohort_hydr%errh2o_growturn_aroot    = (currentCohort%n*ccohort_hydr%errh2o_growturn_aroot + &
         nextCohort%n*ncohort_hydr%errh2o_growturn_aroot)/newn
    ccohort_hydr%errh2o_pheno_aroot       = (currentCohort%n*ccohort_hydr%errh2o_pheno_aroot    + &
         nextCohort%n*ncohort_hydr%errh2o_pheno_aroot)/newn

    ccohort_hydr%is_newly_recruited        = .false.

  end subroutine FuseCohortHydraulics

  ! =====================================================================================
  ! Initialization Routines
  ! =====================================================================================

  subroutine InitHydrCohort(currentSite,currentCohort)

    ! Arguments
    type(ed_site_type), target   :: currentSite
    type(ed_cohort_type), target :: currentCohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr

    if ( hlm_use_planthydro.eq.ifalse ) return
    allocate(ccohort_hydr)
    currentCohort%co_hydr => ccohort_hydr
    call ccohort_hydr%AllocateHydrCohortArrays(currentSite%si_hydr%nlevrhiz)

    ccohort_hydr%is_newly_recruited = .false. 

  end subroutine InitHydrCohort

  ! =====================================================================================
  subroutine DeallocateHydrCohort(currentCohort)

    ! Arguments
    type(ed_cohort_type), target :: currentCohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr

    if ( hlm_use_planthydro.eq.ifalse ) return

    ccohort_hydr => currentCohort%co_hydr
    call ccohort_hydr%DeAllocateHydrCohortArrays()
    deallocate(ccohort_hydr)

    return
  end subroutine DeallocateHydrCohort

  ! =====================================================================================

  subroutine InitHydrSites(sites,bc_in)

    ! Arguments
    type(ed_site_type),intent(inout),target :: sites(:)
    type(bc_in_type),intent(in)             :: bc_in(:)

    ! Locals
    integer :: nsites
    integer :: s
    integer :: j
    integer :: jj
    type(ed_site_hydr_type),pointer :: csite_hydr



    if ( hlm_use_planthydro.eq.ifalse ) return

    ! Initialize any derived hydraulics parameters

    nsites = ubound(sites,1)
    do s=1,nsites
       allocate(csite_hydr)
       sites(s)%si_hydr => csite_hydr
       if ( bc_in(s)%nlevsoil > nlevsoi_hyd_max ) then
          write(fates_log(),*) 'The host land model has defined soil with'
          write(fates_log(),*) bc_in(s)%nlevsoil,' layers, for one of its columns.'
          write(fates_log(),*) 'Fates-hydro temporary array spaces with size'
          write(fates_log(),*) 'nlevsoi_hyd_max = ',nlevsoi_hyd_max,' must be larger'
          write(fates_log(),*) 'see main/FatesHydraulicsMemMod.F90'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Calculate the number of rhizosphere
       ! layers used
       if(ignore_layer1) then
          csite_hydr%i_rhiz_t = 2
          csite_hydr%i_rhiz_b = bc_in(s)%nlevsoil
       else
          csite_hydr%i_rhiz_t = 1
          csite_hydr%i_rhiz_b = bc_in(s)%nlevsoil
       end if
       
       csite_hydr%nlevrhiz = csite_hydr%i_rhiz_b-csite_hydr%i_rhiz_t+1
       call sites(s)%si_hydr%InitHydrSite(numpft,nlevsclass)

       jj=1
       do j=csite_hydr%i_rhiz_t,csite_hydr%i_rhiz_b
          csite_hydr%zi_rhiz(jj)  = bc_in(s)%zi_sisl(j)
          csite_hydr%dz_rhiz(jj) = bc_in(s)%dz_sisl(j) 
          jj=jj+1
       end do
       
   end do

  end subroutine InitHydrSites

  ! ===================================================================================
  subroutine HydrSiteColdStart(sites, bc_in )! , bc_out)


    ! Arguments
    type(ed_site_type),intent(inout),target :: sites(:)
    type(bc_in_type),intent(in)             :: bc_in(:)

    ! Local
    type(ed_site_hydr_type), pointer :: site_hydr
    real(r8) :: smp  ! matric potential temp
    real(r8) :: h2osoi_liqvol ! liquid water content (m3/m3)
    integer :: s
    integer :: j,j_bc
    integer :: nsites
    integer :: nlevrhiz
    class(wrf_type_vg), pointer :: wrf_vg
    class(wkf_type_vg), pointer :: wkf_vg
    class(wrf_type_cch), pointer :: wrf_cch
    class(wkf_type_cch), pointer :: wkf_cch


    nsites = ubound(sites,1)

    do s = 1,nsites

       site_hydr => sites(s)%si_hydr
       nlevrhiz  =  site_hydr%nlevrhiz
       
       do j = 1,nlevrhiz
          j_bc=j+site_hydr%i_rhiz_t-1
          h2osoi_liqvol = min(bc_in(s)%eff_porosity_sl(j_bc), &
               bc_in(s)%h2o_liq_sisl(j_bc)/(site_hydr%dz_rhiz(j)*denh2o))
          
          site_hydr%h2osoi_liqvol_shell(j,1:nshell) = h2osoi_liqvol
          site_hydr%h2osoi_liq_prev(j)              = bc_in(s)%h2o_liq_sisl(j_bc)
       end do
    

       site_hydr%l_aroot_layer(1:site_hydr%nlevrhiz) = 0.0_r8


       ! --------------------------------------------------------------------------------
       ! Initialize water transfer functions
       ! which include both water retention functions (WRFs)
       ! as well as the water conductance (K) functions (WKFs)
       ! But, this is only for soil!
       ! --------------------------------------------------------------------------------
       ! Initialize the Water Retention Functions
       ! -----------------------------------------------------------------------------------

       select case(soil_wrf_type)
       case(van_genuchten_type)
          do j=1,sites(s)%si_hydr%nlevrhiz
             j_bc=j+site_hydr%i_rhiz_t-1
             allocate(wrf_vg)
             site_hydr%wrf_soil(j)%p => wrf_vg
             call wrf_vg%set_wrf_param([alpha_vg, psd_vg, th_sat_vg, th_res_vg])
          end do
       case(campbell_type)
          do j=1,site_hydr%nlevrhiz
             j_bc=j+site_hydr%i_rhiz_t-1
             allocate(wrf_cch)
             site_hydr%wrf_soil(j)%p => wrf_cch
             call wrf_cch%set_wrf_param([bc_in(s)%watsat_sisl(j_bc), &  
                  (-1.0_r8)*bc_in(s)%sucsat_sisl(j_bc)*denh2o*grav_earth*mpa_per_pa*m_per_mm , & 
                  bc_in(s)%bsw_sisl(j_bc)])
          end do
       case(tfs_type)
          write(fates_log(),*) 'TFS water retention curves not available for soil'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

       ! -----------------------------------------------------------------------------------
       ! Initialize the Water Conductance (K) Functions
       ! -----------------------------------------------------------------------------------
       
       select case(soil_wkf_type)
       case(van_genuchten_type)
           do j=1,sites(s)%si_hydr%nlevrhiz
               allocate(wkf_vg)
               site_hydr%wkf_soil(j)%p => wkf_vg
               call wkf_vg%set_wkf_param([alpha_vg, psd_vg, th_sat_vg, th_res_vg, tort_vg])
           end do
        case(campbell_type)
           do j=1,sites(s)%si_hydr%nlevrhiz
              j_bc=j+site_hydr%i_rhiz_t-1
              allocate(wkf_cch)
              site_hydr%wkf_soil(j)%p => wkf_cch
              call wkf_cch%set_wkf_param([bc_in(s)%watsat_sisl(j_bc), &
                   (-1.0_r8)*bc_in(s)%sucsat_sisl(j_bc)*denh2o*grav_earth*mpa_per_pa*m_per_mm , & 
                   bc_in(s)%bsw_sisl(j_bc)])
           end do
       case(tfs_type)
           write(fates_log(),*) 'TFS conductance not used in soil'
           call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

    end do

    ! 
    !!     call UpdateH2OVeg(nsites,sites,bc_out)

    ! --------------------------------------------------------------------------------
    ! All other ed_Hydr_site_type variables are initialized elsewhere:
    !
    ! init_patch() -> UpdateSizeDepRhizHydProps -> shellgeom()
    !         this%v_shell
    !         this%r_node_shell
    !         this%r_out_shell
    !
    ! init_patch() -> UpdateSizeDepRhizHydProps()
    !         this%l_aroot_layer_init
    !         this%l_aroot_1D
    !         this%kmax_upper_shell
    !         this%kmax_lower_shell
    !
    ! hydraulics_bc()
    !         this%supsub_flag
    !         this%errh2o_hyd     =            ! hydraulics_bc
    !         this%dwat_veg       =            ! hydraulics_bc
    !
    ! ed_update_site() -> update_h2oveg()
    !         this%h2oveg
    ! --------------------------------------------------------------------------------

    return
  end subroutine HydrSiteColdStart

  ! =====================================================================================

  subroutine UpdateH2OVeg(nsites,sites,bc_out)

    ! ----------------------------------------------------------------------------------
    ! This subroutine is called following dynamics. After growth has been updated
    ! there needs to be a re-assesment of the how much liquid water is bound in the
    ! plants.  This value is necessary for water balancing in the HLM.
    ! ----------------------------------------------------------------------------------

    ! Arguments
    integer, intent(in)                       :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_out_type), intent(inout)          :: bc_out(nsites)

    ! Locals
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_type), pointer :: currentPatch
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: csite_hydr
    integer :: s
    real(r8) :: balive_patch
    integer :: nstep !number of time steps

    !for debug only
    nstep = get_nstep()

    do s = 1,nsites
       bc_out(s)%plant_stored_h2o_si = 0.0_r8
    end do

    if( hlm_use_planthydro.eq.ifalse ) return

    do s = 1,nsites

       csite_hydr => sites(s)%si_hydr
       csite_hydr%h2oveg = 0.0_r8
       currentPatch => sites(s)%oldest_patch
       do while(associated(currentPatch))         
          currentCohort=>currentPatch%tallest
          do while(associated(currentCohort))
             ccohort_hydr => currentCohort%co_hydr
             !only account for the water for not newly recruit for mass balance
             if(.not.ccohort_hydr%is_newly_recruited) then 
                csite_hydr%h2oveg = csite_hydr%h2oveg + &
                     (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:)) + &
                     ccohort_hydr%th_troot*ccohort_hydr%v_troot + &
                     sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                     denh2o*currentCohort%n
             endif

             currentCohort => currentCohort%shorter
          enddo !cohort
          currentPatch => currentPatch%younger
       enddo !end patch loop

       csite_hydr%h2oveg              = csite_hydr%h2oveg*AREA_INV

       ! Note that h2oveg_dead is incremented wherever we have litter fluxes
       ! and it will be reduced via an evaporation term
       ! growturn_err is a term to accomodate error in growth or turnover. need to be improved for future(CX) 
       bc_out(s)%plant_stored_h2o_si = csite_hydr%h2oveg + csite_hydr%h2oveg_dead - &
            csite_hydr%h2oveg_growturn_err - &
            csite_hydr%h2oveg_pheno_err-&
            csite_hydr%h2oveg_hydro_err

    end do


    return
  end subroutine UpdateH2OVeg

  !=====================================================================================
  subroutine RecruitWUptake(nsites,sites,bc_in,dtime,recruitflag)

    ! ----------------------------------------------------------------------------------
    ! This subroutine is called to calculate the water requirement for newly recruited cohorts
    ! The water update is allocated proportionally to the root biomass, which could be updated
    ! to accomodate the soil moisture and rooting depth for small seedlings (Chonggang XU).
    ! After the root water uptake, is_newly_recruited flag is set to false.
    ! Note, this routine is not accounting for the normal water uptake of new plants
    ! going forward, this routine accounts for the water that needs to be accounted for
    ! as the plants pop into existance.    
    ! ----------------------------------------------------------------------------------

    ! Arguments
    integer, intent(in)                       :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type), intent(in)              :: bc_in(nsites)    
    real(r8), intent(in)                      :: dtime !time (seconds)
    logical, intent(out)                      :: recruitflag      !flag to check if there is newly recruited cohorts

    ! Locals
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_type), pointer :: currentPatch
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: csite_hydr
    integer :: s, j, ft
    integer :: nstep !number of time steps 
    real(r8) :: rootfr !fraction of root in different soil layer
    real(r8) :: recruitw !water for newly recruited cohorts (kg water/m2/s)
    real(r8) :: recruitw_total ! total water for newly recruited cohorts (kg water/m2/s)
    real(r8) :: err !mass error of water for newly recruited cohorts (kg water/m2/s)
    real(r8) :: sumrw_uptake !sum of water take for newly recruited cohorts (kg water/m2/s)
    real(r8) :: sum_l_aroot  !sum of absorbing root lenghts
    recruitflag = .false. 
    do s = 1,nsites 
       csite_hydr => sites(s)%si_hydr
       csite_hydr%recruit_w_uptake = 0.0_r8
       currentPatch => sites(s)%oldest_patch
       recruitw_total = 0.0_r8 
       do while(associated(currentPatch)) 
          currentCohort=>currentPatch%tallest
          do while(associated(currentCohort))
             ccohort_hydr => currentCohort%co_hydr
             ft       = currentCohort%pft
             !-----------------------------------------------------------
             ! recruitment water uptake
             if(ccohort_hydr%is_newly_recruited) then
                recruitflag = .true.
                recruitw =  (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:))    + &
                             ccohort_hydr%th_troot*ccohort_hydr%v_troot                 + &
                             sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                             denh2o*currentCohort%n*AREA_INV/dtime
                recruitw_total = recruitw_total + recruitw
                sum_l_aroot = sum(ccohort_hydr%l_aroot_layer(:))
                do j=1,csite_hydr%nlevrhiz
                   rootfr = ccohort_hydr%l_aroot_layer(j)/sum_l_aroot
                   csite_hydr%recruit_w_uptake(j) = csite_hydr%recruit_w_uptake(j) + &
                        recruitw*rootfr
                end do
                ccohort_hydr%is_newly_recruited = .false.
             endif
             currentCohort=>currentCohort%shorter
          end do !cohort loop	   
          currentPatch => currentPatch%younger
       end do !patch
       !balance check
       sumrw_uptake = sum(csite_hydr%recruit_w_uptake)
       err = recruitw_total - sumrw_uptake 
       if(abs(err)>1.0e-10_r8)then
          do j=1,csite_hydr%nlevrhiz
             csite_hydr%recruit_w_uptake(j) = csite_hydr%recruit_w_uptake(j) + &
                  err*csite_hydr%recruit_w_uptake(j)/sumrw_uptake
          enddo
          write(fates_log(),*) 'math check on recruit water failed.'
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
       endif
    end do ! site loop
    
    !write(fates_log(),*) 'Calculating recruit water'
    !write(fates_log(),*) csite_hydr%recruit_w_uptake


  end subroutine RecruitWUptake

  !=====================================================================================
  subroutine ConstrainRecruitNumber(csite,ccohort, bc_in)

    ! ---------------------------------------------------------------------------
    ! This subroutine constrains the number of plants so that there is enought water 
    !  for newly recruited individuals from the soil
    ! ---------------------------------------------------------------------------

    ! Arguments
    type(ed_site_type), intent(inout), target     :: csite
    type(ed_cohort_type) , intent(inout), target  :: ccohort
    type(bc_in_type)    , intent(in)              :: bc_in 

    ! Locals
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: csite_hydr
    real(r8) :: tmp1
    real(r8) :: watres_local   !minum water content [m3/m3]
    real(r8) :: total_water !total water in rhizosphere at a specific layer (m^3 ha-1)
    real(r8) :: total_water_min !total minimum water in rhizosphere at a specific layer (m^3)
    real(r8) :: rootfr !fraction of root in different soil layer
    real(r8) :: recruitw !water for newly recruited cohorts (kg water/m2/individual)   
    real(r8) :: n, nmin !number of individuals in cohorts
    real(r8) :: sum_l_aroot
    integer :: s, j, ft

    csite_hydr => csite%si_hydr
    ccohort_hydr =>ccohort%co_hydr
    recruitw =  (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:))    + &
         ccohort_hydr%th_troot*ccohort_hydr%v_troot  + &
         sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
         denh2o
    sum_l_aroot = sum(ccohort_hydr%l_aroot_layer(:))
    do j=1,csite_hydr%nlevrhiz
       cohort_recruit_water_layer(j) = recruitw*ccohort_hydr%l_aroot_layer(j)/sum_l_aroot
    end do

    do j=1,csite_hydr%nlevrhiz
       watres_local = csite_hydr%wrf_soil(j)%p%th_from_psi(bc_in%smpmin_si*denh2o*grav_earth*m_per_mm*mpa_per_pa)

       total_water = sum(csite_hydr%v_shell(j,:)*csite_hydr%h2osoi_liqvol_shell(j,:))
       total_water_min = sum(csite_hydr%v_shell(j,:)*watres_local)

       !assumes that only 50% is available for recruit water....
       recruit_water_avail_layer(j)=0.5_r8*max(0.0_r8,total_water-total_water_min)

    end do

    nmin  = 1.0e+36 
    do j=1,csite_hydr%nlevrhiz
       if(cohort_recruit_water_layer(j)>0.0_r8) then
          n = recruit_water_avail_layer(j)/cohort_recruit_water_layer(j)
          nmin = min(n, nmin) 
       endif
    end do
    ccohort%n = min (ccohort%n, nmin) 

  end subroutine ConstrainRecruitNumber


  ! =====================================================================================

  subroutine SavePreviousRhizVolumes(currentSite)

    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(ed_site_hydr_type), pointer    :: csite_hydr

    csite_hydr => currentSite%si_hydr
    csite_hydr%l_aroot_layer_init(:)  = csite_hydr%l_aroot_layer(:)
    csite_hydr%r_node_shell_init(:,:) = csite_hydr%r_node_shell(:,:)
    csite_hydr%v_shell_init(:,:)      = csite_hydr%v_shell(:,:)

    return
  end subroutine SavePreviousRhizVolumes

  ! ======================================================================================

  subroutine UpdateSizeDepRhizVolLenCon(currentSite, bc_in)

    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:


    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in) :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(ed_site_hydr_type), pointer :: csite_hydr
    type(ed_patch_type)  , pointer :: cPatch
    type(ed_cohort_type) , pointer :: cCohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    real(r8)                       :: hksat_s                      ! hksat converted to units of 10^6sec 
                                                                   ! which is equiv to       [kg m-1 s-1 MPa-1]
    integer                        :: j,k                          ! gridcell, soil layer, rhizosphere shell indices
    integer                        :: j_bc                         ! soil layer index of boundary condition
    real(r8)                       :: large_kmax_bound = 1.e4_r8   ! for replacing kmax_bound_shell wherever the 
    ! innermost shell radius is less than the assumed 
    ! absorbing root radius rs1
    ! 1.e-5_r8 from Rudinger et al 1994             	
    integer                        :: nlevrhiz
    integer, parameter :: k_inner = 1   ! innermost rhizosphere shell
    !-----------------------------------------------------------------------

    csite_hydr => currentSite%si_hydr
    nlevrhiz = csite_hydr%nlevrhiz

    ! update cohort-level root length density and accumulate it across cohorts and patches to the column level
    csite_hydr%l_aroot_layer(:)  = 0._r8
    cPatch => currentSite%youngest_patch
    do while(associated(cPatch))
       cCohort => cPatch%tallest
       do while(associated(cCohort))
          ccohort_hydr => cCohort%co_hydr
          csite_hydr%l_aroot_layer(:) = csite_hydr%l_aroot_layer(:) + ccohort_hydr%l_aroot_layer(:)*cCohort%n
          cCohort => cCohort%shorter
       enddo !cohort
       cPatch => cPatch%older
    enddo !patch

    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    do j = 1,nlevrhiz
       ! proceed only if l_aroot_coh has changed
       !       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          call shellGeom( csite_hydr%l_aroot_layer(j), csite_hydr%rs1(j), AREA, csite_hydr%dz_rhiz(j), &
               csite_hydr%r_out_shell(j,:), csite_hydr%r_node_shell(j,:),csite_hydr%v_shell(j,:))
!       end if !has l_aroot_layer changed?
    enddo


    do j = 1,nlevrhiz
       j_bc = j+csite_hydr%i_rhiz_t-1

       ! bc_in%hksat_sisl(j): hydraulic conductivity at saturation (mm H2O /s)
       !
       ! converted from [mm H2O s-1] -> [kg s-1 MPa-1 m-1]
       !
       ! Conversion of Pascals:    1 Pa = 1 kg m-1 s-2
       !
       ! [mm s-1] * 1e-3 [m mm-1] 
       !          * 1 [kg m-1 s-2 Pa-1] 
       !          * 9.8-1 [s2 m-1] 
       !          * 1e6 [Pa MPa-1]
       !                           = [kg s-1 m-1 MPa-1]

       hksat_s = bc_in%hksat_sisl(j_bc) * m_per_mm * 1._r8/grav_earth * pa_per_mpa
       
       ! proceed only if the total absorbing root length (site-level) has changed in this layer
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then

          ! Set the max conductance on the inner shell first.  If the node radius
          ! on the shell is smaller than the root radius, just set the max conductance
          ! to something extremely high.
              
           if( csite_hydr%r_node_shell(j,k_inner) <= csite_hydr%rs1(j) ) then
               csite_hydr%kmax_upper_shell(j,k_inner) = large_kmax_bound
           else
               csite_hydr%kmax_upper_shell(j,k_inner) = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                     log(csite_hydr%r_node_shell(j,k_inner)/csite_hydr%rs1(j))*hksat_s
           end if
           
           csite_hydr%kmax_lower_shell(j,k_inner) = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                 log(csite_hydr%r_out_shell(j,k_inner)/csite_hydr%r_node_shell(j,k_inner) )*hksat_s
           
           do k = 2,nshell
               csite_hydr%kmax_upper_shell(j,k)        = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                     log(csite_hydr%r_node_shell(j,k)/csite_hydr%r_out_shell(j,k-1))*hksat_s
               
               csite_hydr%kmax_lower_shell(j,k)        = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                     log(csite_hydr%r_out_shell(j,k)/csite_hydr%r_node_shell(j,k  ))*hksat_s
           enddo ! loop over rhizosphere shells
           
           


       end if !has l_aroot_layer changed?
    enddo ! loop over soil layers

    return
  end subroutine UpdateSizeDepRhizVolLenCon


  ! =====================================================================================


  subroutine UpdateSizeDepRhizHydProps(currentSite, bc_in )
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:

    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in) :: bc_in


    ! Save current volumes, lenghts and nodes to an "initial"
    ! used to calculate effects in states later on.

    call SavePreviousRhizVolumes(currentSite)

    ! Update the properties of the vegetation-soil hydraulic environment
    ! these are independent on the water state

    call UpdateSizeDepRhizVolLenCon(currentSite, bc_in)


    return
  end subroutine UpdateSizeDepRhizHydProps

  ! =================================================================================

  subroutine UpdateSizeDepRhizHydStates(currentSite, bc_in)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type), intent(in)              :: bc_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: v_rhiz(nlevsoi_hyd_max)                  ! updated volume of all rhizosphere compartments              [m3]
    real(r8) :: r_delta                                  ! change in radius of innermost rhizosphere compartment       [m]
    real(r8) :: dpsidr                                   ! water potential gradient near root surface                  [MPa/m]
    real(r8) :: w_shell_new                              ! updated water volume in rhizosphere compartment             [m3]
    real(r8) :: w_layer_init(nlevsoi_hyd_max)            ! initial water mass by layer                                 [kg]
    real(r8) :: w_layer_interp(nlevsoi_hyd_max)          ! water mass after interpolating to new rhizosphere           [kg]
    real(r8) :: w_layer_new(nlevsoi_hyd_max)             ! water mass by layer after interpolation and fudging         [kg]
    real(r8) :: h2osoi_liq_col_new(nlevsoi_hyd_max)      ! water mass per area after interpolating to new rhizosphere  [kg/m2]
    real(r8) :: s_shell_init(nlevsoi_hyd_max,nshell)     ! initial saturation fraction in rhizosphere compartment      [0-1]
    real(r8) :: s_shell_interp(nlevsoi_hyd_max,nshell)   ! interpolated saturation fraction in rhizosphere compartment [0-1]
    real(r8) :: psi_shell_init(nlevsoi_hyd_max,nshell)   ! initial water potential in rhizosphere compartment          [MPa]
    real(r8) :: psi_shell_interp(nlevsoi_hyd_max,nshell) ! interpolated psi_shell to new r_node_shell                  [MPa]
    real(r8) :: delta_s(nlevsoi_hyd_max)                 ! change in saturation fraction needed to ensure water bal    [0-1]
    real(r8) :: errh2o(nlevsoi_hyd_max)                  ! water budget error after updating                           [kg/m2]
    integer  :: j,k                                      ! gridcell, column, soil layer, rhizosphere shell indicies
    integer  :: j_bc                                     ! level index for boundary conditions
    integer  :: indexc,indexj                            ! column and layer indices where there is a water balance error
    logical  :: found                                    ! flag in search loop
    type(ed_site_hydr_type), pointer :: csite_hydr
    !-----------------------------------------------------------------------

    s_shell_init(:,:)         = 0._r8
    psi_shell_init(:,:)       = 0._r8
    psi_shell_interp(:,:)     = 0._r8
    s_shell_interp(:,:)       = 0._r8

    csite_hydr => currentSite%si_hydr

    if(.false.) then

       do j = 1, csite_hydr%nlevrhiz
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then

             do k = 1,nshell
                psi_shell_init(j,k) = csite_hydr%wrf_soil(j)%p%psi_from_th(csite_hydr%h2osoi_liqvol_shell(j,k)) 
             end do

          end if !has l_aroot_coh changed?
       enddo

       ! interpolate initial psi values by layer and shell
       ! BOC...To-Do: need to constrain psi to be within realistic limits (i.e., < 0)
       do j = 1,csite_hydr%nlevrhiz
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then

             ! fine root length increased, thus shrinking the rhizosphere size
             if(csite_hydr%r_node_shell(j,nshell) < csite_hydr%r_node_shell_init(j,nshell)) then
                r_delta               = csite_hydr%r_node_shell(j,1) - csite_hydr%r_node_shell_init(j,1)
                !dpsidr                = (psi_shell_init(j,2) - psi_shell_init(j,1)) / &
                !                        (csite_hydr%r_node_shell_init(j,2) - csite_hydr%r_node_shell_init(j,1))

                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
                ! HACK for special case of nshell = 1 -- compiler throws error because of index 2 in above line, 
                ! even though at run-time the code should skip over this section: MUST FIX
                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                dpsidr                = (psi_shell_init(j,1) - psi_shell_init(j,1)) / &
                     (csite_hydr%r_node_shell_init(j,1) - csite_hydr%r_node_shell_init(j,1))
                psi_shell_interp(j,1) = dpsidr * r_delta
                do k = 2,nshell
                   r_delta               = csite_hydr%r_node_shell(j,k) - csite_hydr%r_node_shell_init(j,k)
                   dpsidr                = (psi_shell_init(j,k) - psi_shell_init(j,k-1)) / &
                        (csite_hydr%r_node_shell_init(j,k) - csite_hydr%r_node_shell_init(j,k-1))
                   psi_shell_interp(j,k) = dpsidr * r_delta
                enddo
             else                       
                ! fine root length decreased, thus increasing the rhizosphere size
                do k = 1,(nshell-1)
                   r_delta               = csite_hydr%r_node_shell(j,k) - csite_hydr%r_node_shell_init(j,k)
                   dpsidr                = (psi_shell_init(j,k+1) - psi_shell_init(j,k)) / &
                        (csite_hydr%r_node_shell_init(j,k+1) - csite_hydr%r_node_shell_init(j,k))
                   psi_shell_interp(j,k) = dpsidr * r_delta
                enddo
                r_delta               = csite_hydr%r_node_shell(j,nshell) - csite_hydr%r_node_shell_init(j,nshell)
                !dpsidr                = (psi_shell_init(j,nshell) - psi_shell_init(j,nshell-1)) / &
                !                        (csite_hydr%r_node_shell_init(j,nshell) - csite_hydr%r_node_shell_init(j,nshell-1))

                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                ! HACK for special case of nshell = 1 -- compiler throws error because of index nshell-1 in 
                ! above line, even though at run-time the code should skip over this section: MUST FIX
                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                dpsidr                = (psi_shell_init(j,nshell) - psi_shell_init(j,nshell)) / &
                     (csite_hydr%r_node_shell_init(j,nshell) - csite_hydr%r_node_shell_init(j,nshell))

                psi_shell_interp(j,k) = dpsidr * r_delta
             end if
          end if !has l_aroot_coh changed?
       enddo

       ! 1st guess at new s based on interpolated psi
       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1
                 
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
             
             s_shell_interp(j,k) = ( csite_hydr%wrf_soil(j)%p%th_from_psi(psi_shell_interp(j,k)) - bc_in%watres_sisl(j_bc)) / & 
                  (bc_in%watres_sisl(j_bc)+bc_in%watres_sisl(j_bc))

          end if !has l_aroot_coh changed?
       enddo

       ! accumlate water across shells for each layer (initial and interpolated)
       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
             w_layer_init(j)      = 0._r8
             w_layer_interp(j)    = 0._r8
             v_rhiz(j)            = 0._r8
             do k = 1,nshell
                w_layer_init(j)   = w_layer_init(j) + denh2o * &
                     (csite_hydr%v_shell_init(j,k)*csite_hydr%h2osoi_liqvol_shell(j,k) )
                w_layer_interp(j) = w_layer_interp(j) + denh2o * & 
                     (csite_hydr%v_shell(j,k) * &
                     (s_shell_interp(j,k)*(bc_in%watsat_sisl(j_bc)-bc_in%watres_sisl(j_bc))+bc_in%watres_sisl(j_bc)) )
                v_rhiz(j)         = v_rhiz(j) + csite_hydr%v_shell(j,k)
             enddo
          end if !has l_aroot_coh changed?
       enddo

       ! estimate delta_s across all shells needed to ensure total water in each layer doesn't change
       ! BOC...FIX: need to handle special cases where delta_s causes s_shell to go above or below 1 or 0, respectively.
       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
             delta_s(j) = (( w_layer_init(j) - w_layer_interp(j) )/( v_rhiz(j) * denh2o ) - bc_in%watres_sisl(j_bc)) / &
                  (bc_in%watsat_sisl(j_bc)-bc_in%watres_sisl(j_bc))
          end if !has l_aroot_coh changed?
       enddo

       ! update h2osoi_liqvol_shell and h2osoi_liq_shell
       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1
          ! proceed only if l_aroot_coh has changed
          if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
             w_layer_new(j)                = 0._r8
             do k = 1,nshell
                s_shell_interp(j,k)        = s_shell_interp(j,k) + delta_s(j)
                csite_hydr%h2osoi_liqvol_shell(j,k) = s_shell_interp(j,k) * &
                     ( bc_in%watsat_sisl(j_bc)-bc_in%watres_sisl(j_bc) ) + bc_in%watres_sisl(j_bc)
                w_shell_new                = csite_hydr%h2osoi_liqvol_shell(j,k) * &
                     csite_hydr%v_shell(j,k)
                w_layer_new(j)             = w_layer_new(j) + w_shell_new
             enddo
             h2osoi_liq_col_new(j)         = w_layer_new(j)/ v_rhiz(j)
          end if !has l_aroot_coh changed?
       enddo

       ! balance check
       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1
          errh2o(j) = h2osoi_liq_col_new(j) - bc_in%h2o_liq_sisl(j_bc)
          if (abs(errh2o(j)) > 1.e-4_r8) then
             write(fates_log(),*)'WARNING:  water balance error ',&
                  ' updating rhizosphere shells: ',j,errh2o(j)
             write(fates_log(),*)'errh2o= ',errh2o(j), ' [kg/m2]'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       enddo

    end if !nshell > 1

  end subroutine UpdateSizeDepRhizHydStates

  ! ====================================================================================
  
  subroutine BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)

    ! Arguments
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    ! Locals
    integer                                 :: s
    integer                                 :: ifp
    real(r8)                                :: balive_patch
    type(ed_patch_type),pointer             :: cpatch 
    type(ed_cohort_type),pointer            :: ccohort

    do s = 1,nsites

       ifp = 0
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))                 
          ifp=ifp+1

          balive_patch = 0._r8
          ccohort=>cpatch%tallest
          do while(associated(ccohort))
             balive_patch = balive_patch +  &
                  (cCohort%prt%GetState(fnrt_organ, all_carbon_elements) + &
                  cCohort%prt%GetState(sapw_organ, all_carbon_elements) + &
                  cCohort%prt%GetState(leaf_organ, all_carbon_elements))* ccohort%n
             ccohort => ccohort%shorter
          enddo !cohort

          bc_out(s)%btran_pa(ifp) = 0.0_r8
          ccohort=>cpatch%tallest
          do while(associated(ccohort))
             bc_out(s)%btran_pa(ifp) =  bc_out(s)%btran_pa(ifp) + &
                  ccohort%co_hydr%btran * &
                  (cCohort%prt%GetState(fnrt_organ, all_carbon_elements) + &
                  cCohort%prt%GetState(sapw_organ, all_carbon_elements) + &
                  cCohort%prt%GetState(leaf_organ, all_carbon_elements)) * &
                  ccohort%n / balive_patch
             ccohort => ccohort%shorter
          enddo !cohort
          cpatch => cpatch%younger
       enddo !end patch loop
    end do
    return
  end subroutine BTranForHLMDiagnosticsFromCohortHydr

  ! ==========================================================================

  subroutine FillDrainRhizShells(nsites, sites, bc_in, bc_out)
    !
    ! Created by Brad Christoffersen, Jan 2016
    !
    ! !DESCRIPTION:
    ! Parses out mean vertical water fluxes resulting from infiltration,
    ! drainage, and vertical water movement (dwat_kgm2) over radially stratified
    ! rhizosphere shells.
    !
    ! The approach used is heuristic, but based on the principle that water
    ! fluxing out of a layer will preferentially come from rhizosphere
    ! shells with higher water contents/potentials within that layer, and
    ! alternatively, that water fluxing into a layer will preferentially go
    ! into shells with lower water contents/potentials.
    !
    ! This principle is implemented by filling (draining) the rhizosphere
    ! shells in order from the driest (wettest) shell to the wettest (driest).
    ! Each shell is filled (drained) up (down) to the next wettest (driest)
    ! shell until the change in mean layer water (dwat_kgm2) is accounted for.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in)                       :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type), intent(in)              :: bc_in(nsites)
    type(bc_out_type), intent(inout)          :: bc_out(nsites)

    ! Locals
    type(ed_site_hydr_type), pointer :: csite_hydr       ! pointer to site hydraulics object
    real(r8) :: dwat_kgm2                                ! change in layer water content              [kg/m2]
    integer  :: s,j,k                                    ! site, soil layer, rhizosphere shell indicies
    integer  :: i,f,ff,kk                                ! indicies
    integer  :: j_bc                                     ! layer index for matching boundary condition soil layers
    integer  :: indexj                                   ! column and layer indices where there is a water balance error
    integer  :: ordered(nshell) = (/(i,i=1,nshell,1)/)   ! array of rhizosphere indices which have been ordered
    real(r8) :: area_col                                 ! column area                                                    [m2]
    real(r8) :: v_cum                                    ! cumulative shell volume from driest/wettest shell to kth shell [m3]
    real(r8) :: dwat_kg                                  ! water remaining to be distributed across shells                [kg]
    real(r8) :: thdiff                                   ! water content difference between ordered adjacent rhiz shells  [m3 m-3]
    real(r8) :: wdiff                                    ! mass of water represented by thdiff over previous k shells     [kg]
    real(r8) :: errh2o(nlevsoi_hyd_max)                  ! water budget error after updating                              [kg/m2]
    real(r8) :: cumShellH2O                              ! sum of water in all the shells of a specific layer             [kg/m2]
    real(r8) :: h2osoi_liq_shell(nlevsoi_hyd_max,nshell) ! water in the rhizosphere shells                               [kg]
    integer  :: tmp                                      ! temporary
    logical  :: found                                    ! flag in search loop
    !-----------------------------------------------------------------------

    do s = 1,nsites


       ! First step, identify how the liquid water in each layer has changed
       ! since the last time it was updated. This should be due to drainage.
       ! The drainage component should be the total change in liquid water content from the last time
       ! the hydraulics driver was called, and then adding back in the losses due to root uptake
       ! (which was already taken out).

       ! BOC: This was previously in HydrologyDrainage:

       csite_hydr => sites(s)%si_hydr

       do j = 1,csite_hydr%nlevrhiz
          j_bc = j+csite_hydr%i_rhiz_t-1

          cumShellH2O=sum(csite_hydr%h2osoi_liqvol_shell(j,:) *csite_hydr%v_shell(j,:)) * denh2o*AREA_INV
          
          dwat_kgm2 = bc_in(s)%h2o_liq_sisl(j_bc) - cumShellH2O
             
          dwat_kg = dwat_kgm2 * AREA
          
          ! order shells in terms of increasing or decreasing volumetric water content
          ! algorithm same as that used in histFileMod.F90 to alphabetize history tape contents
          if(nshell > 1) then
             do k = nshell-1,1,-1
                do kk = 1,k
                   if (csite_hydr%h2osoi_liqvol_shell(j,ordered(kk)) > &
                        csite_hydr%h2osoi_liqvol_shell(j,ordered(kk+1))) then
                      if (dwat_kg > 0._r8) then  !order increasing
                         tmp           = ordered(kk)
                         ordered(kk)   = ordered(kk+1)
                         ordered(kk+1) = tmp
                      end if
                   else
                      if (dwat_kg < 0._r8) then  !order decreasing
                         tmp           = ordered(kk)
                         ordered(kk)   = ordered(kk+1)
                         ordered(kk+1) = tmp
                      end if
                   end if
                enddo
             enddo
          end if
          
          ! fill shells with water up to the water content of the next-wettest shell, 
          ! in order from driest to wettest (dwat_kg > 0)
          ! ------ OR ------
          ! drain shells' water down to the water content of the next-driest shell, 
          ! in order from wettest to driest (dwat_kg < 0)
          k = 1
          do while ( (dwat_kg /= 0._r8) .and. (k < nshell) )
             thdiff = csite_hydr%h2osoi_liqvol_shell(j,ordered(k+1)) - &
                  csite_hydr%h2osoi_liqvol_shell(j,ordered(k))
             v_cum  = sum(csite_hydr%v_shell(j,ordered(1:k)))
             wdiff  = thdiff * v_cum * denh2o        ! change in h2o [kg / ha] for shells ordered(1:k)
             if(abs(dwat_kg) >= abs(wdiff)) then
                csite_hydr%h2osoi_liqvol_shell(j,ordered(1:k)) = csite_hydr%h2osoi_liqvol_shell(j,ordered(k+1))
                dwat_kg  = dwat_kg - wdiff
             else
                csite_hydr%h2osoi_liqvol_shell(j,ordered(1:k)) = &
                     csite_hydr%h2osoi_liqvol_shell(j,ordered(1:k)) + dwat_kg/denh2o/v_cum
                dwat_kg  = 0._r8
             end if
             k = k + 1
          enddo
          
          if (dwat_kg /= 0._r8) then
             v_cum  = sum(csite_hydr%v_shell(j,ordered(1:nshell)))
             thdiff = dwat_kg / v_cum / denh2o
             do k = nshell, 1, -1
                csite_hydr%h2osoi_liqvol_shell(j,k) = csite_hydr%h2osoi_liqvol_shell(j,k) + thdiff
             end do
          end if
          
          ! m3/m3 * Total volume m3 * kg/m3 = kg
          h2osoi_liq_shell(j,:) = csite_hydr%h2osoi_liqvol_shell(j,:) * &
               csite_hydr%v_shell(j,:) * denh2o
          
          
          errh2o(j) = sum(h2osoi_liq_shell(j,:))*AREA_INV - bc_in(s)%h2o_liq_sisl(j_bc)
          
          if (abs(errh2o(j)) > 1.e-9_r8) then
             write(fates_log(),*)'WARNING:  water balance error in FillDrainRhizShells'
             write(fates_log(),*)'errh2o= ',errh2o(j), ' [kg/m2]'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end do

    end do
    return
  end subroutine FillDrainRhizShells

  ! ====================================================================================

  subroutine hydraulics_bc ( nsites, sites, bc_in, bc_out, dtime)

    ! ----------------------------------------------------------------------------------
    ! added by Brad Christoffersen Jan 2016 for use in ED hydraulics
    !    van Genuchten (1980)-specific functions for the swc (soil water characteristic)
    !    and for the kunsat (unsaturated hydraulic conductivity) curves. Test mod 06/20/2016

    ! resolved the mass-balance bugs and tested Jan, 2018 by C. XU
    !
    ! BOC...for quick implementation avoided JT's abstract interface,
    !    but these should be converted to interfaces in the future
    ! ----------------------------------------------------------------------------------

    !
    ! !DESCRIPTION:
    !s
    ! !USES:
    use FatesUtilsMod  , only : check_var_real

    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    real(r8),intent(in)                     :: dtime

    !
    ! !LOCAL VARIABLES:
    integer :: iv    ! leaf layer
    integer :: ifp   ! index of FATES patch
    integer :: s     ! index of FATES site
    integer :: i     ! shell index
    integer :: j,jj  ! soil layer
    integer :: j_bc  ! soil layer index for boundary conditions
    integer :: k     ! 1D plant-soil continuum array
    integer :: ft    ! plant functional type index
    integer :: sz    ! plant's size class index
    integer :: t     ! previous timesteps (for lwp stability calculation)
    integer :: nstep !number of time steps

    !----------------------------------------------------------------------

    type (ed_patch_type),  pointer     :: cpatch       ! current patch pointer
    type (ed_cohort_type), pointer     :: ccohort      ! current cohort pointer
    type(ed_site_hydr_type), pointer   :: site_hydr    ! site hydraulics pointer
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr ! cohort hydraulics pointer

    ! Local arrays

    ! accumulated water content change over all cohorts in a column   [m3 m-3]
    real(r8)            :: dth_layershell_col(nlevsoi_hyd_max,nshell)

    ! array of soil layer indices which have been ordered
    integer             :: ordered(nlevsoi_hyd_max) = (/(j,j=1,nlevsoi_hyd_max,1)/) 

    ! total absorbing root & rhizosphere conductance (over all shells) by soil layer   [MPa]
    real(r8) :: kbg_layer(nlevsoi_hyd_max)   
    real(r8) :: rootuptake(nlevsoi_hyd_max) ! mass-flux from 1st rhizosphere to absorbing roots            [kg/indiv/layer/step]
    
    real(r8) :: site_runoff         ! If plants are pushing water into saturated soils, we create
                                    ! runoff. This is either banked, or sent to the correct flux pool [kg/m2]
    real(r8) :: aroot_frac_plant    ! The fraction of the total length of absorbing roots contained in one soil layer
                                    ! that are devoted to a single plant
    real(r8) :: wb_err_plant        ! Solve error for a single plant [kg]
    real(r8) :: wb_check_site       ! the water balance error we get from summing fluxes
                                    ! and changes in storage
                                    ! and is just a double check on our error accounting). [kg/m2]
    real(r8) :: dwat_plant          ! change in water mass in the whole plant [kg]
    real(r8) :: qflx_tran_veg_indiv ! individiual transpiration rate [kgh2o indiv-1 s-1]
    real(r8) :: gscan_patch         ! sum of ccohort%gscan across all cohorts within a patch          
    real(r8) :: sapflow             ! mass-flux for the cohort between transporting root and stem  [kg/indiv/step]
    real(r8) :: prev_h2oveg         ! plant water storage at start of timestep (kg/m2)
    real(r8) :: prev_h2osoil        ! soil water storage at start of timestep (kg/m2)
    logical  :: recruitflag         ! flag to check if there is newly recruited cohorts
    real(r8) :: root_flux           ! total water flux into roots [kg/m2]
    real(r8) :: transp_flux         ! total transpiration flux from plants [kg/m2]
    real(r8) :: delta_plant_storage ! change in plant water storage over the step [kg/m2]
    real(r8) :: delta_soil_storage  ! change in soil water storage over the step [kg/m2]
    real(r8) :: sumcheck            ! used to debug mass balance in soil horizon diagnostics
    integer  :: nlevrhiz            ! local for number of rhizosphere levels
    integer  :: sc                  ! size class index
    
    ! ----------------------------------------------------------------------------------
    ! Important note: We are interested in calculating the total fluxes in and out of the
    ! site/column.  Usually, when we do things like this, we acknowledge that FATES
    ! does not consider the bare ground patch.  However, since this routine
    ! calculates "column level" fluxes, we have to factor in that patch-level fluxes
    ! are only accounting for a portion of the area.
    ! ----------------------------------------------------------------------------------

    !For newly recruited cohorts, add the water uptake demand to csite_hydr%recruit_w_uptake
    call RecruitWUptake(nsites,sites,bc_in,dtime,recruitflag)

    !update water storage in veg after incorporating newly recuited cohorts
    if(recruitflag) call UpdateH2OVeg(nsites,sites,bc_out)

    do s = 1, nsites

       site_hydr => sites(s)%si_hydr

       nlevrhiz = site_hydr%nlevrhiz
       
       ! AVERAGE ROOT WATER UPTAKE (BY RHIZOSPHERE SHELL) ACROSS ALL COHORTS WITHIN A COLUMN
       dth_layershell_col(:,:)  = 0._r8
       site_hydr%dwat_veg       = 0._r8
       site_hydr%errh2o_hyd     = 0._r8
       prev_h2oveg    = site_hydr%h2oveg
       prev_h2osoil   = sum(site_hydr%h2osoi_liqvol_shell(:,:) * & 
                        site_hydr%v_shell(:,:)) * denh2o * AREA_INV

       bc_out(s)%qflx_ro_sisl(:) = 0._r8

       ! Zero out diagnotsics that rely on accumulation
       site_hydr%sapflow_scpf(:,:)       = 0._r8
       site_hydr%rootuptake_sl(:)        = 0._r8
       site_hydr%rootuptake0_scpf(:,:)   = 0._r8
       site_hydr%rootuptake10_scpf(:,:)  = 0._r8
       site_hydr%rootuptake50_scpf(:,:)  = 0._r8
       site_hydr%rootuptake100_scpf(:,:) = 0._r8

       
       ! Initialize water mass balancing terms [kg h2o / m2]
       ! --------------------------------------------------------------------------------
       transp_flux          = 0._r8
       root_flux            = 0._r8
       
       ! Initialize the delta in soil water and plant water storage
       ! with the initial condition.
       
       !err_soil = delta_soil_storage - root_flux
       !err_plot = delta_plant_storage - (root_flux - transp_flux)

       ifp = 0
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          ifp = ifp + 1

          ! ----------------------------------------------------------------------------
          ! Objective: Partition the transpiration flux
          ! specfied by the land model to the cohorts. The weighting
          ! factor we use to downscale is the cohort combo term: g_sb_laweight
          ! This term is the stomatal conductance multiplied by total leaf
          ! area.  gscan_patch is the sum over all cohorts, used to normalize.
          ! ----------------------------------------------------------------------------

          gscan_patch   = 0.0_r8
          ccohort=>cpatch%tallest
          do while(associated(ccohort))
             ccohort_hydr => ccohort%co_hydr
             gscan_patch       = gscan_patch + ccohort%g_sb_laweight
             ccohort => ccohort%shorter
          enddo !cohort

          ! The HLM predicted transpiration flux even though no leaves are present?
          if(bc_in(s)%qflx_transp_pa(ifp) > 1.e-10_r8 .and. gscan_patch<nearzero)then
             write(fates_log(),*) 'ERROR in plant hydraulics.'
             write(fates_log(),*) 'The HLM predicted a non-zero total transpiration flux'
             write(fates_log(),*) 'for this patch, yet there is no leaf-area-weighted conductance?'
             write(fates_log(),*) 'transp: ',bc_in(s)%qflx_transp_pa(ifp)
             write(fates_log(),*) 'gscan_patch: ',gscan_patch
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          ccohort=>cpatch%tallest
          do while(associated(ccohort))

             ccohort_hydr => ccohort%co_hydr
             ft       = ccohort%pft

             ! Relative transpiration of this cohort from the whole patch
             ! Note that g_sb_laweight / gscan_patch is the weighting that gives cohort contribution per area
             ! [mm H2O/plant/s]  = [mm H2O/ m2 / s] * [m2 / patch] * [cohort/plant] * [patch/cohort]

             if(ccohort%g_sb_laweight>nearzero) then
                qflx_tran_veg_indiv     = bc_in(s)%qflx_transp_pa(ifp) * cpatch%total_canopy_area * & 
                (ccohort%g_sb_laweight/gscan_patch)/ccohort%n
             else
                qflx_tran_veg_indiv     = 0._r8
             end if

             ! Save the transpiration flux for diagnostics (currently its a constant boundary condition)
             ccohort_hydr%qtop = qflx_tran_veg_indiv*dtime

             transp_flux = transp_flux + (qflx_tran_veg_indiv*dtime)*ccohort%n*AREA_INV

             ! VERTICAL LAYER CONTRIBUTION TO TOTAL ROOT WATER UPTAKE OR LOSS
             !    _____ 
             !   |     |
             !   |leaf |
             !   |_____|
             !      /
             !      \
             !      /
             !    __\__ 
             !   |     |
             !   |stem |
             !   |_____|
             !------/----------------_____---------------------------------  
             !      \               |     |   |    |      |       |        |
             !      /          _/\/\|aroot|   |    |shell | shell | shell  |               layer j-1
             !      \        _/     |_____|   |    | k-1  |   k   |  k+1   |
             !------/------_/--------_____--------------------------------------
             !      \    _/         |     |    |     |       |        |         |
             !    __/__ / _/\/\/\/\/|aroot|    |     | shell | shell  | shell   |          layer j
             !   |     |_/          |_____|    |     |  k-1  |   k    |  k+1    |
             !---|troot|-------------_____----------------------------------------------
             !   |_____|\_          |     |      |      |        |          |           |
             !            \/\/\/\/\/|aroot|      |      | shell  |  shell   |   shell   |  layer j+1
             !                      |_____|      |      |  k-1   |    k     |    k+1    |
             !---------------------------------------------------------------------------


             ! This routine will update the theta values for 1 cohort's flow-path
             ! from leaf to the current soil layer.  This does NOT
             ! update cohort%th_*

             if(use_2d_hydrosolve) then

                 call MatSolve2D(site_hydr,ccohort,ccohort_hydr, &
                       dtime,qflx_tran_veg_indiv, &
                       sapflow,rootuptake(1:nlevrhiz),wb_err_plant,dwat_plant, &
                       dth_layershell_col)
                
             else

                ! ---------------------------------------------------------------------------------
                ! Approach: do nlevsoi_hyd sequential solutions to Richards' equation,
                !           each of which encompass all plant nodes and soil nodes for a given soil layer j,
                !           with the timestep fraction for each layer-specific solution proportional to each 
                ! layer's contribution to the total root-soil conductance
                ! Water potential in plant nodes is updated after each solution
                ! As such, the order across soil layers in which the solution is conducted matters.
                ! For now, the order proceeds across soil layers in order of decreasing root-soil conductance
                ! NET EFFECT: total water removed from plant-soil system remains the same: it 
                !             sums up to total transpiration (qflx_tran_veg_indiv*dtime)
                !             root water uptake in each layer is proportional to each layer's total 
                !             root length density and soil matric potential
                !             root hydraulic redistribution emerges within this sequence when a 
                !             layers have transporting-to-absorbing root water potential gradients of opposite sign
                ! -----------------------------------------------------------------------------------
                
                call OrderLayersForSolve1D(site_hydr, ccohort, ccohort_hydr, ordered, kbg_layer)
                
                call ImTaylorSolve1D(site_hydr,ccohort,ccohort_hydr, &
                                     dtime,qflx_tran_veg_indiv,ordered, kbg_layer, & 
                                     sapflow,rootuptake(1:nlevrhiz), & 
                                     wb_err_plant,dwat_plant, &
                                     dth_layershell_col)

             end if

             ! Remember the error for the cohort
             ccohort_hydr%errh2o  = ccohort_hydr%errh2o + wb_err_plant
             
             ! Update total error in [kg/m2 ground]
             site_hydr%errh2o_hyd = site_hydr%errh2o_hyd + wb_err_plant*ccohort%n*AREA_INV

             ! Accumulate site level diagnostic of plant water change [kg/m2]
             ! (this is zerod)
             site_hydr%dwat_veg   = site_hydr%dwat_veg + dwat_plant*ccohort%n*AREA_INV
        
             ! Update total site-level stored plant water [kg/m2]
             ! (this is not zerod, but incremented)
             site_hydr%h2oveg     = site_hydr%h2oveg + dwat_plant*ccohort%n*AREA_INV

             sc = ccohort%size_class
             
             ! Sapflow diagnostic [kg/ha/s]
             site_hydr%sapflow_scpf(sc,ft) = site_hydr%sapflow_scpf(sc,ft) + sapflow*ccohort%n/dtime

             ! Root uptake per rhiz layer [kg/ha/s]
             site_hydr%rootuptake_sl(1:nlevrhiz) = site_hydr%rootuptake_sl(1:nlevrhiz) + &
                  rootuptake(1:nlevrhiz)*ccohort%n/dtime

             ! Root uptake per pft x size class, over set layer depths [kg/ha/m/s]
             ! These are normalized by depth (in case the desired horizon extends
             ! beyond the actual rhizosphere)

             site_hydr%rootuptake0_scpf(sc,ft) = site_hydr%rootuptake0_scpf(sc,ft) + & 
                  SumBetweenDepths(site_hydr,0._r8,0.1_r8,rootuptake(1:nlevrhiz))*ccohort%n/dtime

             site_hydr%rootuptake10_scpf(sc,ft) = site_hydr%rootuptake10_scpf(sc,ft) + & 
                  SumBetweenDepths(site_hydr,0.1_r8,0.5_r8,rootuptake(1:nlevrhiz))*ccohort%n/dtime

             site_hydr%rootuptake50_scpf(sc,ft) = site_hydr%rootuptake50_scpf(sc,ft) + & 
                  SumBetweenDepths(site_hydr,0.5_r8,1.0_r8,rootuptake(1:nlevrhiz))*ccohort%n/dtime

             site_hydr%rootuptake100_scpf(sc,ft) = site_hydr%rootuptake100_scpf(sc,ft) + & 
                  SumBetweenDepths(site_hydr,1.0_r8,1.e10_r8,rootuptake(1:nlevrhiz))*ccohort%n/dtime
             
             ! ---------------------------------------------------------
             ! Update water potential and frac total conductivity
             ! of plant compartments
             ! ---------------------------------------------------------
             
             call UpdatePlantPsiFTCFromTheta(ccohort,site_hydr)

             ccohort_hydr%btran = wkf_plant(stomata_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(1))
             

             ccohort => ccohort%shorter
          enddo !cohort 

          cpatch => cpatch%younger
      enddo !patch

       ! --------------------------------------------------------------------------------
       ! The cohort level water fluxes are complete, the remainder of this subroutine
       ! is dedicated to doing site level resulting mass balance calculations and checks
       ! --------------------------------------------------------------------------------

       ! Calculate the amount of water fluxing through the roots. It is the sum 
       ! of the change in thr rhizosphere shells.  Note that following this calculation
       ! we may adjust the change in soil water to avoid super-saturation and sub-residual
       ! water contents.  But the pre-adjusted value is the actual amount of root flux.
       ! [kg/m2]

       root_flux = -sum(dth_layershell_col(1:site_hydr%nlevrhiz,:)*site_hydr%v_shell(:,:))*denh2o*AREA_INV

       
       do j=1,site_hydr%nlevrhiz
           j_bc = j+site_hydr%i_rhiz_t-1
          
           ! Update the site-level state variable 
           ! rhizosphere shell water content [m3/m3]
           site_hydr%h2osoi_liqvol_shell(j,:) =  site_hydr%h2osoi_liqvol_shell(j,:) + &
                 dth_layershell_col(j,:)


           bc_out(s)%qflx_soil2root_sisl(j_bc) = &
                 -(sum(dth_layershell_col(j,:)*site_hydr%v_shell(j,:))*denh2o*AREA_INV/dtime) + &
                 site_hydr%recruit_w_uptake(j)
           
           
          ! Save the amount of liquid soil water known to the model after root uptake
          ! This calculation also assumes that 1mm of water is 1kg
          site_hydr%h2osoi_liq_prev(j) = bc_in(s)%h2o_liq_sisl(j_bc) - &
               dtime*bc_out(s)%qflx_soil2root_sisl(j_bc)


          ! We accept that it is possible for gravity to push
          ! water into saturated soils, particularly at night when
          ! transpiration has stopped. In the real world, the water
          ! would be driven out of the layer, although we have no
          ! boundary flux on the rhizospheres in these substeps. To accomodate
          ! this, if soils are pushed beyond saturation minus a small buffer
          ! then we remove that excess, send it to a runoff pool, and 
          ! fix the node's water content to the saturation minus buffer value

          site_runoff          = 0._r8
          if(purge_supersaturation) then
             do i = 1,nshell
                if(site_hydr%h2osoi_liqvol_shell(j,i)>(bc_in(s)%watsat_sisl(j_bc)-thsat_buff)) then
                   
                   ! [m3/m3] * [kg/m3] * [m3/site] * [site/m2] => [kg/m2]
                   site_runoff = site_runoff + & 
                        (site_hydr%h2osoi_liqvol_shell(j,i)-(bc_in(s)%watsat_sisl(j_bc)-thsat_buff)) * &
                        site_hydr%v_shell(j,i)*AREA_INV*denh2o
                   
                   site_hydr%h2osoi_liqvol_shell(j,i) = bc_in(s)%watsat_sisl(j_bc)-thsat_buff
                   
                end if
             end do
             
             bc_out(s)%qflx_ro_sisl(j_bc) = site_runoff/dtime
          end if
      enddo

       
      ! Note that the cohort-level solvers are expected to update
      ! site_hydr%h2oveg

      ! Calculate site total kg's of runoff
      site_runoff = sum(bc_out(s)%qflx_ro_sisl(:))*dtime
      
      delta_plant_storage = site_hydr%h2oveg - prev_h2oveg
       
      delta_soil_storage  = sum(site_hydr%h2osoi_liqvol_shell(:,:) * & 
            site_hydr%v_shell(:,:)) * denh2o * AREA_INV - prev_h2osoil
       
      if(abs(delta_plant_storage - (root_flux - transp_flux)) > 1.e-6_r8 ) then
          write(fates_log(),*) 'Site plant water balance does not close'
          write(fates_log(),*) 'delta plant storage: ',delta_plant_storage,' [kg/m2]'
          write(fates_log(),*) 'integrated root flux: ',root_flux,' [kg/m2]'
          write(fates_log(),*) 'transpiration flux: ',transp_flux,' [kg/m2]'
          write(fates_log(),*) 'end storage: ',site_hydr%h2oveg
          call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      
       if(abs(delta_soil_storage + root_flux + site_runoff) > 1.e-6_r8 ) then
          write(fates_log(),*) 'Site soil water balance does not close'
          write(fates_log(),*) 'delta soil storage: ',delta_soil_storage,' [kg/m2]'
          write(fates_log(),*) 'integrated root flux (pos into root): ',root_flux,' [kg/m2]'
          write(fates_log(),*) 'site runoff: ',site_runoff,' [kg/m2]'
          write(fates_log(),*) 'end storage: ',sum(site_hydr%h2osoi_liqvol_shell(:,:) * & 
               site_hydr%v_shell(:,:)) * denh2o * AREA_INV, & 
               ' [kg/m2]'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if


       !-----------------------------------------------------------------------
       ! mass balance check and pass the total stored vegetation water to HLM
       ! in order for it to fill its balance checks


       ! Compare the integrated error to the site mass balance
       ! error sign is positive towards transpiration overestimation
       ! Loss fluxes should = decrease in storage
       ! (transp_flux + site_runoff) =  -(delta_plant_storage+delta_soil_storage )

       wb_check_site = delta_plant_storage+delta_soil_storage+site_runoff+transp_flux

       if( abs(wb_check_site - site_hydr%errh2o_hyd) > 1.e-10_r8 ) then
           write(fates_log(),*) 'FATES hydro water ERROR balance does not add up [kg/m2]'
           write(fates_log(),*) 'wb_error_site: ',site_hydr%errh2o_hyd
           write(fates_log(),*) 'wb_check_site: ',wb_check_site
           write(fates_log(),*) 'delta_plant_storage: ',delta_plant_storage
           write(fates_log(),*) 'delta_soil_storage: ',delta_soil_storage
           write(fates_log(),*) 'site_runoff: ',site_runoff
           write(fates_log(),*) 'transp_flux: ',transp_flux
       end if

       ! Now check on total error
       if( abs(wb_check_site) > 1.e-6_r8 ) then
           write(fates_log(),*) 'FATES hydro water balance does not add up [kg/m2]'
           write(fates_log(),*) 'site_hydr%errh2o_hyd: ',wb_check_site
           write(fates_log(),*) 'delta_plant_storage: ',delta_plant_storage
           write(fates_log(),*) 'delta_soil_storage: ',delta_soil_storage
           write(fates_log(),*) 'site_runoff: ',site_runoff
           write(fates_log(),*) 'transp_flux: ',transp_flux
       end if


       site_hydr%h2oveg_hydro_err = site_hydr%h2oveg_hydro_err + site_hydr%errh2o_hyd

       bc_out(s)%plant_stored_h2o_si = site_hydr%h2oveg + site_hydr%h2oveg_dead - &
            site_hydr%h2oveg_growturn_err - &
            site_hydr%h2oveg_pheno_err-&
            site_hydr%h2oveg_hydro_err

    enddo !site

    return
  end subroutine Hydraulics_BC

  ! =====================================================================================



  subroutine UpdatePlantKmax(ccohort_hydr,ccohort,csite_hydr)

    ! ---------------------------------------------------------------------------------
    !
    ! This routine sets the maximum conductance of all compartments in the plant, from
    ! leaves, to stem, to transporting root, to the absorbing roots.
    ! These properties are dependent only on the materials (conductivity) and the
    ! geometry of the compartments.
    ! The units of all K_max values are [kg H2O s-1 MPa-1]
    !
    ! There are some different ways to represent overall conductance from node-to-node
    ! throughout the hydraulic system. Universally, all can make use of a system
    ! where we separate the hydraulic compartments of the nodes into the upper (closer
    ! to the sky) and lower (away from the sky) portions of the compartment. It is
    ! possible that due to things like xylem taper, the two portions may have different
    ! conductivity, and therefore differnet conductances.
    !
    ! Assumption 0.  This routine calculates maximum conductivity for 1 plant.
    ! Assumption 1.  The compartment volumes, heights and lengths have all been 
    !                determined, probably called just before this routine.
    !
    ! Steudle, E. Water uptake by roots: effects of water deficit. 
    ! J Exp Bot 51, 1531-1542, doi:DOI 10.1093/jexbot/51.350.1531 (2000).
    ! ---------------------------------------------------------------------------------

    ! Arguments

    type(ed_cohort_hydr_type),intent(inout),target :: ccohort_hydr
    type(ed_cohort_type),intent(in),target         :: ccohort
    type(ed_site_hydr_type),intent(in),target      :: csite_hydr

    ! Locals
    integer :: k                     ! Compartment (node) index
    integer :: j                     ! Soil layer index
    integer :: k_ag                  ! Compartment index for above-ground indexed array
    integer  :: pft                  ! Plant Functional Type index
    real(r8) :: c_sap_dummy          ! Dummy variable (unused) with sapwood carbon [kg]
    real(r8) :: z_lower              ! distance between lower edge and mean petiole height [m]
    real(r8) :: z_upper              ! distance between upper edge and mean petiole height [m]
    real(r8) :: z_node               ! distance between compartment center and mph [m]
    real(r8) :: kmax_lower           ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
    real(r8) :: kmax_node            ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
    real(r8) :: kmax_upper           ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
    real(r8) :: a_sapwood            ! Mean cross section area of sapwood   [m2]
    real(r8) :: rmin_ag              ! Minimum total resistance of all above ground pathways
                                     ! [kg-1 s MPa]
    real(r8) :: kmax_bg              ! Total maximum conductance of all below-ground pathways 
                                     ! from the absorbing roots center nodes to the 
                                     ! transporting root center node
    real(r8) :: rootfr               ! fraction of absorbing root in each soil layer
                                     ! assumes propotion of absorbing root is equal
                                     ! to proportion of total root
    real(r8) :: kmax_layer           ! max conductance between transporting root node
                                     ! and absorbing root node in each layer [kg s-1 MPa-1]
    real(r8) :: surfarea_aroot_layer ! Surface area of absorbing roots in each
                                     ! soil layer [m2]
    real(r8) :: sum_l_aroot          ! sum of plant's total root length
    real(r8),parameter :: taper_exponent = 1._r8/3._r8 ! Savage et al. (2010) xylem taper exponent [-]
    real(r8),parameter :: min_pet_stem_dz = 0.00001_r8  ! Force at least a small difference
                                                       ! in the top of stem and petiole


    pft   = ccohort%pft

    ! Get the cross-section of the plant's sapwood area [m2]
    call bsap_allom(ccohort%dbh,pft,ccohort%canopy_trim,a_sapwood,c_sap_dummy)

    ! Leaf Maximum Hydraulic Conductance
    ! The starting hypothesis is that there is no resistance inside the
    ! leaf, between the petiole and the center of storage.  To override
    ! this, make provisions by changing the kmax to a not-absurdly high 
    ! value.  It is assumed that the conductance in this default case,
    ! is regulated completely by the stem conductance from the stem's
    ! center of storage, to the petiole.

    ccohort_hydr%kmax_petiole_to_leaf = 1.e8_r8


    ! Stem Maximum Hydraulic Conductance

    do k=1, n_hypool_stem

       ! index for "above-ground" arrays, that contain stem and leaf
       ! in one vector
       k_ag = k+n_hypool_leaf

       ! Depth from the petiole to the lower, node and upper compartment edges

       z_lower = ccohort_hydr%z_node_ag(n_hypool_leaf) - ccohort_hydr%z_lower_ag(k_ag)
       z_node  = ccohort_hydr%z_node_ag(n_hypool_leaf) - ccohort_hydr%z_node_ag(k_ag)
       z_upper = max( min_pet_stem_dz,ccohort_hydr%z_node_ag(n_hypool_leaf) - & 
             ccohort_hydr%z_upper_ag(k_ag))


       ! Then we calculate the maximum conductance from each the lower, node and upper 
       ! edges of the compartment to the petiole. The xylem taper factor requires
       ! that the kmax it is scaling is from the point of interest to the mean height
       ! of the petioles.  Then we can back out the conductance over just the path
       ! of the upper and lower compartments, but subtracting them as resistors in
       ! series.

       ! max conductance from upper edge to mean petiole height
       ! If there is no height difference between the upper compartment edge and
       ! the petiole, at least give it some nominal amount to void FPE's
       kmax_upper = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
            xylemtaper(taper_exponent, z_upper) * &
            a_sapwood / z_upper

       ! max conductance from node to mean petiole height
       kmax_node  = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
            xylemtaper(taper_exponent, z_node) * &
            a_sapwood / z_node

       ! max conductance from lower edge to mean petiole height
       kmax_lower = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
            xylemtaper(taper_exponent, z_lower) * &
            a_sapwood / z_lower

       ! Max conductance over the path of the upper side of the compartment
       ccohort_hydr%kmax_stem_upper(k) = (1._r8/kmax_node - 1._r8/kmax_upper)**(-1._r8)

       ! Max conductance over the path on the loewr side of the compartment
       ccohort_hydr%kmax_stem_lower(k) = (1._r8/kmax_lower - 1._r8/kmax_node)**(-1._r8)
       
       if(debug) then
           ! The following clauses should never be true:
           if( (z_lower < z_node) .or. & 
                 (z_node  < z_upper) ) then
               write(fates_log(),*) 'Problem calculating stem Kmax'
               write(fates_log(),*) z_lower, z_node, z_upper
               write(fates_log(),*) kmax_lower*z_lower, kmax_node*z_node, kmax_upper*z_upper
               call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
       end if
       
    enddo

    ! Maximum conductance of the upper compartment in the transporting root
    ! that connects to the lowest stem (btw: z_lower_ag(n_hypool_ag) == 0)

    z_upper = ccohort_hydr%z_lower_ag(n_hypool_leaf)
    z_node  = ccohort_hydr%z_lower_ag(n_hypool_leaf)-ccohort_hydr%z_node_troot

    kmax_node = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
         xylemtaper(taper_exponent, z_node) * &
         a_sapwood / z_node

    kmax_upper = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
         xylemtaper(taper_exponent, z_upper) * &
         a_sapwood / z_upper

    ccohort_hydr%kmax_troot_upper = (1._r8/kmax_node - 1._r8/kmax_upper)**(-1._r8)

    !print*,z_upper,z_node,kmax_upper,kmax_node,ccohort_hydr%kmax_troot_upper


    ! The maximum conductance between the center node of the transporting root
    ! compartment, and the center node of the absorbing root compartment, is calculated
    ! as a residual.  Specifically, we look at the total resistance the plant has in
    ! the stem so far, by adding those resistances in series.
    ! Then we use a parameter to specify what fraction of the resistance
    ! should be below-ground between the transporting root node and the absorbing roots.
    ! After that total is calculated, we then convert to a conductance, and split the
    ! conductance in parallel between root layers, based on the root fraction.
    ! Note* The inverse of max conductance (KMax) is minimum resistance:


    rmin_ag = 1._r8/ccohort_hydr%kmax_petiole_to_leaf + &
         sum(1._r8/ccohort_hydr%kmax_stem_upper(1:n_hypool_stem)) + &
         sum(1._r8/ccohort_hydr%kmax_stem_lower(1:n_hypool_stem)) + &
         1._r8/ccohort_hydr%kmax_troot_upper

    ! Calculate the residual resistance below ground, as a resistor
    ! in series with the existing above ground
    ! Invert to find below-ground kmax
    ! (rmin_ag+rmin_bg)*fr = rmin_ag
    ! rmin_ag + rmin_bg = rmin_ag/fr
    ! rmin_bg = (1/fr-1) * rmin_ag
    !
    ! if kmax_bg = 1/rmin_bg :
    !
    ! kmax_bg = 1/((1/fr-1) * rmin_ag)
    
    kmax_bg = 1._r8/(rmin_ag*(1._r8/EDPftvarcon_inst%hydr_rfrac_stem(pft) - 1._r8))
    

    ! The max conductance of each layer is in parallel, therefore
    ! the kmax terms of each layer, should sum to kmax_bg
    sum_l_aroot = sum(ccohort_hydr%l_aroot_layer(:))
    do j=1,csite_hydr%nlevrhiz

       kmax_layer = kmax_bg*ccohort_hydr%l_aroot_layer(j)/sum_l_aroot

       ! Two transport pathways, in two compartments exist in each layer.
       ! These pathways are connected in serial.
       ! For simplicity, we simply split the resistance between the two.
       ! Mathematically, this results in simply doubling the conductance
       ! and applying to both paths.  Here are the two paths:
       ! 1) is the path between the transporting root's center node, to
       !    the boundary of the transporting root with the boundary of
       !    the absorbing root  (kmax_troot_lower)
       ! 2) is the path between the boundary of the absorbing root and
       !    transporting root, with the absorbing root's center node
       !    (kmax_aroot_upper)

       ccohort_hydr%kmax_troot_lower(j) = 3.0_r8 * kmax_layer
       ccohort_hydr%kmax_aroot_upper(j) = 3.0_r8 * kmax_layer
       ccohort_hydr%kmax_aroot_lower(j) = 3.0_r8 * kmax_layer

    end do

    ! Finally, we calculate maximum radial conductance from the root
    ! surface to its center node.  This transport is not a xylem transport
    ! like the calculations prior to this. This transport is through the
    ! exodermis, cortex, casparian strip and endodermis.  The actual conductance
    ! will possibly depend on the potential gradient (whether out-of the root,
    ! or in-to the root).  So we calculate the kmax's for both cases,
    ! and save them for the final conductance calculation.

    do j=1,csite_hydr%nlevrhiz

       ! Surface area of the absorbing roots for a single plant in this layer [m2]
       surfarea_aroot_layer = 2._r8 * pi_const * & 
            EDPftvarcon_inst%hydr_rs2(ccohort%pft) * ccohort_hydr%l_aroot_layer(j)

       ! Convert from surface conductivity [kg H2O m-2 s-1 MPa-1] to [kg H2O s-1 MPa-1]
       ccohort_hydr%kmax_aroot_radial_in(j) = hydr_kmax_rsurf1 * surfarea_aroot_layer

       ccohort_hydr%kmax_aroot_radial_out(j) = hydr_kmax_rsurf2 * surfarea_aroot_layer

    end do

    return
  end subroutine UpdatePlantKmax

  ! ===================================================================================

  subroutine OrderLayersForSolve1D(site_hydr,cohort,cohort_hydr,ordered, kbg_layer)
    
    ! Arguments (IN)
    type(ed_site_hydr_type), intent(in),target   :: site_hydr
    type(ed_cohort_type), intent(in),target      :: cohort
    type(ed_cohort_hydr_type),intent(in),target  :: cohort_hydr


    ! Arguments (INOUT)
    integer, intent(inout)                       :: ordered(:)
    real(r8), intent(out)                        :: kbg_layer(:)
    
    ! Locals
    
    real(r8) :: kbg_tot                    ! total absorbing root & rhizosphere conductance (over all shells and soil layers  [MPa]
    real(r8) :: psi_inner_shell            ! matric potential of the inner shell, used for calculating
                                           ! which kmax to use when forecasting uptake layer ordering [MPa]
    real(r8) :: psi_aroot                  ! matric potential of absorbing root [MPa]
    real(r8) :: kmax_aroot                 ! max conductance of the absorbing root [kg s-1 Mpa-1]
    real(r8) :: ftc_aroot                  ! fraction of total conductivity of abs root
    real(r8) :: r_bg                       ! total estimated resistance in below ground compartments
                                           ! for each soil layer [s Mpa kg-1] (used to predict order in 1d solve)
    real(r8) :: aroot_frac_plant           ! This is the fraction of absorbing root from one plant
    real(r8) :: kmax_lo                    ! maximum conductance of lower (away from atm) half of path [kg s-1 Mpa-1]
    real(r8) :: kmax_up                    ! maximum conductance of upper (close to atm)  half of path [kg s-1 MPa-1]
    real(r8) :: psi_shell                  ! matric potential of a given shell [-]
    real(r8) :: ftc_shell                  ! fraction of total cond. of a given rhiz shell [-]
    integer  :: tmp                        ! temporarily holds a soil layer index
    integer  :: ft                         ! functional type index of plant
    integer  :: j,jj,k                     ! layer and shell indices
    
    
    kbg_tot      = 0._r8
    kbg_layer(:) = 0._r8

    ft = cohort%pft
    
    do j=1,site_hydr%nlevrhiz

       ! Path is between the absorbing root
       ! and the first rhizosphere shell nodes
       ! Special case. Maximum conductance depends on the 
       ! potential gradient (same elevation, no geopotential
       ! required.
       
       psi_inner_shell = site_hydr%wrf_soil(j)%p%psi_from_th(site_hydr%h2osoi_liqvol_shell(j,1)) 
       
       ! Note, since their is no elevation difference between
       ! the absorbing root and its layer, no need to calc
       ! diff in total, just matric is fine [MPa]
       if(cohort_hydr%psi_aroot(j) < psi_inner_shell) then
          kmax_aroot = cohort_hydr%kmax_aroot_radial_in(j)
       else
          kmax_aroot = cohort_hydr%kmax_aroot_radial_out(j)
       end if
       
       ! Get matric potential [Mpa] of the absorbing root
       psi_aroot = wrf_plant(aroot_p_media,ft)%p%psi_from_th(cohort_hydr%th_aroot(j))
       
       ! Get Fraction of Total Conductivity [-] of the absorbing root
       ftc_aroot = wkf_plant(aroot_p_media,ft)%p%ftc_from_psi(cohort_hydr%psi_aroot(j))

       ! Calculate total effective conductance over path  [kg s-1 MPa-1]
       ! from absorbing root node to 1st rhizosphere shell
       r_bg = 1._r8/(kmax_aroot*ftc_aroot)
       
       ! Path is across the upper an lower rhizosphere comparment
       ! on each side of the nodes. Since there is no flow across the outer
       ! node to the edge, we ignore that last half compartment
       aroot_frac_plant = cohort_hydr%l_aroot_layer(j)/site_hydr%l_aroot_layer(j)
       
       do k = 1,nshell
          
          kmax_up = site_hydr%kmax_upper_shell(j,k)*aroot_frac_plant
          kmax_lo = site_hydr%kmax_lower_shell(j,k)*aroot_frac_plant
          
          psi_shell = site_hydr%wrf_soil(j)%p%psi_from_th(site_hydr%h2osoi_liqvol_shell(j,k))
          
          ftc_shell = site_hydr%wkf_soil(j)%p%ftc_from_psi(psi_shell)
          
          r_bg = r_bg + 1._r8/(kmax_up*ftc_shell)
          if(k<nshell) r_bg = r_bg + 1._r8/(kmax_lo*ftc_shell )
       end do
       
       !! upper bound limited to size()-1 b/c of zero-flux outer boundary condition
       kbg_layer(j)        = 1._r8/r_bg
       kbg_tot             = kbg_tot + kbg_layer(j)

    enddo !soil layer

    
    kbg_layer = kbg_layer/kbg_tot

    ! order soil layers in terms of decreasing volumetric water content
    ! algorithm same as that used in histFileMod.F90 to alphabetize history tape contents
    do j = site_hydr%nlevrhiz-1,1,-1
       do jj = 1,j
          if (kbg_layer(ordered(jj)) <= kbg_layer(ordered(jj+1))) then
             tmp           = ordered(jj)
             ordered(jj)   = ordered(jj+1)
             ordered(jj+1) = tmp
          end if
       enddo
    enddo

    
    return
  end subroutine OrderLayersForSolve1D
  
  ! =================================================================================

  subroutine ImTaylorSolve1D(site_hydr,cohort,cohort_hydr,dtime,q_top, &
       ordered,kbg_layer, sapflow,rootuptake,&
       wb_err_plant,dwat_plant,dth_layershell_col)

    ! -------------------------------------------------------------------------------
    ! Calculate the hydraulic conductances across a list of paths.  The list is a 1D vector, and
    ! the list need not be across the whole path from stomata to the last rhizosphere shell, but
    ! it can only be 1d, which is part of a path through the plant and into 1 soil layer.
    !
    ! Note on conventions:
    ! "Up" upper, refers to the compartment that is closer to the atmosphere
    ! "lo" lower, refers to the compartment that is further from the atmosphere
    ! Weird distinction: since flow from one node to another, will include half of
    ! a compartment on a upper node, and half a compartment of a lower node.  The upp
    ! compartment will be contributing its lower compartment, and the lower node
    ! will be presenting it upper compartment. Yes, confusing, but non-the-less 
    ! accurate.
    ! -------------------------------------------------------------------------------

    ! Arguments (IN)
    type(ed_cohort_type),intent(in),target       :: cohort
    type(ed_cohort_hydr_type),intent(inout),target  :: cohort_hydr
    type(ed_site_hydr_type), intent(in),target   :: site_hydr
    real(r8), intent(in)                         :: dtime
    real(r8), intent(in)                         :: q_top        ! transpiration flux rate at upper boundary [kg -s]
    integer,intent(in)                           :: ordered(:)   ! Layer solution order
    real(r8), intent(in)                         :: kbg_layer(:) ! relative conductance of each layer

    ! Arguments (OUT)

    real(r8),intent(out) :: sapflow                   ! time integrated mass flux between transp-root and stem [kg]
    real(r8),intent(out) :: rootuptake(:)             ! time integrated mass flux between rhizosphere and aroot [kg]
    real(r8),intent(out) :: wb_err_plant              ! total error from the plant, transpiration
                                                      ! should match change in storage [kg]
    real(r8),intent(out) :: dwat_plant                ! Change in plant stored water [kg]
    real(r8),intent(inout) :: dth_layershell_col(:,:) ! accumulated water content change over all cohorts in a column   [m3 m-3])

    ! Locals
    integer :: i              ! node index "i"
    integer :: j              ! path index "j"
    integer :: jj             ! alt path index
    integer :: nsteps         ! number of sub-steps in any given iteration loop, starts at 1 and grows
    integer :: ilayer         ! soil layer index of interest
    integer :: itest          ! node index used for testing and reporting errors
    integer :: ishell         ! rhizosphere shell index of the node
    integer :: ishell_up      ! rhizosphere shell index on the upstream side of flow path (towards soil)
    integer :: ishell_dn      ! rhizosphere shell index on the downstream side of flow path (towards atm)
    integer :: i_up           ! node index on the upstream side of flow path (towards soil)
    integer :: i_dn           ! node index on the downstream side of flow path (towards atm)
    integer :: istep          ! sub-step count index
    integer :: tri_ierr       ! error flag for the tri-diagonal solver 0=passed, 1=failed
    logical :: solution_found ! logical set to true if a solution was found within error tolerance
    real(r8) :: dt_step       ! time [seconds] over-which to calculate solution
    real(r8) :: q_top_eff     ! effective water flux through stomata [kg s-1 plant-1]
    real(r8) :: rootfr_scaler ! Factor to scale down cross-section areas based on what
                              ! fraction of root is in current layer [-]
    real(r8) :: kmax_dn       ! maximum conductance of downstream half of path [kg s-1 Mpa-1]
    real(r8) :: kmax_up       ! maximum conductance of upstream half of path [kg s-1 MPa-1]
    real(r8) :: wb_step_err   ! water balance error over substep [kg]
    real(r8) :: w_tot_beg     ! total plant water prior to solve [kg]
    real(r8) :: w_tot_end     ! total plant water at end of solve [kg]
    real(r8) :: dt_substep    ! timestep length of substeps [s]
    real(r8) :: leaf_water    ! kg of water in the leaf
    real(r8) :: stem_water    ! kg of water in the stem
    real(r8) :: root_water    ! kg of water in the transp and absorbing roots
    real(r8) :: sapflow_lyr   ! sapflow flux [kg] per layer per timestep
    real(r8) :: rootuptake_lyr! rootuptake flux [kg] per layer per timestep 
    real(r8) :: wb_err_layer                ! balance error for the layer [kg/cohort]
    

    real(r8) :: dth_node(n_hypool_tot)          ! change in theta over the timestep
    real(r8) :: th_node_init(n_hypool_tot)      ! "theta" i.e. water content of node [m3 m-3]
                                                ! before the solve
    real(r8) :: th_node(n_hypool_tot)           ! "theta" during the solve (dynamic) [m3 m-3]
    real(r8) :: z_node(n_hypool_tot)            ! elevation of node [m]
    real(r8) :: v_node(n_hypool_tot)            ! volume of the node, ie single plant compartments [m3]
    real(r8) :: psi_node(n_hypool_tot)          ! matric potential on node [Mpa]
    real(r8) :: ftc_node(n_hypool_tot)          ! frac total conductance on node [-]
    real(r8) :: h_node(n_hypool_tot)            ! total potential on node [Mpa]
    real(r8) :: error_arr(n_hypool_tot)         ! array that saves problematic diagnostics for reporting
    real(r8) :: dftc_dtheta_node(n_hypool_tot)  ! deriv FTC w.r.t. theta
    real(r8) :: dpsi_dtheta_node(n_hypool_tot)  ! deriv psi w.r.t. theta
    real(r8) :: k_eff(n_hypool_tot-1)           ! effective (used) conductance over path [kg s-1 MPa-1]
    real(r8) :: a_term(n_hypool_tot-1)          ! "A" term in the tri-diagonal implicit solve [-]
    real(r8) :: b_term(n_hypool_tot-1)          ! "B" term in the tri-diagonal implicit solve [-]
    real(r8) :: k_diag(n_hypool_tot-1)          ! mean time-averaged K over the paths (diagnostic) [kg s-1 Mpa-1]
    real(r8) :: flux_diag(n_hypool_tot-1)       ! time-integrated mass flux over sub-steps [kg]
    real(r8) :: h_diag, psi_diag                ! total and matric potential for error reporting [Mpa]
    real(r8) :: tris_a(n_hypool_tot)            ! left of diagonal terms for tri-diagonal matrix solving delta theta
    real(r8) :: tris_b(n_hypool_tot)            ! center diagonal terms for tri-diagonal matrix solving delta theta
    real(r8) :: tris_c(n_hypool_tot)            ! right of diaongal terms for tri-diagonal matrix solving delta theta
    real(r8) :: tris_r(n_hypool_tot)            ! off (constant coefficients) matrix terms
    real(r8) :: sum_l_aroot                     !
    real(r8) :: aroot_frac_plant                ! This is the fraction of absorbing root from one plant
    real(r8) :: dftc_dpsi                       ! Change in fraction of total conductance wrt change
                                                ! in potential [- MPa-1]
    integer  :: error_code                      ! flag that specifies which check tripped a failed solution
    integer  :: ft                              ! plant functional type
    real(r8) :: q_flow                          ! flow diagnostic [kg]
    real(r8) :: rootfr                          ! rooting fraction of this layer (used for diagnostics)
    ! out of the total absorbing roots from the whole community of plants
    integer  :: iter                      ! iteration count for sub-step loops

    integer, parameter  :: imult    = 3                ! With each iteration, increase the number of substeps
                                                       ! by this much
    integer, parameter  :: max_iter = 20               ! Maximum number of iterations with which we reduce timestep
   
    real(r8), parameter :: max_wb_err      = 1.e-5_r8  ! threshold for water balance error (stop model)   [kg h2o]


    logical, parameter :: no_ftc_radialk = .false.
    logical, parameter :: weight_serial_dt = .true. ! if this is true, and we are not doing spatial parallelism
                                                    ! then we give the fraction of time as a function of how
                                                    ! much conductance the layer has

    associate(pm_node => site_hydr%pm_node)

    ! This is the maximum number of iterations needed for this cohort
    ! (each soil layer has a different number, this saves the max)
    cohort_hydr%iterh1 = 0
    cohort_hydr%iterh2 = 0
    
    ! Initialize plant water error (integrated flux-storage)
    wb_err_plant = 0._r8 

    ! Initialize integrated change in total plant water
    dwat_plant = 0._r8
    
    ! These are diagnostics that must be calculated.
    ! in this routine (uses differentials and actual fluxes)
    ! So we need to zero them, as they are incremented
    ! over the sub-steps
    sapflow = 0._r8
    rootuptake(:) = 0._r8
    
    ft = cohort%pft

    ! Total length of roots per plant for this cohort
    sum_l_aroot = sum(cohort_hydr%l_aroot_layer(:))
    
    ! -----------------------------------------------------------------------------------
    ! As mentioned when calling this routine, we calculate a solution to the flux
    ! equations, sequentially, for the plant and each soil layer.
    ! Go through soil layers in order of decreasing total root-soil conductance
    ! -----------------------------------------------------------------------------------

    do jj=1,site_hydr%nlevrhiz
        
        ilayer = ordered(jj)

        if(do_parallel_stem) then
            ! If we do "parallel" stem
            ! conduits, we integrate
            ! each layer over the whole time, but
            ! reduce the conductance cross section
            ! according to what fraction of root is active
            dt_step = dtime
        else
            if(weight_serial_dt)then
                dt_step = dtime*kbg_layer(ilayer)
            else
                dt_step = dtime/real(site_hydr%nlevrhiz,r8)
            end if
        end if
        
        ! -------------------------------------------------------------------------------
        ! Part 1.  Calculate node quantities:
        !          matric potential: psi_node
        !          fraction of total conductance: ftc_node
        !          total potential (matric + elevatio) h_node
        !          deriv. ftc  wrt  theta: dftc_dtheta_node
        !          deriv. psi  wrt  theta: dpsi_dtheta_node
        ! -------------------------------------------------------------------------------
        
        
        ! This is the fraction of total absorbing root length that a single
        ! plant for this cohort takes up, relative to ALL cohorts at the site. Note:
        ! cohort_hydr%l_aroot_layer(ilayer) is units [m/plant]
        ! site_hydr%l_aroot_layer(ilayer) is units [m/site]
        
        aroot_frac_plant = cohort_hydr%l_aroot_layer(ilayer)/site_hydr%l_aroot_layer(ilayer)
        
        wb_err_layer = 0._r8

        ! If in "spatially parallel" mode, scale down cross section
        ! of flux through top by the root fraction of this layer

        if(do_parallel_stem)then
           rootfr_scaler = cohort_hydr%l_aroot_layer(ilayer)/sum_l_aroot
        else
           rootfr_scaler = 1.0_r8
        end if

        q_top_eff = q_top * rootfr_scaler

        ! For all nodes leaf through rhizosphere
        ! Send node heights and compartment volumes to a node-based array
        
        do i = 1,n_hypool_tot
            
            if (i<=n_hypool_ag) then
                z_node(i)  = cohort_hydr%z_node_ag(i)
                v_node(i)  = cohort_hydr%v_ag(i)
                th_node_init(i) = cohort_hydr%th_ag(i)
            elseif (i==n_hypool_ag+1) then
                z_node(i)  = cohort_hydr%z_node_troot
                v_node(i)  = cohort_hydr%v_troot
                th_node_init(i) = cohort_hydr%th_troot
            elseif (i==n_hypool_ag+2) then
                z_node(i)  = -site_hydr%zi_rhiz(ilayer)
                v_node(i)  = cohort_hydr%v_aroot_layer(ilayer)
                th_node_init(i) = cohort_hydr%th_aroot(ilayer)
            else
                ishell  = i-(n_hypool_ag+2)
                z_node(i)  = -site_hydr%zi_rhiz(ilayer)
                ! The volume of the Rhizosphere for a single plant
                v_node(i)  = site_hydr%v_shell(ilayer,ishell)*aroot_frac_plant
                th_node_init(i) = site_hydr%h2osoi_liqvol_shell(ilayer,ishell)
            end if
        end do
        
        ! Outer iteration loop
        ! This cuts timestep in half and resolve the solution with smaller substeps
        ! This loop is cleared when the model has found a solution
        
        solution_found = .false.
        iter = 0
        do while( .not.solution_found ) 

            ! Gracefully quit if too many iterations have been used
            if(iter>max_iter)then
                call Report1DError(cohort,site_hydr,ilayer,z_node,v_node, & 
                      th_node_init,q_top_eff,dt_step,w_tot_beg,w_tot_end,& 
                      rootfr_scaler,aroot_frac_plant,error_code,error_arr)

                call endrun(msg=errMsg(sourcefile, __LINE__))
            end if

            ! If debugging, then lets re-initialize our diagnostics of
            ! time integrated K and flux across the paths
            if(debug)then
                k_diag    = 0._r8
                flux_diag = 0._r8
            end if

            sapflow_lyr    = 0._r8
            rootuptake_lyr = 0._r8
            
            ! For each attempt, we want to reset theta with the initial value
            th_node(:) = th_node_init(:)
            
            ! Determine how many substeps, and how long they are

            nsteps = max(imult*iter,1)   ! Factor by which we divide through the timestep
                                         ! start with full step (ie dt_fac = 1)
                                         ! Then increase per the "imult" value.
            
            dt_substep = dt_step/real(nsteps,r8) ! This is the sub-stem length in seconds
            
            ! Walk through sub-steps
            do istep = 1,nsteps

                ! Total water mass in the plant at the beginning of this solve [kg h2o]
                w_tot_beg = sum(th_node(:)*v_node(:))*denh2o

                ! Calculate on-node quantities: potential, and derivatives
                do i = 1,n_hypool_plant

                    ! Get matric potential [Mpa]
                    psi_node(i) = wrf_plant(pm_node(i),ft)%p%psi_from_th(th_node(i))

                    ! Get total potential [Mpa]
                    h_node(i) =  mpa_per_pa*denh2o*grav_earth*z_node(i) + psi_node(i)

                    ! Get Fraction of Total Conductivity [-]
                    ftc_node(i) = wkf_plant(pm_node(i),ft)%p%ftc_from_psi(psi_node(i))

                    ! deriv psi wrt theta
                    dpsi_dtheta_node(i) = wrf_plant(pm_node(i),ft)%p%dpsidth_from_th(th_node(i))

                    ! deriv ftc wrt psi

                    dftc_dpsi = wkf_plant(pm_node(i),ft)%p%dftcdpsi_from_psi(psi_node(i))

                    dftc_dtheta_node(i) = dftc_dpsi * dpsi_dtheta_node(i) 

                    ! We have two ways to calculate radial absorbing root conductance
                    ! 1) Assume that water potential does not effect conductance
                    ! 2) The standard FTC function applies

                    if(i==n_hypool_ag+2)then
                        if(no_ftc_radialk) then
                            ftc_node(i)         = 1.0_r8
                            dftc_dtheta_node(i) = 0.0_r8
                        end if
                    end if

                 end do


                ! Same updates as loop above, but for rhizosphere shells
                
                do i = n_hypool_plant+1,n_hypool_tot
                    psi_node(i)         = site_hydr%wrf_soil(ilayer)%p%psi_from_th(th_node(i))
                    h_node(i)           = mpa_per_pa*denh2o*grav_earth*z_node(i) + psi_node(i)
                    ftc_node(i)         = site_hydr%wkf_soil(ilayer)%p%ftc_from_psi(psi_node(i))
                    dpsi_dtheta_node(i) = site_hydr%wrf_soil(ilayer)%p%dpsidth_from_th(th_node(i))
                    dftc_dpsi           = site_hydr%wkf_soil(ilayer)%p%dftcdpsi_from_psi(psi_node(i))
                    dftc_dtheta_node(i) = dftc_dpsi * dpsi_dtheta_node(i) 
                end do

                !--------------------------------------------------------------------------------
                ! Part 2.  Effective conductances over the path-length and Flux terms
                !          over the node-to-node paths
                !--------------------------------------------------------------------------------
                
                ! Path is between the leaf node and first stem node
                ! -------------------------------------------------------------------------------
                
                j        = 1
                i_up     = 2   ! upstream node index
                i_dn     = 1   ! downstream node index
                kmax_dn  = rootfr_scaler*cohort_hydr%kmax_petiole_to_leaf
                kmax_up  = rootfr_scaler*cohort_hydr%kmax_stem_upper(1)

                call GetImTaylorKAB(kmax_up,kmax_dn,        &
                      ftc_node(i_up),ftc_node(i_dn),        & 
                      h_node(i_up),h_node(i_dn),            & 
                      dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                      dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                      k_eff(j),                         &
                      A_term(j),                        & 
                      B_term(j))

                
                ! Path is between stem nodes
                ! -------------------------------------------------------------------------------
                
                do j=2,n_hypool_ag-1

                    i_up = j+1
                    i_dn = j

                    ! "Up" is the "upstream" node, which also uses
                    ! the "upper" side of its compartment for the calculation.
                    ! "dn" is the "downstream" node, which uses the lower
                    ! side of its compartment
                    ! This compartment is the "lower" node, but uses
                    ! the "higher" side of its compartment.

                    kmax_dn  = rootfr_scaler*cohort_hydr%kmax_stem_lower(i_dn-n_hypool_leaf)
                    kmax_up  = rootfr_scaler*cohort_hydr%kmax_stem_upper(i_up-n_hypool_leaf)

                    call GetImTaylorKAB(kmax_up,kmax_dn,        &
                          ftc_node(i_up),ftc_node(i_dn),        & 
                          h_node(i_up),h_node(i_dn),            & 
                          dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                          dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                          k_eff(j),                         &
                          A_term(j),                        & 
                          B_term(j))

                end do

                
                ! Path is between lowest stem and transporting root
                
                j     = n_hypool_ag
                i_up  = j+1
                i_dn  = j
                kmax_dn  = rootfr_scaler*cohort_hydr%kmax_stem_lower(n_hypool_stem)
                kmax_up  = rootfr_scaler*cohort_hydr%kmax_troot_upper

                call GetImTaylorKAB(kmax_up,kmax_dn,        &
                      ftc_node(i_up),ftc_node(i_dn),        & 
                      h_node(i_up),h_node(i_dn),            & 
                      dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                      dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                      k_eff(j),                         &
                      A_term(j),                        & 
                      B_term(j))

                ! Path is between the transporting root 
                ! and the absorbing root for this layer
                ! NOTE: No need to scale by root fraction
                ! even if in parallel mode, already parallel!
                
                j       = n_hypool_ag+1
                i_up    = j+1
                i_dn    = j
                kmax_dn = cohort_hydr%kmax_troot_lower(ilayer) 
                kmax_up = cohort_hydr%kmax_aroot_upper(ilayer)

                call GetImTaylorKAB(kmax_up,kmax_dn,        &
                      ftc_node(i_up),ftc_node(i_dn),        & 
                      h_node(i_up),h_node(i_dn),            & 
                      dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                      dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                      k_eff(j),                         &
                      A_term(j),                        & 
                      B_term(j))

                ! Path is between the absorbing root
                ! and the first rhizosphere shell nodes
                
                j     = n_hypool_ag+2
                i_up  = j+1
                i_dn  = j

                ! Special case. Maximum conductance depends on the 
                ! potential gradient.
                if(h_node(i_up) > h_node(i_dn) ) then
                    kmax_dn = 1._r8/(1._r8/cohort_hydr%kmax_aroot_lower(ilayer) + & 
                              1._r8/cohort_hydr%kmax_aroot_radial_in(ilayer))
                else
                    kmax_dn = 1._r8/(1._r8/cohort_hydr%kmax_aroot_lower(ilayer) + & 
                              1._r8/cohort_hydr%kmax_aroot_radial_out(ilayer))
                end if

                kmax_up = site_hydr%kmax_upper_shell(ilayer,1)*aroot_frac_plant

                call GetImTaylorKAB(kmax_up,kmax_dn,        &
                      ftc_node(i_up),ftc_node(i_dn),        & 
                      h_node(i_up),h_node(i_dn),            & 
                      dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                      dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                      k_eff(j),                         &
                      A_term(j),                        & 
                      B_term(j))

                ! Path is between rhizosphere shells

                do j = n_hypool_ag+3,n_hypool_tot-1

                    i_up = j+1
                    i_dn = j
                    ishell_up = i_up - (n_hypool_tot-nshell)
                    ishell_dn = i_dn - (n_hypool_tot-nshell)

                    kmax_dn = site_hydr%kmax_lower_shell(ilayer,ishell_dn)*aroot_frac_plant
                    kmax_up = site_hydr%kmax_upper_shell(ilayer,ishell_up)*aroot_frac_plant

                    call GetImTaylorKAB(kmax_up,kmax_dn,        &
                          ftc_node(i_up),ftc_node(i_dn),        & 
                          h_node(i_up),h_node(i_dn),            & 
                          dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), &
                          dpsi_dtheta_node(i_up), dpsi_dtheta_node(i_dn), &
                          k_eff(j),                         &
                          A_term(j),                        & 
                          B_term(j))

                end do

                ! -------------------------------------------------------------------------------
                ! Part 3.
                ! Loop through nodes again, build matrix
                ! -------------------------------------------------------------------------------
                
                tris_a(1) = 0._r8
                tris_b(1) = A_term(1) - denh2o*v_node(1)/dt_substep
                tris_c(1) = B_term(1)
                tris_r(1) = q_top_eff - k_eff(1)*(h_node(2)-h_node(1))


                do i = 2,n_hypool_tot-1
                    j = i
                    tris_a(i) = -A_term(j-1)
                    tris_b(i) = A_term(j) - B_term(j-1) - denh2o*v_node(i)/dt_substep
                    tris_c(i) = B_term(j)
                    tris_r(i) = -k_eff(j)*(h_node(i+1)-h_node(i)) + &
                          k_eff(j-1)*(h_node(i)-h_node(i-1))

                end do

                i = n_hypool_tot
                j = n_hypool_tot
                tris_a(i) = -A_term(j-1)
                tris_b(i) = -B_term(j-1) - denh2o*v_node(i)/dt_substep
                tris_c(i) = 0._r8
                tris_r(i) = k_eff(j-1)*(h_node(i)-h_node(i-1))


                ! Calculate the change in theta

                call Hydraulics_Tridiagonal(tris_a, tris_b, tris_c, tris_r, dth_node, tri_ierr)

                if(tri_ierr == 1) then
                    solution_found = .false.
                    error_code = 2
                    error_arr(:) = 0._r8
                    exit
                end if

                ! If we have not broken from the substep loop,
                ! that means this sub-step has been acceptable, and we may 
                ! go ahead and update the water content for the integrator
                
                th_node(:) = th_node(:) + dth_node(:)
                
                ! Mass error (flux - change)
                ! Total water mass in the plant at the beginning of this solve [kg h2o]
                w_tot_end = sum(th_node(:)*v_node(:))*denh2o
                
                wb_step_err = (q_top_eff*dt_substep) - (w_tot_beg-w_tot_end)
                
                if(abs(wb_step_err)>max_wb_step_err .or. any(dth_node(:).ne.dth_node(:)) )then
                    solution_found = .false.
                    error_code = 1
                    error_arr(:) = 0._r8
                    exit
                else
                    ! Note: this is somewhat of a default true. And the sub-steps
                    ! will keep going unless its changed and broken out of
                    ! the loop.
                    solution_found = .true.
                    error_code = 0
                end if

                ! Extra checks
                if( any(th_node(:)<0._r8) ) then
                    solution_found = .false.
                    error_code = 3
                    error_arr(:) = th_node(:)
                    exit
                end if

                ! Calculate new psi for checks
                do i = 1,n_hypool_plant
                    psi_node(i) = wrf_plant(pm_node(i),ft)%p%psi_from_th(th_node(i))
                end do
                do i = n_hypool_plant+1,n_hypool_tot
                    psi_node(i) = site_hydr%wrf_soil(ilayer)%p%psi_from_th(th_node(i))
                end do

                ! Accumulate the water balance error of the layer over the sub-steps
                ! for diagnostic purposes
                ! [kg/m2]
                wb_err_layer = wb_err_layer + wb_step_err

                ! -------------------------------------------------------------------------
                ! Diagnostics
                ! -------------------------------------------------------------------------
                
                ! Sapflow at the base of the tree is the flux rate
                ! between the transporting root node and the first stem node
                ! (note: a path j is between node i and i+1)
                ! [kg] = [kg/s] * [s]
                
                i = n_hypool_ag
                sapflow_lyr = sapflow_lyr + dt_substep * & 
                      (k_eff(i)*(h_node(i+1)-h_node(i)) + &  ! flux at (t) 
                      A_term(i)*dth_node(i)                 + &  ! dq at node i
                      B_term(i)*dth_node(i+1))                   ! dq at node i+1
                
                ! Root uptake is the integrated flux between the first rhizosphere
                ! shell and the absorbing root
                
                i = n_hypool_ag+2
                rootuptake_lyr = rootuptake_lyr + dt_substep * & 
                      (k_eff(i)*(h_node(i+1)-h_node(i)) + &  ! flux at (t) 
                      A_term(i)*dth_node(i)                 + &  ! dq at node i
                      B_term(i)*dth_node(i+1))                   ! dq at node i+1
                
                ! If debug mode is on, lets also track the mass fluxes across each
                ! path, and keep a running average of the effective conductances
                if(debug)then
                   do j=1,n_hypool_tot-1
                      k_diag(j) = k_diag(j) + k_eff(j)*dt_substep/dt_step
                      flux_diag(j) = flux_diag(j) + dt_substep * ( &
                           k_eff(j)*(h_node(j+1)-h_node(j)) + & 
                           A_term(j)*dth_node(j)+ B_term(j)*dth_node(j+1))
                   end do
                end if
                
            end do  ! do istep = 1,nsteps  (substep loop)

            iter=iter+1
            
        end do
    
        ! -----------------------------------------------------------
        ! Do a final check on water balance error sumed over sub-steps
        ! ------------------------------------------------------------
        if ( abs(wb_err_layer) > max_wb_err ) then
            
            write(fates_log(),*)'EDPlantHydraulics water balance error exceeds threshold of = ', max_wb_err
            write(fates_log(),*)'transpiration demand: ', dt_step*q_top_eff,' kg/step/plant'
            
            leaf_water = cohort_hydr%th_ag(1)*cohort_hydr%v_ag(1)*denh2o
            stem_water = sum(cohort_hydr%th_ag(2:n_hypool_ag) * &
                  cohort_hydr%v_ag(2:n_hypool_ag))*denh2o
            root_water = ( cohort_hydr%th_troot*cohort_hydr%v_troot + &
                  sum(cohort_hydr%th_aroot(:)*cohort_hydr%v_aroot_layer(:))) * denh2o
            
            write(fates_log(),*) 'leaf water: ',leaf_water,' kg/plant'
            write(fates_log(),*) 'stem_water: ',stem_water,' kg/plant'
            write(fates_log(),*) 'root_water: ',root_water,' kg/plant'
            write(fates_log(),*) 'LWP: ',cohort_hydr%psi_ag(1)
            write(fates_log(),*) 'dbh: ',cohort%dbh
            write(fates_log(),*) 'pft: ',cohort%pft
            write(fates_log(),*) 'tree lai: ',cohort%treelai,' m2/m2 crown'
            call endrun(msg=errMsg(sourcefile, __LINE__))
        end if


        ! If we have made it to this point, supposedly we have completed the whole time-step
        ! for this cohort x layer combination.  It is now safe to save the delta theta
        ! value and pass it back to the calling routine.  The value passed back is the
        ! change in theta over all sub-steps.
        
        dth_node(:) = th_node(:)-th_node_init(:)


        ! Add the current soil layer's contribution to total
        ! sap and root flux [kg]
        sapflow = sapflow + sapflow_lyr
        rootuptake(ilayer) = rootuptake_lyr
        
        
        ! Record the layer with the most iterations, but only
        ! if it greater than 1. It will default to zero
        ! if no layers took extra iterations.
        if( (real(iter)>cohort_hydr%iterh1) .and. (iter>1) )then
            cohort_hydr%iterlayer = real(ilayer)
        end if
        
        ! Save the number of times we refined our sub-step counts (iterh1)
        cohort_hydr%iterh1 = max(cohort_hydr%iterh1,real(iter))
        ! Save the number of sub-steps we ultimately used
        cohort_hydr%iterh2 = max(cohort_hydr%iterh2,real(nsteps))
        
        ! Update water contents in the relevant plant compartments [m3/m3]
        ! -------------------------------------------------------------------------------

        ! Leaf and above-ground stems
        cohort_hydr%th_ag(1:n_hypool_ag) = cohort_hydr%th_ag(1:n_hypool_ag) + dth_node(1:n_hypool_ag)
        ! Transporting root
        cohort_hydr%th_troot = cohort_hydr%th_troot + dth_node(n_hypool_ag+1)
        ! Absorbing root
        cohort_hydr%th_aroot(ilayer)  = cohort_hydr%th_aroot(ilayer) + dth_node(n_hypool_ag+2)

        ! Change in water per plant [kg/plant]
        dwat_plant = dwat_plant + & 
              (sum(dth_node(1:n_hypool_ag)*cohort_hydr%v_ag(1:n_hypool_ag)) + &
              dth_node(n_hypool_ag+1)*cohort_hydr%v_troot + & 
              dth_node(n_hypool_ag+2)*cohort_hydr%v_aroot_layer(ilayer))*denh2o
        
        ! Remember the error for the cohort
        wb_err_plant = wb_err_plant + wb_err_layer

        ! Save the change in water mass in the rhizosphere. Note that we did
        ! not immediately update the state variables upon completing each
        ! plant-layer solve. We accumulate the difference, and apply them
        ! after all cohort-layers are complete. This allows each cohort
        ! to experience the same water conditions (for good or bad).
        
        if(site_hydr%l_aroot_layer(ilayer)<nearzero)then
            write(fates_log(),*) 'no site level root?'
            write(fates_log(),*) ilayer,site_hydr%l_aroot_layer(ilayer)
            write(fates_log(),*) cohort_hydr%l_aroot_layer(ilayer)
            call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        
        dth_layershell_col(ilayer,:) = dth_layershell_col(ilayer,:) + &
              dth_node((n_hypool_tot-nshell+1):n_hypool_tot) * & 
              cohort_hydr%l_aroot_layer(ilayer) * &
              cohort%n / site_hydr%l_aroot_layer(ilayer)
        
     enddo !soil layer (jj -> ilayer)

    end associate
    return
  end subroutine ImTaylorSolve1D

  ! =====================================================================================

  subroutine Report1DError(cohort, site_hydr, ilayer, z_node, v_node, & 
                           th_node, q_top_eff, dt_step, w_tot_beg, w_tot_end, &
                           rootfr_scaler, aroot_frac_plant, err_code, err_arr)

    ! This routine reports what the initial condition to the 1D solve looks
    ! like, and then quits.

                                                                    ! Arguments (IN)
    type(ed_cohort_type),intent(in),target      :: cohort
    type(ed_site_hydr_type),intent(in), target  :: site_hydr
    integer, intent(in)                         :: ilayer           ! soil layer index of interest
    real(r8), intent(in)                        :: z_node(:)        ! elevation of nodes
    real(r8), intent(in)                        :: v_node(:)        ! volume of nodes
    real(r8), intent(in)                        :: th_node(:)       ! water content of node
    real(r8), intent(in)                        :: dt_step          ! time [seconds] over-which to calculate solution
    real(r8), intent(in)                        :: q_top_eff        ! transpiration flux rate at upper boundary [kg -s]
    real(r8), intent(in)                        :: w_tot_beg        ! total water mass at beginning of step [kg]
    real(r8), intent(in)                        :: w_tot_end        ! total water mass at end of step [kg]
    real(r8), intent(in)                        :: rootfr_scaler    ! What is the root fraction in this layer?
    real(r8), intent(in)                        :: aroot_frac_plant ! What fraction of total absorbring roots
    ! in the soil continuum is from current plant?
    integer, intent(in)                         :: err_code         ! error code
    real(r8), intent(in)                        :: err_arr(:)       ! error diagnostic

    type(ed_cohort_hydr_type),pointer  :: cohort_hydr
    integer :: i
    integer :: ft
    real(r8)              :: leaf_water
    real(r8)              :: stem_water
    real(r8)              :: troot_water
    real(r8)              :: aroot_water
    real(r8), allocatable :: psi_node(:)
    real(r8), allocatable :: h_node(:)

    cohort_hydr => cohort%co_hydr
    ft = cohort%pft

    allocate(psi_node(size(z_node)))
    allocate(h_node(size(z_node)))

    write(fates_log(),*) 'Could not find a stable solution for hydro 1D solve'
    write(fates_log(),*) ''
    write(fates_log(),*) 'error code: ',err_code
    write(fates_log(),*) 'error diag: ',err_arr(:)

    do i = 1,n_hypool_plant
       psi_node(i) =  wrf_plant(site_hydr%pm_node(i),ft)%p%psi_from_th(th_node(i)) 
       h_node(i) =  mpa_per_pa*denh2o*grav_earth*z_node(i) + psi_node(i)
    end do
    do i = n_hypool_plant+1,n_hypool_tot
       psi_node(i) =  site_hydr%wrf_soil(ilayer)%p%psi_from_th(th_node(i)) 
       h_node(i) =  mpa_per_pa*denh2o*grav_earth*z_node(i) + psi_node(i)
    end do


    leaf_water = sum(cohort_hydr%th_ag(1:n_hypool_leaf)* &
         cohort_hydr%v_ag(1:n_hypool_leaf))*denh2o
    stem_water = sum(cohort_hydr%th_ag(n_hypool_leaf+1:n_hypool_ag) * &
         cohort_hydr%v_ag(n_hypool_leaf+1:n_hypool_ag))*denh2o
    troot_water = (cohort_hydr%th_troot*cohort_hydr%v_troot) * denh2o
    aroot_water = sum(cohort_hydr%th_aroot(:)*cohort_hydr%v_aroot_layer(:)) * denh2o
    
    write(fates_log(),*) 'layer: ',ilayer
    write(fates_log(),*) 'wb_step_err = ',(q_top_eff*dt_step) - (w_tot_beg-w_tot_end)
    write(fates_log(),*) 'leaf water: ',leaf_water,' kg/plant'
    write(fates_log(),*) 'stem_water: ',stem_water,' kg/plant'
    write(fates_log(),*) 'troot_water: ',troot_water
    write(fates_log(),*) 'aroot_water: ',aroot_water
    write(fates_log(),*) 'LWP: ',cohort_hydr%psi_ag(1)
    write(fates_log(),*) 'dbh: ',cohort%dbh
    write(fates_log(),*) 'pft: ',cohort%pft
    write(fates_log(),*) 'z nodes: ',z_node(:)
    write(fates_log(),*) 'psi_z: ',h_node(:)-psi_node(:)
    write(fates_log(),*) 'vol,    theta,   H,      kmax-'
    write(fates_log(),*) 'flux:          ', q_top_eff*dt_step 
    write(fates_log(),*) 'l:',v_node(1),th_node(1),h_node(1),psi_node(1)
    write(fates_log(),*) '                      ',cohort_hydr%kmax_stem_upper(1)*rootfr_scaler
    write(fates_log(),*) 's:',v_node(2),th_node(2),h_node(2),psi_node(2)
    write(fates_log(),*) '                      ',1._r8/(1._r8/(cohort_hydr%kmax_stem_lower(1)*rootfr_scaler) + 1._r8/(cohort_hydr%kmax_troot_upper*rootfr_scaler))
    write(fates_log(),*) 't:',v_node(3),th_node(3),h_node(3)
    write(fates_log(),*) '                      ',1._r8/(1._r8/cohort_hydr%kmax_troot_lower(ilayer)+ 1._r8/cohort_hydr%kmax_aroot_upper(ilayer))
    write(fates_log(),*) 'a:',v_node(4),th_node(4),h_node(4)
    write(fates_log(),*) '                   in:',1._r8/(1._r8/cohort_hydr%kmax_aroot_radial_in(ilayer) + &
                                                         1._r8/(site_hydr%kmax_upper_shell(ilayer,1)*aroot_frac_plant)     + & 
                                                         1._r8/cohort_hydr%kmax_aroot_upper(ilayer))
    write(fates_log(),*) '                  out:',1._r8/(1._r8/cohort_hydr%kmax_aroot_radial_out(ilayer) + & 
                                                         1._r8/(site_hydr%kmax_upper_shell(ilayer,1)*aroot_frac_plant)     + & 
                                                         1._r8/cohort_hydr%kmax_aroot_upper(ilayer))
    write(fates_log(),*) 'r1:',v_node(5),th_node(5),h_node(5)
    write(fates_log(),*) '                      ',1._r8/(1._r8/(site_hydr%kmax_lower_shell(ilayer,1)*aroot_frac_plant) + 1._r8/(site_hydr%kmax_upper_shell(ilayer,2)*aroot_frac_plant))
    write(fates_log(),*) 'r2:',v_node(6),th_node(6),h_node(6)
    write(fates_log(),*) '                      ',1._r8/(1._r8/(site_hydr%kmax_lower_shell(ilayer,2)*aroot_frac_plant) + 1._r8/(site_hydr%kmax_upper_shell(ilayer,3)*aroot_frac_plant))
    write(fates_log(),*) 'r3:',v_node(7),th_node(7),h_node(7)
    write(fates_log(),*) '                      ',1._r8/(1._r8/(site_hydr%kmax_lower_shell(ilayer,3)*aroot_frac_plant) + 1._r8/(site_hydr%kmax_upper_shell(ilayer,4)*aroot_frac_plant))
    write(fates_log(),*) 'r4:',v_node(8),th_node(8),h_node(8)
    write(fates_log(),*) '                      ',1._r8/(1._r8/(site_hydr%kmax_lower_shell(ilayer,4)*aroot_frac_plant) + 1._r8/(site_hydr%kmax_upper_shell(ilayer,5)*aroot_frac_plant))
    write(fates_log(),*) 'r5:',v_node(9),th_node(9),h_node(9)
    write(fates_log(),*) 'kmax_aroot_radial_out: ',cohort_hydr%kmax_aroot_radial_out(ilayer)
    write(fates_log(),*) 'surf area of root: ',2._r8 * pi_const * EDPftvarcon_inst%hydr_rs2(ft) * cohort_hydr%l_aroot_layer(ilayer)
    write(fates_log(),*) 'aroot_frac_plant: ',aroot_frac_plant,cohort_hydr%l_aroot_layer(ilayer),site_hydr%l_aroot_layer(ilayer)
    write(fates_log(),*) 'kmax_upper_shell: ',site_hydr%kmax_lower_shell(ilayer,:)*aroot_frac_plant
    write(fates_log(),*) 'kmax_lower_shell: ',site_hydr%kmax_upper_shell(ilayer,:)*aroot_frac_plant
    write(fates_log(),*) ''
    write(fates_log(),*) 'tree lai: ',cohort%treelai,' m2/m2 crown'
    write(fates_log(),*) 'area and area to volume ratios'
    write(fates_log(),*) ''
    write(fates_log(),*) 'a:',v_node(4)
    write(fates_log(),*) '                      ',2._r8 * pi_const * EDPftvarcon_inst%hydr_rs2(ft) * cohort_hydr%l_aroot_layer(ilayer)
    write(fates_log(),*) 'r1:',v_node(5)
    write(fates_log(),*) '                      ',2._r8 * pi_const * site_hydr%r_out_shell(ilayer,1) * cohort_hydr%l_aroot_layer(ilayer)
    write(fates_log(),*) 'r2:',v_node(6)
    write(fates_log(),*) '                      '
    write(fates_log(),*) 'r3:',v_node(7)
    write(fates_log(),*) '                      '
    write(fates_log(),*) 'r4:',v_node(8)
    write(fates_log(),*) '                      '
    write(fates_log(),*) 'r5:',v_node(9)
    
    write(fates_log(),*) 'inner shell kmaxs: ',site_hydr%kmax_lower_shell(:,1)*aroot_frac_plant
    




    deallocate(psi_node)
    deallocate(h_node)


    ! Most likely you will want to end-run after this routine, but maybe not...

    return
  end subroutine Report1DError

  ! =================================================================================

  subroutine GetImTaylorKAB(kmax_up,kmax_dn, &
       ftc_up,ftc_dn, &
       h_up,h_dn, &
       dftc_dtheta_up, dftc_dtheta_dn, &
       dpsi_dtheta_up, dpsi_dtheta_dn, &
       k_eff,   &
       a_term,  & 
       b_term)

    ! -----------------------------------------------------------------------------
    ! This routine will return the effective conductance "K", as well
    ! as two terms needed to calculate the implicit solution (using taylor
    ! first order expansion).  The two terms are generically named A & B.
    ! Thus the name "KAB".  These quantities are specific not to the nodes
    ! themselves, but to the path between the nodes, defined as positive
    ! direction towards atmosphere, from "up"stream side (closer to soil)
    ! and the "d"ow"n" stream side (closer to air)
    ! -----------------------------------------------------------------------------
    ! Arguments
    real(r8),intent(in)    :: kmax_dn, kmax_up               ! max conductance [kg s-1 Mpa-1]
    real(r8),intent(inout) :: ftc_dn, ftc_up                 ! frac total conductance [-]
    real(r8),intent(in)    :: h_dn, h_up                     ! total potential [Mpa]
    real(r8),intent(inout) :: dftc_dtheta_dn, dftc_dtheta_up ! Derivative
                                                             ! of FTC wrt relative water content
    real(r8),intent(in)    :: dpsi_dtheta_dn, dpsi_dtheta_up ! Derivative of matric potential
                                                             ! wrt relative water content
    real(r8),intent(out)   :: k_eff                          ! effective conductance over path [kg s-1 Mpa-1]
    real(r8),intent(out)   :: a_term                         ! "A" term for path (See tech note)
    real(r8),intent(out)   :: b_term                         ! "B" term for path (See tech note)

                                                             ! Locals
    real(r8)               :: h_diff                         ! Total potential difference [MPa]


    ! Calculate difference in total potential over the path [MPa]
    h_diff  = h_up - h_dn

    ! If we do enable "upstream K", then we are saying that 
    ! the fractional loss of conductivity is dictated
    ! by the upstream side of the flow.  In this case,
    ! the change in ftc is only non-zero on that side, and is
    ! zero'd otherwise.

    if(do_upstream_k) then

       if (h_diff>0._r8) then
          ftc_dn       = ftc_up
          dftc_dtheta_dn = 0._r8
       else
          ftc_up         = ftc_dn
          dftc_dtheta_up = 0._r8
       end if

    end if

    ! Calculate total effective conductance over path  [kg s-1 MPa-1]
    k_eff = 1._r8/(1._r8/(ftc_up*kmax_up)+1._r8/(ftc_dn*kmax_dn))

    ! "A" term, which operates on the downstream node (closer to atm)
    a_term = k_eff**2.0_r8 * h_diff * kmax_dn**(-1.0_r8) * ftc_dn**(-2.0_r8) &
         * dftc_dtheta_dn - k_eff * dpsi_dtheta_dn

       
    ! "B" term, which operates on the upstream node (further from atm)
    b_term = k_eff**2.0_r8 * h_diff * kmax_up**(-1.0_r8) * ftc_up**(-2.0_r8) & 
         * dftc_dtheta_up + k_eff * dpsi_dtheta_up
       
    

    return
  end subroutine GetImTaylorKAB

  ! =====================================================================================

  subroutine GetKAndDKDPsi(kmax_dn,kmax_up, &
       h_dn,h_up, &
       ftc_dn,ftc_up, &
       dftc_dpsi_dn, & 
       dftc_dpsi_up, &
       dk_dpsi_dn, &
       dk_dpsi_up, & 
       k_eff)

    ! -----------------------------------------------------------------------------
    ! This routine will return the effective conductance "K", as well
    ! as two terms needed to calculate the implicit solution (using taylor
    ! first order expansion).  The two terms are generically named A & B.
    ! Thus the name "KAB".  These quantities are specific not to the nodes
    ! themselves, but to the path between the nodes, defined as positive
    ! direction from "up"per (closer to atm) and "lo"wer (further from atm).
    ! -----------------------------------------------------------------------------

    real(r8),intent(in)    :: kmax_dn        ! max conductance (downstream) [kg s-1 Mpa-1]
    real(r8),intent(in)    :: kmax_up        ! max conductance (upstream)   [kg s-1 Mpa-1]
    real(r8),intent(in)    :: h_dn           ! total potential (downstream) [MPa]
    real(r8),intent(in)    :: h_up           ! total potential (upstream)   [Mpa]
    real(r8),intent(inout) :: ftc_dn         ! frac total cond (downstream) [-]
    real(r8),intent(inout) :: ftc_up         ! frac total cond (upstream)   [-]
    real(r8),intent(inout) :: dftc_dpsi_dn   ! derivative ftc / theta (downstream)
    real(r8),intent(inout) :: dftc_dpsi_up   ! derivative ftc / theta (upstream)
    
                                             ! of FTC wrt relative water content
    real(r8),intent(out)   :: dk_dpsi_dn     ! change in effective conductance from the
                                             ! downstream pressure node
    real(r8),intent(out)   :: dk_dpsi_up     ! change in effective conductance from the
                                             ! upstream pressure node
    real(r8),intent(out)   :: k_eff          ! effective conductance over path [kg s-1 Mpa-1]

    ! Locals
    real(r8)               :: h_diff                         ! Total potential difference [MPa]
                 ! the effective fraction of total 
                                                             ! conductivity is either governed
                                                             ! by the upstream node, or by both
                                                             ! with a harmonic average
    ! Calculate difference in total potential over the path [MPa]
    h_diff  = h_up - h_dn

    ! If we do enable "upstream K", then we are saying that 
    ! the fractional loss of conductivity is dictated
    ! by the upstream side of the flow.  In this case,
    ! the change in ftc is only non-zero on that side, and is
    ! zero'd otherwise.

    if(do_upstream_k) then

       if (h_diff>0._r8) then
          ftc_dn         = ftc_up
          dftc_dpsi_dn = 0._r8
       else
          ftc_up         = ftc_dn
          dftc_dpsi_up = 0._r8
       end if

    end if

    ! Calculate total effective conductance over path  [kg s-1 MPa-1]
    k_eff = 1._r8/(1._r8/(ftc_up*kmax_up)+1._r8/(ftc_dn*kmax_dn))

    
    dk_dpsi_dn = k_eff**2._r8  * kmax_dn**(-1._r8) * ftc_dn**(-2._r8) * dftc_dpsi_dn

    dk_dpsi_up = k_eff**2._r8  * kmax_up**(-1._r8) * ftc_up**(-2._r8) * dftc_dpsi_up

    

    return
  end subroutine GetKAndDKDPsi
  

  subroutine AccumulateMortalityWaterStorage(csite,ccohort,delta_n)

    ! ---------------------------------------------------------------------------
    ! This subroutine accounts for the water bound in plants that have
    ! just died. This water is accumulated at the site level for all plants
    ! that die.
    ! In another routine, this pool is reduced as water vapor flux, and
    ! passed to the HLM.
    ! ---------------------------------------------------------------------------

    ! Arguments

    type(ed_site_type), intent(inout), target     :: csite
    type(ed_cohort_type) , intent(inout), target  :: ccohort
    real(r8), intent(in)                          :: delta_n ! Loss in number density
    ! for this cohort /ha/day

    real(r8) :: delta_w                                      !water change due to mortality Kg/m2
    ! Locals
    type(ed_site_hydr_type), pointer              :: csite_hydr
    type(ed_cohort_hydr_type), pointer            :: ccohort_hydr

    ccohort_hydr => ccohort%co_hydr
    csite_hydr   => csite%si_hydr
    delta_w =   (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:))      + &
         ccohort_hydr%th_troot*ccohort_hydr%v_troot                  + &
         sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
         denh2o*delta_n*AREA_INV

    csite_hydr%h2oveg_dead = csite_hydr%h2oveg_dead + delta_w


    csite_hydr%h2oveg = csite_hydr%h2oveg - delta_w

    return
  end subroutine AccumulateMortalityWaterStorage

  !-------------------------------------------------------------------------------!

  subroutine RecruitWaterStorage(nsites,sites,bc_out)

    ! ---------------------------------------------------------------------------
    ! This subroutine accounts for the water bound in plants that have
    ! just recruited. This water is accumulated at the site level for all plants
    ! that recruit.
    ! Because this water is taken from the soil in hydraulics_bc, which will not 
    ! be called until the next timestep, this water is subtracted out of
    ! plant_stored_h2o_si to ensure HLM water balance at the beg_curr_day timestep.
    ! plant_stored_h2o_si will include this water when calculated in hydraulics_bc
    ! at the next timestep, when it gets pulled from the soil water.
    ! ---------------------------------------------------------------------------

    ! Arguments
    integer, intent(in)                       :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_out_type), intent(inout)          :: bc_out(nsites)

    ! Locals
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_type), pointer :: currentPatch
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: csite_hydr
    integer :: s

    if( hlm_use_planthydro.eq.ifalse ) return

    do s = 1,nsites

       csite_hydr => sites(s)%si_hydr
       csite_hydr%h2oveg_recruit = 0.0_r8
       currentPatch => sites(s)%oldest_patch 
       do while(associated(currentPatch))
          currentCohort=>currentPatch%tallest
          do while(associated(currentCohort))
             ccohort_hydr => currentCohort%co_hydr
             if(ccohort_hydr%is_newly_recruited) then
                csite_hydr%h2oveg_recruit = csite_hydr%h2oveg_recruit + &
                     (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:)) + &
                     ccohort_hydr%th_troot*ccohort_hydr%v_troot       + &
                     sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                     denh2o*currentCohort%n
             end if
             currentCohort => currentCohort%shorter
          enddo !cohort
          currentPatch => currentPatch%younger
       enddo !end patch loop

       csite_hydr%h2oveg_recruit      = csite_hydr%h2oveg_recruit * AREA_INV

    end do

    return
  end subroutine RecruitWaterStorage

  ! =====================================================================================

  ! =====================================================================================
  ! Utility Functions
  ! =====================================================================================

  subroutine bisect_rootfr(a, b, lower_init, upper_init, xtol, ytol, crootfr, x_new)
    !
    ! !DESCRIPTION: Bisection routine for getting the inverse of the cumulative root
    !  distribution. No analytical soln bc crootfr ~ exp(ax) + exp(bx).
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)     :: a, b        ! pft root distribution constants
    real(r8)      , intent(in)     :: lower_init  ! lower bound of initial x estimate [m]
    real(r8)      , intent(in)     :: upper_init  ! upper bound of initial x estimate [m]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x_new         [m]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for crootfr       [-]
    real(r8)      , intent(in)     :: crootfr     ! cumulative root fraction at x_new [-]
    real(r8)      , intent(out)    :: x_new       ! soil depth                        [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                  ! lower bound x estimate [m]
    real(r8) :: upper                  ! upper bound x estimate [m]
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    !----------------------------------------------------------------------

    lower = lower_init
    upper = upper_init
    f_lo  = zeng2001_crootfr(a, b, lower) - crootfr
    f_hi  = zeng2001_crootfr(a, b, upper) - crootfr
    chg   = upper - lower
    do while(abs(chg) .gt. xtol)
       x_new = 0.5_r8*(lower + upper)
       f_new = zeng2001_crootfr(a, b, x_new) - crootfr
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
    end do
  end subroutine bisect_rootfr

  ! =====================================================================================

  function zeng2001_crootfr(a, b, z, z_max) result(crootfr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: a,b    ! pft parameters
    real(r8) , intent(in) :: z      ! soil depth (m)
    real(r8) , intent(in), optional :: z_max ! max soil depth (m)
    !
    real(r8) :: crootfr_max

    ! !RESULT
    real(r8) :: crootfr            ! cumulative root fraction
    !
    !------------------------------------------------------------------------
    crootfr      = 1._r8 - .5_r8*(exp(-a*z) + exp(-b*z))


    ! If a maximum rooting depth is provided, then
    ! we force everything to sum to unity. We do this by
    ! simply dividing through by the maximum possible
    ! root fraction.

    if(present(z_max))then
       crootfr_max = 1._r8 - .5_r8*(exp(-a*z_max) + exp(-b*z_max))
       crootfr = crootfr/crootfr_max
    end if

    if(debug)then
       if(present(z_max))then
          if((crootfr_max<nearzero) .or. (crootfr_max>1.0_r8) )then
             write(fates_log(),*) 'problem scaling crootfr in zeng2001'
             write(fates_log(),*) 'z_max: ',z_max
             write(fates_log(),*) 'crootfr_max: ',crootfr_max
          end if
       end if
    end if


    return

  end function zeng2001_crootfr

  ! =====================================================================================

  subroutine shellGeom(l_aroot, rs1, area_site, dz, r_out_shell, r_node_shell, v_shell)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil surrounding fine roots remains
    ! the same.
    !
    ! !USES:

    !
    ! !ARGUMENTS:
    real(r8)     , intent(in)             :: l_aroot              ! Total length of absorbing roots
    ! for the whole site, this layer (m)
    real(r8)     , intent(in)             :: rs1                  ! Fine root radius (m)
    real(r8)     , intent(in)             :: area_site            ! Area of site (10,000 m2)
    real(r8)     , intent(in)             :: dz                   ! Width of current soil layer (m)
    real(r8)     , intent(out)            :: r_out_shell(:)       ! Outer radius of each shell (m)
    real(r8)     , intent(out)            :: r_node_shell(:)      ! Radius of the shell's midpoint
    real(r8)     , intent(out)            :: v_shell(:)           ! volume of the rhizosphere shells (m3/ha)
    ! for this layer
    !
    ! !LOCAL VARIABLES:
    integer                        :: k            ! rhizosphere shell indicies
    integer                        :: nshells      ! We don't use the global because of unit testing
    !-----------------------------------------------------------------------

       
    nshells = size(r_out_shell,dim=1)
    
    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    r_out_shell(nshells) = (pi_const*l_aroot/(area_site*dz))**(-0.5_r8)                  ! eqn(8) S98
    if(nshells > 1) then
       do k = 1,nshells-1
          r_out_shell(k)   = rs1*(r_out_shell(nshells)/rs1)**((real(k,r8))/real(nshells,r8))  ! eqn(7) S98
       enddo
    end if
       
    ! set nodal (midpoint) radii of these shells
    ! BOC...not doing this as it requires PFT-specific fine root thickness, but this is at column level
    r_node_shell(1) = 0.5_r8*(rs1 + r_out_shell(1))
    !r_node_shell(1) = 0.5_r8*(r_out_shell(1))
    
    do k = 2,nshells
       r_node_shell(k) = 0.5_r8*(r_out_shell(k-1) + r_out_shell(k))
    enddo
    
    ! update volumes
    do k = 1,nshells
       if(k == 1) then
          v_shell(k) = pi_const*l_aroot*(r_out_shell(k)**2._r8 - rs1**2._r8)
       else
          v_shell(k) = pi_const*l_aroot*(r_out_shell(k)**2._r8 - r_out_shell(k-1)**2._r8)
       end if
    enddo

    return
  end subroutine shellGeom

  ! =====================================================================================

  function xylemtaper(p, dz) result(chi_tapnotap)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: p      ! Savage et al. (2010) taper exponent                                                                [-]
    real(r8) , intent(in) :: dz     ! hydraulic distance from petiole to node of interest                                                [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: atap,btap           ! scaling exponents for total conductance ~ tree size (ratio of stem radius to terminal twig radius)
    real(r8) :: anotap,bnotap       ! same as atap, btap, but not acounting for xylem taper (Savage et al. (2010) p = 0)
    ! NOTE: these scaling exponents were digitized from Fig 2a of Savage et al. (2010)
    ! Savage VM, Bentley LP, Enquist BJ, Sperry JS, Smith DD, Reich PB, von Allmen EI. 2010.
    !    Hydraulic trade-offs and space filling enable better predictions of vascular structure
    !    and function in plants. Proceedings of the National Academy of Sciences 107(52): 22722-22727.
    real(r8) :: lN=0.04_r8          ! petiole length                                                                                     [m]
    real(r8) :: little_n=2._r8      ! number of daughter branches per parent branch, assumed constant throughout tree (self-similarity)  [-]
    real(r8) :: big_n               ! number of branching levels (allowed here to take on non-integer values): increases with tree size  [-]
    real(r8) :: ktap                ! hydraulic conductance along the pathway, accounting for xylem taper                                [kg s-1 MPa-1]
    real(r8) :: knotap              ! hydraulic conductance along the pathway, not accounting for xylem taper                            [kg s-1 MPa-1]
    real(r8) :: num                 ! temporary
    real(r8) :: den                 ! temporary
    !
    ! !RESULT
    real(r8) :: chi_tapnotap        ! ratio of total tree conductance accounting for xylem taper to that without, over interval dz
    !
    !------------------------------------------------------------------------

    anotap  = 7.19903e-13_r8
    bnotap  = 1.326105578_r8
    if (p >= 1.0_r8) then
       btap  = 2.00586217_r8
       atap  = 1.82513E-12_r8
    else if (p >= (1._r8/3._r8) .AND. p < 1._r8) then
       btap  = 1.854812819_r8
       atap  = 6.66908E-13_r8
    else if (p >= (1._r8/6._r8) .AND. p < (1._r8/3._r8)) then
       btap  = 1.628179741_r8
       atap  = 6.58345E-13_r8
    else
       btap  = bnotap
       atap  = anotap
    end if

    num          = 3._r8*log(1._r8 - dz/lN * (1._r8-little_n**(1._r8/3._r8)))
    den          = log(little_n)
    big_n        = num/den - 1._r8
    ktap         = atap   * (little_n**(big_N*  btap/2._r8))
    knotap       = anotap * (little_n**(big_N*bnotap/2._r8))
    chi_tapnotap = ktap / knotap

    return

  end function xylemtaper
  
  ! =====================================================================================
  
  subroutine Hydraulics_Tridiagonal(a, b, c, r, u, ierr)
    !
    ! !DESCRIPTION: An abbreviated version of biogeophys/TridiagonalMod.F90
    !
    ! This solves the form:
    !
    ! a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = r(i)
    !
    ! It assumed that coefficient a(1) and c(N) DNE as there is
    ! no u(0) or u(N-1).
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8), intent(in)    :: a(:)           ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(:)           ! "b" diagonal column of tridiagonal matrix
    real(r8), intent(in)    :: c(:)           ! "c" right off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: r(:)           ! "r" forcing term of tridiagonal matrix
    real(r8), intent(out)   :: u(:)           ! solution
    integer,  intent(out)   :: ierr           ! flag: 0=passed,  1=failed
    !
    ! !LOCAL VARIABLES:
    real(r8) :: bet                           ! temporary
    real(r8) :: gam(10)                       ! temporary
    integer  :: k                             ! index
    integer  :: N                             ! Size of the matrix
    real(r8) :: err                           ! solution error, in units of [m3/m3]
    real(r8) :: rel_err                       ! relative error, normalized by delta theta
    real(r8), parameter :: allowable_rel_err = 0.0001_r8

    !----------------------------------------------------------------------
    N=size(r,dim=1)
    bet = b(1)
    do k=1,N
       if(k == 1) then
          u(k)   = r(k) / bet
       else
          gam(k) = c(k-1) / bet
          bet    = b(k) - a(k) * gam(k)
          u(k)   = (r(k) - a(k)*u(k-1)) / bet
       end if
    enddo

    do k=N-1,1,-1
       u(k)   = u(k) - gam(k+1) * u(k+1)
    enddo

    ! If debug mode, calculate error on the forward solution
    ierr = 0
    if(debug)then
       do k=1,N
          if(k==1)then
             err = abs(r(k) - (b(k)*u(k)+c(k)*u(k+1)))
          elseif(k<N) then
             err = abs(r(k) - (a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)))
          else
             err = abs(r(k) - (a(k)*u(k-1)+b(k)*u(k)))
          end if
          if(abs(u(k))>nearzero)then
              rel_err = abs(err/u(k))
              if( ((rel_err > allowable_rel_err) .and. (err > max_wb_step_err)) .or. & 
                   (err /= err) )then
                  write(fates_log(),*) 'Tri-diagonal solve produced solution with'
                  write(fates_log(),*) 'non-negligable error.'
                  write(fates_log(),*) 'Compartment: ',k
                  write(fates_log(),*) 'Error in forward solution: ',err
                  write(fates_log(),*) 'Estimated delta theta: ',u(k)
                  write(fates_log(),*) 'Rel Error: ',rel_err
                  write(fates_log(),*) 'Reducing time-step'
                  ierr = 1
              end if
          end if
       end do
    end if

  end subroutine Hydraulics_Tridiagonal

  ! =====================================================================================

  subroutine MatSolve2D(site_hydr,cohort,cohort_hydr, &
                        tmx,qtop, &
                        sapflow,rootuptake,wb_err_plant , dwat_plant, & 
                        dth_layershell_site)


    ! ---------------------------------------------------------------------------------
    ! This solution to the plant water flux equations casts all the fluxes through a
    ! cohort, and the rhizosphere shells in ALL layers as a single system of equations.
    ! If thinking of the plant's above ground components as one dimension, and the soil
    ! layers as another, this is a somewhat 2D system (hence "Matrix" in the name).
    ! To improve the quality of the solution and reduce solver error, this also
    ! uses a Newton iteration.  See technical documentation for a full derivation
    ! of the mathematics.  However, in brief, we can describe the flux balance through
    ! any node, considering flux paths labeled j, through that node in set J.
    ! This is an implicit solve, so we balance the change in water mass (defined by
    ! volume V, density rho, and water content theta) with the flux (q) esitmated
    ! at the next time-step q^(t+1).  Note that we continue to solve this equation, using
    ! updated values of water content and pressure (psi), by balancing our fluxes with
    ! the total of previous (theta_p) and remaining (theta_r) water contents.
    !
    !  rho V                 rho V
    !  -----  Del theta_p +  ----- Del theta_r  =  Sum ( q^(t+1) )
    !  Del t                 Del t                  J
    !
    ! The flux at t+1, is simply the current flux (q) and a first order Taylor
    ! expanion (i.e. forward-euler) estimate with the current derivative based
    ! on the current value of theta and psi.
    ! Note also, that the solution is in terms of the matric potential, psi.  This
    ! conversion from theta to psi, requires this derivative (Jacobian) to also
    ! contain not just the rate of change of flux wrt psi, but the change in theta
    ! wrt psi (self term, no cross node terms).
    !
    ! -----------------------------------------------------------------------------------

    
    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    type(ed_site_hydr_type), intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(ed_cohort_hydr_type), target            :: cohort_hydr
    type(ed_cohort_type) , intent(inout), target :: cohort
    real(r8),intent(in)                          :: tmx ! time interval to integrate over [s]
    real(r8),intent(in)                          :: qtop
    real(r8),intent(out) :: sapflow                   ! time integrated mass flux between transp-root and stem [kg]
    real(r8),intent(out) :: rootuptake(:)             ! time integrated mass flux between rhizosphere and aroot [kg]
   
    
    real(r8),intent(out)                         :: wb_err_plant ! total error over plant, transpiration 
    ! should match change in storage [kg/m2]
    real(r8),intent(out)                         :: dwat_plant   ! total change in water mass for the plant [kg]
    real(r8),intent(inout)                       :: dth_layershell_site(:,:)

    integer :: nsteps       ! Number of rounds of attempts we have made
    integer :: i                 ! generic index (sometimes node index)
    integer :: inode             ! node index
    integer :: k                 ! generic node index
    integer :: j, icnx           ! soil layer and connection indices
    integer :: id_dn, id_up      ! Node indices on each side of flux path
    integer :: ishell            ! rhizosphere shell index
     
    integer :: icnv              ! Convergence flag for each solve, see flag definitions
                                 ! below.
                    
    real(r8) :: aroot_frac_plant ! Fraction of rhizosphere this plant "owns"

    real(r8) :: dqflx_dpsi_dn    ! Derivative, change in mass flux per change
                                 ! in matric potential of the down-stream node
                                 ! [kg s-1 Mpa-1]

    real(r8) :: dqflx_dpsi_up    ! Derivative, change in mass flux per change
                                 ! in matric potential of the up-stream node
    ! [kg s-1 Mpa-1]

    real(r8) :: dk_dpsi_dn     ! change in effective conductance from the
                               ! downstream pressure node
    real(r8) :: dk_dpsi_up     ! change in effective conductance from the
                                             ! upstream pressure node

    real(r8) :: residual_amax    ! maximum absolute mass balance residual over all
                                 ! nodes,
                                 ! used for determining convergence. At the point
    
    real(r8) :: rsdx             ! Temporary residual while determining max value


    real(r8) :: rlfx_soil        ! Pressure update reduction factor for soil compartments
    real(r8) :: rlfx_plnt        ! Pressure update reduction factor for plant comparmtents

    real(r8) :: tm               ! Total time integrated after each substep [s]
    real(r8) :: dtime              ! Total time to be integrated this step [s]
    real(r8) :: w_tot_beg     ! total plant water prior to solve [kg]
    real(r8) :: w_tot_end     ! total plant water at end of solve [kg]
    
    real(r8) :: k_eff ! Effective conductivity over the current pathway
                      ! between two nodes.  Factors in fractional
                      ! loss of conductivity on each side of the pathway, and the material maximum
                      ! conductivity on each side  [kg/s/MPa]
    integer :: icnx_ar        ! Connection index of the aroot <-> rhizosphere shell
    
    integer :: nsd               ! node index of highest residual
    integer :: nwtn_iter         ! number of (Newton) iterations on each substep

                                 ! to get a succesfull Newton solve.
    integer :: kshell            ! rhizosphere shell index, 1->nshell
    
    integer :: info
    integer :: nstep             !number of time steps


    ! This is a convergence test.  This is the maximum difference
    ! allowed between the flux balance and the change in storage
    ! on a node. [kg/s] *Note, 1.e-9 = 1 ug/s
    real(r8), parameter :: max_allowed_residual = 1.e-8_r8

    ! Maximum number of times we re-try a round of Newton
    ! iterations, each time decreasing the time-step and
    ! potentially reducing relaxation factors
    integer, parameter :: max_newton_rounds = 10
    
    ! Maximum number of Newton iterations in each round
    integer, parameter :: max_newton_iter = 200

    ! Flag definitions for convergence flag (icnv)
    ! icnv = 1 fail the round due to either wacky math, or
    !          too many Newton iterations
    ! icnv = 2 continue onto next iteration, 
    ! icnv = 3 acceptable solution
    ! icnv = 4 too many failures, aborting
    
    integer, parameter :: icnv_fail_round    = 1
    integer, parameter :: incv_cont_search   = 2
    integer, parameter :: icnv_pass_round    = 3
    integer, parameter :: icnv_complete_fail = 4

    ! Timestep reduction factor when a round of
    ! newton iterations fail. 
    
    real(r8), parameter :: dtime_rf = 0.2_r8

    real(r8), parameter :: rlfx_soil0 = 0.1  ! Initial Pressure update
                                             ! reduction factor for soil compartments
    real(r8), parameter :: rlfx_plnt0 = 0.6  ! Initial Pressure update
                                             ! reduction factor for plant comparmtents


    associate(conn_up      => site_hydr%conn_up, &
              conn_dn      => site_hydr%conn_dn, &
              kmax_up      => site_hydr%kmax_up, &
              kmax_dn      => site_hydr%kmax_dn, &
              q_flux       => site_hydr%q_flux, & 
              residual     => site_hydr%residual, &
              ajac         => site_hydr%ajac, &
              ipiv         => site_hydr%ipiv, & 
              th_node      => site_hydr%th_node, &
              th_node_init => site_hydr%th_node_init, &
              psi_node     => site_hydr%psi_node, &
              pm_node      => site_hydr%pm_node, & 
              ftc_node     => site_hydr%ftc_node, & 
              z_node       => site_hydr%z_node, & 
              v_node       => site_hydr%v_node, &
              dth_node     => site_hydr%dth_node, &
              node_layer   => site_hydr%node_layer, &
              h_node       => site_hydr%h_node, &
              dftc_dpsi_node => site_hydr%dftc_dpsi_node, &
              ft           => cohort%pft)


      !for debug only
      nstep = get_nstep()
      

      ! This NaN's the scratch arrays
      call site_hydr%FlushSiteScratch()

      ! This is the maximum number of iterations needed for this cohort
      ! (each soil layer has a different number, this saves the max)
      cohort_hydr%iterh1 = 0
      cohort_hydr%iterh2 = 0

      ! These are output fluxes from the subroutine, total integrated
      ! mass fluxes [kg] over the time-step. sapflow is the integrated
      ! flux between the transporting root and the 1st stem compartment.
      ! The rootuptake is the integrated flux between the 1st rhizosphere
      ! and absorbing roots
      sapflow = 0._r8
      rootuptake(:) = 0._r8

      ! Chnage in water content, over all substeps [m3/m3]
      dth_node(:) = 0._r8
      
      ! Transfer node heights, volumes and initial water contents for
      ! the transporting root and above ground compartments to the
      ! complete node vector

      do i = 1,n_hypool_ag+n_hypool_troot
         if (i<=n_hypool_ag) then
            z_node(i)  = cohort_hydr%z_node_ag(i)
            v_node(i)  = cohort_hydr%v_ag(i)
            th_node_init(i) = cohort_hydr%th_ag(i)
         elseif (i>n_hypool_ag) then
            z_node(i)  = cohort_hydr%z_node_troot
            v_node(i)  = cohort_hydr%v_troot
            th_node_init(i) = cohort_hydr%th_troot
         end if
      end do
      
      ! Transfer node-heights, volumes and intiial water contents
      ! for below-ground components,
      ! from the cohort structures, into the complete node vector
      i = n_hypool_ag + n_hypool_troot
      
      do j = 1,site_hydr%nlevrhiz

         ! Calculate the fraction of the soil layer
         ! folume that this plant's rhizosphere accounts forPath is across the upper an lower rhizosphere comparment
         ! on each side of the nodes. Since there is no flow across the outer
         ! node to the edge, we ignore that last half compartment
         aroot_frac_plant = cohort_hydr%l_aroot_layer(j)/site_hydr%l_aroot_layer(j)
         
         do k = 1, n_hypool_aroot + nshell
            i = i + 1
            if (k==1) then
               z_node(i)  = -site_hydr%zi_rhiz(j)
               v_node(i)  = cohort_hydr%v_aroot_layer(j)
               th_node_init(i) = cohort_hydr%th_aroot(j)
            else
               kshell  = k-1
               z_node(i)  = -site_hydr%zi_rhiz(j)
               ! The volume of the Rhizosphere for a single plant
               v_node(i)  = site_hydr%v_shell(j,kshell)*aroot_frac_plant
               th_node_init(i) = site_hydr%h2osoi_liqvol_shell(j,kshell)
            end if
         enddo

      enddo

     
      
      ! Total water mass in the plant at the beginning of this solve [kg h2o]
      w_tot_beg = sum(th_node_init(:)*v_node(:))*denh2o

      
      ! Initialize variables and flags that track
      ! the progress of the solve

      tm      = 0
      nsteps  = 0

      outerloop: do while(tm < tmx)

         ! If we are here, then we either are starting the solve,
         ! or, we just completed a solve but did not fully integrate
         ! the time.  Lets update the time-step to be the remainder
         ! of the step.
         dtime = min(tmx*0.01,tmx-tm)
          
         ! Relaxation factors are reset to starting point.
         rlfx_plnt = rlfx_plnt0
         rlfx_soil = rlfx_soil0

         ! Return here if we want to start a new round of Newton
         ! iterations. The previous round was unsucessful either
         ! because it couldn't get a zero residual, or because
         ! a singularity was encountered.
100      continue

         ! Set the current water content as the initial [m3/m3]
         th_node(:) = th_node_init(:)
         
          
         tm = tm + dtime
         nwtn_iter = 0

         ! Return here if you are just continuing the
         ! Newton search for a solution. No need to
         ! update timing information.
200      continue

         nwtn_iter = nwtn_iter + 1

         ! The Jacobian and the residual are incremented,
         ! and the Jacobian is sparse, thus they both need
         ! to be zerod.
         ajac(:,:)   = 0._r8
         residual(:) = 0._r8
        
         
         do k=1,site_hydr%num_nodes

            ! This is the storage gained from previous newton iterations.
            residual(k) = residual(k) + (th_node(k) - th_node_init(k))/dtime*denh2o*v_node(k)
            
            if(pm_node(k) == rhiz_p_media) then

               j = node_layer(k)
               psi_node(k) = site_hydr%wrf_soil(j)%p%psi_from_th(th_node(k))
!!               if ( abs(th_node(k)-site_hydr%wrf_soil(j)%p%th_from_psi(psi_node(k))) > nearzero) then
!!                  print*,'non-reversible WRTs?'
!!                  print*,psi_node(k)
!!                  print*,th_node(k)
!!                  print*,site_hydr%wrf_soil(j)%p%th_from_psi(psi_node(k))
!!                  stop
!!               end if


               
               ! Get total potential [Mpa]
               h_node(k) =  mpa_per_pa*denh2o*grav_earth*z_node(k) + psi_node(k)
               ! Get Fraction of Total Conductivity [-]
               ftc_node(k) = site_hydr%wkf_soil(j)%p%ftc_from_psi(psi_node(k))
               ! deriv ftc wrt psi
               dftc_dpsi_node(k)   = site_hydr%wkf_soil(j)%p%dftcdpsi_from_psi(psi_node(k))

            else

               psi_node(k) = wrf_plant(pm_node(k),ft)%p%psi_from_th(th_node(k))
               ! Get total potential [Mpa]
               h_node(k) =  mpa_per_pa*denh2o*grav_earth*z_node(k) + psi_node(k)
               ! Get Fraction of Total Conductivity [-]
               ftc_node(k) = wkf_plant(pm_node(k),ft)%p%ftc_from_psi(psi_node(k))
               ! deriv ftc wrt psi
               dftc_dpsi_node(k)   = wkf_plant(pm_node(k),ft)%p%dftcdpsi_from_psi(psi_node(k))

            end if

            ! Fill the self-term on the Jacobian's diagonal with the
            ! the change in storage wrt change in psi.
            
            if(pm_node(k) == rhiz_p_media) then
                ajac(k,k) = denh2o*v_node(k)/ &
                      (site_hydr%wrf_soil(j)%p%dpsidth_from_th(th_node(k))*dtime)
            else
                ajac(k,k) = denh2o*v_node(k)/ &
                      (wrf_plant(pm_node(k),ft)%p%dpsidth_from_th(th_node(k))*dtime)
            endif
            
         enddo

!         do i=1,site_hydr%num_nodes
!            print*,i,node_layer(i),pm_node(i),z_node(i),v_node(i),th_node_init(i),psi_node(i),h_node(i)
!         end do
!         stop
         

         ! Calculations of maximum conductance for upstream and downstream sides
         ! of each connection.  This IS dependant on total potential h_node
         ! because of the root-soil radial conductance.

         call SetMaxCondConnections(site_hydr, cohort_hydr, h_node, kmax_dn, kmax_up)
         
         ! calculate boundary fluxes     
         do icnx=1,site_hydr%num_connections

            id_dn = conn_dn(icnx)
            id_up = conn_up(icnx)

            ! The row (first index) of the Jacobian (ajac) represents the
            ! the node for which we are calculating the water balance
            ! The column (second index) of the Jacobian represents the nodes
            ! on which the pressure differentials effect the water balance
            ! of the node of the first index.

            ! This will get the effective K, and may modify FTC depending
            ! on the flow direction

            call GetKAndDKDPsi(kmax_dn(icnx), &
                               kmax_up(icnx), &
                               h_node(id_dn), &
                               h_node(id_up), &
                               ftc_node(id_dn), &
                               ftc_node(id_up), &
                               dftc_dpsi_node(id_dn), & 
                               dftc_dpsi_node(id_up), & 
                               dk_dpsi_dn, &
                               dk_dpsi_up, & 
                               k_eff)

            q_flux(icnx) = k_eff*(h_node(id_up)-h_node(id_dn))
            
            ! See equation (22) in technical documentation
            ! Add fluxes at current time to the residual
            residual(id_dn) = residual(id_dn) - q_flux(icnx)
            residual(id_up) = residual(id_up) + q_flux(icnx)

            ! This is the Jacobian term related to the pressure changes on the down-stream side
            ! and these are applied to both the up and downstream sides (oppositely)
            dqflx_dpsi_dn = -k_eff + (h_node(id_up)-h_node(id_dn)) * dk_dpsi_dn

            ! This is the Jacobian term related to the pressure changes on the up-stream side
            ! and these are applied to both the up and downstream sides (oppositely)
            dqflx_dpsi_up =  k_eff + (h_node(id_up)-h_node(id_dn)) * dk_dpsi_up

            ! Down-stream node's contribution to the down-stream node's Jacobian
            ajac(id_dn,id_dn) = ajac(id_dn,id_dn) + dqflx_dpsi_dn

            ! Down-stream node's contribution to the up-stream node's Jacobian
            ajac(id_up,id_dn) = ajac(id_up,id_dn) - dqflx_dpsi_dn

            ! Up-stream node's contribution to the down-stream node's Jacobian
            ajac(id_dn,id_up) = ajac(id_dn,id_up) + dqflx_dpsi_up

            ! Up-stream node's contribution to the up-stream node's Jacobian
            ajac(id_up,id_up) = ajac(id_up,id_up) - dqflx_dpsi_up

            
         enddo

         ! Add the transpiration flux (known, retrieved from photosynthesis scheme)
         ! to the mass balance on the leaf (1st) node.  This is constant over
         ! the time-step, so no Jacobian term needed (yet)
         
         residual(1) = residual(1) + qtop
         

         ! Start off assuming things will pass, then find numerous
         ! ways to see if it failed
         icnv = icnv_pass_round

         
         ! check residual
         ! if(nstep==15764) print *,'ft,it,residual_amax-',ft,nwtn_iter,residual_amax,'qtop',qtop,psi_node,
         ! 'init-',psi_node_init,'resi-',residual, 'qflux-',q_flux,'v_n',v_node

         ! If we have performed any Newton iterations, then the residual
         ! may reflect a flux that balances (equals) the change in storage. If this is
         ! true, then the residual is zero, and we are done with the sub-step. If it is
         ! not nearly zero, then we must continue our search and perform another solve.
         
         residual_amax = 0._r8
         nsd = 0
         do k = 1, site_hydr%num_nodes
            rsdx = abs(residual(k))
            ! check NaNs
            if( rsdx /= rsdx ) then
               icnv = icnv_fail_round
               exit
            endif
            if( rsdx > residual_amax ) then
               residual_amax = rsdx
               nsd = k
            endif
         enddo

         if(icnv == icnv_fail_round) goto 199
         
         ! If the solution is balanced, none of the residuals
         ! should be very large, and we can ignore another
         ! solve attempt.
         if( residual_amax < max_allowed_residual ) then

             goto 201

         ! In this case, we still have a non-trivially small
         ! residual, yet we have exceeded our iteration cap
         ! Thus we set error flag to 1, which forces a time-step
         ! shortening
         elseif( nwtn_iter > max_newton_iter) then

             icnv = icnv_fail_round
             goto 199


         ! We still have some residual (perhaps this is first step),
         ! have not used too many steps, so we go ahead
         ! and perform a Newton iteration
         else

            ! We wont actually know if we have a good solution
            ! until we complete this step and re-calculate the residual
            ! so we simply flag that we continue the search
            icnv = incv_cont_search

            ! ---------------------------------------------------------------------------
            ! From Lapack documentation
            !
            ! subroutine dgesv(integer 	N  (in),
            !                  integer 	NRHS (in),
            !                  real(r8), dimension( lda, * ) 	A (in/out),
            !                  integer 	LDA (in),
            !                  integer, dimension( * ) 	IPIV (out),
            !                  real(r8), dimension( ldb, * ) 	B (in/out),
            !                  integer 	LDB (in),
            !                  integer 	INFO (out) )
            !
            ! DGESV computes the solution to a real system of linear equations
            ! A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
            ! The LU decomposition with partial pivoting and row interchanges is
            ! used to factor A as   A = P * L * U,
            ! where P is a permutation matrix, L is unit lower triangular, and U is
            ! upper triangular.  The factored form of A is then used to solve the
            ! system of equations A * X = B.
            !
            ! N is the number of linear equations, i.e., the order of the
            ! matrix A.  N >= 0.
            !
            ! NRHS is the number of right hand sides, i.e., the number of columns
            ! of the matrix B.  NRHS >= 0. 
            !
            ! A: 
            ! On entry, the N-by-N coefficient matrix A.
            ! On exit, the factors L and U from the factorization
            ! A = P*L*U; the unit diagonal elements of L are not stored.
            !
            ! LDA is the leading dimension of the array A.  LDA >= max(1,N).
            !
            ! IPIV is the pivot indices that define the permutation matrix P;
            ! row i of the matrix was interchanged with row IPIV(i).
            !
            ! B
            ! On entry, the N-by-NRHS matrix of right hand side matrix B.
            ! On exit, if INFO = 0, the N-by-NRHS solution matrix X.
            !
            ! LDB is the leading dimension of the array B.  LDB >= max(1,N).
            !
            ! INFO:
            ! = 0:  successful exit
            ! < 0:  if INFO = -i, the i-th argument had an illegal value
            ! > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
            !    has been completed, but the factor U is exactly
            !    singular, so the solution could not be computed.
            ! ---------------------------------------------------------------------------
            

            call DGESV(site_hydr%num_nodes,1,ajac,site_hydr%num_nodes,ipiv,residual,site_hydr%num_nodes,info)

            
            if ( info < 0 ) then
               write(fates_log(),*) 'illegal value generated in DGESV() linear'
               write(fates_log(),*) 'system solver, see node: ',-info
               call endrun(msg=errMsg(sourcefile, __LINE__))
            END IF
            if ( info > 0 ) then
               write(fates_log(),*) 'the factorization of linear system in DGESV() generated'
               write(fates_log(),*) 'a singularity at node: ',info
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
            
            ! If info == 0, then
            ! lapack was able to generate a solution.
            ! For A * X = B, 
            ! Where the residual() was B, DGESV() returns
            ! the solution X into the residual array.
            
            ! Update the matric potential of each node.  Since this is a search
            ! we update matric potential as only a fraction of delta psi (residual)
            
            do k = 1, site_hydr%num_nodes
               
               if(pm_node(k) == rhiz_p_media) then
                  psi_node(k) = psi_node(k) + residual(k) * rlfx_soil
                  j = node_layer(k)
             !     print*,'psi:',psi_node(k),residual(k),k,j
                  th_node(k)  = site_hydr%wrf_soil(j)%p%th_from_psi(psi_node(k))
               else
           !       print*,'psi:',psi_node(k),residual(k),k
                  psi_node(k) = psi_node(k) + residual(k) * rlfx_plnt
                  th_node(k) = wrf_plant(pm_node(k),ft)%p%th_from_psi(psi_node(k))
               endif
               
            enddo

!            stop
            
         endif
         
199      continue
         
         if( icnv == icnv_fail_round ) then

            write(fates_log(),'(10x,a)') '---  Convergence Failure  ---'
            write(fates_log(),'(4x,a,1pe11.4,2(a,i6),1pe11.4)') 'Equation Maximum Residual = ', &
                 residual_amax,' Node = ',nsd, 'pft = ',ft, 'qtop: ',qtop

            ! If we have not exceeded our max number
            ! of retrying rounds of Newton iterations, reduce
            ! time and try a new round
            if( nsteps < max_newton_rounds ) then

               tm = tm - dtime
               nsteps = nsteps + 1

               write(fates_log(),*) 'fates hydraulics, MatSolve2D:'
               write(fates_log(),'(4x,a,1pe11.4,1x,2a,1pe11.4,1x,a)')       &
                    'Time Step Reduced From ',dtime,'s',' to ', &
                    min(dtime * dtime_rf,tmx-tm),'s'
               
               dtime = min(dtime * dtime_rf,tmx-tm)
               
               do k = 1,site_hydr%num_nodes
                  th_node(k) = th_node_init(k)
               enddo

               ! Decrease the relaxation factors 
               rlfx_plnt = rlfx_plnt0*(0.9_r8**real(nsteps,r8))
               rlfx_soil = rlfx_soil0*(0.9_r8**real(nsteps,r8))
               
               !
               !---  Number of time step reductions failure: stop simulation  ---
               !
            else
               ! Complete failure to converge even with re-trying
               ! iterations with smaller timestepps and relaxations
               icnv = icnv_complete_fail
            endif
            
         endif

      

         
         if(icnv == icnv_fail_round) then
            goto 100
         elseif(icnv == incv_cont_search) then

            ! THIS MAY BE A GOOD PLACE TO INCREASE
            ! THE RELAXATION FACTORS
            goto 200

         elseif(icnv == icnv_pass_round) then
            dth_node(:) = dth_node(:) + (th_node(:) - th_node_init(:))
            goto 201
         elseif(icnv == icnv_complete_fail) then
             write(fates_log(),*) 'Newton hydraulics solve'
             write(fates_log(),*) 'could not converge on a solution.'
             write(fates_log(),*) 'Perhaps try increasing iteration cap,'
             write(fates_log(),*) 'and decreasing relaxation factors.'
             call endrun(msg=errMsg(sourcefile, __LINE__))
         else
             write(fates_log(),*) 'unhandled failure mode in'
             write(fates_log(),*) 'newton hydraulics solve'
             write(fates_log(),*) 'icnv = ',icnv
             call endrun(msg=errMsg(sourcefile, __LINE__))
         endif


         ! If we have reached this point, we have iterated to
         ! a stable solution (where residual mass balance = 0)
         ! It is possible that we have used a sub-step though,
         ! and need to continue the iteration.

201      continue

         ! Save the number of substeps needed
         cohort_hydr%iterh1 = cohort_hydr%iterh1 + 1

         ! Save the  max number of Newton iterations needed
         cohort_hydr%iterh2 = max(cohort_hydr%iterh2,real(nwtn_iter))

         print*,'Completed a newton solve'
         print*,psi_node(:)
         stop

         ! Save flux diagnostics
         ! ------------------------------------------------------

         sapflow = sapflow + q_flux(n_hypool_ag)*dtime

         do j = 1,site_hydr%nlevrhiz
            ! Connection betwen the 1st rhizosphere and absorbing roots
            icnx_ar = n_hypool_ag + (j-1)*(nshell+1)+2
            rootuptake(j) = q_flux(icnx_ar)*dtime
         enddo


         ! If there are any sub-steps left, we need to update
         ! the initial water content
         th_node_init(:) = th_node(:)
         
      end do outerloop

     

      ! If we have made it here, we have successfully integrated
      ! the water content. Transfer this from scratch space
      ! into the cohort memory structures for plant compartments,
      ! and increment the site-level change in soil moistures

      
      ! Update state variables in plant compartments
      cohort_hydr%th_ag(1:n_hypool_ag) = cohort_hydr%th_ag(1:n_hypool_ag) + dth_node(1:n_hypool_ag)
      cohort_hydr%th_troot             = cohort_hydr%th_troot + dth_node(n_hypool_ag+1)

      ! Change in water per plant [kg/plant]
      dwat_plant = sum(dth_node(1:n_hypool_ag+n_hypool_troot)*v_node(1:n_hypool_ag+n_hypool_troot))*denh2o

      inode = n_hypool_ag+n_hypool_troot
      do j = 1,site_hydr%nlevrhiz
         do k = 1, nshell+1
            inode = inode + 1
            if(k==1) then
               cohort_hydr%th_aroot(j) = cohort_hydr%th_aroot(j)+dth_node(inode)
               dwat_plant = dwat_plant + (dth_node(inode) * v_node(inode))*denh2o
            else
               ishell = k-1
               dth_layershell_site(j,ishell) = dth_layershell_site(j,ishell) + &
                    dth_node(inode) * cohort_hydr%l_aroot_layer(j) * &
                    cohort%n / site_hydr%l_aroot_layer(j)
               
            endif
         enddo
      enddo
      
      ! Total water mass in the plant at the end of this solve [kg h2o]
      w_tot_end = sum(th_node(:)*v_node(:))*denh2o
      
      ! Mass error (flux - change) [kg/m2]
      wb_err_plant = (qtop*dtime)-(w_tot_beg-w_tot_end)

            
     end associate

    return   
  end subroutine MatSolve2D

  ! =====================================================================================
  
  function SumBetweenDepths(site_hydr,depth_t,depth_b,array_in) result(depth_sum)

    ! This function sums the quantity in array_in between depth_t (top)
    ! and depth_b.  It assumes many things. Firstly, that the depth coordinates
    ! for array_in do match site_hydr%zi_rhiz (on rhizosphere layers), and that
    ! those coordinates are positive down.

    type(ed_site_hydr_type), intent(in) :: site_hydr
    real(r8),intent(in)    :: depth_t      ! Top Depth    (positive coordinate)
    real(r8),intent(in)    :: depth_b      ! Bottom depth (positive coordinate)
    real(r8),intent(in)    :: array_in(:)  ! Quantity to be summed (flux?mass?)
    real(r8)               :: depth_sum    ! The summed result we return in units (/depth)
    integer  :: i_rhiz_t                   ! Layer index of top full layer
    integer  :: i_rhiz_b                   ! layer index of bottom full layer
    integer  :: nlevrhiz                   ! Number of rhizosphere layers (not shells)
    real(r8) :: frac                       ! Fraction of partial layer, by depth
    
    i_rhiz_t = count((site_hydr%zi_rhiz-site_hydr%dz_rhiz)<depth_t)+1  ! First layer completely below top depth
    i_rhiz_b = count(site_hydr%zi_rhiz<depth_b)                        ! Last layer completely above bottom depth
    nlevrhiz = size(site_hydr%zi_rhiz)

    depth_sum = 0._r8
    
    ! Trivial, the top depth is deeper than the soil column
    ! return... 0 or nan..?
    if(i_rhiz_t>nlevrhiz) then
       return
    end if
    
    ! Sum all fully encased layers
    if(i_rhiz_b>=i_rhiz_t)then
       depth_sum = depth_sum + sum(array_in(i_rhiz_t:i_rhiz_b)) 
    end if
    
    ! Find fraction contribution from top partial layer (if any)
    if(i_rhiz_t>1) then
       frac = (site_hydr%zi_rhiz(i_rhiz_t-1)-depth_t)/site_hydr%dz_rhiz(i_rhiz_t-1)
       depth_sum = depth_sum + frac*array_in(i_rhiz_t-1)
    end if
    
    ! Find fraction contribution from bottom partial layer (if any)
    if(i_rhiz_b<nlevrhiz) then
       frac = (depth_b-site_hydr%zi_rhiz(i_rhiz_b))/site_hydr%dz_rhiz(i_rhiz_b+1)
       depth_sum = depth_sum + frac*array_in(i_rhiz_b+1)
    end if
    
    depth_sum = depth_sum/(min(depth_b,site_hydr%zi_rhiz(nlevrhiz))-depth_t)
    
  end function SumBetweenDepths
  
  ! =====================================================================================
  
  subroutine SetMaxCondConnections(site_hydr, cohort_hydr, h_node, kmax_dn, kmax_up)
    
    ! -------------------------------------------------------------------------------
    ! This subroutine sets the maximum conductances
    ! on the downstream (towards atm) and upstream (towards
    ! soil) side of each connection. This scheme is somewhat complicated
    ! by the fact that the direction of flow at the root surface impacts
    ! which root surface radial conductance to use, which makes these calculation
    ! dependent on the updating potential in the system, and not just a function
    ! of plant geometry and material properties.
    ! -------------------------------------------------------------------------------
    
    type(ed_site_hydr_type), intent(in),target   :: site_hydr
    type(ed_cohort_hydr_type), intent(in),target :: cohort_hydr
    real(r8),intent(in)  :: h_node(:)        ! Total (matric+height) potential at each node (Mpa)
    real(r8),intent(out) :: kmax_dn(:)       ! Max conductance of downstream sides of connections (kg s-1 MPa-1)
    real(r8),intent(out) :: kmax_up(:)       ! Max conductance of upstream sides of connections   (kg s-1 MPa-1)

    real(r8):: aroot_frac_plant ! Fraction of the cohort's fine-roots
                                ! out of the total in the current layer
    integer :: icnx  ! connection index
    integer :: inode ! node index
    integer :: istem ! stem index
    integer :: k     ! rhizosphere/root index (per level)
    integer :: j     ! soil layer index
    
    kmax_dn(:) = fates_unset_r8
    kmax_up(:) = fates_unset_r8

    ! Set leaf to stem connections (only 1 leaf layer
    ! this will break if we have multiple, as there would
    ! need to be assumptions about which compartment
    ! to connect the leaves to.
    icnx  = 1
    kmax_dn(icnx) = cohort_hydr%kmax_petiole_to_leaf
    kmax_up(icnx) = cohort_hydr%kmax_stem_upper(1)

    ! Stem to stem connections
    do istem = 1,n_hypool_stem-1
       icnx = icnx + 1
       kmax_dn(icnx) = cohort_hydr%kmax_stem_lower(istem)
       kmax_up(icnx) = cohort_hydr%kmax_stem_upper(istem+1)
    enddo

    ! Path is between lowest stem and transporting root
    icnx  = icnx + 1
    kmax_dn(icnx) = cohort_hydr%kmax_stem_lower(n_hypool_stem)
    kmax_up(icnx) = cohort_hydr%kmax_troot_upper

    ! Path is between the transporting root and the absorbing roots
    inode = n_hypool_ag
    do j = 1,site_hydr%nlevrhiz

       aroot_frac_plant = cohort_hydr%l_aroot_layer(j)/site_hydr%l_aroot_layer(j)
       
       do k = 1, n_hypool_aroot + nshell
          icnx = icnx + 1
          inode = inode + 1
          if( k == 1 ) then !troot-aroot
             kmax_dn(icnx) = cohort_hydr%kmax_troot_lower(j) 
             kmax_up(icnx) = cohort_hydr%kmax_aroot_upper(j)

          elseif( k == 2) then ! aroot-soil

             ! Special case. Maximum conductance depends on the 
             ! potential gradient.

             if(h_node(inode) < h_node(inode+1) ) then
                kmax_dn(icnx) = 1._r8/(1._r8/cohort_hydr%kmax_aroot_lower(j) + & 
                     1._r8/cohort_hydr%kmax_aroot_radial_in(j))
             else
                kmax_dn(icnx) = 1._r8/(1._r8/cohort_hydr%kmax_aroot_lower(j) + & 
                     1._r8/cohort_hydr%kmax_aroot_radial_out(j))
             end if
             kmax_up(icnx) = site_hydr%kmax_upper_shell(j,1)*aroot_frac_plant

          else                 ! soil - soil
             kmax_dn(icnx) = site_hydr%kmax_lower_shell(j,k-2)*aroot_frac_plant
             kmax_up(icnx) = site_hydr%kmax_upper_shell(j,k-1)*aroot_frac_plant
          endif
       enddo

    end do


  end subroutine SetMaxCondConnections
  
  ! =====================================================================================

  subroutine InitHydroGlobals()

    ! This routine allocates the Water Transfer Functions (WTFs)
    ! which include both water retention functions (WRFs)
    ! as well as the water conductance (K) functions (WKFs)
    ! But, this is only for plants! These functions have specific
    ! parameters, potentially, for each plant functional type and
    ! each organ (pft x organ), but this can be used globally (across
    ! all sites on the node (machine) to save memory.  These functions
    ! are also applied to soils, but since soil properties vary with
    ! soil layer and location, those functions are bound to the site
    ! structure, and are therefore not "global".

    ! Define
    class(wrf_type_vg), pointer :: wrf_vg
    class(wkf_type_vg), pointer :: wkf_vg
    class(wrf_type_cch), pointer :: wrf_cch
    class(wkf_type_tfs), pointer :: wkf_tfs
    class(wrf_type_tfs), pointer :: wrf_tfs

    integer :: ft            ! PFT index
    integer :: pm            ! plant media index
    integer :: inode         ! compartment node index
    real(r8) :: cap_corr     ! correction for nonzero psi0x (TFS)
    real(r8) :: cap_slp      ! slope of capillary region of curve
    real(r8) :: cap_int      ! intercept of capillary region of curve

    if(hlm_use_planthydro.eq.ifalse) return
    
    ! we allocate from stomata_p_media, which should be zero

    allocate(wrf_plant(stomata_p_media:n_plant_media,numpft))
    allocate(wkf_plant(stomata_p_media:n_plant_media,numpft))
  
    ! -----------------------------------------------------------------------------------
    ! Initialize the Water Retention Functions
    ! -----------------------------------------------------------------------------------

    select case(plant_wrf_type)
    case(van_genuchten_type)
       do ft = 1,numpft
            do pm = 1, n_plant_media
                allocate(wrf_vg)
                wrf_plant(pm,ft)%p => wrf_vg
                call wrf_vg%set_wrf_param([alpha_vg, psd_vg, th_sat_vg, th_res_vg])
            end do
        end do
     case(campbell_type)
        do ft = 1,numpft
           do pm = 1,n_plant_media
              allocate(wrf_cch)
              wrf_plant(pm,ft)%p => wrf_cch
              call wrf_cch%set_wrf_param([EDPftvarcon_inst%hydr_thetas_node(ft,pm), &
                                          EDPftvarcon_inst%hydr_pinot_node(ft,pm), &
                                          9._r8])
           end do
        end do
     case(tfs_type)
        do ft = 1,numpft
           do pm = 1,n_plant_media
              allocate(wrf_tfs)
              wrf_plant(pm,ft)%p => wrf_tfs

              if (pm.eq.leaf_p_media) then   ! Leaf tissue
                 cap_slp    = 0.0_r8
                 cap_int    = 0.0_r8
                 cap_corr   = 1.0_r8
              else               ! Non leaf tissues
                 cap_slp    = (hydr_psi0 - hydr_psicap )/(1.0_r8 - rwccap(pm))  
                 cap_int    = -cap_slp + hydr_psi0    
                 cap_corr   = -cap_int/cap_slp
              end if
              
              call wrf_tfs%set_wrf_param([EDPftvarcon_inst%hydr_thetas_node(ft,pm), &
                                          EDPftvarcon_inst%hydr_resid_node(ft,pm), &
                                          EDPftvarcon_inst%hydr_pinot_node(ft,pm), &
                                          EDPftvarcon_inst%hydr_epsil_node(ft,pm), &
                                          rwcft(pm), & 
                                          cap_corr, &
                                          cap_int, &
                                          cap_slp,real(pm,r8)])
           end do
        end do

    end select

    ! -----------------------------------------------------------------------------------
    ! Initialize the Water Conductance (K) Functions
    ! -----------------------------------------------------------------------------------

    select case(plant_wkf_type)
    case(van_genuchten_type)
        do ft = 1,numpft
            do pm = 1, n_plant_media
                allocate(wkf_vg)
                wkf_plant(pm,ft)%p => wkf_vg
                call wkf_vg%set_wkf_param([alpha_vg, psd_vg, th_sat_vg, th_res_vg, tort_vg])
            end do
           
        end do
    case(campbell_type)
        write(fates_log(),*) 'campbell/clapp-hornberger conductance not used in plants'
        call endrun(msg=errMsg(sourcefile, __LINE__))
    case(tfs_type)
        do ft = 1,numpft
            do pm = 1, n_plant_media
               allocate(wkf_tfs)
               wkf_plant(pm,ft)%p => wkf_tfs
               call wkf_tfs%set_wkf_param([EDPftvarcon_inst%hydr_p50_node(ft,pm), &
                                       EDPftvarcon_inst%hydr_avuln_node(ft,pm)])
            end do
        end do
    end select

    ! There is only 1 stomata conductance hypothesis which uses the p50 and 
    ! vulnerability parameters
    ! -----------------------------------------------------------------------------------

    do ft = 1,numpft
       allocate(wkf_tfs)
       wkf_plant(stomata_p_media,ft)%p => wkf_tfs
       call wkf_tfs%set_wkf_param([EDPftvarcon_inst%hydr_p50_gs(ft), &
                              EDPftvarcon_inst%hydr_avuln_gs(ft)])
    end do

    
    return
  end subroutine InitHydroGlobals

  !! subroutine UpdateLWPMemFLCMin(ccohort_hydr)

  ! This code may be re-introduced at a later date (rgk 08-2019)

  ! SET COHORT-LEVEL BTRAN FOR USE IN NEXT TIMESTEP
  ! first update the leaf water potential memory
  !!    do t=2, numLWPmem
  !!ccohort_hydr%lwp_mem(t-1)      = ccohort_hydr%lwp_mem(t)
  !!end do
  !!ccohort_hydr%lwp_mem(numLWPmem)   = ccohort_hydr%psi_ag(1)
  !!call flc_gs_from_psi(cCohort, ccohort_hydr%psi_ag(1))

  !!refill_rate = -log(0.5)/(ccohort_hydr%refill_days*24._r8*3600._r8)   ! s-1
  !!do k=1,n_hypool_ag
  !!ccohort_hydr%flc_min_ag(k) = min(ccohort_hydr%flc_min_ag(k), ccohort_hydr%flc_ag(k))
  !!if(ccohort_hydr%psi_ag(k) >= ccohort_hydr%refill_thresh .and. &
  !!      ccohort_hydr%flc_ag(k) > ccohort_hydr%flc_min_ag(k)) then   ! then refilling
  !!    ccohort_hydr%flc_min_ag(k) = ccohort_hydr%flc_ag(k) - &
  !!          (ccohort_hydr%flc_ag(k) - ccohort_hydr%flc_min_ag(k))*exp(-refill_rate*dtime)
  !!end if
  !!end do
  !!do k=1,n_hypool_troot
  !!ccohort_hydr%flc_min_troot(k) = min(ccohort_hydr%flc_min_troot(k), ccohort_hydr%flc_troot(k))
  !!if(ccohort_hydr%psi_troot(k) >= ccohort_hydr%refill_thresh .and. &
  !!      ccohort_hydr%flc_troot(k) > ccohort_hydr%flc_min_troot(k)) then   ! then refilling
  !!    ccohort_hydr%flc_min_troot(k) = ccohort_hydr%flc_troot(k) - &
  !!          (ccohort_hydr%flc_troot(k) - ccohort_hydr%flc_min_troot(k))*exp(-refill_rate*dtime)
  !!end if
  !!end do
  !!do j=1,site_hydr%nlevrhiz
  !!ccohort_hydr%flc_min_aroot(j) = min(ccohort_hydr%flc_min_aroot(j), ccohort_hydr%flc_aroot(j))
  !!if(ccohort_hydr%psi_aroot(j) >= ccohort_hydr%refill_thresh .and. &
  !!      ccohort_hydr%flc_aroot(j) > ccohort_hydr%flc_min_aroot(j)) then   ! then refilling
  !!    ccohort_hydr%flc_min_aroot(j) = ccohort_hydr%flc_aroot(j) - &
  !!          (ccohort_hydr%flc_aroot(j) - ccohort_hydr%flc_min_aroot(j))*exp(-refill_rate*dtime)
  !!end if
  !!end do
  !!end subroutine UpdateLWPMemFLCMin



end module FatesPlantHydraulicsMod
