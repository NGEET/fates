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

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!99
   ! (TODO: THE ROW WIDTH ON THIS MODULE ARE TOO LARGE. NAG COMPILERS
   !  WILL FREAK IF LINES ARE TOO LONG.  BEFORE SUBMITTING THIS TO 
   !  MASTER WE NEED TO GO THROUGH AND GET THESE LINES BELOW
   !  100 spaces (for readability), or 130 (for NAG)
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!99

   use FatesGlobals, only      : endrun => fates_endrun
   use FatesGlobals, only      : fates_log

   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : fates_huge
   use FatesConstantsMod, only : denh2o => dens_fresh_liquid_water
   use FatesConstantsMod, only : grav => grav_earth
   use FatesConstantsMod, only : ifalse, itrue
   use FatesConstantsMod, only : pi_const
   use FatesConstantsMod, only : cm2_per_m2
   use FatesConstantsMod, only : g_per_kg

   use EDParamsMod       , only : hydr_kmax_rsurf

   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_patch_type
   use EDTypesMod        , only : ed_cohort_type

   use FatesInterfaceMod  , only : bc_in_type
   use FatesInterfaceMod  , only : bc_out_type

   use FatesInterfaceMod  , only : hlm_use_planthydro

   use FatesAllometryMod, only    : bsap_allom
   use FatesAllometryMod, only    : CrownDepth
   
   use FatesHydraulicsMemMod, only: ed_site_hydr_type
   use FatesHydraulicsMemMod, only: ed_cohort_hydr_type
   use FatesHydraulicsMemMod, only: n_hypool_leaf
   use FatesHydraulicsMemMod, only: n_hypool_tot
   use FatesHydraulicsMemMod, only: n_hypool_stem
   use FatesHydraulicsMemMod, only: numLWPmem
   use FatesHydraulicsMemMod, only: n_hypool_troot
   use FatesHydraulicsMemMod, only: n_hypool_aroot
   use FatesHydraulicsMemMod, only: n_porous_media
   use FatesHydraulicsMemMod, only: nshell
   use FatesHydraulicsMemMod, only: n_hypool_ag
   use FatesHydraulicsMemMod, only: porous_media
   use FatesHydraulicsMemMod, only: cap_slp
   use FatesHydraulicsMemMod, only: cap_int
   use FatesHydraulicsMemMod, only: cap_corr
   use FatesHydraulicsMemMod, only: rwcft
   use FatesHydraulicsMemMod, only: C2B
   use FatesHydraulicsMemMod, only: InitHydraulicsDerived
   use FatesHydraulicsMemMod, only: nlevsoi_hyd_max
   use FatesHydraulicsMemMod, only: cohort_recruit_water_layer
   use FatesHydraulicsMemMod, only: recruit_water_avail_layer

   use PRTGenericMod,          only : all_carbon_elements
   use PRTGenericMod,          only : leaf_organ, fnrt_organ, sapw_organ
   use PRTGenericMod,          only : store_organ, repro_organ, struct_organ
   
   use clm_time_manager  , only : get_step_size, get_nstep

   use EDPftvarcon, only : EDPftvarcon_inst

   ! CIME Globals
   use shr_log_mod , only      : errMsg => shr_log_errMsg
   use shr_infnan_mod   , only : isnan => shr_infnan_isnan


   implicit none

   private
   integer, parameter :: van_genuchten = 1
   integer, parameter :: campbell      = 2
   integer :: iswc = campbell
   
   ! 1=leaf, 2=stem, 3=troot, 4=aroot
   ! Several of these may be better transferred to the parameter file in due time (RGK)

   integer, public :: use_ed_planthydraulics    =  1      ! 0 => use vanilla btran
                                                          ! 1 => use BC hydraulics; 
                                                          ! 2 => use CX hydraulics
   logical, public :: do_dqtopdth_leaf          = .false.  ! should a nonzero dqtopdth_leaf
                                                          ! term be applied to the plant 
                                                          ! hydraulics numerical solution?
   logical, public :: do_dyn_xylemrefill        = .false.  ! should the dynamics of xylem refilling 
                                                          ! (i.e., non-instantaneous) be considered 
                                                          ! within plant hydraulics?
   logical, public :: do_kbound_upstream        = .true.  ! should the hydraulic conductance at the 
                                                          ! boundary between nodes be taken to be a
                                                          ! function of the upstream loss of 
                                                          ! conductivity (flc)?
   logical, public :: do_growthrecruiteffects   = .true. ! should size- or root length-dependent 
                                                          ! hydraulic properties and states be 
                                                          ! updated every day when trees grow or 
                                                          ! when recruitment happens?
   logical,parameter :: debug = .false.                   !flag to report warning in hydro
							  

   character(len=*), parameter, private :: sourcefile = &
         __FILE__

   
   ! We use this parameter as the value for which we set un-initialized values
   real(r8), parameter :: un_initialized = -9.9e32_r8


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
   public :: updateSizeDepTreeHydProps
   public :: updateSizeDepTreeHydStates
   public :: initTreeHydStates
   public :: updateSizeDepRhizHydProps
   public :: updateSizeDepRhizHydStates
   public :: RestartHydrStates
   public :: SavePreviousCompartmentVolumes
   public :: SavePreviousRhizVolumes
   public :: UpdateTreeHydrNodes
   public :: UpdateTreeHydrLenVolCond
   public :: ConstrainRecruitNumber

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
    type(ed_patch_type),pointer       :: cpatch      ! current patch
    type(ed_cohort_type),pointer      :: ccohort     ! current cohort
    type(ed_cohort_hydr_type),pointer :: ccohort_hydr
    integer                           :: s           ! site loop counter

    do s = 1,nsites

       cpatch => sites(s)%oldest_patch
       do while(associated(cpatch))
             
          ccohort => cpatch%shortest
          do while(associated(ccohort))  
             
             ccohort_hydr => ccohort%co_hydr

             ! This calculates node heights
             call UpdateTreeHydrNodes(ccohort_hydr,ccohort%pft,ccohort%hite, &
                                      sites(s)%si_hydr%nlevsoi_hyd,bc_in(s))
                   
             ! This calculates volumes, lengths and max conductances
             call UpdateTreeHydrLenVolCond(ccohort,sites(s)%si_hydr%nlevsoi_hyd,bc_in(s))
             
             ! Since this is a newly initialized plant, we set the previous compartment-size
             ! equal to the ones we just calculated.
             call SavePreviousCompartmentVolumes(ccohort_hydr) 

             ! Set some generic initial values
             ccohort_hydr%refill_days      =  3.0_r8
             ccohort_hydr%lwp_mem(:)       = 0.0_r8
             ccohort_hydr%lwp_stable       = 0.0_r8
             ccohort_hydr%lwp_is_unstable  = .false.
             ccohort_hydr%flc_ag(:)        =  1.0_r8
             ccohort_hydr%flc_troot(:)     =  1.0_r8
             ccohort_hydr%flc_aroot(:)     =  1.0_r8
             ccohort_hydr%flc_min_ag(:)    =  1.0_r8
             ccohort_hydr%flc_min_troot(:)    =  1.0_r8
             ccohort_hydr%flc_min_aroot(:) =  1.0_r8
             ccohort_hydr%refill_thresh    = -0.01_r8
             ccohort_hydr%refill_days      =  3.0_r8

             ccohort => ccohort%taller
          enddo
          
          cpatch => cpatch%younger
       end do
       
       sites(s)%si_hydr%l_aroot_layer_init(:)  = un_initialized
       sites(s)%si_hydr%r_node_shell_init(:,:) = un_initialized
       sites(s)%si_hydr%v_shell_init(:,:)      = un_initialized


       ! Update static quantities related to the rhizosphere
       call UpdateSizeDepRhizVolLenCon(sites(s), bc_in(s))

       ! We update the "initial" values of the rhizosphere after
       ! the previous call to make sure that the conductances are updated
       ! Now we set the prevous to the current so that the water states
       ! are not perturbed
       call SavePreviousRhizVolumes(sites(s), bc_in(s))

    end do
    
    call UpdateH2OVeg(nsites,sites,bc_out)

    return
 end subroutine RestartHydrStates
 
 ! ====================================================================================

 subroutine initTreeHydStates(site_p, cc_p, bc_in)
    
    ! REQUIRED INPUTS:
    !
    !  csite%si_hydr%psisoi_liq_innershell(:)
    !  ccohort_hydr%z_node_troot(:)
    !  ccohort_hydr%z_node_aroot
    !  ccohort_hydr%z_node_ag
    ! 
    ! !DESCRIPTION: 
    !
    ! !USES:

    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target   :: site_p ! current cohort pointer
    type(ed_cohort_type), intent(inout), target :: cc_p ! current cohort pointer
    type(bc_in_type)    , intent(in)            :: bc_in 
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: cCohort
    type(ed_site_type), pointer   :: csite
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    integer  :: j,k,FT                     ! indices
    real(r8) :: dz
    real(r8) :: smp

    cCohort                    => cc_p
    ccohort_hydr               => cCohort%co_hydr
    csite                      => site_p 
    FT                         =  cCohort%pft

    !convert soil water contents to water potential in each soil layer and
    !assign it to the absorbing root (assume absorbing root water potential
    !in equlibrium w/ surrounding soil)
    do j=1, site_p%si_hydr%nlevsoi_hyd

       !call swcVG_psi_from_th(waterstate_inst%h2osoi_liqvol_shell(c,j,1), &
       !   watsat(c,j), watres(c,j), alpha_VG(c,j), n_VG(c,j), m_VG(c,j), l_VG(c,j), &
       !   smp)
       !ccohort_hydr%psi_aroot(j) = smp
       ccohort_hydr%psi_aroot(j) = csite%si_hydr%psisoi_liq_innershell(j)

       call th_from_psi(ft, 4, ccohort_hydr%psi_aroot(j), ccohort_hydr%th_aroot(j), csite%si_hydr, bc_in )
    end do
    
    !initialize plant water potentials at hydrostatic equilibrium (dh/dz = 0)
    !the assumption is made here that initial conditions for soil water will 
    !be in (or at least close to) hydrostatic equilibrium as well, so that
    !it doesn't matter which absorbing root layer the transporting root water
    !potential is referenced to.
    do k=1, n_hypool_troot
       dz = ccohort_hydr%z_node_troot(k) - ccohort_hydr%z_node_aroot(1)
       ccohort_hydr%psi_troot(k) = ccohort_hydr%psi_aroot(1) - 1.e-6_r8*denh2o*grav*dz 
       if (ccohort_hydr%psi_troot(k)>0.0_r8) ccohort_hydr%psi_troot(k) = -0.01_r8
       call th_from_psi(ft, 3, ccohort_hydr%psi_troot(k), ccohort_hydr%th_troot(k), csite%si_hydr, bc_in)
    end do

    !working our way up a tree, assigning water potentials that are in
    !hydrostatic equilibrium with the water potential immediately below
    dz = ccohort_hydr%z_node_ag(n_hypool_ag) - ccohort_hydr%z_node_troot(1)
    ccohort_hydr%psi_ag(n_hypool_ag) = ccohort_hydr%psi_troot(1) - 1.e-6_r8*denh2o*grav*dz
    call th_from_psi(ft, 2, ccohort_hydr%psi_ag(n_hypool_ag), ccohort_hydr%th_ag(n_hypool_ag), csite%si_hydr, bc_in)
    do k=n_hypool_ag-1, 1, -1
       dz = ccohort_hydr%z_node_ag(k) - ccohort_hydr%z_node_ag(k+1)
       ccohort_hydr%psi_ag(k) = ccohort_hydr%psi_ag(k+1) - 1.e-6_r8*denh2o*grav*dz
       if(ccohort_hydr%psi_ag(k)>0.0_r8) ccohort_hydr%psi_ag(k)= -0.01_r8
       call th_from_psi(ft, porous_media(k), ccohort_hydr%psi_ag(k), ccohort_hydr%th_ag(k), csite%si_hydr, bc_in)
    end do
    
    ccohort_hydr%lwp_mem(:)       = ccohort_hydr%psi_ag(1)   ! initializes the leaf water potential memory
    ccohort_hydr%lwp_stable       = ccohort_hydr%psi_ag(1)
    ccohort_hydr%lwp_is_unstable  = .false.                  ! inital value for leaf water potential stability flag
    ccohort_hydr%flc_ag(:)        =  1.0_r8
    ccohort_hydr%flc_troot(:)        =  1.0_r8
    ccohort_hydr%flc_aroot(:)     =  1.0_r8
    ccohort_hydr%flc_min_ag(:)    =  1.0_r8
    ccohort_hydr%flc_min_troot(:)    =  1.0_r8
    ccohort_hydr%flc_min_aroot(:) =  1.0_r8
    ccohort_hydr%refill_thresh    = -0.01_r8
    ccohort_hydr%refill_days      =  3.0_r8
    ccohort_hydr%errh2o_growturn_ag(:)    = 0.0_r8
    ccohort_hydr%errh2o_growturn_troot(:) = 0.0_r8
    ccohort_hydr%errh2o_growturn_aroot(:) = 0.0_r8
    ccohort_hydr%errh2o_pheno_ag(:)       = 0.0_r8
    ccohort_hydr%errh2o_pheno_troot(:)    = 0.0_r8
    ccohort_hydr%errh2o_pheno_aroot(:)    = 0.0_r8
    !ccohort_hydr%th_aroot_prev(:)     =  0.0_r8
    !ccohort_hydr%th_aroot_prev_ucnorr(:)=  0.0_r8
    
    !initialize cohort-level btran
    call flc_gs_from_psi(cCohort, ccohort_hydr%psi_ag(1))
    
  end subroutine initTreeHydStates

  
  ! =====================================================================================

    
  subroutine UpdateTreeHydrNodes(ccohort_hydr,ft,plant_height,nlevsoi_hyd,bc_in)
     
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
     !                     %z_node_troot(:)
     !                     %z_lower_troot(:)
     !                     %z_upper_troot(:)
     !                     %z_node_aroot(:)
     ! --------------------------------------------------------------------------------
     
     ! Arguments
     type(ed_cohort_hydr_type), intent(inout) :: ccohort_hydr 
     integer,intent(in)                       :: ft            ! plant functional type index
     real(r8), intent(in)                     :: plant_height  ! [m]
     integer,intent(in)                       :: nlevsoi_hyd   ! number of soil hydro layers
     type(bc_in_type)       , intent(in)      :: bc_in         ! Boundary Conditions
    
     
     ! Locals
     
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

     roota                      =  EDPftvarcon_inst%roota_par(ft)
     rootb                      =  EDPftvarcon_inst%rootb_par(ft)

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
     
     ! Transporting Root Nodes
     ! in special case where n_hypool_troot = 1, the node depth of the single troot pool
     ! is the depth at which 50% total root distribution is attained
     dcumul_rf                  = 1._r8/real(n_hypool_troot,r8)
     
     do k=1,n_hypool_troot
        cumul_rf                = dcumul_rf*real(k,r8)
        call bisect_rootfr(roota, rootb, 0._r8, 1.E10_r8, &
              0.001_r8, 0.001_r8, cumul_rf, z_cumul_rf)
        z_cumul_rf =  min(z_cumul_rf, abs(bc_in%zi_sisl(nlevsoi_hyd)))
        ccohort_hydr%z_lower_troot(k)   = -z_cumul_rf
        call bisect_rootfr(roota, rootb, 0._r8, 1.E10_r8, &
              0.001_r8, 0.001_r8, cumul_rf-0.5_r8*dcumul_rf, z_cumul_rf)
        z_cumul_rf =  min(z_cumul_rf, abs(bc_in%zi_sisl(nlevsoi_hyd)))
        ccohort_hydr%z_node_troot(k)    = -z_cumul_rf
        call bisect_rootfr(roota, rootb, 0._r8, 1.E10_r8, &
              0.001_r8, 0.001_r8, cumul_rf-1.0_r8*dcumul_rf+1.E-10_r8, z_cumul_rf)
        z_cumul_rf =  min(z_cumul_rf, abs(bc_in%zi_sisl(nlevsoi_hyd)))
        ccohort_hydr%z_upper_troot(k)   = -z_cumul_rf
     enddo
     
     
     ! Absorbing root depth
     ccohort_hydr%z_node_aroot(1:nlevsoi_hyd) = -bc_in%z_sisl(1:nlevsoi_hyd)
     

     ! Shouldn't this be updating the upper and lower values as well?
     ! (RGK 12-2018)
     if(nlevsoi_hyd == 1) then
        ccohort_hydr%z_node_troot(:)    = ccohort_hydr%z_node_aroot(nlevsoi_hyd)
     end if
     
     return
  end subroutine UpdateTreeHydrNodes

  ! =====================================================================================
  
  subroutine SavePreviousCompartmentVolumes(ccohort_hydr)

    type(ed_cohort_hydr_type),intent(inout) :: ccohort_hydr
    
    
    ! Saving the current compartment volumes into an "initial" save-space
    ! allows us to see how the compartments change size when plants
    ! change size and effect water contents
    
    ccohort_hydr%v_ag_init(:)          =  ccohort_hydr%v_ag(:)
    ccohort_hydr%v_troot_init(:)       =  ccohort_hydr%v_troot(:)
    ccohort_hydr%v_aroot_layer_init(:) =  ccohort_hydr%v_aroot_layer(:)
    
    return
  end subroutine SavePreviousCompartmentVolumes
  
  ! =====================================================================================
  
  subroutine updateSizeDepTreeHydProps(currentSite,ccohort,bc_in)


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
    integer                            :: nlevsoi_hyd             ! Number of total soil layers
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    integer                            :: ft

    nlevsoi_hyd                =  currentSite%si_hydr%nlevsoi_hyd
    ccohort_hydr               => ccohort%co_hydr
    ft                         =  ccohort%pft

    ! Save the current vegetation compartment volumes into
    ! a save space so that it can be compared with the updated quantity.
    
    call SavePreviousCompartmentVolumes(ccohort_hydr)
    
    ! This updates all of the z_node positions
    call UpdateTreeHydrNodes(ccohort_hydr,ft,ccohort%hite,nlevsoi_hyd,bc_in)
    
    ! This updates plant compartment volumes, lengths and 
    ! maximum conductances. Make sure for already
    ! initialized vegetation, that SavePreviousCompartment
    ! volumes, and UpdateTreeHydrNodes is called prior to this.
    
    call UpdateTreeHydrLenVolCond(ccohort,nlevsoi_hyd,bc_in)
    
    
  end subroutine updateSizeDepTreeHydProps

  ! =====================================================================================

  subroutine UpdateTreeHydrLenVolCond(ccohort,nlevsoi_hyd,bc_in)
    
      ! -----------------------------------------------------------------------------------
      ! This subroutine calculates three attributes of a plant:
      ! 1) the volumes of storage compartments in the plants
      ! 2) the lenghts of the organs 
      ! 3) the conductances
      ! These and are not dependent on the hydraulic state of the
      ! plant, it is more about the structural characteristics and how much biomass
      ! is present in the different tissues.
      !
      ! Inputs, plant geometries, plant carbon pools, z_node values
      !
      ! -----------------------------------------------------------------------------------

      ! Arguments
      type(ed_cohort_type),intent(inout)  :: ccohort
      integer,intent(in)                  :: nlevsoi_hyd   ! number of soil hydro layers
      type(bc_in_type)       , intent(in) :: bc_in       ! Boundary Conditions
      
      type(ed_cohort_hydr_type),pointer :: ccohort_hydr     ! Plant hydraulics structure
      integer  :: j,k
      integer  :: ft                           ! Plant functional type index 
      real(r8) :: roota                        ! root profile parameter a zeng2001_crootfr
      real(r8) :: rootb                        ! root profile parameter b zeng2001_crootfr
      real(r8) :: leaf_c                       ! Current amount of leaf carbon in the plant                            [kg]
      real(r8) :: fnrt_c                       ! Current amount of fine-root carbon in the plant                       [kg]
      real(r8) :: sapw_c                       ! Current amount of sapwood carbon in the plant                         [kg]
      real(r8) :: struct_c                     ! Current amount of structural carbon in the plant                      [kg]
      real(r8) :: b_canopy_carb                ! total leaf (canopy) biomass in carbon units                           [kgC/indiv]
      real(r8) :: b_canopy_biom                ! total leaf (canopy) biomass in dry wt units                           [kg/indiv]
      real(r8) :: b_woody_carb                 ! total woody biomass in carbon units                                   [kgC/indiv]
      real(r8) :: b_woody_bg_carb              ! belowground woody biomass in carbon units                             [kgC/indiv]
      real(r8) :: b_stem_carb                  ! aboveground stem biomass in carbon units                              [kgC/indiv]
      real(r8) :: b_stem_biom                  ! aboveground stem biomass in dry wt units                              [kg/indiv]
      real(r8) :: b_bg_carb                    ! belowground biomass (coarse + fine roots) in carbon units             [kgC/indiv]
      real(r8) :: b_tot_carb                   ! total individual biomass in carbon units                              [kgC/indiv]
      real(r8) :: v_stem                       ! aboveground stem volume                                               [m3/indiv]
      real(r8) :: z_stem                       ! the height of the plants stem below crown [m]
      real(r8) :: sla                          ! specific leaf area                                                    [cm2/g]
      real(r8) :: v_canopy                     ! total leaf (canopy) volume                                            [m3/indiv]
      real(r8) :: denleaf                      ! leaf dry mass per unit fresh leaf volume                              [kg/m3]     
      real(r8) :: a_sapwood                    ! sapwood area                                                          [m2]
      real(r8) :: a_sapwood_target             ! sapwood cross-section area at reference height, at target biomass     [m2]
      real(r8) :: bsw_target                   ! sapwood carbon, at target                                             [kgC]
      real(r8) :: v_sapwood                    ! sapwood volume                                                        [m3]
      real(r8) :: b_troot_carb                 ! transporting root biomass in carbon units                             [kgC/indiv]
      real(r8) :: b_troot_biom                 ! transporting root biomass in dry wt units                             [kg/indiv]
      real(r8) :: v_troot                      ! transporting root volume                                              [m3/indiv]
      real(r8) :: rootfr                       ! mass fraction of roots in each layer                                  [kg/kg]
      real(r8) :: crown_depth                  ! Depth of the plant's crown [m]
      real(r8) :: kmax_node1_nodekplus1(n_hypool_ag) ! cumulative kmax, petiole to node k+1, 
                                                     ! conduit taper effects excluded   [kg s-1 MPa-1]
      real(r8) :: kmax_node1_lowerk(n_hypool_ag)     ! cumulative kmax, petiole to upper boundary of node k, 
                                                     ! conduit taper effects excluded   [kg s-1 MPa-1]
      real(r8) :: chi_node1_nodekplus1(n_hypool_ag)  ! ratio of cumulative kmax with taper effects 
                                                     ! included to that without  [-]
      real(r8) :: chi_node1_lowerk(n_hypool_ag)      ! ratio of cumulative kmax with taper effects 
                                                     ! included to that without  [-]
      real(r8) :: dz_node1_nodekplus1                ! cumulative distance between canopy 
                                                     ! node and node k + 1                [m]
      real(r8) :: dz_node1_lowerk                    ! cumulative distance between canopy 
                                                     ! node and upper boundary of node k  [m]
      real(r8) :: kmax_treeag_tot                    ! total stem (petiole to transporting root node) 
                                                     ! hydraulic conductance        [kg s-1 MPa-1]
      real(r8) :: kmax_tot                           ! total tree (leaf to root tip) 
                                                     ! hydraulic conductance        [kg s-1 MPa-1]
      
      real(r8),parameter :: taper_exponent = 1._r8/3._r8 ! Savage et al. (2010) xylem taper exponent [-]

      ccohort_hydr => ccohort%co_hydr
      ft           = ccohort%pft

      leaf_c   = ccohort%prt%GetState(leaf_organ, all_carbon_elements)
      sapw_c   = ccohort%prt%GetState(sapw_organ, all_carbon_elements)
      fnrt_c   = ccohort%prt%GetState(fnrt_organ, all_carbon_elements)
      struct_c = ccohort%prt%GetState(struct_organ, all_carbon_elements)

      roota    =  EDPftvarcon_inst%roota_par(ft)
      rootb    =  EDPftvarcon_inst%rootb_par(ft)

      !roota                      =  4.372_r8                           ! TESTING: deep (see Zeng 2001 Table 1)
      !rootb                      =  0.978_r8                           ! TESTING: deep (see Zeng 2001 Table 1)
      !roota                      =  8.992_r8                          ! TESTING: shallow (see Zeng 2001 Table 1)
      !rootb                      =  8.992_r8                          ! TESTING: shallow (see Zeng 2001 Table 1)
      
      if(leaf_c > 0._r8) then


         ! ------------------------------------------------------------------------------
         ! Part 1.  Set the volumes of the leaf, stem and root compartments
         !          and lenghts of the roots
         ! ------------------------------------------------------------------------------

         b_woody_carb               = sapw_c + struct_c
         b_woody_bg_carb            = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(ft)) * b_woody_carb
         b_tot_carb                 = sapw_c + struct_c + leaf_c + fnrt_c
         b_canopy_carb              = leaf_c
         b_bg_carb                  = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(ft)) * b_tot_carb
         b_canopy_biom              = b_canopy_carb * C2B
         
         ! NOTE: SLATOP currently does not use any vertical scaling functions
         ! but that may not be so forever. ie sla = slatop (RGK-082017)
         ! m2/gC * cm2/m2 -> cm2/gC
         sla                        = EDPftvarcon_inst%slatop(ft) * cm2_per_m2 
         
         ! empirical regression data from leaves at Caxiuana (~ 8 spp)
         denleaf                    = -2.3231_r8*sla/C2B + 781.899_r8    
         v_canopy                   = b_canopy_biom / denleaf

         ccohort_hydr%v_ag(1:n_hypool_leaf) = v_canopy / real(n_hypool_leaf,r8)

 
         b_stem_carb  = b_tot_carb - b_bg_carb - b_canopy_carb
         b_stem_biom  = b_stem_carb * C2B                               ! kg DM

         !BOC...may be needed for testing/comparison w/ v_sapwood 
         !    kg  / ( g cm-3 * cm3/m3 * kg/g ) -> m3    
         v_stem       = b_stem_biom / (EDPftvarcon_inst%wood_density(ft)*1.e3_r8  ) 

         ! calculate the sapwood cross-sectional area
         call bsap_allom(ccohort%dbh,ccohort%pft,ccohort%canopy_trim,a_sapwood_target,bsw_target)
         a_sapwood = a_sapwood_target
         
         ! Alternative ways to calculate sapwood cross section
         ! or ....
         ! a_sapwood = a_sapwood_target * ccohort%bsw / bsw_target
         
         !     a_sapwood    = a_leaf_tot / EDPftvarcon_inst%allom_latosa_int(ft)*1.e-4_r8 
         !      m2 sapwood = m2 leaf * cm2 sapwood/m2 leaf *1.0e-4m2
         ! or ...
         !a_sapwood    = a_leaf_tot / ( 0.001_r8 + 0.025_r8 * ccohort%hite ) * 1.e-4_r8

         call CrownDepth(ccohort%hite,crown_depth)
         z_stem       = ccohort%hite - crown_depth
         v_sapwood    = a_sapwood * z_stem
         ccohort_hydr%v_ag(n_hypool_leaf+1:n_hypool_ag) = v_sapwood / n_hypool_stem


         ! Determine belowground biomass as a function of total (sapwood, heartwood, 
         ! leaf, fine root) biomass then subtract out the fine root biomass to get 
         ! coarse (transporting) root biomass
         
         b_troot_carb               = b_woody_bg_carb   
         b_troot_biom               = b_troot_carb * C2B 
         v_troot                    = b_troot_biom / (EDPftvarcon_inst%wood_density(ft)*1.e3_r8)

         !! BOC not sure if/how we should multiply this by the sapwood fraction
         ccohort_hydr%v_troot(:)    = v_troot / n_hypool_troot    

     
         ! Estimate absorbing root total length (all layers)
         ! ------------------------------------------------------------------------------
         ccohort_hydr%l_aroot_tot        = fnrt_c*C2B*EDPftvarcon_inst%hydr_srl(ft)

         ! Estimate absorbing root volume (all layers)
         ! ------------------------------------------------------------------------------
         ccohort_hydr%v_aroot_tot        = pi_const * (EDPftvarcon_inst%hydr_rs2(ft)**2._r8) * &
                                           ccohort_hydr%l_aroot_tot

         
         ! Partition the total absorbing root lengths and volumes into the active soil layers
         ! ------------------------------------------------------------------------------
         if(nlevsoi_hyd == 1) then
            ccohort_hydr%l_aroot_layer(nlevsoi_hyd) = ccohort_hydr%l_aroot_tot
            ccohort_hydr%v_aroot_layer(nlevsoi_hyd) = ccohort_hydr%v_aroot_tot
         else
            do j=1,nlevsoi_hyd
               if(j == 1) then
                  rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j))
               else
                  rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j)) - &
                       zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j-1))
               end if
               ccohort_hydr%l_aroot_layer(j)   = rootfr*ccohort_hydr%l_aroot_tot
               ccohort_hydr%v_aroot_layer(j)   = rootfr*ccohort_hydr%v_aroot_tot
            end do
         end if
         

         ! ------------------------------------------------------------------------------
         ! Part II. Set maximum (size-dependent) hydraulic conductances
         ! ------------------------------------------------------------------------------

         ! first estimate cumulative (petiole to node k) conductances 
         ! without taper as well as the chi taper function
         
         do k=n_hypool_leaf,n_hypool_ag
            dz_node1_lowerk          = ccohort_hydr%z_node_ag(n_hypool_leaf) &
                                     - ccohort_hydr%z_lower_ag(k)
            if(k < n_hypool_ag) then
               dz_node1_nodekplus1   = ccohort_hydr%z_node_ag(n_hypool_leaf) &
                                     - ccohort_hydr%z_node_ag(k+1)
            else
               dz_node1_nodekplus1   = ccohort_hydr%z_node_ag(n_hypool_leaf) &
                                     - ccohort_hydr%z_node_troot(1)
            end if
            kmax_node1_nodekplus1(k) = EDPftvarcon_inst%hydr_kmax_node(ft,2) * a_sapwood / dz_node1_nodekplus1
            kmax_node1_lowerk(k)     = EDPftvarcon_inst%hydr_kmax_node(ft,2) * a_sapwood / dz_node1_lowerk
            chi_node1_nodekplus1(k)  = xylemtaper(taper_exponent, dz_node1_nodekplus1)
            chi_node1_lowerk(k)      = xylemtaper(taper_exponent, dz_node1_lowerk)
            if(.not.do_kbound_upstream) then
               if(crown_depth == 0._r8) then 
                  write(fates_log(),*) 'do_kbound_upstream requires a nonzero canopy depth '
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
            end if
         enddo
         
         
         ! then calculate the conductances at node boundaries as the difference of cumulative conductances
         do k=n_hypool_leaf,n_hypool_ag
            if(k == n_hypool_leaf) then
               ccohort_hydr%kmax_bound(k)    = kmax_node1_nodekplus1(k)  * chi_node1_nodekplus1(k)
               ccohort_hydr%kmax_lower(k)    = kmax_node1_lowerk(k)      * chi_node1_lowerk(k)
            else
               ccohort_hydr%kmax_bound(k)    = ( 1._r8/(kmax_node1_nodekplus1(k)  *chi_node1_nodekplus1(k)  ) - &
                    1._r8/(kmax_node1_nodekplus1(k-1)*chi_node1_nodekplus1(k-1))     ) ** (-1._r8)
               ccohort_hydr%kmax_lower(k)    = ( 1._r8/(kmax_node1_lowerk(k)      *chi_node1_lowerk(k)      ) - &
                    1._r8/(kmax_node1_nodekplus1(k-1)*chi_node1_nodekplus1(k-1))     ) ** (-1._r8)
            end if
            if(k < n_hypool_ag) then
               ccohort_hydr%kmax_upper(k+1)  = ( 1._r8/(kmax_node1_nodekplus1(k)  *chi_node1_nodekplus1(k)  ) - &
                    1._r8/(kmax_node1_lowerk(k)      *chi_node1_lowerk(k)      )     ) ** (-1._r8)
            else if(k == n_hypool_ag) then
               ccohort_hydr%kmax_upper_troot = ( 1._r8/(kmax_node1_nodekplus1(k)  *chi_node1_nodekplus1(k)  ) - &
                    1._r8/(kmax_node1_lowerk(k)      *chi_node1_lowerk(k)      )     ) ** (-1._r8)
            end if

            !!!!!!!!!! FOR TESTING ONLY
            !ccohort_hydr%kmax_bound(:) = 0.02_r8   ! Diurnal lwp variation in coldstart: -0.1 MPa
            ! Diurnal lwp variation in large-tree (50cmDBH) coldstart: less than -0.01 MPa
            !ccohort_hydr%kmax_bound(:) = 0.0016_r8 ! Diurnal lwp variation in coldstart: -0.8 - 1.0 MPa
            ! Diurnal lwp variation in large-tree (50cmDBH) coldstart: -1.5 - 2.0 MPa [seemingly unstable]
            !ccohort_hydr%kmax_bound(:) = 0.0008_r8 ! Diurnal lwp variation in coldstart: -1.5 - 2.0 MPa
            ! Diurnal lwp variation in large-tree (50cmDBH) coldstart: -2.0 - 3.0 MPa [seemingly unstable]
            !ccohort_hydr%kmax_bound(:) = 0.0005_r8 ! Diurnal lwp variation in coldstart: -2.0 - 3.0 MPa and one -5 MPa outlier
            ! Diurnal lwp variation in large-tree (50cmDBH) coldstart: -3.0 - 4.0 MPa and one -10 MPa outlier [Unstable]
            !!!!!!!!!!

         enddo

         ! finally, estimate the remaining tree conductance belowground as a residual
         kmax_treeag_tot = sum(1._r8/ccohort_hydr%kmax_bound(n_hypool_leaf:n_hypool_ag))**(-1._r8)
         kmax_tot        = EDPftvarcon_inst%hydr_rfrac_stem(ft) * kmax_treeag_tot
         ccohort_hydr%kmax_treebg_tot      = ( 1._r8/kmax_tot - 1._r8/kmax_treeag_tot ) ** (-1._r8)
         
         if(nlevsoi_hyd == 1) then
            ccohort_hydr%kmax_treebg_layer(:) = ccohort_hydr%kmax_treebg_tot * &
                                                ccohort%patchptr%rootfr_ft(ft,:)
         else
            do j=1,nlevsoi_hyd
               if(j == 1) then
                  rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j))
               else
                  rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j)) - &
                       zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j-1))
               end if
               ccohort_hydr%kmax_treebg_layer(j) = rootfr*ccohort_hydr%kmax_treebg_tot
            end do
         end if
      end if !check for bleaf
      
    end subroutine UpdateTreeHydrLenVolCond


  ! =====================================================================================
  subroutine updateSizeDepTreeHydStates(currentSite,ccohort)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    use FatesUtilsMod  , only : check_var_real
    use EDTypesMod     , only : dump_cohort
    use EDTypesMod     , only : AREA
    
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
    real(r8) :: th_troot_uncorr(n_hypool_troot) ! uncorrected transporting root water content[m3 m-3]
    real(r8) :: th_aroot_uncorr(currentSite%si_hydr%nlevsoi_hyd)    ! uncorrected absorbing root water content[m3 m-3] 
    real(r8), parameter :: small_theta_num = 1.e-7_r8  ! avoids theta values equalling thr or ths         [m3 m-3]
    integer :: nstep !number of time steps
    !-----------------------------------------------------------------------

    ccohort_hydr => ccohort%co_hydr
    FT      =  cCohort%pft
    
    ! MAYBE ADD A NAN CATCH?  If updateSizeDepTreeHydProps() was not called twice prior to the first
    ! time this routine is called for a new cohort, then v_ag_init(k) will be a nan.
    ! It should be ok, but may be vulnerable if code is changed (RGK 02-2017)

    ! UPDATE WATER CONTENTS (assume water for growth comes from within tissue itself -- apply water mass conservation)
    do k=1,n_hypool_ag
       th_ag_uncorr(k)    = ccohort_hydr%th_ag(k)   * &
                            ccohort_hydr%v_ag_init(k) /ccohort_hydr%v_ag(k)
       ccohort_hydr%th_ag(k) = constrain_water_contents(th_ag_uncorr(k), small_theta_num, ft, k)
    enddo
    do k=1,n_hypool_troot
       th_troot_uncorr(k) = ccohort_hydr%th_troot(k)  	 * &
                                  ccohort_hydr%v_troot_init(k) /ccohort_hydr%v_troot(k)
       ccohort_hydr%th_troot(k) = constrain_water_contents(th_troot_uncorr(k), small_theta_num, ft, 3)
    enddo
    do j=1,currentSite%si_hydr%nlevsoi_hyd
       th_aroot_uncorr(j) = ccohort_hydr%th_aroot(j) * &
                                  ccohort_hydr%v_aroot_layer_init(j)/ccohort_hydr%v_aroot_layer(j)
       ccohort_hydr%th_aroot(j) = constrain_water_contents(th_aroot_uncorr(j), small_theta_num, ft, 4)
       ccohort_hydr%errh2o_growturn_aroot(j) = ccohort_hydr%th_aroot(j) - th_aroot_uncorr(j)
       !call check_var_real(ccohort_hydr%errh2o_growturn_aroot(j),'ccohort_hydr%errh2o_growturn_aroot(j)',err_code)
       !if ((abs(ccohort_hydr%errh2o_growturn_aroot(j)) > 1.0_r8) .or. &
       !     err_code == 1 .or. err_code == 10) then
       !  call dump_cohort(cCohort)
       !end if
    enddo
    
    ! Storing mass balance error
    ! + means water created; - means water destroyed
    ccohort_hydr%errh2o_growturn_ag(:)    = ccohort_hydr%th_ag(:)    - th_ag_uncorr(:)
    ccohort_hydr%errh2o_growturn_troot(:) = ccohort_hydr%th_troot(:) - th_troot_uncorr(:)
    csite_hydr =>currentSite%si_hydr
    csite_hydr%h2oveg_growturn_err = csite_hydr%h2oveg_growturn_err + &
                    (sum(ccohort_hydr%errh2o_growturn_ag(:)*ccohort_hydr%v_ag(:))      + &
                     sum(ccohort_hydr%errh2o_growturn_troot(:)*ccohort_hydr%v_troot(:))   + &
                     sum(ccohort_hydr%errh2o_growturn_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                     denh2o*cCohort%n/AREA
    
    ! UPDATES OF WATER POTENTIALS ARE DONE PRIOR TO RICHARDS' SOLUTION WITHIN FATESPLANTHYDRAULICSMOD.F90
    

  end subroutine updateSizeDepTreeHydStates

! =====================================================================================
  
  function constrain_water_contents(th_uncorr, delta, ft, k) result(th_corr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: th_uncorr ! uncorrected water content (m3 m-3)
    real(r8) , intent(in) :: delta
    integer , intent(in)  :: ft
    integer , intent(in)  :: k
    !
    ! !Local:
    real(r8) :: thr                    ! residual water content (m3 m-3)
    real(r8) :: ths                    ! saturated water content (m3 m-3)
    !
    ! !RESULT
    real(r8) :: th_corr                ! corrected water content
    !
    !------------------------------------------------------------------------
    ths     = EDPftvarcon_inst%hydr_thetas_node(ft,porous_media(k))
    thr     = ths * EDPftvarcon_inst%hydr_resid_node(ft,porous_media(k))
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


     ! BC...PLANT HYDRAULICS - "constants" that change with size. 
     ! Heights are referenced to soil surface (+ = above; - = below)
     ncohort_hydr%z_node_ag          = ocohort_hydr%z_node_ag
     ncohort_hydr%z_node_troot       = ocohort_hydr%z_node_troot
     ncohort_hydr%z_upper_ag         = ocohort_hydr%z_upper_ag
     ncohort_hydr%z_upper_troot      = ocohort_hydr%z_upper_troot
     ncohort_hydr%z_lower_ag         = ocohort_hydr%z_lower_ag
     ncohort_hydr%z_lower_troot      = ocohort_hydr%z_lower_troot
     ncohort_hydr%kmax_upper         = ocohort_hydr%kmax_upper
     ncohort_hydr%kmax_lower         = ocohort_hydr%kmax_lower
     ncohort_hydr%kmax_upper_troot   = ocohort_hydr%kmax_upper_troot
     ncohort_hydr%kmax_bound         = ocohort_hydr%kmax_bound
     ncohort_hydr%kmax_treebg_tot    = ocohort_hydr%kmax_treebg_tot
     ncohort_hydr%v_ag_init          = ocohort_hydr%v_ag_init
     ncohort_hydr%v_ag               = ocohort_hydr%v_ag
     ncohort_hydr%v_troot_init       = ocohort_hydr%v_troot_init
     ncohort_hydr%v_troot            = ocohort_hydr%v_troot
     ncohort_hydr%v_aroot_tot        = ocohort_hydr%v_aroot_tot
     ncohort_hydr%l_aroot_tot        = ocohort_hydr%l_aroot_tot
     ! quantities indexed by soil layer
     ncohort_hydr%z_node_aroot       = ocohort_hydr%z_node_aroot
     ncohort_hydr%kmax_treebg_layer  = ocohort_hydr%kmax_treebg_layer
     ncohort_hydr%v_aroot_layer_init = ocohort_hydr%v_aroot_layer_init
     ncohort_hydr%v_aroot_layer      = ocohort_hydr%v_aroot_layer
     ncohort_hydr%l_aroot_layer      = ocohort_hydr%l_aroot_layer
     
     ! BC PLANT HYDRAULICS - state variables
     ncohort_hydr%th_ag              = ocohort_hydr%th_ag
     ncohort_hydr%th_troot           = ocohort_hydr%th_troot
     ncohort_hydr%psi_ag             = ocohort_hydr%psi_ag
     ncohort_hydr%psi_troot          = ocohort_hydr%psi_troot
     ncohort_hydr%flc_ag             = ocohort_hydr%flc_ag
     ncohort_hydr%flc_troot          = ocohort_hydr%flc_troot
     ncohort_hydr%flc_min_ag         = ocohort_hydr%flc_min_ag

     ncohort_hydr%flc_min_troot      = ocohort_hydr%flc_min_troot

     !refilling status--these are constants are should be moved the fates parameter file(Chonggang XU)
     ncohort_hydr%refill_thresh      = ocohort_hydr%refill_thresh
     ncohort_hydr%refill_days        = ocohort_hydr%refill_days
     ncohort_hydr%btran              = ocohort_hydr%btran

     ncohort_hydr%lwp_mem            = ocohort_hydr%lwp_mem
     ncohort_hydr%lwp_stable         = ocohort_hydr%lwp_stable
     ncohort_hydr%lwp_is_unstable    = ocohort_hydr%lwp_is_unstable
     ncohort_hydr%supsub_flag        = ocohort_hydr%supsub_flag

     ncohort_hydr%iterh1             = ocohort_hydr%iterh1
     ncohort_hydr%iterh2             = ocohort_hydr%iterh2
     ncohort_hydr%errh2o             = ocohort_hydr%errh2o
     ncohort_hydr%errh2o_growturn_ag = ocohort_hydr%errh2o_growturn_ag






     ncohort_hydr%errh2o_pheno_ag    = ocohort_hydr%errh2o_pheno_ag



     ncohort_hydr%errh2o_growturn_troot = ocohort_hydr%errh2o_growturn_troot
     ncohort_hydr%errh2o_pheno_troot    = ocohort_hydr%errh2o_pheno_troot
     ! quantities indexed by soil layer
     ncohort_hydr%th_aroot              = ocohort_hydr%th_aroot
     !ncohort_hydr%th_aroot_prev        = ocohort_hydr%th_aroot_prev
     !ncohort_hydr%th_aroot_prev_uncorr = ocohort_hydr%th_aroot_prev_uncorr
     ncohort_hydr%psi_aroot             = ocohort_hydr%psi_aroot
     ncohort_hydr%flc_aroot             = ocohort_hydr%flc_aroot
     ncohort_hydr%flc_min_aroot         = ocohort_hydr%flc_min_aroot
     
     ncohort_hydr%errh2o_growturn_aroot = ocohort_hydr%errh2o_growturn_aroot
     ncohort_hydr%errh2o_pheno_aroot    = ocohort_hydr%errh2o_pheno_aroot
     
     ! BC PLANT HYDRAULICS - flux terms
     ncohort_hydr%qtop_dt            = ocohort_hydr%qtop_dt
     ncohort_hydr%dqtopdth_dthdt     = ocohort_hydr%dqtopdth_dthdt

     ncohort_hydr%sapflow            = ocohort_hydr%sapflow
     ncohort_hydr%rootuptake         = ocohort_hydr%rootuptake
     ncohort_hydr%rootuptake01       = ocohort_hydr%rootuptake01
     ncohort_hydr%rootuptake02       = ocohort_hydr%rootuptake02
     ncohort_hydr%rootuptake03       = ocohort_hydr%rootuptake03
     ncohort_hydr%rootuptake04       = ocohort_hydr%rootuptake04
     ncohort_hydr%rootuptake05       = ocohort_hydr%rootuptake05
     ncohort_hydr%rootuptake06       = ocohort_hydr%rootuptake06
     ncohort_hydr%rootuptake07       = ocohort_hydr%rootuptake07
     ncohort_hydr%rootuptake08       = ocohort_hydr%rootuptake08
     ncohort_hydr%rootuptake09       = ocohort_hydr%rootuptake09
     ncohort_hydr%rootuptake10       = ocohort_hydr%rootuptake10
     
     ncohort_hydr%is_newly_recruited  = ocohort_hydr%is_newly_recruited

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

     site_hydr => currentSite%si_hydr

     ccohort_hydr => currentCohort%co_hydr
     ncohort_hydr => nextCohort%co_hydr

     ccohort_hydr%th_ag(:)    = (currentCohort%n*ccohort_hydr%th_ag(:)    + &
                                nextCohort%n*ncohort_hydr%th_ag(:))/newn
     ccohort_hydr%th_troot(:)    = (currentCohort%n*ccohort_hydr%th_troot(:)    + &
                                nextCohort%n*ncohort_hydr%th_troot(:))/newn
     ccohort_hydr%th_aroot(:) = (currentCohort%n*ccohort_hydr%th_aroot(:) + &
                                nextCohort%n*ncohort_hydr%th_aroot(:))/newn
     ccohort_hydr%supsub_flag = 0._r8
     ccohort_hydr%iterh1      = 0._r8
     ccohort_hydr%iterh2      = 0._r8

     do k=1,n_hypool_ag
        call psi_from_th(currentCohort%pft, porous_media(k), ccohort_hydr%th_ag(k), &
                         ccohort_hydr%psi_ag(k), site_hydr, bc_in)
        call flc_from_psi(currentCohort%pft, porous_media(k), ccohort_hydr%psi_ag(k), &
                         ccohort_hydr%flc_ag(k), site_hydr, bc_in) 
     end do
     do k=n_hypool_ag+1,n_hypool_ag+n_hypool_troot
        call psi_from_th(currentCohort%pft, 3, ccohort_hydr%th_troot(k-n_hypool_ag), &
                         ccohort_hydr%psi_troot(k-n_hypool_ag), site_hydr, bc_in)
        call flc_from_psi(currentCohort%pft, 3, ccohort_hydr%psi_troot(k-n_hypool_ag), &
                         ccohort_hydr%flc_troot(k-n_hypool_ag), site_hydr, bc_in)
     end do
     do j=1,site_hydr%nlevsoi_hyd
        call psi_from_th(currentCohort%pft, 4, ccohort_hydr%th_aroot(j), &
                         ccohort_hydr%psi_aroot(j), site_hydr, bc_in)
        call flc_from_psi(currentCohort%pft, 4, ccohort_hydr%psi_aroot(j), &
                          ccohort_hydr%flc_aroot(j), site_hydr, bc_in)
     end do
     call flc_gs_from_psi(currentCohort, ccohort_hydr%psi_ag(1))
     ccohort_hydr%qtop_dt        = (currentCohort%n*ccohort_hydr%qtop_dt        + &
                                    nextCohort%n*ncohort_hydr%qtop_dt)/newn
     ccohort_hydr%dqtopdth_dthdt = (currentCohort%n*ccohort_hydr%dqtopdth_dthdt + &
                                    nextCohort%n*ncohort_hydr%dqtopdth_dthdt)/newn
     ccohort_hydr%sapflow        = (currentCohort%n*ccohort_hydr%sapflow        + &
                                    nextCohort%n*ncohort_hydr%sapflow)/newn
     ccohort_hydr%rootuptake     = (currentCohort%n*ccohort_hydr%rootuptake     + &
                                    nextCohort%n*ncohort_hydr%rootuptake)/newn
     ccohort_hydr%rootuptake01   = (currentCohort%n*ccohort_hydr%rootuptake01   + &
                                    nextCohort%n*ncohort_hydr%rootuptake01)/newn
     ccohort_hydr%rootuptake02   = (currentCohort%n*ccohort_hydr%rootuptake02   + &
                                    nextCohort%n*ncohort_hydr%rootuptake02)/newn
     ccohort_hydr%rootuptake03   = (currentCohort%n*ccohort_hydr%rootuptake03   + &
                                    nextCohort%n*ncohort_hydr%rootuptake03)/newn
     ccohort_hydr%rootuptake04   = (currentCohort%n*ccohort_hydr%rootuptake04   + &
                                    nextCohort%n*ncohort_hydr%rootuptake04)/newn
     ccohort_hydr%rootuptake05   = (currentCohort%n*ccohort_hydr%rootuptake05   + &
                                    nextCohort%n*ncohort_hydr%rootuptake05)/newn
     ccohort_hydr%rootuptake06   = (currentCohort%n*ccohort_hydr%rootuptake06   + &
                                    nextCohort%n*ncohort_hydr%rootuptake06)/newn
     ccohort_hydr%rootuptake07   = (currentCohort%n*ccohort_hydr%rootuptake07   + &
                                    nextCohort%n*ncohort_hydr%rootuptake07)/newn
     ccohort_hydr%rootuptake08   = (currentCohort%n*ccohort_hydr%rootuptake08   + &
                                    nextCohort%n*ncohort_hydr%rootuptake08)/newn
     ccohort_hydr%rootuptake09   = (currentCohort%n*ccohort_hydr%rootuptake09   + &
                                    nextCohort%n*ncohort_hydr%rootuptake09)/newn
     ccohort_hydr%rootuptake10   = (currentCohort%n*ccohort_hydr%rootuptake10   + &
                                    nextCohort%n*ncohort_hydr%rootuptake10)/newn

     ccohort_hydr%lwp_mem(:)       = ccohort_hydr%psi_ag(1)
     ccohort_hydr%lwp_stable       = ccohort_hydr%psi_ag(1)
     ccohort_hydr%lwp_is_unstable  = .false.
     ccohort_hydr%flc_min_ag(:)    = (currentCohort%n*ccohort_hydr%flc_min_ag(:)    + &
                                      nextCohort%n*ncohort_hydr%flc_min_ag(:))/newn
     ccohort_hydr%flc_min_troot(:) = (currentCohort%n*ccohort_hydr%flc_min_troot(:) + &
                                      nextCohort%n*ncohort_hydr%flc_min_troot(:))/newn
     ccohort_hydr%flc_min_aroot(:) = (currentCohort%n*ccohort_hydr%flc_min_aroot(:) + &
                                      nextCohort%n*ncohort_hydr%flc_min_aroot(:))/newn
     
     ! need to be migrated to parmeter file (BOC 07/24/2018)
     ccohort_hydr%refill_thresh            = -0.01_r8
     ccohort_hydr%refill_days              =  3.0_r8
     
     ccohort_hydr%errh2o                   = (currentCohort%n*ccohort_hydr%errh2o                   + &
                                              nextCohort%n*ncohort_hydr%errh2o)/newn
     ccohort_hydr%errh2o_growturn_ag(:)    = (currentCohort%n*ccohort_hydr%errh2o_growturn_ag(:)    + &
                                              nextCohort%n*ncohort_hydr%errh2o_growturn_ag(:))/newn
     ccohort_hydr%errh2o_pheno_ag(:)       = (currentCohort%n*ccohort_hydr%errh2o_pheno_ag(:)       + &
                                              nextCohort%n*ncohort_hydr%errh2o_pheno_ag(:))/newn
     ccohort_hydr%errh2o_growturn_troot(:) = (currentCohort%n*ccohort_hydr%errh2o_growturn_troot(:) + &
                                              nextCohort%n*ncohort_hydr%errh2o_growturn_troot(:))/newn
     ccohort_hydr%errh2o_pheno_troot(:)    = (currentCohort%n*ccohort_hydr%errh2o_pheno_troot(:)    + &
                                              nextCohort%n*ncohort_hydr%errh2o_pheno_troot(:))/newn
     ccohort_hydr%errh2o_growturn_aroot(:) = (currentCohort%n*ccohort_hydr%errh2o_growturn_aroot(:) + &
                                              nextCohort%n*ncohort_hydr%errh2o_growturn_aroot(:))/newn
     ccohort_hydr%errh2o_pheno_aroot(:)    = (currentCohort%n*ccohort_hydr%errh2o_pheno_aroot(:)    + &
                                              nextCohort%n*ncohort_hydr%errh2o_pheno_aroot(:))/newn
     
     !ccohort_hydr%th_aroot_prev(:)
     !ccohort_hydr%th_aroot_prev_uncorr(:)

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
    call ccohort_hydr%AllocateHydrCohortArrays(currentSite%si_hydr%nlevsoi_hyd)

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

  subroutine InitHydrSites(sites,bc_in,numpft)

       ! Arguments
       type(ed_site_type),intent(inout),target :: sites(:)
       type(bc_in_type),intent(in)             :: bc_in(:)
       integer,intent(in)                      :: numpft

       ! Locals
       integer :: nsites
       integer :: s
       type(ed_site_hydr_type),pointer :: csite_hydr
       

       if ( hlm_use_planthydro.eq.ifalse ) return
       
       ! Initialize any derived hydraulics parameters
       call InitHydraulicsDerived(numpft)
       
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
          sites(s)%si_hydr%nlevsoi_hyd = bc_in(s)%nlevsoil
          call sites(s)%si_hydr%InitHydrSite()
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
     integer :: j
     integer :: nsites
     integer :: nlevsoil         ! Number of soil layers
     integer :: nlevsoil_hyd     ! Number of hydraulically relevant soil layers
       
     nsites = ubound(sites,1)

     do s = 1,nsites
        site_hydr => sites(s)%si_hydr
        
        nlevsoil     = bc_in(s)%nlevsoil
        nlevsoil_hyd = site_hydr%nlevsoi_hyd

        if ( nlevsoil_hyd == 1) then

           h2osoi_liqvol = min(bc_in(s)%eff_porosity_sl(nlevsoil), &
                bc_in(s)%h2o_liq_sisl(nlevsoil)/(bc_in(s)%dz_sisl(nlevsoil)*denh2o))

           site_hydr%h2osoi_liqvol_shell(nlevsoil_hyd,1:nshell) = h2osoi_liqvol
           site_hydr%h2osoi_liq_prev(nlevsoil_hyd)              = bc_in(s)%h2o_liq_sisl(nlevsoil)
        else
           do j = 1,nlevsoil_hyd
              h2osoi_liqvol = min(bc_in(s)%eff_porosity_sl(j), &
                    bc_in(s)%h2o_liq_sisl(j)/(bc_in(s)%dz_sisl(j)*denh2o))

              site_hydr%h2osoi_liqvol_shell(j,1:nshell) = h2osoi_liqvol
              site_hydr%h2osoi_liq_prev(j)              = bc_in(s)%h2o_liq_sisl(j)
           end do
        end if

        do j = 1, nlevsoil_hyd
           ! Calculate the matric potential on the innner shell (this is used to initialize
           ! xylem and root pressures in new cohorts)
           call swcCampbell_psi_from_th(site_hydr%h2osoi_liqvol_shell(j,1), &
                bc_in(s)%watsat_sisl(j), (-1.0_r8)*bc_in(s)%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                bc_in(s)%bsw_sisl(j), smp)

           site_hydr%psisoi_liq_innershell(j) = smp
           
        end do
        site_hydr%l_aroot_layer(1:site_hydr%nlevsoi_hyd) = 0.0_r8
        
     end do

     ! 
     !!     call UpdateH2OVeg(nsites,sites,bc_out)

     ! --------------------------------------------------------------------------------
     ! All other ed_Hydr_site_type variables are initialized elsewhere:
     !
     ! init_patch() -> updateSizeDepRhizHydProps -> shellgeom()
     !         this%v_shell
     !         this%v_shell_1D
     !         this%r_node_shell
     !         this%r_out_shell
     !         this%r_out_shell_1D
     !         this%r_node_shell_1D
     !
     ! init_patch() -> updateSizeDepRhizHydProps()
     !         this%l_aroot_layer_init
     !         this%l_aroot_1D
     !         this%kmax_upper_shell
     !         this%kmax_bound_shell
     !         this%kmax_lower_shell
     !         this%kmax_upper_shell_1D
     !         this%kmax_bound_shell_1D
     !         this%kmax_lower_shell_1D
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

     use EDTypesMod, only : AREA
     use EDTypesMod     , only : dump_cohort

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
                    sum(ccohort_hydr%th_troot(:)*ccohort_hydr%v_troot(:)) + &
                    sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                    denh2o*currentCohort%n
              endif

              currentCohort => currentCohort%shorter
           enddo !cohort
           currentPatch => currentPatch%younger
        enddo !end patch loop
        
        csite_hydr%h2oveg              = csite_hydr%h2oveg              / AREA

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
     ! This subroutine is called to caluate the water requirement for newly recruited cohorts
     ! The water update is allocated proportionally to the root biomass, which could be updated
     ! to accomodate the soil moisture and rooting depth for small seedlings (Chonggang XU).
     ! After the root water uptake, is_newly_recruited flag is set to false.
     ! ----------------------------------------------------------------------------------

     use EDTypesMod, only : AREA
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
     real(r8) :: roota !root distriubiton parameter a
     real(r8) :: rootb !root distriubiton parameter b
     real(r8) :: rootfr !fraction of root in different soil layer
     real(r8) :: recruitw !water for newly recruited cohorts (kg water/m2/s)
     real(r8) :: recruitw_total ! total water for newly recruited cohorts (kg water/m2/s)
     real(r8) :: err !mass error of water for newly recruited cohorts (kg water/m2/s)
     real(r8) :: sumrw_uptake !sum of water take for newly recruited cohorts (kg water/m2/s)
     
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
	        roota    =  EDPftvarcon_inst%roota_par(ft)
                rootb    =  EDPftvarcon_inst%rootb_par(ft)
	        recruitw =  (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:))    + &
                    sum(ccohort_hydr%th_troot(:)*ccohort_hydr%v_troot(:))             + &
                    sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                    denh2o*currentCohort%n/AREA/dtime
		recruitw_total = recruitw_total + recruitw
                if( csite_hydr%nlevsoi_hyd == 1) then
                    csite_hydr%recruit_w_uptake(1) = csite_hydr%recruit_w_uptake(1)+ &
		                                    recruitw
                else
                  do j=1,csite_hydr%nlevsoi_hyd
                    if(j == 1) then
                       rootfr = zeng2001_crootfr(roota, rootb, bc_in(s)%zi_sisl(j))
                    else
                       rootfr = zeng2001_crootfr(roota, rootb, bc_in(s)%zi_sisl(j)) - &
                          zeng2001_crootfr(roota, rootb, bc_in(s)%zi_sisl(j-1))
                    end if
		    csite_hydr%recruit_w_uptake(j) = csite_hydr%recruit_w_uptake(j) + &
		                                    recruitw*rootfr
                  end do
                end if
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
	   do j=1,csite_hydr%nlevsoi_hyd
	     csite_hydr%recruit_w_uptake(j) = csite_hydr%recruit_w_uptake(j) + &
	         err*csite_hydr%recruit_w_uptake(j)/sumrw_uptake
	   enddo	   
	endif
     end do ! site loop
	
  end subroutine RecruitWUptake	   
  
  !=====================================================================================
  subroutine ConstrainRecruitNumber(csite,ccohort, bc_in)

     ! ---------------------------------------------------------------------------
     ! This subroutine constrains the number of plants so that there is enought water 
     !  for newly recruited individuals from the soil
     ! ---------------------------------------------------------------------------
     use EDTypesMod, only : AREA

     ! Arguments
     type(ed_site_type), intent(inout), target     :: csite
     type(ed_cohort_type) , intent(inout), target  :: ccohort
     type(bc_in_type)    , intent(in)              :: bc_in 

     ! Locals
     type(ed_cohort_hydr_type), pointer :: ccohort_hydr
     type(ed_site_hydr_type), pointer :: csite_hydr
     real(r8) :: tmp1
     real(r8) :: watres_local   !minum water content
     real(r8) :: total_water !total water in rhizosphere at a specific layer (m^3)
     real(r8) :: total_water_min !total minimum water in rhizosphere at a specific layer (m^3)
     real(r8) :: roota !root distriubiton parameter a
     real(r8) :: rootb !root distriubiton parameter b
     real(r8) :: rootfr !fraction of root in different soil layer
     real(r8) :: recruitw !water for newly recruited cohorts (kg water/m2/individual)   
     real(r8) :: n, nmin !number of individuals in cohorts  
     integer :: s, j, ft

     roota                     =  EDPftvarcon_inst%roota_par(ccohort%pft)
     rootb                     =  EDPftvarcon_inst%rootb_par(ccohort%pft)
    
     csite_hydr => csite%si_hydr
     ccohort_hydr =>ccohort%co_hydr
     recruitw =  (sum(ccohort_hydr%th_ag(:)*ccohort_hydr%v_ag(:))    + &
                    sum(ccohort_hydr%th_troot(:)*ccohort_hydr%v_troot(:))  + &
                    sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                    denh2o
     
     do j=1,csite_hydr%nlevsoi_hyd
          if(j == 1) then
                       rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j))
          else
                       rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j)) - &
                      zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j-1))
         end if
	 cohort_recruit_water_layer(j) = recruitw*rootfr
     end do  
     
     do j=1,csite_hydr%nlevsoi_hyd
          select case (iswc)
           case (van_genuchten) 
                 write(fates_log(),*) &
                 'Van Genuchten plant hydraulics is inoperable until further notice'
                 call endrun(msg=errMsg(sourcefile, __LINE__)) 
           case (campbell)
                 call swcCampbell_satfrac_from_psi(bc_in%smpmin_si*denh2o*grav*1.e-9_r8, &
                                    (-1._r8)*bc_in%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                                    bc_in%bsw_sisl(j),     &
                                    tmp1)
                 call swcCampbell_th_from_satfrac(tmp1, &
                               bc_in%watsat_sisl(j),   &
                                    watres_local)
                 
	    
           case default
         end select
         total_water = sum(csite_hydr%v_shell(j,:)*csite_hydr%h2osoi_liqvol_shell(j,:)) * &
                         csite_hydr%l_aroot_layer(j)/&
                         bc_in %dz_sisl(j) 
	 total_water_min = sum(csite_hydr%v_shell(j,:)*watres_local) * &
                         csite_hydr%l_aroot_layer(j)/&
                         bc_in %dz_sisl(j)  		  
	 !assumes that only 50% is available for recruit water....
	 recruit_water_avail_layer(j)=0.5_r8*max(0.0_r8,total_water-total_water_min)
	  
     end do
     
     nmin  = 1.0e+36 
     do j=1,csite_hydr%nlevsoi_hyd
       if(cohort_recruit_water_layer(j)>0.0_r8) then
          n = recruit_water_avail_layer(j)/cohort_recruit_water_layer(j)
          nmin = min(n, nmin) 
       endif
     end do
     ccohort%n = min (ccohort%n, nmin) 

  end subroutine ConstrainRecruitNumber


  ! =====================================================================================

  subroutine SavePreviousRhizVolumes(currentSite, bc_in)

    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in) :: bc_in
    type(ed_site_hydr_type), pointer    :: csite_hydr
    integer                             :: nlevsoi_hyd
    
    csite_hydr => currentSite%si_hydr
    nlevsoi_hyd = csite_hydr%nlevsoi_hyd
     
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
    use EDTypesMod           , only : AREA
    
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
    real(r8)                       :: large_kmax_bound = 1.e4_r8   ! for replacing kmax_bound_shell wherever the 
                                                                   ! innermost shell radius is less than the assumed 
                                                                   ! absorbing root radius rs1
    real(r8)                       :: kmax_root_surf               ! maximum conducitivity for unit root surface
                                                                   ! (kg water/m2 root area/Mpa/s)
                                                                   ! 1.e-5_r8 from Rudinger et al 1994             	
    real(r8)                       :: kmax_root_surf_total         !maximum conducitivity for total root surface(kg water/Mpa/s)
    real(r8)                       :: kmax_soil_total              !maximum conducitivity for total root surface(kg water/Mpa/s)
    integer                        :: nlevsoi_hyd	  

    !-----------------------------------------------------------------------
 
    csite_hydr => currentSite%si_hydr
    nlevsoi_hyd = csite_hydr%nlevsoi_hyd

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
    kmax_root_surf = hydr_kmax_rsurf
    csite_hydr%l_aroot_1D = sum( csite_hydr%l_aroot_layer(:))
    
    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    do j = 1,nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          call shellGeom( csite_hydr%l_aroot_layer(j), csite_hydr%rs1(j), AREA, bc_in%dz_sisl(j), &
                csite_hydr%r_out_shell(j,:), csite_hydr%r_node_shell(j,:),csite_hydr%v_shell(j,:))
       end if !has l_aroot_layer changed?
    enddo
    call shellGeom( csite_hydr%l_aroot_1D, csite_hydr%rs1(1), AREA, sum(bc_in%dz_sisl(1:nlevsoi_hyd)), &
          csite_hydr%r_out_shell_1D(:), csite_hydr%r_node_shell_1D(:), csite_hydr%v_shell_1D(:))
    
    do j = 1,csite_hydr%nlevsoi_hyd
       
       hksat_s = bc_in%hksat_sisl(j) * 1.e-3_r8 * 1/grav * 1.e6_r8

       ! proceed only if the total absorbing root length (site-level) has changed in this layer
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then

          do k = 1,nshell
	     if(k == 1) then
	        kmax_root_surf_total = kmax_root_surf*2._r8*pi_const *csite_hydr%rs1(j)* &
		                       csite_hydr%l_aroot_layer(j)
                if(csite_hydr%r_node_shell(j,k) <= csite_hydr%rs1(j)) then
		   !csite_hydr%kmax_upper_shell(j,k)  = large_kmax_bound
                   !csite_hydr%kmax_bound_shell(j,k)  = large_kmax_bound
                   !csite_hydr%kmax_lower_shell(j,k)  = large_kmax_bound
                   csite_hydr%kmax_upper_shell(j,k)  = kmax_root_surf_total
                   csite_hydr%kmax_bound_shell(j,k)  = kmax_root_surf_total
                   csite_hydr%kmax_lower_shell(j,k)  = kmax_root_surf_total

                else

		   kmax_soil_total = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                         log(csite_hydr%r_node_shell(j,k)/csite_hydr%rs1(j))*hksat_s

                   !csite_hydr%kmax_upper_shell(j,k)  = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
		   !      log(csite_hydr%r_node_shell(j,k)/csite_hydr%rs1(j))*hksat_s
                   !csite_hydr%kmax_bound_shell(j,k)  = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
		   !      log(csite_hydr%r_node_shell(j,k)/csite_hydr%rs1(j))*hksat_s
                   !csite_hydr%kmax_lower_shell(j,k)  = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
		   !      log(csite_hydr%r_node_shell(j,k)/csite_hydr%rs1(j))*hksat_s

                   csite_hydr%kmax_upper_shell(j,k)  = (1._r8/kmax_root_surf_total + &
		                       1._r8/kmax_soil_total)**(-1._r8)    
                   csite_hydr%kmax_bound_shell(j,k)  = (1._r8/kmax_root_surf_total + &
		                       1._r8/kmax_soil_total)**(-1._r8) 
                   csite_hydr%kmax_lower_shell(j,k)  = (1._r8/kmax_root_surf_total + &
		                       1._r8/kmax_soil_total)**(-1._r8)
                end if
		if(j == 1) then
                   if(csite_hydr%r_node_shell(j,k) <= csite_hydr%rs1(j)) then
                     csite_hydr%kmax_upper_shell_1D(k)  = csite_hydr%kmax_upper_shell(1,k)
                     csite_hydr%kmax_bound_shell_1D(k)  = csite_hydr%kmax_bound_shell(1,k)
                     csite_hydr%kmax_lower_shell_1D(k)  = csite_hydr%kmax_lower_shell(1,k)
                   else
                      csite_hydr%kmax_upper_shell_1D(k) = csite_hydr%kmax_upper_shell(1,k)
                      csite_hydr%kmax_bound_shell_1D(k) = csite_hydr%kmax_bound_shell(1,k)
                      csite_hydr%kmax_lower_shell_1D(k) = csite_hydr%kmax_lower_shell(1,k)
                   end if
                end if
             else
                csite_hydr%kmax_upper_shell(j,k)        = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                      log(csite_hydr%r_node_shell(j,k)/csite_hydr%r_out_shell(j,k-1))*hksat_s
                csite_hydr%kmax_bound_shell(j,k)        = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                      log(csite_hydr%r_node_shell(j,k)/csite_hydr%r_node_shell(j,k-1))*hksat_s
                csite_hydr%kmax_lower_shell(j,k)        = 2._r8*pi_const*csite_hydr%l_aroot_layer(j) / &
                      log(csite_hydr%r_out_shell(j,k)/csite_hydr%r_node_shell(j,k  ))*hksat_s
                if(j == 1) then
                   csite_hydr%kmax_upper_shell_1D(k)    = 2._r8*pi_const*csite_hydr%l_aroot_1D / &
                         log(csite_hydr%r_node_shell_1D(k)/csite_hydr%r_out_shell_1D(k-1))*hksat_s
                   csite_hydr%kmax_bound_shell_1D(k)    = 2._r8*pi_const*csite_hydr%l_aroot_1D / &
                         log(csite_hydr%r_node_shell_1D(k)/csite_hydr%r_node_shell_1D(k-1))*hksat_s
                   csite_hydr%kmax_lower_shell_1D(k)    = 2._r8*pi_const*csite_hydr%l_aroot_1D / &
                         log(csite_hydr%r_out_shell_1D( k)/csite_hydr%r_node_shell_1D(k  ))*hksat_s
                end if
             end if
          enddo ! loop over rhizosphere shells
       end if !has l_aroot_layer changed?
    enddo ! loop over soil layers



    return
  end subroutine UpdateSizeDepRhizVolLenCon


  ! =====================================================================================


  subroutine updateSizeDepRhizHydProps(currentSite, bc_in )
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:

    use EDTypesMod           , only : AREA
    
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in) :: bc_in


    ! Save current volumes, lenghts and nodes to an "initial"
    ! used to calculate effects in states later on.
    
    call SavePreviousRhizVolumes(currentSite, bc_in)

    ! Update the properties of the vegetation-soil hydraulic environment
    ! these are independent on the water state
    
    call UpdateSizeDepRhizVolLenCon(currentSite, bc_in)


    return
  end subroutine updateSizeDepRhizHydProps
  
  ! =================================================================================

  subroutine updateSizeDepRhizHydStates(currentSite, bc_in)
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
    real(r8) :: w_shell_new                              ! updated water mass in rhizosphere compartment               [kg]
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
    ! calculate initial s, psi by layer and shell
    do j = 1, csite_hydr%nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          select case (iswc)
          case (van_genuchten)
             do k = 1,nshell
                s_shell_init(j,k)     = (csite_hydr%h2osoi_liqvol_shell(j,k) - bc_in%watres_sisl(j)) / &
                      (bc_in%watsat_sisl(j) - bc_in%watres_sisl(j))
                write(fates_log(),*) 'VG is not available yet'
                call endrun(msg=errMsg(sourcefile, __LINE__))
!                call swcVG_psi_from_satfrac(s_shell_init(j,k),alpha_VG(c,j),n_VG(c,j),m_VG(c,j),l_VG(c,j),psi_shell_init(j,k))
             end do
          case (campbell)
             do k = 1,nshell
                s_shell_init(j,k)     = (csite_hydr%h2osoi_liqvol_shell(j,k) - bc_in%watres_sisl(j)) / &
                      (bc_in%watsat_sisl(j) - bc_in%watres_sisl(j))
                call swcCampbell_psi_from_satfrac( s_shell_init(j,k), &
                      bc_in%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                      bc_in%bsw_sisl(j),psi_shell_init(j,k))
             end do
          case default
             write(fates_log(),*) 'Somehow you picked a PT function that DNE'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end select
       end if !has l_aroot_coh changed?
    enddo
    
    ! interpolate initial psi values by layer and shell
    ! BOC...To-Do: need to constrain psi to be within realistic limits (i.e., < 0)
    do j = 1,csite_hydr%nlevsoi_hyd
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
    do j = 1,csite_hydr%nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          select case (iswc)
          case (van_genuchten)
             do k = 1,nshell
                write(fates_log(),*) 'VG is not available yet'
                call endrun(msg=errMsg(sourcefile, __LINE__))
                !                call swcVG_satfrac_from_psi(psi_shell_interp(j,k), &
                !                alpha_VG(c,j),n_VG(c,j),m_VG(c,j),l_VG(c,j),s_shell_interp(j,k))
             enddo
          case (campbell)
             do k = 1,nshell
                call swcCampbell_satfrac_from_psi(psi_shell_interp(j,k), &
                      (-1._r8)*bc_in%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                      bc_in%bsw_sisl(j), &
                      s_shell_interp(j,k))
             end do
          case default
             write(fates_log(),*) 'Somehow you picked a PT function that DNE'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end select

       end if !has l_aroot_coh changed?
    enddo
    
    ! accumlate water across shells for each layer (initial and interpolated)
    do j = 1,csite_hydr%nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          w_layer_init(j)      = 0._r8
          w_layer_interp(j)    = 0._r8
          v_rhiz(j)            = 0._r8
          do k = 1,nshell
             w_layer_init(j)   = w_layer_init(j) + denh2o/bc_in%dz_sisl(j) * &
                   ( csite_hydr%l_aroot_layer_init(j) * &
                   csite_hydr%v_shell_init(j,k)*csite_hydr%h2osoi_liqvol_shell(j,k) )
             w_layer_interp(j) = w_layer_interp(j) + denh2o/bc_in%dz_sisl(j) * &
                   ( csite_hydr%l_aroot_layer(j)*csite_hydr%v_shell(j,k) * &
                   (s_shell_interp(j,k)*(bc_in%watsat_sisl(j)-bc_in%watres_sisl(j))+bc_in%watres_sisl(j)) )
             v_rhiz(j)         = v_rhiz(j) + csite_hydr%v_shell(j,k)
          enddo
       end if !has l_aroot_coh changed?
    enddo
    
    ! estimate delta_s across all shells needed to ensure total water in each layer doesn't change
    ! BOC...FIX: need to handle special cases where delta_s causes s_shell to go above or below 1 or 0, respectively.
    do j = 1,csite_hydr%nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          delta_s(j) = (( w_layer_init(j) - w_layer_interp(j) )/( v_rhiz(j) * &
                denh2o*csite_hydr%l_aroot_layer(j)/bc_in%dz_sisl(j) ) - bc_in%watres_sisl(j)) / &
                (bc_in%watsat_sisl(j)-bc_in%watres_sisl(j))
       end if !has l_aroot_coh changed?
    enddo
    
    ! update h2osoi_liqvol_shell and h2osoi_liq_shell
    do j = 1,csite_hydr%nlevsoi_hyd
       ! proceed only if l_aroot_coh has changed
       if( csite_hydr%l_aroot_layer(j) /= csite_hydr%l_aroot_layer_init(j) ) then
          w_layer_new(j)                = 0._r8
          do k = 1,nshell
             s_shell_interp(j,k)        = s_shell_interp(j,k) + delta_s(j)
   	     csite_hydr%h2osoi_liqvol_shell(j,k) = s_shell_interp(j,k) * &
               ( bc_in%watsat_sisl(j)-bc_in%watres_sisl(j) ) + bc_in%watres_sisl(j)
	     w_shell_new                = csite_hydr%h2osoi_liqvol_shell(j,k) * &
            csite_hydr%v_shell(j,k) * denh2o
             w_layer_new(j)             = w_layer_new(j) + w_shell_new
          enddo
          h2osoi_liq_col_new(j)         = w_layer_new(j)/( v_rhiz(j)/bc_in%dz_sisl(j) )
       end if !has l_aroot_coh changed?
    enddo
    
    ! balance check
    do j = 1,csite_hydr%nlevsoi_hyd
       ! BOC: PLEASE CHECK UNITS ON h2o_liq_sisl(j) (RGK)
       errh2o(j) = h2osoi_liq_col_new(j) - bc_in%h2o_liq_sisl(j)
       if (abs(errh2o(j)) > 1.e-4_r8) then
          found = .true.
          indexj = j
          if( found ) then
             write(fates_log(),*)'WARNING:  water balance error ',&
                   ' local indexj= ',indexj,&
                   ' errh2o= ',errh2o(indexj)
          end if
       end if
    enddo
   end if !nshell > 1
    
  end subroutine updateSizeDepRhizHydStates


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
                   ccohort%co_hydr%btran(1) * &
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
    use EDtypesMod         , only : AREA
    !
    ! !ARGUMENTS:
    integer, intent(in)                       :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type), intent(in)              :: bc_in(nsites)
    type(bc_out_type), intent(inout)          :: bc_out(nsites)

    ! Locals
    real(r8) :: dwat_kgm2                                            ! change in layer water content              [kg/m2]
    type(ed_site_hydr_type), pointer :: csite_hydr
    integer  :: s,j,k                                  ! site, soil layer, rhizosphere shell indicies
    integer  :: i,f,ff,kk                              ! indicies
    integer  :: indexj                          ! column and layer indices where there is a water balance error
    integer  :: ordered(nshell) = (/(i,i=1,nshell,1)/) ! array of rhizosphere indices which have been ordered
    real(r8) :: area_col                               ! column area                                                    [m2]
    real(r8) :: v_cum                                  ! cumulative shell volume from driest/wettest shell to kth shell [m3]
    real(r8) :: dwat_kg                                ! water remaining to be distributed across shells                [kg]
    real(r8) :: thdiff                                 ! water content difference between ordered adjacent rhiz shells  [m3 m-3]
    real(r8) :: wdiff                                  ! mass of water represented by thdiff over previous k shells     [kg]
    real(r8) :: errh2o(nlevsoi_hyd_max)                ! water budget error after updating                              [kg/m2]
    real(r8) :: cumShellH2O                            ! sum of water in all the shells of a specific layer             [kg/m2]
    real(r8) :: h2osoi_liq_shell(nlevsoi_hyd_max,nshell) !water in the rhizosphere shells                               [kg]
    integer  :: tmp                                    ! temporary
    logical  :: found                                  ! flag in search loop
    !-----------------------------------------------------------------------

    do s = 1,nsites


       ! First step, identify how the liquid water in each layer has changed
       ! since the last time it was updated. This should be due to drainage.
       ! The drainage component should be the total change in liquid water content from the last time
       ! the hydraulics driver was called, and then adding back in the losses due to root uptake
       ! (which was already taken out).

       ! BOC: This was previously in HydrologyDrainage:

       csite_hydr => sites(s)%si_hydr

       do j = 1,csite_hydr%nlevsoi_hyd 
          cumShellH2O=sum(csite_hydr%h2osoi_liqvol_shell(j,:) *csite_hydr%v_shell(j,:)) &
               / bc_in(s)%dz_sisl(j) * csite_hydr%l_aroot_layer(j) * denh2o/AREA
          
          if(csite_hydr%nlevsoi_hyd == 1) then
             dwat_kgm2 = bc_in(s)%h2o_liq_sisl(bc_in(s)%nlevsoil) - cumShellH2O
          else    !  if(csite_hydr%nlevsoi_hyd == bc_in(s)%nlevsoil ) then
             dwat_kgm2 = bc_in(s)%h2o_liq_sisl(j) - cumShellH2O
          end if

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
             v_cum  = sum(csite_hydr%v_shell(j,ordered(1:k))) / &
                      bc_in(s)%dz_sisl(j) * csite_hydr%l_aroot_layer(j)
             wdiff  = thdiff * v_cum * denh2o
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
             v_cum  = sum(csite_hydr%v_shell(j,ordered(1:nshell))) / bc_in(s)%dz_sisl(j) * &
                      csite_hydr%l_aroot_layer(j)
             thdiff = dwat_kg / v_cum / denh2o
             do k = nshell, 1, -1
                csite_hydr%h2osoi_liqvol_shell(j,k) = csite_hydr%h2osoi_liqvol_shell(j,k) + thdiff
             end do
          end if
          
          ! m3/m3 * Total volume m3 * kg/m3 = kg
          h2osoi_liq_shell(j,:) = csite_hydr%h2osoi_liqvol_shell(j,:) * &
               csite_hydr%v_shell(j,:) / bc_in(s)%dz_sisl(j) * csite_hydr%l_aroot_layer(j) * denh2o
          
       enddo
       
       ! balance check
       if(csite_hydr%nlevsoi_hyd .ne. 1) then
          do j = 1,csite_hydr%nlevsoi_hyd
             errh2o(j) = sum(h2osoi_liq_shell(j,:))/AREA - bc_in(s)%h2o_liq_sisl(j)
             
             if (abs(errh2o(j)) > 1.e-9_r8) then
                found = .true.
                indexj = j
                if( found ) then
                   write(fates_log(),*)'WARNING:  water balance error in FillDrainRhizShells',&
                        ' local indexj= ',indexj,&
                        ' errh2o= ',errh2o(indexj)
                end if
             end if
          enddo
       else
          errh2o(csite_hydr%nlevsoi_hyd) = sum(h2osoi_liq_shell(csite_hydr%nlevsoi_hyd,:))/AREA - sum( bc_in(s)%h2o_liq_sisl(:) )
       end if
       
    end do
    return
   end subroutine FillDrainRhizShells

   ! ====================================================================================

  subroutine hydraulics_bc ( nsites, sites,bc_in,bc_out,dtime )
     
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
     use EDTypesMod        , only : AREA
    use FatesUtilsMod  , only : check_var_real
    use EDTypesMod     , only : dump_cohort

     ! ARGUMENTS:
     ! -----------------------------------------------------------------------------------
     integer,intent(in)                      :: nsites
     type(ed_site_type),intent(inout),target :: sites(nsites)
     type(bc_in_type),intent(in)             :: bc_in(nsites)
     type(bc_out_type),intent(inout)         :: bc_out(nsites)
     real(r8),intent(in)                     :: dtime
     
     !
     ! !LOCAL VARIABLES:
     character(len=*), parameter :: sub = 'clm::Hydraulics_bc'
     integer :: iv  ! leaf layer
     integer :: ifp ! index of FATES patch
     integer :: s   ! index of FATES site
     integer :: j,jj! soil layer
     integer :: k   ! 1D plant-soil continuum array
     integer :: ft  ! plant functional type index
     integer :: t   ! previous timesteps (for lwp stability calculation)
     integer :: nstep !number of time steps

     !----------------------------------------------------------------------
     
     type (ed_patch_type),  pointer :: cpatch
     type (ed_cohort_type), pointer :: ccohort
     
     ! hydraulics global constants
     real(r8), parameter :: thresh          = 1.e-7_r8  ! threshold for water balance error (warning only) [mm h2o]
     real(r8), parameter :: thresh_break    = 1.e-4_r8  ! threshold for water balance error (stop model)   [mm h2o]
     real(r8), parameter :: small_theta_num = 1.e-7_r8  ! avoids theta values equalling thr or ths         [m3 m-3]
     
     ! hydraulics timestep adjustments for acceptable water balance error
     integer  :: maxiter        = 1            ! maximum iterations for timestep reduction                       [-]
     integer  :: imult          = 3            ! iteration index multiplier                                      [-]
     real(r8) :: we_area_outer                 ! 1D plant-soil continuum water error                             [kgh2o m-2 individual-1]
     
     ! cohort-specific arrays to hold 1D hydraulics geometric & state variables for entire continuum (leaf,stem,root,soil)
     real(r8) :: z_node(       n_hypool_tot)      ! nodal height of water storage compartments                      [m]
     real(r8) :: z_node_1l(    n_hypool_tot)      ! nodal height of water storage compartments (single-layer soln)  [m]
     real(r8) :: v_node(       n_hypool_tot)      ! volume of water storage compartments                            [m3]
     real(r8) :: v_node_1l(    n_hypool_tot)      ! volume of water storage compartments (single-layer soln)        [m3]
     real(r8) :: psi_node(     n_hypool_tot)      ! water potential in water storage compartments                   [MPa]
     real(r8) :: psi_node_1l(  n_hypool_tot)      ! water potential in water storage compartments (single-layer soln) [MPa]
     real(r8) :: flc_node_1l(  n_hypool_tot)      ! fractional loss of conductivity (single-layer soln)             [-]
     real(r8) :: flc_min_node( n_hypool_tot-nshell)! minimum attained fractional loss of conductivity (for xylem refilling dynamics) [-]
     real(r8) :: dflcdpsi_node_1l(n_hypool_tot)   ! derivative of flc_node_1l wrt psi                               [MPa-1]
     real(r8) :: ths_node(     n_hypool_tot)      ! saturated volumetric water in water storage compartments        [m3 m-3]
     real(r8) :: ths_node_1l(  n_hypool_tot)      ! saturated volumetric water in water storage compartments (single-layer soln) [m3 m-3]
     real(r8) :: thr_node(     n_hypool_tot)      ! residual volumetric water in water storage compartments         [m3 m-3]
     real(r8) :: thr_node_1l(  n_hypool_tot)      ! residual volumetric water in water storage compartments (single-layer soln) [m3 m-3]
     real(r8) :: the_node(     n_hypool_tot)      ! error resulting from supersaturation or below-residual th_node  [m3 m-3]
     real(r8) :: the_node_1l(  n_hypool_tot)      ! like the_node(:) but for specific single soil layer             [m3 m-3]
     real(r8) :: th_node(      n_hypool_tot)      ! volumetric water in water storage compartments                  [m3 m-3]
     real(r8) :: th_node_1l(   n_hypool_tot)      ! volumetric water in water storage compartments (single-layer soln) [m3 m-3]
     real(r8) :: dth_node(     n_hypool_tot)      ! change in volumetric water in water storage compartments        [m3 m-3]
     real(r8) :: dth_node_1l(  n_hypool_tot)      ! like dth_node_1l(:) but for specific single soil layer          [m3 m-3]
     real(r8) :: kmax_bound(   n_hypool_tot)      ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
     real(r8) :: kmax_bound_1l(n_hypool_tot)      ! lower boundary maximum hydraulic conductance of compartments (single-layer soln) [kg s-1 MPa-1]
     real(r8) :: kmax_upper(   n_hypool_tot)      ! maximum hydraulic conductance from node to upper boundary       [kg s-1 MPa-1]
     real(r8) :: kmax_upper_1l(n_hypool_tot)      ! maximum hydraulic conductance from node to upper boundary (single-layer soln)    [kg s-1 MPa-1]
     real(r8) :: kmax_lower(   n_hypool_tot)      ! maximum hydraulic conductance from node to lower boundary       [kg s-1 MPa-1]
     real(r8) :: kmax_lower_1l(n_hypool_tot)      ! maximum hydraulic conductance from node to lower boundary (single-layer soln)    [kg s-1 MPa-1]
     real(r8) :: hdiff_bound_1l( nshell+1)     !
     real(r8) :: k_bound_1l(     nshell+1)     !
     real(r8) :: dhdiffdpsi0_1l( nshell+1)     !
     real(r8) :: dhdiffdpsi1_1l( nshell+1)     !
     real(r8) :: dkbounddpsi0_1l(nshell+1)     !
     real(r8) :: dkbounddpsi1_1l(nshell+1)     !
     real(r8) :: l_aroot_tot_coh               ! total length of absorbing roots across all soil layers (cohort) [m]
     real(r8) :: dwat_veg_coh                  ! total indiv change in stored vegetation water over a timestep   [kg]
     
     ! column-specific arrays to hold rhizosphere geometric & state variables
     real(r8) :: h2osoi_liqvol
     real(r8) :: dz_tot                        ! total soil depth (to bottom of bottom layer)                    [m]
     real(r8) :: l_aroot_tot_col               ! total length of absorbing roots across all soil layers          [m]
     real(r8) :: dth_layershell_col(nlevsoi_hyd_max,nshell) ! accumulated water content change over all cohorts in a column   [m3 m-3]
     real(r8) :: ths_shell_1D(nshell)          ! saturated water content of rhizosphere compartment              [m3 m-3]
     real(r8) :: thr_shell_1D(nshell)          ! residual water content of rhizosphere compartment               [m3 m-3]
     real(r8) :: kmax_bound_shell_1l(nshell)   ! like kmax_bound_shell_1D(:) but for specific single soil layer  [kg s-1 MPa-1]
     real(r8) :: psi_node_shell_1D(nshell)     ! soil matric potential of rhizosphere compartment                [MPa]
     real(r8) :: ths_aroot_1D                  ! saturated water content of 1D representation of fine roots      [m3 m-3]
     real(r8) :: thr_aroot_1D                  ! residual water content of 1D representation of fine roots       [m3 m-3]
     real(r8) :: vtot_aroot_1D                 ! sum of fine root volume across soil layers                      [m3]
     real(r8) :: psi_node_aroot_1D             ! water potential of absorbing root                               [MPa]
     
     ! hydraulics conductances
     real(r8) :: kmax_bound_bylayershell(nlevsoi_hyd_max,nshell) ! maximum conductance at shell boundaries in each layer [kg s-1 MPa-1]
     real(r8) :: kmax_bound_aroot_soil1        ! maximum radial conductance of absorbing roots                   [kg s-1 MPa-1]
     real(r8) :: kmax_bound_aroot_soil2        ! maximum conductance to root surface from innermost rhiz shell   [kg s-1 MPa-1]
     real(r8) :: ksoil_bylayer(nlevsoi_hyd_max)        ! total rhizosphere conductance (over all shells) by soil layer   [MPa]
     real(r8) :: ksoil_tot                     ! total rhizosphere conductance (over all shells and soil layers  [MPa]
     real(r8) :: kbg_layer(nlevsoi_hyd_max)     ! total absorbing root & rhizosphere conductance (over all shells) by soil layer   [MPa]
     real(r8) :: kbg_tot                       ! total absorbing root & rhizosphere conductance (over all shells and soil layers  [MPa]
     real(r8) :: kmax_stem                     ! maximum whole-stem (above troot to leaf) conductance            [kg s-1 MPa-1]
     
     ! hydraulics other
     integer  :: ordered(nlevsoi_hyd_max) = (/(j,j=1,nlevsoi_hyd_max,1)/) ! array of soil layer indices which have been ordered
     real(r8) :: qflx_tran_veg_indiv           ! individiual transpiration rate                                  [kgh2o indiv-1 s-1]
     real(r8) :: qflx_tran_veg_patch_coh
     real(r8) :: gscan_patch                   ! sum of ccohort%gscan across all cohorts within a patch          
     real(r8) :: qtop_dt
     real(r8) :: dqtopdth_dthdt
     real(r8) :: sapflow
     real(r8) :: rootuptake
     real(r8) :: totalrootuptake               !total root uptake per unit area (kg h2o m-2 time step -1) 
     real(r8) :: totaldqtopdth_dthdt
     real(r8) :: totalqtop_dt                  !total transpriation per unit area (kg h2o m-2 time step -1) 
     real(r8) :: total_e                       !mass balance error (kg h2o m-2 time step -1)     
     integer  :: ncoh_col                      ! number of cohorts across all non-veg patches within a column
     real(r8) :: transp_col                    ! Column mean transpiration rate [mm H2O/m2]
                                               ! as defined by the input boundary condition
     real(r8) :: transp_col_check              ! Column mean transpiration rate [mm H2O/m2] as defined
                                               ! by the sum of water fluxes through the cohorts

     real(r8) :: patch_wgt                     ! fraction of current patch relative to the whole site
                                               ! note that this is almost but not quite cpatch%area/AREA
                                               ! as it regards the fraction of canopy area as the relevant
                                               ! area, and assumes that the HLM has it's own patch
                                               ! that is not tracked by FATES which accounts for all
                                               ! non-canopy areas across all patches					       

     real(r8) :: smp                           ! temporary for matric potential (MPa)
     integer  :: tmp
     real(r8) :: tmp1
     real(r8) :: watres_local
     integer  :: pick_1l(nshell+1) = (/(k,k=n_hypool_ag+n_hypool_troot+1,n_hypool_tot,1)/)
     real(r8) :: lwpdiff1, lwpdiff2, Rndiff1, Rndiff2, btran_prev
     logical  :: mono_decr_Rn                  ! flag indicating whether net Radiation is monotonically decreasing
     real(r8) :: refill_rate                   ! rate of xylem refilling  [fraction per unit time; s-1]
     real(r8) :: roota, rootb                  ! parameters for root distribution                                      [m-1]
     real(r8) :: rootfr                        ! root fraction at different soil layers
     real(r8) :: prev_h2oveg                   ! previous time step plant water storage (kg/m2)
     logical  :: recruitflag                   ! flag to check if there is newly recruited cohorts

     type(ed_site_hydr_type), pointer :: site_hydr
     type(ed_cohort_hydr_type), pointer :: ccohort_hydr
     integer  :: err_code = 0

     ! ----------------------------------------------------------------------------------
     ! Important note: We are interested in calculating the total fluxes in and out of the
     ! site/column.  Usually, when we do things like this, we acknowledge that FATES
     ! does not consider the bare ground patch.  However, since this routine
     ! calculates "column level" fluxes, we have to factor in that patch-level fluxes
     ! are only accounting for a portion of the area.
     ! ----------------------------------------------------------------------------------

     ! DEPRECATED: waterstate_inst%psisoi_liq_shell 
     ! Input:  [real(r8) (:,:,:)] soil matric potential (MPa) by layer and rhizosphere shell
     
     !for debug only
     !nstep = get_nstep()
     
     !For newly recruited cohorts, add the water uptake demand to csite_hydr%recruit_w_uptake
     call RecruitWUptake(nsites,sites,bc_in,dtime,recruitflag)
     
     !update water storage in veg after incorporating newly recuited cohorts
     if(recruitflag) call UpdateH2OVeg(nsites,sites,bc_out)
	     
     do s = 1, nsites
          
        site_hydr => sites(s)%si_hydr

        ! AVERAGE ROOT WATER UPTAKE (BY RHIZOSPHERE SHELL) ACROSS ALL COHORTS WITHIN A COLUMN
        dth_layershell_col(:,:) = 0._r8
        site_hydr%dwat_veg       = 0._r8
        site_hydr%errh2o_hyd     = 0._r8
	prev_h2oveg    = site_hydr%h2oveg
        ncoh_col       = 0

        ! Calculate the mean site level transpiration flux
        ! This is usefull to check on mass conservation
        ! of cohort level fluxes
        ! -------------------------------------------------
        ifp = 0
        cpatch => sites(s)%oldest_patch
        transp_col = 0.0_r8
        do while (associated(cpatch))
           ifp = ifp + 1
           patch_wgt = min(1.0_r8,cpatch%total_canopy_area/cpatch%area) * (cpatch%area/AREA)
           transp_col = transp_col +  bc_in(s)%qflx_transp_pa(ifp)*patch_wgt
           cpatch => cpatch%younger
        end do


        ifp = 0
        cpatch => sites(s)%oldest_patch
        do while (associated(cpatch))
           ifp = ifp + 1
           
           ! -----------------------------------------------------------------------------
           ! We apparently want to know the area fraction of
           ! contribution of this patch to the site/column
           ! In the interface:
           ! wt_ed(p) = this%fates(nc)%bc_out(s)%canopy_fraction_pa(ifp)
           ! and then in patch%wtcol(p) = patch%wt_ed(p)
           !
           ! From EDCanopyStructureMode.F90:update_hlm_dynamics():
           ! bc_out(s)%canopy_fraction_pa(ifp) = min(1.0_r8,currentPatch%total_canopy_area/currentPatch%area) * &
           !      (currentPatch%area/AREA)
           ! ----------------------------------------------------------------------------

           patch_wgt = min(1.0_r8,cpatch%total_canopy_area/cpatch%area) * (cpatch%area/AREA)

           ! Total volume transpired from this patch [mm H2O / m2 /s ] * [m2/patch] = [mm H2O / patch / s]
!           qflx_trans_patch_vol = bc_in(s)%qflx_transp_pa(ifp) * (patch_wgt * AREA)

!           do t=2, numLWPmem
!              cpatch_hydr%netRad_mem(t-1) = cpatch_hydr%netRad_mem(t)
!           end do
!           cpatch_hydr%netRad_mem(numLWPmem) = bc_in(s)%swrad_net_pa(ifp) - bc_in(s)%lwrad_net_pa(ifp)

           gscan_patch   = 0.0_r8
           ccohort=>cpatch%tallest
           do while(associated(ccohort))
              ccohort_hydr => ccohort%co_hydr
              gscan_patch       = gscan_patch + ccohort%g_sb_laweight
              if (gscan_patch < 0._r8) then
                 write(fates_log(),*) 'ERROR: negative gscan_patch!'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
              ccohort => ccohort%shorter
           enddo !cohort
           
           ccohort=>cpatch%tallest
           do while(associated(ccohort))
              ccohort_hydr => ccohort%co_hydr
              ft       = ccohort%pft
              ncoh_col = ncoh_col + 1
              ccohort_hydr%qtop_dt         = 0._r8
              ccohort_hydr%dqtopdth_dthdt  = 0._r8
              ccohort_hydr%sapflow         = 0._r8
              ccohort_hydr%rootuptake      = 0._r8
              
              ! Relative transpiration of this cohort from the whole patch
!!              qflx_rel_tran_coh = ccohort%g_sb_laweight/gscan_patch

              qflx_tran_veg_patch_coh      = bc_in(s)%qflx_transp_pa(ifp) * ccohort%g_sb_laweight/gscan_patch

              qflx_tran_veg_indiv          = qflx_tran_veg_patch_coh * cpatch%area* &
	                                     min(1.0_r8,cpatch%total_canopy_area/cpatch%area)/ccohort%n !AREA / ccohort%n
              
              ! [mm H2O/cohort/s] = [mm H2O / patch / s] / [cohort/patch]
!!              qflx_tran_veg_patch_coh      = qflx_trans_patch_vol * qflx_rel_tran_coh


		   
              if(site_hydr%nlevsoi_hyd > 1) then
                 ! BUCKET APPROXIMATION OF THE SOIL-ROOT HYDRAULIC GRADIENT (weighted average across layers)
                 !call map2d_to_1d_shells(soilstate_inst, waterstate_inst, g, c, rs1(c,1), ccohort_hydr%l_aroot_layer*ccohort%n, &
                 !                        (2._r8 * ccohort_hydr%kmax_treebg_layer(:)), ths_shell_1D, &
                 !                        thr_shell_1D, psi_node_shell_1D, &
                 !                        r_out_shell_1D, r_node_shell_1D, v_shell_1D, dz_tot, &
                 !                        ksoil_bylayer, ksoil_tot, kmax_bound_bylayershell)
                 !psi_node(  (n_hypool_tot-nshell+1):n_hypool_tot) = psi_node_shell_1D(:)
                 ! REPRESENTATIVE SINGLE FINE ROOT POOL (weighted average across layers)
                 !call map2d_to_1d_aroot(ft, ccohort, kmax_bound_bylayershell, ths_aroot_1D, thr_aroot_1D, vtot_aroot_1D, psi_node_aroot_1D)
                 !psi_node(n_hypool_ag+n_hypool_troot+1)            = psi_node_aroot_1D
              else if(site_hydr%nlevsoi_hyd == 1) then
                 write(fates_log(),*) 'Single layer hydraulics currently inoperative nlevsoi_hyd==1'
                 call endrun(msg=errMsg(sourcefile, __LINE__))
                 !psi_node(  (n_hypool_tot-nshell+1):n_hypool_tot) = psisoi_liq_shell(c,1,:)
                 !psi_node(  n_hypool_ag+n_hypool_troot+1)            = ccohort_hydr%psi_aroot(1)
                 !flc_min_node(n_hypool_ag+n_hypool_troot+1)          = ccohort_hydr%flc_min_aroot(1)
              end if
              
              ! SET NODE HEIGHTS AND VOLUMES
              z_node(                   1 : n_hypool_ag)           = ccohort_hydr%z_node_ag(:)        ! leaf and stem
              z_node(         (n_hypool_ag+1):(n_hypool_ag+n_hypool_troot)) = ccohort_hydr%z_node_troot(:)        ! transporting root
              z_node((n_hypool_ag+n_hypool_troot+1): n_hypool_tot)          = ccohort_hydr%z_node_aroot(1)     ! absorbing root and rhizosphere shells
              v_node(                   1 : n_hypool_ag)           = ccohort_hydr%v_ag(:)             ! leaf and stem
              v_node(         (n_hypool_ag+1):(n_hypool_ag+n_hypool_troot)) = ccohort_hydr%v_troot(:)             ! transporting root
              if(site_hydr%nlevsoi_hyd == 1) then
                 v_node((n_hypool_ag+n_hypool_troot+1)                    ) = ccohort_hydr%v_aroot_tot         ! absorbing root
                 v_node( (n_hypool_tot-nshell+1): n_hypool_tot)          = &
                       site_hydr%v_shell_1D(:)*ccohort_hydr%l_aroot_tot/sum(bc_in(s)%dz_sisl(:))  ! rhizosphere shells
              end if

		   ! SET SATURATED & RESIDUAL WATER CONTENTS
                   if(site_hydr%nlevsoi_hyd == 1) then
                      ths_node(  (n_hypool_tot-nshell+1):n_hypool_tot) = bc_in(s)%watsat_sisl(1)
		      !! BOC... should the below code exist on HLM side?  watres_col is a new SWC parameter 1
                      !  introduced for the van Genuchten, but does not exist for Campbell SWC.
                      select case (iswc)
                      case (van_genuchten) 
                         write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
                         call endrun(msg=errMsg(sourcefile, __LINE__)) 
                         ! thr_node(  (n_hypool_tot-nshell+1):n_hypool_tot) = bc_in(s)%watres_sisl(1)
                      case (campbell)
                         call swcCampbell_satfrac_from_psi(bc_in(s)%smpmin_si*denh2o*grav*1.e-9_r8, &
                                 (-1._r8)*bc_in(s)%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &
                                 bc_in(s)%bsw_sisl(1),     &
                                 tmp1)
                         call swcCampbell_th_from_satfrac(tmp1, &
                                 bc_in(s)%watsat_sisl(1),   &
                                 watres_local)
                         thr_node(  (n_hypool_tot-nshell+1):n_hypool_tot) = watres_local
                      case default
                      end select
                   end if
		   do k=1,n_hypool_ag+n_hypool_troot+1
		      ths_node(k) = EDPftvarcon_inst%hydr_thetas_node(ft,porous_media(k))
                      thr_node(k) = EDPftvarcon_inst%hydr_thetas_node(ft,porous_media(k)) * &
                                    EDPftvarcon_inst%hydr_resid_node(ft,porous_media(k))
		   enddo

		   ! SET BOUNDARY MAX CONDUCTANCES
		   !! assign cohort-level conductances to the 1D array
		   kmax_bound(                    :             ) = 0._r8
		   kmax_lower(                    :             ) = 0._r8
		   kmax_upper(                    :             ) = 0._r8
		   kmax_bound(                  1 : n_hypool_ag    ) = ccohort_hydr%kmax_bound(:)
		   kmax_upper(                  1 : n_hypool_ag    ) = ccohort_hydr%kmax_upper(:)
		   kmax_lower(                  1 : n_hypool_ag    ) = ccohort_hydr%kmax_lower(:)
		   kmax_upper((        n_hypool_ag+1)              ) = ccohort_hydr%kmax_upper_troot
                   if(site_hydr%nlevsoi_hyd == 1) then
                      !! estimate troot-aroot and aroot-radial components as a residual:
                      !! 25% each of total (surface of aroots to leaves) resistance
                      kmax_bound((        n_hypool_ag+1):(n_hypool_ag+2 )) = 2._r8 * ccohort_hydr%kmax_treebg_tot
                      kmax_lower((        n_hypool_ag+1)              ) = 2._r8 * kmax_bound(n_hypool_ag+1)
                      kmax_upper((        n_hypool_ag+2)              ) = 2._r8 * kmax_bound(n_hypool_ag+1)
                      kmax_lower((        n_hypool_ag+2)              ) = 2._r8 * ccohort_hydr%kmax_treebg_tot
                      kmax_bound_aroot_soil1                         = kmax_bound(n_hypool_ag+2)
                      kmax_bound_aroot_soil2                         = site_hydr%kmax_bound_shell_1D(1) * &
                            ccohort_hydr%l_aroot_tot / site_hydr%l_aroot_1D
                      kmax_bound((        n_hypool_ag+2)              ) = 1._r8/(1._r8/kmax_bound_aroot_soil1 + &
                            1._r8/kmax_bound_aroot_soil2)
                      kmax_bound((n_hypool_tot-nshell+1):(n_hypool_tot-1)) = site_hydr%kmax_bound_shell_1D(2:nshell) * &
                            ccohort_hydr%l_aroot_tot / site_hydr%l_aroot_1D
                      kmax_upper((n_hypool_tot-nshell+1):(n_hypool_tot  )) = site_hydr%kmax_upper_shell_1D(1:nshell) * &
                            ccohort_hydr%l_aroot_tot / site_hydr%l_aroot_1D
                      kmax_lower((n_hypool_tot-nshell+1):(n_hypool_tot  )) = site_hydr%kmax_lower_shell_1D(1:nshell) * &
                            ccohort_hydr%l_aroot_tot / site_hydr%l_aroot_1D
                   end if
		     
                   if(site_hydr%nlevsoi_hyd == 1) then
                      ! CONVERT WATER POTENTIALS TO WATER CONTENTS FOR THE NEW 'BUCKET' 
                      ! RHIZOSPHERE (fine roots and rhizosphere shells)
                      do k = (n_hypool_tot - nshell), n_hypool_tot
                         call th_from_psi(ft, porous_media(k), psi_node(k), th_node(k),site_hydr,bc_in(s))
                      enddo !aroot thru outer rhiz shell
                   end if

		   ! MAP REMAINING WATER CONTENTS (leaf, stem, troot) TO THE 1D ARRAY
		   th_node(               1 : n_hypool_ag          ) = ccohort_hydr%th_ag(:)
		   th_node(     (n_hypool_ag+1):(n_hypool_ag+n_hypool_troot)) = ccohort_hydr%th_troot(:)
		   flc_min_node(          1 : n_hypool_ag          ) = ccohort_hydr%flc_min_ag(:)
		   flc_min_node((n_hypool_ag+1):(n_hypool_ag+n_hypool_troot)) = ccohort_hydr%flc_min_troot(:)
		   
                      mono_decr_Rn = .true.
!		      do t=2, numLWPmem
!                         if((cpatch_hydr%netRad_mem(t) - cpatch_hydr%netRad_mem(t-1)) >= 0._r8) then
                            mono_decr_Rn = .false.
!                            EXIT
!                         end if
!		      end do

                   if(site_hydr%nlevsoi_hyd == 1) then
                      ! 1-D THETA-BASED SOLUTION TO RICHARDS' EQUATION
                      call Hydraulics_1DSolve(ccohort, ft, z_node, v_node, ths_node, &
                            thr_node, kmax_bound, kmax_upper, kmax_lower, &
                            kmax_bound_aroot_soil1, kmax_bound_aroot_soil2, &
                            th_node, flc_min_node, qflx_tran_veg_indiv, &
                            thresh, thresh_break, maxiter, imult, dtime, &
                            dth_node, the_node, we_area_outer, qtop_dt, dqtopdth_dthdt, &
                            sapflow, rootuptake, small_theta_num, &
                            mono_decr_Rn, site_hydr, bc_in(s))

                      ccohort_hydr%errh2o         = we_area_outer             ! kg/m2 ground/individual
                      ccohort_hydr%qtop_dt        = qtop_dt                   
                      ccohort_hydr%dqtopdth_dthdt = dqtopdth_dthdt            
                      ccohort_hydr%sapflow        = sapflow                   
                      ccohort_hydr%rootuptake     = rootuptake                

                      ! UPDATE WATER CONTENT & POTENTIAL IN LEAVES, STEM, AND TROOT (COHORT-LEVEL)  [[NOW THIS IS DONE BELOW]]
                      !do k=1,n_hypool_ag
                      !   ccohort_hydr%th_ag(k)          = th_node(k)
                      !   call psi_from_th(ft, porous_media(k), ccohort_hydr%th_ag(k), ccohort_hydr%psi_ag(k))
                      !enddo
                      !do k=(n_hypool_ag+1),(n_hypool_ag+n_hypool_troot)
                      !   ccohort_hydr%th_troot(k-n_hypool_ag) = th_node(k)
                      !   call psi_from_th(ft, porous_media(k), ccohort_hydr%th_troot(k-n_hypool_ag), ccohort_hydr%psi_troot(k-n_hypool_ag))
                      !enddo
		      ccohort_hydr%th_aroot(1)       = th_node(n_hypool_ag+n_hypool_troot+n_hypool_aroot)
                      call psi_from_th(ft, 4, ccohort_hydr%th_aroot(1), ccohort_hydr%psi_aroot(1), site_hydr, bc_in(s))
		      dwat_veg_coh              = sum(dth_node(1:n_hypool_ag+n_hypool_troot+n_hypool_aroot) * &
                                                  v_node(1:n_hypool_ag+n_hypool_troot+n_hypool_aroot)*denh2o)

		      site_hydr%dwat_veg         = site_hydr%dwat_veg + dwat_veg_coh*ccohort%n/AREA  !*patch_wgt

		      site_hydr%h2oveg           = site_hydr%h2oveg + dwat_veg_coh*ccohort%n/AREA !*patch_wgt
		      !site_hydr%errh2o_hyd       = site_hydr%errh2o_hyd + ccohort_hydr%errh2o*(ccohort%c_area / ccohort%n)/AREA !*patch_wgt
		      site_hydr%errh2o_hyd       = site_hydr%errh2o_hyd + ccohort_hydr%errh2o*ccohort%c_area /AREA !*patch_wgt
		
		      ! ACCUMULATE CHANGE IN SOIL WATER CONTENT OF EACH COHORT TO COLUMN-LEVEL
		      dth_layershell_col(site_hydr%nlevsoi_hyd,:) = dth_layershell_col(site_hydr%nlevsoi_hyd,:) + &
		                                          dth_node((n_hypool_tot-nshell+1):n_hypool_tot) * &
                                                          ccohort_hydr%l_aroot_tot * ccohort%n / site_hydr%l_aroot_1D! * &
							  !patch_wgt
		   else if(site_hydr%nlevsoi_hyd > 1) then
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
                      kbg_layer(:) = 0._r8
                      kbg_tot      = 0._r8
                      do j=1,site_hydr%nlevsoi_hyd
                         z_node_1l((n_hypool_ag+n_hypool_troot+1):(n_hypool_tot)) = bc_in(s)%z_sisl(j)
                         v_node_1l((n_hypool_ag+n_hypool_troot+1)           ) = ccohort_hydr%v_aroot_layer(j)
                         v_node_1l((n_hypool_tot-nshell+1):(n_hypool_tot)) = site_hydr%v_shell(j,:) * &
                               ccohort_hydr%l_aroot_layer(j)/bc_in(s)%dz_sisl(j)
                         kmax_bound_1l(:) = 0._r8 
                         kmax_bound_shell_1l(:) = site_hydr%kmax_bound_shell(j,:) * &
                                                  ccohort_hydr%l_aroot_layer(j) / site_hydr%l_aroot_layer(j)


                         ! transporting-to-absorbing root conductance: factor of 2 means one-half of the total 
                         ! belowground resistance in layer j      
                         kmax_bound_1l((n_hypool_ag+1)) = 2._r8 * ccohort_hydr%kmax_treebg_layer(j)                     
                         ! transporting-to-absorbing root conductance: factor of 2*2 means one-half of the total 
                         ! belowground resistance in layer j, split in half between transporting and absorbing root
		         kmax_lower_1l(n_hypool_ag+1) = 4._r8 * ccohort_hydr%kmax_treebg_layer(j)
                         ! radial absorbing root conductance: factor of 2 means one-half of the 
                         ! total belowground resistance in layer j
                         kmax_bound_aroot_soil1      = 2._r8 * ccohort_hydr%kmax_treebg_layer(j)    
                         ! (root surface)-to-(soil shell#1) conductance
		         kmax_bound_aroot_soil2      = kmax_bound_shell_1l(1) 
                         ! combined (soil shell#1)-to-(absorbing root) conductance

		         kmax_bound_1l(n_hypool_ag+2) = 1._r8/(1._r8/kmax_bound_aroot_soil1 + &
                                                       1._r8/kmax_bound_aroot_soil2)  
		         kmax_upper_1l(n_hypool_ag+2) = kmax_lower_1l(n_hypool_ag+1)
		         kmax_lower_1l(n_hypool_ag+2) = 2._r8 * ccohort_hydr%kmax_treebg_layer(j)
                         ! REMEMBER: kmax_bound_shell_1l defined at the uppper (closer to atmosphere) 
                         ! boundary for each node, while kmax_bound_1l defined at the lower 
                         ! (closer to bulk soil) boundary for each node
                         kmax_bound_1l(n_hypool_tot-nshell+1:n_hypool_tot-1) = kmax_bound_shell_1l(2:nshell)
                         kmax_upper_1l(n_hypool_tot-nshell+1:n_hypool_tot)   = &
                              site_hydr%kmax_upper_shell(j,1:nshell) * &
                              ccohort_hydr%l_aroot_layer(j) / site_hydr%l_aroot_layer(j)
           
                         kmax_lower_1l(n_hypool_tot-nshell+1:n_hypool_tot) = site_hydr%kmax_lower_shell(j,1:nshell) * &
                                  ccohort_hydr%l_aroot_layer(j) / site_hydr%l_aroot_layer(j)

		         th_node_1l(n_hypool_ag+n_hypool_troot+1) = ccohort_hydr%th_aroot(j)

		         th_node_1l(n_hypool_ag+n_hypool_troot+2:n_hypool_tot) = &
                                          site_hydr%h2osoi_liqvol_shell(j,1:nshell)

                         psi_node_1l(     :) = fates_huge
                         flc_node_1l(     :) = fates_huge
                         dflcdpsi_node_1l(:) = fates_huge
                         do k = (n_hypool_ag+n_hypool_troot+1), n_hypool_tot
                            call psi_from_th(ft, porous_media(k), th_node_1l(k), &
                                             psi_node_1l(k),site_hydr, bc_in(s))
                            call flc_from_psi(ft, porous_media(k), psi_node_1l(k), &
                                              flc_node_1l(k), site_hydr, bc_in(s))
                            call dflcdpsi_from_psi(ft, porous_media(k), psi_node_1l(k), &
                                                   dflcdpsi_node_1l(k), site_hydr, bc_in(s))
                         enddo
                         hdiff_bound_1l(  :) = fates_huge
                         k_bound_1l(      :) = fates_huge
                         dhdiffdpsi0_1l(  :) = fates_huge
                         dhdiffdpsi1_1l(  :) = fates_huge
                         dkbounddpsi0_1l( :) = fates_huge
                         dkbounddpsi1_1l( :) = fates_huge
			 
                         ! Get k_bound_1l
                         call boundary_hdiff_and_k(1, z_node_1l(pick_1l), psi_node_1l(pick_1l), & 
                               flc_node_1l(pick_1l), dflcdpsi_node_1l(pick_1l), &
                               kmax_bound_1l(pick_1l), kmax_upper_1l(pick_1l),  &
                               kmax_lower_1l(pick_1l), hdiff_bound_1l, k_bound_1l, dhdiffdpsi0_1l, &
                               dhdiffdpsi1_1l, dkbounddpsi0_1l, dkbounddpsi1_1l, &
                               kmax_bound_aroot_soil1, kmax_bound_aroot_soil2)
                         !! upper bound limited to size()-1 b/c of zero-flux outer boundary condition
                         kbg_layer(j)        = 1._r8/sum(1._r8/k_bound_1l(1:(size(k_bound_1l)-1)))   
                         kbg_tot             = kbg_tot + kbg_layer(j)

		      enddo !soil layer

                      ! order soil layers in terms of decreasing volumetric water content
                      ! algorithm same as that used in histFileMod.F90 to alphabetize history tape contents
                      do j = site_hydr%nlevsoi_hyd-1,1,-1
                         do jj = 1,j
                            if (kbg_layer(ordered(jj)) <= kbg_layer(ordered(jj+1))) then
                               tmp           = ordered(jj)
                               ordered(jj)   = ordered(jj+1)
                               ordered(jj+1) = tmp
                            end if
                         enddo
                      enddo
			 
                      !initialize state variables in leaves to transporting roots
                      z_node_1l     (1:n_hypool_ag+n_hypool_troot)   = z_node(1:n_hypool_ag+n_hypool_troot)
                      v_node_1l     (1:n_hypool_ag+n_hypool_troot)   = v_node(1:n_hypool_ag+n_hypool_troot)
                      ths_node_1l   (1:n_hypool_ag+n_hypool_troot+1) = ths_node(1:n_hypool_ag+n_hypool_troot+1)
                      thr_node_1l   (1:n_hypool_ag+n_hypool_troot+1) = thr_node(1:n_hypool_ag+n_hypool_troot+1)
                      kmax_bound_1l (1:n_hypool_ag)            = kmax_bound(1:n_hypool_ag)
                      kmax_upper_1l (1:n_hypool_ag+1)          = kmax_upper(1:n_hypool_ag+1)
                      kmax_lower_1l (1:n_hypool_ag)            = kmax_lower(1:n_hypool_ag)
                      th_node_1l    (1:n_hypool_ag+n_hypool_troot)   = th_node(1:n_hypool_ag+n_hypool_troot)
                      ccohort_hydr%errh2o                   = 0._r8
  		      ! do j=1,nlevsoi_hyd  ! replace j with ordered(jj) in order 
                      ! to go through soil layers in order of decreasing total root-soil conductance
  		      do jj=1,site_hydr%nlevsoi_hyd

		         !initialize state variables in absorbing roots and rhizosphere shells in each soil layer
                         !z_node_1l(   (         n_hypool_ag+1):(n_hypool_troot         )) = -bc_in(s)%z_sisl(ordered(jj))
                         !! BOC...ad-hoc assume no grav difference bewtween aroot and troot for each layer

                         z_node_1l  (n_hypool_ag+n_hypool_troot+1:n_hypool_tot) = -bc_in(s)%z_sisl(ordered(jj))
		         v_node_1l  (n_hypool_ag+n_hypool_troot+1)           = ccohort_hydr%v_aroot_layer(ordered(jj))
 		         v_node_1l  (n_hypool_tot-nshell+1:n_hypool_tot)  = site_hydr%v_shell(ordered(jj),:) * &
                                                                      ccohort_hydr%l_aroot_layer(ordered(jj))/&
                                                                      bc_in(s)%dz_sisl(ordered(jj))
                         ths_node_1l(n_hypool_tot-nshell+1:n_hypool_tot)  = bc_in(s)%watsat_sisl(ordered(jj))

                         !! BOC... should the below code exist on HLM side?  watres_col is a new 
                         !! SWC parameter introduced for the van Genuchten, but does not exist for Campbell SWC.

                         select case (iswc)
                         case (van_genuchten) 
                            write(fates_log(),*) &
                                 'Van Genuchten plant hydraulics is inoperable until further notice'
                            call endrun(msg=errMsg(sourcefile, __LINE__)) 
                            ! thr_node_1l(   (n_hypool_tot-nshell+1):(n_hypool_tot     )) = bc_in(s)%watres_sisl(ordered(jj))
                         case (campbell)
                            call swcCampbell_satfrac_from_psi(bc_in(s)%smpmin_si*denh2o*grav*1.e-9_r8, &
                                    (-1._r8)*bc_in(s)%sucsat_sisl(ordered(jj))*denh2o*grav*1.e-9_r8, &
                                    bc_in(s)%bsw_sisl(ordered(jj)),     &
                                    tmp1)
                            call swcCampbell_th_from_satfrac(tmp1, &
                                    bc_in(s)%watsat_sisl(ordered(jj)),   &
                                    watres_local)
                            thr_node_1l(   (n_hypool_tot-nshell+1):(n_hypool_tot     )) = watres_local
                         case default
                         end select
		   
                         kmax_bound_shell_1l(:) = site_hydr%kmax_bound_shell(ordered(jj),:) * &
                              ccohort_hydr%l_aroot_layer(ordered(jj)) / site_hydr%l_aroot_layer(ordered(jj))

                         kmax_bound_1l(n_hypool_ag+1) = 2.0_r8 * ccohort_hydr%kmax_treebg_layer(ordered(jj))   

                         ! transporting-to-absorbing root conductance: factor of 2 means 
                         ! one-half of the total belowground resistance in layer j
                         kmax_lower_1l(n_hypool_ag+1) = 4.0_r8 * ccohort_hydr%kmax_treebg_layer(ordered(jj))
                         
                         ! transporting-to-absorbing root conductance: factor of 2*2 means one-half of the total 
                         ! belowground resistance in layer j, split in half between transporting and absorbing root
                         kmax_bound_aroot_soil1    = 2.0_r8 * ccohort_hydr%kmax_treebg_layer(ordered(jj))   

                         ! radial absorbing root conductance: factor of 2 means one-half of 
                         ! the total belowground resistance in layer j
                         kmax_bound_aroot_soil2 = kmax_bound_shell_1l(1)

                         ! (root surface)-to-(soil shell#1) conductance
                         kmax_bound_1l(n_hypool_ag+2) = 1.0_r8 / &
                              (1._r8/kmax_bound_aroot_soil1 + 1._r8/kmax_bound_aroot_soil2)

                         ! combined (soil shell#1)-to-(absorbing root) conductance
                         kmax_upper_1l(n_hypool_ag+2) = kmax_lower_1l(n_hypool_ag+1)
                         kmax_lower_1l(n_hypool_ag+2) = 2.0_r8 * ccohort_hydr%kmax_treebg_layer(ordered(jj))
                         kmax_bound_1l(n_hypool_tot-nshell+1:n_hypool_tot-1) = kmax_bound_shell_1l(2:nshell)

                         ! REMEMBER: kmax_bound_shell_1l defined at the uppper 
                         ! (closer to atmosphere) boundary for each node, while kmax_bound_1l 
                         ! defined at the lower (closer to bulk soil) boundary for each node
                         kmax_upper_1l((n_hypool_tot-nshell+1 ):(n_hypool_tot        )) = &
                               site_hydr%kmax_upper_shell(ordered(jj),1:nshell) * &
                               ccohort_hydr%l_aroot_layer(ordered(jj)) / site_hydr%l_aroot_layer(ordered(jj))
                         kmax_lower_1l((n_hypool_tot-nshell+1 ):(n_hypool_tot        )) = &
                               site_hydr%kmax_lower_shell(ordered(jj),1:nshell) * &
                               ccohort_hydr%l_aroot_layer(ordered(jj)) / site_hydr%l_aroot_layer(ordered(jj))

		         flc_min_node(n_hypool_ag+n_hypool_troot+1)         = ccohort_hydr%flc_min_aroot(ordered(jj))
		         th_node_1l(n_hypool_ag+n_hypool_troot+1)           = ccohort_hydr%th_aroot(ordered(jj))
		         th_node_1l(n_hypool_ag+n_hypool_troot+2:n_hypool_tot) = site_hydr%h2osoi_liqvol_shell(ordered(jj),:)

                         ! the individual-layer Richards' equation solution
                         call Hydraulics_1DSolve(ccohort, ft, &
                               z_node_1l, v_node_1l, ths_node_1l, thr_node_1l, &
                               kmax_bound_1l, kmax_upper_1l, kmax_lower_1l, &
                               kmax_bound_aroot_soil1, kmax_bound_aroot_soil2, &
                               th_node_1l, flc_min_node, qflx_tran_veg_indiv, &
                               thresh, thresh_break, maxiter, imult, &
                               dtime*kbg_layer(ordered(jj))/kbg_tot, &
                               dth_node_1l, the_node_1l, we_area_outer, qtop_dt, &
                               dqtopdth_dthdt, sapflow, rootuptake, small_theta_num, &
                               mono_decr_Rn, site_hydr, bc_in(s))
                   
                         dwat_veg_coh                          = &
                               sum(dth_node_1l(1:n_hypool_ag+n_hypool_troot+n_hypool_aroot)* &
                               v_node_1l(1:n_hypool_ag+n_hypool_troot+n_hypool_aroot)*denh2o)
                         site_hydr%dwat_veg                 = site_hydr%dwat_veg + dwat_veg_coh*ccohort%n/AREA!*patch_wgt
                         site_hydr%h2oveg                   = site_hydr%h2oveg + dwat_veg_coh*ccohort%n/AREA!*patch_wgt
                         ccohort_hydr%errh2o                    = ccohort_hydr%errh2o + we_area_outer                                               
                         !! kg/m2 ground/individual
                         
                         site_hydr%errh2o_hyd               = site_hydr%errh2o_hyd + &
			                                     we_area_outer*ccohort%c_area /AREA!*patch_wgt
                                                             !we_area_outer*(ccohort%c_area / ccohort%n)/AREA!*patch_wgt
                         ccohort_hydr%qtop_dt                   = ccohort_hydr%qtop_dt  + qtop_dt                             ! 
                         ccohort_hydr%dqtopdth_dthdt            = ccohort_hydr%dqtopdth_dthdt + dqtopdth_dthdt                ! 
                         ccohort_hydr%sapflow                   = ccohort_hydr%sapflow + sapflow                              ! 
                         ccohort_hydr%rootuptake                = ccohort_hydr%rootuptake + rootuptake                        ! 
                         SELECT CASE (ordered(jj))  !! select soil layer
                            CASE (1)
                               ccohort_hydr%rootuptake01        = rootuptake
                            CASE (2)
                               ccohort_hydr%rootuptake02        = rootuptake
                            CASE (3)
                               ccohort_hydr%rootuptake03        = rootuptake
                            CASE (4)
                               ccohort_hydr%rootuptake04        = rootuptake
                            CASE (5)
                               ccohort_hydr%rootuptake05        = rootuptake
                            CASE (6)
                               ccohort_hydr%rootuptake06        = rootuptake
                            CASE (7)
                               ccohort_hydr%rootuptake07        = rootuptake
                            CASE (8)
                               ccohort_hydr%rootuptake08        = rootuptake
                            CASE (9)
                               ccohort_hydr%rootuptake09        = rootuptake
                            CASE (10)
                               ccohort_hydr%rootuptake10        = rootuptake
                            CASE DEFAULT
                         end SELECT
	     
		         ! UPDATE WATER CONTENT & POTENTIAL IN AROOT (COHORT-LEVEL)
		         ccohort_hydr%th_aroot(ordered(jj))     = th_node_1l(n_hypool_ag+n_hypool_troot+1)
                         call psi_from_th(ft, porous_media(n_hypool_ag+n_hypool_troot+1), &
                               ccohort_hydr%th_aroot(ordered(jj)), ccohort_hydr%psi_aroot(ordered(jj)), &
                               site_hydr, bc_in(s))
                         call flc_from_psi(ft, porous_media(n_hypool_ag+n_hypool_troot+1), &
                               ccohort_hydr%psi_aroot(ordered(jj)), ccohort_hydr%flc_aroot(ordered(jj)), &
                               site_hydr, bc_in(s))

		         ! ACCUMULATE CHANGE IN SOIL WATER CONTENT OF EACH COHORT TO COLUMN-LEVEL
                         dth_layershell_col(ordered(jj),:) = dth_layershell_col(ordered(jj),:) + &
                               dth_node_1l((n_hypool_tot-nshell+1):n_hypool_tot) * &
                               ccohort_hydr%l_aroot_layer(ordered(jj)) * &
                               ccohort%n / site_hydr%l_aroot_layer(ordered(jj)) !* &
                                                     !patch_wgt
		      enddo !soil layer
		   end if !nlevsoi_hyd > 1
		   
                   ! UPDATE WATER CONTENT & POTENTIAL IN LEAVES, STEM, AND TROOT (COHORT-LEVEL)
                   do k=1,n_hypool_ag
                      if(site_hydr%nlevsoi_hyd == 1) then
                         ccohort_hydr%th_ag(k)          = th_node(k)
                      else
                         ccohort_hydr%th_ag(k)          = th_node_1l(k)
                      endif
                      call psi_from_th(ft, porous_media(k), ccohort_hydr%th_ag(k), &
                            ccohort_hydr%psi_ag(k), site_hydr, bc_in(s) )
                      call flc_from_psi(ft, porous_media(k), ccohort_hydr%psi_ag(k), &
                            ccohort_hydr%flc_ag(k), site_hydr, bc_in(s) ) 
                   enddo
                   do k=(n_hypool_ag+1),(n_hypool_ag+n_hypool_troot)
                      if(site_hydr%nlevsoi_hyd == 1) then
                         ccohort_hydr%th_troot(k-n_hypool_ag) = th_node(k)
                      else
                         ccohort_hydr%th_troot(k-n_hypool_ag) = th_node_1l(k)
                      endif
                      call psi_from_th(ft, porous_media(k), ccohort_hydr%th_troot(k-n_hypool_ag), &
                            ccohort_hydr%psi_troot(k-n_hypool_ag), site_hydr, bc_in(s))
                      call flc_from_psi(ft, porous_media(k), ccohort_hydr%psi_troot(k-n_hypool_ag), &
                            ccohort_hydr%flc_troot(k-n_hypool_ag), site_hydr, bc_in(s))
                   enddo
                   
		   ! SET COHORT-LEVEL BTRAN FOR USE IN NEXT TIMESTEP
                   ! first update the leaf water potential memory
                   do t=2, numLWPmem
                      ccohort_hydr%lwp_mem(t-1)      = ccohort_hydr%lwp_mem(t)
                   end do
                   ccohort_hydr%lwp_mem(numLWPmem)   = ccohort_hydr%psi_ag(1)
                   call flc_gs_from_psi(cCohort, ccohort_hydr%psi_ag(1))
		   
		   refill_rate = -log(0.5)/(ccohort_hydr%refill_days*24._r8*3600._r8)   ! s-1
                   do k=1,n_hypool_ag
                      ccohort_hydr%flc_min_ag(k) = min(ccohort_hydr%flc_min_ag(k), ccohort_hydr%flc_ag(k))
                      if(ccohort_hydr%psi_ag(k) >= ccohort_hydr%refill_thresh .and. &
                            ccohort_hydr%flc_ag(k) > ccohort_hydr%flc_min_ag(k)) then   ! then refilling
                         ccohort_hydr%flc_min_ag(k) = ccohort_hydr%flc_ag(k) - &
                               (ccohort_hydr%flc_ag(k) - ccohort_hydr%flc_min_ag(k))*exp(-refill_rate*dtime)
                      end if
                   end do
                   do k=1,n_hypool_troot
                      ccohort_hydr%flc_min_troot(k) = min(ccohort_hydr%flc_min_troot(k), ccohort_hydr%flc_troot(k))
                      if(ccohort_hydr%psi_troot(k) >= ccohort_hydr%refill_thresh .and. &
                            ccohort_hydr%flc_troot(k) > ccohort_hydr%flc_min_troot(k)) then   ! then refilling
                         ccohort_hydr%flc_min_troot(k) = ccohort_hydr%flc_troot(k) - &
                               (ccohort_hydr%flc_troot(k) - ccohort_hydr%flc_min_troot(k))*exp(-refill_rate*dtime)
                      end if
                   end do
                   do j=1,site_hydr%nlevsoi_hyd
                      ccohort_hydr%flc_min_aroot(j) = min(ccohort_hydr%flc_min_aroot(j), ccohort_hydr%flc_aroot(j))
                      if(ccohort_hydr%psi_aroot(j) >= ccohort_hydr%refill_thresh .and. &
                            ccohort_hydr%flc_aroot(j) > ccohort_hydr%flc_min_aroot(j)) then   ! then refilling
                         ccohort_hydr%flc_min_aroot(j) = ccohort_hydr%flc_aroot(j) - &
                               (ccohort_hydr%flc_aroot(j) - ccohort_hydr%flc_min_aroot(j))*exp(-refill_rate*dtime)
                      end if
                   end do
                   
                   ccohort => ccohort%shorter
                enddo !cohort 

             cpatch => cpatch%younger
          enddo !patch
	  
          ! UPDATE THE COLUMN-LEVEL SOIL WATER CONTENT (LAYER x SHELL)
	  site_hydr%supsub_flag(:) = 999
          do j=1,site_hydr%nlevsoi_hyd
             !! BOC... should the below code exist on HLM side?  watres_col is a new SWC parameter 
             ! introduced for the van Genuchten, but does not exist for Campbell SWC.
             select case (iswc)
             case (van_genuchten) 
                write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
                call endrun(msg=errMsg(sourcefile, __LINE__)) 
                ! watres_local = bc_in(s)%watres_sisl(j)
             case (campbell)
                call swcCampbell_satfrac_from_psi(bc_in(s)%smpmin_si*denh2o*grav*1.e-9_r8, &
                      (-1._r8)*bc_in(s)%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                      bc_in(s)%bsw_sisl(j),     &
                      tmp1)
                call swcCampbell_th_from_satfrac(tmp1, &
                      bc_in(s)%watsat_sisl(j),   &
                      watres_local)
             case default
             end select
             do k=1,nshell
                if ((site_hydr%h2osoi_liqvol_shell(j,k)+dth_layershell_col(j,k)) > &
                      (bc_in(s)%watsat_sisl(j)-small_theta_num)) then
                   site_hydr%supsub_flag(j)           =  k
                   site_hydr%h2osoi_liqvol_shell(j,k) =  bc_in(s)%watsat_sisl(j)-small_theta_num
                else if ((site_hydr%h2osoi_liqvol_shell(j,k)+dth_layershell_col(j,k)) < &
                      (watres_local+small_theta_num)) then
                   site_hydr%supsub_flag(j)           = -k
                   site_hydr%h2osoi_liqvol_shell(j,k) =  watres_local+small_theta_num
                else
                   site_hydr%h2osoi_liqvol_shell(j,k) =  site_hydr%h2osoi_liqvol_shell(j,k) + &
                                                         dth_layershell_col(j,k)
                end if
             enddo

             ! Update the matric potential in the inner-most shell 
             ! (used for setting tissue potentials of new recruits)
             call swcCampbell_psi_from_th(site_hydr%h2osoi_liqvol_shell(j,1), &
                   bc_in(s)%watsat_sisl(j), (-1.0_r8)*bc_in(s)%sucsat_sisl(j)*denh2o*grav*1.e-9_r8, &
                   bc_in(s)%bsw_sisl(j), smp)
             site_hydr%psisoi_liq_innershell(j) = smp


             if(site_hydr%nlevsoi_hyd == 1) then
                bc_out(s)%qflx_soil2root_sisl(1:bc_in(s)%nlevsoil-1) = 0._r8

                ! qflx_rootsoi(c,bc_in(s)%nlevsoil)        = 
                ! -(sum(dth_layershell_col(j,:))*bc_in(s)%dz_sisl(j)*denh2o/dtime)

                bc_out(s)%qflx_soil2root_sisl(bc_in(s)%nlevsoil)     = &
                      -(sum(dth_layershell_col(j,:)*site_hydr%v_shell_1D(:)) * &
                      !site_hydr%l_aroot_1D/bc_in(s)%dz_sisl(j)/AREA*denh2o/dtime)- &   !BOC(10/02/2018)...error - should be plus  
                      site_hydr%l_aroot_1D/bc_in(s)%dz_sisl(j)/AREA*denh2o/dtime)+ &
		      site_hydr%recruit_w_uptake(site_hydr%nlevsoi_hyd)

!                h2osoi_liqvol = min(bc_in(s)%eff_porosity_sl(bc_in(s)%nlevsoil), &
!                     bc_in(s)%h2o_liq_sisl(bc_in(s)%nlevsoil)/(bc_in(s)%dz_sisl(bc_in(s)%nlevsoil)*denh2o))
                
                ! Save the amount of liquid soil water known to the model after root uptake
                site_hydr%h2osoi_liq_prev(site_hydr%nlevsoi_hyd) = bc_in(s)%h2o_liq_sisl(bc_in(s)%nlevsoil) - &
                      dtime*bc_out(s)%qflx_soil2root_sisl(bc_in(s)%nlevsoil)
                
             else
                !qflx_rootsoi(c,j)              = -(sum(dth_layershell_col(j,:))*bc_in(s)%dz_sisl(j)*denh2o/dtime)
                bc_out(s)%qflx_soil2root_sisl(j)               = &
                      -(sum(dth_layershell_col(j,:)*site_hydr%v_shell(j,:)) * &
                      site_hydr%l_aroot_layer(j)/bc_in(s)%dz_sisl(j)/AREA*denh2o/dtime)+ &
		      site_hydr%recruit_w_uptake(j)

                ! h2osoi_liqvol =  min(bc_in(s)%eff_porosity_sl(j), &
                !     bc_in(s)%h2o_liq_sisl(j)/(bc_in(s)%dz_sisl(j)*denh2o))
                
                ! Save the amount of liquid soil water known to the model after root uptake
                ! This calculation also assumes that 1mm of water is 1kg
                site_hydr%h2osoi_liq_prev(j) = bc_in(s)%h2o_liq_sisl(j) - &
                     dtime*bc_out(s)%qflx_soil2root_sisl(j)

	     end if
             
           enddo !site_hydr%nlevsoi_hyd
           !-----------------------------------------------------------------------
           ! mass balance check and pass the total stored vegetation water to HLM
           ! in order for it to fill its balance checks
	   totalrootuptake = 0.0_r8
	   totalqtop_dt = 0.0_r8
	   cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))
	     ccohort=>cpatch%tallest
	     do while(associated(ccohort))
                ccohort_hydr => ccohort%co_hydr
                !totalrootuptake = totalrootuptake + ccohort_hydr%rootuptake* ccohort%n/AREA
		totalqtop_dt= totalqtop_dt+  ccohort_hydr%qtop_dt* ccohort%n/AREA
                ccohort => ccohort%shorter
             enddo !cohort
	     cpatch => cpatch%younger
	   enddo !patch
	   !remove the recruitment water uptake as it has been added to prev_h2oveg 
	   totalrootuptake = sum(bc_out(s)%qflx_soil2root_sisl(:)- &
	                  site_hydr%recruit_w_uptake(:))*dtime
	   
           total_e = site_hydr%h2oveg-(prev_h2oveg + totalrootuptake - totalqtop_dt)
	       
	   site_hydr%h2oveg_hydro_err = site_hydr%h2oveg_hydro_err + total_e
	   
           bc_out(s)%plant_stored_h2o_si = site_hydr%h2oveg + site_hydr%h2oveg_dead - &
                                           site_hydr%h2oveg_growturn_err - &
                                           site_hydr%h2oveg_pheno_err-&
					   site_hydr%h2oveg_hydro_err

           
        enddo !site
      
  end subroutine Hydraulics_BC

  ! =====================================================================================


  subroutine AccumulateMortalityWaterStorage(csite,ccohort,delta_n)

     ! ---------------------------------------------------------------------------
     ! This subroutine accounts for the water bound in plants that have
     ! just died. This water is accumulated at the site level for all plants
     ! that die.
     ! In another routine, this pool is reduced as water vapor flux, and
     ! passed to the HLM.
     ! ---------------------------------------------------------------------------
     use EDTypesMod        , only : AREA	

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
           sum(ccohort_hydr%th_troot(:)*ccohort_hydr%v_troot(:))             + &
           sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
           denh2o*delta_n/AREA
     
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
     use EDTypesMod, only : AREA

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
                       sum(ccohort_hydr%th_troot(:)*ccohort_hydr%v_troot(:)) + &
                       sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)))* &
                       denh2o*currentCohort%n
              end if
              currentCohort => currentCohort%shorter
           enddo !cohort
           currentPatch => currentPatch%younger
        enddo !end patch loop
        
        csite_hydr%h2oveg_recruit      = csite_hydr%h2oveg_recruit      / AREA

     end do
     
     return
  end subroutine RecruitWaterStorage

  
  !-------------------------------------------------------------------------------!
  
  subroutine Hydraulics_1DSolve(cc_p, ft, z_node, v_node, ths_node, thr_node, kmax_bound, &
                                kmax_upper, kmax_lower, kmax_bound_aroot_soil1, kmax_bound_aroot_soil2, &
                                th_node, flc_min_node, qtop, thresh, thresh_break, maxiter, imult, &
                                dtime, dth_node_outer, the_node, we_area_outer, qtop_dt, &
                                dqtopdth_dthdt, sapflow, rootuptake, small_theta_num, mono_decr_Rn, &
                                site_hydr, bc_in)

    use EDTypesMod        , only : AREA				
    !
    ! !DESCRIPTION: 
    !
    !
    ! !ARGUMENTS
    type(ed_cohort_type) , intent(inout), target  :: cc_p                          ! current cohort pointer
    integer  , intent(in)    :: ft                     ! PFT index
    real(r8) , intent(in)    :: z_node(:)              ! nodal height of water storage compartments                      [m]
    real(r8) , intent(in)    :: v_node(:)              ! volume of water storage compartments                            [m3]
    real(r8) , intent(in)    :: ths_node(:)            ! saturated volumetric water in water storage compartments        [m3 m-3]
    real(r8) , intent(in)    :: thr_node(:)            ! residual volumetric water in water storage compartments         [m3 m-3]
    real(r8) , intent(in)    :: kmax_bound(:)          ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8) , intent(in)    :: kmax_upper(:)          ! maximum hydraulic conductance from node to upper boundary       [kg s-1 MPa-1]
    real(r8) , intent(in)    :: kmax_lower(:)          ! maximum hydraulic conductance from node to lower boundary       [kg s-1 MPa-1]
    real(r8) , intent(in)    :: kmax_bound_aroot_soil1 ! maximum radial conductance of absorbing roots                   [kg s-1 MPa-1]
    real(r8) , intent(in)    :: kmax_bound_aroot_soil2 ! maximum conductance to root surface from innermost rhiz shell   [kg s-1 MPa-1]
    real(r8) , intent(inout) :: th_node(:)             ! volumetric water in water storage compartments                  [m3 m-3]
    real(r8) , intent(in)    :: flc_min_node(:)        ! minimum attained fractional loss of conductivity (for xylem refilling dynamics) [-]
    real(r8) , intent(in)    :: qtop                   ! evaporative flux from canopy                                    [kgh2o indiv-1 s-1]
    integer  , intent(in)    :: maxiter                ! maximum iterations for timestep reduction                       [-]
    integer  , intent(in)    :: imult                  ! iteration index multiplier                                      [-]
    real(r8) , intent(in)    :: thresh                 ! threshold for water balance error (warning only)                [mm h2o]
    real(r8) , intent(in)    :: thresh_break           ! threshold for water balance error (stop model)                  [mm h2o]
    real(r8) , intent(in)    :: dtime                  ! timestep size                                                   [s]
    real(r8) , intent(out)   :: dth_node_outer(:)      ! change in volumetric water in water storage compartments        [m3 m-3]
    real(r8) , intent(out)   :: the_node(:)            ! error resulting from supersaturation or below-residual th_node  [m3 m-3]
    real(r8) , intent(out)   :: we_area_outer          ! 1D plant-soil continuum water error                             [kgh2o m-2]
    real(r8) , intent(out)   :: qtop_dt
    real(r8) , intent(out)   :: dqtopdth_dthdt
    real(r8) , intent(out)   :: sapflow
    real(r8) , intent(out)   :: rootuptake
    real(r8) , intent(in)    :: small_theta_num        ! avoids theta values equalling thr or ths                        [m3 m-3]
    logical  , intent(in)    :: mono_decr_Rn           ! flag indicating whether Rn is monotonically decreasing
    type(ed_site_hydr_type), intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(bc_in_type), intent(in)             :: bc_in       ! FATES boundary conditions

    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type) , pointer                :: ccohort                 !
    integer  :: k                             ! 1D plant-soil continuum array
    integer  :: iterh1, iterh2                ! iteration indices                                               [-]
    real(r8) :: w_tot_beg_inner
    real(r8) :: w_tot_end_inner
    real(r8) :: dw_tot_inner
    real(r8) :: w_tot_beg_outer
    real(r8) :: w_tot_end_outer
    real(r8) :: dw_tot_outer
    real(r8) :: we_tot_inner
    real(r8) :: we_area_inner
    real(r8) :: we_vol_inner
    real(r8) :: we_local
    real(r8) :: we_tot_outer
    integer  :: dt_fac                        ! timestep divisor                                                [-]
    real(r8) :: dt_fac2                       ! timestep divisor                                                [-]
    real(r8) :: dt_new                        ! new timestep                                                    [s]
    real(r8) :: th_node_init( n_hypool_tot)      ! initial volumetric water in water storage compartments          [m3 m-3]
    real(r8) :: psi_node(     n_hypool_tot)      ! water potential in water storage compartments                   [MPa]
    real(r8) :: dpsidth_node( n_hypool_tot)      ! derivative of water potential wrt to theta                      [MPa]
    real(r8) :: flc_node(     n_hypool_tot)      ! fractional loss of conductivity at water storage nodes          [-]
    real(r8) :: dflcdpsi_node(n_hypool_tot)      ! derivative of fractional loss of conductivity wrt psi           [MPa-1]
    real(r8) :: k_bound(      n_hypool_tot)      ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
    real(r8) :: q_bound(      n_hypool_tot)      ! lower boundary flux rate                                        [kg s-1]
    real(r8) :: hdiff_bound(  n_hypool_tot)      ! total water potential difference across lower boundary          [MPa-1]
    real(r8) :: dhdiffdpsi0(  n_hypool_tot)      ! derivative of total water potential difference wrt psi above    [-]
    real(r8) :: dhdiffdpsi1(  n_hypool_tot)      ! derivative of total water potential difference wrt psi below    [-]
    real(r8) :: dkbounddpsi0( n_hypool_tot)      ! derivative of lower boundary conductance wrt psi above          [kg s-1 MPa-2]
    real(r8) :: dkbounddpsi1( n_hypool_tot)      ! derivative of lower boundary conductance wrt psi below          [kg s-1 MPa-2]
    real(r8) :: dqbounddpsi0( n_hypool_tot)      ! derivative of lower boundary flux rate wrt psi above            [kg s-1 MPa-1]
    real(r8) :: dqbounddpsi1( n_hypool_tot)      ! derivative of lower boundary flux rate wrt psi below            [kg s-1 MPa-1]
    real(r8) :: dqbounddth0(  n_hypool_tot)      ! derivative of lower boundary flux rate wrt theta above          [kg s-1 m3 m-3]
    real(r8) :: dqbounddth1(  n_hypool_tot)      ! derivative of lower boundary flux rate wrt theta below          [kg s-1 m3 m-3]
    real(r8) :: dth_node_inner(n_hypool_tot)     ! dtheta from inner do while loop                                 [m3 m-3]
    real(r8) :: amx(n_hypool_tot)                ! "a" left off diagonal of tridiagonal matrix                     [kg s-1]
    real(r8) :: bmx(n_hypool_tot)                ! "b" diagonal of tridiagonal matrix                              [kg s-1]
    real(r8) :: cmx(n_hypool_tot)                ! "c" right off diagonal of tridiagonal matrix                    [kg s-1]
    real(r8) :: rmx(n_hypool_tot)                ! "r" forcing term of tridiagonal matrix                          [kg s-1]
    real(r8) :: supsub_flag_node(n_hypool_tot)   ! super saturation or sub residual flag                           [0-no,1-yes]
    real(r8) :: dflcgsdpsi                    ! derivative of stomatal vuln curve wrt to leaf water potential   [MPa-1]
    real(r8) :: dflcgsdth                     ! derivative of stomatal vuln curve wrt to leaf water content     [m-3 m3]
    real(r8) :: dqtopdflcgs                   ! derivative of cohort-level transpiration wrt btran              [kgh2o indiv-1 s-1]
    real(r8) :: dqtopdth_leaf                 ! derivative of transpiration rate wrt to leaf water content      [kgh2o indiv-1 s-1 m-3 m-3]
    real(r8) :: th_prev                       ! temporary                                                       [m3 m-3]
    real(r8) :: dth_prev                      ! temporary                                                       [m3 m-3]
    real(r8) :: dw_total                      ! temporary                                                       [kg]
    real(r8) :: we_k                          ! error for kth node (temporary)                                  [kg]
    real(r8) :: the_k                         ! theta error for kth node (temporary)                            [m3 m-3]
    real(r8) :: supsub_adj_w                  ! water ajustment due to super saturation or sub residual flag 
    logical  :: catch_nan                     ! flag for nan returned from Tridiagaonal
    integer  :: index_nan                     ! highest k index possessing a nan
    integer  :: index_stem
    integer  :: index_aroot
    integer  :: supsub_flag = 0
    integer  :: max_l                         !location of maximum water storage in the array
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
     !----------------------------------------------------------------------

    ccohort           => cc_p
    ccohort_hydr      => ccohort%co_hydr
    dth_node_outer(:) =  0._r8
    we_local          =  0._r8
    supsub_flag_node  =  0._r8
    index_stem        =  n_hypool_leaf + n_hypool_stem
    index_aroot       =  n_hypool_leaf + n_hypool_stem + n_hypool_troot + n_hypool_aroot

    ! STORE INITIAL STATES
    !! in case timestep needs to be chopped in half to balance water
    th_node_init(:) = th_node(:)
		
    ! WATER BALANCE
    w_tot_beg_outer = sum(th_node(:)*v_node(:))*denh2o
    
    ! OUTER DO-WHILE LOOP
    !! cuts timestep in half until all sub-timesteps (inner do-while loop) balance the water budget
    iterh1 = 0
    do while( iterh1 == 0 .or. ((abs(we_local) > thresh .or. supsub_flag /= 0) .and. iterh1 < maxiter) )
       dt_fac = max(imult*iterh1,1)
       dt_fac2 = real(dt_fac,r8)
       dt_new = dtime/dt_fac2

       !! restore initial states for a fresh attempt using new sub-timesteps
       if(iterh1 .gt. 0) then
          th_node(:) = th_node_init(:)
       end if
		     
       ! QUANTITIES OF INTEREST
       qtop_dt        = 0._r8
       dqtopdth_dthdt = 0._r8
       sapflow        = 0._r8
       rootuptake     = 0._r8
	  
       ! INNER DO-WHILE LOOP
       !! repeats, for dt_frac times, the removal of 1/dt_fac * transpiration 
       !! (the top boundary flux condition)
       !! stops and returns to outer loop if at any point the water budget 
       !! doesn't balance, so the timestep can be chopped in half again
       
       iterh2 = 0
       we_local = 0._r8
       supsub_flag = 0
       do while( iterh2 < dt_fac .and. ((abs(we_local) <= thresh .and. &
                 supsub_flag == 0) .or. iterh1 == (maxiter-1)))
          iterh2 = iterh2 + 1

          ! SET DERIVED STATE VARIABLES OVER ALL NODES
          do k = 1, n_hypool_tot
             call psi_from_th(ft, porous_media(k), th_node(k), psi_node(k), site_hydr, bc_in)
             call dpsidth_from_th(ft, porous_media(k), th_node(k), dpsidth_node(k), site_hydr, bc_in)
             call flc_from_psi(ft, porous_media(k), psi_node(k), flc_node(k), site_hydr, bc_in)
             call dflcdpsi_from_psi(ft, porous_media(k), psi_node(k), dflcdpsi_node(k), site_hydr, bc_in)

             if(do_dyn_xylemrefill .and. porous_media(k) <= 4) then
                if(flc_node(k) > flc_min_node(k)) then
                   dflcdpsi_node(k) = 0._r8
		   flc_node(k)      = flc_min_node(k)
                end if
             end if
          enddo
          call dflcgsdpsi_from_psi(ccohort_hydr%psi_ag(1),ft, dflcgsdpsi)
          dflcgsdth   = dflcgsdpsi * dpsidth_node(1)
          dqtopdflcgs = 0.1411985_r8    
          
          !BOC... estimated by trial-and-error: this term multiplied by the maximum value of dflcgsdth gives 150.
          !NEEDED: an estimate for dqtopdflcgs that accounts for variable potential evapotranspiration (efpot).
		  
          ! SET THE DQTOPDTH_LEAF TERM
          if(do_dqtopdth_leaf) then
             dqtopdth_leaf = dqtopdflcgs * dflcgsdth
!             if(mono_decr_Rn .and. ccohort_hydr%psi_ag(1) >= -0.88_r8) then
!                dqtopdth_leaf = 0._r8
!             else if(ccohort_hydr%psi_ag(1) < -0.88_r8) then
!                dqtopdth_leaf = 150._r8
!             else
!                dqtopdth_leaf = 0._r8
!             end if
          else
             dqtopdth_leaf = 0._r8
          end if

          ! SET BOUNDARY PRESSURE DIFFERENCES & CONDUCTANCES
          !! compute water potential differences + conductances and their derivatives wrt water potential
          call boundary_hdiff_and_k((n_hypool_tot-nshell), z_node, psi_node, flc_node, dflcdpsi_node, &
                                    kmax_bound, kmax_upper, kmax_lower, hdiff_bound, k_bound, dhdiffdpsi0, &
				    dhdiffdpsi1, dkbounddpsi0, dkbounddpsi1, &
	                            kmax_bound_aroot_soil1, kmax_bound_aroot_soil2)

          ! SET BOUNDARY FLUX TERMS
          !! compute flux terms and their derivatives wrt water content
          q_bound(1:n_hypool_tot-1)      = -1._r8 * k_bound(1:n_hypool_tot-1) * hdiff_bound(1:n_hypool_tot-1)
          dqbounddpsi0(1:n_hypool_tot-1) = -1._r8 * k_bound(1:n_hypool_tot-1) * dhdiffdpsi0(1:n_hypool_tot-1) - &
                                        dkbounddpsi0(1:n_hypool_tot-1) * hdiff_bound(1:n_hypool_tot-1)
          dqbounddpsi1(1:n_hypool_tot-1) = -1._r8 * k_bound(1:n_hypool_tot-1) * dhdiffdpsi1(1:n_hypool_tot-1) - &
                                        dkbounddpsi1(1:n_hypool_tot-1) * hdiff_bound(1:n_hypool_tot-1)
          dqbounddth0(1:n_hypool_tot-1)  = dqbounddpsi0(1:n_hypool_tot-1) * dpsidth_node(1:n_hypool_tot-1)
          dqbounddth1(1:n_hypool_tot-1)  = dqbounddpsi1(1:n_hypool_tot-1) * dpsidth_node(2:n_hypool_tot)

          !! zero-flux outer soil shell boundary condition
          q_bound(n_hypool_tot)      =    0._r8
          dqbounddpsi0(n_hypool_tot) =    0._r8
          dqbounddpsi1(n_hypool_tot) =    0._r8
          dqbounddth0(n_hypool_tot)  =    0._r8
          dqbounddth1(n_hypool_tot)  =    0._r8

          ! STORE BEGINNING WATER BALANCE
          w_tot_beg_inner = sum(th_node(:)*v_node(:))*denh2o

          ! SET UP TRIDIAGONAL MATRIX
          !! upper (leaf) layer
          k = 1
          rmx(k)    =  qtop - q_bound(k)
          amx(k)    =  0._r8
          bmx(k)    =  dqbounddth0(k) - dqtopdth_leaf - v_node(k)*denh2o/dt_new
          cmx(k)    =  dqbounddth1(k)
          !! intermediate nodes (plant and soil)
          do k=2,(n_hypool_tot-1)
             rmx(k) =  q_bound(k-1) - q_bound(k)
	     amx(k) = -1._r8 * dqbounddth0(k-1)
	     bmx(k) =  dqbounddth0(k) - dqbounddth1(k-1) - v_node(k)*denh2o/dt_new
	     cmx(k) =  dqbounddth1(k)
          enddo
          !! outermost rhizosphere shell
          k = n_hypool_tot
          rmx(k)    =  q_bound(k-1)
          amx(k)    = -1._r8 * dqbounddth0(k-1)
          bmx(k)    = -dqbounddth1(k-1) - v_node(k)*denh2o/dt_new
          cmx(k)    =  0._r8
		      
          ! SOLVE TRIDIAGONAL MATRIX
          call Hydraulics_Tridiagonal(amx, bmx, cmx, rmx, dth_node_inner)
	  
	  ! CATCH NAN VALUES
	  catch_nan = .false.
	  index_nan = 999
	  do k = 1, n_hypool_tot
	     if(isnan(dth_node_inner(k))) then
	        catch_nan = .true.
		index_nan = k
	     end if
	  end do
          if(catch_nan) then
             write(fates_log(),*)'EDPlantHydraulics returns nan at k = ', char(index_nan)
             call endrun(msg=errMsg(sourcefile, __LINE__))
	  end if
	  
	  ! CATCH SUPERSATURATED OR SUB-RESIDUAL WATER CONTENTS
	  ccohort_hydr%supsub_flag = 0._r8
	  supsub_flag = 0
          do k=1,n_hypool_tot
             th_prev     = th_node(k)
	     dth_prev    = dth_node_inner(k)
	     if(     (th_node(k)+dth_node_inner(k)) > (ths_node(k)-small_theta_num)) then
	        supsub_flag         =  k
		supsub_flag_node(k) =  1._r8
	        ccohort_hydr%supsub_flag =  real(k)
		th_node(k)  = ths_node(k)-small_theta_num
	     else if((th_node(k)+dth_node_inner(k)) < (thr_node(k)+small_theta_num)) then
	        supsub_flag         = -k
	        ccohort_hydr%supsub_flag = -real(k)
		supsub_flag_node(k) =  -1._r8
		th_node(k)  = thr_node(k)+small_theta_num
	     else
                th_node(k)  = th_node(k)+dth_node_inner(k)
	     end if
	     dth_node_inner(k) = th_node(k) - th_prev
	     the_node(k) = dth_node_inner(k) - dth_prev
          enddo
	  
	  ! QUANTITIES OF INTEREST
	  qtop_dt         = qtop_dt + &
	                    qtop*dt_new
	  dqtopdth_dthdt  = dqtopdth_dthdt + &
	                    dqtopdth_leaf*dth_node_inner(1)*dt_new
	  sapflow         = sapflow + &
	                    (q_bound(index_stem) + &
	                     dqbounddth0(index_stem)*dth_node_inner(index_stem) + &
		             dqbounddth1(index_stem)*dth_node_inner(index_stem+1))*dt_new
	  rootuptake      = rootuptake + &
	                    (q_bound(index_aroot) + &
	                     dqbounddth0(index_aroot)*dth_node_inner(index_aroot) + &
			     dqbounddth1(index_aroot)*dth_node_inner(index_aroot+1))*dt_new
	  
	  ! UPDATE ERROR TERM
          w_tot_end_inner = sum(th_node(:)*v_node(:))*denh2o
	  dw_tot_inner    = w_tot_end_inner - w_tot_beg_inner
          we_tot_inner    = dw_tot_inner/dt_new + (qtop + dqtopdth_leaf*dth_node_inner(1))
	  we_area_inner   = we_tot_inner/(cCohort%c_area / cCohort%n)
	  we_vol_inner    = we_tot_inner/sum(v_node(:))

          ! we_vol_inner 
          ! different water balance metrics can be chosen here (with an appropriate corresponding thresh)
	  we_local        = we_tot_inner*cCohort%n/AREA     
			
       end do ! loop over sub-timesteps

       iterh1 = iterh1 + 1
		     
    end do ! loop to get a timestep divisor that balances water
    
    ccohort_hydr%iterh1 = real(iterh1)
    ccohort_hydr%iterh2 = real(iterh2)

    ! WATER BALANCE ERROR-HANDLING
    if ( (abs(we_local) > thresh) .and. debug) then
       write(fates_log(),*)'WARNING: plant hydraulics water balance error exceeds threshold of ',&
             thresh
    else if (abs(we_local) > thresh_break) then
       write(fates_log(),*)'EDPlantHydraulics water balance error exceeds threshold of = ', thresh_break
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    ! TOTAL NET WATER BALANCE AND ERROR ADJUST HACK
    w_tot_end_outer   = sum(th_node(:)*v_node(:))*denh2o                            ! kg
    dw_tot_outer      = w_tot_end_outer - w_tot_beg_outer                           ! kg/timestep
    we_tot_outer      = dw_tot_outer + (qtop_dt + dqtopdth_dthdt)                   ! kg/timestep
    we_area_outer     = we_tot_outer/(cCohort%c_area / cCohort%n)                   ! kg/m2 ground/individual
    if(abs(we_tot_outer*cCohort%n)/AREA>1.0e-7_r8) then
      if(debug) then
          write(fates_log(),*)'WARNING: plant hydraulics water balance error exceeds 1.0e-7 and is ajusted for error'
      endif
      !dump the error water to the bin with largest water storage
      max_l  = maxloc(th_node(:)*v_node(:),dim=1)
      th_node(max_l) = th_node(max_l)-  &
                        we_tot_outer/(v_node(max_l)*denh2o)
      th_node(max_l) = min (th_node(max_l),&
                         ths_node(max_l)-small_theta_num) 
      th_node(max_l) = max(th_node(max_l),&
                         thr_node(max_l)+small_theta_num)	  
      w_tot_end_outer   = sum(th_node(:)*v_node(:))*denh2o                            ! kg
      dw_tot_outer      = w_tot_end_outer - w_tot_beg_outer                           ! kg/timestep
      we_tot_outer      = dw_tot_outer + (qtop_dt + dqtopdth_dthdt)                   ! kg/timestep
      we_area_outer     = we_tot_outer/(cCohort%c_area / cCohort%n)                   ! kg/m2 ground/individual   
    end if
    dth_node_outer(:) = th_node(:) - th_node_init(:)

 end subroutine Hydraulics_1DSolve

  !-------------------------------------------------------------------------------!
 subroutine Hydraulics_Tridiagonal(a, b, c, r, u)
    !
    ! !DESCRIPTION: An abbreviated version of biogeophys/TridiagonalMod.F90
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8), intent(in)    :: a(:)           ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(:)           ! "b" diagonal column of tridiagonal matrix
    real(r8), intent(in)    :: c(:)           ! "c" right off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: r(:)           ! "r" forcing term of tridiagonal matrix
    real(r8), intent(out)   :: u(:)           ! solution
    !
    ! !LOCAL VARIABLES:
    real(r8) :: bet                           ! temporary
    real(r8) :: gam(n_hypool_tot)                ! temporary
    integer  :: k                             ! index
    !----------------------------------------------------------------------

    bet = b(1)
    do k=1,n_hypool_tot
       if(k == 1) then
          u(k)   = r(k) / bet
       else
          gam(k) = c(k-1) / bet
	  bet    = b(k) - a(k) * gam(k)
	  u(k)   = (r(k) - a(k)*u(k-1)) / bet
       end if
    enddo
  
    do k=n_hypool_tot-1,1,-1
          u(k)   = u(k) - gam(k+1) * u(k+1)
    enddo
    
 end subroutine Hydraulics_Tridiagonal

  !-------------------------------------------------------------------------------!
  subroutine boundary_hdiff_and_k(k_arootsoil, z_node, psi_node, flc_node, dflcdpsi_node, &
                                  kmax_bound, kmax_upper, kmax_lower, hdiff_bound, k_bound, &
                                  dhdiffdpsi0, dhdiffdpsi1, dkbounddpsi0, dkbounddpsi1, &
				  kmax_bound_aroot_soil1, kmax_bound_aroot_soil2)
    !
    ! !ARGUMENTS
    integer           , intent(in)  :: k_arootsoil            ! index of node where the boundary occurs between root and soil
    real(r8)          , intent(in)  :: z_node(:)              ! height of node                                                  [m]
    real(r8)          , intent(in)  :: psi_node(:)            ! water potential in water storage compartments                   [MPa]
    real(r8)          , intent(in)  :: flc_node(:)            ! fractional loss of conductivity at water storage nodes          [-]
    real(r8)          , intent(in)  :: dflcdpsi_node(:)       ! derivative of fractional loss of conductivity wrt psi           [MPa-1]
    real(r8)          , intent(in)  :: kmax_bound(:)          ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8)          , intent(in)  :: kmax_upper(:)          ! maximum hydraulic conductance from node to upper boundary       [kg s-1 MPa-1]
    real(r8)          , intent(in)  :: kmax_lower(:)          ! maximum hydraulic conductance from node to lower boundary       [kg s-1 MPa-1]
    real(r8)          , intent(out) :: hdiff_bound(:)         ! total water potential difference across lower boundary          [MPa-1]
    real(r8)          , intent(out) :: k_bound(:)             ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
    real(r8)          , intent(out) :: dhdiffdpsi0(:)         ! derivative of total water potential difference wrt psi above    [-]
    real(r8)          , intent(out) :: dhdiffdpsi1(:)         ! derivative of total water potential difference wrt psi below    [-]
    real(r8)          , intent(out) :: dkbounddpsi0(:)        ! derivative of lower boundary conductance wrt psi above          [kg s-1 MPa-2]
    real(r8)          , intent(out) :: dkbounddpsi1(:)        ! derivative of lower boundary conductance wrt psi below          [kg s-1 MPa-2]
    real(r8), optional, intent(in)  :: kmax_bound_aroot_soil1 ! maximum radial conductance of absorbing roots                   [kg s-1 MPa-1]
    real(r8), optional, intent(in)  :: kmax_bound_aroot_soil2 ! maximum conductance to root surface from innermost rhiz shell   [kg s-1 MPa-1]
    !
    ! !LOCAL VARIABLES:
    integer  :: k                                        ! shell index
    real(r8) :: k_bound_aroot_soil1                      ! radial conductance of absorbing roots                           [kg s-1 MPa-1]
    real(r8) :: k_bound_aroot_soil2                      ! conductance to root surface from innermost rhiz shell           [kg s-1 MPa-1]
    real(r8) :: k_lower                                  ! conductance node k to lower boundary                            [kg s-1 MPa-1]
    real(r8) :: k_upper                                  ! conductance node k+1 to upper boundary                          [kg s-1 MPa-1]
    !----------------------------------------------------------------------

    do k = 1, (size(z_node)-1)
       hdiff_bound(k) = 1.e-6_r8*denh2o*grav*(z_node(k) - z_node(k+1)) + &
                        (psi_node(k) - psi_node(k+1))
       if(do_kbound_upstream) then

          ! absorbing root-1st rhizosphere shell boundary. 
          ! Comprised of two distinct conductance terms each with distinct water potentials

          if(k == (k_arootsoil)) then  
             
             k_bound_aroot_soil1 =  kmax_bound_aroot_soil1 * flc_node(k)
             k_bound_aroot_soil2 =  kmax_bound_aroot_soil2 * flc_node(k+1)
             k_bound(k)          =  1._r8/(1._r8/k_bound_aroot_soil1 + 1._r8/k_bound_aroot_soil2)
             dkbounddpsi0(k)     =  ((k_bound(k)/k_bound_aroot_soil1)**2._r8) * & 
                                    kmax_bound_aroot_soil1*dflcdpsi_node(k)
             dkbounddpsi1(k)     =  ((k_bound(k)/k_bound_aroot_soil2)**2._r8) * &
                                    kmax_bound_aroot_soil2*dflcdpsi_node(k+1)
          else
             ! examine direction of water flow; use the upstream node's k for the boundary k.
             ! (as suggested by Ethan Coon, LANL)
             if(hdiff_bound(k) < 0._r8) then
	        k_bound(k)       = kmax_bound(k) * flc_node(k+1)  ! water moving towards atmosphere
	        dkbounddpsi0(k)  = 0._r8
                dkbounddpsi1(k)  = kmax_bound(k) * dflcdpsi_node(k+1)
             else                                           
	        k_bound(k)       = kmax_bound(k) * flc_node(k)    ! water moving towards soil
	        dkbounddpsi0(k)  = kmax_bound(k) * dflcdpsi_node(k)
                dkbounddpsi1(k)  = 0._r8
             end if
          end if
       else
          k_lower                =  kmax_lower(k)   * flc_node(k)
          k_upper                =  kmax_upper(k+1) * flc_node(k+1)
          k_bound(k)             =  1._r8/(1._r8/k_lower + 1._r8/k_upper)
          dkbounddpsi0(k)        =  ((k_bound(k)/k_lower)**2._r8) * kmax_lower(k)  * dflcdpsi_node(k)
          dkbounddpsi1(k)        =  ((k_bound(k)/k_upper)**2._r8) * kmax_upper(k+1)* dflcdpsi_node(k+1)
       end if
       dhdiffdpsi0(k)  =  1.0
       dhdiffdpsi1(k)  = -1.0
    enddo
    k               = size(z_node)
    k_bound(k)      = 0._r8
    dkbounddpsi0(k) = 0._r8
    dkbounddpsi1(k) = 0._r8
  
  end subroutine boundary_hdiff_and_k
  
  !-------------------------------------------------------------------------------!
  subroutine flc_gs_from_psi(cc_p, lwp )
    ! 
    ! !DESCRIPTION: 
    !
    ! !USES:
    !
    ! !ARGUMENTS
    type(ed_cohort_type) , intent(inout), target  :: cc_p ! current cohort pointer
    real(r8)             , intent(in)             :: lwp  !leaf water potential (MPa)
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: cCohort
    integer          :: FT
    !----------------------------------------------------------------------

    cCohort => cc_p
    FT       = cCohort%pft
    ccohort%co_hydr%btran(:) = &
          (1._r8 + (lwp/EDPftvarcon_inst%hydr_p50_gs(ft))**EDPftvarcon_inst%hydr_avuln_gs(ft))**(-1._r8)
    
  end subroutine flc_gs_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine dflcgsdpsi_from_psi(lwp, ft, dflcgsdpsi)
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)             , intent(in)             :: lwp         ! leaf water potential (MPa)
    integer              , intent(in)             :: ft          ! leaf pft
    real(r8)             , intent(out)            :: dflcgsdpsi  ! fractional loss of conductivity  [-]

    !----------------------------------------------------------------------
  
    associate(& 
         avuln_gs => EDPftvarcon_inst%hydr_avuln_gs, &  ! Input: [real(r8) (:) ] stomatal PLC curve: shape parameter                        [-]
         p50_gs   => EDPftvarcon_inst%hydr_p50_gs     & ! Input: [real(r8) (:) ] stomatal PLC curve: water potential at 50% loss of gs,max  [Pa]
         )

    dflcgsdpsi = -1._r8 * (1._r8 + (lwp/p50_gs(FT))**avuln_gs(FT))**(-2._r8) * &
                          avuln_gs(FT)/p50_gs(FT)*(lwp/p50_gs(FT))**(avuln_gs(FT)-1._r8)
	     
    end associate

  end subroutine dflcgsdpsi_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine flc_from_psi(ft, pm, psi_node, flc_node, site_hydr, bc_in )
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity
    
    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: psi_node    ! water potential                  [MPa]
    real(r8)         , intent(out)    :: flc_node    ! fractional loss of conductivity  [-]
    type(ed_site_hydr_type),optional, intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(bc_in_type),optional, intent(in)             :: bc_in       ! FATES boundary conditions

    !
    ! !LOCAL VARIABLES:

    !----------------------------------------------------------------------
  
    associate(& 
         avuln    => EDPftvarcon_inst%hydr_avuln_node , & ! Input: [real(r8) (:,:) ] PLC curve: vulnerability curve shape parameter          [-]
         p50      => EDPftvarcon_inst%hydr_p50_node     & ! Input: [real(r8) (:,:) ] PLC curve: water potential at 50% loss of conductivity  [Pa]
         )
    
    if(pm <= 4) then
       flc_node = 1._r8/(1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))
    else
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__))
!          call unsatkVG_flc_from_psi(psi_node, &
!             site_hydr%alpha_VG(1), &
!	     site_hydr%n_VG(1),     &
!             site_hydr%m_VG(1),     &
!             site_hydr%l_VG(1),     &
!             flc_node)
       case (campbell)
          call unsatkCampbell_flc_from_psi(psi_node, &
             (-1._r8)*bc_in%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &      ! mm * 1e-3 m/mm * 1e3 kg/m3 * 9.8 m/s2 * 1e-6 MPa/Pa = MPa
	     bc_in%bsw_sisl(1),     &
             flc_node)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if
	     
    end associate

  end subroutine flc_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine dflcdpsi_from_psi(ft, pm, psi_node, dflcdpsi_node, site_hydr, bc_in )
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer          , intent(in)     :: ft             ! PFT index
    integer          , intent(in)     :: pm             ! porous media index
    real(r8)         , intent(in)     :: psi_node       ! water potential                  [MPa]
    real(r8)         , intent(out)    :: dflcdpsi_node  ! fractional loss of conductivity  [-]
    type(ed_site_hydr_type),optional, intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(bc_in_type),optional, intent(in)             :: bc_in       ! FATES boundary conditions

    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         avuln    => EDPftvarcon_inst%hydr_avuln_node, & ! Input: [real(r8) (:,:) ] PLC curve: vulnerability curve shape parameter          [-]
         p50      => EDPftvarcon_inst%hydr_p50_node    & ! Input: [real(r8) (:,:) ] PLC curve: water potential at 50% loss of conductivity  [Pa]
         )
    
    if(pm <= 4) then
       dflcdpsi_node = -1._r8 * (1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))**(-2._r8) * &
                                avuln(ft,pm)/p50(ft,pm)*(psi_node/p50(ft,pm))**(avuln(ft,pm)-1._r8)
    else
       select case (iswc)
       case (van_genuchten)    
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__))
          !call unsatkVG_dflcdpsi_from_psi(psi_node, &
          !      site_hydr%alpha_VG(1), &
          !      site_hydr%n_VG(1),     &
          !      site_hydr%m_VG(1),     &
          !      site_hydr%l_VG(1),     &
          !   dflcdpsi_node)
       case (campbell)
          call unsatkCampbell_dflcdpsi_from_psi(psi_node, &
             (-1._r8)*bc_in%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &
	     bc_in%bsw_sisl(1),     &
             dflcdpsi_node)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if
	     
    end associate

  end subroutine dflcdpsi_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine th_from_psi(ft, pm, psi_node, th_node, site_hydr, bc_in)
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to volumetric water contents
    ! !ARGUMENTS
    integer          , intent(in)            :: ft          ! PFT index
    integer          , intent(in)            :: pm          ! porous media index
    real(r8)         , intent(in)            :: psi_node    ! water potential   [MPa]
    real(r8)         , intent(out)           :: th_node     ! water content     [m3 m-3]
    type(ed_site_hydr_type), intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(bc_in_type), intent(in)             :: bc_in       ! FATES boundary conditions

    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                ! lower bound of initial estimate         [m3 m-3]
    real(r8) :: upper                ! upper bound of initial estimate         [m3 m-3]
    real(r8) :: xtol                 ! error tolerance for x-variable          [m3 m-3]
    real(r8) :: ytol                 ! error tolerance for y-variable          [MPa]
    real(r8) :: satfrac              ! soil saturation fraction                [0-1]
    real(r8) :: psi_check
                                              
    !----------------------------------------------------------------------
  
    associate(& 
         thetas   => EDPftvarcon_inst%hydr_thetas_node  , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content
         resid    => EDPftvarcon_inst%hydr_resid_node     & ! Input: [real(r8) (:,:) ] P-V curve: residual water fraction
         )
   
    if(pm <= 4) then

       lower  = thetas(ft,pm)*(resid(ft,pm) + 0.0001_r8)/cap_corr(pm)
       upper  = thetas(ft,pm)
       xtol   = 1.e-16_r8
       ytol   = 1.e-8_r8
       call bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
       call psi_from_th(ft, pm, th_node, psi_check )

       if(psi_check > -1.e-8_r8) then
          write(fates_log(),*)'bisect_pv returned positive value for water potential at pm = ', char(pm)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       
     

    else
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__))
!          call swcVG_satfrac_from_psi(psi_node, &
!                  site_hydr%alpha_VG(1), &
!                  site_hydr%n_VG(1),     &
!                  site_hydr%m_VG(1),     &
!                  site_hydr%l_VG(1),     &
!                  satfrac)
!          call swcVG_th_from_satfrac(satfrac, &
!                  bc_in%watsat_sisl(1),   &
!                  bc_in%watres_sisl(1),   &
!                  th_node)
       case (campbell) 

          call swcCampbell_satfrac_from_psi(psi_node, &
                  (-1._r8)*bc_in%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &
                  bc_in%bsw_sisl(1),     &
                  satfrac)
          call swcCampbell_th_from_satfrac(satfrac, &
                  bc_in%watsat_sisl(1),   &
                  th_node)
       case default
          write(fates_log(),*)  'invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if
	     
    end associate

  end subroutine th_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
    ! 
    ! !DESCRIPTION: Bisection routine for getting the inverse of the plant PV curve.
    !  An analytical solution is not possible because quadratic smoothing functions
    !  are used to remove discontinuities in the PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(inout)  :: lower       ! lower bound of estimate           [m3 m-3]
    real(r8)      , intent(inout)  :: upper       ! upper bound of estimate           [m3 m-3]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x-variable    [m3 m-3]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for y-variable    [MPa]
    real(r8)      , intent(in)     :: psi_node    ! water potential                   [MPa]
    real(r8)      , intent(out)    :: th_node     ! water content                     [m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x_new                  ! new estimate for x in bisection routine
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    integer  :: nitr                   ! number of iterations 
    !----------------------------------------------------------------------
    if(psi_node > 0.0_r8) then
      write(fates_log(),*)'Error: psi_note become positive,&
                           psi_node=',psi_node
      call endrun(msg=errMsg(sourcefile, __LINE__))  
    endif
    call psi_from_th(ft, pm, lower, y_lo)
    call psi_from_th(ft, pm, upper, y_hi)
    f_lo  = y_lo - psi_node
    f_hi  = y_hi - psi_node
    chg   = upper - lower
    nitr = 0
    do while(abs(chg) .gt. xtol .and. nitr < 100)
       x_new = 0.5_r8*(lower + upper)
       call psi_from_th(ft, pm, x_new, y_new)
       f_new = y_new - psi_node
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
       nitr = nitr + 1
    end do
    
    if(nitr .eq. 100)then
        write(fates_log(),*)'Warning: number of iteraction reaches 100 for bisect_pv'
    endif
    
    th_node = x_new
	     
  end subroutine bisect_pv
  
  !-------------------------------------------------------------------------------!
  subroutine psi_from_th(ft, pm, th_node, psi_node, site_hydr, bc_in)
    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content     [m3 m-3]
    real(r8)         , intent(out)    :: psi_node    ! water potential   [MPa]
    type(ed_site_hydr_type), optional, intent(inout),target :: site_hydr        ! ED site_hydr structure
    type(bc_in_type), optional, intent(in)             :: bc_in       ! FATES boundary conditions
    !
    ! !LOCAL VARIABLES:
    real(r8) :: satfrac                  ! saturation fraction [0-1]
    !----------------------------------------------------------------------
  
    if(pm <= 4) then       ! plant

       call tq2(ft, pm, th_node*cap_corr(pm), psi_node)

    else if(pm == 5) then  ! soil

!! NOTE. FIX: The below sidesteps the problem of averaging potentially variable soil hydraulic properties with depth
!!        and simply assigns the bulk soil (bucket) approximation of hydraulic properties as equal to the top soil layer.
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
!          call swcVG_psi_from_th(th_node, &
!                  bc_in%watsat_sisl(1),   &
!                  bc_in%watres_sisl(1),   &
!                  site_hydr%alpha_VG(1), &
!                  site_hydr%n_VG(1),     &
!                  site_hydr%m_VG(1),     &
!                  site_hydr%l_VG(1),     &
!                  psi_node)
       case (campbell)
          call swcCampbell_psi_from_th(th_node, &
                  bc_in%watsat_sisl(1),   &
                  (-1._r8)*bc_in%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &
	          bc_in%bsw_sisl(1),     &
                  psi_node)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

    end if
	     
  end subroutine psi_from_th
  
  !-------------------------------------------------------------------------------!
  subroutine dpsidth_from_th(ft, pm, th_node, y, site_hydr, bc_in) 
    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content                            [m3 m-3]
    real(r8)         , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    type(ed_site_hydr_type), optional,intent(inout) :: site_hydr
    type(bc_in_type), optional,intent(in) :: bc_in
    !
    ! !LOCAL VARIABLES:

    real(r8) :: satfrac                  ! saturation fraction [0-1]
    !----------------------------------------------------------------------
  
    if(pm <= 4) then       ! plant
       call dtq2dth(ft, pm, th_node*cap_corr(pm), y)
    else if(pm == 5) then  ! soil
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
          !call swcVG_dpsidth_from_th(th_node, &
          !        bc_in%watsat_sisl(1),   &
          !        bc_in%watres_sisl(1),   &
          !        site_hydr%alpha_VG(1), &
          !        site_hydr%n_VG(1),     &
          !        site_hydr%m_VG(1),     &
          !        site_hydr%l_VG(1),     &
          !        y)
       case (campbell)
          call swcCampbell_dpsidth_from_th(th_node, &
                  bc_in%watsat_sisl(1),   &
                  (-1._r8)*bc_in%sucsat_sisl(1)*denh2o*grav*1.e-9_r8, &
	          bc_in%bsw_sisl(1),     &
                  y)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if
	     
  end subroutine dpsidth_from_th
  
  !-------------------------------------------------------------------------------!
  subroutine tq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    y = (-y_bq2 + sqrt(y_bq2*y_bq2 - 4._r8*beta2*y_cq2))/(2*beta2)

  end subroutine tq2
  
  !-------------------------------------------------------------------------------!
  subroutine dtq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: dydth_bq2              ! returned derivative from dbq2dth()
    real(r8) :: dydth_cq2              ! returned derivative from dcq2dth()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    call dbq2dth(ft, pm, x, dydth_bq2)
    call dcq2dth(ft, pm, x, dydth_cq2)
    y = 1._r8/(2._r8*beta2)*(-dydth_bq2 + 0.5_r8*((y_bq2*y_bq2 - 4._r8*beta2*y_cq2)**(-0.5_r8)) * &
                                                 (2._r8*y_bq2*dydth_bq2 - 4._r8*beta2*dydth_cq2))

  end subroutine dtq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine bq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = -1._r8*(y_tq1 + y_cavitation)
	     
  end subroutine bq2
  
  !-------------------------------------------------------------------------------!
  subroutine dbq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = -1._r8*(dydth_tq1 + dcavdth)

  end subroutine dbq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine cq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = y_tq1*y_cavitation
	     
  end subroutine cq2
  
  !-------------------------------------------------------------------------------!
  subroutine dcq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq2() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = y_tq1*dcavdth + dydth_tq1*y_cavitation
	     
  end subroutine dcq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine tq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: either calls the elastic region of the PV curve (leaves) or
    ! does a smoothing function for capillary-to-elastic region of the plant PV
    ! curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: y_elastic              ! returned y (psi) value from elasticPV()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call elasticPV(ft, pm, x, y_elastic)
       y = y_elastic
    else if(pm <= 4) then       ! sapwood has a capillary region
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       y = (-y_bq1 - sqrt(y_bq1*y_bq1 - 4._r8*beta1*y_cq1))/(2*beta1)
    end if !porous media

  end subroutine tq1
  
  !-------------------------------------------------------------------------------!
  subroutine dtq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of tq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: dydth_bq1              ! returned derivative from dbq1dth()
    real(r8) :: dydth_cq1              ! returned derivative from dcq1dth()
    real(r8) :: delasticdth            ! returned derivative from delasticPVdth()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call delasticPVdth(ft, pm, x, delasticdth)
       y = delasticdth
    else if(pm <= 4) then       ! sapwood has a capillary region
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       call dbq1dth(ft, pm, x, dydth_bq1)
       call dcq1dth(ft, pm, x, dydth_cq1)
       y = 1._r8/(2._r8*beta1)*(-dydth_bq1 - 0.5_r8*((y_bq1*y_bq1 - 4._r8*beta1*y_cq1)**(-0.5_r8)) * &
                                                    (2._r8*y_bq1*dydth_bq1 - 4._r8*beta1*dydth_cq1))
    end if

  end subroutine dtq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine bq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = -1._r8*(y_capillary + y_elastic)
	     
  end subroutine bq1
  
  !-------------------------------------------------------------------------------!
  subroutine dbq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of bq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = -1._r8*(delasticdth + dcapdth)
	     
  end subroutine dbq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine cq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = y_capillary*y_elastic
	     
  end subroutine cq1
  
  !-------------------------------------------------------------------------------!
  subroutine dcq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = y_elastic*dcapdth + delasticdth*y_capillary
	     
  end subroutine dcq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine cavitationPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    !----------------------------------------------------------------------
  
    call solutepsi(ft, pm, x, y_solute)
    y = y_solute
	     
  end subroutine cavitationPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcavitationPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cavitationPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    y = dsoldth
	     
  end subroutine dcavitationPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine elasticPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    real(r8) :: y_pressure         ! returned y (psi) value from pressurepsi()
    !----------------------------------------------------------------------
  
    call solutepsi(ft, pm, x, y_solute)
    call pressurepsi(ft, pm, x, y_pressure)
    y = y_solute + y_pressure
	     
  end subroutine elasticPV
  
  !-------------------------------------------------------------------------------!
  subroutine delasticPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of elasticPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    real(r8) :: dpressdth         ! returned derivative from dpressurepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    call dpressurepsidth(ft, pm, x, dpressdth)
    y = dsoldth + dpressdth
	     
  end subroutine delasticPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine solutepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes solute water potential (negative) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => EDPftvarcon_inst%hydr_pinot_node,  & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => EDPftvarcon_inst%hydr_thetas_node, & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node   & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
    
      y = pinot(ft,pm)*thetas(ft,pm)*(rwcft(pm) - resid(ft,pm)) / &
            (x - thetas(ft,pm)*resid(ft,pm))
	     
    end associate

  end subroutine solutepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dsolutepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of solutepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => EDPftvarcon_inst%hydr_pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => EDPftvarcon_inst%hydr_thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node       & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
    
      y = -1._r8*thetas(ft,pm)*pinot(ft,pm)*(rwcft(pm) - resid(ft,pm)) / &
            ((x - thetas(ft,pm)*resid(ft,pm))**2._r8)
	     
    end associate

  end subroutine dsolutepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine pressurepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes pressure water potential (positive) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => EDPftvarcon_inst%hydr_pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => EDPftvarcon_inst%hydr_thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node     , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => EDPftvarcon_inst%hydr_epsil_node       & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
      
      y = epsil(ft,pm) * (x - thetas(ft,pm)*rwcft(pm)) / &
            (thetas(ft,pm)*(rwcft(pm)-resid(ft,pm))) - pinot(ft,pm)
	     
    end associate

  end subroutine pressurepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dpressurepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of pressurepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         thetas  => EDPftvarcon_inst%hydr_thetas_node, & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => EDPftvarcon_inst%hydr_epsil_node   & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
      
      y = epsil(ft,pm)/(thetas(ft,pm)*(rwcft(pm) - resid(ft,pm)))
	     
    end associate

  end subroutine dpressurepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine capillaryPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the capillary region of the plant
    !  PV curve (sapwood only)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         thetas    => EDPftvarcon_inst%hydr_thetas_node     & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
   
      y = cap_int(pm) + cap_slp(pm)/thetas(ft,pm)*x
	     
    end associate

  end subroutine capillaryPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcapillaryPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of capillaryPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         thetas    => EDPftvarcon_inst%hydr_thetas_node    & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
   
    y = cap_slp(pm)/thetas(ft,pm)
	     
    end associate

  end subroutine dcapillaryPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine swcVG_satfrac_from_th(th, watsat, watres, satfrac)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns saturation fraction given water content, porosity, and residual water content
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: th       !soil volumetric water content       [m3 m-3]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
  real(r8), intent(out)           :: satfrac  !saturation fraction                 [0-1]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  satfrac = (th - watres)/(watsat - watres)

  end subroutine swcVG_satfrac_from_th

  !-------------------------------------------------------------------------------!
  subroutine swcCampbell_satfrac_from_th(th, watsat, satfrac)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns saturation fraction given water content and porosity
  !
  !USES
  !ARGUMENTS:
  real(r8), intent(in)            :: th       !soil volumetric water content       [m3 m-3]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(out)           :: satfrac  !saturation fraction                 [0-1]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  satfrac = th/watsat

  end subroutine swcCampbell_satfrac_from_th

  !-------------------------------------------------------------------------------!
  subroutine swcVG_psi_from_th(th, watsat, watres, alpha, n, m, l, psi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns water potential given water content and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: th       !volumetric water content       [m3 m-3]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
  real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
  real(r8), intent(in)            :: n        !pore-size distribution index   [-]
  real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
  real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
  real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
  !
  ! !LOCAL VARIABLES:
  real(r8)                        :: satfrac  !saturation fraction            [0-1]
  !------------------------------------------------------------------------------

  call swcVG_satfrac_from_th(th, watsat, watres, satfrac)
  call swcVG_psi_from_satfrac(satfrac, alpha, n, m, l, psi)

  end subroutine swcVG_psi_from_th

  !-------------------------------------------------------------------------------!
  subroutine swcCampbell_psi_from_th(th, watsat, psisat, B, psi)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns water potential given saturation fraction, air-entry pressure and shape parameter
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: th       !volumetric water content       [m3 m-3]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)            :: psisat   !air entry pressure             [MPa]
  real(r8), intent(in)            :: B        !shape parameter                [-]
  real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
  !
  ! !LOCAL VARIABLES:
  real(r8)                        :: satfrac  !saturation fraction            [0-1]
  !------------------------------------------------------------------------------

  call swcCampbell_satfrac_from_th(th, watsat, satfrac)
  call swcCampbell_psi_from_satfrac(satfrac, psisat, B, psi)

  end subroutine swcCampbell_psi_from_th

  !-------------------------------------------------------------------------------!
  subroutine swcVG_psi_from_satfrac(satfrac, alpha, n, m, l, psi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns water potential given saturation fraction and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: satfrac  !saturation fraction            [0-1]
  real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
  real(r8), intent(in)            :: n        !pore-size distribution index   [-]
  real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
  real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
  real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  psi = -1._r8/alpha*(satfrac**(-1._r8/m)-1._r8)**(1._r8/n)

  end subroutine swcVG_psi_from_satfrac

  !-------------------------------------------------------------------------------!
  subroutine swcCampbell_psi_from_satfrac(satfrac, psisat, B, psi)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns water potential given saturation fraction, air-entry pressure and shape parameter
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: satfrac  !saturation fraction            [0-1]
  real(r8), intent(in)            :: psisat   !air entry pressure             [MPa]
  real(r8), intent(in)            :: B        !shape parameter                [-]
  real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  psi = psisat*(satfrac**(-B))

  end subroutine swcCampbell_psi_from_satfrac

  !-------------------------------------------------------------------------------!
  subroutine swcVG_th_from_satfrac(satfrac, watsat, watres, th)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns water content given saturation fraction, porosity and residual water content
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: satfrac  !saturation fraction                 [0-1]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
  real(r8), intent(out)           :: th       !soil volumetric water content       [m3 m-3]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  th = watres + satfrac*(watsat - watres)

  end subroutine swcVG_th_from_satfrac

  !-------------------------------------------------------------------------------!
  subroutine swcCampbell_th_from_satfrac(satfrac, watsat, th)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns water content given saturation fraction and porosity
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: satfrac  !saturation fraction                 [0-1]
  real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(out)           :: th       !soil volumetric water content       [m3 m-3]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  th = satfrac*watsat

  end subroutine swcCampbell_th_from_satfrac

  !-----------------------------------------------------------------------
  subroutine swcVG_satfrac_from_psi(psi, alpha, n, m, l, satfrac)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns saturation fraction given water potential and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: psi      !soil matric potential          [MPa]
  real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
  real(r8), intent(in)            :: n        !pore-size distribution index   [-]
  real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
  real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
  real(r8), intent(out)           :: satfrac  !soil saturation fraction       [0-1]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  satfrac = (1._r8/(1._r8 + (alpha*abs(psi))**n))**m

  end subroutine swcVG_satfrac_from_psi

  !-----------------------------------------------------------------------
  subroutine swcCampbell_satfrac_from_psi(psi, psisat, B, satfrac)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns saturation fraction given water potential and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: psi      !soil matric potential          [MPa]
  real(r8), intent(in)            :: psisat   !air-entry pressure             [MPa]
  real(r8), intent(in)            :: B        !shape parameter                [-]
  real(r8), intent(out)           :: satfrac  !soil saturation fraction       [0-1]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  satfrac = (psi/psisat)**(-1.0_r8/B)

  end subroutine swcCampbell_satfrac_from_psi

  !-----------------------------------------------------------------------
  subroutine swcVG_dpsidth_from_th(th, watsat, watres, alpha, n, m, l, dpsidth)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given water content and SWC parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: th       !volumetric water content                       [m3 m-3]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: satfrac  !saturation fraction            [0-1]
  !------------------------------------------------------------------------------

  call swcVG_satfrac_from_th(th, watsat, watres, satfrac)
  call swcVG_dpsidth_from_satfrac(satfrac, watsat, watres, alpha, n, m, l, dpsidth)

  end subroutine swcVG_dpsidth_from_th

  !-----------------------------------------------------------------------
  subroutine swcCampbell_dpsidth_from_th(th, watsat, psisat, B, dpsidth)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given water content and SWC parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: th       !volumetric water content                       [m3 m-3]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: psisat   !air entry pressure                             [MPa]
  real(r8), intent(in)  :: B        !shape parameter                                [-]
  real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: satfrac  !saturation fraction            [0-1]
  !------------------------------------------------------------------------------

  call swcCampbell_satfrac_from_th(th, watsat, satfrac)
  call swcCampbell_dpsidth_from_satfrac(satfrac, watsat, psisat, B, dpsidth)

  end subroutine swcCampbell_dpsidth_from_th

  !-----------------------------------------------------------------------
  subroutine swcVG_dpsidth_from_satfrac(satfrac, watsat, watres, alpha, n, m, l, dpsidth)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: satfrac  !saturation fraction                            [0-1]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp0    !temporary
  real(r8)              :: temp1    !temporary
  real(r8)              :: temp2    !temporary
  real(r8)              :: temp3    !temporary
  !------------------------------------------------------------------------------

  temp0   = 1._r8/(m*n*alpha*(watsat-watres))
  temp1   = satfrac**(-1._r8/m) - 1._r8
  temp2   = temp1**(1._r8/n - 1._r8)
  temp3   = satfrac**(-1._r8/m - 1._r8)
  dpsidth = temp0*temp2*temp3

  end subroutine swcVG_dpsidth_from_satfrac

  !-----------------------------------------------------------------------
  subroutine swcCampbell_dpsidth_from_satfrac(satfrac, watsat, psisat, B, dpsidth)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: satfrac  !saturation fraction                            [0-1]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
  real(r8), intent(in)  :: B        !shape parameter                                [-]
  real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
  !------------------------------------------------------------------------------

  dpsidth = psisat*(-B)/watsat*(satfrac)**(-B-1._r8)

  end subroutine swcCampbell_dpsidth_from_satfrac

  !-----------------------------------------------------------------------
  subroutine unsatkVG_flc_from_psi(psi, alpha, n, m, l, flc)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns unsaturated hydraulic conductivity 
  ! given water potential and SWC parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: flc      !k/ksat ('fractional loss of conductivity')     [-]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp          !temporary
  real(r8)              :: fac1a         !temporary
  real(r8)              :: fac1b         !temporary
  real(r8)              :: fac1          !temporary
  real(r8)              :: fac2          !temporary
  !------------------------------------------------------------------------------

  temp       = ( alpha*abs(psi)      ) ** (n)
  fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
  fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
  fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
  fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)
  
  flc        =   fac1 * fac2

  end subroutine unsatkVG_flc_from_psi

  !-----------------------------------------------------------------------
  subroutine unsatkCampbell_flc_from_psi(psi, psisat, B, flc)
  !
  ! DESCRIPTION
  ! Campbell (1974) soil water characteristic (retention) curve
  ! returns unsaturated hydraulic conductivity 
  ! given water potential and SWC parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
  real(r8), intent(in)  :: B        !shape parameter                                [-]
  real(r8), intent(out) :: flc      !k/ksat ('fractional loss of conductivity')     [-]
  !------------------------------------------------------------------------------
  
  flc        =  (psi/psisat)**(-2._r8-3._r8/B)

  end subroutine unsatkCampbell_flc_from_psi

  !-----------------------------------------------------------------------
  subroutine unsatkVG_dflcdpsi_from_psi(psi, alpha, n, m, l, dflcdpsi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: dflcdpsi !derivative of k/ksat (flc) wrt psi             [MPa-1]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp          !temporary
  real(r8)              :: fac1a         !temporary
  real(r8)              :: fac1b         !temporary
  real(r8)              :: fac1          !temporary
  real(r8)              :: fac2          !temporary
  real(r8)              :: dtemp         !temporary
  real(r8)              :: dfac1adpsi    !temporary
  real(r8)              :: dfac1bdpsi    !temporary
  real(r8)              :: dfac1dpsi     !temporary
  real(r8)              :: dfac2dpsi     !temporary
  !------------------------------------------------------------------------------

  temp       = ( alpha*abs(psi)      ) ** (n)
  fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
  fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
  fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
  fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)
  
  dtemp      =   n * alpha * ( alpha*abs(psi) ) ** (n-1._r8)
  dfac1adpsi = ( n-1._r8 ) * alpha * ( alpha*abs(psi) ) ** (n-2._r8)
  dfac1bdpsi = ( -1._r8 ) * m * dtemp * ( 1._r8 + temp ) ** (-1._r8*m - 1._r8)
  dfac1dpsi  = ( 2._r8 ) * ( 1._r8 - fac1a*fac1b ) * ( -1._r8*dfac1bdpsi*fac1a - dfac1adpsi*fac1b )
  dfac2dpsi  = ( -0.5_r8 ) * m * dtemp * (1._r8 + temp)**(-0.5_r8*m-1._r8)
  
  dflcdpsi   = ( -1._r8 ) * ( dfac2dpsi*fac1 + dfac1dpsi*fac2 )    ! BOC... mult by -1 because unsatk eqn is based on abs(psi)

  end subroutine unsatkVG_dflcdpsi_from_psi

  !-----------------------------------------------------------------------
  subroutine unsatkCampbell_dflcdpsi_from_psi(psi, psisat, B, dflcdpsi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
  real(r8), intent(in)  :: B        !shape parameter                                [-]
  real(r8), intent(out) :: dflcdpsi !derivative of k/ksat (flc) wrt psi             [MPa-1]
  !------------------------------------------------------------------------------

  dflcdpsi   = psisat*(-2._r8-3._r8/B)*(psi/psisat)**(-3._r8-3._r8/B)

  end subroutine unsatkCampbell_dflcdpsi_from_psi


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
  
  function zeng2001_crootfr(a, b, z) result(crootfr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: a,b    ! pft parameters
    real(r8) , intent(in) :: z      ! soil depth (m)
    !
    ! !RESULT
    real(r8) :: crootfr                         ! cumulative root fraction
    !
    !------------------------------------------------------------------------
    crootfr      = 1._r8 - .5_r8*(exp(-a*z) + exp(-b*z))

    return

  end function zeng2001_crootfr

  ! =====================================================================================

  subroutine shellGeom(l_aroot, rs1, area, dz, r_out_shell, r_node_shell, v_shell)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil surrounding fine roots remains
    ! the same.  
    !
    ! !USES:
    
    !
    ! !ARGUMENTS:
    real(r8)     , intent(in)             :: l_aroot
    real(r8)     , intent(in)             :: rs1
    real(r8)     , intent(in)             :: area
    real(r8)     , intent(in)             :: dz
    real(r8)     , intent(out)            :: r_out_shell(:)
    real(r8)     , intent(out)            :: r_node_shell(:)
    real(r8)     , intent(out)            :: v_shell(:)                 ! volume of a single rhizosphere shell
    !
    ! !LOCAL VARIABLES:
    integer                        :: k                                 ! rhizosphere shell indicies
    !-----------------------------------------------------------------------

    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    r_out_shell(nshell) = (pi_const*l_aroot/(area*dz))**(-0.5_r8)                  ! eqn(8) S98
    if(nshell > 1) then
       do k = 1,nshell-1
          r_out_shell(k)   = rs1*(r_out_shell(nshell)/rs1)**((real(k,r8))/real(nshell,r8))  ! eqn(7) S98
       enddo
    end if

    ! set nodal (midpoint) radii of these shells
    do k = 1,nshell
       if(k == 1) then
          ! BOC...not doing this as it requires PFT-specific fine root thickness, but this is at column level
          ! r_node_shell(k) = 0.5_r8*(rs1 + r_out_shell(k))
          r_node_shell(k) = 0.5_r8*(r_out_shell(k))
       else
          r_node_shell(k) = 0.5_r8*(r_out_shell(k-1) + r_out_shell(k))
       end if
    enddo

    ! update volumes
    do k = 1,nshell
       if(k == 1) then
	  ! BOC...not doing this as it requires PFT-specific fine root thickness but this is at column level
          ! v_shell(k)   = pi*dz*(r_out_shell(k)**2._r8 - rs1**2._r8)
          v_shell(k)     = pi_const*dz*(r_out_shell(k)**2._r8)
       else
          v_shell(k)     = pi_const*dz*(r_out_shell(k)**2._r8 - r_out_shell(k-1)**2._r8)
       end if
    enddo

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
 

end module FatesPlantHydraulicsMod
