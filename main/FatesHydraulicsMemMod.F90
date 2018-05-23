module FatesHydraulicsMemMod

   use FatesConstantsMod, only : r8 => fates_r8
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use FatesConstantsMod, only : itrue,ifalse
   use EDParamsMod      , only : hydr_psi0
   use EDParamsMod      , only : hydr_psicap
   
   implicit none

   ! Number of soil layers for indexing cohort fine root quanitities
   integer, parameter                            :: nlevsoi_hyd = 10

   ! number of distinct types of plant porous media (leaf, stem, troot, aroot)
   integer, parameter                            :: n_porous_media = 4

   integer, parameter                            :: n_hypool_leaf  = 1
   integer, parameter                            :: n_hypool_stem  = 1
   integer, parameter                            :: n_hypool_troot = 1
   integer, parameter                            :: n_hypool_aroot = 1
   integer, parameter                            :: nshell      = 5

   ! number of aboveground plant water storage nodes
   integer, parameter                            :: n_hypool_ag    = n_hypool_leaf+n_hypool_stem

   ! total number of water storage nodes
   integer, parameter                            :: n_hypool_tot   = n_hypool_ag + n_hypool_troot + n_hypool_aroot + nshell

   ! vector indexing the type of porous medium over an arbitrary number of plant pools
   integer, parameter,dimension(n_hypool_tot)       :: porous_media = (/1,2,3,4,5,5,5,5,5/) 

   ! number of previous timestep's leaf water potential to be retained
   integer, parameter                            :: numLWPmem             = 4

   ! mirror of nlevcan, hard-set for simplicity, remove nlevcan_hyd on a rainy day
   ! Note (RGK): uscing nclmax causes annoying circular dependency (this needs EDTypes, EDTypes needs this)
   ! way around that: dynamic allocation, or just keep this, but set the value high
   integer, parameter                            :: nlevcan_hyd = 2                       
                       
   ! Mean fine root radius expected in the bulk soil                
   real(r8), parameter                           :: fine_root_radius_const = 0.001_r8               
   
   ! Constant parameters (for time being, C2B is constant, 
   ! slated for addition to parameter file (RGK 08-2017))
   ! Carbon 2 biomass ratio
   real(r8), parameter                           :: C2B        = 2.0_r8 
              
   ! P-V curve: total RWC @ which elastic drainage begins     [-]        
   real(r8), parameter,dimension(n_porous_media) :: rwcft   = (/1.0_r8,0.958_r8,0.958_r8,0.958_r8/)

   ! P-V curve: total RWC @ which capillary reserves exhausted
   real(r8), parameter,dimension(n_porous_media) :: rwccap  = (/1.0_r8,0.947_r8,0.947_r8,0.947_r8/) 

   ! Derived parameters
   ! ----------------------------------------------------------------------------------------------

   ! P-V curve: slope of capillary region of curve
   real(r8), dimension(n_porous_media)           :: cap_slp                                         

   ! P-V curve: intercept of capillary region of curve
   real(r8), dimension(n_porous_media)           :: cap_int   

   ! P-V curve: correction for nonzero psi0x
   real(r8), dimension(n_porous_media)           :: cap_corr                                        


   type ed_site_hydr_type

      ! Plant Hydraulics
     
     real(r8),allocatable :: v_shell(:,:)           ! Volume of rhizosphere compartment (m) 
     real(r8),allocatable :: v_shell_init(:,:)      ! Previous volume of rhizosphere compartment (m) 
     real(r8),allocatable :: v_shell_1D(:)          ! Volume of rhizosphere compartment (m)
     real(r8),allocatable :: r_node_shell(:,:)      ! Nodal radius of rhizosphere compartment (m)
     real(r8),allocatable :: r_node_shell_init(:,:) ! Previous Nodal radius of rhizosphere compartment (m)
     real(r8),allocatable :: l_aroot_layer(:)       ! Total length (across cohorts) of absorbing
                                                    !  roots by soil layer (m)
     real(r8),allocatable :: l_aroot_layer_init(:)  ! Total length (across cohorts) of absorbing
                                                    !  roots by soil layer (m)
     real(r8),allocatable :: kmax_upper_shell(:,:)  ! Maximum soil hydraulic conductance node k 
                                                    ! to upper (closer to atmosphere) rhiz 
                                                    ! shell boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: kmax_bound_shell(:,:)  ! Maximum soil hydraulic conductance at upper
                                                    ! (closer to atmosphere) rhiz shell 
                                                    ! boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: kmax_lower_shell(:,:)  ! Maximum soil hydraulic conductance node k
                                                    ! to lower (further from atmosphere) 
                                                    ! rhiz shell boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: r_out_shell(:,:)       ! Outer radius of rhizosphere compartment (m)
     real(r8),allocatable :: r_out_shell_1D(:)      ! Outer radius of rhizosphere compartment (m) (USED?)
     real(r8),allocatable :: r_node_shell_1D(:)     ! Nodal radius of rhizosphere compartment (m)

     real(r8),allocatable :: rs1(:)                 ! Mean fine root radius (m) (currently a constant)

     real(r8),allocatable :: kmax_upper_shell_1D(:) ! Maximum soil hydraulic conductance node 
                                                    ! k to upper (closer to atmosphere) rhiz 
                                                    ! shell boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: kmax_bound_shell_1D(:) ! Maximum soil hydraulic conductance at upper 
                                                    ! (closer to atmosphere) rhiz shell
                                                    ! boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: kmax_lower_shell_1D(:) ! Maximum soil hydraulic conductance node 
                                                    ! k to lower (further from atmosphere) rhiz 
                                                    ! shell boundaries (kg s-1 MPa-1)

     integer,allocatable :: supsub_flag(:)          ! index of the outermost rhizosphere shell 
                                                    ! encountering super- or sub-saturation
     real(r8),allocatable :: h2osoi_liqvol_shell(:,:) ! volumetric water in rhizosphere compartment (m3/m3)

     real(r8),allocatable :: h2osoi_liq_prev(:)     ! liquid water mass for the bulk soil layer
                                                    ! defined at the end of the hydraulics sequence
                                                    ! after root water has been extracted.  This should
                                                    ! be equal to the sum of the water over the rhizosphere shells
     
     real(r8),allocatable :: psisoi_liq_innershell(:) ! Matric potential of the inner rhizosphere shell (MPa)
     
     real(r8),allocatable :: recruit_w_uptake(:)    ! recruitment water uptake (kg H2o/m2/s)

    
     real(r8) :: l_aroot_1D                         ! Total (across cohorts) absorbing root 
                                                    ! length across all layers

     real(r8) :: errh2o_hyd                         ! plant hydraulics error summed across 
                                                    ! cohorts to column level (mm)
     real(r8) :: dwat_veg                           ! change in stored water in vegetation
                                                    ! column level (kg)
     real(r8) :: h2oveg                             ! stored water in vegetation (kg/m2)

     real(r8) :: h2oveg_dead                        ! stored water in dead vegetation (kg/m2)

     
     !     Hold Until Van Genuchten is implemented
     ! col inverse of air-entry pressure     [MPa-1]  (for van Genuchten SWC only)
     !     real(r8), allocatable :: alpha_VG(:)  
     ! col pore-size distribution index      [-]      (for van Genuchten SWC only)
     !     real(r8), allocatable :: n_VG(:) 
     ! = 1 - 1/n_VG                          [-]      (for van Genuchten SWC only)   
     !     real(r8), allocatable :: m_VG(:) 
     ! col pore tortuosity parameter         [-]      (for van Genuchten SWC only)    
     !     real(r8), allocatable :: l_VG(:)     
     
  contains
     
     procedure :: InitHydrSite
     
  end type ed_site_hydr_type

  ! This whole structure is actually not used, because netRad_mem() is actually not used
  ! Keeping the code in place in case a patch-level hydraulics variable is desired (RGK 03-2018)

  !type ed_patch_hydr_type
  !   real(r8) ::  netRad_mem(numLWPmem)          ! patch-level net radiation for the previous numLWPmem timesteps [W m-2]
  !end type ed_patch_hydr_type


  type ed_cohort_hydr_type
     
                                                  ! BC...PLANT HYDRAULICS - "constants" that change with size. 
                                                  ! Heights are referenced to soil surface (+ = above; - = below)
     real(r8) ::  z_node_ag(n_hypool_ag)          ! nodal height of aboveground water storage compartments            [m]
     real(r8) ::  z_node_troot(n_hypool_troot)    ! nodal height of belowground water storage compartments            [m]
     real(r8) ::  z_node_aroot(nlevsoi_hyd)       ! nodal height of absorbing root water storage compartments         [m]
     real(r8) ::  z_upper_ag(n_hypool_ag)         ! upper boundary height of aboveground water storage compartments   [m]
     real(r8) ::  z_upper_troot(n_hypool_troot)   ! upper boundary height of belowground water storage compartments   [m]
     real(r8) ::  z_lower_ag(n_hypool_ag)         ! lower boundary height of aboveground water storage compartments   [m]
     real(r8) ::  z_lower_troot(n_hypool_troot)   ! lower boundary height of belowground water storage compartments   [m]
     real(r8) ::  kmax_upper(n_hypool_ag)         ! maximum hydraulic conductance from node to upper boundary         [kg s-1 MPa-1]
     real(r8) ::  kmax_lower(n_hypool_ag)         ! maximum hydraulic conductance from node to lower boundary         [kg s-1 MPa-1]
     real(r8) ::  kmax_upper_troot                ! maximum hydraulic conductance from troot node to upper boundary   [kg s-1 MPa-1]
     real(r8) ::  kmax_bound(n_hypool_ag)         ! maximum hydraulic conductance at lower boundary (canopy to troot) [kg s-1 MPa-1]
     real(r8) ::  kmax_treebg_tot                 ! total belowground tree kmax (troot to surface of absorbing roots) [kg s-1 MPa-1]
     real(r8) ::  kmax_treebg_layer(nlevsoi_hyd)  ! total belowground tree kmax partitioned by soil layer             [kg s-1 MPa-1]
     real(r8) ::  v_ag_init(n_hypool_ag)          ! previous day's volume of aboveground water storage compartments   [m3]
     real(r8) ::  v_ag(n_hypool_ag)               ! volume of aboveground water storage compartments                  [m3]
     real(r8) ::  v_troot_init(n_hypool_troot)    ! previous day's volume of belowground water storage compartments   [m3]
     real(r8) ::  v_troot(n_hypool_troot)         ! volume of belowground water storage compartments                  [m3]
     real(r8) ::  v_aroot_tot                     ! total volume of absorbing roots                                   [m3]
     real(r8) ::  v_aroot_layer_init(nlevsoi_hyd) ! previous day's volume of absorbing roots by soil layer            [m3]
     real(r8) ::  v_aroot_layer(nlevsoi_hyd)      ! volume of absorbing roots by soil layer                           [m3]
     real(r8) ::  l_aroot_tot                     ! total length of absorbing roots                                   [m]
     real(r8) ::  l_aroot_layer(nlevsoi_hyd)      ! length of absorbing roots by soil layer                           [m]

                                                  ! BC PLANT HYDRAULICS - state variables
     real(r8) ::  th_ag(n_hypool_ag)              ! water in aboveground compartments                                 [kgh2o/indiv]
     real(r8) ::  th_troot(n_hypool_troot)        ! water in belowground compartments                                 [kgh2o/indiv]
     real(r8) ::  th_aroot(nlevsoi_hyd)           ! water in absorbing roots                                          [kgh2o/indiv]
     real(r8) ::  lwp_mem(numLWPmem)              ! leaf water potential over the previous numLWPmem timesteps        [MPa]
     real(r8) ::  lwp_stable                      ! leaf water potential just before it became unstable               [MPa]
     logical  ::  lwp_is_unstable                 ! flag for instability of leaf water potential over previous timesteps
     real(r8) ::  psi_ag(n_hypool_ag)             ! water potential in aboveground compartments                       [MPa]
     real(r8) ::  psi_troot(n_hypool_troot)       ! water potential in belowground compartments                       [MPa]
     real(r8) ::  psi_aroot(nlevsoi_hyd)          ! water potential in absorbing roots                                [MPa]
     real(r8) ::  flc_ag(n_hypool_ag)             ! fractional loss of conductivity in aboveground compartments       [-]
     real(r8) ::  flc_troot(n_hypool_troot)       ! fractional loss of conductivity in belowground compartments       [-]
     real(r8) ::  flc_aroot(nlevsoi_hyd)          ! fractional loss of conductivity in absorbing roots                [-]
     real(r8) ::  flc_min_ag(n_hypool_ag)         ! min attained fractional loss of conductivity in 
                                                  ! aboveground compartments (for tracking xylem refilling dynamics) [-]
     real(r8) ::  flc_min_troot(n_hypool_troot)   ! min attained fractional loss of conductivity in 
                                                  ! belowground compartments (for tracking xylem refilling dynamics) [-]
     real(r8) ::  flc_min_aroot(nlevsoi_hyd)      ! min attained fractional loss of conductivity in absorbing roots 
                                                  ! (for tracking xylem refilling dynamics)          [-]
     real(r8) ::  refill_thresh                   ! water potential threshold for xylem refilling to occur            [MPa]
     real(r8) ::  refill_days                     ! number of days required for 50% of xylem refilling to occur       [days]
     real(r8) ::  btran(nlevcan_hyd)              ! leaf water potential limitation on gs                             [0-1]
     real(r8) ::  supsub_flag                     ! k index of last node to encounter supersaturation or 
                                                  ! sub-residual water content  (+ supersaturation; - subsaturation)
     real(r8) ::  iterh1                          ! number of iterations required to achieve tolerable water balance error
     real(r8) ::  iterh2                          ! number of inner iterations
     real(r8) ::  errh2o                          ! total water balance error per unit crown area                     [kgh2o/m2]
     
                                                  ! BC PLANT HYDRAULICS - fluxes
     real(r8) ::  qtop_dt                         ! transpiration boundary condition (+ to atm)                       [kg/indiv/timestep]
     real(r8) ::  dqtopdth_dthdt                  ! transpiration tendency term (+ to atm)                            [kg/indiv/timestep]
                                                  ! NOTE: total transpiration is given by qtop_dt + dqtopdth_dthdt
     real(r8) ::  sapflow                         ! flow at base of tree (+ upward)                                   [kg/indiv/timestep]
     real(r8) ::  rootuptake                      ! net flow into roots (+ into roots)                                [kg/indiv/timestep]
     real(r8) ::  rootuptake01                    ! net flow into roots (+ into roots), soil layer 1                  [kg/indiv/timestep]
     real(r8) ::  rootuptake02                    ! net flow into roots (+ into roots), soil layer 2                  [kg/indiv/timestep]
     real(r8) ::  rootuptake03                    ! net flow into roots (+ into roots), soil layer 3                  [kg/indiv/timestep]
     real(r8) ::  rootuptake04                    ! net flow into roots (+ into roots), soil layer 4                  [kg/indiv/timestep]
     real(r8) ::  rootuptake05                    ! net flow into roots (+ into roots), soil layer 5                  [kg/indiv/timestep]
     real(r8) ::  rootuptake06                    ! net flow into roots (+ into roots), soil layer 6                  [kg/indiv/timestep]
     real(r8) ::  rootuptake07                    ! net flow into roots (+ into roots), soil layer 7                  [kg/indiv/timestep]
     real(r8) ::  rootuptake08                    ! net flow into roots (+ into roots), soil layer 8                  [kg/indiv/timestep]
     real(r8) ::  rootuptake09                    ! net flow into roots (+ into roots), soil layer 9                  [kg/indiv/timestep]
     real(r8) ::  rootuptake10                    ! net flow into roots (+ into roots), soil layer 10                 [kg/indiv/timestep]
                                                  ! BC PLANT HYDRAULICS - flags
     logical ::   is_newly_recuited               !whether the new cohort is newly recuited
  end type ed_cohort_hydr_type
   
 contains

    ! ===================================================================================

    subroutine InitHydrSite(this)
       
       ! Arguments
       class(ed_site_hydr_type),intent(inout) :: this

       allocate(this%v_shell(1:nlevsoi_hyd,1:nshell))         ; this%v_shell = nan
       allocate(this%v_shell_init(1:nlevsoi_hyd,1:nshell))    ; this%v_shell_init = nan
       allocate(this%v_shell_1D(1:nshell))                    ; this%v_shell_1D = nan
       allocate(this%r_node_shell(1:nlevsoi_hyd,1:nshell))    ; this%r_node_shell = nan
       allocate(this%r_node_shell_init(1:nlevsoi_hyd,1:nshell)); this%r_node_shell_init = nan
       allocate(this%r_out_shell(1:nlevsoi_hyd,1:nshell))     ; this%r_out_shell = nan
       allocate(this%l_aroot_layer(1:nlevsoi_hyd))            ; this%l_aroot_layer = nan
       allocate(this%l_aroot_layer_init(1:nlevsoi_hyd))       ; this%l_aroot_layer_init = nan
       allocate(this%kmax_upper_shell(1:nlevsoi_hyd,1:nshell)); this%kmax_upper_shell = nan
       allocate(this%kmax_bound_shell(1:nlevsoi_hyd,1:nshell)); this%kmax_bound_shell = nan
       allocate(this%kmax_lower_shell(1:nlevsoi_hyd,1:nshell)); this%kmax_lower_shell = nan
       allocate(this%r_out_shell_1D(1:nshell))                ; this%r_out_shell_1D = nan
       allocate(this%r_node_shell_1D(1:nshell))               ; this%r_node_shell_1D = nan
       allocate(this%kmax_upper_shell_1D(1:nshell))           ; this%kmax_upper_shell_1D = nan
       allocate(this%kmax_bound_shell_1D(1:nshell))           ; this%kmax_bound_shell_1D = nan
       allocate(this%kmax_lower_shell_1D(1:nshell))           ; this%kmax_lower_shell_1D = nan
       allocate(this%supsub_flag(nlevsoi_hyd))                ; this%supsub_flag = -999
       allocate(this%h2osoi_liqvol_shell(1:nlevsoi_hyd,1:nshell)) ; this%h2osoi_liqvol_shell = nan
       allocate(this%h2osoi_liq_prev(1:nlevsoi_hyd))          ; this%h2osoi_liq_prev = nan
       allocate(this%psisoi_liq_innershell(1:nlevsoi_hyd)); this%psisoi_liq_innershell = nan
       allocate(this%rs1(1:nlevsoi_hyd)); this%rs1(:) = fine_root_radius_const
       allocate(this%recruit_w_uptake(1:nlevsoi_hyd)); this%recruit_w_uptake = nan
       
       this%l_aroot_1D = nan
       this%errh2o_hyd     = nan
       this%dwat_veg       = nan
       this%h2oveg         = nan
       this%h2oveg_dead    = 0.0_r8
       return
    end subroutine InitHydrSite
    
    ! ===================================================================================
    
    subroutine InitHydraulicsDerived()
    
       integer :: k   ! Pool counting index

       do k = 1,n_porous_media
          
          if (k.eq.1) then   ! Leaf tissue
             cap_slp(k)    = 0.0_r8
             cap_int(k)    = 0.0_r8
             cap_corr(k)   = 1.0_r8
          else               ! Non leaf tissues
             cap_slp(k)    = (hydr_psi0 - hydr_psicap )/(1.0_r8 - rwccap(k))  
             cap_int(k)    = -cap_slp(k) + hydr_psi0    
             cap_corr(k)   = -cap_int(k)/cap_slp(k)
          end if
       end do

       return
    end subroutine InitHydraulicsDerived



end module FatesHydraulicsMemMod
