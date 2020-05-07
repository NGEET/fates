module FatesHydraulicsMemMod

   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : fates_unset_r8
   use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=)
   use FatesConstantsMod, only : itrue,ifalse
   use FatesHydroWTFMod,  only : wrf_arr_type
   use FatesHydroWTFMod,  only : wkf_arr_type
   
   implicit none
   private

   logical, parameter, public :: use_2d_hydrosolve = .false.
   
   
   ! Number of soil layers for indexing cohort fine root quanitities
   ! NOTE: The hydraulics code does have some capacity to run a single soil
   ! layer that was developed for comparisons with TFS. However, this has 
   ! not maintained support through the integration with FATES and its
   ! communications with the LSM.  Please do not set nlevsoi_hyd_max
   ! to 1 unless you are developing and testing.

   integer, parameter, public                  :: nlevsoi_hyd_max = 40

   ! number of distinct types of plant porous media (leaf, stem, troot, aroot)
   integer, parameter, public                  :: n_porous_media = 5
   integer, parameter, public                  :: n_plant_media  = 4
   integer, parameter, public                  :: n_hypool_leaf  = 1
   integer, parameter, public                  :: n_hypool_stem  = 1
   integer, parameter, public                  :: n_hypool_troot = 1 ! CANNOT BE CHANGED
   integer, parameter, public                  :: n_hypool_aroot = 1 ! THIS IS "PER-SOIL-LAYER"
   integer, parameter, public                  :: nshell         = 5

   ! number of aboveground plant water storage nodes
   integer, parameter, public                  :: n_hypool_ag    = n_hypool_leaf+n_hypool_stem

   ! total number of water storage nodes
   integer, parameter, public                  :: n_hypool_tot = n_hypool_ag + n_hypool_troot + n_hypool_aroot + nshell
   integer, parameter, public                  :: n_hypool_plant = n_hypool_tot - nshell


   ! vector indexing the type of porous medium over an arbitrary number of plant pools

   integer, parameter, public :: stomata_p_media = 0
   integer, parameter, public :: leaf_p_media  = 1
   integer, parameter, public :: stem_p_media  = 2
   integer, parameter, public :: troot_p_media = 3
   integer, parameter, public :: aroot_p_media = 4
   integer, parameter, public :: rhiz_p_media  = 5

   ! P-V curve: total RWC @ which elastic drainage begins (tfs)     [-]        
   real(r8), parameter, public, dimension(n_plant_media) :: rwcft   = (/1.0_r8,0.958_r8,0.958_r8,0.958_r8/)
   ! P-V curve: total RWC @ which capillary reserves exhausted (tfs)
   real(r8), parameter, public, dimension(n_plant_media) :: rwccap  = (/1.0_r8,0.947_r8,0.947_r8,0.947_r8/) 
   

   ! mirror of nlevcan, hard-set for simplicity, remove nlevcan_hyd on a rainy day
   ! Note (RGK): uscing nclmax causes annoying circular dependency (this needs EDTypes, EDTypes needs this)
   ! way around that: dynamic allocation, or just keep this, but set the value high
   integer, parameter, public                          :: nlevcan_hyd = 2                       
                       
   ! Mean fine root radius expected in the bulk soil                
   real(r8), parameter, public                         :: fine_root_radius_const = 0.0001_r8

   ! Should we ignore the first soil layer and have root layers start on the second?
   logical, parameter, public :: ignore_layer1=.true.
   
   
   ! Derived parameters
   ! ----------------------------------------------------------------------------------------------
   
   !temporatory variables
   real(r8), public :: cohort_recruit_water_layer(nlevsoi_hyd_max)   ! the recruit water requirement for a 
                                                             ! single individual at different layer (kg H2o/m2)
   real(r8), public :: recruit_water_avail_layer(nlevsoi_hyd_max)    ! the recruit water avaibility from soil (kg H2o/m2)

   type, public :: ed_site_hydr_type

      ! Plant Hydraulics
     integer :: i_rhiz_t                            ! Soil layer index of top rhizosphere
     integer :: i_rhiz_b                            ! Soil layer index of bottom rhizospher layer
     integer               :: nlevrhiz              ! Number of rhizosphere levels (vertical layers)
     real(r8), allocatable :: zi_rhiz(:)            ! Depth of the bottom edge of each rhizosphere level [m]
     real(r8), allocatable :: dz_rhiz(:)            ! Width of each rhizosphere level [m]

     real(r8),allocatable :: v_shell(:,:)           ! Volume of rhizosphere compartment (m3) over the
                                                    ! entire site (ha), absolute quantity
     real(r8),allocatable :: v_shell_init(:,:)      ! Previous volume of rhizosphere compartment (m3) 
     real(r8),allocatable :: r_node_shell(:,:)      ! Nodal radius of rhizosphere compartment (m)
     real(r8),allocatable :: r_node_shell_init(:,:) ! Previous Nodal radius of rhizosphere compartment (m)
     real(r8),allocatable :: l_aroot_layer(:)       ! Total length (across cohorts) of absorbing
                                                    !  roots by soil layer (m)
     real(r8),allocatable :: l_aroot_layer_init(:)  ! Total length (across cohorts) of absorbing
                                                    !  roots by soil layer (m)
     real(r8),allocatable :: kmax_upper_shell(:,:)  ! Maximum soil hydraulic conductance node k 
                                                    ! to upper (closer to atmosphere) rhiz 
                                                    ! shell boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: kmax_lower_shell(:,:)  ! Maximum soil hydraulic conductance node k
                                                    ! to lower (further from atmosphere) 
                                                    ! rhiz shell boundaries (kg s-1 MPa-1)
     real(r8),allocatable :: r_out_shell(:,:)       ! Outer radius of rhizosphere compartment (m)


     real(r8),allocatable :: rs1(:)                 ! Mean fine root radius (m) (currently a constant)

     integer, allocatable :: supsub_flag(:)         ! index of the outermost rhizosphere shell 
                                                    ! encountering super- or sub-saturation
     real(r8),allocatable :: h2osoi_liqvol_shell(:,:) ! volumetric water in rhizosphere compartment (m3/m3)

     real(r8),allocatable :: h2osoi_liq_prev(:)     ! liquid water mass for the bulk soil layer
                                                    ! defined at the end of the hydraulics sequence
                                                    ! after root water has been extracted.  This should
                                                    ! be equal to the sum of the water over the rhizosphere shells
                                                    ! [kg/m2]
     
     real(r8),allocatable :: recruit_w_uptake(:)    ! recruitment water uptake (kg H2o/m2/s)

    
     real(r8) :: errh2o_hyd                         ! plant hydraulics error summed across 
                                                    ! cohorts to column level (mm)
     real(r8) :: dwat_veg                           ! change in stored water in vegetation
                                                    ! column level (kg)
     real(r8) :: h2oveg                             ! stored water in vegetation (kg/m2)

     real(r8) :: h2oveg_recruit                     ! stored water in recruits (kg/m2)
     real(r8) :: h2oveg_dead                        ! stored water in dead vegetation (kg/m2)
     real(r8) :: h2oveg_growturn_err                ! error water pool (kg/m2) for increase (growth) or
                                                    !  contraction (turnover) of tissue volumes.
                                                    !  Draw from or add to this pool when
                                                    !  insufficient water available to increase
                                                    !  tissue volume or too much water is
                                                    !  available when tissue volume decreases,
                                                    !  respectively.
     real(r8) :: h2oveg_pheno_err                   ! error water pool (kg/m2) for leaf-on
                                                    !  Draw from or add to this pool when
                                                    !  insufficient plant water available to 
                                                    !  support production of new leaves.
     real(r8) :: h2oveg_hydro_err                   ! error water pool (kg/m2) for hydrodynamics
                                                    !  Draw from or add to this pool when
                                                    !  insufficient plant water available to 
                                                    !  support transpiration


     ! Useful diagnostics
     ! ----------------------------------------------------------------------------------

     real(r8),allocatable ::  sapflow_scpf(:,:)   ! flow at base of tree (+ upward)      [kg/ha/s]
                                                  ! discretized by size x pft

     ! Root uptake per rhiz layer [kg/ha/s]
     real(r8),allocatable :: rootuptake_sl(:)

     ! Root uptake per pft x size class, over set layer depths [kg/ha/m/s]
     ! These are normalized by depth (in case the desired horizon extends
     ! beyond the actual rhizosphere)
     
     real(r8), allocatable :: rootuptake0_scpf(:,:)   ! 0-10 cm
     real(r8), allocatable :: rootuptake10_scpf(:,:)  ! 10-50 cm
     real(r8), allocatable :: rootuptake50_scpf(:,:)  ! 50-100 cm
     real(r8), allocatable :: rootuptake100_scpf(:,:) ! 100+ cm

     
     
     class(wrf_arr_type), pointer :: wrf_soil(:)       ! Water retention function for soil layers
     class(wkf_arr_type), pointer :: wkf_soil(:)       ! Water conductivity (K) function for soil

     ! For the matrix version of the solver we need to define the connection
     ! and type map for the whole system of compartments, from the soil to leaf
     ! as one vector
     
     integer :: num_connections
     integer :: num_nodes
     integer, allocatable :: conn_up(:)
     integer, allocatable :: conn_dn(:)
     integer, allocatable :: pm_node(:)
     integer, allocatable :: node_layer(:)
     integer, allocatable :: ipiv(:)       ! unused, returned from DSEGV
     
     real(r8), allocatable :: residual(:)
     real(r8), allocatable :: ajac(:,:)
     real(r8), allocatable :: th_node_init(:)
     real(r8), allocatable :: th_node(:)
     real(r8), allocatable :: dth_node(:)
     real(r8), allocatable :: h_node(:)
     real(r8), allocatable :: v_node(:)
     real(r8), allocatable :: z_node(:)
     real(r8), allocatable :: psi_node(:)
     real(r8), allocatable :: q_flux(:)
     real(r8), allocatable :: dftc_dpsi_node(:)
     real(r8), allocatable :: ftc_node(:)
     

     real(r8), allocatable :: kmax_up(:)
     real(r8), allocatable :: kmax_dn(:)
     
     
  contains
     
    procedure :: InitHydrSite
    procedure :: SetConnections
    procedure :: FlushSiteScratch
  end type ed_site_hydr_type



  type, public :: ed_cohort_hydr_type


     ! Node heights of compartments [m]
     ! Heights are referenced to soil surface (+ = above; - = below)
     ! Note* The node centers of the absorbing root compartments, are the same
     ! as the soil layer mid-points that they occupy, so no need to save those.
     ! ----------------------------------------------------------------------------------

     real(r8) :: z_node_ag(n_hypool_ag)  ! nodal height of stem and leaf compartments (positive)
     real(r8) :: z_upper_ag(n_hypool_ag) ! height of upper stem and leaf compartment boundaries (positive)
     real(r8) :: z_lower_ag(n_hypool_ag) ! height of lower stem and leaf compartment boundaries (positive)
     real(r8) :: z_node_troot            ! height of transporting root node


     ! Maximum hydraulic conductances  [kg H2O s-1 MPa-1]
     ! ----------------------------------------------------------------------------------

     real(r8) :: kmax_petiole_to_leaf                  ! Max conductance, petiole to leaf
                                                       ! Nominally set to very high value 
     real(r8) :: kmax_stem_upper(n_hypool_stem)        ! Max conductance, upper stem compartments
     real(r8) :: kmax_stem_lower(n_hypool_stem)        ! Max conductance, lower stem compartments
     real(r8) :: kmax_troot_upper                      ! Max conductance, uper portion of the
                                                       ! transporting root
     real(r8),allocatable :: kmax_troot_lower(:)       ! Max conductance in portion of transporting
                                                       ! root compartment that joins each absorbing
                                                       ! root compartment
     real(r8),allocatable :: kmax_aroot_upper(:)       ! Max conductance in the absorbing root
                                                       ! compartment through xylem tissues going
                                                       ! into the transporting root
     real(r8),allocatable :: kmax_aroot_lower(:)       ! Since this pools may actually be a
                                                       ! hybrid that contains transporting
                                                       ! root volume, then we need to factor
                                                       ! in xylem resistance from the absorbing
                                                       ! root edge to the node center

                                                       ! Max conductance in the absorbing
                                                       ! root compartment, radially through the
                                                       ! exodermis, cortex, casparian strip, and 
                                                       ! endodermis, separated for two cases, when:
     real(r8),allocatable :: kmax_aroot_radial_in(:)   ! the potential gradient is positive "into" root
     real(r8),allocatable :: kmax_aroot_radial_out(:)  ! the potential gradient is positive "out of" root


     ! Compartment Volumes and lengths

     real(r8) ::  v_ag_init(n_hypool_ag)          ! previous day's volume of aboveground water storage compartments   [m3]
     real(r8) ::  v_ag(n_hypool_ag)               ! volume of aboveground water storage compartments                  [m3]
     real(r8) ::  v_troot_init                    ! previous day's volume of belowground water storage compartments   [m3]
     real(r8) ::  v_troot                         ! volume of belowground water storage compartments                  [m3]
     real(r8),allocatable :: v_aroot_layer_init(:) ! previous day's volume of absorbing roots by soil layer    [m3]
     real(r8),allocatable :: v_aroot_layer(:)      ! volume of absorbing roots by soil layer                   [m3]
     real(r8),allocatable :: l_aroot_layer(:)      ! length of absorbing roots by soil layer                   [m]
     

     
     ! State variable, relative water content by volume (i.e. "theta")
     real(r8) :: th_ag(n_hypool_ag)              ! water in aboveground compartments                                 [kgh2o/indiv]
     real(r8) :: th_troot                        ! water in belowground compartments                                 [kgh2o/indiv]
     real(r8),allocatable :: th_aroot(:)          ! water in absorbing roots                                          [kgh2o/indiv]
    

     ! Diagnostic, water potential
     real(r8) :: psi_ag(n_hypool_ag)             ! water potential in aboveground compartments                       [MPa]
     real(r8) :: psi_troot                       ! water potential in belowground compartments                       [MPa]
     real(r8),allocatable :: psi_aroot(:)         ! water potential in absorbing roots                                [MPa]

     ! Diagnostic, fraction of total conductivity
     real(r8) :: ftc_ag(n_hypool_ag)              ! ... in above-ground compartments [-]
     real(r8) :: ftc_troot                        ! ... in the transporting root [-]
     real(r8),allocatable :: ftc_aroot(:)         ! ... in the absorbing root [-]


     real(r8) ::  btran                           ! leaf water potential limitation on gs                             [0-1]

     
     real(r8) ::  qtop                            ! mean transpiration flux rate         [kg/cohort/s]
     
     
     ! Variables used for error tracking and flagging
     ! ----------------------------------------------------------------------------------
     
     real(r8) ::  supsub_flag                     ! k index of last node to encounter supersaturation or 
                                                  ! sub-residual water content  (+ supersaturation; - subsaturation)
     real(r8) ::  iterh1                          ! max number of iterations required to achieve tolerable
                                                  ! water balance error (if 1D, associated with iterlayer)
     real(r8) ::  iterh2                          ! number of inner iterations (if 1D, associated with iterlayer)
     real(r8) ::  iterlayer                       ! layer index associated with the highest iterations

     real(r8) ::  errh2o                          ! total water balance error per unit crown area                     [kgh2o/m2]
     real(r8) ::  errh2o_growturn_ag(n_hypool_ag) ! error water pool for increase (growth) or
                                                  !  contraction (turnover) of tissue volumes.
                                                  !  Draw from or add to this pool when
                                                  !  insufficient water available to increase
                                                  !  tissue volume or too much water is
                                                  !  available when tissue volume decreases,
                                                  !  respectively.
     real(r8) ::  errh2o_pheno_ag(n_hypool_ag)    ! error water pool for for leaf-on
                                                  !  Draw from or add to this pool when
                                                  !  insufficient plant water available to 
                                                  !  support production of new leaves.
     real(r8) ::  errh2o_growturn_troot           ! same as errh2o_growturn_ag but for troot pool
     real(r8) ::  errh2o_pheno_troot              ! same as errh2o_pheno_ag but for troot pool
     real(r8) ::  errh2o_growturn_aroot           ! same as errh2o_growturn_ag but for aroot pools
     real(r8) ::  errh2o_pheno_aroot              ! same as errh2o_pheno_ag but for aroot pools


     

    
     ! Other
     ! ----------------------------------------------------------------------------------
     
     logical ::   is_newly_recruited              ! whether the new cohort is newly recruited



     ! ----------------------------------------------------------------------------------
     ! NOT USED, BUT HOLDING FOR FUTURE RE-IMPLEMENTATION
     !real(r8) ::  flc_min_ag(n_hypool_ag)         ! min attained fractional loss of conductivity in 
     !                                             ! aboveground compartments (for tracking xylem refilling dynamics) [-]
     !real(r8) ::  flc_min_troot(n_hypool_troot)   ! min attained fractional loss of conductivity in 
     !                                             ! belowground compartments (for tracking xylem refilling dynamics) [-]
     !real(r8),allocatable ::  flc_min_aroot(:)    ! min attained fractional loss of conductivity in absorbing roots 
     !                                             ! (for tracking xylem refilling dynamics)          [-]
     !real(r8) ::  lwp_mem(numLWPmem)              ! leaf water potential over the previous numLWPmem timesteps        [MPa] 
     !real(r8) ::  lwp_stable                      ! leaf water potential just before it became unstable               [MPa]
     !logical  ::  lwp_is_unstable                 ! flag for instability of leaf water potential over previous timesteps
     !real(r8) ::  refill_thresh                   ! water potential threshold for xylem refilling to occur            [MPa]
     !real(r8) ::  refill_days                     ! number of days required for 50% of xylem refilling to occur       [days]
     ! -----------------------------------------------------------------------------------
  contains
     
     procedure :: AllocateHydrCohortArrays
     procedure :: DeallocateHydrCohortArrays
     
  end type ed_cohort_hydr_type
   
 contains
    
    subroutine AllocateHydrCohortArrays(this,nlevrhiz)
       
       ! Arguments
       class(ed_cohort_hydr_type),intent(inout) :: this
       integer, intent(in)                      :: nlevrhiz

       allocate(this%kmax_troot_lower(1:nlevrhiz))
       allocate(this%kmax_aroot_upper(1:nlevrhiz))
       allocate(this%kmax_aroot_lower(1:nlevrhiz))
       allocate(this%kmax_aroot_radial_in(1:nlevrhiz))
       allocate(this%kmax_aroot_radial_out(1:nlevrhiz))
       allocate(this%v_aroot_layer_init(1:nlevrhiz))
       allocate(this%v_aroot_layer(1:nlevrhiz))
       allocate(this%l_aroot_layer(1:nlevrhiz))
       allocate(this%th_aroot(1:nlevrhiz))
       allocate(this%psi_aroot(1:nlevrhiz))
       allocate(this%ftc_aroot(1:nlevrhiz))

       return
    end subroutine AllocateHydrCohortArrays

    ! ===================================================================================

    subroutine DeallocateHydrCohortArrays(this)

       class(ed_cohort_hydr_type),intent(inout) :: this
       
       deallocate(this%kmax_troot_lower)
       deallocate(this%kmax_aroot_upper)
       deallocate(this%kmax_aroot_lower)
       deallocate(this%kmax_aroot_radial_in)
       deallocate(this%kmax_aroot_radial_out)
       deallocate(this%v_aroot_layer_init)
       deallocate(this%v_aroot_layer)
       deallocate(this%l_aroot_layer)
       deallocate(this%th_aroot)
       deallocate(this%psi_aroot)
       deallocate(this%ftc_aroot)

       return
    end subroutine DeallocateHydrCohortArrays

    ! ===================================================================================

    subroutine InitHydrSite(this,numpft,numlevsclass)
       
       ! Arguments
       class(ed_site_hydr_type),intent(inout) :: this
       integer,intent(in) :: numpft
       integer,intent(in) :: numlevsclass

       associate(nlevrhiz => this%nlevrhiz)

         allocate(this%zi_rhiz(1:nlevrhiz));  this%zi_rhiz(:)  = nan
         allocate(this%dz_rhiz(1:nlevrhiz)); this%dz_rhiz(:) = nan
         allocate(this%v_shell(1:nlevrhiz,1:nshell))         ; this%v_shell = nan
         allocate(this%v_shell_init(1:nlevrhiz,1:nshell))    ; this%v_shell_init = nan
         allocate(this%r_node_shell(1:nlevrhiz,1:nshell))    ; this%r_node_shell = nan
         allocate(this%r_node_shell_init(1:nlevrhiz,1:nshell)); this%r_node_shell_init = nan
         allocate(this%r_out_shell(1:nlevrhiz,1:nshell))     ; this%r_out_shell = nan
         allocate(this%l_aroot_layer(1:nlevrhiz))            ; this%l_aroot_layer = nan
         allocate(this%l_aroot_layer_init(1:nlevrhiz))       ; this%l_aroot_layer_init = nan
         allocate(this%kmax_upper_shell(1:nlevrhiz,1:nshell)); this%kmax_upper_shell = nan
         allocate(this%kmax_lower_shell(1:nlevrhiz,1:nshell)); this%kmax_lower_shell = nan
         allocate(this%supsub_flag(1:nlevrhiz))                ; this%supsub_flag = -999
         allocate(this%h2osoi_liqvol_shell(1:nlevrhiz,1:nshell)) ; this%h2osoi_liqvol_shell = nan
         allocate(this%h2osoi_liq_prev(1:nlevrhiz))          ; this%h2osoi_liq_prev = nan
         allocate(this%rs1(1:nlevrhiz)); this%rs1(:) = fine_root_radius_const
         allocate(this%recruit_w_uptake(1:nlevrhiz)); this%recruit_w_uptake = nan
         
         allocate(this%sapflow_scpf(1:numlevsclass,1:numpft))       ; this%sapflow_scpf = nan
         allocate(this%rootuptake_sl(1:nlevrhiz))                   ; this%rootuptake_sl = nan
         allocate(this%rootuptake0_scpf(1:numlevsclass,1:numpft))   ; this%rootuptake0_scpf = nan
         allocate(this%rootuptake10_scpf(1:numlevsclass,1:numpft))  ; this%rootuptake10_scpf = nan
         allocate(this%rootuptake50_scpf(1:numlevsclass,1:numpft))  ; this%rootuptake50_scpf = nan
         allocate(this%rootuptake100_scpf(1:numlevsclass,1:numpft)) ; this%rootuptake100_scpf = nan
         
         this%errh2o_hyd     = nan
         this%dwat_veg       = nan
         this%h2oveg         = 0.0_r8
         this%h2oveg_recruit = 0.0_r8
         this%h2oveg_dead    = 0.0_r8
         this%h2oveg_growturn_err = 0.0_r8
         this%h2oveg_pheno_err    = 0.0_r8
         this%h2oveg_hydro_err    = 0.0_r8
         
         ! We have separate water transfer functions and parameters
         ! for each soil layer, and each plant compartment type
         allocate(this%wrf_soil(1:nlevrhiz))
         allocate(this%wkf_soil(1:nlevrhiz))
         
         if(use_2d_hydrosolve) then
            
            this%num_connections =  n_hypool_leaf + n_hypool_stem + n_hypool_troot - 1  &
                 + (n_hypool_aroot + nshell) * nlevrhiz
            
            this%num_nodes = n_hypool_leaf + n_hypool_stem + n_hypool_troot  &
                 + (n_hypool_aroot + nshell) * nlevrhiz
            
            ! These are only in the newton-matrix solve
            allocate(this%conn_up(this%num_connections))
            allocate(this%conn_dn(this%num_connections))
            allocate(this%residual(this%num_nodes))
            allocate(this%ajac(this%num_nodes,this%num_nodes))
            allocate(this%th_node_init(this%num_nodes))
            allocate(this%th_node(this%num_nodes))
            allocate(this%dth_node(this%num_nodes))
            allocate(this%h_node(this%num_nodes))
            allocate(this%v_node(this%num_nodes))
            allocate(this%z_node(this%num_nodes))
            allocate(this%psi_node(this%num_nodes))
            allocate(this%q_flux(this%num_connections))
            allocate(this%dftc_dpsi_node(this%num_nodes))
            allocate(this%ftc_node(this%num_nodes))
            allocate(this%pm_node(this%num_nodes))
            allocate(this%ipiv(this%num_nodes))
            allocate(this%node_layer(this%num_nodes))
            
            allocate(this%kmax_up(this%num_connections))
            allocate(this%kmax_dn(this%num_connections))
            
         else
            
            this%num_connections =  n_hypool_leaf + n_hypool_stem + & 
                 n_hypool_troot + n_hypool_aroot + nshell -1 
            
            this%num_nodes = n_hypool_leaf + n_hypool_stem + & 
                   n_hypool_troot + n_hypool_aroot + nshell
            
            allocate(this%conn_up(this%num_connections))
            allocate(this%conn_dn(this%num_connections))
            allocate(this%pm_node(this%num_nodes))
            
            
         end if
         
         call this%SetConnections()
         
          
       end associate
       
       return
    end subroutine InitHydrSite
     
    ! ===================================================================================
    
    subroutine FlushSiteScratch(this)
        class(ed_site_hydr_type),intent(inout) :: this

        if(use_2d_hydrosolve) then
            this%residual(:)       = fates_unset_r8
            this%ajac(:,:)         = fates_unset_r8
            this%th_node_init(:)   = fates_unset_r8
            this%th_node(:)        = fates_unset_r8
            this%dth_node(:)       = fates_unset_r8
            this%h_node(:)         = fates_unset_r8
            this%v_node(:)         = fates_unset_r8
            this%z_node(:)         = fates_unset_r8
            this%psi_node(:)       = fates_unset_r8
            this%ftc_node(:)       = fates_unset_r8
            this%dftc_dpsi_node(:) = fates_unset_r8
!            this%kmax_up(:)        = fates_unset_r8
!            this%kmax_dn(:)        = fates_unset_r8
            this%q_flux(:)         = fates_unset_r8
        end if

    end subroutine FlushSiteScratch

    ! ===================================================================================

    subroutine SetConnections(this)
      
     class(ed_site_hydr_type),intent(inout) :: this
      
     integer :: k, j
     integer :: num_cnxs
     integer :: num_nds
     integer :: nt_ab
     integer :: node_tr_end
     
     num_cnxs = 0
     num_nds = 0
     do k = 1, n_hypool_leaf
        num_cnxs = num_cnxs + 1
        num_nds  = num_nds + 1
        this%conn_dn(num_cnxs) = k           !leaf is the dn, origin, bottom
        this%conn_up(num_cnxs) = k + 1
        this%pm_node(num_nds)  = leaf_p_media
     enddo
     do k = n_hypool_leaf+1, n_hypool_ag
        num_cnxs = num_cnxs + 1
        num_nds  = num_nds + 1
        this%conn_dn(num_cnxs) = k
        this%conn_up(num_cnxs) = k+1
        this%pm_node(num_nds) = stem_p_media
     enddo

     if(use_2d_hydrosolve) then
     
        num_nds     = n_hypool_ag+n_hypool_troot
        node_tr_end = num_nds
        nt_ab       = n_hypool_ag+n_hypool_troot+n_hypool_aroot
        num_cnxs    = n_hypool_ag

        this%pm_node(num_nds) = troot_p_media
        this%node_layer(1:n_hypool_ag) = 0
        this%node_layer(num_nds) = 1
        
        do j = 1,this%nlevrhiz
           do k = 1, n_hypool_aroot + nshell
              num_nds  = num_nds + 1
              num_cnxs = num_cnxs + 1
              this%node_layer(num_nds) = j
              if( k == 1 ) then !troot-aroot
                 !junction node
                 this%conn_dn(num_cnxs) = node_tr_end !absorbing root
                 this%conn_up(num_cnxs) = num_nds
                 this%pm_node(num_nds)  = aroot_p_media
              else
                 this%conn_dn(num_cnxs) = num_nds - 1
                 this%conn_up(num_cnxs) = num_nds
                 this%pm_node(num_nds)  = rhiz_p_media
              endif
           enddo
        end do
     else
        
        this%pm_node(n_hypool_ag+1) = troot_p_media
        this%pm_node(n_hypool_ag+2) = aroot_p_media
        this%pm_node(n_hypool_ag+3:n_hypool_ag+2+nshell) = rhiz_p_media
        
     end if
     
   end subroutine SetConnections
    

end module FatesHydraulicsMemMod
