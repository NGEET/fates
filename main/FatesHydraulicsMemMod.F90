module FatesHydraulicsMemMod

   use FatesConstantsMod, only : r8 => fates_r8
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use FatesConstantsMod, only : itrue,ifalse
   use EDParamsMod      , only : hydr_psi0
   use EDParamsMod      , only : hydr_psicap
   
   implicit none
   private

   ! Number of soil layers for indexing cohort fine root quanitities
   ! NOTE: The hydraulics code does have some capacity to run a single soil
   ! layer that was developed for comparisons with TFS. However, this has 
   ! not maintained support through the integration with FATES and its
   ! communications with the LSM.  Please do not set nlevsoi_hyd_max
   ! to 1 unless you are developing and testing.

   integer, parameter, public                  :: nlevsoi_hyd_max = 40

   ! number of distinct types of plant porous media (leaf, stem, troot, aroot)
   integer, parameter, public                  :: n_porous_media = 4

   integer, parameter, public                  :: n_hypool_leaf  = 1
   integer, parameter, public                  :: n_hypool_stem  = 1
   integer, parameter, public                  :: n_hypool_troot = 1 ! CANNOT BE CHANGED
   integer, parameter, public                  :: n_hypool_aroot = 1 ! THIS IS "PER-SOIL-LAYER"
   integer, parameter, public                  :: nshell         = 5

   ! number of aboveground plant water storage nodes
   integer, parameter, public                  :: n_hypool_ag    = n_hypool_leaf+n_hypool_stem

   ! total number of water storage nodes
   integer, parameter, public                  :: n_hypool_tot = n_hypool_ag + n_hypool_troot + n_hypool_aroot + nshell

   ! vector indexing the type of porous medium over an arbitrary number of plant pools

   integer, parameter, public :: leaf_p_media  = 1
   integer, parameter, public :: stem_p_media  = 2
   integer, parameter, public :: troot_p_media = 3
   integer, parameter, public :: aroot_p_media = 4
   integer, parameter, public :: rhiz_p_media  = 5

   ! This vector holds the identifiers for which porous media type is in the comaprtment
   integer, parameter, public, dimension(n_hypool_tot) :: porous_media = (/leaf_p_media,  & 
                                                                           stem_p_media,  & 
                                                                           troot_p_media, & 
                                                                           aroot_p_media, & 
                                                                           rhiz_p_media,  & 
                                                                           rhiz_p_media,  & 
                                                                           rhiz_p_media,  & 
                                                                           rhiz_p_media,  & 
                                                                           rhiz_p_media /) 

   ! number of previous timestep's leaf water potential to be retained
   integer, parameter, public                          :: numLWPmem             = 4

   ! mirror of nlevcan, hard-set for simplicity, remove nlevcan_hyd on a rainy day
   ! Note (RGK): uscing nclmax causes annoying circular dependency (this needs EDTypes, EDTypes needs this)
   ! way around that: dynamic allocation, or just keep this, but set the value high
   integer, parameter, public                          :: nlevcan_hyd = 2                       
                       
   ! Mean fine root radius expected in the bulk soil                
   real(r8), parameter, public                         :: fine_root_radius_const = 0.001_r8               
   
   ! Constant parameters (for time being, C2B is constant, 
   ! slated for addition to parameter file (RGK 08-2017))
   ! Carbon 2 biomass ratio
   real(r8), parameter, public                         :: C2B        = 2.0_r8 
              
   ! P-V curve: total RWC @ which elastic drainage begins     [-]        
   real(r8), parameter, public, dimension(n_porous_media) :: rwcft   = (/1.0_r8,0.958_r8,0.958_r8,0.958_r8/)

   ! P-V curve: total RWC @ which capillary reserves exhausted
   real(r8), parameter, public, dimension(n_porous_media) :: rwccap  = (/1.0_r8,0.947_r8,0.947_r8,0.947_r8/) 

   ! Derived parameters
   ! ----------------------------------------------------------------------------------------------

   ! P-V curve: slope of capillary region of curve
   real(r8), public, dimension(n_porous_media)           :: cap_slp                                         

   ! P-V curve: intercept of capillary region of curve
   real(r8), public, dimension(n_porous_media)           :: cap_int   

   ! P-V curve: correction for nonzero psi0x
   real(r8), public, dimension(n_porous_media)           :: cap_corr                                        
   
   !temporatory variables
   real(r8), public :: cohort_recruit_water_layer(nlevsoi_hyd_max)   ! the recruit water requirement for a 
                                                             ! single individual at different layer (kg H2o/m2)
   real(r8), public :: recruit_water_avail_layer(nlevsoi_hyd_max)    ! the recruit water avaibility from soil (kg H2o/m2)

   type, public :: ed_site_hydr_type

      ! Plant Hydraulics
     
     integer              :: nlevsoi_hyd            ! The number of soil hydraulic layers
                                                    ! the host model may offer different number of
                                                    ! layers for every site, and hydraulics
                                                    ! may or may not cross that with a simple or
                                                    ! non-simple layering

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

     integer,allocatable :: supsub_flag(:)          ! index of the outermost rhizosphere shell 
                                                    ! encountering super- or sub-saturation
     real(r8),allocatable :: h2osoi_liqvol_shell(:,:) ! volumetric water in rhizosphere compartment (m3/m3)

     real(r8),allocatable :: h2osoi_liq_prev(:)     ! liquid water mass for the bulk soil layer
                                                    ! defined at the end of the hydraulics sequence
                                                    ! after root water has been extracted.  This should
                                                    ! be equal to the sum of the water over the rhizosphere shells
     
     real(r8),allocatable :: psisoi_liq_innershell(:) ! Matric potential of the inner rhizosphere shell (MPa)
     
     
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

     
     ! Variables used for error tracking and flagging
     ! ----------------------------------------------------------------------------------
     
     real(r8) ::  supsub_flag                     ! k index of last node to encounter supersaturation or 
                                                  ! sub-residual water content  (+ supersaturation; - subsaturation)
     real(r8) ::  iterh1                          ! number of iterations required to achieve tolerable water balance error
     real(r8) ::  iterh2                          ! number of inner iterations
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


     
     ! Useful diagnostics
     ! ----------------------------------------------------------------------------------

     real(r8) ::  sapflow                         ! flow at base of tree (+ upward)                                   [kg/indiv/timestep]
     real(r8),allocatable ::  rootuptake(:)       ! net flow into roots (+ into roots)                                [kg/indiv/timestep]
                                                  ! BC PLANT HYDRAULICS - flags

    
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
   
  ! Make public necessary subroutines and functions
  public :: InitHydraulicsDerived

 contains
    
    subroutine AllocateHydrCohortArrays(this,nlevsoil_hydr)
       
       ! Arguments
       class(ed_cohort_hydr_type),intent(inout) :: this
       integer, intent(in)                      :: nlevsoil_hydr

       allocate(this%kmax_troot_lower(1:nlevsoil_hydr))
       allocate(this%kmax_aroot_upper(1:nlevsoil_hydr))
       allocate(this%kmax_aroot_radial_in(1:nlevsoil_hydr))
       allocate(this%kmax_aroot_radial_out(1:nlevsoil_hydr))
       allocate(this%v_aroot_layer_init(1:nlevsoil_hydr))
       allocate(this%v_aroot_layer(1:nlevsoil_hydr))
       allocate(this%l_aroot_layer(1:nlevsoil_hydr))
       allocate(this%th_aroot(1:nlevsoil_hydr))
       allocate(this%psi_aroot(1:nlevsoil_hydr))
       allocate(this%ftc_aroot(1:nlevsoil_hydr))
       allocate(this%rootuptake(1:nlevsoil_hydr))

       return
    end subroutine AllocateHydrCohortArrays

    ! ===================================================================================

    subroutine DeallocateHydrCohortArrays(this)

       class(ed_cohort_hydr_type),intent(inout) :: this
       
       deallocate(this%kmax_troot_lower)
       deallocate(this%kmax_aroot_upper)
       deallocate(this%kmax_aroot_radial_in)
       deallocate(this%kmax_aroot_radial_out)
       deallocate(this%v_aroot_layer_init)
       deallocate(this%v_aroot_layer)
       deallocate(this%l_aroot_layer)
       deallocate(this%th_aroot)
       deallocate(this%psi_aroot)
       deallocate(this%ftc_aroot)
       deallocate(this%rootuptake)

       return
    end subroutine DeallocateHydrCohortArrays

    ! ===================================================================================

    subroutine InitHydrSite(this)
       
       ! Arguments
       class(ed_site_hydr_type),intent(inout) :: this

       associate( nlevsoil_hyd => this%nlevsoi_hyd )
         
         allocate(this%v_shell(1:nlevsoil_hyd,1:nshell))         ; this%v_shell = nan
         allocate(this%v_shell_init(1:nlevsoil_hyd,1:nshell))    ; this%v_shell_init = nan
         allocate(this%r_node_shell(1:nlevsoil_hyd,1:nshell))    ; this%r_node_shell = nan
         allocate(this%r_node_shell_init(1:nlevsoil_hyd,1:nshell)); this%r_node_shell_init = nan
         allocate(this%r_out_shell(1:nlevsoil_hyd,1:nshell))     ; this%r_out_shell = nan
         allocate(this%l_aroot_layer(1:nlevsoil_hyd))            ; this%l_aroot_layer = nan
         allocate(this%l_aroot_layer_init(1:nlevsoil_hyd))       ; this%l_aroot_layer_init = nan
         allocate(this%kmax_upper_shell(1:nlevsoil_hyd,1:nshell)); this%kmax_upper_shell = nan
         allocate(this%kmax_lower_shell(1:nlevsoil_hyd,1:nshell)); this%kmax_lower_shell = nan
         allocate(this%supsub_flag(1:nlevsoil_hyd))                ; this%supsub_flag = -999
         allocate(this%h2osoi_liqvol_shell(1:nlevsoil_hyd,1:nshell)) ; this%h2osoi_liqvol_shell = nan
         allocate(this%h2osoi_liq_prev(1:nlevsoil_hyd))          ; this%h2osoi_liq_prev = nan
         allocate(this%psisoi_liq_innershell(1:nlevsoil_hyd)); this%psisoi_liq_innershell = nan
         allocate(this%rs1(1:nlevsoil_hyd)); this%rs1(:) = fine_root_radius_const
         allocate(this%recruit_w_uptake(1:nlevsoil_hyd)); this%recruit_w_uptake = nan

         this%errh2o_hyd     = nan
         this%dwat_veg       = nan
         this%h2oveg         = 0.0_r8
         this%h2oveg_recruit = 0.0_r8
         this%h2oveg_dead    = 0.0_r8
         this%h2oveg_growturn_err = 0.0_r8
         this%h2oveg_pheno_err    = 0.0_r8
         this%h2oveg_hydro_err    = 0.0_r8
         
       end associate

       return
    end subroutine InitHydrSite
    
    ! ===================================================================================
    
    subroutine InitHydraulicsDerived(numpft)
    
    !use EDPftvarcon,       only : EDPftvarcon_inst
       ! Arguments
       integer,intent(in)                      :: numpft
    
       integer :: k   ! Pool counting index
       integer :: ft

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
       
       do ft=1,numpft
          ! this needs a -999 check (BOC)
          !EDPftvarcon_inst%hydr_pinot_node(ft,:) = EDPftvarcon_inst%hydr_pitlp_node(ft,:) * &
          !                                         EDPftvarcon_inst%hydr_epsil_node(ft,:) / &
          !                                        (EDPftvarcon_inst%hydr_epsil_node(ft,:) - &
          !                                         EDPftvarcon_inst%hydr_pitlp_node(ft,:))
       end do

       return
    end subroutine InitHydraulicsDerived



end module FatesHydraulicsMemMod
