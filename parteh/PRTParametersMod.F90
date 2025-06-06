module PRTParametersMod


  use FatesConstantsMod,     only : r8 => fates_r8
  
  ! This module only holds the parameter definitions for PARTEH and allometry.
  ! This does not hold any of the code used for intiailizing and filling
  ! that data, for that is model dependent (ie FATES may have a different
  ! way than another TBM)
  ! This code does perform checks on parameters.

  type,public ::  prt_param_type

     integer, allocatable :: phen_leaf_habit(:)           ! Leaf phenological habit? Current options include the following:
                                                          !    (actual values defined in FatesConstantsMod.F90)
                                                          !    - ievergreen - evergreen
                                                          !    - ihard_season_decid - obligate cold deciduous (i.e., 
                                                          !      leaves will abscise and flush every winter)
                                                          !    - ihard_stress_decid - obligate drought deciduous (i.e.,
                                                          !      leaves will abscie and flush at least once a year)
                                                          !    - isemi_stress_decid - drought semi-deciduous (i.e., 
                                                          !      partial abscission and flushing are allowed).
     ! Drop fraction for tissues other than leaves (PFT-dependent)
     real(r8), allocatable :: phen_fnrt_drop_fraction(:)  ! Abscission fraction of fine roots
     real(r8), allocatable :: phen_stem_drop_fraction(:)  ! Abscission fraction of stems
     real(r8), allocatable :: phen_drought_threshold(:)   ! For obligate (hard) drought deciduous, this is the threshold
                                                          !    below which plants will abscise leaves, and 
                                                          !    above which plants will flush leaves. For semi-deciduous
                                                          !    plants, this is the threshold below which abscission will
                                                          !    be complete. This depends on the sign. If positive, these
                                                          !    are soil volumetric water content [m3/m3]. If
                                                          !    negative, the values are soil matric potential [mm].
                                                          !    Ignored for non-deciduous plants.
     real(r8), allocatable :: phen_moist_threshold(:)     ! For semi-deciduous, this is the threshold above which flushing 
                                                          !    will be complete.  This depends on the sign. If positive, these
                                                          !    are soil volumetric water content [m3/m3]. If
                                                          !    negative, the values are soil matric potential [mm].
                                                          !    Ignored for non-deciduous plants.
     real(r8), allocatable :: phen_doff_time(:)           ! Minimum number of days that plants must remain leafless before
                                                          !   flushing leaves again.


     ! Growth and Turnover Parameters
     real(r8), allocatable :: senleaf_long_fdrought(:)   ! Multiplication factor for leaf longevity of senescent 
                                                         ! leaves during drought( 1.0 indicates no change)
     real(r8), allocatable :: leaf_long(:,:)             ! Leaf turnover time (longevity) (pft x age-class)
                                                         ! If there is >1 class, it is the longevity from
                                                         ! one class to the next [yr]
                                                         !   For drought-deciduous PFTs, the sum of leaf
                                                         !   longevity across all leaf age classes is also
                                                         !   the maximum length of the growing (i.e., leaves on)
                                                         !   season.
     real(r8), allocatable :: leaf_long_ustory(:,:)      ! As above but for understory trees
     real(r8), allocatable :: root_long(:)               ! root turnover time (longevity) (pft)             [yr]
     real(r8), allocatable :: branch_long(:)             ! Turnover time for branchfall on live trees (pft) [yr]
     real(r8), allocatable :: turnover_nitr_retrans(:,:) ! nitrogen re-translocation fraction (pft x organ)
     real(r8), allocatable :: turnover_phos_retrans(:,:) ! phosphorus re-translocation fraction (pft x organ)
                                                         ! Parameters dimensioned by PFT and leaf age

     ! These vertical n profile scalers affect crown allometry and sla, thus they
     ! are here in the PRT module
     real(r8), allocatable :: leafn_vert_scaler_coeff1(:)  ! Coefficient one for decrease of leaf N through the canopy
     real(r8), allocatable :: leafn_vert_scaler_coeff2(:)  ! Coefficient two for decrease of leaf N through the canopy 
     
     real(r8), allocatable :: grperc(:)                  ! Growth respiration per unit Carbon gained
                                                         ! One value for whole plant
     ! ONLY parteh_mode == 1  [kg/kg]
     !     real(r8), allocatable ::grperc_organ(:,:)     ! Unit growth respiration (pft x organ) [kg/kg]
     !                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !                                                        ! THIS IS NOT READ IN BY THE PARAMETER FILE
     !                                                        ! THIS IS JUST FILLED BY GRPERC.  WE KEEP THIS
     !                                                        ! PARAMETER FOR HYPOTHESIS TESTING (ADVANCED USE)
     !                                                        ! IT HAS THE PRT_ TAG BECAUSE THIS PARAMETER
     !                                                        ! IS USED INSIDE PARTEH, WHILE GRPERC IS APPLIED
     !                                                        ! IN THE LEAF BIOPHYSICS SCHEME
     !                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     real(r8), allocatable :: nitr_stoich_p1(:,:)        ! Parameter 1 for nitrogen stoichiometry (pft x organ) 
     real(r8), allocatable :: phos_stoich_p1(:,:)        ! Parameter 1 for phosphorus stoichiometry (pft x organ) 

     real(r8), allocatable :: nitr_store_ratio(:)        ! This is the ratio of the target nitrogen stored per
                                                         ! target nitrogen that is bound into the tissues
                                                         ! of leaves, fine-roots and sapwood
     
     
     real(r8), allocatable :: phos_store_ratio(:)        ! This is the ratio of the target phosphorus stored per
                                                         ! target phosphorus is bound into the tissues
                                                         ! of leaves, fine-roots and sapwood

     integer, allocatable :: organ_id(:)                 ! Mapping of the organ index in the parameter file, to the
                                                         ! global list of organs found in PRTGenericMod.F90
     real(r8), allocatable :: alloc_priority(:,:)        ! Allocation priority for each organ (pft x organ) [integer 0-6]
     real(r8), allocatable :: cushion(:)                 ! labile carbon storage target as multiple of leaf pool.
     real(r8), allocatable :: leaf_stor_priority(:)      ! leaf turnover vs labile carbon use prioritisation
                                                         ! (1 = lose  leaves, 0 = use store).
     real(r8), allocatable :: dbh_repro_threshold(:)     ! diameter at which mature plants shift allocation
     real(r8), allocatable :: seed_alloc_mature(:)       ! fraction of carbon balance allocated to 
                                                         ! clonal reproduction.
     real(r8), allocatable :: seed_alloc(:)              ! fraction of carbon balance allocated to seeds.
     real(r8), allocatable :: repro_alloc_a(:)           ! ahb added this; sigmoidal shape param relating dbh to seed allocation fraction
     real(r8), allocatable :: repro_alloc_b(:)           ! ahb added this; intercept param relating dbh to seed allocation fraction

     ! Derived parameters

     integer, allocatable :: organ_param_id(:)           ! This is the sparse reverse lookup index map. This is dimensioned
                                                         ! by all the possible organs in parteh, and each index
                                                         ! may point to the index in the parameter file, or will be -1
     
     ! Allometry Parameters
     ! --------------------------------------------------------------------------------------------

     ! Root profile parameters. Note we have separate parameters for those that govern
     ! hydraulics, and those that govern biomass (for decomposition and respiration)

     real(r8), allocatable :: fnrt_prof_mode(:)             ! Fine root profile functional form
     real(r8), allocatable :: fnrt_prof_a(:)                ! Fine root profile scaling parameter A
     real(r8), allocatable :: fnrt_prof_b(:)                ! Fine root profile scaling parameter B

     real(r8), allocatable :: c2b(:)                        ! Carbon to biomass multiplier [kg/kgC]
     real(r8), allocatable :: wood_density(:)               ! wood density  g cm^-3  ...
     integer , allocatable :: woody(:)                      ! Does the plant have wood?      (1=yes, 0=no)
                                                            ! that is occupied by crown
     real(r8), allocatable :: slamax(:)                     ! Maximum specific leaf area of plant (at bottom) [m2/gC]
     real(r8), allocatable :: slatop(:)                     ! Specific leaf area at canopy top [m2/gC]
     real(r8), allocatable :: allom_sai_scaler(:)           ! 
     real(r8), allocatable :: allom_dbh_maxheight(:)        ! dbh at which height growth ceases
     integer , allocatable :: allom_hmode(:)                ! height allometry function type
     integer , allocatable :: allom_lmode(:)                ! maximum leaf allometry function type
     integer , allocatable :: allom_fmode(:)                ! maximum root allometry function type
     integer , allocatable :: allom_amode(:)                ! AGB allometry function type
     integer , allocatable :: allom_cmode(:)                ! Coarse root allometry function type
     integer , allocatable :: allom_smode(:)                ! sapwood allometry function type
     integer , allocatable :: allom_stmode(:)               ! storage allometry functional type 
                                                            !   0 - storage is proportional to maximum leaf biomass 
                                                            !       (considering trimmed)
                                                            !   1 - storage is proportional to maximum leaf biomass 
                                                            !       (untrimmed)
     integer , allocatable :: allom_dmode(:)                ! crown depth allometry function type
                                                            ! (HARD-CODED FOR TIME BEING, RGK 11-2017)
     real(r8), allocatable :: allom_la_per_sa_int(:)        ! Leaf area to sap area conversion, intercept 
                                                            ! (sapwood area / leaf area) [cm2/m2]
     real(r8), allocatable :: allom_la_per_sa_slp(:)        ! Leaf area to sap area conversion, slope 
                                                            ! (sapwood area / leaf area / diameter) [cm2/m2/cm]
     real(r8), allocatable :: allom_l2fr(:)                 ! Fine root biomass per leaf biomass ratio [kgC/kgC]
                                                            ! FOR C-ONLY: this is the static, unchanging ratio
                                                            ! FOR CNP: this is the initial value a cohort starts with
     real(r8), allocatable :: allom_agb_frac(:)             ! Fraction of stem above ground [-]
     real(r8), allocatable :: allom_d2h1(:)                 ! Parameter 1 for d2h allometry (intercept, or "c")
     real(r8), allocatable :: allom_d2h2(:)                 ! Parameter 2 for d2h allometry (slope, or "m")
     real(r8), allocatable :: allom_d2h3(:)                 ! Parameter 3 for d2h allometry (optional)
     real(r8), allocatable :: allom_d2bl1(:)                ! Parameter 1 for d2bl allometry (intercept)
     real(r8), allocatable :: allom_d2bl2(:)                ! Parameter 2 for d2bl allometry (slope)
     real(r8), allocatable :: allom_d2bl3(:)                ! Parameter 3 for d2bl allometry (optional)
     real(r8), allocatable :: allom_blca_expnt_diff(:)      ! Any difference in the exponent between the leaf
                                                            ! biomass and crown area scaling
     real(r8), allocatable :: allom_d2ca_coefficient_max(:) ! upper (savanna) value for crown 
                                                            ! area to dbh coefficient
     real(r8), allocatable :: allom_d2ca_coefficient_min(:) ! lower (closed-canopy forest) value for crown 
                                                            ! area to dbh coefficient
     real(r8), allocatable :: allom_agb1(:)                 ! Parameter 1 for agb allometry
     real(r8), allocatable :: allom_agb2(:)                 ! Parameter 2 for agb allometry
     real(r8), allocatable :: allom_agb3(:)                 ! Parameter 3 for agb allometry
     real(r8), allocatable :: allom_agb4(:)                 ! Parameter 3 for agb allometry

     real(r8), allocatable :: allom_h2cd1(:)                ! Parameter 1 for crown depth allometry. When allom_dmode=1
                                                            !    this is fraction of the height of the plant that is
                                                            !    considered crown (former parameter crown_depth_frac).
     real(r8), allocatable :: allom_h2cd2(:)                ! Exponent for crown depth allometry. Used only when
                                                            !    allom_dmode /= 1.

     real(r8), allocatable :: allom_zroot_max_dbh(:)        ! dbh at which maximum rooting depth saturates (largest possible) [cm]
     real(r8), allocatable :: allom_zroot_max_z(:)          ! the maximum rooting depth defined at dbh = fates_allom_zroot_max_dbh [m]
     real(r8), allocatable :: allom_zroot_min_dbh(:)        ! dbh at which the maximum rooting depth for a recruit is defined [cm]
     real(r8), allocatable :: allom_zroot_min_z(:)          ! the maximum rooting depth defined at dbh = fates_allom_zroot_min_dbh [m]
     real(r8), allocatable :: allom_zroot_k(:)              ! scale coefficient of logistic rooting depth model
     

     ! PID controller parameters
     real(r8), allocatable :: pid_kp(:)                     ! proportion constant in the PID controller for fine-root biomass
     real(r8), allocatable :: pid_ki(:)                     ! integral constant in the PID controller for fine-root biomass
     real(r8), allocatable :: pid_kd(:)                     ! derivative constant in the PID controller for fine-root biomass
     
     real(r8), allocatable :: store_ovrflw_frac(:)          ! For a coupled nutrient enabled simulation with dynamic fine-root biomass,
                                                            ! there will be an excess of at least two of the three species C, N or P.
                                                            ! This specifies how much excess (overflow) is allowed to be retained in storage
                                                            ! beyond the target level before it is either burned (C) or exuded (N or P). The
                                                            ! maximum value is the target * (1+store_ovrflw_frac)
     

     real(r8), allocatable :: nfix_mresp_scfrac(:)            ! Surcharge (as a fraction) to add to maintentance respiration
                                                            ! that is used to pay for N-Fixation
     
  end type prt_param_type

  type(prt_param_type),public :: prt_params          ! Instantiation of the parameter object


  

  
  
end module PRTParametersMod

