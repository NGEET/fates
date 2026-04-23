module FatesParametersInterface

  ! This module simply holds the instantiation
  ! of the generalized FATES parameter file that
  ! is provide by the JSON parser. ie pstruct
  ! This module is intentionaly low dependency
  ! Note also that JSONParameterUtilsMod
  ! only uses shr libraries

  use JSONParameterUtilsMod, only: params_type
  use JSONParameterUtilsMod, only: param_type
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  private
  
  type(params_type) :: pstruct

  ! Parameter indexes
  
  integer, public :: pid_vcmax25top
  integer, public :: pid_cwd_frac
  integer, public :: pid_damage_frac
  integer, public :: pid_damage_bins
  integer, public :: pid_alloc_organ_id
  integer, public :: pid_phen_leaf_habit
  integer, public :: pid_phen_stem_drop_frac
  integer, public :: pid_phen_fnrt_drop_frac
  integer, public :: pid_phen_mindaysoff
  integer, public :: pid_phen_drought_thresh
  integer, public :: pid_phen_moist_thresh
  integer, public :: pid_leaf_slamax
  integer, public :: pid_leaf_slatop
  integer, public :: pid_allom_sai_scaler
  integer, public :: pid_allom_fnrt_prof_a
  integer, public :: pid_allom_fnrt_prof_b
  integer, public :: pid_allom_fnrt_prof_mode
  integer, public :: pid_woody
  integer, public :: pid_wood_density
  integer, public :: pid_recruit_seed_dbh_thresh
  integer, public :: pid_storage_cushion
  integer, public :: pid_store_priority_frac
  integer, public :: pid_turnover_senleaf_fdrought
  integer, public :: pid_turnover_fnrt
  integer, public :: pid_leafn_vert_scaler1
  integer, public :: pid_leafn_vert_scaler2
  integer, public :: pid_recruit_seed_alloc_mature
  integer, public :: pid_recruit_seed_alloc
  integer, public :: pid_trs_repro_alloc_a
  integer, public :: pid_trs_repro_alloc_b
  integer, public :: pid_c2b
  integer, public :: pid_grperc
  integer, public :: pid_allom_dbh_maxh
  integer, public :: pid_allom_hmode
  integer, public :: pid_allom_lmode
  integer, public :: pid_allom_fmode
  integer, public :: pid_allom_amode
  integer, public :: pid_allom_stmode
  integer, public :: pid_allom_cmode
  integer, public :: pid_allom_smode
  integer, public :: pid_allom_dmode
  integer, public :: pid_allom_lapersa_int
  integer, public :: pid_allom_lapersa_slp
  integer, public :: pid_allom_l2fr
  integer, public :: pid_cnp_kp
  integer, public :: pid_cnp_ki
  integer, public :: pid_cnp_kd
  integer, public :: pid_cnp_store_ovrflw_frac
  integer, public :: pid_cnp_nfix1
  integer, public :: pid_allom_d2h1
  integer, public :: pid_allom_d2h2
  integer, public :: pid_allom_d2h3
  integer, public :: pid_allom_agbfrac
  integer, public :: pid_allom_d2l1
  integer, public :: pid_allom_d2l2
  integer, public :: pid_allom_d2l3
  integer, public :: pid_allom_blca_ediff
  integer, public :: pid_allom_d2ca_max
  integer, public :: pid_allom_d2ca_min
  integer, public :: pid_allom_agb1
  integer, public :: pid_allom_agb2
  integer, public :: pid_allom_agb3
  integer, public :: pid_allom_agb4
  integer, public :: pid_allom_h2cd1
  integer, public :: pid_allom_h2cd2
  integer, public :: pid_allom_zroot_maxd
  integer, public :: pid_allom_zroot_maxz
  integer, public :: pid_allom_zroot_mind
  integer, public :: pid_allom_zroot_minz
  integer, public :: pid_allom_zroot_k
  integer, public :: pid_turnover_branch
  integer, public :: pid_cnp_nitr_store_ratio
  integer, public :: pid_cnp_phos_store_ratio
  integer, public :: pid_turnover_leaf_canopy
  integer, public :: pid_turnover_leaf_ustory
  integer, public :: pid_stoich_nitr
  integer, public :: pid_stoich_phos
  integer, public :: pid_alloc_organ_priority
  integer, public :: pid_cnp_nitr_retrans
  integer, public :: pid_cnp_phos_retrans
  
  public :: pstruct
  public :: GetParameterIndices
  public :: Transp2dInt
  public :: Transp2dReal

contains

  subroutine GetParameterIndices()

    ! Query the parameter data structure for the names
    ! of known parameters, and identify the data-structure
    ! index of each, into a named integer constant. These
    ! integer constants will be used for retrieval.
    ! This assignment happens once, here.

    ! Scalar Parameters
    pid_damage_bins = pstruct%GetIndexFromName('fates_history_damage_bin_edges')
    
    pid_vcmax25top = pstruct%GetIndexFromName('fates_leaf_vcmax25top')
    pid_cwd_frac = pstruct%GetIndexFromName('fates_frag_cwd_frac')
    pid_damage_frac = pstruct%GetIndexFromName('fates_damage_frac')
    
    pid_alloc_organ_id = pstruct%GetIndexFromName('fates_alloc_organ_id')
    pid_phen_leaf_habit = pstruct%GetIndexFromName('fates_phen_leaf_habit')
    pid_phen_stem_drop_frac = pstruct%GetIndexFromName('fates_phen_stem_drop_fraction')
    pid_phen_fnrt_drop_frac = pstruct%GetIndexFromName('fates_phen_fnrt_drop_fraction')
    pid_phen_mindaysoff = pstruct%GetIndexFromName('fates_phen_mindaysoff')
    pid_phen_drought_thresh = pstruct%GetIndexFromName('fates_phen_drought_threshold')
    pid_phen_moist_thresh = pstruct%GetIndexFromName('fates_phen_moist_threshold')
    pid_leaf_slamax = pstruct%GetIndexFromName('fates_leaf_slamax')
    pid_leaf_slatop = pstruct%GetIndexFromName('fates_leaf_slatop')
    pid_allom_sai_scaler = pstruct%GetIndexFromName('fates_allom_sai_scaler')
    pid_allom_fnrt_prof_a = pstruct%GetIndexFromName('fates_allom_fnrt_prof_a')
    pid_allom_fnrt_prof_b = pstruct%GetIndexFromName('fates_allom_fnrt_prof_b')
    pid_allom_fnrt_prof_mode = pstruct%GetIndexFromName('fates_allom_fnrt_prof_mode')
    pid_woody = pstruct%GetIndexFromName('fates_woody')
    pid_wood_density = pstruct%GetIndexFromName('fates_wood_density')
    pid_recruit_seed_dbh_thresh = pstruct%GetIndexFromName('fates_recruit_seed_dbh_repro_threshold')
    pid_storage_cushion = pstruct%GetIndexFromName('fates_alloc_storage_cushion')
    pid_store_priority_frac = pstruct%GetIndexFromName('fates_alloc_store_priority_frac')
    pid_turnover_senleaf_fdrought = pstruct%GetIndexFromName('fates_turnover_senleaf_fdrought')
    pid_turnover_fnrt = pstruct%GetIndexFromName('fates_turnover_fnrt')
    pid_leafn_vert_scaler1 = pstruct%GetIndexFromName('fates_leafn_vert_scaler_coeff1')
    pid_leafn_vert_scaler2 = pstruct%GetIndexFromName('fates_leafn_vert_scaler_coeff2')
    pid_recruit_seed_alloc_mature = pstruct%GetIndexFromName('fates_recruit_seed_alloc_mature')
    pid_recruit_seed_alloc = pstruct%GetIndexFromName('fates_recruit_seed_alloc')
    pid_trs_repro_alloc_a = pstruct%GetIndexFromName('fates_trs_repro_alloc_a')
    pid_trs_repro_alloc_b = pstruct%GetIndexFromName('fates_trs_repro_alloc_b')
    pid_c2b = pstruct%GetIndexFromName('fates_c2b')
    pid_grperc = pstruct%GetIndexFromName('fates_grperc')
    pid_allom_dbh_maxh = pstruct%GetIndexFromName('fates_allom_dbh_maxheight')
    pid_allom_hmode = pstruct%GetIndexFromName('fates_allom_hmode')
    pid_allom_lmode = pstruct%GetIndexFromName('fates_allom_lmode')
    pid_allom_fmode = pstruct%GetIndexFromName('fates_allom_fmode')
    pid_allom_amode = pstruct%GetIndexFromName('fates_allom_amode')
    pid_allom_stmode = pstruct%GetIndexFromName('fates_allom_stmode')
    pid_allom_cmode = pstruct%GetIndexFromName('fates_allom_cmode')
    pid_allom_smode = pstruct%GetIndexFromName('fates_allom_smode')
    pid_allom_dmode = pstruct%GetIndexFromName('fates_allom_dmode')
    pid_allom_lapersa_int = pstruct%GetIndexFromName('fates_allom_la_per_sa_int')
    pid_allom_lapersa_slp = pstruct%GetIndexFromName('fates_allom_la_per_sa_slp')
    pid_allom_l2fr = pstruct%GetIndexFromName('fates_allom_l2fr')
    pid_cnp_kp = pstruct%GetIndexFromName('fates_cnp_pid_kp')
    pid_cnp_ki = pstruct%GetIndexFromName('fates_cnp_pid_ki')
    pid_cnp_kd = pstruct%GetIndexFromName('fates_cnp_pid_kd')
    pid_cnp_store_ovrflw_frac = pstruct%GetIndexFromName('fates_cnp_store_ovrflw_frac')
    pid_cnp_nfix1 = pstruct%GetIndexFromName('fates_cnp_nfix1')
    pid_allom_agbfrac = pstruct%GetIndexFromName('fates_allom_agb_frac')
    pid_allom_d2h1 = pstruct%GetIndexFromName('fates_allom_d2h1')
    pid_allom_d2h2 = pstruct%GetIndexFromName('fates_allom_d2h2')
    pid_allom_d2h3 = pstruct%GetIndexFromName('fates_allom_d2h3')
    pid_allom_d2l1 = pstruct%GetIndexFromName('fates_allom_d2bl1')
    pid_allom_d2l2 = pstruct%GetIndexFromName('fates_allom_d2bl2')
    pid_allom_d2l3 = pstruct%GetIndexFromName('fates_allom_d2bl3')
    pid_allom_blca_ediff = pstruct%GetIndexFromName('fates_allom_blca_expnt_diff')
    pid_allom_d2ca_max = pstruct%GetIndexFromName('fates_allom_d2ca_coefficient_max')
    pid_allom_d2ca_min = pstruct%GetIndexFromName('fates_allom_d2ca_coefficient_min')
    pid_allom_agb1 = pstruct%GetIndexFromName('fates_allom_agb1')
    pid_allom_agb2 = pstruct%GetIndexFromName('fates_allom_agb2')
    pid_allom_agb3 = pstruct%GetIndexFromName('fates_allom_agb3')
    pid_allom_agb4 = pstruct%GetIndexFromName('fates_allom_agb4')
    pid_allom_h2cd1 = pstruct%GetIndexFromName('fates_allom_h2cd1')
    pid_allom_h2cd2 = pstruct%GetIndexFromName('fates_allom_h2cd2')
    pid_allom_zroot_maxd = pstruct%GetIndexFromName('fates_allom_zroot_max_dbh')
    pid_allom_zroot_maxz = pstruct%GetIndexFromName('fates_allom_zroot_max_z')
    pid_allom_zroot_mind = pstruct%GetIndexFromName('fates_allom_zroot_min_dbh')
    pid_allom_zroot_minz = pstruct%GetIndexFromName('fates_allom_zroot_min_z')
    pid_allom_zroot_k = pstruct%GetIndexFromName('fates_allom_zroot_k')
    pid_turnover_branch = pstruct%GetIndexFromName('fates_turnover_branch')
    pid_cnp_nitr_store_ratio = pstruct%GetIndexFromName('fates_cnp_nitr_store_ratio')
    pid_cnp_phos_store_ratio = pstruct%GetIndexFromName('fates_cnp_phos_store_ratio')
    pid_turnover_leaf_canopy = pstruct%GetIndexFromName('fates_turnover_leaf_canopy')
    pid_turnover_leaf_ustory = pstruct%GetIndexFromName('fates_turnover_leaf_ustory')
    pid_stoich_nitr = pstruct%GetIndexFromName('fates_stoich_nitr')
    pid_stoich_phos = pstruct%GetIndexFromName('fates_stoich_phos')
    pid_alloc_organ_priority = pstruct%GetIndexFromName('fates_alloc_organ_priority')
    pid_cnp_nitr_retrans = pstruct%GetIndexFromName('fates_cnp_turnover_nitr_retrans')
    pid_cnp_phos_retrans = pstruct%GetIndexFromName('fates_cnp_turnover_phos_retrans')

  end subroutine GetParameterIndices

  subroutine CheckParameters

    type(param_type), pointer :: param(:)
    integer :: corr_id
    real(r8) :: correction
    
    param => pstruct%parameters
    
    ! Apply correction to cwd fractions
    correction = 1._r8 - sum(param(pid_cwd_frac)%r_data_1d(:),dim=1)
    corr_id = maxloc(param(pid_cwd_frac)%r_data_1d(:),dim=1)
    param(pid_cwd_frac)%r_data_1d(corr_id) = param(pid_cwd_frac)%r_data_1d(corr_id) + correction

    
  end subroutine CheckParameters

  
  subroutine Transp2dInt(i_2d_in,i_2d_out)

    ! The FATES JSON parameter files have a legacy from the netcdf 
    ! of having the pfts as the column (inner) indices. At least
    ! that is how it is presented in the text files.
    ! [ [pft1, pft2, pft3],[pf1, pft2, pft3] ]
    !
    ! In fortran, this would have the pft index as the second index
    ! because it goes row x column.
    !
    ! However, we the data arrays in fates use
    ! the PFT as the first dimension... so:
    
    integer :: i_2d_in(:,:)
    integer :: i_2d_out(:,:)
    integer :: i,j

    do i=1,size(i_2d_in,dim=1)
       do j = 1,size(i_2d_in,dim=2)
          i_2d_out(j,i) = i_2d_in(i,j)
       end do
    end do
    
  end subroutine Transp2dInt
  
  subroutine Transp2dReal(r_2d_in,r_2d_out)

    ! The FATES JSON parameter files have a legacy from the netcdf 
    ! of having the pfts as the column (inner) indices. At least
    ! that is how it is presented in the text files.
    ! [ [pft1, pft2, pft3],[pf1, pft2, pft3] ]
    !
    ! In fortran, this would have the pft index as the second index
    ! because it goes row x column.
    !
    ! However, we the data arrays in fates use
    ! the PFT as the first dimension... so:
    
    real(r8) :: r_2d_in(:,:)
    real(r8) :: r_2d_out(:,:)
    integer :: i,j
    
    do i=1,size(r_2d_in,dim=1)
       do j = 1,size(r_2d_in,dim=2)
          r_2d_out(j,i) = r_2d_in(i,j)
       end do
    end do

  end subroutine Transp2dReal
  
end module FatesParametersInterface

