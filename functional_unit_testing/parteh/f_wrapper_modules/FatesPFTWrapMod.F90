! =======================================================================================
!
! This is the wrapper module that provides FATES data structures
!
! =======================================================================================

module EDPftvarcon
  
  use iso_c_binding, only : r8 => c_double
  use iso_c_binding, only : i4 => c_int
  use iso_c_binding, only : c_char
  
  implicit none
  private ! Modules are private by default

  integer,parameter,public :: SHR_KIND_CS = 80                     ! short char
  
  type, public :: EDPftvarcon_inst_type

     real(r8), pointer :: parteh_model(:)              ! The PARTEH model to use

     real(r8), pointer :: prescribed_npp_canopy(:)     ! this is only for the special 
                                                       ! prescribed_physiology_mode
     real(r8), pointer :: prescribed_npp_understory(:) ! this is only for the special 
                                                       ! prescribed_physiology_mode
     real(r8), pointer :: seed_alloc(:)
     real(r8), pointer :: seed_alloc_mature(:)
     real(r8), pointer :: dbh_repro_threshold(:)
     real(r8), pointer :: evergreen(:)
     real(r8), pointer :: season_decid(:)
     real(r8), pointer :: stress_decid(:)
     real(r8), pointer :: woody(:)
     real(r8), pointer :: hgt_min(:)
     real(r8), pointer :: allom_hmode(:)
     real(r8), pointer :: allom_amode(:)
     real(r8), pointer :: allom_lmode(:)
     real(r8), pointer :: allom_smode(:)
     real(r8), pointer :: allom_stmode(:)
     real(r8), pointer :: allom_cmode(:)
     real(r8), pointer :: allom_fmode(:)
     real(r8), pointer :: allom_d2h1(:)
     real(r8), pointer :: allom_d2h2(:)
     real(r8), pointer :: allom_d2h3(:)
     real(r8), pointer :: allom_dbh_maxheight(:)
     real(r8), pointer :: allom_agb1(:)
     real(r8), pointer :: allom_agb2(:)
     real(r8), pointer :: allom_agb3(:)
     real(r8), pointer :: allom_agb4(:)
     real(r8), pointer :: allom_d2bl1(:)
     real(r8), pointer :: allom_d2bl2(:)
     real(r8), pointer :: allom_d2bl3(:)
     real(r8), pointer :: wood_density(:)
     real(r8), pointer :: cushion(:)
     real(r8), pointer :: c2b(:)
     real(r8), pointer :: vcmax25top(:)
     real(r8), pointer :: allom_la_per_sa_int(:)
     real(r8), pointer :: allom_la_per_sa_slp(:)
     real(r8), pointer :: slatop(:)
     real(r8), pointer :: slamax(:)
     real(r8), pointer :: allom_l2fr(:)
     real(r8), pointer :: allom_agb_frac(:)
     real(r8), pointer :: allom_blca_expnt_diff(:)
     real(r8), pointer :: allom_d2ca_coefficient_min(:)
     real(r8), pointer :: allom_d2ca_coefficient_max(:)
     real(r8), pointer :: allom_sai_scaler(:)
     real(r8), pointer :: branch_turnover(:)
     real(r8), pointer :: leaf_long(:)
     real(r8), pointer :: root_long(:)
     real(r8), pointer :: leaf_stor_priority(:)
     real(r8), pointer :: roota_par(:)
     real(r8), pointer :: rootb_par(:)
     real(r8), pointer :: rootprof_beta(:,:)
     


     ! This array matches organ indices in the parameter file
     ! with global indices in PRTGeneric.  The basic global
     ! indices are leaf = 1
     !             fine-root = 2
     !             sapwood  = 3
     !             storage  = 4
     !             reproduction = 5
     !             structural = 6
     ! But, its possible that some organs may be added in
     ! the future, and then all hypotheses will not use the same
     ! set, or some hypotheses will sub-divide.

      ! These arrays hold the stoichiometric parameters
     ! The arrays are dimensioned by PFT X ORGAN
     ! Different formulations may use these parameters differently

     ! Hypothesis 1: Unused                        [na]


     real(r8), pointer :: prt_unit_gr_resp(:,:)
     real(r8), pointer :: prt_nitr_stoich_p1(:,:)
     real(r8), pointer :: prt_nitr_stoich_p2(:,:)
     real(r8), pointer :: prt_phos_stoich_p1(:,:)
     real(r8), pointer :: prt_phos_stoich_p2(:,:)
     real(r8), pointer :: prt_alloc_priority(:,:)
     
     ! THese are new, but not necessarily PARTEH labeled
     real(r8), pointer :: turnover_retrans_mode(:)
     
     real(r8), pointer :: turnover_carb_retrans(:,:)
     real(r8), pointer :: turnover_nitr_retrans(:,:)
     real(r8), pointer :: turnover_phos_retrans(:,:)


  end type EDPftvarcon_inst_type
  
  type, public :: pftptr_var
     real(r8), dimension(:),   pointer :: rp_1d
     real(r8), dimension(:,:), pointer :: rp_2d
     character(len=shr_kind_cs) :: var_name
  end type pftptr_var
  
  type, public :: EDPftvarcon_ptr_type
     type(pftptr_var), allocatable :: var(:)
  end type EDPftvarcon_ptr_type
  
  type(EDPftvarcon_inst_type), public :: EDPftvarcon_inst ! ED ecophysiological constants structure
  type(EDPftvarcon_ptr_type),  public :: EDPftvarcon_ptr  ! Pointer structure for obj-oriented id
  
  integer, public :: numparm     ! Number of different PFT parameters
  integer, public :: num_pft     ! Number of PFTs
  integer, public :: num_organs  ! Number of organs

  ! Make necessary procedures public 
  public :: EDPftvarconPySet
  public :: EDPftvarconAlloc
  
contains

  
  subroutine EDPftvarconPySet(ipft,i2d,rval,name)
    
    implicit none
    ! Arguments
    integer(i4),intent(in) :: ipft
    integer(i4),intent(in) :: i2d         ! Second dimension index
                                          ! if this is >0, use it
    character(kind=c_char,len=*), intent(in) :: name
    real(r8),intent(in) :: rval
    
    ! Locals
    logical :: npfound
    integer :: ip
    integer :: namelen
    
    namelen = len(trim(name))
    
    ip=0
    npfound = .false.
    do ip=1,numparm
       if (trim(name) == trim(EDPftvarcon_ptr%var(ip)%var_name ) ) then
          if(i2d==0) then
             EDPftvarcon_ptr%var(ip)%rp_1d(ipft) = rval
          else
             EDPftvarcon_ptr%var(ip)%rp_2d(ipft,i2d) = rval
          end if
          npfound = .true.
       end if
    end do
    
    if(.not.npfound)then
       print*,"Could not find parameter passed in from python driver"
       print*,"registerred in the fortran wrapper"
       print*,"--",trim(name),"--"
       stop
    end if
    
    ! Performa a check to see if the target array is being filled
    
    if (trim(name) == 'fates_wood_density' ) then
       if (EDPftvarcon_inst%wood_density(ipft) .ne. rval) then
          print*,"F90: POINTER CHECK FAILS:",rval," != ",EDPftvarcon_inst%wood_density(ipft)
          stop
       end if
    end if
    
    return
  end subroutine EDPftvarconPySet
  
  ! ====================================================================================
  
  subroutine EDPftvarconAlloc(numpft_in, numorgans_in)
    
    ! !ARGUMENTS:
    integer(i4), intent(in) :: numpft_in
    integer(i4), intent(in) :: numorgans_in

    ! LOCALS:
    integer                    :: iv   ! The parameter incrementer
    integer, parameter         :: n_beta_dims = 1
    !------------------------------------------------------------------------

    num_pft    = numpft_in
    num_organs = numorgans_in

    allocate( EDPftvarcon_ptr%var (100) ) ! Make this plenty large

    iv=0

    allocate( EDPftvarcon_inst%parteh_model(1:num_pft)); 
    EDPftvarcon_inst%parteh_model (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "parteh_model"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%parteh_model


    allocate( EDPftvarcon_inst%dbh_repro_threshold(1:num_pft)); 
    EDPftvarcon_inst%dbh_repro_threshold (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_dbh_repro_threshold"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%dbh_repro_threshold


    allocate( EDPftvarcon_inst%prescribed_npp_canopy(1:num_pft)); 
    EDPftvarcon_inst%prescribed_npp_canopy (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prescribed_npp_canopy"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%prescribed_npp_canopy

    allocate( EDPftvarcon_inst%prescribed_npp_understory(1:num_pft)); 
    EDPftvarcon_inst%prescribed_npp_understory (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prescribed_npp_understory"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%prescribed_npp_understory

    allocate( EDPftvarcon_inst%seed_alloc(1:num_pft)); 
    EDPftvarcon_inst%seed_alloc (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_seed_alloc"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%seed_alloc
    

    allocate( EDPftvarcon_inst%seed_alloc_mature(1:num_pft)); 
    EDPftvarcon_inst%seed_alloc_mature(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_seed_alloc_mature"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%seed_alloc_mature

    allocate( EDPftvarcon_inst%evergreen(1:num_pft)); 
    EDPftvarcon_inst%evergreen (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_phen_evergreen"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%evergreen

    allocate( EDPftvarcon_inst%season_decid(1:num_pft)); 
    EDPftvarcon_inst%season_decid (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_phen_season_decid"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%season_decid

    allocate( EDPftvarcon_inst%stress_decid(1:num_pft)); 
    EDPftvarcon_inst%stress_decid (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_phen_stress_decid"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%stress_decid


    allocate( EDPftvarcon_inst%woody(1:num_pft)); 
    EDPftvarcon_inst%woody (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_woody"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%woody

    allocate( EDPftvarcon_inst%hgt_min(1:num_pft)); 
    EDPftvarcon_inst%hgt_min (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_recruit_hgt_min"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%hgt_min
    
    allocate( EDPftvarcon_inst%allom_dbh_maxheight(1:num_pft)); 
    EDPftvarcon_inst%allom_dbh_maxheight (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_dbh_maxheight"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_dbh_maxheight

    allocate( EDPftvarcon_inst%allom_hmode(1:num_pft)); 
    EDPftvarcon_inst%allom_hmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_hmode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_hmode
    
    allocate( EDPftvarcon_inst%allom_amode(1:num_pft)); 
    EDPftvarcon_inst%allom_amode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_amode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_amode
    
    allocate( EDPftvarcon_inst%allom_lmode(1:num_pft)); 
    EDPftvarcon_inst%allom_lmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_lmode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_lmode

    allocate( EDPftvarcon_inst%allom_smode(1:num_pft)); 
    EDPftvarcon_inst%allom_smode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_smode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_smode

    allocate( EDPftvarcon_inst%allom_stmode(1:num_pft)); 
    EDPftvarcon_inst%allom_stmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_stmode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_stmode

    allocate( EDPftvarcon_inst%allom_cmode(1:num_pft)); 
    EDPftvarcon_inst%allom_cmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_cmode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_cmode

    allocate( EDPftvarcon_inst%allom_fmode(1:num_pft)); 
    EDPftvarcon_inst%allom_fmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_fmode"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_fmode
     
    allocate( EDPftvarcon_inst%allom_d2h1(1:num_pft)); 
    EDPftvarcon_inst%allom_d2h1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h1"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2h1

    allocate( EDPftvarcon_inst%allom_d2h2(1:num_pft)); 
    EDPftvarcon_inst%allom_d2h2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h2"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2h2

    allocate( EDPftvarcon_inst%allom_d2h3(1:num_pft)); 
    EDPftvarcon_inst%allom_d2h3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h3"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2h3
    
    allocate( EDPftvarcon_inst%allom_agb1(1:num_pft)); 
    EDPftvarcon_inst%allom_agb1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb1"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_agb1

    allocate( EDPftvarcon_inst%allom_agb2(1:num_pft)); 
    EDPftvarcon_inst%allom_agb2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb2"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_agb2

    allocate( EDPftvarcon_inst%allom_agb3(1:num_pft)); 
    EDPftvarcon_inst%allom_agb3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb3"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_agb3

    allocate( EDPftvarcon_inst%allom_agb4(1:num_pft)); 
    EDPftvarcon_inst%allom_agb4(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb4"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_agb4

    allocate( EDPftvarcon_inst%allom_d2bl1(1:num_pft)); 
    EDPftvarcon_inst%allom_d2bl1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl1"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2bl1

    allocate( EDPftvarcon_inst%allom_d2bl2(1:num_pft)); 
    EDPftvarcon_inst%allom_d2bl2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl2"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2bl2

    allocate( EDPftvarcon_inst%allom_d2bl3(1:num_pft)); 
    EDPftvarcon_inst%allom_d2bl3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl3"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2bl3

    allocate( EDPftvarcon_inst%cushion(1:num_pft)); 
    EDPftvarcon_inst%cushion(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_cushion"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%cushion

    allocate( EDPftvarcon_inst%wood_density(1:num_pft)); 
    EDPftvarcon_inst%wood_density(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_wood_density"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%wood_density

    allocate( EDPftvarcon_inst%c2b(1:num_pft)); 
    EDPftvarcon_inst%c2b(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_c2b"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%c2b

    allocate( EDPftvarcon_inst%vcmax25top(1:num_pft)); 
    EDPftvarcon_inst%vcmax25top(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_vcmax25top"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%vcmax25top

    allocate( EDPftvarcon_inst%allom_la_per_sa_int(1:num_pft)); 
    EDPftvarcon_inst%allom_la_per_sa_int(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_la_per_sa_int"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_la_per_sa_int

    allocate( EDPftvarcon_inst%allom_la_per_sa_slp(1:num_pft)); 
    EDPftvarcon_inst%allom_la_per_sa_slp(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_la_per_sa_slp"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_la_per_sa_slp

    allocate( EDPftvarcon_inst%slatop(1:num_pft)); 
    EDPftvarcon_inst%slatop(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_slatop"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%slatop

    allocate( EDPftvarcon_inst%slamax(1:num_pft)); 
    EDPftvarcon_inst%slamax(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_slamax"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%slamax
    
    
    allocate( EDPftvarcon_inst%allom_l2fr(1:num_pft)); 
    EDPftvarcon_inst%allom_l2fr(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_l2fr"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_l2fr

    allocate( EDPftvarcon_inst%allom_agb_frac(1:num_pft)); 
    EDPftvarcon_inst%allom_agb_frac(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb_frac"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_agb_frac
    
    allocate( EDPftvarcon_inst%allom_sai_scaler(1:num_pft)); 
    EDPftvarcon_inst%allom_sai_scaler(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_sai_scaler"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_sai_scaler

    allocate( EDPftvarcon_inst%allom_blca_expnt_diff(1:num_pft)); 
    EDPftvarcon_inst%allom_blca_expnt_diff(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_blca_expnt_diff"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_blca_expnt_diff

    allocate( EDPftvarcon_inst%allom_d2ca_coefficient_min(1:num_pft)); 
    EDPftvarcon_inst%allom_d2ca_coefficient_min(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2ca_coefficient_min"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2ca_coefficient_min

    allocate( EDPftvarcon_inst%allom_d2ca_coefficient_max(1:num_pft)); 
    EDPftvarcon_inst%allom_d2ca_coefficient_max(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2ca_coefficient_max"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%allom_d2ca_coefficient_max

    allocate( EDPftvarcon_inst%branch_turnover(1:num_pft)); 
    EDPftvarcon_inst%branch_turnover(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_branch_turnover"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%branch_turnover

    allocate( EDPftvarcon_inst%leaf_long(1:num_pft)); 
    EDPftvarcon_inst%leaf_long(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_leaf_long"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%leaf_long
    
    allocate( EDPftvarcon_inst%root_long(1:num_pft)); 
    EDPftvarcon_inst%root_long(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_root_long"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%root_long
    
    allocate( EDPftvarcon_inst%leaf_stor_priority(1:num_pft)); 
    EDPftvarcon_inst%leaf_stor_priority(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_leaf_stor_priority"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%leaf_stor_priority

    allocate( EDPftvarcon_inst%roota_par(1:num_pft)); 
    EDPftvarcon_inst%roota_par(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_roota_par"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%roota_par

    allocate( EDPftvarcon_inst%rootb_par(1:num_pft)); 
    EDPftvarcon_inst%rootb_par(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_rootb_par"
    EDPftvarcon_ptr%var(iv)%rp_1d   => EDPftvarcon_inst%rootb_par


    allocate( EDPftvarcon_inst%prt_nitr_stoich_p1(1:num_pft,1:num_organs)); 
    EDPftvarcon_inst%prt_nitr_stoich_p1(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_nitr_stoich_p1"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_nitr_stoich_p1

    
    allocate( EDPftvarcon_inst%prt_phos_stoich_p1(1:num_pft,1:num_organs)); 
    EDPftvarcon_inst%prt_phos_stoich_p1(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_phos_stoich_p1"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_phos_stoich_p1


    allocate( EDPftvarcon_inst%prt_nitr_stoich_p2(1:num_pft,1:num_organs)); 
    EDPftvarcon_inst%prt_nitr_stoich_p2(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_nitr_stoich_p2"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_nitr_stoich_p2

    
    allocate( EDPftvarcon_inst%prt_phos_stoich_p2(1:num_pft,1:num_organs)); 
    EDPftvarcon_inst%prt_phos_stoich_p2(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_phos_stoich_p2"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_phos_stoich_p2


    allocate( EDPftvarcon_inst%prt_unit_gr_resp(1:num_pft,1:num_organs));
    EDPftvarcon_inst%prt_unit_gr_resp(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_unit_gr_resp"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_unit_gr_resp


    allocate( EDPftvarcon_inst%prt_alloc_priority(1:num_pft,1:num_organs));
    EDPftvarcon_inst%prt_alloc_priority(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_prt_alloc_priority"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%prt_alloc_priority
    
    allocate( EDPftvarcon_inst%turnover_retrans_mode(1:num_pft) )
    EDPftvarcon_inst%turnover_retrans_mode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_turnover_retrans_mode"
    EDPftvarcon_ptr%var(iv)%rp_1d => EDPftvarcon_inst%turnover_retrans_mode

    allocate( EDPftvarcon_inst%turnover_carb_retrans(1:num_pft,1:num_organs) )
    EDPftvarcon_inst%turnover_carb_retrans(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_turnover_carb_retrans"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%turnover_carb_retrans

    allocate( EDPftvarcon_inst%turnover_nitr_retrans(1:num_pft,1:num_organs) )
    EDPftvarcon_inst%turnover_nitr_retrans(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_turnover_nitr_retrans"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%turnover_nitr_retrans
    
    allocate( EDPftvarcon_inst%turnover_phos_retrans(1:num_pft,1:num_organs) )
    EDPftvarcon_inst%turnover_phos_retrans(:,:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_turnover_phos_retrans"
    EDPftvarcon_ptr%var(iv)%rp_2d => EDPftvarcon_inst%turnover_phos_retrans


    ! We should gracefully fail if rootprof_beta is requested
    allocate( EDPftvarcon_inst%rootprof_beta(1:num_pft,n_beta_dims)); 
    EDPftvarcon_inst%rootprof_beta(:,:) = nan

    
    numparm = iv

    print*,"F90: ALLOCATED ",numparm," PARAMETERS, FOR ",num_pft," PFTs"
    
    
    return
  end subroutine EDPftvarconAlloc
  
end module EDPftvarcon
