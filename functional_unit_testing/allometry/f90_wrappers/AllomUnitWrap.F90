
! =======================================================================================
!
! This file is an alternative to key files in the fates 
! filesystem. Noteably, we replace fates_r8 and fates_in
! with types that work with "ctypes".  This is
! a key step in working with python
! 
! We also wrap FatesGlobals to reduce the dependancy
! cascade that it pulls in from shr_log_mod.
!
! =======================================================================================

module shr_log_mod

   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int

   contains

   function shr_log_errMsg(source, line) result(ans)
      character(kind=c_char,len=*), intent(in) :: source
      integer(c_int), intent(in) :: line
      character(kind=c_char,len=128) :: ans

      ans = "source: " // trim(source) // " line: "
   end function shr_log_errMsg
   
end module shr_log_mod


module FatesGlobals

   contains

   integer function fates_log()
      fates_log = -1
   end function fates_log
   
   subroutine fates_endrun(msg) 

      implicit none
      character(len=*), intent(in) :: msg    ! string to be printed
      
      stop

   end subroutine fates_endrun

end module FatesGlobals


module EDTypesMod

   use iso_c_binding, only : r8 => c_double

  integer, parameter  :: nclmax = 2
  integer, parameter  :: nlevleaf = 30
  real(r8), parameter :: dinc_ed = 1.0_r8

end module EDTypesMod


module EDPftvarcon
   
   use iso_c_binding, only : r8 => c_double
   use iso_c_binding, only : i4 => c_int
   use iso_c_binding, only : c_char

   integer,parameter :: SHR_KIND_CS = 80                     ! short char

   type, public :: EDPftvarcon_inst_type

   ! VARIABLE-DEFINITIONS-HERE (DO NOT REMOVE THIS LINE, OR MOVE IT)

    real(r8),pointer :: allom_d2h1(:)
    real(r8),pointer :: allom_d2h2(:)
    real(r8),pointer :: allom_d2h3(:)
    real(r8),pointer :: allom_hmode(:)
    real(r8),pointer :: allom_dbh_maxheight(:)
    real(r8),pointer :: allom_agb1(:)
    real(r8),pointer :: allom_agb2(:)
    real(r8),pointer :: allom_agb3(:)
    real(r8),pointer :: allom_agb4(:)
    real(r8),pointer :: wood_density(:)
    real(r8),pointer :: c2b(:)
    real(r8),pointer :: allom_agb_frac(:)
    real(r8),pointer :: allom_amode(:)
    real(r8),pointer :: allom_lmode(:)
    real(r8),pointer :: allom_d2bl1(:)
    real(r8),pointer :: allom_d2bl2(:)
    real(r8),pointer :: allom_d2bl3(:)
    real(r8),pointer :: allom_blca_expnt_diff(:)
    real(r8),pointer :: allom_d2ca_coefficient_min(:)
    real(r8),pointer :: allom_d2ca_coefficient_max(:)
    real(r8),pointer :: slatop(:)
    real(r8),pointer :: slamax(:)
    real(r8),pointer :: allom_sai_scaler(:)
    real(r8),pointer :: allom_smode(:)
    real(r8),pointer :: allom_cmode(:)
    real(r8),pointer :: allom_fmode(:)
    real(r8),pointer :: allom_stmode(:)
    real(r8),pointer :: cushion(:)
    real(r8),pointer :: allom_l2fr(:)
    real(r8),pointer :: allom_la_per_sa_int(:)
    real(r8),pointer :: allom_la_per_sa_slp(:)
    real(r8),pointer :: roota_par(:)
    real(r8),pointer :: rootb_par(:)
    real(r8),pointer :: rootprof_beta(:,:)
    real(r8),pointer :: woody(:)
   end type EDPftvarcon_inst_type
 
  type ptr_var1
     real(r8), dimension(:), pointer :: var_rp
     integer(i4), dimension(:), pointer :: var_ip
     character(len=shr_kind_cs) :: var_name
     integer :: vtype
  end type ptr_var1

  type ptr_var2
     real(r8), dimension(:,:), pointer :: var_rp
     integer(i4), dimension(:,:), pointer :: var_ip
     character(len=shr_kind_cs) :: var_name
     integer :: vtype
  end type ptr_var2

  type EDPftvarcon_ptr_type
     type(ptr_var1), allocatable :: var1d(:)
	 type(ptr_var2), allocatable :: var2d(:)
  end type EDPftvarcon_ptr_type
  

  type(EDPftvarcon_inst_type), public :: EDPftvarcon_inst ! ED ecophysiological constants structure
  type(EDPftvarcon_ptr_type),  public :: EDPftvarcon_ptr  ! Pointer structure for obj-oriented id
  
  integer :: numparm1d ! Number of different PFT parameters
  integer :: numparm2d
  integer :: numpft

  logical, parameter ::  debug = .true.

contains
  

   subroutine EDPftvarconPySet(ipft,rval,ival,name)

      implicit none
      ! Arguments
      integer(i4),intent(in) :: ipft
      character(kind=c_char,len=*), intent(in) :: name
      real(r8),intent(in) :: rval
      integer(i4),intent(in) :: ival
      ! Locals
      logical :: npfound
      integer :: ip
      integer :: namelen
      
      namelen = len(trim(name))

	  if(debug) print*,"F90: ARGS: ",trim(name)," IPFT: ",ipft," RVAL: ",rval," IVAL: ",ival

      ip=0
      npfound = .true.
      do ip=1,numparm1d

         if (trim(name) == trim(EDPftvarcon_ptr%var1d(ip)%var_name ) ) then
            print*,"F90: Found ",trim(name)," in lookup table"
            npfound = .false.
            if(EDPftvarcon_ptr%var1d(ip)%vtype == 1) then ! real
               EDPftvarcon_ptr%var1d(ip)%var_rp(ipft) = rval
            elseif(EDPftvarcon_ptr%var1d(ip)%vtype == 2) then ! integer
               EDPftvarcon_ptr%var1d(ip)%var_ip(ipft) = ival
            else
               print*,"F90: STRANGE TYPE"
               stop
            end if
         end if
      end do

      if(npfound)then
         print*,"F90: The parameter you loaded DNE: ",name(:)
         stop
      end if

      do ip=1,numparm2d
         if (trim(name) == trim(EDPftvarcon_ptr%var2d(ip)%var_name)) then
            print*,"F90: Found ",trim(name)," in lookup table"
			print*,"BUT... WE AVOID USING 2D VARIABLES FOR NOW..."
			print*,"REMOVE THIS TEST"
            stop
         end if
      end do


	  ! Perform a check to see if the target array is being filled
      if (trim(name) == 'fates_allom_d2h1') then
         if (EDPftvarcon_inst%allom_d2h1(ipft) == rval) then
            print*,"F90: POINTER CHECK PASSES:",rval," = ",EDPftvarcon_inst%allom_d2h1(ipft)
         else
            print*,"F90: POINTER CHECK FAILS:",rval," != ",EDPftvarcon_inst%allom_d2h1(ipft)
            stop
         end if
      end if

      if (trim(name) == 'fates_wood_density' ) then
         if (EDPftvarcon_inst%wood_density(ipft) == rval) then
            print*,"F90: POINTER CHECK PASSES:",rval," = ",EDPftvarcon_inst%wood_density(ipft)
         else
            print*,"F90: POINTER CHECK FAILS:",rval," != ",EDPftvarcon_inst%wood_density(ipft)
            stop
         end if
      end if

      return
   end subroutine EDPftvarconPySet


  subroutine EDPftvarconAlloc(numpft_in)
    !

    ! !ARGUMENTS:
    integer(i4), intent(in) :: numpft_in
    ! LOCALS:
    integer                    :: iv1   ! The parameter incrementer
	integer                    :: iv2  
    !------------------------------------------------------------------------

    numpft = numpft_in

    allocate( EDPftvarcon_ptr%var1d(100)) ! Make this plenty large
	allocate( EDPftvarcon_ptr%var2d(100))
	iv1=0
	iv2=0

	! POINTER-SPECIFICATION-HERE (DO NOT REMOVE THIS LINE, OR MOVE IT)

	 allocate(EDPftvarcon_inst%allom_d2h1(1:numpft))
	 EDPftvarcon_inst%allom_d2h1(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2h1"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2h1
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2h2(1:numpft))
	 EDPftvarcon_inst%allom_d2h2(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2h2"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2h2
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2h3(1:numpft))
	 EDPftvarcon_inst%allom_d2h3(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2h3"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2h3
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_hmode(1:numpft))
	 EDPftvarcon_inst%allom_hmode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_hmode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_hmode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_dbh_maxheight(1:numpft))
	 EDPftvarcon_inst%allom_dbh_maxheight(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_dbh_maxheight"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_dbh_maxheight
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_agb1(1:numpft))
	 EDPftvarcon_inst%allom_agb1(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_agb1"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_agb1
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_agb2(1:numpft))
	 EDPftvarcon_inst%allom_agb2(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_agb2"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_agb2
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_agb3(1:numpft))
	 EDPftvarcon_inst%allom_agb3(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_agb3"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_agb3
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_agb4(1:numpft))
	 EDPftvarcon_inst%allom_agb4(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_agb4"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_agb4
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%wood_density(1:numpft))
	 EDPftvarcon_inst%wood_density(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_wood_density"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%wood_density
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%c2b(1:numpft))
	 EDPftvarcon_inst%c2b(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_c2b"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%c2b
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_agb_frac(1:numpft))
	 EDPftvarcon_inst%allom_agb_frac(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_agb_frac"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_agb_frac
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_amode(1:numpft))
	 EDPftvarcon_inst%allom_amode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_amode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_amode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_lmode(1:numpft))
	 EDPftvarcon_inst%allom_lmode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_lmode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_lmode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2bl1(1:numpft))
	 EDPftvarcon_inst%allom_d2bl1(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2bl1"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2bl1
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2bl2(1:numpft))
	 EDPftvarcon_inst%allom_d2bl2(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2bl2"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2bl2
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2bl3(1:numpft))
	 EDPftvarcon_inst%allom_d2bl3(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2bl3"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2bl3
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_blca_expnt_diff(1:numpft))
	 EDPftvarcon_inst%allom_blca_expnt_diff(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_blca_expnt_diff"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_blca_expnt_diff
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2ca_coefficient_min(1:numpft))
	 EDPftvarcon_inst%allom_d2ca_coefficient_min(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2ca_coefficient_min"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2ca_coefficient_min
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_d2ca_coefficient_max(1:numpft))
	 EDPftvarcon_inst%allom_d2ca_coefficient_max(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_d2ca_coefficient_max"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_d2ca_coefficient_max
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%slatop(1:numpft))
	 EDPftvarcon_inst%slatop(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_leaf_slatop"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%slatop
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%slamax(1:numpft))
	 EDPftvarcon_inst%slamax(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_leaf_slamax"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%slamax
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_sai_scaler(1:numpft))
	 EDPftvarcon_inst%allom_sai_scaler(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_sai_scaler"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_sai_scaler
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_smode(1:numpft))
	 EDPftvarcon_inst%allom_smode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_smode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_smode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_cmode(1:numpft))
	 EDPftvarcon_inst%allom_cmode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_cmode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_cmode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_fmode(1:numpft))
	 EDPftvarcon_inst%allom_fmode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_fmode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_fmode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_stmode(1:numpft))
	 EDPftvarcon_inst%allom_stmode(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_stmode"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_stmode
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%cushion(1:numpft))
	 EDPftvarcon_inst%cushion(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_alloc_storage_cushion"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%cushion
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_l2fr(1:numpft))
	 EDPftvarcon_inst%allom_l2fr(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_l2fr"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_l2fr
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_la_per_sa_int(1:numpft))
	 EDPftvarcon_inst%allom_la_per_sa_int(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_la_per_sa_int"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_la_per_sa_int
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%allom_la_per_sa_slp(1:numpft))
	 EDPftvarcon_inst%allom_la_per_sa_slp(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_allom_la_per_sa_slp"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%allom_la_per_sa_slp
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%roota_par(1:numpft))
	 EDPftvarcon_inst%roota_par(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_roota_par"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%roota_par
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%rootb_par(1:numpft))
	 EDPftvarcon_inst%rootb_par(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_rootb_par"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%rootb_par
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

	 allocate(EDPftvarcon_inst%rootprof_beta(1:numpft,1))
	 EDPftvarcon_inst%rootprof_beta(:,:) = nan
	 iv2 = iv2 + 1
	 EDPftvarcon_ptr%var2d(iv2)%var_name = "fates_rootprof_beta"
	 EDPftvarcon_ptr%var2d(iv2)%var_rp   => EDPftvarcon_inst%rootprof_beta
	 EDPftvarcon_ptr%var2d(iv2)%vtype    = 1

	 allocate(EDPftvarcon_inst%woody(1:numpft))
	 EDPftvarcon_inst%woody(:) = nan
	 iv1 = iv1 + 1
	 EDPftvarcon_ptr%var1d(iv1)%var_name = "fates_woody"
	 EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%woody
	 EDPftvarcon_ptr%var1d(iv1)%vtype    = 1

!    allocate( EDPftvarcon_inst%allom_dbh_maxheight   (1:numpft)); EDPftvarcon_inst%allom_dbh_maxheight (:) = nan
!    iv = iv + 1
!    EDPftvarcon_ptr%var1d(iv)%var_name = "fates_allom_dbh_maxheight"
!    EDPftvarcon_ptr%var1d(iv)%var_rp   => EDPftvarcon_inst%allom_dbh_maxheight
!    EDPftvarcon_ptr%var1d(iv)%vtype    = 1

    
    numparm1d = iv1
	numparm2d = iv2


    print*,"F90: ALLOCATED ",numparm1d," PARAMETERS, FOR ",numpft," PFTs"


    return
 end subroutine EDPftvarconAlloc

end module EDPftvarcon
