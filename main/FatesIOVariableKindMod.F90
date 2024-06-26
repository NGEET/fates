module FatesIOVariableKindMod

  use FatesConstantsMod, only : fates_long_string_length
  use FatesGlobals, only : fates_log
  use FatesIODimensionsMod, only : fates_io_dimension_type
  use FatesGlobals          , only : endrun => fates_endrun
  use shr_log_mod           , only : errMsg => shr_log_errMsg

  
  implicit none
  private

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
  ! FIXME(bja, 2016-10) do these need to be strings, or can they be integer enumerations?
  ! FIXME(rgk, 2016-11) these should probably be moved to varkindmod?
  
  character(*), parameter, public :: site_r8 = 'SI_R8'
  character(*), parameter, public :: site_int = 'SI_INT'
  character(*), parameter, public :: site_soil_r8 = 'SI_SOIL_R8'
  character(*), parameter, public :: site_size_pft_r8 = 'SI_SCPF_R8'
  character(*), parameter, public :: site_size_r8 = 'SI_SCLS_R8'
  character(*), parameter, public :: site_coage_pft_r8 = 'SI_CAPF_R8'
  character(*), parameter, public :: site_coage_r8 = 'SI_CACLS_R8'
  character(*), parameter, public :: cohort_r8 = 'CO_R8'
  character(*), parameter, public :: cohort_int = 'CO_INT'
  character(*), parameter, public :: site_pft_r8 = 'SI_PFT_R8'
  character(*), parameter, public :: site_age_r8 = 'SI_AGE_R8'
  character(*), parameter, public :: site_height_r8 = 'SI_HEIGHT_R8'
  character(*), parameter, public :: site_fuel_r8 = 'SI_FUEL_R8'
  character(*), parameter, public :: site_cwdsc_r8 = 'SI_CWDSC_R8'
  character(*), parameter, public :: site_can_r8 = 'SI_CAN_R8'
  character(*), parameter, public :: site_cnlf_r8 = 'SI_CNLF_R8'
  character(*), parameter, public :: site_cdpf_r8 = 'SI_CDPF_R8'
  character(*), parameter, public :: site_cdsc_r8 = 'SI_CDSC_R8'
  character(*), parameter, public :: site_cdam_r8 = 'SI_CDAM_R8'
  character(*), parameter, public :: site_cnlfpft_r8 = 'SI_CNLFPFT_R8'
  character(*), parameter, public :: site_scag_r8 = 'SI_SCAG_R8'
  character(*), parameter, public :: site_scagpft_r8 = 'SI_SCAGPFT_R8'
  character(*), parameter, public :: site_agepft_r8 = 'SI_AGEPFT_R8'
  character(*), parameter, public :: site_agefuel_r8 = 'SI_AGEFUEL_R8'
  character(*), parameter, public :: site_clscpf_r8 = 'SI_CLSCPF_R8'
  character(*), parameter, public :: site_landuse_r8 = 'SI_LANDUSE_R8'
  character(*), parameter, public :: site_lulu_r8 = 'SI_LULU_R8'
  character(*), parameter, public :: site_lupft_r8 = 'SI_LUPFT_R8'
  
  ! Element, and multiplexed element dimensions
  character(*), parameter, public :: site_elem_r8  = 'SI_ELEM_R8'
  character(*), parameter, public :: site_elpft_r8 = 'SI_ELEMPFT_R8'
  character(*), parameter, public :: site_elcwd_r8 = 'SI_ELEMCWD_R8'
  character(*), parameter, public :: site_elage_r8 = 'SI_ELEMAGE_R8'

  ! ------------------------------------------------------------------
  !
  ! History Variable Groups
  !
  ! These are group indices for output variables. We use
  ! these groups to do things like zero-ing and initializing
  !
  ! These groups are updated at the dynamics (daily) step
  ! so they are turned on and off with dimlevel(2)
  !
  ! active when dimlevel(2)>0
  integer, parameter, public :: group_dyna_simple = 1
  integer, parameter, public :: group_nflx_simple = 7
  
  ! active when dimlevel(2)>1
  integer, parameter, public :: group_dyna_complx = 2
  integer, parameter, public :: group_nflx_complx = 8

  ! These groups are updated at the fast step
  ! so they are turned on and off with dimlevel(1)
  !
  ! active when dimlevel(1)>0
  integer, parameter, public :: group_hifr_simple = 3
  integer, parameter, public :: group_hydr_simple = 5

  ! active when dimlevel(1)>1
  integer, parameter, public :: group_hifr_complx = 4
  integer, parameter, public :: group_hydr_complx = 6

  ! -------------------------------------------------------------------
  
  ! NOTE(RGK, 2016) %active is not used yet. Was intended as a check on the HLM->FATES
  ! control parameter passing to ensure all active dimension types received all
  ! dimensioning specifications from the host, but we currently arent using those
  ! passing functions..

  ! This structure is not multi-threaded
  type, public :: fates_io_variable_kind_type
     character(len=fates_long_string_length) :: name ! String labelling this IO type
     integer              :: ndims       ! number of dimensions in this IO type
     integer, allocatable :: dimsize(:)  ! The size of each dimension
     logical, private :: active_
     integer :: dim1_index
     integer :: dim2_index

   contains

     procedure :: Init
     procedure :: set_active
     procedure :: is_active

  end type fates_io_variable_kind_type

  ! Make necessary functions public
  public :: iotype_index

contains

  ! ===================================================================================
  subroutine Init(this, name, num_dims)

    use FatesConstantsMod, only : fates_unset_int
    
    implicit none

    class(fates_io_variable_kind_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: num_dims
    
    this%name = trim(name)
    this%ndims = num_dims
    allocate(this%dimsize(this%ndims))
    this%dimsize(:) = fates_unset_int
    this%active_ = .false.
    this%dim1_index = fates_unset_int
    this%dim2_index = fates_unset_int
    
  end subroutine Init
  
 ! =======================================================================
 subroutine set_active(this)
   implicit none
   class(fates_io_variable_kind_type), intent(inout) :: this
   this%active_ = .true.
 end subroutine set_active

 logical function is_active(this)
   implicit none
   class(fates_io_variable_kind_type), intent(in) :: this
   is_active = this%active_
 end function is_active

  ! ====================================================================================

  function iotype_index(iotype_name, num_dim_kinds, dim_kinds) result(dk_index)

    ! argument
    character(len=*), intent(in) :: iotype_name
    integer, intent(in) :: num_dim_kinds
    type(fates_io_variable_kind_type), intent(in) :: dim_kinds(:)

    ! local
    integer :: dk_index

    do dk_index=1, num_dim_kinds
       if (trim(iotype_name) .eq. trim(dim_kinds(dk_index)%name)) then
          return
       end if
    end do
    write(fates_log(),*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    call endrun(msg=errMsg(sourcefile, __LINE__))

  end function iotype_index
  
end module FatesIOVariableKindMod
