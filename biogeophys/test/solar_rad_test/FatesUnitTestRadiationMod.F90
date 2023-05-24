module FatesUnitTestRadiationMod
  !
  ! DESCRIPTION:
  !		Module to support testing the FATES radiation schemes
  !
  use FatesConstantsMod,        only : r8 => fates_r8
  use FatesUnitTestIOMod,       only : MAX_PATH, MAX_CHAR, OpenFile, logf
  use FatesUnitTestIOMod,       only : GetVar, OpenNCFile, CloseNCFile
  use EDPftvarcon,              only : EDPftvarcon_type
  use FatesParametersInterface, only : fates_parameters_type
  use FatesPatchMod,            only : fates_patch_type

  use shr_kind_mod,             only : SHR_KIND_R8
  use shr_orb_mod,              only : SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use shr_orb_mod,              only : shr_orb_cosz, shr_orb_decl, shr_orb_params

  implicit none
  private

  real, parameter :: PI = 4.0*atan(1.0)   ! Pi (3.1415..)
  real, parameter :: DEG2RAD = PI/180.0   ! Convert from degrees to radians

  type, public :: fates_rad_test
    
    real(r8)                :: lat         ! latitude to simulate [degrees]
    real(r8)                :: lon         ! longitude to simulate [degrees]
    integer                 :: year        ! year to simulate
    integer                 :: jday_start  ! Julian day to start simulation
    real(SHR_KIND_R8)       :: eccen       ! orbital eccentricity
    real(SHR_KIND_R8)       :: mvelp       ! moving vernal equinox long
    real(SHR_KIND_R8)       :: obliqr      ! Earth's obliquity [radians]
    real(SHR_KIND_R8)       :: lambm0      ! mean long of perihelion at vernal equinox [radians]
    real(SHR_KIND_R8)       :: mvelpp      ! moving vernal equinox long of perihelion plus pi [radians]
    integer                 :: num_pft     ! number of pfts to simulate
    integer                 :: num_swb     ! number of shortwave bands to simulate
    real(r8)                :: fcansno     ! fraction of canopy covered by snow [0-1]
    character(len=MAX_PATH) :: patch_file  ! patch data file name
    character(len=MAX_PATH) :: out_file    ! output file name

    contains 

      procedure :: Init
      procedure :: SetRadiationParams
      procedure :: InitPftData
      procedure :: ReadPftNamelist
      procedure :: ReadPatchData
      procedure :: InitPatch
      procedure :: GetOrbitalParams
      procedure :: CalcCosZ

  end type fates_rad_test

  contains
  
  !=======================================================================================

  subroutine Init(this, EDPftvarconInst, rad_nl_file, pft_nl_file)
    !
    !  DESCRIPTION:
    !  Initialize fates_rad_test type
    !

    ! ARGUMENTS:
    class(fates_rad_test),   intent(inout) :: this            ! rad object
    class(EDPftvarcon_type), intent(inout) :: EDPftvarconInst ! PFT variable instance
    character(len=MAX_PATH), intent(in)    :: rad_nl_file     ! radiation namelist file name
    character(len=MAX_PATH), intent(in)    :: pft_nl_file     ! pft namelist name

    ! read in the radiation namelist and set values
    call this%SetRadiationParams(rad_nl_file)

    ! read in pft parameters
    call this%InitPftData(EDPftvarconInst, pft_nl_file)

    ! set these to undefined so we can catch errors later
    this%eccen = SHR_ORB_UNDEF_REAL
    this%mvelp = SHR_ORB_UNDEF_REAL

  end subroutine Init

  !=======================================================================================

  subroutine SetRadiationParams(this, rad_nl_file)
    !
    ! DESCRIPTION:
    !		read in the namelist associated with the radiation unit tests and 
    !   initialize values
    !

    ! ARGUMENTS:
    class(fates_rad_test),   intent(inout) :: this        ! rad object
    character(len=MAX_PATH), intent(in)    :: rad_nl_file ! radiation namelist file name 

    
    ! LOCALS:
    character(len=MAX_PATH) :: patch_file ! patch data file name
    character(len=MAX_PATH) :: out_file   ! output file name
    integer                 :: year, jday ! year and day of year to start simulation
    integer                 :: num_pft    ! number of pfts
    integer                 :: num_swb    ! number of shortwave bands
    real(r8)                :: lat, lon   ! latitude and longitude [radians]
    real(r8)                :: fcansno    ! fraction of canopy covered by snow [0-1]
    character(len=MAX_CHAR) :: msg        ! I/O Error message
    integer                 :: rad_nl     ! unit number for namelist
    integer                 :: ios        ! I/O status

    ! Namelist of radiation parameters
    namelist /radiation/ year, jday, lat, lon, num_pft, num_swb, fcansno, patch_file,    &
      out_file
  
    ! Now read parameters namelist
    rad_nl = OpenFile(trim(rad_nl_file), 'r')
    read(rad_nl, radiation, iostat=ios, iomsg=msg)

    this%year = year
    this%jday_start = jday 
    this%lat = lat*DEG2RAD
    this%lon = lon*DEG2RAD
    this%num_pft = num_pft 
    this%fcansno = fcansno 
    this%patch_file = patch_file 
    this%out_file = out_file
    this%num_pft = num_pft
    this%num_swb = num_swb
    
    if (ios /= 0) then
      ! Problem reading file - tell user.
      write(logf, '(A, I6, A, A)') "Error reading radiation namelist file", ios,         &
        " IOMSG: ", msg
      stop 
    end if
  
    close(rad_nl)
  
  end subroutine SetRadiationParams

  !=======================================================================================

  subroutine InitPftData(this, EDPftvarconInst, pft_nl_file)
    !
    ! DESCRIPTION:
    !		read in the namelist associated with the pft-specific parameters and 
    !   initialize data
    !

    ! ARGUMENTS:
    class(fates_rad_test),   intent(inout) :: this            ! rad object
    character(len=MAX_PATH), intent(in)    :: pft_nl_file     ! pft namelist name
    class(EDPftvarcon_type), intent(inout) :: EDPftvarconInst ! PFT variable instance
  
    ! LOCALS:
    class(fates_parameters_type), allocatable :: fates_params                    ! FATES input parameters
    real(r8)                                  :: clumping_index(this%num_pft)    ! clumping index
    real(r8)                                  :: xl(this%num_pft)                ! leaf-stem orientation index
    real(r8)                                  :: rhol(this%num_pft,this%num_swb) ! leaf reflectance [0-1]
    real(r8)                                  :: rhos(this%num_pft,this%num_swb) ! stem reflectance [0-1]
    real(r8)                                  :: taul(this%num_pft,this%num_swb) ! leaf transmittance [0-1]
    real(r8)                                  :: taus(this%num_pft,this%num_swb) ! stem transmittance [0-1]
  
    ! initialize EDPFTvarcon_inst and fates_params
    call EDPftvarconInst%Init()
    allocate(fates_params)
    call fates_params%Init()
    call EDPftvarconInst%Register(fates_params)
  
    ! read in parameter values
    call this%ReadPftNamelist(pft_nl_file, clumping_index, xl, rhol, rhos, taul, taus)
  
    ! set values
    ! TODO: make this read in from a parameter file, this is hacky
    EDPftvarconInst%clumping_index = clumping_index
    EDPftvarconInst%xl = xl
    EDPftvarconInst%rhol = rhol
    EDPftvarconInst%rhos = rhos
    EDPftvarconInst%taul = taul
    EDPftvarconInst%taus = taus
  
  end subroutine InitPftData

  !=======================================================================================

  subroutine ReadPftNamelist(this, pft_nl_file, clumping_index, xl, rhol, rhos, taul,    &
      taus)
    !
    ! DESCRIPTION:
    !		read in the namelist associated with the pft-specific parameters
    !
    
    ! ARGUMENTS:
    class(fates_rad_test),   intent(inout) :: this                            ! rad object
    character(len=MAX_PATH), intent(in)    :: pft_nl_file                     ! pft namelist name
    real(r8),                intent(out)   :: clumping_index(this%num_pft)    ! clumping index
    real(r8),                intent(out)   :: xl(this%num_pft)                ! leaf-stem orientation index
    real(r8),                intent(out)   :: rhol(this%num_pft,this%num_swb) ! leaf reflectance [0-1]
    real(r8),                intent(out)   :: rhos(this%num_pft,this%num_swb) ! stem reflectance [0-1]
    real(r8),                intent(out)   :: taul(this%num_pft,this%num_swb) ! leaf transmittance [0-1]
    real(r8),                intent(out)   :: taus(this%num_pft,this%num_swb) ! stem transmittance [0-1]
    
    ! LOCALS:
    character(len=MAX_CHAR) :: msg                       ! I/O Error message
    integer                 :: pft_nl                    ! unit number for namelist
    integer                 :: ios                       ! I/O status
    real(r8)                :: leaf_rhonir(this%num_pft) ! leaf NIR reflectance [0-1]
    real(r8)                :: leaf_taunir(this%num_pft) ! leaf NIR transmittance [0-1]
    real(r8)                :: leaf_rhovis(this%num_pft) ! leaf visible reflectance [0-1]
    real(r8)                :: leaf_tauvis(this%num_pft) ! leaf visible transmittance [0-1]
    real(r8)                :: stem_rhonir(this%num_pft) ! stem NIR reflectance [0-1]
    real(r8)                :: stem_taunir(this%num_pft) ! stem NIR transmittance [0-1]
    real(r8)                :: stem_rhovis(this%num_pft) ! stem visible reflectance [0-1]
    real(r8)                :: stem_tauvis(this%num_pft) ! stem visible transmittance [0-1]
    
    ! Namelist of pft-specific parameters
    namelist /params/ clumping_index, xl, leaf_rhonir, leaf_taunir, leaf_rhovis,         &
      leaf_tauvis, stem_rhonir, stem_taunir, stem_rhovis, stem_tauvis
    
    ! Now read parameters namelist
    pft_nl = OpenFile(trim(pft_nl_file), 'r')
    read(pft_nl, params, iostat=ios, iomsg=msg)
    
    if (ios /= 0) then
      ! Problem reading file - tell user.
      write(logf, '(A, I6, A, A)') "Error reading pft namelist file", ios, " IOMSG: ", msg
      stop 
    end if
    
    ! combine arrays
    rhol(:,1) = leaf_rhovis
    rhol(:,2) = leaf_rhonir
    taul(:,1) = leaf_tauvis
    taul(:,2) = leaf_taunir
    rhos(:,1) = stem_rhovis
    rhos(:,2) = stem_rhonir
    taus(:,1) = stem_tauvis
    taus(:,2) = stem_taunir
    
    close(pft_nl)
    
  end subroutine ReadPftNamelist

  !=======================================================================================

  subroutine ReadPatchData(this, canopy_area, elai, esai, nrad)
    !
    ! DESCRIPTION:
    ! Reads and return patch data
    !
  
    ! ARGUMENTS:
    class(fates_rad_test), intent(inout) :: this               ! rad object
    real(r8), allocatable, intent(out)   :: canopy_area(:,:,:) ! canopy area profile
    real(r8), allocatable, intent(out)   :: elai(:,:,:)        ! exposed lai profile
    real(r8), allocatable, intent(out)   :: esai(:,:,:)        ! exposed sai profile
    integer,  allocatable, intent(out)   :: nrad(:,:)          ! number of exposed leaf layers

    ! LOCALS:
    integer :: ncid ! netcdf file unit number

    ! open file
    call OpenNCFile(trim(this%patch_file), ncid)
    
    ! read in data
    call GetVar(ncid, 'can_area', canopy_area)
    call GetVar(ncid, 'elai', elai)
    call GetVar(ncid, 'esai', esai)
    call GetVar(ncid, 'nrad', nrad)

    ! close file
    call CloseNCFile(ncid)

  end subroutine ReadPatchData

  !=======================================================================================

  subroutine InitPatch(this, fatesPatch)
    !
    ! DESCRIPTION:
    ! Reads and return patch data
    !
  
    ! ARGUMENTS:
    class(fates_rad_test),   intent(inout) :: this       ! rad object
    class(fates_patch_type), intent(inout) :: fatesPatch ! patch object

    ! LOCALS:
    real(r8), allocatable :: canopy_area(:,:,:) ! canopy area profile
    real(r8), allocatable :: elai(:,:,:)        ! exposed lai profile
    real(r8), allocatable :: esai(:,:,:)        ! exposed sai profile
    integer,  allocatable :: nrad(:,:)          ! number of exposed leaf layers

    ! PARAMETERS
    integer, parameter :: nlev_soil   = 5 ! this shouldn't really matter

    call this%ReadPatchData(canopy_area, elai, esai, nrad)

    ! initialize the FATES patch and set values
    call fatesPatch%init(this%num_swb, nlev_soil)
    fatesPatch%canopy_area_profile = canopy_area
    fatesPatch%elai_profile = elai 
    fatesPatch%esai_profile = esai
    fatesPatch%nrad = nrad

    ! set the fcansno
    fatesPatch%fcansno = this%fcansno

  end subroutine InitPatch

  !=======================================================================================

  subroutine GetOrbitalParams(this)
    !
    ! DESCRIPTION:
    ! Calculates orbital values
    !
  
    ! ARGUMENTS:
    class(fates_rad_test), intent(inout) :: this ! rad object

    ! LOCALS:
    real(SHR_KIND_R8) :: obliq ! obliquity in degrees

    obliq = SHR_ORB_UNDEF_REAL
    call shr_orb_params(this%year, this%eccen, obliq, this%mvelp, this%obliqr,           &
      this%lambm0, this%mvelpp, .false.)

  end subroutine GetOrbitalParams

  !=======================================================================================

  subroutine CalcCosZ(this, time_step, cosz)
    !
    ! DESCRIPTION:
    ! Calculates cosine of the solar zenith angle.
    ! Assumes 365.0 days/year.
    !
  
    ! ARGUMENTS:
    class(fates_rad_test), intent(inout) :: this      ! rad object
    integer,               intent(in)    :: time_step ! time step
    real(SHR_KIND_R8),     intent(out)   :: cosz      ! cosine of solar zenith angle [radians]

    ! LOCALS:
    real(SHR_KIND_R8) :: cal_day ! Julian cal day (1.xx to 365.xx)
    real(SHR_KIND_R8) :: declin  ! solar declination [radians]
    real(SHR_KIND_R8) :: eccf    ! Earth-sun distance factor (ie. (1/r)**2)
    
    ! calculate proportional day
    cal_day = float(this%jday_start) + float(time_step)*0.02083333_SHR_KIND_R8
    
    ! calculate solar declination
    call shr_orb_decl(cal_day, this%eccen, this%mvelpp, this%lambm0, this%obliqr,        &
      declin, eccf)

    ! calculate cosine of solar zenith angle
    cosz = shr_orb_cosz(cal_day, this%lat, this%lon, declin)

  end subroutine CalcCosZ

end module FatesUnitTestRadiationMod