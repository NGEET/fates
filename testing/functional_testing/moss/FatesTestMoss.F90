program FatesTestMoss

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesMossMod,      only : available_light_under_canopy_and_moss, light_growth_multiplier
  use FatesMossMod,      only : moss_biomass_change_kg_per_m2
  use FatesMossMod,      only : Q_KG_PER_KGMOSS, B_KG_PER_KGMOSS
  use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=), isnan => shr_infnan_isnan

  implicit none

  ! LOCALS:
  integer                            :: i, j              ! looping indices
  integer                            :: n_moss_biomass    ! number of moss biomass levels to test
  integer                            :: n_cla   ! number of cumulative leaf area levels to test
  integer                            :: n_assim   ! number of assimilation levels to test
  real(r8), allocatable              :: cla(:)  ! cumulative leaf area in plot (m2)
  real(r8), allocatable              :: assim(:)  ! moss assimilation (kg/m2 plot)
  real(r8), allocatable              :: moss_biomass(:)  ! moss biomass (kg per m2 plot)
  real(r8), allocatable              :: out_al(:,:) ! output: available light
  real(r8), allocatable              :: out_algf(:,:) ! output: light growth multiplier
  real(r8), allocatable              :: out_resp(:,:) ! output: moss respiration (kg/m2 plot)
  real(r8), allocatable              :: out_mort(:,:) ! output: moss mortality (kg/m2 plot)
  real(r8), allocatable              :: out_biomass_delta(:,:) ! output: change in moss biomass (kg/m2 plot)
  real(r8) :: moss_biomass_after


  ! CONSTANTS:
  logical,          parameter :: debug = .false.
  character(len=*), parameter :: out_file = 'moss_out.nc'    ! output file
  real(r8),         parameter :: min_cla = 0._r8      ! minimum cumulative leaf area to calculate
  real(r8),         parameter :: max_cla = 8000.0_r8     ! maximum cumulative leaf area to calculate
  real(r8),         parameter :: cla_inc = 0.1_r8    ! cumulative leaf area increment to use
  real(r8),         parameter :: min_assim = 0._r8      ! minimum assimilation to calculate
  real(r8),         parameter :: max_assim = 5._r8     ! maximum assimilation to calculate
  real(r8),         parameter :: assim_inc = 0.001_r8    ! assimilation increment to use

  interface

    subroutine WriteMossData(out_file, n_cla, n_assim, n_moss_biomass, cla, assim, moss_biomass, out_al, out_algf, out_resp, out_mort, out_biomass_delta)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n_cla
      integer,          intent(in) :: n_assim
      integer,          intent(in) :: n_moss_biomass
      real(r8),         intent(in) :: cla(:)
      real(r8),         intent(in) :: assim(:)
      real(r8),         intent(in) :: moss_biomass(:)
      real(r8),         intent(in) :: out_al(:,:)
      real(r8),         intent(in) :: out_algf(:,:)
      real(r8),         intent(in) :: out_resp(:,:)
      real(r8),         intent(in) :: out_mort(:,:)
      real(r8),         intent(in) :: out_biomass_delta(:,:)
    end subroutine WriteMossData

  end interface

  ! determine sizes of arrays
  n_cla = int((max_cla - min_cla)/cla_inc + 1)
  n_assim = int((max_assim - min_assim)/assim_inc + 1)
  n_moss_biomass = 7

  ! allocate arrays
  allocate(cla(n_cla))
  allocate(assim(n_assim))
  allocate(moss_biomass(n_moss_biomass))
  allocate(out_al(n_cla, n_moss_biomass))
  allocate(out_algf(n_cla, n_moss_biomass))
  allocate(out_resp(n_assim, n_moss_biomass))
  allocate(out_mort(n_assim, n_moss_biomass))
  allocate(out_biomass_delta(n_assim, n_moss_biomass))

  ! initialize cla array
  do i = 1, n_cla
    cla(i) = min_cla + cla_inc*(i-1)
  end do

  ! initialize assim array
  do i = 1, n_assim
    assim(i) = min_assim + assim_inc*(i-1)
  end do

  ! initialize moss biomass array
  moss_biomass = (/ 0.1_r8, 1._r8, 3._r8, 6._r8, 9._r8, 12._r8, 15._r8 /)

  ! calculate light-related functions
  do i = 1, n_cla
    do j = 1, n_moss_biomass
      out_al(i,j) = available_light_under_canopy_and_moss(cla(i), moss_biomass(j))
      out_algf(i,j) = light_growth_multiplier(cla(i), moss_biomass(j))
    end do
  end do

  ! calculate respiration and mortality
  do i = 1, n_assim
    do j = 1, n_moss_biomass
      call moss_biomass_change_kg_per_m2(Q_KG_PER_KGMOSS, B_KG_PER_KGMOSS, assim(i), moss_biomass(j), out_resp(i,j), out_mort(i,j), moss_biomass_after)
      out_biomass_delta(i,j) = moss_biomass_after - moss_biomass(j)
    end do
  end do

  ! write out data to netcdf file
  call WriteMossData(out_file, n_cla, n_assim, n_moss_biomass, cla, assim, moss_biomass, out_al, out_algf, out_resp, out_mort, out_biomass_delta)

end program FatesTestMoss

! ----------------------------------------------------------------------------------------

subroutine WriteMossData(out_file, n_cla, n_assim, n_moss_biomass, cla, assim, moss_biomass, out_al, out_algf, out_resp, out_mort, out_biomass_delta)
  !
  ! DESCRIPTION:
  ! Writes out data from the moss test
  !
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int
  use FatesConstantsMod,  only : r8 => fates_r8
  implicit none

  ! ARGUMENTS
  character(len=*), intent(in) :: out_file
  integer,          intent(in) :: n_cla
  integer,          intent(in) :: n_assim
  integer,          intent(in) :: n_moss_biomass
  real(r8),         intent(in) :: cla(:)
  real(r8),         intent(in) :: assim(:)
  real(r8),         intent(in) :: moss_biomass(:)
  real(r8),         intent(in) :: out_al(:,:)
  real(r8),         intent(in) :: out_algf(:,:)
  real(r8),         intent(in) :: out_resp(:,:)
  real(r8),         intent(in) :: out_mort(:,:)
  real(r8),         intent(in) :: out_biomass_delta(:,:)

  ! LOCALS:
  integer, parameter   :: n_dims = 3     ! number of dimensions
  integer, parameter   :: n_dim_atts = 2     ! number of attributes for dimensions
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=24)    :: dim_names(n_dims)   ! dimension name(s)
  integer              :: dimIDs(n_dims)      ! dimension ID(s)
  integer              :: claID, assimID, mossbiomassID ! variable ID(s) for dimensions
  integer              :: alID, algfID, respID, mortID, biomassdeltaID

  ! dimension name(s)
  dim_names = [character(len=24) :: 'cla', 'moss_biomass', 'assim']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n_cla, n_moss_biomass, n_assim/), n_dims, dimIDs)

  ! register cumulative leaf area dimension
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'm2', 'Cumulative leaf area across entire plot'], n_dim_atts, claID)

  ! register moss biomass dimension
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'kg/m2 plot', 'Moss biomass (before)'], n_dim_atts, mossbiomassID)

  ! register assimilation dimension
  call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'kg/m2', 'Moss assimilation'], n_dim_atts, assimID)

  ! register out_al
  call RegisterVar(ncid, 'out_al', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'cla moss_biomass', 'units?', 'Available light under canopy and moss'],  &
    3, alID)

  ! register out_algf
  call RegisterVar(ncid, 'out_algf', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'cla moss_biomass', 'unitless?', 'Light growth multiplier'],  &
    3, algfID)

  ! register out_resp
  call RegisterVar(ncid, 'out_resp', (/ dimIDs(3), dimIDs(2) /), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'assim moss_biomass', 'kg/m2', 'Moss respiration'],  &
    3, respID)

  ! register out_mort
  call RegisterVar(ncid, 'out_mort', (/ dimIDs(3), dimIDs(2) /), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'assim moss_biomass', 'kg/m2', 'Moss mortality'],  &
    3, mortID)

  ! register out_biomass_delta
  call RegisterVar(ncid, 'out_biomass_delta', (/ dimIDs(3), dimIDs(2) /), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'assim moss_biomass', 'kg/m2', 'Change in moss biomass'],  &
    3, biomassdeltaID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, claID, cla(:))
  call WriteVar(ncid, mossbiomassID, moss_biomass(:))
  call WriteVar(ncid, assimID, assim(:))
  call WriteVar(ncid, alID, out_al(:,:))
  call WriteVar(ncid, algfID, out_algf(:,:))
  call WriteVar(ncid, respID, out_resp(:,:))
  call WriteVar(ncid, mortID, out_mort(:,:))
  call WriteVar(ncid, biomassdeltaID, out_biomass_delta(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteMossData
