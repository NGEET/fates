program FatesTestMoss

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesMossMod,      only : available_light_under_canopy_and_moss, light_growth_multiplier
  use shr_infnan_mod,    only : nan => shr_infnan_nan, assignment(=), isnan => shr_infnan_isnan

  implicit none

  ! LOCALS:
  integer                            :: i, j              ! looping indices
  integer                            :: n_moss_biomass    ! number of moss biomass levels to test
  integer                            :: n_cla   ! number of cumulative leaf area levels to test
  real(r8), allocatable              :: cla(:)  ! cumulative leaf area in plot (m2)
  real(r8), allocatable              :: moss_biomass(:)  ! moss biomass (kg per m2 plot)
  real(r8), allocatable              :: out_al(:,:) ! output: available light
  real(r8), allocatable              :: out_algf(:,:) ! output: light growth multiplier


  ! CONSTANTS:
  logical,          parameter :: debug = .false.
  character(len=*), parameter :: out_file = 'moss_out.nc'    ! output file
  real(r8),         parameter :: min_cla = 0._r8      ! minimum cumulative leaf area to calculate
  real(r8),         parameter :: max_cla = 8000.0_r8     ! maximum cumulative leaf area to calculate
  real(r8),         parameter :: cla_inc = 0.1_r8    ! cumulative leaf area increment to use

  interface

    subroutine WriteMossData(out_file, n_cla, n_moss_biomass, cla, moss_biomass, out_al, out_algf)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n_cla
      integer,          intent(in) :: n_moss_biomass
      real(r8),         intent(in) :: cla(:)
      real(r8),         intent(in) :: moss_biomass(:)
      real(r8),         intent(in) :: out_al(:,:)
      real(r8),         intent(in) :: out_algf(:,:)
    end subroutine WriteMossData

  end interface

  ! determine sizes of arrays
  n_cla = int((max_cla - min_cla)/cla_inc + 1)
  n_moss_biomass = 9

  ! allocate arrays
  allocate(cla(n_cla))
  allocate(moss_biomass(n_moss_biomass))
  allocate(out_al(n_cla, n_moss_biomass))
  allocate(out_algf(n_cla, n_moss_biomass))

  ! initialize cla array
  do i = 1, n_cla
    cla(i) = min_cla + cla_inc*(i-1)
  end do

  ! initialize moss biomass array: Powers of 10 symmetrical around 10^0
  do i = 1, n_moss_biomass
    moss_biomass(i) = 10._r8**(floor(-real(n_moss_biomass, r8) / 2) + i)
  end do

  ! calculate light-related functions
  do i = 1, n_cla
    do j = 1, n_moss_biomass
      out_al(i,j) = available_light_under_canopy_and_moss(cla(i), moss_biomass(j))
      out_algf(i,j) = light_growth_multiplier(cla(i), moss_biomass(j))
    end do
  end do

  ! write out data to netcdf file
  call WriteMossData(out_file, n_cla, n_moss_biomass, cla, moss_biomass, out_al, out_algf)

end program FatesTestMoss

! ----------------------------------------------------------------------------------------

subroutine WriteMossData(out_file, n_cla, n_moss_biomass, cla, moss_biomass, out_al, out_algf)
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
  integer,          intent(in) :: n_moss_biomass
  real(r8),         intent(in) :: cla(:)
  real(r8),         intent(in) :: moss_biomass(:)
  real(r8),         intent(in) :: out_al(:,:)
  real(r8),         intent(in) :: out_algf(:,:)

  ! LOCALS:
  ! TODO: Set a local parameter for number of dimensions (2)
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=24)    :: dim_names(2)   ! dimension name(s)
  integer              :: dimIDs(2)      ! dimension ID(s)
  integer              :: claID, mossbiomassID ! variable ID(s) for dimensions
  integer              :: alID, algfID

  ! dimension name(s)
  dim_names = [character(len=24) :: 'cla', 'moss_biomass']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n_cla, n_moss_biomass/), 2, dimIDs)

  ! register cumulative leaf area dimension
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'm2', 'Cumulative leaf area across entire plot'], 2, claID)

  ! register moss biomass dimension
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'kg/m2 plot', 'Moss biomass'], 2, mossbiomassID)

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

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, claID, cla(:))
  call WriteVar(ncid, mossbiomassID, moss_biomass(:))
  call WriteVar(ncid, alID, out_al(:,:))
  call WriteVar(ncid, algfID, out_algf(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteMossData
