program FatesUnitTestSF
  !
  ! DESCRIPTION:
  !		Test the FATES SPITFIRE model
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestSFMod, only : ReadDatmData, WriteFireData
  use FatesUnitTestIOMod, only : logf, OpenFile
  use SFFireWeatherMod,   only : fire_weather
  use SFNesterovMod,      only : nesterov_index
  use SFParamsMod,        only : min_precip_thresh

  implicit none

  ! LOCALS:
  class(fire_weather), allocatable :: fireWeather  ! fire weather object
  real(r8), allocatable            :: temp_degC(:) ! daily air temperature [degC]
  real(r8), allocatable            :: precip(:)    ! daily precipitation [mm]
  real(r8), allocatable            :: rh(:)        ! daily relative humidity [%]
  real(r8), allocatable            :: wind(:)      ! daily wind speed [m/s]
  integer,  allocatable            :: time(:)      ! time 
  real(r8), allocatable            :: NI(:)        ! cumulative nesterov index
  integer                          :: n = 365
  integer                          :: i

  ! open log file
  logf = OpenFile("log.txt", mode='rw')

  ! allocate arrays
  allocate(temp_degC(n))
  allocate(precip(n))
  allocate(rh(n))
  allocate(wind(n))
  allocate(time(n))
  allocate(NI(n))

  ! read in DATM data
  call ReadDatmData('BONA_datm.nc', temp_degC, precip, rh, wind)

  ! initialize fire weather
  allocate(nesterov_index :: fireWeather)
  fireWeather%fire_weather_index = 0.0_r8
  
  ! run on time steps
  do i = 1, n
    time(i) = i
    call fireWeather%Calculate(temp_degC(i), precip(i), rh(i))
    NI(i) = fireWeather%fire_weather_index
  end do 

  ! write out data
  call WriteFireData('Fire_unit_out.nc', n, time, temp_degC, precip, rh, NI)

end program FatesUnitTestSF