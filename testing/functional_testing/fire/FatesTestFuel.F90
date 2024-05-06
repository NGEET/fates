program FatesTestFuel
  
  use FatesConstantsMod,           only : r8 => fates_r8 
  use EDTypesMod,                  only : ed_site_type
  use FatesTestFireMod,            only : SetUpFuel, ReadDatmData, WriteFireData
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use SyntheticFuelTypes,          only : fuel_types_array_class
  use SFFireWeatherMod,            only : fire_weather
  use SFNesterovMod,               only : nesterov_index
  use FatesFuelMod,                only : fuel_type
  use FatesFuelClassesMod,         only : nfsc, nfsc_notrunks
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader      ! param reader instance
  type(fuel_types_array_class)                   :: fuel_types_array  ! array of fuel models
  class(fire_weather),               pointer     :: fireWeather       ! fire weather object
  type(fuel_type),                   allocatable :: fuel(:)           ! fuel objects                          
  character(len=:),                  allocatable :: param_file        ! input parameter file
  character(len=:),                  allocatable :: datm_file         ! input DATM driver file
  real(r8),                          allocatable :: temp_degC(:)      ! daily air temperature [degC]
  real(r8),                          allocatable :: precip(:)         ! daily precipitation [mm]
  real(r8),                          allocatable :: rh(:)             ! daily relative humidity [%]
  real(r8),                          allocatable :: wind(:)           ! daily wind speed [m/s]
  real(r8),                          allocatable :: NI(:)             ! Nesterov index
  real(r8),                          allocatable :: fuel_loading(:,:) ! fuel loading [kgC/m2]
  real(r8),                          allocatable :: total_loading(:)  ! total fuel loading [kgC/m2]
  real(r8),                          allocatable :: frac_loading(:,:) ! fractional fuel loading [0-1]
  integer                                        :: i                 ! looping index
  integer                                        :: num_fuel_models   ! number of fuel models to test
  
  ! CONSTANTS:
  integer,          parameter :: n_days = 365             ! number of days to run simulation
  character(len=*), parameter :: out_file = 'fuel_out.nc' ! output file 
  
  ! fuel models to test
  integer, parameter, dimension(3) :: fuel_models = (/102, 183, 164/)
  
  ! number of fuel models to test
  num_fuel_models = size(fuel_models)
  
  ! allocate arrays
  allocate(temp_degC(n_days))
  allocate(precip(n_days))
  allocate(rh(n_days))
  allocate(wind(n_days))
  allocate(NI(n_days))
  allocate(fuel_loading(nfsc_notrunks, num_fuel_models))
  allocate(frac_loading(nfsc_notrunks, num_fuel_models))
  allocate(total_loading(num_fuel_models))
  
  ! read in parameter file name and DATM file from command line
  param_file = command_line_arg(1)
  datm_file = command_line_arg(2)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! read in DATM data
  call ReadDatmData(datm_file, temp_degC, precip, rh, wind)
  
  ! set up fire weather class
  allocate(nesterov_index :: fireWeather)
  call fireWeather%Init()
  
  ! set up fuel objects and calculate loading
  allocate(fuel(num_fuel_models))
  call fuel_types_array%GetFuelModels()
  do i = 1, num_fuel_models
    
    ! uses data from fuel_models to initialize fuel
    call SetUpFuel(fuel(i), fuel_types_array, fuel_models(i))
    
    ! sum up fuel and calculate loading
    call fuel(i)%SumLoading()
    call fuel(i)%CalculateFractionalLoading()
    fuel_loading(:,i) = fuel(i)%loading(:)
    total_loading(i) = fuel(i)%total_loading
    frac_loading(:,i) = fuel(i)%frac_loading(:)
  end do
  
  ! run on time steps
  do i = 1, n_days
    call fireWeather%UpdateIndex(temp_degC(i), precip(i), rh(i), wind(i))
    NI(i) = fireWeather%fire_weather_index
  end do 
  
  ! write out data
  call WriteFireData(out_file, n_days, num_fuel_models, temp_degC, precip, rh, NI,       &
    fuel_loading, frac_loading, total_loading, fuel_models)
  
end program FatesTestFuel