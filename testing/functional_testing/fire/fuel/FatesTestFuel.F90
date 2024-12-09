program FatesTestFuel
  
  use FatesConstantsMod,           only : r8 => fates_r8 
  use EDTypesMod,                  only : ed_site_type
  use FatesTestFireMod,            only : SetUpFuel, ReadDatmData, WriteFireData
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use SyntheticFuelModels,         only : fuel_models_array_class
  use SFFireWeatherMod,            only : fire_weather
  use SFNesterovMod,               only : nesterov_index
  use FatesFuelMod,                only : fuel_type
  use FatesFuelClassesMod,         only : num_fuel_classes
  use SFParamsMod,                 only : SF_val_SAV, SF_val_drying_ratio
  use SFParamsMod,                 only : SF_val_FBD
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader         ! param reader instance
  type(fuel_models_array_class)                  :: fuel_models_array    ! array of fuel models
  class(fire_weather),               pointer     :: fireWeather          ! fire weather object
  type(fuel_type),                   allocatable :: fuel(:)              ! fuel objects                          
  character(len=:),                  allocatable :: param_file           ! input parameter file
  character(len=:),                  allocatable :: datm_file            ! input DATM driver file
  real(r8),                          allocatable :: temp_degC(:)         ! daily air temperature [degC]
  real(r8),                          allocatable :: precip(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: rh(:)                ! daily relative humidity [%]
  real(r8),                          allocatable :: wind(:)              ! daily wind speed [m/s]
  real(r8),                          allocatable :: NI(:)                ! Nesterov index
  real(r8),                          allocatable :: fuel_loading(:,:)    ! fuel loading [kgC/m2]
  real(r8),                          allocatable :: non_trunk_loading(:) ! non-trunk fuel loading [kgC/m2]
  real(r8),                          allocatable :: frac_loading(:,:)    ! fractional fuel loading [0-1]
  real(r8),                          allocatable :: fuel_BD(:)           ! bulk density of fuel [kg/m3]
  real(r8),                          allocatable :: fuel_SAV(:)          ! fuel surface area to volume ratio [/cm]
  real(r8),                          allocatable :: fuel_moisture(:,:)   ! fuel moisture [m3/m3]
  real(r8),                          allocatable :: fuel_MEF(:,:)        ! fuel moisture of extinction [m3/m3
  character(len=100),                allocatable :: fuel_names(:)        ! names of fuel models
  character(len=2),                  allocatable :: carriers(:)          ! carriers of fuel models
  integer                                        :: i, f                 ! looping indices
  integer                                        :: num_fuel_models      ! number of fuel models to test
  
  ! CONSTANTS:
  integer,          parameter :: n_days = 365             ! number of days to run simulation
  character(len=*), parameter :: out_file = 'fuel_out.nc' ! output file 
  
  ! fuel models to test
  integer, parameter, dimension(5) :: fuel_models = (/102, 183, 164, 104, 163/)
  
  ! number of fuel models to test
  num_fuel_models = size(fuel_models)
  
  ! allocate arrays
  allocate(temp_degC(n_days))
  allocate(precip(n_days))
  allocate(rh(n_days))
  allocate(wind(n_days))
  allocate(NI(n_days))
  allocate(fuel_moisture(n_days, num_fuel_models))
  allocate(fuel_MEF(n_days, num_fuel_models))
  allocate(fuel_loading(num_fuel_classes, num_fuel_models))
  allocate(frac_loading(num_fuel_classes, num_fuel_models))
  allocate(fuel_BD(num_fuel_models))
  allocate(fuel_SAV(num_fuel_models))
  allocate(non_trunk_loading(num_fuel_models))
  allocate(fuel_names(num_fuel_models))
  allocate(carriers(num_fuel_models))
  
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
  call fuel_models_array%GetFuelModels()
  do f = 1, num_fuel_models
    
    ! uses data from fuel_models to initialize fuel
    call SetUpFuel(fuel(f), fuel_models_array, fuel_models(f), fuel_names(f), carriers(f))
    
    ! sum up fuel and calculate loading
    call fuel(f)%SumLoading()
    call fuel(f)%CalculateFractionalLoading()
    
    ! calculate geometric properties
    call fuel(f)%AverageBulkDensity_NoTrunks(SF_val_FBD)
    call fuel(f)%AverageSAV_NoTrunks(SF_val_SAV)
    
    ! save values
    fuel_loading(:,f) = fuel(f)%loading(:)
    non_trunk_loading(f) = fuel(f)%non_trunk_loading
    frac_loading(:,f) = fuel(f)%frac_loading(:)
    fuel_BD(f) = fuel(f)%bulk_density_notrunks
    fuel_SAV(f) = fuel(f)%SAV_notrunks
        
  end do
  
  ! run on time steps
  do i = 1, n_days
    call fireWeather%UpdateIndex(temp_degC(i), precip(i), rh(i), wind(i))
    NI(i) = fireWeather%fire_weather_index
    
    ! calculate fuel moisture [m3/m3]
    do f = 1, num_fuel_models
      call fuel(f)%UpdateFuelMoisture(SF_val_SAV, SF_val_drying_ratio, fireWeather)
      fuel_moisture(i, f) = fuel(f)%average_moisture_notrunks
      fuel_MEF(i, f) = fuel(f)%MEF_notrunks
    end do
  end do 
  
  ! write out data
  call WriteFireData(out_file, n_days, num_fuel_models, temp_degC, precip, rh, NI,       &
    fuel_loading, frac_loading, fuel_BD, fuel_SAV, non_trunk_loading, fuel_moisture,         &
    fuel_MEF, fuel_models, carriers)
  
end program FatesTestFuel
