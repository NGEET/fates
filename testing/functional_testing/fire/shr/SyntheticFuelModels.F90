module SyntheticFuelModels
  
  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private
  
  ! Fuel model numbers come from Scott and Burgen (2005) RMRS-GTR-153

  integer, parameter, public, dimension(52) :: all_fuel_models = (/1, 2, 101, 102, 104,  &
                                                      107, 121, 122, 3, 103, 105, 106,   &
                                                      108, 109, 123, 124, 4, 5, 6, 141,  &
                                                      142, 145, 147, 161, 164, 10, 7,    &
                                                      143, 144, 146, 148, 149, 162,      &
                                                      163, 8, 9, 181, 182, 183, 184,     &
                                                      185, 186, 187, 188, 189, 11, 12,   &
                                                      13, 201, 202, 203, 204/)
  
  integer,  parameter :: chunk_size = 10
  real(r8), parameter :: ustons_to_kg = 907.185_r8
  real(r8), parameter :: acres_to_m2 = 4046.86_r8
  real(r8), parameter :: ustons_acre_to_kgC_m2 = ustons_to_kg/acres_to_m2*0.45_r8
  real(r8), parameter :: ft_to_m = 0.3048_r8
  
  ! holds data for fake fuel models that can be used for functional
  ! testing of the FATES fire model
  ! these are taken from the fire behavior fuel models in Scott & Burgan 2005
  type, public :: synthetic_fuel_model
    
    integer            :: fuel_model_index    ! fuel model index
    character(len=2)   :: carrier             ! carrier ('GR', 'GS', etc.)
    character(len=5)   :: fuel_model_code     ! carrier plus fuel model
    character(len=100) :: fuel_model_name     ! long name of fuel model
    real(r8)           :: wind_adj_factor     ! wind adjustment factor
    real(r8)           :: hr1_loading         ! fuel loading for 1 hour fuels [kg/m2]
    real(r8)           :: hr10_loading        ! fuel loading for 10 hour fuels [kg/m2]
    real(r8)           :: hr100_loading       ! fuel loading for 100 hour fuels [kg/m2]
    real(r8)           :: live_herb_loading   ! fuel loading for live herbacious fuels [kg/m2]
    real(r8)           :: live_woody_loading  ! fuel loading for live woody fuels [kg/m2]
    real(r8)           :: fuel_depth          ! fuel bed depth [m]
    contains 
    
      procedure :: InitFuelModel
  
  end type synthetic_fuel_model
    
  ! --------------------------------------------------------------------------------------
  
  ! a class to just hold an array of these fuel models
  type, public :: fuel_models_array_class
  
    type(synthetic_fuel_model), allocatable :: fuel_models(:)  ! array of fuel models
    integer                                 :: num_fuel_models ! number of total fuel models
    
    contains 
      
      procedure :: AddFuelModel
      procedure :: GetFuelModels
      procedure :: FuelModelPosition
  
  end type fuel_models_array_class
  
  ! --------------------------------------------------------------------------------------
  
  contains 
    
  subroutine InitFuelModel(this, fuel_model_index, carrier, fuel_model_name,             &
    wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,        &
    live_woody_loading, fuel_depth)
    !
    ! DESCRIPTION:
    ! Initializes the fuel model with input characteristics
    ! Also converts units as needed
    !
    ! NOTE THE UNITS ON INPUTS
    !
    
    ! ARGUMENTS:
    class(synthetic_fuel_model), intent(inout) :: this
    integer,                     intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),            intent(in)    :: carrier            ! main carrier
    character(len=*),            intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                    intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                    intent(in)    :: hr1_loading        ! loading for 1-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: hr10_loading       ! loading for 10-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: hr100_loading      ! loading for 100-hr fuels [US tons/acre]
    real(r8),                    intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [US tons/acre]
    real(r8),                    intent(in)    :: live_woody_loading ! loading for live woody fuels [US tons/acre]
    real(r8),                    intent(in)    :: fuel_depth         ! fuel bed depth [ft]
        
    this%fuel_model_index = fuel_model_index
    this%carrier = carrier 
    this%fuel_model_name = fuel_model_name
    this%wind_adj_factor = wind_adj_factor
    this%hr1_loading = hr1_loading*ustons_acre_to_kgC_m2 ! convert to kgC/m2
    this%hr10_loading = hr10_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%hr100_loading = hr100_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%live_herb_loading = live_herb_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%live_woody_loading = live_woody_loading*ustons_acre_to_kgC_m2  ! convert to kgC/m2
    this%fuel_depth = fuel_depth*ft_to_m ! convert to m
      
  end subroutine InitFuelModel
  
  ! --------------------------------------------------------------------------------------
    
  subroutine AddFuelModel(this, fuel_model_index, carrier, fuel_model_name,              &
    wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,        &
    live_woody_loading, fuel_depth)
    !
    ! DESCRIPTION:
    ! Adds a fuel model to the dynamic array
    !
    ! NOTE THE UNITS ON INPUTS
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(inout) :: this               ! array of fuel models
    integer,                        intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),               intent(in)    :: carrier            ! main carrier
    character(len=*),               intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                       intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                       intent(in)    :: hr1_loading        ! loading for 1-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: hr10_loading       ! loading for 10-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: hr100_loading      ! loading for 100-hr fuels [US tons/acre]
    real(r8),                       intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [US tons/acre]
    real(r8),                       intent(in)    :: live_woody_loading ! loading for live woody fuels [US tons/acre]
    real(r8),                       intent(in)    :: fuel_depth         ! fuel bed depth [ft]
    
    ! LOCALS:
    type(synthetic_fuel_model)              :: fuel_model         ! fuel model
    type(synthetic_fuel_model), allocatable :: temporary_array(:) ! temporary array to hold data while re-allocating
    
    ! first make sure we have enough space in the array
    if (allocated(this%fuel_models)) then
      ! already allocated to some size
      if (this%num_fuel_models == size(this%fuel_models)) then 
        ! need to add more space
        allocate(temporary_array(size(this%fuel_models) + chunk_size))
        temporary_array(1:size(this%fuel_models)) = this%fuel_models
        call move_alloc(temporary_array, this%fuel_models)
      end if 
      
      this%num_fuel_models = this%num_fuel_models + 1
  
    else 
      ! first element in array 
      allocate(this%fuel_models(chunk_size))
      this%num_fuel_models = 1
    end if 
    
    call fuel_model%InitFuelModel(fuel_model_index, carrier, fuel_model_name,            &
      wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,      &
      live_woody_loading, fuel_depth)
    
    this%fuel_models(this%num_fuel_models) = fuel_model
      
  end subroutine AddFuelModel
  
  ! --------------------------------------------------------------------------------------
  
  integer function FuelModelPosition(this, fuel_model_index)
    !
    ! DESCRIPTION:
    ! Returns the index of a desired fuel model
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(in)  :: this             ! array of fuel models
    integer,                        intent(in)  :: fuel_model_index ! desired fuel model index
        
    ! LOCALS:
    integer :: i ! looping index 
    
    do i = 1, this%num_fuel_models
      if (this%fuel_models(i)%fuel_model_index == fuel_model_index) then
        FuelModelPosition = i
        return
      end if
    end do
    write(*, '(a, i2, a)') "Cannot find the fuel model index ", fuel_model_index, "."
    stop
  
  end function FuelModelPosition
  
  ! --------------------------------------------------------------------------------------
  
  subroutine GetFuelModels(this)
    !
    ! DESCRIPTION:
    ! Returns an array of hard-coded fuel models
    ! these are taken from the fire behavior fuel models in Scott & Burgan 2005
    !
    
    ! ARGUMENTS:
    class(fuel_models_array_class), intent(inout) :: this ! array of fuel models
    
    call this%AddFuelModel(fuel_model_index=1, carrier='GR', fuel_model_name='short grass',   &
      wind_adj_factor=0.36_r8, hr1_loading=0.7_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=2, carrier='GR', fuel_model_name='timber and grass understory', &
      wind_adj_factor=0.36_r8, hr1_loading=2.0_r8, hr10_loading=1.0_r8, hr100_loading=0.5_r8,               &
      live_herb_loading=0.5_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=3, carrier='GR', fuel_model_name='tall grass',    &
      wind_adj_factor=0.44_r8, hr1_loading=3.0_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.5_r8)
      
    call this%AddFuelModel(fuel_model_index=4, carrier='SH', fuel_model_name='chapparal',     &
      wind_adj_factor=0.55_r8, hr1_loading=5.0_r8, hr10_loading=4.0_r8, hr100_loading=2.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=5.0_r8, fuel_depth=6.0_r8)
   
    call this%AddFuelModel(fuel_model_index=5, carrier='SH', fuel_model_name='brush',         &
      wind_adj_factor=0.42_r8, hr1_loading=1.0_r8, hr10_loading=0.5_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=2.0_r8)
    
    call this%AddFuelModel(fuel_model_index=6, carrier='SH', fuel_model_name='dormant brush', &
      wind_adj_factor=0.44_r8, hr1_loading=1.5_r8, hr10_loading=2.5_r8, hr100_loading=2.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.5_r8)

    call this%AddFuelModel(fuel_model_index=7, carrier='SH', fuel_model_name='southern rough', &
      wind_adj_factor=0.44_r8, hr1_loading=1.1_r8, hr10_loading=1.9_r8, hr100_loading=1.0_r8,  &
      live_herb_loading=0.0_r8, live_woody_loading=0.4_r8, fuel_depth=2.5_r8)
      
    call this%AddFuelModel(fuel_model_index=8, carrier='TL', fuel_model_name='compact timber litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.5_r8, hr10_loading=1.0_r8, hr100_loading=2.5_r8,         &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8)
      
    call this%AddFuelModel(fuel_model_index=9, carrier='TL', fuel_model_name='hardwood litter', &
      wind_adj_factor=0.28_r8, hr1_loading=2.9_r8, hr10_loading=0.4_r8, hr100_loading=0.2_r8,   &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8)
    
    call this%AddFuelModel(fuel_model_index=10, carrier='TU', fuel_model_name='timber and litter understorey', &
      wind_adj_factor=0.46_r8, hr1_loading=3.0_r8, hr10_loading=2.0_r8, hr100_loading=5.0_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=11, carrier='SB', fuel_model_name='light slash',  &
      wind_adj_factor=0.36_r8, hr1_loading=1.5_r8, hr10_loading=4.5_r8, hr100_loading=5.5_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=12, carrier='SB', fuel_model_name='medium slash',   &
      wind_adj_factor=0.43_r8, hr1_loading=4.0_r8, hr10_loading=14.0_r8, hr100_loading=16.5_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.3_r8)
      
    call this%AddFuelModel(fuel_model_index=13, carrier='SB', fuel_model_name='heavy slash',    &
      wind_adj_factor=0.46_r8, hr1_loading=7.0_r8, hr10_loading=23.0_r8, hr100_loading=28.1_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=3.0_r8)
    
    call this%AddFuelModel(fuel_model_index=101, carrier='GR', fuel_model_name='short, sparse dry climate grass', &
      wind_adj_factor=0.31_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                     &
      live_herb_loading=0.3_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8)
      
    call this%AddFuelModel(fuel_model_index=102, carrier='GR', fuel_model_name='low load dry climate grass',  &
      wind_adj_factor=0.36_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=1.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=103, carrier='GR', fuel_model_name='low load very coarse humid climate grass',  &
      wind_adj_factor=0.42_r8, hr1_loading=0.1_r8, hr10_loading=0.4_r8, hr100_loading=0.0_r8,                               &
      live_herb_loading=1.5_r8, live_woody_loading=0.0_r8, fuel_depth=2.0_r8)
    
    call this%AddFuelModel(fuel_model_index=104, carrier='GR', fuel_model_name='moderate load dry climate grass',  &
      wind_adj_factor=0.42_r8, hr1_loading=0.3_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                      &
      live_herb_loading=1.9_r8, live_woody_loading=0.0_r8, fuel_depth=2.0_r8)
    
    call this%AddFuelModel(fuel_model_index=105, carrier='GR', fuel_model_name='low load humid climate grass',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.4_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=2.5_r8, live_woody_loading=0.0_r8, fuel_depth=1.5_r8)
    
    call this%AddFuelModel(fuel_model_index=106, carrier='GR', fuel_model_name='moderate load humid climate grass',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                        &
      live_herb_loading=3.4_r8, live_woody_loading=0.0_r8, fuel_depth=1.5_r8)
      
    call this%AddFuelModel(fuel_model_index=107, carrier='GR', fuel_model_name='high load dry climate grass',  &
      wind_adj_factor=0.46_r8, hr1_loading=1.0_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                  &
      live_herb_loading=5.4_r8, live_woody_loading=0.0_r8, fuel_depth=3.0_r8)
    
    call this%AddFuelModel(fuel_model_index=108, carrier='GR', fuel_model_name='high load humid climate grass', &
      wind_adj_factor=0.49_r8, hr1_loading=0.5_r8, hr10_loading=1.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=7.3_r8, live_woody_loading=0.0_r8, fuel_depth=4.0_r8)
      
    call this%AddFuelModel(fuel_model_index=109, carrier='GR', fuel_model_name='very high load humid climate grass-shrub', &
      wind_adj_factor=0.52_r8, hr1_loading=1.0_r8, hr10_loading=1.0_r8, hr100_loading=0.0_r8,                              &
      live_herb_loading=9.0_r8, live_woody_loading=0.0_r8, fuel_depth=5.0_r8)
      
    call this%AddFuelModel(fuel_model_index=121, carrier='GS', fuel_model_name='low load dry climate grass-shrub',  &
      wind_adj_factor=0.35_r8, hr1_loading=0.2_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                       &
      live_herb_loading=0.5_r8, live_woody_loading=0.7_r8, fuel_depth=0.9_r8)
      
    call this%AddFuelModel(fuel_model_index=122, carrier='GS', fuel_model_name='moderate load dry climate grass-shrub',  &
      wind_adj_factor=0.39_r8, hr1_loading=0.5_r8, hr10_loading=0.5_r8, hr100_loading=0.0_r8,                            &
      live_herb_loading=0.6_r8, live_woody_loading=1.0_r8, fuel_depth=1.5_r8)
      
    call this%AddFuelModel(fuel_model_index=123, carrier='GS', fuel_model_name='moderate load humid climate grass-shrub', &
      wind_adj_factor=0.41_r8, hr1_loading=0.3_r8, hr10_loading=0.3_r8, hr100_loading=0.0_r8,                             &
      live_herb_loading=1.5_r8, live_woody_loading=1.3_r8, fuel_depth=1.8_r8)
      
    call this%AddFuelModel(fuel_model_index=124, carrier='GS', fuel_model_name='high load humid climate grass-shrub',  &
      wind_adj_factor=0.42_r8, hr1_loading=1.9_r8, hr10_loading=0.3_r8, hr100_loading=0.1_r8,                          &
      live_herb_loading=3.4_r8, live_woody_loading=7.1_r8, fuel_depth=2.1_r8)
          
    call this%AddFuelModel(fuel_model_index=141, carrier='SH', fuel_model_name='low load dry climate shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=0.3_r8, hr10_loading=0.3_r8, hr100_loading=0.0_r8,                &
      live_herb_loading=0.2_r8, live_woody_loading=1.3_r8, fuel_depth=1.0_r8)
        
    call this%AddFuelModel(fuel_model_index=142, carrier='SH', fuel_model_name='moderate load dry climate shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=1.4_r8, hr10_loading=2.4_r8, hr100_loading=0.8_r8,                     &
      live_herb_loading=0.0_r8, live_woody_loading=3.9_r8, fuel_depth=1.0_r8)
    
    call this%AddFuelModel(fuel_model_index=143, carrier='SH', fuel_model_name='moderate load humid climate shrub', &
      wind_adj_factor=0.44_r8, hr1_loading=0.5_r8, hr10_loading=3.0_r8, hr100_loading=0.0_r8,                       &
      live_herb_loading=0.0_r8, live_woody_loading=6.2_r8, fuel_depth=2.4_r8)
    
    call this%AddFuelModel(fuel_model_index=144, carrier='SH', fuel_model_name='low load humid climate timber-shrub', &
      wind_adj_factor=0.46_r8, hr1_loading=0.9_r8, hr10_loading=1.2_r8, hr100_loading=0.2_r8,                         &
      live_herb_loading=0.0_r8, live_woody_loading=2.6_r8, fuel_depth=3.0_r8)
     
    call this%AddFuelModel(fuel_model_index=145, carrier='SH', fuel_model_name='high load dry climate shrub', &
      wind_adj_factor=0.55_r8, hr1_loading=3.6_r8, hr10_loading=2.1_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=0.0_r8, live_woody_loading=2.9_r8, fuel_depth=6.0_r8)    
    
    call this%AddFuelModel(fuel_model_index=146, carrier='SH', fuel_model_name='low load humid climate shrub', &
      wind_adj_factor=0.42_r8, hr1_loading=2.9_r8, hr10_loading=1.5_r8, hr100_loading=0.0_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=1.4_r8, fuel_depth=2.0_r8)     
      
    call this%AddFuelModel(fuel_model_index=147, carrier='SH', fuel_model_name='very high load dry climate shrub', &
      wind_adj_factor=0.55_r8, hr1_loading=3.5_r8, hr10_loading=5.3_r8, hr100_loading=2.2_r8,                      &
      live_herb_loading=0.0_r8, live_woody_loading=3.4_r8, fuel_depth=6.0_r8)
      
    call this%AddFuelModel(fuel_model_index=148, carrier='SH', fuel_model_name='high load humid climate shrub', &
      wind_adj_factor=0.46_r8, hr1_loading=2.1_r8, hr10_loading=3.4_r8, hr100_loading=0.9_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=4.4_r8, fuel_depth=3.0_r8)  
      
    call this%AddFuelModel(fuel_model_index=149, carrier='SH', fuel_model_name='very high load humid climate shrub', &
      wind_adj_factor=0.5_r8, hr1_loading=4.5_r8, hr10_loading=2.5_r8, hr100_loading=0.0_r8,                         &
      live_herb_loading=1.6_r8, live_woody_loading=7.0_r8, fuel_depth=4.4_r8)  
    
    call this%AddFuelModel(fuel_model_index=161, carrier='TU', fuel_model_name='light load dry climate timber-grass-shrub', &
      wind_adj_factor=0.33_r8, hr1_loading=0.2_r8, hr10_loading=0.9_r8, hr100_loading=1.5_r8,                               &
      live_herb_loading=0.2_r8, live_woody_loading=0.9_r8, fuel_depth=0.6_r8)
        
    call this%AddFuelModel(fuel_model_index=162, carrier='TU', fuel_model_name='moderate load humid climate timber-shrub', &
      wind_adj_factor=0.36_r8, hr1_loading=1.0_r8, hr10_loading=1.8_r8, hr100_loading=1.3_r8,                              &
      live_herb_loading=0.0_r8, live_woody_loading=0.2_r8, fuel_depth=1.0_r8)
  
    call this%AddFuelModel(fuel_model_index=163, carrier='TU', fuel_model_name='moderate load humid climate timber-grass-shrub', &
      wind_adj_factor=0.38_r8, hr1_loading=1.1_r8, hr10_loading=0.2_r8, hr100_loading=0.2_r8,                                    &
      live_herb_loading=0.3_r8, live_woody_loading=0.7_r8, fuel_depth=1.3_r8)
      
    call this%AddFuelModel(fuel_model_index=164, carrier='TU', fuel_model_name='dwarf conifer with understory', &
      wind_adj_factor=0.32_r8, hr1_loading=4.5_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=0.5_r8)
 
    call this%AddFuelModel(fuel_model_index=165, carrier='TU', fuel_model_name='very high load dry climate timber-shrub', &
      wind_adj_factor=0.33_r8, hr1_loading=4.0_r8, hr10_loading=4.0_r8, hr100_loading=3.0_r8,                             &
      live_herb_loading=0.0_r8, live_woody_loading=3.0_r8, fuel_depth=1.0_r8)
    
    call this%AddFuelModel(fuel_model_index=181, carrier='TL', fuel_model_name='low load compact conifer litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.0_r8, hr10_loading=2.2_r8, hr100_loading=3.6_r8,                     &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8)
      
    call this%AddFuelModel(fuel_model_index=182, carrier='TL', fuel_model_name='low load broadleaf litter', &
      wind_adj_factor=0.28_r8, hr1_loading=1.4_r8, hr10_loading=2.3_r8, hr100_loading=2.2_r8,               &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.2_r8)
    
    call this%AddFuelModel(fuel_model_index=183, carrier='TL', fuel_model_name='moderate load conifer litter', &
      wind_adj_factor=0.29_r8, hr1_loading=0.5_r8, hr10_loading=2.2_r8, hr100_loading=2.8_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8)
      
    call this%AddFuelModel(fuel_model_index=184, carrier='TL', fuel_model_name='small downed logs', &
      wind_adj_factor=0.31_r8, hr1_loading=0.5_r8, hr10_loading=1.5_r8, hr100_loading=4.2_r8,       &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8)
      
    call this%AddFuelModel(fuel_model_index=185, carrier='TL', fuel_model_name='high load conifer litter', &
      wind_adj_factor=0.33_r8, hr1_loading=1.2_r8, hr10_loading=2.5_r8, hr100_loading=4.4_r8,              &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.6_r8)
    
    call this%AddFuelModel(fuel_model_index=186, carrier='TL', fuel_model_name='moderate load broadleaf litter', &
      wind_adj_factor=0.29_r8, hr1_loading=2.4_r8, hr10_loading=1.2_r8, hr100_loading=1.2_r8,                    &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8)
      
    call this%AddFuelModel(fuel_model_index=187, carrier='TL', fuel_model_name='large downed logs', &
      wind_adj_factor=0.31_r8, hr1_loading=0.3_r8, hr10_loading=1.4_r8, hr100_loading=8.1_r8,       &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8)
      
    call this%AddFuelModel(fuel_model_index=188, carrier='TL', fuel_model_name='long-needle litter', &
      wind_adj_factor=0.29_r8, hr1_loading=5.0_r8, hr10_loading=1.4_r8, hr100_loading=1.1_r8,        &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8)
      
    call this%AddFuelModel(fuel_model_index=189, carrier='TL', fuel_model_name='very high load broadleaf litter', &
      wind_adj_factor=0.33_r8, hr1_loading=6.7_r8, hr10_loading=3.3_r8, hr100_loading=4.2_r8,                     &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.6_r8)
      
    call this%AddFuelModel(fuel_model_index=201, carrier='SB', fuel_model_name='low load activity fuel', &
      wind_adj_factor=0.36_r8, hr1_loading=1.5_r8, hr10_loading=3.0_r8, hr100_loading=11.1_r8,           &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=202, carrier='SB', fuel_model_name='moderate load activity fuel or low load blowdown', &
      wind_adj_factor=0.36_r8, hr1_loading=4.5_r8, hr10_loading=4.3_r8, hr100_loading=4.0_r8,                                      &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=203, carrier='SB', fuel_model_name='high load activity fuel or moderate load blowdown', &
      wind_adj_factor=0.38_r8, hr1_loading=5.5_r8, hr10_loading=2.8_r8, hr100_loading=3.0_r8,                                       &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.2_r8)
      
    call this%AddFuelModel(fuel_model_index=204, carrier='SB', fuel_model_name='high load blowdown', &
      wind_adj_factor=0.45_r8, hr1_loading=5.3_r8, hr10_loading=3.5_r8, hr100_loading=5.3_r8,        &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=2.7_r8)
      
    end subroutine GetFuelModels
  
  ! --------------------------------------------------------------------------------------
    
end module SyntheticFuelModels
