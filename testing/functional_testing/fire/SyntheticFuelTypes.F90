module SyntheticFuelTypes
  
  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private
  
  integer,  parameter :: chunk_size = 10
  real(r8), parameter :: ustons_to_kg = 907.185_r8
  real(r8), parameter :: acres_to_m2 = 4046.86_r8
  real(r8), parameter :: ft_to_m = 0.3048_r8
  
  ! holds data for fake fuel models that can be used for functional
  ! testing of the FATES fire model
  ! these are taken from the fire behavior fuel models in Scott & Burgan 2005
  type, public :: synthetic_fuel_type
    
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
  
  end type synthetic_fuel_type
    
  ! --------------------------------------------------------------------------------------
  
  ! a class to just hold an array of these fuel models
  type, public :: fuel_types_array_class
  
    type(synthetic_fuel_type), allocatable :: fuel_types(:)  ! array of fuel models
    integer                                :: num_fuel_types ! number of total fuel models
    
    contains 
      
      procedure :: AddFuelModel
      procedure :: GetFuelModels
      procedure :: FuelModelPosition
  
  end type fuel_types_array_class
  
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
    class(synthetic_fuel_type), intent(inout) :: this
    integer,                    intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),           intent(in)    :: carrier            ! main carrier
    character(len=*),           intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                   intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                   intent(in)    :: hr1_loading        ! loading for 1-hr fuels [tons/acre]
    real(r8),                   intent(in)    :: hr10_loading       ! loading for 10-hr fuels [tons/acre]
    real(r8),                   intent(in)    :: hr100_loading      ! loading for 100-hr fuels [tons/acre]
    real(r8),                   intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [tons/acre]
    real(r8),                   intent(in)    :: live_woody_loading ! loading for live woody fuels [tons/acre]
    real(r8),                   intent(in)    :: fuel_depth         ! fuel bed depth [ft]
    
    this%fuel_model_index = fuel_model_index
    this%carrier = carrier 
    this%fuel_model_name = fuel_model_name
    this%wind_adj_factor = wind_adj_factor
    this%hr1_loading = hr1_loading*ustons_to_kg/acres_to_m2*0.45_r8 ! convert to kgC/m2
    this%hr10_loading = hr10_loading*ustons_to_kg/acres_to_m2*0.45_r8  ! convert to kgC/m2
    this%hr100_loading = hr100_loading*ustons_to_kg/acres_to_m2*0.45_r8  ! convert to kgC/m2
    this%live_herb_loading = live_herb_loading*ustons_to_kg/acres_to_m2*0.45_r8  ! convert to kgC/m2
    this%live_woody_loading = live_woody_loading*ustons_to_kg/acres_to_m2*0.45_r8  ! convert to kgC/m2
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
    class(fuel_types_array_class), intent(inout) :: this               ! array of fuel models
    integer,                       intent(in)    :: fuel_model_index   ! fuel model index
    character(len=2),              intent(in)    :: carrier            ! main carrier
    character(len=*),              intent(in)    :: fuel_model_name    ! fuel model long name
    real(r8),                      intent(in)    :: wind_adj_factor    ! wind adjustment factor
    real(r8),                      intent(in)    :: hr1_loading        ! loading for 1-hr fuels [tons/acre]
    real(r8),                      intent(in)    :: hr10_loading       ! loading for 10-hr fuels [tons/acre]
    real(r8),                      intent(in)    :: hr100_loading      ! loading for 100-hr fuels [tons/acre]
    real(r8),                      intent(in)    :: live_herb_loading  ! loading for live herbacious fuels [tons/acre]
    real(r8),                      intent(in)    :: live_woody_loading ! loading for live woody fuels [tons/acre]
    real(r8),                      intent(in)    :: fuel_depth         ! fuel bed depth [ft]
    
    ! LOCALS:
    type(synthetic_fuel_type)              :: fuel_model         ! fuel model
    type(synthetic_fuel_type), allocatable :: temporary_array(:) ! temporary array to hold data while re-allocating
    
    ! first make sure we have enough space in the array
    if (allocated(this%fuel_types)) then
      ! already allocated to some size
      if (this%num_fuel_types == size(this%fuel_types)) then 
        ! need to add more space
        allocate(temporary_array(size(this%fuel_types) + chunk_size))
        temporary_array(1:size(this%fuel_types)) = this%fuel_types
        call move_alloc(temporary_array, this%fuel_types)
      end if 
      
      this%num_fuel_types = this%num_fuel_types + 1
  
    else 
      ! first element in array 
      allocate(this%fuel_types(chunk_size))
      this%num_fuel_types = 1
    end if 
    
    call fuel_model%InitFuelModel(fuel_model_index, carrier, fuel_model_name,            &
      wind_adj_factor, hr1_loading, hr10_loading, hr100_loading, live_herb_loading,      &
      live_woody_loading, fuel_depth)
    
    this%fuel_types(this%num_fuel_types) = fuel_model
      
  end subroutine AddFuelModel
  
  ! --------------------------------------------------------------------------------------
  
  integer function FuelModelPosition(this, fuel_model_index)
    !
    ! DESCRIPTION:
    ! Returns the index of a desired fuel model
    !
    
    ! ARGUMENTS:
    class(fuel_types_array_class), intent(in)  :: this             ! array of fuel models
    integer,                       intent(in)  :: fuel_model_index ! desired fuel model index
        
    ! LOCALS:
    integer :: i ! looping index 
    
    do i = 1, this%num_fuel_types
      if (this%fuel_types(i)%fuel_model_index == fuel_model_index) then
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
    class(fuel_types_array_class), intent(inout) :: this ! array of fuel models
    
    call this%AddFuelModel(fuel_model_index=1, carrier='GR', fuel_model_name='short grass',   &
      wind_adj_factor=0.36_r8, hr1_loading=0.7_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8, &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=2, carrier='GR', fuel_model_name='timber and grass understory', &
      wind_adj_factor=0.36_r8, hr1_loading=2.0_r8, hr10_loading=1.0_r8, hr100_loading=0.5_r8,               &
      live_herb_loading=0.5_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=101, carrier='GR', fuel_model_name='short, sparse dry climate grass', &
      wind_adj_factor=0.31_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                     &
      live_herb_loading=0.3_r8, live_woody_loading=0.0_r8, fuel_depth=0.4_r8)
      
    call this%AddFuelModel(fuel_model_index=102, carrier='GR', fuel_model_name='low load dry climate grass',  &
      wind_adj_factor=0.36_r8, hr1_loading=0.1_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                 &
      live_herb_loading=1.0_r8, live_woody_loading=0.0_r8, fuel_depth=1.0_r8)
      
    call this%AddFuelModel(fuel_model_index=183, carrier='TL', fuel_model_name='moderate load conifer litter', &
      wind_adj_factor=0.29_r8, hr1_loading=0.5_r8, hr10_loading=2.2_r8, hr100_loading=2.8_r8,                  &
      live_herb_loading=0.0_r8, live_woody_loading=0.0_r8, fuel_depth=0.3_r8)
      
    call this%AddFuelModel(fuel_model_index=164, carrier='TU', fuel_model_name='dwarf conifer with understory', &
      wind_adj_factor=0.32_r8, hr1_loading=4.5_r8, hr10_loading=0.0_r8, hr100_loading=0.0_r8,                   &
      live_herb_loading=0.0_r8, live_woody_loading=2.0_r8, fuel_depth=0.5_r8)
      
  end subroutine GetFuelModels
  
  ! --------------------------------------------------------------------------------------
    
end module SyntheticFuelTypes