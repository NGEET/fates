module FatesInterfaceParametersMod
  
  implicit none
  private

  ! Registry keys parameters
  character(len=*), parameter, public :: hlm_fates_decomp_max = 'max_layers_decomp'
  character(len=*), parameter, public :: hlm_fates_soil_level = 'soil_level_number'
  character(len=*), parameter, public :: hlm_fates_decomp_frac_moisture = 'decomp_frac_moisture'
  character(len=*), parameter, public :: hlm_fates_decomp_frac_temperature = 'decomp_frac_temperature'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_cellulose = 'litter_carbon_cellulose'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_lignin = 'litter_carbon_lignin'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_labile = 'litter_carbon_labile'
  
  ! Registry update frequency parameters
  integer, parameter :: registry_update_init = 1       ! variable only needs to be updated during initialization
  integer, parameter :: registry_update_daily = 2      ! variable needs to be updated daily
  integer, parameter :: registry_update_timestep = 3   ! variable needs to be updated at each timestep
  
end module FatesInterfaceParametersMod