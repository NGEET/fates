module FatesInterfaceParametersMod
  
  implicit none
  private

  ! Registry key parameters
  character(len=*), parameter, public :: hlm_fates_decomp= 'decomp_layers'
  character(len=*), parameter, public :: hlm_fates_decomp_max = 'max_decomp_layers'
  character(len=*), parameter, public :: hlm_fates_decomp_thickness= 'decomp_thickness'
  character(len=*), parameter, public :: hlm_fates_decomp_id = 'decomp_id'
  character(len=*), parameter, public :: hlm_fates_rooting_max_depth = 'rooting_max_depth'
  character(len=*), parameter, public :: hlm_fates_soil_level = 'soil_level_number'
  character(len=*), parameter, public :: hlm_fates_decomp_frac_moisture = 'decomp_frac_moisture'
  character(len=*), parameter, public :: hlm_fates_decomp_frac_temperature = 'decomp_frac_temperature'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_cellulose = 'litter_carbon_cellulose'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_lignin = 'litter_carbon_lignin'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_labile = 'litter_carbon_labile'
  character(len=*), parameter, public :: hlm_fates_litter_carbon_total= 'litter_carbon_total'
  character(len=*), parameter, public :: hlm_fates_litter_phosphorus_cellulose = 'litter_phosphorus_cellulose'
  character(len=*), parameter, public :: hlm_fates_litter_phosphorus_lignin = 'litter_phosphorus_lignin'
  character(len=*), parameter, public :: hlm_fates_litter_phosphorus_labile = 'litter_phosphorus_labile'
  character(len=*), parameter, public :: hlm_fates_litter_phosphorus_total= 'litter_phosphorus_total'
  character(len=*), parameter, public :: hlm_fates_litter_nitrogen_cellulose = 'litter_nitrogen_cellulose'
  character(len=*), parameter, public :: hlm_fates_litter_nitrogen_lignin = 'litter_nitrogen_lignin'
  character(len=*), parameter, public :: hlm_fates_litter_nitrogen_labile = 'litter_nitrogen_labile'
  character(len=*), parameter, public :: hlm_fates_litter_nitrogen_total= 'litter_nitrogen_total'

  ! Registry update frequency parameters
  integer, parameter :: registry_update_init = 1       ! variable only needs to be updated during initialization
  integer, parameter :: registry_update_daily = 2      ! variable needs to be updated daily
  integer, parameter :: registry_update_timestep = 3   ! variable needs to be updated at each timestep
  
end module FatesInterfaceParametersMod