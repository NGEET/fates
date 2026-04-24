Module FatesRadiationMemMod

  ! ---------------------------------------------------------------------------
  ! This module is a space that holds data that defines how
  ! FATES in particular uses its radiation schemes.
  !
  ! Alternatively, the TwoStreamMLPEMod is more agnostic.
  ! For instance, TwoStreamMLPEMod makes no assumptions about
  ! which or how many broad bands are used
  !
  ! For now, this module also holds relevant data for Norman radiation
  ! ---------------------------------------------------------------------------

  use FatesConstantsMod, only : r8 => fates_r8
  
  integer, parameter, public :: norman_solver = 1
  integer, parameter, public :: twostr_solver = 2

  integer, parameter, public :: num_rad_stream_types = 2    ! The number of radiation streams used (direct/diffuse)

  integer, parameter, public :: idirect   = 1             ! This is the array index for direct radiation
  integer, parameter, public :: idiffuse  = 2             ! This is the array index for diffuse radiation

  
  ! TODO: we use this cp_maxSWb only because we have a static array q(size=2) of
  ! land-ice abledo for vis and nir.  This should be a parameter, which would
  ! get us on track to start using multi-spectral or hyper-spectral (RGK 02-2017)
  
  integer, parameter, public :: num_swb = 2     ! Number of shortwave bands we use
                                                ! This needs to match what is used in the host model
                                                ! This is visible (1) and near-infrared (2)
  
  integer, parameter, public :: ivis = 1        ! This is the array index for short-wave
                                                ! radiation in the visible spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod

  integer, parameter, public :: inir = 2        ! This is the array index for short-wave
                                                ! radiation in the near-infrared spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod

  integer, parameter, public :: ipar = ivis     ! The photosynthetically active band
                                                ! can be approximated to be equal to the visible band

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: alb_ice(num_swb)   =  (/ 0.80_r8, 0.55_r8 /)

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: rho_snow(num_swb) = (/ 0.80_r8, 0.55_r8 /)

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: tau_snow(num_swb) = (/ 0.01_r8, 0.01_r8 /)


  

  
end Module FatesRadiationMemMod
