module SFNesterovMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  use SFFireWeatherMod,  only : fire_weather
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  
  implicit none
  private

  type, public, extends(fire_weather) :: nesterov_index 

    contains 

      procedure, public :: Init => init_nesterov_fire_weather
      procedure, public :: UpdateIndex => update_nesterov_index
            
  end type nesterov_index

  real(r8), parameter :: min_precip_thresh = 3.0_r8 ! threshold for precipitation above which to zero NI [mm/day]

  contains 

    subroutine init_nesterov_fire_weather(this)
      !
      !  DESCRIPTION:
      !  Initializes class attributes
      
      ! ARGUMENTS
      class(nesterov_index), intent(inout) :: this ! nesterov index extended class

      ! initialize values to 0.0
      this%fire_weather_index   = 0.0_r8
      this%effective_windspeed  = 0.0_r8

    end subroutine init_nesterov_fire_weather

    !-------------------------------------------------------------------------------------

    subroutine update_nesterov_index(this, temp_C, precip, rh, wind)
      !
      !  DESCRIPTION:
      !  Updates Nesterov Index
      
      ! ARGUMENTS
      class(nesterov_index), intent(inout) :: this   ! nesterov index extended class
      real(r8),              intent(in)    :: temp_C ! daily averaged temperature [degrees C]
      real(r8),              intent(in)    :: precip ! daily precipitation [mm]
      real(r8),              intent(in)    :: rh     ! daily relative humidity [%]
      real(r8),              intent(in)    :: wind   ! daily wind speed [m/min]
      
      ! LOCALS:
      real(r8) :: t_dew ! dewpoint temperature [degrees C]

      if (precip > min_precip_thresh) then ! rezero NI if it rains
        this%fire_weather_index = 0.0_r8
      else 
        
        ! Calculate dewpoint temperature
        t_dew = dewpoint(temp_c, rh)
        
        ! Accumulate Nesterov index over fire season. 
        this%fire_weather_index = this%fire_weather_index + calc_nesterov_index(temp_C, t_dew)
      end if 

    end subroutine update_nesterov_index

    !-------------------------------------------------------------------------------------

    real(r8) function calc_nesterov_index(temp_C, t_dew)
      !
      !  DESCRIPTION:
      !  Calculates current day's Nesterov Index for a given input values
      
      ! ARGUMENTS:
      real(r8), intent(in) :: temp_C ! daily averaged temperature [degrees C]
      real(r8), intent(in) :: t_dew  ! daily dewpoint temperature [degrees C]
      
      ! Nesterov 1968. Eq 5, Thonicke et al. 2010
      calc_nesterov_index = (temp_C - t_dew)*temp_C 
      if (calc_nesterov_index < 0.0_r8) calc_nesterov_index = 0.0_r8 ! can't be negative

    end function calc_nesterov_index

    !-------------------------------------------------------------------------------------
    
    real(r8) function dewpoint(temp_C, rh)
      !
      !  DESCRIPTION:
      !  Calculates dewpoint from input air temperature and relative humidity
      !  Uses Equation 8 from Lawrence 2005, https://doi.org/10.1175/BAMS-86-2-225
    
      use FatesConstantsMod, only : dewpoint_a, dewpoint_b
      
      ! ARGUMENTS
      real(r8), intent(in) :: temp_C ! temperature [degrees C]
      real(r8), intent(in) :: rh     ! relative humidity [%]
      
      ! LOCALS
      real(r8) :: yipsolon ! intermediate value for dewpoint calculation
    
      yipsolon = log(max(1.0_r8, rh)/100.0_r8) + (dewpoint_a*temp_C)/(dewpoint_b + temp_C)
      dewpoint = (dewpoint_b*yipsolon)/(dewpoint_a - yipsolon) 
    
    end function dewpoint
    
    !-------------------------------------------------------------------------------------
end module SFNesterovMod