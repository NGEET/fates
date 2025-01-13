module SFFireWeatherMod

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private 

  type, abstract, public :: fire_weather

    real(r8) :: fire_weather_index   ! fire weather index
    real(r8) :: effective_windspeed  ! effective wind speed, corrected for by tree/grass cover [m/min]

    contains 

    procedure(initialize_fire_weather), public, deferred :: Init
    procedure(update_fire_weather),     public, deferred :: UpdateIndex
    procedure,                          public           :: UpdateEffectiveWindSpeed

  end type fire_weather

  abstract interface
    subroutine initialize_fire_weather(this)

      import :: fire_weather 

      class(fire_weather), intent(inout) :: this

    end subroutine initialize_fire_weather

    subroutine update_fire_weather(this, temp_C, precip, rh, wind)

      use FatesConstantsMod, only : r8 => fates_r8

      import :: fire_weather 

      class(fire_weather), intent(inout) :: this
      real(r8),            intent(in)    :: temp_C
      real(r8),            intent(in)    :: precip
      real(r8),            intent(in)    :: rh
      real(r8),            intent(in)    :: wind
      
    end subroutine update_fire_weather
    
   end interface

  contains 

  subroutine UpdateEffectiveWindSpeed(this, wind_speed, tree_fraction, grass_fraction,   &
      bare_fraction)
    !
    !  DESCRIPTION:
    !  Calculates effective wind speed

    ! CONSTANTS:
    real(r8), parameter :: wind_atten_treed = 0.4_r8 ! wind attenuation factor for tree fraction
    real(r8), parameter :: wind_atten_grass = 0.6_r8 ! wind attenuation factor for grass fraction
      
    ! ARGUMENTS
    class(fire_weather), intent(inout) :: this           ! fire weather class
    real(r8),            intent(in)    :: wind_speed     ! wind speed [m/min]
    real(r8),            intent(in)    :: tree_fraction  ! tree fraction [0-1]
    real(r8),            intent(in)    :: grass_fraction ! grass fraction [0-1]
    real(r8),            intent(in)    :: bare_fraction  ! bare ground fraction [0-1]

    this%effective_windspeed = wind_speed*(tree_fraction*wind_atten_treed +              &
      (grass_fraction + bare_fraction)*wind_atten_grass)

  end subroutine UpdateEffectiveWindSpeed

end module SFFireWeatherMod