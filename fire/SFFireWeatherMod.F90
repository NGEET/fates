module SFFireWeatherMod

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private 

  type, abstract, public :: fire_weather
    real(r8) :: fire_weather_index   ! fire weather index
    contains 
    procedure(initialize_fire_weather), public, deferred :: Init
    procedure(update_fire_weather),     public, deferred :: Update
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

end module SFFireWeatherMod