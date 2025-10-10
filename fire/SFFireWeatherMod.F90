module SFFireWeatherMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : ifalse

  implicit none
  private 

  type, abstract, public :: fire_weather

    real(r8) :: fire_weather_index   ! fire weather index
    real(r8) :: effective_windspeed  ! effective wind speed, corrected for by tree/grass cover [m/min]
    integer  :: rx_flag              ! prescribed fire burn window flag [1=burn window present; 0=no burn window]

    contains 

    procedure(initialize_fire_weather), public, deferred :: Init
    procedure(update_fire_weather),     public, deferred :: UpdateIndex
    procedure,                          public           :: UpdateEffectiveWindSpeed
    procedure,                          public           :: UpdateRxfireBurnWindow

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

  subroutine UpdateRxfireBurnWindow(this, rxfire_switch, temp_C, rh, wind, temp_up,      &
    temp_low,rh_up, rh_low, wind_up, wind_low)

    ! ARGUMENTS
    class(fire_weather), intent(inout) :: this           ! fire weather class
    real(r8),            intent(in)    :: temp_C         ! daily averaged temperature [degrees C]
    integer,             intent(in)    :: rxfire_switch  ! whether prescribed fire is turned on  
    real(r8),            intent(in)    :: rh             ! daily relative humidity [%]
    real(r8),            intent(in)    :: wind           ! wind speed [m/min]
    real(r8),            intent(in)    :: temp_up        ! user defined upper bound for temp when define a burn window
    real(r8),            intent(in)    :: temp_low       ! user defined lower bound for temp when define a burn window
    real(r8),            intent(in)    :: rh_up          ! user defined upper bound for relative humidity
    real(r8),            intent(in)    :: rh_low         ! user defined lower bound for relative humidity
    real(r8),            intent(in)    :: wind_up        ! user defined upper bound for wind speed
    real(r8),            intent(in)    :: wind_low       ! user defined lower bound for wind speed

    ! LOCAL VARIABLES
    real(r8) :: t_check  ! intermediate value derived from temp condition check
    real(r8) :: rh_check ! intermediate value derived from RH condition check
    real(r8) :: ws_check ! intermediate value derived from wind speed condition check

    if (rxfire_switch .eq. ifalse) return   
    
    ! check if ambient temperature, relative humidity, and wind speed
    ! are within user defined ranges by comparing current weather
    ! condition to the lower and upper bounds defined. when within range,
    ! it should result in negative value or zero (at the boundary condition)
    ! for each check below 

    t_check = (temp_C - temp_low)*(temp_C - temp_up)
    rh_check = (rh - rh_low)*(rh - rh_up)
    ws_check = (wind - wind_low)*(wind - wind_up)

    if (t_check <= 0.0_r8 .and. rh_check <= 0.0_r8 .and. ws_check <= 0.0_r8) then
      this%rx_flag = 1
    else
      this%rx_flag = 0
    end if

  end subroutine UpdateRxfireBurnWindow
   
end module SFFireWeatherMod
