module SFNesterovMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  use SFFireWeatherMod,  only : fire_weather
  
  implicit none
  private

  type, public, extends(fire_weather) :: nesterov_index 

    contains 

      procedure, public :: Init => init_nesterov_fire_weather
      procedure, public :: UpdateIndex => update_nesterov_index
      procedure         :: calc_nesterov_index 

  end type nesterov_index

  real(r8), parameter :: min_precip_thresh = 3.0_r8 ! threshold for precipitation above which to 0.0 NI

  contains 

    subroutine init_nesterov_fire_weather(this)
      !
      !  DESCRIPTION:
      !  Initializes class attributes
      
      ! ARGUMENTS
      class(nesterov_index), intent(inout) :: this ! nesterov index extended class

      ! initialize values to 0.0
      this%fire_weather_index   = 0.0_r8

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

      if (precip > min_precip_thresh) then ! rezero NI if it rains
        this%fire_weather_index = 0.0_r8
      else 
        ! Accumulate Nesterov index over fire season. 
        this%fire_weather_index = this%fire_weather_index + &
          this%calc_nesterov_index(temp_C, precip, rh)
      end if 

    end subroutine update_nesterov_index

    !-------------------------------------------------------------------------------------

    real(r8) function calc_nesterov_index(this, temp_C, precip, rh)
      !
      !  DESCRIPTION:
      !  Calculates current day's Nesterov Index for given input values

      use SFParamsMod,       only : SF_val_fdi_a, SF_val_fdi_b
      use FatesConstantsMod, only : nearzero
      
      ! ARGUMENTS:
      class(nesterov_index), intent(in) :: this   ! nesterov index extended class
      real(r8),              intent(in) :: temp_C ! daily averaged temperature [degrees C]
      real(r8),              intent(in) :: precip ! daily precipitation [mm]
      real(r8),              intent(in) :: rh     ! daily relative humidity [%]

      ! LOCALS:
      real(r8) :: yipsolon ! intermediate variable for dewpoint calculation
      real(r8) :: dewpoint ! dewpoint

      ! error checking, if temperature equals -1.0*SF_val_fdi_b parameter, we get a divide
      ! by zero
      if (abs(SF_val_fdi_b + temp_C) < nearzero) then
        write(fates_log(), *) 'SF_val_fdi_b: (', SF_val_fdi_b, ') + temp_C: (', temp_C, ') == 0.0 - divide by zero imminent!'
        write(fates_log(), *) 'SF_val_fdi_b should be updated using the parameter file.'
        write(fates_log(), *) 'Otherwise check values for temp_C'
        call endrun(msg=errMsg(sourcefile, __LINE__))
      end if 

      if (precip > min_precip_thresh) then ! NI is 0.0 if it rainsf
        calc_nesterov_index = 0.0_r8
      else 
        ! Calculate dewpoint temperature 
        yipsolon = (SF_val_fdi_a*temp_C)/(SF_val_fdi_b + temp_C) + log(max(1.0_r8, rh)/100.0_r8)

        ! error checking, if SF_val_fdi_a parameter - ypisolon == 0, we get a divide
        ! by zero
        if (abs(SF_val_fdi_a - yipsolon) < nearzero) then
          write(fates_log(), *) 'SF_val_fdi_a: (', SF_val_fdi_a, ') - yipsolon: (', yipsolon, ') == 0.0 - divide by zero imminent!'
          write(fates_log(), *) 'SF_val_fdi_a should be updated using the parameter file.'
          write(fates_log(), *) 'Otherwise check values for yipsolon'
          call endrun(msg=errMsg(sourcefile, __LINE__))
        end if 
        dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a - yipsolon) 
        
        ! Nesterov 1968.  Eq 5, Thonicke et al. 2010
        calc_nesterov_index = (temp_C - dewpoint)*temp_C 
        if (calc_nesterov_index < 0.0_r8) calc_nesterov_index = 0.0_r8 ! can't be negative
      endif

    end function calc_nesterov_index

    !-------------------------------------------------------------------------------------

end module SFNesterovMod