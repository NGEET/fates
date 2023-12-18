module SFNesterovMod

  use FatesConstantsMod, only : r8 => fates_r8
  use SFFireWeatherMod,  only : fire_weather
  use SFParamsMod,       only : min_precip_thresh
  
  implicit none
  private

  type, public, extends(fire_weather) :: nesterov_index 

    contains 

      procedure, public :: Calculate => update_nesterov_index
      procedure         :: calc_nesterov_index 

  end type nesterov_index

  contains 

    subroutine update_nesterov_index(this, temp_C, precip, rh)
      ! ARGUMENTS
      class(nesterov_index), intent(inout) :: this   ! nesterov index extended class
      real(r8),              intent(in)    :: temp_C ! daily averaged temperature [degrees C]
      real(r8),              intent(in)    :: precip ! daily precipitation [mm]
      real(r8),              intent(in)    :: rh     ! daily relative humidity [rh]

      if (precip > min_precip_thresh) then ! rezero NI if it rains
        this%fire_weather_index = 0.0_r8
      else 
        ! Accumulate Nesterov index over fire season. 
        this%fire_weather_index = this%fire_weather_index + &
          this%calc_nesterov_index(temp_C, precip, rh)
      end if 

    end subroutine update_nesterov_index

    !---------------------------------------------------------------------------------------

    real(r8) function calc_nesterov_index(this, temp_C, precip, rh)
      !
      !  DESCRIPTION:
      !  Calculates current day's Nesterov Index for a given input values

      use SFParamsMod, only : SF_val_fdi_a, SF_val_fdi_b, min_precip_thresh
      
      ! ARGUMENTS:
      class(nesterov_index), intent(in) :: this   ! nesterov index extended class
      real(r8),              intent(in) :: temp_C ! daily averaged temperature [degrees C]
      real(r8),              intent(in) :: precip ! daily precipitation [mm]
      real(r8),              intent(in) :: rh     ! daily relative humidity [rh]

      ! LOCALS:
      real(r8) :: yipsolon ! intermediate variable for dewpoint calculation
      real(r8) :: dewpoint ! dewpoint

      if (precip > min_precip_thresh) then ! NI is 0.0 if it rains
        calc_nesterov_index = 0.0_r8
      else 
        ! Calculate dewpoint temperature 
        yipsolon = (SF_val_fdi_a*temp_C)/(SF_val_fdi_b + temp_C) + log(max(1.0_r8, rh)/100.0_r8) 
        dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a - yipsolon) 
        
        ! Nesterov 1968.  Eq 5, Thonicke et al. 2010
        calc_nesterov_index = (temp_C - dewpoint)*temp_C 
        if (calc_nesterov_index < 0.0_r8) calc_nesterov_index = 0.0_r8 ! can't be negative
      endif

    end function calc_nesterov_index

end module SFNesterovMod
