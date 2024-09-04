module FatesFuelMod

  use FatesFuelClassesMod, only : nfsc, nfsc_notrunks, fuel_classes
  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : nearzero
  use SFNesterovMod,       only : nesterov_index
  use SFFireWeatherMod,    only : fire_weather
  use FatesGlobals,        only : fates_log
  use FatesGlobals,        only : endrun => fates_endrun
  use shr_log_mod,         only : errMsg => shr_log_errMsg

  implicit none
  private

  type, public :: fuel_type
    
    real(r8) :: loading(nfsc_notrunks)      ! fuel loading of non-trunks fuel class [kgC/m2]
    real(r8) :: trunk_loading               ! fuel loading of trunk fuel class [kgC/m2]
    real(r8) :: effective_moisture(nfsc)    ! fuel effective moisture all fuel class (moisture/MEF) [m3/m3]
    real(r8) :: frac_loading(nfsc_notrunks) ! fractional loading of non-trunk fuel classes [0-1] 
    real(r8) :: frac_burnt(nfsc)            ! fraction of litter burnt by fire [0-1]
    real(r8) :: total_loading               ! total fuel loading - DOES NOT INCLUDE TRUNKS [kgC/m2]
    real(r8) :: average_moisture            ! weighted average of fuel moisture across non-trunk fuel classes [m3/m3]
    real(r8) :: bulk_density                ! weighted average of bulk density across non-trunk fuel classes [kg/m3]
    real(r8) :: SAV                         ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
    real(r8) :: MEF                         ! weighted average of moisture of extinction across non-trunk fuel classes [m3/m3]

    contains
      
      procedure :: Init
      procedure :: CalculateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity
      procedure :: AverageSAV
      procedure :: BurnFuel

  end type fuel_type
  
  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! just zero everything
      this%loading(1:nfsc_notrunks) = 0.0_r8
      this%trunk_loading = 0.0_r8
      this%frac_loading(1:nfsc_notrunks) = 0.0_r8
      this%frac_burnt(1:nfsc) = 0.0_r8  
      this%effective_moisture(1:nfsc) = 0.0_r8
      this%total_loading = 0.0_r8
      this%average_moisture = 0.0_r8 
      this%bulk_density = 0.0_r8
      this%SAV = 0.0_r8
      this%MEF = 0.0_r8 

    end subroutine Init 

    !-------------------------------------------------------------------------------------

    subroutine CalculateLoading(this, leaf_litter, twig_litter, small_branch_litter,     &
        large_branch_litter, trunk_litter, live_grass)
      ! DESCRIPTION:
      !   Calculates loading for each fuel type

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                ! fuel class
      real(r8),         intent(in)    :: leaf_litter         ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: twig_litter         ! input twig litter [kgC/m2]
      real(r8),         intent(in)    :: small_branch_litter ! input small branch litter [kgC/m2]
      real(r8),         intent(in)    :: large_branch_litter ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: trunk_litter        ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: live_grass          ! input live grass [kgC/m2]

      this%loading(fuel_classes%dead_leaves()) = leaf_litter
      this%loading(fuel_classes%twigs()) = twig_litter
      this%loading(fuel_classes%small_branches()) = small_branch_litter
      this%loading(fuel_classes%large_branches()) = large_branch_litter
      this%loading(fuel_classes%live_grass()) = live_grass
      this%trunk_loading = trunk_litter

    end subroutine CalculateLoading

    !-------------------------------------------------------------------------------------

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sums up the loading

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      this%total_loading = sum(this%loading(1:nfsc_notrunks))

    end subroutine SumLoading

    !-------------------------------------------------------------------------------------
    
    subroutine CalculateFractionalLoading(this)
      ! DESCRIPTION:
      !   Calculates fractional loading for fuel

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      ! sum up loading just in case
      call this%SumLoading()

      if (this%total_loading > nearzero) then
        this%frac_loading(1:nfsc_notrunks) = this%loading(1:nfsc_notrunks)/this%total_loading
      else 
        this%frac_loading(1:nfsc_notrunks) = 0.0_r8
      end if 

    end subroutine CalculateFractionalLoading

    !-------------------------------------------------------------------------------------
    
    subroutine UpdateFuelMoisture(this, sav_fuel, drying_ratio, fireWeatherClass)
      ! DESCRIPTION:
      !   Updates fuel moisture depending on what fire weather class is in use

      ! ARGUMENTS:
      class(fuel_type),    intent(inout) :: this             ! fuel class
      real(r8),            intent(in)    :: sav_fuel(nfsc)   ! surface area to volume ratio of all fuel types [/cm]
      real(r8),            intent(in)    :: drying_ratio     ! drying ratio
      class(fire_weather), intent(in)    :: fireWeatherClass ! fireWeatherClass
      
      real(r8) :: moisture(nfsc)               ! fuel moisture [m3/m3]
      real(r8) :: moisture_of_extinction(nfsc) ! fuel moisture of extinction [m3/m3]
      integer  :: i                            ! looping index
      
      if (this%total_loading > nearzero) then 
        ! calculate fuel moisture [m3/m3] for each fuel class depending on what
        ! fire weather class is in use
        select type (fireWeatherClass)
          class is (nesterov_index)
            call CalculateFuelMoistureNesterov(sav_fuel, drying_ratio,                     &
              fireWeatherClass%fire_weather_index, moisture)
          class default 
            write(fates_log(), *) 'Unknown fire weather class selected.'
            write(fates_log(), *) 'Choose a different fire weather class or upate this subroutine.'
            call endrun(msg=errMsg( __FILE__, __LINE__))
        end select
              
        ! calculate moisture of extinction and fuel effective moisture
        do i = 1, nfsc
          moisture_of_extinction(i) = MoistureOfExtinction(sav_fuel(i))
          this%effective_moisture(i) = moisture(i)/moisture_of_extinction(i)
        end do
        
        ! average fuel moisture across all fuel types except trunks [m3/m3]
        this%average_moisture = sum(this%frac_loading(1:nfsc_notrunks)*moisture(1:nfsc_notrunks))
        
        ! calculate average moisture of extinction across all fuel types except trunks
        this%MEF = sum(this%frac_loading(1:nfsc_notrunks)*moisture_of_extinction(1:nfsc_notrunks))
      else 
        this%effective_moisture(1:nfsc) = 0.0_r8
        this%average_moisture = 0.0_r8 
        this%MEF = 0.0_r8
      end if
       
    end subroutine UpdateFuelMoisture
    
    !-------------------------------------------------------------------------------------
    
    real(r8) function CalculateFractionBurnt(effective_moisture, min_moisture,           &
      mid_moisture, moisture_coeff_low, moisture_slope_low, moisture_coeff_mid,          &
      moisture_slope_mid)
      ! DESCRIPTION:
      !   Calculates fraction burnt of fuel based on input fuel moisture
      
      ! Based on Equation B1 from Thonicke et al. 2010

      ! ARGUMENTS:
      real(r8), intent(in) :: effective_moisture ! effective fuel moisture [m3/m3]
      real(r8), intent(in) :: min_moisture       ! minimum moisture content for initial equation [m3/m3]
      real(r8), intent(in) :: mid_moisture       ! medium moisture content for initial equation [m3/m3]
      real(r8), intent(in) :: moisture_coeff_low ! coefficient for low moisture content [m3/m3]
      real(r8), intent(in) :: moisture_slope_low ! slope for low moisture content [m3/m3]
      real(r8), intent(in) :: moisture_coeff_mid ! coefficient for medium moisture content [m3/m3]
      real(r8), intent(in) :: moisture_slope_mid ! slope for low medium content [m3/m3] 
      
      if (effective_moisture <= min_moisture) then 
        CalculateFractionBurnt = 1.0_r8
      else if (effective_moisture > min_moisture .and. effective_moisture <= mid_moisture) then 
        CalculateFractionBurnt = max(0.0_r8, min(1.0_r8, moisture_coeff_low -            &
          moisture_slope_low*effective_moisture))
      else if (effective_moisture > mid_moisture .and. effective_moisture <= 1.0_r8) then 
        CalculateFractionBurnt = max(0.0_r8, min(1.0_r8, moisture_coeff_mid -            &
          moisture_slope_mid*effective_moisture))
      else 
        CalculateFractionBurnt = 0.0_r8
      end if
      
    end function CalculateFractionBurnt
    
    !-------------------------------------------------------------------------------------
    
    subroutine BurnFuel(this, fuel_consumed)
      ! DESCRIPTION:
      !   Calculates how much fuel burns
      
      ! USES
      use SFParamsMod, only : SF_val_miner_total, SF_val_min_moisture
      use SFParamsMod, only : SF_val_mid_moisture, SF_val_low_moisture_Coeff
      use SFParamsMod, only : SF_val_low_moisture_Slope, SF_val_mid_moisture_Coeff
      use SFParamsMod, only : SF_val_mid_moisture_Slope

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                         ! fuel class
      real(r8),         intent(out)   :: fuel_consumed(nfsc_notrunks) ! amount of fuel consumed in non-trunk litter classes [kgC/m2]
      
      ! LOCALS:
      integer :: i ! looping index for fuel classes
      
      do i = 1, nfsc
        this%frac_burnt = CalculateFractionBurnt(this%effective_moisture(i),             &
          SF_val_min_moisture(i), SF_val_mid_moisture(i), SF_val_low_moisture_Coeff(i),  &
          SF_val_low_moisture_Slope(i), SF_val_mid_moisture_Coeff(i),                    &
          SF_val_mid_moisture_Slope(i))
      end do 
          
      ! we can't ever kill all of the grass
      this%frac_burnt(fuel_classes%live_grass()) = min(0.8_r8,                           &
        this%frac_burnt(fuel_classes%live_grass()))  

      ! reduce burnt amount for mineral content
      this%frac_burnt(1:nfsc) = this%frac_burnt(1:nfsc)*(1.0_r8 - SF_val_miner_total) 
      
      ! calculate amount of fuel burned
      do i = 1, nfsc_notrunks
        fuel_consumed(i) = this%frac_burnt(i)*this%loading(i)
      end do
       
    end subroutine BurnFuel
    
    !-------------------------------------------------------------------------------------
    
    subroutine CalculateFuelMoistureNesterov(sav_fuel, drying_ratio, NI, moisture)
      !
      ! DESCRIPTION:
      !   Updates fuel moisture

      ! ARGUMENTS:
      real(r8), intent(in)  :: sav_fuel(nfsc) ! surface area to volume ratio of all fuel types [/cm]
      real(r8), intent(in)  :: drying_ratio   ! drying ratio
      real(r8), intent(in)  :: NI             ! Nesterov Index
      real(r8), intent(out) :: moisture(nfsc) ! moisture of litter [m3/m3]
      
      ! LOCALS
      integer  :: i         ! looping index
      real(r8) :: alpha_FMC ! intermediate variable for calculating fuel moisture
      
      do i = 1, nfsc
        if (i == fuel_classes%live_grass()) then 
          ! live grass moisture is a function of SAV and changes via Nesterov Index
          ! along the same relationship as the 1 hour fuels
          ! live grass has same SAV as dead grass, but retains more moisture with this calculation
          alpha_FMC = sav_fuel(fuel_classes%twigs())/drying_ratio
        else
          alpha_FMC = sav_fuel(i)/drying_ratio
        end if
        ! Equation 
        moisture(i) = exp(-1.0_r8*alpha_FMC*NI)
      end do
      
    end subroutine CalculateFuelMoistureNesterov
    
    !-------------------------------------------------------------------------------------
    
    real(r8) function MoistureOfExtinction(sav)
      !
      ! DESCRIPTION:
      !   Calculates moisture of extinction based on input surface area to volume ratio
    
      ! MEF (moisure of extinction) depends on compactness of fuel, depth, particle size, wind, slope
      ! Eqn here is eqn 27 from Peterson and Ryan (1986) "Modeling Postfire Conifer Mortality for Long-Range Planning"
      ! MEF: pine needles=0.30 (text near EQ 28 Rothermal 1972)
      ! Table II-1 NFFL mixed fuels models from Rothermal 1983 Gen. Tech. Rep. INT-143 
      ! MEF: short grass=0.12,tall grass=0.25,chaparral=0.20,closed timber litter=0.30,hardwood litter=0.25
      ! Thonicke 2010 SAV values propagated thru P&R86 eqn below gives MEF:tw=0.355, sb=0.44, lb=0.525, tr=0.63, dg=0.248, lg=0.248
      ! Lasslop 2014 Table 1 MEF PFT level:grass=0.2,shrubs=0.3,TropEverGrnTree=0.2,TropDecid Tree=0.3, Extra-trop Tree=0.3

      ! ARGUMENTS:
      real(r8), intent(in) :: sav ! fuel surface area to volume ratio [/cm]
      
      ! CONSTANTS:
      real(r8), parameter :: MEF_a = 0.524_r8
      real(r8), parameter :: MEF_b = 0.066_r8
      
      if (sav <= nearzero) then
        MoistureOfExtinction = 0.0_r8
      else
        MoistureOfExtinction = MEF_a - MEF_b*log(sav)
      end if
    
    end function MoistureOfExtinction
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageBulkDensity(this, bulk_density)
      ! DESCRIPTION:
      !   Calculates average bulk density (not including trunks)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this               ! fuel class
      real(r8),           intent(in)    :: bulk_density(nfsc) ! bulk density of all fuel types [kg/m2]
      
      if (this%total_loading > nearzero) then 
        this%bulk_density = sum(this%frac_loading(1:nfsc_notrunks)*bulk_density(1:nfsc_notrunks))
      else 
        this%bulk_density = 0.0_r8
      end if
    
    end subroutine AverageBulkDensity
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageSAV(this, sav_fuel)
      ! DESCRIPTION:
      !   Calculates average surface area to volume ratio (not including trunks)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this           ! fuel class
      real(r8),           intent(in)    :: sav_fuel(nfsc) ! surface area to volume ratio of all fuel types [/cm]
      
      if (this%total_loading > nearzero) then 
        this%SAV = sum(this%frac_loading(1:nfsc_notrunks)*sav_fuel(1:nfsc_notrunks))
      else 
        this%SAV = 0.0_r8
      end if
    
    end subroutine AverageSAV
    
    !-------------------------------------------------------------------------------------
    
end module FatesFuelMod