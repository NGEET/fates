module FatesFuelMod

  use FatesFuelClassesMod, only : nfsc, fuel_classes
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
    
    real(r8) :: loading(nfsc)             ! fuel loading of non-trunks fuel class [kgC/m2]
    real(r8) :: trunk_loading             ! fuel loading of trunk fuel class [kgC/m2]
    real(r8) :: effective_moisture(nfsc)  ! fuel effective moisture all fuel class (moisture/MEF) [m3/m3]
    real(r8) :: frac_loading(nfsc)        ! fractional loading of non-trunk fuel classes [0-1] 
    real(r8) :: frac_burnt(nfsc)          ! fraction of litter burnt by fire [0-1]
    real(r8) :: total_loading             ! total fuel loading - DOES NOT INCLUDE TRUNKS [kgC/m2]
    real(r8) :: average_moisture          ! weighted average of fuel moisture across non-trunk fuel classes [m3/m3]
    real(r8) :: bulk_density              ! weighted average of bulk density across non-trunk fuel classes [kg/m3]
    real(r8) :: SAV                       ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
    real(r8) :: MEF                       ! weighted average of moisture of extinction across non-trunk fuel classes [m3/m3]

    contains
      
      procedure :: Init
      procedure :: CalculateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity
      procedure :: AverageSAV

  end type fuel_type
  
  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! just zero everything
      this%loading(1:nfsc) = 0.0_r8
      this%trunk_loading = 0.0_r8
      this%frac_loading(1:nfsc) = 0.0_r8
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
      this%loading(fuel_classes%trunks()) = trunk_litter
      this%trunk_loading = trunk_litter

    end subroutine CalculateLoading

    !-------------------------------------------------------------------------------------

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sums up the loading

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! LOCALS:
      integer :: i ! looping index
      
      this%total_loading = 0.0_r8
      do i = 1, nfsc
        if (i /= fuel_classes%trunks()) then 
          this%total_loading = this%total_loading + this%loading(i)
        end if
      end do

    end subroutine SumLoading

    !-------------------------------------------------------------------------------------
    
    subroutine CalculateFractionalLoading(this)
      ! DESCRIPTION:
      !   Calculates fractional loading for fuel

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! LOCALS:
      integer :: i ! looping index

      ! sum up loading just in case
      call this%SumLoading()
      
      if (this%total_loading > nearzero) then
        do i = 1, nfsc 
          if (i /= fuel_classes%trunks()) then 
            this%frac_loading(i) = this%loading(i)/this%total_loading
          end if 
        end do 
      else 
        this%frac_loading(1:nfsc) = 0.0_r8
      end if 

    end subroutine CalculateFractionalLoading

    !-------------------------------------------------------------------------------------
    
    subroutine UpdateFuelMoisture(this, sav_fuel, drying_ratio, fireWeatherClass)
      ! DESCRIPTION:
      !   Updates fuel moisture depending on what fire weather class is in use
      
      use SFParamsMod, only : SF_val_SAV, SF_val_drying_ratio
      
      ! ARGUMENTS:
      class(fuel_type),    intent(inout) :: this             ! fuel class
      real(r8),            intent(in)    :: sav_fuel(nfsc)   ! surface area to volume ratio of all fuel types [/cm]
      real(r8),            intent(in)    :: drying_ratio     ! drying ratio
      class(fire_weather), intent(in)    :: fireWeatherClass ! fireWeatherClass
      
      real(r8) :: moisture(nfsc)               ! fuel moisture [m3/m3]
      real(r8) :: moisture_of_extinction(nfsc) ! fuel moisture of extinction [m3/m3]
      integer  :: i                            ! looping index
      integer  :: tw_sf, dl_sf, lg_sf, lb_sf, tr_sf
      
      tw_sf = fuel_classes%twigs()
      dl_sf = fuel_classes%dead_leaves()
      lg_sf = fuel_classes%live_grass()
      lb_sf = fuel_classes%large_branches()
      tr_sf = fuel_classes%trunks()
      
      if (this%total_loading > nearzero) then 
        ! calculate fuel moisture [m3/m3] for each fuel class depending on what
        ! fire weather class is in use
        select type (fireWeatherClass)
          class is (nesterov_index)
            call CalculateFuelMoistureNesterov(sav_fuel, drying_ratio,                   &
              fireWeatherClass%fire_weather_index, moisture)
          class default 
            write(fates_log(), *) 'Unknown fire weather class selected.'
            write(fates_log(), *) 'Choose a different fire weather class or upate this subroutine.'
            call endrun(msg=errMsg( __FILE__, __LINE__))
        end select
        
        this%average_moisture = 0.0_r8
        this%MEF = 0.0_r8
        do i = 1, nfsc
          ! calculate moisture of extinction and fuel effective moisture
          moisture_of_extinction(i) = MoistureOfExtinction(sav_fuel(i))
          this%effective_moisture(i) = moisture(i)/moisture_of_extinction(i)
          
          ! average fuel moisture  and MEF across all fuel types except trunks [m3/m3]
          if (i /= fuel_classes%trunks()) then 
            this%average_moisture = this%average_moisture + this%frac_loading(i)*moisture(i)
            this%MEF = this%MEF + this%frac_loading(i)*moisture_of_extinction(i)
          end if 
        end do
        
  
        ! this%MEF = sum(this%frac_loading(tw_sf:lb_sf)*moisture_of_extinction(tw_sf:lb_sf))              
        ! this%average_moisture = sum(this%frac_loading(tw_sf:lb_sf)*moisture(tw_sf:lb_sf))
        
        ! this%MEF = this%MEF + sum(this%frac_loading(dl_sf:lg_sf)*moisture_of_extinction(dl_sf:lg_sf))            
        ! this%average_moisture = this%average_moisture + sum(this%frac_loading(dl_sf:lg_sf)*moisture(dl_sf:lg_sf))
        
        ! ! Correct averaging for the fact that we are not using the trunks pool for fire ROS and intensity (5)
        ! ! Consumption of fuel in trunk pool does not influence fire ROS or intensity (Pyne 1996)
        ! if ((1.0_r8 - this%frac_loading(tr_sf)) > nearzero) then
        !   this%MEF = this%MEF*(1.0_r8/(1.0_r8 - this%frac_loading(tr_sf)))
        !   this%average_moisture = this%average_moisture*(1.0_r8/(1.0_r8 - this%frac_loading(tr_sf)))
        ! else
        !   ! somehow the fuel is all trunk. put dummy values from large branches so as not to break things later in code.
        !   this%MEF = moisture_of_extinction(lb_sf)
        !   this%average_moisture = moisture(lb_sf)
        ! endif
        
      else 
        this%effective_moisture(1:nfsc) = 0.0_r8
        this%average_moisture = 0.0000000001_r8 
        this%MEF = 0.0000000001_r8
      end if
       
    end subroutine UpdateFuelMoisture
    
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
    
      ! MEF (moisure of extinction) depends on compactness of fuel, depth, particle size, 
      !  wind, and slope
      ! Equation here is Eq. 27 from Peterson and Ryan (1986) "Modeling Postfire Conifer 
      ! Mortality for Long-Range Planning"
      !
      ! Example MEFs:
      ! pine needles = 0.30 (Rothermal 1972)
      ! short grass = 0.12 (Rothermel 1983; Gen. Tech. Rep. INT-143; Table II-1)
      ! tall grass = 0.24 (Rothermel 1983)
      ! chaparral = 0.20 (Rothermel 1983)
      ! closed timber litter = 0.30 (Rothermel 1983)
      ! hardwood litter = 0.25 (Rothermel 1983)
      ! grass = 0.2 (Lasslop 2014; Table 1)
      ! shrubs = 0.3 (Lasslop 2014; Table 1)
      ! tropical evergreen trees = 0.2 (Lasslop 2014; Table 1)
      ! tropical deciduous trees = 0.3 (Lasslop 2014; Table 1)
      ! extratropical trees = 0.3 (Lasslop 2014; Table 1)
      !
      ! SAV values from Thonicke 2010 give: 
      ! twigs = 0.355, small branches = 0.44, large branches = 0.525, trunks = 0.63
      ! dead leaves = 0.248, live grass = 0.248
      !
    
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
        this%bulk_density = sum(this%frac_loading(1:nfsc)*bulk_density(1:nfsc))
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
        this%SAV = sum(this%frac_loading(1:nfsc)*sav_fuel(1:nfsc))
      else 
        this%SAV = 0.0_r8
      end if
    
    end subroutine AverageSAV
    
    !-------------------------------------------------------------------------------------
    
end module FatesFuelMod