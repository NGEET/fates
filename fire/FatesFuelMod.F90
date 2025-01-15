module FatesFuelMod

  use FatesFuelClassesMod, only : num_fuel_classes, fuel_classes
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
    
    real(r8) :: loading(num_fuel_classes)            ! fuel loading of each fuel class [kgC/m2]
    real(r8) :: effective_moisture(num_fuel_classes) ! fuel effective moisture all fuel class (moisture/MEF) [m3/m3]
    real(r8) :: frac_loading(num_fuel_classes)       ! fractional loading of all fuel classes [0-1] 
    real(r8) :: frac_burnt(num_fuel_classes)         ! fraction of litter burnt by fire [0-1]
    real(r8) :: non_trunk_loading                    ! total fuel loading excluding trunks [kgC/m2]
    real(r8) :: average_moisture_notrunks            ! weighted average of fuel moisture across non-trunk fuel classes [m3/m3]
    real(r8) :: bulk_density_notrunks                ! weighted average of bulk density across non-trunk fuel classes [kg/m3]
    real(r8) :: SAV_notrunks                         ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
    real(r8) :: MEF_notrunks                         ! weighted average of moisture of extinction across non-trunk fuel classes [m3/m3]

    contains
      
      procedure :: Init
      procedure :: Fuse
      procedure :: UpdateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity_NoTrunks
      procedure :: AverageSAV_NoTrunks

  end type fuel_type
  
  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! just zero everything
      this%loading(1:num_fuel_classes) = 0.0_r8
      this%frac_loading(1:num_fuel_classes) = 0.0_r8
      this%frac_burnt(1:num_fuel_classes) = 0.0_r8  
      this%effective_moisture(1:num_fuel_classes) = 0.0_r8
      this%non_trunk_loading = 0.0_r8
      this%average_moisture_notrunks = 0.0_r8 
      this%bulk_density_notrunks = 0.0_r8
      this%SAV_notrunks = 0.0_r8
      this%MEF_notrunks = 0.0_r8 

    end subroutine Init 
    
    !-------------------------------------------------------------------------------------
    
    subroutine Fuse(this, self_area, donor_area, donor_fuel)
      ! DESCRIPTION:
      !   Fuse attributes of this object with another

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this       ! fuel class
      real(r8),         intent(in)    :: self_area  ! area of this fuel class's patch [m2]
      real(r8),         intent(in)    :: donor_area ! area of donor fuel class's patch [m2]
      type(fuel_type),  intent(in)    :: donor_fuel ! donor fuel class
      
      ! LOCALS:
      integer  :: i            ! looping index
      real(r8) :: self_weight  ! weighting of the receiving fuel class
      real(r8) :: donor_weight ! weighting of the donor fuel class
      
      self_weight = self_area/(donor_area + self_area)
      donor_weight = 1.0_r8 - self_weight
      
      do i = 1, num_fuel_classes 
        this%loading(i) = this%loading(i)*self_weight +                                  &
          donor_fuel%loading(i)*donor_weight
        this%frac_loading(i) = this%frac_loading(i)*self_weight +                        &
          donor_fuel%frac_loading(i)*donor_weight
        this%frac_burnt(i) = this%frac_burnt(i)*self_weight +                            &
          donor_fuel%frac_burnt(i)*donor_weight
        this%effective_moisture(i) = this%effective_moisture(i)*self_weight +            &
          donor_fuel%effective_moisture(i)*donor_weight
      end do 
      
      this%non_trunk_loading = this%non_trunk_loading*self_weight +                      &
        donor_fuel%non_trunk_loading*donor_weight
      this%average_moisture_notrunks = this%average_moisture_notrunks*self_weight +                        &
        donor_fuel%average_moisture_notrunks*donor_weight
      this%bulk_density_notrunks = this%bulk_density_notrunks*self_weight +              &
        donor_fuel%bulk_density_notrunks*donor_weight      
      this%SAV_notrunks = this%SAV_notrunks*self_weight + donor_fuel%SAV_notrunks*donor_weight
      this%MEF_notrunks = this%MEF_notrunks*self_weight + donor_fuel%MEF_notrunks*donor_weight    
      
    end subroutine Fuse

    !-------------------------------------------------------------------------------------

    subroutine UpdateLoading(this, leaf_litter, twig_litter, small_branch_litter,     &
        large_branch_litter, trunk_litter, live_grass)
      ! DESCRIPTION:
      !   Updates loading for each fuel type

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

    end subroutine UpdateLoading

    !-------------------------------------------------------------------------------------

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sums up the loading - excludes trunks
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! LOCALS:
      integer :: i ! looping index
      
      this%non_trunk_loading = 0.0_r8
      do i = 1, num_fuel_classes
        if (i /= fuel_classes%trunks()) then 
          this%non_trunk_loading = this%non_trunk_loading + this%loading(i)
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
      
      if (this%non_trunk_loading > nearzero) then
        do i = 1, num_fuel_classes 
          if (i /= fuel_classes%trunks()) then 
            this%frac_loading(i) = this%loading(i)/this%non_trunk_loading
          else 
            this%frac_loading(i) = 0.0_r8
          end if 
        end do 
      else 
        this%frac_loading(1:num_fuel_classes) = 0.0_r8
        this%non_trunk_loading = 0.0_r8
      end if 

    end subroutine CalculateFractionalLoading

    !-------------------------------------------------------------------------------------
    
    subroutine UpdateFuelMoisture(this, sav_fuel, drying_ratio, fireWeatherClass)
      ! DESCRIPTION:
      !   Updates fuel moisture depending on what fire weather class is in use
      
      ! ARGUMENTS:
      class(fuel_type),    intent(inout) :: this                       ! fuel class
      real(r8),            intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8),            intent(in)    :: drying_ratio               ! drying ratio
      class(fire_weather), intent(in)    :: fireWeatherClass           ! fireWeatherClass
      
      real(r8) :: moisture(num_fuel_classes)               ! fuel moisture [m3/m3]
      real(r8) :: moisture_of_extinction(num_fuel_classes) ! fuel moisture of extinction [m3/m3]
      integer  :: i                                        ! looping index
 
      if (this%non_trunk_loading + this%loading(fuel_classes%trunks()) > nearzero) then 
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
        
        this%average_moisture_notrunks = 0.0_r8
        this%MEF_notrunks = 0.0_r8
        do i = 1, num_fuel_classes
          ! calculate moisture of extinction and fuel effective moisture
          moisture_of_extinction(i) = MoistureOfExtinction(sav_fuel(i))
          this%effective_moisture(i) = moisture(i)/moisture_of_extinction(i)
          
          ! average fuel moisture  and MEF across all fuel types except trunks [m3/m3]
          if (i /= fuel_classes%trunks()) then 
            this%average_moisture_notrunks = this%average_moisture_notrunks + this%frac_loading(i)*moisture(i)
            this%MEF_notrunks = this%MEF_notrunks + this%frac_loading(i)*moisture_of_extinction(i)
          end if 
        end do

      else 
        this%effective_moisture(1:num_fuel_classes) = 0.0_r8
        this%average_moisture_notrunks = 0.0_r8
        this%MEF_notrunks = 0.0_r8
      end if
      
    end subroutine UpdateFuelMoisture
    
    !-------------------------------------------------------------------------------------
        
    subroutine CalculateFuelMoistureNesterov(sav_fuel, drying_ratio, NI, moisture)
      !
      ! DESCRIPTION:
      !   Updates fuel moisture

      ! ARGUMENTS:
      real(r8), intent(in)  :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8), intent(in)  :: drying_ratio               ! drying ratio
      real(r8), intent(in)  :: NI                         ! Nesterov Index
      real(r8), intent(out) :: moisture(num_fuel_classes) ! moisture of litter [m3/m3]
      
      ! LOCALS
      integer  :: i         ! looping index
      real(r8) :: alpha_FMC ! intermediate variable for calculating fuel moisture
      
      do i = 1, num_fuel_classes
        if (i == fuel_classes%live_grass()) then 
          ! live grass moisture is a function of SAV and changes via Nesterov Index
          ! along the same relationship as the 1 hour fuels
          ! live grass has same SAV as dead grass, but retains more moisture with this calculation
          alpha_FMC = sav_fuel(fuel_classes%twigs())/drying_ratio
        else
          alpha_FMC = sav_fuel(i)/drying_ratio
        end if
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
      ! pine needles = 0.30 (Rothermel 1972)
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
      
      if (sav <= 0.0_r8) then
        write(fates_log(), *) 'SAV cannot be negative - SAV'
        call endrun(msg=errMsg(__FILE__, __LINE__))
      else
        MoistureOfExtinction = MEF_a - MEF_b*log(sav)
      end if
    
    end function MoistureOfExtinction
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageBulkDensity_NoTrunks(this, bulk_density)
      ! DESCRIPTION:
      !   Calculates average bulk density (not including trunks)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this                           ! fuel class
      real(r8),           intent(in)    :: bulk_density(num_fuel_classes) ! bulk density of all fuel types [kg/m3]
      
      ! LOCALS:
      integer :: i ! looping index
      
      if (this%non_trunk_loading > nearzero) then
        this%bulk_density_notrunks = 0.0_r8
        do i = 1, num_fuel_classes               
          ! average bulk density across all fuel types except trunks 
          if (i /= fuel_classes%trunks()) then 
            this%bulk_density_notrunks = this%bulk_density_notrunks + this%frac_loading(i)*bulk_density(i)
          end if 
        end do
      else 
        this%bulk_density_notrunks = sum(bulk_density(1:num_fuel_classes))/num_fuel_classes
      end if
      
    end subroutine AverageBulkDensity_NoTrunks
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageSAV_NoTrunks(this, sav_fuel)
      ! DESCRIPTION:
      !   Calculates average surface area to volume ratio (not including trunks)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this                       ! fuel class
      real(r8),           intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      
      ! LOCALS:
      integer :: i ! looping index
      
      if (this%non_trunk_loading > nearzero) then
        this%SAV_notrunks = 0.0_r8
        do i = 1, num_fuel_classes               
          ! average bulk density across all fuel types except trunks 
          if (i /= fuel_classes%trunks()) then 
            this%SAV_notrunks = this%SAV_notrunks + this%frac_loading(i)*sav_fuel(i)
          end if 
        end do
      else 
        this%SAV_notrunks = sum(sav_fuel(1:num_fuel_classes))/num_fuel_classes 
      end if
    
    end subroutine AverageSAV_NoTrunks
    
  !---------------------------------------------------------------------------------------
    
end module FatesFuelMod
