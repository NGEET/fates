module FatesFuelMod

  use FatesFuelClassesMod, only : num_fuel_classes, fuel_classes
  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : nearzero
  use SFNesterovMod,       only : nesterov_index
  use SFFireWeatherMod,    only : fire_weather
  use FatesLitterMod,      only : ncwd
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
    real(r8) :: canopy_fuel_load                     ! patch level total canopy fuel load [kg biomass]
    real(r8) :: canopy_base_height                   ! patch level canopy base height at which biomass density > minimum density required for canopy fire spread [m]
    real(r8) :: canopy_bulk_density                  ! patch level canopy fuel bulk density [kg biomass m-3]
    real(r8) :: canopy_water_content                 ! patch level canopy water content [%]


    contains
      
      procedure :: Init
      procedure :: Fuse
      procedure :: UpdateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity_NoTrunks
      procedure :: AverageSAV_NoTrunks
      procedure :: CalculateCanopyFuelLoad
      procedure :: CalculateCanopyBulkDensity
      procedure :: CalculateFuelBurnt
      procedure :: CalculateResidenceTime
      procedure :: CanopyWaterContent

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
      this%canopy_fuel_load = 0.0_r8
      this%canopy_base_height = 0.0_r8
      this%canopy_bulk_density = 0.0_r8
      this%canopy_water_content = 0.0_r8

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
      this%canopy_fuel_load = this%canopy_fuel_load*self_weight +      &
        donor_fuel%canopy_fuel_load*donor_weight
      this%canopy_base_height = this%canopy_base_height*self_weight +   &
        donor_fuel%canopy_base_height*donor_weight
      this%canopy_bulk_density = this%canopy_bulk_density*self_weight +    &
        donor_fuel%canopy_bulk_density*donor_weight
      this%canopy_water_content = this%canopy_water_content*self_weight +    &
        donor_fuel%canopy_water_content*donor_weight
      
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
      class(fuel_type), intent(inout) :: this                       ! fuel class
      real(r8),         intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      
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

    subroutine CalculateCanopyFuelLoad(this, leaf_c, woody_c, &
      cwd_frac_adj, canopy_fuel_1h)
      ! DESCRIPTION:
      ! Calculate canopy fuel load by summing up leaf biomass and 1 hour woody biomass
      ! XLG: I did not change how Sam L. calculate canopy fuel load here except 
      ! switch to the adjusted CWD instead of using the default parameter
      ! Seems only including leaf biomass is recommend by Gruz et al. 2003:
      ! https://www.sierraforestlegacy.org/Resources/Conservation/FireForestEcology/FireScienceResearch/FireEcology/FireEcology-Cruz03.pdf 

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                    ! fuel class
      real(r8),         intent(in)    :: leaf_c                  ! leaf biomass of each cohort [kg biomass]
      real(r8),         intent(in)    :: woody_c                 ! 1-hour woody fuels [kg biomass]
      real(r8),         intent(in)    :: cwd_frac_adj(ncwd)      ! adjusted coarse woody debris fraction [fraction]
      real(r8),         intent(out)   :: canopy_fuel_1h          ! leaf + 1-hour woody biomass of current cohort [kg biomass]

      ! calculate canopy fuel load for FATES crown fire
      canopy_fuel_1h = woody_c * cwd_frac_adj(1) + leaf_c
      this%canopy_fuel_load = this%canopy_fuel_load + canopy_fuel_1h

      
    end subroutine CalculateCanopyFuelLoad


  !---------------------------------------------------------------------------------------
    subroutine CalculateCanopyBulkDensity(this, biom_matrix, max_height)
      ! DESCRIPTION:
      ! Calculate canopy fuel bulk density [kg biomass / m3] by searching for 
      ! the canopy height (canopy base height CBH) at which canopy fuel bed density 
      ! is > minimum canopy fuel density [0.011 kg biomass/m3] for propogating canopy fire; 
      ! CBD is then calculated as canopy fuel load devided by 
      ! crown depth ( CD= max. canopy height - CBH).
      ! XLG: one modification I made compare to Sam L. version is that I recalculate
      ! canopy fuel load after knowing CD by excluding fuels below CBH. As if including
      ! fuels below CBH then we are artificially "condense" canopy fuel bed by deviding total
      ! fuel load by an average, shallower CD. 

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                  ! fuel class
      real(r8),         intent(in)    :: biom_matrix(0:)        ! 1m biomass bin (kg/m3) in the vertical space to sort the canopy fule 
                                                               ! across cohorts, used to search for CBH
      real(r8),         intent(in)    :: max_height            ! the max. cohort height at current patch [m]

      ! Locals:
      real(r8) :: canopy_top_height             ! the highest point at which biomass density is > minimum canopy fuel density [m]
      integer  :: ih                            ! biomass bin counter

      ! Params:
      real(r8), parameter :: min_density_canopy_fuel = 0.011_r8   ! minimum canopy fuel density for propogating canopy fire
                                                                  ! vertically through canopy.
                                                                  ! Scott and Reinhardt 2001 RMRS-RP-29

      ! loop from 1m to 70m to find CBH
      do ih=0,int(max_height)
        if (biom_matrix(ih) > min_density_canopy_fuel) then
          this%canopy_base_height = dble(ih) + 1.0_r8   ! the dimension index of biom_matrix is a ronded-down integer of cohort height
                                                   ! add 1 to be conservative when searching for CBH

          exit
        else                                            ! when this is an open stand that there is no such a height, use max_height
          this%canopy_base_height = max_height
        end if
      end do

      ! now loop from top to bottom to find the highest point at which the minimum bulk density is met
      ! XLG: I modified the way how canopy bulk density is calculated by reducing the maximum canopy height
      ! to the highest point where the minimum bulk density is met. 
  
      do ih = int(max_height), 0, -1
        if (biom_matrix(ih) > min_density_canopy_fuel) then
          canopy_top_height = dble(ih) + 1.0_r8 
          exit
        else
          canopy_top_height = max_height 
        end if
      end do

      ! XLG: We now only calculate canopy bulk density for fuels between
      ! canopy base and top height 
      if ((canopy_top_height - this%canopy_base_height) > nearzero) then
        this%canopy_bulk_density = sum(biom_matrix(int(this%canopy_base_height-1.0_r8):int(canopy_top_height-1.0_r8))) / &
        (canopy_top_height - this%canopy_base_height)
      else
        this%canopy_bulk_density = 0.0_r8
      end if
      

    end subroutine CalculateCanopyBulkDensity

 !---------------------------------------------------------------------------------------

    subroutine CanopyWaterContent(this, co_cwc, co_fuel)
      ! DESCRIPTION:
      ! Calculates live canopy water content  
      ! as sum of fuel load weighted cohort level live fuel moisture content
      !
      ! ARGUMENTS
      class(fuel_type), intent(inout) :: this                  ! fuel class
      real(r8),         intent(in)    :: co_cwc               ! coohort live fuel moisture content [%]
      real(r8),         intent(in)    :: co_fuel               ! cohort canopy fuel, only 1 hour woody + leaf 

      this%canopy_water_content = this%canopy_water_content + &
      co_cwc * co_fuel/this%canopy_fuel_load
      

    end subroutine CanopyWaterContent

 !---------------------------------------------------------------------------------------

   
    subroutine CalculateFuelBurnt(this, fuel_consumed)
      ! DESCRIPTION:
      !   Calculates the fraction and total amount of fuel burnt
      !
      
      use SFParamsMod, only : SF_val_mid_moisture, SF_val_mid_moisture_Coeff
      use SFParamsMod, only : SF_val_mid_moisture_Slope, SF_val_min_moisture
      use SFParamsMod, only : SF_val_low_moisture_Coeff, SF_val_low_moisture_Slope
      use SFParamsMod, only : SF_val_miner_total

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                            ! fuel class
      real(r8),         intent(out)   :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
      
      ! LOCALS:
      real(r8) :: rel_moisture  ! relative moisture of fuel (moist/moisture of extinction) [unitless]
      integer  :: i             ! looping index
      
      ! CONSTANTS:
      real(r8), parameter :: max_grass_frac = 0.8_r8 ! maximum fraction burnt for live grass fuels
      
      this%frac_burnt(:) = 1.0_r8
        
      ! Calculate fraction of litter is burnt for all classes. 
      ! Equation B1 in Thonicke et al. 2010
      do i = 1, num_fuel_classes        
        
        rel_moisture = this%effective_moisture(i)                  
        
        if (rel_moisture <= SF_val_min_moisture(i)) then
          ! very dry litter
          this%frac_burnt(i) = 1.0_r8 
        else if (rel_moisture > SF_val_min_moisture(i) .and. rel_moisture <= SF_val_mid_moisture(i)) then
          ! low to medium moisture
          this%frac_burnt(i) = max(0.0_r8, min(1.0_r8, SF_val_low_moisture_Coeff(i) - &
            SF_val_low_moisture_Slope(i)*rel_moisture))
        else if (rel_moisture > SF_val_mid_moisture(i) .and. rel_moisture <= 1.0_r8) then
          ! medium to high moisture
          this%frac_burnt(i) = max(0.0_r8, min(1.0_r8, SF_val_mid_moisture_Coeff(i) - &
            SF_val_mid_moisture_Slope(i)*rel_moisture))
        else 
          ! very wet litter
          this%frac_burnt(i) = 0.0_r8  
        endif
        
        ! we can't ever kill all of the grass
        if (i == fuel_classes%live_grass()) then
          this%frac_burnt(i) = min(max_grass_frac, this%frac_burnt(i))
        end if
        
        ! reduce fraction burnt based on mineral content
        this%frac_burnt(i) = this%frac_burnt(i)*(1.0_r8 - SF_val_miner_total)
        
        ! calculate fuel consumed
        fuel_consumed(i) = this%frac_burnt(i)*this%loading(i)
      end do

    end subroutine CalculateFuelBurnt
    
    !-------------------------------------------------------------------------------------
    
    subroutine CalculateResidenceTime(this, tau_l)
      !
      !  DESCRIPTION:
      !  Calculates fire residence time, duration of lethal bole heating [min]
      !  This is used for determining cambial kill of woody cohorts
      !
      !  From Peterson & Ryan (1986)
      !
      
      ! ARGUMENTS:
      class(fuel_type), intent(in)  :: this  ! fuel class
      real(r8),         intent(out) :: tau_l ! duration of lethal bole heating [min]
      
      ! LOCALS:
      integer :: i ! looping index
      
      tau_l = 0.0_r8 
      do i = 1, num_fuel_classes
        if (i /= fuel_classes%trunks()) then 
          ! don't include 1000-hr fuels
          ! convert loading from kgC/m2 to g/cm2
          tau_l = tau_l + 39.4_r8*(this%frac_loading(i)*this%non_trunk_loading/0.45_r8/10.0_r8)* &
            (1.0_r8 - ((1.0_r8 - this%frac_burnt(i))**0.5_r8))
        end if 
      end do
      
      ! cap the residence time to 8mins, as suggested by literature survey by P&R (1986)
      tau_l = min(8.0_r8, tau_l) 

    end subroutine CalculateResidenceTime
    
end module FatesFuelMod
