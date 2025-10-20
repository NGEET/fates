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
      real(r8) :: mean_total_SA(num_fuel_classes)      ! mean total surface area per unit fuel cell of each fuel class [m2]
      real(r8) :: weighting_factor(num_fuel_classes)   ! weighting factor for calculating average fuel characteristics [unitless]
      real(r8) :: fuel_weight(num_fuel_classes)        ! weighting factor for calculating average fuel load [unitless]
      real(r8) :: moisture(num_fuel_classes)           ! fuel moisture by fuel class [m3 m-3]
      real(r8) :: effective_moisture(num_fuel_classes) ! fuel effective moisture all fuel class (moisture/MEF) [m3/m3]
      real(r8) :: frac_loading(num_fuel_classes)       ! fractional loading of all fuel classes [0-1]
      real(r8) :: frac_burnt(num_fuel_classes)         ! fraction of litter burnt by fire [0-1]
      real(r8) :: non_trunk_loading                    ! total fuel loading excluding trunks [kgC/m2]
      real(r8) :: weighted_loading_dead                ! this is dead fuel load weighted by weighting_factor given SA:V [kgC/m2]
      real(r8) :: weighted_loading_live                ! this is live fuel load weighted by weighting factor given SA:V [kgC/m2]
      real(r8) :: wf_dead                              ! weighting factor for dead fuel category [unitless]
      real(r8) :: wf_live                              ! weighting factor for live fuel category [unitless]
      real(r8) :: average_moisture_dead                ! weighted average of fuel moisture across all fuel classes for dead fuels [m3/m3]
      real(r8) :: average_moisture_live                ! weighted average of fuel moisture across all fuel classes for live fuels [m3/m3]
      real(r8) :: bulk_density_weighted                ! weighted average of bulk density across all fuel classes [kg/m3]
      real(r8) :: SAV_weighted                         ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
      real(r8) :: MEF_dead                             ! weighted average of moisture of extinction across all dead fuel classes [m3/m3]
      real(r8) :: MEF_live                             ! moisture of extinction for live fuels [m3/m3]

   contains

      procedure :: Init
      procedure :: Fuse
      procedure :: UpdateLoading
      procedure :: CalculateWeightingFactor
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity
      procedure :: AverageSAV
      procedure :: CalculateFuelBurnt
      procedure :: CalculateResidenceTime
      procedure :: FuelLoadWeight

   end type fuel_type

contains

   subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      ! just zero everything
      this%loading(1:num_fuel_classes) = 0.0_r8
      this%mean_total_SA(1:num_fuel_classes) = 0.0_r8
      this%weighting_factor(1:num_fuel_classes) = 0.0_r8
      this%fuel_weight(1:num_fuel_classes) = 0.0_r8
      this%frac_loading(1:num_fuel_classes) = 0.0_r8
      this%frac_burnt(1:num_fuel_classes) = 0.0_r8
      this%moisture(1:num_fuel_classes) = 0.0_r8
      this%effective_moisture(1:num_fuel_classes) = 0.0_r8
      this%non_trunk_loading = 0.0_r8
      this%weighted_loading_dead = 0.0_r8
      this%weighted_loading_live = 0.0_r8
      this%wf_dead = 0.0_r8
      this%wf_live = 0.0_r8
      this%average_moisture_dead = 0.0_r8
      this%average_moisture_live = 0.0_r8
      this%bulk_density_weighted = 0.0_r8
      this%SAV_weighted = 0.0_r8
      this%MEF_dead = 0.0_r8
      this%MEF_live = 0.0_r8

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
         this%mean_total_SA(i) = this%mean_total_SA(i)*self_weight +                      &
            donor_fuel%mean_total_SA(i)*donor_weight
         this%weighting_factor(i) = this%weighting_factor(i)*self_weight +              &
            donor_fuel%weighting_factor(i)*donor_weight
         this%fuel_weight(i) = this%fuel_weight(i)*self_weight +              &
            donor_fuel%fuel_weight(i)*donor_weight
         this%frac_loading(i) = this%frac_loading(i)*self_weight +                        &
            donor_fuel%frac_loading(i)*donor_weight
         this%frac_burnt(i) = this%frac_burnt(i)*self_weight +                            &
            donor_fuel%frac_burnt(i)*donor_weight
         this%moisture(i) = this%moisture(i)*self_weight +            &
            donor_fuel%moisture(i)*donor_weight
         this%effective_moisture(i) = this%effective_moisture(i)*self_weight +            &
            donor_fuel%effective_moisture(i)*donor_weight
      end do

      this%non_trunk_loading = this%non_trunk_loading*self_weight +                      &
         donor_fuel%non_trunk_loading*donor_weight
      this%weighted_loading_dead = this%weighted_loading_dead*self_weight +            &
         donor_fuel%weighted_loading_dead*donor_weight
      this%weighted_loading_live = this%weighted_loading_live*self_weight +            &
         donor_fuel%weighted_loading_live*donor_weight
      this%wf_dead = this%wf_dead*self_weight + donor_fuel%wf_dead*donor_weight
      this%wf_live = this%wf_live*self_weight + donor_fuel%wf_live*donor_weight
      this%average_moisture_dead = this%average_moisture_dead*self_weight +              &
         donor_fuel%average_moisture_dead*donor_weight
      this%average_moisture_live = this%average_moisture_live*self_weight +              &
         donor_fuel%average_moisture_live*donor_weight
      this%bulk_density_weighted = this%bulk_density_weighted*self_weight +              &
         donor_fuel%bulk_density_weighted*donor_weight
      this%SAV_weighted = this%SAV_weighted*self_weight + donor_fuel%SAV_weighted*donor_weight
      this%MEF_dead = this%MEF_dead*self_weight + donor_fuel%MEF_dead*donor_weight
      this%MEF_live = this%MEF_live*self_weight + donor_fuel%MEF_live*donor_weight

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

   subroutine CalculateWeightingFactor(this, sav_fuel, part_dens)
      ! DESCRIPTION:
      ! Calculate total area per fuel cell of each fuel class
      ! and the derived weighting factors

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                       ! fuel class
      real(r8), intent(in)            :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8), intent(in)            :: part_dens                  ! oven-dry particle density [kg m-2]

      ! LOCALS:
      integer  :: i         ! looping index
      real(r8) :: A_dead    ! total surface area of dead fuels
      real(r8) :: A_live    ! total surface area of live fuels

      ! initialize values
      A_dead = 0.0_r8
      A_live = 0.0_r8
      ! calculate mean total area per unit fuel cell of each size EQ 53 in Rothermel 1972
      do i = 1, num_fuel_classes
         this%mean_total_SA(i) = sav_fuel(i)*this%loading(i)/part_dens
         if (i /= fuel_classes%live_grass()) then
            A_dead = A_dead + this%mean_total_SA(i)
         else
            A_live = A_live + this%mean_total_SA(i)
         end if
      end do
      ! calculate weighting factor, EQ. 56 in Rothermel 1972
      do i = 1, num_fuel_classes
         if(i /= fuel_classes%live_grass()) then
            this%weighting_factor(i) = this%mean_total_SA(i)/A_dead
         else
            this%weighting_factor(i) = this%mean_total_SA(i)/A_live
         end if
      end do

      this%wf_dead = A_dead/(A_dead + A_live)
      this%wf_live = A_live/(A_dead + A_live)

   end subroutine CalculateWeightingFactor

   !-------------------------------------------------------------------------------------

   subroutine SumLoading(this)
      ! DESCRIPTION:
      ! Sums up the loading - excludes trunks
      !
      ! Only the 1-h, 10-h and 100-h fuel classes influence fire spread
      ! Rothermel, 1972 (USDA FS GTR INT-115)
      ! Wilson, 1982 (UTINT-289)
      ! Pyne et al., 1996 (Introduction to wildland fire)
      ! XLG 2025-10-16
      ! fuel load is summed using the calculated weighting factors in
      ! a way that if two fuel classes belong to the same group given SAV (see below)
      ! then the two fuel classes share one weighting factor, which is the sum of the weighting_factor
      ! across the two fuel classes, and weight is 0 if SAV < 0.52 cm-1
      !
      ! the six fuel subgroups are:
      ! SAV >= 1200 ft-1 (39.36 cm-1)
      ! 192 ft-1 (6.3 cm-1) <= SAV < 1200 ft-1 (39.36 cm-1)
      ! 96 ft-1 (3.14 cm-1) <= SAV < 192 ft-1 (6.3 cm-1)
      ! 48 ft-1 (1.57 cm-1) <= SAV < 96 ft-1 (3.14 cm-1)
      ! 16 ft-1 (0.52 cm-1) <= SAV < 48 ft-1 (1.57 cm-1)
      ! SAV < 16 ft-1 (0.52 cm-1)
      ! P 14-15 in Albini 1976a , which is refined version for EQs 59-60 in Rothermel 1972

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      ! LOCALS:
      integer    :: i                             ! looping index

      ! ensure weighting factor is updated
      call this%CalculateWeightingFactor(sav_fuel, part_dens)

      ! get fuel load weighting factor given SAV
      call this%FuelLoadWeight()
      this%non_trunk_loading = 0.0_r8
      this%weighted_loading_dead = 0.0_r8
      this%weighted_loading_live = 0.0_r8
      do i = 1, num_fuel_classes
         if (i /= fuel_classes%trunks()) then
            this%non_trunk_loading = this%non_trunk_loading +  &
               this%loading(i)
         end if

         if(i /= fuel_classes%live_grass())then
            this%weighted_loading_dead = this%weighted_loading_dead + &
               this%fuel_weight(i)*this%loading(i)
         else
            this%weighted_loading_live = this%weighted_loading_live + &
               this%fuel_weight(i)*this%loading(i)
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
      real(r8) :: mef_live                                 ! live fuel moisture of extinction [m3/m3]
      integer  :: i                                        ! looping index

      if (this%non_trunk_loading + this%loading(fuel_classes%trunks()) > nearzero) then
         ! calculate fuel moisture [m3/m3] for each fuel class depending on what
         ! fire weather class is in use
         select type (fireWeatherClass)
          class is (nesterov_index)
            call CalculateFuelMoistureNesterov(sav_fuel, drying_ratio,                   &
               fireWeatherClass%fire_weather_index, moisture)
            this%moisture = moisture
          class default
            write(fates_log(), *) 'Unknown fire weather class selected.'
            write(fates_log(), *) 'Choose a different fire weather class or upate this subroutine.'
            call endrun(msg=errMsg( __FILE__, __LINE__))
         end select

         this%average_moisture_dead = 0.0_r8
         this%average_moisture_live = 0.0_r8
         this%MEF_dead = 0.0_r8
         this%MEF_live = 0.0_r8
         do i = 1, num_fuel_classes
            ! effective_moisture is used for determining fractional fuel burnt per fuel class later
            ! to derive this by fuel category (live vs dead) and by fuel sizes
            ! let's still keep this calculation for both dead and live fuels
            moisture_of_extinction(i) = MoistureOfExtinction(sav_fuel(i))
            this%effective_moisture(i) = this%moisture(i)/moisture_of_extinction(i)
            if(i /= fuel_classes%live_grass())then
               ! average fuel moisture  and MEF for dead fuels [m3/m3] EQ. 66 in Rothermel 1972
               ! trunk is not excluded, but it's getting a small weight
               ! MEF_dead is an input parameter in Rothermel 1972. Since we are calculating it for
               ! each fuel size class given SA:V, I am calculating averaged MEF_dead using the
               ! weighting factor 2025-10-16 XLG
               this%average_moisture_dead = this%average_moisture_dead + &
                  this%weighting_factor(i)*this%moisture(i)
               this%MEF_dead = this%MEF_dead + this%weighting_factor(i)*moisture_of_extinction(i)
            end if
         end do

         ! calculate average moisture and MEF for live fuels. We currently only have live grass fuel
         ! The definition of surface fuel bed, however, includes live woody fuel from 1h to 100h
         ! consider including live fine woody fuels later 2025-10-16 XLG

         this%average_moisture_live = this%average_moisture_live + &
            this%weighting_factor(fuel_classes%live_grass())*this%moisture(fuel_classes%live_grass())

         ! this is the Rothermel way of calculating averaged live fuel MEF across all fuel
         ! classes for live fuels, which is then used to calculate live fuel mositure damping coeff
         ! The difference between the SPITFIRE version and original Rothermel model
         ! is that live and dead fuels are always treated separately in the original model
         call LiveFuelMoistureOfExtinction(this%loading, sav_fuel, this%moisture, &
            this%MEF_dead, mef_live)
         this%MEF_live = mef_live
      else
         this%effective_moisture(1:num_fuel_classes) = 0.0_r8
         this%average_moisture_dead = 0.0_r8
         this%average_moisture_live = 0.0_r8
         this%MEF_dead = 0.0_r8
         this%MEF_live = 0.0_r8
      end if

   end subroutine UpdateFuelMoisture

   !-------------------------------------------------------------------------------------

   subroutine CalculateFuelMoistureNesterov(sav_fuel, drying_ratio, NI, moisture)
      !
      ! DESCRIPTION:
      !   Updates fuel moisture

      ! ARGUMENTS:
      real(r8), intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8), intent(in)    :: drying_ratio               ! drying ratio
      real(r8), intent(in)    :: NI                         ! Nesterov Index
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
      ! In Albini 1976a and Rothermel 1972 MEF_dead is an input parameter
      ! MEF_live is then calculated using P16 in Albini 1976a, which is a function of
      ! MEF_dead, dead-to-live load ratio and dead fuel moisture
      ! We'll keep this function for calculating MEF dead for each fuel class and then
      ! calculate the weighted mean MEF_dead given weighting_factor
      ! XLG 2025-10-16

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
   subroutine LiveFuelMoistureOfExtinction(fuel_load, sav_val, fmc, mef_dead, mef_live)
      !
      ! DESCRIPTION
      ! Calculate average MEF for live fuels using P. 16 in Albini 1976a, whiched replaced EQ. 88
      ! in Rothermel 1972; a unit conversion for sav has been applied to convert [cm-1 ----> ft-1 ]
      !
      use FatesLitterMod
      ! ARGUMENTS
      real(r8), intent(in)           :: fuel_load(num_fuel_classes)   ! fuel load [kgC m-2]
      real(r8), intent(in)           :: sav_val(num_fuel_classes)     ! fuel surface area to volume ratio [cm-1]
      real(r8), intent(in)           :: fmc(num_fuel_classes)         ! fuel moisture content by size class [m3/m3]
      real(r8), intent(in)           :: mef_dead                      ! weighted mean of moisture of extinction for dead fuel
      real(r8), intent(out)          :: mef_live                      ! live fuel moisture of extinction
      ! LOCALS
      real(r8)          :: W           ! dead-to-live load ratio
      real(r8)          :: moist_dead  ! dead fuel moisture
      real(r8)          :: w_n, md_n   ! nominator for calculating W and moist_dead
      real(r8)          :: w_d, md_d   ! denominator for calculating W and moist_dead

      integer            :: i           ! looping index

      w_n  = 0.0_r8
      md_n = 0.0_r8
      w_d  = 0.0_r8
      md_d = 0.0_r8
      do i = 1, num_fuel_classes
         if(i /= fuel_classes%live_grass())then
            w_n = w_n + fuel_load(i)*exp(-4.5276_r8/sav_val(i))
            md_n = md_n + fmc(i)*fuel_load(i)*exp(-4.5276_r8/sav_val(i))
            md_d = md_d + fuel_load(i)*exp(-4.5264_r8/sav_val(i))
         else
            w_d = w_d + fuel_load(i)*exp(-16.4042_r8/sav_val(i))
         end if

      end do
      W = w_n/w_d
      moist_dead = md_n/md_d
      mef_live = 2.9_r8 * W * (1.0_r8 - moist_dead/mef_dead) - 0.226_r8

   end subroutine LiveFuelMoistureOfExtinction

   !-------------------------------------------------------------------------------------

   subroutine AverageBulkDensity(this, bulk_density)
      ! DESCRIPTION:
      !   Calculates average bulk density, which is weighted toward fine fuels (1-100-h fuels)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread
      !    Rothermel, 1972 (USDA FS GTR INT-115)
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this                           ! fuel class
      real(r8),           intent(in)    :: bulk_density(num_fuel_classes) ! bulk density of all fuel types [kg/m3]

      ! LOCALS:
      integer     :: i          ! looping index
      real(r8)    :: fbd_dead   ! mean fuel bulk density for dead fuels [kg/m3]
      real(r8)    :: fbd_live   ! mean fuel bulk density for live fuels [kg/m3]

      if (this%non_trunk_loading > nearzero) then
         fbd_dead = 0.0_r8
         fbd_live = 0.0_r8
         do i = 1, num_fuel_classes
            if(i /= fuel_classes%live_grass())then
               fbd_dead = fbd_dead + this%weighting_factor(i)*bulk_density(i)
            else
               fbd_live = fbd_live + this%weighting_factor(i)*bulk_density(i)
            end if

         end do
         ! average bulk density between dead and live fuels using wf_dead and wf_live
         ! same as how SAV is calculated for the entire fuel bed
         this%bulk_density_weighted = this%wf_dead*fbd_dead + this%wf_live*fbd_live
      else
         this%bulk_density_weighted = sum(bulk_density(1:num_fuel_classes))/num_fuel_classes
      end if

   end subroutine AverageBulkDensity

   !-------------------------------------------------------------------------------------

   subroutine AverageSAV(this, sav_fuel)
      ! DESCRIPTION:
      !   Calculates average surface area to volume ratio, weighted toward fine fuels (1-100-h fuels)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread
      !    Rothermel, 1972 (USDA FS GTR INT-115)
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                       ! fuel class
      real(r8),         intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]

      ! LOCALS:
      integer       :: i              ! looping index
      real(r8)      :: sav_dead       ! average SA:V of dead fuels
      real(r8)      :: sav_live       ! average SA:V of live fuels

      if (this%non_trunk_loading > nearzero) then
         sav_dead = 0.0_r8
         sav_live = 0.0_r8
         do i = 1, num_fuel_classes
            ! average bulk density across all fuel types, EQ. 72 in Rothermel 1972
            if (i /= fuel_classes%live_grass()) then
               sav_dead = sav_dead + this%weighting_factor(i)*sav_fuel(i)
            else
               sav_live = sav_live + this%weighting_factor(i)*sav_fuel(i)
            end if
         end do
         ! averaged SA:V between dead and live fuels, EQ 71 in Rothermel 1972
         this%SAV_weighted = this%wf_dead*sav_dead + this%wf_live*sav_live
      else
         this%SAV_weighted = sum(sav_fuel(1:num_fuel_classes))/num_fuel_classes
      end if

   end subroutine AverageSAV

   !---------------------------------------------------------------------------------------

   subroutine CalculateFuelBurnt(this, fuel_consumed)
      ! DESCRIPTION:
      !   Calculates the fraction and total amount of fuel burnt
      !

      use SFParamsMod, only : SF_val_mid_moisture, SF_val_mid_moisture_Coeff
      use SFParamsMod, only : SF_val_mid_moisture_Slope, SF_val_min_moisture
      use SFParamsMod, only : SF_val_low_moisture_Coeff, SF_val_low_moisture_Slope
      use SFParamsMod, only : SF_val_miner_total, SF_val_SAV, SF_val_part_dens

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                            ! fuel class
      real(r8),         intent(out)   :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]

      ! LOCALS:
      real(r8) :: rel_moisture                    ! relative moisture of fuel (moist/moisture of extinction) [unitless]
      integer  :: i                               ! looping index

      ! CONSTANTS:
      real(r8), parameter :: max_grass_frac = 0.8_r8 ! maximum fraction burnt for live grass fuels

      this%frac_burnt(:) = 1.0_r8
      ! get all the weighting factor
      call this%CalculateWeightingFactor(SF_val_SAV, SF_val_part_dens)
      call this%FuelLoadWeight()

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
         ! we don't use the raw loading (which assum a weight of 1 for each fuel class)
         ! switch to using the weighted fuel load per size class 2025-10-16 XLG
         ! see p15, p17 in Albini 1976a for the calculation of reaction intensity based on
         ! weighted fuel load. so we also use the weighted fuel load here just to be consistent
         fuel_consumed(i) = this%frac_burnt(i)*this%loading(i)*this%fuel_weight(i)
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


   !-------------------------------------------------------------------------------------
   subroutine FuelLoadWeight(this)
      !
      ! DESCRIPTION
      ! Assign each fuel class to one fuel group given user-defined SA:V
      ! calculate weighting factor for each fuel class, which will be used
      ! to calculate sum fuel load; we only do this for dead fuel since
      ! dead and live fuels are treated separately and we currently only have
      ! one live fuel, weight should always be 1
      !
      use SFParamsMod, only : SF_val_SAV
      ! ARGUMENTS:
      class(fuel_type), intent(inout)  :: this  ! fuel class
      ! LOCALS:
      integer                  :: ng                             ! sizes
      integer                  :: i, g                           ! looping indices
      integer                  :: group_idx(num_fuel_classes)    ! fuel group assigned
      real(r8), allocatable    :: group_sum(:)                   ! summed weight across fuel sizes that are within the same subgroup
      real(r8), parameter, dimension(5)      :: sav_edg=(/0.53_r8, 1.58_r8, 3.15_r8, 6.30_r8, 39.37_r8/)

      ng = size(sav_edg)
      ! initialize fuel group with group 1 (sav < 0.53)
      group_idx(:) = 1
      ! assign fuel class to one of the fuel groups
      ! for each sav edge value, if sav > sav_edg, increment by 1
      ! at the corresponding location for group_idx
      do g = 1, ng
         where (SF_val_SAV >= sav_edg(g)) group_idx = group_idx + 1
      end do

      allocate(group_sum(ng+1))
      group_sum(:) = 0.0_r8

      ! sum weighting factor across the same subgroup
      ! only using dead fuels now
      ! as this is calculated for dead and live category separately
      ! but dead and live fuels should have same SA:V (which should only be
      ! influenced by size class not dead/live category).
      ! We only have one live grass fuel, which share one SAV value with 1h dead
      do i = 1, num_fuel_classes
         if(i /= fuel_classes%live_grass())then
            group_sum(group_idx(i)) = group_sum(group_idx(i)) + this%weighting_factor(i)
         end if
      end do
      ! map the fuel load weighting factor to each fuel class
      ! if group is 1 (sav < 0.5) then fuel_weight = 0
      ! also live grass fuel has a weight of 1 since we now only have one
      ! fuel size in live fuel category
      do i = 1, num_fuel_classes
         if(group_idx(i) == 1) then
            this%fuel_weight(i) = group_sum(group_idx(i))*0.0_r8
         else if(i == fuel_classes%live_grass())then
            this%fuel_weight(i) = 1.0_r8
         else
            this%fuel_weight(i) = group_sum(group_idx(i))
         end if
      end do

      deallocate(group_sum)

   end subroutine FuelLoadWeight


end module FatesFuelMod
