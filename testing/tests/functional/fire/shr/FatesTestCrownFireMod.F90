
module FatesTestCrownFireMod
   !
   ! DESCRIPTION
   ! modules to support testing crown fire
   !

   use FatesConstantsMod,      only : r8 => fates_r8, nearzero
   use FatesUnitTestIOMod,     only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims
   use FatesUnitTestIOMod,     only : RegisterVar, EndNCDef, WriteVar
   use FatesUnitTestIOMod,     only : type_double, type_char, type_int
   use FatesFuelMod,           only : fuel_type
   use SFParamsMod,            only : SF_val_miner_total, SF_val_part_dens
   use SFEquationsMod,         only : OptimumPackingRatio, ReactionIntensity
   use SFEquationsMod,         only : HeatofPreignition, EffectiveHeatingNumber
   use SFEquationsMod,         only : WindFactor, PropagatingFlux
   use SFEquationsMod,         only : ForwardRateOfSpread

   implicit none
   private

   public :: ROSWrapper
   public :: EffectiveWindWrapper
   public :: WriteCanopyFuelData

   ! ======================================================================================

contains

   subroutine ROSWrapper(fuelType, fire_weather, effect_wind, ROS, i_r)
      !
      ! DESCRIPTION
      ! Calls a sequence of functions to calculate surface fire rate of spread
      !
      ! ARGUMENTS:
      type(fuel_type), intent(in)     :: fuelType     ! fuel type
      real(r8), intent(in)            :: fire_weather ! fire weather index [unitless]
      real(r8), intent(in)            :: effect_wind  ! effective wind speed [m/min]
      real(r8), intent(out)           :: ROS          ! forward rate of spread [m/min]
      real(r8), intent(out)           :: i_r          ! reaction intensity [kJ/m2/min], required for calculating HPA
      !
      ! LOCALS:
      real(r8)                        :: beta         ! packing ratio [unitless]
      real(r8)                        :: beta_op      ! optimum packing ratio [unitless]
      real(r8)                        :: beta_ratio   ! relative packing ratio [unitless]
      real(r8)                        :: xi           ! propagating flux ratio [unitless]
      real(r8)                        :: eps          ! effective heating number [unitless]
      real(r8)                        :: phi_wind     ! wind factor [unitless]
      real(r8)                        :: q_ig         ! heat of pre-ignition [kJ/kg]
      real(r8)                        :: fuel_net     ! fuel load excluding mineral content [kgC/m2]

      ! fraction of fuel array occupied by fuels and optimum packing ratio
      beta = fuelType%bulk_density_notrunks/SF_val_part_dens
      beta_op = OptimumPackingRatio(fuelType%SAV_notrunks)
      ! relative packing ratio
      if (beta_op < nearzero) then
         beta_ratio = 0.0_r8
      else
         beta_ratio = beta/beta_op
      end if

      ! remove mineral content
      fuel_net = fuelType%non_trunk_loading*(1.0_r8 - SF_val_miner_total)

      ! reaction intensity
      i_r = ReactionIntensity(fuel_net/0.45_r8, fuelType%SAV_notrunks, beta_ratio,  &
         fuelType%average_moisture_notrunks, fuelType%MEF_notrunks)

      ! heat of preignition
      q_ig = HeatofPreignition(fuelType%average_moisture_notrunks)

      ! effective heating number
      eps = EffectiveHeatingNumber(fuelType%SAV_notrunks)

      ! wind factor
      phi_wind = WindFactor(effect_wind, beta_ratio, fuelType%SAV_notrunks)

      ! propogation flux
      xi = PropagatingFlux(beta, fuelType%SAV_notrunks)

      ! forward rate of spread [m/min]
      ROS = ForwardRateOfSpread(fuelType%bulk_density_notrunks, eps, q_ig, i_r, xi, phi_wind)

   end subroutine ROSWrapper

   ! ======================================================================================


   subroutine EffectiveWindWrapper(tree_area, grass_area, area, wind, effect_wind)
      !
      ! DESCRIPTION
      ! Call a sequence of calculations to get effective wind speed
      !
      ! ARGUMENTS
      real(r8), intent(in)                  :: tree_area             ! total tree area of the test patch [m2]
      real(r8), intent(in)                  :: grass_area            ! total grass area of the test patch [m2]
      real(r8), intent(in)                  :: area                  ! patch area [m2]
      real(r8), intent(in)                  :: wind                  ! wind speed [m/min]
      real(r8), intent(out)                 :: effect_wind           ! returned effective wind speed [m/min]
      !
      ! LOCALS:
      real(r8)    :: tree_frac          ! tree fraction of the test patch [0-1]
      real(r8)    :: grass_frac         ! grass fraction of the test patch [0-1]
      real(r8)    :: bare_frac          ! bare ground fraction of the test patch [0-1]
      !
      ! CONSTANTS:
      real(r8), parameter   :: wind_atten_tree = 0.4_r8 ! wind attenuation factor for tree fraction
      real(r8), parameter   :: wind_atten_grass = 0.6_r8 ! wind attenuation factor for grass fraction

      tree_frac = tree_area / area
      grass_frac = grass_area / area
      grass_frac = min(grass_frac, 1.0_r8 - tree_frac)
      bare_frac = 1.0_r8 - tree_frac - grass_frac

      effect_wind = wind * (tree_frac*wind_atten_tree +         &
         (grass_frac + bare_frac)*wind_atten_grass )


   end subroutine EffectiveWindWrapper

   ! ======================================================================================

   subroutine WriteCanopyFuelData(out_file, nfuelmods, nstands, nwind, Wind, NI, CWC, &
      CBD, CBH, canopy_fuel_load, ROS_front, FI, FI_init, ROS_actfm10, ROS_critical, ROS_actCI, &
      CFB, ROS_final, FI_final, fuel_models,  patch_types)
      !
      ! DESCRIPTION:
      ! write out data from canopy fuel functional test
      !

      ! ARGUMENTS:
      character(len=*),   intent(in) :: out_file
      integer,            intent(in) :: nfuelmods
      integer,            intent(in) :: nstands
      integer,            intent(in) :: nwind
      real(r8),           intent(in) :: Wind(:)
      real(r8),           intent(in) :: NI(:)
      real(r8),           intent(in) :: CWC(:)
      real(r8),           intent(in) :: CBD(:)
      real(r8),           intent(in) :: CBH(:)
      real(r8),           intent(in) :: canopy_fuel_load(:)
      real(r8),           intent(in) :: ROS_front(:,:,:,:)
      real(r8),           intent(in) :: FI(:,:,:,:)
      real(r8),           intent(in) :: FI_init(:,:)
      real(r8),           intent(in) :: ROS_actfm10(:,:,:,:,:)
      real(r8),           intent(in) :: ROS_critical(:)
      real(r8),           intent(in) :: ROS_actCI(:,:,:,:,:)
      real(r8),           intent(in) :: CFB(:,:,:,:,:)
      real(r8),           intent(in) :: ROS_final(:,:,:,:,:)
      real(r8),           intent(in) :: FI_final(:,:,:,:,:)
      integer,            intent(in) :: fuel_models(:)
      integer,            intent(in) :: patch_types(:)

      ! LOCALS:
      integer           :: ncid         ! netcdf id
      character(len=20) :: dim_names(5) ! dimension names
      integer           :: dimIDs(5)    ! dimension IDs
      integer           :: modID, patchID, windID, niID, cwcID
      integer           :: CBDID, CBHID, cflID
      integer           :: ROS_frontID, ROS_actID,ROS_ciID, ROS_minID, ROS_finID
      integer           :: FIID, FI_initID, FI_finID
      integer           :: CFBID

      ! dimension names
      dim_names = [character(len=20) :: 'wind','fire_weather','canopy_water', &
         'patch_type', 'fuel_model']
      ! open file
      call OpenNCFile(trim(out_file), ncid, 'readwrite')

      ! register dimensions
      call RegisterNCDims(ncid, dim_names, (/nwind, size(NI), size(CWC), &
         nstands, nfuelmods/), 5, dimIDs)

      ! first register dimension variables

      ! register wind speed
      call RegisterVar(ncid, 'wind', dimIDs(1:1), type_int,    &
         [character(len=20) :: 'units', 'long_name' ],    &
         [character(len=150) :: 'm min-1', 'wind speed index'], 2, windID)

      ! register fire weather index
      call RegisterVar(ncid, 'fire_weather', dimIDs(2:2), type_int,    &
         [character(len=20) :: 'units', 'long_name' ],    &
         [character(len=150) :: '', 'fire weather index'], 2, niID)

      ! register canopy water content
      call RegisterVar(ncid, 'canopy_water', dimIDs(3:3), type_int,    &
         [character(len=20) :: 'units', 'long_name' ],    &
         [character(len=150) :: '%', 'canopy water index'], 2, cwcID)

      ! register patch types
      call RegisterVar(ncid, 'patch_type', dimIDs(4:4), type_int,     &
         [character(len=20) :: 'units', 'long_name' ],     &
         [character(len=150) :: '', 'patch type index'], 2, patchID)

      ! register fuel models
      call RegisterVar(ncid, 'fuel_model', dimIDs(5:5), type_int,     &
         [character(len=20) :: 'units', 'long_name' ],     &
         [character(len=150) :: '', 'fuel model index'], 2, modID)

      ! register variables

      ! register actual variables

      ! register CBD
      call RegisterVar(ncid, 'CBD', dimIDs(4:4), type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'patch_type', 'kg m-3', 'canopy bulk density'],  &
         3, CBDID)

      ! register CBH
      call RegisterVar(ncid, 'CBH', dimIDs(4:4), type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'patch_type', 'm', 'canopy base height'],  &
         3, CBHID)

      ! register canopy fuel load
      call RegisterVar(ncid, 'Canopy_fuel', dimIDs(4:4), type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'patch_type', 'kg m-2', 'canopy fuel density'],  &
         3, cflID)

      ! register surface fire spread rate
      call RegisterVar(ncid, 'ROS_front', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type fuel_model', 'm min-1', &
         'surface forward rate of spread'], 3, ROS_frontID)

      ! register surface fire intensity
      call RegisterVar(ncid, 'FI', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type fuel_model', 'kW m-1',  &
         'surface fire intensity'], 3, FIID)

      ! register min FI to initiate crown fire
      call RegisterVar(ncid, 'FI_init', (/dimIDs(4), dimIDs(3)/), type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'patch_type canopy_water', 'kW m-1',         &
         'fire intensity to ignite crown fuel'], 3, FI_initID)

      ! register active crown fire ROS
      call RegisterVar(ncid, 'ROS_active', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(3), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type canopy_water fuel_model',  &
         'm min-1', 'active crown fire ROS assuming FM10'], 3, ROS_actID)

      ! register active crown fire ROS at CI wind speed
      call RegisterVar(ncid, 'ROS_act_ci', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(3), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type canopy_water fuel_model',  &
         'm min-1', 'active crown fire ROS assuming FM10 at crowning index ws'], 3, ROS_ciID)

      ! register critical ROS to sustain active crown fire
      call RegisterVar(ncid, 'ROS_min', dimIDs(4:4), type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'patch_type', 'm min-1', 'min ROS to sustain active crown fire'],  &
         3, ROS_minID)

      ! register crown fraction burnt
      call RegisterVar(ncid, 'CFB', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(3), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type canopy_water fuel_model', '', &
         'crown fraction burnt by crown fire'], 3, CFBID)

      ! register final ROS
      call RegisterVar(ncid, 'ROS_final', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(3), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type canopy_water fuel_model', 'm min-1', &
         'ROS after crown fire update'], 3, ROS_finID)

      ! register final FI
      call RegisterVar(ncid, 'FI_final', (/dimIDs(1), dimIDs(2), dimIDs(4), dimIDs(3), dimIDs(5)/), &
         type_double,        &
         [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
         [character(len=150) :: 'wind fire_weather patch_type canopy_water fuel_model',  &
         'kW m-1', 'FI after crown fire update'],  3, FI_finID)

      ! finish defining variables
      call EndNCDef(ncid)

      ! write out data
      call WriteVar(ncid, patchID, patch_types(:))
      call WriteVar(ncid, modID, fuel_models(:))
      call WriteVar(ncid, windID, Wind(:))
      call WriteVar(ncid, niID, NI(:))
      call WriteVar(ncid, cwcID, CWC(:))
      call WriteVar(ncid, CBDID, CBD(:))
      call WriteVar(ncid, CBHID, CBH(:))
      call WriteVar(ncid, cflID, canopy_fuel_load(:))
      call WriteVar(ncid, ROS_frontID, ROS_front(:,:,:,:))
      call WriteVar(ncid, FIID, FI(:,:,:,:))
      call WriteVar(ncid, FI_initID, FI_init(:,:))
      call WriteVar(ncid, ROS_actID, ROS_actfm10(:,:,:,:,:))
      call WriteVar(ncid, ROS_minID, ROS_critical(:))
      call WriteVar(ncid, ROS_ciID, ROS_actCI(:,:,:,:,:))
      call WriteVar(ncid, CFBID, CFB(:,:,:,:,:))
      call WriteVar(ncid, ROS_finID, ROS_final(:,:,:,:,:))
      call WriteVar(ncid, FI_finID, FI_final(:,:,:,:,:))

      call CloseNCFile(ncid)

   end subroutine WriteCanopyFuelData



end module FatesTestCrownFireMod
