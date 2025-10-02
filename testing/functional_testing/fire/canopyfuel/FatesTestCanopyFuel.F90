program FatesTestCanopyFuel

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use FatesFactoryMod,             only : InitializeGlobals, GetSyntheticPatch
  use SyntheticPatchTypes,         only : synthetic_patch_array_type
  use FatesTestFireMod,            only : SetUpFuel
  use SFFireWeatherMod,            only : fire_weather
  use SFNesterovMod,               only : nesterov_index
  use SyntheticFuelModels,         only : fuel_models_array_class
  use FatesFuelClassesMod,         only : num_fuel_classes, fuel_classes
  use FatesFuelMod,                only : fuel_type
  use FatesTestCrownFireMod,       only : ROSWrapper, EffectiveWindWrapper, WriteCanopyFuelData
  use CrownFireEquationsMod,       only : MaxHeight , BiomassBin
  use CrownFireEquationsMod,       only : PassiveCrownFireIntensity, CrownFireBehaveFM10
  use CrownFireEquationsMod,       only : HeatReleasePerArea, CrownFireCFB
  use CrownFireEquationsMod,       only : CrownFireIntensity
  use SFEquationsMod,              only : ReactionIntensity, FireIntensity
  use FatesAllometryMod,           only : CrownDepth
  use PRTGenericMod,               only : leaf_organ
  use PRTGenericMod,               only : sapw_organ
  use PRTGenericMod,               only : struct_organ
  use PRTGenericMod,               only : carbon12_element
  use PRTParametersMod,            only : prt_params
  use FatesConstantsMod,           only : itrue
  use SFParamsMod,                 only : SF_val_CWD_frac
  use FatesLitterMod,              only : adjust_SF_CWD_frac
  use FatesLitterMod,              only : ncwd
  use SFParamsMod,                 only : SF_val_SAV
  use SFParamsMod,                 only : SF_val_FBD, SF_val_fire_threshold
  use SFParamsMod,                 only : SF_val_part_dens, SF_val_miner_total


  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader ! param reader instance
  type(synthetic_patch_array_type)               :: patch_data   ! array of synthetic patches
  type(fuel_models_array_class)                  :: fuel_models_array    ! array of fuel models
  class(fire_weather),               pointer     :: fireWeather          ! fire weather object
  type(fuel_type),                   allocatable :: fuel(:)              ! fuel objects  
  character(len=:),                  allocatable :: param_file   ! input parameter file
  character(len=100),                allocatable :: fuel_names(:)        ! names of fuel models
  character(len=2),                  allocatable :: carriers(:)      ! carriers of fuel models
  character(len=100),                allocatable :: patch_names(:)       ! names of patch types
  type(fates_patch_type),            pointer     :: patch        ! patch
  type(fates_cohort_type),           pointer     :: cohort       ! cohort
  
  real(r8)                                       :: effect_wind  ! effective wind speed [m/min]
  real(r8)                                       :: max_height   ! max cohort height [m]
  real(r8)                                       :: leaf_c       ! leaf carbon [kg C]
  real(r8)                                       :: sapw_c       ! sap wood carbon [kg C]
  real(r8)                                       :: struct_c     ! strcutural carbon [kg C]
  real(r8)                                       :: woody_c      ! aboveground sap wood + structural carbon [kg C]
  real(r8)                                       :: canopy_fuel_1h ! leaf + 1-hour woody fuel [kg C]
  real(r8),               allocatable            :: biom_matrix(:) ! array to hold biomass at 1m interval [kg biomass]
  real(r8)                                       :: crown_depth  ! crown length of a cohort [m]
  real(r8)                                       :: cbh_co       ! cohort base height, different from path level base height [m]
  real(r8)                                       :: cwd_frac_adj(ncwd) ! adjusted fractional allocation of woody biomass to coarse wood debris pool
  real(r8)                                       :: ROS                 ! surface fire forward rate of spread [m/min]
  real(r8)                                       :: fuel_consumed(num_fuel_classes) ! fuel consumed by fuel class [kgC/m2]
  real(r8)                                       :: tfc_ros             ! overall fuel consumed by spreading fire ignoring trunks [kgC/m2]
  real(r8)                                       :: ROS_active          ! active crown fire ROS using FM10 [m/min]
  real(r8)                                       :: CI                  ! crowning index [m/min]
  real(r8)                                       :: i_r                 ! reaction intensity [kJ/m2/min]
  real(r8)                                       :: HPA                 ! heat release per area [kW/m2]
  real(r8)                                       :: ROS_init            ! ROS for initiating crown fire [m/min]
  real(r8)                                       :: crown_frac_burnt    ! crown fraction burnt by crown fire [0-1]
  real(r8),               allocatable            :: FI(:,:)             ! fire intensity by each fuel model and stand type [kW/m]
  real(r8),               allocatable            :: FI_init(:)          ! initiation fire intensity by stand type [kW/m]
  real(r8),               allocatable            :: canopy_fuel_load(:) ! patch level canopy fuel load by stand type [kg biomass]
  real(r8),               allocatable            :: CBD(:)              ! patch level canopy bulk density [kg biomass / m3]
  real(r8),               allocatable            :: CBH(:)              ! patch level canopy base height [m]
  real(r8),               allocatable            :: ROS_front(:,:)      ! surface fire forward rate of spread by fuel model and stand type [m/min]
  real(r8),               allocatable            :: ROS_actfm10(:,:)    ! active crown fire ROS by drying ratio and wind speed [m/min]
  real(r8),               allocatable            :: ROS_critical(:)     ! critical ROS for active crown fire to occur, by stand type [m/min]
  real(r8),               allocatable            :: ROS_final(:,:)      ! final ROS by stand type and surface fuel model [m/min]
  real(r8),               allocatable            :: FI_final(:,:)       ! final fire intensity by stand type and surface fuel model [kW/m]
  real(r8),               allocatable            :: CFB(:,:)            ! crown fraction burned by stand type and surface fuel model [fraction]
  integer                                        :: p, f, i             ! looping indices
  integer                                        :: num_fuel_models      ! number of fuel models to test
  integer                                        :: num_patch_types      ! number of patch types to test
  integer                                        :: active_crownfire     ! 1 = active crown fire
  integer                                        :: passive_crownfire    ! 1 = passive crown fire

  ! CONSTANTS:
  integer,  parameter :: num_levsoil = 10      ! number of soil layers
  real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]
  real(r8), parameter :: biomass_2_carbon = 0.45_r8  ! biomass to carbon multiplier
  real(r8), parameter :: NI = 5000.0_r8        ! fire weather index to use
  real(r8), parameter :: wind = 500.0_r8       ! wind speed to use [m/min]
  real(r8), parameter :: CWC = 70.0_r8         ! canopy water content [%]
  real(r8), parameter :: drying_ratio = 3000.0_r8
  real(r8), parameter :: wind_atten_tree = 0.4_r8               ! wind attenuation factor for tree fraction
  real(r8), parameter :: wind_atten_grass = 0.6_r8              ! wind attenuation factor for grass fraction
  character(len=*), parameter :: out_file = 'canopyfuel_out.nc' ! output file
  
  ! fuel models and patch types to test
  integer, parameter, dimension(2) :: fuel_models = (/10, 11/)
  integer, parameter, dimension(4) :: patch_ids = (/1, 6, 7, 8/)
  
  ! number of fuel models  and patch types to test
  num_fuel_models = size(fuel_models)
  num_patch_types = size(patch_ids)

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! allocate arrays
  allocate(fuel_names(num_fuel_models))
  allocate(carriers(num_fuel_models))
  allocate(patch_names(num_patch_types))
  allocate(FI(num_patch_types, num_fuel_models))
  allocate(FI_init(num_patch_types))
  allocate(canopy_fuel_load(num_patch_types))
  allocate(CBD(num_patch_types))
  allocate(CBH(num_patch_types))
  allocate(ROS_front(num_patch_types, num_fuel_models))
  allocate(ROS_actfm10(num_patch_types, num_fuel_models))
  allocate(ROS_critical(num_patch_types))
  allocate(ROS_final(num_patch_types, num_fuel_models))
  allocate(FI_final(num_patch_types, num_fuel_models))
  allocate(CFB(num_patch_types, num_fuel_models))

  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()

  ! set up fire weather class
  ! here we just assign a pre-defined value
  allocate(nesterov_index :: fireWeather)
  call fireWeather%Init()
  fireWeather%fire_weather_index = NI

  ! set up fuel objects and calculate loading
  allocate(fuel(num_fuel_models))
  call fuel_models_array%GetFuelModels()

  do f = 1, num_fuel_models
        ! uses data from fuel_models to initialize fuel
        call SetUpFuel(fuel(f), fuel_models_array, fuel_models(f), fuel_names(f), carriers(f))

        ! sum up fuel and calculate loading
        call fuel(f)%SumLoading()
        call fuel(f)%CalculateFractionalLoading()

        ! calculate geometric properties
        call fuel(f)%AverageBulkDensity_NoTrunks(SF_val_FBD)
        call fuel(f)%AverageSAV_NoTrunks(SF_val_SAV)

        ! calculate fuel moisture content
        call fuel(f)%UpdateFuelMoisture(SF_val_SAV, drying_ratio, fireWeather)
  end do

  ! initialize some global data we need
  call InitializeGlobals(step_size)
  
  ! get all the patch data
  call patch_data%GetSyntheticPatchData()
  
  do p = 1, num_patch_types
        i = patch_data%PatchDataPosition(patch_id=patch_ids(p))
        call GetSyntheticPatch(patch_data%patches(i), num_levsoil, patch)
    
        ! update patch tree and grass area
        call patch%UpdateTreeGrassArea()

        ! effective wind speed 
        call EffectiveWindWrapper(patch%total_tree_area, patch%total_grass_area, &
        patch%area, wind, effect_wind)


        ! Calculate ROS for each fuel model
        do f = 1, num_fuel_models
            call ROSWrapper(fuel(f)%bulk_density_notrunks, fuel(f)%SAV_notrunks,    &
            fuel(f)%non_trunk_loading, fuel(f)%average_moisture_notrunks,   &
            fuel(f)%MEF_notrunks, NI, effect_wind, ROS, i_r)
            ROS_front(p,f) = ROS

        end do

        cohort => patch%tallest
        ! initialize max_height
        max_height = 0.0_r8
        ! search for max height
        do while (associated(cohort))
            if (prt_params%woody(cohort%pft) == itrue) then 
                call MaxHeight(cohort%height, max_height)
            end if
            cohort => cohort%shorter
        end do
    
        ! allocate and initialize biom_martix
        allocate(biom_matrix(0:int(max_height)))
        biom_matrix = 0.0_r8

        ! derive canopy fuel load 
        cohort => patch%tallest
        do while (associated(cohort))
            if(prt_params%woody(cohort%pft) == itrue) then
                call CrownDepth(cohort%height, cohort%pft, crown_depth)
                cbh_co = cohort%height - crown_depth
                leaf_c = cohort%prt%GetState(leaf_organ, carbon12_element)
                sapw_c = cohort%prt%GetState(sapw_organ, carbon12_element)
                struct_c = cohort%prt%GetState(struct_organ, carbon12_element)
                woody_c = cohort%n * prt_params%allom_agb_frac(cohort%pft)* &
                (struct_c + sapw_c) / biomass_2_carbon
                leaf_c = cohort%n * leaf_c / biomass_2_carbon
                call adjust_SF_CWD_frac(cohort%dbh, ncwd, SF_val_CWD_frac, cwd_frac_adj)

                ! update canopy fuel load, this has nothing to do with fuel model
                ! but we save it for each fuel model for later use
                do f = 1, num_fuel_models
                    call fuel(f)%CalculateCanopyFuelLoad(leaf_c, woody_c,   &
                    cwd_frac_adj, canopy_fuel_1h)
                end do

                cohort%canopy_fuel_1h = canopy_fuel_1h
        
                ! 1m biomass bin
                call BiomassBin(cbh_co, cohort%height, crown_depth, canopy_fuel_1h, biom_matrix)
            end if     
            cohort => cohort%shorter 
        end do
  
        ! calculate canopy bulk density
        biom_matrix = biom_matrix / patch%area ! kg biomass / m3

        do f = 1, num_fuel_models
            call fuel(f)%CalculateCanopyBulkDensity(biom_matrix, max_height)
        end do
  
        ! save canopy fuel characteristics
        canopy_fuel_load(p) = fuel(1)%canopy_fuel_load
        CBD(p) = fuel(1)%canopy_bulk_density
        CBH(p) = fuel(1)%canopy_base_height

        deallocate(biom_matrix)

        ! calculate initiation surface fire intensity and critical ROS 
        FI_init(p) = PassiveCrownFireIntensity(fuel(1)%canopy_base_height, CWC)
        ROS_critical(p) = 3.0_r8 / CBD(p)

        ! calculate surface fire intensity for each fuel model and 
        ! calculate ROS_active, ROS_final, FI_final, CFB only when FI > FI_init
        do f = 1, num_fuel_models
            call fuel(f)%CalculateFuelBurnt(fuel_consumed)
            tfc_ros = sum(fuel_consumed) - fuel_consumed(fuel_classes%trunks())
            FI(p,f) = FireIntensity(tfc_ros/0.45_r8, ROS_front(p,f)/60.0_r8)

            if(FI(p,f) > SF_val_fire_threshold .and. FI(p,f) > FI_init(p))then
                call CrownFireBehaveFM10(drying_ratio, NI, &
                SF_val_miner_total, SF_val_part_dens, wind, &
                fuel(f)%canopy_bulk_density, ROS_active, CI)
                ROS_actfm10(p,f) = ROS_active

                ! calculate effective wind speed using CI
                call EffectiveWindWrapper(patch%total_tree_area, patch%total_grass_area, &
                patch%area, CI, effect_wind)

                ! calculate ROS_SA, which is the returned ROS
                call ROSWrapper(fuel(f)%bulk_density_notrunks, fuel(f)%SAV_notrunks,    &
                fuel(f)%non_trunk_loading, fuel(f)%average_moisture_notrunks,   &
                fuel(f)%MEF_notrunks, NI, effect_wind, ROS, i_r)

                write (*,*) 'ROS_SA now is: ', ROS
                write (*,*) 'ROS_front now is: ', ROS_front(p,f)

                HPA = HeatReleasePerArea(fuel(f)%SAV_notrunks, i_r)
                ROS_init = (60.0_r8 * FI_init(p)) / HPA
                ! calculate crown fraction burnt
                call CrownFireCFB(ROS_active, ROS_critical(p), ROS_front(p,f), &
                ROS_init, ROS, active_crownfire, passive_crownfire, crown_frac_burnt)
                CFB(p,f) = crown_frac_burnt
            else
                ROS_actfm10(p,f) = 0.0_r8
                crown_frac_burnt = 0.0_r8
                CFB(p,f) = 0.0_r8
            end if

            ! update ROS and FI 
            if(crown_frac_burnt > 0.0_r8)then
                ROS_final(p,f) = ROS_front(p,f) + crown_frac_burnt*(ROS_actfm10(p,f) - ROS_front(p,f))
                FI_final(p,f) = CrownFireIntensity(HPA, fuel(f)%canopy_fuel_load, &
                patch%area, crown_frac_burnt, ROS_final(p,f))
            else 
                ROS_final(p,f) = ROS_front(p,f)
                FI_final(p,f) = FI(p,f)
            end if

        end do

    end do

  ! write out data
  call WriteCanopyFuelData(out_file, num_fuel_models, num_patch_types, CBD, CBH, &
  canopy_fuel_load, ROS_front, FI, FI_init, ROS_actfm10, ROS_critical, CFB,      &
  ROS_final, FI_final, fuel_models, patch_ids)


end program FatesTestCanopyFuel