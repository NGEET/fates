module FatesTestFireMod
  !
  ! DESCRIPTION:
  !		Module to support testing the FATES SPIFTIRE model
  !

  use FatesConstantsMod,   only : r8 => fates_r8 
  use EDTypesMod,          only : ed_site_type
  use FatesPatchMod,       only : fates_patch_type
  use SFNesterovMod,       only : nesterov_index
  use FatesUnitTestIOMod,  only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims
  use FatesUnitTestIOMod,  only : RegisterVar, EndNCDef, WriteVar
  use FatesUnitTestIOMod,  only : type_double, type_int, type_char
  use FatesFuelClassesMod, only : num_fuel_classes
  use SyntheticFuelModels, only : fuel_models_array_class
  use SFParamsMod,         only : SF_val_CWD_frac
  use FatesFuelMod,        only : fuel_type

  implicit none
  private

  public :: SetUpFuel, ReadDatmData, WriteFireData

  contains 

    !=====================================================================================

    subroutine SetUpFuel(fuel, fuel_model_array, fuel_model_index, fuel_name, fuel_carrier)
      !
      ! DESCRIPTION:
      ! Sets up fuel loading
      !
    
      ! ARGUMENTS:
      type(fuel_type),               intent(inout) :: fuel             ! fuel object
      type(fuel_models_array_class), intent(in)    :: fuel_model_array ! array of fuel models
      integer,                       intent(in)    :: fuel_model_index ! fuel model index
      character(len=100),            intent(out)   :: fuel_name        ! name of fuel model
      character(len=2),              intent(out)   :: fuel_carrier     ! fuel carrier for fuel model
      
      ! LOCALS:
      integer  :: i                   ! position of fuel model in array
      real(r8) :: leaf_litter         ! leaf litter [kg/m2] 
      real(r8) :: twig_litter         ! twig litter [kg/m2]
      real(r8) :: small_branch_litter ! small branch litter [kg/m2]
      real(r8) :: large_branch_litter ! large branch litter [kg/m2]
      real(r8) :: grass_litter        ! grass litter [kg/m2]
      

      ! get fuel model position in array
      i = fuel_model_array%FuelModelPosition(fuel_model_index)
      
      ! fuel model data
      leaf_litter = fuel_model_array%fuel_models(i)%hr1_loading
      twig_litter = fuel_model_array%fuel_models(i)%hr10_loading
      
      ! small vs. large branches based on input parameter file
      small_branch_litter = fuel_model_array%fuel_models(i)%hr100_loading*SF_val_CWD_frac(2)/ &
        (SF_val_CWD_frac(2) + SF_val_CWD_frac(3))
      large_branch_litter = fuel_model_array%fuel_models(i)%hr100_loading*SF_val_CWD_frac(3)/ &
        (SF_val_CWD_frac(2) + SF_val_CWD_frac(3))
        
      grass_litter = fuel_model_array%fuel_models(i)%live_herb_loading
      
      fuel_name = fuel_model_array%fuel_models(i)%fuel_model_name
      fuel_carrier = fuel_model_array%fuel_models(i)%carrier
      
      call fuel%UpdateLoading(leaf_litter, twig_litter, small_branch_litter,    &
        large_branch_litter, 0.0_r8, grass_litter)
      
    end subroutine SetUpFuel

    !=====================================================================================
    
    subroutine ReadDatmData(nc_file, temp_degC, precip, rh, wind)
      !
      ! DESCRIPTION:
      ! Reads and returns DATM data
      !
    
      ! ARGUMENTS:
      character(len=*),      intent(in)  :: nc_file      ! netcdf file with DATM data
      real(r8), allocatable, intent(out) :: temp_degC(:) ! daily air temperature [degC]
      real(r8), allocatable, intent(out) :: precip(:)    ! daily precipitation [mm]
      real(r8), allocatable, intent(out) :: rh(:)        ! daily relative humidity [%]
      real(r8), allocatable, intent(out) :: wind(:)      ! daily wind speed [m/s]

      ! LOCALS:
      integer :: ncid ! netcdf file unit number

      ! open file
      call OpenNCFile(trim(nc_file), ncid, 'read')
      
      ! read in data
      call GetVar(ncid, 'temp_degC', temp_degC)
      call GetVar(ncid, 'precip', precip)
      call GetVar(ncid, 'RH', rh)
      call GetVar(ncid, 'wind', wind)

      ! close file
      call CloseNCFile(ncid)

    end subroutine ReadDatmData

    !=====================================================================================

    subroutine WriteFireData(out_file, nsteps, nfuelmods, temp_degC, precip, rh, NI,     &
      loading, frac_loading, fuel_BD, fuel_SAV, non_trunk_loading, fuel_moisture,            &
      fuel_MEF, fuel_models, carriers)
      !
      ! DESCRIPTION:
      ! writes out data from the unit test
      !
    
      ! ARGUMENTS:
      character(len=*),   intent(in) :: out_file
      integer,            intent(in) :: nsteps
      integer,            intent(in) :: nfuelmods
      real(r8),           intent(in) :: temp_degC(:)
      real(r8),           intent(in) :: precip(:)
      real(r8),           intent(in) :: rh(:)
      real(r8),           intent(in) :: NI(:)
      real(r8),           intent(in) :: loading(:,:)
      real(r8),           intent(in) :: frac_loading(:,:)
      real(r8),           intent(in) :: non_trunk_loading(:)
      real(r8),           intent(in) :: fuel_moisture(:,:)
      real(r8),           intent(in) :: fuel_MEF(:,:)
      real(r8),           intent(in) :: fuel_BD(:)
      real(r8),           intent(in) :: fuel_SAV(:)
      integer,            intent(in) :: fuel_models(:)
      character(len=2),   intent(in) :: carriers(:)
      
      ! LOCALS:
      integer, allocatable :: time_index(:) ! array of time index
      integer              :: ncid         ! netcdf id
      integer              :: i            ! looping index
      character(len=20)    :: dim_names(3) ! dimension names
      integer              :: dimIDs(3)    ! dimension IDs
      ! variable IDS
      integer              :: timeID, litterID
      integer              :: modID
      integer              :: tempID, precipID
      integer              :: rhID, NIID, loadingID
      integer              :: frac_loadingID
      integer              :: tot_loadingID
      integer              :: BDID, SAVID
      integer              :: moistID
      integer              :: cID, mefID
      
      ! create pft indices
      allocate(time_index(nsteps))
      do i = 1, nsteps
        time_index(i) = i
      end do
      
      ! dimension names
      dim_names = [character(len=20) :: 'time', 'litter_class', 'fuel_model']

      ! open file
      call OpenNCFile(trim(out_file), ncid, 'readwrite')

      ! register dimensions
      call RegisterNCDims(ncid, dim_names, (/nsteps, num_fuel_classes, nfuelmods/), 3, dimIDs)

      ! first register dimension variables
      
      ! register time
      call RegisterVar(ncid, 'time', dimIDs(1:1), type_int,                             &
        [character(len=20)  :: 'time_origin', 'units', 'calendar', 'long_name'],        &
        [character(len=150) :: '2018-01-01 00:00:00', 'days since 2018-01-01 00:00:00', &
          'gregorian', 'time'],                                                         &
        4, timeID)

      ! register litter class
      call RegisterVar(ncid, 'litter_class', dimIDs(2:2), type_int,       &
        [character(len=20)  :: 'units', 'long_name'],                     &
        [character(len=150) :: '', 'fuel class'], 2, litterID)
        
      ! register fuel models
      call RegisterVar(ncid, 'fuel_model', dimIDs(3:3), type_int,     &
        [character(len=20)  :: 'units', 'long_name'],                 &
        [character(len=150) :: '', 'fuel model index'], 2, modID)

      ! then register actual variables
    
      ! register fuel carriers
      call RegisterVar(ncid, 'carrier', dimIDs(3:3), type_char,            &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
        [character(len=150) :: 'fuel_model_index', '', 'carrier of fuel'],  &
        3, cID)
        
      ! register temperature
      call RegisterVar(ncid, 'temp_degC', dimIDs(1:1), type_double,      &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],     &
        [character(len=150) :: 'time', 'degrees C', 'air temperature'],  &                                                  
        3, tempID)

      ! register precipitation
      call RegisterVar(ncid, 'precip', dimIDs(1:1), type_double,      &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],  &
        [character(len=150) :: 'time', 'mm', 'precipitation'],        &                                                  
        3, precipID)

      ! register relative humidity
      call RegisterVar(ncid, 'RH', dimIDs(1:1), type_double,         &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
        [character(len=150) :: 'time', '%', 'relative humidity'],    &                                                  
        3, rhID)

      ! register Nesterov Index
      call RegisterVar(ncid, 'NI', dimIDs(1:1), type_double,         &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
        [character(len=150) :: 'time', '', 'Nesterov Index'],        &                                                  
        3, NIID)
        
      ! register fuel moisture
      call RegisterVar(ncid, 'fuel_moisture', (/dimIDs(1), dimIDs(3)/), type_double,   &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],                   &
        [character(len=150) :: 'time fuel_model', 'm3 m-3', 'average fuel moisture'],  &                                                  
        3, moistID)
        
      ! register fuel MEF
      call RegisterVar(ncid, 'fuel_MEF', (/dimIDs(1), dimIDs(3)/), type_double,   &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],                   &
        [character(len=150) :: 'time fuel_model', 'm3 m-3', 'average fuel moisture of extinction'],  &                                                  
        3, mefID)
        
        
      ! register fuel loading
      call RegisterVar(ncid, 'fuel_loading', dimIDs(2:3), type_double,                 &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],                   &
        [character(len=150) :: 'litter_class fuel_model', 'kgC m-2', 'fuel loading'],  &                                                  
        3, loadingID)
        
      ! register fractional fuel loading
      call RegisterVar(ncid, 'frac_loading', dimIDs(2:3), type_double,                 &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],                   &
        [character(len=150) :: 'litter_class fuel_model', '', 'fractional loading'],   &                                                    
        3, frac_loadingID)
        
      ! register non-trunk fuel loading
      call RegisterVar(ncid, 'non_trunk_loading', dimIDs(3:3), type_double,    &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
        [character(len=150) :: 'fuel_model', 'kgC m-2', 'total loading'],  &                                                  
        3, tot_loadingID)
        
      ! register fuel bulk density
        call RegisterVar(ncid, 'bulk_density', dimIDs(3:3), type_double,      &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],          &
        [character(len=150) :: 'fuel_model', 'kg m-3', 'fuel bulk density'],  &                                                  
        3, BDID)
        
      ! register fuel SAV
        call RegisterVar(ncid, 'SAV', dimIDs(3:3), type_double,                             &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],                        &
        [character(len=150) :: 'fuel_model', 'cm-1', 'fuel surface area to volume ratio'],  &                                                  
        3, SAVID)
        
      ! finish defining variables
      call EndNCDef(ncid)

      ! write out data
      call WriteVar(ncid, timeID, time_index)
      call WriteVar(ncid, litterID, (/1, 2, 3, 4, 5, 6/))
      call WriteVar(ncid, modID, fuel_models(:))
      call WriteVar(ncid, cID, carriers(:))
      call WriteVar(ncid, tempID, temp_degC(:))
      call WriteVar(ncid, precipID, precip(:))
      call WriteVar(ncid, rhID, rh(:))
      call WriteVar(ncid, NIID, NI(:))
      call WriteVar(ncid, loadingID, loading(:,:))
      call WriteVar(ncid, frac_loadingID, frac_loading(:,:))
      call WriteVar(ncid, tot_loadingID, non_trunk_loading(:))
      call WriteVar(ncid, moistiD, fuel_moisture(:,:))
      call WriteVar(ncid, BDID, fuel_BD(:))
      call WriteVar(ncid, SAVID, fuel_SAV(:))
      call WriteVar(ncid, mefID, fuel_MEF(:,:))
 
      call CloseNCFile(ncid)

    end subroutine WriteFireData

end module FatesTestFireMod
