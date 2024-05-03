module FatesTestFireMod
  !
  ! DESCRIPTION:
  !		Module to support testing the FATES SPIFTIRE model
  !

  use FatesConstantsMod, only : r8 => fates_r8 
  use EDTypesMod,        only : ed_site_type
  use FatesPatchMod,     only : fates_patch_type
  use SFNesterovMod,     only : nesterov_index

  implicit none
  private

  public :: SetUpSite

  contains 

    !=====================================================================================

    subroutine SetUpSite(site)
      !
      ! DESCRIPTION:
      ! Sets up site, patch, litter, and fuel
      ! This only sets up the stuff we actually need for these subroutines
      !
    
      ! ARGUMENTS:
      type(ed_site_type), target :: site ! site object

      ! LOCALS:
      type(fates_patch_type), pointer :: patch
      
      ! set up fire weather class
      allocate(nesterov_index :: site%fireWeather)
      call site%fireWeather%Init()

      ! set up one patch
      allocate(patch)
      call patch%Init(2, 1)

      patch%patchno = 1
      patch%younger => null()
      patch%older   => null()

      site%youngest_patch => patch
      site%oldest_patch   => patch

    end subroutine SetUpSite

    !=====================================================================================

    ! subroutine ReadDatmData(nc_file, temp_degC, precip, rh, wind)
    !   !
    !   ! DESCRIPTION:
    !   ! Reads and returns DATM data
    !   !
    
    !   ! ARGUMENTS:
    !   character(len=*),      intent(in)  :: nc_file      ! netcdf file with DATM data
    !   real(r8), allocatable, intent(out) :: temp_degC(:) ! daily air temperature [degC]
    !   real(r8), allocatable, intent(out) :: precip(:)    ! daily precipitation [mm]
    !   real(r8), allocatable, intent(out) :: rh(:)        ! daily relative humidity [%]
    !   real(r8), allocatable, intent(out) :: wind(:)      ! daily wind speed [m/s]

    !   ! LOCALS:
    !   integer :: ncid ! netcdf file unit number

    !   ! open file
    !   call OpenNCFile(trim(nc_file), ncid, 'read')
      
    !   ! read in data
    !   call GetVar(ncid, 'temp_degC', temp_degC)
    !   call GetVar(ncid, 'precip', precip)
    !   call GetVar(ncid, 'RH', rh)
    !   call GetVar(ncid, 'wind', wind)

    !   ! close file
    !   call CloseNCFile(ncid)

    ! end subroutine ReadDatmData

    !=====================================================================================

    ! subroutine WriteFireData(out_file, nsteps, time_counter, temp_degC, precip, rh, NI,  &
    !   loading, moisture, av_moisture, ros)
    !   !
    !   ! DESCRIPTION:
    !   ! writes out data from the unit test
    !   !
    
    !   ! ARGUMENTS:
    !   character(len=*), intent(in) :: out_file
    !   integer,          intent(in) :: nsteps
    !   integer,          intent(in) :: time_counter(:)
    !   real(r8),         intent(in) :: temp_degC(:)
    !   real(r8),         intent(in) :: precip(:)
    !   real(r8),         intent(in) :: rh(:)
    !   real(r8),         intent(in) :: NI(:)
    !   real(r8),         intent(in) :: loading(:)
    !   real(r8),         intent(in) :: moisture(:,:)
    !   real(r8),         intent(in) :: av_moisture(:)
    !   real(r8),         intent(in) :: ros(:)

    !   ! LOCALS:
    !   integer          :: ncid         ! netcdf id
    !   character(len=8) :: dim_names(2) ! dimension names
    !   integer          :: dimIDs(2)    ! dimension IDs
    !   integer          :: timeID, litterID
    !   integer          :: tempID, precipID, rhID, NIID, loadingID, moistureID
    !   integer          :: av_moistureID, rosID

    !   ! dimension names
    !   dim_names = [character(len=12) :: 'time', 'litter_class']

    !   ! open file
    !   call OpenNCFile(trim(out_file), ncid, 'readwrite')

    !   ! register dimensions
    !   call RegisterNCDims(ncid, dim_names, (/nsteps, nfsc/), 2, dimIDs)

    !   ! register time
    !   call RegisterVar1D(ncid, 'time', dimIDs(1), type_int,                             &
    !     [character(len=20)  :: 'time_origin', 'units', 'calendar', 'long_name'],        &
    !     [character(len=150) :: '2018-01-01 00:00:00', 'days since 2018-01-01 00:00:00', &
    !       'gregorian', 'time'],                                                         &
    !     4, timeID)

    !   ! register litter class
    !   call RegisterVar1D(ncid, 'litter_class', dimIDs(2), type_int,       &
    !     [character(len=20)  :: 'units', 'long_name'],                     &
    !     [character(len=150) :: '', 'litter class'], 2, litterID)

    !   ! register temperature
    !   call RegisterVar1D(ncid, 'temp_degC', dimIDs(1), type_double,      &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'],     &
    !     [character(len=150) :: 'time', 'degrees C', 'air temperature'],  &                                                  
    !     3, tempID)

    !   ! register precipitation
    !   call RegisterVar1D(ncid, 'precip', dimIDs(1), type_double,      &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'],  &
    !     [character(len=150) :: 'time', 'mm', 'precipitation'],        &                                                  
    !     3, precipID)

    !   ! register relative humidity
    !   call RegisterVar1D(ncid, 'RH', dimIDs(1), type_double,         &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    !     [character(len=150) :: 'time', '%', 'relative humidity'],    &                                                  
    !     3, rhID)

    !   ! register Nesterov Index
    !   call RegisterVar1D(ncid, 'NI', dimIDs(1), type_double,         &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    !     [character(len=150) :: 'time', '', 'Nesterov Index'],        &                                                  
    !     3, NIID)

    !   ! register loading
    !   call RegisterVar1D(ncid, 'loading', dimIDs(1), type_double,    &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    !     [character(len=150) :: 'time', 'kg m-2', 'fuel loading'],    &                                                  
    !     3, loadingID)

    !   ! register moisture
    !   call RegisterVar2D(ncid, 'fuel_moisture', dimIDs(1:2), type_double,       &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    !     [character(len=150) :: 'time litter_class', 'm3 m-3', 'fuel moisture'], &                                                  
    !     3, moistureID)

    !   ! register average moisture
    !   call RegisterVar1D(ncid, 'average_moisture', dimIDs(1), type_double,    &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    !     [character(len=150) :: 'time', 'm3 m-3', 'average fuel moisture'],    &                                                  
    !     3, av_moistureID)

    !   ! register ROS
    !   call RegisterVar1D(ncid, 'ros', dimIDs(1), type_double,    &
    !     [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    !     [character(len=150) :: 'time', 'm min-1', 'rate of forward spread'],    &                                                  
    !     3, rosID)

    !   call Check(NF90_ENDDEF(ncid))

    !   call Check(NF90_PUT_VAR(ncid, timeID, time_counter))
    !   call Check(NF90_PUT_VAR(ncid, litterID, (/tw_sf, sb_sf, lb_sf, tr_sf, dl_sf, lg_sf/)))
    !   call Check(NF90_PUT_VAR(ncid, tempID, temp_degC(:)))
    !   call Check(NF90_PUT_VAR(ncid, precipID, precip(:)))
    !   call Check(NF90_PUT_VAR(ncid, rhID, rh(:)))
    !   call Check(NF90_PUT_VAR(ncid, NIID, NI(:)))
    !   call Check(NF90_PUT_VAR(ncid, loadingID, loading(:)))
    !   call Check(NF90_PUT_VAR(ncid, moistureID, moisture(:,:)))
    !   call Check(NF90_PUT_VAR(ncid, av_moistureID, av_moisture(:)))
    !   call Check(NF90_PUT_VAR(ncid, rosID, ros(:)))
 
    !   call CloseNCFile(ncid)

    ! end subroutine WriteFireData

end module FatesTestFireMod