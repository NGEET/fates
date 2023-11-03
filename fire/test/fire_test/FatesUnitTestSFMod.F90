module FatesUnitTestSFMod
  !
  ! DESCRIPTION:
  !		Module to support testing the FATES SPIFTIRE model
  !

  use FatesConstantsMod,  only : r8 => fates_r8 
  use FatesUnitTestIOMod, only : OpenNCFile, CloseNCFile, GetVar, RegisterNCDims
  use FatesUnitTestIOMod, only : RegisterVar1D, type_int, type_double, Check
  use netcdf

  implicit none
  private

  public :: ReadDatmData, WriteFireData

  contains 

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

    subroutine WriteFireData(out_file, nsteps, time_counter, temp_degC, precip, rh, NI)
      !
      ! DESCRIPTION:
      ! writes out data from the unit test
      !
    
      ! ARGUMENTS:
      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: nsteps
      integer,          intent(in) :: time_counter(:)
      real(r8),         intent(in) :: temp_degC(:)
      real(r8),         intent(in) :: precip(:)
      real(r8),         intent(in) :: rh(:)
      real(r8),         intent(in) :: NI(:)

      ! LOCALS:
      integer          :: ncid         ! netcdf id
      character(len=8) :: dim_names(1) ! dimension names
      integer          :: dimIDs(1)    ! dimension IDs
      integer          :: timeID
      integer          :: tempID, precipID, rhID, NIID

      ! dimension names
      dim_names = [character(len=8) :: 'time']

      ! open file
      call OpenNCFile(trim(out_file), ncid, 'readwrite')

      ! register dimensions
      call RegisterNCDims(ncid, dim_names, (/nsteps/), 1, dimIDs)

      ! register time
      call RegisterVar1D(ncid, 'time', dimIDs(1), type_int,                             &
        [character(len=20)  :: 'time_origin', 'units', 'calendar', 'long_name'],        &
        [character(len=150) :: '2018-01-01 00:00:00', 'days since 2018-01-01 00:00:00', &
          'gregorian', 'time'],                                                         &
        4, timeID)

      ! register temperature
      call RegisterVar1D(ncid, 'temp_degC', dimIDs(1), type_double,      &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],     &
        [character(len=150) :: 'time', 'degrees C', 'air temperature'],  &                                                  
        3, tempID)

      ! register precipitation
      call RegisterVar1D(ncid, 'precip', dimIDs(1), type_double,      &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'],  &
        [character(len=150) :: 'time', 'mm', 'precipitation'],        &                                                  
        3, precipID)

      ! register relative humidity
      call RegisterVar1D(ncid, 'RH', dimIDs(1), type_double,       &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
        [character(len=150) :: 'time', '%', 'relative humidity'],    &                                                  
        3, rhID)

      ! register Nesterov Index
      call RegisterVar1D(ncid, 'NI', dimIDs(1), type_double,       &
        [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
        [character(len=150) :: 'time', '', 'Nesterov Index'],    &                                                  
        3, NIID)

      call Check(NF90_ENDDEF(ncid))

      call Check(NF90_PUT_VAR(ncid, timeID, time_counter))
      call Check(NF90_PUT_VAR(ncid, tempID, temp_degC(:)))
      call Check(NF90_PUT_VAR(ncid, precipID, precip(:)))
      call Check(NF90_PUT_VAR(ncid, rhID, rh(:)))
      call Check(NF90_PUT_VAR(ncid, NIID, NI(:)))

      call CloseNCFile(ncid)

    end subroutine WriteFireData

end module FatesUnitTestSFMod