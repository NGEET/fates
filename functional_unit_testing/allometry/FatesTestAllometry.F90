program FatesTestAllometry

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesAllometryMod,           only : h_allom
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader
  character(len=*), parameter        :: param_file = 'fates_params_default.nc'
  character(len=*), parameter        :: out_file = 'allometry_out.nc'
  integer                            :: numpft
  integer                            :: i, j
  integer                            :: numdbh
  real(r8), allocatable              :: dbh(:)       ! diameter at breast height [cm]
  real(r8), allocatable              :: height(:, :) ! height [m]

  ! CONSTANTS:
  real(r8) :: min_dbh = 0.5_r8   ! minimum DBH to calculate [cm]
  real(r8) :: max_dbh = 200.0_r8 ! maximum DBH to calculate [cm]
  real(r8) :: dbh_inc = 0.5_r8   ! DBH increment to use [cm]

  interface

    subroutine WriteAllometryData(out_file, ndbh, numpft, dbh, height)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : RegisterVar1D, WriteVar, RegisterVar2D
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: ndbh, numpft
      real(r8),         intent(in) :: dbh(:)
      real(r8),         intent(in) :: height(:,:)

    end subroutine WriteAllometryData

  end interface
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  numpft = size(prt_params%wood_density, dim=1)

  ! allocate arrays and initialize DBH array
  numdbh = int((max_dbh - min_dbh)/dbh_inc + 1)
  
  allocate(dbh(numdbh))
  allocate(height(numdbh, numpft))
  
  do i = 1, numdbh
    dbh(i) = min_dbh + dbh_inc*(i-1)
  end do

  ! calculate allometries
  do i = 1, numpft
    do j = 1, numdbh
      call h_allom(dbh(j), i, height(j, i))
    end do
  end do

  call WriteAllometryData(out_file, numdbh, numpft, dbh, height)
  
end program FatesTestAllometry

! ----------------------------------------------------------------------------------------

subroutine WriteAllometryData(out_file, numdbh, numpft, dbh, height)
  !
  ! DESCRIPTION:
  ! Writes out data from the allometry test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : RegisterVar1D, WriteVar, RegisterVar2D
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file
  integer,          intent(in) :: numdbh, numpft
  real(r8),         intent(in) :: dbh(:)
  real(r8),         intent(in) :: height(:,:)

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=8)     :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: dbhID, pftID   ! variable IDs for dimensions
  integer              :: heightID

  ! create pft indices
  allocate(pft_indices(numpft))
  do i = 1, numpft
    pft_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'dbh', 'pft']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/numdbh, numpft/), 2, dimIDs)

  ! register dbh
  call RegisterVar1D(ncid, dim_names(1), dimIDs(1), type_double,           &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'cm', 'diameter at breast height'], 2, dbhID)

  ! register pft
  call RegisterVar1D(ncid, dim_names(2), dimIDs(2), type_int,           &
    [character(len=20)  :: 'units', 'long_name'],                     &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)

  ! register height
  call RegisterVar2D(ncid, 'height', dimIDs(1:2), type_double,     &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft dbh', 'm', 'plant height'],      &                                                  
    3, heightID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, dbhID, dbh(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, heightID, height(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteAllometryData