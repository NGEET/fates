program FatesTestAllometry

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesAllometryMod,           only : tree_lai, blmax_allom, carea_allom
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use EDPftvarcon,                 only : EDPftvarcon_inst

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader    ! param reader instance
  character(len=:), allocatable      :: param_file      ! input parameter file
  real(r8),         allocatable      :: dbh(:)          ! dbh [cm]
  real(r8),         allocatable      :: lai(:,:)        ! lai array [m2/m2]
  real(r8),         allocatable      :: leafc(:,:)      ! leaf carbon array [kgC]
  real(r8),         allocatable      :: crown_area(:,:) ! crown area per cohort [m2]
  real(r8)                           :: can_lai(10)     ! canopy lai [m2/m2]
  integer                            :: numpft          ! number of pfts (from parameter file)
  integer                            :: numdbh          ! size of dbh array
  integer                            :: i, j            ! looping indices
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'lai_out.nc' ! output file
  real(r8),         parameter :: min_dbh = 0.5_r8    ! minimum DBH to calculate [cm]
  real(r8),         parameter :: max_dbh = 40.0_r8  ! maximum DBH to calculate [cm]
  real(r8),         parameter :: dbh_inc = 0.5_r8    ! DBH increment to use [cm]

  interface

    subroutine WriteLAIData(out_file, numdbh, numpft, dbh, lai, leaf_c, crown_area)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file             
      integer,          intent(in) :: numdbh             
      integer,          intent(in) :: numpft
      real(r8),         intent(in) :: dbh(:)                
      real(r8),         intent(in) :: lai(:,:)                 
      real(r8),         intent(in) :: leaf_c(:,:)
      real(r8),         intent(in) :: crown_area(:,:)          

    end subroutine WriteLAIData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()

  ! determine sizes of arrays
  numpft = size(prt_params%wood_density, dim=1)
  numdbh = int((max_dbh - min_dbh)/dbh_inc + 1)

  ! allocate arrays
  allocate(dbh(numdbh))
  allocate(leafc(numpft, numdbh))
  allocate(lai(numpft, numdbh))
  allocate(crown_area(numpft, numdbh))
  
  ! initialize this to zero for now
  do i = 1, 10
    can_lai(i) = 0.0_r8
  end do 
    
  do i = 1, numpft
    do j = 1, numdbh
      
      ! initialize dbh
      dbh(j) = min_dbh + dbh_inc*(j-1)
      
      ! calculate leaf c
      call blmax_allom(dbh(j), i, leafc(i,j))
      
      ! calculate crown area
      call carea_allom(dbh(j), 1.0_r8, 1.0_r8, i, 1, crown_area(i,j))
      
      lai(i,j) = tree_lai(leafc(i,j), i, crown_area(i,j), 1.0_r8, 1, can_lai,            &
        EDPftvarcon_inst%vcmax25top(i, 1))
    end do
  end do
  
  ! write out data to netcdf file
  call WriteLAIData(out_file, numdbh, numpft, dbh, lai, leafc, crown_area)

end program FatesTestAllometry

! ----------------------------------------------------------------------------------------

subroutine WriteLAIData(out_file, numdbh, numpft, dbh, lai, leaf_c, crown_area)
  !
  ! DESCRIPTION:
  ! Writes out data from the allometry test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file        ! output file name
  integer,          intent(in) :: numdbh          ! size of dbh array
  integer,          intent(in) :: numpft          ! number of pfts
  real(r8),         intent(in) :: dbh(:)          ! diameter at breast height [cm]
  real(r8),         intent(in) :: lai(:,:)        ! leaf area index [m2/m2]
  real(r8),         intent(in) :: leaf_c(:,:)     ! leaf carbon [kg]
  real(r8),         intent(in) :: crown_area(:,:) ! crown area [m2]

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=8)     :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: dbhID, pftID   ! variable IDs for dimensions
  integer              :: laiID, leafcID
  integer              :: crownID
  
  ! create pft indices
  allocate(pft_indices(numpft))
  do i = 1, numpft
    pft_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'pft', 'dbh']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/numpft, numdbh/), 2, dimIDs)

  ! register pft
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
  
  ! register dbh
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'cm', 'diameter at breast height'], 2, dbhID)

  ! register lai
  call RegisterVar(ncid, 'lai', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft dbh', 'm2/m2', 'leaf area index'],      &
    3, laiID)
    
  ! register leafc
  call RegisterVar(ncid, 'leafc', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft dbh', 'kgC', 'leaf carbon'],      &
    3, leafcID)
    
  ! register crown area
  call RegisterVar(ncid, 'crown_area', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft dbh', 'm2', 'crown area'],      &
    3, crownID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, dbhID, dbh(:))
  call WriteVar(ncid, leafcID, leaf_c(:,:))
  call WriteVar(ncid, laiID, lai(:,:))
  call WriteVar(ncid, crownID, crown_area(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteLAIData