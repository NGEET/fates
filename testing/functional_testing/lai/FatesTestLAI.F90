program FatesTestAllometry

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesAllometryMod,           only : tree_lai
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use EDPftvarcon,                 only : EDPftvarcon_inst

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader ! param reader instance
  character(len=:), allocatable      :: param_file   ! input parameter file
  real(r8),         allocatable      :: lai(:,:)     ! lai array [m2/m2]
  real(r8),         allocatable      :: leafc(:)     ! leaf carbon array [kgC]
  integer                            :: numpft       ! number of pfts (from parameter file)
  integer                            :: numleafc     ! size of leafC array
  integer                            :: i, j         ! looping indices
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'lai_out.nc' ! output file
  real(r8),         parameter :: min_leafc = 0.0_r8  ! minimum leaf C to calculate [kgC]
  real(r8),         parameter :: max_leafc = 15.0_r8 ! maximum leaf C to calculate [kgC]
  real(r8),         parameter :: leafc_inc = 0.1_r8  ! leaf C increment to use [kgC]

  interface

    subroutine WriteLAIData(out_file, numleafc, numpft, lai, leaf_c)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file             
      integer,          intent(in) :: numleafc             
      integer,          intent(in) :: numpft                   
      real(r8),         intent(in) :: lai(:)                 
      real(r8),         intent(in) :: leaf_c(:,:)              

    end subroutine WriteLAIData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()

  ! determine sizes of arrays
  numpft = size(prt_params%wood_density, dim=1)
  numleafc = int((max_leafc - min_leafc)/leafc_inc + 1)

  ! allocate arrays
  allocate(leafc(numleafc))
  allocate(lai(numleafc,numdbh))
  
  ! initialize leaf c array
  do i = 1, numleafc
    leafc(i) = min_leafc + leafc_inc*(i-1)
  end do
  
  ! calculate allometries
  do i = 1, numpft
    do j = 1, numleafc
      lai(i,j) = tree_lai(leaf_c(j), pft, 3.5_r8, 1.0_r8, 1, 0.0_r8,                     &
        EDPftvarcon_inst%vcmax25top(i, 1))
    end do
  end do

  ! write out data to netcdf file
  call WriteLAIData(out_file, numleafc, numpft, lai, leaf_c)

end program FatesTestAllometry

! ----------------------------------------------------------------------------------------

subroutine WriteAllometryData(out_file, numleafc, numpft, lai, leaf_c)
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
  character(len=*), intent(in) :: out_file                 ! output file name
  integer,          intent(in) :: numleafc                 ! size of leafc array
  integer,          intent(in) :: numpft                   ! number of pfts
  real(r8),         intent(in) :: lai(:)                   ! leaf area index [m2/m2]
  real(r8),         intent(in) :: leaf_c(:,:)              ! leaf carbon [kg]

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=8)     :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: leafcID, pftID   ! variable IDs for dimensions
  integer              :: laiID
  ! create pft indices
  allocate(pft_indices(numpft))
  do i = 1, numpft
    pft_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'pft', 'leafc']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/numpft, numleafc/), 2, dimIDs)

  ! register pft
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
  
  ! register leaf C
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'kgC', 'leaf carbon'], 2, leafcID)

  ! register lai
  call RegisterVar(ncid, 'lai', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft leafc', 'm2/m2', 'leaf area index'],      &
    3, laiID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, leafcID, leaf_c(:))
  call WriteVar(ncid, laiID, lai(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteAllometryData