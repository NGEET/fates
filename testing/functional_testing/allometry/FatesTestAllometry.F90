program FatesTestAllometry

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesAllometryMod,           only : h_allom, bagw_allom, blmax_allom
  use FatesAllometryMod,           only : carea_allom, bsap_allom, bbgw_allom
  use FatesAllometryMod,           only : bfineroot, bstore_allom, bdead_allom
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader            ! param reader instance
  character(len=:), allocatable      :: param_file              ! input parameter file
  integer                            :: numpft                  ! number of pfts (from parameter file)
  integer                            :: i, j                    ! looping indices
  integer                            :: numdbh                  ! size of dbh array
  real(r8), allocatable              :: dbh(:)                  ! diameter at breast height [cm]
  real(r8), allocatable              :: height(:, :)            ! height [m]
  real(r8), allocatable              :: bagw(:, :)              ! aboveground woody biomass [kgC]
  real(r8), allocatable              :: blmax(:, :)             ! plant leaf biomass [kgC]
  real(r8), allocatable              :: crown_area(:, :)        ! crown area per cohort [m2]
  real(r8), allocatable              :: sapwood_area(:, :)      ! cross sectional area of sapwood at reference height [m2]
  real(r8), allocatable              :: bsap(:, :)              ! sapwood biomass [kgC]
  real(r8), allocatable              :: bbgw(:, :)              ! belowground woody biomass [kgC]
  real(r8), allocatable              :: fineroot_biomass(:, :)  ! belowground fineroot biomass [kgC]
  real(r8), allocatable              :: bstore(:, :)            ! allometric target storage biomass [kgC]
  real(r8), allocatable              :: bdead(:, :)             ! structural biomass (heartwood/struct) [kgC]
  real(r8), allocatable              :: total_biom_tissues(:,:) ! total biomass calculated as bleaf + bfineroot + bdead + bsap [kgC]
  real(r8), allocatable              :: total_biom_parts(:,:)   ! total biomass calculated as bleaf + bfineroot + agbw + bgbw [kgC]

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'allometry_out.nc'    ! output file
  real(r8),         parameter :: min_dbh = 0.5_r8                 ! minimum DBH to calculate [cm]
  real(r8),         parameter :: max_dbh = 200.0_r8               ! maximum DBH to calculate [cm]
  real(r8),         parameter :: dbh_inc = 0.5_r8                 ! DBH increment to use [cm]

  integer,          parameter :: crown_damage = 1                 ! crown damage
  real(r8),         parameter :: elongation_factor = 1.0_r8       ! elongation factor for stem
  real(r8),         parameter :: elongation_factor_roots = 1.0_r8 ! elongation factor for roots
  real(r8),         parameter :: site_spread = 1.0_r8             ! site spread
  real(r8),         parameter :: canopy_trim = 1.0_r8             ! canopy trim
  real(r8),         parameter :: nplant = 1.0_r8                  ! number of plants per cohort
  real(r8),         parameter :: leaf_to_fineroot = 1.0_r8        ! leaf to fineroot ratio

  interface

    subroutine WriteAllometryData(out_file, ndbh, numpft, dbh, height, bagw, blmax,      &
        crown_area, sapwood_area, bsap, bbgw, fineroot_biomass, bstore, bdead,           &
        total_biom_parts, total_biom_tissues)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: ndbh, numpft
      real(r8),         intent(in) :: dbh(:)
      real(r8),         intent(in) :: height(:,:)
      real(r8),         intent(in) :: bagw(:,:)
      real(r8),         intent(in) :: blmax(:, :)
      real(r8),         intent(in) :: crown_area(:, :)
      real(r8),         intent(in) :: sapwood_area(:, :)
      real(r8),         intent(in) :: bsap(:, :)
      real(r8),         intent(in) :: bbgw(:, :)
      real(r8),         intent(in) :: fineroot_biomass(:, :)
      real(r8),         intent(in) :: bstore(:, :)
      real(r8),         intent(in) :: bdead(:, :)
      real(r8),         intent(in) :: total_biom_parts(:, :)
      real(r8),         intent(in) :: total_biom_tissues(:, :)
    end subroutine WriteAllometryData

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
  allocate(height(numdbh, numpft))
  allocate(bagw(numdbh, numpft))
  allocate(blmax(numdbh, numpft))
  allocate(crown_area(numdbh, numpft))
  allocate(sapwood_area(numdbh, numpft))
  allocate(bsap(numdbh, numpft))
  allocate(bbgw(numdbh, numpft))
  allocate(fineroot_biomass(numdbh, numpft))
  allocate(bstore(numdbh, numpft))
  allocate(bdead(numdbh, numpft))
  allocate(total_biom_parts(numdbh, numpft))
  allocate(total_biom_tissues(numdbh, numpft))

  ! initialize dbh array
  do i = 1, numdbh
    dbh(i) = min_dbh + dbh_inc*(i-1)
  end do

  ! calculate allometries
  do i = 1, numpft
    do j = 1, numdbh
      call h_allom(dbh(j), i, height(j, i))
      call bagw_allom(dbh(j), i, crown_damage, elongation_factor, bagw(j, i))
      call blmax_allom(dbh(j), i, blmax(j, i))
      call carea_allom(dbh(j), nplant, site_spread, i, crown_damage, crown_area(j, i))
      call bsap_allom(dbh(j), i, crown_damage, canopy_trim, elongation_factor,           &
        sapwood_area(j, i), bsap(j, i))
      call bbgw_allom(dbh(j), i, elongation_factor, bbgw(j, i))
      call bfineroot(dbh(j), i, canopy_trim, leaf_to_fineroot, elongation_factor_roots,  &
        fineroot_biomass(j, i))
      call bstore_allom(dbh(j), i, crown_damage, canopy_trim, bstore(j, i))
      call bdead_allom(bagw(j, i), bbgw(j, i), bsap(j, i), i, bdead(j, i))
      total_biom_parts(j, i) = blmax(j, i) + fineroot_biomass(j, i) + bagw(j, i) + bbgw(j, i)
      total_biom_tissues(j, i) = blmax(j, i) + fineroot_biomass(j, i) + bdead(j, i) + bsap(j, i)
    end do
  end do

  ! write out data to netcdf file
  call WriteAllometryData(out_file, numdbh, numpft, dbh, height, bagw, blmax, crown_area, &
    sapwood_area, bsap, bbgw, fineroot_biomass, bstore, bdead, total_biom_parts,         &
    total_biom_tissues)

end program FatesTestAllometry

! ----------------------------------------------------------------------------------------

subroutine WriteAllometryData(out_file, numdbh, numpft, dbh, height, bagw, blmax,        &
  crown_area, sapwood_area, bsap, bbgw, fineroot_biomass, bstore, bdead, total_biom_parts, &
  total_biom_tissues)
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
  integer,          intent(in) :: numdbh                   ! size of dbh array
  integer,          intent(in) :: numpft                   ! number of pfts
  real(r8),         intent(in) :: dbh(:)                   ! diameter at breast height [cm]
  real(r8),         intent(in) :: height(:,:)              ! height [m]
  real(r8),         intent(in) :: bagw(:,:)                ! aboveground biomass [kgC]
  real(r8),         intent(in) :: blmax(:, :)              ! leaf biomass [kgC]
  real(r8),         intent(in) :: crown_area(:, :)         ! crown area [m2]
  real(r8),         intent(in) :: sapwood_area(:, :)       ! sapwood cross-sectional area [m2]
  real(r8),         intent(in) :: bsap(:, :)               ! sapwood biomass [kgC]
  real(r8),         intent(in) :: bbgw(:, :)               ! belowground biomass [kgC]
  real(r8),         intent(in) :: fineroot_biomass(:, :)   ! fineroot biomass [kgC]
  real(r8),         intent(in) :: bstore(:, :)             ! storage biomass [kgC]
  real(r8),         intent(in) :: bdead(:, :)              ! deadwood biomass [kgC]
  real(r8),         intent(in) :: total_biom_parts(:, :)   ! total biomass calculated from parts [kgC]
  real(r8),         intent(in) :: total_biom_tissues(:, :) ! total biomass calculated from tissues [kgC]

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=8)     :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: dbhID, pftID   ! variable IDs for dimensions
  integer              :: heightID, bagwID
  integer              :: blmaxID, c_areaID
  integer              :: sapwoodareaID, bsapID
  integer              :: bbgwID, finerootID
  integer              :: bstoreID, bdeadID
  integer              :: totbiomID1, totbiomID2

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
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'cm', 'diameter at breast height'], 2, dbhID)

  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)

  ! register height
  call RegisterVar(ncid, 'height', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pft dbh', 'm', 'plant height'],      &
    3, heightID)

  ! register aboveground biomass
  call RegisterVar(ncid, 'bagw', dimIDs(1:2), type_double,                     &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                 &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant aboveground woody biomass'], &
    3, bagwID)

  ! register leaf biomass
  call RegisterVar(ncid, 'blmax', dimIDs(1:2), type_double,               &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant maximum leaf biomass'], &
    3, blmaxID)

  ! register crown area
  call RegisterVar(ncid, 'crown_area', dimIDs(1:2), type_double,          &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'pft dbh', 'm2', 'plant crown area per cohort'], &
    3, c_areaID)

  ! register sapwood area
  call RegisterVar(ncid, 'sapwood_area', dimIDs(1:2), type_double,                                 &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                     &
    [character(len=150) :: 'pft dbh', 'm2', 'plant cross section area sapwood at reference height'], &
    3, sapwoodareaID)

  ! register sapwood biomass
  call RegisterVar(ncid, 'bsap', dimIDs(1:2), type_double,           &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant sapwood biomass'], &
    3, bsapID)

  ! register belowground woody biomass
  call RegisterVar(ncid, 'bbgw', dimIDs(1:2), type_double,           &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant belowground woody biomass'], &
    3, bbgwID)

  ! register fineroot biomass
  call RegisterVar(ncid, 'fineroot_biomass', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant fineroot biomass'], &
    3, finerootID)

  ! register storage biomass
  call RegisterVar(ncid, 'bstore', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant storage biomass'], &
    3, bstoreID)

  ! register structural biomass
  call RegisterVar(ncid, 'bdead', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant deadwood (structural/heartwood) biomass'], &
    3, bdeadID)

  ! register total biomass (parts)
  call RegisterVar(ncid, 'total_biomass_parts', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant total biomass calculated from parts'], &
    3, totbiomID1)

  ! register total biomass (tissues)
  call RegisterVar(ncid, 'total_biomass_tissues', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],       &
    [character(len=150) :: 'pft dbh', 'kgC', 'plant total biomass calculated from tissues'], &
    3, totbiomID2)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, dbhID, dbh(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, heightID, height(:,:))
  call WriteVar(ncid, bagwID, bagw(:,:))
  call WriteVar(ncid, blmaxID, blmax(:,:))
  call WriteVar(ncid, c_areaID, crown_area(:,:))
  call WriteVar(ncid, sapwoodareaID, sapwood_area(:,:))
  call WriteVar(ncid, bsapID, bsap(:,:))
  call WriteVar(ncid, bbgwID, bbgw(:,:))
  call WriteVar(ncid, finerootID, fineroot_biomass(:,:))
  call WriteVar(ncid, bstoreID, bstore(:,:))
  call WriteVar(ncid, bdeadID, bdead(:,:))
  call WriteVar(ncid, totbiomID1, total_biom_parts(:,:))
  call WriteVar(ncid, totbiomID2, total_biom_tissues(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteAllometryData