program FatesTestFireMortality
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use PRTParametersMod,            only : prt_params
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader               ! param reader instance
  character(len=:), allocatable      :: param_file                 ! input parameter file
  real(r8),         allocatable      :: tree_diameter(:)           ! tree diameter at breast height [cm]
  real(r8),         allocatable      :: tau_c(:,:)                 ! critical residence time for cambial death [min]
  real(r8),         allocatable      :: mortality(:,:,:)           ! probability of fire mortality
  real(r8),         allocatable      :: crown_kill(:)              ! fraction of crown volume burned
  real(r8),         allocatable      :: tau_r_out(:)               ! relative fire residence time
  real(r8),         allocatable      :: fire_mortality_bySH(:, :)  ! total fire mortality probability (by SH)
  real(r8),         allocatable      :: fire_mortality_bytau(:, :) ! total fire mortality probability (by tau)
  real(r8),         allocatable      :: SH(:)                      ! scorch height [m]
  real(r8),         allocatable      :: tau_l(:)                   ! residence time of fire [min]
  integer                            :: num_pfts                   ! number of pfts (from parameter files)
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'fire_mortality_out.nc' ! output file 
  
  interface

    subroutine TestTauC(num_pfts, tree_diameter, tau_c)

      use FatesConstantsMod, only : r8 => fates_r8 
      use FatesConstantsMod, only : itrue
      use SFEquationsMod,    only : CriticalResidenceTime, BarkThickness
      use EDPftvarcon,       only : EDPftvarcon_inst
      use PRTParametersMod,  only : prt_params
      implicit none
      integer,               intent(in)  :: num_pfts
      real(r8), allocatable, intent(out) :: tree_diameter(:)                
      real(r8), allocatable, intent(out) :: tau_c(:,:) 

    end subroutine TestTauC
    
    subroutine TestMortalityProb(num_pfts, mortality, crown_kill, tau_r_out)
      use FatesConstantsMod, only : r8 => fates_r8
      use FatesConstantsMod, only : itrue
      use SFEquationsMod,    only : CriticalResidenceTime, BarkThickness
      use SFEquationsMod,    only : TotalFireMortality, CrownFireMortality
      use SFEquationsMod,    only : cambial_mort
      use EDPftvarcon,       only : EDPftvarcon_inst
      use PRTParametersMod,  only : prt_params
      implicit none
      integer,               intent(in)  :: num_pfts         
      real(r8), allocatable, intent(out) :: mortality(:,:,:) 
      real(r8), allocatable, intent(out) :: crown_kill(:)
      real(r8), allocatable, intent(out) :: tau_r_out(:)
    end subroutine TestMortalityProb
    
    subroutine TestFireMortalitySensitivity(fire_mortality_bySH, fire_mortality_bytau, &
        SH_out, tau_l_out)
      use FatesConstantsMod, only : r8 => fates_r8
      use SFEquationsMod,    only : CrownFractionBurnt, CambialMortality
      use SFEquationsMod,    only : CrownFireMortality, TotalFireMortality
      implicit none
      real(r8), allocatable, intent(out) :: fire_mortality_bySH(:,:)
      real(r8), allocatable, intent(out) :: fire_mortality_bytau(:,:)
      real(r8), allocatable, intent(out) :: SH_out(:)
      real(r8), allocatable, intent(out) :: tau_l_out(:)
    end subroutine TestFireMortalitySensitivity
    
    subroutine WriteFireMortData(out_file, num_pfts, tree_diameter, tau_c, mortality, &
        crown_kill, tau_r_out, fire_mortality_bySH, SH, fire_mortality_bytau, tau_l)

      use FatesConstantsMod,   only : r8 => fates_r8
      use FatesUnitTestIOMod,  only : OpenNCFile, CloseNCFile, RegisterNCDims
      use FatesUnitTestIOMod,  only : RegisterVar, EndNCDef, WriteVar
      use FatesUnitTestIOMod,  only : type_double
      implicit none
      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: num_pfts
      real(r8),         intent(in) :: tree_diameter(:)
      real(r8),         intent(in) :: tau_c(:,:)
      real(r8),         intent(in) :: mortality(:,:,:) 
      real(r8),         intent(in) :: crown_kill(:)
      real(r8),         intent(in) :: tau_r_out(:)
      real(r8),         intent(in) :: fire_mortality_bySH(:,:)
      real(r8),         intent(in) :: SH(:)
      real(r8),         intent(in) :: fire_mortality_bytau(:,:)
      real(r8),         intent(in) :: tau_l(:)

    end subroutine WriteFireMortData

  end interface
  
  ! read in parameter file name and DATM file from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  num_pfts = size(prt_params%wood_density, dim=1)
  
  ! calculate propagating flux
  call TestTauC(num_pfts, tree_diameter, tau_c)
  
  ! calculate total fire mortality
  call TestMortalityProb(num_pfts, mortality, crown_kill, tau_r_out)
  
  ! test sensitivity to SH and cambial kill
  call TestFireMortalitySensitivity(fire_mortality_bySH, fire_mortality_bytau, SH, tau_l)
  
  ! write output data
  call WriteFireMortData(out_file, num_pfts, tree_diameter, tau_c, mortality, &
    crown_kill, tau_r_out, fire_mortality_bySH, SH, fire_mortality_bytau, tau_l)
    
  ! deallocate arrays
  if (allocated(tree_diameter)) deallocate(tree_diameter)
  if (allocated(tau_c)) deallocate(tau_c)
  if (allocated(mortality)) deallocate(mortality)
  if (allocated(crown_kill)) deallocate(crown_kill)
  if (allocated(tau_r_out)) deallocate(tau_r_out)
  if (allocated(fire_mortality_bySH)) deallocate(fire_mortality_bySH)
  if (allocated(fire_mortality_bytau)) deallocate(fire_mortality_bytau)
  if (allocated(SH)) deallocate(SH)
  if (allocated(tau_l)) deallocate(tau_l)

end program FatesTestFireMortality

!=========================================================================================

subroutine TestFireMortalitySensitivity(fire_mortality_bySH, fire_mortality_bytau, SH_out, &
  tau_l_out)
  !
  ! DESCRIPTION:
  ! Calculates fire mortality for some input tree characteristics and fire behaviors - based
  ! on scorch height
  !
  use FatesConstantsMod, only : r8 => fates_r8
  use SFEquationsMod,    only : CrownFractionBurnt, CambialMortality
  use SFEquationsMod,    only : CrownFireMortality, TotalFireMortality
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: fire_mortality_bySH(:,:)  ! total fire mortality (by SH)
  real(r8), allocatable, intent(out) :: fire_mortality_bytau(:,:) ! total fire mortality (by tau)
  real(r8), allocatable, intent(out) :: SH_out(:)                 ! scorch height of tree [m]
  real(r8), allocatable, intent(out) :: tau_l_out(:)              ! residence times of fire [min]
  
  ! LOCALS:
  integer  :: i, j                  ! looping indices
  real(r8) :: fraction_crown_burned ! fraction of the crown burned
  real(r8) :: cambial_mort          ! cambial mortality
  real(r8) :: bark_scaler           ! cm bark per cm dbh
  real(r8) :: crownfire_mort        ! crown fire mortality
  
  ! CONSTANTS:
  real(r8), parameter, dimension(3) :: SH = (/5.0_r8, 10.0_r8, 20.0_r8/)                                       ! scorch heights of fire [m]
  real(r8), parameter, dimension(3) :: tau_l = (/2.4_r8, 3.1_r8, 6.0_r8/)                                      ! residence times of fire [min]
  real(r8), parameter, dimension(6) :: dbh = (/20.0_r8, 40.0_r8, 20.0_r8, 40.0_r8, 20.0_r8, 40.0_r8/)          ! diameters of trees [cm]
  real(r8), parameter, dimension(6) :: crown_length = (/12.0_r8, 18.5_r8, 13.5_r8, 23.1_r8, 11.9_r8, 18.2_r8/) ! crown depth of trees [m]
  real(r8), parameter, dimension(6) :: bark_thickness = (/1.1_r8, 2.2_r8, 0.9_r8, 1.7_r8, 0.3_r8, 0.6_r8/)     ! bark thickness of trees [cm]
  real(r8), parameter, dimension(6) :: height = (/15.5_r8, 24.4_r8, 17.4_r8, 30.7_r8, 15.1_r8, 24.1_r8/)       ! height of trees [m]
  real(r8), parameter               :: crown_kill_parameter = 0.775_r8 ! parameter for crown kill
  
  allocate(fire_mortality_bySH(size(SH), size(dbh)))
  allocate(fire_mortality_bytau(size(tau_l), size(dbh)))
  allocate(SH_out(size(SH)))
  allocate(tau_l_out(size(tau_l)))
  
  do i = 1, size(SH)
    SH_out(i) = SH(i)
    do j = 1, size(dbh)
      bark_scaler = bark_thickness(j)/dbh(j)
      fraction_crown_burned = CrownFractionBurnt(SH(i), height(j), crown_length(j))
      cambial_mort = CambialMortality(bark_scaler, dbh(j), 6.0_r8)
      crownfire_mort = CrownFireMortality(crown_kill_parameter, fraction_crown_burned)
      fire_mortality_bySH(i,j) = TotalFireMortality(crownfire_mort, cambial_mort)
    end do 
  end do
  
  do i = 1, size(tau_l)
    tau_l_out(i) = tau_l(i)
    do j = 1, size(dbh)
      bark_scaler = bark_thickness(j)/dbh(j)
      fraction_crown_burned = CrownFractionBurnt(10.0_r8, height(j), crown_length(j))
      cambial_mort = CambialMortality(bark_scaler, dbh(j), tau_l(i))
      crownfire_mort = CrownFireMortality(crown_kill_parameter, fraction_crown_burned)
      fire_mortality_bytau(i,j) = TotalFireMortality(crownfire_mort, cambial_mort)
    end do 
  end do
  
end subroutine TestFireMortalitySensitivity

!=========================================================================================

subroutine TestTauC(num_pfts, tree_diameter, tau_c)
  !
  ! DESCRIPTION:
  ! Calculates critical time for cambial kill over a range of diameter values and pfts
  !
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : itrue
  use SFEquationsMod,    only : CriticalResidenceTime, BarkThickness
  use EDPftvarcon,       only : EDPftvarcon_inst
  use PRTParametersMod,  only : prt_params
  
  implicit none
  
  ! ARGUMENTS:
  integer,               intent(in)  :: num_pfts         ! number of pfts
  real(r8), allocatable, intent(out) :: tree_diameter(:) ! tree diameter at breast height [cm]
  real(r8), allocatable, intent(out) :: tau_c(:,:)       ! critical fire residence time for cambial kill [min]

  ! CONSTANTS:
  real(r8), parameter :: dbh_min = 2.5_r8  ! minimum dbh to calculate [cm]
  real(r8), parameter :: dbh_max = 60.0_r8 ! maximum dbh to calculate [cm]
  real(r8), parameter :: dbh_inc = 1.0_r8  ! dbh increment to scale [cm]
  
  ! LOCALS:
  real(r8) :: bark_thickness ! bark thickness [cm]
  integer  :: num_dbh        ! size of dbh array
  integer  :: i, j           ! looping indices
  
  ! allocate arrays
  num_dbh = int((dbh_max - dbh_min)/dbh_inc + 1)
  allocate(tree_diameter(num_dbh))
  allocate(tau_c(num_dbh, num_pfts))
  
  do i = 1, num_dbh
  
    tree_diameter(i) = dbh_min + dbh_inc*(i-1)
    
    do j = 1, num_pfts
      if (prt_params%woody(j) == itrue) then 
        bark_thickness = BarkThickness(EDPftvarcon_inst%bark_scaler(j), tree_diameter(i))
        tau_c(i,j) = CriticalResidenceTime(bark_thickness)
      else 
        tau_c(i,j) = 0.0_r8
      end if
    end do
  end do

end subroutine TestTauC

!=========================================================================================

subroutine TestMortalityProb(num_pfts, mortality, crown_kill, tau_r_out)
  !
  ! DESCRIPTION:
  ! Calculates mortality probability for a range of values of fraction crown killed and 
  ! relative residence time
  !
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : itrue
  use SFEquationsMod,    only : CriticalResidenceTime, BarkThickness
  use SFEquationsMod,    only : TotalFireMortality, CrownFireMortality
  use SFEquationsMod,    only : cambial_mort
  use EDPftvarcon,       only : EDPftvarcon_inst
  use PRTParametersMod,  only : prt_params
  
  implicit none
  
  ! ARGUMENTS:
  integer,               intent(in)  :: num_pfts         ! number of pfts
  real(r8), allocatable, intent(out) :: mortality(:,:,:) ! probability of fire mortality
  real(r8), allocatable, intent(out) :: crown_kill(:)    ! fraction of crown volume burned
  real(r8), allocatable, intent(out) :: tau_r_out(:)     ! relative fire residence times

  ! CONSTANTS:
  real(r8), parameter               :: mort_min = 0.0_r8  ! minimum mortality rate to calculate
  real(r8), parameter               :: mort_max = 1.0_r8  ! maximum mortality rate to calculate
  real(r8), parameter               :: mort_inc = 0.05_r8 ! mortality rates to scale
  real(r8), parameter, dimension(5) :: tau_r = (/0.22_r8, 0.4_r8, 0.66_r8, 1.0_r8, 2.0_r8/) ! relative fire residence times
  
  ! LOCALS:
  real(r8) :: crownfire_mort      ! crown fire mortality
  real(r8) :: cambial_damage_mort ! cambial damage mortality
  integer  :: num_mort            ! size of dbh array
  integer  :: i, j, k             ! looping indices
  
  ! allocate arrays
  num_mort = int((mort_max - mort_min)/mort_inc + 1)
  allocate(mortality(num_mort, num_pfts, size(tau_r)))
  allocate(crown_kill(num_mort))
  allocate(tau_r_out(size(tau_r)))
  
  do i = 1, num_mort
  
    crown_kill(i) = mort_min + mort_inc*(i-1)
    
    do j = 1, num_pfts
      do k = 1, size(tau_r)
        tau_r_out(k) = tau_r(k)
        
        if (prt_params%woody(j) == itrue) then
          cambial_damage_mort = cambial_mort(tau_r(k))
          crownfire_mort = CrownFireMortality(EDPftvarcon_inst%crown_kill(j), crown_kill(i))
          mortality(i,j,k) = TotalFireMortality(crownfire_mort, cambial_damage_mort)
        else 
          mortality(i,j,k) = 0.0_r8
        end if
      end do
    end do
  end do

end subroutine TestMortalityProb

!=========================================================================================

subroutine WriteFireMortData(out_file, num_pfts, tree_diameter, tau_c, mortality, &
  crown_kill, tau_r_out, fire_mortality_bySH, SH, fire_mortality_bytau, tau_l)
  !
  ! DESCRIPTION:
  ! writes out data from the test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, CloseNCFile, RegisterNCDims
  use FatesUnitTestIOMod, only : RegisterVar, EndNCDef, WriteVar
  use FatesUnitTestIOMod, only : type_double, type_int
  
  implicit none
  
  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file
  integer,          intent(in) :: num_pfts
  real(r8),         intent(in) :: tree_diameter(:)
  real(r8),         intent(in) :: tau_c(:,:)
  real(r8),         intent(in) :: mortality(:,:,:)
  real(r8),         intent(in) :: crown_kill(:)
  real(r8),         intent(in) :: tau_r_out(:)
  real(r8),         intent(in) :: fire_mortality_bySH(:,:)
  real(r8),         intent(in) :: SH(:)
  real(r8),         intent(in) :: fire_mortality_bytau(:,:)
  real(r8),         intent(in) :: tau_l(:)
  
  ! LOCALS:
  integer, allocatable :: pft_indices(:)  ! array of pft indices to write out
  integer, allocatable :: tree_indices(:) ! array of tree indices to write out
  integer              :: ncid            ! netcdf id
  character(len=20)    :: dim_names(7)    ! dimension names
  integer              :: dimIDs(7)       ! dimension IDs
  integer              :: i               ! looping index
  integer              :: dbhID, pftID
  integer              :: taurID
  integer              :: taucID
  integer              :: mortID, crownkillID
  integer              :: SHID, firemortbySHID
  integer              :: treeID, firemortbytauID
  integer              :: taulID
  
  ! dimension names
  dim_names = [character(len=20) :: 'dbh', 'pft', 'crown_kill', 'tau_r', 'SH', 'treeID', 'tau_l']
  
  ! create pft indices
  allocate(pft_indices(num_pfts))
  do i = 1, num_pfts
    pft_indices(i) = i
  end do
  
  ! create tree indices
  allocate(tree_indices(size(fire_mortality_bySH, dim=2)))
  do i = 1, size(fire_mortality_bySH, dim=2)
    tree_indices(i) = i
  end do

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/size(tree_diameter), num_pfts, size(crown_kill), &
    size(tau_r_out), size(SH), size(tree_indices), size(tau_l)/), size(dim_names), dimIDs)

  ! first register dimension variables
  
  ! register dbh 
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: 'cm', 'diameter at breast height'], 2, dbhID)
    
  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
    
  ! register crown kill
  call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double,     &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: '', 'fraction crown volume burned'], 2, crownkillID)
  
  ! register tau_r
  call RegisterVar(ncid, dim_names(4), dimIDs(4:4), type_double,       &
    [character(len=20)  :: 'units', 'long_name'],                      &
    [character(len=150) :: '', 'relative fire residence time'], 2, taurID)
    
  ! register scorch height
  call RegisterVar(ncid, dim_names(5), dimIDs(5:5), type_double,       &
    [character(len=20)  :: 'units', 'long_name'],                      &
    [character(len=150) :: 'm', 'scorch height'], 2, SHID)
    
  ! register tree ids
  call RegisterVar(ncid, dim_names(6), dimIDs(6:6), type_int,        &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: '', 'tree indices'], 2, treeID)
    
  ! register tau_l
  call RegisterVar(ncid, dim_names(7), dimIDs(7:7), type_double,     &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: 'min', 'residence time of fire'], 2, taulID)
    
  ! then register actual variables

  ! register tau_c
  call RegisterVar(ncid, 'tau_c', dimIDs(1:2), type_double,                                 &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                            &
    [character(len=150) :: 'pft dbh', 'min', 'critical residence time for cambial death'],  &
    3, taucID)
    
  ! register total mortality
  call RegisterVar(ncid, 'total_mortality', (/dimIDs(3), dimIDs(2), dimIDs(4)/), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                              &
    [character(len=150) :: 'crown_kill pft tau_r', '', 'fire mortality'],                     &
    3, mortID)
    
  ! register total mortality
  call RegisterVar(ncid, 'fire_mortality_bySH', (/dimIDs(5), dimIDs(6)/), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                       &
    [character(len=150) :: 'SH treeID', '', 'fire mortality by SH'],                         &
    3, firemortbySHID)
    
  ! register total mortality
  call RegisterVar(ncid, 'fire_mortality_bytau', (/dimIDs(7), dimIDs(6)/), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                       &
    [character(len=150) :: 'tau_l treeID', '', 'fire mortality by tau'],                         &
    3, firemortbytauID)
      
      
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, dbhID, tree_diameter(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, crownkillID, crown_kill(:))
  call WriteVar(ncid, taurID, tau_r_out(:))
  call WriteVar(ncid, taucID, tau_c(:,:))
  call WriteVar(ncid, mortID, mortality(:,:,:))
  call WriteVar(ncid, SHID, SH(:))
  call WriteVar(ncid, taulID, tau_l(:))
  call WriteVar(ncid, treeID, tree_indices(:))
  call WriteVar(ncid, firemortbySHID, fire_mortality_bySH(:,:))
  call WriteVar(ncid, firemortbytauID, fire_mortality_bytau(:,:))
  
  ! close file
  call CloseNCFile(ncid)

end subroutine WriteFireMortData