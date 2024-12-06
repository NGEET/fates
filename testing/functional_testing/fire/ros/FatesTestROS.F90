program FatesTestROS
  
  use FatesConstantsMod,           only : r8 => fates_r8 
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader          ! param reader instance
  character(len=:), allocatable      :: param_file            ! input parameter file
  real(r8),         allocatable      :: SAV(:)                ! fuel surface area to volume ratio [/cm]
  real(r8),         allocatable      :: beta(:)               ! packing ratio [unitless]
  real(r8),         allocatable      :: propagating_flux(:,:) ! propagating flux [unitless]
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'ros_out.nc' ! output file 
  
  interface

    subroutine TestPropFlux(SAV, propagating_flux, beta)

      use FatesConstantsMod, only : r8 => fates_r8 
      use SFEquationsMod,    only : PropagatingFlux
      implicit none
      real(r8), allocatable, intent(out) :: SAV(:)                
      real(r8), allocatable, intent(out) :: propagating_flux(:,:)
      real(r8), allocatable, intent(out) :: beta(:) 

    end subroutine TestPropFlux
    
    subroutine WriteROSData(out_file, beta, SAV, propagating_flux)

      use FatesConstantsMod, only : r8 => fates_r8
      use FatesUnitTestIOMod,  only : OpenNCFile, CloseNCFile, RegisterNCDims
      use FatesUnitTestIOMod,  only : RegisterVar, EndNCDef, WriteVar
      use FatesUnitTestIOMod,  only : type_double
      implicit none
      character(len=*), intent(in) :: out_file
      real(r8),         intent(in) :: beta(:)
      real(r8),         intent(in) :: SAV(:)
      real(r8),         intent(in) :: propagating_flux(:,:)

    end subroutine WriteROSData

  end interface
  
  ! read in parameter file name and DATM file from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! calculate propagating flux
  call TestPropFlux(SAV, propagating_flux, beta)
  
  ! write output data
  call WriteROSData(out_file, beta, SAV, propagating_flux)
  
  ! deallocate arrays
  if (allocated(propagating_flux)) deallocate(propagating_flux)
  if (allocated(SAV)) deallocate(SAV)
  if (allocated(beta)) deallocate(beta)

end program FatesTestROS

!=========================================================================================

subroutine TestPropFlux(SAV, propagating_flux, beta)
  !
  ! DESCRIPTION:
  ! Calculates propagating flux ratio of a range of SAV and packing ratio values
  !
  use FatesConstantsMod, only : r8 => fates_r8 
  use SFEquationsMod,    only : PropagatingFlux
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: SAV(:)                ! fuel surface area to volume ratio [/cm]
  real(r8), allocatable, intent(out) :: propagating_flux(:,:) ! propagating flux [unitless]
  real(r8), allocatable, intent(out) :: beta(:)               ! packing ratio [unitless]

  ! CONSTANTS:
  real(r8), parameter               :: SAV_min = 0.0_r8   ! minimum SAV to calculate [/cm]
  real(r8), parameter               :: SAV_max = 115.0_r8 ! maximum SAV to calculate [/cm]
  real(r8), parameter               :: SAV_inc = 1.0_r8   ! SAV increment to scale [/cm]
  real(r8), parameter, dimension(4) :: packing_ratio = (/0.02_r8, 0.01_r8, 0.005_r8, 0.001_r8/) ! packing ratios to use [unitless]
  
  ! LOCALS:
  integer :: num_SAV ! size of SAV array
  integer :: i, j    ! looping indices
  
  ! allocate arrays
  num_SAV = int((SAV_max - SAV_min)/SAV_inc + 1)
  allocate(propagating_flux(num_SAV, size(packing_ratio)))
  allocate(SAV(num_SAV))
  allocate(beta(size(packing_ratio)))
  
  do i = 1, num_SAV
  
    SAV(i) = SAV_min + SAV_inc*(i-1)
    
    do j = 1, size(packing_ratio)
      beta(j) = packing_ratio(j)
      propagating_flux(i,j) = PropagatingFlux(packing_ratio(j), SAV(i))
    end do
  end do

end subroutine TestPropFlux

!=========================================================================================

subroutine WriteROSData(out_file, beta, SAV, propagating_flux)
  !
  ! DESCRIPTION:
  ! writes out data from the test
  !
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesUnitTestIOMod,  only : OpenNCFile, CloseNCFile, RegisterNCDims
  use FatesUnitTestIOMod,  only : RegisterVar, EndNCDef, WriteVar
  use FatesUnitTestIOMod,  only : type_double
  
  implicit none
  
  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file
  real(r8),         intent(in) :: beta(:)
  real(r8),         intent(in) :: SAV(:)
  real(r8),         intent(in) :: propagating_flux(:,:)
  
  ! LOCALS:
  integer           :: ncid         ! netcdf id
  character(len=20) :: dim_names(2) ! dimension names
  integer           :: dimIDs(2)    ! dimension IDs
  integer           :: SAVID
  integer           :: betaID
  integer           :: propfluxID
  
  
  ! dimension names
  dim_names = [character(len=20) :: 'SAV_value', 'beta_value']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/size(SAV), size(beta)/), 2, dimIDs)

  ! first register dimension variables
  
  ! register SAV 
  call RegisterVar(ncid, 'SAV', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'units', 'long_name'],           &
    [character(len=150) :: '/cm', 'fuel surface area to volume ratio'], 2, SAVID)
    
  ! register packing ratio
  call RegisterVar(ncid, 'packing_ratio', dimIDs(2:2), type_double,  &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: '', 'packing ratio'], 2, betaID)

  ! then register actual variables

  ! register propagating flux
  call RegisterVar(ncid, 'prop_flux', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
    [character(len=150) :: 'SAV_value beta_value', '', 'propagating flux'],  &
    3, propfluxID)
    
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, SAVID, SAV(:))
  call WriteVar(ncid, betaID, beta(:))
  call WriteVar(ncid, propfluxID, propagating_flux(:,:))
  
  ! close file
  call CloseNCFile(ncid)

end subroutine WriteROSData