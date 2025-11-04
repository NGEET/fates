program FatesTestROS
  
  use FatesConstantsMod,           only : r8 => fates_r8 
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader           ! param reader instance
  character(len=:), allocatable      :: param_file             ! input parameter file
  real(r8),         allocatable      :: reaction_intensity(:)
  real(r8),         allocatable      :: fuel_depth(:)
  real(r8),         allocatable      :: SAV(:)                 ! fuel surface area to volume ratio (for prop flux) [/cm]
  real(r8),         allocatable      :: beta(:)                ! packing ratio [unitless]
  real(r8),         allocatable      :: propagating_flux(:,:)  ! propagating flux [unitless]
  real(r8),         allocatable      :: SAV_values(:)          ! fuel surface area to volume ratio (for reaction vel) [/cm]
  real(r8),         allocatable      :: beta_ratio(:)          ! relative packing ratio [unitless]
  real(r8),         allocatable      :: reaction_velocity(:,:) ! reaction velocity [/min]
  real(r8),         allocatable      :: fuel_moisture(:)       ! fuel moisture [m3/m3]
  real(r8),         allocatable      :: q_ig(:)                ! heat of preignition [kJ/kg]
  real(r8),         allocatable      :: eps(:)                 ! effective heating number [unitless]
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'ros_out.nc' ! output file 
  
  interface
  
    subroutine TestReactionIntensity(fuel_depth, reaction_intensity)
      use FatesConstantsMod, only : r8 => fates_r8
      use SFEquationsMod,    only : ReactionIntensity
      implicit none
      real(r8), allocatable, intent(out) :: fuel_depth(:)        
      real(r8), allocatable, intent(out) :: reaction_intensity(:)
    
    end subroutine TestReactionIntensity

    subroutine TestPropFlux(SAV, propagating_flux, beta)

      use FatesConstantsMod, only : r8 => fates_r8 
      use SFEquationsMod,    only : PropagatingFlux
      implicit none
      real(r8), allocatable, intent(out) :: SAV(:)                
      real(r8), allocatable, intent(out) :: propagating_flux(:,:)
      real(r8), allocatable, intent(out) :: beta(:) 

    end subroutine TestPropFlux
    
    subroutine TestReactionVelocity(SAV, reaction_velocity, beta_ratio)

      use FatesConstantsMod, only : r8 => fates_r8 
      use SFEquationsMod,    only : OptimumReactionVelocity, MaximumReactionVelocity
      implicit none
      real(r8), allocatable, intent(out) :: SAV(:)                
      real(r8), allocatable, intent(out) :: reaction_velocity(:,:)
      real(r8), allocatable, intent(out) :: beta_ratio(:) 

    end subroutine TestReactionVelocity
    
    subroutine TestHeatofPreignition(fuel_moisture, q_ig)
      use FatesConstantsMod, only : r8 => fates_r8 
      use SFEquationsMod,    only : HeatofPreignition
      implicit none
      real(r8), allocatable, intent(out) :: fuel_moisture(:) 
      real(r8), allocatable, intent(out) :: q_ig(:)
    end subroutine TestHeatofPreignition
    
    subroutine TestEffectiveHeatingNumber(eps)
      use FatesConstantsMod, only : r8 => fates_r8 
      use SFEquationsMod,    only : EffectiveHeatingNumber
      implicit none
      
      real(r8), allocatable, intent(out) :: eps(:)
    end subroutine TestEffectiveHeatingNumber
    
    subroutine WriteROSData(out_file, reaction_intensity, fuel_depth, beta, SAV,         &
      propagating_flux, SAV_values, beta_ratio, reaction_velocity, fuel_moisture, q_ig,  &
      eps)

      use FatesConstantsMod, only : r8 => fates_r8
      use FatesUnitTestIOMod,  only : OpenNCFile, CloseNCFile, RegisterNCDims
      use FatesUnitTestIOMod,  only : RegisterVar, EndNCDef, WriteVar
      use FatesUnitTestIOMod,  only : type_double
      implicit none
      character(len=*), intent(in) :: out_file
      real(r8),         intent(in) :: reaction_intensity(:)
      real(r8),         intent(in) :: fuel_depth(:)
      real(r8),         intent(in) :: beta(:)
      real(r8),         intent(in) :: SAV(:)
      real(r8),         intent(in) :: propagating_flux(:,:)
      real(r8),         intent(in) :: SAV_values(:)
      real(r8),         intent(in) :: beta_ratio(:)
      real(r8),         intent(in) :: reaction_velocity(:,:)
      real(r8),         intent(in) :: fuel_moisture(:)
      real(r8),         intent(in) :: q_ig(:)
      real(r8),         intent(in) :: eps(:)

    end subroutine WriteROSData

  end interface
  
  ! read in parameter file name and DATM file from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! calculate reaction intensity
  call TestReactionIntensity(fuel_depth, reaction_intensity)
  
  ! calculate propagating flux
  call TestPropFlux(SAV, propagating_flux, beta)
  
  ! calculate reaction velocity
  call TestReactionVelocity(SAV_values, reaction_velocity, beta_ratio)
  
  ! calculate heat of preignition
  call TestHeatofPreignition(fuel_moisture, q_ig)
  
  ! calculate effective heating number
  call TestEffectiveHeatingNumber(eps)
  
  ! write output data
  call WriteROSData(out_file, reaction_intensity, fuel_depth, beta, SAV,                 &
    propagating_flux, SAV_values, beta_ratio, reaction_velocity, fuel_moisture, q_ig, eps)
  
  ! deallocate arrays
  if (allocated(fuel_depth)) deallocate(fuel_depth)
  if (allocated(reaction_intensity)) deallocate(reaction_intensity)
  if (allocated(propagating_flux)) deallocate(propagating_flux)
  if (allocated(SAV)) deallocate(SAV)
  if (allocated(beta)) deallocate(beta)
  if (allocated(reaction_velocity)) deallocate(reaction_velocity)
  if (allocated(SAV_values)) deallocate(SAV_values)
  if (allocated(beta_ratio)) deallocate(beta_ratio)
  if (allocated(fuel_moisture)) deallocate(fuel_moisture)
  if (allocated(q_ig)) deallocate(q_ig)
  if (allocated(eps)) deallocate(eps)

end program FatesTestROS

!=========================================================================================

subroutine TestReactionIntensity(fuel_depth, reaction_intensity)
  !
  ! DESCRIPTION:
  ! Calculates reaction intensity at different levels of fuel bed depth
  !
  use FatesConstantsMod, only : r8 => fates_r8
  use SFEquationsMod,    only : ReactionIntensity, OptimumPackingRatio
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: fuel_depth(:)         ! fuel depth [m]
  real(r8), allocatable, intent(out) :: reaction_intensity(:) ! reaction intensity [kJ/m2/min]
  
  ! CONSTANTS:
  real(r8), parameter :: fuel_loading = 1.4_r8 ! fuel loading [kg/m2]
  real(r8), parameter :: SAV = 63.0_r8         ! surface area to volume ratio [/cm]
  real(r8), parameter :: moisture = 0.04_r8    ! fuel moisture content [m3/m3]
  real(r8), parameter :: MEF = 0.25_r8         ! moisture of extinction [m3/m3]
  real(r8), parameter :: fd_min = 0.05         ! minimum fuel depth to test [m]
  real(r8), parameter :: fd_max = 2.5          ! maximum fuel depth to test [m]
  real(r8), parameter :: fd_inc = 0.05         ! fuel depth increment [m]
  real(r8), parameter :: part_dens = 513.0_r8  ! particle density [kg/m3]
  
  ! LOCALS:
  real(r8) :: beta_opt     ! optimum packing ratio
  real(r8) :: beta         ! packing ratio
  real(r8) :: beta_ratio   ! relative packing ratio
  real(r8) :: bulk_density ! fuel bulk densith [kg/m3]
  integer  :: num_fd       ! size of fuel depth array
  integer  :: i, j         ! looping indices
  
  ! allocate arrays
  num_fd = int((fd_max - fd_min)/fd_inc + 1)
  allocate(fuel_depth(num_fd))
  allocate(reaction_intensity(num_fd))
  
  beta_opt = OptimumPackingRatio(SAV)
  
  do i = 1, num_fd
  
    fuel_depth(i) = fd_min + fd_inc*(i-1)
    
    bulk_density = fuel_loading/fuel_depth(i)
    beta = bulk_density/part_dens
    beta_ratio = beta/beta_opt
  
    reaction_intensity(i) = ReactionIntensity(fuel_loading, SAV, beta_ratio, moisture,  &
      MEF)
  end do 
end subroutine TestReactionIntensity

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

subroutine TestReactionVelocity(SAV, reaction_velocity, beta_ratio)
  !
  ! DESCRIPTION:
  ! Calculates reaction velocity of a range of SAV and relative packing ratios
  !
  
  use FatesConstantsMod, only : r8 => fates_r8 
  use SFEquationsMod,    only : OptimumReactionVelocity, MaximumReactionVelocity
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: SAV(:)                 ! fuel surface area to volume ratio [/cm]
  real(r8), allocatable, intent(out) :: reaction_velocity(:,:) ! reaction velocity [/min]
  real(r8), allocatable, intent(out) :: beta_ratio(:)          ! relative packing ratio [unitless]

  ! CONSTANTS:
  real(r8), parameter               :: beta_r_min = 0.1_r8 ! minimum beta_ratio to calculate [unitless]
  real(r8), parameter               :: beta_r_max = 5.0_r8 ! maximum beta_ratio to calculate [unitless]
  real(r8), parameter               :: beta_r_inc = 0.1_r8 ! beta_ratio increment to scale [unitless]
  real(r8), parameter, dimension(5) :: SAV_vals = (/100.0_r8, 500.0_r8, 1000.0_r8, 2000.0_r8, 3000.0_r8/) ! SAV to use [/ft]
  
  ! LOCALS:
  integer  :: num_beta_r       ! size of beta_ratio array
  integer  :: i, j             ! looping indices
  real(r8) :: max_reaction_vel ! maximum reaction velocity [/min]
  
  ! allocate arrays
  num_beta_r = int((beta_r_max - beta_r_min)/beta_r_inc + 1)
  allocate(reaction_velocity(num_beta_r, size(SAV_vals)))
  allocate(beta_ratio(num_beta_r))
  allocate(SAV(size(SAV_vals)))
  
  do i = 1, num_beta_r
  
    beta_ratio(i) = beta_r_min + beta_r_inc*(i-1)
    
    do j = 1, size(SAV_vals)
      SAV(j) = SAV_vals(j)/30.48_r8 ! convert from /ft to /cm
      max_reaction_vel =  MaximumReactionVelocity(SAV(j))
      reaction_velocity(i,j) = OptimumReactionVelocity(max_reaction_vel, SAV(j),         &
        beta_ratio(i))
    end do
  end do

end subroutine TestReactionVelocity

!=========================================================================================

subroutine TestHeatofPreignition(fuel_moisture, q_ig)
  !
  ! DESCRIPTION:
  ! Calculates heat of preignition for a range of fuel moisture values
  !
  use FatesConstantsMod, only : r8 => fates_r8 
  use SFEquationsMod,    only : HeatofPreignition
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: fuel_moisture(:) ! fuel moisture [m3/m3]
  real(r8), allocatable, intent(out) :: q_ig(:)          ! heat of preignition [kJ/kg]

  ! CONSTANTS:
  real(r8), parameter :: moist_min = 0.0_r8   ! minimum fuel moisture to calculate [m3/m3]
  real(r8), parameter :: moist_max = 200.0_r8 ! maximum fuel moisture to calculate [m3/m3]
  real(r8), parameter :: moist_inc = 0.5_r8   ! fuel moisture increment to scale [m3/m3]
  
  ! LOCALS:
  integer :: num_moist ! size of fuel moisture array
  integer :: i         ! looping index
  
  ! allocate arrays
  num_moist = int((moist_max - moist_min)/moist_inc + 1)
  allocate(fuel_moisture(num_moist))
  allocate(q_ig(num_moist))
  
  do i = 1, num_moist
    fuel_moisture(i) = moist_min + moist_inc*(i-1)
    q_ig(i) = HeatofPreignition(fuel_moisture(i))
  end do

end subroutine TestHeatofPreignition

!=========================================================================================

subroutine TestEffectiveHeatingNumber(eps)
  !
  ! DESCRIPTION:
  ! Calculates effective heating number for a range of SAV values
  !
  use FatesConstantsMod, only : r8 => fates_r8 
  use SFEquationsMod,    only : EffectiveHeatingNumber
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), allocatable, intent(out) :: eps(:) ! effective heating number [unitless]

  ! CONSTANTS:
  real(r8), parameter :: SAV_min = 0.0_r8   ! minimum SAV to calculate [/cm]
  real(r8), parameter :: SAV_max = 115.0_r8 ! maximum SAV to calculate [/cm]
  real(r8), parameter :: SAV_inc = 1.0_r8   ! SAV increment to scale [/cm]
  
  ! LOCALS:
  real(r8), allocatable :: SAV(:)  ! fuel surface area to volume ratio [/cm]
  integer               :: num_SAV ! size of SAV array
  integer               :: i       ! looping index
  
  ! allocate arrays
  num_SAV = int((SAV_max - SAV_min)/SAV_inc + 1)
  allocate(SAV(num_SAV))
  allocate(eps(num_SAV))
  
  do i = 1, num_SAV
    SAV(i) = SAV_min + SAV_inc*(i-1)
    eps(i) = EffectiveHeatingNumber(SAV(i))
  end do

end subroutine TestEffectiveHeatingNumber

!=========================================================================================

subroutine WriteROSData(out_file, reaction_intensity, fuel_depth, beta, SAV,             &
  propagating_flux, SAV_values, beta_ratio, reaction_velocity, fuel_moisture, q_ig, eps)
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
  real(r8),         intent(in) :: reaction_intensity(:)
  real(r8),         intent(in) :: fuel_depth(:)
  real(r8),         intent(in) :: beta(:)
  real(r8),         intent(in) :: SAV(:)
  real(r8),         intent(in) :: propagating_flux(:,:)
  real(r8),         intent(in) :: SAV_values(:)
  real(r8),         intent(in) :: beta_ratio(:)
  real(r8),         intent(in) :: reaction_velocity(:,:)
  real(r8),         intent(in) :: fuel_moisture(:)
  real(r8),         intent(in) :: q_ig(:)
  real(r8),         intent(in) :: eps(:)
  
  ! LOCALS:
  integer           :: ncid         ! netcdf id
  character(len=20) :: dim_names(6) ! dimension names
  integer           :: dimIDs(6)    ! dimension IDs
  integer           :: fuel_depthID
  integer           :: reactionintenseID
  integer           :: SAVID
  integer           :: betaID
  integer           :: propfluxID
  integer           :: SAVindID
  integer           :: beta_r_ID
  integer           :: reactionvelID
  integer           :: moistID
  integer           :: qigID
  integer           :: epsID
  
  ! dimension names
  dim_names = [character(len=20) :: 'SAV', 'packing_ratio', 'SAV_ind',                   &
    'beta_ratio', 'fuel_moisture', 'fuel_depth']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/size(SAV), size(beta),        &
    size(SAV_values), size(beta_ratio), size(fuel_moisture), size(fuel_depth)/), 6, dimIDs)

  ! first register dimension variables
  
  ! register fuel depth
  call RegisterVar(ncid, 'fuel_depth', dimIDs(6:6), type_double,   &
    [character(len=20)  :: 'units', 'long_name'],           &
    [character(len=150) :: 'm', 'fuel bed depth'], 2, fuel_depthID)
  
  ! register SAV 
  call RegisterVar(ncid, 'SAV', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'units', 'long_name'],           &
    [character(len=150) :: '/cm', 'fuel surface area to volume ratio'], 2, SAVID)
    
  ! register packing ratio
  call RegisterVar(ncid, 'packing_ratio', dimIDs(2:2), type_double,  &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: '', 'packing ratio'], 2, betaID)
    
  ! register SAV values 
  call RegisterVar(ncid, 'SAV_ind', dimIDs(3:3), type_double,   &
    [character(len=20)  :: 'units', 'long_name'],                &
    [character(len=150) :: '/cm', 'fuel surface area to volume ratio'], 2, SAVindID)
    
  ! register relative packing ratio
  call RegisterVar(ncid, 'beta_ratio', dimIDs(4:4), type_double,     &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: '', 'relative packing ratio'], 2, beta_r_ID)
    
  ! register fuel moisture
  call RegisterVar(ncid, 'fuel_moisture', dimIDs(5:5), type_double,     &
    [character(len=20)  :: 'units', 'long_name'],                    &
    [character(len=150) :: 'm3/m3', 'fuel moisture'], 2, moistID)

  ! then register actual variables
  
  ! register reaction intensity
  call RegisterVar(ncid, 'reaction_intensity', dimIDs(6:6), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
    [character(len=150) :: 'fuel_depth', 'kJ/m2/min', 'reaction intensity'],  &
    3, reactionintenseID)

  ! register propagating flux
  call RegisterVar(ncid, 'prop_flux', dimIDs(1:2), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],        &
    [character(len=150) :: 'SAV packing_ratio', '', 'propagating flux'],  &
    3, propfluxID)
    
  ! register reaction velocity 
  call RegisterVar(ncid, 'reaction_velocity', (/dimIDs(4), dimIDs(3)/), type_double,  &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],         &
    [character(len=150) :: 'SAV_ind beta_ratio', '/min', 'reaction velocity'], &
    3, reactionvelID)
    
  ! register heat of preignition
  call RegisterVar(ncid, 'q_ig', dimIDs(5:5), type_double,  &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],         &
    [character(len=150) :: 'fuel_moisture', 'kJ/kg', 'heat of preignition'], &
    3, qigID)
    
  ! register effective heating number
  call RegisterVar(ncid, 'eps', dimIDs(1:1), type_double,  &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],         &
    [character(len=150) :: 'SAV', '', 'effective heating number'], &
    3, epsID)
    
    
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, SAVID, SAV(:))
  call WriteVar(ncid, betaID, beta(:))
  call WriteVar(ncid, propfluxID, propagating_flux(:,:))
  call WriteVar(ncid, SAVindID, SAV_values(:))
  call WriteVar(ncid, beta_r_ID, beta_ratio(:))
  call WriteVar(ncid, reactionvelID, reaction_velocity(:,:))
  call WriteVar(ncid, moistID, fuel_moisture(:))
  call WriteVar(ncid, qigID, q_ig(:))
  call WriteVar(ncid, epsID, eps(:))
  call WriteVar(ncid, fuel_depthID, fuel_depth(:))
  call WriteVar(ncid, reactionintenseID, reaction_intensity(:))
  
  ! close file
  call CloseNCFile(ncid)

end subroutine WriteROSData