program FatesUnitTestRadiation
  !
  ! DESCRIPTION:
  !		Test the FATES radiation schemes
  !

  use FatesUnitTestIOMod,       only : logf, MAX_PATH, open_file, read_patch_data
  use FatesUnitTestOrbitalMod,  only : SHR_KIND_R8
  use FatesUnitTestOrbitalMod,  only : get_orbital_vals, shr_orb_cosz
  use FatesUnitTestOrbitalMod,  only : shr_orb_decl
  use FatesConstantsMod,        only : r8 => fates_r8
  use EDPftvarcon,              only : EDPftvarcon_inst
  use FatesParametersInterface, only : fates_parameters_type
  use FatesPatchMod,            only : fates_patch_type
  use PRTGenericMod,            only : num_elements, element_list, element_pos
  use PRTGenericMod,            only : carbon12_element

  implicit none
    
  ! LOCALS:
  type(fates_patch_type), pointer :: fates_patch        ! patch object
  character(len=MAX_PATH)              :: patch_file         ! patch data file name
  character(len=MAX_PATH)      :: out_file           ! output file name
  integer                      :: year, jday         ! year and day of year to simulate
  integer                      :: num_pft            ! number of pfts
  real(r8)                     :: lat, lon           ! latitude/longitude to simulate [degrees]
  real(r8)                     :: fcansno            ! fraction of canopy covered by snow [0-1]
  real(r8), allocatable        :: canopy_area(:,:,:) ! canopy area profile
  real(r8), allocatable        :: elai(:,:,:)        ! exposed lai profile
  real(r8), allocatable        :: esai(:,:,:)        ! exposed sai profile
  integer,  allocatable        :: nrad(:,:)          ! number of exposed leaf layers

  ! PARAMETERS
  integer, parameter :: numSWB   = 2    ! number of shortwave bands to simulate
  integer, parameter :: nlevsoil = 5    ! this shouldn't really matter
  integer, parameter :: current_tod = 0 ! current time [seconds past 0Z]

  interface

    subroutine read_radiation_namelist(year, jday, lat, lon, num_pft, fcansno, &
        patch_file, out_file)

      use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file
      use FatesConstantsMod,  only : r8 => fates_r8

      implicit none

      character(len=MAX_PATH), intent(out) :: patch_file, out_file 
      integer,                 intent(out) :: year, jday
      integer,                 intent(out) :: num_pft
      real(r8),                intent(out) :: lat, lon
      real(r8),                intent(out) :: fcansno

    end subroutine read_radiation_namelist

    subroutine read_pft_namelist(num_pft, numSWB, clumping_index, xl, rhol)
      
      use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file, logf
      use FatesConstantsMod,  only : r8 => fates_r8

      implicit none
      
      integer,  intent(in)  :: num_pft           
      integer,  intent(in)  :: numSWB                  
      real(r8), intent(out) :: clumping_index(num_pft) 
      real(r8), intent(out) :: xl(num_pft)            
      real(r8), intent(out) :: rhol(num_pft,numSWB)   
 
    end subroutine read_pft_namelist
    subroutine init_pft_data(num_pft, numSWB)
      
      use EDPftvarcon,              only : EDPftvarcon_inst
      use FatesParametersInterface, only : fates_parameters_type
      use FatesConstantsMod,        only : r8 => fates_r8

      implicit none

      integer, intent(in) :: num_pft
      integer, intent(in) :: numSWB
    
    end subroutine init_pft_data

  end interface

  !:...........................................................................:
  
  ! open log file
  logf = open_file("log.txt")

  ! read in namelist to get some runtime parameters
  call read_radiation_namelist(year, jday, lat, lon, num_pft, fcansno,         &
    patch_file, out_file)

  ! initialize pft data
  call init_pft_data(num_pft, numSWB)
    
  ! read in patch data
  call read_patch_data(patch_file, num_pft, canopy_area, elai, esai, nrad)

  ! initialize a patch and set values
  allocate(fates_patch)
  call fates_patch%init(numSWB, nlevsoil)
  fates_patch%canopy_area_profile = canopy_area
  fates_patch%elai_profile = elai 
  fates_patch%esai_profile = esai
  fates_patch%nrad = nrad

  ! call get_orbital_vals(year, logf, eccen, mvelp, obliqr, lambm0, mvelpp)

  ! cosz(:) = 0.0
  ! k_dir(:,:) = 0.0
  
  ! ! for each half-hourly time step in the day
  ! do i = 0, 47
  !     calday = jday + i*0.02083333_SHR_KIND_R8
  !     call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
  !     cosz(i) = shr_orb_cosz(calday, lat, lon, declin)

  !     if (cosz(i) > 0.0_r8) then 
  !         ! call norman radiation scheme
  !         call PatchNormanRadiation(rhol, rhos, taul, taus, xl, clumping_index, &
  !             canopy_area_profile, elai_profile, esai_profile, fcansno,         &
  !             cosz(i), nrad, 2, k_dir(:,i))
  !     end if
  ! end do 

  ! call write_radiation_data(out_file, k_dir, cosz)

  close(logf)

end program FatesUnitTestRadiation

!:.............................................................................:

subroutine read_radiation_namelist(year, jday, lat, lon, num_pft, fcansno,     &
    patch_file, out_file)
  !
  ! DESCRIPTION:
  !		read in the namelist associated with the radiation unit tests and 
  !   initialize values
  !

  use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file, logf
  use FatesConstantsMod,  only : r8 => fates_r8
  
  implicit none

  ! ARGUMENTS:
  character(len=MAX_PATH), intent(out) :: patch_file ! patch data file name
  character(len=MAX_PATH), intent(out) :: out_file   ! output file name
  integer,                 intent(out) :: year, jday ! year and day of year
  integer,                 intent(out) :: num_pft    ! number of pfts
  real(r8),                intent(out) :: lat, lon   ! latitude and longitude [degrees]
  real(r8),                intent(out) :: fcansno    ! fraction of canopy covered by snow [0-1]
  
  ! LOCALS:
  character(len=MAX_PATH) :: rad_nl = 'radiation_nl' ! radiation namelist name
  character(len=MAX_CHAR) :: msg                     ! I/O Error message
  integer                 :: rad_nl_file             ! unit number for namelist
  integer                 :: ios                     ! I/O status
  
  ! Namelist of radiation parameters
  namelist /radiation/ year, jday, lat, lon, num_pft, fcansno, patch_file,     &
    out_file

  ! Now read parameters namelist
  rad_nl_file = open_file(trim(rad_nl), 'r')
  read(rad_nl_file, radiation, iostat=ios, iomsg=msg)
  
  if (ios /= 0) then
    ! Problem reading file - tell user.
    write(logf, '(A, I6, A, A)') "Error reading radiation namelist file",      &
      ios, " IOMSG: ", msg
    stop "Stopped"
  end if

  close(rad_nl_file)

end subroutine read_radiation_namelist

!:.............................................................................:

subroutine read_pft_namelist(num_pft, numSWB, clumping_index, xl, rhol)
!
! DESCRIPTION:
!		read in the namelist associated with the pft-specific parameters
!

use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file, logf
use FatesConstantsMod,  only : r8 => fates_r8

implicit none

! ARGUMENTS:
integer,  intent(in)  :: num_pft                 ! number of pfts
integer,  intent(in)  :: numSWB                  ! number of shortwave bands to simulate
real(r8), intent(out) :: clumping_index(num_pft) ! clumping index
real(r8), intent(out) :: xl(num_pft)             ! leaf-stem orientation index
real(r8), intent(out) :: rhol(num_pft,numSWB)    ! leaf reflectance [0-1]

! LOCALS:
character(len=MAX_PATH) :: pft_nl = 'pft_nlfile'   ! pft namelist name
character(len=MAX_CHAR) :: msg                     ! I/O Error message
integer                 :: pft_nl_file             ! unit number for namelist
integer                 :: ios                     ! I/O status
real(r8)                :: leaf_rhonir(num_pft)    ! leaf NIR reflectance [0-1]
real(r8)                :: leaf_taunir(num_pft)    ! leaf NIR transmittance [0-1]
real(r8)                :: leaf_rhovis(num_pft)    ! leaf visible reflectance [0-1]
real(r8)                :: leaf_tauvis(num_pft)    ! leaf visible transmittance [0-1]
real(r8)                :: stem_rhonir(num_pft)    ! stem NIR reflectance [0-1]
real(r8)                :: stem_taunir(num_pft)    ! stem NIR transmittance [0-1]
real(r8)                :: stem_rhovis(num_pft)    ! stem visible reflectance [0-1]
real(r8)                :: stem_tauvis(num_pft)    ! stem visible transmittance [0-1]

! Namelist of pft-specific parameters
namelist /params/ clumping_index, xl, leaf_rhonir, leaf_taunir, leaf_rhovis,   &
  leaf_tauvis, stem_rhonir, stem_taunir, stem_rhovis, stem_tauvis

! Now read parameters namelist
pft_nl_file = open_file(trim(pft_nl), 'r')
read(pft_nl_file, params, iostat=ios, iomsg=msg)

if (ios /= 0) then
  ! Problem reading file - tell user.
  write(logf, '(A, I6, A, A)') "Error reading pft namelist file",              &
    ios, " IOMSG: ", msg
  stop "Stopped"
end if

! combine arrays
rhol(:,1) = leaf_rhovis
rhol(:,2) = leaf_rhonir

close(pft_nl_file)

end subroutine read_pft_namelist

!:.............................................................................:

subroutine init_pft_data(num_pft, numSWB)
  !
  ! DESCRIPTION:
  !		read in the namelist associated with the pft-specific parameters and 
  !   initialize data
  !

  use EDPftvarcon,              only : EDPftvarcon_inst
  use FatesParametersInterface, only : fates_parameters_type
  use FatesConstantsMod,        only : r8 => fates_r8

  implicit none

  ! ARGUMENTS:
  integer, intent(in) :: num_pft ! number of pfts
  integer, intent(in) :: numSWB  ! number of shortwave bands to simulate

  ! LOCALS:
  class(fates_parameters_type), allocatable :: fates_params            ! FATES input parameters
  real(r8)                                  :: clumping_index(num_pft) ! clumping index
  real(r8)                                  :: xl(num_pft)             ! leaf-stem orientation index
  real(r8)                                  :: rhol(num_pft,numSWB)    ! leaf reflectance [0-1]
  real(r8)                                  :: rhos(num_pft,numSWB)    ! stem reflectance [0-1]
  real(r8)                                  :: taul(num_pft,numSWB)    ! leaf transmittance [0-1]
  real(r8)                                  :: taus(num_pft,numSWB)    ! stem transmittance [0-1]

  ! initialize EDPFTvarcon_inst and fates_params
  call EDPftvarcon_inst%Init()
  allocate(fates_params)
  call fates_params%Init()
  call EDPftvarcon_inst%Register(fates_params)

  ! read in parameter values
  call read_pft_namelist(num_pft, numSWB, clumping_index, xl, rhol)

  ! set values
  ! TODO: make this read in from a parameter file, this is hacky
  EDPftvarcon_inst%clumping_index = clumping_index
  EDPftvarcon_inst%xl = xl
  EDPftvarcon_inst%rhol = rhol
  EDPftvarcon_inst%rhos = rhos
  EDPftvarcon_inst%taul = taul
  EDPftvarcon_inst%taus = taus

end subroutine init_pft_data



