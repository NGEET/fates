program FatesTestEdgeForest

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesEdgeForestMod,          only : gffeb_norm, gffeb_quadratic, is_param_set
  use EDParamsMod, only : ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma,ED_val_edgeforest_gaussian_center
  use EDParamsMod, only : ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma,ED_val_edgeforest_lognormal_center
  use EDParamsMod, only : ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b,ED_val_edgeforest_quadratic_c
  use EDParamsMod, only : ED_val_edgeforest_bin_edges
  use FatesUnitTestParamReaderMod, only :fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use shr_infnan_mod,              only : nan => shr_infnan_nan, assignment(=), isnan => shr_infnan_isnan

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader      ! param reader instance
  character(len=:), allocatable      :: param_file        ! input parameter file
  integer                            :: e, i              ! looping indices
  integer                            :: n_bins            ! number of edge forest bins in the parameter file
  integer                            :: n_pct_nonforest   ! size of pct_nonforest array
  real(r8), allocatable              :: pct_nonforest(:)  ! percent nonforest in site
  real(r8), allocatable              :: frac_in_bin_gaussian(:)  ! output: fraction of forest in a bin with Gaussian fit
  real(r8), allocatable              :: frac_in_bin_lognormal(:) ! output: fraction of forest in a bin with lognormal fit
  real(r8), allocatable              :: frac_in_bin_quadratic(:) ! output: fraction of forest in a bin with quadratic fit
  ! Edge bin parameters
  real(r8) :: amplitude, mu, sigma, a, b, c
  logical  :: bin_found

  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'edge_forest_out.nc'    ! output file
  real(r8),         parameter :: min_pct_nonforest = 0._r8      ! minimum % nonforest to calculate
  real(r8),         parameter :: max_pct_nonforest = 100.0_r8   ! maximum % nonforest to calculate
  real(r8),         parameter :: pct_nonforest_inc = 0.1_r8     ! % nonforest increment to use

  interface

    subroutine WriteEdgeForestData(out_file, n_pct_nonforest, pct_nonforest, frac_in_bin_gaussian, &
        frac_in_bin_lognormal, frac_in_bin_quadratic)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n_pct_nonforest
      real(r8),         intent(in) :: pct_nonforest(:)
      real(r8),         intent(in) :: frac_in_bin_gaussian(:)
      real(r8),         intent(in) :: frac_in_bin_lognormal(:)
      real(r8),         intent(in) :: frac_in_bin_quadratic(:)
    end subroutine WriteEdgeForestData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()

  ! determine sizes of arrays
  n_bins = size(ED_val_edgeforest_bin_edges, dim=1)
  n_pct_nonforest = int((max_pct_nonforest - min_pct_nonforest)/pct_nonforest_inc + 1)

  ! allocate arrays
  allocate(pct_nonforest(n_pct_nonforest))
  allocate(frac_in_bin_gaussian(n_pct_nonforest))
  allocate(frac_in_bin_lognormal(n_pct_nonforest))
  allocate(frac_in_bin_quadratic(n_pct_nonforest))

  ! initialize pct_nonforest array
  do i = 1, n_pct_nonforest
    pct_nonforest(i) = min_pct_nonforest + pct_nonforest_inc*(i-1)
  end do

  ! calculate fraction using parameters of first Gaussian bin, if any
  bin_found = .false.
  do e = 1, n_bins
    if (is_param_set(ED_val_edgeforest_gaussian_amplitude(e))) then
      bin_found = .true.
      amplitude = ED_val_edgeforest_gaussian_amplitude(e)
      mu = ED_val_edgeforest_gaussian_center(e)
      sigma = ED_val_edgeforest_gaussian_sigma(e)
      write(*, '(a, i2, a)') "Gaussian (bin ",e,"):"
      write(*, '(a, E15.6)') "   amplitude: ",amplitude
      write(*, '(a, E15.6)') "          mu: ",mu
      write(*, '(a, E15.6)') "       sigma: ",sigma
      exit
    end if
  end do
  do i = 1, n_pct_nonforest
    if (bin_found) then
      if (pct_nonforest(i) == 0._r8) then
        frac_in_bin_gaussian(i) = 0._r8
      else
        frac_in_bin_gaussian(i) = gffeb_norm(pct_nonforest(i), amplitude, mu, sigma, lognorm=.false.)
      end if
    else
      frac_in_bin_gaussian(i) = nan
    end if
  end do

  ! calculate fraction using parameters of first lognormal bin, if any
  bin_found = .false.
  do e = 1, n_bins
    if (is_param_set(ED_val_edgeforest_lognormal_amplitude(e))) then
      bin_found = .true.
      amplitude = ED_val_edgeforest_lognormal_amplitude(e)
      mu = ED_val_edgeforest_lognormal_center(e)
      sigma = ED_val_edgeforest_lognormal_sigma(e)
      write(*, '(a, i2, a)') "Lognormal (bin ",e,"):"
      write(*, '(a, E15.6)') "   amplitude: ",amplitude
      write(*, '(a, E15.6)') "          mu: ",mu
      write(*, '(a, E15.6)') "       sigma: ",sigma
      exit
    end if
  end do
  do i = 1, n_pct_nonforest
    if (bin_found) then
      if (pct_nonforest(i) == 0._r8) then
        frac_in_bin_lognormal(i) = 0._r8
      else
        frac_in_bin_lognormal(i) = gffeb_norm(pct_nonforest(i), amplitude, mu, sigma, lognorm=.true.)
      end if
    else
      frac_in_bin_lognormal(i) = nan
    end if
  end do

  ! calculate fraction using parameters of first quadratic bin, if any
  bin_found = .false.
  do e = 1, n_bins
    if (is_param_set(ED_val_edgeforest_quadratic_a(e))) then
      bin_found = .true.
      a = ED_val_edgeforest_quadratic_a(e)
      b = ED_val_edgeforest_quadratic_b(e)
      c = ED_val_edgeforest_quadratic_c(e)
      write(*, '(a, i2, a)') "Quadratic (bin ",e,"):"
      write(*, '(a, E15.6)') "   a: ",amplitude
      write(*, '(a, E15.6)') "   b: ",mu
      write(*, '(a, E15.6)') "   c: ",sigma
      exit
    end if
  end do
  do i = 1, n_pct_nonforest
    if (bin_found) then
      if (pct_nonforest(i) == 0._r8) then
        frac_in_bin_quadratic(i) = 0._r8
      else
        frac_in_bin_quadratic(i) = gffeb_quadratic(pct_nonforest(i), a, b, c)
      end if
    else
      frac_in_bin_quadratic(i) = nan
    end if
  end do

  ! write out data to netcdf file
  call WriteEdgeForestData(out_file, n_pct_nonforest, pct_nonforest, frac_in_bin_gaussian, &
    frac_in_bin_lognormal, frac_in_bin_quadratic)

end program FatesTestEdgeForest

! ----------------------------------------------------------------------------------------

subroutine WriteEdgeForestData(out_file, n_pct_nonforest, pct_nonforest, frac_in_bin_gaussian, &
        frac_in_bin_lognormal, frac_in_bin_quadratic)
  !
  ! DESCRIPTION:
  ! Writes out data from the edge forest test
  !
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int
  use FatesConstantsMod,  only : r8 => fates_r8
  implicit none

  ! ARGUMENTS
  character(len=*), intent(in) :: out_file
  integer,          intent(in) :: n_pct_nonforest
  real(r8),         intent(in) :: pct_nonforest(:)
  real(r8),         intent(in) :: frac_in_bin_gaussian(:)
  real(r8),         intent(in) :: frac_in_bin_lognormal(:)
  real(r8),         intent(in) :: frac_in_bin_quadratic(:)

  ! LOCALS:
  ! TODO: Set a local parameter for number of dimensions (1 here; was 2 for Allometry)
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=24)     :: dim_names(1)   ! dimension name(s)
  integer              :: dimIDs(1)      ! dimension ID(s)
  integer              :: pctnonforestID ! variable ID(s) for dimensions
  integer              :: gaussianID, lognormalID, quadraticID

  ! dimension name(s)
  dim_names = [character(len=24) :: 'pct_nonforest']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n_pct_nonforest/), 1, dimIDs)

  ! register percent nonforest
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: '%', 'Percentage of site that is nonforest'], 2, pctnonforestID)

  ! register frac_in_bin_gaussian
  call RegisterVar(ncid, 'frac_in_bin_gaussian', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pct_nonforest', 'unitless', 'Fraction of forest in first bin with Gaussian fit'],  &
    3, gaussianID)

  ! register frac_in_bin_lognormal
  call RegisterVar(ncid, 'frac_in_bin_lognormal', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pct_nonforest', 'unitless', 'Fraction of forest in first bin with lognormal fit'],  &
    3, lognormalID)

  ! register frac_in_bin_quadratic
  call RegisterVar(ncid, 'frac_in_bin_quadratic', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'pct_nonforest', 'unitless', 'Fraction of forest in first bin with quadratic fit'],  &
    3, quadraticID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, pctnonforestID, pct_nonforest(:))
  call WriteVar(ncid, gaussianID, frac_in_bin_gaussian(:))
  call WriteVar(ncid, lognormalID, frac_in_bin_lognormal(:))
  call WriteVar(ncid, quadraticID, frac_in_bin_quadratic(:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteEdgeForestData