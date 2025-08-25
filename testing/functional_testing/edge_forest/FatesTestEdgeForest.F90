program FatesTestEdgeForest

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesUtilsMod,               only : is_param_set
  use FatesEdgeForestMod,          only : gffeb_norm, gffeb_quadratic
  use FatesEdgeForestMod,          only : get_fraction_of_edgeforest_in_each_bin
  use FatesEdgeForestParamsMod, only : ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma,ED_val_edgeforest_gaussian_center
  use FatesEdgeForestParamsMod, only : ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma,ED_val_edgeforest_lognormal_center
  use FatesEdgeForestParamsMod, only : ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b,ED_val_edgeforest_quadratic_c
  use FatesEdgeForestParamsMod, only : ED_val_edgeforest_bin_edges
  use FatesUnitTestParamReaderMod, only :fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use shr_infnan_mod,              only : nan => shr_infnan_nan, assignment(=), isnan => shr_infnan_isnan

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader      ! param reader instance
  character(len=:), allocatable      :: param_file        ! input parameter file
  integer                            :: e, i              ! looping indices
  integer                            :: n_bins            ! number of edge forest bins in the parameter file
  integer                            :: n_frac_forest   ! size of frac_forest array
  real(r8), allocatable              :: frac_forest(:)  ! fraction forest in site
  real(r8), allocatable              :: frac_in_bin_gaussian(:)  ! output: fraction of forest in a bin with Gaussian fit
  real(r8), allocatable              :: frac_in_bin_lognormal(:) ! output: fraction of forest in a bin with lognormal fit
  real(r8), allocatable              :: frac_in_bin_quadratic(:) ! output: fraction of forest in a bin with quadratic fit
  real(r8), allocatable              :: frac_in_every_bin(:,:) ! output: fraction of forest in each bin
  real(r8), allocatable              :: frac_in_every_bin_norm(:,:) ! output: fraction of forest in each bin (normalized)
  ! Edge bin parameters
  real(r8) :: amplitude, mu, sigma, a, b, c
  logical  :: bin_found

  ! CONSTANTS:
  logical,          parameter :: debug = .false.
  character(len=*), parameter :: out_file = 'edge_forest_out.nc'    ! output file
  real(r8),         parameter :: min_frac_forest = 0._r8      ! minimum fraction forest to calculate
  real(r8),         parameter :: max_frac_forest = 1.0_r8   ! maximum fraction forest to calculate
  real(r8),         parameter :: frac_forest_inc = 0.001_r8     ! fraction forest increment to use

  interface

    subroutine WriteEdgeForestData(out_file, n_frac_forest, n_bins, frac_forest, frac_in_bin_gaussian, &
        frac_in_bin_lognormal, frac_in_bin_quadratic, frac_in_every_bin, frac_in_every_bin_norm)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n_frac_forest
      integer,          intent(in) :: n_bins
      real(r8),         intent(in) :: frac_forest(:)
      real(r8),         intent(in) :: frac_in_bin_gaussian(:)
      real(r8),         intent(in) :: frac_in_bin_lognormal(:)
      real(r8),         intent(in) :: frac_in_bin_quadratic(:)
      real(r8),         intent(in) :: frac_in_every_bin(:,:)
      real(r8),         intent(in) :: frac_in_every_bin_norm(:,:)
    end subroutine WriteEdgeForestData

  end interface

  ! read in parameter file name from command line
  param_file = command_line_arg(1)

  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()

  ! determine sizes of arrays
  n_bins = size(ED_val_edgeforest_bin_edges, dim=1)
  n_frac_forest = int((max_frac_forest - min_frac_forest)/frac_forest_inc + 1)

  ! allocate arrays
  allocate(frac_forest(n_frac_forest))
  allocate(frac_in_bin_gaussian(n_frac_forest))
  allocate(frac_in_bin_lognormal(n_frac_forest))
  allocate(frac_in_bin_quadratic(n_frac_forest))
  allocate(frac_in_every_bin(n_frac_forest, n_bins))
  allocate(frac_in_every_bin_norm(n_frac_forest, n_bins))

  ! initialize frac_forest array
  do i = 1, n_frac_forest
    frac_forest(i) = min_frac_forest + frac_forest_inc*(i-1)
  end do

  ! if debugging, print read-in parameters for all bins
  if (debug) then
    do e = 1, n_bins

      ! Gaussian
      if (is_param_set(ED_val_edgeforest_gaussian_amplitude(e))) then
        amplitude = ED_val_edgeforest_gaussian_amplitude(e)
        mu = ED_val_edgeforest_gaussian_center(e)
        sigma = ED_val_edgeforest_gaussian_sigma(e)
        write(*, '(a, i2, a)') "Gaussian (bin ",e,"):"
        write(*, '(a, E15.6)') "   amplitude: ",amplitude
        write(*, '(a, E15.6)') "          mu: ",mu
        write(*, '(a, E15.6)') "       sigma: ",sigma

      ! Lognormal
      else if (is_param_set(ED_val_edgeforest_lognormal_amplitude(e))) then
        bin_found = .true.
        amplitude = ED_val_edgeforest_lognormal_amplitude(e)
        mu = ED_val_edgeforest_lognormal_center(e)
        sigma = ED_val_edgeforest_lognormal_sigma(e)
        write(*, '(a, i2, a)') "Lognormal (bin ",e,"):"
        write(*, '(a, E15.6)') "   amplitude: ",amplitude
        write(*, '(a, E15.6)') "          mu: ",mu
        write(*, '(a, E15.6)') "       sigma: ",sigma

      ! Quadratic
      else if (is_param_set(ED_val_edgeforest_quadratic_a(e))) then
        bin_found = .true.
        a = ED_val_edgeforest_quadratic_a(e)
        b = ED_val_edgeforest_quadratic_b(e)
        c = ED_val_edgeforest_quadratic_c(e)
        write(*, '(a, i2, a)') "Quadratic (bin ",e,"):"
        write(*, '(a, E15.6)') "   a: ",a
        write(*, '(a, E15.6)') "   b: ",b
        write(*, '(a, E15.6)') "   c: ",c
      end if
    end do
  end if

  ! calculate fraction in every bin
  do i = 1, n_frac_forest

    ! Pre-normalization
    call get_fraction_of_edgeforest_in_each_bin(frac_forest(i), n_bins, &
        ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma, ED_val_edgeforest_gaussian_center, &
        ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma, ED_val_edgeforest_lognormal_center, &
        ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b, ED_val_edgeforest_quadratic_c, &
        frac_in_every_bin(i,:), .false.)

    ! Normalized
    call get_fraction_of_edgeforest_in_each_bin(frac_forest(i), n_bins, &
        ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma, ED_val_edgeforest_gaussian_center, &
        ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma, ED_val_edgeforest_lognormal_center, &
        ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b, ED_val_edgeforest_quadratic_c, &
        frac_in_every_bin_norm(i,:), .true.)
  end do

  ! write out data to netcdf file
  call WriteEdgeForestData(out_file, n_frac_forest, n_bins, frac_forest, frac_in_bin_gaussian, &
    frac_in_bin_lognormal, frac_in_bin_quadratic, frac_in_every_bin, frac_in_every_bin_norm)

end program FatesTestEdgeForest

! ----------------------------------------------------------------------------------------

subroutine WriteEdgeForestData(out_file, n_frac_forest, n_bins, frac_forest, frac_in_bin_gaussian, &
        frac_in_bin_lognormal, frac_in_bin_quadratic, frac_in_every_bin, frac_in_every_bin_norm)
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
  use FatesEdgeForestParamsMod, only : ED_val_edgeforest_bin_edges
  implicit none

  ! ARGUMENTS
  character(len=*), intent(in) :: out_file
  integer,          intent(in) :: n_frac_forest
  integer,          intent(in) :: n_bins
  real(r8),         intent(in) :: frac_forest(:)
  real(r8),         intent(in) :: frac_in_bin_gaussian(:)
  real(r8),         intent(in) :: frac_in_bin_lognormal(:)
  real(r8),         intent(in) :: frac_in_bin_quadratic(:)
  real(r8),         intent(in) :: frac_in_every_bin(:,:)
  real(r8),         intent(in) :: frac_in_every_bin_norm(:,:)

  ! LOCALS:
  ! TODO: Set a local parameter for number of dimensions (2)
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=24)     :: dim_names(2)   ! dimension name(s)
  integer              :: dimIDs(2)      ! dimension ID(s)
  integer              :: fracforestID, binID ! variable ID(s) for dimensions
  integer              :: gaussianID, lognormalID, quadraticID
  integer              :: everybinID, everybinnormID

  ! dimension name(s)
  dim_names = [character(len=24) :: 'frac_forest', 'bin']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n_frac_forest, n_bins/), 2, dimIDs)

  ! register fraction forest
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'unitless', 'Fraction of site that is forest'], 2, fracforestID)

  ! register bin
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'unitless', 'FATES edge bin min. distance to nonforest'], 2, binID)

  ! register frac_in_bin_gaussian
  call RegisterVar(ncid, 'frac_in_bin_gaussian', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'frac_forest', 'unitless', 'Fraction of forest in first bin with Gaussian fit'],  &
    3, gaussianID)

  ! register frac_in_bin_lognormal
  call RegisterVar(ncid, 'frac_in_bin_lognormal', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'frac_forest', 'unitless', 'Fraction of forest in first bin with lognormal fit'],  &
    3, lognormalID)

  ! register frac_in_bin_quadratic
  call RegisterVar(ncid, 'frac_in_bin_quadratic', dimIDs(1:1), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'frac_forest', 'unitless', 'Fraction of forest in first bin with quadratic fit'],  &
    3, quadraticID)

  ! register frac_in_every_bin
  call RegisterVar(ncid, 'frac_in_every_bin', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'frac_forest bin', 'unitless', 'Fraction of forest in every bin'],  &
    3, everybinID)

  ! register frac_in_every_bin_norm
  call RegisterVar(ncid, 'frac_in_every_bin_norm', dimIDs(1:2), type_double,   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'], &
    [character(len=150) :: 'frac_forest bin', 'unitless', 'Fraction of forest in every bin (normalized)'],  &
    3, everybinnormID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, fracforestID, frac_forest(:))
  call WriteVar(ncid, binID, ED_val_edgeforest_bin_edges(:))
  call WriteVar(ncid, gaussianID, frac_in_bin_gaussian(:))
  call WriteVar(ncid, lognormalID, frac_in_bin_lognormal(:))
  call WriteVar(ncid, quadraticID, frac_in_bin_quadratic(:))
  call WriteVar(ncid, everybinID, frac_in_every_bin(:,:))
  call WriteVar(ncid, everybinnormID, frac_in_every_bin_norm(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteEdgeForestData
