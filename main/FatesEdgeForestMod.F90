module FatesEdgeForestMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nocomp_bareground
  use FatesConstantsMod, only : nearzero
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use EDTypesMod, only : ed_site_type
  use EDTypesMod, only : AREA
  use FatesPatchMod, only : fates_patch_type
  use FatesEcotypesMod, only : is_patch_forest
  use FatesUtilsMod,    only : is_param_set

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: CalculateTreeGrassAreaSite
  public :: calculate_edgeforest_area
  ! Public for unit testing
  public :: indexx
  public :: get_fraction_of_edgeforest_in_each_bin
  public :: gffeb_norm_numerator
  public :: gffeb_norm_denominator
  public :: gffeb_norm
  public :: gffeb_quadratic
  public :: assign_patch_to_bins

contains

  ! =====================================================================================

  subroutine CalculateTreeGrassAreaSite(csite, tree_fraction, grass_fraction, bare_fraction)
    !
    !  DESCRIPTION:
    !  Calculates total grass, tree, and bare fractions for a site

    use FatesEcotypesMod,   only : is_patch_forest
    use EDParamsMod,        only : forest_tree_fraction_threshold

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: csite  ! site object
    real(r8),           intent(out)   :: tree_fraction  ! total site tree fraction
    real(r8),           intent(out)   :: grass_fraction ! total site grass fraction
    real(r8),           intent(out)   :: bare_fraction  ! total site bare fraction

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch ! patch object

    tree_fraction = 0.0_r8
    grass_fraction = 0.0_r8
    bare_fraction = 1.0_r8 - tree_fraction - grass_fraction

  end subroutine CalculateTreeGrassAreaSite


  subroutine get_number_of_forest_patches(site, n_forest_patches, area_forest_patches, area)
    ! DESCRIPTION
    ! Returns number and area of forest patches at site
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    integer, intent(out) :: n_forest_patches
    real(r8), intent(out) :: area_forest_patches
    real(r8), intent(out) :: area
    !
    ! LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch

    n_forest_patches = 0
    area_forest_patches = 0._r8
    area = 0._r8

  end subroutine get_number_of_forest_patches


  function get_number_of_patches(site) result(n_patches)
    ! DESCRIPTION
    ! Returns number of patches at site
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! RETURN VALUE:
    integer :: n_patches
    !
    ! LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch

    n_patches = 0

  end function get_number_of_patches


  subroutine rank_forest_edge_proximity(site, indices, index_forestpatches_to_allpatches)
    ! DESCRIPTION:
    ! Rank forest patches by their proximity to edge.
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    integer, dimension(:), intent(inout) :: indices  ! Indices to use if you want to sort patches
    integer, dimension(:), intent(inout) :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    !
    ! LOCAL VARIABLES:
    real(r8), dimension(:), allocatable :: array ! Array to be index-sorted.
    integer :: n_forest_patches  ! Number of patches in above arrays
    type(fates_patch_type), pointer  :: currentPatch
    integer :: f  ! index of current forest patch
    integer :: p  ! index of patch

    return
  end subroutine rank_forest_edge_proximity

  function gffeb_norm_numerator(x_in, A, mu, sigma, lognorm)
    real(r8), intent(in) :: x_in
    real(r8), intent(in) :: A      ! Amplitude
    real(r8), intent(in) :: mu     ! Center
    real(r8), intent(in) :: sigma  ! Sigma
    logical,  intent(in) :: lognorm  ! Whether to take log(x)
    real(r8) :: x
    real(r8) :: gffeb_norm_numerator

    gffeb_norm_numerator = 0._r8
  end function gffeb_norm_numerator

  function gffeb_norm_denominator(x, sigma, lognorm)
    use FatesConstantsMod, only : pi => pi_const
    real(r8), intent(in) :: x
    real(r8), intent(in) :: sigma  ! Sigma
    logical,  intent(in) :: lognorm  ! Whether to take log(x)
    real(r8) :: gffeb_norm_denominator

    gffeb_norm_denominator = 1._r8
  end function gffeb_norm_denominator

  function gffeb_norm(x, A, mu, sigma, lognorm)
    real(r8), intent(in) :: x
    real(r8), intent(in) :: A      ! Amplitude
    real(r8), intent(in) :: mu     ! Center
    real(r8), intent(in) :: sigma  ! Sigma
    logical,  intent(in) :: lognorm  ! Whether to take log(x) in numerator
    real(r8) :: gffeb_norm

    gffeb_norm = 0._r8
  end function gffeb_norm

  function gffeb_quadratic(x, a, b, c)
    real(r8), intent(in) :: x
    real(r8), intent(in) :: a, b, c  ! Parameters
    real(r8) :: gffeb_quadratic

    gffeb_quadratic = 0._r8
  end function gffeb_quadratic

  subroutine get_fraction_of_edgeforest_in_each_bin(x, nlevedgeforest, efb_gaussian_amplitudes, efb_gaussian_sigmas, efb_gaussian_centers, efb_lognormal_amplitudes, efb_lognormal_sigmas, efb_lognormal_centers, efb_quadratic_a, efb_quadratic_b, efb_quadratic_c, fraction_forest_in_bin, norm)
    ! DESCRIPTION:
    ! Get the fraction of forest in each bin.
    !
    ! USES
    use FatesConstantsMod, only : pi => pi_const
    !
    ! ARGUMENTS
    real(r8), intent(in)  :: x  ! Independent variable in the fit
    integer, intent(in) :: nlevedgeforest
    real(r8), dimension(:), intent(in) :: efb_gaussian_amplitudes
    real(r8), dimension(:), intent(in) :: efb_gaussian_sigmas
    real(r8), dimension(:), intent(in) :: efb_gaussian_centers
    real(r8), dimension(:), intent(in) :: efb_lognormal_amplitudes
    real(r8), dimension(:), intent(in) :: efb_lognormal_sigmas
    real(r8), dimension(:), intent(in) :: efb_lognormal_centers
    real(r8), dimension(:), intent(in) :: efb_quadratic_a
    real(r8), dimension(:), intent(in) :: efb_quadratic_b
    real(r8), dimension(:), intent(in) :: efb_quadratic_c
    real(r8), dimension(:), intent(out) :: fraction_forest_in_bin
    logical, optional, intent(in) :: norm
    !
    ! LOCAL VARIABLES
    integer :: b  ! Bin index
    real(r8) :: A      ! Amplitude
    real(r8) :: mu     ! Center
    real(r8) :: sigma  ! Sigma
    ! Error checking
    real(r8), parameter :: tol = 1.e-9_r8  ! fraction of total forest area
    real(r8) :: err_chk
    logical :: lognorm
    logical :: do_norm

    fraction_forest_in_bin(:) = 0._r8

  end subroutine get_fraction_of_edgeforest_in_each_bin


  subroutine assign_patch_to_bins(fraction_forest_in_each_bin, area_forest_patches, patch_area, nlevedgeforest, tol, sum_forest_bins_so_far_m2, area_in_edgeforest_bins)
    ! DESCRIPTION
    ! Given one patch in a site, assign its area to edge bin(s).
    !
    ! ARGUMENTS
    real(r8), dimension(:), pointer, intent(in) :: fraction_forest_in_each_bin
    real(r8), intent(in) :: area_forest_patches
    real(r8), intent(in) :: patch_area
    integer,  intent(in) :: nlevedgeforest
    real(r8), intent(in) :: tol
    real(r8), intent(inout) :: sum_forest_bins_so_far_m2
    real(r8), dimension(:), intent(out) :: area_in_edgeforest_bins
    !
    ! LOCAL VARIABLES
    real(r8) :: remaining_to_assign_from_patch_m2
    real(r8) :: remaining_to_assign_to_bin_m2
    integer  :: b
    ! For checks
    real(r8) :: err_chk

    area_in_edgeforest_bins(:) = 0._r8
    sum_forest_bins_so_far_m2 = 0._r8
  end subroutine assign_patch_to_bins


  subroutine assign_patches_to_bins(site, indices, index_forestpatches_to_allpatches, fraction_forest_in_each_bin, n_forest_patches, n_patches, area_forest_patches)
    ! DESCRIPTION
    ! Loops through forest patches from nearest to farthest from edge, assigning their
    ! area to edge bin(s).
    !
    ! USES:
    use FatesInterfaceTypesMod, only : nlevedgeforest
    ! ARGUMENTS
    type(ed_site_type), pointer, intent(in) :: site
    integer, dimension(:), intent(in) :: indices  ! Indices to use if you want to sort patches
    integer, dimension(:), intent(in) :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    real(r8), dimension(:), pointer, intent(in) :: fraction_forest_in_each_bin
    integer, intent(in) :: n_forest_patches  ! Number of forest patches
    integer, intent(in) :: n_patches  ! Number of patches in site
    real(r8), intent(in) :: area_forest_patches
    !
    ! LOCAL VARIABLES
    integer :: f, i, p, b
    type(fates_patch_type), pointer  :: currentPatch
    real(r8) :: sum_forest_bins_so_far_m2
    ! For checks
    real(r8), dimension(nlevedgeforest) :: bin_area_sums
    real(r8), parameter :: tol = 1.e-9_r8  ! m2
    real(r8) :: err_chk
    
    return


  end subroutine assign_patches_to_bins


  subroutine calculate_edgeforest_area(site)
    ! DESCRIPTION:
    ! Loop through forest patches in decreasing order of proximity, calculating the
    ! area of each patch that is in each edge bin.
    !
    ! USES:
    use FatesInterfaceTypesMod, only : nlevedgeforest
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma, ED_val_edgeforest_gaussian_center
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma, ED_val_edgeforest_lognormal_center
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b, ED_val_edgeforest_quadratic_c
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_bin_edges
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! LOCAL VARIABLES:
    integer, dimension(:), allocatable :: indices ! Indices to use if you want to sort patches
    integer, dimension(:), allocatable :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    integer :: n_forest_patches  ! Number of forest patches
    integer :: n_patches  ! Number of patches in site
    real(r8) :: area_forest_patches
    real(r8) :: area
    real(r8) :: frac_forest
    real(r8), dimension(nlevedgeforest), target :: fraction_forest_in_each_bin

    return
  end subroutine calculate_edgeforest_area




  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array.
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx.
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array, index)

    ! Performs an index sort of \texttt{array} and returns the result in
    ! \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !
    ! This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    ! It is not derived from any NR code, but are based on a quicksort routine by
    ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    ! in C, and issued under the GNU General Public License. The conversion to
    ! Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(r8), dimension(:) :: array ! Array to be indexed.
    integer, dimension(:) :: index ! Index of elements of patch_array
    integer :: i

    index = 1

  end subroutine indexx

!==============================================================


  end module FatesEdgeForestMod