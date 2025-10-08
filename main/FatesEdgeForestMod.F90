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
  use FatesEcotypesMod, only : IsPatchForest
  use FatesUtilsMod,    only : is_param_set

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: CalculateTreeGrassAreaSite
  public :: CalculateEdgeForestArea
  public :: ApplyEdgeForestFlamToSite
  ! Public for unit testing
  public :: indexx
  public :: GetFracEdgeForestInEachBin
  public :: GetFracEdgeForestInEachBin_norm_numerator
  public :: GetFracEdgeForestInEachBin_norm_denominator
  public :: GetFracEdgeForestInEachBin_norm
  public :: GetFracEdgeForestInEachBin_quadratic
  public :: AssignPatchToBins
  public :: CalcEdgeForestFlam_1var_1bin
  public :: CalcEdgeForestFlam_1var
  public :: ApplyEdgeForestFlamToPatch_1var
  public :: CheckFlamChangeIntended

contains

  ! =====================================================================================

  subroutine CalculateTreeGrassAreaSite(csite, tree_fraction, grass_fraction, bare_fraction)
    !
    !  DESCRIPTION:
    !  Calculates total grass, tree, and bare fractions for a site

    use FatesEcotypesMod,   only : IsPatchForest
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

    currentPatch => csite%oldest_patch
    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        call currentPatch%UpdateTreeGrassArea()
        tree_fraction = tree_fraction + currentPatch%total_tree_area/AREA
        grass_fraction = grass_fraction + currentPatch%total_grass_area/AREA
        currentPatch%is_forest = IsPatchForest(currentPatch, forest_tree_fraction_threshold)
      end if
      currentPatch => currentPatch%younger
    end do

    ! if cover > 1.0, grasses are under the trees
    grass_fraction = min(grass_fraction, 1.0_r8 - tree_fraction)
    bare_fraction = 1.0_r8 - tree_fraction - grass_fraction

    ! Must come after patch loop with IsPatchForest() call
    call CalculateEdgeForestArea(csite)

  end subroutine CalculateTreeGrassAreaSite


  subroutine GetNumberOfForestPatches(site, n_forest_patches, area_site)
    ! DESCRIPTION
    ! Returns number and area of forest patches at site
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    integer, intent(out) :: n_forest_patches
    real(r8), intent(out) :: area_site
    !
    ! LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch

    n_forest_patches = 0
    site%area_forest_patches = 0._r8
    area_site = 0._r8
    currentPatch => site%youngest_patch
    do while(associated(currentPatch))
       area_site = area_site + currentPatch%area
       if (currentPatch%is_forest) then
          n_forest_patches = n_forest_patches + 1
          site%area_forest_patches = site%area_forest_patches + currentPatch%area
       end if

       currentPatch => currentPatch%older
    enddo

  end subroutine GetNumberOfForestPatches


  subroutine RankForestEdgeProximity(site, indices, index_forestpatches_to_allpatches)
    ! DESCRIPTION:
    ! Rank forest patches by their proximity to edge, using age as a proxy.
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

    ! Skip sites with no forest patches
    n_forest_patches = size(indices)
    if (n_forest_patches == 0) then
       return
    end if

    ! Allocate arrays
    allocate(array(1:n_forest_patches))

    ! Fill arrays
    f = 0
    p = 0
    index_forestpatches_to_allpatches(:) = 0
    currentPatch => site%oldest_patch
    patchloop: do while(associated(currentPatch))
       p = p + 1
       if (.not. currentPatch%is_forest) then
          currentPatch => currentPatch%younger
          cycle
       end if

       f = f + 1
       index_forestpatches_to_allpatches(p) = f

       ! Fill with patch age.
       ! TODO: Add other options. Biomass? Woody biomass?
       array(f) = currentPatch%age

       currentPatch => currentPatch%younger
    end do patchloop

    ! Get indices of sorted forest patches
    call indexx(array, indices)

    ! Clean up
    deallocate(array)
  end subroutine RankForestEdgeProximity

  function GetFracEdgeForestInEachBin_norm_numerator(x_in, A, mu, sigma, lognorm)
    ! DESCRIPTION
    ! Gets numerator at xof normal-like function (Gaussian if lognorm==.true., lognormal otherwise)
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: x_in
    real(r8), intent(in) :: A      ! Amplitude
    real(r8), intent(in) :: mu     ! Center
    real(r8), intent(in) :: sigma  ! Width
    logical,  intent(in) :: lognorm  ! Gaussian function if true, lognormal otherwise
    !
    ! RETURN VALUE:
    real(r8) :: GetFracEdgeForestInEachBin_norm_numerator
    !
    ! LOCAL VARIABLES:
    real(r8) :: x  ! either x_in or its log

    if (lognorm) then
       x = log(x_in)
    else
       x = x_in
    end if

    GetFracEdgeForestInEachBin_norm_numerator = A * exp(-(x - mu)**2 / (2*sigma**2))
  end function GetFracEdgeForestInEachBin_norm_numerator

  function GetFracEdgeForestInEachBin_norm_denominator(x, sigma, lognorm)
    ! DESCRIPTION
    ! Gets denominator at x of normal-like function (Gaussian if lognorm==.true., lognormal otherwise)
    !
    ! ARGUMENTS:
    use FatesConstantsMod, only : pi => pi_const
    real(r8), intent(in) :: x
    real(r8), intent(in) :: sigma  ! Width
    logical,  intent(in) :: lognorm  ! Gaussian function if true, lognormal otherwise
    !
    ! RETURN VALUE:
    real(r8) :: GetFracEdgeForestInEachBin_norm_denominator

    GetFracEdgeForestInEachBin_norm_denominator = sigma * sqrt(2*pi)
    if (lognorm) then
       GetFracEdgeForestInEachBin_norm_denominator = GetFracEdgeForestInEachBin_norm_denominator * x
    end if
  end function GetFracEdgeForestInEachBin_norm_denominator

  function GetFracEdgeForestInEachBin_norm(x, A, mu, sigma, lognorm)
    ! DESCRIPTION
    ! Gets value at x of normal-like function (Gaussian if lognorm==.true., lognormal otherwise)
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: x
    real(r8), intent(in) :: A      ! Amplitude
    real(r8), intent(in) :: mu     ! Center
    real(r8), intent(in) :: sigma  ! Width
    logical,  intent(in) :: lognorm  ! Gaussian function if true, lognormal otherwise
    !
    ! RETURN VALUE:
    real(r8) :: GetFracEdgeForestInEachBin_norm

    GetFracEdgeForestInEachBin_norm = GetFracEdgeForestInEachBin_norm_numerator(x, A, mu, sigma, lognorm) / GetFracEdgeForestInEachBin_norm_denominator(x, sigma, lognorm)
  end function GetFracEdgeForestInEachBin_norm

  function GetFracEdgeForestInEachBin_quadratic(x, a, b, c)
    ! DESCRIPTION
    ! Gets value at x of quadratic function
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: x
    real(r8), intent(in) :: a, b, c  ! Parameters
    !
    ! RETURN VALUE:
    real(r8) :: GetFracEdgeForestInEachBin_quadratic

    GetFracEdgeForestInEachBin_quadratic = a*(x**2) + b*x + c
  end function GetFracEdgeForestInEachBin_quadratic

  subroutine GetFracEdgeForestInEachBin(x, nlevedgeforest, efb_gaussian_amplitudes, efb_gaussian_sigmas, efb_gaussian_centers, efb_lognormal_amplitudes, efb_lognormal_sigmas, efb_lognormal_centers, efb_quadratic_a, efb_quadratic_b, efb_quadratic_c, fraction_forest_in_bin, norm)
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
    logical :: lognorm
    logical :: do_norm
    ! Error checking
    real(r8), parameter :: tol = 1.e-9_r8  ! fraction of total forest area
    real(r8) :: err_chk

    ! Initialize
    fraction_forest_in_bin(:) = 0._r8

    ! If the cell is nearly 0% forest, put any forest in the first edge bin (closest to edge)
    if (x < nearzero) then
       fraction_forest_in_bin(1) = 1._r8
       return
    end if

    ! If the cell is (nearly) 100% forest, it's all "deep forest"
    if (1._r8 - x < nearzero) then
       fraction_forest_in_bin(nlevedgeforest) = 1._r8
       return
    end if

    if (present(norm)) then
       do_norm = norm
    else
       do_norm = .true.
    end if

    binloop: do b = 1, nlevedgeforest

       if (is_param_set(efb_gaussian_amplitudes(b)) .or. is_param_set(efb_lognormal_amplitudes(b))) then
         ! Gaussian or Lognormal
         lognorm = is_param_set(efb_lognormal_amplitudes(b))
         if (lognorm) then
            A = efb_lognormal_amplitudes(b)
            mu = efb_lognormal_centers(b)
            sigma = efb_lognormal_sigmas(b)
         else
            A = efb_gaussian_amplitudes(b)
            mu = efb_gaussian_centers(b)
            sigma = efb_gaussian_sigmas(b)
         end if
         fraction_forest_in_bin(b) = GetFracEdgeForestInEachBin_norm(x, A, mu, sigma, lognorm)

       else if (is_param_set(efb_quadratic_a(b))) then
         ! Quadratic
         fraction_forest_in_bin(b) = GetFracEdgeForestInEachBin_quadratic(x, efb_quadratic_a(b), efb_quadratic_b(b), efb_quadratic_c(b))

       else
         call endrun("Unrecognized bin fit type")
       end if
    end do binloop

    ! Account for fit errors by normalizing to 1
    if (do_norm) then
       fraction_forest_in_bin(:) = fraction_forest_in_bin(:) / sum(fraction_forest_in_bin)

       err_chk = sum(fraction_forest_in_bin) - 1._r8
       if (abs(err_chk) > tol) then
          write(fates_log(),*) "ERROR: bin fractions don't sum to 1;    actual minus expected = ",err_chk
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end if


  end subroutine GetFracEdgeForestInEachBin


  subroutine AssignPatchToBins(fraction_forest_in_each_bin, area_forest_patches, patch_area, nlevedgeforest, tol, sum_forest_bins_so_far_m2, area_in_edgeforest_bins_m2)
    ! DESCRIPTION
    ! Given one patch in a site, assign its area to edge bin(s).
    !
    ! ARGUMENTS
    real(r8), dimension(:), pointer, intent(in) :: fraction_forest_in_each_bin  ! Fraction of site's forest area in each edge bin
    real(r8), intent(in) :: area_forest_patches  ! Total forest area in the site
    real(r8), intent(in) :: patch_area  ! Area of this patch
    integer,  intent(in) :: nlevedgeforest  ! Number of edge forest bins, including "deep forest"
    real(r8), intent(in) :: tol  ! Tolerance for checking total area assigned to bins
    real(r8), intent(inout) :: sum_forest_bins_so_far_m2  ! How much of the site's forest area has been assigned?
    real(r8), dimension(:), intent(out) :: area_in_edgeforest_bins_m2  ! Area of this patch in each edge bin (m2)
    !
    ! LOCAL VARIABLES
    real(r8) :: remaining_to_assign_from_patch_m2  ! How much of this patch's area still needs to be assigned
    real(r8) :: remaining_to_assign_to_bin_m2  ! How much of a given bin's area still needs to be assigned
    integer  :: b
    ! For checks
    real(r8) :: err_chk

    area_in_edgeforest_bins_m2(:) = 0._r8
    remaining_to_assign_from_patch_m2 = patch_area
    binloop: do b = 1, nlevedgeforest

       ! How much area is left for this bin?
       remaining_to_assign_to_bin_m2 = sum(fraction_forest_in_each_bin(1:b))*area_forest_patches - sum_forest_bins_so_far_m2
       if (remaining_to_assign_to_bin_m2 <= 0) then
          cycle
       end if

       ! Assign area
       area_in_edgeforest_bins_m2(b) = min(remaining_to_assign_from_patch_m2, remaining_to_assign_to_bin_m2)
       remaining_to_assign_from_patch_m2 = remaining_to_assign_from_patch_m2 - area_in_edgeforest_bins_m2(b)

       ! Update accounting
       sum_forest_bins_so_far_m2 = sum_forest_bins_so_far_m2 + area_in_edgeforest_bins_m2(b)

       if (remaining_to_assign_from_patch_m2 == 0._r8) then
          exit
       end if
    end do binloop

    ! Check that this patch's complete area was assigned (and no more)
    err_chk = remaining_to_assign_from_patch_m2
    if (abs(err_chk) > tol) then
       write(fates_log(),*) "ERROR: not enough or too much patch area was assigned to bins (check 1); remainder = ",err_chk
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    err_chk = patch_area - sum(area_in_edgeforest_bins_m2)
    if (abs(err_chk) > tol) then
       write(fates_log(),*) "ERROR: not enough or too much patch area was assigned to bins (check 2); remainder = ",err_chk
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
  end subroutine AssignPatchToBins


  subroutine AssignPatchesToBins(site, indices, index_forestpatches_to_allpatches, n_forest_patches, n_patches)
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
    integer, intent(in) :: n_forest_patches  ! Number of forest patches
    integer, intent(in) :: n_patches  ! Number of patches in site
    !
    ! LOCAL VARIABLES
    integer :: f, i, p, b
    type(fates_patch_type), pointer  :: currentPatch
    real(r8) :: sum_forest_bins_so_far_m2
    ! For checks
    real(r8), dimension(nlevedgeforest) :: bin_area_sums
    real(r8), parameter :: tol = 1.e-9_r8  ! m2
    real(r8) :: err_chk

    sum_forest_bins_so_far_m2 = 0._r8
    forestpatchloop: do f = 1, n_forest_patches

       ! Get the i'th patch (which is the f'th forest patch)
       i = indices(f)
       currentPatch => site%oldest_patch
       allpatchloop: do p = 1, n_patches
          if (index_forestpatches_to_allpatches(p) == i) then
             exit
          end if
          currentPatch => currentPatch%younger
       end do allpatchloop
       if ((.not. associated(currentPatch)) .and. (.not. p == n_patches)) then
          write(fates_log(),*) "ERROR: i'th patch (f'th forest patch) not found."
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Make sure this is a forest patch
       if (.not. currentPatch%is_forest) then
          write(fates_log(),*) "ERROR: unexpected non-forest patch in forestpatchloop"
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Assign this patch's area
       call AssignPatchToBins(site%fraction_forest_in_each_bin, site%area_forest_patches, currentPatch%area, nlevedgeforest, tol, sum_forest_bins_so_far_m2, currentPatch%area_in_edgeforest_bins)

    end do forestpatchloop

    ! More checks
    bin_area_sums(:) = 0._r8
    currentPatch => site%oldest_patch
    allpatchloop_check: do while (associated(currentPatch))
       if (currentPatch%is_forest) then

          ! Check that all area of each forest patch is assigned
          err_chk = sum(currentPatch%area_in_edgeforest_bins) - currentPatch%area
          if (abs(err_chk) > tol) then
             write(fates_log(),*) "ERROR: unexpected patch forest bin sum (check 3); actual minus expected = ",err_chk
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if

          ! Accumulate site-wide area in each bin
          binloop_check4a: do b = 1, nlevedgeforest
             bin_area_sums(b) = bin_area_sums(b) + currentPatch%area_in_edgeforest_bins(b) / site%area_forest_patches
          end do binloop_check4a

       end if
       currentPatch => currentPatch%younger
    end do allpatchloop_check
    ! Check that fraction in each bin is what was expected
    binloop_check4b: do b = 1, nlevedgeforest
       err_chk = bin_area_sums(b) - site%fraction_forest_in_each_bin(b)
       if (abs(err_chk) > tol) then
          write(fates_log(),*) "ERROR: unexpected bin sum (check 4); actual minus expected = ",err_chk
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do binloop_check4b
    ! Check that sum of all bin fractions is 1
    err_chk = 1._r8 - sum(site%fraction_forest_in_each_bin(:))
    if (abs(err_chk) > tol) then
       write(fates_log(),*) "ERROR: unexpected bin sum (check 5); actual minus expected = ",err_chk
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
  end subroutine AssignPatchesToBins


  subroutine CalculateEdgeForestArea(site)
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
    type(fates_patch_type), pointer  :: currentPatch
    integer, dimension(:), allocatable :: indices ! Indices to use if you want to sort patches
    integer, dimension(:), allocatable :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    integer :: n_forest_patches  ! Number of forest patches
    integer :: n_patches  ! Number of patches in site
    real(r8) :: area_site
    real(r8) :: frac_forest
    real(r8), dimension(nlevedgeforest), target :: fraction_forest_in_each_bin

    ! Zero out all fractions
    currentPatch => site%oldest_patch
    do while (associated(currentPatch))
      currentPatch%area_in_edgeforest_bins(:) = 0._r8
      currentPatch => currentPatch%younger
    end do

    ! Skip sites with no forest patches
    call GetNumberOfForestPatches(site, n_forest_patches, area_site)
    if (n_forest_patches == 0) then
       return
    end if

    ! Allocate arrays
    allocate(indices(1:n_forest_patches))
    n_patches = site%GetNumberOfPatches()
    allocate(index_forestpatches_to_allpatches(1:n_patches))

    ! Rank forest patches by their proximity to edge
    call RankForestEdgeProximity(site, indices, index_forestpatches_to_allpatches)

    ! Get fraction of forest area in each bin
    frac_forest = site%area_forest_patches / area_site
    call GetFracEdgeForestInEachBin(frac_forest, nlevedgeforest, ED_val_edgeforest_gaussian_amplitude, ED_val_edgeforest_gaussian_sigma, ED_val_edgeforest_gaussian_center, ED_val_edgeforest_lognormal_amplitude, ED_val_edgeforest_lognormal_sigma, ED_val_edgeforest_lognormal_center, ED_val_edgeforest_quadratic_a, ED_val_edgeforest_quadratic_b, ED_val_edgeforest_quadratic_c, fraction_forest_in_each_bin)
    site%fraction_forest_in_each_bin = fraction_forest_in_each_bin

    ! Assign patches to bins
    call AssignPatchesToBins(site, indices, index_forestpatches_to_allpatches, n_forest_patches, n_patches)

    ! Clean up
    deallocate(indices)
    deallocate(index_forestpatches_to_allpatches)
  end subroutine CalculateEdgeForestArea


  elemental function CalcEdgeForestFlam_1var_1bin(mult_factor, add_factor, weather_in) result(weather_out)
    ! DESCRIPTION:
    ! Apply flammability enhancements to one fireWeather variable for one edge bin
    !
    ! USES:
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: mult_factor  ! Multiplicative factor
    real(r8), intent(in) :: add_factor   ! Additive factor
    real(r8), intent(in) :: weather_in   ! Value of weather variable before applying factors
    !
    ! LOCAL VARIABLES:
    real(r8) :: weather_out   ! Value of weather variable after applying factors

    weather_out = (weather_in * mult_factor) + add_factor

  end function CalcEdgeForestFlam_1var_1bin


  subroutine CalcEdgeForestFlam_1var(mult_factors, add_factors, weather_in, weather_out)
    ! DESCRIPTION:
    ! Calculate one fireWeather variable for all edge bins after applying flammability enhancements
    !
    ! USES:
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: mult_factors(:)  ! Multiplicative factors
    real(r8), intent(in) :: add_factors(:)   ! Additive factors
    real(r8), intent(in) :: weather_in   ! Value of weather variable before applying factors
    real(r8), intent(out) :: weather_out(:)  ! Value of weather variable in each bin after applying factors
    !
    ! LOCAL VARIABLES:

    weather_out = CalcEdgeForestFlam_1var_1bin(mult_factors, add_factors, weather_in)

  end subroutine CalcEdgeForestFlam_1var


  subroutine ApplyEdgeForestFlamToPatch_1var(weather_by_edge_bin, patch_area_each_edge_bin, weather_inout)
    ! DESCRIPTION:
    ! Apply enhancements to one fireWeather variable in a patch based on how much of its area is in
    ! each edge bin
    !
    ! USES:
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: weather_by_edge_bin(:)  ! Weather value in each edge bin
    real(r8), intent(in) :: patch_area_each_edge_bin(:)  ! Patch area in each edge bin (unit doesn't matter)
    real(r8), intent(inout) :: weather_inout  ! Value of fireWeather variable in this patch
    !
    ! LOCAL VARIABLES:
    real(r8) :: patch_forest_area  ! Forest area in patch (unit doesn't matter)
    real(r8), allocatable :: patch_weight_each_edge_bin(:)  ! Area weighting of each edge bin in patch

    ! If patch has no or little forest, return early to avoid divide-by-zero
    ! TODO: Or should such patches get the same flammability as nearest edge? It doesn't make much
    ! sense for grassland to have, e.g., lower wind speed than nearest-edge forest just because the
    ! grassland isn't getting the flammability enhancement.
    patch_forest_area = sum(patch_area_each_edge_bin)
    if (patch_forest_area < nearzero) then
      return
    end if

    ! Calculate weight of each edge bin for this patch
    allocate(patch_weight_each_edge_bin(size(patch_area_each_edge_bin)))
    patch_weight_each_edge_bin = patch_area_each_edge_bin(:) / patch_forest_area

    weather_inout = sum(weather_by_edge_bin * patch_weight_each_edge_bin)

    ! Clean up
    deallocate(patch_weight_each_edge_bin)

  end subroutine ApplyEdgeForestFlamToPatch_1var


  function CheckFlamChangeIntended(params_mult, params_add, weather_before, weather_after, tol) result(change_intended)
    ! DESCRIPTION:
    ! Check whether fire weather change was intended. If not but one was found, throw error.
    !
    !
    ! ARGUMENTS:
    real(r8), intent(in) :: params_mult(:)
    real(r8), intent(in) :: params_add(:)
    real(r8), intent(in) :: weather_before
    real(r8), intent(in) :: weather_after
    real(r8), intent(in) :: tol
    !
    ! LOCAL VARIABLES
    real(r8) :: weather_diff
    !
    ! RESULT:
    logical :: change_intended

    change_intended = any(abs(params_mult - 1._r8) > tol) .or. any(abs(params_add) > tol)
    if (.not. change_intended) then
      weather_diff = weather_after - weather_before
      if (abs(weather_diff) > tol) then
        write(fates_log(),*) 'No fire weather change intended but diff ',weather_diff
        call endrun(msg=errMsg(__FILE__, __LINE__))
      end if
    end if
  end function CheckFlamChangeIntended


  subroutine ApplyEdgeForestFlamToSite(site)
    ! DESCRIPTION:
    ! Apply enhancements to one fireWeather variable in a patch based on how much of its area is in
    ! each edge bin
    !
    ! USES:
    use FatesConstantsMod, only: t_water_freeze_k_1atm
    use FatesInterfaceTypesMod, only : nlevedgeforest
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_fireweather_rh_mult, ED_val_edgeforest_fireweather_rh_add
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_fireweather_temp_C_mult, ED_val_edgeforest_fireweather_temp_C_add
    use FatesEdgeForestParamsMod, only : ED_val_edgeforest_fireweather_wind_mult, ED_val_edgeforest_fireweather_wind_add
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! LOCAL VARIABLES:
    real(r8) :: rh, temp_C, wind
    real(r8) :: rh_by_edge_bin(nlevedgeforest)  ! rh value in each edge bin
    real(r8) :: temp_C_by_edge_bin(nlevedgeforest)  ! temp_C value in each edge bin
    real(r8) :: wind_by_edge_bin(nlevedgeforest)  ! wind value in each edge bin
    type(fates_patch_type), pointer :: currentPatch
    real(r8) :: weather_inout
    logical  :: change_intended
    real(r8) :: tol = 1.e-9_r8
    real(r8) :: weather_diff

    ! Get site-level weather values, assuming youngest patch has the same values as all others
    rh = site%youngest_patch%fireWeather%rh
    temp_C = site%youngest_patch%fireWeather%temp_C
    wind = site%youngest_patch%fireWeather%wind

    ! Calculate weather values for each edge bin
    call CalcEdgeForestFlam_1var( &
         ED_val_edgeforest_fireweather_rh_mult, &
         ED_val_edgeforest_fireweather_rh_add, &
         rh, &
         rh_by_edge_bin)
    call CalcEdgeForestFlam_1var( &
         ED_val_edgeforest_fireweather_temp_C_mult, &
         ED_val_edgeforest_fireweather_temp_C_add, &
         temp_C, &
         temp_C_by_edge_bin)
    call CalcEdgeForestFlam_1var( &
         ED_val_edgeforest_fireweather_wind_mult, &
         ED_val_edgeforest_fireweather_wind_add, &
         wind, &
         wind_by_edge_bin)

    ! Update patch values
    currentPatch => site%oldest_patch
    do while(associated(currentPatch))
      ! RH
      weather_inout = currentPatch%fireWeather%rh
      call ApplyEdgeForestFlamToPatch_1var( &
           rh_by_edge_bin, &
           currentPatch%area_in_edgeforest_bins, &
           weather_inout)
      change_intended = CheckFlamChangeIntended(ED_val_edgeforest_fireweather_rh_mult, &
           ED_val_edgeforest_fireweather_rh_add, currentPatch%fireWeather%rh, weather_inout, tol)
      if (change_intended) then
        ! RH can't be negative. In rare cases it can be supersaturated (>100%), so don't check that.
        currentPatch%fireWeather%rh = max(weather_inout, 0._r8)
      end if

      ! Temp
      weather_inout = currentPatch%fireWeather%temp_C
      call ApplyEdgeForestFlamToPatch_1var( &
           temp_C_by_edge_bin, &
           currentPatch%area_in_edgeforest_bins, &
           weather_inout)
      change_intended = CheckFlamChangeIntended(ED_val_edgeforest_fireweather_temp_C_mult, &
           ED_val_edgeforest_fireweather_temp_C_add, currentPatch%fireWeather%temp_C, weather_inout, tol)
      if (change_intended) then
        ! Temperature can't be below absolute zero
        currentPatch%fireWeather%temp_C = max(weather_inout, -t_water_freeze_k_1atm)
      end if

      ! Wind
      weather_inout = currentPatch%fireWeather%wind
      call ApplyEdgeForestFlamToPatch_1var( &
           wind_by_edge_bin, &
           currentPatch%area_in_edgeforest_bins, &
           weather_inout)
      change_intended = CheckFlamChangeIntended(ED_val_edgeforest_fireweather_wind_mult, &
           ED_val_edgeforest_fireweather_wind_add, currentPatch%fireWeather%wind, weather_inout, tol)
      if (change_intended) then
        ! Wind speed can't be negative
        currentPatch%fireWeather%wind = max(weather_inout, 0._r8)
      end if

      currentPatch => currentPatch%younger
    end do

  end subroutine ApplyEdgeForestFlamToSite


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

    if (size(array) /= size(index)) then
       write(fates_log(),*) 'ERROR: INDEXX size mismatch.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    else if (size(array) == 0) then
       write(fates_log(),*) 'ERROR: INDEXX array size 0.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do i=1,size(index)
       index(i)=i
    enddo

    call q_sort_index(array,index,1,size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers,index,left,right)

    !> This is the recursive subroutine actually used by sort_patches.
    !>
    !> This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    !> It is not derived from any NR code, but are based on a quicksort routine by
    !> Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !> in C, and issued under the GNU General Public License. The conversion to
    !> Fortran 90, and modification to do an index sort was done by Ian Rutt.

    implicit none

    real(r8), dimension(:) :: numbers !> Numbers being sorted
    integer, dimension(:) :: index   !> Returned index
    integer :: left, right           !> Limit of sort region

    integer :: ll,rr
    integer :: pv_int,l_hold, r_hold,pivpos
    real(r8) :: pivot

    ll=left
    rr=right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos=index(ll)

    do
       if (.not.(ll < rr)) exit

       do
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr=rr-1
       enddo

       if (ll /= rr) then
          index(ll) = index(rr)
          ll=ll+1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll=ll+1
       enddo

       if (ll /= rr) then
          index(rr) = index(ll)
          rr=rr-1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index


  end module FatesEdgeForestMod
