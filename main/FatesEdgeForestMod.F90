module FatesEdgeForestMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesInterfaceTypesMod, only : num_edge_forest_bins
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use EDTypesMod, only : ed_site_type
  use FatesPatchMod, only : fates_patch_type
  use FatesEcotypesMod, only : is_patch_forest

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: calculate_edge_area
  ! Public for unit testing
  public :: indexx
  public :: get_fraction_of_forest_in_each_bin
  public :: gffeb_lognorm_numerator
  public :: gffeb_lognorm_denominator

contains

  ! =====================================================================================

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
    currentPatch => site%youngest_patch
    do while(associated(currentPatch))
       area = area + currentPatch%area
       if (currentPatch%is_forest) then
          n_forest_patches = n_forest_patches + 1
          area_forest_patches = area_forest_patches + currentPatch%area
       end if

       currentPatch => currentPatch%older
    enddo

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
    currentPatch => site%youngest_patch
    do while(associated(currentPatch))
       n_patches = n_patches + 1
       currentPatch => currentPatch%older
    enddo

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
  end subroutine rank_forest_edge_proximity


  function gffeb_lognorm_numerator(x, A, mu, sigma)
    real(r8), intent(in) :: x
    real(r8), intent(in) :: A      ! Amplitude
    real(r8), intent(in) :: mu     ! Center
    real(r8), intent(in) :: sigma  ! Sigma
    real(r8) :: gffeb_lognorm_numerator

    gffeb_lognorm_numerator = A * exp(-(log(x) - mu)**2 / (2*sigma**2))
  end function gffeb_lognorm_numerator

  function gffeb_lognorm_denominator(x, sigma)
    use FatesConstantsMod, only : pi => pi_const
    real(r8), intent(in) :: x
    real(r8), intent(in) :: sigma  ! Sigma
    real(r8) :: gffeb_lognorm_denominator

    gffeb_lognorm_denominator = sigma * sqrt(2*pi) * x
  end function gffeb_lognorm_denominator

  subroutine get_fraction_of_forest_in_each_bin(x, num_edge_forest_bins, efb_amplitudes, efb_sigmas, efb_centers, efb_decay, fraction_forest_in_bin, norm)
    ! DESCRIPTION:
    ! Get the fraction of forest in each bin.
    !
    ! USES
    use FatesConstantsMod, only : pi => pi_const
    !
    ! ARGUMENTS
    real(r8), intent(in)  :: x  ! Independent variable in the fit
    integer, intent(in) :: num_edge_forest_bins
    real(r8), dimension(:), intent(in) :: efb_amplitudes
    real(r8), dimension(:), intent(in) :: efb_sigmas
    real(r8), dimension(:), intent(in) :: efb_centers
    real(r8), intent(in) :: efb_decay
    real(r8), dimension(:), pointer, intent(in) :: fraction_forest_in_bin
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
    logical :: do_norm

    if (present(norm)) then
       do_norm = norm
    else
       do_norm = .true.
    end if

    binloop: do b = 1, num_edge_forest_bins
       A = efb_amplitudes(b)
       if (b == num_edge_forest_bins) then
          ! Exponential
          fraction_forest_in_bin(b) = A * exp(-x/efb_decay)
       else
          ! Lognormal
         if (x == 0._r8) then
            ! Avoid divide-by-zero
            fraction_forest_in_bin(b) = 0._r8
            cycle
         end if
         mu = efb_centers(b)
         sigma = efb_sigmas(b)
         fraction_forest_in_bin(b) = gffeb_lognorm_numerator(x, A, mu, sigma) / gffeb_lognorm_denominator(x, sigma)
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


  end subroutine get_fraction_of_forest_in_each_bin


  subroutine assign_patches_to_bins(site, indices, index_forestpatches_to_allpatches, fraction_forest_in_each_bin, n_forest_patches, n_patches, area_forest_patches)
    ! DESCRIPTION
    ! Loops through forest patches from nearest to farthest from edge, assigning their
    ! area to edge bin(s).
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
    real(r8) :: remaining_to_assign_from_patch_m2
    real(r8) :: remaining_to_assign_to_bin_m2
    ! For checks
    real(r8), dimension(num_edge_forest_bins) :: bin_area_sums
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
       currentPatch%area_in_edge_forest_bins(:) = 0._r8
       remaining_to_assign_from_patch_m2 = currentPatch%area
       binloop: do b = 1, num_edge_forest_bins

          ! How much area is left for this bin?
          remaining_to_assign_to_bin_m2 = sum(fraction_forest_in_each_bin(1:b))*area_forest_patches - sum_forest_bins_so_far_m2
          if (remaining_to_assign_to_bin_m2 <= 0) then
             cycle
          end if

          ! Assign area
          currentPatch%area_in_edge_forest_bins(b) = min(remaining_to_assign_from_patch_m2, remaining_to_assign_to_bin_m2)
          remaining_to_assign_from_patch_m2 = remaining_to_assign_from_patch_m2 - currentPatch%area_in_edge_forest_bins(b)

          ! Update accounting
          sum_forest_bins_so_far_m2 = sum_forest_bins_so_far_m2 + currentPatch%area_in_edge_forest_bins(b)

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
       err_chk = currentPatch%area - sum(currentPatch%area_in_edge_forest_bins)
       if (abs(err_chk) > tol) then
          write(fates_log(),*) "ERROR: not enough or too much patch area was assigned to bins (check 2); remainder = ",err_chk
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

    end do forestpatchloop

    ! More checks
    bin_area_sums(:) = 0._r8
    currentPatch => site%oldest_patch
    allpatchloop_check: do while (associated(currentPatch))
       if (currentPatch%is_forest) then

          ! Check that all area of each forest patch is assigned
          err_chk = sum(currentPatch%area_in_edge_forest_bins) - currentPatch%area
          if (abs(err_chk) > tol) then
             write(fates_log(),*) "ERROR: unexpected patch forest bin sum (check 3); actual minus expected = ",err_chk
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if

          ! Accumulate site-wide area in each bin
          binloop_check4a: do b = 1, num_edge_forest_bins
             bin_area_sums(b) = bin_area_sums(b) + currentPatch%area_in_edge_forest_bins(b) / area_forest_patches
          end do binloop_check4a

       end if
       currentPatch => currentPatch%younger
    end do allpatchloop_check
    ! Check that fraction in each bin is what was expected
    binloop_check4b: do b = 1, num_edge_forest_bins
       err_chk = bin_area_sums(b) - fraction_forest_in_each_bin(b)
       if (abs(err_chk) > tol) then
          write(fates_log(),*) "ERROR: unexpected bin sum (check 4); actual minus expected = ",err_chk
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do binloop_check4b
  end subroutine assign_patches_to_bins


  subroutine calculate_edge_area(site)
    ! DESCRIPTION:
    ! Loop through forest patches in decreasing order of proximity, calculating the
    ! area of each patch that is in each edge bin.
    !
    ! USES:
    use EDParamsMod, only : ED_val_edgeforest_amplitudes, ED_val_edgeforest_sigmas, ED_val_edgeforest_centers, ED_val_edgeforest_decay
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
    real(r8) :: pct_nonforest
    real(r8), dimension(num_edge_forest_bins), target :: fraction_forest_in_each_bin

    ! Skip sites with no forest patches
    call get_number_of_forest_patches(site, n_forest_patches, area_forest_patches, area)
    if (n_forest_patches == 0) then
       return
    end if

    ! Allocate arrays
    allocate(indices(1:n_forest_patches))
    n_patches = get_number_of_patches(site)
    allocate(index_forestpatches_to_allpatches(1:n_patches))

    ! Get ranks
    call rank_forest_edge_proximity(site, indices, index_forestpatches_to_allpatches)

    ! Get percentage of nonforest area in each bin
    pct_nonforest = 100._r8 * (area - area_forest_patches) / area
    call get_fraction_of_forest_in_each_bin(pct_nonforest, num_edge_forest_bins, ED_val_edgeforest_amplitudes, ED_val_edgeforest_sigmas, ED_val_edgeforest_centers, ED_val_edgeforest_decay, fraction_forest_in_each_bin)

    ! Assign patches to bins
    call assign_patches_to_bins(site, indices, index_forestpatches_to_allpatches, fraction_forest_in_each_bin, n_forest_patches, n_patches, area_forest_patches)

    ! Clean up
    deallocate(indices)
    deallocate(index_forestpatches_to_allpatches)
  end subroutine calculate_edge_area




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