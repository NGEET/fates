module FatesEcotypesMod

  use FatesConstantsMod, only : r8 => fates_r8
  use EDTypesMod, only : ed_site_type
  use FatesPatchMod, only : fates_patch_type

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: is_patch_forest

contains

  ! =====================================================================================

  function does_patch_have_forest_tcthresh(patchptr, forest_tree_fraction_threshold)
  ! DESCRIPTION:
  ! Return boolean: Is this patch "forest"?
  !
  ! ARGUMENTS:
  type(fates_patch_type), intent(in), pointer :: patchptr  ! pointer to patch object
  real(r8), intent(in) :: forest_tree_fraction_threshold ! Tree fraction above which a patch is "forest"
  !
  ! RETURN VALUE
  logical :: does_patch_have_forest_tcthresh
  !
  ! LOCAL VARIABLES
  real(r8) :: tree_fraction = 0._r8

  if (patchptr%area > 0._r8) then
      tree_fraction = patchptr%total_tree_area / patchptr%area
  else
      tree_fraction = 0._r8
  end if

  does_patch_have_forest_tcthresh = tree_fraction > forest_tree_fraction_threshold

  end function does_patch_have_forest_tcthresh


  function does_patch_have_grass_bmthresh(patchptr, grass_biomass_threshold)
  ! DESCRIPTION:
  ! Return boolean: Does this patch have grass biomass above a threshold?
  !
  ! ARGUMENTS:
  type(fates_patch_type), intent(in), pointer :: patchptr  ! pointer to patch object
  real(r8), intent(in) :: grass_biomass_threshold ! Live grass biomass (kgC/m2) above which a patch is considered to "have grass"
  !
  ! RETURN VALUE
  logical :: does_patch_have_grass_bmthresh

  does_patch_have_grass_bmthresh = patchptr%livegrass > grass_biomass_threshold

  end function does_patch_have_grass_bmthresh


  function is_patch_forest(patchptr, forest_tree_fraction_threshold, grass_biomass_threshold)
  ! DESCRIPTION:
  ! Return boolean: Is this patch forest according to tree cover and, optionally, grass biomass?
  !
  ! ARGUMENTS:
  type(fates_patch_type), intent(in), pointer :: patchptr  ! pointer to patch object
  real(r8), intent(in) :: forest_tree_fraction_threshold ! Tree fraction above which a patch is "forest"
  real(r8), intent(in), optional :: grass_biomass_threshold ! Live grass biomass (kgC/m2) above which a patch is considered to "have grass"
  !
  ! RETURN VALUE
  logical :: is_patch_forest

  is_patch_forest = does_patch_have_forest_tcthresh(patchptr, forest_tree_fraction_threshold)
  if (is_patch_forest .and. present(grass_biomass_threshold)) then
     is_patch_forest = does_patch_have_grass_bmthresh(patchptr, grass_biomass_threshold)
  end if

  end function is_patch_forest


end module FatesEcotypesMod
