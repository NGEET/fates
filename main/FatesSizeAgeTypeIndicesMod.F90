module FatesSizeAgeTypeIndicesMod

  use FatesConstantsMod,     only : r8 => fates_r8

  use FatesInterfaceTypesMod,     only : nlevsclass
  use FatesInterfaceTypesMod,     only : nlevage
  use FatesInterfaceTypesMod,     only : nlevheight
  use FatesInterfaceTypesMod,     only : nlevcoage
  use EDTypesMod,                 only : nclmax
  use FatesInterfaceTypesMod,     only : nlevdamage
  use EDParamsMod,           only : ED_val_history_sizeclass_bin_edges
  use EDParamsMod,           only : ED_val_history_ageclass_bin_edges
  use EDParamsMod,           only : ED_val_history_height_bin_edges
  use EDParamsMod,           only : ED_val_history_coageclass_bin_edges
  use EDParamsMod,           only : ED_val_history_damage_bin_edges
  
  implicit none
  private ! Modules are private by default

  ! Make public necessary subroutines and functions
  public :: get_age_class_index
  public :: get_sizeage_class_index
  public :: sizetype_class_index
  public :: get_size_class_index
  public :: get_height_index
  public :: get_sizeagepft_class_index
  public :: get_agepft_class_index
  public :: get_cdamagesize_class_index
  public :: get_cdamagesizepft_class_index
  public :: coagetype_class_index
  public :: get_coage_class_index
  public :: get_agefuel_class_index
  public :: get_layersizetype_class_index
  
contains

  ! =====================================================================================
  
  function get_age_class_index(age) result( patch_age_class ) 

     real(r8), intent(in) :: age
     
     integer :: patch_age_class

     patch_age_class = count(age-ED_val_history_ageclass_bin_edges.ge.0.0_r8)

  end function get_age_class_index

    ! =====================================================================================
  

  function get_layersizetype_class_index(layer,dbh,pft) result(iclscpf)

    ! Get the 1D index for a canopy layer x size x pft triplet
    
    ! Arguments
    integer,intent(in) :: layer
    real(r8),intent(in) :: dbh
    integer,intent(in) :: pft

    integer :: size_class
    integer :: iclscpf
     
    size_class        = get_size_class_index(dbh)
    
    iclscpf = (pft-1)*nlevsclass*nclmax + (size_class-1)*nclmax + layer

    ! FOR ANALYSIS CODE, REVERSE: (assuming indices starting at 1):

    ! pft = ceiling(real(index,r8)/real(nlevsclass*nclmax,r8))
    ! size_class = ceiling(real(index-(pft-1)*nlevsclass*nclmax,r8)/real(nclmax,r8))
    ! layer = index - ((pft-1)*nlevsclass*nclmax + (size_class-1)*nclmax)
    
  end function get_layersizetype_class_index

  ! =====================================================================================

  function get_sizeage_class_index(dbh,age) result(size_by_age_class)
     
     ! Arguments
     real(r8),intent(in) :: dbh
     real(r8),intent(in) :: age

     integer             :: size_class
     integer             :: age_class
     integer             :: size_by_age_class
     
     size_class        = get_size_class_index(dbh)

     age_class         = get_age_class_index(age)
     
     size_by_age_class = (age_class-1)*nlevsclass + size_class

  end function get_sizeage_class_index

  !======================================================================================

  
  function get_cdamagesize_class_index(dbh,cdamage) result(cdamage_by_size_class)
     
     ! Arguments
     real(r8),intent(in) :: dbh
     integer,intent(in) :: cdamage

     integer             :: size_class
     integer             :: cdamage_by_size_class
     
     size_class        = get_size_class_index(dbh)
   
     cdamage_by_size_class = (cdamage-1)*nlevsclass + size_class

  end function get_cdamagesize_class_index


  ! =====================================================================================

  subroutine sizetype_class_index(dbh,pft,size_class,size_by_pft_class)
    
    ! Arguments
    real(r8),intent(in) :: dbh
    integer,intent(in)  :: pft
    integer,intent(out) :: size_class
    integer,intent(out) :: size_by_pft_class
    
    size_class        = get_size_class_index(dbh)
    
    size_by_pft_class = (pft-1)*nlevsclass+size_class

    return
 end subroutine sizetype_class_index

  ! =====================================================================================

  function get_size_class_index(dbh) result(cohort_size_class)

     real(r8), intent(in) :: dbh
     
     integer :: cohort_size_class
     
     cohort_size_class = count(dbh-ED_val_history_sizeclass_bin_edges.ge.0.0_r8)
     
  end function get_size_class_index

  ! =====================================================================================

 subroutine coagetype_class_index(coage,pft,coage_class,coage_by_pft_class)

  ! Arguments
  real(r8),intent(in) :: coage
  integer,intent(in)  :: pft
  integer,intent(out)  :: coage_class
  integer,intent(out)  :: coage_by_pft_class

  coage_class           = get_coage_class_index(coage)

  coage_by_pft_class    = (pft-1)*nlevcoage+coage_class

  return
 end subroutine coagetype_class_index

 ! ========================================================================================

 function get_coage_class_index(coage) result(cohort_coage_class)

   real(r8), intent(in) :: coage

   integer :: cohort_coage_class

   cohort_coage_class = count(coage-ED_val_history_coageclass_bin_edges.ge.0.0_r8)

 end function get_coage_class_index



  ! =====================================================================================

  function get_height_index(height) result(cohort_height_index)

     real(r8), intent(in) :: height
     
     integer :: cohort_height_index
     
     cohort_height_index = count(height-ED_val_history_height_bin_edges.ge.0.0_r8)
     
  end function get_height_index

  ! =====================================================================================

  function get_sizeagepft_class_index(dbh,age,pft) result(size_by_age_by_pft_class)
     
     ! Arguments
     real(r8),intent(in) :: dbh
     real(r8),intent(in) :: age
     integer,intent(in)  :: pft

     integer             :: size_class
     integer             :: age_class
     integer             :: size_by_age_by_pft_class
     
     size_class        = get_size_class_index(dbh)

     age_class         = get_age_class_index(age)
     
     size_by_age_by_pft_class = (age_class-1)*nlevsclass + size_class + &
          (pft-1) * nlevsclass * nlevage

  end function get_sizeagepft_class_index

   ! =====================================================================================

  function get_cdamagesizepft_class_index(dbh,cdamage,pft) result(cdamage_by_size_by_pft_class)
     
     ! Arguments
     real(r8),intent(in) :: dbh
     integer,intent(in) :: cdamage
     integer,intent(in)  :: pft

     integer             :: size_class
     integer             :: cdamage_by_size_by_pft_class
     
     size_class        = get_size_class_index(dbh)
 
     cdamage_by_size_by_pft_class = (cdamage-1)*nlevsclass + size_class + &
          (pft-1) * nlevsclass * nlevdamage

  end function get_cdamagesizepft_class_index

  ! =====================================================================================

  function get_agepft_class_index(age,pft) result(age_by_pft_class)
     
     ! Arguments
     real(r8),intent(in) :: age
     integer,intent(in)  :: pft

     integer             :: age_class
     integer             :: age_by_pft_class
     
     age_class         = get_age_class_index(age)
     
     age_by_pft_class = age_class + (pft-1) * nlevage

  end function get_agepft_class_index

  ! =====================================================================================

  function get_agefuel_class_index(age,fuel) result(age_by_fuel_class)
     
   ! Arguments
   real(r8),intent(in) :: age
   integer,intent(in)  :: fuel

   integer             :: age_class
   integer             :: age_by_fuel_class
   
   age_class         = get_age_class_index(age)
   
   age_by_fuel_class = age_class + (fuel-1) * nlevage

end function get_agefuel_class_index

end module FatesSizeAgeTypeIndicesMod
