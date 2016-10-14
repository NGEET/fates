module FatesConstantsMod
  ! This module is used to define global _immutable_ data. Everything in
  ! this module must have the parameter attribute.

  implicit none

  public

  ! kinds
  integer, parameter :: fates_r8 = selected_real_kind(12) ! 8 byte real

  ! string lengths
  integer, parameter :: fates_avg_flag_length = 3
  integer, parameter :: fates_short_string_length = 32
  integer, parameter :: fates_long_string_length = 199

  integer, parameter :: fates_unset_int = -9999

  ! various magic numbers
  real(fates_r8), parameter ::  fates_special_value = 1.0e36_fates_r8 ! special value for real data, compatible with clm.
  integer, parameter :: fates_int_special_value = -9999    ! keep this negative to avoid conflicts with possible valid values

end module FatesConstantsMod
