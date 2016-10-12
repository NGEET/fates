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

end module FatesConstantsMod
