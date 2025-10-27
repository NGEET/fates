module FatesParametersInterface

  ! This module simply holds the instantiation
  ! of the generalized FATES parameter file that
  ! is provide by the JSON parser. ie pstruct
  ! This module is intentionaly low dependency
  ! Note also that JSONParameterUtilsMod
  ! only uses shr libraries

  use JSONParameterUtilsMod only: params_type

  implicit none
  private
  
  type(params_type) :: pstruct

  public :: pstruct

contains
  
end module FatesParametersInterface

