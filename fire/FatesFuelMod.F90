module FatesFuelMod

  use FatesFuelClassesMod, only : nfsc
  
  implicit none
  private

  integer, parameter :: r8 = selected_real_kind(12)

  type, public :: fuel_type
    real(r8) :: loading(nfsc)      ! fuel loading of each fuel class [kg/m2]
    real(r8) :: total_loading      ! total fuel loading - DOES NOT INCLUDE TRUNKS [kg/m2]
    real(r8) :: frac_loading(nfsc) ! fractional loading of each fuel class [0-1] - TRUNKS SET TO 0.0
    real(r8) :: moisture(nfsc)     ! fuel moisture of each fuel class [m3/m3]
    real(r8) :: average_moisture   ! weighted average of fuel moisture across all fuel classes - DOES NOT INCLUDE TRUNKS [m3/m3]
    real(r8) :: bulk_density       ! weighted average of bulk density across all fuel classes - DOES NOT INCLUDE TRUNKS [kg/m3]
    real(r8) :: SAV                ! weighted average of surface area to volume ratio across all fuel classes - DOES NOT INCLUDE TRUNKS [/cm]
    real(r8) :: MEF                ! weighted average of moisture of extinction across all fuel classes - DOES NOT INCLUDE TRUNKS [m3/m3]

    contains
      procedure :: Init 
  end type fuel_type

  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type) :: this ! fuel class 

      ! just zero everything
      this%loading(:)       = 0.0_r8
      this%total_loading    = 0.0_r8
      this%frac_loading(:)  = 0.0_r8
      this%moisture(:)      = 0.0_r8
      this%average_moisture = 0.0_r8 
      this%bulk_density     = 0.0_r8
      this%SAV              = 0.0_r8
      this%MEF              = 0.0_r8 

    end subroutine Init 

end module FatesFuelMod