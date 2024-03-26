module FatesFuelMod

  use FatesFuelClassesMod, only : nfsc
  use FatesLitterMod,      only : litter_type
  use FatesConstantsMod,   only : nearzero
  
  implicit none
  private

  integer, parameter :: r8 = selected_real_kind(12)

  type, public :: fuel_type
    real(r8) :: loading(nfsc)      ! fuel loading of each fuel class [kgC/m2]
    real(r8) :: moisture(nfsc)     ! fuel moisture of each fuel class [m3/m3]
    real(r8) :: frac_loading(nfsc) ! fractional loading of non-trunk fuel classes [0-1] 
    real(r8) :: total_loading      ! total fuel loading - DOES NOT INCLUDE TRUNKS [kgC/m2]
    real(r8) :: average_moisture   ! weighted average of fuel moisture across non-trunk fuel classes [m3/m3]
    real(r8) :: bulk_density       ! weighted average of bulk density across non-trunk fuel classes [kg/m3]
    real(r8) :: SAV                ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
    real(r8) :: MEF                ! weighted average of moisture of extinction across non-trunk fuel classes [m3/m3]

    contains
      
      procedure :: Init
      procedure :: CalculateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading

  end type fuel_type

  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class 

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

    !-------------------------------------------------------------------------------------

    subroutine CalculateLoading(this, litter, live_grass)
      ! DESCRIPTION:
      !   Calculates loading for each fuel type

      ! ARGUMENTS:
      class(fuel_type),          intent(inout) :: this       ! fuel class
      type(litter_type), target, intent(in)    :: litter     ! litter class
      real(r8),                  intent(in)    :: live_grass ! amount of live grass [kgC/m2]


      this%loading(dl_sf) = sum(litter%leaf_fines(:))
      this%loading(tw_sf:tr_sf) = litter%ag_cwd(:)
      this%loading(lg_sf) = live_grass

    end subroutine CalculateLoading

    !-------------------------------------------------------------------------------------

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sums up the loading

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      this%total_loading = sum(this%loading(:))

    end subroutine SumLoading

    !-------------------------------------------------------------------------------------

    subroutine CalculateFractionalLoading(this)
      ! DESCRIPTION:
      !   Calculates fractional loading

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      ! sum up loading just in case
      call this%SumLoading()

      if (this%total_loading > nearzero) then
        this%frac_loading(:) = this%loading(:)/this%total_loading
      else 
        this%frac_loading(:) = 0.0_r8
      end if 

    end subroutine CalculateFractionalLoading

    !-------------------------------------------------------------------------------------

end module FatesFuelMod