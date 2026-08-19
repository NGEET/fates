module EDTypesMockMod

  use FatesConstantsMod, only: r8 => fates_r8
  use EDTypesMod, only: ed_site_type
  use FatesPatchMod, only: fates_patch_type
  use FatesCohortMod, only: fates_cohort_type
  use FatesHydraulicsMemMod, only : ed_site_hydr_type, ed_cohort_hydr_type
  use PRTAllometricCarbonMod, only : InitPRTGlobalAllometricCarbon
  use EDCohortDynamicsMod, only : InitPRTObject
  use FatesInterfaceTypesMod, only : hlm_parteh_mode
  use PRTGenericMod, only : carbon_only

  implicit none

contains

  subroutine init_mock_site(site)
    type(ed_site_type), intent(inout), target :: site
    
    if (.not. associated(site%si_hydr)) then
       allocate(site%si_hydr)
    end if
    
  end subroutine init_mock_site

  subroutine init_mock_site_patch_cohort(site, patch, cohort)
    type(ed_site_type), pointer, intent(out) :: site
    type(fates_patch_type), pointer, intent(out) :: patch
    type(fates_cohort_type), pointer, intent(out) :: cohort
    
    integer :: i
    
    if (.not. associated(site)) then
      allocate(site)
      site%si_hydr => null()
    end if
    if (.not. associated(patch)) allocate(patch)
    if (.not. associated(cohort)) then
      allocate(cohort)
      cohort%co_hydr => null()
    end if

    site%oldest_patch => patch
    patch%tallest => cohort
    patch%patchno = 1
    patch%younger => null()
    patch%nocomp_pft_label = 1 ! Avoid being treated as bare ground
    cohort%shorter => null()

    call init_mock_site(site)
    
    ! Initialize the PRT global state (can be called safely multiple times in tests)
    hlm_parteh_mode = carbon_only
    call InitPRTGlobalAllometricCarbon()
    
    ! Initialize PRT object for this cohort if not already associated (idempotent for multiple tests)
    if (.not. associated(cohort%prt)) then
       call InitPRTObject(cohort%prt)
    end if

    
    ! Initialize all state variables to positive values to avoid crashes in hydraulics
    if (allocated(cohort%prt%variables)) then
       do i = 1, size(cohort%prt%variables)
          if (associated(cohort%prt%variables(i)%val)) then
             cohort%prt%variables(i)%val(:) = 10.0_r8
          end if
       end do
    end if
    
    if (.not. associated(cohort%co_hydr)) then
       allocate(cohort%co_hydr)
    end if

    ! Cohort density n = 10.0 indiv/m2 represents standard seedling/sapling canopy density
    cohort%n = 10.0_r8
    ! DBH = 10.0 cm represents 10cm diameter sapling trunk
    cohort%dbh = 10.0_r8
    ! Height = 10.0 m matches 10cm DBH allometry
    cohort%height = 10.0_r8
    ! Crown damage fraction = 0.0 represents undamaged intact canopy
    cohort%crowndamage = 0.0_r8
    ! Canopy trim fraction = 1.0 represents full untrimmed leaf canopy
    cohort%canopy_trim = 1.0_r8
    ! Stem efficiency fraction = 1.0 represents undamaged full sapwood conduit capacity
    cohort%efstem_coh = 1.0_r8
    ! Size class index = 1 selects smallest canopy size bin
    cohort%size_class = 1


    
  end subroutine init_mock_site_patch_cohort

end module EDTypesMockMod
