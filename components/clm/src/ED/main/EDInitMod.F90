module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesGlobals              , only : endrun => fates_endrun
  use EDTypesMod                , only : nclmax
  use FatesGlobals              , only : fates_log
  use clm_varctl                , only : use_ed_spit_fire 
  use clm_time_manager          , only : is_restart
  use EDPftvarcon                 , only : EDPftvarcon_inst
  use EDEcophysConType          , only : EDecophyscon
  use EDGrowthFunctionsMod      , only : bdead, bleaf, dbh
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDPatchDynamicsMod        , only : create_patch
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type, area
  use EDTypesMod                , only : ncwd
  use EDTypesMod                , only : nuMWaterMem
  use EDTypesMod                , only : numpft_ed
  
  implicit none
  private

  logical   ::  DEBUG = .false.

  public  :: zero_site
  public  :: init_patches
  public  :: set_site_properties
  private :: init_cohorts

  ! ============================================================================

contains

  ! ============================================================================

  subroutine zero_site( site_in )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) ::  site_in
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    site_in%oldest_patch     => null() ! pointer to oldest patch at the site
    site_in%youngest_patch   => null() ! pointer to yngest patch at the site
    
    ! DISTURBANCE
    site_in%disturbance_rate = 0._r8  ! site level disturbance rates from mortality and fire.
    site_in%dist_type        = 0      ! disturbance dist_type id.
    site_in%total_burn_flux_to_atm = 0._r8 !

    ! PHENOLOGY 
    site_in%status           = 0    ! are leaves in this pixel on or off?
    site_in%dstatus          = 0
    site_in%ED_GDD_site      = nan  ! growing degree days
    site_in%ncd              = nan  ! no chilling days
    site_in%last_n_days(:)   = 999  ! record of last 10 days temperature for senescence model.
    site_in%leafondate       = 999  ! doy of leaf on
    site_in%leafoffdate      = 999  ! doy of leaf off
    site_in%dleafondate      = 999  ! doy of leaf on drought
    site_in%dleafoffdate     = 999  ! doy of leaf on drought
    site_in%water_memory(:)  = nan


    ! SEED
    site_in%seed_bank(:)     = 0._r8

    ! FIRE 
    site_in%acc_ni           = 0.0_r8     ! daily nesterov index accumulating over time. time unlimited theoretically.
    site_in%frac_burnt       = 0.0_r8     ! burn area read in from external file

    ! BGC Balance Checks
    site_in%fates_to_bgc_this_ts = 0.0_r8
    site_in%fates_to_bgc_last_ts = 0.0_r8

    ! termination and recruitment info
    site_in%terminated_nindivs(:,:,:) = 0._r8
    site_in%termination_carbonflux(:) = 0._r8
    site_in%recruitment_rate(:) = 0._r8

    ! demotion/promotion info
    site_in%demotion_rate(:) = 0._r8
    site_in%demotion_carbonflux = 0._r8
    site_in%promotion_rate(:) = 0._r8
    site_in%promotion_carbonflux = 0._r8

  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties( nsites, sites)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    

    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    !
    ! !LOCAL VARIABLES:
    integer  :: s
    real(r8) :: leafon
    real(r8) :: leafoff
    real(r8) :: stat
    real(r8) :: NCD
    real(r8) :: GDD
    real(r8) :: dstat
    real(r8) :: acc_NI
    real(r8) :: watermem
    integer  :: dleafoff
    integer  :: dleafon
    !----------------------------------------------------------------------

    if ( .not. is_restart() ) then
       !initial guess numbers for site condition.
       NCD      = 0.0_r8
       GDD      = 30.0_r8
       leafon   = 100.0_r8
       leafoff  = 300.0_r8
       stat     = 2
       acc_NI   = 0.0_r8
       dstat    = 2
       dleafoff = 300
       dleafon  = 100
       watermem = 0.5_r8

    else ! assignements for restarts

       NCD      = 1.0_r8 ! NCD should be 1 on restart
       GDD      = 0.0_r8
       leafon   = 0.0_r8
       leafoff  = 0.0_r8
       stat     = 1
       acc_NI   = 0.0_r8
       dstat    = 2
       dleafoff = 300
       dleafon  = 100
       watermem = 0.5_r8

    endif

    do s = 1,nsites
       sites(s)%ncd          = NCD
       sites(s)%leafondate   = leafon
       sites(s)%leafoffdate  = leafoff
       sites(s)%dleafoffdate = dleafoff
       sites(s)%dleafondate  = dleafon
       sites(s)%ED_GDD_site  = GDD

       if ( .not. is_restart() ) then
          sites(s)%water_memory(1:numWaterMem) = watermem
       end if

       sites(s)%status = stat
       !start off with leaves off to initialise
       sites(s)%dstatus= dstat
       
       sites(s)%acc_NI     = acc_NI
       sites(s)%frac_burnt = 0.0_r8
       sites(s)%old_stock  = 0.0_r8


    end do

    return
  end subroutine set_site_properties

  ! ============================================================================
  subroutine init_patches( nsites, sites)
    !
    ! !DESCRIPTION:
    !initialize patches on new ground
    !
    ! !USES:
    use EDParamsMod ,  only : ED_val_maxspread
    !
    ! !ARGUMENTS    
    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    !
    ! !LOCAL VARIABLES:
    integer  :: s
    real(r8) :: cwd_ag_local(ncwd)
    real(r8) :: cwd_bg_local(ncwd)
    real(r8) :: spread_local(nclmax)
    real(r8) :: leaf_litter_local(numpft_ed)
    real(r8) :: root_litter_local(numpft_ed)
    real(r8) :: age !notional age of this patch
    type(ed_patch_type), pointer :: newp
    !----------------------------------------------------------------------

    cwd_ag_local(:)      = 0.0_r8 !ED_val_init_litter -- arbitrary value for litter pools. kgC m-2
    cwd_bg_local(:)      = 0.0_r8 !ED_val_init_litter
    leaf_litter_local(:) = 0.0_r8
    root_litter_local(:) = 0.0_r8
    spread_local(:)      = ED_val_maxspread
    age                  = 0.0_r8

    !FIX(SPM,032414) clean this up...inits out of this loop
    do s = 1, nsites

       allocate(newp)

       newp%patchno = 1
       newp%younger => null()
       newp%older   => null()

       sites(s)%youngest_patch => newp
       sites(s)%youngest_patch => newp
       sites(s)%oldest_patch   => newp

       ! make new patch...
       call create_patch(sites(s), newp, age, AREA, &
            spread_local, cwd_ag_local, cwd_bg_local, leaf_litter_local,  &
            root_litter_local) 
       
       call init_cohorts(newp)

    enddo

  end subroutine init_patches

  ! ============================================================================
  subroutine init_cohorts( patch_in )
    !
    ! !DESCRIPTION:
    ! initialize new cohorts on bare ground
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), pointer  :: patch_in
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type),pointer :: temp_cohort
    integer :: cstatus
    integer :: pft
    !----------------------------------------------------------------------

    patch_in%tallest  => null()
    patch_in%shortest => null()

    do pft =  1,numpft_ed !FIX(RF,032414) - turning off veg dynamics

       if(EDecophyscon%initd(pft)>1.0E-7) then

       allocate(temp_cohort) ! temporary cohort

       temp_cohort%pft         = pft
       temp_cohort%n           = EDecophyscon%initd(pft) * patch_in%area
       temp_cohort%hite        = EDecophyscon%hgt_min(pft)
       temp_cohort%dbh         = Dbh(temp_cohort) ! FIX(RF, 090314) - comment out addition of ' + 0.0001_r8*pft   '  - seperate out PFTs a little bit...
       temp_cohort%canopy_trim = 1.0_r8
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + EDPftvarcon_inst%froot_leaf(pft) &
            + EDecophyscon%sapwood_ratio(temp_cohort%pft)*temp_cohort%hite)
       temp_cohort%b           = temp_cohort%balive + temp_cohort%bdead

       if( EDPftvarcon_inst%evergreen(pft) == 1) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft)
          temp_cohort%laimemory = 0._r8
          cstatus = 2
       endif

       if( EDPftvarcon_inst%season_decid(pft) == 1 ) then !for dorment places
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft) !stored carbon in new seedlings.
          if(patch_in%siteptr%status == 2)then 
             temp_cohort%laimemory = 0.0_r8
          else
             temp_cohort%laimemory = Bleaf(temp_cohort)
          endif
          ! reduce biomass according to size of store, this will be recovered when elaves com on.
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%status
       endif

       if ( EDPftvarcon_inst%stress_decid(pft) == 1 ) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft)
          temp_cohort%laimemory = Bleaf(temp_cohort)
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%dstatus
       endif

       if ( DEBUG ) write(fates_log(),*) 'EDInitMod.F90 call create_cohort '

       call create_cohort(patch_in, pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
            temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore, &
            temp_cohort%laimemory,  cstatus, temp_cohort%canopy_trim, 1)

       deallocate(temp_cohort) ! get rid of temporary cohort

       endif

    enddo !numpft

    call fuse_cohorts(patch_in)
    call sort_cohorts(patch_in)

  end subroutine init_cohorts

end module EDInitMod
