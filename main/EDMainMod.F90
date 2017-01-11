module EDMainMod

  ! ===========================================================================
  ! Main ED module.    
  ! ============================================================================

  use shr_kind_mod         , only : r8 => shr_kind_r8
  
  use clm_varctl           , only : iulog
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type  
  use TemperatureType      , only : temperature_type
  use EDCohortDynamicsMod  , only : allocate_live_biomass, terminate_cohorts, fuse_cohorts, sort_cohorts, count_cohorts
  use EDPatchDynamicsMod   , only : disturbance_rates, fuse_patches, spawn_patches, terminate_patches
  use EDPhysiologyMod      , only : canopy_derivs, non_canopy_derivs, phenology, recruitment, trim_canopy
  use SFMainMod            , only : fire_model
  use EDtypesMod           , only : ncwd, numpft_ed, udata
  use EDtypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
  use FatesInterfaceMod    , only : bc_in_type
  use spmdMod              , only : masterproc

  implicit none
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: ed_ecosystem_dynamics
  public  :: ed_update_site
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  
  private :: ed_integrate_state_variables
  private :: ed_total_balance_check
  
  logical :: DEBUG  = .false.
  !
  ! 10/30/09: Created by Rosie Fisher
  !-----------------------------------------------------------------------

contains

  !-------------------------------------------------------------------------------!
  subroutine ed_ecosystem_dynamics(currentSite, bc_in, &
       atm2lnd_inst, &
       temperature_inst)
    !
    ! !DESCRIPTION:
    !  Core of ed model, calling all subsequent vegetation dynamics routines         
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target  :: currentSite
    type(bc_in_type)        , intent(in)             :: bc_in
    type(atm2lnd_type)      , intent(in)             :: atm2lnd_inst
    type(temperature_type)  , intent(in)             :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    !-----------------------------------------------------------------------

    if ( masterproc ) write(iulog,*) 'modelday',bc_in%model_day

    !**************************************************************************
    ! Fire, growth, biogeochemistry. 
    !**************************************************************************

    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file
   
    call ed_total_balance_check(currentSite, 0)
    
    call phenology(currentSite, bc_in )

    call fire_model(currentSite, atm2lnd_inst, temperature_inst)

    ! Calculate disturbance and mortality based on previous timestep vegetation.
    call disturbance_rates(currentSite)

    ! Integrate state variables from annual rates to daily timestep
    call ed_integrate_state_variables(currentSite, bc_in ) 

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organisation 
    !******************************************************************************

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))                 

       ! adds small cohort of each PFT
       call recruitment(0, currentSite, currentPatch)                

       currentPatch => currentPatch%younger
    enddo
       
    call ed_total_balance_check(currentSite,1)

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))

       ! puts cohorts in right order
       call sort_cohorts(currentPatch)            

       ! fuses similar cohorts
       call fuse_cohorts(currentPatch)            

       ! kills cohorts that are too small
       call terminate_cohorts(currentPatch)


       currentPatch => currentPatch%younger
    enddo
   
    call ed_total_balance_check(currentSite,2)

    !*********************************************************************************
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.
    !*********************************************************************************

    ! make new patches from disturbed land
    call spawn_patches(currentSite)       
   
    call ed_total_balance_check(currentSite,3)

    ! fuse on the spawned patches.
    call fuse_patches(currentSite)        
   
    call ed_total_balance_check(currentSite,4)

    ! kill patches that are too small
    call terminate_patches(currentSite)   
   
    call ed_total_balance_check(currentSite,5)

  end subroutine ed_ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine ed_integrate_state_variables(currentSite, bc_in )
    !
    ! !DESCRIPTION:
    ! FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout) :: currentSite
    type(bc_in_type)        , intent(in)   :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class 
    integer  :: ft                    ! Counter for PFT
    real(r8) :: small_no              ! to circumvent numerical errors that cause negative values of things that can't be negative
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    !-----------------------------------------------------------------------

    small_no = 0.0000000000_r8  ! Obviously, this is arbitrary.  RF - changed to zero

    do ft = 1,numpft_ed
       currentSite%dseed_dt(ft) = 0._r8  ! zero the dseed_dt at the site level before looping through patches and adding the fluxes from each patch
    end do
    currentSite%seed_rain_flux(:) = 0._r8  

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))

       currentPatch%age = currentPatch%age + udata%deltat
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(iulog,*) 'negative patch age?',currentPatch%age, &
               currentPatch%patchno,currentPatch%area
       endif

       ! Find the derivatives of the growth and litter processes. 
       call canopy_derivs(currentSite, currentPatch)
       
       ! Update Canopy Biomass Pools
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 

          cohort_biomass_store  = (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)
          currentCohort%dbh    = max(small_no,currentCohort%dbh    + currentCohort%ddbhdt    * udata%deltat )
          currentCohort%balive = currentCohort%balive + currentCohort%dbalivedt * udata%deltat 
          currentCohort%bdead  = max(small_no,currentCohort%bdead  + currentCohort%dbdeaddt  * udata%deltat )
          if ( DEBUG ) then
             write(iulog,*) 'EDMainMod dbstoredt I ',currentCohort%bstore, &
                  currentCohort%dbstoredt,udata%deltat
          end if
          currentCohort%bstore = currentCohort%bstore + currentCohort%dbstoredt * udata%deltat 
          if ( DEBUG ) then
             write(iulog,*) 'EDMainMod dbstoredt II ',currentCohort%bstore, &
                  currentCohort%dbstoredt,udata%deltat
          end if

          if( (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)*currentCohort%n<0._r8)then
            write(iulog,*) 'biomass is negative', currentCohort%n,currentCohort%balive, &
                 currentCohort%bdead,currentCohort%bstore
          endif

          if(abs((currentCohort%balive+currentCohort%bdead+currentCohort%bstore+udata%deltat*(currentCohort%md+ &
               currentCohort%seed_prod)-cohort_biomass_store)-currentCohort%npp_acc) > 1e-8_r8)then
             write(iulog,*) 'issue with c balance in integration', abs(currentCohort%balive+currentCohort%bdead+ &
                  currentCohort%bstore+udata%deltat* &
                 (currentCohort%md+currentCohort%seed_prod)-cohort_biomass_store-currentCohort%npp_acc)
          endif  

          ! THESE SHOULD BE MOVED TO A MORE "VISIBLE" LOCATION (RGK 10-2016)
          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_acc = 0.0_r8
          
          call allocate_live_biomass(currentCohort,1)
  
          currentCohort => currentCohort%taller

       enddo
      
       if ( DEBUG ) then
          write(6,*)'DEBUG18: calling non_canopy_derivs with pno= ',currentPatch%clm_pno
       endif

       call non_canopy_derivs( currentSite, currentPatch, bc_in)

       !update state variables simultaneously according to derivatives for this time period. 

       ! first update the litter variables that are tracked at the patch level
       do c = 1,ncwd
          currentPatch%cwd_ag(c) =  currentPatch%cwd_ag(c) + currentPatch%dcwd_ag_dt(c)* udata%deltat
          currentPatch%cwd_bg(c) =  currentPatch%cwd_bg(c) + currentPatch%dcwd_bg_dt(c)* udata%deltat
       enddo

       do ft = 1,numpft_ed
          currentPatch%leaf_litter(ft) = currentPatch%leaf_litter(ft) + currentPatch%dleaf_litter_dt(ft)* udata%deltat
          currentPatch%root_litter(ft) = currentPatch%root_litter(ft) + currentPatch%droot_litter_dt(ft)* udata%deltat
       enddo

       do c = 1,ncwd
          if(currentPatch%cwd_ag(c)<small_no)then
            write(iulog,*) 'negative CWD_AG', currentPatch%cwd_ag(c),CurrentSite%lat,currentSite%lon
            currentPatch%cwd_ag(c) = small_no
          endif
          if(currentPatch%cwd_bg(c)<small_no)then
            write(iulog,*) 'negative CWD_BG', currentPatch%cwd_bg(c),CurrentSite%lat,CurrentSite%lon
            currentPatch%cwd_bg(c) = small_no
          endif
       enddo

       do ft = 1,numpft_ed
          if(currentPatch%leaf_litter(ft)<small_no)then
            write(iulog,*) 'negative leaf litter numerical error', currentPatch%leaf_litter(ft),CurrentSite%lat,CurrentSite%lon,&
            currentPatch%dleaf_litter_dt(ft),currentPatch%leaf_litter_in(ft),currentPatch%leaf_litter_out(ft),currentpatch%age
            currentPatch%leaf_litter(ft) = small_no
          endif
          if(currentPatch%root_litter(ft)<small_no)then
               write(iulog,*) 'negative root litter numerical error', currentPatch%root_litter(ft), &
               currentPatch%droot_litter_dt(ft)* udata%deltat, &
               CurrentSite%lat,CurrentSite%lon
            currentPatch%root_litter(ft) = small_no
          endif
       enddo

     
       ! update cohort number. This needs to happen after the CWD_input and seed_input calculations as they 
       ! assume the pre-mortality currentCohort%n. 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 
         currentCohort%n = max(small_no,currentCohort%n + currentCohort%dndt * udata%deltat )  
         currentCohort => currentCohort%taller
       enddo

       currentPatch => currentPatch%older

    enddo

    ! at the site level, update the seed bank mass
    do ft = 1,numpft_ed
       currentSite%seed_bank(ft) = currentSite%seed_bank(ft) + currentSite%dseed_dt(ft)*udata%deltat
    enddo

    ! Check for negative values. Write out warning to show carbon balance. 
    do ft = 1,numpft_ed
       if(currentSite%seed_bank(ft)<small_no)then
          write(iulog,*) 'negative seedbank', currentSite%seed_bank(ft)
          currentSite%seed_bank(ft) = small_no
       endif
    enddo


  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite )
    !
    ! !DESCRIPTION:
    ! Calls routines to consolidate the ED growth process.
    ! Canopy Structure to assign canopy layers to cohorts
    ! Canopy Spread to figure out the size of tree crowns
    ! Trim_canopy to figure out the target leaf biomass. 
    ! Extra recruitment to fill empty patches.  
    !
    ! !USES:
    use EDCanopyStructureMod , only : canopy_spread, canopy_structure
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch   
    integer :: cohort_number ! To print out the number of cohorts.  
    integer :: g             ! Counter for sites
    !-----------------------------------------------------------------------

    call canopy_spread(currentSite)

    call ed_total_balance_check(currentSite,6)

    call canopy_structure(currentSite)

    call ed_total_balance_check(currentSite,7)

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       call terminate_cohorts(currentPatch) 

       ! FIX(SPM,040314) why is this needed for BFB restarts? Look into this at some point
       cohort_number = count_cohorts(currentPatch)  
       if ( DEBUG ) then
          write(iulog,*) 'tempCount ',cohort_number
       endif

       ! Note (RF)
       ! This breaks the balance check, but if we leave it out, then 
       ! the first new patch that isn't fused has no cohorts at the end of the spawn process
       ! and so there are radiation errors instead. 
       ! Fixing this would likely require a re-work of how seed germination works which would be tricky. 
       if(currentPatch%countcohorts < 1)then
          !write(iulog,*) 'ED: calling recruitment for no cohorts',currentPatch%siteptr%clmgcell,currentPatch%patchno
          !call recruitment(1, currentSite, currentPatch)
          ! write(iulog,*) 'patch empty',currentPatch%area,currentPatch%age
       endif

       currentPatch => currentPatch%younger    

    enddo

    ! FIX(RF,032414). This needs to be monthly, not annual
    if((udata%time_period == udata%n_sub-1))then 
       write(iulog,*) 'calling trim canopy' 
       call trim_canopy(currentSite)  
    endif

  end subroutine ed_update_site

  !-------------------------------------------------------------------------------!
  subroutine ed_total_balance_check (currentSite, call_index )
    !
    ! !DESCRIPTION:
    ! This routine looks at the carbon in and out of the ED model and compares it to 
    ! the change in total carbon stocks. 
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.  
    ! ed_allsites_inst%flux_out and ed_allsites_inst%flux_in are set where they occur 
    ! in the code. 
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: currentSite
    integer            , intent(in)    :: call_index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: biomass_stock   ! total biomass   in KgC/site
    real(r8) :: litter_stock    ! total litter    in KgC/site
    real(r8) :: seed_stock      ! total seed mass in KgC/site
    real(r8) :: total_stock     ! total ED carbon in KgC/site
    real(r8) :: change_in_stock ! Change since last time we set ed_allsites_inst%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!) 
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site

    ! nb. There is no time associated with these variables 
    ! because this routine can be called between any two 
    ! arbitrary points in code, even if no time has passed. 
    ! Also, the carbon pools are per site/gridcell, so that 
    ! we can account for the changing areas of patches. 

    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
    !-----------------------------------------------------------------------

    change_in_stock = 0.0_r8
    biomass_stock   = 0.0_r8
    litter_stock    = 0.0_r8

    seed_stock   =  sum(currentSite%seed_bank)

    currentPatch => currentSite%oldest_patch 
    do while(associated(currentPatch))

       litter_stock = litter_stock + currentPatch%area * (sum(currentPatch%cwd_ag)+ &
             sum(currentPatch%cwd_bg)+sum(currentPatch%leaf_litter)+sum(currentPatch%root_litter))
       currentCohort => currentPatch%tallest;
       
       do while(associated(currentCohort))
          
          biomass_stock =  biomass_stock + (currentCohort%bdead + currentCohort%balive + &
                currentCohort%bstore) * currentCohort%n
          currentCohort => currentCohort%shorter;
          
       enddo !end cohort loop 

       currentPatch => currentPatch%younger

    enddo !end patch loop

    total_stock     = biomass_stock + seed_stock +litter_stock
    change_in_stock = total_stock - currentSite%old_stock  
    net_flux        = currentSite%flux_in - currentSite%flux_out
    error           = abs(net_flux - change_in_stock)   

    if ( abs(error) > 10e-6 ) then
       write(iulog,*) 'total error: call index: ',call_index, &
                      'in:  ',currentSite%flux_in,   &
                      'out: ',currentSite%flux_out,  &
                      'net: ',net_flux,              &
                      'dstock: ',change_in_stock,    &
                      'error=net_flux-dstock:', error
       write(iulog,*) 'biomass,litter,seeds', biomass_stock,litter_stock,seed_stock
       write(iulog,*) 'lat lon',currentSite%lat,currentSite%lon
    endif

    currentSite%flux_in   = 0.0_r8
    currentSite%flux_out  = 0.0_r8  
    currentSite%old_stock = total_stock

 end subroutine ed_total_balance_check

end module EDMainMod
