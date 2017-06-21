module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : ifalse
  use FatesGlobals              , only : endrun => fates_endrun
  use EDTypesMod                , only : nclmax
  use FatesGlobals              , only : fates_log
  use FatesInterfaceMod         , only : hlm_is_restart
  use EDPftvarcon               , only : EDPftvarcon_inst
  use EDEcophysConType          , only : EDecophyscon
  use EDGrowthFunctionsMod      , only : bdead, bleaf, dbh
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDPatchDynamicsMod        , only : create_patch
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type, area
  use EDTypesMod                , only : ncwd
  use EDTypesMod                , only : nuMWaterMem
  use EDTypesMod                , only : numpft_ed
  use FatesInterfaceMod         , only : bc_in_type
  use EDTypesMod                , only : use_fates_plant_hydro

  ! CIME GLOBALS
  use shr_log_mod               , only : errMsg => shr_log_errMsg

  implicit none
  private

  logical   ::  DEBUG = .false.

  logical, parameter :: do_inv_init = .true.

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

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

    ! diagnostic site-level cwd and litter fluxes
    site_in%CWD_AG_diagnostic_input_carbonflux(:) = 0._r8
    site_in%CWD_BG_diagnostic_input_carbonflux(:) = 0._r8
    site_in%leaf_litter_diagnostic_input_carbonflux(:) = 0._r8
    site_in%root_litter_diagnostic_input_carbonflux(:) = 0._r8

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

    if ( hlm_is_restart == ifalse ) then
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

       if ( hlm_is_restart == ifalse ) then
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
  subroutine init_patches( nsites, sites, bc_in)
    !
    ! !DESCRIPTION:
    ! initialize patches
    ! This may be call a near bare ground initialization, or it may
    ! load patches from an inventory.
     
    !
    ! !USES:
    use shr_file_mod, only      : shr_file_getUnit
    use shr_file_mod, only      : shr_file_freeUnit

    use EDPatchDynamicsMod     , only : dealloc_patch
    use EDParamsMod            , only : ED_val_maxspread
    use FatesPlantHydraulicsMod, only : updateSizeDepRhizHydProps 
    use FatesInventoryInitMod,   only : inv_file_list
    use FatesInventoryInitMod,   only : count_inventory_sites
    use FatesInventoryInitMod,   only : assess_inventory_sites
    use FatesInventoryInitMod,   only : set_inventory_edpatch_type1

    !
    ! !ARGUMENTS    
    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    type(bc_in_type), intent(in)               :: bc_in(nsites)
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
    type(ed_patch_type), pointer :: newpatch
    type(ed_patch_type), pointer :: oldpatch
   

    ! Census Initialization variables
    integer  :: file_unit
    integer  :: nfilesites ! number of sites in the inventory file list
    logical  :: lod       ! logical, file "O"pene"D"
    logical  :: lex       ! logical, file "EX"ists
    integer  :: ios       ! integer, "IO" status
    character(len=512) :: iostr
    logical, parameter :: do_inv_init = .true.
    integer, allocatable            :: inv_format_list(:)
    character(len=256), allocatable :: inv_css_list(:)
    character(len=256), allocatable :: inv_pss_list(:)
    real(r8), allocatable           :: inv_lat_list(:)
    real(r8), allocatable           :: inv_lon_list(:)
    integer                         :: invsite
    integer                         :: ipa   ! Patch index


    ! List out some nominal patch values that are used for Near Bear Ground initializations
    ! as well as initializing inventory
    ! ---------------------------------------------------------------------------------------------
    cwd_ag_local(:)      = 0.0_r8 !ED_val_init_litter -- arbitrary value for litter pools. kgC m-2
    cwd_bg_local(:)      = 0.0_r8 !ED_val_init_litter
    leaf_litter_local(:) = 0.0_r8
    root_litter_local(:) = 0.0_r8
    spread_local(:)      = ED_val_maxspread
    age                  = 0.0_r8
    ! ---------------------------------------------------------------------------------------------
    
    ! ---------------------------------------------------------------------------------------------
    ! Two primary options, either a Near Bear Ground (NBG) or Inventory based cold-start
    ! ---------------------------------------------------------------------------------------------

    if (do_inv_init) then

       ! I. Load the inventory list file, do some file handle checks
       ! ------------------------------------------------------------------------------------------

       file_unit = shr_file_getUnit()
       inquire(file=trim(inv_file_list),exist=lex,opened=lod)
       if( .not.lex ) then   ! The inventory file list DOE
          write(fates_log(), *) 'An inventory Initialization was requested.'
          write(fates_log(), *) 'However the inventory file: ',trim(inv_file_list),' DNE'
          write(fates_log(), *) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if( lod ) then        ! The inventory file should not be open
          write(fates_log(), *) 'The inventory list file is open but should not be.'
          write(fates_log(), *) 'Aborting.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       
       open(unit=file_unit,file=trim(inv_file_list),status='OLD',action='READ',form='FORMATTED')
       rewind(file_unit)

       ! There should be at least 1 line
       read(file_unit,fmt='(A)',iostat=ios) iostr
       read(file_unit,fmt='(A)',iostat=ios) iostr
       if( ios /= 0 ) then
          write(fates_log(), *) 'The inventory file does not contain at least two lines'
          write(fates_log(), *) 'of data, ie a header and 1 site.  Aborting.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       rewind(unit=file_unit)


       ! Count the number of sites that are listed in this file, and allocate storage arrays
       ! ------------------------------------------------------------------------------------------

       nfilesites = count_inventory_sites(file_unit)
       
       allocate(inv_format_list(nfilesites))
       allocate(inv_pss_list(nfilesites))
       allocate(inv_css_list(nfilesites))
       allocate(inv_lat_list(nfilesites))
       allocate(inv_lon_list(nfilesites))
       

       ! Check through the sites that are listed and do some sanity checks
       ! ------------------------------------------------------------------------------------------
       call assess_inventory_sites(file_unit, nfilesites, inv_format_list, &
                                   inv_pss_list, inv_css_list, &
                                   inv_lat_list, inv_lon_list)

       ! We can close the list file now.
       close(file_unit, iostat = ios)
       if( ios /= 0 ) then
          write(fates_log(), *) 'The inventory file needed to be closed, but was still open'
          write(fates_log(), *) 'aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       call shr_file_freeUnit(file_unit)


       ! For each site, identify the most proximal PSS/CSS couplet, read-in the data
       ! allocate linked lists and assign to memory
       do s = 1, nsites
          invsite = &
               minloc( (sites(s)%lat-inv_lat_list(:))**2.0_r8 + (sites(s)%lon-inv_lon_list(:))**2.0_r8 , dim=1)
          
          ! Open the PSS/CSS couplet and initialize the ED data structures.
          ! Lets start withe the PSS
          ! ---------------------------------------------------------------------------------------
          
          file_unit = shr_file_getUnit()
          open(unit=file_unit,file=trim(inv_pss_list(invsite)),status='OLD',action='READ',form='FORMATTED')
          rewind(file_unit)
          read(file_unit,fmt=*) iostr
          print*,"PATCH HEADER:"
          print*,trim(iostr)
          
          ipa = 0
          invpatchloop: do

             allocate(newpatch)

             newpatch%patchno = ipa
             newpatch%younger => null()
             newpatch%older   => null()

             ! This call doesn't do much asside from initializing the patch with
             ! nominal values, NaNs, zero's and allocating some vectors
             call create_patch(sites(s), newpatch, 0.0_r8, 0.0_r8, &
                  spread_local, cwd_ag_local, cwd_bg_local, leaf_litter_local,  &
                  root_litter_local) 

             
             if( inv_format_list(invsite) == 1 ) then
                
                call set_inventory_edpatch_type1(newpatch,file_unit,ipa,ios)
             
             end if
             
             ! If a new line was found in the inventory patch file,
             ! then it will return an IO status flag (ios) of 0
             ! In that case, the patch structure (newpatch) has been filled
             ! with relevant information.
             !
             ! Add it to the site's patch list
             ! ------------------------------------------------------------------------------------
             if(ios==0) then
             
                if(ipa == 0) then
                   ! This is the first patch to be added
                   ! It starts off as the oldest and youngest patch in the list
                   sites(s)%youngest_patch => newpatch
                   sites(s)%oldest_patch   => newpatch
                   oldpatch                => newpatch
                else
                   ! At least for now, we will assume that each subsequent
                   ! patch is a younger one.  We can sort when we are done
                   ! but lets not worry about it immediately
                   newpatch%older                 => oldpatch
                   newpatch%younger               => NULL()
                   sites(s)%youngest_patch        => newpatch
                   oldpatch                       => newpatch
                end if

             ! If a new line was NOT found in the inventory patch file,
             ! then no patch was populated and we should just deallocate the temporary
             ! and move along (tidy up site list and go to the next site)
             else
                
                call dealloc_patch(newpatch)
                deallocate(newpatch)
                exit              ! This should break the do loop

             end if

             
             
             ipa = ipa + 1
          end do invpatchloop

          stop

          
          ! Sort the patch list by age
          ! ---------------------------------------------------------------------------------------


          
          close(file_unit,iostat=ios)
          if( ios /= 0 ) then
             write(fates_log(), *) 'The pss file: ',inv_pss_list(invsite),' could not be closed'
             write(fates_log(), *) 'aborting'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          call shr_file_freeUnit(file_unit)
          stop
       end do

       deallocate(inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list)

    else

       

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

          call init_cohorts(newp, bc_in(s))

          ! This sets the rhizosphere shells based on the plant initialization
          ! The initialization of the plant-relevant hydraulics variables
          ! were set from a call inside of the init_cohorts()->create_cohort() subroutine
          if (use_fates_plant_hydro) then
             call updateSizeDepRhizHydProps(sites(s), bc_in(s))
          end if

       enddo

    end if

  end subroutine init_patches

  ! ============================================================================
  subroutine init_cohorts( patch_in, bc_in)
    !
    ! !DESCRIPTION:
    ! initialize new cohorts on bare ground
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), pointer  :: patch_in
    type(bc_in_type), intent(in)                 :: bc_in
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
       !temp_cohort%n           = 0.5_r8 * 0.0028_r8 * patch_in%area  ! BOC for fixed size runs EDecophyscon%initd(pft) * patch_in%area
       !temp_cohort%hite        = 28.65_r8                            ! BOC translates to DBH of 50cm. EDecophyscon%hgt_min(pft)
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
            temp_cohort%laimemory,  cstatus, temp_cohort%canopy_trim, 1, bc_in)

       deallocate(temp_cohort) ! get rid of temporary cohort

       endif

    enddo !numpft

    call fuse_cohorts(patch_in,bc_in)
    call sort_cohorts(patch_in)

  end subroutine init_cohorts

end module EDInitMod
