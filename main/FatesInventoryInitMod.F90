module FatesInventoryInitMod

   !-----------------------------------------------------------------------------------------------
   ! This module contains the majority of the code used to read forest inventory data and
   ! initialize a simulation from that data.  Some important points:
   ! - This procedure is called through the host model's "cold-start" initialization and not a
   !   restart type of simulation.
   ! - This procedure, if called from a massive grid, is both probably inappropriate, and probably
   !   time consuming.
   ! - This procedure is assumed to be called over a small subset of sites, for instance a single
   !   site, or a small collection of sparse/irregularly spaced group of sites
   !
   ! Created: Ryan Knox June 2017
   ! This code borrows heavily in concept from what is done in ED2.  
   ! At the time of writing this ED2 is unlicensed, and only concepts were borrowed with no direct
   ! code copied.
   !
   !
   ! Update: Jessica Needham October 2023
   ! As discussed in FATES issue #1062 we decided to remove columns not used in FATES from the
   ! PSS and CSS files. 
  !-----------------------------------------------------------------------------------------------

   ! CIME GLOBALS

   use shr_log_mod      , only : errMsg => shr_log_errMsg

   ! FATES GLOBALS
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : pi_const
   use FatesConstantsMod, only : itrue
   use FatesConstantsMod, only : nearzero
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log
   use EDParamsMod      , only : regeneration_model
   use FatesInterfaceTypesMod, only : bc_in_type
   use FatesInterfaceTypesMod, only : hlm_inventory_ctrl_file
   use FatesInterfaceTypesMod, only : nleafage
   use FatesInterfaceTypesMod, only : hlm_current_tod
   use FatesInterfaceTypesMod, only : numpft
   use FatesLitterMod   , only : litter_type
   use EDTypesMod       , only : ed_site_type
   use FatesPatchMod    , only : fates_patch_type
   use FatesCohortMod   , only : fates_cohort_type
   use EDTypesMod       , only : area
   use FatesConstantsMod, only : leaves_on
   use FatesConstantsMod, only : leaves_off
   use FatesConstantsMod, only : ihard_stress_decid
   use FatesConstantsMod, only : isemi_stress_decid
   use PRTGenericMod    , only : num_elements
   use PRTGenericMod    , only : element_list
   use EDTypesMod       , only : phen_cstat_nevercold
   use EDTypesMod       , only : phen_cstat_iscold
   use EDTypesMod       , only : phen_dstat_timeoff
   use EDTypesMod       , only : phen_dstat_moistoff
   use PRTParametersMod , only : prt_params
   use EDPftvarcon      , only : EDPftvarcon_inst
   use FatesInterfaceTypesMod, only : hlm_parteh_mode
   use EDCohortDynamicsMod, only : InitPRTObject
   use PRTGenericMod,       only : prt_carbon_allom_hyp
   use PRTGenericMod,       only : prt_cnp_flex_allom_hyp
   use PRTGenericMod,       only : prt_vartypes
   use PRTGenericMod,       only : leaf_organ
   use PRTGenericMod,       only : fnrt_organ
   use PRTGenericMod,       only : sapw_organ
   use PRTGenericMod,       only : store_organ
   use PRTGenericMod,       only : struct_organ
   use PRTGenericMod,       only : repro_organ
   use PRTGenericMod,       only : carbon12_element
   use PRTGenericMod,       only : nitrogen_element
   use PRTGenericMod,       only : phosphorus_element
   use PRTGenericMod,       only : SetState
   use FatesConstantsMod,   only : primaryland
   use FatesRunningMeanMod, only : ema_lpa
   use PRTGenericMod,       only : StorageNutrientTarget
   use FatesConstantsMod,   only : fates_unset_int
   use EDCanopyStructureMod, only : canopy_summarization, canopy_structure
   use FatesRadiationMemMod, only : num_swb
   implicit none
   private

   ! This derived type is to allow an array of pointers to the LL patch structure
   ! This is different than allocating a vector of patches. This is needed for
   ! quickly matching cohort string identifiers, the indices that match thos identifiers
   ! with a patch.  BY having a vector of patch pointers that lines up with the string
   ! identifier array, this can be done quickly.
   type pp_array
      type(fates_patch_type), pointer :: cpatch
   end type pp_array

   character(len=*), parameter, private :: sourcefile =  __FILE__

   logical, parameter :: debug_inv = .false.       ! Debug flag for devs

   ! String length specifiers
   integer, parameter :: patchname_strlen = 64
   integer, parameter :: cohortname_strlen = 64
   integer, parameter :: line_strlen = 512
   integer, parameter :: path_strlen = 256


   real(r8), parameter :: max_site_adjacency_deg = 0.05_r8 ! The maximum distance in degrees
                                                           ! allowed between a site's coordinate
                                                           ! defined in model memory and a physical
                                                           ! site listed in the file

   logical, parameter :: do_inventory_out = .false.


   public :: initialize_sites_by_inventory

contains

   ! ==============================================================================================

   subroutine initialize_sites_by_inventory(nsites,sites,bc_in)

      ! !USES:
      use shr_file_mod, only        : shr_file_getUnit
      use shr_file_mod, only        : shr_file_freeUnit
      use FatesConstantsMod, only   : nearzero
      use EDPatchDynamicsMod, only  : fuse_patches
      use EDCohortDynamicsMod, only : fuse_cohorts
      use EDPatchDynamicsMod, only  : patch_pft_size_profile

      ! Arguments
      integer,            intent(in)               :: nsites
      type(ed_site_type), intent(inout), target    :: sites(nsites)
      type(bc_in_type),   intent(in)               :: bc_in(nsites)

      ! Locals
      type(ed_site_type), pointer                  :: currentSite
      type(fates_patch_type), pointer                 :: currentpatch
      type(fates_cohort_type), pointer                :: currentcohort
      type(fates_patch_type), pointer                 :: newpatch
      type(fates_patch_type), pointer                 :: olderpatch
      type(fates_patch_type), pointer                 :: head_of_unsorted_patch_list
      type(fates_patch_type), pointer                 :: next_in_unsorted_patch_list
      integer                                      :: sitelist_file_unit   ! fortran file unit for site list
      integer                                      :: pss_file_unit        ! fortran file unit for the pss file
      integer                                      :: css_file_unit        ! fortran file unit for the css file
      integer                                      :: nfilesites           ! number of sites in file list
      logical                                      :: lopen                ! logical, file is open
      logical                                      :: lexist               ! logical, file exists
      integer                                      :: ios                  ! integer, "IO" status
      character(len=line_strlen)                   :: header_str           ! large string for whole lines
      real(r8)                                     :: age_init             ! dummy value for creating a patch
      real(r8)                                     :: area_init            ! dummy value for creating a patch
      integer                                      :: s                    ! site index
      integer                                      :: ipa                  ! patch index
      integer                                      :: iv, ft, ic
      integer                                      :: total_cohorts        ! cohort counter for error checking
      integer,                         allocatable :: inv_format_list(:)   ! list of format specs
      character(len=path_strlen),      allocatable :: inv_css_list(:)      ! list of css file names
      character(len=path_strlen),      allocatable :: inv_pss_list(:)      ! list of pss file names

      real(r8),                        allocatable :: inv_lat_list(:)      ! list of lat coords
      real(r8),                        allocatable :: inv_lon_list(:)      ! list of lon coords
      integer                                      :: invsite              ! index of inventory site
                                                                           ! closest to actual site
      integer                                      :: el                   ! loop counter for number of elements
      character(len=patchname_strlen)              :: patch_name           ! patch ID string in the file
      integer                                      :: npatches             ! number of patches found in PSS
      type(pp_array),                  allocatable :: patch_pointer_vec(:) ! vector of pointers to patch LL
      character(len=patchname_strlen), allocatable :: patch_name_vec(:)    ! vector of patch ID strings
      real(r8)                                     :: basal_area_postf     ! basal area before fusion (m2/ha)
      real(r8)                                     :: basal_area_pref      ! basal area after fusion (m2/ha)
    
      ! I. Load the inventory list file, do some file handle checks
      ! ------------------------------------------------------------------------------------------

      sitelist_file_unit = shr_file_getUnit()


      inquire(file=trim(hlm_inventory_ctrl_file),exist=lexist,opened=lopen)
      if( .not.lexist ) then   ! The inventory file list DNE
         write(fates_log(), *) 'An inventory Initialization was requested.'
         write(fates_log(), *) 'However the inventory file: ',trim(hlm_inventory_ctrl_file),' DNE'
         write(fates_log(), *) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      if( lopen ) then        ! The inventory file should not be open
         write(fates_log(), *) 'The inventory list file is open but should not be.'
         write(fates_log(), *) 'Aborting.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      open(unit=sitelist_file_unit,file=trim(hlm_inventory_ctrl_file),status='OLD',action='READ',form='FORMATTED')
      rewind(sitelist_file_unit)

      ! There should be at least 1 line
      read(sitelist_file_unit,fmt='(A)',iostat=ios) header_str
      read(sitelist_file_unit,fmt='(A)',iostat=ios) header_str
      if( ios /= 0 ) then
         write(fates_log(), *) 'The inventory file does not contain at least two lines'
         write(fates_log(), *) 'of data, ie a header and 1 site.  Aborting.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      rewind(unit=sitelist_file_unit)


      ! Count the number of sites that are listed in this file, and allocate storage arrays
      ! ------------------------------------------------------------------------------------------

      nfilesites = count_inventory_sites(sitelist_file_unit)

      allocate(inv_format_list(nfilesites))
      allocate(inv_pss_list(nfilesites))
      allocate(inv_css_list(nfilesites))
      allocate(inv_lat_list(nfilesites))
      allocate(inv_lon_list(nfilesites))


      ! Check through the sites that are listed and do some sanity checks
      ! ------------------------------------------------------------------------------------------
      call assess_inventory_sites(sitelist_file_unit, nfilesites, inv_format_list, &
            inv_pss_list, inv_css_list, &
            inv_lat_list, inv_lon_list)

      ! We can close the list file now
      close(sitelist_file_unit, iostat = ios)
      if( ios /= 0 ) then
         write(fates_log(), *) 'The inventory file needed to be closed, but was still open'
         write(fates_log(), *) 'aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      call shr_file_freeUnit(sitelist_file_unit)


      ! For each site, identify the most proximal PSS/CSS couplet, read-in the data
      ! allocate linked lists and assign to memory
      do s = 1, nsites
         invsite = &
               minloc( (sites(s)%lat-inv_lat_list(:))**2.0_r8 + &
               (sites(s)%lon-inv_lon_list(:))**2.0_r8 , dim=1)

         ! Do a sanity check on the distance separation between physical site and model site
         if ( sqrt( (sites(s)%lat-inv_lat_list(invsite))**2.0_r8 + &
               (sites(s)%lon-inv_lon_list(invsite))**2.0_r8 ) > max_site_adjacency_deg ) then
            write(fates_log(), *) 'Model site at lat:',sites(s)%lat,' lon:',sites(s)%lon
            write(fates_log(), *) 'has no reasonably proximal site in the inventory site list.'
            write(fates_log(), *) 'Closest is at lat:',inv_lat_list(invsite),' lon:',inv_lon_list(invsite)
            write(fates_log(), *) 'Separation must be less than ',max_site_adjacency_deg,' degrees'
            write(fates_log(), *) 'Exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         ! Open the PSS/CSS couplet and initialize the ED data structures.
         ! Lets start withe the PSS
         ! ---------------------------------------------------------------------------------------

         pss_file_unit = shr_file_getUnit()
         open(unit=pss_file_unit,file=trim(inv_pss_list(invsite)), &
               status='OLD',action='READ',form='FORMATTED')
         rewind(pss_file_unit)
         read(pss_file_unit,fmt=*) header_str

         ! Do one quick pass through just to count lines
         ipa = 0
         countpatchloop: do
            read(pss_file_unit,fmt=*,iostat=ios) header_str
            if(ios/=0) exit
            ipa = ipa + 1
         end do countpatchloop
         rewind(pss_file_unit)
         read(pss_file_unit,fmt=*) header_str

         npatches = ipa
         allocate(patch_name_vec(npatches))
         allocate(patch_pointer_vec(npatches))


         do ipa=1,npatches

            ! This call doesn't do much asside from initializing the patch with
            ! nominal values, NaNs, zero's and allocating some vectors. We should
            ! be able to get the following values from the patch files. But on
            ! the patch creation step, we don't have that information.

            age_init            = 0.0_r8
            area_init           = 0.0_r8
            allocate(newpatch)
            call newpatch%Create(age_init, area_init, primaryland,           &
               fates_unset_int, num_swb, numpft, sites(s)%nlevsoil,         &
               hlm_current_tod, regeneration_model)

            newpatch%patchno = ipa
            newpatch%younger => null()
            newpatch%older   => null()

            if( inv_format_list(invsite) == 1 ) then
               call set_inventory_patch_type1(newpatch,pss_file_unit,ipa,ios,patch_name)
            end if

            ! Add it to the site's patch list
            ! ------------------------------------------------------------------------------------

            patch_name_vec(ipa)           = trim(patch_name)
            patch_pointer_vec(ipa)%cpatch => newpatch

            if(ipa == 1) then
               ! This is the first patch to be added
               ! It starts off as the oldest and youngest patch in the list
               sites(s)%youngest_patch => newpatch
               sites(s)%oldest_patch   => newpatch
            else

               ! Insert this patch into the patch LL
               ! First check the two end-cases

               ! Youngest Patch
               if(newpatch%age<=sites(s)%youngest_patch%age)then
                  newpatch%older                  => sites(s)%youngest_patch
                  newpatch%younger                => null()
                  sites(s)%youngest_patch%younger => newpatch
                  sites(s)%youngest_patch         => newpatch

                  ! Oldest Patch
               else if(newpatch%age>sites(s)%oldest_patch%age)then
                  newpatch%older              => null()
                  newpatch%younger            => sites(s)%oldest_patch
                  sites(s)%oldest_patch%older => newpatch
                  sites(s)%oldest_patch       => newpatch

                  ! Somewhere in the middle
               else
                  currentpatch => sites(s)%youngest_patch
                  do while(associated(currentpatch))
                     olderpatch => currentpatch%older
                     if(associated(currentpatch%older)) then
                        if(newpatch%age >= currentpatch%age .and. &
                              newpatch%age < olderpatch%age) then
                           ! Set the new patches pointers
                           newpatch%older   => currentpatch%older
                           newpatch%younger => currentpatch
                           ! Fix the patch's older pointer
                           currentpatch%older => newpatch
                           ! Fix the older patch's younger pointer
                           olderpatch%younger => newpatch
                        end if
                     end if
                     currentPatch => olderpatch
                  enddo
               end if
            end if
         end do

         close(pss_file_unit,iostat=ios)
         if( ios /= 0 ) then
            write(fates_log(), *) 'The pss file: ',inv_pss_list(invsite),' could not be closed'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         call shr_file_freeUnit(pss_file_unit)

         if(debug_inv) then
            write(fates_log(),*) 'Raw List of Inventory Patches, Age Sorted:'
            currentpatch => sites(s)%youngest_patch
            do while(associated(currentpatch))
               write(fates_log(),*) ' AGE: ',currentpatch%age,' AREA: ',currentpatch%area
               currentPatch => currentpatch%older
            enddo
         end if

         
         ! OPEN THE CSS FILE
         ! ---------------------------------------------------------------------------------------
         css_file_unit = shr_file_getUnit()
         open(unit=css_file_unit,file=trim(inv_css_list(invsite)), &
               status='OLD',action='READ',form='FORMATTED')
         rewind(css_file_unit)
         read(css_file_unit,fmt=*) header_str

         ! Read in each cohort line. Each line is associated with a patch from the PSS
         ! file via a patch name identification string.  We pass the whole site pointer
         ! to this routine, because inside the routine we identify the patch by making
         ! comparisons with patch_name_vec and identifying the patch pointer
         ! from patch_pointer_vec

         invcohortloop: do
            if ( inv_format_list(invsite) == 1 ) then
               call set_inventory_cohort_type1(sites(s),bc_in(s),css_file_unit, &
                     npatches, patch_pointer_vec,patch_name_vec, ios)
            end if
            if ( ios/=0 ) exit
         end do invcohortloop

         close(css_file_unit,iostat=ios)
         if( ios/=0 ) then
            write(fates_log(), *) 'The css file: ',inv_css_list(invsite),' could not be closed'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         call shr_file_freeUnit(css_file_unit)

         deallocate(patch_pointer_vec,patch_name_vec)

         ! Report Basal Area (as a check on if things were read in)
         ! ------------------------------------------------------------------------------
         basal_area_pref = 0.0_r8
         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentcohort => currentpatch%tallest
            do while(associated(currentcohort))
               basal_area_pref = basal_area_pref + &
                    currentcohort%n*0.25*((currentcohort%dbh/100.0_r8)**2.0_r8)*pi_const
               currentcohort => currentcohort%shorter
            end do
            currentPatch => currentpatch%older
         enddo

         write(fates_log(),*) '-------------------------------------------------------'
         write(fates_log(),*) 'Basal Area from inventory, BEFORE fusion'
         write(fates_log(),*) 'Lat: ',sites(s)%lat,' Lon: ',sites(s)%lon
         write(fates_log(),*) basal_area_pref,' [m2/ha]'
         write(fates_log(),*) '-------------------------------------------------------'
                  
         ! Update the patch index numbers and fuse the cohorts in the patches
         ! ----------------------------------------------------------------------------------------
         ipa=1
         total_cohorts = 0

         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentpatch%patchno = ipa
            ipa=ipa+1

            ! Perform Cohort Fusion
            call fuse_cohorts(sites(s), currentpatch,bc_in(s))
            call currentpatch%SortCohorts()

            ! This calculates %num_cohorts
            call currentPatch%CountCohorts()
            total_cohorts = total_cohorts + currentPatch%num_cohorts

            currentPatch => currentpatch%older
         enddo

         if(total_cohorts .eq. 0)then
            write(fates_log(), *) 'The inventory initialization produced no cohorts.'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! Fuse patches
         ! ----------------------------------------------------------------------------------------
         call fuse_patches(sites(s), bc_in(s) )

         ! Report Basal Area (as a check on if things were read in)
         ! ----------------------------------------------------------------------------------------
         !call canopy_structure(sites(s),bc_in(s))
         basal_area_postf = 0.0_r8
         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentcohort => currentpatch%tallest
            do while(associated(currentcohort))
               basal_area_postf = basal_area_postf + &
                     currentcohort%n*0.25*((currentcohort%dbh/100.0_r8)**2.0_r8)*pi_const
               currentcohort => currentcohort%shorter
            end do
            
            currentPatch => currentpatch%older
         enddo



         write(fates_log(),*) '-------------------------------------------------------'
         write(fates_log(),*) 'Basal Area from inventory, AFTER fusion'
         write(fates_log(),*) 'Lat: ',sites(s)%lat,' Lon: ',sites(s)%lon
         write(fates_log(),*) basal_area_postf,' [m2/ha]'
         write(fates_log(),*) '-------------------------------------------------------'

         ! If this is flagged as true, the post-fusion inventory will be written to file
         ! in the run directory.
         
         if(do_inventory_out)then
             call write_inventory_type1(sites(s))
         end if

      end do
      
      deallocate(inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list)

      return
   end subroutine initialize_sites_by_inventory

   ! ==============================================================================================

   function count_inventory_sites(sitelist_file_unit) result(nsites)

      ! Simple function that counts the number of lines in the inventory descriptor file

      ! Arguments
      integer, intent(in)        :: sitelist_file_unit

      ! Locals
      character(len=line_strlen) :: header_str
      character(len=line_strlen) :: site_str
      integer                    :: ios
      integer                    :: nsites


      ! Set the file position to the top of the file
      ! Read in the header line
      ! Read through sites and check coordinates and file existence
      rewind(unit=sitelist_file_unit)
      read(sitelist_file_unit,fmt='(A)') header_str
      nsites = 0
      do
         read(sitelist_file_unit,fmt='(A)',iostat=ios) site_str
         if (ios/=0) exit
         nsites = nsites + 1
      end do

      return
   end function count_inventory_sites

   ! ==============================================================================================

   subroutine assess_inventory_sites(sitelist_file_unit,nsites, inv_format_list, &
         inv_pss_list,inv_css_list, &
         inv_lat_list,inv_lon_list)

      ! -------------------------------------------------------------------------------------------
      ! This subroutine looks through the inventory descriptor file
      ! and line by line reads information about the available inventory
      ! sites, and saves their information (such as location and file path)
      ! to arrays.  This routine also does some simple checks to make
      ! sure it is not reading nonsense
      !
      ! File Format for the inventory site file:
      ! 1 line header
      ! 1 line listing each available inventory site with the following fields:
      ! type     latitude    longitude   pss-name   css-name
      !
      ! The fields for each site are described as follows:
      !
      ! <short-name>    <value-kind>     <description>
      !
      ! type            integer          We will accomodate different file format with different
      !                                  field values as the need arises. format 1 will read in
      !                                  datasets via "set_inventory_patch_type1()",
      !                                  "set_inventory_cohort_type1()"
      !
      ! latitude        float            The geographic latitude coordinate of the site
      ! longitude       float            The geogarphic longitude coordinate of the site
      ! pss-name        string           The full path to the patch descriptor file (PSS)
      ! css-name        string           The full path to the cohort descriptor file (CSS)
      ! -------------------------------------------------------------------------------------------


      ! Arguments
      integer, intent(in)                      :: sitelist_file_unit       ! file unit for sitelist
      integer, intent(in)                      :: nsites                   ! number of inventory sites
      integer, intent(inout)                   :: inv_format_list(nsites)  ! array of formats for each inventory site
      character(len=path_strlen),intent(inout) :: inv_pss_list(nsites)     ! array of pss file paths for each site
      character(len=path_strlen),intent(inout) :: inv_css_list(nsites)     ! array of css file paths for each site
      real(r8),intent(inout)                   :: inv_lat_list(nsites)     ! array of latitudes for each site
      real(r8),intent(inout)                   :: inv_lon_list(nsites)     ! array of longitudes for each site

      ! Locals
      character(len=line_strlen)               :: header_str  ! a string to hold the header information
      character(len=line_strlen)               :: site_str    ! a string to hold each site-line in the file
      integer                                  :: isite       ! site index
      integer                                  :: ios         ! fortran read status flag
      character(len=path_strlen)               :: pss_file
      character(len=path_strlen)               :: css_file
      real(r8)                                 :: site_lat    ! inventory site latitude
      real(r8)                                 :: site_lon    ! site longitude
      integer                                  :: iblnk       ! Index used for string parsing
      integer                                  :: file_format ! format type (1=legacy ED pss/css)
      logical                                  :: lexist      ! file existence flag

      rewind(unit=sitelist_file_unit)
      read(sitelist_file_unit,fmt='(4A)') header_str

      do isite=1,nsites

         ! Read in the whole line
         read(sitelist_file_unit,fmt='(a)',iostat=ios) site_str

         ! Parse the format identifier
         read(site_str,*) file_format

         ! Parse the latitude
         site_str = adjustl(site_str)
         iblnk = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,*) site_lat

         ! Parse the longitude
         site_str = adjustl(site_str)
         iblnk = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,*) site_lon

         ! Parse the pss file name
         site_str = adjustl(site_str)
         iblnk    = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         iblnk    = index(site_str,' ')
         read(site_str(:iblnk),fmt='(1A)') pss_file

         ! Parse the css file name
         site_str = adjustl(site_str)
         iblnk    = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,fmt='(1A)') css_file


         if ( site_lat < -90.0_r8 .or. site_lat > 90.0_r8 ) then
            write(fates_log(), *) 'read invalid latitude: ',site_lat,' from inventory site list'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! Longitude should be converted to positive coordinate

         if( site_lon<0.0_r8 ) site_lon = 360.0_r8 + site_lon

         if ( site_lon < 0.0_r8 .or. site_lon > 360.0_r8 ) then
            write(fates_log(), *) 'read invalid longitude: ',site_lon,' from inventory site list'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         inquire(file=trim(pss_file),exist=lexist)
         if( .not.lexist ) then
            write(fates_log(), *) 'the following pss file could not be found:'
            write(fates_log(), *) trim(pss_file)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         inquire(file=trim(css_file),exist=lexist)
         if( .not.lexist ) then
            write(fates_log(), *) 'the following css file could not be found:'
            write(fates_log(), *) trim(css_file)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! If we have made it to this point, then in all likelihood, the PSS/CSS
         ! File has probably been correctly identified

         inv_format_list(isite) = file_format
         inv_pss_list(isite) = pss_file
         inv_css_list(isite) = css_file
         inv_lat_list(isite) = site_lat
         inv_lon_list(isite) = site_lon

      end do

      return
   end subroutine assess_inventory_sites

   ! ==============================================================================================

   subroutine set_inventory_patch_type1(newpatch,pss_file_unit,ipa,ios,patch_name)

      ! --------------------------------------------------------------------------------------------
      ! This subroutine reads in a line of an inventory patch file (pss)
      ! And populates a new patch with that information.
      ! This routine specifically reads PSS files that are "Type 1" formatted
      !
      ! The file is formatted text, which contains 1 header line to label columns
      ! and then 1 line for each patch containing the following fields:
      !
      ! time	(year)     year of measurement
      ! patch	(string)   patch id string
      ! trk	(integer)  LU type index (0 non-forest, 1 secondary, 2 primary
      ! age	(years)    Time since this patch was disturbed (created)
      ! area	(fraction) Fraction of the site occupied by this patch
      ! --------------------------------------------------------------------------------------------

      use FatesSizeAgeTypeIndicesMod, only: get_age_class_index
      use EDtypesMod, only: AREA

      ! Arguments
      type(fates_patch_type),intent(inout), target   :: newpatch      ! Patch structure
      integer,intent(in)                          :: pss_file_unit ! Self explanatory
      integer,intent(in)                          :: ipa           ! Patch index (line number)
      integer,intent(out)                         :: ios           ! Return flag
      character(len=patchname_strlen),intent(out) :: patch_name    ! unique string identifier of patch

      ! Locals
      type(litter_type),pointer                   :: litt
      integer                                     :: el         ! index for elements
      real(r8)                                    :: p_time     ! Time patch was recorded
      integer                                     :: p_trk      ! Land Use index (see above descriptions)
      character(len=patchname_strlen)             :: p_name     ! unique string identifier of patch
      real(r8)                                    :: p_age      ! Patch age [years]
      real(r8)                                    :: p_area     ! Patch area [fraction]
      integer                                     :: icwd       ! index for counting CWD pools
      integer                                     :: ipft       ! index for counting PFTs
      real(r8)                                    :: pftfrac    ! the inverse of the total number of PFTs

      character(len=30),parameter    :: hd_fmt = &
            '(A5,2X,A20,2X,A4,2X,A5,2X,A17)'
      character(len=47),parameter    :: wr_fmt = &
            '(F5.2,2X,A20,2X,I4,2X,F5.2,2X,F17.14)'

      read(pss_file_unit,fmt=*,iostat=ios) p_time, p_name, p_trk, p_age, p_area

      if (ios/=0) return

      patch_name = trim(p_name)

      if( debug_inv) then
         write(*,fmt=hd_fmt) &
               ' time','               patch',' trk','  age','             area'
         write(*,fmt=wr_fmt) &
               p_time, p_name, p_trk, p_age, p_area
      end if

      ! Fill in the patch's memory structures

      newpatch%age       = p_age
      newpatch%age_class = get_age_class_index(newpatch%age)
      newpatch%area      = p_area*AREA

      ! The non-litter patch soil variables need to be sent to the HLM
      ! This is not being sent as of this message (RGK 06-2017)

      ! ---------------------------------------------------------------------
      ! The litter and CWD could be estimated from at least two methods
      ! 1) after reading in the cohort data, assuming a SS turnover rate
      ! 2) again assuming SS, back out the CWD and Litter flux rates into
      !    the SSC, STSC and FSC pools that balance with their decomp rates
      !    and then use those flux rates, to calculate the CWD and litter
      !    pool sizes based on another SS model where flux out balances with
      !    mortality and litter fluxes into the non-decomposing pool
      !
      ! This is significant science modeling and does not have a simple
      ! first hack solution. (RGK 06-2017)
      ! ----------------------------------------------------------------------

      do el=1,num_elements
         litt => newpatch%litter(el)

         call litt%InitConditions(init_leaf_fines=0._r8, &
              init_root_fines=0._r8, &
              init_ag_cwd=0._r8,     &
              init_bg_cwd=0._r8,     &
              init_seed=0._r8,   &
              init_seed_germ=0._r8)

      end do

      return
   end subroutine set_inventory_patch_type1


   ! ==============================================================================================

   subroutine set_inventory_cohort_type1(csite,bc_in,css_file_unit,npatches, &
                                           patch_pointer_vec,patch_name_vec,ios)

      ! --------------------------------------------------------------------------------------------
      ! This subroutine reads in a line of an inventory cohort file (css)
      ! And populates a new cohort with that information.
      ! This routine specifically reads CSS files that are "Type 1" formatted
      !
      ! The file formatted text, which contains 1 header line to label columns
      ! and then 1 line for each cohort containing the following fields:
      !
      ! FILE FORMAT:
      ! time	(year)     year of measurement
      ! patch	(string)   patch id string associated with this cohort
      ! dbh	(cm)       diameter at breast height. Optional, set height to negative if used
      ! height  (m)        height of vegetation in m. Optional, set dbh to negative if used
      ! pft     (integer)  the plant functional type index (must be consistent with param file)
      ! n 	(/m2)      The plant number density
      ! --------------------------------------------------------------------------------------------

      use FatesAllometryMod         , only : h_allom
      use FatesAllometryMod         , only : h2d_allom
      use FatesAllometryMod         , only : bagw_allom
      use FatesAllometryMod         , only : bbgw_allom
      use FatesAllometryMod         , only : bleaf
      use FatesAllometryMod         , only : bfineroot
      use FatesAllometryMod         , only : bsap_allom
      use FatesAllometryMod         , only : bdead_allom
      use FatesAllometryMod         , only : bstore_allom

      use EDCohortDynamicsMod , only : create_cohort
      use FatesInterfaceTypesMod   , only : numpft

      ! Arguments
      type(ed_site_type),intent(inout), target    :: csite         ! current site
      type(bc_in_type),intent(in)                 :: bc_in         ! boundary conditions
      integer, intent(in)                         :: css_file_unit     ! Self explanatory
      integer, intent(in)                         :: npatches      ! number of patches
      type(pp_array), intent(in)                  :: patch_pointer_vec(npatches)
      character(len=patchname_strlen), intent(in) :: patch_name_vec(npatches)
      integer,intent(out)                         :: ios           ! Return flag

      ! Locals
      class(prt_vartypes), pointer                :: prt_obj
      real(r8)                                    :: c_time        ! Time patch was recorded
      character(len=patchname_strlen)             :: p_name        ! The patch associated with this cohort
      character(len=patchname_strlen)             :: c_name        ! Cohort name
      real(r8)                                    :: c_dbh         ! diameter at breast height (cm)
      real(r8)                                    :: c_height      ! tree height (m)
      integer                                     :: c_pft         ! plant functional type index
      real(r8)                                    :: c_nplant      ! plant density (/m2)
      real(r8)                                    :: site_spread   ! initial guess of site spread
                                                                   ! should be quickly re-calculated
      integer,parameter                           :: rstatus = 0   ! recruit status

      type(fates_patch_type), pointer                :: cpatch        ! current patch pointer
      type(fates_cohort_type), pointer               :: temp_cohort   ! temporary patch (needed for allom funcs)
      integer                                     :: ipa           ! patch idex
      integer                                     :: iage
      integer                                     :: el
      integer                                     :: element_id
      logical                                     :: matched_patch ! check if cohort was matched w/ patch
      real(r8) :: c_agw    ! carbon biomass above ground non-leaf [kgC]
      real(r8) :: c_bgw    ! carbon biomass below ground non-leaf [kgC]
      real(r8) :: c_leaf   ! carbon biomass in leaves [kgC]
      real(r8) :: c_fnrt   ! carbon biomass in fine roots [kgC]
      real(r8) :: c_sapw   ! carbon biomass in sapwood [kgC]
      real(r8) :: c_struct ! carbon biomass in structure [kgC]
      real(r8) :: c_store  ! carbon biomass in storage [kgC]
      real(r8) :: a_sapw   ! area of sapwood at reference height [m2]
      real(r8) :: m_struct ! Generic (any element) mass for structure [kg]
      real(r8) :: m_leaf   ! Generic mass for leaf  [kg]
      real(r8) :: m_fnrt   ! Generic mass for fine-root  [kg]
      real(r8) :: m_sapw   ! Generic mass for sapwood [kg]
      real(r8) :: m_store  ! Generic mass for storage [kg]
      real(r8) :: m_repro  ! Generic mass for reproductive tissues [kg]
      real(r8) :: fnrt_drop_fraction ! Fine-root abscission fraction
      real(r8) :: stem_drop_fraction ! Stem abscission fraction
      integer  :: i_pft, ncohorts_to_create
     
      character(len=35),parameter    :: hd_fmt = &
           '(A7,2X,A20,2X,A5,2X,A6,2X,A4,2X,A9)'
      character(len=43),parameter    :: wr_fmt = &
           '(F7.1,2X,A20,2X,F5.2,2X,F6.2,2X,I4,2X,F9.6)'

      real(r8), parameter :: abnormal_large_nplant = 1000.0_r8  ! Used to catch bad values
      real(r8), parameter :: abnormal_large_dbh    = 500.0_r8   ! I've never heard of a tree > 3m
      real(r8), parameter :: abnormal_large_height = 500.0_r8   ! I've never heard of a tree > 500m tall
      integer,  parameter :: recruitstatus = 0
      logical, parameter :: old_type1_override = .false.

      if(old_type1_override) then
         ! time patch cohort dbh hite pft nplant bdead alive Avgrg 
         read(css_file_unit,fmt=*,iostat=ios) c_time, p_name, c_name, c_dbh, &
              c_height, c_pft, c_nplant
      else
         read(css_file_unit,fmt=*,iostat=ios) c_time, p_name, c_dbh, &
              c_height, c_pft, c_nplant
      end if
         
      if( debug_inv) then
         write(*,fmt=wr_fmt) &
              c_time, p_name, c_dbh, c_height, c_pft, c_nplant
      end if

      if (ios/=0) return

      ! Identify the patch based on the patch_name
      matched_patch = .false.
      do ipa=1,npatches
         if( trim(p_name) == trim(patch_name_vec(ipa))) then
            cpatch => patch_pointer_vec(ipa)%cpatch
            matched_patch = .true.
         end if
      end do

      if(.not.matched_patch)then
         write(fates_log(), *) 'could not match a cohort with a patch'
         write(fates_log(),fmt=hd_fmt) &
            '   time','               patch','  dbh','height',' pft','   nplant'
         write(fates_log(),fmt=wr_fmt) &
               c_time, p_name, c_dbh, c_height, c_pft, c_nplant
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      ! Run some sanity checks on the input data
      ! pft, nplant and dbh are the critical ones in this format specification
      ! -------------------------------------------------------------------------------------------

      if (c_pft > numpft ) then
         write(fates_log(), *) 'inventory pft: ',c_pft
         write(fates_log(), *) 'An inventory cohort file specified a pft index'
         write(fates_log(), *) 'greater than the maximum specified pfts numpft'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_pft < 0 ) then
         write(fates_log(), *) 'inventory pft: ',c_pft
         write(fates_log(), *) 'The inventory produced a cohort with <0 pft index'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_dbh < nearzero .and. c_height < nearzero) then
         write(fates_log(), *) 'inventory dbh: ', c_dbh
         write(fates_log(), *) 'and inventory height: ',c_height
         write(fates_log(), *) 'are both zero or negative. One must be positive.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_dbh > nearzero .and. c_height > nearzero) then
         write(fates_log(), *) 'inventory dbh: ', c_dbh
         write(fates_log(), *) 'and inventory height: ',c_height
         write(fates_log(), *) 'are both positive. One must be zero or negative.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      
      if (c_dbh > abnormal_large_dbh ) then
         write(fates_log(), *) 'inventory dbh: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with very large diameter [cm]'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_height > abnormal_large_height ) then
         write(fates_log(), *) 'inventory height: ', c_height
         write(fates_log(), *) 'The inventory produced a cohort with very large height [m]'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_nplant <=0 ) then
         write(fates_log(), *) 'inventory nplant: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with <= 0 density /m2'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_nplant > abnormal_large_nplant ) then
         write(fates_log(), *) 'inventory nplant: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with very large density /m2'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_pft .eq. 0 ) then
         if(debug_inv)then
            write(fates_log(), *) 'inventory pft: ',c_pft
            write(fates_log(), *) 'SPECIAL CASE TRIGGERED: PFT == 0 and therefore this subroutine'
            write(fates_log(), *) 'will assign a cohort with n = n_orig/numpft to every cohort in range 1 to numpft'
         end if
         ncohorts_to_create = numpft
      else
         ncohorts_to_create = 1
      end if

      do i_pft = 1,ncohorts_to_create
         allocate(temp_cohort)   ! A temporary cohort is needed because we want to make

         if (c_pft .ne. 0 ) then
             ! normal case: assign each cohort to its specified PFT
             temp_cohort%pft         = c_pft
         else
             ! special case, make an identical cohort for each PFT
             temp_cohort%pft         = i_pft
         endif

         temp_cohort%n           = c_nplant * cpatch%area / real(ncohorts_to_create,r8)
         
         temp_cohort%crowndamage = 1  ! assume undamaged 

         if( c_dbh> 0._r8)then
            temp_cohort%dbh         = c_dbh
            call h_allom(c_dbh,temp_cohort%pft,temp_cohort%height)
         else
            temp_cohort%height  = c_height
            call h2d_allom(c_height,temp_cohort%pft,temp_cohort%dbh)
         end if

         temp_cohort%canopy_trim = 1.0_r8

         ! Determine the phenology status and the elongation factors.
         fnrt_drop_fraction = prt_params%phen_fnrt_drop_fraction(temp_cohort%pft)
         stem_drop_fraction = prt_params%phen_stem_drop_fraction(temp_cohort%pft)

         if( prt_params%season_decid(temp_cohort%pft) == itrue .and. &
              any(csite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold])) then
            ! Cold deciduous and season is for leaves off. Set leaf status and 
            ! elongation factors accordingly
            temp_cohort%efleaf_coh = 0.0_r8
            temp_cohort%effnrt_coh = 1._r8 - fnrt_drop_fraction
            temp_cohort%efstem_coh = 1._r8 - stem_drop_fraction

            temp_cohort%status_coh = leaves_off

         elseif ( any(prt_params%stress_decid(temp_cohort%pft) == [ihard_stress_decid,isemi_stress_decid])) then
            ! Drought deciduous.  For the default approach, elongation factor is either
            ! zero (full abscission) or one (fully flushed), but this can also be a
            ! fraction in other approaches. Here we assume that leaves are "on" (i.e.
            ! either fully flushed or growing) if elongation factor is not 0 for the 
            ! initial conditions.
            ! 
            ! For tissues other than leaves, the actual drop fraction is a combination
            ! of the elongation factor (e) and the drop fraction (x), which will ensure
            ! that the remaining tissue biomass will be exactly e when x=1, and exactly
            ! the original biomass when x = 0.
            temp_cohort%efleaf_coh = csite%elong_factor(temp_cohort%pft)
            temp_cohort%effnrt_coh = 1.0_r8 - (1.0_r8 - temp_cohort%efleaf_coh ) * fnrt_drop_fraction
            temp_cohort%efstem_coh = 1.0_r8 - (1.0_r8 - temp_cohort%efleaf_coh ) * stem_drop_fraction

            if (temp_cohort%efleaf_coh > 0.0_r8) then
               ! Assume leaves are growing even if they are not fully flushed.
               temp_cohort%status_coh = leaves_on
            else
               ! Leaves are off (abscissing).
               temp_cohort%status_coh = leaves_off
            end if
         else
            ! Evergreen, or deciduous PFT during the growing season. Assume tissues are fully flushed.
            temp_cohort%efleaf_coh = 1.0_r8
            temp_cohort%effnrt_coh = 1.0_r8
            temp_cohort%efstem_coh = 1.0_r8

            temp_cohort%status_coh = leaves_on
         end if

         call bagw_allom(temp_cohort%dbh,temp_cohort%pft, &
              temp_cohort%crowndamage, temp_cohort%efstem_coh, c_agw)
         ! Calculate coarse root biomass from allometry
         call bbgw_allom(temp_cohort%dbh,temp_cohort%pft, temp_cohort%efstem_coh, c_bgw)

         ! Calculate the leaf biomass (calculates a maximum first, then applies canopy trim
         ! and sla scaling factors)
         call bleaf(temp_cohort%dbh,temp_cohort%pft,temp_cohort%crowndamage,&
              temp_cohort%canopy_trim, temp_cohort%efleaf_coh, c_leaf)
         
         ! Calculate fine root biomass

         temp_cohort%l2fr = prt_params%allom_l2fr(temp_cohort%pft)
         call bfineroot(temp_cohort%dbh,temp_cohort%pft,temp_cohort%canopy_trim,temp_cohort%l2fr, &
              temp_cohort%effnrt_coh, c_fnrt)

         ! Calculate sapwood biomass
         call bsap_allom(temp_cohort%dbh,temp_cohort%pft,temp_cohort%crowndamage, &
              temp_cohort%canopy_trim, temp_cohort%efstem_coh, a_sapw, c_sapw)
         
         call bdead_allom( c_agw, c_bgw, c_sapw, temp_cohort%pft, c_struct )
         call bstore_allom(temp_cohort%dbh, temp_cohort%pft, temp_cohort%crowndamage,temp_cohort%canopy_trim, c_store)

         prt_obj => null()
         call InitPRTObject(prt_obj)

         !  (Keeping as an example)
         ! Allocate running mean functions
         !allocate(temp_cohort%tveg_lpa)
         !call temp_cohort%tveg_lpa%InitRMean(ema_lpa,init_value=cpatch%tveg_lpa%GetMean())
         
         do el = 1,num_elements

            element_id = element_list(el)

            ! If this is carbon12, then the initialization is straight forward
            ! otherwise, we use stoichiometric ratios
            select case(element_id)
            case(carbon12_element)

               m_struct = c_struct
               m_leaf   = c_leaf
               m_fnrt   = c_fnrt
               m_sapw   = c_sapw
               m_store  = c_store
               m_repro  = 0._r8

            case(nitrogen_element)

               ! For inventory runs, initialize nutrient contents half way between max and min stoichiometries
               m_struct = c_struct * &
                    prt_params%nitr_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(struct_organ))

               m_leaf   = c_leaf * &
                    prt_params%nitr_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(leaf_organ))

               m_fnrt   = c_fnrt * &
                    prt_params%nitr_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(fnrt_organ))

               m_sapw   = c_sapw * &
                    prt_params%nitr_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(sapw_organ))

               m_repro  = 0._r8

               m_store = StorageNutrientTarget(temp_cohort%pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)

            case(phosphorus_element)

               m_struct = c_struct * &
                    prt_params%phos_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(struct_organ))

               m_leaf   = c_leaf * &
                    prt_params%phos_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(leaf_organ))

               m_fnrt   = c_fnrt * &
                    prt_params%phos_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(fnrt_organ))
                    
               m_sapw   = c_sapw * &
                    prt_params%phos_stoich_p1(temp_cohort%pft,prt_params%organ_param_id(sapw_organ))

               m_repro  = 0._r8

               m_store = StorageNutrientTarget(temp_cohort%pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)

            end select

            select case(hlm_parteh_mode)
            case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp )

               ! Equally distribute leaf mass into available age-bins
               do iage = 1,nleafage
                  call SetState(prt_obj,leaf_organ, element_id,m_leaf/real(nleafage,r8),iage)
               end do

               call SetState(prt_obj,fnrt_organ, element_id, m_fnrt)
               call SetState(prt_obj,sapw_organ, element_id, m_sapw)
               call SetState(prt_obj,store_organ, element_id, m_store)
               call SetState(prt_obj,struct_organ, element_id, m_struct)
               call SetState(prt_obj,repro_organ, element_id, m_repro)

            case default
               write(fates_log(),*) 'Unspecified PARTEH module during inventory intitialization'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end select

         end do

         call prt_obj%CheckInitialConditions()

         call create_cohort(csite, cpatch, temp_cohort%pft, temp_cohort%n, temp_cohort%height, &
              temp_cohort%coage, temp_cohort%dbh, &
              prt_obj, temp_cohort%efleaf_coh, temp_cohort%effnrt_coh, &
              temp_cohort%efstem_coh, temp_cohort%status_coh, rstatus, &
              temp_cohort%canopy_trim,temp_cohort%c_area, &
              1, temp_cohort%crowndamage, csite%spread, bc_in)

         deallocate(temp_cohort) ! get rid of temporary cohort

      end do
      call cpatch%ValidateCohorts()

      return
    end subroutine set_inventory_cohort_type1

   ! ====================================================================================

   subroutine write_inventory_type1(currentSite)

       ! --------------------------------------------------------------------------------
       ! This subroutine writes the cohort/patch inventory type files in the "type 1"
       ! format.
       ! The files will have a lat/long tag added to their name, and will be
       ! generated in the run folder.
       ! JFN Oct 2023 - updated to get rid of unused ED columns      
       ! --------------------------------------------------------------------------------

       use shr_file_mod, only        : shr_file_getUnit
       use shr_file_mod, only        : shr_file_freeUnit

       ! Arguments
       type(ed_site_type), target :: currentSite

       ! Locals
       type(fates_patch_type), pointer          :: currentpatch
       type(fates_cohort_type), pointer         :: currentcohort

       character(len=128)                    :: pss_name_out         ! output file string
       character(len=128)                    :: css_name_out         ! output file string
       integer                               :: pss_file_out
       integer                               :: css_file_out
       integer                               :: ilat_int,ilat_dec    ! for output string parsing
       integer                               :: ilon_int,ilon_dec    ! for output string parsing
       character(len=32)                     :: patch_str
       character(len=32)                     :: cohort_str
       integer                               :: ipatch
       integer                               :: icohort
       character(len=1)                      :: ilat_sign,ilon_sign

       ! Generate pss/css file name based on the location of the site
       ilat_int = abs(int(currentSite%lat))
       ilat_dec = int(100000*(abs(currentSite%lat) - real(ilat_int,r8)))
       ilon_int = abs(int(currentSite%lon))
       ilon_dec = int(100000*(abs(currentSite%lon) - real(ilon_int,r8)))

       if(currentSite%lat>=0._r8)then
           ilat_sign = 'N'
       else
           ilat_sign = 'S'
       end if
       if(currentSite%lon>=0._r8)then
           ilon_sign = 'E'
       else
           ilon_sign = 'W'
       end if
 
       write(pss_name_out,'(A8,I2.2,A1,I5.5,A1,A1,I3.3,A1,I5.5,A1,A4)') &
             'pss_out_',ilat_int,'.',ilat_dec,ilat_sign,'_',ilon_int,'.',ilon_dec,ilon_sign,'.txt'
       write(css_name_out,'(A8,I2.2,A1,I5.5,A1,A1,I3.3,A1,I5.5,A1,A4)') &
             'css_out_',ilat_int,'.',ilat_dec,ilat_sign,'_',ilon_int,'.',ilon_dec,ilon_sign,'.txt'

       pss_file_out       = shr_file_getUnit()
       css_file_out       = shr_file_getUnit()

       open(unit=pss_file_out,file=trim(pss_name_out), status='UNKNOWN',action='WRITE',form='FORMATTED')
       open(unit=css_file_out,file=trim(css_name_out), status='UNKNOWN',action='WRITE',form='FORMATTED')

       write(pss_file_out,*) 'time patch trk age area'
       write(css_file_out,*) 'time patch dbh height pft nplant'

       ipatch=0
       currentpatch => currentSite%youngest_patch
       do while(associated(currentpatch))
           ipatch=ipatch+1

           write(patch_str,'(A7,i4.4,A)') '<patch_',ipatch,'>'

           write(pss_file_out,*) '0000 ',trim(patch_str),' 2 ',currentPatch%age,currentPatch%area/AREA

           icohort=0
           currentcohort => currentpatch%tallest
           do while(associated(currentcohort))
               icohort=icohort+1
               write(css_file_out,*) '0000 ',trim(patch_str), &
                     currentCohort%dbh,-3.0_r8,currentCohort%pft,currentCohort%n/currentPatch%area

               currentcohort => currentcohort%shorter
           end do
           currentPatch => currentpatch%older
       enddo

       close(css_file_out)
       close(pss_file_out)

       call shr_file_freeUnit(css_file_out)
       call shr_file_freeUnit(pss_file_out)

   end subroutine write_inventory_type1

end module FatesInventoryInitMod
