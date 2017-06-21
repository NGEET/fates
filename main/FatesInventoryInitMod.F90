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
   ! This code borrows heavily in concept from what is done in ED2.  We will also do our best to
   ! maintain compatibility with the PSS/CSS file formats that were used in ED2.
   ! See: https://github.com/EDmodel/ED2/blob/master/ED/src/io/ed_read_ed10_20_history.f90
   ! At the time of writing this ED2 is unlicensed, and only concepts were borrowed with no direct
   ! code copied.
   !-----------------------------------------------------------------------------------------------

   ! CIME GLOBALS
   
   use shr_log_mod , only      : errMsg => shr_log_errMsg

   ! FATES GLOBALS
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log

   use EDTypesMod       , only : ed_site_type, ed_patch_type, ed_cohort_type, area
  
   implicit none
   
   character(len=*), parameter :: inv_file_list = 'inventory_file_list.txt'

   character(len=*), parameter, private :: sourcefile = &
         __FILE__

   logical, parameter          :: debug_inv = .true.

   public :: count_inventory_sites
   public :: inv_file_list
   public :: assess_inventory_sites
   public :: set_inventory_edpatch_type1

contains

   ! ==============================================================================================


   function count_inventory_sites(file_unit) result(nsites)

      integer, intent(in) :: file_unit

      character(len=512) :: header_str
      character(len=512) :: site_str
      integer          :: ios
      real(r8)         :: site_lat
      real(r8)         :: site_lon
      character(len=256) :: pss_file
      character(len=256) :: css_file

      integer          :: nsites


      ! Set the file position to the top of the file
      ! Read in the header line
      ! Read through sites and check coordinates and file existence
      rewind(unit=file_unit)
      read(file_unit,fmt='(A)') header_str
      nsites = 0
      do
         read(file_unit,fmt='(A)',iostat=ios) site_str
         if (ios/=0) exit
         nsites = nsites + 1
      end do

      return
   end function count_inventory_sites

   ! ==============================================================================================

   subroutine assess_inventory_sites(file_unit,nsites, inv_format_list, &
                                     inv_pss_list,inv_css_list, &
                                     inv_lat_list,inv_lon_list)

      integer, intent(in) :: file_unit
      integer, intent(in) :: nsites
      integer, intent(inout)           :: inv_format_list(nsites)
      character(len=256),intent(inout) :: inv_pss_list(nsites)
      character(len=256),intent(inout) :: inv_css_list(nsites)
      real(r8),intent(inout)           :: inv_lat_list(nsites)
      real(r8),intent(inout)           :: inv_lon_list(nsites)

      character(len=512) :: header_str
      character(len=512) :: site_str
      integer            :: isite
      integer            :: ios
      character(len=256) :: pss_file
      character(len=256) :: css_file
      real(r8)           :: site_lat
      real(r8)           :: site_lon
      integer            :: iblnk
      integer            :: file_format
      logical            :: lex
 
      rewind(unit=file_unit)
      read(file_unit,fmt='(4A)') header_str

      do isite=1,nsites

         ! Read in the whole line
         read(file_unit,fmt='(a)',iostat=ios) site_str

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

         if ( site_lon < -180.0_r8 .or. site_lon > 360.0_r8 ) then
            write(fates_log(), *) 'read invalid longitude: ',site_lon,' from inventory site list'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         inquire(file=trim(pss_file),exist=lex)
         if( .not.lex ) then
            write(fates_log(), *) 'the following pss file could not be found:'
            write(fates_log(), *) trim(pss_file)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         inquire(file=trim(css_file),exist=lex)
         if( .not.lex ) then
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
      


   end subroutine assess_inventory_sites

   ! ==============================================================================================

   subroutine set_inventory_edpatch_type1(newpatch,file_unit,ipa,ios)

     ! --------------------------------------------------------------------------------------------
     ! This subroutine reads in a line of an inventory patch file (pss)
     ! And populates a new patch with that information.
     ! This routine specifically reads PSS files that are "Type 1" formatted
     ! 
     ! FILE FORMAT:
     ! time	(year)     year of measurement
     ! patch	(string)   patch id string
     ! trk	(integer)  LU type index (0 non-forest, 1 secondary, 2 primary
     ! age	(years)    Time since this patch was disturbed (created)
     ! area	(fraction) Fraction of the site occupied by this patch
     ! water	(NA)       Water content of soil (NOT USED)
     ! fsc	(kg/m2)    Fast Soil Carbon
     ! stsc	(kg/m2)    Structural Soil Carbon
     ! stsl	(kg/m2)    Structural Soil Lignan
     ! ssc	(kg/m2)    Slow Soil Carbon
     ! psc	(NA)       Passive Soil Carbon (NOT USED)
     ! msn	(kg/m2)    Mineralized Soil Nitrogen
     ! fsn      (kg/m2)    Fast Soil Nitrogen
     ! --------------------------------------------------------------------------------------------

     use EDTypesMod, only: get_age_class_index
     use EDtypesMod, only: AREA
     use EDTypesMod, only: numpft_ed
     use EDTypesMod, only: ncwd
     use SFParamsMod , only : SF_val_CWD_frac
     use EDParamsMod , only : ED_val_ag_biomass

     
     
     ! Arguments
     type(ed_patch_type),intent(inout), target :: newpatch  ! Patch structure
     integer,intent(in)                       :: file_unit ! Self explanatory
     integer,intent(in)                       :: ipa       ! Patch index (line number)
     integer,intent(out)                      :: ios       ! Return flag

     real(r8)                        :: p_time  ! Time patch was recorded
     character(len=64)               :: p_name  ! Unique ID string defining patch
     real(r8)                        :: p_trk   ! Land Use index
                                                ! 0 = Agriculture, 1 = Secondary Forest
                                                ! 2 = Primary Forest, 3 = Forest Plantation
                                                ! 4 = Burnt Patch, 5 = Abandoned (secondary growth)
                                                ! 6 = Logged Forest
     real(r8)                        :: p_age   ! Patch age [years]
     real(r8)                        :: p_area  ! Patch area [fraction]
     real(r8)                        :: p_water ! Patch water (unused)
     real(r8)                        :: p_fsc   ! Patch fast soil carbon
     real(r8)                        :: p_stsc  ! Patch structural soil carbon
     real(r8)                        :: p_stsl  ! Patch structural soil lignans
     real(r8)                        :: p_ssc   ! Patch slow soil carbon
     real(r8)                        :: p_psc   ! Patch P soil carbon
     real(r8)                        :: p_msn   ! Patch mean soil nitrogen
     real(r8)                        :: p_fsn   ! Patch fast soil nitrogen
     integer                         :: icwd    ! index for counting CWD pools
     integer                         :: ipft    ! index for counting PFTs
     real(r8)                        :: pftfrac ! the inverse of the total number of PFTs

     character(len=128),parameter    :: wr_fmt = &
          '(F5.2,2X,A4,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2)'

     real(r8), parameter             :: cwdfrac = 0.95  ! CWD is 95% of structural biomass (GUESS, BAD ONE)
     real(r8), parameter             :: leaffrac = 0.5  ! leaf litter is this fraction of total

     
     read(file_unit,fmt=*,iostat=ios) p_time, p_name, p_trk, p_age, p_area, &
                                      p_water,p_fsc, p_stsc, p_stsl, p_ssc, &
                                      p_psc, p_msn, p_fsn
     
     if (ios/=0) return
     
     if( debug_inv) then
        write(*,fmt=wr_fmt) &
             p_time, p_name, p_trk, p_age, p_area, &
             p_water,p_fsc, p_stsc, p_stsl, p_ssc, &
             p_psc, p_msn, p_fsn
     end if
     
     ! Fill in the patch's memory structures

     newpatch%age       = p_age
     newpatch%age_class = get_age_class_index(newpatch%age)
     newpatch%area      = p_area*AREA

     ! The non-litter patch soil variables need to be sent to the HLM
     ! This is not being sent as of this message (RGK 06-2017)

     ! Estimate CWD and litter pools from p_stsc (twig,s branch,l branch, trunk)

     ! Lets start out assuming that CWD is about 95% of the total structural 
     ! carbon pool from non-living biomass

     do icwd = 1, ncwd
        newpatch%cwd_ag(icwd) = p_stsc*cwdfrac*ED_val_ag_biomass * SF_val_CWD_frac(icwd)
        newpatch%cwd_bg(icwd) = p_stsc*cwdfrac*(1.0_r8 - ED_val_ag_biomass) * SF_val_CWD_frac(icwd)
     end do

     pftfrac = 1.0_r8/dble(numpft_ed)

     do ipft = 1, numpft_ed
        newpatch%leaf_litter(ipft) = p_stsc*(1.0_r8-cwdfrac)*leaffrac*pftfrac
        newpatch%root_litter(ipft) = p_stsc*(1.0_r8-cwdfrac)*(1.0_r8-leaffrac)*pftfrac
     end do

     return

   end subroutine set_inventory_edpatch_type1




end module FatesInventoryInitMod
