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
   use shr_file_mod, only      : shr_file_getUnit
   use shr_log_mod , only      : errMsg => shr_log_errMsg

   ! FATES GLOBALS
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log
  
   implicit none
   
   character(len=*), parameter :: inv_file_list = 'inventory_file_list.txt'

   character(len=*), parameter, private :: sourcefile = &
         __FILE__

   public :: count_inventory_sites
   public :: get_inventory_file_unit
   public :: inv_file_list
   public :: assess_inventory_sites

contains

   



   ! ==============================================================================================

   ! This is is as wrappery as wrapper get
   function get_inventory_file_unit() result(file_unit)
      integer :: file_unit
      file_unit = shr_file_getUnit()
      return
   end function get_inventory_file_unit


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

   subroutine assess_inventory_sites(file_unit,nsites, &
                                     inv_pss_list,inv_css_list, &
                                     inv_lat_list,inv_lon_list)

      integer, intent(in) :: file_unit
      integer, intent(in) :: nsites
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
      logical            :: lex
 
      rewind(unit=file_unit)
      read(file_unit,fmt='(4A)') header_str
      print*,trim(header_str)
      do isite=1,nsites

         ! Read in the whole line
         read(file_unit,fmt='(a)',iostat=ios) site_str

         ! Parse the first entry from the line (latitude)
         read(site_str,*) site_lat

         ! Parse the second entry from the line (longitude)
         site_str = adjustl(site_str)
         iblnk = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,*) site_lon

         site_str = adjustl(site_str)
         iblnk    = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         iblnk    = index(site_str,' ')
         read(site_str(:iblnk),fmt='(1A)') pss_file

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

         inv_pss_list(isite) = pss_file
         inv_css_list(isite) = css_file
         inv_lat_list(isite) = site_lat
         inv_lon_list(isite) = site_lon

      end do
      


   end subroutine assess_inventory_sites


end module FatesInventoryInitMod
