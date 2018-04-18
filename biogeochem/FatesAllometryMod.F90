!===============================================================================
!
! FatesAllometryMod.F90
!
! A library of functions that calculate plant allometry and their
! derivatives.  Most relationships are related to diameter [cm].  All
! derivatives with respect to change in diameter have same units.
! In some cases multiple areguments will be provided, yet
! those arguments in those cases are trivially determined from diameter.
!
! Each function presented in this library is written to also return the
! derivative with respect to diameter, if the argument is provided.
!
! The name convention of the functions follows the form d...  Which
! indicates "diameter to ...".  Allometries for the following variables are
! calculated:
! h:      height [m]
! bagw:   biomass above ground woody tissues [kgC]
!         this is an all encompassing definition of "woody", which 
!         is intended to include all non-leaf woody or fibrous
!         tissues, including sap and structural materials
! blmax:  biomass in leaves when leaves are "on allometry"
!         this also is the "pre-trimmed" value, which is the maximum
!         or potential leaf mass a tree may have [kgC]
! bbgw:   biomass below ground woody tissues [kgC] 
!         similar to AGBW, this essentially encompasses
!         all non-fineroot belowground tissues
! bfrmax: biomass in fine roots when "on allometry" [kgC]
! bsap: biomass in sapwood (above and below) [kgC]
! bdead: biomass (above and below) in the structural pool [kgC]
!
!
! The following function switches are rused
!
! allom_hmode, integer, height allometry function type
! allom_lmode, integer, maximum leaf allometry function type
! allom_rmode, integer, maximum root allometry function type
! allom_amode, integer, AGB allometry function type
! allom_cmode, integer, coarse root allometry function type
! allom_smode, integer, sapwood allometry function type
! allom_stmode, integer, storage allometry function type
!
! The following parameters (traits) are used
!
! wood_density, mean stem wood specific gravity (heart,sap,bark)
! allom_latosa_int, leaf area to sap area ratio, intercept [m2/cm2]
! allom_latosa_slp, leaf area to sap area ratio, slope on diameter [m2/cm2/cm]
! c2b, real, carbon to biomass multiplier (~2.0)
! allom_l2fr, fine root biomass per leaf biomass ratio [kgC/kgC]
! allom_agb_frac, the fraction of stem above ground [-]
! allom_d2h1, parameter 1 for d2h allometry (intercept)
! allom_d2h2, parameter 2 for d2h allometry (slope)
! allom_d2h3, parameter 3 for d2h allometry (optional)
! allom_d2bl1, parameter 1 for d2bl allometry (intercept)
! allom_d2bl2, parameter 2 for d2bl allometry (slope)
! allom_d2bl3, parameter 3 for d2bl allometry (optional)
! allom_agb1
! allom_agb2
! allom_agb3
! allom_dbh_maxheight, dbh at maximum height [cm]
! h_min, the height associated with newly recruited plant [m]
!
! Note - i4 types are expressed explicitly to accomodate unit testing calls
!        to this module
!
! Explanation of pools and their control volumes
!
! ------------------------------------------------------------------------------
!
! total biomass = bleaf + bfineroot + agbw + bgbw  
!                  ... or ...
! total biomass = bleaf + bfineroot + bdead + bsap
!                  ... and ...
! bdead         = agbw + bgbw - bsap
!
! ------------------------------------------------------------------------------
!
! Initial Implementation: Ryan Knox July 2017
!
!===============================================================================

module FatesAllometryMod

  ! If this is a unit-test, these globals will be provided by a wrapper

  use EDPFTvarcon      , only : EDPftvarcon_inst
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : g_per_kg 
  use FatesConstantsMod, only : cm2_per_m2
  use FatesConstantsMod, only : kg_per_Megag
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use FatesGlobals     , only : fates_log
  use FatesGlobals     , only : endrun => fates_endrun
  use EDTypesMod       , only : nlevleaf, dinc_ed

  implicit none

  private
  public :: h2d_allom     ! Generic height to diameter wrapper
  public :: h_allom       ! Generic diameter to height wrapper
  public :: bagw_allom    ! Generic AGWB (above grnd. woody bio) wrapper
  public :: blmax_allom   ! Generic maximum leaf biomass wrapper
  public :: bleaf         ! Generic actual leaf biomass wrapper
  public :: storage_fraction_of_target ! storage as fraction of leaf biomass
  public :: tree_lai      ! Calculate tree-level LAI from actual leaf biomass
  public :: tree_sai      ! Calculate tree-level SAI from target leaf biomass
  public :: bsap_allom    ! Generic sapwood wrapper
  public :: bbgw_allom    ! Generic coarse root wrapper
  public :: bfineroot     ! Generic actual fine root biomass wrapper
  public :: bdead_allom   ! Generic bdead wrapper
  public :: carea_allom   ! Generic crown area wrapper
  public :: bstore_allom  ! Generic maximum storage carbon wrapper
  public :: StructureResetOfDH ! Method to set DBH to sync with structure biomass
  public :: CheckIntegratedAllometries

  logical         , parameter :: verbose_logging = .false.
  character(len=*), parameter :: sourcefile = __FILE__

  ! If testing b4b with older versions, do not remove sapwood
  ! Our old methods with saldarriaga did not remove sapwood from the
  ! bdead pool.  But newer allometries are providing total agb
  ! which includes sapwood. Although a small quantity, it needs to be removed
  ! from the agb pool.
  ! Additionally, our calculation of sapwood biomass may be missing some unit conversions

contains

  ! ============================================================================
  ! Parameter Checks
  ! ============================================================================

  ! Checks to make sure parameters are not within expected ranges for each
  ! functions

  ! Check to make sure Martinez-Cano height cap is not on, or explicitly allowed
   

  ! ===========================================================================
  ! Helper Routines
  ! ===========================================================================

   
   ! ============================================================================
   
  subroutine CheckIntegratedAllometries(dbh,ipft,canopy_trim, &
       bl,bfr,bsap,bstore,bdead, &
       grow_leaf, grow_fr, grow_sap, grow_store, grow_dead, &
       max_err, l_pass)

     ! This routine checks the error on the carbon allocation
     ! integration step.  The integrated quantities should
     ! be a close match on the diagnosed quantities.
     ! We don't have to worry about small accumulating biases,
     ! or small errors because the scheme automatically pushes
     ! all carbon pools towards the allometric value corresponding
     ! to the diameter on each step, prior to performing an integration.

     real(r8),intent(in) :: dbh    ! diameter of plant [cm]
     integer,intent(in)  :: ipft   ! plant functional type index
     real(r8),intent(in) :: canopy_trim ! trimming function
     real(r8),intent(in) :: bl     ! integrated leaf biomass [kgC]
     real(r8),intent(in) :: bfr    ! integrated fine root biomass [kgC]
     real(r8),intent(in) :: bsap   ! integrated sapwood biomass [kgC]
     real(r8),intent(in) :: bstore ! integrated storage biomass [kgC]
     real(r8),intent(in) :: bdead  ! integrated structural biomass [kgc]
     logical,intent(in)  :: grow_leaf ! on-off switch for leaf growth
     logical,intent(in)  :: grow_fr   ! on-off switch for root growth
     logical,intent(in)  :: grow_sap  ! on-off switch for sapwood
     logical,intent(in)  :: grow_store! on-off switch for storage
     logical,intent(in)  :: grow_dead ! on-off switch for structure
     real(r8),intent(in) :: max_err   ! maximum allowable error

     logical,intent(out) :: l_pass   ! Error flag (pass=true,no-pass=false)
     
     real(r8) :: height            ! diagnosed height [m]
     real(r8) :: bl_diag           ! diagnosed leaf biomass [kgC]
     real(r8) :: bfr_diag          ! diagnosed fine-root biomass [kgC]
     real(r8) :: bsap_diag         ! diagnosed sapwood biomass [kgC]
     real(r8) :: bdead_diag        ! diagnosed structural biomass [kgC]
     real(r8) :: bstore_diag       ! diagnosed storage biomass [kgC]
     real(r8) :: bagw_diag         ! diagnosed agbw [kgC]
     real(r8) :: bbgw_diag         ! diagnosed below ground wood [kgC]



     l_pass = .true.  ! Default assumption is that step passed

     if (grow_leaf) then
        call bleaf(dbh,ipft,canopy_trim,bl_diag)
        if( abs(bl_diag-bl) > max_err ) then
           if(verbose_logging) then
              write(fates_log(),*) 'disparity in integrated/diagnosed leaf carbon'
              write(fates_log(),*) 'resulting from the on-allometry growth integration step'
              write(fates_log(),*) 'bl (integrated): ',bl
              write(fates_log(),*) 'bl (diagnosed): ',bl_diag
              write(fates_log(),*) 'relative error: ',abs(bl_diag-bl)/bl_diag
           end if
           l_pass = .false.
        end if
     end if
        
     if (grow_fr) then
        call bfineroot(dbh,ipft,canopy_trim,bfr_diag)
        if( abs(bfr_diag-bfr) > max_err ) then
           if(verbose_logging) then
              write(fates_log(),*) 'disparity in integrated/diagnosed fineroot carbon'
              write(fates_log(),*) 'resulting from the on-allometry growth integration step'
              write(fates_log(),*) 'bfr (integrated): ',bfr
              write(fates_log(),*) 'bfr (diagnosed): ',bfr_diag
              write(fates_log(),*) 'relative error: ',abs(bfr_diag-bfr)/bfr_diag
           end if
           l_pass = .false.
        end if
     end if

     if (grow_sap) then
        call bsap_allom(dbh,ipft,canopy_trim,bsap_diag)
        if( abs(bsap_diag-bsap) > max_err ) then
           if(verbose_logging) then
              write(fates_log(),*) 'disparity in integrated/diagnosed sapwood carbon'
              write(fates_log(),*) 'resulting from the on-allometry growth integration step'
              write(fates_log(),*) 'bsap (integrated): ',bsap
              write(fates_log(),*) 'bsap (diagnosed): ',bsap_diag
              write(fates_log(),*) 'relative error: ',abs(bsap_diag-bsap)/bsap_diag
           end if
           l_pass = .false.
        end if
     end if
        
     if (grow_store) then
        call bstore_allom(dbh,ipft,canopy_trim,bstore_diag)
        if( abs(bstore_diag-bstore) > max_err ) then
           if(verbose_logging) then
              write(fates_log(),*) 'disparity in integrated/diagnosed storage carbon'
              write(fates_log(),*) 'resulting from the on-allometry growth integration step'
              write(fates_log(),*) 'bsap (integrated): ',bstore
              write(fates_log(),*) 'bsap (diagnosed): ',bstore_diag
              write(fates_log(),*) 'relative error: ',abs(bstore_diag-bstore)/bstore_diag
           end if
           l_pass = .false.
        end if
     end if

     if (grow_dead) then
        call bsap_allom(dbh,ipft,canopy_trim,bsap_diag)
        call bagw_allom(dbh,ipft,bagw_diag)
        call bbgw_allom(dbh,ipft,bbgw_diag)
        call bdead_allom( bagw_diag, bbgw_diag, bsap_diag, ipft, bdead_diag )        
        if( abs(bdead_diag-bdead) > max_err ) then
           if(verbose_logging) then
              write(fates_log(),*) 'disparity in integrated/diagnosed structural carbon'
              write(fates_log(),*) 'resulting from the on-allometry growth integration step'
              write(fates_log(),*) 'bdead (integrated): ',bdead
              write(fates_log(),*) 'bdead (diagnosed): ',bdead_diag
              write(fates_log(),*) 'relative error: ',abs(bdead_diag-bdead)/bdead_diag
           end if
           l_pass = .false.
        end if
     end if
        
     return
   end subroutine CheckIntegratedAllometries



  ! ============================================================================
  ! Generic height to diameter interface
  ! ============================================================================

  subroutine h2d_allom(h,ipft,d,dddh)

    real(r8),intent(in)           :: h     ! height of plant [m]
    integer(i4),intent(in)        :: ipft  ! PFT index
    real(r8),intent(out)          :: d     ! plant diameter [cm]
    real(r8),intent(out),optional :: dddh  ! change in diameter per height [cm/m]

    associate(  p1          => EDPftvarcon_inst%allom_d2h1(ipft), &
                p2          => EDPftvarcon_inst%allom_d2h2(ipft), &
                p3          => EDPftvarcon_inst%allom_d2h3(ipft), &
                allom_hmode => EDPftvarcon_inst%allom_hmode(ipft))

      select case(int(allom_hmode))
      case (1) ! Obrien et al. 199X BCI
         call h2d_obrien(h,p1,p2,d,dddh)
      case (2)  ! poorter 2006
         call h2d_poorter2006(h,p1,p2,p3,d,dddh)
      case (3) ! 2 parameter power function
         call h2d_2pwr(h,p1,p2,d,dddh)
      case (4) ! chave 2014
         call h2d_chave2014(h,p1,p2,p3,d,dddh)
      case (5) ! Martinez-Cano
         call h2d_martcano(h,p1,p2,p3,d,dddh)
      case DEFAULT
         write(fates_log(),*) 'An undefined h2d allometry was specified: ',allom_hmode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      
    end associate
    return
  end subroutine h2d_allom

  ! ============================================================================
  ! Generic height interface
  ! ============================================================================

  subroutine h_allom(d,ipft,h,dhdd)

    ! Arguments
    real(r8),intent(in)     :: d     ! plant diameter [cm]
    integer(i4),intent(in)  :: ipft  ! PFT index
    real(r8),intent(out)    :: h     ! plant height [m]
    real(r8),intent(out),optional :: dhdd  ! change in height per diameter [m/cm]
    
    associate( dbh_maxh => EDPftvarcon_inst%allom_dbh_maxheight(ipft), &
               p1       => EDPftvarcon_inst%allom_d2h1(ipft), &
               p2       => EDPftvarcon_inst%allom_d2h2(ipft), &
               p3       => EDPftvarcon_inst%allom_d2h3(ipft), &
               allom_hmode => EDPftvarcon_inst%allom_hmode(ipft))
      
      select case(int(allom_hmode))
      case (1)   ! "obrien"
         call d2h_obrien(d,p1,p2,dbh_maxh,h,dhdd)
      case (2)   ! "poorter06"
         call d2h_poorter2006(d,p1,p2,p3,dbh_maxh,h,dhdd)
      case (3)   ! "2parameter power function h=a*d^b "
         call d2h_2pwr(d,p1,p2,dbh_maxh,h,dhdd)
      case (4)   ! "chave14")
         call d2h_chave2014(d,p1,p2,p3,dbh_maxh,h,dhdd)
      case (5)   ! Martinez-Cano
         call d2h_martcano(d,p1,p2,p3,dbh_maxh,h,dhdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined height allometry was specified: ',allom_hmode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      
    end associate
    return
  end subroutine h_allom

  ! ============================================================================
  ! Generic AGB interface
  ! ============================================================================
  
  subroutine bagw_allom(d,ipft,bagw,dbagwdd)


    real(r8),intent(in)    :: d       ! plant diameter [cm]
    integer(i4),intent(in) :: ipft    ! PFT index
    real(r8),intent(out)   :: bagw    ! biomass above ground woody tissues
    real(r8),intent(out),optional :: dbagwdd  ! change in agbw per diameter [kgC/cm]

    real(r8)               :: h       ! height
    real(r8)               :: dhdd    ! change in height wrt d

    associate( p1           => EDPftvarcon_inst%allom_agb1(ipft), &
               p2           => EDPftvarcon_inst%allom_agb2(ipft), &
               p3           => EDPftvarcon_inst%allom_agb3(ipft), &
               p4           => EDPftvarcon_inst%allom_agb4(ipft), &
               wood_density => EDPftvarcon_inst%wood_density(ipft), &
               c2b          => EDPftvarcon_inst%c2b(ipft), &
               agb_frac     => EDPftvarcon_inst%allom_agb_frac(ipft), &
               allom_amode  => EDPftvarcon_inst%allom_amode(ipft))
      
      select case(int(allom_amode))
      case (1) !"salda")
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_salda(d,h,dhdd,p1,p2,p3,p4,wood_density,c2b,agb_frac,bagw,dbagwdd) 
      case (2) !"2par_pwr")
         ! Switch for woodland dbh->drc
         call d2bagw_2pwr(d,p1,p2,c2b,bagw,dbagwdd)
      case (3) !"chave14") 
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_chave2014(d,h,dhdd,p1,p2,wood_density,c2b,bagw,dbagwdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined AGB allometry was specified: ',allom_amode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      
    end associate
    return
  end subroutine bagw_allom

  ! ============================================================================
  ! Generic diameter to maximum leaf biomass interface
  ! ============================================================================
  
  subroutine blmax_allom(d,ipft,blmax,dblmaxdd)

    real(r8),intent(in)    :: d         ! plant diameter [cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]

    associate( dbh_maxh    => EDPftvarcon_inst%allom_dbh_maxheight(ipft), &
               rho         => EDPftvarcon_inst%wood_density(ipft), &
               c2b         => EDPftvarcon_inst%c2b(ipft),          &
               allom_lmode => EDPftvarcon_inst%allom_lmode(ipft),  &
               p1          => EDPftvarcon_inst%allom_d2bl1(ipft),  &
               p2          => EDPftvarcon_inst%allom_d2bl2(ipft),  &
               p3          => EDPftvarcon_inst%allom_d2bl3(ipft))
      
      select case(int(allom_lmode))
      case(1) !"salda")
         call d2blmax_salda(d,p1,p2,p3,rho,dbh_maxh,c2b,blmax,dblmaxdd)
      case(2) !"2par_pwr")
         call d2blmax_2pwr(d,p1,p2,c2b,blmax,dblmaxdd)
      case(3) ! dh2blmax_2pwr
         call dh2blmax_2pwr(d,p1,p2,dbh_maxh,c2b,blmax,dblmaxdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined leaf allometry was specified: ', &
              allom_lmode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
    end associate
    return
  end subroutine blmax_allom

  ! ============================================================================
  ! Generic crown area allometry wrapper
  ! ============================================================================
  
  subroutine carea_allom(d,nplant,site_spread,ipft,c_area)
     
     real(r8),intent(in)    :: d           ! plant diameter [cm]
     real(r8),intent(in)    :: site_spread ! site level spread factor (crowdedness)
     real(r8),intent(in)    :: nplant      ! number of plants [1/ha]
     integer(i4),intent(in) :: ipft        ! PFT index
     real(r8),intent(out)   :: c_area       ! crown area per plant (m2)

     real(r8)               :: d_eff     ! Effective diameter (cm)
     
     associate( dbh_maxh    => EDPftvarcon_inst%allom_dbh_maxheight(ipft), &
                allom_lmode => EDPftvarcon_inst%allom_lmode(ipft),  &
                d2bl_p2     => EDPftvarcon_inst%allom_d2bl2(ipft),  &
                d2bl_ediff  => EDPftvarcon_inst%allom_blca_expnt_diff(ipft), &
                d2ca_min    => EDPftvarcon_inst%allom_d2ca_coefficient_min(ipft), &
                d2ca_max    => EDPftvarcon_inst%allom_d2ca_coefficient_max(ipft))
       
       select case(int(allom_lmode))
       case(1,3) ! "salda" and "height capped generic two power"
          d_eff = min(d,dbh_maxh)
          call carea_2pwr(d_eff,site_spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max,c_area)
       case(2)   ! "2par_pwr")
          call carea_2pwr(d,site_spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max,c_area)
       case DEFAULT
         write(fates_log(),*) 'An undefined leaf allometry was specified: ', &
               allom_lmode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select

      c_area = c_area * nplant


    end associate
    return
  end subroutine carea_allom

  ! =====================================================================================
        
  subroutine bleaf(d,ipft,canopy_trim,bl,dbldd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target bleaf
    ! based on trimming. Because trimming
    ! is not allometry and rather an emergent property,
    ! this routine is not name-spaced with allom_
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d             ! plant diameter [cm]
    integer(i4),intent(in) :: ipft          ! PFT index
    real(r8),intent(in)    :: canopy_trim   ! trimming function
    real(r8),intent(out)   :: bl            ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dbldd  ! change leaf bio per diameter [kgC/cm]
    
    real(r8) :: blmax
    real(r8) :: dblmaxdd
    
    call blmax_allom(d,ipft,blmax,dblmaxdd)
    
    ! -------------------------------------------------------------------------
    ! Adjust for canopies that have become so deep that their bottom layer is 
    ! not producing any carbon... 
    ! nb this will change the allometry and the effects of this remain untested
    ! RF. April 2014  
    ! -------------------------------------------------------------------------
    
    bl = blmax * canopy_trim
    
    if(present(dbldd))then
       dbldd = dblmaxdd * canopy_trim
    end if
    
    return
  end subroutine bleaf
  
  ! =====================================================================================

  subroutine storage_fraction_of_target(b_leaf, bstore, frac)

    !--------------------------------------------------------------------------------
    ! returns the storage pool as a fraction of its target (only if it is below its target)
    ! used in both the carbon starvation mortlaity scheme as wella s the optional respiration throttling logic
    !--------------------------------------------------------------------------------

    real(r8),intent(in)    :: b_leaf
    real(r8),intent(in)    :: bstore
    real(r8),intent(out)   :: frac

    if( b_leaf > 0._r8 .and. bstore <= b_leaf )then
       frac = bstore/ b_leaf
    else
       frac = 1._r8
    endif

  end subroutine storage_fraction_of_target

  ! =====================================================================================
  
  real(r8) function tree_lai( bl, status_coh, pft, c_area, n )

    ! ============================================================================
    !  LAI of individual trees is a function of the total leaf area and the total canopy area.   
    ! ============================================================================

    real(r8), intent(in) :: bl            ! plant leaf biomass [kg]     
    integer, intent(in)  :: status_coh    ! growth status of plant  (2 = leaves on , 1 = leaves off)
    integer, intent(in)  :: pft
    real(r8), intent(in) :: c_area        ! areal extent of canopy (m2)
    real(r8), intent(in) :: n             ! number of individuals in cohort per 'area' (10000m2 default)

    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC

    if( bl  <  0._r8 .or. pft  ==  0 ) then
       write(fates_log(),*) 'problem in treelai',bl,pft
    endif

    slat = g_per_kg * EDPftvarcon_inst%slatop(pft) ! m2/g to m2/kg
    leafc_per_unitarea = bl/(c_area/n) !KgC/m2
    if(leafc_per_unitarea > 0.0_r8)then
       tree_lai = leafc_per_unitarea * slat  !kg/m2 * m2/kg = unitless LAI 
    else
       tree_lai = 0.0_r8
    endif


    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(tree_lai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much lai' , tree_lai , pft , nlevleaf * dinc_ed
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    return

  end function tree_lai

  ! ============================================================================

  real(r8) function tree_sai( dbh, pft, canopy_trim, c_area, n )

    ! ============================================================================
    !  SAI of individual trees is a function of the target leaf biomass
    ! ============================================================================

    real(r8),intent(in)  :: dbh
    integer, intent(in)  :: pft
    real(r8),intent(in)  :: canopy_trim
    real(r8), intent(in) :: c_area        ! areal extent of canopy (m2)
    real(r8), intent(in) :: n             ! number of individuals in cohort per 'area' (10000m2 default)

    real(r8) :: leafc_per_unitarea ! KgC of target leaf per m2 area of ground.
    real(r8) :: sai_scaler     
    real(r8) :: b_leaf

    sai_scaler = g_per_kg * EDPftvarcon_inst%allom_sai_scaler(pft)  ! m2/g to m2/kg

    call bleaf(dbh,pft,canopy_trim,b_leaf)

    leafc_per_unitarea = b_leaf/(c_area/n) !KgC/m2

    tree_sai = leafc_per_unitarea * sai_scaler !kg/m2 * m2/kg = unitless SAI 

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(tree_sai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much sai' , tree_sai , pft , nlevleaf * dinc_ed
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    return

  end function tree_sai
  
  ! ============================================================================
  ! Generic sapwood biomass interface
  ! ============================================================================

  subroutine bsap_allom(d,ipft,canopy_trim,bsap,dbsapdd)
    
    real(r8),intent(in)           :: d         ! plant diameter [cm]
    integer(i4),intent(in)        :: ipft      ! PFT index
    real(r8),intent(in)           :: canopy_trim
    real(r8),intent(out)          :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per d [kgC/cm]

    real(r8) :: h         ! Plant height [m]
    real(r8) :: dhdd
    real(r8) :: bl
    real(r8) :: dbldd
    real(r8) :: bbgw
    real(r8) :: dbbgwdd
    real(r8) :: bagw
    real(r8) :: dbagwdd
    real(r8) :: bsap_cap  ! cap sapwood so that it is no larger
                          ! than some specified proportion of woody biomass
                          ! should not trip, and only in small plants

    ! Constrain sapwood so that its above ground portion be no larger than 
    ! X% of total woody/fibrous (ie non leaf/fineroot) tissues
    real(r8),parameter :: max_frac = 0.95_r8 

    select case(int(EDPftvarcon_inst%allom_smode(ipft)))
       ! ---------------------------------------------------------------------
       ! Currently both sapwood area proportionality methods use the same
       ! machinery.  The only differences are related to the parameter
       ! checking at the beginning.  For constant proportionality, the slope
       ! of the la:sa to diameter line is zero.
       ! ---------------------------------------------------------------------
    case(1,2) !"constant","dlinear") 

       call h_allom(d,ipft,h,dhdd)
       call bleaf(d,ipft,canopy_trim,bl,dbldd)
       call bsap_dlinear(d,h,dhdd,bl,dbldd,ipft,bsap,dbsapdd)

       ! Perform a capping/check on total woody biomass
       call bagw_allom(d,ipft,bagw,dbagwdd)
       call bbgw_allom(d,ipft,bbgw,dbbgwdd)
       
       ! Force sapwood to be less than a maximum fraction of total biomass
       ! (this comes into play typically in very small plants)
       bsap_cap = max_frac*(bagw+bbgw)
       bsap     = min( bsap_cap,bsap)

       if(present(dbsapdd))then
          if ( bsap  >= bsap_cap ) then
             dbsapdd = max_frac*(dbagwdd+dbbgwdd)
          end if
       end if

    case(9) ! deprecated (9)

       call h_allom(d,ipft,h,dhdd)
       call bleaf(d,ipft,canopy_trim,bl,dbldd)
       call bsap_deprecated(d,h,dhdd,bl,dbldd,ipft,bsap,dbsapdd)

    case DEFAULT
       write(fates_log(),*) 'An undefined sapwood allometry was specified: ', &
            EDPftvarcon_inst%allom_smode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    return
  end subroutine bsap_allom
  
  ! ============================================================================
  ! Generic below ground woody biomass (structure and live/conducting tissues)
  ! (this pool, (like agb and leaves) is assumed to contain all belowground
  ! non-fineroot biomass.
  ! ============================================================================

  subroutine bbgw_allom(d,ipft,bbgw,dbbgwdd)

    real(r8),intent(in)           :: d          ! plant diameter [cm]
    integer(i4),intent(in)        :: ipft       ! PFT index
    real(r8),intent(out)          :: bbgw       ! below ground woody biomass [kgC]
    real(r8),intent(out),optional :: dbbgwdd    ! change bbgw  per diam [kgC/cm]
    
    real(r8)    :: bagw       ! above ground biomass [kgC]
    real(r8)    :: dbagwdd    ! change in agb per diameter [kgC/cm]
    
    select case(int(EDPftvarcon_inst%allom_cmode(ipft)))
    case(1) !"constant")
       call bagw_allom(d,ipft,bagw,dbagwdd)
       call bbgw_const(d,bagw,dbagwdd,ipft,bbgw,dbbgwdd)
    case DEFAULT
       write(fates_log(),*) 'An undefined coarse root allometry was specified: ', &
             EDPftvarcon_inst%allom_cmode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    return
  end subroutine bbgw_allom
 
  ! ============================================================================
  ! Fine root biomass allometry wrapper
  ! ============================================================================
  
  subroutine bfineroot(d,ipft,canopy_trim,bfr,dbfrdd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target fineroot biomass
    ! based on functions that may or may not have prognostic properties. 
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d              ! plant diameter [cm]
    integer(i4),intent(in) :: ipft           ! PFT index
    real(r8),intent(in)    :: canopy_trim    ! trimming function
    real(r8),intent(out)   :: bfr            ! fine root biomass [kgC]
    real(r8),intent(out),optional :: dbfrdd  ! change leaf bio per diameter [kgC/cm]
    
    real(r8) :: blmax      ! maximum leaf biomss per allometry
    real(r8) :: dblmaxdd
    real(r8) :: bfrmax
    real(r8) :: dbfrmaxdd
    real(r8) :: slascaler
    
    select case(int(EDPftvarcon_inst%allom_fmode(ipft)))
    case(1) ! "constant proportionality with TRIMMED target bleaf"
       
       call blmax_allom(d,ipft,blmax,dblmaxdd)
       call bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
       bfr    = bfrmax * canopy_trim
       if(present(dbfrdd))then
          dbfrdd = dbfrmaxdd * canopy_trim
       end if
    case(2) ! "constant proportionality with UNTRIMMED target bleaf"
       
       call blmax_allom(d,ipft,blmax,dblmaxdd)
       call bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
       bfr    = bfrmax
       if(present(dbfrdd))then
          dbfrdd = dbfrmaxdd
       end if

    case DEFAULT 
       write(fates_log(),*) 'An undefined fine root allometry was specified: ', &
            EDPftvarcon_inst%allom_fmode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    
    return
  end subroutine bfineroot


  ! ============================================================================
  ! Storage biomass interface
  ! ============================================================================
  
  subroutine bstore_allom(d,ipft,canopy_trim,bstore,dbstoredd)

     real(r8),intent(in)           :: d            ! plant diameter [cm]
     integer(i4),intent(in)        :: ipft         ! PFT index
     real(r8),intent(in)           :: canopy_trim  ! Crown trimming function [0-1]
     real(r8),intent(out)          :: bstore       ! allometric target storage [kgC]
     real(r8),intent(out),optional :: dbstoredd    ! change storage per cm [kgC/cm]
     
     real(r8) :: bl          ! Allometric target leaf biomass
     real(r8) :: dbldd       ! Allometric target change in leaf biomass per cm
    
     
     ! TODO: allom_stmode needs to be added to the parameter file
     
     associate( allom_stmode => EDPftvarcon_inst%allom_stmode(ipft), &
                cushion      => EDPftvarcon_inst%cushion(ipft) )

       select case(int(allom_stmode))
       case(1) ! Storage is constant proportionality of trimmed maximum leaf
          ! biomass (ie cushion * bleaf)
          
          call bleaf(d,ipft,canopy_trim,bl,dbldd)
          call bstore_blcushion(d,bl,dbldd,cushion,ipft,bstore,dbstoredd)
          
       case DEFAULT 
          write(fates_log(),*) 'An undefined fine storage allometry was specified: ', &
                allom_stmode
          write(fates_log(),*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       
     end associate
     return
  end subroutine bstore_allom
  


  ! ============================================================================
  ! Dead biomass interface
  ! ============================================================================
  
  subroutine bdead_allom(bagw,bbgw,bsap,ipft,bdead,dbagwdd,dbbgwdd,dbsapdd,dbdeaddd)

     real(r8),intent(in)  :: bagw      ! biomass above ground wood (agb) [kgC]
     real(r8),intent(in)  :: bbgw       ! biomass below ground (bgb) [kgC]
     real(r8),intent(in)  :: bsap      ! sapwood biomass [kgC]
     integer(i4),intent(in) :: ipft    ! PFT index
     real(r8),intent(out) :: bdead     ! dead biomass (heartw/struct) [kgC]
     
     real(r8),intent(in),optional  :: dbagwdd    ! change in agb per d [kgC/cm]
     real(r8),intent(in),optional  :: dbbgwdd    ! change in bbgw per d [kgC/cm]
     real(r8),intent(in),optional  :: dbsapdd   ! change in bsap per d [kgC/cm]
     real(r8),intent(out),optional :: dbdeaddd  ! change in bdead per d [kgC/cm]
     
     ! bdead is diagnosed as the mass balance from all other pools
     ! and therefore, no options are necessary
     !
     ! Assumption: We assume that the leaf biomass component of AGB is negligable
     ! and do not assume it is part of the AGB measurement, nor are fineroots part of the
     ! bbgw. Therefore, it is not removed from AGB and BBGW in the calculation of dead mass.

    
    associate( agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft))

      select case(int(EDPftvarcon_inst%allom_amode(ipft)))
      case(1) ! Saldariagga mass allometry originally calculated bdead directly.
              ! we assume proportionality between bdead and bagw
       
         bdead = bagw/agb_fraction 
         if(present(dbagwdd) .and. present(dbdeaddd))then
            dbdeaddd = dbagwdd/agb_fraction
         end if
         
      case(2,3)
         
         bdead = bagw + bbgw - bsap
         if(present(dbagwdd) .and. present(dbbgwdd) .and. &
            present(dbdeaddd) .and. present(dbsapdd) )then
            dbdeaddd = dbagwdd+dbbgwdd-dbsapdd
         end if
         
      case DEFAULT
         
         write(fates_log(),*) 'An undefined AGB allometry was specified: ',&
                              EDPftvarcon_inst%allom_amode(ipft)
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
       
      end select
    end associate
    return
  end subroutine bdead_allom

  ! ============================================================================
  ! Specific bfrmax relationships
  ! ============================================================================
  
  subroutine bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: blmax     ! max leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bfrmax    ! max fine-root root biomass [kgC]
    real(r8),intent(out),optional :: dbfrmaxdd ! change frmax bio per diam [kgC/cm]
    
    associate( l2fr => EDPftvarcon_inst%allom_l2fr(ipft) )
      
      bfrmax = blmax*l2fr
      
      ! dbfr/dd = dbfrmax/dblmax * dblmax/dd
      if(present(dbfrmaxdd))then
         dbfrmaxdd = dblmaxdd*l2fr
      end if
      
    end associate
    return
  end subroutine bfrmax_const

  ! ============================================================================
  ! Specific bbgw relationships
  ! ============================================================================
  
  subroutine bbgw_const(d,bagw,dbagwdd,ipft,bbgw,dbbgwdd)

    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: bagw       ! above ground biomass [kg]
    real(r8),intent(in)    :: dbagwdd    ! change in agb per diameter [kg/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bbgw       ! coarse root biomass [kg]
    real(r8),intent(out),optional :: dbbgwdd    ! change croot bio per diam [kg/cm]

    associate( agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft) )
      
      bbgw = (1.0_r8/agb_fraction-1.0_r8)*bagw 
      
      ! Derivative
      ! dbbgw/dd = dbbgw/dbagw * dbagw/dd
      if(present(dbbgwdd))then
         dbbgwdd = (1.0_r8/agb_fraction-1.0_r8)*dbagwdd
      end if
      
    end associate
    return
 end subroutine bbgw_const

 ! ============================================================================
 ! Specific d2bsap relationships
 ! ============================================================================
  
 subroutine bsap_deprecated(d,h,dhdd,bleaf,dbleafdd,ipft,bsap,dbsapdd)
    

    
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    real(r8),intent(in)    :: bleaf     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dbleafdd  ! change in blmax per diam [kgC/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per diameter [kgC/cm]
    
    real(r8)               :: latosa    ! applied leaf area to sap area 
                                        ! may or may not contain diameter correction
    real(r8)               :: hbl2bsap  ! sapwood biomass per lineal height and kg of leaf
    
    
    associate ( latosa_int => EDPftvarcon_inst%allom_latosa_int(ipft), &
                latosa_slp => EDPftvarcon_inst%allom_latosa_slp(ipft), &
                sla        => EDPftvarcon_inst%slatop(ipft), &
                wood_density => EDPftvarcon_inst%wood_density(ipft), &
                c2b          => EDPftvarcon_inst%c2b(ipft), & 
                agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft) )

      ! notes:
      ! latosa_int units of [/m]
      ! (1/latosa)* slatop*    gtokg    *   cm2tom2     / c2b   * mg2kg  * dens
      ! density (g/cm3 == Mg/m3 b/c  1e6 = 100^3)
      ! [cm2/m2] * [m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3] = [/m]
      !          0.012 * 1000 * (1/10000) / 2 * 1000 * 0.7

      bsap = bleaf * latosa_int * h
      
      if(present(dbsapdd))then
         dbsapdd = latosa_int*(h*dbleafdd + bleaf*dhdd)
      end if

    end associate
    return
  end subroutine bsap_deprecated

  ! ========================================================================

  subroutine bsap_dlinear(d,h,dhdd,bleaf,dbleafdd,ipft,bsap,dbsapdd)
    
    ! -------------------------------------------------------------------------
    ! Calculate sapwood biomass based on leaf area to sapwood area
    ! proportionality.  In this function, the leaftosapwood area is a function
    ! of plant size, see Calvo-Alvarado and Bradley Christoferson
    ! In this case: parameter latosa (from constant proportionality)
    !   is the intercept of the diameter function.
    ! 
    ! Important note: this is above and below-ground sapwood
    !
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    real(r8),intent(in)    :: bleaf     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dbleafdd  ! change in blmax per diam [kgC/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per diameter [kgC/cm]
    
    real(r8)               :: latosa    ! applied leaf area to sap area 
                                        ! may or may not contain diameter correction
    real(r8)               :: hbl2bsap  ! sapwood biomass per lineal height and kg of leaf
    
    
    associate ( latosa_int => EDPftvarcon_inst%allom_latosa_int(ipft), &
                latosa_slp => EDPftvarcon_inst%allom_latosa_slp(ipft), &
                sla        => EDPftvarcon_inst%slatop(ipft), &
                wood_density => EDPftvarcon_inst%wood_density(ipft), &
                c2b          => EDPftvarcon_inst%c2b(ipft), & 
                agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft) )

      ! ------------------------------------------------------------------------
      ! Calculate sapwood biomass per linear height and kgC of leaf [m-1]
      ! Units: 
      ! latosa * slatop*    gtokg    *   cm2tom2     / c2b     * mg2kg  * dens
      ! [cm2/m2]*[m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3]
      !        ->[cm2/gC]
      !                  ->[cm2/kgC]
      !                                ->[m2/kgC]
      !                                             ->[m2/kg]
      !                                                       ->[m2/Mg]
      !                                                                  ->[/m]
      ! ------------------------------------------------------------------------

      latosa   = latosa_int + d*latosa_slp
      hbl2bsap = latosa*sla*g_per_kg*wood_density*kg_per_Megag/(c2b*cm2_per_m2 )
      
      bsap =  hbl2bsap * (h/agb_fraction) * bleaf

      ! Derivative
      ! dbldmaxdd is deriv of blmax wrt dbh (use directives to check oop)
      ! dhdd is deriv of height wrt dbh (use directives to check oop)
      if(present(dbsapdd))then
         dbsapdd = hbl2bsap*(h*dbleafdd + bleaf*dhdd)/agb_fraction
      end if
      
    end associate
    return
  end subroutine bsap_dlinear

  ! ============================================================================
  ! Specific storage relationships
  ! ============================================================================
  
  subroutine bstore_blcushion(d,bl,dbldd,cushion,ipft,bstore,dbstoredd)
     
     ! This discracefully simple subroutine calculates allometric target
     ! storage biomass based on a constant-specified ratio (cushion)
     ! of storage to target allometricc leaf biomass

     real(r8),intent(in)    :: d                  ! plant diameter [cm]
     real(r8),intent(in)    :: bl                 ! plant leaf biomass [kgC]
     real(r8),intent(in)    :: dbldd              ! change in blmax per diam [kgC/cm]
     real(r8),intent(in)    :: cushion            ! simple constant ration bstore/bleaf
     integer(i4),intent(in) :: ipft               ! PFT index
     real(r8),intent(out)   :: bstore             ! plant leaf biomass [kgC]
     real(r8),intent(out),optional :: dbstoredd   ! change leaf bio per diameter [kgC/cm]
     
     
     bstore = bl * cushion
     
     if(present(dbstoredd)) then
        dbstoredd = dbldd * cushion
     end if

     return
  end subroutine bstore_blcushion


  ! ============================================================================
  ! Specific d2blmax relationships
  ! ============================================================================
  
  subroutine d2blmax_salda(d,p1,p2,p3,rho,dbh_maxh,c2b,blmax,dblmaxdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: p1
    real(r8),intent(in)    :: p2
    real(r8),intent(in)    :: p3
    real(r8),intent(in)    :: rho       ! plant's wood specific gravity
    real(r8),intent(in)    :: dbh_maxh  ! dbh at which max height occurs
    real(r8),intent(in)    :: c2b       ! c to biomass multiplier (~2.0)
    
    real(r8),intent(out)   :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out),optional   :: dblmaxdd  ! change leaf bio per diam [kgC/cm]
    
    if( d<dbh_maxh) then
       blmax = p1 * d**p2 * rho**p3
    else
       blmax = p1 * dbh_maxh**p2 * rho**p3
    end if
    
    if(present(dblmaxdd))then
       if( d<dbh_maxh) then
          dblmaxdd = p1*p2 * d**(p2-1.0_r8) * rho**p3
       else
          dblmaxdd = 0.0
       end if
    end if
    
    return
  end subroutine d2blmax_salda
  
  ! ===========================================================================
  
  subroutine d2blmax_2pwr(d,p1,p2,c2b,blmax,dblmaxdd)
    
    ! ======================================================================
    ! This is a power function for leaf biomass from plant diameter.
    !
    ! log(bl) = a2 + b2*log(h)
    ! bl = exp(a2) * h**b2
    ! ======================================================================
    
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    real(r8),intent(in)  :: p1        ! parameter 1
    real(r8),intent(in)  :: p2        ! parameter 2
    real(r8),intent(in)  :: c2b       ! carbon to biomass multiplier
    
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    
    blmax    = p1*d**p2 / c2b
    
    if(present(dblmaxdd))then
       dblmaxdd = p1*p2 *d**(p2-1.0_r8) / c2b
    end if
    
    return
  end subroutine d2blmax_2pwr

  ! ===========================================================================
  
  subroutine dh2blmax_2pwr(d,p1,p2,dbh_maxh,c2b,blmax,dblmaxdd)
    
    ! -------------------------------------------------------------------------
    ! This formulation is very similar to d2blmax_2pwr
    ! The difference is that for very large trees that have reached peak
    ! height, we cap leaf biomass.
    ! --------------------------------------------------------------------------
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    real(r8),intent(in)  :: p1        ! parameter 1
    real(r8),intent(in)  :: p2        ! parameter 2
    real(r8),intent(in)  :: c2b       ! carbon 2 biomass multiplier
    real(r8),intent(in)  :: dbh_maxh  ! dbh at maximum height
    
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    
    blmax    = p1*min(d,dbh_maxh)**p2/c2b
    
    ! If this plant has reached its height cap, then it is not
    ! adding leaf mass.  In this case, dhdd = 0
    if(present(dblmaxdd))then
       if(d>=dbh_maxh)then
          dblmaxdd = 0.0_r8
       else
          dblmaxdd = p1*p2*d**(p2-1.0_r8) / c2b
       end if
    end if
    
    return
  end subroutine dh2blmax_2pwr
  
  ! =========================================================================
  ! Diameter to height (D2H) functions
  ! =========================================================================
  
  subroutine d2h_chave2014(d,p1,p2,p3,dbh_maxh,h,dhdd)
    
    ! "d2h_chave2014"
    ! "d to height via Chave et al. 2014"
    
    ! This function calculates tree height based on tree diameter and the
    ! environmental stress factor "E", as per Chave et al. 2015 GCB
    ! As opposed to previous allometric models in ED, in this formulation
    ! we do not impose a hard cap on tree height.  But, maximum_height
    ! is an important parameter, but instead of imposing a hard limit, in
    ! the new methodology, it will be used to trigger a change in carbon
    ! balance accounting.  Such that a tree that hits its maximum height will
    ! begin to route available NPP into seed and defense respiration.
    !
    ! The stress function is based on the geographic location of the site.  If
    ! a user decides to use Chave2015 allometry, the E factor will be read in
    ! from a global gridded dataset and assigned for each ED patch (note it
    ! will be the same for each ED patch, but this distinction will help in
    ! porting ED into different models (patches are pure ED).  It
    ! assumes that the site is within the pan-tropics, and is a linear function
    ! of climatic water deficit, temperature seasonality and precipitation
    ! seasonality.  See equation 6b of Chave et al.
    ! The relevant equation for height in this function is 6a of the same
    ! manuscript, and is intended to pair with diameter to relate with
    ! structural biomass as per equation 7 (in which H is implicit).
    !
    ! Chave et al. Improved allometric models to estimate the abovegroud
    ! biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
    !
    ! =========================================================================
    
    real(r8),intent(in)  :: d        ! plant diameter [cm]
    real(r8),intent(in)  :: p1       ! parameter a 
    real(r8),intent(in)  :: p2       ! parameter b
    real(r8),intent(in)  :: p3       ! parameter c
    real(r8),intent(in)  :: dbh_maxh ! dbh at maximum height [cm]
    
    real(r8),intent(out) :: h     ! plant height [m]
    real(r8),intent(out),optional :: dhdd  ! change in height per diameter [m/cm]
    real(r8) :: p1e
    
    p1e = p1   ! -eclim (assumed that p1 already has eclim removed)
    if(d>=dbh_maxh) then
       h    = exp( p1e + p2*log(dbh_maxh) + p3*log(dbh_maxh)**2.0 )
    else
       h    = exp( p1e + p2*log(d) + p3*log(d)**2.0 )
    end if
    
    if(present(dhdd))then
       if(d>=dbh_maxh ) then
          dhdd = 0.0_r8
       else
          dhdd = exp(p1e) * ( 2.0_r8*p3*d**(p2-1.0_r8+p3*log(d))*log(d) + &
               p2*d**(p2-1.0_r8+p3*log(d)) )
       end if
    end if
    return
  end subroutine d2h_chave2014

  ! ===========================================================================
  
  subroutine d2h_poorter2006(d,p1,p2,p3,dbh_maxh,h,dhdd)
    
    ! "d2h_poorter2006"
    ! "d to height via Poorter et al. 2006, these routines use natively
    !  asymtotic functions"
    !
    ! Poorter et al calculated height diameter allometries over a variety of
    ! species in Bolivia, including those that could be classified in guilds
    ! that were Partial-shade-tolerants, long-lived pioneers, shade-tolerants
    ! and short-lived pioneers.  There results between height and diameter
    ! found that small stature trees had less of a tendency to asymotote in
    ! height and showed more linear relationships, and the largest stature
    ! trees tended to show non-linear relationships that asymtote.
    !
    ! h = h_max*(1-exp(-a*d**b))
    !
    ! Poorter L, L Bongers and F Bongers.  Architecture of 54 moist-forest tree
    ! species: traits, trade-offs, and functional groups.  Ecology 87(5), 2006.
    !
    ! =========================================================================

    real(r8),intent(in)  :: d     ! plant diameter [cm] 
    real(r8),intent(in)     :: p1       ! parameter a = h_max
    real(r8),intent(in)     :: p2       ! parameter b
    real(r8),intent(in)     :: p3       ! parameter c
    real(r8),intent(in)     :: dbh_maxh ! dbh at maximum height
    real(r8),intent(out) :: h     ! plant height [m]
    real(r8),intent(out),optional :: dhdd  ! change in height per diameter [m/cm]
    
    h = p1*(1.0_r8 - exp(p2*min(d,dbh_maxh)**p3))
    
    !h = h_max - h_max (exp(a*d**b))
    !f(x) = -h_max*exp(g(x))
    !g(x) = a*d**b
    !d/dx f(g(x) = f'(g(x))*g'(x) = -a1*exp(a2*d**a3) * a3*a2*d**(a3-1)

    if(present(dhdd))then
       if( d>=dbh_maxh ) then
          dhdd = 0.0_r8
       else
          dhdd = -p1*exp(p2*d**p3) * p3*p2*d**(p3-1.0_r8)
       end if
    end if

    return
  end subroutine d2h_poorter2006

  ! ===========================================================================
  
  subroutine d2h_2pwr(d,p1,p2,dbh_maxh,h,dhdd)
    
    ! =========================================================================
    ! "d2h_2pwr"
    ! "d to height via 2 parameter power function"
    ! where height h is related to diameter by a linear relationship on the log
    ! transform where log(a) is the intercept and b is the slope parameter.
    !
    ! log(h) = log(a) + b*log(d)
    ! h      = exp(log(a)) * exp(log(d))**b
    ! h      = a*d**b
    !
    ! This functional form is used often in temperate allometries
    ! Therefore, no base reference is cited.  Although, the reader is pointed
    ! towards Dietze et al. 2008, King 1991, Ducey 2012 and many others for
    ! reasonable parameters.  Note that this subroutine is intended only for
    ! trees well below their maximum height, ie initialization.
    !
    ! =========================================================================
    ! From King et al. 1990 at BCI for saplings
    ! log(d) = a + b*log(h)
    ! d = exp(a) * h**b
    ! h = (1/exp(a)) * d**(1/b)
    ! h = p1*d**p2  where p1 = 1/exp(a) = 1.07293  p2 = 1/b = 1.4925
    ! d = (h/p1)**(1/p2)
    ! For T. tuberculata (canopy tree) a = -0.0704, b = 0.67
    ! =========================================================================
    
    ! args
    ! =========================================================================
    ! d: diameter at breast height
    ! p1: the intercept parameter 
    !                       (however exponential of the fitted log trans)
    ! p2: the slope parameter
    ! return:
    ! h: total tree height [m]
    ! =========================================================================

    real(r8),intent(in)     :: d        ! plant diameter [cm]
    real(r8),intent(in)     :: p1       ! parameter a
    real(r8),intent(in)     :: p2       ! parameter b
    real(r8),intent(in)     :: dbh_maxh ! dbh where max height occurs [cm]
    real(r8),intent(out)    :: h        ! plant height [m]
    real(r8),intent(out),optional    :: dhdd     ! change in height per diameter [m/cm]
    
    h    = p1*min(d,dbh_maxh)**p2
    
    if(present(dhdd))then
       if( d>=dbh_maxh ) then
          dhdd = 0.0_r8
       else
          dhdd = (p2*p1)*d**(p2-1.0_r8)
       end if
    end if
    
    return
  end subroutine d2h_2pwr
  
  ! ============================================================================
  
  subroutine d2h_obrien(d,p1,p2,dbh_maxh,h,dhdd)
    
    real(r8),intent(in)    :: d        ! plant diameter [cm]
    real(r8),intent(in)    :: p1       ! parameter a
    real(r8),intent(in)    :: p2       ! parameter b
    real(r8),intent(in)    :: dbh_maxh ! dbh where max height occurs [cm]
    real(r8),intent(out)   :: h        ! plant height [m]
    real(r8),intent(out),optional   :: dhdd     ! change in height per diameter [m/cm]
    
    !p1 = 0.64
    !p2 = 0.37
    h    = 10.0_r8**(log10(min(d,dbh_maxh))*p1+p2)
    
    if(present(dhdd))then
       if(d>=dbh_maxh ) then
          dhdd = 0.0_r8
       else
          dhdd = p1*10.0_r8**p2*d**(p1-1.0_r8)
       end if
    end if
    
    
    return
  end subroutine d2h_obrien
  
  ! ===========================================================================
  
  subroutine d2h_martcano(d,p1,p2,p3,dbh_maxh,h,dhdd)
    
    ! =========================================================================
    ! "d2h_martcano"
    ! "d to height via 3 parameter Michaelis-Menten following work at BCI
    ! by Martinez-Cano et al. 2016
    ! 
    ! h = (a*d**b)/(c+d**b)
    !
    ! h' = [(a*d**b)'(c+d**b) - (c+d**b)'(a*d**b)]/(c+d**b)**2
    ! dhdd = h' = [(ba*d**(b-1))(c+d**b) - (b*d**(b-1))(a*d**b)]/(c+d**b)**2
    !
    ! args
    ! =========================================================================
    ! d: diameter at breast height
    ! h: total tree height [m]
    ! =========================================================================
    
    real(r8),intent(in)  :: d     ! plant diameter [cm]   
    real(r8),intent(in)  :: p1       ! parameter a
    real(r8),intent(in)  :: p2       ! parameter b 
    real(r8),intent(in)  :: p3       ! parameter c
    real(r8),intent(in)  :: dbh_maxh ! diameter at maximum height
    real(r8),intent(out) :: h     ! plant height [m]
    real(r8),intent(out),optional :: dhdd  ! change in height per diameter [m/cm]
    
    h = (p1*min(d,dbh_maxh)**p2)/(p3+min(d,dbh_maxh)**p2)
    
    if(present(dhdd))then
       if(d>=dbh_maxh ) then
          dhdd = 0.0
       else
          dhdd = ((p2*p1*d**(p2-1._r8))*(p3+d**p2) - &
               (p2*d**(p2-1._r8))*(p1*d**p2)        )/ &
               (p3+d**p2)**2._r8
       end if
    end if
    
    return
  end subroutine d2h_martcano
  
  
  ! =========================================================================
  ! Diameter to (2) above-ground biomass
  ! =========================================================================
  
  subroutine dh2bagw_chave2014(d,h,dhdd,p1,p2,wood_density,c2b,bagw,dbagwdd)
    
    ! =========================================================================
    ! This function calculates tree structural biomass from tree diameter,
    ! height and wood density.
    !
    ! Chave et al. Improved allometric models to estimate the abovegroud
    ! biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
    !
    ! Input arguments:
    ! d: Diameter at breast height [cm]
    ! rho:  wood specific gravity (dry wood mass per green volume)
    ! height: total tree height [m]
    ! a1: structural biomass allometry parameter 1 (intercept)
    ! a2: structural biomass allometry parameter 2 (slope)
    ! Output:
    ! bagw:   Total above ground biomass [kgC]
    !
    ! =========================================================================
    
    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: h       ! plant height [m]
    real(r8),intent(in)  :: dhdd    ! change in height wrt diameter
    real(r8),intent(in)  :: p1  ! allometry parameter 1
    real(r8),intent(in)  :: p2  ! allometry parameter 2
    real(r8),intent(in)  :: wood_density
    real(r8),intent(in)  :: c2b
    real(r8),intent(out) :: bagw     ! plant height [m]
    real(r8),intent(out),optional :: dbagwdd  ! change in agb per diameter [kgC/cm]
    
    real(r8) :: dbagwdd1,dbagwdd2,dbagwdd3
    
    bagw   = (p1 * (wood_density*d**2.0_r8*h)**p2)/c2b
    
    
    if(present(dbagwdd))then
       ! Need the the derivative of height to diameter to
       ! solve the derivative of agb with height
       dbagwdd1  = (p1*wood_density**p2)/c2b
       dbagwdd2  = p2*d**(2.0_r8*p2)*h**(p2-1.0_r8)*dhdd
       dbagwdd3  = h**p2*2.0_r8*p2*d**(2.0_r8*p2-1.0_r8)
       dbagwdd   = dbagwdd1*(dbagwdd2 + dbagwdd3)
    end if
    
    return
  end subroutine dh2bagw_chave2014
  
  subroutine d2bagw_2pwr(d,p1,p2,c2b,bagw,dbagwdd)
    
    ! =========================================================================
    ! This function calculates tree above ground biomass according to 2
    ! parameter power functions. (slope and intercepts of a log-log
    ! diameter-agb fit:
    !
    ! These relationships are typical for temperate/boreal plants in North
    ! America.  Parameters are available from Chojnacky 2014 and Jenkins 2003
    !
    ! Note that we are using an effective diameter here, as Chojnacky 2014
    ! and Jenkins use diameter at the base of the plant for "woodland" species
    ! The diameters should be converted prior to this routine if drc.
    !
    ! Input arguments:
    ! diam: effective diameter (d or drc) in cm
    ! FOR SPECIES THAT EXPECT DCM, THIS NEEDS TO BE PRE-CALCULATED!!!!
    ! Output:
    ! agb:   Total above ground biomass [kgC]
    !
    ! =========================================================================
    ! Aceraceae, Betulaceae, Fagaceae and Salicaceae comprised nearly
    ! three-fourths of the hardwood species (Table 3)
    !
    ! Fabaceae and Juglandaceae had specific gravities .0.60 and were
    ! combined, as were Hippocastanaceae and Tilaceae with specific gravities
    ! near 0.30. The remaining 9 families, which included mostly species with
    ! specific gravity 0.45–0.55, were initially grouped to construct a general
    ! hardwood taxon for those families having few published biomass equa-
    ! tions however, 3 warranted separation, leaving 6 families for the general
    ! taxon.
    ! bagw = exp(b0 + b1*ln(diameter))/c2b
    ! =========================================================================

    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: p1  ! allometry parameter 1
    real(r8),intent(in)  :: p2  ! allometry parameter 2
    real(r8),intent(in)  :: c2b     ! carbon to biomass multiplier ~2
    real(r8),intent(out) :: bagw     ! plant height [m]
    real(r8),intent(out),optional :: dbagwdd  ! change in agb per diameter [kgC/cm]
    
    bagw    = (p1 * d**p2)/c2b
    if(present(dbagwdd))then
       dbagwdd = (p2*p1*d**(p2-1.0_r8))/c2b
    end if
    
    return
  end subroutine d2bagw_2pwr
  
  
  subroutine dh2bagw_salda(d,h,dhdd,p1,p2,p3,p4, &
       wood_density,c2b,allom_agb_frac,bagw,dbagwdd)
    
    ! --------------------------------------------------------------------
    ! Calculate stem biomass from height(m) dbh(cm) and wood density(g/cm3)
    ! default params using allometry of J.G. Saldarriaga et al 1988 - Rio Negro
    ! Journal of Ecology vol 76 p938-958  
    ! Saldarriaga 1988 provided calculations on total dead biomass
    ! So here, we calculate total dead, and then remove the below-ground
    ! fraction
    ! --------------------------------------------------------------------
    
    real(r8),intent(in) :: d             ! plant diameter [cm]
    real(r8),intent(in) :: h             ! plant height [m]
    real(r8),intent(in) :: dhdd       ! change in height wrt diameter
    real(r8),intent(in) :: p1         !    = 0.06896_r8
    real(r8),intent(in) :: p2         !    = 0.572_r8
    real(r8),intent(in) :: p3         !    = 1.94_r8
    real(r8),intent(in) :: p4         !    = 0.931_r8
    real(r8),intent(in) :: c2b            ! carbon 2 biomass ratio
    real(r8),intent(in) :: wood_density   
    real(r8),intent(in) :: allom_agb_frac
    real(r8),intent(out) :: bagw     ! plant biomass [kgC/indiv]
    real(r8),intent(out),optional :: dbagwdd  ! change in agb per diameter [kgC/cm]
    
    real(r8) :: term1,term2,term3
    
    
    bagw = allom_agb_frac*p1*(h**p2)*(d**p3)*(wood_density**p4)

    ! Add sapwood calculation to this?
    
    ! bagw     = a1 * h**a2 * d**a3 * r**a4
    ! dbagw/dd = a1*r**a4 * d/dd (h**a2 * d**a3)
    ! dbagw/dd = a1*r**a4 * [ d**a3 *d/dd(h**a2) + h**a2*d/dd(d**a3) ]
    ! dbagw/dd = a1*r**a4 * [ d**a3 * a2*h**(a2-1)dh/dd + h**a2*a3*d**(a3-1)]
    
    ! From code
    !   dbagw/dd =  a3 * a1 *(h**a2)*(d**(a3-1))* (r**a4) + a2*a1*(h**(a2-1))*(d**a3)*(r**a4)*dhdd 
    !   dbagw/dd =  a1*r**a4 * [ d**a3 * a2* h**(a2-1)*dhdd +  a3 * (h**a2)*(d**(a3-1)) ]
    
    
    if(present(dbagwdd)) then

       term1 = allom_agb_frac*p1*(wood_density**p4)
       term2 = (h**p2)*p3*d**(p3-1.0_r8)
       term3  = p2*h**(p2-1.0_r8)*(d**p3)*dhdd
       dbagwdd = term1*(term2+term3)
       
    end if

    return
  end subroutine dh2bagw_salda
  
  ! ============================================================================
  ! height to diameter conversions
  ! Note that these conversions will only back-calculate the actual diameter
  ! for plants that have not started to experience height capping or an
  ! asymptote.  In these cases they may be called effective diameter.
  ! ============================================================================

  subroutine h2d_chave2014(h,p1,p2,p3,de,ddedh)

    
    real(r8),intent(in)  :: h       ! plant height [m]
    real(r8),intent(in)  :: p1
    real(r8),intent(in)  :: p2
    real(r8),intent(in)  :: p3
    
    real(r8),intent(out) :: de      ! effective plant diameter [cm]
    real(r8),intent(out),optional :: ddedh   ! effective change in d per height [cm/m]
    
    real(r8) :: p1e, eroot, dbh1,dhpdd
    
    p1e   = p1 !-eclim  (assumed that p1 already has eclim removed)
    eroot = (-p2 + sqrt(p2**2.0_r8 + 4.0_r8*log(h*exp(-p1e))*p3)) & 
         /(2.0_r8*p3)
    
    de = exp(eroot)
    
    if(present(ddedh))then
       ! Invert the derivative at d without asymtote
       dhpdd = exp(p1e)*( p3*2.0_r8*de**(p2-1.0_r8)*log(de)* &
            exp(p3*log(de)**2) + p2*de**(p2-1.0_r8)* &
            exp(p3*log(de)**2.0_r8) )
       
       ddedh = 1.0_r8/dhpdd
    end if
    
    !    term1 = exp(-p2/(2*p3))
    !    term2 = exp(p2**2/(4*p3**2))
    !    term3 = exp(-p1e/p3)
    !    term4 = h**(1/p3-1.0_r8)/(p3)
    !    d   = term1*term2*term3*term4
    return
  end subroutine h2d_chave2014
  
  ! ============================================================================
  
  subroutine h2d_poorter2006(h,p1,p2,p3,d,dddh)
    
    ! -------------------------------------------------------------------------
    ! Note that the height to diameter conversion in poorter is only necessary
    ! when initializing saplings.  In other methods, height to diameter is
    ! useful when defining the dbh at which point to asymptote, as maximum
    ! height is the user set parameter.  This function should not need to set a
    ! dbh_max parameter for instance, but it may end up doing so anyway, even
    ! if it is not used, do to poor filtering.  The poorter et al. d2h and h2d
    ! functions are already asymptotic, and the parameter governing maximum
    ! height is the p1 parameter.  Note as dbh gets very large, the
    ! exponential goes to zero, and the maximum height approaches p1.
    ! However, the asymptote has a much different shape than the logistic, so
    ! therefore in the Poorter et al functions, we do not set p1 == h_max.
    ! That being said, if an h_max that is greater than p1 is passed to this
    ! function, it will return a complex number. During parameter
    ! initialization, a check will be placed that forces:
    ! h_max = p1*0.98
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)  :: h      ! plant height [m]
    real(r8),intent(in)  :: p1
    real(r8),intent(in)  :: p2
    real(r8),intent(in)  :: p3
    
    real(r8),intent(out) :: d      ! plant diameter [cm]
    real(r8),intent(out),optional :: dddh   ! change in d per height [cm/m]
    
    ! -------------------------------------------------------------------------
    ! h = a1*(1 - exp(a2*d**a3))
    ! h = a1 - a1*exp(a2*d**a3)
    ! a1-h = a1*exp(a2*d**a3)
    ! (a1-h)/a1 = exp(a2*d**a3)
    ! log(1-h/a1) = a2*d**a3
    ! [log(1-h/a1)/a2]**(1/a3) = d
    !
    ! derivative dd/dh
    ! dd/dh = [log((a1-h)/a1)/a2]**(1/a3)'
    !       = (1/a3)*[log((a1-h)/a1)/a2]**(1/a3-1)* [(log(a1-h)-log(a1))/a2]'
    !       = (1/a3)*[log((a1-h)/a1)/a2]**(1/a3-1) * (1/(a2*(h-a1))
    ! dd/dh = -((log(1-h/a1)/a2)**(1/a3-1))/(a2*a3*(a1-h))
    ! -------------------------------------------------------------------------
    
    d = (log(1.0_r8-h/p1)/p2)**(1.0_r8/p3)
    
    if(present(dddh))then
       dddh = -((log(1-h/p1)/p2)**(1.0_r8/p3-1.0_r8))/ &
            (p2*p3*(p1-h))
    end if
    
    return
  end subroutine h2d_poorter2006
  
  ! ============================================================================
  
  subroutine h2d_2pwr(h,p1,p2,d,dddh)
    
    
    real(r8),intent(in)  :: h      ! plant height [m]
    real(r8),intent(in)  :: p1     ! parameter 1
    real(r8),intent(in)  :: p2     ! parameter 2
    
    real(r8),intent(out) :: d      ! plant diameter [cm]
    real(r8),intent(out),optional :: dddh   ! change in d per height [cm/m]
    
    !h = a1*d**a2
    d = (h/p1)**(1.0_r8/p2)
    
    !    d = (1/a1)**(1/a2)*h**(1/a2)
    if(present(dddh)) then
       dddh = (1.0_r8/p2)*(1.0_r8/p1)**(1.0_r8/p2) &
            *h**(1.0_r8/p2-1.0_r8)
    end if
    
    return
  end subroutine h2d_2pwr
  
  ! ============================================================================
  
  subroutine h2d_obrien(h,p1,p2,d,dddh)
    
    real(r8),intent(in)  :: h      ! plant height [m]
    real(r8),intent(in)  :: p1
    real(r8),intent(in)  :: p2
    
    real(r8),intent(out)   :: d      ! plant diameter [cm]
    real(r8),intent(out),optional   :: dddh   ! change in d per height [cm/m]
    
    d      = 10.0_r8**((log10(h)-p2)/p1)
    
    if(present(dddh))then
       dddh = 1.0_r8/(p1*10.0_r8**p2*d**(p1-1.0_r8))
    end if
    
    return
  end subroutine h2d_obrien
  
  ! ============================================================================
  
  subroutine h2d_martcano(h,p1,p2,p3,d,dddh)
    
    ! =========================================================================
    ! "d2h_martcano"
    ! "d to height via 3 parameter Michaelis-Menten following work at BCI
    ! by Martinez-Cano et al. 2016
    ! 
    ! h = (a*d**b)/(c+d**b)
    !
    ! d = [(h*c)/(a-h)]**(1/b)
    ! d = [(h*c)**(1/b)] / [(a-h)**(1/b)]
    ! d' = {[(h*c)**(1/b)]' [(a-h)**(1/b)] -  [(a-h)**(1/b)]'[(h*c)**(1/b)]} /
    !       [(a-h)**(1/b)]**2
    ! dddh = d' = {[(1/b)(h*c)**(1/b-1)] [(a-h)**(1/b)] -  
    !              [(1/b)(a-h)**(1/b-1)] [(h*c)**(1/b)]} /
    !             [(a-h)**(1/b)]**2
    ! 
    ! =========================================================================
    
    real(r8),intent(in)  :: h     ! plant height [m]
    real(r8),intent(in)  :: p1
    real(r8),intent(in)  :: p2
    real(r8),intent(in)  :: p3
    
    real(r8),intent(out)   :: d     ! plant diameter [cm]
    real(r8),intent(out),optional :: dddh  ! change in diameter per height [cm/m]
    
    d = ((h*p3)/(p1-h))**(1._r8/p2)
    
    if(present(dddh))then
       dddh =  (((1._r8/p2)*(h*p3)**(1._r8/p2-1._r8))*((p1-h)**(1._r8/p2)) - & 
            ((1._r8/p2)*(p1-h)**(1._r8/p2-1._r8))* ((h*p3)**(1._r8/p2)) ) / &
            ((p1-h)**(1._r8/p2))**2._r8
    end if
    return
  end subroutine h2d_martcano

  ! =============================================================================
  ! Specific diameter to crown area allometries
  ! =============================================================================

  
  subroutine carea_2pwr(d,spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max,c_area)

     ! ============================================================================
     ! Calculate area of ground covered by entire cohort. (m2)
     ! Function of DBH (cm) canopy spread (m/cm) and number of individuals. 
     ! ============================================================================

     real(r8),intent(in) :: d           ! diameter [cm]
     real(r8),intent(in) :: spread      ! site level relative spread score [0-1]
     real(r8),intent(in) :: d2bl_p2     ! parameter 2 in the diameter->bleaf allometry (exponent)
     real(r8),intent(in) :: d2bl_ediff  ! area difference factor in the diameter-bleaf allometry (exponent)
     real(r8),intent(in) :: d2ca_min    ! minimum diameter to crown area scaling factor
     real(r8),intent(in) :: d2ca_max    ! maximum diameter to crown area scaling factor
     real(r8),intent(out) :: c_area     ! crown area for one plant [m2]
     
     real(r8)            :: crown_area_to_dbh_exponent
     real(r8)            :: spreadterm  ! Effective 2bh to crown area scaling factor
     
     ! default is to use the same exponent as the dbh to bleaf exponent so that per-plant 
     ! canopy depth remains invariant during growth, but allowed to vary via the 
     ! allom_blca_expnt_diff term (which has default value of zero)
     crown_area_to_dbh_exponent = d2bl_p2 + d2bl_ediff
     
     ! ----------------------------------------------------------------------------------
     ! The function c_area is called during the process of canopy position demotion
     ! and promotion. As such, some cohorts are temporarily elevated to canopy positions
     ! that are outside the number of alloted canopy spaces.  Ie, a two story canopy
     ! may have a third-story plant, if only for a moment.  However, these plants
     ! still need to generate a crown area to complete the promotion, demotion process.
     ! So we allow layer index exceedence here and force it down to max.
     ! (rgk/cdk 05/2017)
     ! ----------------------------------------------------------------------------------
     
     ! apply site-level spread elasticity to the cohort crown allometry term
    
     spreadterm = spread * d2ca_max + (1._r8 - spread) * d2ca_min
     
     c_area = spreadterm * d ** crown_area_to_dbh_exponent
     
  end subroutine carea_2pwr
  
  ! =========================================================================

  subroutine StructureResetOfDH( bdead, ipft, canopy_trim, d, h )

     ! =========================================================================
     ! This subroutine estimates the diameter based on the structural biomass
     ! using the allometric functions. Since allometry is specified with diameter
     ! as the independant variable, we must do this through a search algorithm.
     ! Here, we keep searching until the difference between actual structure and
     ! the predicted structure based on the searched diameter is within a tolerance.
     ! T
     ! ============================================================================

     use FatesConstantsMod     , only : calloc_abs_error
     ! Arguments

     real(r8),intent(in)           :: bdead ! actual bdead [kgC]
     integer(i4),intent(in)        :: ipft  ! PFT index
     real(r8),intent(in)           :: canopy_trim
     real(r8),intent(inout)        :: d     ! plant diameter [cm]
     real(r8),intent(out)          :: h     ! plant height
     
     ! Locals
     real(r8)  :: bt_sap,dbt_sap_dd  ! target sap wood at current d
     real(r8)  :: bt_agw,dbt_agw_dd  ! target AG wood at current d
     real(r8)  :: bt_bgw,dbt_bgw_dd  ! target BG wood at current d
     real(r8)  :: bt_dead,dbt_dead_dd ! target struct wood at current d
     real(r8)  :: dd                  ! diameter increment for each step
     real(r8)  :: d_try               ! trial diameter
     real(r8)  :: bt_dead_try         ! trial structure biomasss
     real(r8)  :: dbt_dead_dd_try     ! trial structural derivative
     real(r8)  :: step_frac           ! step fraction
     integer   :: counter 
     real(r8), parameter :: step_frac0  = 0.9_r8
     integer, parameter  :: max_counter = 200

     call bsap_allom(d,ipft,canopy_trim,bt_sap,dbt_sap_dd)
     call bagw_allom(d,ipft,bt_agw,dbt_agw_dd)
     call bbgw_allom(d,ipft,bt_bgw,dbt_bgw_dd)
     call bdead_allom(bt_agw,bt_bgw, bt_sap, ipft, bt_dead, dbt_agw_dd, &
                      dbt_bgw_dd, dbt_sap_dd, dbt_dead_dd)

     ! This calculates a diameter increment based on the difference
     ! in structural mass and the target mass, and sets it to a fraction
     ! of the diameter increment
     counter = 0
     step_frac = step_frac0
     do while( (bdead-bt_dead) > calloc_abs_error )

        ! vulnerable to div0
        dd    = step_frac*(bdead-bt_dead)/dbt_dead_dd
        d_try = d + dd
        
        call bsap_allom(d_try,ipft,canopy_trim,bt_sap,dbt_sap_dd)
        call bagw_allom(d_try,ipft,bt_agw,dbt_agw_dd)
        call bbgw_allom(d_try,ipft,bt_bgw,dbt_bgw_dd)
        call bdead_allom(bt_agw,bt_bgw, bt_sap, ipft, bt_dead_try, dbt_agw_dd, &
              dbt_bgw_dd, dbt_sap_dd, dbt_dead_dd_try)

        ! Prevent overshooting
        if(bt_dead_try > (bdead+calloc_abs_error)) then
           step_frac = step_frac*0.5_r8
        else
           step_frac = step_frac0
           d         = d_try
           bt_dead   = bt_dead_try
           dbt_dead_dd = dbt_dead_dd_try
        end if
        counter = counter + 1
        if (counter>max_counter) then
           write(fates_log(),*) 'Having trouble converging on dbh reset'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
     end do

     call h_allom(d,ipft,h)
     if(counter>10)then
        write(fates_log(),*) 'dbh counter: ',counter
     end if

     ! At this point, the diameter, height and their target structural biomass
     ! should be pretty close to and greater than actual

     return
  end subroutine StructureResetOfDH

  ! =========================================================================
  
  subroutine cspline(x1,x2,y1,y2,dydx1,dydx2,x,y,dydx)
    
    ! ============================================================================
    ! This subroutine performs a cubic spline interpolation between known
    ! endpoints.  The endpoints have known coordinats and slopes
    !
    ! This routine may come in handy if we ever want to start using allometries
    ! that are different for juvenile and adult plants, and then connect
    ! the curves smoothly
    !
    ! ============================================================================
    
    ! Arguments
    
    real(r8),intent(in) :: x1     ! Lower endpoint independent
    real(r8),intent(in) :: x2     ! Upper endpoint independent
    real(r8),intent(in) :: y1     ! Lower endpoint dependent
    real(r8),intent(in) :: y2     ! Upper endpoint dependent
    real(r8),intent(in) :: dydx1  ! Lower endpoint slope
    real(r8),intent(in) :: dydx2  ! Upper endpoint slope
    real(r8),intent(in) :: x      ! Independent
    real(r8),intent(out) :: y     ! Dependent
    real(r8),intent(out) :: dydx  ! Slope
    
    ! Temps
    real(r8) :: t
    real(r8) :: a
    real(r8) :: b
    
    t = (x-x1)/(x2-x1)
    a = dydx1*(x2-x1) - (y2-y1)
    b = -dydx2*(x2-x1) + (y2-y1)
    
    y    = (1.0_r8-t)*y1 + t*y2 + t*(1.0_r8-t)*(a*(1.0_r8-t) + b*t)
    dydx = (y2-y1)/(x2-x1) + (1.0_r8-2.0_r8*t)*(a*(1.0_r8-t)+b*t)/(x2-x1) + t*(1.0_r8-t)*(b-a)/(x2-x1)
    return
  end subroutine cspline
  
end module FatesAllometryMod
