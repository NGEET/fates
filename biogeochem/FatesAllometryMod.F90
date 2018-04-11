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
! h:  height [m]
! bag:  biomass above ground [kgC]  (aka AGB)
! blmax:  biomass in leaves when leaves are "on allometry"
!         this also is the "pre-trimmed" value, which is the maximum
!         or potential leaf mass a tree may have [kgC]
! bcr: biomass in coarse roots [kgC] (belowground sap+dead wood, no fines)
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
!
! OPEN QUESTIONS:
!      SHOULD SAPWOOD ALLOMETRY BE EFFECTED BY TRIMMING, OR BE OFF OF BLMAX?
!      WHAT POOLS DO WE ASSUME ARE CONTAINED in AGB ALLOMETRY?
!      WE ARE NOT EXTENDING SAPWOOD BELOW GROUND, WE PROBABLY SHOULD
!
!
! Carbon Pool Configurations are as follows, and assume a constant proportionality
!      between above and below-ground pools.  Sapwood (bsap) is both above and below
!      ground. Above ground biomass contains above ground dead wood, above ground
!      sapwood and leaves.  Coarse roots (bcr) contain only the below ground
!      deadwood, not the fine roots or the sapwood.  Coarse roots are typically
!      a proportion of above ground deadwood.
!      Leaf biomass, height and above ground biomass typically have non-linear
!      allometry models.  The default for sapwood is the pipe model.
! 
!  We ignore leaf biomass contributions in allometry to AGB:
! 
!  bag   = (bdead+bsap)*agb_frac 
!  bdead = bag - (bsap*agb_frac) + bcr
!  bcr   = bdead * (1-agb_frac) = (bag - (bsap*agb_frac) + bcr)*(1-agb_frac)
!
!
! Initial Implementation: Ryan Knox July 2017
!
!===============================================================================

module FatesAllometryMod

  ! If this is a unit-test, these globals will be provided by a wrapper

  use EDPFTvarcon      , only : EDPftvarcon_inst
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use FatesGlobals     , only : fates_log
  use FatesGlobals     , only : endrun => fates_endrun

  implicit none

  private
  public :: h2d_allom     ! Generic height to diameter wrapper
  public :: h_allom       ! Generic diameter to height wrapper
  public :: bag_allom     ! Generic AGB wrapper
  public :: blmax_allom   ! Generic maximum leaf biomass wrapper
  public :: bleaf         ! Generic actual leaf biomass wrapper
  public :: bsap_allom    ! Generic sapwood wrapper
  public :: bcr_allom     ! Generic coarse root wrapper
  public :: bfineroot     ! Generic actual fine root biomass wrapper
  public :: bdead_allom   ! Generic bdead wrapper
  public :: carea_allom   ! Generic crown area wrapper

  character(len=*), parameter :: sourcefile = __FILE__

  ! If testing b4b with older versions, do not remove sapwood
  ! Our old methods with saldarriaga did not remove sapwood from the
  ! bdead pool.  But newer allometries are providing total agb
  ! which includes sapwood. Although a small quantity, it needs to be removed
  ! from the agb pool.
  ! Additionally, our calculation of sapwood biomass may be missing some unite conversions
  !
  logical,parameter :: test_b4b = .true.

contains

  ! ============================================================================
  ! Parameter Checks
  ! ============================================================================

  ! Checks to make sure parameters are not within expected ranges for each
  ! functions

  ! Check to make sure Martinez-Cano height cap is not on, or explicitly allowed


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
  
  subroutine bag_allom(d,h,ipft,bag,dbagdd)


    real(r8),intent(in)    :: d       ! plant diameter [cm]
    real(r8),intent(in)    :: h       ! plant height [m]
    integer(i4),intent(in) :: ipft    ! PFT index
    real(r8),intent(out)   :: bag     ! plant height [m]
    real(r8),intent(out),optional   :: dbagdd  ! change in agb per diameter [kgC/cm]

    real(r8)               :: hj      ! height (dummy arg)
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
         call h_allom(d,ipft,hj,dhdd)
         call dh2bag_salda(d,h,dhdd,p1,p2,p3,p4,wood_density,c2b,agb_frac,bag,dbagdd) 
      case (2) !"2par_pwr")
         ! Switch for woodland dbh->drc
         call d2bag_2pwr(d,p1,p2,c2b,bag,dbagdd)
      case (3) !"chave14") 
         call h_allom(d,ipft,hj,dhdd)
         call dh2bag_chave2014(d,h,dhdd,p1,p2,wood_density,c2b,bag,dbagdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined AGB allometry was specified: ',allom_amode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      
    end associate
    return
  end subroutine bag_allom

  ! ============================================================================
  ! Generic diameter to maximum leaf biomass interface
  ! ============================================================================
  
  subroutine blmax_allom(d,h,ipft,blmax,dblmaxdd)

    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
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
        
  subroutine bleaf(d,h,ipft,canopy_trim,bl,dbldd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target bleaf
    ! based on trimming. Because trimming
    ! is not allometry and rather an emergent property,
    ! this routine is not name-spaced with allom_
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d             ! plant diameter [cm]
    real(r8),intent(in)    :: h             ! plant height [m]
    integer(i4),intent(in) :: ipft          ! PFT index
    real(r8),intent(in)    :: canopy_trim   ! trimming function
    real(r8),intent(out)   :: bl            ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dbldd  ! change leaf bio per diameter [kgC/cm]
    
    real(r8) :: blmax
    real(r8) :: dblmaxdd
    
    call blmax_allom(d,h,ipft,blmax,dblmaxdd)
    
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
    real(r8) :: blmax
    real(r8) :: dblmaxdd
    real(r8) :: bag
    real(r8) :: dbagdd
    
    select case(int(EDPftvarcon_inst%allom_smode(ipft)))
       ! ---------------------------------------------------------------------
       ! Currently both sapwood area proportionality methods use the same
       ! machinery.  The only differences are related to the parameter
       ! checking at the beginning.  For constant proportionality, the slope
       ! of the la:sa to diameter line is zero.
       ! ---------------------------------------------------------------------
    case(1,2) !"constant","dlinear") 

       call h_allom(d,ipft,h,dhdd)
       if(test_b4b)then
          call bleaf(d,h,ipft,canopy_trim,blmax,dblmaxdd)
       else
          call blmax_allom(d,h,ipft,blmax,dblmaxdd)
       end if
       call bag_allom(d,h,ipft,bag,dbagdd)
       call bsap_dlinear(d,h,dhdd,blmax,dblmaxdd,bag,dbagdd,ipft,bsap,dbsapdd)
    case DEFAULT
       write(fates_log(),*) 'An undefined sapwood allometry was specified: ', &
            EDPftvarcon_inst%allom_smode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    return
  end subroutine bsap_allom
  
  ! ============================================================================
  ! Generic coarse root biomass interface
  ! ============================================================================

  subroutine bcr_allom(d,h,ipft,bcr,dbcrdd)


    real(r8),intent(in)           :: d         ! plant diameter [cm]
    real(r8),intent(in)           :: h         ! plant height [m]
    integer(i4),intent(in)        :: ipft      ! PFT index
    real(r8),intent(out)          :: bcr       ! coarse root biomass [kgC]
    real(r8),intent(out),optional :: dbcrdd    ! change croot bio per diam [kgC/cm]
    
    real(r8)    :: bag       ! above ground biomass [kgC]
    real(r8)    :: dbagdd    ! change in agb per diameter [kgC/cm]
    
    select case(int(EDPftvarcon_inst%allom_cmode(ipft)))
    case(1) !"constant")
       call bag_allom(d,h,ipft,bag,dbagdd)
       call bcr_const(d,bag,dbagdd,ipft,bcr,dbcrdd)
    case DEFAULT
       write(fates_log(),*) 'An undefined coarse root allometry was specified: ', &
            EDPftvarcon_inst%allom_cmode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    return
  end subroutine bcr_allom

  ! ============================================================================
  ! Fine root biomass allometry wrapper
  ! ============================================================================
  
  subroutine bfineroot(d,h,ipft,canopy_trim,bfr,dbfrdd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target fineroot biomass
    ! based on functions that may or may not have prognostic properties. 
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d              ! plant diameter [cm]
    real(r8),intent(in)    :: h              ! plant height [m]
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
    case(1) ! "constant proportionality with bleaf"
       
       call blmax_allom(d,h,ipft,blmax,dblmaxdd)
       call bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
       bfr    = bfrmax * canopy_trim
       if(present(dbfrdd))then
          dbfrdd = dbfrmaxdd * canopy_trim
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
  ! Dead biomass interface
  ! ============================================================================
  
  subroutine bdead_allom(bag,bcr,bsap,ipft,bdead,dbagdd,dbcrdd,dbsapdd,dbdeaddd)


    real(r8),intent(in)  :: bag       ! agb [kgC]
    real(r8),intent(in)  :: bcr       ! coarse root biomass [kgC]
    real(r8),intent(in)  :: bsap      ! sapwood biomass [kgC]
    integer(i4),intent(in) :: ipft    ! PFT index
    real(r8),intent(out) :: bdead     ! dead biomass (heartw/struct) [kgC]
    
    real(r8),intent(in),optional  :: dbagdd    ! change in agb per d [kgC/cm]
    real(r8),intent(in),optional  :: dbcrdd    ! change in croot per d [kgC/cm]
    real(r8),intent(in),optional  :: dbsapdd   ! change in bsap per d [kgC/cm]
    real(r8),intent(out),optional :: dbdeaddd  ! change in bdead per d [kgC/cm]
    
    ! bdead is diagnosed as the mass balance from all other pools
    ! and therefore, no options are necessary
    
    associate( agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft))

      select case(int(EDPftvarcon_inst%allom_amode(ipft)))
      case(1) ! Saldariagga mass allometry originally calculated bdead directly.
              ! we assume proportionality between bdead and bag
       
         bdead = bag/agb_fraction 
         if(present(dbagdd) .and. present(dbdeaddd))then
            dbdeaddd = dbagdd/agb_fraction
         end if
         
      case(2,3)
         
         bdead = bag+bcr-bsap
         if(present(dbagdd) .and. present(dbcrdd) .and. &
            present(dbdeaddd) .and. present(dbsapdd) )then
            dbdeaddd = dbagdd+dbcrdd-dbsapdd
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
  ! Specific bcr relationships
  ! ============================================================================
  
  subroutine bcr_const(d,bag,dbagdd,ipft,bcr,dbcrdd)

    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: bag       ! above ground biomass [kg]
    real(r8),intent(in)    :: dbagdd    ! change in agb per diameter [kg/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bcr       ! coarse root biomass [kg]
    real(r8),intent(out),optional :: dbcrdd    ! change croot bio per diam [kg/cm]

    associate( agb_fraction => EDPftvarcon_inst%allom_agb_frac(ipft) )
      
      ! btot = bag + bcr
      ! bag = btot*agb_fraction
      ! bag/agb_fraction = bag + bcr
      ! bcr = bag*(1/agb_fraction-1)
      bcr = bag*(1.0_r8/agb_fraction-1.0_r8)
      
      ! Derivative
      ! dbcr/dd = dbcr/dbag * dbag/dd
      if(present(dbcrdd))then
         dbcrdd = (1.0_r8/agb_fraction-1.0_r8)*dbagdd
      end if
      
    end associate
    return
  end subroutine bcr_const

  ! ============================================================================
  ! Specific d2bsap relationships
  ! ============================================================================
  
  subroutine bsap_dlinear(d,h,dhdd,blmax,dblmaxdd,bag,dbagdd,ipft,bsap,dbsapdd)
    
    use FatesConstantsMod, only : g_per_kg
    use FatesConstantsMod, only : cm2_per_m2
    use FatesConstantsMod, only : kg_per_Megag
    
    ! -------------------------------------------------------------------------
    ! Calculate sapwood biomass based on leaf area to sapwood area
    ! proportionality.  In this function, the leaftosapwood area is a function
    ! of plant size, see Calvo-Alvarado and Bradley Christoferson
    ! In this case: parameter latosa (from constant proportionality)
    !   is the intercept of the diameter function.
    !
    ! For very small plants, the fraction can get very large, so cap the amount
    ! of sapwood at X! of agb-bleaf
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    real(r8),intent(in)    :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    real(r8),intent(in)    :: bag       ! aboveground biomass [kgC]
    real(r8),intent(in)    :: dbagdd    ! change in agb per diam [kgC/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per diameter [kgC/cm]
    

    real(r8)               :: latosa    ! applied leaf area to sap area 
                                        ! may or may not contain diameter correction
    real(r8)               :: hbl2bsap  ! sapwood biomass per lineal height and kg of leaf
    
    ! Constrain sapwood to be no larger than 75% of total agb
    real(r8),parameter :: max_agbfrac = 0.75_r8 
    
    associate ( latosa_int => EDPftvarcon_inst%allom_latosa_int(ipft), &
                latosa_slp => EDPftvarcon_inst%allom_latosa_slp(ipft), &
                sla        => EDPftvarcon_inst%slatop(ipft), &
                wood_density => EDPftvarcon_inst%wood_density(ipft), &
                c2b          => EDPftvarcon_inst%c2b(ipft))

      ! ------------------------------------------------------------------------
      ! Calculate sapwood biomass per linear height and kgC of leaf [m-1]
      ! Units: 
      ! (1/latosa)* slatop*    gtokg    *   cm2tom2     / c2b   * mg2kg  * dens
      ! [cm2/m2]*[m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3]
      !            ->[cm2/gC]
      !                      ->[cm2/kgC]
      !                                   ->[m2/kgC]
      !                                                 ->[m2/kg]
      !                                                          ->[m2/Mg]
      !                                                                  ->[/m]
      ! ------------------------------------------------------------------------

      if(test_b4b) then
         
         bsap = blmax * latosa_int * h
         
         if(present(dbsapdd))then
            dbsapdd = latosa_int*(h*dblmaxdd + blmax*dhdd)
         end if
         
      else
         
         latosa = latosa_int + d*latosa_slp
         hbl2bsap = sla*g_per_kg*wood_density*kg_per_Megag/(latosa*c2b*cm2_per_m2 )
         
         ! Force sapwood to be less than a maximum fraction of total alive biomass
         ! (this comes into play typically in very small plants)
         bsap = min(max_agbfrac*bag,hbl2bsap * h * blmax)
         
         ! Derivative
         ! dbldmaxdd is deriv of blmax wrt dbh (use directives to check oop)
         ! dhdd is deriv of height wrt dbh (use directives to check oop)
         if(present(dbsapdd))then
            if (bsap<max_agbfrac*bag) then
               dbsapdd = hbl2bsap*(h*dblmaxdd + blmax*dhdd)
            else
               dbsapdd = max_agbfrac*dbagdd
            end if
         end if
      end if
    end associate
    return
  end subroutine bsap_dlinear
  
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
          if(test_b4b)then
             dblmaxdd = p1*p2 * d**(p2) * rho**p3
          else
             dblmaxdd = p1*p2 * d**(p2-1.0_r8) * rho**p3
          end if
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
  
  subroutine dh2bag_chave2014(d,h,dhdd,p1,p2,wood_density,c2b,bag,dbagdd)
    
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
    ! bag:   Total above ground biomass [kgC]
    !
    ! =========================================================================
    
    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: h       ! plant height [m]
    real(r8),intent(in)  :: dhdd    ! change in height wrt diameter
    real(r8),intent(in)  :: p1  ! allometry parameter 1
    real(r8),intent(in)  :: p2  ! allometry parameter 2
    real(r8),intent(in)  :: wood_density
    real(r8),intent(in)  :: c2b
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]
    
    real(r8) :: dbagdd1,dbagdd2,dbagdd3
    
    bag   = (p1 * (wood_density*d**2.0_r8*h)**p2)/c2b
    
    
    if(present(dbagdd))then
       ! Need the the derivative of height to diameter to
       ! solve the derivative of agb with height
       dbagdd1  = (p1*wood_density**p2)/c2b
       dbagdd2  = p2*d**(2.0_r8*p2)*h**(p2-1.0_r8)*dhdd
       dbagdd3  = h**p2*2.0_r8*p2*d**(2.0_r8*p2-1.0_r8)
       dbagdd   = dbagdd1*(dbagdd2 + dbagdd3)
    end if
    
    return
  end subroutine dh2bag_chave2014
  
  subroutine d2bag_2pwr(d,p1,p2,c2b,bag,dbagdd)
    
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
    ! bag = exp(b0 + b1*ln(diameter))/c2b
    ! =========================================================================

    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: p1  ! allometry parameter 1
    real(r8),intent(in)  :: p2  ! allometry parameter 2
    real(r8),intent(in)  :: c2b     ! carbon to biomass multiplier ~2
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]
    
    bag    = (p1 * d**p2)/c2b
    if(present(dbagdd))then
       dbagdd = (p2*p1*d**(p2-1.0_r8))/c2b
    end if
    
    return
  end subroutine d2bag_2pwr
  
  
  subroutine dh2bag_salda(d,h,dhdd,p1,p2,p3,p4, &
       wood_density,c2b,allom_agb_frac,bag,dbagdd)
    
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
    real(r8),intent(out) :: bag     ! plant biomass [kgC/indiv]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]
    
    real(r8) :: term1,term2,term3
    
    
    bag = allom_agb_frac*p1*(h**p2)*(d**p3)*(wood_density**p4)

    ! Add sapwood calculation to this?
    
    ! bag     = a1 * h**a2 * d**a3 * r**a4
    ! dbag/dd = a1*r**a4 * d/dd (h**a2 * d**a3)
    ! dbag/dd = a1*r**a4 * [ d**a3 *d/dd(h**a2) + h**a2*d/dd(d**a3) ]
    ! dbag/dd = a1*r**a4 * [ d**a3 * a2*h**(a2-1)dh/dd + h**a2*a3*d**(a3-1)]
    
    ! From code
    !   dbag/dd =  a3 * a1 *(h**a2)*(d**(a3-1))* (r**a4) + a2*a1*(h**(a2-1))*(d**a3)*(r**a4)*dhdd 
    !   dbag/dd =  a1*r**a4 * [ d**a3 * a2* h**(a2-1)*dhdd +  a3 * (h**a2)*(d**(a3-1)) ]
    
    
    if(present(dbagdd)) then

       term1 = allom_agb_frac*p1*(wood_density**p4)
       term2 = (h**p2)*p3*d**(p3-1.0_r8)
       term3  = p2*h**(p2-1.0_r8)*(d**p3)*dhdd
       dbagdd = term1*(term2+term3)
       
    end if

    return
  end subroutine dh2bag_salda
  
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
 
  
  ! ===========================================================================
  
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
