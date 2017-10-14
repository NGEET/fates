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
! derivative with respect to diameter using a logical switch "dswitch"
! (derivative-switch). With one exception, the h2d function returns the
! change in diameter with respect to height.
!
! The name convention of the functions follows the form d2...  Which
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
! "on allometry" assumes the plant has leaves flushed, and has had
! sufficient carbon to meet maintenance turnover.
!
! The following traits are used:
! allom_hmode, integer, height allometry function type
! allom_lmode, integer, maximum leaf allometry function type
! allom_rmode, integer, maximum root allometry function type
! allom_amode, integer, AGB allometry function type
! allom_cmode, integer, coarse root allometry function type
! allom_smode, integer, sapwood allometry function type
! wood_density, real, mean stem wood specific gravity (heart,sap,bark)
! allom_latosa_int, real, leaf area to sap area ratio, intercept [m2/cm2]
! allom_latosa_slp, real, leaf area to sap area ratio, slope on diameter
! [m2/cm2/cm]
! c2b, real, carbon to biomass ratio (~2.0)
! allom_l2fr, real, fine root biomass per leaf biomass ratio [kgC/kgC]
! allom_agb_fraction, real, the fraction of stem above ground [-]
! allom_d2h1, real, parameter 1 for d2h allometry (intercept)
! allom_d2h2, real, parameter 2 for d2h allometry (slope)
! allom_d2h3, real, parameter 3 for d2h allometry (optional)
! eclim             influence parameter for d2h allometry (potentially not a parameter)
! allom_d2bl1, real, parameter 1 for d2bl allometry (intercept)
! allom_d2bl2, real, parameter 2 for d2bl allometry (slope)
! allom_d2bl3, real, parameter 3 for d2bl allometry (optional)
! allom_agb1
! allom_agb2
! allom_agb3
!
! h_max, real, maximum height of a functional type/group
! h_min, real, the height associated with newly recruited plant [m]
! dbh_min, real, the dbh associated with a newly recruited plant [cm]
! dbh_max, real, the diameter associated with maximum height [cm]
!                diagnosed from maxh using non-asymptotic functions
!
! Note - i4 types are expressed explicitly to accomodate unit testing calls
!        to this module
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

  character(len=*), parameter, private :: sourcefile = &
        __FILE__


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
      case (1) ! chave 2014
         call h2d_chave2014(h,p1,p2,p3,d,dddh)
      case (2)  ! poorter 2006
         call h2d_poorter2006(h,p1,p2,p3,d,dddh)
      case (3) ! 2 parameter power function
         call h2d_2pwr(h,p1,p2,d,dddh)
      case (4) ! Obrien et al. 199X BCI
         call h2d_obrien(h,p1,p2,d,dddh)
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

     ! Locals
     integer                 :: allom_hmode
     real(r8)                :: h_sap
     real(r8)                :: h_ad
     real(r8)                :: dhdd_sap
     real(r8)                :: dhdd_ad
     real(r8)                :: p1
     real(r8)                :: p2
     real(r8)                :: p3
                  
     associate( dbh_maxh => EDPftvarcon_inst%allom_dbh_maxheight(ipft), &
                p1       => EDPftvarcon_inst%allom_d2h1(ipft), &
                p2       => EDPftvarcon_inst%allom_d2h2(ipft), &
                p3       => EDPftvarcon_inst%allom_d2h3(ipft), &
                allom_hmode => EDPftvarcon_inst%allom_hmode(ipft))
       
       select case(int(allom_hmode))
       case (1)   ! "chave14")
          call d2h_chave2014(d,p1,p2,p3,dbh_maxh,h,dhdd)
       case (2)   ! "poorter06"
          call d2h_poorter2006(d,p1,p2,p3,dbh_maxh,h,dhdd)
       case (3)   ! "2parameter power function h=a*d^b "
          call d2h_2pwr(d,p1,p2,dbh_maxh,h,dhdd)
       case (4)   ! "obrien"
          call d2h_obrien(d,p1,p2,dbh_maxh,h,dhdd)
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


    associate( &
          p1 => EDPftvarcon_inst%allom_agb1(ipft), &
          p2 => EDPftvarcon_inst%allom_agb2(ipft), &
          p3 => EDPftvarcon_inst%allom_agb3(ipft), &
          p4 => EDPftvarcon_inst%allom_agb4(ipft), &
          wood_density => EDPftvarcon_inst%wood_density(ipft), &
          c2b => EDPftvarcon_inst%c2b(ipft), &
          agb_frac => EDPftvarcon_inst%allom_agb_frac(ipft), &
          allom_amode => EDPftvarcon_inst%allom_amode(ipft))

      select case(int(allom_amode))
      case (1) !"chave14") 
         call dh2bag_chave2014(d,h,ipft,p1,p2,wood_density,c2b,bag,dbagdd)
      case (2) !"2par_pwr")
         ! Switch for woodland dbh->drc
         call d2bag_2pwr(d,p1,p2,c2b,bag,dbagdd)
      case (3) !"salda")
         call dh2bag_salda(d,h,ipft,p1,p2,p3,p4,wood_density,c2b,agb_frac,bag,dbagdd)
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

    associate( &
          dbh_maxh    => EDPftvarcon_inst%allom_dbh_maxheight(ipft), &
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
         call dh2blmax_2pwr(d,ipft,p1,p2,c2b,blmax,dblmaxdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined leaf allometry was specified: ', &
               allom_lmode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
    end associate
    return
 end subroutine blmax_allom

  ! =====================================================================================

 subroutine bleaf(d,h,ipft,canopy_trim,bl,dbldd)
     
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target bleaf
    ! based on trimming and sla scaling. Because trimming
    ! is not allometry and rather an emergent property,
    ! this routine is not name-spaces with allom_
    ! -------------------------------------------------------------------------

    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(in)    :: canopy_trim ! trimming function
    real(r8),intent(out)   :: bl        ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dbldd  ! change leaf bio per diameter [kgC/cm]

    real(r8) :: blmax
    real(r8) :: dblmaxdd
    real(r8) :: slascaler

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
  
  subroutine bsap_allom(d,ipft,bsap,dbsapdd)

    real(r8),intent(in)           :: d         ! plant diameter [cm]
    integer(i4),intent(in)        :: ipft      ! PFT index
    real(r8),intent(out)          :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per d [kgC/cm]

    real(r8)    :: h         ! plant height [m]
    real(r8)    :: blmax     ! plant leaf biomass [kgC]
    real(r8)    :: dblmaxdd  ! chage in blmax per diam [kgC/cm]
    real(r8)    :: dhdd      ! change in height per diameter [m/cm]

    select case(int(EDPftvarcon_inst%allom_smode(ipft)))
       ! ---------------------------------------------------------------------
       ! Currently both sapwood area proportionality methods use the same
       ! machinery.  The only differences are related to the parameter
       ! checking at the beginning.  For constant proportionality, the slope
       ! of the la:sa to diameter line is zero.
       ! ---------------------------------------------------------------------
    case(1,2) !"constant","dlinear") 
       call h_allom(d,ipft,h,dhdd)
       call blmax_allom(d,h,ipft,blmax,dblmaxdd)
       call bsap_dlinear(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)
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

    call bag_allom(d,h,ipft,bag,dbagdd)

    select case(int(EDPftvarcon_inst%allom_cmode(ipft)))
    case(1) !"constant")
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

  subroutine bdead_allom(bag,bcr,bsap,bdead, &
                         dbagdd,dbcrdd,dbsapdd,dbdeaddd)

    
    real(r8),intent(in)  :: bag       ! agb [kgC]
    real(r8),intent(in)  :: bcr       ! coarse root biomass [kgC]
    real(r8),intent(in)  :: bsap      ! sapwood biomass [kgC]
    real(r8),intent(out) :: bdead     ! dead biomass (heartw/struct) [kgC]
    
    real(r8),intent(in),optional  :: dbagdd    ! change in agb per d [kgC/cm]
    real(r8),intent(in),optional  :: dbcrdd    ! change in croot per d [kgC/cm]
    real(r8),intent(in),optional  :: dbsapdd   ! change in bsap per d [kgC/cm]
    real(r8),intent(out),optional :: dbdeaddd  ! change in bdead per d [kgC/cm]

    ! If testing b4b with older versions, do not remove sapwood
    ! Our old methods with saldarriaga did not remove sapwood from the
    ! bdead pool.  But newer allometries are providing total agb
    ! which includes sapwood. Although a small quantity, it needs to be removed
    ! from the agb pool.
    logical,parameter :: test_b4b = .true.
    
    ! bdead is diagnosed as the mass balance from all other pools
    ! and therefore, no options are necessary
    
    if(test_b4b) then
       bdead = bag+bcr
    else
       bdead = bag+bcr-bsap
    end if
    
    if(present(dbagdd) .and. present(dbcrdd) .and. present(dbsapdd) .and. present(dbdeaddd) )then
       if(test_b4b) then
          dbdeaddd = dbagdd+dbcrdd
       else
          dbdeaddd = dbagdd+dbcrdd-dbsapdd
       end if
    end if
    
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

  subroutine bsap_dlinear(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)


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
    real(r8),intent(in)    :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per diameter [kgC/cm]

    real(r8)               :: latosa    ! m2 leaf area per cm2 sap area
    real(r8)               :: hbl2bsap  ! sapwood biomass per lineal height
                                          ! and kg of leaf
    real(r8)               :: bag       ! aboveground biomass [kgC]
    real(r8)               :: dbagdd    ! change in agb per diam [kgC/cm]

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
      
      latosa = latosa_int + d*latosa_slp

      hbl2bsap = sla*g_per_kg*wood_density*kg_per_Megag/(latosa*c2b*cm2_per_m2 )

      call bag_allom(d,h,ipft,bag,dbagdd)

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
       blmax    = p1 * dbh_maxh**p2 * rho**p3
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
    ! Notes:
    ! From King et al. 1990 at BCI for saplings
    !
    ! log(bl) = a2 + b2*log(h)
    ! bl = exp(a2) * h**b2
    !
    ! and:
    !
    ! log(d) = a1 + b1*log(h)
    ! d = exp(a1) * h**b1
    ! h = (1/exp(a1)) * d**(1/b1)
    !
    ! bl = exp(a2) * ((1/exp(a1)) * d**(1/b1))**b2
    ! bl = exp(a2) * (1/exp(a1))**b2 * d**(b2/b1)
    ! bl = p1 * d ** p2
    ! where: p1 = exp(a2) * (1/exp(a1))**b2 = 0.0201198
    !        p2 = (b2/b1) = 3.1791044
    ! For T. tuberculata (canopy tree):
    ! a1 = -0.0704, b1 = 0.67
    ! a2 = -4.056,  b2 = 2.13
    ! ======================================================================
    
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    real(r8),intent(in)  :: p1        ! parameter 1
    real(r8),intent(in)  :: p2        ! parameter 2
    real(r8),intent(in)  :: c2b       ! carbon to biomass multiplier

    real(r8),intent(out) :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]

    blmax    = p1*d**p2 / c2b

    if(present(dblmaxdd))then
       dblmaxdd = p1*p2 *d**(p2-1.0_r8) / c2b
    end if

    return
 end subroutine d2blmax_2pwr

  ! ===========================================================================
  
  subroutine dh2blmax_2pwr(d,ipft,p1,p2,c2b,blmax,dblmaxdd)
    
    ! -------------------------------------------------------------------------
    ! This formulation is very similar to d2blmax_2pwr
    ! The difference is that for very large trees that have reached peak
    ! height, we calculate an effective diameter to estimate the leaf biomass.
    ! The effective diameter is the diameter that matches the height on
    ! the non-capped portion of the plants height-diameter curve. For pft's
    ! that have naturally asymptotic formulations (Poorter, Michaeles-Menten, 
    ! etc)
    ! This will render the same results as d2blmax_2pwr, as the effective 
    ! diameter equals the actual diameter.  But for hard caps and logistic 
    ! caps, this will prevent trees with huge diameters and non-emergent 
    ! heights to have reasonable leaf levels.
    ! --------------------------------------------------------------------------

    real(r8),intent(in)  :: d         ! plant diameter [cm]
    integer,intent(in)   :: ipft      ! pft index
    real(r8),intent(in)  :: p1        ! parameter 1
    real(r8),intent(in)  :: p2        ! parameter 2
    real(r8),intent(in)  :: c2b       ! carbon 2 biomass multiplier

    real(r8),intent(out) :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]

    real(r8)             :: h         ! plant height
    real(r8)             :: dhdd      ! height to diameter differential
    real(r8)             :: dbh_eff   ! effective diameter
    real(r8)             :: dddh_eff  ! effective diameter to height differential
    real(r8)             :: ddeffdd   ! effective diameter to diameter differential


    ! This call is needed to calculate the rate of change of
    ! the actual h with d
    call h_allom(d,ipft,h,dhdd)
    call h2d_allom(h,ipft,dbh_eff,dddh_eff)
         
    ! This is the rate of change of the effective diameter
    ! with respect to the actual diameter (1.0 in non-height capped)
    ddeffdd  = dddh_eff * dhdd
    blmax    = p1*dbh_eff**p2/c2b

    ! If this plant has reached its height cap, then it is not
    ! adding leaf mass.  In this case, dhdd = 0
    if(present(dblmaxdd))then
       if(dhdd>0.0_r8) then
          dblmaxdd = p1*p2*dbh_eff**(p2-1.0_r8) / c2b * ddeffdd
       else
          dblmaxdd = 0.0_r8
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
    
    !eclim: Chave's climatological influence parameter "E"

    
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
  ! Diameter 2 above-ground biomass
  ! =========================================================================

  subroutine dh2bag_chave2014(d,h,ipft,d2bag1,d2bag2,wood_density,c2b,bag,dbagdd)

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
    integer ,intent(in)  :: ipft    ! plant pft
    real(r8),intent(in)  :: d2bag1  ! allometry parameter 1
    real(r8),intent(in)  :: d2bag2  ! allometry parameter 2
    real(r8),intent(in)  :: wood_density
    real(r8),intent(in)  :: c2b
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]

    real(r8) :: hj,dhdd
    real(r8) :: dbagdd1,dbagdd2,dbagdd3

    bag   = (d2bag1 * (wood_density*d**2.0_r8*h)**d2bag2)/c2b


    if(present(dbagdd))then
       ! Need the the derivative of height to diameter to
       ! solve the derivative of agb with height
       call h_allom(d,ipft,hj,dhdd)
       
       dbagdd1  = (d2bag1*wood_density**d2bag2)/c2b
       dbagdd2  = d2bag2*d**(2.0_r8*d2bag2)*h**(d2bag2-1.0_r8)*dhdd
       dbagdd3  = h**d2bag2*2.0_r8*d2bag2*d**(2.0_r8*d2bag2-1.0_r8)
       dbagdd   = dbagdd1*(dbagdd2 + dbagdd3)
    end if
       
    return
 end subroutine dh2bag_chave2014

 subroutine d2bag_2pwr(d,d2bag1,d2bag2,c2b,bag,dbagdd)

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
    real(r8),intent(in)  :: d2bag1  ! allometry parameter 1
    real(r8),intent(in)  :: d2bag2  ! allometry parameter 2
    real(r8),intent(in)  :: c2b     ! carbon to biomass multiplier ~2
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]
    
    bag    = (d2bag1 * d**d2bag2)/c2b
    if(present(dbagdd))then
       dbagdd = (d2bag2*d2bag1*d**(d2bag2-1.0_r8))/c2b
    end if
    
    return
 end subroutine d2bag_2pwr
  
  
  subroutine dh2bag_salda(d,h,ipft,d2bag1,d2bag2,d2bag3,d2bag4, &
                          wood_density,c2b,allom_agb_frac,bag,dbagdd)

    ! --------------------------------------------------------------------
    ! Calculate stem biomass from height(m) dbh(cm) and wood density(g/cm3)
    ! default params using allometry of J.G. Saldarriaga et al 1988 - Rio Negro
    ! Journal of Ecology vol 76 p938-958  
    ! Saldarriaga 1988 provided calculations on total dead biomass
    ! So here, we calculate total dead, and then call and remove
    ! coarse root and sapwood. We ignore fineroot and leaf 
    ! in the calculations
    ! --------------------------------------------------------------------

    
    real(r8),intent(in)  :: d             ! plant diameter [cm]
    real(r8),intent(in)  :: h             ! plant height [m]
    integer(i4),intent(in) :: ipft        ! PFT index
    real(r8),intent(in) :: d2bag1         !    = 0.06896_r8
    real(r8),intent(in) :: d2bag2         !    = 0.572_r8
    real(r8),intent(in) :: d2bag3         !    = 1.94_r8
    real(r8),intent(in) :: d2bag4         !    = 0.931_r8
    real(r8),intent(in) :: c2b            ! carbon 2 biomass ratio
    real(r8),intent(in) :: wood_density   
    real(r8),intent(in) :: allom_agb_frac

    real(r8),intent(out) :: bag     ! plant biomass [kgC/indiv]
    real(r8),intent(out),optional :: dbagdd  ! change in agb per diameter [kgC/cm]

    real(r8) :: term1,term2,term3,hj,dhdd


    bag = allom_agb_frac*d2bag1*(h**d2bag2)*(d**d2bag3)*(wood_density**d2bag4)
    
    ! Add sapwood calculation to this

    ! bag     = a1 * h**a2 * d**a3 * r**a4
    ! dbag/dd = a1*r**a4 * d/dd (h**a2*d**a3)
    ! dbag/dd = a1*r**a4 * [ d**a3 *d/dd(h**a2) + h**a2*d/dd(d**a3) ]
    ! dbag/dd = a1*r**a4 * [ d**a3 * a2*h**(a2-1)dh/dd + h**a2*a3*d**(a3-1)]

    if(present(dbagdd)) then
       term1 = allom_agb_frac*d2bag1*(wood_density**d2bag4)
       term2 = (h**d2bag2)*d2bag3*d**(d2bag3-1.0_r8)
       
       call h_allom(d,ipft,hj,dhdd)
       term3 = d2bag2*h**(d2bag2-1)*(d**d2bag3)*dhdd
       dbagdd   = term1*(term2+term3)
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

  ! ===========================================================================

  subroutine cspline(x1,x2,y1,y2,dydx1,dydx2,x,y,dydx)

    ! ============================================================================
    ! This subroutine performs a cubic spline interpolation between known
    ! endpoints.  The endpoints have known coordinats and slopes
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
