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
! allom_la_per_sa_int, leaf area per sapwood area, intercept [m2/cm2]
! allom_la_per_sa_slp, leaf area per sapwood area, slope on diameter [m2/cm2/cm]
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

  use PRTParametersMod, only  : prt_params
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : g_per_kg 
  use FatesConstantsMod, only : cm2_per_m2
  use FatesConstantsMod, only : kg_per_Megag
  use FatesConstantsMod, only : calloc_abs_error
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : itrue
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : pi_const
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use FatesGlobals     , only : fates_log
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : FatesWarn,N2S,A2S,I2S
  use EDParamsMod      , only : nlevleaf,dinc_vai,dlower_vai
  use DamageMainMod    , only : GetCrownReduction

  implicit none

  private
  public :: h2d_allom     ! Generic height to diameter wrapper
  public :: h_allom       ! Generic diameter to height wrapper
  public :: bagw_allom    ! Generic AGWB (above grnd. woody bio) wrapper
  public :: blmax_allom   ! Generic maximum leaf biomass wrapper
  public :: bleaf         ! Generic actual leaf biomass wrapper
  public :: storage_fraction_of_target ! storage as fraction of leaf biomass
  public :: bsap_allom    ! Generic sapwood wrapper
  public :: bbgw_allom    ! Generic coarse root wrapper
  public :: bfineroot     ! Generic actual fine root biomass wrapper
  public :: bdead_allom   ! Generic bdead wrapper
  public :: carea_allom   ! Generic crown area wrapper
  public :: bstore_allom  ! Generic maximum storage carbon wrapper
  public :: decay_coeff_vcmax ! vertical canopy decay rate, scaled on vcmax
  public :: ForceDBH      ! Method to set DBH to sync with structure
                          ! or fineroot biomass
  public :: CheckIntegratedAllometries
  public :: CrownDepth
  public :: set_root_fraction  ! Generic wrapper to calculate normalized
                               ! root profiles
  public :: leafc_from_treelai ! Calculate target leaf carbon for a given treelai for SP mode
  public :: VegAreaLayer

  public :: tree_lai_sai       ! LAI and SAI calculations must work together, thus they
                               ! should never be called separately

  logical         , parameter :: verbose_logging = .false.
  character(len=*), parameter :: sourcefile = __FILE__

  
  logical, parameter :: debug = .false.

  character(len=1024) :: warn_msg   ! for defining a warning message
  
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
   
  subroutine CheckIntegratedAllometries(dbh,ipft,crowndamage, &
       canopy_trim, elongf_leaf, elongf_fnrt, elongf_stem, &
       l2fr, bl,bfr,bsap,bstore,bdead, &
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
     integer,intent(in)  :: crowndamage ! crowndamage [1: undamaged, >1 damaged]
     real(r8),intent(in) :: canopy_trim ! trimming function
     real(r8),intent(in) :: elongf_leaf ! Leaf elongation factor
     real(r8),intent(in) :: elongf_fnrt ! Fine-root elongation factor
     real(r8),intent(in) :: elongf_stem ! Stem elongation factor
     real(r8),intent(in) :: l2fr   ! leaf to fine-root biomass multiplier (fr/leaf)
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
     real(r8) :: asap_diag         ! sapwood area (dummy) [m2]
     real(r8) :: bdead_diag        ! diagnosed structural biomass [kgC]
     real(r8) :: bstore_diag       ! diagnosed storage biomass [kgC]
     real(r8) :: bagw_diag         ! diagnosed agbw [kgC]
     real(r8) :: bbgw_diag         ! diagnosed below ground wood [kgC]



     l_pass = .true.  ! Default assumption is that step passed

     if (grow_leaf) then
        call bleaf(dbh,ipft,crowndamage, canopy_trim, elongf_leaf, bl_diag)
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
        call bfineroot(dbh,ipft,canopy_trim,l2fr, elongf_fnrt, bfr_diag)
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
        call bsap_allom(dbh,ipft,crowndamage, canopy_trim, elongf_stem, asap_diag,bsap_diag)
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
        call bstore_allom(dbh,ipft,crowndamage, canopy_trim,bstore_diag)
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
        call bsap_allom(dbh,ipft,crowndamage, canopy_trim, elongf_stem,asap_diag,bsap_diag)
        call bagw_allom(dbh,ipft,crowndamage, elongf_stem, bagw_diag)
        call bbgw_allom(dbh,ipft, elongf_stem,bbgw_diag)
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

    associate(  p1          => prt_params%allom_d2h1(ipft), &
                p2          => prt_params%allom_d2h2(ipft), &
                p3          => prt_params%allom_d2h3(ipft), &
                allom_hmode => prt_params%allom_hmode(ipft))

      select case(allom_hmode)
      case (1) ! O'Brien et al 1995, BCI
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
    
    associate( dbh_maxh => prt_params%allom_dbh_maxheight(ipft), &
               p1       => prt_params%allom_d2h1(ipft), &
               p2       => prt_params%allom_d2h2(ipft), &
               p3       => prt_params%allom_d2h3(ipft), &
               allom_hmode => prt_params%allom_hmode(ipft))
      
      select case(allom_hmode)
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
  
  subroutine bagw_allom(d,ipft,crowndamage, elongf_stem, bagw,dbagwdd)

    use DamageMainMod, only : GetCrownReduction
    use FatesParameterDerivedMod, only : param_derived
    
    real(r8),intent(in)    :: d       ! plant diameter [cm]
    integer(i4),intent(in) :: ipft    ! PFT index
    integer(i4),intent(in) :: crowndamage ! crowndamage [1: undamaged, >1: damaged]
    real(r8),intent(in)    :: elongf_stem ! Stem elongation factor
    real(r8),intent(out)   :: bagw    ! biomass above ground woody tissues
    real(r8),intent(out),optional :: dbagwdd  ! change in agbw per diameter [kgC/cm]

    real(r8)               :: h       ! height
    real(r8)               :: dhdd    ! change in height wrt d
    real(r8)               :: crown_reduction  ! crown reduction from damage
    real(r8)               :: branch_frac ! fraction of aboveground woody biomass in branches
   
    associate( p1           => prt_params%allom_agb1(ipft), &
               p2           => prt_params%allom_agb2(ipft), &
               p3           => prt_params%allom_agb3(ipft), &
               p4           => prt_params%allom_agb4(ipft), &
               wood_density => prt_params%wood_density(ipft), &
               c2b          => prt_params%c2b(ipft), &
               agb_frac     => prt_params%allom_agb_frac(ipft), &
               allom_amode  => prt_params%allom_amode(ipft))

      branch_frac = param_derived%branch_frac(ipft)
      
      select case(allom_amode)
      case (1) !"salda")
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_salda(d,h,dhdd,p1,p2,p3,p4,wood_density,c2b,agb_frac,bagw,dbagwdd) 
      case (2) !"2par_pwr")
         ! Switch for woodland dbh->drc
         call d2bagw_2pwr(d,p1,p2,c2b,bagw,dbagwdd)
      case (3) !"chave14") 
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_chave2014(d,h,dhdd,p1,p2,wood_density,c2b,bagw,dbagwdd)
      case (4) ! 3par_pwr
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_3pwr(d,h,dhdd,p1,p2,p3,wood_density,c2b,bagw,dbagwdd)
      case (5) ! 3par_pwr_grass
         call h_allom(d,ipft,h,dhdd)
         call dh2bagw_3pwr_grass(d,h,dhdd,p1,p2,p3,c2b,bagw,dbagwdd)
      case DEFAULT
         write(fates_log(),*) 'An undefined AGB allometry was specified: ',allom_amode
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      
      ! Potentially reduce AGB based on crown damage (crown_reduction) and/or 
      ! phenology (elongf_stem).
      if(crowndamage > 1) then
         call GetCrownReduction(crowndamage, crown_reduction)
         bagw = elongf_stem * ( bagw - (bagw * branch_frac * crown_reduction) )
         if(present(dbagwdd))then
            dbagwdd = elongf_stem * ( dbagwdd - (dbagwdd * branch_frac * crown_reduction) )
         end if
      else
         bagw = elongf_stem * bagw
         if (present(dbagwdd)) then
            dbagwdd = elongf_stem * dbagwdd
         end if
      end if


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

    real(r8)               :: h       ! height
    real(r8)               :: dhdd    ! change in height wrt d

    associate( dbh_maxh    => prt_params%allom_dbh_maxheight(ipft), &
               rho         => prt_params%wood_density(ipft), &
               slatop      => prt_params%slatop(ipft),       &
               c2b         => prt_params%c2b(ipft),          &
               allom_lmode => prt_params%allom_lmode(ipft),  &
               p1          => prt_params%allom_d2bl1(ipft),  &
               p2          => prt_params%allom_d2bl2(ipft),  &
               p3          => prt_params%allom_d2bl3(ipft))
      
      select case(allom_lmode)
      case(1) !"salda")
         call d2blmax_salda(d,p1,p2,p3,rho,dbh_maxh,c2b,blmax,dblmaxdd)
      case(2) !"2par_pwr")
         call d2blmax_2pwr(d,p1,p2,c2b,blmax,dblmaxdd)
      case(3) ! dh2blmax_2pwr
         call dh2blmax_2pwr(d,p1,p2,dbh_maxh,c2b,blmax,dblmaxdd)
      case(4) ! dh2blmax_3pwr
         call h_allom(d,ipft,h,dhdd)
         call dh2blmax_3pwr(d,h,dhdd,p1,p2,p3,slatop,dbh_maxh,c2b,blmax,dblmaxdd)
      case (5) ! dh2blmax_3pwr_grass
         call h_allom(d,ipft,h,dhdd)
         call dh2blmax_3pwr_grass(d,h,dhdd,p1,p2,p3,dbh_maxh,c2b,blmax,dblmaxdd)
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
  
  subroutine carea_allom(dbh,nplant,site_spread,ipft,crowndamage,c_area,inverse)
     
     real(r8),intent(inout) :: dbh          ! plant diameter at breast (reference) height [cm]
     real(r8),intent(in)    :: site_spread  ! site level spread factor (crowdedness)
     real(r8),intent(in)    :: nplant       ! number of plants [1/ha]
     integer(i4),intent(in) :: ipft         ! PFT index
     integer(i4),intent(in) :: crowndamage  ! crown damage class [1: undamaged, >1: damaged]
     real(r8),intent(inout) :: c_area       ! crown area per cohort (m2)
     logical,optional,intent(in) :: inverse ! if true, calculate dbh from crown area 
                                            ! instead of crown area from dbh

     real(r8)               :: dbh_eff      ! Effective diameter (cm)
     real(r8)               :: height       ! height
     logical                :: do_inverse   ! local copy of the inverse argument
                                            ! defaults to false
     logical                :: capped_allom ! if we are using an allometry that caps
                                            ! crown area at height, we need to make
                                            ! special considerations
     
     associate( dbh_maxh    => prt_params%allom_dbh_maxheight(ipft), &
                allom_lmode => prt_params%allom_lmode(ipft),  &
                d2bl_p2     => prt_params%allom_d2bl2(ipft),  &
                d2bl_ediff  => prt_params%allom_blca_expnt_diff(ipft), &
                d2ca_min    => prt_params%allom_d2ca_coefficient_min(ipft), &
                d2ca_max    => prt_params%allom_d2ca_coefficient_max(ipft))
       
       if( .not. present(inverse) ) then 
          do_inverse = .false.
       else
          do_inverse = inverse
          if (do_inverse) then
             c_area = c_area / nplant
          endif
       endif

       select case(allom_lmode)
       case(1)
          dbh_eff = min(dbh,dbh_maxh)
          call carea_2pwr(dbh_eff,site_spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max, &
               crowndamage,c_area, do_inverse)
          capped_allom = .true.
       case(2)   ! "2par_pwr")
          call carea_2pwr(dbh,site_spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max, & 
               crowndamage, c_area, do_inverse)
          capped_allom = .false.
       case(3,5)
          dbh_eff = min(dbh,dbh_maxh)
          call carea_2pwr(dbh_eff,site_spread,d2bl_p2,d2bl_ediff,d2ca_min,d2ca_max, &
               crowndamage, c_area, do_inverse)
          capped_allom = .true.
       case (4)
          dbh_eff = min(dbh,dbh_maxh)
          call h_allom(dbh,ipft,height)
          call carea_3pwr(dbh_eff,height,ipft,dbh_maxh, site_spread,d2bl_p2, &
               d2bl_ediff, d2ca_min,d2ca_max,crowndamage, c_area, do_inverse)
          capped_allom = .true.
       case DEFAULT
          write(fates_log(),*) 'An undefined leaf allometry was specified: ', &
               allom_lmode
          write(fates_log(),*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       

       if (capped_allom .and. do_inverse) then
          if (dbh_eff .lt. dbh_maxh) then
             dbh = dbh_eff
          else
             ! In this situation, we are signaling to the
             ! calling routine that we we cannot calculate
             ! dbh from crown area, because we have already
             ! hit the area cap, and the two are not proportional
             ! anymore.  hopefully, the calling routine has an alternative
             dbh = fates_unset_r8
          endif
       endif

       c_area = c_area * nplant

    end associate
    return
  end subroutine carea_allom

  ! =====================================================================================
        
  subroutine bleaf(d,ipft,crowndamage,canopy_trim,elongf_leaf,bl,dbldd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target bleaf
    ! based on trimming. Because trimming
    ! is not allometry and rather an emergent property,
    ! this routine is not name-spaced with allom_
    ! -------------------------------------------------------------------------

    use DamageMainMod      , only : GetCrownReduction
    
    real(r8),intent(in)    :: d             ! plant diameter [cm]
    integer(i4),intent(in) :: ipft          ! PFT index
    integer(i4),intent(in) :: crowndamage   ! crown damage class [1: undamaged, >1: damaged]
    real(r8),intent(in)    :: canopy_trim   ! trimming function
    real(r8),intent(in)    :: elongf_leaf   ! Leaf elongation factor (phenology)
    real(r8),intent(out)   :: bl            ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbldd  ! change leaf bio per diameter [kgC/cm]
    
    real(r8) :: blmax
    real(r8) :: dblmaxdd
    real(r8) :: crown_reduction
    
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


    ! Potentially reduce leaf biomass based on crown damage (crown_reduction) and/or
    ! phenology (elongf_leaf).
    if ( crowndamage > 1 ) then

       call  GetCrownReduction(crowndamage, crown_reduction)
       bl = elongf_leaf * bl * (1.0_r8 - crown_reduction)
       if(present(dbldd))then
          dbldd = elongf_leaf * dblmaxdd * canopy_trim * (1.0_r8 - crown_reduction)
       end if
    else
       bl = elongf_leaf * bl
       if (present(dbldd)) then
          dbldd = elongf_leaf * dbldd
       end if
    end if

    return
  end subroutine bleaf
  
  ! =====================================================================================

  subroutine storage_fraction_of_target(c_store_target, c_store, frac)

    !--------------------------------------------------------------------------------
    !    This subroutine returns the ratio between the storage pool and the target 
    ! storage.  This subroutine is used both the carbon starvation mortality scheme 
    ! and the optional respiration throttling. We impose checks so it cannot go negative
    ! due to truncation errors, but this function can return values greater than 1.
    ! 
    !    Fractions exceeding do not impact the default linear carbon starvation model
    ! (mort_cstarvation_model=2), because mortality becomes zero, but they allow carbon 
    ! starvation mortality rates to continue decaying when the exponential carbon 
    ! starvation model is used (mort_cstarvation_model=2).
    ! 
    !    Fraction values above 1 do not impact lowstorage_maintresp_reduction either,
    ! as that routine imposes no reduction once the fraction exceeds 1.
    !--------------------------------------------------------------------------------

    real(r8),intent(in)    :: c_store_target  ! target storage carbon [kg]
    real(r8),intent(in)    :: c_store         ! storage carbon [kg]
    real(r8),intent(out)   :: frac

    frac = max(0._r8, c_store / max( c_store_target, nearzero) )

  end subroutine storage_fraction_of_target

  ! =====================================================================================

  real(r8) function tree_lai( leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top)

    ! -----------------------------------------------------------------------------------
    ! LAI of individual trees is a function of the total leaf area and the total 
    ! canopy area.   
    ! ----------------------------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: leaf_c                    ! plant leaf carbon [kg]
    integer, intent(in)  :: pft                       ! Plant Functional Type index
    real(r8), intent(in) :: c_area                    ! areal extent of canopy (m2)
    real(r8), intent(in) :: nplant                    ! number of individuals in cohort per ha
    integer, intent(in)  :: cl                        ! canopy layer index
    real(r8), intent(in) :: canopy_lai(:)             ! total leaf area index of 
                                                      ! each canopy layer
    real(r8), intent(in) :: vcmax25top                ! maximum carboxylation rate at canopy
                                                      ! top, ref 25C

    ! !LOCAL VARIABLES:
    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC
    real(r8) :: canopy_lai_above   ! total LAI of canopy layer overlying this tree
    real(r8) :: vai_per_lai        ! ratio of vegetation area index (ie. sai+lai) 
                                   ! to lai for individual tree
    real(r8) :: kn                 ! coefficient for exponential decay of 1/sla and 
                                   ! vcmax with canopy depth
    real(r8) :: sla_max            ! Observational constraint on how large sla 
                                   ! (m2/gC) can become
    real(r8) :: leafc_slamax       ! Leafc_per_unitarea at which sla_max is reached
    real(r8) :: clim               ! Upper limit for leafc_per_unitarea in exponential 
                                   ! tree_lai function
    !----------------------------------------------------------------------

    if( leaf_c  <  -1.1_r8*calloc_abs_error .or. pft  ==  0 ) then
       write(fates_log(),*) 'negative leaf carbon in LAI calculation?'
       write(fates_log(),*) 'or.. pft was zero?'
       write(fates_log(),*) 'problem in treelai',leaf_c,pft
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    slat = g_per_kg * prt_params%slatop(pft) ! m2/g to m2/kg
    leafc_per_unitarea = leaf_c/(c_area/nplant) !KgC/m2
    
    if(leafc_per_unitarea > 0.0_r8)then


       if (cl==1) then ! if in we are in the canopy (top) layer)
          canopy_lai_above = 0._r8
       else
          canopy_lai_above = sum(canopy_lai(1:cl-1))
       end if

       ! Coefficient for exponential decay of 1/sla with canopy depth:
       kn = decay_coeff_vcmax(vcmax25top, &
                              prt_params%leafn_vert_scaler_coeff1(pft), &
                              prt_params%leafn_vert_scaler_coeff2(pft))
       
       ! take PFT-level maximum SLA value, even if under a thick canopy (which has units of m2/gC),
       ! and put into units of m2/kgC
       sla_max = g_per_kg*prt_params%slamax(pft)
       ! Leafc_per_unitarea at which sla_max is reached due to exponential sla profile in canopy:
       leafc_slamax = (slat - sla_max * exp(-1.0_r8 * kn * canopy_lai_above)) / &
            (-1.0_r8 * kn * slat * sla_max)
       if(leafc_slamax < 0.0_r8)then
          leafc_slamax = 0.0_r8
       endif

       ! Calculate tree_lai (m2 leaf area /m2 ground) = unitless LAI
       !----------------------------------------------------------------------
       ! If leafc_per_unitarea is less than leafc_slamax,
       ! sla with depth in the canopy will not exceed sla_max.
       ! In this case, we can use an exponential profile for sla throughout the entire canopy.
       ! The exponential profile for sla is given by:
       ! sla(at a given canopy depth) = slat / exp(-kn (canopy_lai_above + tree_lai)
       ! 
       ! We can solve for tree_lai using the above function for the sla profile and first setting 
       ! leafc_per_unitarea = integral of e^(-kn(x + canopy_lai_above)) / slatop
       ! over x = 0 to tree_lai
       ! Then, rearranging the equation to solve for tree_lai.

       if (leafc_per_unitarea <= leafc_slamax)then
          tree_lai = (log(exp(-1.0_r8 * kn * canopy_lai_above) - &
               kn * slat * leafc_per_unitarea) + &
               (kn * canopy_lai_above)) / (-1.0_r8 * kn)

          ! If leafc_per_unitarea becomes too large, tree_lai becomes an imaginary number 
          ! (because the tree_lai equation requires us to take the natural log of something >0)
          ! Thus, we include the following error message in case leafc_per_unitarea becomes too large.
          clim = (exp(-1.0_r8 * kn * canopy_lai_above)) / (kn * slat)
          if (leafc_per_unitarea >= clim) then
             write(fates_log(),*) 'too much leafc_per_unitarea' , leafc_per_unitarea, clim, pft, canopy_lai_above
             write(fates_log(),*) 'Aborting'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

          ! When leafc_per_unitarea is greater than leafc_slamax, 
          ! tree_lai could become so great that the sla profile surpasses sla_max at depth.
          ! In this case, we use the exponential profile to calculate tree_lai until
          ! we reach the maximum allowed value for sla (sla_max).
          ! Then, calculate the remaining tree_lai using a linear function of sla_max and the remaining leafc.
       
       else if(leafc_per_unitarea > leafc_slamax)then
          
          ! Add exponential and linear portions of tree_lai
          ! Exponential term for leafc = leafc_slamax; 
          ! Linear term (static sla = sla_max) for portion of leafc > leafc_slamax
          tree_lai = ((log(exp(-1.0_r8 * kn * canopy_lai_above) - &
               kn * slat * leafc_slamax) + &
               (kn * canopy_lai_above)) / (-1.0_r8 * kn)) + &
               (leafc_per_unitarea - leafc_slamax) * sla_max

          ! if leafc_slamax becomes too large, tree_lai_exp becomes an imaginary number 
          ! (because the tree_lai equation requires us to take the natural log of something >0)
          ! Thus, we include the following error message in case leafc_slamax becomes too large.
          clim = (exp(-1.0_r8 * kn * canopy_lai_above)) / (kn * slat)
          if(leafc_slamax >= clim)then
             write(fates_log(),*) 'too much leafc_slamax' , &
                  leafc_per_unitarea, leafc_slamax, clim, pft, canopy_lai_above
             write(fates_log(),*) 'Aborting'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
       end if ! (leafc_per_unitarea  > leafc_slamax)
    else
       tree_lai = 0.0_r8
    endif ! (leafc_per_unitarea > 0.0_r8)

    return
  end function tree_lai

  ! ============================================================================

  real(r8) function tree_sai(pft, dbh, crowndamage, canopy_trim, elongf_stem, c_area, nplant, &
                             cl, canopy_lai, treelai, vcmax25top, call_id )

    ! ============================================================================
    !  SAI of individual trees is a function of the LAI of individual trees
    ! ============================================================================

    integer, intent(in)  :: pft
    real(r8), intent(in) :: dbh
    integer, intent(in)  :: crowndamage
    real(r8), intent(in) :: canopy_trim        ! trimming function (0-1)
    real(r8), intent(in) :: elongf_stem        ! Elongation factor for stems.
    real(r8), intent(in) :: c_area             ! crown area (m2)
    real(r8), intent(in) :: nplant             ! number of plants
    integer, intent(in)  :: cl                 ! canopy layer index
    real(r8), intent(in) :: canopy_lai(:)      ! total leaf area index of 
                                               ! each canopy layer
    real(r8), intent(in) :: treelai            ! tree LAI for checking purposes only
    real(r8), intent(in) :: vcmax25top         ! maximum carboxylation rate at top of crown
    integer,intent(in)   :: call_id            ! flag specifying where this is called
                                               ! from
    real(r8)             :: h
    real(r8)             :: target_lai
    real(r8)             :: target_bleaf

    ! Assume fully flushed leaves, so stem area index is independent on leaf phenology.
    ! SAI can be downscaled by stem phenology (typically applied to grasses only).
    call bleaf(dbh, pft, crowndamage, canopy_trim, 1.0_r8, target_bleaf)

    target_lai = tree_lai(target_bleaf, pft, c_area, nplant, cl,&
         canopy_lai, vcmax25top) 

    tree_sai   =  elongf_stem * prt_params%allom_sai_scaler(pft) * target_lai

    return
  end function tree_sai

  ! ============================================================================
  
  subroutine tree_lai_sai(leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top, &
                          dbh, crowndamage, canopy_trim, elongf_stem, call_id, & 
                          treelai, treesai)

    ! This is the public wrapper for calling plant lai and sai
    ! The purpose of this wrapper is to ensure that they are called together,
    ! and in the correct order, and that the capping is applied after
    ! the base calculations are performed.
    
    real(r8), intent(in) :: leaf_c                    ! plant leaf carbon [kg]
    integer, intent(in)  :: pft                       ! Plant Functional Type index
    real(r8), intent(in) :: c_area                    ! areal extent of canopy (m2)
    real(r8), intent(in) :: nplant                    ! number of individuals in cohort per ha
    integer, intent(in)  :: cl                        ! canopy layer index
    real(r8), intent(in) :: canopy_lai(:)             ! total leaf area index of 
                                                      ! each canopy layer
    real(r8), intent(in) :: vcmax25top                ! maximum carboxylation rate at canopy
    real(r8), intent(in) :: dbh
    integer, intent(in)  :: crowndamage
    real(r8), intent(in) :: canopy_trim        ! trimming function (0-1)
    real(r8), intent(in) :: elongf_stem        ! Elongation factor for stems.
    integer,intent(in)   :: call_id            ! flag specifying where this is called from
    ! from
    
    real(r8), intent(out) :: treelai         ! plant LAI [m2 leaf area/m2 crown area]
    real(r8), intent(out) :: treesai         ! plant SAI [m2 stem area/m2 crown area]
   
    
    treelai = tree_lai( leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top)
    
    treesai = tree_sai( pft, dbh, crowndamage, canopy_trim, elongf_stem, c_area, nplant, &
                             cl, canopy_lai, treelai, vcmax25top, call_id )

    ! Don't allow lai+sai to exceed the vertical discretization bounds
    if( (treelai + treesai) > (sum(dinc_vai)) )then
       treelai = sum(dinc_vai) * (1._r8 - prt_params%allom_sai_scaler(pft)) - nearzero
       treesai = sum(dinc_vai) * prt_params%allom_sai_scaler(pft) - nearzero
    end if
    
    
    return
  end subroutine tree_lai_sai

  
! =====================================================================================
 
  real(r8) function leafc_from_treelai( treelai, treesai, pft, c_area, nplant, cl, vcmax25top)
 
    ! -----------------------------------------------------------------------------------
    ! Calculates the amount of leaf carbon which is needed to generate a given treelai. 
    ! iss the inverse of the 'tree_lai function. 
    ! ----------------------------------------------------------------------------------
 
    ! !ARGUMENTS
    real(r8), intent(in) :: treelai                    ! desired tree lai m2/m2
    real(r8), intent(in) :: treesai                    ! desired tree lai m2/m2
    integer, intent(in)  :: pft                       ! Plant Functional Type index
    real(r8), intent(in) :: c_area                    ! areal extent of canopy (m2)
    real(r8), intent(in) :: nplant                    ! number of individuals in cohort per ha
    integer, intent(in)  :: cl                        ! canopy layer index
    real(r8), intent(in) :: vcmax25top                ! maximum carboxylation rate at canopy
                                                      ! top, ref 25C
 
    ! !LOCAL VARIABLES:
    real(r8) :: leaf_c                    ! plant leaf carbon [kg]
    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC
    real(r8) :: vai_per_lai        ! ratio of vegetation area index (ie. sai+lai)
                                   ! to lai for individual tree
    real(r8) :: kn                 ! coefficient for exponential decay of 1/sla and
                                   ! vcmax with canopy depth
    real(r8) :: sla_max            ! Observational constraint on how large sla
                                   ! (m2/gC) can become
    real(r8) :: leafc_slamax       ! Leafc_per_unitarea at which sla_max is reached
    real(r8) :: clim               ! Upper limit for leafc_per_unitarea in exponential
                                   ! tree_lai function
    real(r8) :: tree_lai_at_slamax ! lai at which we reach the maximum sla value.
    real(r8) :: leafc_linear_phase ! amount of leaf carbon needed to get to the target treelai
                                   ! when the slamax value has been reached (i.e. deep layers with unchanging sla)

    !----------------------------------------------------------------------
 
    if( treelai  < 0._r8.or. pft  ==  0 ) then
       write(fates_log(),*) 'negative tree lai in leafc_from_treelai?'
       write(fates_log(),*) 'or.. pft was zero?'
       write(fates_log(),*) 'problem in leafc_from_treelai',treelai,pft
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
 
    if(cl>1)then
      write(fates_log(),*) 'in sub-canopy layer in leafc_from_treelai'
      write(fates_log(),*) 'this is not set up to work for lower canopy layers.'
      write(fates_log(),*) 'problem in leafc_from_treelai',cl,pft
      call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    if( (treelai + treesai) > (sum(dinc_vai)) )then
       write(fates_log(),*) 'SP tree lai cannot exceed sum of dinc_vai',treelai,treesai,sum(dinc_vai)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    ! convert PFT-level canopy top and maximum SLA values and convert from m2/gC to m2/kgC
    slat = g_per_kg * prt_params%slatop(pft)
    sla_max = g_per_kg * prt_params%slamax(pft)

    ! Coefficient for exponential decay of 1/sla with canopy depth:
    kn = decay_coeff_vcmax(vcmax25top, &
                           prt_params%leafn_vert_scaler_coeff1(pft), &
                           prt_params%leafn_vert_scaler_coeff2(pft))

    if(treelai > 0.0_r8)then 
       ! Leafc_per_unitarea at which sla_max is reached due to exponential sla profile in canopy:
       leafc_slamax = max(0.0_r8,(slat - sla_max) / (-1.0_r8 * kn * slat * sla_max))

       ! treelai at which we reach maximum sla.
       tree_lai_at_slamax = (log( 1.0_r8- kn * slat * leafc_slamax)) / (-1.0_r8 * kn)

       if(treelai < tree_lai_at_slamax)then
         ! Inversion of the exponential phase calculation of treelai for a given leafc_per_unitarea
         leafc_per_unitarea = (1.0_r8-exp(treelai*(-1.0_r8 * kn)))/(kn*slat)
       else ! we exceed the maxumum sla
 
        ! Add exponential and linear portions of tree_lai
        ! Exponential term for leafc = leafc_slamax;
         leafc_linear_phase = (treelai-tree_lai_at_slamax)/sla_max
         leafc_per_unitarea = leafc_slamax + leafc_linear_phase
       end if
       leafc_from_treelai = leafc_per_unitarea*(c_area/nplant)
    else
       leafc_from_treelai = 0.0_r8 
    endif ! (leafc_per_unitarea > 0.0_r8)
 
    return
  end function leafc_from_treelai
 
  ! =====================================================================================






  ! ============================================================================
  ! Generic sapwood biomass interface
  ! ============================================================================

  subroutine bsap_allom(d,ipft,crowndamage,canopy_trim,elongf_stem, sapw_area,bsap,dbsapdd)

    use DamageMainMod , only : GetCrownReduction
    use FatesParameterDerivedMod, only : param_derived
    
    real(r8),intent(in)           :: d           ! plant diameter [cm]
    integer(i4),intent(in)        :: ipft        ! PFT index
    integer(i4),intent(in)        :: crowndamage ! Crown damage class [1: undamaged, >1: damaged]
    real(r8),intent(in)           :: canopy_trim
    real(r8),intent(in)           :: elongf_stem ! Elongation factor for stems (phenology)
    real(r8),intent(out)          :: sapw_area   ! cross section area of
                                                 ! plant sapwood at reference [m2]
    real(r8),intent(out)          :: bsap        ! sapwood biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd     ! change in sapwood biomass
                                                 ! per d [kgC/cm]

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

    real(r8) :: crown_reduction  ! amount that crown is damage by
    real(r8) :: agb_frac         ! aboveground biomass fraction
    real(r8) :: branch_frac      ! fraction of aboveground woody biomass in branches
    
    ! Constrain sapwood so that its above ground portion be no larger than 
    ! X% of total woody/fibrous (ie non leaf/fineroot) tissues
    real(r8),parameter :: max_frac = 0.95_r8 

    agb_frac = prt_params%allom_agb_frac(ipft)
    branch_frac = param_derived%branch_frac(ipft)
      
    
    select case(prt_params%allom_smode(ipft))
       ! ---------------------------------------------------------------------
       ! Currently only one sapwood allometry model. the slope
       ! of the la:sa to diameter line is zero.
       ! ---------------------------------------------------------------------
    case(1) ! linearly related to leaf area based on target leaf biomass
            ! and slatop (no provisions for slamax)

       !  We assume fully flushed leaves, so sapwood biomass is independent of leaf phenology
       ! (but could be modulated by stem phenology).
       call h_allom(d,ipft,h,dhdd)
       call bleaf(d,ipft,1,canopy_trim,1.0_r8,bl,dbldd)
       call bsap_ltarg_slatop(d,h,dhdd,bl,dbldd,ipft,sapw_area,bsap,dbsapdd)

       ! if trees are damaged reduce bsap by percent crown loss *
       ! fraction of biomass that would be in branches (pft specific)
       if(crowndamage > 1)then

          call GetCrownReduction(crowndamage, crown_reduction)
          bsap = elongf_stem * ( bsap - (bsap * agb_frac *  branch_frac * crown_reduction) )
          if(present(dbsapdd))then
             dbsapdd = elongf_stem * &
                       ( dbsapdd - (dbsapdd * agb_frac * branch_frac * crown_reduction) )
          end if
       else
          bsap = elongf_stem * bsap
          if (present(dbsapdd)) then
             dbsapdd = elongf_stem * dbsapdd
          end if
       end if
       
       
       ! Perform a capping/check on total woody biomass
       call bagw_allom(d,ipft,crowndamage, elongf_stem, bagw,dbagwdd)
       call bbgw_allom(d,ipft, elongf_stem,bbgw,dbbgwdd)
       
       ! Force sapwood to be less than a maximum fraction of total biomass
       ! We omit the sapwood area from this calculation
       ! (this comes into play typically in very small plants)
       bsap_cap = max_frac*(bagw+bbgw)

       if(bsap>bsap_cap) then
          bsap     = bsap_cap
          if(present(dbsapdd))then
             dbsapdd = max_frac*(dbagwdd+dbbgwdd)
          end if
       end if

    case(2) ! this is a 'sapwood' function specifically for grass PFT that do not produce
            ! dead woody biomass. So bsap = bagw. Might remove the bsap and bdead for grass
            ! in the future as there is no need to distinguish the two for grass above- and belowground biomass

       call SapwoodAreaGrass(d,sapw_area)
       call bagw_allom(d,ipft,crowndamage,elongf_stem,bagw,dbagwdd)
       call bbgw_allom(d,ipft, elongf_stem,bbgw,dbbgwdd)
       
       bsap = bagw + bbgw

       ! This is a grass-only functionnal type, no need to run crown-damage effects

       bsap = elongf_stem * bsap
       if (present(dbsapdd))then
          dbsapdd = elongf_stem * (dbagwdd + dbbgwdd)
       end if
       
       if(present(dbsapdd))then
          dbsapdd = dbagwdd + dbbgwdd
       end if

    case DEFAULT
       write(fates_log(),*) 'An undefined sapwood allometry was specified: ', &
            prt_params%allom_smode(ipft)
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

  subroutine bbgw_allom(d,ipft,elongf_stem,bbgw,dbbgwdd)

    real(r8),intent(in)           :: d           ! plant diameter [cm]
    integer(i4),intent(in)        :: ipft        ! PFT index
    real(r8),intent(in)           :: elongf_stem ! Elongation factor for stems (phenology)
    real(r8),intent(out)          :: bbgw        ! below ground woody biomass [kgC]
    real(r8),intent(out),optional :: dbbgwdd     ! change bbgw  per diam [kgC/cm]
    
    real(r8)    :: bagw       ! above ground biomass [kgC]
    real(r8)    :: dbagwdd    ! change in agb per diameter [kgC/cm]
    
    select case(prt_params%allom_cmode(ipft))
    case(1) !"constant")
       ! bbgw not affected by damage so use target allometry no damage. But note that bbgw
       ! is affected by stem phenology (typically applied only to grasses). We do not need
       ! to account for stem phenology in bbgw_const because bbgw will be proportional to
       ! bagw, and bagw is downscaled due to stem phenology.
       call bagw_allom(d,ipft,1, elongf_stem, bagw,dbagwdd)
       call bbgw_const(d,bagw,dbagwdd,ipft,bbgw,dbbgwdd)
    case DEFAULT
       write(fates_log(),*) 'An undefined coarse root allometry was specified: ', &
             prt_params%allom_cmode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select
    return
  end subroutine bbgw_allom
 
  ! ============================================================================
  ! Fine root biomass allometry wrapper
  ! ============================================================================
  
  subroutine bfineroot(d,ipft,canopy_trim,l2fr,elongf_fnrt,bfr,dbfrdd)
    
    ! -------------------------------------------------------------------------
    ! This subroutine calculates the actual target fineroot biomass
    ! based on functions that may or may not have prognostic properties. 
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d             ! plant diameter [cm]
    integer(i4),intent(in) :: ipft          ! PFT index
    real(r8),intent(in)    :: canopy_trim   ! trimming function
    real(r8),intent(in)    :: l2fr          ! leaf to fineroot scaler
                                            ! this is either a PFT parameter
                                            ! constant (when no nutrient model)
                                            ! or dynamic (with nutrient model)
    real(r8),intent(in)    :: elongf_fnrt   ! Elongation factor for fine roots
    real(r8),intent(out)   :: bfr           ! fine root biomass [kgC]
    real(r8),intent(out),optional :: dbfrdd ! change leaf bio per diameter [kgC/cm]
    
    real(r8) :: blmax      ! maximum leaf biomss per allometry
    real(r8) :: dblmaxdd
    real(r8) :: bfrmax
    real(r8) :: dbfrmaxdd
    real(r8) :: slascaler
    
    select case(prt_params%allom_fmode(ipft))
    case(1) ! "constant proportionality with TRIMMED target bleaf"
       
       call blmax_allom(d,ipft,blmax,dblmaxdd)

       bfr = blmax*l2fr*canopy_trim
       
       if(present(dbfrdd))then
          dbfrdd = dblmaxdd*l2fr * canopy_trim
          
       end if
    case(2) ! "constant proportionality with UNTRIMMED target bleaf"
       
       call blmax_allom(d,ipft,blmax,dblmaxdd)

       bfr = blmax*l2fr
       if(present(dbfrdd))then
          dbfrdd = dblmaxdd*l2fr
       end if

    case DEFAULT 
       write(fates_log(),*) 'An undefined fine root allometry was specified: ', &
            prt_params%allom_fmode(ipft)
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select


    ! Reduce fine-root biomass due to phenology.
    bfr = elongf_fnrt * bfr
    if (present(dbfrdd)) then
       dbfrdd = elongf_fnrt * dbfrdd
    end if


    return
  end subroutine bfineroot


  ! ============================================================================
  ! Storage biomass interface
  ! ============================================================================
  
  subroutine bstore_allom(d,ipft,crowndamage, canopy_trim,bstore,dbstoredd)

     real(r8),intent(in)           :: d            ! plant diameter [cm]
     integer(i4),intent(in)        :: ipft         ! PFT index
     integer(i4),intent(in)        :: crowndamage  ! Crowndamage class [1: undamaged, >1: damaged]
     real(r8),intent(in)           :: canopy_trim  ! Crown trimming function [0-1]
     real(r8),intent(out)          :: bstore       ! allometric target storage [kgC]
     real(r8),intent(out),optional :: dbstoredd    ! change storage per cm [kgC/cm]
     
     real(r8) :: bl          ! Allometric target leaf biomass
     real(r8) :: dbldd       ! Allometric target change in leaf biomass per cm
     real(r8) :: blmax       ! Allometric target leaf biomass (UNTRIMMED)
     real(r8) :: dblmaxdd    ! Allometric target change in leaf biomass per cm (UNTRIMMED)
    
     
     associate( allom_stmode => prt_params%allom_stmode(ipft), &
                cushion      => prt_params%cushion(ipft) )

       select case(allom_stmode)
       case(1) ! Storage is constant proportionality of trimmed maximum leaf
          ! biomass (ie cushion * bleaf), and thus leaf phenology is ignored.
          call bleaf(d,ipft, crowndamage, canopy_trim, 1.0_r8, bl, dbldd)
          call bstore_blcushion(d,bl,dbldd,cushion,ipft,bstore,dbstoredd)

       case(2) ! Storage is constant proportionality of untrimmed maximum leaf
          ! biomass (ie cushion * bleaf_max)
          call blmax_allom(d,ipft,blmax,dblmaxdd)
          call bstore_blcushion(d,blmax,dblmaxdd,cushion,ipft,bstore,dbstoredd)

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

    
    associate( agb_fraction => prt_params%allom_agb_frac(ipft))

      select case(prt_params%allom_amode(ipft))
      case(1) ! Saldariagga mass allometry originally calculated bdead directly.
              ! we assume proportionality between bdead and bagw
       
         bdead = bagw/agb_fraction 
         if(present(dbagwdd) .and. present(dbdeaddd))then
            dbdeaddd = dbagwdd/agb_fraction
         end if
         
      case(2,3,4,5)
         
         bdead = bagw + bbgw - bsap
         if(present(dbagwdd) .and. present(dbbgwdd) .and. &
            present(dbdeaddd) .and. present(dbsapdd) )then
            dbdeaddd = dbagwdd+dbbgwdd-dbsapdd
         end if
         
      case DEFAULT
         
         write(fates_log(),*) 'An undefined AGB allometry was specified: ',&
                              prt_params%allom_amode(ipft)
         write(fates_log(),*) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
       
      end select
    end associate
    return
  end subroutine bdead_allom


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

    associate( agb_fraction => prt_params%allom_agb_frac(ipft) )
      
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

  subroutine bsap_ltarg_slatop(d,h,dhdd,bleaf,dbleafdd,ipft, &
        sapw_area,bsap,dbsapdd)
    
    ! -------------------------------------------------------------------------
    ! Calculate sapwood carbon based on its leaf area per sapwood area
    ! proportionality with the plant's target leaf area.
    ! of plant size, see Calvo-Alvarado and Bradley Christoferson.
    ! 
    ! Important note 1: This is above and below-ground sapwood
    ! Important note 2: Since we need continuous calculation of
    !                   sapwood dependent on plant size, we cannot
    !                   use actual leaf area (which is canopy dependent).
    !                   So, this method estimates a leaf area that is
    !                   based only on the specific leaf area (SLA) of 
    !                   the canopy top.
    !
    ! -------------------------------------------------------------------------
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    real(r8),intent(in)    :: bleaf     ! plant leaf target biomass [kgC]
    real(r8),intent(in)    :: dbleafdd  ! change in blmax per diam [kgC/cm]
    integer(i4),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: sapw_area ! area of sapwood crosssection at 
                                        ! reference height [m2]
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dbsapdd   ! change leaf bio per diameter [kgC/cm]
    
    real(r8)               :: la_per_sa  ! effective leaf area per sapwood area
                                         ! [m2/cm]
    real(r8)               :: term1      ! complex term for solving derivative
    real(r8)               :: dterm1_dh  ! deriv of term1 wrt height
    real(r8)               :: dterm1_dd  ! deriv of term1 wrt diameter
    real(r8)               :: hbl2bsap   ! sapwood biomass per lineal height
    
    
    associate ( la_per_sa_int => prt_params%allom_la_per_sa_int(ipft), &
                la_per_sa_slp => prt_params%allom_la_per_sa_slp(ipft), &
                slatop        => prt_params%slatop(ipft), &
                wood_density  => prt_params%wood_density(ipft), &
                c2b           => prt_params%c2b(ipft), & 
                agb_fraction  => prt_params%allom_agb_frac(ipft) )


      ! Calculate sapwood biomass per linear height and kgC of leaf [m-1]
      ! Units: 
      ! Note: wood_density is in units of specific gravity, which is also 
      !       Mg / m3  (megagrams, ie 1000 kg / m3)
      ! 1 /la_per_sa * slatop*     gtokg    *   cm2tom2     / c2b     * mg2kg  * dens
      ! [cm2/m2]     * [m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3]
      !        ->[cm2/gC]
      !                  ->[cm2/kgC]
      !                                ->[m2/kgC]
      !                                             ->[m2/kg]
      !                                                       ->[m2/Mg]
      !                                                                  ->[/m]
      ! ------------------------------------------------------------------------

      ! This is a term that combines unit conversion and specific leaf
      ! area.  This term does not contain the proportionality
      ! between leaf area and sapwood cross-section. This is
      ! because their may be a height dependency, and will effect the 
      ! derivative wrt diameter.
      hbl2bsap   = slatop * g_per_kg * wood_density * kg_per_Megag / (c2b*cm2_per_m2 )

      ! Calculate area. Note that no c2b conversion here, because it is
      ! wood density that is in biomass units, SLA is in units [m2/gC].
      ! [m2]    = [m2/gC] * [kgC] * [gC/kgC] / ( [m2/cm2] * [cm2/m2])
      la_per_sa = la_per_sa_int + h*la_per_sa_slp
      sapw_area = slatop * bleaf * g_per_kg / (la_per_sa*cm2_per_m2 )
      

      ! Note the total depth of the plant is approximated by the 
      ! above ground fraction. This fraction is actually associated
      ! with biomass, but we use it here as well to help us assess
      ! how much sapwood is above and below ground.
      ! total_depth * agb_fraction = height 

      ! Integrate the mass per leaf biomass per depth of the plant
      ! Include above and below ground components
      ! [kgC] = [kgC/kgC/m]   * [kgC]     * [m]
      ! ------------------------------------------------------------------------

!      bsap =  hbl2bsap/(la_per_sa_int + h*la_per_sa_slp) * (h/agb_fraction) * bleaf

      ! "term1" combines two height dependent functions. The numerator is
      ! how sapwood volume scales in the vertical direction.  The denominator
      ! is the leaf_area per sapwood area ratio [m2/cm2], which is height dependent
      ! (for non-zero slope parameters)

      term1 = h/(la_per_sa_int + h*la_per_sa_slp)
      bsap  = (hbl2bsap/agb_fraction) * term1 * bleaf 


      ! dbldmaxdd is deriv of blmax wrt dbh (use directives to check oop)
      ! dhdd is deriv of height wrt dbh (use directives to check oop)
      if(present(dbsapdd))then
         dterm1_dh = la_per_sa_int  / (la_per_sa_int + la_per_sa_slp*h)**2.0_r8
         dterm1_dd = dterm1_dh * dhdd
         dbsapdd  = hbl2bsap/agb_fraction * (bleaf*dterm1_dd + term1 *dbleafdd)
      end if
      
    end associate
    return
  end subroutine bsap_ltarg_slatop



  ! ============================================================================
  ! Area of sap wood cross-section specifically for grass PFT
  ! ============================================================================

  subroutine SapwoodAreaGrass(d,sapw_area)

     !---------------------------------------------------------------------------
     ! This function calculates sapwood cross-sectional area specifically for grass
     ! PFT using basal diameter (cm) of the entire plant as size reference,
     ! assume sapwood area of the entire plant as the sum of the cross-sectional area
     ! of each grass tiller
     ! such that water transport through sapwood can be seen as a collective behavior
     ! of all tillers
     ! No reference. Might update this to more theoretical-based approach once there
     ! is empirical evidence
     !----------------
     ! Input arguments
     !----------------
     ! d                -- basal diameter               [cm]

     !----------------
     ! Output variables
     !----------------
     ! sapw_area        -- sapwood cross-sectional area   [m2]

     !---Arguments
     real(r8), intent(in)              :: d             ! plant basal diameter               [    cm]
     real(r8), intent(out)             :: sapw_area     ! sapwood cross-sectional area       [    m2]

     ! Calculate sapwood cross-sectional area assuming sapwood geometry as a
     ! cylinder and basal diameter is the diameter of the cylinder
     sapw_area = (pi_const * ((d / 2.0_r8)**2.0_r8)) / cm2_per_m2

     return

  end subroutine SapwoodAreaGrass

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
    ! ======================================================================
    
    ! p1 and p2 represent the parameters that govern total beaf dry biomass, 
    ! and the output argument blmax is the leaf carbon only
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    real(r8),intent(in)  :: p1        ! parameter 1  (slope)
    real(r8),intent(in)  :: p2        ! parameter 2  (curvature, exponent)
    real(r8),intent(in)  :: c2b       ! carbon to biomass multiplier (~2)
    
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    
    blmax    = (p1*d**p2) / c2b
    
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
    real(r8),intent(in)  :: p1        ! parameter 1 (slope)
    real(r8),intent(in)  :: p2        ! parameter 2 (curvature, exponent)
    real(r8),intent(in)  :: c2b       ! carbon 2 biomass multiplier
    real(r8),intent(in)  :: dbh_maxh  ! dbh at maximum height
    
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(out),optional :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    
    ! reproduce Saldarriaga:
    ! blmax = p1 * dbh_maxh**p2 * rho**p3
    !       = 0.07 * dbh_maxh**p2 * 0.7*0.55 = (p1 + p2* dbh**p3) / c2b
    !       p1 = 0
    !       p2 = (0.07 * 0.7^0.55)*2 = 0.11506201034678605

    blmax    = (p1*min(d,dbh_maxh)**p2)/c2b
    
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


  ! ===========================================================================


  subroutine dh2blmax_3pwr(d,h,dhdd,p1,p2,p3,slatop,dbh_maxh,c2b,blmax,dblmaxdd)
     !--------------------------------------------------------------------------
     !
     !     This function calculates the maximum leaf biomass from reference 
     ! diameter, plant height and top-of-the-canopy (fully sunlit) SLA. This 
     ! functional form is similar to Lescure et al. (1983) and Longo et al.
     ! (2020), except that it uses SLA as an additional scaler for the 
     ! allometric equation that can have a different exponent from 
     ! (DBH^2 * Height).
     !
     ! -----------------
     !  References
     ! -----------------
     ! Lescure JP, Puig H, Riera B, Leclerc D, Beekman A , Beneteau A. 1983.
     !    La phytomasse epigee d'une foret dense en Guiane Francaise
     !    Acta Oecol.-Oec. Gen. 4: 237-251.
     !    URL http://www.documentation.ird.fr/hor/fdi:010005089
     !
     ! Longo M, Saatchi SS, Keller M, Bowman KW, Ferraz A, Moorcroft PR,
     !    Morton D, Bonal D, Brando P, Burban B et al. 2020. Impacts of
     !    degradation on water, energy, and carbon cycling of the Amazon 
     !    tropical forests. J. Geophys. Res.-Biogeosci. 125:
     !    e2020JG005677. doi:10.1029/2020JG005677.
     !
     ! -----------------
     !  Input arguments
     ! -----------------
     ! d        -- Diameter at breast height          [     cm]
     ! h        -- Total tree height                  [      m]
     ! dhdd     -- Height derivative with dbh         [   m/cm]
     ! p1       -- Parameter 1 (log-intercept)        [     --]
     ! p2       -- Parameter 2 (power, or log-slope)  [     --]
     ! p3       -- Parameter 3 (power, or log-slope)  [     --]
     ! slatop   -- Top-of-canopy specific leaf area   [  m2/gC]
     ! dbh_maxh -- DBH at maximum height              [     cm]
     ! c2b      -- Carbon to biomass multiplier ~ 2   [ kg/kgC]
     !
     ! ------------------
     !  Output arguments
     ! ------------------
     ! blmax    -- Plant leaf biomass                 [    kgC]
     ! dblmaxdd -- Plant leaf biomass derivative      [ kgC/cm]
     !
     ! ------------------
     !   Suggested first guess for parameters, based on Longo et al. (2020) and
     ! corrected to FATES units (first guess based on a very limited leaf area
     ! data set).
     ! ------------------
     ! p1 =  0.000468
     ! p2 =  0.641
     ! p3 = -1.000
     !--------------------------------------------------------------------------


     !--- Arguments
     real(r8), intent(in)            :: d        ! plant diameter                 [    cm]
     real(r8), intent(in)            :: h        ! plant height                   [     m]
     real(r8), intent(in)            :: dhdd     ! Height derivative wrt diameter [  m/cm]
     real(r8), intent(in)            :: p1       ! Log-intercept parameter        [     -]
     real(r8), intent(in)            :: p2       ! Log-slope parameter for size   [     -]
     real(r8), intent(in)            :: p3       ! Log-slope parameter for SLA    [     -]
     real(r8), intent(in)            :: slatop   ! Top canopy specific leaf area  [ m2/gC]
     real(r8), intent(in)            :: c2b      ! Carbon to biomass multiplier   [kg/kgC]
     real(r8), intent(in)            :: dbh_maxh ! dbh at maximum height          [    cm]
     real(r8), intent(out)           :: blmax    ! Leaf biomass                   [   kgC]
     real(r8), intent(out), optional :: dblmaxdd ! Leaf biomass derivative        [kgC/cm]
     !--- Local variables
     real(r8) :: duse
     !---~---



     !--- Cap DBH
     duse = min(d,dbh_maxh)
     !---~---


     !--- Find leaf biomass
     blmax = p1 * (duse*duse*h)**p2 * slatop**p3 / c2b
     !---~---


     !---~---
     !   Compute the leaf biomass derivative with DBH if needed.
     !---~---
     if (present(dblmaxdd)) then
        if (d >= dbh_maxh) then
           !---~---
           !   Leaf area is capped at the maximum DBH. This may be removed in the
           ! future.
           !---~---
           dblmaxdd = 0._r8
           !---~---
        else
           !---~---
           !   Find the leaf biomass derivative, noting that height is actually
           ! a function of DBH.
           !---~---
           dblmaxdd = p2 * blmax * ( 2._r8 / duse + dhdd / h )
           !---~---
        end if
     end if
     !---~---

     return
  end subroutine dh2blmax_3pwr



  ! =========================================================================

  subroutine dh2blmax_3pwr_grass(d,h,dhdd,p1,p2,p3,dbh_maxh,c2b,blmax,dblmaxdd)
    !------------------------------------------------------------------------
    !
    ! This function calculates the maximum leaf biomass using diameter (basal
    ! diameter for grass) and plant height based on grass leaf allometry developed
    ! in Gao et al. 2024
    !
    !-------------------
    ! References
    !-------------------
    ! Gao X., Koven C., and Kueppers L. 2024. Allometric relationships and trade-offs
    ! in 11 common Mediterranean-climate grasses. Ecological Applications.
    ! https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2976

    !------------------
    ! Input arguments
    !------------------
    ! d          -- Basal diameter for grass or any herbaceous plants             [      cm]
    ! h          -- Plant height                                                  [       m]
    ! dhdd       -- Height derivative with dbh                                    [    m/cm]
    ! p1         -- Parameter 1 (log-intercept)                                   [      --]
    ! p2         -- Parameter 2 (log-slope associated with d)                     [      --]
    ! p3         -- Parameter 3 (log-slope associated with h)                     [      --]
    ! dbh_maxh   -- DBH at maximum height                                         [      cm]
    ! c2b        -- Carbon to biomass multiplier ~ 2                              [  kg/kgC]
    !
    !-----------------
    ! Output arguments
    !-----------------
    ! blmax      -- Leaf biomass                                                  [      kgC]
    ! dblmaxdd   -- Leaf biomass derivative                                       [   kgC/cm]
    !
    !-----------------
    ! Default parameters have been updated for three FATES grass PFTs according to Gao et al. 2024
    !---------------------------------------------------------------------------------------------
    
    !-----Arguments
    real(r8), intent(in)              :: d         ! plant diameter                   [    cm]
    real(r8), intent(in)              :: h         ! plant height                     [     m]
    real(r8), intent(in)              :: dhdd      ! height derivative                [  m/cm]
    real(r8), intent(in)              :: p1        ! log-intercept parameter          [     -]
    real(r8), intent(in)              :: p2        ! log-slope associated with d      [     -]
    real(r8), intent(in)              :: p3        ! log-slope associated with h      [     -]
    real(r8), intent(in)              :: c2b       ! carbon to biomass multiplier     [kg/kgC]
    real(r8), intent(in)              :: dbh_maxh  ! diameter at maximum height       [    cm]
    real(r8), intent(out)             :: blmax     ! leaf biomass                     [   kgC]
    real(r8), intent(out), optional   :: dblmaxdd  ! leaf biomass derivative          [kgC/cm]
    !----Local variables
    real(r8)  :: duse


    !----Cap diameter
    duse = min(d, dbh_maxh)

    !----Calculate leaf biomass
    blmax = p1 * duse**p2 * h**p3 / c2b

    !----Calculate leaf biomass derivative if needed

    if (present(dblmaxdd))then
       if(d .ge. dbh_maxh)then
          dblmaxdd = 0._r8
       else
          dblmaxdd = blmax * (p2 / duse + p3 * dhdd / h)
       end if
    end if

    return
  end subroutine dh2blmax_3pwr_grass
  
 
          
    
  ! =========================================================================
  ! Diameter to height (D2H) functions
  ! =========================================================================
  
  subroutine d2h_chave2014(d,p1,p2,p3,dbh_maxh,h,dhdd)
    
    ! "d2h_chave2014"
    ! "d to height via Chave et al. 2014"
    
    ! This function calculates tree height based on tree diameter and the
    ! environmental stress factor "E", as per Chave et al. 2014 GCB
    ! As opposed to previous allometric models in ED, in this formulation
    ! we do not impose a hard cap on tree height.  But, maximum_height
    ! is an important parameter, but instead of imposing a hard limit, in
    ! the new methodology, it will be used to trigger a change in carbon
    ! balance accounting.  Such that a tree that hits its maximum height will
    ! begin to route available NPP into seed and defense respiration.
    !
    ! The stress function is based on the geographic location of the site.  If
    ! a user decides to use Chave2014 allometry, the E factor will be read in
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
    ! p1 =  0.893 - E
    ! p2 =  0.76
    ! p3 = -0.034
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
    !  asymtotic functions (Weibull function)"
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
    ! biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2014.
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
    ! Chave's Paper has p1 = 0.0673, p2 = 0.976
    !
    
    ! =========================================================================
    
    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: h       ! plant height [m]
    real(r8),intent(in)  :: dhdd    ! change in height wrt diameter
    real(r8),intent(in)  :: p1  ! allometry parameter 1
    real(r8),intent(in)  :: p2  ! allometry parameter 2
    real(r8),intent(in)  :: wood_density
    real(r8),intent(in)  :: c2b
    real(r8),intent(out) :: bagw     ! plant aboveground biomass [kgC]
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
  
  ! ============================================================================



  subroutine dh2bagw_3pwr(d,h,dhdd,p1,p2,p3,wood_density,c2b,bagw,dbagwdd)
     !--------------------------------------------------------------------------
     !
     !     This function calculates the maximum above-ground biomass from 
     ! reference diameter, plant height and wood density. This functional form
     ! is an intermediate between Saldarriaga et al. (1988) and Chave et al.
     ! (2014), because the wood-density exponent is independent on the
     ! plant size (DBH^2 * Height) but the diameter and height are scaled 
     ! together through the size function.
     !
     ! -----------------
     !  References
     ! -----------------
     ! Chave J, Rejou-Mechain M, Burquez A, Chidumayo E, Colgan MS, Delitti WB,
     !    Duque A, Eid T, Fearnside PM, Goodman RC et al. 2014. Improved
     !    allometric models to estimate the aboveground biomass of tropical
     !    trees. Glob. Change Biol. 20: 3177-3190. doi:10.1111/gcb.12629.
     !
     ! Saldarriaga JG, West DC, Tharp ML , Uhl C. 1988. Long-term 
     !    chronosequence of forest succession in the upper Rio Negro of
     !    Colombia and Venezuela. J. Ecol. 76: 938-958.
     !    doi:10.2307/2260625.
     !
     ! -----------------
     !  Input arguments
     ! -----------------
     ! d            -- Diameter at breast height           [     cm]
     ! h            -- Total tree height                   [      m]
     ! dhdd         -- Height derivative with dbh          [   m/cm]
     ! p1           -- Parameter 1 (log-intercept)         [     --]
     ! p2           -- Parameter 2 (power, or log-slope)   [     --]
     ! p3           -- Parameter 3 (power, or log-slope)   [     --]
     ! wood_density -- Wood density                        [  g/cm3]
     ! c2b          -- Carbon to biomass multiplier ~ 2    [ kg/kgC]
     !
     ! ------------------
     !  Output arguments
     ! ------------------
     ! bagw         -- Above-ground biomass per individual [    kgC]
     ! dbagwdd      -- Above-ground biomass derivative     [ kgC/cm]
     !
     !--------------------------------------------------------------------------


     !--- Arguments
     real(r8), intent(in)            :: d            ! plant diameter              [    cm]
     real(r8), intent(in)            :: h            ! plant height                [     m]
     real(r8), intent(in)            :: dhdd         ! Height deriv. wrt diameter  [  m/cm]
     real(r8), intent(in)            :: p1           ! Log-intercept parameter     [     -]
     real(r8), intent(in)            :: p2           ! Log-slope parameter (size)  [     -]
     real(r8), intent(in)            :: p3           ! Log-slope parameter (WD)    [     -]
     real(r8), intent(in)            :: wood_density ! Wood density                [ g/cm3]
     real(r8), intent(in)            :: c2b          ! Carbon to biomass factor    [kg/kgC]
     real(r8), intent(out)           :: bagw         ! Above-ground biomass        [   kgC]
     real(r8), intent(out), optional :: dbagwdd      ! AG biomass derivative       [kgC/cm]
     !---~---


     !--- Find above-ground biomass
     bagw = p1 * (d*d*h)**p2 * wood_density**p3 / c2b
     !---~---


     !---~---
     !   Compute the above-ground biomass derivative with DBH if needed, noting that
     ! height is actually a function of DBH.
     !---~---
     if (present(dbagwdd)) then
        dbagwdd = p2 * bagw * ( 2._r8 / d + dhdd / h )
     end if
     !---~---

     return
  end subroutine dh2bagw_3pwr


  ! ============================================================================


  subroutine dh2bagw_3pwr_grass(d,h,dhdd,p1,p2,p3,c2b,bagw,dbagwdd)
    !---------------------------------------------------------------------------
    !
    ! This function calculates aboveground biomass (excluding leaf biomass) using
    ! basal diamerer (cm) and plant height (m) as size references, specifically
    ! for grass or herbaceous plants (can be used for other PFTs if supported by data)
    !
    !----------------
    ! Reference
    !----------------
    ! Gao X., Koven C., and Kueppers L. 2024. Allometric relationships and trade-offs in 11
    ! common Mediterranean-climate grasses. Ecological Applications.
    ! https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2976
    !
    !----------------
    ! Input arguments
    !----------------
    ! d                    -- Basal diameter                      [    cm]
    ! h                    -- Plant height                        [     m]
    ! dhdd                 -- Height derivative with diameter     [  m/cm]
    ! p1                   -- Log-intercept                       [     -]
    ! p2                   -- Log-slope associated with d         [     -]
    ! p3                   -- Log-slope associated with h         [     -]
    ! c2b                  -- Carbon to biomass multiplier        [kg/kgC]
    !
    !----------------
    ! Output variables
    !----------------
    ! bagw                 -- Aboveground biomass                 [   kgC]
    ! dbagwdd              -- Aboveground biomass derivative      [kgC/cm]
    !
    !---------------------------------------------------------------------------


    !----Arguments
    real(r8), intent(in)              :: d               ! plant diameter                      [    cm]
    real(r8), intent(in)              :: h               ! plant height                        [     m]
    real(r8), intent(in)              :: dhdd            ! height derivative w/ diameter       [  m/cm]
    real(r8), intent(in)              :: p1              ! log-intercept                       [     -]
    real(r8), intent(in)              :: p2              ! log-slope associated with d         [     -]
    real(r8), intent(in)              :: p3              ! log-slope associated with h         [     -]
    real(r8), intent(in)              :: c2b             ! biomass to carbon multiplier        [kg/kgC]
    real(r8), intent(out)             :: bagw            ! aboveground biomass excluding leaf  [   kgC]
    real(r8), intent(out),optional    :: dbagwdd         ! aboveground biomass derivative      [kgC/cm]

    !----Calculate aboveground biomass

    bagw = p1 * (d**p2) * (h**p3) / c2b

    !----Compute the aboveground biomass derivative with basal diameter if needed
    if (present(dbagwdd)) then
       dbagwdd = p2 * bagw / d + p3 * bagw * dhdd / h
    end if

    return
  end subroutine dh2bagw_3pwr_grass
  
    

  ! ============================================================================

  

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
    real(r8),intent(in)  :: p1      ! allometry parameter 1
    real(r8),intent(in)  :: p2      ! allometry parameter 2
    real(r8),intent(in)  :: c2b     ! carbon to biomass multiplier ~2
    real(r8),intent(out) :: bagw    ! plant aboveground biomass [kg C]
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


  ! =====================================================================================


  subroutine CrownDepth(height,ipft,crown_depth)

     !--------------------------------------------------------------------------
     !    This routine returns the depth of a plant's crown.  Which is the 
     ! length from the bottom of the crown to the top in the vertical dimension.
     ! 
     !    The original mode (now allom_dmode = 1) used only a fraction of the
     ! plant's height. Alternatively, allom_dmode = 2 uses the same functional
     ! form as Poorter et al. (2006), with an additional constraint to prevent
     ! crown depth to exceed the plant's height.
     !
     ! -----------------
     !  References
     ! -----------------
     ! Poorter L, Bongers L , Bongers F. 2006. Architecture of 54 moist-forest
     !    tree species: traits, trade-offs, and functional groups. Ecology 87:
     !    1289-1301. doi:10.1890/0012-9658(2006)87[1289:AOMTST]2.0.CO;2.
     !
     !--------------------------------------------------------------------------

     real(r8),intent(in)  :: height      ! The height of the plant   [m]
     integer ,intent(in)  :: ipft        ! functional type index
     real(r8),intent(out) :: crown_depth ! The depth of the crown    [m]

     associate( p1          => prt_params%allom_h2cd1(ipft), &
                p2          => prt_params%allom_h2cd2(ipft), &
                allom_dmode => prt_params%allom_dmode(ipft))

        select case (allom_dmode)
        case (1) ! Default, linear relationship with height
           crown_depth = p1 * height
        case (2) ! Power law, akin to Poorter et al. (2006).
           !---~---
           !   Apply the two coefficients, but make sure crown depth does not exceed 
           ! the plant's height.
           !---~---
           crown_depth = min(height, p1 * height ** p2)
           !---~---
        case default
          write(fates_log(),*) 'Invalid settings for crown depth mode for PFT ',ipft,'.'
          write(fates_log(),*) 'Current allom_dmode: ',allom_dmode,'.'
          write(fates_log(),*) 'Valid allom_dmode values: 1 or 2.'
          write(fates_log(),*) 'Aborting'
          call endrun(msg=errMsg(sourcefile, __LINE__))
        end select
     end associate

     return
  end subroutine CrownDepth
  ! =====================================================================================





  ! =============================================================================
  ! Specific diameter to crown area allometries
  ! =============================================================================

  
  subroutine carea_2pwr(dbh,spread,d2bl_p2,d2bl_ediff,d2ca_min, & 
                       d2ca_max,crowndamage,c_area,inverse)

     ! ============================================================================
     ! Calculate area of ground covered by entire cohort. (m2)
     ! Function of DBH (cm) canopy spread (m/cm) and number of individuals. 
     ! ============================================================================

     real(r8),intent(inout) :: dbh      ! diameter at breast (refernce) height [cm]
     real(r8),intent(in) :: spread      ! site level relative spread score [0-1]
     real(r8),intent(in) :: d2bl_p2     ! parameter 2 in the diameter->bleaf allometry (exponent)
     real(r8),intent(in) :: d2bl_ediff  ! area difference factor in the diameter-bleaf allometry (exponent)
     real(r8),intent(in) :: d2ca_min    ! minimum diameter to crown area scaling factor
     real(r8),intent(in) :: d2ca_max    ! maximum diameter to crown area scaling factor
     integer,intent(in)  :: crowndamage ! crowndamage class [1: undamaged, >1: damaged]
     real(r8),intent(inout) :: c_area   ! crown area for one plant [m2]
     logical,intent(in)  :: inverse     ! if true, calculate dbh from crown area rather than its reverse
     
     real(r8)            :: crown_area_to_dbh_exponent
     real(r8)            :: spreadterm  ! Effective 2bh to crown area scaling factor
     real(r8)            :: crown_reduction
     
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
     
     if ( .not. inverse) then
        c_area = spreadterm * dbh ** crown_area_to_dbh_exponent

        if(crowndamage > 1) then
           call GetCrownReduction(crowndamage, crown_reduction)
           c_area = c_area * (1.0_r8 - crown_reduction)
        end if
        
     else
        if(crowndamage > 1) then
           call GetCrownReduction(crowndamage, crown_reduction)
           c_area = c_area/(1.0_r8 - crown_reduction)
        end if
        dbh = (c_area / spreadterm) ** (1./crown_area_to_dbh_exponent)
     endif
     
  end subroutine carea_2pwr


  ! =============================================================================


  subroutine carea_3pwr(dbh,height,ipft,dbh_maxh,spread,dh2bl_p2,dh2bl_ediff, &
                        dh2ca_min,dh2ca_max,crowndamage,c_area,inverse)
     !---~---
     !    Calculate area of ground covered by entire cohort. (m2)
     ! Function of DBH (cm), height (m), canopy spread (m/cm) and number of 
     ! individuals. 
     !---~---

     !--- List of arguments
     real(r8)   , intent(inout) :: dbh         ! Diameter at breast/ref/ height     [   cm]
     real(r8)   , intent(inout) :: height      ! Height                             [    m]
     integer(i4), intent(in)    :: ipft        ! PFT index
     real(r8)   , intent(in)    :: dbh_maxh    ! Minimum DBH at maximum height      [   cm]
     real(r8)   , intent(in)    :: spread      ! site level relative spread score   [  0-1]
     real(r8)   , intent(in)    :: dh2bl_p2    ! Exponent for size (bleaf)          [    -]
     real(r8)   , intent(in)    :: dh2bl_ediff ! Difference in size exponent        [    -]
                                               !    between crown area and bleaf
     real(r8)   , intent(in)    :: dh2ca_min   ! Minimum (closed forest) scaling    [    -]
                                               !    coefficient for crown area
     real(r8)   , intent(in)    :: dh2ca_max   ! Maximum (savannah) scaling         [    -]
                                               !    coefficient for crown area
     integer    , intent(in)    :: crowndamage ! Crown damage class                 [    -]
                                               !    [1: undamaged, >1: damaged]
     real(r8)   , intent(inout) :: c_area      ! crown area for one plant           [   m2]
     logical    , intent(in)    :: inverse     ! If true, calculate dbh from crown
                                               !    area rather than its reverse
     !--- Local variables
     real(r8) :: size            ! Size (Diameter^2 * Height)                       [cm2 m]
     real(r8) :: dh2ca_p1        ! Effective scaling factor (crown area)            [    -]
     real(r8) :: dh2ca_p2        ! Effective exponent (crown area)                  [    -]
     real(r8) :: crown_reduction ! Crown area reduction due to damage.              [    -]
     !---~---


     !---~---
     !   Define the scaling (log-intercept) and exponent (log-slope) parameters for
     ! crown area. The scaling parameter accounts for the site-level spread elasticity.
     ! The exponent is defined in terms of the leaf biomass exponent plus an offset 
     ! parameter (allom_blca_expnt_diff). This is done because the default in FATES is
     ! for both exponents to be same (i.e., allom_blca_expnt_diff = 0.) so the per-plant 
     ! canopy area remains invariant during growth. However, allometric models in general
     ! predict that leaf area grows faster than crown area. 
     !---~---
     dh2ca_p1 = spread * dh2ca_max + (1._r8 - spread) * dh2ca_min
     dh2ca_p2 = dh2bl_p2 + dh2bl_ediff
     !---~---



     !---~---
     !   Decide whether to use DBH and height to find crown area (default) or the
     ! other way round.
     !---~---
     select case (inverse)
     case (.false.)
        !--- Find the maximum area
        size   = dbh * dbh * height
        c_area = dh2ca_p1 * size ** dh2ca_p2
        !---~---

        !--- Reduce area if the crown is damaged.
        if (crowndamage > 1) then
           call GetCrownReduction(crowndamage, crown_reduction)
           c_area = c_area * (1.0_r8 - crown_reduction)
        end if
        !---~---
     case (.true.)
        !--- Reduce area if the crown is damaged.
        if (crowndamage > 1) then
           call GetCrownReduction(crowndamage, crown_reduction)
           c_area = c_area * (1.0_r8 - crown_reduction)
        end if
        !---~---


        !---~---
        !   Find the size, then use a root-finding algorithm to find DBH.
        !---~---
        size = ( c_area / dh2ca_p1 ) ** ( 1.0_r8 / dh2ca_p2 )
        call size2dbh(size,ipft,dbh,dbh_maxh)
        !---~---
     end select
     !---~---

     return
  end subroutine carea_3pwr


  ! =========================================================================

  subroutine set_root_fraction(root_fraction, ft, zi, max_nlevroot)

    !
    ! !DESCRIPTION:
    !  Calculates the fractions of the root biomass in each layer for each pft. 
    !  It assumes an exponential decay.  If the soil depth is shallower than
    !  then exponential attenuation function, then it will normalize
    !  the profile and divide through.
    !
    ! !USES:

    !
    ! !ARGUMENTS
    real(r8),intent(inout) :: root_fraction(:) ! Normalized profile
    integer, intent(in)    :: ft               ! functional typpe
    real(r8),intent(in)    :: zi(0:)            ! Center of depth [m]

    ! The soil may not be active over the soil whole column due to things
    ! like permafrost. If so, compress profile over the maximum depth
    integer,optional, intent(in)   :: max_nlevroot  

    
    ! locals
    real(r8) :: a_par  ! local temporary for "a" parameter
    real(r8) :: b_par  ! ""  "b" parameter
    
    ! Parameters
    !
    ! TO-DO: NEXT TIME WE ROLL OUT A NEW PARAMETER INTERFACE, ADD
    ! PROFILE SWAPPING FLAGS.  OR IF THERE IS NO DEMAND< LEAVE AS IS.
    !
    !
    ! Two context exist 'hydraulic' and 'biomass'.  This allows us to
    ! allow different profiles for how water is drawn from the soil
    ! and different profiles to define the biomass for litter flux.
    ! These two context can currently choose 1 of the following three
    ! methods of defining the profile: 1) A 1 parameter exponential, 2)
    ! a beta profile defined by Jackson et al. and 3) a 2 parameter
    ! exponential.
    ! All methods return a normalized profile.

    integer, parameter :: jackson_beta_profile_type   = 1
    integer, parameter :: exponential_1p_profile_type = 2
    integer, parameter :: exponential_2p_profile_type = 3

    integer :: root_profile_type
    integer :: corr_id(1)        ! This is the bin with largest fraction
                                 ! add/subtract any corrections there
    integer :: nlevroot
    real(r8) :: correction       ! This correction ensures that root fractions
                                 ! sum to 1.0

    !----------------------------------------------------------------------
    
    if(size(zi) .ne. (size(root_fraction)+1)) then
       write(fates_log(),*) 'layer interface array should be 1 larger than'
       write(fates_log(),*) 'root fraction array'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    nlevroot = ubound(zi,1)
    
    ! Set root fraction to zero in all layers, as some may be inactive
    ! and we will only calculate the profiles over those 
    root_fraction(:) = 0._r8

    if(present(max_nlevroot))then
       if(debug .and. max_nlevroot<0)then
          write(fates_log(),*) 'A maximum rooting layer depth <0 was specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       nlevroot = min(max_nlevroot,nlevroot)
    end if
    
    select case(nint(prt_params%fnrt_prof_mode(ft)))
    case ( exponential_1p_profile_type ) 
       call exponential_1p_root_profile(root_fraction(1:nlevroot), zi(0:nlevroot), prt_params%fnrt_prof_a(ft)) 
    case ( jackson_beta_profile_type )
       call jackson_beta_root_profile(root_fraction(1:nlevroot), zi(0:nlevroot), prt_params%fnrt_prof_a(ft))
    case ( exponential_2p_profile_type ) 
       call exponential_2p_root_profile(root_fraction(1:nlevroot), zi(0:nlevroot), & 
             prt_params%fnrt_prof_a(ft),prt_params%fnrt_prof_b(ft))

    case default
       write(fates_log(),*) 'An undefined root profile type was specified'
       write(fates_log(),*) 'Aborting'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select


    correction = 1._r8 - sum(root_fraction)
    corr_id = maxloc(root_fraction)
    root_fraction(corr_id(1)) = root_fraction(corr_id(1)) + correction



    return
  end subroutine set_root_fraction

  ! =====================================================================================
    
  subroutine exponential_2p_root_profile(root_fraction, zi, a, b)

    !
    ! !ARGUMENTS
    real(r8),intent(out) :: root_fraction(:)
    real(r8),intent(in)  :: zi(0:)
    real(r8),intent(in)  :: a     ! Exponential shape parameter a
    real(r8),intent(in)  :: b     ! Exponential shape parameter b

    ! Locals
    integer  :: nlevsoil    ! Number of soil layers
    integer  :: lev         ! soil layer index
    real(r8) :: sum_rootfr  ! sum of root fraction for normalization


    ! Original default parameters: 
    !
    ! broadleaf_evergreen_tropical_tree
    ! needleleaf_evergreen_extratrop_tree
    ! needleleaf_colddecid_extratrop_tree
    ! broadleaf_evergreen_extratrop_tree
    ! broadleaf_hydrodecid_tropical_tree
    ! broadleaf_colddecid_extratrop_tree
    ! broadleaf_evergreen_extratrop_shrub
    ! broadleaf_hydrodecid_extratrop_shrub
    ! broadleaf_colddecid_extratrop_shrub
    ! arctic_c3_grass
    ! cool_c3_grass
    ! c4_grass
    !
    ! a = 7, 7, 7, 7, 6, 6, 7, 7, 7, 11, 11, 11 ;
    ! b = 1, 2, 2, 1, 2, 2, 1.5, 1.5, 1.5, 2, 2, 2 ;


    nlevsoil = ubound(zi,1)

    sum_rootfr = 0.0_r8
    do lev = 1, nlevsoil
       root_fraction(lev) = .5_r8*( &
             exp(-a * zi(lev-1))  &
             + exp(-b * zi(lev-1))  &
             - exp(-a * zi(lev))    &
             - exp(-b * zi(lev)))

       sum_rootfr = sum_rootfr + root_fraction(lev)
    end do

    ! Normalize the root profile
    root_fraction(1:nlevsoil) = root_fraction(1:nlevsoil)/sum_rootfr
    
    return
  end subroutine exponential_2p_root_profile

  ! =====================================================================================
  
  subroutine exponential_1p_root_profile(root_fraction, zi, a)

    !
    ! !ARGUMENTS
    real(r8),intent(out) :: root_fraction(:)
    real(r8),intent(in)  :: zi(0:)
    real(r8),intent(in)  :: a     ! Exponential shape parameter a

    !
    ! LOCAL VARIABLES:
    integer :: lev         ! soil depth layer index
    integer :: nlevsoil    ! number of soil layers
    real(r8) :: depth      ! Depth to middle of layer [m]
    real(r8) :: sum_rootfr ! sum of rooting profile for normalization

    ! Typical default parameter is a = 3.  
    ! how steep profile is
    ! for root C inputs (1/ e-folding depth) (1/m)
    
    nlevsoil = ubound(zi,1)
    
    ! define rooting profile from exponential parameters
    sum_rootfr = 0.0_r8
    do lev = 1,  nlevsoil
       root_fraction(lev) = exp(-a * 0.5*(zi(lev)+zi(lev-1)) )
       sum_rootfr = sum_rootfr + root_fraction(lev)
    end do
    
    ! Normalize the root profile
    root_fraction(1:nlevsoil) = root_fraction(1:nlevsoil)/sum_rootfr
    
    
    return
  end subroutine exponential_1p_root_profile
    
  ! =====================================================================================

  subroutine jackson_beta_root_profile(root_fraction, zi, a)

    ! -----------------------------------------------------------------------------------
    ! use beta distribution parameter from Jackson et al., 1996
    ! -----------------------------------------------------------------------------------
    ! !ARGUMENTS
    real(r8),intent(out) :: root_fraction(:) ! fraction of root mass in each soil layer
    real(r8),intent(in)  :: zi(0:)           ! depth of layer interfaces 0-nlevsoil
    real(r8),intent(in)  :: a                ! Exponential shape parameter a

    !
    ! LOCAL VARIABLES:
    integer :: lev         ! soil depth layer index
    integer :: nlevsoil    ! number of soil layers
    real(r8) :: sum_rootfr ! sum of rooting profile, for normalization 
    

    ! Original defaults in fates, a = 0.976 (all Pfts)

    nlevsoil = ubound(zi,1)
   
    sum_rootfr = 0.0_r8
    do lev = 1, nlevsoil
       root_fraction(lev) = &
             ( a ** ( zi(lev-1)*100._r8) - a ** ( zi(lev)*100._r8) )
       sum_rootfr = sum_rootfr + root_fraction(lev)
    end do
    
    ! Normalize the root profile
    root_fraction(1:nlevsoil) = root_fraction(1:nlevsoil)/sum_rootfr

    return
  end subroutine jackson_beta_root_profile

  ! =====================================================================================

  
  real(r8) function decay_coeff_vcmax(vcmax25top,slope_param,intercept_param)
    
    ! ---------------------------------------------------------------------------------
    ! This function estimates the decay coefficient used to estimate vertical
    ! attenuation of properties in the canopy.
    !
    ! Decay coefficient (kn) is a function of vcmax25top for each pft.
    !
    ! Currently, this decay is applied to vcmax attenuation, SLA (optionally)
    ! and leaf respiration (optionally w/ Atkin)
    !
    ! ---------------------------------------------------------------------------------
    
    !ARGUMENTS

    real(r8),intent(in) :: vcmax25top
    real(r8),intent(in) :: slope_param      ! multiplies vcmax25top
    real(r8),intent(in) :: intercept_param  ! adds to vcmax25top

    
    !LOCAL VARIABLES
    ! -----------------------------------------------------------------------------------
    
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
    ! kn = 0.11. Here, we derive kn from vcmax25 as in Lloyd et al 
    ! (2010) Biogeosciences, 7, 1833-1859
    ! This function is also used to vertically scale leaf maintenance
    ! respiration.
    
    decay_coeff_vcmax = exp(slope_param * vcmax25top - intercept_param)
    
    return
  end function decay_coeff_vcmax
  
  ! =====================================================================================

  subroutine ForceDBH( ipft, crowndamage, canopy_trim, elongf_leaf, elongf_stem, d, h, bdead, bl )

     ! =========================================================================
     ! This subroutine estimates the diameter based on either the structural biomass
     ! (if woody) or the leaf biomass using the allometric 
     ! functions. Since allometry is specified with diameter
     ! as the independent variable, we must do this through a search algorithm.
     ! Here, we keep searching until the difference between actual structure and
     ! the predicted structure based on the searched diameter is within a tolerance.
     ! ============================================================================
  use FatesConstantsMod     , only : calloc_abs_error
     ! Arguments


     integer(i4),intent(in)        :: ipft  ! PFT index
     integer(i4),intent(in)        :: crowndamage ! crowndamage [1: undamaged, >1: damaged]
     real(r8),intent(in)           :: canopy_trim
     real(r8),intent(in)           :: elongf_leaf ! Elongation factor: leaves (phenology)
     real(r8),intent(in)           :: elongf_stem ! Elongation factor: stem (phenology)
     real(r8),intent(inout)        :: d     ! plant diameter [cm]
     real(r8),intent(out)          :: h     ! plant height
     real(r8),intent(in),optional  :: bdead ! Structural biomass
     real(r8),intent(in),optional  :: bl    ! Leaf biomass
   
     
     ! Locals
     real(r8)  :: bt_sap,dbt_sap_dd  ! target sap wood at current d
     real(r8)  :: bt_agw,dbt_agw_dd  ! target AG wood at current d
     real(r8)  :: bt_bgw,dbt_bgw_dd  ! target BG wood at current d
     real(r8)  :: bt_dead,dbt_dead_dd ! target struct wood at current d
     real(r8)  :: bt_leaf,dbt_leaf_dd ! target leaf at current d
     real(r8)  :: at_sap              ! sapwood area (dummy) m2
     real(r8)  :: dd                  ! diameter increment for each step
     real(r8)  :: d_try               ! trial diameter
     real(r8)  :: bt_dead_try         ! trial structure biomasss
     real(r8)  :: dbt_dead_dd_try     ! trial structural derivative
     real(r8)  :: bt_leaf_try         ! trial leaf biomass
     real(r8)  :: dbt_leaf_dd_try     ! trial leaf derivative
     real(r8)  :: step_frac           ! step fraction
     integer   :: counter 
     real(r8), parameter :: step_frac0  = 0.9_r8
     integer, parameter  :: max_counter = 200
  
     
     ! Do reduce "if" calls, we break this call into two parts
     if ( prt_params%woody(ipft) == itrue ) then

        if(.not.present(bdead)) then
           write(fates_log(),*) 'woody plants must use structure for dbh reset'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
        
        call bsap_allom(d,ipft,crowndamage, canopy_trim, elongf_stem,at_sap,bt_sap,dbt_sap_dd)
        call bagw_allom(d,ipft,crowndamage, elongf_stem, bt_agw,dbt_agw_dd)
        call bbgw_allom(d,ipft, elongf_stem,bt_bgw,dbt_bgw_dd)

        call bdead_allom(bt_agw,bt_bgw, bt_sap, ipft, bt_dead, dbt_agw_dd, &
             dbt_bgw_dd, dbt_sap_dd, dbt_dead_dd)

        ! This calculates a diameter increment based on the difference
        ! in structural mass and the target mass, and sets it to a fraction
        ! of the diameter increment
        counter = 0
        step_frac = step_frac0
        do while( (bdead-bt_dead) > calloc_abs_error .and. dbt_dead_dd>0.0_r8)
           
           ! vulnerable to div0
           dd    = step_frac*(bdead-bt_dead)/dbt_dead_dd
           d_try = d + dd
        
           call bsap_allom(d_try,ipft,crowndamage, canopy_trim, elongf_stem,at_sap, &
                bt_sap,dbt_sap_dd)
           call bagw_allom(d_try,ipft,crowndamage, elongf_stem,  bt_agw,dbt_agw_dd)
           call bbgw_allom(d_try,ipft, elongf_stem, bt_bgw,dbt_bgw_dd)


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
        

     else

        if(.not.present(bl)) then
           write(fates_log(),*) 'grasses must use leaf for dbh reset'
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        call bleaf(d,ipft,crowndamage,canopy_trim,elongf_leaf,bt_leaf,dbt_leaf_dd)

        counter = 0
        step_frac = step_frac0
        do while( (bl-bt_leaf) > calloc_abs_error .and. dbt_leaf_dd>0.0_r8)

           dd    = step_frac*(bl-bt_leaf)/dbt_leaf_dd
           d_try = d + dd
           
           call bleaf(d_try,ipft,crowndamage,canopy_trim,elongf_stem,bt_leaf_try,dbt_leaf_dd_try)

           ! Prevent overshooting                                                                                           
           if(bt_leaf_try > (bl+calloc_abs_error)) then
              step_frac = step_frac*0.5_r8
           else
              step_frac = step_frac0
              d         = d_try
              bt_leaf   = bt_leaf_try
              dbt_leaf_dd = dbt_leaf_dd_try
           end if
           counter = counter + 1
           if (counter>max_counter) then
              write(fates_log(),*) 'Having trouble converging on dbh reset'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
        end do

     end if

     call h_allom(d,ipft,h)
     if(counter>20)then
        write(fates_log(),*) 'dbh counter: ',counter,' is woody: ',&
             (prt_params%woody(ipft) == itrue)

        if(prt_params%woody(ipft)==itrue)then
           warn_msg = 'dbh counter: '//trim(I2S(counter))//' is woody'
        else
           warn_msg = 'dbh counter: '//trim(I2S(counter))//' is not woody'
        end if
        call FatesWarn(warn_msg,index=3)
     end if

     

     return
  end subroutine ForceDBH

  ! =========================================================================

  subroutine VegAreaLayer(tree_lai,tree_sai,tree_height,iv,nv,pft,snow_depth, & 
       vai_top,vai_bot, elai_layer,esai_layer,tlai_layer,tsai_layer)

    ! -----------------------------------------------------------------------------------
    ! This routine returns the exposed leaf and stem areas (m2 of leaf and stem) per m2 of
    ! ground inside the crown, for the leaf-layer specified.
    ! -----------------------------------------------------------------------------------

    real(r8),intent(in) :: tree_lai         ! the in-crown leaf area index for the plant
                                            ! [m2 leaf/m2 crown footprint]
    real(r8),intent(in) :: tree_sai         ! the in-crown stem area index for the plant
                                            ! [m2 stem/m2 crown footprint]
    real(r8),intent(in) :: tree_height      ! the height of the plant [m]
    integer,intent(in)  :: iv               ! vegetation layer index
    integer,intent(in)  :: nv               ! this plants total number of veg layers
    integer,intent(in)  :: pft              ! plant functional type index
    real(r8),intent(in) :: snow_depth       ! the depth of snow on the ground [m]
    real(r8),intent(out) :: vai_top
    real(r8),intent(out) :: vai_bot          ! the VAI of the bin top and bottom
    real(r8),intent(out) :: elai_layer       ! exposed leaf area index of the layer
    real(r8),intent(out) :: esai_layer       ! exposed stem area index of the layer
    real(r8),optional,intent(out) :: tlai_layer       ! total leaf area index of the layer
    real(r8),optional,intent(out) :: tsai_layer       ! total stem area index of the layer

                                 ! [m2 of leaf in bin / m2 crown footprint]
    real(r8) :: tree_vai         ! the in-crown veg area index for the plant
    real(r8) :: crown_depth      ! crown depth of the plant [m]
    real(r8) :: frac_crown_depth ! fraction of the crown depth (relative to plant height)
    real(r8) :: fraction_exposed ! fraction of the veg media that is above snow
    real(r8) :: layer_top_height ! Physical height of the layer top relative to ground [m]
    real(r8) :: layer_bot_height ! Physical height of the layer bottom relative to ground [m]
    real(r8) :: tlai,tsai        ! temporary total area indices [m2/m2]
    real(r8) :: fleaf            ! fraction of biomass in layer that is leaf
    real(r8) :: remainder        ! old-method: remainder of biomass in last bin
    integer, parameter :: layer_height_const_depth = 1 ! constant physical depth assumption
    integer, parameter :: layer_height_const_lad   = 2 ! constant leaf area depth assumption
    integer, parameter :: layer_height_method = layer_height_const_depth
    
    tree_vai = tree_lai + tree_sai

    ! Ratio between crown depth and plant height
    call CrownDepth(tree_height,pft,crown_depth)
    frac_crown_depth = crown_depth / tree_height

    if_any_vai: if(tree_vai>0._r8)then

       if(iv==0)then
          vai_top = 0.0
          vai_bot = tree_vai
       else

          if(iv>1)then
             vai_top = dlower_vai(iv) - dinc_vai(iv)
          else
             vai_top = 0._r8
          end if

          if(iv<nv) then
             vai_bot = dlower_vai(iv)
          else
             vai_bot = tree_vai
          end if
       end if

       if(layer_height_method .eq. layer_height_const_depth)then
          if(iv==0)then
             layer_top_height = tree_height
             layer_bot_height = tree_height*(1._r8 - frac_crown_depth)
          else
             layer_top_height = tree_height*(1._r8 - real(iv-1,r8)/real(nv,r8)*frac_crown_depth)
             layer_bot_height = tree_height*(1._r8 - real(iv,r8)/real(nv,r8)*frac_crown_depth)
          end if
       else
          layer_top_height = tree_height*(1._r8 - frac_crown_depth*vai_top/tree_vai)
          layer_bot_height = tree_height*(1._r8 - frac_crown_depth*vai_bot/tree_vai)
       end if

       fraction_exposed =  1._r8 - max(0._r8,(min(1._r8, (snow_depth-layer_bot_height)/(layer_top_height-layer_bot_height))))

       tlai = (vai_bot-vai_top) * tree_lai / tree_vai
       tsai = (vai_bot-vai_top) * tree_sai / tree_vai

       if(present(tlai_layer)) tlai_layer = tlai
       if(present(tsai_layer)) tsai_layer = tsai

       elai_layer = fraction_exposed * tlai
       esai_layer = fraction_exposed * tsai

       ! Update the vai at the bottom to be removed/decreased if there is no exposure
       ! set the vai top and bottoms to the snow layer if below
       !vai_top = min(vai_top,fraction_exposed*tree_vai)

       vai_bot = vai_top + fraction_exposed*(vai_bot-vai_top)

    else

       if(present(tlai_layer)) tlai_layer = 0._r8
       if(present(tsai_layer)) tsai_layer = 0._r8
       elai_layer = 0._r8
       esai_layer = 0._r8
       vai_bot = 0._r8
       vai_top = 0._r8

    end if if_any_vai


    return
  end subroutine VegAreaLayer
  
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
   ! ============================================================================



   ! ============================================================================
   !    This function finds the DBH when size (DBH^2 * Height) is known but we
   ! cannot find DBH analytically due to the non-linear relationship between DBH
   ! and height. This is borrowed from the same approach applied in ED2 for
   ! root finding. It starts with the Newton's method, which should quickly
   ! converge to the solution. In the unlikely case of failure, we use the
   ! Regula Falsi (Illinois) method as a back-up.
   ! ============================================================================
   subroutine size2dbh(size,ipft,dbh,dbh_maxh)
      !--- Arguments.
      real(r8)   , intent(in)    :: size        ! Size (DBH^2 * Height)           [cm2 m]
      integer(i4), intent(in)    :: ipft        ! PFT index                       [    -]
      real(r8)   , intent(inout) :: dbh         ! Diameter at breast height       [   cm]
      real(r8)   , intent(in)    :: dbh_maxh    ! Minimum DBH at maximum height   [   cm]
      !--- Local variables
      real(r8)                   :: hgt         ! Height                          [    m]
      real(r8)                   :: dhgtddbh    ! Height derivative               [ m/cm]
      real(r8)                   :: size_maxh   ! Minimum size at maximum height  [cm2 m]
      real(r8)                   :: deriv       ! Function derivative             [ cm m]
      real(r8)                   :: afun        ! Function value (lower guess)    [cm2 m]
      real(r8)                   :: rfun        ! Function value (RF new guess)   [cm2 m]
      real(r8)                   :: zfun        ! Function value (upper guess)    [cm2 m]
      real(r8)                   :: adbh        ! DBH: lower guess                [   cm]
      real(r8)                   :: rdbh        ! DBH: updated guess (Reg. Falsi) [   cm]
      real(r8)                   :: zdbh        ! DBH: upper guess                [   cm]
      real(r8)                   :: delta       ! Second guess for the RF method  [   cm]
      integer                    :: itn         ! Iteration counter -- Newton     [    -]
      integer                    :: iti         ! Iteration counter -- Reg. Falsi [    -]
      logical                    :: converged   ! Has the solution converged?     [  T|F]
      logical                    :: zside       ! Converging on the upper size?   [  T|F]
      !--- Local constants.
      real(r8) , parameter :: toler =1.0e-12_r8 ! Relative tolerance              [   --]
      integer  , parameter :: maxit_newt = 10   ! Cap in iterations -- Newton     [   --]
      integer  , parameter :: maxit_rf   = 100  ! Cap in iterations -- Reg. Falsi [   --]
      !---~---


      !---~---
      !   Find the maximum size beyond which the height is assumed constant. In this
      ! case, DBH can be determined without the iterative approach.
      !---~---
      call h_allom(dbh_maxh,ipft,hgt)
      size_maxh = dbh_maxh * dbh_maxh * hgt
      if (size >= size_maxh) then
         dbh = sqrt(size/hgt)
         return
      end if
      !---~---


      !--- First guess: use current DBH.
      adbh  = dbh
      call h_allom(adbh,ipft,hgt,dhgtddbh)
      afun  = adbh * adbh * hgt - size
      deriv = 2.0_r8 * adbh * hgt + adbh * adbh * dhgtddbh
      !---~---


      !--- Copy just in case it fails at the first iteration.
      zdbh = adbh
      zfun = afun
      !---~---


      !---~---
      !   Enter the Newton's method loop
      !---~---
      converged = .false.
      newton_loop: do itn = 1, maxit_newt
         !--- If derivative is too flat, go to Regula Falsi
         if ( abs(deriv) < toler) exit newton_loop
         !---~---


         !--- Copy the previous guess.
         adbh = zdbh
         afun = zfun
         !---~---


         !--- Find the new guess, and evaluate the function and derivative.
         zdbh  = adbh - afun / deriv
         call h_allom(zdbh,ipft,hgt,dhgtddbh)
         zfun  = zdbh * zdbh * hgt - size
         deriv = 2.0_r8 * zdbh * hgt + zdbh * zdbh * dhgtddbh
         !---~---

         !--- Check convergence.
         converged = abs(adbh - zdbh) < toler * zdbh
         if (converged) then
            !--- Convergence by iterations.
            dbh = 0.5_r8 * (adbh + zdbh)
            return
            !---~---
         else if (abs(zfun) < nearzero) then
            !--- Convergence by luck.
            dbh = zdbh
            return
            !---~---
         end if
         !---~---
      end do newton_loop
      !---~---



      !---~---
      !   If we have reached this point, then Newton's method has failed. Use Regula
      ! Falsi instead. For this, we must have two guesses whose function evaluation has
      ! opposite signs.
      !---~---
      if (afun * zfun <= -nearzero) then
         !--- We already have two guesses with opposite signs.
         zside = .true.
         !---~---
      else
         !--- Look for another guess with opposite sign.
         if (abs(zfun-afun) < 100._r8 * toler * adbh) then
            delta = 100._r8 * toler * adbh
         else
            delta = max( abs( afun * (zdbh-adbh) / (zfun-afun) ),100._r8 * toler * adbh )
         end if
         !---~---


         !---~---
         !   Try guesses on both sides of the first guess, sending guesses increasingly
         ! further away until we find a good guess.
         !---~---
         zdbh  = adbh + delta
         zside = .false.
         zguess_loop: do iti=1,maxit_rf
            zdbh = adbh + real((-1)**iti * (iti+3)/2,r8) * delta
            call h_allom(zdbh,ipft,hgt)
            zfun  = zdbh * zdbh * hgt - size
            zside = afun * zfun < -nearzero
            if (zside) exit zguess_loop
         end do zguess_loop

         !---~---
         !   Issue an error in case the function failed finding a second guess.
         !---~---
         if (.not. zside) then
            write (unit=*,fmt='(a)')           '---~---'
            write (unit=*,fmt='(a)')           ' Failed finding the second guess:'
            write (unit=*,fmt='(a)')           '---~---'
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Input:   '
            write (unit=*,fmt='(a,1x,es14.7)') ' + size  =',size
            write (unit=*,fmt='(a,1x,es14.7)') ' + dbh   =',dbh
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Current guesses and evaluations:'
            write (unit=*,fmt='(a,1x,es14.7)') ' + adbh  =',adbh
            write (unit=*,fmt='(a,1x,es14.7)') ' + afun  =',afun
            write (unit=*,fmt='(a,1x,es14.7)') ' + zdbh  =',zdbh
            write (unit=*,fmt='(a,1x,es14.7)') ' + zfun  =',zfun
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           '---~---'
            write(fates_log(),*) 'Second guess for Regula Falsi method not found.'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         !---~---
      end if
      !---~---

      !---~---
      !   Proceed to the regula falsi loop.
      !---~---
      regfalsi_loop: do iti=1,maxit_rf

         !--- Update solution.
         rdbh =  ( zfun * adbh - afun * zdbh ) / ( zfun - afun)
         !---~---


         !---~---
         !   Check for convergence. In case it converged, we can exit the sub-routine.
         !---~---
         converged = abs(rdbh - adbh) < toler * max(rdbh,adbh)
         if (converged) exit regfalsi_loop
         !---~---


         !--- Find the new function evaluation.
         call h_allom(rdbh,ipft,hgt)
         rfun = rdbh * rdbh * hgt - size
         !---~---


         !---~---
         !   Define the new searching interval based on the intermediate value theorem.
         !---~---
         if (abs(rfun) < nearzero) then
            !--- Converged by luck.
            converged = .true.
            exit regfalsi_loop
            !---~---
         else if (rfun * afun <= -nearzero ) then
            !--- Guess is between lower and current guess.
            zdbh = rdbh
            zfun = rfun
            !--- If we are updating the upper side again, halve afun (Regula Falsi method).
            if (zside) afun = afun * 0.5_r8
            !--- Flag that we have just updated the upper side.
            zside = .true.
            !---~---
         else
            !--- Guess is between current and upper guess.
            adbh = rdbh
            afun = rfun
            !--- If we are updating the lower side again, halve zfun (Regula Falsi method).
            if (.not. zside) zfun = zfun * 0.5_r8
            !--- Flag that we have just updated the lower side.
            zside = .false.
         end if
      end do regfalsi_loop
      !---~---


      !---~---
      !   Check that the Regula Falsi method indeed converged.
      !---~---
      if (converged) then
         !--- Yes, return the last guess
         dbh = rdbh
         !---~---
      else
         !--- No, report the bad news
         write (unit=*,fmt='(a)')           '---~---'
         write (unit=*,fmt='(a)')           ' Failed finding the solution:'
         write (unit=*,fmt='(a)')           '---~---'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Input:   '
         write (unit=*,fmt='(a,1x,es14.7)') ' + size  =',size
         write (unit=*,fmt='(a,1x,es14.7)') ' + dbh   =',dbh
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Current guesses and evaluations:'
         write (unit=*,fmt='(a,1x,es14.7)') ' + adbh  =',adbh
         write (unit=*,fmt='(a,1x,es14.7)') ' + afun  =',afun
         write (unit=*,fmt='(a,1x,es14.7)') ' + rdbh  =',rdbh
         write (unit=*,fmt='(a,1x,es14.7)') ' + rfun  =',rfun
         write (unit=*,fmt='(a,1x,es14.7)') ' + zdbh  =',zdbh
         write (unit=*,fmt='(a,1x,es14.7)') ' + zfun  =',zfun
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           '---~---'
         write(fates_log(),*) 'Size to DBH routine failed to converge!'
         call endrun(msg=errMsg(sourcefile, __LINE__))
         !---~---
      end if
      !---~---

      return
   end subroutine size2dbh
   ! ============================================================================



end module FatesAllometryMod
