! =======================================================================================
!
! This is the wrapper module that provides callable functions
! so that PARTEH fortran data structures can be intsantiated from python
! Half of the instantiation will occur by binding inherited data structures
! to cohorts, but the other half is the creation of a mapping table,
! of which we have only 1 per instance.  That happens here.
!
! Note: In FATES, the equivalent routine would probably live in FatesInterfaceMod.F90
!
! =======================================================================================

module FatesPARTEHWrapMod

  use PRTAllometricCarbonMod, only : InitPRTGlobalAllometricCarbon
  use PRTAllometricCNPMod, only    : InitPRTGlobalAllometricCNP
  use FatesGlobals        , only : endrun => fates_endrun
  use FatesGlobals        , only : fates_log
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use iso_c_binding, only : r8 => c_double
  use iso_c_binding, only : i4 => c_int
  use iso_c_binding, only : c_char

  implicit none
  private

  character(len=*), parameter, private :: sourcefile = __FILE__

  ! Make necessary procedures public

  public :: SPMapPyset

contains


  subroutine SPMapPyset()  !prt_mode)


    ! Update... Instantiate all of them?
    
!    integer(i4), intent(in) :: prt_mode
    
!    select case(int(prt_mode))
!    case (1)

     ! We actually initialize all hypotheses, since we are intercomparing.


    call InitPRTGlobalAllometricCarbon()

!    case(2)

    call InitPRTGlobalAllometricCNP()
       
!    case DEFAULT
!       write(fates_log(),*) 'You specified an unknown PRT module'
!       write(fates_log(),*) 'Aborting'
!       call endrun(msg=errMsg(sourcefile, __LINE__))
       
       
!    end select
    
  end subroutine SPMapPyset






end module FatesPARTEHWrapMod
