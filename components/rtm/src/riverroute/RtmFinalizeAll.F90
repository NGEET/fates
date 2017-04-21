module RtmFinalizeAll

!
! gateway to other Rtm routines to clean up memory.  Called from:
!
! rof_comp_mct::rof_final_mct
! rof_comp_esmf::rof_final_esmf
!

  implicit none
  private

  public :: RtmFinalizeMemory !  other Rtm routines to clean up memory

contains

  subroutine RtmFinalizeMemory

    !-----------------------------------------------------------------------
    use RunoffMod ,      only : RunoffFinalize
    use RtmMod ,         only : RtmFinalize
    use RtmRestFile,     only : RtmRestFinalize
    use RtmHistFile ,    only : RtmHistFileFinalize
    use RtmTimeManager,  only : timemgr_finalize
    use RtmVar,          only : rtm_active
    !-----------------------------------------------------------------------
    implicit none

    !
    ! deal with clean up of memory for parts of RTM here
    !
    if (rtm_active) then

       call RtmFinalize()
       call RunoffFinalize()
       call RtmRestFinalize()
       call RtmHistFileFinalize()
       !
       ! clean up ESMF clock memory here.  There is no
       ! equivalent for rof_final_esmf.
       ! For further info. please see comments in 
       ! RtmTimeManager.F90::timemgr_finalize
       !
       call timemgr_finalize()

    end if

  end subroutine RtmFinalizeMemory

end module RtmFinalizeAll
