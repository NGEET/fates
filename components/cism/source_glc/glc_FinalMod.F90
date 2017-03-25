!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_FinalMod

!BOP
! !MODULE: glc_FinalMod
! !DESCRIPTION:
!  This module contains the glc finalization method that shuts down glc
!  gracefully (we hope).  It exits the message environment and checks 
!  for successful execution.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_FinalMod.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from POP_FinalMod.F90
!
! !USES:

   use glc_kinds_mod
   use glc_communicate, only: exit_message_environment
   use glc_fields, only: ice_sheet
   use glad_main, only: end_glad
   use glc_constants, only: stdout

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_final

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_final
! !INTERFACE:

 subroutine glc_final()

! !DESCRIPTION:
!  This routine shuts down glc by exiting all relevent environments.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: ierr

!-----------------------------------------------------------------------
!
!  exit glad gracefully
!
!-----------------------------------------------------------------------

   call end_glad(ice_sheet, close_logfile=.false.)

!-----------------------------------------------------------------------
!
!  exit the communication environment
!
!-----------------------------------------------------------------------

   call exit_message_environment(ierr)
   ! note that ierr isn't checked here...

!-----------------------------------------------------------------------
!EOC

 end subroutine glc_final

!***********************************************************************

 end module glc_FinalMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
