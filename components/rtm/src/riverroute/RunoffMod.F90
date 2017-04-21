module RunoffMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing utilities for history file and coupler runoff data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mct_mod
  use RtmVar         , only : iulog, spval, ispval
  use rtm_cpl_indices, only : nt_rtm

! !PUBLIC TYPES:
  implicit none
  private

  public :: runoff_flow
  type runoff_flow
     !    - local initialization
     real(r8), pointer :: lonc(:)   => null()       ! lon of cell
     real(r8), pointer :: latc(:)   => null()     ! lat of cell
     real(r8), pointer :: area(:)   => null()     ! area of cell
     integer , pointer :: gindex(:) => null()     ! global index
     integer , pointer :: dsi(:)    => null()     ! downstream index

     !    - local runtime
     real(r8), pointer :: runoff(:,:)     => null() ! RTM flow (m**3 H2O/s)
     real(r8), pointer :: runofflnd(:,:)  => null() ! runoff masked for land (m**3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)  => null() ! runoff masked for ocn  (m**3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)    => null() ! RTM change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:) => null() ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:) => null() ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:)       => null() ! RTM storage (m**3)
     real(r8), pointer :: volrlnd(:,:)    => null() ! RTM storage masked for land (m**3)
     real(r8), pointer :: fluxout(:,:)    => null() ! RTM cell tracer outlflux (m3/s)
     real(r8), pointer :: fthresh(:)      => null() ! RTM water flood threshold
     real(r8), pointer :: flood(:)        => null() ! RTM water (flood) sent back to clm (mm/s)

     !    - global 
     integer , pointer :: mask(:)         => null() ! mask of cell 0=none, 1=lnd, 2=ocn
     real(r8), pointer :: rlon(:)         => null() ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)         => null() ! rtm latitude list, 1d

     real(r8)          :: totarea          ! global area
     integer           :: numr             ! rtm gdc global number of cells
     integer           :: numrl            ! rtm gdc global number of lnd cells
     integer           :: numro            ! rtm gdc global number of ocn cells

     !    - local
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! local number of cells

     !    - 1d field pointers for history files (currently needed)
     real(r8), pointer :: runofflnd_nt1(:)  => null() 
     real(r8), pointer :: runofflnd_nt2(:)  => null() 
     real(r8), pointer :: runoffocn_nt1(:)  => null() 
     real(r8), pointer :: runoffocn_nt2(:)  => null() 
     real(r8), pointer :: dvolrdtlnd_nt1(:) => null() 
     real(r8), pointer :: dvolrdtlnd_nt2(:) => null() 
     real(r8), pointer :: dvolrdtocn_nt1(:) => null() 
     real(r8), pointer :: dvolrdtocn_nt2(:) => null() 
     real(r8), pointer :: volr_nt1(:)       => null() 
     real(r8), pointer :: volr_nt2(:)       => null() 
  end type runoff_flow
  !
  type (runoff_flow), public :: runoff

  public :: RunoffInit
  public :: RunoffFinalize

contains

  subroutine RunoffInit(begr, endr, numr)

    integer, intent(in) :: begr, endr, numr

    integer :: ier

    allocate(runoff%runoff(begr:endr,nt_rtm),     &
             runoff%dvolrdt(begr:endr,nt_rtm),    &
             runoff%runofflnd(begr:endr,nt_rtm),  &
             runoff%dvolrdtlnd(begr:endr,nt_rtm), &
             runoff%runoffocn(begr:endr,nt_rtm),  &
             runoff%dvolrdtocn(begr:endr,nt_rtm), &
             runoff%area(begr:endr),              &
             runoff%volr(begr:endr,nt_rtm),       &
             runoff%volrlnd(begr:endr,nt_rtm),    &
             runoff%fluxout(begr:endr,nt_rtm),    &
             runoff%lonc(begr:endr),              &
             runoff%latc(begr:endr),              &
             runoff%dsi(begr:endr),               &
             runoff%runofflnd_nt1(begr:endr),     &
             runoff%runofflnd_nt2(begr:endr),     &
             runoff%runoffocn_nt1(begr:endr),     &
             runoff%runoffocn_nt2(begr:endr),     &
             runoff%volr_nt1(begr:endr),          &
             runoff%volr_nt2(begr:endr),          &
             runoff%dvolrdtlnd_nt1(begr:endr),    &
             runoff%dvolrdtlnd_nt2(begr:endr),    &
             runoff%dvolrdtocn_nt1(begr:endr),    &
             runoff%dvolrdtocn_nt2(begr:endr),    &
             runoff%mask(numr),                   &
             runoff%gindex(begr:endr),            &
             runoff%fthresh(begr:endr),           &
             runoff%flood(begr:endr),             &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
       call shr_sys_abort
    end if

    runoff%runoff(:,:)     = 0._r8
    runoff%runofflnd(:,:)  = spval
    runoff%runoffocn(:,:)  = spval
    runoff%dvolrdt(:,:)    = 0._r8
    runoff%dvolrdtlnd(:,:) = spval
    runoff%dvolrdtocn(:,:) = spval
    runoff%volr(:,:)       = 0._r8
    runoff%volrlnd(:,:)    = 0._r8
    runoff%volr_nt1(:)     = 0._r8
    runoff%volr_nt2(:)     = 0._r8
    runoff%gindex(:)       = ispval
    runoff%fthresh(:)      = spval
    runoff%flood(:)        = 0._r8

  end subroutine RunoffInit

  subroutine RunoffFinalize()

    if (associated(runoff%runoff)) deallocate(runoff%runoff)
    if (associated(runoff%dvolrdt)) deallocate(runoff%dvolrdt)
    if (associated(runoff%runofflnd)) deallocate(runoff%runofflnd)
    if (associated(runoff%dvolrdtlnd)) deallocate(runoff%dvolrdtlnd)
    if (associated(runoff%runoffocn)) deallocate(runoff%runoffocn)
    if (associated(runoff%dvolrdtocn)) deallocate(runoff%dvolrdtocn)
    if (associated(runoff%area)) deallocate(runoff%area)
    if (associated(runoff%volr)) deallocate(runoff%volr)
    if (associated(runoff%volrlnd)) deallocate(runoff%volrlnd)
    if (associated(runoff%fluxout)) deallocate(runoff%fluxout)
    if (associated(runoff%lonc)) deallocate(runoff%lonc)
    if (associated(runoff%latc)) deallocate(runoff%latc)
    if (associated(runoff%rlon)) deallocate(runoff%rlon)
    if (associated(runoff%rlat)) deallocate(runoff%rlat)
    if (associated(runoff%dsi)) deallocate(runoff%dsi)
    if (associated(runoff%runofflnd_nt1)) deallocate(runoff%runofflnd_nt1)
    if (associated(runoff%runofflnd_nt2)) deallocate(runoff%runofflnd_nt2)
    if (associated(runoff%runoffocn_nt1)) deallocate(runoff%runoffocn_nt1)
    if (associated(runoff%runoffocn_nt2)) deallocate(runoff%runoffocn_nt2)
    if (associated(runoff%volr_nt1)) deallocate(runoff%volr_nt1)
    if (associated(runoff%volr_nt2)) deallocate(runoff%volr_nt2)
    if (associated(runoff%dvolrdtlnd_nt1)) deallocate(runoff%dvolrdtlnd_nt1)
    if (associated(runoff%dvolrdtlnd_nt2)) deallocate(runoff%dvolrdtlnd_nt2)
    if (associated(runoff%dvolrdtocn_nt1)) deallocate(runoff%dvolrdtocn_nt1)
    if (associated(runoff%dvolrdtocn_nt2)) deallocate(runoff%dvolrdtocn_nt2)
    if (associated(runoff%mask)) deallocate(runoff%mask)
    if (associated(runoff%gindex)) deallocate(runoff%gindex)
    if (associated(runoff%fthresh)) deallocate(runoff%fthresh)
    if (associated(runoff%flood)) deallocate(runoff%flood)

  end subroutine RunoffFinalize

end module RunoffMod
