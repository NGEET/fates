module FatesHydroWTFMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : pa_per_mpa
  use FatesConstantsMod, only : mpa_per_pa
  use FatesConstantsMod, only : mm_per_m
  use FatesConstantsMod, only : m_per_mm
  use FatesConstantsMod, only : denh2o => dens_fresh_liquid_water
  use FatesConstantsMod, only : grav_earth
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : pi_const
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log
  use shr_log_mod      , only : errMsg => shr_log_errMsg

  implicit none
  private

  ! -------------------------------------------------------------------------------------
  ! This module contains all unit (F)unctions associated with (W)ater (T)ransfer.
  ! e.g. WTFs
  ! These are also called "pedotransfer" functions, however, since these
  ! may be applied to xylems, stems, etc, they are not limited to soils (pedo).
  ! -------------------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  real(r8), parameter :: min_ftc = 0.0005_r8   ! Minimum allowed fraction of total conductance
                                               
  ! Bounds on saturated fraction, outside of which we use linear PV or stop flow
  ! In this context, the saturated fraction is defined by the volumetric WC "th"
  ! and the volumetric residual and saturation "th_res" and "th_sat": (th-th_r)/(th_sat-th_res)

  real(r8), parameter :: min_sf_interp = 0.02 ! Linear interpolation below this saturated frac
  real(r8), parameter :: max_sf_interp = 0.95 ! Linear interpolation above this saturated frac

  real(r8), parameter :: quad_a1 = 0.80_r8  ! smoothing factor "A" term
                                            ! in the capillary-elastic region

  real(r8), parameter :: quad_a2 = 0.99_r8  ! Smoothing factor or "A" term in
                                            ! elastic-caviation region


  ! Generic class that can be extended to describe
  ! specific water retention functions

  type, public :: wrf_type
   contains
     procedure :: th_from_psi     => th_from_psi_base
     procedure :: psi_from_th     => psi_from_th_base
     procedure :: dpsidth_from_th => dpsidth_from_th_base
     procedure :: set_wrf_param   => set_wrf_param_base
  end type wrf_type


  ! Generic class that can be extended to describe
  ! water conductance functions

  type, public :: wkf_type
   contains
     procedure :: ftc_from_psi      => ftc_from_psi_base
     procedure :: dftcdpsi_from_psi => dftcdpsi_from_psi_base
     procedure :: set_wkf_param     => set_wkf_param_base
  end type wkf_type

  ! The WRF and WKF types cannot be arrays themselves
  ! we require these holders

  type, public :: wrf_arr_type
     class(wrf_type), pointer :: p
     real(r8) :: th_sat
     real(r8) :: psi_sat
  end type wrf_arr_type

  type, public :: wkf_arr_type
      class(wkf_type), pointer :: p
  end type wkf_arr_type


  ! =====================================================================================
  ! Van Genuchten WTF Definitions
  ! =====================================================================================

  ! Water Retention Function
  type, public, extends(wrf_type) :: wrf_type_vg
     real(r8) :: alpha   ! Inverse air entry parameter         [m3/Mpa]
     real(r8) :: psd     ! Inverse width of pore size distribution parameter
     real(r8) :: th_sat  ! Saturation volumetric water content [m3/m3]
     real(r8) :: th_res  ! Residual volumetric water content   [m3/m3]

   contains
     procedure :: th_from_psi     => th_from_psi_vg
     procedure :: psi_from_th     => psi_from_th_vg
     procedure :: dpsidth_from_th => dpsidth_from_th_vg
     procedure :: set_wrf_param   => set_wrf_param_vg
  end type wrf_type_vg

  ! Water Conductivity Function
  type, public, extends(wkf_type) :: wkf_type_vg
     real(r8) :: alpha   ! Inverse air entry parameter         [m3/Mpa]
     real(r8) :: psd     ! Inverse width of pore size distribution parameter
     real(r8) :: tort    ! Tortuosity parameter (sometimes "l")
     real(r8) :: th_sat  ! Saturation volumetric water content [m3/m3]
     real(r8) :: th_res  ! Residual volumetric water content   [m3/m3]
   contains
     procedure :: ftc_from_psi      => ftc_from_psi_vg
     procedure :: dftcdpsi_from_psi => dftcdpsi_from_psi_vg
     procedure :: set_wkf_param     => set_wkf_param_vg
  end type wkf_type_vg

  ! =====================================================================================
  ! Clapp-Hornberger and Campbell (CCH) water retention and conductivity functions
  ! =====================================================================================

  ! Water Retention Function
  type, public, extends(wrf_type) :: wrf_type_cch
     real(r8) :: th_sat   ! Saturation volumetric water content         [m3/m3]
     real(r8) :: psi_sat  ! Bubbling pressure (potential at saturation) [Mpa]
     real(r8) :: beta     ! Clapp-Hornberger "beta" parameter           [-]
     real(r8) :: psi_max     ! psi where satfrac = max_sf_interp, and use linear
     real(r8) :: dpsidth_max ! deriv wrt theta for psi_max
   contains
     procedure :: th_from_psi     => th_from_psi_cch
     procedure :: psi_from_th     => psi_from_th_cch
     procedure :: dpsidth_from_th => dpsidth_from_th_cch
     procedure :: set_wrf_param   => set_wrf_param_cch
  end type wrf_type_cch

  ! Water Conductivity Function
  type, public, extends(wkf_type) :: wkf_type_cch
     real(r8) :: th_sat   ! Saturation volumetric water content         [m3/m3]
     real(r8) :: psi_sat  ! Bubbling pressure (potential at saturation) [Mpa]
     real(r8) :: beta     ! Clapp-Hornberger "beta" parameter           [-]
   contains
     procedure :: ftc_from_psi      => ftc_from_psi_cch
     procedure :: dftcdpsi_from_psi => dftcdpsi_from_psi_cch
     procedure :: set_wkf_param     => set_wkf_param_cch
  end type wkf_type_cch

  ! =====================================================================================
  ! TFS functions
  ! =====================================================================================

  ! Water Retention Function from TFS
  type, public, extends(wrf_type) :: wrf_type_tfs

     real(r8) :: th_sat   ! Saturation volumetric water content         [m3/m3]
     real(r8) :: th_res   ! Residual volumentric water content          [m3/m3]
     real(r8) :: pinot    ! osmotic potential at full turger            [MPa]
     real(r8) :: epsil    ! bulk elastic modulus                        [MPa]
     real(r8) :: rwc_ft   ! RWC @ full turgor, (elastic drainage begins)[-]
     real(r8) :: cap_corr ! correction for nonzero psi0x
     real(r8) :: cap_int  ! intercept of capillary region of curve
     real(r8) :: cap_slp  ! slope of capillary region of curve
     integer  :: pmedia   ! self describing porous media index

     real(r8) :: psi_max     ! psi matching max_sf_interp where we start linear interp
     real(r8) :: dpsidth_max ! dpsi_dth where we start linear interp
     real(r8) :: psi_min     ! psi matching min_sf_interp
     real(r8) :: dpsidth_min ! dpsi_dth where we start min interp

   contains
     procedure :: th_from_psi     => th_from_psi_tfs
     procedure :: psi_from_th     => psi_from_th_tfs
     procedure :: dpsidth_from_th => dpsidth_from_th_tfs
     procedure :: set_wrf_param   => set_wrf_param_tfs
     procedure :: bisect_pv
  end type wrf_type_tfs

  ! Water Conductivity Function
  type, public, extends(wkf_type) :: wkf_type_tfs
     real(r8) :: p50          ! matric potential at 50% conductivity loss   [Mpa]
     real(r8) :: avuln        ! vulnerability curve parameter
     real(r8) :: th_sat       ! volumetric water content at saturation

   contains
     procedure :: ftc_from_psi      => ftc_from_psi_tfs
     procedure :: dftcdpsi_from_psi => dftcdpsi_from_psi_tfs
     procedure :: set_wkf_param     => set_wkf_param_tfs
  end type wkf_type_tfs


contains

  ! =====================================================================================
  ! Functional definitions follow here
  ! Start off by writing the base types, which ultimately should never be pointed to.
  ! =====================================================================================

  subroutine set_wrf_param_base(this,params_in)
    class(wrf_type)     :: this
    real(r8),intent(in) :: params_in(:)
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end subroutine set_wrf_param_base
  subroutine set_wkf_param_base(this,params_in)
    class(wkf_type)     :: this
    real(r8),intent(in) :: params_in(:)
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end subroutine set_wkf_param_base
  function th_from_psi_base(this,psi) result(th)
    class(wrf_type)     :: this
    real(r8),intent(in) :: psi
    real(r8)            :: th
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function th_from_psi_base
  function psi_from_th_base(this,th) result(psi)
    class(wrf_type)     :: this
    real(r8),intent(in) :: th
    real(r8)            :: psi
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function psi_from_th_base
  function dpsidth_from_th_base(this,th) result(dpsidth)
    class(wrf_type)     :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function dpsidth_from_th_base
  function ftc_from_psi_base(this,psi) result(ftc)
    class(wkf_type)     :: this
    real(r8),intent(in) :: psi
    real(r8)            :: ftc
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function ftc_from_psi_base
  function dftcdpsi_from_psi_base(this,psi) result(dftcdpsi)
    class(wkf_type)     :: this
    real(r8),intent(in) :: psi
    real(r8)            :: dftcdpsi
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function dftcdpsi_from_psi_base

  ! =====================================================================================
  ! Van Genuchten Functions are defined here
  ! =====================================================================================

  subroutine set_wrf_param_vg(this,params_in)

    class(wrf_type_vg) :: this
    real(r8), intent(in) :: params_in(:)

    this%alpha  = params_in(1)
    this%psd    = params_in(2)
    this%th_sat = params_in(3)
    this%th_res = params_in(4)

    return
  end subroutine set_wrf_param_vg

  ! =====================================================================================

  subroutine set_wkf_param_vg(this,params_in)

    class(wkf_type_vg) :: this
    real(r8), intent(in) :: params_in(:)

    this%alpha  = params_in(1)
    this%psd    = params_in(2)
    this%th_sat = params_in(3)
    this%th_res = params_in(4)
    this%tort   = params_in(5)

    return
  end subroutine set_wkf_param_vg

  ! =====================================================================================

  function th_from_psi_vg(this,psi) result(th)

    ! Van Genuchten (1980) calculation of volumetric water content (theta)
    ! from matric potential.

    class(wrf_type_vg)   :: this
    real(r8), intent(in) :: psi           ! Matric potential [MPa]
    real(r8)             :: satfrac       ! Saturated fraction [-]
    real(r8)             :: th            ! Volumetric Water Cont [m3/m3]

    real(r8)             :: psi_interp      ! psi where we start lin interp [Mpa]
    real(r8)             :: th_interp       ! th where we start lin interp
    real(r8)             :: dpsidth_interp  ! change in psi during lin interp (slope)
    real(r8)             :: m               ! pore size distribution param (1/n)

    m   = 1._r8/this%psd

    ! pressure above which we use a linear function
    psi_interp = -(1._r8/this%alpha)*(max_sf_interp**(1._r8/(m-1._r8)) - 1._r8 )**m

    if(psi<psi_interp) then

       ! Saturation fraction
       satfrac = (1._r8 + (-this%alpha*psi)**this%psd)**(-1._r8+m)

       ! convert to volumetric water content
       th = satfrac*(this%th_sat-this%th_res) + this%th_res

    else
       th_interp = max_sf_interp * (this%th_sat-this%th_res) + this%th_res
       dpsidth_interp = this%dpsidth_from_th(th_interp)

       th = th_interp + (psi-psi_interp)/dpsidth_interp

    end if

  end function th_from_psi_vg

  ! =====================================================================================

  function psi_from_th_vg(this,th) result(psi)

    ! Van Genuchten (1980) calculation of matric potential from
    ! volumetric water content (theta).

    class(wrf_type_vg)   :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi            ! matric potential [MPa]
    real(r8)             :: m              ! inverse of psd
    real(r8)             :: satfrac        ! saturated fraction
    real(r8)             :: th_interp      ! theta where we start interpolation
    real(r8)             :: psi_interp     ! psi at interpolation point
    real(r8)             :: dpsidth_interp

    !------------------------------------------------------------------------------------
    ! saturation fraction is the origial equation in vg 1980, we just
    ! need to invert it:
    ! (note "psd" is the pore-size-distribution parameter, equivalent to "n" from the
    ! manuscript.)
    !
    ! satfrac = (1._r8 + (alpha*psi)**psd)**(1._r8/psd-1)
    !
    ! *also modified to accomodate linear pressure regime for super-saturation
    ! -----------------------------------------------------------------------------------

    m   = 1._r8/this%psd
    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)

    if(satfrac>=max_sf_interp) then

       th_interp = max_sf_interp * (this%th_sat-this%th_res) + this%th_res
       dpsidth_interp = this%dpsidth_from_th(th_interp)
       psi_interp = -(1._r8/this%alpha)*(max_sf_interp**(1._r8/(m-1._r8)) - 1._r8 )**m
       psi = psi_interp + dpsidth_interp*(th-th_interp)

    else
       
       psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m


    end if

  end function psi_from_th_vg

  ! =====================================================================================

  function dpsidth_from_th_vg(this,th) result(dpsidth)

    class(wrf_type_vg)  :: this
    real(r8),intent(in) :: th           ! water content
    real(r8)            :: a1           ! parameter intermediary
    real(r8)            :: m1           ! parameter intermediary
    real(r8)            :: m2           ! parameter intermediary
    real(r8)            :: satfrac      ! saturation fraction
    real(r8)            :: dsatfrac_dth ! deriv satfrac wrt theta
    real(r8)            :: dpsidth      ! change in matric potential WRT VWC
    real(r8)            :: th_interp    ! vwc where we start interpolation range

    a1 = 1._r8/this%alpha
    m1 = 1._r8/this%psd
    m2 = 1._r8/(m1-1._r8)

    th_interp = max_sf_interp * (this%th_sat-this%th_res) + this%th_res

    ! Since we apply linear interpolation beyond the max and min saturated fractions
    ! we just cap satfrac at those values and calculate the derivative there
    !!    satfrac = max(min(max_sf_interp,(th-this%th_res)/(this%th_sat-this%th_res)),min_sf_interp)

    if(th>th_interp) then
       satfrac = max_sf_interp
    else
       satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
    end if

    dsatfrac_dth = 1._r8/(this%th_sat-this%th_res)

    ! psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m
    ! psi = -a1 * (satfrac**m2 - 1)** m1
    ! dpsi dth = -(m1)*a1*(satfrac**m2-1)**(m1-1) * m2*(satfrac)**(m2-1)*dsatfracdth

    ! f(x) = satfrac**m2 -1
    ! g(x) = a1*f(x)**m1
    ! dpsidth = g'(f(x)) f'(x)

    dpsidth = -m1*a1*(satfrac**m2 - 1._r8)**(m1-1._r8) * m2*satfrac**(m2-1._r8)*dsatfrac_dth


  end function dpsidth_from_th_vg

  ! =====================================================================================

  function ftc_from_psi_vg(this,psi) result(ftc)

    class(wkf_type_vg)  :: this
    real(r8),intent(in) :: psi
    real(r8)            :: num ! numerator term
    real(r8)            :: den ! denominator term
    real(r8)            :: ftc
    real(r8)            :: psi_eff
    real(r8)            :: m          ! inverse pore size distribution param (1/psd)

    m   = 1._r8/this%psd
    
    if(psi<0._r8) then

       ! VG 1980 assumes a postive pressure convention...
       psi_eff = -psi

       num = (1._r8 - (this%alpha*psi_eff)**(this%psd-1._r8) * &
            (1._r8 + (this%alpha*psi_eff)**this%psd)**(-(1._r8-m)))**2._r8
       den = (1._r8 + (this%alpha*psi_eff)**this%psd)**(this%tort*(1._r8-m))

       ! Make sure this is well behaved
       ftc = min(1._r8,max(min_ftc,num/den))

    else
       ftc = 1._r8

    end if

  end function ftc_from_psi_vg

  ! ====================================================================================

  function dftcdpsi_from_psi_vg(this,psi) result(dftcdpsi)

    ! The derivative of the fraction of total conductivity
    ! Note, this function is fairly complex. To get the derivative
    ! we brake it into terms, see the technical note.

    class(wkf_type_vg) :: this
    real(r8),intent(in) :: psi
    real(r8) :: psi_eff  ! VG 1980 assumed positive convention, so we switch sign
    real(r8) :: t1       ! term 1 in numerator
    real(r8) :: t2       ! term 2 in numerator
    real(r8) :: t3       ! term 3 (denomenator)
    real(r8) :: dt1      ! derivative of term 1
    real(r8) :: dt2      ! derivative of term 2
    real(r8) :: dt3      ! derivative of term 3
    real(r8) :: ftc      ! calculate current ftc to see if we are at min
    real(r8) :: dftcdpsi ! change in frac total cond wrt psi
    real(r8) :: m        ! pore size distribution param (1/psd)
    
    m   = 1._r8/this%psd

    if(psi>=0._r8) then
       dftcdpsi = 0._r8
    else
       psi_eff = -psi  ! switch VG 1980 convention

       ftc = this%ftc_from_psi(psi)

       if(ftc<=min_ftc) then
          dftcdpsi = 0._r8   ! We cap ftc, so derivative is zero
       else

          t1  = (this%alpha*psi_eff)**(this%psd-1._r8)
          dt1 = this%alpha*(this%psd-1._r8)*(this%alpha*psi_eff)**(this%psd-2._r8)

          t2  = (1._r8 + (this%alpha*psi_eff)**this%psd)**(m-1._r8)
          dt2 = (m-1._r8) * &
               (1._r8 + (this%alpha*psi_eff)**this%psd)**(m-2._r8) * &
               this%psd * (this%alpha*psi_eff)**(this%psd-1._r8) * this%alpha

          t3  = (1._r8 + (this%alpha*psi_eff)**this%psd)**(this%tort*( 1._r8-m))
          dt3 = this%tort*(1._r8-m) * &
               (1._r8 + (this%alpha*psi_eff)**this%psd )**(this%tort*(1._r8-m)-1._r8) * &
               this%psd * (this%alpha*psi_eff)**(this%psd-1._r8) * this%alpha

          dftcdpsi = 2._r8*(1._r8-t1*t2)*(t1*dt2 + t2*dt1)/t3 - &
               t3**(-2._r8)*dt3*(1._r8-t1*t2)**2._r8
       end if

    end if

  end function dftcdpsi_from_psi_vg

  ! =====================================================================================
  ! =====================================================================================
  ! Campbell, Clapp-Hornberger Water Retention Functions
  ! =====================================================================================
  ! =====================================================================================

  subroutine set_wrf_param_cch(this,params_in)

    class(wrf_type_cch) :: this
    real(r8), intent(in) :: params_in(:)
    real(r8) :: th_max

    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)

    ! Set DERIVED constants
    ! used for interpolating in extreme ranges
    th_max           = max_sf_interp*this%th_sat-1.e-9_r8
    this%psi_max     = this%psi_from_th(th_max)
    this%dpsidth_max = this%dpsidth_from_th(th_max)

    return
  end subroutine set_wrf_param_cch

  ! =====================================================================================

  subroutine set_wkf_param_cch(this,params_in)

    class(wkf_type_cch) :: this
    real(r8), intent(in) :: params_in(:)

    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)
    return
  end subroutine set_wkf_param_cch

  ! =====================================================================================

  function th_from_psi_cch(this,psi) result(th)

    class(wrf_type_cch)  :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th
    real(r8)             :: satfrac

    if(psi>this%psi_max) then
        ! Linear range for extreme values
        th = max_sf_interp*this%th_sat + (psi-this%psi_max)/this%dpsidth_max
    else
        th = this%th_sat*(psi/this%psi_sat)**(-1.0_r8/this%beta)
    end if
    return
  end function th_from_psi_cch

  ! =====================================================================================

  function psi_from_th_cch(this,th) result(psi)

    class(wrf_type_cch)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi
    real(r8)             :: satfrac

    satfrac = th/this%th_sat
    if(satfrac>max_sf_interp) then
        psi = this%psi_max + this%dpsidth_max*(th-max_sf_interp*this%th_sat)
    else
        psi = this%psi_sat*(th/this%th_sat)**(-this%beta)
    end if

  end function psi_from_th_cch

  ! =====================================================================================

  function dpsidth_from_th_cch(this,th) result(dpsidth)

    class(wrf_type_cch)  :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth

    ! Differentiate:
    ! psi = this%psi_sat*(th/this%th_sat)**(-this%beta)

    dpsidth = -this%beta*this%psi_sat/this%th_sat * (th/this%th_sat)**(-this%beta-1._r8)


  end function dpsidth_from_th_cch

  ! =====================================================================================

  function ftc_from_psi_cch(this,psi) result(ftc)

    class(wkf_type_cch) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: psi_eff
    real(r8)            :: ftc

    ! ftc = (th/th_sat)**(2*b+3)
    !     = (th_sat*(psi/psi_sat)**(-1/b)/th_sat)**(2*b+3)
    !     = ((psi/psi_sat)**(-1/b))**(2*b+3)
    !     = (psi/psi_sat)**(-2-3/b)


    psi_eff = min(psi,this%psi_sat)

    ftc = (psi_eff/this%psi_sat)**(-2._r8-3._r8/this%beta)

  end function ftc_from_psi_cch

  ! ====================================================================================

  function dftcdpsi_from_psi_cch(this,psi) result(dftcdpsi)

    class(wkf_type_cch) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: dftcdpsi ! change in frac total cond wrt psi

    ! Differentiate:
    ! ftc = (psi/this%psi_sat)**(-2._r8-3._r8/this%beta)

    ! Note that if we assume a constant, capped FTC=1.0
    ! at saturation, then the derivative is zero there
    if(psi<this%psi_sat)then
       dftcdpsi = (-2._r8-3._r8/this%beta) / this%psi_sat * &
            (psi/this%psi_sat)**(-3._r8-3._r8/this%beta)
    else
       dftcdpsi = 0._r8
    end if


  end function dftcdpsi_from_psi_cch


  ! =====================================================================================
  ! TFS style functions
  ! =====================================================================================

  subroutine set_wkf_param_tfs(this,params_in)

    class(wkf_type_tfs)  :: this
    real(r8), intent(in) :: params_in(:)

    this%p50    = params_in(1)
    this%avuln  = params_in(2)

    return
  end subroutine set_wkf_param_tfs

  subroutine set_wrf_param_tfs(this,params_in)

    class(wrf_type_tfs)  :: this
    real(r8), intent(in) :: params_in(:)
    real(r8) :: th_max
    real(r8) :: th_min

    this%th_sat   = params_in(1)
    this%th_res   = params_in(2)
    this%pinot    = params_in(3)
    this%epsil    = params_in(4)
    this%rwc_ft   = params_in(5)
    this%cap_corr = params_in(6)
    this%cap_int  = params_in(7)
    this%cap_slp  = params_in(8)
    this%pmedia   = int(params_in(9))

    ! Set DERIVED constants
    ! used for interpolating in extreme ranges
    th_max=max_sf_interp*(this%th_sat-this%th_res)+this%th_res-1.e-9_r8
    th_min=min_sf_interp*(this%th_sat-this%th_res)+this%th_res+1.e-9_r8
    this%psi_max     = this%psi_from_th(th_max)
    this%dpsidth_max = this%dpsidth_from_th(th_max)
    this%psi_min     = this%psi_from_th(th_min)
    this%dpsidth_min = this%dpsidth_from_th(th_min)


    return
  end subroutine set_wrf_param_tfs

  ! =====================================================================================

  function th_from_psi_tfs(this,psi) result(th)

    class(wrf_type_tfs)  :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th

    ! !LOCAL VARIABLES:
    real(r8) :: lower                ! lower bound of initial estimate         [m3 m-3]
    real(r8) :: upper                ! upper bound of initial estimate         [m3 m-3]

    real(r8) :: satfrac              ! soil saturation fraction                [0-1]
    real(r8) :: psi_check

    !----------------------------------------------------------------------


    if(psi>this%psi_max) then

        ! Linear range for extreme values
        th = this%th_res+max_sf_interp*(this%th_sat-this%th_res) + &
              (psi-this%psi_max)/this%dpsidth_max

    elseif(psi<this%psi_min) then

        ! Linear range for extreme values
        th = this%th_res+min_sf_interp*(this%th_sat-this%th_res) + &
              (psi-this%psi_min)/this%dpsidth_min

    else

       ! The bisection scheme performs a search via method of bisection,
       ! we need to define bounds with which to start
       lower  = this%th_res+min_sf_interp*(this%th_sat-this%th_res)-1.e-9_r8
       upper  = this%th_res+max_sf_interp*(this%th_sat-this%th_res)+1.e-9_r8

       call this%bisect_pv(lower, upper, psi, th)
       psi_check = this%psi_from_th(th)
       if( psi_check > -1.e-8_r8) then
          write(fates_log(),*)'bisect_pv returned positive value for water potential?'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

   end if
   return
  end function th_from_psi_tfs

  ! =====================================================================================

  function psi_from_th_tfs(this,th) result(psi)

    class(wrf_type_tfs)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi

    ! locals
    real(r8) :: th_corr        ! corrected vol wc [m3/m3]
    real(r8) :: psi_sol        ! pressure from solute term
    real(r8) :: psi_press      ! pressure from "pressure term"
    real(r8) :: psi_elastic    ! press from elastic
    real(r8) :: psi_capillary  ! press from capillary
    real(r8) :: psi_capelast   ! press from smoothed capillary/elastic
    real(r8) :: psi_cavitation ! press from cavitation
    real(r8) :: b,c            ! quadratic smoothing terms
    real(r8) :: satfrac        ! saturated fraction (between res and sat)

    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)

    if(satfrac>max_sf_interp) then

       psi = this%psi_max + this%dpsidth_max * &
            (th-(max_sf_interp*(this%th_sat-this%th_res)+this%th_res))
       
    elseif(satfrac<min_sf_interp) then
       
       psi = this%psi_min + this%dpsidth_min * &
              (th-(min_sf_interp*(this%th_sat-this%th_res)+this%th_res))
       
    else
       
       th_corr = th * this%cap_corr
       
       ! Perform two rounds of quadratic smoothing, 1st smooth
       ! the elastic and capilary, and then smooth their
       ! combined with the caviation
       
       call solutepsi(th_corr,this%rwc_ft,this%th_sat,this%th_res,this%pinot,psi_sol)
       call pressurepsi(th_corr,this%rwc_ft,this%th_sat,this%th_res,this%pinot,this%epsil,psi_press)
       
       psi_elastic = psi_sol + psi_press
       
       if(this%pmedia == 1) then            ! leaves have no capillary region in their PV curves
          
          psi_capelast = psi_elastic
          
       else if(this%pmedia <= 4) then       ! sapwood has a capillary region
          
          call capillarypsi(th_corr,this%th_sat,this%cap_int,this%cap_slp,psi_capillary)
          
          b = -1._r8*(psi_capillary + psi_elastic)
          c = psi_capillary*psi_elastic
          psi_capelast = (-b - sqrt(b*b - 4._r8*quad_a1*c))/(2._r8*quad_a1)
          
       else
          write(fates_log(),*) 'TFS WRF was called for an inelligable porous media'
          call endrun(msg=errMsg(sourcefile, __LINE__))
          
       end if !porous media
       
       ! Now lets smooth the result of capilary elastic with cavitation
       
       psi_cavitation = psi_sol
       b = -1._r8*(psi_capelast + psi_cavitation)
       c = psi_capelast*psi_cavitation
       
       psi = (-b + sqrt(b*b - 4._r8*quad_a2*c))/(2._r8*quad_a2)
    end if
 
    return
  end function psi_from_th_tfs

  ! =====================================================================================

  function dpsidth_from_th_tfs(this,th) result(dpsidth)

    class(wrf_type_tfs) :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth


    ! locals
    real(r8) :: th_corr        ! corrected vol wc [m3/m3]
    real(r8) :: satfrac        ! saturated fraction (between res and sat)
    real(r8) :: psi_sol
    real(r8) :: psi_press
    real(r8) :: psi_elastic    ! press from elastic
    real(r8) :: psi_capillary  ! press from capillary
    real(r8) :: psi_capelast   ! press from smoothed capillary/elastic
    real(r8) :: psi_cavitation ! press from cavitation
    real(r8) :: b,c            ! quadratic smoothing terms
    real(r8) :: dbdth,dcdth    ! derivs of quad smoohting terms
    real(r8) :: dsol_dth
    real(r8) :: dpress_dth
    real(r8) :: delast_dth
    real(r8) :: dcap_dth
    real(r8) :: dcapelast_dth
    real(r8) :: dcav_dth

    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
    if(satfrac>max_sf_interp) then

        dpsidth = this%dpsidth_max

    elseif(satfrac<min_sf_interp) then

        dpsidth = this%dpsidth_min

    else
       th_corr = th*this%cap_corr
       
       ! Perform two rounds of quadratic smoothing, 1st smooth
       ! the elastic and capilary, and then smooth their
       ! combined with the caviation
       
       call solutepsi(th_corr,this%rwc_ft,this%th_sat,this%th_res,this%pinot,psi_sol)
       call pressurepsi(th_corr,this%rwc_ft,this%th_sat,this%th_res,this%pinot,this%epsil,psi_press)
       
       call dsolutepsidth(th,this%th_sat,this%th_res,this%rwc_ft,this%pinot,dsol_dth)
       call dpressurepsidth(this%th_sat,this%th_res,this%rwc_ft,this%epsil,dpress_dth)
       
       delast_dth = dsol_dth + dpress_dth
       psi_elastic = psi_sol + psi_press
       
       
       if(this%pmedia == 1) then         ! leaves have no capillary region in their PV curves
          
          psi_capelast = psi_elastic
          dcapelast_dth = delast_dth
          
       else if(this%pmedia <= 4) then    ! sapwood has a capillary region
          
          call capillarypsi(th,this%th_sat,this%cap_int,this%cap_slp,psi_capillary)
          
          b = -1._r8*(psi_capillary + psi_elastic)
          c = psi_capillary*psi_elastic
          psi_capelast = (-b - sqrt(b*b - 4._r8*quad_a1*c))/(2._r8*quad_a1)
          
          call dcapillarypsidth(this%cap_slp,this%th_sat,dcap_dth)
          
          dbdth = -1._r8*(delast_dth + dcap_dth)
          dcdth = psi_elastic*dcap_dth + delast_dth*psi_capillary
          
          
          dcapelast_dth = 1._r8/(2._r8*quad_a1) * &
               (-dbdth - 0.5_r8*((b*b - 4._r8*quad_a1*c)**(-0.5_r8)) * &
               (2._r8*b*dbdth - 4._r8*quad_a1*dcdth))
          
       else
          write(fates_log(),*) 'TFS WRF was called for an ineligible porous media'
          call endrun(msg=errMsg(sourcefile, __LINE__))
          
       end if !porous media
       
       ! Now lets smooth the result of capilary elastic with cavitation
       
       psi_cavitation = psi_sol
       
       b = -1._r8*(psi_capelast + psi_cavitation)
       c = psi_capelast*psi_cavitation
       
       dcav_dth = dsol_dth
       
       dbdth = -1._r8*(dcapelast_dth + dcav_dth)
       dcdth = psi_capelast*dcav_dth + dcapelast_dth*psi_cavitation
       
       dpsidth = 1._r8/(2._r8*quad_a2)*(-dbdth + 0.5_r8*((b*b - 4._r8*quad_a2*c)**(-0.5_r8)) * &
            (2._r8*b*dbdth - 4._r8*quad_a2*dcdth))
    end if

    return
  end function dpsidth_from_th_tfs

  ! =====================================================================================

  function ftc_from_psi_tfs(this,psi) result(ftc)

    class(wkf_type_tfs) :: this
    real(r8),intent(in) :: psi   !
    real(r8)            :: ftc
    real(r8)            :: psi_eff

    psi_eff = min(0._r8,psi)

    ftc = max(min_ftc,1._r8/(1._r8 + (psi_eff/this%p50)**this%avuln))

  end function ftc_from_psi_tfs

  ! ====================================================================================

  function dftcdpsi_from_psi_tfs(this,psi) result(dftcdpsi)

    class(wkf_type_tfs) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: ftc      ! differentiated func
    real(r8)            :: fx       ! portion of ftc function
    real(r8)            :: dfx      ! differentiation of portion of func
    real(r8)            :: dftcdpsi ! change in frac total cond wrt psi

    ! Differentiate
    ! ftc = 1._r8/(1._r8 + (psi/this%p50(ft))**this%avuln(ft))

    if(psi>0._r8)then
       dftcdpsi = 0._r8
    else
       ftc = 1._r8/(1._r8 + (psi/this%p50)**this%avuln)
       if(ftc<min_ftc) then
          dftcdpsi = 0._r8
       else
          fx  = 1._r8 + (psi/this%p50)**this%avuln
          dfx = this%avuln*(psi/this%p50)**(this%avuln-1._r8) * (1._r8/this%p50)
          dftcdpsi = -fx**(-2._r8)*dfx
       end if
    end if

  end function dftcdpsi_from_psi_tfs

  ! =====================================================================================

  subroutine solutepsi(th,rwc_ft,th_sat,th_res,pinot,psi)
    !
    ! !DESCRIPTION: computes solute water potential (negative) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS

    real(r8)      , intent(in)     :: th          ! vol wc       [m3 m-3]
    real(r8)      , intent(in)     :: rwc_ft
    real(r8)      , intent(in)     :: th_sat
    real(r8)      , intent(in)     :: th_res
    real(r8)      , intent(in)     :: pinot
    real(r8)      , intent(out)    :: psi         ! water potential   [MPa]

    ! -----------------------------------------------------------------------------------
    ! From eq 8, Christopherson et al:
    !
    ! psi = pino/RWC*, where RWC*=(rwc-rwc_res)/(rwc_ft-rwc_res)
    ! psi = pino * (rwc_ft-rwc_res)/(rwc-rwc_res)
    !
    ! if rwc_res =  th_res/th_sat
    !
    !     = pino * (rwc_ft - th_res/th_sat)/(th/th_sat - th_res/th_sat )
    !     = pino * (th_sat*rwc_ft - th_res)/(th - th_res)
    ! -----------------------------------------------------------------------------------
    
    psi = pinot * (th_sat*rwc_ft - th_res) / (th - th_res)

    return
  end subroutine solutepsi

  ! ====================================================================================

  subroutine dsolutepsidth(th,th_sat,th_res,rwc_ft,pinot,dpsi_dth)

    !
    ! !DESCRIPTION: returns derivative of solutepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)     :: th
    real(r8)      , intent(in)     :: th_sat
    real(r8)      , intent(in)     :: th_res
    real(r8)      , intent(in)     :: rwc_ft
    real(r8)      , intent(in)     :: pinot
    real(r8)      , intent(out)    :: dpsi_dth

    ! -----------------------------------------------------------------------------------
    ! Take derivative of eq 8 (function solutepsi)
    ! psi      =  pinot * (th_sat*rwc_ft - th_res) * (th - th_res)^-1
    ! dpsi_dth = -pinot * (th_sat*rwc_ft - th_res) * (th - th_res)^-2
    ! -----------------------------------------------------------------------------------
    
    dpsi_dth = -1._r8*pinot*(th_sat*rwc_ft - th_res )*(th - th_res)**(-2._r8)

    return
  end subroutine dsolutepsidth

  ! =====================================================================================

  subroutine pressurepsi(th,rwc_ft,th_sat,th_res,pinot,epsil,psi)
    !
    ! !DESCRIPTION: computes pressure water potential (positive) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8) , intent(in)  :: th
    real(r8) , intent(in)  :: rwc_ft
    real(r8) , intent(in)  :: th_sat
    real(r8) , intent(in)  :: th_res
    real(r8) , intent(in)  :: pinot
    real(r8) , intent(in)  :: epsil
    real(r8) , intent(out) :: psi         ! water potential   [MPa]

    psi = epsil * (th - th_sat*rwc_ft) / (th_sat*rwc_ft-th_res) - pinot

    return
  end subroutine pressurepsi


  ! =====================================================================================

  subroutine dpressurepsidth(th_sat,th_res,rwc_ft,epsil,dpsi_dth)
    !
    ! !DESCRIPTION: returns derivative of pressurepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)     :: th_sat
    real(r8)      , intent(in)     :: th_res
    real(r8)      , intent(in)     :: rwc_ft
    real(r8)      , intent(in)     :: epsil
    real(r8)      , intent(out)    :: dpsi_dth       ! derivative of water potential wrt theta  [MPa m3 m-3]

    dpsi_dth = epsil/(th_sat*rwc_ft - th_res)

    return
  end subroutine dpressurepsidth

  ! =====================================================================================

  subroutine capillarypsi(th,th_sat,cap_int,cap_slp,psi)
    !
    ! !DESCRIPTION: computes water potential in the capillary region of the plant
    !  PV curve (sapwood only)
    !
    ! !ARGUMENTS

    real(r8)      , intent(in)     :: th          ! water content     [m3 m-3]
    real(r8)      , intent(in)     :: th_sat
    real(r8)      , intent(in)     :: cap_int
    real(r8)      , intent(in)     :: cap_slp
    real(r8)      , intent(out)    :: psi         ! water potential   [MPa]

    psi = cap_int + th*cap_slp/th_sat

    return
  end subroutine capillarypsi

  ! =====================================================================================

  subroutine dcapillarypsidth(cap_slp,th_sat,y)
    !
    ! !DESCRIPTION: returns derivative of capillaryPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)       , intent(in)     :: cap_slp
    real(r8)       , intent(in)     :: th_sat
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]

    y = cap_slp/th_sat

    return
  end subroutine dcapillarypsidth

  ! =====================================================================================
  
  subroutine bisect_pv(this,lower, upper, psi, th)
    !
    ! !DESCRIPTION: Bisection routine for getting the inverse of the plant PV curve.
    !  An analytical solution is not possible because quadratic smoothing functions
    !  are used to remove discontinuities in the PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    class(wrf_type_tfs)  :: this
    real(r8)      , intent(inout)  :: lower       ! lower bound of estimate           [m3 m-3]
    real(r8)      , intent(inout)  :: upper       ! upper bound of estimate           [m3 m-3]
    real(r8)      , intent(in)     :: psi         ! water potential                   [MPa]
    real(r8)      , intent(out)    :: th          ! water content                     [m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x_new                  ! new estimate for x in bisection routine
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    integer  :: nitr                   ! number of iterations

    !----------------------------------------------------------------------
    real(r8),parameter :: xtol = 1.e-16_r8     ! error tolerance for th          [m3 m-3]
    real(r8),parameter :: ytol = 1.e-8_r8      ! error tolerance for psi         [MPa]

    if(psi  > 0.0_r8) then
        write(fates_log(),*)'Error: psi become positive during pv bisection'
        write(fates_log(),*)'psi: ',psi
      call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    y_lo = this%psi_from_th(lower)
    y_hi = this%psi_from_th(upper)

    f_lo  = y_lo - psi
    f_hi  = y_hi - psi
    chg   = upper - lower

    nitr = 0
    do while(abs(chg) .gt. xtol .and. nitr < 100)
       x_new = 0.5_r8*(lower + upper)
       y_new = this%psi_from_th(x_new)
       f_new = y_new - psi
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
       nitr = nitr + 1
    end do

    if(nitr .eq. 100)then
        write(fates_log(),*)'Warning: number of iteraction reaches 100 for bisect_pv'
    endif

    th = x_new

  end subroutine bisect_pv


end module FatesHydroWTFMod
