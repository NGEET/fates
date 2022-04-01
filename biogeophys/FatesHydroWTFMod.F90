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


  real(r8), parameter :: min_ftc = 0.00001_r8   ! Minimum allowed fraction of total conductance

  ! Bounds on saturated fraction, outside of which we use linear PV or stop flow
  ! In this context, the saturated fraction is defined by the volumetric WC "th"
  ! and the volumetric residual and saturation "th_res" and "th_sat": (th-th_r)/(th_sat-th_res)

  real(r8), parameter :: min_sf_interp = 0.01  ! Linear interpolation below this saturated frac
  real(r8), parameter :: max_sf_interp = 0.998 ! Linear interpolation above this saturated frac

  real(r8), parameter :: quad_a1 = 0.80_r8  ! smoothing factor "A" term
                                            ! in the capillary-elastic region

  real(r8), parameter :: quad_a2 = 0.99_r8  ! Smoothing factor or "A" term in
                                            ! elastic-caviation region


  ! Generic class that can be extended to describe
  ! specific water retention functions

  type, public :: wrf_type

      ! These min and maxes are used only to make the PV functions well
      ! behaved around the endpoints. By doing so we allow the functions
      ! to enter a linear range above the max and below the min, which
      ! should be very close to residual and saturation respectively.
      ! The expectation is that the solvers should never step deeply
      ! into these linear regions, and they only exist reall to handle
      ! strange cases where the solvers overshoot and predict above and below
      ! saturation and residual respectively.

      real(r8) :: psi_max     ! psi matching max_sf_interp where we start linear interp
      real(r8) :: psi_min     ! psi matching min_sf_interp
      real(r8) :: dpsidth_max ! dpsi_dth where we start linear interp
      real(r8) :: dpsidth_min ! dpsi_dth where we start min interp
      real(r8) :: th_min      ! vwc matching min_sf_interp where we start linear interp
      real(r8) :: th_max      ! vwc matching max_sf_interp where we start linear interp

  contains

     procedure :: th_from_psi     => th_from_psi_base
     procedure :: psi_from_th     => psi_from_th_base
     procedure :: dpsidth_from_th => dpsidth_from_th_base
     procedure :: set_wrf_param   => set_wrf_param_base
     procedure :: get_thsat       => get_thsat_base

     ! All brands of WRFs have access to these tools to operate
     ! above and below sat and residual, should they want to
     procedure, non_overridable :: psi_linear_sat
     procedure, non_overridable :: psi_linear_res
     procedure, non_overridable :: th_linear_sat
     procedure, non_overridable :: th_linear_res
     procedure, non_overridable :: set_min_max
     procedure, non_overridable :: get_thmin

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
     real(r8) :: n_vg    ! pore size distribution parameter, psd in original code
     real(r8) :: m_vg    ! m in van Genuchten 1980, also a pore size distribtion parameter , 1-m in original code
     real(r8) :: th_sat  ! Saturation volumetric water content [m3/m3]
     real(r8) :: th_res  ! Residual volumetric water content   [m3/m3]
   contains
     procedure :: th_from_psi     => th_from_psi_vg
     procedure :: psi_from_th     => psi_from_th_vg
     procedure :: dpsidth_from_th => dpsidth_from_th_vg
     procedure :: set_wrf_param   => set_wrf_param_vg
     procedure :: get_thsat       => get_thsat_vg
  end type wrf_type_vg

  ! Water Conductivity Function
  type, public, extends(wkf_type) :: wkf_type_vg
     real(r8) :: alpha   ! Inverse air entry parameter         [m3/Mpa]
     real(r8) :: n_vg    ! pore size distribution parameter
     real(r8) :: m_vg    ! m in van Genuchten 1980, also a pore size distribtion parameter
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
   contains
     procedure :: th_from_psi     => th_from_psi_cch
     procedure :: psi_from_th     => psi_from_th_cch
     procedure :: dpsidth_from_th => dpsidth_from_th_cch
     procedure :: set_wrf_param   => set_wrf_param_cch
     procedure :: get_thsat       => get_thsat_cch
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
  ! Type1 Smooth approximation of Clapp-Hornberger and Campbell (CCH) water retention and conductivity functions
  ! Bisht et al. Geosci. Model Dev., 11, 4085–4102, 2018
  ! =====================================================================================

  ! Water Retention Function
  type, public, extends(wrf_type) :: wrf_type_smooth_cch
     real(r8) :: th_sat   ! Saturation volumetric water content         [m3/m3]
     real(r8) :: psi_sat  ! Bubbling pressure (potential at saturation) [Mpa]
     real(r8) :: beta     ! Clapp-Hornberger "beta" parameter           [-]
     real(r8) :: scch_pu  ! An estimated breakpoint capillary pressure, below which the specified water retention curve is applied. It is also the lower limit when the smoothing function is applied. [Mpa]
     real(r8) :: scch_ps  ! An estimated breakpoint capillary pressure, an upper limit where smoothing funciton is applied. [Mpa]
     real(r8) :: scch_b2  ! constant coefficient of the quadratic term in the smoothing polynomial function [-]
     real(r8) :: scch_b3  ! constant coefficient of the cubic term in the smoothing polynomial function [-]
   contains
     procedure :: th_from_psi     => th_from_psi_smooth_cch
     procedure :: psi_from_th     => psi_from_th_smooth_cch
     procedure :: dpsidth_from_th => dpsidth_from_th_smooth_cch
     procedure :: set_wrf_param   => set_wrf_param_smooth_cch
     procedure :: get_thsat       => get_thsat_smooth_cch
  end type wrf_type_smooth_cch

  ! Water Conductivity Function
  type, public, extends(wkf_type) :: wkf_type_smooth_cch
     real(r8) :: th_sat   ! Saturation volumetric water content         [m3/m3]
     real(r8) :: psi_sat  ! Bubbling pressure (potential at saturation) [Mpa]
     real(r8) :: beta     ! Clapp-Hornberger "beta" parameter           [-]
     real(r8) :: scch_pu  ! An estimated breakpoint capillary pressure, below which the specified water retention curve is applied. It is also the lower limit when the smoothing function is applied. [Mpa]
     real(r8) :: scch_ps  ! An estimated breakpoint capillary pressure, an upper limit where smoothing funciton is applied. [Mpa]
     real(r8) :: scch_b2  ! constant coefficient of the quadratic term in the smoothing polynomial function [-]
     real(r8) :: scch_b3  ! constant coefficient of the cubic term in the smoothing polynomial function [-]
   contains
     procedure :: ftc_from_psi      => ftc_from_psi_smooth_cch
     procedure :: dftcdpsi_from_psi => dftcdpsi_from_psi_smooth_cch
     procedure :: set_wkf_param     => set_wkf_param_smooth_cch
  end type wkf_type_smooth_cch

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
   contains
     procedure :: th_from_psi     => th_from_psi_tfs
     procedure :: psi_from_th     => psi_from_th_tfs
     procedure :: dpsidth_from_th => dpsidth_from_th_tfs
     procedure :: set_wrf_param   => set_wrf_param_tfs
     procedure :: get_thsat       => get_thsat_tfs
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

  ! Generic Functions usable by all
  ! Note that these are linear extrapolations.
  ! They should only be used with the expectation that they will allow
  ! for solutions outside the expected range, with the understanding these
  ! are temporary pertubations, probably through fluctuations in precision
  ! of numerical integration.
  ! ============================================================================

  subroutine set_min_max(this,th_res,th_sat)

      ! This routine uses max_sf_interp and min_sft_interp
      ! to define the bounds of where the linear ranges start and stop

      class(wrf_type)   :: this
      real(r8),intent(in) :: th_res
      real(r8),intent(in) :: th_sat

      this%th_max      = max_sf_interp*(th_sat-th_res)+th_res
      this%th_min      = min_sf_interp*(th_sat-th_res)+th_res
      this%psi_max     = this%psi_from_th(this%th_max-tiny(this%th_max))
      this%dpsidth_max = this%dpsidth_from_th(this%th_max-tiny(this%th_max))
      this%psi_min     = this%psi_from_th(this%th_min+tiny(this%th_min))
      this%dpsidth_min = this%dpsidth_from_th(this%th_min+tiny(this%th_min))

  end subroutine set_min_max

  ! ============================================================================

  function psi_linear_res(this,th) result(psi)

      ! Calculate psi in linear range below residual

      class(wrf_type)   :: this
      real(r8),intent(in)  :: th    ! vol. wat. cont   [m3/m3]
      real(r8)             :: psi   ! Matric potential [MPa]

      psi = this%psi_min + this%dpsidth_min * (th-this%th_min)

  end function psi_linear_res

  ! ===========================================================================

  function get_thmin(this) result(th_min)

    class(wrf_type)   :: this
    real(r8) :: th_min
    
    th_min = this%th_min
    
  end function get_thmin

  ! ===========================================================================
  
  function psi_linear_sat(this,th) result(psi)

      ! Calculate psi in linear range above saturation

      class(wrf_type)   :: this
      real(r8),intent(in)  :: th    ! vol. wat. cont   [m3/m3]
      real(r8)             :: psi   ! Matric potential [MPa]

      psi = this%psi_max + this%dpsidth_max * (th-this%th_max)

  end function psi_linear_sat

  ! ===========================================================================

  function th_linear_sat(this,psi) result(th)

      ! Calculate th from psi in linear range above saturation

      class(wrf_type)   :: this
      real(r8),intent(in)   :: psi   ! Matric potential [MPa]
      real(r8)              :: th    ! vol. wat. cont   [m3/m3]

      th = this%th_max + (psi-this%psi_max)/this%dpsidth_max
  end function th_linear_sat

  ! ===========================================================================

  function th_linear_res(this,psi) result(th)

      ! Calculate th from psi in linear range above saturation

      class(wrf_type)   :: this
      real(r8),intent(in)   :: psi   ! Matric potential [MPa]
      real(r8)              :: th    ! vol. wat. cont   [m3/m3]

      th = this%th_min + (psi-this%psi_min)/this%dpsidth_min

  end function th_linear_res

  ! ===========================================================================

  subroutine set_wrf_param_base(this,params_in)
    class(wrf_type)     :: this
    real(r8),intent(in) :: params_in(:)
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end subroutine set_wrf_param_base
  function get_thsat_base(this) result(th_sat)
    class(wrf_type)     :: this
    real(r8)            :: th_sat
    th_sat = 0._r8
    write(fates_log(),*) 'The base thsat call'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function get_thsat_base
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
    th = 0._r8
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function th_from_psi_base
  function psi_from_th_base(this,th) result(psi)
    class(wrf_type)     :: this
    real(r8),intent(in) :: th
    real(r8)            :: psi
    psi = 0._r8
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function psi_from_th_base
  function dpsidth_from_th_base(this,th) result(dpsidth)
    class(wrf_type)     :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth
    dpsidth = 0._r8
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function dpsidth_from_th_base
  function ftc_from_psi_base(this,psi) result(ftc)
    class(wkf_type)     :: this
    real(r8),intent(in) :: psi
    real(r8)            :: ftc
    ftc = 0._r8
    write(fates_log(),*) 'The base water retention function'
    write(fates_log(),*) 'should never be actualized'
    write(fates_log(),*) 'check how the class pointer was setup'
    call endrun(msg=errMsg(sourcefile, __LINE__))
  end function ftc_from_psi_base
  function dftcdpsi_from_psi_base(this,psi) result(dftcdpsi)
    class(wkf_type)     :: this
    real(r8),intent(in) :: psi
    real(r8)            :: dftcdpsi
    dftcdpsi = 0._r8
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
    this%n_vg   = params_in(2)
    this%m_vg   = params_in(3)    
    this%th_sat = params_in(4)
    this%th_res = params_in(5)

    call this%set_min_max(this%th_res,this%th_sat)

    return
  end subroutine set_wrf_param_vg

  ! =====================================================================================

  subroutine set_wkf_param_vg(this,params_in)

    class(wkf_type_vg) :: this
    real(r8), intent(in) :: params_in(:)

    this%alpha  = params_in(1)
    this%n_vg   = params_in(2)
    this%m_vg   = params_in(3)
    this%th_sat = params_in(4)
    this%th_res = params_in(5)
    this%tort   = params_in(6)

    return
  end subroutine set_wkf_param_vg

  ! =====================================================================================

  function get_thsat_vg(this) result(th_sat)
      class(wrf_type_vg)   :: this
      real(r8) :: th_sat

      th_sat = this%th_sat

  end function get_thsat_vg

  ! =====================================================================================


  
  ! =====================================================================================

  function th_from_psi_vg(this,psi) result(th)

    ! Van Genuchten (1980) calculation of volumetric water content (theta)
    ! from matric potential.

    class(wrf_type_vg)   :: this
    real(r8), intent(in) :: psi            ! Matric potential [MPa]
    real(r8)             :: satfrac        ! Saturated fraction [-]
    real(r8)             :: th             ! Volumetric Water Cont [m3/m3]
    real(r8)             :: dpsidth_interp ! change in psi during lin interp (slope)
    real(r8)             :: m              ! pore size distribution param 1
    real(r8)             :: n              ! pore size distribution param 2

    m = this%m_vg
    n = this%n_vg
    
    if(psi>this%psi_max) then

        ! Linear range for extreme values
        th = this%th_linear_sat(psi)

    elseif(psi<this%psi_min) then

        ! Linear range for extreme values
        th = this%th_linear_res(psi)

    else

       ! Saturation fraction
       satfrac = (1._r8 + (-this%alpha*psi)**n)**(-m)
       
       ! convert to volumetric water content
       th = satfrac*(this%th_sat-this%th_res) + this%th_res


    end if

  end function th_from_psi_vg

  ! =====================================================================================

  function psi_from_th_vg(this,th) result(psi)

    ! Van Genuchten (1980) calculation of matric potential from
    ! volumetric water content (theta).

    class(wrf_type_vg)   :: this
    real(r8), intent(in)  :: th
    real(r8)             :: psi            ! matric potential [MPa]
    real(r8)             :: m
    real(r8)             :: n              
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
    ! satfrac = (1._r8 + (alpha*psi)**psd)**(-m)
    !
    ! *also modified to accomodate linear pressure regime for super-saturation
    ! -----------------------------------------------------------------------------------

    if(th>this%th_max)then

        psi = this%psi_linear_sat(th)

    elseif(th<this%th_min)then

        psi = this%psi_linear_res(th)

    else
       m = this%m_vg
       n = this%n_vg
       satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
       psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(-m)) - 1._r8 )**(1/n)

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
    m1 = 1._r8/this%n_vg
    m2 = -1._r8/this%m_vg

    if(th > this%th_max) then

        dpsidth = this%dpsidth_max

    elseif(th<this%th_min) then

        dpsidth = this%dpsidth_min

    else

        satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
        dsatfrac_dth = 1._r8/(this%th_sat-this%th_res)

        ! psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(-m)) - 1._r8 )**(1/n)
        ! psi = -a1 * (satfrac**m2 - 1)** m1
        ! dpsi dth = -(m1)*a1*(satfrac**m2-1)**(m1-1) * m2*(satfrac)**(m2-1)*dsatfracdth
        ! f(x) = satfrac**m2 -1
        ! g(x) = a1*f(x)**m1
        ! dpsidth = g'(f(x)) f'(x)
        dpsidth = -m1*a1*(satfrac**m2 - 1._r8)**(m1-1._r8) * m2*satfrac**(m2-1._r8)*dsatfrac_dth
        
    end if

  end function dpsidth_from_th_vg

  ! =====================================================================================

  function ftc_from_psi_vg(this,psi) result(ftc)

    class(wkf_type_vg)  :: this
    real(r8),intent(in) :: psi
    real(r8)            :: num ! numerator term
    real(r8)            :: den ! denominator term
    real(r8)            :: ftc
    real(r8)            :: psi_eff
    real(r8)            :: m          ! pore size distribution param ()
    real(r8)            :: n          ! pore size distribution param (psd)

    n   = this%n_vg
    m   = this%m_vg

    if(psi<0._r8) then

       ! VG 1980 assumes a postive pressure convention...
       psi_eff = -psi
       num = (1._r8 - ((this%alpha*psi_eff)**(n) / &
            (1._r8 + (this%alpha*psi_eff)**n))**m)**2._r8
       den = (1._r8 + (this%alpha*psi_eff)**n)**(this%tort*(m))

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
    real(r8) :: n        

    n   =this%n_vg
    m   =this%m_vg

    if(psi>=0._r8) then
       dftcdpsi = 0._r8
    else
       psi_eff = -psi  ! switch VG 1980 convention

       ftc = this%ftc_from_psi(psi)

       if(ftc<=min_ftc) then
          dftcdpsi = 0._r8   ! We cap ftc, so derivative is zero
       else

          t1  = (this%alpha*psi_eff)**(n*m)
          dt1 = this%alpha*(n*m)*(this%alpha*psi_eff)**(n*m-1._r8)

          t2  = (1._r8 + (this%alpha*psi_eff)**n)**(-m)
          dt2 = (-m) * &
               (1._r8 + (this%alpha*psi_eff)**n)**(-m-1._r8) * &
               n * (this%alpha*psi_eff)**(n-1._r8) * this%alpha

          t3  = (1._r8 + (this%alpha*psi_eff)**n)**(this%tort*(m))
          dt3 = this%tort*(m) * &
               (1._r8 + (this%alpha*psi_eff)**n )**(this%tort*(m)-1._r8) * &
               n * (this%alpha*psi_eff)**(n-1._r8) * this%alpha
          
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
    real(r8) :: th_max ! saturated water content

    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)

    ! Set DERIVED constants
    ! used for interpolating in extreme ranges
    this%th_max      = max_sf_interp*this%th_sat
    this%psi_max     = this%psi_from_th(this%th_max-tiny(this%th_max))
    this%dpsidth_max = this%dpsidth_from_th(this%th_max-tiny(this%th_max))
    this%th_min      = fates_unset_r8
    this%psi_min     = fates_unset_r8
    this%dpsidth_min = fates_unset_r8

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

  function get_thsat_cch(this) result(th_sat)
      class(wrf_type_cch)   :: this
      real(r8) :: th_sat

      th_sat = this%th_sat

  end function get_thsat_cch

  ! =====================================================================================

  function th_from_psi_cch(this,psi) result(th)

    class(wrf_type_cch)  :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th

    if(psi>this%psi_max) then
        ! Linear range for extreme values
        th = this%th_max + (psi-this%psi_max)/this%dpsidth_max
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

    if(th>this%th_max) then
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
    if(th>this%th_max) then
        dpsidth = this%dpsidth_max
    else
        dpsidth = -this%beta*this%psi_sat/this%th_sat * (th/this%th_sat)**(-this%beta-1._r8)
    end if

  end function dpsidth_from_th_cch

  ! =====================================================================================

  function ftc_from_psi_cch(this,psi) result(ftc)

    class(wkf_type_cch) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: psi_eff
    real(r8)            :: ftc

    ! th = th_sat * (psi/psi_sat)^(-1/b)

    ! ftc = (th/th_sat)^(2*b+3)
    ! ftc = ( th_sat * (psi/psi_sat)^(-1/b) / th_sat) ^(2*b+3)
    !     = ((psi/psi_sat)^(-1/b))^(2*b+3)
    !     = (psi/psi_sat)^(-2-3/b)


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
  ! =====================================================================================
  ! Type1 Smooth approximation Campbell, Clapp-Hornberger Water Retention Functions
  ! =====================================================================================
  ! =====================================================================================
  subroutine set_wrf_param_smooth_cch(this,params_in)

    class(wrf_type_smooth_cch) :: this
    real(r8), intent(in) :: params_in(:)
    integer  :: styp    ! an option to force constant coefficient of the quadratic
                        ! term 0 (styp = 1) or to force the constant coefficient of
                        ! the cubic term 0 (styp/=2)
    real(r8) :: th_max  ! saturated water content [-]

    ! !LOCAL VARIABLES:
    real(r8) :: pu                ! an estimated breakpoint at which the constant
                                  ! coefficient of the quadratic term (styp=2)
                                  ! or the cubic term (styp/=2) is 0 [Mpa]
    real(r8) :: bcAtPu            ! working local
    real(r8) :: lambdaDeltaPuOnPu ! working local
    real(r8) :: oneOnDeltaPu      ! working local
    real(r8) :: lambda            ! working local, inverse of Clapp and Hornberger "b"
    real(r8) :: alpha             ! working local 
    real(r8) :: ps                ! working local, 90% of entry pressure [Mpa]
    


    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)
    styp = int(params_in(4))


    alpha = -1._r8/this%psi_sat 
    lambda = 1.0_r8/this%beta
    ps = -0.9_r8/alpha
    this%scch_ps     = ps
    ! Choose `pu` that forces `scch_b2 = 0`.
    if(styp == 1) then
      pu               = findGu_SBC_zeroCoeff(lambda, 3, -alpha*ps) / (-alpha)
      this%scch_pu     = pu

      ! Find helper constants.
      bcAtPu            = (-alpha*pu)**(-lambda)
      lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
      oneOnDeltaPu      = 1.d0 / (pu - ps)

      ! Store coefficients for cubic function.
      this%scch_b2 = 0.d0
      this%scch_b3 = (2.d0 - bcAtPu*(2.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu * oneOnDeltaPu
      if( this%scch_b3 <= 0.d0 ) then
         write(fates_log(),*) 'set_wrf_param_smooth_cch b3 <=0',pu,ps,alpha,lambda,oneOnDeltaPu,lambdaDeltaPuOnPu,bcAtPu,this%psi_sat
         call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
    else
      ! Choose `pu` that forces `sbc_b3 = 0`.
      pu               = findGu_SBC_zeroCoeff(lambda, 2, -alpha*ps) / (-alpha)
      this%scch_pu = pu

      ! Find helper constants.
      bcAtPu            = (-alpha*pu)**(-lambda)
      lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
      oneOnDeltaPu      = 1.d0 / (pu - ps)

      ! Store coefficients for cubic function.
      this%scch_b2 = -(3.d0 - bcAtPu*(3.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu* oneOnDeltaPu
      if( this%scch_b2 >= 0.d0 ) then
         write(fates_log(),*) 'set_wrf_param_smooth_cch b2 <= 0'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      this%scch_b3 = 0.d0
     
    endif
    ! Set DERIVED constants
    ! used for interpolating in extreme ranges
    this%th_max      = max_sf_interp*this%th_sat
    this%psi_max     = this%psi_from_th(this%th_max-tiny(this%th_max))
    this%dpsidth_max = this%dpsidth_from_th(this%th_max-tiny(this%th_max))
    this%th_min      = 1.e-8_r8
    this%psi_min     = fates_unset_r8
    this%dpsidth_min = fates_unset_r8

    return
  end subroutine set_wrf_param_smooth_cch



  ! =====================================================================================

  subroutine set_wkf_param_smooth_cch(this,params_in)

    class(wkf_type_smooth_cch) :: this
    real(r8), intent(in) :: params_in(:)
    integer  :: styp              ! an option to force constant coefficient of the
                                  ! quadratic term 0 (styp = 1) or to force the constant
                                  ! coefficient of the cubic term 0 (styp/=2)
    real(r8) :: pu                ! an estimated breakpoint at which the constant
                                  ! coefficient of the quadratic term (styp=2) or
                                  ! the cubic term (styp/=2) is 0 [Mpa]
    real(r8) :: bcAtPu            ! working local
    real(r8) :: lambdaDeltaPuOnPu ! working local
    real(r8) :: oneOnDeltaPu      ! working local
    real(r8) :: lambda            ! working local
    real(r8) :: alpha             ! working local
    real(r8) :: ps                ! working local, 90% of entry pressure [Mpa]

    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)
    styp = int(params_in(4))


    alpha = -1._r8/this%psi_sat 
    lambda = 1.0_r8/this%beta
    ps = -0.9_r8/alpha
    this%scch_ps     = ps
    ! Choose `pu` that forces `scch_b2 = 0`.
    if(styp == 1) then
      pu               = findGu_SBC_zeroCoeff(lambda, 3, -alpha*ps) / (-alpha)
      this%scch_pu     = pu

      ! Find helper constants.
      bcAtPu            = (-alpha*pu)**(-lambda)
      lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
      oneOnDeltaPu      = 1.d0 / (pu - ps)

      ! Store coefficients for cubic function.
      this%scch_b2 = 0.d0
      this%scch_b3 = (2.d0 - bcAtPu*(2.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu * oneOnDeltaPu
      if( this%scch_b3 <= 0.d0 ) then
         write(fates_log(),*) 'set_wrf_param_smooth_cch b3 <=0',pu,ps,alpha,lambda,oneOnDeltaPu,lambdaDeltaPuOnPu,bcAtPu,this%psi_sat
         call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
    else
      ! Choose `pu` that forces `sbc_b3 = 0`.
      pu               = findGu_SBC_zeroCoeff(lambda, 2, -alpha*ps) / (-alpha)
      this%scch_pu = pu

      ! Find helper constants.
      bcAtPu            = (-alpha*pu)**(-lambda)
      lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
      oneOnDeltaPu      = 1.d0 / (pu - ps)

      ! Store coefficients for cubic function.
      this%scch_b2 = -(3.d0 - bcAtPu*(3.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu* oneOnDeltaPu
      if( this%scch_b2 >= 0.d0 ) then
         write(fates_log(),*) 'set_wrf_param_smooth_cch b2 <= 0'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      endif
      this%scch_b3 = 0.d0
     
    endif
    return
  end subroutine set_wkf_param_smooth_cch

  ! =====================================================================================

  function get_thsat_smooth_cch(this) result(th_sat)
      class(wrf_type_smooth_cch)   :: this
      real(r8) :: th_sat
      
      th_sat = this%th_sat
      
  end function get_thsat_smooth_cch
  
  ! =====================================================================================
  
  function th_from_psi_smooth_cch(this,psi) result(th)

    class(wrf_type_smooth_cch)  :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th

    real(r8)                                  :: alpha
    real(r8)                                  :: lambda
    real(r8)                                  :: sat
    real(r8)                                  :: pc
    real(r8)                                  :: ps
    real(r8)                                  :: deltaPc
    real(r8)                                  :: dSe_dpc

    alpha   = -1._r8/this%psi_sat
    lambda  = 1._r8/this%beta
    pc = psi

    if( pc <= this%scch_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       sat      = (-alpha*pc)**(-lambda)
    elseif( pc < this%scch_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - this%scch_ps
       sat      = 1.d0 + deltaPc*deltaPc*(this%scch_b2 + deltaPc*this%scch_b3)
    else
       ! Saturated regime.
       ! Here, `pc >= ps`.
       sat       = 1.d0
    endif
    th = sat * this%th_sat


    return
  end function th_from_psi_smooth_cch

  ! =====================================================================================

  function psi_from_th_smooth_cch(this,th) result(psi)

    class(wrf_type_smooth_cch)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi

    real(r8)                                  :: sat_res
    real(r8)                                  :: alpha
    real(r8)                                  :: lambda
    real(r8)                                  :: Se
    real(r8)                                  :: sat
    real(r8)                                  :: pc
    real(r8)                                  :: xL
    real(r8)                                  :: xc
    real(r8)                                  :: xR
    real(r8)                                  :: resid
    real(r8)                                  :: dx
    real(r8), parameter                       :: relTol = 1.d-9
    integer                                   :: iter
    
    sat_res = 0._r8
    alpha   = -1._r8/this%psi_sat
    lambda  = 1._r8/this%beta

    sat = max(1.e-6,th/this%th_sat)
    if( sat < 1.d0 ) then
       ! Find the `pc` that satisfies the unmodified Brooks-Corey function.
       Se = sat
       pc = -(Se**(-1.d0/lambda)) / alpha
       if( pc > this%scch_pu ) then
          ! Here, solution is in the cubic smoothing regime.
          if( this%scch_b2 == 0.d0 ) then
             ! Note know `b3 > 0`.
             pc = this%scch_ps - ((1.d0 - Se) / this%scch_b3)**(1.d0/3.d0)
          elseif( this%scch_b3 == 0.d0 ) then
             ! Note know `b2 < 0`.
             pc = this%scch_ps - sqrt((Se - 1.d0) / this%scch_b2)
          else
             ! Here, want to solve general cubic
             ! `1 + b2*x^2 + b3*x^3 = Se`
             ! where `x = pc - pu`.
             ! Write as residual function
             ! `r = x^2 * (b2 + b3*x) + (1 - Se)`.
             ! Perform a Newton-Raphson search on `x`.
             ! Have
             ! `dr/dx = x*(2*b2 + 3*b3*x)`
             ! And Newton-Raphson sets
             ! `x[i+1] = x[i] - r[i]/(dr/dx[i])`.
             ! Note that r{0} = 1 - Se > 0.
             ! Therefore maintain the right bracket as having a positive
             ! residual, and the left bracket as having a negative residual.
             ! Note that it is possible, due to numerical effects with `pc`
             ! very close to `pu`, to get an `xL` with a positive residual.
             ! However, in this case also have `xc` very close to `xL`, and
             ! the Newton-Raphson search will converge after a single step.
             ! Therefore do not insert a special test to catch the case here.
             xL = this%scch_pu - this%scch_ps
             xR = 0.d0
             xc = pc - this%scch_ps
             ! write(unit=*, fmt='("SatFunc_SatToPc_SBC: NR search:",
             ! 6(a,g15.6))')  &
             !     ' pu', scch_pu, ' ps', scch_ps,  &
             !     ' xL', xL, ' xR', xR,  &
             !     ' r{xL}', xL*xL*(scch_b2 + scch_b3*xL) +
             !     1.d0 - Se,  &
             !     ' r{xR}', xR*xR*(scch_b2 + scch_b3*xR) +
             !     1.d0 - Se

             iter = 0
             dx = 1.e20_r8 ! something large
             do while( abs(dx) >= -relTol*this%scch_pu )
                
                ! Here, assume:
                ! + Have a bracket on the root, between `xL` and `xR`.
                ! + The residual `r{xL} < 0` and `r{xR} > 0`.
                ! + Have a current guess `xc` at the root.  However, that guess
                !     might not lie in the bracket.

                iter=iter+1
                
                ! Reset `xc` using bisection if necessary.
                if( xc<=xL .or. xc>=xR ) then
                   ! write(unit=*, fmt='("Bisecting")')
                   xc = xL + 0.5d0*(xR - xL)
                endif

                ! Find NR step.
                dx = this%scch_b3 * xc
                resid = xc*xc*(this%scch_b2 + dx) + 1.d0 - Se
                dx = resid / (xc*(2.d0*this%scch_b2 + 3.d0*dx))

                ! Update bracket.
                if( resid > 0.d0 ) then
                   xR = xc
                else
                   xL = xc
                endif

                ! Take the Newton-Raphson step.
                xc = xc - dx
                ! write(unit=*, fmt='(6(a,g15.6))')  &
                !     ' xL', xL, ' xc', xc, ' xR', xR,  &
                !     ' r{xL}', xL*xL*(scch_b2 + scch_b3*xL) +
                !     1.d0 - Se,  &
                !     ' r{xc}', xc*xc*(scch_b2 + scch_b3*xc) +
                !     1.d0 - Se,  &
                !     ' r{xR}', xR*xR*(scch_b2 + scch_b3*xR) +
                !     1.d0 - Se

                if( iter>10000) then
                   write(fates_log(),*) "psi_from_th_smooth_cch iteration not converging"
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
             enddo

             ! Here, have `xc = pc - ps`.
             pc = xc + this%scch_ps
          endif
       endif
    else
       pc = 0.d0
    endif
    psi = pc


  end function psi_from_th_smooth_cch

  ! =====================================================================================

  function dpsidth_from_th_smooth_cch(this,th) result(dpsidth)

    class(wrf_type_smooth_cch)  :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth


    real(r8)                      :: pc
    real(r8)                      :: sat
    real(r8)                      :: dsat_dP
    !
    ! !LOCAL VARIABLES:
    real(r8)                      :: sat_res
    real(r8)                      :: alpha
    real(r8)                      :: lambda
    real(r8)                      :: Se
    real(r8)                      :: deltaPc
    real(r8)                      :: dSe_dpc

    sat_res = 0._r8
    alpha   = -1._r8/this%psi_sat
    lambda  = 1._r8/this%beta

    pc = 1._r8 * this%psi_from_th(th)
    if( pc <= this%scch_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       Se      = (-alpha*pc)**(-lambda)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = -lambda*Se/pc
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
       dpsidth = 1._r8/(dsat_dp * this%th_sat)
    elseif( pc < this%scch_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - this%scch_ps
       Se      = 1.d0 + deltaPc*deltaPc*(this%scch_b2 + deltaPc*this%scch_b3)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = deltaPc*(2*this%scch_b2 + 3*deltaPc*this%scch_b3)
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
       dpsidth = 1._r8/(dsat_dp * this%th_sat)
    else
       ! Saturated regime.
       ! Here, `pc >= ps`.

       dpsidth = this%dpsidth_max
    endif


  end function dpsidth_from_th_smooth_cch

  ! =====================================================================================

  function ftc_from_psi_smooth_cch(this,psi) result(ftc)

    class(wkf_type_smooth_cch) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: ftc

    real(r8)                     :: pc
    real(r8)                     :: kr
    real(r8)                     :: dkr_dP
    !
    real(r8)                     :: sat_res
    real(r8)                     :: alpha
    real(r8)                     :: lambda
    real(r8)                     :: Se
    real(r8)                     :: deltaPc
    real(r8)                     :: dSe_dpc
    real(r8)                     :: dkr_dSe

    pc = psi
    sat_res = 0._r8
    alpha   = -1._r8/this%psi_sat
    lambda  = 1._r8/this%beta

    if( pc <= this%scch_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       Se      = (-alpha*pc)**(-lambda)
       kr      = Se ** (3._r8+2._r8/lambda)

    elseif( pc < this%scch_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - this%scch_ps
       Se      = 1.d0 + deltaPc*deltaPc*(this%scch_b2 + deltaPc*this%scch_b3)
       kr      = Se ** (3._r8+2._r8/lambda)

    else
       ! Saturated regime.
       ! Here, `pc >= ps`.
       kr        = 1.d0
    endif
    ftc = max(kr, min_ftc)


  end function ftc_from_psi_smooth_cch

  ! ====================================================================================

  function dftcdpsi_from_psi_smooth_cch(this,psi) result(dftcdpsi)

    class(wkf_type_smooth_cch) :: this
    real(r8),intent(in) :: psi
    real(r8)            :: dftcdpsi ! change in frac total cond wrt psi

    real(r8)            :: pc
    real(r8)            :: kr
    real(r8)            :: dkr_dP
    !
    real(r8)            :: sat_res
    real(r8)            :: alpha
    real(r8)            :: lambda
    real(r8)            :: Se
    real(r8)            :: deltaPc
    real(r8)            :: dSe_dpc
    real(r8)            :: dkr_dSe

    pc = psi
    sat_res = 0._r8
    alpha   = -1._r8/this%psi_sat
    lambda  = 1._r8/this%beta

    if( pc <= this%scch_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       Se      = (-alpha*pc)**(-lambda)

       dSe_dpc = -lambda*Se/pc

       kr      = Se ** (3.d0 + 2.d0/lambda)

       dkr_dSe = (3.d0 + 2.d0/lambda)*kr/Se
       dkr_dp  = dkr_dSe*dSe_dpc
    elseif( pc < this%scch_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - this%scch_ps
       Se      = 1.d0 + deltaPc*deltaPc*(this%scch_b2 + deltaPc*this%scch_b3)

       dSe_dpc = deltaPc*(2*this%scch_b2 + 3*deltaPc*this%scch_b3)

       kr      = Se ** (2.5d0 + 2.d0/lambda)

       dkr_dSe = (2.5d0 + 2.d0/lambda)*kr/Se
       dkr_dp  = dkr_dSe*dSe_dpc
    else
       ! Saturated regime.
       ! Here, `pc >= ps`.
       kr        = 1.d0
       dkr_dP    = 0.d0
    endif
    dftcdpsi = dkr_dP
    if(kr<=min_ftc) then
          dftcdpsi = 0._r8
    endif


  end function dftcdpsi_from_psi_smooth_cch


  !------------------------------------------------------------------------
  ! Find `pu` that forces a coefficient of the smoothing cubic polynomial to zero.
  ! Bisht et al. Geosci. Model Dev., 11, 4085–4102, 2018, coded in VSFM
  !
  ! Work in terms of multipliers of `pc0`:
  !
  ! + Argument `gs` satisfies `ps = gs*pc0`.
  ! + Return `gu` such that `pu = gu*pc0`.
  !
  ! Argument `AA`:
  !
  ! + To set `b2 = 0`, let `A = 3`.
  ! + To set `b3 = 0`, let `A = 2`.
  !
  real(r8) function findGu_SBC_zeroCoeff(lambda, AA, gs)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    real(r8) , intent(in) :: lambda
    real(r8) , intent(in) :: gs
    integer   , intent(in) :: AA
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: guLeft, gu, guRight  ! Bracketed search.
    real(r8)              :: deltaGu, resid, dr_dGu  ! Newton-Raphson search.
    real(r8)              :: guInv, guToMinusLam, gsOnGu  ! Helper variables.
    real(r8), parameter   :: relTol = 1.d-12
    integer               :: iter
    
    ! Check arguments.
    !   Note this is more for documentation than anything else-- this
    ! fcn should only get used internally, by trusted callers.
    if( lambda<=0.d0 .or. lambda>=2.d0  &
         .or. (AA/=2 .and. AA/=3)  &
         .or. gs>=1.d0 .or. gs<0.d0 ) then
       write(fates_log(),*) 'findGu_SBC_zeroCoeff: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Approach:
    ! + Bracketed Newton-Raphson search.
    ! + Note expect `1 < gu <= gu{gs=0}`.
    ! + Note if this was a critical inner loop, could try solving for
    !     `gui == 1/gu`, rather than for `gu`, in order to avoid division.
    !     Could also try using Ridder's method, since the residual here
    !     has a strong exponential component.

    ! Initialize.
    gu = (AA / (AA + lambda))**(-1.d0/lambda)  ! Solution if `gs = 0`.

    ! Search for root, using bracketed Newton-Raphson.
    !   Not necessary if `gs = 0`.
    if( gs > 0.d0 ) then
       guLeft = 1.d0
       guRight = gu

       ! Test for convergence.
       !   Note this test implicitly also tests `resid == 0`.
       iter = 0
       deltaGu=1.e20_r8 ! something large
       do while( abs(deltaGu) >= relTol*abs(gu) )

          ! Here, assume:
          ! + Have an bracket on the root, between `guLeft` and `guRight`.
          ! + The derivative `dr/d{gu} > 0`.
          ! + The residual `r{guLeft} < 0`, and `r{guRight} > 0`.
          ! + Have a current guess `gu` at the root.  However, that guess
          !     might not lie in the bracket (and does not at first iteration).

          iter=iter+1
          
          ! Reset `gu` using bisection if necessary.
          if( gu<=guLeft .or. gu>=guRight ) then
             gu = guLeft + 0.5d0*(guRight - guLeft)
          endif

          ! Find residual.
          guInv = 1.d0 / gu
          guToMinusLam = gu**(-lambda)  ! Could also do `guInv**lambda`; not sure if any numerical consequences.
          gsOnGu = gs * guInv
          resid = AA - guToMinusLam*(AA + lambda - lambda*gsOnGu)

          ! Update bracket.
          if( resid < 0.d0 ) then
             guLeft = gu
          else
             guRight = gu
          endif

          ! Find next guess using Newton-Raphson's method.
          dr_dGu = (1.d0 + lambda) * (1.d0 - gsOnGu) + (AA - 1)
          dr_dGu = lambda * guToMinusLam * guInv * dr_dGu
          deltaGu = resid / dr_dGu
          ! write(unit=*, fmt='("findGu_SBC_zeroCoeff, NR step: ", 6(a,g15.6))')  &
          !     'guLeft', guLeft, 'gu', gu, 'guRight', guRight,  &
          !     'deltaGu', deltaGu, 'resid', resid, 'dr_dGu', dr_dGu
          gu = gu - deltaGu

          if( iter>10000) then
             write(fates_log(),*) "findGu_SBC_zeroCoeff iteration not converging"
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          
       enddo
    endif

    ! Finish up.
    !   Note assuming the last Newton-Raphson step landed in the bracket,
    ! and had a smaller residual than either of the bracket points.  This
    ! seems a safe enough assumption, compared to cost of tracking residuals.
    findGu_SBC_zeroCoeff = gu
    if(gu /= gu) then
       write(fates_log(),*)'gu = nan in findGu_SBC_zeroCoeff: ',AA,gs
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

  end function findGu_SBC_zeroCoeff

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

    this%th_sat   = params_in(1)
    this%th_res   = params_in(2)
    this%pinot    = params_in(3)
    this%epsil    = params_in(4)
    this%rwc_ft   = params_in(5)
    this%cap_corr = params_in(6)
    this%cap_int  = params_in(7)
    this%cap_slp  = params_in(8)
    this%pmedia   = int(params_in(9))

    call this%set_min_max(this%th_res,this%th_sat)

    return
  end subroutine set_wrf_param_tfs

  ! =====================================================================================

  function get_thsat_tfs(this) result(th_sat)
      class(wrf_type_tfs)   :: this
      real(r8) :: th_sat

      th_sat = this%th_sat

  end function get_thsat_tfs

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

        th = this%th_linear_sat(psi)

    elseif(psi<this%psi_min) then

        ! Linear range for extreme values
        th = this%th_linear_res(psi)

    else

       ! The bisection scheme performs a search via method of bisection,
       ! we need to define bounds with which to start
       lower  = this%th_min-1.e-9_r8
       upper  = this%th_max+1.e-9_r8

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

    if(th>this%th_max)then

        psi = this%psi_linear_sat(th)

   elseif(th<this%th_min)then

       psi = this%psi_linear_res(th)

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

    if(th > this%th_max) then

        dpsidth = this%dpsidth_max

    elseif(th<this%th_min) then

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

    psi_eff = min(-nearzero,psi)

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
    ! psi = pinot/RWC*, where RWC*=(rwc-rwc_res)/(rwc_ft-rwc_res)
    ! psi = pinot * (rwc_ft-rwc_res)/(rwc-rwc_res)
    !
    ! if rwc_res =  th_res/th_sat
    !
    !     = pinot * (rwc_ft - th_res/th_sat)/(th/th_sat - th_res/th_sat )
    !     = pinot * (th_sat*rwc_ft - th_res)/(th - th_res)
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
