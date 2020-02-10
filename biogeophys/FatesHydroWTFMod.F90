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


  real(r8), parameter :: min_ftc = 0.0005_r8

  real(r8), parameter :: min_rwc_interp = 0.02
  real(r8), parameter :: max_rwc_interp = 0.95

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

     real(r8), parameter :: beta2 = 0.99_r8  ! Smoothing factor

   contains
     procedure :: th_from_psi     => th_from_psi_tfs
     procedure :: psi_from_th     => psi_from_th_tfs
     procedure :: dpsidth_from_th => dpsidth_from_th_tfs
     procedure :: set_wrf_param   => set_wrf_param_tfs
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

    real(r8)             :: psi_interp    ! psi where we start lin interp
    real(r8)             :: th_interp     ! th where we start lin interp
    real(r8)             :: dpsidth_interp
    real(r8)             :: m

    m   = 1._r8/this%psd
    
    ! pressure above which we use a linear function
    psi_interp = -(1._r8/this%alpha)*(max_rwc_interp**(1._r8/(m-1._r8)) - 1._r8 )**m

    !    psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m 

    
    if(psi<psi_interp) then
       
       ! Saturation fraction
       satfrac = (1._r8 + (-this%alpha*psi)**this%psd)**(-1._r8+1._r8/this%psd)
       
       ! convert to volumetric water content
       th = satfrac*(this%th_sat-this%th_res) + this%th_res

    else
       th_interp = max_rwc_interp * (this%th_sat-this%th_res) + this%th_res
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
    ! satfrac = (1._r8 + (alpha*psi)**n)**(1._r8/n-1)
    ! we also modify these functions to 
    ! -----------------------------------------------------------------------------------

    m   = 1._r8/this%psd
    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
    
    if(satfrac>=max_rwc_interp) then

       th_interp = max_rwc_interp * (this%th_sat-this%th_res) + this%th_res
       dpsidth_interp = this%dpsidth_from_th(th_interp)
       psi_interp = -(1._r8/this%alpha)*(max_rwc_interp**(1._r8/(m-1._r8)) - 1._r8 )**m
       psi = psi_interp + dpsidth_interp*(th-th_interp)

!!    elseif(satfrac<min_rwc_interp .or. .false.) then
       
 !!      th_interp = min_rwc_interp * (this%th_sat-this%th_res) + this%th_res
 !!      dpsidth_interp = this%dpsidth_from_th(th_interp)
!!       psi_interp = -(1._r8/this%alpha)*(min_rwc_interp**(1._r8/(m-1._r8)) - 1._r8 )**m 
!!       psi = psi_interp + dpsidth_interp*(th-th_interp)

    else

    
       psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m 

       
    end if

  end function psi_from_th_vg

  ! =====================================================================================

  function dpsidth_from_th_vg(this,th) result(dpsidth)

    class(wrf_type_vg)  :: this
    real(r8),intent(in) :: th
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

    th_interp = max_rwc_interp * (this%th_sat-this%th_res) + this%th_res
    
    ! Since we apply linear interpolation beyond the max and min saturated fractions
    ! we just cap satfrac at those values and calculate the derivative there
    !!    satfrac = max(min(max_rwc_interp,(th-this%th_res)/(this%th_sat-this%th_res)),min_rwc_interp)

    if(th>th_interp) then
       satfrac = max_rwc_interp
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

    if(psi<0._r8) then
    
       ! VG 1980 assumes a postive pressure convention...
       psi_eff = -psi

       num = (1._r8 - (this%alpha*psi_eff)**(this%psd-1._r8) * & 
            (1._r8 + (this%alpha*psi_eff)**this%psd)**(-(1._r8-1._r8/this%psd)))**2._r8
       den = (1._r8 + (this%alpha*psi_eff)**this%psd)**(this%tort*(1._r8-1._r8/this%psd))

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
       
          t2  = (1._r8 + (this%alpha*psi_eff)**this%psd)**(1._r8/this%psd-1._r8)
          dt2 = (1._r8/this%psd-1._r8) * & 
               (1._r8 + (this%alpha*psi_eff)**this%psd)**(1._r8/this%psd-2._r8) * & 
               this%psd * (this%alpha*psi_eff)**(this%psd-1._r8) * this%alpha
          
          t3  = (1._r8 + (this%alpha*psi_eff)**this%psd)**(this%tort*( 1._r8-1._r8/this%psd))
          dt3 = this%tort*(1._r8-1._r8/this%psd) * & 
               (1._r8 + (this%alpha*psi_eff)**this%psd )**(this%tort*(1._r8-1._r8/this%psd)-1._r8) * & 
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
    
    this%th_sat  = params_in(1)
    this%psi_sat = params_in(2)
    this%beta    = params_in(3)
    
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

    th = this%th_sat*(psi/this%psi_sat)**(-1.0_r8/this%beta)

  end function th_from_psi_cch

  ! =====================================================================================

  function psi_from_th_cch(this,th) result(psi)
    
    class(wrf_type_cch)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi

    psi = this%psi_sat*(th/this%th_sat)**(-this%beta)

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
    
    this%th_sat  = params_in(1)
    this%th_res  = params_in(2)
    this%pinot   = params_in(3)
    this%epsil   = params_in(4)

    return
  end subroutine set_wrf_param_tfs

  ! =====================================================================================

  function th_from_psi_tfs(this,psi) result(th)
    
    class(wrf_type_tfs)  :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th



  end function th_from_psi_tfs

  ! =====================================================================================

  function psi_from_th_tfs(this,th) result(psi)
    
    class(wrf_type_tfs)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi

    x = th_node*cap_corr(pm)
    ! 
    call bq2(x, y_bq2)
    call cq2(x, y_cq2)
    
    psi = (-y_bq2 + sqrt(y_bq2*y_bq2 - 4._r8*beta2*y_cq2))/(2._r8*this%beta2)


    return
  end function psi_from_th_tfs
  
  ! expanded form of psi from th tfs

  function psi_from_th_tfs(this,th) result(psi)
    
    class(wrf_type_tfs)  :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi

    ! locals
    real(r8) :: th_corr        ! corrected vol wc [m3/m3]
    real(r8) :: psi_elastic    ! press from elastic
    real(r8) :: psi_capillary  ! press from capillary
    real(r8) :: psi_capelast   ! press from smoothed capillary/elastic
    real(r8) :: psi_cavitation ! press from cavitation
    real(r8) :: b,c            ! quadratic smoothing terms

    th_corr = th * this%cap_corr

    ! Perform two rounds of quadratic smoothing, 1st smooth
    ! the elastic and capilary, and then smooth their
    ! combined with the caviation

    call elasticPV(th_corr, psi_elastic)
    
    if(this%pmedia == 1) then            ! leaves have no capillary region in their PV curves
        
        psi_capelast = psi_elastic

    else if(pm <= 4) then       ! sapwood has a capillary region
        
        call capillaryPV(th_corr, psi_capillary)
        
        b = -1._r8*(psi_capillary + psi_elastic)
        c = psi_capillary*psi_elastic
        psi_cap_elas = (-b - sqrt(b*b - 4._r8*this%beta*c))/(2._r8*this%beta)

    else
        write(fates_log(),*) 'TFS WRF was called for an inelligable porous media'
        call endrun(msg=errMsg(sourcefile, __LINE__))

    end if !porous media

    ! Now lets smooth the result of capilary elastic with cavitation

    ! call cavitationPV(th_corr, psi_cavitation)
    ! this is caviation PV:
    call solutepsi(th_corr,this%rwcft,this%th_sat,this%th_res,this%pinot,psi_cavitation)

    b = -1._r8*(psi_cap_elas + psi_cavitation)
    c = psi_cap_elas*psi_cavitation

    psi = (-b + sqrt(b*b - 4._r8*this%beta*c))/(2._r8*this%beta)


    return
  end function psi_from_th_tfs

  


  ! =====================================================================================

  function dpsidth_from_th_tfs(this,th) result(dpsidth)

    class(wrf_type_tfs) :: this
    real(r8),intent(in) :: th
    real(r8)            :: dpsidth


    
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
  ! The following routines are for calculating water retention functions
  ! and their derivatives in the TFS model
  ! =====================================================================================

  subroutine bq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = -1._r8*(y_tq1 + y_cavitation)
	     
  end subroutine bq2
  
  !-------------------------------------------------------------------------------!
  subroutine dbq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = -1._r8*(dydth_tq1 + dcavdth)

  end subroutine dbq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine cq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = y_tq1*y_cavitation
	     
  end subroutine cq2
  
  !-------------------------------------------------------------------------------!
  subroutine dcq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq2() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = y_tq1*dcavdth + dydth_tq1*y_cavitation
	     
  end subroutine dcq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine tq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: either calls the elastic region of the PV curve (leaves) or
    ! does a smoothing function for capillary-to-elastic region of the plant PV
    ! curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: y_elastic              ! returned y (psi) value from elasticPV()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call elasticPV(ft, pm, x, y_elastic)
       y = y_elastic
    else if(pm <= 4) then       ! sapwood has a capillary region

        ! bq1
        call capillaryPV(ft, pm, x, y_capillary)
        call elasticPV(ft, pm, x, y_elastic)
       
        b = -1._r8*(y_capillary + y_elastic)
 
        ! cq1
!        call capillaryPV(ft, pm, x, y_capillary)
!        call elasticPV(ft, pm, x, y_elastic)
        c = y_capillary*y_elastic


!       call bq1(ft, pm, x, y_bq1)
!       call cq1(ft, pm, x, y_cq1)
       y = (-b - sqrt(b*b - 4._r8*beta1*c))/(2*b)

    end if !porous media

  end subroutine tq1
  




  !-------------------------------------------------------------------------------!
  subroutine dtq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of tq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: dydth_bq1              ! returned derivative from dbq1dth()
    real(r8) :: dydth_cq1              ! returned derivative from dcq1dth()
    real(r8) :: delasticdth            ! returned derivative from delasticPVdth()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call delasticPVdth(ft, pm, x, delasticdth)
       y = delasticdth
    else if(pm <= 4) then       ! sapwood has a capillary region
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       call dbq1dth(ft, pm, x, dydth_bq1)
       call dcq1dth(ft, pm, x, dydth_cq1)
       y = 1._r8/(2._r8*beta1)*(-dydth_bq1 - 0.5_r8*((y_bq1*y_bq1 - 4._r8*beta1*y_cq1)**(-0.5_r8)) * &
                                                    (2._r8*y_bq1*dydth_bq1 - 4._r8*beta1*dydth_cq1))
    end if

  end subroutine dtq1dth
  
  !-------------------------------------------------------------------------------!

  !-------------------------------------------------------------------------------!
  subroutine dbq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of bq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = -1._r8*(delasticdth + dcapdth)
	     
  end subroutine dbq1dth

  ! ------------------------------------------------------------------------------

  subroutine bq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = -1._r8*(y_capillary + y_elastic)
	     
  end subroutine bq1
  
  !-------------------------------------------------------------------------------!
  subroutine cq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = y_capillary*y_elastic
	     
  end subroutine cq1
  
  !-------------------------------------------------------------------------------!
  subroutine dcq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = y_elastic*dcapdth + delasticdth*y_capillary
	     
  end subroutine dcq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine cavitationPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    !----------------------------------------------------------------------
 
!    call solutepsi(th,rwcft,th_sat,th_res,pinot,psi)
!    call solutepsi(ft, pm, x, y_solute)
    y = y_solute
	     
  end subroutine cavitationPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcavitationPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cavitationPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    y = dsoldth
	     
  end subroutine dcavitationPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine elasticPV(th,rwcft,th_sat,th_res,pinot,psi)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)  :: th
    real(r8)      , intent(in)  :: rwcft
    real(r8)      , intent(in)  :: th_sat
    real(r8)      , intent(in)  :: th_res
    real(r8)      , intent(in)  :: pinot
    real(r8)      , intent(out) :: psi           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    real(r8) :: y_pressure         ! returned y (psi) value from pressurepsi()
    !----------------------------------------------------------------------

    call solutepsi(th,rwcft,th_sat,th_res,pinot,y_solute)

    call pressurepsi(th,rwcft,th_sat,th_res,pinot,epsil,y_pressure)

    psi = y_solute + y_pressure
	     
  end subroutine elasticPV
  
  !-------------------------------------------------------------------------------!
  subroutine delasticPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of elasticPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    real(r8) :: dpressdth         ! returned derivative from dpressurepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    call dpressurepsidth(ft, pm, x, dpressdth)
    y = dsoldth + dpressdth
	     
  end subroutine delasticPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine solutepsi(th,rwcft,th_sat,th_res,pinot,psi)
    ! 
    ! !DESCRIPTION: computes solute water potential (negative) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS

    real(r8)      , intent(in)     :: th          ! vol wc       [m3 m-3]
    real(r8)      , intent(in)     :: rwcft       ! rel wc       [-]
    real(r8)      , intent(in)     :: th_sat
    real(r8)      , intent(in)     :: th_res
    real(r8)      , intent(in)     :: pinot 
    real(r8)      , intent(out)    :: psi         ! water potential   [MPa]

    psi = pinot*th_sat*(rwcft - th_res) / (th - th_sat*th_res)
	     
    end associate

  end subroutine solutepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dsolutepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of solutepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => EDPftvarcon_inst%hydr_pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => EDPftvarcon_inst%hydr_thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node       & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
    
      y = -1._r8*thetas(ft,pm)*pinot(ft,pm)*(rwcft(pm) - resid(ft,pm)) / &
            ((x - thetas(ft,pm)*resid(ft,pm))**2._r8)
	     
    end associate

  end subroutine dsolutepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine pressurepsi(th,rwcft,th_sat,th_res,pinot,epsil,psi)
    ! 
    ! !DESCRIPTION: computes pressure water potential (positive) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8) , intent(in)  :: th
    real(r8) , intent(in)  :: rwcft       ! rel wc       [-]
    real(r8) , intent(in)  :: th_sat
    real(r8) , intent(in)  :: th_res
    real(r8) , intent(in)  :: pinot 
    real(r8) , intent(in)  :: epsil       
    real(r8) , intent(out) :: psi         ! water potential   [MPa]


    psi = epsil * (th - th_sat*rwcft) / &
            (th_sat*(rwcft-th_res)) - pinot
	     
    end associate

  end subroutine pressurepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dpressurepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of pressurepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         thetas  => EDPftvarcon_inst%hydr_thetas_node, & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 resid   => EDPftvarcon_inst%hydr_resid_node , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => EDPftvarcon_inst%hydr_epsil_node   & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
      
      y = epsil(ft,pm)/(thetas(ft,pm)*(rwcft(pm) - resid(ft,pm)))
	     
    end associate

  end subroutine dpressurepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine capillaryPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the capillary region of the plant
    !  PV curve (sapwood only)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         thetas    => EDPftvarcon_inst%hydr_thetas_node     & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
   
      y = cap_int(pm) + cap_slp(pm)/thetas(ft,pm)*x
	     
    end associate

  end subroutine capillaryPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcapillaryPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of capillaryPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         thetas    => EDPftvarcon_inst%hydr_thetas_node    & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
   
    y = cap_slp(pm)/thetas(ft,pm)
	     
    end associate

  end subroutine dcapillaryPVdth
  
end module FatesHydroWTFMod
