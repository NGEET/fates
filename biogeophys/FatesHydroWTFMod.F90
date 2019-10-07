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


  real(r8), parameter :: min_ftc = 0.005_r8


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
  ! Plant-only fractional loss of conductivity from Chrisoffersen et al. (tfs model)
  ! =====================================================================================

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
!  procedure :: th_from_psi     => th_from_psi
!  procedure :: psi_from_th     => psi_from_th
!  procedure :: dpsidth_from_th => dpsidth_from_th
!  procedure :: set_wrf_param   => set_wrf_param
! 
!  procedure :: set_wkf_param     => set_wkf_param


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

    !satfrac = (1._r8/(1._r8 + (alpha*abs(psi))**n))**m
    ! Saturation fraction
    ! 

    satfrac = (1._r8 + (-this%alpha*psi)**this%psd)**(-1._r8+1._r8/this%psd)
    
    ! convert to volumetric water content
    th = satfrac*(this%th_sat-this%th_res) + this%th_res

  end function th_from_psi_vg

  ! =====================================================================================
  
  function psi_from_th_vg(this,th) result(psi)
    
    ! Van Genuchten (1980) calculation of matric potential from
    ! volumetric water content (theta).

    class(wrf_type_vg)   :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi
    real(r8)             :: m          ! inverse of psd
    real(r8)             :: satfrac    ! saturated fraction
   
    !------------------------------------------------------------------------------------
    ! saturation fraction is the origial equation in vg 1980, we just
    ! need to invert it:
    ! satfrac = (1._r8 + (alpha*psi)**n)**(1._r8/n-1)
    ! -----------------------------------------------------------------------------------
    
    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)

    m   = 1._r8/this%psd
    psi = -(1._r8/this%alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m 

  end function psi_from_th_vg

  ! =====================================================================================

  function dpsidth_from_th_vg(this,th) result(dpsidth)

    class(wrf_type_vg)  :: this
    real(r8),intent(in) :: th
    real(r8)            :: a1       ! parameter intermediary
    real(r8)            :: m1       ! parameter intermediary
    real(r8)            :: m2       ! parameter intermediary
    real(r8)            :: satfrac  ! saturation fraction
    real(r8)            :: dpsidth  ! change in matric potential WRT VWC

    a1 = 1._r8/this%alpha
    m1 = 1._r8/this%psd
    m2 = 1._r8/(m1-1._r8)

    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)

    ! psi = -a1*(satfrac**m2 - 1._r8 )**m1
    ! f(x) = satfrac**m2 -1 
    ! g(x) = a1*f(x)**m1
    ! dpsidth = g'(f(x)) f'(x)

    dpsidth = -(m2/(this%th_sat - this%th_res))*m1*a1*(satfrac**m2 - 1._r8)**(m1-1._r8)

  end function dpsidth_from_th_vg

  ! =====================================================================================

  function ftc_from_psi_vg(this,psi) result(ftc)

    class(wkf_type_vg)  :: this
    real(r8),intent(in) :: psi
    real(r8)            :: num ! numerator term
    real(r8)            :: den ! denominator term
    real(r8)            :: ftc

    num = (1._r8 - (-this%alpha*psi)**(this%psd-1._r8) * & 
         (1._r8 + (-this%alpha*psi)**this%psd)**(-(1._r8-1._r8/this%psd)))**2._r8
    den = (1._r8 + (-this%alpha*psi)**this%psd)**(this%tort*(1._r8-1._r8/this%psd))
    
    ftc = num/den

  end function ftc_from_psi_vg

  ! ====================================================================================

  function dftcdpsi_from_psi_vg(this,psi) result(dftcdpsi)

    ! The derivative of the fraction of total conductivity
    ! Note, this function is fairly complex. To get the derivative
    ! we brake it into terms, and also into numerator and denominator
    ! and then differentiate those by parts
    class(wkf_type_vg) :: this
    real(r8),intent(in) :: psi
    real(r8) :: t1       ! term 1 in numerator
    real(r8) :: t2       ! term 2 in numerator
    real(r8) :: dt1      ! derivative of term 1
    real(r8) :: dt2      ! derivative of term 2
    real(r8) :: num      ! numerator
    real(r8) :: dnum     ! derivative of numerator
    real(r8) :: den      ! denominator
    real(r8) :: dden     ! derivative of denominator
    real(r8) :: dftcdpsi ! change in frac total cond wrt psi
    
    t1  = (-this%alpha*psi)**(this%psd-1._r8)
    dt1 = this%alpha**(this%psd-1._r8)*(this%psd-1._r8)*psi**(this%psd-2._r8)

    t2  = (1._r8 + (-this%alpha*psi)**this%psd)**(-1._r8+1._r8/this%psd)
    dt2 = -(1._r8-1._r8/this%psd) * & 
         (1._r8 + (-this%alpha*psi)**this%psd)**(1._r8/this%psd) * & 
         this%psd*(this%alpha**this%psd)*(-psi)**(this%psd-1._r8)
   
    num  = (1._r8 - t1*t2)**2._r8
    dnum = 2._r8 * (1._r8 - t1*t2) * ( t1*dt2 + t2*dt1 )

    den  = (1._r8 + (-this%alpha*psi)**this%psd)**(this%tort*( 1._r8-1._r8/this%psd))
    dden = (this%tort*( 1._r8-1._r8/this%psd)) * & 
          (1._r8 + (-this%alpha*psi)**this%psd)**(this%tort*( 1._r8-1._r8/this%psd)-1._r8) * & 
          this%alpha**this%psd * this%psd * (-psi)**(this%psd-1._r8)


    dftcdpsi = dnum*den**(-1._r8) - (den**(-2._r8))*dden*num 

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
    
    class(wrf_type_cch)   :: this
    real(r8), intent(in) :: psi
    real(r8)             :: th
    real(r8)             :: satfrac

    satfrac = (psi/this%psi_sat)**(-1.0_r8/this%beta)

    th = satfrac*this%th_sat

  end function th_from_psi_cch

  ! =====================================================================================

  function psi_from_th_cch(this,th) result(psi)
    
    class(wrf_type_cch)   :: this
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

    ! th = this%th_sat*(psi/this%psi_sat)**(-1.0_r8/this%beta)
    ! ftc = ((psi/this%psi_sat)**(-1.0_r8/this%beta))**(2._r8*this%beta+3._r8)
    ! 
    ! Prevent super-saturation from generating unreasonable FTCs

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
  ! Fractional loss of conductivity via TFS style functions
  ! =====================================================================================

  subroutine set_wkf_param_tfs(this,params_in)
    
    class(wkf_type_tfs)  :: this
    real(r8), intent(in) :: params_in(:)
    
    this%th_sat = params_in(1)
    this%p50    = params_in(2)
    this%avuln  = params_in(3)
    
    return
  end subroutine set_wkf_param_tfs

  ! =====================================================================================

  function ftc_from_psi_tfs(this,psi) result(ftc)
    
    class(wkf_type_tfs) :: this
    real(r8),intent(in) :: psi   ! 
    real(r8)            :: ftc
    
    ftc = max(min_ftc,1._r8/(1._r8 + (psi/this%p50)**this%avuln))

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

    ftc = 1._r8/(1._r8 + (psi/this%p50)**this%avuln)
    if(ftc<min_ftc) then
       dftcdpsi = 0._r8
    else
       fx  = 1._r8 + (psi/this%p50)**this%avuln
       dfx = this%avuln*(psi/this%p50)**(this%avuln-1._r8)
       dftcdpsi = -fx**(-2._r8)*dfx
    end if

  end function dftcdpsi_from_psi_tfs
  

  
end module FatesHydroWTFMod
