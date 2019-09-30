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

  implicit none
  private

  ! -------------------------------------------------------------------------------------
  ! This module contains all unit (F)unctions associated with (W)ater (T)ransfer.
  ! e.g. WTFs
  ! These are also called "pedotransfer" functions, however, since these
  ! may be applied to xylems, stems, etc, they are not limited to soils (pedo).
  ! -------------------------------------------------------------------------------------

  ! Generic class that can be extended to describe
  ! specific Pedo-transfer Functions

  type, public :: wtf_type

     ! These base pointers should never be actually called
     procedure :: th_from_psi     => th_from_psi_base
     procedure :: psi_from_th     => psi_from_th_base
     procedure :: dpsidth_from_th => dpsidth_from_th_base
     procedure :: flc_from_th     => flc_from_th_base
     procedure :: dflc_from_th    => dflcdth_from_th_base
     procedure :: set_param       => set_param_base

  end type wtf_type
  
  ! The Van Genuchten Pedo-transfer functions
  type, public, extends(wtf_type) :: wtf_type_vg

     real(r8), allocatable :: alpha   ! Inverse air entry parameter
     real(r8), allocatable :: psd     ! Inverse width of pore size distribution parameter
     real(r8), allocatable :: th_sat  ! Saturation volumetric water content [m3/m3]
     real(r8), allocatable :: th_res  ! Residual volumetric water content   [m3/m3]
     
   contains

     procedure :: th_from_psi     => th_from_psi_vg
     procedure :: psi_from_th     => psi_from_th_vg
     procedure :: dpsidth_from_th => dpsidth_from_th_vg
     procedure :: flc_from_th     => flc_from_th_vg
     procedure :: dflc_from_th    => dflcdth_from_th_vg
     procedure :: set_param       => set_param_vg
     
  end type wtf_type_vg
  

  ! This object holds all of the functions and
  ! parameters for the different porous media types
  type(pft_type), public, pointer :: wtfs(:)
  

contains

  ! =====================================================================================

  subroutine set_param_vg(this,alpha_in,psd_in,th_sat_in,th_res_in)
    
    class(wtf_type_vg) :: this
    real(r8), intent(in) :: alpha_in
    real(r8), intent(in) :: psd_in
    real(r8), intent(in) :: th_sat_in
    real(r8), intent(in) :: th_res_in

    this%alpha  = alpha_in
    this%psd    = psd_in
    this%th_sat = th_sat_in
    this%th_res = th_res_in
    
    return
  end subroutine set_param_vg
  
  ! =====================================================================================
  
  function th_from_psi_vg(psi) result(th)
    
    ! Van Genuchten (1980) calculation of volumetric water content (theta)
    ! from matric potential.
    
    class(wtf_type_vg)   :: this
    real(r8), intent(in) :: psi
    real(r8)             :: satfrac
    real(r8)             :: th

    !satfrac = (1._r8/(1._r8 + (alpha*abs(psi))**n))**m
    ! Saturation fraction

    satfrac = (1._r8 + (this%alpha*psi)**this%psd)**(-1+(1._r8/this%psd))
    
    ! convert to volumetric water content
    th = satfrac*(this%th_sat-this%th_res) + this%th_res

  end function th_from_psi_vg

  ! =====================================================================================

  function psi_from_th_vg(th) result(psi)
    
    ! Van Genuchten (1980) calculation of matric potential from
    ! volumetric water content (theta).

    class(wtf_type_vg)   :: this
    real(r8),intent(in)  :: th
    real(r8)             :: psi
    real(r8)             :: satfrac
   
    !------------------------------------------------------------------------------------
    ! saturation fraction is the origial equation in vg 1980, we just
    ! need to invert it:
    ! satfrac = (1._r8 + (alpha*psi)**n)**(1._r8/n-1)
    ! -----------------------------------------------------------------------------------
    
    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)
    
    psi = (1._r8/this%alpha)*(satfrac**(1._r8/(1._r8/this%psd-1._r8)) - 1._r8 )**(1._r8/this%psd) 

  end function psi_from_th_vg

  ! =====================================================================================

  function dpsidth_from_th_vg(th) result(dpsidth)

    class(wtf_type_vg)  :: this
    real(r8),intent(in) :: th
    real(r8)            :: a1       ! parameter intermediary
    real(r8)            :: m1       ! parameter intermediary
    real(r8)            :: m2       ! parameter intermediary
    real(r8)            :: satfrac  ! saturation fraction

    a1 = 1._r8/this%alpha
    m1 = 1._r8/this%psd
    m2 = 1._r8/(m1-1._r8)

    satfrac = (th-this%th_res)/(this%th_sat-this%th_res)

    ! psi = a1*(satfrac**m2 - 1._r8 )**m1
    ! f(x) = satfrac**m2 -1 
    ! g(x) = a1*f(x)**m1
    ! dpsidth = g'(f(x)) f'(x)

    dpsidth = (m2/(this%th_sat - this%th_res))*m1*a1*(satfrac**m2 - 1._r8)**(m1-1._r8)

  end function dpsidth_from_th_vg


  ! =====================================================================================
 

  function flc_from_th_vg(th) result(flc)

    num = (1._r8 - (this%alpha*psi)**(this%psd-1._r8) * & 
         (1._r8 + (this%alpha*psi)**this%psd)**(-(1._r8-1._r8/this%psd)))**2._r8
    den = (1._r8 + (this%alpha*psi)**this%psd)**(this%lt*(1._r8-1._r8/this%psd))
    
    flc = num/den
    

  end function flc_from_th_vg

  ! ====================================================================================

  function dflcdpsi_from_psi_vg(th) result(dflcdpsi)
    

    ! The derivative of the fraction of total conductivity
    ! Note, this function is fairly complex. To get the derivative
    ! we brake it into terms, and also into numerator and denominator
    ! and then differentiate those by parts
    class(wtf_type_vg) :: this
    real(r8),intent(in) :: psi
    real(r8) :: t1   ! term 1 in numerator
    real(r8) :: t2   ! term 2 in numerator
    real(r8) :: dt1  ! derivative of term 1
    real(r8) :: dt2  ! derivative of term 2
    real(r8) :: num  ! numerator
    real(r8) :: dnum ! derivative of numerator
    real(r8) :: den  ! denominator
    real(r8) :: dden ! derivative of denominator
    real(r8) :: dflcdpsi ! change in frac total cond wrt psi

    
    t1  = (this%alpha*psi)**(this%psd-1._r8)
    dt1 = this%alpha**(this%psd-1._r8)*(this%psd-1._r8)*psi**(this%psd-2._r8)

    t2  = (1._r8 + (this%alpha*psi)**this%psd)**-(1._r8-1._r8/this%psd)
    dt2 = -(1._r8-1._r8/this%psd) * (1._r8 + (this%alpha*psi)**this%psd)**(1._r8/this%psd) * this%psd*(this%alpha**psd)*psi**(this%psd-1._r8)
   
    num  = (1._r8 - t1*t2)**2._r8
    dnum = 2._r8 * (1._r8 - t1*t2) * ( t1*dt2 + t2*dt1 )

    den  = (1._r8 + (this%alpha*psi)**this%psd)**(this%lt*( 1._r8-1._r8/this%psd))
    dden = (this%lt*( 1._r8-1._r8/this%psd)) * & 
          (1._r8 + (this%alpha*psi)**this%psd)**(this%lt*( 1._r8-1._r8/this%psd)-1._r8) * & 
          this%alpha**this%psd * this%psd * psi**(this%psd-1._r8)


    dflcdpsi = dnum*den**-1 + -(den**-2)*dden*num 



  end function dflcdth_from_th_vg



  


  
end module FatesHydroWTFMod
