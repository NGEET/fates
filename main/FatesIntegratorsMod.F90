module FatesIntegratorsMod

   use EDTypesMod          , only : ed_site_type
   use EDTypesMod          , only : ed_patch_type
   use EDTypesMod          , only : ed_cohort_type
   use FatesConstantsMod, only    : r8 => fates_r8

   implicit none
   integer, parameter :: max_states = 20
   
   public :: RKF45
   public :: Euler

contains

  subroutine RKF45(DerivFunction,Y,Ymask,dx,x,ccohort,max_err,Yout,l_pass)

      ! ---------------------------------------------------------------------------------
      ! Runge-Kutta-Fehlerg  4/5 order adaptive explicit integration
      ! 
      ! 
      ! ---------------------------------------------------------------------------------

      ! Arguments
      
      real(r8),intent(in), dimension(:)         :: Y        ! dependent variable (array)
      logical,intent(in), dimension(:)          :: Ymask    ! logical mask defining what is on
      real(r8),intent(in)                       :: dx       ! step size of independent variable
      real(r8),intent(in)                       :: x        ! independent variable (time?)
      type(ed_cohort_type),intent(inout),target :: ccohort  ! Cohort derived type
      real(r8),intent(in)                       :: max_err  ! Maximum allowable error (absolute)
      real(r8),intent(inout), dimension(:)      :: Yout     ! The output vector
      logical,intent(out)                       :: l_pass   ! Was this a successfully step?

      ! Locals
      integer                             :: nY       ! size of Y
      real(r8), dimension(max_states)     :: Ytemp    ! scratch space for the dependent variable
      real(r8)                            :: xtemp
      real(r8), dimension(max_states)     :: K0
      real(r8), dimension(max_states)     :: K1
      real(r8), dimension(max_states)     :: K2
      real(r8), dimension(max_states)     :: K3
      real(r8), dimension(max_states)     :: K4
      real(r8), dimension(max_states)     :: K5
      real(r8)                            :: err45    ! Estimated integrator error
      
      real(r8), parameter :: min_step_fraction = 0.25_r8

      real(r8), parameter :: t1   = 1.0/4.0
      real(r8), parameter :: f1_0 = 1.0/4.0

      real(r8), parameter :: t2   = 3.0/8.0
      real(r8), parameter :: f2_0 = 3.0/32.0
      real(r8), parameter :: f2_1 = 9.0/32.0

      real(r8), parameter :: t3   = 12.0/13.0
      real(r8), parameter :: f3_0 = 1932.0/2197.0
      real(r8), parameter :: f3_1 = -7200.0/2197.0
      real(r8), parameter :: f3_2 = 7296.0/2197.0

      real(r8), parameter :: t4   = 1.0
      real(r8), parameter :: f4_0 = 439.0/216.0
      real(r8), parameter :: f4_1 = -8.0
      real(r8), parameter :: f4_2 = 3680.0/513.0
      real(r8), parameter :: f4_3 = -845.0/4104.0

      real(r8), parameter :: t5   = 0.5
      real(r8), parameter :: f5_0 = -8.0/27.0
      real(r8), parameter :: f5_1 = 2.0
      real(r8), parameter :: f5_2 = -3544.0/2565.0
      real(r8), parameter :: f5_3 = 1859.0/4104.0
      real(r8), parameter :: f5_4 = -11.0/40.0

      real(r8), parameter :: y_0 = 25.0/216.0
      real(r8), parameter :: y_2 = 1408.0/2565.0
      real(r8), parameter :: y_3 = 2197.0/4104.0
      real(r8), parameter :: y_4 = -1.0/5.0

      real(r8), parameter :: z_0 = 16.0/135.0
      real(r8), parameter :: z_2 = 6656.0/12825.0
      real(r8), parameter :: z_3 = 28561.0/56430.0
      real(r8), parameter :: z_4 = -9.0/50.0
      real(r8), parameter :: z_5 = 2.0/55.0
      
      ! Input Functional Argument
      interface
         function DerivFunction(Y,Ymask,x,ccohort) result(dYdx)
              use EDTypesMod          , only : ed_site_type
              use EDTypesMod          , only : ed_patch_type
              use EDTypesMod          , only : ed_cohort_type
              use FatesConstantsMod, only    : r8 => fates_r8
              real(r8),intent(in), dimension(:)        :: Y        ! dependent variable (array)
              logical,intent(in), dimension(:)         :: Ymask    ! logical mask defining what is on
              real(r8),intent(in)                      :: x        ! independent variable (time?)
              type(ed_cohort_type),intent(in),target   :: ccohort  ! Cohort derived type
              real(r8),dimension(lbound(Y,dim=1):ubound(Y,dim=1)) :: dYdx     ! Derivative of dependent variable
          end function DerivFunction
       end interface
       
       nY = size(Y,1)

       ! 0th Step
       Ytemp(1:nY) = Y(1:nY)
       xtemp       = x
       K0(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)

       ! 1st Step
       Ytemp(1:nY) = Y(1:nY) + dx * (f1_0*K0(1:nY))
       xtemp       = x + t1*dx
       K1(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)
       
       ! 2nd Step
       Ytemp(1:nY) = Y(1:nY) + dx * ( f2_0*K0(1:nY) + f2_1*K1(1:nY) )
       xtemp       = x + t2*dx
       K2(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)
       
       ! 3rd Step
       Ytemp(1:nY) = Y(1:nY) + dx * ( f3_0*K0(1:nY) + f3_1*K1(1:nY) + &
                                      f3_2*K2(1:nY))
       xtemp       = x + t3*dx
       K3(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)

       ! 4th Step
       Ytemp(1:nY) = Y(1:nY) + dx * ( f4_0*K0(1:nY) + f4_1*K1(1:nY) + &
                                      f4_2*K2(1:nY) + f4_3*K3(1:nY))
       xtemp       = x + t4*dx
       K4(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)
       
       ! 5th Step
       Ytemp(1:nY) = Y(1:nY) + dx * ( f5_0*K0(1:nY) + f5_1*K1(1:nY) + &
                                      f5_2*K2(1:nY) + f5_3*K3(1:nY) + &
                                      f5_4*K4(1:nY))
       xtemp       = x + t5*dx
       K5(1:nY)    = DerivFunction(Ytemp(1:nY),Ymask,xtemp,ccohort)

       
       ! Evaluate error on the 4/5 steps

       ! 4th order
       Ytemp(1:nY) = Y(1:nY) + dx * ( y_0*K0(1:nY) + y_2*K2(1:nY) + &
                                      y_3*K3(1:nY) + y_4*K4(1:nY) )
       ! 5th order
       Yout(1:nY)  = Y(1:nY) + dx * ( z_0*K0(1:nY) + z_2*K2(1:nY) + &
                                      z_3*K3(1:nY) + z_4*K4(1:nY) + &
                                      z_5*K5(1:nY) )
       
       ! Take the maximum absolute error across all variables
       ! To prevent weirdness set a nominal lower bound
       err45 = maxval(abs(Yout(1:nY)-Ytemp(1:nY)))

       ! --------------------------------------------------------------------------------
       ! Evaluate error and either approve/reject step.
       ! 
       ! Update our estimate of the optimal time-step. We won't update
       ! the current time-step based on this, but we will save this info
       ! to help decide the starting sub-step on the next full step
       ! The equations may be so smooth that the error estimate is so low that it creates
       ! an overflow on the divide, set a lower bound based on max_err.
       ! 1e-5, as an error ratio will shorten the timestep to ~5% of original
       ! --------------------------------------------------------------------------------

       ccohort%ode_opt_step = dx * max(min_step_fraction, &
                                       0.840896 * (max_err/ max(err45,0.00001*max_err))**0.25)

       if(err45 > max_err) then
          l_pass                 = .false.
       else
          l_pass                 = .true.
       end if

       return
    end subroutine RKF45

    ! ===================================================================================

    subroutine Euler(DerivFunction,Y,Ymask,dx,x,ccohort,Yout)

      ! ---------------------------------------------------------------------------------
      ! Simple Euler Integration
      ! ---------------------------------------------------------------------------------

      ! Arguments
      
      real(r8),intent(in), dimension(:)         :: Y        ! dependent variable (array)
      logical,intent(in), dimension(:)          :: Ymask    ! logical mask defining what is on
      real(r8),intent(in)                       :: dx       ! step size of independent variable
      real(r8),intent(in)                       :: x        ! independent variable (time?)
      type(ed_cohort_type),intent(inout),target :: ccohort  ! Cohort derived type
      real(r8),intent(inout), dimension(:)      :: Yout     ! The output vector

      ! Locals
      integer                             :: nY       ! size of Y
      real(r8), dimension(max_states)     :: Ytemp    ! scratch space for the dependent variable
      real(r8)                            :: xtemp
      real(r8), dimension(max_states)     :: dYdx
      
      ! Input Functional Argument
      interface
         function DerivFunction(Y,Ymask,x,ccohort) result(dYdx)
              use EDTypesMod          , only : ed_site_type
              use EDTypesMod          , only : ed_patch_type
              use EDTypesMod          , only : ed_cohort_type
              use FatesConstantsMod, only    : r8 => fates_r8
              real(r8),intent(in), dimension(:)      :: Y        ! dependent variable (array)
              logical,intent(in), dimension(:)       :: Ymask    ! logical mask defining what is on
              real(r8),intent(in)                    :: x        ! independent variable (time?)
              type(ed_cohort_type),intent(in),target :: ccohort  ! Cohort derived type
              real(r8),dimension(lbound(Y,dim=1):ubound(Y,dim=1)) :: dYdx     ! Derivative of dependent variable
          end function DerivFunction
       end interface

       nY = size(Y,1)
       
       dYdx(1:nY)  = DerivFunction(Y(1:nY),Ymask,x,ccohort)
       Yout(1:nY)  = Y(1:nY) + dx * dYdx(1:nY)
       

       return
    end subroutine Euler


 end module FatesIntegratorsMod
