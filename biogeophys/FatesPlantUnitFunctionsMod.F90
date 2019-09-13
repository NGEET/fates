module FatesHydroUnitFunctionsMod

  ! This module contains hydraulics functions that are readily broken down into
  ! unit tests.  These are functions that mostly operate on primitive
  ! arguments, are smaller in scope, and are allowed to access the 
  ! parameter constants EDPftvarcon_inst and params

  use FatesConstants, only : fates_unset_r8
  use EDPftvarcon, only : pft_p => EDPftvarcon_inst
  use EDParamsMod      , only : hydr_psi0
  use EDParamsMod      , only : hydr_psicap
  
  implicit none
  private

  logical, parameter :: debug=.true.
  character(len=*), parameter, private :: sourcefile = __FILE__


  integer, parameter :: van_genuchten = 1
  integer, parameter :: campbell      = 2
  integer, parameter :: iswc = campbell


  ! P-V curve: total RWC @ which elastic drainage begins     [-]        
  real(r8), allocatable :: rwcft(:)   !  = (/1.0_r8,0.958_r8,0.958_r8,0.958_r8/)
  
  ! P-V curve: total RWC @ which capillary reserves exhausted
  real(r8), allocatable :: rwccap(:)  !  = (/1.0_r8,0.947_r8,0.947_r8,0.947_r8/) 

  ! P-V curve: slope of capillary region of curve
  real(r8), allocatable :: cap_slp(:)

  ! P-V curve: intercept of capillary region of curve
  real(r8), allocatable :: cap_int(:)

  ! P-V curve: correction for nonzero psi0x
  real(r8), allocatable :: cap_corr(:)


  public :: Hydraulics_Tridiagonal
  public :: flc_gs_from_psi
  public :: dflcgsdpsi_from_psi
  public :: flc_from_psi
  public :: dflcdpsi_from_psi
  public :: th_from_psi
  public :: psi_from_th
  public :: dpsidth_from_th
  public :: bisect_rootfr
  public :: zeng2001_crootfr
  public :: shellGeom
  public :: xylemtaper
  public :: InitAllocatePlantMedia
  public :: SetPlantMediaParam

contains


  ! =====================================================================================

  subroutine InitAllocatePlantMedia(n_plant_media)
    
    ! We only allocate for plant porous media, we do 
    ! not use these arrays to inform on soil relationships
    integer,intent(in) :: n_plant_media

    allocate(rwcft(n_plant_media))
    allocate(rwccap(n_plant_media))
    allocate(cap_slp(n_plant_media))
    allocate(cap_int(n_plant_media))
    allocate(cap_corr(n_plant_media))

    rwcft(:) = fates_unset_r8
    rwcap(:) = fates_unset_r8
    cap_slp(:) = fates_unset_r8
    cap_int(:) = fates_unset_r8
    cap_corr(:) = fates_unset_r8

    return
  end subroutine InitAllocatePlantMedia

  ! =====================================================================================
  
  subroutine SetPlantMediaParam(pm,rwcft_in,rwcap_in)

    ! To avoid complications that would arise from linking this
    ! module with the FatesHydraulicsMemMod.F90 during unit tests, we
    ! store some of these arrays that are indexed by "porous_media"
    ! as globals in this module.
    
    integer,intent(in)  :: pm      ! porous media index
    real(r8),intent(in) :: rwcft_in  ! rwcft for this pm
    real(r8),intent(in) :: rwcap_in  ! rwcap for this pm

    rwcft(pm)  = rwft_in
    rwccap(pm) = rwcap_in
    
    if (pm.eq.1) then   ! Leaf tissue
       cap_slp(pm)    = 0.0_r8
       cap_int(pm)    = 0.0_r8
       cap_corr(pm)   = 1.0_r8
    else               ! Non leaf tissues
       cap_slp(pm)    = (hydr_psi0 - hydr_psicap )/(1.0_r8 - rwccap(pm))  
       cap_int(pm)    = -cap_slp(pm) + hydr_psi0    
       cap_corr(pm)   = -cap_int(pm)/cap_slp(pm)
    end if
    
    return
  end subroutine SetPlantMediaParam
  
  ! =====================================================================================

  subroutine Hydraulics_Tridiagonal(a, b, c, r, u)
    !
    ! !DESCRIPTION: An abbreviated version of biogeophys/TridiagonalMod.F90
    !
    ! This solves the form:
    !
    ! a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = r(i)
    !
    ! It assumed that coefficient a(1) and c(N) DNE as there is
    ! no u(0) or u(N-1).
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8), intent(in)    :: a(:)           ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(:)           ! "b" diagonal column of tridiagonal matrix
    real(r8), intent(in)    :: c(:)           ! "c" right off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: r(:)           ! "r" forcing term of tridiagonal matrix
    real(r8), intent(out)   :: u(:)           ! solution
    !
    ! !LOCAL VARIABLES:
    real(r8) :: bet                           ! temporary
    real(r8) :: gam(n_hypool_tot)             ! temporary
    integer  :: k                             ! index
    real(r8) :: err                           ! solution error, in units of [m3/m3]
    real(r8) :: rel_err                       ! relative error, normalized by delta theta
    real(r8), parameter :: allowable_rel_err = 0.001_r8
    !    real(r8), parameter :: allowable_err = 1.e-6_r8
    !----------------------------------------------------------------------

    bet = b(1)
    do k=1,n_hypool_tot
       if(k == 1) then
          u(k)   = r(k) / bet
       else
          gam(k) = c(k-1) / bet
          bet    = b(k) - a(k) * gam(k)
          u(k)   = (r(k) - a(k)*u(k-1)) / bet
       end if
    enddo

    do k=n_hypool_tot-1,1,-1
       u(k)   = u(k) - gam(k+1) * u(k+1)
    enddo

    ! If debug mode, calculate error on the forward solution
    if(debug)then
       do k=1,n_hypool_tot

          if(k==1)then
             err = abs(r(k) - (b(k)*u(k)+c(k)*u(k+1)))
          elseif(k<n_hypool_tot) then
             err = abs(r(k) - (a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)))
          else
             err = abs(r(k) - (a(k)*u(k-1)+b(k)*u(k)))
          end if

          rel_err = abs(err/u(k))

          if((rel_err > allowable_rel_err)) then !.and. (err > allowable_err) )then
             write(fates_log(),*) 'Tri-diagonal solve produced solution with'
             write(fates_log(),*) 'non-negligable error.'
             write(fates_log(),*) 'Compartment: ',k
             write(fates_log(),*) 'Error in forward solution: ',err
             write(fates_log(),*) 'Estimated delta theta: ',u(k)
             write(fates_log(),*) 'Rel Error: ',rel_err
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

       end do
    end if

  end subroutine Hydraulics_Tridiagonal

  !===============================================================================!

  function flc_gs_from_psi( lwp, ft ) result( btran )

    ! 
    ! !DESCRIPTION: Calculates fractional loss of conductance
    !               across the stomata (gs).

    !
    ! !ARGUMENTS
    real(r8) , intent(in) :: lwp  !leaf water potential (MPa)
    integer  , intent(in) :: ft
    real(r8)              :: btran

    btran = &
         (1._r8 + & 
         (lwp/pft_p%hydr_p50_gs(ft))**pft_p%hydr_avuln_gs(ft))**(-1._r8)

  end function flc_gs_from_psi

  !===============================================================================!
  
  function dflcgsdpsi_from_psi(lwp, ft) result (dflcgsdpsi)

    ! Calculate the derivative of change in fractional loss of conductivity
    ! WRT matric potential.

    ! !ARGUMENTS
    real(r8), intent(in) :: lwp       ! leaf water potential (MPa)
    integer , intent(in) :: ft        ! leaf pft

    real(r8)           :: dflcgsdpsi  ! fractional loss of conductivity  [-]


    associate(avuln_gs => pft_p%hydr_avuln_gs, &  ! Stomatal PLC curve: shape parameter [-]
              p50_gs   => pft_p%hydr_p50_gs)      ! Stomatal PLC curve: water potential
                                                       ! at 50% loss of gs,max  [Pa]

      
      dflcgsdpsi = -1._r8 * (1._r8 + (lwp/p50_gs(ft))**avuln_gs(ft))**(-2._r8) * &
           avuln_gs(ft)/p50_gs(ft)*(lwp/p50_gs(ft))**(avuln_gs(ft)-1._r8)
      
    end associate
    
  end function dflcgsdpsi_from_psi

  !===============================================================================!
  
  function flc_from_psi(ft, pm, psi_node, suc_sat, bsw) result(flc_node)

    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity

    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: psi_node    ! water potential      [MPa]
    real(r8), optional,intent(in)     :: suc_sat     ! minimum soil suction [mm]
    real(r8), optional,intent(in)     :: bsw         ! col Clapp and Hornberger "b"

    real(r8) :: flc_node                             ! frac loss of conductivity [-]

    associate(& 
         avuln    => pft_p%hydr_avuln_node , & ! PLC curve: vulnerability curve shape parameter          [-]
         p50      => pft_p%hydr_p50_node     & ! PLC curve: water potential at 50% loss of conductivity  [Pa]
         )
      
      if(pm <= 4) then
         flc_node = 1._r8/(1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))
      else
         select case (iswc)
         case (van_genuchten)
            write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
            call endrun(msg=errMsg(sourcefile, __LINE__))
            !          call unsatkVG_flc_from_psi(psi_node, &
            !             site_hydr%alpha_VG(1), &
            !	     site_hydr%n_VG(1),     &
            !             site_hydr%m_VG(1),     &
            !             site_hydr%l_VG(1),     &
            !             flc_node)
         case (campbell)
            call unsatkCampbell_flc_from_psi(psi_node, &
                 -1._r8*suc_sat*denh2o*grav_earth*m_per_mm*mpa_per_pa, &
                 bsw,  &
                 flc_node)
         case default
            write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select
      end if

    end associate

  end function flc_from_psi

  !===============================================================================!

  function dflcdpsi_from_psi(ft, pm, psi_node,  suc_sat, bsw) result(dflcdpsi_node)

    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity

    integer          , intent(in)     :: ft             ! PFT index
    integer          , intent(in)     :: pm             ! porous media index
    real(r8)         , intent(in)     :: psi_node       ! water potential   [MPa]
    real(r8), optional,intent(in)     :: suc_sat     ! minimum soil suction [mm]
    real(r8), optional,intent(in)     :: bsw         ! col Clapp and Hornberger "b"
    real(r8) :: dflcdpsi_node  ! deriv fractional loss of conductivity  [-] 

    associate(& 
         avuln    => pft_p%hydr_avuln_node, & ! vulnerability curve shape parameter          [-]
         p50      => pft_p%hydr_p50_node    & ! water potential at 50% loss of conductivity  [Pa]
         )

      if(pm <= 4) then
         dflcdpsi_node = -1._r8 * (1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))**(-2._r8) * &
              avuln(ft,pm)/p50(ft,pm)*(psi_node/p50(ft,pm))**(avuln(ft,pm)-1._r8)
      else
         select case (iswc)
         case (van_genuchten)    
            write(fates_log(),*) 'Van Genuchten plant hydraulics is '
            write(fates_log(),*) 'inoperable until further notice'
            call endrun(msg=errMsg(sourcefile, __LINE__))
            !call unsatkVG_dflcdpsi_from_psi(psi_node, &
            !      site_hydr%alpha_VG(1), &
            !      site_hydr%n_VG(1),     &
            !      site_hydr%m_VG(1),     &
            !      site_hydr%l_VG(1),     &
            !   dflcdpsi_node)
         case (campbell)
            call unsatkCampbell_dflcdpsi_from_psi(psi_node, &
                 -1._r8*suc_sat*denh2o*grav_earth*m_per_mm*mpa_per_pa, &
                 bsw,     &
                 dflcdpsi_node)
         case default
            write(fates_log(),*) 'ERROR: invalid soil water characteristic '
            write(fates_log(),*) 'function specified, iswc = '//char(iswc)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select
      end if

    end associate

  end function dflcdpsi_from_psi

  !===============================================================================!
  
  function th_from_psi(ft, pm, psi_node, th_sat, suc_sat, bsw) result(th_node)
    
    ! 
    ! Generic function that calls the correct specific functions for converting
    ! plant tissue or soil water potentials to volumetric water contents

    ! !ARGUMENTS
    integer          , intent(in)            :: ft          ! PFT index
    integer          , intent(in)            :: pm          ! porous media index
    real(r8)         , intent(in)            :: psi_node    ! water potential   [MPa]
    real(r8), optional,intent(in)     :: th_sat      ! water content at saturation
    ! (porosity for soil) [m3 m-3]
    real(r8), optional,intent(in)     :: suc_sat     ! minimum soil suction [mm]
    real(r8), optional,intent(in)     :: bsw         ! col Clapp and Hornberger "b"

    real(r8)          :: th_node     ! water content     [m3 m-3]

    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                ! lower bound of initial estimate         [m3 m-3]
    real(r8) :: upper                ! upper bound of initial estimate         [m3 m-3]
    real(r8) :: xtol                 ! error tolerance for x-variable          [m3 m-3]
    real(r8) :: ytol                 ! error tolerance for y-variable          [MPa]
    real(r8) :: satfrac              ! soil saturation fraction                [0-1]
    real(r8) :: psi_check


    associate(& 
         thetas   => pft_p%hydr_thetas_node  , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content
         resid    => pft_p%hydr_resid_node     & ! Input: [real(r8) (:,:) ] P-V curve: residual water fraction
         )

      if(pm <= 4) then

         lower  = thetas(ft,pm)*(resid(ft,pm) + 0.0001_r8)/cap_corr(pm)
         upper  = thetas(ft,pm)
         xtol   = 1.e-16_r8
         ytol   = 1.e-8_r8
         call bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
         psi_check = psi_from_th(ft, pm, th_node)

         if(psi_check > -1.e-8_r8) then
            write(fates_log(),*)'bisect_pv returned positive value for water potential at pm = ', char(pm)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

      else
         select case (iswc)
         case (van_genuchten)
            write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
            call endrun(msg=errMsg(sourcefile, __LINE__))
            !          call swcVG_satfrac_from_psi(psi_node, &
            !                  site_hydr%alpha_VG(1), &
            !                  site_hydr%n_VG(1),     &
            !                  site_hydr%m_VG(1),     &
            !                  site_hydr%l_VG(1),     &
            !                  satfrac)
            !          call swcVG_th_from_satfrac(satfrac, &
            !                  bc_in%watsat_sisl(1),   &
            !                  bc_in%watres_sisl(1),   &
            !                  th_node)
         case (campbell) 

            call swcCampbell_satfrac_from_psi(psi_node, &
                 (-1._r8)*suc_sat*denh2o*grav_earth*1.e-9_r8, &
                 bsw,     &
                 satfrac)
            call swcCampbell_th_from_satfrac(satfrac, &
                 th_sat, & 
                 th_node)
         case default
            write(fates_log(),*)  'invalid soil water characteristic function specified, iswc = '//char(iswc)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select
      end if

    end associate
    
  end function th_from_psi
  
  !===============================================================================!
  
  subroutine bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
    ! 
    ! !DESCRIPTION: Bisection routine for getting the inverse of the plant PV curve.
    !  An analytical solution is not possible because quadratic smoothing functions
    !  are used to remove discontinuities in the PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(inout)  :: lower       ! lower bound of estimate           [m3 m-3]
    real(r8)      , intent(inout)  :: upper       ! upper bound of estimate           [m3 m-3]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x-variable    [m3 m-3]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for y-variable    [MPa]
    real(r8)      , intent(in)     :: psi_node    ! water potential                   [MPa]
    real(r8)      , intent(out)    :: th_node     ! water content                     [m3 m-3]
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

    if(psi_node > 0.0_r8) then
       write(fates_log(),*)'Error: psi_note become positive, psi_node=',psi_node
       call endrun(msg=errMsg(sourcefile, __LINE__))  
    endif
    call psi_from_th(ft, pm, lower, y_lo)
    call psi_from_th(ft, pm, upper, y_hi)
    f_lo  = y_lo - psi_node
    f_hi  = y_hi - psi_node
    chg   = upper - lower
    nitr = 0
    do while(abs(chg) .gt. xtol .and. nitr < 100)
       x_new = 0.5_r8*(lower + upper)
       call psi_from_th(ft, pm, x_new, y_new)
       f_new = y_new - psi_node
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

    th_node = x_new

  end subroutine bisect_pv

  !===============================================================================!

  function psi_from_th(ft, pm, th_node, th_sat, suc_sat, bsw) result(psi_node)

    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content     [m3 m-3]
    real(r8), optional,intent(in)     :: th_sat      ! water content at saturation
    ! (porosity for soil) [m3 m-3]
    real(r8), optional,intent(in)     :: suc_sat     ! minimum soil suction [mm]
    real(r8), optional,intent(in)     :: bsw         ! col Clapp and Hornberger "b" 

    !
    ! !LOCAL VARIABLES:
    real(r8) :: satfrac                  ! saturation fraction [0-1]

    ! Result
    real(r8) :: psi_node    ! water potential   [MPa]


    if(pm <= 4) then       ! plant

       call tq2(ft, pm, th_node*cap_corr(pm), psi_node)

    else if(pm == 5) then  ! soil

       !! NOTE. FIX: The below sidesteps the problem of averaging potentially variable soil hydraulic properties with depth
       !!        and simply assigns the bulk soil (bucket) approximation of hydraulic properties as equal to the top soil layer.
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
          !          call swcVG_psi_from_th(th_node, &
          !                  bc_in%watsat_sisl(1),   &
          !                  bc_in%watres_sisl(1),   &
          !                  site_hydr%alpha_VG(1), &
          !                  site_hydr%n_VG(1),     &
          !                  site_hydr%m_VG(1),     &
          !                  site_hydr%l_VG(1),     &
          !                  psi_node)
       case (campbell)
          call swcCampbell_psi_from_th(th_node,th_sat,               &
               -1._r8*suc_sat*denh2o*grav_earth*m_per_mm*mpa_per_pa, &
               bsw,                                                  &
               psi_node)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

    end if

  end function psi_from_th

  !===============================================================================!
  
  function dpsidth_from_th(ft, pm, th_node, th_sat, suc_sat, bsw) result(dpsidth)
    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content                            [m3 m-3]
    real(r8), optional,intent(in)     :: th_sat      ! water content at saturation
    ! (porosity for soil) [m3 m-3]
    real(r8), optional,intent(in)     :: suc_sat     ! minimum soil suction [mm]
    real(r8), optional,intent(in)     :: bsw         ! col Clapp and Hornberger "b" 
    real(r8)         , intent(out)    :: dpsidth     ! derivative of water potential wrt theta  [MPa m3 m-3]

    !
    ! !LOCAL VARIABLES:

    real(r8) :: satfrac                  ! saturation fraction [0-1]


    if(pm <= 4) then       ! plant
       call dtq2dth(ft, pm, th_node*cap_corr(pm), dpsidth)
    else if(pm == 5) then  ! soil
       select case (iswc)
       case (van_genuchten)
          write(fates_log(),*) 'Van Genuchten plant hydraulics is inoperable until further notice'
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
          !call swcVG_dpsidth_from_th(th_node, &
          !        bc_in%watsat_sisl(1),   &
          !        bc_in%watres_sisl(1),   &
          !        site_hydr%alpha_VG(1), &
          !        site_hydr%n_VG(1),     &
          !        site_hydr%m_VG(1),     &
          !        site_hydr%l_VG(1),     &
          !        y)
       case (campbell)
          call swcCampbell_dpsidth_from_th(th_node, &
               th_sat, & 
               -1._r8*suc_sat*denh2o*grav_earth*m_per_mm*mpa_per_pa, &
               bsw, & 
               dpsidth)
       case default
          write(fates_log(),*) 'ERROR: invalid soil water characteristic function specified, iswc = '//char(iswc)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if

  end function dpsidth_from_th

  !===============================================================================!

  subroutine tq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
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
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------

    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    y = (-y_bq2 + sqrt(y_bq2*y_bq2 - 4._r8*beta2*y_cq2))/(2*beta2)

  end subroutine tq2

  !===============================================================================!

  subroutine dtq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
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
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: dydth_bq2              ! returned derivative from dbq2dth()
    real(r8) :: dydth_cq2              ! returned derivative from dcq2dth()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------

    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    call dbq2dth(ft, pm, x, dydth_bq2)
    call dcq2dth(ft, pm, x, dydth_cq2)
    y = 1._r8/(2._r8*beta2)*(-dydth_bq2 + 0.5_r8*((y_bq2*y_bq2 - 4._r8*beta2*y_cq2)**(-0.5_r8)) * &
         (2._r8*y_bq2*dydth_bq2 - 4._r8*beta2*dydth_cq2))

  end subroutine dtq2dth

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       y = (-y_bq1 - sqrt(y_bq1*y_bq1 - 4._r8*beta1*y_cq1))/(2*beta1)
    end if !porous media

  end subroutine tq1

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

  !===============================================================================!

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

    call solutepsi(ft, pm, x, y_solute)
    y = y_solute

  end subroutine cavitationPV

  !===============================================================================!

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

  !===============================================================================!

  subroutine elasticPV(ft, pm, x, y)
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
    real(r8) :: y_pressure         ! returned y (psi) value from pressurepsi()
    !----------------------------------------------------------------------

    call solutepsi(ft, pm, x, y_solute)
    call pressurepsi(ft, pm, x, y_pressure)
    y = y_solute + y_pressure

  end subroutine elasticPV

  !===============================================================================!

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

  !===============================================================================!

  subroutine solutepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes solute water potential (negative) as a function of
    !  water content for the plant PV curve.
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
         pinot   => pft_p%hydr_pinot_node,  & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pft_p%hydr_thetas_node, & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         resid   => pft_p%hydr_resid_node   & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
      
      y = pinot(ft,pm)*thetas(ft,pm)*(rwcft(pm) - resid(ft,pm)) / &
           (x - thetas(ft,pm)*resid(ft,pm))
      
    end associate
    
  end subroutine solutepsi
  
  !===============================================================================!

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
         pinot   => pft_p%hydr_pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pft_p%hydr_thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         resid   => pft_p%hydr_resid_node       & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )

      y = -1._r8*thetas(ft,pm)*pinot(ft,pm)*(rwcft(pm) - resid(ft,pm)) / &
           ((x - thetas(ft,pm)*resid(ft,pm))**2._r8)
      
    end associate
  end subroutine dsolutepsidth

  !===============================================================================!

  subroutine pressurepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes pressure water potential (positive) as a function of
    !  water content for the plant PV curve.
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
         pinot   => pft_p%hydr_pinot_node     , & ! P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pft_p%hydr_thetas_node    , & ! P-V curve: saturated volumetric water content for node   [m3 m-3]
         resid   => pft_p%hydr_resid_node     , & ! P-V curve: residual fraction                             [-]
         epsil   => pft_p%hydr_epsil_node       & ! P-V curve: bulk elastic modulus                          [MPa]
         )
      
      y = epsil(ft,pm) * (x - thetas(ft,pm)*rwcft(pm)) / &
           (thetas(ft,pm)*(rwcft(pm)-resid(ft,pm))) - pinot(ft,pm)
      
    end associate
  end subroutine pressurepsi
  
  !===============================================================================!

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
         thetas  => pft_p%hydr_thetas_node, & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         resid   => pft_p%hydr_resid_node , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => pft_p%hydr_epsil_node   & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
      
      y = epsil(ft,pm)/(thetas(ft,pm)*(rwcft(pm) - resid(ft,pm)))
      
    end associate
    
  end subroutine dpressurepsidth
  
  !===============================================================================!

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
         thetas    => pft_p%hydr_thetas_node     & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
      
      y = cap_int(pm) + cap_slp(pm)/thetas(ft,pm)*x
      
    end associate
  end subroutine capillaryPV

  !===============================================================================!

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
         thetas    => pft_p%hydr_thetas_node    & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
      
      y = cap_slp(pm)/thetas(ft,pm)
      
    end associate
    
  end subroutine dcapillaryPVdth

  !===============================================================================!

  subroutine swcVG_satfrac_from_th(th, watsat, watres, satfrac)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns saturation fraction given water content, porosity, and residual water content
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: th       !soil volumetric water content       [m3 m-3]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
    real(r8), intent(out)           :: satfrac  !saturation fraction                 [0-1]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    satfrac = (th - watres)/(watsat - watres)

  end subroutine swcVG_satfrac_from_th

  !===============================================================================!

  subroutine swcCampbell_satfrac_from_th(th, watsat, satfrac)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns saturation fraction given water content and porosity
    !
    !USES
    !ARGUMENTS:
    real(r8), intent(in)            :: th       !soil volumetric water content       [m3 m-3]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(out)           :: satfrac  !saturation fraction                 [0-1]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    satfrac = th/watsat

  end subroutine swcCampbell_satfrac_from_th

  !===============================================================================!

  subroutine swcVG_psi_from_th(th, watsat, watres, alpha, n, m, l, psi)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns water potential given water content and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: th       !volumetric water content       [m3 m-3]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
    real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
    real(r8), intent(in)            :: n        !pore-size distribution index   [-]
    real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
    real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
    real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8)                        :: satfrac  !saturation fraction            [0-1]
    !------------------------------------------------------------------------------

    call swcVG_satfrac_from_th(th, watsat, watres, satfrac)
    call swcVG_psi_from_satfrac(satfrac, alpha, n, m, l, psi)

  end subroutine swcVG_psi_from_th

  !===============================================================================!

  subroutine swcCampbell_psi_from_th(th, watsat, psisat, B, psi)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns water potential given saturation fraction, air-entry pressure and shape parameter
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: th       !volumetric water content       [m3 m-3]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)            :: psisat   !air entry pressure             [MPa]
    real(r8), intent(in)            :: B        !shape parameter                [-]
    real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8)                        :: satfrac  !saturation fraction            [0-1]
    !------------------------------------------------------------------------------

    call swcCampbell_satfrac_from_th(th, watsat, satfrac)
    call swcCampbell_psi_from_satfrac(satfrac, psisat, B, psi)

  end subroutine swcCampbell_psi_from_th

  !===============================================================================!

  subroutine swcVG_psi_from_satfrac(satfrac, alpha, n, m, l, psi)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns water potential given saturation fraction and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: satfrac  !saturation fraction            [0-1]
    real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
    real(r8), intent(in)            :: n        !pore-size distribution index   [-]
    real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
    real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
    real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    psi = -1._r8/alpha*(satfrac**(-1._r8/m)-1._r8)**(1._r8/n)

  end subroutine swcVG_psi_from_satfrac

  !===============================================================================!

  subroutine swcCampbell_psi_from_satfrac(satfrac, psisat, B, psi)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns water potential given saturation fraction, air-entry pressure and shape parameter
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: satfrac  !saturation fraction            [0-1]
    real(r8), intent(in)            :: psisat   !air entry pressure             [MPa]
    real(r8), intent(in)            :: B        !shape parameter                [-]
    real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    psi = psisat*(satfrac**(-B))

  end subroutine swcCampbell_psi_from_satfrac

  !===============================================================================!

  subroutine swcVG_th_from_satfrac(satfrac, watsat, watres, th)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns water content given saturation fraction, porosity and residual water content
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: satfrac  !saturation fraction                 [0-1]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)            :: watres   !volumetric residual soil water      [m3 m-3]
    real(r8), intent(out)           :: th       !soil volumetric water content       [m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    th = watres + satfrac*(watsat - watres)

  end subroutine swcVG_th_from_satfrac

  !===============================================================================!

  subroutine swcCampbell_th_from_satfrac(satfrac, watsat, th)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns water content given saturation fraction and porosity
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: satfrac  !saturation fraction                 [0-1]
    real(r8), intent(in)            :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(out)           :: th       !soil volumetric water content       [m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    th = satfrac*watsat

  end subroutine swcCampbell_th_from_satfrac

  !======================================================================-
  subroutine swcVG_satfrac_from_psi(psi, alpha, n, m, l, satfrac)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns saturation fraction given water potential and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: psi      !soil matric potential          [MPa]
    real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
    real(r8), intent(in)            :: n        !pore-size distribution index   [-]
    real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
    real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
    real(r8), intent(out)           :: satfrac  !soil saturation fraction       [0-1]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    satfrac = (1._r8/(1._r8 + (alpha*abs(psi))**n))**m

  end subroutine swcVG_satfrac_from_psi

  !======================================================================-
  subroutine swcCampbell_satfrac_from_psi(psi, psisat, B, satfrac)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns saturation fraction given water potential and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)            :: psi      !soil matric potential          [MPa]
    real(r8), intent(in)            :: psisat   !air-entry pressure             [MPa]
    real(r8), intent(in)            :: B        !shape parameter                [-]
    real(r8), intent(out)           :: satfrac  !soil saturation fraction       [0-1]
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------------

    satfrac = (psi/psisat)**(-1.0_r8/B)

  end subroutine swcCampbell_satfrac_from_psi

  !======================================================================-
  subroutine swcVG_dpsidth_from_th(th, watsat, watres, alpha, n, m, l, dpsidth)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given water content and SWC parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: th       !volumetric water content                       [m3 m-3]
    real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
    real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
    real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
    real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
    real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
    real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: satfrac  !saturation fraction            [0-1]
    !------------------------------------------------------------------------------

    call swcVG_satfrac_from_th(th, watsat, watres, satfrac)
    call swcVG_dpsidth_from_satfrac(satfrac, watsat, watres, alpha, n, m, l, dpsidth)

  end subroutine swcVG_dpsidth_from_th

  !======================================================================-
  subroutine swcCampbell_dpsidth_from_th(th, watsat, psisat, B, dpsidth)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given water content and SWC parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: th       !volumetric water content                       [m3 m-3]
    real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)  :: psisat   !air entry pressure                             [MPa]
    real(r8), intent(in)  :: B        !shape parameter                                [-]
    real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: satfrac  !saturation fraction            [0-1]
    !------------------------------------------------------------------------------

    call swcCampbell_satfrac_from_th(th, watsat, satfrac)
    call swcCampbell_dpsidth_from_satfrac(satfrac, watsat, psisat, B, dpsidth)

  end subroutine swcCampbell_dpsidth_from_th

  !======================================================================-
  subroutine swcVG_dpsidth_from_satfrac(satfrac, watsat, watres, alpha, n, m, l, dpsidth)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given saturation fraction and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: satfrac  !saturation fraction                            [0-1]
    real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
    real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
    real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
    real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
    real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
    real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: temp0    !temporary
    real(r8)              :: temp1    !temporary
    real(r8)              :: temp2    !temporary
    real(r8)              :: temp3    !temporary
    !------------------------------------------------------------------------------

    temp0   = 1._r8/(m*n*alpha*(watsat-watres))
    temp1   = satfrac**(-1._r8/m) - 1._r8
    temp2   = temp1**(1._r8/n - 1._r8)
    temp3   = satfrac**(-1._r8/m - 1._r8)
    dpsidth = temp0*temp2*temp3

  end subroutine swcVG_dpsidth_from_satfrac

  !======================================================================-
  subroutine swcCampbell_dpsidth_from_satfrac(satfrac, watsat, psisat, B, dpsidth)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given saturation fraction and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: satfrac  !saturation fraction                            [0-1]
    real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
    real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
    real(r8), intent(in)  :: B        !shape parameter                                [-]
    real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
    !------------------------------------------------------------------------------

    dpsidth = psisat*(-B)/watsat*(satfrac)**(-B-1._r8)

  end subroutine swcCampbell_dpsidth_from_satfrac

  !======================================================================-
  subroutine unsatkVG_flc_from_psi(psi, alpha, n, m, l, flc)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns unsaturated hydraulic conductivity 
    ! given water potential and SWC parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
    real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
    real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
    real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
    real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
    real(r8), intent(out) :: flc      !k/ksat ('fractional loss of conductivity')     [-]
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: temp          !temporary
    real(r8)              :: fac1a         !temporary
    real(r8)              :: fac1b         !temporary
    real(r8)              :: fac1          !temporary
    real(r8)              :: fac2          !temporary
    !------------------------------------------------------------------------------

    temp       = ( alpha*abs(psi)      ) ** (n)
    fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
    fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
    fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
    fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)

    flc        =   fac1 * fac2

  end subroutine unsatkVG_flc_from_psi

  !======================================================================-
  
  subroutine unsatkCampbell_flc_from_psi(psi, psisat, B, flc)
    !
    ! DESCRIPTION
    ! Campbell (1974) soil water characteristic (retention) curve
    ! returns unsaturated hydraulic conductivity 
    ! given water potential and SWC parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
    real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
    real(r8), intent(in)  :: B        !shape parameter                                [-]
    real(r8), intent(out) :: flc      !k/ksat ('fractional loss of conductivity')     [-]
    !------------------------------------------------------------------------------

    flc        =  (psi/psisat)**(-2._r8-3._r8/B)

  end subroutine unsatkCampbell_flc_from_psi

  !======================================================================-
  subroutine unsatkVG_dflcdpsi_from_psi(psi, alpha, n, m, l, dflcdpsi)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given saturation fraction and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
    real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
    real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
    real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
    real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
    real(r8), intent(out) :: dflcdpsi !derivative of k/ksat (flc) wrt psi             [MPa-1]
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: temp          !temporary
    real(r8)              :: fac1a         !temporary
    real(r8)              :: fac1b         !temporary
    real(r8)              :: fac1          !temporary
    real(r8)              :: fac2          !temporary
    real(r8)              :: dtemp         !temporary
    real(r8)              :: dfac1adpsi    !temporary
    real(r8)              :: dfac1bdpsi    !temporary
    real(r8)              :: dfac1dpsi     !temporary
    real(r8)              :: dfac2dpsi     !temporary
    !------------------------------------------------------------------------------

    temp       = ( alpha*abs(psi)      ) ** (n)
    fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
    fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
    fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
    fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)

    dtemp      =   n * alpha * ( alpha*abs(psi) ) ** (n-1._r8)
    dfac1adpsi = ( n-1._r8 ) * alpha * ( alpha*abs(psi) ) ** (n-2._r8)
    dfac1bdpsi = ( -1._r8 ) * m * dtemp * ( 1._r8 + temp ) ** (-1._r8*m - 1._r8)
    dfac1dpsi  = ( 2._r8 ) * ( 1._r8 - fac1a*fac1b ) * ( -1._r8*dfac1bdpsi*fac1a - dfac1adpsi*fac1b )
    dfac2dpsi  = ( -0.5_r8 ) * m * dtemp * (1._r8 + temp)**(-0.5_r8*m-1._r8)

    dflcdpsi   = ( -1._r8 ) * ( dfac2dpsi*fac1 + dfac1dpsi*fac2 )    ! BOC... mult by -1 because unsatk eqn is based on abs(psi)

  end subroutine unsatkVG_dflcdpsi_from_psi

  !======================================================================-
  subroutine unsatkCampbell_dflcdpsi_from_psi(psi, psisat, B, dflcdpsi)
    !
    ! DESCRIPTION
    ! van Genuchten (1980) soil water characteristic (retention) curve
    ! returns derivative of water water potential with respect to water content
    ! given saturation fraction and shape parameters
    !
    !USES
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
    real(r8), intent(in)  :: psisat   !air-entry pressure                             [MPa]
    real(r8), intent(in)  :: B        !shape parameter                                [-]
    real(r8), intent(out) :: dflcdpsi !derivative of k/ksat (flc) wrt psi             [MPa-1]
    !------------------------------------------------------------------------------

    dflcdpsi   = psisat*(-2._r8-3._r8/B)*(psi/psisat)**(-3._r8-3._r8/B)

  end subroutine unsatkCampbell_dflcdpsi_from_psi


  ! =====================================================================================
  ! Utility Functions
  ! =====================================================================================

  subroutine bisect_rootfr(a, b, lower_init, upper_init, xtol, ytol, crootfr, x_new)
    ! 
    ! !DESCRIPTION: Bisection routine for getting the inverse of the cumulative root
    !  distribution. No analytical soln bc crootfr ~ exp(ax) + exp(bx).
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)     :: a, b        ! pft root distribution constants
    real(r8)      , intent(in)     :: lower_init  ! lower bound of initial x estimate [m]
    real(r8)      , intent(in)     :: upper_init  ! upper bound of initial x estimate [m]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x_new         [m]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for crootfr       [-]
    real(r8)      , intent(in)     :: crootfr     ! cumulative root fraction at x_new [-]
    real(r8)      , intent(out)    :: x_new       ! soil depth                        [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                  ! lower bound x estimate [m]
    real(r8) :: upper                  ! upper bound x estimate [m]
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    !----------------------------------------------------------------------

    lower = lower_init
    upper = upper_init
    f_lo  = zeng2001_crootfr(a, b, lower) - crootfr
    f_hi  = zeng2001_crootfr(a, b, upper) - crootfr
    chg   = upper - lower
    do while(abs(chg) .gt. xtol)
       x_new = 0.5_r8*(lower + upper)
       f_new = zeng2001_crootfr(a, b, x_new) - crootfr
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
    end do
  end subroutine bisect_rootfr

  ! =====================================================================================

  function zeng2001_crootfr(a, b, z, z_max) result(crootfr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: a,b    ! pft parameters
    real(r8) , intent(in) :: z      ! soil depth (m)
    real(r8) , intent(in), optional :: z_max ! max soil depth (m)
    !
    real(r8) :: crootfr_max

    ! !RESULT
    real(r8) :: crootfr            ! cumulative root fraction
    !
    !------------------------------------------------------------------------
    crootfr      = 1._r8 - .5_r8*(exp(-a*z) + exp(-b*z))


    ! If a maximum rooting depth is provided, then
    ! we force everything to sum to unity. We do this by 
    ! simply dividing through by the maximum possible
    ! root fraction.

    if(present(z_max))then
       crootfr_max = 1._r8 - .5_r8*(exp(-a*z_max) + exp(-b*z_max))
       crootfr = crootfr/crootfr_max
    end if

    if(debug)then
       if(present(z_max))then
          if((crootfr_max<nearzero) .or. (crootfr_max>1.0_r8) )then
             write(fates_log(),*) 'problem scaling crootfr in zeng2001'
             write(fates_log(),*) 'z_max: ',z_max
             write(fates_log(),*) 'crootfr_max: ',crootfr_max
          end if
       end if
    end if


    return

  end function zeng2001_crootfr

  ! =====================================================================================

  subroutine shellGeom(l_aroot, rs1, area_site, dz, r_out_shell, r_node_shell, v_shell)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil surrounding fine roots remains
    ! the same.  
    !
    ! !USES:

    !
    ! !ARGUMENTS:
    real(r8)     , intent(in)             :: l_aroot              ! Total length of absorbing roots
    ! for the whole site, this layer (m)
    real(r8)     , intent(in)             :: rs1                  ! Fine root radius (m)
    real(r8)     , intent(in)             :: area_site            ! Area of site (10,000 m2)
    real(r8)     , intent(in)             :: dz                   ! Width of current soil layer (m)
    real(r8)     , intent(out)            :: r_out_shell(:)       ! Outer radius of each shell (m)
    real(r8)     , intent(out)            :: r_node_shell(:)      ! Radius of the shell's midpoint
    real(r8)     , intent(out)            :: v_shell(:)           ! volume of the rhizosphere shells (m3/ha)
    ! for this layer
    !
    ! !LOCAL VARIABLES:
    integer                        :: k                                 ! rhizosphere shell indicies
    !-----------------------------------------------------------------------

    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    r_out_shell(nshell) = (pi_const*l_aroot/(area_site*dz))**(-0.5_r8)                  ! eqn(8) S98
    if(nshell > 1) then
       do k = 1,nshell-1
          r_out_shell(k)   = rs1*(r_out_shell(nshell)/rs1)**((real(k,r8))/real(nshell,r8))  ! eqn(7) S98
       enddo
    end if

    ! set nodal (midpoint) radii of these shells
    ! BOC...not doing this as it requires PFT-specific fine root thickness, but this is at column level
    r_node_shell(1) = 0.5_r8*(rs1 + r_out_shell(1))
    !r_node_shell(1) = 0.5_r8*(r_out_shell(1))

    do k = 2,nshell
       r_node_shell(k) = 0.5_r8*(r_out_shell(k-1) + r_out_shell(k))
    enddo

    ! update volumes
    if(voltype==bcvol)then
       do k = 1,nshell
          if(k == 1) then
             ! BOC...not doing this as it requires PFT-specific fine root thickness but this is at column level
             v_shell(k)   = pi_const*dz*(r_out_shell(k)**2._r8 - rs1**2._r8)

          else
             v_shell(k)   = pi_const*dz*(r_out_shell(k)**2._r8 - r_out_shell(k-1)**2._r8)
          end if
       enddo
    elseif(voltype==rkvol)then
       do k = 1,nshell
          if(k == 1) then
             v_shell(k) = pi_const*l_aroot*(r_out_shell(k)**2._r8 - rs1**2._r8)
          else
             v_shell(k) = pi_const*l_aroot*(r_out_shell(k)**2._r8 - r_out_shell(k-1)**2._r8)
          end if
       enddo
    end if

  end subroutine shellGeom

  ! =====================================================================================

  function xylemtaper(p, dz) result(chi_tapnotap)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: p      ! Savage et al. (2010) taper exponent                                                                [-]
    real(r8) , intent(in) :: dz     ! hydraulic distance from petiole to node of interest                                                [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: atap,btap           ! scaling exponents for total conductance ~ tree size (ratio of stem radius to terminal twig radius)
    real(r8) :: anotap,bnotap       ! same as atap, btap, but not acounting for xylem taper (Savage et al. (2010) p = 0)
    ! NOTE: these scaling exponents were digitized from Fig 2a of Savage et al. (2010)
    ! Savage VM, Bentley LP, Enquist BJ, Sperry JS, Smith DD, Reich PB, von Allmen EI. 2010.
    !    Hydraulic trade-offs and space filling enable better predictions of vascular structure
    !    and function in plants. Proceedings of the National Academy of Sciences 107(52): 22722-22727.
    real(r8) :: lN=0.04_r8          ! petiole length                                                                                     [m]
    real(r8) :: little_n=2._r8      ! number of daughter branches per parent branch, assumed constant throughout tree (self-similarity)  [-]
    real(r8) :: big_n               ! number of branching levels (allowed here to take on non-integer values): increases with tree size  [-]
    real(r8) :: ktap                ! hydraulic conductance along the pathway, accounting for xylem taper                                [kg s-1 MPa-1]
    real(r8) :: knotap              ! hydraulic conductance along the pathway, not accounting for xylem taper                            [kg s-1 MPa-1]
    real(r8) :: num                 ! temporary
    real(r8) :: den                 ! temporary
    !
    ! !RESULT
    real(r8) :: chi_tapnotap        ! ratio of total tree conductance accounting for xylem taper to that without, over interval dz
    !
    !------------------------------------------------------------------------

    anotap  = 7.19903e-13_r8
    bnotap  = 1.326105578_r8
    if (p >= 1.0_r8) then
       btap  = 2.00586217_r8
       atap  = 1.82513E-12_r8
    else if (p >= (1._r8/3._r8) .AND. p < 1._r8) then
       btap  = 1.854812819_r8
       atap  = 6.66908E-13_r8
    else if (p >= (1._r8/6._r8) .AND. p < (1._r8/3._r8)) then
       btap  = 1.628179741_r8
       atap  = 6.58345E-13_r8
    else
       btap  = bnotap
       atap  = anotap
    end if

    num          = 3._r8*log(1._r8 - dz/lN * (1._r8-little_n**(1._r8/3._r8)))
    den          = log(little_n)
    big_n        = num/den - 1._r8
    ktap         = atap   * (little_n**(big_N*  btap/2._r8))
    knotap       = anotap * (little_n**(big_N*bnotap/2._r8))
    chi_tapnotap = ktap / knotap

    return

  end function xylemtaper


end module FatesHydroUnitFunctionsMod
