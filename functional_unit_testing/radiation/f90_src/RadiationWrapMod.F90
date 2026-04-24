module RadiationWrapMod

  use TwoStreamMLPEMod
  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double
  
  implicit none
  public
  save
  
  integer(kind=c_int), parameter :: param_string_length = 32
  
  type(twostream_type) :: twostream
  
  
contains
  
  subroutine InitAllocate(n_layer,n_column)

    integer(kind=c_int), intent(in) :: n_layer
    integer(kind=c_int), intent(in) :: n_column
    
    integer(kind=c_int) :: ican


    call twostream%AllocInitTwoStream((/1,2/),n_layer,n_column)

    
    twostream%n_lyr = n_layer

    do ican = 1,n_layer
       twostream%n_col(ican) = n_column
    end do

    twostream%force_prep = .true.

    call twostream%GetNScel()

    twostream%frac_snow = 0._r8
    twostream%frac_snow_old = 1._r8
    
    print*,"Allocated twostream instance"
    print*," with ",twostream%n_scel," elements"
    
    return
  end subroutine InitAllocate


  subroutine Dealloc()

    call twostream%DeallocTwoStream()
    
  end subroutine Dealloc

  
  subroutine SetRadParam(val,pft,ib,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft
    integer(kind=c_int), intent(in) :: ib

    select case(trim(pname))
    case('rhol')
       rad_params%rhol(ib,pft) = val
    case('rhos')
       rad_params%rhos(ib,pft) = val
    case('taul')
       rad_params%taul(ib,pft) = val
    case('taus')
       rad_params%taus(ib,pft) = val
    case('xl')
       rad_params%xl(pft) = val     
    case('clumping_index')
       rad_params%clumping_index(pft) = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
       
  end subroutine SetRadParam

  ! =============================================================================
  
  subroutine SetGroundSnow(ib,val,pname)

    real(r8), intent(in)            :: val
    integer,  intent(in)            :: ib
    character(kind=c_char,len=*), intent(in)    :: pname

    select case(trim(pname))
    case('albedo_grnd_diff')
       twostream%band(ib)%albedo_grnd_diff = val
    case('albedo_grnd_beam')
       twostream%band(ib)%albedo_grnd_beam = val
       case default
       print*,"An unknown parameter name was sent to ground/snow"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
  end subroutine SetGroundSnow

  ! =============================================================================
  
  subroutine SetupCanopy(ican,icol,pft,area,lai,sai)

    integer(kind=c_int), intent(in) :: ican  ! Canopy layer index
    integer(kind=c_int), intent(in) :: icol  ! Column (pft) position index
    integer(kind=c_int), intent(in) :: pft   ! PFT index
    real(r8), intent(in)            :: area  ! columns fraction of the ground
    real(r8), intent(in)            :: lai   ! LAI
    real(r8), intent(in)            :: sai
    
    
    twostream%scelg(ican,icol)%pft  = pft
    twostream%scelg(ican,icol)%area = area
    twostream%scelg(ican,icol)%lai  = lai
    twostream%scelg(ican,icol)%sai  = sai

    return
  end subroutine SetupCanopy
  
  subroutine WrapCanopyPrep(frac_snow)

    real(kind=r8),intent(in)       :: frac_snow
    
    call twostream%CanopyPrep(frac_snow)
    
  end subroutine WrapCanopyPrep
  
  subroutine WrapZenithPrep(cosz)

    real(kind=r8),intent(in) :: cosz
    
    call twostream%ZenithPrep(cosz)

    return
  end subroutine WrapZenithPrep

  subroutine WrapSetDownwelling(ib,Rbeam_atm,Rdiff_atm)
    
    integer(c_int) :: ib
    real(r8)  :: Rbeam_atm          ! Intensity of beam radiation at top of canopy [W/m2 ground]
    real(r8)  :: Rdiff_atm          ! Intensity of diffuse radiation at top of canopy [W/m2 ground]    

    twostream%band(ib)%Rbeam_atm = Rbeam_atm
    twostream%band(ib)%Rdiff_atm = Rdiff_atm
    
    return
  end subroutine WrapSetDownwelling
    
  
  subroutine WrapSolve(ib,boundary_type,Rbeam_atm,Rdiff_atm, &
       albedo_beam, & 
       albedo_diff, &
       frac_abs_can_beam, &
       frac_abs_can_diff, &
       frac_beam_grnd_beam, &
       frac_diff_grnd_beam, &
       frac_diff_grnd_diff)

    integer(c_int) :: ib
    integer(c_int) :: boundary_type

    real(r8)  :: albedo_beam
    real(r8)  :: albedo_diff
    real(r8)  :: err_solve
    real(r8)  :: err_consv
    real(r8)  :: frac_abs_can_beam
    real(r8)  :: frac_abs_can_diff
    real(r8)  :: frac_beam_grnd_beam
    real(r8)  :: frac_diff_grnd_beam
    real(r8)  :: frac_diff_grnd_diff
    real(r8)  :: Rbeam_atm          ! Intensity of beam radiation at top of canopy [W/m2 ground]
    real(r8)  :: Rdiff_atm          ! Intensity of diffuse radiation at top of canopy [W/m2 ground]

    real(r8) :: taulamb(50)
    real(r8) :: omega(50,50)
    integer  :: ipiv(50)
    
    call twostream%Solve(ib,boundary_type, & 
         Rbeam_atm,Rdiff_atm, &
         taulamb, &
         omega, &
         ipiv, &
         albedo_beam, & 
         albedo_diff, &
         err_consv,   &
         frac_abs_can_beam, &
         frac_abs_can_diff, &
         frac_beam_grnd_beam, &
         frac_diff_grnd_beam, &
         frac_diff_grnd_diff)

    return
  end subroutine WrapSolve

  subroutine WrapGetIntensity(ican,icol,ib,vai,r_diff_dn,r_diff_up,r_beam)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: vai
    real(r8)    :: r_diff_dn
    real(r8)    :: r_diff_up
    real(r8)    :: r_beam
    
    r_diff_dn = twostream%GetRdDn(ican,icol,ib,vai)
    r_diff_up = twostream%GetRdUp(ican,icol,ib,vai)
    r_beam    = twostream%GetRb(ican,icol,ib,vai)

    return
  end subroutine WrapGetIntensity

  subroutine WrapGetAbsRad(ican,icol,ib,vai_top,vai_bot,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: vai_top,vai_bot
    real(r8)    :: Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac,Rb_abs,Rd_abs

    call twostream%GetAbsRad(ican,icol,ib,vai_top,vai_bot,Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac)
    
    return
  end subroutine WrapGetAbsRad
  
  subroutine WrapGetParams(ican,icol,ib,Kb,Kd,om,betad,betab)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: Kb,Kd,om,betad,betab

    Kb    = twostream%scelg(ican,icol)%Kb
    Kd    = twostream%scelg(ican,icol)%Kd
    om    = twostream%band(ib)%scelb(ican,icol)%om
    betad = twostream%band(ib)%scelb(ican,icol)%betad
    betab = twostream%band(ib)%scelb(ican,icol)%betab
    
    return
  end subroutine WrapGetParams

  subroutine WrapForceParams(ican,icol,ib,val,pname)

    ! This will overwrite the 2-stream parameters
    ! that are derived from the fates params

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    
    select case(trim(pname))
    case('Kb')
       twostream%scelg(ican,icol)%Kb = val
    case('Kd')
       twostream%scelg(ican,icol)%Kd = val
    case('om')
       twostream%band(ib)%scelb(ican,icol)%om = val
    case('betab')
       twostream%band(ib)%scelb(ican,icol)%betab = val
    case('betad')
       twostream%band(ib)%scelb(ican,icol)%betad = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
  
  end subroutine WrapForceParams
  
end module RadiationWrapMod
