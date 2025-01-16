module LeafBiophysSuppMod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double

  use LeafBiophysicsMod, only : lb_params
  
  implicit none
  public
  save
  
  integer(kind=c_int), parameter :: param_string_length = 32
  logical, parameter :: debug = .true.
  
contains
  
  subroutine AllocLeafParam(numpft)

    integer, intent(in) :: numpft
    
    allocate(lb_params%c3psn(numpft))
    allocate(lb_params%stomatal_btran_model(numpft))
    allocate(lb_params%agross_btran_model(numpft))
    allocate(lb_params%medlyn_slope(numpft))
    allocate(lb_params%bb_slope(numpft))
    allocate(lb_params%stomatal_intercept(numpft))
    allocate(lb_params%maintresp_leaf_ryan1991_baserate(numpft))
    allocate(lb_params%maintresp_leaf_atkin2017_baserate(numpft))
    allocate(lb_params%maintresp_reduction_curvature(numpft))
    allocate(lb_params%maintresp_reduction_intercept(numpft))
    allocate(lb_params%maintresp_reduction_upthresh(numpft))
    allocate(lb_params%vcmaxha(numpft))
    allocate(lb_params%jmaxha(numpft))
    allocate(lb_params%vcmaxhd(numpft))
    allocate(lb_params%jmaxhd(numpft))
    allocate(lb_params%vcmaxse(numpft))
    allocate(lb_params%jmaxse(numpft))

    return
  end subroutine AllocLeafParam

  subroutine DeallocLeafParam()
    
    deallocate(lb_params%c3psn)
    deallocate(lb_params%stomatal_btran_model)
    deallocate(lb_params%agross_btran_model)
    deallocate(lb_params%medlyn_slope)
    deallocate(lb_params%bb_slope)
    deallocate(lb_params%stomatal_intercept)
    deallocate(lb_params%maintresp_leaf_ryan1991_baserate)
    deallocate(lb_params%maintresp_leaf_atkin2017_baserate)
    deallocate(lb_params%maintresp_reduction_curvature)
    deallocate(lb_params%maintresp_reduction_intercept)
    deallocate(lb_params%maintresp_reduction_upthresh)
    deallocate(lb_params%vcmaxha)
    deallocate(lb_params%jmaxha)
    deallocate(lb_params%vcmaxhd)
    deallocate(lb_params%jmaxhd)
    deallocate(lb_params%vcmaxse)
    deallocate(lb_params%jmaxse)
    
  end subroutine DeallocLeafParam
  
  ! =====================================================================================
  
  subroutine SetLeafParam(val,pft,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft

    select case(trim(pname))
    case('fates_daylength_factor_switch')
       lb_params%dayl_switch = nint(val)
    case('fates_leaf_stomatal_model')
       lb_params%stomatal_model = nint(val)
    case('fates_leaf_stomatal_assim_model')
       lb_params%stomatal_assim_model = nint(val)
    case('fates_leaf_photo_tempsens_model')
       lb_params%photo_tempsens_model = nint(val)
    case('fates_leaf_c3psn')
       lb_params%c3psn(pft) = nint(val)
    case('fates_leaf_stomatal_btran_model')
       lb_params%stomatal_btran_model(pft) = nint(val)
    case('fates_leaf_agross_btran_model')
       lb_params%agross_btran_model(pft) = nint(val)
    case('fates_leaf_stomatal_slope_ballberry')
       lb_params%bb_slope(pft) = val
    case('fates_leaf_stomatal_slope_medlyn')
       lb_params%medlyn_slope(pft) = val
    case('fates_leaf_stomatal_intercept')
       lb_params%stomatal_intercept(pft) = val
    case('fates_maintresp_reduction_curvature')
       lb_params%maintresp_reduction_curvature(pft) = val
    case('fates_maintresp_reduction_intercept')
       lb_params%maintresp_reduction_intercept(pft) = val
    case('fates_maintresp_reduction_upthresh')
       lb_params%maintresp_reduction_upthresh(pft) = val
    case('fates_maintresp_leaf_atkin2017_baserate')
       lb_params%maintresp_leaf_atkin2017_baserate(pft) = val
    case('fates_maintresp_leaf_ryan1991_baserate')
       lb_params%maintresp_leaf_ryan1991_baserate(pft) = val
    case('fates_leaf_vcmaxha')
       lb_params%vcmaxha(pft) = val
    case('fates_leaf_jmaxha')
       lb_params%jmaxha(pft) = val
    case('fates_leaf_vcmaxhd')
       lb_params%vcmaxhd(pft) = val
    case('fates_leaf_jmaxhd')
       lb_params%jmaxhd(pft) = val
    case('fates_leaf_vcmaxse')
       lb_params%vcmaxse(pft) = val
    case('fates_leaf_jmaxse')
       lb_params%jmaxse(pft) = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
       
  end subroutine SetLeafParam
  
  subroutine DumpParams()

    print*,"Leaf Biophys Parameters In Structure lb_params: "
    print*,'fates_daylength_factor_switch: ',lb_params%dayl_switch
    print*,'fates_leaf_stomatal_model: ',lb_params%stomatal_model
    print*,'fates_leaf_stomatal_assim_model: ',lb_params%stomatal_assim_model
    print*,'fates_leaf_photo_tempsens_model: ',lb_params%photo_tempsens_model
    print*,'fates_leaf_c3psn: ',lb_params%c3psn
    print*,'fates_leaf_stomatal_slope_ballberry: ',lb_params%bb_slope
    print*,'fates_leaf_stomatal_slope_medlyn: ',lb_params%medlyn_slope
    print*,'fates_leaf_stomatal_intercept: ',lb_params%stomatal_intercept
    print*,'fates_maintresp_reduction_curvature: ',lb_params%maintresp_reduction_curvature
    print*,'fates_maintresp_reduction_intercept: ',lb_params%maintresp_reduction_intercept
    print*,'fates_maintresp_reduction_upthresh: ',lb_params%maintresp_reduction_upthresh
    print*,'fates_maintresp_leaf_atkin2017_baserate: ',lb_params%maintresp_leaf_atkin2017_baserate
    print*,'fates_maintresp_leaf_ryan1991_baserate: ',lb_params%maintresp_leaf_ryan1991_baserate
    print*,'fates_leaf_vcmaxha: ',lb_params%vcmaxha
    print*,'fates_leaf_jmaxha: ',lb_params%jmaxha
    print*,'fates_leaf_vcmaxhd: ',lb_params%vcmaxhd
    print*,'fates_leaf_jmaxhd: ',lb_params%jmaxhd
    print*,'fates_leaf_vcmaxse: ',lb_params%vcmaxse
    print*,'fates_leaf_jmaxse: ',lb_params%jmaxse
    
  end subroutine DumpParams
  
end module LeafBiophysSuppMod
