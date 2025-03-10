module FatesLeafBiophysParamsMod

  use FatesConstantsMod , only: r8 => fates_r8
  use FatesConstantsMod , only: fates_check_param_set
  use FatesParametersInterface, only : param_string_length
  use FatesGlobals,   only : fates_log
  use FatesGlobals,   only : endrun => fates_endrun
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use LeafBiophysicsMod, only : lb_params,btran_on_gs_gs1,btran_on_ag_none
  use FatesParametersInterface, only : fates_parameters_type
  ! Register the parameters we want the host to provide, and
  ! indicate whether they are fates parameters or host parameters
  ! that need to be synced with host values.
  use FatesParametersInterface, only : fates_parameters_type, param_string_length
  use FatesParametersInterface, only : dimension_name_pft, dimension_shape_1d
  use FatesParametersInterface, only : dimension_shape_scalar, dimension_name_scalar
  use FatesUtilsMod, only : ArrayNint
  
  implicit none
  private ! Modules are private by default
  save

  public :: LeafBiophysRegisterParams
  public :: LeafBiophysReceiveParams
  public :: LeafBiophysReportParams
  
  character(len=*), parameter :: sourcefile = &
        __FILE__

  integer, parameter  :: lower_bound_pft = 1

contains

  ! =====================================================================================

  subroutine LeafBiophysRegisterParams(fates_params)


    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_pft/)
    integer, parameter :: dim_lower_bound(1) = (/ lower_bound_pft /)
    character(len=param_string_length) :: name
    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)


    ! Register PFT dimensioned
    
    name = 'fates_leaf_c3psn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_stomatal_btran_model'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_agross_btran_model'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_leaf_stomatal_slope_ballberry'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_stomatal_slope_medlyn'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)
    
    name = 'fates_leaf_stomatal_intercept'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_reduction_curvature'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_reduction_intercept'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_reduction_upthresh'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_leaf_atkin2017_baserate'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_maintresp_leaf_ryan1991_baserate'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_vcmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxha'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_vcmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxhd'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_vcmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    name = 'fates_leaf_jmaxse'
    call fates_params%RegisterParameter(name=name, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names, lower_bounds=dim_lower_bound)

    return
  end subroutine LeafBiophysRegisterParams

  ! =====================================================================================

  subroutine LeafBiophysReceiveParams(fates_params)

    !use FatesInterfaceTypesMod, only : hlm_daylength_factor_switch
    !use FatesInterfaceTypesMod, only : hlm_stomatal_model
    !use FatesInterfaceTypesMod, only : hlm_stomatal_assim_model
    !use FatesInterfaceTypesMod, only : hlm_photo_tempsens_model
    
    class(fates_parameters_type), intent(inout) :: fates_params
    real(r8), allocatable :: tmpreal(:)  ! Temporary variable to hold floats
    real(r8)              :: tmpscalar
    character(len=param_string_length) :: name

    !lb_params%dayl_switch    = hlm_daylength_factor_switch
    !lb_params%stomatal_model = hlm_stomatal_model
    !lb_params%stomatal_assim_model = hlm_stomatal_assim_model
    !lb_params%photo_tempsens_model = hlm_photo_tempsens_model

    name = 'fates_leaf_c3psn'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(lb_params%c3psn(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,lb_params%c3psn)
    deallocate(tmpreal)

    name = 'fates_leaf_stomatal_btran_model'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(lb_params%stomatal_btran_model(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,lb_params%stomatal_btran_model)
    deallocate(tmpreal)

    name = 'fates_leaf_agross_btran_model'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=tmpreal)
    allocate(lb_params%agross_btran_model(size(tmpreal,dim=1)))
    call ArrayNint(tmpreal,lb_params%agross_btran_model)
    deallocate(tmpreal)
    
    name = 'fates_leaf_stomatal_slope_medlyn'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%medlyn_slope)

    name = 'fates_leaf_stomatal_slope_ballberry'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%bb_slope)

    name = 'fates_leaf_stomatal_intercept'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%stomatal_intercept)

    name = 'fates_maintresp_leaf_ryan1991_baserate'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%maintresp_leaf_ryan1991_baserate)

    name = 'fates_maintresp_leaf_atkin2017_baserate'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%maintresp_leaf_atkin2017_baserate)

    name = 'fates_maintresp_reduction_curvature'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%maintresp_reduction_curvature)

    name = 'fates_maintresp_reduction_intercept'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%maintresp_reduction_intercept)

    name = 'fates_maintresp_reduction_upthresh'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%maintresp_reduction_upthresh)

    name = 'fates_leaf_vcmaxha'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%vcmaxha)

    name = 'fates_leaf_jmaxha'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%jmaxha)

    name = 'fates_leaf_vcmaxhd'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%vcmaxhd)

    name = 'fates_leaf_jmaxhd'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%jmaxhd)

    name = 'fates_leaf_vcmaxse'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%vcmaxse)

    name = 'fates_leaf_jmaxse'
    call fates_params%RetrieveParameterAllocate(name=name, &
         data=lb_params%jmaxse)


    return
  end subroutine LeafBiophysReceiveParams

  ! ====================================================================================

  subroutine LeafBiophysReportParams(is_master)

    ! Argument
    logical, intent(in) :: is_master  ! Only log if this is the master proc

    logical, parameter :: debug_report = .false.
    character(len=32),parameter :: fmt_rout = '(a,F16.8)'
    character(len=32),parameter :: fmt_iout = '(a,I8)'
    
    integer :: npft,ipft

    if(debug_report .and. is_master) then
       write(fates_log(),fmt_iout) 'fates_leaf_c3psn = ',lb_params%c3psn
       write(fates_log(),fmt_iout) 'fates_leaf_stomatal_btran_model = ',lb_params%stomatal_btran_model
       write(fates_log(),fmt_iout) 'fates_leaf_agross_btran_model = ',lb_params%agross_btran_model
       write(fates_log(),fmt_rout) 'fates_leaf_vcmaxha = ',lb_params%vcmaxha
       write(fates_log(),fmt_rout) 'fates_leaf_jmaxha = ',lb_params%jmaxha
       write(fates_log(),fmt_rout) 'fates_leaf_vcmaxhd = ',lb_params%vcmaxhd
       write(fates_log(),fmt_rout) 'fates_leaf_jmaxhd = ',lb_params%jmaxhd
       write(fates_log(),fmt_rout) 'fates_leaf_vcmaxse = ',lb_params%vcmaxse
       write(fates_log(),fmt_rout) 'fates_leaf_jmaxse = ',lb_params%jmaxse
       write(fates_log(),fmt_iout) 'fates_daylength_factor_switch = ',lb_params%dayl_switch
       write(fates_log(),fmt_iout) 'fates_leaf_stomatal_model = ',lb_params%stomatal_model
       write(fates_log(),fmt_iout) 'fates_leaf_stomatal_assim_model = ',lb_params%stomatal_assim_model
       write(fates_log(),fmt_iout) 'fates_leaf_photo_tempsens_model = ',lb_params%photo_tempsens_model
       write(fates_log(),fmt_rout) 'fates_leaf_stomatal_slope_medlyn = ',lb_params%medlyn_slope
       write(fates_log(),fmt_rout) 'fates_leaf_stomatal_slope_ballberry = ',lb_params%bb_slope
       write(fates_log(),fmt_rout) 'fates_leaf_stomatal_intercept = ',lb_params%stomatal_intercept
       write(fates_log(),fmt_rout) 'fates_maintresp_leaf_ryan1991_baserate = ',lb_params%maintresp_leaf_ryan1991_baserate
       write(fates_log(),fmt_rout) 'fates_maintresp_leaf_atkin2017_baserate = ',lb_params%maintresp_leaf_atkin2017_baserate
       write(fates_log(),fmt_rout) 'fates_maintresp_reduction_curvature = ',lb_params%maintresp_reduction_curvature
       write(fates_log(),fmt_rout) 'fates_maintresp_reduction_intercept = ',lb_params%maintresp_reduction_intercept
       write(fates_log(),fmt_rout) 'fates_maintresp_reduction_upthresh = ',lb_params%maintresp_reduction_upthresh
    end if

    
  end subroutine LeafBiophysReportParams

end module FatesLeafBiophysParamsMod
