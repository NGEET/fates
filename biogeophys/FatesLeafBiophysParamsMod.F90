module FatesLeafBiophysParamsMod

  use FatesConstantsMod , only: r8 => fates_r8
  use FatesConstantsMod , only: fates_check_param_set
  use FatesGlobals,   only : fates_log
  use FatesGlobals,   only : endrun => fates_endrun
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use LeafBiophysicsMod, only : lb_params,btran_on_gs_gs1,btran_on_ag_none
  use JSONParameterUtilsMod,only : params_type,param_type
  
  implicit none
  private ! Modules are private by default
  save


  public :: TransferParamsLeafBiophys
  public :: LeafBiophysReportParams
  
  character(len=*), parameter :: sourcefile = &
        __FILE__

contains

  ! =====================================================================================

  subroutine TransferParamsLeafBiophys(pstruct)

    type(params_type) :: pstruct         ! Data structure containing all parameters and dimensions
    type(param_type),pointer :: param_p  ! Pointer to one specific parameter
    integer                  :: numpft
    
    numpft = pstruct%GetDimSizeFromName('fates_pft')
    
    param_p => pstruct%GetParamFromName('fates_leaf_c3psn')
    allocate(lb_params%c3psn(numpft))
    lb_params%c3psn(:) = param_p%i_data_1d(:)

    param_p => pstruct%GetParamFromName('fates_leaf_stomatal_btran_model')
    allocate(lb_params%stomatal_btran_model(numpft))
    lb_params%stomatal_btran_model(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_agross_btran_model')
    allocate(lb_params%agross_btran_model(numpft))
    lb_params%agross_btran_model(:) = param_p%i_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_stomatal_slope_medlyn')
    allocate(lb_params%medlyn_slope(numpft))
    lb_params%medlyn_slope(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_stomatal_slope_ballberry')
    allocate(lb_params%bb_slope(numpft))
    lb_params%bb_slope(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_stomatal_intercept')
    allocate(lb_params%stomatal_intercept(numpft))
    lb_params%stomatal_intercept(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_leaf_ryan1991_baserate')
    allocate(lb_params%maintresp_leaf_ryan1991_baserate(numpft))
    lb_params%maintresp_leaf_ryan1991_baserate(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_leaf_atkin2017_baserate')
    allocate(lb_params%maintresp_leaf_atkin2017_baserate(numpft))
    lb_params%maintresp_leaf_atkin2017_baserate(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_reduction_curvature')
    allocate(lb_params%maintresp_reduction_curvature(numpft))
    lb_params%maintresp_reduction_curvature(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_reduction_intercept')
    allocate(lb_params%maintresp_reduction_intercept(numpft))
    lb_params%maintresp_reduction_intercept(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_maintresp_reduction_upthresh')
    allocate(lb_params%maintresp_reduction_upthresh(numpft))
    lb_params%maintresp_reduction_upthresh(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_vcmaxha')
    allocate(lb_params%vcmaxha(numpft))
    lb_params%vcmaxha(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_jmaxha')
    allocate(lb_params%jmaxha(numpft))
    lb_params%jmaxha(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_vcmaxhd')
    allocate(lb_params%vcmaxhd(numpft))
    lb_params%vcmaxhd(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_jmaxhd')
    allocate(lb_params%jmaxhd(numpft))
    lb_params%jmaxhd(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_vcmaxse')
    allocate(lb_params%vcmaxse(numpft))
    lb_params%vcmaxse(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_jmaxse')
    allocate(lb_params%jmaxse(numpft))
    lb_params%jmaxse(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName('fates_leaf_fnps')
    allocate(lb_params%fnps(numpft))
    lb_params%fnps(:) = param_p%r_data_1d(:)
    
    return
  end subroutine TransferParamsLeafBiophys

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
       write(fates_log(),fmt_rout) 'fates_leaf_fnps = ',lb_params%fnps
       write(fates_log(),fmt_rout) 'nl: electron_transport_model = ',lb_params%electron_transport_model
       write(fates_log(),fmt_iout) 'nl: daylength_factor_switch = ',lb_params%dayl_switch
       write(fates_log(),fmt_iout) 'nl: leaf_stomatal_model = ',lb_params%stomatal_model
       write(fates_log(),fmt_iout) 'nl: leaf_stomatal_assim_model = ',lb_params%stomatal_assim_model
       write(fates_log(),fmt_iout) 'nl: leaf_photo_tempsens_model = ',lb_params%photo_tempsens_model
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
