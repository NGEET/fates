module FatesEdgeForestParamsMod
   !
   ! module that deals with reading the edge forest parameter file
   !
   use FatesConstantsMod,        only : r8 => fates_r8
   use FatesInterfaceTypesMod,  only : nlevedgeforest
   use FatesParametersInterface, only : param_string_length
   use FatesGlobals,             only : fates_log
   use FatesGlobals,             only : endrun => fates_endrun
   use FatesUtilsMod,            only : is_param_set
   use shr_log_mod,              only : errMsg => shr_log_errMsg

   implicit none
   private
   save

   !
   ! this is what the user can use for the actual values
   !
   real(r8),protected,allocatable,public :: ED_val_edgeforest_gaussian_amplitude(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_gaussian_sigma(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_gaussian_center(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_lognormal_amplitude(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_lognormal_sigma(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_lognormal_center(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_quadratic_a(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_quadratic_b(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_quadratic_c(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_rh_mult(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_temp_C_mult(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_wind_mult(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_rh_add(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_temp_C_add(:)
   real(r8),protected,allocatable,public :: ED_val_edgeforest_fireweather_wind_add(:)

   character(len=param_string_length),parameter,public :: ED_name_edgeforest_gaussian_amplitude = "fates_edgeforest_gaussian_amplitude"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_gaussian_sigma = "fates_edgeforest_gaussian_sigma"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_gaussian_center = "fates_edgeforest_gaussian_center"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_lognormal_amplitude = "fates_edgeforest_lognormal_amplitude"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_lognormal_sigma = "fates_edgeforest_lognormal_sigma"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_lognormal_center = "fates_edgeforest_lognormal_center"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_quadratic_a = "fates_edgeforest_quadratic_a"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_quadratic_b = "fates_edgeforest_quadratic_b"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_quadratic_c = "fates_edgeforest_quadratic_c"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_bin_edges = "fates_edgeforest_bin_edges"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_rh_mult = "fates_edgeforest_fireweather_rh_mult"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_temp_C_mult = "fates_edgeforest_fireweather_temp_C_mult"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_wind_mult = "fates_edgeforest_fireweather_wind_mult"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_rh_add = "fates_edgeforest_fireweather_rh_add"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_temp_C_add = "fates_edgeforest_fireweather_temp_C_add"
   character(len=param_string_length),parameter,public :: ED_name_edgeforest_fireweather_wind_add = "fates_edgeforest_fireweather_wind_add"

   character(len=*), parameter, private :: sourcefile =  __FILE__

   public :: EdgeForestRegisterParams
   public :: EdgeForestReceiveParams
   public :: EdgeForestCheckParams

contains

  ! =====================================================================================

  subroutine check_all_unset(value_array, b)
    real(r8), dimension(:), intent(in) :: value_array
    integer, intent(in) :: b
    integer :: i
    do i = 1, size(value_array)
       if (is_param_set(value_array(i))) then
          write(fates_log(),*) 'Multiple fit types found for bin ',b
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end do
  end subroutine check_all_unset

  subroutine EdgeForestCheckParams(is_master)

     ! ----------------------------------------------------------------------------------
     !
     ! This subroutine performs logical checks on user supplied parameters.  It cross
     ! compares various parameters and will fail if they don't make sense.
     ! E.g.: Each bin should have parameters for exactly one fit type.
     ! -----------------------------------------------------------------------------------

     logical, intent(in) :: is_master    ! Only check if this is the master proc

     real(r8) :: gaussian_amplitude
     real(r8) :: gaussian_sigma
     real(r8) :: gaussian_center
     real(r8) :: lognormal_amplitude
     real(r8) :: lognormal_sigma
     real(r8) :: lognormal_center
     real(r8) :: quadratic_a
     real(r8) :: quadratic_b
     real(r8) :: quadratic_c
     integer  :: b

     if (.not. is_master) return

     ! Check each bin
     do b = 1, nlevedgeforest

        gaussian_amplitude = ED_val_edgeforest_gaussian_amplitude(b)
        gaussian_sigma = ED_val_edgeforest_gaussian_sigma(b)
        gaussian_center = ED_val_edgeforest_gaussian_center(b)
        lognormal_amplitude = ED_val_edgeforest_lognormal_amplitude(b)
        lognormal_sigma = ED_val_edgeforest_lognormal_sigma(b)
        lognormal_center = ED_val_edgeforest_lognormal_center(b)
        quadratic_a = ED_val_edgeforest_quadratic_a(b)
        quadratic_b = ED_val_edgeforest_quadratic_b(b)
        quadratic_c = ED_val_edgeforest_quadratic_c(b)

        if (is_param_set(gaussian_amplitude)) then
           ! Has all gaussian parameters
           if (.not. (is_param_set(gaussian_center) .and. is_param_set(gaussian_sigma))) then
              write(fates_log(),*) 'Not all gaussian forest edge parameters found for bin ',b
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           ! Has no other parameters
           call check_all_unset( (/ lognormal_amplitude, lognormal_sigma, lognormal_center, quadratic_a, quadratic_b, quadratic_c /), b )

        else if (is_param_set(lognormal_amplitude)) then
           ! Has all lognormal parameters
           if (.not. (is_param_set(lognormal_center) .and. is_param_set(lognormal_sigma))) then
              write(fates_log(),*) 'Not all lognormal forest edge parameters found for bin ',b
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           ! Has no other parameters
           call check_all_unset( (/ gaussian_amplitude, gaussian_sigma, gaussian_center, quadratic_a, quadratic_b, quadratic_c /), b )

        else if (is_param_set(quadratic_a)) then
           ! Has all quadratic parameters
           if (.not. (is_param_set(quadratic_b) .and. is_param_set(quadratic_c))) then
              write(fates_log(),*) 'Not all quadratic forest edge parameters found for bin ',b
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           ! Has no other parameters
           call check_all_unset( (/ gaussian_amplitude, gaussian_sigma, gaussian_center, lognormal_amplitude, lognormal_sigma, lognormal_center /), b )

        else
           write(fates_log(),*) 'Unrecognized bin fit type for bin ',b
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
     end do

     return
  end subroutine EdgeForestCheckParams

  !-----------------------------------------------------------------------
  subroutine EdgeForestRegisterParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar
    use FatesParametersInterface, only : dimension_name_edgeforest_bins
    use FatesParametersInterface, only : dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_edgeforest(1)= (/dimension_name_edgeforest_bins/)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_gaussian_amplitude, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_gaussian_sigma, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_gaussian_center, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_lognormal_amplitude, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_lognormal_sigma, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_lognormal_center, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_quadratic_a, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_quadratic_b, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_quadratic_c, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_rh_mult, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_temp_C_mult, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_wind_mult, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_rh_add, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_temp_C_add, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

    call fates_params%RegisterParameter(name=ED_name_edgeforest_fireweather_wind_add, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_edgeforest)

 end subroutine EdgeForestRegisterParams

 !-----------------------------------------------------------------------
 subroutine EdgeForestReceiveParams(fates_params)

   use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

   implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_gaussian_amplitude, &
         data=ED_val_edgeforest_gaussian_amplitude)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_gaussian_sigma, &
         data=ED_val_edgeforest_gaussian_sigma)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_gaussian_center, &
         data=ED_val_edgeforest_gaussian_center)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_lognormal_amplitude, &
         data=ED_val_edgeforest_lognormal_amplitude)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_lognormal_sigma, &
         data=ED_val_edgeforest_lognormal_sigma)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_lognormal_center, &
         data=ED_val_edgeforest_lognormal_center)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_quadratic_a, &
         data=ED_val_edgeforest_quadratic_a)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_quadratic_b, &
         data=ED_val_edgeforest_quadratic_b)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_quadratic_c, &
         data=ED_val_edgeforest_quadratic_c)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_bin_edges, &
         data=ED_val_edgeforest_bin_edges)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_rh_mult, &
         data=ED_val_edgeforest_fireweather_rh_mult)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_temp_C_mult, &
         data=ED_val_edgeforest_fireweather_temp_C_mult)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_wind_mult, &
         data=ED_val_edgeforest_fireweather_wind_mult)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_rh_add, &
         data=ED_val_edgeforest_fireweather_rh_add)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_temp_C_add, &
         data=ED_val_edgeforest_fireweather_temp_C_add)

    call fates_params%RetrieveParameterAllocate(name=ED_name_edgeforest_fireweather_wind_add, &
         data=ED_val_edgeforest_fireweather_wind_add)

  end subroutine EdgeForestReceiveParams


end module FatesEdgeForestParamsMod
