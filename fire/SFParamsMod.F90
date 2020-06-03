module SFParamsMod
   !
   ! module that deals with reading the SF parameter file
   !
   use FatesConstantsMod , only: r8 => fates_r8
   use FatesConstantsMod , only: fates_check_param_set
   use EDtypesMod        , only: NFSC
   use FatesLitterMod    , only: ncwd
   use FatesParametersInterface, only : param_string_length
   use FatesGlobals,   only : fates_log
   use FatesGlobals,   only : endrun => fates_endrun
   use shr_log_mod      , only : errMsg => shr_log_errMsg

   implicit none
   private ! Modules are private by default
   save

   !
   ! this is what the user can use for the actual values
   !

   real(r8),protected, public :: SF_val_fdi_a
   real(r8),protected, public :: SF_val_fdi_b
   real(r8),protected, public :: SF_val_fdi_alpha
   real(r8),protected, public :: SF_val_miner_total
   real(r8),protected, public :: SF_val_fuel_energy
   real(r8),protected, public :: SF_val_part_dens
   real(r8),protected, public :: SF_val_miner_damp
   real(r8),protected, public :: SF_val_max_durat
   real(r8),protected, public :: SF_val_durat_slope
   real(r8),protected, public :: SF_val_drying_ratio
   real(r8),protected, public :: SF_val_fire_threshold    ! threshold for fires that spread or go out. kW/m (Pyne 1996)
   real(r8),protected, public :: SF_val_CWD_frac(ncwd)
   real(r8),protected, public :: SF_val_max_decomp(NFSC)
   real(r8),protected, public :: SF_val_SAV(NFSC)
   real(r8),protected, public :: SF_val_FBD(NFSC)
   real(r8),protected, public :: SF_val_min_moisture(NFSC)
   real(r8),protected, public :: SF_val_mid_moisture(NFSC)
   real(r8),protected, public :: SF_val_low_moisture_Coeff(NFSC)
   real(r8),protected, public :: SF_val_low_moisture_Slope(NFSC)
   real(r8),protected, public :: SF_val_mid_moisture_Coeff(NFSC)
   real(r8),protected, public :: SF_val_mid_moisture_Slope(NFSC)

   character(len=param_string_length),parameter :: SF_name_fdi_a = "fates_fire_fdi_a"
   character(len=param_string_length),parameter :: SF_name_fdi_b = "fates_fire_fdi_b"
   character(len=param_string_length),parameter :: SF_name_fdi_alpha = "fates_fire_fdi_alpha"
   character(len=param_string_length),parameter :: SF_name_miner_total = "fates_fire_miner_total"
   character(len=param_string_length),parameter :: SF_name_fuel_energy = "fates_fire_fuel_energy"
   character(len=param_string_length),parameter :: SF_name_part_dens = "fates_fire_part_dens"
   character(len=param_string_length),parameter :: SF_name_miner_damp = "fates_fire_miner_damp"
   character(len=param_string_length),parameter :: SF_name_max_durat = "fates_fire_max_durat"
   character(len=param_string_length),parameter :: SF_name_durat_slope = "fates_fire_durat_slope"
   character(len=param_string_length),parameter :: SF_name_drying_ratio = "fates_fire_drying_ratio"
   character(len=param_string_length),parameter :: SF_name_fire_threshold = "fates_fire_threshold"
   character(len=param_string_length),parameter :: SF_name_CWD_frac = "fates_CWD_frac"
   character(len=param_string_length),parameter :: SF_name_max_decomp = "fates_max_decomp"
   character(len=param_string_length),parameter :: SF_name_SAV = "fates_fire_SAV"
   character(len=param_string_length),parameter :: SF_name_FBD = "fates_fire_FBD"
   character(len=param_string_length),parameter :: SF_name_min_moisture = "fates_fire_min_moisture"
   character(len=param_string_length),parameter :: SF_name_mid_moisture = "fates_fire_mid_moisture"
   character(len=param_string_length),parameter :: SF_name_low_moisture_Coeff = "fates_fire_low_moisture_Coeff"
   character(len=param_string_length),parameter :: SF_name_low_moisture_Slope = "fates_fire_low_moisture_Slope"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_Coeff = "fates_fire_mid_moisture_Coeff"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_Slope = "fates_fire_mid_moisture_Slope"

   character(len=*), parameter, private :: sourcefile = &
         __FILE__


   real(r8), parameter,private :: min_fire_threshold = 0.0001_r8  ! The minimum reasonable fire intensity threshold [kW/m]


   public :: SpitFireRegisterParams
   public :: SpitFireReceiveParams
   public :: SpitFireCheckParams


contains

  ! =====================================================================================

  subroutine SpitFireCheckParams(is_master)

     ! ----------------------------------------------------------------------------------
     !
     ! This subroutine performs logical checks on user supplied parameters.  It cross
     ! compares various parameters and will fail if they don't make sense.  
     ! Examples:
     ! Decomposition rates should not be less than zero or greater than 1
     ! -----------------------------------------------------------------------------------

     logical, intent(in) :: is_master    ! Only log if this is the master proc


     integer :: c      ! debris type loop counter
     integer :: corr_id(1)        ! This is the bin with largest fraction
                                  ! add/subtract any corrections there
     real(r8) :: correction       ! This correction ensures that root fractions
                                  ! sum to 1.0


     if(.not.is_master) return
     
     ! Move these checks to initialization
     do c = 1,nfsc
        if ( SF_val_max_decomp(c) < 0._r8) then
           write(fates_log(),*) 'Decomposition rates should be >0'
           write(fates_log(),*) 'c = ',c,' SF_val_max_decomp(c) = ',SF_val_max_decomp(c)
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end if
     end do

     ! Check if the CWD fraction sums to unity, if it is not wayyy off,
     ! add a small correction to the largest pool. 
     ! This is important for tight mass conservation
     ! checks

     if(abs(1.0_r8 - sum(SF_val_CWD_frac(1:ncwd))) > 1.e-5_r8) then
         write(fates_log(),*) 'The CWD fractions from index 1:4 must sum to unity'
         write(fates_log(),*) 'SF_val_CWD_frac(1:ncwd) = ',SF_val_CWD_frac(1:ncwd)
         write(fates_log(),*) 'error = ',1.0_r8 - sum(SF_val_CWD_frac(1:ncwd))
         call endrun(msg=errMsg(sourcefile, __LINE__))
     else
         correction = 1._r8 - sum(SF_val_CWD_frac(1:ncwd))
         corr_id = maxloc(SF_val_CWD_frac(1:ncwd))
         SF_val_CWD_frac(corr_id(1)) = SF_val_CWD_frac(corr_id(1)) + correction
     end if

     ! Check to see if the fire threshold is above the minimum and set at all
     if(SF_val_fire_threshold < min_fire_threshold .or. & 
           SF_val_fire_threshold >  fates_check_param_set ) then
         write(fates_log(),*) 'The fates_fire_threshold parameter must be set, and > ',min_fire_threshold
         write(fates_log(),*) 'The value is set at :',SF_val_fire_threshold
         write(fates_log(),*) 'Please provide a reasonable value, aborting.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
     end if


     return
  end subroutine SpitFireCheckParams

  !-----------------------------------------------------------------------
  subroutine SpitFireParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    SF_val_fdi_a = nan
    SF_val_fdi_b = nan
    SF_val_fdi_alpha = nan
    SF_val_miner_total = nan
    SF_val_fuel_energy = nan
    SF_val_part_dens = nan
    SF_val_miner_damp = nan
    SF_val_max_durat = nan
    SF_val_durat_slope = nan
    SF_val_drying_ratio = nan
    SF_val_fire_threshold = nan
    SF_val_CWD_frac(:) = nan
    SF_val_max_decomp(:) = nan
    SF_val_SAV(:) = nan
    SF_val_FBD(:) = nan
    SF_val_min_moisture(:) = nan
    SF_val_mid_moisture(:) = nan
    SF_val_low_moisture_Coeff(:) = nan
    SF_val_low_moisture_Slope(:) = nan
    SF_val_mid_moisture_Coeff(:) = nan
    SF_val_mid_moisture_Slope(:) = nan

  end subroutine SpitFireParamsInit

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call SpitFireParamsInit()
    call SpitFireRegisterScalars(fates_params)
    call SpitFireRegisterNCWD(fates_params)
    call SpitFireRegisterNFSC(fates_params)

 end subroutine SpitFireRegisterParams

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveParams(fates_params)
   
   use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar
   
   implicit none
   
    class(fates_parameters_type), intent(inout) :: fates_params
    
    call SpitFireReceiveScalars(fates_params)
    call SpitFireReceiveNCWD(fates_params)
    call SpitFireReceiveNFSC(fates_params)
    
  end subroutine SpitFireReceiveParams

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterScalars(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)
    
    call fates_params%RegisterParameter(name=SF_name_fdi_a, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fdi_b, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fdi_alpha, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_miner_total, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fuel_energy, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_part_dens, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_miner_damp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_max_durat, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_durat_slope, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_drying_ratio, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fire_threshold, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

  end subroutine SpitFireRegisterScalars

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveScalars(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params
    real(r8) :: tmp_real
    

    call fates_params%RetreiveParameter(name=SF_name_fdi_a, &
         data=SF_val_fdi_a)

    call fates_params%RetreiveParameter(name=SF_name_fdi_b, &
         data=SF_val_fdi_b)

    call fates_params%RetreiveParameter(name=SF_name_fdi_alpha, &
         data=SF_val_fdi_alpha)

    call fates_params%RetreiveParameter(name=SF_name_miner_total, &
         data=SF_val_miner_total)

    call fates_params%RetreiveParameter(name=SF_name_fuel_energy, &
         data=SF_val_fuel_energy)

    call fates_params%RetreiveParameter(name=SF_name_part_dens, &
         data=SF_val_part_dens)

    call fates_params%RetreiveParameter(name=SF_name_miner_damp, &
         data=SF_val_miner_damp)

    call fates_params%RetreiveParameter(name=SF_name_max_durat, &
         data=SF_val_max_durat)

    call fates_params%RetreiveParameter(name=SF_name_durat_slope, &
         data=SF_val_durat_slope)

    call fates_params%RetreiveParameter(name=SF_name_drying_ratio, &
         data=SF_val_drying_ratio)

    call fates_params%RetreiveParameter(name=SF_name_fire_threshold, &
         data=SF_val_fire_threshold)




  end subroutine SpitFireReceiveScalars

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterNCWD(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_cwd, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_cwd(1) = (/dimension_name_cwd/)

    call fates_params%RegisterParameter(name=SF_name_CWD_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_cwd)

  end subroutine SpitFireRegisterNCWD

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveNCWD(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetreiveParameter(name=SF_name_CWD_frac, &
         data=SF_val_CWD_frac)

    

  end subroutine SpitFireReceiveNCWD

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterNFSC(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_fsc, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_fsc/)

    call fates_params%RegisterParameter(name=SF_name_SAV, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_FBD, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_min_moisture, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_low_moisture_Coeff, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_low_moisture_Slope, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture_Coeff, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture_Slope, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_max_decomp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

  end subroutine SpitFireRegisterNFSC

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveNFSC(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    
    call fates_params%RetreiveParameter(name=SF_name_SAV, &
         data=SF_val_SAV)

    call fates_params%RetreiveParameter(name=SF_name_FBD, &
         data=SF_val_FBD)

    call fates_params%RetreiveParameter(name=SF_name_min_moisture, &
         data=SF_val_min_moisture)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture, &
         data=SF_val_mid_moisture)

    call fates_params%RetreiveParameter(name=SF_name_low_moisture_Coeff, &
         data=SF_val_low_moisture_Coeff)

    call fates_params%RetreiveParameter(name=SF_name_low_moisture_Slope, &
         data=SF_val_low_moisture_Slope)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture_Coeff, &
         data=SF_val_mid_moisture_Coeff)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture_Slope, &
         data=SF_val_mid_moisture_Slope)

    call fates_params%RetreiveParameter(name=SF_name_max_decomp, &
         data=SF_val_max_decomp)

  end subroutine SpitFireReceiveNFSC
  !-----------------------------------------------------------------------


end module SFParamsMod
