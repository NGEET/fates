module SFParamsMod
   !
   ! module that deals with reading the SF parameter file
   !
   use FatesConstantsMod,        only : r8 => fates_r8
   use FatesConstantsMod,        only : fates_check_param_set
   use FatesFuelClassesMod,      only : num_fuel_classes
   use FatesLitterMod,           only : ncwd
   use FatesGlobals,             only : fates_log
   use FatesGlobals,             only : endrun => fates_endrun
   use shr_log_mod,              only : errMsg => shr_log_errMsg
   use JSONParameterUtilsMod,only : params_type,param_type
   
   implicit none
   private 
   save

   !
   ! this is what the user can use for the actual values
   !
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
   real(r8),protected, public :: SF_val_max_decomp(num_fuel_classes)
   real(r8),protected, public :: SF_val_SAV(num_fuel_classes)
   real(r8),protected, public :: SF_val_FBD(num_fuel_classes)
   real(r8),protected, public :: SF_val_min_moisture(num_fuel_classes)
   real(r8),protected, public :: SF_val_mid_moisture(num_fuel_classes)
   real(r8),protected, public :: SF_val_low_moisture_Coeff(num_fuel_classes)
   real(r8),protected, public :: SF_val_low_moisture_Slope(num_fuel_classes)
   real(r8),protected, public :: SF_val_mid_moisture_Coeff(num_fuel_classes)
   real(r8),protected, public :: SF_val_mid_moisture_Slope(num_fuel_classes)
    ! Prescribed fire relevant parameters
   real(r8),protected, public :: SF_val_rxfire_tpup   ! temperature upper threshold above which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_tplw   ! temperature lower threshold below which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_rhup   ! relative humidity upper threshold above which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_rhlw   ! relative humidity lower threshold below which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_wdup   ! wind speed upper threshold above which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_wdlw   ! wind speed lower threshold below which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_AB     ! prescribed fire burned fraction per day
   real(r8),protected, public :: SF_val_rxfire_min_threshold ! minimum fire energy at or below which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_max_threshold ! maximum fire energy at or above which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_fuel_min      ! minimum fuel load at or below which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_fuel_max      ! maximum fuel load at or above which rx fire is disallowed
   real(r8),protected, public :: SF_val_rxfire_min_frac      ! minimum burnable fraction at site level at or above which rx fire is allowed

   character(len=*), parameter, private :: sourcefile =  __FILE__
   real(r8),         parameter, private :: min_fire_threshold = 0.0001_r8  ! The minimum reasonable fire intensity threshold [kW/m]

   public :: TransferParamsSpitFire
   
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
     do c = 1,num_fuel_classes
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
    SF_val_rxfire_tpup = nan
    SF_val_rxfire_tplw = nan
    SF_val_rxfire_rhup = nan
    SF_val_rxfire_rhlw = nan
    SF_val_rxfire_wdup = nan
    SF_val_rxfire_wdlw = nan
    SF_val_rxfire_AB   = nan
    SF_val_rxfire_min_threshold = nan
    SF_val_rxfire_max_threshold = nan
    SF_val_rxfire_fuel_min = nan
    SF_val_rxfire_fuel_max = nan
    SF_val_rxfire_min_frac = nan

  end subroutine SpitFireParamsInit

  ! ======================================================================

  subroutine TransferParamsSpitFire(pstruct)

    ! -----------------------------------------------------------------------------------
    ! Transfer SpitFire parameter values from the data structure "pstruct"
    ! to named primitive data structures
    ! -----------------------------------------------------------------------------------

    implicit none

    type(params_type) :: pstruct         ! Data structure containing all parameters and dimensions
    type(param_type),pointer :: param_p  ! Pointer to one specific parameter

    call SpitFireParamsInit()

    param_p => pstruct%GetParamFromName("fates_fire_fdi_alpha")
    SF_val_fdi_alpha = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_miner_total")
    SF_val_miner_total = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_fuel_energy")
    SF_val_fuel_total = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_part_dens")
    SF_val_part_dens = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_miner_damp")
    SF_val_miner_damp = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_max_durat")
    SF_val_max_durat = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_durat_slope")
    SF_val_durat_slope = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_drying_ratio")
    SF_val_drying_ratio = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_fire_threshold")
    SF_val_fire_threshold = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_frag_cwd_frac")
    SF_val_CWD_frac(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName("fates_frag_maxdecomp")
    SF_val_max_decomp(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_SAV")
    SF_val_SAV(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_FBD")
    SF_val_FBD(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_min_moisture")
    SF_val_min_moisture(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_mid_moisture")
    SF_val_mid_moisture(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_low_moisture_Coeff")
    SF_val_low_moisture_Coeff(:) = param_p%r_data_1d(:)

    param_p => pstruct%GetParamFromName("fates_fire_low_moisture_Slope")
    SF_val_low_moisture_Slope(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_mid_moisture_Coeff")
    SF_val_mid_moisture_Coeff(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_fire_mid_moisture_Slope")
    SF_val_mid_moisture_Slope(:) = param_p%r_data_1d(:)
    
    param_p => pstruct%GetParamFromName("fates_rxfire_temp_upthreshold")
    SF_val_rxfire_tpup = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_temp_lwthreshold")
    SF_val_rxfire_tplw = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_rh_upthreshold")
    SF_val_rxfire_rhup = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_rh_lwthreshold")
    SF_val_rxfire_rhlw = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_wind_upthreshold")
    SF_val_rxfire_wdup = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_wind_lwthreshold")
    SF_val_rxfire_wdlw = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_AB")
    SF_val_rxfire_AB = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_min_threshold")
    SF_val_rxfire_min_threshold = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_max_threshold")
    SF_val_rxfire_max_threshold = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_fuel_min")
    SF_val_rxfire_fuel_min = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_fuel_max")
    SF_val_rxfire_fuel_max = param_p%r_data_scalar
    
    param_p => pstruct%GetParamFromName("fates_rxfire_min_frac")
    SF_val_rxfire_min_frac = param_p%r_data_scalar
    
  end subroutine TransferParamsSpitFire

end module SFParamsMod
