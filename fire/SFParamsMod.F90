module SFParamsMod
   !
   ! module that deals with reading the SF parameter file
   !
   use shr_kind_mod      , only: r8 => shr_kind_r8
   use EDtypesMod        , only: NLSC,NFSC,NCWD
   use FatesParametersInterface, only : param_string_length

   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   real(r8),protected :: SF_val_fdi_a
   real(r8),protected :: SF_val_fdi_b
   real(r8),protected :: SF_val_fdi_alpha
   real(r8),protected :: SF_val_miner_total
   real(r8),protected :: SF_val_fuel_energy
   real(r8),protected :: SF_val_part_dens
   real(r8),protected :: SF_val_miner_damp
   real(r8),protected :: SF_val_max_durat
   real(r8),protected :: SF_val_durat_slope
   real(r8),protected :: SF_val_alpha_SH
   real(r8),protected :: SF_val_alpha_FMC(NLSC)
   real(r8),protected :: SF_val_CWD_frac(NCWD)
   real(r8),protected :: SF_val_max_decomp(NLSC)
   real(r8),protected :: SF_val_SAV(NFSC)
   real(r8),protected :: SF_val_FBD(NFSC)
   real(r8),protected :: SF_val_min_moisture(NFSC)
   real(r8),protected :: SF_val_mid_moisture(NFSC)
   real(r8),protected :: SF_val_low_moisture_C(NFSC)
   real(r8),protected :: SF_val_low_moisture_S(NFSC)
   real(r8),protected :: SF_val_mid_moisture_C(NFSC)
   real(r8),protected :: SF_val_mid_moisture_S(NFSC)
  
   character(len=param_string_length),parameter :: SF_name_fdi_a = "fdi_a"
   character(len=param_string_length),parameter :: SF_name_fdi_b = "fdi_b"
   character(len=param_string_length),parameter :: SF_name_fdi_alpha = "fdi_alpha"
   character(len=param_string_length),parameter :: SF_name_miner_total = "miner_total"
   character(len=param_string_length),parameter :: SF_name_fuel_energy = "fuel_energy"
   character(len=param_string_length),parameter :: SF_name_part_dens = "part_dens"
   character(len=param_string_length),parameter :: SF_name_miner_damp = "miner_damp"
   character(len=param_string_length),parameter :: SF_name_max_durat = "max_durat"
   character(len=param_string_length),parameter :: SF_name_durat_slope = "durat_slope"
   character(len=param_string_length),parameter :: SF_name_alpha_SH = "alpha_SH"
   character(len=param_string_length),parameter :: SF_name_alpha_FMC = "alpha_FMC"
   character(len=param_string_length),parameter :: SF_name_CWD_frac = "CWD_frac"
   character(len=param_string_length),parameter :: SF_name_max_decomp = "max_decomp"
   character(len=param_string_length),parameter :: SF_name_SAV = "SAV"
   character(len=param_string_length),parameter :: SF_name_FBD = "FBD"
   character(len=param_string_length),parameter :: SF_name_min_moisture = "min_moisture"
   character(len=param_string_length),parameter :: SF_name_mid_moisture = "mid_moisture"
   character(len=param_string_length),parameter :: SF_name_low_moisture_C = "low_moisture_C"
   character(len=param_string_length),parameter :: SF_name_low_moisture_S = "low_moisture_S"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_C = "mid_moisture_C"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_S = "mid_moisture_S"

   public :: SFParamsRead
   public :: SpitFireParamsInit
   public :: SpitFireRegisterParams
   public :: SpitFireReceiveParams
  
contains
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
    SF_val_alpha_SH = nan
    SF_val_alpha_FMC(:) = nan
    SF_val_CWD_frac(:) = nan
    SF_val_max_decomp(:) = nan
    SF_val_SAV(:) = nan
    SF_val_FBD(:) = nan
    SF_val_min_moisture(:) = nan
    SF_val_mid_moisture(:) = nan
    SF_val_low_moisture_C(:) = nan
    SF_val_low_moisture_S(:) = nan
    SF_val_mid_moisture_C(:) = nan
    SF_val_mid_moisture_S(:) = nan

  end subroutine SpitFireParamsInit
  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)

    !call SpitFireParamsInit()
    

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

    call fates_params%RegisterParameter(name=SF_name_alpha_SH, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

 end subroutine SpitFireRegisterParams

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveParams(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    
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

    call fates_params%RetreiveParameter(name=SF_name_alpha_SH, &
         data=SF_val_alpha_SH)

  end subroutine SpitFireReceiveParams
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine SFParamsRead(ncid)
     !
     ! calls to initialize parameter instance and do ncdio read
     !
     use ncdio_pio    , only : file_desc_t
     
     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     call SpitFireParamsInit()
     call SFParamsReadLocal(ncid)

  end subroutine SFParamsRead
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine SFParamsReadLocal(ncid)
     !
     ! read the netcdf file and populate internalInstScalar
     !
     use ncdio_pio         , only : file_desc_t
     use paramUtilMod      , only : readNcdio

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     ! local vars
     character(len=32)  :: subname = 'SFParamsReadLocal::'

     !
     ! call read function
     !

     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_fdi_a, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_fdi_a)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_fdi_b, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_fdi_b)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_fdi_alpha, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_fdi_alpha)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_miner_total, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_miner_total)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_fuel_energy, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_fuel_energy)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_part_dens, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_part_dens)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_miner_damp, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_miner_damp)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_max_durat, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_max_durat)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_durat_slope, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_durat_slope)
     !X!
     !X! call readNcdio(ncid = ncid, &
     !X!    varName=SF_name_alpha_SH, &
     !X!    callingName=subname, &
     !X!    retVal=SF_val_alpha_SH)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_alpha_FMC, &
         callingName=subname, &
         retVal=SF_val_alpha_FMC)

      call readNcdio(ncid = ncid, &
         varName=SF_name_CWD_frac, &
         callingName=subname, &
         retVal=SF_val_CWD_frac)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_max_decomp, &
         callingName=subname, &
         retVal=SF_val_max_decomp)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_SAV, &
         callingName=subname, &
         retVal=SF_val_SAV)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_FBD, &
         callingName=subname, &
         retVal=SF_val_FBD)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_min_moisture, &
         callingName=subname, &
         retVal=SF_val_min_moisture)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture, &
         callingName=subname, &
         retVal=SF_val_mid_moisture)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_low_moisture_C, &
         callingName=subname, &
         retVal=SF_val_low_moisture_C)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_low_moisture_S, &
         callingName=subname, &
         retVal=SF_val_low_moisture_S)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture_C, &
         callingName=subname, &
         retVal=SF_val_mid_moisture_C)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture_S, &
         callingName=subname, &
         retVal=SF_val_mid_moisture_S)
  
  end subroutine SFParamsReadLocal
  !-----------------------------------------------------------------------

end module SFParamsMod
