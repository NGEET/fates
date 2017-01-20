module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  ! ! USES:
  use clm_varctl , only : paramfile, iulog, use_ed, use_cn
  use spmdMod    , only : masterproc
  use fileutils  , only : getfil
  use ncdio_pio  , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio  , only : file_desc_t , ncd_inqdid, ncd_inqdlen

  implicit none
  private
  !
  public :: readParameters
  private :: readFatesParameters
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParameters (nutrient_competition_method, photosyns_inst)
    !
    ! ! USES:
    use CNSharedParamsMod                 , only : CNParamsReadShared
    use CNGapMortalityMod                 , only : readCNGapMortParams                    => readParams
    use CNMRespMod                        , only : readCNMRespParams                      => readParams
    use CNFUNMod                          , only : readCNFUNParams                        => readParams
    use CNPhenologyMod                    , only : readCNPhenolParams                     => readParams
    use SoilBiogeochemCompetitionMod      , only : readSoilBiogeochemCompetitionParams    => readParams
    use SoilBiogeochemNLeachingMod        , only : readSoilBiogeochemNLeachingParams      => readParams
    use SoilBiogeochemNitrifDenitrifMod   , only : readSoilBiogeochemNitrifDenitrifParams => readParams
    use SoilBiogeochemLittVertTranspMod   , only : readSoilBiogeochemLittVertTranspParams => readParams
    use SoilBiogeochemPotentialMod        , only : readSoilBiogeochemPotentialParams      => readParams
    use SoilBiogeochemDecompMod           , only : readSoilBiogeochemDecompParams         => readParams
    use SoilBiogeochemDecompCascadeBGCMod , only : readSoilBiogeochemDecompBgcParams      => readParams
    use SoilBiogeochemDecompCascadeCNMod  , only : readSoilBiogeochemDecompCnParams       => readParams
    use ch4Mod                            , only : readCH4Params                          => readParams
    use NutrientCompetitionMethodMod      , only : nutrient_competition_method_type
    use clm_varctl,                         only : NLFilename_in
    use PhotosynthesisMod                 , only : photosyns_type
    use EDSharedParamsMod                 , only : EDParamsReadShared
    !
    ! !ARGUMENTS:
    type(photosyns_type)                   , intent(in) :: photosyns_inst
    class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    character(len=32)  :: subname = 'readParameters'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: reading CLM '//' parameters '
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    !
    ! Above ground biogeochemistry...
    !
    if (use_cn) then
       call nutrient_competition_method%readParams(ncid)
       call readCNGapMortParams(ncid)
       call readCNMRespParams(ncid)
       call readCNFUNParams(ncid)
       call readCNPhenolParams(ncid)
    end if

    !
    ! Soil biogeochemistry...
    !
    if (use_cn .or. use_ed) then
       call readSoilBiogeochemCompetitionParams(ncid)
       call readSoilBiogeochemDecompBgcParams(ncid)
       call readSoilBiogeochemDecompCnParams(ncid)
       call readSoilBiogeochemDecompParams(ncid)
       call readSoilBiogeochemLittVertTranspParams(ncid)
       call readSoilBiogeochemNitrifDenitrifParams(ncid)
       call readSoilBiogeochemNLeachingParams(ncid)
       call readSoilBiogeochemPotentialParams(ncid)
       call CNParamsReadShared(ncid, NLFilename_in)  ! this is called CN params but really is for the soil biogeochem parameters

       ! FIXME(bja, 2017-01) ED shared params must be read from the
       ! host file, not the fates file to be consistent with the host.
       call EDParamsReadShared(ncid)

       call readCH4Params (ncid)
    end if

    !
    ! Biogeophysics
    !
    call photosyns_inst%ReadParams( ncid )


    !
    call ncd_pio_closefile(ncid)

    call readFatesParameters()

  end subroutine readParameters

  !-----------------------------------------------------------------------
  subroutine readFatesParameters()

    use clm_varctl, only : fates_paramfile
    use shr_kind_mod, only: r8 => shr_kind_r8
    use paramUtilMod, only : readNcdio

    use EDParamsMod, only : FatesRegisterParams, FatesReceiveParams
    use SFParamsMod, only : SFParamsRead, SpitFireRegisterParams, SpitFireReceiveParams
    use FatesParametersInterface, only : fates_parameters_type
    use FatesParametersInterface, only : param_string_length
    use FatesParametersInterface, only : dimension_shape_scalar, dimension_shape_1d, dimension_shape_2d
    use FatesParametersInterface, only : dimension_name_pft

    implicit none

    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    character(len=32)  :: subname = 'readFatesParameters'
    class(fates_parameters_type), allocatable :: fates_params
    integer :: i, num_params, dimension_shape
    real(r8), allocatable :: data(:, :)
    character(len=param_string_length) :: name

    if (use_ed) then
       if (masterproc) then
          write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: CLM reading ED/FATES '//' parameters '
       end if

       call getfil (fates_paramfile, locfn, 0)
       call ncd_pio_openfile (ncid, trim(locfn), 0)
       call ncd_inqdid(ncid, dimension_name_pft, dimid)
       call ncd_inqdlen(ncid, dimid, npft)

       ! read using the old infrastrructure
       call SFParamsRead(ncid)

       ! read parameters with new fates parammeter infrastructure
       allocate(fates_params)
       allocate(data(npft, npft)) ! FIXME(bja, 2017-01) correct? maxpft? needs to be max of all dimensions expected to receive!
       call fates_params%Init()
       call FatesRegisterParams(fates_params)
       call SpitFireRegisterParams(fates_params)
       num_params = fates_params%num_params()
       do i = 1, num_params
          call fates_params%GetMetaData(i, name, dimension_shape)
          select case(dimension_shape)
          case(dimension_shape_scalar)
             call readNcdio(ncid, name, subname, data(1, 1))
             call fates_params%SetData(i, data(1, 1))
          case(dimension_shape_1d)
             call readNcdio(ncid, name, subname, data(:, 1))
             call fates_params%SetData(i, data(:, 1))
          case(dimension_shape_2d)
             call readNcdio(ncid, name, subname, data(:, :))
             call fates_params%SetData(i, data(:, :))
          case default
             ! error, unsupported number of dimensions
          end select
       end do
       call FatesReceiveParams(fates_params)
       call SpitFireReceiveParams(fates_params)
       deallocate(data)
       deallocate(fates_params)
       call ncd_pio_closefile(ncid)
    end if

  end subroutine readFatesParameters

end module readParamsMod
