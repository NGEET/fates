module FatesUnitTestParamReaderMod

  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesConstantsMod         , only : fates_check_param_set
  use FatesParametersInterface,   only : pstruct
  use JSONParameterUtilsMod,      only : JSONRead
  use JSONParameterUtilsMod,  	  only : JSONSetInvalid
  use JSONParameterUtilsMod,  	  only : JSONSetLogInit
  use JSONParameterUtilsMod,  	  only : JSONDumpParameter
  use PRTParametersMod,           only : prt_params
  use FatesInterfaceTypesMod,     only : nleafage
  use FatesParameterDerivedMod,   only : param_derived
  use FatesGlobals              , only : fates_log
  use FatesLeafBiophysParamsMod , only : TransferParamsLeafBiophys
  use EDParamsMod               , only : TransferParamsGeneric
  use SFParamsMod               , only : TransferParamsSpitFire
  use PRTInitParamsFatesMod     , only : TransferParamsPRT
  use EDPftvarcon               , only : TransferParamsPFT

  implicit none
  private
  
  logical, parameter :: debug=.false.

  public :: ReadParameters

contains
  
  ! --------------------------------------------------------------------------------------

  subroutine ReadParameters(param_file)
    !
    ! DESCRIPTION:
    ! Read 'fates_params' parameters from storage
    !
    ! ARGUMENTS:
    character(len=*) :: param_file

    call JSONSetInvalid(fates_check_param_set+10._r8)
    call JSONSetLogInit(fates_log())
    call JSONRead(param_file,pstruct)

    ! Transfer parameters from the json datastructure, into
    ! primitive data structures
    call TransferParamsGeneric(pstruct)
    call TransferParamsSpitFire(pstruct)
    call TransferParamsPRT(pstruct)
    call TransferParamsLeafBiophys(pstruct)
    call TransferParamsPFT(pstruct)
   
    if(debug) call pstruct%ReportAccessCounts()

    nleafage = size(prt_params%leaf_long, dim=2)
  
    ! initialize derived parameters
    call param_derived%Init(size(prt_params%wood_density, dim=1))

  end subroutine ReadParameters

end module FatesUnitTestParamReaderMod
