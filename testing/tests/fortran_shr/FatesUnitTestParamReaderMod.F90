module FatesUnitTestParamReaderMod

  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesConstantsMod,          only : fates_check_param_set
  use FatesParametersInterface,   only : pstruct
  use JSONParameterUtilsMod,      only : JSONRead
  use JSONParameterUtilsMod,      only : JSONSetInvalid
  use JSONParameterUtilsMod,      only : JSONSetLogInit
  use JSONParameterUtilsMod,      only : JSONDumpParameter
  use PRTParametersMod,           only : prt_params
  use FatesInterfaceTypesMod,     only : nleafage
  use FatesInterfaceTypesMod,     only : numpft
  use FatesParameterDerivedMod,   only : param_derived
  use FatesGlobals,               only : fates_log
  use FatesLeafBiophysParamsMod,  only : TransferParamsLeafBiophys
  use EDParamsMod,                only : TransferParamsGeneric
  use SFParamsMod,                only : TransferParamsSpitFire
  use PRTInitParamsFatesMod,      only : TransferParamsPRT
  use PRTInitParamsFatesMod,      only : PRTDerivedParams
  use EDPftvarcon,                only : TransferParamsPFT
  use FatesInterfaceTypesMod,     only : hlm_maintresp_leaf_model
  use FatesConstantsMod,          only : lmrmodel_atkin_etal_2017
  use LeafBiophysicsMod,          only : NegativeRdarkTempC
  use FatesTestLeafPhotoMod,      only : LeafNitrogenContent

  implicit none
  private
  
  ! growth temperature the Atkin et al. (2017) parameter check reports against.
  ! Not a limit the drivers enforce - a PFT whose threshold falls below this is
  ! one whose Rdark will be capped at zero over part of a normal driver run
  real(r8), parameter :: reference_tgrowth_C = 25.0_r8 ! (degrees C)
  

  public :: ReadParameters
  public :: CheckLeafRespParams

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
   
    nleafage = size(prt_params%leaf_long, dim=2)
    numpft = size(prt_params%wood_density, dim=1)
  
    ! initialize derived parameters
    call param_derived%Init(size(prt_params%wood_density, dim=1))
    
    call PRTDerivedParams() 

  end subroutine ReadParameters
  
  ! --------------------------------------------------------------------------------------

  subroutine CheckLeafRespParams()
    !
    ! DESCRIPTION:
    ! Reports the growth temperature above which each PFT's Atkin et al. (2017)
    ! reference dark respiration is capped at zero, and flags any PFT whose
    ! threshold falls below reference_tgrowth_C
    !
    ! Production runs this as part of EDPftvarcon's FatesCheckParams, which the
    ! standalone drivers do not call - it depends on host-model namelist state
    ! no driver sets
    !
    ! A no-op unless the Atkin model is selected, so the calling driver must
    ! have set hlm_maintresp_leaf_model already, and PRTDerivedParams() must
    ! have populated prt_params%organ_param_id (LeafNitrogenContent needs it)
    !

    ! LOCALS:
    integer  :: ipft          ! pft looping index
    real(r8) :: lnc_top       ! leaf N content at the canopy top [gN/m2 leaf]
    real(r8) :: neg_lmr_tempC ! growth temperature at which Rdark is capped at zero [degrees C]

    if (hlm_maintresp_leaf_model /= lmrmodel_atkin_etal_2017) return
    do ipft = 1, numpft
      lnc_top = LeafNitrogenContent(ipft)
      neg_lmr_tempC = NegativeRdarkTempC(ipft, lnc_top)
      if (neg_lmr_tempC < reference_tgrowth_C) then
        write(*,*) '  pft ', ipft, neg_lmr_tempC, ' WARNING: below reference growth temperature'
      end if
    end do
    
  ! --------------------------------------------------------------------------------------

  end subroutine CheckLeafRespParams

end module FatesUnitTestParamReaderMod
