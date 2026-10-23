module EDBtranMod

  !-------------------------------------------------------------------------------------
  ! Description:
  ! 
  ! ------------------------------------------------------------------------------------

  use EDPftvarcon       , only : EDPftvarcon_inst
  use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm 
  use FatesConstantsMod , only : itrue,ifalse,nearzero
  use FatesConstantsMod , only : nocomp_bareground
  use FatesConstantsMod , only : rsnbl_math_prec
  use EDTypesMod        , only : ed_site_type
  use FatesPatchMod,      only : fates_patch_type
  use EDParamsMod,        only : maxpft
  use EDParamsMod,        only : soil_tfrz_thresh
  use FatesCohortMod,     only : fates_cohort_type
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use FatesInterfaceTypesMod , only : bc_in_type, &
       bc_out_type, &
       numpft
  use FatesInterfaceTypesMod , only : hlm_use_planthydro
  use FatesGlobals      , only : fates_log
  use FatesAllometryMod , only : set_root_fraction
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesPlantHydraulicsMod, only : BTranForHLMDiagnosticsFromCohortHydr
  use FatesGlobals, only : FatesWarn,N2S,A2S,I2S

  !
  implicit none
  private

  character(len=512) :: warn_msg  ! Message string for warnings and logging
  logical, parameter :: debug = .false.
  
  public :: RootUptakeFractions
  public :: FatesInternalBtran
  public :: get_active_suction_layers
  public :: check_layer_water

contains 

  ! ====================================================================================

  logical function check_layer_water(h2o_liq_vol, tempk)

    implicit none
    ! Arguments
    real(r8),intent(in) :: h2o_liq_vol
    real(r8),intent(in) :: tempk

    check_layer_water = .false.

    if ( h2o_liq_vol .gt. 0._r8 ) then
       if ( tempk .gt. soil_tfrz_thresh + tfrz) then
          check_layer_water = .true.
       end if
    end if
    return
  end function check_layer_water

  ! =====================================================================================

  subroutine get_active_suction_layers(nsites, sites, bc_in, bc_out)

    ! Arguments

    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    ! !LOCAL VARIABLES:
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    !------------------------------------------------------------------------------

    do s = 1,nsites
       if (bc_in(s)%filter_btran) then
          do j = 1,bc_in(s)%nlevsoil
             bc_out(s)%active_suction_sl(j) = check_layer_water( bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j) )
          end do
       else
          bc_out(s)%active_suction_sl(:) = .false.
       end if
    end do

  end subroutine get_active_suction_layers

  ! =====================================================================================

  subroutine FatesInternalBtran( nsites, sites, bc_in, bc_out)


    ! ---------------------------------------------------------------------------------
    ! Calculate the transpiration wetness function (BTRAN) internal to FATES,
    ! currently this is specific to each patch and PFT. Although, patch%btran_ft
    ! does NOT have any patch specific information used, it should be site. However,
    ! with cohorts having different rooting depths, this ideally should be migrated
    ! to a cohort specific variable.
    !
    ! Boundary conditions in: bc_in(s)%eff_porosity_sl(j)    unfrozen porosity
    !                         bc_in(s)%watsat_sl(j)          porosity
    !                         bc_in(s)%active_uptake_sl(j)   frozen/not frozen
    !                         bc_in(s)%smp_sl(j)             suction
    ! ---------------------------------------------------------------------------------

    ! Arguments

    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type),pointer             :: cpatch ! Current Patch Pointer
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    integer  :: ifp               ! patch vector index for the site
    integer  :: ft                ! plant functional type index
    real(r8) :: smp_node          ! matrix potential
    real(r8) :: suction_lim       ! suction limitation independent
                                  ! of root density
    !------------------------------------------------------------------------------

    associate(                                 &
         smpsc     => EDPftvarcon_inst%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         smpso     => EDPftvarcon_inst%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         )

    do s = 1,nsites

       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))

          ifp = cpatch%patchno
          cpatch%btran_ft(:) = 0._r8
          if_bare: if(cpatch%nocomp_pft_label.ne.nocomp_bareground)then ! only for veg patches

             ! PFT-averaged point level root fraction for extraction purposese.
             ! The cohort's conductance g_sb_laweighted, contains a weighting factor
             ! based on the cohort's leaf area. units: [m/s] * [m2]
             
             do ft = 1,numpft

                call set_root_fraction(sites(s)%rootfrac_scr, ft, sites(s)%zi_soil, &
                     bc_in(s)%max_rooting_depth_index_col ) 

                do j = 1,bc_in(s)%nlevsoil

                   ! Calculations are only relevant where liquid water exists
                   ! see clm_fates%wrap_btran for calculation with CLM/ALM

                   if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then

                      smp_node = max(smpsc(ft), bc_in(s)%smp_sl(j))

                      ! (negative - smaller negative) / (larger negative - smaller negative)
                      
                      suction_lim  = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))* &
                           (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)
                      
                      ! root water uptake is not linearly proportional to root density,
                      ! to allow proper deep root function. Replace with equations from SPA/Newman. FIX(RF,032414)
                      
                      cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + suction_lim * sites(s)%rootfrac_scr(j)
                   end if
                end do
             end do
          end if if_bare

          cpatch => cpatch%younger
       end do
    end do
    end associate
  
  end subroutine FatesInternalBtran

  ! ========================================================================
  
  subroutine RootUptakeFractions(nsites, sites, bc_in, bc_out)

    ! Calculate uptake fractions that are used by the
    ! land model so that we can identify the correct layers
    ! to pull soil moisture from
    
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type),pointer             :: cpatch ! Current Patch Pointer
    type(fates_cohort_type),pointer            :: ccohort ! Current cohort pointer
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    integer  :: ifp               ! patch vector index for the site
    integer  :: ft                ! plant functional type index
    real(r8) :: smp_node          ! matrix potential
    real(r8) :: suction_lim       ! suction limitation to transpiration independent
    ! of root density
    real(r8) :: pftgs(numpft)     ! pft weighted stomatal conductance m/s
    real(r8) :: temprootr
    real(r8) :: sum_pftgs         ! sum of weighted conductances (for normalization)
    real(r8) :: root_relcond(numpft,numlevsoil_max)  ! relative root conductance in each pft x layer
    real(r8) :: num_present_pfts
    real(r8) :: pfts_present(numpft)

    
    associate(                                 &
         smpsc     => EDPftvarcon_inst%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         smpso     => EDPftvarcon_inst%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         )

    do s = 1,nsites

       root_relcond(:,:) = 0._r8
       bc_out(s)%rootr_pasl(:,:) = 0._r8

       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))

          ifp = cpatch%patchno
          
          if_bare: if(cpatch%nocomp_pft_label.ne.nocomp_bareground .and. &
                      cpatch%num_cohorts > 0 )then ! only for veg patches

             do ft = 1,numpft

                call set_root_fraction(sites(s)%rootfrac_scr, ft, sites(s)%zi_soil, &
                     bc_in(s)%max_rooting_depth_index_col ) 

                do j = 1,bc_in(s)%nlevsoil

                   ! Calculations are only relevant where liquid water exists
                   ! see clm_fates%wrap_btran for calculation with CLM/ALM

                   if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then
                      
                      smp_node = max(smpsc(ft), bc_in(s)%smp_sl(j))

                      ! higher conducts more water (0-1)
                      suction_lim  = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))* &
                           (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)

                      root_relcond(ft,j) = sites(s)%rootfrac_scr(j)*suction_lim
                   else
                      root_relcond(ft,j) = 0._r8
                   end if
                   
                end do !j
                
                ! Normalize relative root concuctances to get layer contribution to ET
                do j = 1,bc_in(s)%nlevsoil  
                   if (cpatch%btran_ft(ft)  >  nearzero) then
                      root_relcond(ft,j) = root_relcond(ft,j)/cpatch%btran_ft(ft)
                   else
                      root_relcond(ft,j) = 0._r8
                   end if
                end do

             end do
             
             pfts_present(1:numpft) = 0._r8
             pftgs(1:numpft)=0._r8
             sum_pftgs = 0._r8
             num_present_pfts = 0._r8
             ccohort => cpatch%tallest
             do while(associated(ccohort))
                if(cpatch%btran_ft(ccohort%pft) > nearzero)then
                   if (ccohort%g_sb_laweight > 0._r8) then
                      sum_pftgs = sum_pftgs + ccohort%g_sb_laweight
                      pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%g_sb_laweight
                   end if
                   pfts_present(ccohort%pft) = 1._r8
                end if
                ccohort => ccohort%shorter
             end do

             num_present_pfts = sum(pfts_present(1:numpft))
             
             ! Process the boundary output, this is necessary for calculating the soil-moisture
             ! sink term across the different layers in driver/host.  Photosynthesis will
             ! pass the host a total transpiration for the patch.  This needs rootr to be
             ! distributed over the soil layers.
             
             if (num_present_pfts > 0._r8) then
                ! this is a bit arbitrary threshold but with tiny weight it
                ! does not matter which pft has lower resistance.
                if (sum_pftgs > rsnbl_math_prec) then  
                   pftgs(1:numpft)=pftgs(1:numpft)/sum_pftgs
                else
                   pftgs(1:numpft)=pfts_present(1:numpft)/num_present_pfts
                endif
             else 
                pftgs(1:numpft)=0._r8
             endif

             bc_out(s)%rootr_pasl(ifp,1:bc_in(s)%nlevsoil) = 0._r8
             bc_out(s)%btran_pa(ifp) = 0._r8
             do  ft = 1,numpft
                do j = 1, bc_in(s)%nlevsoil
                   bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                        root_relcond(ft,j) * pftgs(ft)
                enddo
                ! Calculate the BTRAN that is passed back to the HLM
                ! used only for diagnostics. If plant hydraulics is turned off
                ! we are using the patchxpft level btran calculation
                if (hlm_use_planthydro .eq. ifalse) then
                   bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft)  * pftgs(ft)
                end if
             enddo
                   
             temprootr = sum(bc_out(s)%rootr_pasl(ifp,1:bc_in(s)%nlevsoil))

             if(abs(1.0_r8-temprootr) > rsnbl_math_prec .and. temprootr > rsnbl_math_prec)then
                
                warn_msg = 'Normalized root resistances dont sum to 1, s:'//trim(I2S(s))// &
                     ' ifp: '//trim(I2S(ifp))// &
                     ' temprootr: '//trim(N2S(temprootr)//' sum_pftgs: '//trim(N2S(sum_pftgs))// &
                     ' num_present_pfts: '//trim(N2S(num_present_pfts))
                call FatesWarn(warn_msg,index=5)
                
                do j = 1,bc_in(s)%nlevsoil
                   bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j)/temprootr
                enddo
                
             end if
          endif if_bare
          cpatch => cpatch%younger
       end do

    end do
  end associate
  
  if(hlm_use_planthydro.eq.itrue) then
     call BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)
  end if
  
  return
end subroutine RootUptakeFractions


end module EDBtranMod
