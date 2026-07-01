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

  !
  implicit none
  private


  logical, parameter :: debug = .true.
  character(len=*), parameter :: sourcefile = __FILE__
  
  public :: btran_ed
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

  subroutine btran_ed( nsites, sites, bc_in, bc_out)

    use FatesPlantHydraulicsMod, only : BTranForHLMDiagnosticsFromCohortHydr


    ! ---------------------------------------------------------------------------------
    ! Calculate the transpiration wetness function (BTRAN) and the root uptake
    ! distribution (ROOTR).
    ! Boundary conditions in: bc_in(s)%eff_porosity_sl(j)    unfrozen porosity
    !                         bc_in(s)%watsat_sl(j)          porosity
    !                         bc_in(s)%active_uptake_sl(j)   frozen/not frozen
    !                         bc_in(s)%smp_sl(j)             suction
    ! Boundary conditions out: bc_out(s)%rootr_pasl          root uptake distribution
    !                          bc_out(s)%btran_pa            wetness factor
    ! ---------------------------------------------------------------------------------

    ! Arguments

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
    real(r8) :: rresis            ! suction limitation to transpiration independent
    ! of root density
    real(r8) :: pftgs(1:numpft)     ! pft weighted stomatal conductance m/s
    logical  :: valid_pft(1:numpft) ! pft mask
    real(r8) :: num_valid_pfts   ! 
    real(r8) :: temprootr
    real(r8) :: sum_pftgs         ! sum of weighted conductances (for normalization)
    real(r8), allocatable :: root_resis(:,:)  ! Root resistance in each pft x layer
    !------------------------------------------------------------------------------

    associate(                                 &
         smpsc     => EDPftvarcon_inst%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         smpso     => EDPftvarcon_inst%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
         )

    do s = 1,nsites

       allocate(root_resis(numpft,bc_in(s)%nlevsoil))

       bc_out(s)%rootr_pasl(:,:) = 0._r8

       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))

          ifp = cpatch%patchno
          
          if_bare: if(cpatch%nocomp_pft_label.ne.nocomp_bareground .and. (cpatch%num_cohorts > 0))then ! only for veg patches

             valid_pft(1:numpft) = .false.
             pftgs(1:numpft)=0.0_r8
             sum_pftgs = 0.0_r8
             ccohort => cpatch%tallest
             do while(associated(ccohort))
                if (ccohort%g_sb_laweight > rsnbl_math_prec) then
                   sum_pftgs = sum_pftgs + ccohort%g_sb_laweight
                   pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%g_sb_laweight
                   valid_pft(ccohort%pft) = .true.
                endif
                ccohort => ccohort%shorter
             enddo
             num_valid_pfts = 0.0_r8
             do ft = 1,numpft
                if(valid_pft(ft)) num_valid_pfts = num_valid_pfts +1._r8
             enddo
                  
             do ft = 1,numpft
                if (valid_pft(ft)) then
                   if (sum_pftgs < num_valid_pfts * rsnbl_math_prec) then
                         pftgs(ft) = 1._r8/num_valid_pfts ! condition above already means 1 valid pft
                   else
                         pftgs(ft) = pftgs(ft)/sum_pftgs
                   endif
                else
                  pftgs(ft) = 0._r8
                endif
             end do
               
             ! THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR COHORTS (RGK)
             cpatch%btran_ft(:) = 0.0_r8
             root_resis(:,:) = 0.0_r8
             do ft = 1,numpft
                if (valid_pft(ft)) then
                   call set_root_fraction(sites(s)%rootfrac_scr, ft, sites(s)%zi_soil, &
                         bc_in(s)%max_rooting_depth_index_col ) 

                   do j = 1,bc_in(s)%nlevsoil

                      ! Calculations are only relevant where liquid water exists
                      ! see clm_fates%wrap_btran for calculation with CLM/ALM
                      if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then

                         smp_node = max(smpsc(ft), bc_in(s)%smp_sl(j))
      
                         rresis  = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))*               &
                               (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)

                         root_resis(ft,j) = root_resis(ft,j) + rresis * sites(s)%rootfrac_scr(j)

                         ! root water uptake is not linearly proportional to root density,
                         ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)

                         cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + root_resis(ft,j)

                      else
                         root_resis(ft,j) = 0._r8
                      end if
                   end do !j
                else
                   root_resis(ft,1:bc_in(s)%nlevsoil) = 0.0_r8
                   cpatch%btran_ft(ft) = 0.0_r8
                endif ! valid_pfts
             end do ! PFT

               ! remove this check when merging to noresm
               if (debug) then
                  do ft=1,numpft
                  if (sum(root_resis(ft,1:bc_in(s)%nlevsoil)) .ne. cpatch%btran_ft(ft)) then
                     write(fates_log(),*) 'MVD: btran_ft not equal rresist ',cpatch%nocomp_pft_label,ft
                     write(fates_log(),*)' btran sum', cpatch%btran_ft(ft),sum(root_resis(ft,1:bc_in(s)%nlevsoil))
                     write(fates_log(),*)' rootr ',root_resis(ft,1:bc_in(s)%nlevsoil)
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  endif
                  enddo
               endif

             ! PFT-averaged point level root fraction for extraction purposese.
             ! The cohort's conductance g_sb_laweighted, contains a weighting factor
             ! based on the cohort's leaf area. units: [m/s] * [m2]
             do ft = 1,numpft
                ! Normalize root resistances to get layer contribution to ET
                do j = 1,bc_in(s)%nlevsoil  
                   if (cpatch%btran_ft(ft)  >  nearzero) then
                      root_resis(ft,j) = root_resis(ft,j)/cpatch%btran_ft(ft)
                   else
                      root_resis(ft,j) = 0._r8
                   end if
                end do
             enddo
             ! Process the boundary output, this is necessary for calculating the soil-moisture
             ! sink term across the different layers in driver/host.  Photosynthesis will
             ! pass the host a total transpiration for the patch.  This needs rootr to be
             ! distributed over the soil layers.
             bc_out(s)%rootr_pasl(ifp,:) = 0.0_r8
             bc_out(s)%btran_pa(ifp)  = 0.0_r8
             do j = 1, bc_in(s)%nlevsoil
                do ft = 1,numpft
                   bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                             root_resis(ft,j) * pftgs(ft)
                enddo
             enddo
             if (hlm_use_planthydro.eq.ifalse) then
                do ft = 1,numpft
                   bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft) *pftgs(ft)
                end do
             endif

             temprootr = sum(bc_out(s)%rootr_pasl(ifp,1:bc_in(s)%nlevsoil))

             if(abs(1.0_r8-temprootr) > rsnbl_math_prec .and. abs(temprootr) > rsnbl_math_prec)then
                if (debug) then 
                   write(fates_log(),*) 'error with rootr in canopy fluxes',temprootr,sum_pftgs,bc_out(s)%btran_pa(ifp)
                   ! remove this endrun later.
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                endif
                temprootr = abs(temprootr)
                do j = 1,bc_in(s)%nlevsoil
                   bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j)/temprootr
                enddo
                
             end if
          endif if_bare
          cpatch => cpatch%younger
       end do

       deallocate(root_resis)

    end do

    if(hlm_use_planthydro.eq.itrue) then
       call BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)
    end if

  end associate

end subroutine btran_ed


end module EDBtranMod
