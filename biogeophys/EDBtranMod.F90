module EDBtranMod
   
   !-------------------------------------------------------------------------------------
   ! Description:
   ! 
   ! ------------------------------------------------------------------------------------
   
   use EDPftvarcon       , only : EDPftvarcon_inst
   use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm 
   use FatesConstantsMod , only : itrue,ifalse,nearzero
   use EDTypesMod        , only : ed_site_type,       &
                                  ed_patch_type,      &
                                  ed_cohort_type,     &
                                  maxpft
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use FatesInterfaceTypesMod , only : bc_in_type, &
                                  bc_out_type, &
                                  numpft
   use FatesInterfaceTypesMod , only : hlm_use_planthydro
   use FatesGlobals      , only : fates_log
   use FatesAllometryMod , only : set_root_fraction

   !
   implicit none
   private
   
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
       if ( tempk .gt. tfrz-2._r8) then
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
      type(ed_patch_type),pointer             :: cpatch ! Current Patch Pointer
      type(ed_cohort_type),pointer            :: ccohort ! Current cohort pointer
      integer  :: s                 ! site
      integer  :: j                 ! soil layer
      integer  :: ifp               ! patch vector index for the site
      integer  :: ft                ! plant functional type index
      real(r8) :: smp_node          ! matrix potential
      real(r8) :: rresis            ! suction limitation to transpiration independent
                                    ! of root density
      real(r8) :: pftgs(maxpft)     ! pft weighted stomatal conductance m/s
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

           ifp = 0
           cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))                 
              ifp=ifp+1
              
              ! THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR COHORTS (RGK)
              
              do ft = 1,numpft
                  
                  call set_root_fraction(sites(s)%rootfrac_scr, ft, sites(s)%zi_soil ) 

                 cpatch%btran_ft(ft) = 0.0_r8
                 do j = 1,bc_in(s)%nlevsoil
                    
                    ! Calculations are only relevant where liquid water exists
                    ! see clm_fates%wrap_btran for calculation with CLM/ALM
                    
                    if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then
                       
                       smp_node = max(smpsc(ft), bc_in(s)%smp_sl(j))
                       
                       rresis  = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))*               &
                            (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)
                       
                       root_resis(ft,j) = sites(s)%rootfrac_scr(j)*rresis
                       
                       ! root water uptake is not linearly proportional to root density,
                       ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)

                       cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + root_resis(ft,j)
                       
                    else
                       root_resis(ft,j) = 0._r8
                    end if
                    
                 end do !j
                 
                 ! Normalize root resistances to get layer contribution to ET
                 do j = 1,bc_in(s)%nlevsoil  
                    if (cpatch%btran_ft(ft)  >  nearzero) then
                        root_resis(ft,j) = root_resis(ft,j)/cpatch%btran_ft(ft)
                    else
                        root_resis(ft,j) = 0._r8
                    end if
                 end do
                 
              end do !PFT
              
              ! PFT-averaged point level root fraction for extraction purposese.
              ! The cohort's conductance g_sb_laweighted, contains a weighting factor
              ! based on the cohort's leaf area. units: [m/s] * [m2]
              
              pftgs(1:maxpft) = 0._r8
              ccohort => cpatch%tallest
              do while(associated(ccohort))
                 pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%g_sb_laweight
                 ccohort => ccohort%shorter
              enddo
              
              ! Process the boundary output, this is necessary for calculating the soil-moisture
              ! sink term across the different layers in driver/host.  Photosynthesis will
              ! pass the host a total transpiration for the patch.  This needs rootr to be
              ! distributed over the soil layers.
              sum_pftgs = sum(pftgs(1:numpft))

              do j = 1, bc_in(s)%nlevsoil
                 bc_out(s)%rootr_pasl(ifp,j) = 0._r8
                 do ft = 1,numpft
                    if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                       !bit-retart test as a result? FIX(RF,032414)  
                       bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                             root_resis(ft,j) * pftgs(ft)/sum_pftgs
                    else
                       bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                             root_resis(ft,j) * 1._r8/real(numpft,r8)
                    end if
                 enddo
              enddo
              
              ! Calculate the BTRAN that is passed back to the HLM
              ! used only for diagnostics. If plant hydraulics is turned off
              ! we are using the patchxpft level btran calculation
              
              if(hlm_use_planthydro.eq.ifalse) then
                 !weight patch level output BTRAN for the
                 bc_out(s)%btran_pa(ifp) = 0.0_r8
                 do ft = 1,numpft
                    if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                       !bit-retart test as a result? FIX(RF,032414)   
                       bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft)  * pftgs(ft)/sum_pftgs
                    else
                       bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft) * 1./numpft
                    end if
                 enddo
              end if

              temprootr = sum(bc_out(s)%rootr_pasl(ifp,1:bc_in(s)%nlevsoil))

              if(abs(1.0_r8-temprootr) > 1.0e-10_r8 .and. temprootr > 1.0e-10_r8)then
                 write(fates_log(),*) 'error with rootr in canopy fluxes',temprootr,sum_pftgs
                 do j = 1,bc_in(s)%nlevsoil
                    bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j)/temprootr
                 enddo
              end if
              
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
