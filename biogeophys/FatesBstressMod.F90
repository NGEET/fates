module FATESBstressMod
   
   !-------------------------------------------------------------------------------------
   ! Description: calculate the stress impact on transpiration from salinity and sulphide in soils
   !              note that water stress is calculated in EDBtranMod or HYDRO
   ! ------------------------------------------------------------------------------------
   use EDPftvarcon       , only : EDPftvarcon_inst
   !
   implicit none
   private
   
   public :: btran_stress_fates

contains
 ! =====================================================================================

  subroutine btran_stress_fates( nsites, sites, bc_in, bc_out)

    use FatesPlantHydraulicsMod, only : BTranForHLMDiagnosticsFromCohortHydr

      
      ! ---------------------------------------------------------------------------------
      ! Calculate the transpiration wetness function (BTRAN) and the root uptake
      ! distribution (ROOTR).
      ! Boundary conditions in: bc_in(s)%eff_porosity_sl(j)    unfrozen porosity
      !                         bc_in(s)%watsat_sl(j)          porosity
      !                         bc_in(s)%active_uptake_sl(j)   frozen/not frozen
      !                         bc_in(s)%smp_sl(j)             suction
      ! Boundary conditions out: bc_out(s)%rootr_pasl          root uptake distribution
      !                          bc_out(s)%btran_stress_pa            wetness factor
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
      !------------------------------------------------------------------------------
      
      associate(                                 &
            smpsc     => EDPftvarcon_inst%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            smpso     => EDPftvarcon_inst%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            )
        
        do s = 1,nsites

           bc_out(s)%rootr_pasl(:,:) = 0._r8

           ifp = 0
           cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))                 
              ifp=ifp+1
              
              ! THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR COHORTS (RGK)
              
              do ft = 1,numpft
                 cpatch%btran_ft(ft) = 0.0_r8
                 do j = 1,bc_in(s)%nlevsoil
                    
                    ! Calculations are only relevant where liquid water exists
                    ! see clm_fates%wrap_btran for calculation with CLM/ALM
                    
                    if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then
                       
                       smp_node = max(smpsc(ft), bc_in(s)%smp_sl(j))
                       
                       rresis  = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))*               &
                            (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)
                       
                       cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)*rresis
                       
                       ! root water uptake is not linearly proportional to root density,
                       ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                       ! cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)**0.3*rresis_ft(ft,j)/ &
                       ! sum(cpatch%rootfr_ft(ft,1:nlevsoil)**0.3)
                       cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + cpatch%rootr_ft(ft,j)
                       
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
                    end if
                    
                 end do !j
                 
                 ! Normalize root resistances to get layer contribution to ET
                 do j = 1,bc_in(s)%nlevsoil  
                    if (cpatch%btran_ft(ft)  >  0.0_r8) then
                       cpatch%rootr_ft(ft,j) = cpatch%rootr_ft(ft,j)/cpatch%btran_ft(ft)
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
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
                            cpatch%rootr_ft(ft,j) * pftgs(ft)/sum_pftgs
                    else
                       bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                            cpatch%rootr_ft(ft,j) * 1._r8/dble(numpft)
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
        
        end do
           
        if(hlm_use_planthydro.eq.itrue) then
           call BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)
        end if
        
      end associate
      
    end subroutine btran_stress_fates      


end module FATESBsalMod

