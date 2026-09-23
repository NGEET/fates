module EDAccumulateFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period. 
  ! The fluxes are stored per cohort, and the gpp_tstep (etc) fluxes are calculated in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times. 
  ! Rosie Fisher. March 2014. 
  !
  ! !USES:
  use FatesGlobals, only      : endrun => fates_endrun 
  use FatesGlobals, only      : fates_log
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use FatesConstantsMod , only : r8 => fates_r8
  use FatesConstantsMod , only : nocomp_bareground

  implicit none
  private
  !
  public :: AccumulateFluxes_ED
  public :: FatesAccumulatePatchFluxes
  
  logical :: debug = .false.  ! for debugging this module

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------------

  subroutine AccumulateFluxes_ED(nsites, sites, bc_in, bc_out, dt_time)

    !
    ! !DESCRIPTION:
    ! see above
    !
    ! !USES:

    use EDTypesMod        , only : ed_site_type
    use FatesPatchMod,      only : fates_patch_type
    use FatesInterfaceTypesMod , only : bc_in_type,bc_out_type

    !
    ! !ARGUMENTS    
    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type),   intent(in)            :: bc_in(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)
    real(r8),           intent(in)            :: dt_time  ! timestep interval
    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type) , pointer  :: cpatch ! current patch
    integer :: s  ! ed site
    integer :: ifp ! index fates patch
    !----------------------------------------------------------------------

    do s = 1, nsites

       ! Note: Do not attempt to accumulate or log any
       ! heterotrophic respiration fluxes from the HLM here
       ! It is likely this has not been calculated yet (ELM/CLM)
       
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))

          ifp = cpatch%patchno
          
          if(cpatch%nocomp_pft_label.ne.nocomp_bareground)then

             if( bc_in(s)%filter_photo_pa(ifp) == 3 ) then

                call FatesAccumulatePatchFluxes(ifp, sites(s), dt_time)

             end if
          end if

          cpatch => cpatch%younger
       end do
    end do
    return

  end subroutine AccumulateFluxes_ED

  ! ============================================================================
  
  subroutine FatesAccumulatePatchFluxes(ifp, site, dt_time)

    !
    ! !DESCRIPTION:
    ! see above
    !
    ! !USES:

    use EDTypesMod        , only : ed_site_type
    use FatesPatchMod,      only : fates_patch_type
    use FatesCohortMod,     only : fates_cohort_type

    !
    ! !ARGUMENTS    
    integer,            intent(in)            :: ifp
    type(ed_site_type), intent(inout), target :: site
    real(r8),           intent(in)            :: dt_time  ! timestep interval
    !
    ! !LOCAL VARIABLES:
    integer :: ico
    type(fates_cohort_type), pointer  :: ccohort ! current cohort
    type(fates_patch_type) , pointer  :: cpatch ! current patch
    integer :: iv !leaf layer

    !----------------------------------------------------------------------

    cpatch => site%pa_vec(ifp)%p

    ico=0
    ccohort => cpatch%tallest
    do while(associated(ccohort))
       ico=ico+1
       ! Accumulate fluxes from hourly to daily values. 
       ! _tstep fluxes are KgC/indiv/timestep _acc are KgC/indiv/day

       ! g_sb_laweight has an order mis-match, which forces us
       ! to keep it's cohort variable for now (since we need this
       ! variable on restart)
       ccohort%g_sb_laweight = cpatch%coarrays%g_sb_laweight(ico)
       ccohort%gpp_acc  = ccohort%gpp_acc  + cpatch%coarrays%gpp_tstep(ico)
       ccohort%resp_m_acc = ccohort%resp_m_acc + cpatch%coarrays%resp_m_tstep(ico)
       
       ccohort%sym_nfix_daily = ccohort%sym_nfix_daily + cpatch%coarrays%sym_nfix_tstep(ico)
       
       ! weighted mean of D13C by gpp
       if((ccohort%gpp_acc + cpatch%coarrays%gpp_tstep(ico)) .eq. 0.0_r8) then
          ccohort%c13disc_acc = 0.0_r8
       else
          ccohort%c13disc_acc  = ((ccohort%c13disc_acc * ccohort%gpp_acc) + &
               (cpatch%coarrays%c13disc_clm(ico)* cpatch%coarrays%gpp_tstep(ico))) / &
               (ccohort%gpp_acc + cpatch%coarrays%gpp_tstep(ico))
       endif
       
       do iv=1,ccohort%nv
          if(ccohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year. 
             ccohort%year_net_uptake(iv) = 0._r8
          end if
          ccohort%year_net_uptake(iv) = ccohort%year_net_uptake(iv) + cpatch%coarrays%ts_net_uptake(iv,ico)
       enddo
       
       ccohort => ccohort%shorter
    enddo
    return
  end subroutine FatesAccumulatePatchFluxes


  
end module EDAccumulateFluxesMod

