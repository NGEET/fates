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

    use EDTypesMod        , only : ed_site_type, AREA
    use FatesPatchMod,      only : fates_patch_type
    use FatesCohortMod,     only : fates_cohort_type
    use FatesInterfaceTypesMod , only : bc_in_type,bc_out_type
    use EDtypesMod               , only : AREA_INV
    use clm_time_manager  , only : get_curr_days_per_year
    use FatesConstantsMod        , only : sec_per_day
    use FatesConstantsMod        , only : days_per_year    
    use FatesConstantsMod        , only : g_per_kg
    
    !
    ! !ARGUMENTS    
    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type),   intent(in)            :: bc_in(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)
    real(r8),           intent(in)            :: dt_time  ! timestep interval
    !
    ! !LOCAL VARIABLES:
    type(fates_cohort_type), pointer  :: ccohort ! current cohort
    type(fates_patch_type) , pointer  :: cpatch ! current patch
    integer :: iv !leaf layer
    integer :: c  ! clm/alm column
    integer :: s  ! ed site
    integer :: ifp ! index fates patch
    real    :: ind_per_m2
    !----------------------------------------------------------------------
    
    do s = 1, nsites

       ! Note: Do not attempt to accumulate or log any
       ! heterotrophic respiration fluxes from the HLM here
       ! It is likely this has not been calculated yet (ELM/CLM)
       
       cpatch => sites(s)%oldest_patch
       bc_out(s)%npp_acc_site = 0._r8
       bc_out(s)%npp_site = 0._r8
       
       do while (associated(cpatch))

          ifp = cpatch%patchno
          
          if(cpatch%nocomp_pft_label.ne.nocomp_bareground)then

             if( bc_in(s)%filter_photo_pa(ifp) == 3 ) then
                ccohort => cpatch%shortest
                do while(associated(ccohort))
                   ind_per_m2 = ccohort%n * AREA_INV

                   ! Accumulate fluxes from hourly to daily values. 
                   ! _tstep fluxes are KgC/indiv/timestep _acc are KgC/indiv/day
                   ccohort%gpp_acc  = ccohort%gpp_acc  + ccohort%gpp_tstep 
                   ccohort%resp_m_acc = ccohort%resp_m_acc + ccohort%resp_m_tstep

                   ! Make npp_acc variable for the site level to add to the NBP balance check
                   ! Convert from kgC/ind to gC/m2
                   if(.not.ccohort%isnew)then
                      bc_out(s)%npp_acc_site = bc_out(s)%npp_acc_site + &
                     (ccohort%gpp_acc - ccohort%resp_m_acc) &
                      * ind_per_m2 * g_per_kg &
                      -(ccohort%resp_g_acc_hold+ccohort%resp_excess_hold) * & 
                      ind_per_m2 * g_per_kg * dt_time/(days_per_year*sec_per_day)
                      ! gresp is converted from kgC/indiv/year to gC/m2/timestep. 
                   endif 

               ! Net Ecosystem Production [kgC/m2/s]. Use yesterday's growth respiration
               ! This is taken from the NEP history variable calculation.
               ! first add GPP-Rm and convert units from kgC/indiv/timestep to gC/m2/s
               ! then smooth out yesterdays's calculated growth respiration and
               ! convert units from kgC/indiv/year  to gC/m2/s. 
                  if(.not.ccohort%isnew)then
                    bc_out(s)%npp_site = bc_out(s)%npp_site + (ccohort%gpp_tstep-ccohort%resp_m_tstep) &
                        * ind_per_m2 * g_per_kg / dt_time - &
                        (ccohort%resp_g_acc_hold+ccohort%resp_excess_hold) * &
                        ind_per_m2 * g_per_kg / (days_per_year*sec_per_day)
                  endif
                   ccohort%sym_nfix_daily = ccohort%sym_nfix_daily + ccohort%sym_nfix_tstep
                   
                   ! weighted mean of D13C by gpp
                   if((ccohort%gpp_acc + ccohort%gpp_tstep) .eq. 0.0_r8) then
                      ccohort%c13disc_acc = 0.0_r8
                   else
                      ccohort%c13disc_acc  = ((ccohort%c13disc_acc * ccohort%gpp_acc) + &
                           (ccohort%c13disc_clm * ccohort%gpp_tstep)) / &
                           (ccohort%gpp_acc + ccohort%gpp_tstep)
                   endif

                   do iv=1,ccohort%nv
                      if(ccohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year. 
                         ccohort%year_net_uptake(iv) = 0._r8
                      end if
                      ccohort%year_net_uptake(iv) = ccohort%year_net_uptake(iv) + ccohort%ts_net_uptake(iv)
                   enddo

                   ccohort => ccohort%taller
                enddo ! while(associated(ccohort))
             end if
          end if ! not bare ground

          cpatch => cpatch%younger
       end do  ! while(associated(cpatch))

   end do
    return

  end subroutine AccumulateFluxes_ED

end module EDAccumulateFluxesMod

