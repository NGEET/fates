
module EDLoggingMortalityMod

   ! ====================================================================================
   !  Purpose: 1. create logging mortalities: 
   !           (a)logging mortality (cohort level)
   !           (b)collateral mortality (cohort level)
   !           (c)infrastructure mortality (cohort level)
   !           2. move the logged trunk fluxes from live into product pool 
   !           3. move logging-associated mortality fluxes from live to CWD
   !           4. keep carbon balance (in ed_total_balance_check)
   !
   !  Yi Xu
   !  Date: 2017
   ! ====================================================================================

   use FatesConstantsMod , only : r8 => fates_r8
   use EDTypesMod        , only : ed_cohort_type
   use EDTypesMod        , only : ed_patch_type
   use EDTypesMod        , only : ncwd
   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_resources_management_type
   use EDTypesMod        , only : dtype_ilog
   use EDTypesMod        , only : dtype_ifall
   use EDTypesMod        , only : dtype_ifire
   use EDPftvarcon       , only : EDPftvarcon_inst
   use EDParamsMod       , only : logging_event_code
   use EDParamsMod       , only : logging_dbhmin
   use EDParamsMod       , only : logging_collateral_frac 
   use EDParamsMod       , only : logging_direct_frac
   use EDParamsMod       , only : logging_mechanical_frac 
   use FatesInterfaceMod , only : hlm_current_year
   use FatesInterfaceMod , only : hlm_current_month
   use FatesInterfaceMod , only : hlm_current_day
   use FatesInterfaceMod , only : hlm_model_day
   use FatesInterfaceMod , only : hlm_day_of_year
   use FatesConstantsMod , only : itrue,ifalse
   use FatesGlobals      , only : endrun => fates_endrun 
   use FatesGlobals      , only : fates_log
   use shr_log_mod       , only : errMsg => shr_log_errMsg

   implicit none
   private

   logical, protected :: logging_time   ! If true, logging should be 
                                        ! performed during the current time-step

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   public :: LoggingMortality_frac
   public :: logging_litter_fluxes
   public :: logging_time
   public :: IsItLoggingTime

contains

   subroutine IsItLoggingTime(is_master,currentSite)

      ! -------------------------------------------------------------------------------
      ! This subroutine determines if the current dynamics step should enact
      ! the logging module.
      ! This is done by comparing the current model time to the logging event
      ! ids.  If there is a match, it is logging time.
      ! -------------------------------------------------------------------------------
     
      integer, intent(in) :: is_master
      type(ed_site_type), intent(inout), target :: currentSite     ! site structure

      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)
      integer :: log_date  ! Day of month for logging exctracted from event code
      integer :: log_month ! Month of year for logging extraced from event code
      integer :: log_year  ! Year for logging extracted from event code
      character(len=64) :: fmt = '(a,i2.2,a,i2.2,a,i4.4)'

      logging_time = .false.
      icode = int(logging_event_code)

!      if(hlm_use_logging.eq.ifalse) return     ! Don't turn on until fates-clm adds
                                                ! this to the interface (RGK 08-2017)

      if(icode .eq. 1) then
         ! Logging is turned off
         logging_time = .false.

      else if(icode .eq. 2) then
         ! Logging event on the first step
         if( hlm_model_day.eq.1 ) then
            logging_time = .true.
         end if

      else if(icode .eq. 3) then
         ! Logging event every day
         logging_time = .true.

      else if(icode .eq. 4) then
         ! logging event once a month
         if(hlm_current_day.eq.1  ) then
            logging_time = .true.
         end if

      else if(icode < 0 .and. icode > -366) then
         ! Logging event every year on specific day of year
         if(hlm_day_of_year .eq. icode  ) then
            logging_time = .true.
         end if

      else if(icode > 10000 ) then
         ! Specific Event: YYYYMMDD
         log_date  = icode - int(100* floor(real(icode)/100))
         log_year  = floor(real(icode)/10000)
         log_month = floor(real(icode)/100) - log_year*100

         if( hlm_current_day.eq.log_date    .and. &
               hlm_current_month.eq.log_month .and. &
               hlm_current_year.eq.log_year ) then
            logging_time = .true.
         end if
      else 
         ! Bad logging event flag
         write(fates_log(),*) 'An invalid logging code was specified in fates_params'
         write(fates_log(),*) 'Check EDLoggingMortalityMod.F90:IsItLoggingTime()'
         write(fates_log(),*) 'for a breakdown of the valid codes and change'
         write(fates_log(),*) 'fates_logging_event_code in the file accordingly.'
         write(fates_log(),*) 'exiting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      ! Initialize some site level diagnostics that are calculated for each event
      currentSite%resources_management%delta_litter_stock  = 0.0_r8
      currentSite%resources_management%delta_biomass_stock = 0.0_r8
      currentSite%resources_management%delta_individual    = 0.0_r8


      if(logging_time .and. (is_master.eq.itrue) ) then
         write(fates_log(),fmt) 'Logging Event Enacted on date: ', &
               hlm_current_month,'-',hlm_current_day,'-',hlm_current_year
      end if
      return
   end subroutine IsItLoggingTime

   ! ======================================================================================

   subroutine LoggingMortality_frac( site_in, pft_i, dbh, &
         lmort_logging,lmort_collateral,lmort_infra )

      ! Arguments
      type(ed_site_type), intent(inout), target :: site_in          ! site structure
      integer,  intent(in)                      :: pft_i            ! pft index 
      real(r8), intent(in)                      :: dbh              ! diameter at breast height (cm)
      real(r8), intent(out)                     :: lmort_logging    ! logging mortality_rates, share same size class
      real(r8), intent(out)                     :: lmort_collateral ! logging impact mortality_rates, share same size class
      real(r8), intent(out)                     :: lmort_infra      ! infrastructure mortality_rates, share same size class

      ! Locals

      ! Parameters
      real(r8), parameter                       :: adjustment = 1.0 ! adjustment for mortality rates

      if (logging_time) then 
         if(EDPftvarcon_inst%woody(pft_i) == 1)then ! only set logging rates for trees
            ! Log trees whose DBH > 50cm 
            ! Pass logging rates to cohort level 

            if (dbh >= logging_dbhmin ) then
               lmort_logging = logging_direct_frac * adjustment
            else
               lmort_logging = 0.0_r8   ! Added RGK
            end if
            lmort_collateral = logging_collateral_frac * adjustment
            lmort_infra      = logging_mechanical_frac * adjustment
            !damage rates for size class < & > threshold_size need to be specified seperately
         else
            lmort_logging    = 0.0_r8
            lmort_collateral = 0.0_r8
            lmort_infra      = 0.0_r8
         end if
      else 
         lmort_logging    = 0.0_r8
         lmort_collateral = 0.0_r8
         lmort_infra      = 0.0_r8
      end if

   end subroutine LoggingMortality_frac

   ! ============================================================================

   subroutine logging_litter_fluxes(cp_Site, cp_target, new_patch_target, patch_site_areadis)


      !  DESCRIPTION:
      !  Carbon going from ongoing mortality into CWD pools. 
      !  This module includes all the mortalities except fire mortality and logging mortalities
      !  Purpose: 
      !	  1) move logging-associated carbon to CWD and litter pool
      !   2) move the logging trunk from live into product pool 
      !   3) carbon balance check 
      !  E.g,:
      !  Remove trunk of logged trees from litter/CWD
      !  Add other parts of logged trees and all parts of collaterally and mechanically 
      !  damaged trees into CWD/litter  


      !USES:
      use SFParamsMod,  only : SF_val_cwd_frac
      use EDtypesMod,   only : area
      use EDtypesMod,   only : ed_site_type, ed_patch_type,ed_cohort_type
      use EDGrowthFunctionsMod, only : c_area

      ! !ARGUMENTS:
      type(ed_site_type)  , intent(inout), target  :: cp_Site 
      type(ed_patch_type) , intent(inout), target  :: cp_target 
      type(ed_patch_type) , intent(inout), target  :: new_patch_target

      real(r8)            , intent(inout)          :: patch_site_areadis
      !LOCAL VARIABLES:
      real(r8)                                     :: cwd_litter_density
      real(r8)                                     :: litter_area ! area over which to distribute this litter. 
      type(ed_site_type)  , pointer                :: currentSite
      type(ed_patch_type) , pointer                :: currentPatch 
      type(ed_patch_type) , pointer                :: new_patch
      type(ed_cohort_type), pointer                :: currentCohort
      real(r8)                                     :: bcroot ! amount of below ground coarse root per cohort  kgC. (goes into CWD_BG)
      real(r8)                                     :: bstem  ! amount of above ground stem biomass per cohort  kgC.(goes into CWG_AG)
      real(r8)                                     :: dead_tree_density ! # trees killed by all logging per m2
      real(r8)                                     :: damaged_dead_tree_density ! # trees killed by collateral and mechanical damages per m2
      real(r8)                                     :: logging_density ! # logged trees per m2 (logged trunk)
      real(r8)                                     :: np_mult 
                                                      !Fraction of the new patch which came from the current patch
      real(r8)                                     :: flux_out
      real(r8)                                     :: trunk_product_site !kgC/m2
      !debug
      real(r8)                                     :: delta_litter_stock !kgC/m2
      real(r8)                                     :: delta_biomass_stock !kgC/m2 
      real(r8)                                     :: delta_individual 
      integer                                      :: p,c

      currentSite  => cp_Site
      currentPatch => cp_target
      new_patch => new_patch_target

      flux_out=0.0_r8
      trunk_product_site=0.0_r8

      !debug*************************************
      delta_litter_stock=0.0_r8
      delta_biomass_stock=0.0_r8
      delta_individual=0.0_r8
      !******************************************    

      patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate 
      currentCohort => currentPatch%shortest
      do while(associated(currentCohort))       
         p = currentCohort%pft
         if(currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifire) .and. &
              currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifall) )then 
            !mortality is dominant disturbance 

            if(EDPftvarcon_inst%woody(p) == 1)then   
               !woody=1 trees only, woody=0 grass, canopy_layer=0 woody=1 small trees
               !******************************************/
               ! PART 1) Put twig, small branch and large branch of dead trees biomass from logging 
               !           and damage into above ground litter pool.
               !         Put below ground dead trees biomass 
               !         from logging and damage into below ground litter pool
               ! This happens before the plant number was updated, caused by selective logging. 
               ! but the plant number caused by other mortalities have been updated 
               ! in ED_integrate_state_variable

               !************************************/ 
               ! Number of trees that died because of the selective logging, per m2 of ground. 
               ! stem biomass per tree goes to CWD_ag
               ! Total stem biomass = dead biomass + alive biomass + store biomass (bstore) - 
               !                       root biomoass (br+bstore) - leaf biomass(bl) 
               !                    = dead biomass + alive biomass - br - bl 

               ! balive = bl + br + bsw 
               ! stem biomass per tree
               bstem  = (currentCohort%bsw + currentCohort%bdead) * &
                    EDPftvarcon_inst%allom_agb_frac(p)

               ! coarse root biomass per tree
               bcroot = (currentCohort%bsw + currentCohort%bdead) * &
                    (1.0_r8 - EDPftvarcon_inst%allom_agb_frac(p))

               ! density of dead trees per m2. Beware /AREA= 10000.0_r8 
               ! This is site-level density, not patch level 
               dead_tree_density  = (currentCohort%lmort_logging + &
                    currentCohort%lmort_collateral + currentCohort%lmort_infra )* &
                    currentCohort%n / AREA  

               ! direct logged dead tree per m2 Beware /AREA= 10000.0_r8 
               ! This is site-level density, not patch level 
               damaged_dead_tree_density = (currentCohort%lmort_collateral + &
                    currentCohort%lmort_infra )* currentCohort%n / AREA  

               new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + & 
                    dead_tree_density * (currentCohort%bl) * AREA/currentPatch%area * &
                    patch_site_areadis/new_patch%area
               new_patch%root_litter(p) = new_patch%root_litter(p) + &
                    dead_tree_density * (currentCohort%br+currentCohort%bstore)* &
                    AREA/currentPatch%area*patch_site_areadis/new_patch%area

               !KGC/m2 this is patch level density, integrated from site level 
               currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + &
                    dead_tree_density * (currentCohort%bl)*AREA/currentPatch%area
               currentPatch%root_litter(p) = currentPatch%root_litter(p) + &
                    dead_tree_density * (currentCohort%br+currentCohort%bstore) * &
                    AREA/currentPatch%area

               !debug******************************************************
               !kGC/m2 delta litter !kGC/m2
               delta_litter_stock = delta_litter_stock +dead_tree_density * &
                    (currentCohort%bl) + dead_tree_density * &
                    (currentCohort%br+currentCohort%bstore)
               delta_biomass_stock = delta_biomass_stock +dead_tree_density * &
                    (currentCohort%bdead + currentCohort%balive +currentCohort%bstore )
               delta_individual = delta_individual +dead_tree_density
               !***********************************************************

               ! track as diagnostic fluxes
               currentSite%leaf_litter_diagnostic_input_carbonflux(p) = &
                    currentSite%leaf_litter_diagnostic_input_carbonflux(p) + &
                    (currentCohort%bl) *  dead_tree_density

               currentSite%root_litter_diagnostic_input_carbonflux(p) = &
                    currentSite%root_litter_diagnostic_input_carbonflux(p) + &
                    (currentCohort%br+currentCohort%bstore) * dead_tree_density

               !above ground coarse woody debris of twigs, 
               !small branches and large branches are from logged and damaged trees 
               do c = 1,3
                  new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * &
                       SF_val_CWD_frac(c) * bstem *AREA/currentPatch%area * &
                       patch_site_areadis/new_patch%area
                  currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + &
                       dead_tree_density * SF_val_CWD_frac(c) * &
                       bstem *AREA/currentPatch%area

                  !debug***************************************************
                  !kGC/m2 delta litter !kGC/m2
                  delta_litter_stock =delta_litter_stock + dead_tree_density &
                       * SF_val_CWD_frac(c) * bstem 
                  !********************************************************

                  ! track as diagnostic fluxes
                  currentSite%CWD_AG_diagnostic_input_carbonflux(c) = &
                       currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                       SF_val_CWD_frac(c) * bstem * dead_tree_density
               enddo

               ! above ground coarse woody debris of trunk are only from damaged trees.
               !The logged trees of this part are transported off-site. 
               do c = 4,4
                  !This is patch level density, integrated from site level
                  new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + &
                       damaged_dead_tree_density * SF_val_CWD_frac(c) * &
                       bstem *AREA/currentPatch%area * patch_site_areadis/new_patch%area
                  currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + &
                       damaged_dead_tree_density * SF_val_CWD_frac(c) * &
                       bstem*AREA/currentPatch%area

                  !kGC/m2 delta litter !kGC/m2
                  delta_litter_stock =delta_litter_stock + &
                       damaged_dead_tree_density * SF_val_CWD_frac(c) * bstem

                  !*******************************************************************/
                  ! PART 2 Save logged trunk into product tool
                  ! PART 3 Pass product tool into flux out for carbon balance check 
                  !*******************************************************************/
                  logging_density = (currentCohort%lmort_logging * currentCohort%n) / AREA  

                  !kGC/m2
                  flux_out = flux_out + logging_density * SF_val_CWD_frac(c) * bstem
                  trunk_product_site = trunk_product_site + &
                       logging_density * SF_val_CWD_frac(c) * bstem  

                  ! track as diagnostic fluxes
                  currentSite%CWD_AG_diagnostic_input_carbonflux(c) = &
                       currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                       SF_val_CWD_frac(c) * bstem * damaged_dead_tree_density
               enddo

               !below ground coarse woody debris of all sizes are from logged and damaged trees
               do c = 1,ncwd
                  new_patch%cwd_bg(c) = new_patch%cwd_bg(c) + dead_tree_density * &
                       SF_val_CWD_frac(c) * bcroot *AREA/currentPatch%area * &
                       patch_site_areadis/new_patch%area

                  currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + &
                       dead_tree_density * SF_val_CWD_frac(c) * bcroot * &
                       AREA/currentPatch%area

                  !kGC/m2 delta litter !kGC/m2
                  delta_litter_stock =delta_litter_stock + &
                       dead_tree_density * SF_val_CWD_frac(c) * bcroot

                  ! track as diagnostic fluxes
                  currentSite%CWD_BG_diagnostic_input_carbonflux(c) = &
                       currentSite%CWD_BG_diagnostic_input_carbonflux(c) + &
                       SF_val_CWD_frac(c) * bcroot * dead_tree_density
               enddo
            end if
         end if
         currentCohort => currentCohort%taller      
      enddo !currentCohort

      !flux_out has been averaged to site level , not patch level
      currentSite%flux_out = currentSite%flux_out + flux_out*AREA   

      
      currentSite%resources_management%trunk_product_site  = &
           currentSite%resources_management%trunk_product_site + &
           trunk_product_site*AREA

      currentSite%resources_management%delta_litter_stock  = &
           currentSite%resources_management%delta_litter_stock + &
           delta_litter_stock* AREA

      currentSite%resources_management%delta_biomass_stock = &
           currentSite%resources_management%delta_biomass_stock + &
           delta_biomass_stock* AREA

      currentSite%resources_management%delta_individual    = &
           currentSite%resources_management%delta_individual + &
           delta_individual* AREA

      currentCohort => new_patch%shortest
      do while(associated(currentCohort))
         currentCohort%c_area = c_area(currentCohort)
         currentCohort => currentCohort%taller
      enddo

   end subroutine logging_litter_fluxes

end module EDLoggingMortalityMod
