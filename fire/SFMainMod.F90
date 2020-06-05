  module SFMainMod

  ! ============================================================================
  ! All subroutines realted to the SPITFIRE fire routine. 
  ! Code originally developed by Allan Spessa & Rosie Fisher as part of the NERC-QUEST project.  
  ! ============================================================================

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesInterfaceTypesMod     , only : hlm_masterproc ! 1= master process, 0=not master process
  use EDTypesMod            , only : numWaterMem
  use FatesGlobals          , only : fates_log

  use FatesInterfaceTypesMod     , only : bc_in_type

  use EDPftvarcon           , only : EDPftvarcon_inst

  use EDTypesMod            , only : element_pos
  use EDtypesMod            , only : ed_site_type
  use EDtypesMod            , only : ed_patch_type
  use EDtypesMod            , only : ed_cohort_type
  use EDtypesMod            , only : AREA
  use EDtypesMod            , only : DL_SF
  use EDTypesMod            , only : TW_SF
  use EDtypesMod            , only : LB_SF
  use EDtypesMod            , only : LG_SF
  use FatesLitterMod        , only : ncwd
  use EDtypesMod            , only : NFSC
  use EDtypesMod            , only : TR_SF
  use FatesLitterMod        , only : litter_type

  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : SetState
  use FatesInterfaceTypesMod     , only : numpft

  implicit none
  private

  public :: fire_model
  public :: fire_danger_index 
  public :: charecteristics_of_fuel
  public :: rate_of_spread
  public :: ground_fuel_consumption
  public :: wind_effect
  public :: area_burnt_intensity
  public :: crown_scorching
  public :: crown_damage
  public :: cambial_damage_kill
  public :: post_fire_mortality

  integer :: write_SF = 0     ! for debugging
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains

  ! ============================================================================
  !        Area of site burned by fire           
  ! ============================================================================
  subroutine fire_model( currentSite, bc_in)

    use FatesInterfaceTypesMod, only : hlm_use_spitfire

    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in)            :: bc_in
    

    type (ed_patch_type), pointer :: currentPatch

    !zero fire things
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentPatch%frac_burnt = 0.0_r8
       currentPatch%fire       = 0
       currentPatch => currentPatch%older
    enddo

    if(write_SF==1)then
       write(fates_log(),*) 'use_spitfire',hlm_use_spitfire
    endif

    if( hlm_use_spitfire == itrue )then
       call fire_danger_index(currentSite, bc_in)
       call wind_effect(currentSite, bc_in) 
       call charecteristics_of_fuel(currentSite)
       call rate_of_spread(currentSite)
       call ground_fuel_consumption(currentSite)
       call area_burnt_intensity(currentSite)
       call crown_scorching(currentSite)
       call crown_damage(currentSite)
       call cambial_damage_kill(currentSite)
       call post_fire_mortality(currentSite)
    end if

  end subroutine fire_model

  !*****************************************************************
  subroutine  fire_danger_index ( currentSite, bc_in) 

   !*****************************************************************
   ! currentSite%acc_NI is the accumulated Nesterov fire danger index

    use SFParamsMod, only  : SF_val_fdi_a, SF_val_fdi_b
    use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod , only : sec_per_day

    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in)            :: bc_in

    real(r8) :: temp_in_C  ! daily averaged temperature in celcius
    real(r8) :: rainfall   ! daily precip in mm/day
    real(r8) :: rh         ! daily rh 
    
    real(r8) :: yipsolon   !intermediate varable for dewpoint calculation
    real(r8) :: dewpoint   !dewpoint in K 
    real(r8) :: d_NI       !daily change in Nesterov Index. C^2 
    integer  :: iofp       ! index of oldest the fates patch
  
    ! NOTE that the boundary conditions of temperature, precipitation and relative humidity
    ! are available at the patch level. We are currently using a simplification where the whole site
    ! is simply using the values associated with the first patch.
    ! which probably won't have much inpact, unless we decide to ever calculated the NI for each patch.  
    
    iofp = currentSite%oldest_patch%patchno
    
    temp_in_C  = bc_in%t_veg24_pa(iofp) - tfrz
    rainfall   = bc_in%precip24_pa(iofp)*sec_per_day
    rh         = bc_in%relhumid24_pa(iofp)
    
    if (rainfall > 3.0_r8) then !rezero NI if it rains... 
       d_NI = 0.0_r8
       currentSite%acc_NI = 0.0_r8
    else 
       yipsolon = (SF_val_fdi_a* temp_in_C)/(SF_val_fdi_b+ temp_in_C)+log(rh/100.0_r8) 
       dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a-yipsolon) !Standard met. formula
       d_NI = ( temp_in_C-dewpoint)* temp_in_C !follows Nesterov 1968.  Equation 5. Thonicke et al. 2010.
       if (d_NI < 0.0_r8) then !Change in NI cannot be negative. 
          d_NI = 0.0_r8 !check 
       endif
    endif
    currentSite%acc_NI = currentSite%acc_NI + d_NI        !Accumulate Nesterov index over the fire season. 

  end subroutine fire_danger_index


  !*****************************************************************
  subroutine  charecteristics_of_fuel ( currentSite )
  !*****************************************************************

    use SFParamsMod, only  : SF_val_drying_ratio, SF_val_SAV, SF_val_FBD

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort
    type(litter_type), pointer    :: litt_c

    real(r8) alpha_FMC(nfsc)     ! Relative fuel moisture adjusted per drying ratio
    real(r8) fuel_moisture(nfsc) ! Scaled moisture content of small litter fuels. 
    real(r8) MEF(nfsc)           ! Moisture extinction factor of fuels     integer n 

    fuel_moisture(:) = 0.0_r8
    
    

    currentPatch => currentSite%oldest_patch; 
    do while(associated(currentPatch))  

       litt_c => currentPatch%litter(element_pos(carbon12_element))
       
       ! How much live grass is there? 
       currentPatch%livegrass = 0.0_r8 
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))
          if(EDPftvarcon_inst%woody(currentCohort%pft) == 0)then 
             
             currentPatch%livegrass = currentPatch%livegrass + &
                  currentCohort%prt%GetState(leaf_organ, all_carbon_elements) * &
                  currentCohort%n/currentPatch%area

          endif
          currentCohort => currentCohort%shorter
       enddo
       
       ! There are SIX fuel classes
       ! 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
       ! NCWD =4  NFSC = 6
       ! tw_sf = 1, lb_sf = 3, tr_sf = 4, dl_sf = 5, lg_sf = 6,
     
       ! zero fire arrays. 
       currentPatch%fuel_eff_moist = 0.0_r8 
       currentPatch%fuel_bulkd     = 0.0_r8  !this is kgBiomass/m3 for use in rate of spread equations
       currentPatch%fuel_sav       = 0.0_r8 
       currentPatch%fuel_frac(:)   = 0.0_r8 
       currentPatch%fuel_mef       = 0.0_r8
       currentPatch%sum_fuel       = 0.0_r8
       currentPatch%fuel_frac      = 0.0_r8

       if(write_sf == itrue)then
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter1 ',sum(litt_c%leaf_fines(:))
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter2 ',sum(litt_c%ag_cwd(:))
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter3 ',currentPatch%livegrass
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' sum fuel', currentPatch%sum_fuel
       endif

       currentPatch%sum_fuel =  sum(litt_c%leaf_fines(:)) + &
                                sum(litt_c%ag_cwd(:)) + &
                                currentPatch%livegrass
       if(write_SF == itrue)then
          if ( hlm_masterproc == itrue ) write(fates_log(),*) 'sum fuel', currentPatch%sum_fuel,currentPatch%area
       endif
       ! ===============================================
       ! Average moisture, bulk density, surface area-volume and moisture extinction of fuel
       ! ================================================   
                  
       if (currentPatch%sum_fuel > 0.0) then        
          ! Fraction of fuel in litter classes
          currentPatch%fuel_frac(dl_sf)       = sum(litt_c%leaf_fines(:))/ currentPatch%sum_fuel
          currentPatch%fuel_frac(tw_sf:tr_sf) = litt_c%ag_cwd(:) / currentPatch%sum_fuel    

          if(write_sf == itrue)then
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff1 ',currentPatch%fuel_frac
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff2 ',currentPatch%fuel_frac
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff2a ', &
                  lg_sf,currentPatch%livegrass,currentPatch%sum_fuel
          endif

          currentPatch%fuel_frac(lg_sf)       = currentPatch%livegrass       / currentPatch%sum_fuel   
          MEF(1:nfsc)                         = 0.524_r8 - 0.066_r8 * log10(SF_val_SAV(1:nfsc)) 

          !--- weighted average of relative moisture content---
          ! Equation 6 in Thonicke et al. 2010. across twig, small branch, large branch, and dead leaves
          ! dead leaves and twigs included in 1hr pool per Thonicke (2010) 
          ! Calculate fuel moisture for trunks to hold value for fuel consumption
          alpha_FMC(tw_sf:dl_sf)      = SF_val_SAV(tw_sf:dl_sf)/SF_val_drying_ratio
          
          fuel_moisture(tw_sf:dl_sf)  = exp(-1.0_r8 * alpha_FMC(tw_sf:dl_sf) * currentSite%acc_NI) 
 
          if(write_SF == itrue)then
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff3 ',currentPatch%fuel_frac
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'fm ',fuel_moisture
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'csa ',currentSite%acc_NI
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'sfv ',alpha_FMC
          endif
          
          ! live grass moisture is a function of SAV and changes via Nesterov Index
          ! along the same relationship as the 1 hour fuels (live grass has same SAV as dead grass,
          ! but retains more moisture with this calculation.)
          fuel_moisture(lg_sf)        = exp(-1.0_r8 * ((SF_val_SAV(tw_sf)/SF_val_drying_ratio) * currentSite%acc_NI))          
 
          ! Average properties over the first three litter pools (twigs, s branches, l branches) 
          currentPatch%fuel_bulkd     = sum(currentPatch%fuel_frac(tw_sf:lb_sf) * SF_val_FBD(tw_sf:lb_sf))     
          currentPatch%fuel_sav       = sum(currentPatch%fuel_frac(tw_sf:lb_sf) * SF_val_SAV(tw_sf:lb_sf))              
          currentPatch%fuel_mef       = sum(currentPatch%fuel_frac(tw_sf:lb_sf) * MEF(tw_sf:lb_sf))              
          currentPatch%fuel_eff_moist = sum(currentPatch%fuel_frac(tw_sf:lb_sf) * fuel_moisture(tw_sf:lb_sf))         
          if(write_sf == itrue)then
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff4 ',currentPatch%fuel_eff_moist
          endif
          ! Add on properties of dead leaves and live grass pools (5 & 6)
          currentPatch%fuel_bulkd     = currentPatch%fuel_bulkd    + sum(currentPatch%fuel_frac(dl_sf:lg_sf) * SF_val_FBD(dl_sf:lg_sf))      
          currentPatch%fuel_sav       = currentPatch%fuel_sav      + sum(currentPatch%fuel_frac(dl_sf:lg_sf) * SF_val_SAV(dl_sf:lg_sf))
          currentPatch%fuel_mef       = currentPatch%fuel_mef      + sum(currentPatch%fuel_frac(dl_sf:lg_sf) * MEF(dl_sf:lg_sf))            
          currentPatch%fuel_eff_moist = currentPatch%fuel_eff_moist+ sum(currentPatch%fuel_frac(dl_sf:lg_sf) * fuel_moisture(dl_sf:lg_sf))

          ! Correct averaging for the fact that we are not using the trunks pool for fire ROS and intensity (5)
          ! Consumption of fuel in trunk pool does not influence fire ROS or intensity (Pyne 1996)
          currentPatch%fuel_bulkd     = currentPatch%fuel_bulkd     * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_sav       = currentPatch%fuel_sav       * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_mef       = currentPatch%fuel_mef       * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_eff_moist = currentPatch%fuel_eff_moist * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf))) 
     
          ! Pass litter moisture into the fuel burning routine (all fuels: twigs,s branch,l branch,trunk,dead leaves,live grass)
          ! (wo/me term in Thonicke et al. 2010) 
          currentPatch%litter_moisture(tw_sf:lb_sf) = fuel_moisture(tw_sf:lb_sf)/MEF(tw_sf:lb_sf)   
          currentPatch%litter_moisture(tr_sf)       = fuel_moisture(tr_sf)/MEF(tr_sf)
          currentPatch%litter_moisture(dl_sf)       = fuel_moisture(dl_sf)/MEF(dl_sf)
          currentPatch%litter_moisture(lg_sf)       = fuel_moisture(lg_sf)/MEF(lg_sf)
          
       else

          if(write_SF == itrue)then

             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'no litter fuel at all',currentPatch%patchno, &
                  currentPatch%sum_fuel,sum(litt_c%ag_cwd(:)),sum(litt_c%leaf_fines(:))

          endif
          currentPatch%fuel_sav = sum(SF_val_SAV(1:nfsc))/(nfsc) ! make average sav to avoid crashing code. 

          if ( hlm_masterproc == itrue ) write(fates_log(),*) 'problem with spitfire fuel averaging'

          ! FIX(SPM,032414) refactor...should not have 0 fuel unless everything is burnt
          ! off.
          currentPatch%fuel_eff_moist = 0.0000000001_r8 
          currentPatch%fuel_bulkd     = 0.0000000001_r8 
          currentPatch%fuel_frac(:)   = 0.0000000001_r8 
          currentPatch%fuel_mef       = 0.0000000001_r8
          currentPatch%sum_fuel       = 0.0000000001_r8
          currentPatch%fuel_frac      = 0.0000000001_r8

       endif
       ! check values. 
       ! FIX(SPM,032414) refactor...
       if(write_SF == itrue.and.currentPatch%fuel_sav <= 0.0_r8.or.currentPatch%fuel_bulkd <=  &
            0.0_r8.or.currentPatch%fuel_mef <= 0.0_r8.or.currentPatch%fuel_eff_moist <= 0.0_r8)then
            if ( hlm_masterproc == itrue ) write(fates_log(),*) 'problem with spitfire fuel averaging'
       endif 
       
       currentPatch => currentPatch%younger

    enddo !end patch loop
    
  end subroutine charecteristics_of_fuel


  !*****************************************************************
  subroutine  wind_effect ( currentSite, bc_in) 
  !*****************************************************************.

    ! Routine called daily from within ED within a site loop.
    ! Calculates the effective windspeed based on vegetation charecteristics.
    ! currentSite%wind is daily wind converted to m/min for Spitfire units 

    use FatesConstantsMod, only : sec_per_min

    type(ed_site_type) , intent(inout), target :: currentSite
    type(bc_in_type)   , intent(in)            :: bc_in

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    real(r8) :: total_grass_area     ! per patch,in m2
    real(r8) :: tree_fraction        ! site level. no units
    real(r8) :: grass_fraction       ! site level. no units
    real(r8) :: bare_fraction        ! site level. no units 
    integer  :: iofp                 ! index of oldest fates patch


    ! note - this is a patch level temperature, which probably won't have much inpact, 
    ! unless we decide to ever calculated the NI for each patch.  

    iofp = currentSite%oldest_patch%patchno
    currentSite%wind = bc_in%wind24_pa(iofp) * sec_per_min !Convert to m/min for SPITFIRE

    if(write_SF == itrue)then
       if ( hlm_masterproc == itrue ) write(fates_log(),*) 'wind24', currentSite%wind
    endif
    ! --- influence of wind speed, corrected for surface roughness----
    ! --- averaged over the whole grid cell to prevent extreme divergence 
    ! average_wspeed = 0.0_r8   
    tree_fraction = 0.0_r8
    grass_fraction = 0.0_r8
    currentPatch=>currentSite%oldest_patch;  
    do while(associated(currentPatch))
       currentPatch%total_tree_area = 0.0_r8
       total_grass_area = 0.0_r8
       currentCohort => currentPatch%tallest
 
       do while(associated(currentCohort))
          if (debug) write(fates_log(),*) 'SF currentCohort%c_area ',currentCohort%c_area
          if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
             currentPatch%total_tree_area = currentPatch%total_tree_area + currentCohort%c_area
          else
             total_grass_area = total_grass_area + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo
       tree_fraction = tree_fraction + min(currentPatch%area,currentPatch%total_tree_area)/AREA
       grass_fraction = grass_fraction + min(currentPatch%area,total_grass_area)/AREA 
       
       if(debug)then
         write(fates_log(),*) 'SF  currentPatch%area ',currentPatch%area
         write(fates_log(),*) 'SF  currentPatch%total_area ',currentPatch%total_tree_area
         write(fates_log(),*) 'SF  total_grass_area ',tree_fraction,grass_fraction
         write(fates_log(),*) 'SF  AREA ',AREA
       endif
       
       currentPatch => currentPatch%younger
    enddo !currentPatch loop

    !if there is a cover of more than one, then the grasses are under the trees
    grass_fraction = min(grass_fraction,1.0_r8-tree_fraction) 
    bare_fraction = 1.0_r8 - tree_fraction - grass_fraction
    if(write_sf == itrue)then
       if ( hlm_masterproc == itrue ) write(fates_log(),*) 'grass, trees, bare', &
            grass_fraction, tree_fraction, bare_fraction
    endif

    currentPatch=>currentSite%oldest_patch;

    do while(associated(currentPatch))       
       currentPatch%total_tree_area = min(currentPatch%total_tree_area,currentPatch%area)
       ! effect_wspeed in units m/min      
       currentPatch%effect_wspeed = currentSite%wind * (tree_fraction*0.4_r8+(grass_fraction+bare_fraction)*0.6_r8)
      
       currentPatch => currentPatch%younger
    enddo !end patch loop

  end subroutine wind_effect

  !*****************************************************************
  subroutine rate_of_spread ( currentSite ) 
    !*****************************************************************.
    !Routine called daily from within ED within a site loop.
    !Returns the updated currentPatch%ROS_front value for each patch.

    use SFParamsMod, only  : SF_val_miner_total, &
                             SF_val_part_dens,   &
                             SF_val_miner_damp,  &
                             SF_val_fuel_energy
    
    use FatesInterfaceTypesMod, only : hlm_current_day, hlm_current_month

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type), pointer :: currentPatch

    ! Rothermal fire spread model parameters. 
    real(r8) beta,beta_op         ! weighted average of packing ratio (unitless)
    real(r8) ir                   ! reaction intensity (kJ/m2/min)
    real(r8) xi,eps,phi_wind      ! all are unitless
    real(r8) q_ig                 ! heat of pre-ignition (kJ/kg)
    real(r8) reaction_v_opt,reaction_v_max !reaction velocity (per min)!optimum and maximum
    real(r8) moist_damp,mw_weight ! moisture dampening coefficient and ratio fuel moisture to extinction
    real(r8) beta_ratio           ! ratio of beta/beta_op
    real(r8) a_beta               ! dummy variable for product of a* beta_ratio for react_v_opt equation
    real(r8) a,b,c,e              ! function of fuel sav

    logical, parameter :: debug_windspeed = .false. !for debugging
    real(r8),parameter :: q_dry = 581.0_r8          !heat of pre-ignition of dry fuels (kJ/kg) 

    currentPatch=>currentSite%oldest_patch;  

    do while(associated(currentPatch))
              
        ! ---initialise parameters to zero.--- 
       beta_ratio = 0.0_r8; q_ig = 0.0_r8; eps = 0.0_r8;   a = 0.0_r8;   b = 0.0_r8;   c = 0.0_r8;   e = 0.0_r8
       phi_wind = 0.0_r8;   xi = 0.0_r8;   reaction_v_max = 0.0_r8;  reaction_v_opt = 0.0_r8; mw_weight = 0.0_r8
       moist_damp = 0.0_r8;   ir = 0.0_r8; a_beta = 0.0_r8;     
       currentPatch%ROS_front = 0.0_r8

       ! remove mineral content from net fuel load per Thonicke 2010 for ir calculation
       currentPatch%sum_fuel  = currentPatch%sum_fuel * (1.0_r8 - SF_val_miner_total) !net of minerals

       ! ----start spreading---

       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) &
            'SF - currentPatch%fuel_bulkd ',currentPatch%fuel_bulkd
       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) &
            'SF - SF_val_part_dens ',SF_val_part_dens

       ! beta = packing ratio (unitless)
       ! fraction of fuel array volume occupied by fuel or compactness of fuel bed 
       beta = currentPatch%fuel_bulkd / SF_val_part_dens
       
       ! Equation A6 in Thonicke et al. 2010
       ! packing ratio (unitless) 
       beta_op = 0.200395_r8 *(currentPatch%fuel_sav**(-0.8189_r8))

       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta ',beta
       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta_op ',beta_op
       beta_ratio = beta/beta_op   !unitless

       if(write_sf == itrue)then
          if ( hlm_masterproc == itrue ) write(fates_log(),*) 'esf ',currentPatch%fuel_eff_moist
       endif

       ! ---heat of pre-ignition---
       !  Equation A4 in Thonicke et al. 2010
       !  Rothermal EQ12= 250 Btu/lb + 1116 Btu/lb * fuel_eff_moist
       !  conversion of Rothermal (1972) EQ12 in BTU/lb to current kJ/kg 
       !  q_ig in kJ/kg 
       q_ig = q_dry +2594.0_r8 * currentPatch%fuel_eff_moist

       ! ---effective heating number---
       ! Equation A3 in Thonicke et al. 2010.  
       eps = exp(-4.528_r8 / currentPatch%fuel_sav)     
       ! Equation A7 in Thonicke et al. 2010
       b = 0.15988_r8 * (currentPatch%fuel_sav**0.54_r8)
       ! Equation A8 in Thonicke et al. 2010
       c = 7.47_r8 * (exp(-0.8711_r8 * (currentPatch%fuel_sav**0.55_r8)))
       ! Equation A9 in Thonicke et al. 2010. 
       e = 0.715_r8 * (exp(-0.01094_r8 * currentPatch%fuel_sav))

       if (debug) then
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - c ',c
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - currentPatch%effect_wspeed ', &
                                                                         currentPatch%effect_wspeed
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - b ',b
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta_ratio ',beta_ratio
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - e ',e
       endif

       ! Equation A5 in Thonicke et al. 2010
       ! phi_wind (unitless)
       ! convert current_wspeed (wind at elev relevant to fire) from m/min to ft/min for Rothermel ROS eqn
       phi_wind = c * ((3.281_r8*currentPatch%effect_wspeed)**b)*(beta_ratio**(-e))


       ! ---propagating flux----
       ! Equation A2 in Thonicke et al.2010 and Eq. 42 Rothermal 1972
       ! xi (unitless)       
       xi = (exp((0.792_r8 + 3.7597_r8 * (currentPatch%fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
            (192_r8+7.9095_r8 * currentPatch%fuel_sav)      
      
       ! ---reaction intensity----
       ! Equation in table A1 Thonicke et al. 2010. 
       a = 8.9033_r8 * (currentPatch%fuel_sav**(-0.7913_r8))
       a_beta = exp(a*(1.0_r8-beta_ratio))  !dummy variable for reaction_v_opt equation
  
       ! Equation in table A1 Thonicke et al. 2010.
       ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
       ! reaction_v_max = Equation 36 in Rothermal 1972 and Fig 12 
       reaction_v_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (currentPatch%fuel_sav**(-1.5_r8)))
       ! reaction_v_opt =  Equation 38 in Rothermal 1972 and Fig 11
       reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta

       ! mw_weight = relative fuel moisture/fuel moisture of extinction
       ! average values for litter pools (dead leaves, twigs, small and large branches) plus grass
       mw_weight = currentPatch%fuel_eff_moist/currentPatch%fuel_mef
       
       ! Equation in table A1 Thonicke et al. 2010. 
       ! moist_damp is unitless
       moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
            (3.52_r8*(mw_weight**3.0_r8))))

       ! FIX(SPM, 040114) ask RF if this should be an endrun
       ! if(write_SF == itrue)then
       ! write(fates_log(),*) 'moist_damp' ,moist_damp,mw_weight,currentPatch%fuel_eff_moist,currentPatch%fuel_mef
       ! endif
       
       ! ir = reaction intenisty in kJ/m2/min
       ! currentPatch%sum_fuel converted from kgC/m2 to kgBiomass/m2 for ir calculation
       ir = reaction_v_opt*(currentPatch%sum_fuel/0.45_r8)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp 

       ! write(fates_log(),*) 'ir',gamma_aptr,moist_damp,SF_val_fuel_energy,SF_val_miner_damp

       if (((currentPatch%fuel_bulkd) <= 0.0_r8).or.(eps <= 0.0_r8).or.(q_ig <= 0.0_r8)) then
          currentPatch%ROS_front = 0.0_r8
       else ! Equation 9. Thonicke et al. 2010. 
            ! forward ROS in m/min
          currentPatch%ROS_front = (ir*xi*(1.0_r8+phi_wind)) / (currentPatch%fuel_bulkd*eps*q_ig)
          ! write(fates_log(),*) 'ROS',currentPatch%ROS_front,phi_wind,currentPatch%effect_wspeed
          ! write(fates_log(),*) 'ros calcs',currentPatch%fuel_bulkd,ir,xi,eps,q_ig
       endif
       ! Equation 10 in Thonicke et al. 2010
       ! backward ROS from Can FBP System (1992) in m/min
       ! backward ROS wind not changed by vegetation 
       currentPatch%ROS_back = currentPatch%ROS_front*exp(-0.012_r8*currentSite%wind) 

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine  rate_of_spread

  !*****************************************************************
  subroutine  ground_fuel_consumption ( currentSite ) 
  !*****************************************************************
    !returns the  the hypothetic fuel consumed by the fire

    use SFParamsMod, only : SF_val_miner_total, SF_val_min_moisture, &
         SF_val_mid_moisture, SF_val_low_moisture_Coeff, SF_val_low_moisture_Slope, &
         SF_val_mid_moisture_Coeff, SF_val_mid_moisture_Slope

    type(ed_site_type) , intent(in), target :: currentSite
    type(ed_patch_type), pointer    :: currentPatch
    type(litter_type), pointer      :: litt_c           ! carbon 12 litter pool
    
    real(r8) :: moist           !effective fuel moisture
    real(r8) :: tau_b(nfsc)     !lethal heating rates for each fuel class (min) 
    real(r8) :: fc_ground(nfsc) !proportion of fuel consumed

    integer  :: c

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch))
       currentPatch%burnt_frac_litter = 1.0_r8       
       ! Calculate fraction of litter is burnt for all classes. 
       ! Equation B1 in Thonicke et al. 2010---
       do c = 1, nfsc    !work out the burnt fraction for all pools, even if those pools dont exist.         
          moist = currentPatch%litter_moisture(c)                  
          ! 1. Very dry litter
          if (moist <= SF_val_min_moisture(c)) then
             currentPatch%burnt_frac_litter(c) = 1.0_r8  
          endif
          ! 2. Low to medium moistures
          if (moist > SF_val_min_moisture(c).and.moist <= SF_val_mid_moisture(c)) then
             currentPatch%burnt_frac_litter(c) = max(0.0_r8,min(1.0_r8,SF_val_low_moisture_Coeff(c)- &
                  SF_val_low_moisture_Slope(c)*moist)) 
          else
          ! For medium to high moistures. 
             if (moist > SF_val_mid_moisture(c).and.moist <= 1.0_r8) then
                currentPatch%burnt_frac_litter(c) = max(0.0_r8,min(1.0_r8,SF_val_mid_moisture_Coeff(c)- &
                     SF_val_mid_moisture_Slope(c)*moist))
             endif

          endif
          ! Very wet litter        
          if (moist >= 1.0_r8) then !this shouldn't happen? 
             currentPatch%burnt_frac_litter(c) = 0.0_r8  
          endif          
       enddo !c   

       ! we can't ever kill -all- of the grass. 
       currentPatch%burnt_frac_litter(lg_sf) = min(0.8_r8,currentPatch%burnt_frac_litter(lg_sf ))  
       ! reduce burnt amount for mineral content. 
       currentPatch%burnt_frac_litter = currentPatch%burnt_frac_litter * (1.0_r8-SF_val_miner_total) 

       !---Calculate amount of fuel burnt.---    

       litt_c => currentPatch%litter(element_pos(carbon12_element))
       FC_ground(tw_sf:tr_sf) = currentPatch%burnt_frac_litter(tw_sf:tr_sf) * litt_c%ag_cwd(tw_sf:tr_sf)
       FC_ground(dl_sf)       = currentPatch%burnt_frac_litter(dl_sf)   * sum(litt_c%leaf_fines(:))
       FC_ground(lg_sf)       = currentPatch%burnt_frac_litter(lg_sf)   * currentPatch%livegrass      

       ! Following used for determination of cambial kill follows from Peterson & Ryan (1986) scheme 
       ! less empirical cf current scheme used in SPITFIRE which attempts to mesh Rothermel 
       ! and P&R, and while solving potential inconsistencies, actually results in BIG values for 
       ! fire residence time, thus lots of vegetation death!   
       ! taul is the duration of the lethal heating.  
       ! The /10 is to convert from kgC/m2 into gC/cm2, as in the Peterson and Ryan paper #Rosie,Jun 2013
        
       do c = 1,nfsc  
          tau_b(c)   =  39.4_r8 *(currentPatch%fuel_frac(c)*currentPatch%sum_fuel/0.45_r8/10._r8)* &
               (1.0_r8-((1.0_r8-currentPatch%burnt_frac_litter(c))**0.5_r8))  
       enddo
       tau_b(tr_sf)   =  0.0_r8
       ! Cap the residence time to 8mins, as suggested by literature survey by P&R (1986).
       currentPatch%tau_l = min(8.0_r8,sum(tau_b)) 

       !---calculate overall fuel consumed by spreading fire --- 
       ! ignore 1000hr fuels. Just interested in fuels affecting ROS   
       currentPatch%TFC_ROS = sum(FC_ground)-FC_ground(tr_sf)  

       currentPatch=>currentPatch%younger;
    enddo !end patch loop

  end subroutine ground_fuel_consumption

  
  !*****************************************************************
  subroutine  area_burnt_intensity ( currentSite ) 
  !*****************************************************************

    !returns the updated currentPatch%FI value for each patch.

    !currentPatch%FI  avg fire intensity of flaming front during day. Backward ROS plays no role here. kJ/m/s or kW/m.
    !currentSite%FDI  probability that an ignition will start a fire
    !currentSite%NF   number of lighting strikes per day per km2
    !currentPatch%ROS_front  forward ROS (m/min) 
    !currentPatch%TFC_ROS total fuel consumed by flaming front (kgC/m2)

    use FatesInterfaceTypesMod, only : hlm_use_spitfire
    use EDParamsMod,       only : ED_val_nignitions
    use EDParamsMod,       only : cg_strikes    ! fraction of cloud-to-ground ligtning strikes
    use FatesConstantsMod, only : years_per_day
    use SFParamsMod,       only : SF_val_fdi_alpha,SF_val_fuel_energy, &
         SF_val_max_durat, SF_val_durat_slope, SF_val_fire_threshold
    
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type), pointer :: currentPatch

    real(r8) ROS !m/s
    real(r8) W   !kgBiomass/m2
    real(r8) lb               !length to breadth ratio of fire ellipse (unitless)
    real(r8) df               !distance fire has travelled forward in m
    real(r8) db               !distance fire has travelled backward in m
    real(r8) AB               !daily area burnt in m2 per km2
    
    real(r8) size_of_fire !in m2
    real(r8),parameter :: km2_to_m2 = 1000000.0_r8 !area conversion for square km to square m

    !  ---initialize site parameters to zero--- 
    currentSite%frac_burnt = 0.0_r8  

    
    ! Equation 7 from Venevsky et al GCB 2002 (modification of equation 8 in Thonicke et al. 2010) 
    ! FDI 0.1 = low, 0.3 moderate, 0.75 high, and 1 = extreme ignition potential for alpha 0.000337
    currentSite%FDI  = 1.0_r8 - exp(-SF_val_fdi_alpha*currentSite%acc_NI)
    
    !NF = number of lighting strikes per day per km2 scaled by cloud to ground strikes
    currentSite%NF = ED_val_nignitions * years_per_day * cg_strikes

    ! If there are 15  lightning strikes per year, per km2. (approx from NASA product for S.A.) 
    ! then there are 15 * 1/365 strikes/km2 each day 
 

    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch))
       !  ---initialize patch parameters to zero---
       currentPatch%fire       = 0
       currentPatch%FD         = 0.0_r8
       currentPatch%frac_burnt = 0.0_r8
       

       if (currentSite%NF > 0.0_r8) then
          
          ! Equation 14 in Thonicke et al. 2010
          ! fire duration in minutes
          currentPatch%FD = (SF_val_max_durat+1.0_r8) / (1.0_r8 + SF_val_max_durat * &
                            exp(SF_val_durat_slope*currentSite%FDI))
          if(write_SF == itrue)then
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'fire duration minutes',currentPatch%fd
          endif
          !equation 15 in Arora and Boer CTEM model.Average fire is 1 day long.
          !currentPatch%FD = 60.0_r8 * 24.0_r8 !no minutes in a day

           
       ! The feedback between vegetation structure and ellipse size if turned off for now, 
       ! to reduce the positive feedback in the syste,
       ! This will also be investigated by William Hoffmans proposal. 
          !      if (currentPatch%effect_wspeed < 16.67_r8) then !16.67m/min = 1km/hr 
          lb = 1.0_r8
          !     else 
          !FIX(RF,032414) FOR NO GRASS
          !        lb = currentPatch%total_canopy_area/currentPatch%area*(1.0_r8)+(8.729_r8 * &
          ! ((1.0_r8 -(exp(-0.03_r8 * 0.06_r8 * currentPatch%effect_wspeed)))**2.155_r8)) !&
          !&       +currentPatch%fpc_grass*(1.1_r8+((0.06_r8*currentPatch%effect_wspeed)**0.0464))

          !      endif

          !     if (lb > 8.0_r8)then
          !       lb = 8.0_r8  !Constraint Canadian Fire Behaviour System
          !     endif
          ! ---- calculate length of major axis---
          db = currentPatch%ROS_back  * currentPatch%FD !m
          df = currentPatch%ROS_front * currentPatch%FD !m

          ! --- calculate area burnt---
          if(lb > 0.0_r8) then
     
             ! Equation 1 in Thonicke et al. 2010
             ! To Do: Connect here with the Li & Levis GDP fire suppression algorithm. 
             ! Equation 16 in arora and boer model JGR 2005
             ! AB = AB *3.0_r8

             !size of fire = equation 14 Arora and Boer JGR 2005
             size_of_fire = ((3.1416_r8/(4.0_r8*lb))*((df+db)**2.0_r8))

             !AB = daily area burnt = size fires in m2 * num ignitions per day per km2 * prob ignition starts fire
             !AB = m2 per km2 per day
             AB = size_of_fire * currentSite%NF * currentSite%FDI

             !frac_burnt 
             currentPatch%frac_burnt = (min(0.99_r8, AB / km2_to_m2)) * currentPatch%area/area 
             
             if(write_SF == itrue)then
                if ( hlm_masterproc == itrue ) write(fates_log(),*) 'frac_burnt',currentPatch%frac_burnt
             endif

          endif ! lb

         ROS   = currentPatch%ROS_front / 60.0_r8 !m/min to m/sec 
         W     = currentPatch%TFC_ROS / 0.45_r8 !kgC/m2 to kgbiomass/m2          

         ! EQ 15 Thonicke et al 2010
         !units of fire intensity = (kJ/kg)*(kgBiomass/m2)*(m/min)*unitless_fraction
         currentPatch%FI = SF_val_fuel_energy * W * ROS * currentPatch%frac_burnt !kj/m/s, or kW/m
       
         if(write_sf == itrue)then
             if( hlm_masterproc == itrue ) write(fates_log(),*) 'fire_intensity',currentPatch%fi,W,currentPatch%ROS_front
         endif

         !'decide_fire' subroutine 
         if (currentPatch%FI > SF_val_fire_threshold) then !track fires greater than kW/m energy threshold
            currentPatch%fire = 1 ! Fire...    :D
          
         else     
            currentPatch%fire       = 0 ! No fire... :-/
            currentPatch%FD         = 0.0_r8
            currentPatch%frac_burnt = 0.0_r8
         endif         
          
       endif! NF ignitions check

       
       ! accumulate frac_burnt % at site level
       currentSite%frac_burnt = currentSite%frac_burnt + currentPatch%frac_burnt    

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine area_burnt_intensity



  !*****************************************************************
  subroutine  crown_scorching ( currentSite ) 
  !*****************************************************************

    !currentPatch%FI       average fire intensity of flaming front during day.  kW/m.
    !currentPatch%SH(pft)  scorch height for all cohorts of a given PFT on a given patch (m)

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type), pointer  :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    real(r8) ::  tree_ag_biomass ! total amount of above-ground tree biomass in patch. kgC/m2
    real(r8) ::  leaf_c          ! leaf carbon      [kg]
    real(r8) ::  sapw_c          ! sapwood carbon   [kg]
    real(r8) ::  struct_c        ! structure carbon [kg]

    integer  ::  i_pft


    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch)) 
       
       tree_ag_biomass = 0.0_r8
       if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if (EDPftvarcon_inst%woody(currentCohort%pft) == 1) then !trees only
                
                leaf_c = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                sapw_c = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
                
                tree_ag_biomass = tree_ag_biomass + &
                     currentCohort%n * (leaf_c + & 
                     EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))
 
             endif !trees only
             currentCohort=>currentCohort%shorter;
          enddo !end cohort loop

          do i_pft=1,numpft
             if (tree_ag_biomass > 0.0_r8  .and. EDPftvarcon_inst%woody(i_pft) == 1) then 
                
                !Equation 16 in Thonicke et al. 2010 !Van Wagner 1973 EQ8 !2/3 Byram (1959)
                currentPatch%Scorch_ht(i_pft) = EDPftvarcon_inst%fire_alpha_SH(i_pft) * (currentPatch%FI**0.667_r8)
                
                if(write_SF == itrue)then
                   if ( hlm_masterproc == itrue ) write(fates_log(),*) 'currentPatch%SH',currentPatch%Scorch_ht(i_pft)
                endif
             else
                currentPatch%Scorch_ht(i_pft) = 0.0_r8
             endif ! tree biomass
          end do

       endif !fire

       currentPatch => currentPatch%younger;  
    enddo !end patch loop

  end subroutine crown_scorching

  !*****************************************************************
  subroutine  crown_damage ( currentSite )
    !*****************************************************************

    !returns the updated currentCohort%fraction_crown_burned for each tree cohort within each patch.
    !currentCohort%fraction_crown_burned is the proportion of crown affected by fire

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch)) 
       if (currentPatch%fire == 1) then

          currentCohort=>currentPatch%tallest

          do while(associated(currentCohort))  
             currentCohort%fraction_crown_burned = 0.0_r8
             if (EDPftvarcon_inst%woody(currentCohort%pft) == 1) then !trees only
                ! Flames lower than bottom of canopy. 
                ! c%hite is height of cohort
                if (currentPatch%Scorch_ht(currentCohort%pft) < &
                     (currentCohort%hite-currentCohort%hite*EDPftvarcon_inst%crown(currentCohort%pft))) then 
                   currentCohort%fraction_crown_burned = 0.0_r8
                else
                   ! Flames part of way up canopy. 
                   ! Equation 17 in Thonicke et al. 2010. 
                   ! flames over bottom of canopy but not over top.
                   if ((currentCohort%hite > 0.0_r8).and.(currentPatch%Scorch_ht(currentCohort%pft) >=  &
                        (currentCohort%hite-currentCohort%hite*EDPftvarcon_inst%crown(currentCohort%pft)))) then 

                        currentCohort%fraction_crown_burned = (currentPatch%Scorch_ht(currentCohort%pft) - &
                                currentCohort%hite*(1.0_r8 - &
                                EDPftvarcon_inst%crown(currentCohort%pft)))/(currentCohort%hite* &
                                EDPftvarcon_inst%crown(currentCohort%pft)) 

                   else 
                      ! Flames over top of canopy. 
                      currentCohort%fraction_crown_burned =  1.0_r8 
                   endif

                endif
                ! Check for strange values. 
                currentCohort%fraction_crown_burned = min(1.0_r8, max(0.0_r8,currentCohort%fraction_crown_burned))              
             endif !trees only
             !shrink canopy to account for burnt section.     
             !currentCohort%canopy_trim = min(currentCohort%canopy_trim,(1.0_r8-currentCohort%fraction_crown_burned)) 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine crown_damage

  !*****************************************************************
  subroutine  cambial_damage_kill ( currentSite ) 
    !*****************************************************************
    ! routine description.
    ! returns the probability that trees dies due to cambial char
    ! currentPatch%tau_l = duration of lethal stem heating (min). Calculated at patch level.

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    real(r8) :: tau_c !critical time taken to kill cambium (minutes) 
    real(r8) :: bt    !bark thickness in cm.

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch)) 

       if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if (EDPftvarcon_inst%woody(currentCohort%pft) == 1) then !trees only
                ! Equation 21 in Thonicke et al 2010
                bt = EDPftvarcon_inst%bark_scaler(currentCohort%pft)*currentCohort%dbh ! bark thickness. 
                ! Equation 20 in Thonicke et al. 2010. 
                tau_c = 2.9_r8*bt**2.0_r8 !calculate time it takes to kill cambium (min)
                ! Equation 19 in Thonicke et al. 2010
                if ((currentPatch%tau_l/tau_c) >= 2.0_r8) then
                   currentCohort%cambial_mort = 1.0_r8
                else
                   if ((currentPatch%tau_l/tau_c) > 0.22_r8) then
                      currentCohort%cambial_mort = (0.563_r8*(currentPatch%tau_l/tau_c)) - 0.125_r8
                   else
                      currentCohort%cambial_mort = 0.0_r8
                   endif
                endif
             endif !trees 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?

       currentPatch=>currentPatch%younger;

    enddo !end patch loop

  end subroutine cambial_damage_kill

  !*****************************************************************
  subroutine  post_fire_mortality ( currentSite )
  !*****************************************************************

    !  returns the updated currentCohort%fire_mort value for each tree cohort within each patch.
    !  currentCohort%fraction_crown_burned is proportion of crown affected by fire
    !  currentCohort%crownfire_mort  probability of tree post-fire mortality due to crown scorch
    !  currentCohort%cambial_mort  probability of tree post-fire mortality due to cambial char
    !  currentCohort%fire_mort  post-fire mortality from cambial and crown damage assuming two are independent.

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch)) 

       if (currentPatch%fire == 1) then 
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))  
             currentCohort%fire_mort = 0.0_r8
             currentCohort%crownfire_mort = 0.0_r8
             if (EDPftvarcon_inst%woody(currentCohort%pft) == 1) then
                ! Equation 22 in Thonicke et al. 2010. 
                currentCohort%crownfire_mort = EDPftvarcon_inst%crown_kill(currentCohort%pft)*currentCohort%fraction_crown_burned**3.0_r8
                ! Equation 18 in Thonicke et al. 2010. 
                currentCohort%fire_mort = max(0._r8,min(1.0_r8,currentCohort%crownfire_mort+currentCohort%cambial_mort- &
                     (currentCohort%crownfire_mort*currentCohort%cambial_mort)))  !joint prob.   
             else
                currentCohort%fire_mort = 0.0_r8 !Set to zero. Grass mode of death is removal of leaves.
             endif !trees

             currentCohort => currentCohort%shorter

          enddo !end cohort loop
       endif !fire?

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine post_fire_mortality

  ! ============================================================================
end module SFMainMod
