  module SFMainMod

  ! ============================================================================
  ! All subroutines realted to the SPITFIRE fire routine. 
  ! Code originally developed by Allan Spessa & Rosie Fisher as part of the NERC-QUEST project.  
  ! ============================================================================

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : pi_const
  use FatesInterfaceTypesMod     , only : hlm_masterproc ! 1= master process, 0=not master process
  use EDTypesMod            , only : numWaterMem
  use FatesGlobals          , only : fates_log
  use FatesInterfaceTypesMod, only : hlm_spitfire_mode
  use FatesInterfaceTypesMod, only : hlm_sf_nofire_def
  use FatesInterfaceTypesMod, only : hlm_sf_scalar_lightning_def
  use FatesInterfaceTypesMod, only : hlm_sf_successful_ignitions_def
  use FatesInterfaceTypesMod, only : hlm_sf_anthro_ignitions_def
  use FatesInterfaceTypesMod, only : bc_in_type
  
  use EDPftvarcon           , only : EDPftvarcon_inst
  use PRTParametersMod      , only : prt_params
  
  use PRTGenericMod         , only : element_pos
  use EDtypesMod            , only : ed_site_type
  use EDtypesMod            , only : ed_patch_type
  use EDtypesMod            , only : ed_cohort_type
  use EDtypesMod            , only : AREA
  use EDtypesMod            , only : DL_SF
! use EDtypesMod            , only : crown_fire_threshold  ! TODO slevis: NOT used; if end up using, look at SF_val_fire_threshold as template
  use EDTypesMod            , only : TW_SF
  use EDtypesMod            , only : LB_SF
  use EDtypesMod            , only : LG_SF
  use FatesLitterMod        , only : ncwd
  use EDtypesMod            , only : NFSC
  use EDtypesMod            , only : TR_SF
  use FatesLitterMod        , only : litter_type

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
  use FatesAllometryMod,      only : CrownDepth

  implicit none
  private

  public :: fire_model
  public :: fire_danger_index 
  public :: characteristics_of_fuel
  public :: characteristics_of_crown
  public :: rate_of_spread
  public :: ground_fuel_consumption
  public :: wind_effect
  public :: area_burnt_intensity
  public :: active_crown_fire
  public :: crown_scorching
  public :: crown_damage
  public :: cambial_damage_kill
  public :: post_fire_mortality

  ! The following parameter represents one of the values of hlm_spitfire_mode
  ! and more of these appear in subroutine area_burnt_intensity below
  ! NB. The same parameters are set in /src/biogeochem/CNFireFactoryMod
  integer :: write_SF = 0     ! for debugging
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains

  ! ============================================================================
  !        Area of site burned by fire           
  ! ============================================================================
  subroutine fire_model( currentSite, bc_in)

    

    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in)            :: bc_in
    

    type (ed_patch_type), pointer :: currentPatch


    !zero fire things
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentPatch%frac_burnt = 0.0_r8
       currentPatch%FI         = 0.0_r8
       currentPatch%FD         = 0.0_r8
       currentPatch%fire       = 0
       currentPatch%active_crown_fire_flg = 0
       currentPatch => currentPatch%older
    enddo

    if(write_SF==1)then
       write(fates_log(),*) 'spitfire_mode', hlm_spitfire_mode
    endif

    if( hlm_spitfire_mode > hlm_sf_nofire_def )then
       call fire_danger_index(currentSite, bc_in)
       call wind_effect(currentSite, bc_in) 
       call characteristics_of_fuel(currentSite)
       call characteristics_of_crown(currentSite)
       call rate_of_spread(currentSite)
       call ground_fuel_consumption(currentSite)
       call area_burnt_intensity(currentSite, bc_in)
       call active_crown_fire (currentSite)
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
    
    temp_in_C  = currentSite%oldest_patch%tveg24%GetMean() - tfrz
    rainfall   = bc_in%precip24_pa(iofp)*sec_per_day
    rh         = bc_in%relhumid24_pa(iofp)
    
    if (rainfall > 3.0_r8) then !rezero NI if it rains... 
       d_NI = 0.0_r8
       currentSite%acc_NI = 0.0_r8
    else 
       yipsolon = (SF_val_fdi_a* temp_in_C)/(SF_val_fdi_b+ temp_in_C)+log(max(1.0_r8,rh)/100.0_r8) 
       dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a-yipsolon) !Standard met. formula
       d_NI = ( temp_in_C-dewpoint)* temp_in_C !follows Nesterov 1968.  Eq 5, Thonicke et al. 2010.
       if (d_NI < 0.0_r8) then !Change in NI cannot be negative. 
          d_NI = 0.0_r8 !check 
       endif
    endif
    currentSite%acc_NI = currentSite%acc_NI + d_NI        !Accumulate Nesterov index over fire season. 

  end subroutine fire_danger_index


  !*****************************************************************
  subroutine  characteristics_of_fuel ( currentSite)
  !*****************************************************************

    use SFParamsMod, only: SF_val_drying_ratio, SF_val_SAV, SF_val_FBD, &
                           SF_val_miner_total

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort
    type(litter_type), pointer    :: litt_c

    real(r8) alpha_FMC(nfsc)     ! Relative fuel moisture adjusted per drying ratio
    real(r8) fuel_moisture(nfsc) ! Scaled moisture content of small litter fuels
    real(r8) MEF(nfsc)           ! Moisture extinction factor of fuels, integer n

    fuel_moisture(:) = 0.0_r8
    
    

    currentPatch => currentSite%oldest_patch; 
    do while(associated(currentPatch))  

       litt_c => currentPatch%litter(element_pos(carbon12_element))
       
       ! How much live grass is there? 
       currentPatch%livegrass = 0.0_r8 
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))
          if( int(prt_params%woody(currentCohort%pft)) == ifalse)then 
             
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
     

       if(write_sf == itrue)then
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter1 ',sum(litt_c%leaf_fines(:))
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter2 ',sum(litt_c%ag_cwd(:))
          if ( hlm_masterproc == itrue ) write(fates_log(),*) ' leaf_litter3 ',currentPatch%livegrass
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
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'ff2a ', &
                  lg_sf,currentPatch%livegrass,currentPatch%sum_fuel
          endif

          currentPatch%fuel_frac(lg_sf)       = currentPatch%livegrass       / currentPatch%sum_fuel
          
          ! MEF (moisure of extinction) depends on compactness of fuel, depth, particle size, wind, slope
          ! Eq here is Eq 27 from Peterson and Ryan (1986) "Modeling Postfire Conifer Mortality for Long-Range Planning"
          ! but lots of other approaches in use out there...
          ! MEF: pine needles=0.30 (text near Eq 28 Rothermal 1972)
          ! Table II-1 NFFL mixed fuels models from Rothermal 1983 Gen. Tech. Rep. INT-143 
          ! MEF: short grass=0.12,tall grass=0.25,chaparral=0.20,closed timber litter=0.30,hardwood litter=0.25
          ! Thonicke 2010 SAV values propagated thru P&R86 eqn below gives MEF:tw=0.355, sb=0.44, lb=0.525, tr=0.63, dg=0.248, lg=0.248
          ! Lasslop 2014 Table 1 MEF PFT level:grass=0.2,shrubs=0.3,TropEverGrnTree=0.2,TropDecid Tree=0.3, Extra-trop Tree=0.3
          MEF(1:nfsc)                         = 0.524_r8 - 0.066_r8 * log(SF_val_SAV(1:nfsc)) 

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

       endif
       ! check values. 
       ! FIX(SPM,032414) refactor...
       if(write_SF == itrue.and.currentPatch%fuel_sav <= 0.0_r8.or.currentPatch%fuel_bulkd <=  &
            0.0_r8.or.currentPatch%fuel_mef <= 0.0_r8.or.currentPatch%fuel_eff_moist <= 0.0_r8)then
            if ( hlm_masterproc == itrue ) write(fates_log(),*) 'problem with spitfire fuel averaging'
       endif 
       
       ! remove mineral content from net fuel load per Thonicke 2010
       ! for ir calculation in subr. rate_of_spread
       ! slevis moved here because rate_of_spread is now called twice/timestep
       currentPatch%sum_fuel  = currentPatch%sum_fuel * (1.0_r8 - SF_val_miner_total) !net of minerals

       currentPatch => currentPatch%younger

    enddo !end patch loop
    
  end subroutine characteristics_of_fuel

  !****************************************************************
  subroutine  characteristics_of_crown ( currentSite )
  !****************************************************************.  

    !returns the live crown fuel characteristics within each patch.
    ! passive_crown_FI is minimum fire intensity to ignite canopy crown fuel

    use SFParamsMod,    only : SF_VAL_CWD_FRAC

    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    ! ARGUMENTS

    ! LOCAL
    real(r8) ::  crown_depth          ! depth of crown (m)
    real(r8) ::  height_cbb           ! clear branch bole height or crown base height (m)
    real(r8) ::  max_height           ! max cohort on patch (m)
    real(r8) ::  crown_ignite_energy  ! heat yield for crown (kJ/kg)
    real(r8) ::  tree_sapw_struct_c   ! above-ground tree struct and sap biomass in cohort (kgC)
    real(r8) ::  leaf_c                  ! leaf carbon (kgC)
    real(r8) ::  sapw_c                  ! sapwood carbon (kgC)
    real(r8) ::  struct_c                ! structure carbon (kgC)
    real(r8) ::  twig_sapw_struct_c      ! above-ground twig sap and struct in cohort (kgC)
    real(r8) ::  crown_fuel_c            ! biomass of 1 hr fuels (leaves,twigs) in cohort (kg C)
    real(r8) ::  crown_fuel_biomass      ! biomass of crown fuel in cohort (kg biomass)
    real(r8) ::  crown_fuel_per_m        ! crown fuel per 1m section in cohort
    real(r8) ::  height_base_canopy      ! lowest height of fuels in patch to carry fire in crown

    integer  ::  ih                      ! counter

    real, dimension(70):: biom_matrix   ! matrix to track biomass from bottom to 70m
    real(r8),parameter :: min_density_canopy_fuel = 0.011_r8 !min canopy fuel density (kg/m3) sufficient to
                                                             !propogate fire vertically through canopy
                                                             !Scott and Reinhardt 2001 RMRS-RP-29
    real(r8),parameter :: foliar_moist_content = 1.0_r8      !foliar moisture content default 100% Scott & Reinhardt 2001


    !returns the live crown fuel characteristics within each patch.
    ! passive_crown_FI is the required minimum fire intensity to ignite canopy crown fuel

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       !zero Patch level variables
       height_base_canopy                   = 0.0_r8
       currentPatch%canopy_fuel_load        = 0.0_r8
       currentPatch%passive_crown_FI        = 0.0_r8
       currentPatch%canopy_bulk_density     = 0.0_r8
       max_height = 0._r8
       biom_matrix(:) = 0._r8

          currentCohort=>currentPatch%tallest
          do while(associated(currentCohort))

             !zero cohort level variables
             tree_sapw_struct_c                   = 0.0_r8
             leaf_c                               = 0.0_r8
             sapw_c                               = 0.0_r8
             struct_c                             = 0.0_r8
             twig_sapw_struct_c                   = 0.0_r8
             crown_fuel_c                         = 0.0_r8
             crown_fuel_biomass                   = 0.0_r8
             crown_fuel_per_m                     = 0.0_r8

             ! Calculate crown 1hr fuel biomass (leaf, twig sapwood, twig structural biomass)
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees

                call CrownDepth(currentCohort%hite,currentCohort%pft,crown_depth)
                height_cbb   = currentCohort%hite - crown_depth

                !find patch max height for stand canopy fuel
                if (currentCohort%hite > max_height) then
                   max_height = currentCohort%hite
                endif

                leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                sapw_c   = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)

                tree_sapw_struct_c =  currentCohort%n * &
                        (prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))

                twig_sapw_struct_c =  tree_sapw_struct_c * SF_VAL_CWD_frac(1)   !only 1hr fuel

                crown_fuel_c = (currentCohort%n * leaf_c) + twig_sapw_struct_c  !crown fuel (kgC)

                crown_fuel_biomass = crown_fuel_c / 0.45_r8            ! crown fuel (kg biomass)

                crown_fuel_per_m = crown_fuel_biomass / crown_depth    ! kg biomass per m

                !sort crown fuel into bins from bottom to top of crown
                !accumulate across cohorts to find density within canopy 1m sections
                do ih = max(1, nint(height_cbb)), max(1, nint(currentCohort%hite))
                   biom_matrix(ih) = biom_matrix(ih) + crown_fuel_per_m
                end do

                !accumulate available canopy fuel for patch (kg biomass)
                ! use this in CFB (crown fraction burn) calculation and FI final
                currentPatch%canopy_fuel_load = currentPatch%canopy_fuel_load + crown_fuel_biomass  !canopy fuel in patch

             endif !trees only

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop

          biom_matrix(:) = biom_matrix(:) / currentPatch%area    !kg biomass/m3

          !loop from 1m to 70m to find bin with total density = 0.011 kg/m3
          !min canopy fuel density to propogate fire vertically in canopy across patch
          do ih=1,70
             if (biom_matrix(ih) > min_density_canopy_fuel) then
                height_base_canopy = float(ih)
                exit
             end if
          end do

          !canopy_bulk_denisty (kg/m3) for Patch
          if (max_height - height_base_canopy > 0._r8) then
             currentPatch%canopy_bulk_density = sum(biom_matrix) / (max_height - height_base_canopy)
          else
             currentPatch%canopy_bulk_density = 0._r8
          end if

          ! Note: crown_ignition_energy to be calculated based on PFT foliar moisture content from FATES-Hydro
          ! or create foliar moisture % based on BTRAN
          ! Use foliar_moisture(currentCohort%pft) and compute weighted PFT average with Eq 3 Van Wagner 1977
          ! in place of foliar_moist_content parameter

          ! Eq 3 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
          ! h = 460.0 + 25.9*m
          ! h = crown_ignite_energy (kJ/kg), m = foliar moisture content based on dry fuel (%)
          crown_ignite_energy = 460.0 + 25.9 * foliar_moist_content

          ! Crown fuel ignition potential (kW/m), Eq 4 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
          ! FI = (Czh)**3/2 where z=canopy base height,h=heat of crown ignite energy, FI=fire intensity
          ! 0.01 = C, empirical constant Van Wagner 1977 Eq 4 for 6m canopy base height, 100% FMC, FI 2500kW/m
          ! passive_crown_FI = min fire intensity to ignite canopy fuel (kW/m or kJ/m/s)
          currentPatch%passive_crown_FI = (0.01_r8 * height_base_canopy * crown_ignite_energy)**1.5_r8
      
      currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine characteristics_of_crown

  !*****************************************************************
  subroutine  wind_effect ( currentSite, bc_in) 
  !*****************************************************************.

    ! Routine called daily from within ED within a site loop.
    ! Calculates the effective windspeed based on vegetation characteristics.
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
          if( int(prt_params%woody(currentCohort%pft)) == itrue)then
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

  !*******************************************************************
  subroutine rate_of_spread ( currentSite)
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

    ! ARGUMENTS

    ! LOCAL VARIABLES

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
    real(r8) time_r               ! residence time (min)
    real(r8) temp1, temp2

    real(r8),parameter :: q_dry = 581.0_r8             !heat of pre-ignition of dry fuels (kJ/kg)
    real(r8),parameter :: wind_reduce = 0.2_r8         !wind reduction factor (%)


    currentPatch=>currentSite%oldest_patch;  

    do while(associated(currentPatch))

!! clean up this initialise to zero section??
              
        ! ---initialise parameters to zero.--- 
       beta_ratio = 0.0_r8; q_ig = 0.0_r8; eps = 0.0_r8;   a = 0.0_r8;   b = 0.0_r8;   c = 0.0_r8;   e = 0.0_r8
       phi_wind = 0.0_r8;   xi = 0.0_r8;   reaction_v_max = 0.0_r8;  reaction_v_opt = 0.0_r8; mw_weight = 0.0_r8
       moist_damp = 0.0_r8;   ir = 0.0_r8; a_beta = 0.0_r8;     
       currentPatch%ROS_front = 0.0_r8

       ! ----start spreading---

       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) &
            'SF - currentPatch%fuel_bulkd ',currentPatch%fuel_bulkd
       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) &
            'SF - SF_val_part_dens ',SF_val_part_dens

       ! beta = packing ratio (unitless)
       ! fraction of fuel array volume occupied by fuel or compactness of fuel bed 
       beta = currentPatch%fuel_bulkd / SF_val_part_dens
       
       ! Eq A6 in Thonicke et al. 2010
       ! packing ratio (unitless) 
       beta_op = 0.200395_r8 *(currentPatch%fuel_sav**(-0.8189_r8))

       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta ',beta
       if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta_op ',beta_op
       beta_ratio = beta/beta_op   !unitless

       if(write_sf == itrue)then
          if ( hlm_masterproc == itrue ) write(fates_log(),*) 'esf ',currentPatch%fuel_eff_moist
       endif

       ! ---heat of pre-ignition---
       !  Eq A4 in Thonicke et al. 2010, Eq 12 Rothermel 1972
       !  50 Btu/lb + 1116 Btu/lb * fuel_eff_moist
       !  conversion of Rothermel (1972) Eq 12 in BTU/lb to current kJ/kg 
       !  q_ig in kJ/kg 
       q_ig = q_dry + 2594.0_r8 * currentPatch%fuel_eff_moist

       ! ---effective heating number---
       ! Eq A3 in Thonicke et al. 2010.  
       eps = exp(-4.528_r8 / currentPatch%fuel_sav)     
       ! Eq A7 in Thonicke et al. 2010 per Eq 49 Rothermel 1972
       b = 0.15988_r8 * (currentPatch%fuel_sav**0.54_r8)
       ! Eq A8 in Thonicke et al. 2010 per Eq 48 Rothermel 1972 
       c = 7.47_r8 * (exp(-0.8711_r8 * (currentPatch%fuel_sav**0.55_r8)))
       ! Eq A9 in Thonicke et al. 2010. (has typo, using coefficient Eq 50 Rothermel 1972)
       e = 0.715_r8 * (exp(-0.01094_r8 * currentPatch%fuel_sav))

       if (debug) then
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - c ',c
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - currentPatch%effect_wspeed ', &
                                                                         currentPatch%effect_wspeed
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - b ',b
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - beta_ratio ',beta_ratio
          if ( hlm_masterproc == itrue .and.debug) write(fates_log(),*) 'SF - e ',e
       endif

       ! Eq A5 in Thonicke et al. 2010
       ! phi_wind (unitless)
       ! convert current_wspeed (wind at elev relevant to fire) from m/min to ft/min for Rothermel ROS eqn
       phi_wind = c * ((3.281_r8*currentPatch%effect_wspeed)**b)*(beta_ratio**(-e))


       ! ---propagating flux----
       ! Eq A2 in Thonicke et al.2010 and Eq 42 Rothermel 1972
       ! xi (unitless)       
       xi = (exp((0.792_r8 + 3.7597_r8 * (currentPatch%fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
            (192_r8+7.9095_r8 * currentPatch%fuel_sav)      
      
       ! ---reaction intensity----
       ! Eq in table A1 Thonicke et al. 2010. 
       a = 8.9033_r8 * (currentPatch%fuel_sav**(-0.7913_r8))
       a_beta = exp(a*(1.0_r8-beta_ratio))  !dummy variable for reaction_v_opt equation
  
       ! Eq in table A1 Thonicke et al. 2010.
       ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
       ! reaction_v_max = Eq 36 in Rothermel 1972 and Fig 12 
       reaction_v_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (currentPatch%fuel_sav**(-1.5_r8)))
       ! reaction_v_opt =  Eq 38 in Rothermel 1972 and Fig 11
       reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta

       ! mw_weight = relative fuel moisture/fuel moisture of extinction
       ! average values for litter pools (dead leaves, twigs, small and large branches) plus grass
       mw_weight = currentPatch%fuel_eff_moist/currentPatch%fuel_mef
       
       ! Eq in table A1 Thonicke et al. 2010. 
       ! moist_damp is unitless
       moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
            (3.52_r8*(mw_weight**3.0_r8))))

       ! ir = reaction intenisty in kJ/m2/min
       ! currentPatch%sum_fuel converted from kgC/m2 to kgBiomass/m2 for ir calculation
       ir = reaction_v_opt*(currentPatch%sum_fuel/0.45_r8)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp 

       ! write(fates_log(),*) 'ir',gamma_aptr,moist_damp,SF_val_fuel_energy,SF_val_miner_damp

     
       if (currentPatch%fuel_bulkd <= 0.0_r8 .or. eps <= 0.0_r8 .or. &
           q_ig <= 0.0_r8 .or. ir <= 0.0_r8 .or. &
           currentPatch%fuel_sav <= 0.0_r8 .or. xi <= 0.0_r8 .or. &
           c <= 0.0_r8 .or. beta_ratio <= 0.0_r8 .or. b <= 0.0_r8) then
          currentPatch%ROS_front = 0.0_r8
!         currentPatch%ROS_torch = 0.0_r8  ! potentially useful diagnostic
       else ! Eq 9. Thonicke et al. 2010. 
            ! forward ROS in m/min
          currentPatch%ROS_front = (ir*xi*(1.0_r8+phi_wind)) / (currentPatch%fuel_bulkd*eps*q_ig)
          ! write(fates_log(),*) 'ROS',currentPatch%ROS_front,phi_wind,currentPatch%effect_wspeed
          ! write(fates_log(),*) 'ros calcs',currentPatch%fuel_bulkd,ir,xi,eps,q_ig

          ! calculate heat release per unit area (HPA)(kJ/m2), Eq 2 Scott & Reinhardt 2001
          ! and residence time (min), Eq 3 Scott & Reinhardt 2001
          time_r = 12.595_r8 / currentPatch%fuel_sav
          currentPatch%heat_per_area = ir * time_r

          ! calculate torching index based on wind speed and crown fuels 
          ! ROS for crown torch initation (m/min), Eq 18 Scott & Reinhardt 2001
          ! temp1 requires max(0... to avoid raising a negative value to a
          ! fractional temp2 power.
          ! TODO slevis: omitting "- phi_s" after the -1 bc phi_s = 0 for now.
          !              @jkshuman recommended naming it slope_factor.
!         temp1 = max(0._r8, (60._r8 * currentPatch%passive_crown_FI * &
! currentPatch%fuel_bulkd * eps * q_ig / (currentPatch%heat_per_area * ir * xi) - 1._r8) / &
! (c * beta_ratio**-e))
!         temp2 = 1._r8 / b
!         currentPatch%ROS_torch = temp1**temp2 / (54.683_r8 * wind_reduce)
       endif
       ! Eq 10 in Thonicke et al. 2010
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
    real(r8) :: fc_ground(nfsc) !total amount of fuel consumed per area of burned ground (kg C / m2 of burned area)

    integer  :: c

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch))
       currentPatch%burnt_frac_litter(:) = 1.0_r8       
       ! Calculate fraction of litter is burnt for all classes. 
       ! Eq B1 in Thonicke et al. 2010---
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
       currentPatch%burnt_frac_litter(:) = currentPatch%burnt_frac_litter(:) * (1.0_r8-SF_val_miner_total) 

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
  subroutine  area_burnt_intensity ( currentSite, bc_in)
  !*****************************************************************

    !returns the updated currentPatch%FI value for each patch.

    !currentPatch%FI  avg fire intensity of flaming front during day. Backward ROS plays no role here. kJ/m/s or kW/m.
    !currentSite%FDI  probability that an ignition will start a fire
    !currentSite%NF   number of lighting strikes per day per km2
    !currentPatch%ROS_front  forward ROS (m/min) 
    !currentPatch%TFC_ROS total fuel consumed by flaming front (kgC/m2 of burned area)

    use EDParamsMod,       only : ED_val_nignitions
    use EDParamsMod,       only : cg_strikes    ! fraction of cloud-to-ground ligtning strikes
    use FatesConstantsMod, only : years_per_day
    use SFParamsMod,       only : SF_val_fdi_alpha,SF_val_fuel_energy, &
                                  SF_val_max_durat, SF_val_durat_slope, SF_val_fire_threshold
    
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type), pointer :: currentPatch
    type(bc_in_type), intent(in) :: bc_in

    ! ARGUMENTS

    ! LOCAL VARIABLES
    real(r8) ROS                     !rate of spread (m/s)
    real(r8) W                       !available fuel (kgBiomass/m2)
    real(r8) :: tree_fraction_patch  !patch level. no units
    real(r8) df                      !distance fire has travelled forward (m)
    real(r8) db                      !distance fire has travelled backward (m)
    real(r8) AB                      !daily area burnt (m2 per km2)
    real(r8) size_of_fire            !in m2
    real(r8) cloud_to_ground_strikes ! [fraction] depends on hlm_spitfire_mode
    real(r8) anthro_ign_count        ! anthropogenic ignition count/km2/day
    integer :: iofp                  ! index of oldest fates patch
    real(r8), parameter :: pot_hmn_ign_counts_alpha = 0.0035_r8  ! Potential human ignition counts (alpha in Li et al. 2012) (#/person/month)
    real(r8), parameter :: km2_to_m2 = 1000000.0_r8              ! area conversion for square km to square m
    real(r8), parameter :: m_per_min__to__km_per_hour = 0.06_r8  ! convert wind speed from m/min to km/hr
    real(r8), parameter :: forest_grassland_lengthtobreadth_threshold = 0.55_r8 ! tree canopy cover below which to use 
                                                                                ! grassland length-to-breadth eqn
                                                                                ! 0.55 = benchmark forest cover, Staver 2010

    !  ---initialize site parameters to zero--- 
    currentSite%NF_successful = 0._r8
    
    ! Eq 7 from Venevsky et al GCB 2002 (modification of Eqn 8, Thonicke et al. 2010) 
    ! FDI 0.1 = low, 0.3 moderate, 0.75 high, and 1 = extreme ignition potential for alpha 0.000337
    if (hlm_spitfire_mode == hlm_sf_successful_ignitions_def) then
       currentSite%FDI = 1.0_r8           ! READING "SUCCESSFUL IGNITION" DATA
                                          ! force ignition potential to be extreme
       cloud_to_ground_strikes = 1.0_r8   ! cloud_to_ground = 1 = use 100% incoming observed ignitions
    else  ! USING LIGHTNING DATA
       currentSite%FDI  = 1.0_r8 - exp(-SF_val_fdi_alpha*currentSite%acc_NI)
       cloud_to_ground_strikes = cg_strikes
    end if
    
    !NF = number of lighting strikes per day per km2 scaled by cloud to ground strikes
    iofp = currentSite%oldest_patch%patchno
    if (hlm_spitfire_mode == hlm_sf_scalar_lightning_def ) then
       currentSite%NF = ED_val_nignitions * years_per_day * cloud_to_ground_strikes
    else    ! use external daily lightning ignition data
       currentSite%NF = bc_in%lightning24(iofp) * cloud_to_ground_strikes
    end if

    ! If there are 15  lightning strikes per year, per km2. (approx from NASA product for S.A.) 
    ! then there are 15 * 1/365 strikes/km2 each day 
 
    ! Calculate anthropogenic ignitions according to Li et al. (2012)
    ! Add to ignitions by lightning
    if (hlm_spitfire_mode == hlm_sf_anthro_ignitions_def) then
      ! anthropogenic ignitions (count/km2/day)
      !           =  ignitions/person/month * 6.8 * population_density **0.43 /approximate days per month
      anthro_ign_count = pot_hmn_ign_counts_alpha * 6.8_r8 * bc_in%pop_density(iofp)**0.43_r8 / 30._r8
                           
       currentSite%NF = currentSite%NF + anthro_ign_count

    end if

    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch))
       !  ---initialize patch parameters to zero---
       currentPatch%FI         = 0._r8
       currentPatch%fire       = 0
       currentPatch%FD         = 0.0_r8
       currentPatch%frac_burnt = 0.0_r8
       
       if (currentSite%NF > 0.0_r8) then
          
          ! Eq 14 in Thonicke et al. 2010
          ! fire duration in minutes
          currentPatch%FD = (SF_val_max_durat+1.0_r8) / (1.0_r8 + SF_val_max_durat * &
                            exp(SF_val_durat_slope*currentSite%FDI))
          if(write_SF == itrue)then
             if ( hlm_masterproc == itrue ) write(fates_log(),*) 'fire duration minutes',currentPatch%fd
          endif
          !Eq 15 in Arora and Boer CTEM model.Average fire is 1 day long.
          !currentPatch%FD = 60.0_r8 * 24.0_r8 !no minutes in a day

          tree_fraction_patch  = 0.0_r8
          tree_fraction_patch  = currentPatch%total_tree_area/currentPatch%area
       
          if(debug)then
             write(fates_log(),*) 'SF  currentPatch%area ',currentPatch%area
             write(fates_log(),*) 'SF  currentPatch%total_area ',currentPatch%total_tree_area
             write(fates_log(),*) 'SF  patch tree fraction ',tree_fraction_patch
             write(fates_log(),*) 'SF  AREA ',AREA
          endif         
 
          if ((currentPatch%effect_wspeed*m_per_min__to__km_per_hour) < 1._r8) then !16.67m/min = 1km/hr 
             currentPatch%lb = 1.0_r8
          else
             if (tree_fraction_patch > forest_grassland_lengthtobreadth_threshold) then 
                 ! Eq 79 forest fuels (Canadian Forest Fire Behavior Prediction System Ont.Inf.Rep. ST-X-3, 1992)
                 currentPatch%lb = (1.0_r8 + (8.729_r8 * &
  ((1.0_r8 -(exp(-0.03_r8 * m_per_min__to__km_per_hour * currentPatch%effect_wspeed)))**2.155_r8)))
             else ! Eq 80 grass fuels (CFFBPS Ont.Inf.Rep. ST-X-3, 1992, with correction from errata published in 
                  ! Inf.Rep. GLC-X-10 (Bottom et al., 2009) because of typo in CFFBPS Ont.Inf.Rep. ST-X-3, 1992)
                 currentPatch%lb = 1.1_r8 * &
  ((m_per_min__to__km_per_hour * currentPatch%effect_wspeed)**0.464_r8)
             endif
          endif

          !     if (lb > 8.0_r8)then
          !       lb = 8.0_r8  !Constraint Canadian Fire Behaviour System
          !     endif
          ! ---- calculate length of major axis---
          db = currentPatch%ROS_back  * currentPatch%FD !m
          df = currentPatch%ROS_front * currentPatch%FD !m

          ! --- calculate area burnt---
          if (currentPatch%lb > 0.0_r8) then

             ! Eq 1 in Thonicke et al. 2010
             ! To Do: Connect here with the Li & Levis GDP fire suppression algorithm. 
             ! Eq 16 in arora and boer model JGR 2005
             ! AB = AB *3.0_r8

             !size of fire = Eq 14 Arora and Boer JGR 2005 (area of an ellipse)
             size_of_fire = (df + db) * (df + db) * pi_const / (4.0_r8 * currentPatch%lb)

             ! AB = daily area burnt = size fires in m2 * num ignitions per day per km2 * prob ignition starts fire
             ! AB = m2 per km2 per day
             ! the denominator in the units of currentSite%NF is total gridcell area, but since we assume that ignitions 
             ! are equally probable across patches, currentSite%NF is equivalently per area of a given patch
             ! thus AB has units of m2 burned area per km2 patch area per day
             AB = size_of_fire * currentSite%NF * currentSite%FDI

             ! frac_burnt 
             ! just a unit conversion from AB, to become area burned per area patch per day, 
             ! or just the fraction of the patch burned on that day
             currentPatch%frac_burnt = (min(0.99_r8, AB / km2_to_m2))
             
             if(write_SF == itrue)then
                if ( hlm_masterproc == itrue ) write(fates_log(),*) 'frac_burnt',currentPatch%frac_burnt
             endif

          else
             currentPatch%frac_burnt = 0._r8
          endif ! lb

         ROS   = currentPatch%ROS_front / 60.0_r8 !m/min to m/sec for FI calculation
         W     = currentPatch%TFC_ROS / 0.45_r8   !kgC/m2 of burned area to kgbiomass/m2 of burned area

         ! Eq 15 Thonicke et al 2010
         !units of fire intensity = (kJ/kg)*(kgBiomass/m2)*(m/sec)
         currentPatch%FI = SF_val_fuel_energy * W * ROS  !kj/m/s, or kW/m
       
         if(write_sf == itrue)then
             if( hlm_masterproc == itrue ) write(fates_log(),*) 'fire_intensity',currentPatch%fi,W,currentPatch%ROS_front
         endif

         !'decide_fire' subroutine 
         if (currentPatch%FI > SF_val_fire_threshold) then !track fires greater than kW/m energy threshold
            currentPatch%fire = 1 ! Fire...    :D
            !
            currentSite%NF_successful = currentSite%NF_successful + &
                 currentSite%NF * currentSite%FDI * currentPatch%area / area
            !
         else     
            currentPatch%fire       = 0 ! No fire... :-/
            currentPatch%FD         = 0.0_r8
            currentPatch%frac_burnt = 0.0_r8
         endif         
          
       endif! NF ignitions check
       
       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine area_burnt_intensity

  !*****************************************************************
  subroutine  active_crown_fire ( currentSite)
  !*****************************************************************

    !evaluates if there will be an active crown fire based on canopy fuel and rate of spread
    !returns final rate of spread and fire intensity in patch with added fuel from active crown fire.
    !currentCohort%fraction_crown_burned is the proportion of crown affected by fire

    use EDParamsMod, only : active_crown_fire_switch
    use SFParamsMod, only  : SF_val_miner_total, SF_val_part_dens, SF_val_miner_damp, &
                             SF_val_fuel_energy, SF_val_drying_ratio

 
    type(ed_site_type), intent(in), target :: currentSite

    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    ! ARGUMENTS

    ! LOCAL VARIABLES
    ! Active crown Rothermel fire spread model parameters using FM 10
    real(r8) beta,beta_op         ! weighted average of packing ratio (unitless)
    real(r8) ir                   ! reaction intensity (kJ/m2/min)
    real(r8) xi,eps,phi_wind      ! all are unitless
    real(r8) q_ig                 ! heat of pre-ignition (kJ/kg)
    real(r8) reaction_v_opt,reaction_v_max !reaction velocity (per min)!optimum and maximum
    real(r8) moist_damp,mw_weight ! moisture dampening coefficient and ratio fuel moisture to extinction
    real(r8) beta_ratio           ! ratio of beta/beta_op
    real(r8) a_beta               ! dummy variable for product of a* beta_ratio for react_v_opt equation
    real(r8) a,b,c,e              ! function of fuel sav
    real(r8) total_fuel           ! total fuel (kg biomass/m2)
    real(r8) net_fuel             ! net fuel (kg biomass/m2) without minerals
    real(r8) fuel_depth           ! fuel depth (m)
    real(r8) fuel_bd              ! fuel bulk density (kg biomass/m3)
    real(r8) fuel_sav             ! fuels average sav 
    real(r8) fuel_eff_moist       ! fuels effective moisture
    real(r8) fuel_moist1hr        ! moisture 1 hour fuels
    real(r8) fuel_moist10hr       ! moisture 10 hour fuels
    real(r8) fuel_moist100hr      ! moisture 100 hour fuels
    real(r8) fuel_moistlive       ! moisture live fuels
    real(r8) SAV_1hr              ! surface area to volume 1 hour fuels (twigs)
    real(r8) SAV_10hr             ! surface area to volume 10 hour fuels (small branches)
    real(r8) SAV_100hr            ! surface area to volume 100 hour fuels (large branches)
    real(r8) SAV_live             ! surface area to volume live fuels
    real(r8) midflame_wind        ! 40% of open wind speed, Scott & Reinhardt 2001 
    real(r8) db                   ! distance fire has traveld backward (m)
    real(r8) df                   ! distance fire has travelled forward (m)
    real(r8) AB                   ! daily area burnt (m2 per km2)  
    real(r8) size_of_fire         ! in m2
    real(r8) ROS_active           ! actual rate of spread (m/min) using FM 10 fuels
    real(r8) ROS_active_min       ! minimum rate of spread to ignite active crown fire
    real(r8) CI_temp              ! temporary variable to calculate wind_active_min
    real(r8) wind_active_min      ! open windspeed to sustain active crown fire where ROS_SA = ROS_active_min
    real(r8) ROS_SA               ! rate of spread for surface fire with wind_active_min
    real(r8) canopy_frac_burnt    ! fraction of canopy fuels consumed (0, surface fire to 1,active crown fire) 
    real(r8) ROS_final            ! final rate of spread for combined surface and canopy spread (m/min)
    real(r8) FI_final             ! final fireline intensity (kW/m or kJ/m/sec) with canopy consumption 

    real(r8),parameter :: q_dry = 581.0_r8                 !heat of pre-ignition of dry fuels (kJ/kg)
    ! fuel loading, MEF, and depth from Anderson 1982 Aids to determining fuel models for fire behavior
    ! SAV values from BEHAVE model Burgan & Rothermel 1984)
    real(r8),parameter :: fuel_1hr     = 3.01_r8             ! FM 10 1-hr fuel loading (US tons/acre)
    real(r8),parameter :: fuel_10hr    = 2.0_r8              ! FM 10 10-hr fuel loading (US tons/acre)             
    real(r8),parameter :: fuel_100hr   = 5.01_r8             ! FM 10 100-hr fuel loading (US tons/acre)
    real(r8),parameter :: fuel_live    = 2.0_r8              ! FM 10 live fuel loading (US tons/acre)
    real(r8),parameter :: fuel_mef     = 0.25_r8             ! FM 10 moisture of extinction (volumetric)
    real(r8),parameter :: fuel_depth_ft= 1.0_r8              ! FM 10 fuel depth (ft)
    real(r8),parameter :: sav_1hr_ft   = 2000.0_r8           ! FM 10 1-hr SAV (ft2/ft3)
    real(r8),parameter :: sav_10hr_ft  = 109.0_r8            ! FM 10 10-hr SAV (ft2/ft3)             
    real(r8),parameter :: sav_100hr_ft = 30.0_r8             ! FM 10 100-hr SAV (ft2/ft3)
    real(r8),parameter :: sav_live_ft  = 1650.0_r8           ! FM 10 live SAV (ft2/ft3)
    real(r8),parameter :: tonnes_acre_to_kg_m2 = 0.2241701   ! convert tons/acre to kg/m2
    real(r8),parameter :: sqft_cubicft_to_sqm_cubicm = 0.03280844 !convert ft2/ft3 to m2/m3
    real(r8),parameter :: canopy_ignite_energy = 18000_r8    ! heat yield for canopy fuels (kJ/kg)
    real(r8),parameter :: critical_mass_flow_rate = 0.05_r8  ! critical mass flow rate (kg/m2/sec)for crown fire
    real(r8),parameter :: km2_to_m2 = 1000000.0_r8           ! area conversion for square km to square m

    integer  :: passive_canopy_fuel_flg                    ! flag if canopy fuel true for vertical spread


    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))

       if (currentPatch%fire == 1) then
          passive_canopy_fuel_flg = 0         !does patch have canopy fuels for vertical spread?
          ROS_active = 0.0_r8

          ! check initiation of passive crown fire
          if (currentPatch%FI >= currentPatch%passive_crown_FI) then
             passive_canopy_fuel_flg = 1      !enough passive canopy fuels for vertical spread

             ! Calculate rate of spread using FM 10 as in Rothermel 1977 
             ! fuel characteristics 
             total_fuel = (fuel_1hr + fuel_10hr + fuel_100hr + fuel_live) * tonnes_acre_to_kg_m2 

             SAV_1hr   = sav_1hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_10hr  = sav_10hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_100hr = sav_100hr_ft * sqft_cubicft_to_sqm_cubicm
             SAV_live  = sav_live_ft * sqft_cubicft_to_sqm_cubicm

             fuel_moist1hr    = exp(-1.0_r8 * ((SAV_1hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moist10hr   = exp(-1.0_r8 * ((SAV_10hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moist100hr  = exp(-1.0_r8 * ((SAV_100hr/SF_val_drying_ratio) * currentSite%acc_NI))
             fuel_moistlive   = exp(-1.0_r8 * ((SAV_live/SF_val_drying_ratio) * currentSite%acc_NI))

             fuel_depth       = fuel_depth_ft *0.3048           !convert to meters
             fuel_bd          = total_fuel/fuel_depth

             fuel_sav         = SAV_1hr *(fuel_1hr/total_fuel) + SAV_10hr*(fuel_10hr/total_fuel) + & 
                                 SAV_100hr*(fuel_100hr/total_fuel) + SAV_live*(fuel_live/total_fuel)

             fuel_eff_moist = fuel_moist1hr *(fuel_1hr/total_fuel) + fuel_moist10hr*(fuel_10hr/total_fuel) + & 
                               fuel_moist100hr*(fuel_100hr/total_fuel) + fuel_moistlive*(fuel_live/total_fuel)

             ! remove mineral content from net fuel load
             net_fuel = total_fuel * (1.0_r8 - SF_val_miner_total) !net of minerals

             ! ---start spreading---
             !beta = packing ratio (unitless)
             beta = fuel_bd / SF_val_part_dens
             beta_op = 0.200395_r8 *(fuel_sav**(-0.8189_r8))
             beta_ratio = beta/beta_op  

             ! -- heat of pre-ignition --
             q_ig = q_dry + 2594.0_r8 * fuel_eff_moist

             ! ---effective heating number---
             ! Eq A3 in Thonicke et al. 2010.  
             eps = exp(-4.528_r8 / fuel_sav)     
             ! Eq A7 in Thonicke et al. 2010 per Eq 49, Rothermel 1972
             b = 0.15988_r8 * (fuel_sav**0.54_r8)
             ! Eq A8 in Thonicke et al. 2010 per Eq 48, Rothermel 1972 
             c = 7.47_r8 * (exp(-0.8711_r8 * (fuel_sav**0.55_r8))) 
             ! Eq A9 in Thonicke et al. 2010. (typo in Eq A9, using coefficient Eq 50, Rothermel 1972)
             e = 0.715_r8 * (exp(-0.01094_r8 * fuel_sav))

             midflame_wind = currentSite%wind *0.40_r8  !Scott & Reinhardt 2001 40% open wind speed

             ! Eq A5 in Thonicke et al. 2010
             ! convert current_wspeed (wind at elev relevant to fire) from m/min to ft/min for Rothermel ROS eqn
             phi_wind = c * ((3.281_r8*midflame_wind)**b)*(beta_ratio**(-e)) !unitless

             ! ---propagating flux = xi (unitless) 
             ! Eq A2 in Thonicke et al.2010 and Eq 42 Rothermel 1972
             xi = (exp((0.792_r8 + 3.7597_r8 * (fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
                (192_r8+7.9095_r8 * fuel_sav) 

             ! ---reaction intensity----
             ! Eq in table A1 Thonicke et al. 2010. 
             a = 8.9033_r8 * (fuel_sav**(-0.7913_r8))
             a_beta = exp(a*(1.0_r8-beta_ratio))  !dummy variable for reaction_v_opt equation
  
             ! Eq in table A1 Thonicke et al. 2010.
             ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
             ! reaction_v_max = Eq 36 in Rothermel 1972 and Fig 12 
             reaction_v_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (fuel_sav**(-1.5_r8)))
             ! reaction_v_opt =  Eq 38 in Rothermel 1972 and Fig 11
             reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta

             ! mw_weight = relative fuel moisture/fuel moisture of extinction
             mw_weight = fuel_eff_moist/fuel_mef
       
             ! Eq in table A1 Thonicke et al. 2010. (unitless)
             moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
               (3.52_r8*(mw_weight**3.0_r8))))

             ! ir = reaction intenisty in kJ/m2/min
             ! sum_fuel as kgBiomass/m2 for ir calculation
             ir = reaction_v_opt*(net_fuel)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp  
 
             ! actual ROS (m/min) for FM 10 fuels for open windspeed, Eq 8 Scott & Reinhardt 2001
             ROS_active = 3.34_r8 * ((ir*xi*(1.0_r8+phi_wind)) / (fuel_bd * eps * q_ig))

             ! critical min rate of spread (m/min) for active crowning
             ROS_active_min = (critical_mass_flow_rate / fuel_bd) * 60.0_r8

             ! check threshold intensity and rate of spread
             if (active_crown_fire_switch == 1 .and. &
                 ROS_active >= ROS_active_min) then
                currentPatch%active_crown_fire_flg = 1  ! active crown fire ignited
                !ROS_final = ROS_surface+CFB(ROS_active - ROS_surface), Eq 21 Scott & Reinhardt 2001
                !with active crown fire CFB (canopy fraction burned) = 100%
                canopy_frac_burnt = 1.0_r8

             else 
                currentPatch%active_crown_fire_flg = 0  ! only passive crown fire with partial crown burnt

                ! phi_slope is not used yet. consider adding with later
                ! development
                ! calculate open wind speed critical to sustain active crown
                ! fire Eq 20 Scott & Reinhardt
                if (ir > 0._r8 .and. currentPatch%canopy_bulk_density > 0._r8) then
                   CI_temp = ((164.8_r8 * eps * q_ig)/(ir * currentPatch%canopy_bulk_density)) - 1.0_r8
                else
                   CI_temp = 0._r8
                end if

                wind_active_min = 0.0457_r8 * (CI_temp / 0.001612_r8)**0.7_r8

                ! use open wind speed "wind_active_min" for ROS surface fire
                ! where ROS_SA=ROS_active_min
                ROS_SA =  (ir * xi * (1.0_r8 + wind_active_min)) / (fuel_bd * eps * q_ig)

                ! canopy fraction burnt, Eq 28 Scott & Reinhardt Appendix A
                canopy_frac_burnt = max(0._r8, min(1.0_r8, &
                   (currentPatch%ROS_front - ROS_active_min) / (ROS_SA - ROS_active_min)))
                
             endif !check intensity & ROS for active crown fire thresholds

             !ROS_final = ROS_surface+CFB(ROS_active - ROS_surface), Eq 21 Scott & Reinhardt 2001
             ROS_final = currentPatch%ROS_front + &
                canopy_frac_burnt * (ROS_active - currentPatch%ROS_front)

             ! recalculate area burned with new ROS_front value from ROS_final
             ! ---- re-calculate length of major axis for df using new ROS_front value from ROS final---
             db = currentPatch%ROS_back  * currentPatch%FD !(m) 
             df = ROS_final * currentPatch%FD              !(m) update with ROS final 
             
             ! update ROS_front with ROS_final for output variable 
             currentPatch%ROS_front = ROS_final

             ! --- calculate updated area burnt using df from ROS final---
             if (currentPatch%lb > 0.0_r8) then

                ! Eq 1 in Thonicke et al. 2010
                ! To Do: Connect here with the Li & Levis GDP fire suppression algorithm. 
                ! Eq 16 in arora and boer model JGR 2005
                ! AB = AB *3.0_r8

                !size of fire = Eq 14 Arora and Boer JGR 2005 (area of an ellipse)
                size_of_fire = (df + db) * (df + db) * pi_const / (4.0_r8 * currentPatch%lb)

                ! AB = daily area burnt = size fires in m2 * num ignitions per day per km2 * prob ignition starts fire
                ! AB = m2 per km2 per day
                ! the denominator in the units of currentSite%NF is total gridcell area, but since we assume that ignitions 
                ! are equally probable across patches, currentSite%NF is equivalently per area of a given patch
                ! thus AB has units of m2 burned area per km2 patch area per day
                AB = size_of_fire * currentSite%NF * currentSite%FDI

                ! frac_burnt 
                ! just a unit conversion from AB, to become area burned per area patch per day, 
                ! or just the fraction of the patch burned on that day
                currentPatch%frac_burnt = (min(0.99_r8, AB / km2_to_m2))
             
                if(write_SF == itrue)then
                   if ( hlm_masterproc == itrue ) write(fates_log(),*) 'frac_burnt',currentPatch%frac_burnt
                endif

             else  
                currentPatch%frac_burnt = 0.0_r8
             endif ! lb

             !final fireline intensity (kJ/m/sec or kW/m), Eq 22 Scott & Reinhardt 2001
             FI_final = (currentPatch%heat_per_area + &
  currentPatch%canopy_fuel_load * canopy_ignite_energy * canopy_frac_burnt) * &
  currentPatch%ROS_front / 60._r8
             ! update patch FI to adjust according to potential canopy fuel consumed (passive and active)
             currentPatch%FI = FI_final

           endif !check if passive crown fire?
       endif !fire?

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine active_crown_fire


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
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees only

                leaf_c = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                sapw_c = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
                
                tree_ag_biomass = tree_ag_biomass + &
                      currentCohort%n * (leaf_c + & 
                      prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))
             endif !trees only
             currentCohort=>currentCohort%shorter;
          enddo !end cohort loop

          do i_pft=1,numpft
             if (tree_ag_biomass > 0.0_r8  .and. int(prt_params%woody(i_pft)) == itrue) then 
                
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

    real(r8) ::  crown_depth          ! depth of crown (m)
    real(r8) ::  height_cbb           ! clear branch bole height or crown base height (m) for cohort

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       !zero Patch level variables

       if (currentPatch%fire == 1) then

          currentCohort=>currentPatch%tallest

          do while(associated(currentCohort))  
             currentCohort%fraction_crown_burned = 0.0_r8
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees

                ! height_cbb = clear branch bole height at base of crown (m)
                ! inst%crown = crown_depth_frac (PFT)
                call CrownDepth(currentCohort%hite,currentCohort%pft,crown_depth)
                height_cbb   = currentCohort%hite - crown_depth

                ! Equation 17 in Thonicke et al. 2010
                ! flames over bottom of canopy, and potentially over top of
                ! canopy
                if (currentCohort%hite > 0.0_r8 .and. &
                    currentPatch%Scorch_ht(currentCohort%pft) >= height_cbb) then
                   if (currentPatch%active_crown_fire_flg == 0) then
                      currentCohort%fraction_crown_burned = min(1.0_r8, &
                      ((currentPatch%Scorch_ht(currentCohort%pft) - height_cbb) / crown_depth))
                   else  ! active crown fire occurring
                      currentCohort%fraction_crown_burned = 1.0_r8
                   end if
                endif  !SH frac crown burnt calculation
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
             currentCohort%cambial_mort = 0.0_r8
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees only
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
             if ( int(prt_params%woody(currentCohort%pft)) == itrue) then
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
