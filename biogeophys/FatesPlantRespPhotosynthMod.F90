module FATESPlantRespPhotosynthMod
   
   !-------------------------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Calculates the plant respiration and photosynthetic fluxes for the FATES model
   ! This code is similar to and was originally based off of the 'photosynthesis' 
   ! subroutine in the CLM model.
   !
   ! Parameter for activation and deactivation energies were taken from:
   ! Activation energy, from:
   ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
   ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
   ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282
   ! High temperature deactivation, from:
   ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
   ! The factor "c" scales the deactivation to a value of 1.0 at 25C
   ! Photosynthesis and stomatal conductance parameters, from:
   ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
   ! ------------------------------------------------------------------------------------
   
   ! !USES:

   use FatesGlobals, only      : endrun => fates_endrun
   use FatesGlobals, only      : fates_log
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : itrue
   use FatesInterfaceMod, only : hlm_use_planthydro
   use FatesInterfaceMod, only : numpft
   use EDTypesMod, only        : maxpft
   use EDTypesMod, only        : nlevleaf
   use EDTypesMod, only        : nclmax
   
   ! CIME Globals
   use shr_log_mod , only      : errMsg => shr_log_errMsg

   implicit none
   private
   
   public :: FatesPlantRespPhotosynthDrive ! Called by the HLM-Fates interface
   
   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   !-------------------------------------------------------------------------------------
   
   ! maximum stomatal resistance [s/m] (used across several procedures)
   real(r8),parameter :: rsmax0 =  2.e4_r8                    
   
   logical   ::  DEBUG = .false.

contains
  
  !--------------------------------------------------------------------------------------
  
  subroutine FatesPlantRespPhotosynthDrive (nsites, sites,bc_in,bc_out,dtime)

    ! -----------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    ! -----------------------------------------------------------------------------------


    ! !USES:

    use EDPftvarcon         , only : EDPftvarcon_inst 

    use FatesSynchronizedParamsMod , only : FatesSynchronizedParamsInst
    use EDTypesMod        , only : ed_patch_type
    use EDTypesMod        , only : ed_cohort_type
    use EDTypesMod        , only : ed_site_type
    use EDTypesMod        , only : maxpft
    use FatesInterfaceMod , only : hlm_numlevsoil
    use FatesInterfaceMod , only : bc_in_type
    use FatesInterfaceMod , only : bc_out_type
    use EDCanopyStructureMod, only : calc_areaindex
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use FatesConstantsMod, only : umol_per_mmol
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesParameterDerivedMod, only : param_derived
    use EDPatchDynamicsMod, only: set_root_fraction
    use EDParamsMod, only : ED_val_bbopt_c3, ED_val_bbopt_c4, ED_val_base_mr_20
    use FatesAllometryMod, only : bleaf, storage_fraction_of_target


    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    real(r8),intent(in)                     :: dtime


    ! LOCAL VARIABLES:
    ! -----------------------------------------------------------------------------------
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort

    ! -----------------------------------------------------------------------------------
    ! These three arrays hold leaf-level biophysical rates that are calculated
    ! in one loop and then sent to the cohorts in another loop.  If hydraulics are
    ! on, we calculate a unique solution for each level-cohort-layer combination.
    ! If we are not using hydraulics, we calculate a unique solution for each
    ! level-pft-layer combination.  Thus the following three arrays are statically 
    ! allocated for the maximum space of the two cases (numCohortsPerPatch)
    ! The "_z" suffix indicates these variables are discretized at the "leaf_layer"
    ! scale.
    ! Note: For these temporary arrays, we have the leaf layer dimension first
    ! and the canopy layer last. This order is chosen for efficiency. The arrays
    ! such as leaf area that are bound to the patch structure DO NOT follow this order
    ! as they are used in many other parts of the code with different looping, we
    ! are not modifying its order now.
    ! -----------------------------------------------------------------------------------

    ! leaf maintenance (dark) respiration [umol CO2/m**2/s]
    real(r8) :: lmr_z(nlevleaf,maxpft,nclmax)

    ! stomatal resistance [s/m]
    real(r8) :: rs_z(nlevleaf,maxpft,nclmax)    

    ! net leaf photosynthesis averaged over sun and shade leaves. [umol CO2/m**2/s]
    real(r8) :: anet_av_z(nlevleaf,maxpft,nclmax)  
    
    ! Mask used to determine which leaf-layer biophysical rates have been
    ! used already
    logical :: rate_mask_z(nlevleaf,maxpft,nclmax)

    real(r8) :: vcmax_z            ! leaf layer maximum rate of carboxylation 
                                   ! (umol co2/m**2/s)
    real(r8) :: jmax_z             ! leaf layer maximum electron transport rate 
                                   ! (umol electrons/m**2/s)
    real(r8) :: tpu_z              ! leaf layer triose phosphate utilization rate 
                                   ! (umol CO2/m**2/s)
    real(r8) :: kp_z               ! leaf layer initial slope of CO2 response 
                                   ! curve (C4 plants)
   
    real(r8) :: mm_kco2            ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2             ! Michaelis-Menten constant for O2 (Pa)
    real(r8) :: co2_cpoint         ! CO2 compensation point (Pa)
    real(r8) :: btran_eff          ! effective transpiration wetness factor (0 to 1) 
    real(r8) :: bbb                ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8) :: kn(maxpft)         ! leaf nitrogen decay coefficient
    real(r8) :: cf                 ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
    real(r8) :: gb_mol             ! leaf boundary layer conductance (molar form: [umol /m**2/s])
    real(r8) :: ceair              ! vapor pressure of air, constrained (Pa)
    real(r8) :: nscaler            ! leaf nitrogen scaling coefficient
    real(r8) :: leaf_frac          ! ratio of to leaf biomass to total alive biomass
    real(r8) :: laican             ! canopy sum of lai_z
    real(r8) :: vai                ! leaf and steam area in ths layer. 
    real(r8) :: tcsoi              ! Temperature response function for root respiration. 
    real(r8) :: tcwood             ! Temperature response function for wood
    
    real(r8) :: elai               ! exposed LAI (patch scale)
    real(r8) :: live_stem_n        ! Live stem (above-ground sapwood) 
                                   ! nitrogen content (kgN/plant)
    real(r8) :: live_croot_n       ! Live coarse root (below-ground sapwood) 
                                   ! nitrogen content (kgN/plant)
    real(r8) :: froot_n            ! Fine root nitrogen content (kgN/plant)
    real(r8) :: g_sb_leaves        ! Mean combined (stomata+boundary layer) leaf conductance [m/s]
                                   ! over all of the patch's leaves.  The "sb" refers to the combined
                                   ! "s"tomatal and "b"oundary layer.
                                   ! This quantity is relevant on leaf surfaces. It does not
                                   ! have units of /m2 leaf per say, but is implicitly on leaf surfaces
    real(r8) :: r_sb_leaves        ! Mean leaf resistance over all the patch's leaves [s/m]
                                   ! This is the direct reciprocal of g_sb_leaves
    real(r8) :: r_stomata          ! Mean stomatal resistance across all leaves in the patch [s/m]


    real(r8) :: maintresp_reduction_factor  ! factor by which to reduce maintenance respiration when storage pools are low
    real(r8) :: b_leaf             ! leaf biomass kgC
    real(r8) :: frac               ! storage pool as a fraction of target leaf biomass
    real(r8) :: check_elai         ! This is a check on the effective LAI that is calculated
                                   ! over each cohort x layer.
    real(r8) :: cohort_eleaf_area  ! This is the effective leaf area [m2] reported by each cohort
    
    ! -----------------------------------------------------------------------------------
    ! Keeping these two definitions in case they need to be added later
    !
    ! -----------------------------------------------------------------------------------
    !real(r8) :: psncanopy_pa  ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    !real(r8) :: lmrcanopy_pa  ! patch sunlit leaf maintenance respiration rate (umol CO2/m**2/s) 

    integer  :: cl,s,iv,j,ps,ft,ifp ! indices
    integer  :: nv                  ! number of leaf layers
    integer  :: NCL_p               ! number of canopy layers in patch

    ! Parameters
    ! -----------------------------------------------------------------------
    ! Base maintenance respiration rate for plant tissues base_mr_20
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    !
    ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
    ! (gC/gN/s)
    ! ------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Photosynthesis and stomatal conductance parameters, from:
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    ! -----------------------------------------------------------------------------------

    ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    ! For C3 and C4 plants
    ! -----------------------------------------------------------------------------------
    real(r8), dimension(0:1) :: bbbopt 

    associate(  &
         c3psn     => EDPftvarcon_inst%c3psn  , &
         slatop    => EDPftvarcon_inst%slatop , & ! specific leaf area at top of canopy, 
                                                  ! projected area basis [m^2/gC]
         woody     => EDPftvarcon_inst%woody  , & ! Is vegetation woody or not? 
         leafcn    => EDPftvarcon_inst%leafcn , & ! leaf C:N (gC/gN)
         frootcn   => EDPftvarcon_inst%frootcn, & ! froot C:N (gc/gN)   ! slope of BB relationship
         q10       => FatesSynchronizedParamsInst%Q10 )

      bbbopt(0) = ED_val_bbopt_c4
      bbbopt(1) = ED_val_bbopt_c3
      
      do s = 1,nsites

         ! Multi-layer parameters scaled by leaf nitrogen profile.
         ! Loop through each canopy layer to calculate nitrogen profile using
         ! cumulative lai at the midpoint of the layer
         
         ifp = 0
         currentpatch => sites(s)%oldest_patch
         do while (associated(currentpatch))  

            ifp   = ifp+1
            NCL_p = currentPatch%NCL_p
            
            ! Part I. Zero output boundary conditions
            ! ---------------------------------------------------------------------------
            bc_out(s)%rssun_pa(ifp)     = 0._r8
            bc_out(s)%rssha_pa(ifp)     = 0._r8

            g_sb_leaves = 0._r8
            check_elai  = 0._r8

            !psncanopy_pa = 0._r8
            !lmrcanopy_pa = 0._r8

            ! Part II. Filter out patches 
            ! Patch level filter flag for photosynthesis calculations
            ! has a short memory, flags:
            ! 1 = patch has not been called
            ! 2 = patch is currently marked for photosynthesis
            ! 3 = patch has been called for photosynthesis already
            ! ---------------------------------------------------------------------------
            if(bc_in(s)%filter_photo_pa(ifp)==2)then


               ! Part III. Calculate the number of sublayers for each pft and layer.  
               ! And then identify which layer/pft combinations have things in them.  
               ! Output:
               ! currentPatch%ncan(:,:)
               ! currentPatch%canopy_mask(:,:)
               call UpdateCanopyNCanNRadPresent(currentPatch)


               ! Part IV.  Identify some environmentally derived parameters:
               !           These quantities are biologically irrelevant
               !  Michaelis-Menten constant for CO2 (Pa)
               !  Michaelis-Menten constant for O2 (Pa)
               !  CO2 compensation point (Pa)
               !  leaf boundary layer conductance of h20
               !  constrained vapor pressure
               call GetCanopyGasParameters(bc_in(s)%forc_pbot,       & ! in
                                           bc_in(s)%oair_pa(ifp),    & ! in
                                           bc_in(s)%t_veg_pa(ifp),   & ! in
                                           bc_in(s)%tgcm_pa(ifp),    & ! in
                                           bc_in(s)%eair_pa(ifp),    & ! in
                                           bc_in(s)%esat_tv_pa(ifp), & ! in
                                           bc_in(s)%rb_pa(ifp),      & ! in
                                           mm_kco2,                  & ! out              
                                           mm_ko2,                   & ! out
                                           co2_cpoint,               & ! out
                                           cf,                       & ! out
                                           gb_mol,                   & ! out
                                           ceair)                      ! out

               ! Part V.  Pre-process some variables that are PFT dependent
               ! but not environmentally dependent
               ! ------------------------------------------------------------------------

               do ft = 1,numpft

                  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
                  ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al 
                  ! (2010) Biogeosciences, 7, 1833-1859
                  ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
                  
                  if (bc_in(s)%dayl_factor_pa(ifp)  ==  0._r8) then
                     kn(ft) =  0._r8
                  else
                     kn(ft) = exp(0.00963_r8 * EDPftvarcon_inst%vcmax25top(ft) - 2.43_r8)
                  end if
                  
               end do !ft 

               call set_root_fraction(currentPatch,bc_in(s)%zi_sisl)

               ! ------------------------------------------------------------------------
               ! Part VI: Loop over all leaf layers.
               ! The concept of leaf layers is a result of the radiative transfer scheme.
               ! A leaf layer has uniform radiation environment.  Leaf layers are a group
               ! of vegetation surfaces (stems and leaves) which inhabit the same 
               ! canopy-layer "CL", have the same functional type "ft" and within those
               ! two partitions are further partitioned into vertical layers where
               ! downwelling radiation attenuates in order.
               ! In this phase we loop over the leaf layers and calculate the
               ! photosynthesis and respiration of the layer (since all biophysical
               ! properties are homogeneous).  After this step, we can loop through
               ! our cohort list, associate each cohort with its list of leaf-layers
               ! and transfer these quantities to the cohort.
               ! With plant hydraulics, we must realize that photosynthesis and
               ! respiration will be different for leaves of each cohort in the leaf
               ! layers, as they will have there own hydraulic limitations.
               ! NOTE: Only need to flush mask on the number of used pfts, not the whole 
               ! scratch space.
               ! ------------------------------------------------------------------------
               rate_mask_z(:,1:numpft,:) = .false.

               if(currentPatch%countcohorts > 0.0)then   ! Ignore empty patches

                  currentCohort => currentPatch%tallest
                  do while (associated(currentCohort)) ! Cohort loop
                     
                     ! Identify the canopy layer (cl), functional type (ft)
                     ! and the leaf layer (IV) for this cohort
                     ft = currentCohort%pft
                     cl = currentCohort%canopy_layer
                     
                     call bleaf(currentCohort%dbh,currentCohort%pft,currentCohort%canopy_trim,b_leaf)
                     call storage_fraction_of_target(b_leaf, currentCohort%bstore, frac)
                     call lowstorage_maintresp_reduction(frac,currentCohort%pft, &
                          maintresp_reduction_factor)

                     ! are there any leaves of this pft in this layer?
                     if(currentPatch%canopy_mask(cl,ft) == 1)then 
                        
                        if(cl==1)then !are we in the top canopy layer or a shaded layer?
                           laican = 0._r8
                        else

                           laican = sum(currentPatch%canopy_layer_tai(1:cl-1)) 

                        end if
                        
                        ! Loop over leaf-layers
                        do iv = 1,currentCohort%nv
                           
                           ! ------------------------------------------------------------
                           ! If we are doing plant hydro-dynamics (or any run-type
                           ! where cohorts may generate different photosynthetic rates
                           ! of other cohorts in the same canopy-pft-layer combo), 
                           ! we re-calculate the leaf biophysical rates for the 
                           ! cohort-layer combo of interest.  
                           ! but in the vanilla case, we only re-calculate if it has
                           ! not been done yet.
                           ! ------------------------------------------------------------
                           
                           if ( .not.rate_mask_z(iv,ft,cl) .or. (hlm_use_planthydro.eq.itrue) ) then
                              
                              if (hlm_use_planthydro.eq.itrue) then
!                                 write(fates_log(),*) 'hlm_use_planthydro'
!                                 write(fates_log(),*) 'has been set to true.  You have inadvertently'
!                                 write(fates_log(),*) 'turned on a future feature that is not in the'
!                                 write(fates_log(),*) 'FATES codeset yet. Please set this to'
!                                 write(fates_log(),*) 'false and re-compile.'
!                                 call endrun(msg=errMsg(sourcefile, __LINE__))

                                 bbb   = max (bbbopt(nint(c3psn(ft)))*currentCohort%co_hydr%btran(1), 1._r8)
                                 btran_eff = currentCohort%co_hydr%btran(1) 
                              else
                                 bbb   = max (bbbopt(nint(c3psn(ft)))*currentPatch%btran_ft(ft), 1._r8)
                                 btran_eff = currentPatch%btran_ft(ft)
                              end if
                              
                              ! Vegetation area index
                              vai = (currentPatch%elai_profile(cl,ft,iv)+currentPatch%esai_profile(cl,ft,iv)) 
                              if (iv == 1) then
                                 laican = laican + 0.5_r8 * vai
                              else
                                 laican = laican + 0.5_r8 * (currentPatch%elai_profile(cl,ft,iv-1)+ &
                                       currentPatch%esai_profile(cl,ft,iv-1)+vai)
                              end if
                              
                              ! Scale for leaf nitrogen profile
                              nscaler = exp(-kn(ft) * laican)

                              ! Part VII: Calculate dark respiration (leaf maintenance) for this layer
                              call LeafLayerMaintenanceRespiration( param_derived%lmr25top(ft),&  ! in
                                                                    nscaler,                  &  ! in
                                                                    ft,                       &  ! in
                                                                    bc_in(s)%t_veg_pa(ifp),   &  ! in
                                                                    lmr_z(iv,ft,cl))             ! out
                              
                              ! Part VII: Calculate (1) maximum rate of carboxylation (vcmax), 
                              ! (2) maximum electron transport rate, (3) triose phosphate 
                              ! utilization rate and (4) the initial slope of CO2 response curve 
                              ! (C4 plants). Earlier we calculated their base rates as dictated
                              ! by their plant functional type and some simple scaling rules for
                              ! nitrogen limitation baesd on canopy position (not prognostic).
                              ! These rates are the specific rates used in the actual photosynthesis
                              ! calculations that take localized environmental effects (temperature)
                              ! into consideration.

                              call LeafLayerBiophysicalRates(currentPatch%ed_parsun_z(cl,ft,iv), &  ! in
                                                             ft,                                 &  ! in
                                                             EDPftvarcon_inst%vcmax25top(ft),    &  ! in
                                                             param_derived%jmax25top(ft),        &  ! in
                                                             param_derived%tpu25top(ft),         &  ! in
                                                             param_derived%kp25top(ft),          &  ! in
                                                             nscaler,                            &  ! in
                                                             bc_in(s)%t_veg_pa(ifp),             &  ! in
                                                             btran_eff,                          &  ! in
                                                             vcmax_z,                            &  ! out
                                                             jmax_z,                             &  ! out
                                                             tpu_z,                              &  ! out
                                                             kp_z )                                 ! out
                              
                              ! Part IX: This call calculates the actual photosynthesis for the 
                              ! leaf layer, as well as the stomatal resistance and the net assimilated carbon.

                              call LeafLayerPhotosynthesis(currentPatch%f_sun(cl,ft,iv),    &  ! in
                                                        currentPatch%ed_parsun_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_parsha_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_laisun_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_laisha_z(cl,ft,iv), &  ! in
                                                currentPatch%canopy_area_profile(cl,ft,iv), &  ! in
                                                        ft,                                 &  ! in
                                                        vcmax_z,                            &  ! in
                                                        jmax_z,                             &  ! in
                                                        tpu_z,                              &  ! in
                                                        kp_z,                               &  ! in
                                                        bc_in(s)%t_veg_pa(ifp),             &  ! in
                                                        bc_in(s)%esat_tv_pa(ifp),           &  ! in
                                                        bc_in(s)%forc_pbot,                 &  ! in
                                                        bc_in(s)%cair_pa(ifp),              &  ! in
                                                        bc_in(s)%oair_pa(ifp),              &  ! in
                                                        btran_eff,                          &  ! in
                                                        bbb,                                &  ! in
                                                        cf,                                 &  ! in
                                                        gb_mol,                             &  ! in
                                                        ceair,                              &  ! in
                                                        mm_kco2,                            &  ! in
                                                        mm_ko2,                             &  ! in
                                                        co2_cpoint,                         &  ! in
                                                        lmr_z(iv,ft,cl),                    &  ! in
                                                        currentPatch%psn_z(cl,ft,iv),       &  ! out
                                                        rs_z(iv,ft,cl),                     &  ! out
                                                        anet_av_z(iv,ft,cl))                   ! out

                              rate_mask_z(iv,ft,cl) = .true.
                           end if
                        end do

                        ! Zero cohort flux accumulators.
                        currentCohort%npp_tstep  = 0.0_r8
                        currentCohort%resp_tstep = 0.0_r8
                        currentCohort%gpp_tstep  = 0.0_r8
                        currentCohort%rdark      = 0.0_r8
                        currentCohort%resp_m     = 0.0_r8
                        currentCohort%ts_net_uptake = 0.0_r8

                        ! ---------------------------------------------------------------
                        ! Part VII: Transfer leaf flux rates (like maintenance respiration,
                        ! carbon assimilation and conductance) that are defined by the 
                        ! leaf layer (which is area independent, ie /m2) onto each cohort
                        ! (where the rates become per cohort, ie /individual). Most likely
                        ! a sum over layers.
                        ! ---------------------------------------------------------------
                        nv = currentCohort%nv
                        call ScaleLeafLayerFluxToCohort(nv,                                    & !in
                                                        currentPatch%psn_z(cl,ft,1:nv),        & !in
                                                        lmr_z(1:nv,ft,cl),                     & !in
                                                        rs_z(1:nv,ft,cl),                      & !in
                                                        currentPatch%elai_profile(cl,ft,1:nv), & !in
                                                        currentCohort%c_area,                  & !in
                                                        currentCohort%n,                       & !in
                                                        bc_in(s)%rb_pa(ifp),                   & !in
                                                        maintresp_reduction_factor,            & !in
                                                        currentCohort%g_sb_laweight,           & !out
                                                        currentCohort%gpp_tstep,               & !out
                                                        currentCohort%rdark,                   & !out
                                                        cohort_eleaf_area)                       !out
                        
                        ! Net Uptake does not need to be scaled, just transfer directly
                        currentCohort%ts_net_uptake(1:nv) = anet_av_z(1:nv,ft,cl) * umolC_to_kgC

                     else
                        
                        ! In this case, the cohort had no leaves, 
                        ! so no productivity,conductance, transpiration uptake
                        ! or dark respiration
                        cohort_eleaf_area       = 0.0_r8
                        currentCohort%gpp_tstep = 0.0_r8 
                        currentCohort%rdark = 0.0_r8 
                        currentCohort%g_sb_laweight = 0.0_r8 
                        currentCohort%ts_net_uptake(:) = 0.0_r8
                        
                     end if  ! if(currentPatch%canopy_mask(cl,ft) == 1)then
                     

                     ! ------------------------------------------------------------------
                     ! Part VIII: Calculate maintenance respiration in the sapwood and
                     ! fine root pools.
                     ! ------------------------------------------------------------------
                
                     ! Calculate the amount of nitrogen in the above and below ground 
                     ! stem and root pools, used for maint resp
                     ! We are using the fine-root C:N ratio as an approximation for
                     ! the sapwood pools.
                     ! Units are in (kgN/plant)
                     ! ------------------------------------------------------------------
                     live_stem_n = EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * currentCohort%bsw / &
                           frootcn(currentCohort%pft)
                     live_croot_n = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * currentCohort%bsw  / &
                           frootcn(currentCohort%pft)
                     froot_n       = currentCohort%br / frootcn(currentCohort%pft) 
                     
                     
                     !------------------------------------------------------------------------------
                     ! Calculate Whole Plant Respiration
                     ! (this doesn't really need to be in this iteration at all, surely?)
                     ! Response: (RGK 12-2016): I think the positioning of these calls is 
                     ! appropriate as of now.  Maintenance calculations in sapwood and roots
                     ! vary by cohort and with changing temperature at the minimum, and there are
                     ! no sub-pools chopping up those pools any finer that need to be dealt with.
                     !------------------------------------------------------------------------------
                     
                     ! Live stem MR (kgC/plant/s) (above ground sapwood)
                     ! ------------------------------------------------------------------
                     if (woody(ft) == 1) then
                        tcwood = q10**((bc_in(s)%t_veg_pa(ifp)-tfrz - 20.0_r8)/10.0_r8) 
                        ! kgC/s = kgN * kgC/kgN/s
                        currentCohort%livestem_mr  = live_stem_n * ED_val_base_mr_20 * tcwood * maintresp_reduction_factor
                     else
                        currentCohort%livestem_mr  = 0._r8
                     end if
                     
                     
                     ! Fine Root MR  (kgC/plant/s)
                     ! ------------------------------------------------------------------
                     currentCohort%froot_mr = 0._r8
                     do j = 1,hlm_numlevsoil
                        tcsoi  = q10**((bc_in(s)%t_soisno_gl(j)-tfrz - 20.0_r8)/10.0_r8)
                        currentCohort%froot_mr = currentCohort%froot_mr + &
                              froot_n * ED_val_base_mr_20 * tcsoi * currentPatch%rootfr_ft(ft,j) * maintresp_reduction_factor
                     enddo
                     
                     ! Coarse Root MR (kgC/plant/s) (below ground sapwood)
                     ! ------------------------------------------------------------------
                     if (woody(ft) == 1) then
                        currentCohort%livecroot_mr = 0._r8
                        do j = 1,hlm_numlevsoil
                           ! Soil temperature used to adjust base rate of MR
                           tcsoi  = q10**((bc_in(s)%t_soisno_gl(j)-tfrz - 20.0_r8)/10.0_r8)
                           currentCohort%livecroot_mr = currentCohort%livecroot_mr + &
                                 live_croot_n * ED_val_base_mr_20 * tcsoi * &
                                 currentPatch%rootfr_ft(ft,j) * maintresp_reduction_factor
                        enddo
                     else
                        currentCohort%livecroot_mr = 0._r8    
                     end if


                     ! ------------------------------------------------------------------
                     ! Part IX: Perform some unit conversions (rate to integrated) and
                     ! calcualate some fluxes that are sums and nets of the base fluxes
                     ! ------------------------------------------------------------------
                     
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 904 ', currentCohort%resp_m
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 905 ', currentCohort%rdark
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 906 ', currentCohort%livestem_mr
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 907 ', currentCohort%livecroot_mr
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 908 ', currentCohort%froot_mr
                        
                    

                     ! add on whole plant respiration values in kgC/indiv/s-1  
                     currentCohort%resp_m = currentCohort%livestem_mr + &
                                            currentCohort%livecroot_mr + &
                                            currentCohort%froot_mr

                     ! no drought response right now.. something like:
                     ! resp_m = resp_m * (1.0_r8 - currentPatch%btran_ft(currentCohort%pft) * &
                     !                    EDPftvarcon_inst%resp_drought_response(ft))   

                     currentCohort%resp_m = currentCohort%resp_m + currentCohort%rdark
                     
                     ! convert from kgC/indiv/s to kgC/indiv/timestep       
                     currentCohort%resp_m        = currentCohort%resp_m  * dtime
                     currentCohort%gpp_tstep     = currentCohort%gpp_tstep * dtime
                     currentCohort%ts_net_uptake = currentCohort%ts_net_uptake * dtime

                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 911 ', currentCohort%gpp_tstep
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 912 ', currentCohort%resp_tstep
                     if ( DEBUG ) write(fates_log(),*) 'EDPhoto 913 ', currentCohort%resp_m
                     
                     currentCohort%resp_g     = EDPftvarcon_inst%grperc(ft) * &
                                                (max(0._r8,currentCohort%gpp_tstep - &
                                                currentCohort%resp_m))
                     currentCohort%resp_tstep = currentCohort%resp_m + &
                                                currentCohort%resp_g ! kgC/indiv/ts
                     currentCohort%npp_tstep  = currentCohort%gpp_tstep - &
                                                currentCohort%resp_tstep  ! kgC/indiv/ts
                     
                     ! Accumulate the combined conductance (stomatal+leaf boundary layer)
                     ! Note that currentCohort%g_sb_laweight is weighted by the leaf area 
                     ! of each cohort and has units of [m/s] * [m2 leaf]

                     g_sb_leaves  = g_sb_leaves + currentCohort%g_sb_laweight
                     
                     ! Accumulate the total effective leaf area from all cohorts
                     ! in this patch. Normalize by canopy area outside the loop
                     check_elai   = check_elai  + cohort_eleaf_area
                     
                     currentCohort => currentCohort%shorter
                     
                  enddo  ! end cohort loop.   
               end if !count_cohorts is more than zero.
               
               check_elai = check_elai / currentPatch%total_canopy_area
               elai       = calc_areaindex(currentPatch,'elai')

               ! Normalize canopy total conductance by the effective LAI
               ! The value here was integrated over each cohort x leaf layer
               ! and was weighted by m2 of effective leaf area for each layer
               
               if(check_elai>tiny(check_elai)) then
                  
                  ! Normalize the leaf-area weighted canopy conductance
                  ! The denominator is the total effective leaf area in the canopy,
                  ! units of [m/s]*[m2] / [m2] = [m/s]
                  g_sb_leaves = g_sb_leaves / (elai*currentPatch%total_canopy_area)
                  
                  if( g_sb_leaves > (1._r8/rsmax0) ) then 
                     
                     ! Combined mean leaf resistance is the inverse of mean leaf conductance
                     r_sb_leaves  = 1.0_r8/g_sb_leaves
                     
                     if (r_sb_leaves<bc_in(s)%rb_pa(ifp)) then
                        write(fates_log(),*) 'Combined canopy resistance was somehow smaller than'
                        write(fates_log(),*) 'its boundary layer resistance component'
                        write(fates_log(),*) 'r_sb_leaves [s/m]: ',r_sb_leaves
                        write(fates_log(),*) 'bc_in(s)%rb_pa(ifp) [s/m]: ',bc_in(s)%rb_pa(ifp)
                        call endrun(msg=errMsg(sourcefile, __LINE__))
                     end if
                     
                     ! Mean leaf stomatal resistance for all patch leaves
                     r_stomata = (r_sb_leaves - bc_in(s)%rb_pa(ifp))
                                          
                  else
                     
                     ! Here we prevent super high resistances
                     ! and use a nominal value when conductance is low
                     r_stomata = rsmax0
                     
                  end if
                  
                  ! This will be multiplied by scaled by effective LAI in the host model
                  ! when it comes time to calculate a flux rate per unit ground
                  bc_out(s)%rssun_pa(ifp) = r_stomata
                  bc_out(s)%rssha_pa(ifp) = r_stomata
                  
                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_stomata  = cf / r_stomata
                  
               else
                  
                  ! But this will prevent it from using an unintialized value
                  bc_out(s)%rssun_pa(ifp) = rsmax0
                  bc_out(s)%rssha_pa(ifp) = rsmax0

                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_stomata  = cf / rsmax0
                  
               end if
               
               ! This value is used for diagnostics, the molar form of conductance
               ! is what is used in the field usually, so we track that form
               currentPatch%c_lblayer = cf / bc_in(s)%rb_pa(ifp)
               
            end if
            
            currentPatch => currentPatch%younger
            
         end do
         
      end do !site loop
      
    end associate
  end subroutine FatesPlantRespPhotosynthDrive
  
  ! =======================================================================================
  
  subroutine LeafLayerPhotosynthesis(f_sun_lsl,         &  ! in
                                     parsun_lsl,        &  ! in
                                     parsha_lsl,        &  ! in
                                     laisun_lsl,        &  ! in
                                     laisha_lsl,        &  ! in
                                     canopy_area_lsl,   &  ! in
                                     ft,                &  ! in
                                     vcmax,             &  ! in
                                     jmax,              &  ! in
                                     tpu,               &  ! in
                                     co2_rcurve_islope, &  ! in
                                     veg_tempk,         &  ! in
                                     veg_esat,          &  ! in
                                     can_press,         &  ! in
                                     can_co2_ppress,    &  ! in
                                     can_o2_ppress,     &  ! in
                                     btran,             &  ! in
                                     bbb,               &  ! in
                                     cf,                &  ! in
                                     gb_mol,            &  ! in
                                     ceair,             &  ! in
                                     mm_kco2,           &  ! in
                                     mm_ko2,            &  ! in
                                     co2_cpoint,        &  ! in
                                     lmr,               &  ! in
                                     psn_out,           &  ! out
                                     rstoma_out,        &  ! out
                                     anet_av_out)          ! out

    ! ------------------------------------------------------------------------------------
    ! This subroutine calculates photosynthesis and stomatal conductance within each leaf 
    ! sublayer.
    ! A note on naming conventions: As this subroutine is called for every
    ! leaf-sublayer, many of the arguments are specific to that "leaf sub layer"
    ! (LSL), those variables are given a dimension tag "_lsl"
    ! Other arguments or variables may be indicative of scales broader than the LSL.
    ! ------------------------------------------------------------------------------------
    
    use EDPftvarcon       , only : EDPftvarcon_inst
    
    ! Arguments
    ! ------------------------------------------------------------------------------------
    real(r8), intent(in) :: f_sun_lsl         ! 
    real(r8), intent(in) :: parsun_lsl        ! Absorbed PAR in sunlist leaves
    real(r8), intent(in) :: parsha_lsl        ! Absorved PAR in shaded leaves
    real(r8), intent(in) :: laisun_lsl        ! LAI in sunlit leaves
    real(r8), intent(in) :: laisha_lsl        ! LAI in shaded leaves
    real(r8), intent(in) :: canopy_area_lsl   ! 
    integer,  intent(in) :: ft                ! (plant) Functional Type Index
    real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: tpu               ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(r8), intent(in) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
    real(r8), intent(in) :: veg_tempk         ! vegetation temperature
    real(r8), intent(in) :: veg_esat          ! saturation vapor pressure at veg_tempk (Pa)
    
    ! Important Note on the following gas pressures.  This photosynthesis scheme will iteratively
    ! solve for the co2 partial pressure at the leaf surface (ie in the stomata). The reference
    ! point for these input values are NOT within that boundary layer that separates the stomata from
    ! the canopy air space.  The reference point for these is on the outside of that boundary
    ! layer.  This routine, which operates at the leaf scale, makes no assumptions about what the
    ! scale of the refernce is, it could be lower atmosphere, it could be within the canopy
    ! but most likely it is the closest value one can get to the edge of the leaf's boundary
    ! layer.  We use the convention "can_" because a reference point of within the canopy
    ! ia a best reasonable scenario of where we can get that information from.
 
   real(r8), intent(in) :: can_press       ! Air pressure NEAR the surface of the leaf (Pa)
   real(r8), intent(in) :: can_co2_ppress  ! Partial pressure of CO2 NEAR the leaf surface (Pa) 
   real(r8), intent(in) :: can_o2_ppress   ! Partial pressure of O2 NEAR the leaf surface (Pa) 
   real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1) 
   real(r8), intent(in) :: bbb             ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), intent(in) :: cf              ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
   real(r8), intent(in) :: gb_mol          ! leaf boundary layer conductance (umol /m**2/s)
   real(r8), intent(in) :: ceair           ! vapor pressure of air, constrained (Pa)
   real(r8), intent(in) :: mm_kco2         ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), intent(in) :: mm_ko2          ! Michaelis-Menten constant for O2 (Pa)
   real(r8), intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
   real(r8), intent(in) :: lmr             ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
   
   real(r8), intent(out) :: psn_out        ! carbon assimilated in this leaf layer umolC/m2/s
   real(r8), intent(out) :: rstoma_out     ! stomatal resistance (1/gs_lsl) (s/m)
   real(r8), intent(out) :: anet_av_out    ! net leaf photosynthesis (umol CO2/m**2/s) 
                                           ! averaged over sun and shade leaves.  

   ! Locals
   ! ------------------------------------------------------------------------
   integer :: c3c4_path_index    ! Index for which photosynthetic pathway 
                                 ! is active.  C4 = 0,  C3 = 1
   integer :: sunsha             ! Index for differentiating sun and shade
   real(r8) :: gstoma            ! Stomatal Conductance of this leaf layer (m/s)
   real(r8) :: agross            ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8) :: anet              ! net leaf photosynthesis (umol CO2/m**2/s)
   real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
   real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
   real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
   real(r8) :: r1,r2             ! roots of quadratic equation
   real(r8) :: co2_intra_c       ! intracellular leaf CO2 (Pa)
   real(r8) :: co2_intra_c_old   ! intracellular leaf CO2 (Pa) (previous iteration)
   logical  :: loop_continue     ! Loop control variable
   integer  :: niter             ! iteration loop index
   real(r8) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
   real(r8) :: gs                ! leaf stomatal conductance (m/s)
   real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
   real(r8) :: gs_mol_err        ! gs_mol for error check
   real(r8) :: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8) :: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8) :: ap                ! product-limited (C3) or CO2-limited  
                                 ! (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
   real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
   real(r8) :: init_co2_intra_c  ! First guess intracellular co2 specific to C path
   ! Parameters
   ! ------------------------------------------------------------------------
   ! Fraction of light absorbed by non-photosynthetic pigments
   real(r8),parameter :: fnps = 0.15_r8       

   ! empirical curvature parameter for electron transport rate
   real(r8),parameter :: theta_psii = 0.7_r8   
   
   ! First guess on ratio between intracellular co2 and the atmosphere
   ! an iterator converges on actual
   real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
   real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

   ! quantum efficiency, used only for C4 (mol CO2 / mol photons) (index 0)
   real(r8),parameter,dimension(0:1) :: quant_eff = [0.05_r8,0.0_r8]

   ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8),parameter,dimension(0:1) :: theta_cj  = [0.80_r8,0.98_r8]

   ! empirical curvature parameter for ap photosynthesis co-limitation
   real(r8),parameter :: theta_ip = 0.95_r8

   associate( bb_slope  => EDPftvarcon_inst%BB_slope)    ! slope of BB relationship

     ! photosynthetic pathway: 0. = c4, 1. = c3
     c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ft))
     
     if (c3c4_path_index == 1) then
        init_co2_intra_c = init_a2l_co2_c3 * can_co2_ppress
     else
        init_co2_intra_c = init_a2l_co2_c4 * can_co2_ppress
     end if

     ! Part III: Photosynthesis and Conductance
     ! ----------------------------------------------------------------------------------
     
     if ( parsun_lsl <= 0._r8 ) then  ! night time

        anet_av_out = 0._r8
        psn_out     = 0._r8
        rstoma_out  = min(rsmax0, 1._r8/bbb * cf)
        
     else ! day time (a little bit more complicated ...)
        
!        if ( DEBUG ) write(fates_log(),*) 'EDphot 594 ',laisun_lsl
!        if ( DEBUG ) write(fates_log(),*) 'EDphot 595 ',laisha_lsl

        !is there leaf area? - (NV can be larger than 0 with only stem area if deciduous)
        if ( laisun_lsl + laisha_lsl > 0._r8 ) then 

!           if ( DEBUG ) write(fates_log(),*) '600 in laisun, laisha loop '
           
           !Loop aroun shaded and unshaded leaves          
           psn_out     = 0._r8    ! psn is accumulated across sun and shaded leaves. 
           rstoma_out  = 0._r8    ! 1/rs is accumulated across sun and shaded leaves. 
           anet_av_out = 0._r8
           gstoma  = 0._r8
           
           do  sunsha = 1,2      
              ! Electron transport rate for C3 plants. 
              ! Convert par from W/m2 to umol photons/m**2/s using the factor 4.6
              ! Convert from units of par absorbed per unit ground area to par 
              ! absorbed per unit leaf area. 
              
              if(sunsha == 1)then !sunlit
                 if(( laisun_lsl * canopy_area_lsl) > 0.0000000001_r8)then

                    qabs = parsun_lsl / (laisun_lsl * canopy_area_lsl )   
                    qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 
                    
                 else
                    qabs = 0.0_r8
                 end if
              else

                 qabs = parsha_lsl / (laisha_lsl * canopy_area_lsl)  
                 qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 

              end if

              !convert the absorbed par into absorbed par per m2 of leaf, 
              ! so it is consistant with the vcmax and lmr numbers. 
              aquad = theta_psii
              bquad = -(qabs + jmax)
              cquad = qabs * jmax
              call quadratic_f (aquad, bquad, cquad, r1, r2)
              je = min(r1,r2)

              ! Initialize intracellular co2
              co2_intra_c = init_co2_intra_c

              niter = 0
              loop_continue = .true.
              do while(loop_continue)                 
                 ! Increment iteration counter. Stop if too many iterations
                 niter = niter + 1
                 
                 ! Save old co2_intra_c
                 co2_intra_c_old = co2_intra_c
                 
                 ! Photosynthesis limitation rate calculations 
                 if (c3c4_path_index == 1)then    

                    ! C3: Rubisco-limited photosynthesis
                    ac = vcmax * max(co2_intra_c-co2_cpoint, 0._r8) / &
                          (co2_intra_c+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))

                    ! C3: RuBP-limited photosynthesis
                    aj = je * max(co2_intra_c-co2_cpoint, 0._r8) / &
                          (4._r8*co2_intra_c+8._r8*co2_cpoint)
                 
                    ! C3: Product-limited photosynthesis 
                    ap = 3._r8 * tpu
                 
                 else
                    
                    ! C4: Rubisco-limited photosynthesis
                    ac = vcmax

                    ! C4: RuBP-limited photosynthesis
                    if(sunsha == 1)then !sunlit
                       !guard against /0's in the night.
                       if((laisun_lsl * canopy_area_lsl) > 0.0000000001_r8) then   
                          aj = quant_eff(c3c4_path_index) * parsun_lsl * 4.6_r8
                          !convert from per cohort to per m2 of leaf)
                          aj = aj / (laisun_lsl * canopy_area_lsl)
                       else
                          aj = 0._r8
                       end if
                    else
                       aj = quant_eff(c3c4_path_index) * parsha_lsl * 4.6_r8
                       aj = aj / (laisha_lsl * canopy_area_lsl)
                    end if

                    ! C4: PEP carboxylase-limited (CO2-limited)
                    ap = co2_rcurve_islope * max(co2_intra_c, 0._r8) / can_press  
                    
                 end if

                 ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap
                 aquad = theta_cj(c3c4_path_index)
                 bquad = -(ac + aj)
                 cquad = ac * aj
                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 ai = min(r1,r2)

                 aquad = theta_ip
                 bquad = -(ai + ap)
                 cquad = ai * ap
                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 agross = min(r1,r2)

                 ! Net carbon assimilation. Exit iteration if an < 0
                 anet = agross  - lmr
                 if (anet < 0._r8) then
                    loop_continue = .false.
                 end if

                 ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
                 ! With an <= 0, then gs_mol = bbb
                 
                 leaf_co2_ppress = can_co2_ppress- 1.4_r8/gb_mol * anet * can_press 
                 leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
                 aquad = leaf_co2_ppress
                 bquad = leaf_co2_ppress*(gb_mol - bbb) - bb_slope(ft) * anet * can_press
                 cquad = -gb_mol*(leaf_co2_ppress*bbb + &
                                  bb_slope(ft)*anet*can_press * ceair/ veg_esat )

                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 gs_mol = max(r1,r2)
                 
                 ! Derive new estimate for co2_intra_c
                 co2_intra_c = can_co2_ppress - anet * can_press * &
                       (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)

                 ! Check for co2_intra_c convergence. Delta co2_intra_c/pair = mol/mol. 
                 ! Multiply by 10**6 to convert to umol/mol (ppm). Exit iteration if 
                 ! convergence criteria of +/- 1 x 10**-6 ppm is met OR if at least ten 
                 ! iterations (niter=10) are completed
                 
                 if ((abs(co2_intra_c-co2_intra_c_old)/can_press*1.e06_r8 <=  2.e-06_r8) &
                       .or. niter == 5) then
                    loop_continue = .false.
                 end if
              end do !iteration loop
              
              ! End of co2_intra_c iteration.  Check for an < 0, in which case gs_mol = bbb
              if (anet < 0._r8) then
                 gs_mol = bbb
              end if
              
              ! Final estimates for leaf_co2_ppress and co2_intra_c 
              ! (needed for early exit of co2_intra_c iteration when an < 0)
              leaf_co2_ppress = can_co2_ppress - 1.4_r8/gb_mol * anet * can_press
              leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
              co2_intra_c = can_co2_ppress - anet * can_press * &
                            (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)
              
              ! Convert gs_mol (umol /m**2/s) to gs (m/s) and then to rs (s/m)
              gs = gs_mol / cf
              
!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 737 ', psn_out
!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 738 ', agross
!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 739 ', f_sun_lsl

              ! Accumulate total photosynthesis umol/m2 ground/s-1. 
              ! weight per unit sun and sha leaves.
              if(sunsha == 1)then !sunlit       
                 psn_out     = psn_out + agross * f_sun_lsl
                 anet_av_out = anet_av_out + anet * f_sun_lsl
                 gstoma  = gstoma + 1._r8/(min(1._r8/gs, rsmax0)) * f_sun_lsl
              else
                 psn_out = psn_out + agross * (1.0_r8-f_sun_lsl)                 
                 anet_av_out = anet_av_out + anet * (1.0_r8-f_sun_lsl) 
                 gstoma  = gstoma + &
                       1._r8/(min(1._r8/gs, rsmax0)) * (1.0_r8-f_sun_lsl) 
              end if

!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 758 ', psn_out
!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 759 ', agross
!              if ( DEBUG ) write(fates_log(),*) 'EDPhoto 760 ', f_sun_lsl
              
              ! Make sure iterative solution is correct
              if (gs_mol < 0._r8) then
                 write (fates_log(),*)'Negative stomatal conductance:'
                 write (fates_log(),*)'gs_mol= ',gs_mol
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
              
              ! Compare with Ball-Berry model: gs_mol = m * an * hs/leaf_co2_ppress p + b
              hs = (gb_mol*ceair + gs_mol* veg_esat ) / ((gb_mol+gs_mol)*veg_esat )
              gs_mol_err = bb_slope(ft)*max(anet, 0._r8)*hs/leaf_co2_ppress*can_press + bbb
              
              if (abs(gs_mol-gs_mol_err) > 1.e-01_r8) then
                 write (fates_log(),*) 'CF: Ball-Berry error check - stomatal conductance error:'
                 write (fates_log(),*) gs_mol, gs_mol_err
              end if
              
           enddo !sunsha loop

           ! This is the stomatal resistance of the leaf layer
           rstoma_out = 1._r8/gstoma
           
        else
           !No leaf area. This layer is present only because of stems. 
           ! (leaves are off, or have reduced to 0)
           psn_out = 0._r8
           rstoma_out = min(rsmax0, 1._r8/bbb * cf)
           
        end if !is there leaf area? 
        
        
     end if    ! night or day 
   end associate
   return
  end subroutine LeafLayerPhotosynthesis
 
  ! =====================================================================================

  subroutine ScaleLeafLayerFluxToCohort(nv,          & ! in   currentCohort%nv
                                        psn_llz,     & ! in   %psn_z(1:currentCohort%nv,ft,cl)
                                        lmr_llz,     & ! in   lmr_z(1:currentCohort%nv,ft,cl)
                                        rs_llz,      & ! in   rs_z(1:currentCohort%nv,ft,cl)
                                        elai_llz,    & ! in   %elai_profile(cl,ft,1:currentCohort%nv)
                                        c_area,      & ! in   currentCohort%c_area
                                        nplant,      & ! in   currentCohort%n
                                        rb,          & ! in   bc_in(s)%rb_pa(ifp)
                                        maintresp_reduction_factor, & ! in 
                                        g_sb_laweight, & ! out  currentCohort%g_sb_laweight [m/s] [m2-leaf]
                                        gpp,         &   ! out  currentCohort%gpp_tstep
                                        rdark,       &   ! out  currentCohort%rdark
                                        cohort_eleaf_area ) ! out [m2]
   
    ! ------------------------------------------------------------------------------------
    ! This subroutine effectively integrates leaf carbon fluxes over the
    ! leaf layers to give cohort totals.
    ! Some arguments have the suffix "_llz".  This indicates that the vector
    ! is stratefied in the leaf-layer (ll)  dimension, and is a portion of the calling
    ! array which has the "_z" tag, thus "llz".
    ! ------------------------------------------------------------------------------------
    
    use FatesConstantsMod, only : umolC_to_kgC
    use EDTypesMod, only : dinc_ed
    
    ! Arguments
    integer, intent(in)  :: nv               ! number of active leaf layers
    real(r8), intent(in) :: psn_llz(nv)      ! layer photosynthesis rate (GPP) [umolC/m2leaf/s]
    real(r8), intent(in) :: lmr_llz(nv)      ! layer dark respiration rate [umolC/m2leaf/s]
    real(r8), intent(in) :: rs_llz(nv)       ! leaf layer stomatal resistance [s/m]
    real(r8), intent(in) :: elai_llz(nv)     ! exposed LAI per layer [m2 leaf/ m2 pft footprint]
    real(r8), intent(in) :: c_area           ! crown area m2/m2
    real(r8), intent(in) :: nplant           ! indiv/m2
    real(r8), intent(in) :: rb               ! leaf boundary layer resistance (s/m)
    real(r8), intent(in) :: maintresp_reduction_factor  ! factor by which to reduce maintenance respiration
    real(r8), intent(out) :: g_sb_laweight      ! Combined conductance (stomatal + boundary layer) for the cohort 
                                        ! weighted by leaf area [m/s]*[m2]
    real(r8), intent(out) :: gpp        ! GPP (kgC/indiv/s)
    real(r8), intent(out) :: rdark      ! Dark Leaf Respiration (kgC/indiv/s)
    real(r8), intent(out) :: cohort_eleaf_area  ! Effective leaf area of the cohort [m2]
    
    ! GPP IN THIS SUBROUTINE IS A RATE. THE CALLING ARGUMENT IS GPP_TSTEP. AFTER THIS
    ! CALL THE RATE WILL BE MULTIPLIED BY THE INTERVAL TO GIVE THE INTEGRATED QUANT.
    
    ! Locals
    integer  :: il                       ! leaf layer index
    real(r8) :: cohort_layer_eleaf_area  ! the effective leaf area of the cohort's current layer [m2]
    
    cohort_eleaf_area = 0.0_r8
    g_sb_laweight             = 0.0_r8
    gpp               = 0.0_r8
    rdark             = 0.0_r8

    do il = 1, nv        ! Loop over the leaf layers this cohort participates in

    
       ! Cohort's total effective leaf area in this layer [m2]
       ! leaf area index of the layer [m2/m2 ground] * [m2 ground]
       ! elai_llz is the LAI for the whole PFT. Multiplying this by the ground
       ! area this cohort contributes, give the cohort's portion of the leaf
       ! area in this layer
       cohort_layer_eleaf_area = elai_llz(il) * c_area
       
       ! Increment the cohort's total effective leaf area [m2]
       cohort_eleaf_area       = cohort_eleaf_area + cohort_layer_eleaf_area
       
       ! Leaf conductance (stomatal and boundary layer)
       ! This should be the weighted average over the leaf surfaces.
       ! Since this is relevant to the stomata, its weighting should be based
       ! on total leaf area, and not really footprint area
       ! [m/s] * [m2 cohort's leaf layer]
       g_sb_laweight = g_sb_laweight + 1.0_r8/(rs_llz(il)+rb) * cohort_layer_eleaf_area
       
       ! GPP    [umolC/m2leaf/s] * [m2 leaf ] -> [umolC/s]   (This is cohort group sum)
       gpp = gpp + psn_llz(il) * cohort_layer_eleaf_area
       
       ! Dark respiration
       ! [umolC/m2leaf/s] * [m2 leaf]    (This is the cohort group sum)
       rdark = rdark + lmr_llz(il) * cohort_layer_eleaf_area
       
    end do

    ! -----------------------------------------------------------------------------------
    ! We DO NOT normalize g_sb_laweight.
    ! The units that we are passing back are [m/s] * [m2 effective leaf]
    ! We will add these up over the whole patch, and then normalized
    ! by the patch's total leaf area in the calling routine
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Convert dark respiration and GPP from [umol/s] to [kgC/plant/s]
    ! Also, apply the maintenance respiration reduction factor
    ! -----------------------------------------------------------------------------------
    
    rdark     = rdark * umolC_to_kgC * maintresp_reduction_factor / nplant
    gpp       = gpp * umolC_to_kgC / nplant
    
    if ( DEBUG ) then
       write(fates_log(),*) 'EDPhoto 816 ', gpp
       write(fates_log(),*) 'EDPhoto 817 ', psn_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 820 ', nv
       write(fates_log(),*) 'EDPhoto 821 ', elai_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 843 ', rdark
       write(fates_log(),*) 'EDPhoto 873 ', nv
       write(fates_log(),*) 'EDPhoto 874 ', cohort_eleaf_area
    endif
    
    return
  end subroutine ScaleLeafLayerFluxToCohort

  ! =====================================================================================

  function ft1_f(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft1_f

  ! =====================================================================================
  
  function fth_f(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temp function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temp function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temp inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth_f

  ! =====================================================================================

  function fth25_f(hd,se)result(ans)
    !
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES    

     use FatesConstantsMod, only : rgas => rgas_J_K_kmol
     use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temp function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25_f
  
  ! =====================================================================================
  
  subroutine quadratic_f (a, b, c, r1, r2)
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     !==============================================================================!
     ! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
     ! Computing (Cambridge University Press, Cambridge), pp. 145.
     !
     ! !REVISION HISTORY:
     ! 4/5/10: Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
     ! 7/23/16: Copied over from CLM by Ryan Knox
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
     real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
     !
     ! !LOCAL VARIABLES:
     real(r8) :: q                        ! Temporary term for quadratic solution
     !------------------------------------------------------------------------------
    
     if (a == 0._r8) then
        write (fates_log(),*) 'Quadratic solution error: a = ',a
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if
   
     if (b >= 0._r8) then
        q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
     else
        q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
     end if
   
     r1 = q / a
     if (q /= 0._r8) then
        r2 = c / q
     else
        r2 = 1.e36_r8
     end if
     
   end subroutine quadratic_f
   
   ! ====================================================================================

   subroutine quadratic_fast (a, b, c, r1, r2)
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     ! THIS METHOD SIMPLY REMOVES THE DIV0 CHECK AND ERROR REPORTING                !
     !==============================================================================!
     ! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
     ! Computing (Cambridge University Press, Cambridge), pp. 145.
     !
     ! !REVISION HISTORY:
     ! 4/5/10: Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
     ! 7/23/16: Copied over from CLM by Ryan Knox
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
     real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
     !
     ! !LOCAL VARIABLES:
     real(r8) :: q                        ! Temporary term for quadratic solution
     !------------------------------------------------------------------------------
    
   !  if (a == 0._r8) then
   !     write (fates_log(),*) 'Quadratic solution error: a = ',a
   !     call endrun(msg=errMsg(sourcefile, __LINE__))
   !  end if
   
     if (b >= 0._r8) then
        q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
     else
        q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
     end if
   
     r1 = q / a
   !  if (q /= 0._r8) then
     r2 = c / q
   !  else
   !     r2 = 1.e36_r8
   !  end if
     
   end subroutine quadratic_fast


   ! ====================================================================================

   subroutine UpdateCanopyNCanNRadPresent(currentPatch)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates two patch level quanities:
      ! currentPatch%ncan   and
      ! currentPatch%canopy_mask
      !
      ! currentPatch%ncan(:,:) is a two dimensional array that indicates
      ! the total number of leaf layers (including those that are not exposed to light)
      ! in each canopy layer and for each functional type.
      !
      ! currentPatch%nrad(:,:) is a two dimensional array that indicates
      ! the total number of EXPOSED leaf layers, but for all intents and purposes
      ! in the photosynthesis routine, this appears to be the same as %ncan...
      !
      ! currentPatch%canopy_mask(:,:) has the same dimensions, is binary, and
      ! indicates whether or not leaf layers are present (by evaluating the canopy area
      ! profile).
      ! ---------------------------------------------------------------------------------
      
      
      use EDTypesMod , only : ed_patch_type
      use EDTypesMod , only : ed_cohort_type

      ! Arguments
      type(ed_patch_type), target :: currentPatch
      type(ed_cohort_type), pointer :: currentCohort
      
      ! Locals
      integer :: cl  ! Canopy Layer Index
      integer :: ft  ! Function Type Index
      integer :: iv  ! index of the exposed leaf layer for each canopy layer and pft
      
      ! Loop through the cohorts in this patch, associate each cohort with a layer and PFT
      ! and use the cohort's memory of how many layer's it takes up to assign the maximum
      ! of the layer/pft index it is in
      ! ---------------------------------------------------------------------------------

      currentPatch%ncan(:,:) = 0
      ! redo the canopy structure algorithm to get round a 
      ! bug that is happening for site 125, FT13.
      currentCohort => currentPatch%tallest
      do while(associated(currentCohort))
         
         currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft) = &
               max(currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft), &
                   currentCohort%NV)
         
         currentCohort => currentCohort%shorter
         
      enddo !cohort   

      ! NRAD = NCAN ...
      currentPatch%nrad = currentPatch%ncan

      ! Now loop through and identify which layer and pft combo has scattering elements
      do cl = 1,nclmax
         do ft = 1,numpft
            currentPatch%canopy_mask(cl,ft) = 0
            do iv = 1, currentPatch%nrad(cl,ft);
               if(currentPatch%canopy_area_profile(cl,ft,iv) > 0._r8)then
                  currentPatch%canopy_mask(cl,ft) = 1
               end if
            end do !iv     
         enddo !ft
      enddo !cl
      
      return
   end subroutine UpdateCanopyNCanNRadPresent
   
   ! ====================================================================================
   
   subroutine GetCanopyGasParameters(can_press, &
                                     can_o2_partialpress, &
                                     veg_tempk, &
                                     air_tempk, &
                                     air_vpress, &
                                     veg_esat,   &
                                     rb,        &
                                     mm_kco2,   &
                                     mm_ko2,    &
                                     co2_cpoint, &
                                     cf,         &
                                     gb_mol, &
                                     ceair)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates the specific Michaelis Menten Parameters (pa) for CO2
      ! and O2, as well as the CO2 compentation point.
      ! ---------------------------------------------------------------------------------
      
      use FatesConstantsMod, only: umol_per_mol
      use FatesConstantsMod, only: mmol_per_mol
      use FatesConstantsMod, only: umol_per_kmol
      use FatesConstantsMod, only : rgas => rgas_J_K_kmol

      ! Arguments
      real(r8), intent(in) :: can_press           ! Air pressure within the canopy (Pa)
      real(r8), intent(in) :: can_o2_partialpress ! Partial press of o2 in the canopy (Pa)
      real(r8), intent(in) :: veg_tempk           ! The temperature of the vegetation (K)
      real(r8), intent(in) :: air_tempk           ! Temperature of canopy air (K)
      real(r8), intent(in) :: air_vpress          ! Vapor pressure of canopy air (Pa)
      real(r8), intent(in) :: veg_esat            ! Saturated vapor pressure at veg surf (Pa)
      real(r8), intent(in) :: rb                  ! Leaf Boundary layer resistance (s/m)

      real(r8), intent(out) :: mm_kco2       ! Michaelis-Menten constant for CO2 (Pa)
      real(r8), intent(out) :: mm_ko2        !  Michaelis-Menten constant for O2 (Pa)
      real(r8), intent(out) :: co2_cpoint    !  CO2 compensation point (Pa)
      real(r8), intent(out) :: cf            ! conversion factor between molar form and velocity form
                                             ! of conductance and resistance: [umol/m3] 
      real(r8), intent(out) :: gb_mol        ! leaf boundary layer conductance (umol H2O/m**2/s)
      real(r8), intent(out) :: ceair         ! vapor pressure of air, constrained (Pa)
      
      ! Locals 
      real(r8) :: kc25                ! Michaelis-Menten constant for CO2 at 25C (Pa)
      real(r8) :: ko25                ! Michaelis-Menten constant for O2 at 25C (Pa)
      real(r8) :: sco                 ! relative specificity of rubisco
      real(r8) :: cp25                ! CO2 compensation point at 25C (Pa)
      
      ! ---------------------------------------------------------------------------------
      ! Intensive values (per mol of air) 
      ! kc, ko, currentPatch, from: Bernacchi et al (2001) 
      ! Plant, Cell and Environment 24:253-259
      ! ---------------------------------------------------------------------------------
      
      real(r8), parameter :: mm_kc25_umol_per_mol       = 404.9_r8
      real(r8), parameter :: mm_ko25_mmol_per_mol       = 278.4_r8
      real(r8), parameter :: co2_cpoint_umol_per_mol    = 42.75_r8
      
      ! Activation energy, from:
      ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
      ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

      real(r8), parameter :: kcha    = 79430._r8  ! activation energy for kc (J/mol)
      real(r8), parameter :: koha    = 36380._r8  ! activation energy for ko (J/mol)
      real(r8), parameter :: cpha    = 37830._r8  ! activation energy for cp (J/mol)
      
      
      ! Derive sco from currentPatch and O2 using present-day O2 (0.209 mol/mol) and re-calculate
      ! currentPatch to account for variation in O2 using currentPatch = 0.5 O2 / sco
      
      ! FIXME (RGK 11-30-2016 THere are more constants here, but I don't have enough information
      ! about what they are or do, so I can't give them more descriptive names. Someone please
      ! fill this in when possible)
      
      kc25 = ( mm_kc25_umol_per_mol / umol_per_mol ) * can_press
      ko25 = ( mm_ko25_mmol_per_mol / mmol_per_mol ) * can_press
      sco  = 0.5_r8 * 0.209_r8 / (co2_cpoint_umol_per_mol / umol_per_mol )
      cp25 = 0.5_r8 * can_o2_partialpress / sco
      
      if( veg_tempk.gt.150_r8 .and. veg_tempk.lt.350_r8 )then
         mm_kco2       = kc25 * ft1_f(veg_tempk, kcha)
         mm_ko2         = ko25 * ft1_f(veg_tempk, koha)
         co2_cpoint     = cp25 * ft1_f(veg_tempk, cpha)
      else
         mm_kco2    = 1.0_r8
         mm_ko2     = 1.0_r8
         co2_cpoint = 1.0_r8
      end if
      
      ! ---------------------------------------------------------------------------------
      !
      ! cf is the conversion factor between molar form and velocity form
      ! of conductance and resistance: [umol/m3]
      !
      ! i.e.
      ! [m/s] * [umol/m3] -> [umol/m2/s]
      !
      ! Breakdown of the conversion factor: [ umol / m3 ]
      !
      ! Rgas [J /K /kmol]
      ! Air Potential Temperature [ K ]
      ! Canopy Pressure      [ Pa ]
      ! conversion: umol/kmol =  1e9
      !
      ! [ Pa * K * kmol umol/kmol  /  J K ] = [ Pa * umol / J ] 
      ! since: 1 Pa = 1 N / m2
      ! [ Pa * umol / J ] = [ N * umol / J m2 ]
      ! since: 1 J = 1 N * m
      ! [ N * umol / J m2 ] = [ N * umol / N m3 ]
      ! [ umol / m3 ]
      ! 
      ! --------------------------------------------------------------------------------

      cf = can_press/(rgas * air_tempk )*umol_per_kmol
      gb_mol = (1._r8/ rb) * cf           
      
      ! Constrain eair >= 0.05*esat_tv so that solution does not blow up. This ensures
      ! that hs does not go to zero. Also eair <= veg_esat so that hs <= 1
      ceair = min( max(air_vpress, 0.05_r8*veg_esat ),veg_esat )



      return
   end subroutine GetCanopyGasParameters

   ! ====================================================================================
   
   subroutine LeafLayerMaintenanceRespiration(lmr25top_ft, &
                                              nscaler,   &
                                              ft,        &
                                              veg_tempk,     &
                                              lmr)

      use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
      use EDPftvarcon         , only : EDPftvarcon_inst
 
      ! Arguments
      real(r8), intent(in)  :: lmr25top_ft  ! canopy top leaf maint resp rate at 25C 
                                            ! for this pft (umol CO2/m**2/s)
      integer,  intent(in)  :: ft           ! (plant) Functional Type Index
      real(r8), intent(in)  :: nscaler      ! Scale for leaf nitrogen profile
      real(r8), intent(in)  :: veg_tempk    ! vegetation temperature
      real(r8), intent(out) :: lmr          ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
      
      ! Locals
      real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
      
      ! Parameter
      real(r8), parameter :: lmrha = 46390._r8    ! activation energy for lmr (J/mol)
      real(r8), parameter :: lmrhd = 150650._r8   ! deactivation energy for lmr (J/mol)
      real(r8), parameter :: lmrse = 490._r8      ! entropy term for lmr (J/mol/K)
      real(r8), parameter :: lmrc = 1.15912391_r8 ! scaling factor for high 
                                                  ! temperature inhibition (25 C = 1.0)

      ! Part I: Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
      ! ----------------------------------------------------------------------------------
      lmr25 = lmr25top_ft * nscaler 
      
      if ( nint(EDpftvarcon_inst%c3psn(ft)) == 1)then
         lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
               fth_f(veg_tempk, lmrhd, lmrse, lmrc)
      else
         lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
         lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
      end if
      
      ! Any hydrodynamic limitations could go here, currently none
      ! lmr = lmr * (nothing)
      
    end subroutine LeafLayerMaintenanceRespiration
     
   ! ====================================================================================

   subroutine LeafLayerBiophysicalRates( parsun_lsl, &
                                         ft,            &
                                         vcmax25top_ft, &
                                         jmax25top_ft, &
                                         tpu25top_ft, &
                                         co2_rcurve_islope25top_ft, &
                                         nscaler,    &
                                         veg_tempk,      &
                                         btran, &
                                         vcmax, &
                                         jmax, &
                                         tpu, &
                                         co2_rcurve_islope )

      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates the localized rates of several key photosynthesis
      ! rates.  By localized, we mean specific to the plant type and leaf layer,
      ! which factors in leaf physiology, as well as environmental effects.
      ! This procedure should be called prior to iterative solvers, and should
      ! have pre-calculated the reference rates for the pfts before this.
      !
      ! The output biophysical rates are:
      ! vcmax: maximum rate of carboxilation,
      ! jmax: maximum electron transport rate,
      ! tpu: triose phosphate utilization rate and
      ! co2_rcurve_islope: initial slope of CO2 response curve (C4 plants)
      ! ---------------------------------------------------------------------------------

      use EDPftvarcon         , only : EDPftvarcon_inst
      use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm

      ! Arguments
      ! ------------------------------------------------------------------------------
      
      real(r8), intent(in) :: parsun_lsl      ! PAR absorbed in sunlit leaves for this layer
      integer,  intent(in) :: ft              ! (plant) Functional Type Index
      real(r8), intent(in) :: nscaler         ! Scale for leaf nitrogen profile
      real(r8), intent(in) :: vcmax25top_ft   ! canopy top maximum rate of carboxylation at 25C 
                                              ! for this pft (umol CO2/m**2/s)
      real(r8), intent(in) :: jmax25top_ft    ! canopy top maximum electron transport rate at 25C 
                                              ! for this pft (umol electrons/m**2/s)
      real(r8), intent(in) :: tpu25top_ft     ! canopy top triose phosphate utilization rate at 25C 
                                              ! for this pft (umol CO2/m**2/s)
      real(r8), intent(in) :: co2_rcurve_islope25top_ft      ! initial slope of CO2 response curve
                                              ! (C4 plants) at 25C, canopy top, this pft
      real(r8), intent(in) :: veg_tempk           ! vegetation temperature
      real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1) 
                                                    
      real(r8), intent(out) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
      real(r8), intent(out) :: jmax              ! maximum electron transport rate 
                                                 ! (umol electrons/m**2/s)
      real(r8), intent(out) :: tpu               ! triose phosphate utilization rate 
                                                 ! (umol CO2/m**2/s)
      real(r8), intent(out) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
      
      ! Locals
      ! -------------------------------------------------------------------------------
      real(r8) :: vcmax25             ! leaf layer: maximum rate of carboxylation at 25C
                                      ! (umol CO2/m**2/s)
      real(r8) :: jmax25              ! leaf layer: maximum electron transport rate at 25C 
                                      ! (umol electrons/m**2/s)
      real(r8) :: tpu25               ! leaf layer: triose phosphate utilization rate at 25C 
                                      ! (umol CO2/m**2/s)
      real(r8) :: co2_rcurve_islope25 ! leaf layer: Initial slope of CO2 response curve 
                                      ! (C4 plants) at 25C
      
      
      ! Parameters
      ! ---------------------------------------------------------------------------------
      real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
      real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
      real(r8) :: tpuha          ! activation energy for tpu (J/mol)
      real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
      real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
      real(r8) :: tpuhd          ! deactivation energy for tpu (J/mol)
      real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
      real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
      real(r8) :: tpuse          ! entropy term for tpu (J/mol/K)
      real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
      real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)
      real(r8) :: tpuc           ! scaling factor for high temperature inhibition (25 C = 1.0)

      vcmaxha = EDPftvarcon_inst%vcmaxha(FT)
      jmaxha  = EDPftvarcon_inst%jmaxha(FT)
      tpuha   = EDPftvarcon_inst%tpuha(FT)
      
      vcmaxhd = EDPftvarcon_inst%vcmaxhd(FT)
      jmaxhd  = EDPftvarcon_inst%jmaxhd(FT)
      tpuhd   = EDPftvarcon_inst%tpuhd(FT)
      
      vcmaxse = EDPftvarcon_inst%vcmaxse(FT)
      jmaxse  = EDPftvarcon_inst%jmaxse(FT)
      tpuse   = EDPftvarcon_inst%tpuse(FT)

      vcmaxc = fth25_f(vcmaxhd, vcmaxse)
      jmaxc  = fth25_f(jmaxhd, jmaxse)
      tpuc   = fth25_f(tpuhd, tpuse)

      if ( parsun_lsl <= 0._r8) then           ! night time
         vcmax             = 0._r8
         jmax              = 0._r8
         tpu               = 0._r8
         co2_rcurve_islope = 0._r8
      else                                     ! day time
         vcmax25 = vcmax25top_ft * nscaler
         jmax25  = jmax25top_ft * nscaler
         tpu25   = tpu25top_ft * nscaler
         co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler
         
         ! Adjust for temperature
         vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
         jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)
         tpu   = tpu25 * ft1_f(veg_tempk, tpuha) * fth_f(veg_tempk, tpuhd, tpuse, tpuc)
         
         if (nint(EDPftvarcon_inst%c3psn(ft))  /=  1) then
            vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
            vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
            vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
         end if
         !q10 response of product limited psn. 
         co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8) 
      end if
      
      ! Adjust for water limitations 
      vcmax = vcmax * btran
      
      return
    end subroutine LeafLayerBiophysicalRates

    subroutine lowstorage_maintresp_reduction(frac, pft, maintresp_reduction_factor)

      ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
      ! of this is that mortality of plants increases when storage is low because they are not able
      ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
      ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1, 
      ! and the curvature of the function is determined by a parameter.

      ! Uses
      use EDPftvarcon         , only : EDPftvarcon_inst 

      ! Arguments
      ! ------------------------------------------------------------------------------
      real(r8), intent(in) :: frac      ! ratio of storage to target leaf biomass
      integer,  intent(in) :: pft       ! what pft is this cohort?
      real(r8), intent(out) :: maintresp_reduction_factor  ! the factor by which to reduce maintenance respiration

      ! --------------------------------------------------------------------------------
      ! Parameters are at the PFT level:
      ! fates_maintresp_reduction_curvature controls the curvature of this.  
      ! If this parameter is zero, then there is no reduction until the plant dies at storage = 0.
      ! If this parameter is one, then there is a linear reduction in respiration below the storage point.
      ! Intermediate values will give some (concave-downwards) curvature.  
      !
      ! maintresp_reduction_intercept controls the maximum amount of throttling.  
      ! zero means no throttling at any point, so it turns this mechanism off completely and so 
      ! allows an entire cohort to die via negative carbon-induced termination mortality.
      ! one means complete throttling, so no maintenance respiration at all, when out of carbon.
      ! ---------------------------------------------------------------------------------
      
       if( frac .lt. 1._r8 )then
          if ( EDPftvarcon_inst%maintresp_reduction_curvature(pft) .ne. 1._r8 ) then
             maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
                  EDPftvarcon_inst%maintresp_reduction_intercept(pft) * &
                  (1._r8 - EDPftvarcon_inst%maintresp_reduction_curvature(pft)**frac) &
                  / (1._r8-EDPftvarcon_inst%maintresp_reduction_curvature(pft))
          else  ! avoid nan answer for linear case
             maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
                  EDPftvarcon_inst%maintresp_reduction_intercept(pft) * frac
          endif
             
       else
          maintresp_reduction_factor = 1._r8
       endif


    end subroutine lowstorage_maintresp_reduction

 end module FATESPlantRespPhotosynthMod
