module ChecksBalancesMod

   use shr_kind_mod,           only : r8 => shr_kind_r8
   use shr_const_mod,          only : SHR_CONST_CDAY
   use EDtypesMod,             only : ed_site_type,ed_patch_type,ed_cohort_type
   use EDTypesMod,             only : AREA
   use EDTypesMod,             only : site_masscheck_type
   use FatesConstantsMod,      only : g_per_kg
   use PRTGenericMod,          only : all_carbon_elements
   use PRTGenericMod,          only : leaf_organ
   use PRTGenericMod,          only : fnrt_organ
   use PRTGenericMod,          only : sapw_organ
   use PRTGenericMod,          only : store_organ
   use PRTGenericMod,          only : repro_organ
   use PRTGenericMod,          only : struct_organ

   implicit none
   
   private
   public :: SummarizeNetFluxes
   public :: FATES_BGC_Carbon_Balancecheck
   public :: SiteCarbonStock
   public :: InitMassBalanceColdStart

contains
   
  subroutine InitMassBalanceColdStart(sites,bc_in)
    
    ! arguments
    type(ed_site_type), intent(inout), target :: sites(:)
    type(bc_in_type), intent(in)              :: bc_in(:)

    ! locals
    type(site_masscheck_type),pointer :: site_mass
    integer :: il
    integer :: s
    
    do il = 1,num_elements

       do s = 1,size(sites(:),dim=1)
          
          site_mass => sites(s)%mass_balance
          
          call SiteMassStock(currentSite,il,total_stock,biomass_stock,litter_stock,seed_stock)

          site_mass%stock_fates     = total_stock
          site_mass%stock_fates_old = total_stock
          
          select case(element_list(il))
          case(carbon12_type)
             
             site_mass%stock_bgc     = bc_in(s)%tot_somc +  bc_in(s)%tot_litc
             site_mass%stock_bgc_old = site_mass%stock_bgc

          case(nitrogen_type)
             write(fates_log(),*) 'Initial conditions for BGC nitrogen must be provided'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          case(phosphorus_type)
             write(fates_log(),*) 'Initial conditions for BGC nitrogen must be provided'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          case default
             write(fates_log(),*) 'An element type was specified that DNE while'
             write(fates_log(),*) 'passing in the total mass in the BGC model'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end select
          
       end do
       
    end do
    
    return
  end subroutine InitMassBalanceColdStart  



   !------------------------------------------------------------------------
   
   subroutine SummarizeNetFluxes( nsites, sites, bc_in, is_beg_day )
      
      ! Summarize the combined production and decomposition fluxes into net fluxes
      ! This is done on the fast timestep, and to be called after both daily ED calls and 
      ! fast BGC calls.  Does not include summarization of fast-timestep productivity calls 
      ! because these must be summarized prior to daily ED calls
      !
      ! Written by Charlie Koven, Feb 2016
      !
      ! !USES: 
      use FatesInterfaceMod , only : bc_in_type
     
      use EDtypesMod        , only : AREA
      !
      implicit none   
      !
      ! !ARGUMENTS    
      
      integer                                 , intent(in)    :: nsites
      type(ed_site_type)                      , intent(inout), target :: sites(nsites)
      type(bc_in_type)                        , intent(in)    :: bc_in(nsites)
      logical                                 , intent(in)    :: is_beg_day
      
      !
      ! !LOCAL VARIABLES:
      integer :: s
      type (litter_type), pointer    :: litt
      type (ed_patch_type), pointer  :: currentPatch
      type (ed_cohort_type), pointer :: currentCohort
      
      ! in FATES timesteps, because of offset between when ED and BGC reconcile the gain 
      ! and loss of litterfall carbon, (i.e. FATES reconciles it instantly, while BGC 
      ! reconciles it incrementally over the subsequent day) calculate the total 
      ! ED -> BGC flux and keep track of the last day's info for balance checking purposes

      do s = 1,nsites
         sites(s)%hr_timeintegrated  = sites(s)%hr_timeintegrated  + bc_in(s)%tot_het_resp * dtime
         sites(s)%npp_timeintegrated = sites(s)%npp_timeintegrated + sites(s)%npp * dtime
      end do
      
      if ( is_beg_day ) then

         do il = 1,num_elements
            
            do s = 1,nsites
               
               ! order of operations in the next to lines is quite important ;)
               
               site_mass%flux_fates_to_bgc_last = site_mass%flux_fates_to_bgc
               site_mass%flux_fates_to_bgc = 0._r8
               site_mass%flux_fates_to_atm = 0._r8
               site_mass%flux_fates_to_usr = 0._r8
               site_mass%flux_fates_to_hd  = 0._r8
               
               
               currentPatch => sites(s)%oldest_patch
               do while(associated(currentPatch))
                  
                  litt => currentPatch%litter(il)
                  
                  ! Fluxes from FATES to the BGC model
                  ! (litter fragmentation) 
                  ! (exudation and nurtient uptake are handled in different timesteps)
                  
                  site_mass%flux_fates_to_bgc = site_mass%flux_fates_to_bgc + &
                       (sum(litt%ag_cwd_frag(:),dim=1) + &
                       sum(sum(litt%bg_cwd_frag(:,:),dim=1),dim=1) + &
                       sum(litt%leaf_fines_frag(:),dim=1) + &
                       sum(sum(litt%root_fines_frag(:,:),dim=1),dim=1))*currentPatch%area
                  
                  ! Fluxes from fates to the atmosphere 
                  ! (burning and seed outflux)
                  site_mass%flux_fates_to_atm = site_mass%flux_fates_to_atm + &
                       site_mass%burn_flux_to_atm + &
                       site_mass%seed_outflux_atm - &
                       site_mass%seed_influx_atm

                  ! Fluxes from fates to the user source/sinks
                  site_mass%flux_fates_to_usr = site_mass%flux_fates_to_usr - &
                       sum(litt%seed_in_extern(:),dim=1)*currentPatch%area
                  
                  ! Fluxes from fates to human-dimensions (ie wood harvesting)
                  site_mass%flux_fates_to_hd = site_mass%flux_fates_to_hd + &
                       litt%harvesting*curentPatch%area
                  
                  
                  currentPatch => currentPatch%younger
               end do !currentPatch
            end do
         end do
      endif
      
      return
   end subroutine SummarizeNetFluxes
   
   ! ====================================================================================
   
   subroutine FATES_BGC_Carbon_Balancecheck(nsites, sites, bc_in, is_beg_day, dtime, nstep)
      
      ! Integrate in time the fluxes into and out of the ecosystem, and compare these 
      ! on a daily timestep to the chagne in carbon stocks of the ecosystem
      !
      ! Written by Charlie Koven, Feb 2016
      !
      ! !USES: 
      use FatesInterfaceMod , only : bc_in_type
      use EDtypesMod        , only : ed_site_type
      !
      implicit none   
      !
      ! !ARGUMENTS    
      integer                                 , intent(in)    :: nsites
      type(ed_site_type)                      , intent(inout), target :: sites(nsites)
      type(bc_in_type)                        , intent(in)    :: bc_in(nsites)
      logical                                 , intent(in)    :: is_beg_day
      real(r8)                                , intent(in)    :: dtime  ! time-step length (s)
      integer                                 , intent(in)    :: nstep  ! time-step index
      
      ! !LOCAL VARIABLES:
      type(site_site_masscheck_type), pointer :: site_mass
      real(r8) :: error_tolerance = 1.e-6_r8
      integer  :: s
      
      ! TODO: THIS INITIALIZATION SHOULD BE IN AN INITIALIZATION PART OF THE CODE
      ! COLD-START PORTION, NSTEP IS >1 FOR RESTARTS RIGHT? (RGK)

      if (nstep .le. 1) then
         ! when starting up the model, initialize the integrator variables
         
         do il = 1, num_elements
            
            do s = 1,nsites
               
               site_mass => sites(s)%mass_balance(il)
               
               
               sites(s)%totecosysc_old       = sites(s)%totecosysc
               sites(s)%totfatesc_old        = sites(s)%totfatesc
               sites(s)%totbgcc_old          = sites(s)%totbgcc
               sites(s)%hr_timeintegrated    = 0._r8
               sites(s)%npp_timeintegrated   = 0._r8
               !
               ! also initialize the ed-BGC flux variables
               sites(s)%fates_to_bgc_this_ts = 0._r8
               sites(s)%fates_to_bgc_last_ts = 0._r8
               !
               sites(s)%cbal_err_fates = 0._r8
               sites(s)%cbal_err_bgc   = 0._r8
               sites(s)%cbal_err_tot = 0._r8
            end do
         endif
         
  
      
      ! If this is on the dynamics time-step, then we calculate the balance checks
      
      if ( is_beg_day ) then

         do il=1,num_elements
         
            do s = 1,nsites
               
               site_mass => sites(s)%mass_balance(il)
               
               if(element_list(il) .eq. carbon12_type) then
                  
                  site_mass%flux_in_fates = sites(s)%npp_timeintegrated

               else
                  write(fates_log(),*) 'Please add daily integrated nutrient fluxes'
                  write(fates_log(),*) 'to the site level input mass balance'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
                  
               ! Add in litter input fluxes
               currentPatch => sites(s)%oldest_patch 
               do while(associated(currentPatch))
                  litt => currentPatch%litter(il)
                  ! kg/day
                  site_mass%flux_in_fates = site_mass%flux_in_fates + &
                        sum(litt%seed_in_extern(:),dim=1)*currentPatch%area
                  currentPatch => currentPatch%younger
               end do

               site_mass%flux_out_fates = 
               site_mass%flux_out_fates = sites(s)%fates_to_bgc_this_ts*SHR_CONST_CDAY
               
               ! "err_fates", the total mass discrepancy between the change
               !              in the total stock from one day to the next, and
               !              the integrated mass fluxes we record coming and going
               !              over the course of that day

               site_mass%err_fates = &
                    (site_mass%stock_fates - site_mass%old_stock_fates)  & ! This is the change in total stock
                    + site_mass%intgr_flux_in_fates &
                    - site_mass%intgr_flux_out_fates
               


         
            sites(s)%cbal_err_fates = sites(s)%totfatesc - & 
                  sites(s)%totfatesc_old - &
                  (sites(s)%npp_timeintegrated + &
                  sites(s)%tot_seed_rain_flux*SHR_CONST_CDAY - &
                  sites(s)%fates_to_bgc_this_ts*SHR_CONST_CDAY - &
                  sites(s)%fire_c_to_atm*SHR_CONST_CDAY)

            sites(s)%cbal_err_fates = sites(s)%cbal_err_fates / SHR_CONST_CDAY
            
            sites(s)%cbal_err_bgc = sites(s)%totbgcc - &
                  sites(s)%totbgcc_old - &
                  (sites(s)%fates_to_bgc_last_ts*SHR_CONST_CDAY - &
                  sites(s)%hr_timeintegrated)
            sites(s)%cbal_err_bgc = sites(s)%cbal_err_bgc / SHR_CONST_CDAY
            
            sites(s)%cbal_err_tot = sites(s)%totecosysc - &
                  sites(s)%totecosysc_old - &
                  (sites(s)%nbp_integrated + &
                  sites(s)%fates_to_bgc_last_ts*SHR_CONST_CDAY - &
                  sites(s)%fates_to_bgc_this_ts*SHR_CONST_CDAY)
            sites(s)%cbal_err_tot = sites(s)%cbal_err_tot / SHR_CONST_CDAY
            
            ! Send the current to the previous/last
            sites(s)%totecosysc_old = sites(s)%totecosysc
            sites(s)%totfatesc_old  = sites(s)%totfatesc
            sites(s)%totbgcc_old    = sites(s)%totbgcc
            sites(s)%npp_timeintegrated = 0._r8
            sites(s)%hr_timeintegrated = 0._r8

            site_mass%intgr_flux_in_fates = 0._r8
            site_mass%intgr_flux_out_fates = 0._r8

            
         end do
         
      endif
      
      return
   end subroutine FATES_BGC_Carbon_Balancecheck
   
   ! ==============================================================================================

   subroutine SiteMassStock(currentSite,i_element,total_stock,biomass_stock,litter_stock,seed_stock)
     
     type(ed_site_type),intent(inout),target :: currentSite
     integer,intent(in)                      :: i_element    ! This is the element index
                                                             ! in FATES (not the parteh global id)
     real(r8),intent(out)                    :: total_stock
     real(r8),intent(out)                    :: litter_stock
     real(r8),intent(out)                    :: biomass_stock
     real(r8),intent(out)                    :: seed_stock
     integer                                 :: element_id   ! parteh element id
     type(litter_type), pointer              :: litt
     type(ed_patch_type), pointer            :: currentPatch
     type(ed_cohort_type), pointer           :: currentCohort
     
     litter_stock  = 0.0_r8
     biomass_stock = 0.0_r8
     seed_stock    = 0.0_r8

     element_id = element_list(i_element)

     currentPatch => currentSite%oldest_patch 
     do while(associated(currentPatch))

        litt => currentPatch%litter(i_element)

        ! Total non-seed litter in [kg]
        litter_stock = litter_stock + currentPatch%area * &
             (sum(litt%ag_cwd)                  + &
              sum(sum(litt%bg_cwd,dim=1)) + &
              sum(litt%leaf_fines)              + &
              sum(sum(litt%root_fines,dim=1)))

        ! Total mass of viable seeds in [kg]
        seed_stock = seed_stock + currentPatch%area * sum(litt%seed)

        ! Total mass on living plants
        currentCohort => currentPatch%tallest
        do while(associated(currentCohort))
           biomass_stock =  biomass_stock + &
                 (currentCohort%prt%GetState(struct_organ,element_id) + &
                 currentCohort%prt%GetState(sapw_organ,element_id) + &
                 currentCohort%prt%GetState(leaf_organ,element_id) + &
                 currentCohort%prt%GetState(fnrt_organ,element_id) + &
                 currentCohort%prt%GetState(store_organ,element_id) + &
                 currentCohort%prt%GetState(repro_organ,element_id) ) &
                 * currentCohort%n
           currentCohort => currentCohort%shorter
        enddo !end cohort loop 
        currentPatch => currentPatch%younger
     enddo !end patch loop
     
     total_stock = biomass_stock + seed_stock + litter_stock

     return
  end subroutine SiteMassStock
  
  ! =====================================================================================
  
  subroutine CheckLitterPools(currentSite,bc_in)

    ! -----------------------------------------------------------------------------------
    !
    ! This subroutine checks that the litter pools do not have weird values.
    ! Currently, we are only checking for negatives.  
    !
    ! This is not a carbon balance check.
    ! 
    ! We will usually keep this routine turned off, unless we are actively 
    ! debugging.
    !
    ! -----------------------------------------------------------------------------------

    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(bc_in_type), intent(in)               :: bc_in
    
    ! Local variables
    type(ed_patch_type), pointer :: currentPatch
    type(litt_vartype), pointer  :: litt    ! Litter object
    integer :: il                           ! Litter element loop index

    ! We only really run this scheme if we think things are really broken.  
    ! The balance checks should be our first line of defense that are
    ! always on.

    logical, parameter :: do_litter_checks = .true.

    
    if(.not.do_litter_checks) return  

    
    do il = 1, size(currentSite%oldest_patch%litter,dim=1)
       
       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))
          
          litt => currrentPatch%litter(il)
          
          do c = 1,ncwd
             if(litt%ag_cwd(c)<0._r8) then
                write(fates_log(),*) 'In pool: ',c
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative AG CWD: ',litt%ag_cwd(c)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
                
             do ilyr = 1,numlevsoil
                if(litt%bg_cwd(c,ilyr)<0._r8) then
                   write(fates_log(),*) 'In pool: ',c
                   write(fates_log(),*) 'Soil layer: ',ilyr
                   write(fates_log(),*) 'Element id: ',element_id
                   write(fates_log(),*) 'Negative BG CWD: ',litt%bg_cwd(c,ilyr)
                   write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
                
             end do
          end do
          
          do pft = 1,numpft
             
             if(litt%seed(pft)<0._r8) then
                write(fates_log(),*) 'For PFT: ',pft
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative seed pool: ',litt%seed(pft)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             if(litt%leaf_fines(pft)<0._r8)then
                write(fates_log(),*) 'For PFT: ',pft
                write(fates_log(),*) 'Element id: ',element_id
                write(fates_log(),*) 'Negative leaf fine litter: ',litt%leaf_fines(pft)
                write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             
             do ilyr = 1,numlevsoil
                if(litt%root_fines(pft,ilyr)<0._r8)then
                   write(fates_log(),*) 'For PFT: ',pft
                   write(fates_log(),*) 'Soil layer: ',ilyr
                   write(fates_log(),*) 'Element id: ',element_id
                   write(fates_log(),*) 'Negative leaf fine litter: ',litt%root_fines(pft,ilyr)
                   write(fates_log(),*) 'lat/lon: ',currentSite%lat,currentSite%lon
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
             end do
             
          end do
          currentPatch => currentPatch%older
       end do
    end do
    
    return
  end subroutine CheckLitterPools




   
end module ChecksBalancesMod
