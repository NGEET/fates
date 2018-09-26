module PRTLossFluxesMod

  use EDPftvarcon, only : EDPftvarcon_inst
  use PRTGenericMod, only : prt_vartypes
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : fnrt_organ
  use PRTGenericMod, only : sapw_organ
  use PRTGenericMod, only : store_organ
  use PRTGenericMod, only : repro_organ
  use PRTGenericMod, only : struct_organ
  use PRTGenericMod, only : all_carbon_species
  use PRTGenericMod, only : carbon_species    ! This is a vector
  use PRTGenericMod, only : carbon12_species
  use PRTGenericMod, only : nitrogen_species
  use PRTGenericMod, only : phosphorous_species
  use PRTGenericMod, only : un_initialized
  use PRTGenericMod, only : check_initialized
  use PRTGenericMod, only : num_organ_types
  use FatesInterfaceMod, only : hlm_freq_day

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : calloc_abs_error
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log 
  use shr_log_mod      , only : errMsg => shr_log_errMsg

  
  implicit none
  private

  ! -------------------------------------------------------------------------------------
  ! This module hosts two public functions that handle all things
  ! related to loss fluxes.  They broadly cover the two types of turnover;
  ! that which happens as events (storms, deciduous drop, herbivory
  ! fire, etc), and maintenance turnover (constant background) 
  ! of evergreens, and branchfall).
  !
  ! IMPORTANT POINTS! 
  ! Retranslocation is handled by a single
  ! flag that defines the mode for each PFT.  So there
  ! are assumptions here.  A deciduous plant does not
  ! have maintenance leaf and fine-root turnover, and vice
  ! versa.  Therefore, the retranslocation parameter
  ! will have different meanings potentially, for each PFT.
  ! 
  ! Branchfall occurs for each PFT (it may be at a reduced rate,
  ! but it will be called none-the-less).
  !
  ! THIS ROUTINE ONLY DEALS WITH LOSSES OF BIOMASS FROM PLANTS THAT ARE SURVIVING
  ! AN EVENT.  IF A PLANT DIES, THEN THIS ROUTINE DOES NOT HANDLE ITS FLUXES. It
  ! is however likely that an event like fire will kill a portion of a population,
  ! and damage the remaining population, these routines will assist in the latter.
  !
  ! EDPftvarcon_inst%turnover_retrans_mode
  ! -------------------------------------------------------------------------------------

  public :: PRTDeciduousTurnover
  public :: PRTMaintTurnover
  public :: PRTBurnLosses
  public :: PRTPhenologyFlush

contains


  subroutine PRTPhenologyFlush(prt, ipft, organ_id, c_store_transfer_frac)
     
     ! This subroutine is used to flush (leaves) from storage upon bud-burst.
     ! Leaves are somewhat implied here, but the function does allow for other
     ! pools (fine-roots) to be flushed from storage as well.
     
     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id
     real(r8),intent(in) :: c_store_transfer_frac  ! carbon mass fraction 
                                                   ! transferred from storage
   
     integer             :: i_var                  ! variable index
     integer             :: i_sp_var               ! index for all species in
                                                   ! a given organ
     integer             :: i_cvar                 ! carbon variable index
     integer             :: i_pos                  ! spatial position index
     integer             :: i_store                ! storage variable index
     integer             :: spec_id                ! global species identifier
     integer             :: num_sp_vars            ! number of species for this organ
     real(r8)            :: mass_transfer          ! The actual mass
                                                   ! removed from storage
                                                   ! for each pool
     real(r8)            :: target_stoich          ! stoichiometry of species of interest
     real(r8)            :: sp_target              ! target nutrient mass for species
     real(r8)            :: sp_demand              ! nutrient demand for species


     associate(organ_map => prt%prt_instance%organ_map)

       

       
       ! First transfer in carbon
       ! --------------------------------------------------------------------------------
       
       i_cvar = prt%prt_instance%sp_organ_map(organ_id,carbon12_species)

       ! Get the variable id of the storage pool for this species (carbon12)
       i_store = prt%prt_instance%sp_organ_map(store_organ,carbon12_species)

       ! Loop over all of the coordinate ids
       do i_pos = 1,prt%variables(i_cvar)%num_pos
          
          ! Calculate the mass transferred out of storage into the pool of interest
          mass_transfer = prt%variables(i_store)%val(i_pos) * c_store_transfer_frac
          
          ! Increment the c pool of interest
          prt%variables(i_cvar)%net_art(i_pos)   = &
                prt%variables(i_cvar)%net_art(i_pos) + mass_transfer
          
          ! Update the c pool
          prt%variables(i_cvar)%val(i_pos)       = &
                prt%variables(i_cvar)%val(i_pos) + mass_transfer
          
          ! Increment the c pool of interest
          prt%variables(i_store)%net_art(i_pos) = &
                prt%variables(i_store)%net_art(i_pos) - mass_transfer
          
          ! Update the c pool
          prt%variables(i_store)%val(i_pos)     = &
                prt%variables(i_store)%val(i_pos) - mass_transfer
          
          
       end do


       ! Transfer in other species
       ! --------------------------------------------------------------------------------

       ! This is the total number of state variables associated
       ! with this particular organ (ie carbon, nitrogen, phosphorous, ...)
       

       do i_sp_var = 1, organ_map(organ_id)%num_vars
          
          i_var  = organ_map(organ_id)%var_id(i_sp_var)
          
          ! Variable index for the species of interest
          spec_id = prt%prt_instance%state_descriptor(i_var)%spec_id
          
          if ( spec_id .ne. carbon12_species ) then

             ! Get the variable id of the storage pool for this species
             i_store = prt%prt_instance%sp_organ_map(store_organ,spec_id)
             
             ! Calculate the stoichiometry with C for this species
             
             if( spec_id == nitrogen_species ) then
                target_stoich = EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,organ_id)
             else if( spec_id == phosphorous_species ) then
                target_stoich = EDPftvarcon_inst%prt_phos_stoich_p1(ipft,organ_id)
             else
                write(fates_log(),*) ' Trying to calculate nutrient flushing target'
                write(fates_log(),*) ' for species that DNE'
                write(fates_log(),*) ' organ: ',organ_id,' species: ',spec_id
                write(fates_log(),*) 'Exiting'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if


             ! Loop over all of the coordinate ids
             do i_pos = 1,prt%variables(i_var)%num_pos
                
                ! The target quanitity for this species is based on the amount
                ! of carbon
                sp_target = prt%variables(i_cvar)%val(i_pos) * target_stoich

                sp_demand = max(0.0_r8,sp_target - prt%variables(i_var)%val(i_pos))

                ! Assume that all of the storage is transferrable
                mass_transfer = min(sp_demand, prt%variables(i_store)%val(i_pos))

                ! Increment the pool of interest
                prt%variables(i_var)%net_art(i_pos)   = &
                      prt%variables(i_var)%net_art(i_pos) + mass_transfer
                
                ! Update the c pool
                prt%variables(i_var)%val(i_pos)       = &
                      prt%variables(i_var)%val(i_pos) + mass_transfer

                ! Increment the c pool of interest
                prt%variables(i_store)%net_art(i_pos) = &
                      prt%variables(i_store)%net_art(i_pos) - mass_transfer
                
                ! Update the c pool
                prt%variables(i_store)%val(i_pos)     = &
                      prt%variables(i_store)%val(i_pos) - mass_transfer

             
             end do
          
          end if

       end do
       

       
     end associate
     return
  end subroutine PRTPhenologyFlush
  
  ! =====================================================================================

  subroutine PRTBurnLosses(prt, organ_id, mass_fraction)

    ! This subroutine assumes that there is no re-translocation associated
    ! with burn. There is only one destiny for burned mass within
    ! the organ, and that is outside the plant.  
    ! It is also assumed that non PARTEH parts of the code (ie the fire-model)
    ! will decide what to do with the burned mass (i.e. sent it to the litter
    ! pool or send to atmosphere, or.. other?)

    class(prt_vartypes) :: prt
    integer,intent(in)  :: organ_id
    real(r8),intent(in) :: mass_fraction

    integer             :: i_pos        ! position index
    integer             :: i_var        ! index for the variable of interest 
    integer             :: i_sp_var     ! loop counter for all species in this organ
    integer             :: num_sp_vars  ! Loop size for iterating over all species
    integer             :: spec_id      ! Species id of the turnover pool
    real(r8)            :: burned_mass  ! Burned mass of each species, in eahc
                                        ! position, in the organ of interest
     
    associate(organ_map => prt%prt_instance%organ_map)

       ! This is the total number of state variables associated
       ! with this particular organ

       num_sp_vars = organ_map(organ_id)%num_vars
       
       do i_sp_var = 1, num_sp_vars
          
          i_var = organ_map(organ_id)%var_id(i_sp_var)
          
          spec_id = prt%prt_instance%state_descriptor(i_var)%spec_id
          
          ! Loop over all of the coordinate ids
          do i_pos = 1,prt%variables(i_var)%num_pos
             
             ! The mass that is leaving the plant
             burned_mass = mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! Track the amount of mass being burned (+ is amount lost)
             prt%variables(i_var)%burned(i_pos) = prt%variables(i_var)%burned(i_pos) &
                  + burned_mass
             
             ! Update the state of the pool to reflect the mass lost
             prt%variables(i_var)%val(i_pos)      = prt%variables(i_var)%val(i_pos) &
                  - burned_mass
             
          end do
          
       end do
       
     end associate
    end subroutine PRTBurnLosses
    
    ! ===================================================================================


    subroutine PRTDeciduousTurnover(prt,ipft,organ_id,mass_fraction)
     
     ! ---------------------------------------------------------------------------------
     ! Generic subroutine (wrapper) calling specialized routines handling
     ! the turnover of tissues in living plants (non-mortal)
     ! ---------------------------------------------------------------------------------
     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id
     real(r8),intent(in) :: mass_fraction
     
     if ( int(EDPftvarcon_inst%turnover_retrans_mode(ipft)) == 1 ) then
        call DeciduousTurnoverSimpleRetranslocation(prt,ipft,organ_id,mass_fraction)
     else
        write(fates_log(),*) 'A retranslocation mode was specified for deciduous drop'
        write(fates_log(),*) 'that is unknown.'
        write(fates_log(),*) 'turnover_retrans_mode= ',EDPftvarcon_inst%turnover_retrans_mode(ipft)
        write(fates_log(),*) 'pft = ',ipft
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if
     
     return
   end subroutine PRTDeciduousTurnover
   

   ! ====================================================================================

   subroutine DeciduousTurnoverSimpleRetranslocation(prt,ipft,organ_id,mass_fraction)

     ! ---------------------------------------------------------------------------------
     ! Calculate losses due to deciduous turnover.
     ! the turnover of tissues in living plants (non-mortal)
     !
     ! ALERT: NO CODE IS CURRENTLY IN PLACE TO LIMIT THE AMOUNT OF CARBON OR NUTRIENT
     ! CAN BE RE-TRANSLOCATED INTO STORAGE. IT IS POSSIBLE THAT THE MAXIMUM IS BEING
     ! WAY OVER-SHOT.  TO FIX THIS, EACH HYPOTHESIS NEEDS TO HAVE WRAPPER CODE
     ! TO PROVIDE A WAY TO CALCULATE MAXIMUM ALLOWABLE STORAGE.
     !
     ! ---------------------------------------------------------------------------------

     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id
     real(r8),intent(in) :: mass_fraction

     integer             :: i_var     ! index for the variable of interest 
     integer             :: i_sp_var  ! loop counter for all species in this organ

     integer             :: num_sp_vars ! Loop size for iterating over all species
                                        ! in the organ that is turning over
     integer             :: spec_id        ! Species id of the turnover pool
     integer             :: store_var_id   ! Variable id of the storage pool
     integer             :: i_pos          ! position index (spatial)
     real(r8)            :: retrans        ! retranslocated fraction 
     real(r8)            :: turnover_mass
     real(r8)            :: retranslocated_mass
     

     associate(organ_map => prt%prt_instance%organ_map)

       if( (organ_id == store_organ) .or. &
           (organ_id == struct_organ) .or. & 
           (organ_id == sapw_organ)) then
        
          write(fates_log(),*) 'Deciduous turnover (leaf drop, etc)'
          write(fates_log(),*) ' was specified for an unexpected organ'
          write(fates_log(),*) ' organ: ',organ_id
          write(fates_log(),*) 'Exiting'
          call endrun(msg=errMsg(__FILE__, __LINE__))
          
       end if

       ! This is the total number of state variables associated
       ! with this particular organ
       num_sp_vars = organ_map(organ_id)%num_vars
       
       do i_sp_var = 1, num_sp_vars
          
          i_var = organ_map(organ_id)%var_id(i_sp_var)
          
          spec_id = prt%prt_instance%state_descriptor(i_var)%spec_id
          
          if ( any(spec_id == carbon_species) ) then
             retrans = EDPftvarcon_inst%turnover_carb_retrans_p1(ipft,organ_id)
          else if( spec_id == nitrogen_species ) then
             retrans = EDPftvarcon_inst%turnover_nitr_retrans_p1(ipft,organ_id)
          else if( spec_id == phosphorous_species ) then
             retrans = EDPftvarcon_inst%turnover_phos_retrans_p1(ipft,organ_id)
          else
             write(fates_log(),*) 'Please add a new re-translocation clause to your '
             write(fates_log(),*) ' organ x species combination'
             write(fates_log(),*) ' organ: ',leaf_organ,' species: ',spec_id
             write(fates_log(),*) 'Exiting'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          
          ! Get the variable id of the storage pool for this species
          store_var_id = prt%prt_instance%sp_organ_map(store_organ,spec_id)
          
          ! Loop over all of the coordinate ids
          do i_pos = 1,prt%variables(i_var)%num_pos
             
           ! The mass that is leaving the plant
             turnover_mass = (1.0_r8 - retrans) * mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! The mass that is going towards storage
             retranslocated_mass = retrans * mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! Track the amount of mass being turned over (+ is amount lost)
             prt%variables(i_var)%turnover(i_pos) = prt%variables(i_var)%turnover(i_pos) &
                  + turnover_mass
             
             ! Track the amount of mass the is being re-translocated (- is amount lost)
             prt%variables(i_var)%net_art(i_pos)  = prt%variables(i_var)%net_art(i_pos)  &
                  - retranslocated_mass
             
             ! Update the state of the pool to reflect the mass lost
             prt%variables(i_var)%val(i_pos)      = prt%variables(i_var)%val(i_pos) &
                  - (turnover_mass + retranslocated_mass) 
             
             ! Now, since re-translocation is handled by the storage pool, 
             ! we add the re-translocated mass to it
             
             prt%variables(store_var_id)%net_art(i_pos)  = &
                  prt%variables(store_var_id)%net_art(i_pos) + retranslocated_mass
             
             prt%variables(store_var_id)%val(i_pos)  = &
                  prt%variables(store_var_id)%val(i_pos) + retranslocated_mass
             
             
          end do
          
       end do
       
     end associate

     return
   end subroutine DeciduousTurnoverSimpleRetranslocation

   ! ====================================================================================
   
   subroutine PRTMaintTurnover(prt,ipft)
      
      ! ---------------------------------------------------------------------------------
      ! Generic subroutine (wrapper) calling specialized routines handling
      ! the turnover of tissues in living plants (non-mortal)
      ! ---------------------------------------------------------------------------------
      class(prt_vartypes) :: prt
      integer,intent(in)  :: ipft
      
      if ( int(EDPftvarcon_inst%turnover_retrans_mode(ipft)) == 1 ) then
         call MaintTurnoverSimpleRetranslocation(prt,ipft)
      else
         write(fates_log(),*) 'A maintenance/retranslocation mode was specified'
         write(fates_log(),*) 'that is unknown.'
         write(fates_log(),*) 'turnover_retrans_mode= ',EDPftvarcon_inst%turnover_retrans_mode(ipft)
         write(fates_log(),*) 'pft = ',ipft
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if
      
      return
   end subroutine PRTMaintTurnover

   ! ===================================================================================
   
   subroutine MaintTurnoverSimpleRetranslocation(prt,ipft)

      ! ---------------------------------------------------------------------------------
      ! This subroutine removes biomass from all applicable pools due to 
      ! "maintenance turnover".  Maintenance turnover, in this context
      ! is the loss of biomass on living plants, due to continuous turnover. 
      !
      ! Notes:
      ! 1) It is assumed that this is called daily.
      ! 2) This is a completely different thing compared to deciduous leaf drop,
      !    or loss of biomass from the death of the plant.
      ! 3) Since this is maintenance turnover, and not a complete drop of leaves for
      !    deciduous trees, we just re-translocate nutrients (if necessary) from the
      !    leaves and roots that leave (no pun intended), into the leaves and roots that
      !    are still rooted to the plant (pun intended). For deciduous, event-based
      !    phenology, we will re-translocate to the storage pool.
      ! 4) There are currently no reaction costs associated with re-translocation
      ! ---------------------------------------------------------------------------------
      
      class(prt_vartypes) :: prt
      integer,intent(in) :: ipft

      
      integer :: i_var
      integer :: spec_id
      integer :: organ_id
      integer :: num_sp_vars
      integer :: i_pos

      real(r8) :: turnover
      real(r8) :: leaf_turnover
      real(r8) :: fnrt_turnover
      real(r8) :: sapw_turnover
      real(r8) :: store_turnover
      real(r8) :: struct_turnover
      real(r8) :: repro_turnover
      real(r8), dimension(num_organ_types) :: base_turnover   ! A temp for the actual turnover removed from pool
      real(r8) :: retrans    ! A temp for the actual re-translocated mass

      
      num_sp_vars = size(prt%variables,1)

      ! -----------------------------------------------------------------------------------
      ! Calculate the turnover rates (maybe this should be done once in the parameter
      ! check routine. Perhaps generate a rate in parameters derived?
      ! -----------------------------------------------------------------------------------

      base_turnover(:) = un_initialized
      
      if ( EDPftvarcon_inst%branch_turnover(ipft) > nearzero ) then
         base_turnover(sapw_organ)   = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
         base_turnover(struct_organ) = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
         base_turnover(store_organ)  = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
      else
         base_turnover(sapw_organ)   = 0.0_r8
         base_turnover(struct_organ) = 0.0_r8
         base_turnover(store_organ)  = 0.0_r8
      end if

      if ( EDPftvarcon_inst%root_long(ipft) > nearzero ) then
         base_turnover(fnrt_organ) = hlm_freq_day / EDPftvarcon_inst%root_long(ipft)
      else
         base_turnover(fnrt_organ) = 0.0_r8
      end if

      if ( (EDPftvarcon_inst%leaf_long(ipft) > nearzero ) .and. &
           (EDPftvarcon_inst%evergreen(ipft) == 1) ) then
         base_turnover(leaf_organ) = hlm_freq_day / EDPftvarcon_inst%leaf_long(ipft)
      else
         base_turnover(leaf_organ) = 0.0_r8
      endif

      base_turnover(repro_organ)  = 0.0_r8

      do i_var = 1, num_sp_vars
         
         organ_id = prt%prt_instance%state_descriptor(i_var)%organ_id
         spec_id = prt%prt_instance%state_descriptor(i_var)%spec_id

         if ( any(spec_id == carbon_species) ) then
            retrans = EDPftvarcon_inst%turnover_carb_retrans_p1(ipft,organ_id)
         else if( spec_id == nitrogen_species ) then
            retrans = EDPftvarcon_inst%turnover_nitr_retrans_p1(ipft,organ_id)
         else if( spec_id == phosphorous_species ) then
            retrans = EDPftvarcon_inst%turnover_phos_retrans_p1(ipft,organ_id)
         else
            write(fates_log(),*) 'Please add a new re-translocation clause to your '
            write(fates_log(),*) ' organ x species combination'
            write(fates_log(),*) ' organ: ',organ_id,' species: ',spec_id
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         if(base_turnover(organ_id) < check_initialized) then
            write(fates_log(),*) 'A maintenance turnover rate for the organ'
            write(fates_log(),*) ' was not specified....'
            write(fates_log(),*) ' organ: ',organ_id,' species: ',spec_id
            write(fates_log(),*) ' base turnover rate: ',base_turnover(organ_id)
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
         ! Loop over all of the coordinate ids

         if(retrans<0.0 .or. retrans>1.0) then
            write(fates_log(),*) 'Unacceptable retranslocation calculated'
            write(fates_log(),*) ' organ: ',organ_id,' species: ',spec_id
            write(fates_log(),*) ' retranslocation fraction: ',retrans
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         do i_pos = 1,prt%variables(i_var)%num_pos
            
            turnover = (1.0_r8 - retrans) * base_turnover(organ_id) * prt%variables(i_var)%val(i_pos)
      
            prt%variables(i_var)%turnover(i_pos) = prt%variables(i_var)%turnover(i_pos) + turnover
            
            prt%variables(i_var)%val(i_pos) = prt%variables(i_var)%val(i_pos)           - turnover

         end do

      end do
      
      return
   end subroutine MaintTurnoverSimpleRetranslocation





end module PRTLossFluxesMod
