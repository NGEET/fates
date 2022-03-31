module PRTLossFluxesMod

  use EDPftvarcon,   only : EDPftvarcon_inst
  use PRTGenericMod, only : prt_vartypes
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : fnrt_organ
  use PRTGenericMod, only : sapw_organ
  use PRTGenericMod, only : store_organ
  use PRTGenericMod, only : repro_organ
  use PRTGenericMod, only : struct_organ
  use PRTGenericMod, only : carbon_elements_list
  use PRTGenericMod, only : carbon12_element
  use PRTGenericMod, only : carbon13_element
  use PRTGenericMod, only : carbon14_element
  use PRTGenericMod, only : nitrogen_element
  use PRTGenericMod, only : phosphorus_element
  use PRTGenericMod, only : un_initialized
  use PRTGenericMod, only : check_initialized
  use PRTGenericMod, only : num_organ_types
  use PRTGenericMod, only : prt_global
  use FatesConstantsMod, only : years_per_day
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : calloc_abs_error
  use FatesConstantsMod, only : itrue
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log 
  use shr_log_mod      , only : errMsg => shr_log_errMsg

  
  implicit none
  private

  ! -------------------------------------------------------------------------------------
  ! These modules house the public functions that handle all things
  ! related to loss fluxes.  They broadly cover the two types of turnover;
  ! that which happens as events (storms, deciduous drop, herbivory
  ! fire, etc), and maintenance turnover (constant background) 
  ! of evergreens, and branchfall).
  !
  ! IMPORTANT POINTS! 
  ! Retranslocation is handled by a single
  ! flag that defines the mode for each PFT.  So there
  ! are assumptions here.  A deciduous plant does not
  ! have maintenance leaf and fine-root turnover.  An evergreen
  ! plant does not have seasonal or stress induced phenology.
  ! Therefore, the retranslocation parameter
  ! will have different meanings potentially, for each PFT. For evergreens,
  ! it will be the retranslocation during maintenance turnover. For deciduous,
  ! it is during leaf drop.
  !
  ! THIS ROUTINE ONLY DEALS WITH LOSSES OF BIOMASS FROM PLANTS THAT ARE SURVIVING
  ! AN EVENT.  IF A PLANT DIES, THEN THESE ROUTINES DO NOT HANDLE ITS FLUXES. It
  ! is however likely that an event like fire will kill a portion of a population,
  ! and damage the remaining population, these routines will assist in the latter.
  !
  ! EDPftvarcon_inst%turnover_retrans_mode
  ! -------------------------------------------------------------------------------------

  public :: PRTDeciduousTurnover
  public :: PRTMaintTurnover
  public :: PRTBurnLosses
  public :: PRTPhenologyFlush
  public :: PRTReproRelease

contains


  subroutine PRTPhenologyFlush(prt, ipft, organ_id, c_store_transfer_frac)
     
     ! ----------------------------------------------------------------------------------
     ! This subroutine is used to flush (leaves) from storage upon bud-burst.
     ! Leaves are somewhat implied here, but the function does allow for other
     ! pools (fine-roots) to be flushed from storage as well.
     ! ----------------------------------------------------------------------------------
     
     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id
     real(r8),intent(in) :: c_store_transfer_frac  ! carbon mass fraction 
                                                   ! transferred from storage
   
     integer             :: i_var                  ! variable index
     integer             :: i_var_of_organ         ! index for all variables in
                                                   ! a given organ (mostly likely
                                                   ! synonymous with diff elements)
     integer             :: i_cvar                 ! carbon variable index for leaves
                                                   ! or other potential organ of interest
     integer             :: i_pos                  ! spatial position index
     integer             :: i_store                ! storage variable index
     integer             :: i_leaf_pos             ! Flush carbon into a specific
                                                   ! leaf pool (probably 1st?)
     integer             :: i_store_pos            ! position index for net allocation
                                                   ! from retranslocatoin in/out
                                                   ! of storage
     integer             :: element_id             ! global element identifier
     real(r8)            :: mass_transfer          ! The actual mass
                                                   ! removed from storage
                                                   ! for each pool
     real(r8)            :: target_stoich          ! stoichiometry of pool of interest
     real(r8)            :: sp_target              ! target nutrient mass for element
     real(r8)            :: sp_demand              ! nutrient demand for element


     ! We currently only allow the flushing and drop of leaves.
     ! If other organs should be desired (like seasonality of fine-roots)
     ! those parameters and clauses need to be added

     !if(organ_id .ne. leaf_organ) then
     if(organ_id .ne. leaf_organ .AND. EDPftvarcon_inst%woody(ipft) == itrue) then
        write(fates_log(),*) 'Deciduous drop and re-flushing only allowed in leaves'
        write(fates_log(),*) ' leaf_organ: ',leaf_organ
        write(fates_log(),*) ' organ: ',organ_id
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

     if(prt_global%hyp_id .le. 2) then
        i_leaf_pos  = 1             ! also used for sapwood and structural for grass
        i_store_pos = 1             ! hypothesis 1/2 only have
                                    ! 1 storage pool
     else
        write(fates_log(),*) 'You picked a hypothesis that has not defined'
        write(fates_log(),*) ' how and where flushing interacts'
        write(fates_log(),*) ' with the storage pool. specifically, '
        write(fates_log(),*) ' if this hypothesis has multiple storage pools'
        write(fates_log(),*) ' to pull carbon/resources from'
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if
     

     associate(organ_map => prt_global%organ_map)

       ! Flush carbon variables first, as their transfer
       ! rates from storage is dependant on the fraction
       ! passed in by the argument.
       ! After the values are updated, we can then
       ! identify the stoichiometry targets which
       ! govern the nutrient fluxes
       
       do i_var_of_organ = 1, organ_map(organ_id)%num_vars
          
          ! The variable index
          i_var  = organ_map(organ_id)%var_id(i_var_of_organ)
          
          ! The element index of the varible of interest
          element_id = prt_global%state_descriptor(i_var)%element_id
          
          ! This will filter IN all carbon related variables
          if( any(element_id == carbon_elements_list) ) then
             
             ! No hypotheses exist for how to flush carbon isotopes
             ! yet.  Please fill this in.
             if(  (element_id == carbon13_element) .or. &
                  (element_id == carbon14_element) )then
                write(fates_log(),*) ' Phenology flushing routine does not know'
                write(fates_log(),*) ' how to handle carbon isotopes. Please'
                write(fates_log(),*) ' evaluate the code referenced in this message'
                write(fates_log(),*) ' and provide a hypothesis.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if

             ! Get the variable id of the storage pool for this element (carbon12)
             i_store = prt_global%sp_organ_map(store_organ,element_id)


             do i_pos = 1,i_leaf_pos
                
                ! Calculate the mass transferred out of storage into the pool of interest
                mass_transfer = prt%variables(i_store)%val(i_store_pos) * &
                                c_store_transfer_frac
                
                ! Increment the c pool of interest's allocation flux
                prt%variables(i_var)%net_alloc(i_pos)   = &
                     prt%variables(i_var)%net_alloc(i_pos) + mass_transfer
                
                ! Update the c pool
                prt%variables(i_var)%val(i_pos)       = &
                     prt%variables(i_var)%val(i_pos) + mass_transfer
                
                ! Increment the storage pool's allocation flux
                prt%variables(i_store)%net_alloc(i_pos) = &
                     prt%variables(i_store)%net_alloc(i_store_pos) - mass_transfer
                
                ! Update the storage c pool
                prt%variables(i_store)%val(i_pos)     = &
                     prt%variables(i_store)%val(i_store_pos) - mass_transfer
                
                
             end do
          end if
       end do


       ! This is the variable index for leaf carbon
       ! used to calculate the targets for nutrient flushing
       i_cvar = prt_global%sp_organ_map(organ_id,carbon12_element)
       if(i_cvar < 1) then
          write(fates_log(),*) 'Could not determine the carbon var id during flushing'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Transfer in other elements (nutrients)
       ! --------------------------------------------------------------------------------
       
       do i_var_of_organ = 1, organ_map(organ_id)%num_vars
          
          i_var  = organ_map(organ_id)%var_id(i_var_of_organ)
          
          ! Variable index for the element of interest
          element_id = prt_global%state_descriptor(i_var)%element_id
          
          ! This will filter OUT all carbon related elements
          if ( .not. any(element_id == carbon_elements_list)   ) then

             ! Get the variable id of the storage pool for this element
             i_store = prt_global%sp_organ_map(store_organ,element_id)
             
             ! Calculate the stoichiometry with C for this element
             
             if( element_id == nitrogen_element ) then
                target_stoich = EDPftvarcon_inst%prt_nitr_stoich_p1(ipft,organ_id)
             else if( element_id == phosphorus_element ) then
                target_stoich = EDPftvarcon_inst%prt_phos_stoich_p1(ipft,organ_id)
             else
                  write(fates_log(),*) ' Trying to calculate nutrient flushing target'
                  write(fates_log(),*) ' for element that DNE'
                  write(fates_log(),*) ' organ: ',organ_id,' element: ',element_id
                  write(fates_log(),*) 'Exiting'
                  call endrun(msg=errMsg(__FILE__, __LINE__))
             end if

             ! Loop over all of the coordinate ids
             do i_pos = 1,i_leaf_pos
                
                ! The target quanitity for this element is based on the amount
                ! of carbon
                sp_target = prt%variables(i_cvar)%val(i_pos) * target_stoich

                sp_demand = max(0.0_r8,sp_target - prt%variables(i_var)%val(i_pos))

                ! Assume that all of the storage is transferrable
                mass_transfer = min(sp_demand, prt%variables(i_store)%val(i_store_pos))

                ! Increment the pool of interest
                prt%variables(i_var)%net_alloc(i_pos) = &
                prt%variables(i_var)%net_alloc(i_pos) + mass_transfer
                
                ! Update the  pool
                prt%variables(i_var)%val(i_pos) = &
                   prt%variables(i_var)%val(i_pos) + mass_transfer

                ! Increment the store pool allocation diagnostic
                prt%variables(i_store)%net_alloc(i_store_pos) = &
                    prt%variables(i_store)%net_alloc(i_store_pos) - mass_transfer
                
                ! Update the store pool
                prt%variables(i_store)%val(i_store_pos) = &
                    prt%variables(i_store)%val(i_store_pos) - mass_transfer
             
             end do
          
           end if

       end do
       
     end associate
     return
  end subroutine PRTPhenologyFlush
  
  ! =====================================================================================

  subroutine PRTBurnLosses(prt, organ_id, mass_fraction)

    ! ----------------------------------------------------------------------------------
    ! This subroutine assumes that there is no re-translocation associated
    ! with burn. There is only one destiny for burned mass within
    ! the organ, and that is outside the plant.  
    ! It is also assumed that non PARTEH parts of the code (ie the fire-model)
    ! will decide what to do with the burned mass (i.e. sent it to the litter
    ! pool or send to atmosphere, or.. other?)
    ! ----------------------------------------------------------------------------------

    class(prt_vartypes) :: prt
    integer,intent(in)  :: organ_id
    real(r8),intent(in) :: mass_fraction

    integer             :: i_pos          ! position index
    integer             :: i_var          ! index for the variable of interest 
    integer             :: i_var_of_organ ! loop counter for all element in this organ
    integer             :: element_id     ! Element id of the turnover pool
    real(r8)            :: burned_mass    ! Burned mass of each element, in eahc
                                          ! position, in the organ of interest
     
    associate(organ_map => prt_global%organ_map)

       ! This is the total number of state variables associated
       ! with this particular organ

       do i_var_of_organ = 1, organ_map(organ_id)%num_vars
          
          i_var = organ_map(organ_id)%var_id(i_var_of_organ)
          
          element_id = prt_global%state_descriptor(i_var)%element_id
          
          ! Loop over all of the coordinate ids
          do i_pos = 1,prt_global%state_descriptor(i_var)%num_pos
             
             ! The mass that is leaving the plant
             burned_mass = mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! Track the amount of mass being burned (+ is amount lost)
             prt%variables(i_var)%burned(i_pos) = prt%variables(i_var)%burned(i_pos) &
                  + burned_mass
             
             ! Update the state of the pool to reflect the mass lost
             prt%variables(i_var)%val(i_pos)    = prt%variables(i_var)%val(i_pos) &
                  - burned_mass
             
          end do
          
       end do
       
     end associate
  end subroutine PRTBurnLosses
    

  ! =====================================================================================


  subroutine PRTReproRelease(prt, organ_id, element_id, mass_fraction, mass_out)

    ! ----------------------------------------------------------------------------------
    ! This subroutine assumes that there is no re-translocation associated
    ! with the release of reproductive tissues.
    ! We also do not have a special flux for the release of reproductive
    ! tissues.  To not confuse this with turnover, we will provide an output
    ! mass flux, and instead of tracking it, we will just set val0 to val
    ! to prevent mass imbalances.
    ! ----------------------------------------------------------------------------------

    class(prt_vartypes)  :: prt
    integer,intent(in)   :: organ_id
    integer,intent(in)   :: element_id
    real(r8),intent(in)  :: mass_fraction
    real(r8),intent(out) :: mass_out

    integer             :: i_pos        ! position index
    integer             :: i_var        ! index for the variable of interest 

     
    associate(organ_map        => prt_global%organ_map, &
              sp_organ_map     => prt_global%sp_organ_map, &
              state_descriptor => prt_global%state_descriptor)

      ! This is the total number of state variables associated
      ! with this particular organ.
      ! In the future, we may have more finely resolved reproductive
      ! tissues (ie seeds, flowers, etc). but now we just have 1.
     
      if (organ_id .ne. repro_organ) then
         write(fates_log(),*) 'Reproductive tissue releases were called'
         write(fates_log(),*) 'for a non-reproductive organ.'
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

      if (element_id .ne. carbon12_element) then
         write(fates_log(),*) 'Reproductive tissue releases were called for a element other than c12'
         write(fates_log(),*) 'Only carbon seed masses are curently handled.'
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

      ! This is the total number of state variables associated
      ! with this particular organ

      i_var = sp_organ_map(organ_id,element_id)

      ! Reproductive mass leaving the plant
      mass_out = 0.0_r8

      ! Loop over all of the coordinate ids
      do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
         
         ! The mass that is leaving the plant
         mass_out = mass_out + mass_fraction * prt%variables(i_var)%val(i_pos)
             
         ! Update the state of the pool to reflect the mass lost
         prt%variables(i_var)%val(i_pos) = prt%variables(i_var)%val(i_pos) - &
               (mass_fraction * prt%variables(i_var)%val(i_pos))
    
         ! Update the val0 (because we don't give this dedicated flux)
         ! This is somewhat of a hack
         prt%variables(i_var)%val0(i_pos) = prt%variables(i_var)%val(i_pos) - &
               prt%variables(i_var)%net_alloc(i_pos)
         
         
      end do
       
    end associate
  end subroutine PRTReproRelease

  ! ===================================================================================

  subroutine PRTDeciduousTurnover(prt,ipft,organ_id,mass_fraction)
     
     ! ---------------------------------------------------------------------------------
     ! Generic subroutine (wrapper) calling specialized routines handling
     ! the turnover of tissues in living plants (non-mortal)
     ! ---------------------------------------------------------------------------------

     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id      ! see PRTGenericMod for organ list
     real(r8),intent(in) :: mass_fraction ! The fraction of mass in this organ that should
                                          ! leave the indicated organ.
     
     ! We currently only allow the flushing and drop of leaves.
     ! If other organs should be desired (like seasonality of fine-roots)
     ! those parameters and clauses need to be added
     
     !if(organ_id .ne. leaf_organ) then
     if(organ_id .ne. leaf_organ .AND. EDPftvarcon_inst%woody(ipft) == itrue) then
        write(fates_log(),*) 'Deciduous drop and re-flushing only allowed in leaves'
        write(fates_log(),*) ' leaf_organ: ',leaf_organ
        write(fates_log(),*) ' organ: ',organ_id
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

     
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
     ! OVER-SHOT.  TO FIX THIS, EACH HYPOTHESIS NEEDS TO HAVE WRAPPER CODE
     ! TO PROVIDE A WAY TO CALCULATE MAXIMUM ALLOWABLE STORAGE.
     !
     ! ---------------------------------------------------------------------------------

     class(prt_vartypes) :: prt
     integer,intent(in)  :: ipft
     integer,intent(in)  :: organ_id            ! see PRTGenericMod for organ list
     real(r8),intent(in) :: mass_fraction       ! The fraction of mass in this organ that should
                                                ! leave the indicated organ.

     integer             :: i_var               ! index for the variable of interest 
     integer             :: i_var_of_organ      ! loop counter for all element in this organ
     integer             :: element_id          ! Element id of the turnover pool
     integer             :: store_var_id        ! Variable id of the storage pool
     integer             :: i_store_pos         ! Position index for storage
     integer             :: i_pos               ! position index (spatial)
     real(r8)            :: retrans             ! retranslocated fraction 
     real(r8)            :: turnover_mass       ! mass sent to turnover (leaves the plant)
     real(r8)            :: retranslocated_mass ! mass redistributed to storage
     

     associate(organ_map => prt_global%organ_map)

       if((organ_id == store_organ) .or. &
          (organ_id == struct_organ) .or. & 
          (organ_id == sapw_organ)) then	   
           
          if (EDPftvarcon_inst%woody(ipft) == itrue) then        
              write(fates_log(),*) 'Deciduous turnover (leaf drop, etc)'
              write(fates_log(),*) ' was specified for an unexpected organ'
              write(fates_log(),*) ' organ: ',organ_id
              write(fates_log(),*) 'Exiting'
              call endrun(msg=errMsg(__FILE__, __LINE__))        
          end if
	  
       end if

       if(prt_global%hyp_id .le. 2) then
          i_store_pos = 1             ! hypothesis 1/2 only have
                                      ! 1 storage pool
       else
          write(fates_log(),*) 'You picked a hypothesis that has not defined'
          write(fates_log(),*) ' how and where flushing interacts'
          write(fates_log(),*) ' with the storage pool. specifically, '
          write(fates_log(),*) ' if this hypothesis has multiple storage pools'
          write(fates_log(),*) ' to pull carbon/resources from'
          write(fates_log(),*) 'Exiting'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       do i_var_of_organ = 1, organ_map(organ_id)%num_vars
          
          i_var = organ_map(organ_id)%var_id(i_var_of_organ)
          
          element_id = prt_global%state_descriptor(i_var)%element_id
          
          if ( any(element_id == carbon_elements_list) ) then
             retrans = EDPftvarcon_inst%turnover_carb_retrans(ipft,organ_id)
          else if( element_id == nitrogen_element ) then
             retrans = EDPftvarcon_inst%turnover_nitr_retrans(ipft,organ_id)
          else if( element_id == phosphorus_element ) then
             retrans = EDPftvarcon_inst%turnover_phos_retrans(ipft,organ_id)
          else
             write(fates_log(),*) 'Please add a new re-translocation clause to your '
             write(fates_log(),*) ' organ x element combination'
             write(fates_log(),*) ' organ: ',leaf_organ,' element: ',element_id
             write(fates_log(),*) 'Exiting'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          
          ! Get the variable id of the storage pool for this element
          store_var_id = prt_global%sp_organ_map(store_organ,element_id)
          
          ! Loop over all of the coordinate ids
          do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos 
             
           ! The mass that is leaving the plant
             turnover_mass = (1.0_r8 - retrans) * mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! The mass that is going towards storage
             retranslocated_mass = retrans * mass_fraction * prt%variables(i_var)%val(i_pos)
             
             ! Track the amount of mass being turned over (+ is amount lost)
             prt%variables(i_var)%turnover(i_pos) = prt%variables(i_var)%turnover(i_pos) &
                  + turnover_mass
             
             ! Track the amount of mass the is being re-translocated (- is amount lost)
             prt%variables(i_var)%net_alloc(i_pos)  = prt%variables(i_var)%net_alloc(i_pos)  &
                  - retranslocated_mass
             
             ! Update the state of the pool to reflect the mass lost
             prt%variables(i_var)%val(i_pos)      = prt%variables(i_var)%val(i_pos) &
                  - (turnover_mass + retranslocated_mass) 
             
             ! Now, since re-translocation is handled by the storage pool, 
             ! we add the re-translocated mass to it
             
             prt%variables(store_var_id)%net_alloc(i_store_pos)  = &
                  prt%variables(store_var_id)%net_alloc(i_store_pos) + retranslocated_mass
             
             prt%variables(store_var_id)%val(i_store_pos)  = &
                  prt%variables(store_var_id)%val(i_store_pos) + retranslocated_mass

          end do
          
       end do
       
     end associate

     return
   end subroutine DeciduousTurnoverSimpleRetranslocation

   ! ====================================================================================
   
   subroutine PRTMaintTurnover(prt,ipft,is_drought)
      
      ! ---------------------------------------------------------------------------------
      ! Generic subroutine (wrapper) calling specialized routines handling
      ! the turnover of tissues in living plants (non-mortal)
      ! ---------------------------------------------------------------------------------
      class(prt_vartypes) :: prt
      integer,intent(in)  :: ipft
      logical,intent(in)  :: is_drought  ! Is this plant/cohort operating in a drought
                                         ! stress context?
      
      if ( int(EDPftvarcon_inst%turnover_retrans_mode(ipft)) == 1 ) then
         call MaintTurnoverSimpleRetranslocation(prt,ipft,is_drought)
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
   
   subroutine MaintTurnoverSimpleRetranslocation(prt,ipft,is_drought)

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
      
      class(prt_vartypes)  :: prt
      integer, intent(in)  :: ipft
      logical, intent(in)  :: is_drought   ! Is this plant/cohort operating in a drought
                                           ! stress context?
      
      integer  :: i_var            ! the variable index
      integer  :: element_id       ! the element associated w/ each variable
      integer  :: organ_id         ! the organ associated w/ each variable
      integer  :: i_pos            ! spatial position loop counter
      integer  :: aclass_sen_id    ! the index of the leaf age class dimension
                                   ! associated with the senescing pool
      integer  :: ipos_1           ! the first index of the "position"
                                   ! loop to cycle. For leaves, we only
                                   ! generate maintenance fluxes from the last
                                   ! senescing class; all other cases this 
                                   ! is assumed to be 1.
      
      real(r8) :: turnover         ! Actual turnover removed from each
                                   ! pool [kg]
      real(r8) :: retrans          ! A temp for the actual re-translocated mass

      ! A temp for the actual turnover removed from pool
      real(r8), dimension(num_organ_types) :: base_turnover   
      
      ! -----------------------------------------------------------------------------------
      ! Calculate the turnover rates (maybe this should be done once in the parameter
      ! check routine. Perhaps generate a rate in parameters derived?
      ! -----------------------------------------------------------------------------------

      base_turnover(:) = un_initialized

      ! All plants can have branch turnover, if branchfall is nonz-ero,
      ! which will reduce sapwood, structure and storage.
      ! -----------------------------------------------------------------------------------
      
      if ( EDPftvarcon_inst%branch_turnover(ipft) > nearzero ) then
         base_turnover(sapw_organ)   = years_per_day / EDPftvarcon_inst%branch_turnover(ipft)
         base_turnover(struct_organ) = years_per_day / EDPftvarcon_inst%branch_turnover(ipft)
         base_turnover(store_organ)  = years_per_day / EDPftvarcon_inst%branch_turnover(ipft)
      else
         base_turnover(sapw_organ)   = 0.0_r8
         base_turnover(struct_organ) = 0.0_r8
         base_turnover(store_organ)  = 0.0_r8
      end if

      ! All plants are allowed to have fine-root turnover if a non-zero
      ! life-span is selected
      ! ---------------------------------------------------------------------------------
      if ( EDPftvarcon_inst%root_long(ipft) > nearzero ) then
         base_turnover(fnrt_organ) = years_per_day / EDPftvarcon_inst%root_long(ipft)
      else
         base_turnover(fnrt_organ) = 0.0_r8
      end if


      ! The last index of the leaf longevity array contains the turnover
      ! timescale for the senescent pool.
      aclass_sen_id = size(EDPftvarcon_inst%leaf_long(ipft,:))
      
      ! Only evergreens have maintenance turnover (must also change trimming logic
      ! if we want to change this)
      ! -------------------------------------------------------------------------------------
      if ( (EDPftvarcon_inst%leaf_long(ipft,aclass_sen_id) > nearzero ) .and. &
           (EDPftvarcon_inst%evergreen(ipft) == itrue) ) then

         if(is_drought) then
            base_turnover(leaf_organ) = years_per_day / &
                  (EDPftvarcon_inst%leaf_long(ipft,aclass_sen_id) * &
                  EDPftvarcon_inst%senleaf_long_fdrought(ipft) ) 
         else
            base_turnover(leaf_organ) = years_per_day / &
                  EDPftvarcon_inst%leaf_long(ipft,aclass_sen_id)
         end if
      else
         base_turnover(leaf_organ) = 0.0_r8
      endif

      base_turnover(repro_organ)  = 0.0_r8

      do i_var = 1, prt_global%num_vars
         
         organ_id = prt_global%state_descriptor(i_var)%organ_id
         element_id = prt_global%state_descriptor(i_var)%element_id

         if ( any(element_id == carbon_elements_list) ) then
            retrans = EDPftvarcon_inst%turnover_carb_retrans(ipft,organ_id)
         else if( element_id == nitrogen_element ) then
            retrans = EDPftvarcon_inst%turnover_nitr_retrans(ipft,organ_id)
         else if( element_id == phosphorus_element ) then
            retrans = EDPftvarcon_inst%turnover_phos_retrans(ipft,organ_id)
         else
            write(fates_log(),*) 'Please add a new re-translocation clause to your '
            write(fates_log(),*) ' organ x element combination'
            write(fates_log(),*) ' organ: ',organ_id,' element: ',element_id
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         if(base_turnover(organ_id) < check_initialized) then
            write(fates_log(),*) 'A maintenance turnover rate for the organ'
            write(fates_log(),*) ' was not specified....'
            write(fates_log(),*) ' organ: ',organ_id,' element: ',element_id
            write(fates_log(),*) ' base turnover rate: ',base_turnover(organ_id)
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
         ! Loop over all of the coordinate ids

         if(retrans<0.0 .or. retrans>1.0) then
            write(fates_log(),*) 'Unacceptable retranslocation calculated'
            write(fates_log(),*) ' organ: ',organ_id,' element: ',element_id
            write(fates_log(),*) ' retranslocation fraction: ',retrans
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         ! Hypotheses 1 & 2 assume that the leaf pools are statified by age
         ! We only generate turnover from the last (senescing) position
         if((organ_id .eq. leaf_organ)) then
            if (prt_global%hyp_id .le. 2) then
               ipos_1 = prt_global%state_descriptor(i_var)%num_pos 
            else
               write(fates_log(),*) 'Unhandled Leaf maintenance turnover condition'
               write(fates_log(),*) 'for PARTEH hypothesis id: ',prt_global%hyp_id
               call endrun(msg=errMsg(__FILE__, __LINE__))
            end if
         else
            ipos_1 = 1
         end if

         do i_pos = ipos_1, prt_global%state_descriptor(i_var)%num_pos 
            
            turnover = (1.0_r8 - retrans) * base_turnover(organ_id) * prt%variables(i_var)%val(i_pos)
      
            prt%variables(i_var)%turnover(i_pos) = prt%variables(i_var)%turnover(i_pos) + turnover
            
            prt%variables(i_var)%val(i_pos) = prt%variables(i_var)%val(i_pos)           - turnover

         end do

      end do
      
      return
   end subroutine MaintTurnoverSimpleRetranslocation





end module PRTLossFluxesMod
