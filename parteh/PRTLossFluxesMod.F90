module PRTLossFluxesMod

  use EDPftvarcon, only : EDPftvarcon_inst
  use PRTGeneric, only : prt_vartypes
  use PRTGeneric, only : leaf_organ
  use PRTGeneric, only : fnrt_organ
  use PRTGeneric, only : sapw_organ
  use PRTGeneric, only : store_organ
  use PRTGeneric, only : repro_organ
  use PRTGeneric, only : truct_organ
  use PRTGeneric, only : all_carbon_species
  use PRTGeneric, only : carbon12_species
  use PRTGeneric, only : nitrogen_species
  use PRTGeneric, only : phosphorous_species
  use FatesInterfaceMod, only : hlm_freq_day

  ! These public flags specify what kind of event based
  ! turnover is happening
  
  implicit none
  private

  integer, public, parameter :: deciduous_drop_event = 1
  integer, public, parameter :: storm_event          = 2
  integer, public, parameter :: fire_event           = 3
  integer, public, parameter :: herbivory_event      = 4

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
  ! is however likely that an event like fire will kill a portion of a populatoin,
  ! and damage the remaining population, these routines will assist in the latter.
  !
  ! EDPftvarcon_inst%turnover_retrans_mode
  ! -------------------------------------------------------------------------------------

  public :: DeciduousTurnover
  public :: MaintTurnover
  public :: PlantBurnLosses

contains

  
  subroutine PlantBurnLosses(prt, ipft, organ_id, mass_fraction)

    ! This subroutine assumes that there is no re-translocation associated
    ! with burn. There is only one destiny for burned mass within
    ! the organ, and that is outside the plant.  
    ! It is also assumed that non PARTEH parts of the code (ie the fire-model)
    ! will decide what to do with the burned mass (i.e. sent it to the litter
    ! pool or send to atmosphere, or.. other?)

    class(prt_vartypes) :: prt
    integer,intent(in)  :: ipft
    integer,intent(in)  :: organ_id
    real(r8),intent(in) :: mass_fraction

    integer             :: i_pos        ! position index
    integer             :: i_var        ! index for the variable of interest 
    integer             :: i_sp_var     ! loop counter for all species in this organ
    integer             :: num_sp_var   ! Loop size for iterating over all species
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
    
    end subroutine PlantBurnLosses
    
    ! ===================================================================================


    subroutine DeciduousTurnover(prt,ipft,organ_id,mass_fraction)
     
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
   end subroutine DeciduousTurnover
   

   ! ====================================================================================

   subroutine DeciduousTurnoverSimpleRetranslocation(pft,ipft,organ_id,mass_fraction)

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

     integer             :: num_sp_var ! Loop size for iterating over all species
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
   
   subroutine MaintTurnover(prt,ipft)
      
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
   end subroutine MaintTurnover

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

      real(r8) :: base_turnover
      real(r8) :: leaf_turnover
      real(r8) :: fnrt_turnover
      real(r8) :: sapw_turnover
      real(r8) :: store_turnover
      real(r8) :: struct_turnover
      real(r8) :: repro_turnover
      real(r8) :: turnover   ! A temp for the actual turnover removed from pool
      real(r8) :: retrans    ! A temp for the actual re-translocated mass
      
      num_sp_vars = size(prt%variables,1)

      ! -----------------------------------------------------------------------------------
      ! Calculate the turnover rates
      ! -----------------------------------------------------------------------------------
      
      if ( EDPftvarcon_inst%branch_turnover(ipft) > nearzero ) then
         sapw_turnover   = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
         struct_turnover = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
         store_turnover  = hlm_freq_day / EDPftvarcon_inst%branch_turnover(ipft)
      else
         sapw_turnover   = 0.0_r8
         struct_turnover = 0.0_r8
         store_turnover  = 0.0_r8
         
      end if
      if ( EDPftvarcon_inst%root_long(ipft) > nearzero ) then
         fnrt_turnover = hlm_freq_day / EDPftvarcon_inst%root_long(ipft)
      else
         fnrt_turnover = 0.0_r8
      end if
      if ( (EDPftvarcon_inst%leaf_long(ipft) > nearzero ) .and. &
           (EDPftvarcon_inst%evergreen(ipft) == 1) ) then
         leaf_turnover = hlm_freq_day / EDPftvarcon_inst%leaf_long(ipft)
      else
         leaf_turnover = 0.0_r8
      endif

      repro_turnover  = 0.0_r8

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
            write(fates_log(),*) ' organ: ',leaf_organ,' species: ',spec_id
            write(fates_log(),*) 'Exiting'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         ! Loop over all of the coordinate ids

         do i_pos = 1,prt%variables(i_var)%num_pos
            
            turnover = (1.0_r8 - retrans) * base_turnover * prt%variables(i_var)%val(i_pos)
            
            prt%variables(i_var)%turnover(i_pos) = prt%variables(i_var)%turnover(i_pos) + turnover
            
            prt%variables(i_var)%val(i_pos) = prt%variables(i_var)%val(i_pos)           - turnover

         end do

      end do
      
      return
   end subroutine MaintTurnoverSimpleRetranslocation





end module PRTLossFluxesMod
