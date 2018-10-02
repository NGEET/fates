module PRTGenericMod

  ! ------------------------------------------------------------------------------------
  ! Plant Allocation and Reactive Transport (PART) +
  ! Extensible Hypotheses (EH) = PARTEH
  ! 
  ! Non-Specific (Generic) Classes and Functions
  ! This contains the base classes for both the variables and the "instance"
  ! This also contains science relevent procedures that are agnostic of hypothesis
  ! such as maintenance turnover and restranslocation.
  !
  ! Ryan Knox, April 2018
  !
  ! ------------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------------
  ! TO-DO: Impose a parameter check function
  !        reproduction must be priority 0 in CNP
  !
  ! TO-DO: Create a generic mapping table that will list all species
  !        of a specific organ of interest.
  ! ------------------------------------------------------------------------------------

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : calloc_abs_error
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log 
 
  use shr_log_mod      , only : errMsg => shr_log_errMsg
 
  
  implicit none

  logical, parameter :: debug = .true.

  integer, parameter :: maxlen_varname   = 128
  integer, parameter :: maxlen_varsymbol = 32
  integer, parameter :: maxlen_varunits  = 32
  integer, parameter :: len_baseunit     = 6

  ! SEND THESE TO CONSTANTS

  ! We use this parameter as the value for which we set un-initialized values
  real(r8), parameter :: un_initialized = -9.9e32_r8

  ! We use this parameter as the value for which we check un-initialized values
  real(r8), parameter :: check_initialized = -8.8e32_r8


  ! -------------------------------------------------------------------------------------
  ! IMPORTANT!
  ! All species in all organs should be expressed in terms of KILOGRAMS
  ! All rates of change are expressed in terms of kilorams / day
  ! This assumption cannot be broken!
  ! -------------------------------------------------------------------------------------
  
  character(len=len_baseunit), parameter :: mass_unit = 'kg'
  character(len=len_baseunit), parameter :: mass_rate_unit = 'kg/day'
  
  ! -------------------------------------------------------------------------------------
  ! Organ types
  ! These are public indices used to map the organs
  ! in each hypothesis to organs that acknowledged in the calling model
  ! -------------------------------------------------------------------------------------

  integer, parameter :: num_organ_types = 6
  integer, parameter :: all_organs    = 0    ! index for all organs
  integer, parameter :: leaf_organ    = 1    ! index for leaf organs
  integer, parameter :: fnrt_organ    = 2    ! index for fine-root organs
  integer, parameter :: sapw_organ    = 3    ! index for sapwood organs
  integer, parameter :: store_organ   = 4    ! index for storage organs
  integer, parameter :: repro_organ   = 5    ! index for reproductive organs
  integer, parameter :: struct_organ  = 6    ! index for structure (dead) organs

  ! -------------------------------------------------------------------------------------
  ! Species types
  ! These are public indices used to map the species in each hypothesis
  ! to the species that are acknowledged in the calling model
  ! -------------------------------------------------------------------------------------

  integer, parameter :: num_species_types     = 17    ! Total number of unique species
                                                    ! curently recognized by PARTEH
                                                    ! should be max index in list below

  ! The following list of unique public indices should be monotonic, and self-explanatory
  
  integer, parameter :: all_carbon_species  = 0
  integer, parameter :: carbon12_species    = 1
  integer, parameter :: carbon13_species    = 2
  integer, parameter :: carbon14_species    = 3
  integer, parameter :: nitrogen_species    = 4
  integer, parameter :: phosphorous_species = 5
  integer, parameter :: potassium_species   = 6
  integer, parameter :: calcium_species     = 7
  integer, parameter :: magnesium_species   = 8
  integer, parameter :: sulfur_species      = 9
  integer, parameter :: chlorine_species    = 10
  integer, parameter :: iron_species        = 11
  integer, parameter :: manganese_species   = 12
  integer, parameter :: zinc_species        = 13
  integer, parameter :: copper_species      = 14
  integer, parameter :: boron_species       = 15
  integer, parameter :: molybdenum_species  = 16
  integer, parameter :: nickel_species      = 17

  
  ! We have some lists of species or lists of organs, such as
  ! a list of all carbon species.  To keep routines simple
  ! we set a global to the maximum list size for scratch arrays.

  integer, parameter :: max_spec_per_group = 3    ! we may query these lists
                                                  ! carbon species is the biggest list
                                                  ! right now


  ! List of all carbon species, the special index "all_carbon_species"
  ! implies the following list of carbon organs
  
  integer, parameter, dimension(3) :: carbon_species   = &
       [carbon12_species, carbon13_species, carbon14_species]
  

  ! The following index specifies the maximum number of unique variables
  ! that could be described by any unique species x organ combination.  In most
  ! scenarios, this is simply 1. But for example, one may want multiple leaf
  ! layers, each representing carbon 12.   Setting this maximum high
  ! will not have a substantial impact on the memory footprint, and it will
  ! not have an effect on loop sizes because looping bounds are variables.

  integer, parameter :: max_types_per_sp_organ = 1


  ! -------------------------------------------------------------------------------------
  ! This is a generic variable type that can be used to describe all
  ! species x organ variable combinations.
  ! 
  ! It is assumed that over the control period (probably 1 day) that
  ! changes in the current state (val) relative to the value at the start of the
  ! control period (val0), are equal to the time integrated flux terms 
  ! (net_art, turnover, etc)
  !
  ! -------------------------------------------------------------------------------------

  type prt_vartype
     
     real(r8),allocatable :: val(:)       ! Instantaneous state variable           [kg]
     real(r8),allocatable :: val0(:)      ! State variable at the beginning 
                                          ! of the control period                  [kg]
     real(r8),allocatable :: net_art(:)   ! Net change due to allocation/transport [kg]
                                          ! over the control period                [kg]
     real(r8),allocatable :: turnover(:)  ! Losses rate due to turnover            [kg]
                                          ! or, any mass destined for litter
                                          ! over the control period

     real(r8),allocatable :: burned(:)    ! Losses due to burn                     [kg]

!     real(r8),allocatable :: herbiv(:)    ! Losses due to herbivory                [kg]

     ! Placeholder
     ! To save on memory, keep this commented out, or simply
     ! add this only in the extension ... ?
     ! real(r8),allocatable           :: coordinate(:,:)

     integer              :: num_pos      ! Number of pools with own position per species x organ
     
  end type prt_vartype


  ! -------------------------------------------------------------------------------------
  ! Input boundary conditions
  ! -------------------------------------------------------------------------------------

  type prt_bctype
     
     real(r8), pointer :: rval
     integer, pointer  :: ival

  end type prt_bctype


  ! -------------------------------------------------------------------------------------
  ! This generic type defines the object that attaches to the instance of any plant.
  ! ie this is the parent object.
  ! It contains the state variable object: variables
  ! as well as the boundary condition pointers bc_inout, bc_in and bc_out
  !
  ! This object also contains the bulk of the PRT routines, including 
  ! extended (hypothesis specific routines) and generic routines (eg   
  ! routines that can operate on any hypothesis) 
  !
  ! There are procedures that are specialized for each module. And then
  ! there are procedures that are supposed to be generic and should support
  ! all the different modules.
  ! -------------------------------------------------------------------------------------

  type prt_vartypes
     
     type(prt_vartype),allocatable :: variables(:)    ! The state variables
     type(prt_bctype), allocatable :: bc_inout(:)     ! These boundaries may be changed
     type(prt_bctype), allocatable :: bc_in(:)        ! These are protected
     type(prt_bctype), allocatable :: bc_out(:)       ! These are overwritten
     real(r8)                      :: ode_opt_step
     
     ! Note this is allocated only once per node/instance
     ! This really is just a pointer, not an allocatable pointer
     type(prt_instance_type), pointer :: prt_instance
     
  contains
     
     ! These are extendable procedures that have specialized
     ! content in each of the different hypotheses
     procedure :: InitAllocate        => InitAllocateBase
     procedure :: DailyPRT            => DailyPRTBase
     procedure :: FastPRT             => FastPRTBase

     ! These are generic functions that should work on all hypotheses

     procedure, non_overridable :: InitPRTVartype
     procedure, non_overridable :: FlushBCs
     procedure, non_overridable :: InitializeInitialConditions
     procedure, non_overridable :: CheckInitialConditions
     procedure, non_overridable :: RegisterBCIn
     procedure, non_overridable :: RegisterBCOut
     procedure, non_overridable :: RegisterBCInout
     procedure, non_overridable :: GetState
     procedure, non_overridable :: GetTurnover
     procedure, non_overridable :: GetBurned
     procedure, non_overridable :: GetNetART
     procedure, non_overridable :: ZeroRates
     procedure, non_overridable :: CheckMassConservation
     procedure, non_overridable :: DeallocatePRTVartypes
     procedure, non_overridable :: WeightedFusePRTVartypes
     procedure, non_overridable :: CopyPRTVartypes
  end type prt_vartypes

  ! -------------------------------------------------------------------------------------
  ! This next section contains that types that describe the whole instance. These are 
  ! things that map the variable types themselves from one model to the next, or help 
  ! decribe the arbitrary variables.  These are not instanced on every plant, they are 
  ! instanced on every model instance.
  ! -------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------
  ! This type simply packs the names and symbols associated with all
  ! the variables for any given hypothesis
  ! -------------------------------------------------------------------------------------
  
  type :: state_descriptor_type
     character(len=maxlen_varname)   :: longname
     character(len=maxlen_varsymbol) :: symbol
     integer                         :: organ_id    ! global id for organ
     integer                         :: spec_id    ! global id for species
     integer                         :: num_pos    ! number of descrete spatial positions

     ! Also, will probably need flags to define different types of groups that this variable
     ! belongs too, which will control things like fusion, normalization, when to zero, etc...

  end type state_descriptor_type
  
  
  ! This structure packs both the mapping structure and the variable descriptors
  ! --------------------------------------------------------------------------------------
  ! This array should contain the lists of indices to 
  ! the species x organ variable structure that is used to map variables to the outside
  ! world.
  !   
  !              
  !                   | carbon | nitrogen | phosphorous | .... |
  !                   ------------------------------------------
  !    leaf           |        |          |             |      |
  !    fine-root      |        |          |             |      |
  !    sapwood        |        |          |             |      |
  !    storage        |        |          |             |      |
  !    reproduction   |        |          |             |      |
  !    structure      |        |          |             |      |
  !    ....           |        |          |             |      |
  !                   ------------------------------------------
  !     
  ! -------------------------------------------------------------------------------------
  
  ! This type will help us loop through all the different variables associated
  ! with a specific organ type. Since variables are a combination of organ and
  ! species, the number of unique variables is capped at the number of species
  ! per each organ.
  
  type organ_map_type
     integer, dimension(1:num_species_types) :: var_id
     integer                                 :: num_vars
  end type organ_map_type


  type prt_instance_type
     
     ! Note that index 0 is reserved for "all" or "irrelevant"
     character(len=maxlen_varname)                             :: hyp_name

     ! This will list the specific variable ids associated with
     ! each organ
     integer, dimension(0:num_organ_types,0:num_species_types) :: sp_organ_map
     type(state_descriptor_type), allocatable                  :: state_descriptor(:)
     type(organ_map_type), dimension(1:num_organ_types)        :: organ_map

  contains
        
     procedure, non_overridable :: ZeroInstance
     procedure, non_overridable :: InitInstance

  end type prt_instance_type

  
  type(prt_instance_type),pointer :: prt_instance


contains

  ! =====================================================================================
  ! Module Functions and Subroutines
  ! =====================================================================================

   subroutine ZeroInstance(this)
      
     class(prt_instance_type)    :: this
     
     integer :: io ! Organ loop counter
     integer :: is ! Species loop counter
     
     ! First zero out the array
     do io = 1,num_organ_types
        do is = 1,num_species_types
           this%sp_organ_map(io,is)      = 0
           this%organ_map(io)%var_id(is) = 0
        end do
        this%organ_map(io)%num_vars      = 0
     end do
        
     return
  end subroutine ZeroInstance
   
  ! =====================================================================================
  
  subroutine InitInstance(this, var_id, long_name, symbol, organ_id, spec_id, num_pos)
     
     class(prt_instance_type)    :: this
     integer, intent(in)         :: var_id
     character(len=*),intent(in) :: long_name
     character(len=*),intent(in) :: symbol
     integer, intent(in)         :: organ_id
     integer, intent(in)         :: spec_id
     integer, intent(in)         :: num_pos

     ! Set the descriptions and the associated organs/species in the variable's
     ! own array

     this%state_descriptor(var_id)%longname = long_name
     this%state_descriptor(var_id)%symbol   = symbol
     this%state_descriptor(var_id)%organ_id = organ_id
     this%state_descriptor(var_id)%spec_id  = spec_id
     this%state_descriptor(var_id)%num_pos  = num_pos

     ! Set the mapping tables for the external model

     this%sp_organ_map(organ_id,spec_id) = var_id

     ! Set another map that helps to locate all the relevant pools associated
     ! with an organ

     this%organ_map(organ_id)%num_vars = this%organ_map(organ_id)%num_vars + 1
     this%organ_map(organ_id)%var_id(this%organ_map(organ_id)%num_vars) = var_id

        
     return
  end subroutine InitInstance

  ! =====================================================================================

  subroutine InitPRTVartype(this)

    class(prt_vartypes) :: this
    

    ! This subroutine should be the first call whenever a prt_vartype object is
    ! instantiated.  This routine handles the allocation (extended procedure)
    ! and then the initializing of states with bogus information, and then
    ! the flushing of all boundary conditions to null.

    call this%InitAllocate()
    call this%InitializeInitialConditions()
    call this%FlushBCs()


    return
  end subroutine InitPRTVartype

  ! =====================================================================================

  subroutine InitializeInitialConditions(this)

    class(prt_vartypes) :: this

    integer :: num_vars   ! Number of variables
    integer :: i_var      ! Variable index

    num_vars = size(this%variables,1)
    
    do i_var = 1, num_vars
       this%variables(i_var)%val(:)      = un_initialized
       this%variables(i_var)%val0(:)     = un_initialized
       this%variables(i_var)%turnover(:) = un_initialized
       this%variables(i_var)%burned(:)   = un_initialized
       this%variables(i_var)%net_art(:)  = un_initialized
    end do

    
    return
  end subroutine InitializeInitialConditions


  ! =============================================================

  subroutine CheckInitialConditions(this)
    
    ! This subroutine is called for every variable defined in each specific
    ! hypothesis.  The global index for the specific hypothesis' variable
    ! will be provided as the second argument.

    class(prt_vartypes) :: this

    integer :: n_vars     ! Number of variables
    integer :: i_var      ! index for iterating variables
    integer :: n_cor_ids  ! Number of coordinate ids
    integer :: i_cor      ! index for iterating coordinate dimension
    integer :: i_gorgan   ! The global organ id for this variable
    integer :: i_gspecies ! The global species id for this variable

    n_vars = size(this%variables,1)

    do i_var = 1, n_vars

       n_cor_ids = size(this%variables(i_var)%val,1)

       do i_cor = 1, n_cor_ids
       
          if(this%variables(i_var)%val(i_cor) < check_initialized) then

             i_gorgan   = this%prt_instance%state_descriptor(i_var)%organ_id
             i_gspecies = this%prt_instance%state_descriptor(i_var)%spec_id

             write(fates_log(),*)'Not all initial conditions for state variables'
             write(fates_log(),*)' in PRT hypothesis: ',trim(this%prt_instance%hyp_name)
             write(fates_log(),*)' were written out.'
             write(fates_log(),*)' i_var: ',i_var
             write(fates_log(),*)' i_cor: ',i_cor
             write(fates_log(),*)' organ_id:',i_gorgan
             write(fates_log(),*)' species_id',i_gspecies
             write(fates_log(),*)'Exiting'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if

       end do
    end do

    return
  end subroutine CheckInitialConditions
  
  ! =====================================================================================

  subroutine FlushBCs(this)

    ! Boundary conditions are pointers to real's and integers in the calling model.
    ! To flush these, we set all pointers to null
    
    ! Arguments
    class(prt_vartypes) :: this
    
    ! Local
    integer :: num_bc_in
    integer :: num_bc_out
    integer :: num_bc_inout
    integer :: i

    if(allocated(this%bc_in))then
       num_bc_in = size(this%bc_in,1)
       do i = 1,num_bc_in
          this%bc_in(i)%rval => null()
          this%bc_in(i)%ival => null()
       end do
    end if
    
    if(allocated(this%bc_out))then
       num_bc_out = size(this%bc_out,1)
       do i = 1,num_bc_out
          this%bc_out(i)%rval => null()
          this%bc_out(i)%ival => null()
       end do
    end if

    if(allocated(this%bc_inout))then
       num_bc_inout = size(this%bc_inout,1)
       do i = 1,num_bc_inout
          this%bc_inout(i)%rval => null()
          this%bc_inout(i)%ival => null()
       end do
    end if
       
    return
  end subroutine FlushBCs

  ! =====================================================================================

  subroutine RegisterBCIn(this,bc_id, bc_rval, bc_ival )
     
    
    ! This subroutine should be called once when PARTEH
    ! object that is bound to the plant object is first intantiated.
    ! Unless there is some reason the boundary condition pointers
    ! are changing.

    
    ! Input Arguments
    
    class(prt_vartypes)       :: this
    integer,intent(in)        :: bc_id
    real(r8),optional, intent(inout), target :: bc_rval
    integer, optional, intent(inout), target :: bc_ival
    
    if(present(bc_ival)) then
       this%bc_in(bc_id)%ival => bc_ival
    end if
       
    if(present(bc_rval)) then
       this%bc_in(bc_id)%rval => bc_rval
    end if
    

    return
  end subroutine RegisterBCIn
  

  ! =====================================================================================


  subroutine CopyPRTVartypes(this, donor_prt_obj)

    ! Arguments
    class(prt_vartypes)                     :: this
    type(prt_vartypes), intent(in), pointer :: donor_prt_obj

    ! Locals

    integer :: i_var    ! loop iterator for variable objects
    integer :: i_bc     ! loop iterator for boundary pointers

    integer :: n_vars
    integer :: num_bc_in
    integer :: num_bc_inout
    integer :: num_bc_out

    ! Here we copy over all information from a donor_prt_object into the current PRT
    ! object.   It is assumed that the current PRT object
    ! has already bee initialized ( ie. InitAllocate() )
    ! variable val0 is omitted, because it is ephemeral and used only during the
    ! allocation process

    n_vars = size(donor_prt_obj%variables,1)

    do i_var = 1, n_vars
       this%variables(i_var)%val(:)       = donor_prt_obj%variables(i_var)%val(:)
       this%variables(i_var)%val0(:)      = donor_prt_obj%variables(i_var)%val0(:)
       this%variables(i_var)%net_art(:)   = donor_prt_obj%variables(i_var)%net_art(:)
       this%variables(i_var)%turnover(:)  = donor_prt_obj%variables(i_var)%turnover(:)
       this%variables(i_var)%burned(:)    = donor_prt_obj%variables(i_var)%burned(:)
    end do

    if(allocated(this%bc_in))then
       num_bc_in = size(this%bc_in,1)
       do i_bc = 1, num_bc_in
          this%bc_in(i_bc)%ival => donor_prt_obj%bc_in(i_bc)%ival
          this%bc_in(i_bc)%rval => donor_prt_obj%bc_in(i_bc)%rval
       end do
    end if
       
    if(allocated(this%bc_out))then
       num_bc_out = size(this%bc_out,1)
       do i_bc = 1, num_bc_out
          this%bc_out(i_bc)%ival => donor_prt_obj%bc_out(i_bc)%ival
          this%bc_out(i_bc)%rval => donor_prt_obj%bc_out(i_bc)%rval
       end do
    end if

    if(allocated(this%bc_inout))then
       num_bc_inout = size(this%bc_inout,1)
       do i_bc = 1, num_bc_inout
          this%bc_inout(i_bc)%ival => donor_prt_obj%bc_inout(i_bc)%ival
          this%bc_inout(i_bc)%rval => donor_prt_obj%bc_inout(i_bc)%rval
       end do
    end if

    this%ode_opt_step = donor_prt_obj%ode_opt_step

    return
  end subroutine CopyPRTVartypes


  ! =====================================================================================

  subroutine WeightedFusePRTVartypes(this,donor_prt_obj, recipient_fuse_weight, position_id)
    
    ! This subroutine fuses two PRT objects together based on a fusion weighting
    ! assigned for the recipient (the class calling this)

    ! Arguments 
    class(prt_vartypes)                     :: this
    type(prt_vartypes), intent(in), pointer :: donor_prt_obj
    real(r8),intent(in)                     :: recipient_fuse_weight   ! This is the weighting
                                                                       ! for the recipient
    integer,intent(in),optional             :: position_id

    ! Locals
    integer :: n_vars  ! Number of variables
    integer :: i_var   ! Loop iterator over variables
    integer :: pos_id   ! coordinate id (defaults to 1)

    n_vars = size(this%variables,1)

    if(present(position_id)) then
       pos_id = position_id
    else
       pos_id = 1
    end if

    do i_var = 1, n_vars
       this%variables(i_var)%val(pos_id)  = recipient_fuse_weight * this%variables(i_var)%val(pos_id) + &
            (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%val(pos_id)
       
       this%variables(i_var)%val0(pos_id)  = recipient_fuse_weight * this%variables(i_var)%val0(pos_id) + &
            (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%val0(pos_id)

       this%variables(i_var)%net_art(pos_id)      = recipient_fuse_weight * this%variables(i_var)%net_art(pos_id) + &
            (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%net_art(pos_id)
       
       this%variables(i_var)%turnover(pos_id)    = recipient_fuse_weight * this%variables(i_var)%turnover(pos_id) + &
            (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%turnover(pos_id)

       this%variables(i_var)%burned(pos_id)    = recipient_fuse_weight * this%variables(i_var)%burned(pos_id) + &
            (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%burned(pos_id)

    end do

    this%ode_opt_step = recipient_fuse_weight * this%ode_opt_step + &
                        (1.0_r8-recipient_fuse_weight) * donor_prt_obj%ode_opt_step


    return
  end subroutine WeightedFusePRTVartypes

  ! =====================================================================================

  subroutine DeallocatePRTVartypes(this)

    class(prt_vartypes) :: this

    integer             :: n_vars
    integer             :: i_var
    
    ! Check to see if there is any value in these pools?
    ! SHould not deallocate if there is any carbon left

    n_vars = size(this%variables,1)
    
    do i_var = 1, n_vars
       deallocate(this%variables(i_var)%val)
       deallocate(this%variables(i_var)%val0)
       deallocate(this%variables(i_var)%net_art)
       deallocate(this%variables(i_var)%turnover)
       deallocate(this%variables(i_var)%burned)
    end do

    deallocate(this%variables)

    if(allocated(this%bc_in))then
       deallocate(this%bc_in)
    end if
    
    if(allocated(this%bc_out))then
       deallocate(this%bc_out)
    end if

    if(allocated(this%bc_inout))then
       deallocate(this%bc_inout)
    end if

    this%ode_opt_step = -9.0e10_r8
    
    this%prt_instance => null()

    return
  end subroutine DeallocatePRTVartypes

  ! =====================================================================================

  subroutine RegisterBCOut(this,bc_id, bc_rval, bc_ival )
    
    ! Input Arguments
    
    class(prt_vartypes)                      :: this
    integer,intent(in)                       :: bc_id
    real(r8), optional, intent(inout),target :: bc_rval
    integer, optional, intent(inout),target  :: bc_ival
    
    if(present(bc_ival)) then
       this%bc_out(bc_id)%ival => bc_ival
    end if

    if(present(bc_rval)) then
       this%bc_out(bc_id)%rval => bc_rval
    end if

    return
  end subroutine RegisterBCOut

  ! =====================================================================================

  subroutine RegisterBCInOut(this,bc_id, bc_rval, bc_ival )
    
    ! Input Arguments
    
    class(prt_vartypes)                      :: this
    integer,intent(in)                       :: bc_id
    real(r8), optional, intent(inout),target :: bc_rval
    integer, optional, intent(inout),target  :: bc_ival
    
    if(present(bc_ival)) then
       this%bc_inout(bc_id)%ival => bc_ival
    end if

    if(present(bc_rval)) then
       this%bc_inout(bc_id)%rval => bc_rval
    end if

    return
  end subroutine RegisterBCInOut

  
  ! =====================================================================================

  subroutine ZeroRates(this)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine zeros all of the rates of change for our variables.
      ! It also sets the initial value to the current state.
      ! This allows us to make mass conservation checks, where
      ! val - val0 = net_art + turnover
      ! ---------------------------------------------------------------------------------

      class(prt_vartypes) :: this

      integer :: n_vars   ! Number of variables
      integer :: i_var    ! Variable index

      n_vars = size(this%variables,1)
      do i_var = 1,n_vars
         this%variables(i_var)%val0(:)        = this%variables(i_var)%val(:)
         this%variables(i_var)%net_art(:)     = 0.0_r8
         this%variables(i_var)%turnover(:)    = 0.0_r8
         this%variables(i_var)%burned(:)      = 0.0_r8
      end do
      
    end subroutine ZeroRates

   ! ====================================================================================

   subroutine CheckMassConservation(this,ipft,position_id)

     class(prt_vartypes) :: this
     integer, intent(in) :: ipft
     integer, intent(in) :: position_id ! Helps to know where
                                        ! in the call sequence this was called

     integer :: n_vars      ! Number of variables
     integer :: i_var       ! Variable index
     integer :: i_pos       ! Position (coordinate) index

     real(r8) :: err
     real(r8) :: rel_err


     n_vars = size(this%variables,1)
     do i_var = 1,n_vars
        
        do i_pos = 1, this%variables(i_var)%num_pos
           
           err = abs((this%variables(i_var)%val(i_pos) - this%variables(i_var)%val0(i_pos)) - &
                  (this%variables(i_var)%net_art(i_pos) &
                   -this%variables(i_var)%turnover(i_pos) & 
                   -this%variables(i_var)%burned(i_pos) ))
           
           if(this%variables(i_var)%val(i_pos) > nearzero) then
              rel_err = err / this%variables(i_var)%val(i_pos)
           else
              rel_err = 0.0_r8
           end if

           if( abs(err) > calloc_abs_error ) then
              write(fates_log(),*) 'PARTEH mass conservation check failed'
              write(fates_log(),*) ' Change in mass over control period should'
              write(fates_log(),*) ' always equal the integrated fluxes.'
              write(fates_log(),*) ' pft id: ',ipft
              write(fates_log(),*) ' position id: ',position_id
              write(fates_log(),*) ' organ id: ',this%prt_instance%state_descriptor(i_var)%organ_id
              write(fates_log(),*) ' species_id: ',this%prt_instance%state_descriptor(i_var)%spec_id
              write(fates_log(),*) ' position id: ',i_pos
              write(fates_log(),*) ' symbol: ',trim(this%prt_instance%state_descriptor(i_var)%symbol)
              write(fates_log(),*) ' longname: ',trim(this%prt_instance%state_descriptor(i_var)%longname)
              write(fates_log(),*) ' err: ',err,' max error: ',calloc_abs_error
              write(fates_log(),*) ' terms: ', this%variables(i_var)%val(i_pos), &
                                               this%variables(i_var)%val0(i_pos), &
                                               this%variables(i_var)%net_art(i_pos), &
                                               this%variables(i_var)%turnover(i_pos), &
                                               this%variables(i_var)%burned(i_pos)
              write(fates_log(),*) ' Exiting.'
              call endrun(msg=errMsg(__FILE__, __LINE__))
           end if

        end do
     end do

     return
   end subroutine CheckMassConservation

   ! ====================================================================================
   
   function GetState(this, organ_id, species_id, position_id) result(sp_organ_val)
      
      ! This function returns the current amount of mass for
      ! any combination of organ and species.  If a position
      ! is provided, it will us it, but otherwise, it will sum over
      ! all dimensions.  It also can accomodate all_carbon_species, which
      ! will return the mass of all carbon isotopes combined.

      class(prt_vartypes)   :: this
      integer,intent(in)    :: organ_id
      integer,intent(in)    :: species_id
      integer,intent(in),optional :: position_id
      real(r8)              :: sp_organ_val

      integer :: i_pos
      integer :: ispec
      integer :: num_species
      integer,dimension(max_spec_per_group) :: spec_ids 
      integer :: i_var
      
      sp_organ_val = 0.0_r8
      
      if(species_id == all_carbon_species) then
         spec_ids(1:3) = carbon_species(1:3)
         num_species  = 3
      else
         num_species  = 1
         spec_ids(1) = species_id
      end if

      if(present(position_id)) then
         i_pos = position_id
      
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) sp_organ_val = sp_organ_val + this%variables(i_var)%val(i_pos)
         end do

      else
         
         do ispec = 1,num_species
            
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0)then

               do i_pos = 1, this%variables(i_var)%num_pos
                  sp_organ_val = sp_organ_val + this%variables(i_var)%val(i_pos)
               end do
            end if
               
         end do

      end if
      
      return
    end function GetState

    
   ! ====================================================================================

   
    function GetTurnover(this, organ_id, species_id, position_id) result(sp_organ_turnover)
      

      ! THis function is very similar to GetState, with the only difference that it
      ! returns the turnover mass so-far during the period of interest.


      class(prt_vartypes)   :: this
      integer,intent(in)    :: organ_id
      integer,intent(in)    :: species_id
      integer,intent(in),optional :: position_id
      real(r8)              :: sp_organ_turnover

      integer :: i_pos
      integer :: ispec
      integer :: num_species
      integer,dimension(max_spec_per_group) :: spec_ids 
      integer :: i_var
      
      sp_organ_turnover = 0.0_r8
      
      if(species_id == all_carbon_species) then
         spec_ids(1:3) = carbon_species(1:3)
         num_species  = 3
      else
         num_species  = 1
         spec_ids(1) = species_id
      end if

      if(present(position_id)) then
         i_pos = position_id
      
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) sp_organ_turnover = sp_organ_turnover + &
                 this%variables(i_var)%turnover(i_pos)
         end do

      else

         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) then
               do i_pos = 1, this%variables(i_var)%num_pos
                  sp_organ_turnover = sp_organ_turnover + this%variables(i_var)%turnover(i_pos)
               end do
            end if
            
         end do

      end if
      
      return
    end function GetTurnover

    ! =========================================================================
    
    function GetBurned(this, organ_id, species_id, position_id) result(sp_organ_burned)

      ! THis function is very similar to GetTurnover, with the only difference that it
      ! returns the turnover mass so-far during the period of interest.

      class(prt_vartypes)         :: this
      integer,intent(in)          :: organ_id
      integer,intent(in)          :: species_id
      integer,intent(in),optional :: position_id
      real(r8)                    :: sp_organ_burned

      integer :: i_pos
      integer :: ispec
      integer :: num_species
      integer,dimension(max_spec_per_group) :: spec_ids 
      integer :: i_var
      
      sp_organ_burned = 0.0_r8
      
      if(species_id == all_carbon_species) then
         spec_ids(1:3) = carbon_species(1:3)
         num_species  = 3
      else
         num_species  = 1
         spec_ids(1) = species_id
      end if

      if(present(position_id)) then
         i_pos = position_id
      
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) sp_organ_burned = sp_organ_burned + &
                  this%variables(i_var)%burned(i_pos)
         end do

      else
         
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) then
               do i_pos = 1, this%variables(i_var)%num_pos
                  sp_organ_burned = sp_organ_burned + this%variables(i_var)%burned(i_pos)
               end do
            end if
            
         end do

      end if
      
      return
   end function GetBurned

   ! ====================================================================================
   
   function GetNetART(this, organ_id, species_id, position_id) result(sp_organ_netart)
      
      ! THis function is very similar to GetTurnover, with the only difference that it
      ! returns the Net changes due to Allocations Reactions and Transport in that pool

      class(prt_vartypes)         :: this
      integer,intent(in)          :: organ_id
      integer,intent(in)          :: species_id
      integer,intent(in),optional :: position_id
      real(r8)                    :: sp_organ_netart

      integer :: i_pos
      integer :: ispec
      integer :: num_species
      integer,dimension(max_spec_per_group) :: spec_ids 
      integer :: i_var
      
      sp_organ_netart = 0.0_r8
      
      if(species_id == all_carbon_species) then
         spec_ids(1:3) = carbon_species(1:3)
         num_species  = 3
      else
         num_species  = 1
         spec_ids(1) = species_id
      end if

      if(present(position_id)) then
         i_pos = position_id
      
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) sp_organ_netart = sp_organ_netart + &
                  this%variables(i_var)%net_art(i_pos)
         end do

      else
         
         do ispec = 1,num_species
            i_var = this%prt_instance%sp_organ_map(organ_id,spec_ids(ispec))
            if(i_var>0) then
               do i_pos = 1, this%variables(i_var)%num_pos
                  sp_organ_netart = sp_organ_netart + this%variables(i_var)%net_art(i_pos)
               end do
            end if
            
         end do

      end if
      
      return
   end function GetNetART

   ! =====================================================================================

   function GetCoordVal(this, organ_id, species_id ) result(prt_val)
      
      class(prt_vartypes)               :: this
      integer,intent(in)                :: organ_id
      integer,intent(in)                :: species_id
      real(r8)                          :: prt_val 
      
      write(fates_log(),*)'Init must be extended by a child class.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
      
    end function GetCoordVal


   ! ====================================================================================
   
   subroutine InitAllocateBase(this)

      class(prt_vartypes) :: this
      
      write(fates_log(),*)'Init must be extended by a child class.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
      
    end subroutine InitAllocateBase

   ! ====================================================================================

   subroutine DailyPRTBase(this)
      
      class(prt_vartypes) :: this
      
      write(fates_log(),*)'Daily PRT Allocation must be extended'
      call endrun(msg=errMsg(__FILE__, __LINE__))
   
   end subroutine DailyPRTBase

   ! ====================================================================================
   
   subroutine FastPRTBase(this)

      class(prt_vartypes) :: this
      
      write(fates_log(),*)'FastReactiveTransport must be extended by a child class.'
      call endrun(msg=errMsg(__FILE__, __LINE__))

   end subroutine FastPRTBase

   ! ====================================================================================

   
   subroutine SetState(prt,organ_id, species_id, state_val, position_id)

     ! CONSIDER INTERFACING THIS AND CALLING DIFFERENT SUBROUTINES BY POINTER

     class(prt_vartypes) :: prt
     integer,intent(in)  :: organ_id
     integer,intent(in)  :: species_id
     real(r8),intent(in) :: state_val
     integer,intent(in),optional :: position_id

     
     integer :: ispec
     integer :: n_vars
     integer,dimension(max_spec_per_group) :: spec_ids 
     integer :: i_var
     integer :: i_pos
     
     if(species_id == all_carbon_species) then
        write(fates_log(),*) 'You cannot set the state of all isotopes simultaneously.'
        write(fates_log(),*) 'You can only set 1. Exiting.'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if
     
     if( present(position_id) ) then
        i_pos = position_id
     else
        i_pos = 1
     end if
     
     
     i_var = prt%prt_instance%sp_organ_map(organ_id,species_id)
     
     if(i_pos>prt%variables(i_var)%num_pos)then
        write(fates_log(),*) 'A position index was specified that is'
        write(fates_log(),*) 'greater than the allocated position space'
        write(fates_log(),*) ' i_pos: ',i_pos
        write(fates_log(),*) ' num_pos: ',prt%variables(i_var)%num_pos
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if


     if(i_var>0) then
        prt%variables(i_var)%val(i_pos) = state_val
     else
        write(fates_log(),*) 'A mass was sent to PARTEH to over-write'
        write(fates_log(),*) ' a pool with a specie x organ combination. '
        write(fates_log(),*) ' that does not exist.'
        write(fates_log(),*) ' organ_id:',organ_id
        write(fates_log(),*) ' species_id:',species_id
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if
        

     return
   end subroutine SetState

   ! ====================================================================================


   



end module PRTGenericMod
