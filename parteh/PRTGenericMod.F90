module PRTGenericMod

  ! ------------------------------------------------------------------------------------
  ! Plant Allocation and Reactive Transport (PART) +
  ! Extensible Hypotheses (EH) = PARTEH
  ! 
  ! Non-Specific (Generic) Classes and Functions
  ! This contains the base classes for both the variables and the global class
  !
  ! General idea:  PARTEH treats its state variables as objects. Each object
  !  can be mapped to, or associated with:
  !  1) an organ
  !  2) a spatial position associated with that organ
  !  3) a chemical element (ie carbon isotope or nutrient), aka chemical species
  !
  !
  ! THIS ROUTINE SHOULD NOT HAVE TO BE MODIFIED TO ACCOMODATE NEW HYPOTHESES
  ! (in principle ...)
  !
  ! Ryan Knox, April 2018
  ! ------------------------------------------------------------------------------------


  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : calloc_abs_error
  use FatesConstantsMod, only : years_per_day
  use FatesConstantsMod, only : days_per_sec
  use FatesGlobals     , only : endrun => fates_endrun
  use FatesGlobals     , only : fates_log 
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use PRTParametersMod , only : prt_params
  
  implicit none
  private ! Modules are private by default

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  integer, parameter, public :: maxlen_varname   = 128
  integer, parameter, public :: maxlen_varsymbol = 32
  integer, parameter, public :: maxlen_varunits  = 32
  integer, parameter, public :: len_baseunit     = 6


  ! We use this parameter as the value for which we set un-initialized values
  real(r8), parameter, public :: un_initialized = -9.9e32_r8

  ! We use this parameter as the value for which we check un-initialized values
  real(r8), parameter, public :: check_initialized = -8.8e32_r8


  ! -------------------------------------------------------------------------------------
  ! IMPORTANT!
  ! All elements in all organs should be expressed in terms of KILOGRAMS
  ! All rates of change are expressed in terms of kilograms / day
  ! This assumption cannot be broken!
  ! -------------------------------------------------------------------------------------
  
  character(len=len_baseunit), parameter, public :: mass_unit = 'kg'
  character(len=len_baseunit), parameter, public :: mass_rate_unit = 'kg/day'

  ! -------------------------------------------------------------------------------------
  ! Allocation Hypothesis Types
  ! These should each have their own module
  ! -------------------------------------------------------------------------------------

  integer, parameter, public :: prt_carbon_allom_hyp   = 1
  integer, parameter, public :: prt_cnp_flex_allom_hyp = 2

  ! -------------------------------------------------------------------------------------
  ! Organ types
  ! These are public indices used to map the organs
  ! in each hypothesis to organs that acknowledged in the calling model
  ! -------------------------------------------------------------------------------------

  integer, parameter, public :: num_organ_types = 6
  integer, parameter, public :: all_organs    = 0    ! index for all organs
  integer, parameter, public :: leaf_organ    = 1    ! index for leaf organs
  integer, parameter, public :: fnrt_organ    = 2    ! index for fine-root organs
  integer, parameter, public :: sapw_organ    = 3    ! index for sapwood organs
  integer, parameter, public :: store_organ   = 4    ! index for storage organs
  integer, parameter, public :: repro_organ   = 5    ! index for reproductive organs
  integer, parameter, public :: struct_organ  = 6    ! index for structure (dead) organs

  ! -------------------------------------------------------------------------------------
  ! Element types
  ! These are public indices used to map the elements (chem species) in each hypothesis
  ! to the element that are acknowledged in the calling model
  ! -------------------------------------------------------------------------------------

  integer, parameter, public :: num_element_types     = 6    ! Total number of unique element
                                                     ! curently recognized by PARTEH
                                                     ! should be max index in list below

  ! The following list are the unique indices associated with the
  ! element used in each hypothesis.  Note these are just POTENTIAL
  ! element.  At the time of writing this, we are very far away from
  ! creating allocation schemes that even use potassium.
  
  integer, parameter, public :: carbon12_element    = 1
  integer, parameter, public :: carbon13_element    = 2
  integer, parameter, public :: carbon14_element    = 3
  integer, parameter, public :: nitrogen_element    = 4
  integer, parameter, public :: phosphorus_element  = 5
  integer, parameter, public :: potassium_element   = 6


  !  The following elements are just placeholders. In the future
  !  if someone wants to develope an allocation hypothesis
  !  that uses nickel, we can just uncomment it from this list

  !  integer, parameter, public :: calcium_element     = 7
  !  integer, parameter, public :: magnesium_element   = 8
  !  integer, parameter, public :: sulfur_element      = 9
  !  integer, parameter, public :: chlorine_element    = 10
  !  integer, parameter, public :: iron_element        = 11
  !  integer, parameter, public :: manganese_element   = 12
  !  integer, parameter, public :: zinc_element        = 13
  !  integer, parameter, public :: copper_element      = 14
  !  integer, parameter, public :: boron_element       = 15
  !  integer, parameter, public :: molybdenum_element  = 16
  !  integer, parameter, public :: nickel_element      = 17

  
  ! We have some lists of elements or lists of organs, such as
  ! a list of all carbon elements.  To keep routines simple
  ! we set a global to the maximum list size for scratch arrays.

  integer, parameter, public :: max_spec_per_group = 3    ! we may query these lists
                                                  ! the carbon elements are the biggest list
                                                  ! right now


  ! List of all carbon elements, the special index "all_carbon_elements"
  ! implies the following list of carbon organs (NOT USED)
  
  integer, parameter, dimension(3), public :: carbon_elements_list   = &
       [carbon12_element, carbon13_element, carbon14_element]



  ! This is the maximum number of leaf age pools allowed on each plant
  ! (used for allocating scratch space)
  integer, parameter, public :: max_nleafage = 4


  ! This is the minimum allowable L2FR, this is needed so that plants
  ! in the understory don't shrink their roots down so far that
  ! they dissappear and cause numerical issues
  real(r8), parameter, public :: l2fr_min = 0.01_r8
  
  ! -------------------------------------------------------------------------------------
  !
  ! The following is the data structure that holds the state (ie carbon,
  ! nutrients, etc) for each pool of each plant.
  !
  ! For example, this could be the carbon 12 of the leaf pool; its instantaneous state,
  !   and its fluxes.
  !
  ! Note also that these are vectors and not scalars, which indicates that there
  ! may be more than 1 discrete spatial positions.  For instance, there might be multiple
  ! leaf layers or something.
  ! 
  ! Since there are many variables, as well as boundary conditions, this object is 
  ! NESTED in the prt_vartypes  (<---- see the "s" at the end?)  structure that follows.
  ! 
  ! Each object will have a unique index associated with it, it will also be mapped
  ! to a specific organ and element combination.
  ! 
  ! It is assumed that over the control period (probably 1 day) that
  ! changes in the current state (val) relative to the value at the start of the
  ! control period (val0), are equal to the time integrated flux terms 
  ! (net_alloc, turnover, etc)
  !
  ! -------------------------------------------------------------------------------------

  type, public :: prt_vartype
     
     real(r8),pointer :: val(:)       ! Instantaneous state variable           [kg]
     real(r8),allocatable :: val0(:)      ! State variable at the beginning 
                                          ! of the control period                  [kg]
     real(r8),allocatable :: net_alloc(:)   ! Net change due to allocation/transport [kg]
                                          ! over the control period                [kg]
     real(r8),allocatable :: turnover(:)  ! Losses due to turnover                 [kg]
                                          ! or, any mass destined for litter
                                          ! over the control period

     real(r8),allocatable :: burned(:)    ! Losses due to burn                     [kg]
     real(r8),allocatable :: damaged(:)   ! Losses due to damage                   [kg]
     real(r8),allocatable :: herbivory(:) ! Losses due to herbivory                [kg]

     ! Placeholder
     ! To save on memory, keep this commented out, or simply
     ! add this only in the extension ... ?
     ! real(r8),allocatable           :: coordinate(:,:)

  end type prt_vartype


  ! -------------------------------------------------------------------------------------
  ! Input boundary conditions.  These will be allocated as an array for each plant.
  ! This type will also be broken into 3 types of boundary conditions:  input only,
  ! output only, and input-output.
  ! -------------------------------------------------------------------------------------

  type, public :: prt_bctype
     
     real(r8), pointer :: rval
     integer, pointer  :: ival

  end type prt_bctype


  ! -------------------------------------------------------------------------------------
  ! The following is the object that is directly attached to each plant.
  !
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

  type, public :: prt_vartypes
     
     type(prt_vartype),allocatable :: variables(:)    ! The state variables and fluxes
     type(prt_bctype), allocatable :: bc_inout(:)     ! These boundaries may be changed
     type(prt_bctype), allocatable :: bc_in(:)        ! These are protected
     type(prt_bctype), allocatable :: bc_out(:)       ! These are overwritten
     real(r8)                      :: ode_opt_step

   contains
     
     ! These are extendable procedures that have specialized
     ! content in each of the different hypotheses
     
     procedure :: DailyPRT            => DailyPRTBase
     procedure :: FastPRT             => FastPRTBase
     procedure :: DamageRecovery      => DamageRecoveryBase
     procedure :: GetNutrientTarget   => GetNutrientTargetBase
     
     ! These are generic functions that should work on all hypotheses

     procedure, non_overridable :: InitAllocate
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
     procedure, non_overridable :: GetHerbivory
     procedure, non_overridable :: GetNetAlloc
     procedure, non_overridable :: ZeroRates
     procedure, non_overridable :: CheckMassConservation
     procedure, non_overridable :: DeallocatePRTVartypes
     procedure, non_overridable :: WeightedFusePRTVartypes
     procedure, non_overridable :: CopyPRTVartypes
     
     
     procedure :: AgeLeaves  ! This routine may be used generically
                             ! but also leaving the door open for over-rides
     

     
  end type prt_vartypes

  
  ! Global identifiers for which elements we are using (apply mostly to litter)

  integer, public              :: num_elements          ! This is the number of elements in this simulation
                                                        ! e.g. (C,N,P,K, etc)
  integer, allocatable, public :: element_list(:)       ! This vector holds the list of global element identifiers
                                                        ! examples are carbon12_element
                                                        ! nitrogen_element, etc.

  integer, public :: element_pos(num_element_types)     ! This is the reverse lookup
                                                        ! for element types. Pick an element
                                                        ! global index, and it gives you
                                                        ! the position in the element_list

  ! -------------------------------------------------------------------------------------
  ! This next section contains the objects that describe the mapping for each specific
  ! hypothesis. It is also a way to call the descriptions of variables for any
  ! arbitrary hypothesis.
  ! These are things that are globally true, not specific to each plant.
  ! For instance the map just contains the list of variable names, not the values for
  ! each plant.
  ! These are not instanced on every plant, they are just instanced once on every model 
  ! machine or memory space. They should only be initialized once and used
  ! as read-only from that point on.
  ! -------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------
  ! This type simply packs the names and symbols associated with all
  ! the variables for any given hypothesis
  ! -------------------------------------------------------------------------------------
  
  type, public :: state_descriptor_type
     character(len=maxlen_varname)   :: longname
     character(len=maxlen_varsymbol) :: symbol
     integer                         :: organ_id    ! global id for organ
     integer                         :: element_id  ! global id for element
     integer                         :: num_pos     ! number of descrete spatial positions

     ! Also, will probably need flags to define different types of groups that this variable
     ! belongs too, which will control things like fusion, normalization, when to zero, etc...

  end type state_descriptor_type
  


  ! This type will help us loop through all the different variables associated
  ! with a specific organ type. Since variables are a combination of organ and
  ! element, the number of unique variables is capped at the number of elements
  ! per each organ.
  
  type, public :: organ_map_type
     integer, dimension(1:num_element_types) :: var_id
     integer                                 :: num_vars
  end type organ_map_type


  ! This structure packs both the mapping structure and the variable descriptors
  ! -------------------------------------------------------------------------------------
  ! This array should contain the lists of indices to 
  ! the element x organ variable structure that is used to map variables to the outside
  ! world.
  !   
  !              
  !                   | carbon | nitrogen | phosphorus  | .... |
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

  type, public :: prt_global_type
     
     ! Note that index 0 is reserved for "all" or "irrelevant"
     character(len=maxlen_varname)                             :: hyp_name

     ! This is the hypothesis index, used internally for some non-specific
     ! routines where different methods are applied
     integer                                                   :: hyp_id

     ! This will save the specific variable id associated with each organ and element
     integer, dimension(0:num_organ_types,0:num_element_types) :: sp_organ_map

     ! This holds the verbose descriptions of the variables, symbols, names, etc
     type(state_descriptor_type), allocatable                  :: state_descriptor(:)

     ! This will save the list of variable ids associated with each organ. There
     ! are multiple of these because we may have multiple element per organ.
     type(organ_map_type), dimension(1:num_organ_types)        :: organ_map

     ! The number of input boundary conditions
     integer                                                   :: num_bc_in      

     ! The number of output boundary conditions                
     integer                                                   :: num_bc_out
     
     ! The number of combo input-output boundary conditions
     integer                                                   :: num_bc_inout
     
     ! The number of variables set by each hypothesis
     integer                                                   :: num_vars
     

  contains
        
     procedure, non_overridable :: ZeroGlobal
     procedure, non_overridable :: RegisterVarInGlobal

  end type prt_global_type

  
  class(prt_global_type),pointer,public :: prt_global

  ! Make necessary procedures public
  public :: GetCoordVal
  public :: SetState
  public :: StorageNutrientTarget
  
contains

  ! =====================================================================================
  ! Module Functions and Subroutines
  ! =====================================================================================


  subroutine ZeroGlobal(this)
  

     ! This subroutine zero's out the map between variable indexes and the
     ! elements and organs they are associated with.
     ! It also sets the counts of the variables and boundary conditions as 
     ! a nonsense number that will trigger a fail if they are specified later.
     ! This routine must be called 


     class(prt_global_type)    :: this
     
     integer :: io ! Organ loop counter
     integer :: is ! Element loop counter
     
     ! First zero out the array
     do io = 1,num_organ_types
        do is = 1,num_element_types
           this%sp_organ_map(io,is)      = 0
           this%organ_map(io)%var_id(is) = 0
        end do
        this%organ_map(io)%num_vars      = 0
     end do
     
     ! Set the number of boundary conditions as a bogus value
     this%num_bc_in     = -9
     this%num_bc_out    = -9
     this%num_bc_inout  = -9

     ! Set the number of variables to a bogus value. This should be
     ! immediately over-written in the routine that is calling this
     this%num_vars = -9

     return
  end subroutine ZeroGlobal
   
  ! =====================================================================================
  
  subroutine RegisterVarInGlobal(this, var_id, long_name, symbol, organ_id, element_id, num_pos)

     
     ! This subroutine is called for each variable that is defined in each specific hypothesis.
     ! For instance, this is called six times in the carbon only hypothesis,
     ! each time providing names, symbols, associated organs and element for each pool.
     
     class(prt_global_type)      :: this
     integer, intent(in)         :: var_id
     character(len=*),intent(in) :: long_name
     character(len=*),intent(in) :: symbol
     integer, intent(in)         :: organ_id
     integer, intent(in)         :: element_id
     integer, intent(in)         :: num_pos

     ! Set the descriptions and the associated organs/element in the variable's
     ! own array

     this%state_descriptor(var_id)%longname   = long_name
     this%state_descriptor(var_id)%symbol     = symbol
     this%state_descriptor(var_id)%organ_id   = organ_id
     this%state_descriptor(var_id)%element_id = element_id
     this%state_descriptor(var_id)%num_pos    = num_pos

     ! Set the mapping tables for the external model

     this%sp_organ_map(organ_id,element_id) = var_id

     ! Set another map that helps to locate all the relevant pools associated
     ! with an organ

     this%organ_map(organ_id)%num_vars = this%organ_map(organ_id)%num_vars + 1
     this%organ_map(organ_id)%var_id(this%organ_map(organ_id)%num_vars) = var_id

        
     return
  end subroutine RegisterVarInGlobal

  ! =====================================================================================

  subroutine InitPRTVartype(this)

    class(prt_vartypes) :: this
    
    
    ! This subroutine should be the first call whenever a prt_vartype object is
    ! instantiated.  
    ! 
    ! Most likely, this will occur whenever a new plant or cohort is created.
    ! 
    ! This routine handles the allocation (extended procedure)
    ! and then the initializing of states with bogus information, and then
    ! the flushing of all boundary conditions to null.

    call this%InitAllocate()                    ! Allocate memory spaces
    call this%InitializeInitialConditions()     ! Set states to a nan-like starter value
    call this%FlushBCs()                        ! Set all boundary condition pointers 
                                                ! to null


    return
  end subroutine InitPRTVartype

  ! =====================================================================================
  
  subroutine InitAllocate(this)
    
     ! ----------------------------------------------------------------------------------
     ! This initialization is called everytime a plant/cohort
     ! is newly recruited.  Like the name implies, we are just allocating space here.
     ! ----------------------------------------------------------------------------------

     class(prt_vartypes) :: this        

     integer :: i_var    ! Variable loop index
     integer :: num_pos  ! The number of positions for each variable

     ! Allocate the boundar condition arrays and flush them to no-data flags
     ! ----------------------------------------------------------------------------------

     if(prt_global%num_bc_in > 0) then
        allocate(this%bc_in(prt_global%num_bc_in))
     end if

     if(prt_global%num_bc_inout > 0) then
        allocate(this%bc_inout(prt_global%num_bc_inout))
     end if

     if(prt_global%num_bc_out > 0) then
        allocate(this%bc_out(prt_global%num_bc_out))
     end if
     
     ! Allocate the state variables
     allocate(this%variables(prt_global%num_vars))
     
     do i_var = 1, prt_global%num_vars
        
        num_pos = prt_global%state_descriptor(i_var)%num_pos 
        
        allocate(this%variables(i_var)%val(num_pos))
        allocate(this%variables(i_var)%val0(num_pos))
        allocate(this%variables(i_var)%turnover(num_pos))
        allocate(this%variables(i_var)%net_alloc(num_pos))
        allocate(this%variables(i_var)%burned(num_pos))
        allocate(this%variables(i_var)%damaged(num_pos))
        allocate(this%variables(i_var)%herbivory(num_pos))
        
     end do

     
     return
  end subroutine InitAllocate

  ! =====================================================================================

  subroutine InitializeInitialConditions(this)

     ! ----------------------------------------------------------------------------------
     ! This routine sets all PARTEH variables to a nonsense value.
     ! This ensures that a fail is triggered if a value is not initialized correctly.
     ! ----------------------------------------------------------------------------------

    class(prt_vartypes) :: this

    integer :: i_var      ! Variable index

    do i_var = 1, prt_global%num_vars
       this%variables(i_var)%val(:)       = un_initialized
       this%variables(i_var)%val0(:)      = un_initialized
       this%variables(i_var)%turnover(:)  = un_initialized
       this%variables(i_var)%burned(:)    = un_initialized
       this%variables(i_var)%damaged(:)   = un_initialized
       this%variables(i_var)%net_alloc(:) = un_initialized
       this%variables(i_var)%herbivory(:) = un_initialized
    end do

    ! Initialize the optimum step size as very large.
    
    this%ode_opt_step = 1e6_r8
    
    return
  end subroutine InitializeInitialConditions


  ! =============================================================

  subroutine CheckInitialConditions(this)
    
    ! This subroutine makes sure that every variable defined
    ! in the hypothesis has been given an initial value.
    !
    ! This should be called following any blocks where initial
    ! conditions are set. In fates, these calls already
    ! exist and when new hypotheses are added, they will
    ! already be checked if the initial conditions are 
    ! specified in parallel with the other hypotheses.

    class(prt_vartypes) :: this

    integer :: i_var      ! index for iterating variables
    integer :: n_cor_ids  ! Number of coordinate ids
    integer :: i_cor      ! index for iterating coordinate dimension
    integer :: i_organ   ! The global organ id for this variable
    integer :: i_element ! The global element id for this variable

    do i_var = 1, prt_global%num_vars

       n_cor_ids = size(this%variables(i_var)%val,1)

       do i_cor = 1, n_cor_ids
       
          if(this%variables(i_var)%val(i_cor) < check_initialized) then

             i_organ   = prt_global%state_descriptor(i_var)%organ_id
             i_element = prt_global%state_descriptor(i_var)%element_id 
             write(fates_log(),*)'Not all initial conditions for state variables'
             write(fates_log(),*)' in PRT hypothesis: ',trim(prt_global%hyp_name)
             write(fates_log(),*)' were written out.'
             write(fates_log(),*)' i_var: ',i_var
             write(fates_log(),*)' i_cor: ',i_cor
             write(fates_log(),*)' organ_id:',i_organ
             write(fates_log(),*)' element_id',i_element
             write(fates_log(),*)'Exiting'
             call endrun(msg=errMsg(sourcefile, __LINE__))
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

    ! This routine must be called once for each "input only" boundary condition of each 
    ! hypothesis.  
    ! The group of calls only needs to happen once, following InitPRTVartype.
    ! Since we use pointers, we don't need to constantly ask for new boundary conditions
    !
    ! The only complication to this would occur, if the boundary condition variable
    ! that these pointers point to is being disassociated. In that case, one would
    ! need to re-register that boundary condition variable.

    
    ! Input Arguments
    
    class(prt_vartypes)                      :: this
    integer,intent(in)                       :: bc_id
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

  subroutine RegisterBCOut(this,bc_id, bc_rval, bc_ival )

    
    ! This routine is similar to the routine above RegisterBCIn, except this 
    ! is for registering "output only" boundary conditions.


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


    ! This routine is similar to the two routines above, except this 
    ! is for registering "input-output" boundary conditions.
    ! These are conditions that are passed into PARTEH, and are expected
    ! to be updated (or not), and passed back to the host (FATES).

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


  subroutine CopyPRTVartypes(this, donor_prt_obj)

    ! Here we copy over all information from a donor_prt_object into the current PRT
    ! object.   It is assumed that the current PRT object
    ! has already been initialized ( ie. InitAllocate() )
    ! variable val0 is omitted, because it is ephemeral and used only during the
    ! allocation process

    ! Arguments
    class(prt_vartypes)                      :: this
    class(prt_vartypes), intent(in), pointer :: donor_prt_obj

    ! Locals

    integer :: i_var    ! loop iterator for variable objects
    integer :: i_bc     ! loop iterator for boundary pointers

    integer :: num_bc_in
    integer :: num_bc_inout
    integer :: num_bc_out

    do i_var = 1, prt_global%num_vars
       this%variables(i_var)%val(:)       = donor_prt_obj%variables(i_var)%val(:)
       this%variables(i_var)%val0(:)      = donor_prt_obj%variables(i_var)%val0(:)
       this%variables(i_var)%net_alloc(:)   = donor_prt_obj%variables(i_var)%net_alloc(:)
       this%variables(i_var)%turnover(:)  = donor_prt_obj%variables(i_var)%turnover(:)
       this%variables(i_var)%burned(:)    = donor_prt_obj%variables(i_var)%burned(:)
       this%variables(i_var)%damaged(:)   = donor_prt_obj%variables(i_var)%damaged(:)
       this%variables(i_var)%herbivory(:) = donor_prt_obj%variables(i_var)%herbivory(:)
    end do

    this%ode_opt_step = donor_prt_obj%ode_opt_step

    return
  end subroutine CopyPRTVartypes


  ! =====================================================================================

  subroutine WeightedFusePRTVartypes(this,donor_prt_obj, recipient_fuse_weight)
    
    ! This subroutine fuses two PRT objects together based on a fusion weighting
    ! assigned for the recipient (the class calling this)

    ! Arguments 
    class(prt_vartypes)                      :: this
    class(prt_vartypes), intent(in), pointer :: donor_prt_obj
    real(r8),intent(in)                      :: recipient_fuse_weight   ! This is the weighting
                                                                        ! for the recipient
!    integer,intent(in),optional              :: position_id

    ! Locals
    integer :: i_var    ! Loop iterator over variables
    integer :: pos_id   ! coordinate id (defaults to 1, if not position_id)

    ! If no argument regarding positions is supplied, then fuse all positions
    ! sequentially
    do i_var = 1, prt_global%num_vars
       
       do pos_id = 1,prt_global%state_descriptor(i_var)%num_pos
          
          this%variables(i_var)%val(pos_id)  = recipient_fuse_weight * this%variables(i_var)%val(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%val(pos_id)
          
          this%variables(i_var)%val0(pos_id)  = recipient_fuse_weight * this%variables(i_var)%val0(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%val0(pos_id)
          
          this%variables(i_var)%net_alloc(pos_id)      = recipient_fuse_weight * this%variables(i_var)%net_alloc(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%net_alloc(pos_id)
          
          this%variables(i_var)%turnover(pos_id)    = recipient_fuse_weight * this%variables(i_var)%turnover(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%turnover(pos_id)
          
          this%variables(i_var)%burned(pos_id)    = recipient_fuse_weight * this%variables(i_var)%burned(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%burned(pos_id)
          
          this%variables(i_var)%damaged(pos_id)    = recipient_fuse_weight * this%variables(i_var)%damaged(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%damaged(pos_id)
          
          this%variables(i_var)%herbivory(pos_id)    = recipient_fuse_weight * this%variables(i_var)%herbivory(pos_id) + &
                (1.0_r8-recipient_fuse_weight) * donor_prt_obj%variables(i_var)%herbivory(pos_id)
       end do
    end do


    this%ode_opt_step = recipient_fuse_weight * this%ode_opt_step + &
                        (1.0_r8-recipient_fuse_weight) * donor_prt_obj%ode_opt_step


    return
  end subroutine WeightedFusePRTVartypes

  ! =====================================================================================

  subroutine DeallocatePRTVartypes(this)

    ! ---------------------------------------------------------------------------------
    ! Unfortunately ... all plants must die. It is sad, but when this happens
    ! we must also deallocate our memory of them.  Man, thats really is sad. Why
    ! must we also forget them...  Well, anyway,  any time a plant/cohort
    ! is deallocated, we must also deallocate all this memory bound in the PARTEH
    ! data structure.  But on the bright side, there will always be new recruits,
    ! a new generation, to allocate as well.  Life must go on.
    ! I suppose since we are recording their life in the history output, in a way
    ! we are remembering them. I feel better now.
    ! ---------------------------------------------------------------------------------

    class(prt_vartypes) :: this
    integer             :: i_var, istat
    character(len=255)  :: smsg
    
    ! Check to see if there is any value in these pools?
    ! SHould not deallocate if there is any carbon left

    if(allocated(this%variables)) then
       do i_var = 1, prt_global%num_vars
          deallocate( &
               this%variables(i_var)%val, &
               this%variables(i_var)%val0, &
               this%variables(i_var)%net_alloc, &
               this%variables(i_var)%turnover, &
               this%variables(i_var)%burned, &
               this%variables(i_var)%herbivory, &               
               stat=istat, errmsg=smsg )
          if (istat/=0) call endrun(msg='DeallocatePRTVartypes 1 stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
       end do
   
       deallocate(this%variables, stat=istat, errmsg=smsg)
       if (istat/=0) call endrun(msg='DeallocatePRTVartypes 2 stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
    end if

    if(allocated(this%bc_in))then
       deallocate(this%bc_in, stat=istat, errmsg=smsg)
       if (istat/=0) call endrun(msg='DeallocatePRTVartypes bc_in stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
    end if
    
    if(allocated(this%bc_out))then
       deallocate(this%bc_out, stat=istat, errmsg=smsg)
       if (istat/=0) call endrun(msg='DeallocatePRTVartypes bc_out stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
    end if

    if(allocated(this%bc_inout))then
       deallocate(this%bc_inout, stat=istat, errmsg=smsg)
       if (istat/=0) call endrun(msg='DeallocatePRTVartypes bc_inout stat/=0:'//trim(smsg)//errMsg(sourcefile, __LINE__))
    end if

    return
  end subroutine DeallocatePRTVartypes
  
  ! =====================================================================================

  subroutine ZeroRates(this)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine zeros all of the rates of change for our variables.
      ! It also sets the initial value to the current state.
      ! This allows us to make mass conservation checks, where
      ! val - val0 = net_alloc + turnover
      ! 
      ! This subroutine is called each day in FATES, which is the control interval
      ! that we conserve carbon from the allocation and turnover process.
      ! ---------------------------------------------------------------------------------

      class(prt_vartypes) :: this

      integer :: i_var    ! Variable index

      do i_var = 1, prt_global%num_vars
         this%variables(i_var)%val0(:)      = this%variables(i_var)%val(:)
         this%variables(i_var)%net_alloc(:) = 0.0_r8
         this%variables(i_var)%turnover(:)  = 0.0_r8
         this%variables(i_var)%burned(:)    = 0.0_r8
         this%variables(i_var)%herbivory(:)    = 0.0_r8
         this%variables(i_var)%damaged(:)   = 0.0_r8    
      end do
      
    end subroutine ZeroRates

   ! ====================================================================================

   subroutine CheckMassConservation(this,ipft,position_id)

      
     ! ---------------------------------------------------------------------------------
     ! At any time, the sum of fluxes should equal the difference between val and val0.
     ! This routine loops over all variables and ensures this is true.
     ! The final argument is any uniqely identifying index that can be used
     ! to differentiate where in the call sequence a failure in conservation occurs.
     ! ---------------------------------------------------------------------------------

     class(prt_vartypes) :: this    
     integer, intent(in) :: ipft        ! functional type of the plant
     integer, intent(in) :: position_id ! Helps to know where
                                        ! in the call sequence this was called

     integer :: i_var       ! Variable index
     integer :: i_pos       ! Position (coordinate) index

     real(r8) :: err        ! absolute error [kg]
     real(r8) :: rel_err    ! error relative to the pool's size [kg]


     do i_var = 1, prt_global%num_vars

        do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
           
           err = abs((this%variables(i_var)%val(i_pos) - this%variables(i_var)%val0(i_pos)) - &
                  (this%variables(i_var)%net_alloc(i_pos) &
                   -this%variables(i_var)%turnover(i_pos) & 
                   -this%variables(i_var)%burned(i_pos) &
                   -this%variables(i_var)%herbivory(i_pos) &
                   -this%variables(i_var)%damaged(i_pos)))
                  
           if(this%variables(i_var)%val(i_pos) > nearzero ) then
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
              write(fates_log(),*) ' organ id: ',prt_global%state_descriptor(i_var)%organ_id
              write(fates_log(),*) ' element_id: ',prt_global%state_descriptor(i_var)%element_id
              write(fates_log(),*) ' position id: ',i_pos
              write(fates_log(),*) ' symbol: ',trim(prt_global%state_descriptor(i_var)%symbol)
              write(fates_log(),*) ' longname: ',trim(prt_global%state_descriptor(i_var)%longname)
              write(fates_log(),*) ' err: ',err,' max error: ',calloc_abs_error
              write(fates_log(),*) ' terms: ', this%variables(i_var)%val(i_pos), &
                                               this%variables(i_var)%val0(i_pos), &
                                               this%variables(i_var)%net_alloc(i_pos), &
                                               this%variables(i_var)%turnover(i_pos), &
                                               this%variables(i_var)%burned(i_pos), &
                                               this%variables(i_var)%herbivory(i_pos), & 
                                               this%variables(i_var)%damaged(i_pos)
              write(fates_log(),*) ' Exiting.'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if

        end do
     end do

     return
   end subroutine CheckMassConservation

   ! ====================================================================================
   
   function GetState(this, organ_id, element_id, position_id) result(state_val)
      
      ! This function returns the current amount of mass for
      ! any combination of organ and element. **IF** a position
      ! is provided, it will use it, but otherwise, it will sum over
      ! all dimensions.  It also can accomodate all_carbon_element, which
      ! will return the mass of all carbon isotopes combined.
     
      class(prt_vartypes)                   :: this
      integer,intent(in)                    :: organ_id           ! Organ type querried
      integer,intent(in)                    :: element_id         ! Element type querried
      integer,intent(in),optional           :: position_id        ! Position querried
      real(r8)                              :: state_val          ! Mass (value) of state variable [kg]
      integer                               :: i_pos              ! position loop counter
      integer                               :: i_var              ! variable id

      if(present(position_id)) then

         i_pos = position_id
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         state_val = this%variables(i_var)%val(i_pos)
         
      else
         
         state_val = 0._r8
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
            state_val = state_val + this%variables(i_var)%val(i_pos)
         end do
         
      end if
      
      return
    end function GetState

    
   ! ====================================================================================

   
    function GetTurnover(this, organ_id, element_id, position_id) result(turnover_val)
      
     
      ! THis function is very similar to GetState, with the only difference that it
      ! returns the turnover mass so-far during the period of interest.
      !
      ! NOTE: THIS HAS NOTHING TO DO WITH SPECIFYING TURNOVER. THIS IS JUST A QUERY FUNCTION


      class(prt_vartypes)                   :: this
      integer,intent(in)                    :: organ_id           ! Organ type querried
      integer,intent(in)                    :: element_id         ! Element type querried
      integer,intent(in),optional           :: position_id        ! Position querried
      real(r8)                              :: turnover_val       ! Amount (value) of turnover [kg]
      integer                               :: i_pos              ! position loop counter
      integer                               :: i_var              ! variable id
      
      if(present(position_id)) then

         i_pos = position_id
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         turnover_val = this%variables(i_var)%turnover(i_pos)

      else

         turnover_val = 0.0_r8
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
            turnover_val = turnover_val + this%variables(i_var)%turnover(i_pos)
         end do
            

      end if
      
      return
    end function GetTurnover

    ! =========================================================================
    
    function GetBurned(this, organ_id, element_id, position_id) result(burned_val)

      ! THis function is very similar to GetTurnover, with the only difference that it
      ! returns the burned mass so-far during the period of interest.
      
      ! NOTE: THIS HAS NOTHING TO DO WITH SPECIFYING BURNING. THIS IS JUST A QUERY FUNCTION

      class(prt_vartypes)                   :: this
      integer,intent(in)                    :: organ_id           ! Organ type querried
      integer,intent(in)                    :: element_id         ! Element type querried
      integer,intent(in),optional           :: position_id        ! Position querried
      real(r8)                              :: burned_val         ! Amount (value) of burned [kg]
      integer                               :: i_pos              ! position loop counter
      integer                               :: i_var              ! variable id

      if(present(position_id)) then

         i_pos = position_id
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         burned_val = this%variables(i_var)%burned(i_pos)

      else
         
         burned_val = 0.0_r8
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
            burned_val = burned_val + this%variables(i_var)%burned(i_pos)
         end do

      end if
      
      return
   end function GetBurned

    ! =========================================================================
    
    function GetHerbivory(this, organ_id, element_id, position_id) result(herbivory_val)

      ! This function is very very similar to GetBurned
      
      class(prt_vartypes)                   :: this
      integer,intent(in)                    :: organ_id           ! Organ type querried
      integer,intent(in)                    :: element_id         ! Element type querried
      integer,intent(in),optional           :: position_id        ! Position querried
      real(r8)                              :: herbivory_val      ! Amount (value) of herbivory [kg]
      integer                               :: i_pos              ! position loop counter
      integer                               :: i_var              ! variable id

      if(present(position_id)) then

         i_pos = position_id
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         herbivory_val = this%variables(i_var)%herbivory(i_pos)

      else
         
         herbivory_val = 0.0_r8
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
            herbivory_val = herbivory_val + this%variables(i_var)%herbivory(i_pos)
         end do

      end if
      
      return
   end function GetHerbivory

   ! ====================================================================================
   
   function GetNetAlloc(this, organ_id, element_id, position_id) result(val_netalloc)
      
      ! THis function is very similar to GetTurnover, with the only difference that it
      ! returns the Net changes due to Allocations Reactions and Transport in that pool

      ! NOTE: THIS HAS NOTHING TO DO WITH SPECIFYING ALLOCATION/TRANSPORT. 
      ! THIS IS JUST A QUERY FUNCTION

      class(prt_vartypes)                   :: this
      integer,intent(in)                    :: organ_id           ! Organ type querried
      integer,intent(in)                    :: element_id         ! Element type querried
      integer,intent(in),optional           :: position_id        ! Position querried
      real(r8)                              :: val_netalloc       ! Amount (value) of allocation [kg]
      integer                               :: i_pos              ! position loop counter
      integer                               :: i_var              ! variable id
      
      if(present(position_id)) then

         i_pos = position_id
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         val_netalloc = this%variables(i_var)%net_alloc(i_pos)
         
      else

         val_netalloc = 0.0_r8
         i_var = prt_global%sp_organ_map(organ_id,element_id)
         do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos 
            val_netalloc = val_netalloc + this%variables(i_var)%net_alloc(i_pos)
         end do

      end if
      
      return
   end function GetNetAlloc

   ! =====================================================================================

   function GetCoordVal(this, organ_id, element_id ) result(prt_val)
      
      ! This is support code that may be helpful when we have variables in parteh
      ! that have multiple discrete spatial positions.
      

      class(prt_vartypes)               :: this
      integer,intent(in)                :: organ_id
      integer,intent(in)                :: element_id
      real(r8)                          :: prt_val 
      
      prt_val = 0._r8
      write(fates_log(),*)'Init must be extended by a child class.'
      call endrun(msg=errMsg(sourcefile, __LINE__))
      
    end function GetCoordVal

   ! ====================================================================================

   subroutine DailyPRTBase(this,phase)
      
     class(prt_vartypes) :: this
     integer,intent(in)  :: phase  ! We allow this and its children to be broken into phases

     write(fates_log(),*)'Daily PRT Allocation must be extended'
     call endrun(msg=errMsg(sourcefile, __LINE__))
   
   end subroutine DailyPRTBase

   ! ====================================================================================
   
   subroutine DamageRecoveryBase(this)

      class(prt_vartypes) :: this
      
      write(fates_log(),*)'?'
      call endrun(msg=errMsg(sourcefile, __LINE__))

   end subroutine DamageRecoveryBase

   ! ====================================================================================
   ! ====================================================================================

   subroutine FastPRTBase(this)
      
      class(prt_vartypes) :: this
      
      write(fates_log(),*)'FastReactiveTransport must be extended by a child class'
      call endrun(msg=errMsg(sourcefile, __LINE__))
   
   end subroutine FastPRTBase

   !====================================================================================   
   subroutine SetState(prt,organ_id, element_id, state_val, position_id)

     ! This routine should only be called for initalizing the state value
     ! of a plant's pools.  A value is passed in to set the state of 
     ! organ and element couplets, and position id if it is provided.
     ! A select statement will most definitely bracket the call to this
     ! routine.  

     class(prt_vartypes)                   :: prt
     integer,intent(in)                    :: organ_id      ! organ of interest
     integer,intent(in)                    :: element_id    ! element of interest
     real(r8),intent(in)                   :: state_val     ! value to be initialized
     integer,intent(in),optional           :: position_id   ! position of interest
     
     integer                               :: i_element     ! loop counter for elements
     integer                               :: i_var         ! variable loop counter
     integer                               :: i_pos         ! position loop counter
     
     if( present(position_id) ) then
        i_pos = position_id
     else
        i_pos = 1
     end if
     
     i_var = prt_global%sp_organ_map(organ_id,element_id)
     
     if(i_pos > prt_global%state_descriptor(i_var)%num_pos )then
        write(fates_log(),*) 'A position index was specified that is'
        write(fates_log(),*) 'greater than the allocated position space'
        write(fates_log(),*) ' i_pos: ',i_pos
        write(fates_log(),*) ' num_pos: ',prt_global%state_descriptor(i_var)%num_pos
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if


     if(i_var>0) then
        prt%variables(i_var)%val(i_pos) = state_val
     else
        write(fates_log(),*) 'A mass was sent to PARTEH to over-write'
        write(fates_log(),*) ' a pool with a specie x organ combination. '
        write(fates_log(),*) ' that does not exist.'
        write(fates_log(),*) ' organ_id:',organ_id
        write(fates_log(),*) ' element_id:',element_id
        write(fates_log(),*) 'Exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if
        

     return
   end subroutine SetState

   ! ====================================================================================

   subroutine AgeLeaves(this,ipft,icanlayer,period_sec)

     ! -----------------------------------------------------------------------------------
     ! If we have more than one leaf age classification, allow
     ! some leaf biomass to transition to the older classes.
     ! It is assumed this routine is called once per day.
     ! Note that there is NO turnover or loss of mass on the plant in this routine.
     ! We are simply moving portions of leaves from a young bin to the next older, but
     ! we are not moving any mass out of the last (oldest) bin.
     ! -----------------------------------------------------------------------------------

     class(prt_vartypes)              :: this
     integer,intent(in)               :: ipft
     integer,intent(in)               :: icanlayer
     real(r8),intent(in)              :: period_sec  ! Time period over which this routine
                                                     ! is called [seconds] daily=86400
     integer                          :: nleafage
     integer                          :: i_age
     integer                          :: i_var
     integer                          :: el
     integer                          :: element_id
     real(r8)                         :: leaf_age_flux_frac
     real(r8),dimension(max_nleafage) :: leaf_m0
     real(r8)                         :: leaf_long

     do el = 1, num_elements

        element_id = element_list(el)

        ! Global position of leaf variable
        i_var = prt_global%sp_organ_map(leaf_organ,element_id)

        ! Size of the leaf carbon variable (number of age bins)
        nleafage = prt_global%state_descriptor(i_var)%num_pos ! Number of leaf age class

        associate(leaf_m => this%variables(i_var)%val(:))
          
          leaf_m0(1:nleafage) = leaf_m(1:nleafage)
          
          if(nleafage>1) then
             do i_age = 1,nleafage-1
                if (prt_params%leaf_long(ipft,i_age)>nearzero) then

                   if (icanlayer  .eq.  1)  then
                      leaf_long = prt_params%leaf_long(ipft,i_age)
                   else
                      leaf_long = prt_params%leaf_long_ustory(ipft,i_age)
                   end if

                   ! Units: [-] = [sec] * [day/sec] * [years/day] * [1/years]
                   leaf_age_flux_frac = period_sec * days_per_sec * years_per_day / leaf_long
                                   
                   leaf_m(i_age)    = leaf_m(i_age)   - leaf_m0(i_age) * leaf_age_flux_frac
                   leaf_m(i_age+1)  = leaf_m(i_age+1) + leaf_m0(i_age) * leaf_age_flux_frac
                   
                end if
             end do
          end if

          
          ! Update the diagnostic on daily rate of change
          do i_age = 1,nleafage
             this%variables(i_var)%net_alloc(i_age) = &
                  this%variables(i_var)%net_alloc(i_age) + &
                  (leaf_m(i_age) - leaf_m0(i_age))
          end do


        end associate
     end do
     
 end subroutine AgeLeaves


 function GetNutrientTargetBase(this,element_id,organ_id,stoich_mode) result(target_m)
    
    class(prt_vartypes) :: this
    integer, intent(in)           :: element_id
    integer, intent(in)           :: organ_id
    integer, intent(in),optional  :: stoich_mode
    real(r8)                      :: target_m    ! Target amount of nutrient for this organ [kg]

    target_m = 0._r8
    
    write(fates_log(),*)'GetNutrientTargetBase must be extended by a child class.'
    call endrun(msg=errMsg(sourcefile, __LINE__))

    return
  end function GetNutrientTargetBase


  ! ====================================================================================
   
  function StorageNutrientTarget(pft, element_id, leaf_target, fnrt_target, sapw_target, struct_target) result(store_target)

     integer :: pft
     integer :: element_id
     real(r8) :: leaf_target    ! Target leaf nutrient mass [kg]
     real(r8) :: fnrt_target    ! Target fineroot nutrient mass [kg]
     real(r8) :: sapw_target    ! Target sapwood nutrient mass [kg]
     real(r8) :: struct_target  ! Target structural nutrient mass [kg]

     real(r8) :: store_target   ! Output: Target storage nutrient mass [kg]
     
     
     ! -------------------------------------------------------------------------------------
     ! Choice of how nutrient storage target is proportioned to
     !   Each choice makes the nutrient storage proportional the the "in-tissue"
     !   total nitrogen content of 1 or more sets of organs
     ! -------------------------------------------------------------------------------------
     
     integer, parameter :: lfs_store_prop = 1    ! leaf-sapwood proportional storage
     integer, parameter :: lfss_store_prop = 2   ! leaf-fnrt-sapw-struct proportional storage
     integer, parameter :: fnrt_store_prop = 3   ! fineroot proportional storage
     integer, parameter :: cstore_store_prop = 4 ! As a proportion to carbon storage times mean CN
     integer, parameter :: lf_store_prop = 5     ! leaf proportional storage
     integer, parameter :: store_prop = lf_store_prop

     
     select case(element_id)
     case(carbon12_element)
        write(fates_log(),*) 'Cannot call StorageNutrientTarget() for carbon'
        write(fates_log(),*) 'exiting'
        call endrun(msg=errMsg(sourcefile, __LINE__))
        
     case(nitrogen_element)
        
        if (store_prop == lfs_store_prop) then

           store_target  = prt_params%nitr_store_ratio(pft) * (leaf_target + sapw_target)

        elseif (store_prop == lf_store_prop) then
           
           store_target  = prt_params%nitr_store_ratio(pft) * leaf_target
   
        elseif(store_prop==lfss_store_prop) then
           
           store_target  = prt_params%nitr_store_ratio(pft) * (leaf_target + fnrt_target + sapw_target + struct_target)
           
        elseif(store_prop==fnrt_store_prop) then

           store_target  = prt_params%nitr_store_ratio(pft) * fnrt_target

        elseif(store_prop==cstore_store_prop) then

           !call bsap_allom(dbh,ipft,canopy_trim,elongf_stem,sapw_area,target_sapw_c)
           !call bagw_allom(dbh,ipft,elongf_stem,agw_c_target)
           !call bbgw_allom(dbh,ipft,elongf_stem,bgw_c_target)
           !call bdead_allom(agw_c_target,bgw_c_target,target_sapw_c,ipft,target_struct_c)
           !call bleaf(dbh,ipft,canopy_trim, elongf_leaf, target_leaf_c)
           !call bfineroot(dbh,ipft,canopy_trim, l2fr, elongf_fnrt, target_fnrt_c)
           !call bstore_allom(dbh,ipft,canopy_trim, target_store_c)

           ! Strategy, store as much nutrient as needed to match carbon's growth potential
           ! ie, nutrient storage is proportional to carbon storage times plant NC ratio

           ! N_so = a * C_so * NC_p
           ! NC_p = ( (N_so + N_lf + N_fr + N_sa + N_de)/C_tot )
           ! N_so = a * C_so * ( N_so/C_tot)  + a * C_so * (N_lf + N_fr + N_sa + N_de)/C_tot )
           ! N_so = (a * C_so * (N_lf + N_fr + N_sa + N_de)/C_tot ) / ( 1 - a * C_so/C_tot)

           !store_target = (target_store_c * prt_params%nitr_store_ratio(pft) * &
           !     (leaf_target + fnrt_target + sapw_target + struct_target)/total_c_target) / &
           !     ( 1._r8 - target_store_c * prt_params%nitr_store_ratio(pft) / total_c_target )
           write(fates_log(),*)'cstore_store_prop method of calculating target nutrient stores not available'
           call endrun(msg=errMsg(sourcefile, __LINE__))
           
        end if
        
             
     case(phosphorus_element)

        if (store_prop == lfs_store_prop) then
           
           store_target  = prt_params%phos_store_ratio(pft) * (leaf_target + fnrt_target + sapw_target)

        elseif (store_prop == lf_store_prop) then
           
           store_target  = prt_params%phos_store_ratio(pft) * leaf_target
    
        elseif(store_prop==lfss_store_prop) then
           
           store_target  = prt_params%phos_store_ratio(pft) * (leaf_target + fnrt_target + sapw_target + struct_target)
    
        elseif(store_prop==fnrt_store_prop) then
           
           store_target  = prt_params%phos_store_ratio(pft) * fnrt_target

         end if
     end select
     
     
   end function StorageNutrientTarget

  

end module PRTGenericMod
