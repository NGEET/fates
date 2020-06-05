module FatesRestartInterfaceMod


  use FatesConstantsMod,       only : r8 => fates_r8
  use FatesConstantsMod,       only : fates_avg_flag_length
  use FatesConstantsMod,       only : fates_short_string_length
  use FatesConstantsMod,       only : fates_long_string_length
  use FatesConstantsMod,       only : itrue
  use FatesConstantsMod,       only : ifalse
  use FatesConstantsMod,       only : fates_unset_r8
  use FatesConstantsMod,       only : primaryforest
  use FatesGlobals,            only : fates_log
  use FatesGlobals,            only : endrun => fates_endrun
  use FatesIODimensionsMod,    only : fates_io_dimension_type
  use FatesIOVariableKindMod,  only : fates_io_variable_kind_type
  use FatesRestartVariableMod, only : fates_restart_variable_type
  use FatesInterfaceTypesMod,       only : nlevcoage
  use FatesInterfaceTypesMod,       only : bc_in_type 
  use FatesInterfaceTypesMod,       only : bc_out_type
  use FatesInterfaceTypesMod,       only : hlm_use_planthydro
  use FatesInterfaceTypesMod,       only : fates_maxElementsPerSite
  use EDCohortDynamicsMod,     only : UpdateCohortBioPhysRates
  use FatesHydraulicsMemMod,   only : nshell
  use FatesHydraulicsMemMod,   only : n_hypool_ag
  use FatesHydraulicsMemMod,   only : n_hypool_troot
  use FatesHydraulicsMemMod,   only : nlevsoi_hyd_max
  use FatesPlantHydraulicsMod, only : UpdatePlantPsiFTCFromTheta
  use PRTGenericMod,           only : prt_global
  use EDCohortDynamicsMod,     only : nan_cohort
  use EDCohortDynamicsMod,     only : zero_cohort
  use EDCohortDynamicsMod,     only : InitPRTObject
  use EDCohortDynamicsMod,     only : InitPRTBoundaryConditions
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesInterfaceTypesMod,       only : nlevsclass
  use FatesLitterMod,          only : litter_type
  use FatesLitterMod,          only : ncwd
  use FatesLitterMod,          only : ndcmpy
  use PRTGenericMod,           only : prt_global
  use EDTypesMod,              only : num_elements


  ! CIME GLOBALS
  use shr_log_mod       , only : errMsg => shr_log_errMsg


  implicit none
  private ! Modules are private by default

  ! ------------------------------------------------------------
  ! A note on variable naming conventions.
  ! Many variables in this restart IO portion of the code will
  ! follow the conventions:
  !
  ! <use_case>_<description>_<dimension>
  !
  ! For instance we use an index for restart variable "ir_"
  ! to point the object that contains the number of patches per
  ! site "npatch" and this value is relevant to all sites "si"
  ! thus:   ir_npatch_si
  !
  ! We also use associations to the data arrays of restart IO
  ! variables "rio", for example the leaf litter "leaf_litter"
  ! is retrieved for every patch and every functional type "paft"
  ! thus: rio_leaf_litter_paft
  !
  ! si: site dimension
  ! pa: patch dimension
  ! co: cohort dimension
  ! ft: functional type dimension
  ! cl: canopy layer dimension (upper, lower, etc)
  ! ls: layer sublayer dimension (fine discretization of upper,lower)
  ! wm: the number of memory slots for water (currently 10)
  ! -------------------------------------------------------------
  
  
  ! Indices to the restart variable object

  integer :: ir_npatch_si 
  integer :: ir_cd_status_si
  integer :: ir_dd_status_si
  integer :: ir_nchill_days_si
  integer :: ir_ncold_days_si
  integer :: ir_leafondate_si
  integer :: ir_leafoffdate_si
  integer :: ir_dleafondate_si
  integer :: ir_dleafoffdate_si
  integer :: ir_acc_ni_si
  integer :: ir_gdd_si
  integer :: ir_trunk_product_si
  integer :: ir_ncohort_pa
  integer :: ir_canopy_layer_co
  integer :: ir_canopy_layer_yesterday_co
  integer :: ir_canopy_trim_co
  integer :: ir_size_class_lasttimestep_co
  integer :: ir_dbh_co
  integer :: ir_coage_co
  integer :: ir_g_sb_laweight_co
  integer :: ir_height_co
  integer :: ir_laimemory_co
  integer :: ir_sapwmemory_co
  integer :: ir_structmemory_co
  integer :: ir_nplant_co
  integer :: ir_gpp_acc_co
  integer :: ir_npp_acc_co
  integer :: ir_resp_acc_co
  integer :: ir_gpp_acc_hold_co
  integer :: ir_npp_acc_hold_co
  integer :: ir_resp_acc_hold_co
  integer :: ir_bmort_co
  integer :: ir_hmort_co
  integer :: ir_cmort_co
  integer :: ir_frmort_co
  integer :: ir_smort_co
  integer :: ir_asmort_co

  !Logging
  integer :: ir_lmort_direct_co
  integer :: ir_lmort_collateral_co
  integer :: ir_lmort_infra_co

  ! Radiation
  integer :: ir_solar_zenith_flag_pa
  integer :: ir_solar_zenith_angle_pa
  integer :: ir_gnd_alb_dif_pasb
  integer :: ir_gnd_alb_dir_pasb


  integer :: ir_ddbhdt_co
  integer :: ir_resp_tstep_co
  integer :: ir_pft_co
  integer :: ir_status_co
  integer :: ir_isnew_co

  ! Litter
  integer :: ir_agcwd_litt
  integer :: ir_bgcwd_litt
  integer :: ir_leaf_litt
  integer :: ir_fnrt_litt
  integer :: ir_seed_litt
  integer :: ir_seedgerm_litt
  integer :: ir_seed_prod_co
  integer :: ir_livegrass_pa
  integer :: ir_age_pa
  integer :: ir_area_pa
  integer :: ir_agesinceanthrodist_pa
  integer :: ir_patchdistturbcat_pa


  ! Site level
  integer :: ir_watermem_siwm
  integer :: ir_vegtempmem_sitm
  integer :: ir_seed_bank_sift
  integer :: ir_spread_si
  integer :: ir_recrate_sift
  integer :: ir_fmortrate_cano_siscpf
  integer :: ir_fmortrate_usto_siscpf
  integer :: ir_imortrate_siscpf
  integer :: ir_fmortrate_crown_siscpf
  integer :: ir_fmortrate_cambi_siscpf
  integer :: ir_termnindiv_cano_siscpf
  integer :: ir_termnindiv_usto_siscpf
  integer :: ir_growflx_fusion_siscpf
  integer :: ir_demorate_sisc
  integer :: ir_promrate_sisc
  integer :: ir_termcflux_cano_si
  integer :: ir_termcflux_usto_si
  integer :: ir_democflux_si
  integer :: ir_promcflux_si
  integer :: ir_imortcflux_si
  integer :: ir_fmortcflux_cano_si
  integer :: ir_fmortcflux_usto_si
  integer :: ir_cwdagin_flxdg
  integer :: ir_cwdbgin_flxdg
  integer :: ir_leaflittin_flxdg
  integer :: ir_rootlittin_flxdg
  integer :: ir_oldstock_mbal
  integer :: ir_errfates_mbal
  integer :: ir_prt_base     ! Base index for all PRT variables


  ! Hydraulic indices
  integer :: ir_hydro_th_ag_covec
  integer :: ir_hydro_th_troot
  integer :: ir_hydro_th_aroot_covec 
  integer :: ir_hydro_liqvol_shell_si
  integer :: ir_hydro_err_growturn_aroot
  integer :: ir_hydro_err_growturn_ag_covec
  integer :: ir_hydro_err_growturn_troot
  integer :: ir_hydro_recruit_si
  integer :: ir_hydro_dead_si
  integer :: ir_hydro_growturn_err_si
  integer :: ir_hydro_pheno_err_si
  integer :: ir_hydro_hydro_err_si

  ! The number of variable dim/kind types we have defined (static)
  integer, parameter, public :: fates_restart_num_dimensions = 2   !(cohort,column)
  integer, parameter, public :: fates_restart_num_dim_kinds = 4    !(cohort-int,cohort-r8,site-int,site-r8)

  ! integer constants for storing logical data
  integer, parameter, public :: old_cohort = 0
  integer, parameter, public :: new_cohort = 1  

  real(r8), parameter, public :: flushinvalid = -9999.0
  real(r8), parameter, public :: flushzero = 0.0
  real(r8), parameter, public :: flushone  = 1.0
  
  ! Local debug flag
  logical, parameter, public :: debug=.false.

  character(len=*), parameter :: sourcefile = &
       __FILE__

  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type, public :: restart_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: cohort1_index(:) ! maps site index to the HIO cohort 1st position
  end type restart_map_type



  type, public :: fates_restart_interface_type

     type(fates_restart_variable_type),allocatable :: rvars(:)
     integer,private :: num_restart_vars_

     ! Instanteate one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(fates_io_variable_kind_type) :: dim_kinds(fates_restart_num_dim_kinds)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure is
     ! allocated by number of threads. This could be dynamically
     ! allocated, but is unlikely to change...?
     ! Note: history io also instanteates fates_io_dimension_type
     type(fates_io_dimension_type) :: dim_bounds(fates_restart_num_dimensions)
     
     type(restart_map_type), pointer :: restart_map(:)

     integer, private :: cohort_index_, column_index_

   contains
     
     ! public functions
     procedure :: Init
     procedure :: SetThreadBoundsEach
     procedure :: assemble_restart_output_types
     procedure :: initialize_restart_vars
     procedure :: num_restart_vars
     procedure :: column_index
     procedure :: cohort_index
     procedure :: set_restart_vectors
     procedure :: create_patchcohort_structure
     procedure :: get_restart_vectors
     procedure :: update_3dpatch_radiation
     
     ! private work functions
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: set_cohort_index
     procedure, private :: set_column_index
     procedure, private :: flush_rvars
     procedure, private :: define_restart_vars
     procedure, private :: set_restart_var
     procedure, private :: DefinePRTRestartVars
     procedure, private :: GetCohortRealVector
     procedure, private :: SetCohortRealVector
     procedure, private :: RegisterCohortVector

  end type fates_restart_interface_type

  


contains

  ! =====================================================================================
  
  subroutine Init(this, num_threads, fates_bounds)
    
    use FatesIODimensionsMod, only : fates_bounds_type, column, cohort

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_cohort_index(dim_count)
    call this%dim_bounds(dim_count)%Init(cohort, num_threads, &
         fates_bounds%cohort_begin, fates_bounds%cohort_end)

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    ! FIXME(bja, 2016-10) assert(dim_count == FatesIOdimensionsmod::num_dimension_types)

    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%restart_map(num_threads))
    
  end subroutine Init  

  ! ======================================================================

  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)
    
    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index
    
    index = this%cohort_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cohort_begin, thread_bounds%cohort_end)
    
    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)
    
  end subroutine SetThreadBoundsEach

  ! ===================================================================================

  subroutine assemble_restart_output_types(this)
    
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int

    implicit none
   
    class(fates_restart_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(cohort_r8, 1, this%cohort_index())
    call this%set_dim_indices(cohort_int, 1, this%cohort_index())

    call this%set_dim_indices(site_r8, 1, this%column_index())
    call this%set_dim_indices(site_int, 1, this%column_index())

  end subroutine assemble_restart_output_types

 ! ===================================================================================
  
  subroutine set_dim_indices(this, dk_name, idim, dim_index)

    use FatesIOVariableKindMod , only : iotype_index

    implicit none

    ! arguments
    class(fates_restart_interface_type), intent(inout) :: this
    character(len=*), intent(in)     :: dk_name
    integer, intent(in)              :: idim  ! dimension index
    integer, intent(in) :: dim_index


    ! local
    integer :: ityp

    ityp = iotype_index(trim(dk_name), fates_restart_num_dim_kinds, this%dim_kinds)

    ! First check to see if the dimension is allocated
    if (this%dim_kinds(ityp)%ndims < idim) then
       write(fates_log(), *) 'Trying to define dimension size to a dim-type structure'
       write(fates_log(), *) 'but the dimension index does not exist'
       write(fates_log(), *) 'type: ',dk_name,' ndims: ',this%dim_kinds(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if

    if (idim == 1) then
       this%dim_kinds(ityp)%dim1_index = dim_index
    else if (idim == 2) then
       this%dim_kinds(ityp)%dim2_index = dim_index
    end if

    ! With the map, we can set the dimension size
    this%dim_kinds(ityp)%dimsize(idim) = this%dim_bounds(dim_index)%upper_bound - &
         this%dim_bounds(dim_index)%lower_bound + 1

 end subroutine set_dim_indices


  ! =======================================================================

  subroutine set_cohort_index(this, index)
    implicit none
    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%cohort_index_ = index
  end subroutine set_cohort_index
  
  integer function cohort_index(this)
    implicit none
    class(fates_restart_interface_type), intent(in) :: this
    cohort_index = this%cohort_index_
  end function cohort_index
  
  ! =======================================================================

  subroutine set_column_index(this, index)
    implicit none
    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%column_index_ = index
  end subroutine set_column_index
  
  integer function column_index(this)
    implicit none
    class(fates_restart_interface_type), intent(in) :: this
    column_index = this%column_index_
 end function column_index
 
 ! =======================================================================

 subroutine init_dim_kinds_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! CO_R8   : 1D cohort scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    ! CO_INT  : 1D cohort scale integers
    ! SI_INT  : 1D site scale integers
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int
    
    implicit none
    
    ! Arguments
    class(fates_restart_interface_type), intent(inout) :: this

    integer :: index

    ! 1d cohort r8
    index = 1
    call this%dim_kinds(index)%Init(cohort_r8, 1)

    ! 1d Site r8
    index = index + 1
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! cohort int
    index = index + 1
    call this%dim_kinds(index)%Init(cohort_int, 1)

    ! site int
    index = index + 1
    call this%dim_kinds(index)%Init(site_int, 1)

    ! FIXME(bja, 2016-10) assert(index == fates_num_dim_kinds)
  end subroutine init_dim_kinds_maps


  ! ====================================================================================

  integer function num_restart_vars(this)
    
    implicit none

    class(fates_restart_interface_type), intent(in) :: this

    num_restart_vars = this%num_restart_vars_
    
  end function num_restart_vars
  
  ! ====================================================================================
  
  subroutine initialize_restart_vars(this)

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this

   ! Determine how many of the restart IO variables registered in FATES
   ! are going to be allocated
   call this%define_restart_vars(initialize_variables=.false.)

   ! Allocate the list of restart output variable objects
   allocate(this%rvars(this%num_restart_vars()))
   
   ! construct the object that defines all of the IO variables
   call this%define_restart_vars(initialize_variables=.true.)
   
 end subroutine initialize_restart_vars

  ! ======================================================================================

 subroutine flush_rvars(this,nc)
 
   class(fates_restart_interface_type)        :: this
   integer,intent(in)                         :: nc

   integer                                   :: ivar
   type(fates_restart_variable_type),pointer :: rvar
   integer                      :: lb1,ub1,lb2,ub2

   do ivar=1,ubound(this%rvars,1)
      associate( rvar => this%rvars(ivar) )
        call rvar%Flush(nc, this%dim_bounds, this%dim_kinds)
      end associate
   end do
   
 end subroutine flush_rvars

 

 ! ====================================================================================
 
 subroutine define_restart_vars(this, initialize_variables)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF RESTART OUTPUT VARIABLES
    !
    ! Please add any restart variables to this registry. This registry will handle
    ! all variables that can make use of 1D column dimensioned or 1D cohort dimensioned
    ! variables.  Note that restarts are only using 1D vectors in ALM and CLM.  If you
    ! have a multi-dimensional variable that is below the cohort scale, then pack
    ! that variable into a cohort-sized output array by giving it a vtype "cohort_r8"
    ! or "cohort_int".  
    !
    ! Unlike history variables, restarts flush to zero.
    ! ---------------------------------------------------------------------------------
   
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_int, cohort_r8
    implicit none
    
    class(fates_restart_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?
    integer :: ivar
    
    
    ivar=0

    ! -----------------------------------------------------------------------------------
    ! Site level variables
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_PatchesPerSite', vtype=site_int, &
         long_name='Total number of FATES patches per column', units='none', flushval = flushinvalid, &
          hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npatch_si )

    call this%set_restart_var(vname='fates_cold_dec_status', vtype=site_int, &
         long_name='status flag for cold deciduous plants', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cd_status_si )

    call this%set_restart_var(vname='fates_drought_dec_status', vtype=site_int, &
         long_name='status flag for drought deciduous plants', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dd_status_si )

    call this%set_restart_var(vname='fates_chilling_days', vtype=site_int, &
         long_name='chilling day counter', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_nchill_days_si )

    call this%set_restart_var(vname='fates_cold_days', vtype=site_int, &
         long_name='cold day counter', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_ncold_days_si )

    call this%set_restart_var(vname='fates_leafondate', vtype=site_int, &
         long_name='the day of year for leaf on', units='day of year', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leafondate_si )

    call this%set_restart_var(vname='fates_leafoffdate', vtype=site_int, &
         long_name='the day of year for leaf off', units='day of year', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leafoffdate_si )

    call this%set_restart_var(vname='fates_drought_leafondate', vtype=site_int, &
         long_name='the day of year for drought based leaf-on', units='day of year', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dleafondate_si )

    call this%set_restart_var(vname='fates_drought_leafoffdate', vtype=site_int, &
         long_name='the day of year for drought based leaf-off', units='day of year', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dleafoffdate_si )

    call this%set_restart_var(vname='fates_acc_nesterov_id', vtype=site_r8, &
         long_name='a nesterov index accumulator', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_acc_ni_si )
    
    call this%set_restart_var(vname='fates_gdd_site', vtype=site_r8, &
         long_name='growing degree days at each site', units='degC days', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gdd_si )

    call this%set_restart_var(vname='fates_trunk_product_site', vtype=site_r8, &
         long_name='Accumulate trunk product flux at site', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_trunk_product_si )


    ! -----------------------------------------------------------------------------------
    ! Variables stored within cohort vectors
    ! Note: Some of these are multi-dimensional variables in the patch/site dimension
    ! that are collapsed into the cohort vectors for storage and transfer
    ! -----------------------------------------------------------------------------------

    ! This variable may be confusing, because it is a patch level variables
    ! but it is using the cohort IO vector to hold data
    call this%set_restart_var(vname='fates_CohortsPerPatch', vtype=cohort_int, &
         long_name='the number of cohorts per patch', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_ncohort_pa )

    call this%set_restart_var(vname='fates_solar_zenith_flag_pa', vtype=cohort_int, &
         long_name='switch specifying if zenith is positive', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_solar_zenith_flag_pa )
    
    call this%set_restart_var(vname='fates_solar_zenith_angle_pa', vtype=cohort_r8, &
         long_name='the angle of the solar zenith for each patch', units='radians', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_solar_zenith_angle_pa )



    ! 1D cohort Variables
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_seed_prod', vtype=cohort_r8, &
         long_name='fates cohort - seed production', units='kgC/plant', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_seed_prod_co )


    call this%set_restart_var(vname='fates_canopy_layer', vtype=cohort_int, &
         long_name='ed cohort - canopy_layer', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_layer_co )

    call this%set_restart_var(vname='fates_canopy_layer_yesterday', vtype=cohort_r8, &
         long_name='ed cohort - canopy_layer_yesterday', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_layer_yesterday_co )

    call this%set_restart_var(vname='fates_canopy_trim', vtype=cohort_r8, &
         long_name='ed cohort - canopy_trim', units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_trim_co )

    call this%set_restart_var(vname='fates_size_class_lasttimestep', vtype=cohort_int, &
         long_name='ed cohort - size-class last timestep', units='index', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_size_class_lasttimestep_co )

    call this%set_restart_var(vname='fates_dbh', vtype=cohort_r8, &
         long_name='ed cohort - diameter at breast height', units='cm', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dbh_co )

    call this%set_restart_var(vname='fates_coage', vtype=cohort_r8, &
         long_name='ed cohort - age in days', units='days', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_coage_co ) 

    call this%set_restart_var(vname='fates_height', vtype=cohort_r8, &
         long_name='ed cohort - plant height', units='m', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_height_co )

    call this%set_restart_var(vname='fates_laimemory', vtype=cohort_r8, &
         long_name='ed cohort - target leaf biomass set from prev year', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_laimemory_co )

    call this%set_restart_var(vname='fates_sapwmemory', vtype=cohort_r8, &
         long_name='ed cohort - target sapwood biomass set from prev year', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_sapwmemory_co )
    
    call this%set_restart_var(vname='fates_structmemory', vtype=cohort_r8, &
         long_name='ed cohort - target structural biomass set from prev year', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_structmemory_co )
    
    call this%set_restart_var(vname='fates_nplant', vtype=cohort_r8, &
         long_name='ed cohort - number of plants in the cohort', &
         units='/patch', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_nplant_co )

    call this%set_restart_var(vname='fates_gpp_acc', vtype=cohort_r8, &
         long_name='ed cohort - accumulated gpp over dynamics step', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gpp_acc_co )

    call this%set_restart_var(vname='fates_npp_acc', vtype=cohort_r8, &
         long_name='ed cohort - accumulated npp over dynamics step', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_acc_co )

    call this%set_restart_var(vname='fates_resp_acc', vtype=cohort_r8, &
         long_name='ed cohort - accumulated respiration over dynamics step', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_resp_acc_co )

    call this%set_restart_var(vname='fates_gpp_acc_hold', vtype=cohort_r8, &
         long_name='ed cohort - current step gpp', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gpp_acc_hold_co )

    call this%set_restart_var(vname='fates_npp_acc_hold', vtype=cohort_r8, &
         long_name='ed cohort - current step npp', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_acc_hold_co )

    call this%set_restart_var(vname='fates_resp_acc_hold', vtype=cohort_r8, &
         long_name='ed cohort - current step resp', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_resp_acc_hold_co )

    call this%set_restart_var(vname='fates_bmort', vtype=cohort_r8, &
         long_name='ed cohort - background mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bmort_co )

    call this%set_restart_var(vname='fates_hmort', vtype=cohort_r8, &
         long_name='ed cohort - hydraulic mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hmort_co )

    call this%set_restart_var(vname='fates_cmort', vtype=cohort_r8, &
         long_name='ed cohort - carbon starvation mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cmort_co )

    call this%set_restart_var(vname='fates_frmort', vtype=cohort_r8, &
         long_name='ed cohort - freezing mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_frmort_co )

    call this%set_restart_var(vname='fates_smort', vtype=cohort_r8, &
         long_name='ed cohort - senescence mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_smort_co )

    call this%set_restart_var(vname='fates_asmort', vtype=cohort_r8, &
         long_name='ed cohort - age senescence mortality rate', &
         units = '/year', flushval = flushzero, & 
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_asmort_co )

    call this%set_restart_var(vname='fates_lmort_direct', vtype=cohort_r8, &
         long_name='ed cohort - directly logging mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_direct_co )

    call this%set_restart_var(vname='fates_lmort_collateral', vtype=cohort_r8, &
         long_name='ed cohort - collateral mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_collateral_co ) 
  
    call this%set_restart_var(vname='fates_lmort_in', vtype=cohort_r8, &
         long_name='ed cohort - mechanical mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_infra_co ) 

    call this%set_restart_var(vname='fates_ddbhdt', vtype=cohort_r8, &
         long_name='ed cohort - differential: ddbh/dt', &
         units='cm/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_ddbhdt_co )

    call this%set_restart_var(vname='fates_resp_tstep', vtype=cohort_r8, &
         long_name='ed cohort - autotrophic respiration over timestep', &
         units='kgC/indiv/timestep', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_resp_tstep_co )

    call this%set_restart_var(vname='fates_pft', vtype=cohort_int, &
         long_name='ed cohort - plant functional type', units='index', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_pft_co )

    call this%set_restart_var(vname='fates_status_coh', vtype=cohort_int, &
         long_name='ed cohort - plant phenology status', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_status_co )

    call this%set_restart_var(vname='fates_isnew', vtype=cohort_int, &
         long_name='ed cohort - binary flag specifying if a plant has experienced a full day cycle', &
         units='0/1', flushval = flushone, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_isnew_co )

    call this%set_restart_var(vname='fates_gsblaweight',vtype=cohort_r8, &
         long_name='ed cohort - leaf-area weighted total stomatal+blayer conductance', &
         units='[m/s]*[m2]', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_g_sb_laweight_co)

    ! Mixed dimension variables using the cohort vector
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_gnd_alb_dif', vtype=cohort_r8, &
         long_name='ground albedo of diffuse radiation vis and ir', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gnd_alb_dif_pasb )

    call this%set_restart_var(vname='fates_gnd_alb_dir', vtype=cohort_r8, &
         long_name='ground albedo of direct radiation vis and ir', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gnd_alb_dir_pasb )

    call this%set_restart_var(vname='fates_spread', vtype=site_r8, &
         long_name='dynamic ratio of dbh to canopy area, by patch x canopy-layer', &
         units='cm/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_spread_si )

    call this%set_restart_var(vname='fates_livegrass', vtype=cohort_r8, &
         long_name='total AGB from grass, by patch', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_livegrass_pa )

    call this%set_restart_var(vname='fates_age', vtype=cohort_r8, &
         long_name='age of the ED patch', units='yr', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_age_pa )

    call this%set_restart_var(vname='fates_age_since_anthro_dist', vtype=cohort_r8, &
         long_name='age of the ED patch since last anthropogenic disturbance', &
         units='yr', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, &
         index = ir_agesinceanthrodist_pa )

    call this%set_restart_var(vname='fates_patchdistturbcat', vtype=cohort_int, &
         long_name='Disturbance label of patch', units='yr', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_patchdistturbcat_pa )

    call this%set_restart_var(vname='fates_area', vtype=cohort_r8, &
         long_name='are of the ED patch', units='m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_area_pa )


    ! Site Level Diagnostics over multiple nutrients


    ! Patch Level Litter Pools are potentially multi-element

    call this%RegisterCohortVector(symbol_base='fates_ag_cwd', vtype=cohort_r8, &
            long_name_base='above ground CWD',  &
            units='kg/m2', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_agcwd_litt) 

    call this%RegisterCohortVector(symbol_base='fates_bg_cwd', vtype=cohort_r8, &
            long_name_base='below ground CWD',  &
            units='kg/m2', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bgcwd_litt) 

    call this%RegisterCohortVector(symbol_base='fates_leaf_fines', vtype=cohort_r8, &
            long_name_base='above ground leaf litter',  &
            units='kg/m2', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leaf_litt) 

    call this%RegisterCohortVector(symbol_base='fates_fnrt_fines', vtype=cohort_r8, &
            long_name_base='fine root litter',  &
            units='kg/m2', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fnrt_litt) 
    
    call this%RegisterCohortVector(symbol_base='fates_seed', vtype=cohort_r8, &
            long_name_base='seed bank (non-germinated)',  &
            units='kg/m2', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_seed_litt)

    call this%RegisterCohortVector(symbol_base='fates_seedgerm', vtype=cohort_r8, &
           long_name_base='seed bank (germinated)',  &
           units='kg/m2', veclength=num_elements, flushval = flushzero, &
           hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_seedgerm_litt)


    ! Site level flux diagnostics for each element

    call this%RegisterCohortVector(symbol_base='fates_cwdagin', vtype=cohort_r8, &
            long_name_base='Input flux of AG CWD', &
            units='kg/ha', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cwdagin_flxdg)

    call this%RegisterCohortVector(symbol_base='fates_cwdbgin', vtype=cohort_r8, &
            long_name_base='Input flux of BG CWD', &
            units='kg/ha', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cwdbgin_flxdg)

    call this%RegisterCohortVector(symbol_base='fates_leaflittin', vtype=cohort_r8, &
            long_name_base='Input flux of leaf litter', &
            units='kg/ha', veclength=num_elements, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leaflittin_flxdg)

    call this%RegisterCohortVector(symbol_base='fates_rootlittin', vtype=cohort_r8, &
           long_name_base='Input flux of root litter', &
           units='kg/ha', veclength=num_elements, flushval = flushzero, &
           hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_rootlittin_flxdg)

    ! Site level Mass Balance State Accounting

    call this%RegisterCohortVector(symbol_base='fates_oldstock', vtype=site_r8, &
         long_name_base='Previous total mass of all fates state variables', &
         units='kg/ha', veclength=num_elements, flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_oldstock_mbal)
    
    call this%RegisterCohortVector(symbol_base='fates_errfates', vtype=site_r8, &
         long_name_base='Previous total mass of error fates state variables', &
         units='kg/ha', veclength=num_elements, flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_errfates_mbal)
    
    

    ! Only register hydraulics restart variables if it is turned on!
    
    if(hlm_use_planthydro==itrue) then

       if ( fates_maxElementsPerSite < (nshell * nlevsoi_hyd_max) ) then
          write(fates_log(), *) ' Ftes plant hydraulics needs space to store site-level hydraulics info.'
          write(fates_log(), *) ' It uses array spaces typically reserved for cohorts to hold this.'
          write(fates_log(), *) ' However, that space defined by fates_maxElementsPerSite must be larger'
          write(fates_log(), *) ' than the product of maximum soil layers x rhizosphere shells'
          write(fates_log(), *) ' See FatesInterfaceMod.F90 for how this array is set'
          write(fates_log(), *) ' fates_maxElementsPerSite = ',fates_maxElementsPerSite
          write(fates_log(), *) ' nshell = ',nshell
          write(fates_log(), *) ' nlevsoi_hyd_max = ',nlevsoi_hyd_max
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       call this%RegisterCohortVector(symbol_base='fates_hydro_th_ag', vtype=cohort_r8, &
            long_name_base='water in aboveground compartments',  &
            units='kg/plant', veclength=n_hypool_ag, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_th_ag_covec) 
       
       call this%RegisterCohortVector(symbol_base='fates_hydro_th_troot', vtype=cohort_r8, &
            long_name_base='water in transporting roots', &
            units='kg/plant', veclength=n_hypool_troot, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_th_troot) 
       
       call this%RegisterCohortVector(symbol_base='fates_hydro_th_aroot', vtype=cohort_r8, &
            long_name_base='water in absorbing roots',  &
            units='kg/plant', veclength=nlevsoi_hyd_max, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_th_aroot_covec) 

       call this%RegisterCohortVector(symbol_base='fates_hydro_err_aroot', vtype=cohort_r8, &
            long_name_base='error in plant-hydro balance in absorbing roots',  &
            units='kg/plant', veclength=nlevsoi_hyd_max, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_err_growturn_aroot) 

       call this%RegisterCohortVector(symbol_base='fates_hydro_err_ag', vtype=cohort_r8, &
            long_name_base='error in plant-hydro balance above ground',  &
            units='kg/plant', veclength=n_hypool_ag, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_err_growturn_ag_covec) 

       call this%RegisterCohortVector(symbol_base='fates_hydro_err_troot', vtype=cohort_r8, &
            long_name_base='error in plant-hydro balance above ground',  &
            units='kg/plant', veclength=n_hypool_troot, flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_err_growturn_troot) 

       ! Site-level volumentric liquid water content (shell x layer)
       call this%set_restart_var(vname='fates_hydro_liqvol_shell', vtype=cohort_r8, &
            long_name='Volumetric water content of rhizosphere compartments (layerxshell)', &
            units='m3/m3', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_liqvol_shell_si )

       ! Site-level water bound in new recruits
       call this%set_restart_var(vname='fates_hydro_recruit_h2o', vtype=site_r8, &
            long_name='Site level water mass used for new recruits', &
            units='kg', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_recruit_si )
       
       ! Site-level water bound in dead plants
       call this%set_restart_var(vname='fates_hydro_dead_h2o', vtype=site_r8, &
            long_name='Site level water bound in dead plants', &
            units='kg', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_dead_si )
       
       ! Site-level water balance error due to growth/turnover
       call this%set_restart_var(vname='fates_hydro_growturn_err', vtype=site_r8, &
            long_name='Site level error for hydraulics due to growth/turnover', &
            units='kg', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_growturn_err_si )

       ! Site-level water balance error due to phenology?
       call this%set_restart_var(vname='fates_hydro_pheno_err', vtype=site_r8, &
            long_name='Site level error for hydraulics due to phenology', &
            units='kg', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_pheno_err_si )

       ! Site-level water balance error in vegetation
       call this%set_restart_var(vname='fates_hydro_hydro_err', vtype=site_r8, &
            long_name='Site level error for hydrodynamics', &
            units='kg', flushval = flushzero, &
            hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hydro_hydro_err_si )
       
    end if


    !
    ! site x time level vars
    !

    call this%set_restart_var(vname='fates_water_memory', vtype=cohort_r8, &
         long_name='last 10 days of volumetric soil water, by site x day-index', &
         units='m3/m3', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_watermem_siwm )

    call this%set_restart_var(vname='fates_vegtemp_memory', vtype=cohort_r8, &
         long_name='last 10 days of 24-hour vegetation temperature, by site x day-index', &
         units='m3/m3', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_vegtempmem_sitm )
    
    call this%set_restart_var(vname='fates_recrate', vtype=cohort_r8, &
         long_name='fates diagnostics on recruitment', &
         units='indiv/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_recrate_sift)

    call this%set_restart_var(vname='fates_fmortrate_canopy', vtype=cohort_r8, &
         long_name='fates diagnostics on fire mortality canopy', &
         units='indiv/ha/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortrate_cano_siscpf)

    call this%set_restart_var(vname='fates_fmortrate_ustory', vtype=cohort_r8, &
         long_name='fates diagnostics on fire mortality ustory', &
         units='indiv/ha/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortrate_usto_siscpf)

    call this%set_restart_var(vname='fates_imortrate', vtype=cohort_r8, &
         long_name='fates diagnostics on impact mortality', &
         units='indiv/ha/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_imortrate_siscpf)

    call this%set_restart_var(vname='fates_fmortrate_crown', vtype=cohort_r8, &
         long_name='fates diagnostics on crown fire mortality', &
         units='indiv/ha/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortrate_crown_siscpf)

    call this%set_restart_var(vname='fates_fmortrate_cambi', vtype=cohort_r8, &
         long_name='fates diagnostics on fire cambial mortality', &
         units='indiv/ha/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortrate_cambi_siscpf)

    call this%set_restart_var(vname='fates_termn_canopy', vtype=cohort_r8, &
         long_name='fates diagnostics on termin mortality canopy', &
         units='indiv/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_termnindiv_cano_siscpf)

    call this%set_restart_var(vname='fates_termn_ustory', vtype=cohort_r8, &
         long_name='fates diagnostics on term mortality ustory', &
         units='indiv/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_termnindiv_usto_siscpf)

    call this%set_restart_var(vname='fates_growflx_fusion', vtype=cohort_r8, &
         long_name='fates diag: rate of indivs moving via fusion', &
         units='indiv/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_growflx_fusion_siscpf)
    
    call this%set_restart_var(vname='fates_demorate', vtype=cohort_r8, &
         long_name='fates diagnoatic rate of indivs demoted', &
         units='indiv/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_demorate_sisc)

    call this%set_restart_var(vname='fates_promrate', vtype=cohort_r8, &
         long_name='fates diagnostic rate of indivs promoted', &
         units='indiv/ha/da', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_promrate_sisc)

    call this%set_restart_var(vname='fates_imortcflux', vtype=site_r8, &
         long_name='biomass of indivs killed due to impact mort', &
         units='kgC/ha/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_imortcflux_si)
 
   call this%set_restart_var(vname='fates_fmortcflux_canopy', vtype=site_r8, &
         long_name='fates diagnostic biomass of canopy fire', &
         units='gC/m2/sec', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortcflux_cano_si)

    call this%set_restart_var(vname='fates_fmortcflux_ustory', vtype=site_r8, &
         long_name='fates diagnostic biomass of understory fire', &
         units='gC/m2/sec', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmortcflux_usto_si)

    call this%set_restart_var(vname='fates_termcflux_canopy', vtype=site_r8, &
         long_name='fates diagnostic term carbon flux canopy', &
         units='', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index =   ir_termcflux_cano_si )

   call this%set_restart_var(vname='fates_termcflux_ustory', vtype=site_r8, &
         long_name='fates diagnostic term carbon flux understory', &
         units='', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index =   ir_termcflux_usto_si )

   call this%set_restart_var(vname='fates_democflux', vtype=site_r8, &
         long_name='fates diagnostic demotion carbon flux', &
         units='', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index =   ir_democflux_si )

   call this%set_restart_var(vname='fates_promcflux', vtype=site_r8, &
         long_name='fates diagnostic promotion carbon flux ', &
         units='', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index =   ir_promcflux_si )




    ! Register all of the PRT states and fluxes

    ir_prt_base = ivar
    call this%DefinePRTRestartVars(initialize_variables,ivar)
       
 
    
    ! Must be last thing before return
    this%num_restart_vars_ = ivar
    
 end subroutine define_restart_vars
  
 ! =====================================================================================
 
  subroutine DefinePRTRestartVars(this,initialize_variables,ivar)

    ! ----------------------------------------------------------------------------------
    ! PARTEH variables are objects.  These objects 
    ! each are registered to have things like names units and symbols
    ! as part of that object.  Thus, when defining, reading and writing restarts,
    ! instead of manually typing out each variable we want, we just loop through
    ! our list of ojbects.
    !
    ! We do have to loop through the different parts of the objects indepenently.
    ! For instance we can't have one loop that covers the states "val", and
    ! the net allocation and reactive transport flux "net_alloc", so we have to loop
    ! these each separately. As other fluxes are added in the future, they need
    ! their own definition.
    !
    ! Some of the code below is about parsing the strings of these objects
    ! and automatically building the names of the PARTEH output variables
    ! as we go.
    !
    ! Note that parteh variables may or may not be scalars. Each variable's
    ! position gets its own variable in the restart file.  So the variable
    ! name will also parse the string for that position.
    ! -----------------------------------------------------------------------------------


     use FatesIOVariableKindMod, only : cohort_r8

     class(fates_restart_interface_type) :: this
     logical, intent(in)                 :: initialize_variables
     integer,intent(inout)               :: ivar      ! global variable counter
      
     integer                             :: dummy_out ! dummy index for variable
                                                      ! position in global file
     integer                             :: i_var     ! loop counter for prt variables
     integer                             :: i_pos     ! loop counter for discrete position

     character(len=32)  :: symbol_base    ! Symbol name without position or flux type
     character(len=128) :: name_base      ! name without position or flux type
     character(len=4)   :: pos_symbol
     character(len=128) :: symbol
     character(len=256) :: long_name

     do i_var = 1, prt_global%num_vars

        ! The base symbol name
        symbol_base = prt_global%state_descriptor(i_var)%symbol
        
        ! The long name of the variable
        name_base = prt_global%state_descriptor(i_var)%longname

        do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
           
           ! String describing the physical position of the variable
           write(pos_symbol, '(I3.3)') i_pos

           ! Register the instantaneous state variable "val"
           ! ----------------------------------------------------------------------------

           ! The symbol that is written to file
           symbol    = trim(symbol_base)//'_val_'//trim(pos_symbol)

           ! The expanded long name of the variable
           long_name = trim(name_base)//', state var, position:'//trim(pos_symbol)

           call this%set_restart_var(vname=trim(symbol), &
                  vtype=cohort_r8, &
                  long_name=trim(long_name), &
                  units='kg', flushval = flushzero, &
                  hlms='CLM:ALM', initialize=initialize_variables, &
                  ivar=ivar, index = dummy_out ) 

           ! Register the turnover flux variables
           ! ----------------------------------------------------------------------------

           ! The symbol that is written to file
           symbol = trim(symbol_base)//'_turn_'//trim(pos_symbol)

           ! The expanded long name of the variable
           long_name     = trim(name_base)//', turnover, position:'//trim(pos_symbol)
           
           call this%set_restart_var(vname=trim(symbol), &
                 vtype=cohort_r8, &
                 long_name=trim(long_name), &
                 units='kg', flushval = flushzero, &
                 hlms='CLM:ALM', initialize=initialize_variables, &
                 ivar=ivar, index = dummy_out ) 
            


           ! Register the net allocation flux variable
           ! ----------------------------------------------------------------------------
           
           ! The symbol that is written to file
           symbol = trim(symbol_base)//'_net_'//trim(pos_symbol)

           ! The expanded long name of the variable
           long_name     = trim(name_base)//', net allocation/transp, position:'//trim(pos_symbol)

           call this%set_restart_var(vname=trim(symbol), &
                  vtype=cohort_r8, &
                  long_name=trim(long_name), &
                  units='kg', flushval = flushzero, &
                  hlms='CLM:ALM', initialize=initialize_variables, &
                  ivar=ivar, index = dummy_out ) 
           


           ! Register the burn flux variable
           ! ----------------------------------------------------------------------------
           ! The symbol that is written to file
           symbol    = trim(symbol_base)//'_burned_'//trim(pos_symbol)

           ! The expanded long name of the variable
           long_name = trim(name_base)//', burned mass:'//trim(pos_symbol)

           call this%set_restart_var(vname=symbol, &
                 vtype=cohort_r8, &
                 long_name=trim(long_name), &
                 units='kg', flushval = flushzero, &
                 hlms='CLM:ALM', initialize=initialize_variables, &
                 ivar=ivar, index = dummy_out ) 

        end do
     end do
      
     return
  end subroutine DefinePRTRestartVars

  ! =====================================================================================

  subroutine RegisterCohortVector(this,symbol_base, vtype, long_name_base, &
                                  units, veclength, flushval, hlms,   &
                                  initialize, ivar, index) 

       
    ! The basic idea here is that instead of saving cohorts with vector data
    ! as long arrays in the restart file, we give each index of the vector
    ! its own variable.  This helps reduce the size of the restart files
    ! considerably.
    
    
    use FatesIOVariableKindMod, only : cohort_r8
    
    class(fates_restart_interface_type) :: this
    character(*),intent(in) :: symbol_base    ! Symbol name without position
    character(*),intent(in) :: vtype          ! String defining variable type 
    character(*),intent(in) :: long_name_base ! name without position
    character(*),intent(in) :: units          ! units for this variable
    integer,intent(in)      :: veclength      ! length of the vector
    real(r8),intent(in)     :: flushval       ! Value to flush to
    character(*),intent(in) :: hlms           ! The HLMs this works in
    logical, intent(in)     :: initialize     ! Is this registering or counting?
    integer,intent(inout)   :: ivar           ! global variable counter
    integer,intent(out)     :: index          ! The variable index for this variable
    
    ! Local Variables
    character(len=4)        :: pos_symbol     ! vectors need text strings for each position
    character(len=128)      :: symbol         ! symbol  name written to file
    character(len=256)      :: long_name      ! long name written to file
    integer                 :: i_pos          ! loop counter for discrete position
    integer                 :: dummy_index
    

    ! We give each vector its own index that points to the first position
    
    index = ivar + 1
    
    do i_pos = 1, veclength
       
       ! String describing the physical position of the variable
       write(pos_symbol, '(I3.3)') i_pos
       
       ! The symbol that is written to file
       symbol    = trim(symbol_base)//'_vec_'//trim(pos_symbol)
       
       ! The expanded long name of the variable
       long_name = trim(long_name_base)//', position:'//trim(pos_symbol)
       
       call this%set_restart_var(vname=trim(symbol), &
            vtype=vtype, &
            long_name=trim(long_name), &
            units=units, flushval = flushval, &
            hlms='CLM:ALM', initialize=initialize, &
            ivar=ivar, index = dummy_index ) 
       
    end do
    
  end subroutine RegisterCohortVector

  ! =====================================================================================
  
  subroutine GetCohortRealVector(this, state_vector, len_state_vector, &
                                 variable_index_base, co_global_index)
    
    ! This subroutine walks through global cohort vector indices
    ! and pulls from the different associated restart variables
    
    class(fates_restart_interface_type) , intent(inout) :: this
    integer,intent(in)     :: len_state_vector
    real(r8),intent(inout) :: state_vector(len_state_vector)
    integer,intent(in)     :: variable_index_base
    integer,intent(in)     :: co_global_index
    
    integer :: i_pos              ! vector position loop index
    integer :: ir_pos_var         ! global variable index
    
    ir_pos_var = variable_index_base
    do i_pos = 1, len_state_vector
       state_vector(i_pos) = this%rvars(ir_pos_var)%r81d(co_global_index)
       ir_pos_var = ir_pos_var + 1
    end do
    return
 end subroutine GetCohortRealVector
  
  ! =====================================================================================   
  
  subroutine SetCohortRealVector(this, state_vector, len_state_vector, &
                                  variable_index_base, co_global_index)

    ! This subroutine walks through global cohort vector indices
    ! and pushes into the restart arrays the different associated restart variables
    
    class(fates_restart_interface_type) , intent(inout) :: this
    integer,intent(in)   :: len_state_vector
    real(r8),intent(in)  :: state_vector(len_state_vector)
    integer,intent(in)   :: variable_index_base
    integer,intent(in)   :: co_global_index
    
    integer :: i_pos              ! vector position loop index
    integer :: ir_pos_var         ! global variable index
    
    ir_pos_var = variable_index_base
    do i_pos = 1, len_state_vector
       this%rvars(ir_pos_var)%r81d(co_global_index) = state_vector(i_pos)
       ir_pos_var = ir_pos_var + 1
    end do
    return
  end subroutine SetCohortRealVector
  

  ! =====================================================================================

  subroutine set_restart_var(this,vname,vtype,long_name,units,flushval, &
        hlms,initialize,ivar,index)

    use FatesUtilsMod, only : check_hlm_list
    use FatesInterfaceTypesMod, only : hlm_name

    ! arguments
    class(fates_restart_interface_type) :: this
    character(len=*),intent(in)  :: vname
    character(len=*),intent(in)  :: vtype
    character(len=*),intent(in)  :: units 
    real(r8), intent(in)         :: flushval
    character(len=*),intent(in)  :: long_name
    character(len=*),intent(in)  :: hlms
    logical, intent(in)          :: initialize
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    
    type(fates_restart_variable_type),pointer :: rvar
    integer :: ub1,lb1,ub2,lb2    ! Bounds for allocating the var
    integer :: ityp
    
    logical :: use_var
    
    use_var = check_hlm_list(trim(hlms), trim(hlm_name))


    if( use_var ) then
       
       ivar  = ivar+1
       index = ivar    
       
       if( initialize )then
          
          call this%rvars(ivar)%Init(vname, units, long_name, vtype, flushval, &
               fates_restart_num_dim_kinds, this%dim_kinds, this%dim_bounds)

       end if
    else
       
       index = 0
    end if
    
    return
 end subroutine set_restart_var

 ! =====================================================================================

 subroutine set_restart_vectors(this,nc,nsites,sites)

   use FatesInterfaceTypesMod, only : fates_maxElementsPerPatch
   use FatesInterfaceTypesMod, only : numpft
   use EDTypesMod, only : ed_site_type
   use EDTypesMod, only : ed_cohort_type
   use EDTypesMod, only : ed_patch_type
   use EDTypesMod, only : maxSWb
   use EDTypesMod, only : numWaterMem
   use EDTypesMod, only : num_vegtemp_mem

    ! Arguments
    class(fates_restart_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)

    ! Locals
    integer  :: s                         ! The local site index
    type(litter_type), pointer :: litt    ! pointer to patch's litter object

    ! ----------------------------------------------------------------------------------
    ! The following group of integers indicate the positional index (idx)
    ! of variables at different scales inside the I/O arrays (io)
    ! Keep in mind that many of these variables have a composite dimension
    ! at the patch scale.  To hold this memory, we borrow the cohort
    ! vector.  Thus the head of each array points to the first cohort
    ! of each patch. "io_idx_co_1st"
    ! ----------------------------------------------------------------------------------
    integer  :: io_idx_si      ! site index
    integer  :: io_idx_co_1st  ! 1st cohort of each patch
    integer  :: io_idx_co      ! cohort index
    integer  :: io_idx_pa_pft  ! each pft within each patch (pa_pft)
    integer  :: io_idx_pa_cwd  ! each cwd class within each patch (pa_cwd)
    integer  :: io_idx_pa_cwsl ! each cwd x soil layer
    integer  :: io_idx_pa_dcsl ! each decomposability x soil layer
    integer  :: io_idx_pa_dc   ! each decomposability index
    integer  :: io_idx_pa_ib   ! each SW band (vis/ir) per patch (pa_ib)
    integer  :: io_idx_si_wmem ! each water memory class within each site
    integer  :: io_idx_si_lyr_shell ! site - layer x shell index
    integer  :: io_idx_si_scpf ! each size-class x pft index within site
    integer  :: io_idx_si_sc   ! each size-class index within site
    integer  :: io_idx_si_capf ! each cohort age-class x pft index within site
    integer  :: io_idx_si_cacls ! each cohort age class index within site
    integer  :: io_idx_si_cwd  ! each site-cwd index
    integer  :: io_idx_si_pft  ! each site-pft index
    integer  :: io_idx_si_vtmem ! indices for veg-temp memory at site


    ! Some counters (for checking mostly)
    integer  :: totalcohorts   ! total cohort count on this thread (diagnostic)
    integer  :: patchespersite   ! number of patches per site
    integer  :: cohortsperpatch  ! number of cohorts per patch 

    integer  :: ft               ! functional type index
    integer  :: el               ! element loop index
    integer  :: ilyr             ! soil layer index
    integer  :: nlevsoil         ! total soil layers in patch of interest
    integer  :: k,j,i            ! indices to the radiation matrix
    integer  :: ir_prt_var       ! loop counter for var x position
    integer  :: i_var            ! loop counter for PRT variables
    integer  :: i_pos            ! loop counter for discrete PRT positions
    integer  :: i_scls           ! loop counter for size-class
    integer  :: i_cacls          ! loop counter for cohort age class
    integer  :: i_cwd            ! loop counter for cwd
    integer  :: i_pft            ! loop counter for pft

    type(fates_restart_variable_type) :: rvar
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort


    associate( rio_npatch_si           => this%rvars(ir_npatch_si)%int1d, &
           rio_cd_status_si            => this%rvars(ir_cd_status_si)%int1d, &
           rio_dd_status_si            => this%rvars(ir_dd_status_si)%int1d, &
           rio_nchill_days_si          => this%rvars(ir_nchill_days_si)%int1d, &
           rio_ncold_days_si           => this%rvars(ir_ncold_days_si)%int1d, &
           rio_leafondate_si           => this%rvars(ir_leafondate_si)%int1d, &
           rio_leafoffdate_si          => this%rvars(ir_leafoffdate_si)%int1d, &
           rio_dleafondate_si          => this%rvars(ir_dleafondate_si)%int1d, &
           rio_dleafoffdate_si         => this%rvars(ir_dleafoffdate_si)%int1d, &
           rio_acc_ni_si               => this%rvars(ir_acc_ni_si)%r81d, &
           rio_gdd_si                  => this%rvars(ir_gdd_si)%r81d, &
           rio_trunk_product_si        => this%rvars(ir_trunk_product_si)%r81d, &
           rio_ncohort_pa              => this%rvars(ir_ncohort_pa)%int1d, &
           rio_solar_zenith_flag_pa    => this%rvars(ir_solar_zenith_flag_pa)%int1d, &
           rio_solar_zenith_angle_pa   => this%rvars(ir_solar_zenith_angle_pa)%r81d, &
           rio_canopy_layer_co         => this%rvars(ir_canopy_layer_co)%int1d, &
           rio_canopy_layer_yesterday_co    => this%rvars(ir_canopy_layer_yesterday_co)%r81d, &
           rio_canopy_trim_co          => this%rvars(ir_canopy_trim_co)%r81d, &
           rio_seed_prod_co            => this%rvars(ir_seed_prod_co)%r81d, &
           rio_size_class_lasttimestep => this%rvars(ir_size_class_lasttimestep_co)%int1d, &
           rio_dbh_co                  => this%rvars(ir_dbh_co)%r81d, &
           rio_coage_co                => this%rvars(ir_coage_co)%r81d, &
           rio_g_sb_laweight_co        => this%rvars(ir_g_sb_laweight_co)%r81d, &
           rio_height_co               => this%rvars(ir_height_co)%r81d, &
           rio_laimemory_co            => this%rvars(ir_laimemory_co)%r81d, &
           rio_sapwmemory_co           => this%rvars(ir_sapwmemory_co)%r81d, &
           rio_structmemory_co         => this%rvars(ir_structmemory_co)%r81d, &
           rio_nplant_co               => this%rvars(ir_nplant_co)%r81d, &
           rio_gpp_acc_co              => this%rvars(ir_gpp_acc_co)%r81d, &
           rio_npp_acc_co              => this%rvars(ir_npp_acc_co)%r81d, &
           rio_resp_acc_co             => this%rvars(ir_resp_acc_co)%r81d, &
           rio_gpp_acc_hold_co         => this%rvars(ir_gpp_acc_hold_co)%r81d, &
           rio_resp_acc_hold_co        => this%rvars(ir_resp_acc_hold_co)%r81d, &
           rio_npp_acc_hold_co         => this%rvars(ir_npp_acc_hold_co)%r81d, &
           rio_bmort_co                => this%rvars(ir_bmort_co)%r81d, &
           rio_hmort_co                => this%rvars(ir_hmort_co)%r81d, &
           rio_cmort_co                => this%rvars(ir_cmort_co)%r81d, &
           rio_smort_co                => this%rvars(ir_smort_co)%r81d, &
           rio_asmort_co               => this%rvars(ir_asmort_co)%r81d, &
           rio_frmort_co               => this%rvars(ir_frmort_co)%r81d, &
           rio_lmort_direct_co         => this%rvars(ir_lmort_direct_co)%r81d, &
           rio_lmort_collateral_co     => this%rvars(ir_lmort_collateral_co)%r81d, &
           rio_lmort_infra_co          => this%rvars(ir_lmort_infra_co)%r81d, &
           rio_ddbhdt_co               => this%rvars(ir_ddbhdt_co)%r81d, &
           rio_resp_tstep_co           => this%rvars(ir_resp_tstep_co)%r81d, &
           rio_pft_co                  => this%rvars(ir_pft_co)%int1d, &
           rio_status_co               => this%rvars(ir_status_co)%int1d, &
           rio_isnew_co                => this%rvars(ir_isnew_co)%int1d, &
           rio_gnd_alb_dif_pasb        => this%rvars(ir_gnd_alb_dif_pasb)%r81d, &
           rio_gnd_alb_dir_pasb        => this%rvars(ir_gnd_alb_dir_pasb)%r81d, &
           rio_spread_si               => this%rvars(ir_spread_si)%r81d, &
           rio_livegrass_pa            => this%rvars(ir_livegrass_pa)%r81d, &
           rio_age_pa                  => this%rvars(ir_age_pa)%r81d, &
           rio_patchdistturbcat_pa     => this%rvars(ir_patchdistturbcat_pa)%int1d, &           
           rio_agesinceanthrodist_pa   => this%rvars(ir_agesinceanthrodist_pa)%r81d, &           
           rio_area_pa                 => this%rvars(ir_area_pa)%r81d, &
           rio_watermem_siwm           => this%rvars(ir_watermem_siwm)%r81d, &
           rio_vegtempmem_sitm         => this%rvars(ir_vegtempmem_sitm)%r81d, &
           rio_recrate_sift            => this%rvars(ir_recrate_sift)%r81d, &
           rio_fmortrate_cano_siscpf   => this%rvars(ir_fmortrate_cano_siscpf)%r81d, &
           rio_fmortrate_usto_siscpf   => this%rvars(ir_fmortrate_usto_siscpf)%r81d, &
           rio_imortrate_siscpf        => this%rvars(ir_imortrate_siscpf)%r81d, &
           rio_fmortrate_crown_siscpf  => this%rvars(ir_fmortrate_crown_siscpf)%r81d, &
           rio_fmortrate_cambi_siscpf  => this%rvars(ir_fmortrate_cambi_siscpf)%r81d, &
           rio_termnindiv_cano_siscpf  => this%rvars(ir_termnindiv_cano_siscpf)%r81d, &
           rio_termnindiv_usto_siscpf  => this%rvars(ir_termnindiv_usto_siscpf)%r81d, &
           rio_growflx_fusion_siscpf   => this%rvars(ir_growflx_fusion_siscpf)%r81d,  &
           rio_demorate_sisc           => this%rvars(ir_demorate_sisc)%r81d, &
           rio_promrate_sisc           => this%rvars(ir_promrate_sisc)%r81d, &
           rio_termcflux_cano_si       => this%rvars(ir_termcflux_cano_si)%r81d, &
           rio_termcflux_usto_si       => this%rvars(ir_termcflux_usto_si)%r81d, &
           rio_democflux_si            => this%rvars(ir_democflux_si)%r81d, &
           rio_promcflux_si            => this%rvars(ir_promcflux_si)%r81d, &
           rio_imortcflux_si           => this%rvars(ir_imortcflux_si)%r81d, &
           rio_fmortcflux_cano_si      => this%rvars(ir_fmortcflux_cano_si)%r81d, &
           rio_fmortcflux_usto_si      => this%rvars(ir_fmortcflux_usto_si)%r81d)


       totalCohorts = 0
       
       ! ---------------------------------------------------------------------------------
       ! Flush arrays to values defined by %flushval (see registry entry in
       ! subroutine define_history_vars()
       ! ---------------------------------------------------------------------------------
       call this%flush_rvars(nc)
       
       do s = 1,nsites
          
          ! Calculate the offsets
          ! fcolumn is the global column index of the current site.
          ! For the first site, if that site aligns with the first column index
          ! in the clump, than the offset should be be equal to begCohort
          
          io_idx_si      = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)

          io_idx_co      = io_idx_co_1st
          io_idx_pa_ib   = io_idx_co_1st
          io_idx_si_wmem = io_idx_co_1st
          io_idx_si_vtmem = io_idx_co_1st


          ! Hydraulics counters  lyr = hydraulic layer, shell = rhizosphere shell
          io_idx_si_lyr_shell = io_idx_co_1st
          io_idx_si_scpf = io_idx_co_1st
          io_idx_si_sc   = io_idx_co_1st
          io_idx_si_capf = io_idx_co_1st
          io_idx_si_cacls= io_idx_co_1st
          
          ! recruitment rate
          do i_pft = 1,numpft
             rio_recrate_sift(io_idx_co_1st+i_pft-1)   = sites(s)%recruitment_rate(i_pft)
          end do

          do el = 1, num_elements

             io_idx_si_cwd = io_idx_co_1st
             io_idx_si_pft = io_idx_co_1st

             do i_cwd=1,ncwd
                this%rvars(ir_cwdagin_flxdg+el-1)%r81d(io_idx_si_cwd) = sites(s)%flux_diags(el)%cwd_ag_input(i_cwd)
                this%rvars(ir_cwdbgin_flxdg+el-1)%r81d(io_idx_si_cwd) = sites(s)%flux_diags(el)%cwd_bg_input(i_cwd)
                io_idx_si_cwd = io_idx_si_cwd + 1
             end do
             
             do i_pft=1,numpft
                this%rvars(ir_leaflittin_flxdg+el-1)%r81d(io_idx_si_pft) = sites(s)%flux_diags(el)%leaf_litter_input(i_pft)
                this%rvars(ir_rootlittin_flxdg+el-1)%r81d(io_idx_si_pft) = sites(s)%flux_diags(el)%root_litter_input(i_pft)
                io_idx_si_pft = io_idx_si_pft + 1
             end do

             this%rvars(ir_oldstock_mbal+el-1)%r81d(io_idx_si) = sites(s)%mass_balance(el)%old_stock
             this%rvars(ir_errfates_mbal+el-1)%r81d(io_idx_si) = sites(s)%mass_balance(el)%err_fates

          end do


          ! canopy spread term
          rio_spread_si(io_idx_si)   = sites(s)%spread
          
          cpatch => sites(s)%oldest_patch
          
          ! new column, reset num patches
          patchespersite = 0
          
          do while(associated(cpatch))
             
             ! found patch, increment
             patchespersite = patchespersite + 1
             
             ccohort => cpatch%shortest
             
             ! new patch, reset num cohorts
             cohortsperpatch = 0
             
             do while(associated(ccohort))
                
                ! found cohort, increment
                cohortsperpatch = cohortsperpatch + 1
                totalCohorts    = totalCohorts + 1
             
                if ( debug ) then
                   write(fates_log(),*) 'CLTV io_idx_co ', io_idx_co
                   write(fates_log(),*) 'CLTV lowerbound ', lbound(rio_npp_acc_co,1) 
                   write(fates_log(),*) 'CLTV upperbound  ', ubound(rio_npp_acc_co,1)
                endif


                ! Fill output arrays of PRT variables
                ! We just loop through the objects, and reference our members relative
                ! the base index of the PRT variables
                ! -----------------------------------------------------------------------

                ir_prt_var = ir_prt_base
                do i_var = 1, prt_global%num_vars
                   do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos
                      
                      ir_prt_var = ir_prt_var + 1
                      this%rvars(ir_prt_var)%r81d(io_idx_co) = &
                            ccohort%prt%variables(i_var)%val(i_pos)

                      ir_prt_var = ir_prt_var + 1
                      this%rvars(ir_prt_var)%r81d(io_idx_co) = &
                            ccohort%prt%variables(i_var)%turnover(i_pos)
                      
                      ir_prt_var = ir_prt_var + 1
                      this%rvars(ir_prt_var)%r81d(io_idx_co) = &
                            ccohort%prt%variables(i_var)%net_alloc(i_pos)

                      ir_prt_var = ir_prt_var + 1
                      this%rvars(ir_prt_var)%r81d(io_idx_co) = &
                            ccohort%prt%variables(i_var)%burned(i_pos)
                      
                   end do
                end do

                
                if(hlm_use_planthydro==itrue)then
                   
                   ! Load the water contents
                   call this%SetCohortRealVector(ccohort%co_hydr%th_ag,n_hypool_ag, &
                                                 ir_hydro_th_ag_covec,io_idx_co)
                   call this%SetCohortRealVector(ccohort%co_hydr%th_aroot,sites(s)%si_hydr%nlevrhiz, &
                                                 ir_hydro_th_aroot_covec,io_idx_co)

                   this%rvars(ir_hydro_th_troot)%r81d(io_idx_co) = ccohort%co_hydr%th_troot

                   ! Load the error terms
                   call this%setCohortRealVector(ccohort%co_hydr%errh2o_growturn_ag, &
                                                 n_hypool_ag, &
                                                 ir_hydro_err_growturn_ag_covec,io_idx_co)
                   
                   this%rvars(ir_hydro_err_growturn_aroot)%r81d(io_idx_co) = &
                        ccohort%co_hydr%errh2o_growturn_aroot
                        
                   this%rvars(ir_hydro_err_growturn_troot)%r81d(io_idx_co) = &
                        ccohort%co_hydr%errh2o_growturn_troot
                        

                end if

                rio_canopy_layer_co(io_idx_co) = ccohort%canopy_layer
                rio_canopy_layer_yesterday_co(io_idx_co) = ccohort%canopy_layer_yesterday
                rio_canopy_trim_co(io_idx_co)  = ccohort%canopy_trim
                rio_seed_prod_co(io_idx_co)    = ccohort%seed_prod
                rio_size_class_lasttimestep(io_idx_co) = ccohort%size_class_lasttimestep
                rio_dbh_co(io_idx_co)          = ccohort%dbh
                rio_coage_co(io_idx_co)        = ccohort%coage
                rio_height_co(io_idx_co)       = ccohort%hite
                rio_laimemory_co(io_idx_co)    = ccohort%laimemory
                rio_sapwmemory_co(io_idx_co)   = ccohort%sapwmemory
                rio_structmemory_co(io_idx_co) = ccohort%structmemory
                rio_g_sb_laweight_co(io_idx_co)= ccohort%g_sb_laweight

                rio_nplant_co(io_idx_co)       = ccohort%n
                rio_gpp_acc_co(io_idx_co)      = ccohort%gpp_acc
                rio_npp_acc_co(io_idx_co)      = ccohort%npp_acc
                rio_resp_acc_co(io_idx_co)     = ccohort%resp_acc
                rio_gpp_acc_hold_co(io_idx_co) = ccohort%gpp_acc_hold
                rio_resp_acc_hold_co(io_idx_co) = ccohort%resp_acc_hold
                rio_npp_acc_hold_co(io_idx_co) = ccohort%npp_acc_hold

                rio_bmort_co(io_idx_co)        = ccohort%bmort
                rio_hmort_co(io_idx_co)        = ccohort%hmort
                rio_cmort_co(io_idx_co)        = ccohort%cmort
                rio_smort_co(io_idx_co)        = ccohort%smort
                rio_asmort_co(io_idx_co)       = ccohort%asmort
                rio_frmort_co(io_idx_co)       = ccohort%frmort                

                !Logging
                rio_lmort_direct_co(io_idx_co)       = ccohort%lmort_direct
                rio_lmort_collateral_co(io_idx_co)   = ccohort%lmort_collateral
                rio_lmort_infra_co(io_idx_co)        = ccohort%lmort_infra

                rio_ddbhdt_co(io_idx_co)       = ccohort%ddbhdt
                rio_resp_tstep_co(io_idx_co)   = ccohort%resp_tstep
                rio_pft_co(io_idx_co)          = ccohort%pft
                rio_status_co(io_idx_co)       = ccohort%status_coh
                if ( ccohort%isnew ) then
                   rio_isnew_co(io_idx_co)     = new_cohort
                else
                   rio_isnew_co(io_idx_co)     = old_cohort
                endif
                
                if ( debug ) then
                   write(fates_log(),*) 'CLTV offsetNumCohorts II ',io_idx_co, &
                         cohortsperpatch
                endif
             
                io_idx_co = io_idx_co + 1
                
                ccohort => ccohort%taller
                
             enddo ! ccohort do while
             
             !
             ! deal with patch level fields here
             !
             rio_livegrass_pa(io_idx_co_1st)   = cpatch%livegrass
             rio_age_pa(io_idx_co_1st)         = cpatch%age
             rio_patchdistturbcat_pa(io_idx_co_1st)   = cpatch%anthro_disturbance_label
             rio_agesinceanthrodist_pa(io_idx_co_1st) = cpatch%age_since_anthro_disturbance
             rio_area_pa(io_idx_co_1st)        = cpatch%area
             
             ! set cohorts per patch for IO
             rio_ncohort_pa( io_idx_co_1st )   = cohortsperpatch
             
             ! Set zenith angle info
             if ( cpatch%solar_zenith_flag ) then
                rio_solar_zenith_flag_pa(io_idx_co_1st)     = itrue
             else
                rio_solar_zenith_flag_pa(io_idx_co_1st)     = ifalse
             endif
             rio_solar_zenith_angle_pa( io_idx_co_1st) = cpatch%solar_zenith_angle

             if ( debug ) then
                write(fates_log(),*) 'offsetNumCohorts III ' &
                      ,io_idx_co,cohortsperpatch
             endif

             ! --------------------------------------------------------------------------
             ! Send litter to the restart arrays
             ! Each element has its own variable, so we have to make sure 
             ! we keep re-setting this 
             ! --------------------------------------------------------------------------

             do el = 0, num_elements-1
                 
                 io_idx_pa_pft  = io_idx_co_1st
                 io_idx_pa_cwd  = io_idx_co_1st
                 io_idx_pa_cwsl = io_idx_co_1st
                 io_idx_pa_dcsl = io_idx_co_1st
                 io_idx_pa_dc   = io_idx_co_1st
                 
                 litt => cpatch%litter(el+1)

                 do i = 1,numpft
                    this%rvars(ir_seed_litt+el)%r81d(io_idx_pa_pft) = litt%seed(i)
                    this%rvars(ir_seedgerm_litt+el)%r81d(io_idx_pa_pft) = litt%seed_germ(i)
                    io_idx_pa_pft = io_idx_pa_pft + 1
                 end do


                 do i = 1,ndcmpy
                     this%rvars(ir_leaf_litt+el)%r81d(io_idx_pa_dc) = litt%leaf_fines(i)
                     io_idx_pa_dc = io_idx_pa_dc + 1
                     do ilyr=1,sites(s)%nlevsoil
                         this%rvars(ir_fnrt_litt+el)%r81d(io_idx_pa_dcsl) = litt%root_fines(i,ilyr)
                         io_idx_pa_dcsl = io_idx_pa_dcsl + 1
                     end do
                 end do
                 
                 do i = 1,ncwd
                     this%rvars(ir_agcwd_litt+el)%r81d(io_idx_pa_cwd) = litt%ag_cwd(i)
                     io_idx_pa_cwd = io_idx_pa_cwd + 1
                     do ilyr=1,sites(s)%nlevsoil
                         this%rvars(ir_bgcwd_litt+el)%r81d(io_idx_pa_cwsl) = litt%bg_cwd(i,ilyr)
                         io_idx_pa_cwsl = io_idx_pa_cwsl + 1
                     end do
                 end do

             end do

             
             do i = 1,maxSWb
                rio_gnd_alb_dif_pasb(io_idx_pa_ib) = cpatch%gnd_alb_dif(i)
                rio_gnd_alb_dir_pasb(io_idx_pa_ib) = cpatch%gnd_alb_dir(i)
                io_idx_pa_ib = io_idx_pa_ib + 1
             end do

             ! Set the first cohort index to the start of the next patch, increment
             ! by the maximum number of cohorts per patch
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch
             
             ! reset counters so that they are all advanced evenly.
             io_idx_pa_pft  = io_idx_co_1st
             io_idx_pa_cwd  = io_idx_co_1st
             io_idx_pa_ib   = io_idx_co_1st
             io_idx_co      = io_idx_co_1st
             
             if ( debug ) then
                write(fates_log(),*) 'CLTV io_idx_co_1st ', io_idx_co_1st
                write(fates_log(),*) 'CLTV numCohort ', cohortsperpatch
                write(fates_log(),*) 'CLTV totalCohorts ', totalCohorts
             end if
             
             cpatch => cpatch%younger
             
          enddo ! cpatch do while
          

          ! Fill the site level diagnostics arrays
          do i_scls = 1, nlevsclass
             
             do i_pft = 1, numpft
             
                rio_fmortrate_cano_siscpf(io_idx_si_scpf)  = sites(s)%fmort_rate_canopy(i_scls, i_pft)
                rio_fmortrate_usto_siscpf(io_idx_si_scpf)  = sites(s)%fmort_rate_ustory(i_scls, i_pft)
                rio_imortrate_siscpf(io_idx_si_scpf)       = sites(s)%imort_rate(i_scls, i_pft)
                rio_fmortrate_crown_siscpf(io_idx_si_scpf) = sites(s)%fmort_rate_crown(i_scls, i_pft)
                rio_fmortrate_cambi_siscpf(io_idx_si_scpf) = sites(s)%fmort_rate_cambial(i_scls, i_pft)
                rio_termnindiv_cano_siscpf(io_idx_si_scpf) = sites(s)%term_nindivs_canopy(i_scls,i_pft)
                rio_termnindiv_usto_siscpf(io_idx_si_scpf) = sites(s)%term_nindivs_ustory(i_scls,i_pft)
                rio_growflx_fusion_siscpf(io_idx_si_scpf)  = sites(s)%growthflux_fusion(i_scls, i_pft)
                
                io_idx_si_scpf = io_idx_si_scpf + 1
             end do

             rio_demorate_sisc(io_idx_si_sc) = sites(s)%demotion_rate(i_scls)
             rio_promrate_sisc(io_idx_si_sc) = sites(s)%promotion_rate(i_scls)
                
             io_idx_si_sc = io_idx_si_sc + 1
          end do
          
          rio_termcflux_cano_si(io_idx_si)  = sites(s)%term_carbonflux_canopy
          rio_termcflux_usto_si(io_idx_si)  = sites(s)%term_carbonflux_ustory
          rio_democflux_si(io_idx_si)       = sites(s)%demotion_carbonflux
          rio_promcflux_si(io_idx_si)       = sites(s)%promotion_carbonflux
          rio_imortcflux_si(io_idx_si)      = sites(s)%imort_carbonflux
          rio_fmortcflux_cano_si(io_idx_si) = sites(s)%fmort_carbonflux_canopy
          rio_fmortcflux_usto_si(io_idx_si) = sites(s)%fmort_carbonflux_ustory

          rio_cd_status_si(io_idx_si)    = sites(s)%cstatus
          rio_dd_status_si(io_idx_si)    = sites(s)%dstatus
          rio_nchill_days_si(io_idx_si)  = sites(s)%nchilldays
          rio_ncold_days_si(io_idx_si)   = sites(s)%ncolddays
          rio_leafondate_si(io_idx_si)   = sites(s)%cleafondate
          rio_leafoffdate_si(io_idx_si)  = sites(s)%cleafoffdate

          rio_dleafondate_si(io_idx_si)  = sites(s)%dleafondate
          rio_dleafoffdate_si(io_idx_si) = sites(s)%dleafoffdate
          rio_acc_ni_si(io_idx_si)       = sites(s)%acc_NI
          rio_gdd_si(io_idx_si)          = sites(s)%grow_deg_days 
          
          ! Accumulated trunk product
          rio_trunk_product_si(io_idx_si) = sites(s)%resources_management%trunk_product_site
          ! set numpatches for this column

          rio_npatch_si(io_idx_si)  = patchespersite
          
          do i = 1,numWaterMem ! numWaterMem currently 10
             rio_watermem_siwm( io_idx_si_wmem ) = sites(s)%water_memory(i)
             io_idx_si_wmem = io_idx_si_wmem + 1
          end do

          do i = 1, num_vegtemp_mem
             rio_vegtempmem_sitm( io_idx_si_vtmem ) = sites(s)%vegtemp_memory(i)
             io_idx_si_vtmem = io_idx_si_vtmem + 1
          end do

          ! -----------------------------------------------------------------------------
          ! Set site-level hydraulics arrays
          ! -----------------------------------------------------------------------------

          if(hlm_use_planthydro==itrue)then

             ! No associate statements because there is no gaurantee these
             ! are allocated

             this%rvars(ir_hydro_recruit_si)%r81d(io_idx_si) = sites(s)%si_hydr%h2oveg_recruit
             this%rvars(ir_hydro_dead_si)%r81d(io_idx_si) = sites(s)%si_hydr%h2oveg_dead
             this%rvars(ir_hydro_growturn_err_si)%r81d(io_idx_si) = sites(s)%si_hydr%h2oveg_growturn_err
             this%rvars(ir_hydro_pheno_err_si)%r81d(io_idx_si) = sites(s)%si_hydr%h2oveg_pheno_err
             this%rvars(ir_hydro_hydro_err_si)%r81d(io_idx_si) = sites(s)%si_hydr%h2oveg_hydro_err

             ! Hydraulics counters  lyr = hydraulic layer, shell = rhizosphere shell
             do i = 1, sites(s)%si_hydr%nlevrhiz
                ! Loop shells
                do k = 1, nshell
                   this%rvars(ir_hydro_liqvol_shell_si)%r81d(io_idx_si_lyr_shell) = &
                        sites(s)%si_hydr%h2osoi_liqvol_shell(i,k)
                   io_idx_si_lyr_shell = io_idx_si_lyr_shell + 1
                end do
             end do
          end if
          
       enddo
       
       if ( debug ) then
          write(fates_log(),*) 'CLTV total cohorts ',totalCohorts
       end if
       
       return
     end associate
   end subroutine set_restart_vectors

   ! ====================================================================================

   subroutine create_patchcohort_structure(this, nc, nsites, sites, bc_in) 

     ! ----------------------------------------------------------------------------------
     ! This subroutine takes a peak at the restart file to determine how to allocate
     ! memory for the state structure, and then makes those allocations. This
     ! subroutine is called prior to the transfer of the restart vectors into the
     ! linked-list state structure.
     ! ---------------------------------------------------------------------------------

     use EDTypesMod,           only : ed_site_type
     use EDTypesMod,           only : ed_cohort_type
     use EDTypesMod,           only : ed_patch_type
     use EDTypesMod,           only : maxSWb
     use FatesInterfaceTypesMod,    only : fates_maxElementsPerPatch
     
     use EDTypesMod,           only : maxpft
     use EDTypesMod,           only : area
     use EDPatchDynamicsMod,   only : zero_patch
     use EDInitMod,            only : zero_site
     use EDInitMod,            only : init_site_vars
     use EDPatchDynamicsMod,   only : create_patch
     use EDPftvarcon,          only : EDPftvarcon_inst
     use FatesAllometryMod,    only : h2d_allom
     

     ! !ARGUMENTS:
     class(fates_restart_interface_type) , intent(inout) :: this
     integer                     , intent(in)            :: nc
     integer                     , intent(in)            :: nsites
     type(ed_site_type)          , intent(inout), target :: sites(nsites)
     type(bc_in_type)            , intent(in)            :: bc_in(nsites)

     ! local variables
     
     type(ed_patch_type) , pointer     :: newp
     type(ed_cohort_type), pointer     :: new_cohort
     type(ed_cohort_type), pointer     :: prev_cohort
     integer                           :: cohortstatus
     integer                           :: s             ! site index
     integer                           :: idx_pa        ! local patch index
     integer                           :: io_idx_si     ! global site index in IO vector
     integer                           :: io_idx_co_1st ! global cohort index in IO vector
     real(r8)                          :: site_spread   ! site sprea dummy var (0-1)
     integer                           :: fto
     integer                           :: ft
     integer                           :: el            ! element loop counter
     integer, parameter                :: recruitstatus = 0

     ! ----------------------------------------------------------------------------------
     ! We really only need the counts for the number of patches per site
     ! and the number of cohorts per patch. These values tell us how much
     ! space to allocate.
     ! ----------------------------------------------------------------------------------
     
     associate( rio_npatch_si  => this%rvars(ir_npatch_si)%int1d , &
                rio_ncohort_pa => this%rvars(ir_ncohort_pa)%int1d )
            
       do s = 1,nsites
          
          io_idx_si  = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)

          call init_site_vars( sites(s), bc_in(s) )
          call zero_site( sites(s) )

          if ( rio_npatch_si(io_idx_si)<0 .or. rio_npatch_si(io_idx_si) > 10000 ) then
             write(fates_log(),*) 'a column was expected to contain a valid number of patches'
             write(fates_log(),*) '0 is a valid number, but this column seems uninitialized',rio_npatch_si(io_idx_si)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       
          ! Initialize the site pointers to null
          sites(s)%youngest_patch         => null()                 
          sites(s)%oldest_patch           => null()

          do idx_pa = 1,rio_npatch_si(io_idx_si)

             if ( debug ) then
                write(fates_log(),*) 'create patch ',idx_pa
                write(fates_log(),*) 'idx_pa 1-cohortsperpatch : ', rio_ncohort_pa( io_idx_co_1st )
             end if
             
             ! create patch
             allocate(newp)    
             
             ! make new patch
             call create_patch(sites(s), newp, fates_unset_r8, fates_unset_r8, primaryforest )

             ! Initialize the litter pools to zero, these
             ! pools will be populated by looping over the existing patches
             ! and transfering in mass
             do el=1,num_elements
                call newp%litter(el)%InitConditions(init_leaf_fines=fates_unset_r8, &
                     init_root_fines=fates_unset_r8, &
                     init_ag_cwd=fates_unset_r8, &
                     init_bg_cwd=fates_unset_r8, &
                     init_seed=fates_unset_r8, &
                     init_seed_germ=fates_unset_r8)
             end do
             
             ! give this patch a unique patch number
             newp%patchno = idx_pa


             ! Iterate over the number of cohorts
             ! the file says are associated with this patch
             ! we are just allocating space here, so we do 
             ! a simple list filling routine
             
             newp%tallest  => null()
             newp%shortest => null()
             prev_cohort   => null()

             do fto = 1, rio_ncohort_pa( io_idx_co_1st )

                allocate(new_cohort)
                call nan_cohort(new_cohort)  
                call zero_cohort(new_cohort)
                new_cohort%patchptr => newp

                ! If this is the first in the list, it is tallest
                if (.not.associated(newp%tallest)) then
                   newp%tallest => new_cohort
                endif
                
                ! Every cohort's taller is the one that came before
                ! (unless it is first)
                if(associated(prev_cohort)) then
                   new_cohort%taller   => prev_cohort
                   prev_cohort%shorter => new_cohort
                end if

                ! Every cohort added takes over as shortest
                newp%shortest => new_cohort

                ! Initialize the PARTEH object and point to the
                ! correct boundary condition fields
                new_cohort%prt => null()
                call InitPRTObject(new_cohort%prt)
                call InitPRTBoundaryConditions(new_cohort)
                
                
                ! Allocate hydraulics arrays
                if( hlm_use_planthydro.eq.itrue ) then
                   call InitHydrCohort(sites(s),new_cohort)
                end if

                ! Update the previous
                prev_cohort => new_cohort
                
             enddo ! ends loop over fto
             
             !
             ! insert this patch with cohorts into the site pointer.  At this
             ! point just insert the new patch in the youngest position
             !
             if (idx_pa == 1) then ! nothing associated yet. first patch is pointed to by youngest and oldest
                
                if ( debug ) write(fates_log(),*) 'idx_pa = 1 ',idx_pa
                
                sites(s)%youngest_patch         => newp                   
                sites(s)%oldest_patch           => newp                        
                sites(s)%youngest_patch%younger => null()
                sites(s)%youngest_patch%older   => null()
                sites(s)%oldest_patch%younger   => null()
                sites(s)%oldest_patch%older     => null()
                
             else if (idx_pa == 2) then ! add second patch to list
                
                if ( debug ) write(fates_log(),*) 'idx_pa = 2 ',idx_pa
                
                sites(s)%youngest_patch         => newp
                sites(s)%youngest_patch%younger => null()
                sites(s)%youngest_patch%older   => sites(s)%oldest_patch
                sites(s)%oldest_patch%younger   => sites(s)%youngest_patch
                sites(s)%oldest_patch%older     => null()

             else ! more than 2 patches, insert patch into youngest slot
                
                if ( debug ) write(fates_log(),*) 'idx_pa > 2 ',idx_pa
                
                newp%older                      => sites(s)%youngest_patch
                sites(s)%youngest_patch%younger => newp
                newp%younger                    => null()
                sites(s)%youngest_patch         => newp
                
             endif
             
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch

          enddo ! ends loop over idx_pa

       enddo ! ends loop over s
       
     end associate
   end subroutine create_patchcohort_structure
   
   ! ====================================================================================

   subroutine get_restart_vectors(this, nc, nsites, sites)

     use EDTypesMod, only : ed_site_type
     use EDTypesMod, only : ed_cohort_type
     use EDTypesMod, only : ed_patch_type
     use EDTypesMod, only : maxSWb
     use FatesInterfaceTypesMod, only : numpft
     use FatesInterfaceTypesMod, only : fates_maxElementsPerPatch
     use EDTypesMod, only : numWaterMem
     use EDTypesMod, only : num_vegtemp_mem
     use FatesSizeAgeTypeIndicesMod, only : get_age_class_index

     ! !ARGUMENTS:
     class(fates_restart_interface_type) , intent(inout) :: this
     integer                     , intent(in)            :: nc
     integer                     , intent(in)            :: nsites
     type(ed_site_type)          , intent(inout), target :: sites(nsites)


     ! locals
     ! ----------------------------------------------------------------------------------
     ! LL pointers
     type(ed_patch_type),pointer  :: cpatch      ! current patch
     type(ed_cohort_type),pointer :: ccohort     ! current cohort
     type(litter_type), pointer   :: litt        ! litter object on the current patch
     ! loop indices
     integer :: s, i, j, k

     ! ----------------------------------------------------------------------------------
     ! The following group of integers indicate the positional index (idx)
     ! of variables at different scales inside the I/O arrays (io)
     ! Keep in mind that many of these variables have a composite dimension
     ! at the patch scale.  To hold this memory, we borrow the cohort
     ! vector.  Thus the head of each array points to the first cohort
     ! of each patch. "io_idx_co_1st"
     ! ----------------------------------------------------------------------------------
     integer  :: io_idx_si      ! site index
     integer  :: io_idx_co_1st  ! 1st cohort of each patch
     integer  :: io_idx_co      ! cohort index
     integer  :: io_idx_pa_pft  ! each pft within each patch (pa_pft)
     integer  :: io_idx_pa_cwd  ! each cwd class within each patch (pa_cwd)
     integer  :: io_idx_pa_cwsl ! each cwd x soil layer
     integer  :: io_idx_pa_dcsl ! each decomposability x soil layer
     integer  :: io_idx_pa_dc   ! each decomposability index 
     integer  :: io_idx_pa_ib   ! each SW radiation band per patch (pa_ib)
     integer  :: io_idx_si_wmem ! each water memory class within each site
     integer  :: io_idx_si_vtmem ! counter for vegetation temp memory
     integer  :: io_idx_si_lyr_shell ! site - layer x shell index
     integer  :: io_idx_si_scpf ! each size-class x pft index within site
     integer  :: io_idx_si_sc   ! each size-class index within site
     integer  :: io_idx_si_cacls ! each coage class index within site
     integer  :: io_idx_si_capf ! each cohort age class x pft index within site
     integer  :: io_idx_si_cwd
     integer  :: io_idx_si_pft

     ! Some counters (for checking mostly)
     integer  :: totalcohorts   ! total cohort count on this thread (diagnostic)
     integer  :: patchespersite   ! number of patches per site
     integer  :: cohortsperpatch  ! number of cohorts per patch 
     integer  :: el               ! loop counter for elements
     integer  :: nlevsoil         ! number of soil layers
     integer  :: ilyr             ! soil layer loop counter
     integer  :: ir_prt_var       ! loop counter for var x position
     integer  :: i_cwd            ! loop counter for cwd
     integer  :: i_var            ! loop counter for PRT variables
     integer  :: i_pos            ! loop counter for discrete PRT positions
     integer  :: i_pft            ! loop counter for pft
     integer  :: i_scls           ! loop counter for size-clas
     integer  :: i_cacls          ! loop counter for cohort age class

     associate( rio_npatch_si         => this%rvars(ir_npatch_si)%int1d, &
          rio_cd_status_si            => this%rvars(ir_cd_status_si)%int1d, &
          rio_dd_status_si            => this%rvars(ir_dd_status_si)%int1d, &
          rio_nchill_days_si          => this%rvars(ir_nchill_days_si)%int1d, &
          rio_ncold_days_si           => this%rvars(ir_ncold_days_si)%int1d, &
          rio_leafondate_si           => this%rvars(ir_leafondate_si)%int1d, &
          rio_leafoffdate_si          => this%rvars(ir_leafoffdate_si)%int1d, &
          rio_dleafondate_si          => this%rvars(ir_dleafondate_si)%int1d, &
          rio_dleafoffdate_si         => this%rvars(ir_dleafoffdate_si)%int1d, &
          rio_acc_ni_si               => this%rvars(ir_acc_ni_si)%r81d, &
          rio_gdd_si                  => this%rvars(ir_gdd_si)%r81d, &
          rio_trunk_product_si        => this%rvars(ir_trunk_product_si)%r81d, &
          rio_ncohort_pa              => this%rvars(ir_ncohort_pa)%int1d, &
          rio_solar_zenith_flag_pa    => this%rvars(ir_solar_zenith_flag_pa)%int1d, &
          rio_solar_zenith_angle_pa   => this%rvars(ir_solar_zenith_angle_pa)%r81d, &
          rio_canopy_layer_co         => this%rvars(ir_canopy_layer_co)%int1d, &
          rio_canopy_layer_yesterday_co         => this%rvars(ir_canopy_layer_yesterday_co)%r81d, &
          rio_canopy_trim_co          => this%rvars(ir_canopy_trim_co)%r81d, &
          rio_seed_prod_co            => this%rvars(ir_seed_prod_co)%r81d, &
          rio_size_class_lasttimestep => this%rvars(ir_size_class_lasttimestep_co)%int1d, &
          rio_dbh_co                  => this%rvars(ir_dbh_co)%r81d, &
          rio_coage_co                => this%rvars(ir_coage_co)%r81d, & 
          rio_g_sb_laweight_co        => this%rvars(ir_g_sb_laweight_co)%r81d, &
          rio_height_co               => this%rvars(ir_height_co)%r81d, &
          rio_laimemory_co            => this%rvars(ir_laimemory_co)%r81d, &
          rio_sapwmemory_co           => this%rvars(ir_sapwmemory_co)%r81d, &
          rio_structmemory_co         => this%rvars(ir_structmemory_co)%r81d, &
          rio_nplant_co               => this%rvars(ir_nplant_co)%r81d, &
          rio_gpp_acc_co              => this%rvars(ir_gpp_acc_co)%r81d, &
          rio_npp_acc_co              => this%rvars(ir_npp_acc_co)%r81d, &
          rio_resp_acc_co             => this%rvars(ir_resp_acc_co)%r81d, &
          rio_gpp_acc_hold_co         => this%rvars(ir_gpp_acc_hold_co)%r81d, &
          rio_resp_acc_hold_co        => this%rvars(ir_resp_acc_hold_co)%r81d, &
          rio_npp_acc_hold_co         => this%rvars(ir_npp_acc_hold_co)%r81d, &
          rio_bmort_co                => this%rvars(ir_bmort_co)%r81d, &
          rio_hmort_co                => this%rvars(ir_hmort_co)%r81d, &
          rio_cmort_co                => this%rvars(ir_cmort_co)%r81d, &
          rio_smort_co                => this%rvars(ir_smort_co)%r81d, &
          rio_asmort_co               => this%rvars(ir_asmort_co)%r81d, &
          rio_frmort_co               => this%rvars(ir_frmort_co)%r81d, &
          rio_lmort_direct_co         => this%rvars(ir_lmort_direct_co)%r81d, &
          rio_lmort_collateral_co     => this%rvars(ir_lmort_collateral_co)%r81d, &
          rio_lmort_infra_co          => this%rvars(ir_lmort_infra_co)%r81d, &
          rio_ddbhdt_co               => this%rvars(ir_ddbhdt_co)%r81d, &
          rio_resp_tstep_co           => this%rvars(ir_resp_tstep_co)%r81d, &
          rio_pft_co                  => this%rvars(ir_pft_co)%int1d, &
          rio_status_co               => this%rvars(ir_status_co)%int1d, &
          rio_isnew_co                => this%rvars(ir_isnew_co)%int1d, &
          rio_gnd_alb_dif_pasb        => this%rvars(ir_gnd_alb_dif_pasb)%r81d, &
          rio_gnd_alb_dir_pasb        => this%rvars(ir_gnd_alb_dir_pasb)%r81d, &
          rio_spread_si               => this%rvars(ir_spread_si)%r81d, &
          rio_livegrass_pa            => this%rvars(ir_livegrass_pa)%r81d, &
          rio_age_pa                  => this%rvars(ir_age_pa)%r81d, &
          rio_patchdistturbcat_pa     => this%rvars(ir_patchdistturbcat_pa)%int1d,  &
          rio_agesinceanthrodist_pa   => this%rvars(ir_agesinceanthrodist_pa)%r81d, &
          rio_area_pa                 => this%rvars(ir_area_pa)%r81d, &
          rio_watermem_siwm           => this%rvars(ir_watermem_siwm)%r81d, &
          rio_vegtempmem_sitm         => this%rvars(ir_vegtempmem_sitm)%r81d, &
          rio_recrate_sift            => this%rvars(ir_recrate_sift)%r81d, &
          rio_fmortrate_cano_siscpf   => this%rvars(ir_fmortrate_cano_siscpf)%r81d, &
          rio_fmortrate_usto_siscpf   => this%rvars(ir_fmortrate_usto_siscpf)%r81d, &
          rio_imortrate_siscpf        => this%rvars(ir_imortrate_siscpf)%r81d, &
          rio_fmortrate_crown_siscpf  => this%rvars(ir_fmortrate_crown_siscpf)%r81d, &
          rio_fmortrate_cambi_siscpf  => this%rvars(ir_fmortrate_cambi_siscpf)%r81d, &
          rio_termnindiv_cano_siscpf  => this%rvars(ir_termnindiv_cano_siscpf)%r81d, &
          rio_termnindiv_usto_siscpf  => this%rvars(ir_termnindiv_usto_siscpf)%r81d, &
          rio_growflx_fusion_siscpf   => this%rvars(ir_growflx_fusion_siscpf)%r81d,  &
          rio_demorate_sisc           => this%rvars(ir_demorate_sisc)%r81d, &
          rio_promrate_sisc           => this%rvars(ir_promrate_sisc)%r81d, &
          rio_termcflux_cano_si       => this%rvars(ir_termcflux_cano_si)%r81d, &
          rio_termcflux_usto_si       => this%rvars(ir_termcflux_usto_si)%r81d, &
          rio_democflux_si            => this%rvars(ir_democflux_si)%r81d, &
          rio_promcflux_si            => this%rvars(ir_promcflux_si)%r81d, &
          rio_imortcflux_si           => this%rvars(ir_imortcflux_si)%r81d, &
          rio_fmortcflux_cano_si      => this%rvars(ir_fmortcflux_cano_si)%r81d, &
          rio_fmortcflux_usto_si      => this%rvars(ir_fmortcflux_usto_si)%r81d)
     

       totalcohorts = 0
     
       do s = 1,nsites
          
          io_idx_si      = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)
          
          io_idx_co      = io_idx_co_1st
          io_idx_pa_ib   = io_idx_co_1st
          io_idx_si_wmem = io_idx_co_1st
          io_idx_si_vtmem = io_idx_co_1st

          ! Hydraulics counters  lyr = hydraulic layer, shell = rhizosphere shell
          io_idx_si_lyr_shell = io_idx_co_1st

          io_idx_si_scpf = io_idx_co_1st
          io_idx_si_sc   = io_idx_co_1st
          io_idx_si_capf = io_idx_co_1st
          io_idx_si_cacls= io_idx_co_1st
          
          ! read seed_bank info(site-level, but PFT-resolved)
          do i_pft = 1,numpft 
             sites(s)%recruitment_rate(i_pft) = rio_recrate_sift(io_idx_co_1st+i_pft-1)
          enddo


          ! Mass balance and diagnostics across elements at the site level
          do el = 1, num_elements

             io_idx_si_cwd = io_idx_co_1st
             io_idx_si_pft = io_idx_co_1st

             do i_cwd=1,ncwd
                sites(s)%flux_diags(el)%cwd_ag_input(i_cwd) = this%rvars(ir_cwdagin_flxdg+el-1)%r81d(io_idx_si_cwd)
                sites(s)%flux_diags(el)%cwd_bg_input(i_cwd) = this%rvars(ir_cwdbgin_flxdg+el-1)%r81d(io_idx_si_cwd)
                io_idx_si_cwd = io_idx_si_cwd + 1
             end do
             
             do i_pft=1,numpft
                sites(s)%flux_diags(el)%leaf_litter_input(i_pft) = this%rvars(ir_leaflittin_flxdg+el-1)%r81d(io_idx_si_pft)
                sites(s)%flux_diags(el)%root_litter_input(i_pft) = this%rvars(ir_rootlittin_flxdg+el-1)%r81d(io_idx_si_pft)
                io_idx_si_pft = io_idx_si_pft + 1
            end do
            
            sites(s)%mass_balance(el)%old_stock = this%rvars(ir_oldstock_mbal+el-1)%r81d(io_idx_si)
            sites(s)%mass_balance(el)%err_fates = this%rvars(ir_errfates_mbal+el-1)%r81d(io_idx_si)

          end do

          sites(s)%spread = rio_spread_si(io_idx_si) 
          
          ! Perform a check on the number of patches per site
          patchespersite = 0
          
          cpatch => sites(s)%oldest_patch
          do while(associated(cpatch))
             
             patchespersite = patchespersite + 1
             
             ccohort => cpatch%shortest
             
             ! new patch, reset num cohorts
             cohortsperpatch = 0
             
             do while(associated(ccohort))        
                
                ! found cohort, increment
                cohortsperpatch  = cohortsperpatch    + 1
                totalcohorts     = totalcohorts + 1
                
                if ( debug ) then
                   write(fates_log(),*) 'CVTL io_idx_co ',io_idx_co
                endif

                ! Fill PRT state variables with array data
                ! We just loop through the objects, and reference our members relative
                ! the base index of the PRT variables
                ! -----------------------------------------------------------------------

                ir_prt_var = ir_prt_base
                do i_var = 1, prt_global%num_vars
                   do i_pos = 1, prt_global%state_descriptor(i_var)%num_pos 

                      ir_prt_var = ir_prt_var + 1
                      ccohort%prt%variables(i_var)%val(i_pos) = &
                            this%rvars(ir_prt_var)%r81d(io_idx_co)

                      ir_prt_var = ir_prt_var + 1
                      ccohort%prt%variables(i_var)%turnover(i_pos) = &
                            this%rvars(ir_prt_var)%r81d(io_idx_co)

                      ir_prt_var = ir_prt_var + 1
                      ccohort%prt%variables(i_var)%net_alloc(i_pos) = &
                            this%rvars(ir_prt_var)%r81d(io_idx_co)

                      ir_prt_var = ir_prt_var + 1
                      ccohort%prt%variables(i_var)%burned(i_pos) = &
                            this%rvars(ir_prt_var)%r81d(io_idx_co)                      
                   end do
                end do

                !ccohort%vcmax25top          
                !ccohort%jmax25top
                !ccohort%tpu25top          
                !ccohort%kp25top


                ccohort%canopy_layer = rio_canopy_layer_co(io_idx_co)
                ccohort%canopy_layer_yesterday = rio_canopy_layer_yesterday_co(io_idx_co)
                ccohort%canopy_trim  = rio_canopy_trim_co(io_idx_co)
                ccohort%seed_prod    = rio_seed_prod_co(io_idx_co)
                ccohort%size_class_lasttimestep = rio_size_class_lasttimestep(io_idx_co)
                ccohort%dbh          = rio_dbh_co(io_idx_co)
                ccohort%coage        = rio_coage_co(io_idx_co)
                ccohort%g_sb_laweight= rio_g_sb_laweight_co(io_idx_co)
                ccohort%hite         = rio_height_co(io_idx_co)
                ccohort%laimemory    = rio_laimemory_co(io_idx_co)
                ccohort%sapwmemory   = rio_sapwmemory_co(io_idx_co)
                ccohort%structmemory= rio_structmemory_co(io_idx_co)
                ccohort%n            = rio_nplant_co(io_idx_co)
                ccohort%gpp_acc      = rio_gpp_acc_co(io_idx_co)
                ccohort%npp_acc      = rio_npp_acc_co(io_idx_co)
                ccohort%resp_acc     = rio_resp_acc_co(io_idx_co)
                ccohort%gpp_acc_hold = rio_gpp_acc_hold_co(io_idx_co)
                ccohort%resp_acc_hold = rio_resp_acc_hold_co(io_idx_co)
                ccohort%npp_acc_hold = rio_npp_acc_hold_co(io_idx_co)

                ccohort%bmort        = rio_bmort_co(io_idx_co)
                ccohort%hmort        = rio_hmort_co(io_idx_co)
                ccohort%cmort        = rio_cmort_co(io_idx_co)
                ccohort%smort        = rio_smort_co(io_idx_co)
                ccohort%asmort       = rio_asmort_co(io_idx_co)
                ccohort%frmort        = rio_frmort_co(io_idx_co)

                !Logging
                ccohort%lmort_direct       = rio_lmort_direct_co(io_idx_co)
                ccohort%lmort_collateral   = rio_lmort_collateral_co(io_idx_co)
                ccohort%lmort_infra        = rio_lmort_infra_co(io_idx_co)

                ccohort%ddbhdt       = rio_ddbhdt_co(io_idx_co)
                ccohort%resp_tstep   = rio_resp_tstep_co(io_idx_co)
                ccohort%pft          = rio_pft_co(io_idx_co)
                ccohort%status_coh   = rio_status_co(io_idx_co)
                ccohort%isnew        = ( rio_isnew_co(io_idx_co) .eq. new_cohort )

                call UpdateCohortBioPhysRates(ccohort)


                ! Initialize Plant Hydraulics

                if(hlm_use_planthydro==itrue)then
                   
                   ! Load the water contents
                   call this%GetCohortRealVector(ccohort%co_hydr%th_ag,n_hypool_ag, &
                                                 ir_hydro_th_ag_covec,io_idx_co)
                   call this%GetCohortRealVector(ccohort%co_hydr%th_aroot,sites(s)%si_hydr%nlevrhiz, &
                                                 ir_hydro_th_aroot_covec,io_idx_co)
                   
                   ccohort%co_hydr%th_troot = this%rvars(ir_hydro_th_troot)%r81d(io_idx_co)
                   
                   call UpdatePlantPsiFTCFromTheta(ccohort,sites(s)%si_hydr)

                   
                   ccohort%co_hydr%errh2o_growturn_aroot = &
                        this%rvars(ir_hydro_err_growturn_aroot)%r81d(io_idx_co)
                   ccohort%co_hydr%errh2o_growturn_troot = &
                        this%rvars(ir_hydro_err_growturn_troot)%r81d(io_idx_co)

                   call this%GetCohortRealVector(ccohort%co_hydr%errh2o_growturn_ag, &
                                                 n_hypool_ag, &
                                                 ir_hydro_err_growturn_ag_covec,io_idx_co)
                end if
                
                io_idx_co = io_idx_co + 1
             
                ccohort => ccohort%taller
                
             enddo ! current cohort do while

             if(cohortsperpatch .ne. rio_ncohort_pa(io_idx_co_1st)) then
                write(fates_log(),*) 'Number of cohorts per patch during retrieval'
                write(fates_log(),*) 'does not match allocation'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if

             !
             ! deal with patch level fields here
             !
             cpatch%livegrass          = rio_livegrass_pa(io_idx_co_1st)
             cpatch%age                = rio_age_pa(io_idx_co_1st)
             cpatch%anthro_disturbance_label       = rio_patchdistturbcat_pa(io_idx_co_1st)
             cpatch%age_since_anthro_disturbance   = rio_agesinceanthrodist_pa(io_idx_co_1st)
             cpatch%area               = rio_area_pa(io_idx_co_1st)
             cpatch%age_class          = get_age_class_index(cpatch%age)

             ! Set zenith angle info
             cpatch%solar_zenith_flag  = ( rio_solar_zenith_flag_pa(io_idx_co_1st) .eq. itrue )
             cpatch%solar_zenith_angle = rio_solar_zenith_angle_pa(io_idx_co_1st)

             ! set cohorts per patch for IO
             
             if ( debug ) then
                write(fates_log(),*) 'CVTL III ' &
                     ,io_idx_co,cohortsperpatch
             endif
             
             ! --------------------------------------------------------------------------
             ! Pull litter from the restart arrays
             ! Each element has its own variable, so we have to make sure 
             ! we keep re-setting this 
             ! --------------------------------------------------------------------------

             do el = 0, num_elements-1
                 
                 io_idx_pa_pft  = io_idx_co_1st
                 io_idx_pa_cwd  = io_idx_co_1st
                 io_idx_pa_cwsl = io_idx_co_1st
                 io_idx_pa_dcsl = io_idx_co_1st
                 io_idx_pa_dc   = io_idx_co_1st

                 litt => cpatch%litter(el+1)
                 nlevsoil = size(litt%bg_cwd,dim=2)

                 do i = 1,numpft
                     litt%seed(i)       = this%rvars(ir_seed_litt+el)%r81d(io_idx_pa_pft)
                     litt%seed_germ(i)  = this%rvars(ir_seedgerm_litt+el)%r81d(io_idx_pa_pft)
                     io_idx_pa_pft      = io_idx_pa_pft + 1
                  end do

                  do i = 1,ndcmpy
                     litt%leaf_fines(i) = this%rvars(ir_leaf_litt+el)%r81d(io_idx_pa_dc)
                     io_idx_pa_dc       = io_idx_pa_dc + 1
                     do ilyr=1,nlevsoil
                         litt%root_fines(i,ilyr) = this%rvars(ir_fnrt_litt+el)%r81d(io_idx_pa_dcsl)
                         io_idx_pa_dcsl = io_idx_pa_dcsl + 1
                     end do
                 end do
                 
                 do i = 1,ncwd

                     litt%ag_cwd(i) = this%rvars(ir_agcwd_litt+el)%r81d(io_idx_pa_cwd)
                     io_idx_pa_cwd = io_idx_pa_cwd + 1
                     
                     do ilyr=1,nlevsoil
                         litt%bg_cwd(i,ilyr) = this%rvars(ir_bgcwd_litt+el)%r81d(io_idx_pa_cwsl)
                         io_idx_pa_cwsl = io_idx_pa_cwsl + 1
                     end do
                 end do

             end do

             do i = 1,maxSWb
                cpatch%gnd_alb_dif(i) = rio_gnd_alb_dif_pasb(io_idx_pa_ib)
                cpatch%gnd_alb_dir(i) = rio_gnd_alb_dir_pasb(io_idx_pa_ib)
                io_idx_pa_ib = io_idx_pa_ib + 1
             end do

             ! Now increment the position of the first cohort to that of the next
             ! patch
             
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch
             
             ! and max the number of allowed cohorts per patch
             io_idx_pa_pft  = io_idx_co_1st
             io_idx_pa_cwd  = io_idx_co_1st
             io_idx_pa_ib   = io_idx_co_1st
             io_idx_co      = io_idx_co_1st
             
             if ( debug ) then
                write(fates_log(),*) 'CVTL io_idx_co_1st ', io_idx_co_1st
                write(fates_log(),*) 'CVTL cohortsperpatch ', cohortsperpatch
                write(fates_log(),*) 'CVTL totalCohorts ', totalCohorts
             end if
             
             cpatch => cpatch%younger
             
          enddo ! patch do while
          
          if(patchespersite .ne. rio_npatch_si(io_idx_si)) then
             write(fates_log(),*) 'Number of patches per site during retrieval does not match allocation'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          do i = 1,numWaterMem
             sites(s)%water_memory(i) = rio_watermem_siwm( io_idx_si_wmem )
             io_idx_si_wmem = io_idx_si_wmem + 1
          end do

          do i = 1, num_vegtemp_mem
             sites(s)%vegtemp_memory(i) = rio_vegtempmem_sitm( io_idx_si_vtmem )
             io_idx_si_vtmem = io_idx_si_vtmem + 1
          end do

          ! -----------------------------------------------------------------------------
          ! Retrieve site-level hydraulics arrays
          ! Note that Hydraulics structures, their allocations, and the length
          ! declaration nlevsoi_hyd should be allocated early on when the code first
          ! allocates sites (before restart info), and when the soils layer is 
          ! first known.
          ! -----------------------------------------------------------------------------

          if(hlm_use_planthydro==itrue)then

             sites(s)%si_hydr%h2oveg_recruit      = this%rvars(ir_hydro_recruit_si)%r81d(io_idx_si)
             sites(s)%si_hydr%h2oveg_dead         = this%rvars(ir_hydro_dead_si)%r81d(io_idx_si)
             sites(s)%si_hydr%h2oveg_growturn_err = this%rvars(ir_hydro_growturn_err_si)%r81d(io_idx_si)
             sites(s)%si_hydr%h2oveg_pheno_err    = this%rvars(ir_hydro_pheno_err_si)%r81d(io_idx_si)
             sites(s)%si_hydr%h2oveg_hydro_err    = this%rvars(ir_hydro_hydro_err_si)%r81d(io_idx_si)

             ! Hydraulics counters  lyr = hydraulic layer, shell = rhizosphere shell
             do i = 1, sites(s)%si_hydr%nlevrhiz
                ! Loop shells
                do k = 1, nshell
                   sites(s)%si_hydr%h2osoi_liqvol_shell(i,k) = &
                        this%rvars(ir_hydro_liqvol_shell_si)%r81d(io_idx_si_lyr_shell)
                   io_idx_si_lyr_shell = io_idx_si_lyr_shell + 1
                end do
             end do

          end if
          

          ! Fill the site level diagnostics arrays
          do i_scls = 1,nlevsclass
             
             do i_pft = 1, numpft
             
                sites(s)%fmort_rate_canopy(i_scls, i_pft)  = rio_fmortrate_cano_siscpf(io_idx_si_scpf)
                sites(s)%fmort_rate_ustory(i_scls, i_pft)  = rio_fmortrate_usto_siscpf(io_idx_si_scpf)
                sites(s)%imort_rate(i_scls, i_pft)         = rio_imortrate_siscpf(io_idx_si_scpf)
                sites(s)%fmort_rate_crown(i_scls, i_pft)   = rio_fmortrate_crown_siscpf(io_idx_si_scpf)
                sites(s)%fmort_rate_cambial(i_scls, i_pft) = rio_fmortrate_cambi_siscpf(io_idx_si_scpf) 
                sites(s)%term_nindivs_canopy(i_scls,i_pft) = rio_termnindiv_cano_siscpf(io_idx_si_scpf)
                sites(s)%term_nindivs_ustory(i_scls,i_pft) = rio_termnindiv_usto_siscpf(io_idx_si_scpf)
                sites(s)%growthflux_fusion(i_scls, i_pft)  = rio_growflx_fusion_siscpf(io_idx_si_scpf)

                io_idx_si_scpf = io_idx_si_scpf + 1
             end do

             sites(s)%demotion_rate(i_scls)  = rio_demorate_sisc(io_idx_si_sc)
             sites(s)%promotion_rate(i_scls) = rio_promrate_sisc(io_idx_si_sc)
                
             io_idx_si_sc = io_idx_si_sc + 1
          end do

          sites(s)%term_carbonflux_canopy   = rio_termcflux_cano_si(io_idx_si)
          sites(s)%term_carbonflux_ustory   = rio_termcflux_usto_si(io_idx_si)
          sites(s)%demotion_carbonflux      = rio_democflux_si(io_idx_si)
          sites(s)%promotion_carbonflux     = rio_promcflux_si(io_idx_si)
          sites(s)%imort_carbonflux         = rio_imortcflux_si(io_idx_si)
          sites(s)%fmort_carbonflux_canopy  = rio_fmortcflux_cano_si(io_idx_si)
          sites(s)%fmort_carbonflux_ustory  = rio_fmortcflux_usto_si(io_idx_si)

          
          ! Site level phenology status flags

          sites(s)%cstatus        = rio_cd_status_si(io_idx_si)
          sites(s)%dstatus        = rio_dd_status_si(io_idx_si)
          sites(s)%nchilldays     = rio_nchill_days_si(io_idx_si)
          sites(s)%ncolddays      = rio_ncold_days_si(io_idx_si)
          sites(s)%cleafondate    = rio_leafondate_si(io_idx_si)
          sites(s)%cleafoffdate   = rio_leafoffdate_si(io_idx_si)
          sites(s)%dleafondate    = rio_dleafondate_si(io_idx_si)
          sites(s)%dleafoffdate   = rio_dleafoffdate_si(io_idx_si)
          sites(s)%acc_NI         = rio_acc_ni_si(io_idx_si)
          sites(s)%grow_deg_days  = rio_gdd_si(io_idx_si)

          sites(s)%resources_management%trunk_product_site = rio_trunk_product_si(io_idx_si)

       end do

       if ( debug ) then
          write(fates_log(),*) 'CVTL total cohorts ',totalCohorts
       end if
       
     end associate
   end subroutine get_restart_vectors
   
   ! ====================================================================================

   subroutine update_3dpatch_radiation(this, nsites, sites, bc_out)

     ! -------------------------------------------------------------------------
     ! This subroutine populates output boundary conditions related to radiation
     ! called upon restart reads.
     ! -------------------------------------------------------------------------

     use EDTypesMod, only            : ed_site_type
     use EDTypesMod, only            : ed_patch_type
     use EDSurfaceRadiationMod, only : PatchNormanRadiation
     use FatesInterfaceTypesMod, only     : hlm_numSWb

     ! !ARGUMENTS:
     class(fates_restart_interface_type) , intent(inout) :: this
     integer                     , intent(in)            :: nsites
     type(ed_site_type)          , intent(inout), target :: sites(nsites)
     type(bc_out_type)           , intent(inout)         :: bc_out(nsites)

     ! locals
     ! ----------------------------------------------------------------------------------
     type(ed_patch_type),pointer  :: currentPatch  ! current patch
     integer                      :: s             ! site counter
     integer                      :: ib            ! radiation band counter
     integer                      :: ifp           ! patch counter

     do s = 1, nsites
        
        ifp = 0
        currentpatch => sites(s)%oldest_patch
        do while (associated(currentpatch))  
           ifp = ifp+1
           
           currentPatch%f_sun      (:,:,:) = 0._r8
           currentPatch%fabd_sun_z (:,:,:) = 0._r8
           currentPatch%fabd_sha_z (:,:,:) = 0._r8
           currentPatch%fabi_sun_z (:,:,:) = 0._r8
           currentPatch%fabi_sha_z (:,:,:) = 0._r8
           currentPatch%fabd       (:)     = 0._r8
           currentPatch%fabi       (:)     = 0._r8

           ! zero diagnostic radiation profiles
           currentPatch%nrmlzd_parprof_pft_dir_z(:,:,:,:) = 0._r8
           currentPatch%nrmlzd_parprof_pft_dif_z(:,:,:,:) = 0._r8
           currentPatch%nrmlzd_parprof_dir_z(:,:,:)       = 0._r8
           currentPatch%nrmlzd_parprof_dif_z(:,:,:)       = 0._r8
           
           ! -----------------------------------------------------------
           ! When calling norman radiation from the short-timestep
           ! we are passing in boundary conditions to set the following
           ! variables:
           ! currentPatch%solar_zenith_flag     (is there daylight?)
           ! currentPatch%solar_zenith_angle    (what is the value?)
           ! -----------------------------------------------------------
           
           if(currentPatch%solar_zenith_flag)then
              
              bc_out(s)%albd_parb(ifp,:) = 0._r8  ! output HLM
              bc_out(s)%albi_parb(ifp,:) = 0._r8  ! output HLM
              bc_out(s)%fabi_parb(ifp,:) = 0._r8  ! output HLM
              bc_out(s)%fabd_parb(ifp,:) = 0._r8  ! output HLM
              bc_out(s)%ftdd_parb(ifp,:) = 1._r8  ! output HLM
              bc_out(s)%ftid_parb(ifp,:) = 1._r8  ! output HLM
              bc_out(s)%ftii_parb(ifp,:) = 1._r8  ! output HLM
                 
              if (maxval(currentPatch%nrad(1,:))==0)then
                 !there are no leaf layers in this patch. it is effectively bare ground. 
                 ! no radiation is absorbed  
                 bc_out(s)%fabd_parb(ifp,:) = 0.0_r8
                 bc_out(s)%fabi_parb(ifp,:) = 0.0_r8
                 do ib = 1,hlm_numSWb

                    ! REQUIRES A FIX HERE albd vs albi

                    bc_out(s)%albd_parb(ifp,ib) = currentPatch%gnd_alb_dir(ib)
                    bc_out(s)%albd_parb(ifp,ib) = currentPatch%gnd_alb_dif(ib)
                    bc_out(s)%ftdd_parb(ifp,ib)= 1.0_r8
                    bc_out(s)%ftid_parb(ifp,ib)= 1.0_r8
                    bc_out(s)%ftii_parb(ifp,ib)= 1.0_r8
                 enddo
              else
                 
                 call PatchNormanRadiation (currentPatch, &
                      bc_out(s)%albd_parb(ifp,:), &
                      bc_out(s)%albi_parb(ifp,:), &
                      bc_out(s)%fabd_parb(ifp,:), &
                      bc_out(s)%fabi_parb(ifp,:), &
                      bc_out(s)%ftdd_parb(ifp,:), &
                      bc_out(s)%ftid_parb(ifp,:), &
                      bc_out(s)%ftii_parb(ifp,:))
              
              endif ! is there vegetation? 
              
           end if    ! if the vegetation and zenith filter is active
     

           currentPatch => currentPatch%younger
        end do       ! Loop linked-list patches
     enddo           ! Loop Sites
     
     return
   end subroutine update_3dpatch_radiation


 end module FatesRestartInterfaceMod
