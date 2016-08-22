Module HistoryIOMod

  
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use clm_varctl      , only : iulog

  implicit none

  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.  Definitions are not provide, for an explanation of
  ! the variable go to its registry.  (IH_ signifies "index history")
  
  ! Indices to 1D Patch variables

  integer, private :: ih_trimming_pa
  integer, private :: ih_area_plant_pa
  integer, private :: ih_area_treespread_pa
  integer, private :: ih_canopy_spread_pa
  integer, private :: ih_nesterov_fire_danger_pa
  integer, private :: ih_spitfire_ROS_pa
  integer, private :: ih_effect_wspeed_pa
  integer, private :: ih_TFC_ROS_pa
  integer, private :: ih_fire_intensity_pa
  integer, private :: ih_fire_area_pa
  integer, private :: ih_scorch_height_pa
  integer, private :: ih_fire_fuel_bulkd_pa
  integer, private :: ih_fire_fuel_eff_moist_pa
  integer, private :: ih_fire_fuel_sav_pa
  integer, private :: ih_fire_fuel_mef_pa
  integer, private :: ih_sum_fuel_pa
  integer, private :: ih_litter_in_pa
  integer, private :: ih_litter_out_pa
  integer, private :: ih_efpot_pa        ! NA
  integer, private :: ih_rb_pa           ! NA
  integer, private :: ih_daily_temp
  integer, private :: ih_daily_rh
  integer, private :: ih_daily_prec
  integer, private :: ih_seed_bank
  integer, private :: ih_seeds_in
  integer, private :: ih_seed_decay
  integer, private :: ih_seed_germination
  integer, private :: ih_bstore
  integer, private :: ih_bdead
  integer, private :: ih_balive
  integer, private :: ih_bleaf
  integer, private :: ih_biomass
  integer, private :: ih_npp
  integer, private :: ih_gpp
  integer, private :: ih_autotr_resp
  integer, private :: ih_maint_resp
  integer, private :: ih_growth_resp
  
  ! Indices to (patch x pft) variables   (using nlevgrnd as surrogate)

  integer, private :: ih_biomass_pa_pft
  integer, private :: ih_leafbiomass_pa_pft
  integer, private :: ih_storebiomass_pa_pft
  integer, private :: ih_nindivs_pa_pft

  ! Indices to (site) variables

  integer, private :: ih_nep_si
  integer, private :: ih_nep_timeintegrated_si
  integer, private :: ih_npp_timeintegrated_si
  integer, private :: ih_hr_timeintegrated_si
  integer, private :: ih_nbp_si
  integer, private :: ih_npp_si
  integer, private :: ih_fire_c_to_atm_si
  integer, private :: ih_ed_to_bgc_this_edts_si
  integer, private :: ih_ed_to_bgc_last_edts_si
  integer, private :: ih_seed_rain_flux_si
  integer, private :: ih_totecosysc_si
  integer, private :: ih_totecosysc_old_si
  integer, private :: ih_totedc_si
  integer, private :: ih_totedc_old_si
  integer, private :: ih_totbgcc_si
  integer, private :: ih_totbgcc_old_si
  integer, private :: ih_biomass_stock_si
  integer, private :: ih_ed_litter_stock_si
  integer, private :: ih_cwd_stock_si
  integer, private :: ih_seed_stock_si
  integer, private :: ih_cbalance_error_ed_si
  integer, private :: ih_cbalance_error_bgc_si
  integer, private :: ih_cbalance_error_total_si
  integer, private :: ih_ed_npatches_si
  integer, private :: ih_ed_ncohorts_si
  
  ! Indices to (site x scpf) variables

  integer, private :: ih_gpp_si_scpf
  integer, private :: ih_npp_totl_si_scpf
  integer, private :: ih_npp_leaf_si_scpf
  integer, private :: ih_npp_seed_si_scpf
  integer, private :: ih_npp_fnrt_si_scpf
  integer, private :: ih_npp_bgsw_si_scpf
  integer, private :: ih_npp_bgdw_si_scpf
  integer, private :: ih_npp_agsw_si_scpf
  integer, private :: ih_npp_agdw_si_scpf
  integer, private :: ih_npp_stor_si_scpf
  integer, private :: ih_litt_leaf_si_scpf
  integer, private :: ih_litt_fnrt_si_scpf
  integer, private :: ih_litt_sawd_si_scpf
  integer, private :: ih_litt_ddwd_si_scpf
  integer, private :: ih_r_leaf_si_scpf
  integer, private :: ih_r_stem_si_scpf
  integer, private :: ih_r_root_si_scpf
  integer, private :: ih_r_stor_si_scpf

  integer, private :: ih_ddbh_si_scpf
  integer, private :: ih_ba_si_scpf
  integer, private :: ih_np_si_scpf
  integer, private :: ih_m1_si_scpf
  integer, private :: ih_m2_si_scpf
  integer, private :: ih_m3_si_scpf
  integer, private :: ih_m4_si_scpf
  integer, private :: ih_m5_si_scpf


  ! The number of variable dim/kind types we have defined (static)
  integer, parameter                :: n_iovar_dk = 6


  ! This structure is not allocated by thread, but the upper and lower boundaries
  ! of the dimension for each thread is saved in the clump_ entry
  type iovar_dim_type
     character(len=32) :: name            ! This should match the name of the dimension 
     integer :: lb                       ! lower bound
     integer :: ub                       ! upper bound
     integer,allocatable :: clump_lb(:)  ! lower bound of thread's portion of HIO array
     integer,allocatable :: clump_ub(:)  ! upper bound of thread's portion of HIO array
  end type iovar_dim_type
  

  
  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type iovar_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: patch1_index(:) ! maps site index to the HIO patch 1st position
  end type iovar_map_type


   
  ! This structure is not multi-threaded
  type iovar_dimkind_type
     character(len=32)    :: name        ! String labelling this IO type
     integer              :: ndims       ! number of dimensions in this IO type
     integer, allocatable :: dimsize(:)  ! The size of each dimension
     logical              :: active
     type(iovar_dim_type), pointer :: dim1_ptr
     type(iovar_dim_type), pointer :: dim2_ptr
  end type iovar_dimkind_type


  
  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)
  type iovar_def_type
     character(len=32)    :: vname
     character(len=24)    :: units
     character(len=128)   :: long
     character(len=24)    :: vtype
     character(len=1)     :: avgflag
     real(r8)             :: flushval
     type(iovar_dimkind_type),pointer :: iovar_dk_ptr
     ! Pointers (only one of these is allocated per variable)
     real(r8), pointer     :: r81d(:)
     real(r8), pointer     :: r82d(:,:)
     real(r8), pointer     :: r83d(:,:,:)
     integer,  pointer     :: int1d(:)
     integer,  pointer     :: int2d(:,:)
     integer,  pointer     :: int3d(:,:,:)
  end type iovar_def_type


  type, public :: fates_hio_interface_type
     
     ! Instance of the list of history output varialbes
     type(iovar_def_type), pointer :: hvars(:)
     integer                       :: n_hvars
     
     ! Instanteat one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(iovar_dimkind_type), pointer :: iovar_dk(:)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure
     ! is allocated by number of threads
     type(iovar_dim_type) :: iopa_dim
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure
     ! is allocated by number of threads
     type(iovar_dim_type) :: iosi_dim
     
     ! This is a structure that contains the boundaries for the
     ! ground level (includes rock) dimension
     type(iovar_dim_type) :: iogrnd_dim

     ! This is a structure that contains the boundaries for the
     ! number of size-class x pft dimension
     type(iovar_dim_type) :: ioscpf_dim


     type(iovar_map_type), pointer :: iovar_map(:)
     
   contains
     
     procedure, public :: update_history_variables
     procedure, public :: define_history_vars
     procedure, public :: set_history_var
     procedure, public :: init_iovar_dk_maps
     procedure, public :: iotype_index
     procedure, public :: set_dim_ptrs
     procedure, public :: get_hvar_bounds
     procedure, public :: dim_init
     procedure, public :: set_dim_thread_bounds

  end type fates_hio_interface_type
   


contains
  
  
  ! ====================================================================================
  
  subroutine update_history_variables(this,nc,sites,nsites,fcolumn)
    
    ! ---------------------------------------------------------------------------------
    ! This is the main call to update the history IO arrays that are registerred with
    ! the Host Model.
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type,   &
                                     ed_cohort_type, &
                                     ed_patch_type,  &
                                     AREA
    ! Arguments
    class(fates_hio_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    integer                 , intent(in)            :: nsites
    integer                 , intent(in)            :: fcolumn(nsites)
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_scaling_scalar ! ratio of canopy to patch area for counteracting patch scaling

    type(iovar_def_type),pointer :: hvar
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    
    associate( hio_trimming_pa         => this%hvars(ih_trimming_pa)%r81d, &
               hio_area_plant_pa       => this%hvars(ih_area_plant_pa)%r81d, &
               hio_area_treespread_pa  => this%hvars(ih_area_treespread_pa)%r81d, &
               hio_biomass_pa_pft      => this%hvars(ih_biomass_pa_pft)%r82d, &
               hio_leafbiomass_pa_pft  => this%hvars(ih_leafbiomass_pa_pft)%r82d, &
               hio_storebiomass_pa_pft => this%hvars(ih_storebiomass_pa_pft)%r82d, &
               hio_nindivs_pa_pft      => this%hvars(ih_nindivs_pa_pft)%r82d )

    ! ---------------------------------------------------------------------------------
    ! Flush arrays to values defined by %flushval (see registry entry in
    ! subroutine define_history_vars()
    ! ---------------------------------------------------------------------------------

    do ivar=1,ubound(this%hvars,1)
       hvar => this%hvars(ivar)
       call this%get_hvar_bounds(hvar,nc,lb1,ub1,lb2,ub2)
       select case(trim(hvar%iovar_dk_ptr%name))
       case('PA_R8') 
          hvar%r81d(lb1:ub1) = hvar%flushval
       case('SI_R8') 
          hvar%r81d(lb1:ub1) = hvar%flushval
       case('PA_GRND_R8') 
          hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
       case('PA_SCPF_R8') 
          hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
       case('SI_GRND_R8') 
          hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
       case('SI_SCPF_R8') 
          hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
       case('PA_INT')
          hvar%int1d(lb1:ub1) = nint(hvar%flushval)
       case default
          write(iulog,*) 'iotyp undefined while flushing history variables'
          stop
          !end_run
       end select
    end do
    

    ! Perform flushes or initializations over the FATES-only space?
    ! ---------------------------------------------------------------------------------

    
    ! ---------------------------------------------------------------------------------
    ! Loop through the FATES scale hierarchy and fill the history IO arrays
    ! ---------------------------------------------------------------------------------

    do s = 1,nsites
       
       io_si  = this%iovar_map(nc)%site_index(s)
       io_pa1 = this%iovar_map(nc)%patch1_index(s)
       io_soipa = io_pa1-1


       ! Set trimming on the soil patch to 1.0
       hio_trimming_pa(io_soipa) = 1.0_r8

              
       ipa = 0
       cpatch => sites(s)%oldest_patch
       do while(associated(cpatch))
          
          io_pa = io_pa1 + ipa
          
          ccohort => cpatch%shortest
          do while(associated(ccohort))
             
             ft = ccohort%pft
             
             if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                
                ! for quantities that are at the CLM patch level, because of the way 
                ! that CLM patches are weighted for radiative purposes this # density needs 
                ! to be over either ED patch canopy area or ED patch total area, whichever is less
                n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                
                ! for quantities that are natively at column level, calculate plant 
                ! density using whole area
                n_perm2   = ccohort%n/AREA   
                
             else
                n_density = 0.0_r8
                n_perm2   = 0.0_r8
             endif
             
             if(associated(cpatch%tallest))then
                hio_trimming_pa(io_pa) = cpatch%tallest%canopy_trim
             else
                hio_trimming_pa(io_pa) = 0.0_r8
             endif
             
             hio_area_plant_pa(io_pa) = 1.0_r8
             
             if (min(cpatch%total_canopy_area,cpatch%area)>0.0_r8) then
                hio_area_treespread_pa(io_pa) = cpatch%total_tree_area  &
                      / min(cpatch%total_canopy_area,cpatch%area)
             else
                hio_area_treespread_pa(io_pa) = 0.0_r8
             end if
             
             
             hio_biomass_pa_pft(io_pa,ft) = hio_biomass_pa_pft(io_pa,ft) + &
                   n_density * ccohort%b * 1.e3_r8
             
             hio_leafbiomass_pa_pft(io_pa,ft) = hio_leafbiomass_pa_pft(io_pa,ft) + &
                   n_density * ccohort%bl       * 1.e3_r8
             
             hio_storebiomass_pa_pft(io_pa,ft) = hio_storebiomass_pa_pft(io_pa,ft) + &
                   n_density * ccohort%bstore   * 1.e3_r8
             
             hio_nindivs_pa_pft(io_pa,ft) = hio_nindivs_pa_pft(io_pa,ft) + &
                   ccohort%n
             

             ccohort => ccohort%taller
          enddo ! cohort loop
          
          ipa = ipa + 1
          cpatch => cpatch%younger
       end do !patch loop
       
    enddo ! site loop
    
  end associate

    return
  end subroutine update_history_variables
  
  ! ====================================================================================
  
  subroutine define_history_vars(this,callstep,nvar)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF HISTORY OUTPUT VARIABLES
    !
    ! This subroutine is called in two contexts, either in count mode or inialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
    !
    ! Note 1 there are different ways you can flush or initialize the output fields.
    ! If you flush to a native type, (such as zero), the entire slab which covers
    ! indices which may not be relevant to FATES, are flushed to this value.  So
    ! in that case, lakes and crops that are not controlled by FATES will zero'd
    ! and when values are scaled up to the land-grid, the zero's for non FATES will
    ! be included.  This is good and correct if nothing is there.  
    !
    ! But, what if crops exist in the host model and occupy a fraction of the land-surface
    ! shared with natural vegetation? In that case, you want to flush your arrays
    ! with a value that the HLM treats as "do not average"
    ! 
    ! If your HLM makes use of, and you want, INTEGER OUTPUT, pass the flushval as
    ! a real.  The applied flush value will use the NINT() intrinsic function
    ! ---------------------------------------------------------------------------------
    
    class(fates_hio_interface_type)        :: this
    character(len=*),intent(in)            :: callstep  ! are we 'count'ing or 'initializ'ing?
    integer,optional,intent(out)           :: nvar      
    
    integer                                :: ivar
    
    if(.not. (trim(callstep).eq.'count' .or. trim(callstep).eq.'initialize') ) then
       write(iulog,*) 'defining history variables in FATES requires callstep count or initialize'
       ! end_run('MESSAGE')
    end if
    
    ivar=0
    call this%set_history_var(vname='TRIMMING2',units='none',                   &
         long='Degree to which canopy expansion is limited by leaf economics',  &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=1.0_r8, ivar=ivar,   &
         callstep=callstep,index = ih_trimming_pa)
    
    call this%set_history_var(vname='AREA_PLANT2',units='m2',                   &
         long='area occupied by all plants',                                    &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8,ivar=ivar,    &
         callstep=callstep,index = ih_area_plant_pa)
    
    call this%set_history_var(vname='AREA_TREES2',units='m2',                   &
         long='area occupied by woody plants',                                  &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8,ivar=ivar,    &
         callstep=callstep,index = ih_area_treespread_pa)

    call this%set_history_var(vname='CANOPY_SPREAD2',units='0-1',               &
         long='Scaling factor between tree basal area and canopy area',         &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8,ivar=ivar,    &
         callstep=callstep,index = ih_canopy_spread_pa)

    call this%set_history_var(vname='PFTBIOMASS2',units='gC/m2',                &
         long='total PFT level biomass',                                        &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8,        &
         ivar=ivar,callstep=callstep, index = ih_biomass_pa_pft )

    call this%set_history_var(vname='PFTLEAFBIOMASS2', units='gC/m2',           &
         long='total PFT level leaf biomass',                                   &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8,        &
         ivar=ivar,callstep=callstep, index = ih_leafbiomass_pa_pft )

    call this%set_history_var(vname='PFTSTOREBIOMASS2',  units='gC/m2',         &
         long='total PFT level stored biomass',                                 &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8,        &
         ivar=ivar,callstep=callstep, index = ih_storebiomass_pa_pft )

    call this%set_history_var(vname='PFTNINDIVS2',  units='indiv / m2',         &
         long='total PFT level number of individuals',                          &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8,        &
         ivar=ivar,callstep=callstep, index = ih_nindivs_pa_pft )

    ! Must be last thing before return
    if(present(nvar)) nvar = ivar
    
    return
    
  end subroutine define_history_vars
  
  ! =====================================================================================
   
  subroutine set_history_var(this,vname,units,long,avgflag,vtype,hlms,flushval,ivar,callstep,index)


    use FatesUtilsMod, only : check_hlm_list
    use EDTypesMod, only    : cp_hlm_name

    ! arguments
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)  :: vname
    character(len=*),intent(in)  :: units
    character(len=*),intent(in)  :: long
    character(len=*),intent(in)  :: avgflag
    character(len=*),intent(in)  :: vtype
    character(len=*),intent(in)  :: hlms
    real(r8),intent(in)          :: flushval ! IF THE TYPE IS AN INT WE WILL round with NINT
    character(len=*),intent(in)  :: callstep
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    ! locals
    type(iovar_def_type),pointer :: hvar
    integer :: ub1,lb1,ub2,lb2    ! Bounds for allocating the var
    integer :: ityp
    
    if( check_hlm_list(trim(hlms),trim(cp_hlm_name)) ) then
       
       ivar  = ivar+1
       index = ivar    
       
       if(trim(callstep).eq.'initialize')then
          
          hvar => this%hvars(ivar)
          hvar%vname = vname
          hvar%units = units
          hvar%long  = long
          hvar%vtype = vtype
          hvar%avgflag = avgflag
          hvar%flushval = flushval
          
          ityp=this%iotype_index(trim(vtype))
          hvar%iovar_dk_ptr => this%iovar_dk(ityp)
          this%iovar_dk(ityp)%active = .true.
          
          nullify(hvar%r81d)
          nullify(hvar%r82d)
          nullify(hvar%r83d)
          nullify(hvar%int1d)
          nullify(hvar%int2d)
          nullify(hvar%int3d)
          
          call this%get_hvar_bounds(hvar,0,lb1,ub1,lb2,ub2)
          
          select case(trim(vtype))
          case('PA_R8')
             allocate(hvar%r81d(lb1:ub1))
          case('SI_R8')
             allocate(hvar%r81d(lb1:ub1))
          case('PA_GRND_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2))
          case('PA_SCPF_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2))
          case('SI_GRND_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2))
          case('SI_SCPF_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2))
          case default
             write(iulog,*) 'Incompatible vtype passed to set_history_var'
             write(iulog,*) 'vtype = ',trim(vtype),' ?'
             stop
             ! end_run
          end select
          
       end if
    else
       
       index = 0
    end if
    
    return
  end subroutine set_history_var
  
  ! =====================================================================================

  subroutine get_hvar_bounds(this,hvar,thread,lb1,ub1,lb2,ub2)

     class(fates_hio_interface_type) :: this
     type(iovar_def_type),target,intent(in) :: hvar
     integer,intent(in)              :: thread
     integer,intent(out)             :: lb1
     integer,intent(out)             :: ub1
     integer,intent(out)             :: lb2
     integer,intent(out)             :: ub2

     ! local
     integer :: ndims

     lb1 = 0
     ub1 = 0
     lb2 = 0
     ub2 = 0

     ndims = hvar%iovar_dk_ptr%ndims

     ! The thread = 0 case is the boundaries for the whole proc/node
     if (thread==0) then
        lb1 = hvar%iovar_dk_ptr%dim1_ptr%lb
        ub1 = hvar%iovar_dk_ptr%dim1_ptr%ub
        if(ndims>1)then
           lb2 = hvar%iovar_dk_ptr%dim2_ptr%lb
           ub2 = hvar%iovar_dk_ptr%dim2_ptr%ub
        end if
     else
        lb1 = hvar%iovar_dk_ptr%dim1_ptr%clump_lb(thread)
        ub1 = hvar%iovar_dk_ptr%dim1_ptr%clump_ub(thread)
        if(ndims>1)then
           lb2 = hvar%iovar_dk_ptr%dim2_ptr%clump_lb(thread)
           ub2 = hvar%iovar_dk_ptr%dim2_ptr%clump_ub(thread)
        end if
     end if
     
     return
  end subroutine get_hvar_bounds


  ! ====================================================================================
  
  subroutine init_iovar_dk_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! PA_R8   : 1D patch scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! note (RGK) %active is not used yet. Was intended as a check on the HLM->FATES
    ! control parameter passing to ensure all active dimension types received all
    ! dimensioning specifications from the host, but we currently arent using those
    ! passing functions..
    ! ----------------------------------------------------------------------------------
    
    ! Arguments
    class(fates_hio_interface_type) :: this
       
    ! Locals
    integer            :: ityp
    integer, parameter :: unset_int = -999
    
    allocate(this%iovar_dk(n_iovar_dk))

    ! 1d Patch
    ityp = 1
    this%iovar_dk(ityp)%name  = 'PA_R8'
    this%iovar_dk(ityp)%ndims = 1
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.  
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! 1d Site
    ityp = 2
    this%iovar_dk(ityp)%name  = 'SI_R8'
    this%iovar_dk(ityp)%ndims = 1
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! patch x ground
    ityp = 3
    this%iovar_dk(ityp)%name = 'PA_GRND_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! patch x size-class/pft
    ityp = 4
    this%iovar_dk(ityp)%name = 'PA_SCPF_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! site x ground
    ityp = 5
    this%iovar_dk(ityp)%name = 'SI_GRND_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! site x size-class/pft
    ityp = 6
    this%iovar_dk(ityp)%name = 'SI_SCPF_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)


   
    
    
    
    return
  end subroutine init_iovar_dk_maps
  
  ! ===================================================================================
  
  subroutine set_dim_ptrs(this,dk_name,idim,dim_target)
    
    ! arguments
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)     :: dk_name
    integer,intent(in)              :: idim  ! dimension index
    type(iovar_dim_type),target     :: dim_target
    
    
    ! local
    integer                         :: ityp
    
    ityp = this%iotype_index(trim(dk_name))
    
    ! First check to see if the dimension is allocated
    if(this%iovar_dk(ityp)%ndims<idim)then
       write(iulog,*)'Trying to define dimension size to a dim-type structure'
       write(iulog,*)'but the dimension index does not exist'
       write(iulog,*)'type: ',dk_name,' ndims: ',this%iovar_dk(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if
    
    if(idim==1) then
       this%iovar_dk(ityp)%dim1_ptr => dim_target
    elseif(idim==2) then
       this%iovar_dk(ityp)%dim2_ptr => dim_target
    end if

    ! With the map, we can set the dimension size
    this%iovar_dk(ityp)%dimsize(idim) = dim_target%ub - dim_target%lb + 1

    
    return
 end subroutine set_dim_ptrs
  
  ! ====================================================================================
  
  function iotype_index(this,iotype_name) result(ityp)
    
    ! argument
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)     :: iotype_name

    ! local
    integer :: ityp
    
    do ityp=1,n_iovar_dk
       if(trim(iotype_name).eq.trim(this%iovar_dk(ityp)%name))then
          return
       end if
    end do
    write(iulog,*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    !end_run
    
  end function iotype_index
   
  ! =====================================================================================

  subroutine dim_init(this,iovar_dim,dim_name,nthreads,lb_in,ub_in)

    ! arguments
    class(fates_hio_interface_type) :: this
    type(iovar_dim_type),target     :: iovar_dim
    character(len=*),intent(in)     :: dim_name
    integer,intent(in)              :: nthreads
    integer,intent(in)              :: lb_in
    integer,intent(in)              :: ub_in

    allocate(iovar_dim%clump_lb(nthreads))
    allocate(iovar_dim%clump_ub(nthreads))
    
    iovar_dim%name = trim(dim_name)
    iovar_dim%lb = lb_in
    iovar_dim%ub = ub_in

    return
  end subroutine dim_init

  ! =====================================================================================

  subroutine set_dim_thread_bounds(this,iovar_dim,nc,lb_in,ub_in)

    class(fates_hio_interface_type) :: this
    type(iovar_dim_type),target     :: iovar_dim
    integer,intent(in)              :: nc    ! Thread index
    integer,intent(in)              :: lb_in
    integer,intent(in)              :: ub_in
    
    iovar_dim%clump_lb(nc) = lb_in
    iovar_dim%clump_ub(nc) = ub_in

    return
  end subroutine set_dim_thread_bounds

   ! ====================================================================================
   ! DEPRECATED, TRANSITIONAL OR FUTURE CODE SECTION
   ! ====================================================================================

   !subroutine set_fates_hio_str(tag,iotype_name,iostr_val)

!       ! Arguments
!       character(len=*),intent(in)           :: tag
!       character(len=*), optional,intent(in) :: iotype_name
!       integer, optional, intent(in)         :: iostr_val

!       ! local variables
!       logical              :: all_set
!       integer,  parameter  :: unset_int = -999
!       real(r8), parameter  :: unset_double = -999.9
!       integer              :: ityp, idim

!       select case (trim(tag))
!       case('flush_to_unset')
!          write(*,*) ''
!          write(*,*) 'Flushing FATES IO types prior to transfer from host'
!          do ityp=1,ubound(iovar_str,1)
!             iovar_str(ityp)%dimsize = unset_int
!             iovar_str(ityp)%active  = .false.
!          end do

!       case('check_allset')
!          do ityp=1,ubound(iovar_str,1)
!             write(*,*) 'Checking to see if ',iovar_str(ityp)%name,' IO communicators were sent to FATES'
!             if(iovar_str(ityp)%active)then
!                if(iovar_str(ityp)%offset .eq. unset_int) then
!                   write(*,*) 'FATES offset information of IO type:',iovar_str(ityp)%name
!                   write(*,*) 'was never set'
!                   ! end_run('MESSAGE')
!                end if
!                do idim=1,iovar_str(ityp)%ndims
!                   if(iovar_str(ityp)%dimsize(idim) .eq. unset_int) then
!                      write(*,*) 'FATES dimension information of IO type:',iovar_str(ityp)%name
!                      write(*,*) 'was never set'
!                      ! end_run('MESSAGE')
!                   end if
!                end do
!             end if
!          end do
!          write(*,*) 'Checked. All history IO specifications properly sent to FATES.'
!       case default

!          ! Must have two arguments if this is not a check or flush
!          if(present(iostr_val) .and. present(iotype_name))then
!
!             ! Tag in this case is dimsize or offset
!             select case (trim(tag))
!
!             case('offset')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%offset = iostr_val
!                write(*,*) 'Transfering offset for IOTYPE',iotype_name,' to FATES'

!             case('dimsize1')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%dimsize(1) = iostr_val
!                write(*,*) 'Transfering 1st dimension size for IOTYPE',iotype_name,' to FATES'

!             case('dimsize2')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize,1)==1)then
!                   write(iulog,*) 'Transfering second dimensional bound to unallocated space'
!                   write(iulog,*) 'type:',iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(2) = iostr_val
!                write(*,*) 'Transfering 2nd dimension size for IOTYPE',iotype_name,' to FATES'

!             case('dimsize3')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize,1)<3)then
!                   write(iulog,*) 'Transfering third dimensional bound to unallocated space'
!                   write(iulog,*) 'type:',iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(3) = iostr_val
!                write(*,*) 'Transfering 3rd dimension size for IOTYPE',iotype_name,' to FATES'

!             case default
!                write(*,*) 'IO parameter not recognized:',trim(tag)
!                ! end_run
!             end select
!          else
!             write(*,*) 'no value was provided for the tag'
!          end if
!
!       end select
!       return
!     end subroutine set_fates_hio_str



end module HistoryIOMod
