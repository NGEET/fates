Module HistoryIOMod

  
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use clm_varctl      , only : iulog

  implicit none


  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.
  
  integer, private :: ind_hio_trimming_pa
  integer, private :: ind_hio_area_plant_pa
  integer, private :: ind_hio_area_treespread_pa


 
  integer, parameter                :: n_iovar_dk = 2

  ! This structure is allocated by thread, and there are two instances: patch and site
  type iovar_bounds_type
     integer :: lb1
     integer :: ub1
     integer,allocatable :: clump_lb1(:)  ! lower bound of thread's portion of HIO array
     integer,allocatable :: clump_ub1(:)  ! upper bound of thread's portion of HIO array
  end type iovar_bounds_type


  
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
     type(iovar_bounds_type), pointer :: bounds_ptr
  end type iovar_dimkind_type


  
  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)
  type iovar_def_type
     character(len=32)    :: vname
     character(len=24)    :: units
     character(len=128)   :: long
     character(len=16)    :: vtype
     character(len=1)     :: avgflag
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
     type(iovar_bounds_type) :: iopa_bounds
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure
     ! is allocated by number of threads
     type(iovar_bounds_type) :: iosi_bounds
     
     type(iovar_map_type), pointer :: iovar_map(:)
     
   contains
     
     procedure, public :: update_history_variables
     procedure, public :: define_history_vars
     procedure, public :: set_history_var
     procedure, public :: init_iovar_dk_maps
     procedure, public :: iotype_index
     procedure, public :: set_bounds_map_ptrs
     
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
                                     ed_patch_type
    ! Arguments
    class(fates_hio_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    integer                 , intent(in)            :: nsites
    integer                 , intent(in)            :: fcolumn(nsites)
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_s     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1  ! IO array bounds for the calling thread
    integer  :: ivar     ! index of IO variable object vector
    type(ed_patch_type),pointer :: cpatch
    
    ! ---------------------------------------------------------------------------------
    ! Flush arrays to zero
    ! INTERF-TODO: We need to define a flush type, some variables may not want to
    ! average in zero's for patches that are 
    ! ---------------------------------------------------------------------------------
    do ivar=1,ubound(this%hvars,1)
       
       lb1 = this%hvars(ivar)%iovar_dk_ptr%bounds_ptr%clump_lb1(nc)
       ub1 = this%hvars(ivar)%iovar_dk_ptr%bounds_ptr%clump_ub1(nc)
       
       select case(trim(this%hvars(ivar)%iovar_dk_ptr%name))
       case('PA_R8') 
          this%hvars(ivar)%r81d(lb1:ub1) = 0.0_r8
       case('SI_R8') 
          this%hvars(ivar)%r81d(lb1:ub1) = 0.0_r8
       case default
          write(iulog,*) 'iotyp undefined while flushing history variables'
          stop
          !end_run
       end select
       
    end do
    
    ! Perform any special flushes

    lb1 = this%hvars(ind_hio_trimming_pa)%iovar_dk_ptr%bounds_ptr%clump_lb1(nc)
    ub1 = this%hvars(ind_hio_trimming_pa)%iovar_dk_ptr%bounds_ptr%clump_ub1(nc)
    this%hvars(ind_hio_trimming_pa)%r81d(lb1:ub1) = 1.0_r8
    
    ! ---------------------------------------------------------------------------------
    ! Loop through the FATES scale hierarchy and fill the history IO arrays
    ! ---------------------------------------------------------------------------------
    


    do s = 1,nsites
       
       io_s   = this%iovar_map(nc)%site_index(s)
       io_pa1 = this%iovar_map(nc)%patch1_index(s)
       io_soipa = io_pa1-1

       ! TRIMMING2 (soil patch): ind_hio_trimming_pa
       this%hvars(ind_hio_trimming_pa)%r81d(io_soipa) = 1.0_r8
              
       ipa = 0
       cpatch => sites(s)%oldest_patch
       do while(associated(cpatch))
          
          io_pa = io_pa1 + ipa
          
          ! TRIMMING2: ind_hio_trimming_pa
          if(associated(cpatch%tallest))then
             this%hvars(ind_hio_trimming_pa)%r81d(io_pa) = cpatch%tallest%canopy_trim
          else
             this%hvars(ind_hio_trimming_pa)%r81d(io_pa) = 0.0_r8
          endif
          
          ! AREA_PLANT2: ind_hio_area_plant_pa
          this%hvars(ind_hio_area_plant_pa)%r81d(io_pa) = 1._r8
          
          ! AREA_TREES: ind_hio_area_treespread_pa
          if (min(cpatch%total_canopy_area,cpatch%area)>0.0_r8) then
             this%hvars(ind_hio_area_treespread_pa)%r81d(io_pa) = cpatch%total_tree_area  &
                  / min(cpatch%total_canopy_area,cpatch%area)
          else
             this%hvars(ind_hio_area_treespread_pa)%r81d(io_pa) = 0.0_r8
          end if
          
          
          ipa = ipa + 1
          cpatch => cpatch%younger
       end do !patch loop
       
    enddo ! site loop
    
    return
  end subroutine update_history_variables
  
  ! ====================================================================================
  
  subroutine define_history_vars(this,callstep,nvar)
    
    ! ---------------------------------------------------------------------------------
    ! This subroutine is called in two contexts, either in count mode or inialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
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
    call this%set_history_var(vname='TRIMMING2',units='none', &
         long='Degree to which canopy expansion is limited by leaf economics', &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',ivar=ivar,               &
         callstep=callstep,index = ind_hio_trimming_pa)
    
    
    call this%set_history_var(vname='AREA_PLANT2',units='m2', &
         long='area occupied by all plants', &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',ivar=ivar,               &
         callstep=callstep,index = ind_hio_area_plant_pa)
    
    
    call this%set_history_var(vname='AREA_TREES2',units='m2', &
         long='area occupied by woody plants', &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',ivar=ivar,               &
         callstep=callstep,index = ind_hio_area_treespread_pa)


    ! Must be last thing before return
    if(present(nvar)) nvar = ivar
    
    return
    
  end subroutine define_history_vars
  
  ! =====================================================================================
   
  subroutine set_history_var(this,vname,units,long,avgflag,vtype,hlms,ivar,callstep,index)


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
    character(len=*),intent(in)  :: callstep
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    ! locals
    type(iovar_def_type),pointer :: hvar
    integer :: ub1,lb1,ub2,ub3    ! Bounds for allocating the var
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
          
          ityp=this%iotype_index(trim(vtype))
          hvar%iovar_dk_ptr => this%iovar_dk(ityp)
          this%iovar_dk(ityp)%active = .true.
          
          nullify(hvar%r81d)
          nullify(hvar%r82d)
          nullify(hvar%r83d)
          nullify(hvar%int1d)
          nullify(hvar%int2d)
          nullify(hvar%int3d)
          
          lb1 = hvar%iovar_dk_ptr%bounds_ptr%lb1
          ub1 = hvar%iovar_dk_ptr%bounds_ptr%ub1
          
          select case(trim(vtype))
          case('PA_R8')
             allocate(hvar%r81d(lb1:ub1))
          case('SI_R8')
             allocate(hvar%r81d(lb1:ub1))
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
  
  ! ====================================================================================
  
  subroutine init_iovar_dk_maps(this,nclumps)
    
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
    integer,intent(in)              :: nclumps
       
    ! Locals
    integer            :: ityp
    integer, parameter :: unset_int = -999
    
    allocate(this%iovar_dk(n_iovar_dk))
    print*,"1"

    ityp = 1
    this%iovar_dk(ityp)%name  = 'PA_R8'
    print*,"2"
    this%iovar_dk(ityp)%ndims = 1
    print*,"3"
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    print*,"4"
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    print*,"5"
    this%iovar_dk(ityp)%active = .false.  
    print*,"6"
    nullify(this%iovar_dk(ityp)%bounds_ptr)
    print*,"7"
    
    ityp = 2
    this%iovar_dk(ityp)%name  = 'SI_R8'
    this%iovar_dk(ityp)%ndims = 1
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%bounds_ptr)

    ! Allocate bounds associated with patches
    allocate(this%iopa_bounds%clump_lb1(nclumps))
    allocate(this%iopa_bounds%clump_ub1(nclumps))

    ! Allocate bounds associated with sites
    allocate(this%iosi_bounds%clump_lb1(nclumps))
    allocate(this%iosi_bounds%clump_ub1(nclumps))
    
    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%iovar_map(nclumps))
    
    return
  end subroutine init_iovar_dk_maps
  
  ! ===================================================================================
  
  subroutine set_bounds_map_ptrs(this,iovar_dk_name,map_ptr)
    
    ! arguments
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)     :: iovar_dk_name
    type(iovar_bounds_type),target  :: map_ptr
    
    ! local
    integer                         :: ityp
    
    ityp = this%iotype_index(trim(iovar_dk_name))
    
    this%iovar_dk(ityp)%bounds_ptr => map_ptr
    
    ! With the map, we can set the first dimension size
    this%iovar_dk(ityp)%dimsize(1) = this%iovar_dk(ityp)%bounds_ptr%ub1 - &
         this%iovar_dk(ityp)%bounds_ptr%lb1 + 1

    
    return
  end subroutine set_bounds_map_ptrs
  
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
