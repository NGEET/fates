module FatesHistoryVariableType

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only : fates_log
  use FatesHistoryVariableKindMod, only : fates_history_variable_kind_type

  implicit none

  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)

  type fates_history_variable_type
     character(len=32)    :: vname
     character(len=24)    :: units
     character(len=128)   :: long
     character(len=24)    :: use_default ! States whether a variable should be turned
                                         ! on the output files by default (active/inactive)
                                         ! It is a good idea to set inactive for very large
                                         ! or infrequently used output datasets
     character(len=24)    :: vtype
     character(len=1)     :: avgflag
     integer              :: upfreq  ! Update frequency (this is for checks and flushing)
                                     ! 1 = dynamics "dyn" (daily)
                                     ! 2 = production "prod" (prob model tstep)
     real(r8)             :: flushval
     type(fates_history_variable_kind_type),pointer :: iovar_dk_ptr
     ! Pointers (only one of these is allocated per variable)
     real(r8), pointer     :: r81d(:)
     real(r8), pointer     :: r82d(:,:)
     real(r8), pointer     :: r83d(:,:,:)
     integer,  pointer     :: int1d(:)
     integer,  pointer     :: int2d(:,:)
     integer,  pointer     :: int3d(:,:,:)
   contains
     procedure, public :: Init => InitHistoryVariableType
     procedure, public :: Flush => FlushVar
     procedure, private :: GetBounds
  end type fates_history_variable_type

contains

  subroutine InitHistoryVariableType(this, vname, units, long, use_default, &
       vtype, avgflag, flushval, upfreq, n_iovar_dk, iovar_dk)

    use FatesHistoryDimensionMod, only : patch_r8, patch_ground_r8, patch_class_pft_r8, &
         site_r8, site_ground_r8, site_class_pft_r8

    implicit none

    class(fates_history_variable_type), intent(inout) :: this
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: long
    character(len=*), intent(in) :: use_default
    character(len=*), intent(in) :: vtype
    character(len=*), intent(in) :: avgflag
    real(r8), intent(in) :: flushval ! If the type is an int we will round with nint
    integer, intent(in) :: upfreq
    integer, intent(in) :: n_iovar_dk
    type(fates_history_variable_kind_type), intent(in), target :: iovar_dk(:)    

    integer :: ityp
    integer :: lb1, ub1, lb2, ub2
    
    this%vname = vname
    this%units = units
    this%long  = long
    this%use_default = use_default
    this%vtype = vtype
    this%avgflag = avgflag
    this%flushval = flushval
    this%upfreq = upfreq

    nullify(this%r81d)
    nullify(this%r82d)
    nullify(this%r83d)
    nullify(this%int1d)
    nullify(this%int2d)
    nullify(this%int3d)

    ityp = iotype_index(trim(vtype), n_iovar_dk, iovar_dk)
    this%iovar_dk_ptr => iovar_dk(ityp)
    this%iovar_dk_ptr%active = .true.
                
    call this%GetBounds(0, lb1, ub1, lb2, ub2)
          
    ! NOTE(rgk, 2016-09) currently, all array spaces are flushed each
    ! time the update is called. The flush here on the initialization
    ! may be redundant, but will prevent issues in the future if we
    ! have host models where not all threads are updating the HHistory
    ! array spaces.

    select case(trim(vtype))
    case(patch_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(site_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(patch_ground_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(patch_class_pft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_ground_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_class_pft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case default
       write(fates_log(),*) 'Incompatible vtype passed to set_history_var'
       write(fates_log(),*) 'vtype = ',trim(vtype),' ?'
       stop
       ! end_run
    end select
    
  end subroutine InitHistoryVariableType
  
  ! =====================================================================================

  subroutine GetBounds(this, thread, lb1, ub1, lb2, ub2)

     class(fates_history_variable_type), intent(inout) :: this

     integer, intent(in)  :: thread

     integer, intent(out) :: lb1
     integer, intent(out) :: ub1
     integer, intent(out) :: lb2
     integer, intent(out) :: ub2

     ! local
     integer :: ndims

     lb1 = 0
     ub1 = 0
     lb2 = 0
     ub2 = 0

     ndims = this%iovar_dk_ptr%ndims

     ! The thread = 0 case is the boundaries for the whole proc/node
     if (thread==0) then
        lb1 = this%iovar_dk_ptr%dim1_ptr%lower_bound
        ub1 = this%iovar_dk_ptr%dim1_ptr%upper_bound
        if(ndims>1)then
           lb2 = this%iovar_dk_ptr%dim2_ptr%lower_bound
           ub2 = this%iovar_dk_ptr%dim2_ptr%upper_bound
        end if
     else
        lb1 = this%iovar_dk_ptr%dim1_ptr%clump_lower_bound(thread)
        ub1 = this%iovar_dk_ptr%dim1_ptr%clump_upper_bound(thread)
        if(ndims>1)then
           lb2 = this%iovar_dk_ptr%dim2_ptr%clump_lower_bound(thread)
           ub2 = this%iovar_dk_ptr%dim2_ptr%clump_upper_bound(thread)
        end if
     end if
     
   end subroutine GetBounds

  subroutine FlushVar(this, thread)

    use FatesHistoryDimensionMod, only : patch_r8, patch_ground_r8, patch_class_pft_r8, &
         site_r8, site_ground_r8, site_class_pft_r8, patch_int

    implicit none

    class(fates_history_variable_type), intent(inout) :: this
    integer, intent(in) :: thread

    integer :: lb1, ub1, lb2, ub2
    
    call this%GetBounds(thread, lb1, ub1, lb2, ub2)

    select case(trim(this%iovar_dk_ptr%name))
    case(patch_r8) 
       this%r81d(lb1:ub1) = this%flushval
    case(site_r8) 
       this%r81d(lb1:ub1) = this%flushval
    case(patch_ground_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(patch_class_pft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_ground_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_class_pft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(patch_int)
       this%int1d(lb1:ub1) = nint(this%flushval)
    case default
       write(fates_log(),*) 'fates history variable type undefined while flushing history variables'
       stop
       !end_run
    end select
    
  end subroutine FlushVar

  ! ====================================================================================
  
  function iotype_index(iotype_name, n_iovar_dk, iovar_dk) result(ityp)
    
    ! argument
    character(len=*), intent(in) :: iotype_name
    integer, intent(in) :: n_iovar_dk
    type(fates_history_variable_kind_type), intent(in) :: iovar_dk(:)

    ! local
    integer :: ityp
    
    do ityp=1, n_iovar_dk
       if(trim(iotype_name).eq.trim(iovar_dk(ityp)%name))then
          return
       end if
    end do
    write(fates_log(),*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    !end_run
    
  end function iotype_index
   

end module FatesHistoryVariableType
