module JSONParameterUtilsMod

  ! =============================================================================
  ! This module holds the data types and the routines used
  ! for scanning a JSON parameter file, and storing the contents
  ! Assumptions:
  ! 1) Conformity to the JSON standard, see: https://www.json.org/json-en.html
  ! 2) Two lowest level objects are expected: "dimensions" and "parameters"
  ! 3) Each parameter object inside "parameters" should contain:
  !    a) dtype
  !    b) dims
  !    c) long_name
  !    d) units
  !    e) data
  ! 4) The "dims" object in every parameter must reference objects listed in "dimensions"
  !    unless it is defined as "scalar"
  ! 5) All "data" must be encapsulated in [ ] brackets
  ! 6) All parameters must have a dtype that is either "float" "string" or "integer"
  ! 7) Other lowest level objects may exist (such as history) but will not be read
  ! 8) Things like whitespace, spacing, if data is given its own line, line-breaks...
  !    "shouldn't" matter, so long as the amount of whitespace is reasonable.
  !    (If there are 10k blank characters, we may run out of buffer)
  !
  ! Example JSON file:
  ! 
  !{
  ! "attributes": {"history": "05/11/25, First instantation,
  !                 copied from: ../parameter_files/fates_params_default.cdl."},
  ! "dimensions": {
  !  "fates_NCWD": 4,
  !  "fates_history_age_bins": 7,
  !  "fates_history_coage_bins": 2
  !  "fates_pft": 14
  !  },
  ! "parameters": {
  !  "fates_history_ageclass_bin_edges": {
  !    "dtype": "float",
  !    "dims": ["fates_history_age_bins"],
  !    "long_name": "Lower edges for age class bins used in age-resolved patch history output",
  !    "units": "yr",
  !    "data": [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
  !  }
  ! }
  !}
  !
  ! =========================================================================================
  

  use, intrinsic :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END,IOSTAT_EOR
  use shr_sys_mod   , only: shr_sys_abort

  implicit none
  
  private

  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  
  ! Data types
  integer, parameter :: r_scalar_type = 1
  integer, parameter :: r_1d_type     = 2
  integer, parameter :: r_2d_type     = 3
  integer, parameter :: i_scalar_type = 4
  integer, parameter :: i_1d_type     = 5
  integer, parameter :: i_2d_type     = 6
  integer, parameter :: c_solo_type   = 7
  integer, parameter :: c_1d_type     = 8
  integer, parameter :: c_2d_type     = 9

  logical, parameter :: debug = .false.

  integer, parameter :: max_ll = 2048          ! Maximum allowable line-length
  integer, parameter :: max_sl = 512           ! Maximum allowable symbol-length
  integer, parameter :: max_ul = 128           ! Maximum characters for units
  integer, parameter :: max_scr = 1000         ! Maximum number of strings that
                                               ! need to be remembered at once. This
                                               ! needs to be larger than the
                                               ! total number of vector/array members
                                               ! in any single variable/parameter
  
  integer, parameter :: vardata_max_len = 4096 ! Maximum number of characters that would
                                               ! be needed to encapulate all the data
                                               ! inside a variable's {}

  integer, parameter :: too_many_iter = 1000   ! Fail-safe for loops with logical termination
  integer, parameter :: count_phase = 1        ! We read parameters twice so that we can allocate
  integer, parameter :: fill_phase  = 2        ! and then fill 

  integer            :: log_unit = -1
  
  ! Numeric value to represent invalid values (like NaN and Null)
  ! Overridable, see procedure below
  real(r8) :: r_invalid = -1.e36_r8
  
  type,public :: dim_type
     character(len=max_sl) :: name                      ! Name of the dimension
     integer               :: size                      ! Number of values associated with dimension
  end type dim_type
  
  type,public ::  param_type
     character(len=max_sl) :: name                      ! The variable symbol (the dictionary key)
     character(len=max_ul) :: units                     ! Physical units of parameter
     character(len=max_ll) :: long_name                 ! A long descriptive name of the parameter
     integer               :: dtype                     ! Data type, see list above
                                                        ! combines to explain dimensionality
                                                        ! and type of value, such as int/string/float
     character(len=max_sl), allocatable :: dim_names(:) ! These are the indices of the dimensions
                                                        ! associated with this variable
     integer               :: ndims                     ! Number of dimensions,
                                                        ! same as size (dim_names)
     integer                         :: access_count    ! This is used to count how many
                                                        ! times this parameter has been read
                                                        ! on initialization, expecte number
                                                        ! of times is 1, helps for debugging

     ! --- Data holding structures. Each parameter will only use ONE of these ----
     real(r8)                        :: r_data_scalar
     real(r8), allocatable           :: r_data_1d(:)
     real(r8), allocatable           :: r_data_2d(:,:)
     integer                         :: i_data_scalar
     integer,  allocatable           :: i_data_1d(:)
     integer,  allocatable           :: i_data_2d(:,:)
     character(len=128)              :: c_data
     character(len=128), allocatable :: c_data_1d(:)
     character(len=128), allocatable :: c_data_2d(:,:)
     
  end type param_type

  type,public :: params_type
     type(param_type), pointer :: parameters(:)
     type(dim_type), pointer :: dimensions(:)
   contains
     procedure :: GetDimSizeFromName
     procedure :: GetParamFromName
     procedure :: ReportAccessCounts
  end type params_type
  
  public :: ReadJSON
  public :: SetInvalid
  public :: SetLogInit
  public :: DumpParameter
  
contains

  subroutine SetInvalid(r_invalid_in)
    real(r8),intent(in) :: r_invalid_in
    r_invalid = r_invalid_in
  end subroutine SetInvalid
  
  subroutine SetLogInit(log_unit_in)
    integer,intent(in) :: log_unit_in
    log_unit = log_unit_in
  end subroutine SetLogInit

  ! ==============================================================================
  
  subroutine ReadJSON(filename,file_unit,pstruct)

    ! This is esssentially the driver for reading through the input parameter file
    ! There should be two meaningful sections (not including attributes)
    ! 1) dimensions
    ! 2) parameters
    !
    ! All parameters must be scalar or associated with a known dimension.
    !
    ! Each section is started with an open bracket '{' and closed with a closed bracket '}'
    ! For each section here we end with a call
    ! to rewind and go back to the top of the file
    
    character(len=*),intent(in) :: filename
    integer,intent(in)          :: file_unit
    type(params_type)           :: pstruct
    
    ! Local
    integer :: io_status   ! Read status return value from fortran internal

    ! scratch space for storing data as strings
    ! before its copied into data structures
    character(len=max_ll), dimension(max_scr) :: string_scr 
    
    ! Flush the scratch string
    call ClearStringScratch(string_scr,max_scr)

    ! Open the file, it is assumed that an open file unit has been identified
    ! and passed in as an argument
    open(unit=file_unit, file=filename, status='old', &
        action='READ', iostat=io_status)
    
    if (io_status /= 0) THEN
       write(log_unit,*) 'ERROR: Could not open parameter file: ', trim(filename)
       write(log_unit,*) 'IOSTAT value: ', io_status
       call shr_sys_abort()  ! Terminate program gracefully if file cannot be opened
    else
       if(debug) write(log_unit,*) 'Successfully opened ',trim(filename)
    end if

    call CheckRogueBrackets(file_unit,filename)
    
    call GetDimensions(file_unit,pstruct)

    call GetParameters(file_unit,string_scr,pstruct)
    
    close(unit=file_unit, iostat=io_status)
    if (io_status /= 0) THEN
       write(log_unit,*) 'ERROR: Could not close file: ', TRIM(filename)
       write(log_unit,*) 'IOSTAT value: ', io_status
       call shr_sys_abort()  ! Terminate program gracefully if file cannot be opened
    else
       if(debug) write(log_unit,*) 'sucessfully closed ',trim(filename)
    end if
    
  end subroutine ReadJSON

  ! =====================================================================================

  subroutine CheckRogueBrackets(file_unit,filename)

    ! Go through the file and simply make sure that
    ! { and } brackets are not found inside any quotes.

    integer,intent(in)          :: file_unit
    character(len=*)   :: filename
    character(len=1)   :: filechar
    logical            :: inside_quotes
    integer            :: io_status
    character(len=256) :: io_msg
    
    inside_quotes = .false.
    
    do_readfile: do 
       
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          exit do_readfile
       elseif (io_status < 0) then
          cycle do_readfile
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if
       
       if (scan(filechar,'"')>0)then
          inside_quotes = .not.inside_quotes
       end if

       if (scan('{}',filechar)>0 .and. inside_quotes) then
          write(log_unit,*)'A curly bracket "} or {"  was found inside a quoted string in ',trim(filename)
          write(log_unit,*)'This confuses the JSON parser. Please remove brackets from inside all strings.'
          call shr_sys_abort()
       end if

    end do do_readfile

    if(debug)then
       write(log_unit,*) 'Scanned JSON file and found no rogue brackets'
    end if
    
    return
  end subroutine CheckRogueBrackets
  
  ! =====================================================================================
  
  subroutine GetDimensions(file_unit,pstruct)

    integer,intent(in)          :: file_unit
    type(params_type)           :: pstruct

    integer,parameter :: dimdata_len = 4096 ! We read in ALL the dimension data
                                            ! from { to } into one string, this
                                            ! is how big this string is
    character(len=max_sl) :: group_str      ! String buffer used to identify
                                            ! where the start of the "dimensions" group is
    character(len=dimdata_len) :: dimdata_str ! string that holds all dimension text
    character(len=1) :: filechar
    character(len=max_sl) :: symb_str
    character(len=max_sl) :: data_str
    logical :: found_dimtag   ! Have we found the string "dimensions" yet?
    logical :: found_close    ! Have we found the closing bracket?
    logical :: found_dims     !
    logical :: found_scalar
    integer :: i
    integer :: sep_id,beg_id,end_id
    integer :: io_status
    integer :: n_dims
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str
    character(len=256) :: io_msg

    ! Step 0: We start at the very beginning of the file, which
    ! allows us to count file positions and return to those file positions
    rewind(file_unit)
    
    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         "dimensions:"
    ! -----------------------------------------------------------------------------------
    found_dimtag = .false.
    group_str = ''
    i=0
    do_dimtag: do while(.not.found_dimtag)

       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          write(log_unit,*) 'encountered EOF'
          call shr_sys_abort()
       elseif (io_status < 0) then    ! this is most-likely and end-line character
          cycle do_dimtag
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if
       
       i=i+1
       
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,'"dimensions"')>0 .and. index(group_str,'{',.true.)==min(i,max_sl))then
          found_dimtag = .true.
       end if
    end do do_dimtag
    
    ! -----------------------------------------------------------------------------------
    ! Step 2: Read in all the data until the closing bracket
    ! -----------------------------------------------------------------------------------
    found_close = .false.
    dimdata_str = ''
    i=0
    do_close: do while(.not.found_close)

       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          write(log_unit,*) 'encountered EOF'
          call shr_sys_abort()
       elseif (io_status < 0) then
          cycle do_close
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if

       i=i+1
       
       if(i<=dimdata_len)then
          dimdata_str(i:i) = filechar
       else
          write(log_unit,*) 'Ran out of room reading in dimension data, increase dimdata_len'
          call shr_sys_abort()
       end if
       if(index(dimdata_str,'}')>0)then
          found_close = .true.
       end if
    end do do_close

    ! -----------------------------------------------------------------------------------
    ! Step 3: Parse the dimension data string for dimension data
    !         Each dimension should have a colon ":".
    !         Also, scan for the "scalar" parameter
    !         !!! Just count for now
    ! -----------------------------------------------------------------------------------

    found_scalar = .false.
    if( index(trim(dimdata_str(1:dimdata_len)),'"scalar"')>0 ) then
       found_scalar = .true.
    end if
    
    i = 0
    beg_id = 1
    found_dims = .false.
    do while(.not.found_dims)

       sep_id = index(dimdata_str(beg_id:dimdata_len),':')
       if(sep_id==0)then
          write(log_unit,*)'Expected more dimensions in the parameter file...'
          call shr_sys_abort()
       end if

       i = i + 1
       
       end_id = index(dimdata_str(beg_id:dimdata_len),',')
       if(end_id==0)then
          end_id = index(dimdata_str(beg_id:dimdata_len),'}')
          if(end_id==0)then
             write(log_unit,*)'Trouble parsing dim data, expecting closing bracket'
             call shr_sys_abort()
          end if
          found_dims = .true.
       end if

       beg_id = beg_id + end_id
       
    end do
    
    ! -----------------------------------------------------------------------------------
    ! Step 4: Allocate dimensions and fill the name and sizes
    !         add in the scalar dimension if it was not already defined
    ! -----------------------------------------------------------------------------------
    n_dims = i

    if(found_scalar) then
       allocate(pstruct%dimensions(n_dims))
    else
       allocate(pstruct%dimensions(n_dims+1))
    end if
    
    beg_id = 1
    do i = 1,n_dims
       sep_id = index(dimdata_str(beg_id:dimdata_len),':') + beg_id-1 ! Should be left most...
       end_id = index(dimdata_str(beg_id:dimdata_len),',')
       if(end_id==0)then
          end_id = index(dimdata_str(beg_id:dimdata_len),'}')
       end if
       end_id = end_id + beg_id -1
       symb_str = dimdata_str(beg_id:sep_id-1)
       data_str = dimdata_str(sep_id+1:end_id-1)
       call StringToStringOrReal(data_str,tmp_str,tmp_real,is_num=.true.)
       pstruct%dimensions(i)%name = trim(CleanSymbol(symb_str)) !CleanSymbol left adjusts...
       pstruct%dimensions(i)%size = int(tmp_real)
       if(debug)write(log_unit,*) trim(pstruct%dimensions(i)%name),":",pstruct%dimensions(i)%size
       beg_id = end_id+1
    end do

    if(.not.found_scalar)then
       pstruct%dimensions(n_dims+1)%name = 'scalar'
       pstruct%dimensions(n_dims+1)%size = 1
    end if
    
    
    return
  end subroutine GetDimensions

  ! =====================================================================================

  subroutine ReadCharVar(file_unit,phase,var_num,string_scr,found_vartag,pstruct)

    integer,intent(in)          :: file_unit
    integer,intent(in)          :: phase
    integer,intent(inout)       :: var_num
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    logical,intent(out)         :: found_vartag
    type(params_type),optional   :: pstruct
    type(param_type), pointer :: param
    
    character(len=max_sl) :: group_str
    character(len=vardata_max_len) :: vardata_str
    character(len=1) :: filechar
    character(len=max_ll) :: symb_str
    character(len=vardata_max_len) :: data_str
    integer :: n_vec_out
    logical :: found_close
    integer :: i,j,k
    integer :: sep_id,beg_id,end_id
    integer :: io_status
    integer :: vardata_len
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str
    integer :: dimsizes(2)
    character(len=256) :: io_msg
    
    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         the name of the variable":", and save the symbol name
    ! -----------------------------------------------------------------------------------
    found_vartag = .false.
    group_str = ''
    i=0
    do_vartag: do while(.not.found_vartag)

       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          write(log_unit,*) 'encountered EOF'
          call shr_sys_abort()
       elseif (io_status < 0) then
          cycle do_vartag
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if

       i=i+1
       
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,':')>0 .and. index(group_str,'{',.true.)==min(i,max_sl))then
          found_vartag = .true.
       elseif(index(group_str,'}')>0) then
          ! We found a closing bracket before a variable was defined
          ! this means we are out of parameters
          return
       end if
       if(i>too_many_iter)then
          write(log_unit,*)'failed to find variable string or group closing bracket'
          call shr_sys_abort()
       end if
    end do do_vartag

    sep_id = index(group_str,':')
    symb_str = trim(CleanSymbol(group_str(1:sep_id-1)))
    
    ! ---------------------------------------------------------------------------------
    ! Step 2: Advance through the file until the next closing bracket. Save
    !         EVERYTHING in a large string.
    ! ---------------------------------------------------------------------------------

    found_close = .false.
    vardata_str = ''
    i=0
    do_close: do while(.not.found_close)

       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          write(log_unit,*) 'encountered EOF'
          call shr_sys_abort()
       elseif (io_status < 0) then
          cycle do_close
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if

       i=i+1
       
       if(i<=vardata_max_len)then
          vardata_str(i:i) = filechar
       else
          write(log_unit,*) 'Ran out of room reading in variable data'
          write(log_unit,*) 'Increase the size of integer constant: JSONParameterUtilsMod:vardata_max_len'
          call shr_sys_abort()
       end if
       if(index(vardata_str,'{')>0)then
          write(log_unit,*)'An open bracket was found nested inside the current parameter:',trim(symb_str)
          write(log_unit,*)'This is an invalid file format for the parameter file'
          call shr_sys_abort()
       end if
       if(index(vardata_str,'}')>0)then
          found_close = .true.
       end if
    end do do_close
    vardata_len = i
    
    if(phase==count_phase)then
       var_num = var_num + 1
       ! We've advance the file pointer to the closing bracket
       ! if we are just counting, we are done
       return
    end if

    ! ---------------------------------------------------------------------------------
    ! Step 3: Read through that string. Fill data structures.
    ! ---------------------------------------------------------------------------------

    param => pstruct%parameters(var_num)
    param%name = trim(adjustl(symb_str))
    param%access_count = 0

    if(debug) write(log_unit,*) 'Parameter: ',trim(param%name)
    
    call GetMetaString(vardata_str,'"dims"',beg_id,end_id)
    call GetStringVec(vardata_str(beg_id:end_id),string_scr,n_vec_out)
    allocate(param%dim_names(n_vec_out))
    dimsizes(:)=-1
    param%ndims = n_vec_out
    do i = 1,n_vec_out
       param%dim_names(i) = trim(CleanSymbol(string_scr(i)))
       dimsizes(i) = pstruct%GetDimSizeFromName(trim(param%dim_names(i)))
    end do
    call ClearStringScratch(string_scr,n_vec_out)

    
    call GetMetaString(vardata_str,'"dtype"',beg_id,end_id)
    call GetStringVec(vardata_str(beg_id:end_id),string_scr,n_vec_out)
    data_str = trim(string_scr(1))
    param%dtype = -999
    if(trim(CleanSymbol(data_str))=='float') then
       if(param%ndims>1)then
          param%dtype = r_2d_type
       else
          if(param%dim_names(1)=='scalar')then
             param%dtype = r_scalar_type
          else
             param%dtype = r_1d_type
          end if
       end if
    elseif(trim(CleanSymbol(data_str))=='integer')then
       if(param%ndims>1)then
          param%dtype = i_2d_type
       else
          if(param%dim_names(1)=='scalar')then
             param%dtype = i_scalar_type
          else
             param%dtype = i_1d_type
          end if
       end if
    elseif(trim(CleanSymbol(data_str))=='string')then
       if(param%ndims>1)then
          param%dtype = c_2d_type
       else
          if(param%dim_names(1)=='scalar')then
             param%dtype = c_solo_type
          else
             param%dtype = c_1d_type
          end if
       end if
    else
       write(log_unit,*)'could not properly identify data type'
       write(log_unit,*)trim(data_str)
       write(log_unit,*)param%name
       call shr_sys_abort()
    end if
    call ClearStringScratch(string_scr,n_vec_out)
    
    call GetMetaString(vardata_str,'"long_name"',beg_id,end_id)
    call GetStringVec(vardata_str(beg_id:end_id),string_scr,n_vec_out)
    data_str = trim(string_scr(1))
    param%long_name = trim(CleanSymbol(data_str))
    call ClearStringScratch(string_scr,n_vec_out)
    
    call GetMetaString(vardata_str,'"units"',beg_id,end_id)
    call GetStringVec(vardata_str(beg_id:end_id),string_scr,n_vec_out)
    data_str = trim(string_scr(1))
    param%units = trim(CleanSymbol(data_str))
    call ClearStringScratch(string_scr,n_vec_out)
    
    call GetMetaString(vardata_str,'"data"',beg_id,end_id)
    call GetStringVec(vardata_str(beg_id:end_id),string_scr,n_vec_out)
    select case(param%dtype)
    case(r_scalar_type)
       call StringToStringOrReal(string_scr(1),tmp_str,tmp_real,is_num=.true.)
       param%r_data_scalar = tmp_real
    case(i_scalar_type)
       call StringToStringOrReal(string_scr(1),tmp_str,tmp_real,is_num=.true.)
       param%i_data_scalar = int(tmp_real)
    case(c_solo_type)
       call StringToStringOrReal(string_scr(1),tmp_str,tmp_real,is_num=.false.)
       param%c_data = trim(tmp_str)
    case(r_1d_type,i_1d_type,c_1d_type)
       if(n_vec_out.ne.dimsizes(1))then
          write(log_unit,*)'1d parameter size does not match dimension size'
          write(log_unit,*) n_vec_out,dimsizes(1)
          write(log_unit,*) trim(param%name)
          write(log_unit,*) trim(data_str)
          write(log_unit,*) len(trim(data_str))
          call shr_sys_abort()
       end if
       if(param%dtype==r_1d_type)allocate(param%r_data_1d(dimsizes(1)))
       if(param%dtype==i_1d_type)allocate(param%i_data_1d(dimsizes(1)))
       if(param%dtype==c_1d_type)allocate(param%c_data_1d(dimsizes(1)))
       do i=1,dimsizes(1)
          if( any(param%dtype==(/r_1d_type,i_1d_type/)) )then
             call StringToStringOrReal(string_scr(i),tmp_str,tmp_real,is_num=.true.)
             if(param%dtype==r_1d_type)param%r_data_1d(i) = tmp_real
             if(param%dtype==i_1d_type)param%i_data_1d(i) = int(tmp_real)
          else
             call StringToStringOrReal(string_scr(i),tmp_str,tmp_real,is_num=.false.)
             param%c_data_1d(i) = trim(tmp_str)
          end if
       end do
    case(r_2d_type,i_2d_type,c_2d_type)
       if(n_vec_out.ne.(dimsizes(1)*dimsizes(2)))then
          write(log_unit,*)'2d parameter size does not match dimension size'
          write(log_unit,*) n_vec_out,dimsizes(1),dimsizes(2)
          write(log_unit,*) trim(data_str)
          write(log_unit,*) "len:",len(trim(data_str))
          write(log_unit,*) trim(param%name)
          call shr_sys_abort()
       end if
       if(param%dtype==r_2d_type)allocate(param%r_data_2d(dimsizes(1),dimsizes(2)))
       if(param%dtype==i_2d_type)allocate(param%i_data_2d(dimsizes(1),dimsizes(2)))
       if(param%dtype==c_2d_type)allocate(param%c_data_2d(dimsizes(1),dimsizes(2)))
       i=0
       do j=1,dimsizes(1)
          do k=1,dimsizes(2)
             i=i+1
             if( any(param%dtype==(/r_2d_type,i_2d_type/)) )then
                call StringToStringOrReal(string_scr(i),tmp_str,tmp_real,is_num=.true.)
                if(param%dtype==r_2d_type)param%r_data_2d(j,k) = tmp_real
                if(param%dtype==i_2d_type)param%i_data_2d(j,k) = int(tmp_real)
             else
                call StringToStringOrReal(string_scr(i),tmp_str,tmp_real,is_num=.false.)
                param%c_data_2d(j,k) = trim(tmp_str)
             end if
          end do
       end do
    end select
    call ClearStringScratch(string_scr,n_vec_out)
    
    return
  end subroutine ReadCharVar

  ! =====================================================================================
  
  subroutine GetParameters(file_unit,string_scr,pstruct)

    integer,intent(in)                  :: file_unit
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    type(params_type)                   :: pstruct

    integer,parameter :: data_len = 4096
    character(len=max_sl) :: group_str
    character(len=1) :: filechar
    logical :: found_vartag
    integer :: i,j
    integer :: filepos0
    integer :: io_status
    integer :: n_vars
    character(len=256) :: io_msg
    
    ! Step 0: We start at the very beginning of the file, which
    ! allows us to count file positions and return to those file positions
    rewind(file_unit)

    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         "parameters:"
    ! -----------------------------------------------------------------------------------
    found_vartag = .false.
    group_str = ''
    i=0
    do_vartag: do while(.not.found_vartag)

       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          write(log_unit,*) 'encountered EOF'
          call shr_sys_abort()
       elseif (io_status < 0) then
          cycle do_vartag
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if
       
       i=i+1
       
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,'"parameters"')>0 .and. index(group_str,'{',.true.)==max_sl)then
          found_vartag = .true.
       end if
    end do do_vartag

    ! Remember this exact file position, we will need to return
    filepos0 = i
    
    ! -----------------------------------------------------------------------------------
    ! Step 2: Start a loop through parameters, with each one we are identifying
    ! the { and } brackets and saving that in a string.
    ! -----------------------------------------------------------------------------------
    found_vartag = .true.
    n_vars=0
    do while(found_vartag)
       call ReadCharVar(file_unit,count_phase,n_vars,string_scr,found_vartag)
    end do

    if(debug) write(log_unit,*) 'Found ',n_vars,' in the "parameters" group'
    allocate(pstruct%parameters(n_vars))

    call GotoPos(file_unit,filepos0)
    
    do i=1,n_vars
       j = i
       call ReadCharVar(file_unit,fill_phase,j,string_scr,found_vartag,pstruct)
    end do

    return
  end subroutine GetParameters

  ! =====================================================================================
  
  subroutine GetMetaString(string_in,meta_tag,beg_id,end_id)

    ! This routine searches through a larger string to identify
    ! the substring associated with a metadata tag (such as "dims",
    ! "data", "dtype", etc. It assumes that the separator : follows
    ! the provided string, with no other characters aside from whitespace
    ! (if any) between the last character of the tag and the separator.
    ! It will first identify the index position of the
    ! intput string that is the first value after the colon, and the index
    ! position of the output string that closes it off. This is demarcated
    ! by either a comma or closing bracket "}", and that character
    ! may not be inside square brackets or quotes.

    character(len=*) :: string_in
    character(len=*) :: meta_tag
    integer          :: beg_id
    integer          :: end_id
    integer          :: sep_id
    integer          :: i
    logical          :: open_bracket1
    logical          :: open_bracket2
    logical          :: open_quote
    logical          :: found_end
    logical          :: found_tag
    
    if( index(trim(string_in),trim(meta_tag))==0 ) then
       write(log_unit,*) 'metadata tag: ',trim(meta_tag),' was not found'
       write(log_unit,*) 'in data string: ',trim(string_in),'xxxxx'
       call shr_sys_abort()
    else

       ! We are looking for the meta tag. But it is possilbe
       ! that the meta-tag is also somewhere else in the string.
       ! However we want to identify the metatag that is followed
       ! immediately by the separator and/ whitespaces before the
       ! separator.
       ! So, we identify the first position where we find the tag.
       ! then we search through the next characters, if we find
       ! neither a whitespace or a separator, we know we found
       ! a red-herring, so increement and keep going.
     
       beg_id = 1
       end_id = len(trim(string_in))
       found_tag = .false.
       search_tag: do while(.not.found_tag)
          i = index(string_in(beg_id:end_id),trim(meta_tag))+beg_id-1+len(trim(meta_tag))
          if(debug) write(log_unit,*) trim(meta_tag),'  i:',i,'string:xx',string_in(i:end_id)
          do 
             if(i>=end_id)then
                write(log_unit,*)'could not find metatag and separator :',trim(meta_tag)
                write(log_unit,*)'in data string: ',trim(string_in)
                call shr_sys_abort()
             end if
             if(scan(string_in(i:i),':')>0)then
                found_tag = .true.
                sep_id = i+1
                exit search_tag
             end if
             if(scan(string_in(i:i),' ')==0)then
                beg_id = i
                cycle search_tag
             end if
             i=i+1
          end do
       end do search_tag
       
       open_bracket1 = .false.
       open_bracket2 = .false.
       open_quote   = .false.
       found_end    = .false.  ! either a comma or a }
       i=sep_id
       search_close: do i = sep_id,end_id

          if(.not.open_bracket1)then
             if(scan(string_in(i:i),'[')>0) open_bracket1=.true.
          else
             if(scan(string_in(i:i),'[')>0) open_bracket2=.true.
          end if
          if(open_bracket2)then
             if(scan(string_in(i:i),']')>0) open_bracket2=.false.
          else
             if(scan(string_in(i:i),']')>0) open_bracket1=.false.
          end if
          
          if(scan(string_in(i:i),'"')>0) open_quote=.not.open_quote

          if(.not.open_bracket1 .and. .not.open_quote)then
             if(scan(string_in(i:i),',}')>0)then
                found_end = .true.
                exit search_close
             end if
          end if
       end do search_close

       if(.not.found_end)then
          write(log_unit,*)'could not find an end to a metadata string'
          write(log_unit,*)'base string: xxx',trim(string_in),'xxx'
          call shr_sys_abort()
       end if
       
       beg_id = sep_id
       end_id = i

    end if
    
  end subroutine GetMetaString
  
  ! =====================================================================================

  subroutine GetStringVec(string_in,string_vec_out,n_vec_out)

    ! This routine parses a string, does some cleaning of extra charachters
    ! and generally uses commas as separators to fill in a vector of
    ! cleaned strings. It is assumed that the string passed in contains
    ! the data, not the symbol.
    !
    ! Essentially this is everything in a JSON file after a ":"
    ! If its a vector, it probably looks like this:  [0, 1, 2, 5, 10, 20, 50],
    ! If its a scalar, it probably looks like this:  5,
    ! It should also parse vectors of length 1 just fine...
    ! First step is to look for a "["
    
    character(len=*),intent(in)                  :: string_in
    character(len=max_ll), dimension(*),intent(inout) :: string_vec_out 
    integer,intent(out)                               :: n_vec_out

    !character(len=2)                                :: unwanted_chars  = (/!#/)
    integer :: i,j,k,l
    integer :: iv    ! index of vector start
    integer :: ivz   ! index of vector end
    logical :: has_brackets
    
    ! Determine if this is a vector or not
    ! (ie iv=0 is scalar, iv>0 is vector)
    ! Take the right most (in case nested brackets)
    
    iv = index(trim(string_in),'[')
    if(iv>0) then
       has_brackets = .true.
       ivz = index(trim(string_in),']',.true.) ! pos of last closing bracket
       if(ivz==0)then
          write(log_unit,*)'no closing bracket when one was opened?'
          call shr_sys_abort()
       end if
    else
       iv = 1
       has_brackets = .false.
       ivz = len(trim(string_in))
    end if

       
    k = 0 ! Counter for the cleaned string position
    l = 1 ! Counter for the output vector
    
    ! Loop through every character of the original string

    do_str_scan: do i = iv, ivz

       ! We are done parsing if we hit one of two conditions..
       ! 1) We are in vector mode and we encounter the end bracket
       ! 2) We are not in vector mode (no opening bracket) and we encounter a comma
       if(has_brackets .and. i==ivz) exit do_str_scan
       if(.not.has_brackets .and. index(string_in(i:i),',')>0) exit do_str_scan

       ! If we encounter a comma, we go to the next index
       ! and cycle
       if(index(string_in(i:i),",")>0) then
          if(k==0)then
             write(log_unit,*) 'Encountered comma before data in string read, aborting'
             call shr_sys_abort()
          end if
          l = l + 1
          k = 0
          cycle do_str_scan
       else
          
          ! SCAN returns the position of the *first* character from the 
          ! set (unwanted_chars) that is present in the string.
          ! If it returns 0, the character is NOT unwanted.
          !j = scan(string_in(i:i),'!#&^[]')
          j = scan(string_in(i:i),'[]')
          
          if (j == 0) then
             ! This is a character we want to keep
             k = k + 1
             ! Append the character to the cleaned string
             string_vec_out(l)(k:k) = string_in(i:i)
          end if
          
       end if
       
    end do do_str_scan

    n_vec_out = l

  end subroutine GetStringVec

  ! =====================================================================================
  
  subroutine PopString(apstring,newchar)
    
    character(len=*) :: apstring
    character(len=1) :: newchar
    integer :: strlen
    integer :: i
    
    strlen = len(apstring)
    do i = 2,strlen
       apstring(i-1:i-1) = apstring(i:i)
    end do
    apstring(strlen:strlen) = newchar
    
    return
  end subroutine PopString

  ! =====================================================================================
  
  function CleanSymbol(string_in) result(string_out)

    character(len=*) string_in
    character(len=max_sl) :: string_out
    integer :: i,j,k,jj
    
    ! Flush the string
    string_out = ''
    string_in  = adjustl(string_in)
    
    k=0
    do i = 1, len(trim(string_in))

       ! SCAN returns the position of the *first* character from the 
       ! set (unwanted_chars) that is present in the string.
       ! If it returns 0, the character is NOT unwanted.
       !j = scan(string_in(i:i),'!#&^"@/><}{')
       j = scan(string_in(i:i),"'")
       jj = scan(string_in(i:i),'"')
       if (j == 0 .and. jj==0 ) then
          ! This is a character we want to keep
          k = k + 1
          ! Append the character to the cleaned string
          string_out(k:k) = string_in(i:i)
       end if

       ! remove trailing commas
       j = scan(string_out,',',.true.)  ! right most occurance...
       if(j==len(trim(string_out)) .and. j>0)then
          string_out(j:j)=''
       end if
       
    end do
    
  end function CleanSymbol

  ! =====================================================================================
  
  subroutine ClearStringScratch(string_scr,n_vec_in)
    character(len=max_ll), dimension(*),intent(inout) :: string_scr
    integer,intent(in) :: n_vec_in
    integer :: i
    do i = 1,n_vec_in
       string_scr(i) = ''
    end do
    return
  end subroutine ClearStringScratch

  ! =====================================================================================
  
  subroutine GotoPos(file_unit,filepos)
    integer            :: file_unit
    integer,intent(in) :: filepos

    character(len=1) dummychar
    integer :: i         ! file position number index
    integer :: io_status ! i/o status of read
    character(len=256) :: io_msg
    
    ! note that we do not include EOLs as valid file positions...
    
    rewind(file_unit)

    if (filepos > 1) then

       i=0
       do_filepos: do while(i <= filepos)
       
          ! empty read
          read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) dummychar
          if (io_status == IOSTAT_END) then
             write(log_unit,*) 'encountered EOF'
             call shr_sys_abort()
          elseif (io_status < 0) then
             cycle do_filepos
          elseif (io_status>0) then
             write(log_unit,*) 'fatal i/o error! code:', io_status
             write(log_unit,*) 'msg: ',io_msg
             call shr_sys_abort()
          end if
          
          i=i+1
          
       end do do_filepos
    end if
    
  end subroutine GotoPos

  ! =====================================================================================
  
  subroutine StringToStringOrReal(string_in,out_string,out_real,is_num)

    ! This routine takes a string, and determines if it is a number
    ! by trying to read the string into a floating point number
    ! It uses the iostat report to determine if it was possible or not
    ! If the string was not a number, it returns a cleaned version of the
    ! string. It was a number, it returns that as well.
    
    character(len=*),intent(in) :: string_in
    logical,intent(in)          :: is_num  ! This specifies if we think the value is a number
    character(len=max_sl),intent(out) :: out_string
    real(r8),intent(out)              :: out_real

    character(len=max_sl)       :: tmp_str
    character(len=max_sl)       :: tmp_lc_str
    integer  :: io_status
    real(r8) :: tmp_real
    logical  :: is_invalid
    
    tmp_str = trim(CleanSymbol(string_in))

    tmp_lc_str = tmp_str

    call to_lowercase(tmp_lc_str)
    
    ! Check if this is a nan
    is_invalid = index(tmp_lc_str,'nan')>0 .or. index(tmp_lc_str,'null')>0

    if(is_invalid)then
       out_real = r_invalid
       out_string = 'null'
    else
       read (unit=tmp_str,fmt=*,iostat=io_status) tmp_real
       if (io_status == 0) then
          if(is_num)then
             out_real = tmp_real
             out_string = 'real number'
          else
             write(log_unit,*) 'Read string, expected number',trim(CleanSymbol(string_in))
             call shr_sys_abort()
          end if
       elseif(io_status /= 0) then
          if(.not.is_num)then
             out_real = -9999.9_r8
             out_string = trim(CleanSymbol(string_in))
          else
             write(log_unit,*) 'Read string, expected character',trim(CleanSymbol(string_in))
             call shr_sys_abort()
          end if
       end if
    end if
    return
  end subroutine StringToStringOrReal

  ! =====================================================================================

  subroutine to_lowercase(input_string)
    character(len=*), intent(inout) :: input_string
    integer :: i, char_code

    ! Use IACHAR to convert the character to its integer
    ! code. Capital letters have a different code than lower
    ! case.
    
    do i = 1,len(input_string)

       char_code = iachar(input_string(i:i))
        
       ! Check if the character is an uppercase letter (A-Z)
       if (char_code >= iachar('A') .and. char_code <= iachar('Z')) then
          ! Add 32 (the ASCII offset) to convert it to lowercase
          input_string(i:i) = ACHAR(char_code + 32)
       end if
    end do
  end subroutine to_lowercase

  ! =====================================================================================

  function GetDimSizeFromName(this,dim_name) result(dim_size)

    class(params_type) :: this
    character(len=*) :: dim_name
    integer          :: dim_size
    integer          :: i
    
    dim_size = -1
    loop_dims: do i = 1,size(this%dimensions)
       if(index(this%dimensions(i)%name,trim(dim_name))>0)then
          dim_size = this%dimensions(i)%size
          exit loop_dims
       end if
    end do loop_dims

    if(dim_size==-1)then
       write(log_unit,*)'could not find a size for unknown dimension: ',trim(dim_name)
       do i = 1,size(this%dimensions)
          write(log_unit,*) trim(this%dimensions(i)%name),this%dimensions(i)%size
       end do
       call shr_sys_abort()
    end if
  end function GetDimSizeFromName

  ! =====================================================================================

  function GetParamFromName(this,param_name) result(param_ptr)

    class(params_type)       :: this
    character(len=*)         :: param_name
    type(param_type),pointer :: param_ptr
    integer                  :: i
    
    nullify(param_ptr)
    loop_params: do i = 1,size(this%parameters)
       if(trim(param_name)==this%parameters(i)%name)then
          param_ptr=>this%parameters(i)
          this%parameters(i)%access_count = this%parameters(i)%access_count + 1
          return
       end if
    end do loop_params

    write(log_unit,*)'Error finding parameter by name,scanned ',size(this%parameters),' parameters'
    write(log_unit,*)'Cant find: ',trim(param_name)
    call shr_sys_abort()
    
  end function GetParamFromName

  ! =====================================================================================

  subroutine ReportAccessCounts(this)

    class(params_type)       :: this
    integer                  :: i
    
    write(log_unit,*) 'Reporting parameter access counts'
    write(log_unit,*) ''
    write(log_unit,*) 'Parameters accessed more than once:'
    loop_params2: do i = 1,size(this%parameters)
       if(this%parameters(i)%access_count>1)then
          write(log_unit,*) trim(this%parameters(i)%name),', count: ',this%parameters(i)%access_count
       end if
    end do loop_params2
    write(log_unit,*) ''
    write(log_unit,*) 'Parameters that were not accessed:'
    loop_params0: do i = 1,size(this%parameters)
       if(this%parameters(i)%access_count==0)then
          write(log_unit,*) trim(this%parameters(i)%name),', count: ',this%parameters(i)%access_count
       end if
    end do loop_params0
    write(log_unit,*) 'All other parameters were accessed once (expected value)'
       

    return
  end subroutine ReportAccessCounts
  
  ! =====================================================================================
  
  subroutine DumpParameter(param)
    
    type(param_type) :: param
    integer :: i
    integer :: j
    
    write(log_unit,*) '----------------------------'
    write(log_unit,*) 'Parameter: ',trim(param%name)
    write(log_unit,*) '  units: ',trim(param%units)
    write(log_unit,*) '  long_name: ',trim(param%long_name)
    write(log_unit,*) '  dtype: ',param%dtype
    do i=1,param%ndims
       write(log_unit,*) '  dimension: ',i,':',trim(param%dim_names(i))
    end do
    write(log_unit,*) '  data: '
    
    
    select case(param%dtype)
    case(r_scalar_type)
       write(log_unit,'(F15.6)')param%r_data_scalar
    case(i_scalar_type)
       write(log_unit,*)param%i_data_scalar
    case(c_solo_type)
       write(log_unit,*)trim(param%c_data)
    case(r_1d_type)
       do i = 1,size(param%r_data_1d,dim=1)
          write(log_unit,'(F15.6)') param%r_data_1d(i)
       end do
    case(i_1d_type)
       do i = 1,size(param%i_data_1d,dim=1)
          write(log_unit,*) param%i_data_1d(i)
       end do
    case(c_1d_type)
       do i = 1,size(param%c_data_1d,dim=1)
          write(log_unit,*) trim(param%c_data_1d(i))
       end do
    case(r_2d_type)
       do i = 1,size(param%r_data_2d,dim=1)
          do j = 1,size(param%r_data_2d,dim=2)
             write(log_unit,'(F15.6)') param%r_data_2d(i,j)
          end do
       end do
    case(i_2d_type)
       do i = 1,size(param%i_data_2d,dim=1)
          do j = 1,size(param%i_data_2d,dim=2)
             write(log_unit,*) param%i_data_2d(i,j)
          end do
       end do
    case(c_2d_type)
       do i = 1,size(param%c_data_2d,dim=1)
          do j = 1,size(param%c_data_2d,dim=2)
             write(log_unit,*) trim(param%c_data_2d(i,j))
          end do
       end do
    end select
    
  end subroutine DumpParameter

  
end module JSONParameterUtilsMod
