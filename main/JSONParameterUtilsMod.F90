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

  integer, parameter :: too_many_iter = 10000  ! Fail-safe for loops with logical termination

  integer            :: log_unit = -1

  integer, parameter :: ascii_wspace_set(5) = (/9,10,12,13,32/)

  integer, parameter :: scalar_obj = 1
  integer, parameter :: array_obj  = 2
  integer, parameter :: nested_obj = 3
  
  ! Numeric value to represent invalid values (like NaN and Null)
  ! Overridable, see procedure below
  real(r8) :: r_invalid = -1.e36_r8

  ! This large string holds the file contents, minus
  ! all non-standard characters like end-lines and such
  ! This will be deallocated at the end of the parsing
  ! process
  character(len=:), allocatable  :: file_buffer
  integer                        :: nbuffer
  character(len=vardata_max_len) :: obj_buffer

  
  character(len=max_sl),parameter :: search_next_tag = 'search-next-tag'
  
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
     integer               :: access_count              ! This is used to count how many
                                                        ! times this parameter has been read
                                                        ! on initialization, expecte number
                                                        ! of times is 1, helps for debugging

     ! --- Data holding structures. Each parameter will only use ONE of these ----
     real(r8)                           :: r_data_scalar
     real(r8), allocatable              :: r_data_1d(:)
     real(r8), allocatable              :: r_data_2d(:,:)
     integer                            :: i_data_scalar
     integer,  allocatable              :: i_data_1d(:)
     integer,  allocatable              :: i_data_2d(:,:)
     character(len=max_sl)              :: c_data
     character(len=max_sl), allocatable :: c_data_1d(:)
     character(len=max_sl), allocatable :: c_data_2d(:,:)
     
  end type param_type

  type,public :: params_type
     type(param_type), pointer :: parameters(:)
     type(dim_type), pointer :: dimensions(:)
   contains
     procedure :: GetDimSizeFromName
     procedure :: GetParamFromName
     procedure :: ReportAccessCounts
     procedure :: Destroy
  end type params_type
  
  public :: JSONRead
  public :: JSONSetInvalid
  public :: JSONSetLogInit
  public :: JSONDumpParameter
  
contains

  subroutine JSONSetInvalid(r_invalid_in)
    real(r8),intent(in) :: r_invalid_in
    r_invalid = r_invalid_in
  end subroutine JSONSetInvalid
  
  subroutine JSONSetLogInit(log_unit_in)
    integer,intent(in) :: log_unit_in
    log_unit = log_unit_in
  end subroutine JSONSetLogInit

  ! ==============================================================================
  
  subroutine JSONRead(filename,pstruct)

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
    type(params_type)           :: pstruct
    
    ! Local
    integer :: io_status   ! Read status return value from fortran internal
    integer :: file_unit

    ! scratch space for storing data as strings
    ! before its copied into data structures
    character(len=max_ll), dimension(max_scr) :: string_scr 
    
    ! Flush the scratch string
    call ClearStringScratch(string_scr,max_scr)

    ! Open the file, it is assumed that an open file unit has been identified
    ! and passed in as an argument
    open(newunit=file_unit, file=filename, status='old', &
        action='READ', iostat=io_status)
    
    if (io_status /= 0) THEN
       write(log_unit,*) 'ERROR: Could not open parameter file: ', trim(filename)
       write(log_unit,*) 'IOSTAT value: ', io_status
       call shr_sys_abort()  ! Terminate program gracefully if file cannot be opened
    else
       if(debug) write(log_unit,*) 'Successfully opened ',trim(filename)
    end if

    ! Transfer text file into one large string buffer
    call AllocFillBuffer(file_unit)

    close(unit=file_unit, iostat=io_status)
    if (io_status /= 0) then
       write(log_unit,*) 'ERROR: Could not close file: ', trim(filename)
       write(log_unit,*) 'IOSTAT value: ', io_status
       call shr_sys_abort()  ! Terminate program gracefully if file cannot be opened
    else
       if(debug) write(log_unit,*) 'sucessfully closed ',trim(filename)
    end if
    
    call GetDimensions(pstruct)

    call GetParameters(string_scr,pstruct)
    
    return
  end subroutine JSONRead

  ! =====================================================================================

  subroutine AllocFillBuffer(file_unit)

    integer,intent(in) :: file_unit
    character(len=1)   :: filechar
    integer            :: io_status
    character(len=256) :: io_msg
    integer            :: i
    logical            :: inside_dq ! Inside a double quote?
    
    rewind(file_unit)
    inside_dq = .false.
    i = 0
    countchar: do 
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          exit countchar
       elseif (io_status < 0) then
          cycle countchar
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if

       ! Omit white-space that is not inside a double quote
       if(scan(filechar,'"')>0)then
          inside_dq = .not.inside_dq
       end if
       if(.not.(.not.inside_dq .and. any(ascii_wspace_set == iachar(filechar)))) then
          i = i + 1
       end if
    end do countchar

    nbuffer = i
    allocate(character(len=nbuffer) :: file_buffer)

    rewind(file_unit)
    
    i = 0
    fillchar: do 
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status, iomsg=io_msg) filechar
       if (io_status == IOSTAT_END) then
          exit fillchar
       elseif (io_status < 0) then
          cycle fillchar
       elseif (io_status>0) then
          write(log_unit,*) 'fatal i/o error! code:', io_status
          write(log_unit,*) 'msg: ',io_msg
          call shr_sys_abort()
       end if
       ! Omit white-space that is not inside a double quote
       if(scan(filechar,'"')>0)then
          inside_dq = .not.inside_dq
       end if
       if(.not.(.not.inside_dq .and. any(ascii_wspace_set == iachar(filechar)))) then
          i = i + 1
          file_buffer(i:i) = filechar
       end if
       
    end do fillchar
    
  end subroutine AllocFillBuffer
  
  ! ======================================================================================
  
  subroutine GetNextObject(beg_pos,obj_type,next_pos,is_terminal)
    
    ! --------------------------------------------------------------------
    ! Provide a position index in the file string buffer
    ! and this will transfer all of the text
    ! associated with the object into the obj_buffer.
    ! If it is complex object with curly brackets, it will return
    ! everything between them.
    ! If it is a single double quoted string, it will
    ! return the contents of the string with the double quotes
    ! If it is a single value, it will return that value as a string
    ! If it is an array (or nested) it will have the array brackets
    ! included.
    ! It will also return the index position FOLLOWING the terminator
    ! --------------------------------------------------------------------

    integer :: beg_pos     ! input, where to start
    integer :: next_pos    ! Index position following this object's terminator
                           ! If this is a nested object, position after the "}"
                           ! If this is the an array or single value, it
                           ! follows the trailing comma. If this is the
                           ! last single value in the nest, its the
                           ! position of the parent's "}"
    integer :: obj_type    ! This is either a "scalar","array" or "nested"
    logical :: is_terminal ! Is this the last object in it's nest?
                           ! ie is the character after this a "}"
    ! locals
    logical :: open_cb     ! Are we inside a curly bracket?
    logical :: open_dq     ! Are we inside double quotes?
    logical :: open_sb1    ! Are we inside outer square brackets?
    logical :: open_sb2    ! Are we insdie inner square brackets?
    integer :: i           ! object buffer position
    integer :: fpos        ! file string buffer position
    character(len=1) :: filechar ! A single character
    
    ! Clear out the string buffer
    obj_buffer = repeat(' ',vardata_max_len)
    
    open_cb  = .false.
    open_dq  = .false.
    open_sb1 = .false.
    open_sb2 = .false.
    is_terminal = .false.

    ! We assume scalar unless we encounter brackets
    obj_type = scalar_obj
    
    i=1
    fpos = beg_pos
    search_obj: do

       filechar = file_buffer(fpos:fpos)
       
       if(scan(filechar,'"')>0)then
          open_dq = .not.open_dq
       end if

       if_doublequote: if(open_dq) then

          obj_buffer(i:i) = filechar
          i = i + 1
          fpos = fpos + 1
          
       else
          
          if (scan(file_buffer(fpos:fpos),'{')>0) then
             obj_type = nested_obj
             open_cb = .true.
             fpos = fpos + 1
             cycle search_obj
          end if
       
          if_curlybrack: if (open_cb)then
          
             if(scan(filechar,'}')>0)then
                fpos = fpos + 1
                if (scan(file_buffer(fpos:fpos),'}')>0) then
                   is_terminal = .true.
                end if
                exit search_obj
             else
                obj_buffer(i:i) = filechar
                i = i + 1
                fpos = fpos + 1
             end if
             
          else
             
             obj_buffer(i:i) = filechar
             i = i + 1
             fpos = fpos + 1
             
             ! Array stuff, if the hard brackets
             ! are open, then we cannot terminate yet
             if_squarebrack: if(open_sb1) then
                if(open_sb2) then
                   if (scan(filechar,']')>0) then
                      open_sb2 = .false.
                   end if
                else
                   if (scan(filechar,'[')>0) then
                      open_sb2 = .true.
                   end if
                   if (scan(filechar,']')>0) then
                      if (scan(file_buffer(fpos:fpos),'}')>0) then
                         is_terminal = .true.
                      end if
                      exit search_obj
                   end if
                end if
             else
                if (scan(filechar,'[')>0) then
                   obj_type = array_obj
                   open_sb1 = .true.
                end if
                ! Its possible that this object is the last
                ! nested object, so it will be either terminated
                ! by a comma or a curly bracket from it's parent
                if (scan(file_buffer(fpos:fpos),',')>0) then
                   fpos = fpos + 1
                   exit search_obj
                end if
                if (scan(file_buffer(fpos:fpos),'}')>0) then
                   fpos = fpos + 1
                   is_terminal = .true.
                   exit search_obj
                end if
                
             end if if_squarebrack
          end if if_curlybrack
       end if if_doublequote

       if (i>vardata_max_len)then
          write(log_unit,*) 'GetNextObject could not close the object within the buffer space.'
          call shr_sys_abort()
       end if
         
    end do search_obj

    next_pos = fpos
      
    return
  end subroutine GetNextObject
    
  ! =====================================================================================
  
  subroutine GetDimensions(pstruct)

    type(params_type)           :: pstruct

    integer,parameter :: dimdata_len = 4096 ! We read in ALL the dimension data
                                            ! from { to } into one string, this
                                            ! is how big this string is
    character(len=max_sl) :: tagname
    integer :: otype
    integer :: tag_pos            ! tag's starting position
    integer :: obj_pos            ! object's first position
    integer :: dim0_pos           ! position at start of dimensions section
    integer :: next_pos
    logical :: is_terminal
    integer :: i               ! local character read position index
    integer :: n_dims
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str
    
    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer such that the last character
    !         the : separator for "dimensions"
    ! -----------------------------------------------------------------------------------
    tagname = '"dimensions"'
    call FindTag(1, tagname, tag_pos, dim0_pos)

    ! -----------------------------------------------------------------------------------
    ! Step 2: Count the number of dimensions
    ! -----------------------------------------------------------------------------------
    n_dims = 0
    next_pos = dim0_pos
    is_terminal = .false.
    do_dimcount: do while(.not.is_terminal)
       tagname = search_next_tag
       call FindTag(next_pos, tagname, tag_pos, obj_pos)
       call GetNextObject(obj_pos,otype,next_pos,is_terminal)
       n_dims = n_dims+1
       if(debug) write(log_unit,*) trim(tagname),trim(obj_buffer)
    end do do_dimcount

    allocate(pstruct%dimensions(n_dims))

    ! -----------------------------------------------------------------------------------
    ! Step 3: Parse the dimension data string for dimension data
    !         Each dimension should have a colon ":".
    !         Also, scan for the "scalar" parameter
    !         !!! Just count for now
    ! -----------------------------------------------------------------------------------
    i = 0
    next_pos = dim0_pos
    is_terminal = .false.
    do_dimget: do while(.not.is_terminal)
       tagname = search_next_tag
       call FindTag(next_pos, tagname, tag_pos, obj_pos)
       call GetNextObject(obj_pos,otype,next_pos,is_terminal)
       i = i + 1
       pstruct%dimensions(i)%name = trim(CleanSymbol(tagname))
       call StringToStringOrReal(trim(obj_buffer),tmp_str,tmp_real,is_num=.true.)
       pstruct%dimensions(i)%size = int(tmp_real)
    end do do_dimget

    if(debug)then
       write(log_unit,*)'--- Dimensions ---'
       do i = 1,n_dims
          write(log_unit,*)trim(pstruct%dimensions(i)%name),pstruct%dimensions(i)%size
       end do
       write(log_unit,*)''
    end if
    
    return
  end subroutine GetDimensions

  ! =====================================================================================

  subroutine ReadCharVar(pstruct,iparm,obj_pos0,param_name,string_scr)

    type(params_type)     :: pstruct
    integer               :: iparm
    integer               :: obj_pos0
    character(len=max_sl) :: param_name
    character(len=max_ll), dimension(*) :: string_scr
    
    type(param_type), pointer :: param
    character(len=vardata_max_len) :: data_str
    integer :: n_vec_out
    integer :: i,j,k
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str
    integer :: dimsizes(2)
    integer :: tag_pos
    integer :: obj_pos
    integer :: next_pos
    logical :: is_terminal
    integer :: otype
    character(len=max_sl) :: tag  ! The keyword for this object
    
    param => pstruct%parameters(iparm)
    param%name = trim(CleanSymbol(param_name))
    param%access_count = 0


    tag = '"dims"'
    call FindTag(obj_pos0,tag, tag_pos, obj_pos)
    call GetNextObject(obj_pos,otype,next_pos,is_terminal)
    call GetStringVec(obj_buffer,string_scr,n_vec_out)
    allocate(param%dim_names(n_vec_out))
    dimsizes(:)=-1
    param%ndims = n_vec_out
    do i = 1,n_vec_out
       param%dim_names(i) = trim(CleanSymbol(string_scr(i)))
       if(index(trim(param%dim_names(i)),'scalar')>0)then
          dimsizes(i) = 1
       else
          dimsizes(i) = pstruct%GetDimSizeFromName(trim(param%dim_names(i)))
       end if
    end do
    call ClearStringScratch(string_scr,n_vec_out)


    tag = '"dtype"'
    call FindTag(obj_pos0,tag, tag_pos, obj_pos)
    call GetNextObject(obj_pos,otype,next_pos,is_terminal)
    call GetStringVec(obj_buffer,string_scr,n_vec_out)
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

    
    tag = '"long_name"'
    call FindTag(obj_pos0,tag, tag_pos, obj_pos)
    call GetNextObject(obj_pos,otype,next_pos,is_terminal)
    call GetStringVec(obj_buffer,string_scr,n_vec_out)
    data_str = trim(string_scr(1))
    param%long_name = trim(CleanSymbol(data_str))
    call ClearStringScratch(string_scr,n_vec_out)


    tag = '"units"'
    call FindTag(obj_pos0,tag, tag_pos, obj_pos)
    call GetNextObject(obj_pos,otype,next_pos,is_terminal)
    call GetStringVec(obj_buffer,string_scr,n_vec_out)
    data_str = trim(string_scr(1))
    param%units = trim(CleanSymbol(data_str))
    call ClearStringScratch(string_scr,n_vec_out)


    tag = '"data"'
    call FindTag(obj_pos0,tag, tag_pos, obj_pos)
    call GetNextObject(obj_pos,otype,next_pos,is_terminal)
    call GetStringVec(obj_buffer,string_scr,n_vec_out)
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

  ! ====================================================================================

  subroutine FindTag(pos0, tag, fpos, epos)

    ! This procedure finds the location in the character buffer
    ! of the provided tag/keyword name and the separator ":".
    ! This does not tell us what type of object it is, just the
    ! character index position of the trailing separator.
    
    ! Arguments
    integer, intent(in)           :: pos0 ! Index to start with, find
                                          ! the first key/tag after this
    
    character(len=max_sl), intent(inout)  :: tag  ! The keyword for this object
    
    integer, intent(out)          :: epos ! Index of the first character
                                          ! following the separator
                                          ! for the tag, ie the
                                          ! start position of it's object
    integer, intent(out)          :: fpos ! Index of the first character
                                          ! in the tag
    
    integer :: spos     ! index of the separator in the truncated string buffer
    integer :: ipos     ! current index
    integer :: i
    logical :: open_dq  ! Are double quotes currently open?
    integer :: len_tag  ! length of the trimmed string

    
    if(index(trim(tag),trim(search_next_tag))>0)then
       
       tag  = repeat(' ',max_sl)
       ipos = pos0
       i    = 1
       open_dq = .false.
       search_tag: do
          if (open_dq) then
             i = i + 1
             tag(i:i) = file_buffer(ipos:ipos)
             if (scan(file_buffer(ipos:ipos),'"')>0) then
                exit search_tag
             end if
          end if
          if (scan(file_buffer(ipos:ipos),'"')>0) then
             fpos     = ipos
             tag(i:i) = file_buffer(ipos:ipos)
             open_dq  = .true.
          end if
          ipos = ipos + 1
          if(ipos>(vardata_max_len+pos0)) then
             write(log_unit,*) 'FindTag, could not find next tag from position: ',pos0
             call shr_sys_abort()
          end if
       end do search_tag
    else
       
       fpos = index(file_buffer(pos0:nbuffer),trim(tag))
       !if (fpos .ne. index(file_buffer(pos0:nbuffer),trim(tag),.true.)) then
       !   write(log_unit,*) 'FindTagPos, tag:',trim(tag),' is non-unique?'
       !   call shr_sys_abort()
       !end if
       if (fpos==0) then
          write(log_unit,*) 'FindTagPos, tag:',trim(tag),' not found?'
          call shr_sys_abort()
       end if
       fpos = fpos + pos0 - 1
       
    end if
    
    len_tag = len(trim(tag))
    
    ! Find the fist position of the separator, which must follow
    ! the last index of the tag
    spos = index(file_buffer(fpos+len_tag-1:nbuffer),':')

    if (spos==0) then
       write(log_unit,*) 'FindTagPos, tag:',trim(tag),'no separator found?'
       call shr_sys_abort()
    end if

    epos = spos + fpos + len_tag - 1

    return
  end subroutine FindTag
    
  ! =====================================================================================
  
  subroutine GetParameters(string_scr,pstruct)

    character(len=max_ll), dimension(*), intent(inout) :: string_scr ! Internal scratch space
    type(params_type), intent(inout)                   :: pstruct

    character(len=max_sl) :: tagname
    integer               :: tag_pos
    integer               :: param0_pos
    integer,parameter :: data_len = 4096
    integer :: i
    integer :: otype
    integer :: next_pos
    integer :: n_params
    integer :: obj_pos
    logical :: is_terminal

    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer such that the last character
    !         the : separator for "parameters"
    ! -----------------------------------------------------------------------------------
    tagname = '"parameters"'
    call FindTag(1, tagname, tag_pos, param0_pos)

    n_params = 0
    next_pos = param0_pos
    is_terminal = .false.
    do_paramcount: do while(.not.is_terminal)
       tagname = search_next_tag
       call FindTag(next_pos, tagname, tag_pos, obj_pos)
       call GetNextObject(obj_pos,otype,next_pos,is_terminal)
       n_params = n_params + 1
    end do do_paramcount

    if(debug) write(log_unit,*) 'Found ',n_params,' in the "parameters" group'
    allocate(pstruct%parameters(n_params))

    i = 0
    next_pos = param0_pos
    is_terminal = .false.
    do_paramget: do while(.not.is_terminal)
       tagname = search_next_tag
       call FindTag(next_pos, tagname, tag_pos, obj_pos)
       call GetNextObject(obj_pos,otype,next_pos,is_terminal)
       i = i + 1
       if(debug)write(log_unit,*) "-------------------"
       if(debug)write(log_unit,*) trim(tagname),"    ",trim(obj_buffer)
       ! Fill in the parameter data structures using the object buffer
       call ReadCharVar(pstruct,i,obj_pos,tagname,string_scr)
    end do do_paramget
    

    return
  end subroutine GetParameters

  ! =====================================================================================

  subroutine GetStringVec(string_in,string_vec_out,n_vec_out)

    ! This routine parses a string, does some cleaning of extra characters
    ! and generally uses commas as separators to fill in a vector of
    ! cleaned strings. It is assumed that the string passed in contains
    ! the data, not the symbol.
    !
    ! We typically use this to parse data following the separator ":"
    ! If its a vector, it probably looks like this:  [0, 1, 2, 5, 10, 20, 50],
    ! If its a scalar, it probably looks like this:  5,
    ! It should also parse vectors of length 1 just fine...
    ! First step is to look for a "["
    
    character(len=*),intent(in)                       :: string_in
    character(len=max_ll), dimension(*),intent(inout) :: string_vec_out 
    integer,intent(out)                               :: n_vec_out

    integer :: i,k,l
    integer :: iv    ! index of vector start
    integer :: ivz   ! index of vector end
    logical :: has_brackets
    logical :: open_dq
    integer :: len_string
    
    ! Determine if this is a vector or not
    ! (ie iv=0 is scalar, iv>0 is vector)
    ! Take the right most (in case nested brackets)

    len_string = len(trim(string_in))
    open_dq = .false.
    has_brackets = .false.
    do i = 1,len_string
       if(scan(string_in(i:i),'"')>0)then
          open_dq=.not.open_dq
       end if
       if(.not.open_dq)then
          if(.not.has_brackets)then
             has_brackets = scan(string_in(i:i),"[")>0
          end if
       end if
    end do

    if(has_brackets)then
       ivz = index(trim(string_in),']',.true.) ! pos of last closing bracket
       if(ivz==0)then
          write(log_unit,*)'no closing bracket when one was opened?'
          call shr_sys_abort()
       end if
       iv = index(trim(string_in),'[')
    else
       iv  = 1
       ivz = len(trim(string_in))
    end if
    
    k = 0 ! Counter for the cleaned string position
    l = 1 ! Counter for the output vector
    
    ! Loop through every character of the original string
    open_dq = .false.
    do_str_scan: do i = iv, ivz

       ! commas and brackets inside quotes
       ! should not be allowed to act as separators
       if(scan(string_in(i:i),'"')>0)then
          open_dq=.not.open_dq
       end if
       
       ! If we hit a comma, we either finish our scan
       ! because that was the last character, unless this
       ! is an array, then we update our index and cycle
       if(.not.open_dq) then
          if(scan(string_in(i:i),',')>0) then
             if(k==0)then
                write(log_unit,*) 'Encountered comma before data in string read, aborting'
                call shr_sys_abort()
             end if
             if(has_brackets)then
                l = l + 1
                k = 0
                cycle do_str_scan
             else
                exit do_str_scan
             end if
          end if
       end if
          
       ! Otherwise, we save anything other than a bracket character
       ! into our string array
       
       if (scan(string_in(i:i),'[]') == 0) then
          ! This is a character we want to keep
          k = k + 1
          ! Append the character to the cleaned string
          string_vec_out(l)(k:k) = string_in(i:i)
       end if
       
    end do do_str_scan

    n_vec_out = l

  end subroutine GetStringVec

  ! =====================================================================================
  
  function CleanSymbol(string_in) result(string_out)

    character(len=*) string_in
    character(len=max_sl) :: string_out
    integer :: i,k
    logical :: inside_dq  ! inside a double quote
    
    ! Flush the string
    string_out = repeat(' ',max_sl)
    string_in  = adjustl(string_in)
    k=0
    inside_dq = .false.
    do i = 1, len(trim(string_in))

       if( scan(string_in(i:i),'"')>0 ) then
          inside_dq = .not.inside_dq
       else
          if( scan(string_in(i:i),',')==0 .or. inside_dq ) then
             ! This is a character we want to keep
             k = k + 1
             ! Append the character to the cleaned string
             string_out(k:k) = string_in(i:i)
          end if
       end if

    end do

    if(inside_dq)then
       write(log_unit,*) 'Tried to clean a symbol, but it did not close double quotes?'
       call shr_sys_abort()
    end if
    
    
  end function CleanSymbol

  ! =====================================================================================
  
  subroutine ClearStringScratch(string_scr,n_vec_in)
    character(len=max_ll), dimension(*),intent(inout) :: string_scr
    integer,intent(in) :: n_vec_in
    integer :: i
    do i = 1,n_vec_in
       string_scr(i) = repeat(' ',max_ll)
    end do
    return
  end subroutine ClearStringScratch

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

    ! This gets rid of things like double quotes
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

  ! ==================================================================================

  subroutine Destroy(this)

    ! Deallocate the parameter and dimension data structure
    
    class(params_type)       :: this
    type(param_type),pointer :: param_ptr
    integer                  :: i

    nullify(param_ptr)
    loop_params: do i = 1,size(this%parameters)
       param_ptr=>this%parameters(i)
       if(allocated(param_ptr%dim_names))deallocate(param_ptr%dim_names)
       if(allocated(param_ptr%r_data_1d))deallocate(param_ptr%r_data_1d)
       if(allocated(param_ptr%r_data_2d))deallocate(param_ptr%r_data_2d)
       if(allocated(param_ptr%i_data_1d))deallocate(param_ptr%i_data_1d)
       if(allocated(param_ptr%i_data_2d))deallocate(param_ptr%i_data_2d)
       if(allocated(param_ptr%c_data_1d))deallocate(param_ptr%c_data_1d)
       if(allocated(param_ptr%c_data_2d))deallocate(param_ptr%c_data_2d)
    end do loop_params
    if(associated(this%parameters))deallocate(this%parameters)
    if(associated(this%dimensions))deallocate(this%dimensions)

  end subroutine Destroy

    
  ! =====================================================================================
  
  subroutine JSONDumpParameter(param)
    
    type(param_type) :: param
    integer :: i
    integer :: j
    
    write(log_unit,*) '----------------------------'
    write(log_unit,'(A14,A)') 'Parameter: ',trim(param%name)
    write(log_unit,'(A14,A)') '  units: ',trim(param%units)
    write(log_unit,'(A14,A)') '  long_name: ',trim(param%long_name)
    write(log_unit,'(A14,I1)')'  dtype: ',param%dtype
    do i=1,param%ndims
       write(log_unit,'(A14,I1,A,A)') '  dimension: ',i,': ',trim(param%dim_names(i))
    end do
    
    select case(param%dtype)
    case(r_scalar_type)
       write(log_unit,'(A,F15.6)')'   data: ',param%r_data_scalar
    case(i_scalar_type)
       write(log_unit,*)'   data: ',param%i_data_scalar
    case(c_solo_type)
       write(log_unit,*)trim(param%c_data)
    case(r_1d_type)
       do i = 1,size(param%r_data_1d,dim=1)
          write(log_unit,'(A,F15.6)')'   data: ',param%r_data_1d(i)
       end do
    case(i_1d_type)
       do i = 1,size(param%i_data_1d,dim=1)
          write(log_unit,*)'   data: ',param%i_data_1d(i)
       end do
    case(c_1d_type)
       do i = 1,size(param%c_data_1d,dim=1)
          write(log_unit,*)'   data: ',trim(param%c_data_1d(i))
       end do
    case(r_2d_type)
       do i = 1,size(param%r_data_2d,dim=1)
          do j = 1,size(param%r_data_2d,dim=2)
             write(log_unit,'(A,F15.6)')'   data: ',param%r_data_2d(i,j)
          end do
       end do
    case(i_2d_type)
       do i = 1,size(param%i_data_2d,dim=1)
          do j = 1,size(param%i_data_2d,dim=2)
             write(log_unit,*)'   data: ',param%i_data_2d(i,j)
          end do
       end do
    case(c_2d_type)
       do i = 1,size(param%c_data_2d,dim=1)
          do j = 1,size(param%c_data_2d,dim=2)
             write(log_unit,*)'   data: ',trim(param%c_data_2d(i,j))
          end do
       end do
    end select
    
  end subroutine JSONDumpParameter

  
end module JSONParameterUtilsMod
