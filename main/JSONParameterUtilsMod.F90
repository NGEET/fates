module JSONParameterUtilsMod


  ! TO-DO!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Add methods that will check to ignore/swap out brackets }{
  ! that are inside quotes with [].
  ! These brackets are crucial for parsing the file
  ! and if they are in quotes, the code will get confused.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! This module holds the data types and the routines used
  ! for scanning a JSON parameter file, and storing the contents

  use, intrinsic :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END,IOSTAT_EOR

  implicit none
  
  private

  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  
  ! Data types
  integer, parameter :: r_scalar_type = 1
  integer, parameter :: r_1d_type     = 2
  integer, parameter :: r_2d_type     = 3
  integer, parameter :: c_1d_type     = 4

  logical, parameter :: debug = .false.

  integer, parameter :: max_ll = 256    ! Maximum allowable line-length
  integer, parameter :: max_sl = 128    ! Maximum allowable symbol-length
  integer, parameter :: max_ul = 64     ! Maximum characters for units
  integer, parameter :: max_scr = 500   ! Maximum vector size for strings
  integer, parameter :: too_many_iter = 1000

  integer, parameter :: count_phase = 1 ! We read variables twice so that we can allocate
  integer, parameter :: fill_phase  = 2  

  ! Numeric value to represent invalid values (like NaN and Null)
  ! Overridable, see procedure below
  real(r8) :: r_invalid = -1.e36_r8
  
  type,public :: dim_type
     character(len=max_sl) :: name
     integer               :: size
  end type dim_type
  
  type,public ::  param_type
     character(len=max_sl) :: name        ! The variable symbol (the dictionary key)
     character(len=max_ul) :: units
     character(len=max_ll) :: long_name
     integer :: dtype                                    ! Data type, see above
     character(len=max_sl), allocatable :: dim_names(:)  ! These are the indices of the dimensions
                                                         ! associated with this variable
     integer               :: ndims   ! Number of dimensions, same as size (dim_names)
     real(r8)              :: r_data_scalar
     real(r8), allocatable :: r_data_1d(:)
     real(r8), allocatable :: r_data_2d(:,:)
     character(len=128), allocatable :: c_data_1d(:)
  end type param_type

  type,public :: params_type
     type(param_type), pointer :: parameters(:)
     type(dim_type), pointer :: dimensions(:)
   contains
     procedure :: GetDimSizeFromName
     procedure :: GetParamFromName
  end type params_type
  
  public :: ReadJSON
  public :: SetInvalid
  
contains

  subroutine SetInvalid(r_invalid_in)
    real(r8),intent(in) :: r_invalid_in
    r_invalid = r_invalid_in
  end subroutine SetInvalid

  
  subroutine ReadJSON(filename,file_unit,pstruct)

    ! This is esssentially the driver for reading through the input parameter file
    ! There should be two meaningful sections (not including attributes)
    ! 1) dimensions
    ! 2) variables
    !
    ! All variables must be scalar or associated with a known dimension.
    !
    ! Each section is started with an open bracket '{' and closed with a closed bracket '}'
    
    character(len=*),intent(in) :: filename
    integer,intent(in)          :: file_unit
    type(params_type)           :: pstruct
    
    ! Local

    integer :: io_status
    integer :: i

    ! scratch space for storing data as strings
    ! before its copied into data structures
    character(len=max_ll), dimension(max_scr) :: string_scr 
    
    ! Flush the scratch string
    call ClearStringScratch(string_scr,max_scr)
    
    open(unit=file_unit, file=filename, status='old', &
        action='READ', iostat=io_status)
    
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not open parameter file: ', trim(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'successfully opened ',trim(filename)
    end if

    call GetDimensions(file_unit,pstruct)

    call GetVariables(file_unit,string_scr,pstruct)
    
    close(unit=file_unit, iostat=io_status)
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not close file: ', TRIM(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'sucessfully closed ',trim(filename)
    end if
    
  end subroutine ReadJSON

  ! =====================================================================================
  
  subroutine GetDimensions(file_unit,pstruct)

    integer,intent(in)          :: file_unit
    type(params_type)           :: pstruct

    integer,parameter :: scan_len = 128
    integer,parameter :: dimdata_len = 4096
    character(len=max_sl) :: group_str
    character(len=dimdata_len) :: dimdata_str
    character(len=1) :: filechar
    character(len=max_sl) :: symb_str
    character(len=max_sl) :: data_str
    logical :: found_dimtag
    logical :: found_close
    logical :: found_dims
    logical :: found_scalar
    integer :: i
    integer :: sep_id,beg_id,end_id
    integer :: io_status
    integer :: n_dims
    logical :: is_num
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str
    
    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         "dimensions:"
    ! -----------------------------------------------------------------------------------
    found_dimtag = .false.
    group_str = ''
    i=0
    do while(.not.found_dimtag)
       i=i+1
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status) filechar
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,'dimensions')>0 .and. index(group_str,'{',.true.)==max_sl)then
          found_dimtag = .true.
       end if
    end do
    
    ! -----------------------------------------------------------------------------------
    ! Step 2: Read in all the data until the closing bracket
    ! -----------------------------------------------------------------------------------
    found_close = .false.
    dimdata_str = ''
    i=0
    do while(.not.found_close)
       i=i+1
       read(unit=file_unit, fmt='(A1)',ADVANCE='NO',iostat=io_status) filechar
       if(i<=dimdata_len)then
          dimdata_str(i:i) = filechar
       else
          write(*,*) 'Ran out of room reading in dimension data'
       end if
       if(index(dimdata_str,'}')>0)then
          found_close = .true.
       end if
    end do

    ! -----------------------------------------------------------------------------------
    ! Step 3: Parse the dimension data string for dimension data
    !         Each dimension should have a colon ":".
    !         !!! Just count for now
    ! -----------------------------------------------------------------------------------
    i = 0
    beg_id = 1
    found_dims = .false.
    found_scalar = .false.
    do while(.not.found_dims)
       sep_id = index(dimdata_str(beg_id:dimdata_len),':') + beg_id-1 ! Should be left most...
       end_id = index(dimdata_str(beg_id:dimdata_len),',')
       if(end_id==0)then
          end_id = index(dimdata_str(beg_id:dimdata_len),'}')
          if(end_id==0)then
             write(*,*)'Trouble parsing dim data'
             stop
          end if
          found_dims = .true.
       end if

       symb_str = dimdata_str(beg_id:sep_id-1)
       if( trim(CleanSymbol(symb_str))=='scalar') then
          found_scalar = .true.
       end if
       
       end_id = end_id + beg_id -1
       i=i+1
       beg_id = end_id+1
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
       call StringToStringOrReal(data_str,is_num,tmp_str,tmp_real)
       if(.not.is_num)then
          write(*,*)'Failed to read dimension size:'
          write(*,*) trim(symb_str)
          write(*,*) trim(data_str)
          stop
       end if
       pstruct%dimensions(i)%name = trim(CleanSymbol(symb_str))
       pstruct%dimensions(i)%size = int(tmp_real)
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
    type(dim_type), pointer   :: dim
    integer,parameter :: scan_len = 128
    integer,parameter :: vardata_max_len = 4096
    character(len=max_sl) :: group_str
    character(len=vardata_max_len) :: vardata_str
    character(len=1) :: filechar
    character(len=max_sl) :: symb_str
    character(len=max_ll) :: data_str
    integer :: n_vec_out
    logical :: found_close
    logical :: found_vars
    logical :: found_all
    logical :: found_dims
    logical :: found_data
    logical :: found_units
    logical :: found_long
    integer :: i,j,k
    integer :: filepos
    integer :: sep_id,beg_id,end_id,next_sep_id
    integer :: io_status
    integer :: n_vars
    integer :: vardata_len
    logical :: is_num
    real(r8):: tmp_real
    integer :: iter
    character(len=max_sl) :: tmp_str
    integer :: dimsizes(2)
    
    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         the name of the variable":", and save the symbol name
    ! -----------------------------------------------------------------------------------
    found_vartag = .false.
    group_str = ''
    i=0
    do while(.not.found_vartag)
       i=i+1
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status) filechar
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,':')>0 .and. index(group_str,'{',.true.)==min(i,max_sl))then
          found_vartag = .true.
       elseif(index(group_str,'}')>0) then
          ! We found a closing bracket before a variable was defined
          ! this means we are out of variables
          return
       end if
       if(i>too_many_iter)then
          write(*,*)'failed to find variable string or group closing bracket'
          stop
       end if
    end do

    sep_id = index(group_str,':')
    symb_str = trim(CleanSymbol(group_str(1:sep_id-1)))
    
    ! ---------------------------------------------------------------------------------
    ! Step 2: Advance through the file until the next closing bracket. Save
    !         EVERYTHING in a large string.
    ! ---------------------------------------------------------------------------------

    found_close = .false.
    vardata_str = ''
    i=0
    do while(.not.found_close)
       i=i+1
       read(unit=file_unit, fmt='(A1)',ADVANCE='NO',iostat=io_status) filechar
       if(i<=vardata_max_len)then
          vardata_str(i:i) = filechar
       else
          write(*,*) 'Ran out of room reading in dimension data'
          stop
       end if
       if(index(vardata_str,'{')>0)then
          write(*,*)'An open bracket was found nested inside the current parameter:',trim(symb_str)
          write(*,*)'This is an invalid file format for the parameter file'
          stop
       end if
       if(index(vardata_str,'}')>0)then
          found_close = .true.
       end if
    end do
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
    param%name = trim(symb_str)

    beg_id = 1
    found_all   = .false.
    found_dims  = .false.
    found_units = .false.
    found_long  = .false.
    found_data  = .false.
    iter = 0
    do_parse_data: do while(.not.found_all)

       iter = iter+1
       sep_id = index(vardata_str(beg_id:vardata_len),':') + beg_id-1 ! Should be left most...
       next_sep_id = index(vardata_str(sep_id+1:vardata_len),':')
       if(next_sep_id==0)then
          ! This is the last variable, test for closing bracket
          end_id = index(vardata_str(1:vardata_len),'}',.true.)
          if(end_id==0)then
             write(*,*)'Trouble parsing var data'
             stop
          end if
          found_all = .true.
       else
          next_sep_id = next_sep_id + sep_id
          end_id = index(vardata_str(sep_id:next_sep_id),',',.true.) + sep_id-1
       end if

       ! We should not have found the right-most comma before the next separator
       ! or... the closing bracket.  So text between beg_id and sep_id should
       ! have the symbol, and text between sep_id and end_id should have the data

       symb_str = vardata_str(beg_id:sep_id-1)
       data_str = vardata_str(sep_id+1:end_id)

       if( index(trim(symb_str),'"dims"')>0 ) then
          call GetStringVec(data_str,string_scr,n_vec_out)
          allocate(param%dim_names(n_vec_out))
          param%ndims = n_vec_out
          do i = 1,n_vec_out
             param%dim_names(i) = trim(CleanSymbol(string_scr(i)))
             dimsizes(i) = pstruct%GetDimSizeFromName(trim(param%dim_names(i)))
          end do
          call ClearStringScratch(string_scr,n_vec_out)
          found_dims = .true.
       end if
       
       if( index(symb_str,'"long_name"')>0 ) then
          param%long_name = trim(CleanSymbol(data_str))
          found_long = .true.
       end if
       
       if( index(symb_str,'"units"')>0 ) then
          param%units = trim(CleanSymbol(data_str))
          found_units = .true.
       end if
       
       if_data: if( index(symb_str,'"data"')>0 ) then
          call GetStringVec(data_str,string_scr,n_vec_out)
          if(n_vec_out<1)then
             write(*,*) 'parameter data was empty?'
             stop
          end if
          call StringToStringOrReal(string_scr(1),is_num,tmp_str,tmp_real)
          if(is_num)then
             if(param%ndims==1)then
                if(param%dim_names(1)=='scalar')then
                   call StringToStringOrReal(string_scr(1),is_num,tmp_str,tmp_real)
                   param%r_data_scalar = tmp_real
                   param%dtype =  r_scalar_type
                   if(debug)then
                      write(*,*)'-------------------------------'
                      write(*,*)trim(param%name)
                      write(*,*)param%r_data_scalar
                   end if
                else
                   param%dtype =  r_1d_type
                   allocate(param%r_data_1d(dimsizes(1)))
                   if(n_vec_out.ne.dimsizes(1))then
                      write(*,*)'parameter size does not match dimension size'
                      write(*,*) n_vec_out,dimsizes(1)
                      write(*,*) trim(param%name)
                      write(*,*) trim(data_str)
                      write(*,*) len(trim(data_str))
                      stop
                   end if
                   do i=1,dimsizes(1)
                      call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                      param%r_data_1d(i) = tmp_real
                   end do
                   if(debug)then
                      write(*,*)'-------------------------------'
                      write(*,*)trim(param%name)
                      write(*,*)param%r_data_1d(:)
                   end if
                end if
             else
                param%dtype = r_2d_type
                if(n_vec_out.ne.(dimsizes(1)*dimsizes(2)))then
                   write(*,*)'parameter size does not match dimension size'
                   write(*,*) n_vec_out,dimsizes(1),dimsizes(2)
                   write(*,*) trim(param%name)
                   stop
                end if
                allocate(param%r_data_2d(dimsizes(1),dimsizes(2)))
                i=0
                do j=1,dimsizes(1)
                   do k=1,dimsizes(2)
                      i=i+1
                      call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                      param%r_data_2d(j,k) = tmp_real
                   end do
                end do
                if(debug)then
                   write(*,*)'-------------------------------'
                   write(*,*)trim(param%name)
                   write(*,*)param%r_data_2d
                end if
             end if
          else
             ! For string data, just follow the data
             param%dtype = c_1d_type
             allocate(param%c_data_1d(n_vec_out))
             do i=1,n_vec_out
                call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                param%c_data_1d(i) = trim(tmp_str)
             end do
             if(debug)then
                write(*,*)'-------------------------------'
                write(*,*)trim(param%name)
                write(*,*)param%c_data_1d
             end if
          end if
          call ClearStringScratch(string_scr,n_vec_out)
          found_data = .true.
       end if if_data

       ! If we found all the data
       found_all = all( (/ found_data,found_units,found_long,found_dims /) )
       
       if(iter>too_many_iter)then
          write(*,*)'couldnt resolve variables data'
          stop
       end if
       
       beg_id = end_id+1
    end do do_parse_data

    return
  end subroutine ReadCharVar

  ! =====================================================================================
  
  subroutine GetVariables(file_unit,string_scr,pstruct)

    integer,intent(in)                  :: file_unit
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    type(params_type)                   :: pstruct

    
    integer,parameter :: scan_len = 128
    integer,parameter :: data_len = 4096
    character(len=max_sl) :: group_str
    character(len=data_len) :: vardata_str
    character(len=1) :: filechar
    character(len=max_sl) :: symb_str
    character(len=max_sl) :: data_str
    logical :: found_varstag
    logical :: found_vartag
    logical :: found_close
    logical :: found_vars
    integer :: i,j
    integer :: filepos0
    integer :: sep_id,beg_id,end_id
    integer :: io_status
    integer :: n_vars
    logical :: is_num
    real(r8):: tmp_real
    character(len=max_sl) :: tmp_str

    ! Step 0: We start at the very beginning of the file, which
    ! allows us to count file positions and return to those file positions
    rewind(file_unit)

    ! -----------------------------------------------------------------------------------
    ! Step 1: Advance the file's character pointer to the open bracket following
    !         "variables:"
    ! -----------------------------------------------------------------------------------
    found_vartag = .false.
    group_str = ''
    i=0
    do while(.not.found_vartag)
       i=i+1
       read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=io_status) filechar
       if(i<=max_sl)then
          group_str(i:i) = filechar
       else
          call PopString(group_str,filechar)
       end if
       if(index(group_str,'variables')>0 .and. index(group_str,'{',.true.)==max_sl)then
          found_vartag = .true.
       end if
    end do

    ! Remember this exact file position, we will need to return
    filepos0 = i
    
    ! -----------------------------------------------------------------------------------
    ! Step 2: Start a loop through variables, with each one we are identifying
    ! the { and } brackets and saving that in a string.
    ! -----------------------------------------------------------------------------------
    found_vartag = .true.
    n_vars=0
    do while(found_vartag)
       call ReadCharVar(file_unit,count_phase,n_vars,string_scr,found_vartag)
    end do

    if(debug) write(*,*) 'Found ',n_vars,' variables'
    allocate(pstruct%parameters(n_vars))

    call GotoPos(file_unit,filepos0)
    
    do i=1,n_vars
       j = i
       call ReadCharVar(file_unit,fill_phase,j,string_scr,found_vartag,pstruct)
    end do

    return
  end subroutine GetVariables

  ! =====================================================================================

  subroutine GetStringVec(string_in,string_vec_out,n_vec_out)

    ! This routine parses a string, does some cleaning of extra charachters
    ! and generally uses commas as separators to fill in a vector of
    ! cleaned strings. It is assumed that the sting passed in contains
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
    integer :: i,j,k,l,ii
    integer :: iv    ! index of vector start

    ! Determine if this is a vector or not
    ! (ie iv=0 is scalar, iv>0 is vector)
    iv = index(trim(string_in),'[')

    k = 0 ! Counter for the cleaned string position
    l = 1 ! Counter for the output vector
    
    ! Loop through every character of the original string

    do_str_scan: do i = iv+1, len(trim(string_in))

       ! We are done parsing if we hit one of two conditions..
       ! 1) We are in vector mode and we encounter the end bracket
       ! 2) We are not in vector mode (no opening bracket) and we encounter a comma
       if(iv>0 .and. index(string_in(i:i),']')>0) exit do_str_scan
       if(iv==0 .and. index(string_in(i:i),',')>0) exit do_str_scan

       ! If we encounter a comma, we go to the next index
       ! and cycle
       if(index(string_in(i:i),",")>0) then
          if(k==0)then
             write(*,*) 'Encountered comma before data in string read, aborting'
             stop
          end if
          l = l + 1
          k = 0
          cycle do_str_scan
       else
          
          ! SCAN returns the position of the *first* character from the 
          ! set (unwanted_chars) that is present in the string.
          ! If it returns 0, the character is NOT unwanted.
          j = scan(string_in(i:i),'!#&^')
          
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
    integer :: i,j,k
    
    ! Flush the string
    string_out = ''
    string_in  = adjustl(string_in)
    
    k=0
    do i = 1, len(trim(string_in))

       ! SCAN returns the position of the *first* character from the 
       ! set (unwanted_chars) that is present in the string.
       ! If it returns 0, the character is NOT unwanted.
       j = scan(string_in(i:i),'!#&^"@/><}{')
       
       if (j == 0) then
          ! This is a character we want to keep
          k = k + 1
          ! Append the character to the cleaned string
          string_out(k:k) = string_in(i:i)
       end if

       ! remove trailing commas
       j = scan(string_out,',',.true.)  ! right most occurance...
       if(j==len(trim(string_out)) .and. j>0)then
          !write(*,*) j,trim(string_out),len(trim(string_out))
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
    integer :: i        ! file position number index
    integer :: iostat   ! i/o status of read
    
    rewind(file_unit)

    if (filepos > 1) then
       do i = 1, filepos
          ! empty read
          read(unit=file_unit,fmt='(A1)',ADVANCE='NO', iostat=iostat) dummychar
          !write(*,*) dummychar
          if (iostat == 0) then
             ! Success
          else if (iostat == IOSTAT_EOR) then
             ! help debugging and writing lines
          else if(iostat<0)then
             write(*,*) 'Error or EOF reached while skipping lines.'
             stop
          end if
       end do
    end if
    
  end subroutine GotoPos

  ! =====================================================================================
  
  subroutine GotoLine(file_unit,line_number)

    ! THIS ISN'T USED, BUT KEEPING IT HERE IN CASE
    ! This will reset the file read pointer
    ! to read the line of the input argument line_number.
    ! Therefore the next line to be read will be the one
    ! after line number.
    
    integer            :: file_unit
    integer,intent(in) :: line_number

    character(len=256) dummy_line
    integer :: i        ! line number index
    integer :: iostat   ! i/o status of read
    
    rewind(file_unit)

    if (line_number > 1) then
       do i = 1, line_number
          ! empty read
          read(unit=file_unit, fmt='(A)', iostat=iostat) dummy_line
          if (iostat /= 0) then
             write(*,*) 'Error or EOF reached while skipping lines.'
             stop
          end if
       end do

       if(debug) write(*,*) 'gotoline last read: ',trim(dummy_line)
       
    end if
    
  end subroutine GotoLine

  ! =====================================================================================
  
  subroutine StringToStringOrReal(string_in,is_num,out_string,out_real)

    ! This routine takes a string, and determines if it is a number
    ! by trying to read the string into a floating point number
    ! It uses the iostat report to determine if it was possible or not
    ! If the string was not a number, it returns a cleaned version of the
    ! string. It was a number, it returns that as well.
    
    character(len=*),intent(in) :: string_in
    logical,intent(out)         :: is_num
    character(len=max_sl),intent(out) :: out_string
    character(len=max_sl)       :: tmp_str
    character(len=max_sl)       :: tmp_lc_str
    real(r8),intent(out) :: out_real
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
       out_string = 'invalid real'
       is_num = .true.
       
    else
       
       read (unit=tmp_str,fmt=*,iostat=io_status) tmp_real
       if (io_status == 0) then
          out_real = tmp_real
          out_string = 'real number'
          is_num = .true.
       else
          out_real = -9999.9_r8
          out_string = trim(CleanSymbol(string_in))
          is_num = .false.
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
       write(*,*)'could not find a size for unknown dimension: ',trim(dim_name)
       stop
    end if
  end function GetDimSizeFromName

  ! =====================================================================================

  function GetParamFromName(this,param_name) result(param_ptr)

    class(params_type)       :: this
    character(len=*)         :: param_name
    type(param_type),pointer :: param_ptr
    integer                  :: i
    logical                  :: found_param

    nullify(param_ptr)
    found_param = .false.
    loop_params: do i = 1,size(this%parameters)
       if(trim(this%parameters(i)%name).eq.trim(param_name))then
          param_ptr=>this%parameters(i)
          found_param = .true.
       end if
    end do loop_params

    if(.not.found_param)then
       write(*,*)'Error finding parameter by name'
       write(*,*)'Cant find: ',trim(param_name)
       stop
    end if
    
  end function GetParamFromName
    
end module JSONParameterUtilsMod
