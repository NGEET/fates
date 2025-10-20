module FatesJSONMod


  ! TO-DO!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Add methods that will check to ignore/swap out brackets }{
  ! that are inside quotes with [].
  ! These brackets are crucial for parsing the file
  ! and if they are in quotes, the code will get confused.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  
  ! This module holds the data types and the routines used
  ! for scanning a JSON parameter file, and storing the contents

  !use FatesConstantsMod,only :: max_pft
  !use FatesConstantsMod,only :: max_plant_organs
  use FatesConstantsMod,only : r8 => fates_r8
  use, intrinsic :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END,IOSTAT_EOR

  implicit none
  
  public

  ! Data types
  integer, parameter :: r_scalar_type = 1
  integer, parameter :: r_1d_type     = 2
  integer, parameter :: r_2d_type     = 3
  integer, parameter :: c_1d_type     = 4

  logical, parameter :: debug = .false.

  integer, parameter :: max_ll = 256    ! Maximum allowable line-length
  integer, parameter :: max_sl = 128    ! Maximum allowable symbol-length
  integer, parameter :: max_ul = 64     ! Maximum characters for units
  integer,parameter  :: too_many_iter = 1000

  integer, parameter :: count_phase = 1 ! We read variables twice so that we can allocate
  integer, parameter :: fill_phase = 2  
  
  type dim_type
     character(len=60) :: name
     integer           :: size
  end type dim_type
  
  type param_type
     character(len=max_sl) :: name        ! The variable symbol (the dictionary key)
     character(len=max_ul) :: units
     character(len=max_ll) :: long_name
     integer :: dtype                                    ! Data type, see above
     character(len=max_sl), allocatable :: dim_names(:)  ! These are the indices of the dimensions
                                                         ! associated with this variable
     integer :: ndims   ! Number of dimensions, same as size (dim_names)
     real(r8)              :: r_data_scalar
     real(r8), allocatable :: r_data_1d(:)
     real(r8), allocatable :: r_data_2d(:,:)
     character(len=128), allocatable :: c_data_1d(:)
 end type param_type

 
 type(param_type), pointer :: fates_params(:)
 type(dim_type),pointer :: param_dimensions(:)

contains
 
  subroutine ScanJSonGroups(filename,file_unit)

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
    
    ! Local

    integer :: io_status
    integer :: i
    character(len=max_ll) :: line
    character(len=max_ll) :: data_str
    character(len=max_sl) :: symb_str
    integer :: data_str_len
    logical :: finished_scan
    integer :: line_num
    integer :: dim_line0 = -1
    integer :: var_line0 = -1
    logical :: finished_dims
    logical :: finished_vars
    logical :: read_complete
    logical :: is_num
    integer :: sep_id
    integer :: tmp_int
    real(r8) :: tmp_real
    character(len=max_sl) :: tmp_str
    integer :: n_dim
    integer :: n_var
    integer :: iter   ! Used for checking against infinite loops
    integer :: ii
    integer :: n_vec_out
    character(len=max_ll), dimension(500) :: string_scr ! scratch space for storing data before its converted
    character(len=30) :: word_vector(1)
    
    finished_dims = .false.
    finished_vars = .false.
    finished_scan = .false.
    
    ! Flush the scratch string
    do i = 1,200
       string_scr(i)=''
    end do
    
    open(unit=file_unit, file=filename, status='old', &
        action='READ', iostat=io_status)
    
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not open file: ', TRIM(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'successfully opened ',trim(filename)
    end if

    call GetDimensions(file_unit) !,param_dimensions)

    call GetVariables(file_unit,string_scr)
    
    line_num = 0
    do_scan: do while(.not.finished_scan)

       ! Just start reading through the file until we find a section of relevance
       read(unit=file_unit, fmt='(A)',iostat=io_status) line
       if(io_status==IOSTAT_END) then
          write(*,*) 'completed file read'
          finished_scan = .true.
       end if
       line_num = line_num + 1
       

       

         
       ! ------------------------------------------------------------------------------
       ! Dimensions
       ! ------------------------------------------------------------------------------
       
       if_dims: if (.not.finished_dims .and. index(line, '"dimensions":') > 0) then

          ! The line number at the beginning of this section, we need to return
          ! to this position, so remember it
          
          dim_line0 = line_num     
          if(debug) write(*, *) 'dimension section found on line ', line_num, ': ', trim(line)

          ! Loop through lines within the dimensions bracket just to get a count
          n_dim = 0
          iter  = 0
          scan_dimensions1: do
             iter = iter+1
             read(unit=file_unit, fmt='(A)',iostat=io_status) line
             if(io_status==IOSTAT_END) then
                write(*,*) 'pre-mature file read completion'
                stop
             end if
             sep_id = index(line,':')
             if(sep_id>0) then
                ! If the bracket comes before the variable, it means this
                ! variable is in the next group
                if(scan(line,'}')>0 .and. index(line,'}')<sep_id ) then
                   exit scan_dimensions1
                end if
                n_dim = n_dim + 1
             end if
             if(index(line,'}')>0) then
                exit scan_dimensions1
             end if
             if(iter>too_many_iter)then
                write(*,*)'Got lost in the dimensions loop'
                stop
             end if
          end do scan_dimensions1
          
          if(debug) write(*,*) 'found ',n_dim,' dimensions'

          ! Performing allocations
          allocate(param_dimensions(n_dim))
          
          ! Going back to top of the dimensions section
          ! (this sets the position in the file to have
          ! just finished reading the line with "dimensions:"
          line_num = dim_line0
          call GotoLine(file_unit,line_num)
          write(*,*) 'RESET LINE NUMBER: ',line_num
          
          ! Loop through lines within the dimensions bracket and
          ! transfer dimension info to data structure
          n_dim = 0
          iter  = 0
          scan_dimensions2: do 

             read(unit=file_unit, fmt='(A)',iostat=io_status) line
             if(io_status==IOSTAT_END) then
                write(*,*) 'pre-mature file read completion'
                stop
             end if
             line_num = line_num + 1
             iter = iter+1

             sep_id = index(line,':')
             if(sep_id>0) then

                ! If the bracket comes before the variable..
                if(scan(line,'}')>0 .and. index(line,'}')<sep_id ) then
                   exit scan_dimensions2
                end if
                
                n_dim = n_dim + 1

                ! lines should have this format:   '   "fates_NCWD": 4,   '
                ! characters to the left of the colon are the symbol
                ! characters to the right of the colon are the value
                symb_str = line(1:sep_id-1)
                data_str = line(sep_id+1:max_ll)
                data_str_len = max_ll-sep_id

                call StringToStringOrReal(data_str,is_num,tmp_str,tmp_real)

                
                !call GetStringVec(data_str,string_scr,n_vec_out)

                !read (unit=string_scr(1), fmt=*, iostat=io_status) tmp_int
                
                !if(n_vec_out>1)then
                !   write(*,*) 'dimension data should be scalar...'
                !   stop
                !end if

                param_dimensions(n_dim)%name = trim(CleanSymbol(symb_str))
                param_dimensions(n_dim)%size = int(tmp_real)
                
                ! Flush the scratch space again
                do ii=1,n_vec_out
                   string_scr(ii)=''
                end do
             end if
                
             if(index(line,'}')>0) then
                exit scan_dimensions2
             end if
             if(iter>too_many_iter)then
                write(*,*)'Got lost in the dimensions loop'
                stop
             end if
          end do scan_dimensions2

          !if(debug) then
             write(*,*)''
             write(*,*)'--- Reporting dimensions ---'
             do ii=1,n_dim
                write(*,*) param_dimensions(ii)%name,param_dimensions(ii)%size
             end do
             write(*,*)''
          !end if
          
          finished_dims = .true.
       end if if_dims

       !stop
       
       if (.not.finished_vars .and.index(line, '"variables":') > 0) then
          
          if(debug) write(*, *) 'variables found on line ', line_num, ': ', trim(line)
          var_line0 = line_num     
          
          ! Loop through lines within the dimensions bracket just to get a count
          n_var = 0
          iter  = 0
          read_complete = .false.
          scan_vars1: do while(.not.read_complete)
             ! This will read everything within the next level of {brackets}
             call ReadVar(file_unit,count_phase,line_num,n_var,read_complete,string_scr)
          end do scan_vars1
          
          allocate(fates_params(n_var))
          n_var =0

          ! rewind so that "variables":{ was the last read
          line_num = var_line0
          call GotoLine(file_unit,line_num)  
          read_complete = .false.
          scan_vars2: do while(.not.read_complete)
             ! This will read everything within the next level of {brackets}
             call ReadVar(file_unit,fill_phase,line_num,n_var,read_complete,string_scr,fates_params(n_var+1))
          end do scan_vars2
          
          write(*,*) 'Found ',n_var,' variables in the file'
          finished_vars = .true.
          finished_scan = .true.
       end if
          
    end do do_scan

    
    close(unit=file_unit, iostat=io_status)
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not close file: ', TRIM(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'sucessfully closed ',trim(filename)
    end if
    
  end subroutine ScanJSonGroups

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
  
  subroutine GetDimensions(file_unit) !,param_dimensions)

    integer,intent(in)          :: file_unit
    !type(dim_type),dimension(:) :: param_dimensions
    
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
    ! -----------------------------------------------------------------------------------
    n_dims = i

    if(found_scalar) then
       allocate(param_dimensions(n_dims))
    else
       allocate(param_dimensions(n_dims+1))
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
       param_dimensions(i)%name = trim(CleanSymbol(symb_str))
       param_dimensions(i)%size = int(tmp_real)
       write(*,*)i,trim(param_dimensions(i)%name),param_dimensions(i)%size
       beg_id = end_id+1
    end do

    if(.not.found_scalar)then
       param_dimensions(n_dims+1)%name = 'scalar'
       param_dimensions(n_dims+1)%size = 1
    end if
    
    
    return
  end subroutine GetDimensions

  ! =====================================================================================

  subroutine ReadCharVar(file_unit,phase,var_num,string_scr,found_vartag,fates_param)

    integer,intent(in)          :: file_unit
    integer,intent(in)          :: phase
    integer,intent(inout)       :: var_num
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    logical,intent(out)         :: found_vartag
    type(param_type),optional   :: fates_param
    
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

    if(phase==count_phase) var_num = var_num + 1
    
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
       ! We've advance the file pointer to the closing bracket
       ! if we are just counting, we are done
       return
    end if

    ! ---------------------------------------------------------------------------------
    ! Step 3: Read through that string. Fill data structures.
    ! ---------------------------------------------------------------------------------

    fates_param%name = trim(symb_str)
    write(*,*) fates_param%name
    !write(*,*) vardata_str(1:vardata_len)
    
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
          allocate(fates_param%dim_names(n_vec_out))
          fates_param%ndims = n_vec_out
          do i = 1,n_vec_out
             fates_param%dim_names(i) = trim(CleanSymbol(string_scr(i)))
             dimsizes(i) = GetDimSizeFromName(trim(fates_param%dim_names(i)))
          end do
          call ClearStringScratch(string_scr,n_vec_out)
          found_dims = .true.
       end if
       
       if( index(symb_str,'"long_name"')>0 ) then
          fates_param%long_name = trim(CleanSymbol(data_str))
          found_long = .true.
       end if
       
       if( index(symb_str,'"units"')>0 ) then
          fates_param%units = trim(CleanSymbol(data_str))
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
             if(fates_param%ndims==1)then
                if(fates_param%dim_names(1)=='scalar')then
                   call StringToStringOrReal(string_scr(1),is_num,tmp_str,tmp_real)
                   fates_param%r_data_scalar = tmp_real
                   fates_param%dtype =  r_scalar_type
                   if(debug)then
                      write(*,*)'-------------------------------'
                      write(*,*)trim(fates_param%name)
                      write(*,*)fates_param%r_data_scalar
                   end if
                else
                   fates_param%dtype =  r_1d_type
                   allocate(fates_param%r_data_1d(dimsizes(1)))
                   if(n_vec_out.ne.dimsizes(1))then
                      write(*,*)'parameter size does not match dimension size'
                      write(*,*) n_vec_out,dimsizes(1)
                      write(*,*) trim(fates_param%name)
                      write(*,*) trim(data_str)
                      write(*,*) len(trim(data_str))
                      stop
                   end if
                   do i=1,dimsizes(1)
                      call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                      fates_param%r_data_1d(i) = tmp_real
                   end do
                   if(debug)then
                      write(*,*)'-------------------------------'
                      write(*,*)trim(fates_param%name)
                      write(*,*)fates_param%r_data_1d(:)
                   end if
                end if
             else
                fates_param%dtype = r_2d_type
                if(n_vec_out.ne.(dimsizes(1)*dimsizes(2)))then
                   write(*,*)'parameter size does not match dimension size'
                   write(*,*) n_vec_out,dimsizes(1),dimsizes(2)
                   write(*,*) trim(fates_param%name)
                   stop
                end if
                allocate(fates_param%r_data_2d(dimsizes(1),dimsizes(2)))
                i=0
                do j=1,dimsizes(1)
                   do k=1,dimsizes(2)
                      i=i+1
                      call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                      fates_param%r_data_2d(j,k) = tmp_real
                   end do
                end do
                if(debug)then
                   write(*,*)'-------------------------------'
                   write(*,*)trim(fates_param%name)
                   write(*,*)fates_param%r_data_2d
                end if
             end if
          else
             ! For string data, just follow the data
             fates_param%dtype = c_1d_type
             allocate(fates_param%c_data_1d(n_vec_out))
             do i=1,n_vec_out
                call StringToStringOrReal(string_scr(i),is_num,tmp_str,tmp_real)
                fates_param%c_data_1d(i) = trim(tmp_str)
             end do
             if(debug)then
                write(*,*)'-------------------------------'
                write(*,*)trim(fates_param%name)
                write(*,*)fates_param%c_data_1d
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
  
  subroutine GetVariables(file_unit,string_scr)

    integer,intent(in)          :: file_unit
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    
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
    integer :: i
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
    
    if(debug)write(*,*) 'Found variables: ',trim(group_str)

    ! -----------------------------------------------------------------------------------
    ! Step 2: Start a loop through variables, with each one we are identifying
    ! the { and } brackets and saving that in a string.
    ! -----------------------------------------------------------------------------------
    found_vartag = .true.
    n_vars=0
    do while(found_vartag)
       call ReadCharVar(file_unit,count_phase,n_vars,string_scr,found_vartag)
    end do

    write(*,*) 'Found ',n_vars,' variables'
    allocate(fates_params(n_vars))

    call GotoPos(file_unit,filepos0)
    
    do i=1,n_vars
       call ReadCharVar(file_unit,fill_phase,n_vars,string_scr,found_vartag,fates_params(i))
    end do
    
    stop
    return
  end subroutine GetVariables
  
       
  subroutine ReadVar(file_unit,read_phase,line_num,var_num,read_complete,string_scr,fates_param)

    integer, intent(in) :: file_unit
    integer, intent(in) :: read_phase  !1=counting, 2=filling
    integer, intent(inout) :: line_num
    integer, intent(inout) :: var_num
    logical, intent(inout) :: read_complete
    character(len=max_ll), dimension(*) :: string_scr ! Internal scratch space
    type(param_type),intent(inout),optional :: fates_param

    ! Read file information for one variable. Either
    ! store information in a data structure, or
    ! make a count and simply move the linepointer
    ! in the file forward.
    !
    ! We keep reading until we first
    ! 1) identify a variable which both a colon ":"
    !    and an open bracket "{"
    ! 2) and then continue reading until that bracket closes.
    ! We update the line number for each read as well.
    ! If the closing bracket is contained in the same line as the
    ! next variable (ie a colon is found after the closing bracket)
    ! then we subtract a line number and call the GoTo function
    !
    ! Example:
    ! ----------------------------------------------------------------
    !
    ! "fates_allom_d2bl3": {
    !  "dims": ["fates_pft"],
    !  "long_name": "Parameter 3 for d2bl allometry",
    !  "units": "unitless",
    !  "data": [0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55,
    !           0.55, 0.55, 0.3417, 0.3417, 0.9948]
    !  },
    !
    ! ---------------------------------------------------------------- 

    logical :: var_open = .true.
    logical :: found_var = .false.
    logical :: found_ob = .false.
    logical :: found_cb = .false.
    character(len=max_ll) :: line
    character(len=max_ll) :: data_str
    character(len=max_sl) :: symb_str
    integer :: sep_id  ! separator (:) character index
    integer :: ob_id   ! opening bracket character index
    integer :: cb_id   ! closing bracket character index
    integer :: cb2_id  ! duplicate closing bracket if two present
    integer :: pc_id   ! previous comma character index
    integer :: io_status
    integer :: iter
    integer :: line_num0
    integer :: line_num_ob,line_num_cb
    integer :: i
    integer :: n_vec_out
    integer, parameter :: too_many_var_iter = 20 
    
    
    ! First lets just find the open and closed brackets, and
    ! record the line numbers...

    if(debug) write(*,*) 'Starting variable read, phase:',read_phase

    found_cb = .false.
    found_ob = .false.
    
    line_num0=line_num
    iter = 0
    do_var0: do while((.not.found_cb) .or. (.not.found_ob))
       
       read(unit=file_unit, fmt='(A)',iostat=io_status) line
       line_num = line_num+1
       
       ob_id = scan(line,'{')  ! open bracket index
       cb_id = scan(line,'}')  ! closed bracket index
       sep_id = scan(line,':', .TRUE.) ! seperator index (right most)

       ! If this is really just a closed bracket with no separator
       ! then there is no more variables to read...
       if(cb_id>0 .and. sep_id==0 .and.(.not.found_ob) )then
          read_complete = .true.
          return
       end if

       ! The variable (symbol name) will be before the bracket
       if(sep_id>0.and. .not.found_ob)then
          ! if the variable is sharing a line with another,
          ! the last character to prune will be the comma
          ! following the closing bracket. Should only be one...
          pc_id = scan(line,',',.TRUE.)
          symb_str = trim(CleanSymbol(line(pc_id+1:sep_id-1)))
          var_num = var_num + 1
       end if
       if(cb_id>0 .and. .not.found_cb .and. found_ob)then
          line_num_cb = line_num
          found_cb = .true.
       end if
       if(ob_id>0 .and. .not.found_ob)then
          line_num_ob = line_num
          found_ob = .true.
       end if
       
       iter = iter+1
       if(iter>too_many_var_iter)then
          write(*,*)'couldnt resolve variable'
          stop
       end if
    end do do_var0
    
    ! Rewind to the line before the opening bracket
    ! so its the next one
    line_num = line_num_ob-1
    call GotoLine(file_unit,line_num)

    do line_num = line_num_ob,line_num_cb
       read(unit=file_unit, fmt='(A)',iostat=io_status) line

       if(read_phase==fill_phase)then
          ! Metadata/Key-words
          !character(len=max_sl) :: name        ! The variable symbol (the dictionary key)
          !character(len=max_ul) :: units
          !character(len=max_ll) :: long_name
          !integer :: data_type
          !character(len=max_sl), allocatable :: dim_names(:)  ! These are the indices of the dimensions
          ! associated with this variable
          
          ! Data Storage (using separate arrays for heterogeneous types)
          !real(r8) :: r_scalar
          !integer  :: i_scalar
          !real(r8), allocatable :: r_data_1d(:)
          !real(r8), allocatable :: r_data_2d(:,:)
          !integer,  allocatable :: i_data_1d(:)
          !integer,  allocatable :: i_data_2d(:)
          !character(len=128), allocatable :: c_data_1d(:)

          ! Keywords: "dims": ["fates_history_size_bins"]
          !           "long_name": "Lower edges for DBH size class bins used in size-resolved cohort history output"
          !           "units": "cm",
          !           "data": [ ]

          if( index(line,'"dims"')>0 ) then
             sep_id = scan(line,':')
             if(sep_id>0)then
                data_str = line(sep_id+1:max_ll)
                call GetStringVec(data_str,string_scr,n_vec_out)
             else
                write(*,*) 'dims should have a separator :'
                stop
             end if
             allocate(fates_param%dim_names(n_vec_out))
             do i = 1,n_vec_out
                fates_param%dim_names(i) = trim(CleanSymbol(string_scr(i)))
             end do
             call ClearStringScratch(string_scr,n_vec_out)
          end if
          if( index(line,'"long_name"')>0 ) then
             sep_id = scan(line,':')
             if(sep_id>0)then
                data_str = line(sep_id+1:max_ll)
             else
                write(*,*) 'long_names should have a separator :'
                stop
             end if
             fates_param%long_name = trim(CleanSymbol(data_str))
          end if
          if( index(line,'"units"')>0 ) then
             sep_id = scan(line,':')
             if(sep_id>0)then
                data_str = line(sep_id+1:max_ll)
             else
                write(*,*) 'long_names should have a separator :'
                stop
             end if
             fates_param%units = trim(CleanSymbol(data_str))
          end if
          if( index(line,'"data"')>0 ) then
             sep_id = scan(line,':')
             if(sep_id>0)then
                data_str = line(sep_id+1:max_ll)
                call GetStringVec(data_str,string_scr,n_vec_out)
                if(n_vec_out<1)then
                   write(*,*) 'parameter data was empty?'
                   stop
                end if
             else
                write(*,*) 'parameter data should have a separator :'
                stop
             end if
             
         !    call StringToStringOrReal(string_scr(1),is_num,tmp_str,tmp_real)
         !    if(is_num) then
         !       allocate(fates_param%
             
             !   allocate(fates_param%data(n_vec_out))
             !   do i = 1,n_vec_out
            
          !      fates_param%data(i) = trim(CleanSymbol(string_scr(i)))
          !   end do
          !   call ClearStringScratch(string_scr,n_vec_out)
          end if
          
          
          
       end if
    end do

    ! check to see if the next variable's symbol is found in the
    ! current line "line_num_cb", and reverse a line if so..
    cb_id = scan(line,'}')  ! closed bracket index
    sep_id = scan(line,':', .TRUE.) ! seperator index (right most)
    cb2_id = scan(line,'}',.TRUE.) ! if this has two closed brackets, cb_id DNE cb2_id

    ! In this scenario, two closing brackets were detected. Signal
    ! that there are no more variables to read
    if(cb_id.ne.cb2_id)then
       read_complete = .true.
    end if
    
    if(sep_id>cb_id)then
       line_num=line_num_cb-1
       call GotoLine(file_unit,line_num)
    end if
    

    return
  end subroutine ReadVar
  

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

  subroutine ClearStringScratch(string_scr,n_vec_in)
    character(len=max_ll), dimension(*),intent(inout) :: string_scr
    integer,intent(in) :: n_vec_in
    integer :: i
    do i = 1,n_vec_in
       string_scr(i) = ''
    end do
    return
  end subroutine ClearStringScratch

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

  subroutine GotoLine(file_unit,line_number)

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

  ! ===============================================================================================
  
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
    real(r8),intent(out) :: out_real
    integer  :: io_status
    real(r8) :: tmp_real

    tmp_str = trim(CleanSymbol(string_in))
    read (unit=tmp_str,fmt=*,iostat=io_status) tmp_real
    if (io_status == 0) then
       out_real = tmp_real
       out_string = 'invalid'
       is_num = .true.
    else
       out_real = -9999.9_r8
       out_string = trim(CleanSymbol(string_in))
       is_num = .false.
    end if
    return
  end subroutine StringToStringOrReal

  ! =====================================================================================

  function GetDimSizeFromName(dim_name) result(dim_size)

    character(len=*) :: dim_name
    integer          :: dim_size
    integer          :: i

    dim_size = -1
    loop_dims: do i = 1,size(param_dimensions)
       if(index(param_dimensions(i)%name,trim(dim_name))>0)then
          dim_size = param_dimensions(i)%size
          exit loop_dims
       end if
    end do loop_dims

    if(dim_size==-1)then
       write(*,*)'could not find a size for unknown dimension: ',trim(dim_name)
       stop
    end if
  end function GetDimSizeFromName
  
    
end module FatesJSONMod
