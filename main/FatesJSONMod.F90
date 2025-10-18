module FatesJSONMod
  
  ! This module holds the data types and the routines used
  ! for scanning a JSON parameter file, and storing the contents

  !use FatesConstantsMod,only :: max_pft
  !use FatesConstantsMod,only :: max_plant_organs
  use FatesConstantsMod,only : r8 => fates_r8
  use, intrinsic :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END

  implicit none
  
  public
  
  integer, parameter :: i_scalar_type = 0
  integer, parameter :: i_1d_type = 1
  integer, parameter :: i_2d_type = 2
  integer, parameter :: r_scalar_type = 3
  integer, parameter :: r_1d_type = 4
  integer, parameter :: r_2d_type = 5
  integer, parameter :: c_scalar_type = 6
  integer, parameter :: c_1d_type = 7

  logical, parameter :: debug = .true.

  integer, parameter :: max_ll = 256    ! Maximum allowable line-length
  integer, parameter :: max_sl = 128    ! Maximum allowable symbol-length

  type dim_type
     character(len=60) :: name
     integer           :: size
  end type dim_type
  
  type param_type

     ! Metadata/Key-words
     character(len=60) :: name        ! The variable symbol (the dictionary key)
     character(len=60) :: units
     character(len=120) :: long_name
     
     integer :: data_type

     integer, allocatable :: dim_ids(:)  ! These are the indices of the dimensions
                                         ! associated with this variable

    ! Data Storage (using separate arrays for heterogeneous types)
    real(r8) :: r_scalar
    integer  :: i_scalar
    real(r8), allocatable :: r_data_1d(:)
    real(r8), allocatable :: r_data_2d(:,:)
    integer,  allocatable :: i_data_1d(:)
    integer,  allocatable :: i_data_2d(:)
    character(len=128), allocatable :: c_data_1d(:)
    
 end type param_type

 type param_group_type
    character(len=64) :: name   ! Group name
    character(len=64), allocatable :: member_names(:)
    type(param_type), allocatable :: members(:)
 end type param_group_type
 
 type(dim_type),pointer :: param_dimensions(:)
 type(param_group_type), pointer :: group_list(:)

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
    logical :: finished_scan = .false.
    integer :: line_num
    integer :: dim_line0 = -1
    logical :: dimensions_found = .false.
    logical :: variables_found = .false.
    logical :: finished_dims = .false.
    integer :: sep_id
    integer :: tmp_int
    integer :: n_dim
    integer :: iter   ! Used for checking against infinite loops
    integer :: ii
    integer :: n_vec_out
    character(len=max_ll), dimension(200) :: string_scr ! scratch space for storing data before its converted
    character(len=30) :: word_vector(1)

    integer,parameter :: too_many_iter = 1000
    
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
          dimensions_found = .true.
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
          call GotoLine(file_unit,dim_line0)
          line_num = dim_line0
          
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
                
                line_num = line_num+1
                n_dim = n_dim + 1

                ! lines should have this format:   '   "fates_NCWD": 4,   '
                ! characters to the left of the colon are the symbol
                ! characters to the right of the colon are the value
                symb_str = line(1:sep_id-1)
                data_str = line(sep_id+1:max_ll)
                data_str_len = max_ll-sep_id

                call GetStringVec(data_str,data_str_len,string_scr,n_vec_out)

                read (unit=string_scr(1), fmt=*, iostat=io_status) tmp_int
                
                if(n_vec_out>1)then
                   write(*,*) 'dimension data should be scalar...'
                   stop
                end if

                param_dimensions(n_dim)%name = trim(CleanSymbol(symb_str))
                param_dimensions(n_dim)%size = tmp_int
                
                ! Flush the scratch space again
                do ii=1,n_vec_out
                   string_scr(ii)=''
                end do
             end if
                
             if(index(line,'}')>0) then
                line_num = line_num+1
                exit scan_dimensions2
             end if
             if(iter>too_many_iter)then
                write(*,*)'Got lost in the dimensions loop'
                stop
             end if
          end do scan_dimensions2

          if(debug) then
             write(*,*)''
             write(*,*)'--- Reporting dimensions ---'
             do ii=1,n_dim
                write(*,*) param_dimensions(ii)%name,param_dimensions(ii)%size
             end do
             write(*,*)''
          end if
          
          finished_dims = .true.
       end if if_dims
       
       if (index(line, '"variables":') > 0) then
          variables_found = .true.
          if(debug) write(*, *) 'variables found on line ', line_num, ': ', trim(line)
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
  
  !function GetIntScalar(string) result(scalarint)
  ! Convert a sting to an integer and clean off any special characters
  ! end function GetIntScalar

  subroutine GetStringVec(string_in,string_len,string_vec_out,n_vec_out)

    ! This routine parses a string, does some cleaning of extra charachters
    ! and generally uses commas as separators to fill in a vector of
    ! cleaned strings
    !
    ! Essentially this is everything in a JSON file after a ":"
    ! If its a vector, it probably looks like this:  [0, 1, 2, 5, 10, 20, 50],
    ! If its a scalar, it probably looks like this:  5,
    ! First step is to look for a "["
    
    character(len=max_ll),intent(in)                  :: string_in
    integer                                           :: string_len
    character(len=max_ll), dimension(:),intent(inout) :: string_vec_out 
    integer,intent(out)                               :: n_vec_out

    !character(len=2)                                :: unwanted_chars  = (/!#/)
    integer :: i,j,k,l,ii
    integer :: iv    ! index of vector start

    ! Determine if this is a vector or not
    ! (ie iv=0 is scalar, iv>0 is vector)
    iv = index(string_in(1:string_len),'[')

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
       
    end do
    
  end function CleanSymbol

  
  
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

       if(debug) write(*,*) 'gotoline: ',trim(dummy_line)
       
    end if
    
  end subroutine GotoLine

    
  
end module FatesJSONMod
