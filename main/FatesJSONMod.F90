module FatesJSONMod
  
  ! This module holds the data types and the routines used
  ! for scanning a JSON parameter file, and storing the contents

  !use FatesConstantsMod,only :: max_pft
  !use FatesConstantsMod,only :: max_plant_organs
  use FatesConstantsMod,only : r8 => fates_r8
  use, intrinsic :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END

  public
  
  integer, parameter :: i_scalar_type = 0
  integer, parameter :: i_1d_type = 1
  integer, parameter :: i_2d_type = 2
  integer, parameter :: r_scalar_type = 3
  integer, parameter :: r_1d_type = 4
  integer, parameter :: r_2d_type = 5
  integer, parameter :: c_scalar_type = 6
  integer, parameter :: c_1d_type = 7


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
   
    character(len=*),intent(in) :: filename
    integer,intent(in)          :: file_unit
    
    ! Local
    
    integer :: io_status
    integer :: i
    character(len=256) :: line
    logical :: finished_scan
    
    open(unit=file_unit, file=filename, status='old', &
        action='READ', iostat=io_status)
    
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not open file: ', TRIM(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'successfully opened ',trim(filename)
    end if

    finished_scan = .false.
    do while(.not.finished_scan)
       read(unit=file_unit, fmt='(A)',iostat=io_status) line
       if(io_status==IOSTAT_END) then
          write(*,*) 'completed file read'
          finished_scan = .true.
       end if
       !       write(*,*) trim(line),"---"

       ! Lets find the dimensions
       
       
    end do

    
    close(unit=file_unit, iostat=io_status)
    if (io_status /= 0) THEN
       write(*,*) 'ERROR: Could not close file: ', TRIM(filename)
       write(*,*) 'IOSTAT value: ', io_status
       stop  ! Terminate program gracefully if file cannot be opened
    else
       write(*,*) 'sucessfully closed ',trim(filename)
    end if
    
  end subroutine ScanJSonGroups
  
end module FatesJSONMod
