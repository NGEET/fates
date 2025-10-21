program TestJSON

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesJSONMod
  implicit none

  integer :: arg_count   ! Number of arguments provided
  character(len=:),allocatable :: filename

  integer :: filename_len
  integer :: file_unit = 18

  type(params_type) :: pstruct
  
  arg_count = command_argument_count()
  if (arg_count < 1) then
     write(*,*) 'ERROR: Missing required argument.'
     write(*,*) 'Usage: ./program_name <input_filename>'
     stop  ! Stop if the filename argument is not provided
  end if
  
  ! Get the length of the first argument to correctly allocate FILENAME
  call get_command_argument(1, length=filename_len)
  allocate(character(len=filename_len) :: filename)
  
  ! Read the first command line argument (index 1) into FILENAME
  call get_command_argument(1, value=filename)

  call ReadJSON(filename,file_unit,pstruct)

  

  deallocate(filename)

end program TestJSON
