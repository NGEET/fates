program TestJSON

  use FatesConstantsMod, only : r8 => fates_r8
  use JSONParameterUtilsMod
  implicit none

  integer :: arg_count   ! Number of arguments provided
  character(len=:),allocatable :: filename

  integer :: filename_len
  integer :: i

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

  call JSONSetInvalid(1.e-32_r8)
  call JSONSetLogInit(6)
  call JSONRead(filename,pstruct)

  
  do i=1,size(pstruct%parameters)
     call JSONDumpParameter(pstruct%parameters(i))
  end do

  
  ! Assert some reads

  ! Check the character fields on pft names

  !vcmax25top = [50. 62. 39. 61. 58. 58. 62. 54. 54. 38. 54. 86. 78. 78.


  
  

  deallocate(filename)

end program TestJSON
