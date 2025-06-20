program FatesHelloWorld

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none

  ! LOCALS:
  ! if you need parameter file leave this here
  type(fates_unit_test_param_reader) :: param_reader ! param reader instance 
  character(len=:), allocatable      :: param_file   ! input parameter file

  ! CONSTANTS:
  ! if you need an output file update name, otherwise delete
  character(len=*), parameter :: out_file = 'my_output_file.nc' ! output file

  ! Command line arguments read in here -----
  ! uncomment if you want the parameter file or need other arguments
  ! param_file = command_line_arg(1)
  
  ! Read in parameter file here -------
  ! uncomment if you need parameters
  ! call param_reader%Init(param_file)
  ! call param_reader%RetrieveParameters()
  
  print *, "Hello HelloWorld"
  
end program FatesHelloWorld

! ----------------------------------------------------------------------------------------
