module FatesArgumentUtils
  
  implicit none
  private
  
  public :: command_line_arg
  
  contains
    
    function command_line_arg(arg_position)
      ! DESCRIPTION:
      ! Retrieves a single argument from the command line at a given position
      
      ! ARGUMENTS:
      integer,                      intent(in)  :: arg_position     ! argument position
      character(len=:), allocatable             :: command_line_arg ! command argument
      
      ! LOCALS:
      integer :: n_args ! number of arguments received
      integer :: arglen ! argument length
      
      ! get total number of arguments
      n_args = command_argument_count()
      
      if (n_args < arg_position) then
        write(*, '(a, i2, a, i2)') "Incorrect number of arguments: ", n_args, ". Should be at least", arg_position, "."
        stop
      end if 
       
      call get_command_argument(arg_position, length=arglen)
      allocate(character(arglen) :: command_line_arg)
      call get_command_argument(arg_position, value=command_line_arg)
      
    end function command_line_arg
  
end module FatesArgumentUtils
