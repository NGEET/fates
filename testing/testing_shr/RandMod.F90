module RandMod
  !
  !  DESCRIPTION
  !  this module implements a random number generator for testing purposes
  !
  
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  private
  
  integer :: n = 9876 ! default shared seed
  
  public :: random0
  public :: set_seed
  
  contains
  
  !---------------------------------------------------------------------------------------
  
  real(r8) function random0(minl, maxl)
    !
    ! DESCRIPTION:
    !	Generate a pseudorandom number with a uniform distribution
    !	in the range minl <= ran < maxl.
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: minl ! minimum limit for random number
    real(r8), intent(in) :: maxl ! maximum limit for random number

    ! calculate next number
    n = mod(8121*n + 28411, 134456)

    ! generate random value from this number
    random0 = real(n)/134456.0_r8
    random0 = minl + (maxl - minl)*random0

  end function random0
  
  !---------------------------------------------------------------------------------------
  
  subroutine set_seed(seed)
    !
    ! DESCRIPTION
    !	set the random number seed
    !

    ! ARGUMENTS:
    integer, intent(in) :: seed ! seed

    ! set
    n = abs(seed)

  end subroutine set_seed
  
  !---------------------------------------------------------------------------------------

end module RandMod