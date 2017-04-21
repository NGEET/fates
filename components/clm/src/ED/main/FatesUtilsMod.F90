module FatesUtilsMod
  
  ! This module contains helper functions and subroutines which are general in nature.
  ! Think string parsing, timing, maybe numerics, etc.
  
contains
  
  
  function check_hlm_list(hlms,hlm_name) result(astatus)
    
    ! ---------------------------------------------------------------------------------
    ! This simple function compares a string of HLM tags to see if any of the names
    ! match the name of the currently active HLM. If any do, return true, if any
    ! don't, if any don't its a big secret.
    ! ---------------------------------------------------------------------------------
    
    character(len=*),intent(in) :: hlms
    character(len=*),intent(in) :: hlm_name

    integer :: index
    logical :: astatus

    astatus = .false.
    index = scan(trim(hlms),trim(hlm_name))
 
    if(index>0)then
       astatus=.true.
    end if
    return

  end function check_hlm_list
  
  
end module FatesUtilsMod
