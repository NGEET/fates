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

    logical :: astatus
    integer                         :: nargs,ih
    character(len=16),dimension(10) :: args
    
    call parse(hlms,':', args, nargs)
    astatus = .false.
    do ih=1,nargs
       if(trim(args(ih)).eq.trim(hlm_name))then
          astatus = .true.
          return
       end if
    end do
    return
  end function check_hlm_list
  
  
  ! ====================================================================================
  
  
  subroutine parse(str,delims,args,nargs)
    
    ! ----------------------------------------------------------------------------------
    ! Original Code by: George Benthien
    ! Stripped down for simplified use by RGK
    ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    ! the delimiters contained in the string 'delims'. Preceding a delimiter in
    ! 'str' by a backslash (\) makes this particular instance not a delimiter.
    ! The integer output variable nargs contains the number of arguments found.
    ! ---------------------------------------------------------------------------------
    
    character(len=*),intent(in)   :: str
    character(len=*),intent(in)   :: delims
    character(len=len_trim(str))  :: strsav
    character(len=*),dimension(:) :: args
    integer,intent(out)           :: nargs
    integer                       :: i,na,k,lenstr
    
    strsav=str
    na=size(args)
    do i=1,na
       args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(strsav)
    if(lenstr==0) return
    k=0
    do
       if(len_trim(strsav) == 0) exit
       nargs=nargs+1
       call split(strsav,delims,args(nargs))
    end do
    
  end subroutine parse
  
  ! ====================================================================================
  
  subroutine split(str,delims,before)
    
    ! ----------------------------------------------------------------------------------
    ! OriGeorge Benthen
    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the 
    ! found delimiter. 
    ! ----------------------------------------------------------------------------------
    
    character(len=*) :: str,delims,before
    character(len=64) :: strtemp
    character :: ch, cha
    integer   :: lenstr,k,i,iposa,ipos
    
    lenstr=len_trim(str)
    if(lenstr == 0) return        ! string str is empty
    k=0
    before=' '
    do i=1,lenstr
       ch=str(i:i)
       ipos=index(delims,ch)         
       if(ipos == 0) then          ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
       end if
       if(ch /= ' ') then          ! character is a delimiter that is not a space
          strtemp=str(i+1:) 
          str = strtemp
          exit
       end if
       cha=str(i+1:i+1)            ! character is a space delimiter
       iposa=index(delims,cha)
       if(iposa > 0) then          ! next character is a delimiter
          strtemp=str(i+2:) 
          str = strtemp
          exit
       else
          strtemp=str(i+1:) 
          str = strtemp
          exit
       end if
    end do
    
    if(i >= lenstr) str=''
    str=adjustl(str)              ! remove initial spaces
    
    return
    
  end subroutine split
  
end module FatesUtilsMod
