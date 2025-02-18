module FatesUtilsMod
  
  ! This module contains helper functions and subroutines which are general in nature.
  ! Think string parsing, timing, maybe numerics, etc.
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only      : fates_log
  use FatesConstantsMod, only : nearzero
  use FatesGlobals,      only : endrun => fates_endrun
  
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  
  implicit none
  private ! Modules are private by default

  ! Make public necessary subroutines and functions
  public :: check_hlm_list
  public :: check_var_real
  public :: GetNeighborDistance
  public :: FindIndex
  public :: QuadraticRootsNSWC
  public :: QuadraticRootsSridharachary
  public :: ArrayNint
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains
  
  subroutine ArrayNint(realarr,intarr)

    real(r8),intent(in)  :: realarr(:)
    integer,intent(out)  :: intarr(:)
    integer  :: i

    do i = 1,size(realarr,dim=1)
       intarr(i) = nint(realarr(i))
    end do
    
    return
  end subroutine ArrayNint

  ! =============================================================================================
  
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
  
  ! =====================================================================================

  subroutine check_var_real(r8_var, var_name, return_code)

     real(r8),intent(in)         :: r8_var
     character(len=*),intent(in) :: var_name
     integer,intent(out)         :: return_code

     real(r8), parameter :: r8_type  = 1.0
     real(r8), parameter :: overflow  = huge(r8_type)
     real(r8), parameter :: underflow = tiny(r8_type)

     return_code = 0

     ! NaN check
     if (r8_var /= r8_var) then
        write(fates_log(),*) 'NaN detected, ',trim(var_name),': ',r8_var
        return_code = 1
     end if
     
     ! Overflow check (within 100th of max precision)
     if (abs(r8_var) > 0.01*overflow) then
        write(fates_log(),*) 'Nigh overflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 10
     end if

     ! Underflow check (within 100x of min precision)
     if (abs(r8_var) < 100.0_r8*underflow) then
        write(fates_log(),*) 'Nigh underflow detected, ',trim(var_name),': ',r8_var
        return_code = return_code + 100
     end if


  end subroutine check_var_real
  
  !==========================================================================================!
  !   Function to compute the great circle distance between two points: the s suffix denotes !
  ! source point, and f denotes the destination - "forepoint"). The results are given in     !
  ! metres. The formula is intended to be accurate for both small and large distances and    !
  ! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
  ! close to the n*pi/2.                                                                     !
  !------------------------------------------------------------------------------------------!
  real(r8) function GreatCircleDist(slons,slonf,slats,slatf)
  
     use FatesConstantsMod, only : earth_radius_eq &
                                 , rad_per_deg 
     implicit none
    
     !----- Local variables. ----------------------------------------------------------------!
     real(r8), intent(in) :: slons
     real(r8), intent(in) :: slonf
     real(r8), intent(in) :: slats
     real(r8), intent(in) :: slatf
    
     !----- Local variables. ----------------------------------------------------------------!
     real(r8)     :: lons
     real(r8)     :: lonf
     real(r8)     :: lats
     real(r8)     :: latf
     real(r8)     :: dlon
     real(r8)     :: dlat
     real(r8)     :: x
     real(r8)     :: y
     !---------------------------------------------------------------------------------------!
    
     !----- Convert the co-ordinates to double precision and to radians. --------------------!
     lons = slons * rad_per_deg
     lonf = slonf * rad_per_deg
     lats = slats * rad_per_deg
     latf = slatf * rad_per_deg
     dlon = lonf - lons
     dlat = latf - lats
    
     !----- Find the arcs. ------------------------------------------------------------------!
     x    = dsin(lats) * dsin(latf) + dcos(lats) * dcos(latf) * dcos(dlon)
     y    = dsqrt( (dcos(latf)*dsin(dlon)) * (dcos(latf)*dsin(dlon))                         &
                 + (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))                  &
                 * (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)) )
    
     !----- Convert the arcs to actual distance. --------------------------------------------!
     GreatCircleDist = earth_radius_eq*datan2(y,x)
    
     return
 
  end function GreatCircleDist
  
  
  ! ======================================================================================
 
  function GetNeighborDistance(gi,gj,latc,lonc) result(gcd)
   
      integer,  intent(in) :: gi,gj           ! indices of gridcells
      real(r8), intent(in) :: latc(:),lonc(:) ! lat/lon of gridcells
      real(r8) :: gcd
      
      gcd = GreatCircleDist(lonc(gi),lonc(gj), &
                            latc(gi),latc(gj))
   
  end function GetNeighborDistance
   
  ! ====================================================================================== 
 
  function FindIndex(input_string_array,string_to_match) result(array_index)
   
      ! ---------------------------------------------------------------------------------
      ! This simple function is a standin for the intrinsic FINDLOC which is not available
      ! with some compilers such as NAG prior to v7.0.  As with FINDLOC, the
      ! function will return zero if a match is not found.
      !
      ! Limitations compared to FINDLOC:
      !   - Only takes one dimensional arrays
      !   - Only take arrays of characters
      !   - Does not allow masking
      ! ---------------------------------------------------------------------------------

      ! Input and output
      character(len=*), intent(in) :: input_string_array(:)
      character(len=*), intent(in) :: string_to_match
      integer :: array_index
      
      ! Locals
      integer :: i

      ! Initialize return index as zero
      array_index = 0

      ! Loop throught the array and compare strings
      do i = 1, size(input_string_array)
         if (trim(input_string_array(i)) .eq. trim(string_to_match)) then
            array_index = i
         end if
      end do
   
  end function FindIndex

  subroutine QuadraticRootsNSWC(a,b,c,root1,root2,err)

    ! This code is based off of routines from the NSWC Mathematics Subroutine Library
    ! From the NSWC README (https://github.com/jacobwilliams/nswc)
    ! "The NSWC Mathematics Subroutine Library is a collection of Fortran 77 routines
    !  specializing in numerical mathematics collected and developed by the U.S.
    !  Naval Surface Warfare Center.  This software is made available, without cost,
    !  to the general scientific community."
    ! The F77 code was updated to modern fortran by Jacob Williams:
    ! https://jacobwilliams.github.io/polyroots-fortran
    ! The FATES adaptation of this aborts if only imaginary roots are generated


    real(r8),intent(in) :: a , b , c !! coefficients
    real(r8),intent(out) :: root1 ! sr !! real part of first root
    real(r8),intent(out) :: root2 ! lr !! real part of second root
    logical,intent(out)  :: err
    real(r8) :: b1, d, e

    err = .false.
    if ( abs(a)<nearzero ) then
        root2 = 0.0_r8
        if ( b/=0.0_r8 ) root2 = -c/b
        root1 = 0.0_r8
    elseif ( abs(c)>nearzero ) then
        ! compute discriminant avoiding overflow
        b1 = b/2.0_r8
        if ( abs(b1)<abs(c) ) then
            e = a
            if ( c<0.0_r8 ) e = -a
            e = b1*(b1/abs(c)) - e
            d = sqrt(abs(e))*sqrt(abs(c))
        else
            e = 1.0_r8 - (a/b1)*(c/b1)
            d = sqrt(abs(e))*abs(b1)
        endif
        if ( e<0.0_r8 ) then
           ! complex conjugate zeros
           write (fates_log(),*)'error, imaginary roots detected in quadratic solve'
           err = .true.
           call endrun(msg=errMsg(sourcefile, __LINE__))
        else
            ! real zeros
            if ( b1>=0.0_r8 ) d = -d
            root1 = (-b1+d)/a
            root2 = 0.0_r8
            if ( root1/=0.0_r8 ) root2 = (c/root1)/a
        endif
    else
        root2 = 0.0_r8
        root1 = -b/a
    endif

  end subroutine QuadraticRootsNSWC
  
  subroutine QuadraticRootsSridharachary(a,b,c,root1,root2,err)


    real(r8),intent(in) :: a , b , c !! coefficients
    real(r8),intent(out) :: root1 ! sr !! real part of first root
    real(r8),intent(out) :: root2 ! lr !! real part of second root
    logical,intent(out)  :: err
    real(r8) :: d    ! discriminant
    real(r8) :: das  ! sqrt(abs(d))

    err = .false.
    ! If a is 0, then equation is not quadratic, but linear
    if (abs(a) < nearzero ) then
       root2 = 0.0_r8
       if ( abs(b)>nearzero ) root2 = -c/b
       root1 = 0.0_r8
       return
    end if

    d   = b * b - 4._r8 * a * c
    das = sqrt(abs(d))
 
    if (d > nearzero) then

       root1 = (-b + das) / (2._r8 * a)
       root2 = (-b - das) / (2._r8 * a)

    elseif (abs(d) <= nearzero) then

       root1 = -b / (2._r8 * a)
       root2 = root1
    else 

       write (fates_log(),*)'error, imaginary roots detected in quadratic solve'
       err = .true.
       ! Disable this endrun and use the return err to track down
       ! error provenance 
       call endrun(msg=errMsg(sourcefile, __LINE__))
       
    end if
    
  end subroutine QuadraticRootsSridharachary

  ! ====================================================================================== 
end module FatesUtilsMod
