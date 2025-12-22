program FatesGreatCircle

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesUtilsMod,     only : GreatCircleDist

  implicit none

  ! LOCALS:
  real(r8) :: inv_lat_list(5)      ! list of lat coords
  real(r8) :: inv_lon_list(5)      ! list of lon coords
  real(r8) :: site_lat_list(5)      ! list of lat coords
  real(r8) :: site_lon_list(5)      ! list of lon coords
  real(r8) :: delta_site_list(5)
  integer  :: i, s
  integer  :: invsite
  
  inv_lat_list = (/-89._r8, -20._r8, 0._r8, 60._r8, 90._r8/)
  inv_lon_list = (/-170._r8, -20._r8, 120.5_r8, 210._r8, 90._r8/)

  site_lat_list = (/-19._r8, -89._r8, 1._r8, 63._r8, 88._r8/)
  site_lon_list = (/-21._r8, -171._r8, 118.5_r8, 214._r8, 78._r8/)

  do s=1, 5
    do i =1, 5
      delta_site_list(i) = GreatCircleDist(site_lon_list(s), inv_lon_list(i), &
        site_lat_list(s),inv_lat_list(i))
    end do
    invsite = minloc(delta_site_list(:), dim=1)
    write(*,'(A,2(F6.1),A,2(F6.1))') "closest to ", site_lat_list(s), site_lon_list(s), &
      " is: ", inv_lat_list(invsite), inv_lon_list(invsite)
    write(*,'(A,F6.1,A)') "   with distance of:", delta_site_list(invsite)/1000._r8, " km"
  end do
  
end program FatesGreatCircle

! ----------------------------------------------------------------------------------------
