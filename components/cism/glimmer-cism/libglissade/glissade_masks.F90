!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_masks.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This module contains routines for computing various masks used by the Glissade 
! velocity solver.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_masks

    use glimmer_global, only: dp
    use glimmer_physcon, only: rhoi, rhoo
    use glissade_grid_operators     
    use glide_types  ! grounding line options
!    use parallel

    implicit none

    private
    public :: glissade_get_masks, glissade_grounded_fraction

  contains

!****************************************************************************

  subroutine glissade_get_masks(nx,          ny,         &
                                thck,        topg,       &
                                eus,         thklim,     &
                                ice_mask,    floating_mask, &
                                ocean_mask,  land_mask)
                                  

    !----------------------------------------------------------------
    ! Compute various masks for Glissade dycore.
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny                ! number of grid cells in each direction

    ! Default dimensions are meters, but this subroutine will work for
    ! any units as long as thck, topg, eus and thklim have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck,                 &! ice thickness (m)
       topg                   ! elevation of topography (m)

    real(dp), intent(in) ::  &
       eus,                  &! eustatic sea level (m), = 0. by default
       thklim                 ! minimum ice thickness for active cells (m)

    integer, dimension(nx,ny), intent(out) ::  &
       ice_mask               ! = 1 if thck > thklim, else = 0  

    integer, dimension(nx,ny), intent(out), optional ::  &
       floating_mask,        &! = 1 if thck > thklim and ice is floating, else = 0
       ocean_mask,           &! = 1 if topg is below sea level and thk <= thklim, else = 0
       land_mask              ! = 1 if topg is at or above sea level

    !----------------------------------------------------------------
    ! Local arguments
    !----------------------------------------------------------------

    integer :: i, j

    !----------------------------------------------------------------
    ! Compute masks in cells
    !----------------------------------------------------------------

    do j = 1, ny
       do i = 1, nx

          if (thck(i,j) > thklim) then
             ice_mask(i,j) = 1
          else
             ice_mask(i,j) = 0
          endif

          if (present(ocean_mask)) then
             if (topg(i,j) < eus .and. thck(i,j) <= thklim) then
                ocean_mask(i,j) = 1
             else
                ocean_mask(i,j) = 0
             endif
          endif

          if (present(floating_mask)) then
             if (topg(i,j) - eus < (-rhoi/rhoo)*thck(i,j) .and. thck(i,j) > thklim) then
                floating_mask(i,j) = 1
             else
                floating_mask(i,j) = 0
             endif
          endif

          if (present(land_mask)) then
             if (topg(i,j) >= eus) then
                land_mask(i,j) = 1
             else
                land_mask(i,j) = 0
             endif
          endif

       enddo
    enddo

  end subroutine glissade_get_masks

!****************************************************************************

  subroutine glissade_grounded_fraction(nx,          ny,        &
                                        thck,        topg,      &
                                        eus,         ice_mask,  &
                                        whichground, f_ground)

    !----------------------------------------------------------------
    ! Compute fraction of ice that is grounded.
    ! This fraction is computed at vertices based on the thickness and
    !  topography of the four neighboring cell centers.
    !
    ! Three cases, based on the value of whichground:
    ! (0) HO_GROUND_NO_GLP: f_ground = 0 or 1 based on flotation criterion
    ! (1) HO_GROUND_GLP: 0 <= f_ground <= 1 based on grounding-line parameterization
    !                    (similar to that of Pattyn 2006)
    ! (2) HO_GROUND_ALL: f_ground = 1 for all cells with ice
    !----------------------------------------------------------------
    
    !TODO: Apply this subroutine in MISMIP test cases

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny                ! number of grid cells in each direction

    ! Default dimensions are meters, but this subroutine will work for
    ! any units as long as thck and topg have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck,                 &! ice thickness (m)
       topg                   ! elevation of topography (m)

    real(dp), intent(in) :: &
       eus                    ! eustatic sea level (= 0 by default)

    integer, dimension(nx,ny), intent(in) ::   &
       ice_mask               ! = 1 for cells where ice is present (thk > thklim), else = 0

    integer, intent(in) ::   &
       whichground            ! option for computing f_ground

    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       f_ground               ! grounded ice fraction at vertex, 0 <= f_ground <= 1
                              ! set to -1 where vmask = 0

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
           
    integer :: i, j

    integer, dimension(nx-1,ny-1) ::   &
       vmask                  ! = 1 for vertices of cells where ice is present (thk > thklim), else = 0

    real(dp), dimension(nx,ny) :: &
       fpat           ! Pattyn function, -rhoo*(topg-eus) / (rhoi*thck)

    real(dp), dimension(nx-1,ny-1) :: &
       stagfpat       ! fpat interpolated to staggered grid

    real(dp) :: a, b, c, d       ! coefficients in bilinear interpolation
                                 ! f(x,y) = a + b*x + c*y + d*x*y

    real(dp) :: f1, f2, f3, f4   ! fpat at different cell centers

    real(dp) ::  &
       var,     &! combination of fpat terms that determines regions to be integrated
       fpat_v    ! fpat interpolated to vertex

    integer :: nfloat     ! number of grounded vertices of a cell (0 to 4)

    logical, dimension(nx,ny) :: &
       cfloat              ! true if fpat > 1 at cell center, else = false

    logical, dimension(2,2) ::   &
       logvar              ! set locally to float or .not.float, depending on nfloat

    real(dp) ::   &
       f_corner, &         ! fractional area in a corner region of the cell
       f_corner1, f_corner2,  &
       f_trapezoid         ! fractional area in a trapezoidal region of the cell

    logical :: adjacent   ! true if two grounded vertices are adjacent (rather than opposite)

    real(dp), parameter :: &
       eps10 = 1.d-10     ! small number

    !WHL - debug
    integer, parameter :: it = 15, jt = 15

    !----------------------------------------------------------------
    ! Compute ice mask at vertices (= 1 if any surrounding cells have ice)
    !----------------------------------------------------------------

    do j = 1, ny-1
       do i = 1, nx-1
          if (ice_mask(i,j+1)==1 .or. ice_mask(i+1,j+1)==1 .or.   &
              ice_mask(i,j)  ==1 .or. ice_mask(i+1,j)  ==1 ) then
             vmask(i,j) = 1
          else
             vmask(i,j) = 0
          endif
       enddo
    enddo


    ! initialize f_ground
    ! Choose a special non-physical value; this value will be overwritten in all cells with ice
    !TODO Choose a different special value?
!    f_ground(:,:) = -1.d0
    f_ground(:,:) = 9.d0

    select case(whichground)

    case(HO_GROUND_NO_GLP)   ! default: no grounding-line parameterization
                             ! f_ground = 1 if fpat <=1, f_ground = 0 if fpat > 1
                             ! Note: Ice is considered grounded at the GL.

       ! Compute Pattyn function at cell centers

       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then
                fpat(i,j) = -rhoo*(topg(i,j) - eus) / (rhoi*thck(i,j))
             else
                fpat(i,j) = 0.d0
             endif
          enddo
       enddo

       ! Interpolate to staggered mesh

       ! For stagger_margin_in = 1, only ice-covered cells are included in the interpolation.
       ! Will return stagfpat = 0. in ice-free regions

       call glissade_stagger(nx,       ny,         &
                             fpat,     stagfpat,   &
                             ice_mask, stagger_margin_in = 1)

       ! Assume grounded if stagfpat <= 1, else floating

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j)==1) then
                if (stagfpat(i,j) <= 1.d0) then
                   f_ground(i,j) = 1.d0
                else
                   f_ground(i,j) = 0.d0
                endif
             endif
          enddo
       enddo

    case(HO_GROUND_GLP)      ! grounding-line parameterization based on Pattyn (2006, JGR)

       ! Compute Pattyn function at grid cell centers

       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then  ! thck > thklim
                fpat(i,j) = -rhoo*(topg(i,j) - eus) / (rhoi*thck(i,j))
             else
                fpat(i,j) = 0.d0  ! this value is never used
             endif
          enddo
       enddo

       ! Interpolate Pattyn function to staggered mesh
       ! For stagger_margin_in = 1, only ice-covered cells are included in the interpolation.
       ! Returns stagfpat = 0. in ice-free regions

       call glissade_stagger(nx,       ny,         &
                             fpat,     stagfpat,   &
                             ice_mask, stagger_margin_in = 1)

       ! Identify cell centers that are floating

       do j = 1, ny
          do i = 1, nx
             if (fpat(i,j) > 1.d0) then
                cfloat(i,j) = .true.
             else
                cfloat(i,j) = .false.
             endif
          enddo
       enddo

       !WHL - debug
       i = it; j = jt
       print*, 'i, j =', i, j
       print*, 'fpat(i:i+1,j+1):', fpat(i:i+1,j+1)
       print*, 'fpat(i:i+1,j)  :', fpat(i:i+1,j)
       print*, 'cfloat(i:i+1,j+1):', cfloat(i:i+1,j+1)
       print*, 'cfloat(i:i+1,j)  :', cfloat(i:i+1,j)

       ! Loop over vertices, computing f_ground for each vertex with vmask = 1

       do j = 1, ny-1
          do i = 1, nx-1

             if (vmask(i,j) == 1) then  ! ice is present in at least one neighboring cell

                if (ice_mask(i,j+1)==1 .and. ice_mask(i+1,j+1)==1 .and.  &
                    ice_mask(i,j)  ==1 .and. ice_mask(i+1,j)  ==1) then

                   ! ice is present in all 4 neighboring cells; interpolate fpat to find f_ground

                   ! Count the number of floating cells surrounding this vertex

                   nfloat = 0
                   if (cfloat(i,j))     nfloat = nfloat + 1
                   if (cfloat(i+1,j))   nfloat = nfloat + 1
                   if (cfloat(i+1,j+1)) nfloat = nfloat + 1
                   if (cfloat(i,j+1))   nfloat = nfloat + 1

                   !WHL - debug
                   if (i==it .and. j==jt) then
                      print*, ' '
                      print*, 'nfloat =', nfloat
                   endif

                   ! Given nfloat, compute f_ground for each vertex
                   ! First the easy cases...
                
                   if (nfloat == 0) then

                      f_ground(i,j) = 1.d0    ! fully grounded

                   elseif (nfloat == 4) then

                      f_ground(i,j) = 0.d0    ! fully floating

                   ! For the other cases the grounding line runs through the rectangular region 
                   !  around this vertex.
                   ! Using the values at the 4 neighboring cells, we approximate fpat(x,y) as
                   !  a bilinear function f(x,y) = a + bx + cy + dxy over the region.
                   ! To find f_ground, we integrate over the region with f(x,y) <= 1
                   !  (or alternatively, we find f_float = 1 - f_ground by integrating
                   !  over the region with f(x,y) > 1).
                   !  
                   ! There are 3 patterns to consider:
                   ! (1) nfloat = 1 or nfloat = 3 (one cell neighbor is not like the others)
                   ! (2) nfloat = 2 and adjacent cells are floating
                   ! (3) nfloat = 2 and diagonally opposite cells are floating

                   elseif (nfloat == 1 .or. nfloat == 3) then
 
                      if (nfloat==1) then
                         logvar(1:2,1:2) = cfloat(i:i+1,j:j+1)
                      else  ! nfloat = 3
                         logvar(1:2,1:2) = .not.cfloat(i:i+1,j:j+1)
                      endif
                      
                      ! Identify the cell that is not like the others
                      ! (i.e., the only floating cell if nfloat = 1, or the only
                      !  grounded cell if nfloat = 3)
                      !
                      ! Diagrams below are for the case nfloat = 1.
                      ! If nfloat = 3, the F and G labels are switched.

                      if (logvar(1,1)) then      ! no rotation
                         f1 = fpat(i,j)          !   G-----G
                         f2 = fpat(i+1,j)        !   |     |
                         f3 = fpat(i+1,j+1)      !   |     |
                         f4 = fpat(i,j+1)        !   F-----G

                      elseif (logvar(2,1)) then  ! rotate by 90 degrees
                         f4 = fpat(i,j)          !   G-----G
                         f1 = fpat(i+1,j)        !   |     |
                         f2 = fpat(i+1,j+1)      !   |     |
                         f3 = fpat(i,j+1)        !   G-----F

                      elseif (logvar(2,2)) then  ! rotate by 180 degrees
                         f3 = fpat(i,j)          !   G-----F
                         f4 = fpat(i+1,j)        !   |     |
                         f1 = fpat(i+1,j+1)      !   |     |
                         f2 = fpat(i,j+1)        !   G-----G

                      elseif (logvar(1,2)) then  ! rotate by 270 degrees
                         f2 = fpat(i,j)          !   F-----G
                         f3 = fpat(i+1,j)        !   |     |
                         f4 = fpat(i+1,j+1)      !   |     |
                         f1 = fpat(i,j+1)        !   G-----G
                      endif
                      
                      ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
                      ! Note: x is to the right and y is up if the southwest cell is not like the others.
                      !       For the other cases we solve the same problem with x and y rotated.
                      !       The rotations are handled by rotating f1, f2, f3 and f4 above.

                      a = f1
                      b = f2 - f1
                      c = f4 - f1
                      d = f1 + f3 - f2 - f4

                      !WHL - debug
                      if (i==it .and. j==jt) then
                         print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
                         print*, 'a, b, c, d =', a, b, c, d
                      endif

                      ! Compute the fractional area of the corner region 
                      ! (floating if nfloat = 1, grounded if nfloat = 3)
                      !
                      ! Here are the relevant integrals:
                      !
                      ! (1) d /= 0:
                      !     integral_0^x0 {y(x) dx}, where x0   = (1-a)/b
                      !                                    y(x) = (1 - (a+b*x)) / (c+d*x)
                      !     = [bc - ad + d) ln(1 + d(1-a)/(bc)) - (1-a)d] / d^2
                      !
                      ! (2) d = 0:
                      !     integral_0^x0 {y(x) dx}, where x0   = (1-a)/b
                      !                                    y(x) = (1 - (a+b*x)) / c
                      !     = (a-1)(a-1) / (2bc)
                      !
                      ! Note: We cannot have bc = 0, because fpat varies in both x and y

                      if (abs(d) > eps10) then
                         f_corner = ((b*c - a*d + d) * log(1.d0 + d*(1.d0 - a)/(b*c)) - (1.d0 - a)*d) / (d*d)
                      else
                         f_corner = (a - 1.d0)*(a - 1.d0) / (2.d0*b*c)
                      endif

                      if (nfloat==1) then  ! f_corner is the floating area
                         f_ground(i,j) = 1.d0 - f_corner
                      else                 ! f_corner is the grounded area
                         f_ground(i,j) = f_corner
                      endif

                      !WHL - debug
                      if (i==it .and. j==jt) then
                         print*, 'f_corner =', f_corner
                         print*, 'f_ground =', f_ground(i,j)
                      endif

                   elseif (nfloat == 2) then

                      ! first the 4 cases where the 2 grounded cells are adjacent
                      ! We integrate over the trapezoid in the floating part of the cell

                      if (cfloat(i,j) .and. cfloat(i+1,j)) then  ! no rotation
                         adjacent = .true.       !   G-----G
                         f1 = fpat(i,j)          !   |     |
                         f2 = fpat(i+1,j)        !   |     |
                         f3 = fpat(i+1,j+1)      !   |     |
                         f4 = fpat(i,j+1)        !   F-----F

                      elseif (cfloat(i+1,j) .and. cfloat(i+1,j+1)) then  ! rotate by 90 degrees
                         adjacent = .true.       !   G-----F
                         f4 = fpat(i,j)          !   |     |
                         f1 = fpat(i+1,j)        !   |     |
                         f2 = fpat(i+1,j+1)      !   |     |
                         f3 = fpat(i,j+1)        !   G-----F

                      elseif (cfloat(i+1,j+1) .and. cfloat(i,j+1)) then  ! rotate by 180 degrees
                         adjacent = .true.       !   F-----F
                         f3 = fpat(i,j)          !   |     |
                         f4 = fpat(i+1,j)        !   |     |
                         f1 = fpat(i+1,j+1)      !   |     |
                         f2 = fpat(i,j+1)        !   G-----G

                      elseif (cfloat(i,j+1) .and. cfloat(i,j)) then   ! rotate by 270 degrees
                         adjacent = .true.       !   F-----G
                         f2 = fpat(i,j)          !   |     |
                         f3 = fpat(i+1,j)        !   |     |
                         f4 = fpat(i+1,j+1)      !   |     |
                         f1 = fpat(i,j+1)        !   F-----G

                      else   ! the 2 grounded cells are diagonally opposite
                         
                         ! We will integrate assuming the two corner regions lie in the lower left
                         ! and upper right, i.e. one of these patterns:
                         !
                         !   F-----G       G-----F
                         !   |     |       |     |
                         !   |  F  |       |  G  |
                         !   |     |       |     |
                         !   G-----F       F-----G
                         !
                         ! Two other patterns are possible, with corner regions in the lower right
                         ! and upper left; these require a rotation before integrating: 
                         
                         !   G-----F       F-----G
                         !   |     |       |     |
                         !   |  F  |       |  G  |
                         !   |     |       |     |
                         !   F-----G       G-----F
                         !   
                         var = fpat(i+1,j)*fpat(i,j+1) - fpat(i,j)*fpat(i+1,j+1)   &
                             + fpat(i,j) + fpat(i+1,j+1) - fpat(i+1,j) - fpat(i,j+1)
                         if (var >= 0.d0) then   ! we have one of the top two patterns
                            f1 = fpat(i,j)
                            f2 = fpat(i+1,j)
                            f3 = fpat(i+1,j+1)
                            f4 = fpat(i,j+1)
                         else   ! we have one of the bottom two patterns; rotate coordinates by 90 degrees
                            f4 = fpat(i,j)
                            f1 = fpat(i+1,j)
                            f2 = fpat(i+1,j+1)
                            f3 = fpat(i,j+1)
                         endif
                      endif  ! grounded cells are adjacent

                      ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
                      a = f1
                      b = f2 - f1
                      c = f4 - f1
                      d = f1 + f3 - f2 - f4

                      ! Integrate the corner areas

                      !WHL - debug
                      if (i==it .and. j==jt) then
                         print*, 'adjacent =', adjacent
                         print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
                         print*, 'a, b, c, d =', a, b, c, d
                      endif

                      if (adjacent) then

                         ! Compute the area of the floating part of the cell
                         ! Here are the relevant integrals:
                         !
                         ! (1) d /= 0:
                         !     integral_0^1 {y(x) dx}, where y(x) = (1 - (a+b*x)) / (c+d*x)
                         !                                  
                         !     = [bc - ad + d) ln(1 + d/c) - bd] / d^2
                         !
                         ! (2) d = 0:
                         !     integral_0^1 {y(x) dx}, where y(x) = (1 - (a+b*x)) / c
                         !                                 
                         !     = -(2a + b - 2) / (2c)
                         !
                         ! Note: We cannot have c = 0, because the passage of the GL
                         !       through the region from left to right implies variation in y.
                         
                         if (abs(d) > eps10) then
                            f_trapezoid = ((b*c - a*d + d) * log(1 + d/c) - b*d) / (d*d)
                         else
                            f_trapezoid = -(2.d0*a + b - 2.d0) / (2.d0*c)
                         endif

                         f_ground(i,j) = 1.d0 - f_trapezoid

                         !WHL - debug
                         if (i==it .and. j==jt) then
                            print*, 'f_trapezoid =', f_trapezoid
                            print*, 'f_ground =', f_ground(i,j)
                         endif

                      else   ! grounded vertices are diagonally opposite

                         ! bug check: make sure some signs are positive as required by the formulas
                         if (b*c - d*(1.d0-a) < 0.d0) then
                            print*, 'Grounding line error: bc - d(1-a) < 0'
                            stop
                         elseif (b*c < 0.d0) then
                            print*, 'Grounding line error: bc < 0'
                            stop
                         elseif ((b+d)*(c+d) < 0.d0) then
                            print*, 'Grounding line error: (b+d)(c+d) < 0'
                            stop
                         endif

                         ! Compute the combined areas of the two corner regions.
                         ! For the lower left region, the integral is the same as above
                         ! (for the case nfloat = 1 or nfloat = 3, with d /= 0).
                         ! For the upper right region, here is the integral:
                         !
                         !     integral_x1^1 {(1-y(x)) dx}, where x1  = (1-a-c)/(b+d)
                         !                                       y(x) = (1 - (a+b*x)) / (c+d*x)
                         !     = {(bc - ad + d) ln[(bc + d(1-a))/((b+d)(c+d))] + d(a + b + c + d - 1)} / d^2
                         !
                         ! The above integral is valid only if (bc + d(1-a)) > 0.
                         ! If this quantity = 0, then the grounding line lies along two lines,
                         ! x0 = (1-a)/b and y0 = (1-a)/c.
                         ! The lower left area is x0*y0 = (1-a)^2 / (bc).
                         ! The upper right area is (1-x0)*(1-y0) = (a+b-1)(a+c-1) / (bc)
                         !
                         ! Note that this pattern is not possible with d = 0

                         !WHL - debug
                         print*, 'Pattern 3: i, j, bc + d(1-a) =', i, j, b*c + d*(1.d0-a)

                         if (abs(b*c + d*(1.d0-a)) > eps10) then  ! the usual case
                            f_corner1 = ((b*c - a*d + d) * log(1.d0 + d*(1.d0-a)/(b*c)) - (1.d0-a)*d) / (d*d)
                            f_corner2 = ((b*c - a*d + d) * log((b*c + d*(1.d0-a))/((b+d)*(c+d)))  &
                                      + d*(a + b + c + d - 1)) / (d*d)
                         else 
                            f_corner1 = (1.d0 - a)*(1.d0 - a) / (b*c)
                            f_corner2 = (a + b - 1.d0)*(a + c - 1.d0) / (b*c)
                         endif
                         
                         ! Determine whether the central point (1/2,1/2) is grounded or floating.
                         ! (Note: fpat_v /= stagfpat(i,j))
                         ! Then compute the grounded area.
                         ! If the central point is floating, the corner regions are grounded;
                         ! if the central point is grounded, the corner regions are floating.

                         fpat_v = a + 0.5d0*b + 0.5d0*c + 0.25d0*d
                         if (fpat_v > 1.d0) then  ! the central point is floating; corners are grounded
                            f_ground(i,j) = f_corner1 + f_corner2
                         else                     ! the central point is grounded; corners are floating
                            f_ground(i,j) = 1.d0 - (f_corner1 + f_corner2)
                         endif

                         !WHL - debug
                         if (i==it .and. j==jt) then
                            print*, 'fpat_v =', fpat_v
                            print*, 'f_corner1 =', f_corner1
                            print*, 'f_corner2 =', f_corner2
                            print*, 'f_ground =', f_ground(i,j)
                         endif

                      endif  ! adjacent or opposite

                   endif     ! nfloat

                else   ! one or more neighboring cells is ice-free, so bilinear interpolation is not possible
                       ! In this case, set f_ground = 0 or 1 based on stagfpat at vertex

                   if (stagfpat(i,j) <= 1.d0) then
                      f_ground(i,j) = 1.d0
                   else
                      f_ground(i,j) = 0.d0
                   endif

                endif     ! ice_mask = 1 in all 4 neighboring cells
             endif        ! vmask = 1
          enddo           ! i
       enddo              ! j

    case(HO_GROUND_ALL)   ! all vertices with ice-covered neighbors are assumed grounded, 
                          ! regardless of thck and topg

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                f_ground(i,j) = 1.d0
             endif
          enddo
       enddo

    end select

  end subroutine glissade_grounded_fraction

!****************************************************************************

  end module glissade_masks

!****************************************************************************

