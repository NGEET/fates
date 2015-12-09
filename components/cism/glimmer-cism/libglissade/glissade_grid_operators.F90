!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains various grid operators for the Glissade dycore, including routines 
! for computing gradients and interpolating between staggered and unstaggered grids.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_grid_operators

    use glimmer_global, only: dp
    use glimmer_log
    use glide_types  ! HO_GRADIENT_MARGIN_*
    use parallel

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger,  &
              glissade_centered_gradient, glissade_upstream_gradient,  &
              glissade_edge_gradient, glissade_vertical_average

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              ice_mask,     stagger_margin_in)

    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where values are included in the average, else = 0
                                ! Typically ice_mask = 1 where ice is present (or thck > thklim), else = 0

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating (including zeroes where ice is absent)
                                !   may be appropriate when computing stagusrf and stagthck on land
                                ! 1 = use only values where ice is present
                                !   preferable for tracers (e.g., temperature, flwa) and ocean margins

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 1  ! default is to average only over the cells with ice present
    endif

    stagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          stagvar(i,j) = (var(i,j+1) + var(i+1,j+1) + var(i,j) + var(i+1,j)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over cells with ice present (ice_mask = 1)

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = ice_mask(i,j+1)*var(i,j+1) + ice_mask(i+1,j+1)*var(i+1,j+1)  &
                 + ice_mask(i,j)  *var(i,j)   + ice_mask(i+1,j)  *var(i+1,j)
          summask = real(ice_mask(i,j+1) + ice_mask(i+1,j+1) + ice_mask(i,j) + ice_mask(i+1,j), dp)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

  end subroutine glissade_stagger

!----------------------------------------------------------------------------

  subroutine glissade_unstagger(nx,           ny,          &
                                stagvar,      unstagvar,   &
                                vmask,        stagger_margin_in)

    ! Given a variable on the staggered grid (dimension nx-1, ny-1), interpolate
    ! to find values on the staggered grid (dimension nx, ny).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(out) ::    &
       unstagvar                ! unstaggered field, defined at cell centers

    integer, dimension(nx-1,ny-1), intent(in) ::        &
       vmask                    ! = 1 for vertices where the value is used in the average, else = 0
                                ! Note: The user needs to compute this mask in the calling subroutine.
                                !       It will likely be based on the scalar ice mask, but the details are left open.

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating
                                ! 1 = use only values where vmask = 1

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 1  ! default is to average over cells where vmask = 1
    endif

    unstagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          unstagvar(i,j) = (stagvar(i,j) + stagvar(i-1,j) + stagvar(i,j-1) + stagvar(i-1,j-1)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over vertices with vmask = 1

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          sumvar = vmask(i-1,j)  *stagvar(i-1,j)   + vmask(i,j)  *stagvar(i,j)  &
                 + vmask(i-1,j-1)*stagvar(i-1,j-1) + vmask(i,j-1)*stagvar(i,j-1)  
          summask = real(vmask(i-1,j) + vmask(i,j) + vmask(i-1,j-1) + vmask(i,j-1), dp)
          if (summask > 0.d0) unstagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

    ! Fill in halo values
    call parallel_halo(unstagvar)

  end subroutine glissade_unstagger

!****************************************************************************

  subroutine glissade_centered_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        land_mask)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring points and is second-order accurate.
    !
    ! There are several choices for computing gradients at the ice margin:
    ! HO_MARGIN_GRADIENT_ALL = 0: All neighbor values are used to compute the gradient, including 
    !  values in ice-free cells.  This convention is used by Glide, but performs poorly for 
    !  ice shelves with a sudden drop in ice thickness and surface elevation at the margin.
    ! HO_MARGIN_GRADIENT_ICE_LAND = 1: Values in ice-covered and/or land cells are used to compute 
    !  the gradient, but values in ice-free ocean cells are ignored.  Where required values are 
    !  missing, the gradient is set to zero.  This reduces to option (0) for land-based problems 
    !  and (2) for ocean-based problems.
    ! HO_MARGIN_GRADIENT_ICE_ONLY = 2: Only values in ice-covered cells (i.e., cells with thck > thklim) 
    !  are used to compute gradients.  Where required values are missing, the gradient is set to zero.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    ! Since option (1) generally works well at both land and shelf boundaries, it is the default.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: summask, gradient_margin
    integer :: i, j

    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)


    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land mask to compute centered gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

          mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells
    endif

    ! Compute the gradients using info in cells with mask = 1

    do j = 1, ny-1
       do i = 1, nx-1

          summask = mask(i,j) + mask(i+1,j) + mask(i,j+1) + mask(i+1,j+1)

          if (summask == 4) then  ! use info in all four neighbor cells
             df_dx(i,j) = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1)) / (2.d0 * dx)
             df_dy(i,j) = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j)) / (2.d0 * dy)

          else  ! use info only in cells with mask = 1
                ! if info is not available, gradient component = 0

             ! df_dx
             if (mask(i,j)==1 .and. mask(i+1,j)==1) then
                df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
             elseif (mask(i,j+1)==1 .and. mask(i+1,j+1)==1) then
                df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
             endif

             ! df_dy
             if (mask(i,j)==1 .and. mask(i,j+1)==1) then
                df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
             elseif (mask(i+1,j)==1 .and. mask(i+1,j+1)==1) then
                df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
             endif

          endif

       enddo    ! i
    enddo       ! j

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Centered gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_centered_gradient

!****************************************************************************

  subroutine glissade_upstream_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        accuracy_flag_in,        &
                                        land_mask)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    !  compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient can be evaluated at two upstream points (for first-order accuracy) 
    !  or at four upstream points (for second-order accuracy).
    ! Note: Upstream is defined by the direction of higher surface elevation
    !  rather than the direction the flow is coming from (though these are
    !  usually the same).
    !
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       accuracy_flag_in         ! = 1 for 1st order, 2 for 2nd order

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: i, j
    real(dp) :: sum1, sum2
    integer :: gradient_margin, accuracy_flag, summask

    !   First-order upstream gradient at vertex(i,j) is based on two points out of f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !
    !   Second-order gradient is based on four points in the upstream direction

    if (present(accuracy_flag_in)) then
       accuracy_flag = accuracy_flag_in
    else
       accuracy_flag = 2   ! default to second-order
    endif

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land mask to compute upstream gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells

    endif

    if (accuracy_flag == 1) then   ! first-order accurate

       do j = 1, ny-1
          do i = 1, nx-1

             ! Compute gradient only if at least one neighbor is ice-covered
             summask = ice_mask(i,j) + ice_mask(i+1,j) + ice_mask(i,j+1) + ice_mask(i+1,j+1)
           
             if (summask > 0) then

                ! Compute df_dx by taking upstream gradient

                sum1 = f(i+1,j+1) + f(i,j+1)
                sum2 = f(i+1,j) + f(i,j)

                if (sum1 > sum2 .and. mask(i+1,j+1)==1 .and. mask(i,j+1)==1) then
                   df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                elseif (sum1 <= sum2 .and. mask(i+1,j)==1 .and. mask(i,j)==1) then
                   df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                endif

                ! Compute df_dy by taking upstream gradient
             
                sum1 = f(i+1,j+1) + f(i+1,j)
                sum2 = f(i,j+1) + f(i,j)
             
                if (sum1 > sum2 .and. mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then
                   df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                elseif (sum1 <= sum2 .and. mask(i,j+1)==1 .and. mask(i,j)==1) then
                   df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                else
                   df_dy(i,j) = 0.d0
                endif

             endif  ! summask > 0 (mask = 1 in at least one neighbor cell)

          enddo
       enddo

    else    ! second-order accurate

       do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             ! Compute gradient only if at least one neighbor is ice-covered
             summask = ice_mask(i,j) + ice_mask(i+1,j) + ice_mask(i,j+1) + ice_mask(i+1,j+1)
           
             if (summask > 0) then

                ! Compute df_dx by taking upstream gradient
             
                ! determine upstream direction

                sum1 = f(i+1,j+1) + f(i,j+1) + f(i+1,j+2) + f(i,j+2)
                sum2 = f(i+1,j) + f(i,j) + f(i+1,j-1) + f(i,j-1)

                if (sum1 > sum2) then

                   summask = mask(i+1,j+1) + mask(i,j+1) + mask(i+1,j+2) + mask(i,j+2)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dx(i,j) = (1.5d0 * (f(i+1,j+1) - f(i,j+1))     &
                                  - 0.5d0 * (f(i+1,j+2) - f(i,j+2))) / dx
                   elseif (mask(i+1,j+1)==1 .and. mask(i,j+1)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dx: i, j, summask =', i, j, summask
                      df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                   endif

                else  ! sum1 <= sum2

                   summask = mask(i+1,j) + mask(i,j) + mask(i+1,j-1) + mask(i,j-1)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dx(i,j) = (1.5d0 * (f(i+1,j)   - f(i,j))     &
                                  - 0.5d0 * (f(i+1,j-1) - f(i,j-1))) / dx
                   elseif (mask(i+1,j)==1 .and. mask(i,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dx: i, j, summask =', i, j, summask
                      df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                   endif
                   
                endif   ! sum1 > sum2

                ! Compute df_dy by taking upstream gradient

                ! determine upstream direction

                sum1 = f(i+1,j+1) + f(i+1,j) + f(i+2,j+1) + f(i+2,j)
                sum2 = f(i,j+1) + f(i,j) + f(i-1,j+1) + f(i-1,j)
             
                if (sum1 > sum2) then

                   summask = mask(i+1,j+1) + mask(i+1,j) + mask(i+2,j+1) + mask(i+2,j)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dy(i,j) = (1.5d0 * (f(i+1,j+1) - f(i+1,j))     &
                                  - 0.5d0 * (f(i+2,j+1) - f(i+2,j))) / dy
                   elseif (mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dy: i, j, summask =', i, j, summask
                      df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                   endif

                else   ! sum1 <= sum2

                   summask = mask(i,j+1) + mask(i,j) + mask(i-1,j+1) + mask(i-1,j)
                   
                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dy(i,j) = (1.5d0 * (f(i,j+1)   - f(i,j))     &
                                  - 0.5d0 * (f(i-1,j+1) - f(i-1,j))) / dy
                   elseif (mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dy: i, j, summask =', i, j, summask
                      df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                   endif

                endif   ! sum1 > sum2

             endif      ! summask > 0 (mask = 1 in at least one neighbor cell)

          enddo     ! i
       enddo        ! j

       ! fill in halo values
       call staggered_parallel_halo(df_dx)
       call staggered_parallel_halo(df_dy)

    endif   ! 1st or 2nd order accurate

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'upstream df_dx:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'upstream df_dy:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif

  end subroutine glissade_upstream_gradient

!****************************************************************************

  subroutine glissade_edge_gradient(nx,           ny,        &
                                    dx,           dy,        &
                                    f,                       &
                                    df_dx,        df_dy,     &
                                    gradient_margin_in,      &
                                    ice_mask,     land_mask)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) at cell edges (i.e., the C grid):
    ! df_dx at the midpoint of the east edge and df_dy at the midpoint of
    ! the north edge.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell edges

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, set gradient to zero
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, set gradient to zero

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask,     &          ! = 1 where ice is present, else = 0
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: gradient_margin
    integer :: i, j

    !   Gradient at east edge(i,j) is based on f(i:i+1,j)
    !   Gradient at north edge(i,j) is based on f(i,j:j+1)
    ! 
    !   |             |
    !   |   (i,j+1)   |
    !   |             |
    !   |             |
    !   ----df_dy------------------
    !   |             |  
    !   |             |
    !   |   (i,j)   df_dx   (i+1,j)
    !   |             |
    !   |             |
    !   |--------------

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask) .and. present(ice_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land and ice masks to compute edge gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       if (present(ice_mask)) then
          mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells
       else
          call write_log('Must pass in ice mask to compute edge gradient with gradient_margin = 2', GM_FATAL)
       endif

    endif

    ! Compute the gradients using info in cells with mask = 1

    do j = 1, ny-1
       do i = 1, nx-1

          ! df_dx

          if (mask(i,j)==1 .and. mask(i+1,j)==1) then
             df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
          endif

          ! df_dy

          if (mask(i,j)==1 .and. mask(i,j+1)==1) then
             df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
          endif

       enddo    ! i
    enddo       ! j

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Edge gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_edge_gradient

!----------------------------------------------------------------------------

  subroutine glissade_vertical_average(nx,         ny,        &
                                       nz,         sigma,     &
                                       mask,                  &
                                       var,        var_2d)

    !----------------------------------------------------------------
    ! Compute the vertical average of a given variable.
    ! Note: It is assumed that the variable is defined at layer midpoints,
    !       and hence has vertical dimension (nz-1).
    ! Note: This subroutine will work for variables on the staggered
    !       horizontal grid if stagthck is passed in place of thck.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,            &     ! horizontal grid dimensions
       nz                       ! number of vertical levels

    real(dp), dimension(nz), intent(in) ::    &
       sigma                    ! sigma vertical coordinate

    logical, dimension(nx, ny), intent(in) ::    &
       mask                     ! compute var_2d where mask = .true.

    real(dp), dimension(nz-1,nx, ny), intent(in) ::    &
       var                      ! 3D field to be averaged vertically

    real(dp), dimension(nx, ny), intent(out) ::    &
       var_2d                   ! 2D vertically averaged field

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j, k

    do j = 1, ny
       do i = 1, nx

          var_2d(i,j) = 0.d0

          if (mask(i,j)) then
             do k = 1, nz-1
                var_2d(i,j) = var_2d(i,j) + var(k,i,j) * (sigma(k+1) - sigma(k))
             enddo
          endif

       enddo
    enddo

  end subroutine glissade_vertical_average

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
