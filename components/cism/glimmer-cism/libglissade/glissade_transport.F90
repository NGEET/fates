!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_transport.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains drivers for incremental remapping and upwind ice transport.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This version was created from ice_transport_driver in CICE, revision 313, 6 Jan. 2011.
! The repository is here: http://oceans11.lanl.gov/svn/CICE

  module glissade_transport

    use glimmer_global, only: dp
    use glimmer_log
    use glissade_remap, only: glissade_horizontal_remap, make_remap_mask, puny
    use parallel 

    implicit none
    save
    private
    public :: glissade_transport_driver, glissade_check_cfl, ntracer

    logical, parameter ::  &
         prescribed_area = .false.  ! if true, prescribe the area fluxed across each edge

    integer :: &
         ntracer = 1   ! number of tracers to transport (just temperature if ntracer = 1)
                       ! can be set to a different value in glissade.F90

!=======================================================================

  contains

!=======================================================================
!
    subroutine glissade_transport_driver(dt,                   &
                                         dx,       dy,         &
                                         nx,       ny,         &
                                         nlyr,     sigma,      &
                                         ntracer,              &
                                         uvel,     vvel,       &
                                         thck,                 &
                                         acab,     bmlt,       &
                                         temp,     age,        &
                                         upwind_transport_in)


      ! This subroutine solves the transport equations for one timestep
      ! using the conservative remapping scheme developed by John Dukowicz
      ! and John Baumgardner and modified for sea ice by William
      ! Lipscomb and Elizabeth Hunke.
      !
      ! This scheme preserves monotonicity of ice area and tracers.  That is,
      ! it does not produce new extrema.  It is second-order accurate in space,
      ! except where gradients are limited to preserve monotonicity. 
      !
      ! Optionally, the remapping scheme can be replaced with a simple
      ! first-order upwind scheme.
      !
      ! author William H. Lipscomb, LANL
      !
      ! input/output arguments

      real(dp), intent(in) ::     &
         dt,                   &! time step (s)
         dx, dy                 ! gridcell dimensions (m)
                                ! (cells assumed to be rectangular)

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer                ! number of tracers

      real(dp), intent(in), dimension(nlyr+1) :: &
         sigma                  ! layer interfaces in sigma coordinates
                                ! top sfc = 0, bottom sfc = 1

      real(dp), intent(in), dimension(nlyr+1,nx-1,ny-1) :: &
         uvel, vvel             ! horizontal velocity components (m/s)
                                ! (defined at horiz cell corners, vertical interfaces)

      real(dp), intent(inout), dimension(nx,ny) :: &
         thck                   ! ice thickness (m), defined at horiz cell centers

      real(dp), intent(in), dimension(nx,ny) :: &
         acab,    &             ! surface mass balance (m/s)
                                ! (defined at horiz cell centers)
         bmlt                   ! basal melt rate (m/s)
                                ! positive for melting, negative for freeze-on
                                ! (defined at horiz cell centers)

      !NOTE: For the enthalpy scheme, the variable called temp is really the enthalpy.
      ! Note vertical dimension; includes surface and bed.
      ! (Surface and bed values are not transported.)
      real(dp), intent(inout), dimension(0:nlyr+1,nx,ny), optional :: &
         temp                   ! ice temperature (or enthalpy)
                                ! (defined at horiz cell centers, vertical layer midpts)
                                

      !TODO - Compute the ice age tracer; currently ice age = 0
      real(dp), intent(inout), dimension(nlyr,nx,ny), optional :: &
         age                    ! ice age

      logical, intent(in), optional ::  &
         upwind_transport_in    ! if true, do first-order upwind transport

      ! local variables

      integer ::     &
         i, j, k         ,&! cell indices
         ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
         nt                ! tracer index

      real(dp), dimension (nx,ny) ::     &
         thck_mask         ! = 1. if ice is present, = 0. otherwise

      real(dp), dimension (nx-1,ny-1) ::     &
         uvel_layer      ,&! uvel averaged to layer midpoint (m/s)
         vvel_layer        ! vvel averaged to layer midpoint (m/s)

      real(dp), dimension (nx,ny,nlyr) ::     &
         thck_layer        ! ice layer thickness (m)

      real(dp), dimension (nx,ny,ntracer,nlyr) ::     &
         tracer            ! tracer values

      real(dp), dimension (nx,ny) ::     &
         tsfc,            &! surface temperature, temp(0,:,:)
         tbed              ! basal temperature, temp(nlyr+1,:,:)

      integer ::     &
         icells            ! number of cells with ice

      integer, dimension(nx*ny) ::     &
         indxi, indxj      ! compressed i/j indices

      real(dp) ::         &
         sum_acab,        &! total accumulation/ablation
         sum_bmlt          ! total basal melting

      !-------------------------------------------------------------------
      ! If prescribed_area is true, the area of each departure region is
      !  computed in advance (e.g., by taking the divergence of the 
      !  velocity field and passed to locate_triangles.  The departure 
      !  regions are adjusted to obtain the desired area.
      ! If false, edgearea is computed in locate_triangles and passed out.
      !-------------------------------------------------------------------

      real(dp), dimension(nx,ny) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      logical, parameter ::     &
         conservation_check = .true. ! if true, check global conservation

      real(dp) ::     &
         msum_init,      &! initial global ice mass
         msum_final       ! final global ice mass

      real(dp), dimension(ntracer) ::     &
         mtsum_init,     &! initial global ice mass*tracer
         mtsum_final      ! final global ice mass*tracer

      real(dp) :: &
         melt_potential   ! total thickness (m) of additional ice that could be melted
                          ! by available acab/bmlt in columns that are completely melted

      logical ::     &
         errflag          ! true if energy is not conserved

      character(len=100) :: message

      real(dp), dimension (:,:,:), allocatable :: &
         worku            ! work array

      real(dp), dimension(nx,ny) ::      &
         uee, vnn            ! cell edge velocities for upwind transport

      logical ::     &
         upwind_transport    ! if true, do first-order upwind transport

      !-------------------------------------------------------------------
      ! Initialize
      !-------------------------------------------------------------------

      if (present(upwind_transport_in)) then
         upwind_transport = upwind_transport_in
      else
         upwind_transport = .false.
      endif

      errflag = .false.
      melt_potential = 0.d0

      !Note: (ilo,ihi) and (jlo,jhi) are the lower and upper bounds of the local domain
      ! (i.e., grid cells owned by this processor).

      ilo = nhalo + 1
      ihi = nx - nhalo
      jlo = nhalo + 1
      jhi = ny - nhalo

      !-------------------------------------------------------------------
      ! NOTE: Mass and tracer arrays (thck, temp, etc.) must be updated in 
      !       halo cells before this subroutine is called. 
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! Fill thickness and tracer arrays.
      ! Assume that temperature (if present) is tracer 1, and age (if present)
      !  is tracer 2.  Add more tracers as needed.
      ! If no tracers are present, then only the ice thickness is transported.
      !  In this case we define a dummy tracer array, since glissade_horizontal_remap
      !  requires that a tracer array be passed in.
      !-------------------------------------------------------------------

      do k = 1, nlyr

         thck_layer(:,:,k) = thck(:,:) * (sigma(k+1) - sigma(k))

         if (present(temp)) then
            tracer(:,:,1,k) = temp(k,:,:)
         else
            tracer(:,:,1,k) = 0.d0    ! dummy array
         endif

         if (present(age) .and. ntracer >= 2) then
            tracer(:,:,2,k) = age(k,:,:)
         endif

         ! Other tracer fields could be added here

      enddo

    !-------------------------------------------------------------------
    ! Set surface and basal temperatures.
    ! Note: If the enthalpy is passed in, then tsfc and tbed are enthalpy values.
    !-------------------------------------------------------------------

      if (present(temp)) then
         tsfc(:,:) = temp(0,:,:)
         tbed(:,:) = temp(nlyr+1,:,:)
      else
         tsfc(:,:) = 0.d0
         tbed(:,:) = 0.d0
      endif

    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities (optional)
    !-------------------------------------------------------------------

      if (conservation_check) then

         call sum_mass_and_tracers(nx,                ny,              &
                                   nlyr,              ntracer,         &
                                   nhalo,                              &
                                   thck_layer(:,:,:), msum_init,       &
                                   tracer(:,:,:,:),   mtsum_init(:))
      endif

      if (upwind_transport) then

         allocate (worku(nx,ny,0:ntracer))

         do k = 1, nlyr

            ! Average corner velocities at layer interfaces to 
            ! edge velocities at layer midpoints.
      
            do j = jlo, jhi
            do i = ilo-1, ihi   ! include west edge of local domain
               uee(i,j) = 0.25d0 * (uvel(k,  i,j) + uvel(k,  i,j-1)    &
                                  + uvel(k+1,i,j) + uvel(k+1,i,j-1))
            enddo
            enddo


            do j = jlo-1, jhi   ! include south edge of local domain
            do i = ilo, ihi
               vnn(i,j) = 0.25d0 * (vvel(k,  i,j) + vvel(k,  i-1,j)    &
                                  + vvel(k+1,i,j) + vvel(k+1,i-1,j))
            enddo
            enddo

            ! Fill work array for transport

            worku(:,:,0) = thck_layer(:,:,k)
            do nt = 1, ntracer
               worku(:,:,nt) = thck_layer(:,:,k) * tracer(:,:,nt,k)
            enddo

            !-----------------------------------------------------------------
            ! Upwind transport
            !-----------------------------------------------------------------
 
            do nt = 0, ntracer
               call upwind_field (nx,             ny,                  &
                                  ilo, ihi,       jlo, jhi,            &
                                  dx,             dy,                  &
                                  dt,             worku(:,:,nt),       &
                                  uee(:,:),       vnn    (:,:))
            enddo   ! ntracer

            ! Recompute thickness and tracers

            thck_layer(:,:,k) = worku(:,:,0)
            do nt = 1, ntracer
               do j = jlo, jhi
               do i = ilo, ihi
                  if (thck_layer(i,j,k) > puny) then
                     tracer(i,j,nt,k) = worku(i,j,nt) / thck_layer(i,j,k)
                  else
                     tracer(i,j,nt,k) = 0.d0
                  endif
               enddo   ! i
               enddo   ! j
            enddo      ! ntracer

         enddo         ! nlyr

         deallocate (worku)

      else    ! remapping transport

      !-------------------------------------------------------------------
      ! Define a mask: = 1 where ice is present (thck > 0), = 0 otherwise         
      ! The mask is used to prevent tracer values in cells without ice from
      !  being used to compute tracer gradients.
      !-------------------------------------------------------------------

         call make_remap_mask (nx,           ny,                 &
                               ilo, ihi,     jlo, jhi,           &
                               nhalo,        icells,             &
                               indxi(:),     indxj(:),           &
                               thck(:,:),    thck_mask(:,:))

      !-------------------------------------------------------------------
      ! Remap ice thickness and tracers; loop over layers
      !-------------------------------------------------------------------

         do k = 1, nlyr

            ! Average velocities to the midpoint of the layer

            uvel_layer(:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
            vvel_layer(:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))

            edgearea_e(:,:) = 0.d0
            edgearea_n(:,:) = 0.d0

            ! If prescribed_area is true, compute edgearea by taking the divergence
            !  of the velocity field.

            if (prescribed_area) then

               do j = jlo, jhi
               do i = ilo-1, ihi
                  edgearea_e(i,j) = (uvel_layer(i,j) + uvel_layer(i,j-1))     &
                                    * 0.5d0 * dy * dt
               enddo
               enddo
 
               do j = jlo-1, jhi
               do i = ilo, ihi
                  edgearea_n(i,j) = (vvel_layer(i,j) + vvel_layer(i-1,j))     &
                                   * 0.5d0 * dx * dt
               enddo
               enddo

            endif

         !-------------------------------------------------------------------
         ! Main remapping routine: Step ice thickness and tracers forward in time.
         !-------------------------------------------------------------------

            call glissade_horizontal_remap (dt,                                  &
                                            dx,                dy,               &
                                            nx,                ny,               &
                                            ntracer,           nhalo,            &
                                            thck_mask(:,:),    icells,           &
                                            indxi(:),          indxj(:),         &
                                            uvel_layer(:,:),   vvel_layer(:,:),  &
                                            thck_layer(:,:,k), tracer(:,:,:,k),  &
                                            edgearea_e(:,:),   edgearea_n(:,:))

         enddo    ! nlyr

      endif       ! remapping v. upwind transport

      !-------------------------------------------------------------------
      ! Check that mass and mass*tracers are exactly conserved by transport.
      ! Note: Conservation errors will occur if the global domain is open
      !       and ice has left the domain. So depending on the application,
      !       there may or may not be a problem when ice is not conserved.
      !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute new values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         call sum_mass_and_tracers(nx,                ny,              &
                                   nlyr,              ntracer,         &
                                   nhalo,                              &
                                   thck_layer(:,:,:), msum_final,      &
                                   tracer(:,:,:,:),   mtsum_final(:))

         ! Check conservation

         if (main_task) then
            call global_conservation (msum_init,     msum_final,      &
                                      errflag,       melt_potential,  &
                                      ntracer,                        &
                                      mtsum_init,    mtsum_final)
            if (errflag) then
               write(message,*) 'WARNING: Conservation error in glissade_horizontal_remap'
!               call write_log(message,GM_FATAL)      ! uncomment if conservation errors should never happen
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
               write(message,*) 'May be OK if global domain is open'
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif

         endif   ! main_task

      endif      ! conservation_check

      !-------------------------------------------------------------------
      ! Add the mass balance at the surface and bed.
      ! The reason to do this now rather than at the beginning of the
      !  subroutine is that the velocity was computed for the old geometry,
      !  before the addition or loss of new mass at the surface and bed.  
      ! TODO: Rethink this ordering if we move to implicit or semi-implicit
      !       timestepping, where the velocity depends on the new geometry.
      !
      ! We assume here that new ice arrives at the surface with the same 
      !  temperature as the surface.
      ! TODO: Make sure this assumption is consistent with energy
      !       conservation for coupled simulations.
      ! TODO: Pass the melt potential back to the climate model as a heat flux?
      !-------------------------------------------------------------------

      call glissade_add_smb(nx,       ny,         &
                            nlyr,     ntracer,    &
                            nhalo,    dt,         &
                            thck_layer(:,:,:),    &
                            tracer(:,:,:,:),      &
                            acab(:,:),            &
                            tsfc(:,:),            &
                            bmlt(:,:),            &
                            tbed(:,:),            &
                            melt_potential )

      !-------------------------------------------------------------------
      ! Next conservation check: Check that mass is conserved, allowing
      !  for mass gain/loss due to acab/bmlt and for any unused melt potential.
      !
      ! NOTE: There is no tracer conservation check here, because there is no
      !       easy way to correct initial mass*tracer values for acab and bmlt.
      !
      !-------------------------------------------------------------------

      if (conservation_check) then

         ! Correct initial global mass for acab and bmlt

         sum_acab = sum(acab(1+lhalo:nx-uhalo,1+lhalo:ny-uhalo))
         sum_bmlt = sum(bmlt(1+lhalo:nx-uhalo,1+lhalo:ny-uhalo))
         call global_sum(sum_acab)
         call global_sum(sum_bmlt)

         msum_init = msum_init + (sum_acab - sum_bmlt)*dt

         ! Compute new global mass and mass*tracer

         call sum_mass_and_tracers(nx,                ny,              &
                                   nlyr,              ntracer,         &
                                   nhalo,                              &
                                   thck_layer(:,:,:), msum_final,      &
                                   tracer(:,:,:,:),   mtsum_final(:))

         ! Check mass conservation

         if (main_task) then
            call global_conservation (msum_init,     msum_final,      &
                                      errflag,       melt_potential)

            if (errflag) then
               write(message,*) 'WARNING: Conservation error in glissade_add_smb'
!               call write_log(message,GM_FATAL)      ! uncomment if conservation errors should never happen
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif

         endif   ! main_task

      endif      ! conservation_check

      !-------------------------------------------------------------------
      ! Interpolate tracers back to sigma coordinates 
      !-------------------------------------------------------------------

      call glissade_vertical_remap(nx,                ny,      &
                                   nlyr,              ntracer, &
                                   sigma(:),                   &
                                   thck_layer(:,:,:),          &
                                   tracer(:,:,:,:) )

      !-------------------------------------------------------------------
      ! Final conservation check: Check that mass and mass*tracers are exactly 
      !                           conserved by vertical remapping..
      !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute new values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         msum_init = msum_final   ! msum_final computed above after glissade_add_smb
         mtsum_init(:) = mtsum_final(:)   ! mtsum_final computed above after glissade_add_smb
 
         call sum_mass_and_tracers(nx,                ny,              &
                                   nlyr,              ntracer,         &
                                   nhalo,                              &
                                   thck_layer(:,:,:), msum_final,      &
                                   tracer(:,:,:,:),   mtsum_final(:))

         ! Check conservation

         if (main_task) then
            call global_conservation (msum_init,     msum_final,  &
                                      errflag,       0.d0,        &  ! ignore melt potential for this check
                                      ntracer,                    &
                                      mtsum_init,    mtsum_final)
            if (errflag) then
               write(message,*) 'WARNING: Conservation error in glissade_vertical_remap'
!               call write_log(message,GM_FATAL)      ! uncomment if conservation errors should never happen
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif

         endif   ! main_task

      endif      ! conservation_check

      !-------------------------------------------------------------------
      ! Halo updates for thickness and tracer arrays
      !
      ! Note: Cannot pass the full 3D array to parallel_halo, because that
      !       subroutine assumes that k is the first rather than third index. 
      !-------------------------------------------------------------------

      do k = 1, nlyr

         call parallel_halo(thck_layer(:,:,k))

         do nt = 1, ntracer
            call parallel_halo(tracer(:,:,nt,k))
         enddo

      enddo

      !-------------------------------------------------------------------
      ! Recompute thickness, temperature and other tracers.
      !-------------------------------------------------------------------

      thck(:,:) = 0.d0

      do k = 1, nlyr

         thck(:,:) = thck(:,:) + thck_layer(:,:,k)

         if (present(temp)) temp(k,:,:) = tracer(:,:,1,k)
         if (present(age) .and. ntracer >= 2) age(k,:,:) = tracer(:,:,2,k)

         ! Could add more tracer fields here

      enddo

    end subroutine glissade_transport_driver

!=======================================================================

    subroutine glissade_check_cfl(ewn,     nsn,      nlyr,          & 
                                  dew,     dns,      sigma,         &
                                  stagthk, dusrfdew, dusrfdns,      &
                                  uvel,    vvel,     deltat,        &
                                  allowable_dt_adv,  allowable_dt_diff)

      ! Calculate maximum allowable time step based on both 
      ! advective and diffusive CFL limits.
      !
      ! author Matt Hoffman, LANL, March 2014
      !
      ! input/output arguments

      integer, intent(in) ::     &
         ewn, nsn    ! number of cells in the x, y dimensions

      integer, intent(in) ::     &
         nlyr        ! number of vertical layers (layer centers)

      real(dp), intent(in) :: &
         dew, dns    ! grid spacing in x, y (not assumed to be equal here), dimensional m

      real(dp), dimension(:), intent(in) :: &
         sigma       ! vertical coordinate spacing

      real(dp), dimension(:,:), intent(in) :: &
         stagthk     ! thickness on the staggered grid, dimensional m

      real(dp), dimension(:,:), intent(in) :: &
         dusrfdew, dusrfdns    ! slope in x,y directions on the staggered grid, dimensionless m/m

      real(dp), dimension(:,:,:), intent(in) :: &
         uvel, vvel  ! 3-d x,y velocity components on the staggered grid, dimensional m/yr

      real(dp), intent(in) :: &
         deltat      ! model deltat (yrs)

      real(dp), intent(out) :: &
         allowable_dt_adv     ! maximum allowable dt (yrs) based on advective CFL 

      real(dp), intent(out) :: &
         allowable_dt_diff    ! maximum allowable dt (yrs) based on diffusive CFL 

      ! Local variables
      integer :: k
      integer :: xs, xe, ys, ye  ! start and end indices for locally owned cells on the staggered grid in the x and y directions
      real(dp), dimension(nlyr, ewn-1, nsn-1) :: uvel_layer, vvel_layer  ! velocities at layer midpoints, stag. grid
      real(dp), dimension(nlyr, ewn-1, nsn-1) :: flux_layer_ew, flux_layer_ns  ! flux for each layer, stag. grid
      real(dp), dimension(ewn-1, nsn-1) :: flux_ew, flux_ns  ! flux for entire thickness, stag. grid

      real(dp) :: maxuvel, maxvvel, maxvel ! maximum velocity in either direction and in both
      real(dp) :: allowable_dt_diff_here ! temporary calculation at each cell of allowable_dt_diff
      integer :: i, j
      real(dp) :: slopemag  ! the magnitude of the surface slope
      real(dp) :: slopedirx, slopediry  ! the unit vector of the slope direction
      real(dp) :: flux_downslope  ! The component of the flux in the downslope direction
      integer :: ierr ! flag for CFL violation
      integer :: procnum  ! processor on which minimum allowable time step occurs
      integer, dimension(3) :: indices_adv  ! z,x,y indices (stag. grid) of where the min. allow. time step occurs for the advective CFL
      integer, dimension(2) :: indices_diff  ! x and y indices (stag. grid) of where the min. allow. time step occurs for  the  diffusive CFL
      character(len=12)  :: dt_string, xpos_string, ypos_string
      character(len=300) :: message
      ierr = 0

      ! Setup some mesh information - start and end indices for locally owned cells on the staggered grid in the x and y directions
      xs = 1+staggered_lhalo
      xe = ewn - 1 - staggered_uhalo
      ys = 1+staggered_lhalo
      ye = nsn - 1 - staggered_uhalo

      ! ------------------------------------------------------------------------
      ! Advective CFL
      ! TODO use depth-averaged velocity or layer-by-layer (top layer only?), or something else (BB09)?
      ! For now check all layers

      ! Calculate depth-averaged flux and velocity on the B-grid.
      ! The IR code basically uses a B-grid, the FO-Upwind method uses a C-grid.
      ! The B-grid calculation should be more conservative because that is where the
      ! velocities are calculated so there will be no averaging.
      ! (Also, IR is the primary advection method, so make this check most appropriate for that.)

      do k = 1, nlyr
         ! Average velocities to the midpoint of the layer
         uvel_layer(k,:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
         vvel_layer(k,:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))
         ! calculate flux components for this layer
         flux_layer_ew(k,:,:) = uvel_layer(k,:,:) * stagthk(:,:) * (sigma(k+1) - sigma(k))
         flux_layer_ns(k,:,:) = vvel_layer(k,:,:) * stagthk(:,:) * (sigma(k+1) - sigma(k))
      enddo
      flux_ew(:,:) = sum(flux_layer_ew, 1)
      flux_ns(:,:) = sum(flux_layer_ns, 1)

      ! Advective CFL calculation - using all layers. Check locally owned cells only!
      maxuvel = maxval(abs(uvel_layer(:,xs:xe,ys:ye)))
      maxvvel = maxval(abs(vvel_layer(:,xs:xe,ys:ye)))
      ! Determine in which direction the max velocity is - Assuming dx=dy here!
      if (maxuvel > maxvvel) then
!         print *, 'max vel is in uvel'
         maxvel = maxuvel
         indices_adv = maxloc(abs(uvel_layer(:,xs:xe,ys:ye)))
      else
!         print *, 'max vel is in vvel'
         maxvel = maxvvel
         indices_adv = maxloc(abs(vvel_layer(:,xs:xe,ys:ye)))
      endif
      indices_adv(2:3) = indices_adv(2:3) + staggered_lhalo  ! want the i,j coordinates WITH the halo present - we got indices into the slice of owned cells
      ! Finally, determine maximum allowable time step based on advectice CFL condition.
      allowable_dt_adv = dew / (maxvel + 1.0d-20)

      ! ------------------------------------------------------------------------
      ! Diffusive CFL
      ! Estimate diffusivity using the relation that the 2-d flux Q=-D grad h and Q=UH, 
      ! where h is surface elevation, D is diffusivity, U is 2-d velocity vector, and H is thickness
      ! Solving for D = UH/-grad h

      allowable_dt_diff = 1.0d20  ! start with a huge value
      indices_diff(:) = 1 ! Initialize these to something, on the off-chance they never get set... (e.g., no ice on this processor)
      ! Loop over locally-owned cells only!
      do j = ys, ye
         do i = xs, xe
            if (stagthk(i,j) > 0.0d0) then  ! don't bother doing all this for non-ice cells
                ! Find downslope vector
                slopemag = dsqrt(dusrfdew(i,j)**2 + dusrfdns(i,j)**2 + 1.0d-20)
                slopedirx = dusrfdew(i,j) / slopemag
                slopediry = dusrfdns(i,j) / slopemag
                ! Estimate flux in the downslope direction (Flux /dot -slopedir)
                flux_downslope = flux_ew(i,j) * (-1.0d0) * slopedirx + flux_ns(i,j) * (-1.0d0) * slopediry  ! TODO check signs here - they seem ok
                !!! Estimate diffusivity in the downslope direction only
                !!diffu = flux_downslope / slopemag
                !!allowable_dt_diff = 0.5d0 * dew**2 / (diffu + 1.0e-20)  ! Note: assuming diffu is isotropic here.
                ! DCFL: dt = 0.5 * dx**2 / D = 0.5 * dx**2 * slopemag / flux_downslope
                allowable_dt_diff_here = 0.5d0 * dew**2 * slopemag / (flux_downslope + 1.0e-20)  ! Note: assuming diffu is isotropic here.  assuming dx=dy
                if (allowable_dt_diff_here < 0.0d0) allowable_dt_diff_here = 1.0d20 ! ignore negative dt's (upgradient flow due to membrane stresses)
                if (allowable_dt_diff_here < allowable_dt_diff) then
                   allowable_dt_diff = allowable_dt_diff_here
                   indices_diff(1) = i
                   indices_diff(2) = j
                endif
            endif
         enddo
      enddo

      ! Determine location limiting the DCFL
!      print *, 'diffu dt', allowable_dt_diff, indices_diff(1), indices_diff(2)

      ! Optional print of local limiting dt on each procesor
      !print *,'LOCAL ADV DT, POSITION', allowable_dt_adv, indices_adv(2), indices_adv(3)
      !print *,'LOCAL DIFF DT, POSITION', allowable_dt_diff, indices_diff(1), indices_diff(2)

      ! ------------------------------------------------------------------------
      ! Now check for errors

      ! Perform global reduce for advective time step and determine where in the domain it occurs
      call parallel_reduce_minloc(xin=allowable_dt_adv, xout=allowable_dt_adv, xprocout=procnum)

      if (deltat > allowable_dt_adv) then
          ierr = 1  ! Advective CFL violation is a fatal error

          ! Get position of the limiting location - do this only if an error message is needed to avoid 2 MPI comms
          call parallel_globalindex(indices_adv(2), indices_adv(3), indices_adv(2), indices_adv(3))  
          ! Note: This subroutine assumes the scalar grid, but should work fine for the stag grid too
          ! indices_adv now has i,j on the global grid for this proc's location
          call broadcast(indices_adv(2), proc=procnum)
          call broadcast(indices_adv(3), proc=procnum)
          ! indices_adv now has i,j on the global grid for the limiting proc's location

          write(dt_string,'(f12.5)') allowable_dt_adv
          write(xpos_string,'(i12)') indices_adv(2)
          write(ypos_string,'(i12)') indices_adv(3)
          write(message,*) 'Advective CFL violation!  Maximum allowable time step for advective CFL condition is ' &
               // trim(adjustl(dt_string)) // ' yr, limited by global position i=' // trim(adjustl(xpos_string)) // ' j=' //trim(adjustl(ypos_string))
          ! Write a warning first before throwing a fatal error so we can also check the diffusive CFL before aborting
          call write_log(trim(message),GM_WARNING)      
      endif

      ! Perform global reduce for diffusive time step and determine where in the domain it occurs
      call parallel_reduce_minloc(xin=allowable_dt_diff, xout=allowable_dt_diff, xprocout=procnum)

      if (deltat > allowable_dt_diff) then
          ! Get position of the limiting location - do this only if an error message is needed to avoid 2 MPI comms
          call parallel_globalindex(indices_diff(1), indices_diff(2), indices_diff(1), indices_diff(2))  
          ! Note: this subroutine assumes the scalar grid, but should work fine for the stag grid too
          ! indices_diff now has i,j on the global grid for this proc's location
          call broadcast(indices_diff(1), proc=procnum)
          call broadcast(indices_diff(2), proc=procnum)
          ! indices_diff now has i,j on the global grid for the limiting proc's location

          write(dt_string,'(f12.5)') allowable_dt_diff
          write(xpos_string,'(i12)') indices_diff(1)
          write(ypos_string,'(i12)') indices_diff(2)
          write(message,*) 'Diffusive CFL violation!  Maximum allowable time step for diffusive CFL condition is ' &
               // trim(adjustl(dt_string)) // ' yr, limited by global position i=' // trim(adjustl(xpos_string)) // ' j=' //trim(adjustl(ypos_string))
          ! Diffusive CFL violation is just a warning (because it may be overly restrictive as currently formulated)
          call write_log(trim(message),GM_WARNING)    
          write(message,*) '(Note the currently implemented diffusive CFL calculation may be overly restrictive for higher-order dycores.)'
          call write_log(trim(message))
      endif

      ! TODO enable this fatal error after more testing!
      ! Now that we have checked both, throw a fatal error for an ACFL violation
      !if (ierr == 1) then
      !   call write_log('Advective CFL violation is a fatal error.  See log for details.', GM_FATAL)
      !endif

    end subroutine glissade_check_cfl

!=======================================================================

    subroutine sum_mass_and_tracers(nx,         ny,        &
                                    nlyr,       ntracer,   &
                                    nhalo,                 &
                                    thck_layer, msum,      &
                                    tracer,     mtsum)

      ! Compute values of globally conserved quantities.,
      ! Assume gridcells of equal area (dx*dy), ice of uniform density.
      
      ! Input/output arguments

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer,              &! number of tracers
         nhalo                  ! number of halo rows

      real(dp), dimension (nx,ny,nlyr), intent(in) ::     &
         thck_layer             ! ice layer thickness

      real(dp), intent(out) ::   &
         msum                   ! total mass (actually thickness, measured in m)

      real(dp), dimension (nx,ny,ntracer,nlyr), intent(in), optional ::   &
         tracer                 ! tracer values

      real(dp), dimension(ntracer), intent(out), optional ::   &
         mtsum                  ! total mass*tracer

      ! Local arguments

      integer :: i, j, nt

      msum  = 0.d0
      if (present(mtsum)) mtsum(:) = 0.d0

      do j = 1+nhalo, ny-nhalo
         do i = 1+nhalo, nx-nhalo
               
            ! accumulate ice mass and mass*tracers
            ! (actually, accumulate thickness, assuming rhoi*dx*dy is the same for each cell)

            msum = msum + sum(thck_layer(i,j,:))

            if (present(mtsum)) then
               do nt = 1, ntracer
                  mtsum(nt) =  mtsum(nt) + sum(tracer(i,j,nt,:)*thck_layer(i,j,:))
               enddo
            endif

         enddo    ! i
      enddo       ! j

      call global_sum(msum)
      if (present(mtsum)) call global_sum(mtsum)

    end subroutine sum_mass_and_tracers

!=======================================================================
!
    subroutine global_conservation (msum_init,  msum_final,         &
                                    errflag,    melt_potential_in,  &
                                    ntracer,                        &
                                    mtsum_init, mtsum_final)
      !
      ! Check whether values of conserved quantities have changed.
      ! An error probably means that ghost cells are treated incorrectly.
      !
      ! author William H. Lipscomb, LANL
      !
      ! input/output arguments

      real(dp), intent(in) ::     &
         msum_init   ,&! initial global ice mass 
         msum_final    ! final global ice mass

      logical, intent(out) ::     &
         errflag       ! true if there is a conservation error

      real(dp), intent(in), optional :: &
         melt_potential_in   ! total thickness (m) of additional ice that could be melted
                             ! by available acab/bmlt in columns that are completely melted

      integer, intent(in), optional ::  &
         ntracer       ! number of tracers

      real(dp), dimension(:), intent(in), optional :: &
         mtsum_init  ,&! initial global ice mass*tracer
         mtsum_final   ! final global ice mass*tracer

      character(len=100) :: message

      integer ::      &
         nt            ! tracer index

      real(dp) ::         &
         melt_potential,  &! melt_potential_in (if present), else = 0 
         diff              ! difference between initial and final values

      if (present(melt_potential_in)) then
         melt_potential = melt_potential_in
      else
         melt_potential = 0.d0
      endif

      errflag = .false.

      if (msum_init > puny) then
         diff = (msum_final - melt_potential) - msum_init
         if (abs(diff/msum_init) > puny) then
            errflag = .true.
            write (message,*) 'glissade_transport: ice mass conservation error'
            call write_log(message)
            write (message,*) 'Initial global mass =', msum_init
            call write_log(message)
!            write (message,*) 'Final global mass =', msum_final
!            call write_log(message)
!            write (message,*) 'Melt potential =', melt_potential
!            call write_log(message)
            write (message,*) 'Final global mass (adjusted for melt potential) =', msum_final - melt_potential
            call write_log(message)
            write (message,*) 'Fractional error =', abs(diff)/msum_init
            call write_log(message)
         endif
      endif

      if (present(mtsum_init)) then
         do nt = 1, ntracer
            if (abs(mtsum_init(nt)) > puny) then
               diff = mtsum_final(nt) - mtsum_init(nt)
               if (abs(diff/mtsum_init(nt)) > puny) then
                  errflag = .true.
                  write (message,*) 'glissade_transport: mass*tracer conservation error'
                  call write_log(message)
                  write (message,*) 'tracer index =', nt
                  call write_log(message)
                  write (message,*) 'Initial global mass*tracer =', mtsum_init(nt)
                  call write_log(message)
                  write (message,*) 'Final global mass*tracer =', mtsum_final(nt)
                  call write_log(message)
                  write (message,*) 'Fractional difference =', abs(diff)/mtsum_init(nt)
                  call write_log(message)
               endif
            endif
         enddo
      endif                     ! present(mtsum_init)

    end subroutine global_conservation

!----------------------------------------------------------------------

    subroutine glissade_add_smb(nx,         ny,         &
                                nlyr,       ntracer,    &
                                nhalo,      dt,         &
                                thck_layer, tracer,     &
                                acab,       tsfc,       &
                                bmlt,       tbed,       &
                                melt_potential)

      ! Adjust the layer thickness based on the surface and basal mass balance

      ! Input/output arguments

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer,              &! number of tracers
         nhalo                  ! number of halo rows

      real(dp), intent(in) ::   &
         dt                     ! time step (s)

      real(dp), dimension (nx,ny,nlyr), intent(inout) ::     &
         thck_layer             ! ice layer thickness

      real(dp), dimension (nx,ny,ntracer,nlyr), intent(inout) ::     &
         tracer                 ! tracer values

      real(dp), intent(in), dimension(nx,ny) :: &
         acab                   ! surface mass balance (m/s)

      real(dp), intent(in), dimension(nx,ny) :: &
         bmlt                   ! basal melt rate (m/s)
                                ! > 0 for melting, < 0 for freeze-on

      real(dp), intent(out) :: &
         melt_potential   ! total thickness (m) of additional ice that could be melted
                          ! by available acab/bmlt in columns that are completely melted

      ! Note: If enthalpy is passed in, then tsfc and tbed are enthalpies.

      real(dp), intent(in), dimension(nx,ny) :: &
         tsfc                   ! surface temperature (deg C)

      real(dp), intent(in), dimension(nx,ny) :: &
         tbed                   ! basal temperature (deg C)

      ! Local variables

      real(dp), dimension(nx,ny,ntracer,nlyr) ::  &
         thck_tracer       ! thck_layer * tracer

      real(dp) :: sfc_accum, sfc_ablat  ! surface accumulation/ablation, from acab
      real(dp) :: bed_accum, bed_ablat  ! bed accumulation/ablation, from bmlt

      integer :: i, j, k

      character(len=100) :: message

      melt_potential = 0.d0

      do j = 1+nhalo, ny-nhalo
         do i = 1+nhalo, nx-nhalo

            ! initialize accumulation/ablation terms
            sfc_accum = 0.d0
            sfc_ablat = 0.d0
            bed_accum = 0.d0
            bed_ablat = 0.d0
            
            ! Add surface accumulation/ablation to ice thickness
            ! Also modify tracers conservatively.
            ! Assume tracer(:,:,1,:) is temperature and is always present
            ! Asumme tracer(:,:,2,:) is ice age and is optionally present

            if (acab(i,j) > 0.d0) then       ! accumulation, added to layer 1

               sfc_accum = acab(i,j)*dt

               ! temperature of top layer
               if (ntracer >= 1) then
                  thck_tracer(i,j,1,1) = thck_layer(i,j,1) * tracer(i,j,1,1)  &
                                       + sfc_accum * tsfc(i,j)
               endif

               ! ice age (= 0 for new accumulation)
               if (ntracer >= 2) then  
                  thck_tracer(i,j,2,1) = thck_layer(i,j,1) * tracer(i,j,2,1) 
                                     ! + sfc_accum * 0.d0

                  !TODO - Add other tracers here, as needed

               endif

               ! new top layer thickess
               thck_layer(i,j,1) = thck_layer(i,j,1) + sfc_accum

               ! new tracer values in top layer
               tracer(i,j,:,1) = thck_tracer(i,j,:,1) / thck_layer(i,j,1)

            elseif (acab(i,j) < 0.d0) then   ! ablation in one or mor layers            

               ! reduce ice thickness (tracer values will not change)

               sfc_ablat = -acab(i,j)*dt   ! positive by definition

               do k = 1, nlyr
                  if (sfc_ablat > thck_layer(i,j,k)) then
                     sfc_ablat = sfc_ablat - thck_layer(i,j,k)
                     thck_layer(i,j,k) = 0.d0
                     tracer(i,j,:,k) = 0.d0
                  else
                     thck_layer(i,j,k) = thck_layer(i,j,k) - sfc_ablat
                     sfc_ablat = 0.d0
                     exit
                  endif
               enddo

               if (sfc_ablat > 0.d0) then
                  melt_potential = melt_potential + sfc_ablat
               endif

            endif  ! acab > 0

            !TODO - Figure out how to handle excess energy given by melt_potential.
            !       Include in the heat flux passed back to CLM?

            ! Note: It is theoretically possible that we could have residual energy remaining for surface
            ! ablation while ice is freezing on at the bed, in which case the surface ablation should
            ! be subtracted from the bed accumulation.  We ignore this possibility for now.

            if (bmlt(i,j) < 0.d0) then       ! freeze-on, added to lowest layer

               bed_accum = -bmlt(i,j)*dt

               ! temperature of bottom layer
               if (ntracer >= 1) then
                  thck_tracer(i,j,1,nlyr) = thck_layer(i,j,nlyr) * tracer(i,j,1,nlyr)  &
                                          + bed_accum * tbed(i,j)
               endif

               ! ice age (= 0 for new accumulation)
               if (ntracer >= 2) then  
                  thck_tracer(i,j,2,nlyr) = thck_layer(i,j,nlyr) * tracer(i,j,2,nlyr) 
                                        ! + bed_accum * 0.d0

                  !TODO - Add other tracers here, as needed
               endif

               ! new bottom layer thickess
               thck_layer(i,j,nlyr) = thck_layer(i,j,nlyr) + bed_accum

               ! new tracer values in bottom layer
               tracer(i,j,:,nlyr) = thck_tracer(i,j,:,nlyr) / thck_layer(i,j,nlyr)

            elseif (bmlt(i,j) > 0.d0) then   ! basal melting in one or more layers            

               ! reduce ice thickness (tracer values will not change)

               bed_ablat = bmlt(i,j)*dt   ! positive by definition

               do k = nlyr, 1, -1
                  if (bed_ablat > thck_layer(i,j,k)) then
                     bed_ablat = bed_ablat - thck_layer(i,j,k)
                     thck_layer(i,j,k) = 0.d0
                     tracer(i,j,:,k) = 0.d0
                  else
                     thck_layer(i,j,k) = thck_layer(i,j,k) - bed_ablat
                     bed_ablat = 0.d0
                     exit
                  endif
               enddo
  
               if (bed_ablat > 0.d0) then
                  melt_potential = melt_potential + bed_ablat
               endif

            endif  ! bmlt < 0

         enddo   ! i
      enddo      ! j

      call global_sum(melt_potential)

    end subroutine glissade_add_smb

!----------------------------------------------------------------------

    subroutine glissade_vertical_remap(nx,       ny,        &
                                       nlyr,     ntracer,   &
                                       sigma,    hlyr,      &
                                       trcr)
 
    ! Conservative remapping of tracer fields from one set of vertical 
    ! coordinates to another.  The remapping is first-order accurate.
    !
    ! Note: The cost of this subroutine scales as nlyr; 
    !       a previous version scaled as nlyr^2.
    !
    ! TODO - Add a 2nd-order accurate vertical remapping scheme?
    !
    ! Author: William Lipscomb, LANL

    implicit none
 
    ! in-out arguments
 
    integer, intent(in) ::  &
         nx, ny,     &! number of cells in EW and NS directions
         nlyr,       &! number of vertical layers
         ntracer      ! number of tracer fields

    real(dp), dimension (nx, ny, nlyr), intent(inout) ::  &
         hlyr         ! layer thickness

    real(dp), dimension (nlyr+1), intent(in) ::  &
         sigma        ! sigma vertical coordinate (at layer interfaces)

    real(dp), dimension (nx, ny, ntracer, nlyr), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! tracer(k) = value at midpoint of layer k
 
    ! local variables
 
    integer :: i, j, k, k1, k2, nt
 
    real(dp), dimension (nlyr+1) ::  &
         z1,        &! layer interfaces in old coordinate system
                     ! z1(1) = 0. = value at top surface
                     ! z1(k) = value at top of layer k
                     ! z1(nlyr+1) = value at bottom surface (= 1 in sigma coordinates)
         z2          ! layer interfaces in new coordinate system

    real(dp) ::        &
         thck,      &! total thickness
         rthck       ! reciprocal of total thickness
 
    real(dp), dimension(ntracer,nlyr) ::       &
         htsum        ! sum of thickness*tracer in a layer         

    real(dp) :: zlo, zhi, hovlp

    do j = 1, ny
       do i = 1, nx

          !-----------------------------------------------------------------
          ! Compute total thickness and reciprocal thickness
          !-----------------------------------------------------------------

          thck = 0.d0
          do k = 1, nlyr
             thck = thck + hlyr(i,j,k)
          enddo

          if (thck > 0.d0) then
             rthck = 1.d0/thck
          else
             rthck = 0.d0
          endif

          !-----------------------------------------------------------------
          ! Determine vertical coordinate z1, given input layer thicknesses.
          ! These are the coordinates from which we start.
          !-----------------------------------------------------------------

          z1(1) = 0.d0
          do k = 2, nlyr
             z1(k) = z1(k-1) +  hlyr(i,j,k-1)*rthck 
          enddo
          z1(nlyr+1) = 1.d0
                       
          !-----------------------------------------------------------------
          ! Set vertical coordinate z2, given sigma.
          ! These are the coordinates to which we remap in the vertical.
          !-----------------------------------------------------------------

          z2(1) = 0.d0
          do k = 2, nlyr
             z2(k) = sigma(k)
          enddo
          z2(nlyr+1) = 1.d0

          !-----------------------------------------------------------------
          ! Compute new layer thicknesses (z2 coordinates)
          !-----------------------------------------------------------------

          do k = 1, nlyr
             hlyr(i,j,k) = (z2(k+1) - z2(k)) * thck
          enddo

          !-----------------------------------------------------------------
          ! Compute sum of h*T for each new layer (k2) by integrating
          ! over the regions of overlap with old layers (k1).
          !-----------------------------------------------------------------
             
          htsum(:,:) = 0.d0
          k1 = 1
          k2 = 1
          do while (k1 <= nlyr .and. k2 <= nlyr)
             zhi = min (z1(k1+1), z2(k2+1)) 
             zlo = max (z1(k1),   z2(k2))
             hovlp = max (zhi-zlo, 0.d0) * thck
             htsum(:,k2) = htsum(:,k2) +  trcr(i,j,:,k1) * hovlp
             if (z1(k1+1) > z2(k2+1)) then
                k2 = k2 + 1
             else
                k1 = k1 + 1
             endif
          enddo

          do k = 1, nlyr
             if (hlyr(i,j,k) > 0.d0) then
                trcr(i,j,:,k) = htsum(:,k) / hlyr(i,j,k)
             else
                trcr(i,j,:,k) = 0.d0
             endif
          enddo    ! k
          
       enddo       ! i
    enddo          ! j

    end subroutine glissade_vertical_remap

!=======================================================================

    subroutine upwind_field (nx,       ny,         &
                             ilo, ihi, jlo, jhi,   &
                             dx,       dy,         &
                             dt,       phi,        &
                             uee,      vnn)
      !
      ! first-order upwind transport algorithm
      !
      !
      ! Authors: Elizabeth Hunke and William Lipscomb, LANL
      !
      ! input/output arguments

      integer, intent (in) ::     &
         nx, ny             ,&! block dimensions
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real(dp), intent(in) ::         &
         dx, dy             ,&! x and y gridcell dimensions
         dt                   ! time step

      real(dp), dimension(nx,ny), &
         intent(inout) ::                       &
         phi                  ! scalar field

      real(dp), dimension(nx,ny),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities

      ! local variables

      integer ::     &
         i, j                   ! standard indices

      real(dp) ::        &
         upwind, y1, y2, a, h   ! function

      real(dp), dimension(nx,ny) ::  &
         worka, workb

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = 0.5d0*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do j = jlo-1, jhi
      do i = ilo-1, ihi
         worka(i,j)= upwind(phi(i,j),phi(i+1,j),uee(i,j),dy)
         workb(i,j)= upwind(phi(i,j),phi(i,j+1),vnn(i,j),dx)
      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi
         phi(i,j) = phi(i,j) - ( worka(i,j)-worka(i-1,j)      &
                               + workb(i,j)-workb(i,j-1) )    &
                               / (dx*dy)
      enddo
      enddo

    end subroutine upwind_field

!=======================================================================

  end module glissade_transport

!=======================================================================
