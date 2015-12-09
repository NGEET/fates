!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade.F90 - part of the Community Ice Sheet Model (CISM)  
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

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!CLEANUP - glissade.F90
!
! NOTE: MJH Lines that start with !### are ones I have identified to be deleted.
!
! This is a new module, originally copied from glide.F90 (William Lipscomb, June 2012)
! Removed SIA-specific code, leaving only the HO code with remapping transport
! Whenever possible, parallel_halo updates should go in this module rather
!  than at lower levels.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade.f90 - part of the CISM ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade

  ! Driver for Glissade (parallel, higher-order) dynamical core

  use glimmer_global, only: dp
  use glimmer_log
  use glide_types
  use glide_io
  use glide_lithot
  use glimmer_config
  use glissade_test, only: glissade_test_halo, glissade_test_transport

  implicit none

  integer, private, parameter :: dummyunit=99

  logical, parameter :: verbose_glissade = .false.

  ! Change either of the following logical parameters to true to carry out simple tests
  logical, parameter :: test_transport = .false.   ! if true, call test_transport subroutine
  real(dp), parameter :: thk_init = 500.d0         ! initial thickness (m) for test_transport
  logical, parameter :: test_halo = .false.        ! if true, call test_halo subroutine

  !WHL - for trying glissade_therm in place of glissade_temp
!!  logical, parameter :: call_glissade_therm = .false.
  logical, parameter :: call_glissade_therm = .true.

contains

!=======================================================================

! Note: There is no glissade_config subroutine; glide_config works for all dycores.

!=======================================================================

  subroutine glissade_initialise(model)

    ! initialise Glissade model instance

    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo  !TODO - Remove call to init_velo?
    !TODO - Remove glissade_temp
    use glissade_temp, only: glissade_init_temp
    use glissade_therm, only: glissade_init_therm
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glide_thck, only : glide_calclsrf
    use glam_strs2, only : glam_velo_init
    use glimmer_coordinates, only: coordsystem_new
    use glide_grid_operators, only: stagvarb
    use glissade_velo_higher, only: glissade_velo_higher_init
    use glide_diagnostics, only: glide_init_diag
    use felix_dycore_interface, only: felix_velo_init
    use glide_bwater
    use glimmer_paramets, only: thk0

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    !TODO - Is glimmer_version_char sitll needed?
    character(len=100), external :: glimmer_version_char

    integer :: i, j

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)

    ! set up coordinate systems, and change to the parallel values of ewn and nsn

    ! With outflow BCs, scalars in the halos are set to zero.
    if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
       call distributed_grid(model%general%ewn,model%general%nsn,outflow_bc_in=.true.)
    else
       call distributed_grid(model%general%ewn,model%general%nsn)
    endif

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0,               &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.d0, model%numerics%dns/2.d0, &
                                              model%numerics%dew,      model%numerics%dns,      &
                                              model%general%ewn-1,     model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

    ! set uniform basal heat flux (positive down)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)

    ! read first time slice
    call glide_io_readall(model,model)

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)

    select case(model%options%whichrelaxed)
    case(RELAXED_TOPO_INPUT)   ! supplied topography is relaxed
       model%isostasy%relx = model%geometry%topg
    case(RELAXED_TOPO_COMPUTE) ! supplied topography is in equilibrium
                               !TODO - Test the case RELAXED_TOPO_COMPUTE
       call not_parallel(__FILE__,__LINE__)
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model, model)

    ! If a 2D bheatflx field is present in the input file, it will have been written 
    !  to model%temper%bheatflx.  For the case model%options%gthf = 0, we want to use
    !  a uniform heat flux instead.
    ! If no bheatflx field is present in the input file, then we default to the 
    !  prescribed uniform value, model%paramets%geot.

    if (model%options%gthf == GTHF_UNIFORM) then

       ! Check to see if this flux was present in the input file
       ! (by checking whether the flux is nonuniform over the domain)
       if (abs(maxval(model%temper%bheatflx) - minval(model%temper%bheatflx)) > 1.d-6) then  
          call write_log('Setting uniform prescribed geothermal flux')
          call write_log('(Set gthf = 1 to read geothermal flux field from input file)')
       endif

       ! set uniform basal heat flux (positive down)
       model%temper%bheatflx = model%paramets%geot

    endif

    ! initialise glissade components

    ! Update some variables in halo cells
    ! Note: We need thck and artm in halo cells so that temperature will be initialized correctly (if not read from input file).
    !       We do an update here for temp in case temp is read from an input file.
    !       If temp is computed in glissade_init_therm (based on the value of options%temp_init),
    !        then the halos will receive the correct values.
    !TODO - Does anything else need an initial halo update?
    call parallel_halo(model%geometry%thck)
    call parallel_halo(model%climate%artm)
    call parallel_halo(model%temper%temp)

    !TODO - Remove call to init_velo in glissade_initialise?
    !       Most of what's done in init_velo is needed for SIA only, but still need velowk for call to wvelintg
    call init_velo(model)

    !TODO - Remove glissade_init_temp option
    if (call_glissade_therm) then
       call glissade_init_therm(model%options%temp_init,    model%options%is_restart,  &
                                model%general%ewn,          model%general%nsn,         &
                                model%general%upn,                                     &
                                model%numerics%sigma,       model%numerics%stagsigma,  &
                                model%geometry%thck*thk0,                              & ! m
                                model%climate%artm,                                    & ! deg C
                                model%temper%temp)                                       ! deg C
    else
       call glissade_init_temp(model)
    endif

    ! Initialize basal hydrology model, if enabled
    call bwater_init(model)

    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call init_lithot(model)
    end if

    ! Dycore-specific velocity solver initialization
    select case (model%options%whichdycore)
    case ( DYCORE_GLAM )   ! glam finite-difference

       call glam_velo_init(model%general%ewn,    model%general%nsn,  &
                           model%general%upn,                        &
                           model%numerics%dew,   model%numerics%dns, &
                           model%numerics%sigma)

    case ( DYCORE_GLISSADE )   ! glissade finite-element

       call glissade_velo_higher_init

    case ( DYCORE_ALBANYFELIX)

       call felix_velo_init(model)

    end select

    !TODO - Add halo updates of state variables here?

    ! If unstagbeta (i.e., beta on the scalar ice grid) was read from an input file,
    !  then interpolate it to beta on the staggered grid.
    ! NOTE: unstagbeta is initialized to unphys_val = -999.d0, so its maxval will be > 0 only if
    !       the field is read in.
    ! We can make an exception for ISHOM case C; for greater accuracy we set beta in 
    !  subroutine calcbeta instead of interpolating from unstagbeta (one processor only).

    if (maxval(model%velocity%unstagbeta) > 0.d0 .and.   &
               model%options%which_ho_babc /= HO_BABC_ISHOMC) then  ! interpolate to staggered grid
       call write_log('Interpolating beta from unstaggered (unstagbeta) to staggered grid (beta)')
       if (maxval(model%velocity%beta) > 0.0d0 ) then
          call write_log('Warning: the input "beta" field will be overwritten with values interpolated from the input "unstagbeta" field!')
       endif

       call parallel_halo(model%velocity%unstagbeta)    ! fill in halo values
       call stagvarb(model%velocity%unstagbeta,  &      ! interpolate
                     model%velocity%beta,        &
                     model%general%ewn,          &
                     model%general%nsn)

    endif  ! unstagbeta > 0

    ! Note: The basal process option is currently disabled.
    ! initialize basal process module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%dttem)
!!    end if      

    ! calculate mask
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    ! and calculate lower and upper ice surface

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

    ! initialise model diagnostics                                                                                                                 \
                                                                                                                                                  
    call glide_init_diag(model)

    ! optional unit tests

    if (test_halo) then
       call glissade_test_halo (model)
       call parallel_finalise
    endif
     
    if (test_transport) then
       where (model%geometry%thck > model%numerics%thklim)
          model%geometry%thck = thk_init/thk0
       elsewhere
          model%geometry%thck = 0.d0
       endwhere
    endif

    ! Initial solve of calcbwat
    ! TODO: Should call to calcbwat go here or in diagnostic solve routine? Make sure consistent with Glide.
    call calcbwat(model, &
                  model%options%whichbwat, &
                  model%temper%bmlt, &
                  model%temper%bwat, &
                  model%temper%bwatflx, &
                  model%geometry%thck, &
                  model%geometry%topg, &
                  model%temper%temp(model%general%upn,:,:), &
                  GLIDE_IS_FLOAT(model%geometry%thkmask), &
                  model%tempwk%wphi)

  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time)

    ! Perform time-step of an ice model instance with the Glissade dycore

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0, thk0
    use glimmer_scales, only: scale_acab
    use glimmer_physcon, only: scyr
    !TODO - Remove glissade_temp option
    use glissade_temp, only: glissade_temp_driver
    use glissade_therm, only: glissade_therm_driver, glissade_temp2enth, glissade_enth2temp
    use glide_mask, only: glide_set_mask, calc_iareaf_iareag
    use glide_ground, only: glide_marinlim
    use glide_grid_operators
    use isostasy
    use glissade_enthalpy
    use glissade_transport, only: glissade_transport_driver, glissade_check_cfl, ntracer
    use glissade_grid_operators
    use glide_thck, only: glide_calclsrf
    use glide_bwater

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp), intent(in) :: time         ! current time in years

    ! --- Local variables ---

    integer :: sc  ! subcycling index

    ! temporary thck and acab arrays in SI units
    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       thck_unscaled,      &! ice thickness (m)
       acab_unscaled        ! surface mass balance (m/s)

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! temporary bmlt array
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
       bmlt_continuity  ! = bmlt if basal mass balance is included in continuity equation
                        ! else = 0

    logical :: do_upwind_transport  ! logical for whether transport code should do upwind transport or incremental remapping
                                    ! set to true for EVOL_UPWIND, else = false

    integer :: i, j, k
    integer :: nx, ny
    integer :: ewn, nsn, upn
    
    !WHL - debug
    real(dp) :: thck_west, thck_east, dthck, u_west, u_east

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    ! ========================

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%timecounter = model%numerics%timecounter + 1
    model%temper%newtemps = .false.

    ! optional transport test
    ! code execution will end when this is done
    if (test_transport) then
       call glissade_test_transport (model)
       return
    endif

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    !TODO Not sure if this is in the right place.  G1=f(G0,T0) and T1=g(G0,T0)  
    !     If we update G1 now, then we will be doing T1=g(G1,T0).
    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! Vertical diffusion and strain heating only; no advection
    ! ------------------------------------------------------------------------ 

    ! Note: These times have units of years
    !       dttem has scaled units, so multiply by tim0/scyr to convert to years

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%dttem*tim0/scyr)) then

      call t_startf('glissade_therm_driver')

      !TODO - Remove glissade_temp option
      if (call_glissade_therm) then

         if (main_task .and. verbose_glissade) print*, 'Call glissade_therm_driver'

         ! Note: glissade_therm_driver uses SI units
         !       Output arguments are temp, waterfrac and bmlt
         call glissade_therm_driver (model%options%whichtemp,                                      &
                                     model%numerics%dttem*tim0,                                    & ! s
                                     model%general%ewn,          model%general%nsn,                &
                                     model%general%upn,                                            &
                                     model%numerics%idiag_local, model%numerics%jdiag_local,       &
                                     model%numerics%rdiag_local,                                   &
                                     model%numerics%sigma,       model%numerics%stagsigma,         &
                                     model%numerics%thklim*thk0, model%numerics%thklim_temp*thk0,  & ! m
                                     model%geometry%thck*thk0,                                     & ! m
                                     model%geometry%topg*thk0,   model%climate%eus*thk0,           & ! m
                                     model%climate%artm,                                           & ! deg C    
                                     model%temper%bheatflx,      model%temper%bfricflx,            & ! W/m2
                                     model%temper%dissip,                                          & ! deg/s
                                     model%temper%bwat*thk0,                                       & ! m
                                     model%temper%temp,                                            & ! deg C
                                     model%temper%waterfrac,                                       & ! unitless
                                     model%temper%bmlt)                                              ! m/s on output
                                     
         ! convert bmlt from m/s to scaled model units
         model%temper%bmlt = model%temper%bmlt * tim0/thk0
                                     
      else
         if (main_task .and. verbose_glissade) print*, 'Call glissade_temp_driver'
         call glissade_temp_driver(model, model%options%whichtemp)
      endif
      call t_stopf('glissade_therm_driver')

      model%temper%newtemps = .true.

      ! Update basal hydrology, if needed
      call calcbwat( model,                                    &
                     model%options%whichbwat,                  &
                     model%temper%bmlt,                        &
                     model%temper%bwat,                        &
                     model%temper%bwatflx,                     &
                     model%geometry%thck,                      &
                     model%geometry%topg,                      &
                     model%temper%temp(model%general%upn,:,:), &
                     GLIDE_IS_FLOAT(model%geometry%thkmask),   &
                     model%tempwk%wphi)

    end if  ! take a temperature time step

    !------------------------------------------------------------------------ 
    ! Halo updates
    !------------------------------------------------------------------------ 

    call parallel_halo(model%temper%bwat)    !TODO: not sure halo update is needed for bwat

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    ! MJH: This now uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme
    ! WHL: We used to have EVOL_NO_THICKNESS = -1 as a Glide option, used to hold the ice surface elevation fixed during CESM runs.  
    !      This option has been replaced by a Glint option, evolve_ice.
    !      We now have EVOL_NO_THICKESS = 5 as a glam/glissade option.  It is used to hold the ice surface elevation fixed
    !       while allowing temperature to evolve, which can be useful for model spinup.  This option might need more testing.

    select case(model%options%whichevol)

       case(EVOL_INC_REMAP, EVOL_UPWIND, EVOL_NO_THICKNESS) 

       if (model%options%whichevol == EVOL_UPWIND) then
          do_upwind_transport = .true.
       else
          do_upwind_transport = .false.
       endif

       ! Use incremental remapping scheme to transport ice thickness (and temperature too, if whichtemp = TEMP_PROGNOSTIC).
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !       at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.  Do we want
       !       to support temperature evolution without thickness evolution?  If so, then the current implementation may be preferred approach.

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! store old thickness
          thck_old = model%geometry%thck
          stagthck_old = model%geomderv%stagthck
       endif

      call t_startf('inc_remap_driver')

       if (main_task) then
          print *, 'Compute dH/dt'
       endif

       ! Halo updates for velocities, thickness and tracers
       ! Velocity update might be needed if velo was not updated in halo at the end of the previous diagnostic solve
       !  (just to be on the safe side).

      call t_startf('new_remap_halo_upds')

      call staggered_parallel_halo(model%velocity%uvel)
      call staggered_parallel_halo(model%velocity%vvel)

      call parallel_halo(model%geometry%thck)
      call parallel_halo(model%temper%temp)

      if (model%options%whichtemp == TEMP_ENTHALPY) then
         call parallel_halo(model%temper%waterfrac)
      endif

      call t_stopf('new_remap_halo_upds')

      call t_startf('glissade_transport_driver')

      if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
         bmlt_continuity(:,:) = model%temper%bmlt(:,:) * thk0/tim0   ! convert to m/s
      else                                                           ! do not include bmlt in continuity equation
         bmlt_continuity(:,:) = 0.d0
      endif

      ! --- First determine CFL limits ---
      ! Note we are using the subcycled dt here (if subcycling is on).
      ! (see note above about the EVOL_NO_THICKNESS option and how it is affected by a CFL violation)
      ! stagthck, dusrfdew/ns and u/vvel need to be from the previous time step (and are at this point)
      call glissade_check_cfl(model%general%ewn,         model%general%nsn,         model%general%upn-1,      &
                              model%numerics%dew * len0, model%numerics%dns * len0, model%numerics%sigma,     &
                              model%geomderv%stagthck * thk0,                                                 &
                              model%geomderv%dusrfdew*thk0/len0, model%geomderv%dusrfdns*thk0/len0,           &
                              model%velocity%uvel * scyr * vel0, model%velocity%vvel * scyr * vel0,           &
                              model%numerics%dt_transport * tim0 / scyr,                                      &
                              model%numerics%adv_cfl_dt,         model%numerics%diff_cfl_dt )

       ! Call the transport driver.
       ! Note: This subroutine assumes SI units:
       !       * dt (s)
       !       * dew, dns, thck (m)
       !       * uvel, vvel, acab, blmt (m/s)
       !       Since thck has intent(inout), we create and pass a temporary array with units of m.
       ! TODO - Pass ice age as tracer 2

       do sc = 1, model%numerics%subcyc
          if (model%numerics%subcyc > 1 .and. main_task) write(*,*) 'Subcycling transport: Cycle ',sc

          ! temporary in/out arrays in SI units (m)                               
          thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0
          acab_unscaled(:,:) = model%climate%acab(:,:) * thk0/tim0

          if (model%options%whichtemp == TEMP_PROGNOSTIC) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)
                                                                ! Note: We are passing arrays in SI units.

             call glissade_transport_driver(model%numerics%dt_transport * tim0,                   &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%temp(:,:,:),                             &
                                            upwind_transport_in = do_upwind_transport )

             ! convert thck and acab back to scaled units
             model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0
             model%climate%acab(:,:) = acab_unscaled(:,:) / (thk0/tim0)

          elseif (model%options%whichtemp == TEMP_ENTHALPY) then  ! Use IR to transport thickness and enthalpy

             ! Derive enthalpy from temperature and waterfrac
             ! Note: glissade_temp2enth expects SI units
             do j = 1, model%general%nsn 
                do i = 1, model%general%ewn
                   call glissade_temp2enth (model%numerics%stagsigma(1:upn-1),        &
                                            model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j),   &
                                            model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j))
                enddo
             enddo

             ! Transport fields, with enthalpy as a tracer instead of temperature
             call glissade_transport_driver(model%numerics%dt_transport * tim0,                   &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      & 
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%enthalpy(:,:,:),                         & 
                                            upwind_transport_in = do_upwind_transport )

          else  ! Use IR to transport thickness only
                ! Note: In glissade_transport_driver, the ice thickness is transported layer by layer,
                !       which is inefficient if no tracers are being transported.  (It would be more
                !       efficient to transport thickness in one layer only, using a vertically
                !       averaged velocity.)  
                ! Not sure if this option will be used in practice.

             call glissade_transport_driver(model%numerics%dt_transport * tim0,                   &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    & 
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:) ,                                &
                                            upwind_transport_in = do_upwind_transport )

          endif  ! whichtemp

          ! convert thck and acab back to scaled units
          model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0
          model%climate%acab(:,:) = acab_unscaled(:,:) / (thk0/tim0)
         
          if (verbose_glissade) then
             print*, ' '
             print*, 'After glissade_transport_driver:'
             print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
             print*, 'max, min acab (m/yr) =', maxval(model%climate%acab)*scale_acab, minval(model%climate%acab)*scale_acab
             print*, 'max, min artm =', maxval(model%climate%artm), minval(model%climate%artm)
             print*, 'thklim =', model%numerics%thklim * thk0
             print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
             print*, ' '
             print*, 'thck:'
             do j = model%general%nsn, 1, -1
                do i = 1, model%general%ewn
                   write(6,'(f6.0)',advance='no') model%geometry%thck(i,j) * thk0
                enddo
                write(6,*) ' '
             enddo
          endif

          ! Update halos of modified fields

         call t_startf('after_remap_haloupds')

         call parallel_halo(model%geometry%thck)

         if (model%options%whichtemp == TEMP_ENTHALPY) then

             ! Update enthalpy in halo cells
             call parallel_halo(model%temper%enthalpy)

             ! Derive new temperature and waterfrac from enthalpy (will be correct in halo cells)
             ! Note: glissade_enth2temp expects SI units
             do j = 1, model%general%nsn 
                do i = 1, model%general%ewn 
                   call glissade_enth2temp(model%numerics%stagsigma(1:upn-1),                                    &
                                           model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j),   &
                                           model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j))
                enddo
             enddo

         else   ! update temperature in halo cells

            call parallel_halo(model%temper%temp)

         endif    ! TEMP_ENTHALPY

         ! NOTE: Halo updates of other tracers, if present, should go here

         call t_stopf('after_remap_haloupds')

       enddo     ! subcycling

       call t_stopf('glissade_transport_driver')

       call t_stopf('inc_remap_driver')

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

    end select

    ! TODO: Not sure topg should be updated here; should be updated after isostasy
    !       if the isostasy is turned on.
    call parallel_halo(model%geometry%topg)

    !------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    !------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    ! --- Calculate updated mask because marinlim calculation needs a mask.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask)

    !TODO - Look at marinlim more carefully and see which halo updates are necessary.
    !       It appears that marinlim only needs the halo of thkmask for case 5 (which was removed).  
    !
    !       glide_set_mask includes a halo update of model%geometry%thkmask; remove this one?
    !       Do we need a halo update for relx? If so, should this be done at initialization?

    call parallel_halo(model%geometry%thkmask)
    call parallel_halo(model%isostasy%relx)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed depth, 
    ! depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

    call glide_marinlim(model%options%whichmarn, &
                        model%geometry%thck,  &
                        model%isostasy%relx,      &
                        model%geometry%topg,  &
                        model%geometry%thkmask,    &
                        model%numerics%mlimit,     &
                        model%numerics%calving_fraction, &
                        model%climate%eus,         &
                        model%climate%calving,  &
                        model%ground, &
                        model%numerics%dew,    &
                        model%numerics%dns, &
                        model%general%nsn, &
                        model%general%ewn)

    !TODO: Think about what halo updates are needed after glide_marinlim. Just thck and thkmask?

    ! halo updates
    call parallel_halo(model%geometry%thck)    ! Updated halo values of thck are needed below in calc_lsrf

    !TODO - Remove this call to glide_set_mask?
    !       This subroutine is called at the beginning of glissade_velo_driver,
    !        so a call here is not needed for the velo diagnostic solve.
    !       The question is whether it is needed for the isostasy.

    ! --- marinlim adjusts thickness for calved ice.  Therefore the mask needs to be recalculated.
    ! --- This time we want to calculate the optional arguments iarea and ivol because thickness 
    ! --- will not change further during this time step.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    !Note: glide_set_mask includes a halo update of model%geometry%thkmask at end of call
    ! call parallel_halo(model%geometry%thkmask)

    ! --- Calculate global area of ice that is floating and grounded.
    !TODO  May want to calculate iareaf and iareag in glide_write_diag and remove those calculations here.  

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%thkmask,                        &
                            model%geometry%iareaf, model%geometry%iareag)

    !TODO - Are these isostasy calls in the right place?
    ! Consider for a forward Euler time step:
    ! With a relaxing mantle model, topg is a prognostic (time-evolving) variable (I think):
    !      topg1 = f(topg0, thk0, ...) 
    ! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable 
    !(comparable to calculating floatation height of ice in the ocean):
    !      topg1 = f(thk1)
    ! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
    ! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.  
    ! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.  
    ! Also, we should think about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.
    
    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       if (model%numerics%time >= model%isostasy%next_calc) then
          model%isostasy%next_calc = model%isostasy%next_calc + model%isostasy%period
          call isos_icewaterload(model)
          model%isostasy%new_load = .true.
       end if
    end if
   
      ! calculate isostatic adjustment and upper and lower ice surface

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    !TODO - Test the local isostasy schemes in the parallel model.
    !       The elastic lithosphere scheme is not expected to work in parallel.

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call isos_compute(model)
    end if

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !TODO - Any halo updates needed at the end of glissade_tstep?

  end subroutine glissade_tstep

!=======================================================================

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables, in particular the ice velocity.
     ! This is needed at the end of each time step once the prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0, thk0, vis0, tau0, evs0
    use glimmer_physcon, only: scyr
    use glide_thck, only: glide_calclsrf
    use glissade_temp, only: glissade_calcflwa
    use glam_velo, only: glam_velo_driver, glam_basal_friction
    use glissade_velo, only: glissade_velo_driver
    use glide_velo, only: wvelintg
    use glissade_masks, only: glissade_get_masks
    use glissade_therm, only: glissade_interior_dissipation_sia,  &
                              glissade_interior_dissipation_first_order, &
                              glissade_flow_factor
    use glam_grid_operators, only: glam_geometry_derivs
    use felix_dycore_interface, only: felix_velo_driver

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer :: i, j, k
    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,     &! = 1 where thck > thklim, else = 0
         floating_mask  ! = 1 where ice is floating, else = 0

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 1. First part of diagnostic solve: 
    !    Now that advection is done, update geometry- and temperature-related 
    !    diagnostic fields that are needed for the velocity solve.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! ------------------------------------------------------------------------
    ! Halo updates for ice topography and thickness
    !
    ! NOTE: There is an optional argument periodic_offset_ew for topg.
    !       This is for ismip-hom experiments. A positive EW offset means that 
    !        the topography in west halo cells will be raised, and the topography 
    !        in east halo cells will be lowered.  This ensures that the topography
    !        and upper surface elevation are continuous between halo cells
    !        and locally owned cells at the edge of the global domain.
    !       In other cases (anything but ismip-hom), periodic_offset_ew = periodic_offset_ns = 0, 
    !        and this argument will have no effect.
    ! ------------------------------------------------------------------------

    call parallel_halo(model%geometry%thck)
    call parallel_halo(model%geometry%topg, periodic_offset_ew = model%numerics%periodic_offset_ew)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    ! ------------------------------------------------------------------------
    ! Update some geometry derivatives
    ! ------------------------------------------------------------------------
    !Note - The fields computed by glam_geometry_derivs are not required by
    !       the Glissade velocity solver (which computes them internally).  
    !       However, some of the fields (stagthck, dusrfdew and dusrfdns) 
    !       are needed during the next timestep by glissade_temp
    !       if we're doing shallow-ice dissipation.
    !TODO - Replace this glam_geometry_derivs call with calls to Glissade subroutines?
    !       (The glam_velo driver includes its own call to glam_geometry_derivs.) 

    call glam_geometry_derivs(model)

    ! ------------------------------------------------------------------------
    ! Update some masks that are used by Glissade subroutines
    ! ------------------------------------------------------------------------

    call glissade_get_masks(model%general%ewn,   model%general%nsn,     &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,            floating_mask)

    ! ------------------------------------------------------------------------ 
    ! Calculate Glen's A
    !
    ! Notes:
    ! (1) Because flwa is not a restart variable in Glissade, no check is included 
    !      here for whether to calculate it on initial time (as is done in Glide).
    ! (2) We are passing in only vertical elements (1:upn-1) of the temp array,
    !       so that it has the same vertical dimensions as flwa.
    ! (3) The flow enhancement factor is 1 by default.
    ! (4) The waterfrac field is ignored unless whichtemp = TEMP_ENTHALPY.
    ! (5) Inputs and outputs of glissade_flow_factor should have SI units.
    ! ------------------------------------------------------------------------

    call glissade_flow_factor(model%options%whichflwa,            &
                              model%options%whichtemp,            &
                              model%numerics%stagsigma,           &
                              model%geometry%thck * thk0,         &  ! scale to m
                              ice_mask,                           &
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%temper%flwa,                  &  ! Pa^{-n} s^{-1}
                              model%paramets%default_flwa / scyr, &  ! scale to Pa^{-n} s^{-1}
                              model%paramets%flow_enhancement_factor,   &
                              model%temper%waterfrac(:,:,:))

    ! Change flwa to model units (glissade_flow_factor assumes SI units of Pa{-n} s^{-1})
    model%temper%flwa(:,:,:) = model%temper%flwa(:,:,:) / vis0

    !TODO - flwa halo update not needed?
    ! Halo update for flwa
    call parallel_halo(model%temper%flwa)

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then
  
       call write_log('Using uvel, vvel from restart file at initial time')

    else

       ! If this is not a restart or we are not at the initial time, then proceed normally.

       if ( (model%numerics%time == model%numerics%tstart) .and. &
         ( (maxval(abs(model%velocity%uvel)) /= 0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel)) /= 0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value will change the answer only within the tolerance of the nonlinear solve.  
          ! If a user already has a good guess from a previous run, they may wish to start things off with it to speed the initial solution.
          call write_log('Using uvel, vvel from input file as initial guess at initial time.')
          call write_log('If this is not desired, please remove those fields from the input file.')
       endif

       if (main_task) then
          print *, ' '
          print *, 'Compute ice velocities, time =', model%numerics%time
       endif

       !! extrapolate value of mintauf into halos to enforce periodic lateral bcs (only if field covers entire domain)
       if (model%options%which_ho_babc == HO_BABC_YIELD_PICARD) then
          call staggered_parallel_halo_extrapolate(model%basalproc%mintauf)
       endif

       ! Call the appropriate velocity solver

       select case (model%options%whichdycore)
       case ( DYCORE_GLAM )   ! glam finite-difference

          call t_startf('glam_velo_driver')
          call glam_velo_driver(model)
          call t_stopf('glam_velo_driver')

       case ( DYCORE_GLISSADE )   ! glissade finite-element

          call t_startf('glissade_velo_driver')
          call glissade_velo_driver(model)
          call t_stopf('glissade_velo_driver')

       case ( DYCORE_ALBANYFELIX)

          call t_startf('felix_velo_driver')
          call felix_velo_driver(model)
          call t_stopf('felix_velo_driver')

       end select
 
       ! Compute internal heat dissipation
       ! This is used in the prognostic temperature calculation during the next time step.
       ! Note: These glissade subroutines assume SI units on input and output

       model%temper%dissip(:,:,:) = 0.d0

       if (model%options%which_ho_disp == HO_DISP_SIA) then

          call glissade_interior_dissipation_sia(model%general%ewn,              &
                                                 model%general%nsn,              &
                                                 model%general%upn,              &
                                                 model%numerics%stagsigma(:),    &
                                                 ice_mask,                       &
!                                                 model%geomderv%stagthck,     &
!                                                 model%temper%flwa,           &
!                                                 model%geomderv%dusrfdew,     &
!                                                 model%geomderv%dusrfdns,     &
                                                 model%geomderv%stagthck * thk0, & ! scale to m
                                                 model%temper%flwa * vis0,       & ! scale to Pa^{-n} s^{-1}
                                                 model%geomderv%dusrfdew * thk0/len0, & ! scale to m/m
                                                 model%geomderv%dusrfdns * thk0/len0, & ! scale to m/m
                                                 model%temper%dissip)
          
       else    ! first-order dissipation                                                                                                                                                               
          call glissade_interior_dissipation_first_order(model%general%ewn,          &
                                                         model%general%nsn,          &
                                                         model%general%upn,          &
                                                         ice_mask,                   &
                                                         model%stress%tau%scalar * tau0,  &  ! scale to Pa
                                                         model%stress%efvs * evs0,   &  ! scale to Pa s
                                                         model%temper%dissip)
          
       endif    ! which_ho_disp
          
       ! If running Glam, compute the basal friction heat flux
       ! (Glissade computes this flux as part of the velocity solution.)
       
       if (model%options%whichdycore == DYCORE_GLAM) then
          call glam_basal_friction(model%general%ewn,                             &
                                   model%general%nsn,                             &
                                   ice_mask,                                      &
                                   floating_mask,                                 &
                                   model%velocity%uvel(model%general%upn,:,:),    &
                                   model%velocity%vvel(model%general%upn,:,:),    &
                                   model%velocity%btraction(:,:,:),               &
                                   model%temper%bfricflx(:,:) )
       endif
       
       if (this_rank==model%numerics%rdiag_local .and. verbose_glissade) then
          i = model%numerics%idiag_local
          j = model%numerics%jdiag_local
          print*, 'k, dissip (deg/yr):'
          do k = 1, model%general%upn-1
             print*, k, model%temper%dissip(k,i,j)*scyr
          enddo
          print*, 'ubas, vbas =', model%velocity%uvel(model%general%upn,i,j),  &
                                  model%velocity%vvel(model%general%upn,i,j)
          print*, 'btraction =',  model%velocity%btraction(:,i,j)
          print*, 'bfricflx =', model%temper%bfricflx(i,j)
       endif

       if (main_task .and. verbose_glissade) then
          print*, ' '
          print*, 'After glissade velocity solve: uvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%uvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo

          print*, ' '
          print*, 'After glissade velocity solve: vvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%vvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo
       endif  ! main_task & verbose_glissade

    endif     ! is_restart

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! compute the velocity norm and basal velocity (for diagnostic output)

    k = model%general%upn
    model%velocity%ubas(:,:) = model%velocity%uvel(k,:,:)
    model%velocity%vbas(:,:) = model%velocity%vvel(k,:,:)
    model%velocity%velnorm(:,:,:) = sqrt(model%velocity%uvel(:,:,:)**2 + model%velocity%vvel(:,:,:)**2)

    ! Copy uvel and vvel to arrays uvel_extend and vvel_extend.
    ! These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
    ! Thus they are better suited for I/O if we have periodic BC,
    !  where the velocity field we are solving for has global dimensions (nx,ny).
    ! Since uvel and vvel are not defined for i = nx or j = ny, the
    !  uvel_extend and vvel_extend arrays will have values of zero at these points.
    ! But these are halo points, so when we write netCDF I/O it shouldn't matter;
    !  we should have the correct values at physical points.
    
    model%velocity%uvel_extend(:,:,:) = 0.d0
    model%velocity%vvel_extend(:,:,:) = 0.d0

    do j = 1, model%general%nsn-1
       do i = 1, model%general%ewn-1
          model%velocity%uvel_extend(:,i,j) = model%velocity%uvel(:,i,j)
          model%velocity%vvel_extend(:,i,j) = model%velocity%vvel(:,i,j)             
       enddo
    enddo
 
    ! Calculate wvel, assuming grid velocity is 0.
    ! This is calculated relative to ice sheet base, rather than a fixed reference location
    ! Note: This current implementation for wvel only supports whichwvel=VERTINT_STANDARD
    call wvelintg(model%velocity%uvel,                        &
                  model%velocity%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%geometry%thck * 0.0d0,                &  ! Just need a 2d array of all 0's for wgrd
                  model%geometry%thck,                        &
                  model%temper%bmlt,                          &
                  model%velocity%wvel_ho)
    ! Note: halos may be wrong for wvel_ho, but since it is currently only used as an output diagnostic variable, that is OK.

    !TODO - I don't think we need to update ubas, vbas, or velnorm, since these are diagnostic only
    !       Also, I don't think efvs is needed in the halo.
    call staggered_parallel_halo(model%velocity%velnorm)
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)
    call parallel_halo(model%stress%efvs)

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

end module glissade

!=======================================================================
