!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_traction.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"
#include "config.inc"

  module glissade_basal_traction

  !-----------------------------------------------------------------------------
  ! Compute or prescribe the basal traction coefficient 'beta' as required by
  ! the higher-order velocity solver.
  ! 
  ! Note that beta is assumed to be a positive constant.  In earlier versions of
  ! the code it was called 'betasquared'.
  !
  ! The units are Pa/(m/yr) if we assume a linear sliding law of the form
  !    taub_x = -beta * u, taub_y = -beta * v
  !
  ! However, the units are Pa if beta is treated as a till yield stress.
  !
  ! Current options are as follows:
  ! 
  ! [0] constant value of 10 Pa/(m/yr) (useful for debugging)
  ! [1] simple hard-coded pattern (useful for debugging)
  ! [2] treat beta value as a till yield stress (in Pa) using Picard iteration 
  ! [3] linear (inverse) function of basal water depth (bwat) 
  ! [4] very large value for beta to enforce no slip everywhere 
  ! [5] beta field passed in from .nc input file as part of standard i/o
  ! [6] no slip everywhere (using Dirichlet BC rather than large beta)
  ! [7] treat beta value as till yield stress (in Pa) using Newton-type iteration (in devel.)
  ! [8] set beta as prescribed for ISMIP-HOM test C (serial only)
  ! [9] power law that used effective pressure
  ! [10] Coulomb friction law

  ! TODO - Renumber HO_BABC options so that, for example, the no-slip options have small numbers?
  !-----------------------------------------------------------------------------

  use glimmer_paramets, only : dp
  use glimmer_physcon,  only : scyr
  use glimmer_paramets, only : vel0, tau0
  use glimmer_log
  use glide_types
  use parallel,         only : staggered_parallel_halo  
  use glissade_grid_operators

  implicit none

!***********************************************************************

contains

!***********************************************************************

  subroutine calcbeta (whichbabc,                    &
                       dew,           dns,           &
                       ewn,           nsn,           &
                       thisvel,       othervel,      &
                       bwat,          beta_const,    &
                       mintauf,       basal_physics, &
                       flwa_basal,    thck,          &
                       mask,          beta)

  ! subroutine to calculate map of beta sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
   
  ! NOTE: Previously, the input arguments were assumed to be dimensionless
  ! and were rescaled in this routine.  Now the input arguments are
  ! assumed to have the units given below.
     
  use glimmer_paramets, only: len0
  use glimmer_physcon, only: gn
  use parallel, only: nhalo

  implicit none

  ! Input/output arguments

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in)                    :: dew, dns           ! m
  real(dp), intent(in), dimension(:,:)    :: thisvel, othervel  ! basal velocity components (m/yr)
  real(dp), intent(in), dimension(:,:)    :: bwat     ! basal water depth (m)
  real(dp), intent(in), dimension(:,:)    :: mintauf  ! till yield stress (Pa)
  real(dp), intent(in)                    :: beta_const  ! spatially uniform beta (Pa yr/m)
  type(glide_basal_physics), intent(in) :: basal_physics  ! basal physics object
  real(dp), intent(in), dimension(:,:) :: flwa_basal  ! flwa for the basal ice layer
  real(dp), intent(in), dimension(:,:) :: thck  ! ice thickness
  integer, intent(in), dimension(:,:)     :: mask ! staggered grid mask
  real(dp), intent(inout), dimension(:,:) :: beta  ! (Pa yr/m)

!WHL - These masks are no longer used
!  logical, intent(in), dimension(:,:), optional ::  &
!     floating_mask,   &! = 1 for cells where ice is present and is floating
!     ocean_mask        ! = 1 for cells where topography is below sea level and ice is absent
!                       ! Note: These masks live on the scalar grid; beta lives on the staggered grid
   
  ! Local variables

  real(dp) :: smallnum = 1.0d-2  ! m/yr

  integer :: ew, ns

  ! SFP added for making beta a function of basal water flux 
  real(dp), dimension(:,:), allocatable :: unstagbeta
  real(dp) :: C, m

  real(dp) :: Ldomain   ! size of full domain
  real(dp) :: omega     ! frequency of beta field
  real(dp) :: dx, dy
  integer :: ilo, ihi, jlo, jhi  ! limits of beta field for ISHOM C case
  integer :: i, j

  ! variables for power law
  real(dp) :: p, q

  ! variables for Coulomb friction law
  real(dp) :: Coulomb_C   ! friction coefficient
  real(dp) :: lambda_max  ! wavelength of bedrock bumps at subgrid scale
  real(dp) :: m_max       ! maximum bed obstacle slope
  real(dp), dimension(size(beta,1), size(beta,2)) :: big_lambda       ! bed rock characteristics
  integer, dimension(size(thck,1), size(thck,2))  :: imask            ! ice grid mask  1=ice, 0=no ice
  real(dp), dimension(size(beta,1), size(beta,2)) :: flwa_basal_stag  ! flwa for the basal ice layer on the staggered grid


  select case(whichbabc)

    case(HO_BABC_CONSTANT)  ! spatially uniform value; useful for debugging and test cases

!      beta(:,:) = 10.d0       ! This is the default value (Pa yr/m)
      beta(:,:) = beta_const   ! Pa yr/m

      ! If floating and ocean masks are passed in, then set beta to zero for shelf/ocean nodes,
      ! overriding the constant value set above. This allows us to model large regions
      ! (e.g., a whole ice sheet) in a simple but physically sensible way without specifying 
      ! a 2D beta field.
      !
      ! Note: A node must be surrounded by four floating or ocean cells to be considered
      !       a shelf/ocean node.  Nodes along the grounding line retain previous values of beta.

      !if (present(floating_mask) .and. present(ocean_mask)) then
      !   do ns = 1, nsn-1
      !      do ew = 1, ewn-1
      !         if ( (floating_mask(ew,ns  )  ==1 .or. ocean_mask(ew,ns)    ==1)   .and.  &
      !              (floating_mask(ew,ns+1)  ==1 .or. ocean_mask(ew,ns+1)  ==1)   .and.  &
      !              (floating_mask(ew+1,ns)  ==1 .or. ocean_mask(ew+1,ns)  ==1)   .and.  &
      !              (floating_mask(ew+1,ns+1)==1 .or. ocean_mask(ew+1,ns+1)==1) ) then
      !            beta(ew,ns) = 0.d0
      !         endif
      !      enddo
      !   enddo
      !endif

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      beta(:,:) = 1.d4        ! Pa yr/m

      !TODO - Change this loop to work in parallel (set beta on the global grid and scatter to local)
      do ns=5, nsn-5
      do ew=1, ewn-1
        beta(ew,ns) = 100.d0      ! Pa yr/m
      end do
      end do

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force beta to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option will provide the till yield stress as calculate from the basal processes
      !!! submodel. Currently, to enable sliding over plastic till, simple specify the value of "beta" as 
      !!! if it were the till yield stress (in units of Pascals).
      
      beta(:,:) = mintauf(:,:) &                                                         ! plastic yield stress (Pa)
                         / dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 + (smallnum)**2 )   ! velocity components (m/yr)

      !!! since beta is updated here, communicate that info to halos
      call staggered_parallel_halo(beta)

    case(HO_BABC_BETA_BWAT)  ! set value of beta as proportional to value of bwat                                         

      !NOTE: This parameterization has not been scientifically tested.
      !TODO - Test option HO_BABC_BETA_BWAT
      !       Where do these constants come from?
      C = 10.d0   ! Does this assume that bwat is in units of m or dimensionless?
      m = 1.d0

      allocate(unstagbeta(ewn,nsn))

      unstagbeta(:,:) = 200.d0   ! Pa yr/m
                                 ! This setting ensures that the parameterization does nothing.  Remove it?

      where ( bwat > 0.d0 .and. unstagbeta > 200.d0 )
          unstagbeta = C / ( bwat**m )
      endwhere

      ! average beta from unstag grid onto stag grid
      beta = 0.5d0 * ( unstagbeta(1:ewn-1,:) + unstagbeta(2:ewn,:) )
      beta = 0.5d0 * ( unstagbeta(:,1:nsn-1) + unstagbeta(:,2:nsn) )
   
      deallocate(unstagbeta) 

    !Note: This is redundant in that it could be implemented by using HO_BETA_CONST with beta_const = 1.d10
    !      But keeping it for historical reasons since many config files use it

    case(HO_BABC_LARGE_BETA)      ! frozen (u=v=0) ice-bed interface

      beta(:,:) = 1.d10           ! Pa yr/m

    case(HO_BABC_ISHOMC)          ! prescribe according to ISMIP-HOM test C

       !Note: Ideally, beta would be read in from an external netCDF file.
       !      However, this is not possible given that the global velocity grid is smaller
       !       than the ice grid and hence not able to fit the full beta field.
       !      The following code sets beta on the full grid as prescribed by Pattyn et al. (2008).
       !NOTE: This works only in serial!

       Ldomain = (ewn-2*nhalo) * dew   ! size of full domain (must be square)
       omega = 2.d0*pi / Ldomain

       ilo = nhalo
       ihi = ewn-nhalo
       jlo = nhalo
       jhi = nsn-nhalo
       
       ! Prescribe beta as in Pattyn et al., The Cryosphere, 2008
       beta(:,:) = 0.d0
       do j = jlo, jhi
          do i = ilo, ihi
             dx = dew * (i-ilo)
             dy = dns * (j-jlo)
             beta(i,j) = 1000.d0 + 1000.d0 * sin(omega*dx) * sin(omega*dy)
          enddo
       enddo

    case(HO_BABC_EXTERNAL_BETA)   ! use value passed in externally from CISM

      ! scale CISM input value to dimensional units of (Pa yr/m)

       ! beta is initialized to a negative value; we can use that fact to check whether
       ! it has been read correctly from the file
       if (maxval(beta) <= 0.d0) then
          call write_log('ERROR: Trying to use HO_BABC_EXTERNAL_BETA, but all beta values are <= 0,')
          call write_log('which implies that beta could not be read from the input file.')
          call write_log('Make sure that beta is in the cism input file,')
          call write_log('or change which_ho_babc to a different option.')
          call write_log('Invalid value for beta. See log file for details.', GM_FATAL)
       end if

!!      beta(:,:) = beta(:,:) * ( tau0 / vel0 / scyr )   ! already dimensional

      ! this is a check for NaNs, which indicate, and are replaced by no slip

      do ns=1, nsn-1
         do ew=1, ewn-1 
            if( beta(ew,ns) /= beta(ew,ns) )then
               beta(ew,ns) = 1.d10     ! Pa yr/m
            endif
         end do
      end do

    case(HO_BABC_POWERLAW)   ! A power law that uses effective pressure
       ! See Cuffey & Paterson, Physics of Glaciers, 4th Ed. (2010), p. 240, eq. 7.17
       ! This is based on Weertman's classic sliding relation (1957) augmented by the bed-separation index described by Bindschadler (1983)
       !   ub = k Taub^p N^-q
       ! rearranging for Taub gives:
       !   Taub = k^(-1/p) ub^(1/p) N^(q/p)

       ! p and q should be _positive_ exponents
       ! TODO: p, q could be turned into config parameters instead of hard-coded
       ! If p/=1, this is nonlinear in velocity
       ! Cuffey & Paterson recommend p=3 and q=1, and k dependent on thermal & mechanical properties of ice and inversely on bed roughness.   
       p = 3.0d0; q = 1.0d0

       beta = basal_physics%friction_powerlaw_k**(-1.0d0/p) * basal_physics%effecpress_stag**(q/p)    &
              * dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 )**(1.0d0/p-1.0d0)

    case(HO_BABC_COULOMB_FRICTION)

      ! Basal stress representation using coulomb friction law
      ! Coulomb sliding law: Schoof 2005 PRS, eqn. 6.2  (see also Pimentel, Flowers & Schoof 2010 JGR)

      ! Need flwa of the basal layer on the staggered grid
      where (thck > 0.0)
        imask = 1
      elsewhere
        imask = 0
      end where
      call glissade_stagger(ewn,        nsn,               &
                           flwa_basal,  flwa_basal_stag,   &
                           imask,       stagger_margin_in = 1)
      ! TODO Not sure if a halo update is needed on flwa_basal_stag!  I don't think so if nhalo>=2.

      ! Setup parameters needed for the friction law
      m_max = basal_physics%Coulomb_Bump_max_slope  !maximum bed obstacle slope(unitless)
      lambda_max = basal_physics%Coulomb_bump_wavelength ! wavelength of bedrock bumps (m)
      ! biglambda = wavelength of bedrock bumps [m] * flwa [Pa^-n yr^-1] / max bed obstacle slope [dimensionless]
      big_lambda = lambda_max / m_max * flwa_basal_stag
      Coulomb_C = basal_physics%Coulomb_C    ! Basal shear stress factor (Pa (m^-1 y)^1/3)
      !gn                         ! Glen's flaw law from parameter module

      beta = Coulomb_C * basal_physics%effecpress_stag * &
             (dsqrt(thisvel**2 + othervel**2 + smallnum**2))**(1.0d0/gn - 1.0d0) * &
             (                                                                     &
              dsqrt(thisvel**2 + othervel**2 + smallnum**2) +                      &
             basal_physics%effecpress_stag**gn * big_lambda                        &
             )**(-1.0d0/gn)

      ! for numerical stability purposes
      where (beta>1.0d8)
              beta = 1.0d8
      end where

    case default
       ! do nothing

   end select

   ! check for areas where ice is floating and make sure beta in these regions is 0 
   do ns=1, nsn-1
     do ew=1, ewn-1 
       if( GLIDE_IS_FLOAT( mask(ew,ns) ) )then
         beta(ew,ns) = 0.d0
       endif
     end do
   end do

 end subroutine calcbeta

!***********************************************************************

end module glissade_basal_traction

!***********************************************************************
