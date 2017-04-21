!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module glc_indexing
  
  !BOP
  ! !MODULE: glc_indexing

  ! !DESCRIPTION:
  ! Contains information about the indexing of the points owned by each processor.
  !
  ! This includes local indices (translation between (i,j) and a scalar 1..n) and global
  ! indices (unique indices across all procs).
  !
  ! Also contains subroutines for translating arrays between (i,j) and scalar 1..n

#include "shr_assert.h"
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_kind_mod, only: R8=>SHR_KIND_R8

  implicit none
  private
  save

  ! !PUBLIC ROUTINES:
  public :: glc_indexing_init
  public :: local_to_global_indices
  public :: vector_to_spatial
  public :: spatial_to_vector

  ! !PUBLIC MODULE VARIABLES:
  integer, public :: nx       ! number of columns owned by this proc
  integer, public :: ny       ! number of rows owned by this proc
  integer, public :: npts     ! total number of points owned by this proc
  integer, public :: nx_tot   ! total number of columns in full grid (all procs)
  integer, public :: ny_tot   ! total number of rows in full grid (all procs)
  integer, public :: npts_tot ! total number of points in full grid (all procs)

  integer, allocatable, private :: local_indices(:,:)  ! mapping from (i,j) to 1..npts
  integer, allocatable, private :: global_indices(:,:) ! unique indices across all procs (matches indexing on mapping files)
  
contains

  !-----------------------------------------------------------------------
  subroutine glc_indexing_init(params, instance_index)
    !
    ! !DESCRIPTION:
    ! Initialize indices stored here.
    !
    ! Note that the global indexing needs to match the indexing on the SCRIP grid file
    ! that is used to generate GLC mapping files for the coupler.
    !
    ! !USES:
    use glad_main, only : glad_params, glad_get_grid_size, glad_get_grid_indices
    !
    ! !ARGUMENTS:
    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet index
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'glc_indexing_init'
    !-----------------------------------------------------------------------

    call glad_get_grid_size(params, instance_index, &
         ewn = nx, nsn = ny, npts = npts, &
         ewn_tot = nx_tot, nsn_tot = ny_tot, npts_tot = npts_tot)
    
    allocate(local_indices(nx, ny))
    allocate(global_indices(nx, ny))
    
    call glad_get_grid_indices(params, instance_index, global_indices, local_indices)
    
  end subroutine glc_indexing_init

  !-----------------------------------------------------------------------
  function local_to_global_indices() result(gindex)
    !
    ! !DESCRIPTION:
    ! Returns an array that maps local indices to global indices
    !
    ! gindex(n) gives the global index corresponding to local index n
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, allocatable :: gindex(:)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: i, j, n

    character(len=*), parameter :: subname = 'local_to_global_indices'
    !-----------------------------------------------------------------------

    allocate(gindex(npts))
    do j = 1,ny
       do i = 1,nx
          n = local_indices(i,j)
          gindex(n) = global_indices(i,j)
       end do
    end do

  end function local_to_global_indices

  !-----------------------------------------------------------------------
  subroutine vector_to_spatial(arr_vector, arr_spatial)
    !
    ! !DESCRIPTION:
    ! Convert a vector array (1..n) to a spatial array (i,j)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: arr_vector(:)
    real(r8), intent(out) :: arr_spatial(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: i, j, n

    character(len=*), parameter :: subname = 'vector_to_spatial'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(arr_vector) == (/npts/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(arr_spatial) == (/nx, ny/)), errMsg(__FILE__, __LINE__))

    do j = 1, ny
       do i = 1, nx
          n = local_indices(i,j)
          arr_spatial(i,j) = arr_vector(n)
       end do
    end do

  end subroutine vector_to_spatial

  !-----------------------------------------------------------------------
  subroutine spatial_to_vector(arr_spatial, arr_vector)
    !
    ! !DESCRIPTION:
    ! Convert a spatial array (i,j) to a vector array (1..n)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: arr_spatial(:,:)
    real(r8), intent(out) :: arr_vector(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i, j, n

    character(len=*), parameter :: subname = 'spatial_to_vector'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(arr_spatial) == (/nx, ny/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(arr_vector) == (/npts/)), errMsg(__FILE__, __LINE__))

    do j = 1, ny
       do i = 1, nx
          n = local_indices(i,j)
          arr_vector(n) = arr_spatial(i,j)
       end do
    end do

  end subroutine spatial_to_vector


end module glc_indexing
