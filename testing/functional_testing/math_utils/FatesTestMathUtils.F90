program FatesTestQuadSolvers

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesUtilsMod,     only : QuadraticRootsNSWC, QuadraticRootsSridharachary
  use FatesUtilsMod,     only : GetNeighborDistance

  implicit none

  ! CONSTANTS:
  integer, parameter          :: n = 4                    ! number of points to test
  character(len=*), parameter :: out_file = 'quad_out.nc' ! output file

  ! LOCALS:
  integer  :: i                ! looping index
  real(r8) :: a(n), b(n), c(n) ! coefficients for quadratic solvers
  real(r8) :: root1(n)         ! real part of first root of quadratic solver
  real(r8) :: root2(n)         ! real part of second root of quadratic solver
  logical  :: err              ! whether quadratic solver encountered an error

  interface

    subroutine WriteQuadData(out_file, n, a, b, c, root1, root2)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file
      integer,          intent(in) :: n
      real(r8),         intent(in) :: a(:)
      real(r8),         intent(in) :: b(:)
      real(r8),         intent(in) :: c(:)
      real(r8),         intent(in) :: root1(:)
      real(r8),         intent(in) :: root2(:)
    end subroutine WriteQuadData

  end interface

  a = (/1.0_r8, 1.0_r8, 5.0_r8, 1.5_r8/)
  b = (/-2.0_r8, 7.0_r8, 10.0_r8, 3.2_r8/)
  c = (/1.0_r8, 12.0_r8, 3.0_r8, 1.1_r8/)

  do i = 1, n
    call QuadraticRootsNSWC(a(i), b(i), c(i), root1(i), root2(i), err)
  end do

  call WriteQuadData(out_file, n, a, b, c, root1, root2)

end program FatesTestQuadSolvers

! ----------------------------------------------------------------------------------------

subroutine WriteQuadData(out_file, n, a, b, c, root1, root2)
  !
  ! DESCRIPTION:
  ! Writes out data from the quadratic solver test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file ! output file name
  integer,          intent(in) :: n        ! number of points to write out
  real(r8),         intent(in) :: a(:)     ! coefficient a
  real(r8),         intent(in) :: b(:)     ! coefficient b
  real(r8),         intent(in) :: c(:)     ! coefficient c
  real(r8),         intent(in) :: root1(:) ! root1 from quadratic solver
  real(r8),         intent(in) :: root2(:) ! root2 from quadratic solver

  ! LOCALS:
  integer          :: n_index(n)   ! array of pft indices to write out
  integer          :: i            ! looping index
  integer          :: ncid         ! netcdf file id
  character(len=8) :: dim_names(1) ! dimension names
  integer          :: dimIDs(1)    ! dimension IDs
  integer          :: aID, bID, cID
  integer          :: root1ID, root2ID

  ! make index
  do i = 1, n
    n_index(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'n']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/n/), 1, dimIDs)

  ! register a
  call RegisterVar(ncid, 'a', dimIDs(1:1), type_double,  &
    [character(len=20)  :: 'units', 'long_name'],                 &
    [character(len=150) :: '', 'coefficient a'], 2, aID)

  ! register b
    call RegisterVar(ncid, 'b', dimIDs(1:1), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'coefficient b'], 2, bID)

  ! register c
    call RegisterVar(ncid, 'c', dimIDs(1:1), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'coefficient c'], 2, cID)

  ! register root1
    call RegisterVar(ncid, 'root1', dimIDs(1:1), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'root 1'], 2, root1ID)

  ! register root2
    call RegisterVar(ncid, 'root2', dimIDs(1:1), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'root 2'], 2, root2ID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, aID, a(:))
  call WriteVar(ncid, bID, b(:))
  call WriteVar(ncid, cID, c(:))
  call WriteVar(ncid, root1ID, root1(:))
  call WriteVar(ncid, root2ID, root2(:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteQuadData
