      module petscMod
!
#include <petsc/finclude/petsc.h>
        !USES:
        use FatesConstantsMod, only : r8 => fates_r8
        use FatesHydraulicsMemMod, only: num_nodes
        use FatesHydraulicsMemMod, only: num_connections
        use FatesHydraulicsMemMod, only: ed_cohort_hydr_type
        use petscsys
        use petscvec
        use petscmat
        use petscksp
        use petscpc

        implicit none

        private

        Mat, public :: fmat
        Vec, public :: fsol_vec
        Vec, public :: frhs_vec
        KSP, public :: fksp

        integer, allocatable, public :: petsc_offset(:)

        public :: petsc_solver_init, &
                 petsc_put_rhs, &
                 petsc_get_solution, &
                 petsc_solve, &
                 petsc_solver_destroy 

      contains

      subroutine petsc_solver_init(myksp,mat,rhs_vec,sol_vec,conn_dn,conn_up)
!
!
      KSP :: myksp
      Mat :: mat
      Vec :: rhs_vec
      Vec :: sol_vec
      PC :: mypc
!
      integer :: conn_dn(num_connections)
      integer :: conn_up(num_connections)
      integer, allocatable :: d_nnz(:)
      PetscErrorCode :: ierr
      integer :: n, nn, icnx, id_dn, id_up
      

      n = num_nodes
!
! --- create petsc solution and right-hand-side vectors
!
      call VecCreateSeq(PETSC_COMM_SELF,num_nodes,sol_vec,ierr)
      call VecDuplicate(sol_vec,rhs_vec,ierr)
!
! --- petsc vectors are zero-based, create offset array
!
      if( .not.allocated(petsc_offset) ) then
        allocate(petsc_offset(num_nodes))
      endif
      do nn=1,num_nodes
        petsc_offset(nn) = nn-1
      enddo
!
! --- calculate number of nonzeros in each row of jacobian
!
      allocate(d_nnz(num_nodes))
      d_nnz(:) = 1
      do icnx = 1,num_connections
        id_dn = conn_dn(icnx)
        id_up = conn_up(icnx)
        d_nnz(id_dn) = d_nnz(id_dn) + 1
        d_nnz(id_up) = d_nnz(id_up) + 1
      enddo
!
! --- create petsc jacobian matrix
!
        call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,0,d_nnz,mat,ierr)
        deallocate(d_nnz)
!
! --- pick petsc linear solver and set options
!
        call KSPCreate(PETSC_COMM_SELF,myksp,ierr)
        call KSPSetFromOptions(myksp,ierr)
        call KSPSetTolerances(myksp,1.d-6, &
                             PETSC_DEFAULT_REAL, &
                             PETSC_DEFAULT_REAL, &
                             PETSC_DEFAULT_INTEGER,ierr)
        call KSPSetType(myksp,KSPBCGS,ierr)
        call KSPGetPC(myksp,mypc,ierr)
!
! --- pick petsc preconditioner and set options
!
        call PCFactorSetZeroPivot(mypc,1.d-20,ierr)
        call PCSetFromOptions(mypc,ierr)
        call PCSetType(mypc,PCILU,ierr)

      end subroutine petsc_solver_init

      subroutine petsc_put_rhs(array,rhs_vec)
      
        implicit none

        real*8 :: array(*)
        Vec :: rhs_vec
        integer :: n
        real*8, pointer :: vec_ptr(:)
        PetscErrorCode :: ierr

        call VecGetSize(rhs_vec,n,ierr)
        call VecGetArrayF90(rhs_vec,vec_ptr,ierr)
        vec_ptr(1:n) = array(1:n)
        call VecRestoreArrayF90(rhs_vec,vec_ptr,ierr)

      end subroutine petsc_put_rhs

      subroutine petsc_get_solution(array,sol_vec)
      
        implicit none

        real*8 :: array(*)
        Vec :: sol_vec
        integer :: n
        real*8, pointer :: vec_ptr(:) 
        PetscErrorCode :: ierr

        call VecGetSize(sol_vec,n,ierr)
        call VecGetArrayF90(sol_vec,vec_ptr,ierr)
        array(1:n) = vec_ptr(1:n)
        call VecRestoreArrayF90(sol_vec,vec_ptr,ierr)

      end subroutine petsc_get_solution

      subroutine petsc_solve(myksp,mat,rhs_vec,sol_vec)

        implicit none

        KSP :: myksp
        Mat :: mat
        Vec :: rhs_vec
        Vec :: sol_vec
        PetscErrorCode :: ierr
        KSPConvergedReason :: reason_ksp
        PetscViewer :: PV_VIEWER,M_VIEWER,SV_VIEWER,KSP_VIEWER 

#if 1
!        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"prob_vec.dat",PV_VIEWER,IERR )
!        CALL VecView( rhs_vec,PV_VIEWER,IERR )
!        CALL PetscViewerDestroy( PV_VIEWER,IERR )
#endif
        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        if( ierr.ne.0 ) then
          write(*,*) 'petsc: matassemblybegin'
          return
        endif
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
        if( ierr.ne.0 ) then
          write(*,*) 'petsc: matassemblyend'
          return
        endif

!        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"matrix.dat",M_VIEWER,IERR )
!        CALL MatView( mat,M_VIEWER,IERR )
!        CALL PetscViewerDestroy( M_VIEWER,IERR )


        call KSPSetOperators(myksp,mat,mat,ierr)
        if( ierr.ne.0 ) then
          write(*,*) 'petsc: kspsetoperators'
          return
        endif
        call KSPSolve(myksp,rhs_vec,sol_vec,ierr) 
        if( ierr.ne.0 ) then
          write(*,*) 'petsc: kspsolve'
          return
        endif
        call KSPGetConvergedReason(myksp,reason_ksp,ierr)
        if( ierr.ne.0 ) then
          write(*,*) 'petsc: kspgetconvergedreason'
          return
        endif
!        CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,"solu_vec.dat",SV_VIEWER,IERR )
!        CALL VecView( sol_vec,SV_VIEWER,IERR )
!        CALL PetscViewerDestroy( SV_VIEWER,IERR )


      end subroutine petsc_solve

      subroutine petsc_solver_destroy(myksp,mat,rhs_vec,sol_vec)

        implicit none

        KSP :: myksp
        Mat :: mat
        Vec :: rhs_vec
        Vec :: sol_vec
        PetscErrorCode :: ierr

        call KSPDestroy(myksp,ierr)
        call MatDestroy(mat,ierr)
        call VecDestroy(sol_vec,ierr)
        call VecDestroy(rhs_vec,ierr)

        if (allocated(petsc_offset)) deallocate(petsc_offset)

      end subroutine petsc_solver_destroy

      end module petscMod

