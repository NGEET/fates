module LUsolveMod
!
        !USES:
        use FatesConstantsMod, only : r8 => fates_r8
        use FatesHydraulicsMemMod, only: num_nodes
        use FatesGlobals, only      : fates_log
        use FatesGlobals, only      : endrun => fates_endrun
        use shr_log_mod , only      : errMsg => shr_log_errMsg

        implicit none

        character(len=*), parameter, private :: sourcefile = &
         __FILE__

        public :: ludcmp
        public :: lubksb

      contains

      subroutine ludcmp(a,n,indx,d)
!     From numerical recipe
      implicit none
      real(r8), dimension(n,n), intent(inout) :: a
      integer, dimension(n), intent(out) :: indx
      real(r8), intent(out) :: d
!     Given an N×N input matrix a, this routine replaces it by the LU decomposition of a
!     rowwise permutation of itself. On output, a is arranged as in equation (2.3.14);
!     indx is an output vector of length N that records the row permutation effected by the partial pivoting;
!     d is output as ± 1 depending on whether the number of row interchanges was even or odd, 
!     respectively. This routine is used in combination with lubksb to solve linear equations or
!     invert a matrix.  REAL(SP), DIMENSION(size(a,1)) :: vv         
!     vv stores the implicit scaling of each row.
      real(r8), parameter :: tiny=1.0e-20_r8 !A small number.
      real(r8) :: aamax, sum, dum
      real(r8) :: vv(n)
      integer :: j,n,imax,k,i
      d=1.0_r8
!     No row interchanges yet.
      vv=maxval(abs(a),dim=2)
!     Loop over rows to get the implicit scaling information.
      if (any(vv == 0.0)) then
        write(fates_log(),*) 'singular matrix in ludcmp'  !There is a row of zeros.
        call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      vv=1.0_r8/vv
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax) then
          do k=1,N
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j).eq.0.) a(j,j)=tiny
        if (j.ne.n) then
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
!
      end subroutine ludcmp

      subroutine lubksb(a,n,indx,b)
      implicit none
      real(r8), dimension(n,n), intent(in) :: a
      integer     , dimension(n), intent(in) :: indx
      real(r8), dimension(n), intent(inout) :: b
!     Solves the set of N linear equations A·X=B. Here the N×N matrix a is input, not as the original matrix A,
!     but rather as its LU decomposition, determined by the routine ludcmp.
!     indx is input as the permutation vector of length N returned by ludcmp.
!     b is input as the right-hand-side vector B,also of length N, and returns with the solution vector X.
!     a and indx are not modified by this routine and can be left in place for successive calls 
!     with different right-hand sides b. This routine takes into account the possibility that b
!     will begin with many zero elements, so it is efficient for use in matrix inversion.
      integer      :: i,n,ii,ll
      real(r8) :: summ
      ii=0
!     When ii is set to a positive value, it will become the in-dex of the first nonvanishing element of b.
!     Wenowdo the forward substitution, equation (2.3.6). The only new
!     wrinkle is to unscramble the permutation as we go.
      do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0) then
          ii=i
!A nonzero element was encountered, so from now on we will have to do the dot product above.
        end if
        b(i)=summ
     end do
     do i=n,1,-1
!    Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
     end do
     end subroutine lubksb
!
end module LUsolveMod
