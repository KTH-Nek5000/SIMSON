c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine setmat(q,c,d,e,hlam,n,m,md,cim)
c
c     Sets up the system matrix to solve a helmholz equation
c
c     d2u/dx2-hlam*u=f is solved in Chebyshev space
c     If cim: for n coefficients of d2u/dx2 and 2 integration constants,
c     which means that we are solving n+2 equations
c     If tau: for n coefficients of u
c     The first two equations are used for boundary conditions.
c     m matrices of size n by n are created, one for each lambda in hlam
c     q hold the elements of the top rows and c, d and e those of the
c     three diagonals exept elements of the first 2 rows se trid.f
c     In case of tau the first 2 elements of q,c,d and e are not used
c
      implicit none

      logical cim
      integer n,m,md
      real q(md,n+2),c(md,n),d(md,n),e(md,n),hlam(m)

      integer i,j,k
      real rk1,rk2,rk3

      call stopnow(434342)

c
c     Equation 5,...,n for cim and equation 3,...,n-4 for tau.
c
      do i=3,n
         k=i-1
         rk1=-1./real(4*k*(k-1))
         rk2=1./real(2*(k-1)*(k+1))
         rk3=-1./real(4*k*(k+1))
         do j=1,m
            c(j,i)=rk1*hlam(j)
            d(j,i)=1.+rk2*hlam(j)
            e(j,i)=rk3*hlam(j)
         end do
      end do
      do j=1,m
         c(j,3)=2.*c(j,3)
      end do

      if (cim) then
c
c     Setup full rows (equation 1,2) and equation 3, 4
c     and d for equation n-1 and n
c
         do j=1,m
            k=n-2
            d(j,n-1)=1.+hlam(j)*real(2*k)/real(4*k*(k-1)*(k+1))
            k=n-1
            d(j,n)=1.+hlam(j)*real(k+1)/real(4*k*(k-1)*(k+1))
            q(j,1)=1.
            q(j,2)=1.
            q(j,3)=.25
            q(j,4)=-1./12.
            q(j,5)=-7./48.
            q(j,n+1)=-real(n-7)/(real(4*n-4)*real(n-3)*real(n-4))
            q(j,n+2)=.5/(real(n-1)*real(n-2)*real(n-3))
            do i=6,n
               q(j,i)=3./(real(i-5)*real(i-4)*real(i-2)*real(i-1))
            end do
            c(j,1)=-hlam(j)
            d(j,1)=1.
            e(j,1)=0.
            c(j,2)=-hlam(j)
            d(j,2)=1.+hlam(j)*.125
            e(j,2)=-hlam(j)*.125
         end do
      else
         do j=1,m
            d(j,n-1)=1.
            d(j,n)=1.
            e(j,n-3)=0.
            e(j,n-2)=0.
            do i=1,n+2
               q(j,i)=1.
            end do
         end do
      end if

      end subroutine setmat
