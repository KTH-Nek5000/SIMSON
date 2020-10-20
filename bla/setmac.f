c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine setmac(q,c,d,e,hlam,n,m,md,neumann,cc)
c
c     Sets up the matrix system to solve a Helmholtz equation
c     by the Chebyshev tau method.
c     d2u/dx2-hlam*u=f is solved in Chebyshev space for n coefficients of u
c     m matrices of size n by n are created, one for each lambda in hlam
c
c     Used by the pressure solver.
c
      implicit none

      integer neumann
      integer n,m,md,i,j,k
      real q(md,n+2),c(md,n),d(md,n),e(md,n),hlam(m),rk1,rk2,rk3
      real cc

c
c     Setup full rows: equations 1,2
c
      if (neumann.eq.1) then
c
c     Neumann bc for pressure.(Dp=1/re*D2v)
c
         do i=1,n
            do j=1,m
               q(j,i) = real(i-1)**2
            end do
         end do
      else if (neumann.eq.2) then
c
c     Mixed bc for electric potential with conducting walls
c
         do i=1,n
            do j=1,m
               q(j,i) = real(i-1)**2 + cc*hlam(j)
            end do
         end do
      else
c
c     Dirichlet bc.(p=-i*(alpha*D2u+beta*D2w)/Re/(alpha**2+beta**2))
c
         do i=1,n
            do j=1,m
               q(j,i) = 1.
            end do
         end do
      end if
c
c     Setup the three diagonals: equations 1,...,n-6
c
      do i=1,n-6
         k=i+1
         rk1 = .25/(real(k  )*real(k-1))
         rk2 = .5 /(real(k-1)*real(k+1))
         rk3 = .25/(real(k  )*real(k+1))
         do j=1,m
            c(j,i) =   -rk1*hlam(j)
            d(j,i) = 1.+rk2*hlam(j)
            e(j,i) =   -rk3*hlam(j)
         end do
      end do
      do j=1,m
         c(j,1)=2.*c(j,1)
      end do
c
c     Equations n-5,n-4
c
      do i=n-5,n-4
         k=i+1
         rk1 = .25/(real(k  )*real(k-1))
         rk2 = .5 /(real(k-1)*real(k+1))
         do j=1,m
            c(j,i) =   -rk1*hlam(j)
            d(j,i) = 1.+rk2*hlam(j)
            e(j,i) = 0.
         end do
      end do
c
c     Equations n-3,n-2 (these may be not diagonally dominant for a small n)
c
      do i=n-3,n-2
         k=i+1
         rk1 = .25/(real(k)*real(k-1))
         do j=1,m
            c(j,i) = -rk1*hlam(j)
            d(j,i) = 1.
            e(j,i) = 0.
         end do
      end do

      end subroutine setmac
