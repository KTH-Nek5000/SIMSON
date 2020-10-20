c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine setmacrhs(rhs,f,c,d,e,hlam,m,n,nrhs,md,nd)
c
c     Sets up the matrix system and the right handed sides to solve
c     a Helmholtz equation by the Chebyshev tau method.
c
c     d2u/dx2-hlam*u=f is solved in Chebyshev space for n coefficients of u
c     m matrices of size n by n are created, one for each lambda in hlam
c
      implicit none

      integer n,m,md,nd,nrhs,i,j,k
      real rhs(md,nd,nrhs),f(md,nd,nrhs),c(md,n),d(md,n),e(md,n),hlam(m)
      real di,rk1,rk2,rk3,rk4,rk5

      do j=1,m
         c(j,1)=-.25*hlam(j)
      end do
      di=1.
      do i=1,n-3
         rk2=.5/(di*(di+2.))
         rk3=-.25/((di+1.)*(di+2.))
         di=di+1.
         do j=1,m
            c(j,i+1)=rk3*hlam(j)
            d(j,i)=1.+rk2*hlam(j)
            e(j,i)=rk3*hlam(j)
         end do
      end do
c
c     For i=4,...,n-4
c
      do k=1,nrhs
         di=2.
         rk3=1./24.
         do i=2,n-6
            rk1=rk3
            rk2=-.5/(di*(di+2.))
            rk3=.25/((di+1.)*(di+2.))
            di=di+1.
            do j=1,m
               rhs(j,i+2,k)=rk1*f(j,i,k)+rk2*f(j,i+2,k)+rk3*f(j,i+4,k)
            end do
         end do
      end do

      rk1=.25/((di+1.)*(di+2.))
      rk2=-.5/(di*(di+2.))
      rk4=-.5/((di+1.)*(di+3.))
      rk5=.25/((di+2.)*(di+3.))
      di=.25/((di+3.)*(di+4.))

      do k=1,nrhs
         do j=1,m
c
c     i=3
c
            rhs(j,3,k)=.25*f(j,1,k)+f(j,5,k)/24.-f(j,3,k)/6.
c
c     For i=n-3,...,n
c
            rhs(j,n-3,k)=rk3*f(j,n-5,k)+rk2*f(j,n-3,k)
            rhs(j,n-2,k)=rk1*f(j,n-4,k)+rk4*f(j,n-2,k)
            rhs(j,n-1,k)=rk5*f(j,n-3,k)
            rhs(j,n,k)=di*f(j,n-2,k)
         end do
      end do

      end subroutine setmacrhs
