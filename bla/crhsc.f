c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine crhsc(f,bc,m,n,nrhs,md,nd)
c
c     Sets up the right handed sides for the Chebyshev tau method.
c     d2u/dx2-hlam*u=f is solved in terms of u in the Chebyshev space
c
c     Used by the pressure solver.
c
      implicit none
c
      integer m,md,n,nd,nrhs,i,j,k,y,yy
      real f(md,nd,nrhs),bc(md,2,nrhs),rk1,rk2,rk3,rk4,rk5,rk6,tm1,tm2
c
c     For j=3
c
      do k=1,nrhs
         do i=1,m
            tm1=f(i,2,k)
            f(i,2,k)=f(i,3,k)
            f(i,3,k)=.25*f(i,1,k)+f(i,5,k)/24.-f(i,3,k)/6.
            f(i,1,k)=tm1
         end do
      end do
c
c     For j=4,...,n-4
c
      do k=1,nrhs
         do j=4,n-4,2
            y=j-1
            yy=j
            rk1=.25/(real(y)   *real(y-1))
            rk2=.5 /(real(y-1) *real(y+1))
            rk3=.25/(real(y)   *real(y+1))
            rk4=.25/(real(yy)  *real(yy-1))
            rk5=.5 /(real(yy-1)*real(yy+1))
            rk6=.25/(real(yy)  *real(yy+1))
            do i=1,m
               tm1=f(i,j,k)
               tm2=f(i,j+1,k)
               f(i,j  ,k) = rk1*f(i,1,k)+rk3*f(i,j+2,k)-rk2*f(i,j,k)
               f(i,j+1,k) = rk4*f(i,2,k)+rk6*f(i,j+3,k)-rk5*f(i,j+1,k)
               f(i,1  ,k) = tm1
               f(i,2  ,k) = tm2
            end do
         end do
      end do
      do k=1,nrhs
         do i=1,m
c
c     For j=n-3,...,n
c
            tm1=f(i,n-3,k)
            tm2=f(i,n-2,k)
            j=n-4
            f(i,n-3,k)=.25/(real(j)*real(j-1))*f(i,1,k)-
     &           .5/(real(j-1)*real(j+1))*f(i,n-3,k)
            j=n-3
            f(i,n-2,k)=.25/(real(j)*real(j-1))*f(i,2,k)-
     &           .5/(real(j-1)*real(j+1))*f(i,n-2,k)
            f(i,n-1,k)=.25/(real(n-2)*real(n-3))*tm1
            f(i,n,k)=.25/(real(n-1)*real(n-2))*tm2
c
c     For j=1,2(boundary conditions for symmetric and anti-symmetric parts.)
c
            f(i,1,k) = bc(i,1,k)
            f(i,2,k) = bc(i,2,k)
         end do
      end do

      end subroutine crhsc
