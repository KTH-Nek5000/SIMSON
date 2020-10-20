c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine curvefit(x,y,np,n,m,intx,inty,intdy,int,a)
c
c     For more information see numerical recipes "svdfit"
c
      implicit none

      integer nmax,mmax
      parameter(nmax=2000,mmax=10)
      integer n,m,int,np,i,j
      real x(np),y(np),a(nmax),sv(mmax,mmax),su(nmax,mmax),sw(mmax)
      real wmax,thresh,b(nmax)
      real intx(int),inty(int),intdy(int)

      do i=1,n
         do j=1,m
            su(i,j)=x(i)**(j-1)
         end do
         b(i)=y(i)
      end do

      call svdcmp(su,n,m,nmax,mmax,sw,sv)

      wmax=0.0
      do j=1,m
         if (sw(j).gt.wmax) wmax=sw(j)
      end do
      thresh=1e-14*wmax
      do j=1,m
         if (sw(j).lt.thresh) sw(j)=0.0
      end do

      call svbksb(su,sw,sv,n,m,nmax,mmax,b,a)

      do i=1,int
         inty(i)=0.0
         intdy(i)=0.0
         do j=1,m
            inty(i)=inty(i)+a(j)*intx(i)**(j-1)
            intdy(i)=intdy(i)+(j-1)*a(j)*intx(i)**(j-2)
         end do
      end do

      return

      end subroutine curvefit
