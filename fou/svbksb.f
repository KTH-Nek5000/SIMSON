c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
c
c     For more information see numerical recipes pg 57
c
      implicit none

      integer i,j,m,n,jj,mp,np,nmax
      real s
      parameter (nmax=10)
      real u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(nmax)

      do j=1,n
         s=0.0
         if (w(j).ne.0.) then
            do i=1,m
               s=s+u(i,j)*b(i)
            end do
            s=s/w(j)
         end if
         tmp(j)=s
      end do
      do j=1,n
         s=0.0
         do jj=1,n
            s=s+v(j,jj)*tmp(jj)
         end do
         x(j)=s
      end do

      return

      end subroutine svbksb
