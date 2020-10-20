c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine spline(x,y,n,yp1,ypn,y2,u)
c
c     Computes a natural spline
c
      implicit none

      integer n
      real x(n),y(n),y2(n),yp1,ypn,u(n)
c
c     Local variables
c
      integer i,k
      real p,qn,sig,un

      if (yp1.gt..99e30) then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1))-yp1))
      end if
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99e30) then
         qn=0.
         un=0.
      else
         qn=0.5
         un=(3./(x(n)-x(n-1))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))))
      end if
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      end subroutine spline
