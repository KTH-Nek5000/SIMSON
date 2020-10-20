c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine spline(x,y,y2,n,xex,yex,nex)
c
c     Does spline interpolation and evaluates all points
c     with zero derivative
c
      implicit none

      include 'par.f'

      integer i,nex,n
      real x(n),y(n),y2(n),xex(n),yex(n)
      real gam(nx+nyp),sig,p,t,s,a,b

      do i=1,n
         xex(i)=0
         yex(i)=0
      end do
c
c     Make spline
c
      y2(1)=0
      gam(1)=0

      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2
         y2(i)=(sig-1)/p
         gam(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*gam(i-1))/p
      end do
      y2(n)=0

      do i=n-1,1,-1
         y2(i)=y2(i)*y2(i+1)+gam(i)
      end do
c
c     Find points with zero derivative
c
      nex=0
      do i=1,n-1
         a=-y2(i)/(y2(i+1)-y2(i))
         b=(a**2-(6*(y(i+1)-y(i))/(x(i+1)-x(i))**2-y2(i+1)-2*y2(i))
     &        /(3*(y2(i+1)-y2(i))))
         if (b.gt.0) then
            b=sqrt(b)
            t=a+b
            s=a-b
            if (t.lt.1.and.t.gt.0.and.(s.lt.0.or.s.gt.1)) then
               nex=nex+1
               xex(nex)=x(i)+t*(x(i+1)-x(i))
               yex(nex)=t*y(i+1)+(1-t)*y(i)+(x(i+1)-x(i))**2*
     &              ((t**3-t)*y2(i+1)+((1-t)**3-(1-t))*y2(i))/6
            end if

            if (s.lt.1.and.s.gt.0.and.(t.lt.0.or.t.gt.1)) then
               t=s
               nex=nex+1
               xex(nex)=x(i)+t*(x(i+1)-x(i))
               yex(nex)=t*y(i+1)+(1-t)*y(i)+(x(i+1)-x(i))**2*
     &              ((t**3-t)*y2(i+1)+((1-t)**3-(1-t))*y2(i))/6
            end if
         end if
      end do

      return

      end subroutine spline
