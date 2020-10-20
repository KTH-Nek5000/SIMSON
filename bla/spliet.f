c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine spliet(xa,ya,y2a,n,x,y)
c
c     Evaluates natural spline interpolation and extrapolation
c     see also spline.f
c
      implicit none

      integer n
      real xa(n),ya(n),y2a(n),x,y

      integer klo,khi,k
      real a,b,h,yp

      if (n.eq.1) then
         y=ya(1)
         return
      end if
      klo=1
      khi=n

 1    if (khi-klo .gt. 1) then
         k=(khi+klo)/2
         if (xa(k) .gt. x) then
            khi=k
         else
            klo=k
         end if
         goto 1
      end if

      h=xa(khi)-xa(klo)
      if (h .eq. 0.) then
         write(*,*) 'Bad xa input, in routine splint'
         stop
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      if (x.lt.xa(klo)) then
c
c     yprime evaluated at x=xa(klo)
c
         yp=-1./h*ya(klo)+1./h*ya(khi)+
     &        ((2.)*(-1./h)*y2a(klo)+(-1.)*(1./h)*y2a(khi))*
     &        (h**2)/6.
         y=ya(1)+(x-xa(1))*yp
      end if
      if (x.gt.xa(khi)) then
c
c     yprime evaluated at x=xa(khi)
c
         yp=-1./h*ya(klo)+1./h*ya(khi)+
     &        ((-1.)*(-1./h)*y2a(klo)+(2.)*(1./h)*y2a(khi))*
     &        (h**2)/6.
         y=ya(n)+(x-xa(n))*yp
      end if

      end subroutine spliet
