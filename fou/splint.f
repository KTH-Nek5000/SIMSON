c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine splint(x,y,y2,n,intx,inty,intdy,intn)

      implicit none

      integer i,n,intn,k
      real x(n),y(n),y2(n),intx(n),inty(n),intdy(n)
      real h,a,b

      k=1
      do i=1,intn
 1       if (intx(i).gt.x(k+1)) then
            k=k+1
            goto 1
         end if
         h=x(k+1)-x(k)
         a=(x(k+1)-intx(i))/h
         b=(intx(i)-x(k))/h
         inty(i)=a*y(k)+b*y(k+1)+
     &        ((a**3-a)*y2(k)+(b**3-b)*y2(k+1))*h**2/6
         intdy(i)=(y(k+1)-y(k))/h+
     &        ((3*b**2-1)*y2(k+1)-(3*a**2-1)*y2(k))*h/6
      end do

      return

      end subroutine splint
