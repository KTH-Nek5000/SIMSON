c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine dcheb(b,a,n,m,md)
c
c     Calculates the first derivative in chebyshev space
c
c     a is the input array, the result is returned in b
c     n number of chebyshev coefficients
c     m number of derivatives
c     md first index length on data array
c
      implicit none

      integer n,m,md
      real a(md,n),b(md,n)

      integer i,k
c
      do k=1,m
         b(k,n-2)=2.*float(n-2)*a(k,n-1)
         b(k,n-1)=2.*float(n-1)*a(k,n)
      end do
      do i=n-3,1,-1
         do k=1,m
            b(k,i)=b(k,i+2)+2.*float(i)*a(k,i+1)
         end do
      end do
      do k=1,m
         b(k,n)=0.0
         b(k,1)=.5*b(k,1)
      end do

      return

      end subroutine dcheb
