c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine d2cheb(b,a,n,m,md)
c
c     Calculates the second derivative in Chebyshev space
c
c     a is the input array, the result is returned in b
c     n number of chebyshev coefficients
c     m number of derivatives
c     md first index length on data array
c
      implicit none

      integer n,m,md
      real a(md,n),b(md,n),tm1,tm2,tm3

      integer i,k

      tm1=4.*real(n-2)*real(n-3)
      tm2=4.*real(n-1)*real(n-2)
      do k=1,m
         b(k,n-3)=tm1*a(k,n-1)
         b(k,n-2)=tm2*a(k,n)
         b(k,n-1)=0.
         b(k,n)=0.
      end do
      do i=n-4,1,-1
         tm1=2.*real(i+1)/real(i+2)
         tm2=real(i)/real(i+2)
         tm3=4.*real(i)*real(i+1)
         do k=1,m
            b(k,i)=tm1*b(k,i+2)-tm2*b(k,i+4)+tm3*a(k,i+2)
         end do
      end do
      do k=1,m
         b(k,1)=.5*b(k,1)
      end do

      end subroutine d2cheb
