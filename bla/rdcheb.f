c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdcheb(data,n,m,md)
c
c     Calculates first derivative in Chebyshev space
c
c     n number of Chebyshev coefficients
c     m number of derivatives
c     md first index length on data array
c
      implicit none

      integer n,m,md
      real data(md,n)

      integer i,k
      real tmp1
      do k=1,m
         tmp1=data(k,n-2)
         data(k,n-2)=2.*real(n-2)*data(k,n-1)
         data(k,n-1)=2.*real(n-1)*data(k,n)
         data(k,n)=tmp1
      end do
      do i=n-3,1,-1
         do k=1,m
            tmp1=data(k,n)
            data(k,n)=data(k,i)
            data(k,i)=data(k,i+2)+2.*real(i)*tmp1
         end do
      end do
      do k=1,m
         data(k,n)=0.0
         data(k,1)=.5*data(k,1)
      end do

      return

      end subroutine rdcheb
