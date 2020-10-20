c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine icheb(b,a,intcon,n,m,md)
c
c     m integrations of n coefficients in Chebychev space
c
c     The result is returned in b. The integration
c     constant intcon must be given and is returned as the zeroth
c     order Chebyshev coefficient.
c
      implicit none

      integer n,m,md
      real a(md,n),intcon(md),b(md,n),dnm

      integer k,i,nm1

      nm1=n-1
      dnm=.5/real(nm1)
      do i=3,nm1
         do k=1,m
            b(k,i)=(a(k,i-1)-a(k,i+1))*(.5/real(i-1))
         end do
      end do
      do k=1,m
         b(k,n)=a(k,nm1)*dnm
c     c0=2
         b(k,2)=a(k,1)-a(k,3)*.5
         b(k,1)=intcon(k)
      end do

      end subroutine icheb
