c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine corrp(a,b,n,m,md)

      implicit none

      integer n,m,md
      real a(md,n),b(md,n)

      integer i,j

      if (mod(n-1,2).eq.0) then
         do j=1,m
            a(j,1)=a(j,1)-.5*b(j,n)
         end do
         do i=3,n,2
            do j=1,m
               a(j,i-1)=a(j,i-1)-
     &              (1.+real((n-2)**2-(i-2)**2)/real(4*n-4))*b(j,n-1)
               a(j,i)=a(j,i)-b(j,n)
            end do
         end do
      else
         do j=1,m
            a(j,1)=a(j,1)-.5*(1.+real((n-2)**2)/real(4*n-4))*b(j,n-1)
            a(j,2)=a(j,2)-b(j,n)
         end do
         do i=4,n,2
            do j=1,m
               a(j,i-1)=a(j,i-1)-
     &              (1.+real((n-2)**2-(i-2)**2)/real(4*n-4))*b(j,n-1)
               a(j,i)=a(j,i)-b(j,n)
            end do
         end do
      end if

      end subroutine corrp
