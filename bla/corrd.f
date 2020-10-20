c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine corrd(a,n,m,md)

      implicit none

      integer n,m,md
      real a(md,n)
      integer i,j

      if (mod(n-1,2).eq.0) then
         do j=1,m
            a(j,1)=a(j,1)-.5*a(j,n)
         end do
         do i=3,n,2
            do j=1,m
               a(j,i)=a(j,i)-a(j,n)
            end do
         end do
      else
         do i=2,n,2
            do j=1,m
               a(j,i)=a(j,i)-a(j,n)
            end do
         end do
      end if

      end subroutine corrd
