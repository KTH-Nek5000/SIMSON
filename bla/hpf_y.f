c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine hpf_y(data,gewy1,gewy2,nxh,nyp,mbz,diags,rs)
c
c     High-pass filter (I-G)
c     u" = (I-G)*u: A*u" = (A-B)*u
c     low-pass filter G
c     u' = G*u:     A*u' = B*u
c
c     Input and output is in data
c     A temporary storage in rs needs to be provided
c
      implicit none

      integer nxh,nyp,mbz
      real data(nxh*mbz,nyp)
      real gewy1(nyp,5),gewy2(nyp,5)

      real diags(nyp,5)
      real rs(nxh*mbz,nyp)

      integer i,j,l
      integer irep,ifail

      real fact3z,factm3z

      rs = 0.
c
c     Explicit filter
c
c     - Boundary points already zero
c     - Asymmetric stencil (k=2 and k=nyp-1)
c
      do l=1,5
         do i=1,nxh*mbz
            rs(i,2)     = rs(i,2    ) +
     &           (gewy2(2,l)-gewy1(2,l))    *data(i,l)
            rs(i,nyp-1) = rs(i,nyp-1) +
     &           (gewy2(nyp-1,l)-gewy1(nyp-1,l))*data(i,nyp-5+l)
         end do
      end do
c
c     - Symmetric stencil
c
      do l=1,5
         do j=3,nyp-2
            do i=1,nxh*mbz
               rs(i,j) = rs(i,j) +
     &              (gewy2(j,l)-gewy1(j,l))*data(i,j-3+l)
            end do
         end do
      end do
c
c     Implicit filter
c
      fact3z  = -gewy2(2,5)     / gewy2(3,5)
      factm3z = -gewy2(nyp-1,1) / gewy2(nyp-2,1)

      data = rs

      do i=1,nxh*mbz
         rs(i,2)     = rs(i,2)     + rs(i,3)    *fact3z
         rs(i,nyp-1) = rs(i,nyp-1) + rs(i,nyp-2)*factm3z
      end do
c
c     Compute solution (LU decomposition computed in init_filter)
c
      irep = 1
      call pntdiag(nyp,nxh*mbz,nxh*mbz,
     &     diags(1,1),diags(1,2),diags(1,3),diags(1,4),diags(1,5),
     &     rs,data,ifail,irep)
      if (ifail.ne.-1) then
         write(*,*) 'Error in pntdiag. Stop.'
         stop
      end if

      end subroutine hpf_y
