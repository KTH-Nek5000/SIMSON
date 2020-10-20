c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getxzc(pxz,uxz,it,sym,npl,uxzs,nzn)
c
c     Get one xzplane from a sequence on file or in core
c
      implicit none

      include 'par.f'

      integer npl,nzn
      real uxzs(nx,nzn,npl)
      real uxz(nx,nzn),pxz(nx+2,nz)
      real sym
      integer it

      integer x,z

      do z=1,nzn
         do x=1,nx
            uxz(x,z)=uxzs(x,z,it)
         end do
      end do

      do z=1,nzn
         do x=1,nx
            pxz(x,z)=uxz(x,z)
         end do
      end do
c
c     Unfold plane if symmetric
c
      if (nfzsym.eq.1) then
         do z=nz/2+2,nz
            do x=1,nx
               pxz(x,z)=pxz(x,nz+2-z)*sym
            end do
         end do
      end if

      return

      end subroutine getxzc
