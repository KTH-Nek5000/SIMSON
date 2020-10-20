c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************

      subroutine bla2arpck(workd,ur,ui,maxn)
c
c     Transform bla-format to workd(ipntr(1))
c

      implicit none

      include 'par.f'

      integer maxn
      real ur(memnx,memny,memnz,3),ui(memnx,memny,memnz,3)
      real workd(maxn)
      integer i,ll,x,y,z,jj

      ll = 3
      do i=1,ll
         do z=1,nz/nproc
            do y = 1,ny
               do x=1,nx/2
                  jj = (i-1)*nz*ny*nx/nproc +
     &                 (z-1)*ny*nx + (y-1)*nx

                  workd(jj+2*x-1) = ur(x,y,z,i)
                  workd(jj+2*x)  =  ui(x,y,z,i)
               end do
            end do
         end do
      end do

      end subroutine bla2arpck
