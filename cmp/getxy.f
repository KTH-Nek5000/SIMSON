c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************

      subroutine getxy(boxr,boxi,zb,i,ur,ui)
c
c     Get an xy box from ur
c
      implicit none

      include 'par.f'

      integer zb,i
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
c
      integer x,y,z

      do z=zb,zb+mbz-1
         do y=1,nyp
            do x=1,nx/2
               boxr(x,z-zb+1,y)=ur(x,y,z,i)
               boxi(x,z-zb+1,y)=ui(x,y,z,i)
            end do
         end do
      end do

      return

      end subroutine getxy
