c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine zwall(ur,pxy)

      implicit none

      include 'par.f'

      complex ur(memnx,memny,memnz,3)
      complex pxy(nx/2,nyp)
c
      integer y

      call getxyp(pxy,1,1,ur)
      do y=1,nyp
         pxy(1,y)=pxy(1,y)-pxy(1,nyp)
      end do
      call putxyp(pxy,1,1,ur)

      return

      end subroutine zwall
