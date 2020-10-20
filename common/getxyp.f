c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getxyp(pln,z,i,ur)
c
c     Get an xy box from ur
c
      implicit none

      include 'par.f'

      integer z,i
      complex pln(nx/2,nyp)
      complex ur(memnx,memny,memnz,7)

      integer x,y
      do y=1,nyp
         do x=1,nx/2
            pln(x,y)=ur(x,y,z,i)
         end do
      end do

      return

      end subroutine getxyp
