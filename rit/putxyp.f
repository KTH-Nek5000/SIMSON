c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine putxyp(pln,z,i,ur)
c
c     Put an xy box to ur
c
      implicit none

      include 'par.f'

      integer z,i
      complex pln(nx/2,nyp)
      complex ur(memnx,memny,memnz,7)
c
      integer x,y
c
      do y=1,nyp
         do x=1,nx/2
            ur(x,y,z,i)=pln(x,y)
         end do
      end do

      return

      end subroutine putxyp
