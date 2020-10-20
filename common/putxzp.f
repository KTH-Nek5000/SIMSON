c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine putxzp(pln,yp,i,ur)
c
c     Put an xy plane into ur
c
      implicit none

      include 'par.f'

      integer yp,i
      complex pln(nx/2+1,nz)
      complex ur(memnx,memny,memnz,7)

      integer x,z

      do z=1,nzc
         do x=1,nx/2
            ur(x,yp,z,i)=pln(x,z)
         end do
      end do

      return

      end subroutine putxzp
