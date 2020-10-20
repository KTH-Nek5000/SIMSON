c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine putxz(boxr,boxi,yb,i,ur,ui)
c
c     Put an xy box into ur
c
      implicit none

      include 'par.f'

      integer yb,i
      real boxr((nxp/2+1),mby,nzd),boxi((nxp/2+1),mby,nzd)
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
c
      integer x,z,y
c
c     Put box in core
c
        do z=1,nz/2
           do y=yb,min(nyp,yb+mby-1)
              do x=1,nx/2
                 ur(x,y,z,i)=boxr(x,y-yb+1,z)
                 ui(x,y,z,i)=boxi(x,y-yb+1,z)
              end do
           end do
        end do
        if (nfzsym.eq.0) then
           do z=nz/2+1,nz
              do y=yb,min(nyp,yb+mby-1)
                 do x=1,nx/2
                    ur(x,y,z,i)=boxr(x,y-yb+1,z+nzp-nz)
                    ui(x,y,z,i)=boxi(x,y-yb+1,z+nzp-nz)
                 end do
              end do
           end do
        end if

      return

      end subroutine putxz
