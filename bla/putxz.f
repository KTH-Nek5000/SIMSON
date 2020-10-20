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
c     size of box:     nxp/2+1 * mby * nzd
c     copied elements: nx/2    * mby * (2*nz/2)
c
c     Note that the whole box is copied, including the oddball modes!
c
c     Note that ur is on the coarse grid in x and z, i.e. nx/2*nyp*nz
c
      implicit none

      include 'par.f'

      integer yb,i
      real boxr((nxp/2+1),mby,nzd),boxi((nxp/2+1),mby,nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer x,z,y

      if (nproc.gt.1) then
         call stopnow(865786)
      end if

      do z=1,nz/2
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               ur(x,y,z,i)=boxr(x,y-yb+1,z)
               ui(x,y,z,i)=boxi(x,y-yb+1,z)
            end do
         end do
      end do
      if (nfzsym.eq.0) then
c
c     If non-symmetric, put also second half of z's
c
         do z=nz/2+1,nz
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  ur(x,y,z,i)=boxr(x,y-yb+1,z+nzp-nz)
                  ui(x,y,z,i)=boxi(x,y-yb+1,z+nzp-nz)
               end do
            end do
         end do
      end if

      end subroutine putxz
