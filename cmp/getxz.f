c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getxz(boxr,boxi,yb,i,ipad,ur,ui)
c
c     Get an xz box from ur
c     Pad for dealiasing if ipad = 1
c
      implicit none

      include 'par.f'

      integer yb,i,ipad
      real boxr((nxp/2+1),mby,nzd),boxi((nxp/2+1),mby,nzd)
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
c
      integer x,y,z,xy,nzpad
      nzpad=(nzp-nz)*ipad
c
c     Get box from core
c
      do z=1,nz/2
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               boxr(x,y-yb+1,z)=ur(x,y,z,i)
               boxi(x,y-yb+1,z)=ui(x,y,z,i)
            end do
         end do
      end do
      if (nfzsym.eq.0) then
         do z=nz/2+1,nz
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  boxr(x,y-yb+1,z+nzpad)=ur(x,y,z,i)
                  boxi(x,y-yb+1,z+nzpad)=ui(x,y,z,i)
               end do
            end do
         end do
      end if
c
c     Pad with zeros
c
      if (ipad.eq.1) then
         do z=(nz+1)/2+1,min(nzpc,nzp+1-nz/2)
            do xy=1,(nxp/2+1)*min(mby,nyp-yb+1)
               boxr(xy,1,z)=0.0
               boxi(xy,1,z)=0.0
            end do
         end do
         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               do x=nx/2+1,nxp/2+1
                  boxr(x,y,z)=0.0
                  boxi(x,y,z)=0.0
               end do
            end do
         end do
      else
c
c     Odd ball zeroing only
c
         if (mod(nz,2).eq.0) then
            do xy=1,(nxp/2+1)*min(mby,nyp-yb+1)
               boxr(xy,1,nz/2+1)=0.0
               boxi(xy,1,nz/2+1)=0.0
            end do
         end if
         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               boxr(nx/2+1,y,z)=0.0
            end do
         end do
      end if

      return

      end subroutine getxz
