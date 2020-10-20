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
c     Get an xz box from ur in core
c
c     boxr,boxi: returned box (nxp x mby x nzd)
c     yb:        y value of box
c     i:         component
c     ipad:      dealiasing if ipad = 1 (otherwise only oddball modes)
c                (upper part is padded...)
c                if ipad = 1, the data in z has the gap in the middle
c     ur, ur:    source of data
c
c     Pad for dealiasing if ipad = 1
c
      implicit none

      include 'par.f'

      integer yb,i,ipad
      real boxr((nxp/2+1),mby,nzd),boxi((nxp/2+1),mby,nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer x,y,z,xy,nzpad

      if (nproc.gt.1) then
         call stopnow(865785)
      end if

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
c
c     If non-symmetric, get also second half of z's and
c     put it in box according to ipad
c
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
c     Padding
c
      if (ipad.eq.1) then
c
c     Pad upper part/centre with zeros (including oddball)
c
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
c     Oddball zeroing only
c
         if (mod(nz,2).eq.0) then
c
c     Only for even nz, take away centre row
c
            do xy=1,(nxp/2+1)*min(mby,nyp-yb+1)
               boxr(xy,1,nz/2+1)=0.0
               boxi(xy,1,nz/2+1)=0.0
            end do
         end if
c
c     Take away last x row (boxi should already be zero)
c
         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               boxr(nx/2+1,y,z)=0.0
            end do
         end do
      end if

      end subroutine getxz
