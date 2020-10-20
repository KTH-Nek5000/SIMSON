c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine pgmr(u1,xy,mxs,mxe,mys,mye,nx,ny,mbox,gl,gu)
c
c     Prints a pgm raw data file from the data in u1
c     the grid is taken from xy (should be cartesian)
c     plotting is done from mxs to mxe and from mys to mye
c
      integer mxs,mxe,mys,mye,nx,ny,mbox
      real u1(nx,ny),xy(nx,ny,2)
      integer nxm,nym
      parameter (nxm=640,nym=400)
      character*1 pixm(nxm*nym+100)
c
      integer i,j
      integer ix,ixl,ixu,iy,iyl,iyu,nxpix,nypix,ixd,iyd
      real umin,umax,xmin,xmax,ymin,ymax,p,q,gl,gu,ratio
      if (gl.eq.gu) then
c
c     Find range of data
c
         umin=1.E30
         umax=-1.E30
         do i=mys,mye
            do j=mxs,mxe
               umin=min(umin,u1(j,i))
               umax=max(umax,u1(j,i))
            end do
         end do
 2000    continue
         gl=umin
         gu=umax
      end if
c
c     Now find the gridsize
c
      xmin=min(xy(mxs,mys,1),xy(mxe,mys,1))
      xmax=max(xy(mxs,mys,1),xy(mxe,mys,1))
      ymin=min(xy(mxs,mys,2),xy(mxs,mye,2))
      ymax=max(xy(mxs,mys,2),xy(mxs,mye,2))
c
c     The direction of grid
c
      ixd=-1
      if (xy(mxs,mys,1).lt.xy(mxe,mys,1)) ixd=1
      iyd=-1
      if (xy(mxs,mys,2).lt.xy(mxs,mye,2)) iyd=1
      if (mbox.gt.1) ratio=real(mbox)/1000./.6
      if (mbox.lt.0) ratio=1000./real(-mbox)
      if (mbox.eq.1) ratio=real(nxm)/real(nym)
      if (mbox.eq.0) ratio=(xmax-xmin)/(ymax-ymin)
      if (ratio.gt.real(nxm)/real(nym)) then
         nxpix=nxm
         nypix=int(real(nxm)/ratio+.5)
       else
         nxpix=int(real(nym)*ratio+.5)
         nypix=nym
      end if
c
c     Create the pixmap through bilinear interpolation
c
      do i=mys,mye-1
         do j=mxs,mxe-1
            ixl=(xy(j,i,1)-xmin)*(real(nxpix-1)/(xmax-xmin))+1.5
            ixu=(xy(j+1,i,1)-xmin)*(real(nxpix-1)/(xmax-xmin))+1.5
            if (j.lt.mxe-1) ixu=ixu-ixd
            iyl=(xy(j,i,2)-ymin)*(real(nypix-1)/(ymax-ymin))+1.5
            iyu=(xy(j,i+1,2)-ymin)*(real(nypix-1)/(ymax-ymin))+1.5
            if (i.lt.mye-1) iyu=iyu-iyd
            do iy=iyl,iyu,iyd
               do ix=ixl,ixu,ixd
                  p=real(iy-iyl)/real(iyu-iyl+iyd)
                  q=real(ix-ixl)/real(ixu-ixl+ixd)
                  pixm(ix+(nypix-iy)*nxpix+100)=char(min(255,max(0,int(
     &                 ((1.-p)*((1.-q)*u1(j,i)+q*u1(j+1,i))+
     &                 p*((1.-q)*u1(j,i+1)+q*u1(j+1,i+1))-gl)*
     &                 255./(gu-gl)+.5))))
               end do
            end do
         end do
      end do
      call wpgmr(pixm,nxpix,nypix,255)

      return

      end subroutine pgmr
