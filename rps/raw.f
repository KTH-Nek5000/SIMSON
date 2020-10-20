c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine raw(u1,xy,mxs,mxe,mys,mye,nx,ny,mbox,gl,gu,wshow,ifil,
     &     uumin,uumax,nothing)
c
c     Prints a raw data file from the data in u1
c     the grid is taken from xy (should be cartesian)
c     plotting is done from mxs to mxe and from mys to mye
c
      implicit none

      integer mxs,mxe,mys,mye,nx,ny,mbox
      real u1(nx,ny),xy(nx,ny,2),gl,gu,ratio
      integer nxm,nym
      logical wshow

      parameter (nxm=480,nym=304)
      character*1 pixm(nxm*nym)
      logical lagl,nothing

      integer i,j
      integer ix,ixl,ixu,iy,iyl,iyu,nxpix,nypix,ixd,iyd,ifil
      real umin,umax,xmin,xmax,ymin,ymax,p,q
      real uumin,uumax
c
c     Find range of data
c
      lagl=gl.eq.gu
        umin=1.E30
        umax=-1.E30
        do i=mys,mye
           do j=mxs,mxe
              umin=min(umin,u1(j,i))
              umax=max(umax,u1(j,i))
           end do
        end do
c
c        write(*,*) mxs,mxe,mys,mye
c        write(*,*) nx,ny
c        write(*,*) mbox
c        write(*,*) xy(1,1,1),xy(1,1,2)
c        write(*,*) xy(nx,ny,1),xy(nx,ny,2)

c        do i=1,nx
c           do j=1,ny
c              ix = xy(i,j,1)
c              xy(i,j,1) = xy(i,j,2)
c              xy(i,j,2) = ix
c           end do
c        end do
      if (lagl) then
         gl=umin
         gu=umax
      end if
      uumin=min(uumin,umin)
      uumax=max(uumax,umax)
      if (nothing) then
         return
      end if
c
c     Now find the gridsize
c
      xmin=min(xy(mxs,mys,1),xy(mxe,mye,1))
      xmax=max(xy(mxs,mys,1),xy(mxe,mye,1))
      ymin=min(xy(mxs,mys,2),xy(mxe,mye,2))
      ymax=max(xy(mxs,mys,2),xy(mxe,mye,2))
c      write(*,*) 'mm ',umin,umax,xmax,xmin,ymax,ymin
c      write(*,*) mxs,mys,mxe,mye
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
c     write(*,*) i,j,iyl,iyu,iyd,ixl,ixu,ixd
            do iy=iyl,iyu,iyd
               do ix=ixl,ixu,ixd
                  p=real(iy-iyl)/real(iyu-iyl+iyd)
                  q=real(ix-ixl)/real(ixu-ixl+ixd)
                  pixm(ix+(nypix-iy)*nxpix)=char(min(255,max(0,int(
     &                 ((1.-p)*((1.-q)*u1(j,i)+q*u1(j+1,i))+
     &                 p*((1.-q)*u1(j,i+1)+q*u1(j+1,i+1))-gl)*
     &                 255./(gu-gl)+.5))))
c        write(*,*) p,q,ichar(pixm(ix+(nypix-iy)*nxpix))
c          write(*,*) ix,ixl,ixu
c         write(*,*) iy,iyl,iyu
               end do
            end do
         end do
      end do
c      call wraw(pixm,nxpix,nypix,ifil)
      call wpgmr(pixm,nxpix,nypix,ifil)
c
c     Write a show file
c
      wshow = .false.
      if (wshow) then
         open(unit=46,file='show')
         write(46,9000) nypix,nxpix
 9000    format('xmovie -height ',i4,' -width ',i4,' plot*.raw')
         close(unit=46)
      end if
      if (lagl) gl=gu

      return

      end subroutine raw
