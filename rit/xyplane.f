c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine xyplane(plane,mx,my,mxs,mys,grid,z,jvar,sym,
     &xr,xl,mxr,pxz,w,ur,prex,prez,eta)
c
c     Cuts a plane and calculates a grid in the xy-plane
c     possibly with subsampling
c
      implicit none

      include 'par.f'

      integer mx,my,mxs,mys,jvar,z,mxr
      real plane(mx,my),grid(mx,my,2),xr,xl
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15),eta(nyp)
      real sym
c
      integer y,yp,x,xp,xpr

      integer luca1,luca2,luca3
c
c     Cut the plane
c
      do yp=1,nyp,mys
        y=(yp-1)/mys+1
c
c     Note how the plane is turned upside down to compensate for reverse numbering
c
        call getxzp(pxz,nyp+1-yp,jvar,ur,sym)
        call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
        call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
        do xp=1,nx,mxs
          xpr=mod(10*nx+xp+mxr-1,nx)+1
          x=(xp-1)/mxs+1
          plane(x,y)=pxz(xpr,z)
       end do
c
c     Extend periodically in the yx-direction
c
        plane(nx/mxs+1,y)=plane(1,y)
      end do
c
c     Make grid
c
      do yp=1,nyp,mys
         y=(yp-1)/mys+1
         do xp=1,nx+mxs,mxs
            x=(xp-1)/mxs+1
            grid(x,my+1-y,1)=real(xp-nx/2-1)/real(nx)*xl+xr
            grid(x,my+1-y,2)=eta(yp)
            if (grid(x,my+1-y,1).lt.199.0) luca1=x
            if (grid(x,my+1-y,1).lt.250.1) luca2=x
            if (grid(x,my+1-y,2).gt.8.0) luca3=my+1-y
cc          write(*,*)real(xp-nx/2-1)/real(nx)*xl+xr,eta(yp)
         end do
      end do
cc      write(*,*)mxs,mys
      open(unit=29,file='pianoxy.dat')
      do yp=1,nyp
         do xp=1,nx+1
            write(29,2200)grid(xp,yp,1),grid(xp,yp,2),plane(xp,yp)
         end do
      end do
 2200 format(f15.10,f15.10,'  ',e18.12)
      close(29)
c      do x=luca1,luca2
c         do y=1,luca3
c            write(21,*)plane(x,y)
c         end do
c         write(23,*)grid(x,100,1)
c      end do
c      do y=1,luca3
c         write(22,*)eta(nyp+1-y)
c      end do
c
c      write(*,*)luca1,luca2,luca3,eta(nyp+1-luca3),my

      return

      end subroutine xyplane
