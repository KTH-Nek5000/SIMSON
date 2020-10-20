c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine yzplane(plane,my,mz,mys,mzs,grid,x,jvar,sym,
     &zl,pxz,w,ur,prex,prez,eta)
c
c     Cuts a plane and calculates a grid in the yz-plane
c     possibly with subsampling
c
      implicit none

      include 'par.f'

      integer mz,my,mzs,mys,jvar,x
      real plane(mz,my),grid(mz,my,2),zl
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15),eta(nyp)
      real sym
c
      integer y,yp,z,zp,x1
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
         do zp=1,nz,mzs
            z=(zp-1)/mzs+1
            x1=mod(10*nx+x-1,nx)+1
            plane(z,y)=pxz(x1,zp)
         end do
c
c     Extend periodically in the z-direction
c
         plane(nz/mzs+1,y)=plane(1,y)
      end do
c
c     Make grid
c
      do yp=1,nyp,mys
         y=(yp-1)/mys+1
         do zp=1,nz+mzs,mzs
            z=(zp-1)/mzs+1
            grid(z,my+1-y,1)=real(2*zp-nz-2)/real(2*nz)*zl
            grid(z,my+1-y,2)=eta(yp)
         end do
      end do

      open(unit=29,file='pianoyz.dat')
      do yp=1,nyp
         do zp=1,nz+1
            write(29,2200)grid(zp,yp,1),grid(zp,yp,2),plane(zp,yp)
         end do
      end do

 2200 format(f15.10,f15.10,'  ',e18.12)
      close(29)

      return

      end subroutine yzplane
