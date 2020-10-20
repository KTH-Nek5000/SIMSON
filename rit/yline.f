c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine yline(graf,gxax,xmin,xmax,ymin,ymax,x,z,
     &    jvar,sym,my,mys,eta,prex,prez,pxz,w,ur)
c
c     Cuts a line in the y-direction, possibly with subsampling
c
      implicit none

      include 'par.f'

      integer jvar,x,z,my,mys
      real graf(my),gxax(my)
      real xmin,xmax,ymin,ymax
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15),eta(nyp)
      real sym
c
      integer yp,y
c
c     Cut the line
c
      open(unit=38,file='lineay.dat')
      x=mod(10*nx+x-1,nx)+1
      do yp=1,nyp,mys
         y=(yp-1)/mys+1
c
c     Note how the line is turned upside down to account for reverse numbering
c
         call getxzp(pxz,nyp+1-yp,jvar,ur,sym)
         call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
         call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
         graf(y)=pxz(x,z)
         gxax(y)=eta(nyp+1-yp)
         write(38,2200)gxax(y),graf(y)
      end do
      ymin=0.
      ymax=0.
      do y=1,my
         ymin=min(ymin,graf(y))
         ymax=max(ymax,graf(y))
      end do
      ymin=ymin*1.2
      ymax=ymax*1.2
      xmin=eta(nyp)
      xmax=eta(1)
 2200 format(f15.10,'  ',e18.12)
      close(38)

      return

      end subroutine yline
