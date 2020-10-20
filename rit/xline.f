c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine xline(graf,gxax,xmin,xmax,ymin,ymax,y,z,
     &    jvar,sym,mx,mxs,xl,mxr,xr,prex,prez,pxz,w,ur)
c
c     Cuts a line in the x-direction, possibly with subsampling
c
      implicit none

      include 'par.f'

      integer jvar,y,z,mx,mxs,mxr
      real graf(mx),gxax(mx)
      real xmin,xmax,ymin,ymax
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15)
      real sym
      real xr,xl
c
      integer x,xp,xpr

      open(unit=38,file='lineax.dat')
c
c     Cut the line
c
      call getxzp(pxz,y,jvar,ur,sym)
      call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
      call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
      do xp=1,nx,mxs
         xpr=mod(10*nx+xp+mxr-1,nx)+1
         x=(xp-1)/mxs+1
         graf(x)=pxz(xpr,z)
         gxax(x)=real(xp-nx/2-1)/real(nx)*xl+xr
         write(38,2200)gxax(x),graf(x)
      end do
c
c     Extend periodically in the x-direction
c
      graf(nx/mxs+1)=graf(1)
      gxax(nx/mxs+1)=real(nx+mxs-nx/2-1)/real(nx)*xl+xr
      write(38,2200)gxax(nx+1),graf(nx+1)
c
c     Grid
c
      ymin=0.
      ymax=0.
      do x=1,nx/mxs+1
         ymin=min(ymin,graf(x))
         ymax=max(ymax,graf(x))
      end do
      ymin=ymin*1.2
      ymax=ymax*1.2
      xmin=gxax(1)
      xmax=gxax(nx/mxs+1)
 2200 format(f16.10,'  ',e18.12)
      close(38)

      return

      end subroutine xline
