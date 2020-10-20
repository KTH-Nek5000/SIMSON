c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine blth(graf,gxax,xmin,xmax,ymin,ymax,z,jvar,sym,dstar,
     &xr,xl,mxr,npoint,pxy,pxy2,pxz,wxy,ur,prex,prez,prey)
c
c     Calculates boundary layer thicknesses and shape factors
c     for boundary layers
c
      implicit none

      include 'par.f'

      integer jvar,z,mxr,npoint(3)
      real graf(nx+nyp+nz,3),gxax(nx+nyp+nz,3)
      real pxy(nx+1,nyp),pxy2(nx+1,nyp),xr,xl,dstar
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15),prey(nyp*2+15)
      real sym
      real xmin,xmax,ymin,ymax
      real wxy(nx+1,nyp)

      integer y,x,xpr
c
c     Get u(x,y)
c
      do y=1,nyp
         call getxzp(pxz,y,jvar,ur,sym)
         call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
         call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
         do x=1,nx
            xpr=mod(10*nx+x+mxr-1,nx)+1
            pxy(x,y)=pxz(xpr,z)
         end do
         pxy(nx+1,y)=pxy(1,y)
      end do
c
c     Calculate integrands
c
      do y=2,nyp
         do x=1,nx+1
            pxy2(x,y)=(1.-pxy(x,y)/pxy(x,1))*pxy(x,y)/pxy(x,1)
            pxy(x,y)=1.-pxy(x,y)/pxy(x,1)
         end do
      end do
      do x=1,nx+1
         pxy(x,1)=0.0
         pxy2(x,1)=0.0
      end do
c
c     Integrate to find displacement and momentum loss thicknesses
c     integration constants
c
      do x=1,nx+1
         graf(x,1)=0.0
      end do
      call vchbf(pxy,wxy,nyp,nx+1,nx+1,1,prey)
      call icheb(wxy,pxy,graf,nyp,nx+1,nx+1)
      do y=1,nyp
         do x=1,nx+1
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx+1,nx+1,1,prey)
      call vchbf(pxy2,wxy,nyp,nx+1,nx+1,1,prey)
      call icheb(wxy,pxy2,graf,nyp,nx+1,nx+1)
      do y=1,nyp
         do x=1,nx+1
            pxy2(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy2,wxy,nyp,nx+1,nx+1,1,prey)
      xmin=0.
      xmax=0.
      ymin=0.
      ymax=0.
      do x=1,nx+1
         gxax(x,1)=real(x-nx/2-1)/real(nx)*xl+xr
         gxax(x,2)=real(x-nx/2-1)/real(nx)*xl+xr
         gxax(x,3)=real(x-nx/2-1)/real(nx)*xl+xr
         graf(x,1)=(pxy(x,1)-pxy(x,nyp))/dstar*(2./real(nyp-1))
         graf(x,2)=(pxy2(x,1)-pxy2(x,nyp))/dstar*(2./real(nyp-1))
         graf(x,3)=graf(x,1)/graf(x,2)
      end do
      npoint(1)=nx+1
      npoint(2)=nx+1
      npoint(3)=nx+1

      return

      end subroutine blth
