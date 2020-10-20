c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine twox(corr,gxax,xmin,xmax,ymin,ymax,y,jvar,sym,xl,
     &     pxz,w,ur,prex,prez)
c
c     Calculates two point correlations in the x-direction for one xz-plane
c
      implicit none

      include 'par.f'

      integer jvar,y
      real corr(nx/2+1),gxax(nx/2+1),xmin,xmax,ymin,ymax,xl
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15)
      real sym
c
      integer z,x
      real c0
c
      call getxzp(pxz,y,jvar,ur,sym)
c
c     Reset wavenumber zero, ie subtract mean
c
      pxz(1,1)=0.0
      call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
c
c     Calculate the two point correllation in x-fourier space
c
      do z=1,nz
         do x=1,nx+2,2
c
c     Multiply each element by its conjugate
c
            pxz(x,z)=pxz(x,z)*pxz(x,z)+pxz(x+1,z)*pxz(x+1,z)
            pxz(x+1,z)=0.0
         end do
      end do
c
c     Transform to physical space
c
      call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
c
c     Integrate in z-direction and normalize
c
      ymin=0.0
      ymax=1.1
      xmin=0.0
      xmax=xl*.5
      do x=1,nx/2+1
         corr(x)=0.0
         do z=1,nz
            corr(x)=corr(x)+pxz(x,z)
         end do
         if (x.eq.1) c0=corr(1)
         if (c0.eq.0) then
            corr(x)=1.0
         else
            corr(x)=corr(x)/c0
         end if
         gxax(x)=xl*real(x-1)/real(nx)
         if (corr(x).lt.ymin) ymin=corr(x)
      end do
      ymin=ymin*1.1

      return

      end subroutine twox
