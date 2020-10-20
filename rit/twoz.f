c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine twoz(corr,gxax,xmin,xmax,ymin,ymax,y,jvar,sym,zl,
     &     pxz,w,ur,prex,prez)
c
c     Calculates two point correlations in the z-direction for one xz-plane
c
      implicit none

      include 'par.f'

      integer jvar,y
      real corr(nz/2+1),gxax(nz/2+1),xmin,xmax,ymin,ymax,zl
      real pxz(nx+2,nz+2),w(nx+2,nz+2)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15)
      real sym
c
      integer z,x
      real c0
      real prezr(nz+15)
c
      call getxzp(pxz,y,jvar,ur,sym)
c
c     Reset wavenumber zero, ie subtract mean
c
      pxz(1,1)=0.0
c
c     Transform to physical space
c
      call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
      call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
c
c     Prepare real transform in the z-direction
c
      call vrffti(nz,prezr,0)
c
c     Real transform in the z-direction
c
      call vrfftf(pxz,pxz(1,2),w,w(1,2),nz,nx,2*nx+4,1,prezr)
c
c     Calculate the two point correllation in z-fourier space
c
      do z=1,nz+2,2
         do x=1,nx
c
c     Multiply each element by its conjugate
c
            pxz(x,z)=pxz(x,z)*pxz(x,z)+pxz(x,z+1)*pxz(x,z+1)
            pxz(x,z+1)=0.0
         end do
      end do
      call vrfftb(pxz,pxz(1,2),w,w(1,2),nz,nx,2*nx+4,1,prezr)
c
c     Integrate in z-direction and normalize
c
      ymin=0.0
      ymax=1.1
      xmin=0.0
      xmax=zl*.5
      do z=1,nz/2+1
         corr(z)=0.0
         do x=1,nx
            corr(z)=corr(z)+pxz(x,z)
         end do
         if (z.eq.1) c0=corr(1)
         if (c0.eq.0) then
            corr(z)=1.0
         else
            corr(z)=corr(z)/c0
         end if
         gxax(z)=zl*real(z-1)/real(nz)
         if (corr(z).lt.ymin) ymin=corr(z)
      end do
      ymin=ymin*1.1

      return

      end subroutine twoz
