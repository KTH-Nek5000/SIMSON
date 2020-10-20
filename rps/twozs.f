c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine twozs(corr,gxax,xmin,xmax,ymin,ymax,sym,zl,
     &     pxz,w,it,npl,nzn,uxzs,prezr,lcmean,umean)
c
c     Calculates two point correlations in the x-direction for one xz-plane
c
      implicit none

      include 'par.f'

      integer it,npl,nzn
      logical lcmean
      real corr(nz/2+1),gxax(nz/2+1),xmin,xmax,ymin,ymax,zl
      real pxz(nx+2,nz+2),w(nx+2,nz+2)
      real uxzs(nx,nzn,npl)
      real prezr(nz+15)
      real sym,umean

      integer z,x
      real c0

      call getxzc(pxz,w,it,sym,npl,uxzs,nzn)
c
c     Subtract mean
c
      if (.not.lcmean) then
         umean=0.0
         do z=1,nz
            do x=1,nx
               umean=umean+pxz(x,z)
            end do
         end do
         umean=umean/real(nx*nz)
      end if
      do z=1,nz
         do x=1,nx
            pxz(x,z)=pxz(x,z)-umean
         end do
      end do
c
c     Real transform in the z-direction
c
      call vrfftf(pxz,pxz(1,2),w,w(1,2),nz,nx,2*nx+4,1,prezr)
c
c     Calculate the two point correlation in z-fourier space
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
         corr(z)=corr(z)/c0
         gxax(z)=zl*real(z-1)/real(nz)
         if (corr(z).lt.ymin) ymin=corr(z)
      end do
      ymin=ymin*1.1

      return

      end subroutine twozs
