c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine twoxs(corr,gxax,xmin,xmax,ymin,ymax,xl,
     &     pxz,w,it,npl,nzn,sym,uxzs,prex,lcmean,umean)
c
c     Calculates two point correlations in the x-direction for one xz-plane
c     read from a time sequence
c
      implicit none

      include 'par.f'

      integer it,npl,nzn
      logical lcmean
      real corr(nx/2+1),gxax(nx/2+1),xmin,xmax,ymin,ymax,xl
      real pxz(nx+2,nz),w(nx+2,nz)
      real uxzs(nx,nzn,npl)
      real prex(nx+15)
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
c     Fourier transform in the x-direction
c
      call vrfftf(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
c
c     Calculate the two point correlation in x-fourier space
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
         corr(x)=corr(x)/c0
         gxax(x)=xl*real(x-1)/real(nx)
         if (corr(x).lt.ymin) ymin=corr(x)
      end do
      ymin=ymin*1.1

      return

      end subroutine twoxs
