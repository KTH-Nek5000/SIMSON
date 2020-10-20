c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getxzs(plane,uxz,grid,u,uc,ialfa,prex,
     & xs,xr,it,sym,npl,uxzs,xl,zl,nzn,ivar,cpl,iu,fltype,high)
c
c     Get one xyplane from a sequence on file or in core and calculate grid
c
      implicit none

      include 'par.f'

      integer npl,nzn
      integer high
      real uxzs(nx,nzn,npl)
      real plane(nx+1,nz+1),uxz(nx,nzn),grid(nx+1,nz+1,2)
      real u(nx+2,nz)
      complex uc(nx/2+1,nz)
      complex ialfa(nx/2+1)
      real prex(nx+15)
      real xs,xr,xl,zl,sym,cpl
      integer it,iu
      integer ivar,fltype
      real sum

      integer x,z
      real ulam

      ulam=0.
      if ((fltype.eq.1.or.fltype.eq.4).and.ivar.eq.1) ulam=1.-cpl**2
      if ((fltype.eq.2.or.fltype.eq.5).and.ivar.eq.1) ulam=cpl
      if (iu.eq.1) ulam=0.

      do x=1,nx
         do z=1,nzn
            uxz(x,z)=uxzs(x,z,it)
         end do
      end do

      if (high.eq.1) then
         write(*,*) 'Do high pass filtering'
         do x=1,nx
            sum=0.
            do z=1,nzn
               sum=sum+uxz(x,z)
            end do
            sum=sum/real(nzn)
            do z=1,nzn
               uxz(x,z)=uxz(x,z)-sum
            end do
         end do
      end if
      do z=1,nzn
         do x=1,nx
            u(x,z)=uxz(x,z)-ulam
         end do
      end do
c
c     Fourier interpolation shifting of data
c     Fourier transform in the x-direction
c
      call vrfftf(u,u(2,1),plane,plane(2,1),nx,nzn,2,nx+2,prex)
      do z=1,nzn
         do x=1,nx/2+1
            uc(x,z)=uc(x,z)*exp((xr-xs)*ialfa(x))/real(nx)
         end do
         uc(nx/2+1,z)=real(uc(nx/2+1,z))
      end do
c
c     Transform to physical space
c
      call vrfftb(u,u(2,1),plane,plane(2,1),nx,nzn,2,nx+2,prex)
c
c     Copy into plane
c
      do z=1,nzn
         do x=1,nx
            plane(x,z)=u(x,z)
         end do
         plane(nx+1,z)=plane(1,z)
      end do
c
c     Unfold plane if symmetric
c
      if (nfzsym.eq.1) then
         do z=nz/2+2,nz
            do x=1,nx+1
               plane(x,z)=plane(x,nz+2-z)*sym
            end do
         end do
      end if
      do x=1,nx+1
         plane(x,nz+1)=plane(x,1)
      end do
c
c     Make grid
c
      do z=1,nz+1
         do x=1,nx+1
            grid(x,z,1)=real(x-nx/2-1)/real(nx)*xl+xr
            grid(x,z,2)=real(z-nz/2-1)/real(nz)*zl
         end do
      end do

      return

      end subroutine getxzs
