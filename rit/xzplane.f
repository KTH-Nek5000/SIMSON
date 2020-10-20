c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine xzplane(plane,mx,mz,mxs,mzs,grid,y,jvar,sym,
     &xr,xl,zl,mxr,w,ur,prex,prez)
c
c     Cuts a plane and calculates a grid in the xz-plane
c     possibly with subsampling
c
      implicit none

      include 'par.f'

      integer mx,mz,mxs,mzs,jvar,y,mxr
      real plane(mx,mz),grid(mx,mz,2),xr,xl,zl
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15)
      real sym
      integer z,zp,x,xp,xpr
      integer,save :: icount = 0
      character(len=20) fname

      call getxzp(pxz,y,jvar,ur,sym)
      call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
      call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
      do zp=1,nz,mzs
         z=(zp-1)/mzs+1
         do xp=1,nx,mxs
            xpr=mod(10*nx+xp+mxr-1,nx)+1
            x=(xp-1)/mxs+1
            plane(x,z)=pxz(xpr,zp)
         end do
c
c     Extend periodically in the x-direction
c
         plane(nx/mxs+1,z)=plane(1,z)
      end do
      do x=1,mx+1
         plane(x,mz)=plane(x,1)
      end do
c
c     Make grid
c
      do zp=1,nz+mzs,mzs
         z=(zp-1)/mxs+1
         do xp=1,nx+mxs,mxs
            x=(xp-1)/mxs+1
            grid(x,z,1)=real(xp-nx/2-1)/real(nx)*xl+xr
            grid(x,z,2)=real(2*zp-nz-2)/real(2*nz)*zl
         end do
      end do
c
c     write out a plane, either binary or ASCII
c
      if (1.eq.0) then
         open(unit=29,file='pianoxz.dat')
         do zp=1,nz+1
            do xp=1,nx+1
               write(29,2200)grid(xp,zp,1),grid(xp,zp,2),plane(xp,zp)
            end do
         end do
 2200    format(f18.10,f18.10,'  ',e18.12)
         close(29)
      else
         write(fname,'(a,i1,a)') 'plane_xz_',icount,'.dat'
         open(unit=29,file=trim(fname),form='unformatted')
         write(29) nx,nz
         write(29) grid(1:nx,1,1)
         write(29) grid(1,1:nz,2)
         write(29) plane(1:nx,1:nz)
         close(29)
         write(*,*) 'wrote ',trim(fname)
         icount = icount + 1
      end if


      end subroutine xzplane
