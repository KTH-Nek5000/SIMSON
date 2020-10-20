c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine boxcfl(cflp,u2r,u2i,yb,deta,xl,zl,wbci)
c
c     Finds the maximum cflp number in a u2 box (containing u,v,y)
c     cflp is a partial cfl number : abs(u)/dx+...
c
      implicit none

      include 'par.f'

      integer yb,wbci
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real cflp,deta(nyp),xl,zl

      integer nxy,x,y,z,xy
      real cflv((nxp/2+1)*mby),dyi((nxp/2+1)*mby),dxi,dzi,cfll,dyinv

      nxy=(nxp/2+1)*min(mby,nyp-yb+1)
c
c     Compute inverse of grid spacing in y
c
      do y = yb,min(yb+mby-1,nyp)
         dyinv = real(ny-1)/real(nyp-1)/deta(y)
c
c     For y close to the surfaces the estimate of the eigenvalue
c     from pi/dy is too high and is replaced by the below
c     which is the (conservative) estimate of lambda< N^2/16 for N>=16
c     this allows roughly 20 times larger timestep than the regular estimate
c     if the timestep is set by the vertical flow through the boundary
c
         if (scalar.eq.0 .and. wbci.eq.0) then
            if (ny.ge.17) dyinv=min(dyinv,real(ny-1)**2/10./3.14)
         end if
         if (scalar.ne.0 .and. wbci.eq.0) then
            if (ny.ge.17) dyinv=min(dyinv,real(ny-1)**2/10./3.14*4)
         end if
c
c     The new conditions causes the increasingly small dy close
c     to the walls to be ignored and dy is set to the value
c     approx. 4 points off the boundary.
c
         do x=1,nxp/2+1
            xy=x+(nxp/2+1)*(y-yb)
            dyi(xy)=dyinv
            cflv(xy)=0.0
         end do
      end do
c
c     Inverse of grid spacing in x and z
c
      dxi = real(nx-2)/xl
      dzi = real(nz-2)/zl
      if (nz.eq.1) dzi=0.

      do z=1,nzpc
         do xy=1,nxy
            cfll=abs(u2r(xy,z,1))*dxi+abs(u2r(xy,z,2))*dyi(xy)+
     &           abs(u2r(xy,z,3))*dzi
            cflv(xy)=max(cflv(xy),cfll)
            cfll=abs(u2i(xy,z,1))*dxi+abs(u2i(xy,z,2))*dyi(xy)+
     &           abs(u2i(xy,z,3))*dzi
            cflv(xy)=max(cflv(xy),cfll)
         end do
      end do

      cflp=0.
      do xy=1,nxy
         cflp=max(cflp,cflv(xy))
      end do

      end subroutine boxcfl
