c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getxys(plane,uxy,grid,u,uc,ialfa,prex,
     &xs,xr,it,npl,uxys,eta,xl,fltype,dstar,ivar,idudy,prey)
c
c     Get one xyplane from a sequence on file or in core and calculate grid
c
      implicit none

      include 'par.f'

      integer npl,ivar,idudy
      real plane(nx+1,nyp),uxy(nx,nyp),uxys(nx,nyp,npl),grid(nx+1,nyp,2)
      real prey(nyp+15)
      real u(nx+2,nyp)
      complex uc(nx/2+1,nyp)
      complex ialfa(nx/2+1)
      real prex(nx+15)
      real eta(nyp)
      real xs,xr,xl,dstar
      integer it,fltype

      integer x,y
      real u0,uwall

      do y=1,nyp
         do x=1,nx
            uxy(x,y)=uxys(x,y,it)
         end do
      end do
      if (ivar.eq.1.and.idudy.le.1.or.idudy.ge.4) then
         do y=1,nyp
            if (idudy.le.1) then
               u0=0.
               do x=1,nx
                  u0=u0+uxy(x,y)
               end do
               u0=u0/real(nx)
            else
               if (fltype.eq.1.or.fltype.eq.4)u0=1.-eta(y)**2+uxy(1,nyp)
               if (fltype.eq.2.or.fltype.eq.5) u0=eta(y)
            end if
            do x=1,nx
               uxy(x,y)=uxy(x,y)-u0
            end do
         end do
      end if
      if (fltype.ne.2.and.fltype.ne.5.and.ivar.eq.1.and.idudy.eq.2) then
         do y=1,nyp
            uwall=uxy(1,nyp)
            do x=1,nx
               uxy(x,y)=uxy(x,y)-uwall
            end do
         end do
      end if
c
      do y=1,nyp
         do x=1,nx
            u(x,nyp+1-y)=uxy(x,y)
         end do
      end do
c
c     Fourier interpolation shifting of data
c     Fourier transform in the x-direction
c
      call vrfftf(u,u(2,1),plane,plane(2,1),nx,nyp,2,nx+2,prex)
      do y=1,nyp
         do x=1,nx/2+1
            uc(x,y)=uc(x,y)*exp((xr-xs)*ialfa(x))/real(nx)
         end do
         uc(nx/2+1,y)=real(uc(nx/2+1,y))
      end do
c
c     Transform to physical space
c
      call vrfftb(u,u(2,1),plane,plane(2,1),nx,nyp,2,nx+2,prex)
c
c     Copy into plane
c
      do y=1,nyp
         do x=1,nx
            plane(x,y)=u(x,y)
         end do
         plane(nx+1,y)=plane(1,y)
      end do
c
c     Differentiate to make dudy
c
      if (idudy.eq.1.or.idudy.eq.3.or.idudy.eq.5) then
c
c     Obs ! using grid as work storage
c
        call vchbf(plane,grid,nyp,nx+1,nx+1,1,prey)
        call rdcheb(plane,nyp,nx+1,nx+1)
        call vchbb(plane,grid,nyp,nx+1,nx+1,1,prey)
        if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6) then
           do y=1,nyp
              do x=1,nx+1
c
c     Minus sign corrects for derivative of upside down plane
c
                 plane(x,y)=-plane(x,y)*dstar*(2./real(nyp-1))
              end do
           end do
        else
           do y=1,nyp
              do x=1,nx+1
c
c     Minus sign corrects for derivative of upside down plane
c
                 plane(x,y)=-plane(x,y)*(2./real(nyp-1))
              end do
           end do
        end if
c
c     Normalize
c
      end if
c
c     Make grid
c
c      if (fltype.ne.3) then
      do y=1,nyp
         do x=1,nx+1
            grid(x,nyp+1-y,1)=real(x-nx/2-1)/real(nx)*xl+xr
            grid(x,nyp+1-y,2)=eta(y)
         end do
      end do
c      else
c        do 2040 y=1,nyp
c        do 2040 x=1,nx+1
c          grid(x,nyp+1-y,1)=real(x-nx/2-1)/real(nx)*xl+xr
c          grid(x,nyp+1-y,2)=(eta(y)+1.)/dstar
c 2040   continue
c      end if

      return

      end subroutine getxys
