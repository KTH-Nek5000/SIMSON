c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getyzs(plane,uyz,grid,it,sym,npl,uyzs,eta,zl,
     &     nzn,fltype,dstar,ivar,idudy,prey)
c
c     Get one yz-plane from a sequence on file or in core and calculate grid
c
      implicit none

      include 'par.f'

      integer npl,ivar,idudy,nzn
      real plane(nyp,nz+1),uyz(nyp,nz+1),uyzs(nyp,nzn,npl)
      real grid(nyp,nz+1,2)
      real prey(nyp+15)
      real eta(nyp)
      real zl,dstar,sym
      integer it,fltype

      integer y,z
      real u0,uwall

      do y=1,nyp
         do z=1,nzn
            uyz(y,z)=uyzs(y,z,it)
         end do
      end do

      if (ivar.eq.1.and.idudy.le.1.or.idudy.ge.4) then
         do y=1,nyp
            u0=0.
            if (idudy.le.1) then
               do z=1,nzn
                  u0=u0+uyz(y,z)
               end do
               u0=u0/real(nzn)
            else
               if (fltype.eq.1.or.fltype.eq.4)
     &              u0=1.-eta(y)**2+uyz(nyp,1)
               if (fltype.eq.2.or.fltype.eq.5) u0=eta(y)
            end if
            do z=1,nzn
               uyz(y,z)=uyz(y,z)-u0
            end do
         end do
      end if

      if (fltype.ne.2.and.fltype.ne.5.and.ivar.eq.1.and.idudy.eq.2) then
         do y=1,nyp
            uwall=uyz(nyp,1)
            do z=1,nzn
               uyz(y,z)=uyz(y,z)-uwall
            end do
         end do
      end if

      do y=1,nyp
         do z=1,nzn
            plane(nyp+1-y,z)=uyz(y,z)
         end do
      end do
c
c     Differentiate to make dudy
c
      if (idudy.eq.1.or.idudy.eq.3.or.idudy.eq.5) then
c
c     Obs ! using grid as work storage
c
         call vchbf(plane,grid,nyp,nzn,1,nz+1,prey)
         call rdcheb(plane,nyp,nx+1,nx+1)
         call vchbb(plane,grid,nyp,nzn,1,nz+1,prey)
         if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6) then
            do y=1,nyp
               do z=1,nzn
c
c     Minus sign corrects for derivative of upside down plane
c     Normalize
c
                  plane(y,z)=-plane(y,z)*dstar*(2./real(nyp-1))
               end do
            end do
         else
            do y=1,nyp
               do z=1,nzn
c
c     Minus sign corrects for derivative of upside down plane
c     Normalize
c
                  plane(y,z)=-plane(y,z)*(2./real(nyp-1))
               end do
            end do
         end if
      end if
c
c     Unfold plane if symmetric
c
      if (nfzsym.eq.1) then
         do z=nz/2+2,nz
            do y=1,nyp
               plane(y,z)=plane(y,nz+2-z)*sym
            end do
         end do
      end if
      do y=1,nyp
         plane(y,nz+1)=plane(y,1)
      end do
c
c     Make grid
c
c      if (fltype.ne.3) then
      do y=1,nyp
         do z=1,nz+1
            grid(nyp+1-y,z,1)=real(z-nz/2-1)/real(nz)*zl
            grid(nyp+1-y,z,2)=eta(y)
         end do
      end do
c      else
c        do 2040 y=1,nyp
c        do 2040 z=1,nz+1
c          grid(nyp+1-y,z,1)=real(z-nz/2-1)/real(nz)*zl
c          grid(nyp+1-y,z,2)=(eta(y)+1.)/dstar
c 2040   continue
c      end if

      return

      end subroutine getyzs
