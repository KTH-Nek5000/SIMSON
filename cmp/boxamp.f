c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine boxamp(amp,campw,kx,kz,nwave,u2r,u2i,om2r,om2i,yb,
     &                  alfa,beta)
c
c     Accumulates the amplitude from an xz-box
c
c     This version accumulates u(i)**2,om(i)**2
c
      implicit none

      include 'par.f'

      integer yb,nwave,kx(nwave),kz(nwave)
      complex campw(nyp,4,nwave)
      real amp(nyp,20)
      real u2r((nxp/2+1),mby,nzd,3),u2i((nxp/2+1),mby,nzd,3)
      real om2r((nxp/2+1),mby,nzd,3),om2i((nxp/2+1),mby,nzd,3)
      real alfa(nx/2),beta(nz)
c
      integer i,z,y,x,zp,y1
      real c
c
c     amp(y,1) streamwise velocity average without (0,0) comp. (squared)
c     amp(y,2) normal velocity  average (squared)
c     amp(y,3) spanwise velocity average without (0,0) (squared)
c     amp(y,4) streamwise vorticity average without (0,0) comp. (squared)
c     amp(y,5) normal vorticity average (squared)
c     amp(y,6) spanwise vorticity average without (0,0) comp.(squared)
c     amp(y,7) normal vorticity squared over wavenumber square average, no (0,0)
c     amp(y,8) reynolds stress average
c     amp(y,9) mean streamwise disturbance velocity
c     amp(y,10) mean spanwise disturbance velocity
c     amp(y,11) mean streamwise vorticity component
c     amp(y,12) mean spanwise vorticity component (to calculate wall shear)
c     amp(y,13-20) empty
c
c     campw(y,i,wave) complex normal velocity and normal vorticity averages
c     from selected wavenumbers
c
      c=1.
      if (nfzsym.eq.1) c=2.
      do i=1,3
         do y=yb,min(nyp,yb+mby-1)
            y1=y-yb+1
            amp(y,i)=0.
            amp(y,i+3)=0.
            do x=2,nx/2
               amp(y,i)=amp(y,i)
     &              +(u2r(x,y1,1,i)**2+u2i(x,y1,1,i)**2)*2
               amp(y,i+3)=amp(y,i+3)
     &              +(om2r(x,y1,1,i)**2+om2i(x,y1,1,i)**2)*2
            end do
            do z=2,nz/2
               amp(y,i)=amp(y,i)
     &                 +(u2r(1,y1,z,i)**2+u2i(1,y1,z,i)**2)*c
               amp(y,i+3)=amp(y,i+3)
     &              +(om2r(1,y1,z,i)**2+om2i(1,y1,z,i)**2)*c
            end do
            do z=2,nz/2
               do x=2,nx/2
                  amp(y,i)=amp(y,i)
     &                 +(u2r(x,y1,z,i)**2+u2i(x,y1,z,i)**2)*c*2
                  amp(y,i+3)=amp(y,i+3)
     &                 +(om2r(x,y1,z,i)**2+om2i(x,y1,z,i)**2)*c*2
               end do
            end do
            if (nfzsym.eq.0) then
               do z=nzp-nz/2+2,nzp
                  amp(y,i)=amp(y,i)+u2r(1,y1,z,i)**2+u2i(1,y1,z,i)**2
                  amp(y,i+3)=amp(y,i+3)+om2r(1,y1,z,i)**2+
     &                 om2i(1,y1,z,i)**2
               end do
               do z=nzp-nz/2+2,nzp
                  do x=2,nx/2
                     amp(y,i)=amp(y,i)
     &                    +(u2r(x,y1,z,i)**2+u2i(x,y1,z,i)**2)*2
                     amp(y,i+3)=amp(y,i+3)
     &                    +(om2r(x,y1,z,i)**2+om2i(x,y1,z,i)**2)*2
                  end do
               end do
            end if
         end do
      end do
c
c     Accumulate amplitude for omy**2/k**2, Reynolds stress and mean
c
      do y=yb,min(nyp,yb+mby-1)
         y1=y-yb+1
         amp(y,7)=0.
         amp(y,8)=0.
         amp(y,9)=u2r(1,y1,1,1)
         amp(y,10)=u2r(1,y1,1,3)
         amp(y,11)=om2r(1,y1,1,1)
         amp(y,12)=om2r(1,y1,1,3)
         do x=2,nx/2
            amp(y,7)=amp(y,7)+
     &           (om2r(x,y1,1,2)**2+om2i(x,y1,1,2)**2)/alfa(x)**2*2
            amp(y,8)=amp(y,8)+2*(
     &           u2r(x,y1,1,1) * u2r(x,y1,1,2) +
     &           u2i(x,y1,1,1) * u2i(x,y1,1,2) )
         end do
         do z=2,nz/2
            amp(y,7)=amp(y,7)+
     &           (om2r(1,y1,z,2)**2+om2i(1,y1,z,2)**2)/beta(z)**2*c
            amp(y,8)=amp(y,8)+2*c*(
     &           u2r(1,y1,z,1) * u2r(1,y1,z,2) +
     &           u2i(1,y1,z,1) * u2i(1,y1,z,2) )
         end do
         do z=2,nz/2
            do x=2,nx/2
               amp(y,7)=amp(y,7)+(om2r(x,y1,z,2)**2+om2i(x,y1,z,2)**2)/
     &              (alfa(x)**2+beta(z)**2)*c*2.
               amp(y,8)=amp(y,8)+2*c*(
     &              u2r(x,y1,z,1) * u2r(x,y1,z,2) +
     &              u2i(x,y1,z,1) * u2i(x,y1,z,2) )
            end do
         end do
         if (nfzsym.eq.0) then
            do z=nz/2+2,nz
               zp=nzp-nz+z
               amp(y,7)=amp(y,7)+
     &              (om2r(1,y1,zp,2)**2+om2i(1,y1,zp,2)**2)/beta(z)**2
               amp(y,8)=amp(y,8)+
     &              u2r(1,y1,zp,1) * u2r(1,y1,zp,2) +
     &              u2i(1,y1,zp,1) * u2i(1,y1,zp,2)
            end do
            do z=nz/2+2,nz
               do x=2,nx/2
                  zp=nzp-nz+z
                  amp(y,7)=amp(y,7)+(om2r(x,y1,zp,2)**2+
     &                 om2i(x,y1,zp,2)**2)/
     &                 (alfa(x)**2+beta(z)**2)*2.
                  amp(y,8)=amp(y,8)+2*(
     &                 u2r(x,y1,zp,1) * u2r(x,y1,zp,2) +
     &                 u2i(x,y1,zp,1) * u2i(x,y1,zp,2) )
               end do
            end do
         end if
      end do
c
c     Accumulate energy for selected wavenumbers
c
        do i=1,nwave
           if (kx(i).lt.0) goto 1000
 2070      continue
           do y=yb,min(nyp,yb+mby-1)
              x=kx(i)+1
              y1=y-yb+1
              z=kz(i)+1
              if (z.le.0) z=z+nzp
c
c     Note that here we don't sum contributions from negative and postive kz
c
              campw(y,1,i)=cmplx(u2r(x,y1,z,1),u2i(x,y1,z,1))
              campw(y,2,i)=cmplx(u2r(x,y1,z,2),u2i(x,y1,z,2))
              campw(y,3,i)=cmplx(u2r(x,y1,z,3),u2i(x,y1,z,3))
              campw(y,4,i)=cmplx(om2r(x,y1,z,2),om2i(x,y1,z,2))
           end do
        end do
 1000 continue

      return

      end subroutine boxamp
