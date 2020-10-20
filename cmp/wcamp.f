c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wcamp(t,amp,eta,wint,re,dstar,fltype,u0upp,u0low)
c
c     amp(y,1) streamwise velocity average without (0,0) comp. (squared)
c     amp(y,2) normal velocity  average without  (squared)
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
      implicit none

      include 'par.f'

      real t,eta(nyp),wint(nyp),re,dstar
      real amp(nyp,20)
      integer fltype
      real sum(100),hplus,e0,um,omm,dum
      real u0upp,u0low
      integer i,j,y

      u0upp=amp(1,9)
      u0low=amp(nyp,9)
      do j=1,20
         sum(j)=0.
      end do
c
c     This loop makes the amplitude files "compatible" with previous ones
c     by adding in the (0,0) component
c     in addition the turbulence production is calculated
c
      do y=1,nyp
         um=0.0
         omm=0.0
         if (fltype.eq.1.or.fltype.eq.4) um=1.-eta(y)**2+u0upp
         if (fltype.eq.2.or.fltype.eq.5) um=(u0upp-u0low)*.5*eta(y)
         if (fltype.eq.1.or.fltype.eq.4)  omm=2.*eta(y)
         if (fltype.eq.2.or.fltype.eq.5) omm=-(u0upp-u0low)*.5
         dum=-omm
         amp(y,1)=amp(y,1)+(amp(y,9)-um)**2
         amp(y,3)=amp(y,3)+amp(y,10)**2
         amp(y,4)=amp(y,4)+amp(y,11)**2
         amp(y,6)=amp(y,6)+(amp(y,12)-omm)**2
         amp(y,7)=amp(y,7)+(amp(y,9)-um)**2+amp(y,10)**2
         amp(y,8)=-dum*amp(y,8)
         amp(y,9)=(amp(y,9)-um)**2
         amp(y,10)=amp(y,10)**2
      end do
      do i=1,10
         do y=1,nyp
            sum(i)=sum(i)+amp(y,i)*wint(y)
         end do
      end do
c
c     Then scale
c
      do j=1,7
         sum(j)=sqrt(sum(j)*.5)
      end do
      e0=sqrt((sum(9)+sum(10))*.5)
      hplus=sqrt((abs(amp(1,12))+abs(amp(nyp,12)))*re/2.)
      if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0)
     &     hplus=sqrt(abs(amp(nyp,12))*re)
      write(*,*) 't',t/dstar
      write(*,*) 'velocity rms  ',(sum(j),j=1,3)
      write(*,*) 'unorm',sqrt(sum(1)**2+sum(2)**2+sum(3)**2)
      write(*,*) 'vorticity rms  ',(sum(j)*dstar,j=4,6)
      write(*,*) 'omy**2/k2, dUuv, e0, h+ ',sum(7),sum(8),e0,hplus

      return

      end subroutine wcamp
