c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine locsc(utau,uinf,delta,deltar,beta,pxys,pxys2,pxysth,
     &     re,dstar,xl,nux,thdelta,pr,scalarind,fltype,tinf,twall)
c
c     Find local skin friction, displacement thickness and momentum
c     thickness, shape factor
c
      implicit none

      include 'par.f'

      integer nxys2,scalarind
      parameter (nxys2=66)
      real utau(nx),uinf(nx),duinf(nx),delta(nx,7)
      real nux(nx),tinf(nx),twall(nx),thdelta(nx,2)
      real deltar(nx),beta(nx)
      real pxys(nx,nyp,nxys),pxysth(nx,nyp,nxysth)
      real pxys2(nx,nyp,nxys2)
      real re,dstar,xl,pr(scalar)

      integer x,y,fltype
      real pxy(nx,nyp),pxy2(nx,nyp),wxy(nx,nyp),txy(nx,nyp),txy2(nx,nyp)
      real prey(nyp*2+15),mflux

      call vcosti(nyp,prey,0)
c
c     Find skin friction
c
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=pxys(x,y,1)
         end do
      end do
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)
      call rdcheb(pxy,nyp,nx,nx)
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         if (pxy(x,1).gt.0) then
            utau(x)=-sqrt(pxy(x,1)*dstar*(2./real(nyp-1))/re)
         else
            utau(x)=sqrt(-pxy(x,1)*dstar*(2./real(nyp-1))/re)
         end if
         if (fltype.eq.1) then
c
c     For channel flow, take u_centreline
c
            uinf(x)=pxys(x,(nyp+1)/2,1)
         else
c
c     otherwise take uinf
c
            uinf(x)=pxys(x,nyp,1)
         end if
      end do
c
c     Compute bulk velocity
c
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=pxys(x,y,1)
         end do
      end do
      do x=1,nx
         delta(x,1)=0.0
      end do
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)

      mflux = 0.
      do y=1,nyp,2
         mflux=mflux-pxy(1,y)*2./real((y-1)**2-1)
      end do
      mflux = mflux * (2./real(nyp-1))

      call icheb(wxy,pxy,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      x=1
      write(*,*) 'Re_bulk=',
     &     (pxy(x,1)-pxy(x,nyp))*(2./real(nyp-1))*re/2.,mflux*re/2





c
c     Compute displacement and momentum thickness
c
      do y=2,nyp
         do x=1,nx
            pxy2(x,y)=(1.-pxys(x,y,1)/uinf(x))*pxys(x,y,1)/uinf(x)
            pxy(x,y)=1.-pxys(x,y,1)/uinf(x)
         end do
      end do
      do x=1,nx
         pxy(x,1)=0.0
         pxy2(x,1)=0.0
      end do
c
c     Integrate to find displacement and momentum loss thicknesses
c     integration constants
c
      do x=1,nx
         delta(x,1)=0.0
      end do
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      call vchbf(pxy2,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy2,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy2(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy2,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         delta(x,1)=(pxy(x,1)-pxy(x,nyp))/dstar*(2./real(nyp-1))
         delta(x,2)=(pxy2(x,1)-pxy2(x,nyp))/dstar*(2./real(nyp-1))
         delta(x,3)=delta(x,1)/delta(x,2)
         delta(x,4)=delta(x,1)*uinf(x)/utau(x)
      end do
c
c     Compute dtheta/dx
c
      do x=1,nx
         pxy(x,1)=delta(x,2)
      end do
      call ddx(pxy2,pxy,xl)
      do x=1,nx
         delta(x,5) = pxy2(x,1)
      end do
c
c     Compute shear stress difference
c
      do y=1,nyp
         do x=1,nx
            pxy2(x,y)=pxys(x,y,4)**2-pxys(x,y,5)**2
         end do
      end do
      do x=1,nx
         pxy2(x,1)=0.0
      end do
      call vchbf(pxy2,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy2,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy2(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy2,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         pxy(x,1)=(pxy2(x,1)-pxy2(x,nyp))/dstar*(2./real(nyp-1))
      end do
      call ddx(pxy2,pxy,xl)
      do x=1,nx
         delta(x,6)=pxy2(x,1)
      end do
c
c     Compute ration between integrated production and dissipation
c
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=pxys2(x,y,89)
         end do
      end do
c
c     Integrate the production in y
c
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         thdelta(x,1)=(pxy(x,1)-pxy(x,nyp))*(2./real(nyp-1))
      end do

      do y=1,nyp
         do x=1,nx
            pxy(x,y)=pxys2(x,y,84)+pxys2(x,y,88)
         end do
      end do
c
c     Integrate the dissipation in y
c
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy,delta,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         thdelta(x,2)=(pxy(x,1)-pxy(x,nyp))*(2./real(nyp-1))
      end do

      do x=1,nx
         delta(x,7)=-thdelta(x,1)/thdelta(x,2)
      end do
c
c     Find d theta / dy at the wall
c
      if (scalar.gt.0) then
         do y=1,nyp
            do x=1,nx
               txy(x,y)=pxysth(x,y,1)
            end do
         end do
         call vchbf(txy,wxy,nyp,nx,nx,1,prey)
         call rdcheb(txy,nyp,nx,nx)
         call vchbb(txy,wxy,nyp,nx,nx,1,prey)
         do x=1,nx
            tinf(x)=pxysth(x,nyp,1)
            twall(x)=pxysth(x,1,1)
c            write(*,*) x, twall(x),tinf(x)
            nux(x)=(txy(x,1)*dstar*(2./real(nyp-1))/((twall(x)-tinf(x))
     &           *re*pr(scalarind)))
         end do
         do x=1,nx
            nux(x)=abs(nux(x))
         end do
      end if
c
c     Compute displacement and momentum thickness for the scalar
c
      if (scalar.gt.0) then
         do y=2,nyp
            do x=1,nx
               txy2(x,y)=(pxysth(x,y,1)-twall(x))
     &          *(-pxysth(x,y,1)+tinf(x))/(tinf(x)-twall(x))**2
               txy(x,y)=(-pxysth(x,y,1)+tinf(x))/(tinf(x)-twall(x))
            end do
         end do
         do x=1,nx
            txy(x,1)=0.0
            txy2(x,1)=0.0
         end do
c
c     Integrate to find displacement and momentum loss thicknesses for scalar
c     integration constants
c
         do x=1,nx
            thdelta(x,:)=0.0
         end do
         call vchbf(txy,wxy,nyp,nx,nx,1,prey)
         call icheb(wxy,txy,delta,nyp,nx,nx)
         do y=1,nyp
            do x=1,nx
               txy(x,y)=wxy(x,y)
            end do
         end do
         call vchbb(txy,wxy,nyp,nx,nx,1,prey)

         call vchbf(txy2,wxy,nyp,nx,nx,1,prey)
         call icheb(wxy,txy2,delta,nyp,nx,nx)
         do y=1,nyp
            do x=1,nx
               txy2(x,y)=wxy(x,y)
            end do
         end do
         call vchbb(txy2,wxy,nyp,nx,nx,1,prey)
         do x=1,nx
           thdelta(x,1)=(txy(x,1)-txy(x,nyp))/dstar*(2./real(nyp-1))
           thdelta(x,2)=(txy2(x,1)-txy2(x,nyp))/dstar*(2./real(nyp-1))
         end do
      end if

      do y=2,nyp
         do x=1,nx
            pxy(x,y)=pxys2(x,y,43)
         end do
      end do
      do x=1,nx
         pxy(x,1)=0.0
      end do
c
c     Integrate to find Reynolds thickness
c     integration constant
c
      do x=1,nx
         deltar(x)=0.0
      end do
      call vchbf(pxy,wxy,nyp,nx,nx,1,prey)
      call icheb(wxy,pxy,deltar,nyp,nx,nx)
      do y=1,nyp
         do x=1,nx
            pxy(x,y)=wxy(x,y)
         end do
      end do
      call vchbb(pxy,wxy,nyp,nx,nx,1,prey)
      do x=1,nx
         deltar(x)=(pxy(x,1)-pxy(x,nyp))/dstar*(2./real(nyp-1))
      end do
      call ddx(pxy2,pxys(1,1,1),xl)
      do x=1,nx
         duinf(x)=pxy2(x,nyp)
         beta(x)=-1*delta(x,1)/utau(x)/utau(x)*duinf(x)*uinf(x)
      end do

      end subroutine locsc
