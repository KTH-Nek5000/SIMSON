c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine trip(om2r,om2i,yb,xl,xsc,eta,txsc,tysc,ttzc,tx0)

      implicit none
      include 'par.f'
      integer yb
      real om2r(nxp/2+1,mby,nzd,3),om2i(nxp/2+1,mby,nzd,3)
      real txsc,tysc,ttzc(nzp+2,4),tx0
      real xl,xsc
      real eta(nyp)

      integer x,y,z
      real yc,xc11,xc22,ex1(nxp/2),ex2(nxp/2),expy
c
c     Adding a volume force of the form
c     Fv=exp(-((x-tx0)/txsc)**2-(y/tysc)**2)*f(z,t)
c
c     Note that we need a coordinate that runs between -xl/2 and xl/2
c     regardless of the shift xsc, otherwise the force would be turned off
c     abruptly when shifted out of the box
c
c     f(z,t) is computed in gtrip.f
c
c     Compute f(x)=exp(-((x-tx0)/txsc)**2)
c
      do x=1,nxp/2
         xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
         xc11=xc11-int((xc11+xl/2.)/xl)*xl
         ex1(x)=exp(-((xc11-tx0)/txsc)**2)
         xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
         xc22=xc22-int((xc22+xl/2.)/xl)*xl
         ex2(x)=exp(-((xc22-tx0)/txsc)**2)
      end do

      do z=1,nzpc
         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            expy=exp(-(yc/tysc)**2)
            do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
               om2r(x,y,z,2)=om2r(x,y,z,2)+expy*ex1(x)*ttzc(z,1)
               om2i(x,y,z,2)=om2i(x,y,z,2)+expy*ex2(x)*ttzc(z,1)
            end do
         end do
      end do
      end subroutine trip
