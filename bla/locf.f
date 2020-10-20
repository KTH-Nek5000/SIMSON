c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine locf(om2r,om2i,yb,xl,zl,xsc,zsc,eta,tc,
     &     loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,fpds8,
     &     fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,g1,g2,
     &     th2r,th2i,u2r,u2i)

c
c     Localized forcing
c
c==== loctyp=1:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z)*f(t)
c
c     zscale>0   g(x,z)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)
c     zscale<0   g(x,z)=exp(-(x-xloc0)/xscale**2)*cos((z-x*lskew)/zscale*2*pi)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     where step is defined in step.f
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c
c==== loctyp=2:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz).*(g(1,z't),g(2,z't),g(3,z't))**fy(y),fx(x)
c
c     g(1,z't)=cos(zbet*z)*cos(tomeg*t)/(2*tomeg)
c     g(2,z't)=cos(zbet*z)*sin(tomeg*t)
c     g(3,z't)=-sin(zbet*z)*sin(tomeg*t)/(2*zbet)
c
c     fx(x)=step((x-xstart)/xrise)-step((x-xend)/xfall+1)
c
c     fy(y)=step((y-ystart)/yrise)-step((y-yend)/yfall+1)
c
c     where step is defined in step.f
c
c==== loctyp=3:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*fx(x)*f(t)
c
c     xtype=0  f(x)=                 exp(-((x-xloc0)/xscale)**2)
c     xtype=1  f(x)=(x-xloc0)/xscale*exp(-((x-xloc0)/xscale)**2)
c
c     f(t)=sin(tomeg*t)
c
c==== loctyp=4:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-y/yscale)*ft(t)
c
c     f(t)=(step((t-tstart)/tscale))-step((t-tend)/tscale+1))*cos(tomeg*t)
c
c==== loctyp=5:
c
c     Related to loctyp=1
c
c     adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z)*f(t)
c
c
c     zscale>0   g(x,z)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)
c     zscale<0   g(x,z)=cos((z-x*lskew)/zscale*2*pi)*exp(-(x-xloc0)/xscale**2)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     h1(t)=aomeg*sin(tomeg*tc) (note: aomeg is relative amplitude
c     between oscillation and stationary force)
c
c     where step is defined in step.f
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c==== loctyp=6:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*fx(x)*f(t)
c     useful for TS waves and corresponding secondary instability (K/H-type)
c
c     f(x)=exp(-((x-xloc0)/xscale)**2)
c
c     g(z) = cos(2pi/zl)
c
c     f2d(t)   = sin(tomeg*t)
c     f3d(t) = sin(tomeg3D*t)
c
c     F=(0,1,0)*exp(-(yc/yscale)**2)*(amp2d*f2d(t)*f(x) +
c                                     amp3d*f3d(t)*f(x)*g(z) )
c
c==== loctyp=7:
c
c     Adding a localised forcing of the temperature (line source)
c
c==== loctyp=8:
c
c     Approximated an impulse input
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-((y-yloc)/yscale)**2)*fx(x)*f(t)
c
c     xtype=0  f(x)=                 exp(-((x-xloc0)/xscale)**2)
c
c     f(t)=exp(-((t-tstart)/tscale)**2)
c
c==== loctyp=10:
c
c     Variation of loctyp 1:
c     Localized force corresponding to haramonic force described 
c     in Högberg & Henningson, 1998
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z)*f(t)*h1(t)
c
c     zscale>0   g(x,z)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)
c     zscale<0   g(x,z)=exp(-(x-xloc0)/xscale**2)*cos((z-x*lskew)/zscale*2*pi)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     where step is defined in step.f
c
c     h1(t)=cos(tomega*tc)
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c==== loctyp=11:
c
c     Variation of loctyp 5:
c     Localized force corresponding to stationary and time-dependent 
c     force described in Högberg & Henningson, 1998
c
c     Adding a stationary and a time dependen force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z,t)*f(t)
c
c     zscale>0   g(x,z,t)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)+
c                         exp(-(x-xo)/xscale**2-(z/zscale)**2)*h1(t)
c     zscale<0   g(x,z,t)=exp(-(x-xloc0)/xscale**2)*
c              cos((z-x*lskew)/zscale*2*pi)+
c                     exp(-(x-xo)/xscale**2)*cos((z-x*lskew)/zscale*2*pi)*h1(t)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     where step is defined in step.f
c
c     h1(t)=aomeg*cos(tomega*tc)
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c==== loctyp=12:
c
c     Cylinder roughness: F=A*chi(r,y)*f(t)
c     

      implicit none

      include 'par.f'

      integer yb,loctyp,ith
      real om2r(nxp/2+1,mby,nzd,3),om2i(nxp/2+1,mby,nzd,3)
      real th2r(nxp/2+1,mby,nzd,4*scalar),th2i(nxp/2+1,mby,nzd,4*scalar)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real xl,zl,xsc,zsc,tc
      real eta(nyp)
      real g1(nxp/2,nzd),g2(nxp/2,nzd)
      real fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5

      real ampx,ampy,ampz,xscale,yscale,zscale,tscale,xloc0,lskew
      real xtype,xstart,xend,xrise,xfall,ystart,yend,yrise,yfall
      real zbet,tomeg,tstart,tend,xo,aomeg,y0,yscale1
      real amp,yloc0

      real xc1(nxp/2),xc2(nxp/2)
      integer x,y,z
      real xc11,xc22,fx1(nxp/2),fx2(nxp/2),f2x1(nxp/2),f2x2(nxp/2)
      real fy,dfy,f2y,ft,f2t,yc,zc,k2,h1
      real pi
      real amp2d,amp3d,tomeg3d,ft3d
      real rad,lam,lam0
      real radius,h,smooth,x1center,x2center,chi1,chi2
      integer iplane
      parameter (pi = 3.1415926535897932385)

      real mm2r(nzd,3),bbb
c
c     Functions
c
      real,external :: step,dstep,chi

c
c Note that we need a coordinate that runs between -xl/2 and xl/2
c regardless of the shift xsc, otherwise the force would be turned off
c abruptly when shifted out of the box
c

      if (loctyp.eq.1) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5

         if (tscale.gt.0..and.tc.gt.5.*tscale) return

         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do

c
c     Construct g(x,z)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*fx1(x)
                  g2(x,z)=exp(-(zc/zscale)**2)*fx2(x)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)
               end do
            end do
        end if
c
c     Construct f(t)
c
        if (tscale.gt.0.) ft = exp(-(tc/tscale)**2)
        if (tscale.lt.0.) ft = step(-(tc/tscale))
        if (tscale.eq.0.) ft=1.

        do z=1,nzpc
           do y=1,min(mby,nyp-yb+1)
              yc=1.+eta(y+yb-1)
              fy=exp(-(yc/yscale)**2)
              do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                 om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*g1(x,z)
                 om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*g2(x,z)
                 om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*g1(x,z)
                 om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*g2(x,z)
                 om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*g1(x,z)
                 om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*g2(x,z)
              end do
          end do
       end do
       

      end if

      if (loctyp.eq.2) then
c
c     Stellans note: disturbance slightly changed since wiegel and
c     turb simulations.
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xstart=fpds1
         xend=fpds2
         xrise=fpds3
         xfall=fpds4
         ystart=fpds5
         yend=fpds6
         yrise=fpds7
         yfall=fpds8
         zbet=fpdds4
         tomeg=fpdds5
         k2=sqrt(tomeg*tomeg+zbet*zbet)
c
c     Construct fx(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x) = step((xc1(x)-xstart)/xrise)-
     &           step((xc1(x)-xend)/xfall+1)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x) = step((xc2(x)-xstart)/xrise)-
     &           step((xc2(x)-xend)/xfall+1)
         end do
c
c     Construct g(i,z)
c
         do z=1,nzpc
            zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
            g1(1,z)=-tomeg*cos(zbet*zc)*sin(tomeg*tc)/k2
            g1(2,z)=cos(zbet*zc)*cos(tomeg*tc)
            g1(3,z)=-zbet*sin(zbet*zc)*cos(tomeg*tc)/k2
         end do

         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            fy=step((yc-ystart)/yrise)-step((yc-yend)/yfall+1)
            dfy=dstep((yc-ystart)/yrise)/yrise
     &           -dstep((yc-yend)/yfall+1)/yfall
            do z=1,nzpc
               do x=1,nxp/2
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*g1(1,z)*dfy*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*g1(1,z)*dfy*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*g1(2,z)*fy *fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*g1(2,z)*fy *fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*g1(3,z)*dfy*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*g1(3,z)*dfy*fx2(x)
               end do
            end do
         end do
      end if

      if (loctyp.eq.3) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         xtype=fp1
         tomeg=fpdds4

c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=((xc11-xloc0)/xscale)**xtype
     &           *exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=((xc22-xloc0)/xscale)**xtype
     &           *exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft=sin(tc*tomeg)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-(yc/yscale)**2)
c            write(*,*) yb,fy

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*fx2(x)
               end do
            end do
         end do
      end if

      if (loctyp.eq.4) then
         ampz=fpdds1
         yscale=fpds1
         tstart=fpds2
         tend=fpds3
         tscale=fpds4
         tomeg=fpdds2
c
c     Construct f(t)
c
         ft=(step((tc-tstart)/tscale)-step((tc-tend)/tscale+1))*
     &        cos(tc*tomeg)

         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            fy=exp(-yc/yscale)
            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft
               end do
            end do
         end do
      end if

      if (loctyp.eq.5) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5
         xo=fpds6
         aomeg=fpdds4
         tomeg=fpdds5
         y0=fpds7
         yscale1=fpds8

         if (tscale.gt.0..and.tc.gt.5.*tscale) return
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            f2x1(x)=exp(-((xc11-xo)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
            f2x2(x)=exp(-((xc22-xo)/xscale)**2)
         end do
c
c
c     Construct g(x,z)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*fx1(x)
                  g2(x,z)=exp(-(zc/zscale)**2)*fx2(x)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)
               end do
            end do
         end if
c
c     Construct f(t)
c
         if (tscale.gt.0.) ft=exp(-(tc/tscale)**2)
         if (tscale.lt.0.) ft=step(-(tc/tscale))
         if (tscale.eq.0.) ft=1.
c
c     Construct f2(t)
c
         f2t=aomeg*sin(tomeg*tc)

         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               yc=1.+eta(y+yb-1)
               fy=exp(-(yc/yscale)**2)
               f2y=exp(-((yc-y0)/yscale1)**2)
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
               end do
            end do
         end do
      end if

      if (loctyp.eq.6) then
c
c     Rename variables to simplify understanding
c
         amp2d=fpdds1
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         tomeg=fpdds4
         amp3d=fpdds2
         tomeg3D=fpdds3
c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct g(1,z)
c
         do z=1,nzpc
            zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
            g1(1,z) = cos(2*pi/zl*zc)
         end do
c
c     Construct f(t)
c
         ft  = sin(tc*tomeg)
         ft3d = sin(tc*tomeg3D)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-(yc/yscale)**2)

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,2)=om2r(x,y,z,2)+amp2d*fy*ft*fx1(x)+
     &                 amp3d*fy*ft3d*fx1(x)*g1(1,z)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+amp2d*fy*ft*fx2(x)+
     &                 amp3d*fy*ft3d*fx2(x)*g1(1,z)
               end do
            end do
         end do
      end if

      if (loctyp.eq.7) then
c
c     Rename variables to simplify understanding
c
         amp=fpdds1
         xloc0=fpds1
         yloc0=fpds2
         xscale=fpds3
         yscale=fpds4
         tomeg=fpdds2
c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft  = cos(tc*tomeg)

         do ith=1,scalar
            do y=1,min(mby,nyp-yb+1)
               yc = 1.+eta(y+yb-1)
               fy = exp(-((yc-yloc0)/yscale)**2)
               do z=1,nzpc
                  do x=1,nxp/2
                     th2r(x,y,z,4+4*(ith-1))=th2r(x,y,z,4+4*(ith-1))
     &                    +amp*fy*ft*fx1(x)
                     th2i(x,y,z,4+4*(ith-1))=th2i(x,y,z,4+4*(ith-1))
     &                    +amp*fy*ft*fx2(x)
                  end do
               end do
            end do
         end do
      end if

      if (loctyp.eq.8) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         yloc0=fpds6
         tscale=fpds4
         tstart=fpds5
c
c     Construct f(x)
c
         xtype = 0.
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=((xc11-xloc0)/xscale)**xtype
     &           *exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=((xc22-xloc0)/xscale)**xtype
     &           *exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft=exp(-((tc-tstart)/tscale)**2)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-((yc-yloc0)/yscale)**2)

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*fx2(x)
               end do
            end do
         end do
      end if

      if (loctyp.eq.9) then

         do y=1,min(mby,nyp-yb+1)
            yc = eta(y+yb-1)
            do z=1,nzpc
               do x=1,nxp/2

                  lam0 = 0.8

                  rad = sqrt((real(2*x-1)/real(nxp)*xl-xl/2.)**2 +
     &                 (real(z)/real(nzpc)*zl-zl/2.)**2)

                  lam = step( (rad-170.)/20. )

                  om2r(x,y,z,1)=om2r(x,y,z,1)-lam*(u2r(x,y,z,1)-yc)
                  om2r(x,y,z,2)=om2r(x,y,z,2)-lam*(u2r(x,y,z,2))
                  om2r(x,y,z,3)=om2r(x,y,z,3)-lam*(u2r(x,y,z,3))

                  rad = sqrt((real(2*x)/real(nxp)*xl-xl/2.)**2 +
     &                 (real(z)/real(nzpc)*zl-zl/2.)**2)

                  lam = step( (rad-170.)/20. )

                  om2i(x,y,z,1)=om2i(x,y,z,1)-lam*(u2i(x,y,z,1)-yc)
                  om2i(x,y,z,2)=om2i(x,y,z,2)-lam*(u2i(x,y,z,2))
                  om2i(x,y,z,3)=om2i(x,y,z,3)-lam*(u2i(x,y,z,3))

               end do
            end do
         end do
      end if

      if (loctyp.eq.10) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5
         tomeg=fpds6

         if (tscale.gt.0..and.tc.gt.5.*tscale) return

         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct g(x,z)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*fx1(x)
                  g2(x,z)=exp(-(zc/zscale)**2)*fx2(x)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)
               end do
            end do
        end if
c
c     Construct f(t)
c
        if (tscale.gt.0.) ft = exp(-(tc/tscale)**2)
        if (tscale.lt.0.) ft = step(-(tc/tscale))
        if (tscale.eq.0.) ft=1.

        do z=1,nzpc
           do y=1,min(mby,nyp-yb+1)
              yc=1.+eta(y+yb-1)
              fy=exp(-(yc/yscale)**2)
              do x=1,nxp/2
c
c     Construct f2(t)
c
         h1=cos(tomeg*tc)
c
c     Volume force, note that it is nonzero on the wall
c
                 om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*g1(x,z)*h1
                 om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*g2(x,z)*h1
                 om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*g1(x,z)*h1
                 om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*g2(x,z)*h1
                 om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*g1(x,z)*h1
                 om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*g2(x,z)*h1
              end do
           end do
        end do

      end if

      if (loctyp.eq.11) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5
         xo=fpds6
         tomeg=fpdds4
         aomeg=fpdds5
         
         if (tscale.gt.0..and.tc.gt.5.*tscale) return
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            f2x1(x)=exp(-((xc11-xo)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
            f2x2(x)=exp(-((xc22-xo)/xscale)**2)
         end do
c
c     Construct f(t)
c

         if (tscale.gt.0.) ft=exp(-(tc/tscale)**2)
         if (tscale.lt.0.) ft=step(-(tc/tscale))
         if (tscale.eq.0.) ft=1.
c     
c     Construct h1(t)
c
         h1=aomeg*cos(tomeg*tc)
c
c
c     Construct g(x,z,t)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*(fx1(x)+f2x1(x)*h1)
                  g2(x,z)=exp(-(zc/zscale)**2)*(fx2(x)+f2x2(x)*h1)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)+cos((zc-(xc1(x)-xo)*lskew)/zscale*2.*pi)
     &                 *f2x1(x)*h1
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)+cos((zc-(xc2(x)-xo)*lskew)/zscale*2.*pi)
     &                 *f2x2(x)*h1
               end do
            end do
         end if
         
         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               yc=1.+eta(y+yb-1)
               fy=exp(-(yc/yscale)**2)
               do x=1,nxp/2
c     
c     Volume force, note that it is nonzero on the wall
c     
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*g1(x,z)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*g2(x,z)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*g1(x,z)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*g2(x,z)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*g1(x,z)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*g2(x,z)
               end do
            end do
         end do
      end if
c     
c     Large-scale vortices
c
      if (loctyp.eq.12) then
c
c     Rename variables
c
         amp=fpdds1
         ampx=fpdds2
         ampy=fpdds3
         zscale=fpds1
         tscale=fpds2
         
         if (tscale.gt.0..and.tc.gt.5.*tscale) return
         
         if (tscale.gt.0.) ft = exp(-(tc/tscale)**2)
         if (tscale.lt.0.) ft = step(-(tc/tscale))
         if (tscale.eq.0.) ft = 1.
         

         if (1.eq.1) then
c
c     Add the vortex as a force in the whole domain
c     forced1: amp=1, ampx=0, ampy=amplitude
c     forced2: amp=1, ampx=1, ampy=amplitude
c
            bbb = 2*pi/zscale
            do x=1,nxp/2
               xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
               xc11=xc11-int((xc11+xl/2.)/xl)*xl
               xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
               xc22=xc22-int((xc22+xl/2.)/xl)*xl
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  do y=1,min(mby,nyp-yb+1)
                     
                     yc=eta(y+yb-1)

c                     if (abs(yc).gt.0.5) then

c     pressure gradient in x direction
c                     om2r(x,y,z,2)=om2r(x,y,z,2)-amp*(ampx*u2r(x,y,z,2)
c     &                    -ampy/2*cos(zc*bbb)*(1+cos(pi*yc)))
c                     om2r(x,y,z,3)=om2r(x,y,z,3)-amp*(ampx*u2r(x,y,z,3)
c     &                    -ampy/2*pi/bbb*sin(zc*bbb)*sin(pi*yc))
c                     
c                     om2i(x,y,z,2)=om2i(x,y,z,2)-amp*(ampx*u2i(x,y,z,2)
c     &                    -ampy/2*cos(zc*bbb)*(1+cos(pi*yc)))
c                     om2i(x,y,z,3)=om2i(x,y,z,3)-amp*(ampx*u2i(x,y,z,3)
c     &                    -ampy/2*pi/bbb*sin(zc*bbb)*sin(pi*yc))

c     pressure gradient in z direction
                     om2r(x,y,z,2)=om2r(x,y,z,2)-amp*(ampx*u2r(x,y,z,2)
     &                    -ampy/2*cos(xc11*bbb)*(1+cos(pi*yc)))
                     om2r(x,y,z,1)=om2r(x,y,z,1)-amp*(ampx*u2r(x,y,z,1)
     &                    -ampy/2*pi/bbb*sin(xc11*bbb)*sin(pi*yc))
                     
                     om2i(x,y,z,2)=om2i(x,y,z,2)-amp*(ampx*u2i(x,y,z,2)
     &                    -ampy/2*cos(xc22*bbb)*(1+cos(pi*yc)))
                     om2i(x,y,z,1)=om2i(x,y,z,1)-amp*(ampx*u2i(x,y,z,1)
     &                    -ampy/2*pi/bbb*sin(xc22*bbb)*sin(pi*yc))
                     
c                  end if

                  end do
               end do
            end do
         else  
c     
c     This is to force the mean
c     forced3: amp=1, ampx=1 , ampy=amplitude
c
            mm2r = 0
            do z=1,nzpc
               do x=1,nxp/2
                  mm2r(z,1) = mm2r(z,1) + u2r(x,1,z,1)
     &                 + u2i(x,1,z,1)
                  mm2r(z,2) = mm2r(z,2) + u2r(x,1,z,2)
     &                 + u2i(x,1,z,2)
                  mm2r(z,3) = mm2r(z,3) + u2r(x,1,z,3)
     &                 + u2i(x,1,z,3)
               end do
            end do
            mm2r = mm2r/nxp
            bbb = 2*pi/zscale
c
c     the forcing is in principle independent of x
c
            do x=1,nxp/2
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  do y=1,min(mby,nyp-yb+1)
                     
                     yc=eta(y+yb-1)

                     om2r(x,y,z,2)=om2r(x,y,z,2)-amp*(ampx*mm2r(z,2)
     &                    -ampy/2*cos(bbb*zc)*(1+cos(pi*yc)))
                     om2r(x,y,z,3)=om2r(x,y,z,3)-amp*(ampx*mm2r(z,3)
     &                    -ampy/2*pi/bbb*sin(bbb*zc)*sin(pi*yc))
                     
                     om2i(x,y,z,2)=om2i(x,y,z,2)-amp*(ampx*mm2r(z,2)
     &                    -ampy/2*cos(bbb*zc)*(1+cos(pi*yc)))
                     om2i(x,y,z,3)=om2i(x,y,z,3)-amp*(ampx*mm2r(z,3)
     &                    -ampy/2*pi/bbb*sin(bbb*zc)*sin(pi*yc))
                     
                  end do
               end do
            end do
         end if
      end if
c
c     LEBU
c
      if (loctyp.eq.13) then
c
c     Rename variables
c
         amp=fpdds1
         h=fpds1
         radius=fpds2
         smooth=fpds3
         x1center=fpds4
         x2center=fpds5
         iplane=int(fp1+0.5)
         tscale=fpds6

c        if (tscale.gt.0..and.tc.gt.5.*tscale) return

c        if (tscale.gt.0.) ft = exp(-(tc/tscale)**2)
c         if (tscale.lt.0.) ft = step(-(tc/tscale))
         if (tscale.eq.0.) ft = 1.

         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            
            do z=1,nzpc
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
               do y=1,min(mby,nyp-yb+1)
                  yc=1.+eta(y+yb-1)
c
c                 Cylinder base in xz-plane
c
                  if (iplane.eq.0) then
                     chi1 = chi(h,radius,smooth,xc11,zc,yc,
     &                    x1center,x2center)
                     chi2 = chi(h,radius,smooth,xc22,zc,yc,
     &                    x1center,x2center)
c     
c                 Cylinder base in xy-plane
c     
                  elseif (iplane.eq.1) then
                     chi1 = chi(h,radius,smooth,xc11,yc,zc,
     &                    x1center,x2center)
                     chi2 = chi(h,radius,smooth,xc22,yc,zc,
     &                    x1center,x2center)
                  end if
                  
                  om2r(x,y,z,1)=om2r(x,y,z,1)-amp*ft*u2r(x,y,z,1)*chi1
                  om2r(x,y,z,2)=om2r(x,y,z,2)-amp*ft*u2r(x,y,z,2)*chi1
                  om2r(x,y,z,3)=om2r(x,y,z,3)-amp*ft*u2r(x,y,z,3)*chi1
                  
                  om2i(x,y,z,1)=om2i(x,y,z,1)-amp*ft*u2i(x,y,z,1)*chi2
                  om2i(x,y,z,2)=om2i(x,y,z,2)-amp*ft*u2i(x,y,z,2)*chi2
                  om2i(x,y,z,3)=om2i(x,y,z,3)-amp*ft*u2i(x,y,z,3)*chi2

               end do
            end do
         end do
      end if



      if (loctyp.gt.13) then
         write(*,*) 'loctyp = ',loctyp,' not implemented.'
         call stopnow(5433332)
      end if

      end subroutine locf


      real function chi(xLB,yLB,s,x1,y1,x3,x1c,x2c)
c
c Smooth step function:
c
c        / 0                                                 x3>=h+s
c        |                                                   r>=r0+s
c        |
c        | 0.5*(1+cos((r-r0)/s*pi))*0.5*(1+cos((x3-h)/s*pi)  h<x3<h+s
c        |                                                   r0<r<r0+s
c        |
c chi = <  0.5*(1+cos((r-r0)/s*pi))                          x3<=h
c        |                                                   r0<r<r0+s
c        |
c        | 0.5*(1+cos((x3-h)/s*pi))                          h<x3<h+s
c        |                                                   r<=r0
c        |
c        \ 1                                                 x3<=h; r<=r0
c
c
c Coordinate system:
c
c                x3
c                ^
c                |   __|__
c                |  /  .  \
c                |  \__|__/ 
c                |  |  .  | 
c                /--|  |  |------> x1
c               /   |  .  |
c              /    |  |  |
c             /     \__.__/
c            /         |
c           V
c          x2            
c

      implicit none
      
      real s,xLB,yLB,x1,y1,x3,x1c,x2c
      
      real pi
      parameter (pi = 3.1415926535897932385)
      
      if (x1.le.(x1c-xLB/2-s).or.x1.ge.(x1c+xLB/2+s)) then
         chi=0.     
      elseif (x1.gt.(x1c-xLB/2-s).and.x1.lt.(x1c-xLB/2)) then
         if (y1.le.(x2c-yLB/2-s)) then
            chi=0.
         elseif (y1.gt.(x2c-yLB/2-s).and.y1.le.(x2c-yLB/2)) then
            chi=(-cos((x1-x1c+xLB/2+s)/s*pi)/2+0.5)*
     &           (-cos((y1-x2c+yLB/2+s)/s*pi)/2+0.5)
         elseif (y1.gt.(x2c-yLB/2).and.y1.lt.(x2c+yLB/2)) then
            chi=-cos((x1-x1c+xLB/2+s)/s*pi)/2+0.5   
         elseif (y1.ge.(x2c+yLB/2).and.y1.le.(x2c+yLB/2+s)) then
            chi=(-cos((x1-x1c+xLB/2+s)/s*pi)/2+0.5)*
     &           (cos((y1-x2c-yLB/2)/s*pi)/2+0.5)
         else
            chi=0.
         endif
      elseif (x1.gt.(x1c-xLB/2).and.x1.lt.(x1c+xLB/2)) then
         if (y1.le.(x2c-yLB/2-s)) then
            chi=0.
         elseif (y1.gt.(x2c-yLB/2-s).and.y1.le.(x2c-yLB/2)) then
            chi=-cos((y1-x2c+yLB/2+s)/s*pi)/2+0.5
         elseif (y1.gt.(x2c-yLB/2).and.y1.lt.(x2c+yLB/2)) then
            chi=1.
         elseif (y1.ge.(x2c+yLB/2).and.y1.le.(x2c+yLB/2+s)) then
            chi=cos((y1-x2c-yLB/2)/s*pi)/2+0.5
         else
            chi=0.
         endif
      elseif (x1.ge.(x1c+xLB/2).and.x1.lt.(x1c+xLB/2+s)) then
         if (y1.le.(x2c-yLB/2-s)) then
            chi=0.
         elseif (y1.gt.(x2c-yLB/2-s).and.y1.le.(x2c-yLB/2)) then
            chi=(cos((x1-x1c-xLB/2)/s*pi)/2+0.5)*
     &           (-cos((y1-x2c+yLB/2+s)/s*pi)/2+0.5)        
         elseif (y1.gt.(x2c-yLB/2).and.y1.lt.(x2c+yLB/2)) then
            chi=cos((x1-x1c-xLB/2)/s*pi)/2+0.5   
         elseif (y1.ge.(x2c+yLB/2).and.y1.le.(x2c+yLB/2+s)) then
            chi=(cos((x1-x1c-xLB/2)/s*pi)/2+0.5)*
     &           (cos((y1-x2c-yLB/2)/s*pi)/2+0.5)
         else
            chi=0.
         endif
      else
         chi=0.
      endif
      
      end function chi

