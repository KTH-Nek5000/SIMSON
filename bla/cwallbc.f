c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cwallbc(my_node,wallvr,wallvi,wp1,wp2,wpds1,
     &     wpds2,wpds3,wpds4,wpds5,wpds6,wpds7,wpds8,wpds9,wpds10,
     &     wpds11,wpds12,wpdds1,wpdds2,
     &     prex,prez,pres,wr,wi,xl,xsc,zl,zsc,tc,wbci,
     &     amptab,wallur,wallui,wallwr,wallwi,alfa,beta,it)
c
c     wbci=1:
c     Setting the boundary condition at the wall of v to
c     v(y=0) = amp*f(x)*cos(zbet*zc)*sin(tomeg*tc)
c
c     f(x) = step((x-xstart)/xrise)-step((x-xend)/xfall+1)
c
c     Wall oscillations are imposed for wbci=3, 4 and 5.
c
c     Note that wallur/i and wallwr/i are Dv and eta at the wall.
c
c     Note also that for the parallel version wallur etc. are still
c     too large, i.e. they could be of size memnx x memnz instead.
c
      implicit none

      include 'par.f'

      integer my_node,it
      real wallur(nxp/2+1,memnz),wallui(nxp/2+1,memnz)
      real wallvr(nxp/2+1,nzd),wallvi(nxp/2+1,nzd)
      real wallwr(nxp/2+1,memnz),wallwi(nxp/2+1,memnz)

      real wp1,wp2,wpds1,wpds2,wpds3,wpds4,wpds5
      real wpds6,wpds7,wpds8,wpds9
      real wpds10,wpds11,wpds12
      real wpdds1,wpdds2
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15)
      real wr(nxp/2+1,nzd),wi(nxp/2+1,nzd)
      real xl,xsc,zl,zsc,tc
      integer wbci

      real amptab(nxp)
      real amp1,amp2
      real xomeg,xstartc
      real alfa(nx/2*mbz),beta(nz)

      integer x,z,zb
      real zc,tmpr,tmpi
      real xc11,xc22,fx1(nxp/2),fx2(nxp/2),ft
      real amp,damp,xstart,xend,xrise,xfall,zbet,tomeg
      real zstart,zend,zrise,zfall,tstart,tend,trise,tfall
      real step
      real pi
      real dist,distom
      parameter (pi = 3.1415926535897932385)

      if (wbci.eq.1) then
         amp=wp1
         damp=wp2
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         zbet=wpdds1
         tomeg=wpdds2

         if (xend.ge.xl/2.0) then
            xend=xend-xl
            if (xstart.ge.xl/2.0) then
               xstart=xstart-xl
            end if     
         end if       
         
         if (tomeg.eq.0) then
c
c     If omega is zero, apply constant wall velocity
c
            do x=1,nxp/2
               xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
               xc11=xc11-int((xc11+xl/2.)/xl)*xl
               fx1(x)=step((xc11-xstart)/xrise)-
     &              step((xc11-xend)/xfall+1)
               xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
               xc22=xc22-int((xc22+xl/2.)/xl)*xl
               fx2(x)=step((xc22-xstart)/xrise)-
     &              step((xc22-xend)/xfall+1)
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc

                  tmpr=amp*fx1(x)*cos(zbet*zc)
                  tmpi=amp*fx2(x)*cos(zbet*zc)

                  wallvr(x,z)=max(tmpr,tmpr*damp)
                  wallvi(x,z)=max(tmpi,tmpi*damp)
               end do
            end do

         else

            do x=1,nxp/2
               xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
               xc11=xc11-int((xc11+xl/2.)/xl)*xl
               fx1(x)=step((xc11-xstart)/xrise)-
     &              step((xc11-xend)/xfall+1)
               xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
               xc22=xc22-int((xc22+xl/2.)/xl)*xl
               fx2(x)=step((xc22-xstart)/xrise)-
     &              step((xc22-xend)/xfall+1)
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc

                  tmpr=amp*fx1(x)*cos(zbet*zc)*sin(tomeg*tc)
                  tmpi=amp*fx2(x)*cos(zbet*zc)*sin(tomeg*tc)

                  wallvr(x,z)=max(tmpr,tmpr*damp)
                  wallvi(x,z)=max(tmpi,tmpi*damp)
               end do
            end do
         end if

      end if

      if (wbci.eq.2) then
         amp=wp1
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         zstart=wpds5
         zend=wpds6
         zrise=wpds7
         zfall=wpds8
         tstart=wpds9
         tend=wpds10
         trise=wpds11
         tfall=wpds12
         distom=wpdds1

         if (xend.ge.xl/2.0) then
            xend=xend-xl
            if (xstart.ge.xl/2.0) then
               xstart=xstart-xl
            end if     
         end if       

         ft=step((tc-tstart)/trise)-step((tc-tend)/tfall+1)
         dist=amp/100.*sin(distom*tc)
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            fx1(x)=step((xc11-xstart)/xrise)-step((xc11-xend)/xfall+1)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            fx2(x)=step((xc22-xstart)/xrise)-step((xc22-xend)/xfall+1)
            do z=1,nzpc
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
               tmpr=step((zc-zstart)/zrise)-step((zc-zend)/zfall+1)
               wallvr(x,z)=(amp*fx1(x)*tmpr+dist)*ft
               wallvi(x,z)=(amp*fx2(x)*tmpr+dist)*ft
            end do
         end do
      end if


      if (wbci.eq.3) then
c
c     Temporal wall oscillation
c
         amp=wp1
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         tstart=wpds9
         tend=wpds10
         trise=wpds11
         tfall=wpds12
         tomeg=wpdds2

         if (trise.eq.0.0.and.tfall.eq.0.0) then
            ft=1.0
         else
            ft=step((tc-tstart)/trise)-step((tc-tend)/tfall+1.)
         end if
         ft=ft*sin(tomeg*(tc-tstart))

         if (xend.ge.xl/2.0) then
            xend=xend-xl
            if (xstart.ge.xl/2.0) then
               xstart=xstart-xl
            end if     
         end if                  

         do x=1,nxp/2            
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            fx1(x)=step((xc11-xstart)/xrise)-
     &           step((xc11-xend)/xfall+1)

            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            fx2(x)=step((xc22-xstart)/xrise)-
     &           step((xc22-xend)/xfall+1)

            if (xstart.gt.xend) then
               fx1(x)=fx1(x)+1
               fx2(x)=fx2(x)+1
            end if

            if (amp.eq.0.0) then
               amp1=amptab(2*x-1)
               amp2=amptab(2*x)
               do z=1,nzpc
                  wallvr(x,z)=amp1*fx1(x)*ft
                  wallvi(x,z)=amp2*fx2(x)*ft
               end do
            else
               do z=1,nzpc
                  wallvr(x,z)=amp*fx1(x)*ft
                  wallvi(x,z)=amp*fx2(x)*ft
               end do
            endif
         end do
      end if

      if (wbci.eq.4) then
c     
c     Spatial/temporal wall oscillation
c
         amp=wp1
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         tstart=wpds9
         tomeg=wpdds2
         xomeg=wpdds1
c
c     If steady, only impose it at first time step
c
         if (it.gt.1.and.tomeg.eq.0) return
         
         if (xend.ge.xl/2.0) then
            xend=xend-xl
            if (xstart.ge.xl/2.0) then
               xstart=xstart-xl
            end if     
         end if                  
         
         do x=1,nxp/2            
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            
            fx1(x)=step((xc11-xstart)/xrise)-
     &           step((xc11-xend)/xfall+1)
            
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            
            fx2(x)=step((xc22-xstart)/xrise)-
     &           step((xc22-xend)/xfall+1)
            
            if (xstart.gt.xend) then
               fx1(x)=fx1(x)+1
               fx2(x)=fx2(x)+1
            end if
            
            xstartc=xstart
            if (xstart.lt.0) then
               xstartc=xstart+xl
            endif
            if (x.le.nxp/4) then
               xstartc=xstartc-xl
            endif

            if (amp.eq.0.0) then
               amp1=amptab(2*x-1)
               amp2=amptab(2*x)
            else
               amp1 = amp
               amp2 = amp
            endif

            if (tstart.ne.0.) then
c
c     Extended spatial boundary condition
c
               do z=1,nzpc
                  wallvr(x,z)=amp1*fx1(x)*sin((xc11-xstartc)*xomeg)*
     &                 sin(tomeg*(tc-tstart))
                  wallvi(x,z)=amp2*fx2(x)*sin((xc22-xstartc)*xomeg)*
     &                 sin(tomeg*(tc-tstart))
               end do
            else
c
c     Normal spatial boundary condition
c
               do z=1,nzpc
                  wallvr(x,z)=amp1*fx1(x)*sin((xc11-xstartc)*xomeg)
                  wallvi(x,z)=amp2*fx2(x)*sin((xc22-xstartc)*xomeg)
               end do
            endif
         end do
      end if

      if (wbci.eq.5) then
c     
c     Travelling wall oscillation
c
         amp=wp1
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         xomeg=wpdds1
         tstart=wpds9
         tend=wpds10
         trise=wpds11
         tfall=wpds12
         tomeg=wpdds2
         
         if (trise.eq.0.0.and.tfall.eq.0.0) then
            ft=1.0
         else
            ft=step((tc-tstart)/trise)-step((tc-tend)/tfall+1.)
         end if

         if (xend.ge.xl/2.0) then
            xend=xend-xl
            if (xstart.ge.xl/2.0) then
               xstart=xstart-xl
            end if     
         end if                  
         
         do x=1,nxp/2            
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            fx1(x)=step((xc11-xstart)/xrise)-
     &           step((xc11-xend)/xfall+1)
            
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            fx2(x)=step((xc22-xstart)/xrise)-
     &           step((xc22-xend)/xfall+1)
            
            if (xstart.gt.xend) then
               fx1(x)=fx1(x)+1
               fx2(x)=fx2(x)+1
            end if
            
            xstartc=xstart
            if (xstart.lt.0) then
               xstartc=xstart+xl
            endif
            if (x.le.nxp/4) then
               xstartc=xstartc-xl
            endif

            if (amp.eq.0.0) then
               amp1 = ft*amptab(2*x-1)
               amp2 = ft*amptab(2*x)
            else
               amp1 = ft*amp
               amp2 = ft*amp
            endif

            do z=1,nzpc
               wallvr(x,z)=amp1*fx1(x)*
     &              sin( (xc11-xstartc)*xomeg-tomeg*(tc-tstart) )
               wallvi(x,z)=amp2*fx2(x)*
     &              sin( (xc22-xstartc)*xomeg-tomeg*(tc-tstart) )
            end do
         end do
      end if

c
c     Postprocess surface data
c
c     Real to half complex transform in x-direction first
c
      call vrfftf(wallvr,wallvi,wr,wi,nxp,nzpc,1,nxp/2+1,prex)
c
c     Then complex transform in z direction
c
      if (nfzsym.eq.0) then
         call vcfftf(wallvr,wallvi,wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      else
         call vcffts(wallvr,wallvi,wr,nzst,nx/2,nxp/2+1,1,pres)
      end if
c
c     Normalize and change to small grid in z
c
      do x=1,nxp/2+1
         do z=1,nz/2
            wallvr(x,z)=wallvr(x,z)/(nxp*nzp)
            wallvi(x,z)=wallvi(x,z)/(nxp*nzp)
         end do
      end do
      if (nfzsym.eq.0) then
         do x=1,nxp/2+1
            do z=nz/2+1,nz
               wallvr(x,z)=wallvr(x,z+nzp-nz)/real(nxp*nzp)
               wallvi(x,z)=wallvi(x,z+nzp-nz)/real(nxp*nzp)
            end do
         end do
      end if
c
c     Remove odd-ball modes? 
c     Not necessary since linearbl is not executed for 
c     the odd-ball modes.
c
      if (mod(nz,2).eq.0) then
         do x=1,nxp/2+1
            wallvr(x,nz/2+1)=0.0
            wallvi(x,nz/2+1)=0.0
         end do
      end if


      if (wbci.eq.3.or.wbci.eq.4.or.wbci.eq.5) then
c
c     If w is needed rather than v...
c
c     For u, use the commented parts
c     
c     Compute Dv and eta
c     
         do z=1,memnz
            zb = my_node*memnz+z
            do x=1,nx/2
               wallur(x,z) =  beta(zb)*wallvi(x,zb)
               wallui(x,z) = -beta(zb)*wallvr(x,zb)
               wallwr(x,z) =  alfa(x)*wallvi(x,zb)
               wallwi(x,z) = -alfa(x)*wallvr(x,zb)
c
c     For u use the following instead
c
c               wallur(x,z) =  alfa(x)*wallvi(x,zb)
c               wallui(x,z) = -alfa(x)*wallvr(x,zb)
c               wallwr(x,z) = -beta(zb)*wallvi(x,zb)
c               wallwi(x,z) =  beta(zb)*wallvr(x,zb)
            end do
         end do
c
c     In (1,1), i.e. (alfa,beta)=(0,0) the values of u and w
c     and not Dv and eta are put back (not necessary for v)
c
         if (my_node.eq.0) then
            wallur(1,1) = 0.
            wallui(1,1) = 0.
            wallwr(1,1) = wallvr(1,1)      
            wallwi(1,1) = 0.
c
c     For u use the following instead
c
c            wallur(1,1) = wallvr(1,1)
c            wallui(1,1) = 0.
c            wallwr(1,1) = 0.      
c            wallwi(1,1) = 0.
         end if
c
c     Remove temporary vr,vi
c
         do z=1,nz
            do x=1,nx/2
               wallvr(x,z) = 0.
               wallvi(x,z) = 0.
            end do
         end do
      end if

      end subroutine cwallbc
