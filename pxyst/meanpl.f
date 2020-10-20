c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine meanpl(utau,uinf,delta,deltar,beta,xl,zl,re,t,namxys,
     &     nux,thdelta,pr,scalarind,fltype)
c
c     Find local skin friction, displacement thickness and momentum
c     thickness, shape factor
c
      implicit none

      include 'par.f'

      real utau(nx),uinf(nx),delta(nx,7),deltar(nx),beta(nx)
      real nux(nx),thdelta(nx,2)
      real re,xl,zl,t,pr(scalar),st1(scalar),st2(scalar)
      character*80 namxys,xmat

      integer x,nplot,npoint(3),mxr,mbox,xpr,ihard
      integer ixvar,iyvar,i,j
      integer flag
      real px(nx+1,3),py(nx+1,3)
      character*20 xhead
      character*40 yhead
      character*24 xvar,yvar
      character*80 heada,headb,headc
      real xr,xmin,xmax,dx,pxmin,pxmax,cl,ct,rex
      integer scalarind,fltype

      st1(1)=.315
      st2(1)=.0175
      st1(2)=.332
      st2(2)=.0134
      st1(3)=.453
      st2(3)=.014
      st1(4)=.332
      st2(4)=.0095
      st1(5)=.453
      st2(5)=.0096
 1000 continue
      write(*,*) 'Variable on y-axis'
      write(*,*) '0 no more mean plots'
      write(*,*) '1 Retau'
      write(*,*) '2 uinf resp. ucl'
      write(*,*) '3 uinf/utau'
      write(*,*) '4 cf'
      write(*,*) '5 h12'
      write(*,*) '6 displth,momentlossth'
      write(*,*) '7 redisplth,remomentlossth'
      write(*,*) '8 redispl+,remomentloss+'
      write(*,*) '9 Reynolds thickness'
      write(*,*) '10 beta'
      write(*,*) '11 Delta'
      write(*,*) '12 Intermittency gamma (ReX)'
      if (scalar.gt.0) then
         write(*,*) '13 Stanton number'
         write(*,*) '14 2St/Cf'
         write(*,*) '15 thermal displth,momentlossth'
         write(*,*) '16 boundary layer thickness'
         write(*,*) '17 thickness ratios'
      end if
      write(*,*) '18 dtheta/dx'
      write(*,*) '19 int(urms2-vrms2)dy'
      write(*,*) '20 int(P)dy/int(eps_flut)dy'
      read(*,*) iyvar
      if (iyvar.eq.0) return
      write(*,*) 'Variable on x-axis'
      write(*,*) '1 x'
      write(*,*) '2 rex'
      write(*,*) '3 retheta'
      if (iyvar.eq.4.or.iyvar.eq.13.or.iyvar.eq.14) then
         write(*,*) '4 x with solution'
         write(*,*) '5 rex with solution'
      end if
      read(*,*) ixvar
      if (ixvar.gt.3) then
         ixvar = ixvar -3
         flag = 1
      else
         flag = 0
      end if
      write(*,9000) xl
 9000 format(' Give start and end for x ( 0-',F8.3,')')
      read(*,*) xmin,xmax
      write(*,*) 'type of output'
      write(*,*) '0 screen only'
      write(*,*) '1 postscript and screen'
      write(*,*) '2 postscript only'
      write(*,*) '3 tektronix file'
      write(*,*) '4 matlab file'
      read(*,*) ihard
      if (ihard.eq.4) then
      write(*,*) 'Give file to write to'
      read(*,'(a)') xmat
      end if
c
      if (xmin.eq.xmax) then
         xmin=0.
         xmax=xl
      end if
      xr=xl/2.+xmin
      dx=xl/nx
      mxr=int(xr/dx+.5)
      if (xr.lt.0.) mxr=-int((-xr)/dx+.5)
      xr=mxr*dx
      do x=1,nx+1
         xpr=mod(10*nx+x+mxr-1,nx)+1
         px(x,1)=real(x-nx/2-1)/real(nx)*xl+xr
         rex=px(x,1)*re+(re/1.73)**2
         if (fltype.eq.6) then
            if (iyvar.eq.1) py(x,1)=re*utau(xpr)*delta(xpr,1)
         else
            if (iyvar.eq.1) py(x,1)=re*utau(xpr)
         end if
         if (iyvar.eq.2) py(x,1)=uinf(xpr)
         if (iyvar.eq.3) py(x,1)=uinf(xpr)/utau(xpr)
         if (iyvar.eq.4) then
            py(x,1)=2.*utau(xpr)**2/uinf(xpr)**2
            py(x,2)=1.33/sqrt(rex)/2.
            py(x,3)=.0576*(rex)**(-.2)
         end if
         if (iyvar.eq.5) py(x,1)=delta(xpr,3)
         if (iyvar.eq.6) py(x,1)=delta(xpr,1)
         if (iyvar.eq.6) py(x,2)=delta(xpr,2)
         if (iyvar.eq.7) py(x,1)=delta(xpr,1)*re*uinf(xpr)
         if (iyvar.eq.7) py(x,2)=delta(xpr,2)*re*uinf(xpr)
         if (iyvar.eq.8) py(x,1)=delta(xpr,1)*re*utau(xpr)
         if (iyvar.eq.8) py(x,2)=delta(xpr,2)*re*utau(xpr)
         if (iyvar.eq.9) py(x,1)=deltar(xpr)
         if (iyvar.eq.10) py(x,1)=beta(xpr)
         if (iyvar.eq.11) py(x,1)=delta(xpr,4)
         if (iyvar.eq.12) then
            cl=1.33/sqrt(rex)/2.
            ct=.0576*(rex)**(-.2)
            py(x,1)=(cl-2.*utau(xpr)**2/uinf(xpr)**2)/(cl-ct)
         end if
         if (iyvar.eq.13) then
            py(x,1)=nux(xpr)
            py(x,2)=st1(scalarind)/(sqrt(rex)
     &           *pr(scalarind)**(2./3.))
            py(x,3)=st2(scalarind)*pr(scalarind)**(-1./10)
     &           *(.664*sqrt(rex))**(-1./4.)
         end if
         if (iyvar.eq.14) py(x,1)=nux(xpr)/(utau(xpr)**2/uinf(xpr)**2)
         if (iyvar.eq.14) py(x,2)=pr(scalarind)**(-2./3.)
         if (iyvar.eq.14) py(x,3)=1.18
         if (iyvar.eq.15) py(x,1)=thdelta(xpr,1)
         if (iyvar.eq.15) py(x,2)=thdelta(xpr,2)
         if (iyvar.eq.16) py(x,1)=4.92*sqrt(rex)/re
         if (iyvar.eq.16) py(x,2)=4.797*sqrt(rex)
     &        /(re*pr(scalarind)**(1./3.))
         if (iyvar.eq.17) py(x,1)=1.73*sqrt(rex)/re

         if (iyvar.eq.18) py(x,1) = delta(xpr,5)
         if (iyvar.eq.19) py(x,1) = delta(xpr,6)
         if (iyvar.eq.20) then
            py(x,1) = delta(xpr,7)
            py(x,2) = 1
         end if

         if (px(x,1).le.xmax+dx/2.) npoint(1)=x
         if (ixvar.eq.2) px(x,1)= rex
         if (ixvar.eq.3) px(x,1)=delta(xpr,2)*re*uinf(xpr)
         px(x,2)=px(x,1)
         px(x,3)=px(x,1)
      end do
      if (ixvar.eq.1) xhead='x'
      if (ixvar.eq.2) xhead='rex'
      if (ixvar.eq.3) xhead='re mom loss th'
      if (iyvar.eq.1) yhead='friction velocity'
      if (iyvar.eq.2) yhead='u at outer boundary'
      if (iyvar.eq.3) yhead='uinf/utau'
      if (iyvar.eq.4) yhead='skin friction 2*tau/(ro*uinf**2)'
      if (iyvar.eq.5) yhead='shape factor h12'
      if (iyvar.eq.6) yhead='displ, mom loss th'
      if (iyvar.eq.7) yhead='re displ,re mom loss th'
      if (iyvar.eq.8) yhead='retau displ,retau mom loss th'
      if (ixvar.eq.1) xvar='x'
      if (ixvar.eq.2) xvar='rex'
      if (ixvar.eq.3) xvar='re theta'
      if (iyvar.eq.1) yvar='retau'
      if (iyvar.eq.2) yvar='uinf'
      if (iyvar.eq.3) yvar='uinf/utau'
      if (iyvar.eq.4) yvar='cf'
      if (iyvar.eq.5) yvar='h12'
      if (iyvar.eq.6) yvar='dstar, theta'
      if (iyvar.eq.7) yvar='re dstar,re theta'
      if (iyvar.eq.8) yvar='redstar+,retheta+'
      if (iyvar.eq.9) yvar='Reynolds thickness'
      if (iyvar.eq.10) yvar='beta'
      if (iyvar.eq.11) yvar='Delta'
      if (iyvar.eq.13) yvar='Stanton'
      if (iyvar.eq.14) yvar='2St/Cf'
      if (iyvar.eq.15) yvar='thermal displ,mom loss th'
      if (iyvar.eq.16) yvar='boundary layer thickness'
      if (iyvar.eq.17) yvar='theta/theta_0'
      if (iyvar.eq.20) yhead='int(p)dy/int(eps_fluc)dy'
      nplot=1
      if (iyvar.ge.6.and.iyvar.le.8) nplot=2

      if (iyvar.eq.4.and.flag.eq.1) nplot=3
      if (iyvar.eq.13.and.flag.eq.1) nplot=3
      if (iyvar.eq.14.and.flag.eq.1) nplot=3
      if (iyvar.ge.15) nplot=2
      if (iyvar.ge.18) nplot=1
      if (iyvar.ge.20) nplot=2

      write(heada,9510) t,re
 9510   format(' t= ',F7.1,' re= ',F7.1,'$')
      write(headb,9020) xl,zl,nx,ny,nz,namxys
 9020 format(' xl = ',F7.2,' zl= ',F6.2,' ',I4,'x',I3,'x',I3,
     &     ' file ',A20,'$')
      headc=xhead//' vs '//yhead
      npoint(2)=npoint(1)
      npoint(3)=npoint(1)

      mbox=1
      pxmin=0.
      pxmax=0.
      if (ixvar.eq.1) then
         pxmin=xmin
         pxmax=xmax
      end if
      if (ihard.eq.4) then
         write(*,*) nx+1, nplot
         open(unit=24,file=xmat)
         do i=1,nx+1
            write(24,'(SP100ES26.16E3)') px(i,1),(py(i,j),j=1,nplot)
         end do
         ihard=0
      end if
      call rita1a(px,py,pxmin,pxmax,0.,0.,npoint,nplot,nx+1,
     &       xvar,yvar,heada,headb,headc,mbox,1,ihard)

      goto 1000

      end subroutine meanpl
