C ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program rit
c
c     Interactive/batch plotting program
c
c     Can be compiled and run in three modes depending on the flag imode
c     1 interactive
c     2 batch
c     3 input for batch
c     Its value is set below
c
      implicit none

      include 'par.f'
c
c     Mode flag see above
c
      integer imode,id
      parameter (imode=1)
      parameter (id=(4-imode)/2)
c
c     Main storage
c
      complex ur(memnx*id+1-id,memny*id+1-id,memnz*id+1-id,5)
      complex urx(nx/2)
c
c     Plot plane
c
      integer nxy,nxz,nyz
      parameter (nxy=(nx+2)*nyp*id+2-2*id)
      parameter (nxz=(nx+2)*(nz+2)*id+1-id)
      parameter (nyz=nyp*(nz+1)*id+1-id)
      real plane(nxy),grid(nxy,2)
      complex cplane(nxy/2)
      complex pxy(nx/2,nyp)
      real plxz(nxz),grxz(nxz,2),plyz(nyz),gryz(nyz,2)
      equivalence (plane,cplane,plxz,plyz),(grid,grxz,gryz)
      real w1(nxy),w2(nxy),w3(nxy),w4(nxy)
      real w1xz(nxz),w2xz(nxz),w3xz(nxz),w4xz(nxz)
      real w1yz(nyz),w2yz(nyz),w3yz(nyz),w4yz(nyz)
      equivalence (w1,w1xz,w1yz),(w2,w2xz,w2yz),(w3,w3xz,w3yz)
      equivalence (w4,w4xz,w4yz)
      real graf(nx+nyp+nz,3),grafy(nyp,20),grafz(nz,3)
      real gxax(nx+nyp+nz,3),gxaxy(nyp,20),gxaxz(nz,3)
      equivalence(graf,grafy,grafz),(gxax,gxaxy,gxaxz)

      character*1 ans
      character*80 namnin,xmat
c
c     Geometrics
c
      real eta(nyp),wint(nyp)
      real pi,xl,zl,dx,xs,xr,re,t,dstar,xrd
      complex egr(nx/2)
      parameter (pi = 3.1415926535897932385)
c
c     Flow
c
      integer fltype
      real prex(nx+15),prey(nyp*2+15),prez(nz*2+15),alfa(nx/2),beta(nz)
c
c     Plot parameters
c
      logical logax,logc,lpgmr,uwsub,lmat
      integer iplot,ivar,jvar,lvar
      integer ifilt,lfilt
      real sfilt,zfilt,lsfilt,lzfilt
      integer mxs,mys,mzs,mx,my,mz,mkx,mky,mkz,mxr
      real sym
      integer ix,ii,i,yi,xi
      real blev,elev,dlev,xp,yp,zp
      integer mbox,ncont,x,y,z,ni,iopt
      character*80 heada,headb
      character*32 xaxis
      character*7 var,filt
      integer ihard
      integer npoint(20),nplot
      real ymin,ymax,xmin,xmax,sfunc,phiuw
      integer inunit,utunit
c
c     New subtract mean flow variables
c
      character*80 namnin2
c      complex meanfr(memnx*id+1-id,memny*id+1-id,memnz*id+1-id,5)
      complex meanfr(1,1,1,5)
      complex mfrx(nx/2)
      integer readmean

      data lvar,lfilt,lsfilt,lzfilt/0,0,0.,0./
c
c     For 32-bit compiled executables with large size
c
      common /tot/ ur,meanfr

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         rit $Rev$   '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
c
c     Set input unit number and open control file
c
      inunit=5
      if (imode.eq.2) then
         inunit=21
         utunit=6
         open(unit=inunit,file='ritin',status='old')
      end if
      if (imode.eq.3) then
         inunit=5
         utunit=21
         open(unit=utunit,file='ritin',status='unknown')
      end if
c
c     Initialize x,y,z-transforms
c
      call vrffti(nx,prex,0)
      call vcosti(nyp,prey,0)
      call vcffti(nz,prez,0)

      write(*,*) 'Give file to plot from'
      read(inunit,9000) namnin
      if (imode.ge.2) write(utunit,9000) namnin
 9000 format(a80)
      if (imode.ne.3) then
        call rdiscp(ur,re,xl,zl,t,xs,dstar,fltype,
     &       namnin,3,w1,urx)
      end if
      do y=1,nyp
         eta(y)=cos(pi*real(y-1)/real(nyp-1))
      end do
c
c     Subtraction of meanflow field
c
c      write(*,*)
c      write(*,*) 'Subtract another field, e.g. meanflow?'
c     &     //' (1) for yes (0) for no :   '
c      read(*,*) readmean
      readmean = 0
      if (readmean.eq.1) then
         write(*,*) 'Give name of other flow file :   '
         read(*,*) namnin2
         call rdiscp(meanfr,re,xl,zl,0.0,xs,dstar,fltype,
     &        namnin2,3,w1,mfrx)
         call submean(ur,meanfr)
      end if
c
c     Compute integration weights
c
      do y=1,nyp
         wint(y)=-.5
         do ni=1,(nyp+1)/2-2
            wint(y) = wint(y) + cos(pi*real(2*ni*(y-1))/real(nyp-1))/
     &           real((2*ni)**2-1)
         end do
         wint(y) = wint(y) + .5*cos(pi*real(y-1))/real((nyp-1)**2-1)
         wint(y) = -4./real(nyp-1)*wint(y)
         if (y.eq.1.or.y.eq.nyp) wint(y)=wint(y)*.5
      end do
      if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6.and.fltype.le.9
     &     .or.abs(fltype).eq.20)
     &     then
         do y=1,nyp
            eta(y)=(1.+eta(y))/dstar
            wint(y)=wint(y)/dstar
         end do
      end if
      if (imode.ne.3) then
         do x=1,nx/2
            alfa(x)=2.*pi/xl*real(x-1)
         end do
         beta(1)=0.
         do z=2,nz/2+1
            beta(z)=2.*pi/zl*real(z-1)
            beta(nz+2-z)=-2.*pi/zl*real(z-1)
         end do
c
c     Zero the wall value...?
c
c         if (fltype.ne.2.and.fltype.ne.5) call zwall(ur,w1)



         write(*,9010)xl,zl,t,re
 9010    format(' xl= ',F8.2,' zl= ',f8.2,' t= ',F8.1,' re= ',F8.1)
         if (fltype.eq.3) write(*,9015) 2./dstar
 9015    format(' yl= ',F8.2)
         if (fltype.eq.-3) write(*,*) 
     &        'Asymptotic suction boundary layer'
         if (fltype.eq.-2) write(*,*)
     &        'Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.-1) write(*,*) 'Falkner-Skan boundary layer'
         if (fltype.eq.1) write(*,*) 'Poiseuille flow'
         if (fltype.eq.2) write(*,*) 'Couette flow'
         if (fltype.eq.3) write(*,*) 'Blasius boundary layer flow'
         if (fltype.eq.4) write(*,*) 'Spatial Poiseuille flow'
         if (fltype.eq.5) write(*,*) 'Spatial Couette flow'
         if (fltype.eq.6) write(*,*) 'Spatial Blasius boundary layer'
         if (fltype.eq.7) write(*,*)
     &        'Spatial Falkner-Skan boundary layer'
         if (fltype.eq.8)
     &        write(*,*) 'Spatial Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.9) write(*,*) 'Spatial parallel boundary layer'
         if (fltype.eq.-20) write(*,*) 
     &        'Buoyancy-driven boundary layer (temporal)'
         if (fltype.eq.20) write(*,*) 
     &        'Buoyancy-driven boundary layer (spatial)'

      end if
c
c     Choose type of plot
c
 1000 continue
      write(*,*) 'Which type of plot ?'
      write(*,*) '0 Stop'
      write(*,*) '1 xy physical space'
      write(*,*) '2 xz physical space'
      write(*,*) '3 yz physical space'
      write(*,*) '4 xy spectral space'
      write(*,*) '5 xz spectral space'
      write(*,*) '6 yz spectral space'
      write(*,*) '7 Resolution check'
      write(*,*) '8 f(x) for one y,z'
      write(*,*) '9 f(y) for one x,z'
      write(*,*) '10 f(z) for one x,y'
      write(*,*) '11 f(y) for (alfa,beta)'
      write(*,*) '12 f(y) for multiple x, one z'
      write(*,*) '13 f(y) for multiple z, one x'
      write(*,*) '14 Two point correlation, x-direction, for one y'
      write(*,*) '15 Two point correlation, z-direction, for one y'
      write(*,*) '16 f+(y) for one x,z'
      write(*,*) '17 abs(f+(y+)) (halfchannel) for one x,z'
      write(*,*) '18 Boundary layer thickness and shape factor'
      write(*,*) '19 f(z) for one x and averaged in y'
      read(inunit,*) iplot
      if (imode.ge.2) write(utunit,*) iplot
      if (iplot.le.0.or.iplot.ge.20) then
        if (imode.eq.2) close(unit=inunit)
        if (imode.eq.3) close(unit=utunit)
        stop
      end if
c
c     Set compulsory parameters
c
      write(*,*) 'Variable ?'
      write(*,*) '1 u'
      write(*,*) '2 v'
      write(*,*) '3 w'
      write(*,*) '4 omegax'
      write(*,*) '5 omegay'
      write(*,*) '6 omegaz'
      write(*,*) '7 u-umean'
      write(*,*) '8 u-ulam'
      write(*,*) '9 dudy'
      write(*,*) '10 Energy = ((u-umean)**2+v**2+w**2)/2'
      write(*,*) '11 Spectral energy (uhat**2+vhat**2+what**2)/2'
      write(*,*)
     &'12 energy (for res. check) (uthat**2+vthat**2+wthat**2)/2'
      write(*,*) '13 Reynolds stress -(u-umean)*v'
      write(*,*) '14 Total shear -(u-umean)*v+dudy*nu'
      write(*,*) '15 Turbulence production -(u-umean)*v*dU'
      write(*,*) '16 Enstrophy, <omega,omega>'
      write(*,*) '17 Modulus of vorticity, <omega,omega>^(1/2)'
      write(*,*) '18 (u,w) projected'
      write(*,*) '19 Continuity, dui/dxi'
      read(inunit,*) ivar
      if (imode.ge.2) write(utunit,*) ivar
c
c     Default options
c
      if (imode.ne.3) then
         write(headb,9020) xl,zl,nx,ny,nz,namnin
 9020    format(' xl = ',F7.2,' zl= ',F6.2,' ',I4,'x',I3,'x',I3,
     &        ' file ',A20,'$')
      end if
      mxs=1
      mys=1
      mzs=1
      mkx=nx/2
      mky=ny
      logax=.false.
      logc=.true.
      lpgmr=.false.
      ifilt=0
c
c     mbox=1 max plotarea, mbox=0 grid selected proportions
c
      mbox=1
      if (iplot.eq.2) mbox=0
      iopt=0
      ihard=-1
c
c     Shifting to follow disturbance in Poiseuille flow
c     xr is the coordinate of the center of the plot
c     xs is the coordinate in the center of the computational box
c
      if (imode.ne.3) then
        xr=0.
c
c     Shifting for channel flow disabled...
c
        if (fltype.eq.1) xr=xl/2.
c        if (fltype.eq.1) xr=2./3.*t
        if (fltype.eq.3.or.fltype.lt.0) xr=2./3.*t
        if (fltype.eq.3.or.fltype.lt.0) xr=0.
        if (fltype.ge.4) xr=xl/2.
c
c     Special for forced rotating channel
c
c        xr=xl/2-3.

        dx=xl/nx
        mxr=(xr-xs)/dx+.5
        if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
        xr=xs+mxr*dx
      end if
c
c     Get variable dependent parameters
c
      if (ivar.eq.18) then
         write(*,*) 'Give angle of projection (degrees)'
         write(*,*) 'Relative x-axis'
         read(inunit,*) phiuw
         if (imode.ge.2) write(utunit,*) phiuw
         phiuw=phiuw/180.*pi
         uwsub=.false.
         if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.4.or.fltype.eq.5)
     &        then
            write(*,*) 'Subtract laminar profile ? (y/n)'
            read(inunit,9900) ans
            if (imode.ge.2) write(utunit,9900) ans
            uwsub=ans.ne.'N'.and.ans.ne.'n'
         end if
      end if
c
c     Get plot dependent parameters
c
      if (iplot.le.6) then
         write(*,*) 'Number of contours ? (0 to select manually)'
         read(inunit,*) ncont
         if (imode.ge.2) write(utunit,*) ncont
         ncont=-ncont
         if (ncont.eq.0) then
            write(*,*)  'give start,end,spacing for contour levels'
            read(inunit,*) blev,elev,dlev
            if (imode.ge.2) write(utunit,*) blev,elev,dlev
            ncont=(elev-blev)/dlev+1.5
            elev=blev+dlev*real(ncont-1)
            write(*,*) 'endlevel ',elev
         end if
      end if

      if (iplot.eq.11) then
         if (imode.eq.3) then
            write(*,*) 'what alfa and beta ?'
         else
            write(*,9031)
 9031       format(' what alfa and beta ? ')
         end if
         read(inunit,*) x,z
         write(*,*)'letti ',x,z

         if (imode.ge.2) write(utunit,*) xp
      end if

      if (iplot.eq.3.or.iplot.eq.9.or.iplot.eq.10.or.iplot.eq.16
     &     .or.iplot.eq.17.or.iplot.eq.19) then
         if (imode.eq.3) then
            write(*,*) 'What x-value ?'
         else
            write(*,9030) xs-xl/2,xs+xl/2
 9030       format(' What x-value ? (primary range ',F8.2,' -',F8.2,')')
         end if
         read(inunit,*) xp
         if (imode.ge.2) write(utunit,*) xp
         if (imode.ne.3) then
            if (xp-xs.lt.0) then
               x=int((xp-xs)/xl*real(nx)-.5)+nx/2+1
            else
               x=int((xp-xs)/xl*real(nx)+.5)+nx/2+1
            end if

            xp=xl/real(nx)*real(x-nx/2-1)+xs
            write(*,9040) xp,x
 9040       format(' Nearest plane x=',F6.3,' plane number',I4)
         end if
      end if
c
      if (iplot.eq.12) then
         if (imode.eq.3) then
            write(*,*) 'What first x-value ?'
         else
            write(*,9035) xs-xl/2,xs+xl/2
 9035       format(' What first x-value ?',
     &           ' (primary range ',F8.2,' -',F8.2,')')
         end if
         read(inunit,*) xp
         if (imode.ge.2) write(utunit,*) xp
         if (imode.ne.3) then
            if (xp-xs.lt.0) then
               x=int((xp-xs)/xl*real(nx)-.5)+nx/2+1
            else
               x=int((xp-xs)/xl*real(nx)+.5)+nx/2+1
            end if
            xp=xl/real(nx)*real(x-nx/2-1)+xs
            write(*,9045) xp,x
 9045       format(' Nearest plane x=',F8.3,' plane number',I4)
         end if
 1010    if (imode.ne.3) then
            write(*,9048) xl/real(nx)
 9048       format(' One grid step is ',F6.3,'.')
         end if
         write(*,*) ' How many steps per plot ?'
         read(inunit,*) ix
         if (imode.ne.3) then
            if (nx/ix.gt.20) then
               write(*,*) 'More than 20 plots not allowed'
               goto 1010
            end if
         end if
         if (imode.ge.2) write(utunit,*) ix
      end if
c
      if (iplot.eq.2.or.iplot.eq.8.or.iplot.eq.10.or.iplot.eq.14
     &     .or.iplot.eq.15) then
         if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.4.or.fltype.eq.5)
     &        then
            write(*,*)  ' What y-value ? (-1 - 1)'
         else
            write(*,*)  ' What y-value ? (0 - ',.1*int(20./dstar+.5)
         end if
         read(inunit,*) yp
         if (imode.ge.2) write(utunit,*) yp
         if (imode.ne.3) then
            if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.4.or.
     &           fltype.eq.5) then
               y=int(acos(yp)/pi*real(nyp-1)+1.5)
            else
               y=int(acos(yp*dstar-1.)/pi*real(nyp-1)+1.5)
            end if
            yp=eta(y)
            write(*,9055) yp,y
 9055       format(' Nearest plane y=',F6.3,' plane number',I4)
         end if
      end if
c
      if (iplot.eq.1.or.iplot.eq.8.or.iplot.eq.9.or.iplot.eq.12.or.
     &     iplot.eq.16.or.iplot.eq.17.or.iplot.eq.18) then
         if (imode.eq.3) then
            write(*,*) ' What z-value ?'
         else
            write(*,9060) -zl/2,zl/2
 9060       format(' What z-value ? (',F8.2,' -',F8.2,')')
         end if
         read(inunit,*) zp
         if (imode.ge.2) write(utunit,*) zp
         if (imode.ne.3) then
            z=int(zp/zl*real(nz)+.5+real(nz/2+1))
            zp=zl/real(nz)*real(z-nz/2-1)
            write(*,9070) zp,z
 9070       format(' Nearest plane z=',F8.2,' plane number',I4)
         end if
      end if
c
c     Set optional parameters
c
      write(*,*) 'Change default options ? (y/n)'
      read(inunit,9900) ans
      if (imode.ge.2) write(utunit,9900) ans
 9900 format(a1)
      lmat = .false.
      if (ans.ne.'N'.and.ans.ne.'n') then
         iopt=1
         write(*,*) 'Type of output'
         write(*,*) '0 screen only'
         write(*,*) '1 postscript and screen'
         write(*,*) '2 postscript only'
         write(*,*) '3 tektronix file'
         if (iplot.eq.7.or.iplot.eq.8
     &        .or.iplot.eq.9.or.iplot.eq.10
     &        .or.iplot.eq.15) write(*,*)
     &        '4 matlab file and screen'
         if (iplot.le.6) write(*,*) '5 screen and pgmr file'
         if (iplot.le.6) write(*,*) '8 tek file and pgmr file'
         read(inunit,*) ihard
         if (ihard.eq.4) then
            write(*,*) 'Give file to write to'
            read(inunit,9000) xmat
            if (imode.ge.2) write(utunit,9000) xmat
            lmat=.true.
            ihard=ihard-4
         end if
         if (imode.ge.2) write(utunit,*) ihard
         if (ihard.ge.4) then
            lpgmr=.true.
            ihard=ihard-5
         end if
         if (iplot.eq.5) then
            write(*,*)  'Logarithmic axes ? (y/n)'
            read(inunit,9900) ans
            if (imode.ge.2) write(utunit,9900) ans
            logax=ans.ne.'N'.and.ans.ne.'n'
            write(*,*)  'Logarithmic contours ? (y/n)'
            read(inunit,9900) ans
            if (imode.ge.2) write(utunit,9900) ans
            logc=ans.ne.'N'.and.ans.ne.'n'
         end if
c
c     Filter
c
         write(*,*) 'Choose filter'
         write(*,*) '0 no filter'
         write(*,*) '1 gaussian lowpass filter'
         write(*,*) '2 gaussian highpass filter'
         write(*,*) '3 rms-value through highpass filter, squaring'
         write(*,*) '  lowpass filtering and taking the square root'
         write(*,*) '4 gaussian anisotropic lowpass filter'
         write(*,*) '5 gaussian anisotropic highpass filter'
         write(*,*) '6 rms-value through anisotropic filter'
         read(inunit,*) ifilt
         if (imode.ge.2) write(utunit,*) ifilt
         if (ifilt.ge.1) then
            if (ifilt.le.3) write(*,*) 'give filter lengthscale l'
            if (ifilt.ge.4) write(*,*) 'give filter lengthscale lx'
            if (ifilt.eq.2.or.ifilt.eq.5) then
               write(*,*) 'the filter is 1-exp(-(k*l/2pi)**2)'
            else
               write(*,*) 'the filter is exp(-(k*l/2pi)**2)'
            end if
            read(inunit,*) sfilt
            if (imode.ge.2) write(utunit,*) sfilt
         end if
         if (ifilt.ge.4.and.ifilt.le.6) then
            write(*,*) 'give filter lengthscale lz'
            read(inunit,*) zfilt
            if (imode.ge.2) write(utunit,*) zfilt
         end if
c
c     Aspect ratio
c
         write(*,*) 'Choose aspect ratio'
         write(*,*) '0 equal scale in both directions '
         write(*,*) '1 largest possible picture (aspect ratio 1.6)'
         read(inunit,*) mbox
         if (imode.ge.2) write(utunit,*) mbox
c
c     Boxshift
c
         if (iplot.eq.1.or.iplot.eq.2.or.iplot.eq.8.or.iplot.eq.18) then
            write(*,*) 'The center of the plot is now at x=',xr
            write(*,*) 'Give offset (0 for no change)'
            read(*,*) xrd
            xr=xr+xrd
            mxr=(xr-xs)/dx+.5
            if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
            xr=xs+mxr*dx
         end if
      end if
c
c     Set computed parameters
c
      sym=1.
      if (ivar.ge.3.and.ivar.le.5) sym=-1.
      mx=nx/mxs+1
      my=(nyp-1)/mys+1
      mz=nz/mzs+1
      mkx=nx/2
      if (logax) mkx=nx/2-1
      mkz=nz-1
      if (logax) mkz=nz/2-1
      if (ivar.eq.1) var='   u   '
      if (ivar.eq.2) var='   v   '
      if (ivar.eq.3) var='   w   '
      if (ivar.eq.4) var=' omegax'
      if (ivar.eq.5) var=' omegay'
      if (ivar.eq.6) var=' omegaz'
      if (ivar.eq.7) var='u-umean'
      if (ivar.eq.8) var=' u-ulam'
      if (ivar.eq.9) var='  dudy '
      if (ivar.ge.10.and.ivar.le.12) var=' energy'
      if (ivar.eq.13) var='  -uv  '
      if (ivar.eq.14) var='tot sh '
      if (ivar.eq.15) var='t prod '
      if (ivar.eq.16) var='enstr  '
      if (ivar.eq.17) var='mod vor'
      if (ivar.eq.18) var='c*u+s*w'
      if (ivar.eq.19) var='dui/dxi'

      if (ifilt.eq.0) filt=' '
      if (ifilt.eq.1) filt='lo-pass'
      if (ifilt.eq.2) filt='hi-pass'
      if (ifilt.eq.3) filt='rms'
      if (ifilt.eq.4) filt='lo-pass'
      if (ifilt.eq.5) filt='hi-pass'
      if (ifilt.eq.6) filt='rms'
c
c     Create plotdata
c     Compute plot variable
c
      if (ivar.gt.19) then
         write(*,*) 'This variable is not yet implemented'
         goto 1000
      end if
      if (ifilt.gt.6) then
         write(*,*) 'The filter is not yet implemented'
         goto 1000
      end if
      jvar=min(ivar,4)
      if ((lvar.eq.ivar.and.ifilt.eq.lfilt.and.(ifilt.eq.0.or.
     &     (ifilt.le.3.and.lsfilt.eq.sfilt).or.
     &     (lsfilt.eq.sfilt.and.lzfilt.eq.zfilt))).and.ivar.ne.18) then
c
c     Reuse old variable 4
c
         jvar=4
      else
         if (ivar.gt.3) then
            if (imode.ne.3) call cvar(ur,ivar,
     &           plane,w1,w2,w2,w3,w4,w4,
     &           prey,prex,prez,re,eta,alfa,beta,dstar,fltype,
     &           phiuw,uwsub)
            lvar=ivar
            lfilt=0
         end if
c
c     Filter
c
         if (ifilt.ne.0) then
            if (imode.ne.3) call gauss(ur,jvar,sym,ifilt,sfilt,zfilt,
     &           alfa,beta,prex,prez,w1,plane,cplane,w2,w3)
            lfilt=ifilt
            lsfilt=sfilt
            lzfilt=zfilt
            lvar=ivar
            jvar=4
         end if
      end if
c
c     Get plane to plot
c
      if (imode.ne.3) then
         if (iplot.eq.1) call xyplane(plane,mx,my,mxs,mys,grid,z,jvar,
     &        sym,xr,xl,mxr,w1,w2,ur,prex,prez,eta)
         if (iplot.eq.2) call xzplane(plane,mx,mz,mxs,mzs,grid,y,jvar,
     &        sym,xr,xl,zl,mxr,w1,ur,prex,prez)
         if (iplot.eq.3) call yzplane(plane,my,mz,mys,mzs,grid,x,jvar,
     &        sym,zl,w1,w2,ur,prex,prez,eta)
         if (iplot.lt.1.or.iplot.gt.19.or.iplot.eq.4.or.iplot.eq.6.or.
     &        iplot.eq.31.or.iplot.eq.13) then
            write(*,*) 'Not yet implemented'
            goto 1000
         end if
c     if (iplot.eq.4) call xyspec(plane,mkx,mky,grid,ivar,ur)
         if (iplot.eq.5) call xzspec(plane,mkx,mkz,grid,jvar,logax,logc,
     &        sym,wint,xl,zl,w1,ur)
c     if (iplot.eq.6) call yzspec(plane,mky,mkz,grid,ivar,ur)
         if (iplot.eq.7)  call reschk(graf,gxax,xmin,xmax,ymin,ymax,
     &        npoint,prey,jvar,w1,w2,ur)
         if (iplot.eq.8) call xline(graf,gxax,xmin,xmax,ymin,ymax,y,z,
     &        jvar,sym,mx,mxs,xl,mxr,xr,prex,prez,plane,w1,ur)
         if (iplot.eq.9.or.iplot.eq.16.or.iplot.eq.17)
     &        call yline(graf,gxax,xmin,xmax,ymin,ymax,x,z,
     &        jvar,sym,my,mys,eta,prex,prez,w1,w2,ur)
         if (iplot.eq.11)
     &        call yline_spec(grafy,gxaxy,xmin,xmax,ymin,ymax,npoint,
     &        x,z,jvar,eta,ur)
         if (iplot.eq.10) call zline(graf,gxax,imode,iopt,xmin,xmax,
     &        ymin,ymax,x,y,jvar,sym,mz,mzs,zl,prex,prez,plane,w1,ur,
     &        inunit,utunit)
         if (iplot.eq.12) call myline(graf,gxax,npoint,nplot,sfunc,
     &        xmin,xmax,ymin,ymax,x,z,ix,xp,jvar,sym,my,mys,
     &        xl,eta,prex,prez,w1,w2,ur)
         if (iplot.eq.14) call twox(graf,gxax,xmin,xmax,ymin,ymax,
     &        y,jvar,sym,xl,w1,w2,ur,prex,prez)
         if (iplot.eq.15) call twoz(graf,gxax,xmin,xmax,ymin,ymax,
     &        y,jvar,sym,zl,w1,w2,ur,prex,prez)
         if (iplot.eq.16) call scalef(graf,my,ymin,ymax,ivar,
     &        w1,w2,prey,re,ur)
         if (iplot.eq.17) call scaley(graf,gxax,my,xmin,xmax,ymin,ymax,
     &        ivar,w1,w2,prey,re,ur)
         if (iplot.eq.18) call blth(graf,gxax,xmin,xmax,ymin,ymax,z,
     &        jvar,sym,dstar,xr,xl,mxr,npoint,w1,w2,w3,w4,ur,
     &        prex,prez,prey)
         if (iplot.eq.19) then
            do z=1,nzc
               call getxyp(pxy,z,jvar,ur)
               do xi=1,nx/2
                  egr(xi) = 0.
                  do yi=1,nyp
                     egr(xi) = egr(xi) + pxy(xi,yi)*wint(yi)
                  end do
               end do
               do yi=1,nyp
                  do xi=1,nx/2
                     pxy(xi,yi) = egr(xi)*dstar/2.
                  end do
               end do
               call putxyp(pxy,z,4,ur)
            end do
            jvar=4
            y=nyp    ! Any y-position will do
            yp=-99.  ! will be used in the figure
            call zline(graf,gxax,imode,iopt,xmin,xmax,
     &        ymin,ymax,x,y,jvar,sym,mz,mzs,zl,prex,prez,plane,w1,ur,
     &        inunit,utunit)
         end if
      end if
c
c     Create plot
c     Note that for imode=3 cont5 should be called to enable
c     answering questions of subarea selection
c     but no other plot routines
c
      if (iplot.eq.1) then
         if (imode.ne.3) write(heada,9510) var,filt,zp,t,re
 9510    format(1x,A7,1X,A7,', z= ',F6.2,' t= ',F7.1,' re= ',F7.1,'$')
         call cont5(mx,my,plane,grid,iopt,
     &        ncont,blev,elev,'x$','y$',heada,headb,mbox,ihard,
     &        imode,inunit,utunit,lpgmr)
      end if
      if (iplot.eq.2) then
         if (imode.ne.3) write(heada,9520) var,filt,yp,t,re
 9520    format(1x,A7,1x,A7,', y= ',F6.3,' t= ',F7.1,' re= ',F7.1,'$')
         call cont5(mx,mz,plane,grid,iopt,
     &        ncont,blev,elev,'x$','z$',heada,headb,mbox,ihard,
     &        imode,inunit,utunit,lpgmr)
      end if
      if (iplot.eq.3) then
         if (imode.ne.3) write(heada,9530) var,filt,xp,t,re
 9530    format(1x,A7,1X,A7,', x= ',F6.2,' t= ',F7.1,' re= ',F7.1,'$')
         call cont5(mz,my,plane,grid,iopt,
     &        ncont,blev,elev,'z$','y$',heada,headb,mbox,ihard,
     &        imode,inunit,utunit,lpgmr)
      end if
      if (iplot.eq.4) call cont5(mkx,mky,plane,grid,iopt,
     &     ncont,blev,elev,'kx$','ky$',heada,headb,mbox,ihard,
     &     imode,inunit,utunit,lpgmr)
      if (iplot.eq.5) then
         if (imode.ne.3) then
            if (logc) then
               write(heada,9550) var,filt,t,re
 9550          format(' -log ',A7,1X,A7,', integr in y, t= ',
     &              F7.1,' re= ',F7.1,'$')
            else
               write(heada,9555) var,filt,t,re
 9555          format(1x,A7,1X,A7,', integr in y, t= ',F7.1,
     &              ' re= ',F7.1,'$')
            end if
         end if
         call cont5(mkx,mkz,plane,grid,iopt,
     &        ncont,blev,elev,'kx$','kz$',heada,headb,mbox,ihard,
     &        imode,inunit,utunit,lpgmr)
      end if
      if (iplot.eq.6) call cont5(mkz,mky,plane,grid,iopt,
     &     ncont,blev,elev,'kz$','ky$',heada,headb,mbox,ihard,
     &     imode,inunit,utunit,lpgmr)
      if (imode.ne.3) then
         if (iplot.eq.7) then
            write(heada,9570) var,filt,t,re
 9570       format(' log ',A7,1X,A7,', res.check, t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,npoint,3,
     &           nx+nyp+nz,'log n$','log e$',heada,headb,' $',mbox,1,
     &           ihard)
         end if

         if (iplot.eq.11) then
            write(heada,9591) var,filt,x,z,t,re
 9591       format(A7,1X,A7,' alfa= ',I3,' beta= ',I3,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxaxy,grafy,xmin,xmax,ymin,ymax,npoint,3,nyp,
     &           'y$','velocity$',heada,headb,' $',mbox,1,ihard)
         end if

         if (iplot.eq.8) then
            write(heada,9580) var,filt,yp,zp,t,re
 9580       format(A7,1X,A7,' y= ',F6.2,' z= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            if (lmat) then
               open(unit=24,file=xmat)
               do i=1,1
                  do yi=1,nx+1
                     write(24,*)gxax(yi,i),graf(yi,i)
                  end do
               end do
            end if
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,mx,1,nx+nyp+nz,
     &           'x$',var//'$',heada,headb,' $',mbox,1,ihard)
         end if
         if (iplot.eq.9) then
            write(heada,9590) var,filt,xp,zp,t,re
 9590       format(A7,1X,A7,' x= ',F6.2,' z= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            if (lmat) then
               open(unit=24,file=xmat)
               do i=1,1
                  do yi=1,nyp
                     write(24,*)gxax(yi,i),graf(yi,i)
                  end do
               end do
            end if
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,my,1,nx+nyp+nz,
     &           'y$',var//'$',heada,headb,' $',mbox,1,ihard)
         end if
         if (iplot.eq.10.or.iplot.eq.19) then
            write(heada,9600) var,filt,xp,yp,t,re
 9600       format(A7,1X,A7,' x= ',F6.2,' y= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            if (lmat) then
               open(unit=24,file=xmat)
               do i=1,1
                  do yi=1,nz
                     write(24,*)gxax(yi,i),graf(yi,i)
                  end do
               end do
            end if
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,mz,1,nx+nyp+nz,
     &           'z$',var//'$',heada,headb,' $',mbox,1,ihard)
         end if
c
         if (iplot.eq.12) then
            write(heada,9620) var,filt,zp,t,re
 9620       format(A7,1X,A7,' z= ',F6.2,' t= ',F7.1,' re= ',F7.1,'$')
            write(xaxis,9625) sfunc,var
 9625       format('x+',F7.1,'*',A7,'$')
            call rita1a(graf,gxax,xmin,xmax,ymin,ymax,npoint,nplot,my,
     &           xaxis,'y$',heada,headb,' $',mbox,0,ihard)
         end if

         if (iplot.eq.14) then
            write(heada,9640) var,filt,yp,t,re
 9640       format(A7,1X,A7,' y= ',F6.3,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,nx/2+1,1,nx/2+1,
     &           'deltax$','corr$',heada,headb,' $',mbox,1,ihard)
         end if
         if (iplot.eq.15) then
            write(heada,9650) var,filt,yp,t,re
 9650       format(A7,1X,A7,' y= ',F6.3,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,nz/2+1,1,nz/2+1,
     &           'deltaz$','corr$',heada,headb,' $',mbox,1,ihard)

            if (lmat) then
               open(unit=24,file=xmat)
               do   i=1,1
                  do  yi=1,nz/2+1
                     write(24,*)gxax(yi,i),graf(yi,i)
                  end do
               end do
            end if
         end if
         if (iplot.eq.16) then
            write(heada,9660) var,filt,xp,zp,t,re
 9660       format(A7,1X,A7,' x= ',F6.2,' z= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,my,1,nx+nyp+nz,
     &           'y$',var//'+$',heada,headb,' $',mbox,1,ihard)
         end if
         if (iplot.eq.17) then
            write(heada,9670) var,filt,xp,zp,t,re
 9670       format(A7,1X,A7,' x= ',F6.2,' z= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,(my+1)/2,1,
     &           nx+nyp+nz,'y+$',var//'+$',heada,headb,' $',mbox,1,
     &           ihard)
         end if
         if (iplot.eq.18) then
            write(*,*) 'Give end for x'
            read(*,*) xmax
            do 1233 ii=1,npoint(1)
               if (gxax(ii,1).gt.xmax) goto 1234
 1233       continue
            ii=npoint(1)+1
 1234       continue
            npoint(1)=ii-1
            npoint(2)=ii-1
            npoint(3)=ii-1
            xmin=gxax(1,1)
            ymin=0.0
            ymax=3.0
            write(heada,9680) var,filt,yp,zp,t,re
 9680       format(A7,1X,A7,' y= ',F6.2,' z= ',F6.2,' t= ',F7.1,
     &           ' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,npoint,3,
     &           nx+nyp+nz,'x$','d1/d2/h12$',heada,headb,' $',mbox,1,
     &           ihard)
         end if
      end if
c
c     Next plot
c
      goto 1000

      end program rit
