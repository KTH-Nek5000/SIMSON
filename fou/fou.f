c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program fou
c
c     Fourier transform in time from multiple fields
c
      implicit none

      include 'par.f'
c
c     Max number of transform frequencies in time (omega:s) mfr
c
      integer mfr
      parameter (mfr=2)
c
c     Max number of files
c
      integer mfile
      parameter (mfile=2)
c
c     Main storage
c
      complex ur(memnx,memny,memnz,3)
      complex urx(nx/2)
      complex spec(nx,nyp,nzc+2*nfzsym,3,mfr)
      real rspec(nx*2,nyp,nzc+2*nfzsym,3,mfr)
      equivalence (spec,rspec)
      real ws(nx*2*nyp*(nzc+2*nfzsym))
c
c     Geometrics
c
      real eta(nyp),wint(nyp),w1(nx+2,nyp)
      real pi,xl,zl,xs,xr,dx,re,t,dstar,sym
      parameter (pi = 3.1415926535897932385)
c
c     Flow
c
      integer fltype
      real prex(nx+15),prez(nz*2+15),beta(nz)
      real pres(nz+2+15),prea(nz*3/4+15),norm
c
c     Planes
c
      real pxz(nx+2,nz),w(nx+2,nz)
c
c     Temporal transform
c
      real fper,tf(mfile),dt(mfile)
      integer nfr
c
c     Filenames, etc
c
      character*80 namefile,fname(mfile)
      integer nfile
      logical header
c
c     Loop indices etc
c
      integer x,y,z,i,j,k,ni
c
c     Plotting
c
      integer ikfkzm,iloga
      parameter (ikfkzm=20)
      real efu(0:nx,ikfkzm),xx(0:nx,ikfkzm),urms(nx),vrms(nx)
      real wrms(nx),e(nx,ikfkzm),plane(0:nx,nyp),grid(0:nx,nyp,2)
      real xmin,xmax,ymin,ymax,eps,blev,elev,dlev
      integer kfa(ikfkzm),kza(ikfkzm),npoint(ikfkzm),mbox,ivar
      integer nplot,kf,kz,idash,ilog,ikfkz,iii,xpr,mxr,ncont
      character*80 heada,headb,headc,xaxis,yaxis
      character*7 mode(ikfkzm)
      logical logp
c
c     Additional functionality (fouLX)
c
      real giu,xpro,modu,maxmd
      integer xprof,ll,lll,zp,nex
      real yy(nyp,3),umod(nyp,3),maxmod(0:nx,3)
      real planex(0:nz,nyp),gridx(0:nz,nyp,2)
      character ans
      real maxps(0:nx,3),yex(nx+nyp),xex(nx+nyp)
      real y2(nx+nyp),u(nyp)
      real maxamp,minamp,psmax(3),psmin(3)
      real maxampy,minampy,ampy,psmaxy(3),lowy,highy
      real lowz0,highz0,lowzp,highzp,psmaxz(4)
      complex ui
      integer ihard
c      real fase

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         fou $Rev$   '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
c
c     Initialize x,y,z-transforms
c
      ui = (0.,1.)
      call vrffti(nx,prex,0)
      call vcffti(nz,prez,0)
      if (nz.gt.2) then
         call vcosti(nz/2+1,pres,0)
         call vsinti(nz/2-1,prea,0)
      end if
c
c     Zero the accumulation storage for time fourier transforms
c
      do x=1,nx*nyp*(nzc+2*nfzsym)*3*mfr
         spec(x,1,1,1,1)=(0.,0.)
      end do

c      write(*,*) 'Give file with filenames?'
c      read(*,*) namefile

      namefile='list.dat'
      open(unit=13,file=namefile,status='old')
      read(13,*) nfile
      if (nfile.gt.mfile) then
         write(*,*) 'Number of files larger than ',mfile
         write(*,*) 'Increase mfile parameter in fou.f'
         stop
      end if
c
c     Loop for each field for times and names
c
      do i=1,nfile
         header=.true.

         read(13,*) fname(i)

         call rdiscr(ur,re,xl,zl,tf(i),xs,dstar,fltype,
     &        fname(i),3,w1,urx,header)

         if (xs.ne.0.) then
            write(*,*) 'Wall velocity non-zero, xs = ',xs
            write(*,*) tf(i),' ',fname(i)
            stop
         end if
      end do
      close (13)

      if (nfile.eq.1) then
         fper=1.
         dt(1)=1.
      else
         fper=(tf(nfile)-tf(1))/real(nfile-1)*real(nfile)
         do i=2,nfile-1
            dt(i)=(tf(i+1)-tf(i-1))*0.5
         end do
         dt(1)=(tf(2)+fper-tf(nfile))*.5
         dt(nfile)=(tf(1)+fper-tf(nfile-1))*.5
      end if

      write(*,*) 'Tper= ',fper
      write(*,*) '     ',1,tf(1)
      do i=2, nfile
         write(*,*) '     ',i,tf(i),dt(i)
      end do
      write(*,*)
      write(*,*) 'Give number of frequences to calculate'
      write(*,*) 'max ',max(min(mfr,nfile/2),1)
      read(*,*) nfr
c
c     Loop for each field
c
      do i=1,nfile
         header=.false.
         call rdiscr(ur,re,xl,zl,t,xs,dstar,fltype,
     &        fname(i),3,w1,urx,header)
         do y=1,nyp
            eta(y)=cos(pi*real(y-1)/real(nyp-1))
         end do
         do y=1,nyp
            wint(y)=-.5
            do ni=1,(nyp+1)/2-2
               wint(y) = wint(y) + cos(pi*real(2*ni*(y-1))/
     &              real(nyp-1))/real((2*ni)**2-1)
            end do
            wint(y) = wint(y) +.5*cos(pi*real(y-1))/real((nyp-1)**2-1)
            wint(y) = -4./real(nyp-1)*wint(y)
            if (y.eq.1.or.y.eq.nyp) wint(y)=wint(y)*.5
         end do
         if (fltype.gt.5) then
            do y=1,nyp
               eta(y)=(1.+eta(y))/dstar
               wint(y)=wint(y)/dstar
            end do
         end if

         beta(1)=0.
         do z=2,nz/2+1
            beta(z)=2.*pi/zl*real(z-1)
            beta(nz+2-z)=-2.*pi/zl*real(z-1)
         end do
c
c     zwall not necessary if xs=0.0
c
c     if (fltype.eq.4.or.fltype.gt.5) call zwall(ur,w1,udr)
c
         write(*,9010)xl,zl,t,re
 9010    format(' xl= ',F8.2,' zl= ',f8.2,' t= ',F8.1,' re= ',F8.1)
         if (fltype.gt.5) write(*,9015) 2./dstar
 9015    format(' yl= ',F8.2)
         if (fltype.eq.4) write(*,*) 'Poiseuille flow'
         if (fltype.eq.5) write(*,*) 'Couette flow'
         if (fltype.gt.5) write(*,*) 'Boundary layer flow'
c
c     Loop over u,v,w velocites and the y-coordinate when accumulating
c     the frequency-wavenumber spectrum
c
         do k=1,3
            sym=1.
            if (k.eq.3) sym=-1.
            do y=1,nyp
               call getxzp(pxz,y,k,ur,sym)
c
c     Fourier transform to physical space
c
               call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
               call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
               if (fltype.eq.4) then
                  do x=1,nx
                     do z=1,nzc+nfzsym
                        if (k.eq.1) pxz(x,z)=pxz(x,z)-(1.-eta(y)**2)
                     end do
                  end do
               end if

               if (nz.le.2) then
                  norm=1
               elseif (nfzsym.eq.1.and.k.le.2) then
                  norm=real(nz/2+1)
               elseif (nfzsym.eq.1.and.k.eq.3) then
                  norm=real(nz/2-1)
               else
                  norm=nz
               end if
c
c     Accumulate fft in time in spec
c
               do j=1,nfr
                  do z=1,nzc+nfzsym
                     do x=1,nx
                        spec(x,y,z,k,j)=spec(x,y,z,k,j)+
     &                       pxz(x,z)*2./norm*
     &                       exp((0.,1.)*2.*pi*real(j-1)*tf(i)/fper)*
     &                       dt(i)/fper
                     end do
                  end do
               end do
            end do
         end do
      end do
c
c     New functionality: YZ PLOT
c
      write(*,*)' yz plot?'
      read(*,*) ans
      if (ans.eq.'y') then
 3035    continue
         write(*,*) 'What x value?'
         read(*,*) xpro

         if (xpro-xs.lt.0) then
            xprof=int((xpro-xs)/xl*float(nx)-.5)+nx/2+1
         else
            xprof=int((xpro-xs)/xl*float(nx)+.5)+nx/2+1
         end if

         xpro=xl/float(nx)*float(xprof-nx/2-1)+xs
         write(*,*)xpro,xprof,nx,xl

         write(*,*)'Nearest plane x=',xpro,' plane number',xprof
         xprof=mod(10*nx+xprof-1,nx)+1
         write(*,*)'Give frequency'
         read(*,*)kf
         j=kf+1
         write(*,*)'Give velocity component'
         write(*,*)'1  u'
         write(*,*)'2  v'
         write(*,*)'3  w'
         read(*,*)ivar
c
c     Type of contour plot
c
         logp=.false.
         write(*,*) 'Number of contours ? (0 to select manually'
     &        //'<0 gives logplot)'
         read(*,*) ncont
         if (ncont.lt.0) then
            ncont=abs(ncont)
            logp=.true.
         end if
         ncont=-ncont
         if (ncont.eq.0) then
            write(*,*)  'Give start,end,spacing for contour levels'
     &           //'spacing<0 gives logplot'
            read(*,*) blev,elev,dlev
            if (dlev.lt.0) then
               dlev=abs(dlev)
               logp=.true.
            end if
            ncont=int((elev-blev)/dlev+1.5)
            elev=blev+dlev*float(ncont-1)
            write(*,*) 'endlevel ',elev
         end if
c
c     Create plane
c
         do y=1,nyp
            do z=1,nz

               planex(z-1,nyp+1-y)=abs(spec(xprof,y,z,ivar,j))*norm/2
c     LB /norm perche' in dft c'e' un *2/norm di troppo!!!!!
               if (z.eq.7)write(43,*)planex(z-1,nyp+1-y),y,nyp+1-y
            end do
c
c     Extend periodically in the yz-direction
c
            planex(nz,nyp+1-y)=planex(0,nyp+1-y)
         end do
c
c     Grid
c
         do y=1,nyp
            do zp=0,nz
               write(29,*)planex(zp,y)
            end do
         end do

         do y=1,nyp
            do zp=1,nz+1
               z=zp
               gridx(z-1,nyp+1-y,1)=float(2*zp-nz-2)/float(2*nz)*zl
               gridx(z-1,nyp+1-y,2)=eta(y)
            end do
         end do
         mbox=1
 9832    format (' om = ',i4,' Re = ',f6.0,'$')
         write(headb,9832)kf,re
         if (ivar.eq.1) heada='u$'
         if (ivar.eq.2) heada='v$'
         if (ivar.eq.3) heada='w$'
         call contfou(nz+1,nyp,planex,gridx,ncont,blev,elev,'z$','y$',
     &        logp,heada,headb,mbox,1)


         write(*,*)'New yz plot?'
         read(*,*)ans
         if (ans.eq.'y') goto 3035
      end if

      write(*,*)'Max u plot?'
      read(*,*) ans
      if (ans.eq.'y') then
 3081    write(*,*)'Give frequency'
         read(*,*)kf
         j=kf+1
         xr=0.
         if (fltype.eq.4) xr=2./3.*t
         if (fltype.gt.5) xr=2./3.*t
c
c     Special for spatial boundary layer
c
         xr=xl/2.
         dx=xl/nx
         mxr=int((xr-xs)/dx+.5)
         if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
         mxr=int((xr-xs)/dx+.5)
         do ll=1,3
            do x=1,nx
               xpr=mod(10*nx+x+mxr-1,nx)+1
               maxmod(x,ll)=0.
               do y=1,nyp
                  do z=1,nz
                     maxmod(x,ll)=max(abs(spec(xpr,y,z,ll,j)),
     &                    maxmod(x,ll))
c
c     y and z position of max
c
                     if (ll.eq.1) then
                        if (maxmod(x,ll).eq.
     &                      abs(spec(xpr,y,z,ll,j))) then
                           maxps(x,1)=eta(y)
                           maxps(x,2)=float(2*z-nz-2)/float(2*nz)*zl
                           maxps(x,3)=z
                        end if
                     end if
                  end do
               end do
c     _LB vedi sopra per questo norm!!!!
c     _LBLBLBLB             maxmod(x,ll)=maxmod(x,ll)*norm
               xx(x,ll)=float(x)/float(nx)*xl
               if (ll.eq.1) then
c     write(88,*)xx(x,1),maxps(x,1),maxps(x,2),
c     &maxmod(x,1)
               end if
            end do
            xx(0,ll)=0.
            maxmod(0,ll)=maxmod(nx,ll)
         end do


         do x=1,nx
            do y=1,nyp
               xpr=mod(10*nx+x+mxr-1,nx)+1
               z=maxps(x,3)
               u(y)=abs(spec(xpr,y,z,1,j))
            end do
            do y=1,nyp
               yex(y)=0.
               xex(y)=0.
            end do
            call spline(eta,u,y2,nyp,xex,yex,nex)
            maxps(x,1)=0.
            maxmod(x,1)=0.
            do i=1,nex
               if (abs(maxmod(x,1)).lt.abs(yex(i))) then
                  maxps(x,1)=xex(i)
                  maxmod(x,1)=yex(i)
               end if
            end do
            write(86,*)xx(x,1),maxps(x,1),maxps(x,2),
     &           maxmod(x,3)
         end do


         maxmod(0,1)=maxmod(nx,1)
         maxps(0,1)=maxps(nx,1)
         maxps(0,2)=maxps(nx,2)
         xaxis='x$'
         yaxis='u v w$'
         heada='max of amplitude of u along x$'
         write (headb,9817) kf,re
 9817    format (' om = ',i4,' Re = ',f6.0,'$')
         headc='solid: u, dash: v, dot: w$'
         npoint(1)=nx+1
         npoint(2)=nx+1
         npoint(3)=nx+1
         nplot=3
         xmin=xx(0,1)
         xmax=xx(nx,1)
         write(*,*)'Do you want to make a subarea section?'
         read(*,*)ans
         if (ans.eq.'y') then
            write(*,*)'xmin =  '
            read(*,*)xmin
            write(*,*)'xmax =  '
            read(*,*)xmax
         end if

         ymin=0.
         ymax=0.
         mbox=1
         idash=1
         call rita1(xx,maxmod,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &        xaxis,yaxis,heada,headb,headc,mbox,idash)
         do x=0,nx
            write(30,*) xx(x,1),maxmod(x,1),maxmod(x,2),maxmod(x,3)
         end do
c
c     Plot the curve of the maximum over y of the
c     chosen variable and growth rate
c
         call foustr(xx,maxmod,xr,xl,heada,headb,re,mbox,ihard)

         write(*,*)'New max u plot?'
         read(*,*)ans
         if (ans.eq.'y') goto 3081
      end if
c
c     Transform to fourier space in spanwise direction
c
      if (nz.gt.2) then
         do j=1,nfr
            do k=1,3
               if (nfzsym.eq.0) then
                  call vcfftf(rspec(1,1,1,k,j),rspec(2,1,1,k,j),
     &                 ws(1),ws(2),nz,nx*nyp,nx*nyp*2,2,prez)
               else
                  if (k.le.2) then
                     call vcffts(rspec(1,1,1,k,j),rspec(2,1,1,k,j),
     &                    ws(1),nz/2+1,nx*nyp,nx*nyp*2,2,pres)
                  else
                     call vcftaf(rspec(1,1,2,k,j),rspec(2,1,2,k,j),
     &                    ws(1),ws(2),nz/2-1,nx*nyp,nx*nyp*2,2,prea)
                  end if
               end if
            end do
         end do
      end if
c
c     Some more new stuff from fouLX
c
      write(*,*) 'Amplitude plots?'
      read(*,*)ans
      if (ans.eq.'y') then
         xr=0.
         if (fltype.eq.4) xr=2./3.*t
         if (fltype.gt.5) xr=2./3.*t
c
c     Special for spatial boundary layer
c
         xr=xl/2.
         dx=xl/nx
         mxr=int((xr-xs)/dx+.5)
         if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
         mxr=int((xr-xs)/dx+.5)

         do x=1,nx
            xpr=mod(10*nx+x+mxr-1,nx)+1
            maxmd=0.
            xpro=float(x)/float(nx)*xl
c
c     Create plane
c
            j=1
            do y=1,nyp
               do z=1,nz+1
                  planex(z-1,nyp+1-y)=0.
                  planex(z-1,nyp+1-y)=real(spec(xpr,y,1,1,j))/2
               end do
            end do
            do y=1,nyp
               do zp=1,nz+1
                  z=zp
                  giu=float(2*zp-nz-2)/float(2*nz)*zl
                  do i=2,nz/2+1*(1-nfzsym)
                     planex(z-1,nyp+1-y)=planex(z-1,nyp+1-y)+real(
c     _     &                    spec(xpr,y,i,1,j))*
     &                    spec(xpr,y,i,1,j)*
c     _     &                    cos(2*pi/zl*giu*(i-1)+pi*(i-1))
     &                    exp(ui*2*pi/zl*giu*(i-1)+ui*pi*(i-1)))
                  end do
               end do
            end do
            maxamp=-100000.
            minamp=100000.
            ampy=-202020202
            highy=-100000.
            lowy=100000.
            highz0=-100000.
            lowz0=100000.
            highzp=-100000.
            lowzp=100000.
            do y=1,nyp
               maxampy=-100000.
               minampy=100000.
               do z=0,nz
                  maxamp=max(planex(z,y),maxamp)
                  minamp=min(planex(z,y),minamp)
                  maxampy=max(planex(z,y),maxampy)
                  minampy=min(planex(z,y),minampy)
c
c     y and z position of max and min
c
                  if (maxamp.eq.planex(z,y)) then
                     psmax(1)=eta(nyp+1-y)
                     psmax(2)=float(2*z-nz-2)/float(2*nz)*zl
                     psmax(3)=z
                  end if
                  if (minamp.eq.planex(z,y)) then
                     psmin(1)=eta(nyp+1-y)
                     psmin(2)=float(2*z-nz-2)/float(2*nz)*zl
                     psmin(3)=z
                  end if
               end do
               ampy=max((maxampy-minampy)/2.,ampy)
               highy=max(maxampy,highy)
               lowy=min(minampy,lowy)
               if (ampy.eq.(maxampy-minampy)/2) then
                  psmaxy(1)=eta(nyp+1-y)
               end if
               if (highy.eq.maxampy) then
                  psmaxy(2)=eta(nyp+1-y)
               end if
               if (lowy.eq.minampy) then
                  psmaxy(3)=eta(nyp+1-y)
               end if

c               highz0=max(planex(0,y),highz0)
               if (highz0.le.planex(0,y)) then
                  highz0=planex(0,y)
c               if (highz0.eq.planex(0,y)) then
                  psmaxz(1)=eta(nyp+1-y)
               end if
c               lowz0=min(planex(0,y),lowz0)
c               if (lowz0.eq.planex(0,y)) then
               if (lowz0.ge.planex(0,y)) then
                  lowz0=planex(0,y)
                  psmaxz(2)=eta(nyp+1-y)
               end if

               if (highzp.le.planex(nz/2,y)) then
                  highzp=planex(nz/2,y)
                  psmaxz(3)=eta(nyp+1-y)
               end if
               if (lowzp.ge.planex(nz/2,y)) then
                  lowzp=planex(nz/2,y)
                  psmaxz(4)=eta(nyp+1-y)
               end if
            end do
            write(44,*)xpro,(maxamp-minamp)/2.,maxamp,minamp
            write(45,*)xpro,maxamp,minamp,psmax,psmin
            write(46,*)xpro,ampy,psmaxy(1)
            write(47,'(10e18.9)')xpro,ampy,highy,lowy,
     &           psmaxy(2),psmaxy(3)
            write(48,'(15e18.9)')xpro,ampy,highz0,lowz0,
     &           psmaxz(1),psmaxz(2),highzp,lowzp,
     &           psmaxz(3),psmaxz(4)
         end do

         close(44)
         close(45)
         close(46)

 3090    write(*,*) 'What x value?'
         read(*,*) xpro

         if (xpro-xs.lt.0) then
            xprof=int((xpro-xs)/xl*float(nx)-.5)+nx/2+1
         else
            xprof=int((xpro-xs)/xl*float(nx)+.5)+nx/2+1
         end if

         xpro=xl/float(nx)*float(xprof-nx/2-1)+xs
         write(*,*)xpro,xprof,nx,xl

         write(*,*)'Nearest plane x=',xpro,' plane number',xprof
         xprof=mod(10*nx+xprof-1,nx)+1
         logp=.false.
         write(*,*) 'Number of contours ? (0 to select manually
     &        <0 gives logplot)'
         read(*,*) ncont
         if (ncont.lt.0) then
            ncont=abs(ncont)
            logp=.true.
         end if
         ncont=-ncont
         if (ncont.eq.0) then
            write(*,*) 'Give start,end,spacing for contour levels'
     &           ,'spacing<0 gives logplot'
            read(*,*) blev,elev,dlev
            if (dlev.lt.0) then
               dlev=abs(dlev)
               logp=.true.
            end if
            ncont=int((elev-blev)/dlev+1.5)
            elev=blev+dlev*float(ncont-1)
            write(*,*) 'End level ',elev
         end if
c
c     Create plane
c
         j=1
         do y=1,nyp
            do z=1,nz+1
               planex(z-1,nyp+1-y)=0.
               planex(z-1,nyp+1-y)=real(spec(xprof,y,1,1,j))/2
            end do
         end do
         do y=1,nyp
            do zp=1,nz+1
               z=zp
               giu=float(2*zp-nz-2)/float(2*nz)*zl
               do i=2,nz/2+1
                  planex(z-1,nyp+1-y)=planex(z-1,nyp+1-y)+real(
c     _     &                     spec(xprof,y,i,1,j))*
     &                 spec(xprof,y,i,1,j)*
c     _     &                     cos(2*pi/zl*giu*(i-1)+pi*(i-1))
     &                 exp(ui*2*pi/zl*giu*(i-1)+ui*pi*(i-1)))
               end do
            end do
         end do
         maxamp=-100000.
         minamp=100000.

         do y=1,nyp
            do z=0,nz
               maxamp=max(planex(z,y),maxamp)
               minamp=min(planex(z,y),minamp)
c
c     y and z position of max and min
c
               if (maxamp.eq.planex(z,y)) then
                  psmax(1)=eta(nyp+1-y)
                  psmax(2)=float(2*z-nz-2)/float(2*nz)*zl
                  psmax(3)=z
               end if
               if (minamp.eq.planex(z,y)) then
                  psmin(1)=eta(nyp+1-y)
                  psmin(2)=float(2*z-nz-2)/float(2*nz)*zl
                  psmin(3)=z
               end if
            end do
         end do
         write(*,*) maxamp,minamp,psmax,psmin,
     &        (maxamp-minamp)/2
         write(*,*) xpro,(maxamp-minamp)/2

         do y=1,nyp
            do zp=1,nz+1
               z=zp
               gridx(z-1,nyp+1-y,1)=float(2*zp-nz-2)/float(2*nz)*zl
               gridx(z-1,nyp+1-y,2)=eta(y)
            end do
         end do
         call contfou(nz+1,nyp,planex,gridx,ncont,blev,elev,'z$','y$',
     &        logp,heada,headb,mbox,-1)
         do y=1,nyp
            do zp=1,nz+1
               write(29,*)planex(zp-1,y)
            end do
         end do
         write(*,*)'new amplitude?'
         read(*,*)ans
         if (ans.eq.'y') goto 3090
      end if

c
c     Evaluation
c
      eps=1.e-19
 1234 write(*,*) 'Fundamental frequency ',2*pi/fper
      write(*,*) '  maximum integer multiplier kf ',nfr-1
      write(*,*) 'Fundamental wavenumber ',beta(2)
      write(*,*) '  maximum abs integer multiplier kz ',nz/2-1
      write(*,*) ' '
      write(*,*) ' Give number of (frequency,wavenumber)-pairs '
      write(*,*) '                maximum = ',ikfkzm
      write(*,*) '  0 ends the program '
      write(*,*) '  1 gives urms,vrms,wrms for 1 mode'
      write(*,*) '  n gives urms,vrms,wrms for n modes'
      write(*,*) ' -1 gives xy contourplot '
      write(*,*) ' -2 branch I,II of neutral stability curve'
      read(*,*) ikfkz
      if (ikfkz.eq.0) goto 9999
c
c     xr is the coordinate of the center of the plot
c     xs is the coordinate in the center of the computational box
c
      xr=0.
      if (fltype.eq.4) xr=2./3.*t
      if (fltype.gt.5) xr=2./3.*t
c
c     Special for spatial boundary layer
c
      xr=xl/2.
      dx=xl/nx
      mxr=int((xr-xs)/dx+.5)
      if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)

      if (ikfkz.lt.0) then
         if (ikfkz.ne.-2) then
            write(*,*) ' Give  frequency and wavenumber '
            read(*,*) kf,kz
         else
            kf=1
            kz=0
         end if

         j=kf+1
         if (nfzsym.eq.1) kz=abs(kz)
         if (kz.ge.0) z=kz+1
         if (kz.lt.0) z=nz+1-abs(kz)
c
c     Frequency is in j, wavenumber in z
c
         write(*,*) 'Variable ?'
         write(*,*) '1 u'
         write(*,*) '2 v'
         write(*,*) '3 w'
         write(*,*) '4 energy'
         write(*,*) '5 y-profiles'
         read(*,*) ivar
c
c     Some more new stuff from fouLX
c

         if (ivar.eq.5) then
c
c     Plot y-profiles
c
c     First, plot maximum in y over x
c     ATTENTION: coordinates corrected to comply with pxyst
c                looping only over selected spanwise wavenumber
c                scaling with sqrt(2)/2 for rms

            maxmod = 0.
            do ll=1,3
c     do lll=1,nzc
               do lll=z,z
                  do x=1,nx
c     xpr=mod(10*nx+x+mxr-1,nx)+1
                     xpr=mod(10*nx+x+mxr,nx)+1
                     maxmd=0.
                     do y=1,nyp
                        maxmd=max(abs(spec(xpr,y,lll,ll,j)),maxmd)
                     end do
                     xx(x,ll)=float(x)/float(nx)*xl
                     maxmod(x,ll)=maxmod(x,ll)+sqrt(2.)/2.*maxmd
                  end do
               end do

               xx(0,ll)=0
               maxmod(0,ll)=maxmod(nx,ll)
            end do

            xaxis='x$'
            yaxis='u v w$'
            heada='maximum amplitude along x$'
            write (headb,9819) kf,kz,re
 9819       format (' om = ',i4,', beta = ',i4,', Re = ',f6.0,'$')
            headc='solid: u, dash: v, dot: w$'

            npoint(1)=nx+1
            npoint(2)=nx+1
            npoint(3)=nx+1
            nplot=3
            xmin=xx(0,1)
            xmax=xx(nx,1)
            write(*,*) 'Do you want to make a subarea selection?'
            read(*,*) ans
            if (ans.eq.'y') then
               write(*,*)'xmin =  '
               read(*,*)xmin
               write(*,*)'xmax =  '
               read(*,*)xmax
            end if
            ymin=0.
            ymax=0.
            mbox=1
            idash=1
            call rita1(xx,maxmod,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)

            write(*,*) 'Wrote maximum data to fort.82'
            do x=0,nx
               write(82,'(4e18.9)') xx(x,1),(maxmod(x,ll),ll=1,3)
            end do
            close(82)

            write(*,*) 'What x value?'
            read(*,*) xpro

            if (xpro-xs.lt.0) then
               xprof=int((xpro-xs)/xl*float(nx)-.5)+nx/2+1
            else
               xprof=int((xpro-xs)/xl*float(nx)+.5)+nx/2+1
            end if

            xpro=xl/float(nx)*float(xprof-nx/2-1)+xs

            write(*,*)'Nearest plane x=',xpro,' plane number',xprof

            xprof=mod(10*nx+xprof-1,nx)+1

            do ll=1,3
               do y=1,nyp
                  modu=abs(spec(xprof,nyp+1-y,z,ll,j))

cfr          fase=atan(aimag(spec(xprof,nyp+1-y,z,ll,j))/
cfr     &              real(spec(xprof,nyp+1-y,z,ll,j)))
cfr
cfr          fase=fase/pi*180
c          if ((fase.gt.0.).and.(real(spec(xprof,nyp+1-y,z,ll,j))
c     &        .lt.0.)) fase=fase-180
c
c          if ((fase.lt.0.).and.(real(spec(xprof,nyp+1-y,z,ll,j))
c     &        .lt.0.)) fase=fase+180
c
cfr          write(66,*) eta(nyp+1-y),modu,fase
cfr       write(67,*)eta(nyp+1-y),real(spec(xprof,nyp+1-y,z,ll,j))
cfr     &          ,aimag(spec(xprof,nyp+1-y,z,ll,j))

                  umod(y,ll)=modu*sqrt(2.)/2.
                  yy(y,ll)=eta(nyp+1-y)
               end do
            end do

            write(*,*) 'Wrote complete plane to fort.69'
            do y=1,nyp
               write(69,*) eta(nyp+1-y)
            end do
            do y=1,nyp
               modu=real(spec(xprof,nyp+1-y,1,1,j))
               write(69,*) modu/2
            end do
            do lll=2,nzc
               do y=1,nyp
                  modu=real(spec(xprof,nyp+1-y,lll,1,j))
                  write(69,*) modu
               end do
            end do
            close(69)

            xaxis='u,v,w$'
            yaxis='y$'
            heada='wave form in y$'
            write (headb,9719) kf,kz,re
 9719       format (' om = ',i4,', beta = ',i4,', Re = ',f6.0,'$')
            headc='solid: u, dash: v, dot: w$'

            npoint(1)=nyp
            npoint(2)=nyp
            npoint(3)=nyp
            npoint(4)=nyp
            nplot=3
            xmin=0.
            xmax=0.
            ymin=0.
            ymax=yy(nyp,1)
            mbox=1
            idash=1
            call rita1(umod,yy,xmin,xmax,ymin,ymax,npoint,nplot,nyp,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)

            write(*,*) 'Write profile to fort.81'
            do y=1,nyp
               write(81,'(5e18.9)') yy(y,1),(umod(y,ll),ll=1,3)
            end do
            close(81)

            goto 1234
c     end plot y-profiles
         end if
c
c     Get plot dependent parameters
c
         logp=.false.
         if (ikfkz.ne.-2) then
         write(*,*) 'Number of contours ? (0 to select manually
     &<0 gives logplot)'
            read(*,*) ncont
            if (ncont.lt.0) then
               ncont=abs(ncont)
               logp=.true.
            end if
            ncont=-ncont
            if (ncont.eq.0) then
              write(*,*)  'Give start,end,spacing for contour levels'
     &,'spacing<0 gives logplot'
              read(*,*) blev,elev,dlev
              if (dlev.lt.0) then
                dlev=abs(dlev)
                logp=.true.
              end if
              ncont=int((elev-blev)/dlev+1.5)
              elev=blev+dlev*real(ncont-1)
              write(*,*) 'endlevel ',elev
            end if
          end if
c
c     Make plane
c     Note how the plane is turned upside down to compensate
c     for reverse numbering
c
          if (ivar.ne.4.) then
            do y=1,nyp
               do x=1,nx
                  xpr=mod(10*nx+x+mxr-1,nx)+1
                  plane(x,nyp+1-y)=abs(spec(xpr,y,z,ivar,j))
               end do
c
c     Extend periodically in the yx-direction
c
               plane(0,nyp+1-y)=plane(nx,nyp+1-y)
            end do
            if (ivar.eq.1) heada='u$'
            if (ivar.eq.2) heada='v$'
            if (ivar.eq.3) heada='w$'
         else
            do y=1,nyp
               do x=1,nx
                  xpr=mod(10*nx+x+mxr-1,nx)+1
                  plane(x,nyp+1-y)=0.5*(abs(spec(xpr,y,z,1,j))**2
     &          +abs(spec(xpr,y,z,2,j))**2+abs(spec(xpr,y,z,3,j))**2)
               end do
c
c     Extend periodically in the yx-direction
c
               plane(0,nyp+1-y)=plane(nx,nyp+1-y)
            end do
            heada='Energy$'
         end if
         if (logp) then
            do y=1,nyp
               do x=0,nx
                  plane(x,y)=-log10(plane(x,y)+eps)
               end do
            end do
         end if
         if (ikfkz.ne.-2) then
c
c     Make grid
c
            do y=1,nyp
               do x=0,nx
                  grid(x,nyp+1-y,1)=real(x-nx/2)/real(nx)*xl+xr
                  grid(x,nyp+1-y,2)=eta(y)
               end do
            end do
            write (headb,9919) kf,kz,re
            mbox=1

            call contfou(nx+1,nyp,plane,grid,ncont,blev,elev,'x$','y$',
     &           logp,heada,headb,mbox,-1)
            goto 1234
         else
            write (headb,9919) kf,kz,re
            mbox=1
            call foubranch(plane,eta,xr,xl,heada,headb,re,mbox,-1)
            goto 1234
         end if
      end if
c
c     End if ikfkz<0
c
      do iii=1,ikfkz
         write(*,*) ' Give ',iii,' frequency and wavenumber '
         read(*,*) kfa(iii),kza(iii)
      end do

      do iii=1,ikfkz
         kf=kfa(iii)
         kz=kza(iii)
         j=kf+1
         if (nfzsym.eq.1) kz=abs(kz)
         if (kz.ge.0) z=kz+1
         if (kz.lt.0) z=nz+1-abs(kz)

         do x=1,nx
            urms(x)=0.
            vrms(x)=0.
            wrms(x)=0.
            do y=1,nyp
               urms(x)=urms(x)+abs(spec(x,y,z,1,j))**2*wint(y)*0.5
               vrms(x)=vrms(x)+abs(spec(x,y,z,2,j))**2*wint(y)*0.5
               wrms(x)=wrms(x)+abs(spec(x,y,z,3,j))**2*wint(y)*0.5
            end do
            e(x,iii)=0.5*(urms(x)+vrms(x)+wrms(x))
c
c     In order to get non-zero values...
c
c            e(x,iii) = e(x,iii) + eps
c            urms(x)  = urms(x)  + eps
c            vrms(x)  = vrms(x)  + eps
c            wrms(x)  = wrms(x)  + eps
         end do
      end do
c
c     More new stuff from fouLX
c

c___ikfkz=1 one may also get the energy integrated in both y and z
      write(*,*)'Would you like the energy integrated in both y and z?'
      write(*,*)'Valid only if ikfkz=1.'
      read(*,*) ans
      if (ans.eq.'y') then
         write(*,*)'Output file fort.80 contains total energy +'
         write(*,*)'urms, vrms and wrms at the selected spanwise ' //
     &        'wavenumber'
         do x=1,nx
            e(x,1)=0
         end do
         do z=1,nz
            do x=1,nx
               efu(x,1)=0
               efu(x,2)=0
               efu(x,3)=0
               do y=1,nyp
                  efu(x,1)=efu(x,1)+abs(spec(x,y,z,1,j))**2*wint(y)*0.5
                  efu(x,2)=efu(x,2)+abs(spec(x,y,z,2,j))**2*wint(y)*0.5
                  efu(x,3)=efu(x,3)+abs(spec(x,y,z,3,j))**2*wint(y)*0.5
               end do
c               e(x,1)=e(x,1)+sqrt(urms(x))
               e(x,1)=e(x,1)+0.5*(efu(x,1)+efu(x,2)+efu(x,3))
            end do
         end do
      end if

      if (ikfkz.eq.1) then
         do x=1,nx
c
c     Corrected to be consistent with pxyst
c
c     xpr=mod(10*nx+x+mxr-1,nx)+1
            xpr=mod(10*nx+x+mxr,nx)+1
            xx(x,1)=real(x)/real(nx)*xl
            xx(x,2)=xx(x,1)
            xx(x,3)=xx(x,1)
            xx(x,4)=xx(x,1)
            efu(x,1)=e(xpr,1)
            efu(x,2)=0.5*urms(xpr)
            efu(x,3)=0.5*vrms(xpr)
            efu(x,4)=0.5*wrms(xpr)
         end do
         efu(0,1)=efu(nx,1)
         efu(0,2)=efu(nx,2)
         efu(0,3)=efu(nx,3)
         efu(0,4)=efu(nx,4)
         xx(0,1)=0.
         xx(0,2)=0.
         xx(0,3)=0.
         xx(0,4)=0.

         write(*,*) 'Plot data written to fort.80'
         do x=0,nx
            write(80,'(5e18.9)')
     &           xx(x,1),efu(x,1),efu(x,2),efu(x,3),efu(x,4)
         end do
         close(80)

         yaxis='energy$'
         heada='total energy, energy per velocity component$'
         write (headb,9919) kf,kz,re
 9919    format (' om = ',i4,', beta = ',i4,', Re = ',f6.0,'$')
         headc='solid: energy, dash: Eu, dot: Ev, chdsh: Ew $'
         nplot=4
c
         xmin=xx(0,1)
         xmax=xx(nx,1)
         xaxis='x$'
         npoint(1)=nx+1
         npoint(2)=nx+1
         npoint(3)=nx+1
         npoint(4)=nx+1
         idash=1
         mbox=1
c
         write(*,*) 'Set type of energy axis'
         write(*,*) '0 linear'
         write(*,*) '1 logarithmic'
         read(*,*) iloga
         if (iloga.eq.0) then
            ymin=0.
            ymax=0.
            call rita1(xx,efu,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)
         else
            ymin=eps
            ymax=5.
            ilog = 2
            call ritlog(xx,efu,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &           xaxis,yaxis,heada,headb,headc,ilog,idash)
         end if
c
c        Write out data to file
c
         call savelist(1,xx,efu,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &        xaxis,yaxis,heada,headb,headc,2,idash)

         goto 1234
      else
         ymin=eps
         ymax=5.
         do iii=1,ikfkz
            do x=1,nx
               xpr=mod(10*nx+x+mxr-1,nx)+1
               xx(x,iii)=real(x)/real(nx)*xl
               efu(x,iii)=max(e(xpr,iii),ymin)
            end do
            efu(0,iii)=efu(nx,iii)
            xx(0,iii)=0.
            npoint(iii)=nx+1
         end do

         write(*,*) 'Plot data written to fort.80'
         do x=0,nx
            write(80,'(10e18.9)') xx(x,1),(efu(x,i),i=1,ikfkz)
         end do
         close(80)

         yaxis='log(energy)$'
         heada='Energy'
         write (headb,9920) re
 9920    format (' Re = ',f6.0,'$')
         headc=' '
         do iii=1,ikfkz
            write(mode(iii),9921) kfa(iii),kza(iii)
 9921       format (' (',i1,',',i1,')')
         end do
         write(headb,9923) (mode(k),k=1,ikfkz)
 9923    format('solid:',a6,', dash:',a6,', dot:',a6,', chdsh:',10a6)
         nplot=ikfkz

         xmin=xx(0,1)
         xmax=xx(nx,1)
         xaxis='x$'
         idash=1
         ilog=2

         call ritlog(xx,efu,xmin,xmax,ymin,ymax,npoint,nplot,nx+1,
     &        xaxis,yaxis,heada,headb,headc,ilog,idash)

         goto 1234
      end if

 9999 stop

      end program fou
