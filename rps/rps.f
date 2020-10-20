c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program rps
c
c     Simple interactive plane sequence plotting program
c
      implicit none

      include 'par.f'
c
c     Main storage
c
      integer maxmem,npl,nzn
      parameter (maxmem=3500000)
      parameter (nzn=nz-(nz/2-1)*nfzsym)
      parameter (npl=maxmem/(nx*(nyp+nzn)))
c      parameter (npl=3000)
      real uxys(memnx*2,memny,npl)
      real uxzs(memnx*2,nzn,npl)
      real uyzs(memny,nzn,npl)
      equivalence(uxys,uxzs,uyzs)
c
c     Plot plane
c
      integer plsz
      parameter (plsz=(nx+2)*(nz+nyp)+nyp*(nz+1))
      real plane(plsz),grid(plsz,2),plmean(plsz),plrms(plsz)
      real w1(plsz),w2(plsz)
      real xy(nx+1,nyp,2),xz(nx+1,nz+1,2),yz(nyp,nz+1,2)
      equivalence (grid,xy,xz,yz)
c
c     Graphs
c
      real graf(nx+nz),gxax(nx+nz),grafm(nx+nz)
      real fx(nx+1,npl),xc(nx+1,npl)
c
c     Transform preprocessing information
c
      real prex(nx+15),prezr(nz+15),prey(nyp*2+15)
c
      character*1 ans, kans
      character*80 namnin
      character*80 timenam
c
c     Geometrics
c
      real eta(nyp)
      real pi,xl,zl,dstar,xr,t,xs,xss,re,xsa(npl),ta(npl),upan,cleft
      real xe,ys,ye,zs,ze
      complex ialfa(nx/2+1)
      parameter (pi = 3.1415926535897932385)
c
c     Flow
c
      integer fltype
c
c     Plot plane data
c
      integer tpl,mpl,ivar,idxpr,idzpr,npr
      real cpl,xpr,ypr,zpr,dxpr,dzpr,xmoff
c
c     Time interval for processing/displaying
c
      integer nint,itint(npl),ipln
      real tint(npl)
c
c     Probe data
c
      real upr(npl,100),tpr(npl,100),uprm(100)
c
c     Spectral analysis
c
      real esp(npl/2,100,2),fesp(npl/2,100)
      integer nbin
      logical lhanw,lper
c
c     Length data
c
      real dlen(npl),tlen(npl)
      real relth
c
c     Plot parameters
c
      logical lraw,lcmean,chgdef,lsubm,logax
      real gl,gu
      integer ncont,jhard,mbox,iplot,idudy,iu,idash
      real dlev,elev,blev
      real ymin,ymax,xmin,xmax
      real ommin,ommax,emin,emax
      real tf,tl,ti,tnext
      integer itf,itl,it,ifil
      integer mxs,mxe,mye,mys,mze,mzs
      integer npoint(npl+100)
      character*80 heada,headb,headc
      character*12 var
      character*32 yaxis
      real umean
c
c     General
c
      integer x,y,z,i,j
      real sym
c
c     Special
c
      integer iix,k,kk
      real xic
      integer sxe,sxs,sze,szs
      integer high
      real umin,umax
      logical nothing,timefile
      logical ttmean
      real esi(nx+1,npl),esigr(nx+1,npl,2)

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         rps $Rev$   '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Max planes: ',npl
c
c     Initialize x,z-transforms
c
      call vrffti(nx,prex,0)
c
c     Prepare real transform in the z-direction
c
      if (mod(nz,2).ne.1) call vrffti(nz,prezr,0)
c
c     Intialize y-transform
c
      call vcosti(nyp,prey,0)
c
      do y=1,nyp
         eta(y)=cos(pi*real(y-1)/real(nyp-1))
      end do
c
      write(*,*)
      write(*,*) 'Give file to plot from'
      read(*,9000) namnin
 9000 format(a32)
      call rseq(uxys,uxzs,uyzs,re,fltype,xl,zl,dstar,t,xs,namnin,
     &     w1,w1,w1,
     &     xsa,ta,tpl,ivar,cpl,mpl,npl,nzn,
     &     umin,umax)
      if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6) then
         do  y=1,nyp
            eta(y)=(1.+eta(y))/dstar
         end do
      end if
      write(*,9010) mpl,ta(1),ta(mpl)
 9010 format(' read ',i5,' planes, t= ',f8.2,' - ',f8.2)
      write(*,*) 'min/max: ',umin,umax
      do x=1,nx/2+1
         ialfa(x)=(0.,1.)*2.*pi/xl*real(x-1)
      end do
c
 1000 continue
      itf=1
      itl=1
      if (mpl.gt.1)  then
         write(*,*)
         write(*,*) 'New plot'
         write(*,*) 'Give start and end time to display/process'
         write(*,*) 'end < start terminates program'
         read(*,*) tf,tl
         if (tl.lt.tf) stop
         do i=2,mpl
            if (tf.gt..5*(ta(i-1)+ta(i))) itf=i
            if (tl.gt..5*(ta(i-1)+ta(i))) itl=i
         end do
c        write(*,*) 'write time file?'
c        read(*,*) ans
c        timefile = .false.
c        if (ans.eq.'y') then
c           write(*,*) 'give filename'
c           read(*,*) timenam
c           timefile =.true.
c        end if
         timefile=.true.
         timenam = 'time'
      end if
c
c     Default options
c
      write(headb,9040) xl,zl,nx,ny,nz,namnin
9040  format(' xl = ',F7.2,' zl= ',F6.2,' ',I4,'x',I3,'x',I3,
     &       ' file ',A20,'$')
c
c     Create plot
c     Defaults
c
      nothing=.false.
      iplot=1
      jhard=0
      mxs=1
      mxe=nx+1
      mys=1
      mye=nyp
      mzs=1
      mze=nz+1
      high=0
c
c     mbox=1 max plotarea, mbox=0 grid selected proportions
c
      mbox=1
      if (tpl.eq.2) mbox=0
      idudy=0
      if (fltype.ge.4) idudy=2
      iu=0
      ti=0.
      gl=0.
      gu=0.
      chgdef=.false.
      ifil=0
      ommin=0.
      ommax=0.
      emin=0.
      emax=0.
c
c     Shifting to follow disturbance
c
      cleft=0.
      if (fltype.le.3) cleft=-xl/2.
      if (fltype.eq.1) upan=2./3.
      if (fltype.eq.2) upan=0.
      if (fltype.eq.3.and.tpl.eq.3) upan=.55
      if (fltype.eq.3.and.tpl.eq.2) upan=.55
      if (fltype.eq.3.and.tpl.eq.1) upan=.7
      if (fltype.ge.4) upan=0.
c
c     Change defaults
c
      write(*,*) 'Change default options ? (y/n)'
      read(*,9900) ans
 9900 format(a1)
      if (ans.ne.'N'.and.ans.ne.'n') then
         chgdef=.true.
c
         write(*,*) 'Give type of output'
         write(*,*) '-3 tektronix stream file'
         write(*,*) '0 screen'
         write(*,*) '1 screen plots and separate postscriptfiles'
         write(*,*) '4 raw files'
         read(*,*) jhard
         lraw=.false.
         if (jhard.eq.4) then
            lraw=.true.
            jhard=0
         end if
         if (tpl.eq.1) then
            write(*,*) 'Give type of plot'
            write(*,*) '1 raw data'
            write(*,*) '2'
            write(*,*) '3'
            write(*,*) '4 temporal average'
            write(*,*) '5 temporal rms'
            write(*,*) '6 time signal from probe'
            write(*,*) '7 length(time)'
            write(*,*) '8 f(x) for one y or z'
            write(*,*) '9'
            write(*,*) '10 frequency spectra from probe'
            read(*,*) iplot
         end if
         if (tpl.eq.2) then
            write(*,*) 'Give type of plot'
            write(*,*) '1 raw data'
            write(*,*) '2 two point correlation, x-direction'
            write(*,*) '3 two point correlation, z-direction'
            write(*,*) '4 temporal average'
            write(*,*) '5 temporal rms'
            write(*,*) '6 time signal from probe'
            write(*,*) '7 length(time)'
            write(*,*) '8 temporal mean of two point correlation, (x)'
            write(*,*) '9 temporal mean of two point correlation, (z)'
            write(*,*) '10 frequency spectra from probe'
            write(*,*) '11 subtract temporal mean...'
            read(*,*) iplot
            if (iplot.eq.11) then
               write(*,*) 'Run type 4 prior to type 11 with the SAME'
               write(*,*) 'panorating speed!'
               ttmean=.true.
               iplot=1
               open(file='instab',unit=86)
            else
               ttmean=.false.
            end if
            lcmean=iplot.eq.8.or.iplot.eq.9
            if (iplot.eq.8.or.iplot.eq.9) iplot=iplot-6
         end if
         if (tpl.eq.3) then
            write(*,*) 'Give type of plot'
            write(*,*) '1 raw data'
            write(*,*) '2'
            write(*,*) '3'
            write(*,*) '4 temporal average'
            write(*,*) '5 temporal rms'
            write(*,*) '6 time signal from probe'
            write(*,*) '7'
            write(*,*) '8'
            write(*,*) '9'
            write(*,*) '10 frequency spectra from probe'
            read(*,*) iplot
         end if
         mbox=1
         if (iplot.eq.1.or.iplot.eq.4.or.iplot.eq.5) then
            write(*,*) 'Give plot proportions',
     &           ' (0 natural, 1 maximum plotarea)'
            read(*,*) mbox
         end if
         if (tf.ne.tl) then
            write(*,*) 'Give time interval between displayed/processed',
     &           ' planes (0 for every)'
            read(*,*) ti
         end if
         if (tpl.ne.2.and.ivar.eq.1.and.(iplot.eq.1.or.iplot.ge.4)) then
            write(*,*)
     &          'Plot variable (0 u-umean, 1 d(u-umean)dy ,2 u, 3 dudy)'
            if (fltype.eq.1.or.fltype.eq.2.or.
     &          fltype.eq.4.or.fltype.eq.5)
     &           write(*,*) '(4 u-ulam, 5 d(u-ulam)dy)'
            read(*,*) idudy
         end if
         if (tpl.eq.2.and.ivar.eq.1) then
            if (iplot.eq.2.or.iplot.eq.3) then
               iu=1
            else
               write(*,*) 'Plot variable (0 u-ulam, 1 u)'
               read(*,*) iu
               if (iu.eq.0) high=1
            end if
         end if
         if (iplot.ge.4.or.tf.eq.tl) upan=0.
         if (iplot.eq.1.and.tf.ne.tl.and.tpl.ne.3.or.iplot.eq.4) then
            write(*,*) 'Give panorating speed (def ',upan,')'
            read(*,*) upan
         end if
         write(*,*) 'plot anything?'
         read(*,*) ans
         if (ans.eq.'n') then
            nothing=.true.
         else
            nothing=.false.
         end if

         if (iplot.eq.6.or.iplot.eq.10) then
            if (tpl.ne.3) then
               write(*,*) 'Give x-coordinate for first probe'
               read(*,*) xpr
               write(*,*) 'Number of probes'
               read(*,*) npr
               if (npr.gt.1) then
                  write(*,*) 'Distance between probes'
                  read(*,*) dxpr
                  idxpr=int(dxpr/xl*real(nx)+.5)
                  dxpr=idxpr*xl/real(nx)
                  write(*,*) 'Closest value ',dxpr
               end if
            else
               write(*,*) 'Give z-coordinate for first probe'
               read(*,*) zpr
               z=int(zpr/zl*real(nz)+.5+real(nz/2+1))
               zpr=zl/real(nz)*real(z-nz/2-1)
               write(*,9070) zpr,z
               write(*,*) 'Number of probes'
               read(*,*) npr
               if (npr.gt.1) then
                  write(*,*) 'Distance between probes'
                  read(*,*) dzpr
                  idzpr=int(dzpr/zl*real(nz)+.5)
                  dzpr=idzpr*zl/real(nz)
                  write(*,*) 'Closest value ',dzpr
               end if
            end if
            if (tpl.ne.2) then
               write(*,*) 'Give y-coordinate for probe'
               read(*,*) ypr
               if (fltype.eq.1.or.fltype.eq.2
     &              .or.fltype.eq.4.or.fltype.eq.5) then
                  y=int(acos(ypr)/pi*real(nyp-1)+1.5)
               else
                  y=int(acos(ypr*dstar-1.)/pi*real(nyp-1)+1.5)
               end if
               ypr=eta(y)
               write(*,9055) ypr,y
 9055          format(' Nearest point y=',F6.3,' point number',I4)
            end if
            if (tpl.eq.2) then
               write(*,*) 'Give z-coordinate for probe'
               read(*,*) zpr
               z=int(zpr/zl*real(nz)+.5+real(nz/2+1))
               zpr=zl/real(nz)*real(z-nz/2-1)
               write(*,9070) zpr,z
 9070          format(' Nearest point z=',F8.2,' point number',I4)
            end if
c
c     Set the left boundary so that our interesting point is at the center
c     of the box in the x-direction
c
            cleft=xpr
            if (iplot.eq.10) then
               write(*,*) 'periodic equidistant data ? (y/n)'
               read(*,'(a)') ans
               lper=ans.ne.'N'.and.ans.ne.'n'
               lhanw=.false.
               if (.not.lper) then
                  write(*,*) 'use Hanning window ?'
                  read(*,'(a)') ans
                  lhanw=ans.ne.'N'.and.ans.ne.'n'
               end if
               write(*,*) 'how many frequencies to a bin ? ',
     &              '(1 for no binning)'
               read(*,*) nbin
               write(*,*) 'use logarithmic axes ? (y/n)'
               read(*,'(a)') ans
               logax=ans.ne.'N'.and.ans.ne.'n'
            end if
            if (.not.logax.and..not.lhanw)  then
               write(*,*) 'subtract mean over time ? (y/n)'
               read(*,'(a)') ans
               lsubm=ans.ne.'N'.and.ans.ne.'n'
            else
               lsubm=lhanw
            end if
            write(*,*) 'select range ? (y/n)'
            read(*,'(a)') ans
            if (ans.ne.'N'.and.ans.ne.'n') then
               write(*,*) 'Give omegamin,omegamax'
               read(*,*) ommin,ommax
               write(*,*) 'Give emin,emax'
               read(*,*) emin,emax
            end if
         end if
         if (iplot.eq.1.or.iplot.eq.4.or.iplot.eq.5) then
            if (upan.eq.0..and.tpl.ne.3) then
               write(*,*) 'Give the coordinate of the left boundary'
               read(*,*) cleft
            end if
         end if
         if (iplot.eq.7) then
            write(*,*) 'Give threshold relative to max'
            read(*,*) relth
         end if
         if (iplot.eq.8) then
            if (tpl.eq.1) then
               write(*,*) 'Give y-coordinate'
               read(*,*) ypr
               if (fltype.eq.1.or.fltype.eq.2) then
                  y=int(acos(ypr)/pi*real(nyp-1)+1.5)
               else
                  y=int(acos(ypr*dstar-1.)/pi*real(nyp-1)+1.5)
               end if
               ypr=eta(y)
               write(*,9055) ypr,y
            end if
            if (tpl.eq.2) then
               write(*,*) 'Give z-coordinate'
               read(*,*) zpr
               z=int(zpr/zl*real(nz)+.5+real(nz/2+1))
               zpr=zl/real(nz)*real(z-nz/2-1)
               write(*,9070) zpr,z
            end if
         end if
c
c     Set levels for raw file
c
         if (lraw) then
            write(*,*) 'Give upper and lower level for greyscale'
            write(*,*) 'Equal for automatic'
            read(*,*) gl,gu
         end if
      end if
c
c     End changing defaults
c
      if (.not.lraw) then
         if (iplot.le.5.and.(tpl.eq.1.or.(tpl.eq.2.and.
     &        (iplot.eq.1.or.iplot.ge.4)).or.tpl.eq.3)) then
            write(*,*) 'Number of contours ? (0 to select manually)'
            read(*,*) ncont
            ncont=-ncont
            if (ncont.eq.0) then
               write(*,*)  'Give start,end,spacing for contour levels'
               read(*,*) blev,elev,dlev
               ncont=int((elev-blev)/dlev+1.5)
               elev=blev+dlev*real(ncont-1)
               write(*,*) 'endlevel ',elev
            end if
         end if
      end if
      sym=1.
      if (ivar.ge.3.and.ivar.le.5) sym=-1.
      if (ivar.eq.1) then
         if (tpl.eq.1.and.idudy.eq.0) var='u-umean'
         if (tpl.eq.1.and.idudy.eq.1) var='d(u-umean)dy'
         if (tpl.eq.1.and.idudy.eq.2) var='u'
         if (tpl.eq.1.and.idudy.eq.3) var='dudy'
         if (tpl.eq.1.and.idudy.eq.4) var='u-ulam'
         if (tpl.eq.1.and.idudy.eq.5) var='d(u-ulam)dy'
         if (tpl.eq.2.and.iu.eq.0) var='u-ulam'
         if (tpl.eq.2.and.iu.eq.1) var='u'
         if (lsubm) var='u-umean'
      end if
      if (ivar.eq.2) var='v'
      if (ivar.eq.2.and.lsubm) var='v-vmean'
      if (ivar.eq.3) var='w'
      if (ivar.eq.3.and.lsubm) var='w-wmean'
      tnext=ta(itf)
      i=0
      do 2015 it=itf,itl
         t=ta(it)
         if (t.lt.tnext) goto 2015
         i=i+1
         tint(i)=t
         itint(i)=it
         tnext=tnext+ti
 2015 continue
      nint=i

      if (timefile) then
         open(file=timenam,unit=10)
         do i=1,nint
            write(10,*) i-1,itint(i),tint(i)
         end do
         close(10)
      end if
c
c     Calculate total mean of velocity for subtraction
c     when calculating mean of correlations
c
      if (tpl.eq.2.and.lcmean) call tmean(umean,tint,nint,tf,tl,itint,
     &     plane,w1,uxzs,npl,sym,nzn)
      write(*,*) 'Keystroke after each plot?'
      read(*,*) kans
      umin=1.E30
      umax=-1.E30
c
c     Main loop
c
      do ipln=1,nint
c
c     Wait for user reply
c
         if (kans.eq.'y') then
            read(*,*) ans
         end if
         it=itint(ipln)
         t=tint(ipln)
         if (nint.gt.1) then
            if (ipln.eq.nint/100+1) then
               write(*,*)
               write(*,*) 'Processed 1%'
            end if
            if (ipln.eq.nint/10+1) write(*,*) 'Processed 10%'
            if (ipln.eq.nint/2+1) write(*,*) 'Processed 50%'
            if (ipln.eq.nint*9/10+1) write(*,*) 'Processed 90%'
         end if
c
c     xr is the coordinate of the center of the plot
c     xs is the coordinate in the center of the computational box
c
         xs=xsa(it)
         xr=upan*t+xl/2.+cleft
c
c     Create plot data
c
         if (tpl.eq.1) call getxys(plane,w1,grid,w2,w2,ialfa,prex,
     &       xs,xr,it,npl,uxys,eta,xl,fltype,dstar,ivar,idudy,prey)
         if (tpl.eq.2.and.(iplot.eq.1.or.iplot.ge.4))
     &        call getxzs(plane,w1,grid,w2,w2,ialfa,prex,
     &        xs,xr,it,sym,npl,uxzs,xl,zl,nzn,ivar,cpl,iu,fltype,high)
         if (tpl.eq.3) call getyzs(plane,w1,grid,it,sym,npl,uyzs,eta,zl,
     &        nzn,fltype,dstar,ivar,idudy,prey)
         if (tpl.eq.2.and.iplot.eq.2) call twoxs(graf,gxax,xmin,xmax,
     &        ymin,ymax,xl,plane,w1,it,npl,nzn,sym,uxzs,prex,lcmean,
     &        umean)
         if (tpl.eq.2.and.iplot.eq.2.and.lcmean)
     &        call accu(graf,grafm,nx/2+1,tint,nint,tf,tl,ipln)
         if (tpl.eq.2.and.iplot.eq.3) call twozs(graf,gxax,xmin,xmax,
     &        ymin,ymax,sym,zl,plane,w1,it,npl,nzn,uxzs,prezr,lcmean,
     &        umean)
         if (tpl.eq.2.and.iplot.eq.3.and.lcmean)
     &        call accu(graf,grafm,nx/2+1,tint,nint,tf,tl,ipln)
         if (iplot.eq.4.or.iplot.eq.5) call statpl(plane,plmean,plrms,
     &        plane,plmean,plrms,plane,plmean,plrms,nzn,tint,nint,tf,
     &        tl,ipln,tpl,sym)
         if (iplot.eq.6.or.iplot.eq.10) call prob(plane,plane,plane,upr,
     &        y,z,tpl,idxpr,idzpr,npr,ipln,npl,t,tpr,lsubm,uprm,tint,
     &        nint,tf,tl)
         if (iplot.eq.7) call dislen(plane,plane,tpl,
     &        dlen(ipln),tlen(ipln),t,relth,nzn,xl)
         if (iplot.eq.8) call fofx(plane,plane,grid,grid,y,z,tpl,
     &        fx(1,ipln),xc(1,ipln))
c
         if (ipln.eq.1.and.(iplot.eq.1.or.iplot.eq.4.or.iplot.eq.5).and.
     &        chgdef) then
            write(*,*) 'Do you want to make a subarea selection (y/n)'
            read(*,7000) ans
 7000       format(a1)
            if (ans.eq.'Y'.or.ans.eq.'y') then
               if (tpl.ne.3) then
                  write(*,9300)xy(mxs,1,1),xy(mxe,1,1)
 9300             format('Give start and end for x (',F8.3,'-',F8.3,')')
                  read(*,*) xss,xe
                  do i=mxs,mxe-1
                     if (xss.lt..5*(xy(i,1,1)+xy(i+1,1,1))) goto 3300
                  end do
                  i=mxe
 3300             mxs=i
                  do i=mxs,mxe-1
                     if (xe.lt..5*(xy(i,1,1)+xy(i+1,1,1))) goto 4010
                  end do
                  i=mxe
 4010             mxe=i
                  write(*,9310)xy(mxs,1,1),xy(mxe,1,1)
 9310             format(' Nearest values  ',F8.3,' - ',F8.3)
               end if
               if (tpl.eq.1) then
                  write(*,9320)xy(1,mys,2),xy(1,mye,2)
 9320         format(' Give start and end for y  (',F8.4,' - ',F8.4,')')
                  read(*,*) ys,ye
                  do i=mys,mye-1
                     if (ys.lt..5*(xy(1,i,2)+xy(1,i+1,2))) goto 4020
                  end do
                  i=mye
 4020             mys=i
                  do i=mys,mye-1
                     if (ye.lt..5*(xy(1,i,2)+xy(1,i+1,2))) goto 4030
                  end do
                  i=mye
 4030             mye=i
                  write(*,9330) xy(1,mys,2),xy(1,mye,2)
 9330             format(' Nearest values  ',F8.4,' - ',F8.4)
               end if
               if (tpl.eq.2) then
                  write(*,9340)xz(1,mzs,2),xz(1,mze,2)
 9340         format(' Give start and end for z  (',F8.3,' - ',F8.3,')')
                  read(*,*) zs,ze
                  do i=mzs,mze-1
                     if (zs.lt..5*(xz(1,i,2)+xz(1,i+1,2))) goto 4040
                  end do
                  i=mze
 4040             mzs=i
                  do i=mzs,mze-1
                     if (ze.lt..5*(xz(1,i,2)+xz(1,i+1,2))) goto 4050
                  end do
                  i=mze
 4050             mze=i
                  write(*,9350)xz(1,mzs,2),xz(1,mze,2)
 9350             format(' Nearest values ',F8.3,' - ',F8.3)
               end if
               if (tpl.eq.3) then
                  write(*,9320)yz(mys,1,2),yz(mye,1,2)
                  read(*,*) ys,ye
                  do 3021 i=mys,mye-1
                     if (ys.lt..5*(yz(i,1,2)+yz(i+1,1,2))) goto 4021
 3021             continue
                  i=mye
 4021             mys=i
                  do i=mys,mye-1
                     if (ye.lt..5*(yz(i,1,2)+yz(i+1,1,2))) goto 4031
                  end do
                  i=mye
 4031             mye=i
                  write(*,9330) yz(mys,1,2),yz(mye,1,2)
c
                  write(*,9340)yz(1,mzs,1),yz(1,mze,1)
                  read(*,*) zs,ze
                  do i=mzs,mze-1
                     if (zs.lt..5*(yz(1,i,1)+yz(1,i+1,1))) goto 4041
                  end do
                  i=mze
 4041             mzs=i
                  do i=mzs,mze-1
                     if (ze.lt..5*(yz(1,i,1)+yz(1,i+1,1))) goto 4051
                  end do
                  i=mze
 4051             mze=i
                  write(*,9350)yz(1,mzs,1),yz(1,mze,1)
               end if
            end if
         end if

         if (ttmean.and.ipln.eq.1) then
            sxs=1
            szs=1
            sxe=nx+1
            sze=nz+1
            write(*,*) 'Give box around streak'
            write(*,5305)xy(sxs,1,1),xy(sxe,1,1)
 5305       format(' Give start and end for x  (',F8.3,' - ',F8.3,')')
            read(*,*) xss,xe
            do i=sxs,sxe-1
               if (xss.lt..5*(xy(i,1,1)+xy(i+1,1,1))) goto 5300
            end do
            i=sxe
 5300       sxs=i
            do i=sxs,sxe-1
               if (xe.lt..5*(xy(i,1,1)+xy(i+1,1,1))) goto 5210
            end do
            i=sxe
 5210       sxe=i
            write(*,5310) xy(sxs,1,1),xy(sxe,1,1)
 5310       format(' Nearest values  ',F8.3,' - ',F8.3)

            write(*,6340) xz(1,szs,2),xz(1,sze,2)
 6340       format(' Givestart and end for z  (',F8.3,' - ',F8.3,')')
            read(*,*) zs,ze
            do i=szs,sze-1
               if (zs.lt..5*(xz(1,i,2)+xz(1,i+1,2))) goto 6240
            end do
            i=sze
 6240       szs=i
            do i=szs,sze-1
               if (ze.lt..5*(xz(1,i,2)+xz(1,i+1,2))) goto 6250
            end do
            i=sze
 6250       sze=i
            write(*,6350)xz(1,szs,2),xz(1,sze,2)
 6350       format(' Nearest values ',F8.3,' - ',F8.3)

         end if
c
c
c
         if (tpl.eq.1.and.iplot.eq.1) then
            write(heada,9500) var,cpl,t,re
 9500       format(1x,A12,' z= ',F6.2,' t= ',F6.1,' re= ',F7.1,'$')
            call cont6(nx+1,nyp,plane,grid,mxs,mxe,mys,mye,ncont,
     &           blev,elev,'x$','y$',heada,headb,mbox,jhard,lraw,gl,gu,
     &           ipln,ifil,umin,umax,nothing)
         end if
         if (tpl.eq.2.and.iplot.eq.1) then
            if (ttmean) then
               do  z=1,plsz
                  plane(z)=plane(z)-plmean(z)
               end do
               call secinst(plane,grid,esi,sxs,szs,sxe,sze,nx+1,nz+1,
     &              npl,ipln,esigr,t)
            end if
            write(heada,9510) var,cpl,t,re
 9510       format(1x,A12,' y= ',F6.3,' t= ',F6.1,' re= ',F7.1,'$')
            call cont6(nx+1,nz+1,plane,grid,mxs,mxe,mzs,mze,ncont,
     &           blev,elev,'x$','z$',heada,headb,mbox,jhard,lraw,gl,gu,
     &           ipln,ifil,umin,umax,nothing)
         end if
         if (tpl.eq.3.and.iplot.eq.1) then
            if (ans.eq.'p') then
               do k=1,nyp*(nz+1)
                  write(78,*) plane(k),grid(k,1),grid(nyp*(nz+1)+k,1)
               end do
               do k=mys,mye
                  do kk=mzs,mze
                     write(77,*) plane(k+(kk-1)*(nyp)),
     &                    grid(k+(kk-1)*(nyp),1),
     &                    grid(nyp*(nz+1)+k+(kk-1)*(nyp),1)
                  end do
               end do
               write(*,*)'crop', mys,mye,mzs,mze
               close(77)
            end if
            write(heada,9515) var,cpl,t,re

 9515       format(1x,A12,' x= ',F6.2,' t= ',F6.1,' re= ',F7.1,'$')
            call cont6(nyp,nz+1,plane,grid,mys,mye,mzs,mze,ncont,
     &           blev,elev,'z$','y$',heada,headb,mbox,jhard,lraw,gl,gu,
     &           ipln,ifil,umin,umax,nothing)
         end if
         if (tpl.eq.2.and.iplot.eq.2.and..not.lcmean) then
            write(heada,9520) var,cpl,t,re
 9520       format(1x,A12,' y= ',F6.3,' t= ',F6.1,' re= ',F7.1,'$')
            ymin=-0.5
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,nx/2+1,1,nx/2+1,
     &           'deltax$','corr$',heada,headb,' $',mbox,1,jhard)
         end if
         if (tpl.eq.2.and.iplot.eq.3.and..not.lcmean) then
            write(heada,9530) var,cpl,t,re
 9530       format(1x,A12,' y= ',F6.3,' t= ',F6.1,' re= ',F7.1,'$')
            call rita1a(gxax,graf,xmin,xmax,ymin,ymax,nz/2+1,1,nz/2+1,
     &           'deltaz$','corr$',heada,headb,' $',mbox,1,jhard)
         end if
c
c     Next plot
c
      end do
c
c     After all planes have been processed
c
      if (ttmean) then
         write(*,*) 'Energy as function of (x,t)...'
         close(86)
      end if
      if (tpl.eq.2.and.iplot.eq.2.and.lcmean) then
         write(heada,9600) var,cpl,tf,tl,re
 9600    format(1x,A12,' y= ',F6.3,' t= ',F6.1,' - ',F6.1,', re= ',
     &        F7.1,'$')
         ymin=-0.5
         xic=-(grafm(1)+grafm(nx/2+1))*.5
         do iix=1,nx/2+1
            xic=xic+grafm(iix)
         end do
         write(headc,9605) xic*xl/real(nx)
 9605    format(1x,'integral length scale =',F8.4)
         call rita1a(gxax,grafm,xmin,xmax,ymin,ymax,nx/2+1,1,nx/2+1,
     &        'deltax$','corr$',heada,headb,headc,mbox,1,jhard)
      end if
c
      if (tpl.eq.2.and.iplot.eq.3.and.lcmean) then
         write(heada,9610) var,cpl,tf,tl,re
 9610    format(1x,A12,' y= ',F6.3,' t= ',F6.1,' - ',F6.1,', re= ',
     &        F7.1,'$')
         call rita1a(gxax,grafm,xmin,xmax,ymin,ymax,nz/2+1,1,nz/2+1,
     &        'deltaz$','corr$',heada,headb,' $',mbox,1,jhard)
      end if
      if (tpl.eq.1.and.iplot.eq.4) then
         write(heada,9620) var,cpl,tf,tl,re
 9620    format('mean',1x,A12,' z= ',F6.2,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nx+1,nyp,plmean,grid,mxs,mxe,mys,mye,ncont,blev,
     &        elev,'x$','y$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,ifil,ifil,umin,umax,nothing)
      end if
      if (tpl.eq.1.and.iplot.eq.5) then
         write(heada,9630) var,cpl,tf,tl,re
 9630    format('rms',1x,A12,' z= ',F6.2,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nx+1,nyp,plrms,grid,mxs,mxe,mys,mye,ncont,blev,
     &        elev,'x$','y$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,umin,umax,nothing)
      end if
      if (tpl.eq.2.and.iplot.eq.4) then
         write(heada,9640) var,cpl,tf,tl,re
 9640    format('mean',1x,A7,' y= ',F6.3,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nx+1,nz+1,plmean,grid,mxs,mxe,mzs,mze,ncont,blev,
     &        elev,'x$','z$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,umin,umax,nothing)
      end if
      if (tpl.eq.2.and.iplot.eq.5) then
         write(heada,9650) var,cpl,tf,tl,re
 9650    format('rms',1x,A12,' y= ',F6.3,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nx+1,nz+1,plrms,grid,mxs,mxe,mzs,mze,ncont,blev,
     &        elev,'x$','z$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,umin,umax,nothing)
      end if
      if (tpl.eq.3.and.iplot.eq.4) then
         write(heada,9641) var,cpl,tf,tl,re
 9641    format('mean',1x,A7,' x= ',F6.3,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nyp,nz+1,plmean,grid,mys,mye,mzs,mze,ncont,blev,
     &        elev,'x$','z$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,umin,umax,nothing)
      end if
      if (tpl.eq.3.and.iplot.eq.5) then
         write(heada,9651) var,cpl,tf,tl,re
 9651    format('rms',1x,A12,' x= ',F6.3,' t= ',F6.1,' - ',F6.1,
     &        ' re= ',F7.1,'$')
         call cont6(nyp,nz+1,plrms,grid,mys,mye,mzs,mze,ncont,blev,
     &        elev,'x$','z$',heada,headb,mbox,jhard,lraw,gl,gu,ipln,
     &        ifil,umin,umax,nothing)
      end if
      if (iplot.eq.6) then
         if (lsubm) then
            do j=1,npr
               do i=1,nint
                  upr(i,j)=upr(i,j)-uprm(j)
               end do
            end do
         end if
         if (tpl.eq.1) then
            write(heada,9660) var,xpr,xpr+(npr-1)*dxpr,ypr,cpl,re
 9660       format(A12,' x= ',F6.2,' - ',F6.2,' y= ',F6.3,' z= ',
     &           F6.3,' re= ',F7.1,'$')
         end if
         if (tpl.eq.2) then
            write(heada,9660) var,xpr,xpr+(npr-1)*dxpr,cpl,zpr,re
         end if
         if (tpl.eq.3) then
            write(heada,9661) var,zpr,zpr+(npr-1)*dzpr,cpl,ypr,re
 9661       format(A12,' z= ',F6.2,' - ',F6.2,' x= ',F6.3,' y= ',
     &           F6.3,' re= ',F7.1,'$')
         end if
         do i=1,npr
            npoint(i)=nint
         end do
         yaxis=var
         if (npr.gt.1) then
            if (tpl.le.2) then
               call stagpr(upr,npr,npl,nint,dxpr,xpr,xmoff)
            else
               call stagpr(upr,npr,npl,nint,dzpr,zpr,xmoff)
            end if
            write(yaxis,9670) var,xmoff
 9670       format(A12,'+x*',F7.3,'$')
         end if
         call rita1a(tpr,upr,tf,tl,0.,0.,npoint,npr,npl,
     &         't$',yaxis,heada,headb,' $',mbox,0,jhard)
      end if
      if (iplot.eq.7) then
         headc=' Disturbance Length $'
         if (tpl.eq.1) then
            write(heada,9680) var,cpl,re
 9680       format(A12,' z= ',F6.3,' re= ',F7.1,'$')
         else
            write(heada,9690) var,cpl,re
 9690       format(A12,' y= ',F6.3,' re= ',F7.1,'$')
         end if
         mbox=1
         call rita1a(tlen,dlen,0.,0.,0.,0.,nint,1,npl,
     &        't$','length$',headc,heada,headb,mbox,0,jhard)
      end if
      if (iplot.eq.8) then
         headc='Streamwise cross section $'
         if (tpl.eq.1) then
            write(heada,9710) var,tf,tl,ypr,cpl,re
 9710       format(A12,' t= ',F7.3,' - ',F7.3,
     &           ' y= ',F6.3,' z= ',F6.3,' re= ',F7.1,'$')
         else
            write(heada,9710) var,tf,tl,cpl,zpr,re
         end if
         mbox=1
         write(yaxis,9720) var
 9720    format(A12,'$')
         idash=1
         do i=1,nint
            npoint(i)=nx+1
         end do
         call rita1a(xc,fx,0.,0.,0.,0.,npoint,nint,nx+1,
     &        'x$',yaxis,headc,heada,headb,mbox,idash,jhard)
      end if
      if (iplot.eq.10) then
         if (lsubm) then
            do j=1,npr
               do i=1,nint
                  upr(i,j)=upr(i,j)-uprm(j)
               end do
            end do
         end if
         call espec(fesp,esp,upr,npl,npr,tint,nint,tf,tl,
     &        lhanw,nbin,lper)
         if (tpl.eq.1) then
            write(heada,9660) var,xpr,xpr+(npr-1)*dxpr,ypr,cpl,re
         end if
         if (tpl.eq.2) then
            write(heada,9660) var,xpr,xpr+(npr-1)*dxpr,cpl,zpr,re
         end if
         if (tpl.eq.3) then
            write(heada,9661) var,zpr,zpr+(npr-1)*dzpr,cpl,ypr,re
         end if
         do i=1,npr
            npoint(i)=1+(nint/2-nbin)/nbin
            if (logax) npoint(i)=(nint/2-nbin)/nbin
         end do
         yaxis=var
         if (npr.gt.1) then
            if (tpl.le.2) then
               call stagpr(esp,npr,npl/2,nint/2,dxpr,xpr,xmoff)
            else
               call stagpr(esp,npr,npl/2,nint/2,dzpr,zpr,xmoff)
            end if
            write(yaxis,9670) var,xmoff
         end if
         write(63,*) npoint(1)
         do i=2,npoint(1)+1
            write(63,*) fesp(i,1),esp(i,1,1)
         end do
         if (logax) then
            call ritlog(fesp(2,1),esp(2,1,1),
     &           ommin,ommax,emin,emax,npoint,
     &           npr,npl/2,'omega',yaxis,heada,headb,' ',3,0)
         else
            call rita1a(fesp,esp,fesp(1,1),fesp(npoint(1),1),0.,0.,
     &           npoint,npr,npl/2,'omega',yaxis,heada,headb,' ',
     &           mbox,0,jhard)
         end if
      end if

      write(*,*) 'Global min/max ',umin,umax

      goto 1000

      end program rps



      subroutine secinst(plane,grid,esi,sxs,szs,sxe,sze,nx,nz,
     &     npl,ipln,esigr,t)

      implicit none

      integer sxs,sxe,szs,sze,nx,nz,npl
      integer x,z,ipln

      real plane(nx,nz),grid(nx,nz,2)
      real esi(nx,npl),esigr(nx,npl,2)

      real t

      do x=sxs,sxe
         esi(x,ipln)=0.
         esigr(x,ipln,1) = grid(x,szs,1)
         esigr(x,ipln,2) = t
         do z=szs,sze
            esi(x,ipln)=esi(x,ipln)+plane(x,z)**2
c            write(*,*) grid(x,z,1),grid(x,z,2),x,z
         end do
         esi(x,ipln) = esi(x,ipln) / (1.*(sze-szs))
         write(86,'(3e18.9)') esigr(x,ipln,1),esigr(x,ipln,2),
     &        esi(x,ipln)
      end do
      write(86,*)

      end subroutine secinst
