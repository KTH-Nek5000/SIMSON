c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
c
c           PPPPPP    LL        OOOOO   TTTTTTTT     111
c           PP   PP   LL       OO   OO     TT       1111
c           PP   PP   LL       OO   OO     TT         11
c           PPPPPP    LL       OO   OO     TT         11
c           PP        LL       OO   OO     TT         11
c           PP        LLLLLLL   OOOOO      TT        1111
c
c     Plot1, version 1.7, Copyright 1993 Dan S. Henningson
c
c     Plotting package that should run on any UNIX system with Tektronix
c     plotting capability which has a FORTRAN compiler. With the 1.2
c     version a Postscript driver has been included. Below are some
c     simple sample programs that uses most of the user intended
c     routines. Most parametrs should be self explanatory. The terminal
c     types supported so far is tektronix (tek) and xterm (xterm1)
c     in the X environment. There are also a contour plotting routine
c     included that originates at NASA Ames.
c
c     In X environment a hardcopy can be obtained in two ways,
c     1) the command 'xwd | xpr -device ps | lpr -P"printername"' will first
c     change the cursor to a cross, then you point the cursor
c     to the window you want to print and click a mouse button,
c     the window will then be printed on the printer "printername".
c     This gives rather bad resolution.
c     2) By pushing the control key, the middle mouse button and
c     selecting  the COPY option on the xterm menu that appears when
c     you have the cursor in the tektronix window you will
c     create a file in your main directory of the form COPY.......
c     This file contains the Tektronix characters used to produce the
c     plot in the window where you selected COPY. It may be printed out
c     by translating the tektronix file to a postscript file with ps4014
c     command, this can then can be printed out. This gives a very good
c     resolution plot.
c
c     Version 1.2 has a new/better way to obtain hardcopies of the plot.
c     Postscript and Tektronix commands  defining the plot is written
c     to file using the ihard parameter as follows:
c
c     ihard =  0 : Tektronix plot to the screen
c     ihard >  0 : Plot defined either in the tektronix format or the
c                  postscript format is written to file plot.**
c                  where ** are numbers starting at 79 and are incremented
c                  as more plots within the same runs are saved
c     ihard <  0 : Plots are written to unit 79 (usually file fort.79)
c                  with all plots to the same file
c
c     ihard = +-1: Postscript plot saved in file and also sent to screen
c     ihard = +-2: Postscript plot saved in file but NO screen plot
c     ihard = +-3: Tektronix plot saved in file and NO screen plot
c
c     Version 1.3 can obtain essentially arbitrary proportions on the
c     plot area. It is done with the mbox variable in the following way.
c
c     mbox = 0 : Plot proportions given by the input data
c     mbox = 1 : Plot proportions set to maximum (1 : 0.6)
c     mbox > 1 : Plot proportions set to (real(mbox)/1000.) : 0.6
c     mbox < 0 : Plot proportions set to 1 : (real(-mbox)/1000.)
c
c     Version 1.4. Changes include:
c     1. generic function names have been substituted for specific ones
c     2. all character variables with long strings have been made
c        at least 100 characters long
c     3. the opening of files in hunit have been changed to work on alliants
c        as well as the other systems previously supported
c     4. plotting with logarithmic axis is added
c
c     Version 1.5 Changes
c
c     A 'stroke' is inserted after each 500 'draw' in postscript to
c     avoid overflowing the stack in postscript interpreters (printers etc).
c
c     The new handling of string termination does not require the
c     string to be terminated with a '$', rather the last non-blank
c     character is printed. However, to comply with earlier versions
c     and to allow printing of trailing blank, if a '$' character exists
c     in the string its first occurence is taken as the termination.
c     This also means that strings consisting solely of a '$' character
c     works (but does not print).
c
c     Examples (quotation marks are not printed,c is center of line):
c
c     string      prints         prints with v1.4
c                    c                  c
c     'abc$'       'abc'              'abc'
c     'abc  $'    'abc  '            'abc  '
c     '$'            ''             did not work
c     'abc'        'abc'            did not work
c     ''             ''             did not work
c     'abc  '      'abc'            did not work
c
c     ritlog also takes 0 for the ilog argument producing non-log plots
c
c     version 1.6 changes
c
c     fixed bug in ritlog
c
c     version 1.7 changes
c
c     the name of the fortran code is now plot1.f
c     to minimize changes in makefiles
c     when given all zero data the rita1a
c     and rita1 routines now set the interval to -0.0001,.0001
c     fixed slight bug in drawing first line on screen
c
c     Dan Henningson, Cambridge 1993-06-15
c
c
c **********************************************************************
c
c      Below are sample programs that shows how to use the package
c
c      program par1
c      real x(101),y(101)
c      common/hard/ihard
c
c      sample program that plots out a parabola
c      note that the ihard parameter has to be put in common
c      unless a higer level routine like rita1 is used
c
c      ihard=-1
c 1    write(*,*) ' give xmin,xtick,xmax,ymin,ytick,ymax,mbox,idash'
c      read(*,*) xmin,xtick,xmax,ymin,ytick,ymax,mbox,idash
c      do i=1,101
c         x(i)=(i-1)/100.
c         y(i)=x(i)**2
c      end do
c      call xterm1
c      call plarea(xmin,xmax,ymin,ymax,mbox)
c      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
c      call head1('This is a heading,')
c      call head2('this is also a heading')
c      call head3('with this also.')
c      call taxis('xaxis','yaxis')
c      call curve(x,y,101,idash)
c      call plend
c      goto 1
c      end
c
c
c      program par2
c      real x(101,4),y(101,4)
c      integer n(4)
c
c      sample program that plots 4 curves on the same plot
c      and uses the automatic tick mark facility
c
c 1    write(*,*) ' give a,b,c,d,mbox,idash,ihard'
c      read(*,*) a,b,c,d,mbox,idash,ihard
c
c      do i=1,101
c         x(i,1)=(i-1)/100.
c         y(i,1)=x(i,1)
c         x(i,2)=(i-1)/100.
c         y(i,2)=x(i,1)**2
c         x(i,3)=(i-1)/100.
c         y(i,3)=x(i,1)**3
c         x(i,4)=(i-1)/100.
c         y(i,4)=x(i,1)**4
c      end do
c      n(1)=101
c      n(2)=101
c      n(3)=101
c      n(4)=101
c
c      call rita1a(x,y,a,b,c,d,n,4,101,'xaxis','yaxis',
c     &          'test program','hej','svejs',mbox,idash,ihard)
c
c      goto 1
c      end
c
c
c      program tlog
c
c      sample program that produces logplots
c
c      real x(501,4),y(501,4)
c      integer n(4)
c
c 1    write(*,*) ' give a,b,c,d,ilog,idash'
c      read(*,*) a,b,c,d,ilog,idash
c
c      do i=1,501
c         x(i,1)=i
c         y(i,1)=x(i,1)
c         x(i,2)=i
c         y(i,2)=x(i,1)**2
c         x(i,3)=i
c         y(i,3)=x(i,1)**3
c         x(i,4)=i
c         y(i,4)=x(i,1)**4
c      end do
c      n(1)=501
c      n(2)=501
c      n(3)=501
c      n(4)=501
c
c      call ritlog(x,y,a,b,c,d,n,4,501,'xaxis','yaxis',
c     &          'test of log plot','hej','svejs',ilog,idash)
c
c      goto 1
c      end
c
c
c      program conr
c
c      sample program that gives contours of a parabolic cyllinder
c
c      real u1(51,51),xy(51,51,2)
c
c      mbox=1
c      write(*,*) ' give xconst,yconst'
c      read(*,*) xconst,yconst
c      do i=1,51
c         do j=1,51
c            xy(i,j,1)=(i-1)*xconst
c            xy(i,j,2)=(j-1)*yconst
c            u1(i,j)=((i-1)*xconst)**2+((j-1)*yconst)**2-2.
c         end do
c      end do
c      mbox=1
c 1    call cont1(u1,xy,51,51,'xaxis','yaxis','Contour','test',mbox)
c      goto 1
c      end
c
c***********************************************************************
c
c     Higher level plotting routines intended for the user.
c
c***********************************************************************
c
      subroutine rita1(x,y,xmin1,xmax1,ymin1,ymax1,npoint,nplot,ndim,
     &                xaxis,yaxis,heada,headb,headc,mbox,idash)
c
c     If xmax1.eq.xmin1  choose xmax and xmin from data
c     If ymax1.eq.ymin1  choose ymax and ymin from data
c
c     x(ndim,nplot)
c     y(ndim,nplot)
c     nplot            number of plots
c     npoint(i)        array containing number of points in each plot
c     ndim             declared length of first index in x, y array
c     mbox             0: prop. as specified, 1: max plot area
c     idash            dashtype
c     xaxis            name of xaxis
c     heada            first heading
c
c     ihard             0: no hardcopy, just tek screen
c                     +-1: tek screen plus ps hardcopy (unit=79+)
c                     +-2: ps  hardcopy (unit=79+)
c                     +-3: tek hardcopy (unit=79+)
c                      <0: no unit increment for output files
c
      implicit none

      integer i,j,nplot,iterm1,ihard,idash,jdash,mbox,ndim
      real xmin,xmax,xdiff,xtick,xmax1,xmin1
      real ymin,ymax,ydiff,ytick,ymax1,ymin1
      real x(ndim,1),y(ndim,1)
      integer npoint(1)
      character*(*) xaxis,yaxis,heada,headb,headc
      common /hard/ ihard
      ihard=-1
c
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      if (xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do j=1,nplot
            do i=1,npoint(j)
               xmax=max(xmax,x(i,j))
               xmin=min(xmin,x(i,j))
            end do
         end do
         xdiff=xmax-xmin
         if (xmax.ne.0) xmax=xmax+0.1*xdiff
         if (xmin.ne.0) xmin=xmin-0.1*xdiff
         if (xmin.eq.xmax) then
            xmin=xmin-0.0001
            xmax=xmax+0.0001
         end if
      end if
      if (ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do j=1,nplot
            do i=1,npoint(j)
               ymax=max(ymax,y(i,j))
               ymin=min(ymin,y(i,j))
            end do
         end do
         ydiff=ymax-ymin
         if (ymax.ne.0) ymax=ymax+0.1*ydiff
         if (ymin.ne.0) ymin=ymin-0.1*ydiff
         if (ymin.eq.ymax) then
            ymin=ymin-0.0001
            ymax=ymax+0.0001
         end if
      end if
      xtick=0.
      ytick=0.
c      write(*,*) ' termtype: 1=tek, 2=xterm'
c      read(*,*) iterm1
      iterm1=2
 32   if (iterm1.eq.1) call tek
      if (iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      do j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
      end do
      call plend
      return

      end subroutine rita1



      subroutine rita1a(x,y,xmin1,xmax1,ymin1,ymax1,npoint,nplot,ndim,
     &                xaxis,yaxis,heada,headb,headc,mbox,idash,jhard)
c
c     if xmax1.eq.xmin1  choose xmax and xmin from data
c     if ymax1.eq.ymin1  choose ymax and ymin from data
c
c     x(ndim,nplot)
c     y(ndim,nplot)
c     nplot            number of plots
c     npoint(i)        array containing number of points in each plot
c     ndim             declared length of first index in x, y array
c     mbox             0: prop. as specified, 1: max plot area
c     idash            0: no dash, 1: dashed lines
c     xaxis            name of xaxis
c     heada            first heading
c     jhard            hardcopy (only difference with rita1 routine)
c                      see rita1 routine for explanation
c
      implicit none

      integer i,j,nplot,iterm1,ihard,jhard,idash,jdash,mbox,ndim
      real xmin,xmax,xdiff,xtick,xmax1,xmin1
      real ymin,ymax,ydiff,ytick,ymax1,ymin1
      real dy
      real x(ndim,1),y(ndim,1)
      integer npoint(1)
      character*(*) xaxis,yaxis,heada,headb,headc
      common /hard/ ihard
      ihard=jhard
c
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      if (xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do j=1,nplot
            do i=1,npoint(j)
               xmax=max(xmax,x(i,j))
               xmin=min(xmin,x(i,j))
            end do
         end do
         xdiff=xmax-xmin
         if (xmax.ne.0) xmax=xmax+0.1*xdiff
         if (xmin.ne.0) xmin=xmin-0.1*xdiff
         if (xmin.eq.xmax) then
           xmin=xmin-0.0001
           xmax=xmax+0.0001
         end if
      end if
      if (ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do j=1,nplot
            do i=1,npoint(j)
               ymax=max(ymax,y(i,j))
               ymin=min(ymin,y(i,j))
            end do
         end do
         ydiff=ymax-ymin
         if (ymax.ne.0) ymax=ymax+0.1*ydiff
         if (ymin.ne.0) ymin=ymin-0.1*ydiff
         if (ymin.eq.ymax) then
           ymin=ymin-0.0001
           ymax=ymax+0.0001
         end if
      end if
      xtick=0.
      ytick=0.
c      write(*,*) ' termtype: 1=tek, 2=xterm'
c      read(*,*) iterm1
      iterm1=2
 32   if (iterm1.eq.1) call tek
      if (iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      dy=(ymax-ymin)*.01
      do j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
      end do
      call plend
      return

      end subroutine rita1a



      subroutine ritlog(x,y,xmin1,xmax1,ymin1,ymax1,npoint,nplot,ndim,
     &                xaxis,yaxis,heada,headb,headc,ilog,idash)
c
c     If xmax1.eq.xmin1  choose xmax and xmin from data
c     If ymax1.eq.ymin1  choose ymax and ymin from data
c
c     x(ndim,nplot)
c     y(ndim,nplot)
c     nplot            number of plots
c     npoint(i)        array containing number of points in each plot
c     ndim             declared length of first index in x, y array
c     ilog             0: both axes lin, 1: x-axis log, 2: y-axis,
c                      3: both axes log
c     idash            dashtype
c     xaxis            name of xaxis
c     heada            first heading
c
c     ihard             0: no hardcopy, just tek screen
c                     +-1: tek screen plus ps hardcopy (unit=79+)
c                     +-2: ps  hardcopy (unit=79+)
c                     +-3: tek hardcopy (unit=79+)
c                      <0: no unit increment for output files
c
      implicit none

      integer i,j,nplot,iterm1,idash,jdash,ilog,mbox,ihard,ndim
      real xmin,xmax,xmax1,xmin1,xtick
      real ymin,ymax,ymax1,ymin1,ytick
      real x(ndim,1),y(ndim,1)
      integer npoint(1)
      character*(*) xaxis,yaxis,heada,headb,headc
      common /hard/ ihard
      ihard=-1
      mbox=1
c
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      if (xmin.eq.xmax) then
         xmin=1.e20
         xmax=-1.e20
         do j=1,nplot
            do i=1,npoint(j)
               xmax=max(xmax,x(i,j))
               xmin=min(xmin,x(i,j))
            end do
         end do
c         xdiff=xmax-xmin
c         if (xmax.ne.0) xmax=xmax+0.1*xdiff
c         if (xmin.ne.0) xmin=xmin-0.1*xdiff
      end if
      if (ymin.eq.ymax) then
         ymin=1.e20
         ymax=-1.e20
         do j=1,nplot
            do i=1,npoint(j)
               ymax=max(ymax,y(i,j))
               ymin=min(ymin,y(i,j))
            end do
         end do
c         ydiff=ymax-ymin
c         if (ymax.ne.0) ymax=ymax+0.1*ydiff
c         if (ymin.ne.0) ymin=ymin-0.1*ydiff
      end if

      if (ilog.eq.1.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               if (x(i,j).le.0.0) then
                  write(*,*) 'x-array does not contain positive numbers'
                  return
               end if
            end do
         end do
      end if
      if (ilog.eq.2.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               if (y(i,j).le.0.0) then
                  write(*,*) 'y-array does not contain positive numbers'
                  return
               end if
            end do
         end do
      end if

      if (ilog.eq.1.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               x(i,j)=log10(x(i,j))
            end do
         end do
         xmin=log10(xmin)
         xmax=log10(xmax)
      end if
      if (ilog.eq.2.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               y(i,j)=log10(y(i,j))
            end do
         end do
         ymin=log10(ymin)
         ymax=log10(ymax)
      end if
      xtick=0.
      ytick=0.
c      write(*,*) ' termtype: 1=tek, 2=xterm'
c      read(*,*) iterm1
      iterm1=2
 32   if (iterm1.eq.1) call tek
      if (iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call logaxis(xmin,xtick,xmax,ymin,ytick,ymax,ilog)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      do j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
      end do
      call plend
      if (ilog.eq.1.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               x(i,j)=10.**(x(i,j))
            end do
         end do
      end if
      if (ilog.eq.2.or.ilog.eq.3) then
         do j=1,nplot
            do i=1,npoint(j)
               y(i,j)=10.**(y(i,j))
            end do
         end do
      end if
      return

      end subroutine ritlog



      subroutine cont1(u1,xy,nx,ny,xaxis,yaxis,heada,headb,mbox)
c
c     Calculates contours of values in an array u1 with associated
c     coordinates xy
c
      implicit none

      integer ic,is,isnext,idash,nd,n,ncont,ncopt,mbox,iterm,ihard
      integer nx,ny,nxdim,naddim
      real aincr,start,end
      real xmin,xmax
      real ymin,ymax
      parameter(nxdim=50000,naddim=5000)
      real u1(nx,ny),xy(nx,ny,2)
      real acont(naddim),xcont(nxdim),ycont(nxdim)
      integer nlev(naddim),work(nxdim),nad(naddim)
      character text*100
      character*(*) xaxis,yaxis,heada,headb
      common /hard/ ihard
      ihard=-1
c
2     write(*,*) ' # of contours (0=manual)'
      read(*,*) ncont
      ncopt=1
      if (ncont.eq.0) then
          ncopt=2
          write(*,*) ' Give start, end and increment of contours'
          read(*,*) start,end,aincr
          ncont=int((end-start)/aincr+.01)+1
          do n=1,ncont
             acont(n)=start+aincr*real(n-1)
          end do
      end if
c
      call contxx(nx,ny,1,nx,1,ny,xy,u1,ncopt,ncont,acont,
     *         nxdim,xcont,ycont,naddim,nad,nlev,nxdim,work)
c
      if (ncont.eq.1) then
      write(text,1221) ncont,acont(1)
 1221 format(i3,' level at:',f7.3)
      else
      write(text,1222) ncont,acont(1),acont(ncont),acont(2)-acont(1)
 1222 format(i3,' lev., beg:',e11.4,' end:',e11.4,
     &' sp:',e11.4)
      end if
      xmin=xy(1,1,1)
      xmax=xy(nx,ny,1)
      ymin=xy(1,1,2)
      ymax=xy(nx,ny,2)
c      write(*,*) ' termtype (1=tek, 2=xterm)'
c      read(*,*) iterm
      iterm=2
      if (iterm.eq.1) call tek
      if (iterm.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,0.,xmax,ymin,0.,ymax)
      call head1(heada)
      call head2(headb)
      call head3(text)
      call taxis(xaxis,yaxis)
c
      isnext=1
      nd=nad(1)
      do ic=2,nd+1
         is=isnext
         isnext=nad(ic)
         idash=1
         if (acont(nlev(ic-1)).lt.0) idash=2
         call curve(xcont(is),ycont(is),isnext-is,idash)
      end do
      call plend
c      goto 2

      end subroutine cont1



c********************************************************************
c
c     Subroutine included for compatibility with older versions
c     of the plotting package, the same thing can now be done
c     better with rita1 or rita1a using the new version of curve
c
       subroutine rita2(x,y,xmin1,xmax1,ymin1,ymax1,npoint,
     &                xaxis,yaxis,heada,headb,headc,mbox)
c
c     Plots boxes for the points in the coordinate positions
c
c     npoint           number of points in plot
c     mbox             0: prop. as specified, 1: max plot area
c     xaxis            name of xaxis
c     heada            first heading
c
      implicit none

      integer i,iterm1,mbox,npoint,ihard
      real xmin,xmax,xmin1,xmax1,xdiff,xtick
      real ymin,ymax,ymin1,ymax1,ydiff,ytick
      real x(1),y(1)
      character*(*) xaxis,yaxis,heada,headb,headc
      common /hard/ihard
      ihard=-1
c
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      if (xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do i=1,npoint
            xmax=max(xmax,x(i))
            xmin=min(xmin,x(i))
         end do
         xdiff=xmax-xmin
         if (xmax.ne.0) xmax=xmax+0.1*xdiff
         if (xmin.ne.0) xmin=xmin-0.1*xdiff
      end if
      if (ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do i=1,npoint
            ymax=max(ymax,y(i))
            ymin=min(ymin,y(i))
         end do
         ydiff=ymax-ymin
         if (ymax.ne.0) ymax=ymax+0.1*ydiff
         if (ymin.ne.0) ymin=ymin-0.1*ydiff
      end if
       xtick=0.
       ytick=0.
c       write(*,*) ' termtype: 1=tek, 2=xterm'
c       read(*,*) iterm1
      iterm1=2
 32   if (iterm1.eq.1) call tek
      if (iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      do i=1,npoint
         if (y(i).lt.0.5.or.y(i).gt.-1.) call box1(x(i),y(i),0.02)
      end do
      call plend
      return

      end subroutine rita2



      subroutine box1(x,y,al)

      implicit none

      real al,al2,x,y
      al2=al/2.
      call plotr(x+al2,y+al2,0)
      call plotr(x+al2,y-al2,1)
      call plotr(x-al2,y-al2,1)
      call plotr(x-al2,y+al2,1)
      call plotr(x+al2,y+al2,1)
      return

      end subroutine box1



c*******************************************************************
c
c     Intermeditate level routines for making custom plots.
c
c*******************************************************************

      subroutine curve(xr,yr,n,idash)
c
c     Draws a curve with different dashes. idash=1: solid, 2: dash
c     3: dot, 4: chain dash, 5: chain dot, 6: long dash.
c
      implicit none

      integer i,i1,n,isymb,idash
      real h1,hd1,hs1,hd2,hs2
      real x1,x2,xor,xsize,xscale,xoff
      real y1,y2,yor,ysize,yscale,yoff

      real xr(1),yr(1)
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      if (idash.ge.0) then
         call setdash(idash,hd1,hs1,hd2,hs2,h1)
         i1=0
         do i=1,n-1
            x1=xr(i)*xscale+xoff
            y1=yr(i)*yscale+yoff
            x2=xr(i+1)*xscale+xoff
            y2=yr(i+1)*yscale+yoff
            call dash(h1,i1,x1,y1,x2,y2,hd1,hs1,hd2,hs2)
            if (idash.eq.1.and.mod(i,500).eq.0) call plota(x2,y2,0)
         end do
      else
         isymb=abs(idash)
         do i=1,n
            x1=xr(i)*xscale+xoff
            y1=yr(i)*yscale+yoff
            call boxes(x1,y1,0.014,isymb)
         end do
      end if
      return

      end subroutine curve



      subroutine plarea(xmin1,xmax1,ymin1,ymax1,mbox1)
c
c     Defines plotting area and user coordinates
c     mbox=1: maximum ploting area, 2: user specified proportions
c
      implicit none

      integer mbox,mbox1
      real xor,xsize,xscale,xoff,xmax,xmin
      real yor,ysize,yscale,yoff,ymax,ymin
      real xmin1,xmax1
      real ymin1,ymax1
      real sizeb,sizeq
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call area(0.,1.313,0.,1.)
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      mbox=mbox1
cc      xor=0.25
      yor=0.15
      xsize=1.0
      ysize=0.6
      if (mbox.gt.1) xsize=real(mbox)/1000.
      if (mbox.lt.0) ysize=real(-mbox)/1000.
      if (mbox.eq.0) then
         sizeq=xsize/ysize
         sizeb=(xmax-xmin)/(ymax-ymin)
         if (sizeb.ge.sizeq) ysize=ysize*sizeq/sizeb
         if (sizeb.lt.sizeq) xsize=xsize*sizeb/sizeq
      end if
      xor=.75-.5*xsize
      xscale=xsize/(xmax-xmin)
      yscale=ysize/(ymax-ymin)
      xoff=xor-xmin*xscale
      yoff=yor-ymin*yscale
      return

      end subroutine plarea



      subroutine across(xmin,xmax,ymin,ymax)
c
c     Draw axis cross
c
      implicit none

      real xmax,xmin
      real ymax,ymin
      if (xmin.lt.0.0.and.xmax.gt.0.0) then
         call plotr(0.,ymin,0)
         call plotr(0.,ymax,1)
      end if
      if (ymin.lt.0.0.and.ymax.gt.0.0) then
         call plotr(xmin,0.,0)
         call plotr(xmax,0.,1)
      end if
      return

      end subroutine across



      subroutine axis(xmin1,xtick1,xmax1,ymin1,ytick1,ymax1)
c
c     Draw axis with tickmarks at xtick and ytick intervals
c     If xtick or ytick are 0. automatic calculation of nice tick marks
c     are done
c
      implicit none

      integer i,istart,iend,itick
      real x,xa,xor,xsize,xscale,xoff,xtick,charx,xmax,xmin
      real y,ya,yor,ysize,yscale,yoff,ytick,chary,ymax,ymin
      real xmin1,xmax1,xtick1
      real ymin1,ymax1,ytick1
      character*100 num
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xmin=xmin1
      xtick=xtick1
      xmax=xmax1
      ymin=ymin1
      ytick=ytick1
      ymax=ymax1
      call font(2,charx,chary)
      if (xtick.eq.0.) then
         itick=9
         if (xsize.le.0.5) itick=5
         call  conscl(xmin,xmax,itick,xtick)
      end if
      if (ytick.eq.0.) call conscl(ymin,ymax,9,ytick)
      call frame(xmin,xmax,ymin,ymax)
      call across(xmin,xmax,ymin,ymax)
      istart=xmin/xtick
      if (xmin.gt.0.0) istart=istart+1
      iend=xmax/xtick
      if (xmax.lt.0.0) iend=iend-1
      do i=istart,iend
         x=i*xtick
         call plotr(x,ymin,0)
         call plotr(x,ymin+charx/yscale,1)
         call flstr(x,num,xtick)
         call xxa(x,xa)
         call writ(num,xa,yor,charx,chary,1,-1.5,-.5)
         call plotr(x,ymax,0)
         call plotr(x,ymax-charx/yscale,1)
      end do
      istart=ymin/ytick
      if (ymin.gt.0.0) istart=istart+1
      iend=ymax/ytick
      if (ymax.lt.0.0) iend=iend-1
      do i=istart,iend
         y=i*ytick
         call plotr(xmin,y,0)
         call plotr(xmin+charx/xscale,y,1)
         call flstr(y,num,ytick)
         call yya(y,ya)
         call writ(num,xor,ya,charx,chary,1,-.5,-1.2)
         call plotr(xmax,y,0)
         call plotr(xmax-charx/xscale,y,1)
      end do
      return

      end subroutine axis



      subroutine logaxis(xmin1,xtick1,xmax1,ymin1,ytick1,ymax1,ilog)
c
c     Draw logaxis with tickmarks at xtick and ytick intervals
c     If xtick or ytick are 0. automatic calculation of nice tick marks
c     are done, for logaxis tickmarks are also generated automatically
c
      implicit none

      integer i,ilog,istart,iend,itick
      real x,xa,xor,xsize,xscale,xoff,xtick,charx,xmax,xmin,xexp,x1
      real y,ya,yor,ysize,yscale,yoff,ytick,chary,ymax,ymin,yexp,y1
      real xmin1,xmax1,xtick1
      real ymin1,ymax1,ytick1
      character*100 num
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xmin=xmin1
      xtick=xtick1
      xmax=xmax1
      ymin=ymin1
      ytick=ytick1
      ymax=ymax1
      call font(2,charx,chary)
      if (xtick.eq.0..and.ilog.ne.1.and.ilog.ne.3) then
         itick=9
         if (xsize.le.0.5) itick=5
         call  conscl(xmin,xmax,itick,xtick)
      end if
      if (ytick.eq.0..and.ilog.ne.2.and.ilog.ne.3)
     &     call conscl(ymin,ymax,9,ytick)
      call frame(xmin,xmax,ymin,ymax)
c      call across(xmin,xmax,ymin,ymax)

      if (ilog.ne.1.and.ilog.ne.3) then
         istart=xmin/xtick
         if (xmin.gt.0.0) istart=istart+1
         iend=xmax/xtick
         if (xmax.lt.0.0) iend=iend-1
         do i=istart,iend
            x=i*xtick
            call plotr(x,ymin,0)
            call plotr(x,ymin+charx/yscale,1)
            call flstr(x,num,xtick)
            call xxa(x,xa)
            call writ(num,xa,yor,charx,chary,1,-1.5,-.5)
            call plotr(x,ymax,0)
            call plotr(x,ymax-charx/yscale,1)
         end do
      else
         xtick=100000000.
         if (xmin.gt.0.0) istart=int(xmin)+1
         if (xmin.le.0.0) istart=int(xmin)
         if (xmax.gt.0.0) iend=int(xmax)
         if (xmax.lt.0.0) iend=int(xmax)-1

         if (iend-istart.gt.10) goto 110
         do 111 itick=1,4
            x1=istart-1+log10(2.*itick)
            if (x1.lt.xmin) goto 111
            call plotr(x1,ymin,0)
            call plotr(x1,ymin+0.5*charx/yscale,1)
            call plotr(x1,ymax,0)
            call plotr(x1,ymax-0.5*charx/yscale,1)
 111     continue

 110     continue
         do 12 i=istart,iend
            x=i
            xexp=10.**i
            call plotr(x,ymin,0)
            call plotr(x,ymin+charx/yscale,1)
            call flstr(xexp,num,xtick)
            call xxa(x,xa)
            call writ(num,xa,yor,charx,chary,1,-1.5,-.5)
            call plotr(x,ymax,0)
            call plotr(x,ymax-charx/yscale,1)

            if (iend-istart.gt.10) goto 12
            do 112 itick=1,4
               x1=x+log10(2.*itick)
               if (x1.gt.xmax) goto 112
               call plotr(x1,ymin,0)
               call plotr(x1,ymin+0.5*charx/yscale,1)
               call plotr(x1,ymax,0)
               call plotr(x1,ymax-0.5*charx/yscale,1)
 112        continue
 12      continue
      end if

      if (ilog.ne.2.and.ilog.ne.3) then
         istart=ymin/ytick
         if (ymin.gt.0.0) istart=istart+1
         iend=ymax/ytick
         if (ymax.lt.0.0) iend=iend-1
         do i=istart,iend
            y=i*ytick
            call plotr(xmin,y,0)
            call plotr(xmin+charx/xscale,y,1)
            call flstr(y,num,ytick)
            call yya(y,ya)
            call writ(num,xor,ya,charx,chary,1,-.5,-1.2)
            call plotr(xmax,y,0)
            call plotr(xmax-charx/xscale,y,1)
         end do
      else
         ytick=100000000.
         if (ymin.gt.0.0) istart=int(ymin)+1
         if (ymin.le.0.0) istart=int(ymin)
         if (ymax.gt.0.0) iend=int(ymax)
         if (ymax.lt.0.0) iend=int(ymax)-1

         if (iend-istart.gt.10) goto 120
         do 121 itick=1,4
            y1=istart-1+log10(2.*itick)
            if (y1.lt.ymin) goto 121
            call plotr(xmin,y1,0)
            call plotr(xmin+0.5*charx/xscale,y1,1)
            call plotr(xmax,y1,0)
            call plotr(xmax-0.5*charx/xscale,y1,1)
 121     continue

 120     continue
         do 22 i=istart,iend
            y=i
            yexp=10.**i
            call plotr(xmin,y,0)
            call plotr(xmin+charx/xscale,y,1)
            call flstr(yexp,num,ytick)
            call yya(y,ya)
            call writ(num,xor,ya,charx,chary,1,-.5,-1.2)
            call plotr(xmax,y,0)
            call plotr(xmax-charx/xscale,y,1)

            if (iend-istart.gt.10) goto 22
            do 122 itick=1,4
               y1=y+log10(2.*itick)
               if (y1.gt.ymax) goto 122
               call plotr(xmin,y1,0)
               call plotr(xmin+0.5*charx/xscale,y1,1)
               call plotr(xmax,y1,0)
               call plotr(xmax-0.5*charx/xscale,y1,1)
 122        continue
 22      continue
      end if
      return

      end subroutine logaxis



      subroutine plotr(xr,yr,ipen)
c
c     ipen=0: move to user coordinate (x,y)
c     ipen=1: draw line to user coordinate (x,y)
c
      implicit none

      integer ipen
      real x,xr,xor,xsize,xscale,xoff
      real y,yr,yor,ysize,yscale,yoff
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      x=xr*xscale+xoff
      y=yr*yscale+yoff
      call plota(x,y,ipen)
      return

      end subroutine plotr



      subroutine xxa(x,xa)
c
c     Convert from user x to absolute x
c
      implicit none

      real xor,xsize,xscale,xoff
      real yor,ysize,yscale,yoff
      real x,xa
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xa=x*xscale+xoff
      return

      end subroutine xxa



      subroutine yya(y,ya)
c
c     Convert from user y to absolute ya
c
      implicit none

      real xor,xsize,xscale,xoff
      real yor,ysize,yscale,yoff
      real y,ya
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      ya=y*yscale+yoff
      return

      end subroutine yya



      subroutine conscl(amin,amax,ncont,ainc)
c
c     Come up with a "nice" scaling of about ncont values between
c     amin and amax.
c
      implicit none

      integer i,nneed,ncont,imin,imax,nnice,ichar
      real amin,amax,ainc,rnice,rmant,char,diff,temp
      logical rev
      dimension rnice(3)
      data rnice/.1,.2,.5/
      data nnice/3/
c
c     Reverse amin,amax if necessary (change back before return)
c
      rev=amax.lt.amin
      if (rev) then
         temp=amin
         amin=amax
         amax=temp
      end if
c
c     As a first approximation, get the difference, its characteristic and
c     mantissa.
c
      diff= (amax-amin)/(ncont+1)
      if (diff.le.0.) goto 20
      char= log10(diff)+1.
c
c     Round char down and get the mantissa.
c
      if (char.ge.0.) then
         ichar= char
      else
         ichar= char-1.
      end if
      rmant= diff*10.**(-ichar)
c
c     What's the next largest "nice" mantissa?
c
      do i= 1,nnice
         if (rmant.le.rnice(i)) goto 10
      end do
      i= nnice
c
c     Got a guess.  calculate a diff, round amin down.
c
 10   continue
      ainc= rnice(i)*10.**ichar
      imin= amin/ainc
      if (amin.gt.imin*ainc) imin= imin+1
      imax= amax/ainc
      if (amax.lt.imax*ainc) imax= imax-1
      nneed= imax+1-imin
c
c     Are we under?
c
      if (nneed.gt.ncont) then
c
c     Nope. Try the next nice number.
c
         if (i.lt.nnice) then
            i= i+1
         else
            ichar= ichar+1
            i= 1
         end if
         goto 10
      end if
c
 20   continue
c
c     If amin,amax reversed change back and change sign of ainc
c
      if (rev) then
         temp=amin
         amin=amax
         amax=temp
         ainc=-ainc
      end if
      return

      end subroutine conscl



c******************************************************************
c
c     Text and string handling routines
c
c******************************************************************

      subroutine flstr(x,num,xtick)
c
c     Read the floating point number x to string num
c     the size of xtick determines the number of decimal places
c
      implicit none

      integer ilen,idec
      real x,xtick,sci
      character*(*) num
      character*100 form
c
      if (x.eq.0.) then
         ilen=3
      else
         ilen=int(log10(abs(x)))
         if (ilen.le.0) then
            ilen=3
         else
            ilen=ilen+3
         end if
      end if
c
      sci=log10(xtick*1.001)
      idec=-int(sci)
c      idec=-int(log10(xtick*1.001))
      if (idec.lt.0) then
         idec=0
      else if (idec.eq.0.and.xtick.ge.1.) then
         idec=0
      else
         idec=idec+1
      end if
c
      if (abs(sci).le.4) then
         ilen=ilen+idec
         if (ilen.lt.10) write(form,97) ilen,idec
         if (ilen.ge.10.and.idec.lt.10) write(form,98) ilen,idec
         if (idec.ge.10) write(form,99) ilen,idec
 97      format(' (f',i1,'.',i1,') ')
 98      format(' (f',i2,'.',i1,') ')
 99      format(' (f',i2,'.',i2,') ')
      else
         form='(1pe8.1,a1)'
      end if
c      write(*,*) x
      write(num,form) x
      return

      end subroutine flstr



      subroutine head1(text)
c
c     Top heading
c
      implicit none

      real xor,xsize,xscale,xoff,charx
      real yor,ysize,yscale,yoff,chary
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(1,charx,chary)
      call writ(text,xor+xsize/2,0.9,charx,chary,1,-.5,-.5)
      return

      end subroutine head1



      subroutine head3(text)
c
c     Bottom heading
c
      implicit none

      real xor,xsize,xscale,xoff,charx
      real yor,ysize,yscale,yoff,chary
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(text,xor+xsize/2,0.8,charx,chary,1,-.5,-.5)
      return

      end subroutine head3



      subroutine head2(text)
c
c     Middle heading
c
      implicit none

      real xor,xsize,xscale,xoff,charx
      real yor,ysize,yscale,yoff,chary
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(text,xor+xsize/2,0.85,charx,chary,1,-.5,-.5)
      return

      end subroutine head2



      subroutine taxis(xtext,ytext)
c
c     Text on x- and y-axis
c
      implicit none

      real xor,xsize,xscale,xoff,charx
      real yor,ysize,yscale,yoff,chary
      character*(*) xtext,ytext
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(xtext,xor+xsize/2,0.075,charx,chary,1,-.5,-.5)
      call writ(ytext,xor-0.15,yor+ysize/2,charx,chary,2,-.5,-.5)
      return

      end subroutine taxis

c*****************************************************************
c
c     Line type and symbol routines
c
c*****************************************************************

      subroutine dash(h1,i,x1,y1,x2,y2,hd1,hs1,hd2,hs2)
c
c      h1 - distance left of dash/space
c      i  : i < 0 curr. space ; i > 0 curr. dash ; i=0 start
c      hdi/hsi  - length of dashes/spaces
c
      implicit none

      integer i,ipen
      real h1,hd1,hd2,hs1,hs2
      real x,y,xm,ym,x1,x2,y1,y2
      real l
      if (i.eq.0) then
          i=1
          call plota(x1,y1,0)
      end if
      x=x1
      y=y1
 1    l=sqrt((x-x2)**2+(y-y2)**2)
      if (h1.gt.l) then
         ipen=1
         if (i.lt.0) ipen=0
         call plota(x2,y2,ipen)
         ipen=1
         h1=h1-l
      else
         xm=x+h1*(x2-x)/l
         ym=y+h1*(y2-y)/l
         ipen=1
         if (i.lt.0) ipen=0
         call plota(xm,ym,ipen)
         ipen=1
         if (i.eq.1) then
            i=-1
            h1=hs1
         else if (i.eq.-1) then
            i=2
            h1=hd2
         else if (i.eq.2) then
            i=-2
            h1=hs2
         else if (i.eq.-2) then
            i=1
            h1=hd1
         end if
         x=xm
         y=ym
         goto 1
      end if
      return

      end subroutine dash



      subroutine setdash(j,hd1,hs1,hd2,hs2,h1)

      implicit none

      integer j
      real h1,hd1,hs1,hd2,hs2,dmmtek
      dmmtek=0.005714
      if (j.le.1) then
          hd1=dmmtek*1000000.0
          hs1=dmmtek*1.0
          hd2=dmmtek*1.0
          hs2=dmmtek*1.0
      end if
      if (j.eq.2) then
          hd1=dmmtek*1.5
          hs1=dmmtek*1.0
          hd2=dmmtek*1.5
          hs2=dmmtek*1.0
      end if
      if (j.eq.3) then
          hd1=dmmtek*0.3
          hs1=dmmtek*1.25
          hd2=dmmtek*0.3
          hs2=dmmtek*1.25
      end if
      if (j.eq.4) then
          hd1=dmmtek*9.0
          hs1=dmmtek*2.5
          hd2=dmmtek*2.5
          hs2=dmmtek*2.5
      end if
      if (j.eq.5) then
          hd1=dmmtek*7.0
          hs1=dmmtek*2.5
          hd2=dmmtek*0.3
          hs2=dmmtek*2.5
      end if
      if (j.ge.6) then
          hd1=dmmtek*10.0
          hs1=dmmtek*5.0
          hd2=dmmtek*10.0
          hs2=dmmtek*5.0
      end if
      h1=hd1
      return

      end subroutine setdash



      subroutine boxes(x,y,aa,i)

      implicit none

      integer i
      real a,aa,x,y
      if (i.eq.1) then
         a=aa/2.
         call plota(x+a,y+a,0)
         call plota(x+a,y-a,1)
         call plota(x-a,y-a,1)
         call plota(x-a,y+a,1)
         call plota(x+a,y+a,1)
      end if
      if (i.eq.2) then
         a=aa/2.
         call plota(x+a,y+a,0)
         call plota(x,y-a,1)
         call plota(x-a,y+a,1)
         call plota(x+a,y+a,1)
      end if
      if (i.ge.3) then
         a=aa/1.414
         call plota(x,y+a,0)
         call plota(x+a,y,1)
         call plota(x,y-a,1)
         call plota(x-a,y,1)
         call plota(x,y+a,1)
      end if
      return

      end subroutine boxes



      subroutine arrow(xr1,yr1,xr2,yr2,d)
c
c     Draws an arrow between point (xr1,yr1) and (xr2,yr2)
c     coordinates are assumed to be user defined and note that
c     we are back at (xr2,yr2) when the routine exits
c
      implicit none

      real d,dx,dy,d1,xr1,xr2,yr1,yr2
      call plotr(xr1,yr1,0)
      call plotr(xr2,yr2,1)
      dx=xr2-xr1
      dy=yr2-yr1
      d1=sqrt(dx*dx+dy*dy)
      dx=dx/d1
      dy=dy/d1
      call plotr(xr2+(-dy-dx)*d,yr2+(dx-dy)*d,1)
      call plotr(xr2,yr2,1)
      call plotr(xr2+(dy-dx)*d,yr2+(-dx-dy)*d,1)
      call plotr(xr2,yr2,1)
      return

      end subroutine arrow



c*****************************************************************
c
c     Driver dependent routines
c
c*****************************************************************

      subroutine plota(x,y,ipen)
c
c     ipen=0: move to absolute coordinate (x,y)
c     ipen=1: draw line to absolute coordinate (x,y)
c
      implicit none

      integer ipen
      real x,y
      call plott(x,y,ipen)
      call ploth(x,y,ipen)
      return

      end subroutine plota



      subroutine frame(xmin,xmax,ymin,ymax)
c
c     Draw frame of plot, make it thicker for postscript
c
      implicit none

      real xmin,xmax,ymin,ymax
      call linew(2.)
c
      call plotr(xmin,ymin,0)
      call plotr(xmax,ymin,1)
      call plotr(xmax,ymax,1)
      call plotr(xmin,ymax,1)
      call plotr(xmin,ymin,1)
c
      call closep
      call linew(0.4)
c
      return
      end subroutine frame



      subroutine writ(text,x,y,charx,chary,irot,yjust,xjust)
c
c     Write text to plot, note differences between tektronix and postscript
c     text                text string;
c                         the first occurence of a $ terminates the
c                         string. alternatively the last non-blank character
c                         is taken as the length of the string.
c     x                   absolute x coordinate
c     y                   absolute y coordinate
c     charx               width of x character (from font, only used in tek)
c     chary               height of y charactrer (measured in absolute units)
c     irot                1: horizontal text, 2: vertical text
c     yjust               0. -> -1., bottom to to top justification
c     xjust               0. -> -1., left to right justification
c
      implicit none

      integer i,n,irot,iunit
      real x,y,xjust,yjust,charx,chary
      character*(*) text
      character*100 name
      common/unit/iunit
      do i=1,100
         name(i:i)=' '
      end do
      n=index(text,'$')
      if (n.eq.0) then
c
c     Find the last non-blank character
c
         do n=len(text),1,-1
            if (text(n:n).ne.' ') goto 1000
         end do
         n=0
 1000    continue
         n=n+1
      end if
      if (iunit.eq.0) goto 1111
      if (irot.eq.1) then
         call plott(x+n*charx*xjust,y+chary*yjust,0)
         call tmode
         write(iunit,111) text(1:n-1)
 111     format(a)
         call gmode
      end if
      if (irot.eq.2) then
         call tmode
         do i=1,n-1
            call plott(x+charx*yjust,y-n*chary*xjust-(i-1)*chary,0)
            call tmode
            write(iunit,111) text(i:i)
            call gmode
         end do
      end if
 1111 continue
c
      call ploth(x,y,0)
      name='('//text(1:n-1)//') p'
      call  writp(name,irot,yjust,xjust)
c
      return

      end subroutine writ



      subroutine font(ifont,charx,chary)
c
c     4 default font sizes, sizes in absolute units for use in tektronix mode
c
      implicit none

      integer ifont
      real charx,chary
      character*1 c(2)
      c(1)=char(27)
      if (ifont.eq.1) then
         c(2)='8'
         charx=0.0184
         chary=charx*1.4
      end if
      if (ifont.eq.2) then
         c(2)='9'
         charx=0.0151
         chary=charx*1.4
      end if
      if (ifont.eq.3) then
         c(2)=':'
         charx=0.0124
         chary=charx*1.4
      end if
      if (ifont.eq.4) then
         c(2)=';'
         charx=0.0101
         chary=charx*1.4
      end if
      call pointqq(c,2)
c
      call pfont(ifont)
c
      return

      end subroutine font



      subroutine xterm1
c
c     Initializes the tektronix plotting, clears the screen and
c     moves to the tektronix window associated with the xterm
c     window.
c
      implicit none

      character*1 c(6)
      integer iterm
      common/term/iterm
      call hunit
      iterm=1
      c(1)=char(27)
      c(2)='['
      c(3)='?'
      c(4)='3'
      c(5)='8'
      c(6)='h'
      call pointqq(c,6)
      call terminal
      call erase
      call pinit
      return

      end subroutine xterm1



      subroutine tek
c
c     Same as xterm1 except that it does not assume X environment
c
      implicit none

      integer iterm
      common/term/iterm
      call hunit
      iterm=0
      call terminal
      call erase
      call pinit
      return
      end subroutine tek



      subroutine hunit
c
c     Initializes tek and ps units
c     iunit=0 means no tek output
c     kunit=0 means no ps  output
c
      implicit none

      integer ihard,iunit,kunit,munit,ierr
      character*100 name,form
      common/hard/ihard
      common/unit/iunit
      common/punit/kunit
c
      if (ihard.eq.0) then
         iunit=6
         kunit=0
      end if
c
      munit=10
      if (ihard.gt.0) then
c
c         ierr=0
c 1000    if (ierr.gt.0) munit=munit+1
c         if (munit.lt.100) form=' (a5,i2) '
c         if (munit.ge.100) form=' (a5,i3) '
c         write(name,form) 'plot.',munit
c         open(unit=munit,file=name,status='new',err=1000,iostat=ierr)
c
c     New way of opening the plot files (should work on alliant)
c
         ierr=1
 1000    if (ierr.eq.0) munit=munit+1
         if (munit.lt.100) form=' (a5,i2) '
         if (munit.ge.100) form=' (a5,i3) '
         write(name,form) 'plot.',munit
         open(unit=munit,file=name,status='old',err=1010,iostat=ierr)
         close(unit=munit)
         goto 1000
 1010    open(unit=munit,file=name,status='new')
      end if
c
      if (abs(ihard).eq.1) then
         iunit=6
         kunit=munit
      end if
      if (abs(ihard).eq.2) then
         iunit=0
         kunit=munit
      end if
      if (abs(ihard).eq.3) then
         iunit=munit
         kunit=0
      end if
c
      end subroutine hunit



      subroutine plend
c
c     Sets mode back to text and in case of X environment moves
c     back to the xterm window
c
      implicit none

      integer iterm,ihard,iunit,kunit,munit
      character*1 c(2)
      common/term/iterm
      common/hard/ihard
      common/unit/iunit
      common/punit/kunit
      munit=max(iunit,kunit)
      call tmode
      c(1)=char(27)
      c(2)=char(3)
      if (iterm.eq.1) call pointqq(c,2)
      call pend
      if (ihard.gt.0.and.munit.ne.6) close(munit)
      return

      end subroutine plend

c*****************************************************************
c
c     Postscript driver routines
c
c*****************************************************************

      subroutine pinit
c
c     Initialize the postscript file
c     The /p routine takes four arguments in the following manner:
c     irot  yjust  xjust  string  p
c
      integer kunit
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      write(kunit,111)
      write(kunit,112)
 111  format('%!',/,
     &   '% Plot1, version 1.4, Copyright 1991, Dan S. Henningson',/,
     &   'gsave',/,
     &   '/d {lineto} bind def',/,
     &   '/m {moveto} bind def',/,
     &   '/s {stroke} bind def',/,
     &   '/w {setlinewidth} bind def',/,
     &   '/c {closepath} bind def',/,
     &   '/sh {show} bind def')
 112  format('/fs 14 def',/,
     &   '/ff {/Times-Roman findfont fs scalefont setfont} bind def',/,
     &   '/ffs {/Symbol findfont fs scalefont setfont} bind def',/,
     &   '/p {currentpoint stroke moveto gsave dup 5 1 roll',/,
     &   'stringwidth pop mul exch fs mul 3 -1 roll 1 sub 90 mul',/,
     &   'rotate rmoveto show currentpoint',/,
     &   'grestore moveto} bind def',/,
     &   'ff',/,
     &   '.4 w')
 1111 continue
      return
c
c     Put in after 'gsave' line in format 111 for side plot:
c
c     &   '600 0 translate',/,
c     &   '90 rotate',/,
c
c     Also change scale in ploth routine,
c     remove multiplication with 0.76923077
c
      end subroutine pinit



      subroutine pend
c
c     Finalize the postscript file
c
      implicit none

      integer kunit
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      write(kunit,*) 's showpage grestore'
 1111 continue
      return
      end

      subroutine closep
c
c     Postscript close path routine for use in frame
c
      implicit none

      integer kunit
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      write(kunit,*) 'c'
 1111 continue
      return

      end subroutine closep



      subroutine pfont(ifont)
c
c     The 4 postdcript default font sizes
c
      implicit none

      integer kunit,ifont
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      if (ifont.eq.1) write(kunit,*) '/fs 16 def'
      if (ifont.eq.2) write(kunit,*) '/fs 14 def'
      if (ifont.eq.3) write(kunit,*) '/fs 12 def'
      if (ifont.eq.4) write(kunit,*) '/fs 10 def'
      write(kunit,*) 'ff'
 1111 continue
      return

      end subroutine pfont



      subroutine writp(name,irot,yjust,xjust)
c
c     Postscript print text
c
      implicit none

      character*(*) name
      integer kunit,n,irot
      real xjust,yjust
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      n=index(name,') p')
      write(kunit,111) irot,yjust,xjust,name(1:n+2)
 111  format(i2,' ',f5.2,' ',f5.2,' ',100a)
 1111 continue
      return

      end subroutine writp



      subroutine ploth(x,y,ipen)
c
c     Write the line drawing commands to the postscript file
c
      implicit none

      integer ipen,kunit
      real x,xp,yp,y
      common/punit/kunit
      if (kunit.eq.0) goto 1111
c
      if (x.gt.1.5.or.x.lt.-0.2) then
         write(kunit,110)
         goto 113
      end if
      if (y.gt.1.2.or.y.lt.-0.2) then
         write(kunit,110)
         goto 113
      end if
 110  format('s ')
c
      xp=600.*x*0.76923077
      yp=600.*y*0.76923077
      if (ipen.eq.0) then
         write(kunit,111) xp,yp
 111     format('s',' ',f6.1,' ',f6.1,' m')
      else
         write(kunit,112) xp,yp
 112     format(f6.1,' ',f6.1,' d')
      end if
 113  continue
 1111 continue
      return

      end subroutine ploth



      subroutine linew(width)
c
c     Change the postscript line width (1 = 1/72 inch)
c
      implicit none

      integer kunit
      real width
      common/punit/kunit
      if (kunit.eq.0) goto 1111
      write(kunit,111) width
 111  format('s',' ',f6.3,' w')
 1111 continue
      return
      end subroutine linew


c*****************************************************************
c
c     Tektronix driver routines
c
c     Low level routines for ploting on a tektronics terminal.
c     The routines are modifications of a those given to by
c     L. Polvani and G. Flierl, MIT.
c
c******************************************************************

      subroutine terminal

      implicit none

      character*1 leadin
      leadin=char(29)
      call pointqq(leadin,1)
      return
      end subroutine terminal



      subroutine gmode

      implicit none

      character*1 leadin
      leadin=char(29)
      call pointqq(leadin,1)
      return
      end subroutine gmode



      subroutine area(xa,xb,ya,yb)

      implicit none

      real x0,x1,xa,xb,y0,y1,ya,yb
      integer hx0,hy0,lx0,ly0,llx0,lly0
      common/penplot/x0,x1,y0,y1,hy0,ly0,llx0,lly0,hx0,lx0
      x0=xa
      y0=ya
      x1=xb
      y1=yb
      hx0=-1
      hy0=-1
      lx0=-1
      ly0=-1
      llx0=-1
      lly0=-1
      return

      end subroutine area



      subroutine tmode

      implicit none

      character*1 leadout
      leadout=char(31)
      call pointqq(leadout,1)
c      call flush(6)
      return

      end subroutine tmode



      subroutine penup

      implicit none

      character*1 pencom
      pencom=char(29)
      call pointqq(pencom,1)
      return

      end subroutine penup



      subroutine erase

      implicit none

      character*1 ers(4)
      ers(1)=char(31)
      ers(2)=char(27)
      ers(3)=char(12)
      ers(4)=char(29)
      call pointqq(ers,4)
      return

      end subroutine erase



      subroutine plott(x,y,ipen)
c
c     Draws line to absolute coordinate (x,y)
c     Moves to absolute coordinate (x,y) if penup
c     was called just prior to plot call (the next plot call
c     will draw line). Absolute coordinates are defined as:
c     0.00 < x < 1.313, 0.00 < y < 1.00, with the maximum plotting area
c     0.25 < x < 1.250, 0.15 < y < 0.75.
c
      implicit none

      character*1 coord(6)
      integer i,ix,iy,ipen,llx,lly
      real x,x0,x1,y,y0,y1
      integer hx0,hy0,lx0,ly0,llx0,lly0
      integer hx,hy,lx,ly
      common/penplot/x0,x1,y0,y1,hy0,ly0,llx0,lly0,hx0,lx0
      ix=(x-x0)/(x1-x0)*4095+.4999999
      iy=(y-y0)/(y1-y0)*3119+.4999999
      if (ix.lt.0 .or. ix.gt.4095 .or. iy.lt.0 .or. iy.gt.3119) then
         call penup
         return
      end if
      if (ipen.eq.0) call penup
      hy=iy/128
      lly=iy-128*hy
      ly=lly/4
      lly=lly-4*ly
      hx=ix/128
      llx=ix-128*hx
      lx=llx/4
      llx=llx-4*lx
      i=0
      if (hy.ne.hy0) then
         i=i+1
         coord(i)=char(hy+32)
      end if
      if (lly.ne.lly0 .or.llx.ne.llx0) then
         i=i+1
         coord(i)=char(96+4*lly+llx)
         i=i+1
         coord(i)=char(ly+96)
      elseif (ly.ne.ly0 .or. hx.ne.hx0) then
         i=i+1
         coord(i)=char(ly+96)
      end if
      if (hx.ne.hx0) then
         i=i+1
         coord(i)=char(hx+32)
      end if
      i=i+1
      coord(i)=char(lx+64)
      call pointqq(coord,i)
      hy0=hy
      ly0=ly
      hx0=hx
      lx0=lx
      lly0=lly
      llx0=llx
      return

      end subroutine plott



      subroutine pointqq(c,n)

      implicit none

      integer i,n,iunit
      character*1 c(n)
      common /unit/iunit
      if (iunit.eq.0) goto 1111
      do i=1,n
         write(iunit,100) c(i)
      end do
 100  format(a,$)
c      do 10 i=1,n
c 10      j=fputc(iunit,c(i))
c      call flush(iunit)
 1111 continue
      return

      end subroutine pointqq



c*****************************************************************
c
c     Contour plotting routines given to me by researchers at
c     NASA Ames in which a number of changes has been made.
c
c*****************************************************************

      subroutine contxx(idim,jdim,is,ie,js,je,xy,f,ncopt,ncont,acont,
     &     nxdim,xcont,ycont,naddim,nad,nlev,iadim,ia)
c
c     Calculate contour lines for the function f in the region is to ie, js to
c     je.  x,y coordinates corresponding to the grid points are in array xy.
c
c     If ncopt=1, figure our own contour levels, up to ncont of them, using
c     "nice" numbers.  frange finds the function range in the given region, and
c     consc1 computes the contour levels. Note that ncont will be revised
c     downward to correspond to the number of contour levels actually used.
c     The contour levels calculated are returned in the array acont.
c     NEW: zero contour will be removed by shifting half a contour spacing
c
c     If ncopt=2, calculate lines for the ncont contour levels specified in
c     acont.
c
c     The contour lines are returned in arrays xcont, and ycont.  nad(1)
c     gives the number of contour lines, and nad(n) points to the start of the
c     nth line (i.e. nad(n+1) points to one past the end of the nth line).
c     nlev(n) returns the contour level of the nth contour line.
c
c     ia is a scratch array. Try a dimension of 3000.
c
      implicit none

      integer i,j,ie,is,je,js,ncont,iw,np,npt,nlev,ik,ix,iy,ixa,iya
      integer m,n,ima,imb,ii,ncopt,ia,nad,nlinep,icont,idim,jdim
      integer nxdim,naddim,iadim
      real fmin,fmax,acont,aincre,wx,wy,xcont,ycont,wyy,wxx,xyf,xy
      real s,za,f
      dimension xy(idim,jdim,2),f(idim,jdim)
      dimension acont(ncont)
      dimension xcont(nxdim),ycont(nxdim)
      dimension nad(naddim),nlev(naddim),ia(iadim)
c
c     If ncopt=1, figure our own contour levels, up to ncont of them, using
c     "nice" numbers.  frange finds the function range in the given region, and
c     consc1 computes the contour levels.  note that ncont will be revised
c     downward to correspond to the number of contour levels actually used.
c
      if (ncopt.eq.1) then
         call frange(idim,jdim,is,ie,js,je,f,fmin,fmax)
         call consc1(fmin,fmax,ncont,acont)
         if (acont(1)*acont(ncont).le.0.) then
            aincre=acont(2)-acont(1)
            do ii=1,ncont
               acont(ii)=acont(ii)+aincre/2.
            end do
            if (acont(ncont).gt.fmax) ncont=ncont-1
            if (acont(1)-aincre.gt.fmin) then
               ncont=ncont+1
               do ii=ncont,2,-1
                  acont(ii)=acont(ii-1)
               end do
               acont(1)=acont(2)-aincre
            end if
         end if
      end if
c ***
      do 2222 i=1,iadim
2222  ia(i)=0
c ***
      iw=3
c
      nad(1)= 1
      nlinep= 2
c
c     One little check.  if is=ie or js=je, return with no contour lines.
c
      if (is.eq.ie .or. js.eq.je) goto 110
c
c     Loop through each contour level.
c
      do 100 icont= 1,ncont
       za= acont(icont)
      m=0
c
c     Scan points and determine points of ia
c
      do  600  j= js+1,je-1
         imb=0
         do  600  i= is,ie
            if (f(i,j).le.za)  goto  601
            if (imb.ne.1)  goto  600
            m=m+1
            if (m.gt.iadim)  goto  210
            ia(m)=1000*i+j
            imb=0
            goto  600
 601        imb=1
 600  continue
c
c     Search start point on boundary line
c
 101  ima=1
      imb=0
      ixa= is-1
      iya= js
 1    ixa=ixa+1
      if (ixa.eq.ie)  ima=2
      goto  5
 2    iya=iya+1
      if (iya.eq.je)  ima=3
      goto  5
 3    ixa=ixa-1
      if (ixa.eq.is)  ima=4
      goto  5
 4    iya=iya-1
      if (iya.eq.js)  ima=5
 5    if (f(ixa,iya).gt.za)  goto  7
      imb=1
 6    goto  (1,2,3,4,91),ima
 7    if (imb.ne.1)  goto  6
c
c     Determine start point
c
      imb=0
      ix=ixa
      iy=iya
      s=f(ixa,iya)
      goto  (21,11,12,13,51),ima
 11   if (iy.ne.js)  goto  31
      goto  21
 12   if (ix.ne.ie)  goto  41
      goto  31
 13   if (iy.ne.je)  goto  51
      goto  41
 10   ix=ia(n)/1000
      iy=ia(n)-1000*ix
      s=f(ix,iy)
      ia(n)=0
      goto  21
c
c     Process to search plot point
c
 20   iy=iy+1
 21   ix=ix-1
      if (ix.lt.is)  goto  90
      i=1
      if (f(ix,iy).le.za)  goto  52
      s=f(ix,iy)
      goto  31
 30   ix=ix-1
 31   iy=iy-1
      if (iy.lt.js)  goto  90
      i=2
      if (f(ix,iy).le.za)  goto  60
      s=f(ix,iy)
      goto  41
 40   iy=iy-1
 41   ix=ix+1
      if (ix.gt.ie)  goto  90
      i=3
      if (f(ix,iy).le.za)  goto  60
      s=f(ix,iy)
      goto  51
 50   ix=ix+1
 51   iy=iy+1
      i=4
      if (iy.gt.je)  goto  90
      if (f(ix,iy).le.za)  goto  60
      s=f(ix,iy)
      goto  21
 52   if (m.eq.0)  goto  60
      ik=1000*ix+iy+1000
      do  602  j=1,m
      if (ia(j).ne.ik)  goto  602
      ia(j)=0
 602  continue
c
c     Calculate plot point
c
 60   xyf=(za-f(ix,iy))/(s-f(ix,iy))
      goto  (61,62,63,64),i
 61   wxx= xy(ix,iy,1)+xyf*(xy(ix+1,iy,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix+1,iy,2)-xy(ix,iy,2))
      goto  65
 62   wxx= xy(ix,iy,1)+xyf*(xy(ix,iy+1,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix,iy+1,2)-xy(ix,iy,2))
      goto  65
 63   wxx= xy(ix,iy,1)+xyf*(xy(ix-1,iy,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix-1,iy,2)-xy(ix,iy,2))
      goto  65
 64   wxx= xy(ix,iy,1)+xyf*(xy(ix,iy-1,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix,iy-1,2)-xy(ix,iy,2))
c
c     Plot
c
 65   continue
c
c     Decide if plot point equal initial plot point
c
      if (iw.ne.3)  goto  66
      np=1
      nad(nlinep)= nad(nlinep-1)
      nlev(nlinep-1)= icont
      npt=0
      xcont(nad(nlinep))= wxx
      ycont(nad(nlinep))= wyy
      wx=wxx
      wy=wyy
      iw=2
      goto  67
 66   continue
      nad(nlinep)= nad(nlinep)+1
      if (nad(nlinep).gt.nxdim) goto 220
      np=np+1
      xcont(nad(nlinep))= wxx
      ycont(nad(nlinep))= wyy
c     if (np .lt. 200) goto 6602
c
c     call draw2d(pt,np,2,2,0)
c     npt=npt+np
c     np=1
c     pt(1,1)= wxx
c     pt(2,1)= wyy
 6602 if (wxx.ne.wx)  goto  67
      if (wyy.eq.wy)  goto  90
c
c     Determine next process
c
 67   goto  (50,20,30,40),i
 90   iw=3
      nad(nlinep)= nad(nlinep)+1
      if (nad(nlinep).gt.nxdim) goto 220
      if (np.gt.1) then
         nlinep= nlinep+1
         if (nlinep.gt.naddim) goto 230
      end if
c     if (np .gt. 1) call draw2d(pt,np,2,2,0)
      if (ima.ne.5)  goto  6
c
c     Search start point
c
      if (m.eq.0)  goto  92
 91   do  n=1,m
         if (ia(n).ne.0)  goto  10
      end do
 92   continue
c
c     Calculate value of next curve
c
 100  continue
c
 110  continue
      nad(1)= nlinep-2
      return
c
c     Warning - ia array full.
c
 210  continue
      write(6,211) iadim
 211  format('  warning - scratch array ia full in contour routine ',
     c     'contxx.'/
     c     '  picture may be incomplete.  array was dimensioned ',
     c     i5,'.')
      goto 110
c
c     Warning - xcont array full.
c
 220  continue
      write(6,221) nxdim
 221  format('  warning - contour line array xcont full in contour ',
     c     'routine contxx.'/
     c     '  picture may be incomplete.  array was dimensioned ',
     c     i5,'.')
      goto 110
c
c     Warning - nad array full.
c
 230  continue
      write(6,231) naddim
 231  format('  warning - contour line pointer array nad full in ',
     c     'contour routine contxx.'/
     c     '  picture may be incomplete.  array was dimensioned ',
     c     i5,'.')
      goto 110

      end subroutine contxx



      subroutine frange(idim,jdim,is,ie,js,je,f,fmin,fmax)
c
c     Find the range (minimum and maximum) of the function f in the region is to
c     ie, js to je.
c
      implicit none

      integer i,j,ie,is,je,js,idim,jdim
      real f,fmin,fmax
      dimension f(idim,jdim)
c
      fmin= f(is,js)
      fmax= fmin
      do j=js,je
         do i=is,ie
            fmin= min(fmin,f(i,j))
            fmax= max(fmax,f(i,j))
         end do
      end do
      return

      end subroutine frange



      subroutine consc1(amin,amax,ncont,acont)
c
c     Come up with a "nice" scaling of about ncont values between amin and amax.
c     ncont is updated to the number of intervals actually needed.
c
      implicit none

      integer i,ncont,nneed,imin,imax,ichar,nnice
      real acont,ainc,rmant,diff,char,amin,amax,rnice
      dimension acont(ncont)
      dimension rnice(4)
      data rnice/.1,.2,.25,.5/
      data nnice/4/
c
c     As a first approximation, get the difference, its characteristic and
c     mantissa.
c
      diff= (amax-amin)/(ncont+1)
      if (diff.le.0.) goto 20
      char= log10(diff)+1.
c
c     Round char down and get the mantissa.
c
      if (char.ge.0.) then
         ichar= char
      else
         ichar= char-1.
      end if
      rmant= diff*10.**(-ichar)
c
c     What's the next largest "nice" mantissa?
c
      do i= 1,nnice
         if (rmant.le.rnice(i)) goto 10
      end do
      i= nnice
c
c     Got a guess. Calculate a diff, round amin down.
c
   10 continue
      ainc= rnice(i)*10.**ichar
      imin= amin/ainc
      if (amin.gt.imin*ainc) imin= imin+1
      imax= amax/ainc
      if (amax.lt.imax*ainc) imax= imax-1
      nneed= imax+1-imin
c
c     Are we under?
c
      if (nneed.gt.ncont) then
c
c     Nope. Try the next nice number.
c
         if (i.lt.nnice) then
            i= i+1
         else
            ichar= ichar+1
            i= 1
         end if
         goto 10
      end if
c
c     Now just set up the acont array and update ncont.  (set up all (original)
c     ncont of the aconts for the folks back home.)
c
      do 1 i= 1,ncont
         acont(i)= (imin-1+i)*ainc
 1    continue
      ncont= nneed
      goto 30
c
c     All values are the same -- just set up one contour level.
c
 20   continue
      acont(1)= amin
      do i=2,ncont
         acont(i)= 0.
      end do
      ncont= 1
c
   30 continue
c
      return

      end subroutine consc1
