c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program par1
      real x(2000,5),y(2000,5),x1(2000,5),y1(2000,5)
      real xp(8000,5),yp(8000,5)
      character*16 xaxis,yaxis,heada,headb,headc
      integer npoint(5),npoint2(5)
      common/hard/ihard
      real pi
      parameter (pi = 3.1415926535897932385)
c
      open(unit=49,file='uvq2.dat')
      read(49,*) n1
      do i=1,n1
        read(49,*) c,npoint(i)
        do j=1,npoint(i)
          read(49,*) k,x(j,i),y(j,i)
        end do
      end do
      close(unit=49)
      open(unit=49,file='rlsl.dat')
      read(49,*) n3
      do i=1,n3
        read(49,*) c,npoint2(j)
        do j=1,npoint2(j)
          read(49,*) k,x1(j,i),y1(j,i)
        end do
      end do
      close(unit=49)
      nplot=min(n1,n3)
c
c cubical point mid-point insertion
c
      do i=1,nplot
c we start with npoint(i) points and insert three extra points
c into the intervals between point 2 to point npoint(i)-1
c i.e. 3*(npoint(i)-3) extra points
c
       xp(1,i)=y1(1,i)
       xp(2,i)=y1(2,i)
       xp(4*npoint(i)-9,i)=y1(npoint(i),i)
       yp(1,i)=y(1,i)
       yp(2,i)=y(2,i)
       yp(4*npoint(i)-9,i)=y(npoint(i),i)
       do j=3,npoint(i)-1
       xp(4*j-9,i)=
     &    (-7.*y1(j-2,i)+105.*y1(j-1,i)+35.*y1(j,i)-5.*y1(j+1,i))/128.
       xp(4*j-8,i)=(-y1(j-2,i)+9.*y1(j-1,i)+9.*y1(j,i)-y1(j+1,i))/16.
       xp(4*j-7,i)=
     &    (-5.*y1(j-2,i)+35.*y1(j-1,i)+105.*y1(j,i)-7.*y1(j+1,i))/128.
       xp(4*j-6,i)=y1(j,i)
       yp(4*j-9,i)=
     &      (-7.*y(j-2,i)+105.*y(j-1,i)+35.*y(j,i)-5.*y(j+1,i))/128.
       yp(4*j-8,i)=(-y(j-2,i)+9.*y(j-1,i)+9.*y(j,i)-y(j+1,i))/16.
       yp(4*j-7,i)=
     &      (-5.*y(j-2,i)+35.*y(j-1,i)+105.*y(j,i)-7.*y(j+1,i))/128.
       yp(4*j-6,i)=y(j,i)
       end do
       npoint(i)=4*npoint(i)-9
      end do

      xmin=0.
      xmax=150.
      ymin=0.
      ymax=0.20
      yaxis='uv/q2'
      xaxis='r_lams_lam'
      heada=' '
      headb=' '
      headc=' '
      idash=1
      mbox=1
      call rita1(xp,yp,xmin,xmax,ymin,ymax,npoint,nplot,8000,
     &        xaxis,yaxis,heada,headb,headc,mbox,idash)
      end

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
c      do 100 i=1,101
c         x(i,1)=(i-1)/100.
c         y(i,1)=x(i,1)
c         x(i,2)=(i-1)/100.
c         y(i,2)=x(i,1)**2
c         x(i,3)=(i-1)/100.
c         y(i,3)=x(i,1)**3
c         x(i,4)=(i-1)/100.
c         y(i,4)=x(i,1)**4
c 100  continue
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
c      do 100 i=1,501
c         x(i,1)=i
c         y(i,1)=x(i,1)
c         x(i,2)=i
c         y(i,2)=x(i,1)**2
c         x(i,3)=i
c         y(i,3)=x(i,1)**3
c         x(i,4)=i
c         y(i,4)=x(i,1)**4
c 100  continue
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
c      do 20 i=1,51
c         do 20 j=1,51
c            xy(i,j,1)=(i-1)*xconst
c            xy(i,j,2)=(j-1)*yconst
c            u1(i,j)=((i-1)*xconst)**2+((j-1)*yconst)**2-2.
c 20   continue
c
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
c     if xmax1.eq.xmin1  choose xmax and xmin from data
c     if ymax1.eq.ymin1  choose ymax and ymin from data
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
      if(xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do 23 j=1,nplot
         do 23 i=1,npoint(j)
            xmax=max(xmax,x(i,j))
            xmin=min(xmin,x(i,j))
 23      continue
         xdiff=xmax-xmin
         if(xmax.ne.0) xmax=xmax+0.1*xdiff
         if(xmin.ne.0) xmin=xmin-0.1*xdiff
         if(xmin.eq.xmax) then
           xmin=xmin-0.0001
           xmax=xmax+0.0001
         end if
      end if
      if(ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do 24 j=1,nplot
         do 24 i=1,npoint(j)
            ymax=max(ymax,y(i,j))
            ymin=min(ymin,y(i,j))
 24      continue
         ydiff=ymax-ymin
         if(ymax.ne.0) ymax=ymax+0.1*ydiff
         if(ymin.ne.0) ymin=ymin-0.1*ydiff
         if(ymin.eq.ymax) then
           ymin=ymin-0.0001
           ymax=ymax+0.0001
         end if
      end if
      xtick=0.
      ytick=0.
      write(*,*) ' termtype: 1=tek, 2=xterm'
      read(*,*) iterm1
 32   if(iterm1.eq.1) call tek
      if(iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      do 100 j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
 100  continue
      call plend
      return
      end

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
      if(xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do 23 j=1,nplot
         do 23 i=1,npoint(j)
            xmax=max(xmax,x(i,j))
            xmin=min(xmin,x(i,j))
 23      continue
         xdiff=xmax-xmin
         if(xmax.ne.0) xmax=xmax+0.1*xdiff
         if(xmin.ne.0) xmin=xmin-0.1*xdiff
         if(xmin.eq.xmax) then
           xmin=xmin-0.0001
           xmax=xmax+0.0001
         end if
      end if
      if(ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do 24 j=1,nplot
         do 24 i=1,npoint(j)
            ymax=max(ymax,y(i,j))
            ymin=min(ymin,y(i,j))
 24      continue
         ydiff=ymax-ymin
         if(ymax.ne.0) ymax=ymax+0.1*ydiff
         if(ymin.ne.0) ymin=ymin-0.1*ydiff
         if(ymin.eq.ymax) then
           ymin=ymin-0.0001
           ymax=ymax+0.0001
         end if
      end if
      xtick=0.
      ytick=0.
c      write(*,*) ' termtype: 1=tek, 2=xterm'
c      read(*,*) iterm1
      iterm1=2
 32   if(iterm1.eq.1) call tek
      if(iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,xtick,xmax,ymin,ytick,ymax)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      dy=(ymax-ymin)*.01
      do 100 j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
100   continue
      call plend
      return
      end

      subroutine ritlog(x,y,xmin1,xmax1,ymin1,ymax1,npoint,nplot,ndim,
     &                xaxis,yaxis,heada,headb,headc,ilog,idash)
c
c     if xmax1.eq.xmin1  choose xmax and xmin from data
c     if ymax1.eq.ymin1  choose ymax and ymin from data
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
      if(xmin.eq.xmax) then
         xmin=1.e20
         xmax=-1.e20
         do 23 j=1,nplot
         do 23 i=1,npoint(j)
            xmax=max(xmax,x(i,j))
            xmin=min(xmin,x(i,j))
 23      continue
c         xdiff=xmax-xmin
c         if(xmax.ne.0) xmax=xmax+0.1*xdiff
c         if(xmin.ne.0) xmin=xmin-0.1*xdiff
      end if
      if(ymin.eq.ymax) then
         ymin=1.e20
         ymax=-1.e20
         do 24 j=1,nplot
         do 24 i=1,npoint(j)
            ymax=max(ymax,y(i,j))
            ymin=min(ymin,y(i,j))
 24      continue
c         ydiff=ymax-ymin
c         if(ymax.ne.0) ymax=ymax+0.1*ydiff
c         if(ymin.ne.0) ymin=ymin-0.1*ydiff
      end if

      if(ilog.eq.1.or.ilog.eq.3) then
         do 25 j=1,nplot
         do 25 i=1,npoint(j)
            if(x(i,j).le.0.0) then
               write(*,*) 'x-array does not contain positive numbners'
               return
            end if
 25      continue
      end if
      if(ilog.eq.2.or.ilog.eq.3) then
         do 26 j=1,nplot
         do 26 i=1,npoint(j)
            if(y(i,j).le.0.0) then
               write(*,*) 'y-array does not contain positive numbners'
               return
            end if
 26      continue
      end if

      if(ilog.eq.1.or.ilog.eq.3) then
         do 27 j=1,nplot
         do 27 i=1,npoint(j)
            x(i,j)=log10(x(i,j))
 27      continue
         xmin=log10(xmin)
         xmax=log10(xmax)
      end if
      if(ilog.eq.2.or.ilog.eq.3) then
         do 28 j=1,nplot
         do 28 i=1,npoint(j)
            y(i,j)=log10(y(i,j))
 28      continue
         ymin=log10(ymin)
         ymax=log10(ymax)
      end if
      xtick=0.
      ytick=0.
      write(*,*) ' termtype: 1=tek, 2=xterm'
      read(*,*) iterm1
 32   if(iterm1.eq.1) call tek
      if(iterm1.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call logaxis(xmin,xtick,xmax,ymin,ytick,ymax,ilog)
      call head1(heada)
      call head2(headb)
      call head3(headc)
      call taxis(xaxis,yaxis)
      do 100 j=1,nplot
         jdash=j*idash
         call curve(x(1,j),y(1,j),npoint(j),jdash)
 100  continue
      call plend
      if(ilog.eq.1.or.ilog.eq.3) then
         do 34 j=1,nplot
         do 34 i=1,npoint(j)
            x(i,j)=10.**(x(i,j))
 34      continue
      end if
      if(ilog.eq.2.or.ilog.eq.3) then
         do 35 j=1,nplot
         do 35 i=1,npoint(j)
            y(i,j)=10.**(y(i,j))
 35      continue
      end if
      return
      end

      subroutine cont1(u1,xy,nx,ny,xaxis,yaxis,heada,headb,mbox)
c
c     calculates contours of values in an array u1 with associated
c     coordinates xy
c
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
      if(ncont.eq.0) then
          ncopt=2
          write(*,*) ' give start, end and increment of contours'
          read(*,*) start,end,aincr
          ncont=int((end-start)/aincr+.01)+1
          do 100 n=1,ncont
             acont(n)=start+aincr*real(n-1)
100       continue
      end if
c
      call contxx(nx,ny,1,nx,1,ny,xy,u1,ncopt,ncont,acont,
     *         nxdim,xcont,ycont,naddim,nad,nlev,nxdim,work)
c
      if(ncont.eq.1) then
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
      if(iterm.eq.1) call tek
      if(iterm.eq.2) call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,0.,xmax,ymin,0.,ymax)
      call head1(heada)
      call head2(headb)
      call head3(text)
      call taxis(xaxis,yaxis)
c
      isnext=1
      nd=nad(1)
      do 500 ic=2,nd+1
          is=isnext
          isnext=nad(ic)
          idash=1
          if(acont(nlev(ic-1)).lt.0) idash=2
          call curve(xcont(is),ycont(is),isnext-is,idash)
 500  continue
      call plend
c      goto 2
      end

c********************************************************************
c
c     subroutine included for compatibility with older versions
c     of the plotting package, the same thing can now be done
c     better with rita1 or rita1a using the new version of curve
c
       subroutine rita2(x,y,xmin1,xmax1,ymin1,ymax1,npoint,
     &                xaxis,yaxis,heada,headb,headc,mbox)
c
c     plots boxes for the points in the coordinate positions
c
c     npoint           number of points in plot
c     mbox             0: prop. as specified, 1: max plot area
c     xaxis            name of xaxis
c     heada            first heading
c
       real x(1),y(1)
       character*(*) xaxis,yaxis,heada,headb,headc
      common /hard/ihard
      ihard=-1
c
      xmin=xmin1
      xmax=xmax1
      ymin=ymin1
      ymax=ymax1
      if(xmin.eq.xmax) then
         xmax=-1.e10
         xmin=1.e10
         do 23 i=1,npoint
            xmax=max(xmax,x(i))
            xmin=min(xmin,x(i))
 23      continue
         xdiff=xmax-xmin
         if(xmax.ne.0) xmax=xmax+0.1*xdiff
         if(xmin.ne.0) xmin=xmin-0.1*xdiff
      end if
      if(ymin.eq.ymax) then
         ymax=-1.e10
         ymin=1.e10
         do 24 i=1,npoint
            ymax=max(ymax,y(i))
            ymin=min(ymin,y(i))
 24      continue
         ydiff=ymax-ymin
         if(ymax.ne.0) ymax=ymax+0.1*ydiff
         if(ymin.ne.0) ymin=ymin-0.1*ydiff
      end if
       xtick=0.
       ytick=0.
       write(*,*) ' termtype: 1=tek, 2=xterm'
       read(*,*) iterm1
 32   if(iterm1.eq.1) call tek
       if(iterm1.eq.2) call xterm1
       call plarea(xmin,xmax,ymin,ymax,mbox)
       call axis(xmin,xtick,xmax,ymin,ytick,ymax)
       call head1(heada)
       call head2(headb)
       call head3(headc)
       call taxis(xaxis,yaxis)
       do 100 i=1,npoint
          if(y(i).lt.0.5.or.y(i).gt.-1.) call box1(x(i),y(i),0.02)
  100  continue
       call plend
       return
       end

       subroutine box1(x,y,al)
       al2=al/2.
       call plotr(x+al2,y+al2,0)
       call plotr(x+al2,y-al2,1)
       call plotr(x-al2,y-al2,1)
       call plotr(x-al2,y+al2,1)
       call plotr(x+al2,y+al2,1)
       return
       end

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
      real xr(1),yr(1)
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      if(idash.ge.0) then
         call setdash(idash,hd1,hs1,hd2,hs2,h1)
         i1=0
         do 10 i=1,n-1
            x1=xr(i)*xscale+xoff
            y1=yr(i)*yscale+yoff
            x2=xr(i+1)*xscale+xoff
            y2=yr(i+1)*yscale+yoff
            call dash(h1,i1,x1,y1,x2,y2,hd1,hs1,hd2,hs2)
            if(idash.eq.1.and.mod(i,500).eq.0) call plota(x2,y2,0)
 10      continue
      else
         isymb=abs(idash)
         do 20 i=1,n
            x1=xr(i)*xscale+xoff
            y1=yr(i)*yscale+yoff
            call boxes(x1,y1,0.014,isymb)
 20      continue
      end if
      return
      end

      subroutine plarea(xmin1,xmax1,ymin1,ymax1,mbox1)
c
c     Defines plotting area and user coordinates
c     mbox=1: maximum ploting area, 2: user specified proportions
c
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
      if(mbox.gt.1) xsize=real(mbox)/1000.
      if(mbox.lt.0) ysize=real(-mbox)/1000.
      if(mbox.eq.0) then
         sizeq=xsize/ysize
         sizeb=(xmax-xmin)/(ymax-ymin)
         if(sizeb.ge.sizeq) ysize=ysize*sizeq/sizeb
         if(sizeb.lt.sizeq) xsize=xsize*sizeb/sizeq
      end if
      xor=.75-.5*xsize
      xscale=xsize/(xmax-xmin)
      yscale=ysize/(ymax-ymin)
      xoff=xor-xmin*xscale
      yoff=yor-ymin*yscale
      return
      end

      subroutine across(xmin,xmax,ymin,ymax)
c
c     draw axis cross
c
      if(xmin.lt.0.0.and.xmax.gt.0.0) then
         call plotr(0.,ymin,0)
         call plotr(0.,ymax,1)
      end if
      if(ymin.lt.0.0.and.ymax.gt.0.0) then
         call plotr(xmin,0.,0)
         call plotr(xmax,0.,1)
      end if
      return
      end

      subroutine axis(xmin1,xtick1,xmax1,ymin1,ytick1,ymax1)
c
c     draw axis with tickmarks at xtick and ytick intervals
c     if xtick or ytick are 0. automatic calculation of nice tick marks
c     are done
c
      character*100 num
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xmin=xmin1
      xtick=xtick1
      xmax=xmax1
      ymin=ymin1
      ytick=ytick1
      ymax=ymax1
      call font(2,charx,chary)
      if(xtick.eq.0.) then
         itick=9
         if(xsize.le.0.5) itick=5
         call  conscl(xmin,xmax,itick,xtick)
      end if
      if(ytick.eq.0.) call conscl(ymin,ymax,9,ytick)
      call frame(xmin,xmax,ymin,ymax)
      call across(xmin,xmax,ymin,ymax)
      istart=xmin/xtick
      if(xmin.gt.0.0) istart=istart+1
      iend=xmax/xtick
      if(xmax.lt.0.0) iend=iend-1
      do 10 i=istart,iend
         x=i*xtick
         call plotr(x,ymin,0)
         call plotr(x,ymin+charx/yscale,1)
         call flstr(x,num,xtick)
         call xxa(x,xa)
         call writ(num,xa,yor,charx,chary,1,-1.5,-.5)
         call plotr(x,ymax,0)
         call plotr(x,ymax-charx/yscale,1)
 10   continue
      istart=ymin/ytick
      if(ymin.gt.0.0) istart=istart+1
      iend=ymax/ytick
      if(ymax.lt.0.0) iend=iend-1
      do 20 i=istart,iend
         y=i*ytick
         call plotr(xmin,y,0)
         call plotr(xmin+charx/xscale,y,1)
         call flstr(y,num,ytick)
         call yya(y,ya)
         call writ(num,xor,ya,charx,chary,1,-.5,-1.2)
         call plotr(xmax,y,0)
         call plotr(xmax-charx/xscale,y,1)
 20   continue
      return
      end

      subroutine logaxis(xmin1,xtick1,xmax1,ymin1,ytick1,ymax1,ilog)
c
c     draw logaxis with tickmarks at xtick and ytick intervals
c     if xtick or ytick are 0. automatic calculation of nice tick marks
c     are done, for logaxis tickmarks are also generated automatically
c
      character*100 num
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xmin=xmin1
      xtick=xtick1
      xmax=xmax1
      ymin=ymin1
      ytick=ytick1
      ymax=ymax1
      call font(2,charx,chary)
      if(xtick.eq.0..and.ilog.ne.1.and.ilog.ne.3) then
         itick=9
         if(xsize.le.0.5) itick=5
         call  conscl(xmin,xmax,itick,xtick)
      end if
      if(ytick.eq.0..and.ilog.ne.2.and.ilog.ne.3)
     &     call conscl(ymin,ymax,9,ytick)
      call frame(xmin,xmax,ymin,ymax)
c      call across(xmin,xmax,ymin,ymax)

      if(ilog.ne.1.and.ilog.ne.3) then
         istart=xmin/xtick
         if(xmin.gt.0.0) istart=istart+1
         iend=xmax/xtick
         if(xmax.lt.0.0) iend=iend-1
         do 10 i=istart,iend
            x=i*xtick
            call plotr(x,ymin,0)
            call plotr(x,ymin+charx/yscale,1)
            call flstr(x,num,xtick)
            call xxa(x,xa)
            call writ(num,xa,yor,charx,chary,1,-1.5,-.5)
            call plotr(x,ymax,0)
            call plotr(x,ymax-charx/yscale,1)
 10      continue
      else
         xtick=100000000.
         if(xmin.gt.0.0) istart=int(xmin)+1
         if(xmin.le.0.0) istart=int(xmin)
         if(xmax.gt.0.0) iend=int(xmax)
         if(xmax.lt.0.0) iend=int(xmax)-1

         if(iend-istart.gt.10) goto 110
         do 111 itick=1,4
            x1=istart-1+log10(2.*itick)
            if(x1.lt.xmin) goto 111
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

            if(iend-istart.gt.10) goto 12
            do 112 itick=1,4
               x1=x+log10(2.*itick)
               if(x1.gt.xmax) goto 112
               call plotr(x1,ymin,0)
               call plotr(x1,ymin+0.5*charx/yscale,1)
               call plotr(x1,ymax,0)
               call plotr(x1,ymax-0.5*charx/yscale,1)
 112        continue

 12      continue
      end if

      if(ilog.ne.2.and.ilog.ne.3) then
         istart=ymin/ytick
         if(ymin.gt.0.0) istart=istart+1
         iend=ymax/ytick
         if(ymax.lt.0.0) iend=iend-1
         do 20 i=istart,iend
            y=i*ytick
            call plotr(xmin,y,0)
            call plotr(xmin+charx/xscale,y,1)
            call flstr(y,num,ytick)
            call yya(y,ya)
            call writ(num,xor,ya,charx,chary,1,-.5,-1.2)
            call plotr(xmax,y,0)
            call plotr(xmax-charx/xscale,y,1)
 20      continue
      else
         ytick=100000000.
         if(ymin.gt.0.0) istart=int(ymin)+1
         if(ymin.le.0.0) istart=int(ymin)
         if(ymax.gt.0.0) iend=int(ymax)
         if(ymax.lt.0.0) iend=int(ymax)-1

         if(iend-istart.gt.10) goto 120
         do 121 itick=1,4
            y1=istart-1+log10(2.*itick)
            if(y1.lt.ymin) goto 121
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

            if(iend-istart.gt.10) goto 22
            do 122 itick=1,4
               y1=y+log10(2.*itick)
               if(y1.gt.ymax) goto 122
               call plotr(xmin,y1,0)
               call plotr(xmin+0.5*charx/xscale,y1,1)
               call plotr(xmax,y1,0)
               call plotr(xmax-0.5*charx/xscale,y1,1)
 122        continue

 22      continue
      end if
      return
      end

      subroutine plotr(xr,yr,ipen)
c
c     ipen=0: move to user coordinate (x,y)
c     ipen=1: draw line to user coordinate (x,y)
c
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      x=xr*xscale+xoff
      y=yr*yscale+yoff
      call plota(x,y,ipen)
      return
      end

      subroutine xxa(x,xa)
c
c     convert form user x to absolute x
c
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      xa=x*xscale+xoff
      return
      end

      subroutine yya(y,ya)
c
c     convert from user y to absolute ya
c
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      ya=y*yscale+yoff
      return
      end

      subroutine conscl(amin,amax,ncont,ainc)
c
c     come up with a "nice" scaling of about ncont values between
c     amin and amax.
c
      logical rev
      dimension rnice(3)
      data rnice/.1,.2,.5/
      data nnice/3/
c
c     reverse amin,amax if necessary (change back before return)
c
      rev=amax.lt.amin
      if(rev) then
         temp=amin
         amin=amax
         amax=temp
      end if
c
c     as a first approximation, get the difference, its characteristic and
c     mantissa.
c
      diff= (amax-amin)/(ncont+1)
      if (diff.le.0.) goto 20
      char= log10(diff)+1.
c
c     round char down and get the mantissa.
c
      if (char.ge.0.) then
         ichar= char
      else
         ichar= char-1.
      end if
      rmant= diff*10.**(-ichar)
c
c     what's the next largest "nice" mantissa?
c
      do 3 i= 1,nnice
         if (rmant.le.rnice(i)) goto 10
 3    continue
      i= nnice
c
c     got a guess.  calculate a diff, round amin down.
c
 10   continue
      ainc= rnice(i)*10.**ichar
      imin= amin/ainc
      if (amin.gt.imin*ainc) imin= imin+1
      imax= amax/ainc
      if (amax.lt.imax*ainc) imax= imax-1
      nneed= imax+1-imin
c
c     are we under?
c
      if (nneed.gt.ncont) then
c
c     nope.  try the next nice number.
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
c     if amin,amax reversed change back and change sign of ainc
c
      if(rev) then
         temp=amin
         amin=amax
         amax=temp
         ainc=-ainc
      end if
      return
      end

c******************************************************************
c
c     text and string handling routines
c
c******************************************************************

      subroutine flstr(x,num,xtick)
c
c     read the floating point number x to string num
c     the size of xtick determines the number of decimal places
c
      character*(*) num
      character*100 form
c
      if(x.eq.0.) then
         ilen=3
      else
         ilen=int(log10(abs(x)))
         if(ilen.le.0) then
            ilen=3
         else
            ilen=ilen+3
         end if
      end if
c
      sci=log10(xtick*1.001)
      idec=-int(sci)
c      idec=-int(log10(xtick*1.001))
      if(idec.lt.0) then
         idec=0
      else if(idec.eq.0.and.xtick.ge.1.) then
         idec=0
      else
         idec=idec+1
      end if
c
      if(abs(sci).le.4) then
         ilen=ilen+idec
         if(ilen.lt.10) write(form,97) ilen,idec
         if(ilen.ge.10.and.idec.lt.10) write(form,98) ilen,idec
         if(idec.ge.10) write(form,99) ilen,idec
 97      format(' (f',i1,'.',i1,') ')
 98      format(' (f',i2,'.',i1,') ')
 99      format(' (f',i2,'.',i2,') ')
      else
         form='(1pe8.1,a1)'
      end if
c
c      write(*,*) x
      write(num,form) x
c
      return
      end

      subroutine head1(text)
c
c     top heading
c
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(1,charx,chary)
      call writ(text,xor+xsize/2,0.9,charx,chary,1,-.5,-.5)
      return
      end

      subroutine head3(text)
c
c     bottom heading
c
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(text,xor+xsize/2,0.8,charx,chary,1,-.5,-.5)
      return
      end

      subroutine head2(text)
c
c     middle heading
c
      character*(*) text
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(text,xor+xsize/2,0.85,charx,chary,1,-.5,-.5)
      return
      end

      subroutine taxis(xtext,ytext)
c
c     text on x- and y-axis
c
      character*(*) xtext,ytext
      common/draw/xor,yor,xsize,ysize,xoff,yoff,xscale,yscale
      call font(2,charx,chary)
      call writ(xtext,xor+xsize/2,0.075,charx,chary,1,-.5,-.5)
      call writ(ytext,xor-0.15,yor+ysize/2,charx,chary,2,-.5,-.5)
      return
      end

c*****************************************************************
c
c     line type and symbol routines
c
c*****************************************************************

      subroutine dash(h1,i,x1,y1,x2,y2,hd1,hs1,hd2,hs2)
c
c      h1 - distance left of dash/space
c      i  : i < 0 curr. space ; i > 0 curr. dash ; i=0 start
c      hdi/hsi  - length of dashes/spaces
c
      real l
      if(i.eq.0) then
          i=1
          call plota(x1,y1,0)
      end if
      x=x1
      y=y1
1      l=sqrt((x-x2)**2+(y-y2)**2)
      if(h1.gt.l) then
            ipen=1
            if(i.lt.0) ipen=0
            call plota(x2,y2,ipen)
            ipen=1
            h1=h1-l
      else
            xm=x+h1*(x2-x)/l
            ym=y+h1*(y2-y)/l
            ipen=1
            if(i.lt.0) ipen=0
            call plota(xm,ym,ipen)
            ipen=1
            if(i.eq.1) then
                  i=-1
                  h1=hs1
            else if(i.eq.-1) then
                  i=2
                  h1=hd2
            else if(i.eq.2) then
                  i=-2
                  h1=hs2
            else if(i.eq.-2) then
                  i=1
                  h1=hd1
            end if
            x=xm
            y=ym
            goto 1
      end if
      return
      end

      subroutine setdash(j,hd1,hs1,hd2,hs2,h1)
      dmmtek=0.005714
      if(j.le.1) then
          hd1=dmmtek*1000000.0
          hs1=dmmtek*1.0
          hd2=dmmtek*1.0
          hs2=dmmtek*1.0
      end if
      if(j.eq.2) then
          hd1=dmmtek*1.5
          hs1=dmmtek*1.0
          hd2=dmmtek*1.5
          hs2=dmmtek*1.0
      end if
      if(j.eq.3) then
          hd1=dmmtek*0.3
          hs1=dmmtek*1.25
          hd2=dmmtek*0.3
          hs2=dmmtek*1.25
      end if
      if(j.eq.4) then
          hd1=dmmtek*9.0
          hs1=dmmtek*2.5
          hd2=dmmtek*2.5
          hs2=dmmtek*2.5
      end if
      if(j.eq.5) then
          hd1=dmmtek*7.0
          hs1=dmmtek*2.5
          hd2=dmmtek*0.3
          hs2=dmmtek*2.5
      end if
      if(j.ge.6) then
          hd1=dmmtek*10.0
          hs1=dmmtek*5.0
          hd2=dmmtek*10.0
          hs2=dmmtek*5.0
      end if
      h1=hd1
      return
      end

      subroutine boxes(x,y,aa,i)
      if(i.eq.1) then
         a=aa/2.
         call plota(x+a,y+a,0)
         call plota(x+a,y-a,1)
         call plota(x-a,y-a,1)
         call plota(x-a,y+a,1)
         call plota(x+a,y+a,1)
      end if
      if(i.eq.2) then
         a=aa/2.
         call plota(x+a,y+a,0)
         call plota(x,y-a,1)
         call plota(x-a,y+a,1)
         call plota(x+a,y+a,1)
      end if
      if(i.ge.3) then
         a=aa/1.414
         call plota(x,y+a,0)
         call plota(x+a,y,1)
         call plota(x,y-a,1)
         call plota(x-a,y,1)
         call plota(x,y+a,1)
      end if
      return
      end

      subroutine arrow(xr1,yr1,xr2,yr2,d)

c     draws an arrow between point (xr1,yr1) and (xr2,yr2)
c     coordinates are assumed to be user defined and note that
c     we are back at (xr2,yr2) when the routine exits

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
      end

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
      call plott(x,y,ipen)
      call ploth(x,y,ipen)
      return
      end

      subroutine frame(xmin,xmax,ymin,ymax)
c
c     draw frame of plot, make it thicker for postscript
c
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
      end

      subroutine writ(text,x,y,charx,chary,irot,yjust,xjust)
c
c     write text to plot, note differences between tektronix and postscript
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
      character*(*) text
      character*100 name
      common/unit/iunit
      do 8 i=1,100
         name(i:i)=' '
 8    continue
      n=index(text,'$')
      if(n.eq.0) then
c find the last non-blank character
        do 2000 n=len(text),1,-1
          if(text(n:n).ne.' ') goto 1000
 2000   continue
        n=0
 1000   continue
        n=n+1
      end if
      if(iunit.eq.0) goto 1111
      if(irot.eq.1) then
         call plott(x+n*charx*xjust,y+chary*yjust,0)
         call tmode
         write(iunit,111) text(1:n-1)
 111     format(a)
         call gmode
      end if
      if(irot.eq.2) then
         call tmode
         do 10 i=1,n-1
            call plott(x+charx*yjust,y-n*chary*xjust-(i-1)*chary,0)
            call tmode
            write(iunit,111) text(i:i)
            call gmode
 10      continue
      end if
 1111 continue
c
      call ploth(x,y,0)
      name='('//text(1:n-1)//') p'
      call  writp(name,irot,yjust,xjust)
c
      return
      end

      subroutine font(ifont,charx,chary)
c
c     4 default font sizes, sizes in absolute units for use in tektronix mode
c
      character*1 c(2)
      c(1)=char(27)
      if(ifont.eq.1) then
         c(2)='8'
         charx=0.0184
         chary=charx*1.4
      end if
      if(ifont.eq.2) then
         c(2)='9'
         charx=0.0151
         chary=charx*1.4
      end if
      if(ifont.eq.3) then
         c(2)=':'
         charx=0.0124
         chary=charx*1.4
      end if
      if(ifont.eq.4) then
         c(2)=';'
         charx=0.0101
         chary=charx*1.4
      end if
      call pointqq(c,2)
c
      call pfont(ifont)
c
      return
      end

      subroutine xterm1
c
c     initializes the tektronix plotting, clears the screen and
c     moves to the tektronix window associated with the xterm
c     window.
c
      character*1 c(6)
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
      end

      subroutine tek
c
c     same as xterm1 except that it does not assume X environment
c
      common/term/iterm
      call hunit
      iterm=0
      call terminal
      call erase
      call pinit
      return
      end

      subroutine hunit
c
c     initializes tek and ps units
c     iunit=0 means no tek output
c     kunit=0 means no ps  output
c
      character*100 name,form
      common/hard/ihard
      common/unit/iunit
      common/punit/kunit
c
      if(ihard.eq.0) then
         iunit=6
         kunit=0
      end if
c
      munit=79
      if(ihard.gt.0) then
c         ierr=0
c 1000    if(ierr.gt.0) munit=munit+1
c         if(munit.lt.100) form=' (a5,i2) '
c         if(munit.ge.100) form=' (a5,i3) '
c         write(name,form) 'plot.',munit
c         open(unit=munit,file=name,status='new',err=1000,iostat=ierr)
c
c new way of opening the plot files (should work on alliant)
         ierr=1
 1000    if(ierr.eq.0) munit=munit+1
         if(munit.lt.100) form=' (a5,i2) '
         if(munit.ge.100) form=' (a5,i3) '
         write(name,form) 'plot.',munit
         open(unit=munit,file=name,status='old',err=1010,iostat=ierr)
         close(unit=munit)
         goto 1000
 1010    open(unit=munit,file=name,status='new')
      end if
c
      if(abs(ihard).eq.1) then
         iunit=6
         kunit=munit
      end if
      if(abs(ihard).eq.2) then
         iunit=0
         kunit=munit
      end if
      if(abs(ihard).eq.3) then
         iunit=munit
         kunit=0
      end if
c
      end

      subroutine plend
c
c     sets mode back to text and in case of X environment moves
c     back to the xterm window
c
      character*1 c(2)
      common/term/iterm
      common/hard/ihard
      common/unit/iunit
      common/punit/kunit
      munit=max(iunit,kunit)
      call tmode
      c(1)=char(27)
      c(2)=char(3)
      if(iterm.eq.1) call pointqq(c,2)
      call pend
      if(ihard.gt.0.and.munit.ne.6) close(munit)
      return
      end

c*****************************************************************
c
c     Postscript driver routines
c
c*****************************************************************

      subroutine pinit
c
c     initialize the postscript file
c     the /p routine takes four arguments in the following manner:
c     irot  yjust  xjust  string  p
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
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
c     put in after 'gsave' line in format 111 for side plot:
c
c     &   '600 0 translate',/,
c     &   '90 rotate',/,
c
c     also change scale in ploth routine,
c     remove multiplication with 0.76923077
c
      end

      subroutine pend
c
c     finalize the postscript file
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
      write(kunit,*) 's showpage grestore'
 1111 continue
      return
      end

      subroutine closep
c
c     postscript close path routine for use in frame
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
      write(kunit,*) 'c'
 1111 continue
      return
      end

      subroutine pfont(ifont)
c
c     the 4 postdcript default font sizes
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
      if(ifont.eq.1) write(kunit,*) '/fs 16 def'
      if(ifont.eq.2) write(kunit,*) '/fs 14 def'
      if(ifont.eq.3) write(kunit,*) '/fs 12 def'
      if(ifont.eq.4) write(kunit,*) '/fs 10 def'
      write(kunit,*) 'ff'
 1111 continue
      return
      end

      subroutine writp(name,irot,yjust,xjust)
c
c     postscript print text
c
      character*(*) name
      common/punit/kunit
      if(kunit.eq.0) goto 1111
      n=index(name,') p')
      write(kunit,111) irot,yjust,xjust,name(1:n+2)
 111  format(i2,' ',f5.2,' ',f5.2,' ',100a)
 1111 continue
      return
      end

      subroutine ploth(x,y,ipen)
c
c     write the line drawing commands to the postscript file
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
c
      if(x.gt.1.5.or.x.lt.-0.2) then
         write(kunit,110)
         goto 113
      end if
      if(y.gt.1.2.or.y.lt.-0.2) then
         write(kunit,110)
         goto 113
      end if
 110  format('s ')
c
      xp=600.*x*0.76923077
      yp=600.*y*0.76923077
      if(ipen.eq.0) then
         write(kunit,111) xp,yp
 111     format('s',' ',f6.1,' ',f6.1,' m')
      else
         write(kunit,112) xp,yp
 112     format(f6.1,' ',f6.1,' d')
      end if
 113  continue
 1111 continue
      return
      end

      subroutine linew(width)
c
c     change the postscript line width (1 = 1/72 inch)
c
      common/punit/kunit
      if(kunit.eq.0) goto 1111
      write(kunit,111) width
 111  format('s',' ',f6.3,' w')
 1111 continue
      return
      end


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
      character*1 leadin
      leadin=char(29)
      call pointqq(leadin,1)
      return
      end

      subroutine gmode
      character*1 leadin
      leadin=char(29)
      call pointqq(leadin,1)
      return
      end

      subroutine area(xa,xb,ya,yb)
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
      end

      subroutine tmode
      character*1 leadout
      leadout=char(31)
      call pointqq(leadout,1)
c      call flush(6)
      return
      end

      subroutine penup
      character*1 pencom
      pencom=char(29)
      call pointqq(pencom,1)
      return
      end

      subroutine erase
      character*1 ers(4)
      ers(1)=char(31)
      ers(2)=char(27)
      ers(3)=char(12)
      ers(4)=char(29)
      call pointqq(ers,4)
      return
      end

      subroutine plott(x,y,ipen)
c
c     draws line to absolute coordinate (x,y)
c     moves to absolute coordinate (x,y) if penup
c     was called just prior to plot call (the next plot call
c     will draw line). Absolute coordinates are defined as:
c     0.00 < x < 1.313, 0.00 < y < 1.00, with the maximum plotting area
c     0.25 < x < 1.250, 0.15 < y < 0.75.
c
      character*1 coord(6)
      integer hx0,hy0,lx0,ly0,llx0,lly0
      integer hx,hy,lx,ly
      common/penplot/x0,x1,y0,y1,hy0,ly0,llx0,lly0,hx0,lx0
          ix=(x-x0)/(x1-x0)*4095+.4999999
          iy=(y-y0)/(y1-y0)*3119+.4999999
      if(ix.lt.0 .or. ix.gt.4095 .or. iy.lt.0 .or. iy.gt.3119) then
            call penup
            return
      end if
      if(ipen.eq.0) call penup
      hy=iy/128
      lly=iy-128*hy
        ly=lly/4
        lly=lly-4*ly
      hx=ix/128
      llx=ix-128*hx
        lx=llx/4
        llx=llx-4*lx
      i=0
      if(hy.ne.hy0.or.hx.ne.hx0) then
            i=i+1
            coord(i)=char(hy+32)
      end if
        if(lly.ne.lly0 .or.llx.ne.llx0) then
           i=i+1
           coord(i)=char(96+4*lly+llx)
           i=i+1
           coord(i)=char(ly+96)
         elseif(ly.ne.ly0 .or. hx.ne.hx0.or.hy.ne.hy0) then
            i=i+1
            coord(i)=char(ly+96)
      end if
      if(hx.ne.hx0.or.hy.ne.hy0) then
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
      end

      subroutine pointqq(c,n)
      character*1 c(n)
      common /unit/iunit
      if(iunit.eq.0) goto 1111
      do 10 i=1,n
 10      write(iunit,100) c(i)
 100     format(a,$)
c      do 10 i=1,n
c 10      j=fputc(iunit,c(i))
c      call flush(iunit)
 1111 continue
      return
      end

c*****************************************************************
c
c     contour plotting routines given to me by researchers at NASA Ames
c     in which a number of changes has been made.
c
c*****************************************************************

      subroutine contxx(idim,jdim,is,ie,js,je,xy,f,ncopt,ncont,acont,
     &nxdim,xcont,ycont,naddim,nad,nlev,iadim,ia)
c
c   calculate contour lines for the function f in the region is to ie, js to
c   je.  x,y coordinates corresponding to the grid points are in array xy.
c
c   if ncopt=1, figure our own contour levels, up to ncont of them, using
c   "nice" numbers.  frange finds the function range in the given region, and
c   consc1 computes the contour levels.  note that ncont will be revised
c   downward to correspond to the number of contour levels actually used.
c   the contour levels calculated are returned in the array acont.
c   NEW: zero contour will be removed by shifting half a contour spacing
c
c   if ncopt=2, calculate lines for the ncont contour levels specified in
c   acont.
c
c   the contour lines are returned in arrays xcont, and ycont.  nad(1)
c   gives the number of contour lines, and nad(n) points to the start of the
c   nth line (i.e. nad(n+1) points to one past the end of the nth line).
c   nlev(n) returns the contour level of the nth contour line.
c
c   ia is a scratch array.  try a dimension of 3000.
c
      dimension xy(idim,jdim,2),f(idim,jdim)
      dimension acont(ncont)
      dimension xcont(nxdim),ycont(nxdim)
      dimension nad(naddim),nlev(naddim),ia(iadim)
c
c   if ncopt=1, figure our own contour levels, up to ncont of them, using
c   "nice" numbers.  frange finds the function range in the given region, and
c   consc1 computes the contour levels.  note that ncont will be revised
c   downward to correspond to the number of contour levels actually used.
c
      if (ncopt.eq.1) then
       call frange(idim,jdim,is,ie,js,je,f,fmin,fmax)
       call consc1(fmin,fmax,ncont,acont)
       if(acont(1)*acont(ncont).le.0.) then
          aincre=acont(2)-acont(1)
          do 1739 ii=1,ncont
             acont(ii)=acont(ii)+aincre/2.
1739      continue
          if(acont(ncont).gt.fmax) ncont=ncont-1
          if(acont(1)-aincre.gt.fmin) then
            ncont=ncont+1
            do 1738 ii=ncont,2,-1
               acont(ii)=acont(ii-1)
1738        continue
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
c   one little check.  if is=ie or js=je, return with no contour lines.
c
      if (is.eq.ie .or. js.eq.je) goto 110
c
c   loop through each contour level.
c
      do 100 icont= 1,ncont
       za= acont(icont)
      m=0
c ****  scan points and determine points of ia
      do  600  j= js+1,je-1
      imb=0
      do  600  i= is,ie
      if(f(i,j).le.za)  go  to  601
      if(imb.ne.1)  go  to  600
      m=m+1
      if(m.gt.iadim)  go  to  210
      ia(m)=1000*i+j
      imb=0
      go  to  600
 601  imb=1
 600  continue
c****  search start point on boundary line
 101  ima=1
      imb=0
      ixa= is-1
      iya= js
 1    ixa=ixa+1
      if(ixa.eq.ie)  ima=2
      go  to  5
 2    iya=iya+1
      if(iya.eq.je)  ima=3
      go  to  5
 3    ixa=ixa-1
      if(ixa.eq.is)  ima=4
      go  to  5
 4    iya=iya-1
      if(iya.eq.js)  ima=5
 5    if(f(ixa,iya).gt.za)  go  to  7
      imb=1
 6    go  to  (1,2,3,4,91),ima
 7    if(imb.ne.1)  go  to  6
c****  determine start point
      imb=0
      ix=ixa
      iy=iya
      s=f(ixa,iya)
      go  to  (21,11,12,13,51),ima
 11   if(iy.ne.js)  go  to  31
      go  to  21
 12   if(ix.ne.ie)  go  to  41
      go  to  31
 13   if(iy.ne.je)  go  to  51
      go  to  41
 10   ix=ia(n)/1000
      iy=ia(n)-1000*ix
      s=f(ix,iy)
      ia(n)=0
      go  to  21
c****  process to search plot point
 20   iy=iy+1
 21   ix=ix-1
      if(ix.lt.is)  go  to  90
      i=1
      if(f(ix,iy).le.za)  go  to  52
      s=f(ix,iy)
      go  to  31
 30   ix=ix-1
 31   iy=iy-1
      if(iy.lt.js)  go  to  90
      i=2
      if(f(ix,iy).le.za)  go  to  60
      s=f(ix,iy)
      go  to  41
 40   iy=iy-1
 41   ix=ix+1
      if(ix.gt.ie)  go  to  90
      i=3
      if(f(ix,iy).le.za)  go  to  60
      s=f(ix,iy)
      go  to  51
 50   ix=ix+1
 51   iy=iy+1
      i=4
      if(iy.gt.je)  go  to  90
      if(f(ix,iy).le.za)  go  to  60
      s=f(ix,iy)
      go  to  21
 52   if(m.eq.0)  go  to  60
      ik=1000*ix+iy+1000
      do  602  j=1,m
      if(ia(j).ne.ik)  go  to  602
      ia(j)=0
 602  continue
c****  calculate plot point
 60   xyf=(za-f(ix,iy))/(s-f(ix,iy))
      go  to  (61,62,63,64),i
 61   wxx= xy(ix,iy,1)+xyf*(xy(ix+1,iy,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix+1,iy,2)-xy(ix,iy,2))
      go  to  65
 62   wxx= xy(ix,iy,1)+xyf*(xy(ix,iy+1,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix,iy+1,2)-xy(ix,iy,2))
      go  to  65
 63   wxx= xy(ix,iy,1)+xyf*(xy(ix-1,iy,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix-1,iy,2)-xy(ix,iy,2))
      go  to  65
 64   wxx= xy(ix,iy,1)+xyf*(xy(ix,iy-1,1)-xy(ix,iy,1))
      wyy= xy(ix,iy,2)+xyf*(xy(ix,iy-1,2)-xy(ix,iy,2))
c****  plot
65    continue
c**** decide if plot point equal initial plot point
      if(iw.ne.3)  go  to  66
      np=1
      nad(nlinep)= nad(nlinep-1)
      nlev(nlinep-1)= icont
      npt=0
      xcont(nad(nlinep))= wxx
      ycont(nad(nlinep))= wyy
      wx=wxx
      wy=wyy
      iw=2
      go  to  67
66    continue
      nad(nlinep)= nad(nlinep)+1
      if (nad(nlinep).gt.nxdim) goto 220
      np=np+1
      xcont(nad(nlinep))= wxx
      ycont(nad(nlinep))= wyy
c     if (np .lt. 200) go to 6602
c
c     call draw2d(pt,np,2,2,0)
c     npt=npt+np
c     np=1
c     pt(1,1)= wxx
c     pt(2,1)= wyy
6602  if(wxx.ne.wx)  go  to  67
      if(wyy.eq.wy)  go  to  90
c****  determine next process
 67   go  to  (50,20,30,40),i
 90   iw=3
      nad(nlinep)= nad(nlinep)+1
      if (nad(nlinep).gt.nxdim) goto 220
      if (np.gt.1) then
       nlinep= nlinep+1
       if (nlinep.gt.naddim) goto 230
      end if
c     if (np .gt. 1) call draw2d(pt,np,2,2,0)
      if(ima.ne.5)  go  to  6
c****  search start point
      if(m.eq.0)  go  to  92
 91   do  603  n=1,m
      if(ia(n).ne.0)  go  to  10
 603  continue
 92   continue
c****  calculate value of next curve
 100  continue
c
  110 continue
      nad(1)= nlinep-2
      return
c
c   warning - ia array full.
c
  210 continue
      write(6,211) iadim
  211 format('  warning - scratch array ia full in contour routine ',
     c 'contxx.'/
     c       '  picture may be incomplete.  array was dimensioned ',
     c i5,'.')
      goto 110
c
c   warning - xcont array full.
c
  220 continue
      write(6,221) nxdim
  221 format('  warning - contour line array xcont full in contour ',
     c 'routine contxx.'/
     c       '  picture may be incomplete.  array was dimensioned ',
     c i5,'.')
      goto 110
c
c   warning - nad array full.
c
  230 continue
      write(6,231) naddim
  231 format('  warning - contour line pointer array nad full in ',
     c 'contour routine contxx.'/
     c       '  picture may be incomplete.  array was dimensioned ',
     c i5,'.')
      goto 110
      end
      subroutine frange(idim,jdim,is,ie,js,je,f,fmin,fmax)
c
c   find the range (minimum and maximum) of the function f in the region is to
c   ie, js to je.
c
      dimension f(idim,jdim)
c
      fmin= f(is,js)
      fmax= fmin
      do 10 j= js,je
       do 10 i= is,ie
          fmin= min(fmin,f(i,j))
          fmax= max(fmax,f(i,j))
   10       continue
      return
      end
      subroutine consc1(amin,amax,ncont,acont)
c
c   come up with a "nice" scaling of about ncont values between amin and amax.
c   ncont is updated to the number of intervals actually needed.
c
      dimension acont(ncont)
      dimension rnice(4)
      data rnice/.1,.2,.25,.5/
      data nnice/4/
c
c   as a first approximation, get the difference, its characteristic and
c   mantissa.
c
      diff= (amax-amin)/(ncont+1)
      if (diff.le.0.) goto 20
      char= log10(diff)+1.
c
c   round char down and get the mantissa.
c
      if (char.ge.0.) then
       ichar= char
      else
       ichar= char-1.
      end if
      rmant= diff*10.**(-ichar)
c
c   what's the next largest "nice" mantissa?
c
      do 3 i= 1,nnice
       if (rmant.le.rnice(i)) goto 10
3      continue
      i= nnice
c
c   got a guess.  calculate a diff, round amin down.
c
   10 continue
      ainc= rnice(i)*10.**ichar
      imin= amin/ainc
      if (amin.gt.imin*ainc) imin= imin+1
      imax= amax/ainc
      if (amax.lt.imax*ainc) imax= imax-1
      nneed= imax+1-imin
c
c   are we under?
c
      if (nneed.gt.ncont) then
c
c   nope.  try the next nice number.
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
c   now just set up the acont array and update ncont.  (set up all (original)
c   ncont of the aconts for the folks back home.)
c
      do 1 i= 1,ncont
       acont(i)= (imin-1+i)*ainc
1      continue
      ncont= nneed
      goto 30
c
c   all values are the same -- just set up one contour level.
c
   20 continue
      acont(1)= amin
      do 2 i= 2,ncont
         acont(i)= 0.
2     continue
      ncont= 1
c
   30 continue
      return
      end
