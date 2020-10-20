c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pext1
c
c     Plots a  number of components from an extremum file
c
      integer iold,i,j,k,n,nn,iplot,iopt,ivar,nplot,ndim
      integer idash,mbox
      real xmin,xmax,ymin,ymax
      character*80 xaxis,heada,headb,headc
      character*80 namext
      character*10 var
      character*20 yaxis
      integer npoint(5),y,nyp
      real t(2000,5),u(2000,2,7),c(2000,2,7,3),v(2000,5)
      real dum
      integer ifirst,icoor

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         pext1 $Rev$ '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Give the name of the extremum file to be plotted'
      read(*,222) namext
 222  format(a)
      write(*,*) 'Ext format: 1=normal, 2=long'
      read(*,*) iold
      if (iold.eq.2) then
         write(*,*) 'Give nyp'
         read(*,*) nyp
      end if
c
c     Read file namext
c
      open(unit=10,file=namext,status='old')
      do i=1,2000
         read(10,*,end=4,err=4) t(i,1)
         do k=1,2
            do j=1,7
               read(10,*,end=4,err=4) u(i,k,j),c(i,k,j,1)
               read(10,*,end=4,err=4) c(i,k,j,2),c(i,k,j,3)
            end do
         end do
         if (iold.eq.2) then
            do j=1,12*nyp
               read(10,*,end=4,err=4) (dum,y=1,2)
               read(10,*,end=4,err=4) (dum,y=1,2)
            end do
         end if
         do j=2,5
            t(i,j)=t(i,1)
         end do
      end do
      write(*,*) 'The file is too large. Read 2000 time steps.'
      i=2001
 4    n=i-1
      close(unit=10)
      nn=index(namext,' ')
      headb='from file '//namext(1:nn)//'$'

 1    write(*,*) ' Negative value for more options'
      write(*,*) ' 0 = End program'
      write(*,*) ' 1 = Min of u,v,w'
      write(*,*) ' 2 = Max of u,v,w'
      write(*,*) ' 3 = Min of omx,omy,omz'
      write(*,*) ' 4 = Max of omx,omy,omz'
      write(*,*) ' 5 = Min and max, one velocity or vorticity'
      write(*,*) ' 6 = Coordinates of min of one vel or vor'
      write(*,*) ' 7 = Coordinates of max of one vel or vor'
      write(*,*) ' 8 = One coordinate of max and min of one vel or vor'
      read(*,*) iplot

      if (iplot.eq.0) stop
      iopt=1
      if (iplot.lt.0) iopt=-1
      iplot=abs(iplot)
      if (iplot.ge.5) then
         write(*,*) 'Which variable'
         write(*,*) '1 u-ulam'
         write(*,*) '2 v'
         write(*,*) '3 w'
         write(*,*) '4 omx'
         write(*,*) '5 omy'
         write(*,*) '6 omz-omzlam'
         write(*,*) '7 omz'
         read(*,*) ivar
         if (ivar.eq.1) var='u-ulam'
         if (ivar.eq.2) var='v'
         if (ivar.eq.3) var='w'
         if (ivar.eq.4) var='omx'
         if (ivar.eq.5) var='omy'
         if (ivar.eq.6) var='omz-omzlam'
         if (ivar.eq.7) var='omz'
      end if

      if (iplot.ge.1.and.iplot.le.4) then
         do j=1,3
            do i=1,n
               v(i,j)=u(i,1+mod(iplot-1,2),j+3*(iplot/3))
            end do
         end do
         if (iplot.eq.1) yaxis='min velocity$'
         if (iplot.eq.2) yaxis='max velocity$'
         if (iplot.eq.3) yaxis='min vorticity$'
         if (iplot.eq.4) yaxis='max vorticity$'
         heada=yaxis
         if (iplot.le.2) headc='solid: u, dash: v, dot: w$'
         if (iplot.ge.3) headc='solid: omx, dash: omy, dot: omz-omzlam$'
         nplot=3
         ifirst=1
      end if

      if (iplot.eq.5) then
         do j=1,2
            do i=1,n
               v(i,j)=u(i,j,ivar)
            end do
         end do
         yaxis=var//'$'
         heada='min and max of'//var//'$'
         headc='solid: min, dash: max $'
         nplot=2
         ifirst=1
      end if
      if (iplot.eq.6.or.iplot.eq.7) then
         write(*,*) 'Which point to start with (1 is first)'
         read(*,*) ifirst
         do j=1,3,2
            do i=ifirst,n
               v(i+1-ifirst,j)=c(i,iplot-5,ivar,j)
            end do
         end do
         do i=ifirst,n
            v(i+1-ifirst,2)=c(i,iplot-5,ivar,2)*10.
         end do
         yaxis='x,10*y,z$'
         if (iplot.eq.6) heada='coordinates of minimum for '//var//'$'
         if (iplot.eq.7) heada='coordinates of maximum for '//var//'$'
         headc='solid: x, dash: 10*y, dot: z$'
         nplot=3
      end if
      if (iplot.eq.8) then
         write(*,*) 'Which coordinate ?'
         write(*,*) ' 1 x'
         write(*,*) ' 2 y'
         write(*,*) ' 3 z'
         read(*,*) icoor
         do i=1,n
            v(i,1)=c(i,1,ivar,icoor)
            v(i,2)=c(i,2,ivar,icoor)
         end do
         if (icoor.eq.1) then
            yaxis='x'
            heada='x-coordinate of min and max for '//var//'$'
         end if
         if (icoor.eq.2) then
            yaxis='y'
            heada='y-coordinate of min and max for '//var//'$'
         end if
         if (icoor.eq.3) then
            yaxis='z'
            heada='z-coordinate of min and max for '//var//'$'
         end if
         headc='solid: min, dash: max'
         nplot=2
         ifirst=1
      end if

      xmin=t(1,1)
      xmax=t(n+1-ifirst,1)
      do j=1,5
         npoint(j)=n+1-ifirst
      end do
      ymin=0.
      ymax=0.
      if (iopt.eq.-1) then
         write(*,*) 'Give xmin, xmax, ymin, ymax'
         read(*,*) xmin,xmax,ymin,ymax
      end if
      xaxis='t$'
      ndim=2000
      idash=1
      mbox=1

      call rita1(t(ifirst,1),v,xmin,xmax,ymin,ymax,npoint,nplot,ndim,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)
      goto 1

      end program pext1
