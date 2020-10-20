      subroutine savelist(iplot,
     &                x,y,xmin1,xmax1,ymin1,ymax1,npoint,nplot,ndim,
     &                xaxis,yaxis,heada,headb,headc,ilog,idash)
c
c     iplot    1  ritlog
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
      implicit none

      integer iplot
      real x(ndim,1),y(ndim,1)
      integer npoint(1),nplot,ndim,ilog,idash
      real xmin1, xmax1,ymin1,ymax1
      character*(*) xaxis,yaxis,heada,headb,headc

      integer i,j

      open(unit=40,file='plot.dat',status='unknown',
     &     form='formatted')
      write(40,'(3I5)') npoint,nplot,ndim
      write(40,'(4E18.9)') xmin1,xmax1,ymin1,ymax1
      write(40,'(3I5)') iplot,ilog,idash
      write(40,'(a80)') heada
      write(40,'(a80)') headb
      write(40,'(a80)') headc
      write(40,'(a80)') xaxis
      write(40,'(a80)') yaxis

      do i=1,npoint(1)
         write(40,'(5E18.9)')
     &        x(i,1),((y(i+(j-1)*ndim,1)),j=1,nplot)
      end do

      close(40)

      return

      end subroutine savelist
