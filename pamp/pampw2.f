c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pampw2
c
c     Plots a  number of components from an wave-amplitude file
c     assumes new format and handles both with and without y-dependence
c
      implicit none

      integer imax,maxwave,nypm
      parameter(imax=10000)
      parameter(nypm=200,maxwave=18)
      character*80 xaxis,yaxis,heada,headb,headc
      character*32 namamp(5)
      character*40 cwave
c      character*7 mode(maxwave)
      character*6 dash(5)
      character*1 comp(3)
      integer npoint(5),kx(maxwave),kz(maxwave),y
      integer iw(maxwave),i,k,j,ndim,mbox,nplot,l
      integer iopt,iplot,nyp,mwave,nn(5),idash
      integer ntot,ni(5)
      real t(imax,5),v(imax,maxwave),u(4,maxwave,imax,5)
      real alfa0,beta0,k2,pi,re
      real ymin,ymax,xmin,xmax
      complex dum
      logical longlist
c
c     Constants
c
      parameter (pi = 3.1415926535897932385)
c
c     Plotting
c
      real fnorm
      integer inorm
c
      data kx/maxwave*0/,kz/maxwave*0/
      data dash/' solid',', dash',',  dot',', chdh',', chdt'/

      comp(1)='u'
      comp(2)='v'
      comp(3)='w'

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         pampw2 $Rev$'//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Give number of wave files (max 5)'
      read(*,*) nplot
      write(*,*) 'Give the names of the wave files, end with ,CR>'
      do k=1,nplot
         read(*,222) namamp(k)
 222     format(a)
      end do

      write(*,*)
c
c     Read files namamp
c
      do k=1,nplot
         nn(k)=index(namamp(k),' ')
         write(*,222) namamp(k)(1:nn(k)-1)

         open(unit=10,file=namamp(k),status='old')
         do i=1,imax
            read(10,*,end=6) t(i,k),mwave,nyp,re,alfa0,beta0,longlist
            do l=1,mwave
               read(10,*,end=6) kx(l),kz(l),(u(j,l,i,k),j=1,4)
            end do

            if (longlist) then
               do j=1,mwave
                  read(10,*,end=6) (dum,y=1,nyp,2)
                  read(10,*,end=6) (dum,y=1,nyp,2)
                  read(10,*,end=6) (dum,y=1,nyp,2)
                  read(10,*,end=6) (dum,y=1,nyp,2)
               end do
            end if
         end do

         write(*,*) 'The declared maximum is reached: ',i-1
         i=imax+1
 6       ni(k)=i-1
         npoint(k)=ni(k)
         close(unit=10)
      end do

      headb='From files '//namamp(1)(1:nn(1)-1)
      headc=dash(1)
      ntot=nn(1)+10
      do k=2,nplot
         headb=headb(1:ntot)//', '//namamp(k)(1:nn(k)-1)
         headc=headc(1:(k-1)*6)//dash(k)
         ntot=ntot+nn(k)+1
      end do
      headb=headb(1:ntot)//'$'
      headc=headc(1:nplot*6)//'$'

      write(*,*) '1 = Norm with init., 2 = No norm'
      read(*,*) inorm

 1    write(*,*) '0 = End Program'
      do i=1,mwave
         write(*,431) i,kx(i),kz(i)
      end do
 431  format(3i3)
      write(*,*) ' '

      write(*,*) ' Choose iwave '
      read(*,*) iw(1)
      if (iw(1).eq.0) goto 9999

      write(cwave,432) kx(iw(1)),kz(iw(1)),alfa0*kx(iw(1))
     &     ,beta0*kz(iw(1))
 432  format('kx=',i2,', kz=',i2,'; alpha=',f7.3,' beta=',f7.3)

      write(*,*) ' Negative value for more options'
      write(*,*) ' 1 = Energy '
      write(*,*) ' 2 = Energy contribution from u'
      write(*,*) ' 3 = Energy contribution from v'
      write(*,*) ' 4 = Energy contribution from w'
      write(*,*) ' 5 = Energy contribution from Dv'
      write(*,*) ' 6 = Energy contribution from eta'
      write(*,*) ' 7 = Change inorm'
      read(*,*) iplot
      fnorm=1.
      iopt=1
      if (iplot.lt.0) iopt=-1
      iplot=abs(iplot)

      if (iplot.eq.1) then
         do k=1,nplot
            if (inorm.eq.1) then
               fnorm=0.5*(u(1,iw(1),1,k)**2+u(2,iw(1),1,k)**2
     &              +u(3,iw(1),1,k)**2)
               if (fnorm.lt.1e-17) fnorm=1.
               write(*,*) 'fnorm= ',fnorm
            end if
            do i=1,ni(k)
               v(i,k)=0.5*(u(1,iw(1),i,k)**2+u(2,iw(1),i,k)**2
     &              +u(3,iw(1),i,k)**2)/fnorm
            end do
         end do
         yaxis='energy$'
         heada='Energy for '//cwave//'$'
      end if

      if (iplot.ge.2.and.iplot.le.4) then
         do k=1,nplot
            if (inorm.eq.1) then
               fnorm=0.5*(u(iplot-1,iw(1),1,k)**2)
               if (fnorm.lt.1e-17) fnorm=1.
               write(*,*) 'fnorm= ',fnorm
            end if
            do i=1,ni(k)
               v(i,k)=0.5*(u(iplot-1,iw(1),i,k)**2)/fnorm
            end do
         end do
         yaxis='energy$'
         heada='Energy in '//comp(iplot-1)// ' for '//cwave//'$'
      end if

      if (iplot.eq.5) then
         k2=(kx(iw(1))*alfa0)**2+(kz(iw(1))*beta0)**2
         if (k2.eq.0.0) goto 1
         do k=1,nplot
            if (inorm.eq.1) then
               fnorm=0.5*(u(1,iw(1),1,k)**2+u(2,iw(1),1,k)**2
     &              +u(3,iw(1),1,k)**2)/fnorm
     &              -0.5*(u(2,iw(1),1,k)**2+u(4,iw(1),1,k)**2/k2)/fnorm
               if (fnorm.lt.1e-17) fnorm=1.
            end if
            do i=1,ni(k)
               v(i,k)=0.5*(u(1,iw(1),i,k)**2+u(2,iw(1),i,k)**2
     &              +u(3,iw(1),i,k)**2)
     &              -0.5*(u(2,iw(1),i,k)**2+u(4,iw(1),i,k)**2/k2)
            end do
         end do
         yaxis='energy$'
         heada='Energy in Dv for '//cwave//'$'
      end if

      if (iplot.eq.6) then
         k2=(kx(iw(1))*alfa0)**2+(kz(iw(1))*beta0)**2
         if (k2.eq.0.0) goto 1
         do k=1,nplot
            if (inorm.eq.1) then
               fnorm=0.5*u(4,iw(1),1,k)**2/k2
               if (fnorm.lt.1e-17) fnorm=1.
            end if
            do i=1,ni(k)
               v(i,k)=0.5*u(4,iw(1),i,k)**2/k2/fnorm
            end do
         end do
         yaxis='energy$'
         heada='Energy in eta for '//cwave//'$'
      end if

      if (iplot.eq.7) then
         write(*,*) '1 = Norm with init., 2 = No norm'
         read(*,*) inorm
      end if
c
c     Plotting
c
      ndim=imax
      idash=1
      mbox=1
c      xmin=t(1,1)
c      xmax=t(n,1)
      xmin=1.
      xmax=1.
      ymin=0.
      ymax=0.
      if (iopt.eq.-1) then
         write(*,*) 'Give xmin, xmax, ymin, ymax'
         read(*,*) xmin,xmax,ymin,ymax
      end if
      xaxis='t$'
      call rita1(t,v,xmin,xmax,ymin,ymax,npoint,nplot,ndim,
     &        xaxis,yaxis,heada,headb,headc,mbox,idash)

      goto 1

 9999 stop

      end program pampw2
