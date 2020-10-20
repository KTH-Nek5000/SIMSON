c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pampw
c
c     Plots a number of components from an wave-amplitude file
c     assumes new format and handles both with and without y-dependence
c
      implicit none

      integer i,j,k,nplot,iplot,iopt,idash,mbox,ndim
      integer imax,nwave,nypm,kxw,kzw,it1,nyp,iw,n,nn,mwave
      real xmin,xmax,ymin,ymax,dt,re,t1
      real pi
c
c     Constants
c
      parameter (pi = 3.1415926535897932385)
      parameter(imax=1000)
      parameter(nypm=97,nwave=36)

      character*80 xaxis,yaxis,heada,headb,headc
      character*32 namamp
      character*40 cwave
      integer npoint(4),kx(nwave),kz(nwave),y,npoint1(2)
      real t(imax,4),v(imax,4),u(4,nwave,imax)
      real alfa0,beta0,k2,eta(nypm,3)
      complex campw(nypm,2,nwave,imax)
      logical longlist

      data kx/nwave*0/,kz/nwave*0/

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         pampw $Rev$ '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Give the name of the wave file to be plotted'
      read(*,222) namamp
 222  format(a)
c
c     Read file namamp
c
      open(unit=10,file=namamp,status='old')
      do i=1,imax
        read(10,*,end=6) t(i,1),mwave,nyp,re,alfa0,beta0,longlist
        do k=1,mwave
           read(10,*,end=6) kx(k),kz(k),(u(j,k,i),j=1,4)
        end do
        t(i,2)=t(i,1)
        t(i,3)=t(i,1)
        t(i,4)=t(i,1)
        if (longlist) then
           do k=1,mwave
              read(10,*,end=6) (campw(y,1,k,i),y=1,nyp)
              read(10,*,end=6) (campw(y,2,k,i),y=1,nyp)
           end do
        end if
      end do
      write(*,*) 'The declared maximum is reached: ',i-1
      i=imax+1
 6    n=i-1
      close(unit=10)
      do j=1,4
         npoint(j)=n
      end do
      npoint1(1)=nyp
      npoint1(2)=nyp
      nn=index(namamp,' ')
      headb='from file '//namamp(1:nn)//'$'
      do y=1,nyp
         do j=1,3
            eta(y,j)=-cos(pi*(y-1)/float(nyp-1))
         end do
      end do

 1    do i=1,mwave
         write(*,431) i,kx(i),kz(i)
      end do
 431  format(3i3)
      write(*,*) ' Choose iwave'
      read(*,*) iw
      kxw=kx(iw)
      kzw=kz(iw)
      write(cwave,432) kxw,kzw,alfa0*kxw,beta0*kzw
 432  format('kx=',i2,', kz=',i2,'; alpha=',f7.3,' beta=',f7.3)

      write(*,*) ' Negative value for more options'
      write(*,*) ' 1 = energy + contributions from u,v,w'
      write(*,*) ' 2 = energy + contributions from v,Dv,eta'
      write(*,*) ' 3 = rms of u,v,w'
      if (longlist) write(*,*) ' 4 = y-dependence at specified time'
      read(*,*) iplot

      iopt=1
      if (iplot.lt.0) iopt=-1
      iplot=abs(iplot)

      if (iplot.eq.1) then
         do i=1,n
            v(i,1)=0.5*(u(1,iw,i)**2+u(2,iw,i)**2+u(3,iw,i)**2)
            v(i,2)=0.5*u(1,iw,i)**2
            v(i,3)=0.5*u(2,iw,i)**2
            v(i,4)=0.5*u(3,iw,i)**2
            write(53,*) t(i,1),v(i,1)
         end do
         yaxis='energy$'
         heada='Energy for '//cwave//'$'
         headc='solid: energy, dash: u^2/2, dot: v^2/2, chdh: w^2/2$'
         nplot=4
      end if

      if (iplot.eq.2) then
         k2=(kxw*alfa0)**2+(kzw*beta0)**2
         if (k2.eq.0.0) goto 1
         do i=1,n
            v(i,1)=0.5*(u(1,iw,i)**2+u(2,iw,i)**2+u(3,iw,i)**2)
            v(i,2)=0.5*u(2,iw,i)**2
            v(i,3)=v(i,1)-0.5*(u(2,iw,i)**2+u(4,iw,i)**2/k2)
            v(i,4)=0.5*u(4,iw,i)**2/k2
         end do
         yaxis='energy$'
         heada='Energy for '//cwave//'$'
         headc='solid: energy, dash: v^2/2, dot: Dv^2/2, chdh: eta^2/2$'
         nplot=4
      end if

      if (iplot.eq.3) then
         do i=1,n
            v(i,1)=u(1,iw,i)
            v(i,2)=u(2,iw,i)
            v(i,3)=u(3,iw,i)
         end do
         yaxis='rms$'
         heada='Rms velocities for '//cwave//'$'
         headc='solid: u, dash: v, dot: w$'
         nplot=3
      end if

      if (iplot.eq.4) then
         write(*,*) 'Give time, less than tmax: ',t(n,1)
         read(*,*) t1
         dt=t(2,1)-t(1,1)
         it1=int(t1/dt+.01)
         t1=it1*dt
         write(*,*) 'Time choosen: ',t1,it1
         pause
         do y=1,nyp
            v(y,1)=abs(campw(nyp-y+1,1,iw,it1))
            v(y,2)=abs(campw(nyp-y+1,2,iw,it1))
            write(*,*) v(y,1),v(y,2)
         end do
         yaxis='y$'
         heada='Rms vel/vor for '//cwave//'$'
         headc='solid: v, dash: eta,$'
         if (kxw.eq.0.and.kzw.eq.0) headc='solid: u, dash: w,$'
         nplot=2
      end if

      ndim=imax
      idash=1
      mbox=1
      if (iplot.le.3) then
         xmin=t(1,1)
         xmax=t(n,1)
         ymin=0.
         ymax=0.
         if (iopt.eq.-1) then
            write(*,*) 'Give xmin, xmax, ymin, ymax'
            read(*,*) xmin,xmax,ymin,ymax
         end if
         xaxis='t$'
         call rita1(t,v,xmin,xmax,ymin,ymax,npoint,nplot,ndim,
     &        xaxis,yaxis,heada,headb,headc,mbox,idash)
      end if
      if (iplot.eq.4) then
         xmin=0.
         xmax=0.
         ymin=eta(1,1)
         ymax=eta(nyp,1)
         if (iopt.eq.-1) then
            write(*,*) 'Give xmin, xmax, ymin, ymax'
            read(*,*) xmin,xmax,ymin,ymax
         end if
         xaxis='amplitude$'
         call rita1(v,eta,xmin,xmax,ymin,ymax,npoint1,nplot,ndim,
     &        xaxis,yaxis,heada,headb,headc,mbox,idash)
      end if
      goto 1

      end program pampw
