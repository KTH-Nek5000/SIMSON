c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pamp1
c
c     Plots a  number of components from an amplitude file.
c     Works with both old and new amplitude files
c
      implicit none

      integer i,j,nplot,iplot,iopt,idash,mbox,ndim
      integer iold,n,nn
      real xmin,xmax,ymin,ymax,energy
      character*80 xaxis,yaxis,heada,headb,headc
      character*32 namamp
      integer npoint(5),y,nyp
      real t(10000,5),v(10000,5),u(10000,10)
      real dum

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         pamp1 $Rev$ '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Give the name of the amplitude file to be plotted'
      read(*,222) namamp
 222  format(a)
c      write(*,*) 'amp format: -1=real old(6), 0=old(7),
c     &                         1=new(10), 2=new long(10+nyp)'
      write(*,*) 'Amplitude file format: 1=short(11), 2=long(11+nyp)'
      read(*,*) iold
      if (iold.eq.2) then
        write(*,*) 'Give nyp'
        read(*,*) nyp
      end if
c
c     Read file namamp
c
      open(unit=10,file=namamp,status='old')
      do i=1,10000
         read(10,*,end=4) t(i,1),(u(i,j),j=1,10)
c         if (iold.eq.-1) read(10,*,end=4) (u(i,j),j=4,6)
c         if (iold.ge.0)  read(10,*,end=4) (u(i,j),j=4,7)
c         if (iold.ge.1)  read(10,*,end=4) (u(i,j),j=8,10)
         if (iold.eq.2) then
            do j=1,10
               read(10,*,end=4) (dum,y=1,nyp)
            end do
         end if
         do j=2,5
            t(i,j)=t(i,1)
         end do
      end do
      i=10001
 4    n=i-1
      close(unit=10)
      do j=1,5
         npoint(j)=n
      end do
      nn=index(namamp,' ')
      headb='From file '//namamp(1:nn)//'$'

 1    write(*,*) ' Negative value for more options'
      if (iold.le.0) write(*,*) ' 1 = Energy'
      if (iold.ge.1) write(*,*) ' 1 = Energy, with and without mean'
      write(*,*) ' 2 = Energy + contributions from u,v,w'
      write(*,*) ' 3 = Energy + contributions from v,Dv,eta'
      write(*,*) ' 4 = % of energy in u,v,w,Dv,eta'
      write(*,*) ' 5 = rms of u,v,w'
      write(*,*) ' 6 = rms of xi,eta,zeta'
      if (iold.ge.1) then
         write(*,*) ' 7 = dUuv (averaged Reynolds stress)'
         write(*,*) ' 8 = h+ (Reynolds number based on u_tau'
      end if
      read(*,*) iplot

      iopt=1
      if (iplot.lt.0) iopt=-1
      iplot=abs(iplot)

      if (iplot.eq.1) then
         do i=1,n
            v(i,1)=0.5*(u(i,1)**2+u(i,2)**2+u(i,3)**2)
            if (iold.ge.1) v(i,2)=v(i,1)-u(i,9)**2*0.5
         end do
         yaxis='energy$'
         heada='Energy$'
         headc=' $'
         nplot=1
         if (iold.ge.1) nplot=2
      end if
c
      if (iplot.eq.2) then
         do i=1,n
            v(i,1)=0.5*(u(i,1)**2+u(i,2)**2+u(i,3)**2)
            v(i,2)=0.5*u(i,1)**2
            v(i,3)=0.5*u(i,2)**2
            v(i,4)=0.5*u(i,3)**2
         end do
         yaxis='energy$'
         heada='Energy$'
         headc='solid: energy, dash: u^2/2, dot: v^2/2, chdh: w^2/2$'
         nplot=4
      end if

      if (iplot.eq.3) then
         if (iold.eq.-1) goto 1
         do i=1,n
            v(i,1)=0.5*(u(i,1)**2+u(i,2)**2+u(i,3)**2)
            v(i,2)=0.5*u(i,2)**2
            v(i,3)=v(i,1)-0.5*(u(i,2)**2+u(i,7)**2)
            v(i,4)=0.5*u(i,7)**2
         end do
         yaxis='energy$'
         heada='Energy$'
         headc='solid: energy, dash: v^2/2, dot: Dv^2/2, chdh: eta^2/2$'
         nplot=4
      end if

      if (iplot.eq.4) then
         if (iold.eq.-1) goto 1
         do i=1,n
            energy=u(i,1)**2+u(i,2)**2+u(i,3)**2
            v(i,1)=100.*u(i,1)**2/energy
            v(i,2)=100.*u(i,2)**2/energy
            v(i,3)=100.*u(i,3)**2/energy
            v(i,4)=100.*(energy-u(i,2)**2-u(i,7)**2)/energy
            v(i,5)=100.*u(i,7)**2/energy
         end do
         yaxis='% of energy$'
         heada='Contributions to the energy (%)$'
         headc='solid: u, dash: v, dot: w, chdh: Dv, chdt: eta$'
         nplot=5
      end if

      if (iplot.eq.5) then
         do i=1,n
            v(i,1)=u(i,1)
            v(i,2)=u(i,2)
            v(i,3)=u(i,3)
         end do
         yaxis='rms$'
         heada='Rms velocities$'
         headc='solid: u, dash: v, dot: w$'
         nplot=3
      end if

      if (iplot.eq.6) then
         do i=1,n
            v(i,1)=u(i,4)
            v(i,2)=u(i,5)
            v(i,3)=u(i,6)
         end do
         yaxis='rms$'
         heada='Rms vorticities$'
         headc='solid: xi, dash: eta, dot: zeta$'
         nplot=3
      end if

      if (iplot.eq.7) then
         do i=1,n
            v(i,1)=u(i,8)
         end do
         yaxis='dUuv$'
         heada='Rms of mean energy production$'
         headc=' $'
         nplot=1
      end if

      if (iplot.eq.8) then
         do i=1,n
            v(i,1)=u(i,10)
         end do
         yaxis='h+$'
         heada='Reynolds number based on friction velocity$'
         headc=' $'
         nplot=1
      end if

      xmin=t(1,1)
      xmax=t(n,1)
      ymin=0.
      ymax=0.
      if (iopt.eq.-1) then
         write(*,*) 'Give xmin, xmax, ymin, ymax'
         read(*,*) xmin,xmax,ymin,ymax
      end if
      xaxis='t$'
      ndim=10000
      idash=1
      mbox=1

      call rita1(t,v,xmin,xmax,ymin,ymax,npoint,nplot,ndim,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)
      goto 1

      end program pamp1
