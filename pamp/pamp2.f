c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pamp2
c
c     Plots a  number of components from a number of amplitude files.
c     Compares data from a number of files, both old and new formats
c
      implicit none

      integer i,j,k,nplot,iplot,iopt,idash,ntot,mbox,inorm,nmax,ndim
      real xmin,xmax,ymin,ymax,fnorm
      parameter(nmax=10000)
      character*80 xaxis,yaxis,heada,headb,headc
      character*32 namamp(5)
      character*6 dash(5)
      integer npoint(5),nn(5),iold(5),y,nyp
      real t(nmax,5),v(nmax,5),u(nmax,10,5),dum
      data dash/' solid',', dash',',  dot',', chdh',', chdt'/

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         pamp2 $Rev$ '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
      write(*,*) 'Give number of amplitude files (max 5)'
      read(*,*) nplot
      write(*,*) 'Give the names of the amplitude file, end with ,CR>'
      do k=1,nplot
         read(*,222) namamp(k)
 222     format(a)
      end do
      write(*,*)
c
c     Read file namamp
c
      do k=1,nplot
         nn(k)=index(namamp(k),' ')
         write(*,222) namamp(k)(1:nn(k)-1),
     &' 1=short(11), 2=long(11+nyp)'
c     &' old or new? (-1=real old, 0=old, 1=new,  2=new, long)'
         read(*,*) iold(k)
         if (iold(k).eq.2) then
           write(*,*) 'Give nyp'
           read(*,*) nyp
         end if
         open(unit=10,file=namamp(k),status='old')
         do i=1,nmax
            read(10,*,end=4) t(i,k),(u(i,j,k),j=1,10)
c            if (iold(k).eq.-1) read(10,*,end=4) (u(i,j,k),j=4,6)
c            if (iold(k).ge.0)  read(10,*,end=4) (u(i,j,k),j=4,7)
c            if (iold(k).ge.1)  read(10,*,end=4) (u(i,j,k),j=8,10)
            if (iold(k).eq.2) then
              do j=1,10
                read(10,*,end=4) (dum,y=1,nyp)
             end do
            end if
         end do
 4       npoint(k)=i-1
         close(unit=10)
      end do
      headb='from files '//namamp(1)(1:nn(1)-1)
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

 1    write(*,*) ' Negative value for more options'
      write(*,*) ' 1 = Energy'
      write(*,*) ' 2 = u rms'
      write(*,*) ' 3 = v rms'
      write(*,*) ' 4 = w rms'
      write(*,*) ' 5 = xi rms'
      write(*,*) ' 6 = eta rms'
      write(*,*) ' 7 = zeta rms'
      write(*,*) ' 8 = Dv/k rms'
      write(*,*) ' 9 = eta/k rms'
      if (iold(1).ge.1) then
         write(*,*) ' 10 = dUuv (averaged Reynolds stress)'
         write(*,*) ' 11 = h+ (Reynolds number based on u_tau'
      end if
      read(*,*) iplot
      iopt=1
      if (iplot.le.0) iopt=-1
      iplot=abs(iplot)

      if (iplot.eq.1) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1)
     &           fnorm=0.5*(u(1,1,k)**2+u(1,2,k)**2+u(1,3,k)**2)
            do i=1,npoint(k)
               v(i,k)=0.5*(u(i,1,k)**2+u(i,2,k)**2+u(i,3,k)**2)/fnorm
            end do
         end do
         yaxis='energy$'
         heada='Energy$'
      end if

      if (iplot.eq.2) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,1,k)
            do i=1,npoint(k)
               v(i,k)=u(i,1,k)/fnorm
            end do
         end do
         yaxis='u$'
         heada='Rms of u$'
      end if

      if (iplot.eq.3) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,2,k)
            do i=1,npoint(k)
               v(i,k)=u(i,2,k)/fnorm
            end do
         end do
         yaxis='v$'
         heada='Rms of v$'
      end if

      if (iplot.eq.4) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,3,k)
            do i=1,npoint(k)
               v(i,k)=u(i,3,k)/fnorm
            end do
         end do
         yaxis='w$'
         heada='Rms of w$'
      end if

      if (iplot.eq.5) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,4,k)
            do i=1,npoint(k)
               v(i,k)=u(i,4,k)/fnorm
            end do
         end do
         yaxis='xi$'
         heada='Rms of streamwise vorticity$'
      end if

      if (iplot.eq.6) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,5,k)
            do i=1,npoint(k)
               v(i,k)=u(i,5,k)/fnorm
            end do
         end do
         yaxis='eta$'
         heada='Rms of normal vorticity$'
      end if

      if (iplot.eq.7) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,6,k)
            do i=1,npoint(k)
               v(i,k)=u(i,6,k)/fnorm
            end do
         end do
         yaxis='zeta$'
         heada='Rms of spanwise vorticity$'
      end if

      if (iplot.eq.8) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1)
     &           fnorm=sqrt(u(1,1,k)**2+u(1,3,k)**2-u(1,7,k)**2)
            do i=1,npoint(k)
               v(i,k)=sqrt(u(i,1,k)**2+u(i,3,k)**2-u(i,7,k)**2)/fnorm
            end do
         end do
         yaxis='Dv/k$'
         heada='Rms of normal velocity derivative over wavenumber$'
      end if

      if (iplot.eq.9) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,7,k)
            do i=1,npoint(k)
               v(i,k)=u(i,7,k)/fnorm
            end do
         end do
         yaxis='eta/k$'
         heada='Rms of normal vorticity over wavenumber$'
      end if

      if (iplot.eq.10) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,8,k)
            do i=1,npoint(k)
               v(i,k)=u(i,8,k)/fnorm
            end do
         end do
         yaxis='dUuv$'
         heada='Rms of mean energy production$'
      end if

      if (iplot.eq.11) then
         do k=1,nplot
            fnorm=1.
            if (inorm.eq.1) fnorm=u(1,10,k)
            do i=1,npoint(k)
               v(i,k)=u(i,10,k)/fnorm
            end do
         end do
         yaxis='h+$'
         heada='Reynolds number based on friction velocity$'
      end if

      xmin=10000.
      xmax=0.
      do k=1,nplot
         xmin=min(xmin,t(1,k))
         xmax=max(xmax,t(npoint(k),k))
      end do
      ymin=0.
      ymax=0.
      if (iopt.eq.-1) then
         write(*,*) 'Give xmin, xmax, ymin, ymax'
         read(*,*) xmin, xmax, ymin, ymax
      end if
      xaxis='t$'
      ndim=nmax
      idash=1
      mbox=1

      call rita1(t,v,xmin,xmax,ymin,ymax,npoint,nplot,ndim,
     &           xaxis,yaxis,heada,headb,headc,mbox,idash)
      goto 1

      end program pamp2
