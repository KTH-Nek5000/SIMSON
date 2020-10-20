c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program tchebe
c
c     Tests evaluation program chebe2
c
      integer nypmax,m
      parameter (nypmax=2000,m=10)
      real a(m,nypmax),eta(nypmax),wxy(m,nypmax),ai(m)
      integer y,nyp
      real prey(nypmax*2+15)
      real ym,pi
      parameter (pi = 3.1415926535897932385)
      integer j
c
 1000 continue
      write(*,*) 'Give number of points'
      read(*,*) nyp
      write(*,*) 'Give point to evaluate at'
      read(*,*) ym
c
c     Set up eta(y)
c
      do 2000 y=1,nyp
        eta(y)=cos(pi*real(y-1)/real(nyp-1))
 2000 continue
      do 2010 j=1,m
      do 2010 y=1,nyp
        a(j,y)=eta(y)**j
 2010 continue
      call vcosti(nyp,prey,0)
      call vchbf(a,wxy,nyp,m,m,1,prey)
      do 2020 y=1,nyp
      do 2020 j=1,m
        a(j,y)=a(j,y)*(2./real(nyp-1))
 2020 continue
      call chebe2(ai,eta(nyp),eta(1),a,m,nyp,ym)
      do 2030 j=1,m
      write(*,*) ai(j),ym**j
 2030 continue
      goto 1000

      end program tchebe
