c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rwavebl(xl,zl,evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,
     &     ampob,amp2d,prey,dstar,ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &     osmod,osn,osnumb,osymax,osbeta,osre,osomega,osalr,osali,
     &     osur,osui,osvr,osvi,oswr,oswi,osfil,my_node,osdamp)
c
c     Reads in wave disturbances (e.g. OS modes)
c
      implicit none

      include 'par.f'

      real evr(nyp),evi(nyp),eur(nyp),eui(nyp),ewr(nyp),ewi(nyp)
      real eer(nyp),eei(nyp)
      real ev2r(nyp),ev2i(nyp),eu2r(nyp),eu2i(nyp),afw2d,alf2d
      real prey(nyp*2+15)
      real alfaw,afw,dstar,betaw,ampob,amp2d

      real rew,k2w
      complex eigv
      complex eigmo(500),alfac
      real eigmr(500),eigmi(500),eigmsr(500),eigmsi(500)
      real edvr(nyp),edvi(nyp),w(nyp),vmax
      real c1,c2,c3,c4,k2r,k2i,k22
      logical fuleig,osmod
      integer i,y,n,osnumb,vers,numb,modus,osmodus

      real modd

      integer osn, osvers
      logical osdamp
      real yyy(1:400)
      real osymax,osbeta(osnf), osre, osomega(osnf)
      real osalr(osnf), osali(osnf)
      real oser(osnf,nyp), osur(osnf,nyp), osvr(osnf,nyp)
      real oswr(osnf,nyp)
      real osei(osnf,nyp), osui(osnf,nyp), osvi(osnf,nyp)
      real oswi(osnf,nyp)
      real osscale(osnf),re,ymax
      complex osgamma(osnf),eof,eofr
      integer jj
      character*80 osfil

      character*20 time

      real xl,zl,kxmax,kzmax,kymax,kzf,kymin
      integer j
      integer my_node

      real pi
      parameter (pi = 3.1415926535897932385)

      eof=cmplx(12.34,56.78)
      fuleig=.false.


      if (ampob.ne.0.) then
c
c     Reads eigenmode for oblique waves from the file wave.dat
c
         open(unit=10,status='old',file='wave.dat')
         read(10,*) rew
         if (ampob.gt.0.) then
            read(10,*) alfac,betaw
         else
            read(10,*) alfaw,betaw
            alfac=cmplx(alfaw,0.)
         end if
         read(10,*) eigv
         alfac=alfac/dstar
         alfaw=real(alfac)
         betaw=betaw/dstar

         k2w=alfaw**2+betaw**2
         if (ampob.gt.0.) then
            do i=1,500
               read(10,*,end=1000) eigmr(i),eigmi(i),eigmsr(i),eigmsi(i)
            end do
         else
            do i=1,500
               read(10,*,end=1000) eigmo(i)
            end do
         end if
         i=1001
 1000    continue
         n=i-1
         close(unit=10)

         afw=real(eigv)/dstar
c
c     Turn plane upside down
c
         if (ampob.lt.0.) then
            do y=1,min(nyp,n)
               evr(y)=-(-1)**y*real(eigmo(y))
               evi(y)=-(-1)**y*aimag(eigmo(y))
               eer(y)=0.
               eei(y)=0.
            end do
            ampob=abs(ampob)
         else
            do y=1,min(nyp,n)
               evr(y)=-(-1)**y*eigmr(y)
               evi(y)=-(-1)**y*eigmi(y)
               eer(y)=-(-1)**y*eigmsr(y)
               eei(y)=-(-1)**y*eigmsi(y)
            end do
         end if
c
c     Get u=i*(alfa*dv-beta*eta)/k2 and w=i*(alfa*eta+beta*dv)
c
         call dcheb(edvr,evr,nyp,1,1)
         call dcheb(edvi,evi,nyp,1,1)

         if (fuleig) then
c
c     Get the u and w using full complex eigenvalue, alfa.
c
            k2r=alfaw*alfaw-aimag(alfac)*aimag(alfac)+betaw*betaw
            k2i=2.*alfaw*aimag(alfac)
            k22=k2r*k2r+k2i*k2i
            c1=(alfaw*k2i-aimag(alfac)*k2r)/k22
            c2=(alfaw*k2r+aimag(alfac)*k2i)/k22
            c3=betaw*k2i/k22
            c4=betaw*k2r/k22
            do y=1,nyp
               eur(y)=c1*edvr(y)-c2*edvi(y)-(c3*eer(y)+c4*eei(y))/dstar
               eui(y)=c1*edvi(y)+c2*edvr(y)-(c3*eei(y)-c4*eer(y))/dstar
               ewr(y)=(c1*eer(y)-c2*eei(y))/dstar+c3*edvr(y)-c4*edvi(y)
               ewi(y)=(c1*eei(y)+c2*eer(y))/dstar+c3*edvi(y)+c4*edvr(y)
            end do
         else
            do y=1,nyp
               eur(y)=(-alfaw*edvi(y)+betaw*eei(y)/dstar)/k2w
               eui(y)=( alfaw*edvr(y)-betaw*eer(y)/dstar)/k2w
               ewr(y)=(-betaw*edvi(y)-alfaw*eei(y)/dstar)/k2w
               ewi(y)=( betaw*edvr(y)+alfaw*eer(y)/dstar)/k2w
            end do
         end if
c
c     Transform to physical space
c
         call vchbb(evr,w,nyp,1,1,1,prey)
         call vchbb(evi,w,nyp,1,1,1,prey)
         call vchbb(eur,w,nyp,1,1,1,prey)
         call vchbb(eui,w,nyp,1,1,1,prey)
         call vchbb(ewr,w,nyp,1,1,1,prey)
         call vchbb(ewi,w,nyp,1,1,1,prey)
c
c     Normalize the wave to vmax=1
c
         vmax=0.
         do y=1,nyp
            vmax=max(vmax,sqrt(evr(y)**2+evi(y)**2))
         end do
         do y=1,nyp
            evr(y)=evr(y)/vmax*ampob
            evi(y)=evi(y)/vmax*ampob
            eur(y)=eur(y)/vmax*ampob
            eui(y)=eui(y)/vmax*ampob
            ewr(y)=ewr(y)/vmax*ampob
            ewi(y)=ewi(y)/vmax*ampob
c     eta(y)=cos(pi*real(y-1)/real(nyp-1))
c       write(43,*) (eta(y)+1.)/dstar,sqrt(eur(y)**2+eui(y)**2)
c    &  ,sqrt(evr(y)**2+evi(y)**2),sqrt(ewr(y)**2+ewi(y)**2)
c    &        evr(y),evi(y),eur(y),eui(y),ewr(y),ewi(y)
         end do
      else
         do y=1,nyp
            evr(y)=0.
            evi(y)=0.
            eur(y)=0.
            eui(y)=0.
            ewr(y)=0.
            ewi(y)=0.
         end do
      end if

      if (amp2d.ne.0.) then
c
c     Reads eigenmode for 2d wave from the file wave2.d
c
         open(unit=10,status='old',file='wave2.d')
         read(10,*) rew
         read(10,*) alf2d
         read(10,*) eigv
         alfaw=alfaw/dstar
         k2w=alfaw**2

         do i=1,1000
            read(10,*,end=5000) eigmr(i),eigmi(i)
         end do
         i=1001
 5000    continue
         n=i-1
         close(unit=10)

         afw2d=real(eigv)/dstar

c     Turn plane upside down
         do y=1,min(nyp,n)
            ev2r(y)=-(-1)**y*eigmr(y)
            ev2i(y)=-(-1)**y*eigmi(y)
         end do
c
c     Get u
c
         call dcheb(edvr,ev2r,nyp,1,1)
         call dcheb(edvi,ev2i,nyp,1,1)

         if (fuleig) then
c
c     Get the u and w using full complex eigenvalue, alfa.
c
            k2r=alfaw*alfaw-aimag(alfac)*aimag(alfac)+betaw*betaw
            k2i=2.*alfaw*aimag(alfac)
            k22=k2r*k2r+k2i*k2i
            c1=(alfaw*k2i-aimag(alfac)*k2r)/k22
            c2=(alfaw*k2r+aimag(alfac)*k2i)/k22
            do y=1,nyp
               eu2r(y)=c1*edvr(y)-c2*edvi(y)
               eu2i(y)=c1*edvi(y)+c2*edvr(y)
            end do
         else
            do y=1,nyp
               eu2r(y)=-alfaw*edvi(y)/k2w
               eu2i(y)=alfaw*edvr(y)/k2w
            end do
         end if
c
c     Transform to physical space
c
         call vchbb(ev2r,w,nyp,1,1,1,prey)
         call vchbb(ev2i,w,nyp,1,1,1,prey)
         call vchbb(eu2r,w,nyp,1,1,1,prey)
         call vchbb(eu2i,w,nyp,1,1,1,prey)
c
c     Normalize the wave to vmax=1
c
         vmax=0.
         do y=1,nyp
            vmax=max(vmax,sqrt(evr(y)**2+evi(y)**2))
         end do
         do y=1,nyp
            ev2r(y)=ev2r(y)/vmax*amp2d
            ev2i(y)=ev2i(y)/vmax*amp2d
            eu2r(y)=eu2r(y)/vmax*amp2d
            eu2i(y)=eu2i(y)/vmax*amp2d
         end do
      else
         do y=1,nyp
            ev2r(y)=0.
            ev2i(y)=0.
            eu2r(y)=0.
            eu2i(y)=0.
         end do
      end if

C ********************************************************************

      if (osmod) then
c
c     Read the Orr-Sommerfeld/Squire modes for free-stream turbulence
c
         if (my_node.eq.0) then
            write(ios,*) '* Reading file ',trim(osfil)
         end if
         open(unit=10,status='old',file=osfil,
     &        form='unformatted')

c
c     There might be a time-string in the osmodes file
c
         time = 'not saved!'
         read(10,err=7777,end=7777) time
         goto 7778
 7777    continue
         backspace(10)
 7778    continue

         read(10) vers
         if (my_node.eq.0) then
            write(ios,*) '* Version of wave file: ',-vers
         end if
         if (vers.gt.-4) then
            write(ioe,*) 'Impossible to process due to version ',-vers
         end if
         read(10) n, re, ymax
         read(10) osomega(1),osalr(1),osali(1),osbeta(1),osgamma(1)
         read(10) numb,modus

         if (n+1.ne.nyp) then
            write(ioe,*) 'Wrong dimension of ',trim(osfil)
            write(ioe,*) n+1,nyp
            stop
         end if

         if (numb.gt.osnf) then
            write(ioe,*) 'Increase osnf in par.f'
            write(ioe,*) 'Number of eigenfunctions ',numb
            write(ioe,*) 'In bla ',osnf
            stop
         end if
c
c     Maximum wavenumber for the present grid
c
         kxmax = pi/(xl/dstar/nx)
         kzmax = pi/(zl/dstar/nz)
         kymax = pi/
     &        (-(1.+cos(pi*real( (ny-1)/2)/real(nyp-1)))/dstar+
     &        (1.+cos(pi*real((ny-1)/2-1)/real(nyp-1)))/dstar)

         kymin = 2*pi/(2./dstar)
         kymin = 0.
         kzf = pi/(zl/dstar)
         if (my_node.eq.0) then
            write(ios,*) '* Limiting wavenumbers'
            write(ios,*) '  kxmax = ',kxmax
            write(ios,*) '  kymax = ',kymax
            write(ios,*) '  kymin = ',kymin
            write(ios,*) '  kzmax = ',kzmax
         end if
c
c     Start reading the osmodes file again with all the modes
c     in a loop
c
         rewind(10)

         time = 'not saved!'
         read(10,err=8888,end=8888) time
         goto 8889
 8888    continue
         backspace(10)
 8889    continue

         jj = 1
         do j=1, numb
            read(10) osvers
            read(10) osn, osre,  osymax
            read(10) osomega(jj),osalr(jj),osali(jj),
     &           osbeta(jj),osgamma(jj)
            read(10) osnumb,osmodus
            read(10) osscale(jj)
            read(10) (yyy(i),oser(jj,i),osei(jj,i),osur(jj,i),
     &           osui(jj,i),osvr(jj,i),osvi(jj,i),oswr(jj,i),
     &           oswi(jj,i),i=1,osn+1)
            read(10) eofr

            if (eof.ne.eofr) then
               write(ioe,*) 'eof wrong.'
               stop
            end if

            if ((yyy(1).ne.osymax).or.(yyy(osn+1).ne.0.)) then
               write(ioe,*) 'Wall-normal direction wrong.'
               stop
            end if
            call checki(n, osn)
            call checkr(re, osre)
            call checkr(ymax, osymax)
            call checki(numb, osnumb)
            call checki(modus, osmodus)
            call checki(vers, osvers)

c            if (abs(modd(osbeta(jj),kzf)).gt.4.*epsilon(kzf)) then
            if (abs(modd(osbeta(jj),kzf)).gt.1e-6) then
               write(ioe,*) modd(osbeta(jj),kzf),'>',4.*epsilon(kzf)
               write(ioe,*) 'osbeta(jj) = ',osbeta(jj),'(',jj,')'
               write(ioe,*)
     &              'not multiple of fundamental wavenumber ',kzf
            end if
c
c     If the OS/SQ modes contain modes with large decay rates
c     (i.e. large osali) it is inpractical to actually take that
c     large growth in the fringe region (i.e. at positions upstream of
c     the inflow) into account. Then, osali is set to 0:
c     (osdamp = .false.)
c
c     damp1 = osamp
c     damp2 = osamp
c
c     However, if an accurate representation of the modes is necessary,
c     then the full mode with growth/decay should be superimposed:
c     (osdamp = .true.)
c
c     damp1 = osamp*exp(-osali(j)*xc1(x))
c     damp2 = osamp*exp(-osali(j)*xc2(x))
c
c
c     Check the wavenumbers
c
            if (osalr(jj).gt.kxmax .or.
     &           osbeta(jj).gt.kzmax .or.
     &           real(osgamma(jj)).gt.kymax.or.
     &           real(osgamma(jj)).lt.kymin) then
            else
c
c     If wavenumber within range, do the boundary layer scaling and
c     increment jj
c
               osalr(jj)   = osalr(jj)   / dstar
               osali(jj)   = osali(jj)   / dstar
               osbeta(jj)  = osbeta(jj)  / dstar
               osomega(jj) = osomega(jj) / dstar

               jj = jj + 1
            end if
         end do
c
c     Correct number of points and number of modes
c
         osn    = osn+1
         osnumb = jj-1

         close(10)
         if (my_node.eq.0) then
            write(ios,*) '* Total    ',numb,' eigenfunctions.'
            write(ios,*) '* Retained ',osnumb,' eigenfunctions.'
            write(ios,*) '* OSDAMP   ',osdamp
            write(ios,*) '* Generated: ',time
            write(ios,*) '* Eigenmodes finished.'
         end if
c     for bla we need osalr, osali, osbeta, osomega, osur....
      end if

      end subroutine rwavebl


      subroutine checkr(a1,a2)
      implicit none
      real a1,a2
      if (a1.ne.a2) then
         write(*,*) 'Parameters do not agree.'
         write(*,*) a1,a2
         call stopnow(54353)
      end if
      end subroutine checkr


      subroutine checki(a1,a2)
      implicit none
      integer a1,a2
      if (a1.ne.a2) then
         write(*,*) 'Parameters do not agree.'
         write(*,*) a1,a2
         call stopnow(42423)
      end if
      end subroutine checki


      real function modd(a,d)
      implicit none
      real a,d

      modd = a-d*nint(a/d)

      end function modd
