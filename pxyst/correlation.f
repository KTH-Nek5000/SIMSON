c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine correlation(namxys)
c
c     Read correlation file
c
      implicit none

      include 'par.f'

      character*80 head1,head2,head3,namxys
      integer i,nx1,nyp1,nzc1,nfzsym1,fltype,corrnum
      real re,xl,zl,t,dum,dstar,bstart,blength,sumw
      logical pou

      real corrx(mcorr),corry(mcorr)
      real corrxx(mcorr),corryy(mcorr)
      integer corrxi(mcorr),corryi(mcorr)
      real corr_ms(4,mcorr)
      integer corrt(mcorr)
      real corr(nzc,mcorr),dz(nzc)
      integer z,it
      integer j,j1
      character(len=2) cc(17)
      data cc/'uu','vv','ww','pp','uv','uw','up','vw',
     &     'vp','wp','Du','tt','tu','tv','tw','tp','Dt'/

      write(*,*) 'Two-point correlation file'
      open(unit=19,file=namxys,form='unformatted')

      read(19) re,pou,xl,zl,t,dum
      re=-re
      read(19) nx1,nyp1,nzc1,nfzsym1
      read(19) fltype,dstar
      if (fltype.ge.6) read(19) bstart,blength
      read(19) sumw
      read(19) corrnum
c
c     Read the correlation
c
      do i=1,corrnum
c
c     Coordinate info
c
         read(19) corrx(i),corry(i),
     &        corrxi(i),corryi(i),
     &        corrxx(i),corryy(i),corrt(i)
c
c     mean and square
c
         read(19) (corr_ms(z,i),z=1,4)
c
c     And data
c
         read(19) (corr(z,i),z=1,nzc)

         do z=1,nzc
            corr(z,i) = (corr(z,i) - corr_ms(1,i)*corr_ms(2,i))/
     &           sqrt((corr_ms(3,i)-corr_ms(1,i)**2)*
     &           (corr_ms(4,i)-corr_ms(2,i)**2))
         end do

      end do
      close(19)



      write(*,*) '* read ',corrnum,' correlations.'
      write(*,*) 'Summation time ',sumw
      write(*,*) 'Time ',t
      write(*,*) nx,nyp,nzc
      write(*,*) xl,2./dstar,zl


 1000 write(*,'(a5,a10,a10)') '#','x','y'
      do i=1,corrnum
         j = abs(corrt(i))
         if (j.gt.11.and.j.le.11+6*scalar) then
            j=j-11
            j1 = mod(j-1,6)+1
            write(*,'(i5,f12.5,f10.5,a3,i3)') i,corrxx(i),corryy(i),
     &           cc(j1+11),(j-j1)/6+1
         else if (j.gt.11+6*scalar) then
            j=j-11
            j1 = mod(j-1,6)+1
            write(*,'(i5,f12.5,f10.5,a3,i3,i3)') i,corrxx(i),corryy(i),
     &           'TT',mod(j1-1,scalar)+1,(j1-1)/scalar+1
         else if (j.le.11) then
            write(*,'(i5,f12.5,f10.5,a3)') i,corrxx(i),corryy(i),
     &           cc(j)
         else
            write(*,'(i5,f12.5,f10.5,a3)') i,corrxx(i),corryy(i),
     &           '??'
         end if
      end do
      write(*,*) 'Give number of correlation (0 stop, <0 save)'
      read(*,*) it
      if (it.le.0.or.it.gt.corrnum) goto 1010
      do i=1,nzc
         dz(i) = (i-1)*zl/nzc
         write(*,*) i,dz(i),corr(i,it)
      end do

      write(head1,9510) t,re,sumw
 9510 format('Correlation',' t= ',F7.1,' re= ',F7.1,' sum= ',F9.4,'$')

      write(head2,9020) xl,zl,nx,ny,nz,namxys
 9020 format(' xl = ',F7.2,' zl= ',F6.2,' ',I4,'x',I3,'x',I3,
     &     ' file ',A20,'$')

      write(head3,9030) corrxx(it),corryy(it)
 9030 format(' x = ',F7.2,' y = ',F6.2,'$')

      call rita1(dz,corr(1,it),0.,0.,0.,0.,nzc,1,1,
     &     'deltaz$','corr$',head1,head2,head3,
     &     1,1)
      goto 1000

 1010 continue
      if (it.lt.0) then
c
c     Save correlation
c
         write(*,*) 'Saving correlation data'
         open(file='correlation.dat',unit=10)
c         write(10,*) 'dz',(corrxx(it),it=1,corrnum)
         do i=1,nzc
            dz(i) = (i-1)*zl/nzc
            write(10,'(1000e18.9)') dz(i),(corr(i,it),it=1,corrnum)
         end do
         close(10)
      end if

      stop

      end subroutine correlation
