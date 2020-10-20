c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rbla(fbla,dfbla,d2fbla,d3fbla,gbl,
     &     dgbl,d2gbl,th,dth,d2th,nbla,dybla,
     &     rlam,dstar2,mbla,namflo,fltype,m1,pr,scalar,thgrad)

      implicit none

      integer nbla,fltype,mbla,scalar
      real fbla(mbla),dfbla(mbla),d2fbla(mbla),d3fbla(mbla),dybla
      real gbl(mbla),dgbl(mbla),d2gbl(mbla)
      real th(mbla,scalar),dth(mbla,scalar),d2th(mbla,scalar)
      real dstar2
      real rlam,thgrad(scalar)
      character*80 namflo

      real m1(scalar),pr(scalar),m1in(scalar),prin(scalar)
      integer i,scalarin,j
      real a,b,yp,rlamin
c
c     Read header
c
      open(unit=17,file=namflo,status='old',form='unformatted')
      read(17)  a,b,nbla,rlamin,dstar2,scalarin

      if (scalarin.ge.1) then
         read(17) (prin(j),m1in(j),j=1,scalar)
      end if
c
c     Grid spacing
c
      dybla=(b-a)/real(nbla-1)

      write(*,*) 'Flow profile from file ',trim(namflo)
      write(*,*) '  eta (start) = ',a
      write(*,*) '  eta (end)   = ',b
      write(*,*) '  grid points = ',nbla
      write(*,*) '  m           = ',rlamin
      write(*,*) '  scalar      = ',scalarin
      do j=1,scalar
         write(*,*) '    Pr        = ',prin(j),j
         write(*,*) '    m1        = ',m1in(j),j
      end do
      write(*,*) '  d1          = ',dstar2,' (displacement thickness)'
c
c     Check file data
c
      if (nbla.gt.mbla) then
         write(*,*) 'The boundary layer table is too large'
         write(*,*) 'Increase mbla in bla.f to at least ',nbla
         write(*,*) 'and recompile'
         stop
      end if
      if (abs(rlam-rlamin).gt.1.e-14) then
         write(*,*) 'The flow field and the flow profile file'
         write(*,*) 'do not have the same pressure gradient'
         stop
      end if
      if (a.ne.0.) then
         write(*,*) 'Data should start at 0'
         stop
      end if
      if (scalar.gt.scalarin) then
         write(*,*) 'Scalar number not the same:'
         write(*,*) '  table: ',scalarin
         write(*,*) '  bls:   ',scalar
         stop
      end if
      if (scalar.lt.scalarin) then
         write(*,*) 'More scalars available in fsc.dat',scalarin
      end if
      do j=1,scalar
         pr(j) = prin(j)
         m1(j) = m1in(j)
      end do
c
c     Read an equidistant table
c
      if (scalar.eq.0) then
         do i=1,nbla
            read(17) yp,fbla(i),dfbla(i),d2fbla(i),d3fbla(i),
     &           gbl(i),dgbl(i),d2gbl(i)
         end do
      else
         do i=1,nbla
            read(17) yp,fbla(i),dfbla(i),d2fbla(i),d3fbla(i),
     &           gbl(i),dgbl(i),d2gbl(i),
     &           (th(i,j),dth(i,j),d2th(i,j),j=1,scalar)
         end do
c
c     Gradient of theta at the wall
c
         do j=1,scalar
            thgrad(j) = dth(1,j)*dstar2
         end do
      end if
c
c     Postprocess data
c
      if (fltype.eq.6.or.fltype.eq.3) then
c
c     Keep only f and its derivatives
c
         do i=1,nbla
            gbl(i)  = 0.0
            dgbl(i) = 0.0
         end do
      end if

      close(unit=17)

      end subroutine rbla
