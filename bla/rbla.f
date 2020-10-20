c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rbla(fbla,nbla,dybla,
     &     rlam,dstar2,fltype,m1,pr,my_node,x0,re,dstar,thgrad)
c
c     Reads the base flow table 'fsc.dat', which needs to be created by
c     fsc.f (in bls-heat directory). At a later time, fsc.f could be
c     included into bla as well.
c
c     The contents of the fbla array are as follows:
c     fbla(1:nbla, 1)    f
c     fbla(1:nbla, 2)    f'
c     fbla(1:nbla, 3)    f''
c     fbla(1:nbla, 4)    f'''
c     fbla(1:nbla, 5)    g      (only if crossflow active)
c     fbla(1:nbla, 6)    g'     (only if crossflow active)
c     fbla(1:nbla, 7)    g''    (only if crossflow active)
c     fbla(1:nbla, 8)    th     (only if scalar=1)
c     fbla(1:nbla, 9)    th'    (only if scalar=1)
c     fbla(1:nbla,10)    th''   (only if scalar=1)
c
      implicit none

      include 'par.f'

      integer nbla,fltype
      real fbla(mbla,7+3*scalar)
      real dstar2,dybla
      real rlam

      real m1(scalar),pr(scalar),prin(scalar),m1in(scalar)
      integer i,scalarin,j
      real a,b,yp,rlamin,x0,re,dstar,thgrad(scalar)

      integer my_node
      logical recreate

      recreate = .false.
c
c     Read header
c
      open(unit=17,file='fsc.dat',
     &     status='old',form='unformatted',err=1111)
      read(17,err=1111,end=1111) a,b,nbla,rlamin,dstar2,scalarin

      if (scalarin.ge.1) then
         read(17) (prin(j),m1in(j),j=1,scalar)
      end if
c
c     Grid spacing
c
      dybla=(b-a)/real(nbla-1)

      if (my_node.eq.0) then
         write(ios,*) ''
         write(ios,*) 'Reading base-flow profile '
         write(ios,*) '  eta (start) : ',a
         write(ios,*) '  eta (end)   : ',b
         write(ios,*) '  grid points : ',nbla
         write(ios,*) '  m           : ',rlamin
         write(ios,*) '  scalar      : ',scalarin
         do j=1,scalar
            write(ios,*) '    Pr        : ',prin(j),j
            write(ios,*) '    m1        : ',m1in(j),j
         end do
         write(ios,*) '  dstar2      : ',dstar2,
     &        ' (displacement thickness)'
      end if
c
c     Check file data
c
      if (nbla.gt.mbla) then
         if (my_node.eq.0) then
            write(ioe,*) 'The base-flow table is too large.'
            write(ioe,*) 'Increase mbla in par.f to at least ',nbla,
     &           ' (now: ',mbla,')'
            write(*,*) 'and recompile'
         end if
         recreate = .true.
      end if
      if (rlam.ne.rlamin) then
         if (my_node.eq.0) then
            write(*,*) 'The flow field and the flow-profile file'
            write(*,*) 'do not have the same pressure gradient'
         end if
         recreate = .true.
      end if
      if (scalar.eq.1.and.scalarin.eq.0) then
         if (my_node.eq.0) then
            write(ioe,*) 'Scalar flag not the same:'
            write(ioe,*) '  table: ',scalarin
            write(ioe,*) '  bla:   ',scalar
         end if
         recreate = .true.
      end if
      if (scalar.gt.scalarin) then
         if (my_node.eq.0) then
            write(ioe,*) 'Scalar number not the same:'
            write(ioe,*) '  table: ',scalarin
            write(ioe,*) '  bls:   ',scalar
         end if
         recreate = .true.
      end if
      if (scalar.lt.scalarin) then
         if (my_node.eq.0) then
            write(ioe,*) 'More scalars available in fsc.dat',scalarin
         end if
      end if
      do j=1,scalar
         if (pr(j).ne.prin(j)) then
            if (my_node.eq.0) then
               write(ioe,*) 'The flow field and the flow-profile file'
               write(ioe,*) 'do not have the same Prandtl number'
               write(ioe,*) 'Scalar nr: ',j
            end if
            recreate = .true.
         end if
         if (m1(j).ne.m1in(j)) then
            if (my_node.eq.0) then
               write(ioe,*) 'The flow field and the flow-profile file'
               write(ioe,*) 'do not have the same exponent ' //
     &              'for the scalar'
               write(ioe,*) 'Scalar nr: ',j
            end if
            recreate = .true.
         end if
      end do
      if (a.ne.0.) then
         write(ioe,*) 'Data should start at eta=0'
         recreate = .true.
      end if

      goto 1112

 1111 continue
      recreate = .true.
      close(17)

 1112 continue

      if (recreate) then
c
c     Compute similarity solution
c
         if (my_node.eq.0) then
            write(ioe,*)
     &           '** Please provide suitable base flow in fsc.dat'
         end if
         call stopnow(123554)

      else
c
c     Read an equidistant table
c
         do i=1,nbla
            read(17) yp,(fbla(i,j),j=1,7+3*scalar)
         end do

         close(unit=17)

         do j=1,scalar
            thgrad(j) = fbla(1,9+3*(j-1))*dstar2
         end do

      end if
c
c     Postprocess data
c
c     For 2D boundary layers, keep only f and its derivatives
c
      if (fltype.eq.6.or.fltype.eq.3) then
         do i=1,nbla
            fbla(i,6)  = 0.0
            fbla(i,7)  = 0.0
         end do
      end if
c
c     x0=Re/(2*dstar2*2) is the inflow location of a flat plate from the
c     leading edge x=0 in inflow displacement thickness. Here x0 is scaled
c     by dstar=2/h2 (h2 is y-box size in initial displacement thickness,
c     dstar2=1.2168 in the Blasius boundary layer).
c     This is needed for e.g. the correct temporal growth of the
c     boundary layer in the temporal case.
c
      x0   = dstar**2*re*(rlam+1.)/(dstar2**2*2.)
      if (my_node.eq.0) then
         write(ios,*) '  x0          : ',x0/dstar
      end if

      end subroutine rbla
