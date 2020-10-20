c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdiscp(ur,re,xl,zl,t,xs,dstar,fltype,
     &     namnin,m,pln,urx)
c
c     Reads m variables from file namnin and puts into ur
c
      implicit none

      include 'par.f'

      character*80 namnin
      integer m,fltype
      complex ur(memnx,memny,memnz,7)
      real re,xl,zl,t,xs,dstar
      logical pou
      complex urx(nx/2)
      complex pln(nx/2,nyp)

      integer x,y,z,i,j,scalarind
      integer nxin,nypin,nzcin,nfzsin
      real bstart,blength,rlam,spanv
      real m1(scalar),pr(scalar),gr(scalar)
c
c     Read file header
c
      do j=1,scalar
         pr(j) = 0.
         m1(j) = 0.
         gr(j) = 0.
      end do

      open(unit=12,file=namnin,status='old',form='unformatted')
      read(12,err=1009) re,pou,xl,zl,t,xs,(pr(j),m1(j),j=1,scalar)

      goto 1008

 1009 continue
      write(*,*) 'Error in header. Check scalar variable in par.f.'
      write(*,*) 'Also check big/little endian compiler option.'
      stop

 1008 continue


      read(12) nxin,nypin,nzcin,nfzsin
      fltype = 0
      dstar  = 0.0
      read(12,err=1010) fltype,dstar

      if (fltype.eq.-1) read(12) rlam
      if (fltype.eq.-2) read(12) rlam,spanv

      if (fltype.eq.4.or.fltype.eq.5) read(12) bstart,blength
      if (fltype.eq.6) read(12) bstart,blength
      if (fltype.ge.7.and.fltype.le.9) 
     &     read(12) bstart,blength,rlam,spanv
      if (abs(fltype).eq.20) read(12) (gr(j),j=1,scalar)

 1010 if (fltype.ge.-3.and.fltype.le.-1.or.
     &     fltype.ge.1.and.fltype.le.11.or.
     &     abs(fltype).eq.20) then
c
c     valid flow type
c     
      else
         write(*,*) 'Invalid flow type ',fltype,dstar
      end if


      scalarind = 0
      if (scalar.gt.0) then
         write(*,'(a5,3a20)') 'No','Pr  ','m1  ','Gr  '
         do j=1,scalar
            write(*,'(i5,3e20.9)') j,pr(j),m1(j),gr(j)
         end do
         write(*,*) 'Which scalar to read? max=',scalar,
     &        '(0 no scalar)'
         read(*,*) scalarind
      end if




c
c     Rescale coordinates in the boundary layer case
c
      if (fltype.eq.1) dstar=1.
      if (fltype.eq.2) dstar=1.
      if (fltype.eq.4) dstar=1.
      if (fltype.eq.5) dstar=1.
      re=re*dstar
      xl=xl/dstar
      zl=zl/dstar
      t=t/dstar
      xs=xs/dstar
c
c     Check file info
c
      write(*,*)
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &    nfzsin.ne.nfzsym) then
         write(*,*) 'Input file has a size other than program'
         write(*,*) 'File parameters nx,nyp,nzc,nfzsym',
     &        nxin,nypin,nzcin,nfzsin
         write(*,*) 'Program parameters',nx,nyp,nzc,nfzsym
         stop
      end if
      write(*,*) 'Reading data'
c
c     Read in velocity and scalar
c
      do i=1,m+scalar
         do z=1,nzc
            if (z+(i-1)*nzc.eq.nzc*(m+scalar)/100+1)
     &           write(*,*) 'Read 1%'
            if (z+(i-1)*nzc.eq.nzc*(m+scalar)/10+1)
     &           write(*,*) 'Read 10%'
            if (z+(i-1)*nzc.eq.nzc*(m+scalar)/2+1)
     &           write(*,*) 'Read 50%'
            if (z+(i-1)*nzc.eq.nzc*(m+scalar)*9/10+1)
     &           write(*,*) 'Read 90%'
            do y=1,nyp
               read(12,end=1000) urx
               do x=1,nx/2
                  pln(x,y)=urx(x)
               end do
            end do
c
c     Zero odd-ball
c
            if (nfzsym.eq.0.and.z.eq.nzc/2+1.and.mod(nz,2).ne.1) then
               do y=1,nyp
                  do x=1,nx/2
                     pln(x,y)=0.0
                  end do
               end do
            end if
c
c     For now, store the scalar field in the w-component
c     of the velocity field
c
            if (scalarind.gt.0.and.i-3.eq.scalarind) then
               call putxyp(pln,z,3,ur)
            else if (i.le.3) then
               call putxyp(pln,z,i,ur)
            end if

         end do
      end do
      close(unit=12)

      if (scalar.gt.0.and.scalarind.gt.0) then
         write(*,*) 'Scalar = ',scalar,' put into w'
         write(*,*) 'pr     = ',pr(scalarind)
         write(*,*) 'm1     = ',m1(scalarind)
      end if
      return

 1000 write(*,*) 'A premature EOF was found'
      write(*,*) 'Read ',3+((i-1)*nzc+z-1)*nyp+y-1,' records'
      write(*,*) 'Expected ',3+m*nzc*nyp
      write(*,*) 'Output may be incorrect'
      close(unit=12)
c
c     Program could be continued, however it's probably more
c     save to exit
c
c      stop

      end subroutine rdiscp
