c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdiscp(ur,re,xl,zl,t,xs,dstar,fltype,
     &bstart,blength,rlam,spanv,spat,namnin,m,pln,urx)
c
c     Reads m variables from file namnin and puts into ur
c
      implicit none

      include 'par.f'

      character*80 namnin
      integer m,fltype
      complex ur(memnx,memny,memnz,7)
      real re,xl,zl,t,xs,dstar,rlam,spanv,bstart,blength
      logical pou,spat
      complex urx(nx/2)
      complex pln(nx/2,nyp)

      integer x,y,z,i,ic
      integer nxin,nypin,nzcin,nfzsin
      open(unit=12,file=namnin,status='old',form='unformatted')
c
c     Read file header
c
      read(12) re,pou,xl,zl,t,xs
      read(12) nxin,nypin,nzcin,nfzsin
      read(12,*) fltype,dstar
      if (fltype.ne.3.and.(fltype.lt.0.or.fltype.lt.6.or.fltype.gt.9))
     &     then
         write(*,*) 'The input file does not contain'
         write(*,*) 'the correct type of flow'
         stop
      end if
      if (spat.and.fltype.lt.6.or.(.not.spat.and.fltype.ge.6)) then
         write(*,*) 'spat and fltype does not match '
         write(*,*) 'Change spat in bla.i or use other flow field'
         stop
      end if
      if (fltype.eq.-1) read(12) rlam
      if (fltype.eq.-2) read(12) rlam,spanv
      if (fltype.eq.6) then
         read(12) bstart,blength
         rlam=0.0
         spanv=0.0
      end if
      if (fltype.ge.7) read(12) bstart,blength,rlam,spanv
c
c     Check file info
c
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &    nfzsin.ne.nfzsym) then
         write(*,*) 'Input file has a size other than program'
         write(*,*) 'File parameters nx,nyp,nzc,nfzsym',
     &        nxin,nypin,nzcin,nfzsin
         write(*,*) 'Program parameters',nx,nyp,nzc,nfzsym
         stop
      end if
      write(*,*) 'Reading data'
      ic=0
      do i=1,m
         do z=1,nzc
            if (ic.eq.nzc*m/100+1) write(*,*) 'Read 1%'
            if (ic.eq.nzc*m/10+1) write(*,*) 'Read 10%'
            if (ic.eq.nzc*m/2+1) write(*,*) 'Read 50%'
            if (ic.eq.nzc*m*9/10+1) write(*,*) 'Read 90%'
            ic=ic+1
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
            call putxyp(pln,z,i,ur)
         end do
      end do
      close(unit=12)
      return
 1000 write(*,*) 'A premature EOF was found'
      write(*,*) 'Read ',3+((i-1)*nzc+z-1)*nyp+y-1,' records'
      write(*,*) 'Expected ',3+m*nzc*nyp
      write(*,*) 'Output may be incorrect'
      close(unit=12)

      return

      end subroutine rdiscp
