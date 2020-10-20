c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdiscr(ur,re,xl,zl,t,xs,dstar,fltype,
     &namnin,m,pln,urx,header)
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

      logical header
      integer x,y,z,i
      integer nxin,nypin,nzcin,nfzsin
c
c     Read file header
c
      open(unit=12,file=namnin,status='old',form='unformatted')
      read(12) re,pou,xl,zl,t,xs
      read(12) nxin,nypin,nzcin,nfzsin
      fltype=0
      dstar=0.0
      read(12,err=1010) fltype,dstar
      if (fltype.ge.6) read(12)
 1010 if (fltype.lt.4.) then
         write(*,9)' The input file is either a temporal simulation or
     &        an old simulation. '
 9       format(/a/)
         fltype=fltype+3
      end if
c
c     Rescale coordinates in the boundary layer case
c
      if (fltype.lt.6.and.fltype.ne.3) dstar=1.
      re=re*dstar
      xl=xl/dstar
      zl=zl/dstar
      t=t/dstar
      xs=xs/dstar
c
c     Skip reading of velocity field
c
      if (header) goto 1100
c
c     Check file info
c
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &    nfzsin.ne.nfzsym) then
         write(*,*) 'input file has a size other than program'
         write(*,*) 'file parameters nx,nyp,nzc,nfzsym',
     &    nxin,nypin,nzcin,nfzsin
         write(*,*) 'program parameters',nx,nyp,nzc,nfzsym
         stop
      end if
      write(*,*) 'Reading data'
c
c     Zero odd-ball
c
      do i=1,m
         do z=1,nzc
            if (z+(i-1)*nzc.eq.nzc*m/100+1) write(*,*) 'Read 1%'
            if (z+(i-1)*nzc.eq.nzc*m/10+1) write(*,*) 'Read 10%'
            if (z+(i-1)*nzc.eq.nzc*m/2+1) write(*,*) 'Read 50%'
            if (z+(i-1)*nzc.eq.nzc*m*9/10+1) write(*,*) 'Read 90%'
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
 1100 continue
      close(unit=12)

      return

      end subroutine rdiscr
