c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rxys(totxys,totxysth,mxys,mxysth,sumw,re,xl,zl,t,dstar,
     &     fltype,bstart,blength,rlam,spanv,namxys,w,pdim,thdim,pr,
     &     scalarind,mhd_n)
c
c     Reads xy statistics from file
c
      include 'par.f'

      integer mxys,mxysth,nscalar,ith,scalarind
      real totxys(nx,nyp,nxys),totxysth(nx,nyp,nxysth,scalar)
      real re,xl,zl,t,dstar,bstart,blength,rlam,spanv
      real pr(scalar),m1(scalar)
      logical pou
      integer fltype
      real sumw
      real w(nx,nyp)
      character*80 namxys
      integer pdim(nxys,2),thdim(nxysth,3)
      real mhd_n,b0(3)

      integer x,y,i,nfzsin,nxin,nypin,nzcin
      real xs,tot
      character ident

      open(unit=19,file=namxys,status='old',form='unformatted')
c
c     Read header
c
      pr = -999999.
      m1 = -999999.
      mhd_n = -999999.
      b0 = -99999.

      read(19,err=1111) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
c
c     For reading two point correlation...
c
      if (re.lt.0.) then
         close(19)
         return
      end if
c
c     Read some optional parameters
c
      read(19,err=1111) ident
      backspace(19)
      if (ident.eq.'A') then
         read(19) ident,mhd_n,b0
         write(*,*) 'A: ',mhd_n,b0
      else
         write(*,*) 'give mhd_n'
         read(*,*) mhd_n
      end if

 1111 continue
      do i=1,scalar
         write(*,*) i,'Pr = ',pr(i),' m1 = ',m1(i)
      end do
      read(19) nxin,nypin,nzcin,nfzsin
      read(19) fltype,dstar
      if (fltype.lt.0) read(19) rlam
      if (fltype.ge.6) read(19) bstart,blength,rlam,spanv
c
c     Rescale from channel to boundary layer scaling
c
      re=re*dstar
      xl=xl/dstar
      zl=zl/dstar
      t=t/dstar
      xs=xs/dstar
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
c
c     Read number of statistic quantities
c
      mxysth=0
      nscalar = 0
      read(19,err=1122) sumw,mxys,mxysth,nscalar
 1122 continue
      if (mxys.gt.nxys) then
         write(*,*)  re,pou,xl,zl,nxin,nypin,nzcin
         write(*,*) 'The file contains ',mxys,' statistics.'
         write(*,*) 'Unable to read more than ',nxys,' .'
         write(*,*) 'To access all statistics :'
         write(*,*) 'increase nxys in par.f and recompile.'
         write(*,*) 'Additional changes necessary are described'
         write(*,*) 'in comments in pxyst.f.'
         mxys=nxys
      end if

      if (scalar.gt.0.and.nscalar.gt.scalar) then
         write(*,*)  re,pou,xl,zl,nxin,nypin,nzcin
         write(*,*) 'The file contains ',nscalar,'scalars.'
         write(*,*) 'Unable to read more than ',scalar,' .'
         write(*,*) 'To access all the scalars :'
         write(*,*) 'increase  scalar in par.f and recompile.'
         write(*,*) 'Additional changes necessary are described'
         write(*,*) 'in comments in pxyst.f.'
         nscalar=scalar
      end if

      if (mxysth.gt.nxysth.and.scalar.gt.0) then
         write(*,*)  re,pou,xl,zl,nxin,nypin,nzcin
         write(*,*) 'The file contains ',mxysth,' scalar statistics.'
         write(*,*) 'Unable to read more than ',nxysth,' .'
         write(*,*) 'To access all statistics :'
         write(*,*) 'increase nxysth in par.f and recompile.'
         write(*,*) 'Additional changes necessary are described'
         write(*,*) 'in comments in pxyst.f.'
         mxysth=nxysth
      end if
c
c     Read statistics
c
      do i=1,mxys
         read(19) w
         do y=1,nyp
            do x=1,nx
               totxys(x,y,i)=w(x,y)
            end do
         end do
      end do
c
c     Rescale from channel to boundary layer scaling
c
      do i=1,mxys
         do y=1,nyp
            do x=1,nx
               totxys(x,y,i)=totxys(x,y,i)*dstar**(-pdim(i,2))
            end do
         end do
      end do
c
c     IMPORTANT BUG FIX
c     (20070317 Philipp Schlatter)
c     Some of the quantities read from bla are read in with a different
c     length scale than they are used in pxyst
c     (i.e. than what is stored in pdim(:,2)).
c     Example: position 10. from bla: <w^2>        ==> V^2/L^2
c                           in pxyst: sqrt(<w'w'>) ==> V  /L
c
c     This is fixed here for the vorticity rms
c
c     Important: CHECK ALL OTHER QUANTITIES, then we can remove
c     this warning.
c
      do y=1,nyp
         do x=1,nx
            totxys(x,y,10)=totxys(x,y,10)*dstar
            totxys(x,y,11)=totxys(x,y,11)*dstar
            totxys(x,y,12)=totxys(x,y,12)*dstar

            if (nxys.ge.84)
     &           totxys(x,y,84)=totxys(x,y,84)/dstar
         end do
      end do



c
c     Read statistics for scalar
c
      if (scalar.gt.0) then
         do ith=1,nscalar
            do i=1,mxysth
               read(19) w
               do y=1,nyp
                  do x=1,nx
                     totxysth(x,y,i,ith)=w(x,y)
                  end do
               end do
            end do
         end do
c
c     Rescale from channel to boundary layer scaling
c
         do ith=1,nscalar
            do i=1,mxysth
               do y=1,nyp
                  do x=1,nx
                     totxysth(x,y,i,ith)=totxysth(x,y,i,ith)*
     &                    dstar**(-thdim(i,2))
                  end do
               end do
            end do
         end do
      end if



      scalarind = 1
      if (scalar.gt.1) then
         write(*,*) 'which scalar to read? max=',scalar
         read(*,*) scalarind
         write(*,*) 'pr     =',pr(scalarind)
         write(*,*) 'm1     =',m1(scalarind)
      end if




c
c     For umax/umin...
c
      if (mxys.ge.43) then
         do i=43,44
            do y=1,nyp
               do x=1,nx
                  totxys(x,y,i)=totxys(x,y,i)*sumw
               end do
            end do
         end do
      end if

      write(*,*) 'fltype=',fltype

      if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.-3) then
         write(*,*) 'Channel/Couette: average over x? (yes 1, no 0)'
         read(*,*) i
         if (i.eq.1) then
            do i=1,mxys
               do y=1,nyp
                  tot = 0
                  do x=1,nx
                     tot =tot + totxys(x,y,i)
                  end do
                  do x=1,nx
                     totxys(x,y,i)=tot/nx
                  end do
               end do
            end do
            do ith=1,nscalar
               do i=1,mxysth
                  do y=1,nyp
                     tot = 0
                     do x=1,nx
                        tot =tot + totxysth(x,y,i,ith)
                     end do
                     do x=1,nx
                        totxysth(x,y,i,ith)=tot/nx
                     end do
                  end do
               end do
            end do
         end if

      end if

      end subroutine rxys
