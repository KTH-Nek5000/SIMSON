c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rseq(uxys,uxzs,uyzs,re,fltype,
     &   xl,zl,dstar,t,xs,namnin,uxy,uxz,uyz,
     &   xsa,ta,tpl,ivar,cpl,mpl,npl,nzn,
     &   umin,umax)
c
c     Reads file namnin and puts into uxys or uxys
c
      implicit none

      include 'par.f'

      integer npl,nzn
      real uxzs(nx,nzn,npl),uxys(nx,nyp,npl),uyzs(nyp,nzn,npl)
      real xsa(npl),ta(npl)
      real prex(nx+15),prey(2*nyp+15)
      real prez(2*nz+15),prezr(nz+15),pres(2*nz+20),prea(2*nz+15)
      real cpl,dstar,t0
      integer tpl,ivar,mpl
      integer fltype
      character*80 namnin
      real re,xl,zl,t,xs
      logical pou
      real uyz(nyp,nzn)
      real uxz(nx,nzn),uxy(nx,nyp)
      integer imax, jmax
      real umax,umin

      integer x,y,z,i,ir
      integer nxin,nypin,nzcin,nfzsin
      integer inunit
      character*3 var(6)
c
c     Preprocess transforms
c
      call vrffti(nx,prex,0)
      call vcosti(nyp,prey,0)
      if (nfzsym.eq.0) then
         call vcffti(nz,prez,0)
         if (nz.ne.1) call vrffti(nz,prezr,0)
      else
         call vcosti(nz/2+1,pres,0)
         call vsinti(nz/2-1,prea,0)
      end if
      inunit=12
      open(unit=12,file=namnin,status='old',form='unformatted')
c
c     Read file header
c
      read(12) re,pou,xl,zl,t,xs
      read(12) nxin,nypin,nzcin,nfzsin
      read(12,err=1010,end=1010) tpl,ivar,cpl,fltype,dstar
      if (fltype.lt.0.or.(fltype.ge.1.and.fltype.le.9)) goto 1015
 1010 continue
      if (pou) fltype=1
      if (.not.pou) fltype=2
 1015 continue
      write(*,*) 'Start time: ',t/dstar
      if (fltype.lt.-2.or.fltype.gt.9.or.fltype.eq.0) then
         write(*,*) 'Invalid flow type ',fltype
      end if
c
c     Rescale coordinates in the boundary layer case
c
      if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.4.or.fltype.eq.5)
     &     dstar=1.
      re=re*dstar
      xl=xl/dstar
      zl=zl/dstar
      t=t/dstar

      if (1.eq.0) then
         t0=t
         t=0.
      else
         t0 = 0.
      end if

      xs=xs/dstar
c
c     Write fileinfo
c
      write(*,9010)xl,zl,t,re
 9010 format(' xl= ',F8.2,' zl= ',f8.2,' t= ',F8.1,' re= ',F8.1)
      if (fltype.eq.3) write(*,9015) 2./dstar
 9015 format(' yl= ',F8.2)
      if (fltype.eq.-2) write(*,*) 'Falkner-Skan-Cooke boundary layer'
      if (fltype.eq.-1) write(*,*) 'Falkner-Skan boundary layer'
      if (fltype.eq.1) write(*,*) 'Poiseuille flow'
      if (fltype.eq.2) write(*,*) 'Couette flow'
      if (fltype.eq.3) write(*,*) 'Blasius boundary layer flow'
      if (fltype.eq.4) write(*,*) 'spatial Poiseuille flow'
      if (fltype.eq.5) write(*,*) 'spatial Couette flow'
      if (fltype.eq.6) write(*,*) 'spatial Blasius boundary layer'
      if (fltype.eq.7) write(*,*) 'spatial Falkner-Skan boundary layer'
      if (fltype.eq.8)
     &     write(*,*) 'spatial Falkner-Skan-Cooke boundary layer'
      if (fltype.eq.9) write(*,*) 'spatial parallel boundary layer'
      var(1)='u'
      var(2)='v'
      var(3)='w'
      var(4)='omx'
      var(5)='omy'
      var(6)='omz'
      if (tpl.eq.1) write(*,9020) var(ivar),cpl
 9020 format(' xy-plane of ',a,' at z= ',F7.2)
      if (tpl.eq.2) write(*,9030) var(ivar),cpl
 9030 format(' xz-plane of ',a,' at y= ',F6.3)
      if (tpl.eq.3) write(*,9040) var(ivar),cpl
 9040 format(' yz-plane of ',a,' at x= ',F7.2)
c
c     Check file info
c
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &    nfzsin.ne.nfzsym) then
         write(*,*) 'Input file has a size other than program'
         write(*,*) 'File parameters nx,nyp,nzc,nfzsym',
     &    nxin,nypin,nzcin,nfzsin
         write(*,*) 'Program parameters',nx,nyp,nzc,nfzsym
         stop
      end if
      umax = -1000.
      umin = 1000.
      do i=1,npl
         read(12,end=1000,err=1000) ta(i),xsa(i)

         if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6) then
            ta(i)=ta(i)/dstar-t0
            xsa(i)=xsa(i)/dstar
         end if
         if (i.gt.1) then
            if (ta(i).lt.ta(i-1)) goto 1000
         else
c               if (ta(1).lt.t) goto 1000
         end if
         if (tpl.eq.1) read(12,end=1020,err=1020) uxy
         if (tpl.eq.2) read(12,end=1020,err=1020) uxz
         if (tpl.eq.3) read(12,end=1020,err=1020) uyz

         if (tpl.eq.1) then
            do imax=1,nx
               do jmax=1,nyp
                  if (uxy(imax,jmax).gt.100000000) then
                     write(*,*) uxy(imax,jmax)
                  end if
                  if (uxy(imax,jmax).gt.umax) then
                     umax = uxy(imax,jmax)
                  end if
                  if (uxy(imax,jmax).lt.umin) then
                     umin = uxy(imax,jmax)
                  end if
               end do
            end do
         end if
         if (tpl.eq.2) then
            do imax=1,nx
               do jmax=1,nzn
                  if (uxy(imax,jmax).gt.umax) then
                     umax = uxy(imax,jmax)
                  end if
                  if (uxy(imax,jmax).lt.umin) then
                     umin = uxy(imax,jmax)
                  end if
               end do
            end do
         end if
         if (tpl.eq.3) then
            do imax=1,nyp
               do jmax=1,nzn
                  if (uxy(imax,jmax).gt.umax) then
                     umax = uxy(imax,jmax)
                  end if
                  if (uxy(imax,jmax).lt.umin) then
                     umin = uxy(imax,jmax)
                  end if
               end do
            end do
         end if
         if (tpl.eq.1) then
            do y=1,nyp
               do x=1,nx
                  uxys(x,y,i)=uxy(x,y)
               end do
            end do
         end if
         if (tpl.eq.2) then
            do z=1,nzn
               do x=1,nx
                  uxzs(x,z,i)=uxz(x,z)
               end do
            end do
         end if
         if (tpl.eq.3) then
            do y=1,nyp
               do z=1,nzn
                  uyzs(y,z,i)=uyz(y,z)
               end do
            end do
         end if
         ir=nint(10.**int(log10(real(i))))
         if (i.eq.ir.or.i.eq.5*ir)
     &        write(*,*) 'read ',i,'f plane(s), t= ',ta(i),
     &        '(',ta(i)+t0,')'
      end do
      write(*,*) 'Unable to read the whole file.'
      write(*,9000) npl
 9000 format('Read ',I5,' planes.')
      write(*,*) 'To read more, increase npl in MAIN and recompile'
      mpl=npl
      close(unit=12)
      return
 1020 write(*,*) 'Read incomplete plane'
 1000 mpl=i-1
      close(unit=12)

      return

      end subroutine rseq
