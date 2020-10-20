c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program cmp
c
c     Calculates linear combinations of different fields
c
      implicit none

      include 'par.f'
c
c     Main storage
c
      real ur(memnx,memny,memnz,3,2),ui(memnx,memny,memnz,3,2)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real urx(nx)
      real fvr(nx/2,mbz,nyp),fvi(nx/2,mbz,nyp)
      real w3(nx/2,mbz,nyp),updtime
c
c     Amplitudes etc
c
      real alfa(nx/2,mbz),beta(nz),wint(nyp)
      real amp(nyp,20)

      character*80 namnin,ynchar
      character*80 namnut
      character*1 ans, ans2

      real eta(nyp),xl,zl,xs,t
      real xshift,zshift
c
c     Flow
c
      integer fltype
      real re,mflux,dstar
      real u0low,u0upp,w0low,w0upp
      real bstart,blength,rlam,spanv

      integer i,x,y,z,yb,n,lfield
      real c,d
      real pi
      parameter (pi = 3.1415926535897932385)
      logical varsiz

      real prey(nyp*2+15)
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)

      integer ll,j,k
      real ca(nx/2),sa(nx/2),cx(nx/2),sx(nx/2)
      real argz,argx,hr

      varsiz=.true.
c
c     Initialize x-transform
c
      call vrffti(nxp,prex,0)
c
c     Initialize y-transform
c
      call vcosti(nyp,prey,0)
c
c     Initialize z-transform
c
      if (nfzsym.eq.0) then
         call vcffti(nzp,prez,0)
      else
         call vcosti(nzst,pres,0)
         call vsinti(nzat,prea,0)
      end if
c
c     Init
c
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         cmp $Rev$   '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)

      do y=1,nyp
         eta(y)=cos(pi*float(y-1)/float(nyp-1))
      end do

      write(*,*) 'Number of fields to be manipulated ? '
      read(*,*) lfield
      do i=1,lfield
         write(*,*) 'Name of file'
         read(*,9000) namnin
 9000    format(a80)
         n=min(i,2)
         call rdisca(ur,ui,u0low,u0upp,w0low,w0upp,re,mflux,
     &        xl,zl,t,xs,dstar,fltype,bstart,blength,rlam,spanv,
     &        namnin,varsiz,3,n,fvr,fvi,w3,urx)
         write(*,9010)xl/dstar,zl/dstar,t/dstar,re*dstar
 9010    format(' xl= ',F8.2,' zl= ',f8.2,' t= ',F8.1,' re= ',F8.1)
         if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0)
     &        write(*,9015) 2./dstar
 9015    format(' yl= ',F8.2)
         if (fltype.eq.-2) write(*,*)'Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.-1) write(*,*) 'Falkner-Skan boundary layer'
         if (fltype.eq.1) write(*,*) 'Poiseuille flow'
         if (fltype.eq.2) write(*,*) 'Couette flow'
         if (fltype.eq.3) write(*,*) 'Temporal Blasius boundary layer'
         if (fltype.eq.4) write(*,*) 'Spatial Poiseuille flow'
         if (fltype.eq.5) write(*,*) 'Spatial Couette flow'
         if (fltype.eq.6) write(*,*) 'Spatial Blasius boundary layer'
         if (fltype.eq.7)write(*,*)'Spatial Falkner-Skan boundary layer'
         if (fltype.eq.8)
     &        write(*,*) 'Spatial Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.9) write(*,*) 'Spatial parallel boundary layer'
         write(*,*) 'Give multiplier for this disturbance'
         write(*,*) '-99 gives power of field'
         read(*,*) c
         if (c.eq.-99) then
            write(*,*) 'Give power'
            read(*,*) d
            if (abs(d).lt.1) then
               write(*,*) 'Numbers between 0 and -1e-3 will be put to 0'
            end if
            call power(n,d,ur,ui,eta,fltype,
     &           prex,prez,pres,prea)
         else
            call add(n,c,ur,ui,eta,fltype)
         end if
      end do
      write(*,*) 'Force symmetry properties? (y/n)'
      read(*,'(a1)') ans
      if (ans.ne.'n'.and.ans.ne.'N') then
         write(*,*)'Make the field symmetric (s) or anti-symmetric (a)?'
         read(*,'(a1)') ans2
         if (ans.eq.'s') then
            write(*,*) 'Making the field symmetric w.r.t. z-axis...'
            do yb=1,nyp,mby
               call symm_z(ur,ui,yb,u2r,u2i,1)
            end do
         else
            write(*,*)'Making the field anti-symmetric w.r.t. z-axis...'
            do yb=1,nyp,mby
               call symm_z(ur,ui,yb,u2r,u2i,2)
            end do
         end if
      end if

      write(*,*) 'Shift in x/z plane? (y/n)'
      read(*,'(a1)') ans
      if (ans.ne.'n'.and.ans.ne.'N') then
         write(*,*) 'Give shifts in x and z'
         read(*,*) xshift,zshift

         do i=1,nx/2
            alfa(i,1) = 2.*pi/xl*real(i-1)
         end do
         beta(1) = 0.
         do k=2,nz/2+1
            beta(k) = 2.*pi/zl*real(k-1)
            beta(nz+2-k) = -2.*pi/zl*real(k-1)
         end do

         do i=1,nx/2
            argx = -xshift*alfa(i,1)
            cx(i) = cos(argx)
            sx(i) = sin(argx)
         end do

         do k=1,nz
            argz = -zshift*beta(k)
            do i=1,nx/2
               ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz)
               sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz)
            end do
            do ll=1,3+scalar
               do j=1,nyp
                  do i=1,nx/2
                     hr=ur(i,j,k,ll,1)*ca(i)-ui(i,j,k,ll,1)*sa(i)
                     ui(i,j,k,ll,1)=
     &                    ui(i,j,k,ll,1)*ca(i)+ur(i,j,k,ll,1)*sa(i)
                     ur(i,j,k,ll,1)=hr
                  end do
               end do
            end do
         end do

      end if


      write(*,*) 'Calculate disturbance amplitude ? (y/n)'
      read(*,'(a1)') ans
      if (ans.ne.'n'.and.ans.ne.'N') then
c
c     Geometrics
c
         do i=1,mbz
            do x=1,nx/2
               alfa(x,i)=2.*pi/xl*float(x-1)
            end do
         end do
         beta(1)=0.
         do z=2,nz/2+1
            beta(z)=2.*pi/zl*float(z-1)
            beta(nz+2-z)=-2.*pi/zl*float(z-1)
         end do
         do y=1,nyp
            eta(y)=cos(pi*float(y-1)/float(nyp-1))
         end do
         do y=1,nyp
            wint(y)=-.5
            do i=1,(nyp+1)/2-2
               wint(y) = wint(y) +cos(pi*float(2*i*(y-1))/float(nyp-1))/
     &              float((2*i)**2-1)
            end do
            wint(y) =wint(y) + .5*cos(pi*float(y-1))/float((nyp-1)**2-1)
            wint(y) = -4./float(nyp-1)*wint(y)
            if (y.eq.1.or.y.eq.nyp) wint(y)=wint(y)*.5
         end do
c
c     Find omx and omz
c
        call comxz(ur,ui,alfa,beta,prey)
c
c     Accumulate amplitudes
c
        do yb=1,nyp,mby
           if (yb.eq.nyp/100+1) then
              write(*,*)
              write(*,*) 'Processed 1%'
           end if
           if (yb.eq.nyp/10+1) write(*,*) 'Processed 10%'
           if (yb.eq.nyp/2+1) write(*,*) 'Processed 50%'
           if (yb.eq.nyp*9/10+1) write(*,*) 'Processed 90%'
           call namp(amp,ur,ui,yb,alfa,beta,u2r,u2i,om2r,om2i,
     &          vext,cext,xs,xl,zl,prex,prez,pres,prea)
        end do
        call wcamp(t,amp,eta,wint,re,dstar,fltype,u0upp,u0low)
        call wcext(vext,cext,u0low,dstar,eta,fltype,xl)
      end if
      write(*,*) 'Give name for output file'
      read(*,9000) namnut
      write(*,*) 'Change time stamp (y/n)'
      read(*,9000) ynchar
      if (ynchar.eq.'y'.or.ynchar.eq.'Y') then
         write(*,*) 'Give new time stamp'
         read(*,*) updtime
         t=updtime
      end if
c
c     Generate accuracy requirement
c
      call wdisca(ur,ui,re,xl,zl,t,xs,dstar,fltype,bstart,blength,
     &     rlam,spanv,namnut,3,fvr,fvi,urx)
      stop

      end program cmp
