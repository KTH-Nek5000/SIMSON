c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program xys_add
c
c     Adding xy statistics for plotting
c
      implicit none
c
      include 'par.f'
c
      integer z,scalarin
      integer nxin,nypin,nzcin,corrnum,nfzsymin

      real corrx(mcorr),corry(mcorr),dum
      real corrxx(mcorr),corryy(mcorr)
      integer corrxi(mcorr),corryi(mcorr)
      real corr(nzc,mcorr)
      real corr1(nzc)
      real corr_ms(4,mcorr)
      real corr_ms1(4,mcorr)
      integer corrt(mcorr)
      real mhd_n,b0(3)
      real xys(nx,nyp,nxys),xysth(nx,nyp,nxysth,scalar)
      real w(nx,nyp)
c
      character*80 namxys
c
c     Geometrics
c
      real xl,zl,re,t,dstar,pr(scalar),m1(scalar)
c
c     Flow
c
      integer fltype,nn
      logical pou
      real blength,bstart,rlam,spanv,sumw,sum_total
      integer i
      logical corrf
      character*1 ans

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                        xys_add $Rev$'//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)


      corr_ms = 0.
      corr = 0.
      xys = 0.
      xysth = 0.
      sum_total=0.

      corrf=.false.
      write(*,*) 'Compiled statistics:          ',nxys
      write(*,*) 'Compiled statistics (scalar): ',nxysth
      write(*,*) 'Max. correlations:            ',mcorr
      write(*,*)
      write(*,*) 'Add statistic files (y/n)'
      write(*,*) 'or correlation files'
      read(*,'(a1)') ans
      if (ans.eq.'n') then
         corrf=.true.
      end if

 666  write(*,*) 'Do you want to add another file (y/n)'
      read(*,'(a1)') ans
      if (ans.ne.'N'.and.ans.ne.'n') then
         write(*,*) 'Give filename of statistics'
         read(*,'(a)') namxys
         if (corrf) then
c
c     Read correlations
c
            open(unit=19,file=namxys,form='unformatted',status='old')
            read(19,err=1001) re,pou,xl,zl,t,dum,
     &           (pr(i),m1(i),i=1,scalar)
 1001       continue
            re=-re
            read(19) nxin,nypin,nzcin,nfzsymin
            read(19) fltype,dstar
            if (fltype.ge.6) read(19) bstart,blength
c
c     Check file info
c
            if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &           nfzsymin.ne.nfzsym) then
               write(*,*) 'Input file has a size other than program'
               write(*,*) 'File parameters nx,nyp,nzc,nfzsym',
     &              nxin,nypin,nzcin,nfzsymin
               write(*,*) 'program parameters',nx,nyp,nzc,nfzsym
               stop
            end if
            read(19) sumw,scalarin
            read(19) corrnum
            nn=nzc
            do i=1,corrnum
               read(19) corrx(i),corry(i),
     &              corrxi(i),corryi(i),
     &              corrxx(i),corryy(i),corrt(i)
               read(19) (corr_ms1(z,i),z=1,4)
               do z=1,4
                  corr_ms(z,i) = corr_ms(z,i) + corr_ms1(z,i)*sumw
               end do
               read(19) (corr1(z),z=1,nn)
               do z=1,nn
                  corr(z,i)=corr(z,i)+corr1(z)*sumw
               end do
            end do
            close(19)
         else
c
c     Read normal statistics
c
            call  read_xys(xys,xysth,sumw,re,xl,zl,t,dstar,
     &           fltype,bstart,blength,rlam,spanv,namxys,w,pr,m1,
     &           mhd_n,b0)
         end if

         sum_total=sum_total+sumw
         if (corrf) then
            write(*,*) 'Read ',corrnum,' correlations'
            write(*,*) 'Weighted averaging time = ',sumw
         else
            write(*,*) 'Read ',nxys,' statistics'
            write(*,*) 'Read ',nxysth,' scalar statistics'
            write(*,*) 'Weighted averaging time = ',sumw/dstar
         end if
         goto 666
      else
         write(*,*) 'Give filename to store total statistics'
         read(*,'(a)') namxys
         sumw=sum_total
         if (corrf) then
c
c     Write correlations
c
            open(unit=19,file=namxys,form='unformatted')
            write(19) -re,.false.,xl,zl,t,dum,(pr(i),m1(i),i=1,scalar)
            write(19) nx,nyp,nzc,nfzsym
            write(19) fltype,dstar
            if (fltype.ge.6) write(19) bstart,blength
            write(19) sumw,scalarin
            write(19) corrnum
            nn=nzc
            do i=1,corrnum
               write(19) corrx(i),corry(i),
     &              corrxi(i),corryi(i),
     &              corrxx(i),corryy(i),corrt(i)
               write(19) (corr_ms(z,i)*(1./sumw),z=1,4)
               write(19) (corr(z,i)*(1./sumw),z=1,nn)
            end do
            close(19)
         else
c
c     Write normal statistics
c
            call wxys(xys,xysth,sumw,re,xl,zl,t,dstar,fltype,
     &           bstart,blength,rlam,spanv,namxys,w,pr,m1,mhd_n,b0)
            write(*,*) 'Wrote ',nxys,' statistics.'
            write(*,*) 'Wrote ',nxysth,' scalar statistics.'
         end if
      end if

      end program xys_add
