c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine read_xys(xys,xysth,sumw,
     &     re,xl,zl,t,dstar,
     &     fltype,bstart,blength,rlam,spanv,namxys,w,pr,m1,
     &     mhd_n,b0)
c
c     Reads xy statistics from file
c
      implicit none

      include 'par.f'

      integer nnxys,nnxysth,scalarin
      real xys(nx,nyp,nxys)
      real xysth(nx,nyp,nxysth,scalar)
      real re,xl,zl,t,dstar,bstart,blength,rlam,spanv
      logical pou
      integer fltype
      real sumw
      real w(nx,nyp)
      real pr(scalar),m1(scalar)
      character*80 namxys
      character ident
      real mhd_n,b0(3)
      integer x,y,i,nfzsin,nxin,nypin,nzcin,ith
      real xs

      open(unit=19,file=namxys,form='unformatted',status='old')
c
c     Read header
c
      pr =    -999999
      m1 =    -999999
      mhd_n = -999999
      b0 =    -999999.
      read(19,err=1001) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)

      read(19,err=1001) ident
      backspace(19)
      if (ident.eq.'A') then
         read(19) ident,mhd_n,b0
         write(*,*) 'A: ',mhd_n,b0
      else
         write(*,*) 'give mhd_n'
         read(*,*) mhd_n
         write(*,*) 'give b0'
         read(*,*) b0(1),b0(2),b0(3)
      end if



 1001 continue
      read(19) nxin,nypin,nzcin,nfzsin
      read(19) fltype,dstar
      if (fltype.lt.0) read(19) rlam
      if (fltype.ge.6) read(19) bstart,blength,rlam,spanv
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
      read(19,err=1002) sumw,nnxys,nnxysth,scalarin
 1002 continue
      if (nnxys.ne.nxys) then
         write(*,*) 'The file contains ',nnxys,' statistics.'
c         stop
      end if
      if (nnxysth.ne.nxysth) then
         write(*,*) 'The file contains ',nnxysth,' scalar statistics.'
c         stop
      end if
      if (scalar.ne.scalarin) then
         write(*,*) 'The file contains ',scalarin,' scalars.'
         stop
      end if
c
c     Read and add statistics
c
      do i=1,nnxys
        read(19) w
        if (i.eq.43) then
c
c     Max/Min
c
           do y=1,nyp
              do x=1,nx
                 xys(x,y,i)=min(xys(x,y,i),w(x,y)*sumw)
              end do
           end do
        else if (i.eq.44) then
           do y=1,nyp
              do x=1,nx
                 xys(x,y,i)=max(xys(x,y,i),w(x,y)*sumw)
              end do
           end do
        else
c
c     Normal statistics
c
           do y=1,nyp
              do x=1,nx
                 xys(x,y,i)=xys(x,y,i)+w(x,y)*sumw
              end do
           end do
        end if
      end do

      do ith=1,scalar
         do i=1,nnxysth
            read(19) w
c
c     Scalar statistics
c
            do y=1,nyp
               do x=1,nx
                  xysth(x,y,i,ith)=xysth(x,y,i,ith)+w(x,y)*sumw
               end do
            end do
         end do
      end do

      end subroutine read_xys
