c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wxys(xys,xysth,sumw,re,xl,zl,t,dstar,fltype,
     &     bstart,blength,rlam,spanv,namxys,w,pr,m1,mhd_n,b0)
c
c     Writes xy statistics to file
c
      implicit none

      include 'par.f'

      real xys(nx,nyp,nxys)
      real xysth(nx,nyp,nxysth,scalar)
      real re,xl,zl,t,dstar,bstart,blength,rlam,spanv
      integer fltype
      real sumw,pr(scalar),m1(scalar)
      real w(nx,nyp)
      character*32 namxys
      real mhd_n,b0(3)
      integer x,y,i,ith

      open(unit=19,file=namxys,form='unformatted')
      rewind(19)
c
c     Write header
c
      write(19) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar)
      write(19) 'A',mhd_n,b0
      write(19) nx,nyp,nzc,nfzsym
      write(19) fltype,dstar
      if (fltype.lt.0) write(19) rlam
      if (fltype.ge.6) write(19) bstart,blength,rlam,spanv
c
c     Write number of statistic quantities
c
      write(19) sumw,nxys,nxysth,scalar
c
c     Write statistics
c
      do i=1,nxys
c
c     Divide accumulated statistics by sum of weights
c
         do y=1,nyp
            do x=1,nx
               w(x,y)=xys(x,y,i)*(1./sumw)
            end do
         end do
         write(19) w
      end do
      do ith=1,scalar
         do i=1,nxysth
c
c     Divide accumulated statistics by sum of weights
c
            do y=1,nyp
               do x=1,nx
                  w(x,y)=xysth(x,y,i,ith)*(1./sumw)
               end do
            end do
            write(19) w
         end do
      end do
      close(unit=19)

      end subroutine wxys
