c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine dislen(plxy,plxz,tpl,dlen,tlen,t,relth,nzn,xl)
c
c     Calculate the length of the disturbed region
c     this routine recieves one plane at the time
c
      implicit none

      include 'par.f'

      integer nzn,tpl
      real plxy(nx+1,nyp),plxz(nx+1,nz+1),relth,dlen,tlen,t,xl

      integer x,y,z
      real umaxx(nx+1),umax,xfirst,xlast

      tlen=t
c
c     First find the maximum modulus as function of x
c
      umax=0.
      do x=1,nx+1
         umaxx(x)=0.
         if (tpl.eq.1) then
            do y=1,nyp
               umaxx(x)=max(umaxx(x),abs(plxy(x,y)))
            end do
         else
            do z=1,nzn
               umaxx(x)=max(umaxx(x),abs(plxz(x,z)))
            end do
         end if
         umax=max(umaxx(x),umax)
         if (umaxx(x).lt.1E-12) umaxx(x)=0.
      end do
c
c     Now find the first and last points that exceeds the
c     threshold relth of the maximum modulus umax
c     and calculate corresponding coordinate by interpolation
c
      do x=1,nx+1
         if (umaxx(x).gt.umax*relth) goto 1000
      end do
c
c     We should never get here except if all values in the plane =0
c
      dlen=0.
      goto 1010
 1000 continue
c
c     First point that exceeds is x, calculate the coordinate
c     for umaxx(x)=umax*relth
c
      if (x.eq.1) then
         xfirst=xl/real(nx)
      else
         xfirst=(real(x)-(umaxx(x)-umax*relth)/(umaxx(x)-umaxx(x-1)))/
     &        real(nx)*xl
      end if
      do x=nx+1,1,-1
         if (umaxx(x).gt.umax*relth) goto 1020
      end do
c
c     We should never ever get here
c
      write(*,*) 'Internal error in dislen'
      goto 1010
 1020 continue
c
c     Last point that exceeds is x, calculate the coordinate
c     for umaxx(x)=umax*relth
c
      if (x.eq.nx+1) then
         xlast=real(nx+1)/real(nx)*xl
      else
         xlast=(real(x)+(umaxx(x)-umax*relth)/(umaxx(x)-umaxx(x+1)))/
     &        real(nx)*xl
      end if
      dlen=xlast-xfirst
      write(99,*) tlen,dlen
 1010 continue

      return

      end subroutine dislen
