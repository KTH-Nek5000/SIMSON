c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine boxext(vext,cext,u2r,u2i,om2r,om2i,xs,yb,xl,zl)
c
c     Finds the extremum value of each variable and its coordinate in a box
c
      implicit none

      include 'par.f'

      integer yb
      real vext(nyp,2,6),cext(nyp,2,6,2),xs
      real u2r((nxp/2+1),mby,nzd,3),u2i((nxp/2+1),mby,nzd,3)
      real om2r((nxp/2+1),mby,nzd,3),om2i((nxp/2+1),mby,nzd,3)
      real xl,zl

      integer i,y,x,y1
      real xc1(nxp/2+1),xc2(nxp/2+1)
c
c     vext(y,1,i) cext(y,1,1-3,i) the minimum comp i in plane y,
c     Its x,y,z coord
c     vext(y,2,i) cext(y,2,1-3,i) the maximum comp i in plane y,
c     Its x,y,z coord
c     i=1 ... 6 : u,v,w,omegax,omegay,omegaz
c
      do x=1,nxp/2
         xc1(x)=real(2*x-1-nxp/2-1)/real(nxp)*xl+xs
         xc2(x)=real(2*x-nxp/2-1)/real(nxp)*xl+xs
      end do
      do i=1,3
         do y=yb,min(nyp,yb+mby-1)
            y1=y-yb+1
            call plnmin(vext(y,1,i),y1,cext(y,1,i,1),cext(y,1,i,2),
     &           u2r(1,1,1,i),u2i(1,1,1,i),xc1,xc2,zl)
            call plnmax(vext(y,2,i),y1,cext(y,2,i,1),cext(y,2,i,2),
     &           u2r(1,1,1,i),u2i(1,1,1,i),xc1,xc2,zl)
            call plnmin(vext(y,1,i+3),y1,cext(y,1,i+3,1),
     &           cext(y,1,i+3,2),
     &           om2r(1,1,1,i),om2i(1,1,1,i),xc1,xc2,zl)
            call plnmax(vext(y,2,i+3),y1,cext(y,2,i+3,1),
     &           cext(y,2,i+3,2),
     &           om2r(1,1,1,i),om2i(1,1,1,i),xc1,xc2,zl)
         end do
      end do

      end subroutine boxext
