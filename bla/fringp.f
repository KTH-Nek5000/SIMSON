c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine fringp(om2r,om2i,u2r,u2i,xsc,xl,yb,
     &     fstart,fend,bu1,bu2,fring1,fring2)
c
c     Fringe function used for the pressure solver.
c     Forces only to the base flow
c
      implicit none

      include 'par.f'

      integer yb
      real xsc,xl
      real fstart,fend
      real bu1(nxp/2+1,nyp,3),bu2(nxp/2+1,nyp,3)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real om2r(nxp/2+1,mby,nzd,3),om2i(nxp/2+1,mby,nzd,3)

      integer x,y,z,y1,xstart,xend,xendc,xst,xen,i,imin
      real fring1(nxp/2),fring2(nxp/2)
      real fstc,fenc
      real r
      parameter(r=0.499999999999999999)

      fstc=fstart+int((xsc-fstart)/xl+r)*xl
      fenc=fend+int((xsc-fend)/xl+r)*xl
      xstart=int(real(nxp/2)*(fstc-xsc+xl/2.)/xl)+1
      xend=int(real(nxp/2)*(fenc-xsc+xl/2.)/xl)+1

      if (xstart.gt.xend) then
         imin=0
         xendc=nxp/2
      else
         imin=1
         xendc=1
      end if

      do i=imin,1
         xst=xstart**i
         xen=max(xend,xendc**i)
c
c     No wave forcing, damping to base flow in the fringe region
c
         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               y1=y+yb-1
               do x=xst,xen
                  om2r(x,y,z,1)=om2r(x,y,z,1)-
     &                 fring1(x)*(u2r(x,y,z,1)-bu1(x,y1,1))
                  om2i(x,y,z,1)=om2i(x,y,z,1)-
     &                 fring2(x)*(u2i(x,y,z,1)-bu2(x,y1,1))
                  om2r(x,y,z,2)=om2r(x,y,z,2)-
     &                 fring1(x)*(u2r(x,y,z,2)-bu1(x,y1,2))
                  om2i(x,y,z,2)=om2i(x,y,z,2)-
     &                 fring2(x)*(u2i(x,y,z,2)-bu2(x,y1,2))
                  om2r(x,y,z,3)=om2r(x,y,z,3)-
     &                 fring1(x)*(u2r(x,y,z,3)-bu1(x,y1,3))
                  om2i(x,y,z,3)=om2i(x,y,z,3)-
     &                 fring2(x)*(u2i(x,y,z,3)-bu2(x,y1,3))
               end do
            end do
         end do
      end do

      end subroutine fringp
