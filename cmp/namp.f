c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine namp(amp,ur,ui,yb,alfa,beta,u2r,u2i,om2r,om2i,
     &     vext,cext,xs,xl,zl,prex,prez,pres,prea)
c
c     Accumulates amp and extremum values for a xz-box
c     i.e. nby xz-planes
c
      implicit none

      include 'par.f'

      integer mwave
      parameter (mwave=1)
      integer yb,kx(mwave),kz(mwave)
      complex campw(nyp,4,mwave)
      real amp(nyp,20)
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real alfa(nx/2*mby),beta(nz),xl,zl
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real vext(nyp,2,6),cext(nyp,2,6,2),xs
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      integer i,xy,z,zp,y,x,npl
      logical sym
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)
c
      do i=1,3
        call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui)
      end do
      call getxz(om2r(1,1,1),om2i(1,1,1),yb,4,1,ur,ui)
      call getxz(om2r(1,1,3),om2i(1,1,3),yb,5,1,ur,ui)
c
      do z=1,nz/2
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               xy=x+(y-yb)*(nxp/2+1)
               om2r(xy,z,2)=-beta(z)*u2i(xy,z,1)+alfa(x)*u2i(xy,z,3)
               om2i(xy,z,2)=beta(z)*u2r(xy,z,1)-alfa(x)*u2r(xy,z,3)
            end do
         end do
      end do
      if (nfzsym.eq.0) then
         do z=nz/2+1,nz
            zp=nzp-nz+z
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  xy=x+(y-yb)*(nxp/2+1)
                  om2r(xy,zp,2)=-beta(z)*u2i(xy,zp,1)
     &                          +alfa(x)*u2i(xy,zp,3)
                  om2i(xy,zp,2)= beta(z)*u2r(xy,zp,1)
     &                          -alfa(x)*u2r(xy,zp,3)
               end do
            end do
         end do
      end if
      do z=nz/2+1,min(nzpc,nzp+1-nz/2)
         do xy=1,(nxp/2+1)*min(mby,nyp-yb+1)
            om2r(xy,z,2)=0.0
            om2i(xy,z,2)=0.0
         end do
      end do
      do z=1,nzpc
         do y=yb,min(nyp,yb+mby-1)
            do x=nx/2+1,nxp/2+1
               xy=x+(y-yb)*(nxp/2+1)
               om2r(xy,z,2)=0.0
               om2i(xy,z,2)=0.0
            end do
         end do
      end do
c
      call boxamp(amp,campw,kx,kz,0,
     &     u2r,u2i,om2r,om2i,yb,alfa,beta)
      npl=min(mby,nyp-yb+1)
      do i=1,3
         sym=i.le.2
         call fft2db(u2r(1,1,i),u2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
         sym=i.eq.3
         call fft2db(om2r(1,1,i),om2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
      end do
      call boxext(vext,cext,u2r,u2i,om2r,om2i,xs,yb,xl,zl)

      return

      end subroutine namp
