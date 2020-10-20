c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine nonlinp(ur,ui,xsc,yb,rot,spat,fstart,fend,
     &     bu1,bu2,alfa,beta,xl,prex,prez,pres,prea,
     &     u2r,u2i,om2r,om2i,wr,wi,my_node,realg1,realg2,
     &     fring1,fring2)

      implicit none

      include 'par.f'

      logical spat
      integer yb
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real rot,alfa(nx/2*mby),beta(nz)
      real fstart,fend
      real xl,xsc
      real bu1(nxp/2+1,nyp,3),bu2(nxp/2+1,nyp,3)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)
      real fring1(nxp/2),fring2(nxp/2)

      logical sym
      integer i,xy,z,zp,nxy,y,y1,x,nxp2,npl,npp
      real h1u,h2u,h1e,h2e
c
c     MPI
c
      integer my_node
      integer realg1,realg2


      npp=min(nyp,yb+mby-1)
      npl=min(mby,nyp-yb+1)
      nxp2=nxp/2+1
      nxy=nxp2*npl
c
c     Get the velocities (padded)
c
      if (nproc.eq.1) then
         do i=1,3
            call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui)
         end do
c
c     Get the streamwise and spanwise vorticity
c
         call getxz(om2r(1,1,1),om2i(1,1,1),yb,4,1,ur,ui)
         call getxz(om2r(1,1,3),om2i(1,1,3),yb,5,1,ur,ui)
      else
#ifdef MPI
         do i=1,3
            call getpxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui,
     &           realg1,realg2,my_node)
         end do
         call getpxz(om2r(1,1,1),om2i(1,1,1),yb,4,1,ur,ui,
     &        realg1,realg2,my_node)
         call getpxz(om2r(1,1,3),om2i(1,1,3),yb,5,1,ur,ui,
     &        realg1,realg2,my_node)
#endif
      end if
c
c     Compute the wall-normal vorticity
c
      do z=1,nz/2
         do y=yb,npp
            y1=(y-yb)*nxp2
            do x=1,nx/2
               xy=x+y1
               om2r(xy,z,2)=-beta(z)*u2i(xy,z,1)+alfa(x)*u2i(xy,z,3)
               om2i(xy,z,2)= beta(z)*u2r(xy,z,1)-alfa(x)*u2r(xy,z,3)
            end do
         end do
      end do

      if (nfzsym.eq.0) then
c
c     Do it also for the upper half in the z direction
c
         do z=nz/2+1,nz
            zp=nzp-nz+z
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1
                  om2r(xy,zp,2) = -beta(z)*u2i(xy,zp,1) +
     &                 alfa(x)*u2i(xy,zp,3)
                  om2i(xy,zp,2) = beta(z)*u2r(xy,zp,1) -
     &                 alfa(x)*u2r(xy,zp,3)
               end do
            end do
         end do
      end if
c
c     And pad/oddball
c
      do z=(nz+1)/2+1,min(nzpc,nzp+1-(nz+1)/2)
         do xy=1,nxy
            om2r(xy,z,2)=0.0
            om2i(xy,z,2)=0.0
         end do
      end do
      do z=1,nzpc
         do y=yb,npp
            y1=(y-yb)*nxp2
            do x=nx/2+1,nxp2
               xy=x+y1
               om2r(xy,z,2)=0.0
               om2i(xy,z,2)=0.0
            end do
         end do
      end do
c
c     Backward Fourier transform (Fourier synthesis)
c
      do i=1,3
         sym=i.le.2
         call fft2db(u2r(1,1,i),u2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
         sym=i.eq.3
         call fft2db(om2r(1,1,i),om2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
      end do
c
c     We are now in the physical space
c     Calculate the full advection in physical space with or without
c     dealiasing.
c
c     H1,H2,H3 and 2*E=u^2+v^2+w^2
c
      do z=1,nzpc
         do xy=1,nxy
            h1u=u2r(xy,z,2)*(om2r(xy,z,3)+
     &           2.*rot)-u2r(xy,z,3)*om2r(xy,z,2)
            h2u=u2r(xy,z,3)*om2r(xy,z,1)-
     &           u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)

            om2r(xy,z,3)=u2r(xy,z,1)*om2r(xy,z,2)-
     &           u2r(xy,z,2)*om2r(xy,z,1)
            om2r(xy,z,1)=h1u
            om2r(xy,z,2)=h2u

            h1e=u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)-
     &           u2i(xy,z,3)*om2i(xy,z,2)
            h2e=u2i(xy,z,3)*om2i(xy,z,1)-
     &           u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)

            om2i(xy,z,3)=u2i(xy,z,1)*om2i(xy,z,2)-
     &           u2i(xy,z,2)*om2i(xy,z,1)
            om2i(xy,z,1)=h1e
            om2i(xy,z,2)=h2e

            wr(xy,z)=u2r(xy,z,1)*u2r(xy,z,1)+
     &           u2r(xy,z,2)*u2r(xy,z,2)+
     &           u2r(xy,z,3)*u2r(xy,z,3)
            wi(xy,z)=u2i(xy,z,1)*u2i(xy,z,1)+
     &           u2i(xy,z,2)*u2i(xy,z,2)+
     &           u2i(xy,z,3)*u2i(xy,z,3)
         end do
      end do
c
c     Fringe region
c
      if (spat) then
         call fringp(om2r,om2i,u2r,u2i,xsc,xl,yb,
     &        fstart,fend,bu1,bu2,fring1,fring2)
      end if
c
c     Forward Fourier transform (Fourier analysis)
c
      do i=1,3
         sym=i.le.2
         call fft2df(om2r(1,1,i),om2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,u2r,u2i)
      end do
      call fft2df(wr,wi,.true.,npl,prex,prez,pres,prea,u2r,u2i)
c
c     Remove the oddball modes in z on the fine grid
c     (if necessary, i.e. nzp even)
c
      if (mod(nzp,2).eq.0) then
         do i=1,3
            do xy=1,nxy
               om2r(xy,nzp/2+1,i)=0.0
               om2i(xy,nzp/2+1,i)=0.0
            end do
         end do
         do xy=1,nxy
            wr(xy,nzp/2+1)=0.0
            wi(xy,nzp/2+1)=0.0
         end do
      end if
c
c     Calculate d/dx H1 + d/dz H3
c     Real part goes to imaginary part and vice versa
c
      do z=1,nz/2
         do y=yb,npp
            y1=(y-yb)*nxp2
            do x=1,nx/2
               xy=x+y1
               om2i(xy,z,1)=
     &              -alfa(x)*om2i(xy,z,1)-beta(z)*om2i(xy,z,3)
               om2r(xy,z,1)=
     &              alfa(x)*om2r(xy,z,1)+beta(z)*om2r(xy,z,3)
            end do
         end do
      end do

      if (nfzsym.eq.0) then
         do z=nz/2+1,nz
            zp=nzp-nz+z
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1
                  om2i(xy,zp,1)=
     &                 -alfa(x)*om2i(xy,zp,1)-beta(z)*om2i(xy,zp,3)
                  om2r(xy,zp,1)=
     &                 alfa(x)*om2r(xy,zp,1)+beta(z)*om2r(xy,zp,3)
               end do
            end do
         end do
      end if
c
c     Put the planes of the nonlinear term back onto the velocities
c     (here, the truncation to the normal grid occurs)
c     NOTE: THE ODDBALL MODES ARE STILL IN THERE!!!!
c
      if (nproc.eq.1) then
         call putxz(om2i(1,1,1),om2r(1,1,1),yb,4,ur,ui)
         call putxz(om2r(1,1,2),om2i(1,1,2),yb,5,ur,ui)
         call putxz(wr         ,wi         ,yb,8,ur,ui)
      else
#ifdef MPI
         call putpxz(om2i(1,1,1),om2r(1,1,1),yb,4,ur,ui,
     &        realg1,realg2,my_node)
         call putpxz(om2r(1,1,2),om2i(1,1,2),yb,5,ur,ui,
     &        realg1,realg2,my_node)
         call putpxz(wr         ,wi         ,yb,8,ur,ui,
     &        realg1,realg2,my_node)
#endif
      end if

      end subroutine nonlinp
