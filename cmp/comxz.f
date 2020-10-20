c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine comxz(ur,ui,alfa,beta,prey)
c
c     Finds the streamwise and spanwise vorticity omx,omz
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real prey(nyp*2+15)
      real u3r(nx/2,mbz,nyp,3),u3i(nx/2,mbz,nyp,3)
      real om3r(nx/2,mbz,nyp,3),om3i(nx/2,mbz,nyp,3)
      real w3(nx/2,mbz,nyp)
      real alfa(nx/2*mbz),beta(nz)

      integer y,z,zb,i,x,mzb,nxz

      nxz=nx/2*mbz
      do zb=1,nzc,mbz
         do i=1,3
            call getxy(u3r(1,1,1,i),u3i(1,1,1,i),zb,i,ur,ui)
            call vchbf(u3r(1,1,1,i),w3,nyp,nxz,nxz,1,prey)
            call vchbf(u3i(1,1,1,i),w3,nyp,nxz,nxz,1,prey)
            do z=zb,zb+mbz-1
               mzb=z-zb+1
               do y=1,ny
                  do x=1,nx/2
                     u3r(x,mzb,y,i)=u3r(x,mzb,y,i)*(2./float(nyp-1))
                     u3i(x,mzb,y,i)=u3i(x,mzb,y,i)*(2./float(nyp-1))
                  end do
               end do
            end do
         end do
c
c     Construct vorticies
c
         do y=1,ny
            do z=zb,zb+mbz-1
               do x=1,nx/2
                  mzb=z-zb+1
                  om3r(x,mzb,y,2)=-u3i(x,mzb,y,1)*beta(z)
     &                            +u3i(x,mzb,y,3)*alfa(x)
                  om3i(x,mzb,y,2)= u3r(x,mzb,y,1)*beta(z)
     &                            -u3r(x,mzb,y,3)*alfa(x)
               end do
            end do
         end do
         call dcheb(om3r,u3r(1,1,1,3),ny,nxz,nxz)
         call dcheb(om3i,u3i(1,1,1,3),ny,nxz,nxz)
         call dcheb(om3r(1,1,1,3),u3r,ny,nxz,nxz)
         call dcheb(om3i(1,1,1,3),u3i,ny,nxz,nxz)
         do y=1,ny
            do z=zb,zb+mbz-1
               do x=1,nx/2
                  mzb=z-zb+1
                  om3r(x,mzb,y,1) = om3r(x,mzb,y,1)
     &                             +u3i(x,mzb,y,2)*beta(z)
                  om3i(x,mzb,y,1) = om3i(x,mzb,y,1)
     &                             -u3r(x,mzb,y,2)*beta(z)
                  om3r(x,mzb,y,3) =-om3r(x,mzb,y,3)
     &                             -u3i(x,mzb,y,2)*alfa(x)
                  om3i(x,mzb,y,3) =-om3i(x,mzb,y,3)
     &                             +u3r(x,mzb,y,2)*alfa(x)
               end do
            end do
         end do
c
c     Pad, transform and store
c
         if (ny+1.le.nyp) then
            do y=ny+1,max(nyp,ny+1)
               do z=zb,zb+mbz-1
                  do x=1,nx/2
                     mzb=z-zb+1
                     om3r(x,mzb,y,1)=0.0
                     om3i(x,mzb,y,1)=0.0
                     om3r(x,mzb,y,3)=0.0
                     om3i(x,mzb,y,3)=0.0
                  end do
               end do
            end do
         end if
         call vchbb(om3r,w3,nyp,nxz,nxz,1,prey)
         call vchbb(om3i,w3,nyp,nxz,nxz,1,prey)
         call putxy(om3r,om3i,zb,4,ur,ui)
         call vchbb(om3r(1,1,1,3),w3,nyp,nxz,nxz,1,prey)
         call vchbb(om3i(1,1,1,3),w3,nyp,nxz,nxz,1,prey)
         call putxy(om3r(1,1,1,3),om3i(1,1,1,3),zb,5,ur,ui)
      end do

      return

      end subroutine comxz
