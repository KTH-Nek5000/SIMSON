c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine reschk(graf,gxax,xmin,xmax,ymin,ymax,ipoint,
     &  prey,jvar,pxy,w3,ur)

      implicit none

      include 'par.f'

      complex pxy(nx/2,nyp),w3(nx/2,nyp)
      real graf(nx+nyp+nz,3)
      real gxax(nx+nyp+nz,3)
      real prey(nyp*2+15)
      real xmin,xmax,ymin,ymax
      integer ipoint(3)
      integer jvar,ml,zf
      real ur(memnx,memny,memnz,4)

      integer x,y,z,i,j
      real c

      ipoint(1)=nx/2-1
      ipoint(2)=ny-1
      ipoint(3)=nz/2-1
      do i=1,3
         do j=1,ipoint(i)
            graf(j,i)=0.0
            gxax(j,i)=log10(real(j))
         end do
      end do
      do z=1,nzc
         c=1.
         if (z.gt.1.and.nfzsym.eq.1) c=2.
         call getxyp(pxy,z,jvar,ur)
         call vchbf(pxy,w3,nyp,nx,nx,1,prey)
         do y=1,ny
            do x=1,nx/2
               pxy(x,y)=pxy(x,y)*(2./real(nyp-1))
            end do
         end do
c
c     x-resolution
c
         do y=1,ny
            do x=2,nx/2
               graf(x-1,1)=graf(x-1,1)+pxy(x,y)*conjg(pxy(x,y))*c
            end do
         end do
c
c     z-resolution
c
         if (z.gt.1) then
            zf=z
            if (z.gt.nz/2) zf=nz+2-z
            do y=1,ny
               do x=1,nx/2
                  graf(zf-1,3)=graf(zf-1,3)+pxy(x,y)*conjg(pxy(x,y))*c
               end do
            end do
         end if
c
c     y-resolution note that over ny the modes are zero
c
         do y=2,ny
            do x=1,nx/2
               graf(y-1,2)=graf(y-1,2)+pxy(x,y)*conjg(pxy(x,y))*c
            end do
         end do
      end do
      ymin=100.
      ymax=-100.
      do i=1,3
         do j=1,ipoint(i)
            if (graf(j,i).gt.1E-50) then
               graf(j,i)=log10(graf(j,i))
            else
               graf(j,i)=-50.
            end if
            ymin=min(ymin,graf(j,i))
            ymax=max(ymax,graf(j,i))
         end do
      end do
      ymin=ymin-.5
      ymax=ymax+.5
      xmin=0.
      ml=max(nx/2-1,ny-1)
      ml=max(ml,nz/2-1)
      xmax=log10(real(ml))

      return

      end subroutine reschk
