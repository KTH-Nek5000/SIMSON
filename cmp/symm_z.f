c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine symm_z(ur,ui,yb,u2r,u2i,symflag)
c
c     Flips the velocity field about the z=0 plane, adds the flipped
c     field to the original field and divides by two. The product is
c     a field in wave space symmetric about the z=0 plane in the u and v
c     components, and anti-symmetric in the w component.
c
      implicit none

      include 'par.f'

      integer yb
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real u22r((nxp/2+1)*mby,nzd,3),u22i((nxp/2+1)*mby,nzd,3)
      integer i,xy,z,y,x,symflag

      if (symflag.ne.1.and.symflag.ne.2) then
         write(*,*) 'symflag needs to be 1 or 2!'
         stop
      end if

      do i=1,3
        call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui)
      end do
c
c     flip in wave space
c     note that we are not modifying the kz=0 wavenumber
c
      do i=1,2
         do z=2,nzd
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  if (symflag.eq.1) then
c
c     force the field to be symmetric
c
                     xy=x+(y-yb)*(nxp/2+1)
                     u22r(xy,z,i)=u2r(xy,nzd-z+2,i)
                     u22i(xy,z,i)=u2i(xy,nzd-z+2,i)
                  else
c
c     force the field to be anti-symmetric
c
                     xy=x+(y-yb)*(nxp/2+1)
                     u22r(xy,z,i)=-u2r(xy,nzd-z+2,i)
                     u22i(xy,z,i)=-u2i(xy,nzd-z+2,i)
                  end if
               end do
            end do
         end do
      end do
c
c     the z-component of velocity needs to have different symmetry
c
      do z=2,nzd
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               if (symflag.eq.1) then
c
c     force the field to be anti-symmetric
c
                  xy=x+(y-yb)*(nxp/2+1)
                  u22r(xy,z,3)=-u2r(xy,nzd-z+2,3)
                  u22i(xy,z,3)=-u2i(xy,nzd-z+2,3)
               else
c
c     force the field to be symmetric
c
                  xy=x+(y-yb)*(nxp/2+1)
                  u22r(xy,z,3)=u2r(xy,nzd-z+2,3)
                  u22i(xy,z,3)=u2i(xy,nzd-z+2,3)
               end if
            end do
         end do
      end do
c
c     Add also kz=0 component (note changes sign for z-component)
c
      if (symflag.eq.1) then
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               xy=x+(y-yb)*(nxp/2+1)
               u22r(xy,1,1)= u2r(xy,1,1)
               u22i(xy,1,1)= u2i(xy,1,1)
               u22r(xy,1,2)= u2r(xy,1,2)
               u22i(xy,1,2)= u2i(xy,1,2)
               u22r(xy,1,3)=-u2r(xy,1,3)
               u22i(xy,1,3)=-u2i(xy,1,3)
            end do
         end do
      else
         do y=yb,min(nyp,yb+mby-1)
            do x=1,nx/2
               xy=x+(y-yb)*(nxp/2+1)
               u22r(xy,1,1)=-u2r(xy,1,1)
               u22i(xy,1,1)=-u2i(xy,1,1)
               u22r(xy,1,2)=-u2r(xy,1,2)
               u22i(xy,1,2)=-u2i(xy,1,2)
               u22r(xy,1,3)= u2r(xy,1,3)
               u22i(xy,1,3)= u2i(xy,1,3)
            end do
         end do
      end if


c
c     add flipped solution and divide sum by 2
c
      do i=1,3
         do z=1,nzd
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  xy=x+(y-yb)*(nxp/2+1)
                  u2r(xy,z,i)=(u2r(xy,z,i)+u22r(xy,z,i))/2.0
                  u2i(xy,z,i)=(u2i(xy,z,i)+u22i(xy,z,i))/2.0
               end do
            end do
         end do
      end do

      do i=1,3
        call putxz(u2r(1,1,i),u2i(1,1,i),yb,i,ur,ui)
      end do

      return

      end subroutine symm_z
