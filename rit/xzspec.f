c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine xzspec(plane,mkx,mkz,grid,jvar,logax,logc,sym,wint,
     &xl,zl,pxz,ur)
c
c     Integrates in y direction and calculates a grid in the xz spectral plane
c
      implicit none

      include 'par.f'

      integer mkx,mkz,jvar
      logical logax,logc
      real plane(mkx,mkz),grid(mkx,mkz,2),wint(nyp),xl,zl
      complex pxz(nx/2+1,nz)
      complex ur(memnx,memny,memnz,4)
      real sym
c
      integer x,y,z
      real pi
      parameter (pi = 3.141592653589789)
c
      do x=1,mkx
         do z=1,mkz
            plane(x,z)=0.
         end do
      end do

      if (.not.logax) then
c
c     Integration loop
c
         do y=1,nyp
            call getxzp(pxz,y,jvar,ur,sym)
            do z=2,nz/2
               do x=1,mkx
                  plane(x,z-1)=plane(x,z-1)+wint(y)*
     &                 pxz(x,z+nz/2)*conjg(pxz(x,z+nz/2))
                  plane(x,z+nz/2-1)=plane(x,z+nz/2-1)+wint(y)*
     &                 pxz(x,z)*conjg(pxz(x,z))
               end do
            end do
            do x=1,mkx
               plane(x,nz/2)=plane(x,nz/2)+wint(y)*
     &              (pxz(x,1)*conjg(pxz(x,1)))
            end do
         end do
         do z=1,mkz
            do x=1,mkx
               grid(x,z,1)=real(x-1)*2.*pi/xl
               grid(x,z,2)=real(z-nz/2)*2.*pi/zl
            end do
         end do
      else
c
c     Integration loop
c
         do y=1,nyp
            call getxzp(pxz,y,jvar,ur,sym)
            do z=1,mkz
               do x=1,mkx
                  plane(x,z)=plane(x,z)+wint(y)*
     &                 (pxz(x+1,z+1)*conjg(pxz(x+1,z+1))+
     &                 pxz(x+1,nz+1-z)*conjg(pxz(x+1,nz+1-z)))
               end do
            end do
         end do
         do z=1,mkz
            do x=1,mkx
               grid(x,z,1)=log10(real(x)*2.*pi/xl)
               grid(x,z,2)=log10(real(z)*2.*pi/zl)
            end do
         end do
      end if
      if (logc) then
         do z=1,mkz
            do x=1,mkx
               if (plane(x,z).eq.0.) then
                  plane(x,z)=30.
               else
                  plane(x,z)=-log10(abs(plane(x,z)))
               end if
            end do
         end do
      end if

      return

      end subroutine xzspec
