c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine mgrid(plane,pxys,pxys2,pxysth,pxysth2,mxys2,mxysth2,
     &     ivar,ivarth,mx,my,grid,xr,xl,mxr,eta)
c
c     Selects variable ivar or ivarth, shifts data
c     and makes a grid
c
      implicit none

      include 'par.f'

      integer mxys2,mxysth2
      real pxys(nx,nyp,nxys),pxys2(nx,nyp,mxys2)
      real pxysth(nx,nyp,nxysth),pxysth2(nx,nyp,mxysth2)
      integer mx,my,mxr,ivar,ivarth
      real eta(nyp)
      real plane(mx,my),grid(mx,my,2),xr,xl
c
      integer y,x,xpr
      if (ivar.gt.0) then
         do y=1,nyp
          do x=1,nx
             xpr=mod(10*nx+x+mxr-1,nx)+1
             plane(x,y)=pxys(xpr,y,ivar)
          end do
       end do
      else
         if (ivar.lt.0) then
            do y=1,nyp
               do x=1,nx
                  xpr=mod(10*nx+x+mxr-1,nx)+1
                  plane(x,y)=pxys2(xpr,y,-ivar)
               end do
            end do
         else
            if (ivarth.gt.0) then
               do y=1,nyp
                  do x=1,nx
                     xpr=mod(10*nx+x+mxr-1,nx)+1
                     plane(x,y)=pxysth(xpr,y,ivarth)
                  end do
               end do
            else
               do y=1,nyp
                  do x=1,nx
                     xpr=mod(10*nx+x+mxr-1,nx)+1
                     plane(x,y)=pxysth2(xpr,y,-ivarth)
                  end do
               end do
            end if
         end if
      end if
c
c     Extend periodically in the yx-direction
c
      do y=1,nyp
        plane(nx+1,y)=plane(1,y)
      end do
c
c     Make grid
c
      do y=1,my
         do x=1,nx+1
            grid(x,y,1)=real(x-nx/2-1)/real(nx)*xl+xr
            grid(x,y,2)=eta(nyp+1-y)
         end do
      end do

      return

      end subroutine mgrid
