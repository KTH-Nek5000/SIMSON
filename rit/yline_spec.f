c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine yline_spec(graf,gxax,xmin,xmax,ymin,ymax,ipoint,
     &    x,z,jvar,eta,ur)
c
c     Cuts a line in the y-direction, possibly with subsampling
c
      implicit none

      include 'par.f'

      integer jvar,x,z,i,j,zp
      real graf(nyp,20),gxax(nyp,20)
      real xmin,xmax,ymin,ymax
      real eta(nyp)
      complex ur(memnx,memny,memnz,4)
      integer ipoint(3)
c
      write(*,*)'in ',z ,x,jvar

      ymin=1000000.0
      ymax=-100000.0

      zp=z
      if (z.lt.0) zp=nz+z

      ipoint(1)=nyp
      ipoint(2)=nyp
      ipoint(3)=nyp

      do i=1,3
         do j=1,ipoint(i)
            graf(j,i)=0.0
            gxax(j,i)=eta(nyp+1-j)
         end do
      end do

      do j=1,nyp
         graf(j,1)=real(ur(x+1,nyp+1-j,zp+1,jvar))
         graf(j,2)=aimag(ur(x+1,nyp+1-j,zp+1,jvar))
         graf(j,3)=abs(ur(x+1,nyp+1-j,zp+1,jvar))
         write(45,*)ur(x+1,nyp+1-j,zp+1,jvar)
      end do

      do i=1,3
         do j=1,nyp
            write(46,*)gxax(j,i),graf(j,i)
         end do
      end do


      xmin=eta(nyp)
      xmax=eta(1)
      do i=1,3
         do j=1,ipoint(i)
            ymin=min(ymin,graf(j,i))
            ymax=max(ymax,graf(j,i))
         end do
      end do
      ymin=ymin-.05*abs(ymin)
      ymax=ymax+.05*abs(ymax)

      open(unit=29,file='liney_spec.dat')
      do j=1,nyp
         write(29,2200)gxax(j,1),graf(j,1),graf(j,2)
      end do
 2200 format(f15.10,'  ',e20.12,e20.12)
      close(29)

      return

      end subroutine yline_spec
