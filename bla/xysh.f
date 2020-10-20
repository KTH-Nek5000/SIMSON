c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine xysh(boxr,boxi,xs,zs,alfa,beta,zb)
c
c     Shift a Fourier-transformed box so that xs=0 and zs=0
c
      implicit none

      include 'par.f'

      integer zb
      real boxr(nx/2*mbz,nyp),boxi(nx/2*mbz,nyp)
      real alfa(nx/2*mbz),beta(nz)
      real xs,zs

      integer xz,y,x,z
      real hr,arg,bbeta(nx/2*mbz)
c
c     Compute the full beta plane for vectorisation
c
      do z=1,mbz
         do x=1,nx/2
            xz=x+nx/2*(z-1)
            bbeta(xz)=beta(z+zb-1)
         end do
      end do
c
c     Perform the shifting using
c     arg  = - xs*alfa - zs*beta
c     box  = box*exp(i*arg)
c          = box*(cos(arg)+i*sin(arg))
c     boxr = boxr*cos(arg) - boxi*sin(arg)
c     boxi = boxi*cos(arg) + boxr*sin(arg)
c
      do y=1,nyp
         do xz=1,nx/2*mbz
            arg = - xs*alfa(xz) - zs*bbeta(xz)
            hr         = boxr(xz,y)*cos(arg) - boxi(xz,y)*sin(arg)
            boxi(xz,y) = boxi(xz,y)*cos(arg) + boxr(xz,y)*sin(arg)
            boxr(xz,y) = hr
         end do
      end do

      end subroutine xysh
