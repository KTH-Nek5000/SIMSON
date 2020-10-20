c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine ddx(da,a,xl)
c
c     Performs da := ddx (a)
c     data are periodic with period xl
c
      implicit none

      include 'par.f'

      real a(nx,nyp),da(nx,nyp),xl

      real w(nx+2,nyp),w2(nx+2,nyp)
      real prex(nx+15)
      integer x,y
      real alfa(nx/2),pi,h
      parameter (pi = 3.1415926535897932385)
c
c     Initialize transform and wave numbers
c
      call vrffti(nx,prex,0)
      do x=1,nx/2
         alfa(x)=2.*pi/xl*real(x-1)
      end do
c
c     Copy a into work array and normalize for the transform
c
      do y=1,nyp
         do x=1,nx
            w(x,y)=a(x,y)*(1./real(nx))
         end do
      end do
c
c     Fourier transform in x
c
      call vrfftf(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Calculate dwdx as i*alfa*w
c
      do y=1,nyp
         do x=1,nx/2
c
c     Imaginary part
c
            h=w(2*x-1,y)*alfa(x)
c
c     Real part
c
            w(2*x-1,y)=-w(2*x,y)*alfa(x)
            w(2*x,y)=h
         end do
      end do
c
c     Reset odd-ball
c
      do y=1,nyp
         w(nx+1,y)=0.0
      end do
c
c     Transform back to physical space
c
      call vrfftb(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Copy the result to the output array
c
      do y=1,nyp
         do x=1,nx
            da(x,y)=w(x,y)
         end do
      end do

      end subroutine ddx
