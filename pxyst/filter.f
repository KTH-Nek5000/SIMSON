c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine filter(plane,mx,my,ifilt,sfilt,xl)
c
c     Scale variables
c
      implicit none

      include 'par.f'

      integer mx,my,ifilt
      real plane(mx,my),sfilt
c
      real w(nx+2,nyp),w2(nx+2,nyp)
      real alfa(nx/2+1),pi,xl
      real prex(nx+15)
      integer x,y
      parameter (pi = 3.1415926535897932385)
c
c     initialize transforms and wave numbers
c
      write(*,*) 'in filter doing ',ifilt
      call vrffti(nx,prex,0)
      do x=1,nx/2+1
         alfa(x)=2.*pi/xl*real(x-1)
      end do
c
c     Copy the dependent variable  into a work array
c
      do y=1,my
         do x=1,nx
            w(x,y)=plane(x,y)
         end do
      end do
c
c     Transform data to fourier space
c
      write(*,*) 'Calling vrfftf'
      call vrfftf(w,w(2,1),w2,w2(2,1),nx,my,2,nx+2,prex)
c
c     Low pass filter in fourier space and normalize
c
      do y=1,my
         do x=1,nx
            w(x,y)=w(x,y)*
     &           exp(-(alfa((x+1)/2)*sfilt/2.*pi)**2)*(1./real(nx))
         end do
      end do
c
c     Reset odd-ball
c
      do y=1,my
         w(nx+1,y)=0.0
      end do
c
c     Transform the filtered data back to physical space
c
      call vrfftb(w,w(2,1),w2,w2(2,1),nx,my,2,nx+2,prex)
c
c     Copy the result back to the dependent variable
c
      do y=1,my
         do x=1,nx
            plane(x,y)=w(x,y)
         end do
      end do
c
c     Wrap periodic data
c
      do y=1,my
        plane(nx+1,y)=plane(1,y)
      end do

      return

      end subroutine filter
