c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine infilt(totxys,totxysth,sfilt,xl,scalarind)
c
c     Filter input data
c
      implicit none

      include 'par.f'

      real totxys(nx,nyp,nxys),totxysth(nx,nyp,nxysth,scalar)
      real xl,sfilt
c
      real w(nx+2,nyp),w2(nx+2,nyp)
      real alfa(nx/2+1),pi
      real prex(nx+15)
      integer x,y,i,scalarind
      parameter (pi = 3.1415926535897932385)
c
c     Initialize transforms and wave numbers
c
      call vrffti(nx,prex,0)
      do x=1,nx/2+1
         alfa(x)=2.*pi/xl*real(x-1)
      end do
      do i=1,nxys
c
c     Copy the dependent variable  into a work array
c
         do y=1,nyp
            do x=1,nx
               w(x,y)=totxys(x,y,i)
            end do
         end do
c
c     Transform data to fourier space
c
        call vrfftf(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Low pass filter in fourier space and normalize
c
        do y=1,nyp
           do x=1,nx+1
              w(x,y)=w(x,y)*
     &             exp(-(alfa((x+1)/2)*sfilt/2.*pi)**2)*(1./real(nx))
           end do
        end do
c
c     Transform the filtered data back to physical space
c
        call vrfftb(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Copy the result back to the dependent variable
c
        do y=1,nyp
           do x=1,nx
              totxys(x,y,i)=w(x,y)
           end do
        end do
      end do

      do i=1,nxysth
c
c     Copy the dependent variable  into a work array
c
         do y=1,nyp
            do x=1,nx
               w(x,y)=totxysth(x,y,i,scalarind)
            end do
         end do
c
c     Transform data to fourier space
c
         call vrfftf(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Low pass filter in fourier space and normalize
c
         do y=1,nyp
            do x=1,nx+1
               w(x,y)=w(x,y)*
     &              exp(-(alfa((x+1)/2)*sfilt/2.*pi)**2)*(1./real(nx))
            end do
         end do
c
c     Transform the filtered data back to physical space
c
         call vrfftb(w,w(2,1),w2,w2(2,1),nx,nyp,2,nx+2,prex)
c
c     Copy the result back to the dependent variable
c
         do y=1,nyp
            do x=1,nx
               totxysth(x,y,i,scalarind)=w(x,y)
            end do
         end do
      end do

      return

      end subroutine infilt
