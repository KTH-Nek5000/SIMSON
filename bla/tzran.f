c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine tzran(ttz,nzt,seed,prezr)
c
c     Create a function ttz with random dependence of z
c     containing nzt Fourier components
c
      implicit none

      include 'par.f'

      integer seed
      real ttz(nzp+2),nzt
c
c     Local variables
c
      integer z
      real w(nzp+2)
      complex ph
      real prezr(nzp+15)
c
c     Local parameters
c
      complex im
      parameter (im=(0.,1.))
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     External function
c
      real ran2
c
c     Set function to zero
c
      do z=1,nzp
         ttz(z)=0.0
      end do
c
c     Construct ttz with alternating real/imaginary parts
c     up to nzt/2+1 and zeros at positions 2 and nzt/2+2
c
c     For beta=0 (mean) the amplitude is randomized
c     maximum amplitude is sqrt(2)/2
c
      ph = exp(im*2.*pi*ran2(seed))
      ttz(1) = real(ph)*sqrt(2.)
      ttz(2) = 0.0

      if (nfzsym.eq.0) then
c
c     For the non-symmetric case, we take a constant
c     amplitude and randomize the phase
c
         do z=2,nzt/2+1
            ph = exp(im*2.*pi*ran2(seed))
            ttz(2*z-1) = real(ph)
            ttz(2*z) = aimag(ph)
         end do
      else
c
c     For the symmetric case the transform should be real
c     which leads to random amplitude instead of random phase
c     we raise the amplitude to compensate
c
         do z=2,nzt/2+1
            ph = exp(im*2.*pi*ran2(seed))
            ttz(2*z-1) = real(ph)*sqrt(2.)
            ttz(2*z) = 0.0
         end do
      end if
c
c     Transform to real space
c
      if (nzp.gt.1) call vrfftb(ttz,ttz(2),w,w(2),nzp,1,2,1,prezr)

      end subroutine tzran
