c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine gauss(ur,jvar,sym,ifilt,sfilt,zfilt,
     &     alfa,beta,prex,prez,
     &           pxy,pxz,cpxz,pxz2,w)
c
      implicit none

      include 'par.f'

      complex ur(memnx,memny,memnz,4)
      complex pxy(nx/2,nyp)
      real pxz(nx+2,nz),w(nx+2,nz),pxz2(nx/2,nz)
      complex cpxz(nx/2+1,nz)
      real prex(nx+15),prez(nz*2+15),alfa(nx/2),beta(nz)
      real zfilt,sfilt,sym
      integer jvar,ifilt
c
      integer x,y,z
      real pi
      parameter (pi = 3.1415926535897932385)

      if (ifilt.eq.1.or.ifilt.eq.4) then
c
c     Gaussian lowpass filter
c
         do z=1,nzc
            call getxyp(pxy,z,jvar,ur)
            if (ifilt.eq.1) then
               do x=1,nx/2
                  do y=1,nyp
                     pxy(x,y)=pxy(x,y)*
     &                    exp(-(alfa(x)**2+beta(z)**2)*(sfilt/2./pi)**2)
                  end do
               end do
            else
               do x=1,nx/2
                  do y=1,nyp
                     pxy(x,y)=pxy(x,y)*
     &                    exp(-alfa(x)**2*(sfilt/2./pi)**2-
     &                    beta(z)**2*(zfilt/2./pi)**2)
                  end do
               end do
            end if
            call putxyp(pxy,z,4,ur)
         end do
      end if
      if (ifilt.eq.2.or.ifilt.eq.5) then
c
c     Gaussian highpass filter
c
         do z=1,nzc
            call getxyp(pxy,z,jvar,ur)
            if (ifilt.eq.2) then
               do x=1,nx/2
                  do y=1,nyp
                     pxy(x,y)=pxy(x,y)*
     &                    (1.-exp(-(alfa(x)**2+beta(z)**2)*
     &                    (sfilt/2./pi)**2))
                  end do
               end do
            else
               do x=1,nx/2
                  do y=1,nyp
                     pxy(x,y)=pxy(x,y)*
     &                    (1.-exp(-alfa(x)**2*(sfilt/2./pi)**2-
     &                    beta(z)**2*(zfilt/2./pi)**2))
                  end do
               end do
            end if
            call putxyp(pxy,z,4,ur)
         end do
      end if
      if (ifilt.eq.3.or.ifilt.eq.6) then
c
c     rms calculation through gaussian highpass filtering,squaring,
c     lowpass filtering and taking the square root
c     start by precalculating filter
c
         if (ifilt.eq.3) then
            do z=1,nz
               do x=1,nx/2
                  pxz2(x,z)=exp(-(alfa(x)**2+beta(z)**2)*
     &                 (sfilt/2./pi)**2)
               end do
            end do
         else
            do z=1,nz
               do x=1,nx/2
                  pxz2(x,z)=
     &                 exp(-alfa(x)**2*(sfilt/2./pi)**2-
     &                 beta(z)**2*(zfilt/2./pi)**2)
               end do
            end do
         end if
         do y=1,nyp
            call getxzp(pxz,y,jvar,ur,sym)
c
c     Highpass filter
c
            do z=1,nz
               do x=1,nx/2
                  cpxz(x,z)=cpxz(x,z)*(1.-pxz2(x,z))
               end do
               cpxz(nx/2+1,z)=0.0
            end do
            call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
            do z=1,nz
               do x=1,nx
                  pxz(x,z)=pxz(x,z)**2
               end do
            end do
            call vrfftf(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
            call vcfftf(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
c
c     Lowpass filter and scale for transform
c
            do z=1,nz
               do x=1,nx/2
                  cpxz(x,z)=cpxz(x,z)*(1./real(nx*nz))*pxz2(x,z)
               end do
               cpxz(nx/2+1,z)=0.0
            end do
            call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
            do z=1,nz
               do x=1,nx
                  pxz(x,z)=sqrt(max(pxz(x,z),0.))
               end do
            end do
            call vrfftf(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
            call vcfftf(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
            do z=1,nz
               do x=1,nx
                  pxz(x,z)=pxz(x,z)*(1./real(nx*nz))
               end do
            end do
            call putxzp(pxz,y,4,ur)
         end do
      end if

      return

      end subroutine gauss
