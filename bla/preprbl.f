c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine preprbl(alfa,beta,eta,deta,xl,zl,prex,
     &     prey,prez,pres,prea,prexn,prezn,presn,prean,prezr,
     &     wint,wd1,gridx,gridy,gridz,dstar,deltaxyz2)
c
c     Initializes certain preprocessing information for FFTs, grid etc.
c
      implicit none

      include 'par.f'

      real alfa(nx/2,mbz),beta(nz),eta(nyp),deta(nyp),xl,zl
      real prex(nxp+15),prey(nyp*2+15),wint(nyp),wd1(nyp,4)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real prezr(nzp+15)
      real gridx(nx),gridy(nyp),gridz(nz)
      real dstar,dy,deltaxyz2(nyp)

      integer i,x,y,z
      real c(nyp)
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Geometrics
c
      do i=1,mbz
         do x=1,nx/2
            alfa(x,i)=2.*pi/xl*real(x-1)
         end do
      end do
      beta(1)=0.
      do z=2,nz/2+1
         beta(z)=2.*pi/zl*real(z-1)
         beta(nz+2-z)=-2.*pi/zl*real(z-1)
      end do
      do y=1,nyp
         eta(y)=cos(pi*real(y-1)/real(nyp-1))
      end do
c
c     Grid size (in physical coordinates)
c
      do i=1,nx
         gridx(i) = xl/real(nx)*real(i-1)/dstar
      end do
      do i=1,nz
         gridz(i) = zl/real(nz)*real(i-1)/dstar
      end do
      do i=1,nyp
         gridy(i)=(1.+cos(pi*real(i-1)/real(nyp-1)))/dstar
      end do
c
c     Grid spacing (in computational scaling, squared)
c
      do i=1,nyp
         if (i.eq.1) then
            dy = 1.-eta(2)
         else if (i.eq.nyp) then
            dy = 1.-eta(2)
         else
            dy = sqrt((eta(i+1)-eta(i))*(eta(i)-eta(i-1)))
         end if
         deltaxyz2(i) = ( xl/real(nx)*zl/real(nz)*dy )**(2./3.)
      end do
c      deta(1)=(1.-eta(2))/2.
      deta(1)=1.-eta(2)
      deta(nyp)=deta(1)
      do y=2,nyp-1
         deta(y)=(eta(y-1)-eta(y+1))*.5
      end do
      do y=1,nyp
         wint(y)=-.5
         do i=1,(nyp+1)/2-2
            wint(y) = wint(y) + cos(pi*real(2*i*(y-1))/real(nyp-1))/
     &           real((2*i)**2-1)
         end do
         wint(y) = wint(y) + .5*cos(pi*real(y-1))/real((nyp-1)**2-1)
         wint(y) = -4./real(nyp-1)*wint(y)
         if (y.eq.1.or.y.eq.nyp) wint(y)=wint(y)*.5
      end do
c
c     Calculate weight for derivative and second derivative at 1
c
      do y=1,nyp
         wd1(y,1)=0.0
         wd1(y,2)=0.0
         wd1(y,3)=0.0
         wd1(y,4)=0.0
         c(y)=1.
      end do
      c(1)=0.5
      c(nyp)=0.5
      prex(1)=1.
      do i=1,nyp
         prex(1)=-prex(1)
         do y=1,nyp
            wd1(y,1)=wd1(y,1)+real(i-1)**2*c(y)*c(i)*
     &           cos(pi*real((i-1)*(y-1))/real(nyp-1))*
     &           (2./real(nyp-1))
            wd1(y,2)=wd1(y,2)+real(i-1)**2*(real(i-1)**2-1.)/3.*
     &           c(y)*c(i)*
     &           cos(pi*real((i-1)*(y-1))/real(nyp-1))*
     &           (2./real(nyp-1))
c     dv and ddv at y=-1
            wd1(y,3)=wd1(y,3)+real(i-1)**2*c(y)*c(i)*prex(1)*
     &           cos(pi*real((i-1)*(y-1))/real(nyp-1))*
     &           (2./real(nyp-1))
            wd1(y,4)=wd1(y,4)-real(i-1)**2*(real(i-1)**2-1.)/3.*
     &           c(y)*c(i)*
     &           prex(1)*cos(pi*real((i-1)*(y-1))/real(nyp-1))*
     &           (2./real(nyp-1))
         end do
      end do
c
c     Initialize x-transform
c
      call vrffti(nxp,prex,0)
      call vrffti(nx,prexn,0)
c
c     Initialize y-transform
c
      call vcosti(nyp,prey,0)
c
c     Initialize z-transform
c
      if (nfzsym.eq.0) then
         call vcffti(nzp,prez,0)
         call vcffti(nz,prezn,0)
      else
         call vcosti(nzst,pres,0)
         call vsinti(nzat,prea,0)
         call vcosti(nz/2+1,presn,0)
         call vsinti(nz/2-1,prean,0)
      end if
      if (nzp.gt.1) call vrffti(nzp,prezr,0)

      end subroutine preprbl
