c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine pvar(pxys,pxys2,pxysth,pxysth2,nxys2,nxysth2,
     &     totxys,totxysth,xl,zl,eta,re,dstar,pr,scalarind,mhd_n,
     &     firstmax)
c
c     Processes raw data to produce central moments etc
c
      implicit none

      include 'par.f'

      integer nxys2,nxysth2,scalarind
      real totxys(nx,nyp,nxys)
      real pxys(nx,nyp,nxys),pxys2(nx,nyp,nxys2),wxys(nx,nyp,nxys)
      real xl,zl,eta(nyp),re,dstar
      real pr(scalar)
      real totxysth(nx,nyp,nxysth,scalar)
      real pxysth(nx,nyp,nxysth),pxysth2(nx,nyp,nxysth2)
      integer i,j,ivar,kvar,lvar,x,y
      real dx,dz,dy,h1,ym,arg,x1
      real prex(nx+15),prey(nyp*2+15)
      real wxy(nx,nyp),wxy2(nx+2,nyp)
      real wxy3(nx,nyp)
      real urmssh(nx+2,nyp)
      real alfa(nx/2+1),pi
      parameter (pi = 3.1415926535897932385)
      real mhd_n
c
c     Blasius
c
      integer n,ip,im,firstmax
      real et,etap,etam,deta,xx,yy
      parameter (n=1000)
      real etab(0:n),f(0:n),fd(0:n),fdd(0:n)
      real etamax,rmsmin,rmsmax,dummy
c
c     Initialize transforms and wave numbers
c
      call vrffti(nx,prex,0)
      call vcosti(nyp,prey,0)
      do x=1,nx/2+1
         alfa(x)=2.*pi/xl*real(x-1)
      end do
c
c     First turn data upside down
c
      do i=1,nxys
c     this line was replace because the equivalence statement was removed
c     do y=1,(nyp-1)/2
         do y=1,nyp
            do x=1,nx
               dummy = totxys(x,nyp+1-y,i)
               pxys(x,nyp+1-y,i) = totxys(x,y,i)
               pxys(x,y,i) = dummy
            end do
         end do
      end do

      if (scalar.gt.0) then
         do i=1,nxysth
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,i)=totxysth(x,nyp+1-y,i,scalarind)
               end do
            end do
         end do
      end if
c
c     Calculate derivatives of mean velocities
c
      do i=1,3
         call ddx(pxys2(1,1,13+i),pxys(1,1,i),xl)
         call ddy(pxys2(1,1,16+i),pxys(1,1,i),2./dstar)
         call ddx(pxys2(1,1,60+i),pxys2(1,1,13+i),xl)
         call ddy(pxys2(1,1,63+i),pxys2(1,1,16+i),2./dstar)
      end do

c
c    Calaulate derivatives of mean scalar
c
      if (scalar.gt.0) then
         call ddx(pxysth2(1,1,1),pxysth(1,1,1),xl)
         call ddy(pxysth2(1,1,2),pxysth(1,1,1),2./dstar)
      end if
c
c     Calaulate derivatives of mean pressure
c
      if (nxys.ge.31) then
         call ddx(pxys2(1,1,115),pxys(1,1,31),xl)
         call ddy(pxys2(1,1,116),pxys(1,1,31),2./dstar)
      end if
c
c     Total convection for uiui
c
      if (nxys.ge.71) then
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,63)+pxys(x,y,68)+pxys(x,y,70)
               wxys(x,y,2)=pxys(x,y,64)+pxys(x,y,66)+pxys(x,y,71)
            end do
         end do
         call ddx(wxys(1,1,3),wxys(1,1,1),xl)
         call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,125)=(wxys(x,y,3)+wxys(x,y,4))/2.
            end do
         end do
      end if
c
c     Total convection for scalar squared
c
      if (scalar.gt.0.and.nxysth.ge.11) then
         call ddx(wxys(1,1,1),pxysth(1,1,6),xl)
         call ddy(wxys(1,1,2),pxysth(1,1,7),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,26)=(wxys(x,y,1)+wxys(x,y,2))/2.
            end do
         end do
      end if
c
c     Total convection for uiuj (142-147)
c
      if (nxys.ge.72) then
         call ddx(wxys(1,1,1),pxys(1,1,63),xl)
         call ddy(wxys(1,1,2),pxys(1,1,66),2./dstar)
         call ddx(wxys(1,1,3),pxys(1,1,68),xl)
         call ddy(wxys(1,1,4),pxys(1,1,69),2./dstar)
         call ddx(wxys(1,1,5),pxys(1,1,70),xl)
         call ddy(wxys(1,1,6),pxys(1,1,71),2./dstar)
         call ddx(wxys(1,1,7),pxys(1,1,66),xl)
         call ddy(wxys(1,1,8),pxys(1,1,68),2./dstar)
         call ddx(wxys(1,1,9),pxys(1,1,67),xl)
         call ddy(wxys(1,1,10),pxys(1,1,72),2./dstar)
         call ddx(wxys(1,1,11),pxys(1,1,72),xl)
         call ddy(wxys(1,1,12),pxys(1,1,69),2./dstar)
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,141+i)=wxys(x,y,2*i-1)+wxys(x,y,2*i)
               end do
            end do
         end do
      end if
c
c     Total convection for ui-theta (38-40)
c
      if (scalar.gt.0.and.nxysth.ge.33) then
         call ddx(wxys(1,1,1),pxysth(1,1,16),xl)
         call ddx(wxys(1,1,2),pxysth(1,1,19),xl)
         call ddx(wxys(1,1,3),pxysth(1,1,20),xl)
         call ddy(wxys(1,1,4),pxysth(1,1,19),2./dstar)
         call ddy(wxys(1,1,5),pxysth(1,1,17),2./dstar)
         call ddy(wxys(1,1,6),pxysth(1,1,21),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,38)=wxys(x,y,1)+wxys(x,y,4)
               pxysth2(x,y,39)=wxys(x,y,2)+wxys(x,y,5)
               pxysth2(x,y,40)=wxys(x,y,3)+wxys(x,y,6)
            end do
         end do
      end if
c
c     Take away mean component for 3rd and 4th order velocity central momentum
c
      if (nxys.ge.77) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys(x,y,74+i)=pxys(x,y,74+i)-4*pxys(x,y,62+i)
     &                 *pxys(x,y,i)+6*pxys(x,y,3+i)*pxys(x,y,i)**2
     &                 -3*pxys(x,y,i)**4
               end do
            end do
         end do
      end if
      if (nxys.ge.65) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys(x,y,62+i)=pxys(x,y,62+i)+2*pxys(x,y,i)**3
     &                 -3*pxys(x,y,3+i)*pxys(x,y,i)
               end do
            end do
         end do
      end if
c
c     Take away mean component for other triple velocity correlations
c
      if (nxys.ge.72) then
         do y=1,nyp
            do x=1,nx
c     u^2 v
               pxys(x,y,66)=pxys(x,y,66)+2*pxys(x,y,1)**2*pxys(x,y,2)
     &              -2*pxys(x,y,13)*pxys(x,y,1)-pxys(x,y,4)*pxys(x,y,2)
c     u^2 w
               pxys(x,y,67)=pxys(x,y,67)+2*pxys(x,y,1)**2*pxys(x,y,3)
     &           -2*pxys(x,y,14)*pxys(x,y,1)-pxys(x,y,4)*pxys(x,y,3)
c     v^2 u
               pxys(x,y,68)=pxys(x,y,68)+2*pxys(x,y,2)**2*pxys(x,y,1)
     &              -2*pxys(x,y,13)*pxys(x,y,2)-pxys(x,y,5)*pxys(x,y,1)
c     v^2 w
               pxys(x,y,69)=pxys(x,y,69)+2*pxys(x,y,2)**2*pxys(x,y,3)
     &              -2*pxys(x,y,15)*pxys(x,y,2)-pxys(x,y,5)*pxys(x,y,3)
c     w^2 u
               pxys(x,y,70)=pxys(x,y,70)+2*pxys(x,y,3)**2*pxys(x,y,1)
     &              -2*pxys(x,y,14)*pxys(x,y,3)-pxys(x,y,6)*pxys(x,y,1)
c     w^2 v
               pxys(x,y,71)=pxys(x,y,71)+2*pxys(x,y,3)**2*pxys(x,y,2)
     &              -2*pxys(x,y,15)*pxys(x,y,3)-pxys(x,y,6)*pxys(x,y,2)
c     u v w
               pxys(x,y,72)=pxys(x,y,72)-pxys(x,y,13)*pxys(x,y,3)
     &              -pxys(x,y,14)*pxys(x,y,2)-pxys(x,y,15)*pxys(x,y,1)
     &              +2*pxys(x,y,3)*pxys(x,y,2)*pxys(x,y,1)
            end do
         end do
      end if
c
c     Take away mean component for 3rd and 4th order scalar central momentum
c
      if (scalar.gt.0.and.nxysth.ge.35) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,35)=pxysth(x,y,35)-4*pxysth(x,y,34)
     &              *pxysth(x,y,1)+6*pxysth(x,y,2)*pxysth(x,y,1)**2
     &              -3*pxysth(x,y,1)**4
            end do
         end do
      end if
      if (scalar.gt.0.and.nxysth.ge.34) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,34)=pxysth(x,y,34)+2*pxysth(x,y,1)**3
     &              -3*pxysth(x,y,2)*pxysth(x,y,1)
            end do
         end do
      end if
c
c     Take away mean component for velocity^2-scalar
c
      if (scalar.gt.0.and.nxysth.ge.21) then
         do ivar=16,18
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=pxysth(x,y,ivar)
     &                 -pxysth(x,y,1)*pxys(x,y,ivar-12)
     &                 -2*pxys(x,y,ivar-15)*pxysth(x,y,ivar-13)
     &                 +2*pxys(x,y,ivar-15)**2*pxysth(x,y,1)
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,19)=pxysth(x,y,19)-pxysth(x,y,1)*pxys(x,y,13)
     &              -pxysth(x,y,3)*pxys(x,y,2)-pxysth(x,y,4)*pxys(x,y,1)
     &              +2*pxys(x,y,1)*pxys(x,y,2)*pxysth(x,y,1)
               pxysth(x,y,20)=pxysth(x,y,20)-pxysth(x,y,1)*pxys(x,y,14)
     &              -pxysth(x,y,3)*pxys(x,y,3)-pxysth(x,y,5)*pxys(x,y,1)
     &              +2*pxys(x,y,1)*pxys(x,y,3)*pxysth(x,y,1)
               pxysth(x,y,21)=pxysth(x,y,21)-pxysth(x,y,1)*pxys(x,y,15)
     &              -pxysth(x,y,4)*pxys(x,y,3)-pxysth(x,y,5)*pxys(x,y,2)
     &              +2*pxys(x,y,2)*pxys(x,y,3)*pxysth(x,y,1)
            end do
         end do
      end if
c
c     Take away mean component for velocity-scalar^2
c
      if (scalar.gt.0.and.nxysth.ge.8) then
         do ivar=6,8
            do y=1,nyp
              do x=1,nx
                  pxysth(x,y,ivar)=pxysth(x,y,ivar)
     &                 -pxysth(x,y,2)*pxys(x,y,ivar-5)
     &                 -2*pxysth(x,y,1)*pxysth(x,y,ivar-3)
     &                 +2*pxys(x,y,ivar-5)*pxysth(x,y,1)**2
               end do
            end do
         end do
      end if
c
c     Take away mean component for 3rd and 4th order pressure central momentum
c
      if (nxys.ge.79) then
         do y=1,nyp
            do x=1,nx
               pxys(x,y,79)=pxys(x,y,79)-4*pxys(x,y,78)
     &              *pxys(x,y,31)+6*pxys(x,y,32)*pxys(x,y,31)**2
     &              -3*pxys(x,y,31)**4
            end do
         end do
      end if
      if (nxys.ge.78) then
         do y=1,nyp
            do x=1,nx
               pxys(x,y,78)=pxys(x,y,78)+2*pxys(x,y,31)**3
     &              -3*pxys(x,y,32)*pxys(x,y,31)
            end do
         end do
      end if
c
c     Take away mean component and take square root for rms
c
      do j=3,9,6
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys(x,y,i+j)=sqrt(max(0.,pxys(x,y,i+j)-
     &                 pxys(x,y,i+j-3)**2))
               end do
            end do
         end do
      end do
c
c     Take away mean component and take square root for rms of scalar field
c
      if (scalar.gt.0.and.nxysth.ge.2) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,2)=sqrt(max(0.,pxysth(x,y,2)
     &              -pxysth(x,y,1)**2))
            end do
         end do
      end if
c
c     Take away mean component and take square root for rms of pressure
c
      if (nxys.ge.32) then
         do y=1,nyp
            do x=1,nx
               pxys(x,y,32)=sqrt(max(0.,pxys(x,y,32)-pxys(x,y,31)**2))
            end do
         end do
      end if
c
c     Take away mean component for off-diagonal Reynolds-stress
c
      do ivar=13,15
         kvar=1
         lvar=2
         if (ivar.eq.14) lvar=3
         if (ivar.eq.15) kvar=3

         do y=1,nyp
            do x=1,nx
               pxys(x,y,ivar)=pxys(x,y,ivar)-
     &              pxys(x,y,kvar)*pxys(x,y,lvar)
            end do
         end do
      end do
c
c     Take away mean component for scalar-velocity
c
      if (scalar.gt.0.and.nxysth.ge.5) then
         do ivar=3,5
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=pxysth(x,y,ivar)-
     &                 pxysth(x,y,1)*pxys(x,y,ivar-2)
               end do
            end do
         end do
c
c     Calculate the derivatives of scalar-velocity
c
         do i=1,3
            call ddx(pxysth2(1,1,1+2*i),pxysth(1,1,2+i),xl)
            call ddy(pxysth2(1,1,2+2*i),pxysth(1,1,2+i),2./dstar)
         end do
      end if
c
c     Take away mean component for pressure-velocity
c
      if (nxys.ge.35) then
         do ivar=33,35
            do y=1,nyp
               do x=1,nx
                  pxys(x,y,ivar)=pxys(x,y,ivar)-
     &                 pxys(x,y,31)*pxys(x,y,ivar-32)
               end do
            end do
         end do
c
c     Calculate derivatives of pressure-velocity
c
         do i=1,3
            call ddx(pxys2(1,1,32+i),pxys(1,1,32+i),xl)
            call ddy(pxys2(1,1,35+i),pxys(1,1,32+i),2./dstar)
         end do
      end if
c
c     Take away the mean component for pressure-velocity_gradient
c
      if (nxys.ge.74) then
         do y=1,nyp
            do x=1,nx
c     pux
               pxys(x,y,36)=pxys(x,y,36)-
     &              pxys2(x,y,14)*pxys(x,y,31)
c     pvy
               pxys(x,y,37)=pxys(x,y,37)-
     &              pxys2(x,y,18)*pxys(x,y,31)
c     puy
               pxys(x,y,39)=pxys(x,y,39)-
     &              pxys2(x,y,17)*pxys(x,y,31)
c     pwx
               pxys(x,y,41)=pxys(x,y,41)-
     &              pxys2(x,y,16)*pxys(x,y,31)
c     pvx
               pxys(x,y,73)=pxys(x,y,73)-
     &              pxys2(x,y,15)*pxys(x,y,31)
c     pwy
               pxys(x,y,74)=pxys(x,y,74)-
     &              pxys2(x,y,19)*pxys(x,y,31)
            end do
         end do
c
c     Calculate the velocity-pressure_gradient
c
c     upx
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,39)=pxys2(x,y,33)-pxys(x,y,36)
            end do
         end do
c     vpx
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,40)=pxys2(x,y,34)-pxys(x,y,73)
            end do
         end do
c     upy
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,41)=pxys2(x,y,36)-pxys(x,y,39)
            end do
         end do
c     vpy
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,42)=pxys2(x,y,37)-pxys(x,y,37)
            end do
         end do
c     wpx,wpy,upz,vpz,wpz
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,110)=pxys2(x,y,35)-pxys(x,y,41)
               pxys2(x,y,111)=pxys2(x,y,38)-pxys(x,y,74)
               pxys2(x,y,112)=-pxys(x,y,42)
               pxys2(x,y,113)=-pxys(x,y,40)
               pxys2(x,y,114)=-pxys(x,y,38)
            end do
         end do
      end if
c
c     Take away mean component for velocity-scalar_gradient
c
      if (scalar.gt.0.and.nxysth.ge.30) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,22+(i-1)*3)=pxysth(x,y,22+(i-1)*3)
     &                 -pxysth2(x,y,1)*pxys(x,y,i)
                  pxysth(x,y,23+(i-1)*3)=pxysth(x,y,23+(i-1)*3)
     &                 -pxysth2(x,y,2)*pxys(x,y,i)
               end do
            end do
         end do
c
c     Calculate the scalar-velocity_gradient
c
         do i=1,3
            do j=1,2
               do y=1,nyp
                  do x=1,nx
                     pxysth2(x,y,8+(i-1)*3+j)=pxysth2(x,y,2*i+j)
     &                    -pxysth(x,y,21+(i-1)*3+j)
                  end do
               end do
            end do
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,11+(i-1)*3)=-pxysth(x,y,24+(i-1)*3)
               end do
            end do
         end do
      end if
c
c     x-twopoint product, subtract mean
c
      do ivar=16,18
         do y=1,nyp
            do x=1,nx
               pxys(x,y,ivar)=pxys(x,y,ivar)-
     &              pxys(x,y,ivar-15)*pxys(mod(x,nx)+1,y,ivar-15)
            end do
         end do
      end do
c
c     y-two point product
c     move down one point and subtract mean
c
      do ivar=19,21
         do y=2,nyp
            do x=1,nx
               pxys(x,y-1,ivar)=pxys(x,y,ivar)-
     &              pxys(x,y,ivar-18)*pxys(x,y-1,ivar-18)
            end do
         end do
      end do
c
c     z-two point product
c     subtract mean
c
      do ivar=22,24
         do y=1,nyp
            do x=1,nx
               pxys(x,y,ivar)=pxys(x,y,ivar)-pxys(x,y,ivar-21)**2
            end do
         end do
      end do
c
c     epsilon_ij = nu <ui,k uj,k>
c     (note that this definition
c     gives a value of epsilon_ii which is one
c     half of that by the definition in Mansour, Kim & Moin)
c
      do i=1,3
         do y=1,nyp
            do x=1,nx
c     re*epsilon_ii
               pxys(x,y,24+i)=pxys(x,y,24+i)-
     &              pxys2(x,y,13+i)**2-pxys2(x,y,16+i)**2
c     epsilon_ii
               pxys2(x,y,19+i)=pxys(x,y,24+i)/re
            end do
         end do
      end do
      do y=1,nyp
         do x=1,nx
c     re*epsilon_12
            pxys(x,y,28)=pxys(x,y,28)-
     &           pxys2(x,y,14)*pxys2(x,y,15)-
     &           pxys2(x,y,17)*pxys2(x,y,18)
c     epsilon_12
            pxys2(x,y,23)=pxys(x,y,28)/re
c     re*epsilon_13
            pxys(x,y,29)=pxys(x,y,29)-
     &           pxys2(x,y,14)*pxys2(x,y,16)-
     &           pxys2(x,y,17)*pxys2(x,y,19)
c     epsilon_13
            pxys2(x,y,24)=pxys(x,y,29)/re
c     re*epsilon_23
            pxys(x,y,30)=pxys(x,y,30)-
     &           pxys2(x,y,15)*pxys2(x,y,16)-
     &           pxys2(x,y,18)*pxys2(x,y,19)
c     epsilon_23
            pxys2(x,y,25)=pxys(x,y,30)/re
         end do
      end do
c
c     See below:
c     pxys2(x,y,26) = pxys2(x,y,20)+pxys2(x,y,21)+pxys2(x,y,22)
c
c     Viscous dissipation
c
      do y=1,nyp
         do x=1,nx
c     -2/Re<SijSij>
            pxys2(x,y,82)=-2./re*pxys(x,y,45)
c     -2/Re<Sij><Sij>
            pxys2(x,y,83)=-2./re*(pxys(x,y,46)**2+
     &           pxys(x,y,49)**2+pxys(x,y,51)**2+2.*pxys(x,y,47)**2+
     &           2.*pxys(x,y,48)**2+2.*pxys(x,y,50)**2)
c     -2/Re<Sij'Sij'>
            pxys2(x,y,84)=pxys2(x,y,82)-pxys2(x,y,83)
         end do
      end do
c
c     SGS dissipation
c
      do y=1,nyp
         do x=1,nx
c     <tau_ijSij>: total SGS dissipation
            pxys2(x,y,86)=pxys(x,y,53)
c     <tau_ij><Sij> mean SGS dissipation
            pxys2(x,y,87)=
     &           pxys(x,y,46)*pxys(x,y,54)+
     &           pxys(x,y,49)*pxys(x,y,57)+
     &           pxys(x,y,51)*pxys(x,y,59)+
     &           2.*pxys(x,y,47)*pxys(x,y,55)+
     &           2.*pxys(x,y,48)*pxys(x,y,56)+
     &           2.*pxys(x,y,50)*pxys(x,y,58)
c     <tau_ij'Sij'> fluctuating SGS dissipation
            pxys2(x,y,88)=pxys2(x,y,86)-pxys2(x,y,87)
         end do
      end do
c
c     SGS dissipation for scalar
c
      if (scalar.gt.0.and.nxysth.ge.36) then
         do y=1,nyp
            do x=1,nx
c     <sig_i dtheta/dx_i>: total SGS dissipation for scalar
               pxysth2(x,y,95)=pxysth(x,y,36)
c     <sig_i><dtheta/dx_i> mean SGS dissipation for scalar
               pxysth2(x,y,96)=pxysth(x,y,42)*pxysth2(x,y,1)+
     &              pxysth(x,y,43)*pxysth2(x,y,2)
c     <sig_i' dtheta/dx_i'> fluctuating SGS dissipation for scalar
               pxysth2(x,y,97)=pxysth2(x,y,95)-pxysth2(x,y,96)
            end do
         end do
      end if
c
c     SGS transport
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,265) = pxys(x,y,1)*pxys(x,y,54) +
     &           pxys(x,y,2) * pxys(x,y,55) +
     &           pxys(x,y,3) * pxys(x,y,56)
         end do
      end do
      call ddx(wxys(1,1,1),pxys2(1,1,265),xl)

      do y=1,nyp
         do x=1,nx
            pxys2(x,y,265) = pxys(x,y,1)*pxys(x,y,55) +
     &           pxys(x,y,2) * pxys(x,y,57) +
     &           pxys(x,y,3) * pxys(x,y,58)
         end do
      end do
      call ddy(wxys(1,1,2),pxys2(1,1,265),2./dstar)
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,266)=wxys(x,y,1)+wxys(x,y,2)
         end do
      end do

      call ddx(wxys(1,1,1),pxys(1,1,80),xl)
      call ddy(wxys(1,1,2),pxys(1,1,81),2./dstar)
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,265)=wxys(x,y,1)+wxys(x,y,2)
         end do
      end do

      do y=1,nyp
         do x=1,nx
            pxys2(x,y,267) = pxys2(x,y,265) - pxys2(x,y,266)
         end do
      end do

c
c     SGS transport for scalar
c
      if (scalar.gt.0.and.nxysth.ge.44) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,98) = pxysth(x,y,1)*pxysth(x,y,42)
            end do
         end do
         call ddx(wxys(1,1,1),pxysth2(1,1,98),xl)

         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,98) = pxysth(x,y,1)*pxysth(x,y,43)
            end do
         end do
         call ddy(wxys(1,1,2),pxysth2(1,1,98),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,99)=wxys(x,y,1)+wxys(x,y,2)
            end do
         end do

         call ddx(wxys(1,1,1),pxysth(1,1,39),xl)
         call ddy(wxys(1,1,2),pxysth(1,1,40),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,98)=wxys(x,y,1)+wxys(x,y,2)
            end do
         end do

         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,100) = pxysth2(x,y,98) - pxysth2(x,y,99)
            end do
         end do
      end if
c
c     Production
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,89) =
     &           -pxys(x,y,4)**2*pxys(x,y,46)
     &           -pxys(x,y,5)**2*pxys(x,y,49)
     &           -pxys(x,y,6)**2*pxys(x,y,51)
     &           -2.*pxys(x,y,13)*pxys(x,y,47)
     &           -2.*pxys(x,y,14)*pxys(x,y,48)
     &           -2.*pxys(x,y,15)*pxys(x,y,50)

            pxys2(x,y,90) =-pxys2(x,y,89)/(pxys2(x,y,84)+pxys2(x,y,88))

         end do
      end do
c
c     Some derived quantities
c     derivatives of the mean field are computed above
c
c     Turbulent kinetic energy k
c
       do y=1,nyp
          do x=1,nx
             pxys2(x,y,1)=(pxys(x,y,4)**2+pxys(x,y,5)**2+
     &            pxys(x,y,6)**2)/2.
          end do
       end do

       goto 1919

       if (nxys.ge.24) then
c
c     q2ro_ii (see Refinement of the Equation for the
c     determination of the turbulent micro-scale)
c     Jovanovic, Ye & Durst
c     LSTM 349/T 92
c     eq 42 page 10
c
          dx=xl/nx
c     first transform urms**2
          do i=1,3
             do y=1,nyp
                do x=1,nx
                   urmssh(x,y)=pxys(x,y,i+3)**2*(1./real(nx))
                end do
             end do
             call vrfftf(urmssh,urmssh(2,1),
     &            wxy2,wxy2(2,1),nx,nyp,2,nx+2,prex)
c     then interpolate data at x+1/2
             do x=1,nx/2+1
                do y=1,nyp
                   arg=0.5*dx*alfa(x)
                   h1=urmssh(2*x-1,y)*cos(arg)-urmssh(2*x,y)*sin(arg)
                   urmssh(2*x,y)=urmssh(2*x,y)*cos(arg)+
     &                  urmssh(2*x-1,y)*sin(arg)
                   urmssh(2*x-1,y)=h1
                end do
             end do
             call vrfftb(urmssh,urmssh(2,1),
     &            wxy2,wxy2(2,1),nx,nyp,2,nx+2,prex)
             do y=1,nyp
                do x=1,nx
c     -ddx2 of two-point product at x+1/2
                   pxys2(x,y,4+i)=
     &                  (2.*urmssh(x,y)-
     &                  2.*pxys(x,y,15+i))/(dx*dx)
                end do
             end do
          end do
c
c     -ddy2 of two-point product at y+1/2
c     first chebyshev transform urms**2
c
        do 2142 i=1,3
          do 2144 y=1,nyp
          do 2144 x=1,nx
             pxys2(x,nyp+1-y,7+i)=pxys(x,y,3+i)**2*(2./real(nyp-1))
 2144      continue
           call vchbf(pxys2(1,1,7+i),wxy,nyp,nx,nx,1,prey)
 2142   continue
c
c     then interpolate to midpoint of interval put in pxys(,,10+i)
c
        do 2146 i=1,3
        do 2146 y=1,nyp-1
          ym=.5*(eta(nyp+1-y)+eta(nyp+1-y-1))
c         call chebe2(pxys2(1,y,10+i),
c     &        eta(nyp),eta(1),pxys2(1,1,7+i),nx,nyp,ym)
 2146  continue
c
c     then calculate second derivative of two point correlation
c
      do 2150 i=1,3
      do 2150 y=1,nyp-1
      do 2150 x=1,nx
        dy=eta(nyp+1-y-1)-eta(nyp+1-y)
        pxys2(x,y,7+i)=(2.*pxys2(x,y,10+i)-
     &       2.*pxys(x,y,18+i))/(dy*dy)
 2150 continue
c
c     Top point : copy data from point below
c
      do 2152 i=1,3
      do 2152 x=1,nx
        pxys2(x,nyp,7+i)=pxys2(x,nyp-1,7+i)
 2152 continue
c
c     The above gives the correlation between nx,nx+1 and ny,ny+1
c     ie ro_ii is calculated at an intermediate point
c     to move back data to the grid we use linear interpolation :
c
      do 2160 i=1,3
      do 2160 y=1,nyp
      h1=.5*(pxys2(nx,y,4+i)+pxys2(1,y,4+i))
      do 2165 x=nx,2,-1
        pxys2(x,y,4+i)=.5*(pxys2(nx-1,y,4+i)+pxys2(x,y,4+i))
 2165 continue
      pxys2(1,y,4+i)=h1
 2160 continue
c
c     Interpolate in y, note that we leave the value nearest the upper and
c     lower boundary without interpolation
c     (for y=1/2) at the wall
c
      do 2170 i=1,3
      do 2170 y=nyp-1,2,-1
      do 2170 x=1,nx
        pxys2(x,y,7+i)=.5*(pxys2(x,y-1,7+i)+pxys2(x,y,7+i))
 2170 continue
c
c     -ddz2 of two-point product (invariant along z)
c
      dz=zl/nz
      do 2140 i=1,3
      do 2140 y=1,nyp
      do 2140 x=1,nx
        pxys2(x,y,10+i)=
     &       (2.*pxys(x,y,3+i)**2-2.*pxys(x,y,21+i))/(dz*dz)
 2140 continue
c
c     Now we have (minus) the second derivatives calculated
c     compute laplace and ro_ii
c
      do 2180 i=1,3
      do 2180 y=1,nyp
      do 2180 x=1,nx
        pxys2(x,y,1+i)=1./5.*
     &       (pxys2(x,y,4+i)+pxys2(x,y,7+i)+pxys2(x,y,10+i))
 2180 continue
      end if


 1919 continue
c
c     Calculate laplace of urms**2,vrms**2,wrms**2
c
      do i=1,3
         do y=1,nyp
            do x=1,nx
               wxy(x,y)=pxys(x,y,3+i)**2
            end do
         end do
         call ddx(wxy3,wxy,xl)
         call ddx(pxys2(1,1,4+i),wxy3,xl)
         call ddy(wxy3,wxy,2./dstar)
         call ddy(wxy,wxy3,2./dstar)
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,4+i)=pxys2(x,y,4+i)+wxy(x,y)
            end do
         end do
      end do
c
c     Microscale lambda from two-point correlation (eq 126, p21):
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,8)=sqrt(max(0.,2.*pxys2(x,y,1)/
     &           (pxys2(x,y,2)+pxys2(x,y,3)+pxys2(x,y,4))))
         end do
      end do
c
c     Re lambda from two-point correlation
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,9)=re*sqrt(2.*pxys2(x,y,1))*pxys2(x,y,8)
         end do
      end do
c
c     epsilon_ii eq(14) p5  and eq(51) p 11
c     derived from the two-point correlation
c     (note that this definition
c     gives a value of epsilon_ii which is one
c     half of that by the definition in Mansour, Kim & Moin)
c
      do i=1,3
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,9+i)=1./re*(.25*pxys2(x,y,4+i)+
     &              5.*pxys2(x,y,1+i))
            end do
         end do
      end do
c
c     epsilon : epsilon_ii from two-point correlation
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,13)=pxys2(x,y,10)+pxys2(x,y,11)+pxys2(x,y,12)
         end do
      end do
c
c     epsilon : epsilon_ii from derivatives
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,26)=pxys2(x,y,20)+pxys2(x,y,21)+pxys2(x,y,22)
         end do
      end do
c
c     Calculate derivatives of k,eps
c
      do i=1,2
         call ddx(pxys2(1,1,28+i),pxys2(1,1,i+(i-1)*24),xl)
         call ddy(pxys2(1,1,30+i),pxys2(1,1,i+(i-1)*24),2./dstar)
      end do
      call ddx(pxys2(1,1,59),pxys2(1,1,29),xl)
      call ddy(pxys2(1,1,60),pxys2(1,1,31),2./dstar)
c
c     -uv/(2k)
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,27)=-pxys(x,y,13)/(2.*pxys2(x,y,1))
         end do
      end do
c
c     Taylor microscale lambda (based on epsilon from derivatives)
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,28)=sqrt(max(0.,2.*pxys2(x,y,1)/
     &           (re*pxys2(x,y,26)/5.-
     &           1./20.*(pxys2(x,y,5)+pxys2(x,y,6)+pxys2(x,y,7)))))
         end do
      end do

c re_lambda (based on epsilon from derivatives)
c      do 3000 y=1,nyp
c      do 3000 x=1,nx
c        pxys2(x,y,29)=re*sqrt(2.*pxys2(x,y,1))*pxys2(x,y,28)
c 3000 continue
c      end if
c S
c      do 3100 y=1,nyp
c      do 3100 x=1,nx
c        pxys2(x,y,30)=pxys2(x,y,15)+pxys2(x,y,17)
c 3100 continue
c      if (nxys.ge.30) then
c S_lambda
c      do 3200 y=1,nyp
c      do 3200 x=1,nx
c        pxys2(x,y,31)=pxys2(x,y,30)*pxys2(x,y,28)/sqrt(2.*pxys2(x,y,1))
c 3200 continue
c S_lambda*re_lambda
c      do 3300 y=1,nyp
c      do 3300 x=1,nx
c        pxys2(x,y,32)=re*pxys2(x,y,30)*pxys2(x,y,28)**2
c      3300 continue

      do y=1,nyp
         do x=1,nx
            wxy(x,y)=pxys(x,y,4)**2-pxys(x,y,5)**2
            pxys2(x,y,44)=wxy(x,y)*pxys2(x,y,26)/
     &           pxys2(x,y,1)*pxys2(x,y,14)
         end do
      end do
c     <urms2-vrms2>x
      call ddx(pxys2(1,1,43),wxy,xl)
c     <uv>yyy
      call ddy(wxy,pxys(1,1,13),2./dstar)
      call ddy(wxy3,wxy,2./dstar)
      call ddy(pxys2(1,1,45),wxy3,2./dstar)
c     <ww>xx,<ww>yy,<ww>x,<ww>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,6)**2
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,100),wxy3,xl)
      call ddx(pxys2(1,1,46),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,101),wxy3,2./dstar)
      call ddy(pxys2(1,1,99),wxy,2./dstar)
c     <uu>xx,<uu>yy,<uu>x,<uu>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,4)**2
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,53),wxy3,xl)
      call ddx(pxys2(1,1,47),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,55),wxy3,2./dstar)
      call ddy(pxys2(1,1,48),wxy,2./dstar)
c     <vv>xx,<vv>yy,<vv>x,<vv>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,5)**2
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,54),wxy3,xl)
      call ddx(pxys2(1,1,49),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,56),wxy3,2./dstar)
      call ddy(pxys2(1,1,50),wxy,2./dstar)
c     <uv>xx,<uv>yy,<uv>x,<uv>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,13)
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,57),wxy3,xl)
      call ddx(pxys2(1,1,51),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,58),wxy3,2./dstar)
      call ddy(pxys2(1,1,52),wxy,2./dstar)
c     <uw>xx,<uw>yy,<uw>x,<uw>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,14)
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,106),wxy3,xl)
      call ddx(pxys2(1,1,102),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,107),wxy3,2./dstar)
      call ddy(pxys2(1,1,103),wxy,2./dstar)
c     <vw>xx,<vw>yy,<vw>x,<vw>y
      do y=1,nyp
         do x=1,nx
            wxy3(x,y)=pxys(x,y,15)
         end do
      end do
      call ddx(wxy,wxy3,xl)
      call ddx(pxys2(1,1,108),wxy3,xl)
      call ddx(pxys2(1,1,104),wxy,xl)
      call ddy(wxy,wxy3,2./dstar)
      call ddy(pxys2(1,1,109),wxy3,2./dstar)
      call ddy(pxys2(1,1,105),wxy,2./dstar)
c
c     Skewness and flatness for the velocities
c
      if (nxys.ge.65) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,116+i)=pxys(x,y,62+i)/pxys(x,y,3+i)**3
               end do
            end do
         end do
      end if
      if (nxys.ge.77) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,119+i)=pxys(x,y,74+i)/pxys(x,y,3+i)**4
               end do
            end do
c     overwrite the wall point with the closest point to the wall
            do x=1,nx
               pxys2(x,1,119+i) = pxys2(x,2,119+i)
            end do
         end do
      end if
c
c     Skewness and flatness for the scalar
c
      if (scalar.gt.0.and.nxysth.ge.34) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,18)=pxysth(x,y,34)/pxysth(x,y,2)**3
            end do
         end do
         do x=1,nx
            do y=1,2
               pxysth2(x,y,18)=0.
            end do
         end do
      end if
      if (scalar.gt.0.and.nxysth.ge.35) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,19)=pxysth(x,y,35)/pxysth(x,y,2)**4
            end do
         end do
         do x=1,nx
            do y=1,2
               pxysth2(x,y,19)=0.
            end do
         end do
      end if
c
c     Skewness and flatness for the pressure
c
      if (nxys.ge.78) then
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,123)=pxys(x,y,78)/pxys(x,y,32)**3
            end do
         end do
      end if
      if (nxys.ge.79) then
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,124)=pxys(x,y,79)/pxys(x,y,32)**4
            end do
         end do
      end if
c
c     Turbulent Pr
c
      if (scalar.gt.0.and.nxysth.ge.5) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,20)=pxys(x,y,13)*pxysth2(x,y,2)
     &              /(pxysth(x,y,4)*pxys2(x,y,17))
            end do
         end do
      end if
c
c     Turbulent kinetic energy for scalar
c
      if (scalar.gt.0.and.nxysth.ge.2) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,21)=pxysth(x,y,2)**2/2.
            end do
         end do
c
c     (k_theta)x,(k_theta)y,(k_theta)xx,(k_theta)yy
c
         call ddx(pxysth2(1,1,22),pxysth2(1,1,21),xl)
         call ddy(pxysth2(1,1,23),pxysth2(1,1,21),2./dstar)
         call ddx(pxysth2(1,1,24),pxysth2(1,1,22),xl)
         call ddy(pxysth2(1,1,25),pxysth2(1,1,23),2./dstar)
      end if
c
c     Mean convection for uiui
c
      if (nxys.ge.71) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)**2*pxys(x,y,1)
                  wxys(x,y,i+3)=pxys(x,y,i)**2*pxys(x,y,2)
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               wxys(x,y,7)=wxys(x,y,1)+wxys(x,y,2)+wxys(x,y,3)
               wxys(x,y,8)=wxys(x,y,4)+wxys(x,y,5)+wxys(x,y,6)
            end do
         end do
         call ddx(wxys(1,1,1),wxys(1,1,7),xl)
         call ddy(wxys(1,1,2),wxys(1,1,8),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,126)=(wxys(x,y,1)+wxys(x,y,2))/2.
            end do
         end do
c
c     Fluctuating convection for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,127)=pxys(x,y,1)*pxys2(x,y,29)+
     &              pxys(x,y,2)*pxys2(x,y,31)
            end do
         end do
c
c     Fluctuating molecular diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
c     pxys2(x,y,130)=(pxys2(x,y,59)+pxys2(x,y,60))/re
               pxys2(x,y,130)=(pxys2(x,y,46)+pxys2(x,y,47)+pxys2(x,y,48)
     &              +pxys2(x,y,49)+pxys2(x,y,50)+pxys2(x,y,99))/(2.*re)
            end do
         end do
c
c     Mean molecular diffusion for uiui
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)**2
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               wxys(x,y,4)=wxys(x,y,1)+wxys(x,y,2)+wxys(x,y,3)
            end do
         end do
         call ddx(wxys(1,1,5),wxys(1,1,4),xl)
         call ddx(wxys(1,1,6),wxys(1,1,5),xl)
         call ddy(wxys(1,1,7),wxys(1,1,4),2./dstar)
         call ddy(wxys(1,1,8),wxys(1,1,7),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,129)=(wxys(x,y,6)+wxys(x,y,8))/(2.*re)
            end do
         end do
c
c     Total molecular diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,128)=pxys2(x,y,129)+pxys2(x,y,130)
            end do
         end do
c
c     Fluctuating turbulent diffusion for uiui
c
         call ddx(wxys(1,1,1),pxys(1,1,63),xl)
         call ddy(wxys(1,1,2),pxys(1,1,66),2./dstar)
         call ddx(wxys(1,1,3),pxys(1,1,68),xl)
         call ddy(wxys(1,1,4),pxys(1,1,64),2./dstar)
         call ddx(wxys(1,1,5),pxys(1,1,70),xl)
         call ddy(wxys(1,1,6),pxys(1,1,71),2./dstar)
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,132)=pxys2(x,y,132)-wxys(x,y,i)/2.
               end do
            end do
         end do
c
c     Mean turbulent diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1) =pxys(x,y,1)*pxys(x,y,4)**2
               wxys(x,y,3) =pxys(x,y,1)*pxys(x,y,13)
               wxys(x,y,5) =pxys(x,y,2)*pxys(x,y,13)
               wxys(x,y,7) =pxys(x,y,2)*pxys(x,y,5)**2
               wxys(x,y,9) =pxys(x,y,3)*pxys(x,y,14)
               wxys(x,y,11)=pxys(x,y,3)*pxys(x,y,15)
            end do
         end do
         do i=1,9,4
            call ddx(wxys(1,1,i+1),wxys(1,1,i),xl)
            call ddy(wxys(1,1,i+3),wxys(1,1,i+2),2./dstar)
         end do
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,131)=pxys2(x,y,131)-wxys(x,y,2*i)
               end do
            end do
         end do
c
c     Mean velocity-pressure_gradient for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,134)=-(pxys(x,y,1)*pxys2(x,y,115)
     &              +pxys(x,y,2)*pxys2(x,y,116))
            end do
         end do
c
c     Fluctuating velocity-pressure_gradient for uiui
c
         do y=1,nyp
            do x=1,nx
c               pxys2(x,y,135)=-(pxys2(x,y,39)+pxys2(x,y,42)
c     &              +pxys2(x,y,114))
               pxys2(x,y,135)=-(pxys2(x,y,33)+pxys2(x,y,37))
            end do
         end do
c
c     Total velocity-pressure_gradient for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,133)= pxys2(x,y,134)+pxys2(x,y,135)
            end do
         end do
c
c     Fluctuating dissipation for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,138)=-pxys2(x,y,26)
            end do
         end do
c
c     Mean dissipation for uiui
c
         do i=14,19
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,137)=pxys2(x,y,137)-(pxys2(x,y,i)**2)/re
               end do
            end do
         end do
c
c     Total dissipation for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,136)=pxys2(x,y,137)+pxys2(x,y,138)
            end do
         end do
c
c     Fluctuating viscous diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,141)=pxys2(x,y,130)
     &              -pxys2(x,y,84)+pxys2(x,y,138)
            end do
         end do
c
c     Mean viscous diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,140)=pxys2(x,y,129)
     &              -pxys2(x,y,83)+pxys2(x,y,137)
            end do
         end do
c
c     Total viscous diffusion for uiui
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,139)=pxys2(x,y,140)+pxys2(x,y,141)
            end do
         end do
      end if
c
c     Mean convection for scalar squared
c
      if (scalar.gt.0.and.nxysth.ge.11) then
         do i=1,2
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxysth(x,y,1)**2*pxys(x,y,i)
               end do
            end do
         end do
         call ddx(wxys(1,1,3),wxys(1,1,1),xl)
         call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,27)=(wxys(x,y,3)+wxys(x,y,4))/2.
            end do
         end do
c
c     Fluctuating convection for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,28)=pxysth2(x,y,22)*pxys(x,y,1)+
     &              pxysth2(x,y,23)*pxys(x,y,2)
            end do
         end do
c
c     Mean molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxysth(x,y,1)**2
            end do
         end do
         call ddx(wxys(1,1,2),wxys(1,1,1),xl)
         call ddx(wxys(1,1,3),wxys(1,1,2),xl)
         call ddy(wxys(1,1,4),wxys(1,1,1),2./dstar)
         call ddy(wxys(1,1,5),wxys(1,1,4),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,30)=(wxys(x,y,3)+wxys(x,y,5))
     &              /(2.*re*pr(scalarind))
            end do
         end do
c
c     Fluctuating molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,31)=(pxysth2(x,y,24)+pxysth2(x,y,25))
     &              /(re*pr(scalarind))
            end do
         end do
c
c     Total molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,29)=pxysth2(x,y,30)+pxysth2(x,y,31)
            end do
         end do
c
c     Mean turbulent diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxysth(x,y,1)*pxysth(x,y,3)
               wxys(x,y,2)=pxysth(x,y,1)*pxysth(x,y,4)
            end do
         end do
         call ddx(wxys(1,1,3),wxys(1,1,1),xl)
         call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,32)=-(wxys(x,y,3)+wxys(x,y,4))
            end do
         end do
c
c     Fluctuating turbulent diffusion for scalar squared
c
         call ddx(wxys(1,1,1),pxysth(1,1,6),xl)
         call ddy(wxys(1,1,2),pxysth(1,1,7),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,33)=-(wxys(x,y,1)+wxys(x,y,2))/2.
            end do
         end do
c
c     Production for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,34)=-pxysth(x,y,3)*pxysth2(x,y,1)-
     &              pxysth(x,y,4)*pxysth2(x,y,2)
            end do
         end do
c
c     Total dissipation for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,35)=-1./(re*pr(scalarind))
     &              *(pxysth(x,y,9)+pxysth(x,y,10)+pxysth(x,y,11))
            end do
         end do
c
c     Take away mean component for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,37)=pxysth2(x,y,35)+1./(re*pr(scalarind))*
     &              (pxysth2(x,y,1)**2+pxysth2(x,y,2)**2)
            end do
         end do
c
c     Mean dissipation for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,36)=pxysth2(x,y,35)-pxysth2(x,y,37)
            end do
         end do
      end if
c
c     Mean convection for uiuj (148-153)
c
      if (nxys.ge.72.and.nxys2.ge.148) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)**2
               end do
            end do
            call ddx(wxys(1,1,i+3),wxys(1,1,i),xl)
            call ddy(wxys(1,1,i+6),wxys(1,1,i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,147+i)=wxys(x,y,i+3)*pxys(x,y,1)
     &                 +wxys(x,y,i+6)*pxys(x,y,2)
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,1)*pxys(x,y,2)
               wxys(x,y,2)=pxys(x,y,1)*pxys(x,y,3)
               wxys(x,y,3)=pxys(x,y,2)*pxys(x,y,3)
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,i+3),wxys(1,1,i),xl)
            call ddy(wxys(1,1,i+6),wxys(1,1,i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,150+i)=wxys(x,y,i+3)*pxys(x,y,1)
     &                 +wxys(x,y,i+6)*pxys(x,y,2)
               end do
            end do
         end do
c
c     Fluctuating convection for uiuj (154-159)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,154)=pxys2(x,y,53)*pxys(x,y,1)+
     &              pxys2(x,y,55)*pxys(x,y,2)
               pxys2(x,y,155)=pxys2(x,y,54)*pxys(x,y,1)+
     &              pxys2(x,y,56)*pxys(x,y,2)
               pxys2(x,y,156)=pxys2(x,y,100)*pxys(x,y,1)+
     &              pxys2(x,y,101)*pxys(x,y,2)
               pxys2(x,y,157)=pxys2(x,y,57)*pxys(x,y,1)+
     &              pxys2(x,y,58)*pxys(x,y,2)
               pxys2(x,y,158)=pxys2(x,y,106)*pxys(x,y,1)+
     &              pxys2(x,y,107)*pxys(x,y,2)
               pxys2(x,y,159)=pxys2(x,y,108)*pxys(x,y,1)+
     &              pxys2(x,y,109)*pxys(x,y,2)
            end do
         end do
c
c    Mean molecular diffusion for uiuj (166-171)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)**2
               end do
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,i+3),wxys(1,1,i),xl)
            call ddx(wxys(1,1,i+6),wxys(1,1,i+3),xl)
            call ddy(wxys(1,1,i+3),wxys(1,1,i),2./dstar)
            call ddy(wxys(1,1,i+9),wxys(1,1,i+3),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,165+i)=(wxys(x,y,i+6)+wxys(x,y,i+9))/re
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,1)*pxys(x,y,2)
               wxys(x,y,2)=pxys(x,y,1)*pxys(x,y,3)
               wxys(x,y,3)=pxys(x,y,2)*pxys(x,y,3)
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,i+3),wxys(1,1,i),xl)
            call ddx(wxys(1,1,i+6),wxys(1,1,i+3),xl)
            call ddy(wxys(1,1,i+3),wxys(1,1,i),2./dstar)
            call ddy(wxys(1,1,i+9),wxys(1,1,i+3),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,168+i)=(wxys(x,y,i+6)+wxys(x,y,i+9))/re
               end do
            end do
         end do
c
c     Fluctuating molecular diffusion for uiuj (172-177)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,171+i)=pxys2(x,y,4+i)/re
               end do
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,i),pxys(1,1,i+12),xl)
            call ddx(wxys(1,1,i+3),wxys(1,1,i),xl)
            call ddy(wxys(1,1,i),pxys(1,1,i+12),2./dstar)
            call ddy(wxys(1,1,i+6),wxys(1,1,i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,174+i)=(wxys(x,y,i+3)+wxys(x,y,i+6))/re
               end do
            end do
         end do
c
c     Total molecular diffusion for uiuj (160-165)
c
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,159+i)=pxys2(x,y,165+i)+pxys2(x,y,171+i)
               end do
            end do
         end do
c
c     Mean turbulent diffusion for uiuj (178-183)
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,4)**2*pxys(x,y,1)
               wxys(x,y,2)=pxys(x,y,13)*pxys(x,y,1)
               wxys(x,y,3)=pxys(x,y,13)*pxys(x,y,2)
               wxys(x,y,4)=pxys(x,y,5)**2*pxys(x,y,2)
               wxys(x,y,5)=pxys(x,y,14)*pxys(x,y,3)
               wxys(x,y,6)=pxys(x,y,15)*pxys(x,y,3)
               wxys(x,y,7)=pxys(x,y,13)*pxys(x,y,1)
               wxys(x,y,9)=pxys(x,y,4)**2*pxys(x,y,2)
               wxys(x,y,8)=pxys(x,y,5)**2*pxys(x,y,1)
               wxys(x,y,10)=pxys(x,y,13)*pxys(x,y,2)
               wxys(x,y,11)=pxys(x,y,14)*pxys(x,y,1)
               wxys(x,y,13)=pxys(x,y,4)**2*pxys(x,y,3)
               wxys(x,y,12)=pxys(x,y,15)*pxys(x,y,1)
               wxys(x,y,14)=pxys(x,y,13)*pxys(x,y,3)
               wxys(x,y,15)=pxys(x,y,14)*pxys(x,y,2)
               wxys(x,y,17)=pxys(x,y,13)*pxys(x,y,3)
               wxys(x,y,16)=pxys(x,y,15)*pxys(x,y,2)
               wxys(x,y,18)=pxys(x,y,5)**2*pxys(x,y,3)
            end do
         end do
         do i=1,9
            call ddx(wxys(1,1,19+2*i),wxys(1,1,2*i-1),xl)
            call ddy(wxys(1,1,20+2*i),wxys(1,1,2*i),2./dstar)
         end do
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,177+i)=-2.*(wxys(x,y,19+2*i)
     &                 +wxys(x,y,20+2*i))
               end do
            end do
         end do
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,181)=-(wxys(x,y,27)+wxys(x,y,28)+
     &              wxys(x,y,29)+wxys(x,y,30))
               pxys2(x,y,182)=-(wxys(x,y,31)+wxys(x,y,32)+
     &              wxys(x,y,33)+wxys(x,y,34))
               pxys2(x,y,183)=-(wxys(x,y,35)+wxys(x,y,36)+
     &              wxys(x,y,37)+wxys(x,y,38))
            end do
         end do
c
c     Fluctuating turbulent diffusion for uiuj (184-189)
c
         call ddx(wxys(1,1,1),pxys(1,1,63),xl)
         call ddy(wxys(1,1,2),pxys(1,1,66),2./dstar)
         call ddx(wxys(1,1,3),pxys(1,1,68),xl)
         call ddy(wxys(1,1,4),pxys(1,1,64),2./dstar)
         call ddx(wxys(1,1,5),pxys(1,1,70),xl)
         call ddy(wxys(1,1,6),pxys(1,1,71),2./dstar)
         call ddx(wxys(1,1,7),pxys(1,1,66),xl)
         call ddy(wxys(1,1,8),pxys(1,1,68),2./dstar)
         call ddx(wxys(1,1,9),pxys(1,1,67),xl)
         call ddy(wxys(1,1,10),pxys(1,1,72),2./dstar)
         call ddx(wxys(1,1,11),pxys(1,1,72),xl)
         call ddy(wxys(1,1,12),pxys(1,1,69),2./dstar)
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,183+i)=-(wxys(x,y,2*i-1)+wxys(x,y,2*i))
               end do
            end do
         end do
c
c     Production for uiuj (190-195)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,190)=-2.*(pxys(x,y,4)**2*pxys2(x,y,14)+
     &              pxys(x,y,13)*pxys2(x,y,17))
               pxys2(x,y,191)=-2.*(pxys(x,y,5)**2*pxys2(x,y,18)+
     &              pxys(x,y,13)*pxys2(x,y,15))
               pxys2(x,y,192)=-2.*(pxys(x,y,14)*pxys2(x,y,16)+
     &              pxys(x,y,15)*pxys2(x,y,19))
               pxys2(x,y,193)=-(pxys(x,y,13)*pxys2(x,y,14)+
     &              pxys(x,y,5)**2*pxys2(x,y,17)+
     &              pxys(x,y,4)**2*pxys2(x,y,15)+
     &              pxys(x,y,13)*pxys2(x,y,18))
               pxys2(x,y,194)=-(pxys(x,y,14)*pxys2(x,y,14)+
     &              pxys(x,y,15)*pxys2(x,y,17)+
     &              pxys(x,y,4)**2*pxys2(x,y,16)+
     &              pxys(x,y,13)*pxys2(x,y,19))
               pxys2(x,y,195)=-(pxys(x,y,14)*pxys2(x,y,15)+
     &              pxys(x,y,15)*pxys2(x,y,18)+
     &              pxys(x,y,13)*pxys2(x,y,16)+
     &              pxys(x,y,5)**2*pxys2(x,y,19))
            end do
         end do
c
c     Mean velocity pressure_gradient for uiuj (202-207)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,202)=-2.*pxys(x,y,1)*pxys2(x,y,115)
               pxys2(x,y,203)=-2.*pxys(x,y,2)*pxys2(x,y,116)
               pxys2(x,y,204)=0.
               pxys2(x,y,205)=-(pxys(x,y,1)*pxys2(x,y,116)+
     &              pxys(x,y,2)*pxys2(x,y,115))
               pxys2(x,y,206)=-pxys(x,y,3)*pxys2(x,y,115)
               pxys2(x,y,207)=-pxys(x,y,3)*pxys2(x,y,116)
            end do
         end do
c
c     Fluctuating velocity pressure_gradient for uiuj (208-213)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,208)=-2.*pxys2(x,y,39)
               pxys2(x,y,209)=-2.*pxys2(x,y,42)
               pxys2(x,y,210)=-2.*pxys2(x,y,114)
               pxys2(x,y,211)=-(pxys2(x,y,41)+pxys2(x,y,40))
               pxys2(x,y,212)=-(pxys2(x,y,112)+pxys2(x,y,110))
               pxys2(x,y,213)=-(pxys2(x,y,113)+pxys2(x,y,111))
            end do
         end do
c
c     Total velocity pressure_gradient for uiuj (196-201)
c
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,195+i)=pxys2(x,y,201+i)+pxys2(x,y,207+i)
               end do
            end do
         end do
c
c     Mean (pressure-velocity)_gradient for uiuj (219-223)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)*pxys(x,y,31)
               end do
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,2+2*i),wxys(1,1,i),xl)
            call ddy(wxys(1,1,3+2*i),wxys(1,1,i),2./dstar)
         end do
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,219)=-2.*wxys(x,y,4)
               pxys2(x,y,220)=-2.*wxys(x,y,5)
               pxys2(x,y,221)=-(wxys(x,y,6)+wxys(x,y,7))
               pxys2(x,y,222)=-wxys(x,y,8)
               pxys2(x,y,223)=-wxys(x,y,9)
            end do
         end do
c
c     Fluctuating (pressure-velocity)_gradient for uiuj (224-228)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,224)=-2.*pxys2(x,y,33)
               pxys2(x,y,225)=-2.*pxys2(x,y,37)
               pxys2(x,y,226)=-(pxys2(x,y,34)+pxys2(x,y,36))
               pxys2(x,y,227)=-pxys2(x,y,35)
               pxys2(x,y,228)=-pxys2(x,y,38)
            end do
         end do
c
c     Total (pressure-velocity)_gradient for uiuj (214-218) (<ww> is zero)
c
         do i=1,5
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,213+i)=pxys2(x,y,218+i)+pxys2(x,y,223+i)
               end do
            end do
         end do
c
c     Mean pressure-strain rate for uiuj (235-240)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,235)=2*pxys(x,y,46)*pxys(x,y,31)
               pxys2(x,y,236)=2*pxys(x,y,49)*pxys(x,y,31)
               pxys2(x,y,237)=2*pxys(x,y,51)*pxys(x,y,31)
               pxys2(x,y,238)=2*pxys(x,y,47)*pxys(x,y,31)
               pxys2(x,y,239)=2*pxys(x,y,48)*pxys(x,y,31)
               pxys2(x,y,240)=2*pxys(x,y,50)*pxys(x,y,31)
            end do
         end do
c
c     Fluctuating pressure-strain rate for uiuj (241-246)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,241)=2*pxys(x,y,36)
               pxys2(x,y,242)=2*pxys(x,y,37)
               pxys2(x,y,243)=2*pxys(x,y,38)
               pxys2(x,y,244)=pxys(x,y,39)+pxys(x,y,73)
               pxys2(x,y,245)=pxys(x,y,42)+pxys(x,y,41)
               pxys2(x,y,246)=pxys(x,y,40)+pxys(x,y,74)
            end do
         end do
c
c     Total pressure-strain rate for uiuj (229-234)
c
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,228+i)=pxys2(x,y,234+i)+pxys2(x,y,240+i)
               end do
            end do
         end do
c
c     Fluctuating epsilon dissipation for uiuj (259-264)
c
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,258+i)=-2.*pxys2(x,y,19+i)
               end do
            end do
         end do
c
c     Mean epsilon dissipation for uiuj (253-258)
c
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,253)=-2.*(pxys2(x,y,14)**2+pxys2(x,y,17)**2)/re
               pxys2(x,y,254)=-2.*(pxys2(x,y,15)**2+pxys2(x,y,18)**2)/re
               pxys2(x,y,255)=-2.*(pxys2(x,y,16)**2+pxys2(x,y,19)**2)/re
               pxys2(x,y,256)=-2.*(pxys2(x,y,14)*pxys2(x,y,15)+
     &              pxys2(x,y,17)*pxys2(x,y,18))/re
               pxys2(x,y,257)=-2.*(pxys2(x,y,14)*pxys2(x,y,16)+
     &              pxys2(x,y,17)*pxys2(x,y,19))/re
               pxys2(x,y,258)=-2.*(pxys2(x,y,15)*pxys2(x,y,16)+
     &              pxys2(x,y,18)*pxys2(x,y,19))/re
            end do
         end do
c
c     Total epsilon dissipation for uiuj (247-252)
c
         do i=1,6
            do y=1,nyp
               do x=1,nx
                  pxys2(x,y,246+i)=pxys2(x,y,252+i)+pxys2(x,y,258+i)
               end do
            end do
         end do
      end if
c
c     Mean convection for ui-theta (41-43)
c
      if (scalar.gt.0.and.nxysth.ge.33) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)*pxysth(x,y,1)
               end do
            end do
            call ddx(wxys(1,1,3+i),wxys(1,1,i),xl)
            call ddy(wxys(1,1,6+i),wxys(1,1,i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,40+i)=wxys(x,y,3+i)*pxys(x,y,1)+
     &                 wxys(x,y,i+6)*pxys(x,y,2)
               end do
            end do
         end do
c
c     Fluctuating convection for ui-theta (44-46)
c
         do ivar=44,46
            do i=1,2
               do y=1,nyp
                  do x=1,nx
                     pxysth2(x,y,ivar)=pxysth2(x,y,ivar)+
     &                    pxysth2(x,y,2+i)*pxys(x,y,i)
                  end do
               end do
            end do
         end do
c
c     Mean molecular diffusion for ui-theta (50-52)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,4*i-3)=pxysth2(x,y,1)*pxys(x,y,i)
                  wxys(x,y,4*i-2)=pxysth2(x,y,2)*pxys(x,y,i)
                  wxys(x,y,4*i-1)=pxys2(x,y,13+i)*pxysth(x,y,1)
                  wxys(x,y,4*i)=pxys2(x,y,16+i)*pxysth(x,y,1)
               end do
            end do
            call ddx(wxys(1,1,4*i+9),wxys(1,1,4*i-3),xl)
            call ddy(wxys(1,1,4*i+10),wxys(1,1,4*i-2),2./dstar)
            call ddx(wxys(1,1,4*i+11),wxys(1,1,4*i-1),xl)
            call ddy(wxys(1,1,4*i+12),wxys(1,1,4*i),2./dstar)
         end do
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,50)=(wxys(x,y,13)+wxys(x,y,14))
     &              /(re*pr(scalarind))+(wxys(x,y,15)+wxys(x,y,16))/re
               pxysth2(x,y,51)=(wxys(x,y,17)+wxys(x,y,18))
     &              /(re*pr(scalarind))+(wxys(x,y,19)+wxys(x,y,20))/re
               pxysth2(x,y,52)=(wxys(x,y,21)+wxys(x,y,22))
     &              /(re*pr(scalarind))+(wxys(x,y,23)+wxys(x,y,24))/re
            end do
         end do
c
c     Fluctuating molecular diffusion for ui-theta (53-55)
c
         call ddx(wxys(1,1,1 ),pxysth(1,1,22),xl)
         call ddy(wxys(1,1,2 ),pxysth(1,1,23),2./dstar)
         call ddx(wxys(1,1,3 ),pxysth2(1,1,9),xl)
         call ddy(wxys(1,1,4 ),pxysth2(1,1,10),2./dstar)
         call ddx(wxys(1,1,5 ),pxysth(1,1,25),xl)
         call ddy(wxys(1,1,6 ),pxysth(1,1,26),2./dstar)
         call ddx(wxys(1,1,7 ),pxysth2(1,1,12),xl)
         call ddy(wxys(1,1,8 ),pxysth2(1,1,13),2./dstar)
         call ddx(wxys(1,1,9 ),pxysth(1,1,28),xl)
         call ddy(wxys(1,1,10),pxysth(1,1,29),2./dstar)
         call ddx(wxys(1,1,11),pxysth2(1,1,15),xl)
         call ddy(wxys(1,1,12),pxysth2(1,1,16),2./dstar)
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,52+i)=(wxys(x,y,4*i-3)+wxys(x,y,4*i-2))
     &                 /(re*pr(scalarind))+(wxys(x,y,4*i-3)
     &                 +wxys(x,y,4*i))/re
               end do
            end do
         end do
c
c    Total molecular diffusion for ui-theta (47-49)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,46+i)=pxysth2(x,y,49+i)+pxysth2(x,y,52+i)
               end do
            end do
         end do
c
c     Mean turbulent diffusion for ui-theta (56-58)
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,1)*pxysth(x,y,3)+
     &              pxys(x,y,4)**2*pxysth(x,y,1)
               wxys(x,y,2)=pxys(x,y,1)*pxysth(x,y,4)+
     &              pxys(x,y,13)*pxysth(x,y,1)
               wxys(x,y,3)=pxys(x,y,2)*pxysth(x,y,3)+
     &              pxys(x,y,13)*pxysth(x,y,1)
               wxys(x,y,4)=pxys(x,y,2)*pxysth(x,y,4)+
     &              pxys(x,y,5)**2*pxysth(x,y,1)
               wxys(x,y,5)=pxys(x,y,3)*pxysth(x,y,3)+
     &              pxys(x,y,14)*pxysth(x,y,1)
               wxys(x,y,6)=pxys(x,y,3)*pxysth(x,y,4)+
     &              pxys(x,y,15)*pxysth(x,y,1)
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,2*i-1),wxys(1,1,2*i-1),xl)
            call ddy(wxys(1,1,2*i),wxys(1,1,2*i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,55+i)=-(wxys(x,y,2*i-1)+wxys(x,y,2*i))
               end do
            end do
         end do
c
c     Fluctuating turbulent diffusion for ui-theta (59-61)
c
         call ddx(wxys(1,1,1),pxysth(1,1,16),xl)
         call ddy(wxys(1,1,2),pxysth(1,1,19),2./dstar)
         call ddx(wxys(1,1,3),pxysth(1,1,19),xl)
         call ddy(wxys(1,1,4),pxysth(1,1,17),2./dstar)
         call ddx(wxys(1,1,5),pxysth(1,1,20),xl)
         call ddy(wxys(1,1,6),pxysth(1,1,21),2./dstar)
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,58+i)=-(wxys(x,y,2*i-1)+wxys(x,y,2*i))
               end do
            end do
         end do
c
c     Production for ui-theta (62-64)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,62)=-(pxysth(x,y,3)*pxys2(x,y,14)
     &              +pxysth(x,y,4)*pxys2(x,y,17)+pxysth2(x,y,1)
     &              *pxys(x,y,4)**2+pxysth2(x,y,2)*pxys(x,y,13))
               pxysth2(x,y,63)=-(pxysth(x,y,3)*pxys2(x,y,15)
     &              +pxysth(x,y,4)*pxys2(x,y,18)+pxysth2(x,y,1)
     &              *pxys(x,y,13)+pxysth2(x,y,2)*pxys(x,y,5)**2)
               pxysth2(x,y,64)=-(pxysth(x,y,3)*pxys2(x,y,16)
     &              +pxysth(x,y,4)*pxys2(x,y,19)+pxysth2(x,y,1)
     &              *pxys(x,y,14)+pxysth2(x,y,2)*pxys(x,y,15))
            end do
         end do
c
c     Total (scalar-pressure)_gradient for ui-theta (65-66) (<pth>z is zreo)
c
         call ddx(pxysth2(1,1,65),pxysth(1,1,12),xl)
         call ddy(pxysth2(1,1,66),pxysth(1,1,12),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,65)=-pxysth2(x,y,65)
               pxysth2(x,y,66)=-pxysth2(x,y,66)
            end do
         end do
c
c     Take away mean component for pressure-scalar
c
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,12)=pxysth(x,y,12)-pxys(x,y,31)*pxysth(x,y,1)
            end do
         end do
c
c     Fluctuating (scalar-pressure)_gradient for ui-theta (69-70)
c
         call ddx(pxysth2(1,1,69),pxysth(1,1,12),xl)
         call ddy(pxysth2(1,1,70),pxysth(1,1,12),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,69)=-pxysth2(x,y,69)
               pxysth2(x,y,70)=-pxysth2(x,y,70)
            end do
         end do
c
c     Mean (scalar-pressure)_gradient for ui-theta (67-68)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,67)=pxysth2(x,y,65)-pxysth2(x,y,69)
               pxysth2(x,y,68)=pxysth2(x,y,66)-pxysth2(x,y,70)
            end do
         end do
c
c     Total pressure-scalar_gradient for ui-theta (71-73)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,70+i)=pxysth(x,y,12+i)
               end do
            end do
         end do
c
c     Take away mean component for pressure-scalar_gradient
c
         do i=1,2
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,12+i)=pxysth(x,y,12+i)-
     &                 pxys(x,y,31)*pxysth2(x,y,i)
               end do
            end do
         end do
c
c     Mean pressure-scalar_gradient for ui-theta (74-76)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,73+i)=pxysth2(x,y,70+i)-pxysth(x,y,12+i)
               end do
            end do
         end do
c
c     Total scalar-pressure_gradient for ui-theta (77-79)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,77)=pxysth2(x,y,65)+pxysth2(x,y,71)
               pxysth2(x,y,78)=pxysth2(x,y,66)+pxysth2(x,y,72)
               pxysth2(x,y,79)=             0.+pxysth2(x,y,73)
            end do
         end do
c
c     Mean scalar-pressure_gradient for ui-theta (80-82)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,80)=-pxysth(x,y,1)*pxys2(x,y,115)
               pxysth2(x,y,81)=-pxysth(x,y,1)*pxys2(x,y,116)
               pxysth2(x,y,82)=0.
cccccccccccccccccccccccccccccccccccccccccc
c  use this gives exact zero residual
c                     pxysth2(x,y,80)=pxysth2(x,y,67)+pxysth2(x,y,74)
c                     pxysth2(x,y,81)=pxysth2(x,y,68)+pxysth2(x,y,75)
cccccccccccccccccccccccccccccccccccccccccc
            end do
         end do
c
c     Fluctuating scalar-pressure_gradient for ui-theta (83-85)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,82+i)=pxysth2(x,y,76+i)-pxysth2(x,y,79+i)
               end do
            end do
         end do
c
c     Total dissipation for ui-theta (86-88)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,85+i)=-(1/re+1/(re*pr(scalarind)))
     &                 *pxysth(x,y,30+i)
               end do
            end do
         end do
c
c     Mean dissipation for ui-theta (89-91)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,88+i)=-(1/re+1/(re*pr(scalarind)))
     &                 *(pxys2(x,y,13+i)*pxysth2(x,y,1)+pxys2(x,y,16+i)
     &                 *pxysth2(x,y,2))
               end do
            end do
         end do
c
c     Fluctuating dissipation for ui-theta (92-94)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,91+i)=pxysth2(x,y,85+i)-pxysth2(x,y,88+i)
               end do
            end do
         end do
c
c     Take away the mean component for the scalar dissipation
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,30+i)=pxysth(x,y,30+i)-(pxys2(x,y,13+i)
     &                 *pxysth2(x,y,1)+pxys2(x,y,16+i)*pxysth2(x,y,2))
               end do
            end do
         end do
      end if

c     Transform scalar from T to 1-T
      if (scalar.gt.0) then
         write(*,*) 'Transform scalar from T to 1-T ?'
         write(*,*) '1 for yes, 0 for no'
         read(*,*) i
         if (i.eq.0) goto 9000
         if (i.eq.1) continue
      end if
c
c     Calaulate the mean scalar and the derivatives
c
      do y=1,nyp
         do x=1,nx
            pxysth(x,y,1)=1-pxysth(x,y,1)
            pxysth2(x,y,1)=-pxysth2(x,y,1)
            pxysth2(x,y,2)=-pxysth2(x,y,2)
         end do
      end do
c
c     Total convection for scalar squared
c
      if (scalar.gt.0.and.nxysth.ge.11) then
         do i=3,7
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=totxysth(x,nyp+1-y,i,scalarind)
               end do
            end do
         end do

         do y=1,nyp
            do x=1,nx
               wxys(x,y,6)=wxys(x,y,6)+pxys(x,y,1)-2.*wxys(x,y,3)
               wxys(x,y,7)=wxys(x,y,7)+pxys(x,y,2)-2.*wxys(x,y,4)
            end do
         end do
         call ddx(wxys(1,1,1),wxys(1,1,6),xl)
         call ddy(wxys(1,1,2),wxys(1,1,7),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,26)=(wxys(x,y,1)+wxys(x,y,2))/2.
            end do
         end do
      end if
c
c     Total convection for ui-theta (38-40)
c
      if (scalar.gt.0.and.nxysth.ge.33) then
         do i=4,6
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=totxys(x,nyp+1-y,i)
               end do
            end do
         end do
         do i=13,15
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=totxys(x,nyp+1-y,i)
               end do
            end do
         end do
         call ddx(wxys(1,1,1),wxys(1,1,4),xl)
         call ddy(wxys(1,1,2),wxys(1,1,13),2./dstar)
         call ddx(wxys(1,1,3),wxys(1,1,13),xl)
         call ddy(wxys(1,1,4),wxys(1,1,5),2./dstar)
         call ddx(wxys(1,1,5),wxys(1,1,14),xl)
         call ddy(wxys(1,1,6),wxys(1,1,15),2./dstar)
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,37+i)=-pxysth2(x,y,37+i)
     &                 +wxys(x,y,2*i-1)+wxys(x,y,2*i)
               end do
            end do
         end do
      end if
c
c     Fluctuating 3rd and 4th order scalar central momentum
c
      if (scalar.gt.0.and.nxysth.ge.35) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,35)=pxysth(x,y,35)
            end do
         end do
      end if
      if (scalar.gt.0.and.nxysth.ge.34) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,34)=-pxysth(x,y,34)
            end do
         end do
      end if
c
c     Fluctuating velocity^2-scalar
c
      if (scalar.gt.0.and.nxysth.ge.21) then
         do ivar=16,21
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=-pxysth(x,y,ivar)
               end do
            end do
         end do
      end if
c
c     Fluctuating velocity-scalar^2
c
      if (scalar.gt.0.and.nxysth.ge.8) then
         do ivar=6,8
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=pxysth(x,y,ivar)
               end do
            end do
         end do
      end if
c
c     The rms of scalar field
c
      if (scalar.gt.0.and.nxysth.ge.2) then
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,2)=pxysth(x,y,2)
            end do
         end do
      end if
c
c     Fluctuating scalar-velocity
c
      if (scalar.gt.0.and.nxysth.ge.5) then
         do ivar=3,5
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=-pxysth(x,y,ivar)
               end do
            end do
         end do
c
c     Calculate the derivatives of scalar-velocity
c
         do i=1,3
            call ddx(pxysth2(1,1,1+2*i),pxysth(1,1,2+i),xl)
            call ddy(pxysth2(1,1,2+2*i),pxysth(1,1,2+i),2./dstar)
         end do
      end if
c
c     Fluctuating velocity-scalar_gradient
c
      if (scalar.gt.0.and.nxysth.ge.30) then
         do ivar=22,30
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,ivar)=-pxysth(x,y,ivar)
               end do
            end do
         end do
c
c     Fluctuating scalar-velocity_gradient
c
         do i=1,3
            do j=1,2
               do y=1,nyp
                  do x=1,nx
                     pxysth2(x,y,8+(i-1)*3+j)=
     &                    -pxysth2(x,y,8+(i-1)*3+j)
                  end do
               end do
            end do
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,11+(i-1)*3)=-pxysth2(x,y,11+(i-1)*3)
               end do
            end do
         end do
      end if
c
c     Skewness and flatness for the scalar
c
      if (scalar.gt.0.and.nxysth.ge.34) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,18)=-pxysth2(x,y,18)
            end do
         end do
      end if
      if (scalar.gt.0.and.nxysth.ge.35) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,19)=pxysth2(x,y,19)
            end do
         end do
      end if
c
c     Turbulent Pr
c
      if (scalar.gt.0.and.nxysth.ge.5) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,20)=pxys(x,y,13)*pxysth2(x,y,2)
     &              /(pxysth(x,y,4)*pxys2(x,y,17))
c     pxysth2(x,y,20)=pxysth2(x,y,20)
            end do
         end do
      end if
c
c     Turbulent kinetic energy for scalar
c
      if (scalar.gt.0.and.nxysth.ge.2) then
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,21)=pxysth2(x,y,21)
            end do
         end do
c
c     (k_theta)x,(k_theta)y,(k_theta)xx,(k_theta)yy
c
         do ivar=22,26
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,ivar)=pxysth2(x,y,ivar)
               end do
            end do
         end do
      end if
c
c     Mean convection for scalar squared
c
      if (scalar.gt.0.and.nxysth.ge.11) then
         do i=1,2
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxysth(x,y,1)**2*pxys(x,y,i)
               end do
            end do
         end do
         call ddx(wxys(1,1,3),wxys(1,1,1),xl)
         call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,27)=(wxys(x,y,3)+wxys(x,y,4))/2.
            end do
         end do
c
c     Fluctuating convection for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,28)=pxysth2(x,y,28)
            end do
         end do
c
c     Mean molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxysth(x,y,1)**2
            end do
         end do
         call ddx(wxys(1,1,2),wxys(1,1,1),xl)
         call ddx(wxys(1,1,3),wxys(1,1,2),xl)
         call ddy(wxys(1,1,4),wxys(1,1,1),2./dstar)
         call ddy(wxys(1,1,5),wxys(1,1,4),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,30)=(wxys(x,y,3)+wxys(x,y,5))
     &              /(2.*re*pr(scalarind))
            end do
         end do
c
c     Fluctuating molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,31)=pxysth2(x,y,31)
            end do
         end do
c
c     Total molecular diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,29)=pxysth2(x,y,30)+pxysth2(x,y,31)
            end do
         end do
c
c     Mean turbulent diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxysth(x,y,1)*pxysth(x,y,3)
               wxys(x,y,2)=pxysth(x,y,1)*pxysth(x,y,4)
            end do
         end do
         call ddx(wxys(1,1,3),wxys(1,1,1),xl)
         call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,32)=-(wxys(x,y,3)+wxys(x,y,4))
            end do
         end do
c
c     Fluctuating turbulent diffusion for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,33)=pxysth2(x,y,33)
            end do
         end do
c
c     Production for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,34)=pxysth2(x,y,34)
            end do
         end do
c
c     Total dissipation for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,35)=pxysth2(x,y,35)
            end do
         end do
c
c     Fluctuating dissipation for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,37)=pxysth2(x,y,37)
            end do
         end do
c
c     Mean dissipation for scalar squared
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,36)=pxysth2(x,y,36)
            end do
         end do
      end if
c
c     Mean convection for ui-theta (41-43)
c
      if (scalar.gt.0.and.nxysth.ge.33) then
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,i)=pxys(x,y,i)*pxysth(x,y,1)
               end do
            end do
            call ddx(wxys(1,1,3+i),wxys(1,1,i),xl)
            call ddy(wxys(1,1,6+i),wxys(1,1,i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,40+i)=wxys(x,y,3+i)*pxys(x,y,1)+
     &                 wxys(x,y,i+6)*pxys(x,y,2)
               end do
            end do
         end do
c
c     Fluctuating convection for ui-theta (44-46)
c
         do ivar=44,46
            do i=1,2
               do y=1,nyp
                  do x=1,nx
                     pxysth2(x,y,ivar)=-pxysth2(x,y,ivar)
                  end do
               end do
            end do
         end do
c
c     Mean molecular diffusion for ui-theta (50-52)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  wxys(x,y,4*i-3)=pxysth2(x,y,1)*pxys(x,y,i)
                  wxys(x,y,4*i-2)=pxysth2(x,y,2)*pxys(x,y,i)
                  wxys(x,y,4*i-1)=pxys2(x,y,13+i)*pxysth(x,y,1)
                  wxys(x,y,4*i)=pxys2(x,y,16+i)*pxysth(x,y,1)
               end do
            end do
            call ddx(wxys(1,1,4*i+9),wxys(1,1,4*i-3),xl)
            call ddy(wxys(1,1,4*i+10),wxys(1,1,4*i-2),2./dstar)
            call ddx(wxys(1,1,4*i+11),wxys(1,1,4*i-1),xl)
            call ddy(wxys(1,1,4*i+12),wxys(1,1,4*i),2./dstar)
         end do
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,50)=(wxys(x,y,13)+wxys(x,y,14))
     &              /(re*pr(scalarind))+(wxys(x,y,15)
     &              +wxys(x,y,16))/re
               pxysth2(x,y,51)=(wxys(x,y,17)+wxys(x,y,18))
     &              /(re*pr(scalarind))+(wxys(x,y,19)
     &              +wxys(x,y,20))/re
               pxysth2(x,y,52)=(wxys(x,y,21)+wxys(x,y,22))
     &              /(re*pr(scalarind))+(wxys(x,y,23)
     &              +wxys(x,y,24))/re
            end do
         end do
c
c     Fluctuating molecular diffusion for ui-theta (53-55)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,52+i)=-pxysth2(x,y,52+i)
               end do
            end do
         end do
c
c     Total molecular diffusion for ui-theta (47-49)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,46+i)=pxysth2(x,y,49+i)
     &                 +pxysth2(x,y,52+i)
               end do
            end do
         end do
c
c     Mean turbulent diffusion for ui-theta (56-58)
c
         do y=1,nyp
            do x=1,nx
               wxys(x,y,1)=pxys(x,y,1)*pxysth(x,y,3)+
     &              pxys(x,y,4)**2*pxysth(x,y,1)
               wxys(x,y,2)=pxys(x,y,1)*pxysth(x,y,4)+
     &              pxys(x,y,13)*pxysth(x,y,1)
               wxys(x,y,3)=pxys(x,y,2)*pxysth(x,y,3)+
     &              pxys(x,y,13)*pxysth(x,y,1)
               wxys(x,y,4)=pxys(x,y,2)*pxysth(x,y,4)+
     &              pxys(x,y,5)**2*pxysth(x,y,1)
               wxys(x,y,5)=pxys(x,y,3)*pxysth(x,y,3)+
     &              pxys(x,y,14)*pxysth(x,y,1)
               wxys(x,y,6)=pxys(x,y,3)*pxysth(x,y,4)+
     &              pxys(x,y,15)*pxysth(x,y,1)
            end do
         end do
         do i=1,3
            call ddx(wxys(1,1,2*i-1),wxys(1,1,2*i-1),xl)
            call ddy(wxys(1,1,2*i),wxys(1,1,2*i),2./dstar)
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,55+i)=
     &                 -(wxys(x,y,2*i-1)+wxys(x,y,2*i))
               end do
            end do
         end do
c
c     Fluctuating turbulent diffusion for ui-theta (59-61)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,58+i)=-pxysth2(x,y,58+i)
               end do
            end do
         end do
c
c     Production for ui-theta (62-64)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,62)=-(pxysth(x,y,3)*pxys2(x,y,14)
     &              +pxysth(x,y,4)*pxys2(x,y,17)+pxysth2(x,y,1)
     &              *pxys(x,y,4)**2+pxysth2(x,y,2)*pxys(x,y,13))
               pxysth2(x,y,63)=-(pxysth(x,y,3)*pxys2(x,y,15)
     &              +pxysth(x,y,4)*pxys2(x,y,18)+pxysth2(x,y,1)
     &              *pxys(x,y,13)+pxysth2(x,y,2)*pxys(x,y,5)**2)
               pxysth2(x,y,64)=-(pxysth(x,y,3)*pxys2(x,y,16)
     &              +pxysth(x,y,4)*pxys2(x,y,19)+pxysth2(x,y,1)
     &              *pxys(x,y,14)+pxysth2(x,y,2)*pxys(x,y,15))
            end do
         end do
c
c     Total (scalar-pressure)_gradient for ui-theta (65-66) (<pth>z is zreo)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,65)=-pxys2(x,y,115)-pxysth2(x,y,65)
               pxysth2(x,y,66)=-pxys2(x,y,116)-pxysth2(x,y,66)
            end do
         end do
c
c     Fluctuating pressure-scalar
c
         do y=1,nyp
            do x=1,nx
               pxysth(x,y,12)=-pxysth(x,y,12)
            end do
         end do
c
c     Fluctuating (scalar-pressure)_gradient for ui-theta (69-70)
c
         call ddx(pxysth2(1,1,69),pxysth(1,1,12),xl)
         call ddy(pxysth2(1,1,70),pxysth(1,1,12),2./dstar)
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,69)=-pxysth2(x,y,69)
               pxysth2(x,y,70)=-pxysth2(x,y,70)
            end do
         end do
c
c     Mean (scalar-pressure)_gradient for ui-theta (67-68)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,67)=pxysth2(x,y,65)-pxysth2(x,y,69)
               pxysth2(x,y,68)=pxysth2(x,y,66)-pxysth2(x,y,70)
            end do
         end do
c
c     Total pressure-scalar_gradient for ui-theta (71-73)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,70+i)=-pxysth2(x,y,70+i)
               end do
            end do
         end do
c
c     Fluctuaing pressure-scalar_gradient
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,12+i)=-pxysth(x,y,12+i)
               end do
            end do
         end do
c
c     Mean pressure-scalar_gradient for ui-theta (74-76)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,73+i)=pxysth2(x,y,70+i)
     &                 -pxysth(x,y,12+i)
               end do
            end do
         end do
c
c     Total scalar-pressure_gradient for ui-theta (77-79)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,77)=pxysth2(x,y,65)+pxysth2(x,y,71)
               pxysth2(x,y,78)=pxysth2(x,y,66)+pxysth2(x,y,72)
               pxysth2(x,y,79)=             0.+pxysth2(x,y,73)
            end do
         end do
c
c     Mean scalar-pressure_gradient for ui-theta (80-82) (82 is zero)
c
         do y=1,nyp
            do x=1,nx
               pxysth2(x,y,80)=-pxysth(x,y,1)*pxys2(x,y,115)
               pxysth2(x,y,81)=-pxysth(x,y,1)*pxys2(x,y,116)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     use this gives exact zero residual
c               pxysth2(x,y,80)=pxysth2(x,y,67)+pxysth2(x,y,74)
c               pxysth2(x,y,81)=pxysth2(x,y,68)+pxysth2(x,y,75)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            end do
         end do
c
c     Fluctuating scalar-pressure_gradient for ui-theta (83-85)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,82+i)=-pxysth2(x,y,82+i)
               end do
            end do
         end do
c
c     Total dissipation for ui-theta (86-88)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,85+i)=-pxysth2(x,y,85+i)
               end do
            end do
         end do
c
c     Mean dissipation for ui-theta (89-91)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,88+i)=-pxysth2(x,y,88+i)
               end do
            end do
         end do
c
c     Fluctuating dissipation for ui-theta (92-94)
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth2(x,y,91+i)=-pxysth2(x,y,91+i)
               end do
            end do
         end do
c
c     The  scalar dissipation
c
         do i=1,3
            do y=1,nyp
               do x=1,nx
                  pxysth(x,y,30+i)=-pxysth(x,y,30+i)
               end do
            end do
         end do
      end if
 9000 continue
c
c     Dissipation length scale Lk
c
      do  y=1,nyp
         do  x=1,nx
            pxys2(x,y,67)=pxys2(x,y,1)**1.5 / pxys2(x,y,26)
         end do
      end do
c
c     Dissipation length scale Lu
c
      do  y=1,nyp
         do  x=1,nx
            pxys2(x,y,68)=pxys(x,y,4)**3. / pxys2(x,y,20) / 2.
         end do
      end do
c
c     Ratio Lu / Lk
c
      do  y=1,nyp
         do  x=1,nx
            pxys2(x,y,69) = pxys2(x,y,68) / pxys2(x,y,67)
         end do
      end do
c
c     Ratio Lu / Lk (norm)
c
      do  y=1,nyp
         do  x=1,nx
            pxys2(x,y,70) = pxys2(x,y,69) *sqrt(1.5)
         end do
      end do
c
c     Tu (local)
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,73) = sqrt(2.*pxys2(x,y,1)/3./
     &           (pxys(x,y,1)**2.+pxys(x,y,2)**2.+pxys(x,y,3)**2.))
         end do
      end do
c
c     Tu (freestream)
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,74) = sqrt(2.*pxys2(x,y,1)/3.)
         end do
      end do
c
c     isotropy
c
      do y=1,nyp
         do x=1,nx
            rmsmax=-10.e10
            rmsmin=10.e10
            do i=4,6
               rmsmax=max(rmsmax,pxys(x,y,i))
               rmsmin=min(rmsmin,pxys(x,y,i))
            end do
            pxys2(x,y,79) = rmsmax-rmsmin
            pxys2(x,y,80) = pxys2(x,y,79)/pxys2(x,y,74)
         end do
      end do
c
c     99% boundary layer thickness
c
      do x=1,nx
         do y=1,nyp
            if (pxys(x,y,1).gt.0.99) then
               goto 1234
            end if
         end do

 1234    continue
c
c     Do linear interpolation
c
         rmsmax = 2./dstar - eta(y-1)
         rmsmin = 2./dstar - eta(y)

         pxys2(x,1,81) = 1./(pxys(x,y,1)-pxys(x,y-1,1))*
     &        (rmsmin*(0.99-pxys(x,y-1,1)) +
     &        rmsmax*(pxys(x,y,1)-0.99))
c
c     Original result
c         pxys2(x,1,81)=2./dstar-eta(y)
      end do
c
c     99% scalar boundary layer thickness
c
c      if (scalar.gt.0) then
c         do x=1,nx
c            do y=1,nyp
c               if (pxysth(x,y,1).gt.0.99) then
c                  goto 2345
c               end if
c            end do

c 2345       continue
c     do linear interpolation
c            rmsmax = 2./dstar - eta(y-1)
c            rmsmin = 2./dstar - eta(y)
c            pxys2(x,1,81) = 1./(pxysth(x,y,1)-pxysth(x,y-1,1))*
c     &           (rmsmin*(0.99-pxysth(x,y-1,1)) +
c     &           rmsmax*(pxysth(x,y,1)-0.99))
c         end do
c      end if
c
c     Smagorinsky constant
c
      do x=1,nx
         do y=1,nyp
            pxys2(x,y,85) = sqrt(max(0.,pxys(x,y,52)))
         end do
      end do
c
c     max urms
c     Note: the first maximum off the wall is considered
c
c_AM      choose if you want the first maximum or just the maximum.
      firstmax=0
      if (firstmax.eq.1) then
         write(*,*)'WARNING FIRST MAX IS 1!!!! hardcoded. It will
     &calculate global maximum (not first)'
      end if
      do x=1,nx
         pxys2(x,1,71) = 0.
         i = 0
         do y=1,nyp
            if (firstmax.eq.0) then
               if (i.eq.0.and.pxys(x,y,4).gt.pxys2(x,1,71)) then
                  pxys2(x,1,71) = pxys(x,y,4)
                  pxys2(x,1,72)   = 2./dstar - eta(y)
               end if
            elseif (firstmax.eq.1) then
               if (pxys(x,y,4).gt.pxys2(x,1,71)) then
                  pxys2(x,1,71) = pxys(x,y,4)
                  pxys2(x,1,72)   = 2./dstar - eta(y)
               end if
            end if
c     if ((2./dstar-eta(y)).gt.10.) then
c     pxys2(x,1,72) = 0.
c     else
c     pxys2(x,1,72)   = 2./dstar - eta(y)
c     end if
            if (y.gt.1) then
               if (pxys(x,y,4).gt.pxys(x,y+1,4)) then
                  i=1
               end if
            end if
         end do
         do y=2,nyp
            pxys2(x,y,71) = pxys2(x,1,71)
            pxys2(x,y,72) = pxys2(x,1,72)
            pxys2(x,y,81) = pxys2(x,1,81)
         end do
      end do
c
c     Blasius profile
c
      etamax = 20.
      deta = etamax / real(n)
c
c     eta similarity coordinates
c
      do i=0, n
        etab(i) = etamax*real(i)/real(n)
      end do
c
c     Compute Blasius solution
c
      call blas(etab,f,fd,fdd,n)

      do x=1,nx
         if (x.lt.nx/2) x1 = x+nx/2
         if (x.gt.nx/2) x1 = x-nx/2
         xx = re/1.7208**2+ (x1-1)*(xl/real(nx-1))
         do y=1,nyp

            yy = 2./dstar-eta(y)

            et = yy / sqrt(xx) *sqrt(re)
            if (et.gt.etamax) then
               pxys2(x,y,75) = 1.
            else
               im = int(et / deta)
               ip = int(et / deta)+1
               etam = im*deta
               etap = ip*deta
               pxys2(x,y,75) =
     &              1./deta*((etap-et)*fd(im)+(et-etam)*fd(ip))
            end if
         end do
      end do

      do x=1,nx
         xx=0.
         yy=0.
         do y=1,nyp
            pxys2(x,y,76) = pxys(x,y,1) - pxys2(x,y,75)
            xx=max(xx,pxys2(x,y,76))
            yy=min(yy,pxys2(x,y,76))
         end do
         if (xx.ne.0..and.yy.ne.0.) then
            do y=1,nyp
               pxys2(x,y,77) = pxys2(x,y,76) / xx
               pxys2(x,y,78) = -pxys2(x,y,76) / yy
            end do
         end if
      end do
c
c     Compute the full Reynolds stresses
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,91)=sqrt(pxys(x,y,4 )**2 + (pxys(x,y,54)))
            pxys2(x,y,92)=sqrt(pxys(x,y,5 )**2 + (pxys(x,y,57)))
            pxys2(x,y,93)=sqrt(pxys(x,y,6 )**2 + (pxys(x,y,59)))
            pxys2(x,y,94)=     pxys(x,y,13)    + pxys(x,y,55)
         end do
      end do
c
c     Deviations from isotropy plus SGS stresses
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,95)=pxys(x,y,4)**2+pxys(x,y,54)-
     &           1./3.*(pxys(x,y,4)**2+pxys(x,y,5)**2+pxys(x,y,6)**2)
            pxys2(x,y,96)=pxys(x,y,5)**2+pxys(x,y,57)-
     &           1./3.*(pxys(x,y,4)**2+pxys(x,y,5)**2+pxys(x,y,6)**2)
            pxys2(x,y,97)=pxys(x,y,6)**2+pxys(x,y,59)-
     &           1./3.*(pxys(x,y,4)**2+pxys(x,y,5)**2+pxys(x,y,6)**2)
         end do
      end do
c
c     Trace of SGS tensor
c
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,98) = pxys(x,y,54)+pxys(x,y,57)+pxys(x,y,59)
         end do
      end do
c
c     Total shear stress
c
      if (nxys2.ge.268) then
         do y=1,nyp
            do x=1,nx
               pxys2(x,y,268) = -pxys2(x,y,94)+pxys2(x,y,17)/re
            end do
         end do
      end if
c
c     MHD statistics
c
c     Take away mean component and take square root for rms
c
      do y=1,nyp
         do x=1,nx
            pxys(x,y,84)=sqrt(max(0.,pxys(x,y,84)-
     &           pxys(x,y,83)**2))

            pxys(x,y,88)=sqrt(max(0.,pxys(x,y,88)-
     &           pxys(x,y,85)**2))
            pxys(x,y,89)=sqrt(max(0.,pxys(x,y,89)-
     &           pxys(x,y,86)**2))
            pxys(x,y,90)=sqrt(max(0.,pxys(x,y,90)-
     &           pxys(x,y,87)**2))

            pxys(x,y,91)=pxys(x,y,91)-
     &           pxys(x,y,85)*pxys(x,y,86)
            pxys(x,y,92)=pxys(x,y,92)-
     &           pxys(x,y,85)*pxys(x,y,87)
            pxys(x,y,93)=pxys(x,y,93)-
     &           pxys(x,y,86)*pxys(x,y,87)

         end do
      end do
c
c     Budget terms for the MHD terms (269-272)
c
      if (nxys2.ge.274) then
         do y=1,nyp
         do x=1,nx
            pxys2(x,y,269)=-mhd_n*(pxys(x,y,85)**2+pxys(x,y,86)**2+
     &           pxys(x,y,87)**2)
            pxys2(x,y,270)=-mhd_n*(pxys(x,y,88)**2+pxys(x,y,89)**2+
     &           pxys(x,y,90)**2)

         end do
      end do

      do y=1,nyp
         do x=1,nx
            wxys(x,y,1)=pxys(x,y,83)*pxys(x,y,85)
            wxys(x,y,2)=pxys(x,y,83)*pxys(x,y,86)
         end do
      end do
      call ddx(wxys(1,1,3),wxys(1,1,1),xl)
      call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,271) = -mhd_n*(wxys(x,y,3)+wxys(x,y,4))
         end do
      end do

      do y=1,nyp
         do x=1,nx
            wxys(x,y,1)=pxys(x,y,94)-pxys(x,y,83)*pxys(x,y,85)
            wxys(x,y,2)=pxys(x,y,95)-pxys(x,y,83)*pxys(x,y,86)
         end do
      end do
      call ddx(wxys(1,1,3),wxys(1,1,1),xl)
      call ddy(wxys(1,1,4),wxys(1,1,2),2./dstar)
      do y=1,nyp
         do x=1,nx
            pxys2(x,y,272) = -mhd_n*(wxys(x,y,3)+wxys(x,y,4))
         end do
      end do


      do y=1,nyp
         do x=1,nx
c
c     Mean budget
c
            pxys2(x,y,273) = -pxys2(x,y,126)-pxys2(x,y,89)+
     &           pxys2(x,y,131)+pxys2(x,y,140)+pxys2(x,y,134)+
     &           pxys2(x,y,83)+pxys2(x,y,269)+pxys2(x,y,271)
c
c     Fluct budget
c
            pxys2(x,y,274) = -pxys2(x,y,127)+pxys2(x,y,89)+
     &           pxys2(x,y,132)+pxys2(x,y,141)+pxys2(x,y,135)+
     &           pxys2(x,y,84)+pxys2(x,y,270)+pxys2(x,y,272)
         end do
      end do
      end if

      end subroutine pvar



      subroutine blas(eta,f,fd,fdd,n)
      implicit none

      integer i,n
      real eta(0:n),f(0:n),fd(0:n),fdd(0:n)
c
c     Boundary conditions
c
      f(0)=0.
      fd(0)=0.
      fdd(0)=0.33205733621545
c
c     Integration
c
      do i=1,n
         call runge(f(i-1),fd(i-1),fdd(i-1),f(i),fd(i),fdd(i),
     $        eta(i)-eta(i-1))
      end do

      end subroutine blas



      subroutine runge(x,y,z,x1,y1,z1,h)

      implicit none

      real h
      real x,y,z
      real x1,y1,z1
      real f1,f2,f3
      real xt1,yt1,zt1
      real xt2,yt2,zt2
      real xt3,yt3,zt3
      real xt4,yt4,zt4

      f1(x,y,z)=y
      f2(x,y,z)=z
      f3(x,y,z)=-x*z /2.

      xt1=h*f1(x,y,z)
      yt1=h*f2(x,y,z)
      zt1=h*f3(x,y,z)

      xt2=h*f1(x+xt1/2.,y+yt1/2.,z+zt1/2.)
      yt2=h*f2(x+xt1/2.,y+yt1/2.,z+zt1/2.)
      zt2=h*f3(x+xt1/2.,y+yt1/2.,z+zt1/2.)

      xt3=h*f1(x+xt2/2.,y+yt2/2.,z+zt2/2.)
      yt3=h*f2(x+xt2/2.,y+yt2/2.,z+zt2/2.)
      zt3=h*f3(x+xt2/2.,y+yt2/2.,z+zt2/2.)

      xt4=h*f1(x+xt3,y+yt3,z+zt3)
      yt4=h*f2(x+xt3,y+yt3,z+zt3)
      zt4=h*f3(x+xt3,y+yt3,z+zt3)

      x1=x+(xt1+2.*xt2+2.*xt3+xt4)/6.
      y1=y+(yt1+2.*yt2+2.*yt3+yt4)/6.
      z1=z+(zt1+2.*zt2+2.*zt3+zt4)/6.

      end subroutine runge
