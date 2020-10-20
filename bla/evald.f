c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine evald(dj1ur,dj1ui,dj1uh1,dj1uh2,dj1uh3,
     &     dj1omr,dj1omi,dj1omh,
     &     d2vr,d2vi,domyr,domyi,
     &     d2omyr,d2omyi,d2vh,d2vh2,domyh,d2omyh,
     &     bbeta,k2i,alfa,ibc)
c
c     Evaluates the derivatives of the
c     particular and homogeneous solutions at y=1
c     with the nonlinear terms calculated in step 2 as a driving force
c
c     Note that the use of ibc is only to avoid unnecessary
c     parts of the code. The code should work if all
c     lines are executed regardless of the value of ibc
c
      implicit none

      include 'par.f'
c
c     Evaluated derivatives at the boundary y=1
c
      real dj1ur(nx/2*mbz,3,5),dj1ui(nx/2*mbz,3,5)
      real dj1uh1(nx/2*mbz,3,5),dj1uh2(nx/2*mbz,3,5)
      real dj1uh3(nx/2*mbz,3,4)
      real dj1omr(nx/2*mbz,5),dj1omi(nx/2*mbz,5)
      real dj1omh(nx/2*mbz,5)
c
c     Field variables
c
      real d2vr(nx/2*mbz,nyp),d2vi(nx/2*mbz,nyp)
      real domyr(nx/2*mbz,nyp),domyi(nx/2*mbz,nyp)
      real d2omyr(nx/2*mbz,nyp),d2omyi(nx/2*mbz,nyp)
      real d2vh(nx/2*mbz,nyp)
      real d2vh2(nx/2*mbz,nyp)
      real domyh(nx/2*mbz,nyp),d2omyh(nx/2*mbz,nyp)
c
c     Wavenumbers
c
      real bbeta(nx/2*mbz),k2i(nx/2*mbz),alfa(nx/2*mbz)

      integer ibc
c
c     Local variables
c
      integer j,xz,nxz,y
      real fy,ffy

      nxz=nx/2*mbz
c
c     We make use of the following properties of the variables at y=1 :
c
c     vr=vi=dvr=dvi=0
c     vh=0,dvh=2
c     vh2=2,dvh2=0
c     omyr=omyi=0, omyh=2
c
c     As a bc :
c     dj1q(:,i,j+1) is the jth derivative of the q
c     where q is one of vr,vi,vh1,vh2,omr,omi or omh
c     calculate the derivatives of each cheb poly
c     at x=1 :
c
      do j=0,4
         do xz=1,nxz
            dj1ur(xz,2,j+1)=0.0
            dj1ui(xz,2,j+1)=0.0
            dj1uh1(xz,2,j+1)=0.0
            dj1uh2(xz,2,j+1)=0.0
            dj1omr(xz,j+1)=0.0
            dj1omi(xz,j+1)=0.0
            dj1omh(xz,j+1)=0.0
         end do
      end do
      do y=1,ny
         do xz=1,nxz
            dj1ur(xz,2,3)=dj1ur(xz,2,3)+d2vr(xz,y)
            dj1ui(xz,2,3)=dj1ui(xz,2,3)+d2vi(xz,y)
            dj1uh1(xz,2,3)=dj1uh1(xz,2,3)+d2vh(xz,y)
            dj1uh2(xz,2,3)=dj1uh2(xz,2,3)+d2vh2(xz,y)
         end do
      end do
      do xz=1,nxz
        dj1omh(xz,1)=2.0
        dj1uh2(xz,2,1)=2.0
        dj1uh1(xz,2,2)=2.0
      end do
      if (ibc.ne.20.and.ibc.ne.120) then
         do y=1,ny
            fy=real(y-1)
            fy=fy*fy
            ffy=fy*(fy-1.)/3.
            do xz=1,nxz
c
c     d(Tk(y))dy (y=1) = k*k
c
               dj1ur(xz,2,4)=dj1ur(xz,2,4)+d2vr(xz,y)*fy
               dj1ui(xz,2,4)=dj1ui(xz,2,4)+d2vi(xz,y)*fy
               dj1uh1(xz,2,4)=dj1uh1(xz,2,4)+d2vh(xz,y)*fy
               dj1uh2(xz,2,4)=dj1uh2(xz,2,4)+d2vh2(xz,y)*fy
c
c     d2(Tk(y))dy2 (y=1) = 1/3*(k*k)*(k*k-1)
c
               dj1ur(xz,2,5)=dj1ur(xz,2,5)+d2vr(xz,y)*ffy
               dj1ui(xz,2,5)=dj1ui(xz,2,5)+d2vi(xz,y)*ffy
               dj1uh1(xz,2,5)=dj1uh1(xz,2,5)+d2vh(xz,y)*ffy
               dj1uh2(xz,2,5)=dj1uh2(xz,2,5)+d2vh2(xz,y)*ffy
               dj1omr(xz,2)=dj1omr(xz,2)+domyr(xz,y)
               dj1omi(xz,2)=dj1omi(xz,2)+domyi(xz,y)
               dj1omh(xz,2)=dj1omh(xz,2)+domyh(xz,y)
               dj1omr(xz,3)=dj1omr(xz,3)+d2omyr(xz,y)
               dj1omi(xz,3)=dj1omi(xz,3)+d2omyi(xz,y)
               dj1omh(xz,3)=dj1omh(xz,3)+d2omyh(xz,y)
               dj1omr(xz,4)=dj1omr(xz,4)+d2omyr(xz,y)*fy
               dj1omi(xz,4)=dj1omi(xz,4)+d2omyi(xz,y)*fy
               dj1omh(xz,4)=dj1omh(xz,4)+d2omyh(xz,y)*fy
            end do
         end do
      end if
c
c     Now generate the derivatives on u and w
c     u=1/k2(i*alfa*dv-i*beta*omy)
c     w=1/k2(i*beta*dv+i*alfa*omy)
c
      do j=1,4
         do xz=1,nxz
            dj1ur(xz,1,j)=
     &           (-alfa(xz)*dj1ui(xz,2,j+1)+
     &           bbeta(xz)*dj1omi(xz,j))*k2i(xz)
            dj1ui(xz,1,j)=
     &           (alfa(xz)*dj1ur(xz,2,j+1)-
     &           bbeta(xz)*dj1omr(xz,j))*k2i(xz)
            dj1ur(xz,3,j)=
     &           (-bbeta(xz)*dj1ui(xz,2,j+1)-
     &           alfa(xz)*dj1omi(xz,j))*k2i(xz)
            dj1ui(xz,3,j)=
     &           (bbeta(xz)*dj1ur(xz,2,j+1)+
     &           alfa(xz)*dj1omr(xz,j))*k2i(xz)
            dj1uh1(xz,1,j)=alfa(xz)*dj1uh1(xz,2,j+1)*k2i(xz)
            dj1uh1(xz,3,j)=bbeta(xz)*dj1uh1(xz,2,j+1)*k2i(xz)
            dj1uh2(xz,1,j)=alfa(xz)*dj1uh2(xz,2,j+1)*k2i(xz)
            dj1uh2(xz,3,j)=bbeta(xz)*dj1uh2(xz,2,j+1)*k2i(xz)
            dj1uh3(xz,1,j)=-bbeta(xz)*dj1omh(xz,j)*k2i(xz)
            dj1uh3(xz,2,j)=0.0
            dj1uh3(xz,3,j)=alfa(xz)*dj1omh(xz,j)*k2i(xz)
         end do
      end do

      end subroutine evald
