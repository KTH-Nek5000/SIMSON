c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine fshift(ur,ui,u0low,u0upp,w0low,w0upp,
     &     prex,prez,pres,prea,
     &     u2r,u2i,u3r,u3i,wr,wi,xs,xl,zl)
c
c     Find an appropriate shift velocity to maximize the CFL number
c     (i.e. time step) and reset the mean spanwise and streamwise flow
c     velocity accordingly. We insert shift velocities so that min and
c     max balance around zero this is a crude approximation to
c     minimizing the CFL number
c
c     Here, ushift,wshift are determined and the velocity is adapted
c     accordingly; u0low,u0upp,w0low,w0upp are updated.
c
      implicit none

      include 'par.f'

      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real u3r(nx/2*mbz,nyp,3),u3i(nx/2*mbz,nyp,3)
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real wr(nxp/2+1,mby,nzd,3),wi(nxp/2+1,mby,nzd,3)
      real u0low,u0upp,w0low,w0upp
      real xs,xl,zl

      logical sym
      integer yb,i,y
      real umin,umax,wmin,wmax,ushift,wshift

#ifdef MPI
      write(*,*) 'fshift not implemented with MPI'
      call stopnow(45345)
#endif

      do yb=1,nyp,mby
         do i=1,3
            call getxz(u2r(1,1,1,i),u2i(1,1,1,i),yb,i,1,ur,ui)
c
c     u,v symmetric, w antisymmetric if this is a 'symmetric' case
c
            sym=i.le.2
            call fft2db(u2r(1,1,1,i),u2i(1,1,1,i),sym,min(mby,nyp-yb+1),
     &           prex,prez,pres,prea,wr,wi)
         end do
         call boxext(vext,cext,u2r,u2i,u2r,u2i,xs,yb,xl,zl)
      end do

      umin=1.e20
      umax=-1.e20
      wmin=1.e20
      wmax=-1.e20
      do y=1,nyp
         umin=min(umin,vext(y,1,1))
         umax=max(umax,vext(y,2,1))
         wmin=min(wmin,vext(y,1,3))
         wmax=max(wmax,vext(y,2,3))
      end do
c
c     We insert shift velocities so that min and max balance around zero
c     this is a crude approximation to minimizing the CFL number
c
      ushift=(umin+umax)*0.5
      wshift=(wmin+wmax)*0.5
c
c     For symmetric simulation we can't introduce spanwise shifts
c
      if (nfzsym.eq.1) wshift=0.0
c
c     Do not introdce small shifts in the spanwise direction
c
      if (abs(wshift).lt..2) wshift=0.
c
c     Get the mean flow plane
c
      do i=1,3
         call getxy(u3r(1,1,i),u3i(1,1,i),1,i,ur,ui)
      end do
c
c     Reset the mean velocity
c
      do y=1,nyp
         u3r(1,y,1)=u3r(1,y,1)-ushift
         u3r(1,y,3)=u3r(1,y,3)-wshift
      end do
      u0low=u3r(1,nyp,1)
      u0upp=u3r(1,1,1)
      w0low=u3r(1,nyp,3)
      w0upp=u3r(1,1,3)
c
c     Put it back
c
      do i=1,3
         call putxy(u3r(1,1,i),u3i(1,1,i),1,i,ur,ui)
      end do

      end subroutine fshift
