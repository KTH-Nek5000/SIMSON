c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine linearp(ur,ui,zb,zbp,re,alfa,beta,cim,prey,
     &     h3r,h3i,u3r,u3i,om3r,om3i,
     &     f,fvr,fvi,
     &     bis,d2phr,d2phi,phr,phi,
     &     vr,vi,d2vr,d2vi,
     &     hvr,hvi,
     &     q,c,d,e,w3,rv1r,rv1i,wd1,it,dtn)
c
c     Linear part of the pressure solver (see section 4.6. in the
c     manual). Implements the solution to the Poisson/Helmholtz equation
c     with appropriate boundary conditions. Addtionally, recomputes
c     the vorticities which were overwritten in nonlinp.
c
c     Only Neumann boundary conditions are supported right now.
c
c     f is mapped onto fvr/fvi
c     d2phr,d2phr is mapped onto phr,phi and onto bis(1,1,2)
c
      implicit none

      include 'par.f'

      integer nxz
      parameter (nxz=nx/2*mbz)

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      logical cim
      integer neumann
      integer zb,it
      real re,dtn
      real alfa(nxz),beta(nz)
      real prey(nyp*2+15)
      real h3r(nxz,nyp,3),h3i(nxz,nyp,3)
      real u3r(nxz,nyp,3),u3i(nxz,nyp,3)
      real om3r(nxz,nyp,3),om3i(nxz,nyp,3)
      real d2phr(nxz,nyp),d2phi(nxz,nyp)
      real phr(nxz,nyp),phi(nxz,nyp)
      real d2vr(nxz,nyp),d2vi(nxz,nyp)
      real vr(nxz,nyp),vi(nxz,nyp)
      real fvr(nxz,nyp),fvi(nxz,nyp)
      real hvr(nxz,nyp),hvi(nxz,nyp)
      real w3(nx/2,mbz,nyp)
      real q(nxz,ny+2),c(nxz,ny),d(nxz,ny),e(nxz,ny)
      real f(nxz,nyp,5),bis(nxz,nyp,6)
      real rv1r(nx/2,nzc,2),rv1i(nx/2,nzc,2),wd1(nyp,4)
c
c     Local variables
c
      integer x,z,i,xz,y,pl
      real bbeta(nxz),k2(nxz),bc(nxz,2,2),c4
c
c     Evaluated derivatives at the boundary y=1
c
      real temp(nyp),temi(nyp),bc0
      real cc
c
c     MPI
c
      integer zbp
c
c-----[----------------------------------------------------------------]
c
c
c     Don't do it for the (spanwise) oddball mode
c
      cc = 0.
      if (mbz.eq.1.and.zb.eq.(nz+1)/2+1) return
c
c     Set Neumann boundary conditions. Dirichlet is not really implemented.
c
      neumann = 1
c
c     Get wavenumbers
c
      do z=1,mbz
         do x=1,nx/2
            xz=x+nx/2*(z-1)
            bbeta(xz)=beta(z+zb-1)
            k2(xz)=alfa(xz)*alfa(xz)+bbeta(xz)*bbeta(xz)
         end do
      end do
c
c     Get an xy box and chebyshev transform for RHS
c     These terms contain d/dx H1 + d/dz H3 and H2
c
      call getxy(h3r(1,1,1),h3i(1,1,1),zbp,4,ur,ui)
      call getxy(h3r(1,1,2),h3i(1,1,2),zbp,5,ur,ui)

      if (neumann.eq.1) then
c
c     Save vr/vi=H_2 at wall and free stream, and normalize.
c     This will later be used in the boundary conditions.
c
         c4=1./(real(nxp)*real(nzp))
         do xz=1,nxz
            vr(xz,1)=h3r(xz,1,2)*c4
            vr(xz,2)=h3r(xz,nyp,2)*c4
            vi(xz,1)=h3i(xz,1,2)*c4
            vi(xz,2)=h3i(xz,nyp,2)*c4
         end do
      else
         call stopnow(754242)
      end if
c
c     Go into Chebyshev space
c
      c4=2./(real(nxp)*real(nyp-1)*real(nzp))
      do i=1,2
         call vchbf(h3r(1,1,i),w3,nyp,nxz,nxz,1,prey)
         call vchbf(h3i(1,1,i),w3,nyp,nxz,nxz,1,prey)
         do y=1,nyp
            do xz=1,nxz
               h3r(xz,y,i)=h3r(xz,y,i)*c4
               h3i(xz,y,i)=h3i(xz,y,i)*c4
            end do
         end do
      end do
c
c     Get an xy box and chebyshev transform for E
c
      call getxy(hvr,hvi,zbp,8,ur,ui)
c
c     Save complete H_2 and energy E at wall for
c     alpha=beta=0.(D(p+E)=H_2, wavenumber zero)
c
      if (zb.eq.1) then
         do y=1,nyp
            temp(y)=h3r(1,y,2)
         end do
         bc0=hvr(1,nyp)/(2.*real(nxp)*real(nzp))
      end if
c
c     Go with E into Chebyshev space and normalise
c
      call vchbf(hvr,w3,nyp,nxz,nxz,1,prey)
      call vchbf(hvi,w3,nyp,nxz,nxz,1,prey)
      c4=2./(2.*real(nxp)*real(nyp-1)*real(nzp))
      do y=1,nyp
         do xz=1,nxz
            hvr(xz,y)=hvr(xz,y)*c4
            hvi(xz,y)=hvi(xz,y)*c4
         end do
      end do
c
c     Calculate fvr, fvi (=dH2/dy)
c
      call dcheb(fvr,h3r(1,1,2),nyp,nxz,nxz)
      call dcheb(fvi,h3i(1,1,2),nyp,nxz,nxz)
c
c     Solve the Poisson equation for p+E.(E=.5*(u**2+v**2+w**2))
c
c     First, compute dH_i / dx_i, i.e. the RHS
c
      do y=1,nyp
         do xz=1,nxz
            fvr(xz,y)=fvr(xz,y)+h3r(xz,y,1)
            fvi(xz,y)=fvi(xz,y)+h3i(xz,y,1)
         end do
      end do
c
c     Compute the boundary conditions at wall and free stream.
c     Dirichlet bc: k^2*(p+E)=-i*alpha*(H_1+(d2-k^2)u/Re-du/dt)-i*beta*(H_3+
c     (d2-k^2)w/Re-dw/dt) It's not implemented.
c
      if (neumann.eq.1) then
c
c     Neumann bc: d(p+E)=H_2+(d2-k^2)v/Re-dv/dt. v and dv/dt at wall are
c     not zero when there is blowing/suction. Thus they are computed with
c     other terms. dv/dt at freestream is also computed.
c
c     Get an xy box for v
c
         call getxy(u3r(1,1,2),u3i(1,1,2),zbp,2,ur,ui)
c
c     Compute d2v at at wall and freestream by the weighting functions, wd1.
c
         do y=1,2
            do xz=1,nxz
               d2vr(xz,y)=0.
               d2vi(xz,y)=0.
            end do
         end do
         do y=1,nyp
            do xz=1,nxz
               d2vr(xz,1)=d2vr(xz,1)+u3r(xz,y,2)*wd1(y,2)
               d2vr(xz,2)=d2vr(xz,2)+u3r(xz,y,2)*wd1(y,4)
               d2vi(xz,1)=d2vi(xz,1)+u3i(xz,y,2)*wd1(y,2)
               d2vi(xz,2)=d2vi(xz,2)+u3i(xz,y,2)*wd1(y,4)
            end do
         end do
c
c     Add (d2-k^2)v/Re-dv/dt at the wall and in the free stream.
c     rv1r/rv1i is v at the previous time step.
c     dv/dt is approximated by first order differences and set
c     to zero dv/dt=0 for first iteration it=1.
c
c
c     This is for boundary layers
c
         c4=1./dtn*real(1-1/it)
         do xz=1,nxz
            pl=(xz-1)/(nx/2)
            vr(xz,1)=vr(xz,1)+(d2vr(xz,1)-k2(xz)*u3r(xz,1,2))/re
     &           -(u3r(xz,1,2)-rv1r(xz-nx/2*pl,zb+pl,1))*c4
            vr(xz,2)=vr(xz,2)+(d2vr(xz,2)-k2(xz)*u3r(xz,nyp,2))/re
     &           -(u3r(xz,nyp,2)-rv1r(xz-nx/2*pl,zb+pl,2))*c4
            vi(xz,1)=vi(xz,1)+(d2vi(xz,1)-k2(xz)*u3i(xz,1,2))/re
     &           -(u3i(xz,1,2)-rv1i(xz-nx/2*pl,zb+pl,1))*c4
            vi(xz,2)=vi(xz,2)+(d2vi(xz,2)-k2(xz)*u3i(xz,nyp,2))/re
     &           -(u3i(xz,nyp,2)-rv1i(xz-nx/2*pl,zb+pl,2))*c4
         end do
c
c     This is for channels with no-slip b.c.
c
c         do xz=1,nxz
c            pl=(xz-1)/(nx/2)
c            vr(xz,1)=vr(xz,1)+(d2vr(xz,1)-k2(xz)*u3r(xz,1,2))/re
c            vr(xz,2)=vr(xz,2)+(d2vr(xz,2)-k2(xz)*u3r(xz,nyp,2))/re
c            vi(xz,1)=vi(xz,1)+(d2vi(xz,1)-k2(xz)*u3i(xz,1,2))/re
c            vi(xz,2)=vi(xz,2)+(d2vi(xz,2)-k2(xz)*u3i(xz,nyp,2))/re
c         end do
c
c     Even/odd decomposition for bc
c
         do xz=1,nxz
            bc(xz,1,1)=0.5*(vr(xz,1)-vr(xz,2))
            bc(xz,2,1)=0.5*(vr(xz,1)+vr(xz,2))
            bc(xz,1,2)=0.5*(vi(xz,1)-vi(xz,2))
            bc(xz,2,2)=0.5*(vi(xz,1)+vi(xz,2))
         end do
      end if
c
c     Chebyshev tau method(p.130, Canuto)
c     Get the right-handed sides for the Chebyshev tau method
c
c     f contains fvr,fvi (i.e. the RHS)
c
      call crhsc(f,bc,nxz,ny,2,nxz,nyp)
c
c     Set up the system and solve for functions (bis contains pr,pi.)
c
      call setmac(q,c,d,e,k2,ny,nxz,nxz,neumann,cc)
      call tridc(bis(1,1,2),q,c,d,e,f,nxz,ny,2,nxz,nyp)
c
c     Correction for kx=kz=0 (i.e. Poisson equation)
c
      if (zb.eq.1) then
c
c     Integrate D<p+E>=<H_2>. Set the last equation to zero for CTM.
c     Boundary condition is E at the wall.
c
         if (.not.cim) temp(ny)=0.
         call icheb(temi,temp,bc0,nyp,1,1)
c
c     Fix the mean pressure <p> by setting zero at wall.
c
         do y=2,nyp-1,2
            temi(1)=temi(1)+temi(y)-temi(y+1)
         end do
c
c     Alternatively, set <p> zero at the other wall
c
c         do y=2,nyp-1,2
c            temi(1)=temi(1)-temi(y)-temi(y+1)
c         end do
c
c     Or as the mean between both walls
c
c         do y=2,nyp-1,2
c            temi(1)=temi(1)-temi(y+1)
c         end do
c
c     Save mean p+E, i.e. overwrite the solution from the Helmholtz eq.
c
         do y=1,nyp
            phr(1,y)=temi(y)
            phi(1,y)=0.
         end do
      end if
c
c     Now, we have p+E for all wavenumbers.
c
c     Compute the pressure p by subtracting E, i.e. p=(p+E)-E
c
      do y=1,nyp
         do xz=1,nxz
            phr(xz,y)=phr(xz,y)-hvr(xz,y)
            phi(xz,y)=phi(xz,y)-hvi(xz,y)
         end do
      end do
c
c     Pad the dealiasing region in y (if dealiasing in y used)
c
      if (ny+1.le.nyp) then
         do y=ny-3,max(nyp,ny+1)
c         do y=ny+1,max(nyp,ny+1)
            do xz=1,nxz
               d2phr(xz,y) = 0.0
               d2phi(xz,y) = 0.0
            end do
         end do
      end if
c
c     Go back to Fourier/physical space and put back into main storage
c
      call vchbb(d2phr,w3,nyp,nxz,nxz,1,prey)
      call vchbb(d2phi,w3,nyp,nxz,nxz,1,prey)
      call putxy(d2phr,d2phi,zbp,8,ur,ui)
c
c     Reconstruct the vorticities (have been overwritten in nonlinp)
c
      c4=2./real(nyp-1)
      do i=1,3
         call getxy(u3r(1,1,i),u3i(1,1,i),zbp,i,ur,ui)
         call vchbf(u3r(1,1,i),w3,nyp,nxz,nxz,1,prey)
         call vchbf(u3i(1,1,i),w3,nyp,nxz,nxz,1,prey)
         do y=1,ny
            do xz=1,nxz
               u3r(xz,y,i)=u3r(xz,y,i)*c4
               u3i(xz,y,i)=u3i(xz,y,i)*c4
            end do
         end do
      end do

      call dcheb(om3r       ,u3r(1,1,3),ny,nxz,nxz)
      call dcheb(om3r(1,1,3),u3r       ,ny,nxz,nxz)
      call dcheb(om3i       ,u3i(1,1,3),ny,nxz,nxz)
      call dcheb(om3i(1,1,3),u3i       ,ny,nxz,nxz)
      do y=1,ny
         do xz=1,nxz
            om3r(xz,y,1)= om3r(xz,y,1) +bbeta(xz)*u3i(xz,y,2)
            om3r(xz,y,3)=-om3r(xz,y,3) -alfa(xz) *u3i(xz,y,2)
            om3i(xz,y,1)= om3i(xz,y,1) -bbeta(xz)*u3r(xz,y,2)
            om3i(xz,y,3)=-om3i(xz,y,3) +alfa(xz) *u3r(xz,y,2)
         end do
      end do
c
c     Pad the dealiasing region in y (if dealiasing in y used)
c
      if (ny+1.le.nyp) then
         do y=ny+1,max(nyp,ny+1)
            do xz=1,nxz
               om3r(xz,y,1) = 0.0
               om3i(xz,y,1) = 0.0
               om3r(xz,y,3) = 0.0
               om3i(xz,y,3) = 0.0
            end do
         end do
      end if
c
c     Transform back to Fourier/Fourier/physical space
c
      call vchbb(om3r(1,1,1),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3i(1,1,1),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3r(1,1,3),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3i(1,1,3),w3,nyp,nxz,nxz,1,prey)
c
c     Store the boxes in the vorticities
c
      call putxy(om3r(1,1,1),om3i(1,1,1),zbp,4,ur,ui)
      call putxy(om3r(1,1,3),om3i(1,1,3),zbp,5,ur,ui)

      end subroutine linearp
