c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine linearbl(ur,ui,puw,zb,an,bn,anp1,bnp1,re,pr,alfa,beta,
     &     u0low,u0upp,w0low,w0upp,du0upp,px,ibc,cim,icorr,
     &     bu1jr,bu1ji,prey,h3r,h3i,hthr,hthi,pthr,pthi,dthr,dthi,
     &     u3r,u3i,om3r,om3i,th3r,th3i,
     &     f,fvr,fvi,fvh,fomyr,fomyi,pomyr,pomyi,fthr,fthi,
     &     bis,d2phr,d2phi,d2phh,phr,phi,phh,phh2,d2omyr,d2omyi,
     &     d2thr,d2thi,
     &     domyr,domyi,domyh,omyh,d2omyh,
     &     vr,vi,vh,vh2,dvr,dvi,dvh,dvh2,d2v,d2vr,d2vi,d2vh,d2vh2,
     &     pvr,pvi,homyr,homyi,hvr,hvi,q,c,d,e,w3,
     &     wallvr,wallvi,wbci,
     &     wallur,wallui,wallwr,wallwi,v_wall,wlvr,wlvi,
     &     iles,gur,gui,taur,taui,gu3r,gu3i,ggu3r,ggu3i,
     &     gthr,gthi,gsr,gsi,gth3r,gth3i,ggth3r,ggth3i,
     &     gewy1,gewy2,gewy12,gewy22,iord,
     &     diags,diags2,gplane,filtxz,filtxz2,
     &     lpfxz,lpfxz2,zbp,ihighorder,cs,my_node,tbc,
     &     dtheta0_upp,dtheta0_low,d2th_bc,fth_bc,dth3r_bc,dth3i_bc,
     &     dth3r_bc2,dth3i_bc2,theta0_low,theta0_upp,cflux,mflux,vsuc,
     &     wallthr,wallthi)
c
c     Advances all velocites a timestep (substep in case of rk3)
c     with the nonlinear terms calculated in step 2 as a driving force
c
c     Calculates the next timestep partial right hand-sides
c     for an xy-box i.e. nbz xy-planes
c
c     Changes :
c     v1.1 of bla
c     New boundary condition for freestream :
c     dudy=-ku for each velocity component and wavenumber k
c     this also implies dudy=0 (free slip) for the mean component.
c     this enables the freestream boundary condition to be moved
c     to lower y/deltastar
c     v2.0 of bla
c     boundary conditions in the case of spatial simulations
c     dudy=0 for each velocity component and wavenumber k
c     v2.2 of bla
c     new boundary conditions controllable from bla.i via the
c     parameter ibc
c
c     Calculates the next timestep partial right hand-sides
c     for an xy-box i.e. nbz xy-planes
c
c     f:   1         f,fvr
c          2         fvi
c          3         fvh
c          4         fomyr,pomyr
c          5         fomyi,pomyi
c          6         fthr 1
c          7         fthi 1
c          8         fthr 2
c          9         fthi 2
c          ...
c
c     bis: 1         bis,phh2
c          2         d2phr,phr
c          3         d2phi,phi
c          4         d2phh,phh
c          5         d2omyr
c          6         d2omyi
c          7         d2thr 1
c          8         d2thi 1
c          9         d2thr 2
c         10         d2thi 2
c          ...
c
c     d2v:1       d2v,d2vh2
c         2       d2vr,pvr
c         3       d2vi,pvi
c         4       d2vh
c
      implicit none

      include 'par.f'

      integer nxz
      parameter (nxz=nx/2*mbz)
c
c     ipdir: direction of pressure forcing
c            ipdir=1 x, ipdir=2 z
c
      integer ipdir
      parameter (ipdir = 1)

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real puw(ny,2+scalar)
      logical icorr,cim
      integer zb,ibc,my_node
      real an,bn,anp1,bnp1,re,pr(scalar)
      real u0low,u0upp,w0low,w0upp,du0upp,px
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
      real alfa(nxz),beta(nz)
      real prey(nyp*2+15)
      real h3r(nxz,nyp,3),h3i(nxz,nyp,3)
      real hthr(nxz,nyp,scalar),hthi(nxz,nyp,scalar)
      real u3r(nxz,nyp,3),u3i(nxz,nyp,3)
      real th3r(nxz,nyp,scalar),th3i(nxz,nyp,scalar)
      real om3r(nxz,nyp,3),om3i(nxz,nyp,3)
      real d2omyr(nxz,nyp),d2omyi(nxz,nyp)
      real domyr(nxz,nyp),domyi(nxz,nyp)
      real domyh(nxz,nyp),omyh(nxz,nyp)
      real d2omyh(nxz,nyp)
      real pomyr(nxz,nyp),pomyi(nxz,nyp)
      real fomyr(nxz,nyp),fomyi(nxz,nyp)
      real homyr(nxz,nyp),homyi(nxz,nyp)
      real d2phr(nxz,nyp),d2phi(nxz,nyp)
      real d2phh(nxz,nyp)
      real fthr (nxz,nyp,2*scalar),fthi (nxz,nyp,2*scalar) ! part of f
      real d2thr(nxz,nyp,2*scalar),d2thi(nxz,nyp,2*scalar) ! part of bis
      real fact
c
c     Variables for tbc>0
c
      real fth_bc (nxz,nyp,2*scalar)
      real d2th_bc(nxz,nyp,4*scalar)
      real dth3r_bc (nxz,nyp,scalar),dth3i_bc (nxz,nyp,scalar)
      real dth3r_bc2(nxz,nyp,scalar),dth3i_bc2(nxz,nyp,scalar)
      real dtheta0_upp(scalar),dtheta0_low(scalar)
      real theta0_low(scalar),theta0_upp(scalar)

      real phr(nxz,nyp),phi(nxz,nyp)
      real phh(nxz,nyp),phh2(nxz,nyp)
      real pthr(nxz,nyp,scalar),pthi(nxz,nyp,scalar)
      real d2vr(nxz,nyp),d2vi(nxz,nyp)
      real d2vh(nxz,nyp),d2vh2(nxz,nyp)
      real dvr(nxz,nyp),dvi(nxz,nyp)
      real dthr(nxz,nyp,scalar),dthi(nxz,nyp,scalar)
      real dvh(nxz,nyp),dvh2(nxz,nyp)
      real vr(nxz,nyp),vi(nxz,nyp)
      real vh(nxz,nyp),vh2(nxz,nyp)
      real fvr(nxz,nyp),fvi(nxz,nyp),fvh(nxz,nyp)
      real pvr(nxz,nyp),pvi(nxz,nyp)
      real hvr(nxz,nyp),hvi(nxz,nyp)
      real w3(nx/2,mbz,nyp)
      real q(nxz,ny+2),c(nxz,ny),d(nxz,ny),e(nxz,ny)
      real f(nxz,nyp,5+2*scalar),bis(nxz,nyp,6+2*scalar)
      real d2v(nxz,nyp,4)
      real wallvr(nxp/2+1,nzd),wallvi(nxp/2+1,nzd)
      real wallthr(nxp/2+1,nzd,scalar),wallthi(nxp/2+1,nzd,scalar)
      integer wbci
      real vsuc
c
c     Wall roughness and wall oscillations
c
      real wallur(nxp/2+1,memnz),wallui(nxp/2+1,memnz)
      real wlvr  (nxp/2+1,memnz),wlvi  (nxp/2+1,memnz)
      real wallwr(nxp/2+1,memnz),wallwi(nxp/2+1,memnz)
      logical v_wall
c
c     Local variables
c
      integer x,z,i,xz,y,is,pl,tbc(scalar),endi,ith
      real bbeta(nxz),k2(nxz),k2i(nxz),lam2(nxz),lam3(nxz,scalar)
      real k(nxz),k3(nxz,scalar)
      real fuw(ny,3+3*scalar),huw(ny,2+scalar),cuw(2,3+3*scalar)
      real uw(ny,3+3*scalar),duw(ny,3+3*scalar),d2uw(ny,3+scalar)
      real duw1(3+3*scalar),duw_1(3+3*scalar)
      real a11,a12,a21,a22,b1,b2,AA,BB
      real ctheta(nxz,6*scalar),dtheta1(6*scalar,nxz)
      real dtheta_1(6*scalar,nxz)

      real vbcr(nxz,2),vbci(nxz,2)
      real thbcr(nxz,2*scalar),thbci(nxz,2*scalar)
      real dvbcr(nxz,2),dvbci(nxz,2)
      real omybcr(nxz,2),omybci(nxz,2)
      real bc(nxz,2,5+4*scalar)
      real cnoll(nxz,2,5+2*scalar)
      real dvy1r(nxz,2),dvy1i(nxz,2),dvy1h(nxz,2)
      real dvy1h2(nxz,2)
      real cvr(nxz,2),cvi(nxz,2)
      real cvh2(nxz,2),cvh(nxz,2)
      real c1,c2,c3(scalar),c4(scalar)
      real q0(ny+2),c0(ny),d0(ny),e0(ny),bc0(2,3+3*scalar)
c
c     Evaluated derivatives at the boundary y=1
c
      real dj1ur(nxz,3,5),dj1ui(nxz,3,5)
      real dj1uh1(nxz,3,5),dj1uh2(nxz,3,5)
      real dj1uh3(nxz,3,4)
      real dj1omr(nxz,5),dj1omi(nxz,5)
      real dj1omh(nxz,5)
      real bcm(nxz,3,3),bcdr(nxz,3),bcdi(nxz,3)
c
c     LES
c
      integer iles,ll,iord,ihighorder
      real gur (memnx,memny,memnz,5)
      real gui (memnx,memny,memnz,5)
      real gsr (memnx,memny,memnz,scalar)
      real gsi (memnx,memny,memnz,scalar)
      real gthr (memnx,memny,memnz,2*scalar)
      real gthi (memnx,memny,memnz,2*scalar)
      real taur(memnx,memny,memnz,6+2*scalar)
      real taui(memnx,memny,memnz,6+2*scalar)
      real gu3r  (nx/2*mbz,nyp),gu3i (nx/2*mbz,nyp)
      real ggu3r (nx/2*mbz,nyp),ggu3i(nx/2*mbz,nyp)
      real gth3r  (nx/2*mbz,nyp),gth3i  (nx/2*mbz,nyp)
      real ggth3r (nx/2*mbz,nyp),ggth3i (nx/2*mbz,nyp)
      real gplane(nx/2,mbz,nyp)
      real gewy1 (nyp,5),gewy2 (nyp,5)
      real gewy12(nyp,5),gewy22(nyp,5)
      real diags (nyp,5)
      real diags2(nyp,5)
      real filtxz (nx/2,nz),lpfxz (nx/2,nz)
      real filtxz2(nx/2,nz),lpfxz2(nx/2,nz)
      real cs
c
c     Mass flux
c
      real mflux1,mflux2,mflux
      logical cflux
c
c     MPI
c
      integer zbp
c
c     Don't do it for oddball mode
c
      if (mbz.eq.1.and.zb.eq.(nz+1)/2+1) return
c
c     Get wavenumbers
c
      do z=1,mbz
         do x=1,nx/2
            xz=x+nx/2*(z-1)
            bbeta(xz)=beta(z+zb-1)
            k2(xz)=alfa(xz)*alfa(xz)+bbeta(xz)*bbeta(xz)
            lam2(xz)=k2(xz)+2.*re/(an+bn)
            k(xz)=sqrt(k2(xz))
         end do
      end do
      do ith=1,scalar
         do z=1,mbz
            do x=1,nx/2
               xz=x+nx/2*(z-1)
               lam3(xz,ith)=k2(xz)+2.*(pr(ith)*re)/(an+bn)
            end do
         end do
      end do

      if (zb.eq.1) then
         k2i(1)=0.
      else
         k2i(1)=1./k2(1)
      end if
      do xz=2,nxz
         k2i(xz)=1./k2(xz)
      end do
c
c     Get an xy box of the nonlinear term and
c     Chebyshev transform to Chebyshev space
c     and normalise
c
      c1=2./(real(nxp)*real(nyp-1)*real(nzp))
      do i=1,3
         call getxy(h3r(1,1,i),h3i(1,1,i),zbp,i,ur,ui)
         call vchbf(h3r(1,1,i),w3,nyp,nxz,nxz,1,prey)
         call vchbf(h3i(1,1,i),w3,nyp,nxz,nxz,1,prey)
         do y=1,ny
            do xz=1,nxz
               h3r(xz,y,i)=h3r(xz,y,i)*c1
               h3i(xz,y,i)=h3i(xz,y,i)*c1
            end do
         end do
      end do
c
c     Get xy boxes of 6 and 7 (already in Chebyshev space)
c
      call getxy(pvr,pvi,zbp,6,ur,ui)
      call getxy(pomyr,pomyi,zbp,7,ur,ui)

      do ith=1,scalar
c
c     Get xy boxes of 8 (non-linear term)
c     transform and normalise
c
         call getxy(hthr(1,1,ith),hthi(1,1,ith),zbp,
     &        8+pressure+3*(ith-1),ur,ui)
         call vchbf(hthr(1,1,ith),w3,nyp,nxz,nxz,1,prey)
         call vchbf(hthi(1,1,ith),w3,nyp,nxz,nxz,1,prey)
         c1=2./(real(nxp)*real(nyp-1)*real(nzp))
         do y=1,ny
            do xz=1,nxz
               hthr(xz,y,ith)=hthr(xz,y,ith)*c1
               hthi(xz,y,ith)=hthi(xz,y,ith)*c1
            end do
         end do
c
c     Get xy boxes of 10 (partial RHS of theta)
c
         call getxy(pthr(1,1,ith),pthi(1,1,ith),zbp,
     &        10+pressure+3*(ith-1),ur,ui)
      end do
c
c     Wavenumber zero, construct complete fu and fw
c
      c1=2.*re/(an+bn)*an
      if (zb.eq.1) then
         do y=1,ny
            fuw(y,1)=puw(y,1)-c1*h3r(1,y,1)
            fuw(y,2)=puw(y,2)-c1*h3r(1,y,3)
            huw(y,1)=h3r(1,y,1)
            huw(y,2)=h3r(1,y,3)
         end do
c
c     Add pressure gradient (if px is time dependent px=px(t+dt/2) for AB2)
c     actually px=0 for a flat plate boundary layer
c     The force is applied in direction x (ipdir=1) or z (ipdir=2).
c
         fuw(1,ipdir)=fuw(1,ipdir)+2.*re*px         
      end if
c
c     Calculate homy,fomy
c     partial hv=i*alfa*h3(,,1)+i*beta*h3(,,3)
c     homy=i*beta*h3(,,1)-i*alfa*h3(,,3)
c
      do y=1,ny
         do xz=1,nxz
            hvr(xz,y)=-alfa(xz)*h3i(xz,y,1)-bbeta(xz)*h3i(xz,y,3)
            hvi(xz,y)= alfa(xz)*h3r(xz,y,1)+bbeta(xz)*h3r(xz,y,3)
            homyr(xz,y)=-bbeta(xz)*h3i(xz,y,1)+alfa(xz)*h3i(xz,y,3)
            homyi(xz,y)= bbeta(xz)*h3r(xz,y,1)-alfa(xz)*h3r(xz,y,3)
            fomyr(xz,y)=pomyr(xz,y)-c1*homyr(xz,y)
            fomyi(xz,y)=pomyi(xz,y)-c1*homyi(xz,y)
         end do
      end do

      call dcheb(h3r(1,1,1),hvr,ny,nxz,nxz)
      call dcheb(h3i(1,1,1),hvi,ny,nxz,nxz)
c
c     Calculate hv,fv
c     Set up homogeneous right hand side for solution of phihom
c
      do y=1,ny
         do xz=1,nxz
            hvr(xz,y)=-(h3r(xz,y,1)+k2(xz)*h3r(xz,y,2))
            hvi(xz,y)=-(h3i(xz,y,1)+k2(xz)*h3i(xz,y,2))
            fvr(xz,y)=pvr(xz,y)-c1*hvr(xz,y)
            fvi(xz,y)=pvi(xz,y)-c1*hvi(xz,y)
            fvh(xz,y)=0.
         end do
      end do
      do ith=1,scalar
         fact = (2.*pr(ith)*re*an)/(an+bn)
         do y=1,ny
            do xz=1,nxz
               fthr(xz,y,1+(ith-1)*2)=pthr(xz,y,ith)-
     &              fact*hthr(xz,y,ith)
               fthi(xz,y,1+(ith-1)*2)=pthi(xz,y,ith)-
     &              fact*hthi(xz,y,ith)
            end do
         end do
      end do
c
c     For wavenumber zero fv and fomy should be zero
c
      if (zb.eq.1) then
         do ith=1,scalar
            do y=1,ny
               fuw(y,4+3*(ith-1))=puw(y,2+ith)-
     &              (2.*pr(ith)*re*an)/(an+bn)*hthr(1,y,ith)
               huw(y,2+ith)=hthr(1,y,ith)
               fthr(1,y,1+(ith-1)*2)=0.
               fthi(1,y,1+(ith-1)*2)=0.
            end do
         end do
         do y=1,ny
            fvr(1,y)=0.
            fvi(1,y)=0.
            fomyr(1,y)=0.
            fomyi(1,y)=0.
         end do
c
c     For symmetric cases fomy should be zero for beta zero
c     as well as fuw for w
c
         if (nfzsym.eq.1) then
            do y=1,ny
               fuw(y,2)=0.0
               do x=1,nx/2
                  fomyr(x,y)=0.0
                  fomyi(x,y)=0.0
               end do
            end do
            do ith=1,scalar
               do y=1,ny
                  do x=1,nx/2
                     fthr(x,y,1+(ith-1)*2)=0.
                     fthi(x,y,1+(ith-1)*2)=0.
                  end do
               end do
            end do
         end if
      end if
c
c     Get physical boundary conditions for y-even (i=1), y-odd (i=2)
c     note that the even boundary conditions for dv corresponds to the even
c     part of dv, not the even part of v (even v has odd dv)
c     IMPORTANT COMMENT : here only input the boundary conditions for
c     the lower surface (y=-1), the bc's at y=1 are imposed later
c
      do i=1,2
         do xz=1,nxz
            vbcr(xz,i)=0.
            vbci(xz,i)=0.
            dvbcr(xz,i)=0.
            dvbci(xz,i)=0.
            omybcr(xz,i)=0.
            omybci(xz,i)=0.
         end do
      end do
      do ith=1,scalar
         do i=1,2
            do xz=1,nxz
               thbcr(xz,i+2*(ith-1))=0.
               thbci(xz,i+2*(ith-1))=0.
            end do
         end do
      end do
c
c     Wall roughness
c     (Note that wallur/wr are smaller, thus one uses zbp 
c     instead of zb)
c
      if (wbci.eq.-1) then
         do xz=1,nxz
            pl=(xz-1)/(nx/2)
            dvbcr(xz,1)=0.5*wallur(xz-pl*nx/2,zbp+pl) !zbp
            dvbcr(xz,2)=-dvbcr(xz,1)
            dvbci(xz,1)=0.5*wallui(xz-pl*nx/2,zbp+pl)
            dvbci(xz,2)=-dvbci(xz,1)

            omybcr(xz,1)=0.5*wallwr(xz-pl*nx/2,zbp+pl)
            omybcr(xz,2)=-omybcr(xz,1)
            omybci(xz,1)=0.5*wallwi(xz-pl*nx/2,zbp+pl)
            omybci(xz,2)=-omybci(xz,1)
         end do

         if (v_wall) then
            do xz=1,nxz
               pl=(xz-1)/(nx/2)
               vbcr(xz,1)=0.5*wlvr(xz-pl*nx/2,zbp+pl)
               vbcr(xz,2)=-vbcr(xz,1)
               vbci(xz,1)=0.5*wlvi(xz-pl*nx/2,zbp+pl)
               vbci(xz,2)=-vbci(xz,1)
            end do
         end if
      end if

      if (wbci.eq.1.or.wbci.eq.2.or.wbci.eq.-2) then
c
c     The spatial distribution of the inhomogeneous Dirichlet b.c.
c     on the lower wall for the vertical velocity
c
         do xz=1,nxz
            pl=(xz-1)/(nx/2)
            vbcr(xz,1)=0.5*wallvr(xz-pl*nx/2,zb+pl)
            vbcr(xz,2)=-vbcr(xz,1)
            vbci(xz,1)=0.5*wallvi(xz-pl*nx/2,zb+pl)
            vbci(xz,2)=-vbci(xz,1)
         end do
c
c     The spatial distribution of the inhomogeneous Dirichlet b.c.
c     on the lower wall for the scalar
c
         do ith=1,scalar
            do xz=1,nxz
               pl=(xz-1)/(nx/2)
               thbcr(xz,1+2*(ith-1))= 0.5*wallthr(xz-pl*nx/2,zb+pl,ith)
               thbcr(xz,2+2*(ith-1))=-thbcr(xz,1+2*(ith-1))
               thbci(xz,1+2*(ith-1))= 0.5*wallthi(xz-pl*nx/2,zb+pl,ith)
               thbci(xz,2+2*(ith-1))=-thbci(xz,1+2*(ith-1))
            end do
         end do

      end if

      if (wbci.eq.3.or.wbci.eq.4.or.wbci.eq.5) then
c
c     The spatial distribution of the inhomogeneous Dirichlet b.c.
c     on the lower wall for the in-plane velocity
c     We use zb instead of zbp because of the array size of wallur, ...
c
         do xz=1,nxz
            pl=(xz-1)/(nx/2)
            dvbcr(xz,1)=0.5*wallur(xz-pl*nx/2,zbp+pl) 
            dvbcr(xz,2)=-dvbcr(xz,1)
            dvbci(xz,1)=0.5*wallui(xz-pl*nx/2,zbp+pl)
            dvbci(xz,2)=-dvbci(xz,1)
            
            omybcr(xz,1)=0.5*wallwr(xz-pl*nx/2,zbp+pl)
            omybcr(xz,2)=-omybcr(xz,1)
            omybci(xz,1)=0.5*wallwi(xz-pl*nx/2,zbp+pl)
            omybci(xz,2)=-omybci(xz,1)
         end do
      end if
c
c     Calculate omy,phipart,phihom
c     boundary conditions at y=+1/-1 for omy: 0 & 0,  phipart: 0 & 0,
c     and phihom: sym:1&1; antisym:1&-1 => (alfa +/- beta)/2 = 1 for
c     both symmetric and antisymmetric cases.
c
      do i=1,2
         do xz=1,nxz
            bc(xz,i,1)=0.
            bc(xz,i,2)=0.
            bc(xz,i,3)=1.
            bc(xz,i,4)=omybcr(xz,i)
            bc(xz,i,5)=omybci(xz,i)
         end do
      end do
      do ith=1,scalar
         do i=1,2
            do xz=1,nxz
               bc(xz,i,6+4*(ith-1)) = thbcr(xz,i+2*(ith-1))
               bc(xz,i,7+4*(ith-1)) = thbci(xz,i+2*(ith-1))
            end do
         end do
      end do

      if (cim) then
c
c     Chebyshev integration method
c     Solve for second derivatives
c     note that we skip the first component bis(,,1) as this corresponds
c     to the vhom2 component (wave solution)
c     (bis contains d2phr,d2phi,d2phh,d2omyr,d2omyi,
c     f contains fvr,fvi,fvh,fomyr,fomyi)
c
         call setmat(q,c,d,e,lam2,ny,nxz,nxz,cim)
         call trid(bis(1,1,2),cnoll,q,c,d,e,f,bc,nxz,ny+2,5,nxz,nyp,cim)

         do ith=1,scalar
            call setmat(q,c,d,e,lam3(1,ith),ny,nxz,nxz,cim)
            call trid(bis(1,1,7+2*(ith-1)),cnoll(1,1,6+2*(ith-1)),q,c,
     &           d,e,f(1,1,6+2*(ith-1)),
     &           bc(1,1,6+4*(ith-1)),nxz,ny+2,2,nxz,nyp,cim)
         end do
c
c     Integrate phipart twice
c
         call icheb(w3,d2phr,cnoll(1,2,1),ny,nxz,nxz)
         call icheb(phr,w3,cnoll(1,1,1),ny,nxz,nxz)
         call icheb(w3,d2phi,cnoll(1,2,2),ny,nxz,nxz)
         call icheb(phi,w3,cnoll(1,1,2),ny,nxz,nxz)
c
c     Integrate phihom twice, note that this is also omyhom
c
         if (ibc.ne.120.and.ibc.ne.20) then
c
c     Get d2omyhom
c
            do y=1,ny
               do xz=1,nxz
                  d2omyh(xz,y)=d2phh(xz,y)
               end do
            end do
            if (icorr) call corrd2(d2omyh,ny,nxz,nxz)
         end if
         call icheb(domyh,d2phh,cnoll(1,2,3),ny,nxz,nxz)
         call icheb(phh,domyh,cnoll(1,1,3),ny,nxz,nxz)
         if (ibc.ne.120.and.ibc.ne.20) then
            if (icorr) call corrd(domyh,ny,nxz,nxz)
c
c     Get omyhom
c
            do y=1,ny
               do xz=1,nxz
                  omyh(xz,y)=phh(xz,y)
               end do
            end do
         end if
c
c     Integrate d2omy twice into om3(,,2)
c
         call icheb(domyr,d2omyr,cnoll(1,2,4),ny,nxz,nxz)
         call icheb(domyi,d2omyi,cnoll(1,2,5),ny,nxz,nxz)
         call icheb(om3r(1,1,2),domyr,cnoll(1,1,4),ny,nxz,nxz)
         call icheb(om3i(1,1,2),domyi,cnoll(1,1,5),ny,nxz,nxz)
         if (icorr) then
            call corrd2(d2omyr,ny,nxz,nxz)
            call corrd2(d2omyi,ny,nxz,nxz)
            call corrd(domyr,ny,nxz,nxz)
            call corrd(domyi,ny,nxz,nxz)
         end if

         do ith=1,scalar
c
c     Integrate d2th twice into th3(,,1)
c
            call icheb(dthr(1,1,ith),d2thr(1,1,1+(ith-1)*2),
     &           cnoll(1,2,6+2*(ith-1)),ny,nxz,nxz)
            call icheb(dthi(1,1,ith),d2thi(1,1,1+(ith-1)*2),
     &           cnoll(1,2,7+2*(ith-1)),ny,nxz,nxz)
            call icheb(th3r(1,1,ith),dthr(1,1,ith) ,
     &           cnoll(1,1,6+2*(ith-1)),ny,nxz,nxz)
            call icheb(th3i(1,1,ith),dthi(1,1,ith) ,
     &           cnoll(1,1,7+2*(ith-1)),ny,nxz,nxz)
            if (icorr) then
               call corrd2(d2thr(1,1,1+(ith-1)*2),ny,nxz,nxz)
               call corrd2(d2thi(1,1,1+(ith-1)*2),ny,nxz,nxz)
               call corrd(dthr(1,1,ith),ny,nxz,nxz)
               call corrd(dthi(1,1,ith),ny,nxz,nxz)
            end if
c
c     dthr,dthi is finished apart from zero mode
c
         end do

      else
c
c     Chebyshev tau method (Canuto, p.130)
c     get the right-handed sides for the Chebyshev tau method
c     (f contains fvr,fvi,fvh,fomyr,fomyi.)
c     solve for functions(bis contains phr,phi,phh,om3r,om3i.
c
         call setmacrhs(bis(1,1,2),f,c,d,e,lam2,nxz,ny,5,nxz,nyp)
         call trid(bis(1,1,2),cnoll,q,c,d,e,f,bc,nxz,ny,5,nxz,nyp,cim)
c
c     Differentiate phh
c
         call dcheb(domyh,phh,ny,nxz,nxz)
         if (ibc.ne.120.and.ibc.ne.20) then
c
c     Get omyhom
c
            do y=1,ny
               do xz=1,nxz
                  omyh(xz,y)=phh(xz,y)
               end do
            end do
c
c     Get d2omyhom
c
            call dcheb(d2omyh,domyh,ny,nxz,nxz)
         end if
         do y=1,ny
            do xz=1,nxz
               om3r(xz,y,2) = d2omyr(xz,y)
               om3i(xz,y,2) = d2omyi(xz,y)
            end do
         end do
c
c     Differentiate omyr, omyi twice
c
         call dcheb(domyr,om3r(1,1,2),ny,nxz,nxz)
         call dcheb(d2omyr,domyr,ny,nxz,nxz)
         call dcheb(domyi,om3i(1,1,2),ny,nxz,nxz)
         call dcheb(d2omyi,domyi,ny,nxz,nxz)

         do ith=1,scalar
            call setmacrhs(bis(1,1,7+2*(ith-1)),f(1,1,6+2*(ith-1)),
     &           c,d,e,lam3(1,ith),nxz,ny,2,nxz,nyp)
            call trid(bis(1,1,7+2*(ith-1)),cnoll(1,1,6+2*(ith-1)),
     &           q,c,d,e,f(1,1,6+2*(ith-1)),bc(1,1,6+4*(ith-1)),
     &           nxz,ny,2,nxz,nyp,cim)
c
c     Differentiate once to get dth (th3 --> dth)
c
            do y=1,ny
               do xz=1,nxz
                  th3r(xz,y,ith)=d2thr(xz,y,1+(ith-1)*2)
                  th3i(xz,y,ith)=d2thi(xz,y,1+(ith-1)*2)
               end do
            end do
            call dcheb(dthr(1,1,ith),th3r(1,1,ith),ny,nxz,nxz)
            call dcheb(dthi(1,1,ith),th3i(1,1,ith),ny,nxz,nxz)
c
c     dthr,dthi is finished apart from zero mode
c
c
c     If special derivative boundary conditions for scalar are needed
c
            if (tbc(ith).gt.0) then
               do xz=1,nxz
                  ctheta(xz,1+6*(ith-1)) = th3r(xz,1,ith)
                  ctheta(xz,2+6*(ith-1)) = th3i(xz,1,ith)
c
c     Solve for homogeneous solution with bc at -1 and/or +1
c
                  if (tbc(ith).eq.1) then
c
c     Derivate at lower boundary
c
                     bc(xz,1,6+4*(ith-1)) =  1./2.
                     bc(xz,2,6+4*(ith-1)) = -1./2.
                     bc(xz,1,7+4*(ith-1)) =  1./2.
                     bc(xz,2,7+4*(ith-1)) = -1./2.
                  else if (tbc(ith).eq.2) then
c
c     Derivative at upper boundary
c
                     bc(xz,1,6+4*(ith-1)) =  1./2.
                     bc(xz,2,6+4*(ith-1)) =  1./2.
                     bc(xz,1,7+4*(ith-1)) =  1./2.
                     bc(xz,2,7+4*(ith-1)) =  1./2.
                  else if (tbc(ith).eq.3) then
c
c     Two derivatives at lower and upper boundary
c
                     bc(xz,1,6+4*(ith-1)) =  1./2.
                     bc(xz,2,6+4*(ith-1)) = -1./2.
                     bc(xz,1,7+4*(ith-1)) =  1./2.
                     bc(xz,2,7+4*(ith-1)) = -1./2.

                     bc(xz,1,8+4*(ith-1)) =  1./2.
                     bc(xz,2,8+4*(ith-1)) =  1./2.
                     bc(xz,1,9+4*(ith-1)) =  1./2.
                     bc(xz,2,9+4*(ith-1)) =  1./2.
                  end if
                  do y=1,ny
                     fth_bc(xz,y,1+2*(ith-1))=0.
                     fth_bc(xz,y,2+2*(ith-1))=0.
                  end do
               end do
c
c     Solve first system (position 6 and 7)
c
               call setmacrhs(d2th_bc(1,1,1+4*(ith-1)),
     &              fth_bc(1,1,1+2*(ith-1)),c,d,e,lam3(1,ith),
     &              nxz,ny,2,nxz,nyp)
               call trid(d2th_bc(1,1,1+4*(ith-1)),
     &              cnoll(1,1,6+2*(ith-1)),q,c,d,e,
     &              fth_bc(1,1,1+2*(ith-1)),bc(1,1,6+4*(ith-1)),nxz,
     &              ny,2,nxz,nyp,cim)

               call dcheb(dth3r_bc(1,1,ith),d2th_bc(1,1,1+4*(ith-1)),
     &              ny,nxz,nxz)
               call dcheb(dth3i_bc(1,1,ith),d2th_bc(1,1,2+4*(ith-1)),
     &              ny,nxz,nxz)

               do xz=1,nxz
                  ctheta(xz,3+6*(ith-1))=d2th_bc(xz,1,1+4*(ith-1))
                  ctheta(xz,4+6*(ith-1))=d2th_bc(xz,1,2+4*(ith-1))
               end do
               if (tbc(ith).eq.3) then
c
c     Solve second system (positions 8 and 9)
c
                  call setmacrhs(d2th_bc(1,1,3+4*(ith-1)),
     &                 fth_bc(1,1,1+2*(ith-1)),c,d,e,lam3(1,ith),
     &                 nxz,ny,2,nxz,nyp)
                  call trid(d2th_bc(1,1,3+4*(ith-1)),
     &                 cnoll(1,1,6+2*(ith-1)),q,c,d,e,
     &                 fth_bc(1,1,1+2*(ith-1)),bc(1,1,8+4*(ith-1)),
     &                 nxz,ny,2,nxz,nyp,cim)
                  call dcheb(dth3r_bc2(1,1,ith),d2th_bc(1,1,3+4*(ith-1))
     &                 ,ny,nxz,nxz)
                  call dcheb(dth3i_bc2(1,1,ith),d2th_bc(1,1,4+4*(ith-1))
     &                 ,ny,nxz,nxz)
                  do xz=1,nxz
                     ctheta(xz,5+6*(ith-1))=d2th_bc(xz,1,3+4*(ith-1))
                     ctheta(xz,6+6*(ith-1))=d2th_bc(xz,1,4+4*(ith-1))
                  end do
               end if
c
c     Now get the derivative at the wall at y=-1 or y=1
c
               if (tbc(ith).eq.3) then
                  do xz=1,nxz
                     dtheta1 (1+6*(ith-1),xz)=0.0
                     dtheta1 (2+6*(ith-1),xz)=0.0
                     dtheta1 (3+6*(ith-1),xz)=0.0
                     dtheta1 (4+6*(ith-1),xz)=0.0
                     dtheta1 (5+6*(ith-1),xz)=0.0
                     dtheta1 (6+6*(ith-1),xz)=0.0
                     dtheta_1(1+6*(ith-1),xz)=0.0
                     dtheta_1(2+6*(ith-1),xz)=0.0
                     dtheta_1(3+6*(ith-1),xz)=0.0
                     dtheta_1(4+6*(ith-1),xz)=0.0
                     dtheta_1(5+6*(ith-1),xz)=0.0
                     dtheta_1(6+6*(ith-1),xz)=0.0

                     do y=1,ny
c
c     Derivative 1 (upper wall)
c
                        dtheta1(1+6*(ith-1),xz)=dtheta1(1+6*(ith-1),xz)
     &                       +dthr(xz,y,ith)
                        dtheta1(2+6*(ith-1),xz)=dtheta1(2+6*(ith-1),xz)
     &                       +dthi(xz,y,ith)
                        dtheta1(3+6*(ith-1),xz)=dtheta1(3+6*(ith-1),xz)
     &                       +dth3r_bc(xz,y,ith)
                        dtheta1(4+6*(ith-1),xz)=dtheta1(4+6*(ith-1),xz)
     &                       +dth3i_bc(xz,y,ith)
                        dtheta1(5+6*(ith-1),xz)=dtheta1(5+6*(ith-1),xz)
     &                       +dth3r_bc2(xz,y,ith)
                        dtheta1(6+6*(ith-1),xz)=dtheta1(6+6*(ith-1),xz)
     &                       +dth3i_bc2(xz,y,ith)
c
c     Derivative -1 (lower wall)
c
                        dtheta_1(1+6*(ith-1),xz)=dtheta_1(1+6*(ith-1),xz
     &                       )+dthr(xz,y,ith)*(-1)**(y-1)
                        dtheta_1(2+6*(ith-1),xz)=dtheta_1(2+6*(ith-1),xz
     &                       )+dthi(xz,y,ith)*(-1)**(y-1)
                        dtheta_1(3+6*(ith-1),xz)=dtheta_1(3+6*(ith-1),xz
     &                       )+dth3r_bc(xz,y,ith)*(-1)**(y-1)
                        dtheta_1(4+6*(ith-1),xz)=dtheta_1(4+6*(ith-1),xz
     &                       )+dth3i_bc(xz,y,ith)*(-1)**(y-1)
                        dtheta_1(5+6*(ith-1),xz)=dtheta_1(5+6*(ith-1),xz
     &                       )+dth3r_bc2(xz,y,ith)*(-1)**(y-1)
                        dtheta_1(6+6*(ith-1),xz)=dtheta_1(6+6*(ith-1),xz
     &                       )+dth3i_bc2(xz,y,ith)*(-1)**(y-1)
                     end do
                  end do
               else if (tbc(ith).eq.2) then
                  do xz=1,nxz
                     dtheta1(1+6*(ith-1),xz)=0.0
                     dtheta1(2+6*(ith-1),xz)=0.0
                     dtheta1(3+6*(ith-1),xz)=0.0
                     dtheta1(4+6*(ith-1),xz)=0.0
                     do y=1,ny
c
c     Derivative 1 (upper wall)
c
                        dtheta1(1+6*(ith-1),xz)=dtheta1(1+6*(ith-1),xz)
     $                       +dthr(xz,y,ith)
                        dtheta1(2+6*(ith-1),xz)=dtheta1(2+6*(ith-1),xz)
     $                       +dthi(xz,y,ith)
                        dtheta1(3+6*(ith-1),xz)=dtheta1(3+6*(ith-1),xz)
     $                       +dth3r_bc(xz,y,ith)
                        dtheta1(4+6*(ith-1),xz)=dtheta1(4+6*(ith-1),xz)
     $                       +dth3i_bc(xz,y,ith)
                     end do
                  end do
               else if (tbc(ith).eq.1) then
                  do xz=1,nxz
                     dtheta1(1+6*(ith-1),xz)=0.0
                     dtheta1(2+6*(ith-1),xz)=0.0
                     dtheta1(3+6*(ith-1),xz)=0.0
                     dtheta1(4+6*(ith-1),xz)=0.0
                     do y=1,ny
c
c     Derivative -1 (lower wall)
c
                        dtheta1(1+6*(ith-1),xz)=dtheta1(1+6*(ith-1),xz)
     $                       +dthr(xz,y,ith)*(-1)**(y-1)
                        dtheta1(2+6*(ith-1),xz)=dtheta1(2+6*(ith-1),xz)
     $                       +dthi(xz,y,ith)*(-1)**(y-1)
                        dtheta1(3+6*(ith-1),xz)=dtheta1(3+6*(ith-1),xz)
     $                       +dth3r_bc(xz,y,ith)*(-1)**(y-1)
                        dtheta1(4+6*(ith-1),xz)=dtheta1(4+6*(ith-1),xz)
     $                       +dth3i_bc(xz,y,ith)*(-1)**(y-1)
                     end do
                  end do
               end if
c
c     Superimpose to get the prescribed value of the dthetady at y=+1,-1
c
               if (tbc(ith).lt.3) then
c
c     One derivative imposed (tbc=1 or tbc=2)
c
                  do xz=1,nxz
                     ctheta(xz,1+6*(ith-1))=ctheta(xz,1+6*(ith-1))
     $                    -(dtheta1(1+6*(ith-1),xz)/dtheta1(3+6*(ith-1)
     $                    ,xz))*ctheta(xz,3+6*(ith-1))
                     ctheta(xz,2+6*(ith-1))=ctheta(xz,2+6*(ith-1))
     $                    -(dtheta1(2+6*(ith-1),xz)/dtheta1(4+6*(ith-1)
     $                    ,xz))*ctheta(xz,4+6*(ith-1))
                     do y=1,ny
                        dthr(xz,y,ith)=dthr(xz,y,ith) -(dtheta1(1+6*(ith
     $                       -1),xz)/dtheta1(3+6*(ith-1),xz))
     $                       *dth3r_bc(xz,y,ith)
                        dthi(xz,y,ith)=dthi(xz,y,ith) -(dtheta1(2+6*(ith
     $                       -1),xz)/dtheta1(4+6*(ith-1),xz))
     $                       *dth3i_bc(xz,y,ith)

                     end do
                  end do
                  call icheb(th3r(1,1,ith),dthr(1,1,ith),ctheta(1,1+6
     $                 *(ith-1)),ny,nxz,nxz)
                  call icheb(th3i(1,1,ith),dthi(1,1,ith),ctheta(1,2+6
     $                 *(ith-1)),ny,nxz,nxz)
               else if (tbc(ith).eq.3) then
c
c     Two derivatives imposed (tbc=3): System solve necessary
c
                  do xz=1,nxz
                     a11=dtheta1(3+6*(ith-1),xz)
                     a12=dtheta1(5+6*(ith-1),xz)
                     a21=dtheta_1(3+6*(ith-1),xz)
                     a22=dtheta_1(5+6*(ith-1),xz)
                     b1=-dtheta1(1+6*(ith-1),xz)
                     b2=-dtheta_1(1+6*(ith-1),xz)
c
c     Solve and get AA BB
c
                     AA= (a22*b1-b2*a12)/(a22*a11-a12*a21)
                     BB=-(a21*b1-b2*a11)/(a22*a11-a12*a21)

                     ctheta(xz,1+6*(ith-1))=ctheta(xz,1+6*(ith-1))+ AA
     $                    *ctheta(xz,3+6*(ith-1))+BB*ctheta(xz,5+6*(ith
     $                    -1))
                     do y=1,ny
                        dthr(xz,y,ith)=dthr(xz,y,ith)+AA*dth3r_bc(xz,y
     $                       ,ith)+BB*dth3r_bc2(xz,y,ith)
                     end do

                     a11=dtheta1(4+6*(ith-1),xz)
                     a12=dtheta1(6+6*(ith-1),xz)
                     a21=dtheta_1(4+6*(ith-1),xz)
                     a22=dtheta_1(6+6*(ith-1),xz)
                     b1=-dtheta1(2+6*(ith-1),xz)
                     b2=-dtheta_1(2+6*(ith-1),xz)
c
c     Solve and get AA BB
c
                     AA= (a22*b1-b2*a12)/(a22*a11-a12*a21)
                     BB=-(a21*b1-b2*a11)/(a22*a11-a12*a21)
                     ctheta(xz,2+6*(ith-1))=ctheta(xz,2+6*(ith-1))+ AA
     $                    *ctheta(xz,4+6*(ith-1))+BB*ctheta(xz,6+6*(ith
     $                    -1))
                     do y=1,ny
                        dthi(xz,y,ith)=dthi(xz,y,ith)+AA*dth3i_bc(xz,y
     $                       ,ith)+BB*dth3i_bc2(xz,y,ith)
                     end do
                  end do
                  call icheb(th3r(1,1,ith),dthr(1,1,ith),ctheta(1,1+6
     $                 *(ith-1)),ny,nxz,nxz)
                  call icheb(th3i(1,1,ith),dthi(1,1,ith),ctheta(1,2+6
     $                 *(ith-1)),ny,nxz,nxz)
               end if
            end if
         end do
      end if
c
c     Calculate vpart,vhom
c     boundary conditions for vpart, vhom
c     vhom has homogeneous bc, whilst vpart gets the physical bc
c
      do i=1,2
         do xz=1,nxz
            bc(xz,i,1)=1.
            bc(xz,i,2)=vbcr(xz,i)
            bc(xz,i,3)=vbci(xz,i)
            bc(xz,i,4)=0.
         end do
      end do
c
c     Set a homogeneous RHS for vh2
c
      do y=1,ny
         do xz=1,nxz
            phh2(xz,y)=0.
         end do
      end do
      if (cim) then
c
c     Solve for second derivatives, d2v contains d2vh2,d2vr,d2vi,d2vh
c     on exit from trid. bis contains phh2,phr,phi,phh
c
         call setmat(q,c,d,e,k2,ny,nxz,nxz,cim)
         call trid(d2v,cnoll,q,c,d,e,bis,bc,nxz,ny+2,4,nxz,nyp,cim)
         if (icorr) then
            call corrp(phr,d2vr,ny,nxz,nxz)
            call corrp(phi,d2vi,ny,nxz,nxz)
            call corrp(phh,d2vh,ny,nxz,nxz)
            call corrp(phh2,d2vh2,ny,nxz,nxz)
         end if
c
c     Integrate vpart twice
c
         call icheb(dvr,d2vr,cnoll(1,2,2),ny,nxz,nxz)
         call icheb(vr,dvr,cnoll(1,1,2),ny,nxz,nxz)
         call icheb(dvi,d2vi,cnoll(1,2,3),ny,nxz,nxz)
         call icheb(vi,dvi,cnoll(1,1,3),ny,nxz,nxz)
c
c     Integrate vhom twice
c
         call icheb(dvh,d2vh,cnoll(1,2,4),ny,nxz,nxz)
         call icheb(vh,dvh,cnoll(1,1,4),ny,nxz,nxz)
         call icheb(dvh2,d2vh2,cnoll(1,2,1),ny,nxz,nxz)
         call icheb(vh2,dvh2,cnoll(1,1,1),ny,nxz,nxz)
         if (icorr) then
            call corrd2(d2vr,ny,nxz,nxz)
            call corrd(dvr,ny,nxz,nxz)
            call corrd2(d2vi,ny,nxz,nxz)
            call corrd(dvi,ny,nxz,nxz)
            call corrd2(d2vh,ny,nxz,nxz)
            call corrd(dvh,ny,nxz,nxz)
            call corrd2(d2vh2,ny,nxz,nxz)
            call corrd(dvh2,ny,nxz,nxz)
         end if
      else
c
c     Get the right-handed sides for the Chebyshev tau method
c     bis contains phh2,phr,phi,phh.
c     Solve for functions, d2v contains vh2,vr,vi,vh on exit from trid.
c
         call setmacrhs(d2v,bis,c,d,e,k2,nxz,ny,4,nxz,nyp)
         call trid(d2v,cnoll,q,c,d,e,bis,bc,nxz,ny,4,nxz,nyp,cim)
         do y=1,ny
            do xz=1,nxz
               vr(xz,y)=d2vr(xz,y)
               vi(xz,y)=d2vi(xz,y)
               vh(xz,y)=d2vh(xz,y)
               vh2(xz,y)=d2vh2(xz,y)
            end do
         end do
c
c     Differentiate vpart, vhom, vh2 twice.
c
         call dcheb(dvr,vr,ny,nxz,nxz)
         call dcheb(dvi,vi,ny,nxz,nxz)
         call dcheb(dvh,vh,ny,nxz,nxz)
         call dcheb(d2vh,dvh,ny,nxz,nxz)
         call dcheb(dvh2,vh2,ny,nxz,nxz)
      end if
c
c     Superimpose for v,dv/dy,phi
c     Evaluate dvdy at y=1
c
      do is=0,1
         do xz=1,nxz
            dvy1r(xz,is+1)=0.
            dvy1i(xz,is+1)=0.
            dvy1h(xz,is+1)=0.
            dvy1h2(xz,is+1)=0.
         end do
         do y=1+is,ny,2
            do xz=1,nxz
               dvy1r(xz,is+1)=dvy1r(xz,is+1)+dvr(xz,y)
               dvy1i(xz,is+1)=dvy1i(xz,is+1)+dvi(xz,y)
               dvy1h(xz,is+1)=dvy1h(xz,is+1)+dvh(xz,y)
               dvy1h2(xz,is+1)=dvy1h2(xz,is+1)+dvh2(xz,y)
            end do
         end do
c
c     Note that the even parts derivative is odd and vice versa
c     this accounts for the reversing done by 2-is
c
         do xz=1,nxz
            c1=1./dvy1h(xz,is+1)
            cvr(xz,2-is)=(dvbcr(xz,is+1)-dvy1r(xz,is+1))*c1
            cvi(xz,2-is)=(dvbci(xz,is+1)-dvy1i(xz,is+1))*c1
            cvh2(xz,2-is)=-dvy1h2(xz,is+1)*c1
            cvh(xz,2-is)=c1
         end do
      end do

      do is=0,1
         do y=1+is,ny,2
            do xz=1,nxz
               c1=dvh(xz,y)
               dvr(xz,y)=dvr(xz,y)+cvr(xz,2-is)*c1
               dvi(xz,y)=dvi(xz,y)+cvi(xz,2-is)*c1
               dvh2(xz,y)=dvh2(xz,y)+cvh2(xz,2-is)*c1
               c1=d2vh(xz,y)
               d2vr(xz,y)=d2vr(xz,y)+cvr(xz,is+1)*c1
               d2vi(xz,y)=d2vi(xz,y)+cvi(xz,is+1)*c1
               d2vh2(xz,y)=d2vh2(xz,y)+cvh2(xz,is+1)*c1
               c1=vh(xz,y)
               u3r(xz,y,2)=vr(xz,y)+cvr(xz,is+1)*c1
               u3i(xz,y,2)=vi(xz,y)+cvi(xz,is+1)*c1
               vh2(xz,y)=vh2(xz,y)+cvh2(xz,is+1)*c1
c
c     The vh homogeneous solution fulfills the following
c     boundary conditions :
c     vh(y=+-1)=0; laplace vh(y=+1)=2
c     now use vh to find vh1 which satisfies the following boundary
c     conditions : vh1(y=+-1)=0; dvh1(y=-1)=0; dvh(y=+1)=2
c     C is an arbitrary constant.
c     this is done by normalizing the derivative of each part to one
c     at y=1 (this cancels the derivative at y=-1)
c     we will store vh1 in vh's place
c
               vh(xz,y)=vh(xz,y)*cvh(xz,1+is)
               d2vh(xz,y)=d2vh(xz,y)*cvh(xz,1+is)
               dvh(xz,y)=dvh(xz,y)*cvh(xz,2-is)
            end do
         end do
      end do
      if (.not.cim) then
c
c     For Chebyshev-tau method
c
         call dcheb(vr,u3r(1,1,2),ny,nxz,nxz)
         call dcheb(d2vr,vr,ny,nxz,nxz)
         call dcheb(vr,u3i(1,1,2),ny,nxz,nxz)
         call dcheb(d2vi,vr,ny,nxz,nxz)
         call dcheb(vr,vh2,ny,nxz,nxz)
         call dcheb(d2vh2,vr,ny,nxz,nxz)
      end if
c
c     Now find the coefficients for vh1, vh2 and omyh so that
c     we fulfill the selected boundary conditions on u and w.
c     First evaluate the relevant derivatives at y=1 for the particular
c     and homogeneous v and omy solutions
c
      call evald(dj1ur,dj1ui,dj1uh1,dj1uh2,dj1uh3,dj1omr,dj1omi,dj1omh,
     &     d2vr,d2vi,domyr,domyi,d2omyr,d2omyi,d2vh,d2vh2,domyh,d2omyh,
     &     bbeta,k2i,alfa,ibc)
c
c     Then evaluate desired boundary condition selected by ibc
c     for the homogeneous solution and the particular solution
c     and calculate the discrepancy between the bc of the
c     particular solution and the desired ones
c
      call evalbc(bcm,bcdr,bcdi,ibc,zb,
     &     dj1omr,dj1omi,dj1omh,
     &     dj1ur,dj1ui,dj1uh1,dj1uh2,dj1uh3,k,bu1jr,bu1ji,alfa)
c
c     Solve for the coefficients of the homogeneous solutions
c
      call solvbc(bcm,bcdr,bcdi,zb)
c
c     Now finally superimpose the homogeneous solutions on the
c     particular ones v=vp-bcd(,1)*vh-bcd(,2)*vh2, omy=omyp-bcd(,3)*omyh
c
      if (ibc.ne.120.and.ibc.ne.20) then
         do y=1,ny
            do xz=1,nxz
               om3r(xz,y,2)=om3r(xz,y,2)-bcdr(xz,3)*omyh(xz,y)
               om3i(xz,y,2)=om3i(xz,y,2)-bcdi(xz,3)*omyh(xz,y)
               domyr(xz,y)=domyr(xz,y)-bcdr(xz,3)*domyh(xz,y)
               domyi(xz,y)=domyi(xz,y)-bcdi(xz,3)*domyh(xz,y)
            end do
         end do
         if (.not.cim)call dcheb(domyr,om3r(1,1,2),ny,nxz,nxz)
         if (.not.cim)call dcheb(domyi,om3i(1,1,2),ny,nxz,nxz)
      end if

      do y=1,ny
         do xz=1,nxz
            u3r(xz,y,2)=u3r(xz,y,2)-bcdr(xz,1)*vh(xz,y)-
     &           bcdr(xz,2)*vh2(xz,y)
            u3i(xz,y,2)=u3i(xz,y,2)-bcdi(xz,1)*vh(xz,y)-
     &           bcdi(xz,2)*vh2(xz,y)
            dvr(xz,y)=dvr(xz,y)-bcdr(xz,1)*dvh(xz,y)-
     &           bcdr(xz,2)*dvh2(xz,y)
            dvi(xz,y)=dvi(xz,y)-bcdi(xz,1)*dvh(xz,y)-
     &           bcdi(xz,2)*dvh2(xz,y)
         end do
      end do
c
c     If we use the tau method or cim without icorr, getting
c     final phr/phi from phr/phi, phh and cvr/cvi can cause
c     numerical instabilities, specially for a small ny.
c     (Continuity is not fullfilled.)
c     Therefore phr/phi are obtained from Laplacian of normal
c     velocity by definition if we are using the tau method.
c
      if (.not.cim) then
c
c     Chebyshev-tau method
c
         call dcheb(dvr,u3r(1,1,2),ny,nxz,nxz)
         call dcheb(phr,dvr,ny,nxz,nxz)
         call dcheb(dvi,u3i(1,1,2),ny,nxz,nxz)
         call dcheb(phi,dvi,ny,nxz,nxz)
         do y=1,ny
            do xz=1,nxz
               phr(xz,y)=phr(xz,y)-k2(xz)*u3r(xz,y,2)
               phi(xz,y)=phi(xz,y)-k2(xz)*u3i(xz,y,2)
            end do
         end do
      else
c
c     Chebyshev-integration method
c
         do xz=1,nxz
            do is=0,1
               do y=1+is,ny,2
                  c1=phh(xz,y)
                  phr(xz,y)=phr(xz,y)+cvr(xz,is+1)*c1
                  phi(xz,y)=phi(xz,y)+cvi(xz,is+1)*c1
                  phh2(xz,y)=phh2(xz,y)+cvh2(xz,is+1)*c1
                  phh(xz,y)=c1*cvh(xz,1+is)
               end do
            end do
            do y=1,ny
               phr(xz,y)=phr(xz,y)-bcdr(xz,1)*phh(xz,y)-
     &              bcdr(xz,2)*phh2(xz,y)
               phi(xz,y)=phi(xz,y)-bcdi(xz,1)*phh(xz,y)-
     &              bcdi(xz,2)*phh2(xz,y)
            end do
         end do
      end if

      do y=1,ny
         do xz=1,nxz
c
c     Find u,w,omx,omz from dvdy,omy,phi,domydy (this is loop 2800)
c
            u3r(xz,y,1)=(-alfa(xz)*dvi(xz,y)+
     &           bbeta(xz)*om3i(xz,y,2))*k2i(xz)
            u3i(xz,y,1)=(alfa(xz)*dvr(xz,y)-
     &           bbeta(xz)*om3r(xz,y,2))*k2i(xz)
            u3r(xz,y,3)=(-alfa(xz)*om3i(xz,y,2)-
     &           bbeta(xz)*dvi(xz,y))*k2i(xz)
            u3i(xz,y,3)=(alfa(xz)*om3r(xz,y,2)+
     &           bbeta(xz)*dvr(xz,y))*k2i(xz)
            om3r(xz,y,3)=(-bbeta(xz)*domyi(xz,y)+
     &           alfa(xz)*phi(xz,y))*k2i(xz)
            om3i(xz,y,3)=(bbeta(xz)*domyr(xz,y)-
     &           alfa(xz)*phr(xz,y))*k2i(xz)
            om3r(xz,y,1)=(-alfa(xz)*domyi(xz,y)-
     &           bbeta(xz)*phi(xz,y))*k2i(xz)
            om3i(xz,y,1)=(alfa(xz)*domyr(xz,y)+
     &           bbeta(xz)*phr(xz,y))*k2i(xz)
         end do
      end do
c
c     u,w for k=0
c
      if (zb.eq.1) then
         bc0(1,1)=.5*(u0upp+u0low)
         bc0(2,1)=.5*(u0upp-u0low)
         bc0(1,2)=.5*(w0upp+w0low)
         bc0(2,2)=.5*(w0upp-w0low)
c
c     Wall roughness and oscillation
c
         if ((wbci.eq.-1.or.wbci.eq.3.or.wbci.eq.4.or.wbci.eq.5) 
     &        .and.(my_node.eq.0)) then
            bc0(1,1)=.5*(u0upp+wallur(1,1))
            bc0(2,1)=.5*(u0upp-wallur(1,1))
            bc0(1,2)=.5*(w0upp+wallwr(1,1))
            bc0(2,2)=.5*(w0upp-wallwr(1,1))
         end if

         bc0(1,3)=1.
         bc0(2,3)=1.
c
c     Set the scalar boundary conditions (constants)
c
         do ith=1,scalar
            if (tbc(ith).eq.0) then
c
c     Dirichlet at both lower and upper boundary
c
               bc0(1,4+3*(ith-1)) =
     &              1./2.*(theta0_upp(ith)+theta0_low(ith))
               bc0(2,4+3*(ith-1)) =
     &              1./2.*(theta0_upp(ith)-theta0_low(ith))
c
c     In case of jet with scalar, set inhomogeneous Dirichlet b.c. at
c     lower wall
c
               if (wbci.eq.-2) then
                  bc0(1,4+3*(ith-1)) =
     &                 1./2.*(theta0_upp(ith)+wallthr(1,1,ith))
                  bc0(2,4+3*(ith-1)) =
     &                 1./2.*(theta0_upp(ith)-wallthr(1,1,ith))
               end if

            else if (tbc(ith).eq.1) then
c
c     Derivative at -1 (lower boundary), Dirichlet at upper boundary
c
               bc0(1,4+3*(ith-1)) =
     &              1./2.*(theta0_upp(ith)+0)
               bc0(2,4+3*(ith-1)) =
     &              1./2.*(theta0_upp(ith)-0)
               bc0(1,5+3*(ith-1)) =  1./2.   !  t(-1) = 1
               bc0(2,5+3*(ith-1)) = -1./2 .  !  t(+1) = 0
            else if (tbc(ith).eq.2) then
c
c     Derivative at 1 (upper boundary), Dirichlet at lower boundary
c
               bc0(1,4+3*(ith-1)) =
     &              1./2.*(0+theta0_low(ith))
               bc0(2,4+3*(ith-1)) =
     &              1./2.*(0-theta0_low(ith))
               bc0(1,5+3*(ith-1)) =  1./2.   !  t(-1) = 0
               bc0(2,5+3*(ith-1)) =  1./2.   !  t(+1) = 1
            else if (tbc(ith).eq.3) then
c
c     Derivative at -1,+1 (both boundaries)
c
               bc0(1,4+3*(ith-1)) =  0.      !  t(-1) = 0
               bc0(2,4+3*(ith-1)) =  0.      !  t(+1) = 0
               bc0(1,5+3*(ith-1)) =  1./2.   !  t(-1) = 1
               bc0(2,5+3*(ith-1)) = -1./2.   !  t(+1) = 0
               bc0(1,6+3*(ith-1)) =  1./2.   !  t(-1) = 0
               bc0(2,6+3*(ith-1)) =  1./2.   !  t(+1) = 1
            end if
         end do
c
c     Set up homogeneous rhs for the homogenous solutions
c
         do y=1,ny
            fuw(y,3)=0.0
         end do
         do ith=1,scalar
            do y=1,ny
               fuw(y,5+3*(ith-1))=0.0
               fuw(y,6+3*(ith-1))=0.0
            end do
         end do

         if (cim) then
c
c     Chebyshev-integration method
c     Solve for d2u0dy2, d2w0dy2, d2uwhomd2y
c
            call setmat(q0,c0,d0,e0,lam2,ny,1,1,cim)
            call trid(d2uw,cuw,q0,c0,d0,e0,fuw,bc0,1,ny+2,3,1,ny,cim)
c
c     Integrate once for du0dy and dw0dy and duwhomdy
c
            do i=1,3
               call icheb(duw(1,i),d2uw(1,i),cuw(2,i),ny,1,1)
               if (icorr) call corrd(duw(1,i),ny,1,1)
            end do
c
c     Integrate to find u0,w0
c
            if (ibc.eq.0.or.ibc.eq.100) then
               call icheb(uw(1,1),duw(1,1),cuw(1,1),ny,1,1)
               call icheb(uw(1,2),duw(1,2),cuw(1,2),ny,1,1)
            end if

            if (cflux) then
c
c     Constant mass flux not implemented for CIM, but
c     see e.q. in an old cha version.
c
               call stopnow(656456)
            end if

            do ith=1,scalar
               call setmat(q0,c0,d0,e0,lam3(1,ith),ny,1,1,cim)
               call trid(d2uw(1,3+ith),cuw(1,4+3*(ith-1)),q0,c0,d0,e0
     $              ,fuw(1,4+3*(ith-1)),bc0(1,4+3*(ith-1)),1,ny+2,1,1,ny
     $              ,cim)
               call icheb(duw(1,4+3*(ith-1)),d2uw(1,3+ith),cuw(2,4+3
     $              *(ith-1)),ny,1,1)
               if (icorr) call corrd(duw(1,4+3*(ith-1)),ny,1,1)
c
c     For boundary conditions...
c
               call icheb(uw(1,4+3*(ith-1)),duw(1,4+3*(ith-1)),cuw(1,4+3
     $              *(ith-1)),ny,1,1)
            end do

         else
c
c     Get the right-handed sides for the Chebyshev-tau method
c     Solve for u0, w0, uwhom
c
            call setmacrhs(uw,fuw,c0,d0,e0,lam2,1,ny,3,1,ny)
            call trid(uw,cuw,q0,c0,d0,e0,fuw,bc0,1,ny,3,1,ny,cim)

            if (cflux) then
c
c     For constant mass flux, integrate over the channel
c     and then correct the pressure gradient and integrate the
c     equation for wavenumber zero again.
c     This is a roundabout way to solve for the pressure gradient
c     for constant flux rather than implement a new linear solver.
c
c     Note that the direction of the forcing is given by ipdir.
c
c     Get the present mass flux
c
               mflux1 = 0.
               do y=1,ny,2
                  mflux1=mflux1-uw(y,ipdir)*2./real((y-1)**2-1)
               end do
               if (ipdir.eq.1) then
c                  mflux1 = mflux1 - (u0low+u0upp)
               else 
c                  mflux1 = mflux1 - (w0low+w0upp)
               end if
c
c     Double the pressure gradient (i.e. add it once more to fuw(1,1))
c
               fuw(1,ipdir) = fuw(1,ipdir)+2.*re*px
c
c     Solve for u0
c
               call setmacrhs(uw(1,ipdir),fuw(1,ipdir),
     &              c0,d0,e0,lam2,1,ny,1,1,ny)
               call trid(uw(1,ipdir),cuw(1,ipdir),q0,c0,d0,e0,
     &              fuw(1,ipdir),bc0(1,ipdir),1,ny,1,1,ny,cim)
c
c     and integrate this new mass flux
c
               mflux2 = 0.
               do y=1,ny,2
                  mflux2=mflux2-uw(y,ipdir)*2./real((y-1)**2-1)
               end do
               if (ipdir.eq.1) then
c                  mflux2 = mflux2 - (u0low+u0upp)
               else 
c                  mflux2 = mflux2 - (w0low+w0upp)
               end if
c
c     Reset pressure gradient (i.e. subtract both contributions)
c
               fuw(1,ipdir) = fuw(1,ipdir)-4.*re*px
c
c     And calculate the corrected pressure gradient
c
               px=px*(1.+(mflux-mflux1)/(mflux2-mflux1))
c
c     Solve for the final u0
c
               fuw(1,ipdir)=fuw(1,ipdir)+2.*re*px

               call setmacrhs(uw(1,ipdir),fuw(1,ipdir),
     &              c0,d0,e0,lam2,1,ny,1,1,ny)
               call trid(uw(1,ipdir),cuw(1,ipdir),q0,c0,d0,e0,
     &              fuw(1,ipdir),bc0(1,ipdir),1,ny,1,1,ny,cim)

            end if

c
c     To check final mass flux
c
c            mflux2 = 0.
c            do y=1,ny,2
c               mflux2=mflux2-uw(y,1)*2./real((y-1)**2-1)
c            end do
c            write(*,*) 'MFLUX before ',mflux1,' end ',mflux2,
c     &           ' target ',mflux
c     
c     Differentiate u0, w0, uwhom
c
            do i=1,3
               call dcheb(duw(1,i),uw(1,i),ny,1,1)
            end do
            cuw(1,1)=uw(1,1)
            cuw(1,2)=uw(1,2)
            cuw(1,3)=uw(1,3)

            do ith=1,scalar

               if (tbc(ith).lt.3) endi=2
               if (tbc(ith).eq.3) endi=3

               call setmacrhs(uw(1,4+3*(ith-1)),fuw(1,4+3*(ith-1)),c0,d0
     $              ,e0,lam3(1,ith),1,ny,endi,1,ny)
               call trid(uw(1,4+3*(ith-1)),cuw(1,4+3*(ith-1)),q0,c0,d0
     $              ,e0,fuw(1,4+3*(ith-1)),bc0(1,4+3*(ith-1)),1,ny,endi
     $              ,1,ny,cim)
c
c     Solving two systems
c
               do i=4,4+endi-1
                  call dcheb(duw(1,i+3*(ith-1)),uw(1,i+3*(ith-1)),ny,1,1
     $                 )
                  cuw(1,i+3*(ith-1))=uw(1,i+3*(ith-1))
               end do
            end do

         end if



         if (ibc.eq.0.or.ibc.eq.100) then
c
c     Nothing to do since no conditions on derivatives
c
         else if (ibc.eq.150) then
c
c     If one wants u(alpha=0,beta=0)=fixed then one needs to comment
c     out the "first lines" (see !comment this!;
c     then the vorticity is not exactly zero.
c     If one leaves both lines active, one gets omega_z(alpha=0,beta=0)=0,
c     but not constant velocity u(alpha=0,beta=0).
c
c     Now get the derivative in the freestream at y=+1
c
            do i=1,3
               duw1(i)=0.0
               do y=1,ny
                  duw1(i)=duw1(i)+duw(y,i)
               end do
            end do
c
c     Superimpose to get the prescribed value of the dudy and dwdy at y=1
c
            cuw(1,1)=cuw(1,1)-((duw1(1)-du0upp)/duw1(3))*cuw(1,3)!comment this!
            cuw(1,2)=cuw(1,2)-( duw1(2)        /duw1(3))*cuw(1,3)
            do y=1,ny
               duw(y,1)=duw(y,1)-((duw1(1)-du0upp)/duw1(3))*duw(y,3)!comment this!
               duw(y,2)=duw(y,2)-(duw1(2)/duw1(3))*duw(y,3)
            end do
            call icheb(uw(1,1),duw(1,1),cuw(1,1),ny,1,1)!comment this!
            call icheb(uw(1,2),duw(1,2),cuw(1,2),ny,1,1)

         else if (ibc.eq.160) then
c
c     w is a Dirichlet condition, therefore the "second lines"
c     are commented out.
c
c     Now get the derivative in the freestream at y=+1
c
            do i=1,3
               duw1(i)=0.0
               do y=1,ny
                  duw1(i)=duw1(i)+duw(y,i)
               end do
            end do
c
c     Superimpose to get the prescribed value of the dudy and dwdy at y=1
c
            cuw(1,1)=cuw(1,1)-((duw1(1)-du0upp)/duw1(3))*cuw(1,3)
c            cuw(1,2)=cuw(1,2)-( duw1(2)        /duw1(3))*cuw(1,3)
            do y=1,ny
               duw(y,1)=duw(y,1)-((duw1(1)-du0upp)/duw1(3))*duw(y,3)
c               duw(y,2)=duw(y,2)-(duw1(2)/duw1(3))*duw(y,3)
            end do
            call icheb(uw(1,1),duw(1,1),cuw(1,1),ny,1,1)
c            call icheb(uw(1,2),duw(1,2),cuw(1,2),ny,1,1)

         else
c
c     Now get the derivative in the freestream at y=+1
c
            do i=1,3
               duw1(i)=0.0
               do y=1,ny
                  duw1(i)=duw1(i)+duw(y,i)
               end do
            end do
c
c     Superimpose to get the prescribed value of the dudy and dwdy at y=1
c
            cuw(1,1)=cuw(1,1)-((duw1(1)-du0upp)/duw1(3))*cuw(1,3)
            cuw(1,2)=cuw(1,2)-( duw1(2)        /duw1(3))*cuw(1,3)
            do y=1,ny
               duw(y,1)=duw(y,1)-((duw1(1)-du0upp)/duw1(3))*duw(y,3)
               duw(y,2)=duw(y,2)-( duw1(2)        /duw1(3))*duw(y,3)
            end do
            call icheb(uw(1,1),duw(1,1),cuw(1,1),ny,1,1)
            call icheb(uw(1,2),duw(1,2),cuw(1,2),ny,1,1)

         end if

         do ith=1,scalar
            if (tbc(ith).gt.0) then
c
c     Now get the derivative at the wall at y=-1 or y=1
c
               if (tbc(ith).eq.2) then
                  do i=4,5
                     duw1(i+3*(ith-1))=0.0
                     do y=1,ny
c
c     Derivative 1
c
                        duw1(i+3*(ith-1))=duw1(i+3*(ith-1))+duw(y,i+3
     $                       *(ith-1))
                     end do
                  end do
               else if (tbc(ith).eq.1) then
                  do i=4,5
                     duw1(i+3*(ith-1))=0.0
                     do y=1,ny
c
c     Derivative -1
c
                        duw1(i+3*(ith-1))=duw1(i+3*(ith-1))+duw(y,i+3
     $                       *(ith-1))*(-1)**(y-1)
                     end do
                  end do
               else if (tbc(ith).eq.3) then
                  do i=4,6
                     duw1(i+3*(ith-1))=0.0
                     duw_1(i+3*(ith-1))=0.0
                     do y=1,ny
c
c     Derivative 1
c
                        duw1(i+3*(ith-1))=duw1(i+3*(ith-1))+duw(y,i+3
     $                       *(ith-1))
c
c     Derivative -1
c
                        duw_1(i+3*(ith-1))=duw_1(i+3*(ith-1))+duw(y,i+3
     $                       *(ith-1))*(-1)**(y-1)
                     end do
                  end do
               end if
c
c     Superimpose to get the prescribed value of the dthetady at y=-1/y=1
c
               if (tbc(ith).eq.1) then
c
c     Only one derivative...
c
                  cuw(1,4+3*(ith-1))=cuw(1,4+3*(ith-1))-
     &                 ((duw1(4+3*(ith-1))-dtheta0_low(ith))/
     &                 duw1(5+3*(ith-1)))*cuw(1,5+3*(ith-1))
                  do y=1,ny
                     duw(y,4+3*(ith-1))=duw(y,4+3*(ith-1))-
     &                    ((duw1(4+3*(ith-1))-dtheta0_low(ith))/
     &                    duw1(5+3*(ith-1)))*duw(y,5+3*(ith-1))
                  end do
                  call icheb(uw(1,4+3*(ith-1)),duw(1,4+3*(ith-1)),
     &                 cuw(1,4+3*(ith-1)),ny,1,1)
               else if (tbc(ith).eq.2) then
c
c     Only one derivative...
c
                  cuw(1,4+3*(ith-1))=cuw(1,4+3*(ith-1))-
     &                 ((duw1(4+3*(ith-1))-dtheta0_upp(ith))/
     &                 duw1(5+3*(ith-1)))*cuw(1,5+3*(ith-1))
                  do y=1,ny
                     duw(y,4+3*(ith-1))=duw(y,4+3*(ith-1))-
     &                    ((duw1(4+3*(ith-1))-dtheta0_upp(ith))/
     &                    duw1(5+3*(ith-1)))*duw(y,5+3*(ith-1))
                  end do
                  call icheb(uw(1,4+3*(ith-1)),duw(1,4+3*(ith-1)),
     &                 cuw(1,4+3*(ith-1)),ny,1,1)
               else if (tbc(ith).eq.3) then
c
c     Both derivatives...
c
                  a11=duw1(5+3*(ith-1))
                  a12=duw1(6+3*(ith-1))
                  a21=duw_1(5+3*(ith-1))
                  a22=duw_1(6+3*(ith-1))
                  b1=dtheta0_upp(ith)-duw1(4+3*(ith-1))
                  b2=dtheta0_low(ith)-duw_1(4+3*(ith-1))
c
c     Solve and get AA and BB
c
                  AA=(a22*b1-b2*a12)/(a22*a11-a12*a21)
                  BB=-(a21*b1-b2*a11)/(a22*a11-a12*a21)
                  cuw(1,4+3*(ith-1))=cuw(1,4+3*(ith-1))+
     &                 AA*cuw(1,5+3*(ith-1))+BB*cuw(1,6+3*(ith-1))
                  do y=1,ny
                     duw(y,4+3*(ith-1))=duw(y,4+3*(ith-1))+
     &                    AA*duw(y,5+3*(ith-1))+BB*duw(y,6+3*(ith-1))
                  end do
                  call icheb(uw(1,4+3*(ith-1)),duw(1,4+3*(ith-1)),
     &                 cuw(1,4+3*(ith-1)),ny,1,1)
               end if
            end if
         end do
c
c     End test section
c
c     omx0=dw0dy, omz0=-du0dy,
c     Note that the imaginary parts are already zeroed in loop 2800
c
         do y=1,ny
            om3r(1,y,1)= duw(y,2)
            om3r(1,y,3)=-duw(y,1)
            u3r (1,y,1) =  uw(y,1)
            u3r (1,y,3) =  uw(y,2)
c
c     Set the v velocity to zero
c
            u3r (1,y,2) = 0.
            u3i (1,y,2) = 0.
         end do
c
c     Since du/dx and dw/dz are both zero for zero wavenumbers the
c     vertical velocity mean is constant, set it to the mean component
c     of the blowing (or suction)
c
         u3r(1,1,2) = -vsuc
         if (wbci.gt.0 .or. wbci.eq.-2) then
            u3r (1,1,2) = wallvr(1,1)
         endif

         do ith=1,scalar
            do y=1,ny
c
c     Set new solution
c
               th3r(1,y,ith) =  uw(y,4+3*(ith-1))
               dthr(1,y,ith) = duw(y,4+3*(ith-1))
c
c     Probably not needed...
c
               th3i(1,y,ith) = 0.
               dthi(1,y,ith) = 0.
            end do
c
c     dthr,dthi is completely finished with zero mode
c
         end do
c
c     End of k=0 special
c
      end if
c
c     Calculate partial rhs for next step
c
      c1 = 2.*re*(1./(an+bn) + 1./(anp1+bnp1))
      c2 = 2.*re*bnp1/(anp1 + bnp1)
      do ith=1,scalar
         c3(ith) = c1*pr(ith)
         c4(ith) = c2*pr(ith)
      end do
      if (.not.cim) then
c
c     Chebyshev-tau method
c
         c1=2.*re/(anp1+bnp1)
         do ith=1,scalar
            c3(ith)=c1*pr(ith)
         end do
      end if
c
c     First put in all zero order terms
c
      if (zb.eq.1) then
         if (.not.cim) then
            call dcheb(duw(1,1), uw(1,1),ny,1,1)
            call dcheb(fuw(1,1),duw(1,1),ny,1,1)
            call dcheb(duw(1,2), uw(1,2),ny,1,1)
            call dcheb(fuw(1,2),duw(1,2),ny,1,1)

            do ith=1,scalar
               call dcheb(fuw(1,4+3*(ith-1)),duw(1,4+3*(ith-1)),ny,1,1)
            end do
         end if
         do y=1,ny
            puw(y,1)=-fuw(y,1)-c1*u3r(1,y,1)-c2*huw(y,1)
            puw(y,2)=-fuw(y,2)-c1*u3r(1,y,3)-c2*huw(y,2)
            do ith=1,scalar
               puw(y,2+ith)=-fuw(y,4+3*(ith-1))-c3(ith)*th3r(1,y,ith)
     $              -c4(ith)*huw(y,2+ith)
            end do
         end do
      end if
      if (cim) then
         do y=1,ny
            do xz=1,nxz
               pvr(xz,y)=-fvr(xz,y)-c1*phr(xz,y)-c2*hvr(xz,y)
               pvi(xz,y)=-fvi(xz,y)-c1*phi(xz,y)-c2*hvi(xz,y)
               pomyr(xz,y)=-fomyr(xz,y)-c1*om3r(xz,y,2)-c2*homyr(xz,y)
               pomyi(xz,y)=-fomyi(xz,y)-c1*om3i(xz,y,2)-c2*homyi(xz,y)
               do ith=1,scalar
                  pthr(xz,y,ith)=-fthr(xz,y,1+(ith-1)*2)
     &                 -c3(ith)*th3r(xz,y,ith)
     &                 -c4(ith)*hthr(xz,y,ith)
                  pthi(xz,y,ith)=-fthi(xz,y,1+(ith-1)*2)
     &                 -c3(ith)*th3i(xz,y,ith)
     &                 -c4(ith)*hthi(xz,y,ith)
               end do
            end do
         end do
      else
         call dcheb(w3,phr,ny,nxz,nxz)
         call dcheb(fvr,w3,ny,nxz,nxz)
         call dcheb(w3,phi,ny,nxz,nxz)
         call dcheb(fvi,w3,ny,nxz,nxz)
         call dcheb(fomyr,domyr,ny,nxz,nxz)
         call dcheb(fomyi,domyi,ny,nxz,nxz)
         do ith=1,scalar
c
c     fthr,fthi is d2theta
c
            call dcheb(fthr(1,1,1+(ith-1)*2),
     &           dthr(1,1,ith),ny,nxz,nxz)
            call dcheb(fthi(1,1,1+(ith-1)*2),
     &           dthi(1,1,ith),ny,nxz,nxz)
         end do

         do xz=1,nxz
            do ith=1,scalar
               k3(xz,ith)=k2(xz)-c3(ith)
            end do
            k2(xz)=k2(xz)-c1
         end do
         do y=1,ny
            do xz=1,nxz
               pvr(xz,y)=-fvr(xz,y)+k2(xz)*phr(xz,y)-c2*hvr(xz,y)
               pvi(xz,y)=-fvi(xz,y)+k2(xz)*phi(xz,y)-c2*hvi(xz,y)
               pomyr(xz,y)=-fomyr(xz,y)+k2(xz)*om3r(xz,y,2)-
     &              c2*homyr(xz,y)
               pomyi(xz,y)=-fomyi(xz,y)+k2(xz)*om3i(xz,y,2)-
     &              c2*homyi(xz,y)
            end do
         end do
         do ith=1,scalar
            do y=1,ny
               do xz=1,nxz
                  pthr(xz,y,ith)=-fthr(xz,y,1+(ith-1)*2)
     &                 +k3(xz,ith)*th3r(xz,y
     &                 ,ith)-c4(ith)*hthr(xz,y,ith)
                  pthi(xz,y,ith)=-fthi(xz,y,1+(ith-1)*2)
     &                 +k3(xz,ith)*th3i(xz,y
     &                 ,ith)-c4(ith)*hthi(xz,y,ith)
               end do
            end do
         end do
c
c     pthr,pthi is finished
c
      end if
c
c     Pad the dealiasing region in y (if dealiasing in y used)
c
      if (ny+1.le.nyp) then
         do y=ny+1,max(nyp,ny+1)
            do xz=1,nxz
               u3r (xz,y,1) = 0.0
               u3i (xz,y,1) = 0.0
               u3r (xz,y,2) = 0.0
               u3i (xz,y,2) = 0.0
               u3r (xz,y,3) = 0.0
               u3i (xz,y,3) = 0.0
               om3r(xz,y,1) = 0.0
               om3i(xz,y,1) = 0.0
               om3r(xz,y,3) = 0.0
               om3i(xz,y,3) = 0.0
            end do
         end do
         do ith=1,scalar
            do y=ny+1,max(nyp,ny+1)
               do xz=1,nxz
                  th3r(xz,y,ith)=0.0
                  th3i(xz,y,ith)=0.0
               end do
            end do
         end do
      end if
c
c     Transform back to Fourier/Fourier/physical space
c
      call vchbb(u3r(1,1,1) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(u3i(1,1,1) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(u3r(1,1,2) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(u3i(1,1,2) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(u3r(1,1,3) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(u3i(1,1,3) ,w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3r(1,1,1),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3i(1,1,1),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3r(1,1,3),w3,nyp,nxz,nxz,1,prey)
      call vchbb(om3i(1,1,3),w3,nyp,nxz,nxz,1,prey)

      do ith=1,scalar
         call vchbb(th3r(1,1,ith),w3,nyp,nxz,nxz,1,prey)
         call vchbb(th3i(1,1,ith),w3,nyp,nxz,nxz,1,prey)
         call vchbb(dthr(1,1,ith),w3,nyp,nxz,nxz,1,prey)
         call vchbb(dthi(1,1,ith),w3,nyp,nxz,nxz,1,prey)
      end do
c
c     Store the boxes in the velocities, the vorticities and 6,7
c
      call putxy(u3r(1,1,1) ,u3i(1,1,1) ,zbp,1,ur,ui)
      call putxy(u3r(1,1,2) ,u3i(1,1,2) ,zbp,2,ur,ui)
      call putxy(u3r(1,1,3) ,u3i(1,1,3) ,zbp,3,ur,ui)
      call putxy(om3r(1,1,1),om3i(1,1,1),zbp,4,ur,ui)
      call putxy(om3r(1,1,3),om3i(1,1,3),zbp,5,ur,ui)
      call putxy(pvr        ,pvi        ,zbp,6,ur,ui)
      call putxy(pomyr      ,pomyi      ,zbp,7,ur,ui)
      do ith=1,scalar
c
c     Store the boxes of theta, dtheta/dy, prhs in 8-10
c
         call putxy(th3r(1,1,ith),th3i(1,1,ith) ,zbp,
     &        8 +pressure+3*(ith-1),ur,ui)
         call putxy(dthr(1,1,ith),dthi(1,1,ith) ,zbp,
     &        9 +pressure+3*(ith-1),ur,ui)
         call putxy(pthr(1,1,ith),pthi(1,1,ith) ,zbp,
     &        10+pressure+3*(ith-1),ur,ui)
      end do
c
c     LES
c
      if (iles.eq.1.or.iles.eq.3) then
         do ll=1,3
c     Get gu3
            do y=1,nyp
               do xz=1,mbz*nx/2
                  gu3r(xz,y) = u3r(xz,y,ll)
                  gu3i(xz,y) = u3i(xz,y,ll)
               end do
            end do

            if (ihighorder.eq.0) then
c
c     Only two-dimensional filtering
c
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = u3r(xz,y,ll) -
     &                       lpfxz(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll) -
     &                       lpfxz(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.1) then
c
c     Do the multiple filtering via the iterative way
c
               do i=1,iord+1
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
                           ggu3r(xz,y) = filtxz(x,z+zb-1)*gu3r(xz,y)
                           ggu3i(xz,y) = filtxz(x,z+zb-1)*gu3i(xz,y)
                        end do
                     end do
                  end do
                  call lpf_y(ggu3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call lpf_y(ggu3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  do y=1,nyp
                     do xz=1,mbz*nx/2
                        gu3r(xz,y) = gu3r(xz,y)-ggu3r(xz,y)
                        gu3i(xz,y) = gu3i(xz,y)-ggu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.2) then
c
c     Do the multiple filtering via the old way
c     High-pass filter velocities in y
c
               do i=1,iord+1
                  call hpf_y(gu3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call hpf_y(gu3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
               end do
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Store low-pass filtered velocities
c
                        gu3r(xz,y) = u3r(xz,y,ll)-gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll)-gu3i(xz,y)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = u3r(xz,y,ll) -
     &                       lpfxz(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll) -
     &                       lpfxz(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else
               call stopnow(447)
            end if
c
c     Put into gur,gui
c
            call putxy(gu3r,gu3i,zbp,ll,gur,gui)
         end do
      end if


      if (iles.eq.3) then
c
c     Save also the HPF vorticity
c
         do ll=1,2
c
c     Get gu3
c
            do y=1,nyp
               do xz=1,mbz*nx/2
                  gu3r(xz,y) = om3r(xz,y,ll*2-1)
                  gu3i(xz,y) = om3i(xz,y,ll*2-1)
               end do
            end do

            if (ihighorder.eq.0) then
c
c     Only two-dimensional filtering
c
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1) -
     &                       lpfxz(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1) -
     &                       lpfxz(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.1) then
c
c     Do it the iterative way
c
               do i=1,iord+1
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
                           ggu3r(xz,y) = filtxz(x,z+zb-1)*gu3r(xz,y)
                           ggu3i(xz,y) = filtxz(x,z+zb-1)*gu3i(xz,y)
                        end do
                     end do
                  end do
                  call lpf_y(ggu3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call lpf_y(ggu3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  do y=1,nyp
                     do xz=1,mbz*nx/2
                        gu3r(xz,y) = gu3r(xz,y)-ggu3r(xz,y)
                        gu3i(xz,y) = gu3i(xz,y)-ggu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.2) then
c
c     Do it the old way
c     High-pass filter velocities in y
c
               do i=1,iord+1
                  call hpf_y(gu3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call hpf_y(gu3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
               end do
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Store low-pass filtered velocities
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1)-gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1)-gu3i(xz,y)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1) -
     &                       lpfxz(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1) -
     &                       lpfxz(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else
               call stopnow(447)
            end if
c
c     Put into gur,gui
c
            call putxy(gu3r,gu3i,zbp,ll+3,gur,gui)
         end do
      end if





      if (iles.eq.3.and.cs.eq.0) then
c
c     Save the velocities and vorticities with secondary HPF
c
         do ll=1,3
c
c     Get gu3
c
            do y=1,nyp
               do xz=1,mbz*nx/2
                  gu3r(xz,y) = u3r(xz,y,ll)
                  gu3i(xz,y) = u3i(xz,y,ll)
               end do
            end do

            if (ihighorder.eq.0) then
c
c     Only two-dimensional filtering
c
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = u3r(xz,y,ll) -
     &                       lpfxz2(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll) -
     &                       lpfxz2(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.1) then
c
c     Do it the iterative way
c
               do i=1,iord+1
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
                           ggu3r(xz,y) = filtxz2(x,z+zb-1)*gu3r(xz,y)
                           ggu3i(xz,y) = filtxz2(x,z+zb-1)*gu3i(xz,y)
                        end do
                     end do
                  end do
                  call lpf_y(ggu3r,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  call lpf_y(ggu3i,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  do y=1,nyp
                     do xz=1,mbz*nx/2
                        gu3r(xz,y) = gu3r(xz,y)-ggu3r(xz,y)
                        gu3i(xz,y) = gu3i(xz,y)-ggu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.2) then
c
c     Do it the old way
c     High-pass filter velocities in y
c
               do i=1,iord+1
                  call hpf_y(gu3r,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  call hpf_y(gu3i,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
               end do
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Store low-pass filtered velocities
c
                        gu3r(xz,y) = u3r(xz,y,ll)-gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll)-gu3i(xz,y)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = u3r(xz,y,ll) -
     &                       lpfxz2(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = u3i(xz,y,ll) -
     &                       lpfxz2(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else
               call stopnow(447)
            end if
c
c     Put into gur,gui
c
            call putxy(gu3r,gu3i,zbp,ll,taur,taui)
         end do

         do ll=1,2
c
c     Get gu3
c
            do y=1,nyp
               do xz=1,mbz*nx/2
                  gu3r(xz,y) = om3r(xz,y,ll*2-1)
                  gu3i(xz,y) = om3i(xz,y,ll*2-1)
               end do
            end do

            if (ihighorder.eq.0) then
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1) -
     &                       lpfxz2(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1) -
     &                       lpfxz2(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.1) then
c
c     Do it the iterative way
c
               do i=1,iord+1
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
                           ggu3r(xz,y) = filtxz2(x,z+zb-1)*gu3r(xz,y)
                           ggu3i(xz,y) = filtxz2(x,z+zb-1)*gu3i(xz,y)
                        end do
                     end do
                  end do
                  call lpf_y(ggu3r,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  call lpf_y(ggu3i,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  do y=1,nyp
                     do xz=1,mbz*nx/2
                        gu3r(xz,y) = gu3r(xz,y)-ggu3r(xz,y)
                        gu3i(xz,y) = gu3i(xz,y)-ggu3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.2) then
c
c     Do it the old way
c     High-pass filter velocities in y
c
               do i=1,iord+1
                  call hpf_y(gu3r,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
                  call hpf_y(gu3i,gewy12,gewy22,nx/2,nyp,mbz,
     &                 diags2,gplane)
               end do
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Store low-pass filtered velocities
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1)-gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1)-gu3i(xz,y)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gu3r(xz,y) = om3r(xz,y,ll*2-1) -
     &                       lpfxz2(x,z+zb-1)*gu3r(xz,y)
                        gu3i(xz,y) = om3i(xz,y,ll*2-1) -
     &                       lpfxz2(x,z+zb-1)*gu3i(xz,y)
                     end do
                  end do
               end do
            else
               call stopnow(447)
            end if
c
c     Put into gur,gui
c
            call putxy(gu3r,gu3i,zbp,ll+3,taur,taui)
         end do
      end if
c
c     LES for the scalar
c
      if (iles.eq.1.or.iles.eq.3) then
         do ith=1,scalar
c     Get gth3
            do y=1,nyp
               do xz=1,mbz*nx/2
                  gth3r(xz,y) = th3r(xz,y,ith)
                  gth3i(xz,y) = th3i(xz,y,ith)
               end do
            end do

            if (ihighorder.eq.0) then
c
c     Only two-dimensional filtering
c
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gth3r(xz,y) = th3r(xz,y,ith) -
     &                       lpfxz(x,z+zb-1)*gth3r(xz,y)
                        gth3i(xz,y) = th3i(xz,y,ith) -
     &                       lpfxz(x,z+zb-1)*gth3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.1) then
c
c     Do the multiple filtering via the iterative way
c
               do i=1,iord+1
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
                           ggth3r(xz,y) = filtxz(x,z+zb-1)*gth3r(xz,y)
                           ggth3i(xz,y) = filtxz(x,z+zb-1)*gth3i(xz,y)
                        end do
                     end do
                  end do
                  call lpf_y(ggth3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call lpf_y(ggth3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  do y=1,nyp
                     do xz=1,mbz*nx/2
                        gth3r(xz,y) = gth3r(xz,y)-ggth3r(xz,y)
                        gth3i(xz,y) = gth3i(xz,y)-ggth3i(xz,y)
                     end do
                  end do
               end do
            else if (ihighorder.eq.2) then
c
c     Do the multiple filtering via the old way
c     High-pass filter velocities in y
c
               do i=1,iord+1
                  call hpf_y(gth3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
                  call hpf_y(gth3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                 diags,gplane)
               end do
               do y=1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        xz = x+nx/2*(z-1)
c
c     Store low-pass filtered velocities
c
                        gth3r(xz,y) = th3r(xz,y,ith)-gth3r(xz,y)
                        gth3i(xz,y) = th3i(xz,y,ith)-gth3i(xz,y)
c
c     Low-pass filter velocities in x/z and store high-pass filter
c
                        gth3r(xz,y) = th3r(xz,y,ith) -
     &                       lpfxz(x,z+zb-1)*gth3r(xz,y)
                        gth3i(xz,y) = th3i(xz,y,ith) -
     &                       lpfxz(x,z+zb-1)*gth3i(xz,y)
                     end do
                  end do
               end do
            else
               call stopnow(447)
            end if
c
c     Put into gthr,gthi
c
            if (iles.eq.3) then
               call putxy(gth3r,gth3i,zbp,ith,gthr,gthi)
            else
               call putxy(gth3r,gth3i,zbp,ith,gsr,gsi)
            end if
         end do

         if (iles.eq.3) then
c
c     Save also the HPF dtheta/dy
c
            do ith=1,scalar
c
c     Get gth3
c
               do y=1,nyp
                  do xz=1,mbz*nx/2
                     gth3r(xz,y) = dthr(xz,y,ith)
                     gth3i(xz,y) = dthi(xz,y,ith)
                  end do
               end do

               if (ihighorder.eq.0) then
c
c     Only two-dimensional filtering
c
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
c
c     Low-pass filter dtheta/dy in x/z and store high-pass filter
c
                           gth3r(xz,y) = dthr(xz,y,ith) -
     &                          lpfxz(x,z+zb-1)*gth3r(xz,y)
                           gth3i(xz,y) = dthi(xz,y,ith) -
     &                          lpfxz(x,z+zb-1)*gth3i(xz,y)
                        end do
                     end do
                  end do
               else if (ihighorder.eq.1) then
c
c     Do it the iterative way
c
                  do i=1,iord+1
                     do y=1,nyp
                        do z=1,mbz
                           do x=1,nx/2
                              xz = x+nx/2*(z-1)
                              ggth3r(xz,y) = filtxz(x,z+zb-1)
     &                             *gth3r(xz,y)
                              ggth3i(xz,y) = filtxz(x,z+zb-1)
     &                             *gth3i(xz,y)
                           end do
                        end do
                     end do
                     call lpf_y(ggth3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                    diags,gplane)
                     call lpf_y(ggth3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                    diags,gplane)
                     do y=1,nyp
                        do xz=1,mbz*nx/2
                           gth3r(xz,y) = gth3r(xz,y)-ggth3r(xz,y)
                           gth3i(xz,y) = gth3i(xz,y)-ggth3i(xz,y)
                        end do
                     end do
                  end do
               else if (ihighorder.eq.2) then
c
c     Do it the old way
c     High-pass filter scalar in y
c
                  do i=1,iord+1
                     call hpf_y(gth3r,gewy1,gewy2,nx/2,nyp,mbz,
     &                    diags,gplane)
                     call hpf_y(gth3i,gewy1,gewy2,nx/2,nyp,mbz,
     &                    diags,gplane)
                  end do
                  do y=1,nyp
                     do z=1,mbz
                        do x=1,nx/2
                           xz = x+nx/2*(z-1)
c
c     Store low-pass filtered scalar
c
                           gth3r(xz,y) = dthr(xz,y,ith)-gth3r(xz,y)
                           gth3i(xz,y) = dthi(xz,y,ith)-gth3i(xz,y)
c
c     Low-pass filter scalar in x/z and store high-pass filter
c
                           gth3r(xz,y) = dthr(xz,y,ith) -
     &                          lpfxz(x,z+zb-1)*gth3r(xz,y)
                           gth3i(xz,y) = dthi(xz,y,ith) -
     &                          lpfxz(x,z+zb-1)*gth3i(xz,y)
                        end do
                     end do
                  end do
               else
                  call stopnow(447)
               end if
c
c     Put into gthr,gthi
c
               call putxy(gth3r,gth3i,zbp,ith+scalar,gthr,gthi)
            end do
         end if
      end if

      end subroutine linearbl
