c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine les(re,xs,zs,yb,ur,ui,
     &     u0low,w0low,prexn,prezn,presn,prean,
     &     alfa,beta,u2r,u2i,om2i,om2r,wr,wi,
     &     du2r,du2i,
     &     sr,si,lr,li,
     &     mr,mi,br,bi,
     &     sabsr,sabsi,cr,ci,sdr,sdi,
     &     my_node,realg1,realg2,deltaxyz2,
     &     taur,taui,gthr,gthi,
     &     do_stat,xys,ybp,dtn,ineg,iles,cs,
     &     gur,gui,gu2r,gu2i,gom2r,gom2i,gdu2r,gdu2i,
     &     th2r,th2i,lsr,lsi,prt,gall,gridy)

c     Dynamic Smagorinsky model (Germano/Lilly)
c     Computation of the Smagorinsky coefficient C=C_S^2 and
c     the SGS stresses tau_ij
c     Details:
c     - 2D spectral cutoff filter (w_c=pi/2)
c     - clipping at C>0
c     - averaging in z
c
c     Differences to TRANSIT implementation:
c     BLES:    M already includes delta prior to filtering
c              (not relevant for spectral cutoff filter)
c     BLES:    for M, Sij is filtered instead of u
c              (not relevant for spectral cutoff filter)
c     BLES:    delta~^2 = 2^(2*2/3) * delta^2 = 2.519 * delta^2
c     TRANSIT: delta~^2 = 2^(2)     * delta^2 = 4     * delta^2

      implicit none
      include 'par.f'

      real,external :: step
      real gridy(nyp)
      real ydamp
      real xs,zs,re
      integer yb
      real u0low,w0low
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)

      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real du2r((nxp/2+1)*mby,nzd,3,3),du2i((nxp/2+1)*mby,nzd,3,3)

      real gu2r((nxp/2+1)*mby,nzd,3),gu2i((nxp/2+1)*mby,nzd,3)
      real gom2r((nxp/2+1)*mby,nzd,3),gom2i((nxp/2+1)*mby,nzd,3)
      real gdu2r((nxp/2+1)*mby,nzd,3,3),gdu2i((nxp/2+1)*mby,nzd,3,3)

      real th2r((nxp/2+1)*mby,nzd,4,scalar)
      real th2i((nxp/2+1)*mby,nzd,4,scalar)

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real alfa(nx/2*mbz),beta(nz)
      real deltaxyz2(nyp)

      logical gall

      integer nxy,npl
      integer x,y,xy,z,i,j,ll,zz

c     MPI
      integer my_node
      integer realg1,realg2

c     LES
      real taur (memnx,memny,memnz,6+3*scalar)
      real taui (memnx,memny,memnz,6+3*scalar)
      real gur  (memnx,memny,memnz,5)
      real gui  (memnx,memny,memnz,5)
      real gthr (memnx,memny,memnz,2*scalar)
      real gthi (memnx,memny,memnz,2*scalar)
      real sr   ((nxp/2+1)*mby,nzd,6)
      real si   ((nxp/2+1)*mby,nzd,6)
      real lr   ((nxp/2+1)*mby,nzd,6)
      real li   ((nxp/2+1)*mby,nzd,6)
      real lsr   ((nxp/2+1)*mby,nzd,3,scalar)
      real lsi   ((nxp/2+1)*mby,nzd,3,scalar)
      real mr   ((nxp/2+1)*mby,nzd,6)
      real mi   ((nxp/2+1)*mby,nzd,6)
      real br   ((nxp/2+1)*mby,nzd,6)
      real bi   ((nxp/2+1)*mby,nzd,6)
      real sabsr((nxp/2+1)*mby,nzd)
      real sabsi((nxp/2+1)*mby,nzd)
      real cr   ((nxp/2+1)*mby,nzd)
      real ci   ((nxp/2+1)*mby,nzd)
      real cthr ((nxp/2+1)*mby,nzd)
      real cthi ((nxp/2+1)*mby,nzd)
      real sdr  ((nxp/2+1)*mby,nzd)
      real sdi  ((nxp/2+1)*mby,nzd)

      real cnumr,cnumi,cdenr,cdeni,gd2,dr,di,d1r,d1i
      real tsabsr,tsabsi,tsdr,tsdi
      integer ineg,iles,iwale
      real cs

c     statistics
      logical do_stat
      real dtn,c,c1
      integer ybp
      real xys(nx,nyp/nproc+1,nxys)

c     scalar
      integer ith
      real ktr,kti,prt


      npl=min(mby,nyp-yb+1)
      nxy=(nxp/2+1)*npl-(nxp/2+1-nx/2)

      du2r = 0.0
      du2i = 0.0
      u2r  = 0.0
      u2i  = 0.0
      om2r = 0.0
      om2i = 0.0
      th2r = 0.0
      th2i = 0.0

      gdu2r = 0.0
      gdu2i = 0.0
      gu2r  = 0.0
      gu2i  = 0.0
      gom2r = 0.0
      gom2i = 0.0

c
c     If using the WALE model, set iwale to 1
c
      iwale = 0
c
c     Find velocities and vorticities (computed in linearbl/prhs)
c     only needed for classical/dynamic Smagorinsky and
c     dynamic HPF Smagorinsky
c
      if (nproc.eq.1) then
         do i=1,3
            call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,0,ur,ui)
         end do
         call getxz(om2r(1,1,1),om2i(1,1,1),yb,4,0,ur,ui)
         call getxz(om2r(1,1,3),om2i(1,1,3),yb,5,0,ur,ui)
c
c     Get the scalar and its y-derivative
c
         do ith=1,scalar
            call getxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
     &           yb,8+pressure+3*(ith-1),0,ur,ui)
            call getxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
     &           yb,9+pressure+3*(ith-1),0,ur,ui)
         end do
      else
#ifdef MPI
         do i=1,3
            call getpxz(u2r(1,1,i),u2i(1,1,i),yb,i,0,ur,ui,
     &           realg1,realg2,my_node)
         end do
         call getpxz(om2r(1,1,1),om2i(1,1,1),yb,4,0,ur,ui,
     &        realg1,realg2,my_node)
         call getpxz(om2r(1,1,3),om2i(1,1,3),yb,5,0,ur,ui,
     &        realg1,realg2,my_node)
         do ith=1,scalar
            call getpxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
     &           yb,8+pressure+3*(ith-1),0,ur,ui,
     &           realg1,realg2,my_node)
            call getpxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
     &           yb,9+pressure+3*(ith-1),0,ur,ui,
     &           realg1,realg2,my_node)
         end do
#endif
      end if

      if (iles.eq.3) then
c     find HPF velocities and vorticities (computed in linearbl/prhs)
c     used for HPF models
         if (nproc.eq.1) then
            do i=1,3
               call getxz(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,gur,gui)
            end do
            call getxz(gom2r(1,1,1),gom2i(1,1,1),yb,4,0,gur,gui)
            call getxz(gom2r(1,1,3),gom2i(1,1,3),yb,5,0,gur,gui)
c
c     Get the scalar and its y-derivative
c
            do ith=1,scalar
               call getxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
     &              yb,ith,0,gthr,gthi)
               call getxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
     &              yb,ith+scalar,0,gthr,gthi)
            end do
         else
#ifdef MPI
            do i=1,3
               call getpxz(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,gur,gui,
     &              realg1,realg2,my_node)
            end do
            call getpxz(gom2r(1,1,1),gom2i(1,1,1),yb,4,0,gur,gui,
     &           realg1,realg2,my_node)
            call getpxz(gom2r(1,1,3),gom2i(1,1,3),yb,5,0,gur,gui,
     &           realg1,realg2,my_node)
            do ith=1,scalar
               call getpxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
     &              yb,ith,0,gthr,gthi,
     &              realg1,realg2,my_node)
               call getpxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
     &              yb,ith+scalar,0,gthr,gthi,
     &              realg1,realg2,my_node)
            end do
#endif
         end if
      end if

      if (iles.eq.2) then
c     calculate missing vorticity omega_y for the classical quantities
c     only needed for classical Smagorinsky
         do z=1,nzc
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  xy=x+(y-yb)*(nxp/2+1)
                  om2r(xy,z,2)=-beta(z)*u2i(xy,z,1)+alfa(x)*u2i(xy,z,3)
                  om2i(xy,z,2)=beta(z)*u2r(xy,z,1)-alfa(x)*u2r(xy,z,3)
               end do
               om2r(nx/2+1+(y-yb)*(nxp/2+1),z,2)=0.0
            end do
         end do
c
c     Calculate x and z derivatives for scalar field
c
         do ith=1,scalar
            do z=1,nzc
               do y=yb,min(nyp,yb+mby-1)
                  do x=1,nx/2
                     xy=x+(y-yb)*(nxp/2+1)
c
c     Calculate x derivative
c
                     th2r(xy,z,2,ith)=-alfa(x)*
     &                    th2i(xy,z,1,ith)
                     th2i(xy,z,2,ith)= alfa(x)*
     &                    th2r(xy,z,1,ith)
c
c     Calculate z derivative
c
                     th2r(xy,z,4,ith)=-beta(z)*
     &                    th2i(xy,z,1,ith)
                     th2i(xy,z,4,ith)= beta(z)*
     &                    th2r(xy,z,1,ith)
                  end do
                  th2r(nx/2+1+(y-yb)*(nxp/2+1),z,2,ith)=0.0
                  th2r(nx/2+1+(y-yb)*(nxp/2+1),z,4,ith)=0.0
               end do
            end do
         end do

         if (nfzsym.eq.0) then
c
c     Zero oddball (z)
c
            do xy=1,(nxp/2+1)*npl
               om2r(xy,nz/2+1,2)=0.0
               om2i(xy,nz/2+1,2)=0.0
            end do
            do ith=1,scalar
               do xy=1,(nxp/2+1)*npl
                  th2r(xy,nz/2+1,2,ith)=0.0
                  th2i(xy,nz/2+1,2,ith)=0.0
                  th2r(xy,nz/2+1,4,ith)=0.0
                  th2i(xy,nz/2+1,4,ith)=0.0
               end do
            end do
         end if
c
c     Subtract lower wall shifting velocity
c
         do y=1,npl
            xy=1+(y-1)*(nxp/2+1)
            u2r(xy,1,1)=u2r(xy,1,1)-u0low
            u2r(xy,1,3)=u2r(xy,1,3)-w0low
         end do
c
c     Calculate all velocity derivatives
c     horizontal derivatives by multiplikation by alfa or beta
c
         do i=1,3
            do z=1,nzc
               do y=1,npl
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     du2r(xy,z,i,1)=-alfa(x)*u2i(xy,z,i)
                     du2i(xy,z,i,1)= alfa(x)*u2r(xy,z,i)
                     du2r(xy,z,i,3)=-beta(z)*u2i(xy,z,i)
                     du2i(xy,z,i,3)= beta(z)*u2r(xy,z,i)
                  end do
               end do
            end do
         end do
c
c     The wall normal derivatives can be found from :
c     omx  = dwdy - dvdz
c     omz  = dvdx - dudy
c     i.e. :
c     dudy = dvdx - omz
c     dvdy =-dudx - dwdz   (from continuity)
c     dwdy = omx  + dvdz
c
         do z=1,nzc
            do y=1,npl
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  du2r(xy,z,1,2)=-alfa(x)*u2i(xy,z,2) - om2r(xy,z,3)
                  du2i(xy,z,1,2)= alfa(x)*u2r(xy,z,2) - om2i(xy,z,3)
                  du2r(xy,z,2,2)= alfa(x)*u2i(xy,z,1) +
     &                 beta(z)*u2i(xy,z,3)
                  du2i(xy,z,2,2)=-alfa(x)*u2r(xy,z,1) -
     &                 beta(z)*u2r(xy,z,3)
                  du2r(xy,z,3,2)= om2r(xy,z,1)        -
     &                 beta(z)*u2i(xy,z,2)
                  du2i(xy,z,3,2)= om2i(xy,z,1)        +
     &                 beta(z)*u2r(xy,z,2)
               end do
            end do
         end do
c
c     Set oddball modes to zero
c
         do j=1,3
            do i=1,3
               do z=1,nzc
                  do y=1,npl
                     du2r(nx/2+1+(y-1)*(nxp/2+1),z,i,j) = 0.0
                  end do
               end do
            end do
         end do
c
c     Transform to physical space
c
         do i=1,3
c
c     No shifting for LES quantities
c
            if (nfzsym.eq.0) then
               call vcfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
               call vcfftb(om2r(1,1,i),om2i(1,1,i),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
               do j=1,3
                  call vcfftb(du2r(1,1,i,j),du2i(1,1,i,j),wr,wi,nz,
     &                 nxy,(nxp/2+1)*mby,1,prezn)
               end do
            else
               if (i.le.2) then
                  call vcffts(u2r(1,1,i),u2i(1,1,i),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(om2r(1,2,i),om2i(1,2,i),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(du2r(1,1,i,1),du2i(1,1,i,1),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcffts(du2r(1,1,i,2),du2i(1,1,i,2),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(du2r(1,2,i,3),du2i(1,2,i,3),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               else
                  call vcftab(u2r(1,2,i),u2i(1,2,i),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(om2r(1,1,i),om2i(1,1,i),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(du2r(1,2,i,1),du2i(1,2,i,1),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcftab(du2r(1,2,i,2),du2i(1,2,i,2),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(du2r(1,1,i,3),du2i(1,1,i,3),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               end if
            end if
            call vrfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            call vrfftb(om2r(1,1,i),om2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            do j=1,3
               call vrfftb(du2r(1,1,i,j),du2i(1,1,i,j),wr,wi,
     &              nx,nzpc*mby,1,nxp/2+1,prexn)
            end do
         end do
c
c     Maybe some shifting (as above) could be done here as well.
c
         do ith=1,scalar
            do i=1,4
               call vcfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
     &              ,wr,wi,nz,nxy,(nxp/2+1)*mby,1,prezn)
               call vrfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
     &              ,wr,wi,nx,nzpc*mby,1,nxp/2+1,prexn)
            end do
         end do

      else if (iles.eq.3) then
c
c     Calculate missing vorticity omega_y for the HPF quantities
c     only needed for HPF Smagorinsky
c
c     Calculate missing vorticity omega_y
c
         do z=1,nzc
            do y=yb,min(nyp,yb+mby-1)
               do x=1,nx/2
                  xy=x+(y-yb)*(nxp/2+1)
                  gom2r(xy,z,2)=-beta(z)*gu2i(xy,z,1)+
     &                 alfa(x)*gu2i(xy,z,3)
                  gom2i(xy,z,2)=beta(z)*gu2r(xy,z,1)-
     &                 alfa(x)*gu2r(xy,z,3)
               end do
               gom2r(nx/2+1+(y-yb)*(nxp/2+1),z,2)=0.0
            end do
         end do

c
c     Calculate x and z derivatives for scalar field
c
         do ith=1,scalar
            do z=1,nzc
               do y=yb,min(nyp,yb+mby-1)
                  do x=1,nx/2
                     xy=x+(y-yb)*(nxp/2+1)
c
c     Calculate x derivative
c
                     th2r(xy,z,2,ith)=-alfa(x)*
     &                    th2i(xy,z,1,ith)
                     th2i(xy,z,2,ith)= alfa(x)*
     &                    th2r(xy,z,1,ith)
c
c     Calculate z derivative
c
                     th2r(xy,z,4,ith)=-beta(z)*
     &                    th2i(xy,z,1,ith)
                     th2i(xy,z,4,ith)= beta(z)*
     &                    th2r(xy,z,1,ith)
                  end do
                  th2r(nx/2+1+(y-yb)*(nxp/2+1),z,2,ith)=0.0
                  th2r(nx/2+1+(y-yb)*(nxp/2+1),z,4,ith)=0.0
               end do
            end do
         end do

         if (nfzsym.eq.0) then
c
c     Zero oddball (z)
c
            do xy=1,(nxp/2+1)*npl
               gom2r(xy,nz/2+1,2)=0.0
               gom2i(xy,nz/2+1,2)=0.0
            end do
            do ith=1,scalar
               do xy=1,(nxp/2+1)*npl
                  th2r(xy,nz/2+1,2,ith)=0.0
                  th2i(xy,nz/2+1,2,ith)=0.0
                  th2r(xy,nz/2+1,4,ith)=0.0
                  th2i(xy,nz/2+1,4,ith)=0.0
               end do
            end do
         end if
c
c     Subtract lower wall shifting velocity
c     not needed since here we consider HPF quantities
c         do y=1,npl
c            xy=1+(y-1)*(nxp/2+1)
c            gu2r(xy,1,1)=gu2r(xy,1,1)-u0low
c            gu2r(xy,1,3)=gu2r(xy,1,3)-w0low
c         end do
c
c     Calculate all velocity derivatives
c     horizontal derivatives by multiplikation by alfa or beta
         do i=1,3
            do z=1,nzc
               do y=1,npl
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     gdu2r(xy,z,i,1)=-alfa(x)*gu2i(xy,z,i)
                     gdu2i(xy,z,i,1)= alfa(x)*gu2r(xy,z,i)
                     gdu2r(xy,z,i,3)=-beta(z)*gu2i(xy,z,i)
                     gdu2i(xy,z,i,3)= beta(z)*gu2r(xy,z,i)
                  end do
               end do
            end do
         end do
c
c     The wall normal derivatives can be found from :
c     omx  = dwdy - dvdz
c     omz  = dvdx - dudy
c     i.e. :
c     dudy = dvdx - omz
c     dvdy =-dudx - dwdz   (from continuity)
c     dwdy = omx  + dvdz
c
         do z=1,nzc
            do y=1,npl
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  gdu2r(xy,z,1,2)=-alfa(x)*gu2i(xy,z,2) - gom2r(xy,z,3)
                  gdu2i(xy,z,1,2)= alfa(x)*gu2r(xy,z,2) - gom2i(xy,z,3)
                  gdu2r(xy,z,2,2)= alfa(x)*gu2i(xy,z,1) +
     &                 beta(z)*gu2i(xy,z,3)
                  gdu2i(xy,z,2,2)=-alfa(x)*gu2r(xy,z,1) -
     &                 beta(z)*gu2r(xy,z,3)
                  gdu2r(xy,z,3,2)= gom2r(xy,z,1)        -
     &                 beta(z)*gu2i(xy,z,2)
                  gdu2i(xy,z,3,2)= gom2i(xy,z,1)        +
     &                 beta(z)*gu2r(xy,z,2)
               end do
            end do
         end do
c
c     Set oddball modes to zero
c
         do j=1,3
            do i=1,3
               do z=1,nzc
                  do y=1,npl
                     gdu2r(nx/2+1+(y-1)*(nxp/2+1),z,i,j) = 0.0
                  end do
               end do
            end do
         end do
c
c     Transform to physical space
c
         do i=1,3
c
c     no shifting for LES quantities
c
            if (nfzsym.eq.0) then
               call vcfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
               call vcfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
               call vcfftb(gom2r(1,1,i),gom2i(1,1,i),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
               do j=1,3
                  call vcfftb(gdu2r(1,1,i,j),gdu2i(1,1,i,j),wr,wi,nz,
     &                 nxy,(nxp/2+1)*mby,1,prezn)
               end do
            else
               if (i.le.2) then
                  call vcffts(gu2r(1,1,i),gu2i(1,1,i),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcffts(u2r(1,1,i),u2i(1,1,i),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(gom2r(1,2,i),gom2i(1,2,i),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(gdu2r(1,1,i,1),gdu2i(1,1,i,1),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcffts(gdu2r(1,1,i,2),gdu2i(1,1,i,2),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(gdu2r(1,2,i,3),gdu2i(1,2,i,3),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               else
                  call vcftab(gu2r(1,2,i),gu2i(1,2,i),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcftab(u2r(1,2,i),u2i(1,2,i),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(gom2r(1,1,i),gom2i(1,1,i),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  call vcftab(gdu2r(1,2,i,1),gdu2i(1,2,i,1),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcftab(gdu2r(1,2,i,2),gdu2i(1,2,i,2),
     &                 wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  call vcffts(gdu2r(1,1,i,3),gdu2i(1,1,i,3),wr,
     &                 nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               end if
            end if
            call vrfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            call vrfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            call vrfftb(gom2r(1,1,i),gom2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            do j=1,3
               call vrfftb(gdu2r(1,1,i,j),gdu2i(1,1,i,j),wr,wi,
     &              nx,nzpc*mby,1,nxp/2+1,prexn)
            end do
         end do
c
c     Maybe some shifting (as above) could be done here as well.
c
         do ith=1,scalar
            do i=1,4
               call vcfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
     &              ,wr,wi,nz,nxy,(nxp/2+1)*mby,1,prezn)
               call vrfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
     &              ,wr,wi,nx,nzpc*mby,1,nxp/2+1,prexn)
            end do
         end do
      end if

c----------------------------------------------------------------------c
c     we can start computing the LES quantities                        c
c----------------------------------------------------------------------c

c     here we have:
c     u2r,u2i:   velocities  (physical space)
c     om2r,om2i: vorticities (physical space)
c     du2r,du2i: derivatives (physical space)
c     eta:       wall-normal coordinate (1..-1)
c     gridy:     wall-normal coordinate (physical space)
c
c     du2r etc. have to be multiplied by dstar to get the right scaling
c
c     y=nyp is at the wall

      if (scalar.gt.0) then
         if (iles.eq.0.or.iles.eq.1.or.iles.eq.2.or.iles.eq.3) then
         else
            write(*,*) 'NO SGS MODEL FOR SCALAR'
            write(*,*) iles,cs,scalar
            stop
         end if
      end if

      if (iles.eq.2) then
c
c     CLASSICAL EDDY-VISCOSITY MODELS
c
c     Strain rate Sij
c     Note: The scalar gradient is already in th2r,th2i
c
         do y=1,npl
            do z=1,nzc+nfzsym
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)

                  sr(xy,z,1) =     du2r(xy,z,1,1)
                  sr(xy,z,2) = .5*(du2r(xy,z,2,1)+du2r(xy,z,1,2))
                  sr(xy,z,3) = .5*(du2r(xy,z,3,1)+du2r(xy,z,1,3))
                  sr(xy,z,4) =     du2r(xy,z,2,2)
                  sr(xy,z,5) = .5*(du2r(xy,z,2,3)+du2r(xy,z,3,2))
                  sr(xy,z,6) =     du2r(xy,z,3,3)

                  si(xy,z,1) =     du2i(xy,z,1,1)
                  si(xy,z,2) = .5*(du2i(xy,z,2,1)+du2i(xy,z,1,2))
                  si(xy,z,3) = .5*(du2i(xy,z,3,1)+du2i(xy,z,1,3))
                  si(xy,z,4) =     du2i(xy,z,2,2)
                  si(xy,z,5) = .5*(du2i(xy,z,2,3)+du2i(xy,z,3,2))
                  si(xy,z,6) =     du2i(xy,z,3,3)

                  sabsr(xy,z) = sqrt(
     &                 2.*(sr(xy,z,1)**2 + sr(xy,z,4)**2 +
     &                 sr(xy,z,6)**2) +
     &                 4.*(sr(xy,z,2)**2 + sr(xy,z,3)**2 +
     &                 sr(xy,z,5)**2))
                  sabsi(xy,z) = sqrt(
     &                 2.*(si(xy,z,1)**2 + si(xy,z,4)**2 +
     &                 si(xy,z,6)**2) +
     &                 4.*(si(xy,z,2)**2 + si(xy,z,3)**2 +
     &                 si(xy,z,5)**2))

                  if (iwale.eq.1) then
                     sdr(xy,z)=sqrt((
     &                    du2r(xy,z,1,1)*du2r(xy,z,1,2) +
     &                    du2r(xy,z,1,2)*du2r(xy,z,2,2) +
     &                    du2r(xy,z,1,3)*du2r(xy,z,3,2) +
     &                    du2r(xy,z,2,1)*du2r(xy,z,1,1) +
     &                    du2r(xy,z,2,2)*du2r(xy,z,2,1) +
     &                    du2r(xy,z,2,3)*du2r(xy,z,3,1))**2/2. +(
     &                    du2r(xy,z,1,1)*du2r(xy,z,1,3) +
     &                    du2r(xy,z,1,2)*du2r(xy,z,2,3) +
     &                    du2r(xy,z,1,3)*du2r(xy,z,3,3) +
     &                    du2r(xy,z,3,1)*du2r(xy,z,1,1) +
     &                    du2r(xy,z,3,2)*du2r(xy,z,2,1) +
     &                    du2r(xy,z,3,3)*du2r(xy,z,3,1))**2/2. +(
     &                    du2r(xy,z,2,1)*du2r(xy,z,1,3) +
     &                    du2r(xy,z,2,2)*du2r(xy,z,2,3) +
     &                    du2r(xy,z,2,3)*du2r(xy,z,3,3) +
     &                    du2r(xy,z,3,1)*du2r(xy,z,1,2) +
     &                    du2r(xy,z,3,2)*du2r(xy,z,2,2) +
     &                    du2r(xy,z,3,3)*du2r(xy,z,3,2))**2/2.)

                     sdi(xy,z)=sqrt((
     &                    du2i(xy,z,1,1)*du2i(xy,z,1,2) +
     &                    du2i(xy,z,1,2)*du2i(xy,z,2,2) +
     &                    du2i(xy,z,1,3)*du2i(xy,z,3,2) +
     &                    du2i(xy,z,2,1)*du2i(xy,z,1,1) +
     &                    du2i(xy,z,2,2)*du2i(xy,z,2,1) +
     &                    du2i(xy,z,2,3)*du2i(xy,z,3,1))**2/2. +(
     &                    du2i(xy,z,1,1)*du2i(xy,z,1,3) +
     &                    du2i(xy,z,1,2)*du2i(xy,z,2,3) +
     &                    du2i(xy,z,1,3)*du2i(xy,z,3,3) +
     &                    du2i(xy,z,3,1)*du2i(xy,z,1,1) +
     &                    du2i(xy,z,3,2)*du2i(xy,z,2,1) +
     &                    du2i(xy,z,3,3)*du2i(xy,z,3,1))**2/2. +(
     &                    du2i(xy,z,2,1)*du2i(xy,z,1,3) +
     &                    du2i(xy,z,2,2)*du2i(xy,z,2,3) +
     &                    du2i(xy,z,2,3)*du2i(xy,z,3,3) +
     &                    du2i(xy,z,3,1)*du2i(xy,z,1,2) +
     &                    du2i(xy,z,3,2)*du2i(xy,z,2,2) +
     &                    du2i(xy,z,3,3)*du2i(xy,z,3,2))**2/2.)
                  end if
               end do
            end do
         end do

         if (cs.eq.0) then
c     dynamic procedure (Germano/Lilly)
c     Parts of M_ij and L_ij
            do y=1,npl
               do z=1,nzc+nfzsym
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)

                     lr(xy,z,1) = -u2r(xy,z,1)*u2r(xy,z,1)
                     lr(xy,z,2) = -u2r(xy,z,2)*u2r(xy,z,1)
                     lr(xy,z,3) = -u2r(xy,z,3)*u2r(xy,z,1)
                     lr(xy,z,4) = -u2r(xy,z,2)*u2r(xy,z,2)
                     lr(xy,z,5) = -u2r(xy,z,3)*u2r(xy,z,2)
                     lr(xy,z,6) = -u2r(xy,z,3)*u2r(xy,z,3)

                     li(xy,z,1) = -u2i(xy,z,1)*u2i(xy,z,1)
                     li(xy,z,2) = -u2i(xy,z,2)*u2i(xy,z,1)
                     li(xy,z,3) = -u2i(xy,z,3)*u2i(xy,z,1)
                     li(xy,z,4) = -u2i(xy,z,2)*u2i(xy,z,2)
                     li(xy,z,5) = -u2i(xy,z,3)*u2i(xy,z,2)
                     li(xy,z,6) = -u2i(xy,z,3)*u2i(xy,z,3)

                     if (iwale.eq.0) then
                        dr = sabsr(xy,z)*deltaxyz2(yb+y-1)
                        di = sabsi(xy,z)*deltaxyz2(yb+y-1)
                     else
                        dr = sdr(xy,z)**3*deltaxyz2(yb+y-1)
     &                       /(sabsr(xy,z)**5+sqrt(sdr(xy,z))**5)
                        di = sdi(xy,z)**3*deltaxyz2(yb+y-1)
     &                       /(sabsi(xy,z)**5+sqrt(sdi(xy,z))**5)
                     end if

                     mr(xy,z,1) = -sr(xy,z,1)*dr
                     mr(xy,z,2) = -sr(xy,z,2)*dr
                     mr(xy,z,3) = -sr(xy,z,3)*dr
                     mr(xy,z,4) = -sr(xy,z,4)*dr
                     mr(xy,z,5) = -sr(xy,z,5)*dr
                     mr(xy,z,6) = -sr(xy,z,6)*dr

                     mi(xy,z,1) = -si(xy,z,1)*di
                     mi(xy,z,2) = -si(xy,z,2)*di
                     mi(xy,z,3) = -si(xy,z,3)*di
                     mi(xy,z,4) = -si(xy,z,4)*di
                     mi(xy,z,5) = -si(xy,z,5)*di
                     mi(xy,z,6) = -si(xy,z,6)*di

                     br(xy,z,1) = sr(xy,z,1)
                     br(xy,z,2) = sr(xy,z,2)
                     br(xy,z,3) = sr(xy,z,3)
                     br(xy,z,4) = sr(xy,z,4)
                     br(xy,z,5) = sr(xy,z,5)
                     br(xy,z,6) = sr(xy,z,6)

                     bi(xy,z,1) = si(xy,z,1)
                     bi(xy,z,2) = si(xy,z,2)
                     bi(xy,z,3) = si(xy,z,3)
                     bi(xy,z,4) = si(xy,z,4)
                     bi(xy,z,5) = si(xy,z,5)
                     bi(xy,z,6) = si(xy,z,6)

                  end do
               end do
            end do
c
c     Filter M, L, B, u
c
            call filter_cutoff(mr ,mi ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(lr ,li ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(br ,bi ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(u2r,u2i,3,npl,nxy,prexn,prezn,wr,wi)
            if (iwale.eq.1) then
               call filter_cutoff(du2r,du2i,9,npl,nxy,prexn,prezn,wr,wi)
            end if
c
c     Computation of L_ij and M_ij
c
            do y=1,npl
               gd2 = 4.**(2./3.)*deltaxyz2(yb+y-1)
               do z=1,nzc+nfzsym
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)

                     lr(xy,z,1) = lr(xy,z,1) + u2r(xy,z,1)*u2r(xy,z,1)
                     lr(xy,z,2) = lr(xy,z,2) + u2r(xy,z,2)*u2r(xy,z,1)
                     lr(xy,z,3) = lr(xy,z,3) + u2r(xy,z,3)*u2r(xy,z,1)
                     lr(xy,z,4) = lr(xy,z,4) + u2r(xy,z,2)*u2r(xy,z,2)
                     lr(xy,z,5) = lr(xy,z,5) + u2r(xy,z,3)*u2r(xy,z,2)
                     lr(xy,z,6) = lr(xy,z,6) + u2r(xy,z,3)*u2r(xy,z,3)

                     li(xy,z,1) = li(xy,z,1) + u2i(xy,z,1)*u2i(xy,z,1)
                     li(xy,z,2) = li(xy,z,2) + u2i(xy,z,2)*u2i(xy,z,1)
                     li(xy,z,3) = li(xy,z,3) + u2i(xy,z,3)*u2i(xy,z,1)
                     li(xy,z,4) = li(xy,z,4) + u2i(xy,z,2)*u2i(xy,z,2)
                     li(xy,z,5) = li(xy,z,5) + u2i(xy,z,3)*u2i(xy,z,2)
                     li(xy,z,6) = li(xy,z,6) + u2i(xy,z,3)*u2i(xy,z,3)

                     tsabsr = sqrt(
     &                    2.*(br(xy,z,1)**2+br(xy,z,4)**2+
     &                    br(xy,z,6)**2)+
     &                    4.*(br(xy,z,2)**2+br(xy,z,3)**2+
     &                    br(xy,z,5)**2))
                     tsabsi = sqrt(
     &                    2.*(bi(xy,z,1)**2+bi(xy,z,4)**2+
     &                    bi(xy,z,6)**2)+
     &                    4.*(bi(xy,z,2)**2+bi(xy,z,3)**2+
     &                    bi(xy,z,5)**2))

                     if (iwale.eq.1) then
                        tsdr=sqrt((
     &                       du2r(xy,z,1,1)*du2r(xy,z,1,2) +
     &                       du2r(xy,z,1,2)*du2r(xy,z,2,2) +
     &                       du2r(xy,z,1,3)*du2r(xy,z,3,2) +
     &                       du2r(xy,z,2,1)*du2r(xy,z,1,1) +
     &                       du2r(xy,z,2,2)*du2r(xy,z,2,1) +
     &                       du2r(xy,z,2,3)*du2r(xy,z,3,1))**2/2. +(
     &                       du2r(xy,z,1,1)*du2r(xy,z,1,3) +
     &                       du2r(xy,z,1,2)*du2r(xy,z,2,3) +
     &                       du2r(xy,z,1,3)*du2r(xy,z,3,3) +
     &                       du2r(xy,z,3,1)*du2r(xy,z,1,1) +
     &                       du2r(xy,z,3,2)*du2r(xy,z,2,1) +
     &                       du2r(xy,z,3,3)*du2r(xy,z,3,1))**2/2. +(
     &                       du2r(xy,z,2,1)*du2r(xy,z,1,3) +
     &                       du2r(xy,z,2,2)*du2r(xy,z,2,3) +
     &                       du2r(xy,z,2,3)*du2r(xy,z,3,3) +
     &                       du2r(xy,z,3,1)*du2r(xy,z,1,2) +
     &                       du2r(xy,z,3,2)*du2r(xy,z,2,2) +
     &                       du2r(xy,z,3,3)*du2r(xy,z,3,2))**2/2.)

                        tsdi=sqrt((
     &                       du2i(xy,z,1,1)*du2i(xy,z,1,2) +
     &                       du2i(xy,z,1,2)*du2i(xy,z,2,2) +
     &                       du2i(xy,z,1,3)*du2i(xy,z,3,2) +
     &                       du2i(xy,z,2,1)*du2i(xy,z,1,1) +
     &                       du2i(xy,z,2,2)*du2i(xy,z,2,1) +
     &                       du2i(xy,z,2,3)*du2i(xy,z,3,1))**2/2. +(
     &                       du2i(xy,z,1,1)*du2i(xy,z,1,3) +
     &                       du2i(xy,z,1,2)*du2i(xy,z,2,3) +
     &                       du2i(xy,z,1,3)*du2i(xy,z,3,3) +
     &                       du2i(xy,z,3,1)*du2i(xy,z,1,1) +
     &                       du2i(xy,z,3,2)*du2i(xy,z,2,1) +
     &                       du2i(xy,z,3,3)*du2i(xy,z,3,1))**2/2. +(
     &                       du2i(xy,z,2,1)*du2i(xy,z,1,3) +
     &                       du2i(xy,z,2,2)*du2i(xy,z,2,3) +
     &                       du2i(xy,z,2,3)*du2i(xy,z,3,3) +
     &                       du2i(xy,z,3,1)*du2i(xy,z,1,2) +
     &                       du2i(xy,z,3,2)*du2i(xy,z,2,2) +
     &                       du2i(xy,z,3,3)*du2i(xy,z,3,2))**2/2.)
                     end if
                     if (iwale.eq.0) then
                        dr = gd2*tsabsr
                        di = gd2*tsabsi
                     else
                        dr = gd2*tsdr**3/(tsabsr**5+sqrt(tsdr)**5)
                        di = gd2*tsdi**3/(tsabsi**5+sqrt(tsdi)**5)
                     end if

                     mr(xy,z,1) = mr(xy,z,1) + br(xy,z,1)*dr
                     mr(xy,z,2) = mr(xy,z,2) + br(xy,z,2)*dr
                     mr(xy,z,3) = mr(xy,z,3) + br(xy,z,3)*dr
                     mr(xy,z,4) = mr(xy,z,4) + br(xy,z,4)*dr
                     mr(xy,z,5) = mr(xy,z,5) + br(xy,z,5)*dr
                     mr(xy,z,6) = mr(xy,z,6) + br(xy,z,6)*dr

                     mi(xy,z,1) = mi(xy,z,1) + bi(xy,z,1)*di
                     mi(xy,z,2) = mi(xy,z,2) + bi(xy,z,2)*di
                     mi(xy,z,3) = mi(xy,z,3) + bi(xy,z,3)*di
                     mi(xy,z,4) = mi(xy,z,4) + bi(xy,z,4)*di
                     mi(xy,z,5) = mi(xy,z,5) + bi(xy,z,5)*di
                     mi(xy,z,6) = mi(xy,z,6) + bi(xy,z,6)*di

                  end do
               end do
            end do
         end if
      else if (iles.eq.3) then
c
c     HPF models (dynamic or non-dynamic)
c

c
c     Strain rate Sij
c
         do y=1,npl
            do z=1,nzc+nfzsym
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)

                  sr(xy,z,1) =     gdu2r(xy,z,1,1)
                  sr(xy,z,2) = .5*(gdu2r(xy,z,2,1)+gdu2r(xy,z,1,2))
                  sr(xy,z,3) = .5*(gdu2r(xy,z,3,1)+gdu2r(xy,z,1,3))
                  sr(xy,z,4) =     gdu2r(xy,z,2,2)
                  sr(xy,z,5) = .5*(gdu2r(xy,z,2,3)+gdu2r(xy,z,3,2))
                  sr(xy,z,6) =     gdu2r(xy,z,3,3)

                  si(xy,z,1) =     gdu2i(xy,z,1,1)
                  si(xy,z,2) = .5*(gdu2i(xy,z,2,1)+gdu2i(xy,z,1,2))
                  si(xy,z,3) = .5*(gdu2i(xy,z,3,1)+gdu2i(xy,z,1,3))
                  si(xy,z,4) =     gdu2i(xy,z,2,2)
                  si(xy,z,5) = .5*(gdu2i(xy,z,2,3)+gdu2i(xy,z,3,2))
                  si(xy,z,6) =     gdu2i(xy,z,3,3)

                  sabsr(xy,z) = sqrt(
     &                 2.*(sr(xy,z,1)**2 + sr(xy,z,4)**2 +
     &                 sr(xy,z,6)**2) +
     &                 4.*(sr(xy,z,2)**2 + sr(xy,z,3)**2 +
     &                 sr(xy,z,5)**2))
                  sabsi(xy,z) = sqrt(
     &                 2.*(si(xy,z,1)**2 + si(xy,z,4)**2 +
     &                 si(xy,z,6)**2) +
     &                 4.*(si(xy,z,2)**2 + si(xy,z,3)**2 +
     &                 si(xy,z,5)**2))
c                  sabsr(xy,z) = sqrt(2*sr(xy,z,4)**2)
c                  sabsi(xy,z) = sqrt(2*si(xy,z,4)**2)
                  if (iwale.eq.1) then
                     sdr(xy,z)=sqrt((
     &                    gdu2r(xy,z,1,1)*gdu2r(xy,z,1,2) +
     &                    gdu2r(xy,z,1,2)*gdu2r(xy,z,2,2) +
     &                    gdu2r(xy,z,1,3)*gdu2r(xy,z,3,2) +
     &                    gdu2r(xy,z,2,1)*gdu2r(xy,z,1,1) +
     &                    gdu2r(xy,z,2,2)*gdu2r(xy,z,2,1) +
     &                    gdu2r(xy,z,2,3)*gdu2r(xy,z,3,1))**2/2. +(
     &                    gdu2r(xy,z,1,1)*gdu2r(xy,z,1,3) +
     &                    gdu2r(xy,z,1,2)*gdu2r(xy,z,2,3) +
     &                    gdu2r(xy,z,1,3)*gdu2r(xy,z,3,3) +
     &                    gdu2r(xy,z,3,1)*gdu2r(xy,z,1,1) +
     &                    gdu2r(xy,z,3,2)*gdu2r(xy,z,2,1) +
     &                    gdu2r(xy,z,3,3)*gdu2r(xy,z,3,1))**2/2. +(
     &                    gdu2r(xy,z,2,1)*gdu2r(xy,z,1,3) +
     &                    gdu2r(xy,z,2,2)*gdu2r(xy,z,2,3) +
     &                    gdu2r(xy,z,2,3)*gdu2r(xy,z,3,3) +
     &                    gdu2r(xy,z,3,1)*gdu2r(xy,z,1,2) +
     &                    gdu2r(xy,z,3,2)*gdu2r(xy,z,2,2) +
     &                    gdu2r(xy,z,3,3)*gdu2r(xy,z,3,2))**2/2.)

                     sdi(xy,z)=sqrt((
     &                    gdu2i(xy,z,1,1)*gdu2i(xy,z,1,2) +
     &                    gdu2i(xy,z,1,2)*gdu2i(xy,z,2,2) +
     &                    gdu2i(xy,z,1,3)*gdu2i(xy,z,3,2) +
     &                    gdu2i(xy,z,2,1)*gdu2i(xy,z,1,1) +
     &                    gdu2i(xy,z,2,2)*gdu2i(xy,z,2,1) +
     &                    gdu2i(xy,z,2,3)*gdu2i(xy,z,3,1))**2/2. +(
     &                    gdu2i(xy,z,1,1)*gdu2i(xy,z,1,3) +
     &                    gdu2i(xy,z,1,2)*gdu2i(xy,z,2,3) +
     &                    gdu2i(xy,z,1,3)*gdu2i(xy,z,3,3) +
     &                    gdu2i(xy,z,3,1)*gdu2i(xy,z,1,1) +
     &                    gdu2i(xy,z,3,2)*gdu2i(xy,z,2,1) +
     &                    gdu2i(xy,z,3,3)*gdu2i(xy,z,3,1))**2/2. +(
     &                    gdu2i(xy,z,2,1)*gdu2i(xy,z,1,3) +
     &                    gdu2i(xy,z,2,2)*gdu2i(xy,z,2,3) +
     &                    gdu2i(xy,z,2,3)*gdu2i(xy,z,3,3) +
     &                    gdu2i(xy,z,3,1)*gdu2i(xy,z,1,2) +
     &                    gdu2i(xy,z,3,2)*gdu2i(xy,z,2,2) +
     &                    gdu2i(xy,z,3,3)*gdu2i(xy,z,3,2))**2/2.)
                  end if

               end do
            end do
         end do

         if (cs.eq.0) then
c     dynamic procedure (Germano/Lilly)
c     Parts of M_ij and L_ij
            do y=1,npl
               do z=1,nzc+nfzsym
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)

                     lr(xy,z,1) = -u2r(xy,z,1)*u2r(xy,z,1)
                     lr(xy,z,2) = -u2r(xy,z,2)*u2r(xy,z,1)
                     lr(xy,z,3) = -u2r(xy,z,3)*u2r(xy,z,1)
                     lr(xy,z,4) = -u2r(xy,z,2)*u2r(xy,z,2)
                     lr(xy,z,5) = -u2r(xy,z,3)*u2r(xy,z,2)
                     lr(xy,z,6) = -u2r(xy,z,3)*u2r(xy,z,3)

                     li(xy,z,1) = -u2i(xy,z,1)*u2i(xy,z,1)
                     li(xy,z,2) = -u2i(xy,z,2)*u2i(xy,z,1)
                     li(xy,z,3) = -u2i(xy,z,3)*u2i(xy,z,1)
                     li(xy,z,4) = -u2i(xy,z,2)*u2i(xy,z,2)
                     li(xy,z,5) = -u2i(xy,z,3)*u2i(xy,z,2)
                     li(xy,z,6) = -u2i(xy,z,3)*u2i(xy,z,3)

                     if (iwale.eq.0) then
                        dr = sabsr(xy,z)*deltaxyz2(yb+y-1)
                        di = sabsi(xy,z)*deltaxyz2(yb+y-1)
                     else
                        dr = sdr(xy,z)**3*deltaxyz2(yb+y-1)
     &                       /(sabsr(xy,z)**5+sqrt(sdr(xy,z))**5)
                        di = sdi(xy,z)**3*deltaxyz2(yb+y-1)
     &                       /(sabsi(xy,z)**5+sqrt(sdi(xy,z))**5)
                     end if

                     mr(xy,z,1) = -sr(xy,z,1)*dr
                     mr(xy,z,2) = -sr(xy,z,2)*dr
                     mr(xy,z,3) = -sr(xy,z,3)*dr
                     mr(xy,z,4) = -sr(xy,z,4)*dr
                     mr(xy,z,5) = -sr(xy,z,5)*dr
                     mr(xy,z,6) = -sr(xy,z,6)*dr

                     mi(xy,z,1) = -si(xy,z,1)*di
                     mi(xy,z,2) = -si(xy,z,2)*di
                     mi(xy,z,3) = -si(xy,z,3)*di
                     mi(xy,z,4) = -si(xy,z,4)*di
                     mi(xy,z,5) = -si(xy,z,5)*di
                     mi(xy,z,6) = -si(xy,z,6)*di

                     br(xy,z,1) = sr(xy,z,1)
                     br(xy,z,2) = sr(xy,z,2)
                     br(xy,z,3) = sr(xy,z,3)
                     br(xy,z,4) = sr(xy,z,4)
                     br(xy,z,5) = sr(xy,z,5)
                     br(xy,z,6) = sr(xy,z,6)

                     bi(xy,z,1) = si(xy,z,1)
                     bi(xy,z,2) = si(xy,z,2)
                     bi(xy,z,3) = si(xy,z,3)
                     bi(xy,z,4) = si(xy,z,4)
                     bi(xy,z,5) = si(xy,z,5)
                     bi(xy,z,6) = si(xy,z,6)

                  end do
               end do
            end do

c     get the H2 filtered fields
            if (nproc.eq.1) then
               do i=1,3
                  call getxz(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,taur,taui)
               end do
               call getxz(gom2r(1,1,1),gom2i(1,1,1),yb,4,0,taur,taui)
               call getxz(gom2r(1,1,3),gom2i(1,1,3),yb,5,0,taur,taui)
c
c     Get the scalar and its y-derivative
c
               do ith=1,scalar
c                  call getxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
c     &                 yb,6+ith,0,taur,taui)
c                  call getxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
c     &                 yb,6+ith+scalar,0,taur,taui)
               end do
            else
#ifdef MPI
               do i=1,3
                  call getpxz(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,taur,taui,
     &                 realg1,realg2,my_node)
               end do
               call getpxz(gom2r(1,1,1),gom2i(1,1,1),yb,4,0,taur,taui,
     &              realg1,realg2,my_node)
               call getpxz(gom2r(1,1,3),gom2i(1,1,3),yb,5,0,taur,taui,
     &              realg1,realg2,my_node)
               do ith=1,scalar
c                  call getpxz(th2r(1,1,1,ith),th2i(1,1,1,ith),
c     &                 yb,6+ith,0,taur,taui,
c     &                 realg1,realg2,my_node)
c                  call getpxz(th2r(1,1,3,ith),th2i(1,1,3,ith),
c     &                 yb,6+ith+scalar,0,taur,taui,
c     &                 realg1,realg2,my_node)
               end do
#endif
            end if

c     calculate missing vorticity omega_y
            do z=1,nzc
               do y=yb,min(nyp,yb+mby-1)
                  do x=1,nx/2
                     xy=x+(y-yb)*(nxp/2+1)
                     gom2r(xy,z,2)=-beta(z)*gu2i(xy,z,1)+
     &                    alfa(x)*gu2i(xy,z,3)
                     gom2i(xy,z,2)=beta(z)*gu2r(xy,z,1)-
     &                    alfa(x)*gu2r(xy,z,3)
                  end do
                  gom2r(nx/2+1+(y-yb)*(nxp/2+1),z,2)=0.0
               end do
            end do
c
c     Calculate x and z derivatives for scalar field
c
            do ith=1,scalar
               do z=1,nzc
                  do y=yb,min(nyp,yb+mby-1)
                     do x=1,nx/2
                        xy=x+(y-yb)*(nxp/2+1)
c
c     Calculate x derivative
c
c                        th2r(xy,z,2,ith)=-alfa(x)*
c     &                       th2i(xy,z,1,ith)
c                        th2i(xy,z,2,ith)= alfa(x)*
c     &                       th2r(xy,z,1,ith)
c
c     Calculate z derivative
c
c                        th2r(xy,z,4,ith)=-beta(z)*
c     &                       th2i(xy,z,1,ith)
c                        th2i(xy,z,4,ith)= beta(z)*
c     &                       th2r(xy,z,1,ith)
                     end do
c                     th2r(nx/2+1+(y-yb)*(nxp/2+1),z,2,ith)=0.0
c                     th2r(nx/2+1+(y-yb)*(nxp/2+1),z,4,ith)=0.0
                  end do
               end do
            end do

            if (nfzsym.eq.0) then
c     zero oddball (z)
               do xy=1,(nxp/2+1)*npl
                  gom2r(xy,nz/2+1,2)=0.0
                  gom2i(xy,nz/2+1,2)=0.0
               end do
               do ith=1,scalar
                  do xy=1,(nxp/2+1)*npl
c                     th2r(xy,nz/2+1,2,ith)=0.0
c                     th2i(xy,nz/2+1,2,ith)=0.0
c                     th2r(xy,nz/2+1,4,ith)=0.0
c                     th2i(xy,nz/2+1,4,ith)=0.0
                  end do
               end do
            end if

c     subtract lower wall shifting velocity
c     not needed since HPF
c            do y=1,npl
c               xy=1+(y-1)*(nxp/2+1)
c               gu2r(xy,1,1)=gu2r(xy,1,1)-u0low
c               gu2r(xy,1,3)=gu2r(xy,1,3)-w0low
c            end do

c     calculate all velocity derivatives
c     horizontal derivatives by multiplikation by alfa or beta
            do i=1,3
               do z=1,nzc
                  do y=1,npl
                     do x=1,nx/2
                        xy=x+(y-1)*(nxp/2+1)
                        gdu2r(xy,z,i,1)=-alfa(x)*gu2i(xy,z,i)
                        gdu2i(xy,z,i,1)= alfa(x)*gu2r(xy,z,i)
                        gdu2r(xy,z,i,3)=-beta(z)*gu2i(xy,z,i)
                        gdu2i(xy,z,i,3)= beta(z)*gu2r(xy,z,i)
                     end do
                  end do
               end do
            end do

c     the wall normal derivatives can be found from :
c     omx  = dwdy - dvdz
c     omz  = dvdx - dudy
c     i.e. :
c     dudy = dvdx - omz
c     dvdy =-dudx - dwdz   (from continuity)
c     dwdy = omx  + dvdz
            do z=1,nzc
               do y=1,npl
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     gdu2r(xy,z,1,2)=-alfa(x)*gu2i(xy,z,2) -
     &                    gom2r(xy,z,3)
                     gdu2i(xy,z,1,2)= alfa(x)*gu2r(xy,z,2) -
     &                    gom2i(xy,z,3)
                     gdu2r(xy,z,2,2)= alfa(x)*gu2i(xy,z,1) +
     &                    beta(z)*gu2i(xy,z,3)
                     gdu2i(xy,z,2,2)=-alfa(x)*gu2r(xy,z,1) -
     &                    beta(z)*gu2r(xy,z,3)
                     gdu2r(xy,z,3,2)= gom2r(xy,z,1)        -
     &                    beta(z)*gu2i(xy,z,2)
                     gdu2i(xy,z,3,2)= gom2i(xy,z,1)        +
     &                    beta(z)*gu2r(xy,z,2)
                  end do
               end do
            end do
c
c     set oddball modes to zero
c
            do j=1,3
               do i=1,3
                  do z=1,nzc
                     do y=1,npl
                        gdu2r(nx/2+1+(y-1)*(nxp/2+1),z,i,j) = 0.0
                     end do
                  end do
               end do
            end do
c
c     transform to physical space
c
            do i=1,3
c
c     no shifting for LES quantities
c
               if (nfzsym.eq.0) then
                  call vcfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &                 nxy,(nxp/2+1)*mby,1,prezn)
                  call vcfftb(gom2r(1,1,i),gom2i(1,1,i),wr,wi,nz,
     &                 nxy,(nxp/2+1)*mby,1,prezn)
                  do j=1,3
                     call vcfftb(gdu2r(1,1,i,j),gdu2i(1,1,i,j),wr,wi,nz,
     &                    nxy,(nxp/2+1)*mby,1,prezn)
                  end do
               else
                  if (i.le.2) then
                     call vcffts(gu2r(1,1,i),gu2i(1,1,i),wr,
     &                    nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                     call vcftab(gom2r(1,2,i),gom2i(1,2,i),
     &                    wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                     call vcffts(gdu2r(1,1,i,1),gdu2i(1,1,i,1),wr,
     &                    nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                     call vcffts(gdu2r(1,1,i,2),gdu2i(1,1,i,2),wr,
     &                    nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                     call vcftab(gdu2r(1,2,i,3),gdu2i(1,2,i,3),
     &                    wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                  else
                     call vcftab(gu2r(1,2,i),gu2i(1,2,i),
     &                    wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                     call vcffts(gom2r(1,1,i),gom2i(1,1,i),wr,
     &                    nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                     call vcftab(gdu2r(1,2,i,1),gdu2i(1,2,i,1),
     &                    wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                     call vcftab(gdu2r(1,2,i,2),gdu2i(1,2,i,2),
     &                    wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
                     call vcffts(gdu2r(1,1,i,3),gdu2i(1,1,i,3),wr,
     &                    nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
                  end if
               end if
               call vrfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &              nx,nzpc*mby,1,nxp/2+1,prexn)
               call vrfftb(gom2r(1,1,i),gom2i(1,1,i),wr,wi,
     &              nx,nzpc*mby,1,nxp/2+1,prexn)
               do j=1,3
                  call vrfftb(gdu2r(1,1,i,j),gdu2i(1,1,i,j),wr,wi,
     &                 nx,nzpc*mby,1,nxp/2+1,prexn)
               end do
            end do
            do ith=1,scalar
               do i=1,4
c                  call vcfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
c     &                 ,wr,wi,nz,nxy,(nxp/2+1)*mby,1,prezn)
c                  call vrfftb(th2r(1,1,i,ith),th2i(1,1,i,ith)
c     &                 ,wr,wi,nx,nzpc*mby,1,nxp/2+1,prexn)
               end do
            end do

            do y=1,npl
               do z=1,nzc+nfzsym
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     br(xy,z,1) =     gdu2r(xy,z,1,1)
                     br(xy,z,2) = .5*(gdu2r(xy,z,2,1)+gdu2r(xy,z,1,2))
                     br(xy,z,3) = .5*(gdu2r(xy,z,3,1)+gdu2r(xy,z,1,3))
                     br(xy,z,4) =     gdu2r(xy,z,2,2)
                     br(xy,z,5) = .5*(gdu2r(xy,z,2,3)+gdu2r(xy,z,3,2))
                     br(xy,z,6) =     gdu2r(xy,z,3,3)

                     bi(xy,z,1) =     gdu2i(xy,z,1,1)
                     bi(xy,z,2) = .5*(gdu2i(xy,z,2,1)+gdu2i(xy,z,1,2))
                     bi(xy,z,3) = .5*(gdu2i(xy,z,3,1)+gdu2i(xy,z,1,3))
                     bi(xy,z,4) =     gdu2i(xy,z,2,2)
                     bi(xy,z,5) = .5*(gdu2i(xy,z,2,3)+gdu2i(xy,z,3,2))
                     bi(xy,z,6) =     gdu2i(xy,z,3,3)
                  end do
               end do
            end do
c
c     filter M, L, B, u
c
            call filter_cutoff(mr ,mi ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(lr ,li ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(br ,bi ,6,npl,nxy,prexn,prezn,wr,wi)
            call filter_cutoff(u2r,u2i,3,npl,nxy,prexn,prezn,wr,wi)
            if (iwale.eq.1) then
               call filter_cutoff(gdu2r,gdu2i,9,npl,nxy,prexn,prezn,
     &              wr,wi)
            end if
c
c     computation of L_ij and M_ij
c
            do y=1,npl
               gd2 = 4.**(2./3.)*deltaxyz2(yb+y-1)
               do z=1,nzc+nfzsym
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)

                     lr(xy,z,1) = lr(xy,z,1) + u2r(xy,z,1)*u2r(xy,z,1)
                     lr(xy,z,2) = lr(xy,z,2) + u2r(xy,z,2)*u2r(xy,z,1)
                     lr(xy,z,3) = lr(xy,z,3) + u2r(xy,z,3)*u2r(xy,z,1)
                     lr(xy,z,4) = lr(xy,z,4) + u2r(xy,z,2)*u2r(xy,z,2)
                     lr(xy,z,5) = lr(xy,z,5) + u2r(xy,z,3)*u2r(xy,z,2)
                     lr(xy,z,6) = lr(xy,z,6) + u2r(xy,z,3)*u2r(xy,z,3)

                     li(xy,z,1) = li(xy,z,1) + u2i(xy,z,1)*u2i(xy,z,1)
                     li(xy,z,2) = li(xy,z,2) + u2i(xy,z,2)*u2i(xy,z,1)
                     li(xy,z,3) = li(xy,z,3) + u2i(xy,z,3)*u2i(xy,z,1)
                     li(xy,z,4) = li(xy,z,4) + u2i(xy,z,2)*u2i(xy,z,2)
                     li(xy,z,5) = li(xy,z,5) + u2i(xy,z,3)*u2i(xy,z,2)
                     li(xy,z,6) = li(xy,z,6) + u2i(xy,z,3)*u2i(xy,z,3)

                     tsabsr = sqrt(
     &                    2.*(br(xy,z,1)**2+br(xy,z,4)**2+
     &                    br(xy,z,6)**2)+
     &                    4.*(br(xy,z,2)**2+br(xy,z,3)**2+
     &                    br(xy,z,5)**2))
                     tsabsi = sqrt(
     &                    2.*(bi(xy,z,1)**2+bi(xy,z,4)**2+
     &                    bi(xy,z,6)**2)+
     &                    4.*(bi(xy,z,2)**2+bi(xy,z,3)**2+
     &                    bi(xy,z,5)**2))
c                     tsabsr = sqrt(2*br(xy,z,4)**2)
c                     tsabsi = sqrt(2*bi(xy,z,4)**2)

                     if (iwale.eq.1) then
                        tsdr=sqrt((
     &                       gdu2r(xy,z,1,1)*gdu2r(xy,z,1,2) +
     &                       gdu2r(xy,z,1,2)*gdu2r(xy,z,2,2) +
     &                       gdu2r(xy,z,1,3)*gdu2r(xy,z,3,2) +
     &                       gdu2r(xy,z,2,1)*gdu2r(xy,z,1,1) +
     &                       gdu2r(xy,z,2,2)*gdu2r(xy,z,2,1) +
     &                       gdu2r(xy,z,2,3)*gdu2r(xy,z,3,1))**2/2. +(
     &                       gdu2r(xy,z,1,1)*gdu2r(xy,z,1,3) +
     &                       gdu2r(xy,z,1,2)*gdu2r(xy,z,2,3) +
     &                       gdu2r(xy,z,1,3)*gdu2r(xy,z,3,3) +
     &                       gdu2r(xy,z,3,1)*gdu2r(xy,z,1,1) +
     &                       gdu2r(xy,z,3,2)*gdu2r(xy,z,2,1) +
     &                       gdu2r(xy,z,3,3)*gdu2r(xy,z,3,1))**2/2. +(
     &                       gdu2r(xy,z,2,1)*gdu2r(xy,z,1,3) +
     &                       gdu2r(xy,z,2,2)*gdu2r(xy,z,2,3) +
     &                       gdu2r(xy,z,2,3)*gdu2r(xy,z,3,3) +
     &                       gdu2r(xy,z,3,1)*gdu2r(xy,z,1,2) +
     &                       gdu2r(xy,z,3,2)*gdu2r(xy,z,2,2) +
     &                       gdu2r(xy,z,3,3)*gdu2r(xy,z,3,2))**2/2.)

                        tsdi=sqrt((
     &                       gdu2i(xy,z,1,1)*gdu2i(xy,z,1,2) +
     &                       gdu2i(xy,z,1,2)*gdu2i(xy,z,2,2) +
     &                       gdu2i(xy,z,1,3)*gdu2i(xy,z,3,2) +
     &                       gdu2i(xy,z,2,1)*gdu2i(xy,z,1,1) +
     &                       gdu2i(xy,z,2,2)*gdu2i(xy,z,2,1) +
     &                       gdu2i(xy,z,2,3)*gdu2i(xy,z,3,1))**2/2. +(
     &                       gdu2i(xy,z,1,1)*gdu2i(xy,z,1,3) +
     &                       gdu2i(xy,z,1,2)*gdu2i(xy,z,2,3) +
     &                       gdu2i(xy,z,1,3)*gdu2i(xy,z,3,3) +
     &                       gdu2i(xy,z,3,1)*gdu2i(xy,z,1,1) +
     &                       gdu2i(xy,z,3,2)*gdu2i(xy,z,2,1) +
     &                       gdu2i(xy,z,3,3)*gdu2i(xy,z,3,1))**2/2. +(
     &                       gdu2i(xy,z,2,1)*gdu2i(xy,z,1,3) +
     &                       gdu2i(xy,z,2,2)*gdu2i(xy,z,2,3) +
     &                       gdu2i(xy,z,2,3)*gdu2i(xy,z,3,3) +
     &                       gdu2i(xy,z,3,1)*gdu2i(xy,z,1,2) +
     &                       gdu2i(xy,z,3,2)*gdu2i(xy,z,2,2) +
     &                       gdu2i(xy,z,3,3)*gdu2i(xy,z,3,2))**2/2.)
                     end if

                     if (iwale.eq.0) then
                        dr = gd2*tsabsr
                        di = gd2*tsabsi
                     else
                        dr = gd2*tsdr**3/(tsabsr**5+sqrt(tsdr)**5)
                        di = gd2*tsdi**3/(tsabsi**5+sqrt(tsdi)**5)
                     end if

                     mr(xy,z,1) = mr(xy,z,1) + br(xy,z,1)*dr
                     mr(xy,z,2) = mr(xy,z,2) + br(xy,z,2)*dr
                     mr(xy,z,3) = mr(xy,z,3) + br(xy,z,3)*dr
                     mr(xy,z,4) = mr(xy,z,4) + br(xy,z,4)*dr
                     mr(xy,z,5) = mr(xy,z,5) + br(xy,z,5)*dr
                     mr(xy,z,6) = mr(xy,z,6) + br(xy,z,6)*dr

                     mi(xy,z,1) = mi(xy,z,1) + bi(xy,z,1)*di
                     mi(xy,z,2) = mi(xy,z,2) + bi(xy,z,2)*di
                     mi(xy,z,3) = mi(xy,z,3) + bi(xy,z,3)*di
                     mi(xy,z,4) = mi(xy,z,4) + bi(xy,z,4)*di
                     mi(xy,z,5) = mi(xy,z,5) + bi(xy,z,5)*di
                     mi(xy,z,6) = mi(xy,z,6) + bi(xy,z,6)*di

                  end do
               end do
            end do
         end if
      end if

      if (cs.eq.0) then
c
c     computation of C = C_S^2 for the dynamic procedure
c
         do y=1,npl
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               cnumr = 0.
               cnumi = 0.
               cdenr = 0.
               cdeni = 0.
c
c     spanwise average
c
               do z=1,nzc+nfzsym
                  cnumr = cnumr + mr(xy,z,1)*lr(xy,z,1)+
     &                 mr(xy,z,4)*lr(xy,z,4)+mr(xy,z,6)*lr(xy,z,6)+
     &                 2.*(mr(xy,z,2)*lr(xy,z,2)+mr(xy,z,3)*lr(xy,z,3)+
     &                 mr(xy,z,5)*lr(xy,z,5))
                  cnumi = cnumi + mi(xy,z,1)*li(xy,z,1)+
     &                 mi(xy,z,4)*li(xy,z,4)+mi(xy,z,6)*li(xy,z,6)+
     &                 2.*(mi(xy,z,2)*li(xy,z,2)+mi(xy,z,3)*li(xy,z,3)+
     &                 mi(xy,z,5)*li(xy,z,5))
                  cdenr = cdenr + mr(xy,z,1)*mr(xy,z,1)+
     &                 mr(xy,z,4)*mr(xy,z,4)+mr(xy,z,6)*mr(xy,z,6)+
     &                 2.*(mr(xy,z,2)*mr(xy,z,2)+mr(xy,z,3)*mr(xy,z,3)+
     &                 mr(xy,z,5)*mr(xy,z,5))
                  cdeni = cdeni + mi(xy,z,1)*mi(xy,z,1)+
     &                 mi(xy,z,4)*mi(xy,z,4)+mi(xy,z,6)*mi(xy,z,6)+
     &                 2.*(mi(xy,z,2)*mi(xy,z,2)+mi(xy,z,3)*mi(xy,z,3)+
     &                 mi(xy,z,5)*mi(xy,z,5))
               end do

               if (ineg.eq.1) then
c
c     clipping to positive values C>0 if ineg=1
c
                  if (yb.eq.nyp.or.yb.eq.1) then
c
c     set to zero at the wall (nyp) and at the upper boundary (1)
c     here we could include a smooth damping at the upper boundary
c
                     do z=1,nzc+nfzsym
                        cr(xy,z) = 0.
                        ci(xy,z) = 0.
                        cthr(xy,z) = 0.
                        cthi(xy,z) = 0.
                     end do
                  else
                     do z=1,nzc+nfzsym
                        cr(xy,z) = max(.5*cnumr/cdenr,0.)
                        ci(xy,z) = max(.5*cnumi/cdeni,0.)
                        cthr(xy,z) = max(.5*cnumr/cdenr,0.)
                        cthi(xy,z) = max(.5*cnumi/cdeni,0.)
                        if (iwale.eq.1.and.scalar.gt.0) then
c
c     clipping to values 0.3<C<0.5 if ineg=1
c
c                           cthr(xy,z) = min(.3,cr(xy,z))
c                           cthi(xy,z) = min(.3,ci(xy,z))
                        end if
                     end do
                  end if

                  ydamp = 0.
                  if (gridy(yb).ge.ydamp.and.ydamp.gt.0.) then
                     gd2 = 1.-step((gridy(yb)-ydamp)/(gridy(1)-ydamp))
                     do z=1,nzc+nfzsym
                        cr(xy,z) = cr(xy,z)*gd2
                        ci(xy,z) = ci(xy,z)*gd2
                     end do
                  end if

c
c     Van Driest damping
c
c     here A+=5, and plus unit is for case Re_tau =180 in channel
c     note that here the flag is also ydamp, the same as step damping in the free stream
c                  if (ydamp.gt.0.) then
c                     gd2 = 1.-exp(-(gridy(yb)**3/0.275**3))
c                     do z=1,nzc+nfzsym
c                        cr(xy,z) = cr(xy,z)*gd2
c                        ci(xy,z) = ci(xy,z)*gd2
c                        cthr(xy,z) = cthr(xy,z)*gd2
c                        cthi(xy,z) = cthi(xy,z)*gd2
c                     end do
c                  end if

               else if (ineg.eq.2) then
c
c     positive total viscosity if ineg=2
c
                  if (yb.eq.nyp.or.yb.eq.1) then
c
c     set to zero at the wall (nyp) and at the upper boundary (1)
c
                     do z=1,nzc+nfzsym
                        cr(xy,z) = 0.
                        ci(xy,z) = 0.
                     end do
                  else
                     do z=1,nzc+nfzsym
                        cr(xy,z) = .5*cnumr/cdenr
                        ci(xy,z) = .5*cnumi/cdeni
                     end do
                  end if
               else if (ineg.eq.0) then
c
c     no clipping
c
                  do z=1,nzc+nfzsym
                     cr(xy,z) = .5*cnumr/cdenr
                     ci(xy,z) = .5*cnumi/cdeni
                     cthr(xy,z) = .5*cnumr/cdenr
                     cthi(xy,z) = .5*cnumi/cdeni
                  end do
               else
                  call stopnow(54754)
               end if
            end do
         end do

      else
c
c     Fixed model coefficient for non-dynamic procedure
c
         do y=1,npl
            do x=1,nx/2
               do z=1,nzc+nfzsym
                  xy=x+(y-1)*(nxp/2+1)
                  cr(xy,z) = cs**2
                  ci(xy,z) = cs**2
                  cthr(xy,z) = cs**2
                  cthi(xy,z) = cs**2
               end do
            end do
         end do
      end if
c
c     Here we can modify the model coefficient for the
c     eddy-viscosity (HPF/classic), e.g. wall-damping or
c     damping close to the upper wall.
c





c
c     Computation of SGS stresses tau_ij with FFT normalisation
c     (in physical space)
c
      gd2 = 1./real(nx*nz)
      do y=1,npl
         do z=1,nz
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
c
c     Eddy viscsity nu_t
c
               if (iwale.eq.0) then
                  dr = cr(xy,z)*sabsr(xy,z)*deltaxyz2(yb+y-1)
                  di = ci(xy,z)*sabsi(xy,z)*deltaxyz2(yb+y-1)
               else
                  dr = cr(xy,z)*sdr(xy,z)**3*deltaxyz2(yb+y-1)
     &                 /(sabsr(xy,z)**5+sqrt(sdr(xy,z))**5)
                  di = ci(xy,z)*sdi(xy,z)**3*deltaxyz2(yb+y-1)
     &                 /(sabsi(xy,z)**5+sqrt(sdi(xy,z))**5)
               end if
c
c     Apply clippings to eddy viscosity
c
               if (ineg.eq.0.or.ineg.eq.1) then
c
c     Nothing to be done (clipping applied above)
c
               else if (ineg.eq.2) then
c
c     Clip to nu_t+nu>0 and recompute Smagorinsky coefficient
c
                  if (dr.lt.-1./re) then
                     dr = -1./re
                     cr(xy,z) = dr/sabsr(xy,z)/deltaxyz2(yb+y-1)
                  end if
                  if (di.lt.-1./re) then
                     di = -1./re
                     ci(xy,z) = di/sabsi(xy,z)/deltaxyz2(yb+y-1)
                  end if
               else
                  call stopnow(74334)
               end if
c
c     tau_ij (with normalisation)
c
               d1r = -2.*dr*gd2
               d1i = -2.*di*gd2
               do ll=1,6
                  mr(xy,z,ll) = d1r*sr(xy,z,ll)
                  mi(xy,z,ll) = d1i*si(xy,z,ll)
               end do
c
c     Recompute nu_t for the scalar of the WALE model
c
               if (iwale.eq.1) then
                  dr = cthr(xy,z)*sdr(xy,z)**3*deltaxyz2(yb+y-1)
     &                 /(sabsr(xy,z)**5+sqrt(sdr(xy,z))**5)
                  di = cthi(xy,z)*sdi(xy,z)**3*deltaxyz2(yb+y-1)
     &                 /(sabsi(xy,z)**5+sqrt(sdi(xy,z))**5)
               end if
c
c     Constant turbulent Prandtl number for the eddy-diffusivity model
c
               ktr = dr/prt
               kti = di/prt
c
c     SGS force for scalar (with normalisation)
c
               ktr = -ktr*gd2
               kti = -kti*gd2
               do ith=1,scalar
                  do ll=2,4
                     th2r(xy,z,ll,ith) = ktr*th2r(xy,z,ll,ith)
                     th2i(xy,z,ll,ith) = kti*th2i(xy,z,ll,ith)
                  end do
               end do

            end do
         end do
      end do
c
c     Transform to Fourier space, already multiplied by gd2=1/(nx*nz)
c
      do i=1,6
         call vrfftf(mr(1,1,i),mi(1,1,i),wr,wi,
     &        nx,nzpc*mby,1,nxp/2+1,prexn)
         call vcfftf(mr(1,1,i),mi(1,1,i),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)
      end do
      do ith=1,scalar
         do i=2,4
            call vrfftf(th2r(1,1,i,ith),th2i(1,1,i,ith),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            call vcfftf(th2r(1,1,i,ith),th2i(1,1,i,ith),wr,wi,nz,
     &           nxy,(nxp/2+1)*mby,1,prezn)
         end do
      end do
c
c     Zero odd-ball mode of tau
c
      do i=1,6
         do y=1,npl

            z=nz/2+1
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               mr(xy,z,i) = 0.
               mi(xy,z,i) = 0.
            end do

            x=nx/2+1
            xy=x+(y-1)*(nxp/2+1)
            do z=1,nz
               mr(xy,z,i) = 0.
               mi(xy,z,i) = 0.
            end do
         end do
      end do
      do ith=1,scalar
         do i=2,4
            do y=1,npl

               z=nz/2+1
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  th2r(xy,z,i,ith) = 0.
                  th2i(xy,z,i,ith) = 0.
               end do

               x=nx/2+1
               xy=x+(y-1)*(nxp/2+1)
               do z=1,nz
                  th2r(xy,z,i,ith) = 0.
                  th2i(xy,z,i,ith) = 0.
               end do
            end do
         end do
      end do
c
c     Blow up array (could be avoided by changing putxz)
c
      lr=0.
      li=0.
      do ll=1,6
         do y=1,npl
            do z=1,nz/2
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  lr(xy,z,ll) = mr(xy,z,ll)
                  li(xy,z,ll) = mi(xy,z,ll)
               end do
            end do
            do z=nz/2+1,nz
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  zz = z+nzp-nz
                  lr(xy,zz,ll) = mr(xy,z,ll)
                  li(xy,zz,ll) = mi(xy,z,ll)
               end do
            end do
         end do
      end do
c
c     Blow up for the scalar flux
c
      lsr=0.
      lsi=0.
      do ith=1,scalar
         do ll=1,3
            do y=1,npl
               do z=1,nz/2
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     lsr(xy,z,ll,ith) = th2r(xy,z,ll+1,ith)
                     lsi(xy,z,ll,ith) = th2i(xy,z,ll+1,ith)
                  end do
               end do
               do z=nz/2+1,nz
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     zz = z+nzp-nz
                     lsr(xy,zz,ll,ith) = th2r(xy,z,ll+1,ith)
                     lsi(xy,zz,ll,ith) = th2i(xy,z,ll+1,ith)
                  end do
               end do
            end do
         end do
      end do
c
c     Put tau_ij onto global field
c
      if (nproc.eq.1) then
         do i=1,6
            call putxz(lr(1,1,i),li(1,1,i),yb,i,taur,taui)
         end do
         do ith=1,scalar
            do i=1,3
               call putxz(lsr(1,1,i,ith),lsi(1,1,i,ith),
     &              yb,6+i+(ith-1)*3,taur,taui)
            end do
         end do
      else
#ifdef MPI
         do i=1,6
            call putpxz(lr(1,1,i),li(1,1,i),yb,i,taur,taui,
     &           realg1,realg2,my_node)
         end do
         do ith=1,scalar
            do i=1,3
               call putpxz(lsr(1,1,i,ith),lsi(1,1,i,ith),
     &              yb,6+i+(ith-1)*3,taur,taui,
     &              realg1,realg2,my_node)
            end do
         end do
#endif
      end if
c
c-----------------------------------------------------------------------
c
c     Do the statistics
c
      if (.not.do_stat) return

      if (nxys.ge.52) then
c
c     if xs or zs is not zero, the statistics should be moved here
c     using xzsh. This is however not implemented.
c
         if (xs.ne.0.or.zs.ne.0.or.gall) then
            write(*,*) 'no moving walls allowed.'
            write(*,*) 'otherwise the LES statistics are incorrect.'
            call stopnow(549684)
         end if
c
c     Transform tau_ij back to physical space
c     not needed anymore (see below)
c
         do i=1,6
            call vcfftb(mr(1,1,i),mi(1,1,i),wr,wi,nz,
     &           nxy,(nxp/2+1)*mby,1,prezn)
            call vrfftb(mr(1,1,i),mi(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
         end do

         c1=dtn/real(nzc)
         do y=1,npl
            do z=1,nzc+nfzsym
               c=c1
               if (nfzsym.eq.1.and.(z.eq.1.or.z.eq.nzc+1)) c=c1/2.
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
c     Smagorinsky coefficient (position 52)
                  xys(2*x-1,ybp+y-1,52)=xys(2*x-1,ybp+y-1,52)+
     &                 cr(xy,z)*c
                  xys(2*x  ,ybp+y-1,52)=xys(2*x  ,ybp+y-1,52)+
     &                 ci(xy,z)*c

c     <nu_t> (position 60)
c     observe the clipping introduced above
                  if (iwale.eq.0) then
                     xys(2*x-1,ybp+y-1,60)=xys(2*x-1,ybp+y-1,60)+
     &                    cr(xy,z)*sabsr(xy,z)*deltaxyz2(yb+y-1)*c
                     xys(2*x  ,ybp+y-1,60)=xys(2*x  ,ybp+y-1,60)+
     &                    ci(xy,z)*sabsi(xy,z)*deltaxyz2(yb+y-1)*c
                  else
                     xys(2*x-1,ybp+y-1,60)=xys(2*x-1,ybp+y-1,60)+
     &                    cr(xy,z)*sdr(xy,z)**3*deltaxyz2(yb+y-1)*c
     &                    /(sabsr(xy,z)**5+sqrt(sdr(xy,z))**5)
                     xys(2*x  ,ybp+y-1,60)=xys(2*x  ,ybp+y-1,60)+
     &                    ci(xy,z)*sdi(xy,z)**3*deltaxyz2(yb+y-1)*c
     &                    /(sabsi(xy,z)**5+sqrt(sdi(xy,z))**5)
                  end if

c     the computation of <tau_ij S_ij> and <tau_ij> has been moved to
c     boxxys in order to get the correct statistics for models not based
c     on S_ij, since in this routine S_ij is (generally) not available.
c     However, slightly more CPU time is needed as tau has
c     to be transposed once more by getpxz.
c     in case, also uncomment the above FFT.
c
c     <tau_ij S_ij> (position 53)
c                  dr = mr(xy,z,1)*sr(xy,z,1)+mr(xy,z,4)*sr(xy,z,4)+
c     &                 mr(xy,z,6)*sr(xy,z,6)+2.*mr(xy,z,2)*sr(xy,z,2)+
c     &                 2.*mr(xy,z,3)*sr(xy,z,3)+2.*mr(xy,z,5)*sr(xy,z,5)
c                  di = mi(xy,z,1)*si(xy,z,1)+mi(xy,z,4)*si(xy,z,4)+
c     &                 mi(xy,z,6)*si(xy,z,6)+2.*mi(xy,z,2)*si(xy,z,2)+
c     &                 2.*mi(xy,z,3)*si(xy,z,3)+2.*mi(xy,z,5)*si(xy,z,5)
c                  xys(2*x-1,ybp+y-1,53) = xys(2*x-1,ybp+y-1,53) + dr*c
c                  xys(2*x  ,ybp+y-1,53) = xys(2*x  ,ybp+y-1,53) + di*c
c     <tau_11>
c                  xys(2*x-1,ybp+y-1,54) = xys(2*x-1,ybp+y-1,54) +
c     &                 mr(xy,z,1)*c
c                  xys(2*x  ,ybp+y-1,54) = xys(2*x  ,ybp+y-1,54) +
c     &                 mi(xy,z,1)*c
c     <tau_12>
c                  xys(2*x-1,ybp+y-1,55) = xys(2*x-1,ybp+y-1,55) +
c     &                 mr(xy,z,2)*c
c                  xys(2*x  ,ybp+y-1,55) = xys(2*x  ,ybp+y-1,55) +
c     &                 mi(xy,z,2)*c
c     <tau_13>
c                  xys(2*x-1,ybp+y-1,56) = xys(2*x-1,ybp+y-1,56) +
c     &                 mr(xy,z,3)*c
c                  xys(2*x  ,ybp+y-1,56) = xys(2*x  ,ybp+y-1,56) +
c     &                 mi(xy,z,3)*c
c     <tau_22>
c                  xys(2*x-1,ybp+y-1,57) = xys(2*x-1,ybp+y-1,57) +
c     &                 mr(xy,z,4)*c
c                  xys(2*x  ,ybp+y-1,57) = xys(2*x  ,ybp+y-1,57) +
c     &                 mi(xy,z,4)*c
c     <tau_23>
c                  xys(2*x-1,ybp+y-1,58) = xys(2*x-1,ybp+y-1,58) +
c     &                 mr(xy,z,5)*c
c                  xys(2*x  ,ybp+y-1,58) = xys(2*x  ,ybp+y-1,58) +
c     &                 mi(xy,z,5)*c
c     <tau_33> (position 59)
c                  xys(2*x-1,ybp+y-1,59) = xys(2*x-1,ybp+y-1,59) +
c     &                 mr(xy,z,6)*c
c                  xys(2*x  ,ybp+y-1,59) = xys(2*x  ,ybp+y-1,59) +
c     &                 mi(xy,z,6)*c
               end do
            end do
         end do
      end if

      end subroutine les
