c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine nonlinbl(amp,campw,kx,kz,nwave,cflp,vext,cext,iext,
     &     ur,ui,tc,xsc,zsc,xs,yb,rot,spat,fltype,
     &     fstart,fend,nst,
     &     bu1,bu2,
     &     loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,
     &     fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,
     &     tripf,txsc,tysc,ttzc,tx0,
     &     osmod,osnumb,
     &     osur,osui,osvr,osvi,oswr,oswi,
     &     xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &     evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &     ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &     fbla,nbla,dybla,rlam,cdev,re,xblc,
     &     alfa,beta,eta,deta,xl,zl,prex,prez,pres,prea,
     &     it,icfl,iamp,u2r,u2i,om2r,om2i,th2r,th2i,wr,wi,
     &     iles,gur,gui,chi_scaled,prt,
     &     my_node,realg1,realg2,
     &     ampst,streak,betast,omegast,
     &     ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &     tsmoo,tsmst,tsmend,iampst,phist,
     &     isfd,sfd_chi,sfd_delta,sfdur,sfdui,an,bn,wbci,
     &     pert,lin,bom1,bom2,x0,spanv,rv1r,rv1i,do_press_nl,
     &     wtcom,wtser,ctcom,ctser,gsr,gsi,fmhdr,fmhdi,imhd,du2r,du2i,
     &     b0,mhd_n,bf3,bf3u2r,bf3u2i,ybp,bf3tempr,bf3tempi,PINT,
     &     pr,gr,dstar)
c
c     Calculates the nonlinear terms, accumulates amp and cfl for a xz-box
c     i.e. nby times xz-planes
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif
      logical spat
      integer nwave,yb,it,icfl,iamp,kx(nwave),kz(nwave),iext
      complex campw(nyp,4,nwave)
      real cflp,cflp_bf,amp(nyp,20)
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real rot,alfa(nx/2*mby),beta(nz)
      integer fltype,loctyp,nst
      real fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5
      logical tripf
      real txsc,tysc,ttzc(nzp+2,4),tx0
      real fstart,fend
      real evr(nyp),evi(nyp),eur(nyp),eui(nyp),ewr(nyp),ewi(nyp)
      real afw,alfaw,betaw,ampob,amp2d
      real ev2r(nyp),ev2i(nyp),eu2r(nyp),eu2i(nyp),afw2d,alf2d
      real eta(nyp),deta(nyp),xl,zl,tc,xsc,zsc,xs,cdev,re,xblc
      real rlam
      real fbla(mbla,7+3*scalar)
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real du2r((nxp/2+1)*mby,nzd,3,3)
      real du2i((nxp/2+1)*mby,nzd,3,3)
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)
      real x0,spanv
      real d1f,d2f,d3f,d2g,dybla
      integer nbla
c
c     Temporal forcing
c
      real c1,xx,etabl,ybl,h3u
c
c     Linearized-perturbation
c
      logical pert,lin
      real bom1(nxp/2+1,nyp,3),bom2(nxp/2+1,nyp,3)
c
c     Scalar
c
      real th2r((nxp/2+1)*mby,nzd,4*scalar)
      real th2i((nxp/2+1)*mby,nzd,4*scalar)
      integer ith
      real pr(scalar),gr(scalar),dstar
c
c     Pressure
c
      real rv1r(nx/2,nzc,2*pressure)
      real rv1i(nx/2,nzc,2*pressure)
      logical do_press_nl
c
c     Particles
c
      integer p
      real uf_1,uf_2,uf_3,uf_4,uf_5,uf_6,uf_7
      real vf_1,vf_2,vf_3,vf_4,vf_5,vf_6,vf_7
      real wf_1,wf_2,wf_3,wf_4,wf_5,wf_6,wf_7
      real x_left,z_left
      real appo,appo_x
      integer appo1,appo2
      integer i_index,j_index,k_index,ii_index,kk_index
      real PINT(3,npart,3)
      integer nn_index
c
c     Constant
c
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     OS modes
c
      logical osmod
      integer osnumb
      real osur(osnf,nyp), osvr(osnf,nyp), oswr(osnf,nyp)
      real osui(osnf,nyp), osvi(osnf,nyp), oswi(osnf,nyp)
      real cphas1(osnf,nxp/2,nzd),cphas2(osnf,nxp/2,nzd)
      real sphas1(osnf,nxp/2,nzd),sphas2(osnf,nxp/2,nzd)

      real fring1(nxp/2),fring2(nxp/2)
      real xc1(nxp/2),xc2(nxp/2)

      logical sym
      integer i,xy,z,zp,nxy,y,y1,x,nxp2,npl,npp
      real h1u,h2u
c
c     Jet in crossflow
c
      integer wbci
c
c     Streak generation
c
      logical streak
      real iampst,ampst(2),tsmoo(4),tsmst(2),tsmend(2)
      real betast(2),omegast(2),phist
      real uust_r(nyp,180,2),vvst_r(nyp,180,2),wwst_r(nyp,180,2)
      real uust_i(nyp,180,2),vvst_i(nyp,180,2),wwst_i(nyp,180,2)
      integer ndxst
c
c     LES
c
      real gur(memnx,memny,memnz,5)
      real gui(memnx,memny,memnz,5)
      real gsr(memnx,memny,memnz,scalar)
      real gsi(memnx,memny,memnz,scalar)
      integer iles,ll
      real chi_scaled,prt
c
c     MHD
c
      real fmhdr(memnx,memny,memnz,2)
      real fmhdi(memnx,memny,memnz,2)
      integer imhd
      real cr1,cr2,cr3
      real ci1,ci2,ci3
      real b0(3),mhd_n,fact
c
c     SFD
c
      integer isfd
      real sfd_chi,sfd_delta
      real sfdur(memnx,memny,memnz,6)
      real sfdui(memnx,memny,memnz,6)
      real an,bn
c
c     3D baseflow
c
      logical bf3
      integer ybp
      real bf3u2r(nxp/2+1,nzd,nyp/nproc+1,6)
      real bf3u2i(nxp/2+1,nzd,nyp/nproc+1,6)
      real bf3tempr(nxp/2+1,nzd,3)
      real bf3tempi(nxp/2+1,nzd,3)
c
c     MPI
c
      integer my_node
      integer realg1,realg2
#ifdef MPI
      integer ierror
#endif
c

c     Functions
c
      real, external :: cubip

      real ct1,wt1,ct2,wt2
      real wtcom,wtser,ctcom,ctser

c      call mpi_barrier(mpi_comm_world,ierror)
c      call wall_time(t1)
      ct1=0
      wt1=0
      ct2=0
      wt2=0
      call ctim(ct2,wt2)
      ct1=ct2
      wt1=wt2



      if (mby.gt.1.and.pert) then
         write(*,*)'mby>1 and pert----> not implemented'
         call stopnow(1009)
      end if

      npp=min(nyp,yb+mby-1)
      npl=min(mby,nyp-yb+1)
      nxp2=nxp/2+1
      nxy=nxp2*npl
c
c     u2r,u2i,om2r,om2i are on the dealiasing grid
c     middle is padded, i.e. nonzero values:
c     x = 1..nx/2
c          nx/2
c     z = 1..nz/2, nzp-nz/2+2..nzp
c           nz/2      nz/2-1
c     (if ipad=0 z=1..nz/2, nz-nz/2+2)
c                  nz/2       nz/2-1

c
c     Get the velocities (padded)
c
      if (nproc.eq.1) then
         do i=1,3
            call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui)
         end do
c
c     Get the streamwise and spanwise vorticity
c
         call getxz(om2r(1,1,1),om2i(1,1,1),yb,4,1,ur,ui)
         call getxz(om2r(1,1,3),om2i(1,1,3),yb,5,1,ur,ui)

         do ith=1,scalar
c
c     Get the scalar and its derivative
c
            call getxz(th2r(1,1,1+4*(ith-1)),th2i(1,1,1+4*(ith-1)),yb,
     &           8+pressure+3*(ith-1),1,ur,ui)
            call getxz(th2r(1,1,3+4*(ith-1)),th2i(1,1,3+4*(ith-1)),yb,
     &           9+pressure+3*(ith-1),1,ur,ui)
         end do
      else
#ifdef MPI
         do i=1,3
            call getpxz(u2r(1,1,i),u2i(1,1,i),yb,i,1,ur,ui,
     &           realg1,realg2,my_node)
         end do
         call getpxz(om2r(1,1,1),om2i(1,1,1),yb,4,1,ur,ui,
     &        realg1,realg2,my_node)
         call getpxz(om2r(1,1,3),om2i(1,1,3),yb,5,1,ur,ui,
     &        realg1,realg2,my_node)
         do ith=1,scalar
            call getpxz(th2r(1,1,1+4*(ith-1)),th2i(1,1,1+4*(ith-1)),yb,
     &           8+pressure+3*(ith-1),1,ur,ui,
     &           realg1,realg2,my_node)
            call getpxz(th2r(1,1,3+4*(ith-1)),th2i(1,1,3+4*(ith-1)),yb,
     &           9+pressure+3*(ith-1),1,ur,ui,
     &           realg1,realg2,my_node)
         end do
#endif
      end if

c      call mpi_barrier(mpi_comm_world,ierror)
c      call wall_time(t2)
      call ctim(ct2,wt2)
      wtcom = wtcom + (wt2-wt1)
      ctcom = ctcom + (ct2-ct1)
      ct1=ct2
      wt1=wt2
c
c     Save the velocities for MHD
c     for later use when computing the Lorentz force
c
      if (imhd.eq.1) then
         do ll=1,3
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     du2r(xy,z,ll,2)=u2r(xy,z,ll)
                     du2i(xy,z,ll,2)=u2i(xy,z,ll)
                  end do
               end do
            end do
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     du2r(xy,zp,ll,2)=u2r(xy,zp,ll)
                     du2i(xy,zp,ll,2)=u2i(xy,zp,ll)
                  end do
               end do
            end do
         end do
      end if
c
c     Compute the wall-normal vorticity
c
      do z=1,nz/2
         do y=yb,npp
            y1=(y-yb)*nxp2
            do x=1,nx/2
               xy=x+y1
               om2r(xy,z,2)=-beta(z)*u2i(xy,z,1)+alfa(x)*u2i(xy,z,3)
               om2i(xy,z,2)= beta(z)*u2r(xy,z,1)-alfa(x)*u2r(xy,z,3)
            end do
         end do
      end do

      do ith=1,scalar
c
c     Compute the derivatives of theta th2(:,:,1) in 2,3,4
c
         do z=1,nz/2
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1
                  th2r(xy,z,4+4*(ith-1))=-beta(z)*th2i(xy,z,1+4*(ith-1))
                  th2i(xy,z,4+4*(ith-1))= beta(z)*th2r(xy,z,1+4*(ith-1))

                  th2r(xy,z,2+4*(ith-1))=-alfa(x)*th2i(xy,z,1+4*(ith-1))
                  th2i(xy,z,2+4*(ith-1))= alfa(x)*th2r(xy,z,1+4*(ith-1))

               end do
            end do
         end do
      end do

      if (nfzsym.eq.0) then
c
c     Do it also for the upper half in the z-direction
c
         do z=nz/2+1,nz
            zp=nzp-nz+z
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1
                  om2r(xy,zp,2) = -beta(z)*u2i(xy,zp,1) +
     &                 alfa(x)*u2i(xy,zp,3)
                  om2i(xy,zp,2) = beta(z)*u2r(xy,zp,1) -
     &                 alfa(x)*u2r(xy,zp,3)
               end do
            end do
         end do
         do ith=1,scalar
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     th2r(xy,zp,4+4*(ith-1)) = -beta(z)*
     &                    th2i(xy,zp,1+4*(ith-1))
                     th2i(xy,zp,4+4*(ith-1)) =  beta(z)*
     &                    th2r(xy,zp,1+4*(ith-1))

                     th2r(xy,zp,2+4*(ith-1)) = -alfa(x)*
     &                    th2i(xy,zp,1+4*(ith-1))
                     th2i(xy,zp,2+4*(ith-1)) =  alfa(x)*
     &                    th2r(xy,zp,1+4*(ith-1))
                  end do
               end do
            end do
         end do
      end if
c
c     And pad/oddball
c
      do z=(nz+1)/2+1,min(nzpc,nzp+1-(nz+1)/2)
         do xy=1,nxy
            om2r(xy,z,2)=0.0
            om2i(xy,z,2)=0.0
         end do
      end do
      do z=1,nzpc
         do y=yb,npp
            y1=(y-yb)*nxp2
            do x=nx/2+1,nxp2
               xy=x+y1
               om2r(xy,z,2)=0.0
               om2i(xy,z,2)=0.0
            end do
         end do
      end do

      do ith=1,scalar
         do z=(nz+1)/2+1,min(nzpc,nzp+1-(nz+1)/2)
            do xy=1,nxy
               th2r(xy,z,2+4*(ith-1))=0.0
               th2i(xy,z,2+4*(ith-1))=0.0

               th2r(xy,z,4+4*(ith-1))=0.0
               th2i(xy,z,4+4*(ith-1))=0.0
            end do
         end do
         do z=1,nzpc
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=nx/2+1,nxp2
                  xy=x+y1
                  th2r(xy,z,2+4*(ith-1))=0.0
                  th2i(xy,z,2+4*(ith-1))=0.0

                  th2r(xy,z,4+4*(ith-1))=0.0
                  th2i(xy,z,4+4*(ith-1))=0.0
               end do
            end do
         end do
      end do
c
c     Now we have all velocities and vorticities on the
c     fine grid padded/without oddball
c

c
c     Compute/accumulate the amplitudes
c
      if (mod(it-1-nst,iamp).eq.0.and.iamp.gt.0) then
         call boxamp(amp,campw,kx,kz,nwave,
     &        u2r,u2i,om2r,om2i,yb,alfa,beta)
      end if
c
c     Pressure: store v at the wall and in the free stream for
c               Neumann boundary condition at the last full time step
c
c     ATTENTION: MPI
c      if (pressure.eq.1.and.mod(it-1,nst).eq.0) then
      if (pressure.eq.1.and.do_press_nl) then

         if (yb.eq.1.or.yb.eq.(nby-1)*mby+1) then
c
c     i is 1 for yb=1 and i is 2 for yb=(nby-1)*mby+1 (usually nyp)
c
            i = min(2,yb)
            y1 = (nyp-yb)*nxp2*(1-1/yb)
            do z = 1,1/nz+nz/2*(2-nfzsym)
               zp = z+z/(nz/2+1)*(nzp-nz)
               do x=1,nx/2
                  rv1r(x,z,i) = u2r(x+y1,zp,2)
                  rv1i(x,z,i) = u2i(x+y1,zp,2)
c
c     For Dirichlet boundary conditions one would use
c     rv1r(x,z,i)=-alfa(x)*u2i(x+y1,zp,1)-beta(z)*u2i(x+y1,zp,3)
c     rv1i(x,z,i)= alfa(x)*u2r(x+y1,zp,1)+beta(z)*u2r(x+y1,zp,3)
c
               end do
            end do
c            write(*,*) yb,' is on node ',my_node
         end if

#ifdef MPI
         if (nproc.gt.1) then
c
c     Communicate planes if necessary
c
c     yb=1
c
            if (yb.le.nproc) then
               i = 0
c               write(*,*) 'communicate yb=1',my_node,' from ',i
               call mpi_bcast(rv1r(1,1,1),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
               call mpi_bcast(rv1i(1,1,1),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
            end if
c
c     yb=nyp
c
            if (yb.gt.(nyp/nproc)*nproc) then
               i = mod(nyp-1,nproc)
c               write(*,*) 'communicate yb=nyp',my_node,' from ',i
               call mpi_bcast(rv1r(1,1,2),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
               call mpi_bcast(rv1i(1,1,2),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)

            end if
         end if
#endif

      end if
c
c     Backward Fourier transform (Fourier synthesis)
c
      do i=1,3
         sym=i.le.2
         call fft2db(u2r(1,1,i),u2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
         sym=i.eq.3
         call fft2db(om2r(1,1,i),om2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
      end do

      do ith=1,scalar
         do i=1,4
            sym = i.le.3
            call fft2db(th2r(1,1,i+4*(ith-1)),th2i(1,1,i+4*(ith-1)),sym,
     &           npl,prex,prez,pres,prea,wr,wi)
         end do
      end do
c
c     We are now in the physical space (on dealiasing grid)
c     x=1..nxp
c     z=1..nzp
c     y = yb..yb+mby-1 (y=1 is the upper boundary)


c
c_______________________________________________________________________
c    PARTICLES
c
      if (npart.gt.0) then
         do p=1,npart
            appo=pint(2,p,1)
            j_index = int(acos(appo)/pi*real(nyp-1)+1.)
c            j_index = max(j_index,1)
c            j_index = min(j_index,nyp-1)
c            write(*,*) j_index,appo

c
c     Check whether we have a particle in the present plane
c
            if ((j_index.eq.yb).or.(j_index+1.eq.yb)) then
c
c     Compute the indices
c
               nn_index=2
               if(j_index+1.eq.yb) nn_index=3

               appo_x=mod(PINT(1,p,1),xl)

               i_index=modulo(int(appo_x/xl*real(nxp))+nxp/2,nxp)+1
               ii_index = mod(i_index,2)
               i_index = (i_index+1)/2

               x_left=
     &              real(2*i_index-1*ii_index-nxp/2-1)/real(nxp)*xl+xsc
               x_left=x_left-int((x_left+xl/2.)/xl)*xl
               x_left=x_left+xl*(1.-sign(1.,x_left))/2.

               appo=PINT(3,p,1)
               k_index = modulo(int(appo/zl*real(nzd))+nzd/2,nzd)+1
               if(appo.lt.0.) k_index=k_index-1
               if (k_index.eq.0) then
                  k_index=1
               end if
               kk_index=k_index+1
               if (k_index.eq.nzp) kk_index=1
               z_left=zl*real(k_index-nzp/2-1)/real(nzp)+zsc
c
c     Compute fluid velocity in surrounding box
c
               appo=real(ii_index)

               appo1=i_index
               appo2=mod(i_index,nxp/2)+1

               uf_1=u2r(appo1,k_index,1)*appo
     &              +u2i(appo1,k_index,1)*(1.-appo)
               uf_2=u2i(appo1,k_index,1)*appo
     &              +u2r(appo2,k_index,1)*(1.-appo)
               uf_3=u2r(appo1,kk_index,1)*appo
     &              +u2i(appo1,kk_index,1)*(1.-appo)
               uf_4=u2i(appo1,kk_index,1)*appo
     &              +u2r(appo2,kk_index,1)*(1.-appo)

               vf_1=u2r(appo1,k_index,2)*appo
     &              +u2i(appo1,k_index,2)*(1.-appo)
               vf_2=u2i(appo1,k_index,2)*appo
     &              +u2r(appo2,k_index,2)*(1.-appo)
               vf_3=u2r(appo1,kk_index,2)*appo
     &              +u2i(appo1,kk_index,2)*(1.-appo)
               vf_4=u2i(appo1,kk_index,2)*appo
     &              +u2r(appo2,kk_index,2)*(1.-appo)

               wf_1=u2r(appo1,k_index,3)*appo
     &              +u2i(appo1,k_index,3)*(1.-appo)
               wf_2=u2i(appo1,k_index,3)*appo
     &              +u2r(appo2,k_index,3)*(1.-appo)
               wf_3=u2r(appo1,kk_index,3)*appo
     &              +u2i(appo1,kk_index,3)*(1.-appo)
               wf_4=u2i(appo1,kk_index,3)*appo
     &              +u2r(appo2,kk_index,3)*(1.-appo)
c
c     Interpolate 1 and 2 to get 5
c
               uf_5=(uf_2-uf_1)*(appo_x-x_left)/xl*real(nxp)+uf_1
               vf_5=(vf_2-vf_1)*(appo_x-x_left)/xl*real(nxp)+vf_1
               wf_5=(wf_2-wf_1)*(appo_x-x_left)/xl*real(nxp)+wf_1
c
c     Interpolate 3 and 4 to get 6
c
               uf_6=(uf_4-uf_3)*(appo_x-x_left)/xl*real(nxp)+uf_3
               vf_6=(vf_4-vf_3)*(appo_x-x_left)/xl*real(nxp)+vf_3
               wf_6=(wf_4-wf_3)*(appo_x-x_left)/xl*real(nxp)+wf_3
c
c     Interpolate 5 and 6 to get 7
c
               uf_7=(uf_6-uf_5)*(PINT(3,p,1)-z_left)/zl*real(nzp)+uf_5
               vf_7=(vf_6-vf_5)*(PINT(3,p,1)-z_left)/zl*real(nzp)+vf_5
               wf_7=(wf_6-wf_5)*(PINT(3,p,1)-z_left)/zl*real(nzp)+wf_5
c
c     Store interpolated velocity
c
               PINT(1,p,nn_index)=uf_7
               PINT(2,p,nn_index)=vf_7
               PINT(3,p,nn_index)=wf_7
            end if

         end do

      end if
c
c________________________________________________________________________
c

c
c     Compute CFL number
c
      if (mod(it-1,icfl).eq.0) then
c
c     Comment: for perturbation formulation it seems
c     unnecessary to compute the clf based on baseflow at
c     each time step. However it was left that way since:
c     - support time-varying base flows
c     - computational overhead is minimal
c
         if (pert) then
            if (.not.bf3) then
               call boxcflBF(cflp_bf,bu1,bu2,yb,deta,xl,zl)
            else

               do i = 1,3
                  do z = 1,nzd
                     do x = 1,nxp/2+1
                        bf3tempr(x,z,i) = bf3u2r(x,z,ybp,i)
                        bf3tempi(x,z,i) = bf3u2i(x,z,ybp,i)
                     end do
                  end do
               end do

               call boxcfl(cflp_bf, bf3tempr,bf3tempi,
     &              yb,deta,xl,zl,wbci)
            end if
         else
            cflp_bf=0.
         end if
         if (.not.lin) then
            call boxcfl(cflp,u2r,u2i,yb,deta,xl,zl,wbci)
         else
            cflp=0.
         end if
         cflp=cflp+cflp_bf
      end if
c
c     Compute extrema
c
      if (mod(it-1,iext).eq.0.and.iext.gt.0) call
     &     boxext(vext,cext,u2r,u2i,om2r,om2i,xs,yb,xl,zl)
c
c     Compute the nonlinear terms in physical space
c
      call advection(pert,nxy,rot,om2r,om2i,u2r,u2i,
     &     th2r,th2i,lin,bu1,bu2,bom1,bom2,yb,
     &     bf3,bf3u2r,bf3u2i,ybp,pr,gr,fltype,dstar)
c
c     Fringe region for spatial simulations (or temporal simulations
c     with fringe, like fltype 9)
c
      if (spat) then
         call fring(om2r,om2i,u2r,u2i,tc,xsc,zsc,xl,zl,yb,
     &        fstart,fend,bu1,bu2,
     &        osmod,osnumb,
     &        osur,osui,osvr,osvi,oswr,oswi,
     &        evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &        ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &        xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &        th2r,th2i,
     &        ampst,streak,betast,omegast,
     &        ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &        tsmoo,tsmst,tsmend,iampst,phist,pert)
      end if
c
c     Localized volume force
c
      if (loctyp.ge.1) then
         call locf(om2r,om2i,yb,xl,zl,xsc,zsc,eta,tc,loctyp,
     &        fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,fpds8,
     &        fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,wr,wi,
     &        th2r,th2i,u2r,u2i)
      end if
c
c     Trip forcing (on the velocities)
c
      if (tripf) then
         call trip(om2r,om2i,yb,xl,xsc,eta,txsc,tysc,ttzc,tx0)
      end if
c
c     Forward Fourier transform of nonlinear term H_i
c
      do i = 1,3
         sym = i.le.2
         call fft2df(om2r(1,1,i),om2i(1,1,i),sym,npl,
     &        prex,prez,pres,prea,wr,wi)
      end do

      do ith=1,scalar
         sym = .true.
         call fft2df(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1)),sym,
     &        npl,prex,prez,pres,prea,wr,wi)
      end do
c
c     Remove the oddball modes in z on the fine grid
c     (if necessary, i.e. nzp even)
c
      if (mod(nzp,2).eq.0) then
         do i=1,3
            do xy=1,nxy
               om2r(xy,nzp/2+1,i)=0.0
               om2i(xy,nzp/2+1,i)=0.0
            end do
         end do
         do ith=1,scalar
            do xy=1,nxy
               th2r(xy,nzp/2+1,4+4*(ith-1))=0.0
               th2i(xy,nzp/2+1,4+4*(ith-1))=0.0
            end do
         end do
      end if
c
c     Now the nonlinear term is in Fourier/Fourier/Real space.
c     Note that everything is multiplyed by nxp*nzp.
c
c     Adding a volume force to wavenumber zero (mean component)
c     to get the correct boundary layer growth for parallel flows.
c     Note that we are in Fourier-Chebyshev space but have not yet
c     normalized for the Fourier transform.
c     This means that the force has to be multiplied by nxp*nzp.
c     The time-independent forcing term (f1,0,f3) on the right hand side
c     of the Navier-Stokes equations needs to preserve the initial
c     parallel profile if no disturbances are present.
c     f1 = c*dU/dx-1/Re*d2U/d2y =
c        = c*(m/x0)*(x/x0)^(m-1)*f'+c*(x/x0)^m*f''deta/dx-f'''/(2*x0)
c     f3 = -spanv/Re*d2W/d2y=-spanv*(m+1)/(2*x)*d2W/d2y
c     f is similarity solution of
c     f'''+f*f''+2*m/(m+1)*(1-f'*F')=0
c     g''+fg'=0
c     where
c     f'=df/ds
c     s=y*sqrt((m+1)*U_0/(2*nu*x0)), x0 inflow location from the leading edge,
c     scaled by displacement thickness, c a reference
c     speed (usually a group velocity). These terms are easily computed with
c     the nonlinear terms over here. There is another choice, which is time-
c     dependent forcing (see Spalart & Yang, 1987, JFM).
c
      if (.not.pert) then
c
c     Forcing only if temporal boundary layer flows are considered
c
         if (fltype.eq. 3.or.fltype.eq.9.or.
     &       fltype.eq.-1.or.fltype.eq.-2) then
            c1=real(nxp)*real(nzp)
            xx=(xblc/x0)**rlam
            do y=1,npl
               ybl=1.+eta(y+yb-1)
               etabl=ybl*sqrt((rlam+1.)*re*xx/(2.*xblc))
               h1u = cdev*(rlam/x0)*(xblc/x0)**(rlam-1)
               h2u = cdev*0.5*(rlam-1)*xx**1.5*sqrt(0.5*(rlam+1.)*re)*
     &               (1./xblc)**1.5*ybl
               h3u = 0.5*(rlam+1.)/xblc
               d1f=cubip(etabl,fbla(1,2),dybla,nbla)
               d2f=cubip(etabl,fbla(1,3),dybla,nbla)
               d3f=cubip(etabl,fbla(1,4),dybla,nbla)
               om2r(1+(y-1)*(nxp/2+1),1,1)=om2r(1+(y-1)*(nxp/2+1),1,1)+
     &              (d1f*h1u + d2f*h2u - d3f*xx**2*h3u)*c1
               if (fltype.eq.-2.or.fltype.eq.9) then
                  d2g=cubip(etabl,fbla(1,7),dybla,nbla)
                  om2r(1+(y-1)*(nxp/2+1),1,3)=om2r(1+(y-1)*
     &                 (nxp/2+1),1,3)-spanv*h3u*d2g*c1
               end if
            end do
         end if
      end if
c
c     Additional forcing terms
c     Note that these could also be added in the linear step since
c     they do not involve nonlinear products. Thereby, a communication
c     step could be avoided making everything a bit more efficient.
c
c
c-----------------------------------------------------------------
c
c     MHD Lorentz force
c
      if (imhd.eq.1) then
c
c     Reuse du2r array for the forcing term
c
         do ll=1,2
            if (nproc.eq.1) then
               call getxz(du2r(1,1,ll,1),du2i(1,1,ll,1),
     &              yb,ll,1,fmhdr,fmhdi)
            else
#ifdef MPI
               call getpxz(du2r(1,1,ll,1),du2i(1,1,ll,1),yb,ll,1,
     &              fmhdr,fmhdi,realg1,realg2,my_node)
#endif
            end if
         end do
c
c     Add the MHD force to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
         fact = mhd_n*nxp*nzp
         do z=1,nz/2
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1

                  cr1 =  alfa(x)*du2i(xy,z,1,1)
     &                 +du2r(xy,z,2,2)*b0(3)-du2r(xy,z,3,2)*b0(2)
                  ci1 = -alfa(x)*du2r(xy,z,1,1)
     &                 +du2i(xy,z,2,2)*b0(3)-du2i(xy,z,3,2)*b0(2)

                  cr2 = -du2r(xy,z,2,1)
     &                 +du2r(xy,z,3,2)*b0(1)-du2r(xy,z,1,2)*b0(3)
                  ci2 = -du2i(xy,z,2,1)
     &                 +du2i(xy,z,3,2)*b0(1)-du2i(xy,z,1,2)*b0(3)

                  cr3 =  beta(z)*du2i(xy,z,1,1)
     &                 +du2r(xy,z,1,2)*b0(2)-du2r(xy,z,2,2)*b0(1)
                  ci3 = -beta(z)*du2r(xy,z,1,1)
     &                 +du2i(xy,z,1,2)*b0(2)-du2i(xy,z,2,2)*b0(1)
c
c     Lorentz force
c
                  om2r(xy,z,1) = om2r(xy,z,1) +
     &                 fact*(cr2*b0(3) - cr3*b0(2))
                  om2i(xy,z,1) = om2i(xy,z,1) +
     &                 fact*(ci2*b0(3) - ci3*b0(2))

                  om2r(xy,z,2) = om2r(xy,z,2) +
     &                 fact*(cr3*b0(1) - cr1*b0(3))
                  om2i(xy,z,2) = om2i(xy,z,2) +
     &                 fact*(ci3*b0(1) - ci1*b0(3))

                  om2r(xy,z,3) = om2r(xy,z,3) +
     &                 fact*(cr1*b0(2) - cr2*b0(1))
                  om2i(xy,z,3) = om2i(xy,z,3) +
     &                 fact*(ci1*b0(2) - ci2*b0(1))


                  end do
               end do
            end do
c
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1

                  cr1 =  alfa(x)*du2i(xy,zp,1,1)
     &                 +du2r(xy,zp,2,2)*b0(3)-du2r(xy,zp,3,2)*b0(2)
                  ci1 = -alfa(x)*du2r(xy,zp,1,1)
     &                 +du2i(xy,zp,2,2)*b0(3)-du2i(xy,zp,3,2)*b0(2)

                  cr2 = -du2r(xy,zp,2,1)
     &                 +du2r(xy,zp,3,2)*b0(1)-du2r(xy,zp,1,2)*b0(3)
                  ci2 = -du2i(xy,zp,2,1)
     &                 +du2i(xy,zp,3,2)*b0(1)-du2i(xy,zp,1,2)*b0(3)

                  cr3 =  beta(z)*du2i(xy,zp,1,1)
     &                 +du2r(xy,zp,1,2)*b0(2)-du2r(xy,zp,2,2)*b0(1)
                  ci3 = -beta(z)*du2r(xy,zp,1,1)
     &                 +du2i(xy,zp,1,2)*b0(2)-du2i(xy,zp,2,2)*b0(1)
c
c     Lorentz force
c
                  om2r(xy,zp,1) = om2r(xy,zp,1) +
     &                 fact*(cr2*b0(3) - cr3*b0(2))
                  om2i(xy,zp,1) = om2i(xy,zp,1) +
     &                 fact*(ci2*b0(3) - ci3*b0(2))

                  om2r(xy,zp,2) = om2r(xy,zp,2) +
     &                 fact*(cr3*b0(1) - cr1*b0(3))
                  om2i(xy,zp,2) = om2i(xy,zp,2) +
     &                 fact*(ci3*b0(1) - ci1*b0(3))

                  om2r(xy,zp,3) = om2r(xy,zp,3) +
     &                 fact*(cr1*b0(2) - cr2*b0(1))
                  om2i(xy,zp,3) = om2i(xy,zp,3) +
     &                 fact*(ci1*b0(2) - ci2*b0(1))

                  end do
               end do
            end do
      end if
c
c
c--------------------------------------------------------------
c
c     Large-eddy simulation ( LES )
c
      if (iles.eq.1) then
c
c     ADM-RT model:
c     Get filtered velocities
c
         do ll=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui)
            else
#ifdef MPI
               call getpxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     Add the relaxation term to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     om2r(xy,z,ll) = om2r(xy,z,ll) -
     &                    chi_scaled*u2r(xy,z,1)*nxp*nzp
                     om2i(xy,z,ll) = om2i(xy,z,ll) -
     &                    chi_scaled*u2i(xy,z,1)*nxp*nzp
                  end do
               end do
            end do
c
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     om2r(xy,zp,ll) = om2r(xy,zp,ll) -
     &                    chi_scaled*u2r(xy,zp,1)*nxp*nzp
                     om2i(xy,zp,ll) = om2i(xy,zp,ll) -
     &                    chi_scaled*u2i(xy,zp,1)*nxp*nzp
                  end do
               end do
            end do
         end do

         do ith=1,scalar
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi)
            else
#ifdef MPI
               call getpxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     Add the relaxation term to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     th2r(xy,z,4*(ith-1)+4) = th2r(xy,z,4*(ith-1)+4) -
     &                    chi_scaled*u2r(xy,z,1)*nxp*nzp/prt
                     th2i(xy,z,4*(ith-1)+4) = th2i(xy,z,4*(ith-1)+4) -
     &                    chi_scaled*u2i(xy,z,1)*nxp*nzp/prt
                  end do
               end do
            end do
c
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c           This would however break nz=1 (2D) simulations
c
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     th2r(xy,zp,4*(ith-1)+4) = th2r(xy,zp,4*(ith-1)+4) -
     &                    chi_scaled*u2r(xy,zp,1)*nxp*nzp/prt
                     th2i(xy,zp,4*(ith-1)+4) = th2i(xy,zp,4*(ith-1)+4) -
     &                    chi_scaled*u2i(xy,zp,1)*nxp*nzp/prt
                  end do
               end do
            end do
         end do
      else if (iles.eq.2.or.iles.eq.3) then
c
c     (HPF) Eddy-viscosity model:
c     get the SGS force
c
         do ll=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui)
            else
#ifdef MPI
               call getpxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     Add the SGS force to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     om2r(xy,z,ll) = om2r(xy,z,ll) -
     &                    u2r(xy,z,1)*nxp*nzp
                     om2i(xy,z,ll) = om2i(xy,z,ll) -
     &                    u2i(xy,z,1)*nxp*nzp
                  end do
               end do
            end do
c
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     om2r(xy,zp,ll) = om2r(xy,zp,ll) -
     &                    u2r(xy,zp,1)*nxp*nzp
                     om2i(xy,zp,ll) = om2i(xy,zp,ll) -
     &                    u2i(xy,zp,1)*nxp*nzp
                  end do
               end do
            end do
         end do
c
c     (HPF) Eddy-diffusivity model for the scalars
c
         do ith=1,scalar

            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi)
            else
#ifdef MPI
               call getpxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     Add the SGS force to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     th2r(xy,z,4*(ith-1)+4) = th2r(xy,z,4*(ith-1)+4) -
     &                    u2r(xy,z,1)*nxp*nzp
                     th2i(xy,z,4*(ith-1)+4) = th2i(xy,z,4*(ith-1)+4) -
     &                    u2i(xy,z,1)*nxp*nzp
                  end do
               end do
            end do
c
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c           This would however break nz=1 (2D) simulations
c
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     th2r(xy,zp,4*(ith-1)+4) = th2r(xy,zp,4*(ith-1)+4) -
     &                    u2r(xy,zp,1)*nxp*nzp
                     th2i(xy,zp,4*(ith-1)+4) = th2i(xy,zp,4*(ith-1)+4) -
     &                    u2i(xy,zp,1)*nxp*nzp
                  end do
               end do
            end do
         end do
      end if
c
c     SFD (Selective Frequency Damping)
c
c     Memory layout:
c     sfdur(:,:,:,1-3): temporally filtered velocity field
c     sfdur(:,:,:,4-6): unfiltered-filtered velocity field (RHS^n-1)
c
c     u2r(:,:,1): temporally filtered velocity
c     u2r(:,:,2): unfiltered velocity
c     u2r(:,:,3): unfiltered-filtered velocity (RHS^n)
c     u2r(:,:,1): new filtered velocity
c     u2r(:,:,2): RHS^n-1
c
c     Get filtered and unfiltered velocities into u2r(:,:,1) and u2r(:,:,2)
c
      if (isfd.eq.1) then
         do ll=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,sfdur,sfdui)
               call getxz(u2r(1,1,2),u2i(1,1,2),yb,ll,1,   ur,   ui)

            else
#ifdef MPI
               call getpxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,sfdur,sfdui,
     &              realg1,realg2,my_node)
               call getpxz(u2r(1,1,2),u2i(1,1,2),yb,ll,1,ur,ui,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     Add the SFD relaxation term to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,z,3) = u2r(xy,z,2)-u2r(xy,z,1)
                     u2i(xy,z,3) = u2i(xy,z,2)-u2i(xy,z,1)

                     om2r(xy,z,ll) = om2r(xy,z,ll) -
     &                    sfd_chi*u2r(xy,z,3)*nxp*nzp
                     om2i(xy,z,ll) = om2i(xy,z,ll) -
     &                    sfd_chi*u2i(xy,z,3)*nxp*nzp
                  end do
               end do
            end do
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,zp,3) = u2r(xy,zp,2)-u2r(xy,zp,1)
                     u2i(xy,zp,3) = u2i(xy,zp,2)-u2i(xy,zp,1)

                     om2r(xy,zp,ll) = om2r(xy,zp,ll) -
     &                    sfd_chi*u2r(xy,zp,3)*nxp*nzp
                     om2i(xy,zp,ll) = om2i(xy,zp,ll) -
     &                    sfd_chi*u2i(xy,zp,3)*nxp*nzp
                  end do
               end do
            end do
c
c     Do the RK integration (explicit) of the filtered term
c     Add G^n in u2r(:,:,1)
c
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,z,1) = u2r(xy,z,1) +
     &                    an/sfd_delta*u2r(xy,z,3)
                     u2i(xy,z,1) = u2i(xy,z,1) +
     &                    an/sfd_delta*u2i(xy,z,3)
                  end do
               end do
            end do
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,zp,1) = u2r(xy,zp,1) +
     &                    an/sfd_delta*u2r(xy,zp,3)
                     u2i(xy,zp,1) = u2i(xy,zp,1) +
     &                    an/sfd_delta*u2i(xy,zp,3)
                  end do
               end do
            end do
            if (bn.ne.0) then
               ! get G^n-1 onto u2r(:,:,2)
               if (nproc.eq.1) then
                  call getxz(u2r(1,1,2),u2i(1,1,2),yb,ll+3,1,
     &                 sfdur,sfdui)
               else
#ifdef MPI
                  call getpxz(u2r(1,1,2),u2i(1,1,2),yb,ll+3,1,
     &                 sfdur,sfdui,realg1,realg2,my_node)
#endif
               end if
               ! add G^n-1
               do z=1,nz/2
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1
                        u2r(xy,z,1) = u2r(xy,z,1) +
     &                       bn/sfd_delta*u2r(xy,z,2)
                        u2i(xy,z,1) = u2i(xy,z,1) +
     &                       bn/sfd_delta*u2i(xy,z,2)
                     end do
                  end do
               end do
               do z=nz/2+1,nz
                  zp=nzp-nz+z
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1
                        u2r(xy,zp,1) = u2r(xy,zp,1) +
     &                       bn/sfd_delta*u2r(xy,zp,2)
                        u2i(xy,zp,1) = u2i(xy,zp,1) +
     &                       bn/sfd_delta*u2i(xy,zp,2)
                     end do
                  end do
               end do
            end if
c
c     Put G^n and filtered u
c
            if (nproc.eq.1) then
               call putxz(u2r(1,1,1),u2i(1,1,1),yb,ll  ,sfdur,sfdui)
               call putxz(u2r(1,1,3),u2i(1,1,3),yb,ll+3,sfdur,sfdui)
            else
#ifdef MPI
               call putpxz(u2r(1,1,1),u2i(1,1,1),yb,ll  ,sfdur,sfdui,
     &              realg1,realg2,my_node)
               call putpxz(u2r(1,1,3),u2i(1,1,3),yb,ll+3,sfdur,sfdui,
     &              realg1,realg2,my_node)
#endif
            end if
         end do
      end if
c
c     Put the planes of the nonlinear term back onto the velocities
c     (here, the truncation to the normal grid occurs)
c     NOTE: THE ODDBALL MODES ARE STILL IN THERE!!!!
c

c      call mpi_barrier(mpi_comm_world,ierror)
c      call wall_time(t2)
      call ctim(ct2,wt2)
      wtser = wtser + (wt2-wt1)
      ctser = ctser + (ct2-ct1)
      ct1=ct2
      wt1=wt2
      

      if (nproc.eq.1) then
         do i=1,3
            call putxz(om2r(1,1,i),om2i(1,1,i),yb,i,ur,ui)
         end do
         do ith=1,scalar
            call putxz(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1)),yb,
     &           8+pressure+3*(ith-1),ur,ui)
         end do
      else
#ifdef MPI
         do i=1,3
            call putpxz(om2r(1,1,i),om2i(1,1,i),yb,i,ur,ui,
     &           realg1,realg2,my_node)
         end do
         do ith=1,scalar
            call putpxz(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1)),yb,
     &           8+pressure+3*(ith-1),ur,ui,
     &           realg1,realg2,my_node)
         end do
#endif
      end if

c      call mpi_barrier(mpi_comm_world,ierror)
c      call wall_time(t2)
      call ctim(ct2,wt2)
      wtcom = wtcom + (wt2-wt1)
      ctcom = ctcom + (ct2-ct1)
      ct1=ct2
      wt1=wt2
      
      end subroutine nonlinbl
