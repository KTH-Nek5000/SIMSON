c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine updatewallrough(ur,ui,xsc,xst,xlen,delx,
     &     prexn,prey,prezn,presn,prean,alfa,beta,h_rough,xstart,xend,
     &     xrise,xfall,nzr,rzd,tpsi,v_wall,taylor4,my_node,wallur,
     &     wallui,wlvr,wlvi,wallwr,wallwi)
c
c     Calculates the boundary condition at y=0 resulting from
c     surface roughness. The latter is modelled in terms of a Taylor
c     series (projection of the no-slip conditions from the desired
c     bump contour to the "wall" y=0)
c
c     Taylor series is truncated at 1st order or at 4th order (if
c     taylor4==T)
c
c     The projected wall-boundary conditions are calculated using
c     the CURRENT FLOW FIELD (ur,ui). This routine is called
c     every EVERY-th time step, see bla.f
c
c     Called by bla.f (if wbci==-1, updat==T)
c
c     Calls getxz/getpxz,vchbf,vchbb,dcheb,rdcheb,fft2df,fft2db
c
      implicit none

      include 'par.f'

#ifdef MPI
      include 'mpif.h'
#endif

c
c     Parameters
c
      integer nxm,nxn,nzn,lenx
      parameter (nxm=nxp/2+1,nxn=nx/2+1,nzn=nz/2+1)
      parameter (lenx=nx/8)

      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Global variables
c
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer xst,xlen,nzr,my_node

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15)
      real prean(nz*3/4+15),prey(nyp*2+15),w(nxm,nyp)
      real alfa(nx/2*mbz),beta(nz)
      real rzd(nz+2)
      real delx,xsc,tpsi
      real h_rough,xstart,xend,xrise,xfall
      logical v_wall,taylor4
c
c     Output
c
      real wallur(nxm,memnz),wallui(nxm,memnz)
      real wlvr(nxm,memnz),wlvi(nxm,memnz)
      real wallwr(nxm,memnz),wallwi(nxm,memnz)
c
c     Local variables
c
      integer i,x,xsh,xind,y,z,zb,zst,zen,npl,nxy
      integer halfproc
      real delty,step,xc11,xc22,norm
      real tmpr1,tmpr2
c      real zshift1,zshift2

      real fx1,fx2
      real hfunc(lenx,memnz,2)

      real wr(nxm*mby,nzd),wi(nxm*mby,nzd)
      real u2r(nxm*mby,nzd,3),u2i(nxm*mby,nzd,3)
      real urvec(lenx,nyp,memnz),uivec(lenx,nyp,memnz)

#ifdef MPI
      integer nsr,ierror
      real wyr(nxn,nyp),wyi(nxn,nyp)
      real u3r(nxn,nyp),u3i(nxn,nyp)
      real urv(lenx,nyp,nz),uiv(lenx,nyp,nz)
      real rwallur(nxn*mby,nz),rwallui(nxn*mby,nz)
      real rwallvr(nxn*mby,nz),rwallvi(nxn*mby,nz)
      real rwallwr(nxn*mby,nz),rwallwi(nxn*mby,nz)
#endif

      real app1(nyp),app2(nyp),app21(nyp),app22(nyp)
      real app31(nyp),app32(nyp),app41(nyp),app42(nyp)

      real du1wall(lenx,memnz,3),du2wall(lenx,memnz,3)
      real d2u1wall(lenx,memnz,3),d2u2wall(lenx,memnz,3)
      real d3u1wall(lenx,memnz,3),d3u2wall(lenx,memnz,3)
      real d4u1wall(lenx,memnz,3),d4u2wall(lenx,memnz,3)

      real lwallur(nxn*mby,memnz),lwallui(nxn*mby,memnz)
      real lwallvr(nxn*mby,memnz),lwallvi(nxn*mby,memnz)
      real lwallwr(nxn*mby,memnz),lwallwi(nxn*mby,memnz)


      if (my_node.eq.0) then
         write(ios,*)'--> Wall-roughness B.C.s are updated.'
      end if
c
c     Get the flow field and transform it back to physical space
c
      do i=1,3
         if (nproc.eq.1) then
            do y=1,nyp
               npl=min(mby,nyp-y+1)
               nxy=nxm*npl-nxm+nx/2

               call getxz(u2r(1,1,i),u2i(1,1,i),y,i,0,ur,ui)
c
c     Chordwise backward fft
c
               call vrfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,
     &              nx,nzc*mby,1,nxm,prexn)
c
c     Spanwise backward fft
c
               if (nfzsym.eq.0) then
                  call vcfftb(u2r(1,1,i),u2i(1,1,i),wr,wi,nz,
     &                 nxy,nxm*mby,1,prezn)
               else
                  if (i.le.2) then
                     call vcffts(u2r(1,1,i),u2i(1,1,i),wr,
     &                    nzn,nxy,nxm*mby,1,presn)
                  else
                     call vcftab(u2r(1,2,i),u2i(1,2,i),
     &                    wr,wi,nz/2-1,nxy,nxm*mby,1,prean)
                  end if
               end if

               do z=1,nz
                  do x=1,xlen
                     xind=x+xst-1
                     urvec(x,y,z)=u2r(xind,z,i)
                     uivec(x,y,z)=u2i(xind,z,i)
                  end do
               end do
            end do

         else

#ifdef MPI
            do z=1,memnz

               call getfxy(u3r,u3i,z,i,ur,ui)
c
c     Chordwise backward fft
c
               call vrfftb(u3r,u3i,wyr,wyi,nx,nyp,1,nxn,prexn)

               do y=1,nyp
                  do x=1,xlen
                     xind=x+xst-1
                     urvec(x,y,z)=u3r(xind,y)
                     uivec(x,y,z)=u3i(xind,y)
                  end do
               end do
            end do
c
c     The root process collects the "vel.-field stripe", ...
c
            if (nproc.gt.1) then
               nsr=lenx*nyp*memnz
               call mpi_gather(urvec,nsr,mpi_double_precision,
     &              urv,nsr,mpi_double_precision,0,
     &              mpi_comm_world,ierror)
               call mpi_gather(uivec,nsr,mpi_double_precision,
     &              uiv,nsr,mpi_double_precision,0,
     &              mpi_comm_world,ierror)
            end if
c
c     ... performs the spanwise backward fft ...
c
            if (my_node.eq.0) then
               do y=1,nyp

                  do x=1,xlen
                     urv(x,y,nz/2+1)=0.0
                     uiv(x,y,nz/2+1)=0.0
                  end do

                  if (nfzsym.eq.0) then
                     call vcfftb(urv(:,y,:),uiv(:,y,:),wr,wi,
     &                    nz,xlen,lenx,1,prezn)
                  else
                     if (i.le.2) then
                        call vcffts(urv(:,y,:),uiv(:,y,:),wr,
     &                       nzn,xlen,lenx,1,presn)
                     else
                        call vcftab(urv(:,y,:),uiv(:,y,:),wr,wi,
     &                       nz/2-1,xlen,lenx,1,prean)
                     end if
                  end if

               end do

            end if
c
c     ... and scatters the data back to the partner processes.
c
            if (nproc.gt.1) then
               call mpi_scatter(urv,nsr,mpi_double_precision,
     &              urvec,nsr,mpi_double_precision,0,
     &              mpi_comm_world,ierror)
               call mpi_scatter(uiv,nsr,mpi_double_precision,
     &              uivec,nsr,mpi_double_precision,0,
     &              mpi_comm_world,ierror)
            end if

#endif

         end if
c
c     Calculation of the vel. derivatives at y=0
c
         delty = 2./real(nyp-1) ! y step width (internal scaling)

         do z=1,memnz
            do x=1,xlen
               do y=1,nyp
                  app1(y)=urvec(x,y,z)
                  app2(y)=uivec(x,y,z)
               end do

               call vchbf(app1,w,nyp,1,1,1,prey)
               call vchbf(app2,w,nyp,1,1,1,prey)
c
c     Normalise
c
               do y=1,nyp
                  app1(y)=app1(y)*delty
                  app2(y)=app2(y)*delty
               end do

               call rdcheb(app1,nyp,1,1) !1st deriv. of vel. profs.
               call rdcheb(app2,nyp,1,1)

               if (taylor4) then

                  call dcheb(app21,app1,nyp,1,1) ! 2nd derivative
                  call dcheb(app22,app2,nyp,1,1)

                  call dcheb(app31,app21,nyp,1,1) ! 3rd derivative
                  call dcheb(app32,app22,nyp,1,1)

                  call dcheb(app41,app31,nyp,1,1) ! 4th derivative
                  call dcheb(app42,app32,nyp,1,1)


                  call vchbb(app21,w,nyp,1,1,1,prey)
                  call vchbb(app22,w,nyp,1,1,1,prey)

                  call vchbb(app31,w,nyp,1,1,1,prey)
                  call vchbb(app32,w,nyp,1,1,1,prey)

                  call vchbb(app41,w,nyp,1,1,1,prey)
                  call vchbb(app42,w,nyp,1,1,1,prey)


                  d2u1wall(x,z,i)=app21(nyp) !2nd derivative
                  d2u2wall(x,z,i)=app22(nyp)

                  d3u1wall(x,z,i)=app31(nyp) !3rd derivative
                  d3u2wall(x,z,i)=app32(nyp)

                  d4u1wall(x,z,i)=app41(nyp) !4th derivative
                  d4u2wall(x,z,i)=app42(nyp)

               end if


               call vchbb(app1,w,nyp,1,1,1,prey)
               call vchbb(app2,w,nyp,1,1,1,prey)

               du1wall(x,z,i)=app1(nyp) !1st derivative
               du2wall(x,z,i)=app2(nyp)
            end do
         end do

      end do ! of the i-loop
c
c     Initialize the wall B.C.s
c
      lwallur(:,:)=0.
      lwallui(:,:)=0.
      lwallvr(:,:)=0.
      lwallvi(:,:)=0.
      lwallwr(:,:)=0.
      lwallwi(:,:)=0.

      tpsi = 0.! R E M O V E
c
c     The bump function and the projected wall BCs
c
      xsh = xst-0.5*nxn
      if (xsh.le.0) xsh = xsh+nxn

      do z=1,memnz
         zb=my_node*memnz+z
         if (nzr.eq.1) then
            tmpr1=1.
            tmpr2=1.
         else
            tmpr1 = rzd(zb)
            tmpr2 = tmpr1
         end if

         do x=1,xlen
            xind=x+xsh-1

            xc11=xstart+(2*x-2)*delx+xsc
            fx1=step((xc11-xstart)/xrise)-
     &           step((xc11-xend)/xfall+1)

            xc22 = xstart+(2*x-1)*delx+xsc
            fx2=step((xc22-xstart)/xrise)-
     &           step((xc22-xend)/xfall+1)

c            zshift1 = tpsi*(xc11-xstart)
c            zshift2 = tpsi*(xc22-xstart)

            hfunc(x,z,1)=h_rough*fx1*tmpr1
            hfunc(x,z,2)=h_rough*fx2*tmpr2

            lwallur(xind,z)=-hfunc(x,z,1)*du1wall(x,z,1) !"roughness
            lwallui(xind,z)=-hfunc(x,z,2)*du2wall(x,z,1) !projection"

            lwallwr(xind,z)=-hfunc(x,z,1)*du1wall(x,z,3)
            lwallwi(xind,z)=-hfunc(x,z,2)*du2wall(x,z,3)

            if (taylor4) then !Take 2nd, 3rd, 4th deriv. into account
               lwallur(xind,z)=lwallur(xind,z)-
     &              0.5*(hfunc(x,z,1)**2)*d2u1wall(x,z,1)-
     &              (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,z,1)-
     &              (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,z,1)
               lwallui(xind,z)=lwallui(xind,z)-
     &              0.5*(hfunc(x,z,2)**2)*d2u2wall(x,z,1)-
     &              (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,z,1)-
     &              (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,z,1)

               lwallwr(xind,z)=lwallwr(xind,z)-
     &              0.5*(hfunc(x,z,1)**2)*d2u1wall(x,z,3)-
     &              (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,z,3)-
     &              (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,z,3)
               lwallwi(xind,z)=lwallwi(xind,z)-
     &              0.5*(hfunc(x,z,2)**2)*d2u2wall(x,z,3)-
     &              (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,z,3)-
     &              (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,z,3)
            end if

            if (v_wall) then

               lwallvr(xind,z)=-hfunc(x,z,1)*du1wall(x,z,2)
               lwallvi(xind,z)=-hfunc(x,z,2)*du2wall(x,z,2)

               if (taylor4) then !Take 2nd, 3rd, 4th deriv. into account
                  lwallvr(xind,z)=lwallvr(xind,z)-
     &                 0.5*(hfunc(x,z,1)**2)*d2u1wall(x,z,2)-
     &                 (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,z,2)-
     &                 (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,z,2)
                  lwallvi(xind,z)=lwallvi(xind,z)-
     &                 0.5*(hfunc(x,z,2)**2)*d2u2wall(x,z,2)-
     &                 (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,z,2)-
     &                 (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,z,2)
               end if
            end if
         end do
      end do
c
c     Real forward fft in x direction
c
      call vrfftf(lwallur,lwallui,wr,wi,nx,memnz*mby,1,nxn,prexn)
      call vrfftf(lwallwr,lwallwi,wr,wi,nx,memnz*mby,1,nxn,prexn)
      if (v_wall) call vrfftf(lwallvr,lwallvi,wr,wi,nx,memnz*mby,1,
     &     nxn,prexn)

      if (nproc.eq.1) then
c
c     Complex forward fft in z direction
c
         if (nfzsym.eq.0) then
            call vcfftf(lwallur,lwallui,wr,wi,nz,nx/2,nxn*mby,1,prezn)
            call vcfftf(lwallwr,lwallwi,wr,wi,nz,nx/2,nxn*mby,1,prezn)
            if (v_wall) call vcfftf(lwallvr,lwallvi,wr,wi,
     &           nz,nx/2,nxn*mby,1,prezn)
         else
            call vcffts(lwallur,lwallui,wr,nzn,nx/2,
     &           nxn*mby,1,presn)
            call vcftaf(lwallwr,lwallwi,wr, wi,nz/2-1,nx/2,
     &           nxn*mby,1,prean)
            if (v_wall) call vcffts(lwallvr,lwallvi,wr,
     &           nzn,nx/2,nxn*mby,1,presn)
         end if
      end if

#ifdef MPI
c
c     Collect the boundary planes at process 0 to perform the spanwise fft
c
      if (nproc.gt.1) then

         nsr=nxn*memnz
         call mpi_gather(lwallur,nsr,mpi_double_precision,
     &      rwallur,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_gather(lwallui,nsr,mpi_double_precision,
     &      rwallui,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_gather(lwallwr,nsr,mpi_double_precision,
     &      rwallwr,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_gather(lwallwi,nsr,mpi_double_precision,
     &      rwallwi,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         if (v_wall) then
            call mpi_gather(lwallvr,nsr,mpi_double_precision,
     &      rwallvr,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
            call mpi_gather(lwallvi,nsr,mpi_double_precision,
     &      rwallvi,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         end if

         if (my_node.eq.0) then
c
c     Complex forward fft in z direction
c
            if (nfzsym.eq.0) then
               call vcfftf(rwallur,rwallui,wr,wi,nz,nx/2,
     &              nxn*mby,1,prezn)
               call vcfftf(rwallwr,rwallwi,wr,wi,nz,nx/2,
     &              nxn*mby,1,prezn)
               if (v_wall) call vcfftf(rwallvr,rwallvi,wr,wi,
     &              nz,nx/2,nxn*mby,1,prezn)
            else
               call vcffts(rwallur,rwallui,wr,nzn,nx/2,
     &              nxn*mby,1,presn)
               call vcftaf(rwallwr,rwallwi,wr, wi,nz/2-1,nx/2,
     &              nxn*mby,1,prean)
               if (v_wall) call vcffts(rwallvr,rwallvi,wr,
     &              nzn,nx/2,nxn*mby,1,presn)
            end if

         end if
c
c     Distribute the boundary planes among the processes
c
         call mpi_scatter(rwallur,nsr,mpi_double_precision,
     &        lwallur,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_scatter(rwallui,nsr,mpi_double_precision,
     &        lwallui,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_scatter(rwallwr,nsr,mpi_double_precision,
     &        lwallwr,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         call mpi_scatter(rwallwi,nsr,mpi_double_precision,
     &        lwallwi,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         if (v_wall) then
            call mpi_scatter(rwallvr,nsr,mpi_double_precision,
     &        lwallvr,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
            call mpi_scatter(rwallvi,nsr,mpi_double_precision,
     &        lwallvi,nsr,mpi_double_precision,0,mpi_comm_world,ierror)
         end if

      end if ! of the "nproc.gt.1"-loop

#endif
c
c     Normalization
c
      norm=real(nx*nz)

      halfproc=int(0.50001*nproc)
      zen = memnz
      if (nproc.eq.1) then
         halfproc=1
         zen = nz/2
      end if

      if (my_node.lt.halfproc) then

         do z=1,zen
            do x=1,nx/2
               lwallur(x,z)=lwallur(x,z)/norm
               lwallui(x,z)=lwallui(x,z)/norm
               lwallwr(x,z)=lwallwr(x,z)/norm
               lwallwi(x,z)=lwallwi(x,z)/norm
            end do
         end do
c
c     If bump condition for v is also prescribed
c
         if (v_wall) then
            do z=1,zen
               do x=1,nx/2
                  lwallvr(x,z)=lwallvr(x,z)/norm
                  lwallvi(x,z)=lwallvi(x,z)/norm
               end do
            end do
         end if

      end if
c
c     No z symmetry (nzc=nz): Processes holding 2nd z-half
c
      if (nfzsym.eq.0) then
         zst = 1
         if (nproc.eq.1) then
            halfproc=0
            zst = nz/2+1
         end if

         if (my_node.ge.halfproc) then

            do z=zst,memnz
               do x=1,nx/2
                  lwallur(x,z)=lwallur(x,z)/norm
                  lwallui(x,z)=lwallui(x,z)/norm
                  lwallwr(x,z)=lwallwr(x,z)/norm
                  lwallwi(x,z)=lwallwi(x,z)/norm
               end do
            end do
c
c     If bump condition for v is also prescribed
c
            if (v_wall) then
               do z=zst,memnz
                  do x=1,nx/2
                     lwallvr(x,z)=lwallvr(x,z)/norm
                     lwallvi(x,z)=lwallvi(x,z)/norm
                  end do
               end do
            end if

         end if

      end if
c
c     Obtain Dv and eta
c
      do z=1,memnz
         zb=my_node*memnz+z
c
c     Re{Dv},Im{Dv},Re{eta} and Im{eta}
c
         do x=1,nx/2
            wallur(x,z)=alfa(x)*lwallui(x,z)+beta(zb)*lwallwi(x,z)
            wallui(x,z)=-alfa(x)*lwallur(x,z)-beta(zb)*lwallwr(x,z)
            wallwr(x,z)=alfa(x)*lwallwi(x,z)-beta(zb)*lwallui(x,z)
            wallwi(x,z)=-alfa(x)*lwallwr(x,z)+beta(zb)*lwallur(x,z)
         end do
      end do
c
c     Re{v},Im{v}
c
      if (v_wall) then
         do z=1,memnz
            do x=1,nx/2
               wlvr(x,z)=lwallvr(x,z)
               wlvi(x,z)=lwallvi(x,z)
            end do
         end do
      end if
c
c     In (1,1), i.e. (alfa,beta)=(0,0) the values of u and w
c     and not Dv and eta are put back (not necessary for v)
c
      if (my_node .eq. 0) then
         wallur(1,1)=lwallur(1,1)
         wallwr(1,1)=lwallwr(1,1)
      end if


      end  subroutine updatewallrough
