c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getdt(dt,cflmax,rot,deta,xl,zl,prex,prez,pres,prea,
     &     ur,ui,u2r,u2i,bu1,bu2,pert,lin,wr,wi,
     &     realg1,realg2,wbci,bf3,bf3u2r,bf3u2i,my_node,
     &     bf3tempr,bf3tempi)
c
c     Computes CFL and dt for the first time step
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      real dt,cflmax,deta(nyp),xl,zl,rot
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real wr(nxp/2+1,mby,nzd,3),wi(nxp/2+1,mby,nzd,3)
      real cfl,cflp(nby)
      logical sym
      integer yb,i,myb,wbci
      real pi
      parameter (pi = 3.1415926535897932385)

      logical lin,pert
      real cflp1,cflp2
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
c
c     3D baseflow
c
      logical bf3
      integer x,z
      real bf3u2r(nxp/2+1,nzd,nyp/nproc+1,6)
      real bf3u2i(nxp/2+1,nzd,nyp/nproc+1,6)
      real bf3tempr(nxp/2+1,nzd,3)
      real bf3tempi(nxp/2+1,nzd,3)

c
c     MPI
c
      integer my_node,realg1,realg2
      integer yp,nypp
#ifdef MPI
      real cfl1
      integer ierror,ip
#endif
c
c     Calculates CFL/dt from a timestep
c
      if (nproc.eq.1) then
         nypp = nyp
      else
         nypp = nyp/nproc+1
      end if

      do yp=1,nypp

         yb = (yp-1)*nproc+my_node+1

         myb=(yb-1)/mby+1
         do i=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1,i),u2i(1,1,1,i),yb,i,1,ur,ui)
            else
#ifdef MPI
               call mpi_barrier(mpi_comm_world,ierror)
               call getpxz(u2r(1,1,1,i),u2i(1,1,1,i),yb,i,1,ur,ui,
     &              realg1,realg2,my_node)
#endif
            end if
c
c     u,v symmetric, w antisymmetric if this is a 'symmetric' case
c
            sym=i.le.2
            if (yb.le.nyp) then
               call fft2db(u2r(1,1,1,i),u2i(1,1,1,i),sym,
     &              1,prex,prez,pres,prea,wr,wi)
            end if
         end do
         if (yb.le.nyp) then

            if (pert) then
               if (.not.bf3) then
                  call boxcflbf(cflp1,bu1,bu2,yb,deta,xl,zl)
               else

                  do i=1,3
                     do z = 1,nzd
                        do x = 1,nxp/2+1
                           bf3tempr(x,z,i) = bf3u2r(x,z,yp,i)
                           bf3tempi(x,z,i) = bf3u2i(x,z,yp,i)
                        end do
                     end do
                  end do

                  call boxcfl(cflp1,bf3tempr,bf3tempi,
     &                 yb,deta,xl,zl,wbci)
               end if
            else
               cflp1=0.
            end if
            if (.not.lin) then
               call boxcfl(cflp2,u2r,u2i,yb,deta,xl,zl,wbci)
            else
               cflp2=0.
            end if
            cflp(yp)=cflp1+cflp2
         else
            cflp(yp) = 0.
         end if
      end do
c
c     Communicate CFL
c
      cfl=0.
      do i=1,nypp
         cfl=max(cfl,cflp(i))
      end do
#ifdef MPI
      if (my_node.ne.0) then
         call mpi_ssend(cfl,1,mpi_double_precision,
     &        0,1,mpi_comm_world,ierror)
      else
         if (nproc.gt.1) then
            do ip=1,nproc-1
               call mpi_recv(cfl1,1,mpi_double_precision,
     &              ip,1,mpi_comm_world,mpi_status_ignore,ierror)
               if (cfl.lt.cfl1) cfl=cfl1
            end do
         end if
      end if
#endif
      cfl=cfl*pi+2.*abs(rot)
      dt=cflmax/cfl
#ifdef MPI
      if (nproc.gt.1) then
         call mpi_bcast(cfl,1,mpi_double_precision,0,mpi_comm_world,
     &        ierror)
         call mpi_bcast(dt,1,mpi_double_precision,0,mpi_comm_world,
     &        ierror)
      end if
#endif

      end subroutine getdt
