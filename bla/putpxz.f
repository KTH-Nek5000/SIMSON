c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
#ifndef ALLTOALL
      subroutine putpxz(boxr,boxi,yb,i,ur,ui,
     &     realg1,realg2,my_node)
c
c     Put an xz box into ur,ui with MPI communication
c
      implicit none
      include 'par.f'
      include 'mpif.h'

      integer yb,i
      real boxr(nxp/2+1,nzd),boxi(nxp/2+1,nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer my_node,realg1,realg2

      integer x,z
      integer ypt,ypg,npto,npget,ii,yyp
      integer ierror
      integer status1(mpi_status_size),status2(mpi_status_size)
      integer status3(mpi_status_size),status4(mpi_status_size)
      integer req1,req2
      integer itag,tagmult
c
c     This routine shall only be called if more than one
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(38987)
c
c     Copy to coarse grid
c
      do z=nz/2+1,nz
         do x=1,nx/2
            boxr(x,z)=boxr(x,z+nzp-nz)
            boxi(x,z)=boxi(x,z+nzp-nz)
         end do
      end do
c
c     Base index (yyp=yb for node 0)
c
      yyp=yb-my_node
c
c     Do the global communication
c
      tagmult = (nyp/nproc+10)*2
      do ii=1,nproc-1

         npto=my_node-ii
         npget=my_node+ii
         if (npto.lt.0) npto=nproc+npto
         if (npget.ge.nproc) npget=npget-nproc

         ypt=yyp+my_node
         ypg=yyp+npget

         if (ypt.le.nyp) then
c
c     Send x/z-plane y=ypt to npto
c
            itag = ypt*tagmult+i
            call mpi_isend(boxr(1,1+npto*memnz),1,realg2,
     &           npto,itag,mpi_comm_world,req1,ierror)
            itag = ypt*tagmult+i+1
            call mpi_isend(boxi(1,1+npto*memnz),1,realg2,
     &           npto,itag,mpi_comm_world,req2,ierror)
         end if

         if (ypg.le.nyp) then
c
c     Receive x/z-plane from npget and put on z=npget
c
            itag = ypg*tagmult+i
            call mpi_recv(ur(1,ypg,1,i),1,realg1,npget,
     &           itag,mpi_comm_world,status1,ierror)
            itag = ypg*tagmult+i+1
            call mpi_recv(ui(1,ypg,1,i),1,realg1,npget,
     &           itag,mpi_comm_world,status2,ierror)
         end if

         if (ypt.le.nyp) then
c
c     Wait until planes are delivered
c
            call mpi_wait(req1,status3,ierror)
            call mpi_wait(req2,status4,ierror)
         end if

      end do
c      call mpi_barrier(mpi_comm_world,ierror)
      if (yb.le.nyp) then
c
c     Copy received local boxr,boxi onto ur,ui
c
         do z=1,memnz
            do x=1,nx/2
               ur(x,yb,z,i)=boxr(x,z+my_node*memnz)
               ui(x,yb,z,i)=boxi(x,z+my_node*memnz)
            end do
         end do
      end if

      end subroutine putpxz

#else

      subroutine putpxz(boxr,boxi,yb,i,ur,ui,
     &     realg1,realg2,my_node)
c
c     Put an xz box into ur,ui with MPI communication
c
      implicit none
      include 'par.f'
      include 'mpif.h'

      integer yb,i
      real boxr(nxp/2+1,nzd),boxi(nxp/2+1,nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer my_node,realg1,realg2

      integer x,z,yyp,ierror
c
c     This routine shall only be called if more than one
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(38987)
c
c     Copy to coarse grid
c
      do z=nz/2+1,nz
         do x=1,nx/2
            boxr(x,z)=boxr(x,z+nzp-nz)
            boxi(x,z)=boxi(x,z+nzp-nz)
         end do
      end do
c
c     Base index (yyp=yb for node 0)
c
      yyp=yb-my_node
c
c     Do the global communication
c
      call mpi_alltoall(boxr(1,1),1,
     &     realg2,ur(1,yyp,1,i),1,realg1,
     &     mpi_comm_world,ierror)
      call mpi_alltoall(boxi(1,1),1,
     &     realg2,ui(1,yyp,1,i),1,realg1,
     &     mpi_comm_world,ierror)

      end subroutine putpxz
#endif
