c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
#ifndef ALLTOALL
      subroutine getpxz(boxr,boxi,yb,i,ipad,ur,ui,
     &     realg1,realg2,my_node)
c
c     Get an xz box from ur,ui with MPI communication
c     Pad for dealiasing if ipad = 1
c
c     This is (together with putpxz) the most important and time-consuming
c     subroutine when using MPI.
c
c     Here we use MPI_ISEND
c
c     B: Bufferend, returns immediately
c     S: Synchronous (no buffer used)
c     R: Ready (recv needs to be posted)
c     I: Non-blocking, returns immediately
c
c     The different sends:
c
c     MPI_Send will not return until you can use the
c     send buffer. It may or may not block (it is allowed to buffer,
c     either on the sender or receiver side, or to wait for the matching
c     receive).
c
c     MPI_Bsend: Blocking, Buffered. May buffer; returns immediately
c     and you can use the send buffer. A late add-on to the MPI
c     specification. Should be used only when absolutely necessary.
c
c     MPI_Ssend (blocking, synchronous) will not return until
c     matching receive posted
c
c     MPI_Rsend: May be used ONLY if matching receive already posted.
c     User responsible for writing a correct program.
c
c     MPI_Irsend: Performs a nonblocking ready mode send operation.
c
c     MPI_Isend: Nonblocking send. But not necessarily asynchronous
c     (unless MPI_Issend is used).
c     You can NOT reuse the send buffer until either a successful,
c     wait/test or you KNOW that the message has been received
c     (see MPI_Request_free). Note also that while the I refers to
c     immediate, there is no performance requirement on MPI_Isend.
c     An immediate send must return to the user without requiring a
c     matching receive at the destination. An implementation is free to
c     send the data to the destination before returning, as long as the
c     send call does not block waiting for a matching receive.
c     Different strategies of when to send the data offer different
c     performance advantages and disadvantages that will depend on
c     the application.
c
c     MPI_Ibsend: buffered nonblocking
c
c     MPI_Issend: Synchronous nonblocking. Note that a Wait/Test
c     will complete only when the matching receive is posted.
c
c     MPI_Irsend: As with MPI_Rsend, but nonblocking.
c
c
c     The different recvs:
c
c     MPI_Recv: Performs a blocking receive operation. Together with
c     synchrnous send it forms a completely synchronous communication.
c
c     MPI_Irecv: Performs a nonblocking receive operation. This
c     subroutine starts a nonblocking receive and returns a handle
c     to a request object. You can later use the request to query
c     the status of the communication or wait for it to complete.
c     A nonblocking receive call means the system may start writing
c     data into the receive buffer. Once the nonblocking receive
c     operation is called, do not access any part of the receive
c     buffer until the receive is complete.
c

      implicit none

      include 'par.f'
      include 'mpif.h'

      integer yb,i,ipad
      real boxr((nxp/2+1),nzd),boxi((nxp/2+1),nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer realg1,realg2,my_node

      integer x,z,nzpad
      integer yyp,ii,npto,npget,ypt,ypg
      integer ierror
      integer status1(mpi_status_size),status2(mpi_status_size)
      integer status3(mpi_status_size),status4(mpi_status_size)
      integer req1,req2
      integer itag,tagmult
c
c     This routine shall only be called if more than one
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(3898)
c
c     Base index (yyp=yb for node 0)
c
      yyp = yb-my_node
c
c     Do the global communication
c
      tagmult = (nyp/nproc+10)*2
      do ii=1,nproc-1

         npto=my_node-ii
         npget=my_node+ii
         if (npto.lt.0) npto=nproc+npto
         if (npget.ge.nproc) npget=npget-nproc

         ypt=yyp+npto
         ypg=yyp+my_node

         if (ypt.le.nyp) then
c
c     Send x/z-plane y=ypt to npto
c
            itag = ypt*tagmult+i
            call mpi_isend(ur(1,ypt,1,i),1,realg1,npto,itag,
     &           mpi_comm_world,req1,ierror)
            itag = ypt*tagmult+i+1
            call mpi_isend(ui(1,ypt,1,i),1,realg1,npto,itag,
     &           mpi_comm_world,req2,ierror)
         end if
         if (ypg.le.nyp) then
c
c     Receive x/z-plane from npget and put on z=npget
c
            itag = ypg*tagmult+i
            call mpi_recv(boxr(1,1+npget*memnz),1,realg2,
     &           npget,itag,mpi_comm_world,status1,ierror)
            itag = ypg*tagmult+i+1
            call mpi_recv(boxi(1,1+npget*memnz),1,realg2,
     &           npget,itag,mpi_comm_world,status2,ierror)
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
c     Copy to boxr,boxi
c
         do z=1,memnz
            do x=1,nx/2
               boxr(x,z+my_node*memnz)=ur(x,yb,z,i)
               boxi(x,z+my_node*memnz)=ui(x,yb,z,i)
            end do
         end do
         if (nfzsym.eq.0) then
            nzpad=(nzp-nz)*ipad
            do z=nz/2+1,nz
               do x=1,nx/2
                  boxr(x,z+nzpad)=boxr(x,z)
                  boxi(x,z+nzpad)=boxi(x,z)
               end do
            end do
         end if
c
c     Pad with zeros
c
         if (ipad.eq.1) then
            do z=(nz+1)/2+1,min(nzpc,nzp+1-nz/2)
               do x=1,(nxp/2+1)
                  boxr(x,z)=0.0
                  boxi(x,z)=0.0
               end do
            end do
            do z=1,nzpc
               do x=nx/2+1,nxp/2+1
                  boxr(x,z)=0.0
                  boxi(x,z)=0.0
               end do
            end do
         else
c
c     Oddball zeroing only
c
            if (mod(nz,2).eq.0) then
               do x=1,(nxp/2+1)
                  boxr(x,nz/2+1)=0.0
                  boxi(x,nz/2+1)=0.0
               end do
            end if
            do z=1,nzpc
               boxr(nx/2+1,z)=0.0
            end do
         end if
      end if

      end subroutine getpxz

#else

      subroutine getpxz(boxr,boxi,yb,i,ipad,ur,ui,
     &     realg1,realg2,my_node)
c
c     Get an xz box from ur,ui with MPI communication
c     Pad for dealiasing if ipad = 1
c
      implicit none

      include 'par.f'
      include 'mpif.h'

      integer yb,i,ipad
      real boxr((nxp/2+1),nzd),boxi((nxp/2+1),nzd)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer realg1,realg2,my_node

      integer x,z,nzpad,yyp,ierror
c
c     This routine shall only be called if more than one
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(3898)
c
c     Base index (yyp=yb for node 0)
c
      yyp = yb-my_node
c
c     Do the global communication
c
      call mpi_alltoall(ur(1,yyp,1,i),1,
     &     realg1,boxr(1,1),1,realg2,
     &     mpi_comm_world,ierror)
      call mpi_alltoall(ui(1,yyp,1,i),1,
     &     realg1,boxi(1,1),1,realg2,
     &     mpi_comm_world,ierror)


      if (yb.le.nyp) then
c
c     Copy to boxr,boxi
c
         if (nfzsym.eq.0) then
            nzpad=(nzp-nz)*ipad
            if (nzpad.ne.0) then
               do z=nz/2+1,nz
                  do x=1,nx/2
                     boxr(x,z+nzpad)=boxr(x,z)
                     boxi(x,z+nzpad)=boxi(x,z)
                  end do
               end do
            end if
         end if
c
c     Pad with zeros
c
         if (ipad.eq.1) then
            do z=(nz+1)/2+1,min(nzpc,nzp+1-nz/2)
               do x=1,(nxp/2+1)
                  boxr(x,z)=0.0
                  boxi(x,z)=0.0
               end do
            end do
            do z=1,nzpc
               do x=nx/2+1,nxp/2+1
                  boxr(x,z)=0.0
                  boxi(x,z)=0.0
               end do
            end do
         else
c
c     Oddball zeroing only
c
            if (mod(nz,2).eq.0) then
               do x=1,(nxp/2+1)
                  boxr(x,nz/2+1)=0.0
                  boxi(x,nz/2+1)=0.0
               end do
            end if
            do z=1,nzpc
               boxr(nx/2+1,z)=0.0
            end do
         end if
      end if

      end subroutine getpxz
#endif
