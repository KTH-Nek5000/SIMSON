c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wxys(xys,xysth,
     &     totcorr,totcorr_x,corr,corr_x,ncorr,ncorr_x,corrf,corrf_x,
     &     totcorr_ms,totcorr_ms_x,corr_ms,corr_ms_x,my_node,
     &     series,serf,nser,totseries,
     &     sumw,re,xl,zl,t,dstar,fltype,
     &     bstart,blength,rlam,spanv,namxys,wxy,
     &     corrnam,corrnam_x,corrxx,corryy,corrzz,corryy_x,corrwsave,
     &     corrwsave_x,corrx,corry,corrz,corry_x,corrxi,corryi,corrzi,
     &     corryi_x,pr,m1,corrt,corrt_x,nsamp,namser,sercount,mhd_n,b0,
     &     it,ixyss)
c
c     Sum the statistics xys and the correlations corr
c     from all the processors and writes it to a file
c
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      include 'par.f'
      integer y,i,j,x,yb,z,it, ixyss
      real xys     (nx,nyp/nproc+1,nxys)
      real xysth   (nx,nyp/nproc+1,nxysth,scalar)
      real corr(nzc+2,mcorr),corr_x(nx+2,mcorr)
      real totcorr(nzc+2),totcorr_x(nx+2)
      real totcorr_ms(4),corr_ms(4,mcorr)
      real totcorr_ms_x(4),corr_ms_x(4,mcorr)
      logical corrf,corrf_x
      integer ncorr,ncorr_x,ith

      real pr(scalar),m1(scalar)
      integer fltype
      real sumw,sumww
      real wxy(nx,nyp)
      real mhd_n,b0(3)
c
c     Two-point correlation
c
      real corrxx(mcorr),corryy(mcorr)
      real corrzz(mcorr),corryy_x(mcorr)
      real corrx(mcorr),corry(mcorr)
      real corrz(mcorr),corry_x(mcorr)
      integer corrxi(mcorr),corryi(mcorr)
      integer corrzi(mcorr),corryi_x(mcorr)
      real cd(nzc+2),cd_x(nx+2)
      integer corrt(mcorr),corrt_x(mcorr)
      real corrwsave(nzc+15),cw(nzc+2),cwx(nx+2),corrwsave_x(nx+15)
c
c     Time series
c
      real totseries(msamp,0:mser)
      integer nser,nsamp,sercount
      logical serf
      character(len=80) namser

      real re,xl,zl,t,dstar,bstart,blength,rlam,spanv

      character*80 namxys
      character*80 corrnam,corrnam_x
c
c     Time series
c
      real series(msamp,0:mser)
c
c     MPI
c
      integer my_node
#ifdef MPI
      real txys(nx,nyp/nproc+1)
      integer ierror,ip
#endif

c
c     Velocity, pressure and scalar xy-statistics
c     *******************************************
c
      if (my_node.eq.0) then
c
c     Open statistics file
c
         if (nxys .gt. 0) then
            open(unit=19,file=namxys,form='unformatted')
            rewind(19)
c
c     Write header
c
            write(19) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar)
c     Item 'A' added to identify if line is included or not when reading data
            write(19) 'A',mhd_n,b0
            write(19) nx,nyp,nzc,nfzsym
            write(19) fltype,dstar
            if (fltype.lt.0) write(19) rlam
            if (fltype.ge.6) write(19) bstart,blength,rlam,spanv
c
c     Write number of statistic quantities
c
            write(19) sumw,nxys,nxysth,scalar
         end if

      end if

c
c     Loop over the velocity statistics
c
      do i=1,nxys

c
c     Send the individual statistics to processor 0
c

#ifdef MPI
         if (my_node.gt.0) then
            call mpi_ssend(xys(1,1,i),nx*(nyp/nproc+1),
     &           mpi_double_precision,0,my_node+100,
     &           mpi_comm_world,ierror)
         end if
#endif

         if (my_node.eq.0) then
c
c     Receive individual statistics in txys and put on wxy
c
#ifdef MPI
            if (nproc.gt.1) then
               do ip=1,nproc-1
                  call mpi_recv(txys,nx*(nyp/nproc+1),
     &                 mpi_double_precision,
     &                 ip,ip+100,mpi_comm_world,mpi_status_ignore,
     &                 ierror)
                  do y=1,nyp/nproc+1
                     yb=ip+1+(y-1)*nproc
                     if (yb.le.nyp) then
                        do x=1,nx
                           wxy(x,yb)=txys(x,y)
                        end do
                     end if
                  end do
               end do

            end if
#endif
c
c     Put own statistics on wxy
c
            do y=1,nyp/nproc+1
               yb=1+(y-1)*nproc
               if (yb.le.nyp) then
                  do x=1,nx
                     wxy(x,yb)=xys(x,y,i)
                  end do
               end if
            end do
c
c     Divide accumulated statistics by sum of weights
c     and write velocity statistics
c
            do y=1,nyp
               do  x=1,nx
                  wxy(x,y)=wxy(x,y)*(1./sumw)
               end do
            end do
            write(19) wxy

         end if

      end do
c
c     Loop over scalars and scalar statistics
c
      do ith=1,scalar
         do i=1,nxysth

#ifdef MPI
            if (my_node.gt.0) then
               call mpi_ssend(xysth(1,1,i,ith),nx*(nyp/nproc+1),
     &              mpi_double_precision,0,my_node+200,
     &              mpi_comm_world,ierror)
            end if
#endif

            if (my_node.eq.0) then
c
c     Receive individual statistics in txys and put on wxy
c
#ifdef MPI
               if (nproc.gt.1) then
                  do ip=1,nproc-1
                     call mpi_recv(txys,
     &                    nx*(nyp/nproc+1),
     &                    mpi_double_precision,
     &                    ip,ip+200,mpi_comm_world,mpi_status_ignore,
     &                    ierror)
                     do y=1,nyp/nproc+1
                        yb=ip+1+(y-1)*nproc
                        if (yb.le.nyp) then
                           do x=1,nx
                              wxy(x,yb)=txys(x,y)
                           end do
                        end if
                     end do
                  end do
               end if

#endif
c
c     Put own scalar statistics on wxy
c
               do y=1,nyp/nproc+1
                  yb=1+(y-1)*nproc
                  if (yb.le.nyp) then
                     do x=1,nx
                        wxy(x,yb)=xysth(x,y,i,ith)
                     end do
                  end if
               end do
c
c     Divide accumulated statistics by sum of weights
c     and write scalar statistics
c
               do y=1,nyp
                  do  x=1,nx
                     wxy(x,y)=wxy(x,y)*(1./sumw)
                  end do
               end do
               write(19) wxy
            end if
         end do
      end do
c
c     Close statistics file
c
      if (my_node.eq.0) then
         close(unit=19)
      end if
c
c     Two-point correlations
c     **********************
c
      if (corrf) then
         if (my_node.eq.0) then
            open(unit=19,file=corrnam,form='unformatted')
            rewind(19)
c
c     Write header
c
            write(19) -re*dstar,.false.,xl/dstar,zl/dstar,t/dstar,0.,
     &           (pr(i),m1(i),i=1,scalar)
            write(19) nx,nyp,nzc,nfzsym
            write(19) fltype,dstar
            if (fltype.ge.6) write(19) bstart/dstar,blength/dstar
            write(19) sumw/dstar,scalar
            write(19) ncorr

         end if
c
c     Send the individual statistics to processor 0
c     Loop over the correlations and immediately write to disk;
c     This was done to save temporary storage on processor 0.
c
         do i=1,ncorr
c
c     if having saved a whole x plane...
c
            if (corrxi(i).eq.0) then
               sumww = sumw*nx
            else
               sumww = sumw
            end if


            if (nproc.eq.1) then
               do z=1,nzc+2
                  totcorr(z) = corr(z,i)
               end do
               do z=1,4
                  totcorr_ms(z) = corr_ms(z,i)
               end do
            else
#ifdef MPI
               call mpi_reduce(corr(1,i),totcorr,nzc+2,
     &              mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &              ierror)
               call mpi_reduce(corr_ms(1,i),totcorr_ms,4,
     &              mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &              ierror)
#endif
            end if
            if (my_node.eq.0) then
c
c     Postprocess and write the correlations
c
c
c     Do the backward transform
c
               do j=1,nzc+2
                  cd(j)=totcorr(j)
               end do
               call vrfftb(cd(1),cd(2),cw(1),cw(2),
     &              nzc,1,2,1,corrwsave)
c
c     Write coordinate info
c
               write(19) corrx(i),corry(i),
     &              corrxi(i),corryi(i),
     &              corrxx(i),corryy(i),corrt(i)
c
c     Write mean and square
c
               write(19) (totcorr_ms(z)*(1./sumww),z=1,4)
c
c     And write data
c
               write(19) (cd(z)*(1./sumww),z=1,nzc)

            end if
         end do

         if (my_node.eq.0) then
            close(19)
         end if
      endif

      if (corrf_x) then
         if (my_node.eq.0) then
            open(unit=20,file=corrnam_x,form='unformatted')
            rewind(20)
c
c     Write header
c
            write(20) -re*dstar,.false.,xl/dstar,zl/dstar,t/dstar,0.,
     &           (pr(i),m1(i),i=1,scalar)
            write(20) nx,nyp,nzc,nfzsym
            write(20) fltype,dstar
            if (fltype.ge.6) write(20) bstart/dstar,blength/dstar
            write(20) sumw/dstar,scalar
            write(20) ncorr_x

         end if
c
c     Send the individual statistics to processor 0
c     Loop over the correlations and immediately write to disk;
c     This was done to save temporary storage on processor 0.
c
         do i=1,ncorr_x
            if (corrzi(i).eq.0) then
               sumww = sumw*nz
            else
               sumww = sumw
            end if

            if (nproc.eq.1) then
               do x=1,nx+2
                  totcorr_x(x) = corr_x(x,i)
               end do
               do x=1,4
                  totcorr_ms_x(x) = corr_ms_x(x,i)
               end do
            else
#ifdef MPI
               call mpi_reduce(corr_x(1,i),totcorr_x,nx+2,
     &              mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &              ierror)
               call mpi_reduce(corr_ms_x(1,i),totcorr_ms_x,4,
     &              mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &              ierror)
#endif
            end if
            if (my_node.eq.0) then
c
c     Postprocess and write the correlations
c
c
c     Do the backward transform
c
               do j=1,nx+2
                  cd_x(j)=totcorr_x(j)
               end do
               call vrfftb(cd_x(1),cd_x(2),cwx(1),cwx(2),
     &              nx,1,2,1,corrwsave_x)
c
c     Write coordinate info
c
               write(20) corrz(i),corry_x(i),
     &              corrzi(i),corryi_x(i),
     &              corrzz(i),corryy_x(i),corrt_x(i)
c
c     Write mean and square
c
               write(20) (totcorr_ms_x(x)*(1./sumww),x=1,4)
c
c     And write data
c
               write(20) (cd_x(x)*(1./sumww),x=1,nx)

            end if
         end do

         if (my_node.eq.0) then
            close(20)
         end if
      end if
c
c     Time Series
c     ***********
c
      if (serf) then

         if((it-1).ne.(ixyss)) then
            if (sercount.gt.nsamp) then
               if (my_node.eq.0) then
                  write(ioe,*) 'sercount larger than nsamp'
                  write(ioe,*) sercount,nsamp
               end if
               call stopnow(165645)
            end if
         else
            if ((sercount-1).gt.nsamp) then
               if (my_node.eq.0) then
                  write(ioe,*) 'sercount larger than nsamp'
                  write(ioe,*) sercount,nsamp
               end if
               call stopnow(165646)
            end if
         end if
c
c     Collect all series in totseries
c
         if (nproc.eq.1) then
            do i=0,nser
               do j=1,msamp
                  totseries(j,i) = series(j,i)
               end do
            end do
         else
#ifdef MPI
            call mpi_reduce(series(1,1),totseries(1,1),msamp*nser,
     &           mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &           ierror)
            do i=1,sercount
               totseries(i,0) = series(i,0)
            end do
#endif
         end if
c
c     Write series to disk
c
         if (my_node.eq.0) then
            open(unit=19,file=namser,form='formatted',position='append')
            do i=1,sercount
               totseries(i,0) = totseries(i,0)/dstar
               write(19,'(4000e25.16e3)') (totseries(i,j),j=0,nser)
            end do
            close(19)
         end if
      end if

      end subroutine wxys
