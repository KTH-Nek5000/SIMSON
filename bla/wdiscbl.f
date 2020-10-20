c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &     bstart,blength,rlam,spanv,namnut,m,gall,
     &     boxr,boxi,urx,alfa,zs,beta,my_node,gr)
c
c     Writes m variables from ur,ui to file namnut
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      character*80 namnut
      integer m,fltype
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real re,pr(scalar),xl,zl,t,xs,dstar,gr(scalar)
      real bstart,blength,rlam,spanv,m1(scalar)
      logical gall
      real urx(nx)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real alfa(nx/2*mbz),zs,beta(nz)

      integer x,y,z,i,zb,ii
      real uw0low
c
c     MPI
c
      integer my_node,zbp,ip
#ifdef MPI
      integer ierror

      if (nproc.gt.1) call mpi_barrier(mpi_comm_world,ierror)
#endif

      if (my_node.eq.0) then
c
c     Write file header
c     Data is shifted so that xs,zs=0 for all fields
c     Any existing file will be overwritten
c
         open(unit=11,file=namnut,form='unformatted')
         rewind(11)
         if (scalar.ge.1) then
            write(11) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar)
         else
            write(11) re,.false.,xl,zl,t,0.
         end if
         write(11) nx,nyp,nzc,nfzsym
         write(11) fltype,dstar
         if (fltype.eq.-1) write(11) rlam
         if (fltype.eq.-2) write(11) rlam,spanv
         if (fltype.eq.4) write(11) bstart,blength
         if (fltype.eq.5) write(11) bstart,blength
         if (fltype.ge.6.and.fltype.le.9)
     &        write(11) bstart,blength,rlam,spanv
         if (abs(fltype).eq.20) write(11) (gr(i),i=1,scalar)
         if (xs.ne.0.or.zs.ne.0) then
            if (fltype.ne.2.and.fltype.ne.5) then
               write(ios,*) 'Shift of the box by ',xs,zs
            else
               write(ios,*) 'no shift for Couette flow'
            end if
         end if

      end if

      do i=1,m
c
c     Choose which field to write
c
         if (i.ge.4) then
c
c     This is the scalar
c
            ii = 8+pressure+3*(i-4)
         else
c
c     These are the velocities
c
            ii = i
         end if

         if (my_node.eq.0.and.nproc.gt.1) then
c
c     Note that with the present messaging algorithm, the planes
c     on processor 0 are overwritten. We just save it here to disk and
c     read it in at the end of the subroutine.
c
            write(99) ur(:,:,:,ii),ui(:,:,:,ii)
            close(99)
         end if

         do ip=0,nproc-1
            if (my_node.eq.0) then
c
c     Adjust and write content of boxr,boxi
c
               do zbp = 1,memnz
                  zb = zbp+ip*memnz
                  call getxy(boxr,boxi,zbp,ii,ur,ui)
                  if (gall) then
c
c     The data is Galilei transformed so that u0low,w0low=0
c
                     if ((i.eq.1.or.i.eq.3).and.zb.eq.1) then
                        uw0low=boxr(1,1,nyp)
                        do y=1,nyp
                           boxr(1,1,y)=boxr(1,1,y)-uw0low
                        end do
                     end if
                  end if
c
c     Shift data so that xs,zs=0 for all fields
c     Note: No shift for Couette flow!
c
                  if (fltype.ne.2.and.fltype.ne.5) then
                     call xysh(boxr,boxi,xs,zs,alfa,beta,zb)
                  end if
                  do z=zb,zb+mbz-1
                     do y=1,nyp
                        do x=1,nx/2
                           urx(2*x-1)= boxr(x,z-zb+1,y)
                           urx(2*x)  = boxi(x,z-zb+1,y)
                        end do
                        write(11) urx
                     end do
                  end do
               end do
            end if
#ifdef MPI
            if (nproc.gt.1.and.ip.ne.nproc-1) then
               if (my_node.eq.ip+1) then
c
c     Send complete field to processor 0
c
                  call mpi_ssend(ur(1,1,1,ii),memnx*memny*memnz,
     &                 mpi_double_precision,0,ip+10,
     &                 mpi_comm_world,ierror)
                  call mpi_ssend(ui(1,1,1,ii),memnx*memny*memnz,
     &                 mpi_double_precision,0,ip+20,
     &                 mpi_comm_world,ierror)
               end if
               if (my_node.eq.0.and.nproc.gt.1) then
c
c     Receive individual fields from processors >0
c
                  call mpi_recv(ur(1,1,1,ii),memnx*memny*memnz,
     &                 mpi_double_precision,ip+1,ip+10,mpi_comm_world,
     &                 mpi_status_ignore,ierror)
                  call mpi_recv(ui(1,1,1,ii),memnx*memny*memnz,
     &                 mpi_double_precision,ip+1,ip+20,mpi_comm_world,
     &                 mpi_status_ignore,ierror)
               end if
            end if
c
c     Wait until communication is finished
c     (probably no needed, due to ssend above)
c
            if (nproc.gt.1) call mpi_barrier(mpi_comm_world,ierror)
#endif

         end do

         if (my_node.eq.0.and.nproc.gt.1) then
c
c     Read in the saved part on processor 0
c     Note that this is essential otherwise node 0 would lose
c     its part of the velocity field.
c
            read(99) ur(:,:,:,ii),ui(:,:,:,ii)
            close(99)
         end if

      end do

      if (my_node.eq.0) then
c
c     close the file and flush...
c
         call cflush(11)
         close(unit=11)
c
c     If one wants a status file whether a given velocity
c     file has been written, uncomment the following:
c
c         open(unit=11,file=trim(namnut)//'.written')
c         close(unit=11)
      end if

      end subroutine wdiscbl
