c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine stopnow(i)

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer i

#ifdef MPI
      integer my_node,ierror
#endif

#ifdef MPI
      call mpi_comm_rank(mpi_comm_world,my_node,ierror)
      write(*,*) '*** STOP *** at location (node ',my_node,'):',i
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_finalize(ierror)
#else
      write(*,*) '*** STOP *** at location:',i
#endif

      stop

      end subroutine stopnow
