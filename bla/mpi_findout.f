c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine mpi_findout(communicator,all_names,my_name,
     &     all_procs,my_procs,all_threads,my_threads)

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer communicator, all_nodes, ierror, tot_procs, tot_threads
      integer my_procs, my_threads,my_node

#ifdef OPENMP
      integer omp_get_num_procs, omp_get_num_threads
#endif

      character(len=mpi_max_processor_name) my_name
      character(len=mpi_max_processor_name) all_names(0:*)

      integer all_threads(0:*),all_procs(0:*)
      integer idummy,i
c
c     Get to know about the local environment...
c
      call mpi_comm_rank(communicator, my_node, ierror)
      call mpi_comm_size(communicator, all_nodes, ierror)

#ifdef OPENMP
!$omp parallel
!$omp single
      my_threads = omp_get_num_threads()
      my_procs = omp_get_num_procs()
!$omp end single
!$omp end parallel
#else
      my_threads = 1
      my_procs = 1
#endif

      my_name = ' '
      call mpi_get_processor_name(my_name,idummy,ierror)

      if (all_nodes.gt.1) then
         call mpi_allgather(my_procs,1,mpi_integer4,all_procs,
     &        1,mpi_integer4,communicator,ierror)
         call mpi_allgather(my_threads,1,mpi_integer4,all_threads,1,
     &        mpi_integer4,communicator,ierror)
         call mpi_allgather(my_name,mpi_max_processor_name,
     &        mpi_character,
     &        all_names,mpi_max_processor_name,mpi_character,
     &        communicator,ierror)
      else
         all_procs(0) = my_procs
         all_threads(0) = my_threads
         all_names(0) = my_name
      end if
c
c     Correct c-type strings
c
      do i=0,all_nodes-1
         call remove_0(all_names(i),mpi_max_processor_name)
      end do



      if (my_node.eq.0) then
c
c     Check for duplicates in the nodes list
c
c         do i=0, all_nodes-2
c            do j=i+1, all_nodes-1
c               if (all_names(i).eq.all_names(j)) then
c                  write(*,*) 'duplicate node name: ',trim(all_names(i))
c               end if
c            end do
c         end do

c
c     Give welcome message
c
         tot_procs = 0
         tot_threads = 0
         do i=0, all_nodes-1
            tot_procs = tot_procs + all_procs(i)
            tot_threads = tot_threads + all_threads(i)
         end do

         write(*,'(a,i1,a1,i1,a)')
     &        ' starting MPI-',mpi_version,'.',mpi_subversion,' code.'
         write(*,'(a,i4,a,i4,a,i4,a)')
     &        ' using ',all_nodes,' nodes with total ',tot_procs,
     &        ' processors and ',tot_threads,' threads.'
         do i=0,all_nodes-1
            write(*,'(a,i4,a,i4,a,i4,a,a)')
     &           ' node ',i,': procs=',all_procs(i),
     &           ' threads=',all_threads(i),
     &           ' name=',trim(all_names(i))
         end do

         open(file='nodes.out',unit=11,status='unknown')
         write(11,'(a,i1,a1,i1,a)')
     &        ' starting MPI-',mpi_version,'.',mpi_subversion,' code.'
         write(11,'(a,i4,a,i4,a,i4,a)')
     &        ' using ',all_nodes,' nodes with total ',tot_procs,
     &        ' processors and ',tot_threads,' threads.'
         do i=0,all_nodes-1
            write(11,'(a,i4,a,i4,a,i4,a,a)')
     &           ' node ',i,': procs=',all_procs(i),
     &           ' threads=',all_threads(i),
     &           ' name=',trim(all_names(i))
         end do
         close(11)
      end if

      end subroutine mpi_findout


      subroutine remove_0(s,len)

      implicit none

      integer len,i
      character s(len)

      do i=1,len
         if (ichar(s(i)).eq.0) s(i)=' '
      end do

      end subroutine remove_0
