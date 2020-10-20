c************************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine part_ini(part,pos,pseed,xl,zl,fstart,namnin,my_node,
     &                    pint,t,re)

c     this subroutine was removed because MOD(P,0) is throwing an error
c     with gfortran. Please use the code in part_ini_full.f if needed.


      end subroutine part_ini












      subroutine part_write(part,pos,my_node,namnut,t,pint)
      implicit none
      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif
      real part(10,nppart)
      real pos(3,nppart)
      real pint(3,npart,3)
      real t
      integer my_node
      character*80 namnut
      integer p,i
#ifdef MPI
      integer ierror,ip
      integer status1(mpi_status_size)
#endif


c
c     Copy the local data
c
      do p=1,nppart
         do i=1,3
            pint(i,p,1) = pos(i,p)
            pint(i,p,2) = part(i,p)
         end do
         pint(p,1,3) = part(10,p)
      end do

c
c     Collect data from processors
c
#ifdef MPI
      if (my_node.eq.0) then
         if (nproc.gt.0) then
            do ip=1,nproc-1
               call mpi_recv(pint(1,1+ip*nppart,1),nppart*3,
     &              mpi_double_precision,
     &              ip,ip+100,mpi_comm_world,status1,ierror)
               call mpi_recv(pint(1,1+ip*nppart,2),nppart*3,
     &              mpi_double_precision,
     &              ip,ip+200,mpi_comm_world,status1,ierror)
               call mpi_recv(pint(1+ip*nppart,1,3),nppart,
     &              mpi_double_precision,
     &              ip,ip+300,mpi_comm_world,status1,ierror)
            end do
         end if
      else
         call mpi_ssend(pint(1,1,1),nppart*3,
     &        mpi_double_precision,0,my_node+100,
     &        mpi_comm_world,ierror)
         call mpi_ssend(pint(1,1,2),nppart*3,
     &        mpi_double_precision,0,my_node+200,
     &        mpi_comm_world,ierror)
         call mpi_ssend(pint(1,1,3),nppart,
     &        mpi_double_precision,0,my_node+300,
     &        mpi_comm_world,ierror)
      end if
#endif
c
c     Write output
c
      if (my_node.eq.0) then
         open(unit=11,file=trim(namnut)//'.part',form='unformatted')

         write(11) npart,t
         write(11) ((pint(i,p,1),i=1,3),p=1,npart)
         write(11) ((pint(i,p,2),i=1,3),p=1,npart)
         write(11) ((pint(p,i,3),i=1,1),p=1,npart)

         close(11)
      end if


      end subroutine part_write




      real function pran2(idum)
c
c     A simple portable random number generator
c
c     Requires 32-bit integer arithmetic
c     Taken from Numerical Recipes, William Press et al.
c     gives correlation free random numbers but does not have a very large
c     dynamic range, i.e only generates 714025 different numbers
c     for other use consult the above
c     Set idum negative for initialization
c
      implicit none

      integer idum,j

      integer,parameter :: m=714025,ia=1366,ic=150889
      real,parameter :: rm=1./m

      integer,save :: ir(97),iy
      integer,save :: iff = 0

      if (idum.lt.0.or.iff.eq.0) then
c
c     Initialize
c
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
c
c     Generate random number
c
      j=1+(97*iy)/m
      iy=ir(j)
      pran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum

      end function pran2
