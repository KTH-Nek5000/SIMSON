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
      implicit none

      include 'par.f'

      character*80 namnin
      integer p,pseed,j
      real part(10,nppart)
      real pos(3,nppart)
      real,external :: pran2
      integer kode,my_node,npart_in,pp,npartt
      real pint(3,npart)
      real t,t_in
      real re,retau,stp

      real xl,zl,fstart
c
c     Initialise to zero
c
      part = 0.
      pos  = 0.
      pint = 0.
c
c     Check file existence
c
      open(unit=11,file=trim(namnin)//'.part',status='old',
     &     form='unformatted',iostat=kode)
      if (kode.eq.0) then
         if (my_node.eq.0) then
            write(ios,*) 'Initialise particles from file ',
     &           trim(namnin),'.part'
         end if

         read(11) npart_in,t_in
         if (npart_in.ne.npart) then
            if (my_node.eq.0) then
               write(ios,*) 'npart and npart_in do not agree'
            end if
            call stopnow(556453)
         end if
         if (t.ne.t_in) then
            if (my_node.eq.0) then
               write(ios,*) 't and t_in do not agree'
            end if
            call stopnow(556454)
         end if
c
c     Read file
c     (Might be a bit inefficient)
c
         read(11) pint
         do p=1,nppart
            do j=1,3
               pos(j,p) = pint(j,p+nppart*my_node)
            end do
         end do
         read(11) pint
         do p=1,nppart
            do j=1,3
               part(j,p) = pint(j,p+nppart*my_node)
            end do
         end do
         read(11) (pint(1,p),p=1,npart)
         do p=1,nppart
            part(10,p) = pint(1,p+nppart*my_node)
         end do
         close(11)

      else
         if (my_node.eq.0) then
            write(ios,*) 'Initialise particles evenly'

         end if

         if (abs(mod(real(npart)**(1./3.)+0.5,1.)-0.5).gt.0.000001) then
            if (my_node.eq.0) then
               write(ios,*) 'npart in par.f must be equal to int^3'
            end if
            call stopnow(5564545)
         end if

         npartt=npart
         do p=1,nppart
            pos(1,p)=xl*((real(p-1+nppart*my_node)+0.5*real(npart)
     &           **(2./3.))/real(npart)-mod(real(p-1+nppart*my_node),
     &           real(npart)**(2./3.))/real(npartt))
            pos(2,p)=((mod(real(p-1+nppart*my_node),
     &           real(npart)**(2./3.))/real(npart)**(2./3.)
     &           -mod(real(p-1+nppart*my_node),real(npart)**(1./3.))
     &           /real(npart)**(2./3.)+0.5/real(npartt)**(1./3.))-0.5
     &           )*2.
            pos(3,p)=zl*(mod(real(p+nppart*my_node)-0.5,
     &           real(npart)**(1./3.))/real(npartt)**(1./3.)-0.5)

c
c     St+ = tau_p * (u_tau)**2 * Re
c         = tau_p * Re_tau**2/Re
c
c     NOTE: if there are lagrangian particle, they are the last
c     population to evolve
c

            retau = 180.

            pp = p+my_node*nppart

c            if (pp.le.npart/6) then
c               stp = 1.
c            else if (pp.le.2*npart/6) then
c               stp = 5.
c            else if (pp.le.3*npart/6) then
c               stp = 10.
c            else if (pp.le.4*npart/6) then
c               stp = 50.
c            else if (pp.le.5*npart/6) then
c               stp = 100.
c            else
c               stp = 0.
c            end if

            stp=50.

            part(10,p) = stp * re/(retau**2)




            part(1,p)=0.        !pran2(seed)
            part(2,p)=0.
            Part(3,p)=0.
         end do
      end if

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
