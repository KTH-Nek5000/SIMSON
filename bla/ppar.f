c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine ppar(my_node,all_nodes,ompnthreads,ompnproc)
c
c     Prints out all static parameters
c     and performs basic sanity check
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
      integer ierror
#endif
      character(len=80) current_directory
      character(len=20) hostname,fftstring
      logical end
      integer ntmp,stop_now,my_node
      integer all_nodes,ompnthreads,ompnproc
c
c     Only processor zero tests the parameters
c
      if (my_node.gt.0) goto 9999
c
c     Get the local host name
c
      hostname = 'not known'
      call getenv('HOST',hostname)

      current_directory = 'not known'
      call getenv('PWD',current_directory)

      write(ios,*) 'Name of host machine : ',trim(hostname)
      write(ios,*) 'Current directory    : ',trim(current_directory)
c
c     Writing out preprocessor information
c
#ifdef LINUX_INTEL
      write(ios,*) 'Compiled on          : Linux (Intel Fortran)'
#endif
#ifdef LINUX_PGI
      write(ios,*) 'Compiled on          : Linux (PGI)'
#endif
#ifdef IBM
      write(ios,*) 'Compiled on          : IBM (xlf)'
#endif

#ifdef OPENMP
      write(ios,*) 'Compiled             : with OpenMP'
#else
      write(ios,*) 'Compiled             : without OpenMP'
#endif
#ifdef MPI
      write(ios,*) 'Compiled             : with MPI'
#else
      write(ios,*) 'Compiled             : without MPI'
#endif

      call fft_identify(fftstring)
      write(ios,*) 'Linked FFT package   : ',fftstring

      write(ios,*) 'Compiled for         :',nproc,  ' MPI processors'
      write(ios,*) '                      ',nthread,' OMP threads'
      write(ios,*) 'Running on           :',all_nodes,' MPI processors'
      write(ios,*) '                      ',ompnthreads,' OMP threads'
      write(ios,*) '                      ',ompnproc,
     &     ' OMP processors per node'


      include 'latest_revision'
      write(ios,*)
      write(ios,*)
      write(ios,*) '>>>>>>  Parameters read from par.f <<<<<<<'
      write(ios,*) '--------------------------------------------'//
     &           '-----------------------'
      write(ios,*) '            nx,ny,nz : ',nx,ny,nz
      write(ios,*) '      nfxd,nfyd,nfzd : ',nfxd,nfyd,nfzd
      write(ios,*) '              nfzsym : ',nfzsym
      write(ios,*) '             mby,mbz : ',mby,mbz
      write(ios,*) '               nproc : ',nproc
      write(ios,*) '             nthread : ',nthread
      write(ios,*) '         nxys,nxysth : ',nxys,nxysth
      write(ios,*) '    mcorr,mser,msamp : ',mcorr,mser,msamp
      write(ios,*) '            pressure : ',pressure
      write(ios,*) '              scalar : ',scalar
      write(ios,*) '               npart : ',npart
      write(ios,*) '                osnf : ',osnf
      write(ios,*) '                mbla : ',mbla
      write(ios,*) '         nxp,nyp,nzp : ',nxp,nyp,nzp
      write(ios,*) '        nzc,nzd,nzpc : ',nzc,nzd,nzpc
      write(ios,*) '           nzat,nzst : ',nzat,nzst
      write(ios,*) '             nby,nbz : ',nby,nbz
      write(ios,*) '   memnx,memny,memnz : ',memnx,memny,memnz
      write(ios,*) '             memnxyz : ',memnxyz
      write(ios,*) '         mbox2,mbox3 : ',mbox2,mbox3
      write(ios,*) '              nppart : ',nppart
c
c     Make runtime check of parameters
c
      end = .false.
      end = (nfxd.ne.0.and.nfxd.ne.1).or.
     &      (nfyd.ne.0.and.nfyd.ne.1).or.
     &      (nfzd.ne.0.and.nfzd.ne.1).or.
     &      (nfzsym.ne.0.and.nfzsym.ne.1)
      if (end) then
         write(ios,*)'par.f: nfxd, nfyd, nfzd, nfzsym must be 1 or 0'
      end if
      if (mby.lt.1.or.mbz.lt.1.or.nproc.lt.1) then
         write(ios,*)'par.f: mby, mbz, nproc must be bigger than 0.'
         end=.true.
      end if
      if (mod(nzc,mbz).ne.0) then
         write(ios,*) 'par.f: nzc must divisible by mbz'
         end=.true.
      end if
      if (nz.eq.1.and.nfzd.eq.1) then
         write(ios,*) 'par.f: z-dealiasing for 2-d grid is ignored.'
      end if
      if (nz.eq.1.and.nfzsym.eq.1) then
         write(ios,*) 'par.f: Do not use z-symmetry for 2-d grid'
         write(ios,*) 'par.f: Set nfzsym=0 in par.f and recompile'
         end=.true.
      end if
      if (nz.gt.1) then
         call cnxtmp(ntmp,0.,0.,nz,nfzd+nfzsym)
         if (ntmp.ne.nz) then
            if (nfzsym.eq.1.and.nfzd.eq.1) then
               write(ios,*) 'par.f: As nfzsym=1 and nfzd=1'
               write(ios,*)
     &              'nz must divisible by 8 and factorable by 2,3 and 5'
            elseif (nfzsym.eq.1) then
               write(ios,*) 'par.f: As nfzsym=1'
               write(ios,*)
     &              'nz must divisible by 4 and factorable by 2,3 and 5'
            elseif (nfzd.eq.1) then
               write(ios,*) 'par.f: As nfzd=1'
               write(ios,*)
     &              'nz must divisible by 4 and factorable by 2,3 and 5'
            else
               write(ios,*)
     &              'nz must divisible by 2 and factorable by 2,3 and 5'
            end if
            write(ios,*) 'Closest possible nz is',ntmp
            end=.true.
         end if
      end if
      call cnxtmp(ntmp,0.,0.,nx,nfxd)
      if (ntmp.ne.nx) then
         write(ios,*) 'par.f: nx must be even and'//
     &        ' factorable by 2,3 and 5'
         write(ios,*) 'Closest possible nx is',ntmp
         end=.true.
      end if
      call cnxtmp(ntmp,0.,0.,ny-1,nfyd)
      if (ntmp.ne.ny-1) then
         write(ios,*) 'par.f: ny-1 must be even and '
         write(ios,*) 'factorable by 2,3 and 5'
         write(ios,*) 'Closest possible ny-1 is',ntmp
         end=.true.
      end if
      if (mod(nx,2+2*nfxd).ne.0.or.mod(ny-1,2+2*nfyd).ne.0) then
         write(ios,*)
     &        'par.f: nx or ny-1 must be divisible '
         write(ios,*) 'by 2+2*nfxd or 2+2*nfyd.'
         end=.true.
      end if

      if (mod(nz,nproc).ne.0) then
         write(ios,*) 'par.f: nz must be a multiple of nproc.'
         end=.true.
      end if
c
c     The present implementation does not allow for
c     boxes with size greater than 1. This can however be fixed for
c     runs on vector machines.
c
      if (mby.ne.1) then
         write(ios,*) 'par.f: mby must be 1.'
         end=.true.
      end if
      if (mbz.ne.1) then
         write(ios,*) 'par.f: mbz must be 1.'
         end=.true.
      end if

      if (scalar.lt.0) then
         write(ios,*) 'par.f: scalar must be positive.'
         end=.true.
      end if
      if (pressure.lt.0.or.pressure.gt.1) then
         write(ios,*) 'par.f: pressure must be 1 or 0'
         end=.true.
      end if
      if (nfzsym.ne.0) then
         write(ios,*) 'par.f: Symmetry not tested yet (nfzsym.ne.0)'
         end=.true.
      end if

      if (nproc.ne.all_nodes) then
         write(ios,*) 'par.f: Mismatch nproc in par.f and MPI nodes:'
         write(ios,*) 'Compiled for : ',nproc,' MPI processors'
         write(ios,*) 'Running on   : ',all_nodes,' MPI processors'
         end = .true.
      end if

      if (nthread.lt.ompnthreads) then
         write(ios,*) 'par.f: Increase nthreads:'
         write(ios,*) 'Compiler for max.: ',nthread,' OMP threads'
         write(ios,*) 'Running on       : ',ompnthreads,' OMP threads'
         end = .true.
      end if


      stop_now = 0
      if (end) then
         stop_now = 1
      end if

 9999 continue

#ifdef MPI
c
c     Communicate break criterion
c
      if (nproc.gt.1) then
         call mpi_bcast(stop_now,1,mpi_integer4,0,mpi_comm_world,
     &        ierror)
      end if
#endif
      if (stop_now.eq.1) then
         call stopnow(3)
      end if

      end subroutine ppar
