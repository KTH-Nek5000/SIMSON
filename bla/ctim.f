c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine ctim(ctime,wtime)
c
c     Get wall time and CPU time
c
      implicit none

      real wtime,ctime

      call cpu_time(ctime)
      call wall_time(wtime)

      end subroutine ctim
