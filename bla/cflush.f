c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cflush(iunit)
c
c     Compiler-dependent flush of unit iunit
c
      implicit none
      integer iunit

#ifdef IBM
c
c     IBM xlf version
c
      call flush_(iunit)
      return
#endif
c
c     Sun hp sgi dec linux version
c
      call flush(iunit)

      end subroutine cflush
