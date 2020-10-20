c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
c
c -----------------------------------------------------------------------
c  Purpose:    Solves a pentadiagonal LSE by LU-decomposition.
c
c  Usage:      call  pntdiag (n, m, ld, ud2, ud1, hd, od1, od2, rs, x,
c                             ifehl, irep)
c
c  Arguments:  n     : LGS-dimensions
c              m     : vector length
c              ld    : leading dimension
c              ud2   : subsubdiagonal
c              ud1   : subdiagonal
c              hd    : diagonal
c              od1   : superdiagonal
c              od2   : supersuperdiagonal
c              rs    : right-hand side
c              x     : solution
c              ifehl : singularity indicator
c              irep  : repetition flag
c
c  Calls:      pntdiag_lu, pntdiag_e
c
c  Comments:   Overwrites input matrices with LU-decomposition.
c -----------------------------------------------------------------------
c
      subroutine  pntdiag (n, m, ld, ud2, ud1, hd, od1, od2, rs, x,
     &                     ifehl, irep)

      implicit none

      integer n,m,ld
      integer ifehl,irep

      real  ud1(n), ud2(n), hd(n), od1(n), od2(n)
      real  rs(ld,n), x(ld,n)

      ifehl = -1
      if (n .lt. 4) then
         ifehl = 999
         return
      end if

      if (irep.eq.0) then
         call  pntdiag_lu (n, ud2, ud1, hd, od1, od2, ifehl)
         if (ifehl .ne. -1) return
      end if

      call  pntdiag_e (n, ld, m, ud2, ud1, hd, od1, od2,
     &     rs, x)

      end subroutine pntdiag



c-----------------------------------------------------------------------
c  Purpose:    Computes LU-decomposition of a pentadiagonal LSE.
c
c  Usage:      call  pntdiag_lu (n, ud2, ud1, hd, od1, od2, ifehl)
c
c  Arguments:  n     : LGS-dimensions
c              ud2   : subsubdiagonal
c              ud1   : subdiagonal
c              hd    : diagonal
c              od1   : superdiagonal
c              od2   : supersuperdiagonal
c              ifehl : singularity indicator
c
c  Calls:      pntdiag_lu
c
c  Comments:   Overwrites input matrices with LU-decomposition.
c-----------------------------------------------------------------------

      subroutine  pntdiag_lu (n, ud2, ud1, hd, od1, od2, ifehl)

      implicit none

      real  ud1(n), ud2(n), hd(n), od1(n), od2(n)
      real hilf
      integer n,i,ifehl

      ud2(1)   = 0.0
      ud2(2)   = 0.0
      ud1(1)   = 0.0
      od1(n)   = 0.0
      od2(n-1) = 0.0
      od2(n)   = 0.0
      if (hd(1) .eq. 0.0) then
         ifehl = 1
         return
      end if

      hilf   = 1. / hd(1)
      od1(1) = od1(1) * hilf
      od2(1) = od2(1) * hilf
      hd(2)  = hd(2) - ud1(2) * od1(1)
      hilf   = 1. / hd(2)
      od1(2) = (od1(2) - ud1(2) * od2(1)) * hilf
      od2(2) = od2(2) * hilf

      do i = 3 , n-2 , 1
         ud1(i) = ud1(i) - ud2(i) * od1(i-2)
         hd(i)  = hd(i) - ud2(i) * od2(i-2) - ud1(i) * od1(i-1)
         if (hd(i) .eq. 0.0) then
            ifehl = i
            return
         end if
         hilf = 1. / hd(i)
         od1(i) = (od1(i) - ud1(i) * od2(i-1)) * hilf
         od2(i) = od2(i) * hilf
      end do

      ud1(n-1) = ud1(n-1) - ud2(n-1) * od1(n-3)
      hd(n-1)  = hd(n-1) - ud2(n-1) * od2(n-3) - ud1(n-1) * od1(n-2)
      if (hd(n-1) .eq. 0.0) then
         ifehl = n - 1
         return
      end if

      hilf = 1. / hd(n-1)
      od1(n-1) = (od1(n-1) - ud1(n-1) * od2(n-2)) * hilf
      ud1(n)   = ud1(n) - ud2(n) * od1(n-2)
      hd(n)    = hd(n) -ud2(n) * od2(n-2) - ud1(n) * od1(n-1)
      if (hd(n) .eq. 0.0) then
         ifehl = n
         return
      end if

      end subroutine pntdiag_lu






c-----------------------------------------------------------------------
c  Purpose:    Solves a pentadiagonal LSE by LU-decomposition.
c
c  Usage:      call  pntdiag_e
c                    (n, ld, m, ud2, ud1, hd, od1, od2, rs, x,
c                     ifehl, irep)
c
c  Arguments:  n     : LGS-dimensions
c              ld    : leading dimension
c              m     : vector length
c              ud2   : subsubdiagonal
c              ud1   : subdiagonal
c              hd    : diagonal
c              od1   : superdiagonal
c              od2   : supersuperdiagonal
c              rs    : right-hand side
c              x     : solution
c              ifehl : singularity indicator
c              irep  : repetition flag
c
c  Comments:   Needs LU-decomposed matrices as input.
c
c-----------------------------------------------------------------------

      subroutine  pntdiag_e (n, ld, m,  ud2, ud1, hd, od1, od2,
     &                       rs, x)

      implicit none

      real  ud1(n), ud2(n), hd(n), od1(n), od2(n)
      real  rs(ld,n), x(ld,n),hilf
      integer n,ld,m,i,j

      hilf = 1./hd(1)
      do i = 1 , m
         rs(i,1) = rs(i,1)*hilf
      end do
      hilf = 1./hd(2)
      do i = 1 , m
         rs(i,2) = (rs(i,2)-ud1(2)*rs(i,1))*hilf
      end do

      do j = 3 , n
         hilf = 1./hd(j)
         do i = 1 , m
            rs(i,j) = (rs(i,j)-
     &           ud2(j)*rs(i,j-2)-
     &           ud1(j)*rs(i,j-1))*hilf
         end do
      end do

      do i = 1 , m
         x(i,n) = rs(i,n)
         x(i,n-1) = rs(i,n-1)-
     &        od1(n-1)*x(i,n)
      end do

      do j = n-2 , 1 , -1
         do i = 1 , m
            x(i,j) = rs(i,j)-
     &           od1(j)*x(i,j+1)-
     &           od2(j)*x(i,j+2)
         end do
      end do

      end subroutine pntdiag_e
