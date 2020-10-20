c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine comment(iunit)
c
c     Reads from the formatted unit iunit until the next
c     non-comment line (i.e. without starting with an '#')
c
      implicit none

      integer iunit
      character line

      do
         read(iunit,*) line
         if (line.ne.'#') exit
      end do

      backspace(iunit)

      end subroutine comment
