c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wall_time(wt)
c
c     Return wall-clock time as seconds after Jan. 1, 2005.
c     Support for leap year is deactivated.
c
c     By using a 'save' statement, the wall-time after the first
c     call to the subroutine could be computed, but that is not
c     intended with the present subroutine (e.g. the history file)
c
      implicit none

      real wt
      integer val(8),i,shift,day

      integer mon(12,2)
      data mon /
     &     31,28,31,30,31,30,31,31,30,31,30,31,
     &     31,29,31,30,31,30,31,31,30,31,30,31/
c
c     Get current date and time
c     val(1) : year
c     val(2) : month
c     val(3) : day
c     val(4) : difference to GMT
c     val(5) : hour
c     val(6) : minute
c     val(7) : second
c     val(8) : 1/1000 second
c
      call date_and_time(values=val)
c
c     Determine leap year
c
c      if (mod(val(1),4).eq.0) then
c         if (mod(val(1),100).eq.0) then
c            if (mod(val(1),400).eq.0) then
c               shift=2
c            else
c               shift=1
c            end if
c         else
c            shift=2
c         end if
c      else
c         shift = 1
c      end if
c
      shift = 1
c
c     Construct day of the year
c
      day = val(3)-1
      do i=1,val(2)-1
         day=day+mon(i,shift)
      end do
c
c     And compute wall-clock time
c
      wt = (val(1)-2005)*365*86400+
     &     day*86400+val(5)*3600+val(6)*60+val(7)+val(8)/1000.

      end subroutine wall_time
