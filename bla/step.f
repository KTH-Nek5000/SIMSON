c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      real function step(x)
c
c     Smooth step function:
c     x<=0 : step(x) = 0
c     x>=1 : step(x) = 1
c     Non-continuous derivatives at x=0.02 and x=0.98
c
      implicit none

      real x

      if (x.le.0.02) then
         step = 0.0
      else
         if (x.le.0.98) then
            step = 1./( 1. + exp(1./(x - 1.) + 1./x) )
         else
            step = 1.
         end if
      end if

      end function step
