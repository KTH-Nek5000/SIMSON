c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      real function dstep(x)
c
c     Computes first derivative of smooth step function (see step.f)
c
      implicit none

      real x,t1,t5,t7,t9,t11

      if (x.le.0.98.and.x.ge.0.02) then
         t1 = x-1.
         t5 = exp(1./t1+1./x)
         t7 = (1.+t5)*(1.+t5)
         t9 = t1*t1
         t11= x**2
         dstep = 1./t7*(1./t9+1./t11)*t5
      else
         dstep = 0.
      end if

      end function dstep
