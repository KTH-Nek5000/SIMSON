c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      real function cubip(x,ftab,dx,n)
c
c     Interpolates in an equidistant table ftab with
c     n entries spaced at dx starting at x0=0, ending at x1=dx*(n-1)
c     the function is assumed to be converged towards the top
c     so that for input x > x1 the top element in the table is returned
c
      implicit none

      integer n
      real x,ftab(n),dx

      real a,b,c,d,p
      integer m

      if (x.gt.dx*real(n-1)) then
         cubip=ftab(n)
         return
      end if
      if (x.lt.0) then
         write(*,*) 'CUBIP : entry not in table'
         stop
      end if
      m=int(x/dx)
      m=min(max(m,1),n-3)
      p=2.*x/dx-real(2*m+1)
      a=(-3.*ftab(m)+27.*ftab(m+1)+27.*ftab(m+2)-3.*ftab(m+3))/48.
      b=(1.*ftab(m)-27.*ftab(m+1)+27.*ftab(m+2)-1.*ftab(m+3))/48.
      c=(3.*ftab(m)-3.*ftab(m+1)-3.*ftab(m+2)+3.*ftab(m+3))/48.
      d=(-1.*ftab(m)+3.*ftab(m+1)-3.*ftab(m+2)+1.*ftab(m+3))/48.
      cubip=a+b*p+c*p*p+d*p*p*p

      return

      end function cubip
