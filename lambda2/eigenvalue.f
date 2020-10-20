c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine eigenvalue(a,x)
c     gives back the real eigenvalues of a 3x3 real matrix

      implicit none
      real a(3,3)
      real x(3)

      real ca,cb,cc
      real Q,R

      complex AA,BB,cu

      cu = cmplx(0.,1.)

      ca = -(a(1,1)+a(2,2)+a(3,3))
      cb = -(a(1,2)*a(2,1)+a(1,3)*a(3,1)+a(2,3)*a(3,2)-
     &     a(1,1)*a(2,2)-a(1,1)*a(3,3)-a(2,2)*a(3,3))
      cc = -(a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+
     &     a(1,3)*a(2,1)*a(3,2)-a(1,1)*a(2,3)*a(3,2)-
     &     a(1,2)*a(2,1)*a(3,3)+a(1,1)*a(2,2)*a(3,3))


c     the solution to x^3 + A*x^2 + B*x + C = 0.

      Q = (ca**2-3.*cb)/9.
      R = (2.*ca**3-9.*ca*cb+27.*cc)/54.

      AA = -sign(1.,R)*(abs(R)+sqrt(cmplx(R**2-Q**3)))**(1./3.)

      if (AA.eq.0.) then
         BB = 0.
      else
         BB = cmplx(Q)/AA
      end if

      x(1) = (AA+BB) - ca/3.
      x(2) = -.5*(AA+BB)-ca/3.+cu*sqrt(3.)/2.*(AA-BB)
      x(3) = -.5*(AA+BB)-ca/3.-cu*sqrt(3.)/2.*(AA-BB)

      end subroutine eigenvalue
