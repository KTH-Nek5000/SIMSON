c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rkcoeff(anrk,bnrk,cnrk)

      implicit none
c
c     Set the Runge-Kutta time stepping coefficients
c
      real anrk(4,2),bnrk(4,2),cnrk(4,2)

      anrk = 0.
      bnrk = 0.
      cnrk = 0.
c
c     RK3 (3 stage)
c     see Wray, NASA Technical Report AA214a 4, 1987
c
      anrk(1,1)=8./15.
      anrk(2,1)=5./12.
      anrk(3,1)=3./4.

      bnrk(1,1)=0.
      bnrk(2,1)=-17./60.
      bnrk(3,1)=-5./12.

      cnrk(1,1)=0.
      cnrk(2,1)=8./15.
      cnrk(3,1)=2./3.
c
c     RK4 (4 stage) Van der Houwen-type
c
      anrk(1,2)=8./17.
      anrk(2,2)=17./60.
      anrk(3,2)=5./12.
      anrk(4,2)=3./4.

      bnrk(1,2)=0.
      bnrk(2,2)=-15./68.
      bnrk(3,2)=-17./60.
      bnrk(4,2)=-5./12.

      cnrk(1,2)=0.
      cnrk(2,2)=8./17.
      cnrk(3,2)=8./15.
      cnrk(4,2)=2./3.

      end subroutine rkcoeff
