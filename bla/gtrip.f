c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine gtrip(tc,zsc,zl,
     &     tdt,ttzc,nzt,tamps,tampt,it,seed,ntdt,prezr)
c
c     Generate the trip forcing dependence on the z direction and time
c
      implicit none

      include 'par.f'

      integer it,seed,ntdt
      real tdt,ttzc(nzp+2,4),tamps,tampt
      real tc,zsc,zl,nzt,nzt_steady
c
c     Local variables
c
      real prezr(nzp+15)
      integer z,i
      real p,b
      real hr,arg
      real pi
      real w(nzp+2),ttz(nzp+2)
      parameter (pi = 3.1415926535897932385)
c
c     ttzx(z,i):
c     i=1: actual z depencence
c     i=2: steady z dependence
c     i=3: "old" z dependence
c     i=4: "new" z depncence
c     The actual z dependence is computed as a weighted superposition of
c     ttzx(z,3) and ttzx(z,4). If necessary, a new ttzx(z,4) is computed,
c     and ttzx(z,3) is reset with the values of ttzx(z,4).
c
c     Generate the time independent part ttzc(z,2)
c     at first iteration
c
      if (it.eq.1) then
c
c     Get random distribution and rescale
c
         if (nzt.le.0) then
            write(*,*) 'Trip forcing: nzt must be at least 1:',nzt
            call stopnow(45453)
         end if
         if (nzt.gt.nzp) then
            write(*,*) 'Trip forcing: nzt must be .le. than nzp:',nzt
            call stopnow(45453)
         end if
c
c     Usually nzt_steady (for the steady part) is the same as nzt
c     (for the unsteady part), however for some cases
c     it can be different (e.g. "white" noise on top of 2D disturbances)
c         nzt_steady=20
         nzt_steady = nzt

         call tzran(ttzc(1,2),nzt_steady,seed,prezr)
         do z=1,nzpc
            ttzc(z,2)=tamps/real(nzt)*ttzc(z,2)
         end do
         ntdt=-2
      end if
c
c     Generate new time dependent part if necessary
c     to be able to recreate the trip of restarted simulations,
c     loop from ntdt=-1 up to present trip count.
c
      do i=ntdt+1,int(tc/tdt)
         do z=1,nzpc
            ttzc(z,3)=ttzc(z,4)
         end do
c
c     Get random distribution and rescale
c
         call tzran(ttzc(1,4),nzt,seed,prezr)
         do z=1,nzpc
            ttzc(z,4)=tampt/real(nzt)*ttzc(z,4)
         end do
      end do
c
c     Update trip count as actual time divided by time scale
c
      ntdt=int(tc/tdt)
c
c     Generate the z-dependence of the trip
c     as a smooth transition between old and new trip vectors
c     p is varying from 0 to 1 for a given trip count.
c
      p=(tc-real(ntdt)*tdt)/tdt
      b=p*p*(3.-2.*p)
      do z=1,nzpc
         ttzc(z,1)=ttzc(z,2)+(1.-b)*ttzc(z,3)+b*ttzc(z,4)
      end do
c
c     Fourier shift the trip dependence to zsc
c     (you can't use symmetry and shift sideways)
c
      if (nzp.gt.1.and.nfzsym.eq.0) then
         do z=1,nzpc
            ttz(z)=ttzc(z,1)
         end do
         call vrfftf(ttz,ttz(2),w,w(2),nzp,1,2,1,prezr)
         do z=1,nzp/2+1
c     ttz  = ttz*exp(i*arg) = ttz * ( cos(arg) + i*sin(arg) )
c     ttzr = ttzr*cos(arg) - ttzi*sin(arg)
c     ttzi = ttzi*cos(arg) + ttzr*sin(arg)
            arg=zsc/zl*2.*pi*real(z-1)
            hr=ttz(2*z-1)*cos(arg)-ttz(2*z)*sin(arg)
            ttz(2*z)=ttz(2*z)*cos(arg)+ttz(2*z-1)*sin(arg)
            ttz(2*z-1)=hr
         end do
         call vrfftb(ttz,ttz(2),w,w(2),nzp,1,2,1,prezr)
         do z=1,nzpc
            ttzc(z,1)=ttz(z)*(1./real(nzp))
         end do
      end if

      end subroutine gtrip
