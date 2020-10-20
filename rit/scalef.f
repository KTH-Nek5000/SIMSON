c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
       subroutine scalef(graf,n,ymin,ymax,ivar,pxy,w,prey,re,ur)
c
c     Rescales a graf from outer to wall scaling.
c     only the y-axis is recaled according to the type of
c     variable
c
      implicit none

      include 'par.f'

      integer n,ivar
      real graf(n),ymin,ymax,re
      complex ur(memnx,memny,memnz,4)
      complex pxy(nx/2,nyp),w(nx/2,nyp)
      real prey(nyp*2+15)
c
      real utau,hplus,dudym,sc
      integer y,i
c
c     First calculate dudy mean
c
      call getxyp(pxy,1,1,ur)
      call vchbf(pxy,w,nyp,1,nx,1,prey)
      call rdcheb(pxy,nyp,1,nx)
      call vchbb(pxy,w,nyp,1,nx,1,prey)
      do y=1,nyp
        pxy(1,y)=pxy(1,y)*(2./real(nyp-1))
      end do
      dudym=(abs(real(pxy(1,1)))+abs(real(pxy(1,nyp))))*.5
      utau=sqrt(dudym/re)
      hplus=utau*re
      sc=1./utau
      if ((ivar.ge.4.and.ivar.le.6).or.ivar.eq.9) sc=1./(utau*hplus)
      if (ivar.eq.10) sc=1./(utau**2)
      if (ivar.eq.13.or.ivar.eq.14) sc=1./(utau**2)
c
c     Then rescale
c
      ymin=ymin*sc
      ymax=ymax*sc
      do i=1,n
         graf(i)=graf(i)*sc
      end do

      return

      end subroutine scalef
