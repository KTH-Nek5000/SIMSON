c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine fft2df(boxr,boxi,sym,npl,prex,prez,pres,prea,wr,wi)
c
c     Transform npl planes in a box from physical to Fourier space
c     in x & z direction
c
c     If nfzsym = 1 then select transform symmetry according to flag sym
c
      implicit none

      include 'par.f'

      logical sym
      integer npl
      real boxr((nxp/2+1)*mby,nzd),boxi((nxp/2+1)*mby,nzd)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)

      integer nxy
c
c     First real to half complex transform in x-direction
c
      nxy=(nxp/2+1)*npl-(nxp/2+1-nx/2)
      call vrfftf(boxr,boxi,wr,wi,nxp,nzpc*mby,1,nxp/2+1,prex)
c
c     Then complex transform in z-direction
c
      if (nfzsym.eq.0) then
         call vcfftf(boxr,boxi,wr,wi,nzp,nxy,(nxp/2+1)*mby,1,prez)
      else
         if (sym) then
            call vcffts(boxr,boxi,wr,nzst,nxy,(nxp/2+1)*mby,1,pres)
         else
            call vcftaf(boxr(1,2),boxi(1,2),wr,wi,nzat,nxy,
     &           (nxp/2+1)*mby,1,prea)
         end if
      end if

      end subroutine fft2df
