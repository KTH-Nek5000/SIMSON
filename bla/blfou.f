c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine blfou(bu1f,bu2f,bu1,bu2,bu1jr0,bu1ji0,bu1jr,bu1ji,
     &     prex,wbr,wbi)
c
c     Fourier transforms the base flow in the x-direction
c
      implicit none

      include 'par.f'

      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bu1f(nxp/2+1,nyp,3+scalar),bu2f(nxp/2+1,nyp,3+scalar)
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
      real bu1jr0(nxp/2+1,3+scalar,3),bu1ji0(nxp/2+1,3+scalar,3)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real prex(nxp+15)

      integer x,y,i,j
c
c     Copy base flow to Fourier transform array and prenormalize
c
      do i=1,3+scalar
         do y=1,nyp
            do x=1,nxp/2
               bu1f(x,y,i)=bu1(x,y,i)*(1./real(nxp))
               bu2f(x,y,i)=bu2(x,y,i)*(1./real(nxp))
            end do
         end do
      end do
      do i=1,3+scalar
         call vrfftf(bu1f(1,1,i),bu2f(1,1,i),
     &        wbr,wbi,nxp,nyp,1,nxp/2+1,prex)
      end do
      do i=1,3
         do j=1,3+scalar
            do x=1,nxp/2+1
               bu1jr0(x,j,i)=bu1jr(x,j,i)
               bu1ji0(x,j,i)=bu1ji(x,j,i)
            end do
         end do
      end do

      end subroutine blfou
