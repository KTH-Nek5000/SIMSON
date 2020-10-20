c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine blshift(bu1,bu2,bu1jr,bu1ji,bu1jr0,bu1ji0,
     &     bu1f,bu2f,xsc,xl,prex,wbr,wbi)
c
c     Shift the Fourier-transformed base flow and the
c     corresponding boundary derivatives to xs
c
      implicit none

      include 'par.f'

      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bu1f(nxp/2+1,nyp,3+scalar),bu2f(nxp/2+1,nyp,3+scalar)
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
      real bu1jr0(nxp/2+1,3+scalar,3),bu1ji0(nxp/2+1,3+scalar,3)
      real xsc,xl
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real prex(nxp+15)
      real pi
      parameter (pi = 3.1415926535897932385)

      integer x,y,i,j
      real c(nxp/2+1),s(nxp/2+1),arg
c
c     Fourier shift the transformed base flow
c
      do x=1,nxp/2+1
         arg=xsc*2.*pi/xl*real(x-1)
         c(x)=cos(arg)
         s(x)=sin(arg)
      end do
      do i=1,3+scalar
         do y=1,nyp
            do x=1,nxp/2+1
c
c     bu=bu*exp(i*arg)=bu*(cos(arg)+i*sin(arg))
c
               bu1(x,y,i)=bu1f(x,y,i)*c(x)-bu2f(x,y,i)*s(x)
               bu2(x,y,i)=bu2f(x,y,i)*c(x)+bu1f(x,y,i)*s(x)
            end do
         end do
         call vrfftb(bu1(1,1,i),bu2(1,1,i),
     &        wbr,wbi,nxp,nyp,1,nxp/2+1,prex)
      end do
      do i=1,3
         do j=1,3+scalar
            do x=1,nxp/2+1
               bu1jr(x,j,i)=bu1jr0(x,j,i)*c(x)-bu1ji0(x,j,i)*s(x)
               bu1ji(x,j,i)=bu1ji0(x,j,i)*c(x)+bu1jr0(x,j,i)*s(x)
            end do
         end do
      end do

      end subroutine blshift
