c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine solvbc(bcm,bcdr,bcdi,zb)
c
c     Solves for the coffiecients of the homogeneous
c     solutions vh,vh2 and omyh
c     to enforce the desired boundary conditions at y=1
c
      implicit none

      include 'par.f'

      integer zb
      real bcm(nx/2*mbz,3,3),bcdr(nx/2*mbz,3),bcdi(nx/2*mbz,3)

      integer n1,nxz,j,xz
      real tmp,tm

      nxz=nx/2*mbz
      n1=1
      if (zb.eq.1) then
c
c     This is for wavenumber zero
c
         n1=2
         do j=1,3
            bcdr(1,j)=0.0
            bcdi(1,j)=0.0
         end do
      end if
c
c     Solve    bcm (ar,ai) = (bcdr,bcdi) with overwrite of a into bcd
c
c     Second solve with pivoting by row exchange
c
c     Matrix structure   xxx
c                        xxx
c                        xxx
c
      do xz=n1,nxz
         if (abs(bcm(xz,2,3)).gt.abs(bcm(xz,1,3))) then
c
c     Switch row 1 and 2
c
            do j=1,3
               tmp=bcm(xz,1,j)
               bcm(xz,1,j)=bcm(xz,2,j)
               bcm(xz,2,j)=tmp
            end do
            tmp=bcdr(xz,1)
            bcdr(xz,1)=bcdr(xz,2)
            bcdr(xz,2)=tmp
            tmp=bcdi(xz,1)
            bcdi(xz,1)=bcdi(xz,2)
            bcdi(xz,2)=tmp
         end if
      end do
c
c     Null bcm(:,2,3)
c
      do xz=n1,nxz
c
c     Avoid 0/0
c
         if (bcm(xz,2,3).ne.0.) then
            tm=1./bcm(xz,1,3)
            bcm(xz,2,1)=bcm(xz,2,1)-bcm(xz,1,1)*bcm(xz,2,3)*tm
            bcm(xz,2,2)=bcm(xz,2,2)-bcm(xz,1,2)*bcm(xz,2,3)*tm
            bcdr(xz,2)=bcdr(xz,2)-bcdr(xz,1)*bcm(xz,2,3)*tm
            bcdi(xz,2)=bcdi(xz,2)-bcdi(xz,1)*bcm(xz,2,3)*tm
         end if
      end do
c
c     Matrix structure   xxx
c                        xx0
c                        xxx
c
      do xz=n1,nxz
         if (abs(bcm(xz,1,3)).gt.abs(bcm(xz,3,3))) then
c
c     Switch row 1 and 3
c
            do j=1,3
               tmp=bcm(xz,1,j)
               bcm(xz,1,j)=bcm(xz,3,j)
               bcm(xz,3,j)=tmp
            end do
            tmp=bcdr(xz,1)
            bcdr(xz,1)=bcdr(xz,3)
            bcdr(xz,3)=tmp
            tmp=bcdi(xz,1)
            bcdi(xz,1)=bcdi(xz,3)
            bcdi(xz,3)=tmp
         end if
      end do
c
c     Divide row 3 by bcm(:,3,3) and null bcm(:,1,3)
c
      do xz=n1,nxz
         tm=1./bcm(xz,3,3)
         bcm(xz,3,1)=bcm(xz,3,1)*tm
         bcm(xz,3,2)=bcm(xz,3,2)*tm
         bcdr(xz,3)=bcdr(xz,3)*tm
         bcdi(xz,3)=bcdi(xz,3)*tm
         bcm(xz,1,1)=bcm(xz,1,1)-bcm(xz,3,1)*bcm(xz,1,3)
         bcm(xz,1,2)=bcm(xz,1,2)-bcm(xz,3,2)*bcm(xz,1,3)
         bcdr(xz,1)=bcdr(xz,1)-bcdr(xz,3)*bcm(xz,1,3)
         bcdi(xz,1)=bcdi(xz,1)-bcdi(xz,3)*bcm(xz,1,3)
      end do
c
c     Matrix structure   xx0
c                        xx0
c                        xx1
c
c     Switch row 1 and 2
c
      do xz=n1,nxz
         if (abs(bcm(xz,1,2)).gt.abs(bcm(xz,2,2))) then
            do j=1,2
               tmp=bcm(xz,1,j)
               bcm(xz,1,j)=bcm(xz,2,j)
               bcm(xz,2,j)=tmp
            end do
            tmp=bcdr(xz,1)
            bcdr(xz,1)=bcdr(xz,2)
            bcdr(xz,2)=tmp
            tmp=bcdi(xz,1)
            bcdi(xz,1)=bcdi(xz,2)
            bcdi(xz,2)=tmp
         end if
      end do
c
c     Divide row 2 by bcm(:,2,2) and null bcm(:,1,2)
c
      do xz=n1,nxz
         tm=1./bcm(xz,2,2)
         bcm(xz,2,1)=bcm(xz,2,1)*tm
         bcdr(xz,2)=bcdr(xz,2)*tm
         bcdi(xz,2)=bcdi(xz,2)*tm
         bcm(xz,1,1)=bcm(xz,1,1)-bcm(xz,2,1)*bcm(xz,1,2)
         bcdr(xz,1)=bcdr(xz,1)-bcdr(xz,2)*bcm(xz,1,2)
         bcdi(xz,1)=bcdi(xz,1)-bcdi(xz,2)*bcm(xz,1,2)
      end do
c
c     Matrix structure   x00
c                        x10
c                        xx1
c
c     Divide row 1 by bcm(:,1,1)
c
      do xz=n1,nxz
         bcdr(xz,1)=bcdr(xz,1)*(1./bcm(xz,1,1))
         bcdi(xz,1)=bcdi(xz,1)*(1./bcm(xz,1,1))
      end do
c
c     Matrix structure   100
c                        x10
c                        xx1
c
c     Backsubstitution
c
      do xz=n1,nxz
         bcdr(xz,2)=bcdr(xz,2)-bcdr(xz,1)*bcm(xz,2,1)
         bcdi(xz,2)=bcdi(xz,2)-bcdi(xz,1)*bcm(xz,2,1)
         bcdr(xz,3)=bcdr(xz,3)-bcdr(xz,1)*bcm(xz,3,1)-
     &        bcdr(xz,2)*bcm(xz,3,2)
         bcdi(xz,3)=bcdi(xz,3)-bcdi(xz,1)*bcm(xz,3,1)-
     &        bcdi(xz,2)*bcm(xz,3,2)
      end do

      end subroutine solvbc
