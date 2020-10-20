c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine comp_lam(nx,nyp,nz,dur,ar,qr,er,omr,ivort)
c
c     dur(x,y,z,u/v/w,dx/dy/dz)
c
      implicit none

      integer nx,nyp,nz,ivort
      real dur(nx/2+1,nyp,nz,3,3)
      real omr(nx/2+1,nyp,nz,3*ivort)
      real ar(nx/2+1,nyp,nz)
      real qr(nx/2+1,nyp,nz)
      real er(nx/2+1,nyp,nz)
      real sr(3,3),or(3,3),a(3,3),w(3)
      integer ii,jj,kk,i,j,k
      real lambda

!$omp parallel do private(k,j,i,sr,or,ii,jj,kk,a,w,lambda)
      do k=1,nz
         do j=1,nyp
            do i=1,nx/2
c
c     Compute Omega and S
c
               sr(1,1) = dur(i,j,k,1,1)
               sr(2,1) = .5*(dur(i,j,k,2,1)+dur(i,j,k,1,2))
               sr(3,1) = .5*(dur(i,j,k,3,1)+dur(i,j,k,1,3))

               sr(1,2) = .5*(dur(i,j,k,1,2)+dur(i,j,k,2,1))
               sr(2,2) = dur(i,j,k,2,2)
               sr(3,2) = .5*(dur(i,j,k,3,2)+dur(i,j,k,2,3))

               sr(1,3) = .5*(dur(i,j,k,1,3)+dur(i,j,k,3,1))
               sr(2,3) = .5*(dur(i,j,k,2,3)+dur(i,j,k,3,2))
               sr(3,3) = dur(i,j,k,3,3)

               or(1,1) = 0.
               or(2,1) = .5*(dur(i,j,k,2,1)-dur(i,j,k,1,2))
               or(3,1) = .5*(dur(i,j,k,3,1)-dur(i,j,k,1,3))

               or(1,2) = .5*(dur(i,j,k,1,2)-dur(i,j,k,2,1))
               or(2,2) = 0.
               or(3,2) = .5*(dur(i,j,k,3,2)-dur(i,j,k,2,3))

               or(1,3) = .5*(dur(i,j,k,1,3)-dur(i,j,k,3,1))
               or(2,3) = .5*(dur(i,j,k,2,3)-dur(i,j,k,3,2))
               or(3,3) = 0.

               if (1.eq.0) then
c
c     Compute enstrophy
c
c                  er(i,j,k) = 
c     &                 (dur(i,j,k,2,3)-dur(i,j,k,3,2))**2 +
c     &                 (dur(i,j,k,3,1)-dur(i,j,k,1,3))**2 +
c     &                 (dur(i,j,k,1,2)-dur(i,j,k,2,1))**2
                  

c
c     Compute Q
c     
c                  qr(i,j,k) = 0.
c                  do ii=1,3
c                     do kk=1,3
c                        qr(i,j,k) = qr(i,j,k) 
c     &                       + or(ii,kk)*or(ii,kk)
c     &                       - sr(ii,kk)*sr(ii,kk)
c                     end do
c                  end do
               end if

c
c     Compute lambda2
c
               do ii=1,3
                  do kk=1,3
                     a(ii,kk) = 0.
                     do jj=1,3
                        a(ii,kk) = a(ii,kk) +
     &                       sr(ii,jj)*sr(jj,kk) + or(ii,jj)*or(jj,kk)
                     end do
                  end do
               end do

               call eigenvalue(a,w)

               if (w(1).ge.w(2)) then
                  if (w(2).ge.w(3)) then
                     lambda = w(2)
                  else
                     if (w(1).ge.w(3)) then
                        lambda = w(3)
                     else
                        lambda = w(1)
                     end if
                  end if
               else
                  if (w(1).ge.w(3)) then
                     lambda = w(1)
                  else
                     if (w(2).ge.w(3)) then
                        lambda = w(3)
                     else
                        lambda = w(2)
                     end if
                  end if
               end if

               ar(i,j,k) = lambda

            end do
         end do
      end do
c
c     Compute vorticity
c
      if (ivort .gt. 0) then
!$omp parallel do private(k,j,i)
         do k=1,nz
            do j=1,nyp
               do i=1,nx/2

                  omr(i,j,k,1) = dur(i,j,k,3,2) - dur(i,j,k,2,3)
                  omr(i,j,k,2) = dur(i,j,k,1,3) - dur(i,j,k,3,1) 
                  omr(i,j,k,3) = dur(i,j,k,2,1) - dur(i,j,k,1,2)

                  omr(i,j,k,4) = sqrt(omr(i,j,k,1)**2+
     &                 omr(i,j,k,2)**2+omr(i,j,k,3)**2)

               end do
            end do
         end do
      end if

      end subroutine comp_lam
