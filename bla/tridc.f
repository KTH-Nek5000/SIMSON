c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine tridc(x,q,c,d,e,b,m,n,nrhs,md,nd)
c
c     Solves m n*n pentadiagonal systems ax=b with the two top rows full
c     where only a(i,j)<>0 for i+j even.
c
c     Used by the pressure solver.
c
c     this system splits into two tridiagonal systems with the top row full
c     solves by Gaussian elemination from "bottom up".
c     ex N=17(Ny=16) :
c     n=0:  q1  0 q3  0 q5  0 q7  0 q9  0 q11 0 q13 0 q15 0 q17   x1    b1
c     n=1:   0 q2  0 q4  0 q6  0 q8  0 q10 0 q12 0 q14 0 q16 0    x2    b2
c     n=2:  c1  0 d1  0 e1                                        x3    b3
c     n=3:     c2  0 d2  0 e2                                     x4    b4
c     n=4:        c3  0 d3  0 e3                                  x5    b5
c     n=5:           c4  0 d4  0 e4                               x6    b6
c     n=6:              c5  0 d5  0 e5                            x7    b7
c     n=7:                 c6  0 d6  0 e6                         x8    b8
c     n=8:                    c7  0 d7  0 e7                      x9    b9
c     n=9:                       c8  0 d8  0 e8                   x10   b10
c     n=10:                         c9  0 d9  0 e9                x11 = b11
c     n=11:                            c10 0 d10 0 e10            x12   b12
c     n=12:                               c11 0 d11 0 e11         x13   b13
c     n=13:                                  c12 0 d12 0  0       x14   b14
c     n=14:                                     c13 0 d13 0  0    x15   b15
c     n=15:                                        c14 0 d14 0    x16   b16
c     n=16:                                           c15 0 d15   x17   b17
c     Note that q1=0 for Neumann bc.
c     Each system has nrhs right hand sides
c
      implicit none

      integer m,n,nrhs,md,nd,ns,nl,i,j,k
      real q(md,n+2),c(md,n),d(md,n),e(md,n)
      real x(md,nd,nrhs),b(md,nd,nrhs)

      do ns=1,2
         nl=n-mod(n+ns,2)
c
c     Eliminate upper off-diagonal terms
c
c     For n=nl (note that d(j,nl-2)=1, e(j,nl-4)=0.)
c
         do k=1,nrhs
            do j=1,m
               x(j,ns,k)=b(j,ns,k)
               x(j,nl-2,k)=b(j,nl-2,k)
               x(j,nl,k)=b(j,nl,k)
            end do
         end do
c
c     From n=nl-2 to ns+4
c
         do i=nl-2,ns+4,-2
            do j=1,m
               d(j,i-2)=1./d(j,i-2)
               c(j,i-2)=c(j,i-2)*d(j,i-2)
               d(j,i-4)=d(j,i-4)-c(j,i-2)*e(j,i-4)
            end do
            do k=1,nrhs
               do j=1,m
                  x(j,i,k)=x(j,i,k)*d(j,i-2)
                  x(j,i-2,k)=b(j,i-2,k)-x(j,i,k)*e(j,i-4)
               end do
            end do
         end do
c
c     For n=ns+2
c
         do j=1,m
            d(j,ns)=1./d(j,ns)
            c(j,ns)=c(j,ns)*d(j,ns)
         end do
         do k=1,nrhs
            do j=1,m
               x(j,2+ns,k)=x(j,2+ns,k)*d(j,ns)
            end do
         end do
c
c     Note that the right hand side is now in x(i)
c     Eliminate full top row
c
         do i=nl,ns+2,-2
            do j=1,m
               q(j,i-2)=q(j,i-2)-c(j,i-2)*q(j,i)
            end do
         end do
c
c     Caution unrolled loop
ccj     Because the number of unknowns of the even parts are different from
ccj     that of the odd parts, this caution may be needed. But I really do
ccj     not understand why we need this caution. We may not need it.(?)
ccj      do 2070 i=nl,ns+2,-2
         do i=nl,ns+4,-4
            do k=1,nrhs
               do j=1,m
ccj 2070  x(j,ns,k)=x(j,ns,k)-x(j,i,k)*q(j,i)
                  x(j,ns,k)=x(j,ns,k)-x(j,i,k)*q(j,i)-
     &                 x(j,i-2,k)*q(j,i-2)
               end do
            end do
         end do

         if (mod(nl-(ns+2),4).eq.0) then
            do k=1,nrhs
               do j=1,m
                  x(j,ns,k)=x(j,ns,k)-x(j,2+ns,k)*q(j,ns+2)
               end do
            end do
         end if
c
c     Eliminate lower off-diagonal
c
         do k=1,nrhs
            do j=1,m
c
c     It's for kx=kz=0 in order to avoid floating point exception if Neumann
c
               if (abs(q(j,ns)).eq.0.) then
                  x(j,ns,k)=0.
               else
                  x(j,ns,k)=x(j,ns,k)/q(j,ns)
               end if
            end do
         end do
         do k=1,nrhs
            do i=ns+2,nl,2
               do j=1,m
                  x(j,i,k)=x(j,i,k)-c(j,i-2)*x(j,i-2,k)
               end do
            end do
         end do
      end do

      end subroutine tridc
