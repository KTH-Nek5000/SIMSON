c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine trid(x,xc,q,c,d,e,b,bc,m,n,nrhs,md,nd,cim)
c
c     Solves m n*n pentadiagonal systems ax=b with the two top rows full
c     where only a(i,j)<>0 for i+j even
c
c     This system splits into two tridiagonal systems with the top row full
c     ex n=7 :
c        q(1)  0   q(3)  0   q(5)  0   q(7)      xc(1)        bc(1)
c         0   q(2)  0   q(4)  0   q(6)  0        xc(2)        bc(2)
c        c(1)  0   d(1)  0   e(1)                x(1)         b(1)
c             c(2)  0   d(2)  0   e(2)           x(2)    =    b(2)
c                  c(3)  0   d(3)  0   e(3)      x(3)         b(3)
c                       c(4)  0   d(4)  0        x(4)         b(4)
c                            c(5)  0   d(5)      x(5)         b(5)
c     Each system has nrhs right hand sides
c
      implicit none

      integer m,n,nrhs,md,nd
      logical cim
      real q(md,n)
      real c(md,n-2),d(md,n-2),e(md,n-2)
      real x(md,nd,nrhs),b(md,nd,nrhs),xc(md,2,nrhs),bc(md,2,nrhs)

      integer ns,nl,i,j,k

      if (cim) then
         do ns=1,2
            if (ns.eq.1) then
               nl=(n+1)/2*2-1
            else
               nl=n/2*2
            end if
c
c     Set initial solutions
c
            do k=1,nrhs
               do j=1,m
                  x(j,nl-2,k)=b(j,nl-2,k)
                  xc(j,ns,k)=bc(j,ns,k)
               end do
            end do
c
c     Eliminate upper diagonal
c
            do i=nl,ns+4,-2
               do j=1,m
                  d(j,i-2)=1./d(j,i-2)
                  c(j,i-2)=c(j,i-2)*d(j,i-2)
                  d(j,i-4)=d(j,i-4)-c(j,i-2)*e(j,i-4)
               end do
               do k=1,nrhs
                  do j=1,m
                     x(j,i-2,k)=x(j,i-2,k)*d(j,i-2)
                     x(j,i-4,k)=b(j,i-4,k)-x(j,i-2,k)*e(j,i-4)
                  end do
               end do
            end do
            do j=1,m
               d(j,ns)=1./d(j,ns)
               c(j,ns)=c(j,ns)*d(j,ns)
            end do
            do k=1,nrhs
               do j=1,m
                  x(j,ns,k)=x(j,ns,k)*d(j,ns)
               end do
            end do
c
c     Note that the right hand side now is in x
c     eliminate full row
c
            do i=nl,ns+2,-2
               do j=1,m
                  q(j,i-2)=q(j,i-2)-q(j,i)*c(j,i-2)
               end do
            end do
c
c     Caution unrolled loop
c
            do i=nl,ns+4,-4
               do k=1,nrhs
                  do j=1,m
                     xc(j,ns,k)=xc(j,ns,k)-q(j,i)*x(j,i-2,k)-
     &                    q(j,i-2)*x(j,i-4,k)
                  end do
               end do
            end do
            if (mod(nl-(ns+2),4).eq.0) then
               do k=1,nrhs
                  do j=1,m
                     xc(j,ns,k)=xc(j,ns,k)-q(j,ns+2)*x(j,ns,k)
                  end do
               end do
            end if
c
c     Eliminate lower diagonal
c
            do k=1,nrhs
               do j=1,m
                  xc(j,ns,k)=xc(j,ns,k)/q(j,ns)
                  x(j,ns,k)=x(j,ns,k)-c(j,ns)*xc(j,ns,k)
               end do
            end do
            do i=ns+4,nl,2
               do k=1,nrhs
                  do j=1,m
                     x(j,i-2,k)=x(j,i-2,k)-c(j,i-2)*x(j,i-4,k)
                  end do
               end do
            end do
         end do
      else
c
c     Eliminate upper off-diagonal terms
c     from i=n to n-3(note that d(j,n-2)=d(j,n-3)=1, e(j,n-4)=e(j,n-5)=0.)
c
         do j=1,m
            d(j,n-5)=1./d(j,n-5)
            c(j,n-5)=c(j,n-5)*d(j,n-5)
            d(j,n-4)=1./d(j,n-4)
            c(j,n-4)=c(j,n-4)*d(j,n-4)
         end do
         do k=1,nrhs
            do j=1,m
               x(j,n-3,k)=x(j,n-3,k)*d(j,n-5)
               x(j,n-2,k)=x(j,n-2,k)*d(j,n-4)
            end do
         end do
c
c     From i=n-4 to 3
c
         do i=n-6,1,-1
            do j=1,m
               d(j,i)=1./(d(j,i)-c(j,i+2)*e(j,i))
               c(j,i)=c(j,i)*d(j,i)
            end do
         end do
         do k=1,nrhs
            do i=n-4,3,-1
               do j=1,m
                  x(j,i,k)=(x(j,i,k)-x(j,i+2,k)*e(j,i-2))*d(j,i-2)
               end do
            end do
         end do
c
c     Setup and eliminate full top row(d(j,n)=d(j,n-1)=1)
c
         do j=1,m
            d(j,n-3)=1.-c(j,n-3)
            d(j,n-2)=1.-c(j,n-2)
         end do
         do i=n-4,3,-1
            do j=1,m
               d(j,i)=1.-c(j,i)*d(j,i+2)
            end do
         end do
         do j=1,m
            d(j,1)=1./(1.-c(j,1)*d(j,3))
            d(j,2)=1./(1.-c(j,2)*d(j,4))
         end do
         do k=1,nrhs
            do j=1,m
               x(j,1,k)=bc(j,1,k)-x(j,n,k)
               x(j,2,k)=bc(j,2,k)-x(j,n-1,k)
            end do
         end do
         do k=1,nrhs
            do i=5,n-2,2
               do j=1,m
                  x(j,1,k)=x(j,1,k)-x(j,i,k)*d(j,i)
                  x(j,2,k)=x(j,2,k)-x(j,i-1,k)*d(j,i-1)
               end do
            end do
         end do
         do k=1,nrhs
            do j=1,m
               x(j,1,k)=(x(j,1,k)-x(j,3,k)*d(j,3))*d(j,1)
               x(j,2,k)=x(j,2,k)*d(j,2)
            end do
         end do
c
c     Eliminate lower off-diagonal
c
         do k=1,nrhs
            do i=3,n
               do j=1,m
                  x(j,i,k)=x(j,i,k)-c(j,i-2)*x(j,i-2,k)
               end do
            end do
         end do
      end if

      end subroutine trid
