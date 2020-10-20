c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program fsc
c
c     Computes the similarity solution based on user input:
c     Falkner-Skan-Cooke (f,df,d2f,d3f,g,dg,d2g) and
c     passive scalar (th,dth,d2th).
c
c     Writes file 'fsc.dat'
c
      implicit none

      integer runit
      real m,eps,ymax
      integer i,n,scalar,j
      real,allocatable :: pr(:),m1(:)
      real,allocatable :: eta(:)
      real,allocatable :: f0(:),f1(:),f2(:),f3(:)
      real,allocatable :: g0(:),g1(:),g2(:)
      real,allocatable :: t0(:,:),t1(:,:),t2(:,:)

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                  Falkner-Skan-Cooke solver'//
     &     ' $Rev$          *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)

c
c     For interactive use, set runit to 6 and
c     comment out the open statement.
c
c     For batch use, set runit to 10 and uncomment the
c     open statement. [Default behavior]
c
      runit = 10
      open(unit=runit,file='fsc.i')

c
c     Falkner-Skan parameters
c
      write(*,*) 'Give m [0]'
      read(runit,*) m
      write(*,*) m
c
c     Numerical parameters
c
      write(*,*) 'Give resolution [20001]'
      read(runit,*) n
      write(*,*) n
      write(*,*) 'Give eps [1e-15]'
      read(runit,*) eps
      write(*,*) eps
      write(*,*) 'Give ymax [30]'
      read(runit,*) ymax
      write(*,*) ymax
c
c     Temperature parameters
c
      write(*,*) 'Compute scalar [0-N]'
      read(runit,*) scalar
      write(*,*) scalar
      allocate(pr(scalar),m1(scalar))
      do i=1,scalar
         write(*,*) 'Give pr [0.71] No.',i
         read(runit,*) pr(i)
         write(*,*) pr(i)
         write(*,*) 'Give m1 [0] No.',i
         write(*,*) '   m1=0.0 --> Dirichlet b.c. at the wall'
         write(*,*) '   m1=0.5 --> Neumann b.c. at the wall'
c
c     m1 is the exponent of the scalar variation at the wall
c     m1 = 0.0: constant scalar at the wall
c     m1 = 0.5: constant wall gradient (in conjunction with Blasius)
c
         read(runit,*) m1(i)
         write(*,*) m1(i)
      end do
c
c     Allocate arrays
c
      allocate(eta(n))
      allocate(f0(n),f1(n),f2(n),f3(n))
      allocate(g0(n),g1(n),g2(n))
      allocate(t0(n,scalar),t1(n,scalar),t2(n,scalar))
c
c     Compute solution
c
      call fsc_comp(m,pr,m1,ymax,eps,n,
     &     eta,f0,f1,f2,f3,g0,g1,g2,t0,t1,t2,scalar)
c
c     Write out solution in bla format
c
      write(*,*) 'Write solution to fsc.dat'
      open(unit=11,file='fsc.dat',form='unformatted')
c
c     Write out without scalar
c
      if (scalar.eq.0) then
         write(11) 0.,ymax,n,m,ymax-f0(n),scalar
         do i=1,n
            write(11)
     &           eta(i),f0(i),f1(i),f2(i),f3(i),g0(i),g1(i),g2(i)
         end do
      else
c
c     Write out with scalar
c
         write(11) 0.,ymax,n,m,ymax-f0(n),scalar
         write(11) (pr(j),m1(j),j=1,scalar)
         do i=1,n
            write(11)
     &           eta(i),f0(i),f1(i),f2(i),f3(i),g0(i),g1(i),g2(i),
     &           (t0(i,j),t1(i,j),t2(i,j),j=1,scalar)
         end do
      end if
      close(11)
c
c     Write out solution in ASCII format
c
      write(*,*)
      write(*,*) 'Do you want output in text format ? (0) no (1) yes'
      read(runit,*) i
      write(*,*) i
      if (i.ne.0) then
         write(*,*) 'Write solution to fsc.txt'
         open(unit=11,file='fsc.txt',form='formatted')
         if (scalar.eq.0) then
            do i=1,n
               write(11,'(20e26.17e3)') eta(i),f0(i),f1(i),f2(i),f3(i),
     &              g0(i),g1(i),g2(i)
            end do
         else
            do i=1,n
               write(11,'(200e26.17e3)') eta(i),f0(i),f1(i),f2(i),f3(i),
     &              g0(i),g1(i),g2(i),
     &              (t0(i,j),t1(i,j),t2(i,j),j=1,scalar)
            end do
         end if
         close(11)
      end if

      end program fsc


      subroutine fsc_comp(m,pr,m1,ymax,eps,n,
     &     eta,f0,f1,f2,f3,g0,g1,g2,t0,t1,t2,scalar)
c
c     Solves the Falkner-Skan-Cooke similarity equations by
c     a shooting method with Newton-Raphson.
c     Philipp Schlatter, 2006
c
c     Falkner-Skan:        f'''= - f*f'' - beta*(1-f'*f')
c                          f(0) = f'(0) = 0
c                          f'(infty) = 1
c
c     Cooke:               g'' = - f*g'
c                          g(0) = 0
c                          g(infty) = 1
c
c     Scalar:              th'' = - Pr*f*th' + Pr*m1*(2-beta)*f'*th
c                          th(0) = 1
c                          th(infty) = 0
c                          (only computed if scalar=1)
c
c     Parameters:          beta = 2*m/(m+1)
c                          m = beta/(2-beta)
c                          m1 (equivalent to n)
c
c     If the method should not converge, choose a lower box height
c     (e.g. 10 or 5)
c     and rerun. Then -- if converged -- replace the below values for the
c     initial guess s by the value found with the lower box height.
c
c     Another possibility would be to use bisection instead of Newton-Raphson
c     (at the expense of more integrations)
c
      implicit none

      integer i,iter,n,scalar,j
      real ymax,res,p,pr(scalar),eps
      real d0,d1,d1f,d2,d2f,d3,s
      real beta,m,m1(scalar)

      real,allocatable :: f(:,:)
      real,allocatable :: g(:,:)
      real,allocatable :: th(:,:,:)

      real eta(n)
      real f0(n),f1(n),f2(n),f3(n)
      real g0(n),g1(n),g2(n)
      real t0(n,scalar),t1(n,scalar),t2(n,scalar)

      external g_sys

      beta = 2.*m/(m+1.)
      write(*,*) 'm= ',m,' beta= ',beta
      write(*,*) 'eps=',eps,'ymax=',ymax
      write(*,*) 'N=',n
      write(*,*)
c
c     Check m
c
      if (m.lt.-0.0904) then
         write(*,*) 'Minimum m is -0.0904, now: ',m
         stop
      end if
c
c     Allocate arrays
c
      allocate(f (6,n))
      allocate(g (5,n))
      allocate(th(7,n,scalar))
c
c     Create grid
c
      do i=1,n
         eta(i) = ymax * real(i-1)/real(n-1)
      end do
c
c     Newton-Raphson iteration for f
c
      iter = 0
      s = 1.1
      res = 1.
      write(*,*) 'Newton-Raphson iteration for f: s=',s
      do while (res.gt.eps)
         iter=iter+1

         call fsc_f(s,eta,n,f,beta,m1,pr)

         s = s - (f(2,n)-1.) / f(5,n)
         res = abs(f(2,n)-1.)
         write(*,*) iter,s,res

      end do
      write(*,*) 'Number of iterations ',iter,":  f''(0)=",f(3,1)
c
c     Now compute g
c
      g(1,1) = 0.
      g(2,1) = f(3,1)

      g(3,1) = 0.
      g(4,1) = 0.
      g(5,1) = f(3,1)

      do i=2,n
         call rk4(g(1,i-1),g(1,i),eta(i)-eta(i-1),5,
     &        g_sys,beta,m1,pr)
      end do
c
c     Rescale g to match boundary conditions at infinity
c
      p = 1./g(1,n)
      do i=1,n
         g(1,i) = g(1,i) * p
         g(2,i) = g(2,i) * p
      end do
      write(*,*) "g'(0)=",g(2,1)
c
c     Newton-Raphson iteration for th
c
      do j=1,scalar
         iter = 0
         s = 1.
         res = 1.
         write(*,*) 'Newton-Raphson iteration for th',j,': s=',s
         write(*,*) 'm1=',m1(j),' pr=',pr(j)
         do while (res.gt.eps)
            iter=iter+1

            call fsc_th(s,eta,n,th(1,1,j),beta,m1(j),pr(j),f(3,1))

            s = s - (th(1,n,j)) / th(3,n,j)
            res = abs(th(1,n,j))
            write(*,*) iter,s,res
         end do
         write(*,*) 'Number of iterations ',iter,": th'(0)=",th(2,1,j)
      end do
c
c     Compute the different boundary-layer thicknesses
c
      d0 = 0.
      d1 = 0.
      d2 = 0.
      d3 = 0.
      do i=1,n
c
c     99% thickness
c
         if (f(2,i).le.0.99.and.f(2,i+1).gt.0.99) then
            d0 = eta(i)
         end if
c
c     Integration by trapezoidal rule
c
         p=1.
         if (i.eq.1.or.i.eq.n) p=0.5
         d1 = d1 + p*(1.-f(2,i))
         d2 = d2 + p*(1.-f(2,i))*f(2,i)
         d3 = d3 + p*(1.-f(2,i))*f(2,i)**2
      end do
c
c     Rescaling
c
      d1 = d1/real(n-1)*ymax
      d2 = d2/real(n-1)*ymax
      d3 = d3/real(n-1)*ymax
c
c     Alternative values
c
      d1f = ymax-f(1,n)

      write(*,*)
      write(*,*) 'd0  = ',d0
      write(*,*) 'd1  = ',d1
      write(*,*) '    = ',d1f,' (from f)'
      write(*,*) 'd2  = ',d2
      if (m.eq.0) then
         d2f = f(3,1)
         write(*,*) '    = ',d2f,' (from f)'
      end if
      write(*,*) 'd3  = ',d3
      write(*,*) 'H01 = ',d0/d1f
      write(*,*) 'H12 = ',d1/d2
      if (m.eq.0) then
         write(*,*) 'H12 = ',d1f/d2f,' (from f)'
      else
         write(*,*) 'H12 = ',d1f/d2,' (from f)'
      end if
c
c     Now compute f''', g'' and th'' from the RHS
c     and copy all values to new arrays
c
      do i=1,n
         f0(i) = f(1,i)
         f1(i) = f(2,i)
         f2(i) = f(3,i)
         f3(i) = -f(1,i)*f(3,i)-beta*(1-f(2,i)*f(2,i))
         g0(i) = g(1,i)
         g1(i) = g(2,i)
         g2(i) = -g(2,i)*f(1,i)
         do j=1,scalar
            t0(i,j) = th(1,i,j)
            t1(i,j) = th(2,i,j)
            t2(i,j) = -pr(j)*th(2,i,j)*f(1,i)+pr(j)*m1(j)*(2.-beta)*
     &           (th(1,i,j))*f(2,i)
         end do
      end do

      end subroutine fsc_comp


      subroutine fsc_f(s,eta,n,f,beta,m1,pr)

      implicit none

      integer n,i
      real s,beta,pr,m1
      real eta(n),f(6,n)

      external f_sys

      f(1,1) = 0.
      f(2,1) = 0.
      f(3,1) = s

      f(4,1) = 0.
      f(5,1) = 0.
      f(6,1) = 1.

      do i=2,n
         call rk4(f(1,i-1),f(1,i),eta(i)-eta(i-1),6,
     &        f_sys,beta,m1,pr)
      end do

      end subroutine fsc_f


      subroutine fsc_th(s,eta,n,th,beta,m1,pr,d2f0)

      implicit none

      integer n,i
      real s,beta,pr,m1,d2f0
      real eta(n),th(7,n)

      external th_sys

      th(1,1) = 1.
      th(2,1) = s

      th(3,1) = 0.
      th(4,1) = 1.

      th(5,1) = 0.
      th(6,1) = 0.
      th(7,1) = d2f0

      do i=2,n
         call rk4(th(1,i-1),th(1,i),eta(i)-eta(i-1),7,
     &        th_sys,beta,m1,pr)
      end do

      end subroutine fsc_th


      subroutine f_sys(y1,y2,beta,m1,pr)

      implicit none

      real y1(6),y2(6)
      real beta,pr,m1

      y2(1) = y1(2)
      y2(2) = y1(3)
      y2(3) = -y1(1)*y1(3)-beta*(1.-y1(2)*y1(2))

      y2(4) = y1(5)
      y2(5) = y1(6)
      y2(6) = -(y1(1)*y1(6)+y1(3)*y1(4))+2.*beta*y1(5)*y1(2)

      end subroutine f_sys


      subroutine g_sys(y1,y2,beta,m1,pr)

      implicit none

      real y1(5),y2(5)
      real beta,pr,m1

      y2(1) = y1(2)
      y2(2) = -y1(2)*y1(3)

      y2(3) = y1(4)
      y2(4) = y1(5)
      y2(5) = -y1(3)*y1(5)-beta*(1.-y1(4)*y1(4))

      end subroutine g_sys


      subroutine th_sys(y1,y2,beta,m1,pr)

      implicit none

      real y1(7),y2(7)
      real beta,pr,m1

      y2(1) = y1(2)
      y2(2) = -pr*y1(2)*y1(5) + pr*m1*(2.-beta)*y1(1)*y1(6)

      y2(3) = y1(4)
      y2(4) = -pr*y1(5)*y1(4) + pr*m1*(2.-beta)*y1(6)*y1(3)

      y2(5) = y1(6)
      y2(6) = y1(7)
      y2(7) = -y1(5)*y1(7)-beta*(1.-y1(6)*y1(6))

      end subroutine th_sys


      subroutine rk4(y1,y2,h,nn,der,beta,m1,pr)
c
c     Perform one RK4 step with step size h of a n-vector
c     Input:  y1(n)  old step
c             h      step size
c             nn     vector size
c             der    derivative function
c             beta   parameters for der
c             m1     parameters for der
c             pr     parameters for der
c     Output: y2(nn) new step
c
      implicit none

      integer nn

      real y1(nn),y2(nn),y(nn)
      real k1(nn),k2(nn),k3(nn),k4(nn)
      real h,beta,pr,m1

      external der
c
c     First RK substep
c
      y(:) = y1(:)
      call der(y,k1,beta,m1,pr)
      k1(:)=h*k1(:)
c
c     Second RK substep
c
      y(:) = y1(:) + .5*k1(:)
      call der(y,k2,beta,m1,pr)
      k2(:)=h*k2(:)
c
c     Third RK substep
c
      y(:) = y1(:) + .5*k2(:)
      call der(y,k3,beta,m1,pr)
      k3(:)=h*k3(:)
c
c     Fourth RK substep
c
      y(:) = y1(:) + k3(:)
      call der(y,k4,beta,m1,pr)
      k4(:)=h*k4(:)
c
c     Compose solution after full RK step
c
      y2(:) = y1(:) + ( k1(:)+2.*k2(:)+2.*k3(:)+k4(:) ) / 6.

      end subroutine rk4
