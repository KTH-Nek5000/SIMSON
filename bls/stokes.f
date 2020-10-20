c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
       subroutine stokes(v,dv,omy,eta,deta,ny,alfa,beta,is,n)
c
c     Calculates eigenfunctions for Stokes' equations in channel flow
c
c     n eigenvalue number, is=0 symmetric v mode, is=1 antisymmetric v mode
c     is=2 symmetric vorticity mode, is=3 antisymmetric vorticity mode
c
       implicit none

       integer n,ny,is
       real v(ny),dv(ny),eta(ny),omy(ny),deta(ny)
       real alfa,beta,k,k2
       integer y
       real s,gam,gam1,e
       real pi,eps
       parameter (pi = 3.1415926535897932385)
       parameter (eps = 1.e-10)

       k2=alfa*alfa+beta*beta
       if (k2.eq.0.and.is.le.1) return
       k=sqrt(k2)
c
c     For the case k2=0 the routine returns the eigenmode for u or v in omy
c
       if (k2.eq.0) k=1
       if (is.eq.0) then
c
c     Symmetric eigenfunctions for vertical velocity
c     simple iteration to solve -gam*tan(gam)=k*tanh(k)
c     initial guess
c
          s=k*tanh(k)
          gam=(real(n)-.5)*pi
 1000     gam1=gam
          gam=real(n)*pi-atan(s/gam1)
          if (abs((gam1-gam)/gam).gt.eps) goto 1000
c
c     Calculate the eigenfunction for v and its derivative
c
          do y=1,ny
             v(y)=-(exp(k*(eta(y)-1.))+exp(-k*(eta(y)+1.)))/
     &            (1.+exp(-2.*k))*cos(gam)+cos(gam*eta(y))
             dv(y)=-k*(exp(k*(eta(y)-1.))-exp(-k*(eta(y)+1.)))/
     &            (1.+exp(-2.*k))*cos(gam)-gam*sin(gam*eta(y))
          end do
       end if
       if (is.eq.1) then
c
c     Anti-symmetric eigenfunctions for vertical velocity
c     simple iteration to solve gam*tan(gam)=k*tanh(k)
c     initial guess
c
         gam=(real(n)+.5)*pi
         s=tanh(k)/k
 1010    gam1=gam
         gam=real(n)*pi+atan(s*gam1)
         if (abs((gam1-gam)/gam).gt.eps) goto 1010
c
c     Calculate the eigenfunction for v and its derivative
c
         do y=1,ny
            v(y)=-(exp(k*(eta(y)-1.))-exp(-k*(eta(y)+1.)))/
     &           (1.-exp(-2.*k))*sin(gam)+sin(gam*eta(y))
            dv(y)=-k*(exp(k*(eta(y)-1.))+exp(-k*(eta(y)+1.)))/
     &           (1.-exp(-2.*k))*sin(gam)+gam*cos(gam*eta(y))
         end do
      end if
c
c     Normalize vertical velocity modes to give energy density
c
      if (is.eq.0.or.is.eq.1) then
         e=0.
         do y=1,ny
            e=e+deta(y)*(v(y)**2+dv(y)**2/k2)
         end do
         do y=1,ny
            v(y)=v(y)*(1./sqrt(.5*e))
            dv(y)=dv(y)*(1./sqrt(.5*e))
         end do
      end if
c
c     Symmetric normal vorticity mode
c
      if (is.eq.2) then
         gam=(real(n)-.5)*pi
c
c     Calculate normalized mode
c
         do y=1,ny
            omy(y)=sqrt(2.)*k*cos(gam*eta(y))
         end do
      end if
c
c     Anti-symmetric normal vorticity mode
c
      if (is.eq.3) then
         gam=(real(n))*pi
c
c     Calculate normalized mode
c
         do y=1,ny
            omy(y)=sqrt(2.)*k*sin(gam*eta(y))
         end do
      end if

      end subroutine stokes
