c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine jetinit(wallvr, wallvi, vsuc,tc,dstar,jetdiam,jetx,
     &     jetz,jetmag,jet_prof,zsc,xsc,zl,xl
     &     ,prex,prez,pres,wi,wr,it,my_node,wallthr,wallthi)
c
c     jet_prof = 1 :
c     Imposes a parabolic boundary condition at y = 0:
c
c     vjet(x,y) = jetmag*(1-r^2)exp(-(r/0.7**4))
c
c     wherec
c     r = sqrt(jetx**2+jetz**2)/jetdiam.
c
c     jet_prof = 2 :
c     Imposes a top-hat boundary condition at y = 0:
c
c     vjet(x,y) = 0.5*(1-tanh(5*(r-1/r)))
c
c     where,
c     r = sqrt(jetx**2+jetz**2)/(jetdiam*0.9919).
c
c     The amount of blowing is the velocity ratio: jetmag
c     A suction of the average v-velocity is applied,
c     vsuc = -wallvr(1,1),
c     so the net mass-flux is zero.
c
      implicit none

      include 'par.f'

      real wallvr(nxp/2+1,nzd), wallvi(nxp/2+1, nzd)
      real wallthr(nxp/2+1,nzd,scalar), wallthi(nxp/2+1,nzd,scalar)
      real thmag(scalar)
      real vblow, vsuc, tc, dstar, zsc,xsc, zl, xl, step
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15)
      real wr(nxp/2+1,nzd),wi(nxp/2+1,nzd)
      real jetdiam, jetx, jetz, jetmag
      integer jet_prof, ith

      integer my_node, z,x, nxm,it
      real dx,dz
      real xcr, xci, zc, rr, ri, xdistr, xdisti, zdist
      real width, smoothr, smoothi, alpha
      real rrth, rith
c
c     A ramping of the blowing is made to the velocity ratio,
c     so that the cfl-number does not become to large.
c
      vblow = jetmag*step((tc/dstar-0.)/(1.))
c
c     Ramping back...
c
c      vblow = vblow * (1-step (tc/dstar-3.))
c
c     Without ramping...
c
c      vblow = jetmag

c
c     Allowing for later implementation of arbitrary magnitude of the
c     scalar at jet exit, now set to one
c
      thmag = 1.0*step((tc/dstar-0.)/(1.))
c     thmag = 1.0

      dx = xl/real(nxp)
      dz = zl/real(nzp)
      nxm = nxp/2+1
      width = 0.7
      alpha = 0.9919

      do z = 1,nzp
         zc = real(z-nzp/2-1)*dz+zsc
         do x = 1,nxp/2
            xcr = real(2*x-1-nxm)*dx+xsc !-0.5xl*dstar...
            xcr = xcr-int((xcr+xl/2.)/xl)*xl ! +~0.5*xl*dstar
            xci = real(2*x-nxm)*dx+xsc
            xci = xci-int((xci+xl/2.)/xl)*xl

            xdistr = xcr - jetx*dstar
            xdisti = xci - jetx*dstar
            zdist = zc - (jetz-0.5)*zl

            if (jet_prof.eq.1) then
               rr = sqrt(xdistr**2+zdist**2)*2.0/jetdiam
               ri = sqrt(xdisti**2+zdist**2)*2.0/jetdiam

               smoothr = exp(-(rr/width)**4)
               smoothi = exp(-(ri/width)**4)

               wallvr(x,z) = vblow*(1. - rr**2)*smoothr
               wallvi(x,z) = vblow*(1. - ri**2)*smoothi

            elseif (jet_prof.eq.2) then
               rr = sqrt(xdistr**2+zdist**2)*2.0/(jetdiam*alpha)
               ri = sqrt(xdisti**2+zdist**2)*2.0/(jetdiam*alpha)

               wallvr(x,z) = vblow*0.5*(1. - tanh(5.*(rr-1./rr)))
               wallvi(x,z) = vblow*0.5*(1. - tanh(5.*(ri-1./ri)))
            end if

            if (scalar.gt.0) then
               rrth = sqrt(xdistr**2+zdist**2)*2.0/(jetdiam*alpha)
               rith = sqrt(xdisti**2+zdist**2)*2.0/(jetdiam*alpha)
            end if

            do ith=1,scalar
               wallthr(x,z,ith)=thmag(ith)*0.5*
     &              (1.-tanh(5.*(rrth-1./rrth)))
               wallthi(x,z,ith)=thmag(ith)*0.5*
     &              (1.-tanh(5.*(rith-1./rith)))
            end do

         end do
      end do

c
c     Tranform to Fourier space
c     Real to half complex transform in x-direction first
c     Then complex transform in z direction
c
      call vrfftf(wallvr,wallvi,wr,wi,nxp,nzpc,1,nxp/2+1,prex)
      if (nfzsym.eq.0) then
         call vcfftf(wallvr,wallvi,wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      else
         call vcffts(wallvr,wallvi,wr,nzst,nx/2,nxp/2+1,1,pres)
      end if
c
c     Transform scalar b.c. to Fourier space
c
      do ith=1,scalar
         call vrfftf(wallthr(:,:,ith),wallthi(:,:,ith),
     &   wr,wi,nxp,nzpc,1,nxp/2+1,prex)
         call vcfftf(wallthr(:,:,ith),wallthi(:,:,ith),
     &   wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      end do
c
c     Normalize
c
      do x = 1,nxp/2+1
         do z = 1,nz/2
            wallvr(x,z) = wallvr(x,z)/real(nxp*nzp)
            wallvi(x,z) = wallvi(x,z)/real(nxp*nzp)
         end do
      end do
c
c     Scalar
c
      do ith=1,scalar
         do x = 1,nxp/2+1
            do z = 1,nz/2
               wallthr(x,z,ith) = wallthr(x,z,ith)/real(nxp*nzp)
               wallthi(x,z,ith) = wallthi(x,z,ith)/real(nxp*nzp)
            end do
         end do
      end do

      if (nfzsym.eq.0) then
         do x = 1,nxp/2+1
            do z = nz/2+1,nz
               wallvr(x,z) = wallvr(x,z+nzp-nz)/real(nxp*nzp)
               wallvi(x,z) = wallvi(x,z+nzp-nz)/real(nxp*nzp)
            end do
         end do
      end if
c
c     Scalar
c
      do ith=1,scalar
         do x = 1,nxp/2+1
            do z = nz/2+1,nz
               wallthr(x,z,ith) = wallthr(x,z+nzp-nz,ith)/real(nxp*nzp)
               wallthi(x,z,ith) = wallthi(x,z+nzp-nz,ith)/real(nxp*nzp)
            end do
         end do
      end do
c
c     Apply Suction so that net massflux is zero
c
c      vsuc=-wallvr(1,1)
c      wallvr(1,1)=0.
c
c     write out the undeveloped jet magnitude every time-step
c
      if (my_node.eq.0 .and. vblow.lt.jetmag .and. mod(it,4).eq.0) then
         write(*,'(a,F12.3,a)') 'Jet magnitude:', vblow/jetmag*100.,'%'
c         write(*,*) 'VSUC,WALLVI=',vsuc,wallvi(1,1)
      end if

      end subroutine jetinit
